#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ThermoelasticToolpath.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticSourceOwnership.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/GCode/DepositionSection.hpp"
#include "libslic3r/GCode/DepositionProvenance.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include <stdexcept>
#include <cmath>
#include <array>
#include <memory>
#include <limits>
#include <cstddef>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;
namespace {
TEST_CASE("Complete nominal children preserve authenticated decimal endpoints and terminal heights", "[ThermoelasticToolpath]")
{
    GCodeProcessorResult source;
    source.moves.resize(2);
    auto& start=source.moves[0];
    start.type=EMoveType::Travel;
    start.position={.1f,.2f,.2f};
    start.time[0]=1.f;
    auto& end=source.moves[1];
    end.type=EMoveType::Extrude;
    end.position={1.f,-.222222f,.2f};
    end.time[0]=1.f;
    end.width=.4f;
    end.height=.2f;
    end.delta_extruder=.1f;
    end.extrusion_role=erExternalPerimeter;
    auto registry=std::make_shared<DepositionProvenance>();
    end.deposition_source_id=registry->add({1,0,0,0,0,DepositionSource::Height::Transformed});
    registry->seal();
    source.deposition_provenance=registry;
    auto sections=std::make_shared<DepositionSections>(registry->generation());
    DepositionCommandSection command{{start.position.x(),start.position.y(),start.position.z()},
        {end.position.x(),end.position.y(),end.position.z()},.4f,.2f,.05f,0};
    end.deposition_command_section_id=sections->add(command);
    sections->seal();
    source.deposition_sections=sections;
    end.deposition_design_heights=std::array<double,2>{command.start_height,command.end_height};
    end.deposition_design_width=.4f;
    source.print_statistics.modes[0].time=2.;
    const auto path=timed_deposition(source);
    REQUIRE(path.state==ToolpathState::Complete);
    REQUIRE(path.deposition.size()==1);
    const auto complete=resolve_nominal_deposition(source,path,path.total_s,BeadSection::RectangularEnvelope);
    REQUIRE(complete.segments.size()==1);
    const auto& original=path.deposition.front().geometry;
    const auto& retained=complete.segments.front();
    // Exact identity is the live-command authentication contract, not a
    // numerical closeness check on derived geometry.
    CHECK((retained.end_m.array()==original.end_m.array()).all());
    CHECK(retained.end_height_m==original.end_height_m);
    CHECK(authenticate_live_deposited_source(source,{retained,end.deposition_source_id,0,1})>0);
    const auto partial=resolve_nominal_deposition(source,path,1.5,BeadSection::RectangularEnvelope);
    REQUIRE(partial.segments.size()==1);
    CHECK_THAT((partial.segments.front().end_m-(original.start_m+.5*(original.end_m-original.start_m))).norm(),WithinAbs(0.,1e-15));
    REQUIRE(partial.segments.front().end_height_m);
    CHECK_THAT(*partial.segments.front().end_height_m,
        WithinAbs(original.height_m+.5*(*original.end_height_m-original.height_m),1e-15));
}

void populate(GCodeProcessorResult &source)
{
    source.moves.resize(4);
    source.moves[0].position = {0.f, 0.f, .2f};
    source.moves[1].type = EMoveType::Travel;
    source.moves[1].position = {1.f, 0.f, .2f};
    source.moves[1].time[0] = 3.f;
    auto &extrusion = source.moves[3];
    extrusion.type = EMoveType::Extrude;
    extrusion.position = {3.f, 0.f, .4f};
    extrusion.time[0] = 2.f;
    extrusion.width = .4f;
    extrusion.height = .2f;
    extrusion.extruder_id = 2;
    extrusion.extrusion_role = erPerimeter;
    extrusion.fan_speed = 78.f;
    extrusion.temperature = 220.f;
    source.moves[2] = extrusion;
    source.moves[2].internal_only = true;
    source.moves[2].time[0] = 0.f;
    source.moves[2].position = {2.f, 0.f, .3f};
    source.print_statistics.modes[0].time = 5.;
}

void populate_vertical_join(GCodeProcessorResult &source,
    DepositionCommandSection::Transform transform = DepositionCommandSection::Transform::BrickCourse)
{
    populate(source);
    source.filament_diameters = {1.75f, 1.75f, 1.75f};
    auto registry = std::make_shared<DepositionProvenance>();
    const auto owner = registry->add({1, 0, 0, 0, 2, DepositionSource::Height::Transformed,
                                      DepositionSource::Assignment::NominalRole});
    registry->seal();
    source.deposition_provenance = registry;
    auto sections = std::make_shared<DepositionSections>(registry->generation());
    DepositionCommandSection command{{3., 0., double(.4f)}, {3., 0., double(.5f)}, .45f, .3f, .3f, 2};
    command.transform = transform;
    const auto id = sections->add(command);
    sections->seal();
    source.deposition_sections = sections;
    for (float z : {.45f, .5f}) {
        auto move = source.moves.back();
        move.type = EMoveType::Travel;
        move.position.z() = z;
        move.time[0] = .25f;
        move.delta_extruder = .01f;
        // Travel preview fields are not deposition geometry witnesses.
        move.width = move.height = 0.f;
        move.mm3_per_mm = std::numeric_limits<double>::quiet_NaN();
        move.deposition_source_id = owner;
        move.deposition_command_section_id = id;
        move.deposition_design_width = .45f;
        move.deposition_design_heights = std::array<double, 2>{double(.3f), double(.3f)};
        move.deposition_vertical_extrusion = true;
        source.moves.push_back(move);
    }
    source.print_statistics.modes[0].time = 5.5;
}
}

TEST_CASE("Authenticated vertical joins retain extrusion and timing without inventing a footprint", "[ThermoelasticToolpath]")
{
    GCodeProcessorResult source;
    populate_vertical_join(source);
    const auto path = timed_deposition(source);
    REQUIRE(path.state == ToolpathState::Complete);
    REQUIRE(path.deposition.size() == 1);
    REQUIRE(path.unresolved_deposition.size() == 2);
    CHECK_FALSE(path.spatial_deposition_resolved);
    CHECK_THAT(path.total_s, WithinAbs(5.5, 1e-12));
    for (std::size_t i = 0; i < 2; ++i) {
        const auto &join = path.unresolved_deposition[i];
        CHECK(join.reason == UnresolvedSpatialDeposition::Reason::UnresolvedXYFootprint);
        CHECK(join.material == 2);
        CHECK(join.source_id == source.moves[i + 4].deposition_source_id);
        CHECK(join.command_id == source.moves[i + 4].deposition_command_section_id);
        CHECK_THAT(join.start_s, WithinAbs(5. + .25 * i, 1e-12));
        CHECK_THAT(join.end_s, WithinAbs(5.25 + .25 * i, 1e-12));
        CHECK_THAT(join.extruded_filament_mm, WithinAbs(double(.01f), 1e-12));
        CHECK_THAT((join.end_m - join.start_m).head<2>().norm(), WithinAbs(0., 1e-12));
        CHECK_THAT(join.end_m.z() - join.start_m.z(), WithinAbs(.00005, 1e-10));
    }
    CHECK_THROWS_AS(resolve_nominal_deposition(source, path, 5.5, BeadSection::RoundedRectangle), std::invalid_argument);
    CHECK_THROWS_AS(resolve_nominal_deposition(source, path, 5.5, BeadSection::RoundedRectangle,
        StationaryExtrusionAllowance{1e-12, 1e-6, 1e-6}), std::invalid_argument);
    CHECK(timed_deposition(source, 100, 2).failure == ToolpathFailure::SegmentBudget);
    CHECK_THROWS_AS(timed_deposition(source, 100, 100,
        [] { throw std::runtime_error("cancel"); }), std::runtime_error);
}

TEST_CASE("Vertical deposition witnesses reject stale or altered provenance", "[ThermoelasticToolpath]")
{
    GCodeProcessorResult source;
    populate_vertical_join(source);
    const int mutation = GENERATE(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18);
    auto &move = source.moves[4];
    switch (mutation) {
    case 0: source.deposition_provenance.reset(); break;
    case 1: source.deposition_sections.reset(); break;
    case 2: move.deposition_source_id = 0; break;
    case 3: move.deposition_command_section_id = 0; break;
    case 4: move.extruder_id = 1; break;
    case 5: move.position.x() += .001f; break;
    case 6: move.position.z() = .6f; break;
    case 7: move.deposition_design_width = .5f; break;
    case 8: move.deposition_design_heights.reset(); break;
    case 9: move.delta_extruder = 0.f; break;
    case 10: move.internal_only = true; break;
    case 11: move.type = EMoveType::Unretract; break;
    case 12: move.extrusion_role = erInternalInfill; break;
    case 13: populate_vertical_join(source, DepositionCommandSection::Transform::Ordinary); break;
    case 14: {
        auto stale = std::make_shared<DepositionSections>(source.deposition_provenance->generation() + 1);
        stale->seal();
        source.deposition_sections = stale;
        break;
    }
    case 15: move.type = EMoveType::Seam; break;
    case 16: move.type = EMoveType::Extrude; break;
    case 17: move.deposition_design_heights = std::array<double, 2>{.2, .2}; break;
    case 18: move.delta_extruder = std::numeric_limits<float>::quiet_NaN(); break;
    }
    const auto path = timed_deposition(source);
    CHECK(path.state == ToolpathState::InvalidInput);
    CHECK(path.failure == ToolpathFailure::InvalidVerticalDeposition);
    CHECK(path.deposition.empty());
    CHECK(path.unresolved_deposition.empty());
}

TEST_CASE("Authenticated vertical joins require positive attributed duration", "[ThermoelasticToolpath]")
{
    GCodeProcessorResult source;
    populate_vertical_join(source);
    source.moves[4].time[0] = 0.f;
    const auto path = timed_deposition(source);
    CHECK(path.state == ToolpathState::MissingTiming);
    CHECK(path.failure == ToolpathFailure::MissingExtrusionDuration);
    CHECK(path.unresolved_deposition.empty());
}

TEST_CASE("Unmarked vertical travel and unretraction keep their legacy nondeposition behavior", "[ThermoelasticToolpath]")
{
    GCodeProcessorResult source;
    populate_vertical_join(source);
    const auto type = GENERATE(EMoveType::Travel, EMoveType::Unretract);
    for (std::size_t i = 4; i < source.moves.size(); ++i) {
        source.moves[i].deposition_vertical_extrusion = false;
        source.moves[i].type = type;
    }
    const auto path = timed_deposition(source);
    CHECK(path.state == ToolpathState::Complete);
    CHECK(path.deposition.size() == 1);
    CHECK(path.unresolved_deposition.empty());
    CHECK_THAT(path.total_s, WithinAbs(5.5, 1e-12));
}
TEST_CASE("Modeled deposition preserves original nonplanar paths and cumulative move time", "[ThermoelasticToolpath]")
{
    GCodeProcessorResult source;
    populate(source);
    const auto result = timed_deposition(source);
    REQUIRE(result.state == ToolpathState::Complete);
    REQUIRE(result.deposition.size() == 1);
    CHECK_FALSE(result.deposition_phase_timing_resolved);
    const auto &segment = result.deposition.front();
    CHECK_THAT(segment.start_s, WithinAbs(3., 1e-12));
    CHECK_THAT(segment.end_s, WithinAbs(5., 1e-12));
    CHECK_THAT(segment.geometry.start_m.x(), WithinAbs(.001, 1e-12));
    CHECK_THAT(segment.geometry.end_m.z(), WithinAbs(.0004, 1e-10));
    CHECK(segment.geometry.material == 2);
    CHECK(segment.move_index == 3);
    CHECK_THAT(segment.fan_percent, WithinAbs(78., 1e-12));
    CHECK_THAT(segment.nozzle_target_c, WithinAbs(220., 1e-12));
    CHECK_FALSE(segment.commanded_area_m2.has_value());
    source.moves[3].mm3_per_mm = .08;
    const auto with_volume = timed_deposition(source);
    REQUIRE(with_volume.deposition.front().commanded_area_m2.has_value());
    CHECK_THAT(*with_volume.deposition.front().commanded_area_m2, WithinAbs(.08e-6, 1e-14));
    source.moves[3].mm3_per_mm = -1.;
    CHECK(timed_deposition(source).failure == ToolpathFailure::InvalidCommandedArea);
}
TEST_CASE("Missing processor timing never becomes guessed deposition timing", "[ThermoelasticToolpath]")
{
    GCodeProcessorResult source;
    populate(source);
    source.print_statistics.modes[0].time = 10.;
    CHECK(timed_deposition(source).state == ToolpathState::MissingTiming);
    source.print_statistics.modes[0].time = 5.;
    source.moves[3].time[0] = 0.f;
    CHECK(timed_deposition(source).state == ToolpathState::MissingTiming);
    populate(source);
    CHECK(timed_deposition(source, 2).state == ToolpathState::BudgetExceeded);
    CHECK(timed_deposition(source, 10, 0).deposition.empty());
    CHECK_THROWS_AS(timed_deposition(source, 10, 10,
        [] { throw std::runtime_error("cancel"); }), std::runtime_error);
}

TEST_CASE("Anchored deposition uses precise height while retaining its preview witnesses", "[ThermoelasticToolpath]")
{
    GCodeProcessorResult source;
    populate(source);
    const double offset = GENERATE(0., 10.);
    const bool collapsed = GENERATE(false, true);
    const std::array<double, 2> z{30.123 - offset, (collapsed ? 30.123 : 30.143) - offset};
    auto provenance = std::make_shared<DepositionProvenance>();
    const auto owner = provenance->add({1, 0, 0, 0, 2, DepositionSource::Height::Transformed,
                                      DepositionSource::Assignment::NominalRole});
    provenance->seal();
    source.deposition_provenance = provenance;
    auto sections = std::make_shared<DepositionSections>(provenance->generation());
    DepositionCommandSection command{{1., 0., z[0]}, {3., 0., z[1]}, .4, 0., z[1] - z[0], 2};
    command.nominal_lower_z = std::array<double, 2>{z[0], z[0]};
    const auto id = sections->add(command);
    sections->seal();
    source.deposition_sections = sections;
    source.moves[1].position.z() = float(z[0]);
    auto &move = source.moves[3];
    move.delta_extruder = .00001f;
    move.position.z() = float(z[1]);
    move.deposition_source_id = owner;
    move.deposition_command_section_id = id;
    move.deposition_design_heights = std::array<double, 2>{0., z[1] - z[0]};
    move.deposition_design_width = .4f;
    source.deposition_z_generation = provenance->generation();
    source.deposition_z_spans.push_back({3, owner, id, z, {float(z[0]), float(z[1])}});
    const auto path = timed_deposition(source);
    REQUIRE(path.state == ToolpathState::Complete);
    if (collapsed) {
        CHECK(path.deposition.empty());
        REQUIRE(path.unresolved_deposition.size() == 1);
        const auto &unresolved = path.unresolved_deposition.front();
        CHECK(unresolved.reason == UnresolvedSpatialDeposition::Reason::CollapsedEnvelope);
        CHECK(unresolved.source_id == owner);
        CHECK(unresolved.command_id == id);
        CHECK(unresolved.extruded_filament_mm > 0.);
        CHECK_FALSE(path.spatial_deposition_resolved);
        CHECK_THROWS_AS(resolve_nominal_deposition(source, path, path.total_s,
            BeadSection::RoundedRectangle, StationaryExtrusionAllowance{1e-12, 1e-6, 1e-3}), std::invalid_argument);
    } else {
        REQUIRE(path.deposition.size() == 1);
        const auto &geometry = path.deposition.front().geometry;
        CHECK_THAT(geometry.start_m.z(), WithinAbs(z[0] * .001, 1e-18));
        CHECK_THAT(geometry.end_m.z(), WithinAbs(z[1] * .001, 1e-18));
        REQUIRE(geometry.end_height_m);
        CHECK_THAT(geometry.end_m.z() - *geometry.end_height_m, WithinAbs(z[0] * .001, 1e-17));
    }

    const auto corruption = GENERATE(0, 1, 2, 3, 4, 5, 6, 7);
    switch (corruption) {
    case 0: ++source.deposition_z_spans[0].source_id; break;
    case 1: ++source.deposition_z_spans[0].command_id; break;
    case 2: ++source.deposition_z_generation; break;
    case 3: source.deposition_z_spans[0].preview_z_mm[0] += 1.f; break;
    case 4: source.deposition_z_spans[0].z_mm[1] += .1; break;
    case 5: source.deposition_z_spans.clear(); break;
    case 6: move.deposition_source_id = 0; break;
    case 7: source.deposition_z_spans.push_back(source.deposition_z_spans.front()); break;
    }
    CHECK_THROWS_AS(authenticated_deposition_z(source, 3, 1), std::invalid_argument);
    CHECK(timed_deposition(source).state == ToolpathState::InvalidInput);
}

TEST_CASE("Partial tapered deposition integrates design thickness separately from commanded volume", "[ThermoelasticToolpath]")
{
    GCodeProcessorResult source;
    TimedToolpath path;
    path.state = ToolpathState::Complete;
    path.total_s = 2.;
    TimedDeposition entry;
    entry.geometry = {{0.,0.,0.},{.01,0.,.0002},.0004,0.,0,BeadSection::RectangularEnvelope,.0002};
    entry.end_s = 2.;
    entry.commanded_area_m2 = 1e-7;
    path.deposition.push_back(entry);
    const auto partial = resolve_nominal_deposition(source,path,1.,BeadSection::RectangularEnvelope);
    REQUIRE(partial.segments.size() == 1);
    REQUIRE(partial.segments.front().end_height_m);
    CHECK_THAT(*partial.segments.front().end_height_m,WithinAbs(.0001,1e-12));
    CHECK_THAT(partial.nominal_section_path_volume_m3,WithinAbs(.005*.0004*.0001/2.,1e-18));
    CHECK_THAT(partial.commanded_path_volume_m3,WithinAbs(1e-7*std::sqrt(.005*.005+.0001*.0001),1e-18));
    const auto rounded = resolve_nominal_deposition(source,path,1.,BeadSection::RoundedRectangle);
    const double expected = .005*(.0004*.0001/2.-(1.-std::acos(-1.)/4.)*.0001*.0001/3.);
    CHECK_THAT(rounded.nominal_section_path_volume_m3,WithinAbs(expected,1e-18));
}

TEST_CASE("Stationary extrusion retains timing and explicit unresolved spatial mass", "[ThermoelasticToolpath]")
{
    GCodeProcessorResult source;
    populate(source);
    auto point = source.moves.back();
    point.delta_extruder = .00001f;
    point.time[0] = .25f;
    source.moves.push_back(point);
    source.print_statistics.modes[0].time = 5.25;
    const auto result = timed_deposition(source);
    REQUIRE(result.state == ToolpathState::Complete);
    REQUIRE(result.deposition.size() == 1);
    REQUIRE(result.unresolved_deposition.size() == 1);
    CHECK_FALSE(result.spatial_deposition_resolved);
    const auto &deposit = result.unresolved_deposition.front();
    CHECK_THAT(deposit.start_s, WithinAbs(5., 1e-12));
    CHECK_THAT(deposit.end_s, WithinAbs(5.25, 1e-12));
    CHECK_THAT(deposit.extruded_filament_mm, WithinAbs(.00001, 1e-10));
    CHECK(deposit.material == 2);
    CHECK_THAT((deposit.end_m - deposit.start_m).norm(), WithinAbs(0., 1e-12));
}

TEST_CASE("Collapsed envelopes require an explicit bounded active neighbor omission", "[ThermoelasticToolpath]")
{
    GCodeProcessorResult source;
    source.moves.resize(4);
    source.filament_diameters = {1.75f};
    auto registry = std::make_shared<DepositionProvenance>();
    const auto owner = registry->add({1, 0, 0, 0, 0, DepositionSource::Height::Transformed,
                                      DepositionSource::Assignment::NominalRole});
    registry->seal();
    source.deposition_provenance = registry;
    auto sections = std::make_shared<DepositionSections>(registry->generation());
    const double length = GENERATE(.05, .25);
    const std::array<double, 3> a{0., 0., 8.6}, b{float(length), 0., 8.6}, c{float(length + .1), 0., 8.61};
    DepositionCommandSection collapsed{a, b, .4, 0., 0., 0};
    collapsed.nominal_lower_z = std::array<double, 2>{8.6, 8.6};
    DepositionCommandSection ramp{b, c, .4, 0., c[2] - 8.6, 0};
    ramp.nominal_lower_z = std::array<double, 2>{8.6, 8.6};
    const auto collapsed_id = sections->add(collapsed), ramp_id = sections->add(ramp);
    sections->seal();
    source.deposition_sections = sections;
    source.deposition_z_generation = registry->generation();
    source.moves[0].type = EMoveType::Travel;
    source.moves[0].position = Vec3f(0.f, 0.f, 8.6f);
    source.moves[0].time[0] = 1.f;
    for (std::size_t i = 1; i < source.moves.size(); ++i) {
        auto &move = source.moves[i];
        move.type = EMoveType::Extrude;
        move.extrusion_role = erExternalPerimeter;
        move.deposition_source_id = owner;
        move.width = .4f;
        move.height = .2f;
        move.time[0] = 1.f;
        move.mm3_per_mm = .08;
        move.delta_extruder = i == 1 ? .00001f : (i == 2 ? .00021f : 33.26f);
    }
    source.moves[1].position = Vec3f(float(b[0]), 0.f, 8.6f);
    source.moves[2].position = Vec3f(float(c[0]), 0.f, 8.61f);
    source.moves[3].position = Vec3f(float(c[0] + 1000.), 0.f, 8.61f);
    source.moves[1].deposition_command_section_id = collapsed_id;
    source.moves[2].deposition_command_section_id = ramp_id;
    source.moves[1].deposition_design_heights = std::array<double, 2>{0., 0.};
    source.moves[2].deposition_design_heights = std::array<double, 2>{0., c[2] - 8.6};
    source.moves[1].deposition_design_width = source.moves[2].deposition_design_width = .4f;
    source.deposition_z_spans = {{1, owner, collapsed_id, {8.6, 8.6}, {8.6f, 8.6f}},
                                {2, owner, ramp_id, {8.6, 8.61}, {8.6f, 8.61f}}};
    source.print_statistics.modes[0].time = 4.;
    const auto path = timed_deposition(source);
    REQUIRE(path.state == ToolpathState::Complete);
    REQUIRE(path.unresolved_deposition.size() == 1);
    const CollapsedEnvelopeAllowance allowance{1e-12, 1e-11, 1e-6, .1};
    const auto resolve = [&](const auto &p, double end, const CollapsedEnvelopeAllowance &a) {
        return resolve_nominal_deposition(source, p, end, BeadSection::RoundedRectangle, {}, {}, a);
    };
    CHECK_THROWS_AS(resolve_nominal_deposition(source, path, 4., BeadSection::RoundedRectangle), std::invalid_argument);
    if (length > .2) {
        CHECK_THROWS_AS(resolve(path, 4., allowance), std::invalid_argument);
        return;
    }
    const auto nominal = resolve(path, 4., allowance);
    REQUIRE(nominal.collapsed_envelope.count == 1);
    CHECK(nominal.stationary_extrusion.count == 0);
    CHECK(nominal.segments.size() == path.deposition.size());
    CHECK_THAT(nominal.collapsed_envelope.commanded_volume_m3,
        WithinAbs(std::acos(-1.) * .25 * .00175 * .00175 * double(source.moves[1].delta_extruder) * .001, 1e-25));
    CHECK_THROWS_AS(resolve(path, 2.5, allowance), std::invalid_argument);
    auto strict = allowance;
    strict.maximum_neighbor_fraction = .01;
    CHECK_THROWS_AS(resolve(path, 4., strict), std::invalid_argument);
    strict = allowance;
    strict.maximum_span_volume_m3 = 1e-15;
    CHECK_THROWS_AS(resolve(path, 4., strict), std::invalid_argument);
    strict = allowance;
    strict.maximum_total_volume_m3 = 1e-15;
    CHECK_THROWS_AS(resolve(path, 4., strict), std::invalid_argument);
    strict = allowance;
    strict.maximum_total_fraction = 1e-10;
    CHECK_THROWS_AS(resolve(path, 4., strict), std::invalid_argument);
    strict = allowance;
    strict.maximum_neighbor_fraction = .11;
    CHECK_THROWS_AS(resolve(path, 4., strict), std::invalid_argument);
    auto altered = path;
    altered.unresolved_deposition.front().source_id = owner + 1;
    CHECK_THROWS_AS(resolve(altered, 4., allowance), std::invalid_argument);
    altered = path;
    altered.deposition.front().start_s -= .1;
    CHECK_THROWS_AS(resolve(altered, 4., allowance), std::invalid_argument);
    source.deposition_z_spans[1].source_id = owner + 1;
    CHECK_THROWS_AS(resolve(path, 4., allowance), std::invalid_argument);
}
