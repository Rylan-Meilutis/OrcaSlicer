#include "libslic3r/ThermoelasticSectionPreview.hpp"
#include "libslic3r/ThermoelasticGCodeSections.hpp"
#include "libslic3r/ThermoelasticMaterial.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticLayers.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Eigen/Core>
#include <stdexcept>
#include <memory>
#include <optional>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;

TEST_CASE("Layer preview preserves cancellation from deferred ownership callbacks", "[ThermoelasticSectionPreview]")
{
    GCodeProcessorResult source;
    for (unsigned i = 0; i < 4; ++i) {
        GCodeProcessorResult::MoveVertex move{};
        move.type = i % 2 ? EMoveType::Extrude : EMoveType::Travel;
        move.position = {i % 2 ? 3.f : 1.f, 1.f, 1.f};
        move.width = move.height = 1.f;
        move.mm3_per_mm = 1.;
        move.delta_extruder = i % 2 ? .4f : 0.f;
        move.extruder_id = i / 2;
        move.extrusion_role = erExternalPerimeter;
        move.time = {1.f, 1.f};
        source.moves.push_back(move);
    }
    source.print_statistics.modes[0].time = 4.;
    source.print_statistics.modes[1].time = 4.;
    GCodeSectionScenario scenario(BeadSection::RectangularEnvelope);
    ElasticityMatrix c = ElasticityMatrix::Zero();
    c.diagonal().head<3>().setConstant(1000.);
    c.diagonal().tail<3>().setConstant(500.);
    scenario.materials = {{c}, {c}};
    scenario.filament_to_material = {0, 1};
    scenario.initial_basis_spacing_m = {.001, .0005, .001};
    scenario.study.max_section_height_m = .001;
    scenario.study.relative_tolerance = .05;
    for (unsigned mode = 0; mode < 3; ++mode) {
        unsigned calls = 0;
        const auto owner = [&](double z) -> std::optional<SectionOwnership> {
            SectionOwnership result{z, {}};
            const auto cancel = [&](const auto&) -> std::vector<DepositedSection::Region> {
                ++calls;
                throw ResourceLimit("cancel deferred ownership");
            };
            if (mode == 0) result.resolve_conflict = cancel;
            else if (mode == 1) result.resolve_sampled_conflict = cancel;
            else result.resolve_partition = [&](const SampledSectionConflict&) -> std::shared_ptr<const AuthenticatedSectionPartition> {
                ++calls;
                throw ResourceLimit("cancel deferred partition");
            };
            return result;
        };
        CHECK_THROWS_AS(preview_section_layer(source, 0, scenario, 2000000, 2000000, 200000, owner), ResourceLimit);
        CHECK(calls == 1);
    }
}

TEST_CASE("Static layer previews sample only selected deposition with explicit prescribed strain", "[ThermoelasticSectionPreview]")
{
    GCodeProcessorResult source;
    for (unsigned i = 0; i < 4; ++i) {
        GCodeProcessorResult::MoveVertex move{};
        move.type = i % 2 ? EMoveType::Extrude : EMoveType::Travel;
        move.position = {i < 2 ? (i ? 3.f : 1.f) : (i == 2 ? 3.f : 1.f), 1.f, i < 2 ? 1.f : 2.f};
        move.width = move.height = 1.f;
        move.mm3_per_mm = 1.;
        move.delta_extruder = i % 2 ? .4f : 0.f;
        move.layer_id = i / 2;
        move.extrusion_role = erExternalPerimeter;
        move.time = {1.f, 1.f};
        source.moves.push_back(move);
    }
    source.print_statistics.modes[0].time = 4.;
    source.print_statistics.modes[1].time = 4.;
    GCodeSectionScenario scenario(BeadSection::RectangularEnvelope);
    ElasticityMatrix c = ElasticityMatrix::Zero();
    c.diagonal().head<3>().setConstant(1000.);
    c.diagonal().tail<3>().setConstant(500.);
    scenario.materials = {{c}};
    scenario.materials[0].eigenstrain[2] = -.001;
    scenario.filament_to_material = {0};
    scenario.initial_basis_spacing_m = {.001,.0005,.001};
    scenario.study.max_section_height_m = .001;
    scenario.study.relative_tolerance = .05;
    const auto first = preview_section_layer(source, 0, scenario);
    INFO(first.detail);
    REQUIRE(first.state == SectionLayerPreview::State::Ready);
    REQUIRE(first.samples.size() == 1);
    CHECK(first.samples.front().id == 1);
    CHECK_THAT(first.observation_time_s, Catch::Matchers::WithinAbs(2., 1e-10));
    CHECK_THAT(first.samples.front().displacement_m.z(), Catch::Matchers::WithinAbs(-.0000005, 1e-12));
    const auto second = preview_section_layer(source, 1, scenario);
    INFO(second.detail);
    REQUIRE(second.state == SectionLayerPreview::State::Ready);
    REQUIRE(second.samples.size() == 1);
    CHECK(second.samples.front().id == 3);
    CHECK_THAT(second.observation_time_s, Catch::Matchers::WithinAbs(4., 1e-10));
    CHECK_THAT(second.samples.front().displacement_m.z(), Catch::Matchers::WithinAbs(-.0000015, 1e-12));
    CHECK(preview_section_layer(source, 10, scenario).state == SectionLayerPreview::State::MissingLayer);
    scenario.study.max_basis_levels = 1;
    const auto incomplete = preview_section_layer(source, 0, scenario);
    CHECK(incomplete.state == SectionLayerPreview::State::Unconverged);
    CHECK(incomplete.samples.empty());
    CHECK_THROWS_AS(preview_section_layer(source, 0, scenario, 2000000, 2000000, 200000, {}, {},
        [] { throw ResourceLimit("cancel"); }), ResourceLimit);
}
