#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/GCode/DepositionProvenance.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/ThermoelasticSourceDiagnostics.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <memory>
#include <stdexcept>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;

TEST_CASE("Source geometry diagnostics distinguish physical slopes from display vertices", "[ThermoelasticSourceDiagnostics]")
{
    GCodeProcessorResult result;
    auto registry = std::make_shared<DepositionProvenance>();
    const auto id = registry->add({0, 0, 0, 0, 0, DepositionSource::Height::Transformed});
    registry->seal();
    result.deposition_provenance = registry;
    const auto add               = [&](EMoveType type, float x, float z, float height, bool internal = false) {
        GCodeProcessorResult::MoveVertex move;
        move.type                 = type;
        move.position             = Eigen::Vector3f(x, 0.f, z);
        move.height               = height;
        move.width                = .4f;
        move.deposition_source_id = type == EMoveType::Extrude ? id : 0;
        move.internal_only        = internal;
        result.moves.push_back(move);
    };
    add(EMoveType::Travel, 0, .2f, .2f);
    add(EMoveType::Extrude, .5f, .8f, .2f, true);
    add(EMoveType::Seam, .5f, .9f, .2f);
    add(EMoveType::Extrude, 1, .2f, .2f);
    add(EMoveType::Extrude, 2, .3f, .1f);
    add(EMoveType::Extrude, 2, .3f, .1f);
    add(EMoveType::Extrude, 2, .4f, .1f);
    add(EMoveType::Unretract, 2, .4f, .1f);
    result.moves.back().delta_extruder = 1.f;
    const auto diagnostics             = source_geometry_diagnostics(result);
    REQUIRE(diagnostics.state == SourceGeometryDiagnostics::State::Complete);
    CHECK(diagnostics.visited_moves == result.moves.size());
    REQUIRE(diagnostics.sources.size() == 2);
    const auto& source = diagnostics.sources[id];
    CHECK(source.extrusion_moves == 4);
    CHECK(source.horizontal_moves == 1);
    CHECK(source.sloped_moves == 1);
    CHECK(source.stationary_or_unresolved_moves == 2);
    CHECK(source.height_varies());
    REQUIRE(source.top_z_m);
    CHECK_THAT(source.top_z_m->minimum, Catch::Matchers::WithinAbs(.0002, 1e-10));
    CHECK_THAT(source.top_z_m->maximum, Catch::Matchers::WithinAbs(.0004, 1e-10));
    CHECK(diagnostics.sources[0].stationary_or_unresolved_moves == 1);
    // Diagnostics never rewrite the conservative generator qualification.
    CHECK(registry->source(id)->height == DepositionSource::Height::Transformed);
}

TEST_CASE("Source diagnostics require bounded verified provenance and propagate cancellation", "[ThermoelasticSourceDiagnostics]")
{
    GCodeProcessorResult result;
    CHECK(source_geometry_diagnostics(result).state == SourceGeometryDiagnostics::State::MissingRegistry);
    auto registry = std::make_shared<DepositionProvenance>();
    registry->add({});
    registry->seal();
    result.deposition_provenance = registry;
    result.moves.emplace_back();
    result.moves.back().type                 = EMoveType::Extrude;
    result.moves.back().deposition_source_id = 1;
    CHECK(source_geometry_diagnostics(result).sources[1].stationary_or_unresolved_moves == 1);
    CHECK_THROWS_AS(source_geometry_diagnostics(result, 2000001), std::invalid_argument);
    CHECK_THROWS_AS(source_geometry_diagnostics(result, 0), std::invalid_argument);
    CHECK_THROWS_AS(source_geometry_diagnostics(result, 10, [] { throw std::runtime_error("cancel"); }), std::runtime_error);
    result.moves.back().extruder_id = 1;
    CHECK_THROWS_AS(source_geometry_diagnostics(result), std::invalid_argument);
    result.moves.back().extruder_id          = 0;
    result.moves.back().deposition_source_id = 2;
    CHECK_THROWS_AS(source_geometry_diagnostics(result), std::invalid_argument);
    result.moves.emplace_back();
    CHECK_THROWS_AS(source_geometry_diagnostics(result, 1), ResourceLimit);
}
