#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ThermomechanicalLayerPreview.hpp"
#include "libslic3r/ThermomechanicalGCode.hpp"
#include "libslic3r/ThermomechanicalSimulation.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/Point.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <stdexcept>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;
namespace {
void two_layers(GCodeProcessorResult& source)
{
    GCodeProcessorResult::MoveVertex travel{}, deposit{};
    travel.type            = EMoveType::Travel;
    travel.position        = Vec3f(1.f, 1.f, .2f);
    travel.time[0]         = 1.f;
    deposit                = travel;
    deposit.type           = EMoveType::Extrude;
    deposit.extrusion_role = erExternalPerimeter;
    deposit.width          = .4f;
    deposit.height         = .2f;
    deposit.mm3_per_mm     = .08f;
    deposit.delta_extruder = .04f;
    deposit.position.x()   = 2.f;
    source.moves           = {travel, deposit};
    travel.position        = Vec3f(2.f, 1.f, .4f);
    travel.layer_id        = 1;
    deposit.position       = Vec3f(1.f, 1.f, .4f);
    deposit.layer_id       = 1;
    source.moves.push_back(travel);
    source.moves.push_back(deposit);
    source.print_statistics.modes[0].time = 4.;
}
GCodeThermalScenario explicit_scenario()
{
    GCodeThermalScenario scenario;
    scenario.timing               = TimingScenario::AcknowledgedModeledIntervals;
    scenario.mechanical_materials = {{{0, 0, 0}, 1000., 0.}};
    // Synthetic verification data, not a default filament calibration. Pure Z
    // contraction is an affine stress-free solution despite the bonded bed.
    Eigen::Matrix3d expansion        = Eigen::Matrix3d::Zero();
    expansion(2, 2)                  = 1e-5;
    scenario.thermal_materials       = {{1., 1e6, 360., {{250., 400., expansion}}}};
    scenario.deposited_temperature_k = {350.};
    scenario.boundaries              = {{0., {}, 0., 300.}};
    scenario.spacing_m               = {.0002, .0002, .0001};
    return scenario;
}
LayerPreviewLimits preview_limits()
{
    LayerPreviewLimits limits;
    limits.study.relative_tolerance                                  = .2;
    limits.study.run.simulation.thermal.max_step_temperature_error_k = .001;
    return limits;
}
} // namespace

TEST_CASE("Selected layer maps use that observation and exclude future deposited paths", "[ThermomechanicalLayerPreview]")
{
    GCodeProcessorResult source;
    two_layers(source);
    const auto scenario = explicit_scenario();
    const auto limits   = preview_limits();
    const auto first    = preview_mechanical_layer(source, 0, scenario, limits);
    INFO(first.detail);
    REQUIRE(first.state == LayerMechanicalPreview::State::Ready);
    CHECK(first.integrated_discretization_qualified);
    CHECK_THAT(first.observation_time_s, Catch::Matchers::WithinAbs(2., 1e-12));
    REQUIRE(first.samples.size() == 1);
    CHECK(first.samples.front().move_index == 1);
    const auto& point = first.samples.front();
    CHECK_THAT(point.displacement_m.z(), Catch::Matchers::WithinAbs(-1e-4 * point.position_m.z(), 1e-12));
    CHECK_THAT(point.von_mises_pa, Catch::Matchers::WithinAbs(0., 1e-8));
    REQUIRE_FALSE(first.bed_samples.empty());
    for (const auto& bed : first.bed_samples) {
        CHECK(bed.tributary_area_m2 > 0.);
        CHECK_THAT(bed.tensile_traction_pa, Catch::Matchers::WithinAbs(0., 1e-8));
    }
    const auto second = preview_mechanical_layer(source, 1, scenario, limits);
    INFO(second.detail);
    REQUIRE(second.state == LayerMechanicalPreview::State::Ready);
    CHECK_THAT(second.observation_time_s, Catch::Matchers::WithinAbs(4., 1e-12));
    REQUIRE(second.samples.size() == 2);
    CHECK(second.samples.back().move_index == 3);
    CHECK(second.samples.back().displacement_m.z() < first.samples.front().displacement_m.z());
}

TEST_CASE("Unavailable selected layer studies publish neither zero stress nor stale samples", "[ThermomechanicalLayerPreview]")
{
    GCodeProcessorResult source;
    two_layers(source);
    auto scenario      = explicit_scenario();
    auto limits        = preview_limits();
    const auto missing = preview_mechanical_layer(source, 9, scenario, limits);
    CHECK(missing.state == LayerMechanicalPreview::State::MissingLayer);
    CHECK(missing.samples.empty());
    limits.max_samples = 1;
    const auto budget  = preview_mechanical_layer(source, 1, scenario, limits);
    CHECK(budget.state == LayerMechanicalPreview::State::BudgetExceeded);
    CHECK(budget.samples.empty());
    CHECK(budget.bed_samples.empty());
    limits                = preview_limits();
    limits.study.max_runs = 1;
    const auto unfinished = preview_mechanical_layer(source, 0, scenario, limits);
    CHECK(unfinished.state == LayerMechanicalPreview::State::BudgetExceeded);
    CHECK_FALSE(unfinished.integrated_discretization_qualified);
    CHECK(unfinished.samples.empty());
    limits = preview_limits();
    scenario.thermal_materials.clear();
    const auto unsupported = preview_mechanical_layer(source, 0, scenario, limits);
    CHECK(unsupported.state == LayerMechanicalPreview::State::UnsupportedInput);
    CHECK(unsupported.samples.empty());
    CHECK_THROWS_AS(preview_mechanical_layer(source, 0, explicit_scenario(), limits, {}, [] { throw ResourceLimit("cancel"); }),
                    ResourceLimit);
    CHECK_THROWS_AS(preview_mechanical_layer(source, 0, explicit_scenario(), limits, [](unsigned, bool) { throw ResourceLimit("cancel"); }),
                    ResourceLimit);
}
