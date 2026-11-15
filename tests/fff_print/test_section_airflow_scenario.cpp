#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/SectionAirflowScenario.hpp"
#include "libslic3r/SectionGCodeHistory.hpp"
#include "libslic3r/SectionThermalAssembly.hpp"
#include "libslic3r/SectionIntegration.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/ThermalDeposition.hpp"
#include "libslic3r/ThermalAirflow.hpp"
#include "libslic3r/ThermalAirflowTiming.hpp"
#include "libslic3r/ThermalHistory.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticSections.hpp"
#include "libslic3r/ThermoelasticToolpath.hpp"
#include "libslic3r/Utils.hpp"
#include "test_utils.hpp"
#include <boost/nowide/fstream.hpp>
#include <Eigen/Core>
#include <cstddef>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Slic3r::Thermoelastic;
using namespace Slic3r::ThermalHistory;
using Catch::Matchers::WithinAbs;
using BoundaryElement = Slic3r::ThermalHistory::Element;

namespace {
SectionGCodeScenario scenario()
{
    SectionGCodeScenario result(BeadSection::RectangularEnvelope);
    result.acknowledge_modeled_timing = true;
    result.activation_times_s = {2., 10.};
    result.exposure_times_s = {1., 6., 12.};
    result.cooldown_s = 3.;
    return result;
}
const AirflowExposure empty_law = [](double, const SectionThermalAssembly&, const AirflowState&) {
    return std::vector<BoundaryElement>{};
};
} // namespace

TEST_CASE("Explicit airflow exposure preserves unknown states and uses the last equal-time command", "[SectionAirflowScenario]")
{
    auto study = scenario();
    AirflowState first, second, last;
    first.pwm[0] = 128.;
    second.dock_automatic = true;
    second.unassigned_fan_command = true;
    last.unsupported_firmware_command = true;
    TimedAirflowTrace trace{{{3, 3., first}, {4, 3., second}, {9, 9., last}}, 10., true};
    AirflowState seen;
    double seen_time = -1.;
    apply_timed_airflow(study, trace, [&](double time, const SectionThermalAssembly&, const AirflowState& state) {
        seen = state;
        seen_time = time;
        return std::vector<BoundaryElement>{};
    });
    CHECK(study.activation_times_s == std::vector<double>{2., 10.});
    CHECK(study.exposure_times_s == std::vector<double>{1., 3., 6., 9., 12.});
    SectionThermalAssembly assembly;
    study.exposure(0., assembly);
    CHECK(seen == AirflowState{});
    study.exposure(2.99, assembly);
    CHECK(seen == AirflowState{});
    study.exposure(3., assembly);
    CHECK(seen == second);
    CHECK_FALSE(seen.pwm[0]); // Never interpolate/fill an unknown command.
    study.exposure(8., assembly);
    CHECK(seen == second);
    trace.events.clear(); // Callback owns its trace independently.
    study.exposure(13., assembly);
    CHECK(seen == last);
    CHECK_THAT(seen_time, WithinAbs(13., 1e-12));
    CHECK_THROWS_AS(study.exposure(13.01, assembly), std::invalid_argument);
}

TEST_CASE("Initial airflow commands and empty traces retain explicit boundary semantics", "[SectionAirflowScenario]")
{
    auto study = scenario();
    AirflowState on, seen;
    on.pwm[0] = 255.;
    const auto capture = [&](double, const SectionThermalAssembly&, const AirflowState& state) {
        seen = state;
        return std::vector<BoundaryElement>{};
    };
    apply_timed_airflow(study, {{{1, 0., on}}, 10., true}, capture);
    SectionThermalAssembly assembly;
    study.exposure(0., assembly);
    CHECK(seen == on);
    study.exposure(13., assembly);
    CHECK(seen == on);
    apply_timed_airflow(study, {{}, 10., true}, capture);
    study.exposure(7., assembly);
    CHECK(seen == AirflowState{});
    CHECK(study.exposure_times_s == std::vector<double>{0., 1., 6., 12.});
}

TEST_CASE("Airflow scenario attachment rejects invalid timing and command state atomically", "[SectionAirflowScenario]")
{
    auto study = scenario();
    bool original_called = false;
    study.exposure = [&](double, const SectionThermalAssembly&) {
        original_called = true;
        return std::vector<BoundaryElement>{};
    };
    CHECK_THROWS_AS(apply_timed_airflow(study, {{}, 9., true}, empty_law), std::invalid_argument);
    CHECK_THROWS_AS(apply_timed_airflow(study, {{}, 10., false}, empty_law), std::invalid_argument);
    CHECK_THROWS_AS(apply_timed_airflow(study, {{{1, 3., {}}, {2, 2., {}}}, 10., true}, empty_law), std::invalid_argument);
    CHECK_THROWS_AS(apply_timed_airflow(study, {{{2, 3., {}}, {1, 4., {}}}, 10., true}, empty_law), std::invalid_argument);
    CHECK_THROWS_AS(apply_timed_airflow(study, {{{1, 11., {}}}, 10., true}, empty_law), std::invalid_argument);
    AirflowState invalid;
    invalid.pwm[0] = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(apply_timed_airflow(study, {{{1, 3., invalid}}, 10., true}, empty_law), std::invalid_argument);
    CHECK_THROWS_AS(apply_timed_airflow(study, {{}, 10., true}, {}), std::invalid_argument);
    CHECK(study.exposure_times_s == std::vector<double>{1., 6., 12.});
    study.exposure(0., SectionThermalAssembly{});
    CHECK(original_called);
}

TEST_CASE("Airflow scenario attachment bounds schedules and propagates cancellation", "[SectionAirflowScenario]")
{
    auto study = scenario();
    const TimedAirflowTrace trace{{{1, 3., {}}, {2, 4., {}}}, 10., true};
    CHECK_THROWS_AS(apply_timed_airflow(study, trace, empty_law, 1), std::invalid_argument);
    CHECK_THROWS_AS(apply_timed_airflow(study, trace, empty_law, 100000, [] {
        throw std::runtime_error("cancel");
    }), std::runtime_error);
    study.activation_times_s.clear();
    study.exposure_times_s.clear();
    for (unsigned i = 1; i <= 10000; ++i) study.activation_times_s.push_back(i * .001);
    CHECK_THROWS_AS(apply_timed_airflow(study, trace, empty_law), std::invalid_argument);
    CHECK(study.exposure_times_s.empty());
}

TEST_CASE("Parsed fan changes drive explicit section cooling without changing deposited geometry", "[SectionAirflowScenario]")
{
    using namespace Slic3r;
    const bool was_bbl = GCodeProcessor::s_IsBBLPrinter;
    const ScopeGuard restore([was_bbl] { GCodeProcessor::s_IsBBLPrinter = was_bbl; });
    GCodeProcessor::s_IsBBLPrinter = false;
    const std::vector<std::string> lines{
        "M106 S0", "G90", "M83", ";TYPE:Outer wall", ";WIDTH:2", ";HEIGHT:1",
        "G1 X1 Y1 Z1 F600", "G1 X3 Y1 E1.2 F600", "G1 X3 Y2 F600",
        "M106 S255", "G1 X3 Y1 F600", "G1 X5 Y1 E1.2 F600"};
    ScopedTemporaryFile input(".gcode");
    AirflowTimeline airflow({{0, FanRole::Part}}, FirmwareFanPolicy::GenericM106);
    GCodeReader reader;
    {
        boost::nowide::ofstream stream(input.string());
        for (std::size_t i = 0; i < lines.size(); ++i) {
            stream << lines[i] << '\n';
            reader.parse_line(lines[i], [&](auto&, const auto& line) { airflow.consume(line, i + 1); });
        }
        REQUIRE(stream.good());
    }
    GCodeProcessor processor;
    processor.apply_config(FullPrintConfig{});
    processor.process_file(input.string());
    const auto& source = processor.get_result();
    const auto path = timed_deposition(source);
    REQUIRE(path.state == ToolpathState::Complete);
    REQUIRE(path.deposition.size() == 2);
    const auto trace = time_airflow_events(source, airflow.events());
    REQUIRE(trace.events.size() == 2);
    REQUIRE(trace.events.back().time_s > path.deposition.front().end_s);
    REQUIRE(trace.events.back().time_s < path.total_s);
    SectionGCodeScenario study(BeadSection::RectangularEnvelope);
    study.acknowledge_modeled_timing = true;
    // Deliberately synthetic constants and boundary coefficients: this tests
    // coupling, not a polymer calibration or a PWM-to-convection relation.
    study.thermal_materials = {{Eigen::Matrix3d::Identity() * .2, 2e6, 400.}};
    study.expansion = {{{{250., 450., Eigen::Matrix3d::Identity() * 1e-5}}, Eigen::Matrix3d::Identity(), 400.}};
    study.mechanical_materials = {{Eigen::Matrix<double, 6, 6>::Identity() * 1000.}};
    study.filament_to_material = {0};
    study.spacing_m = {.001, .001, .001};
    study.activation_times_s = {path.deposition.front().end_s, path.total_s};
    study.cooldown_s = .2;
    study.deposition_temperatures_k = [](double) { return std::vector<double>{400.}; };
    apply_timed_airflow(study, trace, empty_law);
    const auto insulated = simulate_section_gcode_history(source, study);
    const auto switch_law = [](double, const SectionThermalAssembly& assembly, const AirflowState& state) {
        std::vector<BoundaryElement> boundaries;
        if (!state.pwm[0] || *state.pwm[0] != 255.) return boundaries;
        const auto add = [&](const auto& face) {
            const auto robin = SectionIntegration::thermal_boundary_element(face.area_moments, 20., 300., assembly.quadrature);
            BoundaryElement boundary;
            boundary.nodes = face.nodes;
            boundary.conductance = robin.conductance;
            boundary.load = robin.load;
            boundary.boundary = true;
            boundary.kind = BoundaryKind::PartFan;
            boundaries.push_back(boundary);
        };
        for (const auto& face : assembly.horizontal_faces)
            if (face.kind != SectionHorizontalFaceKind::Bed) add(face);
        for (const auto& face : assembly.vertical_faces) add(face);
        return boundaries;
    };
    apply_timed_airflow(study, trace, switch_law);
    const auto cooled = simulate_section_gcode_history(source, study);
    REQUIRE(insulated.thermal.state == Slic3r::ThermalHistory::State::Converged);
    REQUIRE(cooled.thermal.state == Slic3r::ThermalHistory::State::Converged);
    REQUIRE(insulated.mechanical.has_value());
    REQUIRE(cooled.mechanical.has_value());
    CHECK(cooled.activated_snapshots == insulated.activated_snapshots);
    CHECK(cooled.activated_snapshots == 2);
    CHECK(cooled.exposure_changes == 1);
    CHECK_THAT(cooled.occupied_volume_m3, WithinAbs(insulated.occupied_volume_m3, 1e-18));
    CHECK_THAT(cooled.thermal.deposited_energy_j, WithinAbs(insulated.thermal.deposited_energy_j, 1e-8));
    CHECK_THAT(insulated.thermal.energy_balance_error_j, WithinAbs(0., 1e-8));
    CHECK_THAT(cooled.thermal.energy_balance_error_j, WithinAbs(0., 1e-8));
    REQUIRE(cooled.thermal.temperatures_k.size() == insulated.thermal.temperatures_k.size());
    const auto temperature_sum = [](const auto& result) {
        return std::accumulate(result.thermal.temperatures_k.begin(), result.thermal.temperatures_k.end(), 0.);
    };
    CHECK(temperature_sum(cooled) < temperature_sum(insulated));
    CHECK(cooled.thermal.boundary_energy_by_kind_j[static_cast<std::size_t>(BoundaryKind::PartFan)] < 0.);
    CHECK(cooled.mechanical->elastic_energy_j > insulated.mechanical->elastic_energy_j);
    CHECK_THAT(trace.total_s, WithinAbs(path.total_s, 1e-12));
    CHECK(timed_deposition(source).deposition.size() == path.deposition.size());
}
