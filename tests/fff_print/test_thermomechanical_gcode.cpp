#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ThermomechanicalGCode.hpp"
#include "libslic3r/ThermomechanicalSimulation.hpp"
#include "libslic3r/ThermoelasticToolpath.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticLayers.hpp"
#include "libslic3r/ThermoelasticGCodeSections.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/ThermalHistory.hpp"
#include "libslic3r/Utils.hpp"
#include "test_utils.hpp"
#include <Eigen/Core>
#include <boost/nowide/fstream.hpp>
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <optional>
#include <string>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;

TEST_CASE("Parsed G-code feeds bounded peeling qualification without inventing bed contact", "[ThermomechanicalGCode][GCodePeeling]")
{
    const bool was_bbl = GCodeProcessor::s_IsBBLPrinter;
    const ScopeGuard restore([was_bbl] { GCodeProcessor::s_IsBBLPrinter = was_bbl; });
    GCodeProcessor::s_IsBBLPrinter = false;
    ScopedTemporaryFile input(".gcode");
    {
        boost::nowide::ofstream stream(input.string());
        stream << "G90\nM83\n;TYPE:Outer wall\n;WIDTH:0.4\n;HEIGHT:0.2\n"
            "G1 X1 Y1 Z0.2 F600\nG1 X2 Y1 E0.04 F600\n";
        REQUIRE(stream.good());
    }
    GCodeProcessor processor;
    processor.apply_config(FullPrintConfig{});
    processor.process_file(input.string());
    const auto path = timed_deposition(processor.get_result());
    REQUIRE(path.state == ToolpathState::Complete);
    REQUIRE(path.deposition.size() == 1);
    GCodeThermalScenario scenario;
    scenario.timing = TimingScenario::AcknowledgedModeledIntervals;
    scenario.mechanical_materials = {{{0,0,0},1000.,0.}};
    scenario.thermal_materials = {{1.,1e6,350.,{{250.,400.,Eigen::Matrix3d::Identity()*1e-5}}}};
    scenario.deposited_temperature_k = {350.};
    scenario.boundaries = {{0.,{},0.,300.}};
    scenario.cooldown_s = 1.;
    const auto &bead = path.deposition.front().geometry;
    scenario.spacing_m = {bead.width_m/2.,bead.width_m/2.,bead.height_m/2.};
    const std::vector<double> times{path.total_s,path.total_s+1.};
    const std::vector<GCodeReleaseRegion> regions{{{.0015,0.},{.01,.01},0,0.,0.,300.}};
    GCodePeelingLimits limits;
    limits.solve.time.run.simulation.thermal.max_step_temperature_error_k = .001;
    const auto study = study_gcode_peeling(processor.get_result(),scenario,times,regions,limits);
    REQUIRE(study.solve.discretization_qualified.has_value());
    CHECK(study.solve.completed_levels == 3);
    CHECK(study.solve.completed_runs == 18);
    CHECK(study.solve.discretization_qualified->history.releases.empty());
    CHECK(study.stationary_extrusion.count == 0);
    CHECK_THAT(study.solve.discretization_qualified->history.final_snapshot.thermal.energy_balance_error_j,WithinAbs(0.,1e-8));
    GCodeProcessorResult stationary;
    stationary.moves = processor.get_result().moves;
    stationary.print_statistics = processor.get_result().print_statistics;
    stationary.filament_diameters = processor.get_result().filament_diameters;
    auto point = stationary.moves.back();
    point.type = EMoveType::Extrude;
    point.internal_only = false;
    point.delta_extruder = 1e-9f;
    point.time[0] = .25f;
    ++point.gcode_id;
    stationary.moves.push_back(point);
    stationary.print_statistics.modes[0].time += .25;
    CHECK_THROWS_AS(simulate_gcode_deposition(stationary,scenario),std::invalid_argument);
    auto approximated = scenario;
    approximated.stationary_extrusion_allowance = StationaryExtrusionAllowance{1e-12,1e-6,1e-3};
    const auto accepted = simulate_gcode_deposition(stationary,approximated);
    CHECK(accepted.stationary_extrusion.count == 1);
    CHECK(accepted.stationary_extrusion.commanded_volume_m3 > 0.);
    const auto original = simulate_gcode_deposition(processor.get_result(),scenario);
    CHECK_THAT(accepted.occupied_volume_m3,WithinAbs(original.occupied_volume_m3,1e-20));
    CHECK_THAT(accepted.commanded_path_volume_m3,WithinAbs(original.commanded_path_volume_m3,1e-20));
    stationary.moves.back().delta_extruder = 1e-7f;
    CHECK_THROWS_AS(simulate_gcode_deposition(stationary,approximated),std::invalid_argument);
    stationary.moves.back().delta_extruder = 1e-9f;
    stationary.moves.back().position.z() += .0001f;
    CHECK_THROWS_AS(simulate_gcode_deposition(stationary,approximated),std::invalid_argument);
    stationary.moves.back().position = point.position;
    approximated.stationary_extrusion_allowance->maximum_volume_m3 = 1e-19;
    CHECK_THROWS_AS(simulate_gcode_deposition(stationary,approximated),std::invalid_argument);
    CHECK_THROWS_AS(study_gcode_peeling(processor.get_result(),scenario,times,{regions.front(),regions.front()},limits),std::invalid_argument);
    auto empty = regions;
    empty.front().minimum_m.x() = .1;
    empty.front().maximum_m.x() = .2;
    CHECK_THROWS_AS(study_gcode_peeling(processor.get_result(),scenario,times,empty,limits),std::invalid_argument);
    auto bounded = limits;
    bounded.max_segments = 0;
    CHECK_THROWS_AS(study_gcode_peeling(processor.get_result(),scenario,times,regions,bounded),ResourceLimit);
    scenario.observation_time_s = path.total_s;
    CHECK_THROWS_AS(study_gcode_peeling(processor.get_result(),scenario,times,regions,limits),std::invalid_argument);
    scenario.observation_time_s.reset();
    scenario.timing = TimingScenario::RejectModeledIntervals;
    CHECK_THROWS_AS(study_gcode_peeling(processor.get_result(),scenario,times,regions,limits),std::invalid_argument);
}

TEST_CASE("Seam visualization does not steal subsequent extrusion timing", "[ThermomechanicalGCode]")
{
    const bool was_bbl = GCodeProcessor::s_IsBBLPrinter;
    const ScopeGuard restore([was_bbl] { GCodeProcessor::s_IsBBLPrinter = was_bbl; });
    GCodeProcessor::s_IsBBLPrinter = false;
    ScopedTemporaryFile input(".gcode");
    {
        boost::nowide::ofstream stream(input.string());
        stream << "G90\nM83\n;TYPE:Outer wall\n;WIDTH:0.4\n;HEIGHT:0.2\n"
            "G1 X1 Y1 Z0.2 F600\nG1 X3 Y1 E0.08\nG1 X3 Y3 E0.08\n"
            "G1 X1 Y3 E0.08\nG1 X1 Y1 E0.08\n;TYPE:Inner wall\nG1 X2 Y1 E0.04\n";
        REQUIRE(stream.good());
    }
    GCodeProcessor processor;
    processor.apply_config(FullPrintConfig{});
    processor.process_file(input.string());
    bool seam_found = false;
    for (const auto &move : processor.get_result().moves) {
        if (move.type == EMoveType::Seam) {
            seam_found = true;
            CHECK_THAT(move.time[0], WithinAbs(0., 1e-12));
        }
        if (move.type == EMoveType::Extrude && !move.internal_only) CHECK(move.time[0] > 0.f);
    }
    REQUIRE(seam_found);
    const auto path = timed_deposition(processor.get_result());
    REQUIRE(path.state == ToolpathState::Complete);
    REQUIRE(path.deposition.size() == 5);
    CHECK_THAT(path.deposition.back().geometry.start_m.x(), WithinAbs(.001, 1e-12));
}

TEST_CASE("Parsed deposition evaluates only explicitly acknowledged thermal scenarios", "[ThermomechanicalGCode]")
{
    const bool was_bbl = GCodeProcessor::s_IsBBLPrinter;
    const ScopeGuard restore([was_bbl] { GCodeProcessor::s_IsBBLPrinter = was_bbl; });
    GCodeProcessor::s_IsBBLPrinter = false;
    GCodeProcessor processor;
    processor.apply_config(FullPrintConfig{});
    ScopedTemporaryFile input(".gcode");
    {
        boost::nowide::ofstream stream(input.string());
        stream << "G90\nM83\n;TYPE:Outer wall\n;WIDTH:0.4\n;HEIGHT:0.2\n"
            "G1 X1 Y1 Z0.2 F600\nG1 X3 Y1 E0.08 F600\n";
        REQUIRE(stream.good());
    }
    processor.process_file(input.string());
    GCodeThermalScenario scenario;
    scenario.mechanical_materials = {{{0,0,0},1000.,0.}};
    scenario.thermal_materials = {{1.,1e6,350.,{{250.,400.,Eigen::Matrix3d::Identity()*1e-5}}}};
    scenario.deposited_temperature_k = {350.};
    scenario.boundaries = {{0.,{},0.,300.}};
    const auto path = timed_deposition(processor.get_result());
    REQUIRE(path.state == ToolpathState::Complete);
    REQUIRE(path.deposition.size() == 1);
    const auto &bead = path.deposition.front().geometry;
    // The generic-file parser derives outer-wall rectangular width from the
    // actual E volume, including its documented 1.05 diameter adjustment.
    // Do not assume a WIDTH hint overrides that physical extrusion amount.
    REQUIRE_FALSE(processor.get_result().filament_diameters.empty());
    CHECK_THAT(processor.get_result().filament_diameters[0], WithinAbs(1.75, 1e-6));
    const double radius = 1.05 * 1.75 / 2.;
    const double expected_width_m = .08 * std::acos(-1.) * radius * radius / (2. * .2) * .001;
    CHECK_THAT(bead.width_m, WithinAbs(expected_width_m, 1e-9));
    CHECK_THAT(bead.height_m, WithinAbs(.0002, 1e-9));
    scenario.spacing_m = {bead.width_m / 2., bead.width_m / 2., bead.height_m / 2.};
    CHECK_THROWS_AS(simulate_gcode_deposition(processor.get_result(), scenario), std::invalid_argument);
    scenario.timing = TimingScenario::AcknowledgedModeledIntervals;
    GCodeSimulationLimits excessive;
    excessive.simulation.thermal.solve.max_cells = 200001;
    CHECK_THROWS_AS(simulate_gcode_deposition(processor.get_result(), scenario, excessive), ResourceLimit);
    const auto result = simulate_gcode_deposition(processor.get_result(), scenario);
    REQUIRE(result.simulation.thermal.state == ThermalHistory::State::Converged);
    REQUIRE(result.simulation.mechanical.has_value());
    CHECK(result.uses_modeled_intervals);
    CHECK_FALSE(result.physical_deposition_timing_resolved);
    CHECK(result.uses_nominal_first_hit_cell_activation);
    CHECK(result.segments_without_commanded_volume == 0);
    const double commanded_volume = .08 * std::acos(-1.) * std::pow(1.75 / 2., 2) * 1e-9;
    CHECK_THAT(result.commanded_path_volume_m3, Catch::Matchers::WithinRel(commanded_volume, 1e-6));
    CHECK_THAT(result.nominal_section_path_volume_m3 / result.commanded_path_volume_m3,
        Catch::Matchers::WithinRel(1.05 * 1.05, 1e-6));
    CHECK(result.occupied_volume_m3 > 0.);
    for (const double temperature : result.simulation.thermal.temperatures_k)
        CHECK_THAT(temperature, WithinAbs(350., 1e-7));
    CHECK_THAT(result.simulation.mechanical->equilibrium.elastic_energy_j, WithinAbs(0., 1e-15));
    // A field from one successful solve is not sufficient: independently
    // refine time integration, then require two successive grid comparisons.
    GCodeStudyLimits study_limits;
    CHECK_THROWS_AS(study_gcode_deposition(processor.get_result(), scenario, study_limits), std::invalid_argument);
    study_limits.run.simulation.thermal.max_step_temperature_error_k = 1e-3;
    study_limits.relative_tolerance = .2;
    const auto study = study_gcode_deposition(processor.get_result(), scenario, study_limits);
    INFO("completed runs=" << study.completed_runs << " levels=" << study.completed_levels
        << " thermal state=" << (study.failed_thermal_state ? static_cast<int>(*study.failed_thermal_state) : -1));
    REQUIRE(study.state == GCodeStudyState::Converged);
    REQUIRE(study.qualified.has_value());
    CHECK(study.completed_levels == 3);
    CHECK(study.completed_runs == 6);
    CHECK(study.qualified->simulation.thermal.temporal_error_controlled);
    auto history_scenario = scenario;
    history_scenario.cooldown_s = 1.;
    const std::vector<double> observations{path.total_s + .25, path.total_s + .6, path.total_s + 1.};
    auto history_limits = study_limits;
    // Cooling and reheating need more integration work than the insulated
    // reference above, especially on the finest grid. Keep the same memory
    // and accuracy limits; give this small transient an explicit work budget.
    history_limits.run.simulation.thermal.max_cell_steps = 100000000;
    const auto history = study_gcode_history(processor.get_result(), history_scenario, observations, history_limits);
    INFO("history state=" << static_cast<int>(history.state) << " runs=" << history.completed_runs);
    REQUIRE(history.state == GCodeStudyState::Converged);
    REQUIRE(history.qualified_samples.size() == observations.size());
    CHECK(history.completed_observations == observations.size());
    CHECK(history.completed_runs == 18);
    REQUIRE(history.peak_energy_sample.has_value());
    CHECK(*history.peak_energy_sample == 0);
    CHECK_THAT(history.qualified_samples.front().elastic_energy_j, WithinAbs(0., 1e-15));
    CHECK_THAT(history.qualified_samples.back().elastic_energy_j, WithinAbs(0., 1e-15));
    CHECK_FALSE(history.failed_observation.has_value());
    // A strong transient does not converge on these three grids even when
    // integration has enough work budget. Never publish its sampled peak as
    // qualified just because each individual linear solve succeeds.
    auto transient_scenario = history_scenario;
    transient_scenario.boundaries = {{0., {{1000.,300.}}, 0.,300.},
        {path.total_s + .3, {{1000.,350.}}, 0.,300.}};
    const auto transient_history = study_gcode_history(processor.get_result(), transient_scenario,
        observations, history_limits);
    CHECK(transient_history.state == GCodeStudyState::Unconverged);
    CHECK(transient_history.qualified_samples.empty());
    CHECK_FALSE(transient_history.peak_force_sample.has_value());
    CHECK_FALSE(transient_history.peak_energy_sample.has_value());
    const auto truncated_history = study_gcode_history(processor.get_result(), history_scenario,
        observations, history_limits, 7);
    CHECK(truncated_history.state == GCodeStudyState::BudgetExceeded);
    CHECK(truncated_history.completed_runs == 7);
    CHECK(truncated_history.completed_observations == 1);
    CHECK(truncated_history.failed_observation == std::optional<std::size_t>{1});
    CHECK(truncated_history.qualified_samples.empty());
    CHECK_FALSE(truncated_history.peak_energy_sample.has_value());
    CHECK_FALSE(truncated_history.peak_force_sample.has_value());
    CHECK_THROWS_AS(study_gcode_history(processor.get_result(), scenario, {1.,1.}, study_limits), std::invalid_argument);
    auto conflicting_schedule = history_scenario;
    conflicting_schedule.observation_time_s = observations.front();
    CHECK_THROWS_AS(study_gcode_history(processor.get_result(), conflicting_schedule, observations, study_limits), std::invalid_argument);
    CHECK_THROWS_AS(study_gcode_history(processor.get_result(), history_scenario, observations, study_limits, 64, {},
        [] { throw std::runtime_error("cancelled history"); }), std::runtime_error);
    auto insufficient = study_limits;
    insufficient.run.simulation.thermal.solve.max_cells = 1;
    const auto limited = study_gcode_deposition(processor.get_result(), scenario, insufficient);
    CHECK(limited.state == GCodeStudyState::BudgetExceeded);
    CHECK_FALSE(limited.qualified.has_value());
    insufficient = study_limits;
    insufficient.relative_tolerance = 1e-10;
    insufficient.contact_area_absolute_tolerance_m2 = 0.;
    insufficient.heat_absolute_tolerance_j = 0.;
    const auto unresolved = study_gcode_deposition(processor.get_result(), scenario, insufficient);
    CHECK(unresolved.state == GCodeStudyState::Unconverged);
    CHECK_FALSE(unresolved.qualified.has_value());
    auto cooling = scenario;
    cooling.cooldown_s = 1.;
    cooling.boundaries = {{0., {{500.,300.}}, 0.,300.}};
    auto coarse_time = study_limits;
    // The nominal cap is inactive across short deposition-event intervals.
    // The study must still refine the actual accepted intervals.
    coarse_time.run.simulation.thermal.max_step_temperature_error_k = 1000.;
    coarse_time.relative_tolerance = 1e-10;
    coarse_time.heat_absolute_tolerance_j = 0.;
    const auto time_unresolved = study_gcode_deposition(processor.get_result(), cooling, coarse_time);
    CHECK(time_unresolved.state == GCodeStudyState::Unconverged);
    CHECK(time_unresolved.completed_runs == 2);
    CHECK(time_unresolved.completed_levels == 0);
    CHECK_FALSE(time_unresolved.qualified.has_value());
    CHECK_THROWS_AS(study_gcode_deposition(processor.get_result(), scenario, study_limits, {},
        [] { throw std::runtime_error("cancelled study"); }), std::runtime_error);
    scenario.observation_time_s = .5 * (path.deposition.front().start_s + path.deposition.front().end_s);
    const auto partial = simulate_gcode_deposition(processor.get_result(), scenario);
    REQUIRE(partial.simulation.mechanical.has_value());
    CHECK(partial.simulation.active_voxel_indices.size() < result.simulation.active_voxel_indices.size());
    CHECK(partial.simulation.thermal.deposited_energy_j < result.simulation.thermal.deposited_energy_j);
    CHECK_THAT(partial.observation_time_s, WithinAbs(*scenario.observation_time_s, 1e-12));
    CHECK_THAT(partial.commanded_path_volume_m3,
        Catch::Matchers::WithinRel(.5 * result.commanded_path_volume_m3, 1e-6));
    CHECK_THAT(partial.nominal_section_path_volume_m3,
        Catch::Matchers::WithinRel(.5 * result.nominal_section_path_volume_m3, 1e-6));
    scenario.observation_time_s = path.total_s + 1.;
    CHECK_THROWS_AS(simulate_gcode_deposition(processor.get_result(), scenario), std::invalid_argument);
    scenario.observation_time_s.reset();
    // Explicit hypothetical rounded section, not an inference about the
    // generic parser's display width or this file's physical deposition.
    auto rounded_scenario = scenario;
    rounded_scenario.section = BeadSection::RoundedRectangle;
    const auto rounded = simulate_gcode_deposition(processor.get_result(), rounded_scenario);
    CHECK_THAT(rounded.commanded_path_volume_m3, WithinAbs(result.commanded_path_volume_m3, 1e-20));
    CHECK(rounded.nominal_section_path_volume_m3 < result.nominal_section_path_volume_m3);
    GCodeProcessorResult with_point;
    with_point.moves = processor.get_result().moves;
    with_point.print_statistics = processor.get_result().print_statistics;
    auto point = with_point.moves.back();
    point.type = EMoveType::Extrude;
    point.internal_only = false;
    point.delta_extruder = .00001f;
    point.time[0] = .25f;
    with_point.moves.push_back(point);
    with_point.print_statistics.modes[0].time += .25;
    CHECK_THROWS_AS(simulate_gcode_deposition(with_point, scenario), std::invalid_argument);
    scenario.observation_time_s = path.total_s - 1e-6;
    CHECK_NOTHROW(simulate_gcode_deposition(with_point, scenario));
    scenario.observation_time_s.reset();
    scenario.deposited_temperature_k.clear();
    CHECK_THROWS_AS(simulate_gcode_deposition(processor.get_result(), scenario), std::invalid_argument);
}

TEST_CASE("A supplied export exposes complete modeled deposition timing", "[.][ThermomechanicalGCodeAudit]")
{
    const char *path = std::getenv("ORCA_THERMOMECHANICAL_GCODE");
    if (!path) SKIP("Set ORCA_THERMOMECHANICAL_GCODE for a read-only export audit.");
    REQUIRE(std::filesystem::file_size(path) <= 256u * 1024u * 1024u);
    GCodeProcessor processor;
    processor.process_file(path);
    const auto timed = timed_deposition(processor.get_result(), 2000000, 2000000);
    std::size_t nonzero = 0;
    for (const auto &move : processor.get_result().moves) nonzero += move.time[0] > 0.f;
    for (std::size_t i = 1; i < processor.get_result().moves.size(); ++i) {
        const auto &move = processor.get_result().moves[i];
        if (!move.internal_only && move.type == EMoveType::Extrude && move.time[0] <= 0.f) {
            UNSCOPED_INFO("First missing extrusion duration: index=" << i << " line=" << move.gcode_id
                << " role=" << static_cast<int>(move.extrusion_role) << " deltaE=" << move.delta_extruder
                << " width=" << move.width << " height=" << move.height
                << " distance=" << (move.position - processor.get_result().moves[i-1].position).norm());
            break;
        }
    }
    CAPTURE(static_cast<int>(timed.state), timed.total_s, timed.deposition.size(),
        processor.get_result().moves.size(), nonzero, processor.get_result().print_statistics.modes[0].time);
    CAPTURE(static_cast<int>(timed.failure));
    CAPTURE(timed.unresolved_deposition.size(), timed.spatial_deposition_resolved);
    for (const auto &point : timed.unresolved_deposition) {
        INFO("Unresolved deposition line=" << point.gcode_id << " material=" << point.material
            << " time=" << point.start_s << " E_mm=" << point.extruded_filament_mm
            << " from=" << point.start_m.transpose() << " to=" << point.end_m.transpose()
            << " role=" << static_cast<int>(point.role));
        CHECK(point.extruded_filament_mm > 0.);
    }
    if (timed.failing_move_index) {
        const auto index = *timed.failing_move_index;
        const auto &move = processor.get_result().moves[index];
        CAPTURE(index, move.gcode_id, move.width, move.height, move.fan_speed, move.temperature,
            move.time[0], move.delta_extruder, move.position.x(), move.position.y(), move.position.z());
        if (timed.previous_physical_move_index) {
            const auto prior = *timed.previous_physical_move_index;
            const auto &previous = processor.get_result().moves[prior];
            CAPTURE(prior, static_cast<int>(previous.type), previous.position.x(), previous.position.y(), previous.position.z());
            REQUIRE(timed.state == ToolpathState::Complete);
        }
        REQUIRE(timed.state == ToolpathState::Complete);
    }
    REQUIRE(timed.state == ToolpathState::Complete);
    REQUIRE_FALSE(timed.deposition.empty());
    CHECK_FALSE(timed.deposition_phase_timing_resolved);
    if (std::getenv("ORCA_THERMOMECHANICAL_LAYERS")) {
        // Geometry-only audit: unresolved point extrusions are counted above,
        // not silently admitted into a thermal/mechanical prediction here.
        std::vector<DepositedSegment> segments;
        segments.reserve(timed.deposition.size());
        std::size_t sloped = 0;
        double maximum_slope_z = 0.;
        unsigned first_sloped_line = 0;
        const char *height_limit = std::getenv("ORCA_THERMOMECHANICAL_LAYERS_Z");
        const double maximum_z = height_limit ? std::stod(height_limit) : 10.;
        const char *material_limit = std::getenv("ORCA_THERMOMECHANICAL_LAYERS_MATERIAL");
        const auto selected_material = material_limit ?
            std::optional<std::size_t>(std::stoull(material_limit)) : std::nullopt;
        // Optional single-material *geometry* hypothesis, never a claim that
        // different filament slots have equal mechanical/thermal properties.
        const bool one_material = std::getenv("ORCA_THERMOMECHANICAL_LAYERS_ONE_MATERIAL") != nullptr;
        const auto bead_section = std::getenv("ORCA_THERMOMECHANICAL_ROUNDED_SECTIONS") ?
            BeadSection::RoundedRectangle : BeadSection::RectangularEnvelope;
        for (const auto &entry : timed.deposition) {
            if (std::max(entry.geometry.start_m.z(), entry.geometry.end_m.z()) > maximum_z) continue;
            if (selected_material && entry.geometry.material != *selected_material) continue;
            const double dz = std::abs(entry.geometry.end_m.z() - entry.geometry.start_m.z());
            if (dz > 0.) {
                if (sloped++ == 0) first_sloped_line = entry.gcode_id;
                maximum_slope_z = std::max(maximum_slope_z, dz);
            }
            segments.push_back(entry.geometry);
            segments.back().section = bead_section;
            if (one_material) segments.back().material = 0;
        }
        UNSCOPED_INFO("Sloped segments=" << sloped << " maximum dz=" << maximum_slope_z
            << " first G-code line=" << first_sloped_line);
        CAPTURE(sloped, maximum_slope_z, first_sloped_line);
        REQUIRE_FALSE(segments.empty());
        if (std::getenv("ORCA_THERMOMECHANICAL_SECTION_SWEEP")) {
            // Two midplanes per nominal 0.2 mm layer, for bounded full-height
            // geometry admission only. This is not converged volume quadrature
            // or a force solve, and does not omit sloped deposition segments.
            double top = 0.;
            for (const auto &segment : segments)
                top = std::max({top, segment.start_m.z(), segment.end_m.z()});
            std::vector<double> planes;
            for (std::size_t i = 0; (double(i) + .5) * .0001 < top; ++i)
                planes.push_back((double(i) + .5) * .0001);
            std::size_t peak_points = 0, nonempty = 0;
            double sampled_volume = 0.;
            const auto visits = [&] {
                try {
                    return visit_deposition_sections(segments, planes, 0.,
                        [&](std::size_t, const auto &section) {
                            peak_points = std::max(peak_points, section.boundary_points);
                            nonempty += !section.regions.empty();
                            sampled_volume += section.area_m2 * .0001;
                        });
                } catch (const LayerMaterialOverlap &overlap) {
                    const auto diagnostic = diagnose_section_overlap(timed, overlap, bead_section, 0.);
                    WARN("Overlap evidence state=" << static_cast<int>(diagnostic.state)
                        << " polygon checks=" << diagnostic.polygon_checks);
                    for (const auto &entry : diagnostic.contributors)
                        WARN("Overlap contributor: line=" << entry.gcode_id << " move=" << entry.move_index
                            << " filament=" << entry.filament << " role=" << static_cast<int>(entry.role)
                            << " area_mm2=" << entry.overlap_area_m2 * 1e6);
                    WARN("Section sweep rejected material overlap: z_mm=" << overlap.bottom_m * 1000.
                        << " material=" << overlap.material << " overlap_mm2=" << overlap.area_m2 * 1e6
                        << " completed planes=" << nonempty << "; no completed geometry or force result.");
                    throw;
                }
            }();
            WARN("Geometry-only streaming sweep: single-material hypothesis=" << one_material
                << " planes=" << planes.size() << " nonempty=" << nonempty
                << " bead visits=" << visits << " peak points=" << peak_points
                << " sampled volume_mm3=" << sampled_volume * 1e9
                << "; no mechanical or quadrature qualification.");
            CHECK(nonempty > 0);
            CHECK(sampled_volume > 0.);
        } else if (const char *section_height = std::getenv("ORCA_THERMOMECHANICAL_SECTION_Z")) {
            const auto section = section_deposition(segments, std::stod(section_height), 0.);
            WARN("Geometry-only section: z=" << section_height << " single-material hypothesis=" << one_material
                << " points=" << section.boundary_points << " area_mm2=" << section.area_m2 * 1e6
                << "; no mechanical qualification.");
            CHECK_FALSE(section.regions.empty());
            CHECK(section.area_m2 > 0.);
        } else {
        REQUIRE(maximum_slope_z < 1e-9);
        const auto layers = [&] {
            try {
                return layered_deposition(segments, 0.);
            } catch (const LayerMaterialOverlap &overlap) {
                CAPTURE(overlap.bottom_m, overlap.top_m, overlap.area_m2, overlap.material);
                REQUIRE(overlap.area_m2 <= 0.);
                throw;
            }
        }();
        CAPTURE(layers.slabs.size(), layers.boundary_points, layers.occupied_volume_m3);
        WARN("Geometry-only snapshot: single-material hypothesis=" << one_material
            << " segments=" << segments.size() << " slabs=" << layers.slabs.size()
            << " points=" << layers.boundary_points << " volume_mm3=" << layers.occupied_volume_m3 * 1e9
            << "; not a mechanical qualification or a complete multi-material body.");
        CHECK_FALSE(layers.slabs.empty());
        CHECK(layers.occupied_volume_m3 > 0.);
        }
    }
    // This is an extraction audit only, not a physical material/adhesion test.
    if (std::getenv("ORCA_THERMOMECHANICAL_PREFLIGHT")) {
        // Deliberately synthetic properties: this checks geometry admission and
        // bounded work only, never reports a material/adhesion prediction.
        std::size_t materials = 0;
        double minimum_width = timed.deposition.front().geometry.width_m;
        double minimum_height = timed.deposition.front().geometry.height_m;
        unsigned thinnest_line = 0;
        for (const auto &entry : timed.deposition) {
            materials = std::max(materials,entry.geometry.material+1);
            minimum_width = std::min(minimum_width,entry.geometry.width_m);
            if (entry.geometry.height_m<minimum_height) {
                minimum_height = entry.geometry.height_m;
                thinnest_line = entry.gcode_id;
            }
        }
        CAPTURE(minimum_width,minimum_height,thinnest_line);
        GCodeThermalScenario scenario;
        scenario.timing = TimingScenario::AcknowledgedModeledIntervals;
        scenario.mechanical_materials.assign(materials,{{0,0,0},1000.,0.});
        scenario.thermal_materials.assign(materials,{1.,1e6,350.,{{250.,400.,Eigen::Matrix3d::Zero()}}});
        scenario.deposited_temperature_k.assign(materials,350.);
        scenario.boundaries = {{0.,{},0.,300.}};
        scenario.spacing_m = {.5*minimum_width,.5*minimum_width,.5*minimum_height};
        scenario.stationary_extrusion_allowance = StationaryExtrusionAllowance{1e-12,1e-6,1e-3};
        GCodeSimulationLimits bounded;
        bounded.max_segments = 2000000;
        bounded.simulation.thermal.solve.max_cells = 200000;
        const auto admission = [&] {
            try {
                return simulate_gcode_deposition(processor.get_result(),scenario,bounded);
            } catch (const DepositionBudget &budget) {
                CAPTURE(budget.what(),budget.occupied_cells,budget.cell_visits,budget.segment_index);
                if (budget.segment_index<timed.deposition.size()) {
                    const auto &entry = timed.deposition[budget.segment_index];
                    CAPTURE(entry.gcode_id,entry.start_s,entry.geometry.start_m.z());
                    CHECK(budget.occupied_cells>0);
                }
                throw;
            }
        };
        CHECK_THROWS_AS(admission(),ResourceLimit);
    }
}

TEST_CASE("An observation excludes future bead conflicts and preserves original activation fractions", "[ThermomechanicalGCode]")
{
    GCodeProcessorResult source;
    source.moves.resize(4);
    source.moves[0].position = {3.f,0.f,.2f};
    auto &first = source.moves[1];
    first.type = EMoveType::Extrude;
    first.position = {4.f,0.f,.2f};
    first.width = .4f; first.height = .2f; first.time[0] = 1.f;
    source.moves[2].type = EMoveType::Travel;
    source.moves[2].position = {1.f,0.f,.2f};
    source.moves[2].time[0] = 1.f;
    source.moves[3] = first;
    source.moves[3].extruder_id = 1;
    source.moves[3].time[0] = 3.f;
    source.print_statistics.modes[0].time = 5.;
    GCodeThermalScenario scenario;
    scenario.timing = TimingScenario::AcknowledgedModeledIntervals;
    scenario.mechanical_materials = {{{0,0,0},1000.,0.},{{0,0,0},1000.,0.}};
    const ThermalMechanicalMaterial material{1.,1e6,350.,{{250.,400.,Eigen::Matrix3d::Identity()*1e-5}}};
    scenario.thermal_materials = {material,material};
    scenario.deposited_temperature_k = {350.,350.};
    scenario.boundaries = {{0.,{{100.,300.}},0.,300.}};
    scenario.spacing_m = {.00019,.00019,.00009};
    scenario.observation_time_s = 3.; // Current line has reached x=2, not x=4.
    const auto prefix = simulate_gcode_deposition(source,scenario);
    REQUIRE(prefix.simulation.mechanical.has_value());
    // Independently truncate the input move itself. Identical cooling requires
    // converting clipped local fractions back to the original 3-second move.
    source.moves[3].position.x() = 2.f;
    source.moves[3].time[0] = 1.f;
    source.print_statistics.modes[0].time = 3.;
    const auto reference = simulate_gcode_deposition(source,scenario);
    REQUIRE(prefix.simulation.thermal.temperatures_k.size() == reference.simulation.thermal.temperatures_k.size());
    CHECK_THAT(prefix.simulation.thermal.final_energy_j, WithinAbs(reference.simulation.thermal.final_energy_j,1e-10));
    for (std::size_t i=0;i<prefix.simulation.thermal.temperatures_k.size();++i)
        CHECK_THAT(prefix.simulation.thermal.temperatures_k[i], WithinAbs(reference.simulation.thermal.temperatures_k[i],1e-7));
    source.moves[3].position.x() = 4.f;
    source.moves[3].time[0] = 3.f;
    source.print_statistics.modes[0].time = 5.;
    scenario.observation_time_s.reset();
    CHECK_THROWS_AS(simulate_gcode_deposition(source,scenario),std::runtime_error);
    source.moves[3].extruder_id = 0; // Remove the conflict; future cells still cost memory.
    GCodeSimulationLimits bounded;
    bounded.simulation.thermal.solve.max_cells = prefix.simulation.active_voxel_indices.size();
    scenario.observation_time_s = 3.;
    CHECK_NOTHROW(simulate_gcode_deposition(source,scenario,bounded));
    scenario.observation_time_s.reset();
    CHECK_THROWS_AS(simulate_gcode_deposition(source,scenario,bounded),ResourceLimit);
}
