#include "ThermomechanicalGCode.hpp"
#include "ThermomechanicalSimulation.hpp"
#include "ThermalHistory.hpp"
#include "ThermoelasticToolpath.hpp"
#include "ThermoelasticMesh.hpp"
#include "GCode/GCodeProcessor.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
namespace {
struct PreparedDeposition {
    DepositionGrid grid;
    std::vector<CellActivation> activation;
    GCodeSimulationResult result;
};
PreparedDeposition prepare_gcode_deposition(const GCodeProcessorResult &source,
    const GCodeThermalScenario &scenario, const GCodeSimulationLimits &limits,
    const std::function<void()> &checkpoint)
{
    if (limits.max_moves == 0 || limits.max_moves > 2000000 || limits.max_segments == 0 ||
        limits.max_segments > 2000000 || limits.max_cell_visits == 0 || limits.max_cell_visits > 2000000 ||
        limits.simulation.thermal.solve.max_cells == 0 || limits.simulation.thermal.solve.max_cells > 200000 ||
        limits.simulation.max_mechanical_nodes == 0 || limits.simulation.max_mechanical_nodes > 12000 ||
        limits.simulation.max_schedule_entries == 0 || limits.simulation.max_schedule_entries > 10000000 ||
        scenario.boundaries.size() > limits.simulation.thermal.max_events ||
        scenario.boundaries.size() > 10000 ||
        scenario.mechanical_materials.size() > 256 || scenario.thermal_materials.size() > 256 ||
        scenario.deposited_temperature_k.size() > 256)
        throw ResourceLimit("G-code thermal scenario exceeds bounded integration limits.");
    if (scenario.timing != TimingScenario::AcknowledgedModeledIntervals)
        throw std::invalid_argument("Thermal simulation requires explicit acknowledgement of modeled timing.");
    if (scenario.mechanical_materials.empty() ||
        scenario.mechanical_materials.size() != scenario.thermal_materials.size() ||
        scenario.mechanical_materials.size() != scenario.deposited_temperature_k.size() ||
        !std::isfinite(scenario.cooldown_s) || scenario.cooldown_s < 0.)
        throw std::invalid_argument("Incomplete thermal material scenario.");
    for (const double temperature : scenario.deposited_temperature_k)
        if (!std::isfinite(temperature) || temperature <= 0.)
            throw std::invalid_argument("Invalid deposited polymer temperature.");
    const auto path = timed_deposition(source, limits.max_moves, limits.max_segments, checkpoint);
    if (path.state == ToolpathState::BudgetExceeded) throw ResourceLimit("Timed toolpath budget exceeded.");
    if (path.state != ToolpathState::Complete || path.deposition.empty())
        throw std::invalid_argument("Complete timed deposition is unavailable.");
    const double end = scenario.observation_time_s.value_or(path.total_s + scenario.cooldown_s);
    if (!std::isfinite(end) || end < 0. || end > path.total_s + scenario.cooldown_s)
        throw std::invalid_argument("Observation time is outside the modeled timeline.");
    PreparedDeposition prepared;
    auto &result = prepared.result;
    const auto nominal = resolve_nominal_deposition(source, path, end, scenario.section,
        scenario.stationary_extrusion_allowance, checkpoint);
    const auto &segments = nominal.segments;
    const auto &retained_fractions = nominal.retained_fractions;
    result.commanded_path_volume_m3 = nominal.commanded_path_volume_m3;
    result.nominal_section_path_volume_m3 = nominal.nominal_section_path_volume_m3;
    result.segments_without_commanded_volume = nominal.segments_without_commanded_volume;
    result.stationary_extrusion = nominal.stationary_extrusion;
    prepared.grid = voxelize_deposition(segments, scenario.mechanical_materials, scenario.spacing_m,
        scenario.bed_z_m, limits.simulation.thermal.solve.max_cells, checkpoint, limits.max_cell_visits);
    const auto &grid = prepared.grid;
    if (grid.deposition.size() != grid.voxels.size())
        throw std::invalid_argument("Missing cell deposition provenance.");
    auto &activation = prepared.activation;
    activation.reserve(grid.voxels.size());
    for (std::size_t i = 0; i < grid.voxels.size(); ++i) {
        if (checkpoint && i % 1024 == 0) checkpoint();
        const auto &hit = grid.deposition[i];
        if (hit.segment_index >= segments.size() || !std::isfinite(hit.path_fraction) ||
            hit.path_fraction < 0. || hit.path_fraction > 1. ||
            grid.voxels[i].material >= scenario.deposited_temperature_k.size())
            throw std::invalid_argument("Invalid cell deposition provenance.");
        const auto &entry = path.deposition[hit.segment_index];
        activation.push_back({entry.start_s + hit.path_fraction * retained_fractions[hit.segment_index] * (entry.end_s - entry.start_s),
            scenario.deposited_temperature_k[grid.voxels[i].material]});
    }
    result.occupied_volume_m3 = double(grid.voxels.size()) * grid.spacing_m.prod();
    result.modeled_print_time_s = path.total_s;
    result.observation_time_s = end;
    return prepared;
}
} // namespace

GCodeSimulationResult simulate_gcode_deposition(const GCodeProcessorResult &source,
    const GCodeThermalScenario &scenario, const GCodeSimulationLimits &limits,
    const std::function<void()> &checkpoint)
{
    auto prepared = prepare_gcode_deposition(source,scenario,limits,checkpoint);
    prepared.result.simulation = simulate_deposition_grid(prepared.grid,prepared.activation,scenario.thermal_materials,
        scenario.boundaries,scenario.intermaterial_resistance_m2_k_w,prepared.result.observation_time_s,
        limits.simulation, checkpoint);
    return std::move(prepared.result);
}

GCodePeelingStudy study_gcode_peeling(const GCodeProcessorResult &source,
    const GCodeThermalScenario &scenario, const std::vector<double> &times,
    const std::vector<GCodeReleaseRegion> &regions, const GCodePeelingLimits &limits,
    const std::function<void()> &checkpoint)
{
    if (scenario.observation_time_s || times.size()<2 || times.size()>4 || regions.empty() || regions.size()>32)
        throw std::invalid_argument("Invalid G-code peeling observation schedule or regions.");
    for (std::size_t i=0; i<times.size(); ++i)
        if (!std::isfinite(times[i]) || times[i]<0. || (i && times[i]<=times[i-1]))
            throw std::invalid_argument("G-code peeling observations must be finite and strictly increasing.");
    for (const auto &region : regions)
        if (!region.minimum_m.allFinite() || !region.maximum_m.allFinite() ||
            (region.minimum_m.array()>=region.maximum_m.array()).any() ||
            (region.material && *region.material>=scenario.mechanical_materials.size()) ||
            !std::isfinite(region.fracture_energy_j_m2) || region.fracture_energy_j_m2<0. ||
            !std::isfinite(region.gap_coefficient_w_m2_k) || region.gap_coefficient_w_m2_k<0. ||
            !std::isfinite(region.gap_temperature_k) || region.gap_temperature_k<=0.)
            throw std::invalid_argument("Invalid G-code release-region properties.");
    auto bounded = scenario;
    bounded.observation_time_s = times.back();
    GCodeSimulationLimits input_limits;
    input_limits.simulation = limits.solve.time.run.simulation;
    input_limits.max_moves = limits.max_moves;
    input_limits.max_segments = limits.max_segments;
    input_limits.max_cell_visits = limits.max_cell_visits;
    auto prepared = prepare_gcode_deposition(source,bounded,input_limits,checkpoint);
    std::vector<GridReleasePatch> patches;
    for (const auto &region : regions)
        patches.push_back({{},region.fracture_energy_j_m2,region.gap_coefficient_w_m2_k,region.gap_temperature_k});
    for (std::size_t i=0; i<prepared.grid.voxels.size(); ++i) {
        if (checkpoint && i%128==0) checkpoint();
        const auto &voxel = prepared.grid.voxels[i];
        if (voxel.index[2]!=0) continue;
        const Eigen::Vector2d center = prepared.grid.origin_m.head<2>() + prepared.grid.spacing_m.head<2>().cwiseProduct(
            Eigen::Vector2d(double(voxel.index[0])+.5,double(voxel.index[1])+.5));
        bool selected = false;
        for (std::size_t p=0; p<regions.size(); ++p) {
            const auto &region = regions[p];
            if ((region.material && *region.material!=voxel.material) ||
                (center.array()<region.minimum_m.array()).any() || (center.array()>=region.maximum_m.array()).any()) continue;
            if (selected) throw std::invalid_argument("G-code release regions overlap on deposited bed contact.");
            selected = true;
            patches[p].voxels.push_back(i);
        }
    }
    for (const auto &patch : patches)
        if (patch.voxels.empty()) throw std::invalid_argument("G-code release region selects no deposited bed contact.");
    return {study_peeling_space(prepared.grid,prepared.activation,scenario.thermal_materials,scenario.boundaries,
        scenario.intermaterial_resistance_m2_k_w,times,patches,limits.solve,checkpoint),prepared.result.stationary_extrusion};
}

GCodeStudyResult study_gcode_deposition(const GCodeProcessorResult &source,
    const GCodeThermalScenario &scenario, const GCodeStudyLimits &limits,
    const std::function<void(unsigned, bool)> &progress,
    const std::function<void()> &checkpoint)
{
    const std::array<double, 6> absolute{
        limits.force_absolute_tolerance_n, limits.elastic_energy_absolute_tolerance_j,
        limits.heat_absolute_tolerance_j, limits.heat_absolute_tolerance_j,
        limits.heat_absolute_tolerance_j, limits.contact_area_absolute_tolerance_m2};
    if (limits.spatial_levels < 3 || limits.spatial_levels > 4 || limits.max_runs == 0 || limits.max_runs > 8 ||
        !std::isfinite(limits.relative_tolerance) || limits.relative_tolerance <= 0. || limits.relative_tolerance >= 1. ||
        !std::isfinite(limits.run.simulation.thermal.max_step_temperature_error_k) ||
        limits.run.simulation.thermal.max_step_temperature_error_k <= 0. ||
        std::any_of(absolute.begin(), absolute.end(), [](double value) { return !std::isfinite(value) || value < 0.; }))
        throw std::invalid_argument("Invalid thermal-mechanical convergence criteria.");
    using Measures = std::array<double, 6 + ThermalHistory::boundary_kind_count>;
    const auto measures = [](const GCodeSimulationResult &run) -> Measures {
        const auto &thermal = run.simulation.thermal;
        const auto &bed = *run.simulation.mechanical;
        Measures values{bed.tensile_force_n, bed.equilibrium.elastic_energy_j,
            thermal.deposited_energy_j, thermal.final_energy_j, thermal.boundary_energy_j,
            bed.contact_area_total_m2};
        std::copy(thermal.boundary_energy_by_kind_j.begin(), thermal.boundary_energy_by_kind_j.end(), values.begin() + 6);
        return values;
    };
    const auto close = [&](const Measures &a, const Measures &b) {
        for (std::size_t i = 0; i < a.size(); ++i)
            if (!std::isfinite(a[i]) || !std::isfinite(b[i]) ||
                std::abs(a[i] - b[i]) > (i < absolute.size() ? absolute[i] : limits.heat_absolute_tolerance_j) +
                    limits.relative_tolerance * std::max(std::abs(a[i]), std::abs(b[i])))
                return false;
        return true;
    };
    GCodeStudyResult result;
    GCodeThermalScenario refined = scenario;
    Measures previous{};
    unsigned stable = 0;
    for (unsigned level = 0; level < limits.spatial_levels; ++level) {
        if (checkpoint) checkpoint();
        Measures coarse_time{};
        double coarse_interval_s = 0.;
        for (unsigned temporal = 0; temporal < 2; ++temporal) {
            if (result.completed_runs >= limits.max_runs) {
                result.state = GCodeStudyState::BudgetExceeded;
                return result;
            }
            if (progress) progress(level, temporal != 0);
            auto run_limits = limits.run;
            if (temporal) {
                run_limits.simulation.thermal.max_step_s *= .5;
                // Event boundaries or adaptive control may make the configured
                // cap inactive. Refine the intervals actually taken, not just
                // a nominal cap which would reproduce the same integration.
                if (coarse_interval_s > 0.)
                    run_limits.simulation.thermal.max_step_s = std::min(
                        run_limits.simulation.thermal.max_step_s, .5 * coarse_interval_s);
                run_limits.simulation.thermal.max_step_temperature_error_k *= .25;
            }
            GCodeSimulationResult run;
            try {
                run = simulate_gcode_deposition(source, refined, run_limits, checkpoint);
            } catch (const ResourceLimit &) {
                result.state = GCodeStudyState::BudgetExceeded;
                return result;
            }
            ++result.completed_runs;
            if (run.simulation.thermal.state == ThermalHistory::State::BudgetExceeded) {
                result.state = GCodeStudyState::BudgetExceeded;
                return result;
            }
            if (run.simulation.thermal.state != ThermalHistory::State::Converged ||
                !run.simulation.thermal.temporal_error_controlled || !run.simulation.mechanical) {
                result.state = GCodeStudyState::SolverFailure;
                result.failed_thermal_state = run.simulation.thermal.state;
                return result;
            }
            const Measures current = measures(run);
            if (!temporal) {
                coarse_time = current;
                coarse_interval_s = run.simulation.thermal.maximum_accepted_interval_s;
                continue;
            }
            if (!close(coarse_time, current)) return result;
            ++result.completed_levels;
            stable = level != 0 && close(previous, current) ? stable + 1 : 0;
            previous = current;
            if (stable >= 2) {
                result.state = GCodeStudyState::Converged;
                result.qualified = std::move(run);
                return result;
            }
        }
        refined.spacing_m *= .5;
    }
    return result;
}

GCodeHistoryStudy study_gcode_history(const GCodeProcessorResult &source,
    const GCodeThermalScenario &scenario, const std::vector<double> &times_s,
    const GCodeStudyLimits &limits, unsigned max_total_runs,
    const std::function<void(std::size_t, unsigned, bool)> &progress,
    const std::function<void()> &checkpoint)
{
    if (times_s.empty() || times_s.size() > 16 || max_total_runs == 0 || max_total_runs > 128 ||
        scenario.observation_time_s || limits.max_runs == 0 || limits.max_runs > 8)
        throw std::invalid_argument("Invalid G-code history study schedule or budget.");
    for (std::size_t i = 0; i < times_s.size(); ++i)
        if (!std::isfinite(times_s[i]) || times_s[i] <= 0. || (i != 0 && times_s[i] <= times_s[i - 1]))
            throw std::invalid_argument("G-code history observations must be finite and strictly increasing.");
    GCodeHistoryStudy result;
    std::vector<GCodeHistorySample> samples;
    for (std::size_t i = 0; i < times_s.size(); ++i) {
        if (checkpoint) checkpoint();
        result.failed_observation = i;
        if (result.completed_runs >= max_total_runs) {
            result.state = GCodeStudyState::BudgetExceeded;
            return result;
        }
        auto observation = scenario;
        observation.observation_time_s = times_s[i];
        auto observation_limits = limits;
        observation_limits.max_runs = std::min(limits.max_runs, max_total_runs - result.completed_runs);
        // The field is destroyed at the end of this iteration; only demands
        // survive. A later failure must not publish an apparently complete peak.
        const auto study = study_gcode_deposition(source, observation, observation_limits,
            [&](unsigned level, bool finer_time) { if (progress) progress(i, level, finer_time); }, checkpoint);
        result.completed_runs += study.completed_runs;
        if (study.state != GCodeStudyState::Converged || !study.qualified || !study.qualified->simulation.mechanical) {
            result.state = study.state == GCodeStudyState::Converged ? GCodeStudyState::SolverFailure : study.state;
            return result;
        }
        const auto &bed = *study.qualified->simulation.mechanical;
        samples.push_back({times_s[i], bed.tensile_force_n, bed.equilibrium.elastic_energy_j});
        ++result.completed_observations;
    }
    result.peak_force_sample = static_cast<std::size_t>(std::max_element(samples.begin(), samples.end(),
        [](const auto &a, const auto &b) { return a.tensile_force_n < b.tensile_force_n; }) - samples.begin());
    result.peak_energy_sample = static_cast<std::size_t>(std::max_element(samples.begin(), samples.end(),
        [](const auto &a, const auto &b) { return a.elastic_energy_j < b.elastic_energy_j; }) - samples.begin());
    result.qualified_samples = std::move(samples);
    result.failed_observation.reset();
    result.state = GCodeStudyState::Converged;
    return result;
}
} // namespace Slic3r::Thermoelastic
