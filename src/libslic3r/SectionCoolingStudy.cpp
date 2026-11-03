#include "SectionCoolingStudy.hpp"
#include "GCode/GCodeProcessor.hpp"
#include "SectionGCodeHistory.hpp"
#include "ThermalHistory.hpp"
#include "ThermoelasticConflictOwnership.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <stdexcept>

namespace Slic3r::Thermoelastic {
std::array<SectionCoolingMetrics, 2> compare_section_cooling(
    const GCodeProcessorResult& source, const SectionGCodeScenario& common,
    const std::array<SectionCoolingCase, 2>& cases, const SectionGCodeLimits& limits,
    const std::function<void()>& checkpoint)
{
    if (checkpoint) checkpoint();
    if (limits.ownership.interface_variant==NativeInterfaceVariant::Strict &&
        (limits.ownership.allow_single_sided_native_boundary || limits.ownership.allow_dual_native_boundary))
        throw std::invalid_argument("Section cooling boundary sensitivity requires an explicit non-strict ownership policy.");
    for (const auto& study : cases) {
        if (study.label.empty() || study.label.size() > 1024 || !study.exposure ||
            common.activation_times_s.empty() || common.activation_times_s.size() > 10000 ||
            study.exposure_times_s.size() > 10000-common.activation_times_s.size())
            throw std::invalid_argument("Invalid section cooling case.");
        double previous = -1.;
        for (double time : study.exposure_times_s) {
            if (!std::isfinite(time) || time < 0. || time <= previous ||
                time > common.activation_times_s.back()+common.cooldown_s)
                throw std::invalid_argument("Invalid section cooling exposure schedule.");
            previous = time;
        }
    }
    if (cases[0].label == cases[1].label)
        throw std::invalid_argument("Section cooling cases need distinct labels.");
    std::array<SectionCoolingMetrics, 2> output;
    for (std::size_t i = 0; i < cases.size(); ++i) {
        if (checkpoint) checkpoint();
        auto scenario = common;
        scenario.exposure = cases[i].exposure;
        scenario.exposure_times_s = cases[i].exposure_times_s;
        const auto history = simulate_section_gcode_history(source, scenario, limits, checkpoint);
        if (history.thermal.state != ThermalHistory::State::Converged || !history.mechanical)
            throw std::invalid_argument("Section cooling comparison requires complete converged solves.");
        auto& result = output[i];
        result.label = cases[i].label;
        for (const auto& force : history.mechanical->generalized_bed_reaction_n) {
            if (!force.allFinite()) throw std::invalid_argument("Nonfinite section cooling reaction.");
            result.tensile_bed_demand_n += std::max(0., -force.z());
        }
        result.elastic_energy_j = history.mechanical->elastic_energy_j;
        result.deposited_energy_j = history.thermal.deposited_energy_j;
        result.final_energy_j = history.thermal.final_energy_j;
        result.boundary_energy_j = history.thermal.boundary_energy_j;
        result.energy_balance_error_j = history.thermal.energy_balance_error_j;
        result.occupied_volume_m3 = history.occupied_volume_m3;
        result.excluded_nominal_volume_m3 = history.excluded_nominal_volume_m3;
        result.modeled_print_time_s = history.modeled_print_time_s;
        result.modeled_observation_time_s = history.modeled_print_time_s + common.cooldown_s;
        result.stationary_extrusion = history.stationary_extrusion;
        result.collapsed_envelope = history.collapsed_envelope;
        result.excluded_segments = history.excluded_segments;
        result.activated_snapshots = history.activated_snapshots;
        result.exposure_changes = history.exposure_changes;
        result.uses_fixed_interface_reference = history.uses_fixed_interface_reference;
        result.allows_single_sided_native_boundary = history.allows_single_sided_native_boundary;
        result.allows_dual_native_boundary = history.allows_dual_native_boundary;
        result.native_interface_variant = history.native_interface_variant;
        result.deposited_overlap_policy = history.deposited_overlap_policy;
        for (const double value : {result.tensile_bed_demand_n, result.elastic_energy_j, result.deposited_energy_j,
                 result.final_energy_j, result.boundary_energy_j, result.energy_balance_error_j,
                 result.occupied_volume_m3, result.excluded_nominal_volume_m3, result.modeled_print_time_s,
                 result.modeled_observation_time_s,
                 result.stationary_extrusion.commanded_volume_m3, result.collapsed_envelope.commanded_volume_m3})
            if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite section cooling metric.");
    }
    const auto& a = output[0];
    const auto& b = output[1];
    if (a.uses_fixed_interface_reference!=b.uses_fixed_interface_reference ||
        a.allows_single_sided_native_boundary!=b.allows_single_sided_native_boundary ||
        a.allows_dual_native_boundary!=b.allows_dual_native_boundary ||
        a.native_interface_variant!=b.native_interface_variant ||
        a.deposited_overlap_policy!=b.deposited_overlap_policy)
        throw std::invalid_argument("Section cooling cases changed their material ownership policy.");
    // Same source and deterministic geometry/deposition callbacks must yield
    // identical inventories. Never present omission or source changes as cooling.
    if (a.occupied_volume_m3 != b.occupied_volume_m3 || a.excluded_nominal_volume_m3 != b.excluded_nominal_volume_m3 ||
        a.deposited_energy_j != b.deposited_energy_j || a.modeled_print_time_s != b.modeled_print_time_s ||
        a.excluded_segments != b.excluded_segments || a.activated_snapshots != b.activated_snapshots ||
        a.stationary_extrusion.count != b.stationary_extrusion.count ||
        a.stationary_extrusion.commanded_volume_m3 != b.stationary_extrusion.commanded_volume_m3 ||
        a.collapsed_envelope.count != b.collapsed_envelope.count ||
        a.collapsed_envelope.commanded_volume_m3 != b.collapsed_envelope.commanded_volume_m3)
        throw std::invalid_argument("Section cooling cases changed geometry, deposition heat or omitted material.");
    if (checkpoint) checkpoint();
    return output;
}
} // namespace Slic3r::Thermoelastic
