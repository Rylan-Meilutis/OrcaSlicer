#pragma once
#include "SectionGCodeHistory.hpp"
#include "SectionThermalAssembly.hpp"
#include "ThermalHistory.hpp"
#include "ThermoelasticToolpath.hpp"
#include "ThermoelasticConflictOwnership.hpp"
#include "ThermoelasticExactSections.hpp"
#include <array>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace Slic3r::Thermoelastic {
struct SectionCoolingCase {
    // Identifies an explicit assumed boundary law, not measured fan behavior.
    std::string label;
    std::vector<double> exposure_times_s;
    std::function<std::vector<ThermalHistory::Element>(double, const SectionThermalAssembly&)> exposure;
};
struct SectionCoolingMetrics {
    std::string label;
    // Positive tensile generalized nodal reactions, not local peel traction
    // or an independently integrated adhesive-interface force.
    double tensile_bed_demand_n = 0., elastic_energy_j = 0.;
    double deposited_energy_j = 0., final_energy_j = 0., boundary_energy_j = 0., energy_balance_error_j = 0.;
    double occupied_volume_m3 = 0., excluded_nominal_volume_m3 = 0., modeled_print_time_s = 0.;
    double modeled_observation_time_s = 0.; // Print duration plus declared cooldown.
    StationaryExtrusionSummary stationary_extrusion;
    CollapsedEnvelopeSummary collapsed_envelope;
    std::size_t excluded_segments = 0, activated_snapshots = 0, exposure_changes = 0;
    // Retain the common material-location assumptions, not just inventories.
    bool uses_fixed_interface_reference = false;
    bool allows_single_sided_native_boundary = false, allows_dual_native_boundary = false;
    NativeInterfaceVariant native_interface_variant = NativeInterfaceVariant::Strict;
    DepositedOverlapPolicy deposited_overlap_policy = DepositedOverlapPolicy::None;
};
// Exactly two sequential simulations of the SAME source and common scenario.
// Only exposure callbacks/schedules are replaced. Common callbacks must be
// deterministic and side-effect free. Per-case limits imply at most twice
// their work budget; only one full simulation field is retained at a time.
// Returns final-time diagnostics only: not peak demand, an improvement verdict,
// convergence study, source authentication, or an adhesion safety assessment.
// Exceptions/cancellation/nonconvergence publish no partial pair.
std::array<SectionCoolingMetrics, 2> compare_section_cooling(
    const GCodeProcessorResult& source, const SectionGCodeScenario& common,
    const std::array<SectionCoolingCase, 2>& cases, const SectionGCodeLimits& per_case_limits = {},
    const std::function<void()>& checkpoint = {});
} // namespace Slic3r::Thermoelastic
