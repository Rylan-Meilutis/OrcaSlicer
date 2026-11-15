#pragma once
#include "ThermalAirflow.hpp"
#include "ThermalAirflowTiming.hpp"
#include "ThermalHistory.hpp"
#include <cstddef>
#include <functional>
#include <vector>

namespace Slic3r::Thermoelastic {
struct SectionGCodeScenario;
struct SectionThermalAssembly;
using AirflowExposure = std::function<std::vector<ThermalHistory::Element>(
    double, const SectionThermalAssembly&, const ThermalHistory::AirflowState&)>;

// Attach modeled command timing to an EXPLICIT caller-supplied boundary law.
// Caller must establish that trace and scenario refer to the same final source;
// matching total duration alone does not authenticate this association.
// Replaces exposure, merges changes into exposure_times_s only, and leaves
// deposition activation/materials unchanged. No PWM-to-convection conversion,
// material inference, measured-temperature claim or Auto qualification occurs.
// State is initially unknown, right-continuous (last equal-time event wins),
// and retained through cooldown. Existing boundary times are preserved.
// Invalid input/cancellation leaves the scenario unchanged.
void apply_timed_airflow(SectionGCodeScenario& scenario,
    const ThermalHistory::TimedAirflowTrace& trace, AirflowExposure exposure,
    std::size_t max_events = 100000, const std::function<void()>& checkpoint = {});
} // namespace Slic3r::Thermoelastic
