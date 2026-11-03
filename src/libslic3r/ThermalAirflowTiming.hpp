#pragma once

#include "ThermalAirflow.hpp"
#include <cstddef>
#include <functional>
#include <vector>

namespace Slic3r {
struct GCodeProcessorResult;
namespace ThermalHistory {

struct TimedAirflowEvent {
    std::size_t line_id;
    double time_s;
    AirflowState state;
};
struct AirflowTimingLimits {
    std::size_t max_moves = 2000000, max_segments = 250000, max_events = 100000;
};
struct TimedAirflowTrace {
    std::vector<TimedAirflowEvent> events;
    double total_s = 0.;
    bool uses_modeled_timing = true;
};

// Caller MUST supply events from the same final G-code stream as the processor
// result, using its one-based line IDs. This API cannot authenticate that source
// association; filenames or coincident line counts do not establish it.
//
// An event is placed after all non-internal records attributed to its line.
// Thus a synthetic event sharing an arc's line follows all its children;
// ordinary fan commands occupy their own lines. Equal-line events retain input
// order. Before-first/after-last events map to zero/total respectively.
//
// The complete timed_deposition admission is required. Travel and attributed
// waits/tool changes (including legacy seam timing) contribute to the total;
// display-only internal vertices do not. These are MODELED line boundaries,
// not measured firmware execution times or resolved heater/deposition phases.
// Unknown/automatic airflow states are preserved, never converted into h or
// assumed to mean zero cooling. No Auto/physical qualification is provided.
TimedAirflowTrace time_airflow_events(const GCodeProcessorResult& source,
    const std::vector<AirflowEvent>& events, const AirflowTimingLimits& limits = {},
    const std::function<void()>& checkpoint = {});

} // namespace ThermalHistory
} // namespace Slic3r
