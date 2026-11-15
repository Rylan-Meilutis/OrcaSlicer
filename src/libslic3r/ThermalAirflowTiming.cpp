#include "ThermalAirflowTiming.hpp"
#include "ThermalAirflow.hpp"

#include "GCode/GCodeProcessor.hpp"
#include "PrintConfig.hpp"
#include "ThermoelasticToolpath.hpp"
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

namespace Slic3r::ThermalHistory {
TimedAirflowTrace time_airflow_events(const GCodeProcessorResult& source,
    const std::vector<AirflowEvent>& events, const AirflowTimingLimits& limits,
    const std::function<void()>& checkpoint)
{
    if (!limits.max_moves || limits.max_moves > 2000000 || !limits.max_segments ||
        limits.max_segments > 2000000 || !limits.max_events || limits.max_events > 100000 ||
        events.size() > limits.max_events)
        throw std::invalid_argument("Invalid airflow timing budget.");
    for (std::size_t i = 0; i < events.size(); ++i) {
        if (checkpoint && i % 1024 == 0) checkpoint();
        if (!events[i].line_id || (i && events[i].line_id < events[i - 1].line_id))
            throw std::invalid_argument("Airflow events require ordered one-based line IDs.");
        for (const auto pwm : events[i].state.pwm)
            if (pwm && (!std::isfinite(*pwm) || *pwm < 0. || *pwm > 255.))
                throw std::invalid_argument("Invalid commanded airflow PWM.");
    }
    TimedAirflowTrace result;
    {
        // Release the temporary geometry before allocating the event output.
        const auto timing = Thermoelastic::timed_deposition(source, limits.max_moves, limits.max_segments, checkpoint);
        if (timing.state != Thermoelastic::ToolpathState::Complete)
            throw std::invalid_argument("Airflow timing requires complete modeled move timing.");
        result.total_s = timing.total_s;
    }
    result.events.reserve(events.size());
    constexpr auto mode = static_cast<std::size_t>(PrintEstimatedStatistics::ETimeMode::Normal);
    std::size_t event = 0, previous_line = 0;
    double elapsed = 0.;
    for (std::size_t i = 0; i < source.moves.size(); ++i) {
        if (checkpoint && i % 1024 == 0) checkpoint();
        const auto& move = source.moves[i];
        if (move.internal_only) continue;
        if (move.gcode_id < previous_line)
            throw std::invalid_argument("Physical move line IDs are not ordered.");
        previous_line = move.gcode_id;
        while (event < events.size() && events[event].line_id < move.gcode_id) {
            if (checkpoint && event % 1024 == 0) checkpoint();
            result.events.push_back({events[event].line_id, elapsed, events[event].state});
            ++event;
        }
        elapsed += move.time[mode];
    }
    while (event < events.size()) {
        if (checkpoint && event % 1024 == 0) checkpoint();
        result.events.push_back({events[event].line_id, elapsed, events[event].state});
        ++event;
    }
    return result;
}
} // namespace Slic3r::ThermalHistory
