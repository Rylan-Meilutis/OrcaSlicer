#include "SectionAirflowScenario.hpp"
#include "SectionGCodeHistory.hpp"
#include "SectionThermalAssembly.hpp"
#include "ThermalAirflow.hpp"
#include "ThermalAirflowTiming.hpp"
#include "ThermalHistory.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iterator>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
void apply_timed_airflow(SectionGCodeScenario& scenario,
    const ThermalHistory::TimedAirflowTrace& trace, AirflowExposure exposure,
    std::size_t max_events, const std::function<void()>& checkpoint)
{
    if (checkpoint) checkpoint();
    if (!exposure || !max_events || max_events > 100000 || trace.events.size() > max_events ||
        !trace.uses_modeled_timing || !scenario.acknowledge_modeled_timing ||
        !std::isfinite(trace.total_s) || trace.total_s <= 0. ||
        !std::isfinite(scenario.cooldown_s) || scenario.cooldown_s < 0. || scenario.cooldown_s > 86400. ||
        scenario.activation_times_s.empty() || scenario.activation_times_s.size() > 10000 ||
        scenario.exposure_times_s.size() > 10000 - scenario.activation_times_s.size() ||
        scenario.activation_times_s.back() != trace.total_s)
        throw std::invalid_argument("Invalid explicit airflow scenario timing or budget.");
    const double end = trace.total_s + scenario.cooldown_s;
    if (!std::isfinite(end)) throw std::invalid_argument("Invalid airflow scenario end time.");
    double previous = 0.;
    for (double time : scenario.activation_times_s) {
        if (!std::isfinite(time) || time <= previous)
            throw std::invalid_argument("Invalid airflow scenario activation schedule.");
        previous = time;
    }
    previous = -1.;
    for (double time : scenario.exposure_times_s) {
        if (!std::isfinite(time) || time < 0. || time <= previous || time > end)
            throw std::invalid_argument("Invalid airflow scenario boundary schedule.");
        previous = time;
    }
    auto times = scenario.exposure_times_s;
    std::vector<double> changes;
    changes.reserve(trace.events.size());
    for (std::size_t i = 0; i < trace.events.size(); ++i) {
        if (checkpoint && i % 1024 == 0) checkpoint();
        const auto& event = trace.events[i];
        if (!event.line_id || !std::isfinite(event.time_s) || event.time_s < 0. || event.time_s > trace.total_s ||
            (i && (event.time_s < trace.events[i - 1].time_s || event.line_id < trace.events[i - 1].line_id)))
            throw std::invalid_argument("Invalid ordered modeled airflow events.");
        for (const auto pwm : event.state.pwm)
            if (pwm && (!std::isfinite(*pwm) || *pwm < 0. || *pwm > 255.))
                throw std::invalid_argument("Invalid explicit airflow command state.");
        if (changes.empty() || changes.back() != event.time_s) changes.push_back(event.time_s);
    }
    times.insert(times.end(), changes.begin(), changes.end());
    std::sort(times.begin(), times.end());
    times.erase(std::unique(times.begin(), times.end()), times.end());
    if (times.size() > 10000 - scenario.activation_times_s.size())
        throw std::invalid_argument("Explicit airflow boundary schedule exceeds its budget.");
    // Capture an independent immutable copy. The callback may outlive the trace.
    std::function<std::vector<ThermalHistory::Element>(double, const SectionThermalAssembly&)> boundary =
        [events = trace.events, law = std::move(exposure), end](double time, const SectionThermalAssembly& assembly) {
            if (!std::isfinite(time) || time < 0. || time > end)
                throw std::invalid_argument("Airflow exposure query is outside its scenario.");
            const auto after = std::upper_bound(events.begin(), events.end(), time,
                [](double t, const ThermalHistory::TimedAirflowEvent& event) { return t < event.time_s; });
            const ThermalHistory::AirflowState state = after == events.begin() ?
                ThermalHistory::AirflowState{} : std::prev(after)->state;
            return law(time, assembly, state);
        };
    if (checkpoint) checkpoint();
    scenario.exposure = std::move(boundary);
    scenario.exposure_times_s = std::move(times);
}
} // namespace Slic3r::Thermoelastic
