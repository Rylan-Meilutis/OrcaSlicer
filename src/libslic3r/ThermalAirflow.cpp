#include "ThermalAirflow.hpp"
#include "GCodeReader.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <utility>

namespace Slic3r::ThermalHistory {
bool AirflowState::operator==(const AirflowState &other) const
{
    return pwm == other.pwm && dock_automatic == other.dock_automatic &&
           unassigned_fan_command == other.unassigned_fan_command &&
           unsupported_firmware_command == other.unsupported_firmware_command;
}

AirflowTimeline::AirflowTimeline(std::map<unsigned, FanRole> channels, FirmwareFanPolicy policy,
                                 std::size_t max_events)
    : m_channels(std::move(channels)), m_policy(policy), m_max_events(max_events)
{
    if (max_events == 0) throw std::invalid_argument("Airflow event budget is empty.");
    if (policy != FirmwareFanPolicy::GenericM106 && policy != FirmwareFanPolicy::RmeLegacy &&
        policy != FirmwareFanPolicy::RmeDockCooldown)
        throw std::invalid_argument("Invalid firmware fan policy.");
    if (policy == FirmwareFanPolicy::RmeLegacy) m_state.dock_automatic = false;
    for (const auto &[channel, role] : m_channels)
        if (channel > 255 || role >= FanRole::Count)
            throw std::invalid_argument("Invalid airflow channel mapping.");
}

void AirflowTimeline::consume(const GCodeReader::GCodeLine &line, std::size_t line_id)
{
    if (m_last_line && line_id <= *m_last_line)
        throw std::invalid_argument("Airflow commands must be chronological.");
    auto next = m_state;
    const bool rme = m_policy != FirmwareFanPolicy::GenericM106;
    const auto unknown = [&](FanRole role) { next.pwm[static_cast<unsigned>(role)].reset(); };
    if (line.cmd_is("M106") || line.cmd_is("M107")) {
        float channel = 0.;
        const bool valid_channel = !line.has('P') ||
            (line.has_value('P', channel) && std::isfinite(channel) &&
             channel >= 0. && channel <= 255. && std::floor(channel) == channel);
        const auto found = valid_channel ? m_channels.find(static_cast<unsigned>(channel)) : m_channels.end();
        // Buddy's M106 handler only changes PWM when S, A or R is present.
        // N/G alone changes the chamber regulator, not its manual PWM request.
        const bool no_pwm_request = rme && line.cmd_is("M106") &&
            !line.has('S') && !line.has('A') && !line.has('R');
        if (no_pwm_request) {
            if ((channel == 3. || channel == 4.) && (line.has('N') || line.has('G')))
                unknown(FanRole::Exhaust);
        } else if (found == m_channels.end()) {
            next.unassigned_fan_command = true;
            // An unrecognized channel may alias a known physical fan.
            for (auto &value : next.pwm) value.reset();
        } else {
            float speed = line.cmd_is("M107") ? 0.f : 255.f;
            const bool valid_speed = line.cmd_is("M107") || !line.has('S') || line.has_value('S', speed);
            const bool direct_dock = rme && channel == 6. && found->second == FanRole::Dock;
            const bool direct_part = rme && channel == 0. && found->second == FanRole::Part;
            // M106 A without S keeps the current print-fan request on Buddy.
            const bool retain_part = direct_part && line.cmd_is("M106") && line.has('A') && !line.has('S');
            if (!valid_speed || !std::isfinite(speed) || speed < 0. || speed > 255. ||
                std::floor(speed) != speed ||
                ((!direct_dock && !direct_part) && (line.has('R') || line.has('T') || line.has('N') || line.has('G'))) ||
                (rme && found->second == FanRole::Dock && next.dock_automatic != false) ||
                (rme && found->second == FanRole::Exhaust))
                unknown(found->second);
            else if (!retain_part)
                next.pwm[static_cast<unsigned>(found->second)] = speed;
        }
    } else if (rme) {
        const auto command = line.cmd();
        const bool tool_change = command.size() > 1 && command.front() == 'T' &&
            std::all_of(command.begin() + 1, command.end(), [](char c) { return c >= '0' && c <= '9'; });
        if (line.cmd_is("M979")) {
            if (m_policy == FirmwareFanPolicy::RmeLegacy) {
                next.unsupported_firmware_command = true;
            } else {
                float enabled = -1.;
                if (line.has_value('S', enabled) && (enabled == 0.f || enabled == 1.f))
                    next.dock_automatic = enabled == 1.f;
                else
                    next.dock_automatic.reset();
                // Mode selection alone does not establish actual dock PWM.
                unknown(FanRole::Dock);
            }
        } else if ((tool_change || line.cmd_is("P0")) && next.dock_automatic != false) {
            unknown(FanRole::Dock);
        } else if (line.cmd_is("M141") || line.cmd_is("M191")) {
            unknown(FanRole::Exhaust);
        }
    }
    if (!(next == m_state)) {
        if (m_events.size() >= m_max_events)
            throw std::runtime_error("Airflow event budget exceeded; thermal history is unresolved.");
        m_events.push_back({line_id, next});
        m_state = next;
    }
    m_last_line = line_id;
}
} // namespace Slic3r::ThermalHistory
