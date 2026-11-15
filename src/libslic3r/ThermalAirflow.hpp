#pragma once

#include "GCodeReader.hpp"
#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <optional>
#include <vector>

namespace Slic3r::ThermalHistory {

enum class FanRole : unsigned { Part, Auxiliary, Dock, Exhaust, Count };
// Select from known installed firmware capabilities, not the printer name alone.
enum class FirmwareFanPolicy { GenericM106, RmeLegacy, RmeDockCooldown };
struct AirflowState {
    // Commanded PWM only, never measured airflow or a heat-transfer coefficient.
    // Missing values require a boundary-condition scenario, not zero cooling.
    std::array<std::optional<double>, static_cast<unsigned>(FanRole::Count)> pwm;
    std::optional<bool> dock_automatic; // Unknown until explicit mode selection.
    bool unassigned_fan_command = false;
    bool unsupported_firmware_command = false;
    bool operator==(const AirflowState &other) const;
};
struct AirflowEvent {
    std::size_t line_id;
    AirflowState state;
};

// Explicit printer-specific P-channel mapping, including aliases if appropriate.
// No universal P1/P2/P3 meaning is assumed. This decoder handles M106/M107;
// macros and closed-loop controllers require runtime uncertainty. The RME
// policies distinguish direct dock PWM from the optional M979 cooldown protocol.
class AirflowTimeline {
public:
    explicit AirflowTimeline(std::map<unsigned, FanRole> channels,
                             FirmwareFanPolicy policy, std::size_t max_events = 100000);
    void consume(const GCodeReader::GCodeLine &line, std::size_t line_id);
    const AirflowState &state() const { return m_state; }
    const std::vector<AirflowEvent> &events() const { return m_events; }
private:
    std::map<unsigned, FanRole> m_channels;
    FirmwareFanPolicy m_policy;
    std::size_t m_max_events;
    std::optional<std::size_t> m_last_line;
    AirflowState m_state;
    std::vector<AirflowEvent> m_events;
};

} // namespace Slic3r::ThermalHistory
