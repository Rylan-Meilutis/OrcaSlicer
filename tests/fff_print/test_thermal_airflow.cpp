#include <catch2/catch_test_macros.hpp>
#include "libslic3r/ThermalAirflow.hpp"
#include "libslic3r/GCodeReader.hpp"
#include <cstddef>
#include <stdexcept>
#include <string>

using namespace Slic3r;
using namespace Slic3r::ThermalHistory;

TEST_CASE("Thermal airflow keeps fan channels separate and runtime dock cooling unresolved", "[ThermalAirflow]")
{
    AirflowTimeline trace({{0,FanRole::Part}, {2,FanRole::Auxiliary}, {3,FanRole::Exhaust}, {6,FanRole::Dock}}, FirmwareFanPolicy::RmeDockCooldown);
    GCodeReader reader;
    std::size_t id = 0;
    const auto consume = [&](const std::string &text) {
        reader.parse_line(text, [&](auto &, const auto &line) { trace.consume(line, ++id); });
    };
    for (const auto &pwm : trace.state().pwm) CHECK_FALSE(pwm.has_value());
    consume("M106 S200");
    consume("M106 P2 S64");
    consume("M106 P3 S90");
    CHECK_FALSE(trace.state().pwm[3].has_value()); // Chamber/filter controller may override it.
    consume("M106 P6 S255");
    CHECK_FALSE(trace.state().dock_automatic.has_value());
    CHECK_FALSE(trace.state().pwm[2].has_value());
    CHECK(trace.state().pwm[0].value() == 200.);
    consume("M107 P2");
    CHECK(trace.state().pwm[0].value() == 200.);
    CHECK(trace.state().pwm[1].value() == 0.);
    consume("M979 S1 C180");
    CHECK(trace.state().dock_automatic.value());
    CHECK_FALSE(trace.state().pwm[2].has_value());
    consume("M106 P6 S0");
    CHECK_FALSE(trace.state().pwm[2].has_value());
    consume("M979 S0");
    consume("M106 P6 S255");
    CHECK(trace.state().pwm[2].has_value());
    consume("T2");
    CHECK(trace.state().pwm[2].value() == 255.); // Explicit manual mode survives tool changes.
    consume("M141 S50");
    CHECK_FALSE(trace.state().pwm[3].has_value());
    consume("M106 P6 S255");
    consume("P0");
    CHECK(trace.state().pwm[2].value() == 255.);
    consume("M979 S1 C180");
    consume("T3");
    CHECK_FALSE(trace.state().pwm[2].has_value());
    consume("M106 P8 S255");
    CHECK(trace.state().unassigned_fan_command);
    for (const auto &pwm : trace.state().pwm) CHECK_FALSE(pwm.has_value());
}

TEST_CASE("Thermal airflow applies only declared aliases and bounds event storage", "[ThermalAirflow]")
{
    AirflowTimeline trace({{0,FanRole::Part},{1,FanRole::Part}}, FirmwareFanPolicy::GenericM106, 2);
    GCodeReader reader;
    std::size_t id = 0;
    const auto consume = [&](const std::string &text) {
        reader.parse_line(text, [&](auto &, const auto &line) { trace.consume(line, ++id); });
    };
    consume("M106 P1 S100 ; M106 P0 S0");
    CHECK(trace.state().pwm[0].value() == 100.);
    consume("M106 S100"); // Same state does not grow the event list.
    CHECK(trace.events().size() == 1);
    consume("M106 S256");
    CHECK_FALSE(trace.state().pwm[0].has_value());
    CHECK_THROWS_AS(consume("M107"), std::runtime_error);
    CHECK_FALSE(trace.state().pwm[0].has_value());
}

TEST_CASE("Legacy RME dock PWM is direct and unsupported cooldown cannot imply automation", "[ThermalAirflow]")
{
    AirflowTimeline trace({{0,FanRole::Part},{3,FanRole::Exhaust},{4,FanRole::Exhaust},{6,FanRole::Dock}}, FirmwareFanPolicy::RmeLegacy);
    GCodeReader reader;
    std::size_t id = 0;
    const auto consume = [&](const std::string &text) {
        reader.parse_line(text, [&](auto &, const auto &line) { trace.consume(line, ++id); });
    };
    CHECK(trace.state().dock_automatic == false);
    CHECK_FALSE(trace.state().pwm[2].has_value());
    consume("M106 P6 S100");
    CHECK(trace.state().pwm[2].value() == 100.);
    consume("M106 P6"); // No S/A/R means no action on Buddy, unlike generic M106.
    consume("T4");
    consume("P0");
    CHECK(trace.state().pwm[2].value() == 100.);
    consume("M106 P6 R T3"); // Legacy P6 ignores R/T and defaults S to 255.
    CHECK(trace.state().pwm[2].value() == 255.);
    consume("M979 S1 C180");
    CHECK(trace.state().unsupported_firmware_command);
    CHECK(trace.state().dock_automatic == false);
    CHECK(trace.state().pwm[2].value() == 255.);
    consume("M107 P6");
    CHECK(trace.state().pwm[2].value() == 0.);
    consume("M106 S123");
    consume("M106");
    consume("M106 A");
    CHECK(trace.state().pwm[0].value() == 123.);
    consume("M106 P3 N40 G10");
    consume("M106 P4 S0");
    consume("M191 S55");
    CHECK_FALSE(trace.state().pwm[3].has_value());
    consume("M106 P6 S1.5"); // Reject malformed byte commands instead of inventing rounding.
    CHECK_FALSE(trace.state().pwm[2].has_value());
}
