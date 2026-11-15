#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "libslic3r/PrintConfig.hpp"
#include "test_helpers.hpp"
#include <string>

using namespace Slic3r;

TEST_CASE("Dock-fan commands require explicit firmware capability and preserve manual selection", "[DockFan]")
{
    const bool supported = GENERATE(false, true);
    const bool automatic = GENERATE(false, true);
    const int seconds = GENERATE(180, 300);
    auto config = Test::multifilament_config(1, {
        {"support_rme_indx_dock_fan_cooldown", supported ? "1" : "0"},
        {"rme_indx_dock_fan_cooldown", automatic ? "1" : "0"},
        {"rme_indx_dock_fan_cooldown_time", std::to_string(seconds)},
        {"machine_start_gcode", "; dock-test-start\nG90\nM83\n"},
        {"machine_end_gcode", "M106 P6 S0\n"},
        {"layer_height", "0.2"}, {"initial_layer_print_height", "0.2"},
        {"automatic_warp_mitigation", "0"}
    });
    const auto gcode = Test::slice({Test::cube(1.)}, config);
    const auto command = gcode.find("\nM979 ");
    if (!supported) {
        CHECK(command == std::string::npos);
    } else {
        REQUIRE(command != std::string::npos);
        const std::string expected = automatic ? "M979 S1 C" + std::to_string(seconds) : "M979 S0";
        CHECK(gcode.compare(command + 1, expected.size(), expected) == 0);
        CHECK(command < gcode.find("; dock-test-start"));
        CHECK(gcode.find("\nM979 ", command + 1) == std::string::npos);
    }
    // Firmware owns asynchronous cooldown; no extra dwell or manual fan rewrite.
    CHECK(gcode.find("\nG4 ") == std::string::npos);
    CHECK(gcode.find("\nM106 P6 S0\n") != std::string::npos);
}
