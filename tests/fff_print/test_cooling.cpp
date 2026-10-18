#include <catch2/catch_all.hpp>

#include <catch2/catch_test_macros.hpp>
#include "test_helpers.hpp"
#include "libslic3r/GCode/CoolingBuffer.hpp"

#include <algorithm>
#include <cstddef>
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Config.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/libslic3r.h"
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::Test;

TEST_CASE("Bridge temperature is independent of arcs and restores normal extrusion", "[Cooling][FeatureTemperature]")
{
    const bool fast = GENERATE(false, true);
    const std::string output = slice({TestMesh::bridge}, {
        {"filament_fast_temperature_changes", fast},
        {"filament_bridge_temperature", 201},
        {"filament_arc_temperature", 202},
        {"nozzle_temperature", 215},
        {"arc_overhang_enabled", false}
    });
    const auto target = output.find("M104 S201");
    CHECK((target != std::string::npos) == fast);
    CHECK(output.find("M104 S202") == std::string::npos);
    if (fast)
        CHECK(output.find("M104 S215", target) != std::string::npos);
}

TEST_CASE("Arc fill uses its own temperature rather than the straight bridge target", "[Cooling][FeatureTemperature]")
{
    const std::string output = slice({TestMesh::bridge}, {
        {"filament_fast_temperature_changes", true},
        {"filament_bridge_temperature", 201},
        {"filament_arc_temperature", 202},
        {"nozzle_temperature", 215},
        {"arc_overhang_enabled", true},
        {"arc_overhang_bridges", true},
        {"arc_overhang_bridge_distance", 0.}
    });
    const auto target = output.find("M104 S202");
    REQUIRE(target != std::string::npos);
    CHECK(output.find("M104 S215", target) != std::string::npos);
}

TEST_CASE("Calibration ignores feature temperature overrides", "[Cooling][FeatureTemperature]")
{
    Print print;
    init_and_process_print({TestMesh::bridge}, print, {
        {"filament_fast_temperature_changes", true},
        {"filament_bridge_temperature", 201},
        {"filament_arc_temperature", 202},
        {"filament_overhang_temperature", 203},
        {"arc_overhang_enabled", false}
    });
    print.calib_mode() = CalibMode::Calib_Temp_Tower;
    const auto output = gcode(print);
    CHECK(output.find("M104 S201") == std::string::npos);
    CHECK(output.find("M104 S202") == std::string::npos);
    CHECK(output.find("M104 S203") == std::string::npos);
}

TEST_CASE("Steep wall temperature does not require overhang fan or speed controls", "[Cooling][FeatureTemperature]")
{
    const std::string output = slice({TestMesh::overhang}, {
        {"filament_fast_temperature_changes", true},
        {"filament_overhang_temperature", 203},
        {"nozzle_temperature", 215},
        {"enable_overhang_speed", false},
        {"enable_overhang_bridge_fan", false},
        {"arc_overhang_enabled", false}
    });
    const auto target = output.find("M104 S203");
    REQUIRE(target != std::string::npos);
    CHECK(output.find("M104 S215", target) != std::string::npos);
}

// The fan is held off for the first close_fan_the_first_x_layers layers, so an explicit
// fan-off command is emitted.
TEST_CASE("Fan is held off for the initial layers", "[Cooling]")
{
    const std::string gcode = slice({ cube(20) }, {
        { "cooling",                      true },
        { "close_fan_the_first_x_layers", 5 },
    });
    CHECK(gcode.find("M106 S0") != std::string::npos);
}

// The cooling pass resolves and strips its internal speed placeholders; none leak into
// the final G-code.
TEST_CASE("Cooling consumes its internal speed markers", "[Cooling]")
{
    const std::string gcode = slice({ cube(20) }, { { "layer_height", 0.2 } });
    CHECK(gcode.find(";_EXTRUDE_SET_SPEED") == std::string::npos);
}

TEST_CASE("Hull-line mitigation limits abrupt layer-time reductions", "[Cooling][HullLine]")
{
    CHECK_THAT(CoolingBuffer::hull_line_target_time(40.f, 20.f, 25.f),
               Catch::Matchers::WithinAbs(32.f, EPSILON));
    CHECK_THAT(CoolingBuffer::hull_line_target_time(40.f, 35.f, 25.f),
               Catch::Matchers::WithinAbs(35.f, EPSILON));
    CHECK_THAT(CoolingBuffer::hull_line_target_time(0.f, 20.f, 25.f),
               Catch::Matchers::WithinAbs(20.f, EPSILON));
}

TEST_CASE("Arc fitting remains active with overhang fan control", "[Cooling][Regression]")
{
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({
        { "enable_arc_fitting",             true },
        { "enable_overhang_bridge_fan",     true },
        { "enable_overhang_speed",          true },
        { "initial_layer_print_height",     0.2 },
        { "layer_height",                   0.2 },
        { "slow_down_for_layer_cooling",    false },
    });
    config.set_key_value("fan_max_speed", new ConfigOptionFloats{20.0});
    config.set_key_value("fan_min_speed", new ConfigOptionFloats{20.0});
    config.set_key_value("overhang_fan_speed", new ConfigOptionInts{100});
    config.set_key_value("overhang_fan_threshold", new ConfigOptionEnumsGeneric{Overhang_threshold_2_4});
    config.set_key_value("layer_change_gcode", new ConfigOptionString{";TEST_LAYER=[layer_num]"});

    const auto has_arc_after_first_layer = [](const std::string &gcode) {
        std::istringstream input(gcode);
        std::string line;
        size_t layers_seen = 0;
        while (std::getline(input, line)) {
            if (line.rfind(";TEST_LAYER=", 0) == 0)
                ++layers_seen;
            else if (layers_seen > 1 && (line.rfind("G2 ", 0) == 0 || line.rfind("G3 ", 0) == 0))
                return true;
        }
        return false;
    };
    const auto has_arc_with_overhang_fan = [](const std::string &gcode) {
        std::istringstream input(gcode);
        std::string line;
        bool overhang_fan_active = false;
        bool outer_wall = false;
        while (std::getline(input, line)) {
            if (line.rfind("M106 S255", 0) == 0)
                overhang_fan_active = true;
            else if (line.rfind("M106 S0", 0) == 0 || line.rfind("M107", 0) == 0)
                overhang_fan_active = false;
            else if (line.rfind(";TYPE:", 0) == 0)
                outer_wall = line == ";TYPE:Outer wall";
            else if (overhang_fan_active && outer_wall && (line.rfind("G2 ", 0) == 0 || line.rfind("G3 ", 0) == 0))
                return true;
        }
        return false;
    };

    const std::string with_overhang_speed = slice({ make_cylinder(25.0, 10.0, 2.0 * PI / 72.0) }, config);
    CHECK(has_arc_after_first_layer(with_overhang_speed));

    config.set_deserialize_strict({ { "enable_overhang_speed", false } });
    const std::string without_overhang_speed = slice({ make_cylinder(25.0, 10.0, 2.0 * PI / 72.0) }, config);
    CHECK(has_arc_after_first_layer(without_overhang_speed));

    const std::string uniform_overhang_gcode = slice({ make_sphere(50.0, 2.0 * PI / 72.0) }, config);
    CHECK(has_arc_with_overhang_fan(uniform_overhang_gcode));
}

TEST_CASE("Overhang fan transitions do not depend on overhang speed", "[Cooling][Regression]")
{
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({
        { "bridge_speed",                   2.0 },
        { "enable_arc_fitting",             false },
        { "enable_overhang_bridge_fan",     true },
        { "enable_overhang_speed",          false },
        { "initial_layer_print_height",     0.3 },
        { "inner_wall_speed",              30.0 },
        { "layer_height",                   0.3 },
        { "outer_wall_speed",              30.0 },
        { "overhang_1_4_speed",            "30" },
        { "overhang_2_4_speed",            "29" },
        { "overhang_3_4_speed",             "6" },
        { "overhang_4_4_speed",             "3" },
        { "slow_down_for_layer_cooling",    false },
        { "slowdown_for_curled_perimeters", false },
    });
    config.set_key_value("fan_max_speed", new ConfigOptionFloats{20.0});
    config.set_key_value("fan_min_speed", new ConfigOptionFloats{20.0});
    config.set_key_value("overhang_fan_speed", new ConfigOptionInts{100});
    config.set_key_value("overhang_fan_threshold", new ConfigOptionEnumsGeneric{Overhang_threshold_2_4});
    config.set_key_value("layer_change_gcode", new ConfigOptionString{";TEST_LAYER_Z=[layer_z]"});

    const auto fan_commands = [](const std::string &gcode) {
        std::vector<std::pair<std::string, std::string>> commands;
        std::istringstream input(gcode);
        std::string layer;
        std::string line;
        while (std::getline(input, line)) {
            if (line.rfind(";TEST_LAYER_Z=", 0) == 0)
                layer = line;
            else if (!layer.empty() && (line.rfind("M106", 0) == 0 || line.rfind("M107", 0) == 0))
                commands.emplace_back(layer, line);
        }
        return commands;
    };
    const auto feedrates = [](const std::string &gcode) {
        std::vector<std::string> values;
        std::istringstream input(gcode);
        std::string word;
        while (input >> word)
            if (!word.empty() && word.front() == 'F')
                values.push_back(word);
        return values;
    };

    constexpr double sphere_radius = 50.0; // 100 mm diameter.
    const std::string without_speed_gcode = slice({make_sphere(sphere_radius, PI / 24.0)}, config);
    config.set_deserialize_strict({{"enable_overhang_speed", true}});
    const std::string with_speed_gcode = slice({make_sphere(sphere_radius, PI / 24.0)}, config);

    const auto without_speed_fan = fan_commands(without_speed_gcode);
    const auto with_speed_fan = fan_commands(with_speed_gcode);
    const auto without_speed_feedrates = feedrates(without_speed_gcode);
    const auto with_speed_feedrates = feedrates(with_speed_gcode);

    REQUIRE_FALSE(without_speed_fan.empty());
    REQUIRE(std::any_of(without_speed_fan.begin(), without_speed_fan.end(),
                        [](const auto &command) { return command.second.find("S255") != std::string::npos; }));
    REQUIRE(with_speed_feedrates != without_speed_feedrates);
    CHECK(with_speed_fan == without_speed_fan);
}
