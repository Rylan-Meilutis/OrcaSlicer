#include <catch2/catch_all.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/Print.hpp"
#include "libslic3r/calib.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include "test_helpers.hpp"
#include "libslic3r/GCode/CoolingBuffer.hpp"
#include "libslic3r/GCode/FanMover.hpp"
#include "libslic3r/GCode.hpp"
#include "libslic3r/GCodeWriter.hpp"

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

TEST_CASE("Custom fan changes restore the active cooling request", "[Cooling][FanResume]")
{
    const bool overhang = GENERATE(false, true);
    const int normal_fan = GENERATE(0, 100);
    PrintConfig config;
    config.fan_min_speed.values = {double(normal_fan)};
    config.fan_max_speed.values = {double(normal_fan)};
    config.reduce_fan_stop_start_freq.values = {true};
    config.close_fan_the_first_x_layers.values = {0};
    config.initial_layer_fan_speed.values = {-1};
    config.full_fan_speed_layer.values = {0};
    config.overhang_fan_speed.values = {80};
    config.enable_overhang_bridge_fan.values = {true};
    config.slow_down_for_layer_cooling.values = {false};
    GCode generator;
    generator.apply_print_config(config);
    generator.writer().set_extruders({0});
    CoolingBuffer cooling(generator);
    std::string input = overhang ? ";_OVERHANG_FAN_START\n" : "";
    input += "M106 S200\n;_FORCE_RESUME_FAN_SPEED\nG1 X1 F600\n";
    const std::string output = cooling.process_layer(std::move(input), 2, true);
    const auto custom = output.find("M106 S200");
    REQUIRE(custom != std::string::npos);
    const std::string expected = normal_fan == 100 ? "M106 S255" : overhang ? "M106 S204" : "M106 S0";
    const auto restored = output.find(expected, custom);
    REQUIRE(restored != std::string::npos);
    CHECK(restored < output.find("G1 X1"));
    CHECK(output.find("_FORCE_RESUME_FAN_SPEED") == std::string::npos);
}

TEST_CASE("Tool boundaries restore incoming cooling after selecting the tool", "[Cooling][FanResume][Toolchange]")
{
    const bool first_layer = GENERATE(false, true);
    const bool custom_fan = GENERATE(false, true);
    const bool same_fan = GENERATE(false, true);
    const bool continuing_feature = GENERATE(false, true);
    CAPTURE(first_layer, custom_fan, same_fan, continuing_feature);
    const auto dynamic = multifilament_config(2, {
        {"filament_type", "ASA;TPU"},
        {"fan_min_speed", same_fan ? "0,0" : "0,40"},
        {"fan_max_speed", same_fan ? "0,0" : "0,40"},
        {"reduce_fan_stop_start_freq", "1,1"},
        {"close_fan_the_first_x_layers", "1,1"},
        {"initial_layer_fan_speed", "-1,-1"}, {"full_fan_speed_layer", "0,0"},
        {"overhang_fan_speed", "80,80"}, {"enable_overhang_bridge_fan", "0,1"},
        {"slow_down_for_layer_cooling", "0,0"}});
    PrintConfig config;
    config.apply(dynamic, true);
    GCode generator;
    generator.apply_print_config(config);
    generator.writer().set_extruders({0, 1});
    CoolingBuffer cooling(generator);
    const std::string cleaning = custom_fan ?
        ";_FAN_MOVER_SKIP_START\nM106 S200\n;_FAN_MOVER_SKIP_END\n;_FORCE_RESUME_FAN_SPEED\n" : "";
    // Consecutive feature paths may share a START across tools. The incoming
    // profile decides whether it applies; END always clears the shared scope.
    const std::string end = continuing_feature ? "" : ";_OVERHANG_FAN_END\n";
    std::string input = ";_OVERHANG_FAN_START\nT1\n" + cleaning +
        end + "G1 X1 E1 F600\nG1 X9 E2 F60\n;_OVERHANG_FAN_START\nT0\n" + cleaning +
        ";_OVERHANG_FAN_END\nG1 X2 E2 F600\n";
    const auto cooled = cooling.process_layer(std::move(input), first_layer ? 0 : 2, true);
    FanMover mover(generator.writer(), 3.f, false, true, false, 0.f);
    const std::string output = mover.process_gcode(cooled, true);
    INFO(output);
    const auto tpu = output.find("T1\n");
    const auto asa = output.find("T0\n");
    REQUIRE(tpu != std::string::npos);
    REQUIRE(asa != std::string::npos);
    const std::string expected_tpu = first_layer || same_fan ? "M106 S0" : "M106 S102";
    // Each incoming tool receives a command immediately after T, not before it.
    CHECK(output.compare(tpu + 3, expected_tpu.size(), expected_tpu) == 0);
    CHECK(output.compare(asa + 3, 7, "M106 S0") == 0);
    if (custom_fan) {
        const auto tpu_custom = output.find("M106 S200", tpu);
        const auto asa_custom = output.find("M106 S200", asa);
        REQUIRE(tpu_custom != std::string::npos);
        REQUIRE(asa_custom != std::string::npos);
        const auto active_tpu = continuing_feature && !first_layer ? "M106 S204" : expected_tpu;
        CHECK(output.find(active_tpu, tpu_custom) < output.find("G1 X1"));
        CHECK(output.find("M106 S0", asa_custom) < output.find("G1 X2"));
    }
    if (continuing_feature && !first_layer)
        CHECK(output.find("M106 S204", tpu) < output.find("G1 X1"));
}

TEST_CASE("Fan look-ahead and kickstarts stay on their side of tool changes", "[Cooling][FanMover][Toolchange]")
{
    GCodeWriter writer;
    PrintConfig config;
    writer.apply_print_config(config);
    SECTION("Identical incoming fan commands survive the tool boundary") {
        FanMover mover(writer, 3.f, false, true, false, 0.f);
        const std::string output = mover.process_gcode(
            "M106 S102\nG1 X10 F600\nT1 S1 L0 D0\nM106 S102\nG1 X20 F600\n", true);
        const auto tool = output.find("T1 S1 L0 D0\n");
        REQUIRE(tool != std::string::npos);
        CHECK(output.find("M106 S102") < tool);
        const auto incoming = output.find("M106 S102", tool);
        REQUIRE(incoming != std::string::npos);
        CHECK(incoming < output.find("G1 X20"));
    }
    SECTION("An unfinished kickstart restores the outgoing fan before selection") {
        FanMover mover(writer, 0.f, false, true, false, 1.f);
        const std::string output = mover.process_gcode("M106 S51\nT1\nM106 S0\nG1 X10 F600\n", true);
        const auto tool = output.find("T1\n");
        REQUIRE(tool != std::string::npos);
        const auto outgoing = output.rfind("M106", tool);
        REQUIRE(outgoing != std::string::npos);
        CHECK(output.compare(outgoing, 8, "M106 S51") == 0);
        CHECK(output.find("M106 S0", tool) < output.find("G1 X10"));
    }
    SECTION("Unknown incoming fan state never kickstarts a fan-off command") {
        FanMover mover(writer, 3.f, false, true, false, 1.f);
        const std::string output = mover.process_gcode("T1\nM106 S0\nG1 X10 F600\n", true);
        CHECK(output.find("M106 S255") == std::string::npos);
        const auto tool = output.find("T1\n");
        REQUIRE(tool != std::string::npos);
        CHECK(output.find("M106 S0", tool) < output.find("G1 X10"));
    }
}

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
