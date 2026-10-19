#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <iterator>
#include "libslic3r/PlaceholderParser.hpp"
#include <catch2/catch_all.hpp>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Config.hpp"
#include "libslic3r/Extruder.hpp"
#include "libslic3r/Geometry.hpp"
#include "libslic3r/libslic3r.h"
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "nlohmann/json.hpp"

#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/GCodeWriter.hpp"
#include "libslic3r/Exception.hpp"
#include "libslic3r/GCode.hpp"
#include "libslic3r/GCode/SpoolManagerMetadata.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/ModelArrange.hpp"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>

#include "test_helpers.hpp"
#include <cmath>
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include <algorithm>
#include <limits>
#include "libslic3r/GCode/BeltKinematics.hpp"
#include "libslic3r/BeltTransform.hpp"
#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Arrange.hpp"

using namespace Slic3r;
using namespace Slic3r::Test;

TEST_CASE("RME startup plans use deposition temperatures and local bounds", "[GCodeWriter][RMEOutput]")
{
    const bool adaptive_pa = GENERATE(false, true);
    const bool inside_bed = GENERATE(false, true);
    auto config = multifilament_config(1, {
        {"filament_type", "HIPS"}, {"nozzle_temperature_initial_layer", "250"},
        {"nozzle_temperature", "240"}, {"skirt_loops", 0}, {"brim_type", "no_brim"},
        {"printable_area", "0x0,200x0,200x200,0x200"},
        {"layer_height", 0.2}, {"initial_layer_print_height", 0.2}
    });
    config.set_key_value("machine_start_gcode", new ConfigOptionString(
        std::string("; RME_DEPOSITION_BOUNDS=1\n; RME_TEMPERATURE_ADAPTIVE_PA=") +
        (adaptive_pa ? "1" : "0") +
        "\nM976 A 0:0:HIPS:250\nM104 S170\nM109 S250\nM555 X1 Y1 W1 H1\nG29 P1\n"));
    auto mesh = cube(5.);
    mesh.translate(inside_bed ? 80. : -20., 80., 0.);
    Print print;
    Model model;
    std::vector<TriangleMesh> meshes;
    meshes.push_back(std::move(mesh));
    init_print(std::move(meshes), print, model, config, nullptr, false);
    const std::string output = gcode(print);
    CHECK(output.find(adaptive_pa ? "M976 A 0:0:HIPS:250/240\n" : "M976 A 0:0:HIPS:250\n") != std::string::npos);
    CHECK(output.find("HIPS:250/170") == std::string::npos);
    const auto probe = output.find("G29 P1\n");
    REQUIRE(probe != std::string::npos);
    const auto bounds = output.rfind("M555 ", probe);
    REQUIRE(bounds != std::string::npos);
    CHECK(output.find('\n', bounds) == probe - 1);
    CHECK(output.substr(bounds, probe - bounds).find("X1 Y1 W1 H1") == std::string::npos);
    GCodeReader reader;
    reader.parse_buffer(output.substr(bounds, probe - bounds), [inside_bed](GCodeReader &, const GCodeReader::GCodeLine &line) {
        float width = 0, height = 0;
        REQUIRE(line.has_value('W', width));
        REQUIRE(line.has_value('H', height));
        if (inside_bed) {
            CHECK(width >= 5.);
            CHECK(width < 6.);
            CHECK(height >= 5.);
            CHECK(height < 6.);
        } else {
            CHECK_THAT(width, Catch::Matchers::WithinAbs(200., 1e-4));
            CHECK_THAT(height, Catch::Matchers::WithinAbs(200., 1e-4));
        }
    });
}

TEST_CASE("Smart bed temperature ignores elevated interface material", "[GCodeWriter][BedTemperature]")
{
    const std::string mode = GENERATE("by_highest_temp", "smart_bed_contact");
    const bool interface_on_bed = GENERATE(false, true);
    auto config = multifilament_config(2, {
        {"bed_temperature_formula", mode}, {"enable_support", true},
        {"support_type", "normal(auto)"}, {"support_filament", interface_on_bed ? 2 : 1},
        {"dont_support_bridges", false},
        {"support_interface_filament", 2}, {"support_interface_top_layers", 2},
        {"enable_prime_tower", false}, {"skirt_loops", 0}, {"brim_type", "no_brim"}
    });
    CHECK(config.option("bed_temperature_formula")->serialize() == mode);
    const BedType bed = config.opt_enum<BedType>("curr_bed_type");
    config.set_key_value(get_bed_temp_1st_layer_key(bed), new ConfigOptionInts{60, 85});
    config.set_key_value(get_bed_temp_key(bed), new ConfigOptionInts{55, 80});
    config.set_key_value("machine_start_gcode", new ConfigOptionString(
        "; BED_PLACEHOLDER {bed_temperature_initial_layer_single}\n"
        "; BED_VECTOR {first_layer_bed_temperature[1]} {bed_temperature[1]}\n"));
    Print print;
    Model model;
    init_print({TestMesh::overhang}, print, model, config);
    const std::string output = gcode(print);
    const auto first = print.get_slice_used_filaments(true);
    const auto all = print.get_slice_used_filaments(false);
    REQUIRE(std::find(all.begin(), all.end(), 1u) != all.end());
    const bool hotter_on_bed = std::find(first.begin(), first.end(), 1u) != first.end();
    REQUIRE(hotter_on_bed == interface_on_bed);
    const int expected_first = hotter_on_bed ? 85 : 60;
    const int expected_normal = mode == "by_highest_temp" || hotter_on_bed ? 80 : 55;
    CHECK(output.find("; BED_PLACEHOLDER " + std::to_string(expected_first)) != std::string::npos);
    if (mode == "smart_bed_contact")
        CHECK(output.find("; BED_VECTOR " + std::to_string(expected_first) + " " +
            std::to_string(expected_normal)) != std::string::npos);
    else
        CHECK(output.find("; BED_VECTOR 85 80") != std::string::npos);
    std::vector<int> bed_commands;
    GCodeReader reader;
    reader.parse_buffer(output, [&](GCodeReader &, const GCodeReader::GCodeLine &line) {
        float target;
        if ((line.cmd_is("M140") || line.cmd_is("M190")) && line.has_value('S', target) && target > 0)
            bed_commands.push_back(int(target));
    });
    REQUIRE_FALSE(bed_commands.empty());
    CHECK(bed_commands.front() == expected_first);
    CHECK(bed_commands.back() == expected_normal);
}

TEST_CASE("Preview heater targets preserve inactive tool temperatures and bed cooldown", "[GCodeWriter][TemperaturePreview]")
{
    GCodeProcessor processor;
    processor.initialize("temperature-preview.gcode");
    processor.process_buffer("G90\nM83\nM140 S60\nM104 T0 S210\nM109 T7 S240\nG1 X10 F600\n"
        "M190 R55\nM109 T0 R205\nG1 X20\nM104 T7 S0\nM140 S0\nG1 X30\n");
    const auto &result = processor.get_result();
    const auto *first = result.temperature_targets_at(7);
    REQUIRE(first != nullptr);
    REQUIRE(first->bed.has_value());
    CHECK_THAT(*first->bed, Catch::Matchers::WithinAbs(60., 0.01));
    REQUIRE(first->tools.count(0) == 1);
    REQUIRE(first->tools.count(7) == 1);
    CHECK_THAT(first->tools.at(0), Catch::Matchers::WithinAbs(210., 0.01));
    CHECK_THAT(first->tools.at(7), Catch::Matchers::WithinAbs(240., 0.01));
    const auto *cooled = result.temperature_targets_at(10);
    REQUIRE(cooled != nullptr);
    CHECK_THAT(*cooled->bed, Catch::Matchers::WithinAbs(55., 0.01));
    CHECK_THAT(cooled->tools.at(0), Catch::Matchers::WithinAbs(205., 0.01));
    CHECK_THAT(cooled->tools.at(7), Catch::Matchers::WithinAbs(240., 0.01));
    const auto *off = result.temperature_targets_at(13);
    REQUIRE(off != nullptr);
    CHECK_THAT(*off->bed, Catch::Matchers::WithinAbs(0., 0.01));
    CHECK_THAT(off->tools.at(7), Catch::Matchers::WithinAbs(0., 0.01));
    // Scrubbing backwards must not show the final shutdown temperatures.
    CHECK_THAT(*result.temperature_targets_at(7)->bed, Catch::Matchers::WithinAbs(60., 0.01));
    GCodeProcessorResult copy;
    copy = result;
    REQUIRE(copy.temperature_targets_at(7) != nullptr);
    CHECK_THAT(copy.temperature_targets_at(7)->tools.at(7), Catch::Matchers::WithinAbs(240., 0.01));
    processor.reset();
    CHECK(processor.get_result().temperature_targets_at(13) == nullptr);
}

TEST_CASE("Printer model overrides affect validator metadata without changing the profile",
          "[GCodeWriter][PrinterModel][Regression]")
{
    const bool skip_config = GENERATE(false, true);
    const std::string override_model = GENERATE(std::string(), std::string("  PrusaCoreOneINDX  "));
    std::string start = "; PROFILE_MODEL={printer_model}\nM862.3 P\"COREONEINDX\"\n";
    for (size_t i = 0; i < 1100; ++i)
        start += "; padding\n";
    Print print;
    Model model;
    init_print({cube(2.)}, print, model, {
        {"printer_model", "Prusa CORE One HF"},
        {"gcode_printer_model", override_model},
        {"gcode_skip_config_block", skip_config},
        {"machine_start_gcode", start},
        {"skirt_loops", 0}, {"brim_type", "no_brim"}
    });
    const std::string output = gcode(print);
    const std::string expected = override_model.empty() ? "Prusa CORE One HF" : "PrusaCoreOneINDX";
    CHECK(print.config().printer_model.value == "Prusa CORE One HF");
    CHECK(print.full_print_config().opt_string("printer_model") == "Prusa CORE One HF");
    CHECK(output.find("; PROFILE_MODEL=Prusa CORE One HF") != std::string::npos);
    CHECK(output.find("M862.3 P\"COREONEINDX\"") != std::string::npos);
    if (!skip_config || !override_model.empty()) {
        const auto position = output.rfind("; printer_model = " + expected + "\n");
        REQUIRE(position != std::string::npos);
        // The actual validator searches only the final 1000 lines.
        CHECK(std::count(output.begin() + position, output.end(), '\n') < 1000);
    } else
        CHECK(output.find("; printer_model = ") == std::string::npos);
    if (!override_model.empty())
        CHECK(output.find("; printer_model = Prusa CORE One HF\n") == std::string::npos);
}

TEST_CASE("Printer model metadata rejects embedded G-code lines",
          "[GCodeWriter][PrinterModel][Regression]")
{
    Print print;
    Model model;
    init_print({cube(2.)}, print, model, {{"gcode_printer_model", "INDX\nM112"}});
    CHECK(print.validate().opt_key == "gcode_printer_model");
    CHECK_THROWS_AS(gcode(print), SlicingError);
}

TEST_CASE("Unused configured tools stay unused in exported startup metadata",
          "[GCodeWriter][INDX][Regression]")
{
    const bool two_tools = GENERATE(false, true);
    auto config = multifilament_config(8, {{"enable_prime_tower", false},
        {"single_extruder_multi_material", true}, {"skirt_loops", 0},
        {"brim_type", "no_brim"}});
    config.set_key_value("filament_density", new ConfigOptionFloats(std::vector<double>(8, 1.24)));
    std::string start = "G90\nM83\n; USED_TOOLS = ";
    for (unsigned int tool = 0; tool < 8; ++tool)
        start += "{if is_extruder_used[" + std::to_string(tool) + "]}" +
                 std::to_string(tool) + ",{endif}";
    config.set_key_value("machine_start_gcode", new ConfigOptionString(start + "\n"));
    const std::string output = slice_with_object_overrides({cube(3), cube(3)}, config,
        {{{"extruder", two_tools ? 1 : 8}}, {{"extruder", 8}}});
    REQUIRE(output.find(two_tools ? "; USED_TOOLS = 0,7,\n" : "; USED_TOOLS = 7,\n") != std::string::npos);
    const std::string prefix = "; filament used [g] = ";
    const size_t begin = output.find(prefix);
    REQUIRE(begin != std::string::npos);
    std::istringstream values(output.substr(begin + prefix.size()));
    for (unsigned int tool = 0; tool < 8; ++tool) {
        double grams = -1.;
        REQUIRE(bool(values >> grams));
        CAPTURE(tool, two_tools);
        if (tool == 7 || (two_tools && tool == 0)) CHECK(grams > 0.);
        else CHECK(grams == 0.);
        if (tool < 7) {
            char comma = 0;
            values >> comma;
            REQUIRE(comma == ',');
        }
    }
}

TEST_CASE("INDX startup primes in the bucket and ejects before printing", "[GCodeWriter][INDX][Regression]")
{
    const bool bundled_profile = GENERATE(false, true);
    const bool tower = GENERATE(false, true);
    const bool heated_chamber = GENERATE(false, true);
    const unsigned int tool_mask = GENERATE(1u, 128u, 129u, 18u, 255u);
    std::vector<unsigned char> used(8, 0);
    std::string tool_list;
    int initial_tool = -1, tool_count = 0;
    for (int tool = 0; tool < 8; ++tool) {
        if ((tool_mask & (1u << tool)) == 0)
            continue;
        used[tool] = 1;
        if (initial_tool < 0) initial_tool = tool;
        if (tool_count++) tool_list += ',';
        tool_list += std::to_string(tool);
    }
    CAPTURE(tool_mask, tower, heated_chamber);
    auto config = multifilament_config(8, {{"enable_prime_tower", tower},
        {"outer_wall_line_width", 0.45}, {"line_width", 0.45}, {"layer_height", 0.2}});
    PlaceholderParser parser;
    parser.apply_config(config);
    parser.set("filament_type", new ConfigOptionStrings(std::vector<std::string>(8, heated_chamber ? "ABS" : "PLA")));
    parser.set("chamber_temperature", new ConfigOptionInts(std::vector<int>(8, heated_chamber ? 45 : 0)));
    parser.set("chamber_minimal_temperature", new ConfigOptionInts(std::vector<int>(8, heated_chamber ? 35 : 0)));
    parser.set("initial_tool", initial_tool);
    parser.set("is_extruder_used", new ConfigOptionBools(used));
    parser.set("first_layer_print_min", new ConfigOptionFloats({100., 100.}));
    parser.set("first_layer_print_max", new ConfigOptionFloats({120., 120.}));
    parser.set("adaptive_bed_mesh_min", new ConfigOptionFloats({99., 99.}));
    parser.set("adaptive_bed_mesh_max", new ConfigOptionFloats({121., 121.}));
    parser.set("print_bed_max", new ConfigOptionFloats({250., 220.}));
    parser.set("first_layer_bed_temperature", new ConfigOptionInts(std::vector<int>(8, 60)));
    parser.set("filament_minimal_purge_on_wipe_tower", new ConfigOptionFloats(std::vector<double>(8, 15.)));
    DynamicConfig outputs;
    outputs.set_key_value("e_retracted", new ConfigOptionFloats(std::vector<double>(8, 0.)));
    PlaceholderParser::ContextData context;
    context.global_config = std::make_unique<DynamicConfig>();
    const auto path = bundled_profile ?
        boost::filesystem::path(PROFILES_DIR) / "Prusa" / "machine" / "Prusa CORE One INDX RME 0.4 nozzle.json" :
        boost::filesystem::path(TEST_DATA_DIR) / "../../resources/gcode/indx_rme_start.gcode";
    std::ifstream input(path.string());
    REQUIRE(input.good());
    const std::string start = bundled_profile ? nlohmann::json::parse(input).at("machine_start_gcode").get<std::string>() :
        std::string(std::istreambuf_iterator<char>(input), {});
    const std::string output = parser.process(start,
                                             initial_tool, nullptr, &outputs, &context);
    CHECK(output.find('{') == std::string::npos);
    REQUIRE(output.find("\nG427 T" + tool_list + " R2 P3 ;") != std::string::npos);
    bool inside = false, pellet = false, wiped = false;
    size_t exits = 0, calibrations = 0, homes = 0, xy_homes = 0, z_homes = 0;
    size_t pa_calibrations = 0, tool_selections = 0, bed_waits = 0;
    bool bed_heating = false;
    GCodeReader reader;
    reader.parse_buffer(output, [&](GCodeReader &, const GCodeReader::GCodeLine &line) {
        float sequence = 0.f;
        if (line.cmd_is("M140") || line.cmd_is("M190")) {
            float target = 0.f;
            const bool has_target = line.has_value('S', target) || line.has_value('R', target);
            REQUIRE(has_target);
            CHECK(target > 0.f); // Never interrupt bed heating during startup.
            bed_heating = true;
            if (line.cmd_is("M190")) {
                CHECK(pa_calibrations == 1);
                ++bed_waits;
            }
        }
        if (line.cmd_is("M976")) {
            CHECK(z_homes == 1);
            CHECK(calibrations == 0);
            CHECK(bed_heating);
            ++pa_calibrations;
        }
        if (line.cmd_is("G427")) {
            CHECK(pa_calibrations == 1);
            CHECK(bed_waits == 1);
            ++calibrations;
        }
        if (!line.cmd().empty() && line.cmd().front() == 'T') {
            // One pickup for Z homing, then one print-tool restoration after
            // PA and offsets. No extra pickup between the calibration stages.
            if (tool_selections > 0) {
                CHECK(pa_calibrations == 1);
                CHECK(calibrations == 1);
                CHECK(line.cmd() == "T" + std::to_string(initial_tool));
            }
            ++tool_selections;
        }
        if (line.cmd_is("G28")) {
            ++homes;
            if (line.raw().rfind("G28 XY", 0) == 0) ++xy_homes;
            if (line.raw().rfind("G28 Z", 0) == 0) {
                CHECK(xy_homes == 1);
                ++z_homes;
            }
        }
        if (line.cmd_is("M870")) {
            CHECK(xy_homes == 1);
            CHECK(z_homes == 0);
        }
        if (line.cmd_is("G12") && line.has_value('S', sequence)) {
            if (sequence == 90) inside = true;
            if (sequence == 1) {
                CHECK(inside);
                wiped = true;
            }
            if (sequence == 30) {
                CHECK(wiped);
                wiped = false;
                pellet = false;
            }
            if (sequence == 91) {
                CHECK_FALSE(pellet);
                inside = false;
                ++exits;
            }
        }
        if (line.cmd_is("G1") && line.has_e() && line.e() > 0) {
            CHECK(inside); // No bed purge line or extrusion outside the cleaner.
            pellet = true;
            wiped = false;
        }
    });
    CHECK(exits == size_t(tool_count == 1 || tower));
    CHECK(calibrations == 1);
    CHECK(pa_calibrations == 1);
    CHECK(tool_selections == 2);
    CHECK(homes == 2);
    CHECK(xy_homes == 1);
    CHECK(z_homes == 1);
    CHECK_FALSE(inside);
}

TEST_CASE("INDX bucket templates eject pellets and emit evaluated filament commands", "[GCodeWriter][INDX][Regression]")
{
    const bool tower = GENERATE(false, true);
    const auto read_template = [](const char *filename) {
        std::ifstream input(std::string(TEST_DATA_DIR) + "/../../resources/gcode/" + filename);
        REQUIRE(input.good());
        return std::string(std::istreambuf_iterator<char>(input), {});
    };
    auto config = multifilament_config(2, {
        {"enable_prime_tower", tower}, {"single_extruder_multi_material", true},
        {"wipe_tower_x", "20"}, {"wipe_tower_y", "20"},
        {"layer_height", 0.2}, {"initial_layer_print_height", 0.2},
        {"skirt_loops", 0}, {"brim_type", "no_brim"}, {"enable_pressure_advance", "0,0"}
    });
    config.set_key_value("machine_start_gcode", new ConfigOptionString("{global retract_toolchange = 8}\nG90\nM83\n"));
    config.set_key_value("change_filament_gcode", new ConfigOptionString(read_template("indx_rme_toolchange.gcode")));
    config.set_key_value("filament_start_gcode", new ConfigOptionStrings({
        "M572 S{if nozzle_diameter[0]==0.4}0.036{else}0{endif} ; evaluated filament zero\n",
        "M572 S{if nozzle_diameter[0]==0.4}0.026{else}0{endif} ; evaluated filament one\n"
    }));
    const std::string output = slice_with_object_overrides(
        {cube(5), cube(5)}, config, {{{"extruder", 1}}, {{"extruder", 2}}});
    bool inside_cleaner = false;
    bool positive_e_since_eject = false;
    bool wiped = false;
    size_t cleaner_exits = 0, filament_commands = 0;
    GCodeReader reader;
    reader.parse_buffer(output, [&](GCodeReader &r, const GCodeReader::GCodeLine &line) {
        if (line.cmd().empty())
            return; // The footer intentionally stores unevaluated profile templates.
        CHECK(line.raw().find('{') == std::string::npos);
        CHECK(line.raw().find('[') == std::string::npos);
        if (line.cmd_is("M572")) {
            ++filament_commands;
            float pressure = 0.f;
            REQUIRE(line.has_value('S', pressure));
            CHECK((std::abs(pressure - 0.036) < 1e-5 || std::abs(pressure - 0.026) < 1e-5));
        }
        if (line.cmd_is("G12")) {
            float sequence = 0.f;
            REQUIRE(line.has_value('S', sequence));
            if (sequence == 90) inside_cleaner = true;
            if (sequence == 1) {
                CHECK(inside_cleaner);
                wiped = true;
            }
            if (sequence == 30) {
                CHECK(inside_cleaner);
                CHECK(wiped);
                wiped = false;
                positive_e_since_eject = false;
            }
            if (sequence == 91) {
                CHECK(inside_cleaner);
                CHECK_FALSE(positive_e_since_eject);
                inside_cleaner = false;
                ++cleaner_exits;
            }
        }
        if (inside_cleaner && line.cmd_is("G1") && line.has_e() && line.e() > 0) {
            positive_e_since_eject = true;
            wiped = false;
        }
    });
    CHECK(cleaner_exits > 0);
    CHECK(filament_commands >= 2);
    CHECK_FALSE(inside_cleaner);
}

TEST_CASE("SpoolManager metadata updates marked filament notes", "[GCodeWriter][SpoolManager]")
{
    const std::string tail =
        "; filament_type = PLA;PETG;ABS\n"
        "; filament used [mm] = 10, 20\n"
        "; filament_notes = [sm_name=Old]; [sm_name = ]; ordinary note\n";

    const std::vector<SpoolManagerMetadata::Filament> filaments{
        {"Purple PLA", "PLA", "#800080", "Purple"},
        {"Blue PETG", "PETG", "#0000FF", "Blue"},
        {"Ignored", "ABS", "#000000", "Black"}
    };
    const std::string updated = SpoolManagerMetadata::update_gcode_tail(tail, filaments);

    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring("; filament used [mm] = 10, 20, 0"));
    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring(
                            "; filament_notes = [sm_name = Purple PLA]; [sm_name = Blue PETG]; ordinary note"));
    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring(
                            "; spool_manager_filament_names = Purple PLA;Blue PETG;Ignored"));
    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring(
                            "; spool_manager_filament_materials = PLA;PETG;ABS"));
    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring(
                            "; spool_manager_filament_colors = #800080;#0000FF;#000000"));
}

TEST_CASE("SpoolManager embeds explicit metadata without legacy note markers", "[GCodeWriter][SpoolManager]")
{
    const std::string tail =
        "; filament_type = PLA\n"
        "; filament used [mm] = 10\n"
        "; filament_notes = ordinary note\n";

    const std::string updated = SpoolManagerMetadata::update_gcode_tail(
        tail, {{"Purple PLA", "PLA", "#800080", "Purple"}});
    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring(
                            "; filament_notes = ordinary note [sm_name = Purple PLA]"));
    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring("; spool_manager_filament_names = Purple PLA"));
}

TEST_CASE("SpoolManager embeds explicit metadata without filament notes", "[GCodeWriter][SpoolManager]")
{
    const std::string tail =
        "; filament_type = PLA\n"
        "; filament used [mm] = 10\n";

    const std::string updated = SpoolManagerMetadata::update_gcode_tail(
        tail, {{"Purple PLA", "PLA", "#800080", "Purple"}});
    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring(
                            "; filament_notes = [sm_name = Purple PLA]"));
    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring("; spool_manager_filament_names = Purple PLA"));
    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring("; spool_manager_filament_materials = PLA"));
    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring("; spool_manager_filament_colors = #800080"));
}

TEST_CASE("Spoolman metadata uses the Nozzle Filament Validator stable identifier",
          "[GCodeWriter][SpoolManager]")
{
    const std::string tail =
        "; filament_type = PLA\n"
        "; filament used [mm] = 10\n"
        "; filament_notes = Keep this [sm_name = old value]\n";
    SpoolManagerMetadata::Filament filament;
    filament.name = "Galaxy Black";
    filament.material = "PLA";
    filament.spool_id = "42";
    filament.provider = "Spoolman";

    const std::string updated = SpoolManagerMetadata::update_gcode_tail(tail, {filament});
    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring(
                            "; filament_notes = Keep this [sm_name = spoolman:42]"));
    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring(
                            "; spool_manager_filament_spool_ids = 42"));
    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring(
                            "; spool_manager_filament_providers = Spoolman"));
    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring(
                            "; spool_manager_filament_validation_ids = spoolman:42"));
}

TEST_CASE("RME compatibility metadata uses its provider-qualified validator identifier",
          "[GCodeWriter][SpoolManager]")
{
    const std::string tail =
        "; filament_type = PLA;PETG\n"
        "; filament used [mm] = 10, 20\n"
        "; filament_notes = [sm_name = old];[sm_name = old]\n";
    SpoolManagerMetadata::Filament internal;
    internal.name = "RME PLA";
    internal.spool_id = "4";
    internal.provider = "RME compatibility";
    internal.inventory_provider = "internal";
    SpoolManagerMetadata::Filament delegated;
    delegated.name = "External PETG";
    delegated.spool_id = "91";
    delegated.provider = "RME compatibility";
    delegated.inventory_provider = "Spoolman Cloud";

    const std::string updated = SpoolManagerMetadata::update_gcode_tail(
        tail, {internal, delegated});
    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring(
                            "; filament_notes = [sm_name = rme:internal:4];[sm_name = rme:spoolman-cloud:91]"));
    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring(
                            "; spool_manager_filament_spool_ids = 4;91"));
    CHECK_THAT(updated, Catch::Matchers::ContainsSubstring(
                            "; spool_manager_filament_validation_ids = rme:internal:4;rme:spoolman-cloud:91"));

    // Upload retry/reprocessing must replace the metadata block, not append a
    // second set of identifiers that a host could parse ambiguously.
    const std::string repeated = SpoolManagerMetadata::update_gcode_tail(updated, {internal, delegated});
    const std::string key = "; spool_manager_filament_validation_ids = ";
    REQUIRE(repeated.find(key) != std::string::npos);
    CHECK(repeated.find(key, repeated.find(key) + key.size()) == std::string::npos);
}

TEST_CASE("SpoolManager selected spool parser preserves tool slot indices", "[GCodeWriter][SpoolManager]")
{
    const std::string response = R"({
        "allSpools": [],
        "selectedSpools": [
            {"displayName":"Galaxy PLA","vendor":"Prusament","material":"PLA","color":"#123456","colorName":"Galaxy"},
            null,
            {"displayName":"Safety PETG","material":"PETG","color":"#FF8800","colorName":"Orange"}
        ]
    })";
    std::vector<SpoolManagerMetadata::Filament> slots;
    std::string error;

    REQUIRE(SpoolManagerMetadata::parse_selected_spools(
        response, slots, error, "RME compatibility"));
    REQUIRE(slots.size() == 3);
    CHECK(slots[0].name == "Galaxy PLA");
    CHECK(slots[0].material == "PLA");
    CHECK(slots[0].color == "#123456");
    CHECK(slots[0].vendor == "Prusament");
    CHECK(slots[1].name.empty());
    CHECK(slots[2].name == "Safety PETG");
}

TEST_CASE("SpoolManager selected spool parser rejects an unassigned machine", "[GCodeWriter][SpoolManager]")
{
    std::vector<SpoolManagerMetadata::Filament> slots;
    std::string error;

    CHECK_FALSE(SpoolManagerMetadata::parse_selected_spools(
        R"({"allSpools":[],"selectedSpools":[null,null]})", slots, error));
    CHECK_FALSE(error.empty());
    CHECK(slots.empty());
}

TEST_CASE("Spoolman tool assignments resolve against its spool inventory", "[GCodeWriter][SpoolManager]")
{
    const std::string response = R"({
        "selectedSpoolIds":{"0":{"spoolId":"42"},"2":{"spoolId":"73"}},
        "spools":[
            {"id":42,"filament":{"name":"Galaxy Black","material":"PLA","color_hex":"112233",
             "vendor":{"id":8,"name":"Prusament"}}},
            {"id":73,"filament":{"name":"Signal Orange","material":"PETG","color_hex":"FF6600",
             "vendor":{"name":"Polymaker"}}}
        ]
    })";
    std::vector<SpoolManagerMetadata::Filament> slots;
    std::string error;

    REQUIRE(SpoolManagerMetadata::parse_selected_spools(response, slots, error));
    REQUIRE(slots.size() == 3);
    CHECK(slots[0].name == "Galaxy Black");
    CHECK(slots[0].vendor == "Prusament");
    CHECK(slots[0].material == "PLA");
    CHECK(slots[0].color == "112233");
    CHECK(slots[0].spool_id == "42");
    CHECK(slots[0].provider == "Spoolman");
    CHECK(slots[1].name.empty());
    CHECK(slots[2].name == "Signal Orange");
}

TEST_CASE("OctoPrint spool inventory exposes every upstream roll for profile mapping",
          "[GCodeWriter][SpoolManager]")
{
    const std::string response = R"({
        "selectedSpoolIds":{"0":{"spoolId":"42"}},
        "spools":[
            {"id":42,"filament":{"name":"Galaxy Black","material":"PLA",
             "vendor":{"name":"Prusament"}}},
            {"id":73,"filament":{"name":"Signal Orange","material":"PETG",
             "vendor":{"name":"Polymaker"}}}
        ]
    })";
    std::vector<SpoolManagerMetadata::Filament> inventory;
    std::string error;

    REQUIRE(SpoolManagerMetadata::parse_spool_inventory(
        response, inventory, error, "Spoolman"));
    REQUIRE(inventory.size() == 2);
    CHECK(inventory[0].spool_id == "42");
    CHECK(inventory[0].provider == "Spoolman");
    CHECK(inventory[1].spool_id == "73");
    CHECK(inventory[1].material == "PETG");
}

TEST_CASE("legacy SpoolManager inventory is available without a tool assignment",
          "[GCodeWriter][SpoolManager]")
{
    const std::string response = R"({
        "allSpools":[
            {"spoolId":"9","displayName":"Workshop Black","vendor":"Generic","material":"ABS"}
        ],
        "selectedSpools":[null]
    })";
    std::vector<SpoolManagerMetadata::Filament> inventory;
    std::string error;

    REQUIRE(SpoolManagerMetadata::parse_spool_inventory(response, inventory, error));
    REQUIRE(inventory.size() == 1);
    CHECK(inventory.front().name == "Workshop Black");
    CHECK(inventory.front().provider == "SpoolManager");
}

TEST_CASE("generic OctoPrint filament providers preserve tool indices", "[GCodeWriter][SpoolManager]")
{
    const std::string response = R"({
        "provider":"RME compatibility",
        "data":{"tools":{
            "1":{"spool":{"id":"roll-b","name":"Blue TPU","material":"TPU",
                  "color":"#0000FF","manufacturer":"Example Filament"}},
            "3":null
        }}
    })";
    std::vector<SpoolManagerMetadata::Filament> slots;
    std::string error;

    REQUIRE(SpoolManagerMetadata::parse_selected_spools(response, slots, error));
    REQUIRE(slots.size() == 4);
    CHECK(slots[0].name.empty());
    CHECK(slots[1].name == "Blue TPU");
    CHECK(slots[1].vendor == "Example Filament");
    CHECK(slots[1].spool_id == "roll-b");
    CHECK(slots[1].provider == "RME compatibility");
    CHECK(slots[3].name.empty());
}

TEST_CASE("RME compatibility supplies its provider name when omitted",
          "[GCodeWriter][SpoolManager]")
{
    const std::string response = R"({"data":{"tools":{"0":{
        "spool":{"id":"rme-roll-7","name":"RME PLA","material":"PLA"}
    }}}})";
    std::vector<SpoolManagerMetadata::Filament> slots;
    std::string error;

    REQUIRE(SpoolManagerMetadata::parse_selected_spools(
        response, slots, error, "RME compatibility"));
    REQUIRE(slots.size() == 1);
    CHECK(slots[0].spool_id == "rme-roll-7");
    CHECK(slots[0].provider == "RME compatibility");
}

TEST_CASE("RME compatibility preserves the inventory provider namespace",
          "[GCodeWriter][SpoolManager]")
{
    const std::string response = R"({"schema":"rme-filament-report-v1","provider":"spoolmanager","data":{"tools":[{
        "spool_id":"27","name":"RME delegated PLA","material":"PLA"
    }]}})";
    std::vector<SpoolManagerMetadata::Filament> slots;
    std::string error;

    REQUIRE(SpoolManagerMetadata::parse_selected_spools(response, slots, error));
    REQUIRE(slots.size() == 1);
    CHECK(slots[0].spool_id == "27");
    CHECK(slots[0].provider == "RME compatibility");
    CHECK(slots[0].inventory_provider == "spoolmanager");
}

TEST_CASE("FilamentManager selections expose nested spool profiles", "[GCodeWriter][SpoolManager]")
{
    const std::string response = R"({"selections":[
        {"tool":0,"spool":{"id":9,"name":"Workshop Black","profile":{
            "vendor":"Generic","material":"ABS","color":"#202020"}}},
        {"tool":1,"spool":null}
    ]})";
    std::vector<SpoolManagerMetadata::Filament> slots;
    std::string error;

    REQUIRE(SpoolManagerMetadata::parse_selected_spools(response, slots, error));
    REQUIRE(slots.size() == 2);
    CHECK(slots[0].name == "Workshop Black");
    CHECK(slots[0].material == "ABS");
    CHECK(slots[0].vendor == "Generic");
    CHECK(slots[0].spool_id == "9");
    CHECK(slots[0].provider == "FilamentManager");
    CHECK(slots[1].name.empty());
}

TEST_CASE("OctoPrint spool mappings prefer a roll over material and default", "[GCodeWriter][SpoolManager]")
{
    SpoolManagerMetadata::Filament filament;
    filament.spool_id = "42";
    filament.provider = "Spoolman";
    filament.vendor = "Prusament";
    filament.material = "PLA";

    CHECK(SpoolManagerMetadata::mapped_profile_name(
              filament,
              {"spoolman:42=Exact Roll PLA", "42=Fallback Roll PLA"},
              {"prusament|pla=Brand PLA", "pla=Generic PLA"},
              "Default PLA") == "Exact Roll PLA");
    filament.spool_id = "99";
    CHECK(SpoolManagerMetadata::mapped_profile_name(
              filament, {}, {" PRUSAMENT | PLA = Brand PLA"}, "Default PLA") == "Brand PLA");
    filament.vendor = "Unknown";
    CHECK(SpoolManagerMetadata::mapped_profile_name(
              filament, {}, {}, "Default PLA") == "Default PLA");
}

TEST_CASE("Belt machine coordinates retain a non-45-degree slicing angle", "[GCodeWriter][belt]")
{
    PrintConfig config;
    config.belt_printer.value               = true;
    config.belt_slice_rotation.value        = BeltRotationAxis::X;
    config.belt_slice_rotation_angle.value  = 30.;
    config.gcode_remap_x.value              = RemapAxis::PosX;
    config.gcode_remap_y.value              = RemapAxis::PosZ;
    config.gcode_remap_z.value              = RemapAxis::PosY;

    GCodeWriter writer;
    install_belt_kinematics(writer, config);
    writer.set_axis_remap(int(config.gcode_remap_x.value),
                          int(config.gcode_remap_y.value),
                          int(config.gcode_remap_z.value));

    // Start with a point in the unrotated model frame, then feed the writer the
    // same rotated coordinate produced by the pre-slice mesh transform. The
    // back-transform must recover the model point before the axis swap and
    // machine-frame shear/scale are applied.
    const Vec3d model(4., 10., 3.);
    Transform3d forward = BeltTransformPipeline::build_forward_transform(config);
    const Vec3d machine = writer.kinematics().to_machine(forward * model);

    // The conventional X-tilt remap produces (x, z, y). At 30 degrees the
    // gantry coordinate is z/sin(30) and belt travel is y + z*cot(30).
    // The complementary tan/inv-cos formulas accidentally used by the unified
    // transform are indistinguishable at 45 degrees, but fail this case.
    REQUIRE_THAT(machine.x(), Catch::Matchers::WithinAbs(4., 1e-9));
    REQUIRE_THAT(machine.y(), Catch::Matchers::WithinAbs(3. / std::sin(Geometry::deg2rad(30.)), 1e-9));
    REQUIRE_THAT(machine.z(), Catch::Matchers::WithinAbs(10. + 3. / std::tan(Geometry::deg2rad(30.)), 1e-9));
}

// Arrange on a 500x500 bed, which keeps coordinates small while still covering large printers.
static void arrange_objects_on_test_bed(Model &model, const DynamicPrintConfig &config)
{
    const BoundingBox bed{Point::new_scale(0., 0.), Point::new_scale(500., 500.)};
    arrange_objects(model, bed, ArrangeParams{scaled(min_object_distance(config))});
}

SCENARIO("set_speed emits values with fixed-point output.", "[GCodeWriter]") {

    GIVEN("GCodeWriter instance") {
        GCodeWriter writer;
        WHEN("set_speed is called to set speed to 99999.123") {
            THEN("Output string is G1 F99999.123") {
                REQUIRE_THAT(writer.set_speed(99999.123), Catch::Matchers::Equals("G1 F99999.123\n"));
            }
        }
        WHEN("set_speed is called to set speed to 1") {
            THEN("Output string is G1 F1") {
                REQUIRE_THAT(writer.set_speed(1.0), Catch::Matchers::Equals("G1 F1\n"));
            }
        }
        WHEN("set_speed is called to set speed to 203.200022") {
            THEN("Output string is G1 F203.2") {
                REQUIRE_THAT(writer.set_speed(203.200022), Catch::Matchers::Equals("G1 F203.2\n"));
            }
        }
        WHEN("set_speed is called to set speed to 203.200522") {
            THEN("Output string is G1 F203.201") {
                REQUIRE_THAT(writer.set_speed(203.200522), Catch::Matchers::Equals("G1 F203.201\n"));
            }
        }
    }
}

SCENARIO("z_hop lifts the nozzle when a lift is requested", "[GCodeWriter]") {
    GIVEN("A writer with the nozzle parked at Z = 10") {
        GCodeWriter writer;
        std::vector<unsigned int> extruder_ids { 0 };
        writer.set_extruders(extruder_ids);
        writer.set_extruder(0);
        writer.travel_to_z(10.0);

        WHEN("z_hop is 1 and an eager lift is requested") {
            writer.config.z_hop.values = { 1.0 };
            std::string gcode = writer.eager_lift(LiftType::NormalLift);
            THEN("a Z move up by z_hop is emitted") {
                REQUIRE_THAT(gcode, Catch::Matchers::ContainsSubstring("Z11"));
            }
        }
        WHEN("z_hop is 0") {
            writer.config.z_hop.values = { 0.0 };
            std::string gcode = writer.eager_lift(LiftType::NormalLift);
            THEN("no lift is emitted") {
                REQUIRE(gcode.empty());
            }
        }
    }
}

SCENARIO("Origin manipulation", "[GCodeWriter]") {
	Slic3r::GCode gcodegen;
	WHEN("set_origin to (10,0)") {
    	gcodegen.set_origin(Vec2d(10,0));
    	REQUIRE(gcodegen.origin() == Vec2d(10, 0));
    }
	WHEN("set_origin to (10,0) and translate by (5, 5)") {
		gcodegen.set_origin(Vec2d(10,0));
		gcodegen.set_origin(gcodegen.origin() + Vec2d(5, 5));
		THEN("origin returns reference to point") {
    		REQUIRE(gcodegen.origin() == Vec2d(15,5));
    	}
    }
}

TEST_CASE("A cached config slot is looked up again whenever its key changes", "[GCodeWriter]")
{
    GCode::ConfigIndexCache cache;
    int lookups = 0;
    auto slot = [&](int filament, size_t layer, size_t generation) {
        return cache.get(filament, layer, generation, [&] { ++lookups; return filament * 100 + int(layer) * 10 + int(generation); });
    };
    REQUIRE(slot(1, 2, 3) == 123);
    REQUIRE(slot(1, 2, 3) == 123);
    REQUIRE(lookups == 1);
    REQUIRE(slot(4, 2, 3) == 423);
    REQUIRE(slot(4, 5, 3) == 453);
    REQUIRE(slot(4, 5, 6) == 456);
    REQUIRE(lookups == 4);
}

// Verify that emit_machine_limits_to_gcode emits the correct max value across
// used extruders (regression for commit b4ee665: "Emit max value of machine
// limit among used extruders").
TEST_CASE("Machine envelope emits max limit among used extruders", "[GCodeWriter]")
{
    SECTION("Single extruder emits its configured values") {
        const std::string gcode = Slic3r::Test::slice({ cube(20) }, {
            { "emit_machine_limits_to_gcode", "1" },
            { "gcode_flavor",                "marlin2" },
            { "gcode_comments",              "1" },
            { "machine_start_gcode",         "" },
            { "layer_height",                "0.2" },
            { "initial_layer_print_height",  "0.2" },
            { "initial_layer_line_width",    "0" },
            { "z_hop",                       "0" },
            // stride-2 options: (normal, silent)
            { "machine_max_acceleration_x",           "500,600" },
            { "machine_max_acceleration_y",           "700,800" },
            { "machine_max_acceleration_z",           "100,200" },
            { "machine_max_acceleration_e",           "5000,6000" },
            { "machine_max_acceleration_extruding",   "1200,1300" },
            { "machine_max_acceleration_retracting",  "1400,1500" },
            { "machine_max_acceleration_travel",      "1600,1700" },
            // stride-2 options: (normal, silent)
            { "machine_max_speed_x",  "100,100" },
            { "machine_max_speed_y",  "110,110" },
            { "machine_max_speed_z",  "10,10" },
            { "machine_max_speed_e",  "50,50" },
            { "machine_max_jerk_x",   "8,8" },
            { "machine_max_jerk_y",   "9,9" },
            { "machine_max_jerk_z",   "0.4,0.4" },
            { "machine_max_jerk_e",   "5,5" },
            { "machine_max_junction_deviation", "0.02,0.03" },
        });

        THEN("M201 uses the normal acceleration values") {
            REQUIRE(gcode.find("M201 X500 Y700 Z100 E5000") != std::string::npos);
        }
        THEN("M203 uses the speed values") {
            REQUIRE(gcode.find("M203 X100 Y110 Z10 E50") != std::string::npos);
        }
        THEN("M204 (Marlin 2) uses extruding / retracting / travel") {
            REQUIRE(gcode.find("M204 P1200 R1400 T1600") != std::string::npos);
        }
        THEN("M205 uses the jerk values") {
            REQUIRE(gcode.find("M205 X8.00 Y9.00 Z0.40 E5.00") != std::string::npos);
        }
        THEN("M205 J uses the junction deviation") {
            REQUIRE(gcode.find("M205 J0.020") != std::string::npos);
        }
    }

    SECTION("Legacy Marlin flavor emits correct format") {
        const std::string gcode = Slic3r::Test::slice({ cube(20) }, {
            { "emit_machine_limits_to_gcode", "1" },
            { "gcode_flavor",                "marlin" },
            { "gcode_comments",              "1" },
            { "machine_start_gcode",         "" },
            { "layer_height",                "0.2" },
            { "initial_layer_print_height",  "0.2" },
            { "initial_layer_line_width",    "0" },
            { "z_hop",                       "0" },
            // All machine limits must be provided — defaults are empty vectors.
            { "machine_max_acceleration_x",           "500,600" },
            { "machine_max_acceleration_y",           "500,600" },
            { "machine_max_acceleration_z",           "500,600" },
            { "machine_max_acceleration_e",           "5000,6000" },
            { "machine_max_acceleration_extruding",   "1200,1300" },
            { "machine_max_acceleration_retracting",  "1400,1500" },
            { "machine_max_acceleration_travel",      "1600,1700" },
            { "machine_max_speed_x",  "100,100" },
            { "machine_max_speed_y",  "110,110" },
            { "machine_max_speed_z",  "10,10" },
            { "machine_max_speed_e",  "50,50" },
            { "machine_max_jerk_x",   "8,8" },
            { "machine_max_jerk_y",   "9,9" },
            { "machine_max_jerk_z",   "0.4,0.4" },
            { "machine_max_jerk_e",   "5,5" },
            { "machine_max_junction_deviation", "0.02,0.03" },
        });

        THEN("Legacy Marlin: M204 travel_acc = extruding_acc") {
            // gcfMarlinLegacy uses extruding acc for travel too
            REQUIRE(gcode.find("M204 P1200 R1400 T1200") != std::string::npos);
        }
        THEN("Legacy Marlin: M205 uses mm/sec format") {
            REQUIRE(gcode.find("M205 X8.00 Y9.00 Z0.40 E5.00") != std::string::npos);
        }
    }

    SECTION("Multi extruder - max of used extruders is emitted") {
        // Build config with 2 extruders that have *different* machine limits.
        // Extruder 1 has higher values; the emitted G-code must use the max.
        DynamicPrintConfig config = DynamicPrintConfig::full_print_config();

        // Print basics
        config.set_key_value("emit_machine_limits_to_gcode", new ConfigOptionBool(true));
        config.set_key_value("gcode_flavor",                 new ConfigOptionEnum<GCodeFlavor>(gcfMarlinFirmware));
        config.set_key_value("gcode_comments",               new ConfigOptionBool(true));
        config.set_key_value("machine_start_gcode",          new ConfigOptionString(""));
        config.set_key_value("layer_height",                 new ConfigOptionFloat(0.2));
        config.set_key_value("initial_layer_print_height",   new ConfigOptionFloat(0.2));
        config.set_key_value("initial_layer_line_width",     new ConfigOptionFloatOrPercent(0, false));
        config.set_key_value("z_hop",                        new ConfigOptionFloats({0}));
        // Print objects sequentially so each uses its own extruder without
        // wipe-tower / tool-change complexity.
        config.set_key_value("print_sequence",               new ConfigOptionEnum<PrintSequence>(PrintSequence::ByObject));

        // 2 extruders
        config.set_key_value("nozzle_diameter",          new ConfigOptionFloats({0.4, 0.4}));
        config.set_key_value("printer_extruder_id",      new ConfigOptionInts({1, 2}));
        config.set_key_value("printer_extruder_variant", new ConfigOptionStrings({"Direct Drive Standard", "Direct Drive Standard"}));
        config.set_key_value("filament_diameter",        new ConfigOptionFloats({1.75, 1.75}));
        config.set_key_value("filament_colour",          new ConfigOptionStrings({"#FF0000", "#00FF00"}));
        config.set_key_value("filament_type",            new ConfigOptionStrings({"PLA", "PLA"}));
        // filament_map maps filament slot index (1-based) → logical extruder ID (1-based).
        // Default [1] maps everything to extruder 0. Need [1, 2] for two distinct extruders.
        // fmmManual prevents auto-computation from overwriting the explicit mapping.
        config.option<ConfigOptionEnum<FilamentMapMode>>("filament_map_mode", true)->value = fmmManual;
        config.set_key_value("filament_map",             new ConfigOptionInts({1, 2}));
        config.set_key_value("default_filament_colour",  new ConfigOptionStrings({"#FF0000", "#00FF00"}));
        config.set_key_value("nozzle_temperature",       new ConfigOptionInts({210, 210}));
        config.set_key_value("nozzle_temperature_range_low",  new ConfigOptionInts({190, 190}));
        config.set_key_value("nozzle_temperature_range_high", new ConfigOptionInts({240, 240}));
        // flush_volumes_matrix must be filament_count^2 * heads_count entries.
        // 2 filaments * 2 * 1 head = 4 entries (all zero — flush volumes not tested here).
        config.set_key_value("flush_multiplier",     new ConfigOptionFloats({1}));
        config.set_key_value("flush_volumes_matrix", new ConfigOptionFloats({0, 0, 0, 0}));

        // Machine limits: extruder 0 low, extruder 1 high
        // Stride-2 (normal, silent pairs): e0_n, e0_s, e1_n, e1_s
        config.set_key_value("machine_max_acceleration_x",          new ConfigOptionFloats({500, 0, 1000, 0}));
        config.set_key_value("machine_max_acceleration_y",          new ConfigOptionFloats({700, 0, 1100, 0}));
        config.set_key_value("machine_max_acceleration_z",          new ConfigOptionFloats({100, 0, 300, 0}));
        config.set_key_value("machine_max_acceleration_e",          new ConfigOptionFloats({5000, 0, 8000, 0}));
        config.set_key_value("machine_max_acceleration_extruding",  new ConfigOptionFloats({1200, 0, 2200, 0}));
        config.set_key_value("machine_max_acceleration_retracting", new ConfigOptionFloats({1400, 0, 2400, 0}));
        config.set_key_value("machine_max_acceleration_travel",     new ConfigOptionFloats({1600, 0, 2600, 0}));
        config.set_key_value("machine_max_speed_x",  new ConfigOptionFloats({100, 0, 200, 0}));
        config.set_key_value("machine_max_speed_y",  new ConfigOptionFloats({110, 0, 210, 0}));
        config.set_key_value("machine_max_speed_z",  new ConfigOptionFloats({10, 0, 30, 0}));
        config.set_key_value("machine_max_speed_e",  new ConfigOptionFloats({50, 0, 80, 0}));
        config.set_key_value("machine_max_jerk_x",   new ConfigOptionFloats({8, 0, 12, 0}));
        config.set_key_value("machine_max_jerk_y",   new ConfigOptionFloats({9, 0, 13, 0}));
        config.set_key_value("machine_max_jerk_z",   new ConfigOptionFloats({0.4, 0, 0.6, 0}));
        config.set_key_value("machine_max_jerk_e",   new ConfigOptionFloats({5, 0, 10, 0}));
        config.set_key_value("machine_max_junction_deviation", new ConfigOptionFloats({0.02, 0, 0.05, 0}));

        // Model: two objects assigned to different extruders
        Model model;
        auto* obj1 = model.add_object();
        obj1->add_volume(cube(20));
        obj1->add_instance();
        // obj1 uses default extruder=1 (0-based index 0)

        auto* obj2 = model.add_object();
        obj2->add_volume(cube(20));
        obj2->add_instance();
        obj2->config.set_key_value("extruder", new ConfigOptionInt(2)); // 0-based index 1

        Print print;
        arrange_objects_on_test_bed(model, config);
        for (auto* mo : model.objects) {
            mo->ensure_on_bed();
            print.auto_assign_extruders(mo);
        }

        print.apply(model, config);
        print.validate();
        print.set_status_silent();
        print.process();

        std::string gcode = Slic3r::Test::gcode(print);

        THEN("M201 contains max (extruder 1's) acceleration values") {
            REQUIRE(gcode.find("M201 X1000 Y1100 Z300 E8000") != std::string::npos);
        }
        THEN("M203 contains max speed values") {
            REQUIRE(gcode.find("M203 X200 Y210 Z30 E80") != std::string::npos);
        }
        THEN("M204 contains max extruding / retracting / travel") {
            REQUIRE(gcode.find("M204 P2200 R2400 T2600") != std::string::npos);
        }
        THEN("M205 contains max jerk values") {
            REQUIRE(gcode.find("M205 X12.00 Y13.00 Z0.60 E10.00") != std::string::npos);
        }
        THEN("M205 contains max m_max_junction_deviation ") {
            REQUIRE(gcode.find("M205 J0.050") != std::string::npos);
        }
    }
}

// Verify that the EXTRUDER_LIMIT macro (GCodeWriter.cpp) correctly:
//  1) Uses the active extruder's specific limit when filament() is known.
//  2) Falls back to the maximum of all extruder limits when filament() is nullptr.
//
// These two behaviours were introduced in:
//  - "Use per-extruder motion limit" (1ab34a7454)
//  - "Use max limit when current extruder is unknown" (b7240ab1c6)
TEST_CASE("EXTRUDER_LIMIT per-extruder clamping and max fallback", "[GCodeWriter]")
{
    // --- Build config with 2 extruders that have different machine limits ---
    // Extruder 0: low limits
    // Extruder 1: high limits
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();

    config.set_key_value("emit_machine_limits_to_gcode", new ConfigOptionBool(true));
    config.set_key_value("gcode_flavor",                 new ConfigOptionEnum<GCodeFlavor>(gcfMarlinFirmware));
    config.set_key_value("gcode_comments",               new ConfigOptionBool(true));
    config.set_key_value("machine_start_gcode",          new ConfigOptionString(""));
    config.set_key_value("layer_height",                 new ConfigOptionFloat(0.2));
    config.set_key_value("initial_layer_print_height",   new ConfigOptionFloat(0.2));
    config.set_key_value("initial_layer_line_width",     new ConfigOptionFloatOrPercent(0, false));
    config.set_key_value("z_hop",                        new ConfigOptionFloats({0}));
    config.set_key_value("print_sequence",               new ConfigOptionEnum<PrintSequence>(PrintSequence::ByObject));

    // 2 extruders, 2 filaments
    config.set_key_value("nozzle_diameter",          new ConfigOptionFloats({0.4, 0.4}));
    config.set_key_value("printer_extruder_id",      new ConfigOptionInts({1, 2}));
    config.set_key_value("printer_extruder_variant", new ConfigOptionStrings({"Direct Drive Standard", "Direct Drive Standard"}));
    config.set_key_value("filament_diameter",        new ConfigOptionFloats({1.75, 1.75}));
    config.set_key_value("filament_colour",          new ConfigOptionStrings({"#FF0000", "#00FF00"}));
    config.set_key_value("filament_type",            new ConfigOptionStrings({"PLA", "PLA"}));
    config.option<ConfigOptionEnum<FilamentMapMode>>("filament_map_mode", true)->value = fmmManual;
    config.set_key_value("filament_map",             new ConfigOptionInts({1, 2}));
    config.set_key_value("default_filament_colour",  new ConfigOptionStrings({"#FF0000", "#00FF00"}));
    config.set_key_value("nozzle_temperature",       new ConfigOptionInts({210, 210}));
    config.set_key_value("nozzle_temperature_range_low",  new ConfigOptionInts({190, 190}));
    config.set_key_value("nozzle_temperature_range_high", new ConfigOptionInts({240, 240}));
    config.set_key_value("flush_multiplier",     new ConfigOptionFloats({1}));
    config.set_key_value("flush_volumes_matrix", new ConfigOptionFloats({0, 0, 0, 0}));

    // --- Machine limits (stride-2: e0_n, e0_s, e1_n, e1_s) ---
    // Extruder 0 has LOW limits, Extruder 1 has HIGH limits.
    config.set_key_value("machine_max_acceleration_x",          new ConfigOptionFloats({500, 0, 1000, 0}));
    config.set_key_value("machine_max_acceleration_y",          new ConfigOptionFloats({500, 0, 1000, 0}));
    config.set_key_value("machine_max_acceleration_z",          new ConfigOptionFloats({100, 0, 200, 0}));
    config.set_key_value("machine_max_acceleration_e",          new ConfigOptionFloats({5000, 0, 5000, 0}));
    config.set_key_value("machine_max_acceleration_extruding",  new ConfigOptionFloats({500, 0, 2000, 0}));
    config.set_key_value("machine_max_acceleration_retracting", new ConfigOptionFloats({600, 0, 2000, 0}));
    config.set_key_value("machine_max_acceleration_travel",     new ConfigOptionFloats({700, 0, 2500, 0}));
    config.set_key_value("machine_max_speed_x",  new ConfigOptionFloats({100, 0, 200, 0}));
    config.set_key_value("machine_max_speed_y",  new ConfigOptionFloats({110, 0, 210, 0}));
    config.set_key_value("machine_max_speed_z",  new ConfigOptionFloats({10, 0, 30, 0}));
    config.set_key_value("machine_max_speed_e",  new ConfigOptionFloats({50, 0, 80, 0}));
    config.set_key_value("machine_max_jerk_x",   new ConfigOptionFloats({5, 0, 15, 0}));
    config.set_key_value("machine_max_jerk_y",   new ConfigOptionFloats({6, 0, 16, 0}));
    config.set_key_value("machine_max_jerk_z",   new ConfigOptionFloats({0.4, 0, 0.8, 0}));
    config.set_key_value("machine_max_jerk_e",   new ConfigOptionFloats({3, 0, 8, 0}));
    config.set_key_value("machine_max_junction_deviation", new ConfigOptionFloats({0.02, 0, 0.08, 0}));

    // --- Print acceleration: 1500 mm/s² ---
    // Exceeds extruder 0's limit (500) → should be clamped to 500.
    // Does NOT exceed extruder 1's limit (2000) → passes through as 1500.
    config.set_key_value("default_acceleration",        new ConfigOptionFloats({1500, 1500}));
    config.set_key_value("outer_wall_acceleration",     new ConfigOptionFloats({1500, 1500}));
    config.set_key_value("inner_wall_acceleration",     new ConfigOptionFloats({1500, 1500}));
    config.set_key_value("top_surface_acceleration",    new ConfigOptionFloats({1500, 1500}));
    config.set_key_value("initial_layer_acceleration",  new ConfigOptionFloats({1500, 1500}));
    config.set_key_value("travel_acceleration",         new ConfigOptionFloats({1500, 1500}));

    // Model: two objects assigned to different extruders
    Model model;
    auto* obj1 = model.add_object();
    obj1->add_volume(cube(20));
    obj1->add_instance();

    auto* obj2 = model.add_object();
    obj2->add_volume(cube(20));
    obj2->add_instance();
    obj2->config.set_key_value("extruder", new ConfigOptionInt(2)); // 0-based index 1

    Print print;
    arrange_objects_on_test_bed(model, config);
    for (auto* mo : model.objects) {
        mo->ensure_on_bed();
        print.auto_assign_extruders(mo);
    }

    print.apply(model, config);
    print.validate();
    print.set_status_silent();
    print.process();

    std::string gcode = Slic3r::Test::gcode(print);

    SECTION("Preamble: max limit among used extruders") {
        THEN("M201 uses max (extruder 1's) acceleration values") {
            REQUIRE(gcode.find("M201 X1000 Y1000 Z200 E5000") != std::string::npos);
        }
        THEN("M204 uses max extruding/retracting/travel") {
            REQUIRE(gcode.find("M204 P2000 R2000 T2500") != std::string::npos);
        }
        THEN("M205 uses max jerk values") {
            REQUIRE(gcode.find("M205 X15.00 Y16.00 Z0.80 E8.00") != std::string::npos);
        }
    }

    SECTION("Preamble: EXTRUDER_LIMIT falls back to max when no filament is active") {
        // set_junction_deviation() is called during preamble with no active filament.
        // EXTRUDER_LIMIT(m_max_junction_deviation) → filament() == nullptr → max of all (0.08).
        THEN("M205 J uses max junction deviation") {
            REQUIRE(gcode.find("M205 J0.080") != std::string::npos);
        }
    }

    SECTION("Print: extruder 0 acceleration clamped to its specific limit") {
        // Extruder 0 machine limit = 500. Print accel = 1500 > 500 → clamped to 500.
        THEN("M204 P500 appears (extruder 0 clamped)") {
            REQUIRE(gcode.find("M204 P500") != std::string::npos);
        }
        THEN("M204 T700 appears (extruder 0 travel clamped)") {
            REQUIRE(gcode.find("M204 T700") != std::string::npos);
        }
    }

    SECTION("Print: extruder 1 acceleration NOT clamped to extruder 0's limit") {
        // Extruder 1 machine limit = 2000. Print accel = 1500 < 2000 → not clamped.
        THEN("M204 P1500 appears (extruder 1 not clamped to 500)") {
            REQUIRE(gcode.find("M204 P1500") != std::string::npos);
        }
    }
}

TEST_CASE("Acceleration and velocity limit commands print their values in general notation", "[GCodeWriter]")
{
    enum class Command { Print, Travel, KlipperLimits };
    struct Case
    {
        GCodeFlavor              flavor;
        Command                  command;
        unsigned int             acceleration;
        double                   jerk;
        bool                     comments;
        std::vector<std::string> present;
        std::vector<std::string> absent;
    };
    // accel_to_decel_factor is 50%, so ACCEL_TO_DECEL is half the acceleration.
    const Case c = GENERATE(values<Case>({
        {gcfKlipper, Command::KlipperLimits, 2000000, 25. / 3., false,
         {"SET_VELOCITY_LIMIT ACCEL=2000000 ", "ACCEL_TO_DECEL=1e+06 ", "SQUARE_CORNER_VELOCITY=8.33333\n"}, {}},
        {gcfKlipper, Command::KlipperLimits, 12345, 0., false, {"ACCEL=12345 ", "ACCEL_TO_DECEL=6172.5\n"}, {"SQUARE_CORNER_VELOCITY"}},
        {gcfKlipper, Command::KlipperLimits, 0, 0.25, true, {"SQUARE_CORNER_VELOCITY=0.25 ", "; adjust VELOCITY_LIMIT"}, {"ACCEL"}},
        {gcfKlipper, Command::Print, 3001, 0., true, {"ACCEL=3001 ", "ACCEL_TO_DECEL=1500.5 ", "; adjust ACCEL_TO_DECEL", "; adjust acceleration"}, {}},
        {gcfMarlinFirmware, Command::Print, 2500, 0., false, {"M204 P2500\n"}, {}},
        {gcfMarlinFirmware, Command::Travel, 7000, 0., false, {"M204 T7000\n"}, {}},
        {gcfRepRapFirmware, Command::Travel, 7000, 0., true, {"M204 T7000 ", "; adjust acceleration"}, {}},
        {gcfMarlinLegacy, Command::Print, 2500, 0., false, {"M204 S2500\n"}, {}},
        {gcfRepetier, Command::Print, 2500, 0., false, {"M201 X2500 Y2500\n"}, {}},
        {gcfRepetier, Command::Travel, 7000, 0., false, {"M202 X7000 Y7000\n"}, {}},
    }));

    struct CommentGuard
    {
        bool saved = GCodeWriter::full_gcode_comment;
        ~CommentGuard() { GCodeWriter::full_gcode_comment = saved; }
    } comment_guard;
    GCodeWriter::full_gcode_comment = c.comments;

    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_key_value("gcode_flavor", new ConfigOptionEnum<GCodeFlavor>(c.flavor));
    config.set_key_value("accel_to_decel_enable", new ConfigOptionBool(true));
    config.set_key_value("accel_to_decel_factor", new ConfigOptionPercent(50));
    for (const char *limit : {"machine_max_acceleration_extruding", "machine_max_acceleration_travel", "machine_max_acceleration_x",
                              "machine_max_acceleration_y", "machine_max_jerk_x", "machine_max_jerk_y"}) {
        std::vector<double> &values = config.option<ConfigOptionFloats>(limit)->values;
        std::fill(values.begin(), values.end(), 0.);
    }
    PrintConfig print_config;
    print_config.apply(config, true);
    GCodeWriter writer;
    writer.apply_print_config(print_config);

    const std::string line = c.command == Command::Print         ? writer.set_print_acceleration(c.acceleration) :
                             c.command == Command::Travel        ? writer.set_travel_acceleration(c.acceleration) :
                                                                   writer.set_accel_and_jerk(c.acceleration, c.jerk);
    INFO(line);
    for (const std::string &token : c.present)
        CHECK_THAT(line, Catch::Matchers::ContainsSubstring(token));
    for (const std::string &token : c.absent)
        CHECK_THAT(line, !Catch::Matchers::ContainsSubstring(token));
}

SCENARIO("Extruder reads the injected config column", "[GCodeWriter][H2C]") {
    GIVEN("A writer whose per-variant arrays hold three columns for two filaments") {
        GCodeWriter writer;
        // Column layout after a migrating regroup: filament 0 -> column 0, filament 1 ->
        // columns 1 (its first variant) and 2 (its second variant).
        writer.config.retraction_length.values   = {0.8, 0.5, 1.2};
        writer.config.z_hop.values               = {0.4, 0.6, 0.9};
        writer.config.retraction_speed.values    = {30., 40., 50.};
        writer.config.filament_flow_ratio.values = {0.98, 1.0, 1.02};
        // Filament-indexed arrays keep one entry per filament.
        writer.config.filament_diameter.values   = {1.75, 1.75};
        writer.set_extruders({0, 1});
        writer.toolchange(1, 1);
        Extruder *fil = writer.filament();
        REQUIRE(fil != nullptr);
        REQUIRE(fil->id() == 1);
        const double crossection = 1.75 * 1.75 * 0.25 * PI;

        WHEN("no column has been injected") {
            THEN("the getters read the filament id's column") {
                REQUIRE(fil->config_index() == 1);
                REQUIRE_THAT(fil->retraction_length(), Catch::Matchers::WithinAbs(0.5, 1e-9));
                REQUIRE_THAT(fil->retract_lift(), Catch::Matchers::WithinAbs(0.6, 1e-9));
                REQUIRE(fil->retract_speed() == 40);
                REQUIRE_THAT(fil->e_per_mm3(), Catch::Matchers::WithinRel(1.0 / crossection, 1e-9));
            }
        }
        WHEN("the second variant column is injected") {
            fil->set_config_index(2);
            THEN("the getters follow the column and the flow cache is rescaled") {
                REQUIRE(fil->config_index() == 2);
                REQUIRE_THAT(fil->retraction_length(), Catch::Matchers::WithinAbs(1.2, 1e-9));
                REQUIRE_THAT(fil->retract_lift(), Catch::Matchers::WithinAbs(0.9, 1e-9));
                REQUIRE(fil->retract_speed() == 50);
                REQUIRE_THAT(fil->e_per_mm3(), Catch::Matchers::WithinRel(1.02 / crossection, 1e-9));
            }
            THEN("filament-indexed reads keep using the filament id") {
                REQUIRE_THAT(fil->filament_diameter(), Catch::Matchers::WithinAbs(1.75, 1e-9));
            }
        }
        WHEN("a negative index is injected") {
            fil->set_config_index(2);
            fil->set_config_index(-1);
            THEN("resolution resets to the filament id") {
                REQUIRE(fil->config_index() == 1);
                REQUIRE_THAT(fil->retraction_length(), Catch::Matchers::WithinAbs(0.5, 1e-9));
                REQUIRE_THAT(fil->e_per_mm3(), Catch::Matchers::WithinRel(1.0 / crossection, 1e-9));
            }
        }
    }
}

// Numeric argument of every line starting with `prefix`, in file order.
static std::vector<int> collect_line_args(const std::string &gcode, const std::string &prefix)
{
    std::vector<int> values;
    std::istringstream stream(gcode);
    std::string line;
    while (std::getline(stream, line))
        if (line.compare(0, prefix.size(), prefix) == 0)
            values.push_back(std::atoi(line.c_str() + int(prefix.size())));
    return values;
}

static int count_lines_with_prefix(const std::string &gcode, const std::string &prefix)
{
    return (int) collect_line_args(gcode, prefix).size();
}

// A toolchange ordinal sequence is healthy when it advances by exactly one per
// change block; a change-less prime-tower visit must not consume an ordinal.
static bool ordinals_consecutive(const std::vector<int> &values)
{
    for (size_t i = 1; i < values.size(); ++i)
        if (values[i] != values[i - 1] + 1)
            return false;
    return true;
}

SCENARIO("Toolchange emission and prefix per printer kind", "[GCodeWriter][H2C]") {
    GIVEN("A dual-extruder writer with two filaments") {
        GCodeWriter writer;
        writer.config.filament_diameter.values = {1.75, 1.75};
        writer.set_extruders({0, 1});

        WHEN("the printer is a BBL machine") {
            writer.set_is_bbl_machine(true);
            THEN("the toolchange prefix is the plain T command") {
                REQUIRE_THAT(writer.toolchange_prefix(), Catch::Matchers::Equals("T"));
            }
            THEN("toolchange emits a single M1020 with the nozzle id") {
                const std::string gcode = writer.toolchange(1, 0);
                REQUIRE_THAT(gcode, Catch::Matchers::ContainsSubstring("M1020 S1 H0"));
                REQUIRE_THAT(gcode, !Catch::Matchers::StartsWith("T1"));
            }
            THEN("the other filament and nozzle emit their own ids") {
                REQUIRE_THAT(writer.toolchange(0, 1), Catch::Matchers::ContainsSubstring("M1020 S0 H1"));
            }
            THEN("an unresolved nozzle keeps the literal -1 convention") {
                REQUIRE_THAT(writer.toolchange(1, -1), Catch::Matchers::ContainsSubstring("M1020 S1 H-1"));
            }
        }
        WHEN("the printer is a BBL machine with manual filament change") {
            writer.set_is_bbl_machine(true);
            writer.config.manual_filament_change.value = true;
            THEN("the manual tag wins over the M1020 form") {
                REQUIRE_THAT(writer.toolchange_prefix(), Catch::Matchers::StartsWith(";"));
                const std::string gcode = writer.toolchange(1, 0);
                REQUIRE_THAT(gcode, Catch::Matchers::ContainsSubstring(writer.toolchange_prefix() + "1"));
                REQUIRE_THAT(gcode, !Catch::Matchers::ContainsSubstring("M1020"));
            }
        }
        WHEN("the printer is not a BBL machine") {
            THEN("toolchange keeps the plain T command") {
                REQUIRE_THAT(writer.toolchange_prefix(), Catch::Matchers::Equals("T"));
                const std::string gcode = writer.toolchange(1, 0);
                REQUIRE_THAT(gcode, Catch::Matchers::StartsWith("T1"));
                REQUIRE_THAT(gcode, !Catch::Matchers::ContainsSubstring("M1020"));
            }
        }
    }
}

// Shared dual-extruder printer config for the toolchange-count scenarios below.
static DynamicPrintConfig dual_extruder_toolchange_config()
{
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_key_value("gcode_flavor",               new ConfigOptionEnum<GCodeFlavor>(gcfMarlinFirmware));
    config.set_key_value("emit_machine_limits_to_gcode", new ConfigOptionBool(false));
    config.set_key_value("machine_start_gcode",        new ConfigOptionString(""));
    config.set_key_value("layer_height",               new ConfigOptionFloat(0.2));
    config.set_key_value("initial_layer_print_height", new ConfigOptionFloat(0.2));
    config.set_key_value("initial_layer_line_width",   new ConfigOptionFloatOrPercent(0, false));
    config.set_key_value("z_hop",                      new ConfigOptionFloats({0., 0.}));
    // The change block carries both a real toolchange command and the ordinal
    // placeholder the stock profiles feed to the firmware.
    config.set_key_value("change_filament_gcode",
                         new ConfigOptionString("T[next_filament_id]\nM620 O{toolchange_count + 1}\n"));

    // 2 extruders, one filament each (manual map so nothing regroups them).
    config.set_key_value("nozzle_diameter",          new ConfigOptionFloats({0.4, 0.4}));
    config.set_key_value("printer_extruder_id",      new ConfigOptionInts({1, 2}));
    config.set_key_value("printer_extruder_variant", new ConfigOptionStrings({"Direct Drive Standard", "Direct Drive Standard"}));
    config.set_key_value("filament_diameter",        new ConfigOptionFloats({1.75, 1.75}));
    config.set_key_value("filament_colour",          new ConfigOptionStrings({"#FF0000", "#00FF00"}));
    config.set_key_value("default_filament_colour",  new ConfigOptionStrings({"#FF0000", "#00FF00"}));
    config.set_key_value("filament_type",            new ConfigOptionStrings({"PLA", "PLA"}));
    config.option<ConfigOptionEnum<FilamentMapMode>>("filament_map_mode", true)->value = fmmManual;
    config.set_key_value("filament_map",             new ConfigOptionInts({1, 2}));
    config.set_key_value("nozzle_temperature",       new ConfigOptionInts({210, 210}));
    config.set_key_value("nozzle_temperature_range_low",  new ConfigOptionInts({190, 190}));
    config.set_key_value("nozzle_temperature_range_high", new ConfigOptionInts({240, 240}));
    config.set_key_value("flush_multiplier",     new ConfigOptionFloats({1}));
    config.set_key_value("flush_volumes_matrix", new ConfigOptionFloats({0, 140, 140, 0}));
    // Inside the 200x200 test bed; the default y, 220, is not, and generation rejects that.
    config.set_key_value("wipe_tower_x",         new ConfigOptionFloats({50.}));
    config.set_key_value("wipe_tower_y",         new ConfigOptionFloats({50.}));
    return config;
}

SCENARIO("Change blocks carry consecutive toolchange ordinals without a duplicate command", "[GCodeWriter][H2C]") {
    GIVEN("Two sequentially printed objects on different extruders of a BBL machine") {
        DynamicPrintConfig config = dual_extruder_toolchange_config();
        config.set_key_value("print_sequence", new ConfigOptionEnum<PrintSequence>(PrintSequence::ByObject));

        Model model;
        auto *obj1 = model.add_object();
        obj1->add_volume(cube(20));
        obj1->add_instance();
        auto *obj2 = model.add_object();
        obj2->add_volume(cube(20));
        obj2->add_instance();
        obj2->config.set_key_value("extruder", new ConfigOptionInt(2));

        auto slice_to_gcode = [&]() {
            Print print;
            print.is_BBL_printer() = true;
            arrange_objects_on_test_bed(model, config);
            for (auto *mo : model.objects) {
                mo->ensure_on_bed();
                print.auto_assign_extruders(mo);
            }
            print.apply(model, config);
            print.validate();
            print.set_status_silent();
            print.process();
            return Slic3r::Test::gcode(print);
        };

        WHEN("the change block already changes the tool") {
            const std::string gcode = slice_to_gcode();
            const std::vector<int> ordinals = collect_line_args(gcode, "M620 O");
            THEN("each change block advances the ordinal by exactly one, without inflation") {
                REQUIRE(!ordinals.empty());
                REQUIRE(ordinals_consecutive(ordinals));
                REQUIRE(ordinals.front() <= 3);
            }
            THEN("the writer's own command is suppressed as a duplicate") {
                REQUIRE(count_lines_with_prefix(gcode, "M1020") == 0);
                REQUIRE(count_lines_with_prefix(gcode, "T1") >= 1);
            }
        }
        WHEN("the change block does not change the tool itself") {
            config.set_key_value("change_filament_gcode",
                                 new ConfigOptionString("M620 O{toolchange_count + 1}\n"));
            const std::string gcode = slice_to_gcode();
            const std::vector<int> ordinals = collect_line_args(gcode, "M620 O");
            THEN("the writer's toolchange survives and carries a nozzle id") {
                REQUIRE(count_lines_with_prefix(gcode, "M1020 S1 H") >= 1);
            }
            THEN("the ordinal sequence stays consecutive") {
                REQUIRE(!ordinals.empty());
                REQUIRE(ordinals_consecutive(ordinals));
                REQUIRE(ordinals.front() <= 3);
            }
        }
    }
}

SCENARIO("Prime-tower visits without a filament change do not advance the toolchange ordinal", "[GCodeWriter][H2C]") {
    GIVEN("A print whose only filament change happens far above the bed") {
        DynamicPrintConfig config = dual_extruder_toolchange_config();
        config.set_key_value("enable_prime_tower", new ConfigOptionBool(true));

        // Filament 2 is used only above z=6, so every tower layer below it is a
        // change-less visit — the exact geometry that used to inflate the ordinal.
        Model model;
        auto *obj = model.add_object();
        obj->add_volume(cube(10));
        obj->add_instance();
        DynamicPrintConfig range_config;
        range_config.set_key_value("extruder", new ConfigOptionInt(2));
        // Every layer range must carry a layer_height (see layer_height_profile_from_ranges).
        range_config.set_key_value("layer_height", new ConfigOptionFloat(0.2));
        obj->layer_config_ranges[{6.0, 10.0}].assign_config(std::move(range_config));

        Print print;
        print.is_BBL_printer() = true;
        arrange_objects_on_test_bed(model, config);
        for (auto *mo : model.objects) {
            mo->ensure_on_bed();
            print.auto_assign_extruders(mo);
        }
        print.apply(model, config);
        print.validate();
        print.set_status_silent();
        print.process();
        const std::string gcode = Slic3r::Test::gcode(print);

        WHEN("the print is exported") {
            const std::vector<int> ordinals = collect_line_args(gcode, "M620 O");
            THEN("the prime-tower toolchange path was exercised") {
                REQUIRE_THAT(gcode, Catch::Matchers::ContainsSubstring("CP TOOLCHANGE START"));
            }
            THEN("dozens of change-less tower layers consume no ordinal") {
                REQUIRE(!ordinals.empty());
                REQUIRE(ordinals_consecutive(ordinals));
                REQUIRE(ordinals.front() <= 3);
            }
            THEN("no duplicate toolchange command follows the change block") {
                REQUIRE(count_lines_with_prefix(gcode, "M1020") == 0);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Real-profile toolchange coverage, targeted. The all-vendors sweep in test_profile_slicing.cpp now
// slices a two-colour cube per printer, so it already expands every shipped change_filament_gcode with
// each printer's DEFAULT extruder variants (both the single-nozzle append_tcr and dual-nozzle set_extruder
// paths). What that sweep can't reach is a variant-conditional branch the defaults never select — H2D's
// change gcode has an `== "Direct Drive TPU High Flow"` block. This scenario forces that branch by handing
// the extruders distinct kits, so an unregistered placeholder inside it still throws "Variable does not
// exist" here instead of only in the field.
// ---------------------------------------------------------------------------

// Two 20mm cubes on separate extruders of a BBL machine, printed by object so exactly
// one real toolchange fires and drives the change_filament_gcode. Returns the g-code.
static std::string slice_two_object_bbl(DynamicPrintConfig &config)
{
    config.set_key_value("print_sequence", new ConfigOptionEnum<PrintSequence>(PrintSequence::ByObject));

    Model model;
    auto *obj1 = model.add_object();
    obj1->add_volume(cube(20));
    obj1->add_instance();
    auto *obj2 = model.add_object();
    obj2->add_volume(cube(20));
    obj2->add_instance();
    obj2->config.set_key_value("extruder", new ConfigOptionInt(2));

    Print print;
    print.is_BBL_printer() = true;
    arrange_objects_on_test_bed(model, config);
    for (auto *mo : model.objects) {
        mo->ensure_on_bed();
        print.auto_assign_extruders(mo);
    }
    print.apply(model, config);
    print.validate();
    print.set_status_silent();
    print.process();
    return Slic3r::Test::gcode(print);
}

// The real change_filament_gcode of a shipped "<printer> 0.4 nozzle" machine profile.
// The profile does not state it inline any more: it names the template carrying it in
// `include`, and the loader layers that template under the preset. Follow the same list
// - so a profile that stops naming one fails here instead of silently slicing G-code it
// no longer ships.
static std::string shipped_change_filament_gcode(const std::string &printer)
{
    const std::string machine_dir   = std::string(PROFILES_DIR) + "/BBL/machine/";
    const std::string machine_path  = machine_dir + "Bambu Lab " + printer + " 0.4 nozzle.json";
    const std::string template_name = "Bambu Lab " + printer + " 0.4 nozzle template change_filament_gcode";
    // PROFILES_DIR is an absolute path baked in at build time; a sparse test checkout
    // without resources/ leaves it missing. Skip rather than dereference a config that
    // never loaded - this is the only fff_print test that reads a shipped profile.
    if (!boost::filesystem::exists(machine_path))
        SKIP("shipped profile not present in this checkout: " << machine_path);

    auto load = [](const std::string &file, std::map<std::string, std::string> &key_values) {
        DynamicPrintConfig        config;
        ConfigSubstitutionContext substitutions{ForwardCompatibilitySubstitutionRule::Enable};
        std::string               reason;
        // false: `inherits` and `include` stay out of the config and land in key_values,
        // for the loader to resolve - neither is a slicing setting.
        config.load_from_json(file, substitutions, false, key_values, reason);
        // Fail loudly on a malformed/renamed profile instead of null-dereferencing in opt_string.
        INFO("profile: " << file << (reason.empty() ? "" : ("  load reason: " + reason)));
        return config;
    };

    std::map<std::string, std::string> machine_values;
    load(machine_path, machine_values);
    REQUIRE(machine_values.count("include") == 1);
    const nlohmann::json includes = nlohmann::json::parse(machine_values["include"]);
    REQUIRE(std::find(includes.begin(), includes.end(), nlohmann::json(template_name)) != includes.end());

    std::map<std::string, std::string> template_values;
    const DynamicPrintConfig           config = load(machine_dir + template_name + ".json", template_values);
    REQUIRE(config.has("change_filament_gcode"));
    return config.opt_string("change_filament_gcode");
}

SCENARIO("Toolchange gcode resolves old/new_extruder_variant from printer_extruder_variant", "[GCodeWriter][H2C]")
{
    GIVEN("a BBL dual-extruder print whose change gcode reads the extruder-variant placeholders") {
        DynamicPrintConfig config = dual_extruder_toolchange_config();
        // A distinctive variant that can only reach the g-code through printer_extruder_variant.
        // Both entries carry it so the assertion is independent of which physical extruder the
        // emitted change routes through.
        config.set_key_value("printer_extruder_variant",
                             new ConfigOptionStrings({"Direct Drive TPU High Flow", "Direct Drive TPU High Flow"}));
        config.set_key_value("change_filament_gcode", new ConfigOptionString(
            "; VARIANT old={old_extruder_variant} new={new_extruder_variant}\nT[next_filament_id]\n"));

        WHEN("the print is sliced") {
            const std::string gcode = slice_two_object_bbl(config);
            THEN("both placeholders resolve to the printer_extruder_variant value") {
                // The resolved line is the proof: an unresolved token or a parser throw would
                // prevent this exact line from being emitted. (A negative "{token}" check is
                // unreliable — the g-code's trailing config dump echoes the raw template.)
                REQUIRE_THAT(gcode, Catch::Matchers::ContainsSubstring(
                    "; VARIANT old=Direct Drive TPU High Flow new=Direct Drive TPU High Flow"));
            }
        }
    }
}

SCENARIO("Global current-tool placeholders resolve in a context with no local injection", "[GCodeWriter][H2C]")
{
    GIVEN("a BBL dual-extruder print whose before_layer_change_gcode reads the current-tool placeholders") {
        DynamicPrintConfig config = dual_extruder_toolchange_config();
        // before_layer_change is one of the contexts that inject NO current_* into their local config
        // (unlike change_filament / machine_end / layer_change), so these placeholders can only resolve
        // through the GLOBAL parser vars published at each toolchange (and at the initial set_extruder).
        // Pre-fix current_filament_id / current_extruder_id / current_nozzle_id were undefined here and the
        // whole slice threw a PlaceholderParserError — the same failure mode X2D's layer_change hit.
        config.set_key_value("before_layer_change_gcode", new ConfigOptionString(
            "; GVAR fid={current_filament_id} eid={current_extruder_id} nid={current_nozzle_id}\n"));

        WHEN("the print is sliced (obj1 on filament 0, obj2 on filament 1)") {
            std::string gcode;
            REQUIRE_NOTHROW(gcode = slice_two_object_bbl(config));
            THEN("all three globals resolve to the CORRECT active-tool values on both sides of the change") {
                // Assert the FULL resolved marker, not just no-throw: obj1 prints on filament 0 (extruder 0,
                // nozzle 0) and obj2 on filament 1 (extruder 1, nozzle 1). Locking every field means a
                // stale or wrong global (e.g. obj2 still reading fid=0, or a mismatched extruder/nozzle id)
                // fails here — this is the value guard that replaces the old "throws on undefined" canary.
                // Only the emitted before_layer_change lines carry resolved values; the trailing config dump
                // keeps the raw "{current_filament_id}" template, so these are unambiguous.
                REQUIRE_THAT(gcode, Catch::Matchers::ContainsSubstring("; GVAR fid=0 eid=0 nid=0"));
                REQUIRE_THAT(gcode, Catch::Matchers::ContainsSubstring("; GVAR fid=1 eid=1 nid=1"));
            }
        }
    }
}

SCENARIO("Shipped dual-nozzle change_filament_gcode resolves during a real slice", "[GCodeWriter][H2C][Profiles]")
{
    const std::string printer = GENERATE(std::string("H2C"), std::string("H2D"), std::string("H2D Pro"), std::string("X2D"));

    GIVEN("the real " + printer + " change_filament_gcode driving a BBL dual-extruder slice") {
        DynamicPrintConfig config = dual_extruder_toolchange_config();
        config.set_key_value("change_filament_gcode", new ConfigOptionString(shipped_change_filament_gcode(printer)));
        // H2D's gcode branches on the extruder variant; give the extruders distinct kits so the
        // "Direct Drive TPU High Flow" branch is reachable.
        config.set_key_value("printer_extruder_variant",
                             new ConfigOptionStrings({"Direct Drive Standard", "Direct Drive TPU High Flow"}));
        // Extruder-indexed machine rates the stock gcode divides by (default size 1); size to 2 extruders.
        config.set_key_value("hotend_cooling_rate", new ConfigOptionFloatsNullable({2.0, 2.0}));
        config.set_key_value("hotend_heating_rate", new ConfigOptionFloatsNullable({2.0, 2.0}));

        THEN("every placeholder resolves (no undefined-variable throw) and the change block runs") {
            std::string gcode;
            REQUIRE_NOTHROW(gcode = slice_two_object_bbl(config));
            // A resolved marker only the emitted change block produces (the trailing config
            // dump keeps the raw "{filament_type[...]}" template), so this confirms the real
            // change_filament_gcode was expanded, not merely echoed.
            REQUIRE_THAT(gcode, Catch::Matchers::ContainsSubstring("set_filament_type:PLA"));
        }
    }
}


SCENARIO("set_pressure_advance emits nothing for negative PA", "[GCodeWriter][PressureAdvance]") {
    GIVEN("A default GCodeWriter") {
        GCodeWriter writer;
        THEN("Negative PA returns empty regardless of firmware flavor") {
            writer.config.gcode_flavor.value = gcfKlipper;
            REQUIRE(writer.set_pressure_advance(-1.0).empty());
            writer.config.gcode_flavor.value = gcfRepRapFirmware;
            REQUIRE(writer.set_pressure_advance(-0.001).empty());
            writer.config.gcode_flavor.value = gcfMarlinFirmware;
            REQUIRE(writer.set_pressure_advance(-100.0).empty());
        }
    }
}

SCENARIO("set_pressure_advance emits Klipper form with optional EXTRUDER=extruder<N>", "[GCodeWriter][PressureAdvance]") {
    GIVEN("A Klipper-flavored GCodeWriter") {
        GCodeWriter writer;
        writer.config.gcode_flavor.value = gcfKlipper;

        WHEN("set_pressure_advance is called without a tool index") {
            std::string out = writer.set_pressure_advance(0.05);
            THEN("Output contains SET_PRESSURE_ADVANCE ADVANCE=0.05 with no EXTRUDER qualifier") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("SET_PRESSURE_ADVANCE"));
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("ADVANCE=0.05; Override"));
                REQUIRE_THAT(out, !Catch::Matchers::ContainsSubstring("EXTRUDER="));
            }
        }
        WHEN("set_pressure_advance is called with tool=0") {
            std::string out = writer.set_pressure_advance(0.04, 0);
            THEN("Output targets EXTRUDER=extruder (no trailing index)") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("EXTRUDER=extruder;"));
                REQUIRE_THAT(out, !Catch::Matchers::ContainsSubstring("EXTRUDER=extruder0"));
            }
        }
        WHEN("set_pressure_advance is called with tool=3") {
            std::string out = writer.set_pressure_advance(0.06, 3);
            THEN("Output targets EXTRUDER=extruder3") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("EXTRUDER=extruder3"));
            }
        }
    }
}

SCENARIO("set_pressure_advance emits RepRapFirmware form with D<N>", "[GCodeWriter][PressureAdvance]") {
    GIVEN("An RRF-flavored GCodeWriter") {
        GCodeWriter writer;
        writer.config.gcode_flavor.value = gcfRepRapFirmware;

        WHEN("set_pressure_advance is called without a tool index") {
            std::string out = writer.set_pressure_advance(0.07);
            THEN("Output keeps the historical D0 rather than depending on the selected tool") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M572 D0 S0.07"));
            }
        }
        WHEN("set_pressure_advance is called with tool=0") {
            std::string out = writer.set_pressure_advance(0.08, 0);
            THEN("Output contains D0, same as the tool-less form") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M572 D0 S0.08"));
            }
        }
        WHEN("set_pressure_advance is called with tool=2") {
            std::string out = writer.set_pressure_advance(0.09, 2);
            THEN("Output contains D2") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M572 D2 S0.09"));
            }
        }
    }
}

SCENARIO("set_pressure_advance emits Marlin 2.x form with optional T<N>", "[GCodeWriter][PressureAdvance]") {
    GIVEN("A Marlin 2-flavored GCodeWriter") {
        GCodeWriter writer;
        writer.config.gcode_flavor.value = gcfMarlinFirmware;

        WHEN("set_pressure_advance is called without a tool index") {
            std::string out = writer.set_pressure_advance(0.10);
            THEN("Output is bare M900 K... with no T qualifier") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M900 K0.1; Override"));
                REQUIRE_THAT(out, !Catch::Matchers::ContainsSubstring(" T"));
            }
        }
        WHEN("set_pressure_advance is called with tool=0") {
            std::string out = writer.set_pressure_advance(0.11, 0);
            THEN("Output contains T0") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M900 K0.11 T0"));
            }
        }
        WHEN("set_pressure_advance is called with tool=1") {
            std::string out = writer.set_pressure_advance(0.12, 1);
            THEN("Output contains T1") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M900 K0.12 T1"));
            }
        }
    }
}

SCENARIO("set_pressure_advance emits Marlin Legacy form without tool qualifier even when tool index is supplied",
         "[GCodeWriter][PressureAdvance]") {
    GIVEN("A Marlin Legacy-flavored GCodeWriter") {
        GCodeWriter writer;
        writer.config.gcode_flavor.value = gcfMarlinLegacy;

        WHEN("set_pressure_advance is called without a tool index") {
            std::string out = writer.set_pressure_advance(0.05);
            THEN("Output is bare M900 K...") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M900 K0.05"));
                REQUIRE_THAT(out, !Catch::Matchers::ContainsSubstring(" T"));
            }
        }
        WHEN("set_pressure_advance is called with tool=2 (a hypothetical IMEX secondary)") {
            std::string out = writer.set_pressure_advance(0.06, 2);
            THEN("Output is still bare M900 — Marlin Legacy has no per-tool LA") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M900 K0.06"));
                REQUIRE_THAT(out, !Catch::Matchers::ContainsSubstring(" T2"));
            }
        }
    }
}

SCENARIO("set_temperature per-flavor command routing", "[GCodeWriter][Temperature]") {
    GIVEN("temperature=210, no tool index, no wait") {
        WHEN("flavor is Marlin 2") {
            std::string out = GCodeWriter::set_temperature(210, gcfMarlinFirmware, false, -1, std::string());
            THEN("output is M104 S210 with no tool qualifier") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M104 S210"));
                REQUIRE_THAT(out, !Catch::Matchers::ContainsSubstring(" T"));
                REQUIRE_THAT(out, !Catch::Matchers::ContainsSubstring("M109"));
            }
        }
        WHEN("flavor is RepRapFirmware") {
            std::string out = GCodeWriter::set_temperature(210, gcfRepRapFirmware, false, -1, std::string());
            THEN("output is G10 S210 (M104 is deprecated on RRF)") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("G10 S210"));
                REQUIRE_THAT(out, !Catch::Matchers::ContainsSubstring("M104"));
            }
        }
        WHEN("flavor is Mach3 or Machinekit") {
            std::string mach3      = GCodeWriter::set_temperature(210, gcfMach3,      false, -1, std::string());
            std::string machinekit = GCodeWriter::set_temperature(210, gcfMachinekit, false, -1, std::string());
            THEN("output uses P-prefix for the value instead of S") {
                REQUIRE_THAT(mach3,      Catch::Matchers::ContainsSubstring("M104 P210"));
                REQUIRE_THAT(machinekit, Catch::Matchers::ContainsSubstring("M104 P210"));
                REQUIRE_THAT(mach3,      !Catch::Matchers::ContainsSubstring("S210"));
            }
        }
    }
}

SCENARIO("set_temperature wait=true handling per firmware", "[GCodeWriter][Temperature]") {
    WHEN("flavor is Marlin 2 with wait") {
        std::string out = GCodeWriter::set_temperature(210, gcfMarlinFirmware, true, -1, std::string());
        THEN("output is M109 S210 (blocking wait)") {
            REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M109 S210"));
            REQUIRE_THAT(out, !Catch::Matchers::ContainsSubstring("M104"));
        }
    }
    WHEN("flavor is MakerWare or Sailfish with wait") {
        std::string mw = GCodeWriter::set_temperature(210, gcfMakerWare, true, -1, std::string());
        std::string sf = GCodeWriter::set_temperature(210, gcfSailfish,  true, -1, std::string());
        THEN("output is empty — these flavors don't support blocking waits") {
            REQUIRE(mw.empty());
            REQUIRE(sf.empty());
        }
    }
    WHEN("flavor is Teacup with wait") {
        std::string out = GCodeWriter::set_temperature(210, gcfTeacup, true, -1, std::string());
        THEN("output emits M104 + a separate M116 poll (Teacup doesn't support M109)") {
            REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M104 S210"));
            REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M116"));
            REQUIRE_THAT(out, !Catch::Matchers::ContainsSubstring("M109"));
        }
    }
    WHEN("flavor is RepRapFirmware with wait") {
        std::string out = GCodeWriter::set_temperature(210, gcfRepRapFirmware, true, -1, std::string());
        THEN("output emits G10 + M116 (same poll pattern as Teacup)") {
            REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("G10 S210"));
            REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M116"));
        }
    }
}

SCENARIO("set_temperature per-tool qualifier routing for IMEX secondary carriages",
         "[GCodeWriter][Temperature]") {
    // IMEX secondary tools never go through a tool-change, so layer-change temperature
    // for them is emitted via the tool-qualified static set_temperature overload.
    GIVEN("temperature=220, tool=2, no wait") {
        WHEN("flavor is Marlin 2") {
            std::string out = GCodeWriter::set_temperature(220, gcfMarlinFirmware, false, 2, std::string());
            THEN("output contains T2 qualifier") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M104 S220 T2"));
            }
        }
        WHEN("flavor is RepRapFirmware") {
            std::string out = GCodeWriter::set_temperature(220, gcfRepRapFirmware, false, 2, std::string());
            THEN("output uses P-prefix for tool (RRF convention), not T") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("G10 S220 P2"));
                REQUIRE_THAT(out, !Catch::Matchers::ContainsSubstring(" T2"));
            }
        }
        WHEN("flavor is Klipper") {
            std::string out = GCodeWriter::set_temperature(220, gcfKlipper, false, 1, std::string());
            THEN("output contains T1 qualifier (Klipper layer-change temperature uses T)") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M104 S220 T1"));
            }
        }
    }
}

SCENARIO("set_temperature instance overload forces tool=-1 on single-extruder writers",
         "[GCodeWriter][Temperature]") {
    // Guards against spuriously emitting `T0` on printers that only have one extruder.
    // The instance overload discards the tool argument when !multiple_extruders.
    GIVEN("A default GCodeWriter (multiple_extruders=false)") {
        GCodeWriter writer;
        writer.config.gcode_flavor.value = gcfMarlinFirmware;

        WHEN("set_temperature is called with tool=2") {
            std::string out = writer.set_temperature(210, false, 2);
            THEN("output has no T qualifier despite the caller passing tool=2") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M104 S210"));
                REQUIRE_THAT(out, !Catch::Matchers::ContainsSubstring(" T"));
            }
        }
    }
    GIVEN("A GCodeWriter with multiple_extruders=true (not SEMM)") {
        GCodeWriter writer;
        writer.config.gcode_flavor.value = gcfMarlinFirmware;
        writer.multiple_extruders = true;

        WHEN("set_temperature is called with tool=2") {
            std::string out = writer.set_temperature(210, false, 2);
            THEN("tool argument passes through — output contains T2") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M104 S210 T2"));
            }
        }
    }
}

SCENARIO("set_pressure_advance emits BBL M900 L1000 M10 regardless of tool index",
         "[GCodeWriter][PressureAdvance]") {
    GIVEN("A BBL-flagged GCodeWriter (the flag overrides firmware flavor routing)") {
        GCodeWriter writer;
        writer.set_is_bbl_machine(true);
        // Flavor intentionally set to something other than the BBL branch to prove the flag wins.
        writer.config.gcode_flavor.value = gcfMarlinFirmware;

        WHEN("set_pressure_advance is called without a tool index") {
            std::string out = writer.set_pressure_advance(0.05);
            THEN("Output is the BBL-specific M900 Kx L1000 M10 form") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M900 K0.05 L1000 M10"));
            }
        }
        WHEN("set_pressure_advance is called with a tool index") {
            std::string out = writer.set_pressure_advance(0.05, 2);
            THEN("BBL output is unchanged — no per-tool qualifier is emitted on BBL printers") {
                REQUIRE_THAT(out, Catch::Matchers::ContainsSubstring("M900 K0.05 L1000 M10"));
                REQUIRE_THAT(out, !Catch::Matchers::ContainsSubstring(" T2"));
                REQUIRE_THAT(out, !Catch::Matchers::ContainsSubstring("EXTRUDER="));
            }
        }
    }
}

TEST_CASE("Custom G-code motion limits are restored before generated moves", "[GCodeWriter]")
{
    const std::string gcode = Slic3r::Test::slice({ cube(20) }, {
        { "gcode_flavor",                "marlin" },
        { "gcode_comments",              "1" },
        { "machine_start_gcode",         "" },
        { "layer_change_gcode",          "M204 S5000\nm205 x5 y5\n" },
        { "layer_height",                "0.2" },
        { "initial_layer_print_height",  "0.2" },
        { "initial_layer_line_width",    "0" },
        { "z_hop",                       "0" },
        { "default_acceleration",        "6000" },
        { "initial_layer_acceleration",  "6000" },
        { "outer_wall_acceleration",     "6000" },
        { "inner_wall_acceleration",     "0" },
        { "default_jerk",                "8" },
        { "initial_layer_jerk",          "8" },
        { "outer_wall_jerk",             "8" },
        { "inner_wall_jerk",             "0" },
    });

    const size_t custom_gcode_pos = gcode.find("m205 x5 y5");
    REQUIRE(custom_gcode_pos != std::string::npos);
    REQUIRE(gcode.find("M204 S6000 ; adjust acceleration", custom_gcode_pos) != std::string::npos);
    REQUIRE(gcode.find("M205 X8 Y8 ; adjust jerk", custom_gcode_pos) != std::string::npos);
}

// Regression test for the belt-printer "illegal gantry move at print start" bug.
//
// On a belt printer the layer-change z-hop is deferred (lazy_lift) and consumed
// by the first travel_to_xyz, whose NormalLift branch lifts in place via
// _travel_to_z(). On a normal printer _travel_to_z emits a Z-only move, but in
// belt mode Z is coupled to Y/X, so _travel_to_z re-emits the current m_pos
// through the belt shear. At print start (and after custom gcode)
// is_current_position_clear() is false and m_pos.xy is still the uninitialised
// origin (0,0), which shears into machine (X=bed_max, Y=layer_z) — a move far up
// the gantry, e.g. "G1 X95 Y168.19 Z237.857". The fix guards that lift on
// is_current_position_clear(), mirroring the SlopeLift branch.
SCENARIO("Belt: the first travel does not lift through the uninitialised origin", "[GCodeWriter][belt]")
{
    GIVEN("A fresh belt-kinematics GCodeWriter configured for an X-tilt 45 degree belt") {
        // Machine-frame + slicer->world back-transform config (X tilt, 45 deg).
        PrintConfig belt_config;
        belt_config.belt_printer.value               = true;
        belt_config.belt_slice_rotation.value        = BeltRotationAxis::X;
        belt_config.belt_slice_rotation_angle.value  = 45.0;
        belt_config.belt_frame_tilt_decouple.value   = false;
        belt_config.belt_frame_tilt_angle.value      = 45.0;

        GCodeWriter writer;
        install_belt_kinematics(writer, belt_config);

        std::vector<unsigned int> extruder_ids { 0 };
        writer.set_extruders(extruder_ids);
        writer.set_extruder(0);
        // travel_speed became per-extruder (ConfigOptionFloatsNullable) upstream.
        writer.config.travel_speed.values       = { 100.0 };
        writer.config.z_hop.values              = { 0.4 };
        writer.config.retract_lift_above.values = { 0.0 };
        writer.config.retract_lift_below.values = { 0.0 };

        // A fresh writer has not established its planar position yet — this is the
        // precondition that made the origin leak into the first move.
        REQUIRE_FALSE(writer.is_current_position_clear());

        WHEN("a layer-change z-hop is pending and we travel to the first object point") {
            // Defer a z-hop, exactly as a retract on layer change leaves it.
            writer.lazy_lift(LiftType::NormalLift);

            // First object point in slicing coordinates: a near-belt point (y ~= -z)
            // so its transformed gantry Y is small (~1mm). The bogus origin lift, in
            // contrast, would shear to machine Y ~= nominal_z.
            const double nominal_z = 100.0;
            std::string gcode = writer.travel_to_xyz(Vec3d(10.0, -(nominal_z - 1.0), nominal_z));

            THEN("no emitted move flies up the gantry; machine Y stays near the part") {
                double max_y = std::numeric_limits<double>::lowest();
                GCodeReader reader;
                reader.parse_buffer(gcode, [&max_y](GCodeReader &, const GCodeReader::GCodeLine &line) {
                    if (line.cmd_is("G1") && line.has(Y))
                        max_y = std::max(max_y, double(line.y()));
                });
                // The destination shears to machine Y ~= 1mm. The old origin-lift bug
                // produced a separate move at machine Y ~= nominal_z (100mm), so any
                // Y well above the part means the origin leaked into a move.
                REQUIRE(max_y > 0.0);   // the destination move was emitted and parsed
                REQUIRE(max_y < 10.0);  // ... and nothing flew up the gantry
            }
        }
    }
}

// Regression test for the belt-printer "phantom extrusion line from Y=0" bug.
//
// GCodeProcessor::store_move_vertex pins a move's stored Z to the first-layer
// height while m_processing_start_custom_gcode is set (the start G-code "prepare"
// stage), because on a normal printer the toolhead Z there is not yet a real print
// height. On a belt printer that override is wrong: Z is written explicitly and the
// designed-view back-transform couples machine Z into the rendered model Y (the
// belt tilt mixes the height and belt-feed axes). Overriding it back-transforms the
// last prepare-stage move (the unretract right before the first extrusion) to
// model Y ~= 0, and libvgcode then draws a phantom extrusion segment from Y ~= 0 to
// the first real toolpath — rendered in the first extrusion role's color. The fix
// keeps the real Z for belt printers (gated on belt_tilt_angle). Here we assert the
// prepare-stage move keeps its real Z so it can no longer leak to Y ~= 0.
SCENARIO("Belt: start-gcode prepare-stage moves keep their real Z", "[GCode][belt]")
{
    // Belt printers are non-Bambu, so the G-code uses the "compatible" reserved
    // tags ("TYPE:" for the extrusion role). The processor selects the tag table
    // from the static s_IsBBLPrinter flag, so mirror the belt-printer setting here
    // (saved/restored so test ordering stays unaffected).
    struct BBLPrinterGuard {
        bool prev = GCodeProcessor::s_IsBBLPrinter;
        BBLPrinterGuard()  { GCodeProcessor::s_IsBBLPrinter = false; }
        ~BBLPrinterGuard() { GCodeProcessor::s_IsBBLPrinter = prev; }
    } bbl_guard;

    GIVEN("A belt G-code whose start sequence travels to a high machine Z before the first extrusion") {
        // The leading "; belt_slice_rotation_angle = 45" header sets belt_tilt_angle
        // (parsed before the body), enabling the belt code path. ;TYPE:Custom before
        // any G1 turns on the prepare stage; ;TYPE:Outer wall turns it off, exactly
        // as a sliced belt print is laid out.
        const std::string gcode =
            "; belt_slice_rotation_angle = 45\n"
            "G90\n"
            "G21\n"
            "M83\n"
            ";TYPE:Custom\n"
            "G1 E-1.5 F2100\n"            // retract at the (0,0,0) origin
            "G1 X45 Y0.3 Z50 F12000\n"   // travel to the approach point (prepare stage)
            "G1 E1.5 F1800\n"            // unretract in place (prepare stage)
            ";TYPE:Outer wall\n"
            "G1 X46 Y0.3 Z50 E0.05\n";   // first extrusion, same Z as the approach

        GCodeProcessor processor;
        processor.process_buffer(gcode);
        const GCodeProcessorResult& result = processor.get_result();

        THEN("the belt code path is active") {
            REQUIRE_THAT(result.belt_tilt_angle, Catch::Matchers::WithinAbs(45.0, 1e-4));
        }

        WHEN("locating the first extrusion and the move that precedes it") {
            size_t first_extrude = result.moves.size();
            for (size_t i = 0; i < result.moves.size(); ++i)
                if (result.moves[i].type == EMoveType::Extrude) { first_extrude = i; break; }

            THEN("an extrusion and a preceding move exist") {
                REQUIRE(first_extrude < result.moves.size());
                REQUIRE(first_extrude > 0);
            }

            THEN("the preceding prepare-stage move shares the extrusion's real Z (no leak to Y=0)") {
                const float extrude_z = result.moves[first_extrude].position.z();
                const float prev_z    = result.moves[first_extrude - 1].position.z();
                // The first extrusion is at the real Z=50; before the fix the
                // prepare-stage move's Z was pinned to the first-layer height
                // (0 here) instead, which back-transforms to model Y ~= 0 and
                // produces the phantom extrusion segment.
                REQUIRE_THAT(extrude_z, Catch::Matchers::WithinAbs(50.0, 1e-3));
                REQUIRE_THAT(prev_z,    Catch::Matchers::WithinAbs(50.0, 1e-3));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Regression tests for the two latent bugs the MachineKinematics refactor
// preserved deliberately and the follow-up commit fixed.
// ---------------------------------------------------------------------------

// Bug 1. _travel_to_z() emits full XYZ whenever the mapping must emit every
// axis, and it builds that point from m_pos. While the position is unknown,
// m_pos.xy is the uninitialised origin, which a reverse remap maps to the far
// corner of the bed. Belt kinematics guarded this; a Cartesian writer with an
// axis remap did not, and would command a rapid across the whole bed.
static void configure_lift_writer(GCodeWriter &writer)
{
    std::vector<unsigned int> extruder_ids { 0 };
    writer.set_extruders(extruder_ids);
    writer.set_extruder(0);
    writer.config.travel_speed.values       = { 100.0 };
    writer.config.travel_speed_z.values     = { 100.0 };
    writer.config.z_hop.values              = { 0.4 };
    writer.config.retract_lift_above.values = { 0.0 };
    writer.config.retract_lift_below.values = { 0.0 };
}

// Largest X word in a chunk of emitted G-code, or lowest() if none.
static double max_emitted_x(const std::string &gcode)
{
    double max_x = std::numeric_limits<double>::lowest();
    GCodeReader reader;
    reader.parse_buffer(gcode, [&max_x](GCodeReader &, const GCodeReader::GCodeLine &line) {
        if (line.cmd_is("G1") && line.has(X))
            max_x = std::max(max_x, double(line.x()));
    });
    return max_x;
}

static size_t count_g1(const std::string &gcode)
{
    size_t n = 0;
    GCodeReader reader;
    reader.parse_buffer(gcode, [&n](GCodeReader &, const GCodeReader::GCodeLine &line) {
        if (line.cmd_is("G1")) ++n;
    });
    return n;
}

SCENARIO("Axis remap: no lift is commanded through the uninitialised origin", "[GCodeWriter][remap]")
{
    // Reverse X: machine X = build_vol_max.x - logical X, so the uninitialised
    // origin maps to the far edge of the bed and is unmistakable in the output.
    const double bed_x = 250.0;

    GIVEN("a writer with a reverse-X remap and an unknown current position") {
        GCodeWriter writer;
        configure_lift_writer(writer);
        writer.set_axis_remap(6, 1, 2);
        writer.set_build_volume_max(Vec3d(bed_x, 250.0, 250.0));
        REQUIRE(writer.kinematics().must_emit_all_axes());
        REQUIRE_FALSE(writer.is_current_position_clear());

        WHEN("a z-hop is pending and we travel to the first point") {
            writer.lazy_lift(LiftType::NormalLift);
            const std::string gcode = writer.travel_to_xyz(Vec3d(10.0, 10.0, 5.0));

            THEN("nothing is commanded at the image of the origin") {
                // The destination maps to machine X = 250 - 10 = 240; the bogus
                // origin lift would have mapped to machine X = 250.
                REQUIRE(max_emitted_x(gcode) < bed_x - 1.0);
            }
            THEN("only the destination move is emitted") {
                REQUIRE(count_g1(gcode) == 1);
            }
        }
    }

    GIVEN("the same writer once its position is known") {
        GCodeWriter writer;
        configure_lift_writer(writer);
        writer.set_axis_remap(6, 1, 2);
        writer.set_build_volume_max(Vec3d(bed_x, 250.0, 250.0));
        writer.travel_to_xyz(Vec3d(20.0, 20.0, 5.0));
        REQUIRE(writer.is_current_position_clear());

        WHEN("a z-hop is pending and we travel again") {
            writer.lazy_lift(LiftType::NormalLift);
            const std::string gcode = writer.travel_to_xyz(Vec3d(30.0, 30.0, 5.0));

            THEN("the separate lift move is still emitted") {
                // Suppression must be pinned to the unknown position, not to the
                // presence of a remap.
                REQUIRE(count_g1(gcode) == 2);
            }
        }
    }

    GIVEN("an identity-mapping writer with an unknown position") {
        GCodeWriter writer;
        configure_lift_writer(writer);
        REQUIRE_FALSE(writer.kinematics().must_emit_all_axes());
        REQUIRE_FALSE(writer.is_current_position_clear());

        WHEN("a z-hop is pending and we travel to the first point") {
            writer.lazy_lift(LiftType::NormalLift);
            const std::string gcode = writer.travel_to_xyz(Vec3d(10.0, 10.0, 5.0));

            THEN("behaviour is unchanged: the lift is still emitted") {
                // Three moves, not two: with no remap and an unknown position the
                // destination is emitted as a separate XY move followed by its own
                // Z move, on top of the lift. That split is the pre-existing
                // identity-mapping path and must not change.
                REQUIRE(count_g1(gcode) == 3);
            }
        }
    }
}

SCENARIO("Axis remap: eager_lift does not lift, or record a lift, at an unknown position",
         "[GCodeWriter][remap]")
{
    GIVEN("a writer with a reverse-X remap and an unknown current position") {
        GCodeWriter writer;
        configure_lift_writer(writer);
        writer.set_axis_remap(6, 1, 2);
        writer.set_build_volume_max(Vec3d(250.0, 250.0, 250.0));
        REQUIRE_FALSE(writer.is_current_position_clear());

        WHEN("an eager lift is requested") {
            const std::string lift = writer.eager_lift(LiftType::NormalLift);

            THEN("no move is emitted") {
                REQUIRE(lift.empty());
            }
            THEN("no lift is recorded, so unlift does not descend from it") {
                // If m_lifted had been set while nothing was commanded, unlift()
                // would emit a descent from a height the machine never reached.
                REQUIRE(writer.unlift().empty());
            }
        }
    }

    GIVEN("an identity-mapping writer with an unknown position") {
        GCodeWriter writer;
        configure_lift_writer(writer);

        WHEN("an eager lift is requested") {
            const std::string lift = writer.eager_lift(LiftType::NormalLift);

            THEN("behaviour is unchanged: the lift is emitted and can be undone") {
                REQUIRE_FALSE(lift.empty());
                REQUIRE_FALSE(writer.unlift().empty());
            }
        }
    }
}

// Bug 2. extrude_arc_to_xy() emits G2/G3 with logical X/Y and I/J and never
// consulted the mapping. An arc is only representable when logical X and Y reach
// the machine unchanged -- which is a narrower question than "is the remap the
// identity", because a mapping that only touches Z leaves every emitted word alone.
SCENARIO("Arc support is decided by whether the mapping leaves X and Y alone", "[GCodeWriter][remap]")
{
    GIVEN("a Cartesian writer") {
        GCodeWriter writer;

        THEN("the identity mapping supports arcs") {
            REQUIRE(writer.kinematics().supports_arc_moves());
        }
        THEN("a Z-only negation still supports arcs") {
            // (+X, +Y, -Z): non-identity, but X, Y, I and J are all untouched.
            writer.set_axis_remap(0, 1, 5);
            REQUIRE(writer.kinematics().must_emit_all_axes());
            REQUIRE(writer.kinematics().supports_arc_moves());
        }
        THEN("a Z-only reversal still supports arcs") {
            writer.set_axis_remap(0, 1, 8);
            REQUIRE(writer.kinematics().supports_arc_moves());
        }
        THEN("swapping X and Y does not support arcs") {
            writer.set_axis_remap(1, 0, 2);
            REQUIRE_FALSE(writer.kinematics().supports_arc_moves());
        }
        THEN("the X-tilt style (x, z, y) remap does not support arcs") {
            writer.set_axis_remap(0, 2, 1);
            REQUIRE_FALSE(writer.kinematics().supports_arc_moves());
        }
    }

    GIVEN("a belt writer") {
        PrintConfig belt_config;
        belt_config.belt_printer.value              = true;
        belt_config.belt_slice_rotation.value       = BeltRotationAxis::X;
        belt_config.belt_slice_rotation_angle.value = 45.0;

        GCodeWriter writer;
        install_belt_kinematics(writer, belt_config);

        THEN("arcs are never supported, because the frame shears") {
            REQUIRE_FALSE(writer.kinematics().supports_arc_moves());
        }
    }
}

SCENARIO("An unrepresentable arc degrades to its chord rather than emitting a wrong G2/G3",
         "[GCodeWriter][remap]")
{
    auto emitted_commands = [](const std::string &gcode) {
        std::vector<std::string> cmds;
        GCodeReader reader;
        reader.parse_buffer(gcode, [&cmds](GCodeReader &, const GCodeReader::GCodeLine &line) {
            if (! line.cmd().empty()) cmds.emplace_back(line.cmd());
        });
        return cmds;
    };

    GIVEN("an identity-mapping writer") {
        GCodeWriter writer;
        configure_lift_writer(writer);

        WHEN("an arc is extruded") {
            const std::string gcode = writer.extrude_arc_to_xy(
                Vec2d(10.0, 0.0), Vec2d(5.0, 0.0), 0.0, /*is_ccw=*/true, "", /*force_no_extrusion=*/true);

            THEN("it is still a G3") {
                const auto cmds = emitted_commands(gcode);
                REQUIRE(cmds.size() == 1);
                REQUIRE(cmds.front() == "G3");
            }
        }
    }

    GIVEN("a writer whose mapping swaps X and Y") {
        GCodeWriter writer;
        configure_lift_writer(writer);
        writer.set_axis_remap(1, 0, 2);

        WHEN("an arc is extruded") {
            const std::string gcode = writer.extrude_arc_to_xy(
                Vec2d(10.0, 0.0), Vec2d(5.0, 0.0), 0.0, /*is_ccw=*/true, "", /*force_no_extrusion=*/true);

            THEN("no arc is emitted; it is approximated with linear moves") {
                const auto cmds = emitted_commands(gcode);
                REQUIRE(! cmds.empty());
                for (const auto &c : cmds)
                    REQUIRE(c == "G1");
            }
        }
    }

    // The first version of this test used dE = 0 with force_no_extrusion, which
    // hid a real bug: the capability check sat AFTER filament()->extrude(dE), so
    // the fallback into extrude_to_xy() advanced E twice. Extrusion accounting has
    // to be asserted with a positive dE.
    GIVEN("a writer whose mapping cannot express arcs, extruding a real amount") {
        GCodeWriter writer;
        configure_lift_writer(writer);
        writer.set_axis_remap(1, 0, 2);
        const double dE = 1.5;
        // used_filament() accumulates across moves; E() is reset per line in
        // relative-E mode, so it would only show the last segment.
        const double used_before = writer.filament()->used_filament();

        WHEN("an arc carrying that extrusion is emitted") {
            const std::string gcode = writer.extrude_arc_to_xy(
                Vec2d(10.0, 0.0), Vec2d(5.0, 0.0), dE, /*is_ccw=*/true, "", /*force_no_extrusion=*/false);

            THEN("exactly dE is accounted for, not twice dE") {
                REQUIRE_THAT(writer.filament()->used_filament() - used_before,
                             Catch::Matchers::WithinAbs(dE, 1e-6));
            }
            THEN("no G2/G3 survives") {
                REQUIRE(gcode.find("G2") == std::string::npos);
                REQUIRE(gcode.find("G3") == std::string::npos);
            }
        }
    }

    GIVEN("a writer whose mapping CAN express arcs, extruding a real amount") {
        GCodeWriter writer;
        configure_lift_writer(writer);
        const double dE = 1.5;
        // used_filament() accumulates across moves; E() is reset per line in
        // relative-E mode, so it would only show the last segment.
        const double used_before = writer.filament()->used_filament();

        WHEN("an arc carrying that extrusion is emitted") {
            const std::string gcode = writer.extrude_arc_to_xy(
                Vec2d(10.0, 0.0), Vec2d(5.0, 0.0), dE, /*is_ccw=*/true, "", /*force_no_extrusion=*/false);

            THEN("it is still a single arc and accounts for dE once") {
                REQUIRE(emitted_commands(gcode).size() == 1);
                REQUIRE_THAT(writer.filament()->used_filament() - used_before,
                             Catch::Matchers::WithinAbs(dE, 1e-6));
            }
        }
    }
}

TEST_CASE("Percent accelerations resolve against the option they are a percentage of", "[GCodeWriter]")
{
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({
        { "gcode_flavor",                       "marlin" },
        { "machine_max_acceleration_extruding", "20000,20000" },
        { "default_acceleration",               "4000" },
        { "initial_layer_acceleration",         "0" },
        { "outer_wall_acceleration",            "3000" },
        { "bridge_acceleration",                "50%" },
        { "sparse_infill_acceleration",         "25%" },
        { "internal_solid_infill_acceleration", "60%" },
        { "sparse_infill_density",              "20%" },
    });
    // get_abs_value_at() resolves each percentage through the ratio_over in the config definitions.
    const std::map<std::string, int> expected = {
        { "Bridge",                int(config.get_abs_value_at("bridge_acceleration", 0)) },
        { "Sparse infill",         int(config.get_abs_value_at("sparse_infill_acceleration", 0)) },
        { "Internal solid infill", int(config.get_abs_value_at("internal_solid_infill_acceleration", 0)) },
    };
    REQUIRE(expected.at("Bridge") == 1500);
    REQUIRE(expected.at("Sparse infill") == 1000);
    REQUIRE(expected.at("Internal solid infill") == 2400);

    std::map<std::string, std::set<int>> accelerations_by_role;
    std::string role;
    int         acceleration = 0;
    GCodeReader reader;
    reader.parse_buffer(Slic3r::Test::slice({ TestMesh::bridge }, config), [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
        float value;
        if (boost::starts_with(line.raw(), ";TYPE:"))
            role = line.raw().substr(6);
        else if (line.cmd_is("M204") && line.has_value('S', value))
            acceleration = int(value);
        else if (line.extruding(self) && line.dist_XY(self) > 0)
            accelerations_by_role[role].insert(acceleration);
    });
    for (const auto &[role_name, value] : expected) {
        INFO(role_name);
        REQUIRE(accelerations_by_role[role_name] == std::set<int>{ value });
    }
}

TEST_CASE("Initial filament start G-code tolerates a list shorter than the filament count", "[GCodeWriter]")
{
    // Two filaments but only one filament_start_gcode entry, and the only object prints on
    // filament 2, so the BBL export starts on a filament past the end of the list.
    DynamicPrintConfig config = multifilament_config(2);
    config.set_key_value("filament_start_gcode", new ConfigOptionStrings({"; filament start"}));

    Model model;
    auto *obj = model.add_object();
    obj->add_volume(cube(20));
    obj->add_instance();
    obj->config.set_key_value("extruder", new ConfigOptionInt(2));

    Print print;
    print.is_BBL_printer() = true;
    arrange_objects_on_test_bed(model, config);
    for (auto *mo : model.objects) {
        mo->ensure_on_bed();
        print.auto_assign_extruders(mo);
    }
    print.apply(model, config);
    print.validate();

    REQUIRE_NOTHROW(Slic3r::Test::gcode(print));
}
