#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "test_helpers.hpp"
#include "libslic3r/Config.hpp"
#include "libslic3r/CoolingAssistant.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/PrintBase.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/WarpMitigation.hpp"
#include "libslic3r/Format/bbs_3mf.hpp"
#include "libslic3r/Preset.hpp"
#include "libslic3r/Semver.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polyline.hpp"
#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/Utils.hpp"
#include "test_utils.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace Slic3r;
using Catch::Matchers::WithinAbs;

TEST_CASE("Automatic ASA cooling preserves a supplied project's extrusion geometry", "[WarpMitigation][.CoolingProject]")
{
    const char *path = std::getenv("ORCA_WARP_REGRESSION_3MF");
    if (path == nullptr) SKIP("Set ORCA_WARP_REGRESSION_3MF to a local ASA project");
    const auto previous_resources = resources_dir();
    set_resources_dir(std::string(TEST_DATA_DIR) + "/../../resources");
    const ScopeGuard restore_resources([&] { set_resources_dir(previous_resources); });
    ScopedSlic3rTemporaryDir temporary_dir("asa_cooling_project");
    DynamicPrintConfig imported;
    ConfigSubstitutionContext substitutions{ForwardCompatibilitySubstitutionRule::Enable};
    PlateDataPtrs plates;
    std::vector<Preset *> presets;
    bool is_bbl = false;
    Semver version;
    auto model = Model::read_from_file(path, &imported, &substitutions,
        LoadStrategy::LoadModel | LoadStrategy::LoadConfig | LoadStrategy::AddDefaultInstances,
        &plates, &presets, &is_bbl, &version);
    release_PlateData_list(plates);
    for (Preset *preset : presets) delete preset;
    auto config = DynamicPrintConfig::full_print_config();
    config.apply(imported);
    const auto risk = assess_warp_risk(config, model);
    REQUIRE_FALSE(risk.asa_cooling_slots.empty());
    // Hold the existing process mitigations constant while toggling Auto.
    apply_warp_mitigations(config, risk);
    config.set_deserialize_strict("automatic_warp_mitigation", "0");
    Print original;
    original.apply(model, config);
    original.process();
    const auto geometry = [](const Print &print) {
        std::vector<Points> paths;
        double volume = 0.;
        for (const auto *object : print.objects())
            for (const auto *layer : object->layers())
                for (const auto *region : layer->regions())
                    for (const auto *collection : {&region->perimeters, &region->fills}) {
                        volume += collection->total_volume();
                        for (const auto &line : collection->as_polylines())
                            paths.push_back(line.points);
                    }
        return std::make_pair(volume, paths);
    };
    const auto before = geometry(original);
    config.set_deserialize_strict("automatic_warp_mitigation", "1");
    original.apply(model, config);
    for (const auto *object : original.objects()) {
        CHECK(object->is_step_done(posPerimeters));
        CHECK(object->is_step_done(posInfill));
    }
    original.process();
    const auto after = geometry(original);
    REQUIRE_FALSE(before.second.empty());
    CHECK(before.second == after.second);
    CHECK_THAT(before.first, WithinAbs(after.first, 1e-9));
    for (size_t slot : risk.asa_cooling_slots)
        CHECK(original.config().reduce_fan_stop_start_freq.get_at(slot) == config.option<ConfigOptionBools>("reduce_fan_stop_start_freq")->get_at(slot));
    if (const char *output = std::getenv("ORCA_WARP_REGRESSION_GCODE")) {
        const std::string gcode = Test::gcode(original);
        REQUIRE_FALSE(gcode.empty());
        std::ofstream stream(output);
        REQUIRE(stream.is_open());
        stream << gcode;
        stream.close();
        REQUIRE(stream.good());
    }
}

TEST_CASE("Automatic ASA mitigation preserves profile fans and deposited paths", "[WarpMitigation][Cooling]")
{
    const bool short_layers = GENERATE(false, true);
    auto config = Test::multifilament_config(1, {
        {"filament_type", "ASA"}, {"brim_type", "outer_only"}, {"brim_width", "10"},
        {"brim_object_gap", "0.1"}, {"layer_height", "0.2"},
        {"initial_layer_print_height", "0.2"}, {"bottom_shell_thickness", "1"},
        {"curr_bed_type", "Textured PEI Plate"}, {"textured_plate_temp", "110"},
        {"chamber_temperature", "50"}, {"reduce_fan_stop_start_freq", "1"},
        {"fan_min_speed", "20"}, {"fan_max_speed", "25"},
        {"close_fan_the_first_x_layers", "0"}, {"full_fan_speed_layer", "0"},
        {"initial_layer_fan_speed", "-1"}, {"slow_down_for_layer_cooling", "0"},
        {"slow_down_layer_time", short_layers ? "100000" : "0"},
        {"fan_cooling_layer_time", short_layers ? "100001" : "0"},
        {"enable_overhang_speed", "1"}, {"slowdown_for_curled_perimeters", "1"},
        {"automatic_warp_mitigation", "0"}});
    const auto original = Test::slice({make_cube(80.,70.,1.)}, config);
    config.set_deserialize_strict("automatic_warp_mitigation", "1");
    const auto mitigated = Test::slice({make_cube(80.,70.,1.)}, config);
    const auto moves = [](const std::string &gcode) {
        std::vector<std::string> result;
        std::istringstream input(gcode);
        std::string line;
        while (std::getline(input, line))
            if (line.rfind("G0 ", 0) == 0 || line.rfind("G1 ", 0) == 0 ||
                line.rfind("G2 ", 0) == 0 || line.rfind("G3 ", 0) == 0)
            {
                // Auto mode may reduce first-layer speed, but must not change
                // coordinates, extrusion amounts or their order.
                std::istringstream words(line);
                std::string word, geometry;
                while (words >> word)
                    if (word.front() != 'F') geometry += word + ' ';
                result.push_back(geometry);
            }
        return result;
    };
    CHECK(moves(original) == moves(mitigated));
    if (short_layers) {
        CHECK(original.find("M106 S63") != std::string::npos);
        CHECK(mitigated.find("M106 S63") != std::string::npos);
    } else {
        CHECK(original.find("M106 S51") != std::string::npos);
        CHECK(mitigated.find("M106 S51") != std::string::npos);
        CHECK(mitigated.find("M106 S0") != std::string::npos);
    }
    CHECK(config.option<ConfigOptionBools>("reduce_fan_stop_start_freq")->get_at(0));
    REQUIRE(apply_material_cooling(config, config, config));
    const auto profile_cooled = Test::slice({make_cube(80.,70.,1.)}, config);
    CHECK(moves(mitigated) == moves(profile_cooled));
    CHECK(profile_cooled.find("M106 S51") == std::string::npos);
    CHECK(profile_cooled.find("M106 S63") == std::string::npos);
}

TEST_CASE("Automatic ASA mitigation retains bridge fan commands", "[WarpMitigation][Cooling]")
{
    auto config = Test::multifilament_config(1, {
        {"filament_type", "ASA"}, {"automatic_warp_mitigation", "1"},
        {"curr_bed_type", "Textured PEI Plate"}, {"textured_plate_temp", "110"},
        {"chamber_temperature", "50"}, {"layer_height", "0.2"},
        {"bottom_shell_thickness", "1"}, {"enable_overhang_bridge_fan", "1"},
        {"overhang_fan_speed", "100"}, {"internal_bridge_fan_speed", "100"},
        {"close_fan_the_first_x_layers", "0"}, {"full_fan_speed_layer", "0"},
        {"initial_layer_fan_speed", "-1"}, {"fan_min_speed", "20"}, {"fan_max_speed", "25"}});
    Print print;
    Model model;
    auto bridge = Test::mesh(Test::TestMesh::bridge);
    bridge.scale(Vec3f(4.f, 4.f, 1.f));
    REQUIRE(apply_material_cooling(config, config, config));
    Test::init_print({bridge}, print, model, config);
    CHECK_THAT(print.config().fan_max_speed.get_at(0), WithinAbs(0., 1e-9));
    const auto gcode = Test::gcode(print);
    CHECK(gcode.find("M106 S255") != std::string::npos);
}

TEST_CASE("Automatic first layer adhesion slows screened ASA without changing later layers", "[WarpMitigation][FirstLayer]")
{
    const std::string material = GENERATE(std::string("ASA"), std::string("TPU"));
    const bool automatic = GENERATE(false, true);
    const double span = GENERATE(2., 80.);
    auto config = Test::multifilament_config(1, {
        {"filament_type", material}, {"automatic_warp_mitigation", automatic ? "1" : "0"},
        {"brim_type", "no_brim"}, {"layer_height", "0.2"},
        {"initial_layer_print_height", "0.2"}, {"bottom_shell_thickness", "1"},
        {"curr_bed_type", "Textured PEI Plate"}, {"textured_plate_temp", "110"},
        {"chamber_temperature", "50"}, {"initial_layer_speed", "45"},
        {"initial_layer_infill_speed", "45"}, {"outer_wall_speed", "60"},
        {"inner_wall_speed", "60"}, {"small_perimeter_threshold", "0"},
        {"slow_down_layers", "0"}, {"slow_down_for_layer_cooling", "0"},
        {"filament_max_volumetric_speed", "100"}, {"gcode_comments", "1"}});
    const std::string output = Test::slice({make_cube(span,span,.6)}, config);
    GCodeReader reader;
    size_t first_moves = 0, later_moves = 0;
    reader.parse_buffer(output, [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
        if (!line.extruding(self) || line.dist_XY(self) <= 0. ||
            std::string(line.comment()).find("perimeter") == std::string::npos)
            return;
        const double speed = line.new_F(self) / 60.;
        if (line.new_Z(self) < .3) {
            ++first_moves;
            CHECK_THAT(speed, WithinAbs(automatic && material == "ASA" && span == 80. ? 33.75 : 45., .01));
        } else {
            ++later_moves;
            CHECK_THAT(speed, WithinAbs(60., .01));
        }
    });
    CHECK(first_moves > 0);
    CHECK(later_moves > 0);
}

TEST_CASE("Automatic warp mitigation changes effective adhesion without changing the project", "[WarpMitigation]")
{
    auto config = Test::multifilament_config(3, {
        {"filament_type", "ASA;ASA;TPU"}, {"brim_type", "painted"},
        {"brim_width", "0"}, {"brim_object_gap", "0.2"},
        {"automatic_warp_mitigation", "1"}});
    config.set_key_value("reduce_fan_stop_start_freq", new ConfigOptionBools({true, true, true}));
    Model model;
    Print print;
    Test::init_print({make_cube(80., 70., 1.)}, print, model, config);
    REQUIRE(assess_warp_risk(config, model).elevated());
    REQUIRE(print.objects().size() == 1);
    CHECK(print.config().reduce_fan_stop_start_freq.get_at(0));
    CHECK(print.config().reduce_fan_stop_start_freq.get_at(1));
    CHECK(print.config().reduce_fan_stop_start_freq.get_at(2));
    CHECK(config.option<ConfigOptionBools>("reduce_fan_stop_start_freq")->get_at(0));
    CHECK(print.objects().front()->config().brim_type.value == btOuterOnly);
    CHECK_THAT(print.objects().front()->config().brim_width.value,
        WithinAbs(assess_warp_risk(config, model).brim_mm(), 1e-6));
    CHECK(config.option<ConfigOptionEnum<BrimType>>("brim_type")->value == btPainted);
    CHECK_THAT(config.opt_float("brim_width"), WithinAbs(0., 1e-6));
    print.process();
    CHECK(print.has_brim());
    REQUIRE_FALSE(print.skirt_brim_groups().empty());
    CHECK_FALSE(print.skirt_brim_groups().front().brims.empty());
    std::vector<StringObjectException> mitigated_warnings;
    print.validate(&mitigated_warnings);
    CHECK(std::any_of(mitigated_warnings.begin(), mitigated_warnings.end(), [](const auto &warning) {
        return warning.opt_key == "automatic_warp_mitigation" && !warning.is_warning;
    }));

    config.set_key_value("automatic_warp_mitigation", new ConfigOptionBool(false));
    print.apply(model, config);
    CHECK(print.config().reduce_fan_stop_start_freq.get_at(0));
    CHECK(print.objects().front()->config().brim_type.value == btPainted);
    CHECK_THAT(print.objects().front()->config().brim_width.value, WithinAbs(0., 1e-6));
    std::vector<StringObjectException> warnings;
    print.validate(&warnings);
    CHECK(std::any_of(warnings.begin(), warnings.end(), [](const auto &warning) {
        return warning.opt_key == "automatic_warp_mitigation" && warning.is_warning;
    }));
}

TEST_CASE("Unused susceptible filament does not change PLA adhesion", "[WarpMitigation]")
{
    auto config = Test::multifilament_config(2, {
        {"filament_type", "PLA;ASA"}, {"brim_type", "no_brim"},
        {"curr_bed_type", "Textured PEI Plate"}, {"textured_plate_temp", "60;100"},
        {"chamber_temperature", "20;20"}, {"bottom_shell_thickness", "1"},
        {"brim_width", "0"}, {"automatic_warp_mitigation", "1"}});
    Model model;
    Print print;
    Test::init_print({Test::cube(10.)}, print, model, config);
    CHECK_FALSE(assess_warp_risk(config, model).elevated());
    REQUIRE(print.objects().size() == 1);
    CHECK(print.objects().front()->config().brim_type.value == btNoBrim);
    std::vector<StringObjectException> warnings;
    print.validate(&warnings);
    CHECK(std::none_of(warnings.begin(), warnings.end(), [](const auto &warning) {
        return warning.opt_key == "automatic_warp_mitigation";
    }));
}

TEST_CASE("Warp risk remains with a brim or raft and clears when the gradient disappears", "[WarpMitigation]")
{
    const int raft_layers = GENERATE(0, 3);
    auto config = Test::multifilament_config(1, {
        {"filament_type", "ASA"}, {"brim_type", "outer_only"},
        {"brim_width", "12"}, {"brim_object_gap", "0"},
        {"raft_layers", raft_layers}, {"automatic_warp_mitigation", "0"}});
    Model model;
    Print print;
    Test::init_print({make_cube(80., 70., 1.)}, print, model, config);
    std::vector<StringObjectException> warnings;
    print.validate(&warnings);
    CHECK(std::count_if(warnings.begin(), warnings.end(), [](const auto &warning) {
        return warning.opt_key == "automatic_warp_mitigation";
    }) == 1);

    config.set_deserialize_strict("chamber_temperature", "150");
    print.apply(model, config);
    warnings.clear();
    print.validate(&warnings);
    CHECK(std::none_of(warnings.begin(), warnings.end(), [](const auto &warning) {
        return warning.opt_key == "automatic_warp_mitigation";
    }));
}

TEST_CASE("New material families apply geometry based mitigation without editing saved settings", "[WarpMitigation]")
{
    const auto type = GENERATE("PLA", "PP", "PA", "PA6", "PA11", "PA12", "POM", "PPS", "PEEK", "PPS-CF", "PA-GF");
    auto config = Test::multifilament_config(1, {
        {"filament_type", type}, {"brim_type", "no_brim"}, {"brim_width", "0"},
        {"curr_bed_type", "Textured PEI Plate"}, {"textured_plate_temp", "100"},
        {"chamber_temperature", "20"}, {"nozzle_temperature", "400"},
        {"bottom_shell_thickness", "1"}, {"layer_height", "0.2"},
        {"automatic_warp_mitigation", "1"}});
    Model model;
    Print print;
    Test::init_print({make_cube(80., 70., 1.)}, print, model, config);
    const auto risk = assess_warp_risk(config, model);
    CHECK(risk.unknown_materials.empty());
    CHECK(risk.elevated());
    CHECK(risk.crystallization_contraction_mm > 0.);
    REQUIRE(print.objects().size() == 1);
    CHECK(print.objects().front()->config().brim_type.value == btOuterOnly);
    CHECK(config.option<ConfigOptionEnum<BrimType>>("brim_type")->value == btNoBrim);
    config.set_deserialize_strict("automatic_warp_mitigation", "0");
    print.apply(model, config);
    CHECK(print.objects().front()->config().brim_type.value == btNoBrim);
}
