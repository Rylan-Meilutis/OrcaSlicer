#include <catch2/catch_all.hpp>

#include <algorithm>

#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/PrintConfigConstants.hpp"
#include "libslic3r/LocalesUtils.hpp"
#include "libslic3r/Preset.hpp"
#include "libslic3r/SequentialGantryGeometry.hpp"
#include "libslic3r/Utils.hpp"

#include "test_utils.hpp"

#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>

#include <boost/filesystem.hpp>
#include <boost/nowide/fstream.hpp>
#include <nlohmann/json.hpp>

#include <sstream>

using namespace Slic3r;

TEST_CASE("New strength and overhang options preserve existing print defaults", "[Config]")
{
    const DynamicPrintConfig config = DynamicPrintConfig::full_print_config();

    CHECK_FALSE(config.opt_bool("seam_start_on_inner_wall"));
    CHECK_FALSE(config.opt_bool("arc_overhang_enabled"));
    CHECK_FALSE(config.opt_bool("bridge_overhang_before_walls"));
    CHECK(config.opt_bool("arc_overhang_bridges"));
    CHECK(config.opt_bool("arc_overhang_overhangs"));
    CHECK(config.opt_bool("arc_overhang_recursive_fill"));
    CHECK_THAT(config.opt_float("arc_overhang_overlap"), Catch::Matchers::WithinAbs(0.0, EPSILON));
    const auto *arc_overhang_flow = config.option<ConfigOptionPercent>("arc_overhang_flow_ratio");
    REQUIRE(arc_overhang_flow != nullptr);
    CHECK_THAT(arc_overhang_flow->value, Catch::Matchers::WithinAbs(100.0, EPSILON));
    DynamicPrintConfig legacy_arc_flow = config;
    REQUIRE_NOTHROW(legacy_arc_flow.set_deserialize_strict("arc_overhang_flow_ratio", "0.8"));
    const auto *legacy_arc_overhang_flow = legacy_arc_flow.option<ConfigOptionPercent>("arc_overhang_flow_ratio");
    REQUIRE(legacy_arc_overhang_flow != nullptr);
    CHECK_THAT(legacy_arc_overhang_flow->value, Catch::Matchers::WithinAbs(80.0, EPSILON));
    DynamicPrintConfig percent_arc_flow = config;
    REQUIRE_NOTHROW(percent_arc_flow.set_deserialize_strict("arc_overhang_flow_ratio", "80%"));
    const auto *percent_arc_overhang_flow = percent_arc_flow.option<ConfigOptionPercent>("arc_overhang_flow_ratio");
    REQUIRE(percent_arc_overhang_flow != nullptr);
    CHECK_THAT(percent_arc_overhang_flow->value, Catch::Matchers::WithinAbs(80.0, EPSILON));
    CHECK_THAT(config.opt_float("arc_overhang_speed"), Catch::Matchers::WithinAbs(5.0, EPSILON));
    CHECK_THAT(config.opt_float("arc_overhang_stabilization_speed"), Catch::Matchers::WithinAbs(5.0, EPSILON));
    const auto *arc_overhang_cooling = config.option<ConfigOptionPercents>("arc_overhang_cooling");
    REQUIRE(arc_overhang_cooling != nullptr);
    CHECK_THAT(arc_overhang_cooling->get_at(0), Catch::Matchers::WithinAbs(100.0, EPSILON));
    DynamicPrintConfig legacy_config = config;
    REQUIRE_NOTHROW(legacy_config.set_deserialize_strict("arc_overhang_cooling", "73%"));
    const auto *legacy_arc_overhang_cooling = legacy_config.option<ConfigOptionPercents>("arc_overhang_cooling");
    REQUIRE(legacy_arc_overhang_cooling != nullptr);
    CHECK_THAT(legacy_arc_overhang_cooling->get_at(0), Catch::Matchers::WithinAbs(73.0, EPSILON));
    CHECK(config.opt_int("arc_overhang_layers") == 1);
    CHECK(config.opt_int("arc_overhang_overhang_speed_layers") == 0);
    CHECK(config.opt_int("arc_overhang_bridge_speed_layers") == 0);
    CHECK_THAT(config.opt_float("arc_overhang_bridge_distance"), Catch::Matchers::WithinAbs(0.0, EPSILON));
    CHECK_THAT(config.opt_float("arc_overhang_min_overhang_distance"), Catch::Matchers::WithinAbs(0.0, EPSILON));
    const auto *inner_walls_flow = config.option<ConfigOptionPercent>("inner_walls_flow_ratio");
    REQUIRE(inner_walls_flow != nullptr);
    CHECK_THAT(inner_walls_flow->value, Catch::Matchers::WithinAbs(100.0, EPSILON));
    DynamicPrintConfig legacy_third_wall_config = config;
    REQUIRE_NOTHROW(legacy_third_wall_config.set_deserialize_strict("third_wall_flow_ratio", "1.25"));
    const auto *legacy_inner_walls_flow =
        legacy_third_wall_config.option<ConfigOptionPercent>("inner_walls_flow_ratio");
    REQUIRE(legacy_inner_walls_flow != nullptr);
    CHECK_THAT(legacy_inner_walls_flow->value, Catch::Matchers::WithinAbs(125.0, EPSILON));
    DynamicPrintConfig percent_inner_walls_config = config;
    REQUIRE_NOTHROW(percent_inner_walls_config.set_deserialize_strict("inner_walls_flow_ratio", "125%"));
    const auto *percent_inner_walls_flow =
        percent_inner_walls_config.option<ConfigOptionPercent>("inner_walls_flow_ratio");
    REQUIRE(percent_inner_walls_flow != nullptr);
    CHECK_THAT(percent_inner_walls_flow->value, Catch::Matchers::WithinAbs(125.0, EPSILON));
    DynamicPrintConfig legacy_stagger_flow_config = config;
    REQUIRE_NOTHROW(legacy_stagger_flow_config.set_deserialize_strict(
        "staggered_perimeter_flow_ratio", "105%"));
    const auto *legacy_stagger_flow =
        legacy_stagger_flow_config.option<ConfigOptionPercent>("inner_walls_flow_ratio");
    REQUIRE(legacy_stagger_flow != nullptr);
    CHECK_THAT(legacy_stagger_flow->value, Catch::Matchers::WithinAbs(105.0, EPSILON));
    CHECK_FALSE(config.opt_bool("staggered_perimeters"));
    CHECK(config.opt_enum<PerimeterLayeringMode>("perimeter_layering") ==
          PerimeterLayeringMode::Standard);
    CHECK(config.opt_enum<TopSurfaceZMode>("top_surface_z_mode") ==
          TopSurfaceZMode::Disabled);
    CHECK_THAT(config.opt_float("smooth_outer_wall_layer_height"),
               Catch::Matchers::WithinAbs(0.1, EPSILON));
    CHECK_THAT(config.opt_float("interlocking_wall_amplitude"),
               Catch::Matchers::WithinAbs(0.08, EPSILON));
    CHECK_THAT(config.opt_float("interlocking_wall_wavelength"),
               Catch::Matchers::WithinAbs(6., EPSILON));
    CHECK_THAT(config.opt_float("interlocking_wall_resolution"),
               Catch::Matchers::WithinAbs(0.5, EPSILON));
    CHECK(config.opt_bool("staggered_perimeters_inner_only"));
    const auto *stagger_offset = config.option<ConfigOptionPercent>("staggered_perimeter_offset");
    REQUIRE(stagger_offset != nullptr);
    CHECK_THAT(stagger_offset->value, Catch::Matchers::WithinAbs(50.0, EPSILON));
    CHECK_FALSE(config.opt_bool("nonplanar_infill"));
    CHECK_THAT(config.opt_float("nonplanar_infill_amplitude"),
               Catch::Matchers::WithinAbs(0.08, EPSILON));
    CHECK_THAT(config.opt_float("nonplanar_infill_wavelength"),
               Catch::Matchers::WithinAbs(6., EPSILON));
    CHECK_THAT(config.opt_float("nonplanar_infill_resolution"),
               Catch::Matchers::WithinAbs(1., EPSILON));
    CHECK_FALSE(config.opt_bool("fuzzy_skin_top_surface"));
    CHECK(config.opt_int("support_interface_top_temperature") == 0);
    CHECK_FALSE(config.opt_bool("slow_down_layer_above_dissimilar_support_interface"));
    CHECK(config.opt_int("hull_line_extra_perimeters") == 0);
    CHECK_THAT(config.opt_float("hull_line_perimeter_expansion"),
               Catch::Matchers::WithinAbs(0., EPSILON));
    CHECK(config.opt_enum<LocalizedShrinkageStrategy>("localized_shrinkage_strategy") ==
          LocalizedShrinkageStrategy::Custom);
    CHECK_THAT(config.opt_float("localized_shrinkage_infill_wall_gap"),
               Catch::Matchers::WithinAbs(0., EPSILON));
    CHECK_THAT(config.opt_float("localized_shrinkage_section_width"),
               Catch::Matchers::WithinAbs(0., EPSILON));
    CHECK_THAT(config.opt_float("localized_shrinkage_section_spacing"),
               Catch::Matchers::WithinAbs(10., EPSILON));
    CHECK_THAT(config.opt_float("localized_shrinkage_perforation_diameter"),
               Catch::Matchers::WithinAbs(0., EPSILON));
    CHECK_THAT(config.opt_float("localized_shrinkage_perforation_spacing"),
               Catch::Matchers::WithinAbs(3., EPSILON));
    const auto *interface_speed =
        config.option<ConfigOptionFloatOrPercent>("dissimilar_support_interface_speed");
    REQUIRE(interface_speed != nullptr);
    CHECK(interface_speed->percent);
    CHECK_THAT(interface_speed->value, Catch::Matchers::WithinAbs(50., EPSILON));
    CHECK(config.opt_int("dissimilar_support_interface_speed_layers") == 0);
    const ConfigOptionDef *spool_sync = print_config_def.get("sync_spool_manager_filament_names");
    REQUIRE(spool_sync != nullptr);
    const auto *spool_sync_default = dynamic_cast<const ConfigOptionBool *>(spool_sync->default_value.get());
    REQUIRE(spool_sync_default != nullptr);
    CHECK_FALSE(spool_sync_default->value);
    const ConfigOptionDef *embed_spool_names =
        print_config_def.get("embed_spool_manager_filament_names");
    REQUIRE(embed_spool_names != nullptr);
    const auto *embed_spool_names_default =
        dynamic_cast<const ConfigOptionBool *>(embed_spool_names->default_value.get());
    REQUIRE(embed_spool_names_default != nullptr);
    CHECK(embed_spool_names_default->value);
    const ConfigOptionDef *spool_sync_mode = print_config_def.get("spool_manager_sync_mode");
    REQUIRE(spool_sync_mode != nullptr);
    const auto *spool_sync_mode_default =
        dynamic_cast<const ConfigOptionEnum<SpoolManagerSyncMode> *>(spool_sync_mode->default_value.get());
    REQUIRE(spool_sync_mode_default != nullptr);
    CHECK(spool_sync_mode_default->value == smsmColorsAndProfiles);
    const ConfigOptionDef *filament_plugin_endpoint =
        print_config_def.get("octoprint_filament_plugin_endpoint");
    REQUIRE(filament_plugin_endpoint != nullptr);
    const auto *filament_plugin_endpoint_default =
        dynamic_cast<const ConfigOptionString *>(filament_plugin_endpoint->default_value.get());
    REQUIRE(filament_plugin_endpoint_default != nullptr);
    CHECK(filament_plugin_endpoint_default->value.empty());
    REQUIRE(print_config_def.get("octoprint_spool_profile_mappings") != nullptr);
    REQUIRE(print_config_def.get("octoprint_material_profile_mappings") != nullptr);
    REQUIRE(print_config_def.get("octoprint_default_filament_profile") != nullptr);
    CHECK(config.opt_int("max_filament_colors") == 0);
    const std::vector<std::string> &printer_options = Preset::printer_options();
    CHECK(std::find(printer_options.begin(), printer_options.end(), "max_filament_colors") !=
          printer_options.end());
    const std::vector<std::string> &physical_printer_options = PhysicalPrinter::printer_options();
    CHECK(std::find(physical_printer_options.begin(), physical_printer_options.end(), "printhost_apikey") !=
          physical_printer_options.end());
    CHECK(std::find(physical_printer_options.begin(), physical_printer_options.end(),
                    "sync_spool_manager_filament_names") != physical_printer_options.end());
    CHECK(std::find(physical_printer_options.begin(), physical_printer_options.end(),
                    "embed_spool_manager_filament_names") != physical_printer_options.end());
    CHECK(std::find(physical_printer_options.begin(), physical_printer_options.end(),
                    "spool_manager_sync_mode") != physical_printer_options.end());
    CHECK(std::find(physical_printer_options.begin(), physical_printer_options.end(),
                    "octoprint_filament_plugin_endpoint") != physical_printer_options.end());
    CHECK(std::find(physical_printer_options.begin(), physical_printer_options.end(),
                    "octoprint_spool_profile_mappings") != physical_printer_options.end());
    CHECK(std::find(physical_printer_options.begin(), physical_printer_options.end(),
                    "octoprint_material_profile_mappings") != physical_printer_options.end());
    CHECK(std::find(physical_printer_options.begin(), physical_printer_options.end(),
                    "octoprint_default_filament_profile") != physical_printer_options.end());
    CHECK(config.opt_string("sequential_print_gantry_geometry").empty());
    CHECK(config.opt_string("sequential_print_gantry_model").empty());
    CHECK_THAT(config.opt_float("nonplanar_toolhead_clearance_angle"),
               Catch::Matchers::WithinAbs(0., EPSILON));
    CHECK_FALSE(config.opt_bool("nonplanar_top_surface"));
    CHECK_THAT(config.opt_float("nonplanar_top_surface_max_angle"),
               Catch::Matchers::WithinAbs(45., EPSILON));
    CHECK(config.opt_int("nonplanar_top_surface_layers") == 5);
    CHECK_THAT(config.opt_float("nonplanar_top_surface_resolution"),
               Catch::Matchers::WithinAbs(0.2, EPSILON));
    CHECK_THAT(config.opt_float("nonplanar_top_surface_min_height"),
               Catch::Matchers::WithinAbs(0.05, EPSILON));
    CHECK_FALSE(config.opt_bool("support_ironing_nonplanar"));
    CHECK_THAT(config.opt_float("support_ironing_nonplanar_max_angle"),
               Catch::Matchers::WithinAbs(45., EPSILON));
    CHECK_THAT(config.opt_float("support_ironing_nonplanar_resolution"),
               Catch::Matchers::WithinAbs(0.2, EPSILON));
    const std::vector<std::string> &print_options = Preset::print_options();
    for (const char *key : {"perimeter_layering", "top_surface_z_mode", "smooth_outer_wall_layer_height",
                            "interlocking_wall_amplitude", "interlocking_wall_wavelength",
                            "interlocking_wall_resolution", "nonplanar_top_surface", "nonplanar_top_surface_max_angle",
                            "nonplanar_top_surface_layers", "nonplanar_top_surface_resolution",
                            "nonplanar_top_surface_min_height",
                            "nonplanar_infill", "nonplanar_infill_amplitude",
                            "nonplanar_infill_wavelength", "nonplanar_infill_resolution",
                            "support_ironing_nonplanar",
                            "support_ironing_nonplanar_max_angle",
                            "support_ironing_nonplanar_resolution",
                            "localized_shrinkage_strategy",
                            "localized_shrinkage_infill_wall_gap",
                            "localized_shrinkage_section_width",
                            "localized_shrinkage_section_spacing",
                            "localized_shrinkage_perforation_diameter",
                            "localized_shrinkage_perforation_spacing",
                            "slow_down_layer_above_dissimilar_support_interface",
                            "dissimilar_support_interface_speed",
                            "dissimilar_support_interface_speed_layers"}) {
        CHECK(std::find(print_options.begin(), print_options.end(), key) != print_options.end());
    }
    CHECK(std::find(printer_options.begin(), printer_options.end(), "nonplanar_toolhead_clearance_angle") !=
          printer_options.end());
    const std::vector<std::string> &filament_options = Preset::filament_options();
    CHECK(std::find(filament_options.begin(), filament_options.end(), "hull_line_mitigation") !=
          filament_options.end());
    CHECK(std::find(filament_options.begin(), filament_options.end(), "hull_line_max_layer_time_variation") !=
          filament_options.end());
    const auto *hull_line_mitigation =
        config.option<ConfigOptionBools>("hull_line_mitigation");
    const auto *hull_line_variation =
        config.option<ConfigOptionPercents>("hull_line_max_layer_time_variation");
    REQUIRE(hull_line_mitigation != nullptr);
    REQUIRE(hull_line_variation != nullptr);
    REQUIRE_FALSE(hull_line_mitigation->values.empty());
    REQUIRE_FALSE(hull_line_variation->values.empty());
    CHECK_FALSE(hull_line_mitigation->values.front());
    CHECK_THAT(hull_line_variation->values.front(),
               Catch::Matchers::WithinAbs(25., EPSILON));
}

TEST_CASE("Legacy perimeter layering options preserve compatible brick and nonplanar features", "[Config][PerimeterLayering]")
{
    DynamicPrintConfig legacy_brick;
    legacy_brick.set_key_value("staggered_perimeters", new ConfigOptionBool(true));
    legacy_brick.handle_legacy_composite();
    CHECK(legacy_brick.opt_enum<PerimeterLayeringMode>("perimeter_layering") ==
          PerimeterLayeringMode::Brick);

    DynamicPrintConfig full_legacy_brick = DynamicPrintConfig::full_print_config();
    full_legacy_brick.set_key_value("staggered_perimeters", new ConfigOptionBool(true));
    full_legacy_brick.handle_legacy_composite();
    CHECK(full_legacy_brick.opt_enum<PerimeterLayeringMode>("perimeter_layering") ==
          PerimeterLayeringMode::Brick);

    DynamicPrintConfig legacy_nonplanar;
    legacy_nonplanar.set_key_value("nonplanar_top_surface", new ConfigOptionBool(true));
    legacy_nonplanar.handle_legacy_composite();
    CHECK(legacy_nonplanar.opt_enum<PerimeterLayeringMode>("perimeter_layering") ==
          PerimeterLayeringMode::Standard);
    CHECK(legacy_nonplanar.opt_enum<TopSurfaceZMode>("top_surface_z_mode") ==
          TopSurfaceZMode::NonplanarTopSurface);

    DynamicPrintConfig legacy_nonplanar_selector;
    legacy_nonplanar_selector.set_key_value("perimeter_layering",
        new ConfigOptionEnum<PerimeterLayeringMode>(PerimeterLayeringMode::Nonplanar));
    legacy_nonplanar_selector.handle_legacy_composite();
    CHECK(legacy_nonplanar_selector.opt_enum<PerimeterLayeringMode>("perimeter_layering") ==
          PerimeterLayeringMode::Standard);
    CHECK(legacy_nonplanar_selector.opt_enum<TopSurfaceZMode>("top_surface_z_mode") ==
          TopSurfaceZMode::NonplanarTopSurface);
    CHECK(legacy_nonplanar_selector.opt_bool("nonplanar_top_surface"));

    DynamicPrintConfig selected_nonplanar = DynamicPrintConfig::full_print_config();
    selected_nonplanar.set_key_value("top_surface_z_mode",
        new ConfigOptionEnum<TopSurfaceZMode>(TopSurfaceZMode::NonplanarTopSurface));
    selected_nonplanar.handle_legacy_composite();
    CHECK(selected_nonplanar.opt_enum<TopSurfaceZMode>("top_surface_z_mode") ==
          TopSurfaceZMode::NonplanarTopSurface);
    CHECK(selected_nonplanar.opt_bool("nonplanar_top_surface"));
    CHECK_FALSE(selected_nonplanar.opt_bool("zaa_enabled"));

    DynamicPrintConfig selected_contouring = DynamicPrintConfig::full_print_config();
    selected_contouring.set_key_value("top_surface_z_mode",
        new ConfigOptionEnum<TopSurfaceZMode>(TopSurfaceZMode::ZContouring));
    selected_contouring.handle_legacy_composite();
    CHECK(selected_contouring.opt_enum<TopSurfaceZMode>("top_surface_z_mode") ==
          TopSurfaceZMode::ZContouring);
    CHECK(selected_contouring.opt_bool("zaa_enabled"));
    CHECK_FALSE(selected_contouring.opt_bool("nonplanar_top_surface"));

    DynamicPrintConfig selected_hybrid = DynamicPrintConfig::full_print_config();
    selected_hybrid.set_key_value("top_surface_z_mode",
        new ConfigOptionEnum<TopSurfaceZMode>(
            TopSurfaceZMode::NonplanarWithZContouringFallback));
    selected_hybrid.handle_legacy_composite();
    CHECK(selected_hybrid.opt_enum<TopSurfaceZMode>("top_surface_z_mode") ==
          TopSurfaceZMode::NonplanarWithZContouringFallback);
    CHECK(selected_hybrid.opt_bool("zaa_enabled"));
    CHECK(selected_hybrid.opt_bool("nonplanar_top_surface"));

    DynamicPrintConfig serialized_hybrid = DynamicPrintConfig::full_print_config();
    REQUIRE_NOTHROW(serialized_hybrid.set_deserialize_strict(
        "top_surface_z_mode", "nonplanar_with_z_contouring_fallback"));
    CHECK(serialized_hybrid.opt_enum<TopSurfaceZMode>("top_surface_z_mode") ==
          TopSurfaceZMode::NonplanarWithZContouringFallback);

    DynamicPrintConfig conflicting;
    conflicting.set_key_value("staggered_perimeters", new ConfigOptionBool(true));
    conflicting.set_key_value("nonplanar_top_surface", new ConfigOptionBool(true));
    conflicting.handle_legacy_composite();
    CHECK(conflicting.opt_enum<PerimeterLayeringMode>("perimeter_layering") ==
          PerimeterLayeringMode::Brick);

    PrintRegionConfig compatible;
    compatible.perimeter_layering.value = PerimeterLayeringMode::Brick;
    compatible.staggered_perimeters.value = true;
    compatible.nonplanar_top_surface.value = true;
    CHECK(brick_perimeters_enabled(compatible));
    CHECK(nonplanar_perimeters_enabled(compatible));

    compatible.perimeter_layering.value = PerimeterLayeringMode::InterlockingWalls;
    CHECK_FALSE(brick_perimeters_enabled(compatible));
    CHECK(nonplanar_perimeters_enabled(compatible));
}

TEST_CASE("Bundled Prusa sequential gantry models match their printer notes", "[Config][Arrange]")
{
    struct ResourcesDirGuard {
        std::string original = resources_dir();
        ~ResourcesDirGuard() { set_resources_dir(original); }
    } resources_dir_guard;
    const boost::filesystem::path source_root =
        boost::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
    set_resources_dir((source_root / "resources").string());

    const std::array<std::pair<const char*, const char*>, 8> cases{{
        {"PRINTER_MODEL_COREONE", "prusa3d_coreone_gantry.stl"},
        {"PRINTER_MODEL_COREONE_L", "prusa3d_coreoneL_gantry.stl"},
        {"PRINTER_MODEL_COREONE_INDX\nSEQ_ARRANGE_MODEL_COREONE_INDX", "prusa3d_coreone_indx_gantry.stl"},
        {"PRINTER_MODEL_MK4IS", "prusa3d_mk4_gantry.stl"},
        {"PRINTER_MODEL_MK4S", "prusa3d_mk4s_gantry.stl"},
        {"PRINTER_MODEL_MK3", "prusa3d_mk3s_gantry.stl"},
        {"PRINTER_MODEL_XLIS", "prusa3d_xl_gantry.stl"},
        {"PRINTER_MODEL_MINI", "prusa3d_mini_gantry.stl"}
    }};

    for (const auto& [notes, filename] : cases) {
        DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
        config.set_key_value("printer_notes", new ConfigOptionString(notes));
        const SequentialGantryGeometry geometry = load_sequential_gantry_geometry(config);
        CAPTURE(notes);
        REQUIRE_FALSE(geometry.empty());
        CHECK(boost::filesystem::path(geometry.model_path).filename().string() == filename);
        CHECK(boost::filesystem::exists(geometry.model_path));
    }
}

TEST_CASE("Custom sequential gantry geometry is machine-independent", "[Config][Arrange]")
{
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_key_value("sequential_print_gantry_geometry", new ConfigOptionString(R"({
        "gantry_model_filename": "custom_gantry.stl",
        "slices": [
            {"height": "2", "type": "convex", "polygons": ["-30,-40; 30,-40; 30,40; -30,40"]},
            {"height": "18", "type": "box", "polygons": ["-200,-10; 200,-10; 200,10; -200,10"]}
        ]
    })"));

    const SequentialGantryGeometry geometry = load_sequential_gantry_geometry(config);
    REQUIRE_FALSE(geometry.empty());
    CHECK_THAT(geometry.conservative_clearance_radius(), Catch::Matchers::WithinAbs(50., EPSILON));
    CHECK_THAT(geometry.first_box_height(), Catch::Matchers::WithinAbs(18., EPSILON));
    CHECK(geometry.model_path == "custom_gantry.stl");
}

SCENARIO("Generic config validation performs as expected.", "[Config]") {
    GIVEN("A config generated from default options") {
        Slic3r::DynamicPrintConfig config = Slic3r::DynamicPrintConfig::full_print_config();
        WHEN( "outer_wall_line_width is set to 250%, a valid value") {
            config.set_deserialize_strict("outer_wall_line_width", "250%");
            THEN( "The config is read as valid.") {
                REQUIRE(config.validate().empty());
            }
        }
        WHEN( "outer_wall_line_width is set to -10, an invalid value") {
            config.set("outer_wall_line_width", -10);
            THEN( "Validate returns error") {
                REQUIRE_FALSE(config.validate().empty());
            }
        }

        WHEN( "wall_loops is set to -10, an invalid value") {
            config.set("wall_loops", -10);
            THEN( "Validate returns error") {
                REQUIRE_FALSE(config.validate().empty());
            }
        }
    }
}

SCENARIO("Config accessor functions perform as expected.", "[Config]") {
    GIVEN("A config generated from default options") {
        Slic3r::DynamicPrintConfig config = Slic3r::DynamicPrintConfig::full_print_config();
        WHEN("A boolean option is set to a boolean value") {
            REQUIRE_NOTHROW(config.set("gcode_comments", true));
            THEN("The underlying value is set correctly.") {
                REQUIRE(config.opt<ConfigOptionBool>("gcode_comments")->getBool() == true);
            }
        }
        WHEN("A boolean option is set to a string value representing a 0 or 1") {
            CHECK_NOTHROW(config.set_deserialize_strict("gcode_comments", "1"));
            THEN("The underlying value is set correctly.") {
                REQUIRE(config.opt<ConfigOptionBool>("gcode_comments")->getBool() == true);
            }
        }
        WHEN("A boolean option is set to a string value representing something other than 0 or 1") {
            THEN("A BadOptionTypeException exception is thrown.") {
                REQUIRE_THROWS_AS(config.set("gcode_comments", "Z"), BadOptionTypeException);
            }
            AND_THEN("Value is unchanged.") {
                REQUIRE(config.opt<ConfigOptionBool>("gcode_comments")->getBool() == false);
            }
        }
        WHEN("A boolean option is set to an int value") {
            THEN("A BadOptionTypeException exception is thrown.") {
                REQUIRE_THROWS_AS(config.set("gcode_comments", 1), BadOptionTypeException);
            }
        }
        WHEN("A numeric option is set from serialized string") {
            config.set_deserialize_strict("raft_layers", "20");
            THEN("The underlying value is set correctly.") {
                REQUIRE(config.opt<ConfigOptionInt>("raft_layers")->getInt() == 20);
            }
        }
	WHEN("An integer-based option is set through the integer interface") {
	    config.set("raft_layers", 100);
	    THEN("The underlying value is set correctly.") {
		REQUIRE(config.opt<ConfigOptionInt>("raft_layers")->getInt() == 100);
	    }
        }
        WHEN("An floating-point option is set through the integer interface") {
            config.set("max_bridge_length", 10);
            THEN("The underlying value is set correctly.") {
                REQUIRE(config.opt<ConfigOptionFloat>("max_bridge_length")->getFloat() == 10.0);
            }
        }
        WHEN("A floating-point option is set through the double interface") {
            config.set("max_bridge_length", 5.5);
            THEN("The underlying value is set correctly.") {
                REQUIRE(config.opt<ConfigOptionFloat>("max_bridge_length")->getFloat() == 5.5);
            }
        }
        WHEN("An integer-based option is set through the double interface") {
            THEN("A BadOptionTypeException exception is thrown.") {
                REQUIRE_THROWS_AS(config.set("top_shell_layers", 5.5), BadOptionTypeException);
            }
        }
        WHEN("A numeric option is set to a non-numeric value.") {
	    auto prev_value = config.opt<ConfigOptionFloat>("max_bridge_length")->getFloat();
            THEN("A BadOptionTypeException exception is thrown.") {
                REQUIRE_THROWS_AS(config.set_deserialize_strict("max_bridge_length", "zzzz"), BadOptionValueException);
            }
            THEN("The value does not change.") {
                REQUIRE(config.opt<ConfigOptionFloat>("max_bridge_length")->getFloat() == prev_value);
            }
        }
        WHEN("A string option is set through the string interface") {
            config.set("machine_end_gcode", "100");
            THEN("The underlying value is set correctly.") {
                REQUIRE(config.opt<ConfigOptionString>("machine_end_gcode")->value == "100");
            }
        }
        WHEN("A string option is set through the integer interface") {
            config.set("machine_end_gcode", 100);
            THEN("The underlying value is set correctly.") {
                REQUIRE(config.opt<ConfigOptionString>("machine_end_gcode")->value == "100");
            }
        }
        WHEN("A string option is set through the double interface") {
            config.set("machine_end_gcode", 100.5);
            THEN("The underlying value is set correctly.") {
                REQUIRE(config.opt<ConfigOptionString>("machine_end_gcode")->value == float_to_string_decimal_point(100.5));
            }
        }
        WHEN("A float or percent is set as a percent through the string interface.") {
            config.set_deserialize_strict("initial_layer_line_width", "100%");
            THEN("Value and percent flag are 100/true") {
                auto tmp = config.opt<ConfigOptionFloatOrPercent>("initial_layer_line_width");
                REQUIRE(tmp->percent == true);
                REQUIRE(tmp->value == 100);
            }
        }
        WHEN("A float or percent is set as a float through the string interface.") {
            config.set_deserialize_strict("initial_layer_line_width", "100");
            THEN("Value and percent flag are 100/false") {
                auto tmp = config.opt<ConfigOptionFloatOrPercent>("initial_layer_line_width");
                REQUIRE(tmp->percent == false);
                REQUIRE(tmp->value == 100);
            }
        }
        WHEN("A float or percent is set as a float through the int interface.") {
            config.set("initial_layer_line_width", 100);
            THEN("Value and percent flag are 100/false") {
                auto tmp = config.opt<ConfigOptionFloatOrPercent>("initial_layer_line_width");
                REQUIRE(tmp->percent == false);
                REQUIRE(tmp->value == 100);
            }
        }
        WHEN("A float or percent is set as a float through the double interface.") {
            config.set("initial_layer_line_width", 100.5);
            THEN("Value and percent flag are 100.5/false") {
                auto tmp = config.opt<ConfigOptionFloatOrPercent>("initial_layer_line_width");
                REQUIRE(tmp->percent == false);
                REQUIRE(tmp->value == 100.5);
            }
        }
        WHEN("A numeric vector is set from serialized string") {
	    config.set_deserialize_strict("temperature_vitrification", "10,20");
            THEN("The underlying value is set correctly.") {
                CHECK(config.opt<ConfigOptionInts>("temperature_vitrification")->get_at(0) == 10);
                CHECK(config.opt<ConfigOptionInts>("temperature_vitrification")->get_at(1) == 20);
            }
        }
	// FIXME: Design better accessors for vector elements
	// The following isn't supported and probably shouldn't be:
	// WHEN("An integer-based vector option is set through the integer interface") {
	//     config.set("temperature_vitrification", 100);
	//     THEN("The underlying value is set correctly.") {
	// 	REQUIRE(config.opt<ConfigOptionInts>("temperature_vitrification")->get_at(0) == 100);
	//     }
        // }
	WHEN("An integer-based vector option is set through the set_key_value interface") {
	    config.set_key_value("temperature_vitrification", new ConfigOptionInts{10,20});
	    THEN("The underlying value is set correctly.") {
                CHECK(config.opt<ConfigOptionInts>("temperature_vitrification")->get_at(0) == 10);
                CHECK(config.opt<ConfigOptionInts>("temperature_vitrification")->get_at(1) == 20);
	    }
        }
        WHEN("An invalid option is requested during set.") {
            THEN("A BadOptionTypeException exception is thrown.") {
                REQUIRE_THROWS_AS(config.set("deadbeef_invalid_option", 1), UnknownOptionException);
                REQUIRE_THROWS_AS(config.set("deadbeef_invalid_option", 1.0), UnknownOptionException);
                REQUIRE_THROWS_AS(config.set("deadbeef_invalid_option", "1"), UnknownOptionException);
                REQUIRE_THROWS_AS(config.set("deadbeef_invalid_option", true), UnknownOptionException);
            }
        }

        WHEN("An invalid option is requested during get.") {
            THEN("A UnknownOptionException exception is thrown.") {
                REQUIRE_THROWS_AS(config.option_throw<ConfigOptionString>("deadbeef_invalid_option", false), UnknownOptionException);
                REQUIRE_THROWS_AS(config.option_throw<ConfigOptionFloat>("deadbeef_invalid_option", false), UnknownOptionException);
                REQUIRE_THROWS_AS(config.option_throw<ConfigOptionInt>("deadbeef_invalid_option", false), UnknownOptionException);
                REQUIRE_THROWS_AS(config.option_throw<ConfigOptionBool>("deadbeef_invalid_option", false), UnknownOptionException);
            }
        }
        WHEN("An invalid option is requested during opt.") {
            THEN("A UnknownOptionException exception is thrown.") {
                REQUIRE_THROWS_AS(config.option_throw<ConfigOptionString>("deadbeef_invalid_option", false), UnknownOptionException);
                REQUIRE_THROWS_AS(config.option_throw<ConfigOptionFloat>("deadbeef_invalid_option", false), UnknownOptionException);
                REQUIRE_THROWS_AS(config.option_throw<ConfigOptionInt>("deadbeef_invalid_option", false), UnknownOptionException);
                REQUIRE_THROWS_AS(config.option_throw<ConfigOptionBool>("deadbeef_invalid_option", false), UnknownOptionException);
            }
        }

        WHEN("getX called on an unset option.") {
            THEN("The default is returned.") {
                REQUIRE(config.opt_float("layer_height") == INITIAL_LAYER_HEIGHT);
                REQUIRE(config.opt_int("raft_layers") == INITIAL_RAFT_LAYERS);
                REQUIRE(config.opt_bool("reduce_crossing_wall") == INITIAL_REDUCE_CROSSING_WALL);
            }
        }

        WHEN("opt_float called on an option that has been set.") {
            config.set("layer_height", INITIAL_LAYER_HEIGHT*2);
            THEN("The set value is returned.") {
                REQUIRE(config.opt_float("layer_height") == INITIAL_LAYER_HEIGHT*2);
            }
        }
    }
}

SCENARIO("Config ini load/save interface", "[Config]") {
    WHEN("new_from_ini is called") {
		Slic3r::DynamicPrintConfig config;
		std::string path = std::string(TEST_DATA_DIR) + "/test_config/new_from_ini.ini";
		config.load_from_ini(path, ForwardCompatibilitySubstitutionRule::Disable);
        THEN("Config object contains ini file options.") {
			REQUIRE(config.option_throw<ConfigOptionStrings>("filament_colour", false)->values.size() == 1);
			REQUIRE(config.option_throw<ConfigOptionStrings>("filament_colour", false)->values.front() == "#ABCD");
        }
    }
}

TEST_CASE("Flush-volume warning predicate respects used filament transitions", "[Config][Regression]")
{
    const std::vector<double> multipliers = {1.0};

    SECTION("Single used filament does not trigger warning with zero transition entries")
    {
        const std::vector<double> matrix = {
            0.0, 0.0,
            0.0, 0.0
        };
        const std::vector<int> used_filaments = {1};

        REQUIRE_FALSE(has_zero_flush_volume_for_used_filaments(matrix, multipliers, used_filaments));
    }

    SECTION("Two used filaments trigger warning when transition flush entry is zero")
    {
        const std::vector<double> matrix = {
            0.0, 0.0,
            0.0, 0.0
        };
        const std::vector<int> used_filaments = {1, 2};

        REQUIRE(has_zero_flush_volume_for_used_filaments(matrix, multipliers, used_filaments));
    }

    SECTION("Two used filaments do not trigger warning when transitions are non-zero")
    {
        const std::vector<double> matrix = {
            0.0, 280.0,
            280.0, 0.0
        };
        const std::vector<int> used_filaments = {1, 2};

        REQUIRE_FALSE(has_zero_flush_volume_for_used_filaments(matrix, multipliers, used_filaments));
    }

    SECTION("Zero multiplier still triggers warning when multiple filaments are used")
    {
        const std::vector<double> matrix = {
            0.0, 280.0,
            280.0, 0.0
        };
        const std::vector<double> zero_multiplier = {0.0};
        const std::vector<int> used_filaments = {1, 2};

        REQUIRE(has_zero_flush_volume_for_used_filaments(matrix, zero_multiplier, used_filaments));
    }
}

// TODO: https://github.com/SoftFever/OrcaSlicer/issues/11269 - Is this test still relevant? Delete if not.
// It was failing so at least "nozzle_type" and "extruder_printable_area" could not be serialized
// and an exception was thrown, but "nozzle_type" has been around for at least 3 months now.
// So maybe this test and the serialization logic in Config.?pp should be deleted if it doesn't get used.
SCENARIO("DynamicPrintConfig serialization", "[Config]") {
    WHEN("DynamicPrintConfig is serialized and deserialized") {
        FullPrintConfig full_print_config;
        DynamicPrintConfig cfg;
        cfg.apply(full_print_config, false);

        std::string serialized;
        // try {
            std::ostringstream ss;
            cereal::BinaryOutputArchive oarchive(ss);
            oarchive(cfg);
            serialized = ss.str();
        // } catch (const std::runtime_error & /* e */) {
        //     // e.what();
        // }
	CAPTURE(serialized.length());

        THEN("Config object contains ini file options.") {
            DynamicPrintConfig cfg2;
            // try {
                std::stringstream ss(serialized);
                cereal::BinaryInputArchive iarchive(ss);
                iarchive(cfg2);
            // } catch (const std::runtime_error & /* e */) {
            //     // e.what();
            // }
	    CAPTURE(cfg.diff_report(cfg2));
            REQUIRE(cfg == cfg2);
        }
    }
}

SCENARIO("update_non_diff_values_to_base_config preserves child vectors when child has more extruders than parent",
         "[Config][Variant]") {
    GIVEN("A 2-extruder child printer config inheriting from a 1-extruder parent") {
        Slic3r::DynamicPrintConfig child;
        Slic3r::DynamicPrintConfig parent;

        child.set_key_value("nozzle_diameter",           new Slic3r::ConfigOptionFloats({0.4, 0.4}));
        child.set_key_value("printer_extruder_id",       new Slic3r::ConfigOptionInts({1, 2}));
        child.set_key_value("printer_extruder_variant",  new Slic3r::ConfigOptionStrings({"Direct Drive Standard", "Direct Drive Standard"}));
        child.set_key_value("retraction_length",         new Slic3r::ConfigOptionFloats({1.5, 1.5}));

        parent.set_key_value("nozzle_diameter",          new Slic3r::ConfigOptionFloats({0.4}));
        parent.set_key_value("printer_extruder_id",      new Slic3r::ConfigOptionInts({1}));
        parent.set_key_value("printer_extruder_variant", new Slic3r::ConfigOptionStrings({"Direct Drive Standard"}));
        parent.set_key_value("retraction_length",        new Slic3r::ConfigOptionFloats({0.8}));

        const Slic3r::t_config_option_keys keys = {
            "retraction_length", "printer_extruder_id", "printer_extruder_variant"
        };
        const std::set<std::string> different_keys = {
            "retraction_length", "printer_extruder_id", "printer_extruder_variant"
        };

        WHEN("update_non_diff_values_to_base_config is called") {
            std::string id_name  = "printer_extruder_id";
            std::string var_name = "printer_extruder_variant";
            child.update_non_diff_values_to_base_config(
                parent, keys, different_keys, id_name, var_name,
                Slic3r::printer_options_with_variant_1,
                Slic3r::printer_options_with_variant_2);

            THEN("printer_extruder_id retains size 2") {
                REQUIRE(child.option<Slic3r::ConfigOptionInts>("printer_extruder_id")->values.size() == 2);
            }
            THEN("printer_extruder_variant retains size 2") {
                REQUIRE(child.option<Slic3r::ConfigOptionStrings>("printer_extruder_variant")->values.size() == 2);
            }
            THEN("retraction_length retains size 2") {
                REQUIRE(child.option<Slic3r::ConfigOptionFloats>("retraction_length")->values.size() == 2);
            }
            THEN("printer_extruder_id values are preserved for both extruders") {
                auto* pe_id = child.option<Slic3r::ConfigOptionInts>("printer_extruder_id");
                REQUIRE(pe_id->values.size() == 2);
                REQUIRE(pe_id->values[0] == 1);
                REQUIRE(pe_id->values[1] == 2);
            }
        }
    }
}

SCENARIO("update_diff_values_to_child_config tolerates legacy machine-limit vector sizes",
         "[Config][Variant]") {
    // Regression: loading a user printer preset that inherits a non-BBL multi-extruder base and
    // overrides stride-2 machine limits used to throw in ConfigOptionVector::set_only_diff
    // ("invalid diff_index size"). The base's machine-limit vectors get length-extended by the
    // nozzle count while it carries no printer_extruder_variant, so the base length (nozzles*2)
    // no longer matches variant_index.size()*2. The throw was caught upstream and DELETED the
    // user's preset file. The merge must instead degrade gracefully.
    GIVEN("A 4-nozzle parent with stride-2 limits extended to nozzles*2 but no printer_extruder_variant") {
        Slic3r::DynamicPrintConfig parent;
        Slic3r::DynamicPrintConfig child;

        parent.set_key_value("nozzle_diameter",
            new Slic3r::ConfigOptionFloats({0.4, 0.4, 0.4, 0.4}));
        parent.set_key_value("machine_max_acceleration_x",
            new Slic3r::ConfigOptionFloats({25000, 25000, 25000, 25000, 25000, 25000, 25000, 25000}));

        // Child user preset declares 4 extruder variants and overrides the machine limit.
        child.set_key_value("printer_extruder_id",
            new Slic3r::ConfigOptionInts({1, 2, 3, 4}));
        child.set_key_value("printer_extruder_variant",
            new Slic3r::ConfigOptionStrings({"Direct Drive Standard", "Direct Drive Standard",
                                             "Direct Drive Standard", "Direct Drive Standard"}));
        child.set_key_value("machine_max_acceleration_x",
            new Slic3r::ConfigOptionFloats({8000, 8000, 8000, 8000, 8000, 8000, 8000, 8000}));

        WHEN("update_diff_values_to_child_config merges the child overrides") {
            std::string id_name  = "printer_extruder_id";
            std::string var_name = "printer_extruder_variant";

            THEN("it does not throw on the legacy size mismatch") {
                REQUIRE_NOTHROW(parent.update_diff_values_to_child_config(
                    child, id_name, var_name,
                    Slic3r::printer_options_with_variant_1,
                    Slic3r::printer_options_with_variant_2));

                AND_THEN("the child's overridden machine limit is preserved") {
                    auto* mx = parent.option<Slic3r::ConfigOptionFloats>("machine_max_acceleration_x");
                    REQUIRE(mx != nullptr);
                    REQUIRE(mx->values.size() >= 2);
                    REQUIRE_THAT(mx->values[0], Catch::Matchers::WithinAbs(8000.0, 1e-6));
                    REQUIRE_THAT(mx->values[1], Catch::Matchers::WithinAbs(8000.0, 1e-6));
                }
            }
        }
    }
}

// SCENARIO("DynamicPrintConfig JSON serialization", "[Config]") {
//     WHEN("DynamicPrintConfig is serialized and deserialized") {
// 	auto now = std::chrono::high_resolution_clock::now();
// 	auto timestamp = now.time_since_epoch().count();
// 	std::stringstream ss;
// 	ss << "catch_test_serialization_" << timestamp << ".json";
// 	std::string filename = (fs::temp_directory_path() / ss.str()).string();

// TODO: Finish making a unit test for JSON serialization
//         FullPrintConfig full_print_config;
//         DynamicPrintConfig cfg;
//         cfg.apply(full_print_config, false);

//         std::string serialized;
//         try {
//             std::ostringstream ss;
//             cereal::BinaryOutputArchive oarchive(ss);
//             oarchive(cfg);
//             serialized = ss.str();
//         } catch (const std::runtime_error & /* e */) {
//             // e.what();
//         }
// 	CAPTURE(serialized.length());

//         THEN("Config object contains ini file options.") {
//             DynamicPrintConfig cfg2;
//             try {
//                 std::stringstream ss(serialized);
//                 cereal::BinaryInputArchive iarchive(ss);
//                 iarchive(cfg2);
//             } catch (const std::runtime_error & /* e */) {
//                 // e.what();
//             }
// 	    CAPTURE(cfg.diff_report(cfg2));
//             REQUIRE(cfg == cfg2);
//         }
//     }
// }

TEST_CASE("save_to_json round-trips plugin capability references as strings", "[Config][plugins]") {
    ScopedTemporaryFile tmp(".json");
    const std::vector<std::string> refs = {
        "local_plugin;;inset",
        "cloud_plugin;550e8400-e29b-41d4-a716-446655440000;inset"
    };

    std::unique_ptr<DynamicPrintConfig> config_ptr(
        DynamicPrintConfig::new_from_defaults_keys({"slicing_pipeline_plugin"}));
    DynamicPrintConfig config = std::move(*config_ptr);
    config.option<ConfigOptionStrings>("slicing_pipeline_plugin", true)->values = refs;
    config.save_to_json(tmp.string(), "test_preset", "User", "1.0.0.0");

    nlohmann::json j;
    {
        boost::nowide::ifstream ifs(tmp.string());
        ifs >> j;
    }
    REQUIRE(j["slicing_pipeline_plugin"] == nlohmann::json(refs));
    CHECK_FALSE(j.contains("plugins"));

    DynamicPrintConfig reloaded = DynamicPrintConfig::full_print_config();
    ConfigSubstitutionContext substitutions(ForwardCompatibilitySubstitutionRule::Disable);
    std::map<std::string, std::string> key_values;
    std::string reason;
    REQUIRE(reloaded.load_from_json(tmp.string(), substitutions, true, key_values, reason) == 0);
    CHECK(reason.empty());
    CHECK(reloaded.option<ConfigOptionStrings>("slicing_pipeline_plugin")->values == refs);
}

TEST_CASE("nonplanar and perimeter layering settings survive process preset JSON", "[Config][Preset][Nonplanar]")
{
    ScopedTemporaryFile tmp(".json");
    const std::vector<std::string> keys = {
        "perimeter_layering", "top_surface_z_mode", "staggered_perimeters",
        "staggered_perimeters_inner_only", "staggered_perimeter_offset",
        "smooth_outer_wall_layer_height", "interlocking_wall_amplitude",
        "interlocking_wall_wavelength", "interlocking_wall_resolution",
        "nonplanar_top_surface", "nonplanar_top_surface_max_angle",
        "nonplanar_top_surface_layers", "nonplanar_top_surface_resolution",
        "nonplanar_top_surface_min_height", "nonplanar_infill",
        "nonplanar_infill_amplitude", "nonplanar_infill_wavelength",
        "nonplanar_infill_resolution", "support_ironing_nonplanar",
        "support_ironing_nonplanar_max_angle",
        "support_ironing_nonplanar_resolution", "localized_shrinkage_strategy",
        "localized_shrinkage_infill_wall_gap",
        "localized_shrinkage_section_width",
        "localized_shrinkage_section_spacing",
        "localized_shrinkage_perforation_diameter",
        "localized_shrinkage_perforation_spacing", "zaa_enabled",
        "zaa_minimize_perimeter_height", "zaa_dont_alternate_fill_direction",
        "zaa_min_z", "bridge_overhang_before_walls"
    };
    for (const std::string &key : keys)
        CHECK(std::find(Preset::print_options().begin(), Preset::print_options().end(), key) !=
              Preset::print_options().end());

    std::unique_ptr<DynamicPrintConfig> config_ptr(
        DynamicPrintConfig::new_from_defaults_keys(keys));
    DynamicPrintConfig config = std::move(*config_ptr);
    config.set_key_value("perimeter_layering",
        new ConfigOptionEnum<PerimeterLayeringMode>(PerimeterLayeringMode::Brick));
    config.set_key_value("top_surface_z_mode",
        new ConfigOptionEnum<TopSurfaceZMode>(
            TopSurfaceZMode::NonplanarWithZContouringFallback));
    config.set_key_value("staggered_perimeters", new ConfigOptionBool(true));
    config.set_key_value("staggered_perimeters_inner_only", new ConfigOptionBool(false));
    config.set_key_value("staggered_perimeter_offset", new ConfigOptionPercent(42.));
    config.set_key_value("nonplanar_top_surface", new ConfigOptionBool(true));
    config.set_key_value("zaa_enabled", new ConfigOptionBool(true));
    config.set_key_value("nonplanar_top_surface_max_angle", new ConfigOptionFloat(37.));
    config.set_key_value("nonplanar_top_surface_layers", new ConfigOptionInt(4));
    config.set_key_value("nonplanar_top_surface_resolution", new ConfigOptionFloat(0.15));
    config.set_key_value("nonplanar_top_surface_min_height", new ConfigOptionFloat(0.04));
    config.set_key_value("nonplanar_infill", new ConfigOptionBool(true));
    config.set_key_value("nonplanar_infill_amplitude", new ConfigOptionFloat(0.07));
    config.set_key_value("nonplanar_infill_wavelength", new ConfigOptionFloat(5.));
    config.set_key_value("nonplanar_infill_resolution", new ConfigOptionFloat(0.75));
    config.save_to_json(tmp.string(), "nonplanar_process", "User", "1.0.0.0");

    nlohmann::json json;
    {
        boost::nowide::ifstream ifs(tmp.string());
        ifs >> json;
    }
    for (const std::string &key : keys)
        CHECK(json.contains(key));

    DynamicPrintConfig reloaded = DynamicPrintConfig::full_print_config();
    ConfigSubstitutionContext substitutions(ForwardCompatibilitySubstitutionRule::Disable);
    std::map<std::string, std::string> key_values;
    std::string reason;
    REQUIRE(reloaded.load_from_json(tmp.string(), substitutions, true, key_values, reason) == 0);
    CHECK(reason.empty());
    CHECK(reloaded.opt_enum<PerimeterLayeringMode>("perimeter_layering") ==
          PerimeterLayeringMode::Brick);
    CHECK(reloaded.opt_enum<TopSurfaceZMode>("top_surface_z_mode") ==
          TopSurfaceZMode::NonplanarWithZContouringFallback);
    CHECK(reloaded.opt_bool("staggered_perimeters"));
    CHECK_FALSE(reloaded.opt_bool("staggered_perimeters_inner_only"));
    CHECK_THAT(reloaded.option<ConfigOptionPercent>("staggered_perimeter_offset")->value,
               Catch::Matchers::WithinAbs(42., EPSILON));
    CHECK(reloaded.opt_bool("nonplanar_top_surface"));
    CHECK(reloaded.opt_bool("zaa_enabled"));
    CHECK(reloaded.opt_bool("nonplanar_infill"));
    CHECK_THAT(reloaded.opt_float("nonplanar_infill_amplitude"),
               Catch::Matchers::WithinAbs(0.07, EPSILON));
}

TEST_CASE("nonplanar clearance settings survive printer preset JSON", "[Config][Preset][Nonplanar]")
{
    ScopedTemporaryFile tmp(".json");
    const std::vector<std::string> keys = {
        "nonplanar_toolhead_clearance_angle", "sequential_print_gantry_geometry",
        "sequential_print_gantry_model"
    };
    for (const std::string &key : keys)
        CHECK(std::find(Preset::printer_options().begin(),
                        Preset::printer_options().end(), key) !=
              Preset::printer_options().end());

    std::unique_ptr<DynamicPrintConfig> config_ptr(
        DynamicPrintConfig::new_from_defaults_keys(keys));
    DynamicPrintConfig config = std::move(*config_ptr);
    config.set_key_value("nonplanar_toolhead_clearance_angle",
                         new ConfigOptionFloat(42.));
    config.save_to_json(tmp.string(), "nonplanar_printer", "User", "1.0.0.0");

    nlohmann::json json;
    {
        boost::nowide::ifstream ifs(tmp.string());
        ifs >> json;
    }
    for (const std::string &key : keys)
        CHECK(json.contains(key));

    DynamicPrintConfig reloaded = DynamicPrintConfig::full_print_config();
    ConfigSubstitutionContext substitutions(
        ForwardCompatibilitySubstitutionRule::Disable);
    std::map<std::string, std::string> key_values;
    std::string reason;
    REQUIRE(reloaded.load_from_json(tmp.string(), substitutions, true,
                                    key_values, reason) == 0);
    CHECK(reason.empty());
    CHECK_THAT(reloaded.opt_float("nonplanar_toolhead_clearance_angle"),
               Catch::Matchers::WithinAbs(42., EPSILON));
}

TEST_CASE("save_to_json writes the same document to a stream as to a file", "[Config]") {
    DynamicPrintConfig config;
    config.set_key_value("layer_height", new ConfigOptionFloat(0.2));
    config.set_key_value("wall_loops", new ConfigOptionInt(3));
    config.set_key_value("filament_type", new ConfigOptionStrings({ "PLA", "PETG" }));
    config.set_key_value("machine_start_gcode", new ConfigOptionString("G28\nG1 Z5"));

    ScopedTemporaryFile tmp(".json");
    config.save_to_json(tmp.string(), "test_preset", "User", "1.0.0.0");
    std::string file_contents;
    {
        boost::nowide::ifstream ifs(tmp.string());
        file_contents.assign(std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>());
    }
    // The file format: one tab per nesting level and a trailing newline.
    REQUIRE_FALSE(file_contents.empty());
    CHECK(file_contents.rfind("{\n\t\"", 0) == 0);
    CHECK(file_contents.back() == '\n');

    std::ostringstream strict, replaced;
    config.save_to_json(strict, "test_preset", "User", "1.0.0.0");
    config.save_to_json(replaced, "test_preset", "User", "1.0.0.0", true);
    CHECK(strict.str() == file_contents);
    CHECK(replaced.str() == file_contents);
    CHECK(nlohmann::json::parse(strict.str())["machine_start_gcode"] == "G28\nG1 Z5");
}

TEST_CASE("save_to_json replaces invalid UTF-8 in a stream only when asked", "[Config]") {
    DynamicPrintConfig config;
    config.set_key_value("machine_start_gcode", new ConfigOptionString("G28 ; \xff"));

    std::ostringstream strict, replaced;
    CHECK_THROWS_AS(config.save_to_json(strict, "test_preset", "User", "1.0.0.0"), nlohmann::json::type_error);
    REQUIRE_NOTHROW(config.save_to_json(replaced, "test_preset", "User", "1.0.0.0", true));
    CHECK(nlohmann::json::parse(replaced.str())["machine_start_gcode"] == "G28 ; \xEF\xBF\xBD");
}

TEST_CASE("save_to_json leaves an existing file untouched when the config cannot be serialized", "[Config]") {
    DynamicPrintConfig config;
    config.set_key_value("machine_start_gcode", new ConfigOptionString("G28 ; \xff"));

    ScopedTemporaryFile tmp(".json");
    {
        boost::nowide::ofstream ofs(tmp.string());
        ofs << "previous";
    }
    CHECK_THROWS_AS(config.save_to_json(tmp.string(), "test_preset", "User", "1.0.0.0"), nlohmann::json::type_error);

    boost::nowide::ifstream ifs(tmp.string());
    const std::string contents((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    CHECK(contents == "previous");
}

TEST_CASE("plugin capability references survive string-map serialization", "[Config][plugins]") {
    const std::vector<std::string> refs = {
        "master_plugin;;header-stamp",
        "Sample Plugin;1f998ea9-0183-4cc5-957f-4eef659ba4e6;G-code Benchmark (.py)"
    };

    DynamicPrintConfig original = DynamicPrintConfig::full_print_config();
    original.option<ConfigOptionStrings>("slicing_pipeline_plugin", true)->values = refs;

    std::map<std::string, std::string> serialized{
        {"slicing_pipeline_plugin", original.option<ConfigOptionStrings>("slicing_pipeline_plugin")->serialize()}
    };
    CHECK(serialized["slicing_pipeline_plugin"].find("\"master_plugin;;header-stamp\"") != std::string::npos);

    DynamicPrintConfig reloaded = DynamicPrintConfig::full_print_config();
    reloaded.load_string_map(serialized, ForwardCompatibilitySubstitutionRule::Disable);

    CHECK(reloaded.option<ConfigOptionStrings>("slicing_pipeline_plugin")->values == refs);
}

TEST_CASE("parse_capability_ref parses local and cloud references", "[Config][plugin]") {
    const auto local = Slic3r::parse_capability_ref("local_plugin;;post_process");
    REQUIRE(local.has_value());
    CHECK(local->name == "local_plugin");
    CHECK(local->capability_name == "post_process");
    CHECK(local->uuid.empty());

    const auto cloud = Slic3r::parse_capability_ref(
        "cloud_plugin;550e8400-e29b-41d4-a716-446655440000;post_process");
    REQUIRE(cloud.has_value());
    CHECK(cloud->name == "cloud_plugin");
    CHECK(cloud->capability_name == "post_process");
    CHECK(cloud->uuid == "550e8400-e29b-41d4-a716-446655440000");
}

TEST_CASE("parse_capability_ref rejects malformed input", "[Config][plugin]") {
    CHECK_FALSE(Slic3r::parse_capability_ref("").has_value());
    CHECK_FALSE(Slic3r::parse_capability_ref("plugin").has_value());
    CHECK_FALSE(Slic3r::parse_capability_ref("plugin;uuid").has_value());
    CHECK_FALSE(Slic3r::parse_capability_ref(";;capability").has_value());
    CHECK_FALSE(Slic3r::parse_capability_ref(";uuid;capability").has_value());
    CHECK_FALSE(Slic3r::parse_capability_ref("plugin;;").has_value());
    CHECK_FALSE(Slic3r::parse_capability_ref("plugin;uuid;").has_value());
}

namespace {
// Installs a stub capability resolver that echoes the capability type into the reference, so tests
// can assert each plugin-backed option resolved with its own ConfigOptionDef::plugin_type. Resets
// the global resolver on teardown -- tests run in random order and other cases assert the
// no-resolver behavior (an absent "plugins" manifest).
struct PluginResolverFixture {
    PluginResolverFixture() {
        ConfigBase::set_resolve_capability_fn([](const std::string& name, const std::string& type) {
            return name.empty() ? std::string() : name + ";;" + type;
        });
    }
    ~PluginResolverFixture() { ConfigBase::set_resolve_capability_fn(nullptr); }
};
} // namespace

TEST_CASE_METHOD(PluginResolverFixture,
    "update_plugin_manifest derives references generically from plugin-backed options",
    "[Config][plugins]") {
    // Both scalar (printer_agent) and vector (slicing_pipeline_plugin) options opt in via a non-empty
    // ConfigOptionDef::plugin_type (is_plugin_backed) and are resolved with it -- there is no hardcoded
    // per-option switch. printer_agent in particular relies on its plugin_type metadata being wired up
    // (it is edited via a dedicated widget, not the plugin_picker).
    std::unique_ptr<DynamicPrintConfig> config_ptr(DynamicPrintConfig::new_from_defaults_keys(
        {"slicing_pipeline_plugin", "printer_agent"}));
    DynamicPrintConfig config = std::move(*config_ptr);
    config.option<ConfigOptionStrings>("slicing_pipeline_plugin", true)->values = {"sp"};
    config.option<ConfigOptionString>("printer_agent", true)->value            = "agent";

    config.update_plugin_manifest();
    const std::vector<std::string> manifest = config.option<ConfigOptionStrings>("plugins")->values;

    using Catch::Matchers::VectorContains;
    REQUIRE_THAT(manifest, VectorContains(std::string("sp;;slicing-pipeline")));
    REQUIRE_THAT(manifest, VectorContains(std::string("agent;;printer-connection")));
    CHECK(manifest.size() == 2);
}

TEST_CASE_METHOD(PluginResolverFixture,
    "update_plugin_manifest de-duplicates references and skips unset options",
    "[Config][plugins]") {
    std::unique_ptr<DynamicPrintConfig> config_ptr(DynamicPrintConfig::new_from_defaults_keys(
        {"slicing_pipeline_plugin", "printer_agent"}));
    DynamicPrintConfig config = std::move(*config_ptr);
    config.option<ConfigOptionStrings>("slicing_pipeline_plugin", true)->values = {"x", "x"};  // duplicate
    // printer_agent stays at its default empty value -> contributes nothing to the manifest.

    config.update_plugin_manifest();
    const std::vector<std::string> manifest = config.option<ConfigOptionStrings>("plugins")->values;

    CHECK(manifest == std::vector<std::string>{"x;;slicing-pipeline"});
}

TEST_CASE("H2C/A2L-era multi-nozzle and pre-heat config keys exist", "[config]") {
    // Foundation keys backing H2C 6-nozzle cluster grouping, the pre-heat/pre-cool time
    // model, and wipe-tower nozzle-change handling. Defaults must keep existing
    // single-nozzle printers behaving identically.
    Slic3r::DynamicPrintConfig config = Slic3r::DynamicPrintConfig::full_print_config();

    // Printer / per-extruder options
    REQUIRE(config.option<ConfigOptionIntsNullable>("extruder_max_nozzle_count") != nullptr);
    REQUIRE(config.option<ConfigOptionIntsNullable>("extruder_max_nozzle_count")->values == std::vector<int>{1});
    REQUIRE(config.option<ConfigOptionBool>("enable_pre_heating") != nullptr);
    REQUIRE(config.option<ConfigOptionBool>("enable_pre_heating")->value == false);
    REQUIRE(config.option<ConfigOptionFloatsNullable>("hotend_cooling_rate") != nullptr);
    REQUIRE(config.option<ConfigOptionFloatsNullable>("hotend_heating_rate") != nullptr);
    REQUIRE(config.option<ConfigOptionFloat>("machine_hotend_change_time") != nullptr);
    REQUIRE(config.option<ConfigOptionFloat>("machine_prepare_compensation_time") != nullptr);

    // Filament pre-cooling / ramming / nozzle-change (nc) options
    REQUIRE(config.option<ConfigOptionIntsNullable>("filament_pre_cooling_temperature") != nullptr);
    REQUIRE(config.option<ConfigOptionIntsNullable>("filament_pre_cooling_temperature_nc") != nullptr);
    REQUIRE(config.option<ConfigOptionFloatsNullable>("filament_preheat_temperature_delta") != nullptr);
    REQUIRE(config.option<ConfigOptionFloatsNullable>("filament_retract_length_nc") != nullptr);
    REQUIRE(config.option<ConfigOptionFloats>("filament_change_length_nc") != nullptr);
    REQUIRE(config.option<ConfigOptionFloats>("filament_prime_volume_nc") != nullptr);
    REQUIRE(config.option<ConfigOptionFloatsNullable>("filament_ramming_travel_time") != nullptr);
    REQUIRE(config.option<ConfigOptionFloatsNullable>("filament_ramming_travel_time_nc") != nullptr);
    REQUIRE(config.option<ConfigOptionFloatsNullable>("filament_ramming_volumetric_speed") != nullptr);
    REQUIRE(config.option<ConfigOptionFloatsNullable>("filament_ramming_volumetric_speed_nc") != nullptr);

    // Spot-check defaults that must not alter existing behavior.
    REQUIRE(config.option<ConfigOptionFloatsNullable>("filament_retract_length_nc")->values == std::vector<double>{10.});
    REQUIRE(config.option<ConfigOptionFloats>("filament_prime_volume_nc")->values == std::vector<double>{60.});
    REQUIRE(config.option<ConfigOptionIntsNullable>("filament_pre_cooling_temperature_nc")->values == std::vector<int>{0});
    REQUIRE(config.option<ConfigOptionFloatsNullable>("filament_ramming_volumetric_speed")->values == std::vector<double>{-1});
}

SCENARIO("ConfigOptionVector::set_to_index with stride=1 copies values correctly", "[Config][set_to_index]") {
    GIVEN("A destination vector and a source vector with 3 values") {
        Slic3r::ConfigOptionFloats dest({0.0});
        Slic3r::ConfigOptionFloats src({10.0, 20.0, 30.0});
        std::vector<int> variant_index = {0, 1, 2};
        int stride = 1;

        WHEN("set_to_index is called with stride=1") {
            dest.set_to_index(&src, variant_index, stride);

            THEN("The destination contains the source values") {
                REQUIRE(dest.values.size() == 3);
                REQUIRE(dest.values[0] == 10.0);
                REQUIRE(dest.values[1] == 20.0);
                REQUIRE(dest.values[2] == 30.0);
            }
        }
    }

    GIVEN("A destination vector and a source vector with subset mapping") {
        Slic3r::ConfigOptionFloats dest({0.0});
        Slic3r::ConfigOptionFloats src({100.0, 200.0, 300.0});
        std::vector<int> variant_index = {1, 2};
        int stride = 1;

        WHEN("set_to_index maps only indices 1 and 2") {
            dest.set_to_index(&src, variant_index, stride);

            THEN("Only the mapped values are copied, default fills the others") {
                REQUIRE(dest.values.size() == 2);
                REQUIRE(dest.values[0] == 200.0);
                REQUIRE(dest.values[1] == 300.0);
            }
        }
    }
}

SCENARIO("ConfigOptionVector::set_to_index with stride=2 copies grouped values correctly", "[Config][set_to_index]") {
    GIVEN("A destination vector and a source vector with stride=2 (e.g., nozzle groups)") {
        // Source has 4 groups of 2 values each: (10,11), (20,21), (30,31), (40,41)
        Slic3r::ConfigOptionFloats dest({0.0});
        Slic3r::ConfigOptionFloats src({10.0, 11.0, 20.0, 21.0, 30.0, 31.0, 40.0, 41.0});
        int stride = 2;

        WHEN("set_to_index maps groups 0, 1, 3") {
            std::vector<int> variant_index = {0, 1, 3};
            dest.set_to_index(&src, variant_index, stride);

            THEN("The destination has 3 groups (6 values) mapped correctly") {
                REQUIRE(dest.values.size() == 6);
                // Group 0: (10, 11)
                REQUIRE(dest.values[0] == 10.0);
                REQUIRE(dest.values[1] == 11.0);
                // Group 1: (20, 21)
                REQUIRE(dest.values[2] == 20.0);
                REQUIRE(dest.values[3] == 21.0);
                // Group 3: (40, 41)
                REQUIRE(dest.values[4] == 40.0);
                REQUIRE(dest.values[5] == 41.0);
            }
        }
    }

    GIVEN("A destination and a single-group source") {
        Slic3r::ConfigOptionFloats dest({0.0});
        // Source has 1 group of 2 values
        Slic3r::ConfigOptionFloats src({50.0, 60.0});
        int stride = 2;

        WHEN("set_to_index maps group 0 from a single-group source") {
            std::vector<int> variant_index = {0};
            dest.set_to_index(&src, variant_index, stride);

            THEN("The destination contains the single group correctly") {
                REQUIRE(dest.values.size() == 2);
                REQUIRE(dest.values[0] == 50.0);
                REQUIRE(dest.values[1] == 60.0);
            }
        }
    }
}

SCENARIO("ConfigOptionVector::set_to_index handles empty dest_index", "[Config][set_to_index]") {
    GIVEN("A destination and source with stride=2") {
        Slic3r::ConfigOptionFloats dest({0.0});
        Slic3r::ConfigOptionFloats src({10.0, 11.0, 20.0, 21.0});
        std::vector<int> variant_index = {};
        int stride = 2;

        WHEN("set_to_index is called with an empty index vector") {
            dest.set_to_index(&src, variant_index, stride);

            THEN("The destination is resized to 0") {
                REQUIRE(dest.values.size() == 0);
            }
        }
    }
}

SCENARIO("ConfigOptionVector::set_to_index handles nil values in source", "[Config][set_to_index]") {
    GIVEN("A source with a nil group (stride=2)") {
        Slic3r::ConfigOptionFloatsNullable dest({0.0});
        Slic3r::ConfigOptionFloatsNullable src({10.0, 11.0,
            Slic3r::ConfigOptionFloatsNullable::nil_value(), Slic3r::ConfigOptionFloatsNullable::nil_value(),
            30.0, 31.0});
        int stride = 2;

        WHEN("set_to_index maps all groups including the nil one") {
            std::vector<int> variant_index = {0, 1, 2};
            dest.set_to_index(&src, variant_index, stride);

            THEN("Non-nil groups are copied and the nil group keeps the default") {
                REQUIRE(dest.values.size() == 6);
                // Group 0: (10, 11) — copied
                REQUIRE(dest.values[0] == 10.0);
                REQUIRE(dest.values[1] == 11.0);
                // Group 1: nil — keeps default (the front value = 10.0)
                REQUIRE(dest.values[2] == 10.0);
                REQUIRE(dest.values[3] == 10.0);
                // Group 2: (30, 31) — copied
                REQUIRE(dest.values[4] == 30.0);
                REQUIRE(dest.values[5] == 31.0);
            }
        }
    }
}

SCENARIO("ConfigOptionVector::set_to_index handles out-of-bounds dest_index", "[Config][set_to_index]") {
    GIVEN("A source with only 2 groups (4 values) but dest_index references group 3") {
        Slic3r::ConfigOptionFloats dest({0.0});
        Slic3r::ConfigOptionFloats src({10.0, 11.0, 20.0, 21.0}); // 2 groups of stride 2
        int stride = 2;

        WHEN("set_to_index maps group 3 which is out of bounds") {
            std::vector<int> variant_index = {0, 3}; // group 3 is out of range
            dest.set_to_index(&src, variant_index, stride);

            THEN("Group 0 is copied, group 3 falls back to default without crashing") {
                REQUIRE(dest.values.size() == 4);
                // Group 0: (10, 11) — copied
                REQUIRE(dest.values[0] == 10.0);
                REQUIRE(dest.values[1] == 11.0);
                // Group 3: out of bounds — keeps default (10.0 = src.values.front())
                REQUIRE(dest.values[2] == 10.0);
                REQUIRE(dest.values[3] == 10.0);
            }
        }
    }
}

SCENARIO("ConfigOptionVector::set_to_index handles negative dest_index values", "[Config][set_to_index]") {
    GIVEN("A destination and source with a negative entry in dest_index") {
        // The dest is initially empty, so resize fills all slots with src.values.front().
        Slic3r::ConfigOptionFloats dest;
        Slic3r::ConfigOptionFloats src({100.0, 101.0, 200.0, 201.0});
        int stride = 2;

        WHEN("set_to_index maps group 0 and a negative index") {
            std::vector<int> variant_index = {-1, 0};
            dest.set_to_index(&src, variant_index, stride);

            THEN("The negative index is skipped, the valid group is copied") {
                REQUIRE(dest.values.size() == 4);
                // Position 0 (variant_index[0] = -1): skipped, keeps default fill
                // from resize (src.values.front() = 100.0, applied to all new elements)
                REQUIRE(dest.values[0] == 100.0);
                REQUIRE(dest.values[1] == 100.0);
                // Position 1 (variant_index[1] = 0): copied from group 0 of src
                REQUIRE(dest.values[2] == 100.0);
                REQUIRE(dest.values[3] == 101.0);
            }
        }
    }
}

SCENARIO("ConfigOptionVector::set_to_index handles single-element groups with stride=1", "[Config][set_to_index]") {
    GIVEN("A destination re-mapping one variant index with a stride=1 source") {
        // Simulates the PrintObject.cpp code path: stride=1, variant_index={1}
        Slic3r::ConfigOptionFloats dest({99.0, 99.0, 99.0, 99.0}); // pre-sized for 4 extruders
        Slic3r::ConfigOptionFloats src({0.5, 0.6, 0.7, 0.8});      // 4 extruder values
        std::vector<int> variant_index = {1}; // only extruder 1 is active
        int stride = 1;

        WHEN("set_to_index is called") {
            dest.set_to_index(&src, variant_index, stride);

            THEN("Only the mapped value is copied, rest are defaulted") {
                REQUIRE(dest.values.size() == 1);
                REQUIRE(dest.values[0] == 0.6);
            }
        }
    }
}

SCENARIO("ConfigOptionVector::set_to_index throws on incompatible type", "[Config][set_to_index]") {
    GIVEN("A Floats destination and an Ints source") {
        Slic3r::ConfigOptionFloats dest({0.0});
        Slic3r::ConfigOptionInts src({1, 2, 3});
        std::vector<int> variant_index = {0};
        int stride = 1;

        WHEN("set_to_index is called with mismatched types") {
            THEN("A ConfigurationError is thrown") {
                REQUIRE_THROWS_AS(dest.set_to_index(&src, variant_index, stride), Slic3r::ConfigurationError);
            }
        }
    }
}

TEST_CASE("read_cli applies valid values and collects non-option arguments", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--nozzle-temperature", "210,190", "--reduce-crossing-wall=1", "model.3mf"};
    REQUIRE(config.read_cli(5, argv, &extra, &keys));
    REQUIRE(config.opt<ConfigOptionInts>("nozzle_temperature")->values == std::vector<int>{210, 190});
    REQUIRE(config.opt<ConfigOptionBool>("reduce_crossing_wall")->value);
    REQUIRE(extra == t_config_option_keys{"model.3mf"});
    REQUIRE(keys == t_config_option_keys{"nozzle_temperature", "reduce_crossing_wall"});
}

TEST_CASE("read_cli rejects nil for a non-nullable vector option", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--nozzle-temperature", "nil"};
    REQUIRE_FALSE(config.read_cli(3, argv, &extra, &keys));
}

TEST_CASE("read_cli rejects an invalid boolean value", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--reduce-crossing-wall=maybe"};
    REQUIRE_FALSE(config.read_cli(2, argv, &extra, &keys));
}

TEST_CASE("read_cli accepts the common spellings of a boolean value", "[Config]") {
    const auto [text, expected] = GENERATE(table<const char*, bool>({
        {"--reduce-crossing-wall=1", true},
        {"--reduce-crossing-wall=true", true},
        {"--reduce-crossing-wall=Yes", true},
        {"--reduce-crossing-wall=on", true},
        {"--reduce-crossing-wall=enabled", true},
        {"--reduce-crossing-wall=TRUE", true},
        {"--reduce-crossing-wall=oN", true},
        {"--reduce-crossing-wall=0", false},
        {"--reduce-crossing-wall=false", false},
        {"--reduce-crossing-wall=No", false},
        {"--reduce-crossing-wall=off", false},
        {"--reduce-crossing-wall=disabled", false},
        {"--reduce-crossing-wall=FALSE", false},
        {"--reduce-crossing-wall=DiSaBlEd", false},
    }));

    DYNAMIC_SECTION(text) {
        Slic3r::DynamicPrintConfig config;
        t_config_option_keys extra, keys;
        const char* argv[] = {"orca-slicer", text};
        REQUIRE(config.read_cli(2, argv, &extra, &keys));
        REQUIRE(config.opt<ConfigOptionBool>("reduce_crossing_wall")->value == expected);
    }
}

TEST_CASE("read_cli accepts the common boolean spellings inside a bools vector", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--filament-soluble=true,no,1"};
    REQUIRE(config.read_cli(2, argv, &extra, &keys));
    REQUIRE(config.opt<ConfigOptionBools>("filament_soluble")->values == std::vector<unsigned char>{1, 0, 1});
}

TEST_CASE("read_cli trims whitespace around boolean spellings", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--reduce-crossing-wall= true ", "--filament-soluble= true , no ,1"};
    REQUIRE(config.read_cli(3, argv, &extra, &keys));
    REQUIRE(config.opt<ConfigOptionBool>("reduce_crossing_wall")->value);
    REQUIRE(config.opt<ConfigOptionBools>("filament_soluble")->values == std::vector<unsigned char>{1, 0, 1});
}

TEST_CASE("read_cli normalizes boolean spellings when a bools vector is repeated", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--filament-soluble=true", "--filament-soluble=off"};
    REQUIRE(config.read_cli(3, argv, &extra, &keys));
    REQUIRE(config.opt<ConfigOptionBools>("filament_soluble")->values == std::vector<unsigned char>{1, 0});
}

TEST_CASE("read_cli keeps nil alongside boolean spellings in a nullable bools vector", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--enable-overhang-speed=nil,yes,off"};
    REQUIRE(config.read_cli(2, argv, &extra, &keys));
    auto* opt = config.opt<ConfigOptionBoolsNullable>("enable_overhang_speed");
    REQUIRE(opt != nullptr);
    REQUIRE(opt->values.size() == 3);
    REQUIRE(opt->is_nil(0));
    REQUIRE(opt->values[1] == 1);
    REQUIRE(opt->values[2] == 0);
}

TEST_CASE("read_cli rejects an empty item inside a bools vector", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--filament-soluble=true,,1"};
    REQUIRE_FALSE(config.read_cli(2, argv, &extra, &keys));
}

TEST_CASE("read_cli rejects an unknown spelling next to a valid one in a bools vector", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--filament-soluble=true,affirmative"};
    REQUIRE_FALSE(config.read_cli(2, argv, &extra, &keys));
}

// The normalization lives in read_cli's boolean branches, so options of other types keep the
// value verbatim - a path named "on" or a colour named "true" must not turn into "1".
TEST_CASE("read_cli leaves boolean spellings alone for non-boolean options", "[Config]") {
    SECTION("string option") {
        Slic3r::DynamicPrintAndCLIConfig config;
        t_config_option_keys extra, keys;
        const char* argv[] = {"orca-slicer", "--logfile=true"};
        REQUIRE(config.read_cli(2, argv, &extra, &keys));
        REQUIRE(config.opt<ConfigOptionString>("logfile")->value == "true");
    }
    SECTION("strings vector option") {
        Slic3r::DynamicPrintConfig config;
        t_config_option_keys extra, keys;
        const char* argv[] = {"orca-slicer", "--filament-colour=on;off"};
        REQUIRE(config.read_cli(2, argv, &extra, &keys));
        REQUIRE(config.opt<ConfigOptionStrings>("filament_colour")->values == std::vector<std::string>{"on", "off"});
    }
}

TEST_CASE("read_cli treats a bare boolean flag as true without consuming the next argument", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--reduce-crossing-wall", "model.3mf"};
    REQUIRE(config.read_cli(3, argv, &extra, &keys));
    REQUIRE(config.opt<ConfigOptionBool>("reduce_crossing_wall")->value);
    REQUIRE(extra == t_config_option_keys{"model.3mf"});
}

TEST_CASE("read_cli rejects an invalid scalar numeric value", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--top-shell-layers", "several"};
    REQUIRE_FALSE(config.read_cli(3, argv, &extra, &keys));
}

TEST_CASE("read_cli appends values when a vector option is repeated", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--nozzle-temperature", "210", "--nozzle-temperature", "190,200"};
    REQUIRE(config.read_cli(5, argv, &extra, &keys));
    REQUIRE(config.opt<ConfigOptionInts>("nozzle_temperature")->values == std::vector<int>{210, 190, 200});
    // the key is recorded once, on first use
    REQUIRE(keys == t_config_option_keys{"nozzle_temperature"});
}

TEST_CASE("read_cli parses a bools vector given in the --flag=values form", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--filament-soluble=1,0,1"};
    REQUIRE(config.read_cli(2, argv, &extra, &keys));
    REQUIRE(config.opt<ConfigOptionBools>("filament_soluble")->values == std::vector<unsigned char>{1, 0, 1});
}

TEST_CASE("read_cli rejects an invalid value inside a bools vector", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--filament-soluble=1,maybe"};
    REQUIRE_FALSE(config.read_cli(2, argv, &extra, &keys));
}

TEST_CASE("read_cli appends true for a bare bools vector flag", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--filament-soluble"};
    REQUIRE(config.read_cli(2, argv, &extra, &keys));
    REQUIRE(config.opt<ConfigOptionBools>("filament_soluble")->values == std::vector<unsigned char>{1});
}

TEST_CASE("read_cli splits a strings vector on semicolons and unescapes quoted items", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--filament-colour", "#FF0000;\"a\\nb\";#00FF00"};
    REQUIRE(config.read_cli(3, argv, &extra, &keys));
    auto& values = config.opt<ConfigOptionStrings>("filament_colour")->values;
    REQUIRE(values == std::vector<std::string>{"#FF0000", "a\nb", "#00FF00"});
}

TEST_CASE("read_cli rejects a strings vector with an unterminated quote", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--filament-colour", "\"oops"};
    REQUIRE_FALSE(config.read_cli(3, argv, &extra, &keys));
}

TEST_CASE("read_cli parses a points vector in the NxM coordinate form", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--printable-area", "0x0,200x0,200x200,0x200"};
    REQUIRE(config.read_cli(3, argv, &extra, &keys));
    auto& points = config.opt<ConfigOptionPoints>("printable_area")->values;
    REQUIRE(points.size() == 4);
    REQUIRE_THAT(points[1].x(), Catch::Matchers::WithinAbs(200.0, 1e-9));
    REQUIRE_THAT(points[1].y(), Catch::Matchers::WithinAbs(0.0, 1e-9));
    REQUIRE_THAT(points[3].x(), Catch::Matchers::WithinAbs(0.0, 1e-9));
    REQUIRE_THAT(points[3].y(), Catch::Matchers::WithinAbs(200.0, 1e-9));
}

// logfile is a CLI-only option, so it needs the config type whose def pulls in cli_misc_config_def.
TEST_CASE("read_cli stores the log file path as a string", "[Config]") {
    Slic3r::DynamicPrintAndCLIConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--logfile", "orca.log"};
    REQUIRE(config.read_cli(3, argv, &extra, &keys));
    REQUIRE(config.opt<ConfigOptionString>("logfile")->value == "orca.log");
}

TEST_CASE("read_cli accepts nil entries for a nullable vector option", "[Config]") {
    Slic3r::DynamicPrintConfig config;
    t_config_option_keys extra, keys;
    const char* argv[] = {"orca-slicer", "--filament-retraction-length", "nil,2.5"};
    REQUIRE(config.read_cli(3, argv, &extra, &keys));
    auto* opt = config.opt<ConfigOptionFloatsNullable>("filament_retraction_length");
    REQUIRE(opt != nullptr);
    REQUIRE(opt->values.size() == 2);
    REQUIRE(opt->is_nil(0));
    REQUIRE_FALSE(opt->is_nil(1));
    REQUIRE_THAT(opt->values[1], Catch::Matchers::WithinAbs(2.5, 1e-9));
}

// get_at() returns values.front() for an out-of-range index, so calling it on an empty vector
// option is UB. filament_id and filament_is_support are unpopulated on a CLI from-scratch slice.
TEST_CASE("get_filament_type treats empty vector options as absent", "[Config][Filament]")
{
    DynamicPrintConfig config;
    std::string displayed;

    SECTION("an empty filament_type yields no type at all")
    {
        config.set_key_value("filament_type", new ConfigOptionStrings());
        REQUIRE(config.get_filament_type(displayed, 0) == "");
    }

    SECTION("an empty filament_is_support falls back to the plain filament type")
    {
        config.set_key_value("filament_type", new ConfigOptionStrings({"PETG"}));
        config.set_key_value("filament_is_support", new ConfigOptionBools());
        REQUIRE(config.get_filament_type(displayed, 0) == "PETG");
        REQUIRE(displayed == "PETG");
    }

    SECTION("a support filament with an empty filament_id resolves from the type alone")
    {
        config.set_key_value("filament_type", new ConfigOptionStrings({"PLA"}));
        config.set_key_value("filament_is_support", new ConfigOptionBools({true}));
        config.set_key_value("filament_id", new ConfigOptionStrings());
        REQUIRE(config.get_filament_type(displayed, 0) == "PLA-S");
        REQUIRE(displayed == "Sup.PLA");
    }

    SECTION("a populated filament_id still selects the support type by id")
    {
        config.set_key_value("filament_type", new ConfigOptionStrings({"PETG"}));
        config.set_key_value("filament_is_support", new ConfigOptionBools({true}));
        config.set_key_value("filament_id", new ConfigOptionStrings({"GFS00"}));
        REQUIRE(config.get_filament_type(displayed, 0) == "PLA-S");
        REQUIRE(displayed == "Sup.PLA");
    }
}

namespace {

// min_object_distance reads exactly these three options.
DynamicPrintConfig spacing_config(PrinterTechnology tech, PrintSequence seq, double clearance_radius)
{
    DynamicPrintConfig c;
    c.set_key_value("printer_technology", new ConfigOptionEnum<PrinterTechnology>(tech));
    c.set_key_value("print_sequence", new ConfigOptionEnum<PrintSequence>(seq));
    c.set_key_value("extruder_clearance_radius", new ConfigOptionFloat(clearance_radius));
    return c;
}

} // namespace

TEST_CASE("min_object_distance floors object spacing per print sequence", "[Config]")
{
    struct Case
    {
        std::string       description;
        PrinterTechnology tech;
        PrintSequence     sequence;
        double            clearance_radius;
        double            expected;
    };

    auto c = GENERATE(values<Case>({
        {"sequential FFF takes a clearance radius above the floor", ptFFF, PrintSequence::ByObject, 12., 12.},
        {"sequential FFF holds the floor at the radius",            ptFFF, PrintSequence::ByObject,  6.,  6.},
        {"sequential FFF holds the floor below the radius",         ptFFF, PrintSequence::ByObject,  4.,  6.},
        {"layered FFF ignores the clearance radius",                ptFFF, PrintSequence::ByLayer,  12.,  6.},
        {"SLA is a flat 6mm",                                       ptSLA, PrintSequence::ByObject, 12.,  6.},
        {"SLA ignores the print sequence too",                      ptSLA, PrintSequence::ByLayer,  12.,  6.},
    }));

    DYNAMIC_SECTION(c.description)
    {
        CHECK_THAT(min_object_distance(spacing_config(c.tech, c.sequence, c.clearance_radius)),
                   Catch::Matchers::WithinAbs(c.expected, 1e-9));
    }
}

TEST_CASE("min_object_distance yields no floor when an FFF config lacks the options", "[Config]")
{
    // Missing options yield 0 rather than an error, so a caller gets no floor at all.
    SECTION("no clearance radius") {
        DynamicPrintConfig c;
        c.set_key_value("printer_technology", new ConfigOptionEnum<PrinterTechnology>(ptFFF));
        c.set_key_value("print_sequence", new ConfigOptionEnum<PrintSequence>(PrintSequence::ByObject));
        CHECK_THAT(min_object_distance(c), Catch::Matchers::WithinAbs(0., 1e-9));
    }

    SECTION("no print sequence") {
        DynamicPrintConfig c;
        c.set_key_value("printer_technology", new ConfigOptionEnum<PrinterTechnology>(ptFFF));
        c.set_key_value("extruder_clearance_radius", new ConfigOptionFloat(12.));
        CHECK_THAT(min_object_distance(c), Catch::Matchers::WithinAbs(0., 1e-9));
    }

    SECTION("nothing at all") {
        CHECK_THAT(min_object_distance(DynamicPrintConfig{}), Catch::Matchers::WithinAbs(0., 1e-9));
    }

    SECTION("an unset printer technology is treated as FFF") {
        DynamicPrintConfig c;
        c.set_key_value("print_sequence", new ConfigOptionEnum<PrintSequence>(PrintSequence::ByObject));
        c.set_key_value("extruder_clearance_radius", new ConfigOptionFloat(12.));
        CHECK_THAT(min_object_distance(c), Catch::Matchers::WithinAbs(12., 1e-9));
    }
}

TEST_CASE("Static print configs compare, order and hash by their option values", "[Config]")
{
    // PrintObjectConfig comes from PRINT_CONFIG_CLASS_DEFINE; PrintConfig combines MachineEnvelopeConfig
    // and GCodeConfig through PRINT_CONFIG_CLASS_DERIVED_DEFINE. Both generate hash(), operator==,
    // operator< and the option registration from the same option list. The hash inequalities use fixed
    // inputs, so they are deterministic; they check that hash() covers the changed option.
    SECTION("default-constructed configs are equal and find their options by key")
    {
        PrintObjectConfig a, b;
        REQUIRE(a == b);
        REQUIRE(a.hash() == b.hash());
        REQUIRE_FALSE(a < b);
        REQUIRE_FALSE(b < a);
        REQUIRE(a.optptr("layer_height") == &a.layer_height);
        REQUIRE(a.optptr("brim_object_gap") == &a.brim_object_gap);
    }

    SECTION("one differing option makes the configs unequal and orders them")
    {
        PrintObjectConfig a, b;
        b.layer_height.value = a.layer_height.value + 0.05;
        REQUIRE(a != b);
        REQUIRE(a.hash() != b.hash());
        REQUIRE(a < b);
        REQUIRE_FALSE(b < a);
    }

    SECTION("ordering is decided by the first option in declaration order that differs")
    {
        PrintObjectConfig a, b;
        a.brim_object_gap.value = b.brim_object_gap.value + 1.0;  // declared first
        a.layer_height.value    = b.layer_height.value - 0.05;    // declared later, points the other way
        REQUIRE(b < a);
        REQUIRE_FALSE(a < b);
    }

    SECTION("a derived config sees differences in its parents and in its own options")
    {
        PrintConfig a, b;
        REQUIRE(a == b);
        REQUIRE(a.hash() == b.hash());

        b.gcode_flavor.value = b.gcode_flavor.value == gcfMarlinLegacy ? gcfKlipper : gcfMarlinLegacy;  // GCodeConfig parent
        REQUIRE(a != b);
        REQUIRE(a.hash() != b.hash());

        PrintConfig c, d;
        d.skirt_distance.value = c.skirt_distance.value + 1.0;  // PrintConfig's own list
        REQUIRE(c != d);
        REQUIRE(c.hash() != d.hash());
        REQUIRE(c.optptr("skirt_distance") == &c.skirt_distance);
        REQUIRE(c.optptr("gcode_flavor") == &c.gcode_flavor);
    }
}
