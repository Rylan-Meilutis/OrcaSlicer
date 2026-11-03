#include <catch2/catch_test_macros.hpp>
#include "libslic3r/Format/ProjectStrategy.hpp"
#include <catch2/catch_message.hpp>
#include "libslic3r/Config.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/ConnectedBodies.hpp"
#include "libslic3r/ExtrusionEntityCollection.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Surface.hpp"
#include "libslic3r/libslic3r.h"
#include "libslic3r/FilamentAssignment.hpp"
#include "libslic3r/Format/bbs_3mf.hpp"
#include "libslic3r/Format/SlicedGCode.hpp"
#include "libslic3r/GCode.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/GCode/ToolMapping.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Preset.hpp"
#include "libslic3r/PresetBundle.hpp"
#include "libslic3r/Semver.hpp"
#include "libslic3r/Utils.hpp"
#include "libslic3r/miniz_extension.hpp"
#include <miniz.h>
#include <boost/nowide/fstream.hpp>
#include "test_helpers.hpp"
#include "test_utils.hpp"
#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <ios>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace Slic3r;

TEST_CASE("Reassigning an imported uniform part changes its sliced filament", "[ToolMapping][FilamentAssignment]")
{
    auto config = Test::multifilament_config(2, {
        {"filament_map", "1,2"}, {"nozzle_diameter", "0.4,0.4"},
        {"single_extruder_multi_material", "0"}, {"enable_prime_tower", "0"},
        {"automatic_warp_mitigation", "0"}});
    Test::size_flush_to_nozzles(config);
    Model model;
    Print print;
    Test::init_print({Test::cube(2.)}, print, model, config);
    auto &object = *model.objects.front();
    object.config.set("extruder", 2);
    auto &part = object.volumes.front()->config;
    for (const char *key : {"outer_wall_filament_id", "inner_wall_filament_id", "sparse_infill_filament_id",
                           "internal_solid_filament_id", "top_surface_filament_id", "bottom_surface_filament_id"})
        part.set(key, 1);
    print.apply(model, config);
    REQUIRE(print.get_object(0)->printing_region(0).config().outer_wall_filament_id.value == 1);
    clear_uniform_feature_filament_overrides(part);
    print.apply(model, config);
    const auto &region = print.get_object(0)->printing_region(0).config();
    CHECK(region.outer_wall_filament_id.value == 2);
    CHECK(region.inner_wall_filament_id.value == 2);
    CHECK(region.sparse_infill_filament_id.value == 2);
    CHECK(region.internal_solid_filament_id.value == 2);
    CHECK(region.top_surface_filament_id.value == 2);
    CHECK(region.bottom_surface_filament_id.value == 2);
}

TEST_CASE("A sliced independent-tool job retains its motions during physical tool dispatch", "[ToolMapping][GCode]")
{
    auto config = Test::multifilament_config(2, {
        {"nozzle_diameter", "0.4,0.4"},
        {"single_extruder_multi_material", "0"},
        {"filament_map", "1,2"},
        {"filament_self_index", "1,2"},
        {"gcode_flavor", "marlin2"},
        {"support_rme_indx_dock_fan_cooldown", "1"},
        {"rme_indx_dock_fan_cooldown", "1"},
        {"layer_height", "0.2"},
        {"initial_layer_print_height", "0.2"},
        {"enable_prime_tower", "0"},
        {"skirt_loops", "0"},
        {"machine_start_gcode", "M151 W255 D0 T500\nT{initial_tool}\nG90\nM83\n"},
        {"machine_end_gcode", "M400\n"},
        {"change_filament_gcode", "T{next_extruder}\nM106 S{255 / 100 * fan_max_speed[next_extruder]}\n"},
        {"automatic_warp_mitigation", "0"}
    });
    config.set_key_value("filament_start_gcode", new ConfigOptionStrings{"", ""});
    config.set_key_value("filament_end_gcode", new ConfigOptionStrings{"", ""});
    config.set_key_value("filament_ids", new ConfigOptionStrings{"fixture", "fixture"});
    Test::size_flush_to_nozzles(config);
    Model model;
    Print print;
    const std::vector<std::vector<ConfigBase::SetDeserializeItem>> overrides{
        {{"extruder", "1"}}, {{"extruder", "2"}}};
    Test::init_print({Test::cube(2.), Test::cube(2.)}, print, model, config, &overrides);
    GCodeProcessorResult result;
    const auto original = Test::gcode(print, &result);
    REQUIRE(result.dispatch_data);
    REQUIRE(result.dispatch_data->valid);
    REQUIRE_FALSE(result.dispatch_data->templates.empty());
    const auto slicing_timestamp = print.step_state_with_timestamp(psSlicingFinished).timestamp;

    // The materials/nozzles have identical physical settings but distinct
    // preview colours, which must follow their physical assignment too.
    // No Print::process or GCode::do_export occurs below.
    auto target = print.full_print_config();
    auto &colours = target.option<ConfigOptionStrings>("filament_colour")->values;
    std::swap(colours[0], colours[1]);
    std::istringstream preflight(original);
    const auto plan = prepare_tool_mapped_gcode(preflight, *result.dispatch_data, target, {1, 0},
        "; filament_type = PLA;PLA\n");
    std::istringstream source(original);
    std::ostringstream destination;
    apply_tool_mapped_gcode(source, destination, plan);
    const auto mapped = destination.str();
    CHECK(mapped.find("\nM979 S1 C180 ") != std::string::npos);
    CHECK(mapped.find("M151 W255 D0 T500\nT1\n") != std::string::npos);
    CHECK(print.step_state_with_timestamp(psSlicingFinished).timestamp == slicing_timestamp);
    const auto motions = [](const std::string &gcode) {
        std::vector<std::string> result;
        std::istringstream stream(gcode);
        std::string line;
        while (std::getline(stream, line))
            if (line.rfind("G0 ", 0) == 0 || line.rfind("G1 ", 0) == 0 ||
                line.rfind("G2 ", 0) == 0 || line.rfind("G3 ", 0) == 0)
                result.push_back(line);
        return result;
    };
    CHECK(motions(mapped) == motions(original));
    CHECK_FALSE(motions(original).empty());

    Model mapped_model(model);
    mapped_model.objects[0]->config.set("extruder", 2);
    mapped_model.objects[1]->config.set("extruder", 1);
    Print mapped_print;
    mapped_print.set_status_silent();
    mapped_print.apply(mapped_model, target);
    REQUIRE(equivalent_tool_mapping_config(print.full_print_config(), mapped_print.full_print_config(), {1,0}));
    // Config-only fixtures can omit identifiers that a PresetBundle normally
    // supplies. Reject this before any mapped output/package file is created,
    // rather than dereferencing an empty metadata vector during packaging.
    auto missing_ids = target;
    missing_ids.set_key_value("filament_ids", new ConfigOptionStrings{});
    Print incomplete;
    incomplete.set_status_silent();
    incomplete.apply(mapped_model, missing_ids);
    REQUIRE(equivalent_tool_mapping_config(incomplete.full_print_config(), incomplete.full_print_config(), {0,1}));
    CHECK_THROWS_AS(validate_tool_mapped_sliced_gcode_3mf(incomplete, incomplete, result, {0,1}), std::invalid_argument);
    ScopedTemporaryFile raw(".gcode"), package(".gcode.3mf");
    {
        boost::nowide::ofstream stream(raw.string(), std::ios::binary);
        stream << mapped;
        REQUIRE(stream.good());
    }
    const auto original_volumes = result.print_statistics.total_volumes_per_extruder;
    const auto original_grouping = result.nozzle_group_result;
    REQUIRE(result.warnings.empty());
    REQUIRE_NOTHROW(validate_tool_mapped_sliced_gcode_3mf(print, mapped_print, result, {1,0}));
    REQUIRE(store_tool_mapped_sliced_gcode_3mf(package.string(), print, mapped_print, result, {1,0}, raw.string()));
    CHECK(result.print_statistics.total_volumes_per_extruder == original_volumes);
    CHECK(result.nozzle_group_result == original_grouping);
    CHECK(print.step_state_with_timestamp(psSlicingFinished).timestamp == slicing_timestamp);
    CHECK_FALSE(mapped_print.is_step_done(psSlicingFinished));
    mz_zip_archive archive{};
    REQUIRE(open_zip_reader(&archive, package.string()));
    const ScopeGuard close([&] { close_zip_reader(&archive); });
    const auto entry = [&](const char *name) {
        std::size_t size = 0;
        void *data = mz_zip_reader_extract_file_to_heap(&archive, name, &size, 0);
        REQUIRE(data != nullptr);
        const std::string text(static_cast<const char *>(data), size);
        mz_free(data);
        return text;
    };
    CHECK(entry("Metadata/plate_1.gcode") == mapped);
    const auto metadata = entry("Metadata/slice_info.config");
    // Source material0 is first in the statistics map, but becomes physical1.
    const auto second_slot = metadata.find("<filament id=\"2\"");
    const auto first_slot = metadata.find("<filament id=\"1\"");
    REQUIRE(second_slot != std::string::npos);
    REQUIRE(first_slot != std::string::npos);
    CHECK(second_slot < first_slot);
    CHECK(metadata.find("<nozzle id=\"1\" extruder_id=\"2\"") != std::string::npos);
    CHECK(metadata.find("<nozzle id=\"0\" extruder_id=\"1\"") != std::string::npos);
    result.warnings.push_back({});
    CHECK_THROWS_AS(validate_tool_mapped_sliced_gcode_3mf(print, mapped_print, result, {1,0}), std::invalid_argument);
    CHECK_THROWS_AS(store_tool_mapped_sliced_gcode_3mf(package.string(), print, mapped_print, result, {1,0}, raw.string()), std::invalid_argument);
}

// Opt-in local regression: imports the user's real configuration and model,
// but never slices, exports, or copies a large G-code file.
TEST_CASE("An imported project retains applied settings across physical tool dispatch", "[ToolMapping][.ToolMappingProject]")
{
    const char *path = std::getenv("ORCA_TOOL_MAPPING_PROJECT");
    if (!path) SKIP("Set ORCA_TOOL_MAPPING_PROJECT to a local three-material, eight-tool 3MF");
    const auto previous_resources = resources_dir();
    set_resources_dir(std::string(TEST_DATA_DIR) + "/../../resources");
    const ScopeGuard restore_resources([&] { set_resources_dir(previous_resources); });
    ScopedSlic3rTemporaryDir temporary_dir("tool_mapping_project");
    DynamicPrintConfig imported;
    ConfigSubstitutionContext substitutions{ForwardCompatibilitySubstitutionRule::Enable};
    PlateDataPtrs plates;
    std::vector<Preset *> presets;
    const ScopeGuard release_import([&] {
        release_PlateData_list(plates);
        for (Preset *preset : presets) delete preset;
    });
    bool is_bbl = false;
    Semver version;
    auto model = Model::read_from_file(path, &imported, &substitutions,
        LoadStrategy::LoadModel | LoadStrategy::LoadConfig | LoadStrategy::AddDefaultInstances,
        &plates, &presets, &is_bbl, &version);
    REQUIRE_FALSE(model.objects.empty());
    if (const char *reassign = std::getenv("ORCA_TOOL_MAPPING_REASSIGN"); reassign && std::string_view(reassign) == "1") {
        // Exercise explicit whole-part reassignment without modifying the user's
        // saved project. Existing mixed feature assignments remain untouched.
        for (auto *object : model.objects) {
            clear_uniform_feature_filament_overrides(object->config);
            for (auto *volume : object->volumes)
                if (volume->is_model_part())
                    clear_uniform_feature_filament_overrides(volume->config);
        }
    }
    Preset::normalize(imported);
    PresetBundle bundle;
    bundle.load_config_model(path, imported, version);
    std::vector<int> tools{1, 2, 4};
    if (const char *mapping = std::getenv("ORCA_TOOL_MAPPING_TOOLS")) {
        tools.clear();
        std::istringstream values(mapping);
        std::string value;
        while (std::getline(values, value, ',')) {
            std::size_t consumed = 0;
            const int tool = std::stoi(value, &consumed);
            REQUIRE(consumed == value.size());
            REQUIRE(tool >= 0);
            REQUIRE(tool < bundle.get_printer_extruder_count());
            tools.push_back(tool);
        }
        REQUIRE_FALSE(tools.empty());
        REQUIRE(std::set<int>(tools.begin(), tools.end()).size() == tools.size());
    } else {
        REQUIRE(bundle.filament_presets.size() == 3);
        REQUIRE(bundle.get_printer_extruder_count() == 8);
    }
    REQUIRE(bundle.filament_presets.size() == tools.size());
    DynamicPrintConfig plate;
    if (!plates.empty()) plate = plates.front()->config;
    auto source_config = bundle.full_config(false);
    source_config.apply(plate);
    Print source;
    source.set_status_silent();
    source.apply(model, source_config);
    Model mapped_model(model);
    auto mapped_config = bundle.tool_mapped_config(mapped_model, tools, plate);
    Print mapped;
    mapped.set_status_silent();
    mapped.apply(mapped_model, mapped_config);
    const auto &a = source.full_print_config();
    const auto &b = mapped.full_print_config();
    REQUIRE(equivalent_tool_mapping_config(a, b, tools));
    CHECK_FALSE(source.is_step_done(psSlicingFinished));
    CHECK_FALSE(mapped.is_step_done(psSlicingFinished));
    if (const char *slice = std::getenv("ORCA_TOOL_MAPPING_SLICE"); slice && std::string_view(slice) == "1") {
        // One logical slice supplies the in-memory certificate. Older exported
        // G-code cannot prove its macro inputs. The helper's only temporary
        // G-code file is removed before this in-memory remapping check.
        GCodeProcessorResult generated;
        const auto original = Test::gcode(source, &generated);
        if (std::getenv("ORCA_PROJECT_CONNECTIVITY")) {
            Model native_model(model);
            auto native_config = source_config;
            for (const char *key : {"interlocking_beam", "rooting", "rooting_xy"}) {
                native_config.set_deserialize_strict(key, "0");
                for (auto *object : native_model.objects) {
                    object->config.set(key, false);
                    for (auto *volume : object->volumes)
                        volume->config.set(key, false);
                }
            }
            Print native;
            native.set_status_silent();
            native.apply(native_model, native_config);
            Print raw;
            raw.set_status_silent();
            raw.apply(model, source_config);
            REQUIRE(native.objects().size() == source.objects().size());
            for (size_t object = 0; object < source.objects().size(); ++object) {
                native.get_object(object)->slice();
                raw.get_object(object)->slice();
                const auto &before = native.get_object(object)->layers();
                const auto &after = source.get_object(object)->layers();
                REQUIRE(before.size() == after.size());
                for (size_t z = 0; z < before.size(); ++z) {
                    ExPolygons native_outline, structural_outline;
                    for (const auto *region : before[z]->regions())
                        append(native_outline, to_expolygons(region->slices.surfaces));
                    for (const auto *region : raw.get_object(object)->get_layer(z)->regions())
                        append(structural_outline, to_expolygons(region->slices.surfaces));
                    CAPTURE(object, z);
                    CHECK(offset_ex(diff_ex(native_outline, structural_outline), -scaled(.001)).empty());
                    CHECK(offset_ex(diff_ex(structural_outline, native_outline), -scaled(.001)).empty());
                }
                for (bool toolpaths : {false, true}) {
                for (unsigned tool : {1u, 2u, 3u}) {
                    CAPTURE(object, tool, toolpaths);
                    std::vector<ExPolygons> deposits(after.size());
                    std::vector<const ExPolygons *> layers;
                    for (size_t z = 0; z < after.size(); ++z) {
                        Polygons coverage;
                        const auto *layer = toolpaths ? after[z] : raw.get_object(object)->get_layer(z);
                        for (const auto *region : layer->regions())
                            if (region->region().extruder(frExternalPerimeter) == tool) {
                                if (toolpaths) {
                                    region->perimeters.polygons_covered_by_width(coverage, 0.f);
                                    region->fills.polygons_covered_by_width(coverage, 0.f);
                                } else
                                    append(coverage, to_polygons(region->slices.surfaces));
                            }
                        deposits[z] = union_ex(coverage);
                        layers.push_back(&deposits[z]);
                    }
                    size_t count = 0;
                    const auto bodies = connected_bodies(layers, count);
                    std::vector<bool> anchored(count, false);
                    std::vector<bool> printable(count, false);
                    for (size_t z = 0; z < before.size(); ++z) {
                        ExPolygons original_material;
                        for (const auto *region : before[z]->regions())
                            if (region->region().extruder(frExternalPerimeter) == tool)
                                append(original_material, to_expolygons(region->slices.surfaces));
                        for (size_t i = 0; i < deposits[z].size(); ++i) {
                            // Polygon unions can leave sub-micron numerical
                            // slivers; these are not disconnected extrusions.
                            printable[bodies[z][i]] = printable[bodies[z][i]] ||
                                !offset_ex(ExPolygons{deposits[z][i]}, -scaled(.001)).empty();
                            if (!intersection_ex(ExPolygons{deposits[z][i]}, original_material).empty())
                                anchored[bodies[z][i]] = true;
                        }
                    }
                    size_t disconnected = 0;
                    for (size_t body = 0; body < count; ++body)
                        disconnected += printable[body] && !anchored[body];
                    if (disconnected) {
                        std::vector<double> volumes(count, 0.);
                        std::vector<size_t> courses(count, 0);
                        for (size_t z = 0; z < deposits.size(); ++z)
                            for (size_t i = 0; i < deposits[z].size(); ++i) {
                                const auto body = bodies[z][i];
                                volumes[body] += deposits[z][i].area() * SCALING_FACTOR * SCALING_FACTOR * after[z]->height;
                                ++courses[body];
                            }
                        for (size_t body = 0; body < count; ++body)
                            if (printable[body] && !anchored[body])
                                WARN("Disconnected object=" << object << " tool=" << tool << " volume_mm3=" << volumes[body]
                                    << " polygon_courses=" << courses[body]);
                    }
                    CHECK(disconnected == 0);
                }
                }
            }
        }
        // Preview identity stays in project filament space, independently of
        // the physical nozzle map used by dispatch.
        std::set<unsigned int> preview_filaments;
        bool preview_colors_match_filaments = true;
        for (const auto &move : generated.moves)
            if (move.type == EMoveType::Extrude) {
                preview_filaments.insert(move.extruder_id);
                preview_colors_match_filaments &= move.cp_color_id == move.extruder_id;
            }
        std::set<unsigned int> expected_filaments;
        for (unsigned int material = 0; material < tools.size(); ++material)
            expected_filaments.insert(material);
        CHECK(preview_filaments == expected_filaments);
        CHECK(preview_colors_match_filaments);
        REQUIRE(generated.dispatch_data);
        REQUIRE(generated.dispatch_data->sealed);
        REQUIRE(equivalent_tool_mapping_config(source.full_print_config(), mapped.full_print_config(), tools));
        const auto timestamp = source.step_state_with_timestamp(psSlicingFinished).timestamp;
        std::string full_config;
        GCode::append_full_config(mapped, full_config);
        GCodeDispatchPlan plan;
        {
            std::istringstream input(original);
            plan = prepare_tool_mapped_gcode(input, *generated.dispatch_data, b, tools, full_config);
        }
        const auto motion_fingerprint = [](std::string_view gcode) {
            std::uint64_t hash = 14695981039346656037ULL, bytes = 0;
            while (!gcode.empty()) {
                const auto end = gcode.find('\n');
                const auto size = end == std::string_view::npos ? gcode.size() : end + 1;
                const auto line = gcode.substr(0, size);
                if (line.size() >= 3 && line[0] == 'G' && line[1] >= '0' && line[1] <= '3' && line[2] == ' ') {
                    for (unsigned char value : line) { hash ^= value; hash *= 1099511628211ULL; }
                    bytes += line.size();
                }
                gcode.remove_prefix(size);
            }
            return std::make_pair(hash, bytes);
        };
        std::ostringstream output;
        {
            std::istringstream input(original);
            apply_tool_mapped_gcode(input, output, plan);
        }
        const auto fingerprint = motion_fingerprint(original);
        REQUIRE(fingerprint.second > 0);
        CHECK(motion_fingerprint(output.str()) == fingerprint);
        CHECK(source.step_state_with_timestamp(psSlicingFinished).timestamp == timestamp);
        CHECK_FALSE(mapped.is_step_done(psSlicingFinished));
    }
}
