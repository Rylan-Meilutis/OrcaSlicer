#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include "libslic3r/Config.hpp"
#include "libslic3r/Format/bbs_3mf.hpp"
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
#include "test_helpers.hpp"
#include "test_utils.hpp"
#include <cstdint>
#include <cstdlib>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace Slic3r;

TEST_CASE("A sliced independent-tool job retains its motions during physical tool dispatch", "[ToolMapping][GCode]")
{
    auto config = Test::multifilament_config(2, {
        {"nozzle_diameter", "0.4,0.4"},
        {"single_extruder_multi_material", "0"},
        {"gcode_flavor", "marlin2"},
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

    // The materials/nozzles are identical in this fixture; only their physical
    // assignment changes. No Print::process or GCode::do_export occurs below.
    auto target = print.full_print_config();
    std::istringstream preflight(original);
    const auto plan = prepare_tool_mapped_gcode(preflight, *result.dispatch_data, target, {1, 0},
        "; filament_type = PLA;PLA\n");
    std::istringstream source(original);
    std::ostringstream destination;
    apply_tool_mapped_gcode(source, destination, plan);
    const auto mapped = destination.str();
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
}

// Opt-in local regression: imports the user's real configuration and model,
// but never slices, exports, or copies a large G-code file.
TEST_CASE("An imported three-material project retains applied settings across eight-tool dispatch", "[ToolMapping][.ToolMappingProject]")
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
    Preset::normalize(imported);
    PresetBundle bundle;
    bundle.load_config_model(path, imported, version);
    REQUIRE(bundle.filament_presets.size() == 3);
    REQUIRE(bundle.get_printer_extruder_count() == 8);
    DynamicPrintConfig plate;
    if (!plates.empty()) plate = plates.front()->config;
    auto source_config = bundle.full_config(false);
    source_config.apply(plate);
    Print source;
    source.set_status_silent();
    source.apply(model, source_config);
    const std::vector<int> tools{1, 2, 4};
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
