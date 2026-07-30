#include "ProfileSourceManager.hpp"

#include "Http.hpp"
#include "libslic3r/AppConfig.hpp"
#include "libslic3r/Config.hpp"
#include "libslic3r/Preset.hpp"
#include "libslic3r/PresetBundle.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Utils.hpp"
#include "libslic3r_version.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <ctime>
#include <iomanip>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <unordered_map>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/nowide/fstream.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/uuid/detail/sha1.hpp>
#include <miniz/miniz.h>
#include <nlohmann/json.hpp>

namespace fs = boost::filesystem;
using json = nlohmann::json;

namespace Slic3r {
namespace {

constexpr const char *CONFIG_KEY = "profile_sources";

std::string sanitize(std::string value)
{
    for (char &ch : value)
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '-' && ch != '_')
            ch = '_';
    while (value.find("__") != std::string::npos)
        boost::replace_all(value, "__", "_");
    boost::trim_if(value, boost::is_any_of("_"));
    return value.empty() ? "source" : value;
}

std::string download(const std::string &url, std::string &error)
{
    std::string body;
    Http::get(url)
        .timeout_connect(15)
        .timeout_max(120)
        .size_limit(512 * 1024 * 1024)
        .tls_verify(true)
        .on_complete([&](std::string response, unsigned status) {
            if (status >= 200 && status < 300)
                body = std::move(response);
            else
                error = "HTTP " + std::to_string(status);
        })
        .on_error([&](std::string, std::string message, unsigned status) {
            error = message.empty() ? "HTTP " + std::to_string(status) : std::move(message);
        })
        .perform_sync();
    return body;
}

bool safe_archive_name(const std::string &name)
{
    const fs::path path(name);
    if (path.is_absolute())
        return false;
    for (const fs::path &part : path)
        if (part == "..")
            return false;
    return true;
}

bool extract_zip(const std::string &body, const fs::path &destination, std::string &error)
{
    mz_zip_archive archive {};
    if (!mz_zip_reader_init_mem(&archive, body.data(), body.size(), 0)) {
        error = "The downloaded file is not a valid ZIP archive.";
        return false;
    }
    const int count = static_cast<int>(mz_zip_reader_get_num_files(&archive));
    for (int i = 0; i < count; ++i) {
        mz_zip_archive_file_stat stat {};
        if (!mz_zip_reader_file_stat(&archive, i, &stat))
            continue;
        const std::string name = stat.m_filename;
        if (!safe_archive_name(name)) {
            error = "The archive contains an unsafe path.";
            mz_zip_reader_end(&archive);
            return false;
        }
        const fs::path output = destination / fs::path(name);
        if (mz_zip_reader_is_file_a_directory(&archive, i)) {
            fs::create_directories(output);
            continue;
        }
        fs::create_directories(output.parent_path());
        const std::string encoded_output = encode_path(output.string().c_str());
        if (!mz_zip_reader_extract_to_file(&archive, i, encoded_output.c_str(), 0)) {
            error = "Could not extract " + name;
            mz_zip_reader_end(&archive);
            return false;
        }
    }
    mz_zip_reader_end(&archive);
    return true;
}

std::string archive_url(std::string url)
{
    boost::trim_right_if(url, boost::is_any_of("/"));
    constexpr const char *github_prefix = "https://github.com/";
    if (boost::starts_with(url, github_prefix) && !boost::ends_with(url, ".zip")) {
        std::string repository = url.substr(std::char_traits<char>::length(github_prefix));
        boost::trim_right_if(repository, boost::is_any_of("/"));
        if (boost::ends_with(repository, ".git"))
            repository.resize(repository.size() - 4);
        // GitHub's zipball endpoint follows the repository's configured
        // default branch. Forks do not consistently use main or master.
        if (std::count(repository.begin(), repository.end(), '/') == 1)
            return "https://api.github.com/repos/" + repository + "/zipball";
    }
    return url;
}

using Section = std::map<std::string, std::string>;
using Sections = std::map<std::string, Section>;

Sections read_sections(const fs::path &path, std::string &error)
{
    boost::property_tree::ptree tree;
    try {
        boost::nowide::ifstream stream(path.string());
        boost::property_tree::read_ini(stream, tree);
    } catch (const std::exception &e) {
        error = e.what();
        return {};
    }
    Sections result;
    for (const auto &section : tree)
        for (const auto &item : section.second)
            result[section.first][item.first] = item.second.data();
    return result;
}

Section read_flat_ini(const fs::path &path)
{
    boost::nowide::ifstream stream(path.string());
    Section result;
    std::string line;
    bool saw_section = false;
    while (std::getline(stream, line)) {
        boost::trim(line);
        if (line.empty() || line.front() == '#' || line.front() == ';')
            continue;
        if (line.front() == '[') {
            saw_section = true;
            continue;
        }
        // A bundle is handled by read_sections(). Do not accidentally combine
        // all of its sections into one profile if its profile types are unsupported.
        if (saw_section)
            continue;
        const size_t separator = line.find('=');
        if (separator == std::string::npos)
            continue;
        std::string key = line.substr(0, separator);
        std::string value = line.substr(separator + 1);
        boost::trim(key);
        boost::trim(value);
        if (!key.empty())
            result[key] = value;
    }
    return result;
}

std::string detect_flat_type(const fs::path &path, const Section &input)
{
    if (auto explicit_type = input.find("ini_type"); explicit_type != input.end()) {
        const std::string value = boost::to_lower_copy(explicit_type->second);
        if (value == "print") return "print";
        if (value == "filament") return "filament";
        if (value == "printer" || value == "machine") return "printer";
    }
    for (fs::path parent = path.parent_path(); !parent.empty();) {
        const std::string directory = boost::to_lower_copy(parent.filename().string());
        if (directory == "print" || directory == "process") return "print";
        if (directory == "filament") return "filament";
        if (directory == "printer" || directory == "machine") return "printer";
        const fs::path next = parent.parent_path();
        if (next == parent)
            break;
        parent = next;
    }
    const auto count = [&](std::initializer_list<const char *> keys) {
        return std::count_if(keys.begin(), keys.end(), [&](const char *key) { return input.count(key) != 0; });
    };
    const size_t printer = count({"nozzle_diameter", "bed_shape", "gcode_flavor", "printer_technology",
                                  "max_print_height", "start_gcode"});
    const size_t filament = count({"filament_type", "filament_diameter", "temperature",
                                   "filament_max_volumetric_speed", "filament_colour"});
    const size_t print = count({"layer_height", "perimeters", "fill_density", "fill_pattern",
                                "perimeter_speed", "infill_speed"});
    const size_t best = std::max({printer, filament, print});
    if (best < 2)
        return {};
    return best == printer ? "printer" : best == filament ? "filament" : "print";
}

Section flatten(const std::string &name, const Sections &sections, std::set<std::string> &visiting,
                std::map<std::string, Section> &cache, size_t &skipped)
{
    if (auto found = cache.find(name); found != cache.end())
        return found->second;
    if (!visiting.insert(name).second) {
        ++skipped;
        return {};
    }
    Section result;
    auto found = sections.find(name);
    if (found == sections.end()) {
        visiting.erase(name);
        ++skipped;
        return result;
    }
    auto inherited = found->second.find("inherits");
    if (inherited != found->second.end()) {
        std::vector<std::string> parents;
        unescape_strings_cstyle(inherited->second, parents);
        for (const std::string &parent : parents) {
            const std::string prefix = name.substr(0, name.find(':') + 1);
            Section base = flatten(prefix + parent, sections, visiting, cache, skipped);
            result.insert(base.begin(), base.end());
        }
    }
    for (const auto &[key, value] : found->second)
        if (key != "inherits")
            result[key] = value;
    visiting.erase(name);
    return cache[name] = std::move(result);
}

const std::unordered_map<std::string, std::string> &key_map()
{
    // Translation coverage follows the GPLv3 SuperSlicer/PrusaSlicer to
    // OrcaSlicer converter by theophile, adapted here to Orca's live option
    // definitions so obsolete source fields are rejected safely.
    static const std::unordered_map<std::string, std::string> map {
        {"arc_fitting", "enable_arc_fitting"},
        {"perimeters", "wall_loops"},
        {"top_solid_layers", "top_shell_layers"},
        {"bottom_solid_layers", "bottom_shell_layers"},
        {"top_solid_min_thickness", "top_shell_thickness"},
        {"bottom_solid_min_thickness", "bottom_shell_thickness"},
        {"fill_density", "sparse_infill_density"},
        {"fill_pattern", "sparse_infill_pattern"},
        {"solid_fill_pattern", "internal_solid_infill_pattern"},
        {"top_fill_pattern", "top_surface_pattern"},
        {"bottom_fill_pattern", "bottom_surface_pattern"},
        {"perimeter_speed", "inner_wall_speed"},
        {"external_perimeter_speed", "outer_wall_speed"},
        {"small_perimeter_speed", "small_perimeter_speed"},
        {"infill_speed", "sparse_infill_speed"},
        {"solid_infill_speed", "internal_solid_infill_speed"},
        {"top_solid_infill_speed", "top_surface_speed"},
        {"gap_fill_speed", "gap_infill_speed"},
        {"support_material_speed", "support_speed"},
        {"support_material_interface_speed", "support_interface_speed"},
        {"bridge_speed_internal", "internal_bridge_speed"},
        {"first_layer_speed", "initial_layer_speed"},
        {"first_layer_infill_speed", "initial_layer_infill_speed"},
        {"extrusion_width", "line_width"},
        {"perimeter_extrusion_width", "inner_wall_line_width"},
        {"external_perimeter_extrusion_width", "outer_wall_line_width"},
        {"infill_extrusion_width", "sparse_infill_line_width"},
        {"solid_infill_extrusion_width", "internal_solid_infill_line_width"},
        {"top_infill_extrusion_width", "top_surface_line_width"},
        {"first_layer_extrusion_width", "initial_layer_line_width"},
        {"support_material_extrusion_width", "support_line_width"},
        {"skirts", "skirt_loops"},
        {"skirt_distance", "skirt_distance"},
        {"skirt_height", "skirt_height"},
        {"brim_separation", "brim_object_gap"},
        {"wipe_tower", "enable_prime_tower"},
        {"wipe_tower_width", "prime_tower_width"},
        {"wipe_tower_brim_width", "prime_tower_brim_width"},
        {"wipe_tower_no_sparse_layers", "wipe_tower_no_sparse_layers"},
        {"avoid_crossing_perimeters", "reduce_crossing_wall"},
        {"avoid_crossing_perimeters_max_detour", "max_travel_detour_distance"},
        {"only_retract_when_crossing_perimeters", "reduce_infill_retraction"},
        {"external_perimeters_first", "external_perimeters_first"},
        {"infill_first", "infill_first"},
        {"first_layer_height", "initial_layer_print_height"},
        {"first_layer_acceleration", "initial_layer_acceleration"},
        {"solid_infill_acceleration", "internal_solid_infill_acceleration"},
        {"top_solid_infill_acceleration", "top_surface_acceleration"},
        {"infill_acceleration", "sparse_infill_acceleration"},
        {"perimeter_acceleration", "inner_wall_acceleration"},
        {"external_perimeter_acceleration", "outer_wall_acceleration"},
        {"fill_angle", "infill_direction"},
        {"infill_overlap", "infill_wall_overlap"},
        {"solid_infill_below_area", "minimum_sparse_infill_area"},
        {"infill_every_layers", "infill_combination"},
        {"first_layer_size_compensation", "elefant_foot_compensation"},
        {"first_layer_size_compensation_layers", "elefant_foot_compensation_layers"},
        {"xy_size_compensation", "xy_contour_compensation"},
        {"xy_inner_size_compensation", "xy_hole_compensation"},
        {"output_filename_format", "filename_format"},
        {"complete_objects", "print_sequence"},
        {"perimeter_generator", "wall_generator"},
        {"overhangs", "detect_overhang_wall"},
        {"thin_walls", "detect_thin_wall"},
        {"extra_perimeters_odd_layers", "alternate_extra_wall"},
        {"extra_perimeters_on_overhangs", "extra_perimeters_on_overhangs"},
        {"extra_perimeters_overhangs", "extra_perimeters_on_overhangs"},
        {"only_one_perimeter_first_layer", "only_one_wall_first_layer"},
        {"only_one_perimeter_top", "only_one_wall_top"},
        {"small_perimeter_min_length", "small_perimeter_threshold"},
        {"slice_closing_radius", "slice_closing_radius"},
        {"slicing_mode", "slicing_mode"},
        {"resolution", "resolution"},
        {"seam_position", "seam_position"},
        {"seam_gap", "seam_gap"},
        {"staggered_inner_seams", "staggered_inner_seams"},
        {"support_material", "enable_support"},
        {"support_material_angle", "support_angle"},
        {"support_material_enforce_layers", "enforce_support_layers"},
        {"support_material_spacing", "support_base_pattern_spacing"},
        {"support_material_contact_distance", "support_top_z_distance"},
        {"support_material_bottom_contact_distance", "support_bottom_z_distance"},
        {"support_material_bottom_interface_layers", "support_interface_bottom_layers"},
        {"support_material_interface_layers", "support_interface_top_layers"},
        {"support_material_interface_spacing", "support_interface_spacing"},
        {"support_material_interface_contact_loops", "support_interface_loop_pattern"},
        {"support_material_buildplate_only", "support_on_build_plate_only"},
        {"support_material_threshold", "support_threshold_angle"},
        {"support_material_pattern", "support_base_pattern"},
        {"support_material_interface_pattern", "support_interface_pattern"},
        {"support_material_xy_spacing", "support_object_xy_distance"},
        {"support_material_layer_height", "independent_support_layer_height"},
        {"support_material_extruder", "support_filament"},
        {"support_material_interface_extruder", "support_interface_filament"},
        {"infill_extruder", "sparse_infill_filament"},
        {"solid_infill_extruder", "solid_infill_filament"},
        {"perimeter_extruder", "wall_filament"},
        {"first_layer_extruder", "first_layer_filament"},
        {"bridge_flow_ratio", "bridge_flow"},
        {"fill_top_flow_ratio", "top_solid_infill_flow_ratio"},
        {"first_layer_flow_ratio", "bottom_solid_infill_flow_ratio"},
        {"bridge_overlap_min", "bridge_density"},
        {"dont_support_bridges", "bridge_no_support"},
        {"machine_limits_usage", "emit_machine_limits_to_gcode"},
        {"remaining_times", "disable_m73"},
        {"draft_shield", "draft_shield"},
        {"spiral_vase", "spiral_mode"},
        {"fuzzy_skin_point_dist", "fuzzy_skin_point_distance"},
        {"fuzzy_skin_thickness", "fuzzy_skin_thickness"},
        {"init_z_rotate", "preferred_orientation"},
        {"retract_length", "retraction_length"},
        {"retract_speed", "retraction_speed"},
        {"deretract_speed", "deretraction_speed"},
        {"retract_lift", "z_hop"},
        {"retract_before_travel", "retraction_minimum_travel"},
        {"retract_layer_change", "retract_when_changing_layer"},
        {"retract_restart_extra", "retract_restart_extra"},
        {"retract_length_toolchange", "retract_length_toolchange"},
        {"retract_restart_extra_toolchange", "retract_restart_extra_toolchange"},
        {"retract_before_wipe", "retract_before_wipe"},
        {"wipe", "wipe"},
        {"retract_lift_top", "retract_lift_enforce"},
        {"before_layer_gcode", "before_layer_change_gcode"},
        {"toolchange_gcode", "change_filament_gcode"},
        {"layer_gcode", "layer_change_gcode"},
        {"feature_gcode", "change_extrusion_role_gcode"},
        {"end_gcode", "machine_end_gcode"},
        {"start_gcode", "machine_start_gcode"},
        {"pause_print_gcode", "machine_pause_gcode"},
        {"bed_shape", "printable_area"},
        {"max_print_height", "printable_height"},
        {"machine_max_feedrate_x", "machine_max_speed_x"},
        {"machine_max_feedrate_y", "machine_max_speed_y"},
        {"machine_max_feedrate_z", "machine_max_speed_z"},
        {"machine_max_feedrate_e", "machine_max_speed_e"},
        {"temperature", "nozzle_temperature"},
        {"first_layer_temperature", "nozzle_temperature_initial_layer"},
        {"bed_temperature", "hot_plate_temp"},
        {"first_layer_bed_temperature", "hot_plate_temp_initial_layer"},
        {"extrusion_multiplier", "filament_flow_ratio"},
        {"filament_colour", "default_filament_colour"},
        {"filament_deretract_speed", "filament_deretraction_speed"},
        {"filament_retract_before_travel", "filament_retraction_minimum_travel"},
        {"filament_retract_length", "filament_retraction_length"},
        {"filament_retract_lift", "filament_z_hop"},
        {"filament_retract_speed", "filament_retraction_speed"},
        {"filament_retract_layer_change", "filament_retract_when_changing_layer"},
        {"max_fan_speed", "fan_max_speed"},
        {"min_fan_speed", "fan_min_speed"},
        {"min_print_speed", "slow_down_min_speed"},
        {"slowdown_below_layer_time", "slow_down_layer_time"},
        {"disable_fan_first_layers", "close_fan_the_first_x_layers"},
        {"fan_below_layer_time", "fan_cooling_layer_time"},
        {"bridge_fan_speed", "overhang_fan_speed"},
        {"external_perimeter_fan_speed", "overhang_fan_threshold"},
        {"fan_always_on", "reduce_fan_stop_start_freq"},
        {"start_filament_gcode", "filament_start_gcode"},
        {"end_filament_gcode", "filament_end_gcode"},
        {"filament_max_volumetric_speed", "filament_max_volumetric_speed"}
    };
    return map;
}

bool is_percent(const std::string &value)
{
    return !value.empty() && value.back() == '%';
}

std::optional<double> number(const std::string &value)
{
    try {
        size_t consumed = 0;
        const double parsed = std::stod(value, &consumed);
        if (consumed == value.size() || (consumed + 1 == value.size() && value.back() == '%'))
            return parsed;
    } catch (...) {
    }
    return std::nullopt;
}

std::string scalar(double value)
{
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(4) << value;
    std::string result = stream.str();
    while (result.size() > 1 && result.back() == '0')
        result.pop_back();
    if (!result.empty() && result.back() == '.')
        result.pop_back();
    return result;
}

std::string percent_to_absolute(const std::string &value, const std::string &base)
{
    if (!is_percent(value))
        return value;
    const auto percentage = number(value);
    const auto base_value = number(base);
    return percentage && base_value && !is_percent(base) ? scalar(*base_value * *percentage / 100.0) : value;
}

std::string converted_value(const std::string &source_key, const std::string &value,
                            const Section &input, const std::string &type)
{
    static const std::unordered_map<std::string, std::string> filament_types {
        {"PET", "PETG"}, {"FLEX", "TPU"}, {"NYLON", "PA"}
    };
    static const std::unordered_map<std::string, std::string> patterns {
        {"rectilinear", "zig-zag"}, {"rectilineargapfill", "zig-zag"},
        {"rectiwithperimeter", "zig-zag"}, {"sawtooth", "zig-zag"},
        {"scatteredrectilinear", "zig-zag"}, {"monotoniclines", "monotonicline"},
        {"smooth", "monotonic"}, {"smoothhilbert", "hilbertcurve"},
        {"smoothtriple", "triangles"}, {"stars", "tri-hexagon"},
        {"concentricgapfill", "concentric"}, {"monotonicgapfill", "monotonic"}
    };
    static const std::unordered_map<std::string, std::string> seams {
        {"cost", "nearest"}, {"random", "random"}, {"allrandom", "random"},
        {"aligned", "aligned"}, {"contiguous", "aligned"}, {"rear", "back"}, {"nearest", "nearest"}
    };
    static const std::unordered_map<std::string, std::string> flavors {
        {"klipper", "klipper"}, {"marlin", "marlin"}, {"marlin2", "marlin2"},
        {"reprap", "reprapfirmware"}, {"reprapfirmware", "reprapfirmware"},
        {"repetier", "reprapfirmware"}, {"smoothie", "reprapfirmware"},
        {"mach3", "reprapfirmware"}, {"machinekit", "reprapfirmware"},
        {"makerware", "reprapfirmware"}, {"sailfish", "reprapfirmware"},
        {"teacup", "reprapfirmware"}, {"sprinter", "reprapfirmware"},
        {"no-extrusion", "reprapfirmware"}
    };
    if (source_key == "filament_type") {
        auto found = filament_types.find(value);
        return found == filament_types.end() ? value : found->second;
    }
    if (source_key == "fill_pattern" || source_key == "solid_fill_pattern" ||
        source_key == "top_fill_pattern" || source_key == "bottom_fill_pattern") {
        auto found = patterns.find(value);
        return found == patterns.end() ? value : found->second;
    }
    if (source_key == "seam_position") {
        auto found = seams.find(value);
        return found == seams.end() ? value : found->second;
    }
    if (source_key == "gcode_flavor") {
        auto found = flavors.find(value);
        return found == flavors.end() ? value : found->second;
    }
    if (source_key == "complete_objects")
        return value == "0" ? "by layer" : "by object";
    if (source_key == "machine_limits_usage")
        return value == "emit_to_gcode" ? "1" : "0";
    if (source_key == "remaining_times")
        return value == "0" ? "1" : "0";
    if (source_key == "draft_shield")
        return value == "disabled" ? "0" : value == "enabled" ? "1" : value;
    if (source_key == "support_material_layer_height")
        return value == "0" ? "0" : "1";
    if (source_key == "infill_every_layers")
        return value == "0" ? "0" : "1";
    if (source_key == "output_filename_format") {
        std::string result = value;
        std::replace(result.begin(), result.end(), '[', '{');
        std::replace(result.begin(), result.end(), ']', '}');
        return result;
    }
    if (source_key == "bridge_flow_ratio" || source_key == "fill_top_flow_ratio" ||
        source_key == "first_layer_flow_ratio") {
        const auto parsed = number(value);
        return is_percent(value) && parsed ? scalar(std::min(2.0, *parsed / 100.0)) : value;
    }
    if (source_key == "external_perimeter_fan_speed") {
        const auto parsed = number(value);
        return parsed && *parsed >= 0 ? scalar(*parsed) + "%" : "0%";
    }
    if (source_key == "external_perimeter_speed" || source_key == "first_layer_speed") {
        auto base = input.find("perimeter_speed");
        return base == input.end() ? value : percent_to_absolute(value, base->second);
    }
    if (source_key == "solid_infill_speed") {
        auto base = input.find("infill_speed");
        return base == input.end() ? value : percent_to_absolute(value, base->second);
    }
    if (source_key == "first_layer_infill_speed") {
        auto base = input.find("infill_speed");
        return base == input.end() ? value : percent_to_absolute(value, base->second);
    }
    if (source_key == "support_material_interface_speed") {
        auto base = input.find("support_material_speed");
        return base == input.end() ? value : percent_to_absolute(value, base->second);
    }
    if (source_key == "support_material_bottom_contact_distance" && value == "0") {
        auto top = input.find("support_material_contact_distance");
        return top == input.end() ? value : top->second;
    }
    if (source_key == "filament_max_volumetric_speed" && number(value).value_or(0.0) <= 0.0) {
        static const std::unordered_map<std::string, std::string> defaults {
            {"PLA", "15"}, {"PET", "10"}, {"PETG", "10"}, {"ABS", "12"}, {"ASA", "12"},
            {"FLEX", "3.2"}, {"TPU", "3.2"}, {"NYLON", "12"}, {"PA", "12"},
            {"PVA", "12"}, {"PC", "12"}, {"HIPS", "8"}, {"PP", "8"}
        };
        auto material = input.find("filament_type");
        auto fallback = material == input.end() ? defaults.end() : defaults.find(material->second);
        return fallback == defaults.end() ? "8" : fallback->second;
    }
    // Prusa print flow is profile-wide, while Orca's similarly named option
    // under filament is the only compatible destination.
    if (source_key == "extrusion_multiplier" && type != "filament")
        return value;
    return value == "nil" ? std::string() : value;
}

bool is_vector_option(const std::string &key, const std::string &preset_type)
{
    const ConfigOptionDef *definition = print_config_def.get(key);
    if (definition == nullptr || (definition->type & coVectorType) == 0)
        return false;
    if (preset_type == "filament")
        return true;
    static const std::set<std::string> profile_vectors {
        "compatible_printers", "compatible_prints", "post_process",
        "printable_area", "thumbnails", "default_filament_profile",
        "machine_max_acceleration_e", "machine_max_acceleration_extruding",
        "machine_max_acceleration_retracting", "machine_max_acceleration_travel",
        "machine_max_acceleration_x", "machine_max_acceleration_y", "machine_max_acceleration_z",
        "machine_max_speed_e", "machine_max_speed_x", "machine_max_speed_y", "machine_max_speed_z",
        "machine_max_jerk_e", "machine_max_jerk_x", "machine_max_jerk_y", "machine_max_jerk_z",
        "machine_min_extruding_rate", "machine_min_travel_rate"
    };
    return profile_vectors.count(key) != 0;
}

void add_value(json &output, const std::string &key, const std::string &value,
               const std::string &preset_type)
{
    if (is_vector_option(key, preset_type)) {
        std::vector<std::string> values;
        unescape_strings_cstyle(value, values);
        if (values.size() <= 1 && value.find(';') != std::string::npos)
            boost::split(values, value, boost::is_any_of(";"), boost::token_compress_off);
        const bool literal_string = key.find("gcode") != std::string::npos ||
                                    key.find("notes") != std::string::npos ||
                                    key == "post_process";
        if (values.size() <= 1 && !literal_string && value.find(',') != std::string::npos)
            boost::split(values, value, boost::is_any_of(","), boost::token_compress_off);
        if (values.empty())
            values.push_back(value);
        output[key] = values;
    } else {
        output[key] = value;
    }
}

std::string preset_type(const std::string &section)
{
    if (boost::starts_with(section, "print:"))
        return "process";
    if (boost::starts_with(section, "filament:"))
        return "filament";
    if (boost::starts_with(section, "printer:"))
        return "machine";
    return {};
}

std::string preset_name(const std::string &section)
{
    const size_t colon = section.find(':');
    return colon == std::string::npos ? section : section.substr(colon + 1);
}

bool write_preset(const fs::path &root, const std::string &section_name, const Section &input,
                  ProfileSourceSyncResult &result)
{
    const std::string type = preset_type(section_name);
    const std::string name = preset_name(section_name);
    if (type.empty() || name.empty() || name.front() == '*')
        return true;

    json output {
        {"type", type},
        {"name", name},
        {"from", "bundle"},
        {"version", SoftFever_VERSION},
        {"instantiation", "true"}
    };
    for (const auto &[source_key, value] : input) {
        if (source_key == "inherits" || source_key == "alias" || source_key == "renamed_from")
            continue;
        if (value == "nil")
            continue;
        auto mapped = key_map().find(source_key);
        std::string key = mapped == key_map().end() ? source_key : mapped->second;
        if (source_key == "extrusion_multiplier")
            key = type == "filament" ? "filament_flow_ratio" : "print_flow_ratio";
        if (!print_config_def.has(key)) {
            ++result.skipped_options;
            continue;
        }
        const std::string converted = converted_value(source_key, value, input, type);
        if (!converted.empty())
            add_value(output, key, converted, type);
        if (source_key == "bed_temperature" || source_key == "first_layer_bed_temperature") {
            const bool initial = source_key == "first_layer_bed_temperature";
            for (const char *plate : {"cool_plate_temp", "eng_plate_temp", "hot_plate_temp", "textured_plate_temp"}) {
                const std::string plate_key = std::string(plate) + (initial ? "_initial_layer" : "");
                if (print_config_def.has(plate_key))
                    add_value(output, plate_key, converted, type);
            }
        }
    }
    if (type == "process") {
        const bool outer_first = input.count("external_perimeters_first") &&
                                 input.at("external_perimeters_first") != "0";
        const bool infill_first = input.count("infill_first") && input.at("infill_first") != "0";
        output["wall_infill_order"] = infill_first ?
            (outer_first ? "infill/outer wall/inner wall" : "infill/inner wall/outer wall") :
            (outer_first ? "outer wall/inner wall/infill" : "inner wall/outer wall/infill");

        const bool ironing = input.count("ironing") && input.at("ironing") != "0";
        output["ironing_type"] = ironing && input.count("ironing_type") ?
            input.at("ironing_type") : "no ironing";

        const bool dynamic_overhang = input.count("enable_dynamic_overhang_speeds") &&
                                      input.at("enable_dynamic_overhang_speeds") != "0";
        output["enable_overhang_speed"] = dynamic_overhang ? "1" : "0";
        if (dynamic_overhang) {
            const std::array<const char *, 4> source {
                "overhang_speed_3", "overhang_speed_2", "overhang_speed_1", "overhang_speed_0"
            };
            const std::array<const char *, 4> destination {
                "overhang_1_4_speed", "overhang_2_4_speed", "overhang_3_4_speed", "overhang_4_4_speed"
            };
            for (size_t i = 0; i < source.size(); ++i)
                if (auto found = input.find(source[i]); found != input.end())
                    output[destination[i]] = found->second;
        }
    }
    if (type == "filament")
        output["filament_id"] = ProfileSourceManager::make_id(name, section_name);

    const char *subdir = type == "process" ? PRESET_PRINT_NAME :
                         type == "filament" ? PRESET_FILAMENT_NAME : PRESET_PRINTER_NAME;
    fs::create_directories(root / subdir);
    fs::path file = root / subdir / (sanitize(name) + ".json");
    boost::nowide::ofstream stream(file.string());
    stream << output.dump(1, '\t') << '\n';
    if (!stream.good())
        return false;
    if (type == "process") ++result.processes;
    else if (type == "filament") ++result.filaments;
    else ++result.printers;
    return true;
}

std::vector<fs::path> selected_ini_files(const fs::path &root)
{
    std::map<fs::path, std::vector<fs::path>> by_directory;
    if (!fs::exists(root))
        return {};
    for (fs::recursive_directory_iterator it(root), end; it != end; ++it) {
        if (!fs::is_regular_file(it->path()) || it->path().extension() != ".ini")
            continue;
        by_directory[it->path().parent_path()].push_back(it->path());
    }
    std::vector<fs::path> result;
    for (const auto &[directory, files] : by_directory) {
        fs::path newest;
        boost::optional<Semver> newest_version;
        for (const fs::path &file : files) {
            const auto version = Semver::parse(file.stem().string());
            if (version && (!newest_version || *newest_version < *version)) {
                newest = file;
                newest_version = version;
            }
        }
        if (newest_version)
            result.push_back(newest);
        else
            result.insert(result.end(), files.begin(), files.end());
    }
    return result;
}

ProfileSourceSyncResult convert_prusa_tree(const fs::path &input_root, const fs::path &output_root)
{
    ProfileSourceSyncResult result;
    const auto files = selected_ini_files(input_root);
    if (files.empty()) {
        result.error = "No PrusaSlicer INI profiles were found in the source.";
        return result;
    }
    for (const fs::path &file : files) {
        std::string parse_error;
        Sections sections = read_sections(file, parse_error);
        if (!parse_error.empty())
            sections.clear();
        const bool has_profiles = std::any_of(sections.begin(), sections.end(),
            [](const auto &entry) { return !preset_type(entry.first).empty(); });
        if (!has_profiles) {
            const Section flat = read_flat_ini(file);
            const std::string flat_type = detect_flat_type(file, flat);
            if (!flat_type.empty()) {
                const std::string name = flat.count("profile_name") ?
                    flat.at("profile_name") : file.stem().string();
                if (!write_preset(output_root, flat_type + ":" + name, flat, result)) {
                    result.error = "Could not write converted preset " + name;
                    return result;
                }
            }
            continue;
        }
        std::set<std::string> visiting;
        std::map<std::string, Section> cache;
        for (const auto &[name, section] : sections) {
            if (preset_type(name).empty())
                continue;
            Section merged = flatten(name, sections, visiting, cache, result.skipped_options);
            if (!write_preset(output_root, name, merged, result)) {
                result.error = "Could not write converted preset " + preset_name(name);
                return result;
            }
        }
    }
    if (result.printers + result.filaments + result.processes == 0)
        result.error = "No compatible printer, filament, or process profiles were found in the source.";
    return result;
}

} // namespace

ProfileSourceManager::ProfileSourceManager(AppConfig &config) : m_config(config) {}

ProfileSourceSyncResult ProfileSourceManager::convert_prusa_profiles(const std::string &input_root,
                                                                      const std::string &output_root)
{
    return convert_prusa_tree(fs::path(input_root), fs::path(output_root));
}

std::vector<ProfileSource> ProfileSourceManager::sources() const
{
    const std::string serialized = m_config.get(CONFIG_KEY);
    if (serialized.empty())
        return default_sources();
    std::vector<ProfileSource> result;
    try {
        for (const json &entry : json::parse(serialized)) {
            ProfileSource source;
            source.id        = entry.value("id", "");
            source.name      = entry.value("name", "");
            source.url       = entry.value("url", "");
            source.format    = entry.value("format", "prusa") == "orca" ? ProfileSource::Format::Orca : ProfileSource::Format::Prusa;
            source.last_sync = entry.value("last_sync", 0LL);
            source.enabled   = entry.value("enabled", false);
            if (!source.id.empty() && !source.url.empty())
                result.push_back(std::move(source));
        }
    } catch (...) {
        return default_sources();
    }
    return result;
}

void ProfileSourceManager::set_sources(const std::vector<ProfileSource> &sources)
{
    json output = json::array();
    for (const ProfileSource &source : sources)
        output.push_back({{"id", source.id}, {"name", source.name}, {"url", source.url},
                          {"format", source.format == ProfileSource::Format::Orca ? "orca" : "prusa"},
                          {"last_sync", source.last_sync}, {"enabled", source.enabled}});
    m_config.set(CONFIG_KEY, output.dump());
    m_config.save();
}

void ProfileSourceManager::add(const ProfileSource &source)
{
    auto current = sources();
    auto found = std::find_if(current.begin(), current.end(), [&](const ProfileSource &item) { return item.id == source.id; });
    if (found == current.end())
        current.push_back(source);
    else
        *found = source;
    set_sources(current);
}

bool ProfileSourceManager::remove(const std::string &id, std::string &error)
{
    auto current = sources();
    current.erase(std::remove_if(current.begin(), current.end(), [&](const ProfileSource &source) { return source.id == id; }), current.end());
    set_sources(current);
    boost::system::error_code ec;
    const std::string user = m_config.get("preset_folder").empty() ? DEFAULT_USER_FOLDER_NAME : m_config.get("preset_folder");
    fs::remove_all(fs::path(data_dir()) / PRESET_USER_DIR / user / PRESET_LOCAL_DIR / ("profile_source_" + sanitize(id)), ec);
    if (ec)
        error = ec.message();
    return !ec;
}

ProfileSourceSyncResult ProfileSourceManager::sync(const ProfileSource &source)
{
    ProfileSourceSyncResult result;
    const std::string user = m_config.get("preset_folder").empty() ? DEFAULT_USER_FOLDER_NAME : m_config.get("preset_folder");
    const fs::path target = fs::path(data_dir()) / PRESET_USER_DIR / user / PRESET_LOCAL_DIR / ("profile_source_" + sanitize(source.id));
    const fs::path staging = target.string() + ".staging";
    boost::system::error_code ec;
    fs::remove_all(staging, ec);
    fs::create_directories(staging, ec);
    if (ec) {
        result.error = ec.message();
        return result;
    }

    std::string error;
    std::string body = download(archive_url(source.url), error);
    if (!error.empty() && boost::starts_with(source.url, "https://github.com/") &&
        !boost::ends_with(source.url, ".zip")) {
        error.clear();
        std::string head_url = source.url;
        boost::trim_right_if(head_url, boost::is_any_of("/"));
        body = download(head_url + "/archive/HEAD.zip", error);
    }
    if (!error.empty()) {
        result.error = error;
        return result;
    }

    fs::path extracted = staging / "download";
    fs::create_directories(extracted);
    if (source.format == ProfileSource::Format::Prusa) {
        if (boost::ends_with(boost::to_lower_copy(source.url), ".ini")) {
            boost::nowide::ofstream stream((extracted / "source.ini").string(), std::ios::binary);
            stream.write(body.data(), body.size());
        } else if (!extract_zip(body, extracted, result.error)) {
            return result;
        }
        result = convert_prusa_tree(extracted, staging);
        if (!result.success())
            return result;
    } else {
        if (!extract_zip(body, extracted, result.error))
            return result;
        for (fs::recursive_directory_iterator it(extracted), end; it != end; ++it) {
            if (!fs::is_regular_file(it->path()) || it->path().extension() != ".json")
                continue;
            const std::string parent = it->path().parent_path().filename().string();
            if (parent != PRESET_PRINT_NAME && parent != PRESET_FILAMENT_NAME && parent != PRESET_PRINTER_NAME)
                continue;
            fs::create_directories(staging / parent);
            // Native repositories contain one machine/process/filament
            // directory per vendor, and common files often share names.
            // Local bundles are flat, so retain the vendor in the destination
            // filename to prevent one fork or vendor from silently replacing
            // another. Preset inheritance uses the JSON `name`, not filenames.
            const std::string vendor = sanitize(it->path().parent_path().parent_path().filename().string());
            const fs::path destination = staging / parent /
                (vendor + "__" + it->path().filename().string());
            fs::copy_file(it->path(), destination, fs::copy_option::overwrite_if_exists);
            if (parent == PRESET_PRINT_NAME) ++result.processes;
            else if (parent == PRESET_FILAMENT_NAME) ++result.filaments;
            else ++result.printers;
        }
    }

    if (result.printers + result.filaments + result.processes == 0) {
        result.error = "No compatible printer, filament, or process profiles were found in the source.";
        return result;
    }

    BundleMetadata metadata;
    metadata.id            = "profile_source_" + sanitize(source.id);
    metadata.name          = source.name;
    metadata.version       = "1.0.0";
    metadata.description   = "Profiles synchronized from " + source.url;
    metadata.author        = source.format == ProfileSource::Format::Prusa ? "PrusaSlicer profile source" : "OrcaSlicer profile source";
    metadata.imported_time = metadata.updated_time = std::time(nullptr);
    if (!metadata.save_to_json((staging / PRESET_BUNDLE_METADATA).string())) {
        result.error = "Could not save source metadata.";
        return result;
    }

    const fs::path backup = target.string() + ".backup";
    fs::remove_all(backup, ec);
    ec.clear();
    if (fs::exists(target)) {
        fs::rename(target, backup, ec);
        if (ec) {
            result.error = ec.message();
            return result;
        }
    }
    fs::rename(staging, target, ec);
    if (ec) {
        result.error = ec.message();
        boost::system::error_code rollback_error;
        if (fs::exists(backup))
            fs::rename(backup, target, rollback_error);
        return result;
    }
    fs::remove_all(backup, ec);
    auto current = sources();
    for (ProfileSource &item : current)
        if (item.id == source.id)
            item.last_sync = std::time(nullptr);
    set_sources(current);
    return result;
}

std::vector<ProfileSource> ProfileSourceManager::stale_enabled_sources(long long max_age_seconds) const
{
    const long long now = std::time(nullptr);
    std::vector<ProfileSource> result;
    for (const ProfileSource &source : sources())
        if (source.enabled && (source.last_sync == 0 || now - source.last_sync >= max_age_seconds))
            result.push_back(source);
    return result;
}

std::vector<ProfileSource> ProfileSourceManager::prusa_sources()
{
    return {
        {"prusa-fff", "Prusa Research FFF", "https://github.com/prusa3d/PrusaSlicer-settings-prusa-fff", ProfileSource::Format::Prusa},
        {"non-prusa-fff", "PrusaSlicer community FFF", "https://github.com/prusa3d/PrusaSlicer-settings-non-prusa-fff", ProfileSource::Format::Prusa}
    };
}

std::vector<ProfileSource> ProfileSourceManager::default_sources()
{
    std::vector<ProfileSource> result = prusa_sources();
    const std::vector<ProfileSource> native_sources {
        {"orca", "OrcaSlicer", "https://github.com/OrcaSlicer/OrcaSlicer", ProfileSource::Format::Orca},
        {"bambu", "Bambu Studio", "https://github.com/bambulab/BambuStudio", ProfileSource::Format::Orca},
        {"creality", "Creality Print", "https://github.com/CrealityOfficial/CrealityPrint", ProfileSource::Format::Orca},
        {"qidi", "QIDI Studio", "https://github.com/QIDITECH/QIDIStudio", ProfileSource::Format::Orca},
        {"superslicer", "SuperSlicer", "https://github.com/supermerill/SuperSlicer", ProfileSource::Format::Prusa}
    };
    result.insert(result.end(), native_sources.begin(), native_sources.end());
    return result;
}

std::string ProfileSourceManager::make_id(const std::string &name, const std::string &url)
{
    boost::uuids::detail::sha1 sha;
    const std::string value = name + '\n' + url;
    sha.process_bytes(value.data(), value.size());
    unsigned int digest[5] {};
    sha.get_digest(digest);
    std::ostringstream output;
    output << sanitize(name) << '_';
    output << std::hex << digest[0] << digest[1];
    return output.str();
}

} // namespace Slic3r
