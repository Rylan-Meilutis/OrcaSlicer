#include "ToolMapping.hpp"
#include "../Config.hpp"
#include "../PlaceholderParser.hpp"
#include "../Preset.hpp"
#include "../PrintConfig.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <ios>
#include <istream>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace Slic3r {

GCodeTemplateSnapshot GCodeTemplateSnapshot::capture(const std::string &name, const std::string &text,
    unsigned current_tool, const PlaceholderParser &parser, const DynamicConfig *overrides,
    const DynamicConfig &outputs, const PlaceholderParser::ContextData &context,
    std::shared_ptr<const GCodeTemplateEnvironment> environment)
{
    GCodeTemplateSnapshot snapshot;
    snapshot.name = name;
    snapshot.text = text;
    snapshot.current_tool = current_tool;
    if (!environment || (parser.external_config() && environment->external_config != *parser.external_config())) {
        auto owned = std::make_shared<GCodeTemplateEnvironment>();
        owned->parser_config = parser.config();
        if (parser.external_config()) owned->external_config = *parser.external_config();
        environment = std::move(owned);
    }
    snapshot.environment = std::move(environment);
    for (const auto &key : parser.config().keys()) {
        const auto *value = parser.config().option(key);
        const auto *base = snapshot.environment->parser_config.option(key);
        if (!base || *value != *base) snapshot.parser_changes.set_key_value(key, value->clone());
    }
    for (const auto &key : snapshot.environment->parser_config.keys())
        if (!parser.config().has(key)) snapshot.removed_parser_keys.push_back(key);
    if (overrides) snapshot.overrides = *overrides;
    snapshot.outputs_before = outputs;
    if (context.global_config) snapshot.globals_before = *context.global_config;
    snapshot.rng_before = context.rng;
    snapshot.inactive_before = context.inactive_global_variables;
    snapshot.legacy_machine_mode_index = context.legacy_machine_mode_index;
    return snapshot;
}

void GCodeTemplateSnapshot::finish(const std::string &output, const DynamicConfig &outputs,
                                  const PlaceholderParser::ContextData &context)
{
    rendered = output;
    outputs_after = outputs;
    if (context.global_config) globals_after = *context.global_config;
    rng_after = context.rng;
    inactive_after = context.inactive_global_variables;
}

namespace {
bool mapped_rows_equal(const ConfigOptionVectorBase &source, const ConfigOptionVectorBase &target,
                       const std::vector<int> &tools, const std::vector<int> *source_slots = nullptr)
{
    if (source.type() != target.type() || source.empty() || target.empty()) return false;
    const auto a = source.vserialize();
    const auto b = target.vserialize();
    for (std::size_t i = 0; i < tools.size(); ++i) {
        const auto si = source.size() == 1 ? 0 : (source_slots ? std::size_t((*source_slots)[i]) : i);
        const auto ti = target.size() == 1 ? 0 : std::size_t(tools[i]);
        if (si >= a.size() || ti >= b.size() || a[si] != b[ti]) return false;
    }
    return true;
}

DynamicConfig relocate_template_config(const DynamicConfig &source, const DynamicPrintConfig *target,
                                      const std::vector<int> &tools, std::size_t capacity)
{
    DynamicConfig result(source);
    if (target)
        for (const auto &key : source.keys())
            if (const auto *option = target->option(key)) result.set_key_value(key, option->clone());

    // Print's parser contains effective material retraction overrides, while
    // full_print_config retains the underlying nozzle defaults. Rebuild that
    // same effective context for the destination slots; copying the defaults
    // alone would change the macro's e_retracted output and writer state.
    if (target) {
        DynamicPrintConfig effective;
        t_config_option_keys changed;
        const auto *mapping = target->option<ConfigOptionInts>("filament_map");
        for (const auto &key : print_config_def.extruder_retract_keys()) {
            if (!source.has(key)) continue;
            const auto *base = target->option(key);
            const auto *material = dynamic_cast<const ConfigOptionVectorBase *>(target->option("filament_" + key));
            if (!base || !material || material->empty()) continue;
            std::vector<int> slots(material->size());
            for (std::size_t i = 0; i < slots.size(); ++i)
                slots[i] = mapping && i < mapping->size() ? mapping->values[i] - 1 : int(i);
            compute_filament_override_value(key, base, base, material, *target, changed, effective, slots);
        }
        result.apply(effective);
    }

    // These are runtime variables rather than settings. Do not translate
    // arbitrary integers: layer numbers, temperatures, times and literal
    // parameters in custom commands are not tool IDs.
    for (const char *key : {"initial_tool", "initial_extruder", "initial_filament_id", "current_extruder",
                            "current_filament_id", "previous_extruder", "next_extruder", "next_filament_id",
                            "filament_extruder_id",
                            "initial_no_support_tool", "initial_no_support_extruder", "initial_no_support_filament_id",
                            "initial_no_support_filament"}) {
        if (const auto *value = source.option<ConfigOptionInt>(key)) {
            if (value->value >= 0 && std::size_t(value->value) < tools.size())
                result.set_key_value(key, new ConfigOptionInt(tools[std::size_t(value->value)]));
        }
    }
    for (const auto &pair : {std::pair<const char *, const char *>{"current_hotend", "current_filament_id"},
                            {"current_extruder_id", "current_filament_id"}, {"next_hotend", "next_filament_id"},
                            {"initial_no_support_hotend", "initial_no_support_filament_id"}}) {
        const auto *physical = source.option<ConfigOptionInt>(pair.first);
        if (!physical) continue;
        const auto *logical = source.option<ConfigOptionInt>(pair.second);
        if (!logical && (std::string_view(pair.first) == "current_extruder_id" || std::string_view(pair.first) == "current_hotend"))
            logical = source.option<ConfigOptionInt>("filament_extruder_id");
        if (!logical && std::string_view(pair.first) == "next_hotend") logical = source.option<ConfigOptionInt>("next_extruder");
        if (!logical) throw std::runtime_error("A physical hotend placeholder has no logical material context.");
        if (logical->value < 0) continue;
        if (std::size_t(logical->value) >= tools.size()) throw std::runtime_error("Invalid logical material for a hotend placeholder.");
        result.set_key_value(pair.first, new ConfigOptionInt(tools[std::size_t(logical->value)]));
    }
    for (const char *key : {"first_non_support_tools", "first_non_support_filaments", "first_tools", "first_filaments"})
        if (const auto *list = source.option<ConfigOptionInts>(key)) {
            auto values = list->values;
            for (int &value : values)
                if (value >= 0 && std::size_t(value) < tools.size()) value = tools[std::size_t(value)];
            result.set_key_value(key, new ConfigOptionInts(values));
        }
    if (source.has("first_non_support_hotend")) {
        const auto *logical = source.option<ConfigOptionInts>("first_non_support_filaments");
        if (!logical) throw std::runtime_error("A hotend list has no logical material context.");
        auto values = logical->values;
        for (int &value : values) {
            if (value < 0) continue;
            if (std::size_t(value) >= tools.size()) throw std::runtime_error("Invalid logical material in a hotend list.");
            value = tools[std::size_t(value)];
        }
        result.set_key_value("first_non_support_hotend", new ConfigOptionInts(values));
    }
    for (const char *key : {"is_extruder_used", "e_position", "e_retracted", "e_restart_extra",
                            "extruded_volume", "extruded_weight", "tool_init",
                            "first_layer_temperature", "temperature", "first_layer_bed_temperature",
                            "bed_temperature", "flush_volumetric_speeds", "flush_temperatures", "retract_length"}) {
        const auto *values = dynamic_cast<const ConfigOptionVectorBase *>(source.option(key));
        if (!values || values->empty()) continue;
        const bool runtime_state = std::string_view(key) == "is_extruder_used" ||
            std::string_view(key) == "e_position" || std::string_view(key) == "e_retracted" ||
            std::string_view(key) == "e_restart_extra" || std::string_view(key) == "extruded_volume" ||
            std::string_view(key) == "extruded_weight" || std::string_view(key) == "tool_init";
        const auto count = runtime_state ? std::max(capacity, values->size()) : capacity;
        if (count > 256) throw std::runtime_error("Runtime tool array exceeds dispatch budget.");
        std::unique_ptr<ConfigOptionVectorBase> relocated(static_cast<ConfigOptionVectorBase *>(values->clone()));
        relocated->resize(count, values);
        // Unassigned runtime tool slots have no extrusion state and are unused.
        std::unique_ptr<ConfigOption> zero(values->clone());
        if (zero->type() == coFloats || zero->type() == coInts || zero->type() == coBools) {
            zero->deserialize("0");
            for (std::size_t i = 0; i < count; ++i) relocated->set_at(zero.get(), i, 0);
        } else {
            throw std::runtime_error("Unsupported runtime array in tool-mapped template.");
        }
        for (std::size_t i = 0; i < tools.size(); ++i) {
            if (runtime_state && i >= values->size()) continue; // Unused tools have no extrusion state.
            const auto index = !runtime_state && values->size() == 1 ? 0 : i;
            if (index >= values->size()) throw std::runtime_error("Incomplete runtime tool array in template.");
            relocated->set_at(values, std::size_t(tools[i]), index);
        }
        result.set_key_value(key, relocated.release());
    }
    // Overrides may supply layer-specific temperatures/flush settings instead
    // of the preset values. Relocate those rows rather than replacing them by
    // the target's ordinary (non-overridden) settings.
    if (!target) {
        for (const auto &key : Preset::filament_options()) {
            const auto *values = dynamic_cast<const ConfigOptionVectorBase *>(source.option(key));
            if (!values || values->empty()) continue;
            std::unique_ptr<ConfigOptionVectorBase> relocated(static_cast<ConfigOptionVectorBase *>(values->clone()));
            relocated->resize(capacity, values);
            for (std::size_t i = 0; i < capacity; ++i) relocated->set_at(values, i, 0);
            for (std::size_t i = 0; i < tools.size(); ++i) {
                const auto index = values->size() == 1 ? 0 : i;
                if (index >= values->size()) throw std::runtime_error("Incomplete material override in tool template.");
                relocated->set_at(values, std::size_t(tools[i]), index);
            }
            result.set_key_value(key, relocated.release());
        }
    }
    if (const auto *mapping = source.option<ConfigOptionInts>("filament_extruder_id")) {
        std::vector<int> values(capacity, 0);
        for (std::size_t i = 0; i < capacity; ++i) values[i] = int(i);
        for (std::size_t i = 0; i < tools.size(); ++i) {
            if (i >= mapping->size() || mapping->values[i] < 0 || std::size_t(mapping->values[i]) >= tools.size())
                throw std::runtime_error("Unsupported physical extruder map in tool template.");
            values[std::size_t(tools[i])] = tools[std::size_t(mapping->values[i])];
        }
        result.set_key_value("filament_extruder_id", new ConfigOptionInts(values));
    }
    return result;
}
}

std::string replay_tool_mapped_template(const GCodeTemplateSnapshot &snapshot,
    const DynamicPrintConfig &target, const std::vector<int> &tools)
{
    if (!snapshot.environment || tools.empty() || snapshot.current_tool >= tools.size())
        throw std::invalid_argument("Invalid custom-template tool mapping.");
    if (snapshot.text.find("nozzle_id") != std::string::npos || snapshot.text.find("nozzle_volume_types") != std::string::npos ||
        snapshot.text.find("filament_nozzle_map") != std::string::npos)
        throw std::runtime_error("Logical nozzle-group templates require normal regeneration.");
    std::set<int> assigned;
    for (int tool : tools)
        if (tool < 0 || tool >= 256 || !assigned.insert(tool).second)
            throw std::invalid_argument("Custom-template tools must be non-negative and one-to-one.");
    std::size_t capacity = std::size_t(*assigned.rbegin()) + 1;
    if (const auto *types = target.option<ConfigOptionStrings>("filament_type")) capacity = std::max(capacity, types->size());
    if (capacity > 256) throw std::runtime_error("Tool-dispatch slot budget exceeded.");

    DynamicConfig parser_config(snapshot.environment->parser_config);
    parser_config.apply(snapshot.parser_changes);
    for (const auto &key : snapshot.removed_parser_keys) parser_config.erase(key);
    const auto evaluate = [&](bool mapped) {
        const auto external = mapped ? relocate_template_config(snapshot.environment->external_config, &target, tools, capacity) : snapshot.environment->external_config;
        PlaceholderParser parser(&external);
        parser.config_writable() = mapped ? relocate_template_config(parser_config, &target, tools, capacity) : parser_config;
        const auto overrides = mapped ? relocate_template_config(snapshot.overrides, nullptr, tools, capacity) : snapshot.overrides;
        auto outputs = mapped ? relocate_template_config(snapshot.outputs_before, nullptr, tools, capacity) : snapshot.outputs_before;
        PlaceholderParser::ContextData context;
        context.rng = snapshot.rng_before;
        context.inactive_global_variables = snapshot.inactive_before;
        context.legacy_machine_mode_index = snapshot.legacy_machine_mode_index;
        if (mapped && !context.legacy_machine_mode_index)
            context.legacy_machine_mode_index = snapshot.current_tool;
        context.global_config = std::make_unique<DynamicConfig>(mapped ?
            relocate_template_config(snapshot.globals_before, nullptr, tools, capacity) : snapshot.globals_before);
        const std::string output = parser.process(snapshot.text,
            mapped ? unsigned(tools[snapshot.current_tool]) : snapshot.current_tool, &overrides, &outputs, &context);
        const auto expected_outputs = mapped ? relocate_template_config(snapshot.outputs_after, nullptr, tools, capacity) : snapshot.outputs_after;
        const auto expected_globals = mapped ? relocate_template_config(snapshot.globals_after, nullptr, tools, capacity) : snapshot.globals_after;
        if (outputs != expected_outputs || *context.global_config != expected_globals ||
            context.rng != snapshot.rng_after || context.inactive_global_variables != snapshot.inactive_after) {
            std::string changed;
            const auto describe = [&](const DynamicConfig &actual, const DynamicConfig &expected, const char *prefix) {
                std::set<std::string> keys;
                for (const auto &key : actual.keys()) keys.insert(key);
                for (const auto &key : expected.keys()) keys.insert(key);
                for (const auto &key : keys) {
                    const auto *a = actual.option(key), *b = expected.option(key);
                    if (!a || !b || *a != *b) changed += std::string(prefix) + key + " ";
                }
            };
            describe(outputs, expected_outputs, "output:");
            describe(*context.global_config, expected_globals, "global:");
            if (context.rng != snapshot.rng_after) changed += "random-generator ";
            if (context.inactive_global_variables != snapshot.inactive_after) changed += "inactive-globals ";
            throw std::runtime_error("Tool mapping changes custom-script state in " + snapshot.name + ": " + changed.substr(0, 1024));
        }
        if (!mapped && output != snapshot.rendered)
            throw std::runtime_error("Custom-script inputs changed; regenerate the print job.");
        return output;
    };
    evaluate(false);
    const auto output = evaluate(true);
    const auto motion_state_lines = [](const std::string &text) {
        std::vector<std::string> lines;
        std::istringstream stream(text);
        std::string line;
        while (std::getline(stream, line)) {
            line.erase(line.find(';') == std::string::npos ? line.size() : line.find(';'));
            const auto begin = line.find_first_not_of(" \t\r");
            if (begin == std::string::npos) continue;
            const auto end = line.find_first_of(" \t\r", begin);
            const auto command = line.substr(begin, end == std::string::npos ? end : end - begin);
            static const std::set<std::string> state_commands{
                "G0", "G1", "G2", "G3", "G4", "G10", "G11", "G12", "G17", "G18", "G19", "G20", "G21",
                "G27", "G28", "G29", "G90", "G91", "G92", "P0", "M17", "M18", "M82", "M83", "M84",
                "M201", "M203", "M204", "M205", "M220", "M221", "M400", "M569", "M593", "M906"};
            if (state_commands.count(command)) lines.push_back(line.substr(begin));
        }
        return lines;
    };
    // The parser's output variables are not the only way a macro can move the
    // machine. Unless a future final-state interpreter certifies an alternate
    // purge path, retain the original macro motions/control modes too.
    const auto mapped_motions = motion_state_lines(output), original_motions = motion_state_lines(snapshot.rendered);
    if (mapped_motions != original_motions)
        throw std::runtime_error("Tool mapping changes custom-script motion in " + snapshot.name + "; regenerate the print job.");
    return output;
}

bool equivalent_tool_mapping_config(const DynamicPrintConfig &source, const DynamicPrintConfig &target,
                                    const std::vector<int> &tools)
{
    if (tools.empty()) return false;
    std::set<int> assigned;
    for (int tool : tools)
        if (tool < 0 || tool >= 256 || !assigned.insert(tool).second) return false;

    // This path handles one independent physical nozzle per logical material.
    // Layer-dependent nozzle maps and virtual mixtures cannot be represented
    // by substituting a single physical tool number.
    const auto valid_map = [](const DynamicPrintConfig &config, std::size_t required, bool identity) {
        if (const auto *dynamic = config.option<ConfigOptionBool>("enable_filament_dynamic_map"); dynamic && dynamic->value)
            return false;
        if (const auto *map = config.option<ConfigOptionInts>("filament_map")) {
            if (map->size() < required) return false;
            for (std::size_t i = 0; i < map->size(); ++i)
                if (map->values[i] <= 0 || (identity && map->values[i] != int(i + 1))) return false;
        }
        if (const auto *mixed = config.option<ConfigOptionBools>("filament_is_mixed"))
            if (std::any_of(mixed->values.begin(), mixed->values.end(), [](bool value) { return value; })) return false;
        if (const auto *nozzles = config.option<ConfigOptionInts>("filament_nozzle_map")) {
            const bool unresolved = std::all_of(nozzles->values.begin(), nozzles->values.end(), [](int value) { return value == 0; });
            if (!unresolved) {
                const auto *extruders = config.option<ConfigOptionInts>("filament_map");
                if (!extruders || nozzles->size() != extruders->size()) return false;
                // The static single-nozzle-per-extruder engine writes this
                // derived cache after slicing. It is not a different nozzle
                // assignment from the all-zero pre-slice placeholder.
                for (std::size_t i = 0; i < nozzles->size(); ++i)
                    if (nozzles->values[i] != extruders->values[i] - 1) return false;
            }
        }
        // Multiple nozzle-volume variants per material need a column map,
        // not the one material-row-per-slot projection supported here.
        if (const auto *indices = config.option<ConfigOptionInts>("filament_self_index")) {
            if (indices->size() != required) return false;
            for (std::size_t i = 0; i < indices->size(); ++i)
                if (indices->values[i] != int(i + 1)) return false;
        }
        return true;
    };
    std::size_t capacity = std::size_t(*assigned.rbegin()) + 1;
    if (const auto *map = target.option<ConfigOptionInts>("filament_map")) capacity = std::max(capacity, map->size());
    if (capacity > 256) return false;
    if (!valid_map(source, tools.size(), false) || !valid_map(target, capacity, true)) return false;
    std::vector<int> source_slots(tools.size());
    const auto *source_map = source.option<ConfigOptionInts>("filament_map");
    const auto *source_nozzles = source.option<ConfigOptionFloats>("nozzle_diameter");
    for (std::size_t i = 0; i < tools.size(); ++i) {
        source_slots[i] = source_map ? source_map->values[i] - 1 : int(i);
        if (source_nozzles && std::size_t(source_slots[i]) >= source_nozzles->size()) return false;
    }

    // The printer arrays remain in physical-slot order. Even if the complete
    // printer configuration is unchanged, moving to a different row may change
    // nozzle diameter, offsets, retract or extrusion limits.
    for (const std::string &key : Preset::nozzle_options()) {
        const auto *a = dynamic_cast<const ConfigOptionVectorBase *>(source.option(key));
        const auto *b = dynamic_cast<const ConfigOptionVectorBase *>(target.option(key));
        if ((a == nullptr) != (b == nullptr)) return false;
        if (a) {
            const auto *source_nozzles = source.option<ConfigOptionFloats>("nozzle_diameter");
            const auto *target_nozzles = target.option<ConfigOptionFloats>("nozzle_diameter");
            if ((source_nozzles && a->size() != 1 && a->size() != source_nozzles->size()) ||
                (target_nozzles && b->size() != 1 && b->size() != target_nozzles->size())) return false;
        }
        if (a && !mapped_rows_equal(*a, *b, tools, &source_slots)) return false;
    }

    std::set<std::string> material_keys(Preset::filament_options().begin(), Preset::filament_options().end());
    for (const char *key : {"filament_settings_id", "filament_ids", "filament_colour", "filament_colour_type", "filament_multi_colour", "filament_volume_map",
                            "filament_is_mixed", "filament_mixed_components", "filament_mixed_sublayer_ratios",
                            "filament_mixed_gradient", "filament_mixed_gradient_range", "filament_mixed_gradient_curve",
                            "filament_mixed_gradient_per_part"})
        material_keys.insert(key);
    auto roles = project_filament_role_keys();
    roles.emplace_back("wipe_tower_filament");
    const std::set<std::string> role_keys(roles.begin(), roles.end());
    std::set<std::string> keys;
    for (const auto &key : source.keys()) keys.insert(key);
    for (const auto &key : target.keys()) keys.insert(key);
    for (const std::string &key : keys) {
        const ConfigOption *a = source.option(key);
        const ConfigOption *b = target.option(key);
        // Host credentials and already-resolved UI bindings never enter the
        // generated path. The concrete role settings are compared below.
        static const std::set<std::string> dispatch_only_keys{
            "print_host", "print_host_webui", "printhost_apikey", "printhost_cafile", "printhost_user",
            "printhost_password", "printhost_port", "project_filament_bindings", "project_filament_roles",
            "filament_map_mode", "compatible_machine_expression_group", "compatible_process_expression_group",
            "inherits_group", "different_settings_to_system", "print_compatible_printers"};
        if (dispatch_only_keys.count(key)) continue;
        const bool indexed = key == "flush_volumes_matrix" || key == "flush_volumes_vector" ||
            key == "first_layer_print_sequence" || key == "other_layers_print_sequence";
        if (a && b && *a == *b && !material_keys.count(key) && !role_keys.count(key) && !indexed) continue;
        if (!a || !b || a->type() != b->type()) return false;
        if (key == "filament_map" || key == "filament_nozzle_map" || key == "filament_self_index") continue; // Validated above.
        if (key == "flush_volumes_matrix" || key == "flush_volumes_vector") {
            const auto *av = dynamic_cast<const ConfigOptionFloats *>(a);
            const auto *bv = dynamic_cast<const ConfigOptionFloats *>(b);
            if (!av || !bv) return false;
            if (av->empty() && bv->empty()) continue;
            const bool matrix = key == "flush_volumes_matrix";
            const auto old_stride = tools.size() * (matrix ? tools.size() : 2);
            const auto new_stride = capacity * (matrix ? capacity : 2);
            if (av->size() % old_stride || bv->size() % new_stride || av->size() / old_stride != bv->size() / new_stride)
                return false;
            for (std::size_t block = 0; block < av->size() / old_stride; ++block)
                for (std::size_t i = 0; i < tools.size(); ++i)
                    for (std::size_t j = 0; j < (matrix ? tools.size() : 2); ++j)
                        if (av->values[block * old_stride + i * (matrix ? tools.size() : 2) + j] !=
                            bv->values[block * new_stride + std::size_t(tools[i]) * (matrix ? capacity : 2) + (matrix ? std::size_t(tools[j]) : j)])
                            return false;
            continue;
        }
        if (key == "first_layer_print_sequence" || key == "other_layers_print_sequence") {
            const auto *av = dynamic_cast<const ConfigOptionInts *>(a);
            const auto *bv = dynamic_cast<const ConfigOptionInts *>(b);
            if (!av || !bv || av->size() != bv->size()) return false;
            std::size_t stride = av->size();
            const bool ranges = key == "other_layers_print_sequence";
            if (ranges) {
                const auto *count = source.option<ConfigOptionInt>("other_layers_print_sequence_nums");
                if (!count || count->value <= 0) { if (*a == *b) continue; return false; }
                if (av->size() % std::size_t(count->value)) return false;
                stride = av->size() / std::size_t(count->value);
                if (stride <= 2) return false;
            }
            for (std::size_t i = 0; i < av->size(); ++i) {
                const auto value = av->values[i];
                if (ranges && i % stride < 2) { if (bv->values[i] != value) return false; continue; }
                if (value < 0 || std::size_t(value) > tools.size() || bv->values[i] != (value == 0 ? 0 : tools[std::size_t(value - 1)] + 1))
                    return false;
            }
            continue;
        }
        if (material_keys.count(key)) {
            const auto *av = dynamic_cast<const ConfigOptionVectorBase *>(a);
            const auto *bv = dynamic_cast<const ConfigOptionVectorBase *>(b);
            if (av && bv) {
                if ((av->size() > 1 && av->size() != tools.size()) || (bv->size() > 1 && bv->size() != capacity))
                    return false; // Variant/grouped arrays need their own dimensional mapping.
                if ((av->empty() && bv->empty()) || mapped_rows_equal(*av, *bv, tools)) continue;
            } else if (*a == *b) continue;
            return false;
        }
        if (role_keys.count(key)) {
            const auto *av = dynamic_cast<const ConfigOptionInt *>(a);
            const auto *bv = dynamic_cast<const ConfigOptionInt *>(b);
            if (!av || !bv || av->value < 0 || std::size_t(av->value) > tools.size()) return false;
            if (bv->value == (av->value == 0 ? 0 : tools[std::size_t(av->value - 1)] + 1)) continue;
        }
        return false;
    }
    return true;
}

bool can_reuse_identity_tool_mapping(const DynamicPrintConfig &source, const DynamicPrintConfig &target,
                                    const std::vector<int> &tools)
{
    if (tools.empty()) return false;
    for (std::size_t i = 0; i < tools.size(); ++i)
        if (tools[i] < 0 || std::size_t(tools[i]) != i) return false;
    // The exact comparison retains support for existing non-independent jobs;
    // the projection only relaxes unused physical-slot material padding.
    if (source == target) return true;
    if (const auto *map = source.option<ConfigOptionInts>("filament_map"))
        for (std::size_t i = 0; i < tools.size(); ++i)
            if (i >= map->size() || map->values[i] != int(i + 1)) return false;
    return equivalent_tool_mapping_config(source, target, tools);
}

namespace {
void fingerprint_bytes(std::uint64_t &hash, std::string_view bytes)
{
    for (unsigned char byte : bytes) { hash ^= byte; hash *= 1099511628211ULL; }
}

std::size_t unsigned_token(std::string_view token)
{
    std::size_t value = 0;
    const auto parsed = std::from_chars(token.data(), token.data() + token.size(), value);
    if (token.empty() || parsed.ec != std::errc{} || parsed.ptr != token.data() + token.size())
        throw std::runtime_error("Unsupported tool operand; regenerate the print job.");
    return value;
}

std::string terminated(std::string text)
{
    if (!text.empty() && text.back() != '\n') text += '\n';
    return text;
}

std::string retain_processed_macro_edits(const std::string &original, const std::string &processed,
                                         const std::string &replayed)
{
    const auto lines = [](std::string_view text) {
        std::vector<std::string_view> result;
        while (!text.empty()) {
            const auto newline = text.find('\n');
            const auto length = newline == std::string_view::npos ? text.size() : newline + 1;
            result.push_back(text.substr(0, length));
            text.remove_prefix(length);
        }
        return result;
    };
    const auto before = lines(original), after = lines(replayed);
    if (before.size() != after.size())
        throw std::runtime_error("Mapped custom script changes line structure; regenerate the print job.");

    struct Replacements {
        std::vector<std::string_view> lines;
        std::size_t consumed = 0;
    };
    std::map<std::string_view, Replacements> edits;
    for (std::size_t i = 0; i < before.size(); ++i)
        if (before[i] != after[i]) edits.try_emplace(before[i]);

    // Account for every occurrence, including unchanged copies of a changed
    // line. Exact counts and original occurrence order are required: a deleted,
    // rewritten or additional matching line makes relocation ambiguous.
    for (std::size_t i = 0; i < before.size(); ++i)
        if (auto edit = edits.find(before[i]); edit != edits.end())
            edit->second.lines.push_back(after[i]);

    // The final stream is generator-sealed. Keep ALL of its other bytes:
    // FanMover can relocate fans and split G1 moves; timing inserts M73; the
    // deposition-bounds pass updates M555. Reconstructing the whole template
    // would discard those transformations. Motion/writer-state equality was
    // checked on the original and replayed expansions before reaching here.
    std::string result;
    result.reserve(processed.size());
    for (const auto line : lines(processed)) {
        if (auto edit = edits.find(line); edit != edits.end()) {
            auto &replacement = edit->second;
            if (replacement.consumed >= replacement.lines.size())
                throw std::runtime_error("A changed custom-script line is ambiguous after processing.");
            result += replacement.lines[replacement.consumed++];
        } else {
            result += line;
        }
    }
    for (const auto &edit : edits)
        if (edit.second.consumed != edit.second.lines.size())
            throw std::runtime_error("A changed custom-script line was removed or rewritten after processing.");
    return result;
}
}

void seal_gcode_dispatch(std::istream &source, GCodeDispatchData &data, const std::function<void()> &checkpoint)
{
    data.sealed = false;
    data.source_bytes = 0;
    data.source_fingerprint = 14695981039346656037ULL;
    std::array<char, 64 * 1024> buffer;
    while (true) {
        if (checkpoint) checkpoint();
        source.read(buffer.data(), std::streamsize(buffer.size()));
        const auto count = source.gcount();
        if (count > 0) {
            data.source_bytes += std::uint64_t(count);
            fingerprint_bytes(data.source_fingerprint, std::string_view(buffer.data(), std::size_t(count)));
        }
        if (source.bad() || (source.fail() && !source.eof())) throw std::runtime_error("Cannot certify generated G-code.");
        if (source.eof()) break;
    }
    data.sealed = true;
}

GCodeDispatchPlan prepare_tool_mapped_gcode(std::istream &source, const GCodeDispatchData &data,
    const DynamicPrintConfig &target, const std::vector<int> &tools, const std::string &full_config,
    const std::function<void(std::uint64_t)> &progress)
{
    if (!data.valid || !data.sealed || tools.empty()) throw std::runtime_error("Missing tool-dispatch certificate; regenerate the print job.");
    for (const char *key : {"post_process", "slicing_pipeline_plugin"})
        if (const auto *scripts = target.option<ConfigOptionStrings>(key))
            for (const auto &script : scripts->values)
                if (!script.empty()) throw std::runtime_error("Post-processing requires normal tool-dispatch regeneration.");
    for (const char *key : {"is_imex", "single_extruder_multi_material", "enable_prime_tower", "enable_pre_heating", "enable_filament_dynamic_map"})
        if (const auto *enabled = target.option<ConfigOptionBool>(key); enabled && enabled->value)
            throw std::runtime_error("This printer/toolpath mode requires normal tool-dispatch regeneration.");
    const auto *flavor = target.option<ConfigOptionEnum<GCodeFlavor>>("gcode_flavor");
    if (!flavor || (flavor->value != gcfMarlinLegacy && flavor->value != gcfMarlinFirmware))
        throw std::runtime_error("This G-code dialect requires normal tool-dispatch regeneration.");
    std::set<int> assigned;
    for (int tool : tools)
        if (tool < 0 || tool >= 256 || !assigned.insert(tool).second) throw std::invalid_argument("Invalid physical tool mapping.");
    GCodeDispatchPlan plan;
    std::set<std::size_t> seen;
    std::size_t retained_bytes = 0;
    constexpr std::size_t memory_limit = 32 * 1024 * 1024;
    const auto add_edit = [&](std::uint64_t offset, std::string expected, std::string replacement) {
        if (expected == replacement) return;
        const auto bytes = expected.size() + replacement.size() + sizeof(GCodeDispatchEdit);
        if (bytes > memory_limit - retained_bytes)
            throw std::runtime_error("Tool-dispatch edit budget exceeded; regenerate the print job.");
        retained_bytes += bytes;
        plan.edits.push_back({offset, std::move(expected), std::move(replacement)});
    };
    const auto read_line = [&](std::string &line) {
        line.clear();
        // Bound even malformed/user-edited files without a newline. std::getline
        // alone would allocate an arbitrarily large line before checking it.
        char ch;
        while (source.get(ch)) {
            line += ch;
            if (line.size() > 1024 * 1024) throw std::runtime_error("G-code line exceeds dispatch budget.");
            if (ch == '\n') break;
        }
        if (source.bad() || (source.fail() && !source.eof())) throw std::runtime_error("Reading G-code dispatch input failed.");
        if (line.empty()) return false;
        fingerprint_bytes(plan.source_fingerprint, line);
        plan.source_bytes += line.size();
        if (progress && (plan.source_bytes % (64 * 1024) < line.size())) progress(plan.source_bytes);
        return true;
    };
    const auto mapped_id = [&](std::string_view token, unsigned base) {
        const auto value = unsigned_token(token);
        if (value < base || value - base >= tools.size())
            throw std::runtime_error("A native command addresses an uncaptured tool; regenerate the print job.");
        return std::to_string(std::uint64_t(tools[value - base]) + base);
    };
    std::string line;
    while (read_line(line)) {
        const auto offset = plan.source_bytes - line.size();
        constexpr std::string_view begin = "; ORCA_DISPATCH_BEGIN ";
        if (line.compare(0, begin.size(), begin) == 0) {
            std::string_view id_text(line.data() + begin.size(), line.size() - begin.size());
            if (!id_text.empty() && id_text.back() == '\n') id_text.remove_suffix(1);
            const auto id = unsigned_token(id_text);
            if (id >= data.templates.size() || !seen.insert(id).second)
                throw std::runtime_error("Invalid or duplicate custom-script certificate.");
            const auto body_offset = plan.source_bytes;
            const std::string end = "; ORCA_DISPATCH_END " + std::to_string(id) + "\n";
            std::string body;
            bool ended = false;
            while (read_line(line)) {
                if (line == end) { ended = true; break; }
                if (body.size() + line.size() > memory_limit) throw std::runtime_error("Custom-script span exceeds dispatch budget.");
                body += line;
            }
            if (!ended)
                throw std::runtime_error("A custom script was changed after slicing; regenerate the print job.");
            auto replayed = terminated(replay_tool_mapped_template(data.templates[id], target, tools));
            try {
                replayed = retain_processed_macro_edits(terminated(data.templates[id].rendered), body, replayed);
            } catch (const std::runtime_error &error) {
                throw std::runtime_error("Template " + data.templates[id].name + ": " + error.what());
            }
            add_edit(body_offset, std::move(body), std::move(replayed));
            continue;
        }
        if (line == "; CONFIG_BLOCK_START\n") {
            const auto body_offset = plan.source_bytes;
            std::string body;
            bool ended = false;
            while (read_line(line)) {
                if (line == "; CONFIG_BLOCK_END\n") { ended = true; break; }
                if (body.size() + line.size() > memory_limit) throw std::runtime_error("Configuration exceeds dispatch budget.");
                body += line;
            }
            if (!ended || full_config.empty()) throw std::runtime_error("Missing mapped configuration metadata.");
            // The mapped Print has not been processed. Its all-zero derived
            // nozzle cache must describe the physical identity slots in the
            // exported metadata, just as normal mapped generation would.
            std::istringstream config_lines(full_config);
            std::string replacement, config_line;
            while (std::getline(config_lines, config_line)) {
                if (config_line.compare(0, 24, "; filament_nozzle_map = ") == 0) {
                    const auto *mapping = target.option<ConfigOptionInts>("filament_map");
                    if (!mapping || mapping->empty()) throw std::runtime_error("Missing mapped nozzle metadata.");
                    auto slots = mapping->values;
                    for (std::size_t i = 0; i < slots.size(); ++i) {
                        if (slots[i] != int(i + 1)) throw std::runtime_error("Non-identity physical nozzle metadata.");
                        --slots[i];
                    }
                    config_line = "; filament_nozzle_map = " + ConfigOptionInts(slots).serialize();
                }
                replacement += config_line + '\n';
                if (replacement.size() > memory_limit) throw std::runtime_error("Configuration exceeds dispatch budget.");
            }
            add_edit(body_offset, std::move(body), std::move(replacement));
            continue;
        }
        if (line.compare(0, 12, "; filament: ") == 0) {
            std::string replacement = "; filament: ";
            std::istringstream fields(line.substr(12));
            std::string field;
            bool first = true;
            while (std::getline(fields, field, ',')) {
                const auto a = field.find_first_not_of(" \t\r\n");
                const auto b = field.find_last_not_of(" \t\r\n");
                if (a == std::string::npos) throw std::runtime_error("Invalid filament metadata.");
                if (!first) replacement += ',';
                replacement += mapped_id(std::string_view(field).substr(a, b - a + 1), 1);
                first = false;
            }
            add_edit(offset, line, replacement + '\n');
            continue;
        }
        bool metadata = false;
        for (const char *key : {"filament_density", "filament_diameter"}) {
            const std::string prefix = std::string("; ") + key + ": ";
            if (line.compare(0, prefix.size(), prefix) == 0) {
                const auto *option = target.option(key);
                if (!option) throw std::runtime_error("Missing mapped filament metadata.");
                add_edit(offset, line, prefix + option->serialize() + '\n');
                metadata = true;
                break;
            }
        }
        if (metadata) continue;
        for (const char *prefix : {"; filament used [mm] = ", "; filament used [cm3] = ", "; filament used [g] = ", "; filament cost = "}) {
            const std::string_view p(prefix);
            if (line.compare(0, p.size(), p) != 0) continue;
            std::vector<std::string> values(std::size_t(*assigned.rbegin()) + 1, "0");
            std::istringstream fields(line.substr(p.size()));
            std::string field;
            std::size_t index = 0;
            while (std::getline(fields, field, ',')) {
                if (index >= tools.size()) throw std::runtime_error("Filament statistics refer to an uncaptured tool.");
                const auto a = field.find_first_not_of(" \t\r\n");
                const auto b = field.find_last_not_of(" \t\r\n");
                if (a == std::string::npos) throw std::runtime_error("Invalid filament statistics.");
                values[std::size_t(tools[index++])] = field.substr(a, b - a + 1);
            }
            std::string replacement(prefix);
            for (std::size_t i = 0; i < values.size(); ++i) replacement += (i ? ", " : "") + values[i];
            add_edit(offset, line, replacement + '\n');
            metadata = true;
            break;
        }
        if (metadata) continue;

        std::string_view code(line);
        code = code.substr(0, code.find(';'));
        const auto first = code.find_first_not_of(" \t\r\n");
        if (first == std::string_view::npos) continue;
        code.remove_prefix(first);
        const auto command_end = code.find_first_of(" \t\r\n");
        const auto command = code.substr(0, command_end);
        if (command.size() > 1 && command.front() == 'T') {
            add_edit(offset + first + 1, std::string(command.substr(1)), mapped_id(command.substr(1), 0));
            continue;
        }
        // Commands outside captured custom spans are generator-owned. Keep a
        // small explicit dialect: future tool-bearing commands fail closed.
        static const std::set<std::string_view> unchanged_commands{
            "G0", "G1", "G2", "G3", "G4", "G10", "G11", "G17", "G18", "G19", "G20", "G21",
            "G28", "G29", "G90", "G91", "G92", "M73", "M82", "M83", "M106", "M107", "M140",
            "M141", "M190", "M191", "M201", "M203", "M204", "M205", "M400", "M486"};
        char operand = '\0';
        if (command == "M104" || command == "M109" || command == "M900") operand = 'T';
        else if (command == "M572") operand = 'D';
        else if (!unchanged_commands.count(command))
            throw std::runtime_error("Unsupported native command in fast tool dispatch: " + std::string(command));
        if (command_end == std::string_view::npos) continue;
        std::size_t pos = command_end;
        while (pos < code.size()) {
            pos = code.find_first_not_of(" \t\r\n", pos);
            if (pos == std::string_view::npos) break;
            const auto end = code.find_first_of(" \t\r\n", pos);
            const auto token = code.substr(pos, end == std::string_view::npos ? end : end - pos);
            if (token.front() == operand)
                add_edit(offset + first + pos + 1, std::string(token.substr(1)), mapped_id(token.substr(1), 0));
            else if ((token.front() == 'T' && command != "M204" && command != "M205") ||
                     (command == "G10" && token.front() == 'P'))
                throw std::runtime_error("Unsupported native tool operand; regenerate the print job.");
            if (end == std::string_view::npos) break;
            pos = end;
        }
    }
    if (seen.size() != data.templates.size()) throw std::runtime_error("A custom-script certificate is missing from the final G-code.");
    if (plan.source_bytes != data.source_bytes || plan.source_fingerprint != data.source_fingerprint)
        throw std::runtime_error("G-code changed after generation; regenerate the print job.");
    if (progress) progress(plan.source_bytes);
    return plan;
}

void apply_tool_mapped_gcode(std::istream &source, std::ostream &target, const GCodeDispatchPlan &plan,
    const std::function<void(std::uint64_t)> &progress)
{
    std::array<char, 64 * 1024> buffer;
    std::uint64_t consumed = 0, fingerprint = 14695981039346656037ULL;
    const auto transfer = [&](std::uint64_t end, bool write) {
        if (end < consumed || end > plan.source_bytes) throw std::runtime_error("Invalid dispatch edit range.");
        while (consumed < end) {
            const auto count = std::streamsize(std::min<std::uint64_t>(buffer.size(), end - consumed));
            source.read(buffer.data(), count);
            if (source.gcount() != count) throw std::runtime_error("Source G-code changed during dispatch.");
            fingerprint_bytes(fingerprint, std::string_view(buffer.data(), std::size_t(count)));
            if (write) target.write(buffer.data(), count);
            consumed += std::uint64_t(count);
            if (!target) throw std::runtime_error("Writing mapped G-code failed.");
            if (progress) progress(consumed);
        }
    };
    for (const auto &edit : plan.edits) {
        transfer(edit.offset, true);
        // The whole-file fingerprint also covers unchanged motion, not just
        // replacement spans. Publish the private destination only after it
        // matches the read-only preflight snapshot.
        transfer(edit.offset + edit.expected.size(), false);
        target.write(edit.replacement.data(), std::streamsize(edit.replacement.size()));
    }
    transfer(plan.source_bytes, true);
    if (source.peek() != std::char_traits<char>::eof() || fingerprint != plan.source_fingerprint)
        throw std::runtime_error("Source G-code changed during dispatch.");
    target.flush();
    if (!target) throw std::runtime_error("Writing mapped G-code failed.");
}

void remap_gcode_bindings(std::istream &source, std::ostream &target,
                         const std::vector<GCodeToolBinding> &bindings,
                         const std::vector<int> &tools,
                         const std::function<void(std::uint64_t)> &progress)
{
    std::set<int> used;
    if (tools.empty()) throw std::invalid_argument("Empty physical tool mapping.");
    for (int tool : tools)
        if (tool < 0 || !used.insert(tool).second)
            throw std::invalid_argument("Tool mapping must be non-negative and one-to-one.");
    std::uint64_t previous_end = 0;
    for (const auto &binding : bindings) {
        if (binding.logical_tool >= tools.size() || binding.index_base > 1 ||
            binding.logical_tool > std::numeric_limits<std::size_t>::max() - binding.index_base)
            throw std::invalid_argument("Invalid G-code tool binding.");
        const auto width = std::to_string(binding.logical_tool + binding.index_base).size();
        if (binding.offset < previous_end || binding.offset > std::numeric_limits<std::uint64_t>::max() - width)
            throw std::invalid_argument("G-code bindings overlap or are not ordered.");
        previous_end = binding.offset + width;
    }
    std::array<char, 64 * 1024> buffer;
    std::uint64_t consumed = 0;
    char previous = '\0';
    const auto report = [&] {
        if (!target) throw std::runtime_error("Writing remapped G-code failed.");
        if (progress) progress(consumed);
    };
    const auto copy_to = [&](std::uint64_t end) {
        while (consumed < end) {
            const auto count = std::streamsize(std::min<std::uint64_t>(buffer.size(), end - consumed));
            source.read(buffer.data(), count);
            if (source.gcount() != count)
                throw std::runtime_error("G-code binding points beyond the source file.");
            target.write(buffer.data(), count);
            previous = buffer[std::size_t(count - 1)];
            consumed += std::uint64_t(count);
            report();
        }
    };
    report();
    for (const auto &binding : bindings) {
        copy_to(binding.offset);
        if ((previous >= '0' && previous <= '9') || previous == '.' || previous == '-' || previous == '+')
            throw std::runtime_error("G-code binding starts inside a numeric value.");
        const std::string expected = std::to_string(binding.logical_tool + binding.index_base);
        source.read(buffer.data(), std::streamsize(expected.size()));
        if (source.gcount() != std::streamsize(expected.size()) ||
            !std::equal(expected.begin(), expected.end(), buffer.begin()))
            throw std::runtime_error("G-code tool binding no longer matches the source file.");
        consumed += expected.size();
        previous = expected.back();
        const auto following = source.peek();
        if ((following >= '0' && following <= '9') || following == '.')
            throw std::runtime_error("G-code binding ends inside a numeric value.");
        // Promote before adding one: a one-based slot must not overflow int.
        target << std::to_string(std::uint64_t(tools[binding.logical_tool]) + binding.index_base);
        report();
    }
    while (true) {
        source.read(buffer.data(), std::streamsize(buffer.size()));
        const auto count = source.gcount();
        if (count > 0) {
            target.write(buffer.data(), count);
            consumed += std::uint64_t(count);
            report();
        }
        if (source.bad() || (source.fail() && !source.eof()))
            throw std::runtime_error("Reading G-code for tool remapping failed.");
        if (source.eof()) break;
    }
    target.flush();
    report();
}

} // namespace Slic3r
