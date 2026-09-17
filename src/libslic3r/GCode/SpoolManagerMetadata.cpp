#include "SpoolManagerMetadata.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <functional>
#include <regex>
#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/nowide/fstream.hpp>
#include "../Utils.hpp"

namespace fs = boost::filesystem;

namespace Slic3r::SpoolManagerMetadata {

namespace {

constexpr size_t tail_line_count = 1000;
constexpr size_t copy_buffer_size = 64 * 1024;

std::string scalar_string(const nlohmann::json &value)
{
    if (value.is_string())
        return value.get<std::string>();
    if (value.is_number_integer())
        return std::to_string(value.get<long long>());
    if (value.is_number_unsigned())
        return std::to_string(value.get<unsigned long long>());
    return {};
}

std::string object_string(const nlohmann::json &object, std::initializer_list<const char *> keys)
{
    if (!object.is_object())
        return {};
    for (const char *key : keys) {
        const auto value = object.find(key);
        if (value != object.end()) {
            const std::string result = scalar_string(*value);
            if (!result.empty())
                return result;
        }
    }
    return {};
}

Filament filament_from_spool(const nlohmann::json &spool, const std::string &provider)
{
    if (!spool.is_object())
        return {};

    const auto nested_spool = spool.find("spool");
    if (nested_spool != spool.end() && nested_spool->is_object())
        return filament_from_spool(*nested_spool, provider);
    const auto filament_it = spool.find("filament");
    const auto profile_it = spool.find("profile");
    const nlohmann::json &filament = filament_it != spool.end() && filament_it->is_object() ? *filament_it :
                                     profile_it != spool.end() && profile_it->is_object() ? *profile_it : spool;
    std::string vendor;
    const auto vendor_it = filament.find("vendor");
    if (vendor_it != filament.end()) {
        vendor = vendor_it->is_object() ? object_string(*vendor_it, {"name", "displayName"})
                                        : scalar_string(*vendor_it);
    }
    if (vendor.empty())
        vendor = object_string(spool, {"vendor", "manufacturer"});

    Filament result;
    result.name       = object_string(spool, {"displayName", "display_name", "name"});
    if (result.name.empty())
        result.name   = object_string(filament, {"name", "displayName", "display_name"});
    result.material   = object_string(filament, {"material", "type", "filament_type"});
    result.color      = object_string(filament, {"color", "color_hex", "colour", "colour_hex"});
    result.color_name = object_string(filament, {"colorName", "color_name", "colourName", "colour_name"});
    result.vendor     = std::move(vendor);
    result.spool_id   = object_string(spool, {"spoolId", "spool_id", "id"});
    result.provider   = provider;
    result.inventory_provider = object_string(spool, {"provider", "inventory_provider"});
    if (result.name.empty() && (!result.vendor.empty() || !result.material.empty())) {
        result.name = result.vendor;
        if (!result.name.empty() && !result.material.empty())
            result.name += ' ';
        result.name += result.material;
    }
    return result;
}

const nlohmann::json *find_spool(const nlohmann::json &spools, const std::string &id)
{
    if (!spools.is_array() || id.empty())
        return nullptr;
    for (const nlohmann::json &spool : spools)
        if (spool.is_object() && object_string(spool, {"spoolId", "spool_id", "id"}) == id)
            return &spool;
    return nullptr;
}

bool finish_slots(std::vector<Filament> &slots, std::string &error)
{
    if (!slots.empty() && std::any_of(slots.begin(), slots.end(), [](const Filament &filament) {
            return !filament.name.empty() || !filament.spool_id.empty();
        }))
        return true;
    slots.clear();
    error = "The OctoPrint filament provider has no spools assigned to its tools or slots.";
    return false;
}

std::string sanitize_metadata_value(std::string value)
{
    std::replace_if(value.begin(), value.end(), [](char ch) {
        return ch == ']' || ch == ';' || ch == '\r' || ch == '\n';
    }, ' ');
    return value;
}

std::string normalized_identifier_component(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    std::replace_if(value.begin(), value.end(), [](unsigned char ch) {
        return std::isspace(ch) != 0;
    }, '-');
    return sanitize_metadata_value(std::move(value));
}

std::string validation_spool_name(const Filament &filament)
{
    if (filament.spool_id.empty())
        return filament.name;

    const std::string provider = normalized_identifier_component(filament.provider);
    if (provider == "spoolman")
        return "spoolman:" + sanitize_metadata_value(filament.spool_id);
    if (provider.find("rme") != std::string::npos) {
        std::string inventory_provider = normalized_identifier_component(filament.inventory_provider);
        if (inventory_provider.empty() || inventory_provider.find("rme") != std::string::npos)
            inventory_provider = "internal";
        return "rme:" + inventory_provider + ':' + sanitize_metadata_value(filament.spool_id);
    }
    return filament.name;
}

std::string metadata_line(const char *key, const std::vector<Filament> &filaments,
                          const std::function<std::string(const Filament &)> &value)
{
    std::string output = std::string("; ") + key + " = ";
    for (size_t index = 0; index < filaments.size(); ++index) {
        if (index > 0)
            output += ';';
        output += sanitize_metadata_value(value(filaments[index]));
    }
    return output + '\n';
}

std::pair<std::string, std::streamoff> read_tail(const fs::path &path)
{
    boost::nowide::ifstream input(path.string(), std::ios::binary);
    if (!input)
        return {{}, -1};

    input.seekg(0, std::ios::end);
    const std::streamoff file_size = input.tellg();
    std::streamoff start = file_size;
    size_t newline_count = 0;
    std::vector<char> buffer(copy_buffer_size);

    while (start > 0 && newline_count <= tail_line_count) {
        const std::streamoff read_size = std::min<std::streamoff>(start, buffer.size());
        start -= read_size;
        input.seekg(start);
        input.read(buffer.data(), read_size);
        newline_count += std::count(buffer.begin(), buffer.begin() + input.gcount(), '\n');
    }

    input.clear();
    input.seekg(start);
    std::string tail(static_cast<size_t>(file_size - start), '\0');
    input.read(tail.data(), tail.size());

    if (newline_count > tail_line_count) {
        size_t skip = newline_count - tail_line_count;
        size_t offset = 0;
        while (skip > 0 && offset < tail.size()) {
            if (tail[offset++] == '\n')
                --skip;
        }
        start += static_cast<std::streamoff>(offset);
        tail.erase(0, offset);
    }

    return {std::move(tail), start};
}

} // namespace

bool parse_selected_spools(const std::string &response, std::vector<Filament> &slots,
                           std::string &error,
                           const std::string &fallback_provider)
{
    try {
        const nlohmann::json root = nlohmann::json::parse(response);
        slots.clear();

        const nlohmann::json *payload = &root;
        const auto data = root.find("data");
        if (data != root.end() && data->is_object())
            payload = &*data;
        const std::string reported_provider = object_string(root, {"provider"});
        const bool rme_filament_report = object_string(root, {"schema"}) == "rme-filament-report-v1";
        const auto provider_or = [&](const char *native_provider) {
            // In the RME report, the top-level provider names the backing
            // inventory (for example `spoolmanager`), not the HTTP backend.
            if (rme_filament_report)
                return fallback_provider.empty() ? std::string("RME compatibility") : fallback_provider;
            if (!reported_provider.empty())
                return reported_provider;
            return fallback_provider.empty() ? std::string(native_provider) :
                                               fallback_provider;
        };

        // Legacy OllisGit SpoolManager and the generic compatibility format
        // may return the complete assigned spool objects directly.
        auto selected = payload->find("selectedSpools");
        if (selected == payload->end())
            selected = payload->find("selected_spools");
        if (selected != payload->end() && selected->is_array()) {
            slots.reserve(selected->size());
            for (const nlohmann::json &entry : *selected)
                slots.emplace_back(entry.is_null() ? Filament{} :
                    filament_from_spool(entry, provider_or("SpoolManager")));
            return finish_slots(slots, error);
        }

        // OllisGit FilamentManager returns one selection per tool and nests
        // its material data under spool.profile.
        const auto selections = payload->find("selections");
        if (selections != payload->end() && selections->is_array()) {
            for (size_t position = 0; position < selections->size(); ++position) {
                const nlohmann::json &selection = (*selections)[position];
                size_t index = position;
                if (selection.is_object()) {
                    const std::string tool = object_string(selection, {"tool", "toolIdx", "tool_index"});
                    if (!tool.empty()) {
                        try { index = static_cast<size_t>(std::stoul(tool)); }
                        catch (...) {}
                    }
                }
                if (slots.size() <= index)
                    slots.resize(index + 1);
                slots[index] = filament_from_spool(
                    selection, provider_or("FilamentManager"));
            }
            return finish_slots(slots, error);
        }

        const auto spools_it = payload->find("spools");
        const nlohmann::json empty_spools = nlohmann::json::array();
        const nlohmann::json &spools = spools_it != payload->end() ? *spools_it : empty_spools;

        // The maintained OctoPrint Spoolman plugin stores tool assignments as
        // selectedSpoolIds and reports the Spoolman inventory separately.
        auto ids = payload->find("selectedSpoolIds");
        if (ids == payload->end())
            ids = payload->find("selected_spool_ids");
        if (ids != payload->end() && ids->is_object()) {
            size_t slot_count = 0;
            for (const auto &[tool, assignment] : ids->items()) {
                try { slot_count = std::max(slot_count, static_cast<size_t>(std::stoul(tool) + 1)); }
                catch (...) { continue; }
            }
            slots.resize(slot_count);
            for (const auto &[tool, assignment] : ids->items()) {
                size_t index;
                try { index = static_cast<size_t>(std::stoul(tool)); }
                catch (...) { continue; }
                const std::string id = assignment.is_object()
                    ? object_string(assignment, {"spoolId", "spool_id", "id"}) : scalar_string(assignment);
                if (const nlohmann::json *spool = find_spool(spools, id))
                    slots[index] = filament_from_spool(
                        *spool, provider_or("Spoolman"));
                else {
                    slots[index].spool_id = id;
                    slots[index].provider = provider_or("Spoolman");
                }
            }
            return finish_slots(slots, error);
        }

        // Provider-neutral response for RME compatibility and third-party
        // plugins: tools may be an array or an object keyed by tool number.
        const auto tools = payload->find("tools");
        if (tools != payload->end() && (tools->is_array() || tools->is_object())) {
            const auto append_tool = [&](size_t index, const nlohmann::json &tool) {
                if (slots.size() <= index)
                    slots.resize(index + 1);
                const nlohmann::json *source = &tool;
                if (tool.is_object()) {
                    const auto nested = tool.find("spool");
                    if (nested != tool.end() && nested->is_object())
                        source = &*nested;
                    else if (const nlohmann::json *spool = find_spool(
                                 spools, object_string(tool, {"spoolId", "spool_id", "id"})))
                        source = spool;
                }
                slots[index] = filament_from_spool(
                    *source, provider_or("OctoPrint"));
                if (tool.is_object()) {
                    const std::string inventory_provider = object_string(
                        tool, {"provider", "inventory_provider"});
                    if (!inventory_provider.empty())
                        slots[index].inventory_provider = inventory_provider;
                    else if (rme_filament_report)
                        slots[index].inventory_provider = reported_provider;
                    if (slots[index].spool_id.empty())
                        slots[index].spool_id = object_string(tool, {"spoolId", "spool_id", "id"});
                }
            };
            if (tools->is_array()) {
                for (size_t index = 0; index < tools->size(); ++index)
                    append_tool(index, (*tools)[index]);
            } else {
                for (const auto &[tool, value] : tools->items()) {
                    try { append_tool(static_cast<size_t>(std::stoul(tool)), value); }
                    catch (...) { continue; }
                }
            }
            return finish_slots(slots, error);
        }

        error = "The OctoPrint filament provider response did not contain recognized tool assignments.";
        return false;
    } catch (const std::exception &exception) {
        error = exception.what();
        slots.clear();
        return false;
    }
    return true;
}

bool parse_spool_inventory(const std::string &response, std::vector<Filament> &inventory,
                           std::string &error,
                           const std::string &fallback_provider)
{
    try {
        const nlohmann::json root = nlohmann::json::parse(response);
        const nlohmann::json *payload = &root;
        if (const auto data = root.find("data"); data != root.end() && data->is_object())
            payload = &*data;

        const std::string reported_provider = object_string(root, {"provider"});
        const bool rme_report = object_string(root, {"schema"}) == "rme-filament-report-v1";
        const auto provider_or = [&](const char *native_provider) {
            if (rme_report)
                return fallback_provider.empty() ? std::string("RME compatibility") : fallback_provider;
            if (!reported_provider.empty())
                return reported_provider;
            return fallback_provider.empty() ? std::string(native_provider) : fallback_provider;
        };

        const nlohmann::json *spools = nullptr;
        for (const char *key : {"allSpools", "all_spools", "spools"}) {
            const auto found = payload->find(key);
            if (found != payload->end() && found->is_array()) {
                spools = &*found;
                break;
            }
        }
        if (spools == nullptr) {
            error = "The OctoPrint filament provider did not report a spool inventory.";
            inventory.clear();
            return false;
        }

        const char *native_provider = payload->contains("allSpools") || payload->contains("all_spools") ?
            "SpoolManager" : "Spoolman";
        inventory.clear();
        inventory.reserve(spools->size());
        for (const nlohmann::json &entry : *spools) {
            Filament filament = filament_from_spool(entry, provider_or(native_provider));
            if (!filament.name.empty() || !filament.spool_id.empty())
                inventory.emplace_back(std::move(filament));
        }
        if (inventory.empty()) {
            error = "The OctoPrint filament provider reported an empty spool inventory.";
            return false;
        }
        return true;
    } catch (const std::exception &exception) {
        error = exception.what();
        inventory.clear();
        return false;
    }
}

std::string mapped_profile_name(const Filament &filament,
                                const std::vector<std::string> &spool_mappings,
                                const std::vector<std::string> &material_mappings,
                                const std::string &default_profile)
{
    const auto normalized = [](std::string value) {
        std::transform(value.begin(), value.end(), value.begin(),
                       [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
        value.erase(std::remove_if(value.begin(), value.end(),
                                   [](unsigned char ch) { return std::isspace(ch) != 0; }),
                    value.end());
        return value;
    };
    const auto lookup = [&normalized](const std::vector<std::string> &mappings,
                                      const std::vector<std::string> &keys) {
        for (const std::string &mapping : mappings) {
            const size_t separator = mapping.find('=');
            if (separator == std::string::npos)
                continue;
            const std::string key = normalized(mapping.substr(0, separator));
            if (std::find(keys.begin(), keys.end(), key) == keys.end())
                continue;
            std::string profile = mapping.substr(separator + 1);
            const size_t first = profile.find_first_not_of(" \t");
            const size_t last = profile.find_last_not_of(" \t");
            return first == std::string::npos ? std::string() : profile.substr(first, last - first + 1);
        }
        return std::string();
    };

    const std::string spool_id = normalized(filament.spool_id);
    if (!spool_id.empty()) {
        std::vector<std::string> keys{spool_id};
        const std::string provider = normalized(filament.provider);
        if (!provider.empty())
            keys.insert(keys.begin(), provider + ':' + spool_id);
        if (std::string profile = lookup(spool_mappings, keys); !profile.empty())
            return profile;
    }

    const std::string material = normalized(filament.material);
    if (!material.empty()) {
        std::vector<std::string> keys{material};
        const std::string vendor = normalized(filament.vendor);
        if (!vendor.empty())
            keys.insert(keys.begin(), vendor + '|' + material);
        if (std::string profile = lookup(material_mappings, keys); !profile.empty())
            return profile;
    }
    return default_profile;
}

std::string update_gcode_tail(const std::string &gcode_tail, const std::vector<Filament> &filaments)
{
    static const std::regex filament_type_re(R"((^|\n); filament_type = ([^\r\n]+))");
    static const std::regex filament_used_re(R"((^|\n); filament used \[mm\] = ([^\r\n]+))");
    static const std::regex filament_notes_re(R"((^|\n); filament_notes = ([^\r\n]+))");
    static const std::regex spool_name_re(R"(\[\s*sm_name\s*=\s*[^\]]*\])");

    std::string output = gcode_tail;
    static const std::regex metadata_re(
        R"((^|\n); spool_manager_filament_(names|materials|colors|spool_ids|providers|validation_ids) = [^\r\n]*)");
    output = std::regex_replace(output, metadata_re, "$1");

    size_t filament_count = 0;
    std::smatch type_match;
    if (std::regex_search(gcode_tail, type_match, filament_type_re)) {
        filament_count = 1 + static_cast<size_t>(
            std::count(type_match[2].first, type_match[2].second, ';'));
    }

    std::smatch used_match;
    if (filament_count > 0 && std::regex_search(output, used_match, filament_used_re)) {
        std::vector<std::string> usage;
        std::stringstream usage_stream(used_match[2].str());
        for (std::string value; std::getline(usage_stream, value, ',');)
            usage.emplace_back(value);
        usage.resize(std::max(usage.size(), filament_count), " 0");

        std::string replacement = used_match[1].str() + "; filament used [mm] =";
        for (size_t index = 0; index < usage.size(); ++index)
            replacement += (index == 0 ? " " : ",") + usage[index];
        output.replace(static_cast<size_t>(used_match.position()), static_cast<size_t>(used_match.length()), replacement);
    }

    std::smatch notes_match;
    std::string generated_notes_metadata;
    if (std::regex_search(output, notes_match, filament_notes_re)) {
        std::vector<std::string> notes;
        std::stringstream notes_stream(notes_match[2].str());
        for (std::string note; std::getline(notes_stream, note, ';');)
            notes.emplace_back(std::move(note));
        notes.resize(std::max(notes.size(), std::min(filament_count, filaments.size())));

        for (size_t index = 0; index < notes.size() && index < filaments.size(); ++index) {
            const std::string validation_name = validation_spool_name(filaments[index]);
            if (validation_name.empty())
                continue;
            const std::string marker = "[sm_name = " + sanitize_metadata_value(validation_name) + "]";
            if (std::regex_search(notes[index], spool_name_re))
                notes[index] = std::regex_replace(notes[index], spool_name_re, marker,
                                                  std::regex_constants::format_first_only);
            else {
                if (!notes[index].empty() && !std::isspace(static_cast<unsigned char>(notes[index].back())))
                    notes[index] += ' ';
                notes[index] += marker;
            }
        }

        std::string replacement = notes_match[1].str() + "; filament_notes = ";
        for (size_t index = 0; index < notes.size(); ++index) {
            if (index > 0)
                replacement += ';';
            replacement += notes[index];
        }
        output.replace(static_cast<size_t>(notes_match.position()), static_cast<size_t>(notes_match.length()), replacement);
    } else if (!filaments.empty()) {
        generated_notes_metadata = "; filament_notes = ";
        const size_t note_count = filament_count == 0 ? filaments.size() :
                                  std::min(filament_count, filaments.size());
        for (size_t index = 0; index < note_count; ++index) {
            if (index > 0)
                generated_notes_metadata += ';';
            const std::string validation_name = validation_spool_name(filaments[index]);
            if (!validation_name.empty())
                generated_notes_metadata += "[sm_name = " + sanitize_metadata_value(validation_name) + "]";
        }
        generated_notes_metadata += '\n';
    }

    const std::string metadata =
        generated_notes_metadata +
        metadata_line("spool_manager_filament_names", filaments,
                      [](const Filament &filament) { return filament.name; }) +
        metadata_line("spool_manager_filament_materials", filaments,
                      [](const Filament &filament) { return filament.material; }) +
        metadata_line("spool_manager_filament_colors", filaments,
                      [](const Filament &filament) { return filament.color; }) +
        metadata_line("spool_manager_filament_spool_ids", filaments,
                      [](const Filament &filament) { return filament.spool_id; }) +
        metadata_line("spool_manager_filament_providers", filaments,
                      [](const Filament &filament) { return filament.provider; }) +
        metadata_line("spool_manager_filament_validation_ids", filaments,
                      [](const Filament &filament) { return validation_spool_name(filament); });

    size_t insertion_position = output.size();
    if (std::regex_search(output, notes_match, filament_notes_re))
        insertion_position = static_cast<size_t>(notes_match.position() + notes_match.length(1));
    else if (std::regex_search(output, type_match, filament_type_re))
        insertion_position = static_cast<size_t>(type_match.position() + type_match.length(1));
    else if (!output.empty() && output.back() != '\n') {
        output += '\n';
        insertion_position = output.size();
    }
    output.insert(insertion_position, metadata);
    return output;
}

bool update_gcode_file(const fs::path &path, const std::vector<Filament> &filaments, std::string &error)
{
    const auto [tail, tail_offset] = read_tail(path);
    if (tail_offset < 0) {
        error = "Unable to read G-code file";
        return false;
    }

    const std::string updated_tail = update_gcode_tail(tail, filaments);
    if (updated_tail == tail)
        return true;

    const fs::path temporary_path = path.parent_path() / fs::unique_path(path.filename().string() + ".%%%%-%%%%.tmp");
    boost::nowide::ifstream input(path.string(), std::ios::binary);
    boost::nowide::ofstream output(temporary_path.string(), std::ios::binary | std::ios::trunc);
    if (!input || !output) {
        error = "Unable to create temporary G-code file";
        return false;
    }

    std::vector<char> buffer(copy_buffer_size);
    std::streamoff remaining = tail_offset;
    while (remaining > 0) {
        const std::streamsize amount = static_cast<std::streamsize>(
            std::min<std::streamoff>(remaining, buffer.size()));
        input.read(buffer.data(), amount);
        const std::streamsize read = input.gcount();
        if (read <= 0)
            break;
        output.write(buffer.data(), read);
        remaining -= read;
    }
    output.write(updated_tail.data(), static_cast<std::streamsize>(updated_tail.size()));
    input.close();
    output.close();

    if (remaining != 0 || !output) {
        fs::remove(temporary_path);
        error = "Unable to rewrite G-code metadata";
        return false;
    }

    if (const std::error_code ec = rename_file(temporary_path.string(), path.string())) {
        fs::remove(temporary_path);
        error = ec.message();
        return false;
    }
    return true;
}

} // namespace Slic3r::SpoolManagerMetadata
