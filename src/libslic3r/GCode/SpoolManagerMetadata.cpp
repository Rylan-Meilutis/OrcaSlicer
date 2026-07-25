#include "SpoolManagerMetadata.hpp"

#include <algorithm>
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

std::string sanitize_metadata_value(std::string value)
{
    std::replace_if(value.begin(), value.end(), [](char ch) {
        return ch == ']' || ch == ';' || ch == '\r' || ch == '\n';
    }, ' ');
    return value;
}

std::string metadata_line(const char *key, const std::vector<Filament> &filaments,
                          const std::function<const std::string &(const Filament &)> &value)
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

bool parse_selected_spools(const std::string &response, std::vector<Filament> &slots, std::string &error)
{
    try {
        const nlohmann::json root = nlohmann::json::parse(response);
        const auto selected = root.find("selectedSpools");
        if (selected == root.end() || !selected->is_array()) {
            error = "The OctoPrint SpoolManager response did not contain selectedSpools.";
            return false;
        }

        slots.clear();
        slots.reserve(selected->size());
        bool has_assigned_spool = false;
        for (const nlohmann::json &entry : *selected) {
            if (entry.is_null()) {
                slots.emplace_back();
                continue;
            }
            if (!entry.is_object()) {
                error = "The OctoPrint SpoolManager selectedSpools list contained an invalid slot.";
                slots.clear();
                return false;
            }

            Filament spool{
                entry.value("displayName", ""),
                entry.value("material", ""),
                entry.value("color", ""),
                entry.value("colorName", ""),
                entry.value("vendor", "")
            };
            has_assigned_spool |= !spool.name.empty();
            slots.emplace_back(std::move(spool));
        }

        if (slots.empty() || !has_assigned_spool) {
            error = "OctoPrint SpoolManager has no spools assigned to its tools or slots.";
            slots.clear();
            return false;
        }
    } catch (const std::exception &exception) {
        error = exception.what();
        slots.clear();
        return false;
    }
    return true;
}

std::string update_gcode_tail(const std::string &gcode_tail, const std::vector<Filament> &filaments)
{
    static const std::regex filament_type_re(R"((^|\n); filament_type = ([^\r\n]+))");
    static const std::regex filament_used_re(R"((^|\n); filament used \[mm\] = ([^\r\n]+))");
    static const std::regex filament_notes_re(R"((^|\n); filament_notes = ([^\r\n]+))");
    static const std::regex spool_name_re(R"(\[\s*sm_name\s*=\s*[^\]]*\])");

    std::string output = gcode_tail;
    static const std::regex metadata_re(
        R"((^|\n); spool_manager_filament_(names|materials|colors) = [^\r\n]*(\r?\n|$))");
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
    if (std::regex_search(output, notes_match, filament_notes_re)) {
        std::vector<std::string> notes;
        std::stringstream notes_stream(notes_match[2].str());
        for (std::string note; std::getline(notes_stream, note, ';');)
            notes.emplace_back(std::move(note));

        for (size_t index = 0; index < notes.size() && index < filaments.size(); ++index) {
            if (filaments[index].name.empty() || !std::regex_search(notes[index], spool_name_re))
                continue;
            notes[index] = std::regex_replace(
                notes[index], spool_name_re, "[sm_name = " + sanitize_metadata_value(filaments[index].name) + "]",
                std::regex_constants::format_first_only);
        }

        std::string replacement = notes_match[1].str() + "; filament_notes = ";
        for (size_t index = 0; index < notes.size(); ++index) {
            if (index > 0)
                replacement += ';';
            replacement += notes[index];
        }
        output.replace(static_cast<size_t>(notes_match.position()), static_cast<size_t>(notes_match.length()), replacement);
    }

    const std::string metadata =
        metadata_line("spool_manager_filament_names", filaments,
                      [](const Filament &filament) -> const std::string & { return filament.name; }) +
        metadata_line("spool_manager_filament_materials", filaments,
                      [](const Filament &filament) -> const std::string & { return filament.material; }) +
        metadata_line("spool_manager_filament_colors", filaments,
                      [](const Filament &filament) -> const std::string & { return filament.color; });

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
