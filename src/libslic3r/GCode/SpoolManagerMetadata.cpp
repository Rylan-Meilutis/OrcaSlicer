#include "SpoolManagerMetadata.hpp"

#include <algorithm>
#include <fstream>
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

std::string sanitize_spool_name(std::string name)
{
    std::replace_if(name.begin(), name.end(), [](char ch) {
        return ch == ']' || ch == '\r' || ch == '\n';
    }, ' ');
    return name;
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

std::string update_gcode_tail(const std::string &gcode_tail, const std::vector<std::string> &spool_names)
{
    static const std::regex filament_type_re(R"((^|\n); filament_type = ([^\r\n]+))");
    static const std::regex filament_used_re(R"((^|\n); filament used \[mm\] = ([^\r\n]+))");
    static const std::regex filament_notes_re(R"((^|\n); filament_notes = ([^\r\n]+))");
    static const std::regex spool_name_re(R"(\[\s*sm_name\s*=\s*[^\]]*\])");

    std::smatch notes_match;
    if (!std::regex_search(gcode_tail, notes_match, filament_notes_re))
        return gcode_tail;

    std::string output = gcode_tail;

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

    if (!std::regex_search(output, notes_match, filament_notes_re))
        return output;

    std::vector<std::string> notes;
    std::stringstream notes_stream(notes_match[2].str());
    for (std::string note; std::getline(notes_stream, note, ';');)
        notes.emplace_back(std::move(note));

    for (size_t index = 0; index < notes.size() && index < spool_names.size(); ++index) {
        if (spool_names[index].empty() || !std::regex_search(notes[index], spool_name_re))
            continue;
        notes[index] = std::regex_replace(
            notes[index], spool_name_re, "[sm_name = " + sanitize_spool_name(spool_names[index]) + "]",
            std::regex_constants::format_first_only);
    }

    std::string replacement = notes_match[1].str() + "; filament_notes = ";
    for (size_t index = 0; index < notes.size(); ++index) {
        if (index > 0)
            replacement += ';';
        replacement += notes[index];
    }
    output.replace(static_cast<size_t>(notes_match.position()), static_cast<size_t>(notes_match.length()), replacement);
    return output;
}

bool update_gcode_file(const fs::path &path, const std::vector<std::string> &spool_names, std::string &error)
{
    const auto [tail, tail_offset] = read_tail(path);
    if (tail_offset < 0) {
        error = "Unable to read G-code file";
        return false;
    }

    const std::string updated_tail = update_gcode_tail(tail, spool_names);
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
