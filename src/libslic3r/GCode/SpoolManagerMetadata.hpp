#pragma once

#include <string>
#include <vector>

#include <boost/filesystem/path.hpp>
#include <nlohmann/json.hpp>

namespace Slic3r::SpoolManagerMetadata {

struct Filament
{
    std::string name;
    std::string material;
    std::string color;
    std::string color_name;
};

inline void to_json(nlohmann::json &json, const Filament &filament)
{
    json = {
        {"name", filament.name},
        {"material", filament.material},
        {"color", filament.color},
        {"color_name", filament.color_name}
    };
}

inline void from_json(const nlohmann::json &json, Filament &filament)
{
    filament.name       = json.value("name", "");
    filament.material   = json.value("material", "");
    filament.color      = json.value("color", "");
    filament.color_name = json.value("color_name", "");
}

// Parses the SpoolManager loadSpoolsByQuery response. The allSpools list is
// intentionally used so the user may map any loaded project filament.
bool parse_spools(const std::string &response, std::vector<Filament> &spools, std::string &error);

// Replaces [sm_name=...] markers in the serialized filament notes and pads
// per-filament usage metadata to the number of configured filaments. Explicit
// name/material/color arrays are also embedded for host-side job validation.
std::string update_gcode_tail(const std::string &gcode_tail, const std::vector<Filament> &filaments);

// Updates only the trailing slicer configuration block and preserves the rest
// of the potentially large G-code file without loading it all into memory.
bool update_gcode_file(const boost::filesystem::path &path, const std::vector<Filament> &filaments,
                       std::string &error);

} // namespace Slic3r::SpoolManagerMetadata
