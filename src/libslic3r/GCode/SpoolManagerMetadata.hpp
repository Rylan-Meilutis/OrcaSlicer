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
    std::string vendor;
    std::string spool_id;
    std::string provider;
    // Inventory namespace reported by provider-neutral backends.  RME uses
    // this to distinguish internal rolls from delegated inventory providers.
    std::string inventory_provider;
};

inline void to_json(nlohmann::json &json, const Filament &filament)
{
    json = {
        {"name", filament.name},
        {"material", filament.material},
        {"color", filament.color},
        {"color_name", filament.color_name},
        {"vendor", filament.vendor},
        {"spool_id", filament.spool_id},
        {"provider", filament.provider},
        {"inventory_provider", filament.inventory_provider}
    };
}

inline void from_json(const nlohmann::json &json, Filament &filament)
{
    filament.name       = json.value("name", "");
    filament.material   = json.value("material", "");
    filament.color      = json.value("color", "");
    filament.color_name = json.value("color_name", "");
    filament.vendor     = json.value("vendor", "");
    filament.spool_id   = json.value("spool_id", "");
    filament.provider   = json.value("provider", "");
    filament.inventory_provider = json.value("inventory_provider", "");
}

// Parses tool assignments reported by the legacy SpoolManager plugin, the
// OctoPrint Spoolman plugin, or a provider implementing the generic
// { data: { tools, spools } } compatibility response. Empty slots are
// preserved so each vector index continues to match the OctoPrint tool index.
bool parse_selected_spools(const std::string &response, std::vector<Filament> &slots,
                           std::string &error,
                           const std::string &fallback_provider = {});

// Resolves an explicitly configured profile name. Mapping entries use
// "provider:spool-id=preset" (or "spool-id=preset") for a roll and
// "vendor|material=preset" (or "material=preset") for a material family.
// Exact-roll mappings take priority over material mappings and the default.
std::string mapped_profile_name(const Filament &filament,
                                const std::vector<std::string> &spool_mappings,
                                const std::vector<std::string> &material_mappings,
                                const std::string &default_profile);

// Replaces [sm_name=...] markers in the serialized filament notes and pads
// per-filament usage metadata to the number of configured filaments. Explicit
// name/material/color arrays are also embedded for host-side job validation.
std::string update_gcode_tail(const std::string &gcode_tail, const std::vector<Filament> &filaments);

// Updates only the trailing slicer configuration block and preserves the rest
// of the potentially large G-code file without loading it all into memory.
bool update_gcode_file(const boost::filesystem::path &path, const std::vector<Filament> &filaments,
                       std::string &error);

} // namespace Slic3r::SpoolManagerMetadata
