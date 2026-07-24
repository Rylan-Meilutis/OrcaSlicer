#pragma once

#include <string>
#include <vector>

#include <boost/filesystem/path.hpp>

namespace Slic3r::SpoolManagerMetadata {

// Replaces [sm_name=...] markers in the serialized filament notes and pads
// per-filament usage metadata to the number of configured filaments.
std::string update_gcode_tail(const std::string &gcode_tail, const std::vector<std::string> &spool_names);

// Updates only the trailing slicer configuration block and preserves the rest
// of the potentially large G-code file without loading it all into memory.
bool update_gcode_file(const boost::filesystem::path &path, const std::vector<std::string> &spool_names,
                       std::string &error);

} // namespace Slic3r::SpoolManagerMetadata
