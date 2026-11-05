#pragma once

#include <string>
#include <vector>

namespace Slic3r {
class Print;
struct GCodeProcessorResult;

// Package an already generated job without referring back to editable GUI
// presets/plate caches, which may use different logical material numbering.
bool store_sliced_gcode_3mf(const std::string &path, const Print &print,
                          GCodeProcessorResult &result, const std::string &gcode_path);
// Read-only metadata preflight, before allocating the mapped output file.
// Throws when packaging requires regeneration; does not certify G-code bytes.
void validate_tool_mapped_sliced_gcode_3mf(const Print &source_print, const Print &mapped_print,
    const GCodeProcessorResult &source_result, const std::vector<int> &tools);
// Package a separately certified fast-remapped file without copying/relabeling
// the logical preview. Source statistics remain tied to logical materials;
// mapped_print supplies the mapped model/config, not unsliced statistics.
// Throws when metadata cannot safely be projected (including opaque warnings).
// The caller must certify and apply the GCodeDispatchPlan before calling this.
bool store_tool_mapped_sliced_gcode_3mf(const std::string &path, const Print &source_print,
    const Print &mapped_print, GCodeProcessorResult &source_result,
    const std::vector<int> &tools, const std::string &mapped_gcode_path);
}
