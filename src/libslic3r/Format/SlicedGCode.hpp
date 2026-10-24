#pragma once

#include <string>

namespace Slic3r {
class Print;
struct GCodeProcessorResult;

// Package an already generated job without referring back to editable GUI
// presets/plate caches, which may use different logical material numbering.
bool store_sliced_gcode_3mf(const std::string &path, const Print &print,
                          GCodeProcessorResult &result, const std::string &gcode_path);
}
