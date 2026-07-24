#pragma once

#include <string>
#include <vector>

#include "Polygon.hpp"

namespace Slic3r {

class ConfigBase;

struct SequentialGantrySlice
{
    double   height = 0.;
    bool     is_box = false;
    Polygons polygons;
};

struct SequentialGantryGeometry
{
    std::string                        model_path;
    std::vector<SequentialGantrySlice> slices;

    bool empty() const { return slices.empty(); }
    double conservative_clearance_radius() const;
    double first_box_height() const;
};

// Loads a custom machine definition when configured, otherwise matches the
// printer notes against the bundled Prusa-style geometry catalog.
SequentialGantryGeometry load_sequential_gantry_geometry(const ConfigBase &config);

} // namespace Slic3r
