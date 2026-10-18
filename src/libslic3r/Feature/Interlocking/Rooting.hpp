#pragma once

#include "libslic3r/ExPolygon.hpp"
#include <functional>

namespace Slic3r {
class PrintObject;

// Snapshot native material contacts before beam interlocking creates artificial
// ones. Roots attach below or alongside the beam band when present, preserving
// the beam band and reserving protective skins around the buried branches.
class Rooting {
public:
    explicit Rooting(const PrintObject &object);
    size_t apply(PrintObject &object, const std::function<void()> &cancel) const;
private:
    size_t apply_xy(PrintObject &object, const std::function<void()> &cancel,
                   const std::vector<std::vector<ExPolygons>> &beams) const;
    std::vector<std::vector<ExPolygons>> m_slices; // layer, region
};
}
