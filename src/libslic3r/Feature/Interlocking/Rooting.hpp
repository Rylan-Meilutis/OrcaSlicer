#pragma once

#include "libslic3r/ExPolygon.hpp"
#include <functional>

namespace Slic3r {
class PrintObject;

// Snapshot native material contacts before beam interlocking creates artificial
// ones. Roots attach below the lowest beam course when present, preserving
// the beam band and reserving protective skins around the buried branches.
class Rooting {
public:
    explicit Rooting(const PrintObject &object);
    size_t apply(PrintObject &object, const std::function<void()> &cancel) const;
private:
    std::vector<std::vector<ExPolygons>> m_slices; // layer, region
};
}
