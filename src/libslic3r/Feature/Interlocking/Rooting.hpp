#pragma once

#include "libslic3r/ExPolygon.hpp"
#include <functional>

namespace Slic3r {
class PrintObject;

// Snapshot native material contacts before beam interlocking creates artificial
// ones. Application after beams rejects conflicting roots instead of overwriting
// other interlocks or making disconnected, clipped branches.
class Rooting {
public:
    explicit Rooting(const PrintObject &object);
    size_t apply(PrintObject &object, const std::function<void()> &cancel) const;
private:
    std::vector<std::vector<ExPolygons>> m_slices; // layer, region
};
}
