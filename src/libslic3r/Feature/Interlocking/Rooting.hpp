#pragma once

#include "libslic3r/ExPolygon.hpp"
#include <functional>

namespace Slic3r {
class PrintObject;

// Snapshot native material contacts before beam interlocking creates artificial
// ones. Accepted roots reserve local envelopes after beams, keeping their floor
// and connecting stems intact while leaving unrelated interlocks unchanged.
class Rooting {
public:
    explicit Rooting(const PrintObject &object);
    size_t apply(PrintObject &object, const std::function<void()> &cancel) const;
private:
    std::vector<std::vector<ExPolygons>> m_slices; // layer, region
};
}
