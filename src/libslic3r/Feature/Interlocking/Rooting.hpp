#pragma once

#include "libslic3r/ExPolygon.hpp"
#include <cstddef>
#include <vector>
#include <functional>

namespace Slic3r {
class PrintObject;

// Snapshot native material contacts before beam interlocking creates artificial
// ones. Roots attach below or alongside the beam band when present, preserving
// the beam band and reserving protective skins around the buried branches.
class Rooting {
public:
    // Called on the slicing thread with completed/total layers (or region pairs
    // for XY roots). Geometry workers never invoke UI progress callbacks.
    using Progress = std::function<void(bool xy, size_t completed, size_t total)>;
    struct ContactResult {
        size_t region;
        size_t layer;
        ExPolygon contact;
        bool placed = false;
    };
    explicit Rooting(const PrintObject &object);
    size_t apply(PrintObject &object, const std::function<void()> &cancel, const Progress &progress = {},
                 std::vector<ContactResult> *z_contacts = nullptr) const;
private:
    size_t apply_xy(PrintObject &object, const std::function<void()> &cancel,
                   const std::vector<std::vector<ExPolygons>> &beams, const Progress &progress) const;
    std::vector<std::vector<ExPolygons>> m_slices; // layer, region
};
}
