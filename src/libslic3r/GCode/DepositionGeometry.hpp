#pragma once

#include "DepositionMesh.hpp"
#include <cstdint>
#include <functional>
#include <iosfwd>

namespace Slic3r::DepositionMesh {

enum class GeometryFormat { Step, StepTessellated, Obj, Stl };
struct GeometryOptions {
    // Each segment or joined path is a closed, faceted nominal bead. Adjacent beads
    // may overlap: no Boolean union, bonding assumption or gap filling occurs.
    unsigned angular_samples = 16;
    // Join exactly collinear, contiguous planar runs with identical bead
    // dimensions and material identities. This changes no occupied volume.
    bool merge_collinear = true;
    // Maximum centreline displacement in mm. Zero disables approximate path
    // compaction. This is separate from angular cross-section faceting error.
    double path_tolerance_mm = 0.;
    // Join gentle planar runs into closed path meshes instead of repeating
    // rounded caps at every sample. The additional join displacement is
    // bounded by path_tolerance_mm. Other paths can still overlap.
    bool join_paths = false;
    std::uint64_t max_output_bytes = 512ULL * 1024 * 1024;
};

// Simplifies only contiguous planar paths with identical dimensions and
// extruder/color identities. Bounds the two-sided centreline Hausdorff error;
// it is not a Boolean union or a material-interface reconstruction.
Snapshot compact_paths(const Snapshot&, double tolerance_mm,
                       const std::function<void()>& checkpoint = {});

// Streams one bounded path mesh at a time, independent of whole-print mesh complexity.
// STEP uses millimetres and one product per segment or joined path (arcs
// are already tessellated). OBJ uses millimetres and one named object per path.
// OBJ coordinates are rounded to 1 nm (maximum 0.5 nm per axis), with compact
// relative vertex indices. STEP retains its independent 13-digit precision.
// Step is AP214 planar B-rep; StepTessellated is compact AP242 indexed geometry
// requiring tessellated-STEP support in the receiving application.
// Stl streams separate closed beads into one binary STL (not a Boolean union),
// and requires a seekable stream for its triangle-count header.
// Output is private until successful return; a failure can leave partial data.
Report export_geometry(const Snapshot &, std::ostream &, GeometryFormat,
                       const GeometryOptions & = {},
                       const std::function<void(int)> &progress = {},
                       const std::function<void()> &checkpoint = {});

} // namespace Slic3r::DepositionMesh
