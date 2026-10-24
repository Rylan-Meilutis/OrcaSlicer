#pragma once

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <stdexcept>
#include <vector>

namespace Slic3r {
struct GCodeProcessorResult;
namespace DepositionMesh {

struct Options {
    double voxel_size_mm = .1;
    std::size_t max_memory_bytes = 512ULL * 1024 * 1024;
    std::uint64_t max_output_bytes = 512ULL * 1024 * 1024;
    std::size_t max_segments = 2000000;
    std::uint64_t max_sample_visits = 250000000;
    unsigned tile_size_cells = 64;
    std::size_t max_tiles = 262144;
    bool reduce_horizontal_faces = true;
};
struct Segment {
    Eigen::Vector3f start, end;
    float width, height;
};
struct Snapshot {
    std::vector<Segment> segments;
    // Retraction recovery and stationary purging cannot be distinguished or
    // reconstructed as spatial solids from extruder-only commands.
    std::size_t unrepresented_extruder_only_moves = 0;
    // Positive-E travel/wipe records include cleaner purging and vertical
    // extrusion. The processor deliberately does not assign them bead geometry.
    std::size_t unrepresented_motion_extrusion_moves = 0;
    // Explicit extrusion vertices coincident with the previous XY position
    // (e.g. a rounded arc endpoint followed by a tiny positive-E command).
    std::size_t unrepresented_stationary_extrusion_moves = 0;
};
struct Report {
    std::uint64_t triangle_count = 0, output_bytes = 0;
    double voxel_size_mm = 0.;
    std::uint64_t sample_visits = 0;
};
// Structured diagnostics avoid losing the offending sampled geometry when an
// export aborts. They are bounded independently of mesh or G-code size.
struct DegenerateTriangleError : std::runtime_error {
    std::array<std::array<float,3>,3> vertices;
    bool planar_reduced;
    DegenerateTriangleError(const std::array<std::array<float,3>,3> &, bool);
};

// A small, owned snapshot, independent of preview visibility and result lifetime.
// Arcs are already tessellated by GCodeProcessor. Includes every extrusion role.
Snapshot snapshot(const GCodeProcessorResult &, const Options & = {},
                  const std::function<void()> &checkpoint = {});

// Nominal rounded-bead union, not CAD or measured printed geometry. Each bead is
// segment + horizontal disc((width-height)/2) + sphere(height/2), shifted down
// height/2. Sloped paths remain sloped. Processor width/height are nominal and
// do not reproduce varying scarf flow or measured deposition cross-sections.
// The sampling
// resolution is explicit; unresolved gaps are not guaranteed to survive. No
// adaptivity or silent coarsening. STL has no material/interface information.
// Output must be a seekable private stream: publish only on successful return.
// OpenVDB meshing itself cannot be interrupted; budgets are checked beforehand.
Report export_stl(const Snapshot &, std::ostream &, const Options & = {},
                  const std::function<void(int)> &progress = {},
                  const std::function<void()> &checkpoint = {});

} // namespace DepositionMesh
} // namespace Slic3r
