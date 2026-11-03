#pragma once
#include "Point.hpp"
#include "SectionIntegration.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include "libslic3r.h"
#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

namespace Slic3r::Thermoelastic {
// Internal exposed portion of a certified prism edge. Material lies to the
// left of a->b; begin/end are fractions of that edge. Contacts are subtracted
// in exact native edge-coordinate/Z intervals, including partial-height ones.
struct VerticalSectionPatch {
    Point a, b;
    double begin, end;
    coord_t bottom, top;
    Eigen::Vector2d outward_normal;
};
// piece must refer to an entry in topology.pieces(), not a detached copy.
std::vector<VerticalSectionPatch> exposed_vertical_section_patches(const SectionTopology& topology,
    const SectionPiece& piece, std::size_t max_patches, const std::function<void(std::size_t)>& charge);

// Numerical neighborhood primitive for a certified stream's current band.
// The caller owns connectivity/partition certification; this does not certify
// geometry. piece must be an entry in the supplied partition.
std::vector<VerticalSectionPatch> partition_vertical_section_patches(const std::vector<SectionPiece>& partition,
    const SectionPiece& piece, std::size_t max_patches, const std::function<void(std::size_t)>& charge);

// Visit in partition order, sharing an exact contact index across the band.
void visit_vertical_section_patches(const std::vector<SectionPiece>& partition, std::size_t max_patches,
    const std::function<void(const SectionPiece&, const std::vector<VerticalSectionPatch>&)>& visitor,
    const std::function<void(std::size_t)>& charge);

// Assign a boundary lying on a grid face only to its material-side cell.
// Oblique intersections use parametric coordinates, not rounded new vertices.
std::optional<SectionIntegration::Moments> vertical_patch_cell_moments(const VerticalSectionPatch& patch,
    const std::array<coord_t,3>& origin, const std::array<coord_t,3>& size);
// XY-only counterpart with identical parametric clipping/material-side rules.
std::optional<SectionIntegration::LineMoments> vertical_patch_cell_xy_moments(const VerticalSectionPatch& patch,
    const std::array<coord_t,3>& origin, const std::array<coord_t,3>& size);
} // namespace Slic3r::Thermoelastic
