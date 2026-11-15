#pragma once
#include "ExPolygon.hpp"
#include <array>
#include <cstddef>
#include <functional>
#include <vector>

namespace Slic3r::Thermoelastic {
namespace Detail {
// Internal certificate shared with regression tests. Normalizes orientation;
// exact-zero triangles are harmless only when the remaining directed chain
// still matches the complete original boundary. No new vertices are allowed.
std::size_t certify_boundary_triangulation(const ExPolygon& polygon,
    std::vector<std::array<std::size_t, 3>>& triangles, std::size_t max_work,
    const std::function<void()>& checkpoint = {});
}
struct SectionBoundaryError {
    // Euclidean coordinate allowance includes rounding of the offset band.
    unsigned uncertainty_nm = 0;
    double added_area_bound_m2 = 0., removed_area_bound_m2 = 0.;
    std::size_t work = 0;
};
// Explicit native-grid representation uncertainty, never a force-error bound.
// The symmetric difference must lie in the rounded requested boundary band,
// or be certified exactly inside its reported (band+3nm) Euclidean allowance;
// for nonzero representation changes, corresponding components and holes must
// survive erosion and dilation. No area-only acceptance, topology repair,
// material overlap, or gap filling.
// Exact equality of the original occupied sets adds no representation error;
// input validation and island-survival guards still apply in that case.
// Sub-grid differences absent from rounded Boolean output require a separate
// exact whole-region band, component/hole and outward-area certificate.
SectionBoundaryError qualify_section_boundary_error(const ExPolygons& original,
    const ExPolygons& represented, unsigned band_grid_units,
    std::size_t max_points, std::size_t max_work, const std::function<void()>& checkpoint = {},
    const std::vector<const ExPolygons*>& material_regions = {});
} // namespace Slic3r::Thermoelastic
