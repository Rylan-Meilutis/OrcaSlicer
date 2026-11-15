#pragma once
#include "ExPolygon.hpp"
#include "libslic3r.h"
#include <Eigen/Core>
#include <array>
#include <functional>

namespace Slic3r::Thermoelastic {
struct SectionTraceExtrema {
    double minimum, maximum;
    Eigen::Vector2d minimum_local, maximum_local;
};
// Exact extrema of a scalar bilinear Q1 field on a polygonal cut footprint,
// up to floating-point arithmetic. Includes edge-interior stationary points:
// testing polygon vertices alone is insufficient on non-axis-aligned edges.
// Corner order is (0,0), (1,0), (0,1), (1,1); returned coordinates are local.
SectionTraceExtrema section_trace_extrema(const ExPolygons& footprint,
    const std::array<coord_t,2>& origin, const std::array<coord_t,2>& size,
    const std::array<double,4>& corner_values, const std::function<void()>& checkpoint = {});
}
