#pragma once

#include "ExPolygon.hpp"
#include "Polyline.hpp"
#include "libslic3r.h"
#include <Eigen/Core>
#include <array>

namespace Slic3r::SectionIntegration {

// LumpedDiagonal retains consistent capacity, but lumps transverse gradient
// weights and Robin operators. It supports only grid-aligned diagonal
// conductivity, preserves affine conduction, and requires spatial refinement.
enum class ThermalQuadrature { Consistent, LumpedDiagonal };

// Internal quadrature for certified, disjoint polygons already clipped to a
// background cell. Coordinates are native nanometres; returned moments use SI
// volume and normalized cell coordinates. This is not a topology certificate.
using Moments = std::array<double, 27>;
using LineMoments = std::array<double, 9>;
// Change normalized coordinates to offset + scale*old_coordinate. Physical
// integration measure is unchanged (there is no Jacobian multiplier).
Moments rebase_moments(const Moments& moments, const Eigen::Vector3d& offset, const Eigen::Vector3d& scale);
Moments volume_moments(const ExPolygons& polygons, const std::array<coord_t, 3>& start,
    const std::array<coord_t, 3>& size, coord_t bottom, coord_t top);
// Caller supplies only exposed horizontal faces, subtracting covered/contact
// regions first. Holes are preserved; moments carry SI area, not volume.
Moments horizontal_surface_moments(const ExPolygons& polygons, const std::array<coord_t, 3>& start,
    const std::array<coord_t, 3>& size, coord_t z);
// Only exposed edge spans, not every contour of every slab: internal tile and
// material interfaces must never receive an external convection condition.
Moments vertical_surface_moments(const Polyline& edge, const std::array<coord_t, 3>& start,
    const std::array<coord_t, 3>& size, coord_t bottom, coord_t top);
// A single edge clipped parametrically to a cell, avoiding native-coordinate
// rounding at oblique grid intersections. XY endpoints and Z limits are in
// normalized cell coordinates; size_m supplies the physical metric.
Moments vertical_segment_moments(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
    const Eigen::Vector3d& size_m, double bottom, double top);
// Separable vertical-face integration: XY moments carry SI length and can be
// reused for different Z intervals without changing any geometric boundary.
LineMoments vertical_segment_xy_moments(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
    const Eigen::Vector2d& size_m);
Moments extrude_line_moments(const LineMoments& xy, double height_m, double bottom, double top);

struct ThermalElement {
    Eigen::Matrix<double, 8, 8> capacity = Eigen::Matrix<double, 8, 8>::Zero();
    Eigen::Matrix<double, 8, 8> conduction = Eigen::Matrix<double, 8, 8>::Zero();
};
// Consistent Q1 capacity and conduction integrated over occupied volume only.
// Node bits encode high X/Y/Z. Call separately for each physical component and
// material; sharing a background cell does not authorize bonding across voids.
ThermalElement thermal_element(const Moments& moments, const Eigen::Vector3d& size_m,
    double volumetric_capacity_j_m3_k, const Eigen::Matrix3d& conductivity_w_m_k,
    ThermalQuadrature quadrature = ThermalQuadrature::Consistent);

struct ThermalBoundaryElement {
    Eigen::Matrix<double, 8, 8> conductance = Eigen::Matrix<double, 8, 8>::Zero();
    Eigen::Matrix<double, 8, 1> load = Eigen::Matrix<double, 8, 1>::Zero();
};
// Explicit Robin condition h*(T_reservoir-T_surface); no PWM conversion.
ThermalBoundaryElement thermal_boundary_element(const Moments& surface_moments,
    double coefficient_w_m2_k, double reservoir_temperature_k,
    ThermalQuadrature quadrature = ThermalQuadrature::Consistent);

} // namespace Slic3r::SectionIntegration
