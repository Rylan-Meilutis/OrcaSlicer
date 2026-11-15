#pragma once

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <functional>
#include <vector>

namespace Slic3r::ThermalHistory {
struct ConductingTetrahedron
{
    std::array<std::size_t, 4> nodes;
    Eigen::Matrix3d conductivity_w_m_k;
    double volumetric_heat_capacity_j_m3_k;
};
struct ThermalHomogenizationResult
{
    Eigen::Matrix3d conductivity_w_m_k;
    double volumetric_heat_capacity_j_m3_k = 0.; // Per TOTAL cell volume.
    double solid_volume_fraction           = 0.;
    // Errors normalized by the affine Voigt bound, including zero-transport islands.
    double relative_symmetry_error = 0.;
    double maximum_energy_error    = 0.;
};

// Explicit solid tetrahedra with insulating voids, perfect shared-face thermal
// contact and matched periodic boundary traces. Three affine temperature-gradient
// solves. No air conductivity, interface resistance or material data is inferred.
// Requires a conforming, non-overlapping tetrahedral mesh; admission checks do
// not constitute a general mesh-intersection or hanging-node validator.
// Disconnected periodic components have independent additive-temperature gauges;
// isolated islands store heat but contribute no macroscopic steady conduction.
ThermalHomogenizationResult homogenize_periodic_conduction(const std::vector<Eigen::Vector3d>& nodes_m,
                                                           const std::vector<ConductingTetrahedron>& elements,
                                                           const Eigen::Vector3d& minimum_m,
                                                           const Eigen::Vector3d& maximum_m,
                                                           std::size_t max_nodes                         = 1024,
                                                           const std::function<void(unsigned)>& progress = {},
                                                           const std::function<void()>& checkpoint       = {});
} // namespace Slic3r::ThermalHistory
