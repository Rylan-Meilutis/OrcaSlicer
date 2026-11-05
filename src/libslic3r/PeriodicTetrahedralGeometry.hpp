#pragma once

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <functional>
#include <vector>

namespace Slic3r::Thermoelastic {
struct PeriodicTetrahedralGeometry
{
    std::vector<std::size_t> node_classes;
    std::vector<double> volumes_m3;
    std::size_t class_count = 0;
    double cell_volume_m3   = 0.;
    double solid_volume_m3  = 0.;
};
// Shared geometric admission for thermal and mechanical periodic cell problems.
// Exact opposing traces, no coincident nodes or degenerate/duplicate tetrahedra.
// The caller must supply a conforming, non-overlapping mesh. General element
// intersection and interior hanging-node detection are not performed here.
PeriodicTetrahedralGeometry validate_periodic_tetrahedral_geometry(const std::vector<Eigen::Vector3d>& nodes_m,
                                                                   const std::vector<std::array<std::size_t, 4>>& elements,
                                                                   const Eigen::Vector3d& minimum_m,
                                                                   const Eigen::Vector3d& maximum_m,
                                                                   std::size_t max_nodes,
                                                                   const std::function<void()>& checkpoint = {});
} // namespace Slic3r::Thermoelastic
