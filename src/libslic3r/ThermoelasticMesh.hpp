#pragma once

#include "ThermoelasticSolver.hpp"
#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

namespace Slic3r::Thermoelastic {

class ResourceLimit : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct Voxel {
    std::array<int, 3> index;
    double young_pa;
    double poisson;
    Eigen::Matrix<double, 6, 1> eigenstrain = Eigen::Matrix<double, 6, 1>::Zero();
    std::size_t material = 0; // retained by deposition rasterization and refinement
};

struct Mesh {
    std::vector<Eigen::Vector3d> nodes_m;
    std::vector<Element> elements;
    std::vector<std::array<std::size_t, 3>> bed_faces;
};

struct DepositedSegment {
    Eigen::Vector3d start_m, end_m; // nozzle path; z is the top of the bead
    double width_m, height_m;
    std::size_t material;
};

struct DepositionGrid {
    Eigen::Vector3d origin_m;
    Eigen::Vector3d spacing_m;
    std::vector<Voxel> voxels;
};

// Top Z of each voxel's face-connected component. Separate bed-contact islands
// must not acquire another object's thermal-gradient length.
std::vector<double> component_top_z(const DepositionGrid &grid,
    const std::function<void()> &checkpoint = {});

// Rasterize deposited paths, not the solid CAD envelope: infill gaps remain
// empty. Horizontal beads use a capsule footprint and constant layer height.
// Material prototypes supply E, nu and eigenstrain; their indices are ignored.
// Conflicting materials in one cell and unresolved bead thickness are errors,
// never an implicit average or a reason to fill the whole object's interior.
DepositionGrid voxelize_deposition(const std::vector<DepositedSegment> &segments,
    const std::vector<Voxel> &materials, const Eigen::Vector3d &spacing_m,
    double bed_z_m, std::size_t max_cells, const std::function<void()> &checkpoint = {},
    std::size_t max_cell_visits = 2000000);

// Six conforming tetrahedra per occupied voxel. Missing voxels remain voids.
// Only face-connected voxels share nodes: diagonal/edge contact must not create
// artificial welds between separate printed islands. Index z=0 touches the bed.
Mesh mesh_voxels(const std::vector<Voxel> &voxels, const Eigen::Vector3d &origin_m,
                 const Eigen::Vector3d &spacing_m, std::size_t max_nodes,
                 const std::function<void()> &checkpoint = {});

struct RefinementResult {
    Mesh mesh;
    BedResult bed;
    unsigned subdivisions = 0;
    bool converged = false;
    bool resource_limited = false;
    double force_relative_change = 0.;
    double energy_relative_change = 0.;
};

// Refine the same occupied volume and prescribed strain field, rather than
// changing the geometry between comparisons. Both integrated tensile demand
// and strain energy must stabilize over two successive refinements. This checks
// discretization only: it does not validate material, thermal or adhesion data.
// Unconverged results must not be used to select automatic mitigations.
RefinementResult solve_refined_voxels(const std::vector<Voxel> &voxels,
    const Eigen::Vector3d &origin_m, const Eigen::Vector3d &spacing_m,
    std::size_t max_nodes, unsigned max_subdivisions, double relative_tolerance,
    const std::function<void(unsigned)> &progress = {},
    const std::function<void()> &checkpoint = {});

} // namespace Slic3r::Thermoelastic
