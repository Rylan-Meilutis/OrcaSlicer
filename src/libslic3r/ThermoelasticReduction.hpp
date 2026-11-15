#pragma once

#include "ThermoelasticSolver.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

namespace Slic3r::Thermoelastic {
class ReductionBudget : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};
struct InterpolationWeight
{
    std::size_t node;
    double weight;
};
struct ReductionLimits
{
    std::size_t max_fine_nodes       = 100000;
    std::size_t max_elements         = 200000;
    std::size_t max_coarse_nodes     = 12000;
    std::size_t max_free_dofs        = 6000;
    std::size_t max_assembly_entries = 5000000;
    std::size_t max_assembly_work    = 20000000;
};

struct ReducedResult
{
    // Displacement/stress and constrained-node forces recovered on fine geometry.
    // relative_residual measures P^T*r, NOT the full fine residual. The recovered
    // forces are not qualified bed-contact demands until reduction convergence
    // against a richer basis/fine reference has been demonstrated.
    Result recovered;
    // P^T*r on constrained coarse DOFs, including contributions from fine nodes
    // that are not themselves constrained. Not a local bed-traction distribution.
    std::vector<Eigen::Vector3d> generalized_reaction_n;
    double fine_free_residual_norm_n = 0.;
    std::size_t free_reduced_dofs    = 0;
};

struct LineBasis
{
    std::vector<Eigen::Vector3d> nodes_m;
    std::vector<std::vector<InterpolationWeight>> rows;
};
// Interpolate only along uninterrupted axis-aligned tetrahedral edges. Distinct
// node IDs are never welded, even at coincident coordinates. Pinned nodes remain
// unit rows. The caller must provide a conforming occupied-grid mesh.
LineBasis occupied_line_basis(const std::vector<Eigen::Vector3d>& nodes_m,
                              const std::vector<Element>& elements,
                              const std::vector<std::size_t>& pinned_nodes,
                              unsigned axis,
                              unsigned stride,
                              const ReductionLimits& limits           = {},
                              const std::function<void()>& checkpoint = {});

// Caller-supplied kinematic subspace u_fine=P*u_coarse. Geometry, material,
// eigenstrain and integration remain on the ORIGINAL fine tetrahedra. No
// homogenization, geometry merging or automatic coarsening occurs here.
// Each row has 1..8 positive weights and reproduces constant and coordinate
// fields. Each basis has connected positive support within one fine component. Each
// constrained fine node must be an exact unit row, so constraints are eliminated
// exactly rather than approximately imposed. Unused or rank-deficient bases
// are rejected. The reported equilibrium residual is in the reduced subspace;
// fine stress/reactions are recovered without asserting fine equilibrium.
ReducedResult solve_reduced(const std::vector<Eigen::Vector3d>& fine_nodes_m,
                            const std::vector<Element>& elements,
                            const std::vector<Constraint>& constraints,
                            const std::vector<Eigen::Vector3d>& coarse_nodes_m,
                            const std::vector<std::vector<InterpolationWeight>>& prolongation,
                            const ReductionLimits& limits           = {},
                            const std::function<void()>& checkpoint = {});
} // namespace Slic3r::Thermoelastic
