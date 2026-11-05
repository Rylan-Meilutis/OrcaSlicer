#pragma once

#include "ThermoelasticSolver.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

namespace Slic3r::Thermoelastic {

struct HomogenizationResult {
    // Engineering strain order xx,yy,zz,xy,yz,xz; stress shear is not doubled.
    Eigen::Matrix<double, 6, 6> stiffness_pa;
    double solid_volume_fraction = 0.;
    double relative_symmetry_error = 0.;
    double maximum_hill_mandel_error = 0.;
};

// Conforming axis-aligned periodic RVE, including explicit voids. Opposite
// solid boundary triangulations must match exactly in their tangential
// coordinates. Duplicate/coincident nodes are rejected rather than welded.
// Six affine-periodic engineering-strain solves; zero input eigenstrain only.
// No density-scaled modulus, thermal coefficient or weld property is inferred.
HomogenizationResult homogenize_periodic(const std::vector<Eigen::Vector3d> &nodes_m,
    const std::vector<Element> &elements, const Eigen::Vector3d &minimum_m,
    const Eigen::Vector3d &maximum_m, std::size_t max_nodes = 1024,
    const std::function<void(unsigned)> &progress = {},
    const std::function<void()> &checkpoint = {});

struct EigenstrainHomogenizationResult {
    HomogenizationResult elastic;
    Eigen::Matrix<double, 6, 1> stress_offset_pa;
    Eigen::Matrix<double, 6, 1> effective_eigenstrain;
    double zero_macro_strain_energy_density_j_m3 = 0.;
    double residual_energy_density_j_m3 = 0.;
};
// Same periodic geometry assumptions, with explicit per-element eigenstrain.
// For engineering macro strain E, average stress is C*E + stress_offset_pa;
// energy per TOTAL cell volume is .5*E.dot(C*E) + E.dot(stress_offset_pa)
// + zero_macro_strain_energy_density_j_m3. Residual energy at the effective
// free strain is retained: a heterogeneous body need not become stress-free.
// No temperature law, interface calibration or macroscopic applicability is inferred.
// Progress indices 0..5 are stiffness loads; 6 is the eigenstrain load.
// Bounded inputs are snapshotted before any callback to keep these loads coherent.
EigenstrainHomogenizationResult homogenize_periodic_eigenstrain(
    const std::vector<Eigen::Vector3d> &nodes_m, const std::vector<Element> &elements,
    const Eigen::Vector3d &minimum_m, const Eigen::Vector3d &maximum_m,
    std::size_t max_nodes = 1024, const std::function<void(unsigned)> &progress = {},
    const std::function<void()> &checkpoint = {});

struct HybridPeriodicResult {
    // Macro elastic field; its elastic_energy_j excludes the constant microenergy.
    Result equilibrium;
    double residual_microenergy_j = 0.;
    double total_energy_j = 0.;
    double homogenized_volume_m3 = 0.;
};
// Explicit conforming macro mesh, not an automatic repeated-region classifier.
// Each optional entry assigns a closure to that element, replacing its C and
// eigenstrain. Unassigned elements and ALL connectivity/restraints are preserved.
// Caller must justify closure applicability and interfaces; shared nodes imply
// perfect bonding. No void, bed contact or adhesion certificate is inferred.
// All inputs are bounded and frozen before callbacks; no partial field returns.
HybridPeriodicResult solve_hybrid_periodic(const std::vector<Eigen::Vector3d>& nodes_m,
    const std::vector<Element>& elements, const std::vector<Constraint>& constraints,
    const std::vector<std::optional<std::size_t>>& element_closures,
    const std::vector<EigenstrainHomogenizationResult>& closures,
    std::size_t max_nodes = 2048, const std::function<void()>& checkpoint = {});

} // namespace Slic3r::Thermoelastic
