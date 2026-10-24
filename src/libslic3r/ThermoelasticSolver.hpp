#pragma once

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <functional>
#include <vector>

namespace Slic3r::Thermoelastic {

// SI throughout: coordinates/displacements in m, E/stress in Pa, reactions in N.
// Linear tetrahedral elements with prescribed, element-constant eigenstrain.
// This solves equilibrium, not a thermal history or a bed adhesion/fracture law.
struct Element {
    std::array<std::size_t, 4> nodes;
    double young_pa;
    double poisson;
    // Engineering strain order: xx, yy, zz, xy, yz, xz.
    Eigen::Matrix<double, 6, 1> eigenstrain = Eigen::Matrix<double, 6, 1>::Zero();
};

struct Constraint {
    std::size_t node;
    unsigned axis;
    double displacement_m = 0.;
};

struct Result {
    std::vector<Eigen::Vector3d> displacement_m;
    std::vector<Eigen::Vector3d> reaction_n;
    // Stress is element-wise; do not smooth across voids or material interfaces.
    std::vector<Eigen::Matrix<double, 6, 1>> stress_pa;
    double elastic_energy_j = 0.;
    double relative_residual = 0.;
};

// Constraints must remove every component's rigid-body modes. Invalid, singular
// and numerically unresolved systems throw instead of returning a safe-looking
// zero field. The optional checkpoint may throw to cancel; no partial result is
// returned. Called between assembly/recovery batches and around factorization.
Result solve(const std::vector<Eigen::Vector3d> &nodes_m,
             const std::vector<Element> &elements,
             const std::vector<Constraint> &constraints,
             const std::function<void()> &checkpoint = {});

struct BedResult {
    Result equilibrium;
    // Lumped boundary area and traction at each mesh node. A positive normal
    // demand means the bed must pull down on the part to retain contact.
    std::vector<double> contact_area_m2;
    std::vector<double> tensile_traction_pa;
    double tensile_force_n = 0.;
    double contact_area_total_m2 = 0.;
};

// Fully bonded, rigid, horizontal bed scenario. Contact triangles must be real
// exterior element faces at bed_z_m, not a bounding rectangle spanning holes.
// This computes the required restraint, NOT the available adhesion or a pass/
// fail prediction. Sharp-edge traction is mesh-dependent and must be checked
// for convergence; integrated demand is retained separately from nodal peaks.
BedResult solve_bonded_bed(const std::vector<Eigen::Vector3d> &nodes_m,
                          const std::vector<Element> &elements,
                          const std::vector<std::array<std::size_t, 3>> &contact_faces,
                          double bed_z_m,
                          const std::function<void()> &checkpoint = {});

} // namespace Slic3r::Thermoelastic
