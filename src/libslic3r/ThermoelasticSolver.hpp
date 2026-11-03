#pragma once

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <functional>
#include <optional>
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
    // Optional world-frame anisotropic constitutive matrix. When present it
    // replaces young_pa/poisson; eigenstrain remains in world coordinates.
    std::optional<Eigen::Matrix<double, 6, 6>> elasticity_pa;
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

// u(node) = affine_displacement(node) + fluctuation(node_class). Classes are
// tied exactly by elimination; only the anchor class's translation is fixed.
// This bounded RVE primitive does not create periodic geometry or penalty terms.
Result solve_affine_periodic(const std::vector<Eigen::Vector3d> &nodes_m,
    const std::vector<Element> &elements, const std::vector<std::size_t> &node_classes,
    const std::vector<Eigen::Vector3d> &affine_displacement_m, std::size_t anchor_node,
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

struct PartialBedResult {
    BedResult bonded;
    std::vector<double> released_area_m2;
    // Signed vertical displacement of original footprint nodes: positive lifts
    // away from the bed. Non-footprint entries are zero, not surface samples.
    std::vector<double> footprint_gap_m;
    double released_area_total_m2 = 0.;
    double maximum_lift_m = 0.;
    double maximum_penetration_m = 0.;
    // Traction-free released faces do not enforce unilateral bed contact.
    // Penetration beyond numerical roundoff makes that prescribed scenario
    // inadmissible, requiring a contact solve rather than a safe-looking lift.
    bool nonpenetrating = true;
};

// Prescribed contact loss, NOT predicted peeling. The caller supplies the full
// original footprint and indices of faces which remain bonded. Crack-front
// nodes shared with a bonded face remain constrained. Unrestrained components
// throw as in solve(); no artificial pins or adhesion capacity are introduced.
PartialBedResult solve_partially_bonded_bed(
    const std::vector<Eigen::Vector3d> &nodes_m, const std::vector<Element> &elements,
    const std::vector<std::array<std::size_t, 3>> &footprint_faces,
    const std::vector<std::size_t> &bonded_face_indices, double bed_z_m,
    const std::function<void()> &checkpoint = {});

struct BedContactLimits {
    unsigned max_iterations = 64;
    double gap_tolerance_m = 1e-12;
    double reaction_tolerance_n = 1e-8;
};
struct ContactBedResult {
    PartialBedResult bed;
    // Only released nodes in compressive contact; positive means upward bed
    // reaction. Other entries are zero. Small negative roundoff is retained.
    std::vector<double> normal_reaction_n;
    unsigned iterations = 0;
};
// Prescribed bond loss with frictionless, unilateral contact on the released
// footprint: gap >= 0, normal reaction >= 0, gap*reaction = 0, within supplied
// tolerances. No friction, adhesive healing or inferred fracture resistance.
// Retained faces remain fixed in all axes. A bounded active-set iteration
// rejects cycling/nonconvergence and detached rigid-body modes; never returns
// a partially converged field. Linearized contact uses the original footprint.
ContactBedResult solve_partially_bonded_contact(
    const std::vector<Eigen::Vector3d> &nodes_m, const std::vector<Element> &elements,
    const std::vector<std::array<std::size_t, 3>> &footprint_faces,
    const std::vector<std::size_t> &bonded_face_indices, double bed_z_m,
    const BedContactLimits &limits = {}, const std::function<void()> &checkpoint = {});

} // namespace Slic3r::Thermoelastic
