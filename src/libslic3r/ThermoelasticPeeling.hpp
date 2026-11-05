#pragma once
#include <Eigen/Core>

#include "ThermoelasticSolver.hpp"
#include <array>
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

namespace Slic3r::Thermoelastic {
struct ReleaseDemand {
    double released_area_m2 = 0.;
    double elastic_energy_j = 0.;
    // Average energy released per newly released area, not a local crack-tip
    // value or an assumed material/plate fracture resistance.
    double incremental_energy_release_j_m2 = 0.;
    double retained_tensile_force_n = 0.;
    double maximum_lift_m = 0.;
};

// Prescribed nested contact-loss scenarios at ONE fixed thermal/material state.
// With no applied mechanical loads, stored elastic energy is the potential
// energy. Decreasing it per unit released area gives a finite-increment driving
// demand. This is not peel propagation or a pass/fail adhesion calculation.
// Each entry lists the footprint face indices to retain, strictly fewer than
// the preceding entry (the initial state is fully bonded). Penetration,
// unrestrained islands and budget failures reject the entire study. Supplying
// contact_limits explicitly enables frictionless compressive re-contact on
// released faces; omitting it retains the traction-free release model.
std::vector<ReleaseDemand> bed_release_demands(
    const std::vector<Eigen::Vector3d> &nodes_m, const std::vector<Element> &elements,
    const std::vector<std::array<std::size_t, 3>> &footprint_faces,
    const std::vector<std::vector<std::size_t>> &retained_faces,
    double bed_z_m, const std::function<void()> &checkpoint = {},
    std::optional<BedContactLimits> contact_limits = std::nullopt);

struct BedReleasePatch {
    std::vector<std::size_t> faces;
    // Explicit interface fracture-work scenario, NOT bulk tensile strength.
    double fracture_energy_j_m2 = 0.;
};
struct BedReleaseEvolutionLimits {
    BedContactLimits contact;
    unsigned max_steps = 16, max_equilibrium_solves = 128;
    double energy_tolerance_j = 1e-12;
    double energy_relative_tolerance = 1e-6;
};
struct BedReleaseStep {
    std::size_t patch = 0;
    double newly_released_area_m2 = 0.;
    double released_elastic_energy_j = 0., fracture_work_j = 0.;
    double maximum_lift_m = 0.;
};
struct BedReleaseEvolution {
    PartialBedResult equilibrium;
    std::vector<std::size_t> retained_faces;
    std::vector<BedReleaseStep> steps;
    unsigned equilibrium_solves = 0;
};
// Conditional, fixed-temperature finite-patch energy study. From each current
// state, evaluate all remaining disjoint patches with unilateral bed contact;
// irreversibly release the patch with greatest positive (elastic energy drop
// minus fracture work). Return only when none of the supplied patches has a
// favorable increment beyond the energy tolerance. This is NOT stability
// against arbitrary cracks, cohesive initiation or a calibrated PEI prediction.
// Mesh/patch refinement and explicit interface properties remain necessary.
// Work is globally bounded, including contact iterations. Cancellation, rigid
// detached components, contact failure and exhausted budgets throw without
// returning a partial evolution as a completed result.
// For continuation on the SAME mesh/patch partition, pass the preceding
// retained_faces as initially_bonded_faces. Entirely released patches are
// skipped; partially bonded patches are rejected. Compressive re-contact does
// not restore adhesion. This transfers contact damage only, not viscoelastic
// history, deposition activation or qualification of time sampling.
BedReleaseEvolution evolve_bed_release(
    const std::vector<Eigen::Vector3d> &nodes_m, const std::vector<Element> &elements,
    const std::vector<std::array<std::size_t, 3>> &footprint_faces,
    const std::vector<std::size_t> &initially_bonded_faces,
    const std::vector<BedReleasePatch> &patches, double bed_z_m,
    const BedReleaseEvolutionLimits &limits = {}, const std::function<void()> &checkpoint = {});
} // namespace Slic3r::Thermoelastic
