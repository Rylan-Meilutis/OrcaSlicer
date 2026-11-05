#include "ThermoelasticPeeling.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticSolver.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
std::vector<ReleaseDemand> bed_release_demands(
    const std::vector<Eigen::Vector3d> &nodes_m, const std::vector<Element> &elements,
    const std::vector<std::array<std::size_t, 3>> &footprint_faces,
    const std::vector<std::vector<std::size_t>> &retained_faces,
    double bed_z_m, const std::function<void()> &checkpoint,
    std::optional<BedContactLimits> contact_limits)
{
    if (nodes_m.size() > 12000 || elements.size() > 100000 || footprint_faces.size() > 400000 || retained_faces.size() > 16)
        throw ResourceLimit("Bed release study exceeds its resource budget.");
    if (retained_faces.empty()) throw std::invalid_argument("Bed release study has no release steps.");
    std::vector<bool> previous(footprint_faces.size(), true);
    for (const auto &faces : retained_faces) {
        if (checkpoint) checkpoint();
        if (faces.empty() || faces.size() >= static_cast<std::size_t>(std::count(previous.begin(), previous.end(), true)))
            throw std::invalid_argument("Bed release steps must strictly reduce retained contact.");
        std::vector<bool> next(footprint_faces.size(), false);
        for (const auto face : faces) {
            if (face >= previous.size() || !previous[face] || next[face])
                throw std::invalid_argument("Bed release steps must be nested distinct footprint subsets.");
            next[face] = true;
        }
        previous.swap(next);
    }
    std::vector<std::size_t> all_faces(footprint_faces.size());
    std::iota(all_faces.begin(), all_faces.end(), std::size_t(0));
    const auto solve_state = [&](const std::vector<std::size_t> &faces) {
        if (contact_limits)
            return solve_partially_bonded_contact(nodes_m, elements, footprint_faces, faces,
                bed_z_m, *contact_limits, checkpoint).bed;
        return solve_partially_bonded_bed(nodes_m, elements, footprint_faces, faces, bed_z_m, checkpoint);
    };
    const auto initial = solve_state(all_faces);
    double energy = initial.bonded.equilibrium.elastic_energy_j;
    double area = 0.;
    std::vector<ReleaseDemand> demands;
    demands.reserve(retained_faces.size());
    for (const auto &faces : retained_faces) {
        if (checkpoint) checkpoint();
        const auto released = solve_state(faces);
        if (!released.nonpenetrating)
            throw std::runtime_error("Released surface penetrates the bed; contact analysis is unresolved.");
        const double next_energy = released.bonded.equilibrium.elastic_energy_j;
        const double delta_area = released.released_area_total_m2 - area;
        const double delta_energy = energy - next_energy;
        // Relaxing displacement restraints cannot increase the minimum stored
        // energy. Do not silently clamp a materially inconsistent solve to zero.
        const double tolerance = 1e-8 * std::max(energy, next_energy);
        if (!(delta_area > 0.) || delta_energy < -tolerance)
            throw std::runtime_error("Bed release energy is numerically inconsistent.");
        const double demand = std::max(0., delta_energy) / delta_area;
        if (!std::isfinite(demand)) throw std::runtime_error("Bed release demand overflow.");
        demands.push_back({released.released_area_total_m2, next_energy, demand,
                           released.bonded.tensile_force_n, released.maximum_lift_m});
        area = released.released_area_total_m2;
        energy = next_energy;
    }
    return demands;
}

BedReleaseEvolution evolve_bed_release(
    const std::vector<Eigen::Vector3d> &nodes_m, const std::vector<Element> &elements,
    const std::vector<std::array<std::size_t, 3>> &footprint_faces,
    const std::vector<std::size_t> &initially_bonded_faces,
    const std::vector<BedReleasePatch> &patches, double bed_z_m,
    const BedReleaseEvolutionLimits &limits, const std::function<void()> &checkpoint)
{
    if (nodes_m.size() > 12000 || elements.size() > 100000 || footprint_faces.size() > 400000 || patches.size() > 32)
        throw ResourceLimit("Bed release evolution exceeds its geometry budget.");
    if (patches.empty() || limits.max_steps == 0 || limits.max_steps > 32 || limits.max_equilibrium_solves == 0 ||
        limits.max_equilibrium_solves > 256 || limits.contact.max_iterations == 0 || limits.contact.max_iterations > 128 ||
        !std::isfinite(limits.energy_tolerance_j) || limits.energy_tolerance_j <= 0. ||
        !std::isfinite(limits.energy_relative_tolerance) || limits.energy_relative_tolerance <= 0. ||
        limits.energy_relative_tolerance >= 1.)
        throw std::invalid_argument("Invalid bed release evolution limits.");
    std::vector<bool> initially_bonded(footprint_faces.size(), false);
    for (const auto face : initially_bonded_faces) {
        if (face >= footprint_faces.size() || initially_bonded[face])
            throw std::invalid_argument("Invalid initial bonded footprint.");
        initially_bonded[face] = true;
    }
    std::vector<std::size_t> owner(footprint_faces.size(), patches.size());
    std::vector<bool> released(patches.size(), false);
    for (std::size_t p = 0; p < patches.size(); ++p) {
        if (checkpoint) checkpoint();
        const auto &patch = patches[p];
        if (patch.faces.empty() || !std::isfinite(patch.fracture_energy_j_m2) || patch.fracture_energy_j_m2 < 0.)
            throw std::invalid_argument("Invalid bed release patch or fracture energy.");
        std::size_t bonded_count = 0;
        for (const auto face : patch.faces) {
            if (face >= owner.size() || owner[face] != patches.size())
                throw std::invalid_argument("Release patches must be disjoint footprint subsets.");
            bonded_count += initially_bonded[face] ? 1 : 0;
            owner[face] = p;
        }
        // Carry the same patch identities across thermal snapshots without
        // re-bonding released faces or charging their fracture work twice.
        // A changed partition needs an explicit history transfer, not an
        // implicit interpretation of a partially released patch.
        if (bonded_count != 0 && bonded_count != patch.faces.size())
            throw std::invalid_argument("A release patch cannot be partially bonded at initialization.");
        released[p] = bonded_count == 0;
    }
    BedReleaseEvolution result;
    const auto solve_state = [&](const std::vector<std::size_t> &faces) {
        if (checkpoint) checkpoint();
        if (result.equilibrium_solves >= limits.max_equilibrium_solves)
            throw ResourceLimit("Bed release equilibrium budget exhausted.");
        auto contact_limits = limits.contact;
        contact_limits.max_iterations = std::min(contact_limits.max_iterations,
            limits.max_equilibrium_solves - result.equilibrium_solves);
        auto state = solve_partially_bonded_contact(nodes_m, elements, footprint_faces,
            faces, bed_z_m, contact_limits, checkpoint);
        result.equilibrium_solves += state.iterations;
        return std::move(state.bed);
    };
    result.retained_faces = initially_bonded_faces;
    result.equilibrium = solve_state(result.retained_faces);
    for (;;) {
        std::optional<PartialBedResult> best;
        std::vector<std::size_t> best_faces;
        BedReleaseStep best_step;
        double best_excess = 0.;
        const double energy = result.equilibrium.bonded.equilibrium.elastic_energy_j;
        for (std::size_t p = 0; p < patches.size(); ++p) {
            if (released[p]) continue;
            std::vector<std::size_t> retained;
            retained.reserve(result.retained_faces.size());
            for (const auto face : result.retained_faces)
                if (owner[face] != p) retained.push_back(face);
            auto trial = solve_state(retained);
            const double area = trial.released_area_total_m2 - result.equilibrium.released_area_total_m2;
            const double next_energy = trial.bonded.equilibrium.elastic_energy_j;
            const double energy_drop = energy - next_energy;
            const double work = patches[p].fracture_energy_j_m2 * area;
            if (!(area > 0.) || !std::isfinite(energy_drop) || !std::isfinite(work) ||
                // An affine, stress-free contraction can have nonzero
                // displacement but only roundoff-sized stored energy. Use the
                // explicit absolute energy tolerance there, not a relative
                // comparison against an almost-zero denominator.
                energy_drop < -(limits.energy_tolerance_j + 1e-8 * std::max(energy, next_energy)))
                throw std::runtime_error("Bed release trial energy is inconsistent.");
            const double excess = energy_drop - work;
            const double margin = limits.energy_tolerance_j + limits.energy_relative_tolerance *
                std::max({energy, next_energy, work});
            if (excess > margin && excess > best_excess) {
                best_excess = excess;
                best_step = {p, area, energy_drop, work, trial.maximum_lift_m};
                best = std::move(trial);
                best_faces = std::move(retained);
            }
        }
        if (!best) return result;
        if (result.steps.size() >= limits.max_steps)
            throw ResourceLimit("Bed release step budget exhausted before the patch study completed.");
        released[best_step.patch] = true;
        result.steps.push_back(best_step);
        result.retained_faces = std::move(best_faces);
        result.equilibrium = std::move(*best);
    }
}
} // namespace Slic3r::Thermoelastic
