#include "ThermoelasticHomogenization.hpp"
#include "PeriodicTetrahedralGeometry.hpp"
#include "ThermoelasticSolver.hpp"
#include "ThermoelasticMaterial.hpp"
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace Slic3r::Thermoelastic {

HomogenizationResult homogenize_periodic(const std::vector<Eigen::Vector3d> &nodes_m,
    const std::vector<Element> &elements, const Eigen::Vector3d &minimum_m,
    const Eigen::Vector3d &maximum_m, std::size_t max_nodes,
    const std::function<void(unsigned)> &progress, const std::function<void()> &checkpoint)
{
    if (elements.size() > 32768) throw std::invalid_argument("Oversized periodic RVE.");
    std::vector<std::array<std::size_t,4>> connectivity;
    connectivity.reserve(elements.size());
    for (const auto &element : elements) {
        if (!element.eigenstrain.isZero(0.))
            throw std::invalid_argument("Stiffness homogenization requires zero eigenstrain.");
        connectivity.push_back(element.nodes);
    }
    const auto geometry = validate_periodic_tetrahedral_geometry(nodes_m, connectivity,
        minimum_m, maximum_m, max_nodes, checkpoint);
    const auto &node_classes = geometry.node_classes;
    const auto &volumes = geometry.volumes_m3;
    const double cell_volume = geometry.cell_volume_m3, solid_volume = geometry.solid_volume_m3;

    HomogenizationResult result;
    result.stiffness_pa.setZero();
    result.solid_volume_fraction = solid_volume / cell_volume;
    constexpr double amplitude = .001;
    std::vector<Eigen::Vector3d> affine(nodes_m.size());
    for (unsigned load = 0; load < 6; ++load) {
        if (checkpoint) checkpoint();
        if (progress) progress(load);
        Eigen::Matrix3d strain = Eigen::Matrix3d::Zero();
        if (load < 3) strain(load, load) = amplitude;
        else {
            constexpr unsigned first[]{0, 1, 0}, second[]{1, 2, 2};
            strain(first[load-3], second[load-3]) = amplitude / 2.;
            strain(second[load-3], first[load-3]) = amplitude / 2.;
        }
        for (std::size_t i = 0; i < nodes_m.size(); ++i) affine[i] = strain * (nodes_m[i] - minimum_m);
        const auto equilibrium = solve_affine_periodic(nodes_m, elements, node_classes, affine, 0, checkpoint);
        Eigen::Matrix<double, 6, 1> average = Eigen::Matrix<double, 6, 1>::Zero();
        for (std::size_t i = 0; i < elements.size(); ++i) average += volumes[i] * equilibrium.stress_pa[i] / cell_volume;
        result.stiffness_pa.col(load) = average / amplitude;
        const double macro_energy = .5 * cell_volume * amplitude * average[load];
        const double scale = std::max(std::abs(macro_energy), std::abs(equilibrium.elastic_energy_j));
        const double error = scale == 0. ? 0. : std::abs(macro_energy - equilibrium.elastic_energy_j) / scale;
        if (!std::isfinite(error) || error > 1e-7)
            throw std::runtime_error("Periodic RVE fails energy consistency.");
        result.maximum_hill_mandel_error = std::max(result.maximum_hill_mandel_error, error);
    }
    const double norm = result.stiffness_pa.norm();
    result.relative_symmetry_error = norm == 0. ? 0. : (result.stiffness_pa - result.stiffness_pa.transpose()).norm() / norm;
    if (!result.stiffness_pa.allFinite() || result.relative_symmetry_error > 1e-7)
        throw std::runtime_error("Periodic effective stiffness is not symmetric.");
    result.stiffness_pa = (.5 * (result.stiffness_pa + result.stiffness_pa.transpose())).eval();
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> eigen(result.stiffness_pa);
    if (eigen.info() != Eigen::Success || eigen.eigenvalues().minCoeff() <= 1e-10 * eigen.eigenvalues().maxCoeff())
        throw std::runtime_error("Periodic effective stiffness has an unresolved mechanism.");
    return result;
}

EigenstrainHomogenizationResult homogenize_periodic_eigenstrain(
    const std::vector<Eigen::Vector3d> &source_nodes, const std::vector<Element> &source_elements,
    const Eigen::Vector3d &source_minimum, const Eigen::Vector3d &source_maximum,
    std::size_t max_nodes, const std::function<void(unsigned)> &progress,
    const std::function<void()> &checkpoint)
{
    if (max_nodes == 0 || max_nodes > 2048 || source_nodes.empty() || source_nodes.size() > max_nodes ||
        source_elements.empty() || source_elements.size() > 32768)
        throw std::invalid_argument("Oversized or empty periodic RVE.");
    // Freeze all source inputs before invoking caller callbacks: all seven loads
    // must refer to one material/geometry state even if a callback edits its source.
    const auto nodes_m = source_nodes;
    const auto elements = source_elements;
    const Eigen::Vector3d minimum_m = source_minimum, maximum_m = source_maximum;
    if (checkpoint) checkpoint();
    // The solver owns no load override; this bounded copy leaves caller loads intact.
    std::vector<Element> unloaded = elements;
    std::vector<std::array<std::size_t, 4>> connectivity;
    connectivity.reserve(elements.size());
    for (std::size_t i = 0; i < elements.size(); ++i) {
        if (checkpoint && i % 1024 == 0) checkpoint();
        if (!elements[i].eigenstrain.allFinite())
            throw std::invalid_argument("Non-finite periodic eigenstrain.");
        unloaded[i].eigenstrain.setZero();
        connectivity.push_back(elements[i].nodes);
    }
    EigenstrainHomogenizationResult result;
    result.elastic = homogenize_periodic(nodes_m, unloaded, minimum_m, maximum_m,
        max_nodes, progress, checkpoint);
    const auto geometry = validate_periodic_tetrahedral_geometry(nodes_m, connectivity,
        minimum_m, maximum_m, max_nodes, checkpoint);
    if (checkpoint) checkpoint();
    if (progress) progress(6);
    const std::vector<Eigen::Vector3d> affine(nodes_m.size(), Eigen::Vector3d::Zero());
    const auto equilibrium = solve_affine_periodic(nodes_m, elements, geometry.node_classes,
        affine, 0, checkpoint);
    result.stress_offset_pa.setZero();
    for (std::size_t i = 0; i < elements.size(); ++i) {
        if (checkpoint && i % 1024 == 0) checkpoint();
        result.stress_offset_pa += geometry.volumes_m3[i] / geometry.cell_volume_m3 * equilibrium.stress_pa[i];
    }
    result.zero_macro_strain_energy_density_j_m3 = equilibrium.elastic_energy_j / geometry.cell_volume_m3;
    const Eigen::LDLT<Eigen::Matrix<double, 6, 6>> factor(result.elastic.stiffness_pa);
    if (factor.info() != Eigen::Success || !factor.isPositive())
        throw std::runtime_error("Periodic thermal closure has unresolved stiffness.");
    result.effective_eigenstrain = factor.solve(-result.stress_offset_pa);
    const auto residual = (result.elastic.stiffness_pa * result.effective_eigenstrain + result.stress_offset_pa).eval();
    const double stress_scale = result.elastic.stiffness_pa.norm() * result.effective_eigenstrain.norm() +
        result.stress_offset_pa.norm();
    const double released_energy = -.5 * result.stress_offset_pa.dot(result.effective_eigenstrain);
    const double energy_scale = std::abs(result.zero_macro_strain_energy_density_j_m3) + std::abs(released_energy);
    result.residual_energy_density_j_m3 = result.zero_macro_strain_energy_density_j_m3 - released_energy;
    if (!result.stress_offset_pa.allFinite() || !result.effective_eigenstrain.allFinite() ||
        !residual.allFinite() || !std::isfinite(stress_scale) || !std::isfinite(energy_scale) ||
        !std::isfinite(result.residual_energy_density_j_m3) || residual.norm() > 1e-10 * stress_scale ||
        result.residual_energy_density_j_m3 < -4096. * std::numeric_limits<double>::epsilon() * energy_scale)
        throw std::runtime_error("Periodic thermal closure fails energy or equilibrium consistency.");
    // Only cancellation-level subtraction error is clamped, not a physical deficit.
    result.residual_energy_density_j_m3 = std::max(0., result.residual_energy_density_j_m3);
    if (checkpoint) checkpoint();
    return result;
}

HybridPeriodicResult solve_hybrid_periodic(const std::vector<Eigen::Vector3d>& source_nodes,
    const std::vector<Element>& source_elements, const std::vector<Constraint>& source_constraints,
    const std::vector<std::optional<std::size_t>>& source_assignments,
    const std::vector<EigenstrainHomogenizationResult>& source_closures,
    std::size_t max_nodes, const std::function<void()>& checkpoint)
{
    if (max_nodes == 0 || max_nodes > 2048 || source_nodes.empty() || source_nodes.size() > max_nodes ||
        source_elements.empty() || source_elements.size() > 32768 ||
        source_constraints.size() > 3 * source_nodes.size() || source_assignments.size() != source_elements.size() ||
        source_closures.size() > 256)
        throw std::invalid_argument("Hybrid periodic mesh exceeds bounded input limits.");
    const auto nodes = source_nodes;
    auto elements = source_elements;
    const auto constraints = source_constraints;
    const auto assignments = source_assignments;
    const auto closures = source_closures;
    if (checkpoint) checkpoint();
    for (const auto& closure : closures) {
        const auto stiffness = validated_elasticity(closure.elastic.stiffness_pa);
        const auto offset = (stiffness * closure.effective_eigenstrain).eval();
        const double expected = closure.residual_energy_density_j_m3 + .5 * closure.effective_eigenstrain.dot(offset);
        const double energy_scale = std::abs(expected) + std::abs(closure.zero_macro_strain_energy_density_j_m3);
        const double stress_scale = offset.norm() + closure.stress_offset_pa.norm();
        if (!closure.effective_eigenstrain.allFinite() || !closure.stress_offset_pa.allFinite() ||
            !offset.allFinite() || !std::isfinite(expected) || !std::isfinite(energy_scale) ||
            !std::isfinite(stress_scale) || !std::isfinite(closure.elastic.solid_volume_fraction) ||
            closure.elastic.solid_volume_fraction <= 0. || closure.elastic.solid_volume_fraction > 1. ||
            !std::isfinite(closure.residual_energy_density_j_m3) || closure.residual_energy_density_j_m3 < 0. ||
            (offset + closure.stress_offset_pa).norm() > 1e-10 * stress_scale ||
            std::abs(expected - closure.zero_macro_strain_energy_density_j_m3) > 1e-10 * energy_scale)
            throw std::invalid_argument("Inconsistent hybrid periodic closure.");
    }
    HybridPeriodicResult result;
    for (std::size_t i = 0; i < elements.size(); ++i) {
        if (checkpoint && i % 256 == 0) checkpoint();
        if (!assignments[i]) continue;
        if (*assignments[i] >= closures.size()) throw std::invalid_argument("Unknown hybrid periodic closure.");
        auto& element = elements[i];
        for (const auto node : element.nodes)
            if (node >= nodes.size() || !nodes[node].allFinite())
                throw std::invalid_argument("Invalid hybrid periodic element node.");
        Eigen::Matrix3d edges;
        for (unsigned axis = 0; axis < 3; ++axis)
            edges.col(axis) = nodes[element.nodes[axis + 1]] - nodes[element.nodes[0]];
        const double volume = std::abs(edges.determinant()) / 6.;
        if (!std::isfinite(volume) || volume <= 0.) throw std::invalid_argument("Invalid hybrid periodic volume.");
        const auto& closure = closures[*assignments[i]];
        element.elasticity_pa = closure.elastic.stiffness_pa;
        element.eigenstrain = closure.effective_eigenstrain;
        result.homogenized_volume_m3 += volume;
        result.residual_microenergy_j += volume * closure.residual_energy_density_j_m3;
    }
    if (!std::isfinite(result.homogenized_volume_m3) || !std::isfinite(result.residual_microenergy_j))
        throw std::invalid_argument("Non-finite hybrid periodic energy or volume.");
    result.equilibrium = solve(nodes, elements, constraints, checkpoint);
    result.total_energy_j = result.equilibrium.elastic_energy_j + result.residual_microenergy_j;
    if (!std::isfinite(result.total_energy_j)) throw std::runtime_error("Non-finite hybrid periodic total energy.");
    if (checkpoint) checkpoint();
    return result;
}

} // namespace Slic3r::Thermoelastic
