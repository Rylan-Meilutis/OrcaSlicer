#include "ThermalHomogenization.hpp"
#include "PeriodicTetrahedralGeometry.hpp"
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace Slic3r::ThermalHistory {
ThermalHomogenizationResult homogenize_periodic_conduction(const std::vector<Eigen::Vector3d>& nodes,
                                                           const std::vector<ConductingTetrahedron>& elements,
                                                           const Eigen::Vector3d& minimum,
                                                           const Eigen::Vector3d& maximum,
                                                           std::size_t max_nodes,
                                                           const std::function<void(unsigned)>& progress,
                                                           const std::function<void()>& checkpoint)
{
    // Bound admission before any per-element allocation.
    if (max_nodes == 0 || max_nodes > 2048 || nodes.size() > max_nodes || elements.size() > 32768)
        throw std::invalid_argument("Oversized thermal periodic cell.");
    std::vector<std::array<std::size_t, 4>> connectivity;
    connectivity.reserve(elements.size());
    for (const auto& element : elements)
        connectivity.push_back(element.nodes);
    const auto geometry = Thermoelastic::validate_periodic_tetrahedral_geometry(nodes, connectivity, minimum, maximum, max_nodes,
                                                                                checkpoint);
    using Gradient      = Eigen::Matrix<double, 3, 4>;
    std::vector<Gradient> gradients;
    std::vector<Eigen::Matrix4d> local_matrices;
    gradients.reserve(elements.size());
    local_matrices.reserve(elements.size());
    std::vector<std::size_t> parent(geometry.class_count);
    std::iota(parent.begin(), parent.end(), std::size_t(0));
    const auto root = [&parent](std::size_t id) {
        while (parent[id] != id) {
            parent[id] = parent[parent[id]];
            id         = parent[id];
        }
        return id;
    };
    ThermalHomogenizationResult result;
    Eigen::Matrix3d voigt = Eigen::Matrix3d::Zero();
    result.conductivity_w_m_k.setZero();
    result.solid_volume_fraction = geometry.solid_volume_m3 / geometry.cell_volume_m3;
    for (std::size_t i = 0; i < elements.size(); ++i) {
        if (checkpoint && i % 128 == 0)
            checkpoint();
        const auto& element = elements[i];
        const auto& k       = element.conductivity_w_m_k;
        const double scale  = k.cwiseAbs().maxCoeff();
        if (!k.allFinite() || !std::isfinite(scale) || scale <= 0. || (k - k.transpose()).cwiseAbs().maxCoeff() > 1e-12 * scale ||
            Eigen::LLT<Eigen::Matrix3d>(k / scale).info() != Eigen::Success || !std::isfinite(element.volumetric_heat_capacity_j_m3_k) ||
            element.volumetric_heat_capacity_j_m3_k <= 0.)
            throw std::invalid_argument("Invalid periodic thermal material.");
        Eigen::Matrix3d jacobian;
        for (unsigned axis = 0; axis < 3; ++axis)
            jacobian.col(axis) = nodes[element.nodes[axis + 1]] - nodes[element.nodes[0]];
        Gradient gradient;
        gradient.rightCols<3>() = jacobian.inverse().transpose();
        gradient.col(0)         = -gradient.rightCols<3>().rowwise().sum();
        gradients.push_back(gradient);
        local_matrices.push_back(geometry.volumes_m3[i] * gradient.transpose() * k * gradient);
        if (!local_matrices.back().allFinite())
            throw std::invalid_argument("Thermal element matrix overflow.");
        result.volumetric_heat_capacity_j_m3_k += geometry.volumes_m3[i] / geometry.cell_volume_m3 *
                                                  element.volumetric_heat_capacity_j_m3_k;
        voigt += geometry.volumes_m3[i] / geometry.cell_volume_m3 * k;
        for (unsigned n = 1; n < 4; ++n)
            parent[root(geometry.node_classes[element.nodes[n]])] = root(geometry.node_classes[element.nodes[0]]);
    }
    if (!std::isfinite(result.volumetric_heat_capacity_j_m3_k))
        throw std::invalid_argument("Thermal capacity overflow.");
    const std::size_t missing = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> unknown(geometry.class_count, missing);
    std::size_t count = 0;
    for (std::size_t i = 0; i < geometry.class_count; ++i)
        if (root(i) != i)
            unknown[i] = count++;
    using Sparse = Eigen::SparseMatrix<double>;
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(16 * elements.size());
    for (std::size_t i = 0; i < elements.size(); ++i) {
        if (checkpoint && i % 128 == 0)
            checkpoint();
        for (unsigned a = 0; a < 4; ++a)
            for (unsigned b = 0; b < 4; ++b) {
                const auto row = unknown[geometry.node_classes[elements[i].nodes[a]]];
                const auto col = unknown[geometry.node_classes[elements[i].nodes[b]]];
                if (row != missing && col != missing)
                    triplets.emplace_back(row, col, local_matrices[i](a, b));
            }
    }
    Sparse matrix(count, count);
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    Eigen::SimplicialLDLT<Sparse> factor;
    if (count) {
        factor.compute(matrix);
        if (factor.info() != Eigen::Success || (factor.vectorD().array() <= 0.).any())
            throw std::runtime_error("Thermal periodic cell is singular.");
    }
    for (unsigned load = 0; load < 3; ++load) {
        if (checkpoint)
            checkpoint();
        if (progress)
            progress(load);
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(count);
        for (std::size_t i = 0; i < elements.size(); ++i) {
            if (checkpoint && i % 128 == 0)
                checkpoint();
            Eigen::Vector4d affine;
            for (unsigned a = 0; a < 4; ++a)
                affine[a] = nodes[elements[i].nodes[a]][load] - minimum[load];
            const Eigen::Vector4d force = -local_matrices[i] * affine;
            for (unsigned a = 0; a < 4; ++a) {
                const auto row = unknown[geometry.node_classes[elements[i].nodes[a]]];
                if (row != missing)
                    rhs[row] += force[a];
            }
        }
        Eigen::VectorXd fluctuation = Eigen::VectorXd::Zero(count);
        if (count) {
            fluctuation = factor.solve(rhs);
            if (factor.info() != Eigen::Success || !fluctuation.allFinite() ||
                (matrix * fluctuation - rhs).norm() > 1e-8 * std::max(rhs.norm(), std::numeric_limits<double>::min()))
                throw std::runtime_error("Thermal periodic solve did not converge.");
        }
        Eigen::Vector3d flux = Eigen::Vector3d::Zero();
        double energy        = 0.;
        for (std::size_t i = 0; i < elements.size(); ++i) {
            if (checkpoint && i % 128 == 0)
                checkpoint();
            Eigen::Vector4d temperature;
            for (unsigned a = 0; a < 4; ++a) {
                const auto node = elements[i].nodes[a], row = unknown[geometry.node_classes[node]];
                temperature[a] = nodes[node][load] - minimum[load] + (row == missing ? 0. : fluctuation[row]);
            }
            const Eigen::Vector3d gradient = gradients[i] * temperature;
            const Eigen::Vector3d response = elements[i].conductivity_w_m_k * gradient;
            flux += geometry.volumes_m3[i] / geometry.cell_volume_m3 * response;
            energy += geometry.volumes_m3[i] * gradient.dot(response);
        }
        result.conductivity_w_m_k.col(load) = flux;
        const double macro                  = geometry.cell_volume_m3 * flux[load];
        // Normalize by the affine (Voigt) energy bound as well: an insulated
        // island's exact transport is zero, so relative-to-zero noise is not a
        // meaningful energy-consistency test.
        const double scale = std::max({std::abs(macro), std::abs(energy), geometry.cell_volume_m3 * voigt(load, load)});
        const double error = scale == 0. ? 0. : std::abs(macro - energy) / scale;
        if (!std::isfinite(error) || error > 1e-7)
            throw std::runtime_error("Thermal periodic energy mismatch.");
        result.maximum_energy_error = std::max(result.maximum_energy_error, error);
    }
    const auto& k                  = result.conductivity_w_m_k;
    const double norm              = std::max(k.norm(), voigt.norm());
    result.relative_symmetry_error = norm == 0. ? 0. : (k - k.transpose()).norm() / norm;
    if (!k.allFinite() || result.relative_symmetry_error > 1e-7)
        throw std::runtime_error("Thermal periodic tensor is not symmetric.");
    result.conductivity_w_m_k = .5 * (k + k.transpose()).eval();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(result.conductivity_w_m_k);
    if (eig.info() != Eigen::Success || eig.eigenvalues().minCoeff() < -1e-10 * norm)
        throw std::runtime_error("Thermal periodic tensor is not passive.");
    return result;
}
} // namespace Slic3r::ThermalHistory
