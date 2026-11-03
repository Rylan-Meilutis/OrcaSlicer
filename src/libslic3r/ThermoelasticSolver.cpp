#include "ThermoelasticSolver.hpp"
#include "ThermoelasticMaterial.hpp"
#include "ThermoelasticReduction.hpp"

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
#include <utility>

namespace Slic3r::Thermoelastic {
namespace {
using StrainMatrix = Eigen::Matrix<double, 6, 12>;
using Elasticity = Eigen::Matrix<double, 6, 6>;

struct ElementMatrices {
    StrainMatrix b = StrainMatrix::Zero();
    Elasticity d = Elasticity::Zero();
    double volume;
};

ElementMatrices matrices(const std::vector<Eigen::Vector3d> &nodes, const Element &e)
{
    if ((!e.elasticity_pa && (!std::isfinite(e.young_pa) || e.young_pa <= 0. ||
        !std::isfinite(e.poisson) || e.poisson <= -1. || e.poisson >= 0.5)) || !e.eigenstrain.allFinite())
        throw std::invalid_argument("Invalid thermoelastic material or eigenstrain.");
    for (std::size_t n : e.nodes)
        if (n >= nodes.size()) throw std::invalid_argument("Invalid tetrahedron node index.");
    Eigen::Matrix3d jacobian;
    for (int i = 0; i < 3; ++i)
        jacobian.col(i) = nodes[e.nodes[i + 1]] - nodes[e.nodes[0]];
    const double determinant = jacobian.determinant();
    const double scale = jacobian.cwiseAbs().maxCoeff();
    if (scale == 0. || !std::isfinite(determinant) ||
        std::abs(determinant) <= 1e-12 * scale * scale * scale)
        throw std::invalid_argument("Degenerate tetrahedron in thermoelastic mesh.");
    ElementMatrices result;
    result.volume = std::abs(determinant) / 6.;
    // Shape-function gradients in physical coordinates. Orientation is allowed
    // to be either sign; volume is positive and gradients retain their signs.
    Eigen::Matrix<double, 4, 3> gradients;
    gradients.bottomRows<3>() = jacobian.inverse();
    gradients.row(0) = -gradients.bottomRows<3>().colwise().sum();
    for (int i = 0; i < 4; ++i) {
        const double x = gradients(i, 0), y = gradients(i, 1), z = gradients(i, 2);
        result.b(0, 3*i) = x;
        result.b(1, 3*i+1) = y;
        result.b(2, 3*i+2) = z;
        result.b(3, 3*i) = y; result.b(3, 3*i+1) = x;
        result.b(4, 3*i+1) = z; result.b(4, 3*i+2) = y;
        result.b(5, 3*i) = z; result.b(5, 3*i+2) = x;
    }
    if (e.elasticity_pa) result.d = validated_elasticity(*e.elasticity_pa);
    else {
        const double mu = e.young_pa / (2. * (1. + e.poisson));
        const double lambda = e.young_pa * e.poisson / ((1. + e.poisson) * (1. - 2. * e.poisson));
        result.d.topLeftCorner<3, 3>().setConstant(lambda);
        result.d.diagonal().head<3>().array() += 2. * mu;
        result.d.diagonal().tail<3>().setConstant(mu);
    }
    if (!result.d.allFinite() || !result.b.allFinite())
        throw std::invalid_argument("Non-finite thermoelastic element coefficients.");
    return result;
}
} // namespace

namespace {
Result solve_mapped(const std::vector<Eigen::Vector3d> &nodes_m, const std::vector<Element> &elements,
    const std::vector<Constraint> &constraints, const std::vector<std::size_t> &node_classes,
    const std::vector<Eigen::Vector3d> &offsets, const std::function<void()> &checkpoint)
{
    if (nodes_m.empty() || elements.empty() || nodes_m.size() > std::size_t(std::numeric_limits<int>::max() / 3))
        throw std::invalid_argument("Empty or oversized thermoelastic mesh.");
    for (const auto &node : nodes_m)
        if (!node.allFinite()) throw std::invalid_argument("Non-finite thermoelastic node.");
    const int dofs = int(nodes_m.size() * 3);
    Eigen::VectorXd displacement = Eigen::VectorXd::Zero(dofs);
    if (!offsets.empty()) {
        if (offsets.size() != nodes_m.size()) throw std::invalid_argument("Invalid affine offset count.");
        for (std::size_t i = 0; i < offsets.size(); ++i) {
            if (!offsets[i].allFinite()) throw std::invalid_argument("Non-finite affine offset.");
            displacement.segment<3>(3*i) = offsets[i];
        }
    }
    std::vector<int> representatives(dofs);
    if (!node_classes.empty() && node_classes.size() != nodes_m.size())
        throw std::invalid_argument("Invalid periodic node class count.");
    for (std::size_t i = 0; i < nodes_m.size(); ++i) {
        const auto representative = node_classes.empty() ? i : node_classes[i];
        if (representative >= nodes_m.size()) throw std::invalid_argument("Invalid periodic node class.");
        for (unsigned axis = 0; axis < 3; ++axis) representatives[3*i+axis] = int(3*representative+axis);
    }
    Eigen::VectorXd prescribed = Eigen::VectorXd::Zero(dofs);
    std::vector<bool> fixed(dofs, false);
    for (const auto &c : constraints) {
        if (c.node >= nodes_m.size() || c.axis > 2 || !std::isfinite(c.displacement_m))
            throw std::invalid_argument("Invalid thermoelastic constraint.");
        const int original = int(3 * c.node + c.axis);
        const int i = representatives[original];
        const double value = c.displacement_m - displacement[original];
        if (fixed[i] && prescribed[i] != value)
            throw std::invalid_argument("Conflicting thermoelastic constraints.");
        fixed[i] = true;
        prescribed[i] = value;
    }
    std::vector<int> free_index(dofs, -1);
    int free_count = 0;
    std::vector<int> class_free_index(dofs, -1);
    std::vector<bool> original_fixed(dofs, false);
    for (int i = 0; i < dofs; ++i) {
        const int representative = representatives[i];
        original_fixed[i] = fixed[representative];
        if (fixed[representative]) displacement[i] += prescribed[representative];
        else {
            if (class_free_index[representative] < 0) class_free_index[representative] = free_count++;
            free_index[i] = class_free_index[representative];
        }
    }
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(free_count);
    std::vector<Eigen::Triplet<double>> entries;
    // Assemble only the reduced system: prescribed displacements are eliminated
    // exactly, not enforced using arbitrary penalty springs.
    for (std::size_t n = 0; n < elements.size(); ++n) {
        if (checkpoint && n % 128 == 0) checkpoint();
        const auto &e = elements[n];
        const auto m = matrices(nodes_m, e);
        const Eigen::Matrix<double, 12, 12> k = m.volume * m.b.transpose() * m.d * m.b;
        const Eigen::Matrix<double, 12, 1> f = m.volume * m.b.transpose() * m.d * e.eigenstrain;
        if (!k.allFinite() || !f.allFinite())
            throw std::runtime_error("Thermoelastic element assembly overflowed.");
        for (int i = 0; i < 12; ++i) {
            const int row = int(3 * e.nodes[i / 3] + i % 3);
            if (original_fixed[row]) continue;
            rhs[free_index[row]] += f[i];
            for (int j = 0; j < 12; ++j) {
                const int col = int(3 * e.nodes[j / 3] + j % 3);
                rhs[free_index[row]] -= k(i, j) * displacement[col];
                if (!original_fixed[col] && k(i, j) != 0.) entries.emplace_back(free_index[row], free_index[col], k(i, j));
            }
        }
    }
    Result result;
    if (free_count != 0) {
        Eigen::SparseMatrix<double> k(free_count, free_count);
        k.setFromTriplets(entries.begin(), entries.end());
        if (checkpoint) checkpoint();
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(k);
        if (checkpoint) checkpoint();
        if (solver.info() != Eigen::Success || !solver.vectorD().allFinite() ||
            solver.vectorD().minCoeff() <= 1e-12 * solver.vectorD().cwiseAbs().maxCoeff())
            throw std::runtime_error("Thermoelastic system is singular or ill-conditioned; check component restraints and mesh quality.");
        const Eigen::VectorXd u = solver.solve(rhs);
        if (solver.info() != Eigen::Success || !u.allFinite())
            throw std::runtime_error("Thermoelastic solve failed.");
        const double scale = rhs.norm() + (k * u).norm();
        result.relative_residual = scale > 0. ? (k * u - rhs).norm() / scale : 0.;
        if (!std::isfinite(result.relative_residual) || result.relative_residual > 1e-8)
            throw std::runtime_error("Thermoelastic equilibrium residual exceeds tolerance.");
        for (int i = 0; i < dofs; ++i)
            if (!original_fixed[i]) displacement[i] += u[free_index[i]];
    }
    result.displacement_m.resize(nodes_m.size());
    result.reaction_n.assign(nodes_m.size(), Eigen::Vector3d::Zero());
    result.stress_pa.reserve(elements.size());
    for (std::size_t i = 0; i < nodes_m.size(); ++i)
        result.displacement_m[i] = displacement.segment<3>(3 * i);
    for (std::size_t n = 0; n < elements.size(); ++n) {
        if (checkpoint && n % 128 == 0) checkpoint();
        const auto &e = elements[n];
        const auto m = matrices(nodes_m, e);
        Eigen::Matrix<double, 12, 1> u;
        for (int i = 0; i < 4; ++i) u.segment<3>(3*i) = result.displacement_m[e.nodes[i]];
        const Eigen::Matrix<double, 6, 1> elastic_strain = m.b * u - e.eigenstrain;
        const Eigen::Matrix<double, 6, 1> stress = m.d * elastic_strain;
        if (!stress.allFinite()) throw std::runtime_error("Non-finite thermoelastic stress.");
        result.stress_pa.push_back(stress);
        result.elastic_energy_j += 0.5 * m.volume * elastic_strain.dot(stress);
        const Eigen::Matrix<double, 12, 1> force = m.volume * m.b.transpose() * stress;
        for (int i = 0; i < 12; ++i)
            if (original_fixed[3 * e.nodes[i/3] + i%3] || !node_classes.empty()) result.reaction_n[e.nodes[i/3]][i%3] += force[i];
    }
    if (!std::isfinite(result.elastic_energy_j))
        throw std::runtime_error("Non-finite thermoelastic energy.");
    if (checkpoint) checkpoint();
    return result;
}
} // namespace

Result solve(const std::vector<Eigen::Vector3d> &nodes_m, const std::vector<Element> &elements,
             const std::vector<Constraint> &constraints, const std::function<void()> &checkpoint)
{
    return solve_mapped(nodes_m, elements, constraints, {}, {}, checkpoint);
}

ReducedResult solve_reduced(const std::vector<Eigen::Vector3d>& nodes,
                            const std::vector<Element>& elements,
                            const std::vector<Constraint>& constraints,
                            const std::vector<Eigen::Vector3d>& coarse,
                            const std::vector<std::vector<InterpolationWeight>>& interpolation,
                            const ReductionLimits& limits,
                            const std::function<void()>& checkpoint)
{
    if (nodes.empty() || elements.empty() || coarse.empty() || interpolation.size() != nodes.size() || limits.max_fine_nodes == 0 ||
        limits.max_fine_nodes > 100000 || constraints.size() > 3 * nodes.size() || limits.max_elements == 0 ||
        limits.max_elements > 200000 || limits.max_coarse_nodes == 0 ||
        limits.max_coarse_nodes > 12000 || limits.max_free_dofs == 0 ||
        limits.max_free_dofs > 6000 || limits.max_assembly_entries == 0 || limits.max_assembly_entries > 5000000 ||
        limits.max_assembly_work == 0 || limits.max_assembly_work > 20000000)
        throw std::invalid_argument("Invalid or oversized mechanical reduction.");
    if (nodes.size() > limits.max_fine_nodes || elements.size() > limits.max_elements || coarse.size() > limits.max_coarse_nodes)
        throw ReductionBudget("Mechanical reduction input budget exceeded.");
    const auto missing = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> parent(nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        parent[i] = i;
        if (!nodes[i].allFinite())
            throw std::invalid_argument("Invalid reduction coordinates.");
    }
    const auto root = [&parent](std::size_t node) {
        while (parent[node] != node) {
            parent[node] = parent[parent[node]];
            node         = parent[node];
        }
        return node;
    };
    std::vector<bool> used(nodes.size(), false);
    for (std::size_t i = 0; i < elements.size(); ++i) {
        if (checkpoint && i % 128 == 0)
            checkpoint();
        for (auto node : elements[i].nodes) {
            if (node >= nodes.size())
                throw std::invalid_argument("Invalid reduction element index.");
            used[node] = true;
        }
        for (unsigned j = 1; j < 4; ++j)
            parent[root(elements[i].nodes[j])] = root(elements[i].nodes[0]);
    }
    if (std::find(used.begin(), used.end(), false) != used.end())
        throw std::invalid_argument("Unused fine reduction node.");
    double length = 0.;
    for (const auto& node : nodes)
        length = std::max(length, (node - nodes.front()).norm());
    if (!std::isfinite(length) || length <= 0.)
        throw std::invalid_argument("Invalid reduction extent.");
    for (const auto& node : coarse)
        if (!node.allFinite())
            throw std::invalid_argument("Invalid coarse reduction coordinates.");
    std::vector<std::size_t> component(coarse.size(), missing);
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        if (checkpoint && i % 128 == 0)
            checkpoint();
        const auto& row = interpolation[i];
        if (row.empty() || row.size() > 8)
            throw std::invalid_argument("Invalid prolongation row width.");
        double sum                 = 0.;
        Eigen::Vector3d reproduced = Eigen::Vector3d::Zero();
        for (std::size_t j = 0; j < row.size(); ++j) {
            const auto& entry = row[j];
            if (entry.node >= coarse.size() || !std::isfinite(entry.weight) || entry.weight <= 0. || entry.weight > 1.)
                throw std::invalid_argument("Invalid interpolation weight.");
            for (std::size_t k = 0; k < j; ++k)
                if (row[k].node == entry.node)
                    throw std::invalid_argument("Duplicate interpolation basis.");
            if (component[entry.node] != missing && component[entry.node] != root(i))
                throw std::invalid_argument("Reduction basis spans disconnected components.");
            component[entry.node] = root(i);
            sum += entry.weight;
            reproduced += entry.weight * (coarse[entry.node] - nodes.front());
        }
        if (std::abs(sum - 1.) > 1e-12 || (reproduced - (nodes[i] - nodes.front())).norm() > 1e-10 * length)
            throw std::invalid_argument("Prolongation does not reproduce affine displacement.");
    }
    if (std::find(component.begin(), component.end(), missing) != component.end())
        throw std::invalid_argument("Unused coarse reduction node.");
    // A globally connected body may contain separated fork/slit faces. Each
    // individual basis support must also connect through positive-weight nodes,
    // not borrow a connection elsewhere through zero-weight nodes.
    std::vector<std::size_t> offsets(nodes.size() + 1, 0);
    for (std::size_t i = 0; i < nodes.size(); ++i)
        offsets[i + 1] = offsets[i] + interpolation[i].size();
    std::vector<std::size_t> support_parent(offsets.back()), first(coarse.size(), missing);
    for (std::size_t i = 0; i < support_parent.size(); ++i)
        support_parent[i] = i;
    const auto support_root = [&support_parent](std::size_t id) {
        while (support_parent[id] != id) {
            support_parent[id] = support_parent[support_parent[id]];
            id                 = support_parent[id];
        }
        return id;
    };
    std::size_t work = 0;
    for (std::size_t e = 0; e < elements.size(); ++e) {
        if (checkpoint && e % 128 == 0)
            checkpoint();
        for (unsigned a = 0; a < 4; ++a)
            for (unsigned b = a + 1; b < 4; ++b) {
                const auto na = elements[e].nodes[a], nb = elements[e].nodes[b];
                for (std::size_t i = 0; i < interpolation[na].size(); ++i)
                    for (std::size_t j = 0; j < interpolation[nb].size(); ++j) {
                        if (work++ >= limits.max_assembly_work)
                            throw ReductionBudget("Reduction support work budget exceeded.");
                        if (interpolation[na][i].node == interpolation[nb][j].node)
                            support_parent[support_root(offsets[na] + i)] = support_root(offsets[nb] + j);
                    }
            }
    }
    for (std::size_t i = 0; i < nodes.size(); ++i)
        for (std::size_t j = 0; j < interpolation[i].size(); ++j) {
            const auto master = interpolation[i][j].node, representative = support_root(offsets[i] + j);
            if (first[master] != missing && first[master] != representative)
                throw std::invalid_argument("Reduction basis has disconnected positive support.");
            first[master] = representative;
        }
    Eigen::VectorXd prescribed = Eigen::VectorXd::Zero(3 * coarse.size());
    std::vector<bool> fixed(3 * coarse.size(), false), fine_fixed(3 * nodes.size(), false);
    std::size_t constraint_count = 0;
    for (const auto& constraint : constraints) {
        if (checkpoint && constraint_count++ % 128 == 0) checkpoint();
        if (constraint.node >= nodes.size() || constraint.axis > 2 || !std::isfinite(constraint.displacement_m))
            throw std::invalid_argument("Invalid reduced constraint.");
        const auto& row = interpolation[constraint.node];
        if (row.size() != 1 || row[0].weight != 1.)
            throw std::invalid_argument("Constrained node must have an exact unit prolongation row.");
        const auto index = 3 * row[0].node + constraint.axis;
        if (fixed[index] && prescribed[index] != constraint.displacement_m)
            throw std::invalid_argument("Conflicting reduced constraints.");
        fixed[index]                                      = true;
        fine_fixed[3 * constraint.node + constraint.axis] = true;
        prescribed[index]                                 = constraint.displacement_m;
    }
    std::vector<int> free_index(fixed.size(), -1);
    int free_count = 0;
    for (std::size_t i = 0; i < fixed.size(); ++i)
        if (!fixed[i])
            free_index[i] = free_count++;
    if (std::size_t(free_count) > limits.max_free_dofs)
        throw ReductionBudget("Reduced solve exceeds its free-DOF budget.");
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(free_count);
    std::vector<Eigen::Triplet<double>> entries;
    for (std::size_t n = 0; n < elements.size(); ++n) {
        if (checkpoint && n % 64 == 0)
            checkpoint();
        const auto& element = elements[n];
        const auto m        = matrices(nodes, element);
        std::vector<std::size_t> masters;
        for (auto node : element.nodes)
            for (const auto& entry : interpolation[node])
                masters.push_back(entry.node);
        std::sort(masters.begin(), masters.end());
        masters.erase(std::unique(masters.begin(), masters.end()), masters.end());
        const std::size_t dofs = 3 * masters.size();
        if (dofs * dofs > limits.max_assembly_work - work)
            throw ReductionBudget("Reduction assembly work budget exceeded.");
        work += dofs * dofs;
        Eigen::MatrixXd b = Eigen::MatrixXd::Zero(6, dofs);
        for (unsigned vertex = 0; vertex < 4; ++vertex)
            for (const auto& entry : interpolation[element.nodes[vertex]]) {
                const auto column = std::lower_bound(masters.begin(), masters.end(), entry.node) - masters.begin();
                b.block<6, 3>(0, 3 * column) += entry.weight * m.b.block<6, 3>(0, 3 * vertex);
            }
        const Eigen::MatrixXd k = m.volume * b.transpose() * m.d * b;
        const Eigen::VectorXd f = m.volume * b.transpose() * m.d * element.eigenstrain;
        if (!k.allFinite() || !f.allFinite())
            throw std::runtime_error("Reduced element assembly overflowed.");
        for (std::size_t i = 0; i < dofs; ++i) {
            const auto row = 3 * masters[i / 3] + i % 3;
            if (fixed[row])
                continue;
            rhs[free_index[row]] += f[i];
            for (std::size_t j = 0; j < dofs; ++j) {
                const auto col = 3 * masters[j / 3] + j % 3;
                rhs[free_index[row]] -= k(i, j) * prescribed[col];
                if (!fixed[col] && k(i, j) != 0.) {
                    if (entries.size() >= limits.max_assembly_entries)
                        throw ReductionBudget("Reduction matrix storage budget exceeded.");
                    entries.emplace_back(free_index[row], free_index[col], k(i, j));
                }
            }
        }
    }
    Result result;
    Eigen::VectorXd coarse_displacement = prescribed;
    if (free_count) {
        Eigen::SparseMatrix<double> k(free_count, free_count);
        k.setFromTriplets(entries.begin(), entries.end());
        if (checkpoint)
            checkpoint();
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> factor;
        factor.compute(k);
        if (checkpoint)
            checkpoint();
        if (factor.info() != Eigen::Success || !factor.vectorD().allFinite() ||
            factor.vectorD().minCoeff() <= 1e-12 * factor.vectorD().cwiseAbs().maxCoeff())
            throw std::runtime_error("Reduced system is singular or ill-conditioned.");
        const Eigen::VectorXd solution = factor.solve(rhs);
        if (factor.info() != Eigen::Success || !solution.allFinite())
            throw std::runtime_error("Reduced solve failed.");
        const double scale       = rhs.norm() + (k * solution).norm();
        result.relative_residual = scale > 0. ? (k * solution - rhs).norm() / scale : 0.;
        if (!std::isfinite(result.relative_residual) || result.relative_residual > 1e-8)
            throw std::runtime_error("Reduced equilibrium residual exceeds tolerance.");
        for (std::size_t i = 0; i < fixed.size(); ++i)
            if (!fixed[i])
                coarse_displacement[i] = solution[free_index[i]];
    }
    result.displacement_m.resize(nodes.size(), Eigen::Vector3d::Zero());
    result.reaction_n.resize(nodes.size(), Eigen::Vector3d::Zero());
    result.stress_pa.reserve(elements.size());
    std::vector<Eigen::Vector3d> free_residual(nodes.size(), Eigen::Vector3d::Zero());
    std::vector<Eigen::Vector3d> generalized_reaction(coarse.size(), Eigen::Vector3d::Zero());
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        if (checkpoint && i % 128 == 0)
            checkpoint();
        for (const auto& entry : interpolation[i])
            result.displacement_m[i] += entry.weight * coarse_displacement.segment<3>(3 * entry.node);
    }
    for (std::size_t n = 0; n < elements.size(); ++n) {
        if (checkpoint && n % 128 == 0)
            checkpoint();
        const auto& element = elements[n];
        const auto m        = matrices(nodes, element);
        Eigen::Matrix<double, 12, 1> u;
        for (unsigned i = 0; i < 4; ++i)
            u.segment<3>(3 * i) = result.displacement_m[element.nodes[i]];
        const Eigen::Matrix<double, 6, 1> strain = m.b * u - element.eigenstrain;
        const Eigen::Matrix<double, 6, 1> stress = m.d * strain;
        if (!stress.allFinite())
            throw std::runtime_error("Reduced stress overflowed.");
        result.stress_pa.push_back(stress);
        result.elastic_energy_j += .5 * m.volume * strain.dot(stress);
        const Eigen::Matrix<double, 12, 1> force = m.volume * m.b.transpose() * stress;
        for (unsigned i = 0; i < 12; ++i) {
            if (fine_fixed[3 * element.nodes[i / 3] + i % 3])
                result.reaction_n[element.nodes[i / 3]][i % 3] += force[i];
            else
                free_residual[element.nodes[i / 3]][i % 3] += force[i];
            for (const auto& entry : interpolation[element.nodes[i / 3]])
                if (fixed[3 * entry.node + i % 3])
                    generalized_reaction[entry.node][i % 3] += entry.weight * force[i];
        }
    }
    if (!std::isfinite(result.elastic_energy_j))
        throw std::runtime_error("Reduced energy overflowed.");
    if (checkpoint)
        checkpoint();
    double residual_squared = 0.;
    for (const auto& force : free_residual)
        residual_squared += force.squaredNorm();
    if (!std::isfinite(residual_squared))
        throw std::runtime_error("Fine residual overflowed.");
    return {std::move(result), std::move(generalized_reaction), std::sqrt(residual_squared), std::size_t(free_count)};
}

Result solve_affine_periodic(const std::vector<Eigen::Vector3d>& nodes_m,
                             const std::vector<Element>& elements,
                             const std::vector<std::size_t>& node_classes,
                             const std::vector<Eigen::Vector3d>& affine_displacement_m,
                             std::size_t anchor_node,
                             const std::function<void()>& checkpoint)
{
    if (nodes_m.size() > 2048 || elements.size() > 32768 || anchor_node >= nodes_m.size() || node_classes.size() != nodes_m.size() ||
        affine_displacement_m.size() != nodes_m.size())
        throw std::invalid_argument("Invalid or oversized periodic mechanical system.");
    // Conservative dense-fill/workspace allowance, not a promise that sparse
    // factorization always remains sparse. Check before allocating triplets.
    const std::size_t dofs = 3 * nodes_m.size();
    if (dofs*dofs*32 + elements.size()*144*32 > 512ULL*1024*1024)
        throw std::runtime_error("Periodic mechanical workspace budget exceeded.");
    std::vector<Constraint> constraints;
    for (unsigned axis = 0; axis < 3; ++axis)
        constraints.push_back({anchor_node, axis, affine_displacement_m[anchor_node][axis]});
    return solve_mapped(nodes_m, elements, constraints, node_classes, affine_displacement_m, checkpoint);
}
namespace {
std::vector<double> footprint_areas(const std::vector<Eigen::Vector3d> &nodes_m,
                          const std::vector<Element> &elements,
                          const std::vector<std::array<std::size_t, 3>> &contact_faces,
                          double bed_z_m, const std::function<void()> &checkpoint)
{
    if (!std::isfinite(bed_z_m) || contact_faces.empty())
        throw std::invalid_argument("A bonded-bed scenario requires a finite bed plane and contact faces.");
    using Face = std::array<std::size_t, 3>;
    std::map<Face, unsigned> incidence;
    for (std::size_t i = 0; i < elements.size(); ++i) {
        if (checkpoint && i % 128 == 0) checkpoint();
        for (unsigned omit = 0; omit < 4; ++omit) {
            Face face;
            unsigned j = 0;
            for (unsigned k = 0; k < 4; ++k)
                if (k != omit) face[j++] = elements[i].nodes[k];
            std::sort(face.begin(), face.end());
            if (++incidence[face] > 2)
                throw std::invalid_argument("Non-manifold thermoelastic mesh.");
        }
    }
    std::vector<double> areas(contact_faces.size(), 0.);
    std::set<Face> seen;
    for (std::size_t i = 0; i < contact_faces.size(); ++i) {
        if (checkpoint && i % 128 == 0) checkpoint();
        Face face = contact_faces[i];
        std::sort(face.begin(), face.end());
        const auto entry = incidence.find(face);
        if (!seen.insert(face).second || entry == incidence.end() || entry->second != 1)
            throw std::invalid_argument("Bed contact must contain distinct exterior tetrahedron faces.");
        for (std::size_t n : face)
            if (n >= nodes_m.size() || !nodes_m[n].allFinite() || std::abs(nodes_m[n].z() - bed_z_m) > 1e-9)
                throw std::invalid_argument("Bed contact face is not on the specified bed plane.");
        const double area = 0.5 * (nodes_m[face[1]] - nodes_m[face[0]]).cross(
            nodes_m[face[2]] - nodes_m[face[0]]).norm();
        if (!std::isfinite(area) || area <= 0.)
            throw std::invalid_argument("Degenerate bed contact face.");
        areas[i] = area;
    }
    return areas;
}

BedResult solve_bed_areas(const std::vector<Eigen::Vector3d> &nodes_m,
                         const std::vector<Element> &elements,
                         const std::vector<std::array<std::size_t, 3>> &faces,
                         const std::vector<double> &areas,
                         const std::function<void()> &checkpoint)
{
    BedResult result;
    result.contact_area_m2.assign(nodes_m.size(), 0.);
    result.tensile_traction_pa.assign(nodes_m.size(), 0.);
    for (std::size_t i = 0; i < faces.size(); ++i) {
        if (checkpoint && i % 128 == 0) checkpoint();
        for (const auto n : faces[i]) result.contact_area_m2[n] += areas[i] / 3.;
        result.contact_area_total_m2 += areas[i];
    }
    std::vector<Constraint> constraints;
    for (std::size_t n = 0; n < nodes_m.size(); ++n)
        if (result.contact_area_m2[n] > 0.)
            for (unsigned axis = 0; axis < 3; ++axis) constraints.push_back({n, axis});
    // Components with no bonded contact retain rigid-body modes and are
    // rejected by solve(), not artificially pinned to manufacture a result.
    result.equilibrium = solve(nodes_m, elements, constraints, checkpoint);
    for (std::size_t n = 0; n < nodes_m.size(); ++n) {
        const double area = result.contact_area_m2[n];
        if (area > 0.) {
            const double tension = std::max(0., -result.equilibrium.reaction_n[n].z());
            result.tensile_traction_pa[n] = tension / area;
            result.tensile_force_n += tension;
        }
    }
    return result;
}
} // namespace

BedResult solve_bonded_bed(const std::vector<Eigen::Vector3d> &nodes_m,
                          const std::vector<Element> &elements,
                          const std::vector<std::array<std::size_t, 3>> &contact_faces,
                          double bed_z_m, const std::function<void()> &checkpoint)
{
    const auto areas = footprint_areas(nodes_m, elements, contact_faces, bed_z_m, checkpoint);
    return solve_bed_areas(nodes_m, elements, contact_faces, areas, checkpoint);
}

PartialBedResult solve_partially_bonded_bed(
    const std::vector<Eigen::Vector3d> &nodes_m, const std::vector<Element> &elements,
    const std::vector<std::array<std::size_t, 3>> &footprint_faces,
    const std::vector<std::size_t> &bonded_face_indices, double bed_z_m,
    const std::function<void()> &checkpoint)
{
    auto bonded_areas = footprint_areas(nodes_m, elements, footprint_faces, bed_z_m, checkpoint);
    std::vector<bool> retained(footprint_faces.size(), false);
    for (const auto i : bonded_face_indices) {
        if (i >= retained.size() || retained[i]) throw std::invalid_argument("Invalid or duplicate bonded footprint face index.");
        retained[i] = true;
    }
    if (bonded_face_indices.empty())
        throw std::runtime_error("Released bed scenario has no restraints; rigid-body motion is unresolved.");
    PartialBedResult result;
    result.released_area_m2.assign(nodes_m.size(), 0.);
    result.footprint_gap_m.assign(nodes_m.size(), 0.);
    for (std::size_t i = 0; i < footprint_faces.size(); ++i) {
        if (checkpoint && i % 128 == 0) checkpoint();
        if (retained[i]) continue;
        for (const auto n : footprint_faces[i]) result.released_area_m2[n] += bonded_areas[i] / 3.;
        result.released_area_total_m2 += bonded_areas[i];
        bonded_areas[i] = 0.;
    }
    result.bonded = solve_bed_areas(nodes_m, elements, footprint_faces, bonded_areas, checkpoint);
    double coordinate_scale = std::abs(bed_z_m);
    for (std::size_t n = 0; n < nodes_m.size(); ++n) {
        if (result.released_area_m2[n] == 0. && result.bonded.contact_area_m2[n] == 0.) continue;
        coordinate_scale = std::max(coordinate_scale, nodes_m[n].norm());
        const double gap = result.bonded.equilibrium.displacement_m[n].z();
        result.footprint_gap_m[n] = gap;
        result.maximum_lift_m = std::max(result.maximum_lift_m, gap);
        result.maximum_penetration_m = std::max(result.maximum_penetration_m, -gap);
    }
    const double roundoff = 128. * std::numeric_limits<double>::epsilon() * coordinate_scale;
    result.nonpenetrating = result.maximum_penetration_m <= roundoff;
    return result;
}
ContactBedResult solve_partially_bonded_contact(
    const std::vector<Eigen::Vector3d> &nodes_m, const std::vector<Element> &elements,
    const std::vector<std::array<std::size_t, 3>> &footprint_faces,
    const std::vector<std::size_t> &bonded_face_indices, double bed_z_m,
    const BedContactLimits &limits, const std::function<void()> &checkpoint)
{
    if (nodes_m.size() > 12000 || elements.size() > 100000 || footprint_faces.size() > 400000 ||
        limits.max_iterations == 0 || limits.max_iterations > 128 ||
        !std::isfinite(limits.gap_tolerance_m) || limits.gap_tolerance_m <= 0. ||
        !std::isfinite(limits.reaction_tolerance_n) || limits.reaction_tolerance_n <= 0.)
        throw std::invalid_argument("Invalid bounded bed-contact solve limits.");
    ContactBedResult result;
    result.bed = solve_partially_bonded_bed(nodes_m, elements, footprint_faces,
        bonded_face_indices, bed_z_m, checkpoint);
    std::vector<bool> active(nodes_m.size(), false);
    for (unsigned iteration = 0; iteration < limits.max_iterations; ++iteration) {
        if (checkpoint) checkpoint();
        ++result.iterations;
        bool changed = false;
        const auto &equilibrium = result.bed.bonded.equilibrium;
        for (std::size_t n = 0; n < nodes_m.size(); ++n) {
            if (result.bed.released_area_m2[n] == 0. || result.bed.bonded.contact_area_m2[n] > 0.) continue;
            // The rigid plate can push upwards, but cannot pull a released
            // node down. Remove such restraints instead of silently rebonding.
            if (active[n] && equilibrium.reaction_n[n].z() < -limits.reaction_tolerance_n) {
                active[n] = false;
                changed = true;
            } else if (!active[n] && equilibrium.displacement_m[n].z() < -limits.gap_tolerance_m) {
                active[n] = true;
                changed = true;
            }
        }
        if (!changed) {
            auto &bed = result.bed;
            bed.maximum_lift_m = bed.maximum_penetration_m = 0.;
            bed.bonded.tensile_force_n = 0.;
            result.normal_reaction_n.assign(nodes_m.size(), 0.);
            for (std::size_t n = 0; n < nodes_m.size(); ++n) {
                const double area = bed.bonded.contact_area_m2[n];
                if (area == 0. && bed.released_area_m2[n] == 0.) continue;
                const double gap = equilibrium.displacement_m[n].z();
                bed.footprint_gap_m[n] = gap;
                bed.maximum_lift_m = std::max(bed.maximum_lift_m, gap);
                bed.maximum_penetration_m = std::max(bed.maximum_penetration_m, -gap);
                if (area > 0.) {
                    const double tension = std::max(0., -equilibrium.reaction_n[n].z());
                    bed.bonded.tensile_traction_pa[n] = tension / area;
                    bed.bonded.tensile_force_n += tension;
                } else if (active[n]) {
                    result.normal_reaction_n[n] = equilibrium.reaction_n[n].z();
                }
            }
            bed.nonpenetrating = bed.maximum_penetration_m <= limits.gap_tolerance_m;
            if (!bed.nonpenetrating) throw std::runtime_error("Bed contact remains penetrating.");
            return result;
        }
        if (iteration + 1 == limits.max_iterations) break;
        std::vector<Constraint> constraints;
        for (std::size_t n = 0; n < nodes_m.size(); ++n) {
            if (result.bed.bonded.contact_area_m2[n] > 0.)
                for (unsigned axis = 0; axis < 3; ++axis) constraints.push_back({n, axis});
            else if (active[n])
                constraints.push_back({n, 2});
        }
        result.bed.bonded.equilibrium = solve(nodes_m, elements, constraints, checkpoint);
    }
    throw std::runtime_error("Bed contact active set did not converge within its iteration budget.");
}
} // namespace Slic3r::Thermoelastic
