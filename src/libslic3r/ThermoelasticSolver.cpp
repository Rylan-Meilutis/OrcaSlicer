#include "ThermoelasticSolver.hpp"

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
    if (!std::isfinite(e.young_pa) || e.young_pa <= 0. ||
        !std::isfinite(e.poisson) || e.poisson <= -1. || e.poisson >= 0.5 || !e.eigenstrain.allFinite())
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
    const double mu = e.young_pa / (2. * (1. + e.poisson));
    const double lambda = e.young_pa * e.poisson / ((1. + e.poisson) * (1. - 2. * e.poisson));
    result.d.topLeftCorner<3, 3>().setConstant(lambda);
    result.d.diagonal().head<3>().array() += 2. * mu;
    result.d.diagonal().tail<3>().setConstant(mu);
    if (!result.d.allFinite() || !result.b.allFinite())
        throw std::invalid_argument("Non-finite thermoelastic element coefficients.");
    return result;
}
} // namespace

Result solve(const std::vector<Eigen::Vector3d> &nodes_m, const std::vector<Element> &elements,
             const std::vector<Constraint> &constraints, const std::function<void()> &checkpoint)
{
    if (nodes_m.empty() || elements.empty() || nodes_m.size() > std::size_t(std::numeric_limits<int>::max() / 3))
        throw std::invalid_argument("Empty or oversized thermoelastic mesh.");
    for (const auto &node : nodes_m)
        if (!node.allFinite()) throw std::invalid_argument("Non-finite thermoelastic node.");
    const int dofs = int(nodes_m.size() * 3);
    Eigen::VectorXd displacement = Eigen::VectorXd::Zero(dofs);
    std::vector<bool> fixed(dofs, false);
    for (const auto &c : constraints) {
        if (c.node >= nodes_m.size() || c.axis > 2 || !std::isfinite(c.displacement_m))
            throw std::invalid_argument("Invalid thermoelastic constraint.");
        const int i = int(3 * c.node + c.axis);
        if (fixed[i] && displacement[i] != c.displacement_m)
            throw std::invalid_argument("Conflicting thermoelastic constraints.");
        fixed[i] = true;
        displacement[i] = c.displacement_m;
    }
    std::vector<int> free_index(dofs, -1);
    int free_count = 0;
    for (int i = 0; i < dofs; ++i)
        if (!fixed[i]) free_index[i] = free_count++;
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
            if (fixed[row]) continue;
            rhs[free_index[row]] += f[i];
            for (int j = 0; j < 12; ++j) {
                const int col = int(3 * e.nodes[j / 3] + j % 3);
                if (fixed[col]) rhs[free_index[row]] -= k(i, j) * displacement[col];
                else if (k(i, j) != 0.) entries.emplace_back(free_index[row], free_index[col], k(i, j));
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
            if (!fixed[i]) displacement[i] = u[free_index[i]];
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
            if (fixed[3 * e.nodes[i/3] + i%3]) result.reaction_n[e.nodes[i/3]][i%3] += force[i];
    }
    if (!std::isfinite(result.elastic_energy_j))
        throw std::runtime_error("Non-finite thermoelastic energy.");
    if (checkpoint) checkpoint();
    return result;
}
BedResult solve_bonded_bed(const std::vector<Eigen::Vector3d> &nodes_m,
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
    BedResult result;
    result.contact_area_m2.assign(nodes_m.size(), 0.);
    result.tensile_traction_pa.assign(nodes_m.size(), 0.);
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
        for (std::size_t n : face) result.contact_area_m2[n] += area / 3.;
        result.contact_area_total_m2 += area;
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
} // namespace Slic3r::Thermoelastic
