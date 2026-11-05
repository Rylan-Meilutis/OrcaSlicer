#include "ThermalHistory.hpp"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

namespace Slic3r::ThermalHistory {
namespace {
bool positive(double value) { return std::isfinite(value) && value > 0.; }
bool nonnegative(double value) { return std::isfinite(value) && value >= 0.; }
double dot(const std::vector<double> &a, const std::vector<double> &b)
{
    double value = 0.;
    for (std::size_t i = 0; i < a.size(); ++i) value += a[i] * b[i];
    return value;
}
} // namespace

bool valid_element(const Element& element, std::size_t count)
{
    if (!element.conductance.allFinite() || !element.load.allFinite() ||
        !element.conductance.isApprox(element.conductance.transpose(),1e-12) ||
        static_cast<std::size_t>(element.kind) >= boundary_kind_count) return false;
    const double scale = element.conductance.cwiseAbs().maxCoeff();
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,8,8>> eigen(element.conductance);
    if (eigen.info() != Eigen::Success || eigen.eigenvalues().minCoeff() < -1e-12*scale) return false;
    if (!element.boundary && (element.load.cwiseAbs().maxCoeff() != 0. ||
        element.conductance.rowwise().sum().cwiseAbs().maxCoeff() > 1e-12*scale)) return false;
    if (element.boundary && (element.conductance.minCoeff() < -1e-12*scale ||
        element.load.minCoeff() < -1e-12*element.load.cwiseAbs().maxCoeff())) return false;
    for (unsigned a = 0; a < 8; ++a) {
        if (element.nodes[a] >= count) return false;
        for (unsigned b = 0; b < a; ++b) if (element.nodes[a] == element.nodes[b]) return false;
    }
    return true;
}

Result step(const std::vector<Cell> &cells, const std::vector<Connection> &connections,
            const std::vector<Boundary> &boundaries, double seconds,
            const Limits &limits, const std::function<void()> &checkpoint,
            const std::vector<Element>& elements)
{
    Result result;
    if (checkpoint) checkpoint();
    if (cells.size() > limits.max_cells || connections.size() > limits.max_connections ||
        boundaries.size() > limits.max_boundaries || elements.size() > limits.max_elements) {
        result.state = State::BudgetExceeded;
        return result;
    }
    if (cells.empty() || !positive(seconds) || !positive(limits.relative_tolerance) ||
        limits.relative_tolerance >= 1.) return result;
    const std::size_t count = cells.size();
    std::vector<double> capacity(count), diagonal(count), rhs(count), x(count);
    for (std::size_t i = 0; i < count; ++i) {
        if (i % 4096 == 0 && checkpoint) checkpoint();
        const auto &cell = cells[i];
        if (!positive(cell.volume_m3) || !positive(cell.volumetric_heat_capacity_j_m3_k) ||
            !positive(cell.temperature_k) || !std::isfinite(cell.source_w)) return result;
        capacity[i] = cell.volume_m3 * cell.volumetric_heat_capacity_j_m3_k;
        diagonal[i] = capacity[i] / seconds;
        rhs[i] = diagonal[i] * cell.temperature_k + cell.source_w;
        x[i] = cell.temperature_k;
        if (!positive(diagonal[i]) || !positive(capacity[i]) || !std::isfinite(rhs[i])) return result;
    }
    for (std::size_t i = 0; i < connections.size(); ++i) {
        if (i % 4096 == 0 && checkpoint) checkpoint();
        const auto &edge = connections[i];
        if (edge.first >= count || edge.second >= count || edge.first == edge.second ||
            !nonnegative(edge.conductance_w_k)) return result;
        diagonal[edge.first] += edge.conductance_w_k;
        diagonal[edge.second] += edge.conductance_w_k;
    }
    for (std::size_t i = 0; i < boundaries.size(); ++i) {
        if (i % 4096 == 0 && checkpoint) checkpoint();
        const auto &boundary = boundaries[i];
        if (boundary.cell >= count || !nonnegative(boundary.conductance_w_k) ||
            !positive(boundary.ambient_temperature_k) ||
            static_cast<std::size_t>(boundary.kind) >= boundary_kind_count) return result;
        diagonal[boundary.cell] += boundary.conductance_w_k;
        rhs[boundary.cell] += boundary.conductance_w_k * boundary.ambient_temperature_k;
    }
    for (std::size_t i = 0; i < count; ++i)
        if (!positive(diagonal[i]) || !std::isfinite(rhs[i])) return result;

    for (const auto& element : elements) {
        if (checkpoint) checkpoint();
        if (!valid_element(element,count)) return result;
        for (unsigned a = 0; a < 8; ++a) {
            diagonal[element.nodes[a]] += element.conductance(a,a);
            rhs[element.nodes[a]] += element.load[a];
        }
    }
    for (std::size_t i = 0; i < count; ++i)
        if (!positive(diagonal[i]) || !std::isfinite(rhs[i])) return result;

    const auto multiply = [&](const std::vector<double> &in, std::vector<double> &out) {
        for (std::size_t i = 0; i < count; ++i) out[i] = diagonal[i] * in[i];
        for (const auto &edge : connections) {
            out[edge.first] -= edge.conductance_w_k * in[edge.second];
            out[edge.second] -= edge.conductance_w_k * in[edge.first];
        }
        for (const auto& element : elements)
            for (unsigned a = 0; a < 8; ++a)
                for (unsigned b = 0; b < 8; ++b)
                    if (a != b) out[element.nodes[a]] += element.conductance(a,b) * in[element.nodes[b]];
    };
    std::vector<double> residual(count), direction(count), product(count), preconditioned(count);
    multiply(x, product);
    for (std::size_t i = 0; i < count; ++i) {
        residual[i] = rhs[i] - product[i];
        direction[i] = preconditioned[i] = residual[i] / diagonal[i];
    }
    const double rhs_norm = std::sqrt(dot(rhs, rhs));
    if (!positive(rhs_norm)) return result;
    const auto convergence_error = [&] {
        double error = std::sqrt(dot(residual,residual))/rhs_norm;
        // A global norm alone can ignore tiny occupied cut cells beside a
        // large heat capacity, even when their temperatures are far from the
        // solution. Also require a locally scaled equation residual at every
        // node. No temperature is clamped or small cell discarded.
        for (std::size_t i=0;i<count;++i) {
            const double scale=std::abs(rhs[i])+diagonal[i]*std::abs(x[i]);
            const double local=scale>0. ? std::abs(residual[i])/scale :
                (residual[i]==0. ? 0. : std::numeric_limits<double>::infinity());
            error=std::max(error,local);
        }
        return error;
    };
    double rz = dot(residual, preconditioned);
    result.relative_residual = convergence_error();
    while (result.relative_residual > limits.relative_tolerance && result.iterations < limits.max_iterations) {
        if (checkpoint) checkpoint();
        multiply(direction, product);
        const double denominator = dot(direction, product);
        if (!positive(denominator) || !positive(rz)) break;
        const double alpha = rz / denominator;
        for (std::size_t i = 0; i < count; ++i) {
            x[i] += alpha * direction[i];
            residual[i] -= alpha * product[i];
            preconditioned[i] = residual[i] / diagonal[i];
        }
        const double next_rz = dot(residual, preconditioned);
        const double beta = next_rz / rz;
        for (std::size_t i = 0; i < count; ++i) direction[i] = preconditioned[i] + beta * direction[i];
        rz = next_rz;
        ++result.iterations;
        result.relative_residual = convergence_error();
    }
    // Check the true residual, not just the recursively updated CG residual.
    multiply(x, product);
    for (std::size_t i = 0; i < count; ++i) residual[i] = rhs[i] - product[i];
    result.relative_residual = convergence_error();
    if (!std::isfinite(result.relative_residual) || result.relative_residual > limits.relative_tolerance) {
        result.state = State::Unconverged;
        return result;
    }
    for (std::size_t i = 0; i < count; ++i) {
        if (!positive(x[i])) {
            result.rejected_temperature=RejectedTemperature{i,x[i]};
            return result;
        }
        result.stored_energy_change_j += capacity[i] * (x[i] - cells[i].temperature_k);
        result.source_energy_j += seconds * cells[i].source_w;
    }
    for (const auto &boundary : boundaries) {
        const double energy = seconds * boundary.conductance_w_k *
                              (boundary.ambient_temperature_k - x[boundary.cell]);
        result.boundary_energy_j += energy;
        result.boundary_energy_by_kind_j[static_cast<std::size_t>(boundary.kind)] += energy;
    }
    for (const auto& element : elements) if (element.boundary) {
        Eigen::Matrix<double,8,1> local;
        for (unsigned a = 0; a < 8; ++a) local[a] = x[element.nodes[a]];
        const double energy = seconds * (element.load-element.conductance*local).sum();
        result.boundary_energy_j += energy;
        result.boundary_energy_by_kind_j[static_cast<std::size_t>(element.kind)] += energy;
    }
    result.energy_balance_error_j = result.stored_energy_change_j - result.source_energy_j - result.boundary_energy_j;
    if (!std::isfinite(result.energy_balance_error_j)) return result;
    result.temperatures_k = std::move(x);
    result.state = State::Converged;
    return result;
}
} // namespace Slic3r::ThermalHistory
