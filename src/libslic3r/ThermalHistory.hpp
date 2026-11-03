#pragma once

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

namespace Slic3r::ThermalHistory {

// Only deposited, active cells belong in the system. The caller adds deposited
// mass and its initial energy, and updates boundary exposure between steps.
struct Cell {
    double volume_m3 = 0.;
    double volumetric_heat_capacity_j_m3_k = 0.;
    double temperature_k = 0.;
    double source_w = 0.;
};

struct Connection {
    std::size_t first = 0, second = 0;
    double conductance_w_k = 0.;
};

enum class BoundaryKind { Ambient, Bed, PartFan, AuxiliaryFan, DockFan, Exhaust, Mixed };
constexpr std::size_t boundary_kind_count = 7;
using BoundaryEnergy = std::array<double, boundary_kind_count>;
struct Boundary {
    std::size_t cell = 0;
    double conductance_w_k = 0.; // Explicit effective h*A, not fan PWM.
    double ambient_temperature_k = 0.;
    BoundaryKind kind = BoundaryKind::Ambient;
};

// Integrated Q1 operator on active nodes. Internal conduction has zero load
// and row sums; boundary operators include a reservoir load and heat channel.
// This numerical input does not certify geometric connectivity or exposure.
struct Element {
    std::array<std::size_t, 8> nodes{};
    Eigen::Matrix<double, 8, 8> conductance = Eigen::Matrix<double, 8, 8>::Zero();
    Eigen::Matrix<double, 8, 1> load = Eigen::Matrix<double, 8, 1>::Zero();
    bool boundary = false;
    BoundaryKind kind = BoundaryKind::Ambient;
};
bool valid_element(const Element& element, std::size_t active_nodes);

struct Limits {
    std::size_t max_cells = 200000;
    std::size_t max_connections = 1000000;
    std::size_t max_boundaries = 1000000;
    unsigned max_iterations = 2000;
    double relative_tolerance = 1e-10;
    std::size_t max_elements = 100000;
};

enum class State { Converged, InvalidInput, BudgetExceeded, Unconverged };
struct RejectedTemperature {
    std::size_t node;
    double kelvin;
};
struct Result {
    State state = State::InvalidInput;
    std::vector<double> temperatures_k;
    unsigned iterations = 0;
    // Maximum of the global relative L2 residual and each equation's residual
    // scaled by |rhs_i| + diagonal_i*|temperature_i|. Both must converge.
    double relative_residual = 0.;
    double stored_energy_change_j = 0.;
    double source_energy_j = 0.;
    double boundary_energy_j = 0.; // Positive into the cells, evaluated at end of step.
    // Signed heat exchange, not a causal attribution of stress to a fan. The
    // explicit boundary conductances must already account for combined airflow.
    BoundaryEnergy boundary_energy_by_kind_j{};
    double energy_balance_error_j = 0.;
    // Failure diagnostic only; an inadmissible solve still publishes no field.
    std::optional<RejectedTemperature> rejected_temperature;
};

// One backward-Euler step with fixed positive heat capacities/conductances.
// Temperature-dependent properties require caller-controlled re-evaluation and
// convergence. Cancellation exceptions propagate; failure publishes no field.
Result step(const std::vector<Cell> &cells, const std::vector<Connection> &connections,
            const std::vector<Boundary> &boundaries, double seconds,
            const Limits &limits = {}, const std::function<void()> &checkpoint = {},
            const std::vector<Element>& elements = {});

} // namespace Slic3r::ThermalHistory
