#pragma once
#include "ThermalHistory.hpp"
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

namespace Slic3r::ThermalHistory {
struct AddedMass {
    std::size_t cell = 0;
    double volume_m3 = 0., volumetric_heat_capacity_j_m3_k = 0., temperature_k = 0.;
};
// Cell IDs are append order unless an explicit coalescence reindexes them.
// An event deposits NEW mass, never reheats
// an existing cell by silently resetting its temperature. Boundary replacement
// explicitly removes covered surfaces and supplies current effective h*A.
struct DepositionEvent {
    double time_s = 0.;
    std::vector<Cell> deposited_cells;
    std::vector<Connection> added_connections;
    std::vector<Boundary> exposed_boundaries;
    // Q1 mass-lumped activation on shared nodes. Add heat/capacity, never reset
    // existing temperatures. Caller owns nonoverlap and deposition provenance.
    std::vector<AddedMass> added_mass;
    std::vector<Element> added_elements;
    std::vector<Element> exposed_elements; // Replaces prior boundary operators.
    // Optional old-node -> new-node dense map, applied BEFORE this event's
    // additions. Only merges are permitted: every old node contributes all
    // its mass, heat and source to one new node. New IDs must cover [0,N).
    // Caller certifies that joined nodes represent the same physical location
    // after contact, never merely membership in the final printed component.
    // All event indices refer to the new numbering; old operators are remapped.
    std::optional<std::vector<std::size_t>> coalesced_nodes;
};
struct HistoryLimits {
    Limits solve;
    // max_steps and max_cell_steps include rejected error-control trials.
    std::size_t max_events = 10000, max_steps = 100000;
    std::size_t max_cell_steps = 10000000;
    double max_step_s = 1.;
    // Zero retains fixed-step backward Euler. Positive enables step doubling:
    // the maximum absolute cell-temperature difference between one full step
    // and two half steps must not exceed this local error tolerance.
    double max_step_temperature_error_k = 0.;
    // Includes rejected trials, independently of the active-node work budget.
    std::size_t max_element_steps = 10000000;
    std::size_t max_coalescence_work = 10000000;
};
struct HistoryResult {
    State state = State::InvalidInput;
    std::vector<double> temperatures_k;
    std::size_t steps = 0; // Accepted time intervals; adaptive intervals use two half steps.
    std::size_t solver_steps = 0, rejected_steps = 0;
    // Accepted interval before the adaptive pair of half steps, if enabled.
    // Zero means no active-cell time interval was integrated.
    double maximum_accepted_interval_s = 0.;
    double maximum_step_temperature_error_k = 0.;
    bool temporal_error_controlled = false;
    std::size_t completed_observations = 0;
    double deposited_energy_j = 0.; // Sum C*T at activation, Kelvin reference.
    double source_energy_j = 0., boundary_energy_j = 0.;
    BoundaryEnergy boundary_energy_by_kind_j{}; // Accepted steps only.
    double final_energy_j = 0., energy_balance_error_j = 0.;
    std::optional<RejectedTemperature> rejected_temperature;
};
struct HistoryObservations {
    // Strictly increasing, at most 4096 times in [0,end_time]. Observations at
    // an event time include ALL deposition/boundary events at that time.
    std::vector<double> times_s;
    // Borrowed active cells in the latest event's numbering; consume synchronously.
    // Coalescence changes IDs, so the producer must retain its geometry mapping.
    // Snapshots are
    // provisional until evolve returns Converged, and are not convergence
    // certificates. The caller must discard derived output if history fails.
    std::function<void(double, const std::vector<Cell>&)> observe;
};
// Ordered explicit physical times, not G-code line numbers or inferred PWM
// cooling. No inactive cells exist in the system. No partial field on failure.
HistoryResult evolve(const std::vector<DepositionEvent> &events, double end_time_s,
    const HistoryLimits &limits = {}, const std::function<void()> &checkpoint = {},
    const HistoryObservations& observations = {});

// Pull one event at a time instead of retaining every exposed-boundary snapshot.
// The producer returns false at end and owns its own construction-work budget.
// Events are consumed in order; cancellation/producer exceptions propagate.
HistoryResult evolve_generated(const std::function<bool(DepositionEvent &)> &next,
    double end_time_s, const HistoryLimits &limits = {},
    const std::function<void()> &checkpoint = {}, const HistoryObservations& observations = {});
} // namespace Slic3r::ThermalHistory
