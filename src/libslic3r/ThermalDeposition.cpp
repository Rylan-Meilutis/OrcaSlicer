#include "ThermalDeposition.hpp"
#include "ThermalHistory.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

namespace Slic3r::ThermalHistory {
HistoryResult evolve(const std::vector<DepositionEvent> &events, double end_time_s,
    const HistoryLimits &limits, const std::function<void()> &checkpoint, const HistoryObservations& observations)
{
    if (events.size() > limits.max_events) {
        HistoryResult result;
        result.state = State::BudgetExceeded;
        return result;
    }
    std::size_t cursor = 0;
    return evolve_generated([&](DepositionEvent &event) {
        if (cursor == events.size()) return false;
        event = events[cursor++];
        return true;
    }, end_time_s, limits, checkpoint, observations);
}

HistoryResult evolve_generated(const std::function<bool(DepositionEvent &)> &next, double end_time_s,
    const HistoryLimits &limits, const std::function<void()> &checkpoint, const HistoryObservations& observations)
{
    HistoryResult result;
    if (!next || !std::isfinite(end_time_s) || end_time_s < 0. ||
        !std::isfinite(limits.max_step_s) || limits.max_step_s <= 0. ||
        !std::isfinite(limits.max_step_temperature_error_k) || limits.max_step_temperature_error_k < 0.) return result;
    if (observations.times_s.size() > 4096) { result.state = State::BudgetExceeded; return result; }
    if (!observations.times_s.empty() && !observations.observe) return result;
    double previous_observation = -1.;
    for (const double at : observations.times_s) {
        if (!std::isfinite(at) || at < 0. || at > end_time_s || at <= previous_observation) return result;
        previous_observation = at;
    }
    std::vector<Cell> cells;
    std::vector<Connection> connections;
    std::vector<Boundary> boundaries;
    std::vector<Element> elements;
    std::size_t internal_elements = 0, element_work = 0;
    std::size_t coalescence_work = 0;
    double time = 0.;
    std::size_t work = 0;
    double proposed_step = limits.max_step_s;
    const auto solve_step = [&](const std::vector<Cell> &input, double dt) {
        Result solved;
        if (result.solver_steps >= limits.max_steps || input.size() > limits.max_cell_steps - work ||
            elements.size() > limits.max_element_steps-element_work) {
            solved.state = State::BudgetExceeded;
            return solved;
        }
        ++result.solver_steps;
        work += input.size();
        element_work += elements.size();
        solved=step(input, connections, boundaries, dt, limits.solve, checkpoint, elements);
        if (solved.rejected_temperature) result.rejected_temperature=solved.rejected_temperature;
        return solved;
    };
    const auto advance = [&](double until) {
        if (!std::isfinite(until) || until < time || until > end_time_s) return State::InvalidInput;
        if (cells.empty()) { time = until; return State::Converged; }
        while (time < until) {
            if (checkpoint) checkpoint();
            const double dt = std::min(proposed_step, until - time);
            if (time + dt <= time) return State::InvalidInput;
            auto solved = solve_step(cells, dt);
            if (solved.state != State::Converged) return solved.state;
            if (limits.max_step_temperature_error_k > 0.) {
                // The two solves use relative durations with fixed coefficients
                // between events. Their midpoint need not be representable on
                // the absolute clock (e.g. adjacent floating-point event times).
                if (dt * .5 <= 0.) return State::Unconverged;
                const auto first = solve_step(cells, dt * .5);
                if (first.state != State::Converged) return first.state;
                auto midpoint = cells;
                for (std::size_t i = 0; i < cells.size(); ++i)
                    midpoint[i].temperature_k = first.temperatures_k[i];
                auto second = solve_step(midpoint, dt * .5);
                if (second.state != State::Converged) return second.state;
                double error = 0.;
                for (std::size_t i = 0; i < cells.size(); ++i)
                    error = std::max(error, std::abs(solved.temperatures_k[i] - second.temperatures_k[i]));
                if (error > limits.max_step_temperature_error_k) {
                    ++result.rejected_steps;
                    proposed_step = dt * .5;
                    continue; // Neither temperature nor heat accounting is committed.
                }
                result.maximum_step_temperature_error_k = std::max(result.maximum_step_temperature_error_k, error);
                second.source_energy_j += first.source_energy_j;
                second.boundary_energy_j += first.boundary_energy_j;
                for (std::size_t kind = 0; kind < boundary_kind_count; ++kind)
                    second.boundary_energy_by_kind_j[kind] += first.boundary_energy_by_kind_j[kind];
                solved = std::move(second);
                proposed_step = error < limits.max_step_temperature_error_k * .25 ?
                    std::min(limits.max_step_s, dt * 2.) : dt;
            }
            for (std::size_t i = 0; i < cells.size(); ++i) cells[i].temperature_k = solved.temperatures_k[i];
            result.source_energy_j += solved.source_energy_j;
            result.boundary_energy_j += solved.boundary_energy_j;
            for (std::size_t kind = 0; kind < boundary_kind_count; ++kind)
                result.boundary_energy_by_kind_j[kind] += solved.boundary_energy_by_kind_j[kind];
            ++result.steps;
            result.maximum_accepted_interval_s = std::max(result.maximum_accepted_interval_s, dt);
            time += dt;
        }
        return State::Converged;
    };
    const auto observe_before = [&](double until, bool inclusive) {
        while (result.completed_observations < observations.times_s.size()) {
            const double at = observations.times_s[result.completed_observations];
            if (at > until || (!inclusive && at == until)) break;
            const auto state = advance(at);
            if (state != State::Converged) return state;
            if (checkpoint) checkpoint();
            observations.observe(at,cells);
            ++result.completed_observations;
        }
        return State::Converged;
    };
    std::size_t event_count = 0;
    while (true) {
        if (checkpoint) checkpoint();
        DepositionEvent event;
        if (!next(event)) break;
        if (event_count++ >= limits.max_events) { result.state = State::BudgetExceeded; return result; }
        // Validate the incoming event clock before advancing any observations.
        if (!std::isfinite(event.time_s) || event.time_s < time || event.time_s > end_time_s) {
            result.state = State::InvalidInput; return result;
        }
        result.state = observe_before(event.time_s,false);
        if (result.state != State::Converged) return result;
        result.state = advance(event.time_s);
        if (result.state != State::Converged) return result;
        if (event.coalesced_nodes) {
            const auto& mapping = *event.coalesced_nodes;
            if (mapping.size() != cells.size() || mapping.empty()) {
                result.state = State::InvalidInput; return result;
            }
            const auto charge = [&](std::size_t count) {
                if (count > limits.max_coalescence_work-coalescence_work) return false;
                coalescence_work += count;
                return true;
            };
            if (!charge(cells.size()) || !charge(connections.size()) || !charge(internal_elements*8)) {
                result.state = State::BudgetExceeded; return result;
            }
            std::size_t count = 0;
            for (const auto id : mapping) {
                if (id >= cells.size()) { result.state = State::InvalidInput; return result; }
                count = std::max(count,id+1);
            }
            struct Accumulated { double volume = 0., capacity = 0., heat = 0., source = 0.; };
            std::vector<Accumulated> sums(count);
            for (std::size_t i = 0; i < cells.size(); ++i) {
                if (checkpoint) checkpoint();
                const auto& cell = cells[i];
                auto& sum = sums[mapping[i]];
                const double capacity = cell.volume_m3*cell.volumetric_heat_capacity_j_m3_k;
                sum.volume += cell.volume_m3;
                sum.capacity += capacity;
                sum.heat += capacity*cell.temperature_k;
                sum.source += cell.source_w;
            }
            std::vector<Cell> merged;
            merged.reserve(count);
            for (const auto& sum : sums) {
                const double density = sum.capacity/sum.volume, temperature = sum.heat/sum.capacity;
                if (!std::isfinite(sum.volume) || sum.volume <= 0. || !std::isfinite(density) || density <= 0. ||
                    !std::isfinite(temperature) || temperature <= 0. || !std::isfinite(sum.source)) {
                    result.state = State::InvalidInput; return result;
                }
                merged.push_back({sum.volume,density,temperature,sum.source});
            }
            for (std::size_t i = 0; i < internal_elements; ++i) {
                if (checkpoint) checkpoint();
                auto& element = elements[i];
                for (auto& node : element.nodes) node = mapping[node];
                // Coalescing distinct corners within an element is a geometry
                // collapse, not an admissible contact join.
                if (!valid_element(element,count)) { result.state = State::InvalidInput; return result; }
            }
            for (auto& link : connections) {
                link.first = mapping[link.first];
                link.second = mapping[link.second];
            }
            connections.erase(std::remove_if(connections.begin(),connections.end(),[](const Connection& link) {
                return link.first == link.second;
            }),connections.end());
            cells = std::move(merged);
        }
        if (event.deposited_cells.size() > limits.solve.max_cells - cells.size() ||
            event.added_connections.size() > limits.solve.max_connections - connections.size() ||
            event.exposed_boundaries.size() > limits.solve.max_boundaries ||
            event.added_mass.size() > limits.solve.max_cells ||
            event.added_elements.size() > limits.solve.max_elements-internal_elements ||
            event.exposed_elements.size() > limits.solve.max_elements-internal_elements-event.added_elements.size()) {
            result.state = State::BudgetExceeded; return result;
        }
        for (const auto &cell : event.deposited_cells) {
            const double capacity = cell.volume_m3 * cell.volumetric_heat_capacity_j_m3_k;
            if (!std::isfinite(capacity) || capacity <= 0. || cell.volume_m3 <= 0. ||
                cell.volumetric_heat_capacity_j_m3_k <= 0. || !std::isfinite(cell.temperature_k) ||
                cell.temperature_k <= 0. || !std::isfinite(cell.source_w)) {
                result.state = State::InvalidInput; return result;
            }
            result.deposited_energy_j += capacity * cell.temperature_k;
            cells.push_back(cell);
        }
        for (const auto& addition : event.added_mass) {
            const double capacity = addition.volume_m3 * addition.volumetric_heat_capacity_j_m3_k;
            if (addition.cell >= cells.size() || !std::isfinite(capacity) || capacity <= 0. ||
                addition.volume_m3 <= 0. || addition.volumetric_heat_capacity_j_m3_k <= 0. ||
                !std::isfinite(addition.temperature_k) || addition.temperature_k <= 0.) {
                result.state = State::InvalidInput; return result;
            }
            auto& cell = cells[addition.cell];
            const double old_capacity = cell.volume_m3 * cell.volumetric_heat_capacity_j_m3_k;
            const double total_capacity = old_capacity+capacity;
            const double volume = cell.volume_m3+addition.volume_m3;
            const double temperature = cell.temperature_k*(old_capacity/total_capacity)+
                addition.temperature_k*(capacity/total_capacity);
            const double volumetric_capacity = total_capacity/volume;
            if (!std::isfinite(total_capacity) || !std::isfinite(volume) || !std::isfinite(temperature) ||
                !std::isfinite(volumetric_capacity) || volumetric_capacity <= 0. || temperature <= 0.) {
                result.state = State::InvalidInput; return result;
            }
            cell = {volume,volumetric_capacity,temperature,cell.source_w};
            result.deposited_energy_j += capacity*addition.temperature_k;
        }
        elements.resize(internal_elements);
        for (const auto& element : event.added_elements) {
            if (element.boundary || !valid_element(element,cells.size())) {
                result.state = State::InvalidInput; return result;
            }
            elements.push_back(element);
        }
        internal_elements = elements.size();
        for (const auto& element : event.exposed_elements) {
            if (!element.boundary || !valid_element(element,cells.size())) {
                result.state = State::InvalidInput; return result;
            }
            elements.push_back(element);
        }
        for (const auto &link : event.added_connections) {
            if (link.first >= cells.size() || link.second >= cells.size() || link.first == link.second ||
                !std::isfinite(link.conductance_w_k) || link.conductance_w_k < 0.) {
                result.state = State::InvalidInput; return result;
            }
            connections.push_back(link);
        }
        for (const auto &boundary : event.exposed_boundaries)
            if (boundary.cell >= cells.size() || !std::isfinite(boundary.conductance_w_k) ||
                boundary.conductance_w_k < 0. || !std::isfinite(boundary.ambient_temperature_k) ||
                boundary.ambient_temperature_k <= 0. ||
                static_cast<std::size_t>(boundary.kind) >= boundary_kind_count) {
                result.state = State::InvalidInput; return result;
            }
        boundaries = event.exposed_boundaries;
    }
    result.state = observe_before(end_time_s,true);
    if (result.state != State::Converged) return result;
    result.state = advance(end_time_s);
    if (result.state != State::Converged) return result;
    for (const auto &cell : cells)
        result.final_energy_j += cell.volume_m3 * cell.volumetric_heat_capacity_j_m3_k * cell.temperature_k;
    result.energy_balance_error_j = result.final_energy_j - result.deposited_energy_j
        - result.source_energy_j - result.boundary_energy_j;
    if (!std::isfinite(result.energy_balance_error_j)) { result.state = State::InvalidInput; return result; }
    for (const auto &cell : cells) result.temperatures_k.push_back(cell.temperature_k);
    result.temporal_error_controlled = limits.max_step_temperature_error_k > 0.;
    return result;
}
} // namespace Slic3r::ThermalHistory
