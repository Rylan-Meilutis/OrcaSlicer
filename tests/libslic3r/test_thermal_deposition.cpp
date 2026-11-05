#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ThermalDeposition.hpp"
#include "libslic3r/ThermalHistory.hpp"
#include <cstddef>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace Slic3r::ThermalHistory;
using Catch::Matchers::WithinAbs;

TEST_CASE("Inadmissible thermal temperatures report the failed node without publishing a field", "[ThermalDeposition]")
{
    const std::vector<Cell> cells{{1.,1.,400.,-1000.}};
    const auto single=step(cells,{}, {},1.);
    REQUIRE(single.state==State::InvalidInput);
    CHECK(single.temperatures_k.empty());
    REQUIRE(single.rejected_temperature.has_value());
    CHECK(single.rejected_temperature->node==0);
    CHECK_THAT(single.rejected_temperature->kelvin,WithinAbs(-600.,1e-8));
    DepositionEvent event;
    event.deposited_cells=cells;
    const auto history=evolve({event},1.);
    REQUIRE(history.state==State::InvalidInput);
    CHECK(history.temperatures_k.empty());
    REQUIRE(history.rejected_temperature.has_value());
    CHECK(history.rejected_temperature->node==0);
    CHECK_THAT(history.rejected_temperature->kelvin,WithinAbs(-600.,1e-8));
}

TEST_CASE("Thermal contact joins preserve separate histories until physical contact", "[ThermalDeposition]")
{
    DepositionEvent initial{0.,{{1.,1.,400.,0.},{1.,3.,300.,0.}},{},{}};
    DepositionEvent joined;
    joined.time_s = 1.;
    joined.coalesced_nodes = std::vector<std::size_t>{0,0};
    joined.added_mass.push_back({0,1.,1.,500.});
    HistoryObservations observations;
    observations.times_s = {.5,1.};
    std::vector<std::vector<Cell>> snapshots;
    observations.observe = [&](double,const std::vector<Cell>& cells) { snapshots.push_back(cells); };
    const auto result = evolve({initial,joined},1.,{},{},observations);
    REQUIRE(result.state == State::Converged);
    REQUIRE(snapshots.size() == 2);
    REQUIRE(snapshots[0].size() == 2);
    CHECK_THAT(snapshots[0][0].temperature_k,WithinAbs(400.,1e-12));
    CHECK_THAT(snapshots[0][1].temperature_k,WithinAbs(300.,1e-12));
    REQUIRE(snapshots[1].size() == 1);
    CHECK_THAT(snapshots[1][0].temperature_k,WithinAbs(360.,1e-12));
    CHECK_THAT(snapshots[1][0].volume_m3,WithinAbs(3.,1e-12));
    CHECK_THAT(result.final_energy_j,WithinAbs(1800.,1e-9));
    CHECK_THAT(result.energy_balance_error_j,WithinAbs(0.,1e-9));
}

TEST_CASE("Thermal contact coalescence rejects incomplete maps and exhausted work", "[ThermalDeposition]")
{
    DepositionEvent initial{0.,{{1.,1.,400.,0.},{1.,3.,300.,0.}},{},{}};
    DepositionEvent joined;
    joined.time_s = 1.;
    for (const std::vector<std::size_t>& invalid : std::vector<std::vector<std::size_t>>{{},{0},{1,1},{0,2}}) {
        joined.coalesced_nodes = invalid;
        const auto result = evolve({initial,joined},1.);
        CHECK(result.state == State::InvalidInput);
        CHECK(result.temperatures_k.empty());
    }
    joined.coalesced_nodes = std::vector<std::size_t>{0,0};
    HistoryLimits limits;
    limits.max_coalescence_work = 1;
    CHECK(evolve({initial,joined},1.,limits).state == State::BudgetExceeded);
}

TEST_CASE("Thermal joins preserve conduction and source power under node renumbering", "[ThermalDeposition]")
{
    DepositionEvent initial{0.,{{1.,1.,400.,1.},{1.,1.,300.,3.},{1.,2.,300.,0.}},
        {{0,2,1.},{1,2,1.},{0,1,0.}},{}};
    DepositionEvent joined;
    joined.coalesced_nodes = std::vector<std::size_t>{0,0,1};
    const auto result = evolve({initial,joined},1.);
    REQUIRE(result.state == State::Converged);
    REQUIRE(result.temperatures_k.size() == 2);
    CHECK_THAT(result.temperatures_k[0],WithinAbs(1004./3.,1e-8));
    CHECK_THAT(result.temperatures_k[1],WithinAbs(952./3.,1e-8));
    CHECK_THAT(result.source_energy_j,WithinAbs(4.,1e-12));
    CHECK_THAT(result.energy_balance_error_j,WithinAbs(0.,1e-8));
}

TEST_CASE("Thermal coalescence remaps existing element nodes without collapsing corners", "[ThermalDeposition]")
{
    DepositionEvent initial;
    initial.deposited_cells.assign(16,Cell{1.,1.,400.,0.});
    for (unsigned half = 0; half < 2; ++half) {
        Element element;
        for (unsigned node = 0; node < 8; ++node) element.nodes[node] = half*8+node;
        initial.added_elements.push_back(element);
    }
    DepositionEvent joined;
    joined.time_s = 1.;
    joined.coalesced_nodes = std::vector<std::size_t>{0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7};
    const auto result = evolve({initial,joined},2.);
    REQUIRE(result.state == State::Converged);
    CHECK(result.temperatures_k.size() == 8);
    CHECK_THAT(result.final_energy_j,WithinAbs(6400.,1e-9));
    joined.coalesced_nodes = std::vector<std::size_t>(16,0);
    CHECK(evolve({initial,joined},2.).state == State::InvalidInput);
}

TEST_CASE("Thermal mass activation rejects invalid additions without publishing a field", "[ThermalDeposition]")
{
    DepositionEvent event{0.,{{1.,1.,300.,0.}},{},{}};
    event.added_mass.push_back({1,1.,1.,400.});
    CHECK(evolve({event},0.).state == State::InvalidInput);
    event.added_mass[0].cell = 0;
    event.added_mass[0].volume_m3 = -1.;
    const auto result = evolve({event},0.);
    CHECK(result.state == State::InvalidInput);
    CHECK(result.temperatures_k.empty());
}

TEST_CASE("Thermal histories validate element events without requiring a subsequent solve", "[ThermalDeposition]")
{
    DepositionEvent event;
    event.deposited_cells.assign(8,Cell{1.,1.,300.,0.});
    Element element;
    for (unsigned i = 0; i < 8; ++i) element.nodes[i] = i;
    element.nodes[7] = 8;
    event.added_elements.push_back(element);
    CHECK(evolve({event},0.).state == State::InvalidInput);
    event.added_elements[0].nodes[7] = 7;
    HistoryLimits limits;
    limits.max_element_steps = 0;
    CHECK(evolve({event},1.,limits).state == State::BudgetExceeded);
    CHECK(evolve({event},0.,limits).state == State::Converged);
    event.added_elements[0].boundary = true;
    CHECK(evolve({event},0.).state == State::InvalidInput);
}

TEST_CASE("Adaptive thermal history accounts only accepted heat exchange by source", "[ThermalDeposition]")
{
    // A boundary change switches from bed heating to dock-fan cooling. Trial
    // steps are deliberately rejected; their heat must not enter the ledger.
    const std::vector<DepositionEvent> events{
        {0., {{1., 1., 300., 0.}}, {}, {{0, 2., 400., BoundaryKind::Bed}}},
        {1., {}, {}, {{0, 2., 300., BoundaryKind::DockFan}}}};
    HistoryLimits limits;
    limits.max_step_s = 1.;
    limits.max_step_temperature_error_k = .01;
    double switch_temperature = 0.;
    HistoryObservations observations;
    observations.times_s = {1.};
    observations.observe = [&](double, const std::vector<Cell>& cells) {
        switch_temperature = cells.at(0).temperature_k;
    };
    const auto result = evolve(events, 2., limits, {}, observations);
    REQUIRE(result.state == State::Converged);
    REQUIRE(result.rejected_steps > 0);
    CHECK_THAT(result.boundary_energy_by_kind_j[static_cast<std::size_t>(BoundaryKind::Bed)],
        WithinAbs(switch_temperature - 300., 1e-8));
    CHECK_THAT(result.boundary_energy_by_kind_j[static_cast<std::size_t>(BoundaryKind::DockFan)],
        WithinAbs(result.temperatures_k.at(0) - switch_temperature, 1e-8));
    double total = 0.;
    for (const double energy : result.boundary_energy_by_kind_j) total += energy;
    CHECK_THAT(total, WithinAbs(result.boundary_energy_j, 1e-9));
    CHECK_THAT(result.energy_balance_error_j, WithinAbs(0., 1e-8));
}

TEST_CASE("Thermal history rejects an unknown airflow channel even without a time step", "[ThermalDeposition]")
{
    const auto result = evolve({{0., {{1., 1., 300., 0.}}, {},
        {{0, 1., 300., static_cast<BoundaryKind>(99)}}}}, 0.);
    CHECK(result.state == State::InvalidInput);
    CHECK(result.temperatures_k.empty());
}

TEST_CASE("Thermal observations split integration and include simultaneous deposition events", "[ThermalDeposition]")
{
    const std::vector<DepositionEvent> events{
        {0.,{{1.,1.,300.,1.}},{},{}},
        {1.,{{1.,1.,400.,0.}},{},{}},
        {1.,{{1.,1.,500.,0.}},{},{}}};
    HistoryObservations observations;
    observations.times_s = {0.,.5,1.,1.5};
    std::size_t calls = 0;
    observations.observe = [&](double time,const std::vector<Cell>& cells) {
        CHECK_THAT(time,WithinAbs(observations.times_s.at(calls++),1e-15));
        REQUIRE(cells.size() == (time < 1. ? 1 : 3));
        CHECK_THAT(cells[0].temperature_k,WithinAbs(300.+time,1e-10));
        if (time >= 1.) CHECK_THAT(cells[2].temperature_k,WithinAbs(500.,1e-10));
    };
    const auto result = evolve(events,1.5,{}, {},observations);
    REQUIRE(result.state == State::Converged);
    CHECK(calls == 4);
    CHECK(result.completed_observations == 4);
    CHECK_THAT(result.energy_balance_error_j,WithinAbs(0.,1e-10));
}

TEST_CASE("Thermal observations reject invalid schedules and propagate consumer failure", "[ThermalDeposition]")
{
    const std::vector<DepositionEvent> events{{0.,{{1.,1.,300.,1.}},{},{}}};
    HistoryObservations observations{{.5},[](double,const std::vector<Cell>&) { throw std::runtime_error("cancel mechanics"); }};
    CHECK_THROWS_AS(evolve(events,1.,{}, {},observations),std::runtime_error);
    observations.observe = [](double,const std::vector<Cell>&) {};
    for (const auto& schedule : {std::vector<double>{.5,.5},std::vector<double>{-.1},std::vector<double>{2.},
        std::vector<double>{std::numeric_limits<double>::quiet_NaN()}}) {
        observations.times_s = schedule;
        const auto invalid = evolve(events,1.,{}, {},observations);
        CHECK(invalid.state == State::InvalidInput);
        CHECK(invalid.completed_observations == 0);
        CHECK(invalid.temperatures_k.empty());
    }
    observations.times_s.assign(4097,.5);
    CHECK(evolve(events,1.,{}, {},observations).state == State::BudgetExceeded);
    observations.times_s = {.5};
    auto invalid_events = events;
    invalid_events.push_back({.75,{},{},{}});
    invalid_events.push_back({.6,{},{},{}});
    const auto failed = evolve(invalid_events,1.,{}, {},observations);
    CHECK(failed.state == State::InvalidInput);
    CHECK(failed.completed_observations == 1);
    CHECK(failed.temperatures_k.empty());
}

TEST_CASE("Adaptive thermal integration retains heat across adjacent clock values", "[ThermalDeposition]")
{
    const double start = 1.;
    const double end = std::nextafter(start, 2.);
    const double duration = end - start;
    REQUIRE(start + duration * .5 == start);
    HistoryLimits limits;
    limits.max_step_temperature_error_k = 1e-5;
    // C=duration J/K and Q=1 W make this tiny interval add exactly 1 K.
    // Dropping it or treating its midpoint as unrepresentable loses real heat.
    const auto result = evolve({{start, {{duration, 1., 300., 1.}}, {}, {}}}, end, limits);
    REQUIRE(result.state == State::Converged);
    REQUIRE(result.temperatures_k.size() == 1);
    CHECK_THAT(result.temperatures_k[0], WithinAbs(301., 1e-10));
    CHECK(result.steps == 1);
    CHECK(result.solver_steps == 3);
    CHECK_THAT(result.source_energy_j, WithinAbs(duration, duration * 1e-10));
    CHECK_THAT(result.energy_balance_error_j, WithinAbs(0., duration * 1e-10));
}

TEST_CASE("Thermal time error control converges to analytical cooling without committing rejected heat", "[ThermalDeposition]")
{
    const std::vector<DepositionEvent> events{{0., {{1., 1., 400., 0.}}, {}, {{0, 1., 300.}}}};
    const auto fixed = evolve(events, 1.);
    REQUIRE(fixed.state == State::Converged);
    CHECK_FALSE(fixed.temporal_error_controlled);
    CHECK_THAT(fixed.maximum_accepted_interval_s, WithinAbs(1., 1e-12));
    CHECK_THAT(fixed.temperatures_k[0], WithinAbs(350., 1e-9));
    HistoryLimits limits;
    limits.max_step_temperature_error_k = 1e-4;
    const auto adaptive = evolve(events, 1., limits);
    REQUIRE(adaptive.state == State::Converged);
    REQUIRE(adaptive.temperatures_k.size() == 1);
    CHECK(adaptive.temporal_error_controlled);
    CHECK(adaptive.rejected_steps > 0);
    CHECK(adaptive.maximum_accepted_interval_s > 0.);
    CHECK(adaptive.maximum_accepted_interval_s < fixed.maximum_accepted_interval_s);
    CHECK(adaptive.solver_steps == 3 * (adaptive.steps + adaptive.rejected_steps));
    CHECK(adaptive.maximum_step_temperature_error_k <= limits.max_step_temperature_error_k);
    CHECK_THAT(adaptive.temperatures_k[0], WithinAbs(300. + 100. * std::exp(-1.), .03));
    CHECK_THAT(adaptive.energy_balance_error_j, WithinAbs(0., 1e-8));
    CHECK_THAT(adaptive.boundary_energy_j, WithinAbs(adaptive.temperatures_k[0] - 400., 1e-8));
    limits.max_steps = 3; // One rejected trial consumes the complete work budget.
    const auto exhausted = evolve(events, 1., limits);
    CHECK(exhausted.state == State::BudgetExceeded);
    CHECK(exhausted.temperatures_k.empty());
    CHECK_FALSE(exhausted.temporal_error_controlled);
    CHECK_THAT(exhausted.maximum_accepted_interval_s, WithinAbs(0., 1e-12));
    CHECK_THAT(exhausted.boundary_energy_j, WithinAbs(0., 1e-12));
    limits.max_step_temperature_error_k = std::numeric_limits<double>::quiet_NaN();
    CHECK(evolve(events, 1., limits).state == State::InvalidInput);
}

TEST_CASE("Adaptive thermal steps stop at airflow and mass activation events", "[ThermalDeposition]")
{
    HistoryLimits limits;
    limits.max_step_temperature_error_k = 1e-5;
    const auto result = evolve({
        {0., {{1., 1., 400., 0.}}, {}, {}},
        {.25, {}, {}, {{0, 1., 300., BoundaryKind::Exhaust}}},
        {.5, {{1., 1., 450., 0.}}, {}, {{0, 1., 300.}, {1, 1., 300.}}}}, 1., limits);
    REQUIRE(result.state == State::Converged);
    REQUIRE(result.temperatures_k.size() == 2);
    CHECK_THAT(result.temperatures_k[0], WithinAbs(300. + 100. * std::exp(-.75), .01));
    CHECK_THAT(result.temperatures_k[1], WithinAbs(300. + 150. * std::exp(-.5), .01));
    CHECK_THAT(result.deposited_energy_j, WithinAbs(850., 1e-12));
    CHECK_THAT(result.energy_balance_error_j, WithinAbs(0., 1e-8));
}

TEST_CASE("Thermal deposition retains activation order and accounts for new mass energy", "[ThermalDeposition]")
{
    const std::vector<DepositionEvent> events{
        {0., {{1., 2., 400., 0.}}, {}, {{0, 1., 300.}}},
        {2., {{1., 2., 500., 0.}}, {{0, 1, 1.}}, {}}};
    auto limits = HistoryLimits{};
    limits.max_step_s = 2.;
    const auto result = evolve(events, 4., limits);
    REQUIRE(result.state == State::Converged);
    REQUIRE(result.temperatures_k.size() == 2);
    // First cell cools alone to 350 K; second does not exist until t=2.
    // Covered convection is removed, then the pair conducts in isolation.
    CHECK_THAT(result.temperatures_k[0], WithinAbs(400., 1e-9));
    CHECK_THAT(result.temperatures_k[1], WithinAbs(450., 1e-9));
    CHECK_THAT(result.deposited_energy_j, WithinAbs(1800., 1e-9));
    CHECK_THAT(result.boundary_energy_j, WithinAbs(-100., 1e-9));
    CHECK_THAT(result.energy_balance_error_j, WithinAbs(0., 1e-9));
    CHECK(result.steps == 2);
}

TEST_CASE("Streaming deposition agrees with stored events and bounds its producer", "[ThermalDeposition]")
{
    const std::vector<DepositionEvent> events{
        {0., {{1., 1., 400., 0.}}, {}, {}},
        {1., {{1., 1., 300., 0.}}, {{0, 1, 1.}}, {}},
        {2., {}, {}, {{1, 1., 280., BoundaryKind::Exhaust}}}};
    std::size_t cursor = 0;
    const auto next = [&](DepositionEvent &event) {
        if (cursor == events.size()) return false;
        event = events[cursor++];
        return true;
    };
    const auto streamed = evolve_generated(next, 3.);
    const auto stored = evolve(events, 3.);
    REQUIRE(streamed.state == State::Converged);
    REQUIRE(streamed.temperatures_k.size() == stored.temperatures_k.size());
    for (std::size_t i = 0; i < streamed.temperatures_k.size(); ++i)
        CHECK_THAT(streamed.temperatures_k[i], WithinAbs(stored.temperatures_k[i], 1e-12));
    CHECK_THAT(streamed.energy_balance_error_j, WithinAbs(0., 1e-9));
    HistoryLimits limits;
    limits.max_events = 1;
    cursor = 0;
    const auto failed = evolve_generated(next, 3., limits);
    CHECK(failed.state == State::BudgetExceeded);
    CHECK(failed.temperatures_k.empty());
    CHECK(cursor == 2);
    CHECK(evolve_generated({}, 3.).state == State::InvalidInput);
    CHECK_THROWS_AS(evolve_generated([](DepositionEvent &) -> bool {
        throw std::runtime_error("producer cancelled");
    }, 3.), std::runtime_error);
}

TEST_CASE("Explicit airflow events split thermal steps at physical times", "[ThermalDeposition]")
{
    const auto result = evolve({{0., {{1., 1., 400., 0.}}, {}, {}},
        {.5, {}, {}, {{0, 1., 300., BoundaryKind::Exhaust}}}}, 1.);
    REQUIRE(result.state == State::Converged);
    CHECK(result.steps == 2);
    CHECK_THAT(result.temperatures_k[0], WithinAbs(300. + 100. / 1.5, 1e-9));
    CHECK_THAT(result.energy_balance_error_j, WithinAbs(0., 1e-9));
}

TEST_CASE("Thermal deposition rejects inactive references and bounded work exhaustion", "[ThermalDeposition]")
{
    const DepositionEvent first{0., {{1., 1., 400., 0.}}, {}, {}};
    CHECK(evolve({first, {-1., {}, {}, {}}}, 2.).state == State::InvalidInput);
    CHECK(evolve({first, {1., {}, {{0, 1, 1.}}, {}}}, 2.).state == State::InvalidInput);
    CHECK(evolve({first, {1., {}, {}, {{1, 1., 300.}}}}, 2.).state == State::InvalidInput);
    HistoryLimits limits;
    limits.max_cell_steps = 1;
    const auto failed = evolve({first}, 2., limits);
    CHECK(failed.state == State::BudgetExceeded);
    CHECK(failed.temperatures_k.empty());
    limits.max_events = 0;
    CHECK(evolve({first}, 2., limits).state == State::BudgetExceeded);
    CHECK_THROWS_AS(evolve({first}, 2., {}, [] { throw std::runtime_error("cancel"); }), std::runtime_error);
}
