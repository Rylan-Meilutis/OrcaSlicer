#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ThermalHistory.hpp"
#include <Eigen/Core>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace Slic3r::ThermalHistory;
using Catch::Matchers::WithinAbs;

TEST_CASE("Thermal element operators reject invalid connectivity and nonphysical conduction", "[ThermalHistory]")
{
    const std::vector<Cell> cells(8,Cell{1.,1.,300.,0.});
    Element element;
    for (unsigned i = 0; i < 8; ++i) element.nodes[i] = i;
    element.conductance.setIdentity();
    // A row-sum sink must be an explicit boundary, not internal conduction.
    CHECK(step(cells,{}, {},1.,{}, {},{element}).state == State::InvalidInput);
    element.boundary = true;
    element.load.setConstant(300.);
    CHECK(step(cells,{}, {},1.,{}, {},{element}).state == State::Converged);
    element.load[0] = -300.;
    CHECK(step(cells,{}, {},1.,{}, {},{element}).state == State::InvalidInput);
    element.load[0] = 300.;
    element.conductance(0,0) = -1.;
    CHECK(step(cells,{}, {},1.,{}, {},{element}).state == State::InvalidInput);
    element.conductance.setIdentity();
    element.conductance(0,1) = .5;
    CHECK(step(cells,{}, {},1.,{}, {},{element}).state == State::InvalidInput);
    element.conductance.setIdentity();
    element.nodes[1] = 0;
    const auto invalid = step(cells,{}, {},1.,{}, {},{element});
    CHECK(invalid.state == State::InvalidInput);
    CHECK(invalid.temperatures_k.empty());
    element.nodes[1] = 1;
    Limits limits;
    limits.max_elements = 0;
    CHECK(step(cells,{}, {},1.,limits, {},{element}).state == State::BudgetExceeded);
}

TEST_CASE("Thermal cells conserve energy under internal conduction", "[ThermalHistory]")
{
    const auto result = step({{1., 2., 400., 0.}, {1., 2., 300., 0.}}, {{0, 1, 1.}}, {}, 1.);
    REQUIRE(result.state == State::Converged);
    CHECK_THAT(result.temperatures_k[0], WithinAbs(375., 1e-9));
    CHECK_THAT(result.temperatures_k[1], WithinAbs(325., 1e-9));
    CHECK_THAT(result.stored_energy_change_j, WithinAbs(0., 1e-9));
    CHECK_THAT(result.energy_balance_error_j, WithinAbs(0., 1e-9));
}

TEST_CASE("Small thermal regions cool independently of unrelated large heat capacities", "[ThermalHistory]")
{
    for (const double scale : {1.,1e-12,1e-20}) {
        const auto result=step({{1.,1.,400.,0.},{scale,1.,400.,0.}}, {},
            {{1,scale,300.,BoundaryKind::Ambient}},1.);
        REQUIRE(result.state==State::Converged);
        REQUIRE(result.temperatures_k.size()==2);
        CHECK_THAT(result.temperatures_k[0],WithinAbs(400.,1e-8));
        CHECK_THAT(result.temperatures_k[1],WithinAbs(350.,1e-8));
    }
}

TEST_CASE("Small connected thermal regions satisfy their local heat balance", "[ThermalHistory]")
{
    for (const double scale : {1.,1e-12,1e-20}) {
        const auto result=step({{1.,1.,400.,0.},{scale,1.,400.,0.}}, {{0,1,scale}},
            {{1,scale,300.,BoundaryKind::Ambient}},1.);
        REQUIRE(result.state==State::Converged);
        REQUIRE(result.temperatures_k.size()==2);
        const double large=(1200.+700.*scale)/(3.+2.*scale);
        CHECK_THAT(result.temperatures_k[0],WithinAbs(large,1e-8));
        CHECK_THAT(result.temperatures_k[1],WithinAbs((700.+large)/3.,1e-8));
    }
}

TEST_CASE("Thermal source and independent airflow boundaries balance deposited energy", "[ThermalHistory]")
{
    const std::vector<Boundary> boundaries{{0, 1., 300., BoundaryKind::PartFan},
        {0, 2., 310., BoundaryKind::AuxiliaryFan}, {0, 3., 320., BoundaryKind::DockFan}};
    const auto result = step({{1., 2., 400., 10.}}, {}, boundaries, 1.);
    REQUIRE(result.state == State::Converged);
    CHECK_THAT(result.temperatures_k[0], WithinAbs((800. + 10. + 300. + 620. + 960.) / 8., 1e-9));
    CHECK_THAT(result.source_energy_j, WithinAbs(10., 1e-9));
    for (const auto &boundary : boundaries)
        CHECK_THAT(result.boundary_energy_by_kind_j[static_cast<std::size_t>(boundary.kind)],
            WithinAbs(boundary.conductance_w_k *
                (boundary.ambient_temperature_k - result.temperatures_k[0]), 1e-9));
    CHECK_THAT(result.boundary_energy_by_kind_j[static_cast<std::size_t>(BoundaryKind::Bed)], WithinAbs(0., 1e-12));
    CHECK_THAT(result.energy_balance_error_j, WithinAbs(0., 1e-9));
    const auto isolated = step({{1., 2., 300., 10.}}, {}, {}, 2.);
    REQUIRE(isolated.state == State::Converged);
    CHECK_THAT(isolated.temperatures_k[0], WithinAbs(310., 1e-9));
}

TEST_CASE("Implicit thermal cooling converges toward the analytical transient", "[ThermalHistory]")
{
    const double exact = 300. + 100. * std::exp(-1.);
    const auto simulate = [](unsigned steps) {
        std::vector<Cell> cells{{1., 2., 400., 0.}};
        for (unsigned i = 0; i < steps; ++i) {
            const auto result = step(cells, {}, {{0, 1., 300.}}, 2. / steps);
            REQUIRE(result.state == State::Converged);
            cells[0].temperature_k = result.temperatures_k[0];
        }
        return cells[0].temperature_k;
    };
    const double coarse = std::abs(simulate(10) - exact);
    const double fine = std::abs(simulate(20) - exact);
    CHECK(fine < coarse * .55);
    CHECK(fine > coarse * .45);
    CHECK_THAT(simulate(200), WithinAbs(exact, .1));
}

TEST_CASE("Enclosure exhaust removes heat through an explicit air reservoir", "[ThermalHistory]")
{
    // Synthetic two-reservoir control: no calibrated fan curve or chamber size.
    // Part exchanges with enclosure air; only the air exchanges with the room.
    const std::vector<Cell> cells{{1., 2., 400., 0.}, {1., 2., 350., 0.}};
    const auto sealed = step(cells, {{0, 1, 1.}}, {}, 1.);
    const auto exhausted = step(cells, {{0, 1, 1.}}, {{1, 2., 300., BoundaryKind::Exhaust}}, 1.);
    REQUIRE(sealed.state == State::Converged);
    REQUIRE(exhausted.state == State::Converged);
    CHECK(exhausted.temperatures_k[0] < sealed.temperatures_k[0]);
    CHECK(exhausted.temperatures_k[1] < sealed.temperatures_k[1]);
    CHECK(exhausted.boundary_energy_j < 0.);
    CHECK_THAT(exhausted.energy_balance_error_j, WithinAbs(0., 1e-9));
}

TEST_CASE("Inactive space contributes no phantom thermal mass", "[ThermalHistory]")
{
    // Caller activates a new cell at its deposition temperature between steps.
    const auto before = step({{1., 1., 300., 0.}}, {}, {}, 1.);
    REQUIRE(before.state == State::Converged);
    const auto after = step({{1., 1., before.temperatures_k[0], 0.}, {1., 1., 500., 0.}}, {{0, 1, 1.}}, {}, 1.);
    REQUIRE(after.state == State::Converged);
    CHECK_THAT(after.temperatures_k[0] + after.temperatures_k[1], WithinAbs(800., 1e-9));
    CHECK_THAT(after.energy_balance_error_j, WithinAbs(0., 1e-9));
}

TEST_CASE("Thermal solving rejects invalid inputs and honors budgets and cancellation", "[ThermalHistory]")
{
    const std::vector<Cell> cells{{1., 1., 400., 0.}};
    CHECK(step(cells, {}, {}, 0.).state == State::InvalidInput);
    CHECK(step(cells, {{0, 1, 1.}}, {}, 1.).state == State::InvalidInput);
    CHECK(step(cells, {}, {{0, -1., 300.}}, 1.).state == State::InvalidInput);
    CHECK(step(cells, {}, {{0, 1., 300., static_cast<BoundaryKind>(-1)}}, 1.).state == State::InvalidInput);
    CHECK(step({{1., 1., std::numeric_limits<double>::quiet_NaN(), 0.}}, {}, {}, 1.).state == State::InvalidInput);
    Limits limits;
    limits.max_cells = 0;
    CHECK(step(cells, {}, {}, 1., limits).state == State::BudgetExceeded);
    limits = {};
    limits.max_iterations = 0;
    const auto unconverged = step(cells, {}, {{0, 1., 300.}}, 1., limits);
    CHECK(unconverged.state == State::Unconverged);
    CHECK(unconverged.temperatures_k.empty());
    CHECK_THROWS_AS(step(cells, {}, {}, 1., {}, [] { throw std::runtime_error("cancelled"); }), std::runtime_error);
}
