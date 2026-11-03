#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ThermomechanicalSimulation.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticMaterial.hpp"
#include "libslic3r/ThermoelasticPeeling.hpp"
#include "libslic3r/ThermoelasticSolver.hpp"
#include "libslic3r/ThermoelasticThermalStrain.hpp"
#include "libslic3r/ThermalHistory.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;

namespace {
ThermalMechanicalMaterial material()
{
    return {1.,1.,300.,{{250.,400.,Eigen::Matrix3d::Identity()*1e-5}}};
}
BoundaryScenario insulated()
{
    return {0.,{},0.,300.};
}
}

TEST_CASE("Deposition simulation qualifies reduced mechanics separately and abstains without a reference", "[ThermomechanicalSimulation]")
{
    DepositionGrid grid{Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), {}};
    std::vector<CellActivation> activation;
    for (int x = 0; x < 6; ++x) for (int z = 0; z < 2; ++z) {
        grid.voxels.push_back({{x, 0, z}, 1000., .3});
        activation.push_back({0., 320. + 10. * z});
    }
    const auto fine = simulate_deposition_grid(grid, activation, {material()}, {insulated()}, 0., 0.);
    REQUIRE(fine.mechanical.has_value());
    SimulationLimits limits;
    limits.reduction.emplace();
    limits.reduction->displacement_relative_tolerance = 1e-12;
    limits.reduction->energy_relative_tolerance = 1e-12;
    const auto reduced = simulate_deposition_grid(grid, activation, {material()}, {insulated()}, 0., 0., limits);
    REQUIRE(reduced.reduction.has_value());
    CHECK_FALSE(reduced.mechanical.has_value());
    CHECK(reduced.reduction->state == ReductionStudy::State::FineReference);
    REQUIRE(reduced.reduction->latest.has_value());
    for (std::size_t node = 0; node < fine.mesh.nodes_m.size(); ++node)
        CHECK_THAT((reduced.reduction->latest->recovered.displacement_m[node] - fine.mechanical->equilibrium.displacement_m[node]).norm(), WithinAbs(0., 1e-10));
    limits.reduction->strides = {4};
    const auto unresolved = simulate_deposition_grid(grid, activation, {material()}, {insulated()}, 0., 0., limits);
    CHECK(unresolved.reduction->state == ReductionStudy::State::Unreconciled);
    limits.reduction->strides = {16, 8};
    limits.reduction->displacement_relative_tolerance = 1.;
    limits.reduction->energy_relative_tolerance = 1.;
    const auto one_comparison = simulate_deposition_grid(grid, activation, {material()}, {insulated()}, 0., 0., limits);
    CHECK(one_comparison.reduction->state == ReductionStudy::State::Unreconciled);
    limits.reduction->strides = {16, 8, 4};
    const auto stabilized = simulate_deposition_grid(grid, activation, {material()}, {insulated()}, 0., 0., limits);
    // Strides16 and8 produce the same anchors on this short body, so they
    // cannot manufacture two independent refinement comparisons.
    CHECK(stabilized.reduction->state == ReductionStudy::State::Unreconciled);
    CHECK_FALSE(stabilized.mechanical.has_value());
    DepositionGrid long_grid{Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), {}};
    std::vector<CellActivation> long_activation;
    for (int x = 0; x < 24; ++x) {
        long_grid.voxels.push_back({{x, 0, 0}, 1000., .3});
        long_activation.push_back({0., 300.});
    }
    const auto stable = simulate_deposition_grid(long_grid, long_activation, {material()}, {insulated()}, 0., 0., limits);
    CHECK(stable.reduction->state == ReductionStudy::State::Stabilized);
    REQUIRE(stable.reduction->solved_strides.size() == 3);
    CHECK_FALSE(stable.mechanical.has_value());
    limits.reduction->solve.max_coarse_nodes = 1;
    const auto bounded = simulate_deposition_grid(grid, activation, {material()}, {insulated()}, 0., 0., limits);
    CHECK(bounded.reduction->state == ReductionStudy::State::BudgetExceeded);
    CHECK_FALSE(bounded.reduction->latest.has_value());
    limits.reduction->solve = {};
    limits.reduction->solve.max_fine_nodes = 1;
    const auto fine_budget = simulate_deposition_grid(grid, activation, {material()}, {insulated()}, 0., 0., limits);
    CHECK(fine_budget.reduction->state == ReductionStudy::State::BudgetExceeded);
}

TEST_CASE("Peeling spatial qualification preserves voids and material volume while rejecting unstable mechanics", "[ThermomechanicalSimulation][PeelingSpaceStudy]")
{
    DepositionGrid grid{{-10.,-20.,30.}, Eigen::Vector3d::Ones(),
        {{{0,0,0},1000.,0.}, {{1,0,0},1000.,0.}, {{3,0,0},1000.,0.}}};
    grid.voxels[2].material = 1;
    auto first = material();
    first.stress_free_temperature_k = 350.;
    auto second = first;
    second.volumetric_heat_capacity_j_m3_k = 3.;
    second.conductivity_w_m_k = 2.;
    const std::vector<CellActivation> activation(3,{0.,350.});
    const std::vector<GridReleasePatch> patches{{{1},0.,0.,300.}};
    PeelingSpaceStudyLimits limits;
    limits.time.run.simulation.thermal.max_step_temperature_error_k = .001;
    limits.time.run.simulation.thermal.max_step_s = .25;
    const auto stable = study_peeling_space(grid, activation, {first,second}, {insulated()}, 0., {0.,1.}, patches, limits);
    REQUIRE(stable.discretization_qualified.has_value());
    CHECK(stable.completed_levels == 3);
    CHECK(stable.completed_runs == 18);
    CHECK(stable.completed_observations == 60);
    const auto &qualified = *stable.discretization_qualified;
    CHECK(qualified.subdivisions == 4);
    CHECK(qualified.grid.voxels.size() == 192);
    CHECK_THAT(qualified.grid.spacing_m.x(), WithinAbs(.25,1e-14));
    CHECK_THAT((qualified.grid.origin_m-grid.origin_m).norm(), WithinAbs(0.,1e-14));
    std::size_t first_count=0, second_count=0;
    for (const auto &voxel : qualified.grid.voxels) {
        // The original empty x=2 column stays empty; no coarse solid envelope.
        const bool in_gap = voxel.index[0] >= 8 && voxel.index[0] < 12;
        CHECK_FALSE(in_gap);
        if (voxel.material == 0) ++first_count;
        else if (voxel.material == 1) ++second_count;
    }
    CHECK(first_count == 128);
    CHECK(second_count == 64);
    CHECK_THAT(qualified.history.final_snapshot.thermal.deposited_energy_j, WithinAbs(1750.,1e-8));
    CHECK_THAT(qualified.history.final_release.equilibrium.bonded.contact_area_total_m2, WithinAbs(3.,1e-10));
    CHECK(qualified.history.releases.empty());

    // Constant temperatures qualify in time, but restrained nonzero strain
    // does not have a grid-independent coarse tetrahedral solution.
    DepositionGrid strained{Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),
        {{{0,0,0},1000.,0.},{{1,0,0},1000.,0.}}};
    first.stress_free_temperature_k = 400.;
    auto strict = limits;
    strict.time.relative_tolerance = 1e-8;
    strict.time.force_tolerance_n = 1e-14;
    strict.time.energy_tolerance_j = 1e-14;
    const auto unresolved = study_peeling_space(strained, {{0.,350.},{0.,350.}}, {first}, {insulated()},
        0., {0.,1.}, {{{1},1e9,0.,300.}}, strict);
    CHECK_FALSE(unresolved.discretization_qualified.has_value());
    CHECK(unresolved.completed_levels == 2);
    // Pure vertical free contraction has an exact affine displacement on
    // every grid, so qualification must also admit a nonzero physical field.
    auto axial = first;
    axial.expansion.front().expansion_per_k.setZero();
    axial.expansion.front().expansion_per_k(2,2) = 1e-5;
    const auto contraction = study_peeling_space(strained, {{0.,350.},{0.,350.}}, {axial}, {insulated()},
        0., {0.,1.}, patches, limits);
    REQUIRE(contraction.discretization_qualified.has_value());
    const auto &contracted = contraction.discretization_qualified->history;
    CHECK(contracted.releases.empty());
    const auto &displacement = contracted.final_release.equilibrium.bonded.equilibrium.displacement_m;
    REQUIRE(displacement.size() == contracted.final_snapshot.mesh.nodes_m.size());
    for (std::size_t i=0; i<displacement.size(); ++i) {
        const Eigen::Vector3d exact{0.,0.,-.0005*contracted.final_snapshot.mesh.nodes_m[i].z()};
        CHECK_THAT((displacement[i]-exact).norm(),WithinAbs(0.,1e-10));
    }
    auto bounded = limits;
    bounded.max_history_runs = 1;
    CHECK_THROWS_AS(study_peeling_space(grid, activation, {first,second}, {insulated()},
        0., {0.,1.}, patches, bounded), ResourceLimit);
    bounded = limits;
    bounded.time.run.simulation.thermal.solve.max_cells = 24;
    first.stress_free_temperature_k = 350.;
    CHECK_THROWS_AS(study_peeling_space(grid, activation, {first,second}, {insulated()},
        0., {0.,1.}, patches, bounded), ResourceLimit);
}

TEST_CASE("Peeling time qualification rejects changing onset and requires independent refinements", "[ThermomechanicalSimulation][PeelingTimeStudy]")
{
    DepositionGrid grid{Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(),
        {{{0,0,0},1000.,0.}, {{1,0,0},1000.,0.}}};
    const std::vector<CellActivation> activation{{0.,350.},{0.,350.}};
    auto properties = material();
    properties.stress_free_temperature_k = 400.;
    const std::vector<GridReleasePatch> patches{{{1},0.,0.,300.}};
    PeelingTimeStudyLimits limits;
    limits.run.simulation.thermal.max_step_temperature_error_k = .001;
    limits.run.simulation.thermal.max_step_s = .25;
    const auto stable = study_peeling_time(grid, activation, {properties}, {insulated()}, 0., {0.,1.}, patches, limits);
    REQUIRE(stable.temporally_qualified.has_value());
    CHECK(stable.completed_runs == 6);
    CHECK(stable.completed_observations == 20);
    CHECK(stable.temporally_qualified->samples.size() == 5);
    CHECK(stable.temporally_qualified->releases.size() == 1);
    CHECK(stable.equilibrium_solves <= limits.max_total_equilibrium_solves);
    // With no strain, the same study must qualify no release, not invent a
    // fracture solely because the assumed interface resistance is zero.
    properties.stress_free_temperature_k = 350.;
    const auto relaxed = study_peeling_time(grid, activation, {properties}, {insulated()}, 0., {0.,1.}, patches, limits);
    REQUIRE(relaxed.temporally_qualified.has_value());
    CHECK(relaxed.temporally_qualified->releases.empty());
    const BoundaryScenario bed{0.,{},1.,300.};
    const auto cooled = simulate_deposition_grid(grid, activation, {properties}, {bed}, 0., 1.);
    const auto demand = bed_release_demands(cooled.mesh.nodes_m, cooled.mesh.elements,
        cooled.mesh.bed_faces, {{0,1}}, 0., {}, BedContactLimits{});
    const std::vector<GridReleasePatch> delayed{{{1},.5*demand.front().incremental_energy_release_j_m2,0.,300.}};
    limits.release_time_tolerance_s = .001;
    const auto unresolved = study_peeling_time(grid, activation, {properties}, {bed}, 0., {0.,1.}, delayed, limits);
    CHECK_FALSE(unresolved.temporally_qualified.has_value());
    // The thermal pair agreed; the second observation schedule exposed an
    // earlier release. Equal final released area is not sufficient evidence.
    CHECK(unresolved.completed_runs == 4);
    auto loose_thermal = limits;
    loose_thermal.run.simulation.thermal.max_step_s = 1.;
    loose_thermal.run.simulation.thermal.max_step_temperature_error_k = 100.;
    loose_thermal.temperature_tolerance_k = 1e-6;
    const auto unresolved_thermal = study_peeling_time(grid, activation, {properties}, {bed}, 0., {0.,1.},
        {{{1},1e9,0.,300.}}, loose_thermal);
    CHECK_FALSE(unresolved_thermal.temporally_qualified.has_value());
    CHECK(unresolved_thermal.completed_runs == 2);
    limits.max_history_runs = 5;
    CHECK_THROWS_AS(study_peeling_time(grid, activation, {properties}, {insulated()}, 0., {0.,1.}, patches, limits), ResourceLimit);
    limits.max_history_runs = 6;
    limits.max_total_observations = 19;
    CHECK_THROWS_AS(study_peeling_time(grid, activation, {properties}, {insulated()}, 0., {0.,1.}, patches, limits), ResourceLimit);
    limits.max_total_observations = 64;
    limits.max_total_equilibrium_solves = 1;
    CHECK_THROWS_AS(study_peeling_time(grid, activation, {properties}, {insulated()}, 0., {0.,1.}, patches, limits), ResourceLimit);
    limits.max_total_equilibrium_solves = 256;
    limits.run.simulation.thermal.max_step_temperature_error_k = 0.;
    CHECK_THROWS_AS(study_peeling_time(grid, activation, {properties}, {insulated()}, 0., {0.,1.}, patches, limits), std::invalid_argument);
}

TEST_CASE("Predicted patch release feeds back into later thermal snapshots", "[ThermomechanicalSimulation][PeelingFeedback]")
{
    DepositionGrid grid{Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(),
        {{{0,0,0},1000.,0.}, {{1,0,0},1000.,0.}}};
    const std::vector<CellActivation> activation{{0.,350.},{0.,350.}};
    auto properties = material();
    properties.stress_free_temperature_k = 400.;
    const BoundaryScenario bed{0.,{},1.,300.};
    // Zero fracture work deliberately isolates feedback bookkeeping; these
    // synthetic properties are not an ASA/PEI adhesion calibration.
    const std::vector<GridReleasePatch> patches{{{1},0.,0.,300.}};
    const auto result = simulate_peeling_history(grid, activation, {properties}, {bed}, 0., {0.,1.}, patches);
    REQUIRE(result.samples.size() == 2);
    CHECK(result.samples[0].new_releases == 1);
    CHECK(result.samples[0].pre_release_energy_j > result.samples[0].elastic_energy_j);
    CHECK(result.samples[1].new_releases == 0);
    REQUIRE(result.releases.size() == 1);
    CHECK(result.releases.front().voxel == 1);
    CHECK_THAT(result.releases.front().time_s, WithinAbs(0., 0.));
    CHECK_FALSE(result.final_snapshot.mechanical.has_value());
    CHECK_FALSE(result.final_snapshot.partial_mechanical.has_value());
    CHECK(result.final_release.steps.empty());
    CHECK_THAT(result.final_release.equilibrium.released_area_total_m2, WithinAbs(1.,1e-12));
    REQUIRE(result.final_snapshot.thermal.temperatures_k.size() == 2);
    CHECK_THAT(result.final_snapshot.thermal.temperatures_k[0], WithinAbs(4350./13.,1e-8));
    CHECK_THAT(result.final_snapshot.thermal.temperatures_k[1], WithinAbs(4450./13.,1e-8));
    CHECK_THAT(result.final_snapshot.thermal.energy_balance_error_j, WithinAbs(0.,1e-8));
    const auto prescribed = simulate_deposition_grid(grid, activation, {properties}, {bed}, 0., 1., {}, {}, result.releases);
    REQUIRE(prescribed.partial_mechanical.has_value());
    CHECK_THAT(result.final_release.equilibrium.bonded.equilibrium.elastic_energy_j,
        WithinAbs(prescribed.partial_mechanical->bed.bonded.equilibrium.elastic_energy_j,1e-12));
    // Extra observations without new release do not reset thermal integration.
    const auto sampled = simulate_peeling_history(grid, activation, {properties}, {bed}, 0., {0.,.5,1.}, patches);
    CHECK(sampled.releases.size() == 1);
    for (std::size_t i=0; i<2; ++i)
        CHECK_THAT(sampled.final_snapshot.thermal.temperatures_k[i],
            WithinAbs(result.final_snapshot.thermal.temperatures_k[i],1e-10));
    const auto strong = simulate_peeling_history(grid, activation, {properties}, {bed}, 0., {0.,1.}, {{{1},1e9,0.,300.}});
    CHECK(strong.releases.empty());
    for (const auto temperature : strong.final_snapshot.thermal.temperatures_k)
        CHECK_THAT(temperature, WithinAbs(330.,1e-8));
    auto initially_relaxed = properties;
    initially_relaxed.stress_free_temperature_k = 350.;
    const auto cooled = simulate_deposition_grid(grid, activation, {initially_relaxed}, {bed}, 0., 1.);
    const auto demand = bed_release_demands(cooled.mesh.nodes_m, cooled.mesh.elements,
        cooled.mesh.bed_faces, {{0,1}}, 0., {}, BedContactLimits{});
    REQUIRE(demand.front().incremental_energy_release_j_m2 > 0.);
    const auto delayed = simulate_peeling_history(grid, activation, {initially_relaxed}, {bed}, 0., {0.,1.,2.},
        {{{1},.5*demand.front().incremental_energy_release_j_m2,0.,300.}});
    REQUIRE(delayed.samples.size() == 3);
    CHECK(delayed.samples[0].new_releases == 0);
    CHECK(delayed.samples[1].new_releases == 1);
    CHECK(delayed.samples[2].new_releases == 0);
    REQUIRE(delayed.releases.size() == 1);
    CHECK_THAT(delayed.releases.front().time_s, WithinAbs(1.,1e-14));
    CHECK_THAT(delayed.final_snapshot.thermal.temperatures_k[0], WithinAbs(4170./13.,1e-8));
    CHECK_THAT(delayed.final_snapshot.thermal.temperatures_k[1], WithinAbs(4230./13.,1e-8));
    auto growing_grid = grid;
    growing_grid.voxels.push_back({{0,0,1},1000.,0.});
    const std::vector<CellActivation> growing_activation{{0.,350.},{0.,350.},{.5,350.}};
    const auto growing = simulate_peeling_history(growing_grid, growing_activation, {properties}, {bed}, 0., {0.,1.}, patches);
    CHECK(growing.releases.size() == 1);
    CHECK(growing.final_snapshot.active_voxel_indices.size() == 3);
    CHECK(growing.final_snapshot.mesh.elements.size() == 18);
    CHECK_THAT(growing.final_release.equilibrium.released_area_total_m2, WithinAbs(1.,1e-12));
    PeelingHistoryLimits limited;
    limited.max_equilibrium_solves = 1;
    CHECK_THROWS_AS(simulate_peeling_history(grid, activation, {properties}, {bed}, 0., {0.,1.}, patches, limited), ResourceLimit);
    CHECK_THROWS_AS(simulate_peeling_history(grid, activation, {properties}, {bed}, 0., {0.,0.}, patches), std::invalid_argument);
    CHECK_THROWS_AS(simulate_peeling_history(grid, activation, {properties}, {bed}, 0., {0.,1.},
        {{{1},0.,0.,300.},{{1},0.,0.,300.}}), std::invalid_argument);
    CHECK_THROWS_AS(simulate_peeling_history(grid, activation, {properties}, {bed}, 0., {0.,1.}, patches, {},
        [] { throw std::runtime_error("cancelled peeling feedback"); }), std::runtime_error);
}

TEST_CASE("Timed bed release changes heat exchange without resetting temperature or restoring adhesion", "[ThermomechanicalSimulation][BedReleaseHistory]")
{
    DepositionGrid grid{Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(),
        {{{0,0,0},1000.,0.}, {{1,0,0},1000.,0.}}};
    const std::vector<CellActivation> activation{{0.,350.},{0.,350.}};
    const BoundaryScenario bed{0.,{},1.,300.};
    const std::vector<BedReleaseEvent> release{{1.,1,0.,300.}};
    const auto result = simulate_deposition_grid(grid, activation, {material()}, {bed}, 0., 2., {}, {}, release);
    REQUIRE(result.thermal.state == Slic3r::ThermalHistory::State::Converged);
    REQUIRE(result.partial_mechanical.has_value());
    CHECK_FALSE(result.mechanical.has_value());
    CHECK(result.partial_mechanical->bed.nonpenetrating);
    // C=1, bed G=2/3, shared-face G=1. Both cells reach 330 K
    // at t=1. Then only the first retains bed exchange: solving the
    // two-cell backward-Euler equations gives 4170/13 and 4230/13 K.
    REQUIRE(result.thermal.temperatures_k.size() == 2);
    CHECK_THAT(result.thermal.temperatures_k[0], WithinAbs(4170./13., 1e-8));
    CHECK_THAT(result.thermal.temperatures_k[1], WithinAbs(4230./13., 1e-8));
    CHECK_THAT(result.thermal.deposited_energy_j, WithinAbs(700., 1e-8));
    CHECK_THAT(result.thermal.energy_balance_error_j, WithinAbs(0., 1e-8));
    CHECK_THAT(result.partial_mechanical->bed.released_area_total_m2, WithinAbs(1., 1e-12));
    CHECK_THAT(result.partial_mechanical->bed.bonded.contact_area_total_m2, WithinAbs(1., 1e-12));
    const auto before = simulate_deposition_grid(grid, activation, {material()}, {bed}, 0., .5, {}, {}, release);
    CHECK(before.mechanical.has_value());
    CHECK_FALSE(before.partial_mechanical.has_value());
    const auto at_release = simulate_deposition_grid(grid, activation, {material()}, {bed}, 0., 1., {}, {}, release);
    REQUIRE(at_release.partial_mechanical.has_value());
    CHECK_THAT(at_release.thermal.temperatures_k[1], WithinAbs(330., 1e-8));
    // Later global bed settings must not overwrite a released face's explicit
    // gap exchange. Even with zero strain/closed contact, adhesion stays lost.
    auto no_strain = material();
    no_strain.expansion.front().expansion_per_k.setZero();
    auto changed_bed = bed;
    changed_bed.time_s = 1.5;
    changed_bed.bed_coefficient_w_m2_k = 2.;
    const auto closed = simulate_deposition_grid(grid, activation, {no_strain}, {bed,changed_bed}, 0., 2., {}, {}, release);
    REQUIRE(closed.partial_mechanical.has_value());
    CHECK_THAT(closed.partial_mechanical->bed.maximum_lift_m, WithinAbs(0., 1e-14));
    CHECK_THAT(closed.partial_mechanical->bed.released_area_total_m2, WithinAbs(1., 1e-12));
    CHECK_THAT(closed.thermal.temperatures_k[0], WithinAbs(3500./11., 1e-8));
    CHECK_THAT(closed.thermal.temperatures_k[1], WithinAbs(3572./11., 1e-8));
    const auto conducting_gap = simulate_deposition_grid(grid, activation, {material()}, {bed},
        0., 2., {}, {}, {{1.,1,1.,300.}});
    REQUIRE(conducting_gap.partial_mechanical.has_value());
    for (const auto temperature : conducting_gap.thermal.temperatures_k)
        CHECK_THAT(temperature, WithinAbs(318., 1e-8));
    // Releasing a later-deposited cell uses grid identity, not thermal append ID.
    const auto reordered = simulate_deposition_grid(grid, {{1.,350.},{0.,350.}}, {material()}, {bed},
        0., 2., {}, {}, {{1.,0,0.,300.}});
    REQUIRE(reordered.partial_mechanical.has_value());
    CHECK(reordered.active_voxel_indices == (std::vector<std::size_t>{1,0}));
    for (std::size_t n=0; n<reordered.mesh.nodes_m.size(); ++n)
        if (reordered.mesh.nodes_m[n].x() == 0.)
            CHECK_THAT(reordered.partial_mechanical->bed.bonded.contact_area_m2[n], WithinAbs(0., 1e-14));
    CHECK_THROWS_AS(simulate_deposition_grid(grid, activation, {material()}, {bed}, 0., 2., {}, {},
        {{1.,2,0.,300.}}), std::invalid_argument);
    CHECK_THROWS_AS(simulate_deposition_grid(grid, {{0.,350.},{2.,350.}}, {material()}, {bed}, 0., 2., {}, {},
        release), std::invalid_argument);
    CHECK_THROWS_AS(simulate_deposition_grid(grid, activation, {material()}, {bed}, 0., 2., {}, {},
        {{1.,1,0.,300.},{2.,1,0.,300.}}), std::invalid_argument);
    CHECK_THROWS_AS(simulate_deposition_grid(grid, activation, {material()}, {bed}, 0., 2., {}, {},
        {{1.,1,-1.,300.}}), std::invalid_argument);
    CHECK_THROWS_AS(simulate_deposition_grid(grid, activation, {material()}, {bed}, 0., 2., {}, {},
        {{1.,1,0.,300.},{.5,0,0.,300.}}), std::invalid_argument);
    CHECK_THROWS_AS(simulate_deposition_grid(grid, activation, {material()}, {bed}, 0., 2., {}, {},
        {{std::numeric_limits<double>::infinity(),1,0.,300.}}), std::invalid_argument);
    auto elevated = grid;
    elevated.voxels[1].index = {0,0,1};
    CHECK_THROWS_AS(simulate_deposition_grid(elevated, activation, {material()}, {bed}, 0., 2., {}, {},
        release), std::invalid_argument);
    SimulationLimits bounded;
    bounded.thermal.max_events = 1;
    CHECK_THROWS_AS(simulate_deposition_grid(grid, activation, {material()}, {bed}, 0., 2., bounded, {},
        release), ResourceLimit);
    SimulationLimits reduced;
    reduced.reduction.emplace();
    reduced.reduction->displacement_relative_tolerance = .01;
    reduced.reduction->energy_relative_tolerance = .01;
    CHECK_THROWS_AS(simulate_deposition_grid(grid, activation, {material()}, {bed}, 0., 2., reduced, {},
        release), std::invalid_argument);
}

TEST_CASE("Directional surface exchange selects physical faces without changing the mechanical footprint", "[ThermomechanicalSimulation]")
{
    using Slic3r::ThermalHistory::BoundaryKind;
    DepositionGrid grid{Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(),
        {{{0,0,0},1000.,0.}, {{2,0,0},1000.,0.}}};
    SurfaceExchange rear{1.,300.,BoundaryKind::Exhaust};
    rear.face_mask = 1u << 3; // +Y only.
    rear.bounds = SurfaceExchange::Bounds{{0.,0.,0.},{1.,1.,1.}};
    const BoundaryScenario directional{0.,{rear},0.,300.};
    const auto result = simulate_deposition_grid(grid,{{0.,350.},{0.,350.}}, {material()}, {directional},0.,1.);
    const auto baseline = simulate_deposition_grid(grid,{{0.,350.},{0.,350.}}, {material()}, {insulated()},0.,1.);
    // One selected face: G=1/(.5+1)=2/3; C=1; implicit Euler gives 330 K.
    CHECK_THAT(result.thermal.temperatures_k[0],WithinAbs(330.,1e-8));
    CHECK_THAT(result.thermal.temperatures_k[1],WithinAbs(350.,1e-8));
    CHECK(result.mesh.bed_faces == baseline.mesh.bed_faces);
    REQUIRE(result.mesh.nodes_m.size() == baseline.mesh.nodes_m.size());
    REQUIRE(result.mesh.elements.size() == baseline.mesh.elements.size());
    for (std::size_t i=0; i<result.mesh.nodes_m.size(); ++i)
        CHECK_THAT((result.mesh.nodes_m[i]-baseline.mesh.nodes_m[i]).norm(),WithinAbs(0.,0.));
    for (std::size_t i=0; i<result.mesh.elements.size(); ++i)
        CHECK(result.mesh.elements[i].nodes == baseline.mesh.elements[i].nodes);

    auto overlap = directional;
    overlap.exposed.push_back(rear);
    overlap.exposed[0].temperature_k = 290.;
    overlap.exposed[1].temperature_k = 310.;
    const auto parallel = simulate_deposition_grid(grid,{{0.,350.},{0.,350.}}, {material()}, {overlap},0.,1.);
    // Two selected mechanisms: total h=2, shared G=1, weighted reservoir 300 K.
    CHECK_THAT(parallel.thermal.temperatures_k[0],WithinAbs(325.,1e-8));
    overlap.exposed[0].face_mask = 0;
    CHECK_THROWS_AS(simulate_deposition_grid(grid,{{0.,350.},{0.,350.}}, {material()}, {overlap},0.,1.),std::invalid_argument);
    overlap = directional;
    overlap.exposed[0].bounds->maximum_m.x() = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(simulate_deposition_grid(grid,{{0.,350.},{0.,350.}}, {material()}, {overlap},0.,1.),std::invalid_argument);
}

TEST_CASE("Directional exchange disappears when its selected face is buried", "[ThermomechanicalSimulation]")
{
    DepositionGrid grid{Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),
        {{{0,0,0},1000.,0.},{{1,0,0},1000.,0.}}};
    SurfaceExchange surface{1.,300.};
    surface.face_mask = 1u << 1; // +X of the first cell only.
    surface.bounds = SurfaceExchange::Bounds{{1.,0.,0.},{1.,1.,1.}};
    const BoundaryScenario boundary{0.,{surface},0.,300.};
    const auto result = simulate_deposition_grid(grid,{{0.,350.},{1.,330.}}, {material()}, {boundary},0.,2.);
    REQUIRE(result.thermal.state == Slic3r::ThermalHistory::State::Converged);
    CHECK_THAT(result.thermal.temperatures_k[0],WithinAbs(330.,1e-8));
    CHECK_THAT(result.thermal.temperatures_k[1],WithinAbs(330.,1e-8));
    CHECK_THAT(result.thermal.energy_balance_error_j,WithinAbs(0.,1e-8));
}

TEST_CASE("Spatial deposition conserves mass and heat across cell activation", "[ThermomechanicalSimulation]")
{
    DepositionGrid grid{Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),{{{0,0,0},1000.,0.},{{1,0,0},1000.,0.}}};
    const auto result=simulate_deposition_grid(grid,{{0.,300.},{1.,350.}},{material()},{insulated()},0.,2.);
    REQUIRE(result.thermal.state==Slic3r::ThermalHistory::State::Converged);
    REQUIRE(result.mechanical.has_value());
    REQUIRE(result.thermal.temperatures_k.size()==2);
    CHECK_THAT(result.thermal.temperatures_k[0],WithinAbs(950./3.,1e-8));
    CHECK_THAT(result.thermal.temperatures_k[1],WithinAbs(1000./3.,1e-8));
    CHECK_THAT(result.thermal.deposited_energy_j,WithinAbs(650.,1e-8));
    CHECK_THAT(result.thermal.energy_balance_error_j,WithinAbs(0.,1e-8));
    grid.voxels[1].index[0]=2; // Air gap is not a conductive connection.
    const auto gap=simulate_deposition_grid(grid,{{0.,300.},{1.,350.}},{material()},{insulated()},0.,2.);
    CHECK_THAT(gap.thermal.temperatures_k[0],WithinAbs(300.,1e-8));
    CHECK_THAT(gap.thermal.temperatures_k[1],WithinAbs(350.,1e-8));
    grid.voxels[1].index=grid.voxels[0].index;
    CHECK_THROWS_AS(simulate_deposition_grid(grid,{{0.,300.},{1.,350.}},{material()},{insulated()},0.,2.),std::invalid_argument);
}

TEST_CASE("Deposited path grids couple explicit fan exchange to mechanical snapshots", "[ThermomechanicalSimulation]")
{
    const std::vector<Voxel> prototypes{{{0,0,0},1000.,0.}};
    const auto grid=voxelize_deposition({{{0.,0.,.001},{.002,0.,.001},.001,.001,0}},prototypes,
                                      {.0005,.0005,.0005},0.,1000);
    std::vector<CellActivation> activation(grid.voxels.size(),{0.,350.});
    auto properties=material();
    properties.volumetric_heat_capacity_j_m3_k=1e6;
    properties.stress_free_temperature_k=350.;
    const BoundaryScenario fan{0.,{{100.,300.,Slic3r::ThermalHistory::BoundaryKind::PartFan}},0.,300.};
    const auto hot=simulate_deposition_grid(grid,activation,{properties},{insulated()},0.,1.);
    const auto cooled=simulate_deposition_grid(grid,activation,{properties},{fan},0.,1.);
    REQUIRE(hot.mechanical.has_value());
    REQUIRE(cooled.mechanical.has_value());
    CHECK_THAT(hot.mechanical->equilibrium.elastic_energy_j,WithinAbs(0.,1e-15));
    CHECK(cooled.mechanical->equilibrium.elastic_energy_j>0.);
    CHECK(cooled.thermal.final_energy_j<hot.thermal.final_energy_j);
    CHECK_THAT(cooled.thermal.energy_balance_error_j,WithinAbs(0.,1e-7));
}

TEST_CASE("Spatial thermal scenarios remove covered faces and retain explicit interface resistance", "[ThermomechanicalSimulation]")
{
    DepositionGrid grid{Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),{{{0,0,0},1000.,0.},{{1,0,0},1000.,0.}}};
    const BoundaryScenario convection{0.,{{1.,300.,Slic3r::ThermalHistory::BoundaryKind::Exhaust}},0.,300.};
    const auto covered=simulate_deposition_grid(grid,{{0.,350.},{0.,350.}},{material()},{convection},0.,1.);
    // Four non-bed faces, each G=A/(dx/(2*k)+1/h)=2/3 W/K.
    CHECK_THAT(covered.thermal.temperatures_k[0],WithinAbs(3450./11.,1e-8));
    grid.voxels[1].material=1;
    const auto interface=simulate_deposition_grid(grid,{{0.,300.},{0.,350.}},{material(),material()},{insulated()},1.,1.);
    CHECK_THAT(interface.thermal.temperatures_k[0],WithinAbs(312.5,1e-8));
    SimulationLimits limits;
    limits.max_schedule_entries=1;
    CHECK_THROWS_AS(simulate_deposition_grid(grid,{{0.,300.},{0.,350.}},{material(),material()},{insulated()},0.,1.,limits),ResourceLimit);
    limits.max_schedule_entries=1000000;
    limits.max_mechanical_nodes=100000000;
    CHECK_THROWS_AS(simulate_deposition_grid(grid,{{0.,300.},{0.,350.}},{material(),material()},{insulated()},0.,1.,limits),ResourceLimit);
    CHECK_THROWS_AS(simulate_deposition_grid(grid,{{0.,300.},{0.,350.}},{material(),material()},{insulated()},0.,1.,{},[]{throw std::runtime_error("cancel");}),std::runtime_error);
}

TEST_CASE("Surface heat exchange includes cell center conduction for convection and bed contact", "[ThermomechanicalSimulation]")
{
    DepositionGrid grid{Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),{{{0,0,0},1000.,0.}}};
    const BoundaryScenario all_faces{0.,{{1.,300.,Slic3r::ThermalHistory::BoundaryKind::Ambient}},1.,300.};
    // Six faces, G=2/3 each; C=1 J/K, backward Euler dt=1 gives (350+4*300)/5.
    const auto result=simulate_deposition_grid(grid,{{0.,350.}},{material()},{all_faces},0.,1.);
    REQUIRE(result.thermal.state==Slic3r::ThermalHistory::State::Converged);
    CHECK_THAT(result.thermal.temperatures_k[0],WithinAbs(310.,1e-8));
    const BoundaryScenario bed_only{0.,{},1.,300.};
    const auto bed=simulate_deposition_grid(grid,{{0.,350.}},{material()},{bed_only},0.,1.);
    CHECK_THAT(bed.thermal.temperatures_k[0],WithinAbs(330.,1e-8));
}

TEST_CASE("Simulation keeps world stiffness distinct from local expansion orientation", "[ThermomechanicalSimulation]")
{
    DepositionGrid grid{Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),{{{0,0,0},1000.,0.}}};
    auto properties=material();
    properties.expansion[0].expansion_per_k.diagonal()<<1e-5,2e-5,3e-5;
    properties.local_to_world<<0.,-1.,0.,1.,0.,0.,0.,0.,1.;
    ElasticityMatrix local=ElasticityMatrix::Identity();
    local.diagonal()<<100.,200.,300.,40.,50.,60.;
    const ElasticityMatrix world=rotate_elasticity(local,properties.local_to_world);
    properties.elasticity_pa=world; // Rotation is explicit, not applied twice.
    const auto result=simulate_deposition_grid(grid,{{0.,350.}},{properties},{insulated()},0.,0.);
    REQUIRE_FALSE(result.mesh.elements.empty());
    REQUIRE(result.mesh.elements[0].elasticity_pa.has_value());
    CHECK_THAT((*result.mesh.elements[0].elasticity_pa-world).norm(),WithinAbs(0.,1e-12));
    CHECK_THAT(result.mesh.elements[0].eigenstrain[0],WithinAbs(.001,1e-15));
    CHECK_THAT(result.mesh.elements[0].eigenstrain[1],WithinAbs(.0005,1e-15));
}

TEST_CASE("Concurrent surface exchanges share one cell conduction resistance", "[ThermomechanicalSimulation]")
{
    using Slic3r::ThermalHistory::BoundaryKind;
    DepositionGrid grid{Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),{{{0,0,0},1000.,0.}}};
    const BoundaryScenario separate{0.,{{1.,290.,BoundaryKind::PartFan},{1.,310.,BoundaryKind::Exhaust}},0.,300.};
    const BoundaryScenario combined{0.,{{2.,300.,BoundaryKind::Ambient}},0.,300.};
    const auto a=simulate_deposition_grid(grid,{{0.,350.}},{material()},{separate},0.,1.);
    const auto b=simulate_deposition_grid(grid,{{0.,350.}},{material()},{combined},0.,1.);
    REQUIRE(a.thermal.state==Slic3r::ThermalHistory::State::Converged);
    // Five exposed faces, each G=1/(.5+1/2)=1; bed is insulated.
    CHECK_THAT(a.thermal.temperatures_k[0],WithinAbs(1850./6.,1e-8));
    CHECK_THAT(a.thermal.temperatures_k[0],WithinAbs(b.thermal.temperatures_k[0],1e-8));
    CHECK_THAT(a.thermal.boundary_energy_j,WithinAbs(b.thermal.boundary_energy_j,1e-8));
    CHECK_THAT(a.thermal.boundary_energy_by_kind_j[static_cast<std::size_t>(BoundaryKind::Mixed)],
        WithinAbs(a.thermal.boundary_energy_j,1e-8));
    CHECK_THAT(a.thermal.boundary_energy_by_kind_j[static_cast<std::size_t>(BoundaryKind::Ambient)],WithinAbs(0.,1e-12));
    CHECK_THAT(b.thermal.boundary_energy_by_kind_j[static_cast<std::size_t>(BoundaryKind::Ambient)],
        WithinAbs(b.thermal.boundary_energy_j,1e-8));
}

TEST_CASE("A larger thermally generated contraction gradient increases prescribed release demand", "[ThermomechanicalSimulation]")
{
    using Slic3r::ThermalHistory::BoundaryKind;
    DepositionGrid grid{Eigen::Vector3d::Zero(),{.001,.001,.0005},{}};
    for(int x=0;x<8;++x) for(int y=0;y<2;++y) for(int z=0;z<2;++z)
        grid.voxels.push_back({{x,y,z},2.9e9,0.});
    auto properties=material();
    properties.conductivity_w_m_k=.2;
    properties.volumetric_heat_capacity_j_m3_k=1e6;
    properties.stress_free_temperature_k=350.;
    // Contraction only in the strip direction isolates the bending sign.
    properties.expansion[0].expansion_per_k.setZero();
    properties.expansion[0].expansion_per_k(0,0)=1e-4;
    const std::vector<CellActivation> activation(grid.voxels.size(),{0.,350.});
    const BoundaryScenario mild{0.,{{250.,330.,BoundaryKind::PartFan}},1000.,350.};
    const BoundaryScenario strong{0.,{{250.,310.,BoundaryKind::PartFan}},1000.,350.};
    const auto a=simulate_deposition_grid(grid,activation,{properties},{mild},0.,1.);
    const auto b=simulate_deposition_grid(grid,activation,{properties},{strong},0.,1.);
    REQUIRE(a.mechanical.has_value());
    REQUIRE(b.mechanical.has_value());
    CHECK(a.thermal.temperatures_k[1]<a.thermal.temperatures_k[0]);
    std::vector<std::size_t> retained;
    for(std::size_t i=0;i<a.mesh.bed_faces.size();++i) {
        double x=0.;
        for(const auto node : a.mesh.bed_faces[i]) x+=a.mesh.nodes_m[node].x()/3.;
        if(x<.001) retained.push_back(i);
    }
    const auto demand_a=bed_release_demands(a.mesh.nodes_m,a.mesh.elements,a.mesh.bed_faces,{retained},0.);
    const auto demand_b=bed_release_demands(b.mesh.nodes_m,b.mesh.elements,b.mesh.bed_faces,{retained},0.);
    CHECK(demand_a[0].maximum_lift_m>0.);
    CHECK(demand_a[0].incremental_energy_release_j_m2>0.);
    // Linear heat transfer doubles temperature departures; elasticity then
    // doubles displacement and quadruples energy, without fitting adhesion.
    CHECK_THAT(demand_b[0].maximum_lift_m,WithinAbs(2.*demand_a[0].maximum_lift_m,1e-10));
    CHECK_THAT(demand_b[0].incremental_energy_release_j_m2,
        WithinAbs(4.*demand_a[0].incremental_energy_release_j_m2,1e-6));
}
