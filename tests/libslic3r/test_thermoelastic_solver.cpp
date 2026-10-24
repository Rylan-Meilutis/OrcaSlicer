#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include "libslic3r/ThermoelasticSolver.hpp"
#include "libslic3r/ThermoelasticField.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;

namespace {
const std::vector<Eigen::Vector3d> nodes{{0., 0., 0.}, {0.01, 0., 0.}, {0., 0.01, 0.}, {0., 0., 0.01}};
Element thermal_element()
{
    Element e{{0, 1, 2, 3}, 2.9e9, 0.35};
    e.eigenstrain.head<3>().setConstant(-0.005);
    return e;
}
std::vector<Constraint> fixed_nodes()
{
    std::vector<Constraint> c;
    for (std::size_t i = 0; i < nodes.size(); ++i)
        for (unsigned axis = 0; axis < 3; ++axis) c.push_back({i, axis, 0.});
    return c;
}
} // namespace

TEST_CASE("Deposited mechanical paths preserve infill gaps and reject unresolved material interfaces", "[ThermoelasticSolver]")
{
    const std::vector<Voxel> materials{{{0,0,0}, 2.9e9, .35}, {{0,0,0}, 1e8, .35}};
    DepositedSegment first{{0.,0.,.0002}, {.002,0.,.0002}, .0004, .0002, 0};
    auto second = first;
    second.start_m.y() = second.end_m.y() = .002;
    const Eigen::Vector3d spacing(.0001,.0001,.0001);
    const auto grid = voxelize_deposition({first, second}, materials, spacing, 0., 10000);
    REQUIRE_FALSE(grid.voxels.empty());
    for (const auto &voxel : grid.voxels) {
        const double y = grid.origin_m.y() + (voxel.index[1] + .5) * spacing.y();
        const bool in_gap = y > .0003 && y < .0017;
        CHECK_FALSE(in_gap);
        CHECK(voxel.index[2] < 2);
    }
    std::swap(first.start_m, first.end_m);
    const auto reversed = voxelize_deposition({first, second}, materials, spacing, 0., 10000);
    REQUIRE(reversed.voxels.size() == grid.voxels.size());
    for (std::size_t i = 0; i < grid.voxels.size(); ++i)
        CHECK(reversed.voxels[i].index == grid.voxels[i].index);
    auto conflict = first;
    conflict.material = 1;
    CHECK_THROWS_AS(voxelize_deposition({first, conflict}, materials, spacing, 0., 10000), std::runtime_error);
    CHECK_THROWS_AS(voxelize_deposition({first}, materials, {.001,.001,.0001}, 0., 10000), std::invalid_argument);
    CHECK_THROWS_AS(voxelize_deposition({first}, materials, spacing, 0., 1), std::runtime_error);
    CHECK_THROWS_AS(voxelize_deposition(std::vector<DepositedSegment>(20, first),
        materials, spacing, 0., 10000, {}, 1000), ResourceLimit);
    CHECK_THROWS_AS(voxelize_deposition({first}, materials, spacing, 0., 10000,
        [] { throw std::runtime_error("cancelled"); }), std::runtime_error);
}

TEST_CASE("Spatial stress lookup preserves tension and does not extrapolate into empty space", "[ThermoelasticSolver]")
{
    const auto e = thermal_element();
    const auto result = solve(nodes, {e}, fixed_nodes());
    const Field field(nodes, {e}, result);
    const auto sample = field.sample({.001, .001, .001});
    REQUIRE(sample.has_value());
    CHECK_THAT(sample->maximum_principal_pa, WithinAbs(e.young_pa * .005 / (1. - 2. * e.poisson), 1e-6));
    // Hydrostatic tension has zero von Mises stress; it is not zero tensile
    // stress and must not be presented as a bed-adhesion safety margin.
    CHECK_THAT(sample->von_mises_pa, WithinAbs(0., 1e-7));
    CHECK_THAT(sample->displacement_m.norm(), WithinAbs(0., 1e-15));
    CHECK_FALSE(field.sample({.009, .009, .009}).has_value());
    CHECK_FALSE(field.sample({.02, 0., 0.}).has_value());
    CHECK(field.sample({0., 0., 0.}).has_value());
    auto wrong = result;
    wrong.stress_pa.clear();
    CHECK_THROWS_AS(Field(nodes, {e}, wrong), std::invalid_argument);
}

TEST_CASE("Mechanical voxel meshing shares faces without welding disconnected islands", "[ThermoelasticSolver]")
{
    Voxel first{{0,0,0}, 2.9e9, .35};
    first.eigenstrain.head<3>().setConstant(-.005);
    Voxel second = first;
    second.index = {1,0,0};
    const auto adjacent = mesh_voxels({first, second}, Eigen::Vector3d::Zero(), {.001,.001,.001}, 100);
    CHECK(adjacent.nodes_m.size() == 12);
    CHECK(adjacent.elements.size() == 12);
    const auto bonded = solve_bonded_bed(adjacent.nodes_m, adjacent.elements, adjacent.bed_faces, 0.);
    CHECK_THAT(bonded.contact_area_total_m2, WithinAbs(2e-6, 1e-16));
    second.index = {1,1,0};
    const auto diagonal = mesh_voxels({first, second}, Eigen::Vector3d::Zero(), {.001,.001,.001}, 100);
    CHECK(diagonal.nodes_m.size() == 16);
    REQUIRE_NOTHROW(solve_bonded_bed(diagonal.nodes_m, diagonal.elements, diagonal.bed_faces, 0.));
    second.index = {1,0,1};
    const auto floating = mesh_voxels({first, second}, Eigen::Vector3d::Zero(), {.001,.001,.001}, 100);
    CHECK_THROWS_AS(solve_bonded_bed(floating.nodes_m, floating.elements, floating.bed_faces, 0.), std::runtime_error);
    CHECK_THROWS_AS(mesh_voxels({first, first}, Eigen::Vector3d::Zero(), {.001,.001,.001}, 100), std::invalid_argument);
    CHECK_THROWS_AS(mesh_voxels({first}, Eigen::Vector3d::Zero(), {.001,.001,.001}, 4), std::runtime_error);
}

TEST_CASE("A remote connection does not weld the touching edges of a mechanical slit", "[ThermoelasticSolver]")
{
    std::vector<Voxel> voxels;
    for (const auto &index : {std::array<int,3>{1,1,0}, {1,0,0}, {2,0,0},
                             {3,0,0}, {3,1,0}, {3,2,0}, {2,2,0}})
        voxels.push_back({index, 2.9e9, .35});
    const auto mesh = mesh_voxels(voxels, Eigen::Vector3d::Zero(), {.001,.001,.001}, 100);
    // The ends touch along this vertical edge, but have no common face.
    // Their remote connection around the opening must not close the slit.
    for (double z : {0., .001}) {
        std::size_t count = 0;
        for (const auto &node : mesh.nodes_m)
            if ((node - Eigen::Vector3d(.002,.002,z)).norm() < 1e-12) ++count;
        CHECK(count == 2);
    }
    REQUIRE_NOTHROW(solve_bonded_bed(mesh.nodes_m, mesh.elements, mesh.bed_faces, 0.));
}

TEST_CASE("Mechanical convergence requires repeated refinements of the same volume", "[ThermoelasticSolver]")
{
    const std::vector<Voxel> voxels{{{0,0,0}, 2.9e9, .35}};
    std::vector<unsigned> progress;
    const auto converged = solve_refined_voxels(voxels, Eigen::Vector3d::Zero(), {.001,.001,.001},
        1000, 8, .05, [&progress](unsigned subdivision) { progress.push_back(subdivision); });
    CHECK(converged.converged);
    CHECK(converged.subdivisions == 4);
    CHECK(progress == std::vector<unsigned>{1,2,4});
    CHECK_THAT(converged.bed.contact_area_total_m2, WithinAbs(1e-6, 1e-16));
    CHECK_THAT(converged.bed.equilibrium.elastic_energy_j, WithinAbs(0., 1e-16));
    const auto coarse = solve_refined_voxels(voxels, Eigen::Vector3d::Zero(), {.001,.001,.001}, 1000, 1, .05);
    CHECK_FALSE(coarse.converged);
    CHECK(coarse.subdivisions == 1);
    CHECK_FALSE(coarse.resource_limited);
    std::vector<Voxel> block;
    for (int z = 0; z < 2; ++z)
        for (int y = 0; y < 2; ++y)
            for (int x = 0; x < 2; ++x) block.push_back({{x,y,z}, 2.9e9, .35});
    const auto limited = solve_refined_voxels(block, Eigen::Vector3d::Zero(), {.001,.001,.001}, 27, 4, .05);
    CHECK_FALSE(limited.converged);
    CHECK(limited.resource_limited);
    CHECK_THROWS_AS(solve_refined_voxels(voxels, Eigen::Vector3d::Zero(), {.001,.001,.001},
        1000, 4, 0.), std::invalid_argument);
    CHECK_THROWS_AS(solve_refined_voxels(voxels, Eigen::Vector3d::Zero(), {.001,.001,.001},
        1000, 4, .05, {}, [] { throw std::runtime_error("cancelled"); }), std::runtime_error);
}

TEST_CASE("Separate mechanical components retain their own thermal gradient height", "[ThermoelasticSolver]")
{
    DepositionGrid grid{Eigen::Vector3d::Zero(), {.001,.001,.001},
        {{{0,0,0}, 2.9e9, .35}, {{0,0,1}, 2.9e9, .35}, {{10,0,0}, 2.9e9, .35}}};
    auto tops = component_top_z(grid);
    REQUIRE(tops.size() == 3);
    CHECK_THAT(tops[0], WithinAbs(.002, 1e-12));
    CHECK_THAT(tops[2], WithinAbs(.001, 1e-12));
    grid.voxels.push_back({{0,0,2}, 2.9e9, .35});
    tops = component_top_z(grid);
    CHECK_THAT(tops[0], WithinAbs(.003, 1e-12));
    CHECK_THAT(tops[2], WithinAbs(.001, 1e-12));
}

TEST_CASE("Unrestrained thermal contraction has zero elastic stress", "[ThermoelasticSolver]")
{
    // Six independent constraints remove rigid-body motion without restraining
    // uniform contraction. Exercise either tetrahedron orientation.
    auto e = thermal_element();
    const std::vector<Constraint> c{{0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2}};
    for (int orientation = 0; orientation < 2; ++orientation) {
        const auto result = solve(nodes, {e}, c);
        for (std::size_t i = 0; i < nodes.size(); ++i)
            CHECK_THAT((result.displacement_m[i] + 0.005 * nodes[i]).norm(), WithinAbs(0., 1e-14));
        CHECK_THAT(result.stress_pa[0].norm(), WithinAbs(0., 1e-6));
        CHECK_THAT(result.elastic_energy_j, WithinAbs(0., 1e-15));
        CHECK(result.relative_residual < 1e-12);
        std::swap(e.nodes[1], e.nodes[2]);
    }
}

TEST_CASE("Fully restrained thermal contraction matches hydrostatic elastic stress and energy", "[ThermoelasticSolver]")
{
    const auto e = thermal_element();
    const auto result = solve(nodes, {e}, fixed_nodes());
    const double sigma = e.young_pa * 0.005 / (1. - 2. * e.poisson);
    for (int i = 0; i < 3; ++i) CHECK_THAT(result.stress_pa[0][i], WithinAbs(sigma, 1e-6));
    CHECK_THAT(result.stress_pa[0].tail<3>().norm(), WithinAbs(0., 1e-8));
    const double volume = 0.01 * 0.01 * 0.01 / 6.;
    CHECK_THAT(result.elastic_energy_j, WithinAbs(0.5 * volume * 3. * sigma * 0.005, 1e-12));
    Eigen::Vector3d total = Eigen::Vector3d::Zero();
    for (const auto &force : result.reaction_n) total += force;
    CHECK_THAT(total.norm(), WithinAbs(0., 1e-10));
    CHECK(result.reaction_n[1].x() > 0.);
}

TEST_CASE("Affine displacement patch recovers engineering shear stress", "[ThermoelasticSolver]")
{
    auto e = thermal_element();
    e.eigenstrain.setZero();
    auto c = fixed_nodes();
    // u_x = gamma*y, all other components zero.
    for (auto &constraint : c)
        if (constraint.axis == 0) constraint.displacement_m = 0.001 * nodes[constraint.node].y();
    const auto result = solve(nodes, {e}, c);
    CHECK_THAT(result.stress_pa[0][3], WithinAbs(0.001 * e.young_pa / (2. * (1. + e.poisson)), 1e-7));
    CHECK_THAT(result.stress_pa[0].head<3>().norm(), WithinAbs(0., 1e-8));
}

TEST_CASE("Thermoelastic analysis rejects unsupported and invalid systems", "[ThermoelasticSolver]")
{
    auto e = thermal_element();
    CHECK_THROWS_AS(solve(nodes, {e}, {}), std::runtime_error);
    auto orphan = nodes;
    orphan.emplace_back(1., 1., 1.);
    CHECK_THROWS_AS(solve(orphan, {e}, fixed_nodes()), std::runtime_error);
    e.nodes[3] = e.nodes[2];
    CHECK_THROWS_AS(solve(nodes, {e}, fixed_nodes()), std::invalid_argument);
    e = thermal_element();
    e.poisson = 0.5;
    CHECK_THROWS_AS(solve(nodes, {e}, fixed_nodes()), std::invalid_argument);
    e = thermal_element();
    auto c = fixed_nodes();
    c.push_back({0, 0, 1.});
    CHECK_THROWS_AS(solve(nodes, {e}, c), std::invalid_argument);
    CHECK_THROWS_AS(solve(nodes, {e}, fixed_nodes(), [] { throw std::runtime_error("cancelled"); }), std::runtime_error);
}

TEST_CASE("A through-thickness eigenstrain gradient bends a spatial cantilever", "[ThermoelasticSolver]")
{
    // Six conforming tetrahedra per brick, all sharing the same body diagonal.
    // nu=0 isolates bending from Poisson contraction. Euler-Bernoulli theory
    // gives tip displacement -gradient*L^2/2 for this stress-free curvature.
    const double length = 0.02, height = 0.002, gradient = 0.5;
    std::vector<double> errors;
    // Linear tetrahedra are overly stiff in bending on very coarse meshes;
    // require four/eight elements through thickness and verify convergence.
    for (int resolution : {4, 8}) {
        const int nx = 10 * resolution, ny = resolution, nz = resolution;
        const auto id = [=](int x, int y, int z) { return std::size_t((z * (ny + 1) + y) * (nx + 1) + x); };
        std::vector<Eigen::Vector3d> mesh;
        std::vector<Constraint> constraints;
        for (int z = 0; z <= nz; ++z)
            for (int y = 0; y <= ny; ++y)
                for (int x = 0; x <= nx; ++x) {
                    mesh.emplace_back(length*x/nx, height*y/ny, height*z/nz);
                    if (x == 0)
                        for (unsigned axis = 0; axis < 3; ++axis) constraints.push_back({id(x,y,z), axis});
                }
        std::vector<Element> elements;
        for (int z = 0; z < nz; ++z)
            for (int y = 0; y < ny; ++y)
                for (int x = 0; x < nx; ++x) {
                    const std::size_t v[] = {id(x,y,z), id(x+1,y,z), id(x,y+1,z), id(x+1,y+1,z),
                        id(x,y,z+1), id(x+1,y,z+1), id(x,y+1,z+1), id(x+1,y+1,z+1)};
                    for (const auto &tet : {std::array<int,4>{0,1,3,7}, {0,3,2,7}, {0,2,6,7},
                                           {0,6,4,7}, {0,4,5,7}, {0,5,1,7}}) {
                        Element e{{v[tet[0]], v[tet[1]], v[tet[2]], v[tet[3]]}, 2.9e9, 0.};
                        double z_mean = 0.;
                        for (std::size_t n : e.nodes) z_mean += mesh[n].z() / 4.;
                        e.eigenstrain[0] = gradient * (z_mean - height / 2.);
                        elements.push_back(e);
                    }
                }
        const auto result = solve(mesh, elements, constraints);
        const double expected = -gradient * length * length / 2.;
        const double tip = result.displacement_m[id(nx,ny/2,nz/2)].z();
        CHECK(tip < 0.);
        CHECK_THAT(tip, WithinAbs(expected, 0.2 * std::abs(expected)));
        CHECK(result.relative_residual < 1e-8);
        errors.push_back(std::abs(tip - expected));
    }
    CHECK(errors.back() < errors.front());
}

TEST_CASE("Bonded bed demand uses real contact faces and balances normal reactions", "[ThermoelasticSolver]")
{
    const std::vector<Eigen::Vector3d> mesh{
        {0.,0.,0.}, {.01,0.,0.}, {0.,.01,0.}, {.01,.01,0.},
        {0.,0.,.002}, {.01,0.,.002}, {0.,.01,.002}, {.01,.01,.002}};
    std::vector<Element> elements;
    for (const auto &tet : {std::array<std::size_t,4>{0,1,3,7}, {0,3,2,7}, {0,2,6,7},
                           {0,6,4,7}, {0,4,5,7}, {0,5,1,7}}) {
        Element e{tet, 2.9e9, .35};
        double z = 0.;
        for (std::size_t n : tet) z += mesh[n].z() / 4.;
        e.eigenstrain[0] = e.eigenstrain[1] = -.005 * z / .002;
        elements.push_back(e);
    }
    // Only a sign/balance test: this deliberately coarse brick is not a
    // converged force prediction for a printed part.
    const auto result = solve_bonded_bed(mesh, elements, {{0,1,3}, {0,3,2}}, 0.);
    CHECK_THAT(result.contact_area_total_m2, WithinAbs(0.0001, 1e-15));
    CHECK_THAT(result.contact_area_m2[7], WithinAbs(0., 1e-15));
    CHECK(result.tensile_force_n > 0.);
    double integrated = 0., normal_balance = 0.;
    for (std::size_t n = 0; n < mesh.size(); ++n) {
        integrated += result.tensile_traction_pa[n] * result.contact_area_m2[n];
        normal_balance += result.equilibrium.reaction_n[n].z();
        CHECK(result.tensile_traction_pa[n] >= 0.);
    }
    CHECK_THAT(integrated, WithinAbs(result.tensile_force_n, 1e-12));
    CHECK_THAT(normal_balance, WithinAbs(0., 1e-8));
    const auto reversed = solve_bonded_bed(mesh, elements, {{3,1,0}, {2,3,0}}, 0.);
    CHECK_THAT(reversed.tensile_force_n, WithinAbs(result.tensile_force_n, 1e-12));
    CHECK_THROWS_AS(solve_bonded_bed(mesh, elements, {}, 0.), std::invalid_argument);
    CHECK_THROWS_AS(solve_bonded_bed(mesh, elements, {{0,1,3}, {3,1,0}}, 0.), std::invalid_argument);
    CHECK_THROWS_AS(solve_bonded_bed(mesh, elements, {{0,1,7}}, 0.), std::invalid_argument);
    CHECK_THROWS_AS(solve_bonded_bed(mesh, elements, {{0,1,3}}, 0.001), std::invalid_argument);
}
