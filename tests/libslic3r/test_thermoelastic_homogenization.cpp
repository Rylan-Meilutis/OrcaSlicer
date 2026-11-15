#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticSolver.hpp"
#include "libslic3r/ThermoelasticHomogenization.hpp"
#include "libslic3r/PeriodicTetrahedralGeometry.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstddef>
#include <array>
#include <limits>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <vector>

using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
constexpr double contraction = -.005;

Mesh perforated_mesh()
{
    std::vector<Voxel> voxels;
    for (int x = 0; x < 6; ++x)
        for (int y = 0; y < 4; ++y) {
            if (x > 0 && x < 5 && y > 0 && y < 3) continue;
            Voxel voxel{{x, y, 0}, 2.e9, .35};
            voxel.eigenstrain.head<3>().setConstant(contraction);
            voxels.push_back(voxel);
        }
    return mesh_voxels(voxels, Eigen::Vector3d::Zero(), {.001,.001,.001}, 1000);
}

Eigen::Matrix3d rotation(double angle)
{
    Eigen::Matrix3d result;
    result << std::cos(angle), -std::sin(angle), 0.,
              std::sin(angle),  std::cos(angle), 0.,
              0., 0., 1.;
    return result;
}

std::size_t node_at(const Mesh &mesh, const Eigen::Vector3d &point)
{
    for (std::size_t i = 0; i < mesh.nodes_m.size(); ++i)
        if ((mesh.nodes_m[i] - point).norm() < 1e-14) return i;
    throw std::runtime_error("Missing analytical fixture anchor.");
}

Eigen::Matrix3d stress_tensor(const Eigen::Matrix<double, 6, 1> &stress)
{
    Eigen::Matrix3d result;
    result << stress[0], stress[3], stress[5],
              stress[3], stress[1], stress[4],
              stress[5], stress[4], stress[2];
    return result;
}
} // namespace

TEST_CASE("Rotating a perforated isotropic body does not invent directional free thermal contraction", "[ThermoelasticHomogenization]")
{
    const double angle = GENERATE(0., .5235987755982988, 1.5707963267948966);
    auto mesh = perforated_mesh();
    const auto origin = node_at(mesh, {0.,0.,0.});
    const auto along_x = node_at(mesh, {.006,0.,0.});
    const auto along_y = node_at(mesh, {0.,.004,0.});
    const auto transform = rotation(angle);
    for (auto &node : mesh.nodes_m) node = (transform * node).eval();
    // Exactly six scalar constraints remove rigid motion without preventing
    // the known affine free-contraction solution, including after rotation.
    const unsigned transverse = std::abs(mesh.nodes_m[along_x].x()) >=
        std::abs(mesh.nodes_m[along_x].y()) ? 1 : 0;
    std::vector<Constraint> constraints{{origin,0}, {origin,1}, {origin,2},
        {along_x, transverse, contraction * mesh.nodes_m[along_x][transverse]},
        {along_x,2}, {along_y,2}};
    const auto result = solve(mesh.nodes_m, mesh.elements, constraints);
    CHECK_THAT(result.elastic_energy_j, WithinAbs(0., 1e-16));
    for (std::size_t i = 0; i < mesh.nodes_m.size(); ++i) {
        CHECK_THAT((result.displacement_m[i] - contraction * mesh.nodes_m[i]).norm(), WithinAbs(0., 1e-13));
        CHECK_THAT(result.reaction_n[i].norm(), WithinAbs(0., 1e-8));
    }
    for (const auto &stress : result.stress_pa)
        CHECK_THAT(stress.norm(), WithinAbs(0., 1e-4));
}

TEST_CASE("Rotating a thermal body together with its restraint preserves energy and rotates stress", "[ThermoelasticHomogenization]")
{
    auto mesh = perforated_mesh();
    std::vector<Constraint> constraints;
    for (std::size_t i = 0; i < mesh.nodes_m.size(); ++i)
        if (std::abs(mesh.nodes_m[i].x()) < 1e-14)
            for (unsigned axis = 0; axis < 3; ++axis) constraints.push_back({i, axis});
    const auto reference = solve(mesh.nodes_m, mesh.elements, constraints);
    REQUIRE(reference.elastic_energy_j > 1e-12);
    const auto transform = rotation(.731);
    for (auto &node : mesh.nodes_m) node = (transform * node).eval();
    const auto rotated = solve(mesh.nodes_m, mesh.elements, constraints);
    CHECK_THAT(rotated.elastic_energy_j, WithinRel(reference.elastic_energy_j, 1e-9));
    for (std::size_t i = 0; i < mesh.nodes_m.size(); ++i) {
        CHECK_THAT((rotated.displacement_m[i] - transform * reference.displacement_m[i]).norm(), WithinAbs(0., 1e-13));
        CHECK_THAT((rotated.reaction_n[i] - transform * reference.reaction_n[i]).norm(), WithinAbs(0., 1e-8));
    }
    for (std::size_t i = 0; i < mesh.elements.size(); ++i)
        CHECK_THAT((stress_tensor(rotated.stress_pa[i]) -
            transform * stress_tensor(reference.stress_pa[i]) * transform.transpose()).norm(), WithinAbs(0., 1e-4));
}

TEST_CASE("Periodic homogenization recovers homogeneous isotropic stiffness and work", "[ThermoelasticHomogenization]")
{
    constexpr double young = 2.e9, poisson = .3;
    const auto mesh = mesh_voxels({Voxel{{0,0,0}, young, poisson}}, Eigen::Vector3d::Zero(), {.001,.001,.001}, 100);
    const auto result = homogenize_periodic(mesh.nodes_m, mesh.elements, Eigen::Vector3d::Zero(), {.001,.001,.001});
    const double mu = young / (2. * (1. + poisson));
    const double lambda = young * poisson / ((1. + poisson) * (1. - 2. * poisson));
    Eigen::Matrix<double, 6, 6> expected = Eigen::Matrix<double, 6, 6>::Zero();
    expected.topLeftCorner<3, 3>().setConstant(lambda);
    expected.diagonal().head<3>().array() += 2. * mu;
    expected.diagonal().tail<3>().setConstant(mu);
    CHECK_THAT((result.stiffness_pa - expected).norm(), WithinAbs(0., 1e-5));
    CHECK_THAT(result.solid_volume_fraction, WithinAbs(1., 1e-14));
    CHECK(result.maximum_hill_mandel_error < 1e-10);
    CHECK(result.relative_symmetry_error < 1e-10);
}

TEST_CASE("Periodic layered materials recover series and parallel stiffness without density scaling", "[ThermoelasticHomogenization]")
{
    // Equal layer thickness and nu=0 decouple the analytical normal cases.
    const auto mesh = mesh_voxels({Voxel{{0,0,0}, 2.e9, 0.}, Voxel{{1,0,0}, 1.e9, 0.}},
        Eigen::Vector3d::Zero(), {.001,.001,.001}, 100);
    const auto result = homogenize_periodic(mesh.nodes_m, mesh.elements, Eigen::Vector3d::Zero(), {.002,.001,.001});
    CHECK_THAT(result.stiffness_pa(0,0), WithinRel(2. / (1. / 2.e9 + 1. / 1.e9), 1e-10));
    CHECK_THAT(result.stiffness_pa(1,1), WithinRel(1.5e9, 1e-10));
    CHECK_THAT(result.stiffness_pa(2,2), WithinRel(1.5e9, 1e-10));
    CHECK_THAT(result.stiffness_pa(3,3), WithinRel(2. / (1. / 1.e9 + 1. / .5e9), 1e-10));
    CHECK_THAT(result.stiffness_pa(4,4), WithinRel(.75e9, 1e-10));
    CHECK(result.maximum_hill_mandel_error < 1e-10);
}

TEST_CASE("Periodic RVE preserves voids and rejects unmatched geometry and free islands", "[ThermoelasticHomogenization]")
{
    std::vector<Voxel> voxels;
    for (int x = 0; x < 3; ++x)
        for (int y = 0; y < 3; ++y)
            if (x != 1 || y != 1) voxels.push_back({{x,y,0}, 2.e9, .3});
    auto mesh = mesh_voxels(voxels, Eigen::Vector3d::Zero(), {.001,.001,.001}, 1000);
    const auto result = homogenize_periodic(mesh.nodes_m, mesh.elements, Eigen::Vector3d::Zero(), {.003,.003,.001});
    CHECK_THAT(result.solid_volume_fraction, WithinAbs(8. / 9., 1e-13));
    CHECK(result.stiffness_pa(0,0) > 0.);
    CHECK(result.maximum_hill_mandel_error < 1e-7);
    CHECK_THROWS_AS(homogenize_periodic(mesh.nodes_m, mesh.elements, Eigen::Vector3d::Zero(),
        {.003,.003,.001}, 1000, {}, [] { throw std::runtime_error("cancelled"); }), std::runtime_error);
    const auto boundary = node_at(mesh, {.003,0.,0.});
    mesh.nodes_m[boundary].y() = .00001;
    CHECK_THROWS_AS(homogenize_periodic(mesh.nodes_m, mesh.elements, Eigen::Vector3d::Zero(),
        {.003,.003,.001}), std::invalid_argument);
    const auto island = mesh_voxels({Voxel{{0,0,0}, 2.e9, .3}}, {.001,.001,.001}, {.001,.001,.001}, 100);
    CHECK_THROWS_AS(homogenize_periodic(island.nodes_m, island.elements, Eigen::Vector3d::Zero(),
        {.003,.003,.003}), std::runtime_error);
}

TEST_CASE("Periodic eigenstrain closure recovers a homogeneous free strain", "[ThermoelasticHomogenization]")
{
    auto mesh = mesh_voxels({Voxel{{0,0,0}, 2.e9, .3}}, Eigen::Vector3d::Zero(), {.001,.001,.001}, 100);
    Eigen::Matrix<double, 6, 1> strain;
    strain << -.003, -.002, .001, .0004, -.0002, .0001;
    for (auto &element : mesh.elements) element.eigenstrain = strain;
    const auto result = homogenize_periodic_eigenstrain(mesh.nodes_m, mesh.elements,
        Eigen::Vector3d::Zero(), {.001,.001,.001});
    CHECK_THAT((result.effective_eigenstrain - strain).norm(), WithinAbs(0., 1e-14));
    CHECK_THAT(result.residual_energy_density_j_m3, WithinAbs(0., 1e-9));
    CHECK_THAT((result.stress_offset_pa + result.elastic.stiffness_pa * strain).norm(), WithinAbs(0., 1e-7));
    CHECK_THROWS_AS(homogenize_periodic(mesh.nodes_m, mesh.elements,
        Eigen::Vector3d::Zero(), {.001,.001,.001}), std::invalid_argument);
    CHECK_THROWS_AS(homogenize_periodic_eigenstrain(mesh.nodes_m, mesh.elements,
        Eigen::Vector3d::Zero(), {.001,.001,.001}, 4), std::invalid_argument);
    CHECK_THROWS_AS(homogenize_periodic_eigenstrain(mesh.nodes_m, mesh.elements,
        Eigen::Vector3d::Zero(), {.001,.001,.001}, 100, {},
        [] { throw std::runtime_error("cancelled"); }), std::runtime_error);
    CHECK_THROWS_AS(homogenize_periodic_eigenstrain(mesh.nodes_m, mesh.elements,
        Eigen::Vector3d::Zero(), {.001,.001,.001}, 100,
        [](unsigned load) { if (load == 6) throw std::runtime_error("load cancelled"); }), std::runtime_error);
    mesh.elements.front().eigenstrain[0] = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(homogenize_periodic_eigenstrain(mesh.nodes_m, mesh.elements,
        Eigen::Vector3d::Zero(), {.001,.001,.001}), std::invalid_argument);
}

TEST_CASE("Periodic voids change stiffness without scaling uniform free eigenstrain", "[ThermoelasticHomogenization]")
{
    std::vector<Voxel> voxels;
    Eigen::Matrix<double, 6, 1> strain;
    strain << -.003, -.002, -.001, .0002, -.0003, .0001;
    for (int x = 0; x < 3; ++x)
        for (int y = 0; y < 3; ++y) {
            if (x == 1 && y == 1) continue;
            Voxel voxel{{x,y,0}, 2.e9, .3};
            voxel.eigenstrain = strain;
            voxels.push_back(voxel);
        }
    const auto mesh = mesh_voxels(voxels, Eigen::Vector3d::Zero(), {.001,.001,.001}, 1000);
    const auto closure = homogenize_periodic_eigenstrain(mesh.nodes_m, mesh.elements,
        Eigen::Vector3d::Zero(), {.003,.003,.001});
    CHECK_THAT(closure.elastic.solid_volume_fraction, WithinAbs(8. / 9., 1e-13));
    CHECK_THAT((closure.effective_eigenstrain - strain).norm(), WithinAbs(0., 1e-13));
    CHECK_THAT(closure.residual_energy_density_j_m3, WithinAbs(0., 1e-8));
    CHECK_THAT((closure.stress_offset_pa + closure.elastic.stiffness_pa * strain).norm(), WithinAbs(0., 1e-6));
}

TEST_CASE("Periodic layered eigenstrain retains incompatible energy and predicts mixed loading", "[ThermoelasticHomogenization]")
{
    Voxel first{{0,0,0}, 2.e9, 0.}, second{{1,0,0}, 1.e9, 0.};
    first.eigenstrain[1] = -.004;
    second.eigenstrain[1] = -.001;
    const auto mesh = mesh_voxels({first, second}, Eigen::Vector3d::Zero(), {.001,.001,.001}, 100);
    const Eigen::Vector3d maximum(.002,.001,.001);
    const auto closure = homogenize_periodic_eigenstrain(mesh.nodes_m, mesh.elements,
        Eigen::Vector3d::Zero(), maximum);
    const double effective = (2.e9 * -.004 + 1.e9 * -.001) / 3.e9;
    CHECK_THAT(closure.effective_eigenstrain[1], WithinRel(effective, 1e-10));
    const double minimum_energy = .25 * (2.e9 * std::pow(effective + .004, 2) +
        1.e9 * std::pow(effective + .001, 2));
    CHECK_THAT(closure.residual_energy_density_j_m3, WithinRel(minimum_energy, 1e-10));
    REQUIRE(closure.residual_energy_density_j_m3 > 0.);
    std::vector<std::array<std::size_t, 4>> connectivity;
    for (const auto &element : mesh.elements) connectivity.push_back(element.nodes);
    const auto geometry = validate_periodic_tetrahedral_geometry(mesh.nodes_m, connectivity,
        Eigen::Vector3d::Zero(), maximum, 100);
    Eigen::Matrix<double, 6, 1> macro;
    macro << .001, -.002, .003, .0004, -.0007, .0002;
    Eigen::Matrix3d tensor;
    tensor << macro[0], macro[3] / 2., macro[5] / 2.,
        macro[3] / 2., macro[1], macro[4] / 2.,
        macro[5] / 2., macro[4] / 2., macro[2];
    std::vector<Eigen::Vector3d> affine;
    for (const auto &node : mesh.nodes_m) affine.push_back(tensor * node);
    const auto fine = solve_affine_periodic(mesh.nodes_m, mesh.elements, geometry.node_classes, affine, 0);
    Eigen::Matrix<double, 6, 1> average = Eigen::Matrix<double, 6, 1>::Zero();
    for (std::size_t i = 0; i < mesh.elements.size(); ++i)
        average += geometry.volumes_m3[i] / geometry.cell_volume_m3 * fine.stress_pa[i];
    CHECK_THAT((average - closure.elastic.stiffness_pa * macro - closure.stress_offset_pa).norm(),
        WithinAbs(0., 1e-6));
    const double predicted_energy = .5 * macro.dot(closure.elastic.stiffness_pa * macro) +
        macro.dot(closure.stress_offset_pa) + closure.zero_macro_strain_energy_density_j_m3;
    CHECK_THAT(fine.elastic_energy_j / geometry.cell_volume_m3, WithinRel(predicted_energy, 1e-10));
}

TEST_CASE("Periodic closure freezes source geometry and loading before callbacks", "[ThermoelasticHomogenization]")
{
    auto mesh = mesh_voxels({Voxel{{0,0,0}, 2.e9, .3}}, Eigen::Vector3d::Zero(), {.001,.001,.001}, 100);
    for (auto &element : mesh.elements) element.eigenstrain.head<3>().setConstant(-.002);
    Eigen::Vector3d minimum = Eigen::Vector3d::Zero(), maximum(.001,.001,.001);
    const auto reference = homogenize_periodic_eigenstrain(mesh.nodes_m, mesh.elements, minimum, maximum);
    bool mutated = false;
    const auto frozen = homogenize_periodic_eigenstrain(mesh.nodes_m, mesh.elements, minimum, maximum, 100, {}, [&] {
        if (mutated) return;
        mutated = true;
        mesh.nodes_m.clear();
        mesh.elements.clear();
        minimum.setConstant(1.);
        maximum.setZero();
    });
    REQUIRE(mutated);
    CHECK_THAT((frozen.elastic.stiffness_pa - reference.elastic.stiffness_pa).norm(), WithinAbs(0., 1e-8));
    CHECK_THAT((frozen.effective_eigenstrain - reference.effective_eigenstrain).norm(), WithinAbs(0., 1e-14));
    CHECK_THAT((frozen.stress_offset_pa - reference.stress_offset_pa).norm(), WithinAbs(0., 1e-8));
    CHECK_THAT(frozen.zero_macro_strain_energy_density_j_m3,
        WithinRel(reference.zero_macro_strain_energy_density_j_m3, 1e-12));
}

TEST_CASE("Hybrid periodic assignment preserves explicit geometry and accounts for microenergy", "[ThermoelasticHomogenization]")
{
    Voxel left{{0,0,0}, 2.e9, 0.}, right{{1,0,0}, 1.e9, 0.};
    left.eigenstrain[1] = -.004;
    right.eigenstrain[1] = -.001;
    const auto rve = mesh_voxels({left,right}, Eigen::Vector3d::Zero(), {.001,.001,.001}, 100);
    const auto closure = homogenize_periodic_eigenstrain(rve.nodes_m, rve.elements,
        Eigen::Vector3d::Zero(), {.002,.001,.001});
    // Disconnected explicit cubes retain their gap and independent bed restraints.
    auto mesh = mesh_voxels({left, Voxel{{3,0,0}, 3.e9, .2}}, Eigen::Vector3d::Zero(), {.001,.001,.001}, 100);
    std::vector<Constraint> constraints;
    for (std::size_t i = 0; i < mesh.nodes_m.size(); ++i)
        if (mesh.nodes_m[i].z() < 1e-14)
            for (unsigned axis = 0; axis < 3; ++axis) constraints.push_back({i, axis});
    std::vector<std::optional<std::size_t>> assignments(mesh.elements.size());
    auto explicit_elements = mesh.elements;
    for (std::size_t i = 0; i < 6; ++i) {
        assignments[i] = 0;
        explicit_elements[i].elasticity_pa = closure.elastic.stiffness_pa;
        explicit_elements[i].eigenstrain = closure.effective_eigenstrain;
    }
    const auto expected = solve(mesh.nodes_m, explicit_elements, constraints);
    const auto result = solve_hybrid_periodic(mesh.nodes_m, mesh.elements, constraints, assignments, {closure});
    CHECK_THAT(result.homogenized_volume_m3, WithinRel(1e-9, 1e-12));
    CHECK_THAT(result.residual_microenergy_j, WithinRel(1e-9 * closure.residual_energy_density_j_m3, 1e-12));
    CHECK_THAT(result.total_energy_j, WithinRel(expected.elastic_energy_j + result.residual_microenergy_j, 1e-12));
    REQUIRE(result.equilibrium.displacement_m.size() == expected.displacement_m.size());
    for (std::size_t i = 0; i < expected.displacement_m.size(); ++i) {
        CHECK_THAT((result.equilibrium.displacement_m[i] - expected.displacement_m[i]).norm(), WithinAbs(0., 1e-15));
        CHECK_THAT((result.equilibrium.reaction_n[i] - expected.reaction_n[i]).norm(), WithinAbs(0., 1e-10));
    }
    CHECK_THROWS_AS(solve_hybrid_periodic(mesh.nodes_m, mesh.elements, constraints, assignments, {closure}, 4), std::invalid_argument);
    CHECK_THROWS_AS(solve_hybrid_periodic(mesh.nodes_m, mesh.elements, constraints, assignments, {closure}, 100,
        [] { throw std::runtime_error("cancelled"); }), std::runtime_error);
    auto invalid = closure;
    invalid.residual_energy_density_j_m3 = -1.;
    CHECK_THROWS_AS(solve_hybrid_periodic(mesh.nodes_m, mesh.elements, constraints, assignments, {invalid}), std::invalid_argument);
    invalid = closure;
    invalid.stress_offset_pa[0] = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(solve_hybrid_periodic(mesh.nodes_m, mesh.elements, constraints, assignments, {invalid}), std::invalid_argument);
    invalid = closure;
    invalid.zero_macro_strain_energy_density_j_m3 += 1.;
    CHECK_THROWS_AS(solve_hybrid_periodic(mesh.nodes_m, mesh.elements, constraints, assignments, {invalid}), std::invalid_argument);
    const auto frozen = solve_hybrid_periodic(mesh.nodes_m, mesh.elements, constraints, assignments, {closure}, 100,
        [&] { mesh.elements.clear(); });
    CHECK_THAT(frozen.total_energy_j, WithinRel(result.total_energy_j, 1e-12));
    mesh.elements = explicit_elements;
    assignments.front() = 1;
    CHECK_THROWS_AS(solve_hybrid_periodic(mesh.nodes_m, mesh.elements, constraints, assignments, {closure}), std::invalid_argument);
    assignments.pop_back();
    CHECK_THROWS_AS(solve_hybrid_periodic(mesh.nodes_m, mesh.elements, constraints, assignments, {closure}), std::invalid_argument);
}

TEST_CASE("A porous core hybrid coupon reports its error against resolved walls and voids",
    "[ThermoelasticHomogenization][.HybridThermalCoupon]")
{
    // Diagnostic only: repeating cells is not itself proof of scale separation
    // or convergence. No material/contact calibration or Auto ranking follows.
    const int repeats = GENERATE(1, 3);
    const int fine_xy_cells = 6 + 3 * repeats;
    const auto started = std::chrono::steady_clock::now();
    const auto checkpoint = [&] {
        if (std::chrono::steady_clock::now() - started > std::chrono::seconds(30))
            throw std::runtime_error("Hybrid coupon exceeds its diagnostic time budget.");
    };
    constexpr double young = 2.e9, poisson = .3, pitch = .001;
    const double core_maximum = (3 + 3 * repeats) * pitch;
    std::vector<Voxel> fine_voxels;
    std::array<EigenstrainHomogenizationResult, 3> closure;
    for (int z = 0; z < 3; ++z) {
        const double strain = -.001 * (z + 1);
        std::vector<Voxel> periodic_voxels;
        for (int x = 0; x < 3; ++x)
            for (int y = 0; y < 3; ++y) {
                Voxel voxel{{x, y, 0}, young, poisson};
                voxel.eigenstrain.head<3>().setConstant(strain);
                if (x != 1 || y != 1) periodic_voxels.push_back(voxel);
            }
        const auto rve = mesh_voxels(periodic_voxels, Eigen::Vector3d::Zero(), {pitch,pitch,pitch}, 499, checkpoint);
        closure[z] = homogenize_periodic_eigenstrain(rve.nodes_m, rve.elements,
            Eigen::Vector3d::Zero(), {3 * pitch,3 * pitch,pitch}, 499, {}, checkpoint);
        for (int x = 0; x < fine_xy_cells; ++x)
            for (int y = 0; y < fine_xy_cells; ++y) {
                if (x >= 3 && x < fine_xy_cells - 3 && y >= 3 && y < fine_xy_cells - 3 &&
                    (x - 3) % 3 == 1 && (y - 3) % 3 == 1) continue;
                Voxel voxel{{x, y, z}, young, poisson};
                voxel.eigenstrain.head<3>().setConstant(strain);
                fine_voxels.push_back(voxel);
            }
    }
    const auto fine_started = std::chrono::steady_clock::now();
    const auto fine = mesh_voxels(fine_voxels, Eigen::Vector3d::Zero(), {pitch,pitch,pitch}, 1100, checkpoint);
    REQUIRE(fine.nodes_m.size() <= 1100);
    const std::vector<EigenstrainHomogenizationResult> closures(closure.begin(), closure.end());
    const auto restraints = [&](const Mesh &mesh) {
        std::vector<Constraint> constraints;
        for (std::size_t i = 0; i < mesh.nodes_m.size(); ++i) {
            const auto &point = mesh.nodes_m[i];
            // Constrain only nodes on the exterior-wall bottom in both models.
            // The coarse core face has only shared wall-corner nodes: its
            // interpolation can over-restrain the interior despite this rule.
            // This is part of the measured error, not a bed-contact certificate.
            if (std::abs(point.z()) < 1e-14 &&
                (point.x() <= 3 * pitch + 1e-14 || point.x() >= core_maximum - 1e-14 ||
                 point.y() <= 3 * pitch + 1e-14 || point.y() >= core_maximum - 1e-14))
                for (unsigned axis = 0; axis < 3; ++axis) constraints.push_back({i, axis});
        }
        return constraints;
    };
    const auto fine_constraints = restraints(fine);
    const auto fine_result = solve(fine.nodes_m, fine.elements, fine_constraints, checkpoint);
    const double fine_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - fine_started).count();
    const double closure_seconds = std::chrono::duration<double>(fine_started - started).count();
    const auto measures = [](const Result &result, const std::vector<Constraint> &constraints, double offset) {
        std::array<double, 3> values{result.elastic_energy_j + offset, 0., 0.};
        for (const auto &constraint : constraints)
            if (constraint.axis == 2) values[1] += std::max(0., -result.reaction_n[constraint.node].z());
        for (const auto &displacement : result.displacement_m)
            values[2] = std::max(values[2], displacement.norm());
        return values;
    };
    const auto reference = measures(fine_result, fine_constraints, 0.);
    constexpr const char *names[]{"energy_J", "tensile_generalized_reaction_N", "maximum_nodal_displacement_m"};
    std::array<double, 3> previous{};
    bool have_previous = false;
    for (const int subdivision : {1, 2, 3}) {
        const auto macro_started = std::chrono::steady_clock::now();
        const int xy_cells = (repeats + 2) * subdivision;
        const double xy_pitch = 3 * pitch / subdivision;
        std::vector<Voxel> macro_voxels;
        for (int z = 0; z < 3; ++z)
            for (int x = 0; x < xy_cells; ++x)
                for (int y = 0; y < xy_cells; ++y) {
                    Voxel voxel{{x, y, z}, young, poisson};
                    voxel.eigenstrain.head<3>().setConstant(-.001 * (z + 1));
                    macro_voxels.push_back(voxel);
                }
        auto macro = mesh_voxels(macro_voxels, Eigen::Vector3d::Zero(), {xy_pitch,xy_pitch,pitch}, 1100, checkpoint);
        REQUIRE(macro.nodes_m.size() <= 1100);
        std::vector<std::optional<std::size_t>> assignments(macro.elements.size());
        for (std::size_t i = 0; i < macro.elements.size(); ++i) {
            const auto& element = macro.elements[i];
            Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
            for (const auto node : element.nodes) centroid += macro.nodes_m[node] / 4.;
            if (centroid.x() > 3 * pitch && centroid.x() < core_maximum &&
                centroid.y() > 3 * pitch && centroid.y() < core_maximum) {
                const auto layer = static_cast<std::size_t>(centroid.z() / pitch);
                REQUIRE(layer < closure.size());
                assignments[i] = layer;
            }
        }
        const auto macro_constraints = restraints(macro);
        const auto macro_result = solve_hybrid_periodic(macro.nodes_m, macro.elements, macro_constraints,
            assignments, closures, 1100, checkpoint);
        const double macro_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - macro_started).count();
        const auto approximation = measures(macro_result.equilibrium, macro_constraints, macro_result.residual_microenergy_j);
        std::cout << "Hybrid coupon: prescribed strain, fixed porous closure, no convergence verdict; fine_nodes="
                  << fine.nodes_m.size() << " macro_nodes=" << macro.nodes_m.size()
                  << " macro_xy_pitch_m=" << xy_pitch << " repeats_per_xy_side=" << repeats
                  << " fine_mesh_and_solve_s=" << fine_seconds << " macro_mesh_and_solve_s=" << macro_seconds
                  << " closure_and_fixture_setup_s=" << closure_seconds << '\n';
        for (std::size_t i = 0; i < reference.size(); ++i) {
            REQUIRE(std::isfinite(reference[i]));
            REQUIRE(std::isfinite(approximation[i]));
            REQUIRE(reference[i] > 0.);
            REQUIRE(approximation[i] > 0.);
            std::cout << names[i] << " fine=" << reference[i] << " macro=" << approximation[i]
                      << " relative_error=" << std::abs(approximation[i] - reference[i]) / reference[i];
            if (have_previous)
                std::cout << " relative_change_from_previous=" << std::abs(approximation[i] - previous[i]) / previous[i];
            std::cout << '\n';
        }
        previous = approximation;
        have_previous = true;
    }
    checkpoint();
}
