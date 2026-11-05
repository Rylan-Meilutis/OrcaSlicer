#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ThermoelasticReduction.hpp"
#include "libslic3r/ThermoelasticSolver.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <vector>

using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
struct Basis
{
    std::vector<Eigen::Vector3d> nodes;
    std::vector<std::vector<InterpolationWeight>> rows;
};
// This fixture-specific interpolation does not merge or fill any geometry.
// Bed nodes remain explicit, and absent endpoints retain a fine node.
Basis shell_basis(const Mesh& mesh, int stride)
{
    using Index = std::array<int, 3>;
    std::map<Index, std::size_t> fine, coarse;
    std::vector<Index> indices;
    for (std::size_t i = 0; i < mesh.nodes_m.size(); ++i) {
        Index index;
        for (unsigned axis = 0; axis < 3; ++axis)
            index[axis] = int(std::lround(mesh.nodes_m[i][axis] * 1000.));
        fine[index] = i;
        indices.push_back(index);
    }
    Basis result;
    result.rows.resize(indices.size());
    const auto add = [&](const Index& index) {
        const auto found = coarse.emplace(index, result.nodes.size());
        if (found.second)
            result.nodes.push_back(mesh.nodes_m[fine.at(index)]);
        return found.first->second;
    };
    for (std::size_t i = 0; i < indices.size(); ++i) {
        const auto index = indices[i];
        auto left = index, right = index;
        left[0]  = index[0] / stride * stride;
        right[0] = left[0] + stride;
        if (index[2] == 0 || index[0] % stride == 0 || !fine.count(left) || !fine.count(right)) {
            result.rows[i] = {{add(index), 1.}};
        } else {
            const double fraction = double(index[0] - left[0]) / stride;
            result.rows[i]        = {{add(left), 1. - fraction}, {add(right), fraction}};
        }
    }
    return result;
}
Mesh printed_shell()
{
    std::vector<Voxel> voxels;
    for (int x = 0; x < 12; ++x)
        for (int y = 0; y < 8; ++y)
            for (int z = 0; z < 3; ++z) {
                if (z > 0 && x > 0 && x < 11 && y > 0 && y < 7)
                    continue;
                Voxel voxel{{x, y, z}, 2e9, .3};
                voxel.eigenstrain.head<3>().setConstant(-.001 * (z + .5));
                voxels.push_back(voxel);
            }
    return mesh_voxels(voxels, Eigen::Vector3d::Zero(), {.001, .001, .001}, 3000);
}
} // namespace

TEST_CASE("Occupied line bases preserve bed constraints and do not interpolate across cavity gaps", "[ThermoelasticReduction]")
{
    const auto mesh = printed_shell();
    std::vector<std::size_t> pinned;
    for (std::size_t node = 0; node < mesh.nodes_m.size(); ++node)
        if (mesh.nodes_m[node].z() == 0.) pinned.push_back(node);
    const auto basis = occupied_line_basis(mesh.nodes_m, mesh.elements, pinned, 0, 4);
    CHECK(basis.nodes_m.size() < mesh.nodes_m.size());
    for (auto node : pinned) {
        REQUIRE(basis.rows[node].size() == 1);
        CHECK_THAT(basis.rows[node][0].weight, WithinAbs(1., 0.));
    }
    for (std::size_t node = 0; node < mesh.nodes_m.size(); ++node) {
        Eigen::Vector3d reproduced = Eigen::Vector3d::Zero();
        for (const auto& entry : basis.rows[node]) {
            reproduced += entry.weight * basis.nodes_m[entry.node];
            // Upper interior Y lines contain only the two separated side walls.
            if (mesh.nodes_m[node].z() > .001 && mesh.nodes_m[node].y() > .001 && mesh.nodes_m[node].y() < .007)
                CHECK(std::abs(basis.nodes_m[entry.node].x() - mesh.nodes_m[node].x()) <= .0011);
        }
        CHECK_THAT((reproduced - mesh.nodes_m[node]).norm(), WithinAbs(0., 1e-15));
    }
    std::vector<Constraint> constraints;
    for (auto node : pinned) for (unsigned axis = 0; axis < 3; ++axis) constraints.push_back({node, axis});
    CHECK_NOTHROW(solve_reduced(mesh.nodes_m, mesh.elements, constraints, basis.nodes_m, basis.rows));
    CHECK_THROWS_AS(occupied_line_basis(mesh.nodes_m, mesh.elements, pinned, 3, 4), std::invalid_argument);
    CHECK_THROWS_AS(occupied_line_basis(mesh.nodes_m, mesh.elements, pinned, 0, 4, {}, [] { throw std::runtime_error("cancel"); }), std::runtime_error);
}

TEST_CASE("Galerkin reduction converges to a fine printed shell without filling its cavity", "[ThermoelasticReduction]")
{
    const auto mesh = printed_shell();
    std::vector<Constraint> constraints;
    for (std::size_t i = 0; i < mesh.nodes_m.size(); ++i)
        if (mesh.nodes_m[i].z() == 0. && (mesh.nodes_m[i].x() < .008 || mesh.nodes_m[i].y() < .004))
            for (unsigned axis = 0; axis < 3; ++axis)
                constraints.push_back({i, axis});
    const auto fine        = solve(mesh.nodes_m, mesh.elements, constraints);
    double previous_energy = 1e300;
    for (const int stride : {4, 2, 1}) {
        const auto basis   = shell_basis(mesh, stride);
        const auto reduced = solve_reduced(mesh.nodes_m, mesh.elements, constraints, basis.nodes, basis.rows);
        const auto& result = reduced.recovered;
        INFO("stride=" << stride << " fine nodes=" << mesh.nodes_m.size() << " basis=" << basis.nodes.size());
        CHECK(result.elastic_energy_j >= fine.elastic_energy_j * (1. - 1e-9));
        CHECK(result.elastic_energy_j <= previous_energy * (1. + 1e-9));
        previous_energy = result.elastic_energy_j;
        double error = 0., magnitude = 0.;
        for (std::size_t i = 0; i < mesh.nodes_m.size(); ++i) {
            error += (result.displacement_m[i] - fine.displacement_m[i]).squaredNorm();
            magnitude += fine.displacement_m[i].squaredNorm();
        }
        INFO("relative displacement error=" << std::sqrt(error / magnitude)
                                            << " relative energy excess=" << result.elastic_energy_j / fine.elastic_energy_j - 1.);
        if (stride == 2) {
            CHECK(basis.nodes.size() < mesh.nodes_m.size());
            CHECK(std::sqrt(error / magnitude) < .2);
            CHECK(result.elastic_energy_j / fine.elastic_energy_j - 1. < .2);
        }
        if (stride == 1) {
            CHECK_THAT(reduced.fine_free_residual_norm_n, WithinAbs(0., 1e-7));
            CHECK_THAT(result.elastic_energy_j, WithinRel(fine.elastic_energy_j, 1e-9));
            CHECK_THAT(std::sqrt(error / magnitude), WithinAbs(0., 1e-9));
            for (std::size_t i = 0; i < mesh.nodes_m.size(); ++i)
                CHECK_THAT((result.reaction_n[i] - fine.reaction_n[i]).norm(), WithinAbs(0., 1e-8));
        }
        for (const auto& c : constraints)
            CHECK_THAT(result.displacement_m[c.node][c.axis], WithinAbs(c.displacement_m, 1e-15));
    }
}

TEST_CASE("Galerkin reduction rejects invalid bases and bounded work", "[ThermoelasticReduction]")
{
    const auto mesh = printed_shell();
    auto basis      = shell_basis(mesh, 2);
    const std::vector<Constraint> fixed{{0, 0}, {0, 1}, {0, 2}};
    const std::vector<Constraint> oversized(3 * mesh.nodes_m.size() + 1, Constraint{0,0});
    CHECK_THROWS_AS(solve_reduced(mesh.nodes_m, mesh.elements, oversized, basis.nodes, basis.rows), std::invalid_argument);
    auto invalid              = basis;
    invalid.rows[0][0].weight = .9;
    CHECK_THROWS_AS(solve_reduced(mesh.nodes_m, mesh.elements, fixed, invalid.nodes, invalid.rows), std::invalid_argument);
    invalid = basis;
    invalid.nodes[0].x() += .0001;
    CHECK_THROWS_AS(solve_reduced(mesh.nodes_m, mesh.elements, fixed, invalid.nodes, invalid.rows), std::invalid_argument);
    const auto interpolated = std::find_if(basis.rows.begin(), basis.rows.end(), [](const auto& row) { return row.size() > 1; });
    REQUIRE(interpolated != basis.rows.end());
    const std::vector<Constraint> inexact{{std::size_t(interpolated - basis.rows.begin()), 0}};
    CHECK_THROWS_AS(solve_reduced(mesh.nodes_m, mesh.elements, inexact, basis.nodes, basis.rows), std::invalid_argument);
    ReductionLimits limits;
    limits.max_assembly_work = 1;
    CHECK_THROWS_AS(solve_reduced(mesh.nodes_m, mesh.elements, fixed, basis.nodes, basis.rows, limits), std::invalid_argument);
    CHECK_THROWS_AS(solve_reduced(mesh.nodes_m, mesh.elements, fixed, basis.nodes, basis.rows, {},
                                  [] { throw std::runtime_error("cancel"); }),
                    std::runtime_error);
}

TEST_CASE("Galerkin bases cannot bridge disconnected fine bodies", "[ThermoelasticReduction]")
{
    const auto mesh = mesh_voxels({{{0, 0, 0}, 1000., 0.}, {{1, 1, 0}, 1000., 0.}}, Eigen::Vector3d::Zero(), {.001, .001, .001}, 100);
    // Merge coincident basis positions intentionally; the fine mesh preserves
    // separate corner identities and no face-area contact exists.
    const auto basis = shell_basis(mesh, 1);
    CHECK_THROWS_AS(solve_reduced(mesh.nodes_m, mesh.elements, {}, basis.nodes, basis.rows), std::invalid_argument);
}

TEST_CASE("Galerkin bases cannot tie separated supports through a distant shell connection", "[ThermoelasticReduction]")
{
    const auto mesh = printed_shell();
    auto basis      = shell_basis(mesh, 1);
    std::vector<std::size_t> selected;
    for (std::size_t i = 0; i < mesh.nodes_m.size(); ++i)
        if (std::abs(mesh.nodes_m[i].x() - .006) < 1e-12 && std::abs(mesh.nodes_m[i].z() - .003) < 1e-12 &&
            (mesh.nodes_m[i].y() == 0. || std::abs(mesh.nodes_m[i].y() - .008) < 1e-12))
            selected.push_back(i);
    REQUIRE(selected.size() == 2);
    const Eigen::Vector3d midpoint = .5 * (mesh.nodes_m[selected[0]] + mesh.nodes_m[selected[1]]);
    const auto shared              = basis.nodes.size();
    basis.nodes.push_back(midpoint);
    for (auto node : selected) {
        const auto own   = basis.rows[node][0].node;
        basis.nodes[own] = 2. * mesh.nodes_m[node] - midpoint;
        basis.rows[node] = {{own, .5}, {shared, .5}};
    }
    // Affine reproduction and global component membership both hold, but this
    // shared basis has two disconnected positive-support islands across the cavity.
    CHECK_THROWS_AS(solve_reduced(mesh.nodes_m, mesh.elements, {}, basis.nodes, basis.rows), std::invalid_argument);
}

TEST_CASE("Generalized reduced reactions include unconstrained fine residual contributions", "[ThermoelasticReduction]")
{
    Voxel left{{0, 0, 0}, 1000., 0.}, right{{1, 0, 0}, 1000., 0.};
    left.eigenstrain[0]  = .001;
    right.eigenstrain[0] = .003;
    const auto mesh      = mesh_voxels({left, right}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), 100);
    Basis basis;
    std::map<std::array<int, 3>, std::size_t> masters;
    for (const auto& point : mesh.nodes_m)
        if (point.x() != 1.) {
            masters[{int(point.x()), int(point.y()), int(point.z())}] = basis.nodes.size();
            basis.nodes.push_back(point);
        }
    std::vector<Constraint> constraints;
    for (std::size_t i = 0; i < mesh.nodes_m.size(); ++i) {
        const auto& point = mesh.nodes_m[i];
        const int y = int(point.y()), z = int(point.z());
        if (point.x() == 1.)
            basis.rows.push_back({{masters.at({0, y, z}), .5}, {masters.at({2, y, z}), .5}});
        else
            basis.rows.push_back({{masters.at({int(point.x()), y, z}), 1.}});
        if (point.x() == 0.)
            for (unsigned axis = 0; axis < 3; ++axis)
                constraints.push_back({i, axis});
    }
    const auto reduced = solve_reduced(mesh.nodes_m, mesh.elements, constraints, basis.nodes, basis.rows);
    double raw = 0., generalized = 0.;
    for (const auto& force : reduced.recovered.reaction_n)
        raw += force.x();
    for (const auto& force : reduced.generalized_reaction_n)
        generalized += force.x();
    CHECK_THAT(raw, WithinAbs(-1., 1e-10));
    CHECK_THAT(generalized, WithinAbs(0., 1e-10));
    CHECK(reduced.fine_free_residual_norm_n > .1);
    CHECK(reduced.recovered.relative_residual < 1e-10);
    for (std::size_t i = 0; i < mesh.nodes_m.size(); ++i)
        if (mesh.nodes_m[i].x() == 2.)
            CHECK_THAT(reduced.recovered.displacement_m[i].x(), WithinAbs(.004, 1e-12));
    for (std::size_t i = 0; i < mesh.elements.size(); ++i)
        CHECK_THAT(reduced.recovered.stress_pa[i][0], WithinAbs(i < 6 ? 1. : -1., 1e-10));
    Basis identity;
    identity.nodes = mesh.nodes_m;
    for (std::size_t i = 0; i < mesh.nodes_m.size(); ++i)
        identity.rows.push_back({{i, 1.}});
    const auto exact = solve_reduced(mesh.nodes_m, mesh.elements, constraints, identity.nodes, identity.rows);
    CHECK_THAT(exact.fine_free_residual_norm_n, WithinAbs(0., 1e-10));
    for (const auto& stress : exact.recovered.stress_pa)
        CHECK_THAT(stress.norm(), WithinAbs(0., 1e-10));
}
