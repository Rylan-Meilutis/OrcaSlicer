#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ThermalHomogenization.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <stdexcept>
#include <vector>

using namespace Slic3r::ThermalHistory;
using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;

namespace {
std::vector<ConductingTetrahedron> thermal_elements(const Mesh& mesh)
{
    std::vector<ConductingTetrahedron> elements;
    for (const auto& element : mesh.elements)
        elements.push_back({element.nodes, 2. * Eigen::Matrix3d::Identity(), 1000.});
    return elements;
}
} // namespace

TEST_CASE("Periodic thermal solid recovers its anisotropic tensor and heat capacity", "[ThermalHomogenization]")
{
    const auto mesh = mesh_voxels({{{0, 0, 0}, 1000., 0.}}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), 100);
    auto elements   = thermal_elements(mesh);
    Eigen::Matrix3d expected;
    expected << 2., .2, .1, .2, 3., .3, .1, .3, 4.;
    for (auto& element : elements)
        element.conductivity_w_m_k = expected;
    const auto result = homogenize_periodic_conduction(mesh.nodes_m, elements, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones());
    CHECK_THAT((result.conductivity_w_m_k - expected).norm(), WithinAbs(0., 1e-12));
    CHECK_THAT(result.volumetric_heat_capacity_j_m3_k, WithinAbs(1000., 1e-10));
    CHECK_THAT(result.solid_volume_fraction, WithinAbs(1., 1e-14));
    CHECK(result.maximum_energy_error < 1e-12);
}

TEST_CASE("Periodic thermal laminate recovers series and parallel conduction", "[ThermalHomogenization]")
{
    const auto mesh = mesh_voxels({{{0, 0, 0}, 1000., 0.}, {{1, 0, 0}, 1000., 0.}}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), 100);
    auto elements   = thermal_elements(mesh);
    for (std::size_t i = 6; i < elements.size(); ++i) {
        elements[i].conductivity_w_m_k              = 4. * Eigen::Matrix3d::Identity();
        elements[i].volumetric_heat_capacity_j_m3_k = 3000.;
    }
    const auto result        = homogenize_periodic_conduction(mesh.nodes_m, elements, Eigen::Vector3d::Zero(), {2., 1., 1.});
    Eigen::Matrix3d expected = Eigen::Matrix3d::Zero();
    expected.diagonal() << 8. / 3., 3., 3.;
    CHECK_THAT((result.conductivity_w_m_k - expected).norm(), WithinAbs(0., 1e-11));
    CHECK_THAT(result.volumetric_heat_capacity_j_m3_k, WithinAbs(2000., 1e-10));
    CHECK(result.maximum_energy_error < 1e-10);
}

TEST_CASE("Periodic thermal voids do not acquire fictitious conductivity or mass", "[ThermalHomogenization]")
{
    const auto rod    = mesh_voxels({{{0, 1, 1}, 1000., 0.}, {{1, 1, 1}, 1000., 0.}, {{2, 1, 1}, 1000., 0.}}, Eigen::Vector3d::Zero(),
                                    Eigen::Vector3d::Ones(), 100);
    const auto result = homogenize_periodic_conduction(rod.nodes_m, thermal_elements(rod), Eigen::Vector3d::Zero(), {3., 3., 3.});
    Eigen::Matrix3d expected = Eigen::Matrix3d::Zero();
    expected(0, 0)           = 2. / 9.;
    CHECK_THAT((result.conductivity_w_m_k - expected).norm(), WithinAbs(0., 1e-12));
    CHECK_THAT(result.volumetric_heat_capacity_j_m3_k, WithinAbs(1000. / 9., 1e-10));
    const auto island   = mesh_voxels({{{1, 1, 1}, 1000., 0.}}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), 100);
    const auto isolated = homogenize_periodic_conduction(island.nodes_m, thermal_elements(island), Eigen::Vector3d::Zero(), {3., 3., 3.});
    CHECK_THAT(isolated.conductivity_w_m_k.norm(), WithinAbs(0., 1e-12));
    CHECK_THAT(isolated.volumetric_heat_capacity_j_m3_k, WithinAbs(1000. / 27., 1e-10));
}

TEST_CASE("Periodic thermal admission rejects invalid geometry and material and respects cancellation", "[ThermalHomogenization]")
{
    auto mesh     = mesh_voxels({{{0, 0, 0}, 1000., 0.}}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), 100);
    auto elements = thermal_elements(mesh);
    CHECK_THROWS_AS(homogenize_periodic_conduction(mesh.nodes_m, elements, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), 4),
                    std::invalid_argument);
    CHECK_THROWS_AS(homogenize_periodic_conduction(mesh.nodes_m, elements, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), 100, {},
                                                   [] { throw std::runtime_error("cancel"); }),
                    std::runtime_error);
    elements[0].conductivity_w_m_k(0, 0) = -1.;
    CHECK_THROWS_AS(homogenize_periodic_conduction(mesh.nodes_m, elements, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones()),
                    std::invalid_argument);
    elements = thermal_elements(mesh);
    mesh.nodes_m.push_back({.5, .5, .5});
    CHECK_THROWS_AS(homogenize_periodic_conduction(mesh.nodes_m, elements, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones()),
                    std::invalid_argument);
    mesh.nodes_m.pop_back();
    mesh.nodes_m[0].y() = .1;
    CHECK_THROWS_AS(homogenize_periodic_conduction(mesh.nodes_m, elements, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones()),
                    std::invalid_argument);
}
