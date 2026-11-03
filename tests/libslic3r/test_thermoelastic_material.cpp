#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ThermoelasticMaterial.hpp"
#include "libslic3r/ThermoelasticSolver.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;

TEST_CASE("Isotropic stiffness is invariant under material rotation", "[ThermoelasticMaterial]")
{
    ElasticityMatrix c = ElasticityMatrix::Zero();
    c.topLeftCorner<3,3>().setConstant(120.);
    c.diagonal().head<3>().array() += 80.;
    c.diagonal().tail<3>().setConstant(40.);
    const Eigen::Matrix3d q = Eigen::AngleAxisd(.73, Eigen::Vector3d(1.,2.,3.).normalized()).toRotationMatrix();
    CHECK_THAT((rotate_elasticity(c,q) - c).norm(), WithinAbs(0., 1e-10));
}

TEST_CASE("Orthotropic stiffness rotates normal and engineering shear components", "[ThermoelasticMaterial]")
{
    ElasticityMatrix c = ElasticityMatrix::Zero();
    c.diagonal() << 100.,200.,300.,40.,50.,60.;
    Eigen::Matrix3d q;
    q << 0.,-1.,0., 1.,0.,0., 0.,0.,1.;
    ElasticityMatrix expected = ElasticityMatrix::Zero();
    expected.diagonal() << 200.,100.,300.,40.,60.,50.;
    CHECK_THAT((rotate_elasticity(c,q) - expected).norm(), WithinAbs(0., 1e-12));
    const Eigen::Matrix3d arbitrary = Eigen::AngleAxisd(.73, Eigen::Vector3d(1.,2.,3.).normalized()).toRotationMatrix();
    CHECK_THAT((rotate_elasticity(rotate_elasticity(c,arbitrary),arbitrary.transpose()) - c).norm(), WithinAbs(0., 1e-10));
}

TEST_CASE("Explicit elasticity passes an affine tetrahedral strain patch", "[ThermoelasticMaterial]")
{
    const std::vector<Eigen::Vector3d> nodes{{0.,0.,0.},{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
    Element element{{0,1,2,3}, 0., 0.}; // E/nu ignored only when explicit tensor exists.
    ElasticityMatrix c = ElasticityMatrix::Identity() * 100.;
    c(0,1) = c(1,0) = 20.;
    element.elasticity_pa = c;
    element.eigenstrain[0] = .002;
    Eigen::Matrix3d gradient;
    gradient << .01,.004,.003, .002,.02,.005, .001,.007,.03;
    std::vector<Constraint> constraints;
    for (std::size_t n = 0; n < nodes.size(); ++n) {
        const Eigen::Vector3d displacement = gradient * nodes[n];
        for (unsigned axis = 0; axis < 3; ++axis) constraints.push_back({n,axis,displacement[axis]});
    }
    Eigen::Matrix<double,6,1> strain;
    strain << .01,.02,.03,.006,.012,.004;
    const Eigen::Matrix<double,6,1> elastic = strain - element.eigenstrain;
    const auto result = solve(nodes,{element},constraints);
    CHECK_THAT((result.stress_pa[0] - c * elastic).norm(), WithinAbs(0., 1e-12));
    CHECK_THAT(result.elastic_energy_j, WithinAbs(elastic.dot(c * elastic) / 12., 1e-12));
    (*element.elasticity_pa)(0,0) = -1.;
    CHECK_THROWS_AS(solve(nodes,{element},constraints), std::invalid_argument);
}

TEST_CASE("Material tensors reject nonpassive data and improper rotations", "[ThermoelasticMaterial]")
{
    ElasticityMatrix c = ElasticityMatrix::Identity();
    c(0,1) = .1;
    CHECK_THROWS_AS(validated_elasticity(c), std::invalid_argument);
    c.setIdentity(); c(0,0) = 0.;
    CHECK_THROWS_AS(validated_elasticity(c), std::invalid_argument);
    c(0,0) = -1.;
    CHECK_THROWS_AS(validated_elasticity(c), std::invalid_argument);
    c(0,0) = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(validated_elasticity(c), std::invalid_argument);
    c.setIdentity();
    Eigen::Matrix3d q = Eigen::Matrix3d::Identity();
    q(0,0) = -1.;
    CHECK_THROWS_AS(rotate_elasticity(c,q), std::invalid_argument);
    q(0,0) = 2.;
    CHECK_THROWS_AS(rotate_elasticity(c,q), std::invalid_argument);
    q(0,0) = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(rotate_elasticity(c,q), std::invalid_argument);
}
