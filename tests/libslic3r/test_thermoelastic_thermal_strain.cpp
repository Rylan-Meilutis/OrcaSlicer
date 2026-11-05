#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ThermoelasticThermalStrain.hpp"
#include "libslic3r/ThermoelasticSolver.hpp"
#include "libslic3r/ThermalHistory.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;

TEST_CASE("Trilinear thermal conversion retains interior cross terms and rejects law breakpoints", "[ThermoelasticThermalStrain]")
{
    const Eigen::Matrix3d alpha = Eigen::Vector3d(1e-5,2e-5,3e-5).asDiagonal();
    const auto orientation = Eigen::AngleAxisd(.4,Eigen::Vector3d::UnitY()).toRotationMatrix();
    const std::vector<ExpansionInterval> intervals{{250.,300.,alpha},{300.,400.,2.*alpha}};
    std::array<double,8> temperatures;
    for (unsigned corner = 0; corner < 8; ++corner)
        temperatures[corner] = 270.+((corner & 3) == 3 ? 8. : 0.)+(corner == 7 ? 5. : 0.);
    const auto field = trilinear_thermal_eigenstrain(intervals,orientation,375.,temperatures);
    const Eigen::Vector3d p(.2,.4,.7);
    Eigen::Matrix<double,6,1> interpolated = Eigen::Matrix<double,6,1>::Zero();
    for (unsigned corner = 0; corner < 8; ++corner) {
        double weight = 1.;
        for (unsigned axis = 0; axis < 3; ++axis)
            weight *= corner & (1u << axis) ? p[axis] : 1.-p[axis];
        interpolated += weight*field.at_corners[corner];
    }
    const auto exact = thermal_eigenstrain(intervals,orientation,375.,270.+8.*p.x()*p.y()+5.*p.prod());
    CHECK_THAT((interpolated-exact).norm(),WithinAbs(0.,1e-15));
    temperatures[7] = 310.;
    CHECK_THROWS_AS(trilinear_thermal_eigenstrain(intervals,orientation,375.,temperatures),std::invalid_argument);
    temperatures[7] = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(trilinear_thermal_eigenstrain(intervals,orientation,375.,temperatures),std::invalid_argument);
}

TEST_CASE("Thermal strain integrates temperature intervals and rotates tensor shear", "[ThermoelasticThermalStrain]")
{
    Eigen::Matrix3d alpha = Eigen::Matrix3d::Zero();
    alpha.diagonal() << 1e-5,2e-5,3e-5;
    const std::vector<ExpansionInterval> intervals{{250.,300.,alpha},{300.,350.,2.*alpha}};
    const auto strain = thermal_eigenstrain(intervals, Eigen::Matrix3d::Identity(), 325.,275.);
    CHECK_THAT(strain[0], WithinAbs(-75e-5,1e-15));
    CHECK_THAT(strain[1], WithinAbs(-150e-5,1e-15));
    CHECK_THAT((thermal_eigenstrain(intervals,Eigen::Matrix3d::Identity(),275.,325.)+strain).norm(),WithinAbs(0.,1e-15));
    const Eigen::Matrix3d q = Eigen::AngleAxisd(.7,Eigen::Vector3d::UnitZ()).toRotationMatrix();
    const Eigen::Matrix3d expected = q * (-75. * alpha) * q.transpose();
    const auto rotated = thermal_eigenstrain(intervals,q,325.,275.);
    CHECK_THAT(rotated[3],WithinAbs(2.*expected(0,1),1e-15));
    CHECK_THAT(rotated[0],WithinAbs(expected(0,0),1e-15));
    CHECK_THAT(thermal_eigenstrain(intervals,q,300.,300.).norm(),WithinAbs(0.,1e-15));
}

TEST_CASE("Thermal strain rejects gaps extrapolation and nonphysical orientation", "[ThermoelasticThermalStrain]")
{
    const Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
    const auto alpha = (1e-5 * identity).eval();
    CHECK_THROWS_AS(thermal_eigenstrain({},identity,300.,290.),std::invalid_argument);
    CHECK_THROWS_AS(thermal_eigenstrain({{280.,320.,alpha}},identity,350.,300.),std::invalid_argument);
    CHECK_THROWS_AS(thermal_eigenstrain({{280.,290.,alpha},{300.,320.,alpha}},identity,280.,320.),std::invalid_argument);
    CHECK_THROWS_AS(thermal_eigenstrain({{280.,310.,alpha},{300.,320.,alpha}},identity,280.,320.),std::invalid_argument);
    CHECK_THROWS_AS(thermal_eigenstrain({{280.,320.,alpha}},2.*identity,300.,290.),std::invalid_argument);
    Eigen::Matrix3d asymmetric = alpha;
    asymmetric(0,1) = 1e-5;
    CHECK_THROWS_AS(thermal_eigenstrain({{280.,320.,asymmetric}},identity,300.,290.),std::invalid_argument);
}

TEST_CASE("Affine temperature fields preserve rotated expansion and world coordinates", "[ThermoelasticThermalStrain]")
{
    Eigen::Matrix3d alpha = Eigen::Matrix3d::Zero();
    alpha.diagonal() << 1e-5,2e-5,3e-5;
    const std::vector<ExpansionInterval> intervals{{250.,300.,alpha},{300.,400.,2.*alpha}};
    const auto q = Eigen::AngleAxisd(.7,Eigen::Vector3d::UnitY()).toRotationMatrix();
    const Eigen::Vector3d lo(.1,.2,.3), hi(.2,.3,.4), anchor(.15,.25,.35), gradient(10.,-20.,30.);
    const auto field = linear_thermal_eigenstrain(intervals,q,375.,anchor,275.,gradient,lo,hi);
    // The stress-free reference crosses the expansion-law change, while the
    // entire current spatial field remains in the lower-temperature interval.
    for (unsigned corner = 0; corner < 8; ++corner) {
        Eigen::Vector3d point;
        for (unsigned axis = 0; axis < 3; ++axis) point[axis] = corner&(1u<<axis) ? hi[axis] : lo[axis];
        const auto exact = thermal_eigenstrain(intervals,q,375.,275.+gradient.dot(point-anchor));
        CHECK_THAT((field.at_world_zero+field.gradient_per_m*point-exact).norm(),WithinAbs(0.,1e-15));
    }
    const auto uniform = linear_thermal_eigenstrain(intervals,q,375.,anchor,300.,Eigen::Vector3d::Zero(),lo,hi);
    CHECK_THAT(uniform.gradient_per_m.norm(),WithinAbs(0.,1e-15));
    CHECK_THAT((uniform.at_world_zero-thermal_eigenstrain(intervals,q,375.,300.)).norm(),WithinAbs(0.,1e-15));
    CHECK_THROWS_AS(linear_thermal_eigenstrain(intervals,q,375.,anchor,300.,gradient,lo,hi),std::invalid_argument);
    CHECK_THROWS_AS(linear_thermal_eigenstrain(intervals,q,375.,anchor,251.,gradient,lo,hi),std::invalid_argument);
    CHECK_THROWS_AS(linear_thermal_eigenstrain(intervals,q,375.,hi+Eigen::Vector3d::Ones(),275.,gradient,lo,hi),std::invalid_argument);
    CHECK_THROWS_AS(linear_thermal_eigenstrain(intervals,q,375.,anchor,275.,gradient,hi,lo),std::invalid_argument);
    auto invalid = gradient;
    invalid.z() = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(linear_thermal_eigenstrain(intervals,q,375.,anchor,275.,invalid,lo,hi),std::invalid_argument);
    // Adjacent intervals with the same tensor need no subdivision.
    CHECK_NOTHROW(linear_thermal_eigenstrain({{250.,300.,alpha},{300.,400.,alpha}},q,375.,anchor,300.,gradient,lo,hi));
}

TEST_CASE("Explicit fan and exhaust heat exchange changes constrained thermal stress", "[ThermoelasticThermalStrain]")
{
    namespace Heat = Slic3r::ThermalHistory;
    // Synthetic constant properties verify coupling, not a PWM-to-airflow law.
    const std::vector<Heat::Cell> cells{{1.,1.,350.,0.}};
    const auto fan = Heat::step(cells,{},{{0,1.,300.,Heat::BoundaryKind::PartFan}},1.);
    const auto both = Heat::step(cells,{},{{0,1.,300.,Heat::BoundaryKind::PartFan},
                                        {0,1.,300.,Heat::BoundaryKind::Exhaust}},1.);
    REQUIRE(fan.state == Heat::State::Converged);
    REQUIRE(both.state == Heat::State::Converged);
    CHECK_THAT(fan.temperatures_k[0],WithinAbs(325.,1e-8));
    CHECK_THAT(both.temperatures_k[0],WithinAbs(950./3.,1e-8));
    const std::vector<Eigen::Vector3d> nodes{{0.,0.,0.},{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
    std::vector<Constraint> fixed;
    for (std::size_t n=0;n<nodes.size();++n)
        for(unsigned axis=0;axis<3;++axis) fixed.push_back({n,axis,0.});
    Element element{{0,1,2,3},0.,0.};
    Eigen::Matrix<double,6,6> stiffness = Eigen::Matrix<double,6,6>::Identity()*1000.;
    stiffness(1,1)=2000.;
    element.elasticity_pa=stiffness;
    const std::vector<ExpansionInterval> expansion{{250.,350.,Eigen::Matrix3d::Identity()*1e-5}};
    element.eigenstrain=thermal_eigenstrain(expansion,Eigen::Matrix3d::Identity(),350.,fan.temperatures_k[0]);
    const auto fan_stress=solve(nodes,{element},fixed).stress_pa[0];
    CHECK_THAT(fan_stress[0],WithinAbs(.25,1e-9));
    CHECK_THAT(fan_stress[1],WithinAbs(.5,1e-9));
    element.eigenstrain=thermal_eigenstrain(expansion,Eigen::Matrix3d::Identity(),350.,both.temperatures_k[0]);
    const auto combined_stress=solve(nodes,{element},fixed).stress_pa[0];
    CHECK(combined_stress[0]>fan_stress[0]);
    CHECK_THAT(combined_stress[0],WithinAbs(1./3.,1e-9));
}
