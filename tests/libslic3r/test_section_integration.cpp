#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/SectionIntegration.hpp"
#include "libslic3r/ThermalDeposition.hpp"
#include "libslic3r/ThermalHistory.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Polyline.hpp"
#include "libslic3r/libslic3r.h"
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::SectionIntegration;
using Catch::Matchers::WithinAbs;

namespace {
ExPolygon rectangle(coord_t x0, coord_t y0, coord_t x1, coord_t y1)
{
    return ExPolygon(Polygon{{x0,y0},{x1,y0},{x1,y1},{x0,y1}});
}
const std::array<coord_t,3> origin{0,0,0}, size{1000000000,1000000000,1000000000};
}

TEST_CASE("Separated vertical moments match independent Gaussian integration for every monomial", "[SectionIntegration]")
{
    const Eigen::Vector3d metric(.003,.002,.001);
    const std::array<double,3> nodes{.5-.5*std::sqrt(3./5.),.5,.5+.5*std::sqrt(3./5.)};
    const std::array<double,3> weights{5./18.,4./9.,5./18.};
    const std::array<std::array<Eigen::Vector2d,2>,6> edges{{
        {Eigen::Vector2d(.2,.8),Eigen::Vector2d(.9,.1)},
        {Eigen::Vector2d(0.,.3),Eigen::Vector2d(1.,.3)},
        {Eigen::Vector2d(.3,0.),Eigen::Vector2d(.3,1.)},
        {Eigen::Vector2d(.3,.3),Eigen::Vector2d(.3,.3)},
        {Eigen::Vector2d(0.,1e-12),Eigen::Vector2d(1e-12,0.)},
        {Eigen::Vector2d(1.-1e-12,0.),Eigen::Vector2d(1.,1e-12)}}};
    for (const auto& edge : edges) for (const bool reverse : {false,true}) {
        const Eigen::Vector2d a = edge[reverse ? 1 : 0];
        const Eigen::Vector2d b = edge[reverse ? 0 : 1];
        const auto xy = vertical_segment_xy_moments(a,b,metric.head<2>());
        for (const auto interval : {std::array<double,2>{0.,1.},std::array<double,2>{.2,.7},
                                   std::array<double,2>{.999999999,1.}}) {
            const auto moments = extrude_line_moments(xy,metric.z(),interval[0],interval[1]);
            const double area = (b-a).cwiseProduct(metric.head<2>()).norm()*metric.z()*(interval[1]-interval[0]);
            for (unsigned p = 0; p < 3; ++p)
                for (unsigned q = 0; q < 3; ++q)
                    for (unsigned r = 0; r < 3; ++r) {
                        double expected = 0.;
                        for (unsigned i = 0; i < 3; ++i)
                            for (const double s : {-1.,1.}) {
                                const Eigen::Vector2d point = a+nodes[i]*(b-a);
                                const double z = interval[0]+(interval[1]-interval[0])*(.5+s/(2.*std::sqrt(3.)));
                                expected += area*weights[i]*.5*std::pow(point.x(),int(p))*
                                    std::pow(point.y(),int(q))*std::pow(z,int(r));
                            }
                        CHECK_THAT(moments[p+3*q+9*r],WithinAbs(expected,area*1e-13));
                        if (expected > 0.)
                            CHECK_THAT(moments[p+3*q+9*r],Catch::Matchers::WithinRel(expected,2e-12));
                    }
        }
    }
}

TEST_CASE("Polygon surface moments match analytic triangle integrals in either winding", "[SectionIntegration]")
{
    constexpr std::array<double,7> factorial{1.,1.,2.,6.,24.,120.,720.};
    for (const bool reverse : {false,true}) {
        ExPolygon triangle(Polygon{{0,0},{size[0],0},{0,size[1]}});
        if (reverse) triangle.contour.reverse();
        const auto surface=horizontal_surface_moments({triangle},origin,size,size[2]/2);
        const auto volume=volume_moments({triangle},origin,size,0,size[2]);
        for (unsigned p=0;p<3;++p)
            for (unsigned q=0;q<3;++q)
                for (unsigned r=0;r<3;++r) {
                    const double xy=factorial[p]*factorial[q]/factorial[p+q+2];
                    CHECK_THAT(surface[p+3*q+9*r],WithinAbs(xy*std::pow(.5,int(r)),1e-14));
                    CHECK_THAT(volume[p+3*q+9*r],WithinAbs(xy/(r+1),1e-14));
                }
    }
}

TEST_CASE("Tiny translated surface triangles retain physical convection operators", "[SectionIntegration]")
{
    for (const coord_t offset : {coord_t(0),coord_t(500000000),coord_t(999999990)}) {
        ExPolygon triangle(Polygon{{offset,offset},{offset+3,offset},{offset,offset+3}});
        const auto moments=horizontal_surface_moments({triangle},origin,size,size[2]/2);
        CHECK_THAT(moments[0],Catch::Matchers::WithinRel(4.5e-18,1e-10));
        const auto robin=thermal_boundary_element(moments,20.,300.);
        ThermalHistory::Element element;
        for (unsigned i=0;i<8;++i) element.nodes[i]=i;
        element.boundary=true;
        element.conductance=robin.conductance;
        element.load=robin.load;
        CHECK(ThermalHistory::valid_element(element,8));
    }
}

TEST_CASE("Diagonal lumped thermal quadrature preserves affine conduction and temperature bounds", "[SectionIntegration]")
{
    for (const auto& footprint : {rectangle(100000000,200000000,800000000,900000000),
                                ExPolygon(Polygon{{100000000,200000000},{800000000,300000000},{400000000,900000000}})}) {
        const auto moments=volume_moments({footprint},origin,size,0,size[2]);
        const Eigen::Matrix3d conductivity=Eigen::Vector3d(1.,2.,3.).asDiagonal();
        const auto reference=thermal_element(moments,Eigen::Vector3d::Ones(),1.,conductivity);
        const auto lumped=thermal_element(moments,Eigen::Vector3d::Ones(),1.,conductivity,ThermalQuadrature::LumpedDiagonal);
        CHECK_THAT((lumped.capacity-reference.capacity).norm(),WithinAbs(0.,1e-14));
        for (unsigned axis=0;axis<3;++axis) {
            Eigen::Matrix<double,8,1> coordinate;
            for (unsigned i=0;i<8;++i) coordinate[i]=(i>>axis)&1u;
            CHECK_THAT(((lumped.conduction-reference.conduction)*coordinate).norm(),WithinAbs(0.,1e-12));
        }
        std::vector<ThermalHistory::Cell> cells;
        ThermalHistory::Element internal;
        internal.conductance=lumped.conduction;
        for (unsigned i=0;i<8;++i) {
            internal.nodes[i]=i;
            cells.push_back({lumped.capacity.row(i).sum(),1.,400.,0.});
            CHECK_THAT(lumped.conduction.row(i).sum(),WithinAbs(0.,1e-12));
            for (unsigned j=0;j<8;++j) if (i!=j) CHECK(lumped.conduction(i,j)<=0.);
        }
        std::vector<ThermalHistory::Element> elements{internal};
        for (const bool bed : {false,true}) {
            const auto face=horizontal_surface_moments({footprint},origin,size,bed ? 0 : size[2]);
            const auto robin=thermal_boundary_element(face,10.,bed ? 400. : 300.,ThermalQuadrature::LumpedDiagonal);
            CHECK_THAT(robin.conductance.sum(),Catch::Matchers::WithinRel(10.*face[0],1e-12));
            CHECK_THAT(robin.load.sum(),Catch::Matchers::WithinRel(10.*face[0]*(bed ? 400. : 300.),1e-12));
            for (unsigned i=0;i<8;++i) if (bool(i&4u)==bed) {
                CHECK_THAT(robin.conductance(i,i),WithinAbs(0.,0.));
                CHECK_THAT(robin.load[i],WithinAbs(0.,0.));
            }
            ThermalHistory::Element boundary;
            boundary.nodes=internal.nodes;
            boundary.conductance=robin.conductance;
            boundary.load=robin.load;
            boundary.boundary=true;
            elements.push_back(boundary);
        }
        const auto cooled=ThermalHistory::step(cells,{}, {},100.,{}, {},elements);
        REQUIRE(cooled.state==ThermalHistory::State::Converged);
        for (const double temperature : cooled.temperatures_k) {
            CHECK(temperature>=300.-1e-8);
            CHECK(temperature<=400.+1e-8);
        }
        CHECK_THAT(cooled.energy_balance_error_j,WithinAbs(0.,1e-8));
        Eigen::Matrix3d rotated=conductivity;
        rotated(0,1)=rotated(1,0)=.1;
        CHECK_THROWS_AS(thermal_element(moments,Eigen::Vector3d::Ones(),1.,rotated,ThermalQuadrature::LumpedDiagonal),
                        std::invalid_argument);
        CHECK_THROWS_AS(thermal_element(moments,Eigen::Vector3d::Ones(),1.,conductivity,static_cast<ThermalQuadrature>(9)),
                        std::invalid_argument);
    }
}

TEST_CASE("Diagonal lumped conduction converges for a three dimensional insulated heat mode", "[SectionIntegration]")
{
    // On a unit insulated cube, cos(pi*x)cos(pi*y)cos(pi*z) decays
    // with rate (kx+ky+kz)*pi^2 for unit heat capacity. This checks
    // curved, mixed-axis fields, not only the affine patch condition.
    const double pi=std::acos(-1.);
    const double rate=6.*pi*pi;
    double previous_error=1.;
    for (const unsigned divisions : {4u,8u,16u}) {
        const unsigned side=divisions+1;
        const std::size_t count=std::size_t(side)*side*side;
        const auto index=[side](unsigned x,unsigned y,unsigned z) {
            return std::size_t(x)+std::size_t(side)*(y+std::size_t(side)*z);
        };
        const coord_t native_step=size[0]/divisions;
        const std::array<coord_t,3> cell_size{native_step,native_step,native_step};
        const auto moments=volume_moments({rectangle(0,0,native_step,native_step)},origin,cell_size,0,native_step);
        const auto element=thermal_element(moments,Eigen::Vector3d::Constant(1./divisions),1.,
            Eigen::Vector3d(1.,2.,3.).asDiagonal(),ThermalQuadrature::LumpedDiagonal);
        std::vector<double> mode(count),mass(count,0.),action(count,0.);
        for (unsigned z=0;z<side;++z) for (unsigned y=0;y<side;++y) for (unsigned x=0;x<side;++x)
            mode[index(x,y,z)]=std::cos(pi*x/divisions)*std::cos(pi*y/divisions)*std::cos(pi*z/divisions);
        for (unsigned z=0;z<divisions;++z) for (unsigned y=0;y<divisions;++y) for (unsigned x=0;x<divisions;++x) {
            std::array<std::size_t,8> nodes;
            Eigen::Matrix<double,8,1> local;
            for (unsigned i=0;i<8;++i) {
                nodes[i]=index(x+(i&1u),y+((i>>1)&1u),z+((i>>2)&1u));
                local[i]=mode[nodes[i]];
            }
            const Eigen::Matrix<double,8,1> flux=element.conduction*local;
            for (unsigned i=0;i<8;++i) {
                mass[nodes[i]]+=element.capacity.row(i).sum();
                action[nodes[i]]+=flux[i];
            }
        }
        double squared_error=0.,squared_reference=0.;
        for (std::size_t i=0;i<count;++i) {
            const double expected=rate*mass[i]*mode[i];
            squared_error+=std::pow(action[i]-expected,2);
            squared_reference+=expected*expected;
        }
        const double error=std::sqrt(squared_error/squared_reference);
        CAPTURE(divisions,error);
        CHECK(error<.06);
        if (divisions>4) CHECK(error<.27*previous_error);
        previous_error=error;
    }
}

TEST_CASE("Moment rebasing preserves independently integrated physical moments", "[SectionIntegration]")
{
    const ExPolygons footprint{rectangle(310000000,220000000,610000000,720000000)};
    const auto original=volume_moments(footprint,origin,size,150000000,850000000);
    const std::array<coord_t,3> shifted{300000000,200000000,100000000}, smaller{400000000,600000000,800000000};
    Eigen::Vector3d offset,scale;
    for (unsigned axis=0;axis<3;++axis) {
        offset[axis]=-double(shifted[axis])/double(smaller[axis]);
        scale[axis]=double(size[axis])/double(smaller[axis]);
    }
    const auto transformed=rebase_moments(original,offset,scale);
    const auto direct=volume_moments(footprint,shifted,smaller,150000000,850000000);
    for (unsigned i=0;i<27;++i) CHECK_THAT(transformed[i],WithinAbs(direct[i],1e-13));
}

TEST_CASE("Concave surface moments are independent of the fan origin and winding", "[SectionIntegration]")
{
    const Polygon shape{{0,0},{1000000000,0},{1000000000,250000000},
                        {250000000,250000000},{250000000,1000000000},{0,1000000000}};
    for (std::size_t start=0;start<shape.points.size();++start)
        for (const bool reverse : {false,true}) {
            auto ring=shape;
            std::rotate(ring.points.begin(),ring.points.begin()+start,ring.points.end());
            if (reverse) ring.reverse();
            const auto moments=horizontal_surface_moments({ExPolygon(ring)},origin,size,size[2]/2);
            for (unsigned p=0;p<3;++p)
                for (unsigned q=0;q<3;++q)
                    for (unsigned r=0;r<3;++r) {
                        // Disjoint bottom strip and upper left strip.
                        const double xy=(std::pow(.25,int(q+1))+
                            std::pow(.25,int(p+1))*(1.-std::pow(.25,int(q+1))))/double((p+1)*(q+1));
                        CHECK_THAT(moments[p+3*q+9*r],WithinAbs(xy*std::pow(.5,int(r)),1e-14));
                    }
        }
}

TEST_CASE("Polygon moments preserve holes under translated anisotropic coordinates", "[SectionIntegration]")
{
    const std::array<coord_t,3> start{17000000,-31000000,7000000}, extent{4000000,8000000,2000000};
    for (const bool reverse_outer : {false,true}) for (const bool reverse_hole : {false,true}) {
        auto polygon=rectangle(start[0],start[1],start[0]+extent[0],start[1]+extent[1]);
        polygon.holes.push_back(rectangle(start[0]+extent[0]/4,start[1]+extent[1]/4,
            start[0]+3*extent[0]/4,start[1]+3*extent[1]/4).contour);
        if (reverse_outer) polygon.contour.reverse();
        if (reverse_hole) polygon.holes[0].reverse();
        const auto actual=horizontal_surface_moments({polygon},start,extent,start[2]+extent[2]/2);
        for (unsigned p=0;p<3;++p)
            for (unsigned q=0;q<3;++q)
                for (unsigned r=0;r<3;++r) {
                    const auto hole_integral=[](unsigned power) {
                        return (std::pow(.75,int(power+1))-std::pow(.25,int(power+1)))/(power+1);
                    };
                    const double xy=1./((p+1)*(q+1))-hole_integral(p)*hole_integral(q);
                    CHECK_THAT(actual[p+3*q+9*r],Catch::Matchers::WithinRel(
                        xy*double(extent[0])*double(extent[1])*1e-18*std::pow(.5,int(r)),1e-12));
                }
    }
}

TEST_CASE("Section capacity integrates the consistent trilinear mass matrix", "[SectionIntegration]")
{
    const auto moments = volume_moments({rectangle(0,0,size[0],size[1])}, origin, size, 0, size[2]);
    const auto element = thermal_element(moments, Eigen::Vector3d::Ones(), 216., Eigen::Matrix3d::Identity());
    for (unsigned a = 0; a < 8; ++a)
        for (unsigned b = 0; b < 8; ++b) {
            double expected = 1.;
            for (unsigned axis = 0; axis < 3; ++axis)
                if ((a & (1u << axis)) == (b & (1u << axis))) expected *= 2.;
            CHECK_THAT(element.capacity(a,b), WithinAbs(expected,1e-12));
        }
    CHECK_THAT(element.capacity.sum(), WithinAbs(216.,1e-12));
    CHECK_THAT((element.conduction * Eigen::Matrix<double,8,1>::Ones()).norm(), WithinAbs(0.,1e-12));
    CHECK(element.capacity.llt().info() == Eigen::Success);
}

TEST_CASE("Thermal area moments preserve translated tiny holes", "[SectionIntegration]")
{
    const auto integrate=[](coord_t offset) {
        auto polygon=rectangle(offset,offset,offset+100,offset+100);
        polygon.holes.push_back(Polygon{{offset+2,offset+2},{offset+2,offset+5},{offset+5,offset+2}});
        return horizontal_surface_moments({polygon},{offset,offset,0},{100,100,100},50);
    };
    const auto reference=integrate(0), translated=integrate(1000000000);
    CHECK_THAT(reference[0],Catch::Matchers::WithinRel((10000.-4.5)*1e-18,1e-12));
    for (std::size_t i=0;i<reference.size();++i)
        CHECK_THAT(translated[i],Catch::Matchers::WithinRel(reference[i],1e-12));
}

TEST_CASE("Thin upper volume bands retain all positive polynomial moments", "[SectionIntegration]")
{
    const auto actual=volume_moments({rectangle(0,0,size[0],size[1])},origin,size,size[2]-1,size[2]);
    // Integrate (1-u)^r over the one-native-unit interval 0 <= u <= t.
    const double t=1./double(size[2]);
    const std::array<double,3> z{t,t-t*t/2.,t-t*t+t*t*t/3.};
    for (unsigned p=0;p<3;++p)
        for (unsigned q=0;q<3;++q)
            for (unsigned r=0;r<3;++r) {
                const double expected=z[r]/((p+1)*(q+1));
                CHECK(actual[p+3*q+9*r]>0.);
                CHECK_THAT(actual[p+3*q+9*r],Catch::Matchers::WithinRel(expected,1e-12));
            }
}

TEST_CASE("Section conduction preserves holes and anisotropic linear temperature gradients", "[SectionIntegration]")
{
    auto polygon = rectangle(0,0,size[0],size[1]);
    polygon.holes.push_back(rectangle(250000000,250000000,750000000,750000000).contour);
    const auto moments = volume_moments({polygon}, origin, size, 0, size[2]);
    Eigen::Matrix3d conductivity;
    conductivity << 2.,.3,.1, .3,3.,.2, .1,.2,4.;
    const auto element = thermal_element(moments,Eigen::Vector3d::Ones(),2.,conductivity);
    const Eigen::Vector3d gradient(3.,-2.,4.);
    Eigen::Matrix<double,8,1> temperature = Eigen::Matrix<double,8,1>::Zero();
    for (unsigned node = 0; node < 8; ++node)
        for (unsigned axis = 0; axis < 3; ++axis)
            if (node & (1u << axis)) temperature[node] += gradient[axis];
    CHECK_THAT(moments[0],WithinAbs(.75,1e-12));
    CHECK_THAT(element.capacity.sum(),WithinAbs(1.5,1e-12));
    CHECK_THAT(temperature.dot(element.conduction * temperature),
        WithinAbs(.75 * gradient.dot(conductivity * gradient),1e-10));
    CHECK_THAT((element.conduction-element.conduction.transpose()).norm(),WithinAbs(0.,1e-12));
}

TEST_CASE("Section thermal integration is additive across deposition subintervals", "[SectionIntegration]")
{
    const ExPolygons footprint{rectangle(0,0,500000000,size[1])};
    const auto whole = thermal_element(volume_moments(footprint,origin,size,0,size[2]),
        Eigen::Vector3d::Ones(),2.,Eigen::Matrix3d::Identity());
    const auto lower = thermal_element(volume_moments(footprint,origin,size,0,250000000),
        Eigen::Vector3d::Ones(),2.,Eigen::Matrix3d::Identity());
    const auto upper = thermal_element(volume_moments(footprint,origin,size,250000000,size[2]),
        Eigen::Vector3d::Ones(),2.,Eigen::Matrix3d::Identity());
    CHECK_THAT((whole.capacity-lower.capacity-upper.capacity).norm(),WithinAbs(0.,1e-12));
    CHECK_THAT((whole.conduction-lower.conduction-upper.conduction).norm(),WithinAbs(0.,1e-12));
    CHECK_THAT(whole.capacity.sum(),WithinAbs(1.,1e-12));
}

TEST_CASE("Section thermal integration rejects invalid constitutive properties", "[SectionIntegration]")
{
    const auto moments = volume_moments({rectangle(0,0,size[0],size[1])},origin,size,0,size[2]);
    CHECK_THROWS_AS(thermal_element(moments,Eigen::Vector3d::Zero(),1.,Eigen::Matrix3d::Identity()),std::invalid_argument);
    CHECK_THROWS_AS(thermal_element(moments,Eigen::Vector3d::Ones(),-1.,Eigen::Matrix3d::Identity()),std::invalid_argument);
    CHECK_THROWS_AS(thermal_element(moments,Eigen::Vector3d::Ones(),1.,-Eigen::Matrix3d::Identity()),std::invalid_argument);
    auto invalid = moments;
    invalid[2] = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(thermal_element(invalid,Eigen::Vector3d::Ones(),1.,Eigen::Matrix3d::Identity()),std::invalid_argument);
}

TEST_CASE("Section surface heat exchange excludes holes and inactive opposite face nodes", "[SectionIntegration]")
{
    auto polygon = rectangle(0,0,size[0],size[1]);
    polygon.holes.push_back(rectangle(250000000,250000000,750000000,750000000).contour);
    const auto face = horizontal_surface_moments({polygon},origin,size,size[2]);
    const auto boundary = thermal_boundary_element(face,10.,300.);
    const Eigen::Matrix<double,8,1> temperature = Eigen::Matrix<double,8,1>::Constant(350.);
    CHECK_THAT((boundary.load-boundary.conductance*temperature).sum(),WithinAbs(-375.,1e-9));
    CHECK_THAT(boundary.load.head<4>().norm(),WithinAbs(0.,1e-12));
    CHECK_THAT(boundary.conductance.topRows<4>().norm(),WithinAbs(0.,1e-12));
    CHECK_THAT(boundary.conductance.sum(),WithinAbs(7.5,1e-12));
    CHECK_THROWS_AS(thermal_boundary_element(face,-1.,300.),std::invalid_argument);
}

TEST_CASE("Section thermal elements recover the analytical bed to air conduction solution", "[SectionIntegration]")
{
    const std::array<coord_t,3> cell_size{10000000,10000000,10000000};
    const ExPolygons footprint{rectangle(0,0,cell_size[0],cell_size[1])};
    const auto element = thermal_element(volume_moments(footprint,origin,cell_size,0,cell_size[2]),
        Eigen::Vector3d::Constant(.01),1e6,.2*Eigen::Matrix3d::Identity());
    const auto bed = thermal_boundary_element(horizontal_surface_moments(footprint,origin,cell_size,0),1000.,383.);
    const auto air = thermal_boundary_element(horizontal_surface_moments(footprint,origin,cell_size,cell_size[2]),10.,323.);
    const Eigen::Matrix<double,8,8> operator_matrix = element.conduction+bed.conductance+air.conductance;
    const Eigen::Matrix<double,8,1> temperature = operator_matrix.ldlt().solve(bed.load+air.load);
    const double flux = (383.-323.)/(1./1000.+.01/.2+1./10.);
    for (unsigned node = 0; node < 8; ++node)
        CHECK_THAT(temperature[node],WithinAbs(node & 4 ? 323.+flux/10. : 383.-flux/1000.,1e-9));
    const double bed_heat = (bed.load-bed.conductance*temperature).sum();
    const double air_heat = (air.load-air.conductance*temperature).sum();
    CHECK_THAT(bed_heat,WithinAbs(flux*.0001,1e-10));
    CHECK_THAT(bed_heat+air_heat,WithinAbs(0.,1e-10));
    // Capacity belongs only to polymer, not the bounding air reservoir.
    CHECK_THAT(element.capacity.sum(),WithinAbs(1.,1e-12));
}

TEST_CASE("Vertical section heat exchange integrates oblique exposed edges without bounding box area", "[SectionIntegration]")
{
    const Polyline edge{{0,0},{size[0],size[1]}};
    const auto face = vertical_surface_moments(edge,origin,size,250000000,750000000);
    const double area = .5 * std::sqrt(2.);
    CHECK_THAT(face[0],WithinAbs(area,1e-12));
    const auto boundary = thermal_boundary_element(face,10.,300.);
    CHECK_THAT(boundary.conductance.sum(),WithinAbs(10.*area,1e-12));
    CHECK_THAT(boundary.load.sum(),WithinAbs(3000.*area,1e-9));
    // x=y=t along the oblique face: integral x*y dA = area/3.
    CHECK_THAT(face[4],WithinAbs(area/3.,1e-12));
    const Polyline reversed{{size[0],size[1]},{0,0}};
    const auto reverse_face = vertical_surface_moments(reversed,origin,size,250000000,750000000);
    for (std::size_t i = 0; i < face.size(); ++i) CHECK_THAT(reverse_face[i],WithinAbs(face[i],1e-12));
}

TEST_CASE("Section transient conduction conserves heat without a boundary reservoir", "[SectionIntegration]")
{
    const auto moments = volume_moments({rectangle(0,0,size[0],size[1])},origin,size,0,size[2]);
    const auto element = thermal_element(moments,Eigen::Vector3d::Ones(),2.,Eigen::Matrix3d::Identity());
    Eigen::Matrix<double,8,1> temperature;
    for (unsigned node = 0; node < 8; ++node) temperature[node] = node & 4 ? 400. : 300.;
    const double initial_heat = (element.capacity*temperature).sum();
    const double dt = .1;
    const Eigen::Matrix<double,8,8> matrix = element.capacity/dt+element.conduction;
    const auto solve = matrix.ldlt();
    double variance = (temperature.array()-350.).square().sum();
    for (unsigned step = 0; step < 20; ++step) {
        temperature = solve.solve(element.capacity*temperature/dt).eval();
        CHECK_THAT((element.capacity*temperature).sum(),WithinAbs(initial_heat,1e-9));
        const double next_variance = (temperature.array()-350.).square().sum();
        CHECK(next_variance < variance);
        variance = next_variance;
    }
}

TEST_CASE("Global section history activates shared nodes without erasing deposited heat", "[SectionIntegration]")
{
    namespace TH = Slic3r::ThermalHistory;
    const auto element = thermal_element(volume_moments({rectangle(0,0,size[0],size[1])},origin,size,0,size[2]),
        Eigen::Vector3d::Ones(),2.,Eigen::Matrix3d::Identity());
    TH::Element lower, upper;
    lower.conductance = upper.conductance = element.conduction;
    TH::DepositionEvent first, second;
    second.time_s = 1.;
    for (unsigned node = 0; node < 8; ++node) {
        lower.nodes[node] = node;
        upper.nodes[node] = node+4;
        const double volume = element.capacity.row(node).sum()/2.;
        first.deposited_cells.push_back({volume,2.,300.,0.});
        if (node < 4) second.added_mass.push_back({node+4,volume,2.,500.});
        else second.deposited_cells.push_back({volume,2.,500.,0.});
    }
    first.added_elements.push_back(lower);
    second.added_elements.push_back(upper);
    TH::HistoryObservations observations;
    observations.times_s = {1.};
    observations.observe = [](double, const std::vector<TH::Cell>& cells) {
        REQUIRE(cells.size() == 12);
        for (unsigned node = 0; node < 12; ++node)
            CHECK_THAT(cells[node].temperature_k,WithinAbs(node < 4 ? 300. : (node < 8 ? 400. : 500.),1e-9));
    };
    TH::HistoryLimits limits;
    limits.max_step_temperature_error_k = .01;
    const auto result = TH::evolve({first,second},20.,limits,{},observations);
    REQUIRE(result.state == TH::State::Converged);
    CHECK_THAT(result.deposited_energy_j,WithinAbs(1600.,1e-9));
    CHECK_THAT(result.final_energy_j,WithinAbs(1600.,1e-7));
    CHECK_THAT(result.energy_balance_error_j,WithinAbs(0.,1e-7));
    CHECK(result.temperatures_k.front() > 300.);
    CHECK(result.temperatures_k.back() < 500.);
}

TEST_CASE("Global section history recovers conduction across shared element faces", "[SectionIntegration]")
{
    namespace TH = Slic3r::ThermalHistory;
    const ExPolygons footprint{rectangle(0,0,size[0],size[1])};
    const auto volume = thermal_element(volume_moments(footprint,origin,size,0,size[2]),
        Eigen::Vector3d::Ones(),1.,Eigen::Matrix3d::Identity());
    TH::DepositionEvent event;
    for (unsigned node = 0; node < 12; ++node)
        event.deposited_cells.push_back({node >= 4 && node < 8 ? .25 : .125,1.,350.,0.});
    for (unsigned cell = 0; cell < 2; ++cell) {
        TH::Element conduction;
        for (unsigned node = 0; node < 8; ++node) conduction.nodes[node] = node+4*cell;
        conduction.conductance = volume.conduction;
        event.added_elements.push_back(conduction);
        TH::Element boundary;
        boundary.nodes = conduction.nodes;
        boundary.boundary = true;
        boundary.kind = cell ? TH::BoundaryKind::Exhaust : TH::BoundaryKind::Bed;
        const auto surface = thermal_boundary_element(horizontal_surface_moments(footprint,origin,size,cell ? size[2] : 0),
            1.,cell ? 300. : 400.);
        boundary.conductance = surface.conductance;
        boundary.load = surface.load;
        event.exposed_elements.push_back(boundary);
    }
    const auto result = TH::evolve({event},100.);
    REQUIRE(result.state == TH::State::Converged);
    REQUIRE(result.temperatures_k.size() == 12);
    for (unsigned node = 0; node < 12; ++node)
        CHECK_THAT(result.temperatures_k[node],WithinAbs(node < 4 ? 375. : (node < 8 ? 350. : 325.),1e-6));
    CHECK_THAT(result.energy_balance_error_j,WithinAbs(0.,1e-6));
    CHECK(result.boundary_energy_by_kind_j[static_cast<std::size_t>(TH::BoundaryKind::Bed)] > 0.);
    CHECK(result.boundary_energy_by_kind_j[static_cast<std::size_t>(TH::BoundaryKind::Exhaust)] < 0.);
}

TEST_CASE("Section deposition removes covered convection operators from subsequent heat steps", "[SectionIntegration]")
{
    namespace TH = Slic3r::ThermalHistory;
    const ExPolygons footprint{rectangle(0,0,size[0],size[1])};
    const auto solid = thermal_element(volume_moments(footprint,origin,size,0,size[2]),
        Eigen::Vector3d::Ones(),1.,Eigen::Matrix3d::Identity());
    const auto surface = thermal_boundary_element(horizontal_surface_moments(footprint,origin,size,size[2]),10.,300.);
    TH::DepositionEvent first, cover;
    cover.time_s = 1.;
    TH::Element conduction, boundary;
    for (unsigned node = 0; node < 8; ++node) {
        conduction.nodes[node] = boundary.nodes[node] = node;
        first.deposited_cells.push_back({.125,1.,400.,0.});
    }
    conduction.conductance = solid.conduction;
    boundary.boundary = true;
    boundary.kind = TH::BoundaryKind::PartFan;
    boundary.conductance = surface.conductance;
    boundary.load = surface.load;
    first.added_elements.push_back(conduction);
    first.exposed_elements.push_back(boundary);
    // Explicitly insulated after this exposure event. No old boundary survives.
    double heat_at_cover = 0.;
    TH::HistoryObservations observations;
    observations.times_s = {1.};
    observations.observe = [&](double,const std::vector<TH::Cell>& cells) {
        for (const auto& cell : cells)
            heat_at_cover += cell.volume_m3*cell.volumetric_heat_capacity_j_m3_k*cell.temperature_k;
    };
    const auto result = TH::evolve({first,cover},3.,{}, {},observations);
    REQUIRE(result.state == TH::State::Converged);
    CHECK(heat_at_cover < 400.);
    CHECK_THAT(result.final_energy_j,WithinAbs(heat_at_cover,1e-8));
    CHECK_THAT(result.boundary_energy_by_kind_j[static_cast<std::size_t>(TH::BoundaryKind::PartFan)],
        WithinAbs(heat_at_cover-400.,1e-8));
}
