#include "SectionIntegration.hpp"
#include "ExPolygon.hpp"
#include "Point.hpp"
#include "Polygon.hpp"
#include "Polyline.hpp"
#include "libslic3r.h"
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace Slic3r::SectionIntegration {
namespace {
std::array<double,9> polygon_moments(const Polygon& ring, coord_t x0, coord_t y0, coord_t dx, coord_t dy)
{
    // Signed fan triangles avoid subtracting cell-scale edge integrals for a
    // tiny translated surface. The native determinant preserves its area.
    // A Duffy map turns each triangle into a square; the polynomial including
    // its Jacobian has degree at most five in u and four in v, so 3x3 Gauss
    // integration is exact. Signed triangles also handle concave rings.
    const std::array<double,3> nodes{.5-.5*std::sqrt(3./5.),.5,.5+.5*std::sqrt(3./5.)};
    constexpr std::array<double,3> weights{5./18.,4./9.,5./18.};
    std::array<double,9> result{};
    if (ring.points.size()<3) return result;
    const auto& origin=ring.points.front();
    const double x=double(origin.x()-x0)/double(dx), y=double(origin.y()-y0)/double(dy);
    for (std::size_t i=1;i+1<ring.points.size();++i) {
        const auto a=ring.points[i]-origin, b=ring.points[i+1]-origin;
        const boost::multiprecision::int128_t determinant=
            boost::multiprecision::int128_t(a.x())*b.y()-boost::multiprecision::int128_t(a.y())*b.x();
        if (determinant==0) continue;
        const double scale=determinant.convert_to<double>()/double(dx)/double(dy);
        const double ax=double(a.x())/double(dx), ay=double(a.y())/double(dy);
        const double bx=double(b.x())/double(dx), by=double(b.y())/double(dy);
        for (unsigned u=0;u<3;++u)
            for (unsigned v=0;v<3;++v) {
                const double s=nodes[u], t=(1.-s)*nodes[v];
                const double sx=x+s*ax+t*bx, sy=y+s*ay+t*by;
                const std::array<double,3> xp{1.,sx,sx*sx}, yp{1.,sy,sy*sy};
                const double weight=scale*weights[u]*weights[v]*(1.-s);
                for (unsigned q=0;q<3;++q)
                    for (unsigned p=0;p<3;++p) result[p+3*q]+=weight*xp[p]*yp[q];
            }
    }
    return result;
}
unsigned product_index(unsigned a, unsigned b)
{
    return (a & 1) + (b & 1) + 3 * (((a >> 1) & 1) + ((b >> 1) & 1)) +
        9 * (((a >> 2) & 1) + ((b >> 2) & 1));
}
std::array<double,9> area_moments(const ExPolygons& polygons, const std::array<coord_t,3>& start,
    const std::array<coord_t,3>& size)
{
    std::array<double,9> result{};
    const auto add = [&](const Polygon& ring, double sign) {
        // Global double area can lose the sign of a native sliver/hole.
        boost::multiprecision::int128_t twice_area=0;
        for (std::size_t i=0;i<ring.points.size();++i) {
            const auto& a=ring.points[i];
            const auto& b=ring.points[(i+1)%ring.points.size()];
            twice_area+=boost::multiprecision::int128_t(a.x())*b.y()-boost::multiprecision::int128_t(a.y())*b.x();
        }
        sign *= twice_area >= 0 ? 1. : -1.;
        const auto moments=polygon_moments(ring,start[0],start[1],size[0],size[1]);
        for (unsigned i=0;i<9;++i) result[i]+=sign*moments[i];
    };
    for (const auto& polygon : polygons) {
        add(polygon.contour,1.);
        for (const auto& hole : polygon.holes) add(hole,-1.);
    }
    return result;
}
Eigen::Matrix<double, 8, 8> basis_polynomials()
{
    Eigen::Matrix<double, 8, 8> basis;
    for (unsigned node = 0; node < 8; ++node)
        for (unsigned monomial = 0; monomial < 8; ++monomial) {
            double coefficient = 1.;
            for (unsigned axis = 0; axis < 3; ++axis) {
                const bool high = node & (1u << axis), power = monomial & (1u << axis);
                coefficient *= high ? (power ? 1. : 0.) : (power ? -1. : 1.);
            }
            basis(monomial, node) = coefficient;
        }
    return basis;
}
}

Moments rebase_moments(const Moments& moments, const Eigen::Vector3d& offset, const Eigen::Vector3d& scale)
{
    if (!offset.allFinite() || !scale.allFinite() || (scale.array()<=0.).any())
        throw std::invalid_argument("Invalid moment coordinate transform.");
    Moments result=moments;
    unsigned stride=1;
    for (unsigned axis=0;axis<3;++axis,stride*=3) {
        const long double a=scale[axis], b=offset[axis];
        for (unsigned i=0;i<27;++i) if ((i/stride)%3==0) {
            const long double m0=result[i],m1=result[i+stride],m2=result[i+2*stride];
            result[i+stride]=double(a*m1+b*m0);
            result[i+2*stride]=double(a*a*m2+2*a*b*m1+b*b*m0);
        }
    }
    for (const double value : result)
        if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite transformed moment.");
    return result;
}

Moments volume_moments(const ExPolygons& polygons, const std::array<coord_t, 3>& start,
    const std::array<coord_t, 3>& size, coord_t bottom, coord_t top)
{
    Moments result{};
    const double lo = double(bottom - start[2]) / double(size[2]), hi = double(top - start[2]) / double(size[2]);
    const double volume_scale = double(size[0]) * double(size[1]) * double(size[2]) * 1e-27;
    const auto xy = area_moments(polygons,start,size);
    const double thickness=double(top-bottom)/double(size[2]);
    const std::array<double,3> z{thickness,thickness*(hi+lo)/2.,thickness*(hi*hi+hi*lo+lo*lo)/3.};
    for (unsigned p = 0; p < 3; ++p)
        for (unsigned q = 0; q < 3; ++q) {
            for (unsigned r = 0; r < 3; ++r)
                result[p + 3 * q + 9 * r] = xy[p+3*q] * volume_scale * z[r];
        }
    return result;
}

Moments horizontal_surface_moments(const ExPolygons& polygons, const std::array<coord_t, 3>& start,
    const std::array<coord_t, 3>& size, coord_t z)
{
    Moments result{};
    const double normalized_z = double(z-start[2])/double(size[2]);
    const auto xy = area_moments(polygons,start,size);
    const std::array<double,3> zp{1.,normalized_z,normalized_z*normalized_z};
    for (unsigned p = 0; p < 3; ++p)
        for (unsigned q = 0; q < 3; ++q) {
            for (unsigned r = 0; r < 3; ++r)
                result[p+3*q+9*r] = xy[p+3*q] * double(size[0]) * double(size[1]) * 1e-18 * zp[r];
        }
    return result;
}

Moments vertical_surface_moments(const Polyline& edge, const std::array<coord_t,3>& start,
    const std::array<coord_t,3>& size, coord_t bottom, coord_t top)
{
    Moments result{};
    const double lo = double(bottom-start[2])/double(size[2]), hi = double(top-start[2])/double(size[2]);
    const Eigen::Vector3d metric(double(size[0])*1e-9,double(size[1])*1e-9,double(size[2])*1e-9);
    for (std::size_t i = 1; i < edge.points.size(); ++i) {
        const auto& a = edge.points[i-1];
        const auto& b = edge.points[i];
        const auto part = vertical_segment_moments(
            {double(a.x()-start[0])/double(size[0]),double(a.y()-start[1])/double(size[1])},
            {double(b.x()-start[0])/double(size[0]),double(b.y()-start[1])/double(size[1])},metric,lo,hi);
        for (std::size_t j = 0; j < result.size(); ++j) result[j] += part[j];
    }
    return result;
}

Moments vertical_segment_moments(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
    const Eigen::Vector3d& size, double lo, double hi)
{
    return extrude_line_moments(vertical_segment_xy_moments(a,b,size.head<2>()),size.z(),lo,hi);
}

LineMoments vertical_segment_xy_moments(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
    const Eigen::Vector2d& size)
{
    LineMoments result{};
    const Eigen::Vector2d delta = b-a;
    const double scale = delta.cwiseProduct(size).norm();
    // With t in [-1/2,1/2], odd powers integrate to zero, E[t^2]=1/12
    // and E[t^4]=1/80. Expand the nine monomials once rather than running
    // nested polynomial products and divisions for every exposed edge.
    const Eigen::Vector2d midpoint = (a+b)*.5;
    const double x = midpoint.x(), y = midpoint.y();
    const double dx2 = delta.x()*delta.x(), dy2 = delta.y()*delta.y();
    const double cross = delta.x()*delta.y()/12.;
    const double xx = x*x+dx2/12., yy = y*y+dy2/12.;
    result = {1.,x,xx,y,x*y+cross,y*xx+2.*x*cross,
        yy,x*yy+2.*y*cross,
        x*x*y*y+(x*x*dy2+y*y*dx2)/12.+4.*x*y*cross+dx2*dy2/80.};
    for (auto& value : result) value *= scale;
    return result;
}

Moments extrude_line_moments(const LineMoments& xy, double height, double lo, double hi)
{
    Moments result{};
    // Factored differences avoid cancellation for very thin native Z bands.
    const double thickness = height*(hi-lo);
    const std::array<double,3> z{thickness,thickness*(hi+lo)/2.,
        thickness*(hi*hi+hi*lo+lo*lo)/3.};
    for (unsigned r = 0; r < 3; ++r)
        for (unsigned i = 0; i < 9; ++i) result[i+9*r] = xy[i]*z[r];
    return result;
}

ThermalElement thermal_element(const Moments& moments, const Eigen::Vector3d& size,
    double capacity, const Eigen::Matrix3d& conductivity, ThermalQuadrature quadrature)
{
    if (quadrature!=ThermalQuadrature::Consistent && quadrature!=ThermalQuadrature::LumpedDiagonal)
        throw std::invalid_argument("Invalid thermal quadrature.");
    if (quadrature==ThermalQuadrature::LumpedDiagonal)
        for (unsigned a=0;a<3;++a) for (unsigned b=0;b<3;++b)
            if (a!=b && conductivity(a,b)!=0.)
                throw std::invalid_argument("Lumped thermal quadrature requires grid-aligned diagonal conductivity.");
    if (!size.allFinite() || (size.array() <= 0.).any() || !std::isfinite(capacity) || capacity <= 0. ||
        !conductivity.allFinite() || !conductivity.isApprox(conductivity.transpose(), 1e-12) ||
        conductivity.llt().info() != Eigen::Success || moments[0] <= 0. ||
        !std::all_of(moments.begin(), moments.end(), [](double v) { return std::isfinite(v); }))
        throw std::invalid_argument("Invalid section thermal element properties.");
    const auto basis = basis_polynomials();
    std::array<Eigen::Matrix<double, 3, 8>, 8> gradient;
    for (auto& g : gradient) g.setZero();
    for (unsigned node = 0; node < 8; ++node)
        for (unsigned monomial = 0; monomial < 8; ++monomial) {
            for (unsigned axis = 0; axis < 3; ++axis)
                if (monomial & (1u << axis))
                    gradient[monomial ^ (1u << axis)](axis, node) += basis(monomial,node) / size[axis];
        }
    ThermalElement result;
    for (unsigned a = 0; a < 8; ++a)
        for (unsigned b = 0; b < 8; ++b) {
            const double integral = moments[product_index(a, b)];
            result.capacity.noalias() += (capacity * integral) * basis.row(a).transpose() * basis.row(b);
            result.conduction.noalias() += integral * gradient[a].transpose() * conductivity * gradient[b];
        }
    if (quadrature==ThermalQuadrature::LumpedDiagonal) {
        // Transverse gradient lumping replaces squared Q1 transverse weights
        // by their linear weights. For diagonal conductivity this gives twelve
        // nonnegative edge conductances, preserving the exact affine action.
        // Mixed modes change: this option still needs spatial refinement.
        const auto original=result.conduction;
        result.conduction.setZero();
        for (unsigned axis=0;axis<3;++axis) {
            Eigen::Matrix<double,8,1> coordinate;
            for (unsigned node=0;node<8;++node) coordinate[node]=(node>>axis)&1u;
            const auto force=(original*coordinate).eval();
            for (unsigned low=0;low<8;++low) if (!(low&(1u<<axis))) {
                const unsigned high=low|(1u<<axis);
                const double conductance=.5*(force[high]-force[low]);
                if (!std::isfinite(conductance) || conductance<0.)
                    throw std::invalid_argument("Unresolved lumped thermal edge conductance.");
                result.conduction(low,low)+=conductance;
                result.conduction(high,high)+=conductance;
                result.conduction(low,high)-=conductance;
                result.conduction(high,low)-=conductance;
            }
        }
        for (unsigned axis=0;axis<3;++axis) {
            Eigen::Matrix<double,8,1> coordinate;
            for (unsigned node=0;node<8;++node) coordinate[node]=(node>>axis)&1u;
            if (((result.conduction-original)*coordinate).norm()>1e-10*original.norm())
                throw std::invalid_argument("Lumped thermal quadrature does not preserve affine conduction.");
        }
    }
    if (!result.capacity.allFinite() || !result.conduction.allFinite())
        throw std::invalid_argument("Section thermal element overflow.");
    return result;
}

ThermalBoundaryElement thermal_boundary_element(const Moments& moments, double coefficient, double temperature,
                                                ThermalQuadrature quadrature)
{
    if (quadrature!=ThermalQuadrature::Consistent && quadrature!=ThermalQuadrature::LumpedDiagonal)
        throw std::invalid_argument("Invalid thermal boundary quadrature.");
    if (!std::isfinite(coefficient) || coefficient < 0. || !std::isfinite(temperature) || temperature <= 0. ||
        moments[0] <= 0. || !std::all_of(moments.begin(),moments.end(),[](double v) { return std::isfinite(v); }))
        throw std::invalid_argument("Invalid section thermal boundary properties.");
    const auto basis = basis_polynomials();
    ThermalBoundaryElement result;
    for (unsigned a = 0; a < 8; ++a) {
        result.load.noalias() += (coefficient * temperature * moments[product_index(a,0)]) * basis.row(a).transpose();
        for (unsigned b = 0; b < 8; ++b)
            result.conductance.noalias() += (coefficient * moments[product_index(a,b)]) *
                basis.row(a).transpose() * basis.row(b);
    }
    if (quadrature==ThermalQuadrature::LumpedDiagonal) {
        // Transform first moments directly to Q1 weights, pairing matching
        // terms before multiplication. Summing an already expanded 8x8 matrix
        // loses exact zero weights on a face at a background-cell boundary.
        Eigen::Matrix<double,8,1> weights;
        for (unsigned node=0;node<8;++node) weights[node]=moments[product_index(node,0)];
        for (unsigned axis=0;axis<3;++axis)
            for (unsigned low=0;low<8;++low) if (!(low&(1u<<axis))) weights[low]-=weights[low|(1u<<axis)];
        if (weights.minCoeff()<0.) throw std::invalid_argument("Unresolved lumped thermal boundary weight.");
        result.conductance=(coefficient*weights).asDiagonal();
        result.load=coefficient*temperature*weights;
    }
    if (!result.load.allFinite() || !result.conductance.allFinite())
        throw std::invalid_argument("Section thermal boundary overflow.");
    return result;
}
} // namespace Slic3r::SectionIntegration
