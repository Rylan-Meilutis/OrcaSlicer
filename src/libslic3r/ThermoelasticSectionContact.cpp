#include "ThermoelasticSectionContact.hpp"
#include "ExPolygon.hpp"
#include "Point.hpp"
#include "Polygon.hpp"
#include "libslic3r.h"
#include <Eigen/Core>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>

namespace Slic3r::Thermoelastic {
SectionTraceExtrema section_trace_extrema(const ExPolygons& footprint,
    const std::array<coord_t,2>& origin, const std::array<coord_t,2>& size,
    const std::array<double,4>& values, const std::function<void()>& checkpoint)
{
    for (unsigned axis = 0; axis < 2; ++axis)
        if (size[axis] <= 0 || size[axis] > 10000000000LL || std::abs(double(origin[axis])) > 1e10)
            throw std::invalid_argument("Invalid section trace cell.");
    for (double value : values)
        if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite section trace field.");
    const double a = values[0], b = values[1]-a, c = values[2]-a, d = values[3]-values[1]-values[2]+a;
    if (!std::isfinite(b) || !std::isfinite(c) || !std::isfinite(d))
        throw std::invalid_argument("Unresolved section trace polynomial.");
    SectionTraceExtrema result{};
    bool initialized = false;
    const auto consider = [&](const Eigen::Vector2d& p) {
        const double value = a+b*p.x()+c*p.y()+d*p.x()*p.y();
        if (!std::isfinite(value)) throw std::runtime_error("Nonfinite section trace value.");
        if (!initialized || value < result.minimum) { result.minimum = value; result.minimum_local = p; }
        if (!initialized || value > result.maximum) { result.maximum = value; result.maximum_local = p; }
        initialized = true;
    };
    const auto local = [&](const Point& p) {
        Eigen::Vector2d q;
        for (unsigned axis = 0; axis < 2; ++axis) {
            if (p[axis] < origin[axis] || p[axis] > origin[axis]+size[axis])
                throw std::invalid_argument("Section trace extends outside its cell.");
            q[axis] = double(p[axis]-origin[axis])/double(size[axis]);
        }
        return q;
    };
    std::size_t points = 0;
    const auto visit = [&](const Polygon& ring) {
        if (ring.points.size() < 3 || ring.points.size() > 2000000-points)
            throw std::invalid_argument("Invalid or excessive section trace boundary.");
        points += ring.points.size();
        for (std::size_t i = 0; i < ring.points.size(); ++i) {
            if (checkpoint) checkpoint();
            const Eigen::Vector2d p = local(ring.points[i]);
            const Eigen::Vector2d delta = local(ring.points[(i+1)%ring.points.size()])-p;
            consider(p);
            const double quadratic = d*delta.x()*delta.y();
            const double linear = b*delta.x()+c*delta.y()+d*(p.x()*delta.y()+p.y()*delta.x());
            if (quadratic != 0.) {
                const double t = -linear/(2.*quadratic);
                if (t > 0. && t < 1.) consider(p+t*delta);
            }
        }
    };
    // A nonconstant bilinear polynomial is harmonic, hence has no strict
    // interior extremum. Boundaries (including holes) suffice.
    for (const auto& polygon : footprint) {
        visit(polygon.contour);
        for (const auto& hole : polygon.holes) visit(hole);
    }
    if (!initialized) throw std::invalid_argument("Empty section trace footprint.");
    return result;
}
}
