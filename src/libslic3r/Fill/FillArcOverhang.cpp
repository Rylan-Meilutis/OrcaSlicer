#include "FillArcOverhang.hpp"

#include "../ClipperUtils.hpp"
#include "../Geometry.hpp"

#include <algorithm>
#include <cmath>

namespace Slic3r {

void FillArcOverhang::_fill_surface_single(const FillParams              &params,
                                           unsigned int,
                                           const std::pair<float, Point> &direction,
                                           ExPolygon                      expolygon,
                                           Polylines                     &polylines_out)
{
    if (expolygon.empty() || params.density <= 0.f)
        return;

    const BoundingBox bbox = get_extents(expolygon);
    const coord_t spacing = std::max<coord_t>(1, scale_(this->spacing) / params.density);
    const Vec2d axis(std::cos(direction.first), std::sin(direction.first));
    const Point midpoint = bbox.center();
    const double ray_length = bbox.size().cast<double>().norm() + 4. * spacing;

    Polyline ray;
    ray.points.emplace_back(midpoint - (axis * ray_length).cast<coord_t>());
    ray.points.emplace_back(midpoint + (axis * ray_length).cast<coord_t>());
    Polylines clipped_ray = intersection_pl(ray, expolygon);
    if (clipped_ray.empty())
        return;

    Point anchor = clipped_ray.front().first_point();
    double anchor_projection = anchor.cast<double>().dot(axis);
    for (const Polyline &segment : clipped_ray) {
        for (const Point &point : {segment.first_point(), segment.last_point()}) {
            const double projection = point.cast<double>().dot(axis);
            if (projection < anchor_projection) {
                anchor = point;
                anchor_projection = projection;
            }
        }
    }

    // Keep the center just inside the bridge surface so the smallest arc is
    // anchored by a useful length instead of degenerating to a point.
    const Point center = anchor + (axis * (0.5 * spacing)).cast<coord_t>();
    const double max_radius = std::max((bbox.min - center).cast<double>().norm(),
                              std::max((bbox.max - center).cast<double>().norm(),
                              std::max((Point(bbox.min.x(), bbox.max.y()) - center).cast<double>().norm(),
                                       (Point(bbox.max.x(), bbox.min.y()) - center).cast<double>().norm()))) + spacing;
    const double resolution = std::max<double>(scale_(params.resolution), SCALED_EPSILON);

    Polylines arcs;
    for (double radius = 0.5 * spacing; radius <= max_radius; radius += spacing) {
        const double angle_step = std::clamp(2. * std::acos(std::max(-1., 1. - resolution / radius)),
                                             Geometry::deg2rad(2.), Geometry::deg2rad(15.));
        Polyline circle;
        const size_t segments = std::max<size_t>(24, size_t(std::ceil(2. * M_PI / angle_step)));
        circle.points.reserve(segments + 1);
        for (size_t i = 0; i <= segments; ++i) {
            const double angle = 2. * M_PI * double(i) / double(segments);
            circle.points.emplace_back(center.x() + coord_t(std::lround(radius * std::cos(angle))),
                                       center.y() + coord_t(std::lround(radius * std::sin(angle))));
        }
        append(arcs, intersection_pl(circle, expolygon));
    }

    std::sort(arcs.begin(), arcs.end(), [&center](const Polyline &lhs, const Polyline &rhs) {
        return (lhs.first_point() - center).squaredNorm() < (rhs.first_point() - center).squaredNorm();
    });
    append(polylines_out, std::move(arcs));
}

} // namespace Slic3r
