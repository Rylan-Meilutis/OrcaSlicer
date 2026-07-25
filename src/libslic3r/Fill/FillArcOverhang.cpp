#include "FillArcOverhang.hpp"

#include "../ClipperUtils.hpp"
#include "../Geometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

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
    const double resolution = std::max<double>(scale_(params.resolution), SCALED_EPSILON);

    // Emit each family from its supported center outwards. A family stops when
    // the printable half-circle reaches another boundary; continuing concentric
    // circles beyond that point creates paths whose endpoints begin in midair.
    auto generate_family = [&](const Point &family_center, const Point &start_from, double first_radius,
                               const ExPolygon &region,
                               Polylines &arcs, Polylines &outer_arcs) {
        Point previous = start_from;
        const BoundingBox region_bbox = get_extents(region);
        const double family_max_radius =
            std::max((region_bbox.min - family_center).cast<double>().norm(),
                     std::max((region_bbox.max - family_center).cast<double>().norm(),
                              std::max((Point(region_bbox.min.x(), region_bbox.max.y()) - family_center).cast<double>().norm(),
                                       (Point(region_bbox.max.x(), region_bbox.min.y()) - family_center).cast<double>().norm()))) +
            spacing;

        for (double radius = first_radius; radius <= family_max_radius; radius += spacing) {
            const double angle_step = std::clamp(2. * std::acos(std::max(-1., 1. - resolution / radius)),
                                                 Geometry::deg2rad(2.), Geometry::deg2rad(15.));
            Polyline circle;
            const size_t segments = std::max<size_t>(24, size_t(std::ceil(2. * M_PI / angle_step)));
            circle.points.reserve(segments + 1);
            for (size_t i = 0; i <= segments; ++i) {
                const double angle = 2. * M_PI * double(i) / double(segments);
                circle.points.emplace_back(family_center.x() + coord_t(std::lround(radius * std::cos(angle))),
                                           family_center.y() + coord_t(std::lround(radius * std::sin(angle))));
            }

            Polylines ring = intersection_pl(circle, region);
            ring.erase(std::remove_if(ring.begin(), ring.end(), [spacing](const Polyline &path) {
                return path.length() < 0.25 * spacing;
            }), ring.end());
            if (ring.empty())
                break;

            // Keep every next path on the endpoint side where the previous path
            // finished. This prevents a cross-over travel before each larger arc.
            Polylines ordered;
            ordered.reserve(ring.size());
            while (!ring.empty()) {
                size_t nearest_idx = 0;
                bool reverse = false;
                double nearest_distance = std::numeric_limits<double>::max();
                for (size_t idx = 0; idx < ring.size(); ++idx) {
                    const double first_distance = (ring[idx].first_point() - previous).cast<double>().squaredNorm();
                    const double last_distance = (ring[idx].last_point() - previous).cast<double>().squaredNorm();
                    if (std::min(first_distance, last_distance) < nearest_distance) {
                        nearest_distance = std::min(first_distance, last_distance);
                        nearest_idx = idx;
                        reverse = last_distance < first_distance;
                    }
                }
                Polyline path = std::move(ring[nearest_idx]);
                ring.erase(ring.begin() + nearest_idx);
                if (reverse)
                    path.reverse();
                previous = path.last_point();
                ordered.emplace_back(std::move(path));
            }

            const double clipped_length = std::accumulate(
                ordered.begin(), ordered.end(), 0.,
                [](double length, const Polyline &path) { return length + path.length(); });
            outer_arcs = ordered;
            append(arcs, std::move(ordered));

            // A supported arc family needs at least a half circle. Once another
            // boundary clips it below that amount, use this last arc as a parent
            // for recursive child families instead of extending unsupported rings.
            if (clipped_length < 0.9 * M_PI * radius)
                break;
        }
    };

    Polylines arcs;
    Polylines frontier;
    generate_family(center, center, 0.5 * spacing, expolygon, arcs, frontier);

    const bool recursive_fill = params.config != nullptr && params.config->arc_overhang_recursive_fill.value;
    if (recursive_fill && !arcs.empty()) {
        constexpr size_t max_branches = 128;
        for (size_t branch = 0; branch < max_branches && !frontier.empty(); ++branch) {
            const ExPolygons covered = union_ex(offset(arcs, float(0.52 * spacing)));
            ExPolygons remaining = diff_ex(ExPolygons{expolygon}, covered);
            // Polygonized line offsets leave narrow radial slivers between
            // otherwise adjacent arcs. They are not printable child regions.
            remaining = opening_ex(remaining, float(0.75 * spacing));
            remaining.erase(std::remove_if(remaining.begin(), remaining.end(), [spacing](const ExPolygon &part) {
                return std::abs(area(part)) < double(spacing) * double(spacing);
            }), remaining.end());
            if (remaining.empty())
                break;

            auto largest = std::max_element(remaining.begin(), remaining.end(), [](const ExPolygon &lhs, const ExPolygon &rhs) {
                return std::abs(area(lhs)) < std::abs(area(rhs));
            });
            ExPolygons branch_regions = intersection_ex(
                offset_ex(ExPolygons{*largest}, float(0.75 * spacing)), ExPolygons{expolygon});
            if (branch_regions.empty())
                break;
            auto branch_region = std::max_element(branch_regions.begin(), branch_regions.end(),
                [](const ExPolygon &lhs, const ExPolygon &rhs) {
                    return std::abs(area(lhs)) < std::abs(area(rhs));
                });
            const Lines remaining_boundary = to_lines(*branch_region);
            Point child_center;
            double nearest_distance = std::numeric_limits<double>::max();
            // A child may branch from any already printed arc, not only from the
            // last root ring. Concave surfaces and holes often leave an uncovered
            // pocket beside an earlier parent arc.
            for (const Polyline &parent : arcs) {
                for (const Point &point : parent.points) {
                    for (const Line &line : remaining_boundary) {
                        const double distance = line.distance_to_squared(point);
                        if (distance < nearest_distance) {
                            nearest_distance = distance;
                            child_center = point;
                        }
                    }
                }
            }
            if (nearest_distance > 2.25 * double(spacing) * double(spacing))
                break;

            Polylines child_arcs;
            Polylines child_frontier;
            // The covered mask extends just beyond half a line spacing. Start
            // children slightly farther out so their first arc overlaps the
            // parent while still intersecting the uncovered pocket.
            generate_family(child_center, arcs.back().last_point(), 0.55 * spacing,
                            *branch_region, child_arcs, child_frontier);
            if (child_arcs.empty())
                break;
            append(arcs, std::move(child_arcs));
            append(frontier, std::move(child_frontier));
        }
    }

    append(polylines_out, std::move(arcs));
}

} // namespace Slic3r
