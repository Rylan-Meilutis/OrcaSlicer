#include "FillArcOverhang.hpp"

#include "../ClipperUtils.hpp"
#include "../Geometry.hpp"

#include <algorithm>
#include <array>
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
    Vec2d axis(std::cos(direction.first), std::sin(direction.first));
    const Point midpoint = bbox.center();
    const double ray_length = bbox.size().cast<double>().norm() + 4. * spacing;

    Polyline ray;
    ray.points.emplace_back(midpoint - (axis * ray_length).cast<coord_t>());
    ray.points.emplace_back(midpoint + (axis * ray_length).cast<coord_t>());
    Polylines clipped_ray = intersection_pl(ray, expolygon);
    if (clipped_ray.empty())
        return;

    Point anchor;
    bool has_supported_anchor = false;
    auto supported_anchor_on_ray = [&params](const Polyline &candidate_ray, const Point &ray_midpoint,
                                             Point &supported_anchor) {
        if (params.arc_anchor_regions == nullptr || params.arc_anchor_regions->empty())
            return false;
        const Polylines supported_ray = intersection_pl(Polylines{candidate_ray}, *params.arc_anchor_regions);
        double nearest_center_distance = std::numeric_limits<double>::max();
        bool found = false;
        for (const Polyline &segment : supported_ray) {
            for (const Point &point : {segment.first_point(), segment.last_point()}) {
                const double distance = (point - ray_midpoint).cast<double>().squaredNorm();
                if (distance < nearest_center_distance) {
                    supported_anchor = point;
                    nearest_center_distance = distance;
                    found = true;
                }
            }
        }
        return found;
    };
    has_supported_anchor = supported_anchor_on_ray(ray, midpoint, anchor);
    if (!has_supported_anchor && params.arc_anchor_regions != nullptr &&
        !params.arc_anchor_regions->empty()) {
        // FillBase expresses bridge directions in its line-pattern convention,
        // which is perpendicular to the physical span for some bridge surfaces.
        // Prefer the orthogonal axis when it is the one that actually reaches a
        // lower-layer anchor.
        const Vec2d alternate_axis(-axis.y(), axis.x());
        Polyline alternate_ray;
        alternate_ray.points.emplace_back(midpoint - (alternate_axis * ray_length).cast<coord_t>());
        alternate_ray.points.emplace_back(midpoint + (alternate_axis * ray_length).cast<coord_t>());
        Polylines alternate_clipped_ray = intersection_pl(alternate_ray, expolygon);
        Point alternate_anchor;
        if (!alternate_clipped_ray.empty() &&
            supported_anchor_on_ray(alternate_ray, midpoint, alternate_anchor)) {
            axis = alternate_axis;
            ray = std::move(alternate_ray);
            clipped_ray = std::move(alternate_clipped_ray);
            anchor = alternate_anchor;
            has_supported_anchor = true;
        }
    }
    if (!has_supported_anchor) {
        anchor = clipped_ray.front().first_point();
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
    } else if ((midpoint - anchor).cast<double>().dot(axis) < 0.) {
        // A one-sided overhang may be supported from either side. Always grow
        // from the selected supported wall into the surface interior.
        axis = -axis;
    }
    const Vec2d normal(-axis.y(), axis.x());

    // Keep the center just inside the bridge surface so the smallest arc is
    // anchored by a useful length instead of degenerating to a point.
    const Point center = anchor + (axis * (0.5 * spacing)).cast<coord_t>();
    const double resolution = std::max<double>(scale_(params.resolution), SCALED_EPSILON);
    const bool recursive_fill = params.config != nullptr && params.config->arc_overhang_recursive_fill.value;

    // Emit each family from its supported anchor outwards. Starting the source
    // circle on the anchor side is important for a closed first ring: clipping
    // cannot otherwise move its arbitrary first vertex onto the perimeter.
    auto generate_family = [&](const Point &family_anchor, const Point &family_center, double first_radius,
                               double radius_limit, const ExPolygon &region, size_t max_paths, Polylines &arcs) {
        Point previous = family_anchor;
        const Vec2d anchor_vector = (family_anchor - family_center).cast<double>();
        const double start_angle = anchor_vector.squaredNorm() > 0. ?
            std::atan2(anchor_vector.y(), anchor_vector.x()) : std::atan2(axis.y(), axis.x()) + M_PI;
        const BoundingBox region_bbox = get_extents(region);
        double family_max_radius =
            std::max((region_bbox.min - family_center).cast<double>().norm(),
                     std::max((region_bbox.max - family_center).cast<double>().norm(),
                              std::max((Point(region_bbox.min.x(), region_bbox.max.y()) - family_center).cast<double>().norm(),
                                       (Point(region_bbox.max.x(), region_bbox.min.y()) - family_center).cast<double>().norm()))) +
            spacing;
        if (radius_limit > 0.)
            family_max_radius = std::min(family_max_radius, radius_limit);

        for (double radius = first_radius; radius <= family_max_radius; radius += spacing) {
            const double angle_step = std::clamp(2. * std::acos(std::max(-1., 1. - resolution / radius)),
                                                 Geometry::deg2rad(2.), Geometry::deg2rad(15.));
            Polyline circle;
            const size_t segments = std::max<size_t>(24, size_t(std::ceil(2. * M_PI / angle_step)));
            circle.points.reserve(segments + 1);
            for (size_t i = 0; i <= segments; ++i) {
                const double angle = start_angle + 2. * M_PI * double(i) / double(segments);
                circle.points.emplace_back(family_center.x() + coord_t(std::lround(radius * std::cos(angle))),
                                           family_center.y() + coord_t(std::lround(radius * std::sin(angle))));
            }

            Polylines ring = intersection_pl(circle, region);
            ring.erase(std::remove_if(ring.begin(), ring.end(), [spacing, recursive_fill](const Polyline &path) {
                if (path.length() < 0.25 * spacing)
                    return true;
                if (!recursive_fill || path.points.size() < 3)
                    return false;

                const Point &first = path.first_point();
                const Point &last = path.last_point();
                const double chord_length = (last - first).cast<double>().norm();
                if (chord_length < 4. * spacing)
                    return false;

                const Line chord(first, last);
                double max_sagitta_squared = 0.;
                for (const Point &point : path.points)
                    max_sagitta_squared = std::max(max_sagitta_squared, chord.distance_to_squared(point));
                // A large circle clipped to a narrow strip is mathematically an
                // arc but prints like a straight bridge. Leave that area for a
                // smaller supported child family instead.
                return max_sagitta_squared < std::pow(1.25 * spacing, 2);
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

            const bool path_limit_reached = ordered.size() >= max_paths - std::min(max_paths, arcs.size());
            if (path_limit_reached)
                ordered.resize(max_paths - arcs.size());
            append(arcs, std::move(ordered));
            if (path_limit_reached)
                break;
        }
    };

    Polylines arcs;
    const std::array<Point, 4> bbox_corners{
        bbox.min, bbox.max, Point(bbox.min.x(), bbox.max.y()), Point(bbox.max.x(), bbox.min.y())};
    double axis_min = std::numeric_limits<double>::max();
    double axis_max = std::numeric_limits<double>::lowest();
    double normal_min = std::numeric_limits<double>::max();
    double normal_max = std::numeric_limits<double>::lowest();
    for (const Point &corner : bbox_corners) {
        const Vec2d position = corner.cast<double>();
        axis_min = std::min(axis_min, position.dot(axis));
        axis_max = std::max(axis_max, position.dot(axis));
        normal_min = std::min(normal_min, position.dot(normal));
        normal_max = std::max(normal_max, position.dot(normal));
    }
    const double axis_span = axis_max - axis_min;
    const double normal_span = normal_max - normal_min;
    const bool narrow_region = axis_span <= 8. * spacing && normal_span > 1.5 * axis_span;
    const bool narrow_chord_region = normal_span <= 12. * spacing && axis_span > 1.5 * normal_span;
    const double narrow_radius_limit =
        narrow_region ? std::max(1.5 * spacing, axis_span + spacing) :
        narrow_chord_region ? std::max(1.5 * spacing, normal_span) : 0.;

    // The primary family always covers the complete radial extent. Recursive
    // fill is only responsible for details that concentric primary arches cannot
    // reach, such as pockets behind holes and concave boundaries.
    generate_family(anchor, center, 0.5 * spacing, narrow_radius_limit, expolygon,
                    std::numeric_limits<size_t>::max(), arcs);

    if (narrow_chord_region && has_supported_anchor && params.arc_anchor_regions != nullptr) {
        // When the narrow bridge is supported at both ends, grow a second tight
        // family from the far wall. The two supported families meet near the
        // middle instead of turning the last few rings into shallow full-span
        // chords.
        const Polylines supported_ray = intersection_pl(Polylines{ray}, *params.arc_anchor_regions);
        Point opposite_anchor = anchor;
        double opposite_distance = 0.;
        for (const Polyline &segment : supported_ray) {
            for (const Point &point : {segment.first_point(), segment.last_point()}) {
                const double distance = (point - anchor).cast<double>().norm();
                if (distance > opposite_distance) {
                    opposite_distance = distance;
                    opposite_anchor = point;
                }
            }
        }
        if (opposite_distance > 2. * spacing) {
            Vec2d opposite_axis = (midpoint - opposite_anchor).cast<double>();
            if (opposite_axis.squaredNorm() > 0.) {
                opposite_axis.normalize();
                const Point opposite_center =
                    opposite_anchor + (opposite_axis * (0.5 * spacing)).cast<coord_t>();
                generate_family(opposite_anchor, opposite_center, 0.5 * spacing, narrow_radius_limit,
                                expolygon, std::numeric_limits<size_t>::max(), arcs);
            }
        }
    }

    if (narrow_region && has_supported_anchor) {
        // A long, narrow bridge should be made from several locally anchored
        // arches, not one very large circle whose shallow segment behaves like
        // an ordinary sagging bridge. Overlap neighboring families slightly so
        // this primary pass still covers the strip when recursive refinement is
        // disabled.
        const double family_pitch = std::max(2. * spacing, 1.25 * narrow_radius_limit);
        const size_t family_count = size_t(std::ceil(normal_span / family_pitch));
        for (size_t family_idx = 1; family_idx <= family_count; ++family_idx) {
            for (const double sign : {-1., 1.}) {
                const Point family_midpoint =
                    midpoint + (normal * (sign * family_idx * family_pitch)).cast<coord_t>();
                const double normal_projection = family_midpoint.cast<double>().dot(normal);
                if (normal_projection < normal_min || normal_projection > normal_max)
                    continue;
                Polyline family_ray;
                family_ray.points.emplace_back(family_midpoint - (axis * ray_length).cast<coord_t>());
                family_ray.points.emplace_back(family_midpoint + (axis * ray_length).cast<coord_t>());
                Point family_anchor;
                if (!supported_anchor_on_ray(family_ray, family_midpoint, family_anchor))
                    continue;
                Vec2d family_axis = axis;
                if ((family_midpoint - family_anchor).cast<double>().dot(family_axis) < 0.)
                    family_axis = -family_axis;
                const Point family_center =
                    family_anchor + (family_axis * (0.5 * spacing)).cast<coord_t>();
                generate_family(family_anchor, family_center, 0.5 * spacing, narrow_radius_limit,
                                expolygon, std::numeric_limits<size_t>::max(), arcs);
            }
        }
    }

    Lines printed_lines;
    const auto segment_intersection = [](const Line &first, const Line &second, Point &intersection) {
        const Vec2d p = first.a.cast<double>();
        const Vec2d q = second.a.cast<double>();
        const Vec2d r = (first.b - first.a).cast<double>();
        const Vec2d s = (second.b - second.a).cast<double>();
        const auto cross = [](const Vec2d &lhs, const Vec2d &rhs) {
            return lhs.x() * rhs.y() - lhs.y() * rhs.x();
        };
        const double denominator = cross(r, s);
        if (std::abs(denominator) <= std::numeric_limits<double>::epsilon())
            return false;
        const double first_position = cross(q - p, s) / denominator;
        const double second_position = cross(q - p, r) / denominator;
        constexpr double position_epsilon = 1e-9;
        if (first_position < -position_epsilon || first_position > 1. + position_epsilon ||
            second_position < -position_epsilon || second_position > 1. + position_epsilon)
            return false;
        const Vec2d position = p + std::clamp(first_position, 0., 1.) * r;
        intersection = Point(coord_t(std::lround(position.x())), coord_t(std::lround(position.y())));
        return true;
    };
    const auto trim_at_printed_path =
        [spacing, &printed_lines, &segment_intersection](Polyline path, bool require_anchor) {
        Polyline trimmed;
        if (path.points.size() < 2)
            return trimmed;

        if (require_anchor) {
            double start_distance_squared = std::numeric_limits<double>::max();
            for (const Line &printed : printed_lines)
                start_distance_squared = std::min(
                    start_distance_squared, printed.distance_to_squared(path.first_point()));
            if (start_distance_squared > std::pow(1.1 * spacing, 2))
                return trimmed;
        }

        trimmed.points.emplace_back(path.first_point());
        const double intersection_epsilon_squared =
            double(SCALED_EPSILON) * double(SCALED_EPSILON);
        for (size_t point_idx = 1; point_idx < path.points.size(); ++point_idx) {
            const Line segment(path.points[point_idx - 1], path.points[point_idx]);
            Point nearest_intersection;
            double nearest_distance_squared = std::numeric_limits<double>::max();
            for (const Line &printed : printed_lines) {
                Point intersection;
                if (!segment_intersection(segment, printed, intersection))
                    continue;
                const double distance_squared =
                    (intersection - segment.a).cast<double>().squaredNorm();
                // A path may begin on an existing extrusion. Any subsequent
                // contact terminates it before the nozzle crosses that line.
                if (point_idx == 1 && distance_squared <= intersection_epsilon_squared)
                    continue;
                if (distance_squared < nearest_distance_squared) {
                    nearest_distance_squared = distance_squared;
                    nearest_intersection = intersection;
                }
            }

            if (nearest_distance_squared < std::numeric_limits<double>::max()) {
                if (nearest_intersection != trimmed.last_point())
                    trimmed.points.emplace_back(nearest_intersection);
                break;
            }
            trimmed.points.emplace_back(segment.b);
        }

        if (trimmed.length() < 0.25 * spacing)
            trimmed.points.clear();
        return trimmed;
    };

    if (recursive_fill && !arcs.empty()) {
        // Recursive mode may create more than one primary family. Resolve their
        // overlaps before computing uncovered pockets; otherwise child arches
        // inherit crossed geometry and cannot repair the resulting gaps.
        Polylines noncrossing_arcs;
        noncrossing_arcs.reserve(arcs.size());
        for (Polyline &arc : arcs) {
            Polyline trimmed = trim_at_printed_path(std::move(arc), false);
            if (trimmed.empty())
                continue;
            append(printed_lines, to_lines(trimmed));
            noncrossing_arcs.emplace_back(std::move(trimmed));
        }
        arcs = std::move(noncrossing_arcs);

        // Lay down the largest useful family first, then refine the pockets left
        // beside it in generations. Subtracting only each generation's new
        // coverage avoids the quadratic whole-model union that made the original
        // recursive implementation cancel slicing on complex models.
        const ExPolygons covered = union_ex(offset(arcs, float(0.52 * spacing)));
        ExPolygons remaining = diff_ex(ExPolygons{expolygon}, covered);
        auto prune_remaining = [spacing](ExPolygons &regions) {
            // Polygonized line offsets leave narrow radial slivers between
            // otherwise adjacent arcs. Remove only sub-line-width noise: wider
            // wedges are printable and must remain available to child arches.
            regions = opening_ex(regions, float(0.4 * spacing));
            regions.erase(std::remove_if(regions.begin(), regions.end(), [spacing](const ExPolygon &part) {
                return std::abs(area(part)) < 0.25 * double(spacing) * double(spacing);
            }), regions.end());
            std::sort(regions.begin(), regions.end(), [](const ExPolygon &lhs, const ExPolygon &rhs) {
                return std::abs(area(lhs)) > std::abs(area(rhs));
            });
        };
        prune_remaining(remaining);

        // Recursive fill is a refinement pass. These limits keep pathological
        // meshes bounded while allowing smaller arches to grow from the arches
        // produced by the preceding generation.
        constexpr size_t max_generations = 32;
        constexpr size_t max_child_regions = 256;
        constexpr size_t max_child_paths = 8192;
        size_t child_path_count = 0;
        size_t child_region_count = 0;
        for (size_t generation = 0;
             generation < max_generations && !remaining.empty() &&
             child_region_count < max_child_regions && child_path_count < max_child_paths;
             ++generation) {
            const size_t parent_arc_count = arcs.size();
            const size_t parent_point_count = std::accumulate(
                arcs.begin(), arcs.begin() + parent_arc_count, size_t(0),
                [](size_t count, const Polyline &path) { return count + path.points.size(); });
            const size_t point_stride = std::max<size_t>(1, (parent_point_count + 8191) / 8192);
            Polylines generation_arcs;

            const size_t regions_this_generation =
                std::min(remaining.size(), max_child_regions - child_region_count);
            for (size_t region_idx = 0;
                 region_idx < regions_this_generation &&
                 child_path_count + generation_arcs.size() < max_child_paths;
                 ++region_idx, ++child_region_count) {
                ExPolygons branch_regions = intersection_ex(
                    offset_ex(ExPolygons{remaining[region_idx]}, float(0.55 * spacing)), ExPolygons{expolygon});
                if (branch_regions.empty())
                    continue;
                auto branch_region = std::max_element(branch_regions.begin(), branch_regions.end(),
                    [](const ExPolygon &lhs, const ExPolygon &rhs) {
                        return std::abs(area(lhs)) < std::abs(area(rhs));
                    });
                // Locate the parent beside the actual uncovered pocket. The
                // expanded branch may share the model's outer contour, which can
                // otherwise select an unrelated parent point and produce only a
                // tiny clipped circle at the wrong side of the region.
                const Lines remaining_boundary = to_lines(remaining[region_idx]);
                Point child_anchor;
                Point pocket_boundary;
                double nearest_distance = std::numeric_limits<double>::max();
                size_t point_idx = 0;
                for (size_t parent_idx = 0; parent_idx < parent_arc_count; ++parent_idx) {
                    const Polyline &parent = arcs[parent_idx];
                    for (const Point &point : parent.points) {
                        if (point_idx++ % point_stride != 0)
                            continue;
                        for (const Line &line : remaining_boundary) {
                            const double distance = line.distance_to_squared(point);
                            if (distance < nearest_distance) {
                                nearest_distance = distance;
                                child_anchor = point;
                                line.distance_to_squared(point, &pocket_boundary);
                            }
                        }
                    }
                }
                if (nearest_distance > 2.25 * double(spacing) * double(spacing))
                    continue;

                Polylines child_arcs;
                // The covered mask extends just beyond half a line spacing. Start
                // children slightly farther out, with the circle center directed
                // into the pocket. Its first ring passes through child_anchor, so
                // extrusion begins on an already printed parent rather than in
                // the unsupported pocket.
                Vec2d child_direction = (pocket_boundary - child_anchor).cast<double>();
                if (child_direction.squaredNorm() == 0.)
                    child_direction = (remaining[region_idx].contour.centroid() - child_anchor).cast<double>();
                if (child_direction.squaredNorm() == 0.)
                    continue;
                child_direction.normalize();
                const double child_radius = 0.55 * spacing;
                const Point child_center =
                    child_anchor + (child_direction * child_radius).cast<coord_t>();
                generate_family(child_anchor, child_center, child_radius, 0., *branch_region,
                                max_child_paths - child_path_count - generation_arcs.size(), child_arcs);
                for (Polyline &child_arc : child_arcs) {
                    Polyline trimmed = trim_at_printed_path(std::move(child_arc), true);
                    if (trimmed.empty())
                        continue;
                    append(printed_lines, to_lines(trimmed));
                    generation_arcs.emplace_back(std::move(trimmed));
                }
            }

            if (generation_arcs.empty())
                break;

            child_path_count += generation_arcs.size();
            const ExPolygons newly_covered = union_ex(offset(generation_arcs, float(0.52 * spacing)));
            append(arcs, std::move(generation_arcs));
            remaining = diff_ex(remaining, newly_covered);
            prune_remaining(remaining);
        }
    }

    append(polylines_out, std::move(arcs));
}

} // namespace Slic3r
