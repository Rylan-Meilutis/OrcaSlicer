#include "FillArcOverhang.hpp"

#include "../AABBTreeLines.hpp"
#include "../ClipperUtils.hpp"
#include "../Geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <unordered_map>

namespace Slic3r {

namespace {

struct ArcGridCell {
    int64_t x;
    int64_t y;

    bool operator==(const ArcGridCell &other) const { return x == other.x && y == other.y; }
};

struct ArcGridCellHash {
    size_t operator()(const ArcGridCell &cell) const
    {
        const size_t x_hash = std::hash<int64_t>{}(cell.x);
        const size_t y_hash = std::hash<int64_t>{}(cell.y);
        return x_hash ^ (y_hash + 0x9e3779b9 + (x_hash << 6) + (x_hash >> 2));
    }
};

// Recursive arc fill may contain thousands of short segments. Comparing every
// new segment with the complete set of previously printed segments made the
// overlap-removal pass quadratic. This grid keeps the same exact intersection
// tests while limiting them to lines whose bounding boxes are nearby.
class PrintedArcIndex {
public:
    PrintedArcIndex(Lines &lines, coord_t cell_size)
        : m_lines(lines), m_cell_size(std::max<coord_t>(1, cell_size))
    {}

    void add(const Polyline &path)
    {
        for (size_t point_idx = 1; point_idx < path.points.size(); ++point_idx)
            add(Line(path.points[point_idx - 1], path.points[point_idx]));
    }

    template<typename Visitor>
    void visit(const Line &query, coord_t padding, Visitor &&visitor)
    {
        visit_indices(
            query, padding,
            [this, &visitor](size_t line_idx) {
                visitor(m_lines[line_idx]);
            });
    }

private:
    int64_t cell(coord_t coordinate) const
    {
        if (coordinate >= 0)
            return coordinate / m_cell_size;
        return -int64_t((uint64_t(-(coordinate + 1)) + 1 + m_cell_size - 1) / m_cell_size);
    }

    template<typename Visitor>
    void visit_indices(const Line &query, coord_t padding, Visitor &&visitor)
    {
        if (++m_query_id == 0) {
            std::fill(m_seen.begin(), m_seen.end(), 0);
            ++m_query_id;
        }

        const int64_t min_x = cell(std::min(query.a.x(), query.b.x()) - padding);
        const int64_t max_x = cell(std::max(query.a.x(), query.b.x()) + padding);
        const int64_t min_y = cell(std::min(query.a.y(), query.b.y()) - padding);
        const int64_t max_y = cell(std::max(query.a.y(), query.b.y()) + padding);
        for (int64_t x = min_x; x <= max_x; ++x) {
            for (int64_t y = min_y; y <= max_y; ++y) {
                const auto found = m_cells.find({x, y});
                if (found == m_cells.end())
                    continue;
                for (const size_t line_idx : found->second) {
                    if (m_seen[line_idx] == m_query_id)
                        continue;
                    m_seen[line_idx] = m_query_id;
                    visitor(line_idx);
                }
            }
        }
    }

    void add(const Line &line)
    {
        const size_t line_idx = m_lines.size();
        m_lines.emplace_back(line);
        m_seen.emplace_back(0);

        const int64_t min_x = cell(std::min(line.a.x(), line.b.x()));
        const int64_t max_x = cell(std::max(line.a.x(), line.b.x()));
        const int64_t min_y = cell(std::min(line.a.y(), line.b.y()));
        const int64_t max_y = cell(std::max(line.a.y(), line.b.y()));
        for (int64_t x = min_x; x <= max_x; ++x)
            for (int64_t y = min_y; y <= max_y; ++y)
                m_cells[{x, y}].emplace_back(line_idx);
    }

    Lines                                                   &m_lines;
    coord_t                                                  m_cell_size;
    std::unordered_map<ArcGridCell, std::vector<size_t>,
                       ArcGridCellHash>                       m_cells;
    std::vector<uint32_t>                                    m_seen;
    uint32_t                                                 m_query_id{0};
};

} // namespace

void FillArcOverhang::_fill_surface_single(const FillParams              &params,
                                           unsigned int,
                                           const std::pair<float, Point> &direction,
                                           ExPolygon                      expolygon,
                                           Polylines                     &polylines_out)
{
    if (expolygon.empty() || params.density <= 0.f)
        return;
    const BoundingBox bbox = get_extents(expolygon);
    const coord_t line_width = std::max<coord_t>(
        1, scale_(params.flow.width() > 0.f ? params.flow.width() : this->spacing));
    const coord_t nominal_spacing =
        std::max<coord_t>(1, coord_t(std::lround(scale_(this->spacing) / params.density)));
    const double overlap_percent =
        params.config == nullptr ? 0. :
        std::clamp(params.config->arc_overhang_overlap.value, 0., 99.);
    const coord_t line_overlap = coord_t(std::lround(
        double(line_width) * overlap_percent / 100.));
    const coord_t spacing = std::max<coord_t>(1, nominal_spacing - line_overlap);
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
    const double anchor_embed_depth = std::max(1.25 * double(line_width), double(spacing));
    const auto supported_anchor_in_regions =
        [anchor_embed_depth, &expolygon](const Polyline &candidate_ray,
                             const Point &ray_midpoint,
                             const ExPolygons *anchor_regions,
                             Point &supported_anchor) {
        if (anchor_regions == nullptr || anchor_regions->empty())
            return false;
        const ExPolygons usable_anchor_regions =
            intersection_ex(*anchor_regions, ExPolygons{expolygon});
        if (usable_anchor_regions.empty())
            return false;
        const Polylines supported_ray =
            intersection_pl(Polylines{candidate_ray},
                            usable_anchor_regions);
        double nearest_center_distance = std::numeric_limits<double>::max();
        bool found = false;
        for (const Polyline &segment : supported_ray) {
            if (segment.points.size() < 2)
                continue;
            const Point &first = segment.first_point();
            const Point &last  = segment.last_point();
            const double first_distance = (first - ray_midpoint).cast<double>().squaredNorm();
            const double last_distance  = (last - ray_midpoint).cast<double>().squaredNorm();
            const Point &transition = first_distance < last_distance ? first : last;
            const Point &interior   = first_distance < last_distance ? last : first;
            const Vec2d support_direction = (interior - transition).cast<double>();
            const double available_depth = support_direction.norm();
            if (available_depth <= SCALED_EPSILON)
                continue;

            const double distance = std::min(first_distance, last_distance);
            if (distance >= nearest_center_distance)
                continue;

            // A circle tangent to the support boundary has only point contact
            // and starts like an unsupported extrusion. Move the seed into the
            // retained wall/lower-layer footprint so the first arc has a solid
            // supported lead-in before it grows into free air.
            const double embed_depth = std::min(anchor_embed_depth, 0.8 * available_depth);
            const Vec2d embedded =
                transition.cast<double>() + support_direction * (embed_depth / available_depth);
            supported_anchor = Point(coord_t(std::lround(embedded.x())),
                                     coord_t(std::lround(embedded.y())));
            nearest_center_distance = distance;
            found = true;
        }
        return found;
    };
    const auto select_anchor_axis =
        [&](const ExPolygons *anchor_regions) {
        if (supported_anchor_in_regions(
                ray, midpoint, anchor_regions, anchor))
            return true;

        // FillBase expresses bridge directions in its line-pattern convention,
        // which is perpendicular to the physical span for some bridge surfaces.
        // Check the orthogonal axis against the same preferred support class
        // before considering a broader fallback mask.
        const Vec2d alternate_axis(-axis.y(), axis.x());
        Polyline alternate_ray;
        alternate_ray.points.emplace_back(midpoint - (alternate_axis * ray_length).cast<coord_t>());
        alternate_ray.points.emplace_back(midpoint + (alternate_axis * ray_length).cast<coord_t>());
        Polylines alternate_clipped_ray = intersection_pl(alternate_ray, expolygon);
        Point alternate_anchor;
        if (!alternate_clipped_ray.empty() &&
            supported_anchor_in_regions(
                alternate_ray, midpoint, anchor_regions,
                alternate_anchor)) {
            axis = alternate_axis;
            ray = std::move(alternate_ray);
            clipped_ray = std::move(alternate_clipped_ray);
            anchor = alternate_anchor;
            return true;
        }
        return false;
    };
    has_supported_anchor =
        select_anchor_axis(params.arc_root_anchor_regions);
    if (!has_supported_anchor)
        has_supported_anchor =
            select_anchor_axis(params.arc_anchor_regions);
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

    // Match the reference construction: all primary circles share the
    // selected point on the supported start region as their center. The anchor
    // selection above already embeds that point into the preceding layer, so a
    // one-width first circle has a real supported runway without shifting the
    // entire concentric family away from its support.
    const Point center = anchor;
    const double resolution = std::max<double>(scale_(params.resolution), SCALED_EPSILON);
    const bool recursive_fill = params.config != nullptr && params.config->arc_overhang_recursive_fill.value;
    std::vector<BoundingBox> obstacle_region_bboxes;
    if (params.arc_obstacle_regions != nullptr) {
        obstacle_region_bboxes.reserve(
            params.arc_obstacle_regions->size());
        for (const ExPolygon &obstacle :
             *params.arc_obstacle_regions)
            obstacle_region_bboxes.emplace_back(get_extents(obstacle));
    }

    // Emit each family from its supported anchor outwards. Starting the source
    // circle on the anchor side is important for a closed first ring: clipping
    // cannot otherwise move its arbitrary first vertex onto the perimeter.
    auto generate_family = [&](const Point &family_anchor, const Point &family_center, double first_radius,
                               double radius_limit, const ExPolygon &region, size_t max_paths,
                               bool reject_shallow_paths, Polylines &arcs,
                               ExPolygons *swept_area,
                               Polylines *outer_frontier) {
        Point previous = family_anchor;
        const Vec2d anchor_vector = (family_anchor - family_center).cast<double>();
        const double start_angle = anchor_vector.squaredNorm() > 0. ?
            std::atan2(anchor_vector.y(), anchor_vector.x()) : std::atan2(axis.y(), axis.x()) + M_PI;
        // Bounding-box corners may be far outside a narrow or diagonal pocket.
        // Circles beyond the farthest polygon vertex cannot intersect that
        // pocket (apart from their chord tolerance), so avoid generating and
        // clipping guaranteed-empty rings.
        double family_max_radius_squared = 0.;
        const auto include_polygon_radius =
            [&family_center, &family_max_radius_squared](
                const Polygon &polygon) {
                for (const Point &point : polygon.points)
                    family_max_radius_squared = std::max(
                        family_max_radius_squared,
                        (point - family_center)
                            .cast<double>()
                            .squaredNorm());
            };
        include_polygon_radius(region.contour);
        for (const Polygon &hole : region.holes)
            include_polygon_radius(hole);
        double family_max_radius =
            std::sqrt(family_max_radius_squared) + resolution +
            SCALED_EPSILON;
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
            ring.erase(std::remove_if(ring.begin(), ring.end(), [spacing, line_width, reject_shallow_paths](const Polyline &path) {
                // Reject sub-bead clipping remnants before they can become
                // anchors for a child family. Removing them only in the final
                // sanitizer leaves later paths depending on geometry which is
                // no longer emitted.
                if (path.length() <
                    std::max(0.25 * double(spacing),
                             0.75 * double(line_width)))
                    return true;
                if (!reject_shallow_paths || path.points.size() < 3)
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
                continue;

            if (swept_area != nullptr) {
                Polygon disk;
                disk.points.assign(circle.points.begin(),
                                   circle.points.end() - 1);
                *swept_area = intersection_ex(
                    ExPolygons{ExPolygon(std::move(disk))},
                    ExPolygons{region});
            }
            if (outer_frontier != nullptr)
                *outer_frontier = ring;

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
    // A 0.4 mm nozzle commonly leaves bridge slots around 4-5 mm wide. At
    // that width, clipping an expanding circle still produces long chords
    // that are almost straight, even though the slot is unmistakably narrow
    // compared with its length. Keep using bounded-radius translated arches
    // through sixteen pitches; the aspect-ratio check prevents this from
    // changing compact or generally wide overhangs.
    const bool narrow_region = axis_span <= 16. * spacing && normal_span > 1.5 * axis_span;
    const bool narrow_chord_region = normal_span <= 32. * spacing &&
                                     axis_span > 1.5 * normal_span;
    const double narrow_radius_limit =
        narrow_region ? std::max(1.5 * spacing, axis_span + spacing) :
        narrow_chord_region ? std::max(1.5 * spacing, normal_span) : 0.;

    const bool chained_primary = (narrow_region || narrow_chord_region) &&
                                 has_supported_anchor &&
                                 narrow_radius_limit > 0.;
    if (chained_primary) {
        const bool has_retained_obstacle_paths =
            params.arc_obstacle_paths != nullptr &&
            !params.arc_obstacle_paths->empty();
        // Very large circles clipped to a narrow span become nearly straight
        // bridge lines. Tile the span with parallel, fixed-radius circular
        // arches instead. Translating one arch by one line width preserves its
        // curvature and cannot cross the preceding arch.
        auto generate_parallel_arcs = [&](const Point &start_anchor, Vec2d family_axis,
                                          double target_distance, double phase_offset) {
            if (family_axis.squaredNorm() == 0.)
                return;
            family_axis.normalize();
            Vec2d family_normal(-family_axis.y(), family_axis.x());
            double transverse_min = std::numeric_limits<double>::max();
            double transverse_max = std::numeric_limits<double>::lowest();
            for (const Point &corner : bbox_corners) {
                const double projection =
                    corner.cast<double>().dot(family_normal);
                transverse_min = std::min(transverse_min, projection);
                transverse_max = std::max(transverse_max, projection);
            }
            const double half_width = 0.5 * (transverse_max - transverse_min);
            const double radius = std::max(half_width + spacing,
                                           1.25 * half_width);
            const double half_angle = std::asin(
                std::clamp(half_width / radius, 0., 1.));
            const double sagitta =
                radius - std::sqrt(std::max(0., radius * radius -
                                                 half_width * half_width));
            const double angle_step = std::clamp(
                2. * std::acos(std::max(-1., 1. - resolution / radius)),
                Geometry::deg2rad(2.), Geometry::deg2rad(15.));
            const size_t segments = std::max<size_t>(
                8, size_t(std::ceil(2. * half_angle / angle_step)));
            const double start_projection =
                start_anchor.cast<double>().dot(family_axis);
            const double normal_midpoint =
                0.5 * (transverse_min + transverse_max);
            // Parallel arches are translations of the same curve, so their
            // centerline pitch may use the configured arc spacing directly.
            // This also preserves the arc-overhang overlap option: increasing
            // overlap reduces spacing and therefore emits more arches.
            const double path_pitch = double(spacing);
            phase_offset = std::clamp(
                phase_offset, 0.,
                std::nextafter(path_pitch, 0.));
            // Keep every generated arch inside the target span. Overshooting
            // by the sagitta made clipping turn the last several arches into
            // nearly straight fragments (or two-point lines) at the far end.
            const double usable_distance =
                std::max(
                    0.,
                    target_distance - sagitta -
                        phase_offset);
            const size_t path_count = std::max<size_t>(
                1, size_t(std::floor(usable_distance / path_pitch)) + 1);
            Point previous = start_anchor;
            for (size_t path_idx = 0;
                 path_idx < path_count;
                 ++path_idx) {
                const double peak_projection =
                    start_projection +
                    sagitta + phase_offset +
                    path_idx * path_pitch;
                const double center_projection = peak_projection - radius;
                Polyline arch;
                arch.points.reserve(segments + 1);
                for (size_t segment_idx = 0; segment_idx <= segments; ++segment_idx) {
                    const double angle =
                        -half_angle + 2. * half_angle *
                                          double(segment_idx) / double(segments);
                    const Vec2d position =
                        family_axis * (center_projection + radius * std::cos(angle)) +
                        family_normal * (normal_midpoint + radius * std::sin(angle));
                    arch.points.emplace_back(
                        coord_t(std::lround(position.x())),
                        coord_t(std::lround(position.y())));
                }

                Polylines clipped = intersection_pl(Polylines{arch}, expolygon);
                if (clipped.empty())
                    continue;
                size_t nearest_idx = 0;
                bool reverse = false;
                double nearest_distance = std::numeric_limits<double>::max();
                for (size_t idx = 0; idx < clipped.size(); ++idx) {
                    const double first_distance =
                        (clipped[idx].first_point() - previous).cast<double>().squaredNorm();
                    const double last_distance =
                        (clipped[idx].last_point() - previous).cast<double>().squaredNorm();
                    if (std::min(first_distance, last_distance) < nearest_distance) {
                        nearest_distance = std::min(first_distance, last_distance);
                        nearest_idx = idx;
                        reverse = last_distance < first_distance;
                    }
                }
                Polyline path = std::move(clipped[nearest_idx]);
                if (path_idx == 0) {
                    reverse =
                        (path.last_point() - start_anchor).cast<double>().squaredNorm() <
                        (path.first_point() - start_anchor).cast<double>().squaredNorm();
                }
                if (reverse)
                    path.reverse();
                previous = path.last_point();
                arcs.emplace_back(std::move(path));
            }
        };

        if (narrow_region) {
            // The bridge direction crosses the short dimension. Growing one
            // ever-larger circle along the long dimension eventually produces
            // almost straight clipped chords. Print a curved transverse arch at
            // the supported midpoint, then advance in both directions. The
            // validation pass below retains these only when both ends attach or
            // the complete arch follows an already deposited bead.
            Vec2d positive_chain_axis = normal;
            if (positive_chain_axis.squaredNorm() == 0.)
                positive_chain_axis = Vec2d(-axis.y(), axis.x());
            positive_chain_axis.normalize();
            double positive_distance = 0.;
            double negative_distance = 0.;
            for (const Point &corner : bbox_corners) {
                const double distance =
                    (corner - anchor).cast<double>().dot(positive_chain_axis);
                positive_distance = std::max(positive_distance, distance);
                negative_distance = std::max(negative_distance, -distance);
            }
            // When a retained wall exists the supported concentric frontier
            // below is translated instead. Avoid constructing this generic
            // one-ended family only to discard it afterwards.
            if (!has_retained_obstacle_paths) {
                generate_parallel_arcs(
                    anchor, positive_chain_axis,
                    positive_distance, 0.);
                generate_parallel_arcs(
                    anchor, -positive_chain_axis,
                    negative_distance, 0.);
            }
        } else {
            Point opposite_anchor = anchor;
            double opposite_distance = 0.;
            const ExPolygons *opposite_anchor_regions =
                params.arc_root_anchor_regions != nullptr &&
                !params.arc_root_anchor_regions->empty() ?
                    params.arc_root_anchor_regions :
                    params.arc_anchor_regions;
            const Polylines support_on_axis =
                intersection_pl(Polylines{ray}, *opposite_anchor_regions);
            for (const Polyline &segment : support_on_axis) {
                for (const Point &point : {segment.first_point(), segment.last_point()}) {
                    const double distance = (point - anchor).cast<double>().norm();
                    if (distance > opposite_distance) {
                        opposite_distance = distance;
                        opposite_anchor = point;
                    }
                }
            }

            Vec2d chain_axis = (midpoint - anchor).cast<double>();
            if (chain_axis.squaredNorm() == 0.)
            chain_axis = axis;
            if (opposite_distance > std::max(2. * spacing, 0.5 * axis_span)) {
                // Meet in the middle when both sides are supported. Each half is
                // independently printable and its first extrusion starts on its
                // own perimeter wall. Align the second half to the first
                // family's pitch so their last arches meet without leaving a
                // one-pitch center gap or depositing two nearly coincident
                // arches.
                chain_axis.normalize();
                const Vec2d family_normal(
                    -chain_axis.y(), chain_axis.x());
                double transverse_min =
                    std::numeric_limits<double>::max();
                double transverse_max =
                    std::numeric_limits<double>::lowest();
                for (const Point &corner : bbox_corners) {
                    const double projection =
                        corner.cast<double>().dot(
                            family_normal);
                    transverse_min =
                        std::min(transverse_min, projection);
                    transverse_max =
                        std::max(transverse_max, projection);
                }
                const double half_width =
                    0.5 * (transverse_max -
                           transverse_min);
                const double radius =
                    std::max(
                        half_width + spacing,
                        1.25 * half_width);
                const double sagitta =
                    radius -
                    std::sqrt(std::max(
                        0.,
                        radius * radius -
                            half_width * half_width));
                double opposite_phase =
                    std::fmod(
                        opposite_distance -
                            2. * sagitta,
                        double(spacing));
                if (opposite_phase < 0.)
                    opposite_phase += spacing;

                generate_parallel_arcs(
                    anchor, chain_axis,
                    0.5 * opposite_distance, 0.);
                generate_parallel_arcs(opposite_anchor,
                                       (midpoint - opposite_anchor).cast<double>(),
                                       0.5 * opposite_distance,
                                       opposite_phase);
            } else {
                chain_axis.normalize();
                double target_distance = 0.;
                for (const Point &corner : bbox_corners)
                    target_distance = std::max(
                        target_distance,
                        (corner - anchor).cast<double>().dot(chain_axis));
                generate_parallel_arcs(
                    anchor, chain_axis,
                    target_distance, 0.);
            }
        }

        if (has_retained_obstacle_paths) {
            // A translated narrow family may initially have only its departure
            // end on a retained wall. Seed it with the reference concentric
            // family instead of accepting that one-ended extrusion. Validation
            // prints these rings first; deferred translated arches are retried
            // afterwards and may then attach to the completed backbone.
            Polylines translated_arcs = std::move(arcs);
            arcs.clear();
            generate_family(
                anchor, center, 0.5 * double(spacing), 0., expolygon,
                std::numeric_limits<size_t>::max(), true, arcs,
                nullptr, nullptr);
            if (narrow_region && !arcs.empty()) {
                // The generic translated arches span the slot from one wall
                // to the other. With support on only one wall their first
                // member necessarily ends in free air and is correctly
                // rejected below, leaving the long ends of the slot empty.
                //
                // Instead translate the two halves of the largest accepted
                // concentric arch. Each new half is one bead pitch beside its
                // predecessor, starts on the retained wall, and retains the
                // same substantial curvature. This builds outward in both
                // directions without crossings or a one-ended seed.
                size_t frontier_idx = 0;
                double frontier_extent =
                    std::numeric_limits<double>::lowest();
                const double anchor_projection =
                    anchor.cast<double>().dot(axis);
                for (size_t path_idx = 0;
                     path_idx < arcs.size(); ++path_idx) {
                    const double first_projection =
                        arcs[path_idx].first_point()
                            .cast<double>()
                            .dot(axis);
                    const double last_projection =
                        arcs[path_idx].last_point()
                            .cast<double>()
                            .dot(axis);
                    // A usable backbone half leaves and returns to the
                    // retained wall. Outer clipped circle fragments with one
                    // endpoint on the far edge are precisely the shallow,
                    // one-ended paths this construction is replacing.
                    if (std::abs(first_projection -
                                 anchor_projection) >
                            double(line_width) ||
                        std::abs(last_projection -
                                 anchor_projection) >
                            double(line_width))
                        continue;
                    double maximum =
                        std::numeric_limits<double>::lowest();
                    double minimum =
                        std::numeric_limits<double>::max();
                    for (const Point &point :
                         arcs[path_idx].points) {
                        const double projection =
                            point.cast<double>().dot(axis);
                        maximum = std::max(maximum, projection);
                        minimum = std::min(minimum, projection);
                    }
                    if (maximum - minimum >
                        frontier_extent) {
                        frontier_extent = maximum - minimum;
                        frontier_idx = path_idx;
                    }
                }

                const Polyline &frontier = arcs[frontier_idx];
                size_t peak_idx = 0;
                double peak_projection =
                    std::numeric_limits<double>::lowest();
                for (size_t point_idx = 0;
                     point_idx < frontier.points.size();
                     ++point_idx) {
                    const double projection =
                        frontier.points[point_idx]
                            .cast<double>()
                            .dot(axis);
                    if (projection > peak_projection) {
                        peak_projection = projection;
                        peak_idx = point_idx;
                    }
                }

                if (peak_idx > 0 &&
                    peak_idx + 1 <
                        frontier.points.size()) {
                    std::array<Polyline, 2> halves;
                    halves[0].points.insert(
                        halves[0].points.end(),
                        frontier.points.begin(),
                        frontier.points.begin() +
                            peak_idx + 1);
                    halves[1].points.insert(
                        halves[1].points.end(),
                        frontier.points.begin() +
                            peak_idx,
                        frontier.points.end());

                    for (Polyline &half : halves) {
                        const double mean_normal =
                            std::accumulate(
                                half.points.begin(),
                                half.points.end(), 0.,
                                [&normal](double sum,
                                          const Point &point) {
                                    return sum +
                                        point.cast<double>()
                                            .dot(normal);
                                }) /
                            double(half.points.size());
                        const double center_normal =
                            center.cast<double>().dot(normal);
                        const double side =
                            mean_normal >= center_normal ?
                                1. : -1.;
                        if ((half.last_point() - anchor)
                                .cast<double>()
                                .squaredNorm() <
                            (half.first_point() - anchor)
                                .cast<double>()
                                .squaredNorm())
                            half.reverse();

                        const size_t maximum_steps =
                            size_t(std::ceil(
                                normal_span /
                                double(spacing))) + 2;
                        for (size_t step = 1;
                             step <= maximum_steps;
                             ++step) {
                            const Vec2d offset =
                                side * double(step) *
                                double(spacing) * normal;
                            Polyline shifted = half;
                            for (Point &point :
                                 shifted.points) {
                                const Vec2d position =
                                    point.cast<double>() +
                                    offset;
                                point = Point(
                                    coord_t(std::lround(
                                        position.x())),
                                    coord_t(std::lround(
                                        position.y())));
                            }
                            Polylines clipped =
                                intersection_pl(
                                    Polylines{
                                        std::move(shifted)},
                                    expolygon);
                            if (clipped.empty())
                                break;
                            auto longest =
                                std::max_element(
                                    clipped.begin(),
                                    clipped.end(),
                                    [](const Polyline &lhs,
                                       const Polyline &rhs) {
                                        return lhs.length() <
                                               rhs.length();
                                    });
                            if (longest->length() <
                                0.25 * spacing)
                                break;
                            if ((longest->last_point() -
                                 anchor)
                                    .cast<double>()
                                    .squaredNorm() <
                                (longest->first_point() -
                                 anchor)
                                    .cast<double>()
                                    .squaredNorm())
                                longest->reverse();
                            arcs.emplace_back(
                                std::move(*longest));
                        }
                    }
                } else {
                    append(arcs,
                           std::move(translated_arcs));
                }
            } else {
                append(arcs, std::move(translated_arcs));
            }
        }
    } else {
        // Outside narrow spans, one complete directly anchored family supplies
        // the large-arc backbone before any recursive detail is considered.
        generate_family(anchor, center, spacing, 0., expolygon,
                        std::numeric_limits<size_t>::max(), false, arcs,
                        nullptr, nullptr);
    }

    Lines printed_lines;
    PrintedArcIndex printed_index(printed_lines, 4 * spacing);
    if (params.arc_prior_paths != nullptr) {
        for (const Polyline &path : *params.arc_prior_paths)
            printed_index.add(path);
    }
    Lines obstacle_lines;
    PrintedArcIndex obstacle_index(obstacle_lines, 4 * spacing);
    if (params.arc_obstacle_paths != nullptr) {
        for (const Polyline &obstacle : *params.arc_obstacle_paths)
            obstacle_index.add(obstacle);
    }
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
    const coord_t gcode_coordinate_step = scale_(0.001);
    const auto rounded_proper_intersection =
        [gcode_coordinate_step](const Line &first,
                                const Line &second,
                                Point &intersection) {
            const auto rounded =
                [gcode_coordinate_step](const Point &point) {
                    return Point(
                        coord_t(std::llround(
                            double(point.x()) /
                            double(gcode_coordinate_step))) *
                            gcode_coordinate_step,
                        coord_t(std::llround(
                            double(point.y()) /
                            double(gcode_coordinate_step))) *
                            gcode_coordinate_step);
                };
            const Point first_a = rounded(first.a);
            const Point first_b = rounded(first.b);
            const Point second_a = rounded(second.a);
            const Point second_b = rounded(second.b);
            const Vec2d p = first_a.cast<double>();
            const Vec2d q = second_a.cast<double>();
            const Vec2d r = (first_b - first_a).cast<double>();
            const Vec2d s = (second_b - second_a).cast<double>();
            const auto cross =
                [](const Vec2d &lhs, const Vec2d &rhs) {
                    return lhs.x() * rhs.y() -
                           lhs.y() * rhs.x();
                };
            const double denominator = cross(r, s);
            if (std::abs(denominator) <=
                std::numeric_limits<double>::epsilon())
                return false;
            const double first_position =
                cross(q - p, s) / denominator;
            const double second_position =
                cross(q - p, r) / denominator;
            constexpr double endpoint_epsilon = 1e-8;
            if (first_position <= endpoint_epsilon ||
                first_position >= 1. - endpoint_epsilon ||
                second_position <= endpoint_epsilon ||
                second_position >= 1. - endpoint_epsilon)
                return false;
            const Vec2d position = p + first_position * r;
            intersection = Point(
                coord_t(std::lround(position.x())),
                coord_t(std::lround(position.y())));
            return true;
        };

    // A retained perimeter is both a collision boundary and a valid support
    // bead. Clip out only its protected center while preserving overlap onto
    // both sides of the bead. Unlike splitting at a shared centerline point,
    // the resulting distinct endpoints cannot be chained back together across
    // the perimeter when extrusion entities are ordered for G-code.
    const auto split_at_obstacle_paths =
        [spacing, &params, &obstacle_region_bboxes](Polyline path) {
        if (path.points.size() < 2)
            return Polylines{};
        if (params.arc_obstacle_regions == nullptr ||
            params.arc_obstacle_regions->empty())
            return Polylines{std::move(path)};

        const BoundingBox path_bbox = get_extents(path);
        ExPolygons relevant_obstacles;
        for (size_t obstacle_idx = 0;
             obstacle_idx < obstacle_region_bboxes.size();
             ++obstacle_idx) {
            if (path_bbox.overlap(
                    obstacle_region_bboxes[obstacle_idx]))
                relevant_obstacles.emplace_back(
                    (*params.arc_obstacle_regions)[obstacle_idx]);
        }
        if (relevant_obstacles.empty())
            return Polylines{std::move(path)};

        Polylines fragments = diff_pl(
            Polylines{path}, relevant_obstacles);
        const auto path_position = [&path](const Point &point) {
            double best_distance_squared = std::numeric_limits<double>::max();
            double best_position = 0.;
            double cumulative_length = 0.;
            for (size_t point_idx = 1; point_idx < path.points.size(); ++point_idx) {
                const Vec2d start = path.points[point_idx - 1].cast<double>();
                const Vec2d delta =
                    (path.points[point_idx] - path.points[point_idx - 1]).cast<double>();
                const double length_squared = delta.squaredNorm();
                if (length_squared <= 0.)
                    continue;
                const double segment_position = std::clamp(
                    (point.cast<double>() - start).dot(delta) / length_squared,
                    0., 1.);
                const Vec2d nearest = start + segment_position * delta;
                const double distance_squared =
                    (point.cast<double>() - nearest).squaredNorm();
                if (distance_squared < best_distance_squared) {
                    best_distance_squared = distance_squared;
                    best_position =
                        cumulative_length + segment_position * std::sqrt(length_squared);
                }
                cumulative_length += std::sqrt(length_squared);
            }
            return best_position;
        };

        std::vector<std::pair<double, Polyline>> ordered;
        ordered.reserve(fragments.size());
        for (Polyline &fragment : fragments) {
            if (fragment.length() < 0.25 * spacing)
                continue;
            double start_position = path_position(fragment.first_point());
            double end_position = path_position(fragment.last_point());
            if (end_position < start_position) {
                fragment.reverse();
                std::swap(start_position, end_position);
            }
            ordered.emplace_back(start_position, std::move(fragment));
        }
        std::sort(ordered.begin(), ordered.end(),
                  [](const auto &lhs, const auto &rhs) {
                      return lhs.first < rhs.first;
                  });

        fragments.clear();
        fragments.reserve(ordered.size());
        for (auto &entry : ordered)
            fragments.emplace_back(std::move(entry.second));
        return fragments;
    };

    // A child circle can be clipped into a curve which never crosses its
    // parent centerline but follows it closely enough to deposit a second bead
    // over the first one. Centerline-intersection tests cannot detect that
    // case. Reject a sustained near-coincident run while still allowing the
    // short endpoint contact needed to anchor one arc to another.
    const auto trim_leading_sustained_retrace =
        [line_width, spacing](Polyline &path,
                              PrintedArcIndex &occupied_index) {
        if (path.points.size() < 2)
            return;

        const coord_t duplicate_clearance =
            std::max<coord_t>(
                SCALED_EPSILON,
                coord_t(std::lround(std::max(
                    0.20 * double(line_width),
                    0.55 * double(spacing)))));
        const double sample_step =
            std::max<double>(SCALED_EPSILON, 0.15 * line_width);
        const double maximum_duplicate_run = 0.65 * double(line_width);
        double duplicate_run = 0.;
        double path_position = 0.;
        for (size_t point_idx = 1;
             point_idx < path.points.size(); ++point_idx) {
            const Vec2d start =
                path.points[point_idx - 1].cast<double>();
            const Vec2d delta =
                (path.points[point_idx] -
                 path.points[point_idx - 1]).cast<double>();
            const double length = delta.norm();
            if (length <= 0.)
                continue;
            const size_t samples = std::max<size_t>(
                1, size_t(std::ceil(length / sample_step)));
            const double interval = length / double(samples);
            for (size_t sample_idx = 0;
                 sample_idx < samples; ++sample_idx) {
                const Vec2d position =
                    start + delta *
                                ((double(sample_idx) + 0.5) /
                                 double(samples));
                const Point sample(
                    coord_t(std::lround(position.x())),
                    coord_t(std::lround(position.y())));
                double nearest_distance_squared =
                    std::numeric_limits<double>::max();
                occupied_index.visit(
                    Line(sample, sample), duplicate_clearance,
                    [&sample, &delta, &nearest_distance_squared](
                        const Line &printed) {
                        constexpr double minimum_parallel_dot_squared =
                            0.82;
                        const Vec2d printed_delta =
                            (printed.b - printed.a).cast<double>();
                        const double direction_norm_squared =
                            delta.squaredNorm() *
                            printed_delta.squaredNorm();
                        const double direction_dot =
                            delta.dot(printed_delta);
                        if (direction_norm_squared <= 0. ||
                            direction_dot * direction_dot <
                                minimum_parallel_dot_squared *
                                    direction_norm_squared)
                            return;
                        nearest_distance_squared = std::min(
                            nearest_distance_squared,
                            printed.distance_to_squared(sample));
                    });
                if (nearest_distance_squared <
                    double(duplicate_clearance) *
                        double(duplicate_clearance)) {
                    duplicate_run += interval;
                } else {
                    if (duplicate_run >= maximum_duplicate_run) {
                        // Retain the useful portion after the child leaves its
                        // parent bead. The transition point is still within one
                        // bead width of that parent and therefore remains a
                        // supported extrusion start.
                        path.clip_start(path_position);
                    }
                    return;
                }
                path_position += interval;
            }
        }
        if (duplicate_run >= maximum_duplicate_run)
            path.points.clear();
    };

    const auto trim_sustained_retrace =
        [&trim_leading_sustained_retrace](
            Polyline &path, PrintedArcIndex &occupied_index) {
        trim_leading_sustained_retrace(path, occupied_index);
        if (path.points.size() < 2)
            return;

        // A fragment may converge onto a previously emitted family at its far
        // anchor. Preserve a crossing or sharply converging bead contact, but
        // remove a terminal run which follows the earlier centerline in
        // parallel. Reversing lets the same leading-run logic find the exact
        // departure point without introducing a second clipping algorithm.
        path.reverse();
        trim_leading_sustained_retrace(path, occupied_index);
        if (path.points.size() >= 2)
            path.reverse();
    };

    const auto trim_at_printed_path =
        [spacing, line_width, &params,
         &printed_index, &obstacle_index, &segment_intersection,
         &rounded_proper_intersection](Polyline path,
                                       bool require_anchor) {
        Polyline trimmed;
        if (path.points.size() < 2)
            return trimmed;

        if (require_anchor) {
            // Anchoring is a bead-contact requirement, not a path-pitch
            // requirement. Spacing may exceed the extrusion width for free-air
            // arcs; using it here can accept two centerlines whose deposited
            // beads still have an air gap between them.
            const coord_t anchor_distance = line_width;
            const bool has_retained_perimeters =
                params.arc_obstacle_paths != nullptr &&
                !params.arc_obstacle_paths->empty();
            const auto is_anchored =
                [&params, &printed_index, &obstacle_index,
                 anchor_distance, has_retained_perimeters](const Point &point) {
                double printed_distance_squared =
                    std::numeric_limits<double>::max();
                const auto measure_printed_distance =
                    [&printed_distance_squared, &point](const Line &printed) {
                    printed_distance_squared = std::min(
                        printed_distance_squared,
                        printed.distance_to_squared(point));
                };
                printed_index.visit(
                    Line(point, point), anchor_distance,
                    measure_printed_distance);

                double obstacle_distance_squared =
                    std::numeric_limits<double>::max();
                const auto measure_obstacle_distance =
                    [&obstacle_distance_squared, &point](const Line &obstacle) {
                    obstacle_distance_squared = std::min(
                        obstacle_distance_squared,
                        obstacle.distance_to_squared(point));
                };
                obstacle_index.visit(
                    Line(point, point), anchor_distance,
                    measure_obstacle_distance);

                const double maximum_distance_squared =
                    double(anchor_distance) * double(anchor_distance);
                if (printed_distance_squared <= maximum_distance_squared ||
                    obstacle_distance_squared <= maximum_distance_squared)
                    return true;

                // Root regions are the intersection with the actual previous
                // layer. They remain a valid physical runway even when this
                // layer also retains perimeter walls. Refusing that support
                // made the first closed ring fail its arbitrary start-point
                // check; every larger ring then lost its predecessor and most
                // of the overhang surface was left empty.
                if (params.arc_root_anchor_regions != nullptr &&
                    std::any_of(
                        params.arc_root_anchor_regions->begin(),
                        params.arc_root_anchor_regions->end(),
                        [&point](const ExPolygon &region) {
                            return region.contains(point);
                        }))
                    return true;

                // A lower-layer support region is sufficient when there is no
                // same-layer wall to attach to. When retained perimeters do
                // exist, however, accepting the broad geometric mask alone can
                // seed a narrow bridge near its midpoint in free space. In that
                // case require physical bead contact with a perimeter or with
                // an arc that has already been printed.
                if (has_retained_perimeters)
                    return false;
                return params.arc_anchor_regions != nullptr &&
                       std::any_of(params.arc_anchor_regions->begin(),
                                   params.arc_anchor_regions->end(),
                                   [&point](const ExPolygon &region) {
                                       return region.contains(point);
                                   });
            };

            if (!is_anchored(path.first_point())) {
                if (!is_anchored(path.last_point()))
                    return trimmed;
                path.reverse();
            }
        }

        trimmed.points.emplace_back(path.first_point());
        const double intersection_epsilon_squared =
            double(SCALED_EPSILON) * double(SCALED_EPSILON);
        for (size_t point_idx = 1; point_idx < path.points.size(); ++point_idx) {
            const Line segment(path.points[point_idx - 1], path.points[point_idx]);
            Point nearest_intersection;
            double nearest_distance_squared = std::numeric_limits<double>::max();
            const auto record_occupied_intersection =
                [&](const Line &occupied) {
                // Recursive children may start a few integer coordinates apart
                // on the same parent tangent. Nearly parallel first chords can
                // then overlap for most of a bead before their numerically
                // unstable line intersection is detected farther downstream.
                // Treat that duplicate lead-in as occupied instead of laying
                // two extrusions on top of one another.
                if (point_idx == 1 &&
                    occupied.distance_to_squared(segment.a) <=
                        intersection_epsilon_squared &&
                    occupied.distance_to_squared(segment.b) <=
                        std::pow(0.10 * line_width, 2) &&
                    segment.length() >= 0.25 * spacing) {
                    nearest_distance_squared = 0.;
                    nearest_intersection = segment.a;
                    return;
                }
                Point intersection;
                if (!rounded_proper_intersection(
                        segment, occupied, intersection)) {
                    if (!segment_intersection(
                            segment, occupied, intersection))
                        return;
                }
                const double distance_squared =
                    (intersection - segment.a).cast<double>().squaredNorm();
                // A path may begin on an existing extrusion. Any subsequent
                // contact terminates it before the nozzle crosses that line.
                if (point_idx == 1 && distance_squared <= intersection_epsilon_squared)
                    return;
                if (distance_squared < nearest_distance_squared) {
                    nearest_distance_squared = distance_squared;
                    nearest_intersection = intersection;
                }
            };
            printed_index.visit(
                segment, SCALED_EPSILON, record_occupied_intersection);
            // Polygonized perimeter cores are deliberately narrower than the
            // deposited bead so an arc can overlap its edge for anchoring.
            // They are not, however, a complete collision test: short wall
            // segments and safety-offset quantization may leave no core at all.
            // Check the retained centerlines as well and terminate at the
            // first contact instead of allowing the arc to cross a wall that
            // will already have been printed on this layer.
            obstacle_index.visit(
                segment, SCALED_EPSILON, record_occupied_intersection);

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
    const auto is_shallow_arc = [spacing](const Polyline& path) {
        if (path.points.size() < 2)
            return false;
        const Line chord(path.first_point(), path.last_point());
        // A two-point extrusion is a line, not an arch. It may only survive
        // below when its entire length is already supported.
        if (path.points.size() == 2)
            return true;
        if (chord.length() < 2. * spacing)
            return false;
        // The paper calls out large-radius fragments as unstable because they
        // behave like straight bridges. Require ten percent bow on short
        // nozzle-scale spans, increasing to a full pitch on longer chords.
        // Applying the full-pitch requirement at only four pitches discarded
        // the complete 2.4 mm bridge family even though those arches still
        // have substantial curvature.
        const double minimum_sagitta = std::min(double(spacing), 0.10 * chord.length());
        double max_sagitta_squared   = 0.;
        for (const Point& point : path.points)
            max_sagitta_squared = std::max(max_sagitta_squared, chord.distance_to_squared(point));
        return max_sagitta_squared < minimum_sagitta * minimum_sagitta;
    };
    const auto tighten_shallow_arc =
        [spacing, resolution, &expolygon,
         &center](const Polyline &path) {
            Polyline tightened;
            if (path.points.size() < 2)
                return tightened;

            const Point &first = path.first_point();
            const Point &last = path.last_point();
            const Vec2d chord =
                (last - first).cast<double>();
            const double chord_length = chord.norm();
            if (chord_length < 2. * spacing)
                return tightened;

            const Vec2d midpoint =
                0.5 * (first.cast<double>() +
                       last.cast<double>());
            Vec2d normal(-chord.y(), chord.x());
            normal /= chord_length;

            // Preserve the side on which the original large circle bowed.
            // When clipping leaves only a two-point chord, use the side away
            // from the primary center, which is the unsupported corner the
            // rejected fragment was intended to cover.
            double signed_original_sagitta = 0.;
            for (const Point &point : path.points) {
                const double signed_distance =
                    (point.cast<double>() - midpoint)
                        .dot(normal);
                if (std::abs(signed_distance) >
                    std::abs(signed_original_sagitta))
                    signed_original_sagitta =
                        signed_distance;
            }
            if (std::abs(signed_original_sagitta) <=
                SCALED_EPSILON)
                signed_original_sagitta =
                    (midpoint - center.cast<double>())
                            .dot(normal) >= 0. ?
                        1. : -1.;
            const double side =
                signed_original_sagitta >= 0. ? 1. : -1.;

            // Match the curvature test with a small safety margin. Bending
            // every rejected lane much farther than necessary makes adjacent
            // replacements fan apart and opens wedge-shaped holes between
            // them.
            const double minimum_sagitta =
                std::min(
                    double(spacing),
                    0.10 * chord_length);
            const double sagitta = std::min(
                0.25 * chord_length,
                1.05 * minimum_sagitta);
            const double radius =
                chord_length * chord_length /
                    (8. * sagitta) +
                0.5 * sagitta;
            const double center_distance =
                radius - sagitta;
            const Vec2d arc_center =
                midpoint - side * center_distance * normal;
            const double start_angle = std::atan2(
                first.y() - arc_center.y(),
                first.x() - arc_center.x());
            const double end_angle = std::atan2(
                last.y() - arc_center.y(),
                last.x() - arc_center.x());
            double sweep = end_angle - start_angle;
            while (sweep <= -M_PI)
                sweep += 2. * M_PI;
            while (sweep > M_PI)
                sweep -= 2. * M_PI;

            // The selected minor arc must bow to the same side as the
            // rejected source fragment.
            const double middle_angle =
                start_angle + 0.5 * sweep;
            const Vec2d arc_middle =
                arc_center +
                radius * Vec2d(std::cos(middle_angle),
                               std::sin(middle_angle));
            if ((arc_middle - midpoint).dot(normal) *
                    side <
                0.) {
                sweep += sweep > 0. ?
                    -2. * M_PI : 2. * M_PI;
            }
            if (std::abs(sweep) > M_PI + 1e-6)
                return tightened;

            const double angle_step = std::clamp(
                2. * std::acos(std::max(
                    -1.,
                    1. - resolution / radius)),
                Geometry::deg2rad(2.),
                Geometry::deg2rad(15.));
            const size_t segments =
                std::max<size_t>(
                    3,
                    size_t(std::ceil(
                        std::abs(sweep) /
                        angle_step)));
            tightened.points.reserve(segments + 1);
            for (size_t segment_idx = 0;
                 segment_idx <= segments;
                 ++segment_idx) {
                const double angle =
                    start_angle +
                    sweep * double(segment_idx) /
                        double(segments);
                tightened.points.emplace_back(
                    coord_t(std::lround(
                        arc_center.x() +
                        radius * std::cos(angle))),
                    coord_t(std::lround(
                        arc_center.y() +
                        radius * std::sin(angle))));
            }
            tightened.points.front() = first;
            tightened.points.back() = last;

            // Do not create a more curved path by pushing it outside the
            // bridge surface. A clipped replacement would lose one of its
            // supported endpoints and belongs in recursive handling instead.
            const Polylines inside =
                intersection_pl(
                    Polylines{tightened}, expolygon);
            const double inside_length =
                std::accumulate(
                    inside.begin(), inside.end(), 0.,
                    [](double total,
                       const Polyline &part) {
                        return total + part.length();
                    });
            if (inside_length <
                tightened.length() -
                    0.10 * double(spacing))
                tightened.points.clear();
            return tightened;
        };
    const auto is_shallow_child_arc =
        [spacing, &is_shallow_arc](const Polyline &path) {
        if (is_shallow_arc(path))
            return true;
        const Line chord(path.first_point(), path.last_point());
        if (chord.length() < 4. * spacing)
            return false;
        double max_sagitta_squared = 0.;
        for (const Point &point : path.points)
            max_sagitta_squared = std::max(
                max_sagitta_squared,
                chord.distance_to_squared(point));
        // Recursive families are the replacement for the large family's
        // nearly straight remnants. A long child must bow by at least one
        // complete pitch or another exposed center must fill this region.
        return max_sagitta_squared <
               double(spacing) * double(spacing);
    };
    const auto is_fully_supported = [&params, spacing](const Polyline& path) {
        if (params.arc_anchor_regions == nullptr || params.arc_anchor_regions->empty())
            return false;
        const Polylines unsupported =
            diff_pl(Polylines{path}, *params.arc_anchor_regions);
        return std::all_of(
            unsupported.begin(), unsupported.end(),
            [spacing](const Polyline &fragment) {
                return fragment.length() < 0.25 * spacing;
            });
    };
    // Polygonized circles approximate the requested radius within the fill
    // resolution. Include that geometric tolerance, but require genuine bead
    // overlap rather than accepting another complete pitch of unsupported
    // travel.
    const coord_t lateral_support_distance =
        std::min(line_width, spacing) +
        coord_t(std::ceil(resolution));
    const auto point_has_root_support =
        [&params](const Point &point) {
            return params.arc_root_anchor_regions != nullptr &&
                   std::any_of(
                       params.arc_root_anchor_regions->begin(),
                       params.arc_root_anchor_regions->end(),
                       [&point](const ExPolygon &region) {
                           return region.contains(point);
                       });
        };
    const auto point_has_lateral_support =
        [lateral_support_distance, &params, &printed_index,
         &obstacle_index](const Point &point) {
        const auto point_in_regions =
            [&point](const ExPolygons *regions) {
                return regions != nullptr &&
                       std::any_of(regions->begin(), regions->end(),
                                   [&point](const ExPolygon &region) {
                                       return region.contains(point);
                                   });
            };
        if (point_in_regions(params.arc_root_anchor_regions) ||
            point_in_regions(params.arc_anchor_regions))
            return true;

        double nearest_distance_squared =
            std::numeric_limits<double>::max();
        const auto measure_distance =
            [&nearest_distance_squared, &point](const Line &support) {
                nearest_distance_squared = std::min(
                    nearest_distance_squared,
                    support.distance_to_squared(point));
            };
        printed_index.visit(
            Line(point, point), lateral_support_distance, measure_distance);
        obstacle_index.visit(
            Line(point, point), lateral_support_distance, measure_distance);
        return nearest_distance_squared <=
               double(lateral_support_distance) *
                   double(lateral_support_distance);
    };
    const coord_t supported_lead_length =
        std::max<coord_t>(1, coord_t(std::ceil(0.75 * double(line_width))));
    const coord_t support_sample_step =
        std::max<coord_t>(1, coord_t(std::ceil(0.20 * double(line_width))));
    const auto has_supported_origin = [&point_has_lateral_support, supported_lead_length, support_sample_step](const Polyline& path) {
        // A single coincident vertex is not an anchor. Require a useful initial
        // length to remain over the previous-layer footprint or alongside an
        // already deposited bead, so extrusion pressure and adhesion are
        // established before the arc enters free air.
        if (path.points.size() < 2)
            return false;
        coord_t remaining = supported_lead_length;
        for (size_t point_idx = 1;
             point_idx < path.points.size() && remaining > 0; ++point_idx) {
            const Point &start = path.points[point_idx - 1];
            const Point &end = path.points[point_idx];
            const Vec2d delta = (end - start).cast<double>();
            const double segment_length = delta.norm();
            if (segment_length <= 0.)
                continue;
            const double checked_length =
                std::min<double>(segment_length, remaining);
            const size_t samples = std::max<size_t>(
                1, size_t(std::ceil(
                       checked_length / double(support_sample_step))));
            for (size_t sample_idx = 0; sample_idx <= samples; ++sample_idx) {
                const double distance = checked_length * double(sample_idx) / double(samples);
                const Vec2d position  = start.cast<double>() + delta * (distance / segment_length);
                const Point sample(coord_t(std::lround(position.x())), coord_t(std::lround(position.y())));
                if (!point_has_lateral_support(sample))
                    return false;
            }
            remaining -= coord_t(std::lround(checked_length));
        }
        // A path shorter than the target lead is acceptable only because every
        // point of that complete path was checked above.
        return true;
    };
    const auto has_supported_ends = [&has_supported_origin](const Polyline& path) {
        if (!has_supported_origin(path))
            return false;
        Polyline reversed = path;
        reversed.reverse();
        return has_supported_origin(reversed);
    };
    const auto has_supported_endpoint = [&has_supported_origin, &point_has_lateral_support](const Polyline& path) {
        // A primary family already has a completed concentric backbone. Its
        // departure still needs a runway, while bead contact at the arrival is
        // enough to close the arc onto a perimeter or a completed sibling.
        return has_supported_origin(path) && point_has_lateral_support(path.last_point());
    };
    const auto is_fully_laterally_supported = [&point_has_lateral_support, support_sample_step](const Polyline& path) {
        if (path.points.size() < 2)
            return false;
        for (size_t point_idx = 1; point_idx < path.points.size(); ++point_idx) {
            const Point& start          = path.points[point_idx - 1];
            const Point& end            = path.points[point_idx];
            const Vec2d delta           = (end - start).cast<double>();
            const double segment_length = delta.norm();
            if (segment_length <= 0.)
                continue;
            const size_t samples = std::max<size_t>(1, size_t(std::ceil(segment_length / double(support_sample_step))));
            for (size_t sample_idx = 0; sample_idx <= samples; ++sample_idx) {
                const Vec2d position = start.cast<double>() + delta * (double(sample_idx) / double(samples));
                const Point sample(
                    coord_t(std::lround(position.x())),
                    coord_t(std::lround(position.y())));
                if (!point_has_lateral_support(sample))
                    return false;
            }
        }
        return true;
    };
    const auto point_has_emitted_bead_support =
        [lateral_support_distance, &params, &printed_index,
         &obstacle_index](const Point &point) {
        if (params.arc_root_anchor_regions != nullptr &&
            std::any_of(
                params.arc_root_anchor_regions->begin(),
                params.arc_root_anchor_regions->end(),
                [&point](const ExPolygon &region) {
                    return region.contains(point);
                }))
            return true;

        double nearest_distance_squared =
            std::numeric_limits<double>::max();
        const auto measure_distance =
            [&nearest_distance_squared, &point](const Line &support) {
                nearest_distance_squared = std::min(
                    nearest_distance_squared,
                    support.distance_to_squared(point));
            };
        printed_index.visit(
            Line(point, point), lateral_support_distance,
            measure_distance);
        obstacle_index.visit(
            Line(point, point), lateral_support_distance,
            measure_distance);
        return nearest_distance_squared <=
               double(lateral_support_distance) *
                   double(lateral_support_distance);
    };
    const auto is_fully_emitted_bead_supported =
        [&point_has_emitted_bead_support,
         support_sample_step](const Polyline &path) {
        if (path.points.size() < 2)
            return false;
        for (size_t point_idx = 1;
             point_idx < path.points.size(); ++point_idx) {
            const Point &start = path.points[point_idx - 1];
            const Point &end   = path.points[point_idx];
            const Vec2d delta  = (end - start).cast<double>();
            const double segment_length = delta.norm();
            if (segment_length <= 0.)
                continue;
            const size_t samples = std::max<size_t>(
                1, size_t(std::ceil(
                       segment_length /
                       double(support_sample_step))));
            for (size_t sample_idx = 0;
                 sample_idx <= samples; ++sample_idx) {
                const Vec2d position =
                    start.cast<double>() +
                    delta * (double(sample_idx) /
                             double(samples));
                if (!point_has_emitted_bead_support(
                        Point(
                            coord_t(std::lround(position.x())),
                            coord_t(std::lround(position.y())))))
                    return false;
            }
        }
        return true;
    };
    const coord_t stable_branch_lead_length =
        std::max<coord_t>(
            1,
            coord_t(std::ceil(
                1.50 * double(line_width))));
    const auto has_stable_emitted_branch_lead =
        [&point_has_emitted_bead_support,
         support_sample_step,
         stable_branch_lead_length](const Polyline &path) {
            if (path.points.size() < 2)
                return false;
            coord_t remaining = stable_branch_lead_length;
            for (size_t point_idx = 1;
                 point_idx < path.points.size() &&
                 remaining > 0;
                 ++point_idx) {
                const Point &start =
                    path.points[point_idx - 1];
                const Point &end =
                    path.points[point_idx];
                const Vec2d delta =
                    (end - start).cast<double>();
                const double segment_length = delta.norm();
                if (segment_length <= 0.)
                    continue;
                const double checked_length =
                    std::min<double>(
                        segment_length, remaining);
                const size_t samples =
                    std::max<size_t>(
                        1,
                        size_t(std::ceil(
                            checked_length /
                            double(support_sample_step))));
                for (size_t sample_idx = 0;
                     sample_idx <= samples;
                     ++sample_idx) {
                    const Vec2d position =
                        start.cast<double>() +
                        delta *
                            (checked_length *
                             double(sample_idx) /
                             double(samples) /
                             segment_length);
                    if (!point_has_emitted_bead_support(
                            Point(
                                coord_t(std::lround(position.x())),
                                coord_t(std::lround(position.y())))))
                        return false;
                }
                remaining -= coord_t(
                    std::lround(checked_length));
            }
            return true;
        };
    const double minimum_free_air_path_length =
        1.5 * double(line_width);
    const auto is_short_free_air_path =
        [minimum_free_air_path_length,
         &is_fully_emitted_bead_supported](const Polyline &path) {
        return path.length() < minimum_free_air_path_length &&
               !is_fully_emitted_bead_supported(path);
    };
    const auto split_from_supported_point =
        [&point_has_lateral_support, &point_has_root_support,
         &params, spacing](Polyline path) {
        const bool first_supported =
            point_has_lateral_support(path.first_point());
        const bool last_supported =
            point_has_lateral_support(path.last_point());
        const bool has_root_context =
            params.arc_root_anchor_regions != nullptr &&
            !params.arc_root_anchor_regions->empty();
        const bool root_endpoint_supported =
            point_has_root_support(path.first_point()) ||
            point_has_root_support(path.last_point());
        if (path.points.size() < 3 ||
            (first_supported && last_supported &&
             (!has_root_context || root_endpoint_supported)))
            return Polylines{std::move(path)};

        const auto find_supported_idx =
            [&path](const auto &is_supported) {
                for (size_t point_idx = 1;
                     point_idx + 1 < path.points.size(); ++point_idx)
                    if (is_supported(path.points[point_idx]))
                        return point_idx;
                return path.points.size();
            };
        // When the preceding layer's actual perimeter footprint is known,
        // start the root family there. A broader supported-slice mask is valid
        // for later arcs, but choosing it first can place the initial bead one
        // complete line width away from the perimeter that must establish the
        // free-air shape.
        size_t supported_idx =
            has_root_context ?
                find_supported_idx(point_has_root_support) :
                path.points.size();
        if (supported_idx == path.points.size())
            supported_idx =
                find_supported_idx(point_has_lateral_support);
        if (supported_idx == path.points.size())
            return Polylines{std::move(path)};

        Polylines result;
        if (path.is_closed()) {
            // Closed rings have no meaningful endpoint. Rotate their arbitrary
            // seam to the point that is physically supported.
            Polyline rotated;
            rotated.points.reserve(path.points.size());
            for (size_t point_idx = supported_idx;
                 point_idx + 1 < path.points.size(); ++point_idx)
                rotated.points.emplace_back(path.points[point_idx]);
            for (size_t point_idx = 0; point_idx <= supported_idx; ++point_idx)
                rotated.points.emplace_back(path.points[point_idx]);
            if (rotated.length() >= 0.25 * spacing)
                result.emplace_back(std::move(rotated));
            return result;
        }

        // An open arc that reaches support in its middle is printable as two
        // paths leaving that shared anchor in opposite directions. Keeping it
        // as one path would force one half to begin in free air.
        Polyline forward;
        forward.points.insert(
            forward.points.end(),
            path.points.begin() + supported_idx, path.points.end());
        if (forward.length() >= 0.25 * spacing)
            result.emplace_back(std::move(forward));

        Polyline backward;
        backward.points.reserve(supported_idx + 1);
        for (size_t point_idx = supported_idx + 1; point_idx > 0; --point_idx)
            backward.points.emplace_back(path.points[point_idx - 1]);
        if (backward.length() >= 0.25 * spacing)
            result.emplace_back(std::move(backward));
        return result;
    };
    const auto family_path_radius =
        [](const Polyline &path, const Point &family_center) {
        double radius = 0.;
        for (const Point &point : path.points)
            radius = std::max(
                radius,
                (point - family_center).cast<double>().norm());
        return radius;
    };
    const auto select_outer_frontier = [&family_path_radius, spacing](const Polylines& family, const Point& family_center) {
        Polylines frontier;
        double outer_radius = 0.;
        for (const Polyline &path : family)
            outer_radius = std::max(
                outer_radius,
                family_path_radius(path, family_center));
        for (const Polyline &path : family) {
            if (family_path_radius(path, family_center) >= outer_radius - 0.25 * spacing)
                frontier.emplace_back(path);
        }
        return frontier;
    };
    const auto swept_disk = [resolution](const Point& family_center, double radius, const ExPolygon& region) {
        ExPolygons swept;
        if (radius <= 0.)
            return swept;
        const double angle_step = std::clamp(2. * std::acos(std::max(-1., 1. - resolution / radius)), Geometry::deg2rad(2.),
                                             Geometry::deg2rad(15.));
        const size_t segments   = std::max<size_t>(24, size_t(std::ceil(2. * M_PI / angle_step)));
        Polygon disk;
        disk.points.reserve(segments);
        for (size_t segment_idx = 0; segment_idx < segments; ++segment_idx) {
            const double angle = 2. * M_PI * double(segment_idx) / double(segments);
            disk.points.emplace_back(family_center.x() + coord_t(std::lround(radius * std::cos(angle))),
                                     family_center.y() + coord_t(std::lround(radius * std::sin(angle))));
        }
        return intersection_ex(ExPolygons{ExPolygon(std::move(disk))}, ExPolygons{region});
    };

    double nominal_primary_outer_radius = 0.;
    for (const Polyline& path : arcs)
        nominal_primary_outer_radius = std::max(nominal_primary_outer_radius, family_path_radius(path, center));

    const auto extend_to_nearby_support = [&expolygon, &printed_lines, &obstacle_lines, line_width, spacing, &point_has_lateral_support,
                                           &segment_intersection](Polyline& path) {
        if (path.points.size() < 2 || point_has_lateral_support(path.last_point()))
            return true;

        const Point end                 = path.last_point();
        double nearest_distance_squared = std::numeric_limits<double>::max();
        Vec2d nearest_point;
        const auto consider_support = [&end, &nearest_distance_squared, &nearest_point](const Line& support) {
            const Vec2d start           = support.a.cast<double>();
            const Vec2d delta           = (support.b - support.a).cast<double>();
            const double length_squared = delta.squaredNorm();
            if (length_squared <= 0.)
                return;
            const double position         = std::clamp((end.cast<double>() - start).dot(delta) / length_squared, 0., 1.);
            const Vec2d projected         = start + position * delta;
            const double distance_squared = (end.cast<double>() - projected).squaredNorm();
            if (distance_squared < nearest_distance_squared) {
                nearest_distance_squared = distance_squared;
                nearest_point            = projected;
            }
        };
        for (const Line& support : obstacle_lines)
            consider_support(support);
        for (const Line& support : printed_lines)
            consider_support(support);

        const double nearest_distance = std::sqrt(nearest_distance_squared);
        if (!std::isfinite(nearest_distance) || nearest_distance > 2.25 * double(spacing))
            return false;

        const Vec2d direction         = (nearest_point - end.cast<double>()) / std::max(nearest_distance, 1.);
        const double extension_length = std::max(0., nearest_distance - 0.50 * double(line_width));
        const Vec2d extended_position = end.cast<double>() + extension_length * direction;
        Point extended_end(coord_t(std::lround(extended_position.x())), coord_t(std::lround(extended_position.y())));
        if (extended_end == end)
            return false;

        const Line extension_line(end, extended_end);
        double nearest_intersection_squared = (extended_end - end).cast<double>().squaredNorm();
        const auto stop_at_first_support    = [&extension_line, &end, &extended_end, &nearest_intersection_squared,
                                            &segment_intersection](const Line& support) {
            Point intersection;
            if (!segment_intersection(extension_line, support, intersection))
                return;
            const double distance_squared = (intersection - end).cast<double>().squaredNorm();
            if (distance_squared <= double(SCALED_EPSILON) * double(SCALED_EPSILON) || distance_squared >= nearest_intersection_squared)
                return;
            nearest_intersection_squared = distance_squared;
            extended_end                 = intersection;
        };
        for (const Line& support : obstacle_lines)
            stop_at_first_support(support);
        for (const Line& support : printed_lines)
            stop_at_first_support(support);

        const Polyline extension(Points{end, extended_end});
        const Polylines inside     = intersection_pl(Polylines{extension}, expolygon);
        const double inside_length = std::accumulate(inside.begin(), inside.end(), 0.,
                                                     [](double total, const Polyline& part) { return total + part.length(); });
        if (inside_length < extension.length() - 0.10 * double(line_width))
            return false;

        path.points.emplace_back(extended_end);
        return point_has_lateral_support(path.last_point());
    };

    if (!arcs.empty()) {
        // Resolve intersections for every mode. Complex clipping may split even
        // a primary family into fragments that cross after ordering; allowing
        // those crossings produces duplicate extrusion and nozzle collisions.
        // Recursive modes additionally rely on this clean backbone when finding
        // their remaining uncovered pockets.
        Polylines pending_arcs;
        pending_arcs.reserve(arcs.size());
        const bool validate_primary_anchors = has_supported_anchor ||
                                              (params.arc_obstacle_paths != nullptr && !params.arc_obstacle_paths->empty()) ||
                                              (params.arc_prior_paths != nullptr && !params.arc_prior_paths->empty());
        const bool require_supported_ends = params.arc_obstacle_paths != nullptr && !params.arc_obstacle_paths->empty();
        for (Polyline& arc : arcs)
            append(pending_arcs, split_at_obstacle_paths(std::move(arc)));

        Polylines noncrossing_arcs;
        noncrossing_arcs.reserve(pending_arcs.size());
        // A concentric ring may be split into several disconnected pieces by
        // holes and retained walls. Some pieces are not printable when first
        // encountered but become supported after another piece reaches them
        // from a perimeter. Retry deferred primary pieces before resorting to
        // small recursive families. This preserves the large-arc backbone and
        // prevents printable regions from being permanently discarded merely
        // because of their initial traversal order.
        while (!pending_arcs.empty()) {
            Polylines deferred;
            deferred.reserve(pending_arcs.size());
            const size_t accepted_before = noncrossing_arcs.size();
            for (Polyline& fragment : pending_arcs) {
                Polylines anchored_fragments = chained_primary ? Polylines{std::move(fragment)} :
                                                                 split_from_supported_point(std::move(fragment));
                for (Polyline& anchored_fragment : anchored_fragments) {
                    Polyline trimmed = trim_at_printed_path(anchored_fragment, validate_primary_anchors);
                    if (trimmed.empty()) {
                        deferred.emplace_back(std::move(anchored_fragment));
                        continue;
                    }
                    // A short straight lead-in is harmless while it remains fully
                    // on an existing wall or lower-layer footprint. Never keep a
                    // shallow first path once it enters free air: the bridge must
                    // leave its anchor as an arch so later paths inherit a stable
                    // curved foundation.
                    if (validate_primary_anchors && is_shallow_arc(trimmed) && !is_fully_supported(trimmed) &&
                        !is_fully_laterally_supported(trimmed)) {
                        Polyline tightened =
                            tighten_shallow_arc(trimmed);
                        if (!tightened.empty())
                            tightened =
                                trim_at_printed_path(
                                    std::move(tightened),
                                    validate_primary_anchors);
                        if (tightened.empty() ||
                            (is_shallow_arc(tightened) &&
                             !is_fully_supported(tightened) &&
                             !is_fully_laterally_supported(
                                 tightened))) {
                            deferred.emplace_back(
                                std::move(
                                    anchored_fragment));
                            continue;
                        }
                        trimmed = std::move(tightened);
                    }
                    if (validate_primary_anchors &&
                        is_short_free_air_path(trimmed)) {
                        deferred.emplace_back(
                            std::move(anchored_fragment));
                        continue;
                    }
                    if (validate_primary_anchors && chained_primary &&
                        require_supported_ends &&
                        is_shallow_child_arc(trimmed) &&
                        !is_fully_supported(trimmed) &&
                        !is_fully_laterally_supported(trimmed)) {
                        deferred.emplace_back(
                            std::move(anchored_fragment));
                        continue;
                    }
                    if (validate_primary_anchors && require_supported_ends && !has_supported_endpoint(trimmed) &&
                        !extend_to_nearby_support(trimmed)) {
                        deferred.emplace_back(std::move(anchored_fragment));
                        continue;
                    }
                    if (validate_primary_anchors && !require_supported_ends && !has_supported_origin(trimmed)) {
                        deferred.emplace_back(std::move(anchored_fragment));
                        continue;
                    }
                    trim_sustained_retrace(
                        trimmed, printed_index);
                    if (trimmed.length() <
                        std::max(
                            0.25 * double(spacing),
                            0.75 * double(line_width)))
                        continue;
                    printed_index.add(trimmed);
                    noncrossing_arcs.emplace_back(std::move(trimmed));
                }
            }
            if (noncrossing_arcs.size() == accepted_before)
                break;
            pending_arcs = std::move(deferred);
        }
        arcs = std::move(noncrossing_arcs);
    }
    const size_t primary_arc_count = arcs.size();
    if (recursive_fill && !arcs.empty()) {
        // Follow the breadth-first construction used by the reference arc
        // overhang implementation: finish the directly anchored concentric
        // family, subtract its deposited area, then grow complete child
        // families from already printed parent arcs.  In particular, never
        // expand a whole uncovered pocket.  Doing so made a child's clipping
        // region overlap unrelated arcs and perimeters, allowing geometrically
        // valid circles to turn into crossing G-code after clipping.
        // The reference algorithm treats a one-width-pitch concentric family
        // as a filled clipped disk before applying breadth-first children.
        // Retain any substantial holes caused by our additional collision and
        // anchor validation, while filtering the bead-pitch slivers that are
        // already filled by the family's deposited width.
        ExPolygons primary_covered;
        if (!chained_primary) {
            primary_covered = swept_disk(center, nominal_primary_outer_radius, expolygon);
        } else {
            primary_covered = union_ex(offset(arcs, float(0.52 * line_width)));
        }
        ExPolygons actual_uncovered = diff_ex(ExPolygons{expolygon}, union_ex(offset(arcs, float(0.52 * line_width))));
        // Collapse isolated bead-pitch seams before adding actual coverage
        // defects to the reference swept-disk work queue. Keeping the complete
        // centerline-shaped residual for every rejected fragment makes the
        // breadth-first search recursively chase thousands of sub-bead lanes.
        actual_uncovered     = opening_ex(actual_uncovered, float(0.15 * spacing));
        ExPolygons remaining = diff_ex(ExPolygons{expolygon}, primary_covered);
        append(remaining, std::move(actual_uncovered));
        remaining            = union_ex(remaining);
        auto prune_remaining = [spacing](ExPolygons& regions) {
            // A line-width-pitch family is continuous deposited material even
            // though polygonizing each centerline independently leaves
            // quantized sub-bead slivers. Do not recursively fill those
            // numerical seams; retain only residuals wide enough to represent
            // an actual hole in the bottom surface.
            regions.erase(std::remove_if(regions.begin(), regions.end(),
                                         [spacing](const ExPolygon& part) {
                                             if (std::abs(area(part)) < 0.02 * double(spacing) * double(spacing))
                                                 return true;
                                             // Use opening only as a width
                                             // test. Replacing the residual by
                                             // its opened shape removes narrow
                                             // necks which may lead to a large
                                             // real uncovered lobe.
                                             return opening_ex(ExPolygons{part}, float(0.15 * spacing)).empty();
                                         }),
                          regions.end());
            std::sort(regions.begin(), regions.end(),
                      [](const ExPolygon& lhs, const ExPolygon& rhs) { return std::abs(area(lhs)) > std::abs(area(rhs)); });
        };
        prune_remaining(remaining);

        constexpr size_t max_child_paths     = 8192;
        constexpr size_t max_parents         = 2048;
        constexpr double max_child_radius_mm = 15.;
        size_t child_path_count              = 0;

        struct RecursiveParent
        {
            Polylines frontier;
            Point center;
        };
        std::vector<RecursiveParent> parents;
        parents.reserve(256);
        // The reference implementation queues the completed outer ring of a
        // family. Native collision and anchor validation can also leave a
        // substantial hole inside that nominal disk, however. Queue any
        // accepted primary component that actually borders such a residual
        // region; queuing every inner ring unconditionally repeatedly
        // regenerated already covered space and made slicing quadratic.
        const ExPolygons parent_seed_regions = offset_ex(remaining, float(1.10 * line_width));
        for (const Polyline& path : arcs) {
            const bool borders_remaining = std::any_of(path.points.begin(), path.points.end(), [&parent_seed_regions](const Point& point) {
                return std::any_of(parent_seed_regions.begin(), parent_seed_regions.end(),
                                   [&point](const ExPolygon& region) { return region.contains(point); });
            });
            if (borders_remaining)
                parents.push_back({Polylines{path}, center});
        }
        if (parents.empty()) {
            for (const Polyline& path : select_outer_frontier(arcs, center))
                parents.push_back({Polylines{path}, center});
        }
        if (params.arc_obstacle_paths != nullptr) {
            for (const Polyline& path : *params.arc_obstacle_paths)
                parents.push_back({Polylines{path}, path.first_point()});
        }
        if (params.arc_prior_paths != nullptr) {
            for (const Polyline& path : *params.arc_prior_paths)
                parents.push_back({Polylines{path}, path.first_point()});
        }
        Lines surface_boundary_lines = to_lines(expolygon);
        const AABBTreeLines::LinesDistancer<Line> surface_boundary_index(std::move(surface_boundary_lines));

        size_t parent_idx = 0;
        while (parent_idx < parents.size() && parent_idx < max_parents && !remaining.empty() && child_path_count < max_child_paths) {
            Lines parent_lines = to_lines(parents[parent_idx].frontier);
            if (parent_lines.empty()) {
                ++parent_idx;
                continue;
            }
            const AABBTreeLines::LinesDistancer<Line> parent_index(std::move(parent_lines));
            bool spawned_child = false;

            for (size_t region_idx = 0; region_idx < remaining.size() && child_path_count < max_child_paths; ++region_idx) {
                // Port the reference get_farthest_point rule: select an
                // exposed point on the current parent which is still adjacent
                // to the unfilled region and farthest from the model boundary.
                Point child_anchor;
                Point child_region_point;
                double nearest_distance           = std::numeric_limits<double>::max();
                double farthest_boundary_distance = std::numeric_limits<double>::lowest();
                for (const Line& boundary : to_lines(remaining[region_idx])) {
                    const double boundary_length = std::max<double>(boundary.length(), 1.);
                    const size_t samples         = std::max<size_t>(1, size_t(std::ceil(boundary_length / std::max<double>(spacing, 1.))));
                    for (size_t sample_idx = 0; sample_idx <= samples; ++sample_idx) {
                        const double t     = double(sample_idx) / double(samples);
                        const Vec2d sample = boundary.a.cast<double>() + t * (boundary.b - boundary.a).cast<double>();
                        const Point point(coord_t(std::lround(sample.x())), coord_t(std::lround(sample.y())));
                        const auto [distance, line_idx, nearest_point] = parent_index.distance_from_lines_extra<false>(point);
                        if (distance > 1.50 * line_width)
                            continue;
                        const Point parent_point(coord_t(std::lround(nearest_point.x())), coord_t(std::lround(nearest_point.y())));
                        const auto [boundary_distance, boundary_line_idx,
                                    nearest_boundary_point] = surface_boundary_index.distance_from_lines_extra<false>(parent_point);
                        if (boundary_distance > farthest_boundary_distance ||
                            (boundary_distance == farthest_boundary_distance && distance < nearest_distance)) {
                            farthest_boundary_distance = boundary_distance;
                            nearest_distance           = distance;
                            child_region_point         = point;
                            child_anchor               = parent_point;
                        }
                    }
                }
                if (nearest_distance > 1.50 * line_width)
                    continue;

                Vec2d remaining_direction = (child_region_point - child_anchor).cast<double>();
                if (remaining_direction.squaredNorm() == 0.)
                    remaining_direction = axis;
                else
                    remaining_direction.normalize();

                Polylines accepted_family;
                Point child_center           = child_anchor;
                double accepted_outer_radius = 0.;
                // The postprocessor uses a 2 mm center offset, then retries
                // tight regions with zero offset. Keep a bead-width midpoint
                // attempt because it gives small details a supported runway.
                const std::array<double, 3> center_offsets{scale_(2.), double(line_width), 0.};
                for (const double center_offset : center_offsets) {
                    child_center              = child_anchor - Point(coord_t(std::lround(center_offset * remaining_direction.x())),
                                                                     coord_t(std::lround(center_offset * remaining_direction.y())));
                    const double first_radius = center_offset + double(line_width) / 1.5;
                    const double radius_limit = std::max(first_radius, std::min(scale_(max_child_radius_mm),
                                                                                farthest_boundary_distance - 0.5 * double(spacing)));
                    Polylines child_arcs;
                    generate_family(child_anchor, child_center, first_radius, radius_limit, remaining[region_idx],
                                    max_child_paths - child_path_count, false, child_arcs, nullptr, nullptr);

                    // Process paths in radius order so every accepted ring can
                    // support the next. A clipped component with only one
                    // supported end is omitted, while independently printable
                    // siblings remain available to fill the layer.
                    size_t child_idx = 0;
                    while (child_idx < child_arcs.size()) {
                        const double group_radius = family_path_radius(child_arcs[child_idx], child_center);
                        size_t group_end          = child_idx + 1;
                        while (group_end < child_arcs.size() &&
                               std::abs(family_path_radius(child_arcs[group_end], child_center) - group_radius) < 0.25 * spacing)
                            ++group_end;

                        Polylines candidate_group;
                        bool crossing_group = false;
                        for (size_t path_idx = child_idx; path_idx < group_end; ++path_idx) {
                            Polylines fragments = split_at_obstacle_paths(std::move(child_arcs[path_idx]));
                            if (fragments.empty())
                                continue;
                            for (Polyline& fragment : fragments) {
                                Polyline trimmed = trim_at_printed_path(std::move(fragment), true);
                                if (trimmed.empty() ||
                                    is_short_free_air_path(trimmed) ||
                                    is_shallow_child_arc(trimmed))
                                    continue;
                                if (!has_stable_emitted_branch_lead(
                                        trimmed))
                                    trimmed.reverse();
                                if (!has_stable_emitted_branch_lead(
                                        trimmed))
                                    continue;
                                // Reversing an independently clipped child can
                                // turn its formerly harmless trailing contact
                                // into a duplicate lead along the parent.
                                // Re-run proximity trimming in the actual print
                                // direction before accepting the path.
                                trim_sustained_retrace(
                                    trimmed, printed_index);
                                if (trimmed.length() <
                                        std::max(
                                            0.25 * double(spacing),
                                            0.75 * double(line_width)))
                                    continue;
                                candidate_group.emplace_back(
                                    std::move(trimmed));
                            }
                        }
                        if (candidate_group.empty())
                            break;
                        // Components from one clipped circle are validated
                        // before any of them enters printed_index. An
                        // independently printable sibling must not be discarded
                        // merely because another clipped component has only one
                        // supported end. Reject the radius only if two of the
                        // retained siblings actually cross.
                        for (size_t lhs_idx = 0; lhs_idx < candidate_group.size() && !crossing_group; ++lhs_idx) {
                            for (size_t rhs_idx = lhs_idx + 1; rhs_idx < candidate_group.size() && !crossing_group; ++rhs_idx) {
                                const Polyline& lhs = candidate_group[lhs_idx];
                                const Polyline& rhs = candidate_group[rhs_idx];
                                for (size_t lhs_segment = 1; lhs_segment < lhs.points.size() && !crossing_group; ++lhs_segment) {
                                    for (size_t rhs_segment = 1; rhs_segment < rhs.points.size(); ++rhs_segment) {
                                        Point intersection;
                                        if (!segment_intersection(Line(lhs.points[lhs_segment - 1], lhs.points[lhs_segment]),
                                                                  Line(rhs.points[rhs_segment - 1], rhs.points[rhs_segment]), intersection))
                                            continue;
                                        const bool lhs_endpoint = intersection == lhs.first_point() || intersection == lhs.last_point();
                                        const bool rhs_endpoint = intersection == rhs.first_point() || intersection == rhs.last_point();
                                        if (!lhs_endpoint || !rhs_endpoint) {
                                            crossing_group = true;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                        if (crossing_group)
                            break;

                        // Siblings from one radius may physically support each
                        // other's trailing lead. Resolve that dependency in a
                        // printable order: only add a path after both of its
                        // ends are supported by previously deposited material.
                        // If no candidate can advance, the remaining cyclic or
                        // one-ended group is intentionally left for another
                        // family.
                        Polylines accepted_group;
                        while (!candidate_group.empty()) {
                            bool accepted_candidate = false;
                            for (size_t candidate_idx = 0; candidate_idx < candidate_group.size(); ++candidate_idx) {
                                if (!has_supported_ends(candidate_group[candidate_idx]))
                                    continue;
                                printed_index.add(candidate_group[candidate_idx]);
                                accepted_group.emplace_back(std::move(candidate_group[candidate_idx]));
                                candidate_group.erase(candidate_group.begin() + candidate_idx);
                                accepted_candidate = true;
                                break;
                            }
                            if (!accepted_candidate)
                                break;
                        }
                        if (accepted_group.empty())
                            break;
                        append(accepted_family, std::move(accepted_group));
                        accepted_outer_radius = group_radius;
                        child_idx             = group_end;
                    }
                    if (!accepted_family.empty())
                        break;
                }
                if (accepted_family.empty())
                    continue;

                Polylines accepted_frontier = select_outer_frontier(accepted_family, child_center);
                // Use the reference filled-disk model when the deposited beads
                // cover it apart from sub-bead pitch seams. If rejecting a
                // one-ended component leaves a substantial portion absent,
                // subtract only actual bead footprints so another supported
                // family can fill that hole.
                ExPolygons ideal_swept_area =
                    swept_disk(child_center, accepted_outer_radius,
                               remaining[region_idx]);
                ExPolygons deposited_area =
                    union_ex(offset(
                        accepted_family,
                        float(0.52 * line_width)));
                const ExPolygons substantial_missing =
                    opening_ex(
                        diff_ex(ideal_swept_area, deposited_area),
                        float(0.30 * spacing));
                ExPolygons accepted_swept_area =
                    substantial_missing.empty() ?
                        std::move(ideal_swept_area) :
                        std::move(deposited_area);

                append(arcs, accepted_family);
                child_path_count += accepted_family.size();
                remaining = diff_ex(remaining, accepted_swept_area);
                prune_remaining(remaining);
                if (!accepted_frontier.empty() && parents.size() < max_parents)
                    parents.push_back({std::move(accepted_frontier), child_center});
                spawned_child = true;
                // Re-evaluate this same parent against the newly exposed
                // remaining-space boundary. This is the reference BFS branch
                // loop and is essential when one parent arc feeds two pockets.
                break;
            }

            if (!spawned_child)
                ++parent_idx;
        }

    }
    // Families are generated and refined independently, and a complex bridge
    // surface may contain hundreds of clipped recursive fragments.  Validate
    // the final combined sequence once more before returning it.  This catches
    // intersections between separately accepted families (including families
    // inherited from an earlier fill surface) and avoids emitting tiny
    // nozzle-scale fragments which otherwise become isolated blobs in G-code.
    Lines final_lines;
    PrintedArcIndex final_index(final_lines, 4 * spacing);
    if (params.arc_prior_paths != nullptr) {
        for (const Polyline &path : *params.arc_prior_paths)
            final_index.add(path);
    }
    if (params.arc_obstacle_paths != nullptr) {
        for (const Polyline &path : *params.arc_obstacle_paths)
            final_index.add(path);
    }
    Polylines sanitized_arcs;
    sanitized_arcs.reserve(arcs.size());
    size_t sanitized_primary_count = 0;
    // Keep three-quarters of a bead of printable centerline. Short recursive
    // fragments at this scale still close visible nozzle-width pockets on
    // narrow bridges; the collision and sustained-retrace passes below reject
    // fragments which would merely redraw a perimeter or an earlier arc.
    const double minimum_arc_length =
        std::max(0.25 * double(spacing),
                 0.75 * double(line_width));
    const coord_t final_collision_tolerance =
        std::max<coord_t>(
            SCALED_EPSILON, gcode_coordinate_step);
    const bool has_final_support_context =
        (params.arc_root_anchor_regions != nullptr &&
         !params.arc_root_anchor_regions->empty()) ||
        (params.arc_anchor_regions != nullptr &&
         !params.arc_anchor_regions->empty()) ||
        (params.arc_prior_paths != nullptr &&
         !params.arc_prior_paths->empty()) ||
        (params.arc_obstacle_paths != nullptr &&
         !params.arc_obstacle_paths->empty());
    const auto is_fully_supported_by_final_sequence =
        [&params, &final_index, lateral_support_distance,
         support_sample_step](const Polyline &path) {
            const auto point_supported =
                [&params, &final_index,
                 lateral_support_distance](const Point &point) {
                    const auto in_regions =
                        [&point](const ExPolygons *regions) {
                            return regions != nullptr &&
                                   std::any_of(
                                       regions->begin(), regions->end(),
                                       [&point](const ExPolygon &region) {
                                           return region.contains(point);
                                       });
                        };
                    if (in_regions(params.arc_root_anchor_regions) ||
                        in_regions(params.arc_anchor_regions))
                        return true;

                    double nearest_distance_squared =
                        std::numeric_limits<double>::max();
                    final_index.visit(
                        Line(point, point), lateral_support_distance,
                        [&nearest_distance_squared,
                         &point](const Line &support) {
                            nearest_distance_squared = std::min(
                                nearest_distance_squared,
                                support.distance_to_squared(point));
                        });
                    return nearest_distance_squared <=
                           double(lateral_support_distance) *
                               double(lateral_support_distance);
                };

            for (size_t point_idx = 1;
                 point_idx < path.points.size(); ++point_idx) {
                const Point &start = path.points[point_idx - 1];
                const Point &end   = path.points[point_idx];
                const Vec2d delta  = (end - start).cast<double>();
                const double length = delta.norm();
                if (length <= 0.)
                    continue;
                const size_t samples = std::max<size_t>(
                    1, size_t(std::ceil(
                           length / double(support_sample_step))));
                for (size_t sample_idx = 0;
                     sample_idx <= samples; ++sample_idx) {
                    const Vec2d position =
                        start.cast<double>() +
                        delta * (double(sample_idx) /
                                 double(samples));
                    if (!point_supported(Point(
                            coord_t(std::lround(position.x())),
                            coord_t(std::lround(position.y())))))
                        return false;
                }
            }
            return true;
        };
    for (size_t arc_idx = 0; arc_idx < arcs.size(); ++arc_idx) {
        Polyline &arc = arcs[arc_idx];
        if (arc.points.size() < 2)
            continue;

        Polyline sanitized;
        sanitized.points.emplace_back(arc.first_point());
        bool collided = false;
        for (size_t point_idx = 1;
             point_idx < arc.points.size() && !collided; ++point_idx) {
            const Line segment(
                arc.points[point_idx - 1], arc.points[point_idx]);
            Point nearest_intersection;
            double nearest_distance_squared =
                std::numeric_limits<double>::max();
            final_index.visit(
                segment, final_collision_tolerance,
                [&](const Line &occupied) {
                    Point intersection;
                    // The rounded proper-intersection helper alone misses
                    // acute crossings after coordinates are quantized to
                    // G-code resolution. At this final pass every actual
                    // segment intersection must terminate the later path.
                    if (!segment_intersection(
                            segment, occupied, intersection))
                        return;
                    const double distance_squared =
                        (intersection - segment.a)
                            .cast<double>()
                            .squaredNorm();
                    // Contact at the departure point is the intended physical
                    // anchor. Any later contact terminates the extrusion
                    // before its centerline crosses deposited material.
                    if (point_idx == 1 &&
                        distance_squared <=
                            double(SCALED_EPSILON) *
                                double(SCALED_EPSILON))
                        return;
                    if (distance_squared <
                        nearest_distance_squared) {
                        nearest_distance_squared = distance_squared;
                        nearest_intersection = intersection;
                    }
                });

            if (nearest_distance_squared <
                std::numeric_limits<double>::max()) {
                Point safe_endpoint = nearest_intersection;
                const Vec2d approach =
                    (nearest_intersection -
                     sanitized.last_point())
                        .cast<double>();
                const double approach_length =
                    approach.norm();
                // The G-code writer rounds XY to 0.001 mm. Ending exactly at
                // an analytic intersection may consequently round a few
                // microns past the occupied centerline. Stop one output step
                // early; the deposited beads still overlap while their
                // centerlines cannot cross after serialization.
                if (approach_length >
                    double(gcode_coordinate_step)) {
                    const Vec2d position =
                        nearest_intersection.cast<double>() -
                        approach *
                            (double(gcode_coordinate_step) /
                             approach_length);
                    safe_endpoint = Point(
                        coord_t(std::lround(position.x())),
                        coord_t(std::lround(position.y())));
                }
                if (safe_endpoint !=
                    sanitized.last_point())
                    sanitized.points.emplace_back(
                        safe_endpoint);
                collided = true;
            } else {
                sanitized.points.emplace_back(segment.b);
            }
        }

        if (sanitized.length() < minimum_arc_length)
            continue;
        // The collision pass may shorten an otherwise curved arc until the
        // surviving free-air portion is effectively a straight bridge. Keep
        // short nozzle-scale gap fillers, but do not reintroduce the shallow
        // spans rejected during family generation.
        if (is_shallow_arc(sanitized) &&
            has_final_support_context &&
            !is_fully_supported_by_final_sequence(sanitized))
            continue;
        final_index.add(sanitized);
        sanitized_arcs.emplace_back(std::move(sanitized));
        if (arc_idx < primary_arc_count)
            ++sanitized_primary_count;
    }
    arcs = std::move(sanitized_arcs);

    // Clipping a concentric family around holes and retained walls can
    // interleave disconnected corner fragments. Defer a fragment until its
    // complete bead is beside already deposited material. This is a stable
    // retry pass: the large primary family remains ahead of recursive detail,
    // and only paths which originally appeared too early move later. Cyclic
    // components need one bridge seed; choose the seed with the shortest
    // unsupported run, then resume dependency ordering.
    if (arcs.size() > 1) {
        Lines ordering_lines;
        PrintedArcIndex ordering_index(ordering_lines, 4 * spacing);
        if (params.arc_prior_paths != nullptr) {
            for (const Polyline &path : *params.arc_prior_paths)
                ordering_index.add(path);
        }
        if (params.arc_obstacle_paths != nullptr) {
            for (const Polyline &path : *params.arc_obstacle_paths)
                ordering_index.add(path);
        }

        const coord_t ordering_support_distance = line_width;
        const double ordering_support_distance_squared =
            double(ordering_support_distance) *
            double(ordering_support_distance);
        const auto point_in_regions =
            [](const Point &point, const ExPolygons *regions) {
                return regions != nullptr &&
                       std::any_of(
                           regions->begin(), regions->end(),
                           [&point](const ExPolygon &region) {
                               return region.contains(point);
                           });
            };
        const auto point_supported =
            [&params, &ordering_index,
             &point_in_regions,
             ordering_support_distance,
             ordering_support_distance_squared](const Point &point) {
                if (point_in_regions(
                        point, params.arc_root_anchor_regions) ||
                    point_in_regions(
                        point, params.arc_anchor_regions))
                    return true;
                double nearest_distance_squared =
                    std::numeric_limits<double>::max();
                ordering_index.visit(
                    Line(point, point),
                    ordering_support_distance,
                    [&nearest_distance_squared,
                     &point](const Line &support) {
                        nearest_distance_squared =
                            std::min(
                                nearest_distance_squared,
                                support.distance_to_squared(point));
                    });
                return nearest_distance_squared <=
                       ordering_support_distance_squared;
            };
        const auto has_root_supported_origin =
            [&params, &point_in_regions, support_sample_step,
             supported_lead_length](const Polyline &path) {
                if (params.arc_root_anchor_regions == nullptr ||
                    params.arc_root_anchor_regions->empty() ||
                    path.points.size() < 2)
                    return false;
                coord_t remaining = supported_lead_length;
                for (size_t point_idx = 1;
                     point_idx < path.points.size() &&
                     remaining > 0; ++point_idx) {
                    const Point &start = path.points[point_idx - 1];
                    const Point &end   = path.points[point_idx];
                    const Vec2d delta  = (end - start).cast<double>();
                    const double length = delta.norm();
                    if (length <= 0.)
                        continue;
                    const double checked_length =
                        std::min<double>(length, remaining);
                    const size_t samples = std::max<size_t>(
                        1, size_t(std::ceil(
                               checked_length /
                               double(support_sample_step))));
                    for (size_t sample_idx = 0;
                         sample_idx <= samples; ++sample_idx) {
                        const Vec2d position =
                            start.cast<double>() +
                            delta *
                                (checked_length *
                                 double(sample_idx) /
                                 double(samples) / length);
                        if (!point_in_regions(
                                Point(
                                    coord_t(std::lround(position.x())),
                                    coord_t(std::lround(position.y()))),
                                params.arc_root_anchor_regions))
                            return false;
                    }
                    remaining -= coord_t(
                        std::lround(checked_length));
                }
                return true;
            };
        const auto longest_unsupported_run =
            [&point_supported,
             support_sample_step](const Polyline &path) {
                double longest = 0.;
                double current = 0.;
                for (size_t point_idx = 1;
                     point_idx < path.points.size(); ++point_idx) {
                    const Point &start = path.points[point_idx - 1];
                    const Point &end   = path.points[point_idx];
                    const Vec2d delta  = (end - start).cast<double>();
                    const double length = delta.norm();
                    if (length <= 0.)
                        continue;
                    const size_t samples = std::max<size_t>(
                        1, size_t(std::ceil(
                               length /
                               double(support_sample_step))));
                    for (size_t sample_idx = 0;
                         sample_idx <= samples; ++sample_idx) {
                        const Vec2d position =
                            start.cast<double>() +
                            delta *
                                (double(sample_idx) /
                                 double(samples));
                        if (point_supported(Point(
                                coord_t(std::lround(position.x())),
                                coord_t(std::lround(position.y()))))) {
                            current = 0.;
                        } else if (sample_idx > 0) {
                            current += length / double(samples);
                            longest = std::max(longest, current);
                        }
                    }
                }
                return longest;
            };
        const auto fully_supported =
            [&point_supported,
             support_sample_step](const Polyline &path) {
                for (size_t point_idx = 1;
                     point_idx < path.points.size(); ++point_idx) {
                    const Point &start = path.points[point_idx - 1];
                    const Point &end   = path.points[point_idx];
                    const Vec2d delta  = (end - start).cast<double>();
                    const double length = delta.norm();
                    if (length <= 0.)
                        continue;
                    const size_t samples = std::max<size_t>(
                        1, size_t(std::ceil(
                               length /
                               double(support_sample_step))));
                    for (size_t sample_idx = 0;
                         sample_idx <= samples; ++sample_idx) {
                        const Vec2d position =
                            start.cast<double>() +
                            delta *
                                (double(sample_idx) /
                                 double(samples));
                        if (!point_supported(Point(
                                coord_t(std::lround(position.x())),
                                coord_t(std::lround(position.y())))))
                            return false;
                    }
                }
                return true;
            };
        const auto order_group =
            [&ordering_index, &longest_unsupported_run,
             &fully_supported, &point_supported,
             &has_root_supported_origin]
            (Polylines pending, Polylines &ordered,
             bool seed_original_root) {
                if (seed_original_root && !pending.empty()) {
                    size_t root_idx =
                        std::numeric_limits<size_t>::max();
                    bool reverse_root = false;
                    for (size_t path_idx = 0;
                         path_idx < pending.size(); ++path_idx) {
                        if (has_root_supported_origin(
                                pending[path_idx])) {
                            root_idx = path_idx;
                            break;
                        }
                        Polyline reversed = pending[path_idx];
                        reversed.reverse();
                        if (has_root_supported_origin(reversed)) {
                            root_idx = path_idx;
                            reverse_root = true;
                            break;
                        }
                    }
                    if (root_idx !=
                        std::numeric_limits<size_t>::max()) {
                        if (reverse_root)
                            pending[root_idx].reverse();
                        ordering_index.add(pending[root_idx]);
                        ordered.emplace_back(
                            std::move(pending[root_idx]));
                        pending.erase(
                            pending.begin() + root_idx);
                    }
                }
                while (!pending.empty()) {
                    Polylines deferred;
                    deferred.reserve(pending.size());
                    bool progressed = false;
                    for (Polyline &path : pending) {
                        if (!fully_supported(path)) {
                            deferred.emplace_back(std::move(path));
                            continue;
                        }
                        if (!ordered.empty() &&
                            (path.last_point() -
                             ordered.back().last_point())
                                    .cast<double>()
                                    .squaredNorm() <
                                (path.first_point() -
                                 ordered.back().last_point())
                                    .cast<double>()
                                    .squaredNorm())
                            path.reverse();
                        ordering_index.add(path);
                        ordered.emplace_back(std::move(path));
                        progressed = true;
                    }
                    if (deferred.empty())
                        break;
                    if (!progressed) {
                        size_t seed_idx =
                            std::numeric_limits<size_t>::max();
                        bool reverse_seed = false;
                        double shortest_run =
                            std::numeric_limits<double>::max();
                        for (size_t path_idx = 0;
                             path_idx < deferred.size(); ++path_idx) {
                            const bool first_supported =
                                point_supported(
                                    deferred[path_idx].first_point());
                            if (!first_supported)
                                continue;
                            const bool last_supported =
                                point_supported(
                                    deferred[path_idx].last_point());
                            const double run =
                                longest_unsupported_run(
                                    deferred[path_idx]);
                            if (run < shortest_run) {
                                shortest_run = run;
                                seed_idx = path_idx;
                                reverse_seed =
                                    last_supported &&
                                    !ordered.empty() &&
                                    (deferred[path_idx].last_point() -
                                     ordered.back().last_point())
                                            .cast<double>()
                                            .squaredNorm() <
                                        (deferred[path_idx].first_point() -
                                         ordered.back().last_point())
                                            .cast<double>()
                                            .squaredNorm();
                            }
                        }
                        if (seed_idx ==
                                std::numeric_limits<size_t>::max() &&
                            !ordered.empty()) {
                            for (size_t path_idx = 0;
                                 path_idx < deferred.size(); ++path_idx) {
                                if (!point_supported(
                                        deferred[path_idx].last_point()) ||
                                    (deferred[path_idx].last_point() -
                                     ordered.back().last_point())
                                            .cast<double>()
                                            .squaredNorm() >
                                        (deferred[path_idx].first_point() -
                                         ordered.back().last_point())
                                            .cast<double>()
                                            .squaredNorm())
                                    continue;
                                const double run =
                                    longest_unsupported_run(
                                        deferred[path_idx]);
                                if (run < shortest_run) {
                                    shortest_run = run;
                                    seed_idx = path_idx;
                                    reverse_seed = true;
                                }
                            }
                        }
                        // Generation already requires an anchored origin, so
                        // this fallback is only for degenerate coordinate
                        // quantization where neither endpoint survives the
                        // final support-distance test.
                        if (seed_idx ==
                            std::numeric_limits<size_t>::max()) {
                            seed_idx = 0;
                            for (size_t path_idx = 0;
                                 path_idx < deferred.size(); ++path_idx) {
                                const double run =
                                    longest_unsupported_run(
                                        deferred[path_idx]);
                                if (run < shortest_run) {
                                    shortest_run = run;
                                    seed_idx = path_idx;
                                }
                            }
                            if (!ordered.empty() &&
                                (deferred[seed_idx].last_point() -
                                 ordered.back().last_point())
                                        .cast<double>()
                                        .squaredNorm() <
                                    (deferred[seed_idx].first_point() -
                                     ordered.back().last_point())
                                        .cast<double>()
                                        .squaredNorm())
                                reverse_seed = true;
                        }
                        if (reverse_seed)
                            deferred[seed_idx].reverse();
                        ordering_index.add(deferred[seed_idx]);
                        ordered.emplace_back(
                            std::move(deferred[seed_idx]));
                        deferred.erase(
                            deferred.begin() + seed_idx);
                    }
                    pending = std::move(deferred);
                }
            };

        Polylines primary;
        primary.reserve(sanitized_primary_count);
        for (size_t path_idx = 0;
             path_idx < sanitized_primary_count; ++path_idx)
            primary.emplace_back(std::move(arcs[path_idx]));
        Polylines recursive;
        recursive.reserve(arcs.size() - sanitized_primary_count);
        for (size_t path_idx = sanitized_primary_count;
             path_idx < arcs.size(); ++path_idx)
            recursive.emplace_back(std::move(arcs[path_idx]));

        Polylines ordered;
        ordered.reserve(arcs.size());
        order_group(std::move(primary), ordered, true);
        sanitized_primary_count = ordered.size();
        order_group(std::move(recursive), ordered, false);
        // Standalone fill callers do not provide the previous-layer support
        // map. Preserve their historical same-side traversal without
        // pretending that either endpoint is mechanically preferred. In a
        // real print the explicit support data above remains authoritative.
        if (params.arc_prior_paths == nullptr &&
            params.arc_obstacle_paths == nullptr &&
            params.arc_root_anchor_regions == nullptr &&
            params.arc_anchor_regions == nullptr) {
            for (size_t path_idx = 1;
                 path_idx < ordered.size(); ++path_idx) {
                if ((ordered[path_idx].last_point() -
                     ordered[path_idx - 1].last_point())
                            .cast<double>()
                            .squaredNorm() <
                    (ordered[path_idx].first_point() -
                     ordered[path_idx - 1].last_point())
                            .cast<double>()
                            .squaredNorm())
                    ordered[path_idx].reverse();
            }
        }
        arcs = std::move(ordered);
    }

    // Dependency ordering may reverse a fragment after its generation-time
    // validation. Validate only the leading contact again in the exact order
    // that will reach G-code. Do not seed this index with perimeter obstacles:
    // sustained attachment to a wall is an anchor, while retracing an earlier
    // arc is duplicate extrusion.
    Lines emitted_arc_lines;
    PrintedArcIndex emitted_arc_index(
        emitted_arc_lines, 4 * spacing);
    Polylines emitted_arcs;
    emitted_arcs.reserve(arcs.size());
    for (size_t arc_idx = 0; arc_idx < arcs.size(); ++arc_idx) {
        Polyline &arc = arcs[arc_idx];
        trim_sustained_retrace(arc, emitted_arc_index);
        if (arc.length() < minimum_arc_length)
            continue;
        // Proximity trimming happens in final print order and can shorten a
        // valid recursive circle into a long, nearly straight chord. Such a
        // remnant no longer provides the curved free-air support promised by
        // the pattern; leave its residual for another child family instead of
        // emitting a sag-prone bridge line.
        if (arc_idx >= sanitized_primary_count &&
            is_shallow_child_arc(arc))
            continue;
        emitted_arc_index.add(arc);
        emitted_arcs.emplace_back(std::move(arc));
    }
    arcs = std::move(emitted_arcs);

    if (params.arc_prior_paths != nullptr)
        append(*params.arc_prior_paths, arcs);
    append(polylines_out, std::move(arcs));
}

} // namespace Slic3r
