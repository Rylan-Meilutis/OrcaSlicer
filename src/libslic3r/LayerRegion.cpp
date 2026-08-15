#include "Layer.hpp"
#include "BridgeDetector.hpp"
#include "ClipperUtils.hpp"
#include "Geometry.hpp"
#include "PerimeterGenerator.hpp"
#include "Point.hpp"
#include "Print.hpp"
#include "Surface.hpp"
#include "BoundingBox.hpp"
#include "SVG.hpp"
#include "Algorithm/RegionExpansion.hpp"

#include <string>
#include <map>
#include <limits>

#include <boost/log/trivial.hpp>
#include <boost/algorithm/clamp.hpp>

namespace Slic3r {

namespace {

bool path_is_fully_buried(const ExtrusionPath &path,
                          const ExPolygons &upper_coverage)
{
    const Polygons footprint = path.polygons_covered_by_width(
        float(SCALED_EPSILON));
    return !footprint.empty() && !upper_coverage.empty() && diff_ex(
        footprint, upper_coverage, ApplySafetyOffset::No).empty();
}

bool path_is_stagger_wall(const ExtrusionPath &path, int inset_idx,
                          bool inner_only)
{
    const bool alternating_inner = inset_idx > 0 && inset_idx % 2 == 1 &&
        path.role() == erPerimeter;
    const bool outer = !inner_only && inset_idx == 0 &&
        (path.role() == erExternalPerimeter || path.role() == erPerimeter);
    return (alternating_inner || outer) && path.polyline.points.size() >= 2;
}

bool path_is_stagger_candidate(const ExtrusionPath &path, int inset_idx,
                               bool inner_only)
{
    return path_is_stagger_wall(path, inset_idx, inner_only) &&
           !path.z_contoured;
}

bool path_can_be_staggered(const ExtrusionPath &path, int inset_idx,
                           bool inner_only,
                           const ExPolygons &current_coverage,
                           const ExPolygons &upper_coverage)
{
    if (!path_is_stagger_candidate(path, inset_idx, inner_only))
        return false;

    if (!inner_only && inset_idx == 0) {
        // An external bead intentionally straddles the model contour, so the
        // inner-wall "whole bead is inside the next slice" test can never
        // accept it. Instead require the slice boundary to be unchanged in
        // the complete bead corridor. This admits vertical structural walls
        // while rejecting top edges and sloped walls where either slice gains
        // or loses material beside the path.
        const Polygons footprint = path.polygons_covered_by_width(
            float(SCALED_EPSILON));
        if (footprint.empty() || current_coverage.empty() ||
            upper_coverage.empty())
            return false;
        const ExPolygons changed = xor_ex(
            current_coverage, upper_coverage, ApplySafetyOffset::No);
        return changed.empty() || intersection_ex(
            footprint, changed, ApplySafetyOffset::No).empty();
    }

    // Raising an inner wall is safe only when its complete deposited bead is
    // buried by the next model layer. Checking the centerline alone lets the
    // outside half of a bead emerge on sloped walls and beside top surfaces,
    // which makes the brick pattern visible and changes the model envelope.
    return path_is_fully_buried(path, upper_coverage);
}

void stagger_path(ExtrusionPath &path, int inset_idx, coord_t z_offset)
{
    path.inset_idx = inset_idx;
    for (Point3 &point : path.polyline.points)
        point.z() += z_offset;
    path.z_contoured = true;
    path.staggered_perimeter = true;
}

void set_staggered_transition_height(ExtrusionPath &path, coord_t z_offset,
                                     bool entering_course)
{
    const double original_height = std::max(0.01, double(path.height));
    const double offset = unscale<double>(z_offset);
    const double transition_height = entering_course ?
        original_height + offset : std::max(0.01, original_height - offset);
    path.mm3_per_mm *= transition_height / original_height;
    path.height = float(transition_height);
    // The short exit bead returns the course to the nominal layer grid. It is
    // ordinary inner-wall geometry, not a raised brick bead; keeping the brick
    // preview tag here made the final exposed perimeter look bricked even
    // though this half-height course is the required realignment layer.
    path.staggered_perimeter = entering_course;
    path.staggered_transition = true;
    path.z_contoured = true;
}

ExtrusionPaths add_staggered_course_transitions(
    const ExtrusionPath &source, int inherited_inset_idx, coord_t z_offset,
    bool inner_only, const ExPolygons &lower_staggered_coverage)
{
    const int inset_idx = source.inset_idx >= 0 ? source.inset_idx : inherited_inset_idx;
    if (!path_is_stagger_wall(source, inset_idx, inner_only) ||
        source.nonplanar_surface ||
        lower_staggered_coverage.empty()) {
        ExtrusionPath result = source;
        result.inset_idx = inset_idx;
        // With no lower brick course, every newly raised path is the entrance
        // bead and must span the otherwise missing half layer below it.
        if (result.staggered_perimeter)
            set_staggered_transition_height(result, z_offset, true);
        return {std::move(result)};
    }

    const Polyline source_polyline = source.polyline.to_polyline();
    Polylines fragments = intersection_pl(
        Polylines{source_polyline}, lower_staggered_coverage);
    append(fragments, diff_pl(Polylines{source_polyline}, lower_staggered_coverage));
    if (fragments.size() <= 1) {
        ExtrusionPath result = source;
        result.inset_idx = inset_idx;
        const bool lower_course_present = !intersection_pl(
            Polylines{source_polyline}, lower_staggered_coverage).empty();
        if (result.staggered_perimeter && !lower_course_present)
            set_staggered_transition_height(result, z_offset, true);
        else if (!result.staggered_perimeter && lower_course_present)
            set_staggered_transition_height(result, z_offset, false);
        return {std::move(result)};
    }
    restore_source_path_order(source_polyline, fragments);

    ExtrusionPaths result;
    result.reserve(fragments.size());
    for (const Polyline &fragment : fragments) {
        if (fragment.points.size() < 2)
            continue;
        result.emplace_back(Polyline3(fragment), source);
        ExtrusionPath &path = result.back();
        path.inset_idx = inset_idx;
        const Point &sample = fragment.points[fragment.points.size() / 2];
        const bool lower_course_present = std::any_of(
            lower_staggered_coverage.begin(), lower_staggered_coverage.end(),
            [&sample](const ExPolygon &polygon) { return polygon.contains(sample); });
        if (path.staggered_perimeter && !lower_course_present)
            set_staggered_transition_height(path, z_offset, true);
        else if (!path.staggered_perimeter && lower_course_present)
            set_staggered_transition_height(path, z_offset, false);
    }
    return result.empty() ? ExtrusionPaths{source} : result;
}

void ramp_staggered_endpoint(ExtrusionPath &path, bool at_start,
                             coord_t z_offset)
{
    Points3 &points = path.polyline.points;
    if (!path.staggered_perimeter || points.size() < 2 || z_offset <= 0)
        return;
    if (std::none_of(points.begin(), points.end(), [z_offset](const Point3 &point) {
            return point.z() >= z_offset;
        }))
        return;

    // This path contains a deliberate Brick-course entry/exit ramp. Keep that
    // provenance even if a later course-boundary pass classifies the bead as
    // an ordinary (non-raised) wall for preview purposes. Otherwise G-code
    // mistakes the variable Z for a non-planar finishing surface and enables
    // toolhead-clearance lifts for every subsequent path on the layer.
    path.staggered_transition = true;

    if (!at_start)
        path.polyline.reverse();
    const Point3 start = points.front();
    const Point3 next = points[1];
    const double segment_length =
        (next.to_point() - start.to_point()).cast<double>().norm();
    if (segment_length > 1.) {
        const double ramp_length = std::min(
            segment_length * 0.9, double(scale_(std::max(0.1f, path.width))));
        const auto interpolate = [&start, &next, segment_length](double distance,
                                                                coord_t z) {
            const double ratio = distance / segment_length;
            return Point3(
                coord_t(std::llround(double(start.x()) +
                    (double(next.x()) - double(start.x())) * ratio)),
                coord_t(std::llround(double(start.y()) +
                    (double(next.y()) - double(start.y())) * ratio)), z);
        };
        const coord_t base_z = start.z() - z_offset;
        points.front().z() = base_z;
        // Explicit half-height point makes the transition visible to the
        // motion planner and prevents an abrupt Z-only jump at the boundary
        // between nominal and brick-offset fragments.
        points.insert(points.begin() + 1,
                      interpolate(0.5 * ramp_length, base_z + z_offset / 2));
        if (ramp_length < segment_length * 0.9)
            points.insert(points.begin() + 2,
                          interpolate(ramp_length, base_z + z_offset));
        path.polyline.fitting_result.clear();
    }
    if (!at_start)
        path.polyline.reverse();
}

void smooth_staggered_transitions(ExtrusionPaths &paths, coord_t z_offset,
                                  bool closed)
{
    if (paths.size() < 2)
        return;
    const size_t transition_count = closed ? paths.size() : paths.size() - 1;
    for (size_t idx = 0; idx < transition_count; ++idx) {
        ExtrusionPath &before = paths[idx];
        ExtrusionPath &after = paths[(idx + 1) % paths.size()];
        const auto is_raised = [](const ExtrusionPath &path) {
            return path.staggered_perimeter && std::any_of(
                path.polyline.points.begin(), path.polyline.points.end(),
                [](const Point3 &point) { return point.z() > 0; });
        };
        const bool before_raised = is_raised(before);
        const bool after_raised = is_raised(after);
        if (before_raised == after_raised ||
            before.last_point3().to_point() != after.first_point3().to_point())
            continue;
        if (before_raised)
            ramp_staggered_endpoint(before, false, z_offset);
        else
            ramp_staggered_endpoint(after, true, z_offset);
    }
}

ExtrusionPaths stagger_supported_fragments(const ExtrusionPath &source,
                                           int inherited_inset_idx,
                                           coord_t z_offset,
                                           bool inner_only,
                                           const ExPolygons &current_coverage,
                                           const ExPolygons &upper_coverage)
{
    const int inset_idx = source.inset_idx >= 0 ? source.inset_idx : inherited_inset_idx;
    if (!path_is_stagger_candidate(source, inset_idx, inner_only))
        return {source};
    if (path_can_be_staggered(
            source, inset_idx, inner_only, current_coverage, upper_coverage)) {
        ExtrusionPath result = source;
        stagger_path(result, inset_idx, z_offset);
        return {std::move(result)};
    }

    // Clip the centerline against coverage inset by the bead radius. This
    // partitions a path that crosses a top-surface boundary without allowing
    // any raised bead fragment to protrude beyond the next layer. Intersected
    // and subtracted fragments are restored to source order before extrusion.
    const coord_t clearance = scale_(0.5 * source.width) + coord_t(SCALED_EPSILON);
    const ExPolygons safe_centerline_area = shrink_ex(upper_coverage, clearance);
    if (safe_centerline_area.empty())
        return {source};

    const Polyline source_polyline = source.polyline.to_polyline();
    Polylines fragments = intersection_pl(Polylines{source_polyline}, safe_centerline_area);
    append(fragments, diff_pl(Polylines{source_polyline}, safe_centerline_area));
    if (fragments.size() <= 1)
        return {source};
    restore_source_path_order(source_polyline, fragments);

    ExtrusionPaths result;
    result.reserve(fragments.size());
    for (const Polyline &fragment : fragments) {
        if (fragment.points.size() < 2)
            continue;
        result.emplace_back(Polyline3(fragment), source);
        ExtrusionPath &path = result.back();
        path.inset_idx = inset_idx;
        if (path_can_be_staggered(
                path, inset_idx, inner_only, current_coverage, upper_coverage))
            stagger_path(path, inset_idx, z_offset);
    }
    return result.empty() ? ExtrusionPaths{source} : result;
}

coord_t entity_max_relative_z(const ExtrusionEntity &entity)
{
    const auto path_max_z = [](const ExtrusionPath &path) {
        coord_t maximum = 0;
        for (const Point3 &point : path.polyline.points)
            maximum = std::max(maximum, point.z());
        return maximum;
    };
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity))
        return path_max_z(*path);
    if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
        coord_t maximum = 0;
        for (const ExtrusionPath &path : multipath->paths)
            maximum = std::max(maximum, path_max_z(path));
        return maximum;
    }
    if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
        coord_t maximum = 0;
        for (const ExtrusionPath &path : loop->paths)
            maximum = std::max(maximum, path_max_z(path));
        return maximum;
    }
    if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        coord_t maximum = 0;
        for (const ExtrusionEntity *child : collection->entities)
            maximum = std::max(maximum, entity_max_relative_z(*child));
        return maximum;
    }
    return 0;
}

void start_closed_course_at_lowest_path(ExtrusionPaths &paths)
{
    if (paths.size() < 2)
        return;
    const auto minimum_z = [](const ExtrusionPath &path) {
        coord_t minimum = std::numeric_limits<coord_t>::max();
        for (const Point3 &point : path.polyline.points)
            minimum = std::min(minimum, point.z());
        return minimum;
    };
    const auto first = std::min_element(paths.begin(), paths.end(),
        [&minimum_z](const ExtrusionPath &left, const ExtrusionPath &right) {
            const coord_t left_z = minimum_z(left);
            const coord_t right_z = minimum_z(right);
            if (left_z != right_z)
                return left_z < right_z;
            // At the same Z, start in the conventional wall. This keeps the
            // seam and prime on supported geometry before entering a raised
            // interlocking segment.
            return !left.staggered_perimeter && right.staggered_perimeter;
        });
    std::rotate(paths.begin(), first, paths.end());
}

void interlock_inner_wall_path(ExtrusionPath &path, int inset_idx,
                               double depth, double wavelength,
                               double resolution, int layer_id)
{
    if (path.polyline.points.size() < 2 || depth <= EPSILON ||
        wavelength <= EPSILON || resolution <= EPSILON)
        return;

    double total_length = 0.;
    for (const Line3 &line : path.polyline.lines())
        total_length += (line.b.to_point() - line.a.to_point()).cast<double>().norm() *
                        SCALING_FACTOR;
    if (total_length <= 2. * resolution)
        return;

    const double phase = (layer_id & 1) ? M_PI : 0.;
    const double transition_length = std::min(0.5 * wavelength, 0.25 * total_length);
    Points3 sampled;
    sampled.reserve(size_t(std::ceil(total_length / resolution)) + 2);
    sampled.emplace_back(path.polyline.points.front().x(),
                         path.polyline.points.front().y(), coord_t(0));

    double distance = 0.;
    for (const Line3 &line : path.polyline.lines()) {
        const Vec2d delta = (line.b.to_point() - line.a.to_point()).cast<double>();
        const double length = delta.norm() * SCALING_FACTOR;
        const size_t segments = std::max<size_t>(1, size_t(std::ceil(length / resolution)));
        for (size_t segment = 1; segment <= segments; ++segment) {
            const double t = double(segment) / double(segments);
            const Vec2d xy = line.a.to_point().cast<double>() + t * delta;
            const double along = distance + t * length;
            const double edge_distance = std::min(along, total_length - along);
            const double u = transition_length <= EPSILON ? 1. :
                std::clamp(edge_distance / transition_length, 0., 1.);
            const double envelope = u * u * (3. - 2. * u);
            // Keep the entire course at or below its nominal plane. This
            // creates alternating interlock without entering the next layer's
            // motion envelope or requiring a cross-layer dependency reorder.
            const double wave = 0.5 + 0.5 * std::sin(
                2. * M_PI * along / wavelength + phase);
            sampled.emplace_back(coord_t(std::llround(xy.x())),
                                 coord_t(std::llround(xy.y())),
                                 coord_t(-scale_(envelope * depth * wave)));
        }
        distance += length;
    }
    sampled.back().z() = 0;
    path.polyline.points = std::move(sampled);
    path.polyline.remove_duplicate_points();
    path.polyline.fitting_result.clear();
    path.inset_idx = inset_idx;
    path.z_contoured = true;
    path.nonplanar_interlocking_wall = true;
}

void apply_nonplanar_interlocking_walls(ExtrusionEntity &entity,
                                        int inherited_inset_idx,
                                        double depth, double wavelength,
                                        double resolution, int layer_id,
                                        const ExPolygons &current_coverage,
                                        const ExPolygons &upper_coverage)
{
    const int inset_idx = entity.inset_idx >= 0 ? entity.inset_idx : inherited_inset_idx;
    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
        if (path_can_be_staggered(*path, inset_idx, true,
                                  current_coverage, upper_coverage))
            interlock_inner_wall_path(*path, inset_idx, depth, wavelength,
                                      resolution, layer_id);
    } else if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity)) {
        for (ExtrusionPath &path : multipath->paths)
            apply_nonplanar_interlocking_walls(path, inset_idx, depth, wavelength,
                                               resolution, layer_id, current_coverage,
                                               upper_coverage);
    } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
        for (ExtrusionPath &path : loop->paths)
            apply_nonplanar_interlocking_walls(path, inset_idx, depth, wavelength,
                                               resolution, layer_id, current_coverage,
                                               upper_coverage);
    } else if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
        for (ExtrusionEntity *child : collection->entities)
            apply_nonplanar_interlocking_walls(*child, inset_idx, depth, wavelength,
                                               resolution, layer_id, current_coverage,
                                               upper_coverage);
    }
}

void apply_staggered_perimeters(ExtrusionEntity &entity, int inherited_inset_idx,
                                coord_t z_offset, bool inner_only,
                                const ExPolygons &current_coverage,
                                const ExPolygons &upper_coverage,
                                const ExPolygons &lower_staggered_coverage)
{
    const int inset_idx = entity.inset_idx >= 0 ? entity.inset_idx : inherited_inset_idx;
    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
        if (path_can_be_staggered(
                *path, inset_idx, inner_only, current_coverage, upper_coverage)) {
            stagger_path(*path, inset_idx, z_offset);
            const bool continues_lower_course = !intersection_pl(
                Polylines{path->polyline.to_polyline()}, lower_staggered_coverage).empty();
            if (!continues_lower_course)
                set_staggered_transition_height(*path, z_offset, true);
        } else if (path_is_stagger_candidate(*path, inset_idx, inner_only) &&
                   !intersection_pl(Polylines{path->polyline.to_polyline()},
                                    lower_staggered_coverage).empty()) {
            path->inset_idx = inset_idx;
            set_staggered_transition_height(*path, z_offset, false);
        }
    } else if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity)) {
        // Arachne may split one wall into ordinary and overhang paths. Do not
        // let one exposed segment disable brick staggering for every buried
        // segment in the same multipath. G-code already supports a different
        // Z on consecutive Polyline3 paths and emits the connecting Z move.
        ExtrusionPaths paths;
        for (const ExtrusionPath &path : multipath->paths) {
            ExtrusionPaths supported = stagger_supported_fragments(
                path, inset_idx, z_offset, inner_only, current_coverage,
                upper_coverage);
            for (const ExtrusionPath &fragment : supported)
                append(paths, add_staggered_course_transitions(
                    fragment, inset_idx, z_offset, inner_only,
                    lower_staggered_coverage));
        }
        smooth_staggered_transitions(paths, z_offset, false);
        multipath->paths = std::move(paths);
    } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
        // The same loop may contain independently classified overhang and
        // perimeter sections. Raise only the buried ordinary sections; the
        // exposed portions stay at nominal Z and retain the model envelope.
        ExtrusionPaths paths;
        for (const ExtrusionPath &path : loop->paths) {
            ExtrusionPaths supported = stagger_supported_fragments(
                path, inset_idx, z_offset, inner_only, current_coverage,
                upper_coverage);
            for (const ExtrusionPath &fragment : supported)
                append(paths, add_staggered_course_transitions(
                    fragment, inset_idx, z_offset, inner_only,
                    lower_staggered_coverage));
        }
        // A closed course has no geometric beginning, so choose its lowest
        // supported fragment as the emission start. This retains one
        // continuous loop while preventing a seam prime at the raised height
        // followed by a descent through an already deposited brick bead.
        start_closed_course_at_lowest_path(paths);
        smooth_staggered_transitions(paths, z_offset, true);
        loop->paths = std::move(paths);
    } else if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
        for (ExtrusionEntity *child : collection->entities)
            apply_staggered_perimeters(
                *child, inset_idx, z_offset, inner_only, current_coverage,
                upper_coverage,
                lower_staggered_coverage);
        // Resolve the physical Z dependency before applying the configured
        // wall order: lower courses must exist before a neighboring raised
        // wall is deposited. Stable ordering preserves the user's inner /
        // outer wall order among entities at the same relative height.
        std::stable_sort(collection->entities.begin(), collection->entities.end(),
            [](const ExtrusionEntity *left, const ExtrusionEntity *right) {
                return entity_max_relative_z(*left) < entity_max_relative_z(*right);
            });
    }
}

coord_t entity_min_relative_z(const ExtrusionEntity &entity)
{
    const auto path_min_z = [](const ExtrusionPath &path) {
        coord_t minimum = 0;
        for (const Point3 &point : path.polyline.points)
            minimum = std::min(minimum, point.z());
        return minimum;
    };
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity))
        return path_min_z(*path);
    coord_t minimum = 0;
    if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
        for (const ExtrusionPath &path : multipath->paths)
            minimum = std::min(minimum, path_min_z(path));
    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
        for (const ExtrusionPath &path : loop->paths)
            minimum = std::min(minimum, path_min_z(path));
    } else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        for (const ExtrusionEntity *child : collection->entities)
            minimum = std::min(minimum, entity_min_relative_z(*child));
    }
    return minimum;
}

bool smooth_outer_loop_is_eligible(const ExtrusionLoop &loop,
                                   const ExPolygons &changed_footprint)
{
    if (loop.paths.empty())
        return false;
    Polylines centerlines;
    for (const ExtrusionPath &path : loop.paths) {
        if (path.role() != erExternalPerimeter || path.z_contoured ||
            path.polyline.points.size() < 2)
            return false;
        centerlines.emplace_back(path.polyline.to_polyline());
    }
    if (changed_footprint.empty())
        return true;
    const coord_t clearance = scale_(0.5 * loop.paths.front().width) +
        coord_t(SCALED_EPSILON);
    return intersection_pl(centerlines, offset_ex(changed_footprint, clearance)).empty();
}

ExtrusionEntitiesPtr smooth_outer_wall_courses(
    const ExtrusionLoop &source, double layer_height, double target_height)
{
    // Choose the evenly dividing course height closest to the request. Merely
    // rounding layer_height / target_height selects one pass for e.g.
    // 0.20 / 0.14, even though two 0.10 mm courses are the closer result.
    const double ratio = layer_height / target_height;
    const int lower_passes = std::max(1, int(std::floor(ratio)));
    const int upper_passes = std::max(1, int(std::ceil(ratio)));
    const int passes =
        std::abs(layer_height / lower_passes - target_height) <=
                std::abs(layer_height / upper_passes - target_height) ?
            lower_passes : upper_passes;
    if (passes < 2)
        return {};
    const double course_height = layer_height / double(passes);
    ExtrusionEntitiesPtr courses;
    courses.reserve(size_t(passes));
    for (int pass = 0; pass < passes; ++pass) {
        auto *loop = static_cast<ExtrusionLoop *>(source.clone());
        const coord_t relative_z = scale_(
            -layer_height + (double(pass) + 1.) * course_height);
        for (ExtrusionPath &path : loop->paths) {
            for (Point3 &point : path.polyline.points)
                point.z() += relative_z;
            path.height = float(course_height);
            path.mm3_per_mm /= double(passes);
            path.z_contoured = true;
            path.nonplanar_clearance_validated = true;
            path.smooth_outer_wall = true;
            path.set_reverse();
        }
        courses.emplace_back(loop);
    }
    return courses;
}

void apply_smooth_outer_wall_courses(ExtrusionEntityCollection &collection,
                                     double layer_height, double target_height,
                                     const ExPolygons &changed_footprint)
{
    for (size_t idx = 0; idx < collection.entities.size(); ++idx) {
        ExtrusionEntity *entity = collection.entities[idx];
        if (auto *nested = dynamic_cast<ExtrusionEntityCollection *>(entity)) {
            apply_smooth_outer_wall_courses(
                *nested, layer_height, target_height, changed_footprint);
        } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(entity);
                   loop != nullptr &&
                   smooth_outer_loop_is_eligible(*loop, changed_footprint)) {
            ExtrusionEntitiesPtr courses = smooth_outer_wall_courses(
                *loop, layer_height, target_height);
            if (!courses.empty()) {
                delete entity;
                collection.entities.erase(collection.entities.begin() + idx);
                collection.entities.insert(collection.entities.begin() + idx,
                    courses.begin(), courses.end());
                idx += courses.size() - 1;
            }
        }
    }
    // Independent sub-courses must be deposited before any neighboring path
    // at the owning layer's nominal Z. This dependency takes precedence over
    // the selected wall sequence; stable ordering preserves that sequence for
    // entities sharing the same minimum Z.
    std::stable_sort(collection.entities.begin(), collection.entities.end(),
        [](const ExtrusionEntity *left, const ExtrusionEntity *right) {
            return entity_min_relative_z(*left) < entity_min_relative_z(*right);
        });
}

} // namespace

Flow LayerRegion::flow(FlowRole role) const
{
    return this->flow(role, m_layer->height);
}

Flow LayerRegion::flow(FlowRole role, double layer_height) const
{
    return m_region->flow(*m_layer->object(), role, layer_height, m_layer->id() == 0);
}

Flow LayerRegion::bridging_flow(FlowRole role, bool thick_bridge) const
{
    const PrintRegion       &region         = this->region();
    const PrintRegionConfig &region_config  = region.config();
    const PrintObject       &print_object   = *this->layer()->object();
    Flow bridge_flow;
    // Here this->extruder(role) - 1 may underflow to MAX_INT, but then the get_at() will fall back to zero'th element, so everything is all right.
    auto nozzle_diameter = float(print_object.print()->config().nozzle_diameter.get_at(region.extruder(role) - 1));
    const ConfigOptionFloatOrPercent& bridge_width_opt = region_config.bridge_line_width;
    const double                      bridge_width      = bridge_width_opt.get_abs_value(nozzle_diameter);
    const bool                        has_bridge_width  = bridge_width > 0.;
    const double                      bridge_flow_ratio = region_config.bridge_flow;

    if (thick_bridge) {
        // The old Slic3r way (different from all other slicers): Use rounded extrusions.
        // Get the configured nozzle_diameter for the extruder associated to the flow role requested.
        float thread_diameter = has_bridge_width ? float(bridge_width) : nozzle_diameter;
        if (bridge_flow_ratio > 0.)
            thread_diameter *= float(sqrt(bridge_flow_ratio));
        bridge_flow = Flow::bridging_flow(thread_diameter, nozzle_diameter);
    } else {
        // The same way as other slicers: Use normal extrusions. Apply bridge_flow while maintaining the original spacing.
        Flow base_flow = this->flow(role);
        if (has_bridge_width)
            base_flow = Flow(float(bridge_width), base_flow.height(), nozzle_diameter);
        bridge_flow = base_flow.with_flow_ratio(bridge_flow_ratio);
    }
    return bridge_flow;

}

// Fill in layerm->fill_surfaces by trimming the layerm->slices by the cummulative layerm->fill_surfaces.
void LayerRegion::slices_to_fill_surfaces_clipped()
{
    // Note: this method should be idempotent, but fill_surfaces gets modified 
    // in place. However we're now only using its boundaries (which are invariant)
    // so we're safe. This guarantees idempotence of prepare_infill() also in case
    // that combine_infill() turns some fill_surface into VOID surfaces.
    // Collect polygons per surface type.
    std::array<SurfacesPtr, size_t(stCount)> by_surface;
    for (Surface &surface : this->slices.surfaces)
        by_surface[size_t(surface.surface_type)].emplace_back(&surface);
    // Trim surfaces by the fill_boundaries.
    this->fill_surfaces.surfaces.clear();
    for (size_t surface_type = 0; surface_type < size_t(stCount); ++ surface_type) {
        const SurfacesPtr &this_surfaces = by_surface[surface_type];
        if (! this_surfaces.empty())
            this->fill_surfaces.append(intersection_ex(this_surfaces, this->fill_expolygons), SurfaceType(surface_type));
    }
}

void LayerRegion::make_perimeters(const SurfaceCollection &slices, const LayerRegionPtrs &compatible_regions, SurfaceCollection* fill_surfaces, ExPolygons* fill_no_overlap)
{
    this->perimeters.clear();
    this->thin_fills.clear();

    const PrintConfig       &print_config  = this->layer()->object()->print()->config();
    const PrintRegionConfig &region_config = this->region().config();
    const PrintObjectConfig& object_config = this->layer()->object()->config();
    // This needs to be in sync with PrintObject::_slice() slicing_mode_normal_below_layer!
    bool spiral_mode = print_config.spiral_mode &&
        //FIXME account for raft layers.
        (this->layer()->id() >= size_t(region_config.bottom_shell_layers.value) &&
         this->layer()->print_z >= region_config.bottom_shell_thickness - EPSILON);

    double model_rotation_rad = 0.0;
    if (region_config.align_infill_direction_to_model) {
        auto m = this->layer()->object()->trafo().matrix();
        model_rotation_rad = std::atan2((double)m(1, 0), (double)m(0, 0));
    }

    PerimeterGenerator g(
        // input:
        &slices,
        &compatible_regions,
        this->layer()->height,
        this->layer()->slice_z,
        this->flow(frPerimeter),
        &region_config,
        &this->layer()->object()->config(),
        &print_config,
        spiral_mode,
        model_rotation_rad,
        
        // output:
        &this->perimeters,
        &this->thin_fills,
        fill_surfaces,
        //BBS
        fill_no_overlap
    );
    
    if (this->layer()->lower_layer != nullptr)
        // Cummulative sum of polygons over all the regions.
        g.lower_slices = &this->layer()->lower_layer->lslices;
    if (this->layer()->upper_layer != NULL)
        g.upper_slices = &this->layer()->upper_layer->lslices;

    int region_id = this->region().print_object_region_id();
    if (this->layer()->upper_layer != NULL)
        g.upper_slices_same_region = &this->layer()->upper_layer->get_region(region_id)->slices;

    g.layer_id              = (int)this->layer()->id();
    g.ext_perimeter_flow    = this->flow(frExternalPerimeter);
    // Arc-overhang bottom surfaces replace unsupported wall segments. Keep any
    // remaining bridge perimeter anchors at ordinary perimeter flow: applying
    // thick/overextruded bridge flow to them makes the free edge sag before the
    // arc family can bond to it.
    g.overhang_flow         = region_config.arc_overhang_enabled ?
                                  this->flow(frPerimeter) :
                                  this->bridging_flow(frPerimeter, object_config.thick_bridges);
    g.solid_infill_flow     = this->flow(frSolidInfill);

    if (this->layer()->object()->config().wall_generator.value == PerimeterGeneratorType::Arachne && !spiral_mode)
        g.process_arachne();
    else
        g.process_classic();

    if (brick_perimeters_enabled(region_config) && this->layer()->lower_layer != nullptr &&
        !spiral_mode) {
        const double offset_ratio = region_config.staggered_perimeter_offset.get_abs_value(1.);
        const coord_t z_offset = scale_(this->layer()->height * offset_ratio);
        if (z_offset > 0) {
            // Aggregate object coverage is intentional: a wall is buried even
            // when the next layer covering it belongs to another print region.
            const ExPolygons upper_coverage = this->layer()->upper_layer == nullptr ?
                ExPolygons{} : this->layer()->upper_layer->lslices;
            ExPolygons lower_staggered_coverage;
            // Perimeters are generated in parallel across layers, so course
            // continuity must be derived from slice geometry rather than from
            // whether the lower layer's path entities happen to be ready. A
            // previous brick course can exist only where the same buried wall
            // corridor exists in both the current and lower slices. Layer zero
            // is deliberately conventional, making layer one an entry course.
            if (this->layer()->lower_layer->lower_layer != nullptr) {
                const coord_t clearance = scale_(0.5 * this->flow(frPerimeter).width()) +
                    coord_t(SCALED_EPSILON);
                lower_staggered_coverage = intersection_ex(
                    shrink_ex(this->layer()->lslices, clearance),
                    shrink_ex(this->layer()->lower_layer->lslices, clearance));
            }
            apply_staggered_perimeters(
                this->perimeters, -1, z_offset,
                region_config.staggered_perimeters_inner_only.value,
                this->layer()->lslices,
                upper_coverage,
                lower_staggered_coverage);
        }
    }
    if (region_config.perimeter_layering.value ==
            PerimeterLayeringMode::SmoothOuterWall &&
        this->layer()->lower_layer != nullptr && !spiral_mode) {
        const double target_height = region_config.smooth_outer_wall_layer_height.value;
        if (target_height > EPSILON && this->layer()->height > target_height + EPSILON) {
            const ExPolygons changed_footprint = union_ex(
                diff_ex(this->layer()->lslices,
                        this->layer()->lower_layer->lslices,
                        ApplySafetyOffset::No),
                diff_ex(this->layer()->lower_layer->lslices,
                        this->layer()->lslices,
                        ApplySafetyOffset::No));
            apply_smooth_outer_wall_courses(
                this->perimeters, this->layer()->height, target_height,
                changed_footprint);
        }
    }
    if ((region_config.perimeter_layering.value ==
             PerimeterLayeringMode::InterlockingWalls ||
         region_config.perimeter_layering.value ==
             PerimeterLayeringMode::Nonplanar) &&
        this->layer()->lower_layer != nullptr &&
        this->layer()->upper_layer != nullptr && !spiral_mode) {
        const double depth = std::min(
            region_config.interlocking_wall_amplitude.value,
            0.45 * this->layer()->height);
        apply_nonplanar_interlocking_walls(
            this->perimeters, -1, depth,
            region_config.interlocking_wall_wavelength.value,
            region_config.interlocking_wall_resolution.value,
            int(this->layer()->id()), this->layer()->lslices,
            this->layer()->upper_layer->lslices);
    }
}

#if 1

// Extract surfaces of given type from surfaces, extract fill (layer) thickness of one of the surfaces.
static ExPolygons fill_surfaces_extract_expolygons(Surfaces &surfaces, std::initializer_list<SurfaceType> surface_types, double &thickness)
{
    size_t cnt = 0;
    for (const Surface &surface : surfaces)
        if (std::find(surface_types.begin(), surface_types.end(), surface.surface_type) != surface_types.end()) {
            ++cnt;
            thickness = surface.thickness;
        }
    if (cnt == 0)
        return {};

    ExPolygons out;
    out.reserve(cnt);
    for (Surface &surface : surfaces)
        if (std::find(surface_types.begin(), surface_types.end(), surface.surface_type) != surface_types.end())
            out.emplace_back(std::move(surface.expolygon));
    return out;
}

struct ExpansionZone
{
    ExPolygons                           expolygons;
    Algorithm::RegionExpansionParameters parameters;
    bool                                 expanded_into = false;
};

// Cache for detecting bridge orientation and merging regions with overlapping expansions.
struct Bridge {
    ExPolygon expolygon;
    uint32_t group_id;
    std::vector<Algorithm::RegionExpansionEx>::const_iterator bridge_expansion_begin;
    std::optional<double> angle{std::nullopt};
};

// Group the bridge surfaces by overlaps.
uint32_t group_id(std::vector<Bridge> &bridges, uint32_t src_id) {
    uint32_t group_id = bridges[src_id].group_id;
    while (group_id != src_id) {
        src_id = group_id;
        group_id = bridges[src_id].group_id;
    }
    bridges[src_id].group_id = group_id;
    return group_id;
};

std::vector<Bridge> get_grouped_bridges(
    ExPolygons&& bridge_expolygons,
    const std::vector<Algorithm::RegionExpansionEx>& bridge_expansions
) {
    using namespace Algorithm;

    std::vector<Bridge> result;
    {
        result.reserve(bridge_expansions.size());
        uint32_t group_id = 0;
        using std::move_iterator;
        for (ExPolygon& expolygon : bridge_expolygons)
            result.push_back({ std::move(expolygon), group_id ++, bridge_expansions.end() });
    }


    // Detect overlaps of bridge anchors inside their respective shell regions.
    // bridge_expansions are sorted by boundary id and source id.
    for (auto expansion_iterator = bridge_expansions.begin(); expansion_iterator != bridge_expansions.end();) {
        auto boundary_region_begin = expansion_iterator;
        auto boundary_region_end = std::find_if(
            next(expansion_iterator),
            bridge_expansions.end(),
            [&](const RegionExpansionEx& expansion){
                return expansion.boundary_id != expansion_iterator->boundary_id;
            }
        );

        // Cache of bboxes per expansion boundary.
        std::vector<BoundingBox> bounding_boxes;
        bounding_boxes.reserve(std::distance(boundary_region_begin, boundary_region_end));
        std::transform(
            boundary_region_begin,
            boundary_region_end,
            std::back_inserter(bounding_boxes),
            [](const RegionExpansionEx& expansion){
                return get_extents(expansion.expolygon.contour);
            }
        );

        // For each bridge anchor of the current source:
        for (;expansion_iterator != boundary_region_end; ++expansion_iterator) {
            auto candidate_iterator = std::next(expansion_iterator);
            for (;candidate_iterator != boundary_region_end; ++candidate_iterator) {
                const BoundingBox& current_bounding_box{
                    bounding_boxes[expansion_iterator - boundary_region_begin]
                };
                const BoundingBox& candidate_bounding_box{
                    bounding_boxes[candidate_iterator - boundary_region_begin]
                };
                if (
                    expansion_iterator->src_id != candidate_iterator->src_id
                    && current_bounding_box.overlap(candidate_bounding_box)
                    // One may ignore holes, they are irrelevant for intersection test.
                    && !intersection(expansion_iterator->expolygon.contour, candidate_iterator->expolygon.contour).empty()
                ) {
                    // The two bridge regions intersect. Give them the same (lower) group id.
                    uint32_t id  = group_id(result, expansion_iterator->src_id);
                    uint32_t id2 = group_id(result, candidate_iterator->src_id);
                    if (id < id2)
                        result[id2].group_id = id;
                    else
                        result[id].group_id = id2;
                }
            }
        }
    }
    return result;
}

void detect_bridge_directions(
    const Algorithm::WaveSeeds& bridge_anchors,
    std::vector<Bridge>& bridges,
    const std::vector<ExpansionZone>& expansion_zones
) {
    if (expansion_zones.empty()) {
        throw std::runtime_error("At least one expansion zone must exist!");
    }
    auto it_bridge_anchor = bridge_anchors.begin();
    for (uint32_t bridge_id = 0; bridge_id < uint32_t(bridges.size()); ++ bridge_id) {
        Bridge &bridge = bridges[bridge_id];
        Polygons anchor_areas;
        int32_t last_anchor_id = -1;
        for (; it_bridge_anchor != bridge_anchors.end() && it_bridge_anchor->src == bridge_id; ++ it_bridge_anchor) {
            if (last_anchor_id != int(it_bridge_anchor->boundary)) {
                last_anchor_id = int(it_bridge_anchor->boundary);

                unsigned start_index{};
                unsigned end_index{};
                for (const ExpansionZone& expansion_zone: expansion_zones) {
                    end_index += expansion_zone.expolygons.size();
                    if (last_anchor_id < static_cast<int64_t>(end_index)) {
                        append(anchor_areas, to_polygons(expansion_zone.expolygons[last_anchor_id - start_index]));
                        break;
                    }
                    start_index += expansion_zone.expolygons.size();
                }
            }
        }
        Lines lines{to_lines(diff_pl(to_polylines(bridge.expolygon), expand(anchor_areas, float(SCALED_EPSILON))))};
        auto [bridging_dir, unsupported_dist] = detect_bridging_direction(lines, to_polygons(bridge.expolygon));
        bridge.angle = M_PI + std::atan2(bridging_dir.y(), bridging_dir.x());

        if constexpr (false) {
            coordf_t    stroke_width = scale_(0.06);
            BoundingBox bbox         = get_extents(anchor_areas);
            bbox.merge(get_extents(bridge.expolygon));
            bbox.offset(scale_(1.));
            ::Slic3r::SVG
                svg(debug_out_path(("bridge" + std::to_string(*bridge.angle) + "_" /* + std::to_string(this->layer()->bottom_z())*/).c_str()),
                bbox);
            svg.draw(bridge.expolygon, "cyan");
            svg.draw(lines, "green", stroke_width);
            svg.draw(anchor_areas, "red");
        }
    }
}

Surfaces merge_bridges(
    std::vector<Bridge>& bridges,
    const std::vector<Algorithm::RegionExpansionEx>& bridge_expansions,
    const float closing_radius
) {
    for (auto it = bridge_expansions.begin(); it != bridge_expansions.end(); ) {
        bridges[it->src_id].bridge_expansion_begin = it;
        uint32_t src_id = it->src_id;
        for (++ it; it != bridge_expansions.end() && it->src_id == src_id; ++ it) ;
    }

    Surfaces result;
    for (uint32_t bridge_id = 0; bridge_id < uint32_t(bridges.size()); ++ bridge_id) {
        if (group_id(bridges, bridge_id) == bridge_id) {
            // Head of the group.
            Polygons acc;
            for (uint32_t bridge_id2 = bridge_id; bridge_id2 < uint32_t(bridges.size()); ++ bridge_id2)
                if (group_id(bridges, bridge_id2) == bridge_id) {
                    append(acc, to_polygons(std::move(bridges[bridge_id2].expolygon)));
                    auto it_bridge_expansion = bridges[bridge_id2].bridge_expansion_begin;
                    assert(it_bridge_expansion == bridge_expansions.end() || it_bridge_expansion->src_id == bridge_id2);
                    for (; it_bridge_expansion != bridge_expansions.end() && it_bridge_expansion->src_id == bridge_id2; ++ it_bridge_expansion)
                        append(acc, to_polygons(it_bridge_expansion->expolygon));
                }
            //FIXME try to be smart and pick the best bridging angle for all?
            if (!bridges[bridge_id].angle) {
                assert(false && "Bridge angle must be pre-calculated!");
            }
            Surface templ{ stBottomBridge, {} };
            templ.bridge_angle = bridges[bridge_id].angle ? *bridges[bridge_id].angle : -1;
            //NOTE: The current regularization of the shells can create small unasigned regions in the object (E.G. benchy)
            // without the following closing operation, those regions will stay unfilled and cause small holes in the expanded surface.
            // look for narrow_ensure_vertical_wall_thickness_region_radius filter.
            ExPolygons final = closing_ex(acc, closing_radius);
            // without safety offset, artifacts are generated (GH #2494)
            // union_safety_offset_ex(acc)
            for (ExPolygon &ex : final)
                result.emplace_back(templ, std::move(ex));
        }
    }
    return result;
}

struct ExpansionResult {
    Algorithm::WaveSeeds anchors;
    std::vector<Algorithm::RegionExpansionEx> expansions;
};

ExpansionResult expand_expolygons(
    const ExPolygons& expolygons,
    std::vector<ExpansionZone>& expansion_zones
) {
    using namespace Algorithm;
    WaveSeeds bridge_anchors;
    std::vector<RegionExpansionEx> bridge_expansions;

    unsigned processed_bridges_count = 0;
    for (ExpansionZone& expansion_zone : expansion_zones) {
        WaveSeeds seeds{wave_seeds(
            expolygons,
            expansion_zone.expolygons,
            expansion_zone.parameters.tiny_expansion,
            true
        )};
        std::vector<RegionExpansionEx> expansions{propagate_waves_ex(
            seeds,
            expansion_zone.expolygons,
            expansion_zone.parameters
        )};

        for (WaveSeed &seed : seeds)
            seed.boundary += processed_bridges_count;
        for (RegionExpansionEx &expansion : expansions)
            expansion.boundary_id += processed_bridges_count;

        expansion_zone.expanded_into = ! expansions.empty();

        append(bridge_anchors, std::move(seeds));
        append(bridge_expansions, std::move(expansions));

        processed_bridges_count += expansion_zone.expolygons.size();
    }
    return {bridge_anchors, bridge_expansions};
}

// Extract bridging surfaces from "surfaces", expand them into "shells" using expansion_params,
// detect bridges.
// Trim "shells" by the expanded bridges.
Surfaces expand_bridges_detect_orientations(
    Surfaces &surfaces,
    std::vector<ExpansionZone>& expansion_zones,
    const float closing_radius
)
{
    using namespace Slic3r::Algorithm;

    double thickness;
    ExPolygons bridge_expolygons = fill_surfaces_extract_expolygons(surfaces, {stBottomBridge}, thickness);
    if (bridge_expolygons.empty())
        return {};

    // Calculate bridge anchors and their expansions in their respective shell region.
    ExpansionResult expansion_result{expand_expolygons(
        bridge_expolygons,
        expansion_zones
    )};

    std::vector<Bridge> bridges{get_grouped_bridges(
        std::move(bridge_expolygons),
        expansion_result.expansions
    )};
    bridge_expolygons.clear();

    std::sort(expansion_result.anchors.begin(), expansion_result.anchors.end(), Algorithm::lower_by_src_and_boundary);
    detect_bridge_directions(expansion_result.anchors, bridges, expansion_zones);

    // Merge the groups with the same group id, produce surfaces by merging source overhangs with their newly expanded anchors.
    std::sort(expansion_result.expansions.begin(), expansion_result.expansions.end(), [](auto &l, auto &r) {
        return l.src_id < r.src_id || (l.src_id == r.src_id && l.boundary_id < r.boundary_id);
    });
    Surfaces out{merge_bridges(bridges, expansion_result.expansions, closing_radius)};

    // Clip by the expanded bridges.
    for (ExpansionZone& expansion_zone : expansion_zones)
        if (expansion_zone.expanded_into)
            expansion_zone.expolygons = diff_ex(expansion_zone.expolygons, out);
    return out;
}

Surfaces expand_merge_surfaces(
    Surfaces &surfaces,
    SurfaceType surface_type,
    std::vector<ExpansionZone>& expansion_zones,
    const float closing_radius,
    const double bridge_angle = -1
)
{
    using namespace Slic3r::Algorithm;

    double thickness;
    ExPolygons src = fill_surfaces_extract_expolygons(surfaces, {surface_type}, thickness);
    if (src.empty())
        return {};

    unsigned processed_expolygons_count = 0;
    std::vector<RegionExpansion> expansions;
    for (ExpansionZone& expansion_zone : expansion_zones) {
        std::vector<RegionExpansion> zone_expansions = propagate_waves(src, expansion_zone.expolygons, expansion_zone.parameters);
        expansion_zone.expanded_into = !zone_expansions.empty();

        for (RegionExpansion &expansion : zone_expansions)
            expansion.boundary_id += processed_expolygons_count;

        processed_expolygons_count += expansion_zone.expolygons.size();
        append(expansions, std::move(zone_expansions));
    }

    std::vector<ExPolygon> expanded = merge_expansions_into_expolygons(std::move(src), std::move(expansions));
    //NOTE: The current regularization of the shells can create small unasigned regions in the object (E.G. benchy)
    // without the following closing operation, those regions will stay unfilled and cause small holes in the expanded surface.
    // look for narrow_ensure_vertical_wall_thickness_region_radius filter.
    expanded = closing_ex(expanded, closing_radius);
    // Trim the zones by the expanded expolygons.
    for (ExpansionZone& expansion_zone : expansion_zones)
        if (expansion_zone.expanded_into)
            expansion_zone.expolygons = diff_ex(expansion_zone.expolygons, expanded);

    Surface templ{ surface_type, {} };
    templ.bridge_angle = bridge_angle;
    Surfaces out;
    out.reserve(expanded.size());
    for (auto &expoly : expanded)
        out.emplace_back(templ, std::move(expoly));
    return out;
}

void LayerRegion::process_external_surfaces(const Layer *lower_layer, const Polygons *lower_layer_covered)
{
    using namespace Slic3r::Algorithm;

#ifdef SLIC3R_DEBUG_SLICE_PROCESSING
    export_region_fill_surfaces_to_svg_debug("4_process_external_surfaces-initial");
#endif /* SLIC3R_DEBUG_SLICE_PROCESSING */

    // Width of the perimeters.
    float shell_width = 0;
    float expansion_min = 0;
    if (int num_perimeters = this->region().config().wall_loops; num_perimeters > 0) {
        Flow external_perimeter_flow = this->flow(frExternalPerimeter);
        Flow perimeter_flow          = this->flow(frPerimeter);
        shell_width  = 0.5f * external_perimeter_flow.scaled_width() + external_perimeter_flow.scaled_spacing();
        shell_width += perimeter_flow.scaled_spacing() * (num_perimeters - 1);
        expansion_min = perimeter_flow.scaled_spacing();
    } else {
        // TODO: Maybe there is better solution when printing with zero perimeters, but this works reasonably well, given the situation
        shell_width   = float(SCALED_EPSILON);
        expansion_min = float(SCALED_EPSILON);;
    }

    // Scaled expansions of the respective external surfaces.
    float                           expansion_top           = shell_width * sqrt(2.);
    float                           expansion_bottom        = expansion_top;
    float                           expansion_bottom_bridge = expansion_top;
    // Expand by waves of expansion_step size (expansion_step is scaled), but with no more steps than max_nr_expansion_steps.
    const float                     expansion_step          = scaled<float>(0.1);
    // Don't take more than max_nr_steps for small expansion_step.
    static constexpr const size_t   max_nr_expansion_steps  = 5;
    // Radius (with added epsilon) to absorb empty regions emering from regularization of ensuring, viz  const float narrow_ensure_vertical_wall_thickness_region_radius = 0.5f * 0.65f * min_perimeter_infill_spacing;
    const float closing_radius = 0.55f * 0.65f * 1.05f * this->flow(frSolidInfill).scaled_spacing();

    // Expand the top / bottom / bridge surfaces into the shell thickness solid infills.
    double     layer_thickness;
    ExPolygons shells = union_ex(fill_surfaces_extract_expolygons(this->fill_surfaces.surfaces, { stInternalSolid }, layer_thickness));
    ExPolygons sparse = union_ex(fill_surfaces_extract_expolygons(this->fill_surfaces.surfaces, {stInternal}, layer_thickness));
    ExPolygons top_expolygons = union_ex(fill_surfaces_extract_expolygons(this->fill_surfaces.surfaces, {stTop}, layer_thickness));
    const auto expansion_params_into_sparse_infill = RegionExpansionParameters::build(expansion_min, expansion_step, max_nr_expansion_steps);
    const auto expansion_params_into_solid_infill  = RegionExpansionParameters::build(expansion_bottom_bridge, expansion_step, max_nr_expansion_steps);

    std::vector<ExpansionZone> expansion_zones{
        ExpansionZone{std::move(shells), expansion_params_into_solid_infill},
        ExpansionZone{std::move(sparse), expansion_params_into_sparse_infill},
        ExpansionZone{std::move(top_expolygons), expansion_params_into_solid_infill},
    };

    SurfaceCollection bridges;
    {
        BOOST_LOG_TRIVIAL(trace) << "Processing external surface, detecting bridges. layer" << this->layer()->print_z;
        // ORCA: Relative/Align Bridge Angle
        const auto  &region_config    = this->region().config();
        const double custom_angle_deg = region_config.bridge_angle.value;
        const bool   relative_angle   = region_config.relative_bridge_angle.value;
        const double custom_angle_rad = Geometry::deg2rad(custom_angle_deg);

        double align_offset_rad = 0.0;
        if (region_config.align_infill_direction_to_model) {
            auto m = this->layer()->object()->trafo().matrix();
            align_offset_rad = std::atan2((double)m(1, 0), (double)m(0, 0));
        }

        bridges.surfaces = (custom_angle_deg > 0.0 && !relative_angle) ?
            expand_merge_surfaces(this->fill_surfaces.surfaces, stBottomBridge, expansion_zones, closing_radius, custom_angle_rad + align_offset_rad) :
            expand_bridges_detect_orientations(this->fill_surfaces.surfaces, expansion_zones, closing_radius);
        if (custom_angle_deg > 0.0 && relative_angle) {
            for (Surface &bridge_surface : bridges.surfaces) {
                if (bridge_surface.bridge_angle >= 0)
                    bridge_surface.bridge_angle += custom_angle_rad;
            }
        }
        BOOST_LOG_TRIVIAL(trace) << "Processing external surface, detecting bridges - done";
#ifdef SLIC3R_DEBUG_SLICE_PROCESSING
        {
            static int iRun = 0;
            bridges.export_to_svg(debug_out_path("bridges-after-grouping-%d.svg", iRun++).c_str(), true);
        }
#endif
    }

    this->fill_surfaces.remove_types({stTop});
    {
        Surface top_templ(stTop, {});
        top_templ.thickness = layer_thickness;
        this->fill_surfaces.append(std::move(expansion_zones.back().expolygons), top_templ);
    }

    expansion_zones.pop_back();

    expansion_zones.at(0).parameters = RegionExpansionParameters::build(expansion_bottom, expansion_step, max_nr_expansion_steps);
    Surfaces bottoms = expand_merge_surfaces(this->fill_surfaces.surfaces, stBottom, expansion_zones, closing_radius);

    expansion_zones.at(0).parameters = RegionExpansionParameters::build(expansion_top, expansion_step, max_nr_expansion_steps);
    Surfaces tops = expand_merge_surfaces(this->fill_surfaces.surfaces, stTop, expansion_zones, closing_radius);

    // turn too small internal regions into solid regions according to the user setting
    if (!this->layer()->object()->print()->config().spiral_mode && this->region().config().sparse_infill_density.value > 0) {
        // scaling an area requires two calls!
        double min_area = scale_(scale_(this->region().config().minimum_sparse_infill_area.value));
        ExPolygons small_regions{};
        expansion_zones[1].expolygons.erase(std::remove_if(expansion_zones[1].expolygons.begin(), expansion_zones[1].expolygons.end(), [min_area, &small_regions](ExPolygon& ex_polygon) {
            if (ex_polygon.area() <= min_area) {
                small_regions.push_back(ex_polygon);
                return true;
            }
            return false;
        }), expansion_zones[1].expolygons.end());

        if (!small_regions.empty()) {
            expansion_zones[0].expolygons = union_ex(expansion_zones[0].expolygons, small_regions);
        }
    }

//    this->fill_surfaces.remove_types({ stBottomBridge, stBottom, stTop, stInternal, stInternalSolid });
    this->fill_surfaces.clear();
    unsigned zones_expolygons_count = 0;
    for (const ExpansionZone& zone : expansion_zones)
        zones_expolygons_count += zone.expolygons.size();
    reserve_more(this->fill_surfaces.surfaces, zones_expolygons_count + bridges.size() + bottoms.size() + tops.size());
    {
        Surface solid_templ(stInternalSolid, {});
        solid_templ.thickness = layer_thickness;
        this->fill_surfaces.append(std::move(expansion_zones[0].expolygons), solid_templ);
    }
    {
        Surface sparse_templ(stInternal, {});
        sparse_templ.thickness = layer_thickness;
        this->fill_surfaces.append(std::move(expansion_zones[1].expolygons), sparse_templ);
    }
    this->fill_surfaces.append(std::move(bridges.surfaces));
    this->fill_surfaces.append(std::move(bottoms));
    this->fill_surfaces.append(std::move(tops));

#ifdef SLIC3R_DEBUG_SLICE_PROCESSING
    export_region_fill_surfaces_to_svg_debug("4_process_external_surfaces-final");
#endif /* SLIC3R_DEBUG_SLICE_PROCESSING */
}
#else

//#define EXTERNAL_SURFACES_OFFSET_PARAMETERS ClipperLib::jtMiter, 3.
//#define EXTERNAL_SURFACES_OFFSET_PARAMETERS ClipperLib::jtMiter, 1.5
#define EXTERNAL_SURFACES_OFFSET_PARAMETERS ClipperLib::jtSquare, 0.

void LayerRegion::process_external_surfaces(const Layer *lower_layer, const Polygons *lower_layer_covered)
{
    const bool      has_infill = this->region().config().sparse_infill_density.value > 0.;
    //BBS
    auto nozzle_diameter = this->region().nozzle_dmr_avg(this->layer()->object()->print()->config());
    const float margin = float(scale_(EXTERNAL_INFILL_MARGIN));
    const float bridge_margin = std::min(float(scale_(BRIDGE_INFILL_MARGIN)), float(scale_(nozzle_diameter * BRIDGE_INFILL_MARGIN / 0.4)));

    // BBS
    const PrintObjectConfig& object_config = this->layer()->object()->config();

#ifdef SLIC3R_DEBUG_SLICE_PROCESSING
    export_region_fill_surfaces_to_svg_debug("3_process_external_surfaces-initial");
#endif /* SLIC3R_DEBUG_SLICE_PROCESSING */

    // 1) Collect bottom and bridge surfaces, each of them grown by a fixed 3mm offset
    // for better anchoring.
    // Bottom surfaces, grown.
    Surfaces                    bottom;
    // Bridge surfaces, initialy not grown.
    Surfaces                    bridges;
    // Top surfaces, grown.
    Surfaces                    top;
    // Internal surfaces, not grown.
    Surfaces                    internal;
    // Areas, where an infill of various types (top, bottom, bottom bride, sparse, void) could be placed.
    Polygons                    fill_boundaries = to_polygons(this->fill_expolygons);
    Polygons  					lower_layer_covered_tmp;

    // Collect top surfaces and internal surfaces.
    // Collect fill_boundaries: If we're slicing with no infill, we can't extend external surfaces over non-existent infill.
    // This loop destroys the surfaces (aliasing this->fill_surfaces.surfaces) by moving into top/internal/fill_boundaries!

    {
        // Voids are sparse infills if infill rate is zero.
        Polygons voids;

        double max_grid_area = -1;
        if (this->layer()->lower_layer != nullptr)
            max_grid_area = this->layer()->lower_layer->get_sparse_infill_max_void_area();
        for (const Surface &surface : this->fill_surfaces.surfaces) {
            if (surface.is_top()) {
                // Collect the top surfaces, inflate them and trim them by the bottom surfaces.
                // This gives the priority to bottom surfaces.
                if (max_grid_area < 0 || surface.expolygon.area() < max_grid_area)
                    surfaces_append(top, offset_ex(surface.expolygon, margin, EXTERNAL_SURFACES_OFFSET_PARAMETERS), surface);
                else
                    //BBS: Don't need to expand too much in this situation. Expand 3mm to eliminate hole and 1mm for contour
                    surfaces_append(top, intersection_ex(offset(surface.expolygon.contour, margin / 3.0, EXTERNAL_SURFACES_OFFSET_PARAMETERS),
                                                         offset_ex(surface.expolygon, margin, EXTERNAL_SURFACES_OFFSET_PARAMETERS)), surface);
            } else if (surface.surface_type == stBottom || (surface.surface_type == stBottomBridge && lower_layer == nullptr)) {
                // Grown by 3mm.
                surfaces_append(bottom, offset_ex(surface.expolygon, margin, EXTERNAL_SURFACES_OFFSET_PARAMETERS), surface);
            } else if (surface.surface_type == stBottomBridge) {
                if (! surface.empty())
                    bridges.emplace_back(surface);
            }
            if (surface.is_internal()) {
            	assert(surface.surface_type == stInternal || surface.surface_type == stInternalSolid);
            	if (! has_infill && lower_layer != nullptr)
            		polygons_append(voids, surface.expolygon);
            	internal.emplace_back(std::move(surface));
            }
        }
        if (! has_infill && lower_layer != nullptr && ! voids.empty()) {
        	// Remove voids from fill_boundaries, that are not supported by the layer below.
            if (lower_layer_covered == nullptr) {
            	lower_layer_covered = &lower_layer_covered_tmp;
            	lower_layer_covered_tmp = to_polygons(lower_layer->lslices);
            }
            if (! lower_layer_covered->empty())
            	voids = diff(voids, *lower_layer_covered);
            fill_boundaries = diff(fill_boundaries, voids);
        }
    }

#if 0
    {
        static int iRun = 0;
        bridges.export_to_svg(debug_out_path("bridges-before-grouping-%d.svg", iRun ++), true);
    }
#endif

    if (bridges.empty())
    {
        fill_boundaries = union_safety_offset(fill_boundaries);
    } else
    {
        // 1) Calculate the inflated bridge regions, each constrained to its island.
        ExPolygons               fill_boundaries_ex = union_safety_offset_ex(fill_boundaries);
        std::vector<Polygons>    bridges_grown;
        std::vector<BoundingBox> bridge_bboxes;

#ifdef SLIC3R_DEBUG_SLICE_PROCESSING
        {
            static int iRun = 0;
            SVG svg(debug_out_path("3_process_external_surfaces-fill_regions-%d.svg", iRun ++).c_str(), get_extents(fill_boundaries_ex));
            svg.draw(fill_boundaries_ex);
            svg.draw_outline(fill_boundaries_ex, "black", "blue", scale_(0.05)); 
            svg.Close();
        }

//        export_region_fill_surfaces_to_svg_debug("3_process_external_surfaces-initial");
#endif /* SLIC3R_DEBUG_SLICE_PROCESSING */
 
        {
            // Bridge expolygons, grown, to be tested for intersection with other bridge regions.
            std::vector<BoundingBox> fill_boundaries_ex_bboxes = get_extents_vector(fill_boundaries_ex);
            bridges_grown.reserve(bridges.size());
            bridge_bboxes.reserve(bridges.size());
            for (size_t i = 0; i < bridges.size(); ++ i) {
                // Find the island of this bridge.
                const Point pt = bridges[i].expolygon.contour.points.front();
                int idx_island = -1;
                for (int j = 0; j < int(fill_boundaries_ex.size()); ++ j)
                    if (fill_boundaries_ex_bboxes[j].contains(pt) && 
                        fill_boundaries_ex[j].contains(pt)) {
                        idx_island = j;
                        break;
                    }
                // Grown by 3mm.
                //BBS: eliminate too narrow area to avoid generating bridge on top layer when wall loop is 1
                //Polygons polys = offset(bridges[i].expolygon, bridge_margin, EXTERNAL_SURFACES_OFFSET_PARAMETERS);
                Polygons polys = offset2({ bridges[i].expolygon }, -scale_(nozzle_diameter * 0.1), bridge_margin, EXTERNAL_SURFACES_OFFSET_PARAMETERS);
                if (idx_island == -1) {
				    BOOST_LOG_TRIVIAL(trace) << "Bridge did not fall into the source region!";
                } else {
                    // Found an island, to which this bridge region belongs. Trim it,
                    polys = intersection(polys, fill_boundaries_ex[idx_island]);
                }
                bridge_bboxes.push_back(get_extents(polys));
                bridges_grown.push_back(std::move(polys));
            }
        }

        // 2) Group the bridge surfaces by overlaps.
        std::vector<size_t> bridge_group(bridges.size(), (size_t)-1);
        size_t n_groups = 0; 
        for (size_t i = 0; i < bridges.size(); ++ i) {
            // A grup id for this bridge.
            size_t group_id = (bridge_group[i] == size_t(-1)) ? (n_groups ++) : bridge_group[i];
            bridge_group[i] = group_id;
            // For all possibly overlaping bridges:
            for (size_t j = i + 1; j < bridges.size(); ++ j) {
                if (! bridge_bboxes[i].overlap(bridge_bboxes[j]))
                    continue;
                if (intersection(bridges_grown[i], bridges_grown[j]).empty())
                    continue;
                // The two bridge regions intersect. Give them the same group id.
                if (bridge_group[j] != size_t(-1)) {
                    // The j'th bridge has been merged with some other bridge before.
                    size_t group_id_new = bridge_group[j];
                    for (size_t k = 0; k < j; ++ k)
                        if (bridge_group[k] == group_id)
                            bridge_group[k] = group_id_new;
                    group_id = group_id_new;
                }
                bridge_group[j] = group_id;
            }
        }

        // 3) Merge the groups with the same group id, detect bridges.
        {
			BOOST_LOG_TRIVIAL(trace) << "Processing external surface, detecting bridges. layer" << this->layer()->print_z << ", bridge groups: " << n_groups;
            for (size_t group_id = 0; group_id < n_groups; ++ group_id) {
                size_t n_bridges_merged = 0;
                size_t idx_last = (size_t)-1;
                for (size_t i = 0; i < bridges.size(); ++ i) {
                    if (bridge_group[i] == group_id) {
                        ++ n_bridges_merged;
                        idx_last = i;
                    }
                }
                if (n_bridges_merged == 0)
                    // This group has no regions assigned as these were moved into another group.
                    continue;
                // Collect the initial ungrown regions and the grown polygons.
                ExPolygons  initial;
                Polygons    grown;
                for (size_t i = 0; i < bridges.size(); ++ i) {
                    if (bridge_group[i] != group_id)
                        continue;
                    initial.push_back(std::move(bridges[i].expolygon));
                    polygons_append(grown, bridges_grown[i]);
                }
                // detect bridge direction before merging grown surfaces otherwise adjacent bridges
                // would get merged into a single one while they need different directions
                // also, supply the original expolygon instead of the grown one, because in case
                // of very thin (but still working) anchors, the grown expolygon would go beyond them
                // ORCA: Relative/Align Bridge Angle
                const auto &region_config   = this->region().config();
                const double custom_angle_deg = region_config.bridge_angle.value;
                const bool   relative_angle   = region_config.relative_bridge_angle.value;
                const double custom_angle_rad = Geometry::deg2rad(custom_angle_deg);

                double align_offset_rad = 0.0;
                if (region_config.align_infill_direction_to_model) {
                    auto m = this->layer()->object()->trafo().matrix();
                    align_offset_rad = std::atan2((double)m(1, 0), (double)m(0, 0));
                }

                if (custom_angle_deg > 0.0 && !relative_angle) {
                    bridges[idx_last].bridge_angle = custom_angle_rad + align_offset_rad;
                } else {
                    auto [bridging_dir, unsupported_dist] = detect_bridging_direction(to_polygons(initial), to_polygons(lower_layer->lslices));
                    bridges[idx_last].bridge_angle = PI + std::atan2(bridging_dir.y(), bridging_dir.x());
                    if (custom_angle_deg > 0.0 && relative_angle)
                        bridges[idx_last].bridge_angle += custom_angle_rad;
                }

                /*
                BridgeDetector bd(initial, lower_layer->lslices, this->bridging_flow(frInfill, object_config.thick_bridges).scaled_width());
                #ifdef SLIC3R_DEBUG
                printf("Processing bridge at layer %zu:\n", this->layer()->id());
                #endif
                //BBS: use 0 as custom angle to enable auto detection all the time
                double custom_angle = Geometry::deg2rad(this->region().config().bridge_angle.value);
                if(custom_angle > 0)
                        bridges[idx_last].bridge_angle = custom_angle;
				else if (bd.detect_angle(custom_angle)) {
                    bridges[idx_last].bridge_angle = bd.angle;
                    if (this->layer()->object()->has_support()) {
//                        polygons_append(this->bridged, bd.coverage());
                        append(this->unsupported_bridge_edges, bd.unsupported_edges());
                    }
				} else if (custom_angle > 0) {
					// Bridge was not detected (likely it is only supported at one side). Still it is a surface filled in
					// using a bridging flow, therefore it makes sense to respect the custom bridging direction.
					bridges[idx_last].bridge_angle = custom_angle;
				}
                */
                // without safety offset, artifacts are generated (GH #2494)
                surfaces_append(bottom, union_safety_offset_ex(grown), bridges[idx_last]);
            }

            fill_boundaries = to_polygons(fill_boundaries_ex);
			BOOST_LOG_TRIVIAL(trace) << "Processing external surface, detecting bridges - done";
		}

    #if 0
        {
            static int iRun = 0;
            bridges.export_to_svg(debug_out_path("bridges-after-grouping-%d.svg", iRun ++), true);
        }
    #endif
    }

    Surfaces new_surfaces;
    {
        // Merge top and bottom in a single collection.
        surfaces_append(top, std::move(bottom));
        // Intersect the grown surfaces with the actual fill boundaries.
        Polygons bottom_polygons = to_polygons(bottom);
        for (size_t i = 0; i < top.size(); ++ i) {
            Surface &s1 = top[i];
            if (s1.empty())
                continue;
            Polygons polys;
            polygons_append(polys, to_polygons(std::move(s1)));
            for (size_t j = i + 1; j < top.size(); ++ j) {
                Surface &s2 = top[j];
                if (! s2.empty() && surfaces_could_merge(s1, s2)) {
                    polygons_append(polys, to_polygons(std::move(s2)));
                    s2.clear();
                }
            }
            if (s1.is_top())
                // Trim the top surfaces by the bottom surfaces. This gives the priority to the bottom surfaces.
                polys = diff(polys, bottom_polygons);
            surfaces_append(
                new_surfaces,
                // Don't use a safety offset as fill_boundaries were already united using the safety offset.
                intersection_ex(polys, fill_boundaries),
                s1);
        }
    }
    
    // Subtract the new top surfaces from the other non-top surfaces and re-add them.
    Polygons new_polygons = to_polygons(new_surfaces);
    for (size_t i = 0; i < internal.size(); ++ i) {
        Surface &s1 = internal[i];
        if (s1.empty())
            continue;
        Polygons polys;
        polygons_append(polys, to_polygons(std::move(s1)));
        for (size_t j = i + 1; j < internal.size(); ++ j) {
            Surface &s2 = internal[j];
            if (! s2.empty() && surfaces_could_merge(s1, s2)) {
                polygons_append(polys, to_polygons(std::move(s2)));
                s2.clear();
            }
        }
        ExPolygons new_expolys = diff_ex(polys, new_polygons);
        polygons_append(new_polygons, to_polygons(new_expolys));
        surfaces_append(new_surfaces, std::move(new_expolys), s1);
    }
    
    this->fill_surfaces.surfaces = std::move(new_surfaces);

#ifdef SLIC3R_DEBUG_SLICE_PROCESSING
    export_region_fill_surfaces_to_svg_debug("3_process_external_surfaces-final");
#endif /* SLIC3R_DEBUG_SLICE_PROCESSING */
}
#endif

void LayerRegion::prepare_fill_surfaces()
{
#ifdef SLIC3R_DEBUG_SLICE_PROCESSING
    export_region_slices_to_svg_debug("2_prepare_fill_surfaces-initial");
    export_region_fill_surfaces_to_svg_debug("2_prepare_fill_surfaces-initial");
#endif /* SLIC3R_DEBUG_SLICE_PROCESSING */ 

    /*  Note: in order to make the psPrepareInfill step idempotent, we should never
        alter fill_surfaces boundaries on which our idempotency relies since that's
        the only meaningful information returned by psPerimeters. */
    
    bool spiral_mode = this->layer()->object()->print()->config().spiral_mode;

    // if no solid layers are requested, turn top/bottom surfaces to internal
    if (! spiral_mode && this->region().config().top_shell_layers == 0) {
        for (Surface &surface : this->fill_surfaces.surfaces)
            if (surface.is_top())
                //BBS
                //surface.surface_type = this->layer()->object()->config().infill_only_where_needed ? stInternalVoid : stInternal;
                surface.surface_type = PrintObject::infill_only_where_needed ? stInternalVoid : stInternal;
    }
    if (this->region().config().bottom_shell_layers == 0) {
        for (Surface &surface : this->fill_surfaces.surfaces)
            if (surface.is_bottom()) // (surface.surface_type == stBottom)
                surface.surface_type = stInternal;
    }

    if (!spiral_mode && fabs(this->region().config().sparse_infill_density.value - 100.) < EPSILON) {
        // Turn all internal sparse infill into solid infill, if sparse_infill_density is 100%
        for (Surface &surface : this->fill_surfaces.surfaces)
            if (surface.surface_type == stInternal)
                surface.surface_type = stInternalSolid;
    }

#ifdef SLIC3R_DEBUG_SLICE_PROCESSING
    export_region_slices_to_svg_debug("2_prepare_fill_surfaces-final");
    export_region_fill_surfaces_to_svg_debug("2_prepare_fill_surfaces-final");
#endif /* SLIC3R_DEBUG_SLICE_PROCESSING */
}

double LayerRegion::infill_area_threshold() const
{
    double ss = this->flow(frSolidInfill).scaled_spacing();
    return ss*ss;
}

void LayerRegion::trim_surfaces(const Polygons &trimming_polygons)
{
#ifndef NDEBUG
    for (const Surface &surface : this->slices.surfaces)
        assert(surface.surface_type == stInternal);
#endif /* NDEBUG */
	this->slices.set(intersection_ex(this->slices.surfaces, trimming_polygons), stInternal);
}

void LayerRegion::elephant_foot_compensation_step(const float elephant_foot_compensation_perimeter_step, const Polygons &trimming_polygons)
{
#ifndef NDEBUG
    for (const Surface &surface : this->slices.surfaces)
        assert(surface.surface_type == stInternal);
#endif /* NDEBUG */
    Polygons tmp = intersection(this->slices.surfaces, trimming_polygons);
    append(tmp, diff(this->slices.surfaces, opening(this->slices.surfaces, elephant_foot_compensation_perimeter_step)));
    this->slices.set(union_ex(tmp), stInternal);
}

void LayerRegion::export_region_slices_to_svg(const char *path) const
{
    BoundingBox bbox;
    for (Surfaces::const_iterator surface = this->slices.surfaces.begin(); surface != this->slices.surfaces.end(); ++surface)
        bbox.merge(get_extents(surface->expolygon));
    Point legend_size = export_surface_type_legend_to_svg_box_size();
    Point legend_pos(bbox.min(0), bbox.max(1));
    bbox.merge(Point(std::max(bbox.min(0) + legend_size(0), bbox.max(0)), bbox.max(1) + legend_size(1)));

    SVG svg(path, bbox);
    const float transparency = 0.5f;
    for (Surfaces::const_iterator surface = this->slices.surfaces.begin(); surface != this->slices.surfaces.end(); ++surface)
        svg.draw(surface->expolygon, surface_type_to_color_name(surface->surface_type), transparency);
    for (Surfaces::const_iterator surface = this->fill_surfaces.surfaces.begin(); surface != this->fill_surfaces.surfaces.end(); ++surface)
        svg.draw(surface->expolygon.lines(), surface_type_to_color_name(surface->surface_type));
    export_surface_type_legend_to_svg(svg, legend_pos);
    svg.Close();
}

// Export to "out/LayerRegion-name-%d.svg" with an increasing index with every export.
void LayerRegion::export_region_slices_to_svg_debug(const char *name) const
{
    static std::map<std::string, size_t> idx_map;
    size_t &idx = idx_map[name];
    this->export_region_slices_to_svg(debug_out_path("LayerRegion-slices-%s-%d.svg", name, idx ++).c_str());
}

void LayerRegion::export_region_fill_surfaces_to_svg(const char *path) const
{
    BoundingBox bbox;
    for (Surfaces::const_iterator surface = this->fill_surfaces.surfaces.begin(); surface != this->fill_surfaces.surfaces.end(); ++surface)
        bbox.merge(get_extents(surface->expolygon));
    Point legend_size = export_surface_type_legend_to_svg_box_size();
    Point legend_pos(bbox.min(0), bbox.max(1));
    bbox.merge(Point(std::max(bbox.min(0) + legend_size(0), bbox.max(0)), bbox.max(1) + legend_size(1)));

    SVG svg(path, bbox);
    const float transparency = 0.5f;
    for (const Surface &surface : this->fill_surfaces.surfaces) {
        svg.draw(surface.expolygon, surface_type_to_color_name(surface.surface_type), transparency);
        svg.draw_outline(surface.expolygon, "black", "blue", scale_(0.05)); 
    }
    export_surface_type_legend_to_svg(svg, legend_pos);
    svg.Close();
}

// Export to "out/LayerRegion-name-%d.svg" with an increasing index with every export.
void LayerRegion::export_region_fill_surfaces_to_svg_debug(const char *name) const
{
    static std::map<std::string, size_t> idx_map;
    size_t &idx = idx_map[name];
    this->export_region_fill_surfaces_to_svg(debug_out_path("LayerRegion-fill_surfaces-%s-%d.svg", name, idx ++).c_str());
}

void LayerRegion::simplify_entity_collection(ExtrusionEntityCollection* entity_collection)
{
    for (size_t i = 0; i < entity_collection->entities.size(); i++) {
        if (ExtrusionEntityCollection* collection = dynamic_cast<ExtrusionEntityCollection*>(entity_collection->entities[i]))
            this->simplify_entity_collection(collection);
        else if (ExtrusionPath* path = dynamic_cast<ExtrusionPath*>(entity_collection->entities[i]))
            this->simplify_path(path);
        else if (ExtrusionMultiPath* multipath = dynamic_cast<ExtrusionMultiPath*>(entity_collection->entities[i]))
            this->simplify_multi_path(multipath);
        else if (ExtrusionLoop* loop = dynamic_cast<ExtrusionLoop*>(entity_collection->entities[i]))
            this->simplify_loop(loop);
        else
            throw Slic3r::InvalidArgument("Invalid extrusion entity supplied to simplify_entity_collection()");
    }
}

void LayerRegion::simplify_path(ExtrusionPath* path)
{
    const auto print_config = this->layer()->object()->print()->config();
    const bool spiral_mode = print_config.spiral_mode;
    const bool enable_arc_fitting = print_config.enable_arc_fitting;
    const auto scaled_resolution = scaled<double>(print_config.resolution.value);

    // Arc-overhang paths have already been clipped against previously emitted
    // paths. Preserve fitting data approved by the generator, but never fit or
    // simplify one here: doing so without the generator's obstacle context may
    // bow a segment outside its validated corridor.
    if (is_arc_fill(path->role())) {
        return;
    }

    if (enable_arc_fitting &&
        !spiral_mode) {
        if (path->role() == erInternalInfill)
            path->simplify_by_fitting_arc(SCALED_SPARSE_INFILL_RESOLUTION);
        else
            path->simplify_by_fitting_arc(scaled_resolution);
    } else {
        path->simplify(scaled_resolution);
    }
}

void LayerRegion::simplify_multi_path(ExtrusionMultiPath* multipath)
{
    const auto print_config = this->layer()->object()->print()->config();
    const bool spiral_mode = print_config.spiral_mode;
    const bool enable_arc_fitting = print_config.enable_arc_fitting;
    const auto scaled_resolution = scaled<double>(print_config.resolution.value);

    for (size_t i = 0; i < multipath->paths.size(); ++i) {
        if (is_arc_fill(multipath->paths[i].role())) {
            continue;
        }
        if (enable_arc_fitting &&
            !spiral_mode) {
            if (multipath->paths[i].role() == erInternalInfill)
                multipath->paths[i].simplify_by_fitting_arc(SCALED_SPARSE_INFILL_RESOLUTION);
            else
                multipath->paths[i].simplify_by_fitting_arc(scaled_resolution);
        } else {
            multipath->paths[i].simplify(scaled_resolution);
        }
    }
}

void LayerRegion::simplify_loop(ExtrusionLoop* loop)
{
    const auto print_config = this->layer()->object()->print()->config();
    const bool spiral_mode = print_config.spiral_mode;
    const bool enable_arc_fitting = print_config.enable_arc_fitting;
    const auto scaled_resolution = scaled<double>(print_config.resolution.value);

    for (size_t i = 0; i < loop->paths.size(); ++i) {
        if (is_arc_fill(loop->paths[i].role())) {
            continue;
        }
        if (enable_arc_fitting &&
            !spiral_mode) {
            if (loop->paths[i].role() == erInternalInfill)
                loop->paths[i].simplify_by_fitting_arc(SCALED_SPARSE_INFILL_RESOLUTION);
            else
                loop->paths[i].simplify_by_fitting_arc(scaled_resolution);
        } else {
            loop->paths[i].simplify(scaled_resolution);
        }
    }
}

}
 
