#include "ClipperUtils.hpp"
#include "Arachne/WallToolPaths.hpp"
#include "ExtrusionEntity.hpp"
#include "ExtrusionEntityCollection.hpp"
#include "Flow.hpp"
#include "Geometry.hpp"
#include "Layer.hpp"
#include "NonplanarSurface.hpp"
#include "Print.hpp"
#include "SequentialGantryGeometry.hpp"
#include "ShortestPath.hpp"
#include "SLA/IndexedMesh.hpp"
#include "TriangleMesh.hpp"
#include "Utils.hpp"
#include "VariableWidth.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cmath>
#include <map>
#include <numeric>
#include <optional>
#include <queue>
#include <type_traits>

#include <boost/log/trivial.hpp>

#include <tbb/blocked_range.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/parallel_for.h>

namespace Slic3r {

Polylines nonplanar_wall_centerlines(const ExtrusionEntity &entity, int inset)
{
    Polylines result;
    const auto collect = [&](const auto &self, const ExtrusionEntity &item, int parent_inset) -> void {
        if (const auto *path = dynamic_cast<const ExtrusionPath *>(&item)) {
            const int path_inset = parent_inset >= 0 ? parent_inset : path->inset_idx;
            if (is_perimeter(path->role()) && path_inset == inset && path->polyline.points.size() >= 2) {
                Polyline line;
                line.points.reserve(path->polyline.points.size());
                for (const Point3 &point : path->polyline.points)
                    line.points.push_back(point.to_point());
                line.remove_duplicate_points();
                if (line.points.size() >= 2)
                    result.emplace_back(std::move(line));
            }
            return;
        }
        const int owner_inset = item.inset_idx >= 0 ? item.inset_idx : parent_inset;
        if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&item)) {
            for (const ExtrusionPath &path : loop->paths)
                self(self, path, owner_inset);
        } else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&item)) {
            for (const ExtrusionPath &path : multipath->paths)
                self(self, path, owner_inset);
        } else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&item)) {
            for (const ExtrusionEntity *child : collection->entities)
                if (child != nullptr)
                    self(self, *child, owner_inset);
        }
    };
    collect(collect, entity, -1);
    return result;
}

bool filter_nonplanar_course_paths(ExtrusionEntitiesPtr &paths,
                                   const std::vector<uint8_t> &retained)
{
    if (paths.size() != retained.size())
        return false;
    // A failed collision check is not proof that a wall is redundant. Test
    // the entire selection before deleting anything, including repair fill.
    for (size_t i = 0; i < paths.size(); ++i)
        if (!retained[i] && (paths[i] == nullptr || is_perimeter(paths[i]->role())))
            return false;
    for (size_t i = paths.size(); i-- > 0;)
        if (!retained[i]) {
            delete paths[i];
            paths.erase(paths.begin() + i);
        }
    return true;
}

size_t nonplanar_solid_course_begin(
    const std::vector<NonplanarCourseProfile> &courses,
    double minimum_surface_z, double maximum_surface_z,
    size_t minimum_layers, double minimum_thickness)
{
    if (courses.empty())
        return 0;
    size_t begin = courses.size() -
        std::min(courses.size(), std::max(size_t(1), minimum_layers));
    const auto &top = courses.back();
    while (begin > 0) {
        const auto &below = courses[begin - 1];
        const double high_thickness = top.plane_z - below.plane_z;
        const double low_thickness = high_thickness +
            (top.blend - below.blend) * (minimum_surface_z - maximum_surface_z);
        // Thickness is affine in surface Z, so the two extrema bound the
        // entire patch. Counting nominal owners overestimates thin low-side
        // beads and can leave less solid material than the profile requests.
        if (std::min(low_thickness, high_thickness) + EPSILON >= minimum_thickness)
            break;
        --begin;
    }
    return begin;
}

double nonplanar_wall_course_z(double planar_z, double projected_z,
    double surface_z, double course_z, bool fully_selected_path)
{
    if (fully_selected_path)
        return course_z;
    const double surface_delta = surface_z - planar_z;
    // Zero displacement is the planar attachment, not full selection. In
    // particular, a surface crossing the owner plane cannot justify moving
    // that shared endpoint onto a different transition course.
    const double influence = std::abs(surface_delta) <= EPSILON ? 0. :
        std::clamp((projected_z - planar_z) / surface_delta, 0., 1.);
    return planar_z + influence * (course_z - planar_z);
}

std::vector<std::pair<double, double>> nonplanar_segment_overlap(
    const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &d,
    double xy_tolerance, double z_tolerance)
{
    // Most pairs have disjoint bead bounds. Reject them before constructing
    // the closest-point intervals (and their temporary vectors).
    for (int axis = 0; axis < 3; ++axis) {
        const double tolerance = axis == 2 ? z_tolerance : xy_tolerance;
        if (std::min(a[axis], b[axis]) > std::max(c[axis], d[axis]) + tolerance ||
            std::max(a[axis], b[axis]) < std::min(c[axis], d[axis]) - tolerance)
            return {};
    }
    const Vec3d r = b - a;
    const Vec3d s = d - c;
    const double length_squared = s.head<2>().squaredNorm();
    // u(t) is the support-segment parameter. First restrict it to the Z
    // band, then minimize XY distance within that interval. Checking Z only
    // at the unconstrained XY-nearest point misses real contact on a slope.
    struct Affine {
        double origin, slope;
        double at(double t) const { return origin + slope * t; }
    };
    const Affine zero{0., 0.}, one{1., 0.};
    const Affine nearest{
        length_squared > 1e-16 ? (a - c).head<2>().dot(s.head<2>()) / length_squared : 0.,
        length_squared > 1e-16 ? r.head<2>().dot(s.head<2>()) / length_squared : 0.};
    const bool sloped_support = std::abs(s.z()) > 1e-16;
    Affine z_low = zero, z_high = one;
    if (sloped_support) {
        z_low = {(a.z() - c.z() - z_tolerance) / s.z(), r.z() / s.z()};
        z_high = {(a.z() - c.z() + z_tolerance) / s.z(), r.z() / s.z()};
        if (s.z() < 0.)
            std::swap(z_low, z_high);
    }
    // All changes in the active clamp bound occur at intersections of these
    // affine functions. Within each interval the closest feasible point is
    // affine, and its XY contact inequality is an exact quadratic.
    const std::array<Affine, 5> bounds{zero, one, nearest, z_low, z_high};
    std::vector<double> cuts{0., 1.};
    for (size_t i = 0; i < bounds.size(); ++i)
        for (size_t j = i + 1; j < bounds.size(); ++j) {
            const double slope = bounds[i].slope - bounds[j].slope;
            if (std::abs(slope) <= 1e-16)
                continue;
            const double t = (bounds[j].origin - bounds[i].origin) / slope;
            if (t > 0. && t < 1.)
                cuts.push_back(t);
        }
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    std::vector<std::pair<double, double>> intervals;
    for (size_t idx = 1; idx < cuts.size(); ++idx) {
        double begin = cuts[idx - 1], end = cuts[idx];
        const double middle = 0.5 * (begin + end);
        const Affine lower = z_low.at(middle) > 0. ? z_low : zero;
        const Affine upper = z_high.at(middle) < 1. ? z_high : one;
        if (lower.at(middle) > upper.at(middle))
            continue;
        const Affine closest = nearest.at(middle) < lower.at(middle) ? lower :
                               nearest.at(middle) > upper.at(middle) ? upper : nearest;
        const Vec3d origin = a - c - s * closest.origin;
        const Vec3d direction = r - s * closest.slope;
        const double aa = direction.head<2>().squaredNorm();
        const double bb = 2. * origin.head<2>().dot(direction.head<2>());
        const double cc = origin.head<2>().squaredNorm() - xy_tolerance * xy_tolerance;
        if (aa <= 1e-16) {
            if (cc > 0.)
                continue;
        } else {
            const double discriminant = bb * bb - 4. * aa * cc;
            if (discriminant < 0.)
                continue;
            const double root = std::sqrt(discriminant);
            begin = std::max(begin, (-bb - root) / (2. * aa));
            end = std::min(end, (-bb + root) / (2. * aa));
        }
        if (!sloped_support) {
            if (std::abs(r.z()) <= 1e-16) {
                if (std::abs(a.z() - c.z()) > z_tolerance)
                    continue;
            } else {
                const double t0 = (c.z() - a.z() - z_tolerance) / r.z();
                const double t1 = (c.z() - a.z() + z_tolerance) / r.z();
                begin = std::max(begin, std::min(t0, t1));
                end = std::min(end, std::max(t0, t1));
            }
        }
        if (end > begin + 1e-12) {
            if (!intervals.empty() && begin <= intervals.back().second + 1e-12)
                intervals.back().second = end;
            else
                intervals.emplace_back(begin, end);
        }
    }
    return intervals;
}

namespace {

bool make_nonplanar_path(LayerRegion &region, const sla::IndexedMesh &mesh,
                         const std::vector<uint8_t> &selected_facets, ExtrusionPath &path,
                         const ExPolygons *allowed_area = nullptr,
                         double maximum_drape_height_override = 0.)
{
    const bool is_outer_perimeter = path.role() == erExternalPerimeter ||
                                    (is_perimeter(path.role()) && path.inset_idx == 0);
    const bool is_wall = is_perimeter(path.role());
    if ((path.role() != erTopSolidInfill && path.role() != erIroning &&
         !is_wall) ||
        path.polyline.points.size() < 2)
        return false;

    const PrintRegionConfig &config = region.region().config();
    const Layer             &layer  = *region.layer();
    const double resolution = std::max(0.05, config.nonplanar_top_surface_resolution.value);
    // Start just inside the material. A ray originating exactly on the slice
    // plane may numerically miss the same boundary triangle for one adjacent
    // raster endpoint and hit it for the next, creating an artificial Z cliff.
    const double source_z   = layer.slice_z + mesh.ground_level() - 0.01;
    // This path is emitted with the current layer, before any later layer can
    // support it. Therefore it may follow the exposed surface only through the
    // thickness already supplied by the current shell layer. Using the number
    // of top shell layers here allowed an early perimeter to climb several
    // millimetres onto a future part of the surface and print in mid-air.
    // Multi-layer non-planar surfaces are formed by adjacent supported shell
    // cutouts; no individual cutout may skip the material below it.
    const double maximum_drape_height = maximum_drape_height_override > 0. ?
        maximum_drape_height_override :
        std::max(double(path.height), layer.height);

    struct SurfaceSample {
        Vec2d  xy;
        double target_z {0.};
        int    face {-1};
        bool   upward {false};
    };

    auto is_allowed = [allowed_area](const Vec2d &xy) {
        return allowed_area == nullptr || Geometry::contains(*allowed_area,
            Point(scale_(xy.x()), scale_(xy.y())));
    };

    auto query_surface = [&mesh, &selected_facets, &layer, source_z,
                          &is_allowed,
                          is_outer_perimeter, probe_distance =
                              std::max(0.05, 0.5 * double(path.width))]
                         (const Vec2d &xy, const Vec2d &lateral) {
        SurfaceSample sample;
        sample.xy = xy;

        // Most paths on most layers are remote from every selected surface.
        // Reject those samples in 2D before entering IndexedMesh's AABB tree.
        // The domain is expanded when constructed, so boundary walls and
        // their half-bead lateral probes remain eligible.
        if (!is_allowed(xy))
            return sample;

        auto query_at = [&](const Vec2d &probe) {
            SurfaceSample candidate;
            candidate.xy = xy;
            const sla::IndexedMesh::hit_result hit = mesh.query_ray_hit(
                Vec3d(probe.x(), probe.y(), source_z), Vec3d::UnitZ());
            if (hit.is_hit()) {
                candidate.face     = hit.face();
                candidate.upward   = size_t(hit.face()) < selected_facets.size() &&
                                     selected_facets[size_t(hit.face())] != 0;
                candidate.target_z = hit.position().z() - mesh.ground_level();
            }
            return candidate;
        };

        sample = query_at(xy);
        // An external perimeter lies on the intersection between the shallow
        // upper face and its steep side wall. A ray cast exactly through that
        // shared edge may consistently return the side-wall triangle, making
        // perimeter-only surfaces such as Benchy's bow, window rims and small
        // posts invisible to the non-planar discovery pass. Probe half a bead
        // to either side of the centerline and retain the selected upper face,
        // while keeping the emitted XY point on the dimensional contour. Both
        // sides are tested because loop winding is not stable after clipping.
        if (is_outer_perimeter && !sample.upward && lateral.squaredNorm() > EPSILON) {
            const Vec2d offset = lateral.normalized() * probe_distance;
            const std::array<SurfaceSample, 2> candidates = {
                query_at(xy + offset), query_at(xy - offset)};
            for (const SurfaceSample &candidate : candidates) {
                if (!candidate.upward)
                    continue;
                if (!sample.upward ||
                    std::abs(candidate.target_z - layer.print_z) <
                        std::abs(sample.target_z - layer.print_z))
                    sample = candidate;
            }
        }
        sample.xy = xy;
        return sample;
    };

    auto evaluate_surface = [&](const SurfaceSample &surface) {
        double z_offset = 0.;
        bool   sample_is_eligible = false;
        if (surface.face >= 0 && surface.upward && is_allowed(surface.xy)) {
            const double target_z = surface.target_z;
            {
                // Top-surface source regions are already disjoint XY cutouts.
                // Draping only the samples nearest to this nominal layer chops
                // one raster into layer-height fragments and recreates a seam
                // at every Z boundary. Project the complete source path onto
                // the connected surface instead. Its XYZ extrusion may cross
                // several nominal layers; reachability is validated later
                // against the actual toolhead geometry.
                const double height_delta = target_z - layer.print_z;
                if (height_delta <= maximum_drape_height + EPSILON &&
                    height_delta >= -double(path.height) - EPSILON) {
                    z_offset = height_delta;
                    sample_is_eligible = true;
                }
            }
        }
        return std::pair<double, bool>{sample_is_eligible ? z_offset : 0., sample_is_eligible};
    };

    Pointf3s          sampled;
    std::vector<char> eligible;

    const Points3 &points = path.polyline.points;
    for (size_t segment_idx = 0; segment_idx + 1 < points.size(); ++segment_idx) {
        const Vec2d a = unscale(points[segment_idx].to_point());
        const Vec2d b = unscale(points[segment_idx + 1].to_point());
        const Vec2d delta = b - a;
        const Vec2d lateral(-delta.y(), delta.x());
        const int segment_count = std::max(1, int(std::ceil(delta.norm() / resolution)));
        const int first_sample = segment_idx == 0 ? 0 : 1;

        std::vector<SurfaceSample> surface_samples(size_t(segment_count) + 1);
        surface_samples.front() = query_surface(a, lateral);
        surface_samples.back()  = query_surface(b, lateral);

        // A triangle and its XY projection are convex. If two samples hit the
        // same face, every requested sample between them lies on that face and
        // its height is exactly linear. Preserve the configured output point
        // spacing while avoiding an AABB-tree raycast for every 0.2 mm point.
        auto populate_surface_samples = [&](auto &&self, int first, int last) -> void {
            if (last - first <= 1)
                return;
            const SurfaceSample &lhs = surface_samples[size_t(first)];
            const SurfaceSample &rhs = surface_samples[size_t(last)];
            if (lhs.face >= 0 && lhs.face == rhs.face) {
                for (int idx = first + 1; idx < last; ++idx) {
                    const double ratio = double(idx - first) / double(last - first);
                    SurfaceSample &sample = surface_samples[size_t(idx)];
                    sample.xy       = a + delta * (double(idx) / double(segment_count));
                    sample.target_z = lhs.target_z + (rhs.target_z - lhs.target_z) * ratio;
                    sample.face     = lhs.face;
                    sample.upward   = lhs.upward;
                }
                return;
            }
            const int middle = first + (last - first) / 2;
            surface_samples[size_t(middle)] = query_surface(
                a + delta * (double(middle) / double(segment_count)), lateral);
            self(self, first, middle);
            self(self, middle, last);
        };
        populate_surface_samples(populate_surface_samples, 0, segment_count);

        for (int sample_idx = first_sample; sample_idx <= segment_count; ++sample_idx) {
            const SurfaceSample &surface = surface_samples[size_t(sample_idx)];
            const auto [z_offset, sample_is_eligible] = evaluate_surface(surface);
            sampled.emplace_back(surface.xy.x(), surface.xy.y(), z_offset);
            eligible.push_back(sample_is_eligible);
        }
    }

    if (sampled.size() < 2)
        return false;

    // Fill generators extend centerlines to the edge of their clipped area.
    // The endpoint may consequently sit a few scaled units outside the mesh
    // and miss the raycast even though the bead is part of the same sloped top
    // surface. Continue the nearest valid height across those end samples;
    // do not bridge internal misses, which may represent real unsafe geometry.
    const auto first_eligible = std::find(eligible.begin(), eligible.end(), true);
    const auto last_eligible  = std::find(eligible.rbegin(), eligible.rend(), true);
    if (first_eligible != eligible.end() && last_eligible != eligible.rend()) {
        const size_t first_idx = size_t(std::distance(eligible.begin(), first_eligible));
        const size_t last_idx  = eligible.size() - 1 -
            size_t(std::distance(eligible.rbegin(), last_eligible));
        for (size_t idx = 0; idx < first_idx && is_allowed(sampled[idx].head<2>()); ++idx) {
            sampled[idx].z() = sampled[first_idx].z();
            eligible[idx] = true;
        }
        for (size_t idx = last_idx + 1; idx < eligible.size() &&
             is_allowed(sampled[idx].head<2>()); ++idx) {
            sampled[idx].z() = sampled[last_idx].z();
            eligible[idx] = true;
        }
    }

    // An isolated mesh hit cannot establish a surface direction. Flatten it
    // rather than creating a small Z spike between planar samples.
    for (size_t idx = 0; idx < eligible.size(); ++idx) {
        const bool has_eligible_neighbor =
            (idx > 0 && eligible[idx - 1]) || (idx + 1 < eligible.size() && eligible[idx + 1]);
        if (eligible[idx] && !has_eligible_neighbor) {
            eligible[idx] = false;
            sampled[idx].z() = 0.;
        }
    }

    // A selected perimeter remains one slicer-owned wall. Blend the surface
    // displacement into the nominal planar course over the same distance
    // profile at both boundaries. Measured from the boundary, this profile is
    // direction independent, so planar->non-planar is exactly the reverse of
    // non-planar->planar without maintaining two algorithms.
    //
    // A completely selected closed loop has no eligibility boundary. Its seam
    // is the transition boundary: start and finish at the planar course and
    // blend onto the projected surface between them. Without this case the
    // wall was topologically closed but still approached as a separate raised
    // operation, which produced the visible recess and weak attachment seen on
    // small rings and bosses.
    const bool has_planar_samples = std::any_of(
        eligible.begin(), eligible.end(), [](char value) { return value == 0; });
    const bool has_nonplanar_samples = std::any_of(
        eligible.begin(), eligible.end(), [](char value) { return value != 0; });
    const bool closed_wall = is_wall && sampled.size() > 2 &&
        (sampled.front().head<2>() - sampled.back().head<2>()).norm() <= EPSILON;
    if ((has_planar_samples && has_nonplanar_samples) ||
        (closed_wall && !has_planar_samples)) {
        const double transition_length = std::max(2., 4. * double(path.width));
        std::vector<double> before(eligible.size(), 0.);
        std::vector<double> after(eligible.size(), 0.);
        const bool seam_is_transition = closed_wall && !has_planar_samples;
        if (eligible.front() && !seam_is_transition)
            before.front() = transition_length;
        for (size_t idx = 1; idx < eligible.size(); ++idx)
            if (eligible[idx] && eligible[idx - 1])
                before[idx] = std::min(
                    transition_length,
                    before[idx - 1] +
                        (sampled[idx].head<2>() -
                         sampled[idx - 1].head<2>()).norm());
        if (eligible.back() && !seam_is_transition)
            after.back() = transition_length;
        for (size_t idx = eligible.size() - 1; idx > 0; --idx)
            if (eligible[idx - 1] && eligible[idx])
                after[idx - 1] = std::min(
                    transition_length,
                    after[idx] +
                        (sampled[idx].head<2>() -
                         sampled[idx - 1].head<2>()).norm());
        for (size_t idx = 0; idx < eligible.size(); ++idx) {
            if (!eligible[idx])
                continue;
            const double t = std::clamp(
                std::min(before[idx], after[idx]) / transition_length, 0., 1.);
            const double smooth = t * t * (3. - 2. * t);
            sampled[idx].z() *= smooth;
        }
    }

    // Do not impose a process-level surface-angle ceiling here. G-code emission
    // validates every variable-Z segment against the printer's detailed
    // toolhead geometry or its fallback clearance cone and flattens paths that
    // cannot be reached safely.

    const bool adjusted = std::any_of(sampled.begin(), sampled.end(),
        [](const Vec3d &point) { return std::abs(point.z()) > EPSILON; });
    if (!adjusted)
        return false;

    Polyline3 polyline;
    polyline.points.reserve(sampled.size());
    for (const Vec3d &point : sampled)
        polyline.append(Point3(scale_(point.x()), scale_(point.y()), scale_(point.z())));

    // Ray sampling determines where the mesh surface changes; it does not mean
    // every sample must become a G-code vertex. A triangle is planar, so most
    // consecutive samples are collinear in XYZ. Removing points whose 3D
    // deviation is below five microns preserves the sampled surface while
    // avoiding thousands of redundant moves and clearance checks per layer.
    polyline.simplify(scale_(std::min(0.005, resolution * 0.025)));
    if (polyline.points.size() < 2)
        return false;

    path.polyline          = std::move(polyline);
    path.z_contoured       = true;
    path.nonplanar_surface = true;
    path.nonplanar_schedule_owned = false;
    // A wall may be one directed part of a loop or multipath. Reversing that
    // part independently disconnects it from both neighbors even though its
    // unordered XY endpoints still appear to form a closed graph. Only free
    // raster paths may choose their own high-to-low direction here.
    if (!is_wall) {
        if (path.polyline.points.front().z() < path.polyline.points.back().z())
            path.reverse();
        path.set_reverse();
    }
    return true;
}

std::optional<double> nonplanar_height(const ExtrusionEntity &entity)
{
    Points3 points;
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity)) {
        if (!path->nonplanar_surface)
            return std::nullopt;
        points = path->polyline.points;
    } else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
        for (const ExtrusionPath &path : multipath->paths) {
            if (!path.nonplanar_surface)
                continue;
            append(points, path.polyline.points);
        }
    } else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        std::optional<double> maximum;
        for (const ExtrusionEntity *child : collection->entities)
            if (const std::optional<double> height = nonplanar_height(*child)) {
                maximum = maximum ? std::max(*maximum, *height) : height;
            }
        return maximum;
    }
    if (points.empty())
        return std::nullopt;
    return unscale_(std::max_element(points.begin(), points.end(),
        [](const Point3 &lhs, const Point3 &rhs) { return lhs.z() < rhs.z(); })->z());
}

std::vector<ExtrusionPath> split_planar_and_nonplanar(const ExtrusionPath &path)
{
    if (!path.nonplanar_surface || path.polyline.points.size() < 2)
        return {path};

    std::vector<ExtrusionPath> result;
    Points3 chunk{path.polyline.points.front()};
    const coord_t z_epsilon = scale_(EPSILON);
    bool chunk_is_nonplanar = std::abs(path.polyline.points.front().z()) > z_epsilon;

    for (size_t idx = 1; idx < path.polyline.points.size(); ++idx) {
        const Point3 &point = path.polyline.points[idx];
        const bool segment_is_nonplanar =
            std::abs(path.polyline.points[idx - 1].z()) > z_epsilon ||
            std::abs(point.z()) > z_epsilon;
        if (segment_is_nonplanar != chunk_is_nonplanar && chunk.size() >= 2) {
            ExtrusionPath split = path;
            split.polyline.points = std::move(chunk);
            split.z_contoured = chunk_is_nonplanar;
            split.nonplanar_surface = chunk_is_nonplanar;
            result.emplace_back(std::move(split));
            chunk = {path.polyline.points[idx - 1]};
        }
        chunk_is_nonplanar = segment_is_nonplanar;
        chunk.push_back(point);
    }

    if (chunk.size() >= 2) {
        ExtrusionPath split = path;
        split.polyline.points = std::move(chunk);
        split.z_contoured = chunk_is_nonplanar;
        split.nonplanar_surface = chunk_is_nonplanar;
        result.emplace_back(std::move(split));
    }
    return result;
}

void restore_candidate_path_as_planar(ExtrusionPath &path)
{
    for (Point3 &point : path.polyline.points)
        point.z() = 0;
    path.polyline.fitting_result.clear();
    path.z_contoured = false;
    path.nonplanar_surface = false;
    path.nonplanar_schedule_owned = false;
    path.nonplanar_transition = false;
    path.nonplanar_clearance_validated = false;
}

bool projected_wall_loop_is_closed(const ExtrusionLoop &loop)
{
    constexpr double join_tolerance_mm = 0.05;
    const double join_tolerance_squared =
        std::pow(double(scale_(join_tolerance_mm)), 2.);

    if (loop.paths.empty())
        return false;
    for (size_t i = 0; i < loop.paths.size(); ++i) {
        const ExtrusionPath &path = loop.paths[i];
        const ExtrusionPath &next = loop.paths[(i + 1) % loop.paths.size()];
        if (path.polyline.points.size() < 2 || next.polyline.points.size() < 2 ||
            (path.last_point3() - next.first_point3()).cast<double>().squaredNorm() >
                join_tolerance_squared)
            return false;
    }
    return true;
}

void process_entity(LayerRegion &region, const sla::IndexedMesh &mesh,
                    const std::vector<uint8_t> &selected_facets, ExtrusionEntity &entity,
                    const ExPolygons *allowed_area = nullptr)
{
    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
        make_nonplanar_path(region, mesh, selected_facets, *path, allowed_area);
    } else if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity)) {
        const ExtrusionPaths original = multipath->paths;
        for (ExtrusionPath &path : multipath->paths) {
            path.inset_idx = multipath->inset_idx;
            make_nonplanar_path(region, mesh, selected_facets, path, allowed_area);
        }
        if (nonplanar_height(*multipath)) {
            const bool all_wall_paths_nonplanar = !multipath->paths.empty() &&
                std::all_of(multipath->paths.begin(), multipath->paths.end(),
                    [](const ExtrusionPath &path) {
                        return is_perimeter(path.role()) && path.nonplanar_surface;
                    });
            // A fully projected open multipath is still an isolated wall
            // operation. Mixed multipaths retain their planar entry/exit and
            // are the intended in-place transition representation.
            if (all_wall_paths_nonplanar &&
                multipath->first_point() != multipath->last_point()) {
                multipath->paths = original;
                return;
            }
            if (multipath->first_point3().z() > multipath->last_point3().z())
                multipath->reverse();
            multipath->set_reverse();
        }
    } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
        const ExtrusionPaths original = loop->paths;
        for (ExtrusionPath &path : loop->paths) {
            path.inset_idx = loop->inset_idx;
            make_nonplanar_path(region, mesh, selected_facets, path, allowed_area);
        }
        // nonplanar_height() describes sortable free paths, not loops.
        if (std::any_of(loop->paths.begin(), loop->paths.end(),
                [](const ExtrusionPath &path) { return path.nonplanar_surface; }) &&
            !projected_wall_loop_is_closed(*loop))
            loop->paths = original;
    } else if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
        ExtrusionEntitiesPtr rebuilt;
        rebuilt.reserve(collection->entities.size());
        for (ExtrusionEntity *child : collection->entities) {
            if (auto *path = dynamic_cast<ExtrusionPath *>(child)) {
                std::optional<ExtrusionPath> original;
                if (is_perimeter(path->role()))
                    original = *path;
                make_nonplanar_path(region, mesh, selected_facets, *path, allowed_area);
                std::vector<ExtrusionPath> split =
                    split_planar_and_nonplanar(*path);
                // A bare open perimeter path has no owning loop which can
                // connect it to the rest of the wall graph.  Projecting the
                // whole fragment therefore creates a detached raised line
                // which is approached with a retract/travel (Arachne may
                // produce these for narrow residual wall regions).  Mixed
                // paths are safe: the planar pieces and both transitions are
                // retained below in one ExtrusionMultiPath.  Fully projected
                // walls are safe only when their original entity is closed;
                // open fragments remain conventional until a graph-level
                // neighbor can own both endpoints.
                if (is_perimeter(path->role()) && split.size() == 1 &&
                    split.front().nonplanar_surface &&
                    split.front().first_point() != split.front().last_point()) {
                    *path = *original;
                    split = {*path};
                }
                if (is_perimeter(path->role()) && split.size() > 1) {
                    auto *connected = new ExtrusionMultiPath(std::move(split));
                    connected->inset_idx = path->inset_idx;
                    connected->set_reverse();
                    rebuilt.emplace_back(connected);
                } else {
                    for (ExtrusionPath &part : split)
                        rebuilt.emplace_back(part.clone_move());
                }
                delete child;
            } else {
                process_entity(region, mesh, selected_facets, *child, allowed_area);
                rebuilt.emplace_back(child);
            }
        }
        collection->entities = std::move(rebuilt);
        std::stable_sort(collection->entities.begin(), collection->entities.end(),
            [](const ExtrusionEntity *lhs, const ExtrusionEntity *rhs) {
                const std::optional<double> lhs_height = nonplanar_height(*lhs);
                const std::optional<double> rhs_height = nonplanar_height(*rhs);
                if (!lhs_height)
                    return bool(rhs_height);
                if (!rhs_height)
                    return false;
                // Print the high side first, then descend. Printing low paths
                // first makes the nozzle body sweep through those beads while
                // reaching the neighboring high raster.
                return *lhs_height > *rhs_height;
            });
        if (std::any_of(collection->entities.begin(), collection->entities.end(),
                [](const ExtrusionEntity *child) { return bool(nonplanar_height(*child)); }))
            collection->no_sort = true;
    }
}

void contour_existing_perimeter_graph(
    LayerRegion &region, const sla::IndexedMesh &mesh,
    const std::vector<uint8_t> &selected_facets, ExtrusionEntity &entity,
    const ExPolygons &allowed_area, double maximum_drape_height,
    std::vector<ExtrusionPath *> &adjusted_paths)
{
    const auto contour_path = [&](ExtrusionPath &path) {
        if (!is_perimeter(path.role()) || path.nonplanar_surface ||
            path.nonplanar_clearance_validated)
            return;
        if (make_nonplanar_path(region, mesh, selected_facets, path,
                                &allowed_area, maximum_drape_height))
            adjusted_paths.emplace_back(&path);
    };

    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
        contour_path(*path);
    } else if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity)) {
        for (ExtrusionPath &path : multipath->paths) {
            path.inset_idx = multipath->inset_idx;
            contour_path(path);
        }
    } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
        const size_t first_adjusted = adjusted_paths.size();
        const ExtrusionPaths original = loop->paths;
        for (ExtrusionPath &path : loop->paths) {
            path.inset_idx = loop->inset_idx;
            contour_path(path);
        }
        // Preserve the perimeter generator's loop topology. A partial mesh
        // projection may modify several constituent paths, but it may never
        // turn their connected loop into detached raised fragments.
        if (adjusted_paths.size() != first_adjusted &&
            !projected_wall_loop_is_closed(*loop)) {
            adjusted_paths.resize(first_adjusted);
            loop->paths = original;
        }
    } else if (auto *collection =
                   dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
        for (ExtrusionEntity *child : collection->entities)
            contour_existing_perimeter_graph(
                region, mesh, selected_facets, *child, allowed_area,
                maximum_drape_height, adjusted_paths);
    }
}

void collect_nonplanar_top_coverage(const ExtrusionEntity &entity, Polygons &coverage)
{
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity)) {
        // A shallow component is sometimes represented only by its outer wall
        // at a particular slice (small rims, posts and the narrowing bow are
        // common examples). Its bead footprint is sufficient to discover the
        // connected selected surface patch; the patch projection, not this
        // narrow seed, defines the eventual replacement skin. Do not use the
        // whole closed-loop interior here because that can seed unrelated flat
        // or steep geometry merely enclosed by the same perimeter.
        const bool is_outer_outline = is_perimeter(path->role()) && path->inset_idx == 0;
        if (path->nonplanar_surface &&
            (path->role() == erTopSolidInfill || is_outer_outline))
            path->polygons_covered_by_spacing(coverage, float(scale_(0.02)));
    } else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
        for (const ExtrusionPath &path : multipath->paths)
            collect_nonplanar_top_coverage(path, coverage);
    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
        for (const ExtrusionPath &path : loop->paths)
            collect_nonplanar_top_coverage(path, coverage);
    } else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        for (const ExtrusionEntity *child : collection->entities)
            collect_nonplanar_top_coverage(*child, coverage);
    }
}

void collect_nonplanar_outline_areas(const ExtrusionEntity &entity, Polygons &areas)
{
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity)) {
        if (path->nonplanar_surface && is_perimeter(path->role()) &&
            path->inset_idx == 0 && !path->polyline.points.empty() &&
            path->polyline.points.front() == path->polyline.points.back()) {
            Polygon polygon;
            polygon.points = to_points(path->polyline.points);
            if (!polygon.points.empty())
                polygon.points.pop_back();
            if (polygon.is_valid())
                areas.emplace_back(std::move(polygon));
        }
    } else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
        for (const ExtrusionPath &path : multipath->paths)
            collect_nonplanar_outline_areas(path, areas);
    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
        const bool owns_nonplanar_outline = loop->inset_idx == 0 &&
            std::any_of(loop->paths.begin(), loop->paths.end(), [](const ExtrusionPath &path) {
                return path.nonplanar_surface && is_perimeter(path.role());
            });
        if (owns_nonplanar_outline) {
            Polygon polygon = loop->polygon();
            if (polygon.is_valid())
                areas.emplace_back(std::move(polygon));
        }
    } else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        for (const ExtrusionEntity *child : collection->entities)
            collect_nonplanar_outline_areas(*child, areas);
    }
}

void append_without_nonplanar_top_overlap(const ExtrusionEntity &entity,
                                          const ExPolygons &top_coverage,
                                          const BoundingBox &coverage_bounds,
                                          ExtrusionEntitiesPtr &out,
                                          bool bridges_and_sparse_only = false,
                                          bool preserve_accepted_nonplanar = true,
                                          const std::function<ExPolygons(const ExtrusionPath &)> *coverage_at_path = nullptr,
                                          const std::map<const ExtrusionPath *, ExtrusionPaths> *path_replacements = nullptr)
{
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity)) {
        if (path_replacements != nullptr) {
            const auto replacement = path_replacements->find(path);
            if (replacement != path_replacements->end()) {
                for (const ExtrusionPath &fragment : replacement->second)
                    out.emplace_back(fragment.clone());
                return;
            }
        }
        if (coverage_at_path != nullptr) {
            const ExPolygons coverage = (*coverage_at_path)(*path);
            if (coverage.empty())
                out.emplace_back(path->clone());
            else
                append_without_nonplanar_top_overlap(entity, coverage,
                    get_extents(coverage), out, bridges_and_sparse_only,
                    preserve_accepted_nonplanar);
            return;
        }
        if (preserve_accepted_nonplanar && path->nonplanar_surface &&
            path->nonplanar_clearance_validated) {
            // This is an already accepted surface from an earlier, lower run.
            // A later candidate may inspect the same nominal layer while
            // transitioning another feature; never carve the completed skin.
            out.emplace_back(path->clone());
        } else if (bridges_and_sparse_only &&
                   path->role() != erInternalInfill &&
                   !is_bridge_or_arc_fill(path->role())) {
            out.emplace_back(path->clone());
        } else if (path->polyline.points.empty() ||
                   !path->polyline.to_polyline().bounding_box().overlap(
                       coverage_bounds)) {
            // Most replacement footprints are local while a layer may contain
            // hundreds of unrelated wall and infill paths. Clipper still
            // builds and sweeps a complete edge graph when two inputs have
            // disjoint bounds. Reject those pairs here; the exact subtraction
            // below remains authoritative for the small overlapping subset.
            out.emplace_back(path->clone());
        } else {
            // The replacement footprint owns every model extrusion in this
            // height band, not only top fill and walls. A sloping skin often
            // intersects paths classified as internal bridge or sparse fill
            // on an earlier nominal layer. Retaining those paths produces
            // same-Z crossings with the draped skin and makes an otherwise
            // valid surface unprintable. This helper is only called for an
            // object's fills/perimeters and with the shell-specific coverage,
            // so clipping all conventional roles here does not affect support
            // or material outside the replaced volume. Provisional projected
            // outlines are clipped too: the accepted Arachne wall stack now
            // replaces them, and retaining an old outer segment produced the
            // same-Z crossings seen on small curved Bunny surfaces.
            ExtrusionEntityCollection clipped;
            path->subtract_expolygons(top_coverage, &clipped);
            for (auto fragment = clipped.entities.begin();
                 fragment != clipped.entities.end();) {
                if (auto *fragment_path = dynamic_cast<ExtrusionPath *>(*fragment)) {
                    fragment_path->nonplanar_replacement_remainder = true;
                    // Exact polygon subtraction may leave alternating
                    // sub-bead tips of the old planar top raster along the
                    // accepted skin boundary. They are neither a printable
                    // planar surface nor useful buried support: hiding the
                    // replacement exposes them as floating dashes. Buried
                    // solid/sparse fill and perimeter fragments are retained
                    // because they may be the validated foundation or lateral
                    // anchor for the replacement.
                    if (fragment_path->role() == erTopSolidInfill &&
                        unscale<double>(fragment_path->length()) + EPSILON <
                            std::max(0.05, double(fragment_path->width))) {
                        delete *fragment;
                        fragment = clipped.entities.erase(fragment);
                        continue;
                    }
                }
                ++fragment;
            }
            out.insert(out.end(),
                       std::make_move_iterator(clipped.entities.begin()),
                       std::make_move_iterator(clipped.entities.end()));
            clipped.entities.clear();
        }
    } else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
        const bool owns_accepted_surface = !multipath->paths.empty() && std::all_of(
            multipath->paths.begin(), multipath->paths.end(),
            [](const ExtrusionPath &path) {
                return path.nonplanar_surface &&
                       path.nonplanar_clearance_validated;
            });
        // Variable-flow fill has the same continuity contract as walls.
        // Keeping each child unchanged but putting it into a collection
        // turns thickness samples into separate extrusion/retraction events.
        if (preserve_accepted_nonplanar && owns_accepted_surface)
            out.emplace_back(multipath->clone());
        else {
            auto *rebuilt = new ExtrusionEntityCollection;
            rebuilt->no_sort = true;
            rebuilt->inset_idx = multipath->inset_idx;
            for (const ExtrusionPath &path : multipath->paths)
                append_without_nonplanar_top_overlap(
                    path, top_coverage, coverage_bounds, rebuilt->entities,
                    bridges_and_sparse_only, preserve_accepted_nonplanar, coverage_at_path, path_replacements);
            // Clipping may change only part of a connected path. Reassemble
            // consecutive retained pieces at exact XYZ endpoints; only an
            // actual removed interval should introduce a new print move.
            ExtrusionEntitiesPtr connected;
            ExtrusionPaths chain;
            const auto flush = [&]() {
                if (chain.empty())
                    return;
                if (chain.size() == 1)
                    connected.emplace_back(new ExtrusionPath(std::move(chain.front())));
                else {
                    auto *joined = new ExtrusionMultiPath(std::move(chain));
                    joined->inset_idx = multipath->inset_idx;
                    connected.emplace_back(joined);
                }
                chain.clear();
            };
            for (ExtrusionEntity *child : rebuilt->entities) {
                auto *path = dynamic_cast<ExtrusionPath *>(child);
                if (path == nullptr || path->polyline.points.size() < 2) {
                    flush();
                    connected.emplace_back(child);
                    continue;
                }
                if (!chain.empty() && chain.back().last_point3() != path->first_point3())
                    flush();
                chain.emplace_back(std::move(*path));
                delete child;
            }
            flush();
            rebuilt->entities = std::move(connected);
            if (rebuilt->empty())
                delete rebuilt;
            else
                out.emplace_back(rebuilt);
        }
    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
        const bool owns_accepted_wall = !loop->paths.empty() && std::all_of(
            loop->paths.begin(), loop->paths.end(),
            [](const ExtrusionPath &path) {
                return path.nonplanar_surface &&
                       path.nonplanar_clearance_validated &&
                       is_perimeter(path.role());
            });
        if (preserve_accepted_nonplanar && owns_accepted_wall)
            out.emplace_back(loop->clone());
        else {
            // Removing a section of a closed perimeter necessarily makes the
            // remainder open. Keep those pieces in one unsortable owner with
            // the original inset metadata: flattening them into the parent
            // lets nearest-path scheduling silently override wall_sequence.
            auto *rebuilt = new ExtrusionEntityCollection;
            rebuilt->no_sort = true;
            rebuilt->inset_idx = loop->inset_idx;
            for (const ExtrusionPath &path : loop->paths)
                append_without_nonplanar_top_overlap(
                    path, top_coverage, coverage_bounds, rebuilt->entities,
                    bridges_and_sparse_only, preserve_accepted_nonplanar, coverage_at_path, path_replacements);
            // A layer-wide clipping pass also visits unaffected loops. Do
            // not turn those into independently scheduled paths: retaining
            // every original segment means the original closed wall and its
            // seam policy still apply.
            const bool unchanged = rebuilt->entities.size() == loop->paths.size() &&
                std::equal(rebuilt->entities.begin(), rebuilt->entities.end(), loop->paths.begin(),
                    [](const ExtrusionEntity *entity, const ExtrusionPath &original) {
                        const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
                        return path != nullptr && path->polyline.points == original.polyline.points;
                    });
            if (unchanged) {
                delete rebuilt;
                out.emplace_back(loop->clone());
            } else if (rebuilt->empty())
                delete rebuilt;
            else
                out.emplace_back(rebuilt);
        }
    } else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        auto *rebuilt = new ExtrusionEntityCollection;
        rebuilt->no_sort = collection->no_sort;
        rebuilt->inset_idx = collection->inset_idx;
        for (const ExtrusionEntity *child : collection->entities)
            append_without_nonplanar_top_overlap(
                *child, top_coverage, coverage_bounds, rebuilt->entities,
                bridges_and_sparse_only, preserve_accepted_nonplanar, coverage_at_path, path_replacements);
        if (rebuilt->empty())
            delete rebuilt;
        else
            out.emplace_back(rebuilt);
    }
}

void remove_nonplanar_top_overlap(ExtrusionEntityCollection &collection,
                                  const ExPolygons &top_coverage,
                                  bool bridges_and_sparse_only = false,
                                  bool preserve_accepted_nonplanar = true,
                                  const std::function<ExPolygons(const ExtrusionPath &)> *coverage_at_path = nullptr,
                                  const std::map<const ExtrusionPath *, ExtrusionPaths> *path_replacements = nullptr)
{
    if (top_coverage.empty() && coverage_at_path == nullptr)
        return;
    const BoundingBox coverage_bounds = get_extents(top_coverage);
    ExtrusionEntitiesPtr rebuilt;
    rebuilt.reserve(collection.entities.size());
    for (const ExtrusionEntity *entity : collection.entities)
        append_without_nonplanar_top_overlap(
            *entity, top_coverage, coverage_bounds, rebuilt,
            bridges_and_sparse_only,
            preserve_accepted_nonplanar, coverage_at_path, path_replacements);
    collection.clear();
    collection.entities = std::move(rebuilt);
}

bool prune_empty_extrusion_paths(ExtrusionEntity &entity)
{
    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity))
        return path->polyline.points.size() >= 2;
    if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity)) {
        multipath->paths.erase(
            std::remove_if(multipath->paths.begin(), multipath->paths.end(),
                [](const ExtrusionPath &path) {
                    return path.polyline.points.size() < 2;
                }),
            multipath->paths.end());
        return !multipath->paths.empty();
    }
    if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
        loop->paths.erase(
            std::remove_if(loop->paths.begin(), loop->paths.end(),
                [](const ExtrusionPath &path) {
                    return path.polyline.points.size() < 2;
                }),
            loop->paths.end());
        return !loop->paths.empty();
    }
    if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
        auto &entities = collection->entities;
        for (auto it = entities.begin(); it != entities.end();) {
            if (prune_empty_extrusion_paths(**it)) {
                ++it;
            } else {
                delete *it;
                it = entities.erase(it);
            }
        }
        return !collection->empty();
    }
    return true;
}

int perimeter_inset(const ExtrusionEntity &entity)
{
    if (entity.inset_idx >= 0)
        return entity.inset_idx;
    if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity))
        return multipath->paths.empty() ? -1 : int(multipath->paths.front().inset_idx);
    if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity))
        return loop->paths.empty() ? -1 : int(loop->paths.front().inset_idx);
    if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity))
        for (const ExtrusionEntity *child : collection->entities)
            if (const int inset = perimeter_inset(*child); inset >= 0)
                return inset;
    return -1;
}

void restore_nonplanar_wall_sequence(ExtrusionEntityCollection &collection,
                                     WallSequence wall_sequence,
                                     int wall_loops,
                                     bool above_first_layer)
{
    for (ExtrusionEntity *entity : collection.entities)
        if (auto *nested = dynamic_cast<ExtrusionEntityCollection *>(entity))
            restore_nonplanar_wall_sequence(
                *nested, wall_sequence, wall_loops, above_first_layer);
    const bool sandwich = wall_sequence == WallSequence::InnerOuterInner &&
                          wall_loops >= 3 && above_first_layer;
    const auto key = [wall_sequence, sandwich](const ExtrusionEntity *entity) {
        const int inset = perimeter_inset(*entity);
        if (inset < 0)
            return std::pair<int, int>{1, 0};
        if (sandwich) {
            if (inset >= 2)
                return std::pair<int, int>{0, -inset};
            return std::pair<int, int>{0, inset == 0 ? 0 : 1};
        }
        const bool outside_first = wall_sequence == WallSequence::OuterInner ||
                                   wall_sequence == WallSequence::InnerOuterInner;
        return std::pair<int, int>{0, outside_first ? inset : -inset};
    };
    std::stable_sort(collection.entities.begin(), collection.entities.end(),
        [&](const ExtrusionEntity *left, const ExtrusionEntity *right) {
            return key(left) < key(right);
        });
    collection.no_sort = true;
}

bool paths_have_proper_same_z_crossing(
    const ExtrusionPath &first, double first_owner_z,
    const ExtrusionPath &second, double second_owner_z,
    double z_tolerance, double xy_tolerance = 0.05)
{
    // The writer rounds XY/Z to 0.001 mm. On steep, nearly tangent paths that
    // quantization changes the interpolated height at the intersection by
    // several microns. Keep final graph arbitration outside the emitted
    // 0.05 mm collision band instead of stopping exactly on its boundary.
    const double effective_z_tolerance = z_tolerance + 0.01;
    const Points3 &first_points = first.polyline.points;
    const Points3 &second_points = second.polyline.points;
    // The test is symmetric. Build and sort the broad-phase index for the
    // shorter path; feature courses are often compared against a much longer
    // roof raster, and rebuilding that larger index for every wall dominated
    // non-planar validation time.
    if (second_points.size() > first_points.size())
        return paths_have_proper_same_z_crossing(
            second, second_owner_z, first, first_owner_z, z_tolerance,
            xy_tolerance);
    struct SegmentBounds {
        size_t index;
        double min_x;
        double max_x;
        double min_y;
        double max_y;
        double start_z;
        double end_z;
    };
    std::vector<SegmentBounds> second_segments;
    second_segments.reserve(second_points.size() > 0 ?
                                second_points.size() - 1 : 0);
    for (size_t idx = 1; idx < second_points.size(); ++idx) {
        const Vec2d start = unscale(second_points[idx - 1].to_point());
        const Vec2d end = unscale(second_points[idx].to_point());
        second_segments.push_back({
            idx, std::min(start.x(), end.x()), std::max(start.x(), end.x()),
            std::min(start.y(), end.y()), std::max(start.y(), end.y()),
            second_owner_z + unscale_(second_points[idx - 1].z()),
            second_owner_z + unscale_(second_points[idx].z())});
    }
    std::sort(second_segments.begin(), second_segments.end(),
        [](const SegmentBounds &left, const SegmentBounds &right) {
            return left.min_x < right.min_x;
        });
    for (size_t first_idx = 1; first_idx < first_points.size(); ++first_idx) {
        const Vec2d a = unscale(first_points[first_idx - 1].to_point());
        const Vec2d b = unscale(first_points[first_idx].to_point());
        const Vec2d first_delta = b - a;
        const double first_min_x = std::min(a.x(), b.x());
        const double first_max_x = std::max(a.x(), b.x());
        const double first_min_y = std::min(a.y(), b.y());
        const double first_max_y = std::max(a.y(), b.y());
        const double first_start_z = first_owner_z +
            unscale_(first_points[first_idx - 1].z());
        const double first_end_z = first_owner_z +
            unscale_(first_points[first_idx].z());
        const double first_min_z = std::min(first_start_z, first_end_z);
        const double first_max_z = std::max(first_start_z, first_end_z);
        for (const SegmentBounds &bounds : second_segments) {
            // G-code coordinates are rounded and the deposited bead has a
            // finite width. Keep the broad phase consistent with the exact
            // near-collinear test below: centerlines separated by only a few
            // microns may become the same emitted line and must not escape
            // collision arbitration at this early bounds check.
            if (bounds.min_x > first_max_x + xy_tolerance)
                break;
            if (bounds.max_x < first_min_x - xy_tolerance ||
                bounds.min_y > first_max_y + xy_tolerance ||
                bounds.max_y < first_min_y - xy_tolerance)
                continue;
            const double second_start_z = bounds.start_z;
            const double second_end_z = bounds.end_z;
            if (std::max(first_min_z,
                         std::min(second_start_z, second_end_z)) >
                    std::min(first_max_z,
                             std::max(second_start_z, second_end_z)) +
                        effective_z_tolerance)
                continue;
            const size_t second_idx = bounds.index;
            const Vec2d c = unscale(second_points[second_idx - 1].to_point());
            const Vec2d d = unscale(second_points[second_idx].to_point());
            const Vec2d second_delta = d - c;
            const double denominator =
                first_delta.x() * second_delta.y() -
                first_delta.y() * second_delta.x();
            const double first_length = first_delta.norm();
            const double second_length = second_delta.norm();
            if (first_length <= EPSILON || second_length <= EPSILON)
                continue;
            if (std::abs(denominator) <=
                    1e-4 * first_length * second_length) {
                const double offset_cross =
                    first_delta.x() * (c.y() - a.y()) -
                    first_delta.y() * (c.x() - a.x());
                if (std::abs(offset_cross) > xy_tolerance * first_length)
                    continue;
                const double first_length_squared = first_delta.squaredNorm();
                const double second_length_squared = second_delta.squaredNorm();
                const double second_start_on_first =
                    (c - a).dot(first_delta) / first_length_squared;
                const double second_end_on_first =
                    (d - a).dot(first_delta) / first_length_squared;
                // Preserve genuine endpoint bonding, but leave enough guard
                // for the 0.001 mm XY quantization performed by the G-code
                // writer. A crossing at roughly one percent of a short
                // segment may move from just outside to just inside after
                // rounding, so final graph arbitration has to be slightly
                // more conservative than the emitted-path audit.
                constexpr double endpoint_guard = 0.005;
                const double overlap_start = std::max(
                    endpoint_guard,
                    std::min(second_start_on_first, second_end_on_first));
                const double overlap_end = std::min(
                    1. - endpoint_guard,
                    std::max(second_start_on_first, second_end_on_first));
                if ((overlap_end - overlap_start) * first_length <= 0.05)
                    continue;
                const auto z_difference = [&](double first_t) {
                    const Vec2d position = a + first_t * first_delta;
                    const double second_t =
                        (position - c).dot(second_delta) /
                        second_length_squared;
                    if (second_t <= endpoint_guard ||
                        second_t >= 1. - endpoint_guard)
                        return std::optional<double>{};
                    const double first_z = first_start_z +
                        first_t * (first_end_z - first_start_z);
                    const double second_z = second_start_z +
                        second_t * (second_end_z - second_start_z);
                    return std::optional<double>{first_z - second_z};
                };
                const std::optional<double> start_difference =
                    z_difference(overlap_start);
                const std::optional<double> end_difference =
                    z_difference(overlap_end);
                // Along collinear XY spans both Z profiles are linear. Their
                // closest approach is therefore at an overlap endpoint or at
                // the zero crossing of their difference. Sampling only the
                // midpoint misses the common ramp-vs-nearly-flat case where
                // the beads cross closer to one end of the shared span.
                if ((start_difference &&
                     std::abs(*start_difference) < effective_z_tolerance) ||
                    (end_difference &&
                     std::abs(*end_difference) < effective_z_tolerance) ||
                    (start_difference && end_difference &&
                     *start_difference * *end_difference <= 0.))
                    return true;
                continue;
            }
            const Vec2d offset = c - a;
            const double first_t =
                (offset.x() * second_delta.y() -
                 offset.y() * second_delta.x()) / denominator;
            const double second_t =
                (offset.x() * first_delta.y() -
                 offset.y() * first_delta.x()) / denominator;
            // A fill endpoint meeting a wall is the intended structural bond.
            // Only a centerline which enters and leaves the wall is a proper
            // crossing and therefore a duplicate/collision candidate.
            constexpr double endpoint_guard = 0.005;
            if (first_t <= endpoint_guard ||
                first_t >= 1. - endpoint_guard ||
                second_t <= endpoint_guard ||
                second_t >= 1. - endpoint_guard)
                continue;
            const double first_z = first_start_z +
                first_t * (first_end_z - first_start_z);
            const double second_z = second_start_z +
                second_t * (second_end_z - second_start_z);
            if (std::abs(first_z - second_z) < effective_z_tolerance)
                return true;
        }
    }
    return false;
}

bool bounds_overlap_with_tolerance(const BoundingBox &first,
                                   const BoundingBox &second,
                                   coord_t tolerance)
{
    return first.min.x() <= second.max.x() + tolerance &&
           second.min.x() <= first.max.x() + tolerance &&
           first.min.y() <= second.max.y() + tolerance &&
           second.min.y() <= first.max.y() + tolerance;
}

bool path_has_proper_same_z_self_crossing(const ExtrusionPath &path,
                                          double owner_z,
                                          double z_tolerance)
{
    const Points3 &points = path.polyline.points;
    if (points.size() < 4)
        return false;
    const bool closed = points.front() == points.back();
    const coord_t xy_tolerance = scale_(0.02);
    const coord_t scaled_z_tolerance = scale_(z_tolerance + 0.01);
    for (size_t first_idx = 1; first_idx < points.size(); ++first_idx) {
        const Point3 &first_start = points[first_idx - 1];
        const Point3 &first_end = points[first_idx];
        const coord_t first_min_x = std::min(first_start.x(), first_end.x());
        const coord_t first_max_x = std::max(first_start.x(), first_end.x());
        const coord_t first_min_y = std::min(first_start.y(), first_end.y());
        const coord_t first_max_y = std::max(first_start.y(), first_end.y());
        const coord_t first_min_z = std::min(first_start.z(), first_end.z());
        const coord_t first_max_z = std::max(first_start.z(), first_end.z());
        for (size_t second_idx = first_idx + 2;
             second_idx < points.size(); ++second_idx) {
            // The first and last segments of a closed loop share its seam and
            // are adjacent, not a crossing.  All other non-adjacent segment
            // intersections represent a path which enters its own deposited
            // bead and must be rejected before G-code ordering can hide it.
            if (closed && first_idx == 1 && second_idx + 1 == points.size())
                continue;
            const Point3 &second_start = points[second_idx - 1];
            const Point3 &second_end = points[second_idx];
            // Most pairs in a long raster or perimeter are remote. Avoid
            // constructing two ExtrusionPath objects and rebuilding the
            // generic segment index unless their emitted XY and Z envelopes
            // can actually touch. This check uses the same tolerances as the
            // exact routine below and therefore changes only cost, not the
            // collision decision.
            if (std::min(second_start.x(), second_end.x()) >
                    first_max_x + xy_tolerance ||
                std::max(second_start.x(), second_end.x()) <
                    first_min_x - xy_tolerance ||
                std::min(second_start.y(), second_end.y()) >
                    first_max_y + xy_tolerance ||
                std::max(second_start.y(), second_end.y()) <
                    first_min_y - xy_tolerance ||
                std::min(second_start.z(), second_end.z()) >
                    first_max_z + scaled_z_tolerance ||
                std::max(second_start.z(), second_end.z()) <
                    first_min_z - scaled_z_tolerance)
                continue;
            ExtrusionPath first(path);
            first.polyline.points = {first_start, first_end};
            ExtrusionPath second(path);
            second.polyline.points = {second_start, second_end};
            if (paths_have_proper_same_z_crossing(
                    first, owner_z, second, owner_z, z_tolerance))
                return true;
        }
    }
    return false;
}

// Join open, surface-following wall fragments to the planar wall graph they
// replace.  Non-planar top skin is generated after conventional perimeters,
// but an open boundary fragment is not a second shell: it is an edge of the
// same wall.  Keeping the two entity lists independent leaves a retract/travel
// at each boundary and, more importantly, creates two closed or weakly
// touching contours.  This routine consumes the clipped planar fragments and
// emits connected ExtrusionMultiPaths in the generated wall order.
//
// Return false rather than committing a partially connected graph.  The
// caller retains an unmodified copy of the layer and atomically restores the
// conventional slice in that case.
[[maybe_unused]] bool integrate_nonplanar_perimeters(
    ExtrusionEntityCollection &perimeters, ExtrusionEntitiesPtr &generated)
{
    struct Segment {
        ExtrusionPath path;
        bool nonplanar {false};
        size_t source_index {size_t(-1)};
    };
    struct CandidateJoin {
        size_t nonplanar_endpoint {0};
        size_t planar_endpoint {0};
        double distance_squared {0.};
    };

    ExtrusionEntityCollection flattened = perimeters.flatten(false);
    std::vector<Segment> segments;
    std::vector<size_t> planar_segment_by_entity(flattened.entities.size(), size_t(-1));
    std::vector<size_t> nonplanar_segment_by_entity(generated.size(), size_t(-1));

    for (size_t entity_idx = 0; entity_idx < flattened.entities.size(); ++entity_idx) {
        const auto *path = dynamic_cast<const ExtrusionPath *>(flattened.entities[entity_idx]);
        if (path == nullptr || !is_perimeter(path->role()) ||
            path->polyline.points.size() < 2)
            continue;
        planar_segment_by_entity[entity_idx] = segments.size();
        segments.push_back(Segment{*path, false, entity_idx});
    }
    const size_t planar_count = segments.size();
    for (size_t entity_idx = 0; entity_idx < generated.size(); ++entity_idx) {
        const auto *path = dynamic_cast<const ExtrusionPath *>(generated[entity_idx]);
        if (path == nullptr || !path->nonplanar_surface ||
            !is_perimeter(path->role()) || path->polyline.points.size() < 2 ||
            path->first_point() == path->last_point())
            continue;
        nonplanar_segment_by_entity[entity_idx] = segments.size();
        segments.push_back(Segment{*path, true, entity_idx});
    }
    if (segments.size() == planar_count)
        return true;

    const size_t endpoint_count = 2 * segments.size();
    std::vector<size_t> joined_to(endpoint_count, size_t(-1));
    std::vector<CandidateJoin> candidates;
    for (size_t nonplanar_idx = planar_count;
         nonplanar_idx < segments.size(); ++nonplanar_idx) {
        const Segment &nonplanar = segments[nonplanar_idx];
        const double join_distance = std::clamp(
            1.75 * double(nonplanar.path.width), 0.35, 0.9);
        const double maximum_distance_squared =
            double(scale_(join_distance)) * double(scale_(join_distance));
        for (size_t planar_idx = 0; planar_idx < planar_count; ++planar_idx) {
            const Segment &planar = segments[planar_idx];
            if (planar.path.role() != nonplanar.path.role() ||
                planar.path.inset_idx != nonplanar.path.inset_idx)
                continue;
            for (size_t nonplanar_end = 0; nonplanar_end < 2; ++nonplanar_end) {
                const Point np = nonplanar_end == 0 ?
                    nonplanar.path.first_point() : nonplanar.path.last_point();
                for (size_t planar_end = 0; planar_end < 2; ++planar_end) {
                    const Point pp = planar_end == 0 ?
                        planar.path.first_point() : planar.path.last_point();
                    const double distance_squared = (np - pp).squaredNorm();
                    if (distance_squared <= maximum_distance_squared)
                        candidates.push_back(CandidateJoin{
                            2 * nonplanar_idx + nonplanar_end,
                            2 * planar_idx + planar_end,
                            distance_squared});
                }
            }
        }
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const CandidateJoin &left, const CandidateJoin &right) {
                  return left.distance_squared < right.distance_squared;
              });
    for (const CandidateJoin &candidate : candidates) {
        if (joined_to[candidate.nonplanar_endpoint] != size_t(-1) ||
            joined_to[candidate.planar_endpoint] != size_t(-1))
            continue;
        joined_to[candidate.nonplanar_endpoint] = candidate.planar_endpoint;
        joined_to[candidate.planar_endpoint] = candidate.nonplanar_endpoint;
    }
    for (size_t segment_idx = planar_count; segment_idx < segments.size(); ++segment_idx)
        if (joined_to[2 * segment_idx] == size_t(-1) ||
            joined_to[2 * segment_idx + 1] == size_t(-1)) {
            std::array<double, 2> nearest{
                std::numeric_limits<double>::max(),
                std::numeric_limits<double>::max()};
            for (size_t end = 0; end < 2; ++end) {
                const Point endpoint = end == 0 ?
                    segments[segment_idx].path.first_point() :
                    segments[segment_idx].path.last_point();
                for (size_t planar_idx = 0; planar_idx < planar_count; ++planar_idx) {
                    if (segments[planar_idx].path.role() !=
                            segments[segment_idx].path.role() ||
                        segments[planar_idx].path.inset_idx !=
                            segments[segment_idx].path.inset_idx)
                        continue;
                    nearest[end] = std::min({
                        nearest[end],
                        unscale<double>((endpoint -
                            segments[planar_idx].path.first_point())
                                .cast<double>().norm()),
                        unscale<double>((endpoint -
                            segments[planar_idx].path.last_point())
                                .cast<double>().norm())});
                }
            }
            BOOST_LOG_TRIVIAL(debug)
                << "Rejecting nonplanar wall graph: open wall inset "
                << segments[segment_idx].path.inset_idx
                << " has an endpoint without a planar-wall match; nearest="
                << nearest[0] << ',' << nearest[1]
                << " mm, planar candidates=" << planar_count;
            return false;
        }

    std::vector<uint8_t> visited(segments.size(), 0);
    std::vector<uint8_t> consumed_planar(flattened.entities.size(), 0);
    for (size_t seed = planar_count; seed < segments.size(); ++seed) {
        if (visited[seed])
            continue;

        // Discover the complete alternating planar/non-planar component.
        std::vector<size_t> component;
        std::vector<size_t> pending{seed};
        while (!pending.empty()) {
            const size_t segment_idx = pending.back();
            pending.pop_back();
            if (std::find(component.begin(), component.end(), segment_idx) != component.end())
                continue;
            component.push_back(segment_idx);
            for (size_t end = 0; end < 2; ++end) {
                const size_t neighbor_endpoint = joined_to[2 * segment_idx + end];
                if (neighbor_endpoint != size_t(-1))
                    pending.push_back(neighbor_endpoint / 2);
            }
        }

        size_t start_segment = component.front();
        size_t start_end = 0;
        bool found_open_start = false;
        for (size_t segment_idx : component)
            for (size_t end = 0; end < 2; ++end)
                if (joined_to[2 * segment_idx + end] == size_t(-1)) {
                    start_segment = segment_idx;
                    start_end = end;
                    found_open_start = true;
                    break;
                }
        if (!found_open_start) {
            const auto planar_it = std::find_if(
                component.begin(), component.end(),
                [planar_count](size_t idx) { return idx < planar_count; });
            if (planar_it != component.end())
                start_segment = *planar_it;
        }

        auto *combined = new ExtrusionMultiPath;
        size_t segment_idx = start_segment;
        size_t enter_end = start_end;
        while (!visited[segment_idx]) {
            visited[segment_idx] = 1;
            ExtrusionPath oriented = segments[segment_idx].path;
            if (enter_end == 1)
                oriented.reverse();
            combined->paths.emplace_back(std::move(oriented));
            if (!segments[segment_idx].nonplanar)
                consumed_planar[segments[segment_idx].source_index] = 1;

            const size_t exit_endpoint = 2 * segment_idx + (1 - enter_end);
            const size_t neighbor_endpoint = joined_to[exit_endpoint];
            if (neighbor_endpoint == size_t(-1))
                break;
            const size_t neighbor_segment = neighbor_endpoint / 2;
            const size_t neighbor_end = neighbor_endpoint % 2;
            const Point3 &from = combined->paths.back().last_point3();
            const Point3 &to = neighbor_end == 0 ?
                segments[neighbor_segment].path.first_point3() :
                segments[neighbor_segment].path.last_point3();
            if (from != to) {
                ExtrusionPath connector = segments[neighbor_segment].path;
                connector.polyline.points = {from, to};
                connector.polyline.fitting_result.clear();
                connector.nonplanar_surface = true;
                connector.nonplanar_schedule_owned = false;
                connector.nonplanar_transition = true;
                connector.nonplanar_clearance_validated = true;
                combined->paths.emplace_back(std::move(connector));
            }
            segment_idx = neighbor_segment;
            enter_end = neighbor_end;
        }
        combined->inset_idx = segments[seed].path.inset_idx;
        combined->set_reverse();

        const auto first_nonplanar = std::find_if(
            component.begin(), component.end(),
            [planar_count](size_t idx) { return idx >= planar_count; });
        assert(first_nonplanar != component.end());
        const size_t output_idx = segments[*first_nonplanar].source_index;
        delete generated[output_idx];
        generated[output_idx] = combined;
        for (size_t component_idx : component)
            if (component_idx >= planar_count) {
                const size_t source_idx = segments[component_idx].source_index;
                if (source_idx != output_idx) {
                    delete generated[source_idx];
                    generated[source_idx] = nullptr;
                }
            }
    }

    ExtrusionEntitiesPtr remaining_planar;
    remaining_planar.reserve(flattened.entities.size());
    for (size_t idx = 0; idx < flattened.entities.size(); ++idx) {
        if (!consumed_planar[idx])
            remaining_planar.emplace_back(flattened.entities[idx]->clone());
    }
    perimeters.clear();
    if (!remaining_planar.empty()) {
        // LayerRegion::perimeters stores top-level island collections; several
        // downstream consumers intentionally cast each child to a collection.
        // flatten() above is only an integration workspace, so restore that
        // container contract instead of leaking leaf paths into the top level.
        auto *remaining = new ExtrusionEntityCollection;
        remaining->no_sort = flattened.no_sort;
        remaining->entities = std::move(remaining_planar);
        perimeters.entities.emplace_back(remaining);
    }
    perimeters.no_sort = flattened.no_sort;

    ExtrusionEntitiesPtr remaining_generated;
    remaining_generated.reserve(generated.size());
    for (size_t idx = 0; idx < generated.size(); ++idx) {
        if (generated[idx] != nullptr)
            remaining_generated.emplace_back(generated[idx]);
    }
    generated = std::move(remaining_generated);
    return true;
}

} // namespace

void Layer::make_nonplanar_top_surfaces(
    const sla::IndexedMesh &mesh,
    const std::vector<std::vector<uint8_t>> &selected_facets,
    const std::vector<std::vector<ExPolygons>> &selected_projections)
{
    // Mark exposed top skin and its outermost contour as candidate seeds. The
    // path projector accepts an outer wall only when the selected upper face
    // is within one physical layer of that wall, so vertical side walls below
    // an overhanging roof remain planar. Including the contour is necessary
    // for perimeter-only details (small rims, holes, wheels, windows, posts)
    // that have no top-solid raster from which to discover their surface.
    for (LayerRegion *region : this->regions()) {
        if (!nonplanar_perimeters_enabled(region->region().config()))
            continue;
        const int region_id = region->region().print_region_id();
        if (region_id < 0 || size_t(region_id) >= selected_facets.size() ||
            selected_facets[size_t(region_id)].empty())
            continue;
        const ExPolygons *allowed_area =
            size_t(region_id) < selected_projections.size() &&
            this->id() < selected_projections[size_t(region_id)].size() ?
                &selected_projections[size_t(region_id)][this->id()] : nullptr;
        process_entity(*region, mesh, selected_facets[size_t(region_id)], region->fills,
                       allowed_area);
        process_entity(*region, mesh, selected_facets[size_t(region_id)], region->perimeters,
                       allowed_area);
    }

    // Candidate marking must be non-destructive. Earlier code immediately
    // subtracted the provisional footprint from top-solid fill and expected
    // consolidation to restore it later. A rejected patch could reset flags
    // on surviving paths, but it could not resurrect paths already deleted by
    // that subtraction, leaving internal-solid infill exposed as the visible
    // top layer. Consolidation below first restores these intact candidates,
    // then carves only the bead footprint of a complete patch which passed
    // support and toolhead-clearance validation.
}

namespace {

bool entity_has_nonplanar_surface(const ExtrusionEntity &entity)
{
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity))
        return path->nonplanar_surface && !path->nonplanar_feature_transition;
    if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity))
        return std::any_of(multipath->paths.begin(), multipath->paths.end(),
            [](const ExtrusionPath &path) {
                return path.nonplanar_surface && !path.nonplanar_feature_transition;
            });
    if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity))
        return std::any_of(loop->paths.begin(), loop->paths.end(),
            [](const ExtrusionPath &path) {
                return path.nonplanar_surface && !path.nonplanar_feature_transition;
            });
    if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity))
        return std::any_of(collection->entities.begin(), collection->entities.end(),
            [](const ExtrusionEntity *child) { return entity_has_nonplanar_surface(*child); });
    return false;
}

bool entity_is_entirely_nonplanar(const ExtrusionEntity &entity)
{
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity))
        return path->nonplanar_surface && !path->nonplanar_feature_transition;
    if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity))
        return !multipath->paths.empty() &&
            std::all_of(multipath->paths.begin(), multipath->paths.end(),
                [](const ExtrusionPath &path) {
                    return path.nonplanar_surface && !path.nonplanar_feature_transition;
                });
    if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity))
        return !loop->paths.empty() &&
            std::all_of(loop->paths.begin(), loop->paths.end(),
                [](const ExtrusionPath &path) {
                    return path.nonplanar_surface && !path.nonplanar_feature_transition;
                });
    if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity))
        return !collection->entities.empty() &&
            std::all_of(collection->entities.begin(), collection->entities.end(),
                [](const ExtrusionEntity *child) { return entity_is_entirely_nonplanar(*child); });
    return false;
}

ExtrusionPath *clone_path_with_metadata(const ExtrusionPath &path)
{
    auto *copy = new ExtrusionPath(path);
    copy->inset_idx = path.inset_idx;
    return copy;
}

void split_nonplanar_entity(const ExtrusionEntity &entity,
                            ExtrusionEntitiesPtr &planar,
                            ExtrusionEntitiesPtr &nonplanar)
{
    if (!entity_has_nonplanar_surface(entity)) {
        ExtrusionEntity *copy = entity.clone();
        copy->inset_idx = entity.inset_idx;
        planar.emplace_back(copy);
        return;
    }
    if (entity_is_entirely_nonplanar(entity)) {
        ExtrusionEntity *copy = entity.clone();
        copy->inset_idx = entity.inset_idx;
        nonplanar.emplace_back(copy);
        return;
    }

    if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
        for (const ExtrusionPath &path : multipath->paths)
            (path.nonplanar_surface ? nonplanar : planar).emplace_back(
                clone_path_with_metadata(path));
    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
        // A mixed loop is already open at the planar/non-planar boundary. Keep
        // its pieces as paths instead of falsely reconnecting them at one Z.
        for (const ExtrusionPath &path : loop->paths)
            (path.nonplanar_surface ? nonplanar : planar).emplace_back(
                clone_path_with_metadata(path));
    } else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        auto *remaining = new ExtrusionEntityCollection;
        remaining->no_sort = collection->no_sort;
        remaining->inset_idx = collection->inset_idx;
        for (const ExtrusionEntity *child : collection->entities)
            split_nonplanar_entity(*child, remaining->entities, nonplanar);
        if (remaining->empty())
            delete remaining;
        else
            planar.emplace_back(remaining);
    }
}

[[maybe_unused]] void extract_nonplanar_entities(
    ExtrusionEntityCollection &collection, ExtrusionEntitiesPtr &nonplanar)
{
    ExtrusionEntitiesPtr remaining;
    remaining.reserve(collection.entities.size());
    for (const ExtrusionEntity *entity : collection.entities)
        split_nonplanar_entity(*entity, remaining, nonplanar);
    collection.clear();
    collection.entities = std::move(remaining);
}

void offset_entity_z(ExtrusionEntity &entity, coord_t delta)
{
    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
        for (Point3 &point : path->polyline.points)
            point.z() += delta;
    } else if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity)) {
        for (ExtrusionPath &path : multipath->paths)
            offset_entity_z(path, delta);
    } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
        for (ExtrusionPath &path : loop->paths)
            offset_entity_z(path, delta);
    } else if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
        for (ExtrusionEntity *child : collection->entities)
            offset_entity_z(*child, delta);
    }
}

void set_nonplanar_transition(ExtrusionEntity &entity, bool transition)
{
    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
        if (path->nonplanar_surface)
            path->nonplanar_transition = transition;
    } else if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity)) {
        for (ExtrusionPath &path : multipath->paths)
            set_nonplanar_transition(path, transition);
    } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
        for (ExtrusionPath &path : loop->paths)
            set_nonplanar_transition(path, transition);
    } else if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
        for (ExtrusionEntity *child : collection->entities)
            set_nonplanar_transition(*child, transition);
    }
}

void make_planar_solid_foundation(ExtrusionEntity &entity,
                                  const ExPolygons &foundation_material)
{
    const auto make_path = [&foundation_material](ExtrusionPath &path) {
        for (Point3 &point : path.polyline.points)
            point.z() = 0;
        if (is_perimeter(path.role())) {
            // The final-surface outline is morphed toward this real planar
            // course, but a rounded corner may still finish a few hundredths
            // beyond the supporting bead envelope. Clamp only unsupported
            // points to that envelope; already valid dimensional centerlines
            // remain untouched.
            const ExPolygons bead_support = offset_ex(
                foundation_material,
                scale_(0.5 * double(path.width) + 0.02));
            if (!bead_support.empty())
                for (Point3 &point : path.polyline.points)
                    if (!Geometry::contains(
                            bead_support, point.to_point())) {
                        const Point supported = projection_onto(
                            bead_support, point.to_point());
                        point.x() = supported.x();
                        point.y() = supported.y();
                    }
            // Course zero is the last fully horizontal perimeter below the
            // staircase. Preserve its structural role and inset, but keep it
            // conventional. Calling a fixed-Z copy a non-planar transition
            // made the G-code look like a post-processed shell and caused it
            // to be scheduled independently from the planar foundation it
            // replaces. The following aligned course is the first true
            // transition and carries the variable-Z metadata.
            path.z_contoured = false;
            // Keep scheduler ownership so this foundation is emitted before
            // the first raised course. It is deliberately not Z-contoured or
            // a transition, and G-code/preview therefore expose its original
            // structural perimeter role. Dropping surface ownership here
            // moved the foundation back into ordinary island ordering, after
            // the transition it was meant to support.
            path.nonplanar_surface = true;
            path.nonplanar_schedule_owned = true;
            path.nonplanar_transition = false;
            path.nonplanar_clearance_validated = false;
            path.polyline.fitting_result.clear();
            return;
        }
        path.set_extrusion_role(erSolidInfill);
        path.z_contoured = false;
        path.nonplanar_surface = false;
        path.nonplanar_schedule_owned = false;
        path.nonplanar_transition = false;
        path.nonplanar_clearance_validated = false;
        path.polyline.fitting_result.clear();
    };
    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
        make_path(*path);
    } else if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity)) {
        for (ExtrusionPath &path : multipath->paths)
            make_path(path);
    } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
        for (ExtrusionPath &path : loop->paths)
            make_path(path);
    } else if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
        for (ExtrusionEntity *child : collection->entities)
            make_planar_solid_foundation(*child, foundation_material);
    }
}

void set_nonplanar_local_course_height(ExtrusionPath &path,
                                       double local_height)
{
    if (path.height <= EPSILON || path.width <= EPSILON)
        return;
    const double old_area = Flow(path.width, path.height, 0.f).mm3_per_mm();
    const double new_area = Flow(
        path.width, float(local_height), 0.f).mm3_per_mm();
    if (old_area > EPSILON)
        path.mm3_per_mm *= new_area / old_area;
    path.height = float(local_height);
}

void adjust_nonplanar_stack_flow(
    ExtrusionEntity &entity, double source_owner_z, double planar_support_z,
    double blend_step, double minimum_height, double maximum_height)
{
    const auto adjust_path = [&](ExtrusionPath &path) {
        if (path.polyline.points.size() < 2)
            return;
        double weighted_span = 0.;
        double total_length = 0.;
        for (size_t idx = 1; idx < path.polyline.points.size(); ++idx) {
            const Point3 &first = path.polyline.points[idx - 1];
            const Point3 &last = path.polyline.points[idx];
            const double length = unscale_(
                (last.to_point() - first.to_point()).cast<double>().norm());
            if (length <= EPSILON)
                continue;
            const double first_surface_z =
                source_owner_z + unscale_(first.z());
            const double last_surface_z =
                source_owner_z + unscale_(last.z());
            weighted_span += length * 0.5 *
                (std::max(0., first_surface_z - planar_support_z) +
                 std::max(0., last_surface_z - planar_support_z));
            total_length += length;
        }
        if (total_length <= EPSILON)
            return;
        const double local_height = std::clamp(
            blend_step * weighted_span / total_length,
            minimum_height, maximum_height);
        set_nonplanar_local_course_height(path, local_height);
    };
    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
        adjust_path(*path);
    } else if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity)) {
        for (ExtrusionPath &path : multipath->paths)
            adjust_path(path);
    } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
        for (ExtrusionPath &path : loop->paths)
            adjust_path(path);
    } else if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
        for (ExtrusionEntity *child : collection->entities)
            adjust_nonplanar_stack_flow(
                *child, source_owner_z, planar_support_z, blend_step,
                minimum_height, maximum_height);
    }
}

void flatten_interlocking_wall_foundation(
    ExtrusionEntity &entity, const ExPolygons &surface_area)
{
    const auto flatten_path = [&surface_area](ExtrusionPath &path) {
        if (!path.nonplanar_interlocking_wall ||
            std::none_of(path.polyline.points.begin(),
                         path.polyline.points.end(),
                         [&surface_area](const Point3 &point) {
                             return Geometry::contains(
                                 surface_area, point.to_point());
                         }))
            return;
        // Buried sine-wave walls deliberately dip below their nominal plane.
        // The course immediately below a mesh-following surface is instead
        // its dimensional foundation: retaining a trough here can increase
        // the first surface-course gap by the full interlock amplitude. End
        // the buried interlock one course early and hand the unchanged XY
        // perimeter graph to the flat-to-non-planar transition above it.
        for (Point3 &point : path.polyline.points)
            point.z() = 0;
        path.nonplanar_interlocking_wall = false;
        path.z_contoured = false;
        path.polyline.fitting_result.clear();
    };
    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity))
        flatten_path(*path);
    else if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity))
        for (ExtrusionPath &path : multipath->paths)
            flatten_path(path);
    else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity))
        for (ExtrusionPath &path : loop->paths)
            flatten_path(path);
    else if (auto *collection =
                 dynamic_cast<ExtrusionEntityCollection *>(&entity))
        for (ExtrusionEntity *child : collection->entities)
            flatten_interlocking_wall_foundation(*child, surface_area);
}

size_t planar_support_layer_index(const std::vector<Layer *> &layers,
                                  double requested_z)
{
    // Use the first printable course at or above the requested bottom of the
    // transition. Rounding down adds up to one complete planar layer to the
    // height that the requested number of non-planar courses must span. The
    // resulting course-to-course distance can then exceed the configured
    // layer height, which appears as a floating finishing shell in preview
    // and produces weak or absent bonding on the printer. Rounding upward
    // keeps the complete blend inside the requested shell-depth budget. XY
    // anchoring below still projects the first course onto this real slice,
    // so the transition begins on printable planar geometry.
    auto lower = std::lower_bound(
        layers.begin(), layers.end(), requested_z - EPSILON,
        [](const Layer *layer, double z) { return layer->print_z < z; });
    return lower == layers.end() ? layers.size() - 1 :
        size_t(std::distance(layers.begin(), lower));
}

double planar_support_layer_z(const std::vector<Layer *> &layers,
                              double requested_z)
{
    return layers[planar_support_layer_index(layers, requested_z)]->print_z;
}

size_t planar_transition_anchor_layer_index(
    const std::vector<Layer *> &layers, double minimum_surface_z,
    double transition_depth)
{
    // A transition cannot eliminate stair stepping when every XY sample is
    // anchored to surface_z - shell_depth: that height field is the original
    // staircase, merely copied below the finishing surface.  Use one real,
    // horizontal slice for the complete connected patch instead.  The last
    // horizontal course below the requested transition depth is the final
    // ordinary course before this surface begins adopting any lateral step.
    // Round upward so the first generated course remains within one nominal
    // layer of its retained foundation.
    auto first_surface_layer = std::lower_bound(
        layers.begin(), layers.end(),
        minimum_surface_z - transition_depth - EPSILON,
        [](const Layer *layer, double z) { return layer->print_z < z; });
    return first_surface_layer == layers.end() ? layers.size() - 1 :
        size_t(std::distance(layers.begin(), first_surface_layer));
}

struct NonplanarPatch;
void destroy_entities(ExtrusionEntitiesPtr &entities);
std::optional<double> patch_surface_z_at_xy(
    const NonplanarPatch &patch, const sla::IndexedMesh &mesh,
    const Point &point);
std::optional<double> patch_surface_z_near_xy(
    const NonplanarPatch &patch, const sla::IndexedMesh &mesh,
    const Point &point, double maximum_distance);
size_t owner_layer_for_z(const std::vector<Layer *> &layers, size_t first_layer,
                         size_t last_layer, double required_z);

bool blend_wall_z_from_planar_support(
    ExtrusionPath &path, const std::vector<Layer *> &layers,
    const NonplanarPatch &patch, const sla::IndexedMesh &mesh,
    coordf_t destination_owner_z, size_t transition_anchor_layer,
    double planar_layer_height, double minimum_course_height,
    double blend_step, double surface_reference_z, double surface_depth,
    double blend)
{
    // Preserve the exact XY centerline, inset, width, role, brick offset and
    // perimeter ordering produced by the ordinary slicer. Only the height
    // field is introduced gradually. Consequently every course is the same
    // continuous perimeter it would have been in a planar slice, rather than
    // a separately generated concentric shell.
    // A connected perimeter may also contain an ordinary planar portion
    // outside this surface. Keep that portion on its owner course, while the
    // contoured part uses exactly the same support-to-surface interpolation
    // as the raster stack. Scaling the old relative Z about the nominal owner
    // made adjacent wall and raster courses describe different height fields
    // and cross at the same physical Z.
    Points3 transitioned = path.polyline.points;
    std::vector<std::optional<double>> sampled_surface_z;
    sampled_surface_z.reserve(transitioned.size());
    for (const Point3 &point : transitioned) {
        const std::optional<double> surface_z = patch_surface_z_near_xy(
            patch, mesh, point.to_point(),
            0.5 * double(path.width) + 0.05);
        sampled_surface_z.emplace_back(surface_z);
    }
    const bool closed_surface_loop = transitioned.size() > 2 &&
        transitioned.front().to_point() == transitioned.back().to_point() &&
        std::all_of(sampled_surface_z.begin(), sampled_surface_z.end(),
                    [](const std::optional<double> &z) { return z.has_value(); });
    // Arachne also represents a wholly selected wall as connected open width
    // segments. If its projected surface happens to cross the owner plane,
    // that zero offset is not a planar attachment. Recognize full projection
    // from the complete path before resolving the ambiguous individual point.
    const bool fully_selected_path = closed_surface_loop || std::equal(
        transitioned.begin(), transitioned.end(), sampled_surface_z.begin(),
        [&](const Point3 &point, const std::optional<double> &surface_z) {
            return surface_z && std::abs(destination_owner_z + unscale_(point.z()) -
                (*surface_z - mesh.ground_level())) <= 1e-5;
        });
    for (size_t point_idx = 0; point_idx < transitioned.size(); ++point_idx) {
        Point3 &point = transitioned[point_idx];
        // The generic contouring pass eases a fully selected closed loop down
        // to its nominal owner Z at the seam so a stand-alone raised loop can
        // be approached.  Accepted surface stacks already have the dense flat
        // foundation and every preceding wall course below them. Repeating
        // that seam dip on every course creates several unsupported points at
        // nominal stair height and breaks what should be one continuous wall.
        // Remove it only when every point belongs to this surface; mixed
        // planar/non-planar loops retain their in-path entry and exit blend.
        if (point.z() == 0 && !sampled_surface_z[point_idx])
            continue;
        const double original_absolute_z =
            destination_owner_z + unscale_(point.z());
        const std::optional<double> &mesh_z = sampled_surface_z[point_idx];
        if (!mesh_z)
            return false;
        const double surface_z = *mesh_z - mesh.ground_level();
        // Every transition course starts at the same real planar foundation
        // and progressively adopts the mesh height field.  The resulting
        // local course height is adjusted separately to match this blend;
        // retaining nominal-height flow here would overfill the low side of
        // a broad roof where the courses are intentionally closer together.
        // Match the retained slicer wall to the same parallel surface course
        // as the generated skin.  Interpolating from a single global anchor Z
        // caused every nominal perimeter between the patch low and high points
        // to become a separately raised closed ring.  `surface_depth` is the
        // authoritative physical offset of this course from the mesh.
        const double transition_plane_z =
            surface_reference_z - surface_depth;
        const double fully_transitioned_z = transition_plane_z +
            blend * (surface_z - surface_reference_z);
        // contour_path_to_surface() already constructed the symmetric
        // planar<->surface entry profile along this exact perimeter. Preserve
        // that profile while moving its fully selected portion onto the
        // patch-wide transition stack. Treating every non-zero sample as
        // fully selected turned the first sample after a planar endpoint into
        // an abrupt Z jump and could exceed the safe cone despite a shallow
        // target surface.
        // A closed loop wholly owned by this patch is already connected to
        // itself and to the preceding wall course.  Its provisional seam dip
        // is an approach aid for a stand-alone contoured loop, not part of
        // the accepted transition stack.  Reapplying that dip independently
        // on each nominal outline makes the slightly different outlines cut
        // across one another near the seam.  Use the authoritative shell
        // blend over the complete closed loop.  Open/mixed fragments still
        // keep the symmetric in-path planar<->surface influence so their
        // endpoints join the untouched planar perimeter without a Z step.
        const double absolute_z = nonplanar_wall_course_z(
            destination_owner_z, original_absolute_z, surface_z,
            fully_transitioned_z, fully_selected_path);
        point.z() = scale_(absolute_z - destination_owner_z);
    }
    path.polyline.points = std::move(transitioned);
    path.polyline.remove_duplicate_points();
    path.polyline.fitting_result.clear();

    // Candidate contouring must not destroy the selected structural wall
    // mode: the patch may still be rejected or only replace a local section.
    // Once this exact wall is accepted into the surface stack, however, it
    // must have one authoritative Z schedule. Restore nominal course flow and
    // remove brick/interlocking offsets only from the accepted segment.
    if (path.staggered_perimeter || path.nonplanar_interlocking_wall ||
        std::abs(double(path.height) - planar_layer_height) > EPSILON) {
        const double old_height = std::max(0.01, double(path.height));
        path.mm3_per_mm *= planar_layer_height / old_height;
        path.height = float(planar_layer_height);
        path.staggered_perimeter = false;
        path.nonplanar_interlocking_wall = false;
    }

    // Parallel courses retain the process layer height.  The XYZ segment
    // length accounts for the extra material needed by sloped motion.
    set_nonplanar_local_course_height(path, planar_layer_height);
    (void) layers;
    (void) transition_anchor_layer;
    (void) minimum_course_height;
    (void) blend_step;
    (void) blend;

    // mm3_per_mm remains the original planar/brick value. G-code emission
    // measures each segment in XYZ, so extrusion rises naturally as the
    // blended Z slope increases without inventing a second flow multiplier or
    // overriding the existing inner-wall flow option.
    return true;
}

coord_t quantized_support_inset(coord_t inset)
{
    const coord_t quantum = scale_(0.02);
    return inset <= 0 ? coord_t(0) :
        coord_t(((inset + quantum - 1) / quantum) * quantum);
}

bool blend_shell_from_planar_support(
    ExtrusionEntity &entity, const std::vector<Layer *> &layers,
    const NonplanarPatch &patch, const sla::IndexedMesh &mesh,
    coordf_t source_owner_z, coordf_t destination_owner_z,
    size_t transition_anchor_layer, double total_shell_depth,
    coord_t fill_support_inset, double surface_reference_z,
    double surface_depth, double blend, bool use_local_support,
    const std::map<std::pair<size_t, coord_t>, ExPolygons> &support_cache,
    const std::map<std::pair<size_t, int>, Polylines> &wall_support_cache)
{
    // Variable-width Arachne paths may differ by only a few microns. Building
    // a separate Clipper offset for every such width made the support cache
    // scale with path count (hundreds of offsets on the Benchy roof). Round
    // outward to a small conservative grid: sharing the next-larger inset can
    // only shrink the accepted support, while reducing the cache to the small
    // number of physically meaningful bead/inset combinations.
    const auto blend_path = [&](ExtrusionPath &path) {
        const bool perimeter = is_perimeter(path.role());
        const Flow path_flow(path.width, path.height, 0.f);
        // A wall in the transition stack must originate from the matching
        // centerline contour of the lower planar shell. Merely checking that
        // its final XY point lies somewhere inside the lower solid leaves all
        // wall clones on the final footprint, producing a vertical stack of
        // self-contained rings. Include the wall inset when constructing the
        // support contour so outer and inner walls retain their separation.
        const coord_t support_inset = quantized_support_inset(perimeter ?
            scale_(std::max(0., 0.5 * double(path.width) - 0.02)) +
                coord_t(std::max(0, path.inset_idx)) *
                    path_flow.scaled_spacing() :
            fill_support_inset);
        for (Point3 &point : path.polyline.points) {
            double surface_z = source_owner_z + unscale_(point.z());
            const Point original_xy = point.to_point();
            // A developable roof must originate from one real, horizontal
            // source course. Selecting a support layer from the
            // final surface height independently at each XY point reproduced
            // the planar staircase inside course zero: a single generated
            // perimeter could be pulled toward several unrelated slice
            // contours and become a set of floating, laterally shifted rings.
            // It also meant the advertised "flat foundation" was not flat in
            // topology even though all of its points shared one Z value.
            // Morph the complete path from the transition anchor selected for
            // this patch.  Subsequent aligned courses release that one lateral
            // correction with `blend`, so they remain bonded while converging
            // continuously on the actual mesh contour. A genuinely doubly
            // curved crown is different: there is no single corresponding
            // lower contour, so its closed tracks retain local parallel-
            // surface support selected from the buried surface offset.
            const size_t support_idx = use_local_support ?
                owner_layer_for_z(
                    layers, transition_anchor_layer, layers.size() - 1,
                    surface_z - total_shell_depth) :
                transition_anchor_layer;

            // A Z-only blend leaves the complete transition stack on the
            // final surface's XY footprint. At a sloped outside edge that
            // footprint may lie beyond the lower planar slice, producing a
            // row of floating endpoints which variable-width gap fill cannot
            // make structural. Anchor the bottom shell to the real supporting
            // slice and continuously release that lateral correction as the
            // stack approaches the model surface. Interior points already
            // supported by the planar slice remain stationary in XY.
            if (blend < 1. - EPSILON) {
                const auto key = std::make_pair(support_idx, support_inset);
                auto cached = support_cache.find(key);
                assert(cached != support_cache.end());
                bool inset_collapsed = cached->second.empty();
                if (inset_collapsed && perimeter) {
                    // Arachne can fill a narrow annulus with a widened medial
                    // wall: inset_idx * spacing is no longer its distance from
                    // the boundary. Preserve that centerline when its complete
                    // bead fits; never project it onto the raw model edge.
                    cached = support_cache.find({support_idx, quantized_support_inset(
                        scale_(std::max(0., 0.5 * double(path.width) - 0.02)))});
                    assert(cached != support_cache.end());
                }
                const Point surface_xy = point.to_point();
                // Perimeters always morph from the actual lower contour,
                // even when the final point happens to fall inside the solid
                // support polygon. Fill points only need correction when they
                // extend beyond printable support. Limit the correction to a
                // locally reachable boundary so interior raster lines do not
                // collapse onto a remote contour.
                std::optional<Point> native_support_xy;
                if (perimeter) {
                    const auto walls = wall_support_cache.find({support_idx, std::max(0, path.inset_idx)});
                    if (walls == wall_support_cache.end() || walls->second.empty())
                        return false;
                    double closest = std::numeric_limits<double>::infinity();
                    for (const Polyline &line : walls->second) {
                        const Point candidate = surface_xy.projection_onto(line);
                        const double distance = (candidate - surface_xy).cast<double>().squaredNorm();
                        if (distance < closest) {
                            closest = distance;
                            native_support_xy = candidate;
                        }
                    }
                }
                if (native_support_xy || !cached->second.empty()) {
                    const Point support_xy = native_support_xy ? *native_support_xy :
                        projection_onto(cached->second, surface_xy);
                    const bool outside_support =
                        !Geometry::contains(cached->second, surface_xy);
                    // A shallow roof moves its contour laterally by only the
                    // run implied by this finishing stack.  At an organic
                    // crown, however, the nearest boundary of the much wider
                    // lower slice may be several millimetres away. Pulling a
                    // small crown ring all the way to that remote boundary
                    // creates a false overhang, crosses the surrounding
                    // staircase, and can make an otherwise printable top cap
                    // fail clearance. Morph supported walls only when their
                    // matching lower contour is locally reachable during this
                    // transition; unsupported endpoints are always anchored.
                    const double boundary_distance =
                        unscale_((support_xy - surface_xy).cast<double>().norm());
                    const double maximum_local_morph =
                        total_shell_depth + 0.5 * double(path.width);
                    // An unsupported bottom-course point may never remain on
                    // the exposed footprint merely because its nearest real
                    // support contour is farther away than the preferred
                    // morph distance. That was the source of the long random
                    // floating scanline ends at the Benchy roof/eave. Anchor
                    // every genuinely outside point; the aligned-path
                    // crossing and clearance passes below reject a fold which
                    // cannot be connected safely. Apply the distance cap only
                    // to already-supported perimeter points, where moving to
                    // a remote organic contour would be unnecessary.
                    if (outside_support ||
                        (perimeter && (native_support_xy || !inset_collapsed) &&
                         boundary_distance <= maximum_local_morph + EPSILON)) {
                        point.x() = coord_t(std::llround(
                            (1. - blend) * double(support_xy.x()) +
                            blend * double(surface_xy.x())));
                        point.y() = coord_t(std::llround(
                            (1. - blend) * double(support_xy.y()) +
                            blend * double(surface_xy.y())));
                    }
                }
            }
            // XY morphing changes the location at which the common height
            // field must be evaluated. Reusing the height from the final
            // contour at another XY tilts the roof differently from attached
            // features, which already query the field at their actual XY.
            if (blend > EPSILON && point.to_point() != original_xy) {
                const auto projected_z = patch_surface_z_at_xy(patch, mesh, point.to_point());
                if (!projected_z)
                    return false;
                surface_z = *projected_z - mesh.ground_level();
            }
            // Physical courses are parallel offsets from the selected mesh
            // surface.  `surface_depth` used to be ignored here; all courses
            // were instead interpolated from one global anchor plane.  That
            // made the high side of a shallow roof climb through dozens of
            // synthetic Z levels and left isolated wall rings in their wake.
            const double transition_plane_z =
                surface_reference_z - surface_depth;
            const double transitioned_z = transition_plane_z +
                blend * (surface_z - surface_reference_z);
            point.z() = scale_(transitioned_z - destination_owner_z);
        }
        path.polyline.remove_duplicate_points();
        path.polyline.fitting_result.clear();
        return true;
    };
    bool succeeded = true;
    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
        succeeded = blend_path(*path);
    } else if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity)) {
        for (ExtrusionPath &path : multipath->paths)
            succeeded &= blend_path(path);
    } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
        for (ExtrusionPath &path : loop->paths)
            succeeded &= blend_path(path);
    } else if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
        for (ExtrusionEntity *child : collection->entities)
            succeeded &= blend_shell_from_planar_support(
                *child, layers, patch, mesh, source_owner_z, destination_owner_z,
                transition_anchor_layer, total_shell_depth,
                fill_support_inset, surface_reference_z, surface_depth, blend,
                use_local_support, support_cache, wall_support_cache);
    }
    return succeeded;
}

ExtrusionEntitiesPtr clone_nonplanar_shell(
    const ExtrusionEntitiesPtr &source, const std::vector<Layer *> &layers,
    const NonplanarPatch &patch, const sla::IndexedMesh &mesh,
    coordf_t source_owner_z, coordf_t destination_owner_z,
    size_t transition_anchor_layer, double total_shell_depth,
    coord_t fill_support_inset, double surface_reference_z,
    double surface_depth, double blend, bool use_local_support,
    bool transition,
    const std::map<std::pair<size_t, coord_t>, ExPolygons> &support_cache,
    const std::map<std::pair<size_t, int>, Polylines> &wall_support_cache)
{
    ExtrusionEntitiesPtr result;
    result.reserve(source.size());
    for (const ExtrusionEntity *entity : source) {
        ExtrusionEntity *copy = entity->clone();
        if (!blend_shell_from_planar_support(
            *copy, layers, patch, mesh, source_owner_z, destination_owner_z,
        transition_anchor_layer, total_shell_depth,
        fill_support_inset, surface_reference_z, surface_depth, blend,
        use_local_support, support_cache, wall_support_cache)) {
            delete copy;
            destroy_entities(result);
            return {};
        }
        set_nonplanar_transition(*copy, transition);
        result.emplace_back(copy);
    }
    return result;
}

void collect_shell_support_insets(const ExtrusionEntity &entity,
                                  coord_t fill_support_inset,
                                  std::set<coord_t> &insets)
{
    const auto collect_path = [&](const ExtrusionPath &path) {
        if (!is_perimeter(path.role())) {
            insets.insert(quantized_support_inset(fill_support_inset));
            return;
        }
        const Flow path_flow(path.width, path.height, 0.f);
        insets.insert(quantized_support_inset(
            scale_(std::max(0., 0.5 * double(path.width) - 0.02))));
        insets.insert(quantized_support_inset(
            scale_(std::max(0., 0.5 * double(path.width) - 0.02)) +
            coord_t(std::max(0, path.inset_idx)) * path_flow.scaled_spacing()));
    };
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity))
        collect_path(*path);
    else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity))
        for (const ExtrusionPath &path : multipath->paths)
            collect_path(path);
    else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity))
        for (const ExtrusionPath &path : loop->paths)
            collect_path(path);
    else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity))
        for (const ExtrusionEntity *child : collection->entities)
            collect_shell_support_insets(*child, fill_support_inset, insets);
}

bool region_has_nonplanar_surface(const LayerRegion &region)
{
    return entity_has_nonplanar_surface(region.perimeters) ||
           entity_has_nonplanar_surface(region.fills);
}

void include_nonplanar_z_range(const ExtrusionEntity &entity, coord_t layer_z,
                               coord_t &minimum_z, coord_t &maximum_z)
{
    auto include_path = [layer_z, &minimum_z, &maximum_z](const ExtrusionPath &path) {
        if (!path.nonplanar_surface || path.nonplanar_feature_transition)
            return;
        for (const Point3 &point : path.polyline.points) {
            const coord_t absolute_z = layer_z + point.z();
            minimum_z = std::min(minimum_z, absolute_z);
            maximum_z = std::max(maximum_z, absolute_z);
        }
    };
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity)) {
        include_path(*path);
    } else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
        for (const ExtrusionPath &path : multipath->paths)
            include_path(path);
    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
        for (const ExtrusionPath &path : loop->paths)
            include_path(path);
    } else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        for (const ExtrusionEntity *child : collection->entities)
            include_nonplanar_z_range(*child, layer_z, minimum_z, maximum_z);
    }
}

void restore_nonplanar_as_planar(ExtrusionEntity &entity)
{
    auto restore_path = [](ExtrusionPath &path) {
        if (!path.nonplanar_surface || path.nonplanar_feature_transition ||
            path.nonplanar_clearance_validated)
            return;
        for (Point3 &point : path.polyline.points)
            point.z() = 0;
        path.z_contoured = false;
        path.nonplanar_surface = false;
        path.nonplanar_schedule_owned = false;
        path.nonplanar_transition = false;
        path.nonplanar_clearance_validated = false;
    };
    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
        restore_path(*path);
    } else if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity)) {
        for (ExtrusionPath &path : multipath->paths)
            restore_path(path);
    } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
        for (ExtrusionPath &path : loop->paths)
            restore_path(path);
    } else if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
        for (ExtrusionEntity *child : collection->entities)
            restore_nonplanar_as_planar(*child);
    }
}

void append_ordered_entities(ExtrusionEntityCollection &destination,
                             ExtrusionEntitiesPtr &&entities,
                             coord_t z_delta)
{
    if (entities.empty())
        return;

    auto *ordered = new ExtrusionEntityCollection;
    ordered->no_sort = true;
    ordered->entities = std::move(entities);
    for (ExtrusionEntity *entity : ordered->entities)
        offset_entity_z(*entity, z_delta);
    destination.entities.emplace_back(ordered);
    destination.no_sort = true;
}

struct NonplanarPatch {
    ExPolygons           projection;
    std::vector<uint8_t> facets;
    // Dense membership is retained for O(1) ray-hit filtering.  Iterating it
    // is prohibitively expensive on detailed meshes, though: a Benchy patch
    // may own a few hundred faces in a 225k-face model, and shell-band commit
    // visits that patch for many layers.  Keep the compact member list as the
    // iteration source so cost follows patch complexity rather than complete
    // model tessellation.
    std::vector<size_t>  face_indices;
    double               minimum_z {0.};
    double               maximum_z {0.};
    bool                 top_cap {false};
};

std::optional<double> patch_surface_z_at_xy(
    const NonplanarPatch &patch, const sla::IndexedMesh &mesh, const Point &point)
{
    return project_nonplanar_patch_height(mesh, patch.projection, patch.facets, point);
}

} // namespace

std::optional<double> project_nonplanar_patch_height(
    const sla::IndexedMesh &mesh, const ExPolygons &projection,
    const std::vector<uint8_t> &facets, const Point &point)
{
    if (projection.empty())
        return std::nullopt;
    // Clamp only samples outside the selected surface. projection_onto()
    // always returns a boundary point, including for a point already inside
    // the polygon. Applying it unconditionally made every interior raster
    // sample inherit the Z of its nearest roof edge or feature opening while
    // retaining the original emitted XY. The result was a set of false
    // plateaus which could sit several layer heights away from the mesh.
    const bool on_surface = Geometry::contains(projection, point);
    const Point sample = on_surface ? point : projection_onto(projection, point);
    const bool in_opening = !on_surface && std::any_of(
        projection.begin(), projection.end(), [&](const ExPolygon &polygon) {
            return polygon.contour.contains(point);
        });
    const Point target = in_opening ? point : sample;
    auto selected_patch_hit = [&facets, &mesh, &target](const Point &candidate)
        -> std::optional<double> {
        const std::vector<sla::IndexedMesh::hit_result> hits = mesh.query_ray_hits(
            Vec3d(unscale_(candidate.x()), unscale_(candidate.y()),
                  mesh.ground_level() - mesh.ground_level_offset() - 0.02),
            Vec3d::UnitZ());
        // At a continuing boss/chimney, the nearest upward hit may belong to
        // the feature rather than the selected roof. Search the complete hit
        // stack and accept the first facet owned by this connected patch.
        // Choosing a later selected hit is incorrect for multi-valued meshes:
        // a distinct upper sheet may share this XY projection while the
        // printable surface for this patch is the lower selected sheet.
        for (const sla::IndexedMesh::hit_result &hit : hits)
            if (hit.is_hit() && size_t(hit.face()) < facets.size() &&
                facets[size_t(hit.face())]) {
                if (candidate == target)
                    return hit.position().z();
                // An attached feature removes roof triangles from its own
                // footprint. Copying the nearest edge Z makes that feature's
                // course a plateau rather than a continuation of the roof.
                // Extend the selected facet plane to the actual sample XY;
                // do not extend beyond the patch's exterior boundary.
                const Vec3d normal = mesh.normal_by_face_id(hit.face());
                if (normal.z() <= EPSILON)
                    return std::nullopt;
                const Vec3d position = hit.position();
                return position.z() -
                    (normal.x() * (unscale_(target.x()) - position.x()) +
                     normal.y() * (unscale_(target.y()) - position.y())) / normal.z();
            }
        return std::nullopt;
    };
    if (const std::optional<double> z = selected_patch_hit(sample))
        return z;

    // A fused continuing feature leaves no roof triangles beneath its own
    // footprint. projection_onto() consequently lands on the shared patch
    // boundary, where an exact ray can miss both adjacent triangles because
    // of floating-point edge rules. Probe a few microns around that boundary;
    // the facet filter above ensures that only the actual selected surface,
    // never the feature wall, supplies the height field.
    static constexpr std::array<std::pair<double, double>, 8> directions {{
        {1., 0.}, {-1., 0.}, {0., 1.}, {0., -1.},
        {0.7071067811865476, 0.7071067811865476},
        {-0.7071067811865476, 0.7071067811865476},
        {0.7071067811865476, -0.7071067811865476},
        {-0.7071067811865476, -0.7071067811865476}
    }};
    for (const double radius : {0.005, 0.02, 0.05})
        for (const auto &[dx, dy] : directions) {
            const Point candidate(
                sample.x() + scale_(radius * dx),
                sample.y() + scale_(radius * dy));
            if (const std::optional<double> z = selected_patch_hit(candidate))
                return z;
        }
    return std::nullopt;
}

namespace {

std::optional<double> patch_surface_z_near_xy(
    const NonplanarPatch &patch, const sla::IndexedMesh &mesh,
    const Point &point, double maximum_distance)
{
    if (const std::optional<double> z =
            patch_surface_z_at_xy(patch, mesh, point))
        return z;
    if (patch.projection.empty())
        return std::nullopt;
    // The dimensional centerline of an outside wall is roughly half a bead
    // beyond the projected upper face. Query the corresponding face boundary
    // before treating it as an unrelated planar wall. This mirrors the
    // lateral probe used during candidate contouring and prevents the surface
    // stack from leaving its eave at nominal stair height merely because an
    // exact edge ray returned the side-wall triangle.
    const Point boundary = projection_onto(patch.projection, point);
    const double distance = unscale_(
        (boundary - point).cast<double>().norm());
    return distance <= maximum_distance + EPSILON ?
        patch_surface_z_at_xy(patch, mesh, boundary) : std::nullopt;
}

// Build the transition from a real solid course generated by the ordinary
// slicer. Its scanline phase, clipping, wall contact and variable-width flow
// therefore match material which is actually deposited below it. Projecting
// a newly generated final-surface raster and cloning it downward cannot offer
// that guarantee: the final footprint and phase may have no bead beneath them
// at the flat transition foundation.
} // namespace

ExtrusionEntitiesPtr project_nonplanar_source_course(
    LayerRegion &source_region, const sla::IndexedMesh &mesh,
    const ExPolygons &surface_projection, const std::vector<uint8_t> &facets,
    const ExPolygons &course_domain,
    coordf_t destination_owner_z, double maximum_drape_height,
    const NonplanarCourseProjection *course)
{
    ExtrusionEntitiesPtr result;
    if (course != nullptr &&
        (!std::isfinite(course->profile.plane_z) || !std::isfinite(course->maximum_surface_z) ||
         !std::isfinite(course->profile.blend) || course->profile.blend < 0. || course->profile.blend > 1.))
        return result;
    const double resolution = std::max(
        0.05, source_region.region().config().nonplanar_top_surface_resolution.value);
    const double source_z = source_region.layer()->print_z;
    const auto collect = [&](auto &&self, const ExtrusionEntity &entity) -> void {
        if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity)) {
            if (course != nullptr ?
                (!is_infill(path->role()) || is_arc_fill(path->role()) || path->role() == erIroning) :
                (path->role() != erTopSolidInfill && path->role() != erSolidInfill))
                return;
            const Polyline source = path->polyline.to_polyline();
            Polylines fragments = intersection_pl(Polylines{source}, course_domain);
            restore_source_path_order(source, fragments);
            for (const Polyline &fragment : fragments) {
                if (fragment.points.size() < 2 ||
                    fragment.length() <= SCALED_EPSILON)
                    continue;
                ExtrusionPath projected(*path);
                projected.polyline.points.clear();
                projected.polyline.points.reserve(fragment.points.size());
                for (const Point &point : fragment.points)
                    projected.polyline.points.emplace_back(
                        point.x(), point.y(), coord_t(0));
                if (course == nullptr)
                    projected.set_extrusion_role(erTopSolidInfill);
                projected.z_contoured = false;
                projected.nonplanar_surface = false;
                projected.nonplanar_schedule_owned = false;
                projected.nonplanar_transition = course != nullptr && course->profile.blend < 1.;
                projected.nonplanar_clearance_validated = false;

                // This is a complete solid course, not a perimeter entering
                // and leaving a selected patch.  make_nonplanar_path() also
                // applies its path-direction entry/exit easing to ineligible
                // samples.  On a raster clipped around a chimney that turns
                // the clipped endpoint back down to source_z; cloning the
                // result then repeats the same flat endpoint through several
                // transition courses and finally creates a near-vertical
                // extrusion.  Project every sample of this already-clipped
                // course directly onto the selected surface instead.  The
                // transition is between supported courses, while each course
                // remains one coherent height field.
                Points3 draped;
                bool complete = true;
                const auto emit_draped = [&]() {
                    if (draped.size() < 2) {
                        draped.clear();
                        return;
                    }
                    ExtrusionPath section(projected);
                    section.polyline.points = std::move(draped);
                    section.polyline.simplify(scale_(
                        std::min(0.005, resolution * 0.025)));
                    if (section.polyline.points.size() >= 2) {
                        section.z_contoured = true;
                        section.nonplanar_surface = true;
                        section.nonplanar_schedule_owned = false;
                        section.polyline.fitting_result.clear();
                        result.emplace_back(section.clone_move());
                    }
                    draped.clear();
                };
                const Points &points = fragment.points;
                for (size_t segment_idx = 0;
                     complete && segment_idx + 1 < points.size();
                     ++segment_idx) {
                    const Vec2d a = unscale(points[segment_idx]);
                    const Vec2d b = unscale(points[segment_idx + 1]);
                    const Vec2d delta = b - a;
                    const int sample_count = std::max(
                        1, int(std::ceil(delta.norm() / resolution)));
                    for (int sample_idx = segment_idx == 0 ? 0 : 1;
                         sample_idx <= sample_count; ++sample_idx) {
                        const Vec2d xy = a + delta *
                            (double(sample_idx) / double(sample_count));
                        const Point sample(scale_(xy.x()), scale_(xy.y()));
                        const std::optional<double> surface_z =
                            project_nonplanar_patch_height(mesh, surface_projection, facets, sample);
                        // Mesh coordinates are centered independently of the
                        // printer's layer Z. Compare and emit in print space,
                        // just as the wall and shell projection paths do.
                        const double target_surface_z = surface_z ?
                            *surface_z - mesh.ground_level() : 0.;
                        const double print_z = course == nullptr ? target_surface_z :
                            course->profile.plane_z + course->profile.blend *
                                (target_surface_z - course->maximum_surface_z);
                        if (!surface_z ||
                            print_z < source_z - double(path->height) - EPSILON ||
                            print_z > source_z + maximum_drape_height + EPSILON) {
                            // Only domain clipping may create endpoints.
                            // Keeping either side of a failed projection
                            // creates unanchored scanline ends in the middle
                            // of the patch. Leave this fragment conventional;
                            // other complete scanlines remain candidates.
                            complete = false;
                            break;
                        }
                        draped.emplace_back(
                            sample.x(), sample.y(),
                            coord_t(scale_(print_z - destination_owner_z)));
                    }
                }
                if (complete)
                    emit_draped();
            }
        } else if (const auto *multipath =
                       dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
            for (const ExtrusionPath &path : multipath->paths)
                self(self, path);
        } else if (const auto *loop =
                       dynamic_cast<const ExtrusionLoop *>(&entity)) {
            for (const ExtrusionPath &path : loop->paths)
                self(self, path);
        } else if (const auto *collection =
                       dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
            for (const ExtrusionEntity *child : collection->entities)
                self(self, *child);
        }
    };
    collect(collect, source_region.fills);
    return result;
}

// Remove sub-bead fill remnants, not dimensional wall fragments. At this
// stage Arachne's variable-width walls are still a flat list of paths; their
// owning loops/multipaths are assembled later. A short wall here is therefore
// not evidence of an isolated extrusion and deleting it opens the wall.
//
// Remnants are removed and the later deposited-coverage calculation leaves
// the corresponding conventional material in place. This is preferable to
// rejecting an otherwise broad printable surface because Arachne produced a
// microscopic taper at one corner; only material covered by surviving paths
// is removed from the ordinary slice during commit.
void remove_short_nonplanar_fill_paths(ExtrusionEntitiesPtr &entities)
{
    entities.erase(
        std::remove_if(entities.begin(), entities.end(),
            [](ExtrusionEntity *entity) {
                const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
                if (path == nullptr || is_perimeter(path->role()) ||
                    path->polyline.points.size() < 2)
                    return false;
                const double minimum_run = std::max(
                    0.40, double(path->width));
                if (unscale<double>(path->length()) + EPSILON >= minimum_run)
                    return false;
                BOOST_LOG_TRIVIAL(debug)
                    << "Dropping isolated nonplanar extrusion length="
                    << unscale<double>(path->length())
                    << " width=" << path->width
                    << " role=" << int(path->role());
                delete entity;
                return true;
            }),
        entities.end());
}

namespace {

ExPolygons continuing_feature_footprint(
    const std::vector<Layer *> &layers, size_t region_idx, size_t surface_owner,
    const NonplanarPatch &patch, double contact_distance)
{
    if (surface_owner + 1 >= layers.size() ||
        region_idx >= layers[surface_owner + 1]->regions().size())
        return {};

    // A feature is eligible only when it continues above the finishing skin
    // and physically meets that skin. Remote towers inside the toolhead
    // clearance envelope remain ordinary collision obstacles; only connected
    // geometry such as a chimney, boss, or rib receives sloped transition
    // courses.
    const ExPolygons contact_band = offset_ex(
        patch.projection, scale_(std::max(0.05, contact_distance)));
    const ExPolygons above = to_expolygons(
        layers[surface_owner + 1]->regions()[region_idx]->slices.surfaces);
    ExPolygons result;
    for (const ExPolygon &component : above) {
        if (!intersection_ex(ExPolygons{component}, contact_band).empty())
            result.push_back(component);
    }
    return union_ex(result);
}

void close_near_complete_wall_ring(ExtrusionPath &path)
{
    if (!is_perimeter(path.role()) || path.polyline.points.size() < 3)
        return;
    const coord_t seam_length = path.polyline.points.back().distance_to(
        path.polyline.points.front());
    if (seam_length > SCALED_EPSILON &&
        seam_length <= scale_(double(path.width) + 0.05) &&
        path.polyline.to_polyline().length() >= scale_(5. * double(path.width)))
        path.polyline.points.push_back(path.polyline.points.front());
}

void close_near_complete_wall_ring(ExtrusionPaths &paths)
{
    if (paths.empty() ||
        !std::all_of(paths.begin(), paths.end(), [](const ExtrusionPath &path) {
            return is_perimeter(path.role()) &&
                   !path.polyline.points.empty();
        }))
        return;

    const Point3 &first = paths.front().polyline.points.front();
    ExtrusionPath &last_path = paths.back();
    const coord_t seam_length =
        last_path.polyline.points.back().distance_to(first);
    const double maximum_width = std::accumulate(
        paths.begin(), paths.end(), 0.,
        [](double width, const ExtrusionPath &path) {
            return std::max(width, double(path.width));
        });
    double component_length = 0.;
    for (const ExtrusionPath &path : paths)
        component_length += path.polyline.to_polyline().length();

    // A gap no wider than one bead is a seam in an otherwise complete ring,
    // not a model opening. Larger openings retain their sliced topology.
    if (seam_length > SCALED_EPSILON &&
        seam_length <= scale_(maximum_width + 0.05) &&
        component_length >= scale_(5. * maximum_width))
        last_path.polyline.points.push_back(first);
}

bool warp_feature_path_to_surface(
    ExtrusionPath &path, const ExPolygons &feature_footprint,
    const NonplanarPatch &patch, const sla::IndexedMesh &mesh,
    double owner_z, double base_z, double surface_reference_z, double blend,
    bool before_current_layer, size_t feature_course,
    const std::vector<Layer *> *transition_layers = nullptr,
    size_t transition_anchor_layer = size_t(-1))
{
    if (path.nonplanar_surface || path.polyline.points.empty())
        return false;

    Points3 warped = path.polyline.points;
    for (Point3 &point : warped) {
        const std::optional<double> surface_z =
            patch_surface_z_at_xy(patch, mesh, point.to_point());
        if (!surface_z)
            return false;
        const double absolute_surface_z = *surface_z - mesh.ground_level();
        // The continuing feature must use the exact same per-XY height field
        // as the finishing shell it depends on. Using one support Z derived
        // from the patch high point made every chimney course converge on the
        // same roof surface: the low side of adjacent walls ended up only
        // 0.03 mm apart while retaining nominal layer flow. Derive the local
        // planar anchor from the sampled surface, just as
        // blend_shell_from_planar_support() does for the roof. Successive
        // feature courses then remain a printable stack while progressively
        // adopting the roof slope.
        double target_z = base_z + blend *
            (absolute_surface_z - surface_reference_z);
        // Transport the structural wall course as well as its nominal plane.
        // Discarding the brick offset while keeping its role and flow lowers
        // the first feature wall into the preceding raised course (notably
        // the chimney at a 50% brick offset). These continuing feature walls
        // are not the flat finishing cap: preserve their original interlock
        // until the ordinary source course itself returns to nominal Z.
        if (path.staggered_perimeter || path.staggered_transition)
            target_z += unscale_(point.z());
        // Real continuing-feature courses pass transition_layers so their
        // caller supplies a different base_z for every extracted source
        // course. Keep that physical course offset and add the same blended
        // roof displacement used by the surface stack. Resetting all courses
        // to the common anchor here collapses chimney/boss rings at the low
        // edge; projecting every course fully onto the surface moves the
        // stack ahead of its deposited roof support. The common expression
        // above preserves both separation and dependency order.
        (void) transition_layers;
        (void) transition_anchor_layer;
        point.z() = scale_(target_z - owner_z);
    }
    path.polyline.points = std::move(warped);
    close_near_complete_wall_ring(path);
    path.polyline.fitting_result.clear();
    const auto [minimum_z, maximum_z] = std::minmax_element(
        path.polyline.points.begin(), path.polyline.points.end(),
        [](const Point3 &left, const Point3 &right) {
            return left.z() < right.z();
        });
    const bool variable_z = minimum_z != path.polyline.points.end() &&
        maximum_z != path.polyline.points.end() &&
        maximum_z->z() - minimum_z->z() >= scale_(0.001);
    // A locally horizontal part of a dependency course is still ordered with
    // that course, but it is an ordinary wall rather than a non-planar
    // transition. Keeping the transition flag on these fixed-Z bore rings
    // produced the isolated loops visible in preview and falsely claimed a Z
    // blend where none existed. Scheduler ownership preserves dependency
    // order without changing its structural wall role.
    path.z_contoured = variable_z;
    path.nonplanar_surface = variable_z;
    path.nonplanar_schedule_owned = !variable_z;
    path.nonplanar_transition = variable_z;
    path.nonplanar_feature_transition = true;
    path.nonplanar_feature_course = feature_course;
    path.nonplanar_clearance_validated = false;
    path.nonplanar_before_current_layer = before_current_layer;
    return true;
}

size_t warp_feature_entity_to_surface(
    ExtrusionEntity &entity, const ExPolygons &feature_footprint,
    const NonplanarPatch &patch, const sla::IndexedMesh &mesh,
    double owner_z, double base_z, double surface_reference_z, double blend,
    bool before_current_layer, size_t feature_course,
    const std::vector<Layer *> *transition_layers = nullptr,
    size_t transition_anchor_layer = size_t(-1))
{
    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity))
        return warp_feature_path_to_surface(
            *path, feature_footprint, patch, mesh, owner_z, base_z,
            surface_reference_z, blend, before_current_layer,
            feature_course, transition_layers, transition_anchor_layer) ? 1 : 0;
    size_t count = 0;
    if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity)) {
        for (ExtrusionPath &path : multipath->paths) {
            if (multipath->inset_idx >= 0 && is_perimeter(path.role()))
                path.inset_idx = multipath->inset_idx;
            count += warp_feature_path_to_surface(
                path, feature_footprint, patch, mesh, owner_z, base_z,
                surface_reference_z, blend, before_current_layer,
                feature_course, transition_layers, transition_anchor_layer);
        }
        close_near_complete_wall_ring(multipath->paths);
    } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
        for (ExtrusionPath &path : loop->paths) {
            if (loop->inset_idx >= 0 && is_perimeter(path.role()))
                path.inset_idx = loop->inset_idx;
            count += warp_feature_path_to_surface(
                path, feature_footprint, patch, mesh, owner_z, base_z,
                surface_reference_z, blend, before_current_layer,
                feature_course, transition_layers, transition_anchor_layer);
        }
        close_near_complete_wall_ring(loop->paths);
    } else if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
        for (ExtrusionEntity *child : collection->entities)
            count += warp_feature_entity_to_surface(
                *child, feature_footprint, patch, mesh, owner_z, base_z,
                surface_reference_z, blend, before_current_layer,
                feature_course, transition_layers, transition_anchor_layer);
    }
    return count;
}

size_t count_warpable_feature_paths(const ExtrusionEntity &entity)
{
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity))
        return !path->nonplanar_surface && !path->polyline.points.empty();
    size_t count = 0;
    if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
        for (const ExtrusionPath &path : multipath->paths)
            count += !path.nonplanar_surface && !path.polyline.points.empty();
    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
        for (const ExtrusionPath &path : loop->paths)
            count += !path.nonplanar_surface && !path.polyline.points.empty();
    } else if (const auto *collection =
                   dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        for (const ExtrusionEntity *child : collection->entities)
            count += count_warpable_feature_paths(*child);
    }
    return count;
}

void normalize_fixed_feature_transition(ExtrusionEntity &entity)
{
    const auto normalize_path = [](ExtrusionPath &path) {
        if (!path.nonplanar_feature_transition ||
            path.nonplanar_leveling_transition ||
            path.polyline.points.empty())
            return;
        const auto [minimum_z, maximum_z] = std::minmax_element(
            path.polyline.points.begin(), path.polyline.points.end(),
            [](const Point3 &left, const Point3 &right) {
                return left.z() < right.z();
            });
        if (maximum_z->z() - minimum_z->z() >= scale_(0.001))
            return;
        // A ring may be horizontal while intentionally owned by a later
        // nominal layer for dependency scheduling. Its constant point Z is
        // then a real physical offset, not disposable transition metadata.
        // Clearing z_contoured made the writer ignore that offset and emit
        // the ring at m_nominal_z, where it collided with a later course.
        const bool offset_from_owner =
            std::abs(minimum_z->z()) >= scale_(0.001);
        path.z_contoured = offset_from_owner;
        path.nonplanar_surface = offset_from_owner;
        path.nonplanar_schedule_owned = true;
        path.nonplanar_transition = false;
    };
    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
        normalize_path(*path);
    } else if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity)) {
        for (ExtrusionPath &path : multipath->paths)
            normalize_path(path);
    } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
        for (ExtrusionPath &path : loop->paths)
            normalize_path(path);
    } else if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
        for (ExtrusionEntity *child : collection->entities)
            normalize_fixed_feature_transition(*child);
    }
}

bool feature_path_touches(const ExtrusionPath &path,
                          const ExPolygons &membership_area)
{
    // A path merely crossing a continuing feature still belongs to the
    // supporting surface. Defer only paths fully owned by the feature; this
    // prevents a chimney footprint from capturing long roof/deck paths.
    if (path.nonplanar_surface || path.polyline.points.empty())
        return false;
    const Polyline planar = path.polyline.to_polyline();
    const double total_length = planar.length();
    if (total_length <= SCALED_EPSILON)
        return false;

    // Vertex counts are not a measure of geometric ownership. A planar infill
    // chain may have many closely spaced vertices while it curves around a
    // chimney and only two vertices on a long diagonal across the roof. The
    // old majority-of-vertices test consequently captured and warped the
    // complete roof chain, deleting its foundation and producing crossings
    // with retained infill. Classify by clipped physical length instead.
    double owned_length = 0.;
    for (const Polyline &fragment : intersection_pl(Polylines{planar}, membership_area))
        owned_length += fragment.length();

    // The membership area already includes half a bead of cross-section
    // tolerance. Requiring a strong length majority retains slightly changing
    // feature rings while rejecting supporting paths that merely traverse the
    // feature on their way across the surrounding surface.
    constexpr double minimum_owned_fraction = 0.75;
    return owned_length >= minimum_owned_fraction * total_length;
}

bool connected_feature_paths_touch(const std::vector<ExtrusionPath> &paths,
                                   const ExPolygons &membership_area)
{
    // Arachne represents one continuous wall as a loop or multipath whose
    // children may change width and role. Classify that connected entity by
    // aggregate physical length. Applying the threshold independently to
    // each child cuts a ring at every short width transition and the G-code
    // scheduler subsequently inserts retracts through the middle wall.
    double total_length = 0.;
    double owned_length = 0.;
    for (const ExtrusionPath &path : paths) {
        if (path.nonplanar_surface || path.polyline.points.empty())
            continue;
        const Polyline planar = path.polyline.to_polyline();
        const double path_length = planar.length();
        if (path_length <= SCALED_EPSILON)
            continue;
        total_length += path_length;
        for (const Polyline &fragment :
             intersection_pl(Polylines{planar}, membership_area))
            owned_length += fragment.length();
    }
    constexpr double minimum_owned_fraction = 0.75;
    return total_length > SCALED_EPSILON &&
           owned_length >= minimum_owned_fraction * total_length;
}

bool entity_has_feature_path(const ExtrusionEntity &entity,
                             const ExPolygons &membership_area)
{
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity))
        return feature_path_touches(*path, membership_area) &&
               is_perimeter(path->role());
    if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity))
        return connected_feature_paths_touch(multipath->paths, membership_area) &&
               std::any_of(multipath->paths.begin(), multipath->paths.end(),
                   [](const ExtrusionPath &path) {
                       return is_perimeter(path.role());
                   });
    if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity))
        return connected_feature_paths_touch(loop->paths, membership_area) &&
               std::any_of(loop->paths.begin(), loop->paths.end(),
                   [](const ExtrusionPath &path) {
                       return is_perimeter(path.role());
                   });
    if (const auto *collection =
            dynamic_cast<const ExtrusionEntityCollection *>(&entity))
        return std::any_of(collection->entities.begin(), collection->entities.end(),
            [&membership_area](const ExtrusionEntity *child) {
                return entity_has_feature_path(*child, membership_area);
            });
    return false;
}

void split_feature_entity(const ExtrusionEntity &entity, const ExPolygons &footprint,
                          ExtrusionEntitiesPtr &remaining, ExtrusionEntitiesPtr &deferred)
{
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity)) {
        (feature_path_touches(*path, footprint) ? deferred : remaining)
            .emplace_back(clone_path_with_metadata(*path));
    } else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
        // A multipath is connected by contract. Keep it atomic even when its
        // Arachne children have different widths or roles.
        (connected_feature_paths_touch(multipath->paths, footprint) ?
             deferred : remaining).emplace_back(entity.clone());
    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
        // Do not dismantle a closed wall into independently sortable paths.
        (connected_feature_paths_touch(loop->paths, footprint) ?
             deferred : remaining).emplace_back(entity.clone());
    } else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        auto *kept = new ExtrusionEntityCollection;
        kept->no_sort = collection->no_sort;
        kept->inset_idx = collection->inset_idx;
        for (const ExtrusionEntity *child : collection->entities)
            split_feature_entity(*child, footprint, kept->entities, deferred);
        if (kept->empty())
            delete kept;
        else
            remaining.emplace_back(kept);
    }
}

void extract_feature_entities(ExtrusionEntityCollection &collection,
                              const ExPolygons &footprint,
                              ExtrusionEntitiesPtr &deferred)
{
    // Arachne may emit one wall ring as several top-level paths rather than
    // an ExtrusionLoop. Once any such fragment is owned by the attached
    // feature, follow exact endpoint connectivity and extract the complete
    // inset. Classifying each fragment independently leaves alternating arcs
    // on the conventional course and creates the missing middle chimney wall.
    std::vector<size_t> standalone;
    for (size_t idx = 0; idx < collection.entities.size(); ++idx)
        if (const auto *path =
                dynamic_cast<const ExtrusionPath *>(collection.entities[idx]);
            path != nullptr && is_perimeter(path->role()) &&
            path->polyline.points.size() >= 2)
            standalone.push_back(idx);
    std::vector<uint8_t> extracted(collection.entities.size(), 0);
    const coord_t join_tolerance = scale_(0.05);
    const auto endpoint_connected = [join_tolerance](
            const ExtrusionPath &left, const ExtrusionPath &right) {
        return std::min({
            left.first_point().distance_to(right.first_point()),
            left.first_point().distance_to(right.last_point()),
            left.last_point().distance_to(right.first_point()),
            left.last_point().distance_to(right.last_point())}) <= join_tolerance;
    };
    for (size_t seed : standalone) {
        if (extracted[seed])
            continue;
        const auto &seed_path =
            *dynamic_cast<const ExtrusionPath *>(collection.entities[seed]);
        if (!feature_path_touches(seed_path, footprint))
            continue;
        std::vector<size_t> component{seed};
        extracted[seed] = 1;
        for (size_t pending = 0; pending < component.size(); ++pending) {
            const auto &path = *dynamic_cast<const ExtrusionPath *>(
                collection.entities[component[pending]]);
            for (size_t candidate : standalone) {
                if (extracted[candidate])
                    continue;
                const auto &other = *dynamic_cast<const ExtrusionPath *>(
                    collection.entities[candidate]);
                if (other.inset_idx == path.inset_idx &&
                    endpoint_connected(path, other)) {
                    extracted[candidate] = 1;
                    component.push_back(candidate);
                }
            }
        }
        ExtrusionPaths connected;
        connected.reserve(component.size());
        for (size_t idx : component)
            connected.emplace_back(*dynamic_cast<const ExtrusionPath *>(
                collection.entities[idx]));
        chain_and_reorder_extrusion_paths(connected);
        close_near_complete_wall_ring(connected);
        auto *multipath = new ExtrusionMultiPath(std::move(connected));
        multipath->inset_idx = seed_path.inset_idx;
        multipath->set_reverse();
        deferred.emplace_back(multipath);
    }

    ExtrusionEntitiesPtr remaining;
    remaining.reserve(collection.entities.size());
    for (size_t idx = 0; idx < collection.entities.size(); ++idx)
        if (!extracted[idx])
            split_feature_entity(
                *collection.entities[idx], footprint, remaining, deferred);
    collection.clear();
    collection.entities = std::move(remaining);
}

std::optional<std::pair<double, double>> feature_surface_z_range(
    const ExPolygons &footprint, const NonplanarPatch &patch,
    const sla::IndexedMesh &mesh)
{
    double minimum = std::numeric_limits<double>::max();
    double maximum = std::numeric_limits<double>::lowest();
    for (const ExPolygon &polygon : footprint)
        for (size_t contour_idx = 0; contour_idx <= polygon.holes.size(); ++contour_idx) {
            const Polygon &contour = polygon.contour_or_hole(contour_idx);
            for (const Point &point : contour.points)
                if (const std::optional<double> z =
                        patch_surface_z_at_xy(patch, mesh, point)) {
                    minimum = std::min(minimum, *z - mesh.ground_level());
                    maximum = std::max(maximum, *z - mesh.ground_level());
                }
        }
    if (minimum > maximum)
        return std::nullopt;
    return std::pair<double, double>{minimum, maximum};
}

void restore_deferred_feature(ExtrusionEntity &entity)
{
    auto restore_path = [](ExtrusionPath &path) {
        for (Point3 &point : path.polyline.points)
            point.z() = 0;
        path.z_contoured = false;
        path.nonplanar_surface = false;
        path.nonplanar_schedule_owned = false;
        path.nonplanar_transition = false;
        path.nonplanar_clearance_validated = false;
        path.nonplanar_before_current_layer = false;
        path.nonplanar_feature_transition = false;
        path.nonplanar_feature_course = 0;
        path.nonplanar_leveling_transition = false;
    };
    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
        restore_path(*path);
    } else if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity)) {
        for (ExtrusionPath &path : multipath->paths)
            restore_path(path);
    } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
        for (ExtrusionPath &path : loop->paths)
            restore_path(path);
    } else if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
        for (ExtrusionEntity *child : collection->entities)
            restore_deferred_feature(*child);
    }
}

void mark_nonplanar_leveling_transition(ExtrusionEntity &entity)
{
    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
        if (path->nonplanar_feature_transition)
            path->nonplanar_leveling_transition = true;
    } else if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity)) {
        for (ExtrusionPath &path : multipath->paths)
            if (path.nonplanar_feature_transition)
                path.nonplanar_leveling_transition = true;
    } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
        for (ExtrusionPath &path : loop->paths)
            if (path.nonplanar_feature_transition)
                path.nonplanar_leveling_transition = true;
    } else if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
        for (ExtrusionEntity *child : collection->entities)
            mark_nonplanar_leveling_transition(*child);
    }
}

std::vector<NonplanarPatch> build_nonplanar_patches(
    const sla::IndexedMesh &mesh, const std::vector<uint8_t> &selected_facets)
{
    const indexed_triangle_set &its = *mesh.get_triangle_mesh();
    const std::vector<Vec3i32> neighbors = its_face_neighbors_par(its);
    std::vector<uint8_t> visited(selected_facets.size(), 0);
    std::vector<NonplanarPatch> result;
    std::queue<size_t> pending;

    for (size_t seed = 0; seed < selected_facets.size(); ++seed) {
        if (!selected_facets[seed] || visited[seed])
            continue;

        NonplanarPatch patch;
        patch.facets.assign(selected_facets.size(), 0);
        patch.minimum_z = std::numeric_limits<double>::max();
        patch.maximum_z = std::numeric_limits<double>::lowest();
        Polygons triangles;
        pending.push(seed);
        visited[seed] = 1;
        while (!pending.empty()) {
            const size_t face = pending.front();
            pending.pop();
            patch.facets[face] = 1;
            patch.face_indices.push_back(face);
            const Vec3i32 &indices = its.indices[face];
            Polygon triangle;
            triangle.points.reserve(3);
            for (int vertex = 0; vertex < 3; ++vertex) {
                const Vec3f &point = its.vertices[size_t(indices[vertex])];
                triangle.points.emplace_back(scale_(point.x()), scale_(point.y()));
                patch.minimum_z = std::min(patch.minimum_z, double(point.z()));
                patch.maximum_z = std::max(patch.maximum_z, double(point.z()));
            }
            if (triangle.area() < 0)
                triangle.reverse();
            triangles.emplace_back(std::move(triangle));
            for (int edge = 0; edge < 3; ++edge) {
                const int neighbor = neighbors[face][edge];
                if (neighbor >= 0 && selected_facets[size_t(neighbor)] &&
                    !visited[size_t(neighbor)]) {
                    visited[size_t(neighbor)] = 1;
                    pending.push(size_t(neighbor));
                }
            }
        }

        patch.projection = union_ex(triangles);
        if (!patch.projection.empty())
            result.emplace_back(std::move(patch));
    }
    return result;
}

std::vector<Vec3d> clip_vertices_at_z(const std::vector<Vec3d> &input,
                                      double limit, bool keep_above);

NonplanarPatch clip_patch_to_z_range(const NonplanarPatch &patch,
                                     const sla::IndexedMesh &mesh,
                                     double minimum_z, double maximum_z)
{
    NonplanarPatch clipped;
    clipped.facets.assign(patch.facets.size(), 0);
    clipped.minimum_z = std::numeric_limits<double>::max();
    clipped.maximum_z = std::numeric_limits<double>::lowest();
    clipped.top_cap = patch.top_cap;

    const indexed_triangle_set &its = *mesh.get_triangle_mesh();
    Polygons polygons;
    for (const size_t face : patch.face_indices) {
        const Vec3i32 &indices = its.indices[face];
        std::vector<Vec3d> vertices{
            its.vertices[size_t(indices[0])].cast<double>(),
            its.vertices[size_t(indices[1])].cast<double>(),
            its.vertices[size_t(indices[2])].cast<double>()};
        vertices = clip_vertices_at_z(vertices, minimum_z, true);
        vertices = clip_vertices_at_z(vertices, maximum_z, false);
        if (vertices.size() < 3)
            continue;

        Polygon polygon;
        polygon.points.reserve(vertices.size());
        for (const Vec3d &vertex : vertices) {
            polygon.points.emplace_back(scale_(vertex.x()), scale_(vertex.y()));
            clipped.minimum_z = std::min(clipped.minimum_z, vertex.z());
            clipped.maximum_z = std::max(clipped.maximum_z, vertex.z());
        }
        if (polygon.area() < 0)
            polygon.reverse();
        if (!polygon.is_valid())
            continue;
        clipped.facets[face] = 1;
        clipped.face_indices.push_back(face);
        polygons.emplace_back(std::move(polygon));
    }
    clipped.projection = union_ex(polygons);
    return clipped;
}

std::vector<NonplanarPatch> active_nonplanar_patches(
    const std::vector<NonplanarPatch> &candidates,
    const sla::IndexedMesh &mesh,
    const ExPolygons &seed_coverage, const std::vector<uint8_t> &emitted_facets,
    double active_minimum_z, double active_maximum_z,
    double first_layer_height, double last_layer_height)
{
    std::vector<NonplanarPatch> result;
    for (const NonplanarPatch &patch : candidates) {
        if (patch.maximum_z + EPSILON < active_minimum_z ||
            patch.minimum_z - EPSILON > active_maximum_z)
            continue;

        // Seeds are inset extrusion centerlines, not the mesh boundary. On
        // a shallow roof their Z extrema omit the first/last bead of the
        // surface. Clipping there creates an artificial wall inside the roof
        // which has no matching wall on the flat foundation. Use the complete
        // connected component when it fits within the seed courses' height
        // envelope. Longer organic slopes still use a local band. This only
        // selects a candidate; support and toolhead checks remain mandatory.
        const bool complete_component =
            patch.minimum_z >= active_minimum_z - first_layer_height - EPSILON &&
            patch.maximum_z <= active_maximum_z + last_layer_height + EPSILON;
        NonplanarPatch active = complete_component ? patch : clip_patch_to_z_range(
            patch, mesh, std::max(patch.minimum_z, active_minimum_z),
            std::min(patch.maximum_z, active_maximum_z));
        if (active.projection.empty())
            continue;

        // A clipped band is emitted atomically, so only a facet participating
        // in this band suppresses it. A previously emitted lower band must not
        // hide a disjoint top cap from a later run.
        bool emitted = false;
        for (const size_t face : active.face_indices)
            if (face < emitted_facets.size() && emitted_facets[face]) {
                emitted = true;
                break;
            }
        if (emitted)
            continue;

        // XY overlap alone is insufficient: on models such as Benchy, a low
        // hull/deck path may project beneath the cabin roof. Associate a patch
        // with the path run that actually occupies its surface-height range.
        if (!intersection_ex(active.projection, seed_coverage).empty()) {
            BOOST_LOG_TRIVIAL(debug)
                << "Activated nonplanar surface band z=[" << active.minimum_z
                << ',' << active.maximum_z << "] from component z=["
                << patch.minimum_z << ',' << patch.maximum_z << ']';
            result.emplace_back(std::move(active));
        }
    }
    return result;
}

bool patch_is_accessible(const NonplanarPatch &patch,
                         const std::vector<Layer *> &layers, size_t region_idx,
                         double ground, double clearance_angle_degrees)
{
    const double minimum_z = patch.minimum_z - ground;
    const double maximum_z = patch.maximum_z - ground;
    const double angle = Geometry::deg2rad(std::clamp(clearance_angle_degrees, 1., 89.));
    std::vector<ExPolygons> surface_bands(layers.size());
    ExPolygons remaining = patch.projection;
    size_t first = layers.size();
    size_t last = 0;

    // Partition the projection by the ordinary top step that currently owns
    // each part of the curved surface. This is the collision construction used
    // by the reference Slic3r implementation: the steps are support geometry,
    // while the offset between them models the toolhead clearance cone.
    for (size_t idx = layers.size(); idx-- > 0;) {
        const Layer &layer = *layers[idx];
        if (layer.slice_z > maximum_z + layer.height)
            continue;
        if (layer.slice_z + layer.height < minimum_z)
            break;
        const ExPolygons upper = idx + 1 < layers.size() ?
            to_expolygons(layers[idx + 1]->regions()[region_idx]->slices.surfaces) :
            ExPolygons{};
        const ExPolygons current_top = diff_ex(
            to_expolygons(layer.regions()[region_idx]->slices.surfaces), upper,
            ApplySafetyOffset::No);
        surface_bands[idx] = intersection_ex(
            remaining, current_top, ApplySafetyOffset::Yes);
        if (!surface_bands[idx].empty()) {
            remaining = diff_ex(
                remaining, surface_bands[idx], ApplySafetyOffset::Yes);
            first = std::min(first, idx);
            last = std::max(last, idx);
        }
    }
    if (first == layers.size())
        return false;

    double collision_area = 0.;
    for (size_t high = first; high <= last; ++high) {
        const ExPolygons high_geometry = to_expolygons(
            layers[high]->regions()[region_idx]->slices.surfaces);
        for (size_t low = first; low < high; ++low) {
            if (surface_bands[low].empty())
                continue;
            const double height = layers[high]->slice_z - layers[low]->slice_z;
            const coord_t radius = scale_(height / std::tan(angle));
            if (radius <= 0)
                continue;
            const ExPolygons envelope = offset_ex(surface_bands[low], radius,
                ClipperLib::jtSquare, 100000.);
            const ExPolygons collisions = intersection_ex(
                high_geometry, diff_ex(envelope, patch.projection));
            for (const ExPolygon &collision : collisions)
                collision_area += std::abs(collision.area()) *
                    SCALING_FACTOR * SCALING_FACTOR;
            // Ignore sub-bead tessellation flecks, but reject a patch before
            // toolpath generation when a printable amount of geometry is in
            // the clearance envelope. Never flatten it later in G-code.
            if (collision_area > 0.5)
                return false;
        }
    }
    return true;
}

const ExtrusionPath *first_nonplanar_path(const ExtrusionEntity &entity, ExtrusionRole role)
{
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity))
        return path->nonplanar_surface && path->role() == role ? path : nullptr;
    if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
        for (const ExtrusionPath &path : multipath->paths)
            if (path.nonplanar_surface && path.role() == role)
                return &path;
    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
        for (const ExtrusionPath &path : loop->paths)
            if (path.nonplanar_surface && path.role() == role)
                return &path;
    } else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        for (const ExtrusionEntity *child : collection->entities)
            if (const ExtrusionPath *path = first_nonplanar_path(*child, role))
                return path;
    }
    return nullptr;
}

const ExtrusionPath *first_path_with_role(
    const ExtrusionEntity &entity, ExtrusionRole role)
{
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity))
        return path->role() == role ? path : nullptr;
    if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
        for (const ExtrusionPath &path : multipath->paths)
            if (path.role() == role)
                return &path;
    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
        for (const ExtrusionPath &path : loop->paths)
            if (path.role() == role)
                return &path;
    } else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        for (const ExtrusionEntity *child : collection->entities)
            if (const ExtrusionPath *path = first_path_with_role(*child, role))
                return path;
    }
    return nullptr;
}

const ExtrusionPath *first_path_with_role(
    const std::vector<Layer *> &layers, size_t first_layer, size_t last_layer,
    size_t region_idx, ExtrusionRole role)
{
    for (size_t layer_idx = first_layer; layer_idx <= last_layer; ++layer_idx) {
        const LayerRegion &region = *layers[layer_idx]->regions()[region_idx];
        for (const ExtrusionEntity *entity : region.perimeters)
            if (const ExtrusionPath *path = first_path_with_role(*entity, role))
                return path;
        for (const ExtrusionEntity *entity : region.fills)
            if (const ExtrusionPath *path = first_path_with_role(*entity, role))
                return path;
    }
    return nullptr;
}

const ExtrusionPath *first_nonplanar_path(
    const std::vector<Layer *> &layers, size_t first_layer, size_t last_layer,
    size_t region_idx, ExtrusionRole role)
{
    for (size_t layer_idx = first_layer; layer_idx <= last_layer; ++layer_idx) {
        const LayerRegion &region = *layers[layer_idx]->regions()[region_idx];
        if (const ExtrusionPath *path = first_nonplanar_path(region.perimeters, role))
            return path;
        if (const ExtrusionPath *path = first_nonplanar_path(region.fills, role))
            return path;
    }
    return nullptr;
}

void destroy_entities(ExtrusionEntitiesPtr &entities)
{
    for (ExtrusionEntity *entity : entities)
        delete entity;
    entities.clear();
}

struct WeightedPatchGradient {
    Vec2d gradient;
    double weight {0.};
    Point centroid;
};

struct PatchRasterCache {
    std::vector<int8_t> facet_clearance;
    ExPolygons reachable_projection;
    std::vector<WeightedPatchGradient> reachable_gradients;
    bool domain_ready {false};
    bool doubly_curved {false};
};

ExtrusionEntitiesPtr make_patch_raster(
    const NonplanarPatch &patch, const sla::IndexedMesh &mesh,
    const ExtrusionPath &path_template, coordf_t owner_z, double resolution,
    const SequentialGantryGeometry &gantry, InfillPattern top_surface_pattern,
    double surface_depth = 0.,
    const ExPolygons *continuing_feature = nullptr,
    const ExtrusionPath *outline_template = nullptr,
    PatchRasterCache *shared_cache = nullptr,
    const ExtrusionPath *inner_outline_template = nullptr,
    size_t outline_count = 1,
    const Arachne::WallToolPathsParams *arachne_params = nullptr,
    bool emit_outlines = true)
{
    ExtrusionEntitiesPtr output;
    const coord_t spacing = Flow(path_template.width, path_template.height, 0.f).scaled_spacing();
    if (spacing <= 0 || patch.projection.empty())
        return output;

    const double ground = mesh.ground_level();
    PatchRasterCache local_cache;
    PatchRasterCache &cache = shared_cache == nullptr ? local_cache : *shared_cache;
    if (cache.facet_clearance.size() != patch.facets.size()) {
        cache.facet_clearance.assign(patch.facets.size(), -1);
        cache.reachable_projection.clear();
        cache.reachable_gradients.clear();
        cache.domain_ready = false;
        cache.doubly_curved = false;
    }
    std::vector<int8_t> &facet_clearance = cache.facet_clearance;
    // Include the boundary in the physical occupancy test.  Clipper's point
    // containment is intentionally strict, but a toolhead vertex touching a
    // model edge is still a collision.  A 0.01 mm expansion is below the path
    // sampling resolution and prevents a toolhead whose reach lands exactly
    // on a tessellation edge from incorrectly clearing the whole facet.
    const ExPolygons clearance_projection = offset_ex(
        patch.projection, scale_(0.01));
    const BoundingBox patch_bbox = get_extents(clearance_projection);
    // Adjacent mesh facets share vertices, so the finite-toolhead check asks
    // about the same translated toolhead samples many times. Point-in-polygon
    // testing against a detailed roof outline is substantially more expensive
    // than the tangent-plane arithmetic around it. Cache those exact scaled
    // points in a concurrent map because facet clearance is evaluated by TBB.
    // A benign duplicate insertion may race, but every thread computes the
    // same immutable containment result.
    tbb::concurrent_unordered_map<Point, uint8_t, PointHash>
        projection_contains_cache;
    auto projection_contains = [&clearance_projection, &patch_bbox,
                                &projection_contains_cache](const Point &point) {
        if (!patch_bbox.contains(point))
            return false;
        const auto cached = projection_contains_cache.find(point);
        if (cached != projection_contains_cache.end())
            return cached->second != 0;
        // The connected facet union is already the exact XY domain needed by
        // this finite-patch clearance test, including holes around continuing
        // features. A 3D mesh ray produced the same answer but traversed the
        // complete AABB tree for every toolhead vertex on every rejected
        // facet, dominating non-planar analysis on detailed meshes.
        const bool contained = Geometry::contains(clearance_projection, point);
        projection_contains_cache.emplace(point, contained ? 1 : 0);
        return contained;
    };
    auto facet_is_clear = [&patch, &mesh, &gantry, &path_template,
                           &facet_clearance, &projection_contains](size_t face) {
        if (gantry.empty())
            return true;
        int8_t &cached = facet_clearance[face];
        if (cached >= 0)
            return cached != 0;

        const Vec3d normal = mesh.normal_by_face_id(int(face));
        if (normal.z() <= EPSILON)
            return (cached = 0) != 0;
        const Vec2d gradient(-normal.x() / normal.z(), -normal.y() / normal.z());
        const double slope = gradient.norm();
        const double nozzle_reach = gantry.maximum_xy_reach(0.001);
        const double support_contact = std::max(1e-4, double(path_template.height)) +
            (nozzle_reach + 0.5 * double(path_template.width)) * slope;
        bool tangent_plane_is_clear = true;
        for (const SequentialGantrySlice &slice : gantry.slices) {
            // The obstacle supporting this shell is the preceding draped skin,
            // one bead height below the requested surface. Apply that real
            // vertical separation to every hotend cross-section; the nozzle
            // slice additionally receives its lateral bead-contact allowance.
            const double allowed_height = slice.height <= EPSILON ?
                support_contact : slice.height + double(path_template.height);
            for (const Polygon &polygon : slice.polygons) {
                for (const Point &vertex : polygon.points) {
                    const Vec2d offset = unscale(vertex);
                    // Height of the tangent surface beneath this point of the
                    // fixed-orientation toolhead. This directional support
                    // test derives the actual approach envelope from every
                    // gantry slice instead of imposing one scalar angle.
                    if (gradient.dot(offset) > allowed_height + 1e-4) {
                        tangent_plane_is_clear = false;
                        break;
                    }
                }
                if (!tangent_plane_is_clear)
                    break;
            }
            if (!tangent_plane_is_clear)
                break;
        }
        if (tangent_plane_is_clear)
            return (cached = 1) != 0;

        // The tangent-plane test above is a cheap conservative bound, but it
        // treats the local surface as if it extended infinitely beneath the
        // complete toolhead. Near a finite roof this rejects facets merely
        // because a distant heater-block corner would intersect that
        // imaginary continuation, producing triangular planar/non-planar
        // islands. Recheck rejected facets against the actual connected mesh
        // patch. Toolhead samples outside the patch are empty space, not a
        // collision.
        const indexed_triangle_set &its = *mesh.get_triangle_mesh();
        if (face >= its.indices.size())
            return (cached = 0) != 0;
        const Vec3i32 &indices = its.indices[face];
        std::array<Vec3d, 4> anchors;
        anchors[3] = Vec3d::Zero();
        for (size_t idx = 0; idx < 3; ++idx) {
            anchors[idx] = its.vertices[size_t(indices[int(idx)])].cast<double>();
            anchors[3] += anchors[idx] / 3.;
        }
        for (const Vec3d &anchor : anchors)
            for (const SequentialGantrySlice &slice : gantry.slices) {
                const double allowed_height = slice.height <= EPSILON ?
                    support_contact : slice.height + double(path_template.height);
                // If even the highest point of this finite patch is below the
                // toolhead slice, no XY sample from that slice can collide.
                // This eliminates nearly all expensive surface queries for
                // the large upper hotend cross-sections.
                if (patch.maximum_z <= anchor.z() + allowed_height + 1e-4)
                    continue;
                for (const Polygon &polygon : slice.polygons)
                    for (const Point &vertex : polygon.points) {
                        const Vec2d candidate = anchor.head<2>() + unscale(vertex);
                        const Point candidate_scaled(scale_(candidate.x()), scale_(candidate.y()));
                        // A sample outside the connected finishing patch is
                        // empty space. Within it, the facet's tangent plane is
                        // the exact surface over that triangle and remains the
                        // same conservative directional clearance test used
                        // above. Restricting that test to the finite patch is
                        // what removes false heater-block collisions at a roof
                        // edge, without launching a full-mesh ray for every
                        // facet/toolhead vertex pair.
                        if (projection_contains(candidate_scaled) &&
                            gradient.dot(unscale(vertex)) > allowed_height + 1e-4) {
                            return (cached = 0) != 0;
                        }
                    }
            }
        return (cached = 1) != 0;
    };

    // Build one coherent reachable domain from the facets allowed by the
    // detailed toolhead. Steep boundary fillets may legitimately be outside
    // the safe cone while the broad roof inside them is printable. The finite
    // patch test above prevents distant toolhead corners from cutting false
    // diagonal holes through this domain; the later all-depth check remains
    // atomic for every raster generated from it.
    const indexed_triangle_set &its = *mesh.get_triangle_mesh();
    if (!cache.domain_ready) {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, patch.face_indices.size()),
            [&](const tbb::blocked_range<size_t> &range) {
                for (size_t index = range.begin(); index < range.end(); ++index)
                    facet_is_clear(patch.face_indices[index]);
            });
        Polygons rejected_triangles;
        for (const size_t face : patch.face_indices) {
            const Vec3i32 &indices = its.indices[face];
            Polygon triangle;
            triangle.points.reserve(3);
            for (int vertex = 0; vertex < 3; ++vertex) {
                const Vec3f &point = its.vertices[size_t(indices[vertex])];
                triangle.points.emplace_back(scale_(point.x()), scale_(point.y()));
            }
            if (triangle.area() < 0)
                triangle.reverse();
            if (triangle.is_valid()) {
                if (!facet_is_clear(face)) {
                    rejected_triangles.emplace_back(std::move(triangle));
                    continue;
                }
                const Vec3d normal = mesh.normal_by_face_id(int(face));
                const double weight = std::abs(triangle.area());
                if (normal.z() > EPSILON && weight > 0.)
                    cache.reachable_gradients.push_back({
                        Vec2d(-normal.x() / normal.z(), -normal.y() / normal.z()),
                        weight,
                        Point((triangle.points[0].x() + triangle.points[1].x() +
                               triangle.points[2].x()) / 3,
                              (triangle.points[0].y() + triangle.points[1].y() +
                               triangle.points[2].y()) / 3)});
            }
        }
        // patch.projection is already the union of every selected triangle.
        // Most practical roofs reject no facets, so unioning hundreds of
        // thousands of the same mesh triangles again for toolhead filtering
        // was pure work. Subtract only the usually tiny rejected set; the
        // common all-reachable case becomes a cheap copy of the existing
        // projection.
        cache.reachable_projection = rejected_triangles.empty() ?
            patch.projection :
            diff_ex(patch.projection, union_ex(rejected_triangles));

        // A few edge facets can clear a bulky toolhead simply because most of
        // its footprint hangs beyond the model.  Emitting those isolated
        // slivers does not remove the visible staircase and, on rounded ears
        // and similar features, cuts holes out of an otherwise sound planar
        // skin.  Require the collision-safe domain to cover a meaningful part
        // of the connected candidate.  This remains shape-independent: broad
        // roofs with a rejected fillet keep nearly all of their area, while a
        // surface which the gantry can approach only at its rim falls back
        // atomically to the conventional layers.
        constexpr double minimum_coherent_coverage = 0.55;
        const double candidate_area = std::abs(area(patch.projection));
        const double reachable_area = std::abs(area(cache.reachable_projection));
        if (candidate_area > 0. &&
            reachable_area <= minimum_coherent_coverage * candidate_area) {
            BOOST_LOG_TRIVIAL(debug)
                << "Rejecting nonplanar candidate: collision-safe coverage "
                << reachable_area / candidate_area
                << " is below the coherent-surface threshold";
            cache.reachable_projection.clear();
            cache.reachable_gradients.clear();
        }
        cache.domain_ready = true;
    }
    const ExPolygons &reachable_projection = cache.reachable_projection;
    const std::vector<WeightedPatchGradient> &reachable_gradients =
        cache.reachable_gradients;
    const double patch_projection_area = std::abs(area(patch.projection));
    const bool complete_reachable_domain = patch_projection_area > 0. &&
        std::abs(std::abs(area(reachable_projection)) - patch_projection_area) <=
            std::max(1., patch_projection_area * 1e-6);
    // The replacement owns its boundary explicitly. This is important for
    // dimensional accuracy at an outside edge and for small negative details
    // such as holes and window openings. Previously the fill was extended
    // into a continuing feature to hide the missing outline; that made the
    // roof raster and the feature perimeter extrude over the same interface.
    // A draped contour on each side now owns exactly one bead, while the fill
    // is inset to meet (rather than cross) the contour.
    const ExtrusionPath *boundary_template = outline_template;
    const double boundary_width = boundary_template == nullptr ? 0. :
        std::max(0., double(boundary_template->width));
    ExPolygons solid_feature_footprint;
    if (continuing_feature != nullptr) {
        // Preserve the feature's complete slice topology. Converting every
        // ExPolygon to only its outer contour filled all of its holes before
        // subtracting the feature from the roof raster. On a hollow or
        // compound continuation (the Benchy chimney/cabin is a representative
        // case), that incorrectly assigned exposed roof inside those holes to
        // the feature and left a broad unextruded band. A genuine bore is
        // still protected: it is absent from reachable_projection, so keeping
        // the hole here cannot generate a roof path across modeled empty
        // space.
        solid_feature_footprint = union_ex(*continuing_feature);
    }
    const auto project_line = [&](const Points &points, const ExtrusionPath &source,
                                  bool closed) -> ExtrusionPath * {
        if (points.size() < 2)
            return nullptr;
        Points sampled = points;
        if (closed && sampled.front() != sampled.back())
            sampled.push_back(sampled.front());
        Polyline3 projected;
        for (size_t segment_idx = 0; segment_idx + 1 < sampled.size(); ++segment_idx) {
            const Vec2d a = unscale(sampled[segment_idx]);
            const Vec2d b = unscale(sampled[segment_idx + 1]);
            const Vec2d delta = b - a;
            const int count = std::max(1, int(std::ceil(delta.norm() / resolution)));
            for (int sample_idx = segment_idx == 0 ? 0 : 1;
                 sample_idx <= count; ++sample_idx) {
                const Vec2d xy = a + delta * (double(sample_idx) / count);
                const Point sample_point(scale_(xy.x()), scale_(xy.y()));
                const std::optional<double> surface_z =
                    patch_surface_z_at_xy(patch, mesh, sample_point);
                if (!surface_z)
                    return nullptr;
                projected.append(Point3(Vec3crd(
                    sample_point.x(), sample_point.y(),
                    scale_(*surface_z - ground - owner_z - surface_depth))));
            }
        }
        projected.simplify(scale_(std::min(0.005, resolution * 0.025)));
        const double maximum_retrace_squared =
            double(scale_(0.1)) * double(scale_(0.1));
        for (size_t point_idx = 1;
             point_idx + 1 < projected.points.size();) {
            if ((projected.points[point_idx - 1].to_point() -
                 projected.points[point_idx + 1].to_point()).squaredNorm() <=
                    double(scale_(0.01)) * double(scale_(0.01)) &&
                (projected.points[point_idx].to_point() -
                 projected.points[point_idx - 1].to_point()).squaredNorm() <=
                    maximum_retrace_squared) {
                projected.points.erase(
                    projected.points.begin() + point_idx,
                    projected.points.begin() + point_idx + 2);
                if (point_idx > 1)
                    --point_idx;
            } else {
                ++point_idx;
            }
        }
        if (projected.points.size() < 2)
            return nullptr;
        auto *path = new ExtrusionPath(std::move(projected), source);
        path->inset_idx = source.inset_idx;
        path->z_contoured = true;
        path->nonplanar_surface = true;
        path->nonplanar_schedule_owned = false;
        path->nonplanar_transition = surface_depth > EPSILON;
        path->set_reverse();
        return path;
    };

    // Generate in the projected 2D domain, then deform the complete result
    // back to the model surface (Ahlers/QuickCurve construction).
    // The boundary centerline is half a perimeter width inside the model. The
    // first fill centerline is another half perimeter plus half fill width
    // inward, so their deposited bead envelopes touch without intersecting.
    ExPolygons centerline_area;
    // Gap repair is generated before the slicer-owned wall graph is committed.
    // Reserve the bead envelope of the equivalent Arachne outlines now, so a
    // later repair path cannot mistake the perimeter band for uncovered top
    // skin and extrude through an inner/outer wall at the same Z.
    Polygons reserved_wall_coverage;
    bool generated_arachne_outlines = false;
    if (boundary_template != nullptr && outline_count > 0 &&
        arachne_params != nullptr) {
        std::map<unsigned int, size_t> projected_outline_paths_by_inset;
        const ExtrusionPath &inner_template = inner_outline_template == nullptr ?
            *boundary_template : *inner_outline_template;
        const coord_t outer_width = scale_(boundary_template->width);
        const coord_t inner_width = scale_(inner_template.width);
        // A partially reachable surface must not acquire a synthetic closed
        // perimeter around the collision boundary. Generate the real model
        // wall graph first and retain only its reachable edges; these open
        // fragments are spliced into the surviving planar graph at commit.
        ExPolygons wall_domain = complete_reachable_domain ?
            reachable_projection : patch.projection;
        if (!solid_feature_footprint.empty()) {
            // The continuing feature owns the interface wall. Fill its
            // footprint only for Arachne's boundary analysis so the roof does
            // not generate a second set of concentric walls around the same
            // chimney/boss opening. The solid footprint is removed from the
            // raster below, so this does not close a modeled bore.
            append(wall_domain, solid_feature_footprint);
            wall_domain = union_ex(wall_domain);
        }
        Polygons wall_outline = to_polygons(wall_domain);
        Arachne::WallToolPaths wall_tool_paths(
            wall_outline, outer_width, inner_width,
            outline_count, 0, path_template.height, *arachne_params);
        const std::vector<Arachne::VariableWidthLines> &wall_groups =
            wall_tool_paths.getToolPaths();
        for (const Arachne::VariableWidthLines &wall_group : wall_groups)
            for (const Arachne::ExtrusionLine &wall : wall_group) {
                if (wall.empty())
                    continue;
                const bool external = wall.inset_idx == 0;
                ExtrusionPath wall_source = external ?
                    *boundary_template : inner_template;
                wall_source.set_extrusion_role(
                    external ? erExternalPerimeter : erPerimeter);
                wall_source.inset_idx = wall.inset_idx;
                const Flow wall_flow(
                    wall_source.width, wall_source.height, 0.f);
                // Arachne deliberately tapers the ends of thin-feature walls
                // down to hairline junctions.  That is useful for fitting a
                // planar bead into a sharp corner, but a draped boundary must
                // cover and bond to the model envelope at every sampled Z.
                // Sub-nozzle tips left visible gaps at the Benchy roof eaves
                // even though their centerlines reached the mesh boundary.
                // Clamp this non-planar copy to the configured printable bead
                // minimum; do not alter Arachne's ordinary planar output.
                Arachne::ExtrusionLine printable_wall = wall;
                const coord_t configured_minimum_width = scale_(
                    std::max(0., double(arachne_params->min_bead_width)));
                const coord_t physical_minimum_width = scale_(
                    0.75 * double(wall_source.width));
                const coord_t minimum_width = std::max(
                    configured_minimum_width, physical_minimum_width);
                for (Arachne::ExtrusionJunction &junction :
                     printable_wall.junctions)
                    junction.w = std::max(junction.w, minimum_width);
                ExtrusionPaths planar_paths;
                extrusion_paths_append(
                    planar_paths, printable_wall, wall_source.role(), wall_flow);
                for (ExtrusionPath &planar_path : planar_paths) {
                    planar_path.inset_idx = wall.inset_idx;
                    const Polylines wall_fragments = complete_reachable_domain ?
                        Polylines{planar_path.polyline.to_polyline()} :
                        intersection_pl(
                            Polylines{planar_path.polyline.to_polyline()},
                            offset_ex(reachable_projection, scale_(0.01)));
                    for (const Polyline &fragment : wall_fragments)
                        if (ExtrusionPath *projected = project_line(
                                fragment.points, planar_path, false)) {
                            projected->polygons_covered_by_width(
                                reserved_wall_coverage, 0.f);
                            // Perimeters are contoured in
                            // Layer::make_nonplanar_top_surfaces(), where they
                            // retain their seam and their place in the original
                            // wall loop.  This Arachne pass is still needed to
                            // derive the variable-width inner contour for the
                            // skin, but emitting its outlines here creates a
                            // second, smaller set of closed rings which is not
                            // connected to the sliced wall graph.
                            if (emit_outlines)
                            {
                                ++projected_outline_paths_by_inset[
                                    projected->inset_idx];
                                output.emplace_back(projected);
                            }
                            else
                                delete projected;
                        }
                }
            }
        std::ostringstream outline_summary;
        for (const auto &[inset, count] : projected_outline_paths_by_inset)
            outline_summary << inset << ':' << count << ' ';
        BOOST_LOG_TRIVIAL(debug)
            << "Projected nonplanar Arachne outlines requested="
            << outline_count << " emitted=" << emit_outlines
            << " paths_by_inset=" << outline_summary.str();
        ExPolygons inner_contour = union_ex(wall_tool_paths.getInnerContour());
        if (!inner_contour.empty()) {
            // WallToolPaths::getInnerContour() is useful for detecting that
            // Arachne produced a viable wall set, but its contour may include
            // the outside side of an open/variable-width group.  Insetting
            // that result by only half a top-fill spacing let roof scanlines
            // begin on the outer wall and cross every requested wall loop.
            // Derive the fill ownership from the same deposited spacings as
            // the planar perimeter generator instead.  This intentionally
            // leaves a small bonding overlap, while making it impossible for
            // a top-fill centerline to penetrate the generated wall stack.
            // Derive fill ownership from the bead envelopes that actually
            // survived projection. Using outline_count here assumed that all
            // requested Arachne walls were emitted. At a mesh/feature
            // boundary an individual wall may be clipped or rejected while
            // the coherent surface remains valid; insetting for that missing
            // wall then leaves a real unextruded band. Start with the domain
            // in which a fill bead remains inside the model, then remove the
            // actual projected wall envelopes plus half a fill spacing. This
            // preserves the intended 0.02 mm bonding overlap and cannot cross
            // a retained wall merely because fewer walls fitted locally.
            centerline_area = offset_ex(
                reachable_projection, -coord_t(std::llround(0.5 * double(spacing))));
            const ExPolygons projected_wall_material =
                union_ex(reserved_wall_coverage);
            if (!projected_wall_material.empty()) {
                const coord_t wall_to_fill_clearance = std::max<coord_t>(
                    0, coord_t(std::llround(0.5 * double(spacing))) -
                           scale_(0.02));
                centerline_area = diff_ex(
                    centerline_area,
                    offset_ex(projected_wall_material, wall_to_fill_clearance),
                    ApplySafetyOffset::No);
            }
        }
        generated_arachne_outlines = true;
    }

    coord_t fill_inset = coord_t(std::llround(0.5 * double(spacing)));
    if (boundary_template != nullptr && !generated_arachne_outlines) {
        const ExPolygons model_outline_centerlines = offset_ex(
            patch.projection, -scale_(0.5 * boundary_width));
        const ExPolygons reachable_wall_domain = offset_ex(
            reachable_projection, scale_(0.01));
        for (const ExPolygon &polygon : model_outline_centerlines) {
            Polylines contours{polygon.contour.split_at_first_point()};
            for (const Polygon &hole : polygon.holes)
                contours.emplace_back(hole.split_at_first_point());
            const Polylines clipped_contours =
                intersection_pl(contours, reachable_wall_domain);
            for (const Polyline &fragment : clipped_contours)
                if (ExtrusionPath *path = project_line(
                        fragment.points, *boundary_template, false)) {
                    path->polygons_covered_by_width(
                        reserved_wall_coverage, 0.f);
                    if (emit_outlines)
                        output.emplace_back(path);
                    else
                        delete path;
                }
        }
        const coord_t boundary_spacing = Flow(
            boundary_template->width, boundary_template->height, 0.f).scaled_spacing();
        // The outline remains half its nominal width inside the model for
        // dimensional accuracy. Place the first top-skin centerline using the
        // two deposited spacings, not their nominal widths, then add a tiny
        // bonding overlap. Width-based placement left a visible unextruded
        // band at the non-planar outline; a larger overlap made the two bead
        // envelopes collide instead.
        fill_inset = scale_(0.5 * boundary_width) +
            (boundary_spacing + spacing) / 2 - scale_(0.02);
    }
    if (!generated_arachne_outlines)
        centerline_area = offset_ex(reachable_projection, -fill_inset);

    // A deferred chimney/boss course owns its own wall footprint.  Do not
    // run the roof raster underneath that wall: the two paths would be
    // emitted at nearly identical Z and their deposited beads would overlap.
    // Stop the top-skin centerline one half fill width before the feature's
    // modeled boundary, retaining 0.02 mm of bead-envelope contact for a
    // sealed interface.  The surface-following feature course then supplies
    // the structural wall without either a gap or duplicate extrusion.
    if (!solid_feature_footprint.empty()) {
        const coord_t feature_contact_margin = std::max<coord_t>(
            0, scale_(0.5 * double(path_template.width) - 0.02));
        centerline_area = diff_ex(
            centerline_area,
            offset_ex(solid_feature_footprint, feature_contact_margin),
            ApplySafetyOffset::No);
    }
    if (centerline_area.empty())
        return output;
    // Classify the patch from its area-weighted slope field. A plane or a
    // developable roof has gradient variation in at most one tangent
    // direction, and is best covered by rows along that principal direction.
    // A hill/dome varies materially in both directions; concentric offsets
    // follow that two-dimensional curvature without imposing a single global
    // staircase direction. An explicit concentric top-surface selection also
    // remains authoritative.
    const ExPolygons classification_interior = offset_ex(
        reachable_projection,
        -scale_(2. * std::max(double(path_template.width), resolution)));
    const auto use_for_classification = [&classification_interior](
                                            const WeightedPatchGradient &sample) {
        return classification_interior.empty() ||
               Geometry::contains(classification_interior, sample.centroid);
    };
    Vec2d mean_gradient = Vec2d::Zero();
    double gradient_weight = 0.;
    for (const WeightedPatchGradient &sample : reachable_gradients) {
        if (!use_for_classification(sample))
            continue;
        mean_gradient += sample.gradient * sample.weight;
        gradient_weight += sample.weight;
    }
    if (gradient_weight > 0.)
        mean_gradient /= gradient_weight;
    Eigen::Matrix2d gradient_covariance = Eigen::Matrix2d::Zero();
    for (const WeightedPatchGradient &sample : reachable_gradients) {
        if (!use_for_classification(sample))
            continue;
        const Vec2d centered = sample.gradient - mean_gradient;
        gradient_covariance += sample.weight * centered * centered.transpose();
    }
    if (gradient_weight > 0.)
        gradient_covariance /= gradient_weight;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> gradient_solver(
        gradient_covariance);
    const Vec2d gradient_variance = gradient_solver.eigenvalues();
    // Rounded edge/fillet facets on an otherwise developable roof contribute
    // a small second eigenvalue. Require both a meaningful absolute spread
    // and a strongly two-dimensional ratio so those boundary facets do not
    // turn a broad roof into a concentric hill. A dome has comparable
    // variance in both axes and remains comfortably above these thresholds.
    constexpr double minimum_two_axis_variance = 0.0025;
    const bool doubly_curved = gradient_variance.x() > minimum_two_axis_variance &&
        gradient_variance.x() > 0.35 * gradient_variance.y();
    cache.doubly_curved = doubly_curved;
    BOOST_LOG_TRIVIAL(debug)
        << "Nonplanar surface pattern eigenvalues=" << gradient_variance.transpose()
        << " mean_gradient=" << mean_gradient.transpose()
        << " doubly_curved=" << doubly_curved;
    const bool use_concentric = top_surface_pattern == ipConcentric || doubly_curved;

    Polylines clipped;
    double scanline_first_row = 0.;
    double scanline_row_spacing = unscale_(spacing);
    Vec2d row_direction = mean_gradient.norm() > 0.01 ?
        mean_gradient.normalized() :
        gradient_solver.eigenvectors().col(1).normalized();
    if (!std::isfinite(row_direction.x()) || !std::isfinite(row_direction.y()) ||
        row_direction.squaredNorm() < 0.5)
        row_direction = Vec2d::UnitX();
    if (row_direction.x() < -EPSILON ||
        (std::abs(row_direction.x()) <= EPSILON && row_direction.y() < 0.))
        row_direction = -row_direction;
    const Vec2d row_normal(-row_direction.y(), row_direction.x());

    if (use_concentric) {
        ExPolygons inset = centerline_area;
        while (!inset.empty()) {
            append(clipped, to_polylines(inset));
            inset = offset_ex(inset, -spacing);
        }
    } else {
        const BoundingBox bbox = get_extents(centerline_area);
        const std::array<Vec2d, 4> corners{
            unscale(bbox.min),
            Vec2d(unscale_(bbox.max.x()), unscale_(bbox.min.y())),
            unscale(bbox.max),
            Vec2d(unscale_(bbox.min.x()), unscale_(bbox.max.y()))};
        double minimum_row = std::numeric_limits<double>::max();
        double maximum_row = std::numeric_limits<double>::lowest();
        double minimum_along = std::numeric_limits<double>::max();
        double maximum_along = std::numeric_limits<double>::lowest();
        for (const Vec2d &corner : corners) {
            minimum_row = std::min(minimum_row, row_normal.dot(corner));
            maximum_row = std::max(maximum_row, row_normal.dot(corner));
            minimum_along = std::min(minimum_along, row_direction.dot(corner));
            maximum_along = std::max(maximum_along, row_direction.dot(corner));
        }
        const double spacing_mm = unscale_(spacing);
        // Fit the scanline lattice to both sides of the printable centerline
        // span. Merely centering a nominal-spacing lattice still leaves up to
        // half a spacing *in addition to* the wall-to-fill inset already
        // represented by centerline_area. On a sloped eave that produced a
        // real wedge between the innermost draped wall and the first skin row
        // (the Benchy roof is a representative case). Use one more interval
        // when necessary and distribute it evenly, so spacing never exceeds
        // the requested value and both edge rows remain bonded to the wall.
        // Keep the rows a tiny distance inside the clipping boundary to avoid
        // zero-length tangencies at sharp corners.
        const double row_span = std::max(0., maximum_row - minimum_row);
        const double edge_inset = std::min(0.02, 0.5 * row_span);
        const double usable_row_span = std::max(0., row_span - 2. * edge_inset);
        const size_t interval_count = usable_row_span <= EPSILON ? 0 :
            std::max<size_t>(1, size_t(std::ceil(usable_row_span / spacing_mm)));
        const size_t row_count = interval_count + 1;
        const double fitted_spacing = interval_count == 0 ? 0. :
            usable_row_span / double(interval_count);
        const double first_row = minimum_row + edge_inset;
        scanline_first_row = first_row;
        scanline_row_spacing = fitted_spacing;
        Polylines scanlines;
        scanlines.reserve(row_count);
        for (size_t row_idx = 0; row_idx < row_count; ++row_idx) {
            const double row = first_row + fitted_spacing * double(row_idx);
            const Vec2d first = row_normal * row +
                row_direction * (minimum_along - spacing_mm);
            const Vec2d last = row_normal * row +
                row_direction * (maximum_along + spacing_mm);
            scanlines.emplace_back(Point(scale_(first.x()), scale_(first.y())),
                                   Point(scale_(last.x()), scale_(last.y())));
        }
        clipped = intersection_pl(scanlines, centerline_area);
    }

    // Clipper does not promise spatial ordering for the fragments it returns.
    // Order complete scanline rows before alternating them; otherwise the
    // connector pass sees unrelated roof fragments as neighbors and cannot
    // form continuous serpentine paths.
    const auto row_index = [&row_normal, scanline_first_row,
                            scanline_row_spacing](const Polyline &line) {
        const Vec2d point = unscale(line.points.front());
        return scanline_row_spacing <= EPSILON ? int64_t(0) :
            std::llround((row_normal.dot(point) - scanline_first_row) /
                         scanline_row_spacing);
    };
    const auto line_along_min = [&row_direction](const Polyline &line) {
        return std::min(row_direction.dot(unscale(line.points.front())),
                        row_direction.dot(unscale(line.points.back())));
    };
    if (!use_concentric)
        std::sort(clipped.begin(), clipped.end(), [&row_index, &line_along_min](
                  const Polyline &left, const Polyline &right) {
            const int64_t left_row = row_index(left);
            const int64_t right_row = row_index(right);
            return left_row < right_row ||
                   (left_row == right_row && line_along_min(left) < line_along_min(right));
        });
    bool reverse_row = false;
    for (size_t row_begin = 0; !use_concentric && row_begin < clipped.size();) {
        size_t row_end = row_begin + 1;
        const int64_t current_row = row_index(clipped[row_begin]);
        while (row_end < clipped.size() && row_index(clipped[row_end]) == current_row)
            ++row_end;
        if (reverse_row) {
            std::reverse(clipped.begin() + row_begin, clipped.begin() + row_end);
            for (size_t idx = row_begin; idx < row_end; ++idx)
                clipped[idx].reverse();
        }
        reverse_row = !reverse_row;
        row_begin = row_end;
    }

    for (Polyline &line : clipped) {
        if (line.points.size() < 2)
            continue;
        Polyline3 projected;
        bool valid = true;
        for (size_t segment_idx = 0; segment_idx + 1 < line.points.size(); ++segment_idx) {
            const Vec2d a = unscale(line.points[segment_idx]);
            const Vec2d b = unscale(line.points[segment_idx + 1]);
            const Vec2d delta = b - a;
            const int count = std::max(1, int(std::ceil(delta.norm() / resolution)));
            for (int sample = segment_idx == 0 ? 0 : 1; sample <= count; ++sample) {
                const Vec2d xy = a + delta * (double(sample) / count);
                const Point sample_point(scale_(xy.x()), scale_(xy.y()));
                const std::optional<double> surface_z =
                    patch_surface_z_at_xy(patch, mesh, sample_point);
                if (!surface_z) {
                    valid = false;
                    break;
                }
                projected.append(Point3(scale_(xy.x()), scale_(xy.y()),
                    scale_(*surface_z - ground - owner_z - surface_depth)));
            }
            if (!valid)
                break;
        }
        if (!valid || projected.points.size() < 2)
            continue;
        projected.simplify(scale_(std::min(0.005, resolution * 0.025)));
        if (projected.points.size() < 2)
            continue;
        ExtrusionPath path(std::move(projected), path_template);
        path.z_contoured = true;
        path.nonplanar_surface = true;
        path.nonplanar_schedule_owned = false;
        path.nonplanar_transition = surface_depth > EPSILON;
        path.set_reverse();
        // Alternate neighboring scanlines so short boundary connectors can be
        // retained as a continuous serpentine after clearance validation.
        // Every shell depth is supported by the preceding draped shell, so an
        // ascending line is no less supported than its descending neighbor.
        output.emplace_back(path.clone_move());
    }

    // Clipping the regular raster against an Arachne inner contour aligns
    // centerlines, not the envelopes of the deposited beads.  At a curved or
    // variable-width wall the round cap of a clipped raster line may therefore
    // leave a small triangular opening even though the nominal centerline
    // overlap is correct.  A larger blanket overlap is not safe here: it makes
    // the same paths cross on tighter curves and around deferred features.
    //
    // Measure the coverage of the paths we actually generated and repair only
    // narrow residuals.  This is the same medial-axis/variable-width principle
    // used by planar Arachne gap fill, but every repair path is projected back
    // onto the selected surface below.  Consequently the repair bonds the top
    // skin to its walls without filling modeled openings or adding a second
    // full-density raster.
    Polygons deposited;
    for (const ExtrusionEntity *entity : output)
        entity->polygons_covered_by_width(deposited, 0.f);
    append(deposited, reserved_wall_coverage);

    ExPolygons repair_target = offset_ex(reachable_projection, -scale_(0.01));
    if (!solid_feature_footprint.empty()) {
        // The deferred feature course owns this contact band.  Its walls are
        // emitted later in dependency order, so treating the band as a roof
        // gap would put a second extrusion underneath the feature.
        const coord_t feature_contact_margin = std::max<coord_t>(
            0, scale_(0.5 * double(path_template.width) - 0.02));
        repair_target = diff_ex(
            repair_target,
            offset_ex(solid_feature_footprint, feature_contact_margin),
            ApplySafetyOffset::No);
    }

    ExPolygons uncovered = diff_ex(
        repair_target, union_ex(deposited), ApplySafetyOffset::No);
    if (!uncovered.empty()) {
        const double minimum_gap_width =
            0.2 * double(spacing) * (1. - INSET_OVERLAP_TOLERANCE);
        const double maximum_gap_width = 2. * double(spacing);
        // Feed the actual uncovered domain to Arachne's medial-axis filter.
        // The former opening/difference pre-filter erased tapered wall-to-skin
        // bands before medial_axis() could classify their local width.  On a
        // sloped Benchy roof that left a 0.2--0.3 mm strip with no bead at the
        // eave even though both neighboring paths were otherwise valid.  The
        // medial-axis width limits below already reject numerical slivers and
        // broad regions, and unlike a global morphological opening they retain
        // a band whose width varies along a curved/variable-width wall.
        ExPolygons narrow_gaps = std::move(uncovered);

        ThickPolylines medial_lines;
        for (ExPolygon &gap : narrow_gaps) {
            gap.douglas_peucker(SCALED_RESOLUTION * 0.1);
            gap.medial_axis(minimum_gap_width, maximum_gap_width, &medial_lines);
        }

        // Short wedge-shaped anchors are useful here, unlike optional planar
        // gap fill where very short paths are commonly filtered for cosmetic
        // reasons.  Reject only numerical remnants that cannot produce a
        // meaningful motion at the configured projection resolution.
        const double minimum_anchor_length = scale_(
            std::max(0.05, std::min(0.2, resolution)));
        medial_lines.erase(
            std::remove_if(medial_lines.begin(), medial_lines.end(),
                [minimum_anchor_length](const ThickPolyline &line) {
                    return line.length() < minimum_anchor_length;
                }),
            medial_lines.end());

        ExtrusionEntitiesPtr planar_repairs;
        const Flow repair_flow(
            path_template.width, path_template.height, 0.f);
        variable_width(
            medial_lines, path_template.role(), repair_flow, planar_repairs);
        size_t projected_repairs = 0;
        for (ExtrusionEntity *entity : planar_repairs) {
            const auto append_projected = [&](const ExtrusionPath &planar_path) {
                if (ExtrusionPath *projected = project_line(
                        to_points(planar_path.polyline.points),
                        planar_path, false)) {
                    output.emplace_back(projected);
                    ++projected_repairs;
                }
            };
            if (const auto *path = dynamic_cast<const ExtrusionPath *>(entity)) {
                append_projected(*path);
            } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(entity)) {
                for (const ExtrusionPath &path : loop->paths)
                    append_projected(path);
            }
            delete entity;
        }
        BOOST_LOG_TRIVIAL(debug)
            << "Added " << projected_repairs
            << " variable-width nonplanar surface anchors";
    }
    return output;
}

std::vector<Vec3d> clip_vertices_at_z(const std::vector<Vec3d> &input,
                                      double limit, bool keep_above)
{
    std::vector<Vec3d> output;
    if (input.empty())
        return output;
    auto inside = [limit, keep_above](const Vec3d &point) {
        return keep_above ? point.z() + EPSILON >= limit : point.z() <= limit + EPSILON;
    };
    Vec3d previous = input.back();
    bool previous_inside = inside(previous);
    for (const Vec3d &current : input) {
        const bool current_inside = inside(current);
        if (current_inside != previous_inside) {
            const double dz = current.z() - previous.z();
            if (std::abs(dz) > EPSILON) {
                const double ratio = std::clamp((limit - previous.z()) / dz, 0., 1.);
                output.emplace_back(previous + ratio * (current - previous));
            }
        }
        if (current_inside)
            output.emplace_back(current);
        previous = current;
        previous_inside = current_inside;
    }
    return output;
}

ExPolygons patch_shell_band_at_layer(const NonplanarPatch &patch,
                                     const sla::IndexedMesh &mesh,
                                     double layer_z, double shell_depth,
                                     double layer_height,
                                     double upper_clearance = 0.)
{
    const indexed_triangle_set &its = *mesh.get_triangle_mesh();
    const double ground = mesh.ground_level();
    // Conventional paths are removed only where this nominal Z plane passes
    // through the replacement shell. This avoids the previous full-column
    // cutout, which deleted solid support below shallow surfaces.
    // Apply the same full-course clearance below the nominal plane. On a
    // descending part of the shell a generated transition may sit just below
    // the next conventional course (for example 37.50 mm versus a 37.60 mm
    // bridge). A half-layer lower bound retained that course and let it cross
    // the draped path from above. The course a complete layer below remains as
    // physical support.
    const double lower_surface_z = ground + layer_z - layer_height;
    // Keep the conventional course immediately below the lowest draped
    // shell. Extending this upper bound by another complete layer selected
    // precisely that support course: after it was clipped away, the first
    // generated roof pass could be 0.4 mm (and at an outside stair, much
    // farther) above the nearest previously deposited bead. The lower bound
    // still removes a stair course whose bead reaches into a descending
    // generated path; the upper bound describes only the requested shell
    // depth, leaving its dependency in place.
    const double upper_surface_z =
        ground + layer_z + shell_depth + upper_clearance;
    // Commit visits every potentially affected owner layer so it can remove
    // conventional paths crossing an accepted shell.  Most of those planes
    // are well outside this patch's mesh Z extent.  Reject them before
    // clipping any triangles or invoking Clipper; this is an exact interval
    // test and therefore cannot change the generated replacement area.
    if (upper_surface_z + EPSILON < patch.minimum_z ||
        lower_surface_z > patch.maximum_z + EPSILON)
        return {};
    Polygons polygons;
    for (const size_t face : patch.face_indices) {
        const Vec3i32 &indices = its.indices[face];
        std::vector<Vec3d> vertices{
            its.vertices[size_t(indices[0])].cast<double>(),
            its.vertices[size_t(indices[1])].cast<double>(),
            its.vertices[size_t(indices[2])].cast<double>()};
        vertices = clip_vertices_at_z(vertices, lower_surface_z, true);
        vertices = clip_vertices_at_z(vertices, upper_surface_z, false);
        if (vertices.size() < 3)
            continue;
        Polygon polygon;
        polygon.points.reserve(vertices.size());
        for (const Vec3d &vertex : vertices)
            polygon.points.emplace_back(scale_(vertex.x()), scale_(vertex.y()));
        if (polygon.area() < 0)
            polygon.reverse();
        if (polygon.is_valid())
            polygons.emplace_back(std::move(polygon));
    }
    return union_ex(polygons);
}

void collect_generated_material_at_z(
    const ExtrusionEntity &entity, double owner_z, double layer_z,
    double z_tolerance, bool include_planar_foundation, Polygons &coverage,
    bool round_caps = false)
{
    const auto collect_path = [&](const ExtrusionPath &path) {
        if ((!path.nonplanar_surface &&
             !(include_planar_foundation &&
               path.nonplanar_clearance_validated)) ||
            path.polyline.points.size() < 2)
            return;
        for (size_t point_idx = 1;
             point_idx < path.polyline.points.size(); ++point_idx) {
            const Point3 &first = path.polyline.points[point_idx - 1];
            const Point3 &last = path.polyline.points[point_idx];
            const double first_z = owner_z + unscale_(first.z());
            const double last_z = owner_z + unscale_(last.z());
            const double dz = last_z - first_z;
            double begin = 0.;
            double end = 1.;
            if (std::abs(dz) <= EPSILON) {
                if (std::abs(first_z - layer_z) > z_tolerance)
                    continue;
            } else {
                const double first_t =
                    (layer_z - z_tolerance - first_z) / dz;
                const double last_t =
                    (layer_z + z_tolerance - first_z) / dz;
                begin = std::max(0., std::min(first_t, last_t));
                end = std::min(1., std::max(first_t, last_t));
                if (end <= begin + EPSILON)
                    continue;
            }
            const auto interpolate = [&](double ratio) {
                return Point3(
                    coord_t(std::llround(
                        double(first.x()) + ratio * double(last.x() - first.x()))),
                    coord_t(std::llround(
                        double(first.y()) + ratio * double(last.y() - first.y()))),
                    coord_t(0));
            };
            ExtrusionPath clipped(path);
            clipped.polyline.points = {interpolate(begin), interpolate(end)};
            if (clipped.polyline.length() > SCALED_EPSILON) {
                if (round_caps)
                    polygons_append(coverage, offset(clipped.polyline.to_polyline(),
                        float(scale_(path.width / 2.)), ClipperLib::jtRound,
                        scale_(0.005), ClipperLib::etOpenRound));
                else
                    clipped.polygons_covered_by_width(coverage, 0.f);
            }
        }
    };

    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity)) {
        collect_path(*path);
    } else if (const auto *multipath =
                   dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
        for (const ExtrusionPath &path : multipath->paths)
            collect_path(path);
    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
        for (const ExtrusionPath &path : loop->paths)
            collect_path(path);
    } else if (const auto *collection =
                   dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        for (const ExtrusionEntity *child : collection->entities)
            collect_generated_material_at_z(
                *child, owner_z, layer_z, z_tolerance,
                include_planar_foundation, coverage, round_caps);
    }
}

ExPolygons generated_material_at_z(
    const std::vector<ExtrusionEntitiesPtr> &generated_by_layer,
    const std::vector<Layer *> &layers, double layer_z, double z_tolerance,
    bool include_planar_foundation = false)
{
    Polygons coverage;
    for (size_t owner = 0;
         owner < generated_by_layer.size() && owner < layers.size(); ++owner)
        for (const ExtrusionEntity *entity : generated_by_layer[owner])
            collect_generated_material_at_z(
                *entity, layers[owner]->print_z, layer_z, z_tolerance,
                include_planar_foundation, coverage);
    return union_ex(coverage);
}

void collect_deposited_coverage(const ExtrusionEntity &entity,
                                Polygons &coverage)
{
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity)) {
        path->polygons_covered_by_width(coverage, 0.f);
    } else if (const auto *multipath =
                   dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
        for (const ExtrusionPath &path : multipath->paths)
            path.polygons_covered_by_width(coverage, 0.f);
    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
        for (const ExtrusionPath &path : loop->paths)
            path.polygons_covered_by_width(coverage, 0.f);
    } else if (const auto *collection =
                   dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        for (const ExtrusionEntity *child : collection->entities)
            collect_deposited_coverage(*child, coverage);
    }
}

bool path_has_supported_run(const ExtrusionPath &path,
                            const ExPolygons &support)
{
    if (path.polyline.points.size() < 2 || support.empty())
        return false;
    const ExPolygons contact = offset_ex(
        support, scale_(0.5 * double(path.width) + 0.05));
    const double maximum_unsupported_run =
        std::max(0.8, 2. * double(path.width));
    double unsupported_run = 0.;
    for (size_t idx = 1; idx < path.polyline.points.size(); ++idx) {
        const Point3 &first = path.polyline.points[idx - 1];
        const Point3 &last = path.polyline.points[idx];
        const Point midpoint(
            coord_t((int64_t(first.x()) + int64_t(last.x())) / 2),
            coord_t((int64_t(first.y()) + int64_t(last.y())) / 2));
        if (Geometry::contains(contact, midpoint)) {
            unsupported_run = 0.;
            continue;
        }
        unsupported_run += unscale_(
            (last.to_point() - first.to_point()).cast<double>().norm());
        if (unsupported_run > maximum_unsupported_run)
            return false;
    }
    return true;
}

void collect_extrusion_path_references(
    const ExtrusionEntity &entity,
    std::vector<const ExtrusionPath *> &paths)
{
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity)) {
        if (path->polyline.points.size() >= 2)
            paths.emplace_back(path);
    } else if (const auto *multipath =
                   dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
        for (const ExtrusionPath &path : multipath->paths)
            collect_extrusion_path_references(path, paths);
    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
        for (const ExtrusionPath &path : loop->paths)
            collect_extrusion_path_references(path, paths);
    } else if (const auto *collection =
                   dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        for (const ExtrusionEntity *child : collection->entities)
            collect_extrusion_path_references(*child, paths);
    }
}

} // namespace

static bool nonplanar_paths_are_supported(
    const ExtrusionEntity &entity,
    const std::vector<const ExtrusionPath *> &support_paths,
    double entity_owner_z, double support_owner_z,
    bool allow_short_unsupported_runs, bool walls_only)
{
    double unsupported_run = 0.;
    std::optional<Point3> previous_endpoint;
    const auto path_is_supported = [
        &support_paths, entity_owner_z,
        support_owner_z,
        allow_short_unsupported_runs, walls_only, &unsupported_run,
        &previous_endpoint](const ExtrusionPath &path) {
        // Solid foundation raster may bridge ordinary sparse infill just like
        // conventional internal-solid/top fill. The failure being guarded
        // here is a generated structural wall with no wall or solid material
        // under it; applying bead-by-bead support to raster fill would reject
        // every legitimate roof over sparse infill.
        if (walls_only && !is_perimeter(path.role()))
            return true;
        if (path.polyline.points.size() < 2 || support_paths.empty())
            return false;
        if (!previous_endpoint || *previous_endpoint != path.first_point3())
            unsupported_run = 0.;
        previous_endpoint = path.last_point3();
        // This is the foundation of a generated replacement stack, not an
        // ordinary overhang wall. Require every deposited segment to contact
        // the preceding course. Arachne represents one loop as many short
        // paths; applying the normal short-overhang allowance independently
        // to each fragment allowed all of them to pass and routing later
        // joined them into a several-millimetre floating ring. Validate the
        // actual bead envelope here, before fragments can be chained.
        const double maximum_unsupported_run =
            std::max(0.8, 2. * double(path.width));
        for (size_t idx = 1; idx < path.polyline.points.size(); ++idx) {
            const Point3 &first = path.polyline.points[idx - 1];
            const Point3 &last = path.polyline.points[idx];
            // Test the entire segment, not only its midpoint. A short support
            // crossing the middle of a long wall must not support both ends.
            // Shift the query into the middle of the permitted vertical gap
            // so the shared interval solver can enforce asymmetric Z bounds.
            const double minimum_gap = 0.04;
            const double maximum_gap = double(path.height) + 0.06;
            const double gap_center = 0.5 * (minimum_gap + maximum_gap);
            const double gap_tolerance = 0.5 * (maximum_gap - minimum_gap);
            Vec3d a(unscale_(first.x()), unscale_(first.y()),
                    entity_owner_z + unscale_(first.z()) - gap_center);
            Vec3d b(unscale_(last.x()), unscale_(last.y()),
                    entity_owner_z + unscale_(last.z()) - gap_center);
            std::vector<std::pair<double, double>> contact;
            bool complete_contact = false;
            for (const ExtrusionPath *support : support_paths) {
                const double contact_distance = 0.5 *
                    (double(path.width) + double(support->width)) + 0.02;
                const Points3 &points = support->polyline.points;
                for (size_t support_idx = 1; support_idx < points.size(); ++support_idx) {
                    const Point3 &start = points[support_idx - 1];
                    const Point3 &end = points[support_idx];
                    auto intervals = nonplanar_segment_overlap(a, b,
                        Vec3d(unscale_(start.x()), unscale_(start.y()),
                              support_owner_z + unscale_(start.z())),
                        Vec3d(unscale_(end.x()), unscale_(end.y()),
                              support_owner_z + unscale_(end.z())),
                        contact_distance, gap_tolerance);
                    complete_contact = intervals.size() == 1 &&
                        intervals.front().first == 0. && intervals.front().second == 1.;
                    if (complete_contact)
                        break;
                    append(contact, intervals);
                }
                if (complete_contact)
                    break;
            }
            if (complete_contact) {
                unsupported_run = 0.;
                continue;
            }
            std::sort(contact.begin(), contact.end());
            const double length = (b.head<2>() - a.head<2>()).norm();
            double covered_until = 0.;
            const auto gap_is_allowed = [&](double fraction) {
                unsupported_run += fraction * length;
                return allow_short_unsupported_runs ?
                    unsupported_run <= maximum_unsupported_run + EPSILON :
                    unsupported_run <= EPSILON;
            };
            for (const auto &[begin, end] : contact) {
                if (end <= covered_until)
                    continue;
                if (begin > covered_until && !gap_is_allowed(begin - covered_until))
                    return false;
                unsupported_run = 0.;
                covered_until = end;
            }
            if (!gap_is_allowed(1. - covered_until))
                return false;
        }
        return true;
    };

    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity))
        return path_is_supported(*path);
    if (const auto *multipath =
            dynamic_cast<const ExtrusionMultiPath *>(&entity))
        return std::all_of(
            multipath->paths.begin(), multipath->paths.end(),
            path_is_supported);
    if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity))
        return std::all_of(loop->paths.begin(), loop->paths.end(),
                           path_is_supported);
    if (const auto *collection =
            dynamic_cast<const ExtrusionEntityCollection *>(&entity))
        return std::all_of(
            collection->entities.begin(), collection->entities.end(),
            [&support_paths, entity_owner_z,
             support_owner_z,
             allow_short_unsupported_runs, walls_only](const ExtrusionEntity *child) {
                return child != nullptr &&
                    nonplanar_paths_are_supported(
                        *child, support_paths,
                        entity_owner_z, support_owner_z,
                        allow_short_unsupported_runs, walls_only);
            });
    return false;
}

bool nonplanar_walls_are_supported(
    const ExtrusionEntity &entity,
    const std::vector<const ExtrusionPath *> &support_paths,
    double entity_owner_z, double support_owner_z,
    bool allow_short_unsupported_runs)
{
    return nonplanar_paths_are_supported(entity, support_paths,
        entity_owner_z, support_owner_z, allow_short_unsupported_runs, true);
}

bool nonplanar_course_is_supported(
    const ExtrusionEntity &entity,
    const std::vector<const ExtrusionPath *> &support_paths,
    double entity_owner_z, double support_owner_z)
{
    return nonplanar_paths_are_supported(entity, support_paths,
        entity_owner_z, support_owner_z, false, false);
}

bool native_foundation_supports_course(
    const ExtrusionEntitiesPtr &course, const LayerRegion &foundation,
    double owner_z, double foundation_z, double minimum_height, double maximum_height)
{
    std::vector<const ExtrusionPath *> support_paths;
    collect_extrusion_path_references(foundation.perimeters, support_paths);
    collect_extrusion_path_references(foundation.fills, support_paths);
    // A raised brick/interlocking path is an obstacle, not simply a path to
    // omit from the support list. Reuse requires a genuinely flat native
    // course; mixed-height courses need an explicit realignment first.
    if (std::any_of(support_paths.begin(), support_paths.end(),
        [](const ExtrusionPath *path) {
            return path->nonplanar_surface || std::any_of(
                path->polyline.points.begin(), path->polyline.points.end(),
                [](const Point3 &point) { return point.z() != 0; });
        }))
        return false;
    if (course.empty() || support_paths.empty())
        return false;
    for (const ExtrusionEntity *entity : course) {
        if (entity == nullptr)
            return false;
        std::vector<const ExtrusionPath *> paths;
        collect_extrusion_path_references(*entity, paths);
        if (paths.empty())
            return false;
        for (const ExtrusionPath *path : paths) {
            if (path->polyline.points.size() < 2)
                return false;
            ExtrusionPath segment(*path);
            for (size_t i = 1; i < path->polyline.points.size(); ++i) {
                const Point3 &a = path->polyline.points[i - 1];
                const Point3 &b = path->polyline.points[i];
                const double low = owner_z + unscale_(std::min(a.z(), b.z())) - foundation_z;
                const double high = owner_z + unscale_(std::max(a.z(), b.z())) - foundation_z;
                if (low < minimum_height - EPSILON || high > maximum_height + EPSILON)
                    return false;
                segment.polyline.points = {a, b};
                segment.height = low;
                // The foundation was verified flat above. The first course's
                // local thickness is therefore z(s) - foundation_z, not the
                // minimum thickness held constant along the whole segment.
                // Check its contact footprint on that plane after bounding
                // both endpoint heights. This changes only the validation
                // probe; course_flow emits the actual varying thickness.
                for (Point3 &point : segment.polyline.points)
                    point.z() = 0;
                if (!nonplanar_course_is_supported(segment, support_paths,
                        foundation_z + low, foundation_z))
                    return false;
            }
        }
    }
    return true;
}

namespace {

ExPolygons retained_nonplanar_wall_material_at_z(
    const std::vector<Layer *> &layers, size_t region_idx, double layer_z,
    double z_tolerance)
{
    Polygons coverage;
    for (const Layer *layer : layers) {
        if (region_idx >= layer->regions().size())
            continue;
        collect_generated_material_at_z(
            layer->regions()[region_idx]->perimeters,
            layer->print_z, layer_z, z_tolerance, false, coverage);
    }
    return union_ex(coverage);
}

size_t owner_layer_for_z(const std::vector<Layer *> &layers, size_t first_layer,
                         size_t last_layer, double required_z)
{
    size_t owner = std::min(first_layer, layers.size() - 1);
    while (owner < layers.size() && layers[owner]->print_z + EPSILON < required_z)
        ++owner;
    return std::min(owner, std::min(last_layer, layers.size() - 1));
}

std::vector<PrintedToolpathObstacle> printed_model_obstacles(
    const std::vector<Layer *> &layers, size_t last_layer,
    const ExPolygons *local_surface = nullptr)
{
    std::vector<PrintedToolpathObstacle> obstacles;
    obstacles.reserve(last_layer + 1);
    for (size_t layer_idx = 0; layer_idx <= last_layer && layer_idx < layers.size(); ++layer_idx) {
        ExPolygons footprint_ex = layers[layer_idx]->lslices;
        if (local_surface != nullptr)
            footprint_ex = diff_ex(footprint_ex, *local_surface, ApplySafetyOffset::No);
        Polygons footprint = to_polygons(footprint_ex);
        if (!footprint.empty())
            obstacles.emplace_back(layers[layer_idx]->print_z, std::move(footprint));
    }
    return obstacles;
}

template<class Obstacles>
bool generated_path_is_clear(
    ExtrusionPath &path, coordf_t owner_z,
    const SequentialGantryGeometry &gantry,
    const ToolheadClearanceFallback &fallback,
    const Obstacles &remote_obstacles,
    std::string *failure = nullptr)
{
    std::vector<Vec3d> points;
    points.reserve(path.polyline.points.size());
    for (const Point3 &point : path.polyline.points)
        points.emplace_back(unscale_(point.x()), unscale_(point.y()),
                            owner_z + unscale_(point.z()));
    if (points.size() < 2)
        return false;

    for (size_t idx = 1; idx < points.size(); ++idx) {
        double support_contact = std::max(1e-4, double(path.height));
        const double xy_distance = (points[idx].head<2>() - points[idx - 1].head<2>()).norm();
        if (xy_distance > EPSILON) {
            const double slope = std::abs(points[idx].z() - points[idx - 1].z()) / xy_distance;
            const double nozzle_reach = gantry.empty() ? fallback.radius :
                gantry.maximum_xy_reach(0.001);
            support_contact += (nozzle_reach + 0.5 * double(path.width)) * slope;
        }
        // Layer slices represent the whole local patch at one nominal Z. On a
        // shallow surface that overstates the height of its staircase by
        // several layers and makes every valid draped line collide with its
        // own support. The raster is already constrained to the selected mesh
        // facets and generated at every requested support depth, so clearance
        // envelopes are checked against geometry outside that connected patch.
        // Nearby chimneys and walls remain in remote_obstacles and are rejected.
        // Detailed gantry geometry contains the nozzle envelope at height
        // zero.  Running the fallback cone first therefore repeats the most
        // expensive segment/obstacle query and may reject paths that the
        // printer-specific geometry explicitly permits.  An empty geometry
        // still selects the conservative fallback inside the shared checker.
        const ToolheadCollision collision = check_toolhead_clearance(
                points[idx - 1], points[idx], gantry, remote_obstacles,
                fallback, support_contact);
        if (collision) {
            if (failure != nullptr) {
                std::ostringstream message;
                message << "segment=[" << points[idx - 1].transpose()
                        << "]-[" << points[idx].transpose()
                        << "] obstacle=" << collision.obstacle_index
                        << " gantry_slice=" << collision.slice_index
                        << " contact_allowance=" << support_contact;
                *failure = message.str();
            }
            return false;
        }
    }
    path.nonplanar_clearance_validated = true;
    return true;
}

template<class Obstacles>
const PrintedToolpathObstacle &obstacle_at(const Obstacles &obstacles,
                                           size_t index)
{
    using Element = std::decay_t<decltype(obstacles[index])>;
    if constexpr (std::is_pointer_v<Element>)
        return *obstacles[index];
    else
        return obstacles[index];
}

template<class Obstacles>
bool validate_deferred_feature(
    ExtrusionEntity &entity, coordf_t owner_z,
    const SequentialGantryGeometry &gantry,
    const ToolheadClearanceFallback &fallback,
    const Obstacles &remote_obstacles)
{
    auto validate_path = [&](ExtrusionPath &path) {
        if (generated_path_is_clear(
                path, owner_z, gantry, fallback, remote_obstacles))
            return true;
        ExtrusionPath nozzle_probe(path);
        const SequentialGantryGeometry nozzle_only;
        const bool nozzle_clear = generated_path_is_clear(
            nozzle_probe, owner_z, nozzle_only, fallback, remote_obstacles);
        coord_t minimum_z = std::numeric_limits<coord_t>::max();
        coord_t maximum_z = std::numeric_limits<coord_t>::lowest();
        for (const Point3 &point : path.polyline.points) {
            minimum_z = std::min(minimum_z, point.z());
            maximum_z = std::max(maximum_z, point.z());
        }
        BOOST_LOG_TRIVIAL(debug)
            << "Deferred feature clearance rejected "
            << (nozzle_clear ? "by gantry" : "by nozzle")
            << " at owner_z=" << owner_z
            << " path_z=[" << unscale_(minimum_z) << ','
            << unscale_(maximum_z) << "]";
        if (!nozzle_clear) {
            std::vector<Vec3d> points;
            points.reserve(path.polyline.points.size());
            for (const Point3 &point : path.polyline.points)
                points.emplace_back(unscale_(point.x()), unscale_(point.y()),
                                    owner_z + unscale_(point.z()));
            for (size_t idx = 1; idx < points.size(); ++idx) {
                const SequentialGantryGeometry nozzle_only;
                double support_contact = std::max(1e-4, double(path.height));
                const double xy_distance =
                    (points[idx].head<2>() - points[idx - 1].head<2>()).norm();
                if (xy_distance > EPSILON) {
                    const double slope = std::abs(
                        points[idx].z() - points[idx - 1].z()) / xy_distance;
                    const double nozzle_reach = gantry.empty() ? fallback.radius :
                        gantry.maximum_xy_reach(0.001);
                    support_contact +=
                        (nozzle_reach + 0.5 * double(path.width)) * slope;
                }
                const ToolheadCollision collision = check_toolhead_clearance(
                    points[idx - 1], points[idx], nozzle_only,
                    remote_obstacles, fallback, support_contact);
                if (collision) {
                    const PrintedToolpathObstacle &obstacle =
                        obstacle_at(remote_obstacles, collision.obstacle_index);
                    BOOST_LOG_TRIVIAL(debug)
                        << "Deferred feature nozzle obstacle "
                        << collision.obstacle_index << " top_z="
                        << obstacle.top_z
                        << " segment_z=[" << points[idx - 1].z() << ','
                        << points[idx].z() << "] segment_xy=["
                        << points[idx - 1].x() << ',' << points[idx - 1].y()
                        << " -> " << points[idx].x() << ',' << points[idx].y()
                        << "] obstacle_bbox=["
                        << unscale_(obstacle.bbox.min.x())
                        << ','
                        << unscale_(obstacle.bbox.min.y())
                        << " -> "
                        << unscale_(obstacle.bbox.max.x())
                        << ','
                        << unscale_(obstacle.bbox.max.y())
                        << ']';
                    break;
                }
            }
        }
        return false;
    };
    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity))
        return validate_path(*path);
    if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity))
        return std::all_of(
            multipath->paths.begin(), multipath->paths.end(),
            [&](ExtrusionPath &path) {
                return validate_path(path);
            });
    if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity))
        return std::all_of(
            loop->paths.begin(), loop->paths.end(),
            [&](ExtrusionPath &path) {
                return validate_path(path);
            });
    if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity))
        return std::all_of(
            collection->entities.begin(), collection->entities.end(),
            [&](ExtrusionEntity *child) {
                return validate_deferred_feature(
                    *child, owner_z, gantry, fallback, remote_obstacles);
            });
    return false;
}

} // namespace

std::vector<NonplanarCourseProfile> plan_nonplanar_courses(
    double anchor_z, double minimum_surface_z, double maximum_surface_z,
    size_t minimum_courses, double minimum_height, double maximum_height)
{
    const double rise = minimum_surface_z - anchor_z;
    const double span = maximum_surface_z - minimum_surface_z;
    if (!std::isfinite(rise) || !std::isfinite(span) || rise <= 0. || span < 0. ||
        !std::isfinite(minimum_height) || !std::isfinite(maximum_height) ||
        minimum_height <= 0. || maximum_height < minimum_height)
        return {};
    // smoothstep's maximum derivative is 1.5. This conservative bound
    // covers every XY height in the patch, not just its highest point.
    const double required_steps = std::ceil((rise + 1.5 * span) / maximum_height);
    const double available_steps = std::floor(rise / minimum_height);
    if (required_steps > available_steps ||
        double(minimum_courses) > available_steps ||
        required_steps >= double(std::numeric_limits<size_t>::max()))
        return {};
    const size_t count = std::max<size_t>({1, minimum_courses, size_t(required_steps)});
    std::vector<NonplanarCourseProfile> result;
    if (count >= result.max_size())
        return {};
    result.reserve(count + 1);
    for (size_t idx = 0; idx <= count; ++idx) {
        const double t = double(idx) / double(count);
        const double blend = t * t * (3. - 2. * t);
        result.push_back({anchor_z + t * rise + blend * span, blend});
    }
    return result;
}

std::vector<ExtrusionPath> split_nonplanar_path_by_height(
    const ExtrusionPath &path, const std::vector<double> &segment_heights)
{
    assert(path.polyline.points.size() == segment_heights.size() + 1);
    ExtrusionPaths result;
    for (size_t idx = 0; idx < segment_heights.size(); ++idx) {
        const double height = segment_heights[idx];
        assert(std::isfinite(height) && height > 0.);
        if (!result.empty() && std::abs(double(result.back().height) - height) < 1e-6) {
            result.back().polyline.points.push_back(path.polyline.points[idx + 1]);
            continue;
        }
        // Copy metadata, not the entire source polyline for every thickness
        // change. The latter makes finely sampled variable-flow paths quadratic.
        result.emplace_back(Polyline3(Points3{path.polyline.points[idx], path.polyline.points[idx + 1]}), path);
        ExtrusionPath &part = result.back();
        set_nonplanar_local_course_height(part, height);
    }
    for (ExtrusionPath &part : result)
        normalize_fixed_feature_transition(part);
    return result;
}

void consolidate_nonplanar_top_surfaces(
    const std::vector<Layer *> &layers, const sla::IndexedMesh &mesh,
    const std::vector<std::vector<uint8_t>> &selected_facets,
    const std::function<void(const NonplanarProgress &)> &progress)
{
    if (layers.size() < 2)
        return;

    const PrintConfig &print_config = layers.front()->object()->print()->config();
    SequentialGantryGeometry gantry = load_sequential_gantry_geometry(print_config);
    const double nozzle_diameter = print_config.nozzle_diameter.values.empty() ? 0.4 :
        *std::max_element(print_config.nozzle_diameter.values.begin(),
                          print_config.nozzle_diameter.values.end());
    gantry.set_nozzle_tip_diameter(nozzle_diameter);
    const ToolheadClearanceFallback fallback = make_nonplanar_clearance_cone(print_config);

    size_t region_count = std::numeric_limits<size_t>::max();
    for (const Layer *layer : layers)
        region_count = std::min(region_count, layer->regions().size());

    for (size_t region_idx = 0; region_idx < region_count; ++region_idx) {
        const auto report_progress = [&, region_idx](
            double layer_position, NonplanarProgressStage stage,
            size_t current, size_t total) {
            if (progress)
                progress({
                    (double(region_idx) +
                     layer_position / double(layers.size())) /
                        double(region_count),
                    stage, current, total});
        };
        report_progress(0., NonplanarProgressStage::SurfaceTopology,
                        region_idx + 1, region_count);
        if (region_idx >= selected_facets.size() || selected_facets[region_idx].empty())
            continue;
        // Mesh connectivity and the projected union depend only on the
        // selected surface facets, not on the nominal layer run currently
        // being considered. Build them once instead of repeating the same
        // neighbor search and large Clipper union for every candidate run.
        const std::vector<NonplanarPatch> region_patches =
            build_nonplanar_patches(mesh, selected_facets[region_idx]);
        std::vector<uint8_t> emitted_facets(selected_facets[region_idx].size(), 0);
        struct EmittedPatchProjection {
            ExPolygons projection;
            double minimum_z;
            double maximum_z;
        };
        std::vector<EmittedPatchProjection> emitted_patch_projections;
        size_t layer_idx = 1;
        while (layer_idx < layers.size()) {
            if (!region_has_nonplanar_surface(*layers[layer_idx]->regions()[region_idx])) {
                report_progress(
                    double(layer_idx + 1),
                    NonplanarProgressStage::SurfaceTopology,
                    region_idx + 1, region_count);
                ++layer_idx;
                continue;
            }

            const size_t run_begin = layer_idx;
            size_t run_end = run_begin;
            while (run_end + 1 < layers.size() &&
                   region_has_nonplanar_surface(*layers[run_end + 1]->regions()[region_idx]))
                ++run_end;

            coord_t minimum_z = std::numeric_limits<coord_t>::max();
            coord_t maximum_z = std::numeric_limits<coord_t>::lowest();
            for (size_t idx = run_begin; idx <= run_end; ++idx) {
                const LayerRegion &region = *layers[idx]->regions()[region_idx];
                const coord_t layer_z = scale_(layers[idx]->print_z);
                // Eligibility belongs to the exposed skin. A nearby curved
                // outer wall can have a larger Z range even when the deck fill
                // itself is effectively flat.
                include_nonplanar_z_range(region.fills, layer_z, minimum_z, maximum_z);
                include_nonplanar_z_range(region.perimeters, layer_z, minimum_z, maximum_z);
            }
            const coord_t minimum_span = scale_(layers[run_end]->regions()[region_idx]
                ->region().config().nonplanar_top_surface_min_height.value);
            BOOST_LOG_TRIVIAL(debug) << "Nonplanar candidate run region=" << region_idx
                                     << " layers=" << run_begin << '-' << run_end
                                     << " z_span=" << unscale_(maximum_z - minimum_z)
                                     << " required_span=" << unscale_(minimum_span);
            if (minimum_z == std::numeric_limits<coord_t>::max() ||
                maximum_z - minimum_z < minimum_span) {
                // Facet components are filtered before path generation, but a
                // small edge of such a component may leak onto an otherwise
                // flat deck. Apply the configured height threshold to the
                // actual printable run as well, then restore the original
                // planar top paths instead of emitting a cosmetic Z wobble.
                for (size_t idx = run_begin; idx <= run_end; ++idx) {
                    LayerRegion &region = *layers[idx]->regions()[region_idx];
                    restore_nonplanar_as_planar(region.perimeters);
                    restore_nonplanar_as_planar(region.fills);
                }
                layer_idx = run_end + 1;
                continue;
            }

            Polygons seed_polygons;
            for (size_t source_idx = run_begin; source_idx <= run_end; ++source_idx) {
                collect_nonplanar_top_coverage(
                    layers[source_idx]->regions()[region_idx]->fills, seed_polygons);
                collect_nonplanar_top_coverage(
                    layers[source_idx]->regions()[region_idx]->perimeters, seed_polygons);
            }
            const ExPolygons seed_coverage = union_ex(seed_polygons);

            const ExtrusionPath *top_template_ptr = first_nonplanar_path(
                layers, run_begin, run_end, region_idx, erTopSolidInfill);
            if (top_template_ptr == nullptr)
                top_template_ptr = first_path_with_role(
                    layers, run_begin, run_end, region_idx, erTopSolidInfill);
            const ExtrusionPath *outline_template_ptr = first_path_with_role(
                layers, run_begin, run_end, region_idx, erExternalPerimeter);
            const ExtrusionPath *inner_outline_template_ptr = first_path_with_role(
                layers, run_begin, run_end, region_idx, erPerimeter);
            const ExtrusionPath *ironing_template_ptr = first_nonplanar_path(
                layers, run_begin, run_end, region_idx, erIroning);
            if (ironing_template_ptr == nullptr)
                ironing_template_ptr = first_path_with_role(
                    layers, run_begin, run_end, region_idx, erIroning);
            const std::optional<ExtrusionPath> outline_template =
                outline_template_ptr == nullptr ? std::nullopt :
                std::optional<ExtrusionPath>(*outline_template_ptr);
            const std::optional<ExtrusionPath> inner_outline_template =
                inner_outline_template_ptr == nullptr ? std::nullopt :
                std::optional<ExtrusionPath>(*inner_outline_template_ptr);
            std::optional<ExtrusionPath> top_template = top_template_ptr == nullptr ?
                std::nullopt : std::optional<ExtrusionPath>(*top_template_ptr);
            const bool top_template_from_outline =
                !top_template && outline_template.has_value();
            if (top_template_from_outline) {
                // Supply flow for the complete wall/skin generator even when
                // ordinary slicing fills a narrow rim entirely with walls.
                // Its projected domain retains holes; buried courses still
                // come from real source geometry and require support below.
                top_template = *outline_template;
                top_template->set_extrusion_role(erTopSolidInfill);
                top_template->inset_idx = 0;
            }
            const std::optional<ExtrusionPath> ironing_template = ironing_template_ptr == nullptr ?
                std::nullopt : std::optional<ExtrusionPath>(*ironing_template_ptr);

            std::vector<NonplanarPatch> patches = active_nonplanar_patches(
                region_patches, mesh, seed_coverage, emitted_facets,
                unscale_(minimum_z) + mesh.ground_level(),
                unscale_(maximum_z) + mesh.ground_level(),
                layers[run_begin]->height, layers[run_end]->height);
            BOOST_LOG_TRIVIAL(debug) << "Nonplanar candidate run has " << patches.size()
                                     << " connected surface patches";
            const PrintRegionConfig &region_config =
                layers[run_end]->regions()[region_idx]->region().config();
            std::optional<Arachne::WallToolPathsParams> arachne_params;
            if (layers.front()->object()->config().wall_generator.value ==
                PerimeterGeneratorType::Arachne) {
                arachne_params = Arachne::make_paths_params(
                    int(run_end), layers.front()->object()->config(), print_config);
                arachne_params->is_top_or_bottom_layer = true;
            }
            const size_t outline_count = size_t(std::max(
                1, region_config.wall_loops.value));
            const double resolution = std::max(0.05,
                region_config.nonplanar_top_surface_resolution.value);
            const size_t shell_count = size_t(std::max(1,
                region_config.nonplanar_top_surface_layers.value));
            const double shell_layer_height = std::max(
                double(top_template ? top_template->height : layers[run_end]->height), 0.01);
            // The configured value counts the visible non-planar courses.
            // One additional buried course is generated one nominal height
            // below them.  It is the surface-following equivalent of internal
            // solid infill: without it, the low side of a broad roof may have
            // no deposited bead below the first finishing course even though
            // the mesh itself is solid there.  Unlike the former global flat
            // foundation, this support course is a parallel surface offset,
            // so every adjacent course keeps the configured physical height.
            // `nonplanar_top_surface_layers` is the number of physical surface
            // courses, including the exposed course.  The former reconstruction
            // treated it as a minimum and added one course for every nominal Z
            // plane crossed by the complete patch.  A five-layer Benchy roof
            // consequently became 28--33 generated courses and replaced deep
            // structural geometry with repeated projected rings.  A non-planar
            // skin is instead a bounded stack of parallel surface offsets; the
            // conventional slice below it remains its support.
            const size_t physical_shell_count = std::max<size_t>(1, shell_count);
            const double total_shell_depth = shell_layer_height *
                double(physical_shell_count);

            // A complete connected organic surface is always attempted first.
            // If it cannot pass dependency, support or clearance validation,
            // retry only its highest coherent cap. This is deliberately a
            // one-time fallback: recursively shaving a surface until a few
            // lines happen to pass would recreate unsafe cosmetic fragments.
            // Two finishing-stack depths retain enough conventional material
            // beneath the cap for a gradual transition, while the 1 mm floor
            // gives fine-layer profiles a visibly useful smoothing region.
            patches.reserve(std::max<size_t>(1, 2 * patches.size()));
            auto enqueue_top_cap = [&](const NonplanarPatch &source,
                                       size_t source_idx,
                                       const char *reason) {
                if (source.top_cap)
                    return;
                const double cap_height = std::max(1., 2. * total_shell_depth);
                if (source.maximum_z - source.minimum_z <= cap_height + EPSILON)
                    return;
                NonplanarPatch cap = clip_patch_to_z_range(
                    source, mesh, source.maximum_z - cap_height,
                    source.maximum_z);
                if (cap.projection.empty() ||
                    cap.maximum_z - cap.minimum_z + EPSILON <
                        unscale_(minimum_span))
                    return;
                cap.top_cap = true;
                BOOST_LOG_TRIVIAL(debug)
                    << "Retrying nonplanar patch " << source_idx
                    << " as top cap z=[" << cap.minimum_z << ','
                    << cap.maximum_z << "] after " << reason;
                patches.emplace_back(std::move(cap));
            };

            struct GeneratedShell {
                size_t owner_layer {0};
                double surface_depth {0.};
                double blend {1.};
                ExtrusionEntitiesPtr paths;
            };
            struct FeatureCourse {
                size_t owner_layer {0};
                double surface_reference_z {0.};
                double surface_depth {0.};
                double blend {1.};
            };
            struct AcceptedPatch {
                size_t patch_index {0};
                ExPolygons coverage;
                ExPolygons feature_footprint;
                std::vector<FeatureCourse> feature_courses;
                size_t surface_owner {0};
                size_t transition_anchor_layer {0};
                size_t feature_source_begin {0};
                double first_shell_blend {1.};
                bool owns_generated_boundary {false};
                bool before_current_layer {false};
                bool native_foundation {false};
            };
            std::vector<ExtrusionEntitiesPtr> generated_by_layer(layers.size());
            std::vector<AcceptedPatch> accepted_patches;
            size_t first_operation_layer = run_begin;
            size_t last_operation_layer = run_end;

            for (size_t patch_idx = 0; patch_idx < patches.size(); ++patch_idx) {
                const NonplanarPatch &patch = patches[patch_idx];
                if (!top_template)
                    continue;

                // Facet ownership alone does not prevent a second connected
                // component or queued top-cap fallback from projecting onto
                // an already accepted physical roof. Its flat foundation is
                // then emitted inside the first patch's transition stack,
                // creating unsupported full-width paths and duplicate top
                // fill. Reject only substantial overlap in both XY and Z;
                // vertically separate stacked surfaces remain independent.
                const double projected_patch_area =
                    std::abs(area(patch.projection));
                const bool projection_already_owned = std::any_of(
                    emitted_patch_projections.begin(),
                    emitted_patch_projections.end(),
                    [&patch, projected_patch_area](
                        const EmittedPatchProjection &emitted) {
                        if (patch.minimum_z > emitted.maximum_z + EPSILON ||
                            emitted.minimum_z > patch.maximum_z + EPSILON)
                            return false;
                        const double emitted_area =
                            std::abs(area(emitted.projection));
                        const double shared_area = std::abs(area(intersection_ex(
                            patch.projection, emitted.projection)));
                        return shared_area >
                            0.25 * std::min(projected_patch_area, emitted_area);
                    });
                if (projection_already_owned) {
                    BOOST_LOG_TRIVIAL(debug)
                        << "Skipping nonplanar patch " << patch_idx
                        << ": its projected surface is already owned by an accepted patch";
                    continue;
                }

                const BoundingBox patch_bounds = get_extents(patch.projection);
                BOOST_LOG_TRIVIAL(debug)
                    << "Analyzing nonplanar patch " << patch_idx
                    << " bbox=[" << unscale_(patch_bounds.min.x()) << ','
                    << unscale_(patch_bounds.min.y()) << "]-["
                    << unscale_(patch_bounds.max.x()) << ','
                    << unscale_(patch_bounds.max.y()) << "] z=["
                    << patch.minimum_z << ',' << patch.maximum_z
                    << "] area=" << std::abs(area(patch.projection)) *
                                           SCALING_FACTOR * SCALING_FACTOR;

                const double minimum_surface_z =
                    patch.minimum_z - mesh.ground_level();
                size_t transition_anchor_layer =
                    planar_transition_anchor_layer_index(
                        layers, minimum_surface_z, total_shell_depth);
                const double maximum_surface_z =
                    patch.maximum_z - mesh.ground_level();
                double minimum_course_height = 0.01;
                double maximum_course_height = 0.8 * nozzle_diameter;
                for (double height : print_config.min_layer_height.values)
                    if (height > 0.)
                        minimum_course_height = std::max(minimum_course_height, height);
                for (double height : print_config.max_layer_height.values)
                    if (height > 0.)
                        maximum_course_height = std::min(maximum_course_height, height);
                for (double diameter : print_config.nozzle_diameter.values)
                    if (diameter > 0.)
                        maximum_course_height = std::min(maximum_course_height, 0.8 * diameter);
                std::vector<NonplanarCourseProfile> course_plan;
                // Start earlier when the low side cannot hold the courses
                // required by the high side. Foundation geometry/support is
                // still validated before any original material is removed.
                for (;;) {
                    course_plan = plan_nonplanar_courses(
                        layers[transition_anchor_layer]->print_z,
                        minimum_surface_z, maximum_surface_z, physical_shell_count,
                        minimum_course_height, maximum_course_height);
                    if (!course_plan.empty() || transition_anchor_layer == 0)
                        break;
                    --transition_anchor_layer;
                }
                if (course_plan.empty()) {
                    enqueue_top_cap(patch, patch_idx, "no printable transition course spacing");
                    continue;
                }
                const double transition_anchor_z = layers[transition_anchor_layer]->print_z;
                const size_t transition_shell_count = course_plan.size() - 1;
                BOOST_LOG_TRIVIAL(debug)
                    << "Nonplanar patch " << patch_idx
                    << " transition anchor layer=" << transition_anchor_layer
                    << " z=" << transition_anchor_z
                    << " requested finishing shells=" << shell_count
                    << " transition steps=" << transition_shell_count
                    << " surface z=[" << minimum_surface_z << ','
                    << maximum_surface_z << ']';

                report_progress(
                    double(run_begin), NonplanarProgressStage::SurfaceProjection,
                    patch_idx + 1, patches.size());

                std::vector<GeneratedShell> shells;
                bool native_foundation = false;
                // Occupancy is fixed before support/clearance planning. Keep
                // the aligned wall/raster indices for projection, but never
                // let a row omitted from printing act as deposited support.
                std::map<const ExtrusionEntity *, bool> course_occupancy;
                const auto is_emitted = [&](const ExtrusionEntity *entity) {
                    const auto found = course_occupancy.find(entity);
                    return found == course_occupancy.end() || found->second;
                };
                // Course zero is a real flat foundation. The remaining
                // courses interpolate from that course to the surface.
                shells.reserve(transition_shell_count + 1);
                const size_t surface_support_begin = nonplanar_solid_course_begin(
                    course_plan, minimum_surface_z, maximum_surface_z,
                    size_t(std::max(1, region_config.top_shell_layers.value)),
                    region_config.top_shell_thickness.value);
                for (size_t shell_idx = 0;
                     shell_idx <= transition_shell_count; ++shell_idx) {
                    const double blend = course_plan[shell_idx].blend;
                    const double surface_depth = maximum_surface_z - course_plan[shell_idx].plane_z;
                    // Solid versus sparse fill does not change the height
                    // field. Switching here to parallel offsets can put the
                    // first solid course below an already deposited course.
                    const double required_z =
                        maximum_surface_z - surface_depth;
                    const size_t owner_layer = owner_layer_for_z(
                        layers,
                        shell_idx == 0 ? transition_anchor_layer :
                            std::min(transition_anchor_layer + 1,
                                     layers.size() - 1),
                        layers.size() - 1, required_z);
                    GeneratedShell shell;
                    shell.owner_layer = owner_layer;
                    shell.surface_depth = surface_depth;
                    shell.blend = blend;
                    shells.emplace_back(std::move(shell));
                }

                const size_t surface_owner = shells.back().owner_layer;
                ExPolygons feature_footprint = continuing_feature_footprint(
                    layers, region_idx, surface_owner, patch,
                    2. * std::max(double(top_template->width), resolution));

                // Dependency deferral is intended for a compact continuation
                // through the finishing skin (a chimney, boss, post, or rib),
                // not for another major body that merely belongs to the same
                // connected slice component. Excluding a body comparable to
                // the patch itself can make an otherwise colliding surface
                // look printable, only for hundreds of dependent paths to
                // fail after the ordinary shell has already been replaced.
                // Area is scale-independent here because only the ratio is
                // used. The deliberately generous bound admits real attached
                // features while keeping the acceptance decision atomic for
                // large hull/cabin/storage-box continuations.
                constexpr double maximum_deferred_feature_area_ratio = 0.35;
                const double patch_area = std::abs(area(patch.projection));
                const double feature_area = std::abs(area(feature_footprint));
                if (patch_area <= 0. ||
                    feature_area > maximum_deferred_feature_area_ratio * patch_area) {
                    if (!feature_footprint.empty())
                        BOOST_LOG_TRIVIAL(debug)
                            << "Nonplanar patch " << patch_idx
                            << " ignores non-compact continuation (area ratio "
                            << feature_area / std::max(patch_area, 1.) << ')';
                    // A footprint larger than the surface is the surrounding
                    // hull/cabin continuing normally through adjacent layers,
                    // not an attached chimney or boss.  Treating it as a
                    // dependency falsely rejected the Benchy stern (ratios of
                    // 1.0--10.5 were observed) and could enqueue a smaller,
                    // disconnected cap.  Leave ordinary continuation in the
                    // normal schedule; clearance validation below remains the
                    // authority for whether the complete patch is printable.
                    feature_footprint.clear();
                }

                // Owner layers depend only on surface depth, so determine the
                // feature interface before generating any raster. Previously
                // every shell was ray-projected, destroyed, and projected a
                // second time when a continuing feature existed. Generating
                // the final aligned stack once removes the dominant duplicate
                // work on roofs with chimneys, posts, bosses, or ribs.
                // All transition shells share the same XY paths and sampled
                // model surface. Project the finishing shell once, then blend
                // its clones from the nearest supporting planar courses into
                // the final surface height field. The bottom course therefore
                // rests on the staircase instead of floating as a parallel
                // offset; successive courses progressively remove that stair
                // step. Clearance is still validated for every clone below,
                // so this avoids duplicate ray queries without weakening the
                // collision model.
                PatchRasterCache raster_cache;
                GeneratedShell &surface_shell = shells.back();
                ExPolygons source_course_domain = patch.projection;
                if (!feature_footprint.empty()) {
                    source_course_domain = diff_ex(
                        source_course_domain,
                        offset_ex(feature_footprint, scale_(
                            std::max(0., 0.5 * double(top_template->width) -
                                             0.02))),
                        ApplySafetyOffset::No);
                }
                ExtrusionEntitiesPtr projected_source_course =
                    project_nonplanar_source_course(
                        *layers[transition_anchor_layer]->regions()[region_idx],
                        mesh, patch.projection, patch.facets, source_course_domain,
                        layers[surface_shell.owner_layer]->print_z,
                        maximum_surface_z - transition_anchor_z +
                            shell_layer_height);
                const bool has_real_source_course =
                    !projected_source_course.empty();
                constexpr size_t minimum_coherent_surface_paths = 4;
                // Reusing an ordinary solid course also requires retaining
                // its native wall graph. Otherwise the generated Arachne
                // boundary owns the replacement, subject to support and
                // clearance validation of the complete stack.
                const bool surface_owns_boundary =
                    outline_template.has_value() && !has_real_source_course;
                // Preserve the ordinary slicer's raster phase, clipping and
                // wall anchoring.  The Ahlers pipeline first creates a normal
                // coherent fill/perimeter graph and only then projects that
                // graph pointwise onto the mesh.  Preferring a newly generated
                // final-surface raster here inverted those stages: its contour
                // was subsequently cloned down to the anchor and manufactured
                // flat loops where no planar extrusion existed (the isolated
                // Benchy stern and corrupted roof regression).
                //
                // A generated raster remains a fallback only when this patch
                // truly has no ordinary solid source course.  Such a fallback
                // must still pass the complete deposited-support validation
                // below; it is never allowed to stand in for a missing source
                // on a normal top surface.
                if (surface_owns_boundary) {
                    // An exposed surface needs one topology owner for its
                    // final walls and skin.  Reusing only clipped solid-infill
                    // source paths leaves the finishing course without any
                    // perimeter, while separately contouring nominal stair
                    // walls creates the detached rings this pass is replacing.
                    // Generate the complete Arachne course once and use its
                    // exact indexed clones for the bounded stack below.
                    surface_shell.paths = make_patch_raster(
                        patch, mesh, *top_template,
                        layers[surface_shell.owner_layer]->print_z, resolution,
                        gantry, region_config.top_surface_pattern.value, 0.,
                        feature_footprint.empty() ? nullptr : &feature_footprint,
                        outline_template ? &*outline_template : nullptr,
                        &raster_cache,
                        inner_outline_template ? &*inner_outline_template : nullptr,
                        outline_count,
                        arachne_params ? &*arachne_params : nullptr,
                        true);
                } else if (has_real_source_course) {
                    surface_shell.paths = std::move(projected_source_course);
                    BOOST_LOG_TRIVIAL(debug)
                        << "Nonplanar patch " << patch_idx << " reuses "
                        << surface_shell.paths.size()
                        << " ordinary-slicer source paths";
                } else if (!top_template_from_outline) {
                    surface_shell.paths = make_patch_raster(
                        patch, mesh, *top_template,
                        layers[surface_shell.owner_layer]->print_z, resolution,
                        gantry, region_config.top_surface_pattern.value, 0.,
                        feature_footprint.empty() ? nullptr : &feature_footprint,
                        outline_template ? &*outline_template : nullptr,
                        &raster_cache,
                        inner_outline_template ? &*inner_outline_template : nullptr,
                        outline_count,
                        arachne_params ? &*arachne_params : nullptr,
                        false);
                }
                // A complete generated boundary does not consume the optional
                // solid source projection. Release those temporary entities.
                destroy_entities(projected_source_course);
                if (surface_owns_boundary && raster_cache.doubly_curved) {
                    // A dome has no single perimeter direction or matching
                    // flat lower contour. Its concentric top tracks should
                    // meet the ordinary curved walls, while cloning a newly
                    // generated final crown outline down to one planar anchor
                    // creates self-intersecting rings. Keep complete generated
                    // boundaries for developable roofs and narrow slopes, but
                    // let a genuinely two-axis crown use its existing sliced
                    // wall graph and regenerate only the adaptive skin.
                    destroy_entities(surface_shell.paths);
                    surface_shell.paths = make_patch_raster(
                        patch, mesh, *top_template,
                        layers[surface_shell.owner_layer]->print_z, resolution,
                        gantry, region_config.top_surface_pattern.value, 0.,
                        feature_footprint.empty() ? nullptr : &feature_footprint,
                        outline_template ? &*outline_template : nullptr,
                        &raster_cache,
                        inner_outline_template ? &*inner_outline_template : nullptr,
                        outline_count,
                        arachne_params ? &*arachne_params : nullptr,
                        false);
                }
                if (!surface_owns_boundary && !surface_shell.paths.empty() &&
                    surface_shell.paths.size() < minimum_coherent_surface_paths) {
                    destroy_entities(surface_shell.paths);
                    surface_shell.paths.clear();
                    // A sparse ordinary source course is useful evidence that
                    // the patch is a real top surface, but two isolated rows
                    // cannot describe a coherent dome/crown. Fall back to the
                    // complete adaptive raster instead of falling through
                    // with an empty course. The resulting full patch still
                    // has to pass the same atomic support, crossing and
                    // toolhead-clearance validation below; no sparse fragment
                    // is retained merely because it happened to be clear.
                    if (!top_template_from_outline)
                        surface_shell.paths = make_patch_raster(
                            patch, mesh, *top_template,
                            layers[surface_shell.owner_layer]->print_z,
                            resolution, gantry,
                            region_config.top_surface_pattern.value, 0.,
                            feature_footprint.empty() ? nullptr :
                                &feature_footprint,
                            outline_template ? &*outline_template : nullptr,
                            &raster_cache,
                            inner_outline_template ?
                                &*inner_outline_template : nullptr,
                            outline_count,
                            arachne_params ? &*arachne_params : nullptr,
                            surface_owns_boundary);
                }
                remove_short_nonplanar_fill_paths(surface_shell.paths);
                const bool generated_surface_owns_boundary = std::any_of(
                    surface_shell.paths.begin(), surface_shell.paths.end(),
                    [](const ExtrusionEntity *entity) {
                        const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
                        return path != nullptr && is_perimeter(path->role());
                    });
                report_progress(
                    double(run_begin), NonplanarProgressStage::BoundaryRepair,
                    patch_idx + 1, patches.size());
                size_t path_count = surface_shell.paths.size();
                // A two- or three-bead sliver is not a coherent finishing
                // skin. It cannot visibly remove stair stepping, yet replacing
                // several nominal layers around it creates a disproportionate
                // boundary and has produced wall-envelope overlaps on small
                // organic details and general torture-test geometry. Require
                // enough adjacent paths to span several nozzle widths; small
                // holes and bosses remain ordinary, dimensionally accurate
                // perimeters until a complete surrounding surface can be
                // generated atomically.
                if (path_count > 0 && path_count < minimum_coherent_surface_paths) {
                    // A narrow rim can be completely covered by fewer than
                    // four variable-width walls. Count alone is not a coverage
                    // test. Admit a complete generated graph, preserving mesh
                    // holes, but still reject sparse isolated raster fragments.
                    // All support/clearance/topology checks below still apply.
                    Polygons beads;
                    for (const ExtrusionEntity *entity : surface_shell.paths)
                        entity->polygons_covered_by_width(beads, 0.f);
                    const ExPolygons missing = diff_ex(source_course_domain,
                        offset_ex(union_ex(beads), scale_(0.05)));
                    if (area(missing) > 0.05 * area(source_course_domain)) {
                        BOOST_LOG_TRIVIAL(warning)
                            << "Rejecting nonplanar patch " << patch_idx << ": "
                            << path_count << " paths leave more than 5% of the finishing domain uncovered";
                        destroy_entities(surface_shell.paths);
                        continue;
                    }
                }
                const Flow top_flow(
                    top_template->width, top_template->height, 0.f);
                const Flow outer_flow(
                    outline_template ? outline_template->width : top_template->width,
                    outline_template ? outline_template->height : top_template->height,
                    0.f);
                const Flow inner_flow(
                    inner_outline_template ? inner_outline_template->width :
                                             outer_flow.width(),
                    inner_outline_template ? inner_outline_template->height :
                                             outer_flow.height(),
                    0.f);
                // The lower transition raster is supported inside the same
                // complete wall stack as ordinary top fill. Anchoring it only
                // half a fill width inside the raw slice allowed the cloned
                // raster to extend through one or more lower perimeter
                // centerlines. The final course still has blend one and keeps
                // the exact surface footprint generated above.
                const coord_t fill_support_inset =
                    coord_t(std::llround(0.5 * double(outer_flow.scaled_width()))) +
                    coord_t(outline_count > 0 ? outline_count - 1 : 0) *
                        inner_flow.scaled_spacing() +
                    (inner_flow.scaled_spacing() + top_flow.scaled_spacing()) / 2;
                if (path_count != 0) {
                    // The buried course is generated from the transition
                    // anchor, not by copying the exposed boundary downward.
                    // Features such as holes and chimney walls may not exist
                    // at every buried depth, so cloning the final perimeter
                    // creates unsupported rings.  The source-course morph
                    // preserves the ordinary slicer's topology while the
                    // following courses progressively adopt the surface.
                    // Each course owns independent path copies and reads only
                    // immutable slice geometry, so build them concurrently.
                    // The prior serial loop repeated the same expensive
                    // polygon projection 10--20 times on a shallow roof and
                    // dominated the multi-minute Benchy slice.
                    // Every course uses the same anchor slice and only a
                    // small set of line-width/inset combinations. Build those
                    // expensive Clipper offsets once. Previously each parallel
                    // clone rebuilt them, multiplying the roof cost by the
                    // 10--20 transition courses visible in a shallow slope.
                    std::set<coord_t> support_insets;
                    std::set<int> wall_insets;
                    for (const ExtrusionEntity *entity : surface_shell.paths)
                        if (is_perimeter(entity->role()))
                            wall_insets.insert(std::max(0, entity->inset_idx));
                    std::map<std::pair<size_t, int>, Polylines> wall_support_cache;
                    for (const ExtrusionEntity *entity : surface_shell.paths)
                        collect_shell_support_insets(
                            *entity, fill_support_inset, support_insets);
                    std::map<std::pair<size_t, coord_t>, ExPolygons>
                        support_cache;
                    // Support projection is queried only near this patch, but
                    // offsetting the complete object slice makes a local roof
                    // pay for every hull, cabin, text and multicolour contour.
                    // On the Benchy this single cache-building loop dominated
                    // minutes of CPU time. Crop before offsetting, with enough
                    // margin that the artificial window edge is farther away
                    // than any permitted transition morph.
                    const coord_t maximum_support_inset = support_insets.empty() ?
                        0 : *support_insets.rbegin();
                    BoundingBox support_bounds = patch_bounds;
                    support_bounds.offset(scale_(
                        total_shell_depth +
                        unscale_(maximum_support_inset) +
                        std::max(double(top_template->width), resolution) + 1.));
                    const ExPolygons support_window{
                        ExPolygon(support_bounds.polygon())};
                    for (size_t support_idx = transition_anchor_layer;
                         support_idx <= surface_owner; ++support_idx) {
                        for (const int inset : wall_insets)
                            wall_support_cache.emplace(std::make_pair(support_idx, inset),
                                nonplanar_wall_centerlines(
                                    layers[support_idx]->regions()[region_idx]->perimeters, inset));
                        const ExPolygons local_slices = intersection_ex(
                            layers[support_idx]->lslices, support_window);
                        for (const coord_t inset : support_insets) {
                            ExPolygons support = inset > 0 ?
                                offset_ex(local_slices, -inset) : local_slices;
                            support_cache.emplace(
                                std::make_pair(support_idx, inset),
                                std::move(support));
                        }
                    }
                    tbb::parallel_for(
                        tbb::blocked_range<size_t>(0, shells.size() - 1),
                        [&](const tbb::blocked_range<size_t> &range) {
                            for (size_t shell_idx = range.begin();
                                 shell_idx < range.end(); ++shell_idx) {
                                GeneratedShell &shell = shells[shell_idx];
                                shell.paths = clone_nonplanar_shell(
                                    surface_shell.paths, layers, patch, mesh,
                                    layers[surface_shell.owner_layer]->print_z,
                                    layers[shell.owner_layer]->print_z,
                                    transition_anchor_layer, total_shell_depth,
                                    fill_support_inset, maximum_surface_z,
                                    shell.surface_depth, shell.blend,
                                    raster_cache.doubly_curved, true,
                                    support_cache, wall_support_cache);
                            }
                        });

                    if (std::any_of(shells.begin(), shells.end(),
                            [](const GeneratedShell &shell) { return shell.paths.empty(); })) {
                        for (GeneratedShell &shell : shells)
                            destroy_entities(shell.paths);
                        enqueue_top_cap(patch, patch_idx, "incomplete transition surface projection");
                        continue;
                    }

                    // This legacy raster selection is only an occupancy
                    // proposal, not a bridge plan. Validate the actual sparse
                    // proposal below; native per-course fill/bridge generation
                    // must replace it when a dense handoff lacks support.
                    const double density = std::clamp(
                        double(region_config.sparse_infill_density.value), 0., 100.);
                    const size_t stride = density <= EPSILON ?
                        std::numeric_limits<size_t>::max() :
                        std::max<size_t>(1, size_t(std::ceil(100. / density)));
                    size_t fill_rank = 0;
                    for (size_t idx = 0; idx < surface_shell.paths.size(); ++idx) {
                        const auto *path = dynamic_cast<const ExtrusionPath *>(surface_shell.paths[idx]);
                        const bool wall = path != nullptr && is_perimeter(path->role());
                        const bool structural = path != nullptr && !wall &&
                            unscale<double>(path->length()) + EPSILON >=
                                std::max(1.0, 2.5 * double(path->width));
                        const bool sparse = structural &&
                            stride != std::numeric_limits<size_t>::max() && fill_rank % stride == 0;
                        if (structural)
                            ++fill_rank;
                        for (size_t course = 0; course < shells.size(); ++course)
                            course_occupancy.emplace(shells[course].paths[idx],
                                wall || course >= surface_support_begin || sparse);
                    }

                    // Course zero is the last fully horizontal structural
                    // course. Keep its fill in the ordinary planar schedule
                    // and preserve only scheduler ownership needed by its
                    // generated wall graph. The first path advertised as a
                    // non-planar transition is therefore shells[1], which is
                    // both supported by this complete course and actually
                    // changes Z along the surface.
                    for (ExtrusionEntity *entity : shells.front().paths)
                        make_planar_solid_foundation(
                            *entity,
                            layers[shells.front().owner_layer]->lslices);

                    // Course zero is useful only if it is a real foundation,
                    // not merely a final-surface contour flattened onto an
                    // earlier owner layer. Validate it against material that
                    // will already have been deposited on the immediately
                    // preceding planar course. The broad Benchy stern exposed
                    // this distinction: its generated ring was collision-free
                    // in empty space but the nearest matching wall was 2--8 mm
                    // below, so every later transition course floated as an
                    // isolated shell. Reject the complete patch before any
                    // conventional geometry is carved; never salvage a few
                    // individually supported fragments.
                    bool foundation_supported = transition_anchor_layer > 0;
                    if (foundation_supported) {
                        const LayerRegion &support_region =
                            *layers[transition_anchor_layer - 1]
                                 ->regions()[region_idx];
                        std::vector<const ExtrusionPath *> support_paths;
                        collect_extrusion_path_references(
                            support_region.perimeters, support_paths);
                        collect_extrusion_path_references(
                            support_region.fills, support_paths);
                        foundation_supported = std::all_of(
                            shells.front().paths.begin(),
                            shells.front().paths.end(),
                            [&support_paths, &layers, &raster_cache,
                             transition_anchor_layer](
                                const ExtrusionEntity *entity) {
                                return entity != nullptr &&
                                    nonplanar_walls_are_supported(
                                        *entity, support_paths,
                                        layers[transition_anchor_layer]->print_z,
                                        layers[transition_anchor_layer - 1]
                                            ->print_z,
                                        raster_cache.doubly_curved);
                            });
                    }
                    if (!foundation_supported && shells.size() > 1) {
                        // A native bridge course can be a valid foundation
                        // even when a copied finishing outline has no wall
                        // directly underneath it. Reuse the native course
                        // only when every first angled bead contacts it.
                        const LayerRegion &anchor =
                            *layers[transition_anchor_layer]->regions()[region_idx];
                        const double owner_z = layers[shells[1].owner_layer]->print_z;
                        const double anchor_z = layers[transition_anchor_layer]->print_z;
                        native_foundation = native_foundation_supports_course(
                            shells[1].paths, anchor, owner_z, anchor_z,
                            minimum_course_height, maximum_course_height);
                        foundation_supported = native_foundation;
                        if (native_foundation)
                            BOOST_LOG_TRIVIAL(debug) << "Nonplanar patch " << patch_idx
                                << " reuses native foundation at z=" << anchor_z;
                    }
                    if (!foundation_supported) {
                        BOOST_LOG_TRIVIAL(debug)
                            << "Rejecting nonplanar patch " << patch_idx
                            << ": its flat foundation is not continuously "
                               "bonded to the preceding planar course";
                        for (GeneratedShell &shell : shells)
                            destroy_entities(shell.paths);
                        enqueue_top_cap(
                            patch, patch_idx,
                            "unsupported flat transition foundation");
                        continue;
                    }

                    // Local course thickness is assigned when the connected
                    // paths are assembled below. XYZ length alone does not
                    // account for changing separation between the courses.

                    // XY morphing is evaluated independently for each
                    // Arachne fragment. Near a variable-width edge, two
                    // neighboring fragments may reach the same XY at
                    // different blend rates; an intermediate course can then
                    // sit more than one physical layer above the nearest bead
                    // even though the complete stack has enough courses.
                    // Constrain the emitted stack sequentially against the
                    // material actually present in the preceding aligned
                    // shell. This is a Z-only correction: XY topology,
                    // perimeter order, and the final surface target remain
                    // unchanged whenever the planned course is supported.
                    struct ShellSupportSegment {
                        Vec2d first;
                        Vec2d last;
                        double first_z;
                        double last_z;
                    };
                    constexpr double shell_support_radius = 0.4;
                    constexpr double shell_support_bucket = 1.0;
                    for (size_t shell_idx = 1;
                         shell_idx < shells.size(); ++shell_idx) {
                        const GeneratedShell &previous = shells[shell_idx - 1];
                        GeneratedShell &current = shells[shell_idx];
                        std::map<std::pair<int64_t, int64_t>,
                                 std::vector<ShellSupportSegment>> support_buckets;
                        for (const ExtrusionEntity *entity : previous.paths) {
                            if (!is_emitted(entity))
                                continue;
                            const auto *path =
                                dynamic_cast<const ExtrusionPath *>(entity);
                            if (path == nullptr)
                                continue;
                            const Points3 &points = path->polyline.points;
                            for (size_t idx = 1; idx < points.size(); ++idx) {
                                const Vec2d first =
                                    unscale(points[idx - 1].to_point())
                                        .cast<double>();
                                const Vec2d last =
                                    unscale(points[idx].to_point())
                                        .cast<double>();
                                const ShellSupportSegment segment{
                                    first, last,
                                    layers[previous.owner_layer]->print_z +
                                        unscale_(points[idx - 1].z()),
                                    layers[previous.owner_layer]->print_z +
                                        unscale_(points[idx].z())};
                                const int64_t min_x = int64_t(std::floor(
                                    (std::min(first.x(), last.x()) -
                                     shell_support_radius) /
                                    shell_support_bucket));
                                const int64_t max_x = int64_t(std::floor(
                                    (std::max(first.x(), last.x()) +
                                     shell_support_radius) /
                                    shell_support_bucket));
                                const int64_t min_y = int64_t(std::floor(
                                    (std::min(first.y(), last.y()) -
                                     shell_support_radius) /
                                    shell_support_bucket));
                                const int64_t max_y = int64_t(std::floor(
                                    (std::max(first.y(), last.y()) +
                                     shell_support_radius) /
                                    shell_support_bucket));
                                for (int64_t x = min_x; x <= max_x; ++x)
                                    for (int64_t y = min_y; y <= max_y; ++y)
                                        support_buckets[{x, y}].push_back(segment);
                            }
                        }
                        const double current_owner_z =
                            layers[current.owner_layer]->print_z;
                        for (ExtrusionEntity *entity : current.paths) {
                            if (!is_emitted(entity))
                                continue;
                            auto *path = dynamic_cast<ExtrusionPath *>(entity);
                            if (path == nullptr)
                                continue;
                            bool changed = false;
                            for (Point3 &point : path->polyline.points) {
                                const Vec2d sample =
                                    unscale(point.to_point()).cast<double>();
                                const int64_t bucket_x = int64_t(std::floor(
                                    sample.x() / shell_support_bucket));
                                const int64_t bucket_y = int64_t(std::floor(
                                    sample.y() / shell_support_bucket));
                                double nearest_distance_squared =
                                    std::numeric_limits<double>::infinity();
                                double support_z = 0.;
                                const auto found =
                                    support_buckets.find({bucket_x, bucket_y});
                                if (found != support_buckets.end())
                                    for (const ShellSupportSegment &segment :
                                         found->second) {
                                        const Vec2d delta =
                                            segment.last - segment.first;
                                        const double squared =
                                            delta.squaredNorm();
                                        const double ratio =
                                            squared <= EPSILON ? 0. :
                                            std::clamp(
                                                (sample - segment.first)
                                                        .dot(delta) /
                                                    squared,
                                                0., 1.);
                                        const double distance_squared =
                                            (sample - (segment.first +
                                                       ratio * delta))
                                                .squaredNorm();
                                        if (distance_squared >=
                                            nearest_distance_squared)
                                            continue;
                                        nearest_distance_squared =
                                            distance_squared;
                                        support_z = segment.first_z +
                                            ratio * (segment.last_z -
                                                     segment.first_z);
                                    }
                                if (nearest_distance_squared >
                                    shell_support_radius *
                                        shell_support_radius)
                                    continue;
                                const double current_z = current_owner_z +
                                    unscale_(point.z());
                                const double maximum_z =
                                    support_z + maximum_course_height;
                                if (current_z > maximum_z + EPSILON) {
                                    point.z() = scale_(maximum_z -
                                                       current_owner_z);
                                    changed = true;
                                }
                            }
                            if (changed)
                                path->polyline.fitting_result.clear();
                        }
                    }
                    // Every generated course is an indexed clone of the
                    // finishing path.  Enforce spacing on those corresponding
                    // samples in both directions: the forward pass keeps each
                    // buried course clear of its support, while the backward
                    // pass keeps the immutable mesh-following surface clear
                    // by moving only buried courses downward.  This avoids a
                    // late pairwise repair which can collapse unrelated
                    // contours onto the exposed roof.
                    constexpr double preferred_generated_course_clearance =
                        0.078;
                    for (size_t path_idx = 0;
                         path_idx < surface_shell.paths.size(); ++path_idx) {
                        std::vector<ExtrusionPath *> course_paths;
                        course_paths.reserve(shells.size());
                        bool aligned = true;
                        size_t point_count = 0;
                        for (GeneratedShell &shell : shells) {
                            auto *path = dynamic_cast<ExtrusionPath *>(
                                shell.paths[path_idx]);
                            if (path == nullptr ||
                                (!course_paths.empty() &&
                                 path->polyline.points.size() != point_count)) {
                                aligned = false;
                                break;
                            }
                            point_count = path->polyline.points.size();
                            course_paths.emplace_back(path);
                        }
                        if (!aligned || point_count == 0)
                            continue;
                        for (size_t point_idx = 0;
                             point_idx < point_count; ++point_idx) {
                            const double foundation_z =
                                layers[shells.front().owner_layer]->print_z +
                                unscale_(course_paths.front()
                                             ->polyline.points[point_idx].z());
                            const double surface_z =
                                layers[shells.back().owner_layer]->print_z +
                                unscale_(course_paths.back()
                                             ->polyline.points[point_idx].z());
                            const double generated_course_clearance = std::min(
                                preferred_generated_course_clearance,
                                std::max(0., surface_z - foundation_z) /
                                    double(shells.size() - 1));
                            for (size_t shell_idx = 1;
                                 shell_idx + 1 < shells.size(); ++shell_idx) {
                                const double previous_z =
                                    layers[shells[shell_idx - 1].owner_layer]
                                        ->print_z +
                                    unscale_(course_paths[shell_idx - 1]
                                                 ->polyline.points[point_idx]
                                                 .z());
                                Point3 &point = course_paths[shell_idx]
                                                    ->polyline.points[point_idx];
                                const double owner_z =
                                    layers[shells[shell_idx].owner_layer]
                                        ->print_z;
                                const double current_z =
                                    owner_z + unscale_(point.z());
                                if (current_z < previous_z +
                                                    generated_course_clearance)
                                    point.z() = scale_(
                                        previous_z + generated_course_clearance -
                                        owner_z);
                            }
                            for (size_t shell_idx = shells.size() - 1;
                                 shell_idx > 1; --shell_idx) {
                                const double current_z =
                                    layers[shells[shell_idx].owner_layer]
                                        ->print_z +
                                    unscale_(course_paths[shell_idx]
                                                 ->polyline.points[point_idx]
                                                 .z());
                                Point3 &previous = course_paths[shell_idx - 1]
                                                       ->polyline.points[point_idx];
                                const double previous_owner_z =
                                    layers[shells[shell_idx - 1].owner_layer]
                                        ->print_z;
                                const double previous_z = previous_owner_z +
                                    unscale_(previous.z());
                                if (previous_z > current_z -
                                                   generated_course_clearance)
                                    previous.z() = scale_(
                                        current_z - generated_course_clearance -
                                        previous_owner_z);
                            }
                        }
                        for (ExtrusionPath *path : course_paths)
                            path->polyline.fitting_result.clear();
                    }
                    // Candidate source paths may already carry transition
                    // metadata from the marking pass. The authoritative last
                    // course is the exposed surface, not another buried
                    // transition, so reset its preview/ordering role here.
                    for (ExtrusionEntity *entity : surface_shell.paths)
                        set_nonplanar_transition(*entity, false);
                }
                // Arachne boundary repair around a clipped feature can emit
                // more than tiny crossing tips. On the Benchy roof it emitted
                // one narrow, full-height boundary strand across thirteen
                // otherwise valid raster rows. Validate the complete path
                // graph at physical Z and greedily remove the strand involved
                // in the most crossings (then the lower-volume strand for a
                // tie). Apply identical removals to every aligned transition
                // shell so course topology remains consistent.
                std::vector<std::vector<BoundingBox>> shell_path_bounds(
                    shells.size(), std::vector<BoundingBox>(path_count));
                for (size_t shell_idx = 0; shell_idx < shells.size(); ++shell_idx)
                    for (size_t path_idx = 0; path_idx < path_count; ++path_idx)
                        if (const auto *path =
                                dynamic_cast<const ExtrusionPath *>(
                                    shells[shell_idx].paths[path_idx]);
                            path != nullptr && path->polyline.points.size() >= 2)
                            shell_path_bounds[shell_idx][path_idx] =
                                path->polyline.to_polyline().bounding_box();
                std::vector<uint8_t> retained_raster(path_count, uint8_t(1));
                std::vector<size_t> crossing_counts(path_count, 0);
                for (;;) {
                    std::fill(crossing_counts.begin(), crossing_counts.end(), 0);
                    size_t crossing_pairs = 0;
                    for (size_t left = 0; left < path_count; ++left) {
                        if (!retained_raster[left])
                            continue;
                        for (size_t right = left + 1; right < path_count; ++right) {
                            if (!retained_raster[right])
                                continue;
                            bool pair_crosses = false;
                            // The surface course may be collision-free while
                            // the same two aligned paths cross after their
                            // lower transition courses are inset and blended
                            // toward the planar support.  Validate every
                            // physical course which will actually be emitted.
                            for (size_t shell_idx = 0;
                                 shell_idx < shells.size(); ++shell_idx) {
                                const GeneratedShell &shell = shells[shell_idx];
                                if (!is_emitted(shell.paths[left]) || !is_emitted(shell.paths[right]))
                                    continue;
                                if (!shell_path_bounds[shell_idx][left].overlap(
                                        shell_path_bounds[shell_idx][right]))
                                    continue;
                                const auto *left_path =
                                    dynamic_cast<const ExtrusionPath *>(
                                        shell.paths[left]);
                                const auto *right_path =
                                    dynamic_cast<const ExtrusionPath *>(
                                        shell.paths[right]);
                                // make_patch_raster() owns only the finishing
                                // skin here (outlines are deliberately not
                                // emitted).  Variable-width Arachne repairs
                                // may retain a gap-fill-derived extrusion role
                                // even though the G-code processor correctly
                                // classifies the resulting non-planar path as
                                // top surface.  Role-filtering therefore hid
                                // exactly the long boundary-repair strand this
                                // pass is intended to reject.
                                if (left_path == nullptr || right_path == nullptr)
                                    continue;
                                if (paths_have_proper_same_z_crossing(
                                        *left_path,
                                        layers[shell.owner_layer]->print_z,
                                        *right_path,
                                        layers[shell.owner_layer]->print_z,
                                        0.05)) {
                                    BOOST_LOG_TRIVIAL(debug)
                                        << "Nonplanar aligned path crossing shell="
                                        << shell_idx << " left=" << left
                                        << " inset=" << left_path->inset_idx
                                        << " role=" << int(left_path->role())
                                        << " right=" << right
                                        << " inset=" << right_path->inset_idx
                                        << " role=" << int(right_path->role());
                                    pair_crosses = true;
                                    break;
                                }
                            }
                            if (!pair_crosses)
                                continue;
                            ++crossing_counts[left];
                            ++crossing_counts[right];
                            ++crossing_pairs;
                        }
                    }
                    if (crossing_pairs == 0)
                        break;
                    size_t victim = size_t(-1);
                    size_t victim_crossing_count = 0;
                    double victim_volume = std::numeric_limits<double>::max();
                    bool victim_is_wall = true;
                    for (size_t path_idx = 0; path_idx < path_count; ++path_idx) {
                        if (!retained_raster[path_idx] ||
                            crossing_counts[path_idx] == 0)
                            continue;
                        const auto *path = dynamic_cast<const ExtrusionPath *>(
                            surface_shell.paths[path_idx]);
                        if (path == nullptr)
                            continue;
                        const bool path_is_wall = is_perimeter(path->role());
                        const double volume = path->mm3_per_mm *
                            unscale<double>(path->length());
                        // The generated Arachne wall is the dimensional
                        // boundary for a fully replaced top surface. Never
                        // sacrifice it while any crossing skin or repair
                        // strand remains. Selecting the global maximum first
                        // used to delete an inner outline that crossed five
                        // independent raster strands, even though removing
                        // those five non-wall strands preserved the requested
                        // wall count and allowed the later gap repair to
                        // replace only the unsafe local skin. Within the same
                        // structural class, remove the path resolving the most
                        // crossings and use deposited volume as the final tie
                        // breaker.
                        if ((victim == size_t(-1)) ||
                            (victim_is_wall && !path_is_wall) ||
                            (victim_is_wall == path_is_wall &&
                             (crossing_counts[path_idx] >
                                  victim_crossing_count ||
                              (crossing_counts[path_idx] ==
                                   victim_crossing_count &&
                               volume < victim_volume)))) {
                            victim = path_idx;
                            victim_crossing_count = crossing_counts[path_idx];
                            victim_volume = volume;
                            victim_is_wall = path_is_wall;
                        }
                    }
                    if (victim == size_t(-1))
                        break;
                    if (const auto *path =
                            dynamic_cast<const ExtrusionPath *>(
                                surface_shell.paths[victim]))
                        BOOST_LOG_TRIVIAL(debug)
                            << "Removing aligned nonplanar path=" << victim
                            << " inset=" << path->inset_idx
                            << " role=" << int(path->role())
                            << " wall=" << is_perimeter(path->role())
                            << " crossing_count=" << victim_crossing_count;
                    retained_raster[victim] = 0;
                }
                // A single projected Arachne line can also fold across
                // itself after its XY contour is morphed toward the planar
                // foundation. Pairwise validation cannot see that. Remove
                // the aligned path from every shell when any physical course
                // self-intersects; leaving it creates a real double extrusion
                // even though the path count and endpoint graph look valid.
                for (size_t path_idx = 0; path_idx < path_count; ++path_idx) {
                    if (!retained_raster[path_idx])
                        continue;
                    for (size_t shell_idx = 0;
                         shell_idx < shells.size(); ++shell_idx) {
                        const auto *path = dynamic_cast<const ExtrusionPath *>(
                            shells[shell_idx].paths[path_idx]);
                        if (path == nullptr ||
                            !path_has_proper_same_z_self_crossing(
                                *path,
                                layers[shells[shell_idx].owner_layer]->print_z,
                                0.05))
                            continue;
                        BOOST_LOG_TRIVIAL(debug)
                            << "Removing self-crossing aligned nonplanar path="
                            << path_idx << " shell=" << shell_idx
                            << " inset=" << path->inset_idx
                            << " role=" << int(path->role());
                        retained_raster[path_idx] = 0;
                        break;
                    }
                }
                bool wall_graph_preserved = true;
                for (GeneratedShell &shell : shells)
                    if (!filter_nonplanar_course_paths(shell.paths, retained_raster)) {
                        wall_graph_preserved = false;
                        break;
                    }
                if (!wall_graph_preserved) {
                    // These are still temporary courses. No ordinary wall or
                    // fill may be carved for a candidate missing an inset.
                    for (GeneratedShell &shell : shells)
                        destroy_entities(shell.paths);
                    enqueue_top_cap(patch, patch_idx,
                        "surface projection would remove a required perimeter");
                    continue;
                }
                path_count = surface_shell.paths.size();
                std::map<unsigned int, size_t> retained_outlines_by_inset;
                for (const ExtrusionEntity *entity : surface_shell.paths)
                    if (const auto *path =
                            dynamic_cast<const ExtrusionPath *>(entity);
                        path != nullptr && is_perimeter(path->role()))
                        ++retained_outlines_by_inset[path->inset_idx];
                std::ostringstream retained_outline_summary;
                for (const auto &[inset, count] : retained_outlines_by_inset)
                    retained_outline_summary << inset << ':' << count << ' ';
                BOOST_LOG_TRIVIAL(debug)
                    << "Retained nonplanar Arachne outlines paths_by_inset="
                    << retained_outline_summary.str();
                // Narrow wall-dominated candidates now use the same complete
                // generated graph and support stack as broad skins above.
                // Do not revive the legacy independently projected-loop route:
                // it emitted isolated rings without a connected foundation.
                const bool perimeter_only_patch = false;
                for (size_t shell_idx = 0; shell_idx < shells.size(); ++shell_idx)
                    BOOST_LOG_TRIVIAL(debug)
                        << "Nonplanar patch " << patch_idx << " shell "
                        << shell_idx << " owner_layer=" << shells[shell_idx].owner_layer
                        << " raster_paths=" << shells[shell_idx].paths.size();
                if (path_count == 0) {
                    BOOST_LOG_TRIVIAL(debug) << "Rejecting nonplanar patch " << patch_idx
                                             << ": no printable raster paths";
                    for (GeneratedShell &shell : shells)
                    {
                        destroy_entities(shell.paths);
                    }
                    enqueue_top_cap(patch, patch_idx, "surface reachability rejection");
                    continue;
                }
                // Connected geometry above the skin (a chimney, boss, rib,
                // etc.) is printed as a dependent sloped course after the
                // skin. It must therefore be removed from the obstacle field
                // before validating the skin. Performing this classification
                // after clearance made the dependency itself reject the roof,
                // so the deferral code could never run.
                // Raster centerlines terminate inside the patch but the
                // deposited bead and nozzle footprint extend by half a line
                // width. Treat that narrow adjoining band as part of the
                // replacement boundary; otherwise the first and last raster
                // rows collide with the conventional outline they are meant
                // to meet. The small bound deliberately cannot hide a remote
                // tower or other independent obstacle.
                const coord_t replacement_margin =
                    scale_(0.5 * top_template->width + 0.05);
                ExPolygons clearance_exclusion = offset_ex(
                    patch.projection, replacement_margin);
                if (!feature_footprint.empty()) {
                    append(clearance_exclusion, offset_ex(
                        feature_footprint, replacement_margin));
                    clearance_exclusion = union_ex(clearance_exclusion);
                }

                // Subtract the connected replacement domain from each nominal
                // layer once. Pointer prefixes are cheap views into this stable
                // cache and avoid rebuilding identical Clipper differences for
                // every transition shell.
                const size_t last_shell_owner = std::max_element(
                    shells.begin(), shells.end(),
                    [](const GeneratedShell &left, const GeneratedShell &right) {
                        return left.owner_layer < right.owner_layer;
                    })->owner_layer;
                const std::vector<PrintedToolpathObstacle> patch_obstacles =
                    printed_model_obstacles(
                        layers, last_shell_owner, &clearance_exclusion);
                const auto obstacle_prefix = [&patch_obstacles, &layers](size_t owner_layer) {
                    std::vector<const PrintedToolpathObstacle *> prefix;
                    const double maximum_z = layers[owner_layer]->print_z + EPSILON;
                    prefix.reserve(patch_obstacles.size());
                    for (const PrintedToolpathObstacle &obstacle : patch_obstacles) {
                        if (obstacle.top_z > maximum_z)
                            break;
                        prefix.emplace_back(&obstacle);
                    }
                    return prefix;
                };
                std::vector<std::optional<std::vector<const PrintedToolpathObstacle *>>>
                    remote_obstacle_cache(layers.size());

                // A finishing raster is accepted only when the same complete
                // line is printable at every requested shell depth. This
                // guarantees that the visible skin has continuous support
                // beneath it instead of a one-layer bridge or a floating end.
                std::vector<uint8_t> clear_at_every_depth(path_count, 1);
                // Finishing surfaces are emitted after the ordinary geometry
                // on their owner layer.  A connected continuing feature is
                // removed from the obstacle field above and is emitted with
                // the surface as a deferred course, so its mere presence must
                // not force the complete patch into the early schedule.
                // Accepting a path merely because it clears before the current
                // layer creates a raised surface followed by ordinary walls
                // and fill, requiring repeated clearance lifts and defeating
                // the purpose of a final non-planar skin.  Such a path is not
                // safe in the finishing schedule and the atomic patch must
                // fall back to its original planar geometry instead.
                const bool patch_requires_early_emission = false;
                for (size_t shell_idx = 0; shell_idx < shells.size(); ++shell_idx) {
                    GeneratedShell &shell = shells[shell_idx];
                    report_progress(
                        double(run_begin) +
                            double(run_end - run_begin + 1) *
                                (double(shell_idx) /
                                 double(std::max<size_t>(1, shells.size()))),
                        NonplanarProgressStage::ShellClearance,
                        shell_idx + 1, shells.size());
                    if (!remote_obstacle_cache[shell.owner_layer])
                        remote_obstacle_cache[shell.owner_layer] =
                            obstacle_prefix(shell.owner_layer);
                    const auto &remote_obstacles =
                        *remote_obstacle_cache[shell.owner_layer];
                    tbb::parallel_for(
                        tbb::blocked_range<size_t>(0, shell.paths.size()),
                        [&](const tbb::blocked_range<size_t> &range) {
                            for (size_t path_idx = range.begin(); path_idx < range.end(); ++path_idx) {
                                if (!clear_at_every_depth[path_idx] || !is_emitted(shell.paths[path_idx]))
                                    continue;
                                auto *path = dynamic_cast<ExtrusionPath *>(shell.paths[path_idx]);
                                if (path == nullptr) {
                                    clear_at_every_depth[path_idx] = 0;
                                    continue;
                                }
                                const bool clear_after_current_layer = generated_path_is_clear(
                                    *path, layers[shell.owner_layer]->print_z,
                                    gantry, fallback, remote_obstacles);
                                if (clear_after_current_layer)
                                    continue;

                                // The connected replacement domain has already
                                // been removed from this obstacle set. A
                                // remaining current-layer collision therefore
                                // belongs to unrelated geometry and cannot be
                                // repaired by printing the finishing skin
                                // early. Reject the complete surface instead.
                                path->nonplanar_clearance_validated = false;
                                clear_at_every_depth[path_idx] = 0;
                            }
                        });
                }

                report_progress(
                    double(run_begin) +
                        double(run_end - run_begin + 1) *
                            (double(patch_idx + 1) /
                             double(std::max<size_t>(1, patches.size()))),
                    NonplanarProgressStage::ShellClearance,
                    shells.size(), shells.size());

                size_t clear_path_count = size_t(std::count(
                    clear_at_every_depth.begin(), clear_at_every_depth.end(), uint8_t(1)));
                if (generated_surface_owns_boundary &&
                    clear_path_count + 1 == path_count) {
                    const auto rejected = std::find(
                        clear_at_every_depth.begin(),
                        clear_at_every_depth.end(), uint8_t(0));
                    const size_t rejected_idx = size_t(std::distance(
                        clear_at_every_depth.begin(), rejected));
                    const auto *rejected_path = dynamic_cast<const ExtrusionPath *>(
                        shells.back().paths[rejected_idx]);
                    const BoundingBox rejected_bounds = rejected_path == nullptr ?
                        BoundingBox{} :
                        rejected_path->polyline.to_polyline().bounding_box();
                    const Vec2d rejected_size = rejected_path == nullptr ?
                        Vec2d::Constant(std::numeric_limits<double>::max()) :
                        unscale(rejected_bounds.size());
                    // Boundary repair may leave one short tapered remnant
                    // between two printable Arachne anchors. It is not a roof
                    // course or wall and cannot be emitted at any transition
                    // depth. Removing the same sub-bead fragment from every
                    // aligned course preserves topology; the later deposited-
                    // coverage pass still rejects the patch if its neighbors
                    // do not close the resulting bead envelope.
                    if (rejected_path != nullptr &&
                        rejected_path->role() == erTopSolidInfill &&
                        std::max(rejected_size.x(), rejected_size.y()) <= 0.35) {
                        BOOST_LOG_TRIVIAL(debug)
                            << "Dropping one unreachable sub-bead nonplanar "
                               "repair fragment from patch " << patch_idx
                            << " bbox=" << rejected_size.transpose();
                        for (GeneratedShell &shell : shells) {
                            delete shell.paths[rejected_idx];
                            shell.paths.erase(shell.paths.begin() + rejected_idx);
                        }
                        clear_at_every_depth.erase(
                            clear_at_every_depth.begin() + rejected_idx);
                        --path_count;
                        clear_path_count = path_count;
                    }
                }
                if (clear_path_count != path_count) {
                    std::map<std::pair<int, size_t>, size_t> rejected_by_role;
                    for (size_t path_idx = 0; path_idx < path_count; ++path_idx)
                        if (!clear_at_every_depth[path_idx])
                            if (const auto *path = dynamic_cast<const ExtrusionPath *>(
                                    shells.back().paths[path_idx]))
                                ++rejected_by_role[{int(path->role()), path->inset_idx}];
                    BOOST_LOG_TRIVIAL(debug) << "Rejecting nonplanar patch " << patch_idx
                                             << ": " << clear_path_count << " of "
                                             << path_count
                                             << " paths clear at every shell depth; "
                                                "partial roofs are not emitted";
                    for (const auto &[role_and_inset, count] : rejected_by_role)
                        BOOST_LOG_TRIVIAL(debug)
                            << "  rejected " << count << " path(s), role="
                            << ExtrusionEntity::role_to_string(
                                   ExtrusionRole(role_and_inset.first))
                            << ", inset=" << role_and_inset.second;
                    for (GeneratedShell &shell : shells)
                    {
                        destroy_entities(shell.paths);
                    }
                    enqueue_top_cap(
                        patch, patch_idx, "toolhead clearance rejection");
                    continue;
                }

                const auto course_flow = [&](size_t shell_idx, const ExtrusionPath &source) {
                    if (shell_idx == 0)
                        return ExtrusionPaths{source};
                    ExtrusionPath sampled(source);
                    sampled.polyline.points = {source.first_point3()};
                    std::vector<double> heights;
                    const auto &current = course_plan[shell_idx];
                    const auto &previous = course_plan[shell_idx - 1];
                    for (size_t idx = 1; idx < source.polyline.points.size(); ++idx) {
                        const Point3 &a = source.polyline.points[idx - 1];
                        const Point3 &b = source.polyline.points[idx];
                        const double length = unscale_((b - a).cast<double>().norm());
                        const size_t steps = std::max<size_t>(1, size_t(std::ceil(length / resolution)));
                        for (size_t step = 1; step <= steps; ++step) {
                            const Vec3d midpoint = a.cast<double>() +
                                ((double(step) - 0.5) / double(steps)) * (b - a).cast<double>();
                            const Point xy(coord_t(std::llround(midpoint.x())),
                                           coord_t(std::llround(midpoint.y())));
                            const auto mesh_z = patch_surface_z_at_xy(patch, mesh, xy);
                            if (!mesh_z)
                                throw SlicingError("Could not determine non-planar course thickness");
                            const double surface_z = *mesh_z - mesh.ground_level();
                            const double height = current.plane_z - previous.plane_z +
                                (current.blend - previous.blend) * (surface_z - maximum_surface_z);
                            if (height < minimum_course_height - 1e-4 ||
                                height > maximum_course_height + 1e-4)
                                throw SlicingError("Non-planar course thickness exceeds its planned bounds");
                            const Vec3d endpoint = a.cast<double>() +
                                (double(step) / double(steps)) * (b - a).cast<double>();
                            const Point3 point(coord_t(std::llround(endpoint.x())),
                                               coord_t(std::llround(endpoint.y())),
                                               coord_t(std::llround(endpoint.z())));
                            if (point == sampled.polyline.points.back())
                                continue;
                            sampled.polyline.points.push_back(point);
                            heights.push_back(height);
                        }
                    }
                    return split_nonplanar_path_by_height(sampled, heights);
                };

                // This is the last gate before touching the native wall
                // graph. Test the emitted occupancy with the same local
                // thickness used by the writer, not the unthinned finishing
                // raster or its provisional nominal layer height. A missing
                // support row requires bridge planning; relabeling it as a
                // bridge here would not establish physical anchors.
                bool courses_supported = true;
                std::vector<std::vector<const ExtrusionPath *>> course_supports(shells.size());
                std::vector<double> support_owner_z(shells.size(), 0.);
                for (size_t course = 1; course < shells.size() && courses_supported; ++course) {
                    auto &support = course_supports[course];
                    double &support_z = support_owner_z[course];
                    support_z = layers[shells[course - 1].owner_layer]->print_z;
                    if (native_foundation && course == 1) {
                        const LayerRegion &anchor = *layers[transition_anchor_layer]->regions()[region_idx];
                        collect_extrusion_path_references(anchor.perimeters, support);
                        collect_extrusion_path_references(anchor.fills, support);
                        support_z = layers[transition_anchor_layer]->print_z;
                    } else {
                        for (const ExtrusionEntity *entity : shells[course - 1].paths)
                            if (is_emitted(entity))
                                collect_extrusion_path_references(*entity, support);
                    }
                    for (const ExtrusionEntity *entity : shells[course].paths) {
                        if (!is_emitted(entity))
                            continue;
                        const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
                        if (path == nullptr) {
                            courses_supported = false;
                            break;
                        }
                        const ExtrusionPaths parts = course_flow(course, *path);
                        courses_supported = !parts.empty() && std::all_of(parts.begin(), parts.end(),
                            [&](const ExtrusionPath &part) {
                                return nonplanar_course_is_supported(part, support,
                                    layers[shells[course].owner_layer]->print_z, support_z);
                            });
                        if (!courses_supported) {
                            BOOST_LOG_TRIVIAL(debug) << "Rejecting nonplanar patch " << patch_idx
                                << ": emitted course " << course << " has no deposited support for "
                                << ExtrusionEntity::role_to_string(path->role())
                                << "; a native fill/bridge plan is required";
                            break;
                        }
                    }
                }
                if (!courses_supported) {
                    for (GeneratedShell &shell : shells)
                        destroy_entities(shell.paths);
                    enqueue_top_cap(patch, patch_idx, "unsupported emitted course occupancy");
                    continue;
                }

                // The wall is not regenerated from the patch outline. Keep
                // the original perimeter graph that was contoured during the
                // layer pass and validate those exact, connected path pieces
                // against the same obstacle set as the skin. This preserves
                // seam placement and planar/non-planar topology established
                // by the perimeter generator.
                std::vector<std::pair<ExtrusionPath *, size_t>> retained_walls;
                const double retained_wall_width = std::max({
                    double(top_template->width),
                    outline_template ? double(outline_template->width) : 0.,
                    inner_outline_template ?
                        double(inner_outline_template->width) : 0.});
                // Only the perimeter band of an exposed surface belongs to
                // its non-planar wall graph.  Selecting the complete patch
                // interior also picked up ordinary stair boundaries and
                // feature chords which merely crossed the same XY domain;
                // those paths were then elevated as stand-alone loops.  The
                // band is deep enough for every configured wall inset while
                // preserving the exact Arachne centerlines and their normal
                // ordering.  A path which enters the band remains one mixed
                // planar/non-planar entity, so make_nonplanar_path() supplies
                // its symmetric in-path transition instead of appending a
                // second contour.
                const coord_t wall_selection_outer =
                    scale_(0.5 * retained_wall_width + 0.05);
                const coord_t wall_selection_inner = scale_(
                    (std::max(1, region_config.wall_loops.value) + 0.75) *
                        retained_wall_width +
                    0.05);
                ExPolygons wall_selection = diff_ex(
                    offset_ex(patch.projection, wall_selection_outer),
                    offset_ex(patch.projection, -wall_selection_inner),
                    ApplySafetyOffset::No);
                if (!feature_footprint.empty()) {
                    // The dependency scheduler below owns continuing
                    // chimney/boss/post walls. Contouring those paths here as
                    // ordinary roof walls and then extracting the original
                    // planar feature course later creates two coincident wall
                    // stacks. Keep one topology owner: the roof pass handles
                    // only the surrounding surface; feature deferral warps
                    // the original perimeter graph exactly once.
                    wall_selection = diff_ex(
                        wall_selection,
                        offset_ex(feature_footprint, scale_(0.05)),
                        ApplySafetyOffset::No);
                }
                // Retained slicer-owned walls must use the same course stack
                // as the replacement skin. Otherwise every nominal stair
                // course is projected independently and the first wall is
                // already sitting on a stepped contour. The lowest requested
                // course has blend zero and is therefore moved onto the real
                // planar support Z; later courses progressively introduce the
                // original mesh-derived relative Z. A layer can own more
                // than one shallow shell course; its existing perimeter can
                // represent only one physical wall. Preserve the lowest
                // owner's minimum blend so the transition always has a flat
                // anchor. Later owners use their largest blend, while the
                // generated skin supplies any intermediate support course.
                std::map<size_t, std::pair<double, double>> wall_course_by_owner;
                for (const GeneratedShell &shell : shells) {
                    if (native_foundation && shell.owner_layer == transition_anchor_layer)
                        continue;
                    auto [position, inserted] = wall_course_by_owner.emplace(
                        shell.owner_layer,
                        std::pair<double, double>{shell.surface_depth, shell.blend});
                    if (!inserted) {
                        // One ordinary wall graph can represent only one
                        // physical course. Keep the upper/final offset for a
                        // shared owner; generated skin paths supply the buried
                        // intermediate course.
                        if (shell.surface_depth < position->second.first)
                            position->second = {shell.surface_depth, shell.blend};
                    }
                }
                // Candidate discovery covers only nominal layers whose paths
                // happen to intersect the selected surface. A wide shallow
                // surface may need transition owners before or after that
                // run. Project the ordinary slicer-owned perimeter graph on
                // every real transition owner as well. This keeps outer and
                // inner walls, seams, ordering and dimensional XY contours
                // intact; only their Z field is blended below. Generating new
                // patch outlines here would create isolated concentric rings
                // and duplicate walls.
                std::vector<ExtrusionPath *> added_transition_walls;
                const size_t first_wall_owner = wall_course_by_owner.begin()->first;
                const size_t last_wall_owner = wall_course_by_owner.rbegin()->first;
                // A rejected attempt must restore the complete original wall
                // entities, not merely zero their Z offsets. Blending changes
                // flow, role and brick metadata, and contouring may split loops.
                // Leaving those changes behind corrupts the planar fallback
                // and the next patch attempted on the same owner.
                std::vector<std::pair<ExtrusionEntityCollection *, ExtrusionEntityCollection>> wall_backups;
                if (!generated_surface_owns_boundary)
                    for (const auto &[owner, course] : wall_course_by_owner) {
                        auto &perimeters = layers[owner]->regions()[region_idx]->perimeters;
                        wall_backups.emplace_back(&perimeters, perimeters);
                    }
                ScopeGuard restore_wall_graph([&wall_backups]() {
                    for (auto &[perimeters, original] : wall_backups)
                        *perimeters = std::move(original);
                });
                const double maximum_drape_height = std::max(
                    shell_layer_height,
                    maximum_surface_z - transition_anchor_z +
                        shell_layer_height);
                if (!generated_surface_owns_boundary)
                    for (const auto &[source_idx, course] :
                         wall_course_by_owner) {
                        // Course zero is the ordinary slicer's last complete
                        // horizontal perimeter.  Keep that original entity
                        // planar and let the first blended course connect to
                        // it on the following owner. Marking a generated or
                        // projected fixed-Z copy as "Non-planar transition"
                        // produced an isolated closed loop before any Z
                        // transition occurred and could clip one of the
                        // configured chimney insets from the real wall graph.
                        if (course.second <= EPSILON)
                            continue;
                        LayerRegion &source =
                            *layers[source_idx]->regions()[region_idx];
                        contour_existing_perimeter_graph(
                            source, mesh, patch.facets, source.perimeters,
                            wall_selection, maximum_drape_height,
                            added_transition_walls);
                    }
                bool retained_wall_projection_valid = true;
                const auto collect_retained_walls = [&](auto &&self,
                                                        ExtrusionEntity &entity,
                                                        size_t source_idx) -> void {
                    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
                        if (!path->nonplanar_surface ||
                            !is_perimeter(path->role()) ||
                            path->nonplanar_clearance_validated)
                            return;
                        const auto [minimum_z, maximum_z] = std::minmax_element(
                            path->polyline.points.begin(),
                            path->polyline.points.end(),
                            [](const Point3 &left, const Point3 &right) {
                                return left.z() < right.z();
                            });
                        // A wall which rounds to one fixed emitted Z is not a
                        // surface-following transition.  Such candidates were
                        // skipped by belongs_to_patch below (all relative Z
                        // values were zero), retained their provisional
                        // nonplanar metadata, and leaked out as detached
                        // closed loops several millimetres above their last
                        // real support. Preserve the slicer's authoritative
                        // owner course and classification instead.
                        if (minimum_z == path->polyline.points.end() ||
                            maximum_z == path->polyline.points.end() ||
                            maximum_z->z() - minimum_z->z() < scale_(0.001)) {
                            restore_candidate_path_as_planar(*path);
                            return;
                        }
                        const bool belongs_to_patch = std::any_of(
                            path->polyline.points.begin(),
                            path->polyline.points.end(),
                            [&wall_selection](const Point3 &point) {
                                return point.z() != 0 && Geometry::contains(
                                    wall_selection, point.to_point());
                            });
                        const auto course = wall_course_by_owner.find(source_idx);
                        if (belongs_to_patch && course != wall_course_by_owner.end()) {
                            if (!blend_wall_z_from_planar_support(
                                    *path, layers, patch, mesh,
                                    layers[source_idx]->print_z,
                                    transition_anchor_layer,
                                    layers[source_idx]->height,
                                    std::clamp(
                                        region_config
                                            .nonplanar_top_surface_min_height.value,
                                        0.01, double(layers[source_idx]->height)),
                                    1. / double(transition_shell_count),
                                    maximum_surface_z,
                                    course->second.first,
                                    course->second.second)) {
                                BOOST_LOG_TRIVIAL(debug)
                                    << "Nonplanar patch " << patch_idx
                                    << " wall projection failed: owner=" << source_idx
                                    << " inset=" << path->inset_idx
                                    << " role=" << int(path->role());
                                retained_wall_projection_valid = false;
                                return;
                            }
                            // A roof eave can begin outside every preceding
                            // wall even though the interior finishing raster
                            // has a dense transition foundation. It is a real
                            // cantilever/bridge, not a supported non-planar
                            // wall. Preserve the connected slicer-owned path,
                            // but retain the conventional overhang policy so
                            // it receives the appropriate speed, cooling and
                            // preview classification instead of masquerading
                            // as a normally supported perimeter.
                            if (source_idx > 0) {
                                const ExPolygons lower_support = offset_ex(
                                    layers[source_idx - 1]->lslices,
                                    scale_(0.5 * double(path->width) + 0.02));
                                if (std::any_of(
                                        path->polyline.points.begin(),
                                        path->polyline.points.end(),
                                        [&lower_support](const Point3 &point) {
                                            return !Geometry::contains(
                                                lower_support,
                                                point.to_point());
                                        }))
                                    path->set_extrusion_role(
                                        erOverhangPerimeter);
                            }
                            path->nonplanar_transition =
                                course->second.first > EPSILON;
                            retained_walls.emplace_back(path, source_idx);
                        }
                    } else if (auto *multipath =
                                   dynamic_cast<ExtrusionMultiPath *>(&entity)) {
                        for (ExtrusionPath &path : multipath->paths)
                            self(self, path, source_idx);
                    } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
                        for (ExtrusionPath &path : loop->paths)
                            self(self, path, source_idx);
                    } else if (auto *collection =
                                   dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
                        for (ExtrusionEntity *child : collection->entities)
                            self(self, *child, source_idx);
                    }
                };
                if (!generated_surface_owns_boundary)
                    for (size_t source_idx = first_wall_owner;
                         source_idx <= last_wall_owner; ++source_idx)
                        if (wall_course_by_owner.count(source_idx) != 0)
                            collect_retained_walls(
                                collect_retained_walls,
                                layers[source_idx]->regions()[region_idx]->perimeters,
                                source_idx);
                // The ordinary perimeter graph intentionally keeps long
                // straight segments. Once adjacent courses acquire different
                // Z profiles, however, checking only their original vertices
                // is insufficient: a ramp can pass through a nearly flat
                // supporting segment between its endpoints. Sample the exact
                // slicer-owned centerline (without changing XY geometry) at
                // the configured non-planar resolution before enforcing
                // course separation. This is the same error-control principle
                // used by mesh-projected raster paths and keeps the following
                // support pass valid over segment interiors as well as nodes.
                const double wall_sample_length = std::max(0.1, resolution);
                for (const auto &[path, unused_owner] : retained_walls) {
                    (void) unused_owner;
                    if (path->polyline.points.size() < 2)
                        continue;
                    Points3 sampled;
                    sampled.reserve(path->polyline.points.size());
                    sampled.emplace_back(path->polyline.points.front());
                    for (size_t point_idx = 1;
                         point_idx < path->polyline.points.size(); ++point_idx) {
                        const Point3 &start =
                            path->polyline.points[point_idx - 1];
                        const Point3 &end = path->polyline.points[point_idx];
                        const double xy_length = unscale_(
                            (end.to_point() - start.to_point())
                                .cast<double>().norm());
                        const size_t subdivisions = start.z() == end.z() ? 1 :
                            std::max<size_t>(
                                1, size_t(std::ceil(
                                       xy_length / wall_sample_length)));
                        for (size_t sample_idx = 1;
                             sample_idx <= subdivisions; ++sample_idx) {
                            const double ratio = double(sample_idx) /
                                double(subdivisions);
                            sampled.emplace_back(Point3(
                                coord_t(std::llround(
                                    double(start.x()) +
                                    ratio * double(end.x() - start.x()))),
                                coord_t(std::llround(
                                    double(start.y()) +
                                    ratio * double(end.y() - start.y()))),
                                coord_t(std::llround(
                                    double(start.z()) +
                                    ratio * double(end.z() - start.z())))));
                        }
                    }
                    path->polyline.points = std::move(sampled);
                    path->polyline.remove_duplicate_points();
                    path->polyline.fitting_result.clear();
                }

                // The provisional planar<->surface blend is constructed on
                // each ordinary perimeter independently.  Seams and slightly
                // different layer contours can therefore give corresponding
                // points different blend influence, allowing a later course
                // to dip below an earlier one near the handoff.  Preserve the
                // original connected XY graph, but enforce monotonic physical
                // Z wherever matching wall courses overlap.  The mesh surface
                // remains a hard upper bound, so this cannot grow a wall past
                // the requested model geometry.
                // Keep more than the 0.05 mm rounded-centerline collision
                // tolerance between successive projected wall courses.  The
                // correction is evaluated at resampled vertices, while two
                // slightly different Arachne contours may intersect between
                // those vertices; using the configured minimum course height
                // leaves enough margin for that interpolation and G-code
                // coordinate rounding without moving the exposed surface.
                const double minimum_wall_course_separation = std::max(
                    0.0601,
                    std::clamp(
                        region_config.nonplanar_top_surface_min_height.value,
                        0.01, shell_layer_height));
                for (size_t high_idx = 0;
                     high_idx < retained_walls.size(); ++high_idx) {
                    auto &[high_path, high_owner_idx] = retained_walls[high_idx];
                    const double high_owner_z = layers[high_owner_idx]->print_z;
                    const double lateral_tolerance =
                        0.75 * std::max(double(high_path->width), retained_wall_width);
                    const double lateral_tolerance_squared =
                        lateral_tolerance * lateral_tolerance;
                    size_t previous_owner_idx = size_t(-1);
                    for (const auto &[candidate, candidate_owner_idx] :
                         retained_walls)
                        if (candidate_owner_idx < high_owner_idx &&
                            candidate->inset_idx == high_path->inset_idx &&
                            (previous_owner_idx == size_t(-1) ||
                             candidate_owner_idx > previous_owner_idx))
                            previous_owner_idx = candidate_owner_idx;
                    if (previous_owner_idx == size_t(-1))
                        continue;
                    for (Point3 &high_point : high_path->polyline.points) {
                        if (high_point.z() == 0)
                            continue;
                        const Vec2d high_xy = unscale(high_point.to_point());
                        double supporting_z = std::numeric_limits<double>::lowest();
                        for (size_t low_idx = 0; low_idx < retained_walls.size(); ++low_idx) {
                            const auto &[low_path, low_owner_idx] =
                                retained_walls[low_idx];
                            if (low_owner_idx != previous_owner_idx ||
                                low_path->inset_idx != high_path->inset_idx)
                                continue;
                            const Points3 &low_points = low_path->polyline.points;
                            for (size_t segment_idx = 1;
                                 segment_idx < low_points.size(); ++segment_idx) {
                                const Vec2d first =
                                    unscale(low_points[segment_idx - 1].to_point());
                                const Vec2d last =
                                    unscale(low_points[segment_idx].to_point());
                                const Vec2d delta = last - first;
                                const double length_squared = delta.squaredNorm();
                                const double ratio = length_squared <= EPSILON ? 0. :
                                    std::clamp((high_xy - first).dot(delta) /
                                                   length_squared,
                                               0., 1.);
                                const Vec2d nearest = first + ratio * delta;
                                if ((nearest - high_xy).squaredNorm() >
                                    lateral_tolerance_squared)
                                    continue;
                                const double low_z =
                                    layers[low_owner_idx]->print_z +
                                    unscale_(low_points[segment_idx - 1].z()) +
                                    ratio * unscale_(
                                        low_points[segment_idx].z() -
                                        low_points[segment_idx - 1].z());
                                supporting_z = std::max(supporting_z, low_z);
                            }
                        }
                        if (!std::isfinite(supporting_z))
                            continue;
                        const std::optional<double> mesh_z =
                            patch_surface_z_at_xy(patch, mesh, high_point.to_point());
                        if (!mesh_z)
                            continue;
                        const double surface_z = *mesh_z - mesh.ground_level();
                        const double current_z = high_owner_z +
                            unscale_(high_point.z());
                        const double corrected_z = std::min(
                            surface_z,
                            std::max(current_z,
                                     supporting_z +
                                         minimum_wall_course_separation));
                        high_point.z() = scale_(corrected_z - high_owner_z);
                    }
                    high_path->polyline.remove_duplicate_points();
                    high_path->polyline.fitting_result.clear();
                }

                // The forward pass above prevents a later wall course from
                // dipping into an earlier one, but that is only half of the
                // support constraint.  At the low edge of a shallow roof the
                // finishing course is fixed by the mesh while its preceding
                // blended course may still be more than one printable layer
                // below it.  Do not pull the finishing surface away from the
                // model.  Propagate its required support height backwards
                // through matching slicer-owned wall centerlines instead.
                // Processing owners from high to low makes the constraint
                // reach the flat foundation in one pass.
                // Leave margin for simplification and rounded G-code Z
                // coordinates; the physical acceptance limit is one nominal
                // course plus 0.06 mm, not an exact floating-point boundary.
                constexpr double maximum_wall_course_allowance = 0.02;
                std::vector<size_t> wall_order(retained_walls.size());
                std::iota(wall_order.begin(), wall_order.end(), size_t(0));
                std::stable_sort(
                    wall_order.begin(), wall_order.end(),
                    [&retained_walls](size_t lhs, size_t rhs) {
                        return retained_walls[lhs].second >
                               retained_walls[rhs].second;
                    });
                for (size_t high_wall_idx : wall_order) {
                    auto &[high_path, high_owner_idx] =
                        retained_walls[high_wall_idx];
                    size_t previous_owner_idx = size_t(-1);
                    for (const auto &[candidate, candidate_owner_idx] :
                         retained_walls)
                        if (candidate_owner_idx < high_owner_idx &&
                            candidate->inset_idx == high_path->inset_idx &&
                            (previous_owner_idx == size_t(-1) ||
                             candidate_owner_idx > previous_owner_idx))
                            previous_owner_idx = candidate_owner_idx;
                    if (previous_owner_idx == size_t(-1))
                        continue;

                    const double high_owner_z =
                        layers[high_owner_idx]->print_z;
                    const double maximum_separation =
                        std::max(double(high_path->height), shell_layer_height) +
                        maximum_wall_course_allowance;
                    const double lateral_tolerance = 0.75 * std::max(
                        double(high_path->width), retained_wall_width);
                    const double lateral_tolerance_squared =
                        lateral_tolerance * lateral_tolerance;
                    for (const Point3 &high_point : high_path->polyline.points) {
                        const Vec2d high_xy = unscale(high_point.to_point());
                        const double high_z = high_owner_z +
                            unscale_(high_point.z());
                        ExtrusionPath *nearest_path = nullptr;
                        size_t nearest_segment = 0;
                        double nearest_distance_squared =
                            std::numeric_limits<double>::max();
                        for (auto &[low_path, low_owner_idx] : retained_walls) {
                            if (low_owner_idx != previous_owner_idx ||
                                low_path->inset_idx != high_path->inset_idx)
                                continue;
                            const Points3 &low_points = low_path->polyline.points;
                            for (size_t segment_idx = 1;
                                 segment_idx < low_points.size(); ++segment_idx) {
                                const Vec2d first = unscale(
                                    low_points[segment_idx - 1].to_point());
                                const Vec2d last = unscale(
                                    low_points[segment_idx].to_point());
                                const Vec2d delta = last - first;
                                const double length_squared =
                                    delta.squaredNorm();
                                const double ratio = length_squared <= EPSILON ?
                                    0. : std::clamp(
                                        (high_xy - first).dot(delta) /
                                            length_squared,
                                        0., 1.);
                                const double distance_squared =
                                    (first + ratio * delta - high_xy)
                                        .squaredNorm();
                                if (distance_squared < nearest_distance_squared) {
                                    nearest_distance_squared = distance_squared;
                                    nearest_path = low_path;
                                    nearest_segment = segment_idx;
                                }
                            }
                        }
                        if (nearest_path == nullptr ||
                            nearest_distance_squared > lateral_tolerance_squared)
                            continue;

                        const double low_owner_z =
                            layers[previous_owner_idx]->print_z;
                        const double required_support_z =
                            high_z - maximum_separation;
                        for (size_t endpoint :
                             {nearest_segment - 1, nearest_segment}) {
                            Point3 &low_point =
                                nearest_path->polyline.points[endpoint];
                            const double current_z = low_owner_z +
                                unscale_(low_point.z());
                            if (current_z >= required_support_z - EPSILON)
                                continue;
                            const std::optional<double> mesh_z =
                                patch_surface_z_near_xy(
                                    patch, mesh, low_point.to_point(),
                                    0.5 * double(nearest_path->width) + 0.05);
                            if (!mesh_z)
                                continue;
                            const double surface_z =
                                *mesh_z - mesh.ground_level();
                            const double corrected_z = std::min(
                                {surface_z,
                                 high_z - minimum_wall_course_separation,
                                 required_support_z});
                            if (corrected_z > current_z)
                                low_point.z() =
                                    scale_(corrected_z - low_owner_z);
                        }
                    }
                }
                for (const auto &[path, owner_idx] : retained_walls) {
                    // Back-propagating support can raise only a pair of
                    // vertices on a long source segment. Spread that change
                    // along the existing centerline so it cannot create a
                    // sharp Z kink. Raising the lower neighbor preserves the
                    // support constraint; the selected mesh is the hard cap.
                    constexpr double maximum_transition_gradient = 0.8;
                    const double owner_z = layers[owner_idx]->print_z;
                    for (size_t pass = 0;
                         pass < path->polyline.points.size(); ++pass) {
                        bool changed = false;
                        for (size_t segment_idx = 1;
                             segment_idx < path->polyline.points.size();
                             ++segment_idx) {
                            Point3 &first =
                                path->polyline.points[segment_idx - 1];
                            Point3 &last = path->polyline.points[segment_idx];
                            const double xy = unscale_(
                                (last.to_point() - first.to_point())
                                    .cast<double>().norm());
                            if (xy <= EPSILON)
                                continue;
                            const double first_z =
                                owner_z + unscale_(first.z());
                            const double last_z = owner_z + unscale_(last.z());
                            Point3 *lower = first_z < last_z ? &first : &last;
                            const double high_z = std::max(first_z, last_z);
                            const double required_z = high_z -
                                maximum_transition_gradient * xy;
                            const double low_z = std::min(first_z, last_z);
                            if (low_z >= required_z - EPSILON)
                                continue;
                            const std::optional<double> mesh_z =
                                patch_surface_z_near_xy(
                                    patch, mesh, lower->to_point(),
                                    0.5 * double(path->width) + 0.05);
                            if (!mesh_z)
                                continue;
                            const double corrected_z = std::min(
                                required_z,
                                *mesh_z - mesh.ground_level());
                            if (corrected_z > low_z + EPSILON) {
                                lower->z() = scale_(corrected_z - owner_z);
                                changed = true;
                            }
                        }
                        if (!changed)
                            break;
                    }
                    path->polyline.remove_duplicate_points();
                    path->polyline.fitting_result.clear();
                }

                bool retained_walls_clear = retained_wall_projection_valid;
                for (const auto &[path, source_idx] : retained_walls) {
                    if (!remote_obstacle_cache[source_idx])
                        remote_obstacle_cache[source_idx] = obstacle_prefix(source_idx);
                    std::string clearance_failure;
                    if (generated_path_is_clear(
                            *path, layers[source_idx]->print_z,
                            gantry, fallback, *remote_obstacle_cache[source_idx],
                            &clearance_failure))
                        continue;
                    BOOST_LOG_TRIVIAL(debug)
                        << "Nonplanar patch " << patch_idx
                        << " wall clearance failed: owner=" << source_idx
                        << " inset=" << path->inset_idx
                        << " role=" << int(path->role())
                        << ' ' << clearance_failure;
                    retained_walls_clear = false;
                    break;
                }
                if (perimeter_only_patch && retained_walls_clear) {
                    // A compact crown may span many nominal owners.  Testing
                    // every projected ring against the single flat anchor
                    // footprint made a much later ring look supported merely
                    // because its XY happened to overlap the old hull.  The
                    // resulting path could start several millimetres above
                    // the most recently deposited material (the detached
                    // Benchy stern loops were the real-world example).
                    //
                    // Validate each retained slicer course against the course
                    // which will actually precede its owner.  This is still
                    // the normal perimeter overhang policy -- short runs may
                    // cantilever -- but an entire closed post-process ring can
                    // no longer borrow support from an unrelated old layer.
                    std::map<size_t, std::vector<const ExtrusionPath *>>
                        support_by_owner;
                    for (const auto &[path, owner_idx] : retained_walls) {
                        if (owner_idx == 0) {
                            retained_walls_clear = false;
                            break;
                        }
                        auto found = support_by_owner.find(owner_idx);
                        if (found == support_by_owner.end()) {
                            const LayerRegion &support_region =
                                *layers[owner_idx - 1]->regions()[region_idx];
                            std::vector<const ExtrusionPath *> support_paths;
                            collect_extrusion_path_references(
                                support_region.perimeters, support_paths);
                            collect_extrusion_path_references(
                                support_region.fills, support_paths);
                            found = support_by_owner.emplace(
                                owner_idx, std::move(support_paths)).first;
                        }
                        // A compact feature is not exempt from the physical
                        // support rule.  Coverage polygons discard the
                        // relative Z of a previously contoured support path,
                        // so a ring at the same XY could be accepted even
                        // when the closest deposited bead was millimetres
                        // below it.  This was the source of the detached
                        // Benchy flag-hole/chimney rings. Validate against
                        // the real 3D centerlines just like a surface patch.
                        if (!nonplanar_walls_are_supported(
                                *path, found->second,
                                layers[owner_idx]->print_z,
                                layers[owner_idx - 1]->print_z,
                                raster_cache.doubly_curved)) {
                            retained_walls_clear = false;
                            break;
                        }
                    }
                    for (size_t path_idx = 0;
                         retained_walls_clear && path_idx < retained_walls.size();
                         ++path_idx) {
                        const auto &[path, owner_idx] = retained_walls[path_idx];
                        const Points3 &points = path->polyline.points;
                        for (size_t point_idx = 1; point_idx < points.size(); ++point_idx) {
                            const double xy = unscale_(
                                (points[point_idx].to_point() -
                                 points[point_idx - 1].to_point())
                                    .cast<double>().norm());
                            const double dz = std::abs(unscale_(
                                points[point_idx].z() -
                                points[point_idx - 1].z()));
                            if (xy > EPSILON && dz / xy <= 1.)
                                continue;
                            retained_walls_clear = false;
                            break;
                        }
                        for (size_t other_idx = path_idx + 1;
                             retained_walls_clear &&
                             other_idx < retained_walls.size(); ++other_idx) {
                            const auto &[other, other_owner_idx] =
                                retained_walls[other_idx];
                            if (paths_have_proper_same_z_crossing(
                                    *path, layers[owner_idx]->print_z,
                                    *other, layers[other_owner_idx]->print_z,
                                    0.05))
                                retained_walls_clear = false;
                        }
                    }
                    if (!retained_walls_clear)
                        BOOST_LOG_TRIVIAL(debug)
                            << "Rejecting compact perimeter-only nonplanar patch "
                            << patch_idx
                            << ": its wall courses cross or exceed the safe Z gradient";
                }
                // The first contoured wall course is printed at the flat
                // transition anchor.  Clearance alone cannot prove that it
                // has material underneath: on a narrowing bow it may clear
                // every obstacle while spanning millimetres beyond the last
                // deposited wall.  Validate the actual preceding extrusion
                // footprint before replacing any conventional geometry.  A
                // short unsupported run remains governed by normal overhang
                // behavior; a long run rejects the patch atomically.
                if (!perimeter_only_patch && retained_walls_clear &&
                    !raster_cache.doubly_curved) {
                    std::map<size_t, std::vector<const ExtrusionPath *>>
                        support_by_owner;
                    for (const auto &[path, owner_idx] : retained_walls) {
                        if (owner_idx == 0) {
                            retained_walls_clear = false;
                            break;
                        }
                        auto found = support_by_owner.find(owner_idx);
                        if (found == support_by_owner.end()) {
                            const LayerRegion &support_region =
                                *layers[owner_idx - 1]->regions()[region_idx];
                            std::vector<const ExtrusionPath *> support_paths;
                            collect_extrusion_path_references(
                                support_region.perimeters, support_paths);
                            collect_extrusion_path_references(
                                support_region.fills, support_paths);
                            found = support_by_owner.emplace(
                                owner_idx, std::move(support_paths)).first;
                        }
                        if (!nonplanar_walls_are_supported(
                                *path, found->second,
                                layers[owner_idx]->print_z,
                                layers[owner_idx - 1]->print_z,
                                raster_cache.doubly_curved)) {
                            retained_walls_clear = false;
                            break;
                        }
                    }
                    if (!retained_walls_clear)
                        BOOST_LOG_TRIVIAL(debug)
                            << "Rejecting nonplanar patch " << patch_idx
                            << ": a connected wall course has no continuous "
                               "deposited support on its preceding owner";
                }
                if (!retained_walls_clear) {
                    BOOST_LOG_TRIVIAL(debug)
                        << "Rejecting nonplanar patch " << patch_idx
                        << ": its connected perimeter graph does not clear the toolhead";
                    for (GeneratedShell &shell : shells)
                    {
                        destroy_entities(shell.paths);
                    }
                    for (const auto &[path, source_idx] : retained_walls) {
                        (void) source_idx;
                        path->nonplanar_clearance_validated = false;
                        restore_candidate_path_as_planar(*path);
                    }
                    for (ExtrusionPath *path : added_transition_walls) {
                        path->nonplanar_clearance_validated = false;
                        restore_candidate_path_as_planar(*path);
                    }
                    enqueue_top_cap(
                        patch, patch_idx, "connected perimeter clearance rejection");
                    continue;
                }

                // A nominal stair perimeter may have been selected and
                // projected even though the finishing raster replaces the
                // interior it encloses. Keeping that closed wall produces an
                // inner rectangular box through the roof: multiple raster
                // rows cross its centerline at the same physical Z. Preserve
                // only retained wall paths which meet the skin at endpoints
                // or remain outside it. Conflicting paths stay provisional so
                // the commit pass clips them from the replacement volume.
                std::vector<uint8_t> retained_wall_crosses_raster(
                    retained_walls.size(), 0);
                bool retained_wall_topology_conflict = false;
                std::vector<BoundingBox> retained_wall_bounds;
                std::vector<std::pair<double, double>> retained_wall_z_bounds;
                retained_wall_bounds.reserve(retained_walls.size());
                retained_wall_z_bounds.reserve(retained_walls.size());
                for (const auto &[path, owner] : retained_walls) {
                    retained_wall_bounds.emplace_back(
                        path->polyline.to_polyline().bounding_box());
                    double minimum_z = std::numeric_limits<double>::max();
                    double maximum_z = std::numeric_limits<double>::lowest();
                    for (const Point3 &point : path->polyline.points) {
                        const double z = layers[owner]->print_z +
                            unscale_(point.z());
                        minimum_z = std::min(minimum_z, z);
                        maximum_z = std::max(maximum_z, z);
                    }
                    retained_wall_z_bounds.emplace_back(minimum_z, maximum_z);
                }
                // When two successive nominal wall outlines cross after
                // projection, keep the earlier course as physical support and
                // let the aligned generated shell own only the conflicting
                // later course. Removing both sides of the pair leaves the
                // generated wall floating; retaining both deposits duplicate
                // material at the intersection. This one-sided replacement
                // is the variable-layer analogue of terminating a local
                // course as the model contour narrows.
                for (size_t wall_idx = 0;
                     wall_idx < retained_walls.size(); ++wall_idx)
                    for (size_t other_idx = wall_idx + 1;
                         other_idx < retained_walls.size(); ++other_idx) {
                        const auto &[wall, source_idx] =
                            retained_walls[wall_idx];
                        const auto &[other, other_source_idx] =
                            retained_walls[other_idx];
                        const size_t later_idx = source_idx > other_source_idx ?
                            wall_idx : other_idx;
                        if (retained_wall_crosses_raster[later_idx])
                            continue;
                        if (std::max(source_idx, other_source_idx) -
                                std::min(source_idx, other_source_idx) >
                                    transition_shell_count ||
                            !retained_wall_bounds[wall_idx].overlap(
                                retained_wall_bounds[other_idx]) ||
                            retained_wall_z_bounds[wall_idx].first >
                                retained_wall_z_bounds[other_idx].second + 0.05 ||
                            retained_wall_z_bounds[other_idx].first >
                                retained_wall_z_bounds[wall_idx].second + 0.05)
                            continue;
                        if (!paths_have_proper_same_z_crossing(
                                *wall, layers[source_idx]->print_z,
                                *other, layers[other_source_idx]->print_z,
                                0.05, 0.02))
                            continue;
                        retained_wall_crosses_raster[later_idx] = 1;
                        retained_wall_topology_conflict = true;
                    }
                for (size_t wall_idx = 0; wall_idx < retained_walls.size(); ++wall_idx) {
                    const auto &[wall, source_idx] = retained_walls[wall_idx];
                    for (const GeneratedShell &shell : shells) {
                        if (retained_wall_crosses_raster[wall_idx])
                            break;
                        for (const ExtrusionEntity *entity : shell.paths) {
                            const auto *raster = dynamic_cast<const ExtrusionPath *>(entity);
                            if (raster == nullptr || raster->role() != erTopSolidInfill)
                                continue;
                            if (paths_have_proper_same_z_crossing(
                                    *wall, layers[source_idx]->print_z,
                                    *raster, layers[shell.owner_layer]->print_z,
                                    0.05, 0.02)) {
                                retained_wall_crosses_raster[wall_idx] = 1;
                                retained_wall_topology_conflict = true;
                                break;
                            }
                        }
                    }
                }
                // The retained paths are the slicer's complete dimensional
                // wall graph.  A crossing is not repairable by dropping one
                // wall or relabeling it as planar: either choice commits a
                // hybrid patch with a missing perimeter or an unchanged
                // collision.  Keep patch replacement atomic and restore the
                // conventional slice when its entire wall graph cannot be
                // connected to the finishing skin without an overlap.
                if (retained_wall_topology_conflict) {
                    BOOST_LOG_TRIVIAL(debug)
                        << "Rejecting nonplanar patch " << patch_idx
                        << ": projected wall topology intersects another "
                           "wall course or finishing raster";
                    for (GeneratedShell &shell : shells)
                        destroy_entities(shell.paths);
                    for (const auto &[path, source_idx] : retained_walls) {
                        (void) source_idx;
                        path->nonplanar_clearance_validated = false;
                        restore_candidate_path_as_planar(*path);
                    }
                    for (ExtrusionPath *path : added_transition_walls) {
                        path->nonplanar_clearance_validated = false;
                        restore_candidate_path_as_planar(*path);
                    }
                    enqueue_top_cap(
                        patch, patch_idx, "connected wall topology rejection");
                    continue;
                }
                // A continuing feature receives its matching sloped course
                // immediately after the finishing skin. Dependency-aware
                // deferral removes the ordinary feature course from the local
                // obstacle field, while unrelated geometry remains a hard
                // rejection instead of moving the surface to an early slot.
                // Route the surviving fragments by their actual endpoints.
                // Clipping around a chimney or cabin can create more than one
                // fragment per raster row, for which a simple alternating row
                // order still produces many full-width travels. A nearest-end
                // route is deterministic, generic, and is applied identically
                // to every transition shell so corresponding paths remain
                // aligned through the complete non-planar skin.
                std::vector<size_t> ordered_path_indices;
                ordered_path_indices.reserve(clear_path_count);
                if (native_foundation) {
                    // Course corrections and wall routing run after the
                    // initial support test. Recheck the geometry we will emit
                    // and ensure its planned thickness still matches the
                    // native flat foundation before committing replacement.
                    const GeneratedShell &first = shells[1];
                    const double owner_z = layers[first.owner_layer]->print_z;
                    bool supported = native_foundation_supports_course(
                        first.paths, *layers[transition_anchor_layer]->regions()[region_idx],
                        owner_z, layers[transition_anchor_layer]->print_z,
                        minimum_course_height, maximum_course_height);
                    for (const ExtrusionEntity *entity : first.paths) {
                        std::vector<const ExtrusionPath *> paths;
                        collect_extrusion_path_references(*entity, paths);
                        for (const ExtrusionPath *path : paths)
                            for (const Point3 &point : path->polyline.points) {
                                const auto mesh_z = patch_surface_z_at_xy(patch, mesh, point.to_point());
                                if (!mesh_z || std::abs(owner_z + unscale_(point.z()) -
                                    (course_plan[1].plane_z + course_plan[1].blend *
                                     (*mesh_z - mesh.ground_level() - maximum_surface_z))) > 0.001)
                                    supported = false;
                            }
                    }
                    if (!supported) {
                        for (GeneratedShell &shell : shells)
                            destroy_entities(shell.paths);
                        enqueue_top_cap(patch, patch_idx, "native foundation does not support the corrected first course");
                        continue;
                    }
                }
                std::vector<uint8_t> path_unvisited = clear_at_every_depth;
                Point current_endpoint;
                bool have_current_endpoint = false;
                while (ordered_path_indices.size() < clear_path_count) {
                    size_t closest_idx = size_t(-1);
                    bool reverse_closest = false;
                    double closest_distance = std::numeric_limits<double>::max();
                    for (size_t path_idx = 0; path_idx < path_count; ++path_idx) {
                        if (!path_unvisited[path_idx])
                            continue;
                        const auto *path = dynamic_cast<const ExtrusionPath *>(
                            shells.back().paths[path_idx]);
                        if (path == nullptr)
                            continue;
                        if (!have_current_endpoint) {
                            closest_idx = path_idx;
                            break;
                        }
                        const double forward_distance =
                            (current_endpoint - path->first_point3().to_point()).squaredNorm();
                        const double reverse_distance =
                            (current_endpoint - path->last_point3().to_point()).squaredNorm();
                        const bool reverse = reverse_distance < forward_distance;
                        const double distance = reverse ? reverse_distance : forward_distance;
                        if (distance < closest_distance) {
                            closest_distance = distance;
                            closest_idx = path_idx;
                            reverse_closest = reverse;
                        }
                    }
                    if (closest_idx == size_t(-1))
                        break;
                    if (reverse_closest)
                        for (GeneratedShell &shell : shells)
                            if (auto *path = dynamic_cast<ExtrusionPath *>(
                                    shell.paths[closest_idx]))
                                path->polyline.reverse();
                    path_unvisited[closest_idx] = 0;
                    ordered_path_indices.push_back(closest_idx);
                    const auto *path = dynamic_cast<const ExtrusionPath *>(
                        shells.back().paths[closest_idx]);
                    current_endpoint = path->last_point3().to_point();
                    have_current_endpoint = true;
                }

                // make_patch_raster() produces outlines and top fill as one
                // aligned shell so they can be clearance-validated
                // atomically. Restore the same user-selected ordering used by
                // the conventional perimeter generator before scheduling the
                // shell. The nearest-end pass above remains the stable
                // tie-breaker for disconnected fragments at the same inset.
                //
                // Inner/Outer/Inner uses the normal Orca sandwich sequence:
                // walls from the third inward are deposited inside-out, then
                // the outer wall, and finally the first inner wall. With fewer
                // than three walls it has the conventional Outer/Inner
                // fallback. Infill-first is orthogonal to wall sequence and
                // is honored here as it is for planar regions.
                WallSequence wall_sequence = region_config.wall_sequence.value;
                if (region_config.wall_loops.value >= 3 &&
                    region_config.inner_walls_flow_ratio.get_abs_value(1.) >
                        1.0 + EPSILON &&
                    surface_owner > 0)
                    wall_sequence = WallSequence::InnerOuterInner;
                const bool sandwich_sequence =
                    wall_sequence == WallSequence::InnerOuterInner &&
                    region_config.wall_loops.value >= 3 && surface_owner > 0;
                const bool infill_first = region_config.is_infill_first.value;
                const auto path_order_key = [&](size_t path_idx) {
                    const auto *path = dynamic_cast<const ExtrusionPath *>(
                        shells.back().paths[path_idx]);
                    const bool wall = path != nullptr &&
                        (path->role() == erExternalPerimeter ||
                         path->role() == erPerimeter);
                    const int family = wall == infill_first ? 1 : 0;
                    if (!wall)
                        return std::pair<int, int>{family, 0};
                    const int inset = std::max(0, path->inset_idx);
                    if (sandwich_sequence) {
                        if (inset >= 2)
                            return std::pair<int, int>{family, -inset};
                        return std::pair<int, int>{family, inset == 0 ? 0 : 1};
                    }
                    if (wall_sequence == WallSequence::OuterInner ||
                        wall_sequence == WallSequence::InnerOuterInner)
                        return std::pair<int, int>{family, inset};
                    return std::pair<int, int>{family, -inset};
                };
                std::stable_sort(
                    ordered_path_indices.begin(), ordered_path_indices.end(),
                    [&](size_t left, size_t right) {
                        return path_order_key(left) < path_order_key(right);
                    });

                // Only the exposed course is a top skin. Every lower generated
                // course is a geometric transition and must retain the normal
                // structural density. Repeating the final top raster for the
                // configured transition count turns the whole transition into
                // solid plastic, consumes excessive material and carves away
                // ordinary infill which should remain between its sparse rows.
                // The same row phase is used on every deep course so each
                // retained bead has direct support. Conventional material is
                // carved from the model only after these actual emitted paths
                // are known, so omitted rows leave the ordinary sliced infill
                // in place rather than creating holes.
                const size_t solid_shell_begin = surface_support_begin;

                // Join neighboring, collision-safe scanlines into continuous
                // serpentine paths. The connector is accepted only when it
                // remains inside the same surface patch and clears unrelated
                // geometry at every shell depth. This removes a retract,
                // clearance lift, and full-width travel for each raster line
                // without weakening the generic gantry check.
                std::vector<uint8_t> merge_after(
                    ordered_path_indices.size() > 0 ? ordered_path_indices.size() - 1 : 0,
                    uint8_t(0));
                const coord_t maximum_connector = scale_(
                    1.5 * std::max(double(top_template->width), resolution));
                const double maximum_connector_squared =
                    double(maximum_connector) * double(maximum_connector);
                const double maximum_connector_slope = std::tan(
                    std::clamp(double(region_config.nonplanar_top_surface_max_angle.value),
                               0., 89.) * M_PI / 180.);
                for (size_t order_idx = 0;
                     order_idx + 1 < ordered_path_indices.size(); ++order_idx) {
                    const size_t path_idx = ordered_path_indices[order_idx];
                    const size_t next_path_idx = ordered_path_indices[order_idx + 1];
                    bool connector_is_clear = true;
                    for (GeneratedShell &shell : shells) {
                        if (!is_emitted(shell.paths[path_idx]) || !is_emitted(shell.paths[next_path_idx]))
                            continue;
                        auto *before = dynamic_cast<ExtrusionPath *>(shell.paths[path_idx]);
                        auto *after = dynamic_cast<ExtrusionPath *>(shell.paths[next_path_idx]);
                        if (before == nullptr || after == nullptr ||
                            before->role() != after->role() ||
                            before->inset_idx != after->inset_idx ||
                            std::abs(before->height - after->height) > EPSILON ||
                            before->first_point3() == before->last_point3() ||
                            after->first_point3() == after->last_point3() ||
                            (before->last_point3().to_point() -
                             after->first_point3().to_point()).squaredNorm() >
                                maximum_connector_squared) {
                            connector_is_clear = false;
                            break;
                        }
                        // Variable-width Arachne pieces still form one wall.
                        // Preserve their individual flow in a multipath, but
                        // only bridge endpoint-scale gaps, not separate walls.
                        if (std::abs(before->width - after->width) > EPSILON ||
                            std::abs(before->mm3_per_mm - after->mm3_per_mm) > EPSILON) {
                            const Vec3d gap = (after->first_point3() - before->last_point3())
                                .cast<double>() * SCALING_FACTOR;
                            if (!is_perimeter(before->role()) || gap.norm() > 0.05 ||
                                std::abs(gap.z()) > maximum_connector_slope * gap.head<2>().norm() + EPSILON) {
                                connector_is_clear = false;
                                break;
                            }
                        }
                        const Line connector(
                            before->last_point3().to_point(),
                            after->first_point3().to_point());
                        if (!std::any_of(
                                patch.projection.begin(), patch.projection.end(),
                                [&connector](const ExPolygon &polygon) {
                                    return polygon.contains(connector);
                                })) {
                            connector_is_clear = false;
                            break;
                        }
                        ExtrusionPath connector_path = *before;
                        connector_path.polyline.points = {
                            before->last_point3(), after->first_point3()};
                        connector_path.nonplanar_clearance_validated = false;
                        const size_t course = size_t(&shell - shells.data());
                        if (course > 0 && is_emitted(before) && is_emitted(after)) {
                            const ExtrusionPaths parts = course_flow(course, connector_path);
                            if (!std::all_of(parts.begin(), parts.end(), [&](const ExtrusionPath &part) {
                                    return nonplanar_course_is_supported(part, course_supports[course],
                                        layers[shell.owner_layer]->print_z, support_owner_z[course]);
                                })) {
                                connector_is_clear = false;
                                break;
                            }
                        }
                        if (native_foundation && &shell == &shells[1] &&
                            !native_foundation_supports_course(
                                ExtrusionEntitiesPtr{&connector_path},
                                *layers[transition_anchor_layer]->regions()[region_idx],
                                layers[shell.owner_layer]->print_z,
                                layers[transition_anchor_layer]->print_z,
                                minimum_course_height, maximum_course_height)) {
                            connector_is_clear = false;
                            break;
                        }
                        // The connector is an extrusion, not a travel.  A
                        // short connector can remain inside the patch and
                        // clear the toolhead while still cutting across one
                        // or more already routed raster rows at essentially
                        // the same physical Z.  This happened on the Benchy
                        // roof where one nominally short serpentine turn
                        // crossed thirteen top-surface rows.  Validate the
                        // connector against the complete shell, excluding
                        // only the two paths it joins; endpoint contact with
                        // either neighbor is intentionally permitted by
                        // paths_have_proper_same_z_crossing().
                        for (size_t other_idx = 0;
                             other_idx < shell.paths.size(); ++other_idx) {
                            if (other_idx == path_idx ||
                                other_idx == next_path_idx)
                                continue;
                            const auto *other = dynamic_cast<const ExtrusionPath *>(
                                shell.paths[other_idx]);
                            if (other != nullptr &&
                                paths_have_proper_same_z_crossing(
                                    connector_path,
                                    layers[shell.owner_layer]->print_z,
                                    *other,
                                    layers[shell.owner_layer]->print_z,
                                    0.05)) {
                                connector_is_clear = false;
                                break;
                            }
                        }
                        if (!connector_is_clear)
                            break;
                        const auto &connector_obstacles =
                            *remote_obstacle_cache[shell.owner_layer];
                        if (!generated_path_is_clear(
                                connector_path, layers[shell.owner_layer]->print_z,
                                gantry, fallback, connector_obstacles)) {
                            connector_is_clear = false;
                            break;
                        }
                    }
                    merge_after[order_idx] = connector_is_clear;
                }

                Polygons patch_coverage_polygons;
                for (size_t path_idx = 0; path_idx < path_count; ++path_idx)
                    if (clear_at_every_depth[path_idx]) {
                        auto *finishing_path = dynamic_cast<ExtrusionPath *>(
                            shells.back().paths[path_idx]);
                        finishing_path->polygons_covered_by_spacing(
                            patch_coverage_polygons, float(scale_(0.02)));
                    }
                // A merged serpentine extrudes the short segment between two
                // neighboring raster rows as well. Include those accepted
                // connectors in the replacement footprint before carving the
                // conventional shell. Otherwise internal solid/bridge paths
                // may survive underneath a connector and be printed first;
                // a later transition shell then descends through that bead.
                // Shell rasters are aligned at every depth, so the finishing
                // shell's XY connector footprint covers the complete stack.
                for (size_t order_idx = 0; order_idx < merge_after.size(); ++order_idx) {
                    if (!merge_after[order_idx])
                        continue;
                    const auto *before = dynamic_cast<const ExtrusionPath *>(
                        shells.back().paths[ordered_path_indices[order_idx]]);
                    const auto *after = dynamic_cast<const ExtrusionPath *>(
                        shells.back().paths[ordered_path_indices[order_idx + 1]]);
                    if (before == nullptr || after == nullptr)
                        continue;
                    ExtrusionPath connector = *before;
                    connector.polyline.points = {
                        before->last_point3(), after->first_point3()};
                    connector.polygons_covered_by_spacing(
                        patch_coverage_polygons, float(scale_(0.02)));
                }
                for (size_t shell_idx = 0; shell_idx < shells.size(); ++shell_idx) {
                    GeneratedShell &shell = shells[shell_idx];
                    if (native_foundation && shell_idx == 0)
                        continue;
                    ExtrusionPath *chain = nullptr;
                    ExtrusionPaths wall_parts;
                    std::optional<size_t> previous_order_idx;
                    const auto flush_chain = [&]() {
                        if (chain == nullptr)
                            return;
                        chain->nonplanar_clearance_validated = true;
                        wall_parts.emplace_back(std::move(*chain));
                        delete chain;
                        ExtrusionPaths parts;
                        for (const ExtrusionPath &part : wall_parts)
                            append(parts, course_flow(shell_idx, part));
                        wall_parts.clear();
                        if (parts.empty())
                            throw SlicingError("Non-planar course has no printable path");
                        const int inset = parts.front().inset_idx;
                        if (is_perimeter(parts.front().role()) &&
                            parts.front().first_point3() == parts.back().last_point3()) {
                            auto *loop = new ExtrusionLoop(std::move(parts));
                            loop->inset_idx = inset;
                            generated_by_layer[shell.owner_layer].emplace_back(loop);
                        } else if (parts.size() > 1) {
                            auto *multipath = new ExtrusionMultiPath(std::move(parts));
                            multipath->inset_idx = inset;
                            generated_by_layer[shell.owner_layer].emplace_back(multipath);
                        } else {
                            generated_by_layer[shell.owner_layer].emplace_back(new ExtrusionPath(std::move(parts.front())));
                        }
                        chain = nullptr;
                        first_operation_layer = std::min(
                            first_operation_layer, shell.owner_layer);
                        last_operation_layer = std::max(
                            last_operation_layer, shell.owner_layer);
                    };
                    for (size_t order_idx = 0;
                         order_idx < ordered_path_indices.size(); ++order_idx) {
                        const size_t path_idx = ordered_path_indices[order_idx];
                        auto *path = dynamic_cast<ExtrusionPath *>(shell.paths[path_idx]);
                        shell.paths[path_idx] = nullptr;
                        if (!is_emitted(path)) {
                            flush_chain();
                            previous_order_idx.reset();
                            delete path;
                            continue;
                        }
                        if (shell_idx + 1 < shells.size() &&
                            !is_perimeter(path->role()))
                            path->set_extrusion_role(
                                shell_idx < solid_shell_begin ?
                                    erInternalInfill : erSolidInfill);
                        if (chain != nullptr && previous_order_idx &&
                            *previous_order_idx + 1 == order_idx &&
                            merge_after[*previous_order_idx]) {
                            const size_t original_point_count =
                                chain->polyline.points.size();
                            chain->polyline.points.insert(
                                chain->polyline.points.end(),
                                path->polyline.points.begin(),
                                path->polyline.points.end());
                            chain->polyline.fitting_result.clear();
                            // Individual raster rows and their short joining
                            // connector may each be valid while the assembled
                            // serpentine still folds across an earlier row.
                            // Validate the entity that will actually reach the
                            // G-code writer.  If the join is unsafe, retain all
                            // extrusion by splitting the chain at this point;
                            // routing will insert a travel between the two
                            // genuinely disconnected paths.
                            ExtrusionPath whole_wall;
                            const ExtrusionPath *validation_path = chain;
                            if (!wall_parts.empty()) {
                                whole_wall = wall_parts.front();
                                for (size_t part = 1; part < wall_parts.size(); ++part)
                                    whole_wall.polyline.points.insert(whole_wall.polyline.points.end(),
                                        wall_parts[part].polyline.points.begin() + 1,
                                        wall_parts[part].polyline.points.end());
                                whole_wall.polyline.points.insert(whole_wall.polyline.points.end(),
                                    chain->polyline.points.begin() + 1, chain->polyline.points.end());
                                validation_path = &whole_wall;
                            }
                            if (path_has_proper_same_z_self_crossing(
                                    *validation_path,
                                    layers[shell.owner_layer]->print_z,
                                    0.05)) {
                                chain->polyline.points.resize(
                                    original_point_count);
                                flush_chain();
                                chain = path;
                            } else if (std::abs(chain->width - path->width) > EPSILON ||
                                       std::abs(chain->mm3_per_mm - path->mm3_per_mm) > EPSILON) {
                                chain->polyline.points.resize(original_point_count);
                                const Point3 endpoint = chain->last_point3();
                                chain->nonplanar_clearance_validated = true;
                                wall_parts.emplace_back(std::move(*chain));
                                delete chain;
                                chain = path;
                                if (endpoint != chain->first_point3())
                                    chain->polyline.points.insert(chain->polyline.points.begin(), endpoint);
                                chain->polyline.fitting_result.clear();
                            } else {
                                delete path;
                            }
                        } else {
                            flush_chain();
                            chain = path;
                        }
                        previous_order_idx = order_idx;
                    }
                    flush_chain();
                }
                for (GeneratedShell &shell : shells)
                {
                    for (ExtrusionEntity *entity : shell.paths)
                        delete entity;
                }

                if (perimeter_only_patch)
                    for (size_t wall_idx = 0;
                         wall_idx < retained_walls.size(); ++wall_idx) {
                        const auto &[path, source_idx] = retained_walls[wall_idx];
                        (void) source_idx;
                        if (retained_wall_crosses_raster[wall_idx])
                            continue;
                        path->nonplanar_clearance_validated = true;
                        path->polygons_covered_by_spacing(
                            patch_coverage_polygons, float(scale_(0.02)));
                    }
                ExPolygons patch_coverage = union_ex(patch_coverage_polygons);
                if (!patch_coverage.empty()) {
                    restore_wall_graph.reset();
                    if (!perimeter_only_patch) {
                        for (size_t wall_idx = 0;
                             wall_idx < retained_walls.size(); ++wall_idx) {
                            const auto &[path, source_idx] = retained_walls[wall_idx];
                            if (retained_wall_crosses_raster[wall_idx])
                                continue;
                            path->nonplanar_clearance_validated = true;
                            path->polygons_covered_by_spacing(
                                patch_coverage_polygons, float(scale_(0.02)));
                        }
                        patch_coverage = union_ex(patch_coverage_polygons);
                    }
                    BOOST_LOG_TRIVIAL(debug) << "Accepted nonplanar patch " << patch_idx
                                             << " with "
                                             << std::count(clear_at_every_depth.begin(),
                                                           clear_at_every_depth.end(), uint8_t(1))
                                             << " continuous raster paths";
                    if (ironing_template && !perimeter_only_patch) {
                        const size_t owner_layer = shells.back().owner_layer;
                        ExtrusionEntitiesPtr ironing_paths = make_patch_raster(
                            patch, mesh, *ironing_template,
                            layers[owner_layer]->print_z, resolution, gantry,
                            region_config.top_surface_pattern.value, 0.,
                            feature_footprint.empty() ? nullptr : &feature_footprint);
                        if (!remote_obstacle_cache[owner_layer])
                            remote_obstacle_cache[owner_layer] =
                                obstacle_prefix(owner_layer);
                        for (ExtrusionEntity *entity : ironing_paths) {
                            auto *path = dynamic_cast<ExtrusionPath *>(entity);
                            const auto &ironing_remote_obstacles =
                                *remote_obstacle_cache[owner_layer];
                            // make_patch_raster() derives ironing from the
                            // exact same reachable patch as the complete,
                            // already accepted finishing shell. Re-testing its
                            // dense narrow beads against a polygonized union of
                            // the top rasters introduced scaled-coordinate gaps
                            // and rejected every valid path. Patch membership
                            // is guaranteed by construction; retain the actual
                            // toolhead-clearance test for each ironing path.
                            if (path != nullptr && generated_path_is_clear(
                                    *path, layers[owner_layer]->print_z,
                                    gantry, fallback, ironing_remote_obstacles)) {
                                path->nonplanar_before_current_layer =
                                    patch_requires_early_emission;
                                generated_by_layer[owner_layer].emplace_back(entity);
                                last_operation_layer = std::max(
                                    last_operation_layer, owner_layer);
                            } else {
                                delete entity;
                            }
                        }
                    }
                    std::vector<FeatureCourse> feature_courses;
                    size_t feature_source_begin = surface_owner + 1;
                    if (!feature_footprint.empty()) {
                        const ExPolygons feature_membership_area = offset_ex(
                            feature_footprint,
                            scale_(0.5 * std::max({
                                nozzle_diameter, double(top_template->width),
                                resolution}) + 0.05));
                        for (size_t source_idx = transition_anchor_layer + size_t(native_foundation);
                             source_idx <= surface_owner; ++source_idx) {
                            const ExtrusionEntityCollection &perimeters =
                                layers[source_idx]->regions()[region_idx]->perimeters;
                            if (!std::any_of(
                                    perimeters.entities.begin(),
                                    perimeters.entities.end(),
                                    [&feature_membership_area](
                                        const ExtrusionEntity *entity) {
                                        return entity_has_feature_path(
                                            *entity, feature_membership_area);
                                    }))
                                continue;
                            feature_source_begin = source_idx;
                            break;
                        }
                    }
                    // An attached chimney/boss participates in the exact same
                    // physical course stack as the roof beneath it. A separate
                    // interpolation from the patch's global low/high Z range
                    // assigned the feature rings to later owner layers than
                    // the roof boundary they were meant to meet. The roof was
                    // consequently emitted first, unsupported beside the
                    // still-planar feature, and the delayed ring then crossed
                    // it. Reuse each generated shell's authoritative owner,
                    // depth and blend so both sides of the interface are
                    // deposited together at the same physical height.
                    // A mesh contact may begin part-way through the roof's
                    // transition stack. Starting before its first real sliced
                    // perimeter makes the extracted course empty and causes
                    // the complete feature run (including its outer wall) to
                    // be discarded. Use every available physical feature
                    // course, but never synthesize one below that first real
                    // cross-section.
                    const size_t available_feature_courses =
                        feature_source_begin <= surface_owner ?
                            surface_owner - feature_source_begin + 1 : 0;
                    const size_t feature_shell_count = std::min(
                        shells.size() - 1, available_feature_courses);
                    feature_courses.reserve(feature_shell_count);
                    const size_t first_feature_shell = std::max<size_t>(
                        1, shells.size() - feature_shell_count);
                    for (size_t feature_idx = 0;
                         feature_idx < feature_shell_count; ++feature_idx) {
                        const GeneratedShell &shell =
                            shells[first_feature_shell + feature_idx];
                        feature_courses.push_back({
                            shell.owner_layer, maximum_surface_z,
                            shell.surface_depth,
                            shell.blend});
                    }
                    if (!native_foundation && transition_anchor_layer > 0)
                        flatten_interlocking_wall_foundation(
                            layers[transition_anchor_layer - 1]
                                ->regions()[region_idx]->perimeters,
                            wall_selection);
                    accepted_patches.push_back({
                        patch_idx, std::move(patch_coverage),
                        std::move(feature_footprint), std::move(feature_courses),
                        surface_owner, transition_anchor_layer,
                        feature_source_begin,
                        course_plan[1].blend,
                        generated_surface_owns_boundary,
                        patch_requires_early_emission, native_foundation});
                    emitted_patch_projections.push_back({
                        patch.projection, patch.minimum_z, patch.maximum_z});
                    for (const size_t face : patch.face_indices)
                        emitted_facets[face] = 1;
                }
            }

            if (accepted_patches.empty()) {
                // Restore the original supported slice when no complete patch
                // is reachable. Partial flattening or trimming at G-code time
                // would leave flagged variable-Z paths in otherwise planar
                // geometry and can produce the same holes as a partial roof.
                for (size_t source_idx = run_begin; source_idx <= run_end;
                     ++source_idx) {
                    LayerRegion &source =
                        *layers[source_idx]->regions()[region_idx];
                    restore_nonplanar_as_planar(source.fills);
                    restore_nonplanar_as_planar(source.perimeters);
                }
                layer_idx = run_end + 1;
                continue;
            }

            // Restore the ordinary shell first, then carve only the local
            // surface-height band occupied by the generated draped shell.
            // This preserves solid infill below the transition and removes
            // every same-height perimeter/fill that could cross a draped bead.
            // Both operands below originally described centerline ownership:
            // a conventional wall whose centerline was just outside either
            // polygon could still put half of its deposited bead into the
            // variable-Z shell. Expand the accepted bead envelope and the
            // locally active height band by the largest participating bead
            // radius, leaving 0.02 mm of intentional contact. Clipping path
            // centerlines against that Minkowski sum removes deep overlaps
            // without opening a seam at the replacement boundary.
            const double maximum_conventional_width = std::max({
                double(top_template->width),
                outline_template ? double(outline_template->width) : 0.,
                inner_outline_template ? double(inner_outline_template->width) : 0.});
            // The surviving conventional bead must reach the generated bead
            // envelope at every replacement boundary.  Top-cap fallbacks
            // previously added another 0.20 mm of empty clearance here. That
            // avoided envelope overlap, but it also isolated the cap from the
            // surrounding shell and could leave a visible/weak lip. Clip the
            // conventional centerline by one half of its largest possible
            // width minus a small bonding overlap, identically for complete
            // patches and local top caps. The generated variable-width
            // anchors already close residuals on the non-planar side.
            // End-on joints need more than a numerical overlap to survive
            // path simplification and to render as a genuinely connected
            // bead envelope. This remains far below one line width and is
            // applied only where the conventional path is clipped against
            // the accepted replacement boundary, so it cannot create a
            // second parallel wall or fill crossing.
            constexpr double boundary_bonding_overlap = 0.05;
            const coord_t replacement_bead_margin = scale_(std::max(
                0., 0.5 * maximum_conventional_width -
                        boundary_bonding_overlap));
            std::vector<coord_t> replacement_bead_margins;
            replacement_bead_margins.reserve(accepted_patches.size());
            std::vector<ExPolygons> accepted_bead_clearance;
            accepted_bead_clearance.reserve(accepted_patches.size());
            for (const AcceptedPatch &accepted : accepted_patches) {
                replacement_bead_margins.emplace_back(
                    replacement_bead_margin);
                accepted_bead_clearance.emplace_back(
                    offset_ex(accepted.coverage,
                              replacement_bead_margin));
            }

            // Only the paths that formed this candidate run are provisional.
            // Earlier accepted runs may live below this one (for example a
            // Benchy roof below the chimney cap); restoring every preceding
            // layer here silently flattened and untagged that already-valid
            // non-planar surface when a later candidate was processed.
            for (size_t source_idx = run_begin; source_idx <= run_end;
                 ++source_idx) {
                LayerRegion &source =
                    *layers[source_idx]->regions()[region_idx];
                restore_nonplanar_as_planar(source.fills);
                restore_nonplanar_as_planar(source.perimeters);
            }
            // Snapshot complete continuing-feature wall graphs before the
            // height-band replacement pass clips conventional stair walls.
            // Deferral is intentionally performed after shell commit, so
            // extracting only from the post-clipped graph can otherwise lose
            // one or more configured insets. The live graph is still
            // extracted later to avoid duplicates; these snapshots provide
            // the atomic, pre-clip source for projection and rollback.
            std::vector<std::vector<ExtrusionEntitiesPtr>>
                feature_course_snapshots(
                    accepted_patches.size(),
                    std::vector<ExtrusionEntitiesPtr>(layers.size()));
            for (size_t accepted_idx = 0;
                 accepted_idx < accepted_patches.size(); ++accepted_idx) {
                const AcceptedPatch &accepted = accepted_patches[accepted_idx];
                if (accepted.feature_footprint.empty() ||
                    accepted.feature_source_begin > accepted.surface_owner)
                    continue;
                const ExPolygons membership = offset_ex(
                    accepted.feature_footprint,
                    scale_(0.5 * std::max({
                        nozzle_diameter, double(top_template->width), resolution}) +
                        0.05));
                for (size_t source_idx = accepted.feature_source_begin;
                     source_idx <= accepted.surface_owner; ++source_idx) {
                    const ExtrusionEntityCollection &perimeters =
                        layers[source_idx]->regions()[region_idx]->perimeters;
                    ExtrusionEntityCollection snapshot_source;
                    snapshot_source.no_sort = perimeters.no_sort;
                    snapshot_source.inset_idx = perimeters.inset_idx;
                    for (const ExtrusionEntity *entity : perimeters.entities)
                        snapshot_source.entities.emplace_back(entity->clone());
                    extract_feature_entities(
                        snapshot_source, membership,
                        feature_course_snapshots[accepted_idx][source_idx]);
                }
            }
            std::vector<std::pair<const ExtrusionPath *, double>> finishing_paths;
            for (size_t owner = 0; owner < generated_by_layer.size(); ++owner)
                for (const ExtrusionEntity *entity : generated_by_layer[owner]) {
                    std::vector<const ExtrusionPath *> paths;
                    collect_extrusion_path_references(*entity, paths);
                    for (const ExtrusionPath *path : paths)
                        if (path->nonplanar_surface && !path->nonplanar_transition &&
                            !path->nonplanar_feature_transition &&
                            (path->role() == erTopSolidInfill || is_perimeter(path->role())))
                            finishing_paths.emplace_back(path, layers[owner]->print_z);
                }
            // Material replacement follows the exact height of the first
            // generated blend.  Clearance retains the wider physical shell
            // band used by the gantry checks, but it must never be reused to
            // decide which solid course is the transition foundation.
            for (size_t source_idx = 1;
                 source_idx <= last_operation_layer; ++source_idx) {
                LayerRegion &source = *layers[source_idx]->regions()[region_idx];
                const bool is_transition_support_layer = std::any_of(
                    accepted_patches.begin(), accepted_patches.end(),
                    [source_idx](const AcceptedPatch &accepted) {
                        return (accepted.transition_anchor_layer > 0 &&
                                accepted.transition_anchor_layer - 1 == source_idx) ||
                               (accepted.native_foundation &&
                                accepted.transition_anchor_layer == source_idx);
                    });
                ExPolygons replacement_fill_area;
                ExPolygons replacement_wall_area;
                ExPolygons replacement_clearance_area;
                for (size_t accepted_idx = 0;
                     accepted_idx < accepted_patches.size(); ++accepted_idx) {
                    const AcceptedPatch &accepted = accepted_patches[accepted_idx];
                    ExPolygons clearance_band = patch_shell_band_at_layer(
                        patches[accepted.patch_index], mesh,
                        layers[source_idx]->print_z, total_shell_depth,
                        layers[source_idx]->height,
                        layers[source_idx]->height);
                    if (!clearance_band.empty()) {
                        clearance_band = offset_ex(
                            clearance_band,
                            replacement_bead_margins[accepted_idx]);
                        append(replacement_clearance_area,
                               intersection_ex(
                                   clearance_band,
                                   accepted_bead_clearance[accepted_idx]));
                    }
                    // Do not remove every conventional wall throughout the
                    // patch's global low-to-high owner range.  A sloped course
                    // crosses many nominal Z planes, but most of those planar
                    // walls are its staircase support rather than duplicate
                    // boundary material.  The exact 3D bead-at-Z calculation
                    // below removes only centerlines which physically share a
                    // generated course.  The former broad subtraction deleted
                    // all three wall insets below the bounded roof stack and
                    // left its first course floating by several millimetres.
                }
                // Material replacement must be driven by the paths that will
                // actually be deposited, not by the complete mesh projection
                // of an idealized first shell.  On a shallow roof the latter
                // spans XY columns where the first generated course is still
                // one or more layers above the local staircase.  Carving that
                // theoretical band deleted those intervening source courses
                // and left the low edge of the Benchy roof floating by
                // 0.4--1.8 mm.  Project the emitted 3D bead envelope onto this
                // exact conventional center Z instead.  This retains every
                // course below the transition and removes only material that
                // would physically share the replacement course.
                constexpr double generated_centerline_z_tolerance = 0.05;
                append(replacement_fill_area, offset_ex(
                    generated_material_at_z(
                        generated_by_layer, layers,
                        layers[source_idx]->print_z,
                        generated_centerline_z_tolerance, true),
                    replacement_bead_margin));
                append(replacement_wall_area, offset_ex(
                    generated_material_at_z(
                        generated_by_layer, layers,
                        layers[source_idx]->print_z,
                        generated_centerline_z_tolerance),
                    replacement_bead_margin));
                // Retained perimeter graphs are contoured in place rather
                // than copied into generated_by_layer. Their physical bead
                // envelope is nevertheless replacement material at the same
                // absolute Z. Include it when clipping conventional wall and
                // top-fill paths, or a later source course can cross directly
                // through the accepted roof/chimney transition while the
                // non-planar wall itself is correctly preserved below.
                const ExPolygons retained_wall_material =
                    retained_nonplanar_wall_material_at_z(
                        layers, region_idx, layers[source_idx]->print_z,
                        generated_centerline_z_tolerance);
                append(replacement_fill_area, offset_ex(
                    retained_wall_material, replacement_bead_margin));
                append(replacement_wall_area, offset_ex(
                    retained_wall_material, replacement_bead_margin));
                // The final mesh-following skin replaces the rounded-up
                // planar cap above it, too. A symmetric center-Z comparison
                // misses that cap when its nozzle plane lies more than the
                // join tolerance above the true surface. Hybrid contouring
                // subsequently lowers the leftover cap into the finished
                // wall. Use a one-sided layer interval for finishing paths
                // only; never remove the supporting course below the skin.
                Polygons superseded_cap;
                const double half_height = 0.5 * layers[source_idx]->height;
                for (const auto &[path, owner_z] : finishing_paths)
                    collect_generated_material_at_z(*path, owner_z,
                        layers[source_idx]->print_z - half_height,
                        half_height, false, superseded_cap);
                const ExPolygons cap_coverage = offset_ex(
                    union_ex(superseded_cap), replacement_bead_margin);
                append(replacement_wall_area, cap_coverage);
                append(replacement_fill_area, cap_coverage);
                replacement_fill_area = union_ex(replacement_fill_area);
                replacement_wall_area = union_ex(replacement_wall_area);
                replacement_clearance_area =
                    union_ex(replacement_clearance_area);
                // The course immediately below a patch's flat anchor is the
                // support that was used by its bead-to-bead acceptance test.
                // Never carve it after acceptance. Generated-material-at-Z
                // may include a retained provisional wall descending into
                // this nominal plane; using that later path as a replacement
                // mask deleted the already validated support and turned an
                // initially sound foundation into a floating G-code ring.
                // It is one normal layer below the generated foundation, so
                // preserving the complete conventional course is the same
                // dependency rule used by ordinary top-solid printing.
                if (!is_transition_support_layer &&
                    !replacement_fill_area.empty())
                    remove_nonplanar_top_overlap(
                        source.fills, replacement_fill_area);
                if (!replacement_clearance_area.empty()) {
                    // Clearance is a toolhead validation envelope, not a
                    // material-removal envelope. In the planar Benchy control
                    // slice, an internal bridge one course below the first
                    // draped roof pass is its actual foundation. Removing all
                    // bridge/sparse paths throughout this wider band carved
                    // that foundation away after the exact replacement pass
                    // above had deliberately retained it. Conventional fill
                    // is therefore clipped only by replacement_fill_area;
                    // generated_path_is_clear() remains responsible for the
                    // larger gantry/nozzle envelope.
                    // A draped shell is owned by the accepted surface, not by
                    // one nominal source layer.  Its lowest points may
                    // descend through a conventional stair course whose
                    // owner differs from every generated shell owner (for
                    // example the Benchy bow at 18.8 mm beneath a shell owned
                    // by 19.0 mm).  Testing only generated_by_layer[source]
                    // retained that stair wall and let the later XYZ raster
                    // cross it at essentially the same Z. Remove conventional
                    // wall centerlines throughout the same local height band
                    // as fill. Material outside replacement_area is preserved,
                    // leaving the planar wall fragments attached to the
                    // connected retained boundary instead of closing beneath
                    // the generated skin.
                    // Perimeters are deposited material, so clip them only
                    // where a generated wall actually replaces their bead.
                    // The larger clearance band is for sparse/bridge paths
                    // and toolhead travel; applying it to walls deleted the
                    // planar foundation beneath the first draped course and
                    // left real Benchy roof walls 0.4--0.6 mm unsupported.
                    // The exact first-shell footprint keeps the immediately
                    // preceding wall course while still removing stair walls
                    // that share or overlie the generated shell's physical Z.
                    if (!is_transition_support_layer &&
                        !replacement_wall_area.empty())
                        remove_nonplanar_top_overlap(
                            source.perimeters, replacement_wall_area);
                }
                WallSequence wall_sequence = region_config.wall_sequence.value;
                if (region_config.wall_loops.value >= 3 &&
                    region_config.inner_walls_flow_ratio.get_abs_value(1.) >
                        1.0 + EPSILON &&
                    source_idx > 0)
                    wall_sequence = WallSequence::InnerOuterInner;
                restore_nonplanar_wall_sequence(
                    source.perimeters, wall_sequence,
                    region_config.wall_loops.value, source_idx > 0);
            }
            for (size_t destination_idx = first_operation_layer;
                 destination_idx <= last_operation_layer; ++destination_idx)
            {
                report_progress(
                    double(destination_idx), NonplanarProgressStage::Commit,
                    destination_idx - first_operation_layer + 1,
                    last_operation_layer - first_operation_layer + 1);
                LayerRegion &destination =
                    *layers[destination_idx]->regions()[region_idx];
                ExtrusionEntitiesPtr generated_walls;
                ExtrusionEntitiesPtr generated_fills;
                for (ExtrusionEntity *entity : generated_by_layer[destination_idx]) {
                    bool wall = false;
                    if (const auto *path = dynamic_cast<const ExtrusionPath *>(entity))
                        wall = is_perimeter(path->role());
                    else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(entity))
                        wall = !multipath->paths.empty() && std::all_of(
                            multipath->paths.begin(), multipath->paths.end(),
                            [](const ExtrusionPath &path) { return is_perimeter(path.role()); });
                    else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(entity))
                        wall = !loop->paths.empty() && std::all_of(
                            loop->paths.begin(), loop->paths.end(),
                            [](const ExtrusionPath &path) {
                                return is_perimeter(path.role());
                            });
                    if (wall)
                        generated_walls.emplace_back(entity);
                    else
                        generated_fills.emplace_back(entity);
                }
                generated_by_layer[destination_idx].clear();
                append_ordered_entities(
                    destination.perimeters, std::move(generated_walls), 0);
                append_ordered_entities(
                    destination.fills, std::move(generated_fills), 0);
                WallSequence wall_sequence = region_config.wall_sequence.value;
                if (region_config.wall_loops.value >= 3 &&
                    region_config.inner_walls_flow_ratio.get_abs_value(1.) >
                        1.0 + EPSILON && destination_idx > 0)
                    wall_sequence = WallSequence::InnerOuterInner;
                restore_nonplanar_wall_sequence(
                    destination.perimeters, wall_sequence,
                    region_config.wall_loops.value, destination_idx > 0);
            }

            // A feature meeting a draped skin cannot already be standing at
            // its conventional layer heights when the skin arrives. Defer its
            // intersecting courses, print the first courses on the exact same
            // height field as the replacement shells, then consume the
            // remaining courses above the final skin while continuously
            // reducing the slope. The last deferred course is horizontal at
            // the destination layer, so ordinary slicing resumes without a Z
            // discontinuity. Connectivity and height alone drive this; there
            // are no model- or printer-specific feature names.
            for (size_t accepted_patch_idx = 0;
                 accepted_patch_idx < accepted_patches.size();
                 ++accepted_patch_idx) {
                const AcceptedPatch &accepted =
                    accepted_patches[accepted_patch_idx];
                if (accepted.feature_footprint.empty())
                    continue;
                report_progress(
                    double(run_end), NonplanarProgressStage::FeatureScheduling,
                    accepted.patch_index + 1, patches.size());
                const NonplanarPatch &patch = patches[accepted.patch_index];
                const std::optional<std::pair<double, double>> surface_range =
                    feature_surface_z_range(
                        accepted.feature_footprint, patch, mesh);
                if (!surface_range)
                    continue;

                const double local_surface_z = surface_range->second;
                // Build the feature-membership tolerance once. Previously it
                // was offset again for every Arachne fragment, which made a
                // real Benchy spend most of dependency deferral inside
                // identical ClipperOffset calls.
                const ExPolygons feature_membership_area = offset_ex(
                    accepted.feature_footprint,
                    scale_(0.5 * std::max({
                        nozzle_diameter, double(top_template->width), resolution}) +
                        0.05));
                size_t first_source = accepted.feature_source_begin;
                // The main replacement pass may consume a nominal stair
                // perimeter which was present when the dependency was first
                // classified. Re-evaluate the surviving slicer-owned graph
                // before extraction and skip only leading empty courses. A
                // missing course in the middle remains an atomic rejection.
                while (first_source <= accepted.surface_owner) {
                    const ExtrusionEntityCollection &perimeters =
                        layers[first_source]->regions()[region_idx]->perimeters;
                    if (std::any_of(
                            perimeters.entities.begin(),
                            perimeters.entities.end(),
                            [&feature_membership_area](
                                const ExtrusionEntity *entity) {
                                return entity_has_feature_path(
                                    *entity, feature_membership_area);
                            }))
                        break;
                    ++first_source;
                }

                struct DeferredCourse {
                    size_t source_layer {0};
                    size_t destination_layer {0};
                    bool synthetic {false};
                    ExtrusionEntitiesPtr perimeters;
                    ExtrusionEntitiesPtr fills;
                };
                const size_t available_source_courses =
                    first_source <= accepted.surface_owner ?
                        accepted.surface_owner - first_source + 1 : 0;
                const size_t surface_course_count = std::min(
                    accepted.feature_courses.size(), available_source_courses);
                const size_t feature_course_offset =
                    accepted.feature_courses.size() - surface_course_count;
                if (surface_course_count == 0 ||
                    surface_course_count > accepted.surface_owner + 1)
                    continue;
                // Use feature cross-sections from the physical low side of
                // the transition, not the last N global layers below the
                // patch's highest owner. On a sloped roof an angled post or
                // flag socket changes XY as Z rises; taking its upper rings
                // and moving them down shifts the whole continuation away
                // from the roof attachment. Starting at first_source keeps
                // each warped course aligned with the model geometry at the
                // height it replaces.
                const size_t surface_source_begin = first_source;
                if (surface_source_begin + surface_course_count >
                    accepted.surface_owner + 1)
                    continue;
                const double resume_plane_z =
                    layers[accepted.surface_owner]->print_z;
                // The attachment may sit on the low side of a roof whose
                // global owner is selected by the high side. Build as many
                // physical-height courses as are needed to reach that owner
                // plane; assigning one course to each later global layer
                // would leave a large unsupported jump at the attachment.
                const size_t leveling_course_count = size_t(std::max(1.,
                    std::ceil(std::max(0., resume_plane_z - surface_range->first) /
                              shell_layer_height)));

                const auto extract_course = [&](size_t source_idx,
                                                size_t destination_idx) {
                    DeferredCourse course;
                    course.source_layer = source_idx;
                    course.destination_layer = destination_idx;
                    LayerRegion &source =
                        *layers[source_idx]->regions()[region_idx];
                    extract_feature_entities(
                        source.perimeters, feature_membership_area,
                        course.perimeters);
                    const ExtrusionEntitiesPtr &snapshot =
                        feature_course_snapshots[accepted_patch_idx][source_idx];
                    if (!snapshot.empty()) {
                        destroy_entities(course.perimeters);
                        course.perimeters.reserve(snapshot.size());
                        for (const ExtrusionEntity *entity : snapshot)
                            course.perimeters.emplace_back(entity->clone());
                    }
                    // Fill remains on its conventional course and is clipped
                    // later only where a successfully validated transformed
                    // wall actually occupies the same physical Z. Moving all
                    // fill intersecting this footprint would capture complete
                    // roof rasters and distort the supporting surface.
                    return course;
                };

                // Courses below the first replacement shell would leave the
                // continuing feature standing in the toolhead envelope before
                // the sloped skin arrives. Keep them recoverable until all
                // replacement courses pass projection and collision checks;
                // on success the underlying solid roof owns this volume.
                std::vector<DeferredCourse> superseded_courses;
                superseded_courses.reserve(
                    accepted.surface_owner + 1 - surface_course_count);
                for (size_t source_idx = first_source;
                     source_idx < surface_source_begin; ++source_idx)
                    superseded_courses.emplace_back(
                        extract_course(source_idx, source_idx));

                std::vector<DeferredCourse> courses;
                courses.reserve(surface_course_count + leveling_course_count);
                bool complete_course_run = true;
                for (size_t course_idx = 0;
                     course_idx < surface_course_count; ++course_idx) {
                    report_progress(
                        double(run_end), NonplanarProgressStage::FeatureProjection,
                        course_idx + 1,
                        surface_course_count + leveling_course_count);
                    DeferredCourse course = extract_course(
                        surface_source_begin + course_idx,
                        accepted.feature_courses[
                            feature_course_offset + course_idx].owner_layer);
                    if (course.perimeters.empty())
                        complete_course_run = false;
                    courses.emplace_back(std::move(course));
                }
                // The synthetic leveling run replaces the remaining nominal
                // courses up to the global owner plane. Remove those original
                // paths atomically as well; leaving them in place recreates a
                // shifted duplicate wall beside the correctly anchored
                // feature transition.
                for (size_t source_idx =
                         surface_source_begin + surface_course_count;
                     source_idx <= accepted.surface_owner; ++source_idx)
                    superseded_courses.emplace_back(
                        extract_course(source_idx, source_idx));
                if (complete_course_run && !courses.empty()) {
                    // Leveling courses are synthetic continuations of the
                    // last real feature cross-section. Normal feature paths
                    // on the following global layer remain untouched and are
                    // therefore the exact handoff target after flattening.
                    const DeferredCourse &leveling_template = courses.back();
                    for (size_t leveling_idx = 0;
                         leveling_idx < leveling_course_count; ++leveling_idx) {
                        DeferredCourse course;
                        course.source_layer = accepted.surface_owner;
                        course.destination_layer = accepted.surface_owner;
                        course.synthetic = true;
                        course.perimeters.reserve(
                            leveling_template.perimeters.size());
                        for (const ExtrusionEntity *entity :
                             leveling_template.perimeters)
                            course.perimeters.emplace_back(entity->clone());
                        course.fills.reserve(leveling_template.fills.size());
                        for (const ExtrusionEntity *entity :
                             leveling_template.fills)
                            course.fills.emplace_back(entity->clone());
                        courses.emplace_back(std::move(course));
                    }
                }

                // Preserve the original entities when the sliced feature has
                // a discontinuity or does not supply one course per shell.
                // Such a feature is not a safe candidate for deferral.
                if (!complete_course_run) {
                    for (DeferredCourse &course : superseded_courses) {
                        LayerRegion &source =
                            *layers[course.source_layer]->regions()[region_idx];
                        append_ordered_entities(
                            source.perimeters, std::move(course.perimeters), 0);
                        append_ordered_entities(
                            source.fills, std::move(course.fills), 0);
                    }
                    for (DeferredCourse &course : courses) {
                        LayerRegion &source =
                            *layers[course.source_layer]->regions()[region_idx];
                        append_ordered_entities(
                            source.perimeters, std::move(course.perimeters), 0);
                        append_ordered_entities(
                            source.fills, std::move(course.fills), 0);
                    }
                    BOOST_LOG_TRIVIAL(debug)
                        << "Nonplanar patch " << accepted.patch_index
                        << " kept its continuing feature planar: incomplete course run";
                    continue;
                }

                size_t warped_count = 0;
                bool projection_succeeded = true;
                for (size_t course_idx = 0;
                     course_idx < surface_course_count; ++course_idx) {
                    DeferredCourse &source = courses[course_idx];
                    const FeatureCourse &destination =
                        accepted.feature_courses[
                            feature_course_offset + course_idx];
                    const double owner_z = layers[destination.owner_layer]->print_z;
                    size_t course_warped = 0;
                    size_t expected_paths = 0;
                    for (ExtrusionEntity *entity : source.perimeters) {
                        expected_paths += count_warpable_feature_paths(*entity);
                        course_warped += warp_feature_entity_to_surface(
                            *entity, accepted.feature_footprint, patch, mesh,
                            owner_z,
                            destination.surface_reference_z -
                                destination.surface_depth,
                            destination.surface_reference_z,
                            destination.blend,
                            accepted.before_current_layer, course_idx + 1,
                            &layers, accepted.transition_anchor_layer);
                    }
                    for (ExtrusionEntity *entity : source.fills) {
                        expected_paths += count_warpable_feature_paths(*entity);
                        course_warped += warp_feature_entity_to_surface(
                            *entity, accepted.feature_footprint, patch, mesh,
                            owner_z,
                            destination.surface_reference_z -
                                destination.surface_depth,
                            destination.surface_reference_z,
                            destination.blend,
                            accepted.before_current_layer, course_idx + 1,
                            &layers, accepted.transition_anchor_layer);
                    }
                    projection_succeeded &= expected_paths > 0 &&
                        course_warped == expected_paths;
                    if (course_warped != expected_paths)
                        BOOST_LOG_TRIVIAL(debug)
                            << "Nonplanar patch " << accepted.patch_index
                            << " projected " << course_warped << " of "
                            << expected_paths << " paths in feature course "
                            << course_idx
                            << " from layer " << source.source_layer
                            << " to layer " << source.destination_layer;
                    warped_count += course_warped;
                }

                for (size_t leveling_idx = 0;
                     leveling_idx < leveling_course_count; ++leveling_idx) {
                    report_progress(
                        double(run_end), NonplanarProgressStage::FeatureProjection,
                        surface_course_count + leveling_idx + 1,
                        surface_course_count + leveling_course_count);
                    DeferredCourse &source =
                        courses[surface_course_count + leveling_idx];
                    const double ratio = double(leveling_idx + 1) /
                        double(leveling_course_count);
                    const double blend = 1. - ratio;
                    const double owner_z = resume_plane_z;
                    // Advance by approximately one physical layer while the
                    // slope fades. The final synthetic course is horizontal
                    // on the current owner plane; the untouched conventional
                    // feature on the next layer is then exactly one ordinary
                    // layer above it.
                    const double base_z = local_surface_z +
                        ratio * (resume_plane_z - local_surface_z);
                    size_t course_warped = 0;
                    size_t expected_paths = 0;
                    for (ExtrusionEntity *entity : source.perimeters) {
                        expected_paths += count_warpable_feature_paths(*entity);
                        course_warped += warp_feature_entity_to_surface(
                            *entity, accepted.feature_footprint, patch, mesh,
                            owner_z, base_z, local_surface_z, blend,
                            accepted.before_current_layer,
                            surface_course_count + leveling_idx + 1);
                        mark_nonplanar_leveling_transition(*entity);
                    }
                    for (ExtrusionEntity *entity : source.fills) {
                        expected_paths += count_warpable_feature_paths(*entity);
                        course_warped += warp_feature_entity_to_surface(
                            *entity, accepted.feature_footprint, patch, mesh,
                            owner_z, base_z, local_surface_z, blend,
                            accepted.before_current_layer,
                            surface_course_count + leveling_idx + 1);
                        mark_nonplanar_leveling_transition(*entity);
                    }
                    projection_succeeded &= expected_paths > 0 &&
                        course_warped == expected_paths;
                    if (course_warped != expected_paths)
                        BOOST_LOG_TRIVIAL(debug)
                            << "Nonplanar patch " << accepted.patch_index
                            << " projected " << course_warped << " of "
                            << expected_paths << " paths in leveling course "
                            << leveling_idx
                            << " from layer " << source.source_layer;
                    warped_count += course_warped;
                }

                // Adjacent feature cross-sections may shift laterally as a
                // boss emerges through a slope. Independent interpolation
                // can then make their centerlines cross at nearly the same
                // absolute Z even though they represent successive physical
                // courses. Constrain each course against the preceding one
                // at matching insets. This preserves the sliced XY contour
                // while maintaining one real bead height of vertical order.
                const auto collect_course_paths =
                    [](auto &&self, ExtrusionEntity &entity,
                       std::vector<ExtrusionPath *> &paths) -> void {
                        if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
                            if (is_perimeter(path->role()))
                                paths.push_back(path);
                        } else if (auto *multipath =
                                       dynamic_cast<ExtrusionMultiPath *>(&entity)) {
                            for (ExtrusionPath &path : multipath->paths)
                                self(self, path, paths);
                        } else if (auto *loop =
                                       dynamic_cast<ExtrusionLoop *>(&entity)) {
                            for (ExtrusionPath &path : loop->paths)
                                self(self, path, paths);
                        } else if (auto *collection =
                                       dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
                            for (ExtrusionEntity *child : collection->entities)
                                self(self, *child, paths);
                        }
                    };
                for (size_t course_idx = 1;
                     course_idx < courses.size(); ++course_idx) {
                    DeferredCourse &previous = courses[course_idx - 1];
                    DeferredCourse &current = courses[course_idx];
                    std::vector<ExtrusionPath *> previous_paths;
                    std::vector<ExtrusionPath *> current_paths;
                    for (ExtrusionEntity *entity : previous.perimeters)
                        collect_course_paths(
                            collect_course_paths, *entity, previous_paths);
                    for (ExtrusionEntity *entity : current.perimeters)
                        collect_course_paths(
                            collect_course_paths, *entity, current_paths);
                    const double previous_owner_z =
                        layers[previous.destination_layer]->print_z;
                    const double current_owner_z =
                        layers[current.destination_layer]->print_z;
                    for (ExtrusionPath *path : current_paths) {
                        bool changed = false;
                        for (Point3 &point : path->polyline.points) {
                            const Vec2d sample =
                                unscale(point.to_point()).cast<double>();
                            double supporting_z =
                                -std::numeric_limits<double>::infinity();
                            for (const ExtrusionPath *support : previous_paths) {
                                if (support->inset_idx != path->inset_idx)
                                    continue;
                                for (size_t point_idx = 1;
                                     point_idx < support->polyline.points.size();
                                     ++point_idx) {
                                    const Point3 &first =
                                        support->polyline.points[point_idx - 1];
                                    const Point3 &last =
                                        support->polyline.points[point_idx];
                                    const Vec2d a =
                                        unscale(first.to_point()).cast<double>();
                                    const Vec2d b =
                                        unscale(last.to_point()).cast<double>();
                                    const Vec2d delta = b - a;
                                    const double squared = delta.squaredNorm();
                                    const double ratio = squared <= EPSILON ? 0. :
                                        std::clamp((sample - a).dot(delta) /
                                                       squared,
                                                   0., 1.);
                                    if ((sample - (a + ratio * delta)).norm() >
                                        0.5 * (double(path->width) +
                                               double(support->width)) + 0.05)
                                        continue;
                                    const double z = previous_owner_z +
                                        unscale_(first.z()) + ratio *
                                            unscale_(last.z() - first.z());
                                    supporting_z = std::max(supporting_z, z);
                                }
                            }
                            if (!std::isfinite(supporting_z))
                                continue;
                            const double current_z = current_owner_z +
                                unscale_(point.z());
                            // This is an ordering/contact constraint, not a
                            // second layer-height generator. A full nominal
                            // layer offset can push a curved ring through the
                            // finishing raster. Stay just beyond the 0.05 mm
                            // same-Z collision tolerance while leaving the
                            // model-projected surface profile authoritative.
                            constexpr double minimum_course_clearance = 0.0501;
                            const double minimum_z = supporting_z +
                                minimum_course_clearance;
                            if (current_z + EPSILON < minimum_z) {
                                point.z() = scale_(minimum_z - current_owner_z);
                                changed = true;
                            }
                        }
                        if (changed)
                            path->polyline.fitting_result.clear();
                    }
                }

                // Keep the ordered feature stack on the wall side of the
                // already generated finishing skin. The small course-order
                // correction above can otherwise lift a curved wall by a few
                // microns through a nearby raster. Clamp only against a top
                // path above the wall and within one combined bead radius.
                struct OwnedTopPath {
                    const ExtrusionPath *path;
                    double owner_z;
                };
                std::vector<OwnedTopPath> owned_top_paths;
                const auto collect_owned_top_paths =
                    [&](auto &&self, const ExtrusionEntity &entity,
                        double owner_z) -> void {
                        if (const auto *path =
                                dynamic_cast<const ExtrusionPath *>(&entity)) {
                            if (path->nonplanar_surface &&
                                path->role() == erTopSolidInfill)
                                owned_top_paths.push_back({path, owner_z});
                        } else if (const auto *multipath =
                                       dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
                            for (const ExtrusionPath &path : multipath->paths)
                                self(self, path, owner_z);
                        } else if (const auto *loop =
                                       dynamic_cast<const ExtrusionLoop *>(&entity)) {
                            for (const ExtrusionPath &path : loop->paths)
                                self(self, path, owner_z);
                        } else if (const auto *collection =
                                       dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
                            for (const ExtrusionEntity *child : collection->entities)
                                self(self, *child, owner_z);
                        }
                    };
                for (size_t owner = 0; owner < layers.size(); ++owner)
                    collect_owned_top_paths(
                        collect_owned_top_paths,
                        layers[owner]->regions()[region_idx]->fills,
                        layers[owner]->print_z);
                const auto separate_wall_from_top_skin =
                    [&owned_top_paths](ExtrusionPath &wall,
                                       double owner_z) {
                        bool changed = false;
                        for (size_t wall_idx = 1;
                             wall_idx < wall.polyline.points.size(); ++wall_idx) {
                            Point3 &wall_first = wall.polyline.points[wall_idx - 1];
                            Point3 &wall_last = wall.polyline.points[wall_idx];
                            const Vec2d wall_start =
                                unscale(wall_first.to_point()).cast<double>();
                            const Vec2d wall_delta =
                                unscale(wall_last.to_point()).cast<double>() -
                                wall_start;
                            for (const OwnedTopPath &top : owned_top_paths)
                                for (size_t top_idx = 1;
                                     top_idx < top.path->polyline.points.size();
                                     ++top_idx) {
                                    const Point3 &top_first =
                                        top.path->polyline.points[top_idx - 1];
                                    const Point3 &top_last =
                                        top.path->polyline.points[top_idx];
                                    const Vec2d top_start =
                                        unscale(top_first.to_point()).cast<double>();
                                    const Vec2d top_delta =
                                        unscale(top_last.to_point()).cast<double>() -
                                        top_start;
                                    const double denominator =
                                        wall_delta.x() * top_delta.y() -
                                        wall_delta.y() * top_delta.x();
                                    if (std::abs(denominator) <= EPSILON)
                                        continue;
                                    const Vec2d offset = top_start - wall_start;
                                    const double wall_t =
                                        (offset.x() * top_delta.y() -
                                         offset.y() * top_delta.x()) /
                                        denominator;
                                    const double top_t =
                                        (offset.x() * wall_delta.y() -
                                         offset.y() * wall_delta.x()) /
                                        denominator;
                                    if (wall_t <= 0. || wall_t >= 1. ||
                                        top_t <= 0. || top_t >= 1.)
                                        continue;
                                    const double wall_z = owner_z +
                                        unscale_(wall_first.z()) + wall_t *
                                            unscale_(wall_last.z() - wall_first.z());
                                    const double top_z = top.owner_z +
                                        unscale_(top_first.z()) + top_t *
                                            unscale_(top_last.z() - top_first.z());
                                    constexpr double skin_clearance = 0.051;
                                    if (top_z <= wall_z + EPSILON ||
                                        top_z - wall_z >= skin_clearance)
                                        continue;
                                    const coord_t correction = scale_(
                                        skin_clearance - (top_z - wall_z));
                                    wall_first.z() -= correction;
                                    wall_last.z() -= correction;
                                    changed = true;
                                }
                        }
                        if (changed)
                            wall.polyline.fitting_result.clear();
                        return changed;
                    };
                for (DeferredCourse &course : courses) {
                    std::vector<ExtrusionPath *> wall_paths;
                    for (ExtrusionEntity *entity : course.perimeters)
                        collect_course_paths(
                            collect_course_paths, *entity, wall_paths);
                    const double owner_z =
                        layers[course.destination_layer]->print_z;
                    for (ExtrusionPath *wall : wall_paths) {
                        bool changed = false;
                        for (Point3 &point : wall->polyline.points) {
                            const Vec2d sample =
                                unscale(point.to_point()).cast<double>();
                            const double wall_z = owner_z + unscale_(point.z());
                            double nearest_skin_z =
                                std::numeric_limits<double>::infinity();
                            for (const OwnedTopPath &top : owned_top_paths) {
                                for (size_t point_idx = 1;
                                     point_idx < top.path->polyline.points.size();
                                     ++point_idx) {
                                    const Point3 &first =
                                        top.path->polyline.points[point_idx - 1];
                                    const Point3 &last =
                                        top.path->polyline.points[point_idx];
                                    const Vec2d a =
                                        unscale(first.to_point()).cast<double>();
                                    const Vec2d b =
                                        unscale(last.to_point()).cast<double>();
                                    const Vec2d delta = b - a;
                                    const double squared = delta.squaredNorm();
                                    const double ratio = squared <= EPSILON ? 0. :
                                        std::clamp((sample - a).dot(delta) /
                                                       squared,
                                                   0., 1.);
                                    if ((sample - (a + ratio * delta)).norm() >
                                        0.5 * (double(wall->width) +
                                               double(top.path->width)) + 0.02)
                                        continue;
                                    const double skin_z = top.owner_z +
                                        unscale_(first.z()) + ratio *
                                            unscale_(last.z() - first.z());
                                    if (skin_z > wall_z + EPSILON)
                                        nearest_skin_z = std::min(
                                            nearest_skin_z, skin_z);
                                }
                            }
                            if (!std::isfinite(nearest_skin_z) ||
                                wall_z <= nearest_skin_z - 0.0501 + EPSILON)
                                continue;
                            point.z() = scale_(
                                nearest_skin_z - 0.0501 - owner_z);
                            changed = true;
                        }
                        for (size_t wall_idx = 1;
                             wall_idx < wall->polyline.points.size();
                             ++wall_idx) {
                            Point3 &wall_first =
                                wall->polyline.points[wall_idx - 1];
                            Point3 &wall_last =
                                wall->polyline.points[wall_idx];
                            const Vec2d wall_start =
                                unscale(wall_first.to_point()).cast<double>();
                            const Vec2d wall_delta =
                                unscale(wall_last.to_point()).cast<double>() -
                                wall_start;
                            for (const OwnedTopPath &top : owned_top_paths) {
                                for (size_t top_idx = 1;
                                     top_idx < top.path->polyline.points.size();
                                     ++top_idx) {
                                    const Point3 &top_first =
                                        top.path->polyline.points[top_idx - 1];
                                    const Point3 &top_last =
                                        top.path->polyline.points[top_idx];
                                    const Vec2d top_start =
                                        unscale(top_first.to_point()).cast<double>();
                                    const Vec2d top_delta =
                                        unscale(top_last.to_point()).cast<double>() -
                                        top_start;
                                    const double denominator =
                                        wall_delta.x() * top_delta.y() -
                                        wall_delta.y() * top_delta.x();
                                    if (std::abs(denominator) <= EPSILON)
                                        continue;
                                    const Vec2d offset = top_start - wall_start;
                                    const double wall_t =
                                        (offset.x() * top_delta.y() -
                                         offset.y() * top_delta.x()) /
                                        denominator;
                                    const double top_t =
                                        (offset.x() * wall_delta.y() -
                                         offset.y() * wall_delta.x()) /
                                        denominator;
                                    if (wall_t <= 0. || wall_t >= 1. ||
                                        top_t <= 0. || top_t >= 1.)
                                        continue;
                                    const double wall_z = owner_z +
                                        unscale_(wall_first.z()) + wall_t *
                                            unscale_(wall_last.z() -
                                                     wall_first.z());
                                    const double top_z = top.owner_z +
                                        unscale_(top_first.z()) + top_t *
                                            unscale_(top_last.z() - top_first.z());
                                    constexpr double skin_clearance = 0.051;
                                    if (top_z <= wall_z + EPSILON ||
                                        top_z - wall_z >= skin_clearance)
                                        continue;
                                    const coord_t correction = scale_(
                                        skin_clearance - (top_z - wall_z));
                                    wall_first.z() -= correction;
                                    wall_last.z() -= correction;
                                    changed = true;
                                }
                            }
                        }
                        if (changed)
                            wall->polyline.fitting_result.clear();
                        // Keep the complete feature wall just inside the
                        // finishing skin. G-code simplification may remove
                        // micron-scale local correction vertices; a uniform
                        // 0.01 mm bonding inset survives simplification while
                        // remaining far below geometric/extrusion tolerance.
                        for (Point3 &point : wall->polyline.points)
                            point.z() -= scale_(0.01);
                    }
                }

                // Boundary-contact corrections above may flatten a locally
                // horizontal bore or boss ring after its initial projection.
                // Normalize metadata only after those corrections so preview
                // and scheduling never advertise a fixed-Z closed loop as a
                // non-planar transition.
                for (DeferredCourse &course : courses) {
                    for (ExtrusionEntity *entity : course.perimeters)
                        normalize_fixed_feature_transition(*entity);
                    for (ExtrusionEntity *entity : course.fills)
                        normalize_fixed_feature_transition(*entity);
                }

                // Validate the transformed courses against the same generic
                // gantry model used for the skin. The local patch is excluded
                // from the obstacle set because it supports these courses;
                // unrelated walls and features remain collision obstacles.
                ExPolygons local_transition_area = patch.projection;
                // The cross-section of a continuing boss/chimney changes
                // slightly through the deferred height range. Excluding only
                // its first post-surface slice leaves a thin part of an older
                // ring in the nominal-layer obstacle map, even though that
                // ring is the local support for the transformed course. Grow
                // only the feature footprint by one bead radius; remote walls
                // and every part of the model beyond that local allowance
                // remain collision obstacles.
                append(local_transition_area, offset_ex(
                    accepted.feature_footprint,
                    scale_(0.5 * nozzle_diameter + 0.05)));
                local_transition_area = union_ex(local_transition_area);
                // Every deferred course sees the same immutable model with a
                // different last-printed-layer prefix. Building that model per
                // course repeated all layer/transition Clipper differences and
                // dominated feature clearance on real meshes. Construct it
                // once, then cache pointer prefixes exactly as the finishing
                // shell validator above does.
                const std::vector<PrintedToolpathObstacle> feature_obstacles =
                    printed_model_obstacles(
                        layers, layers.size() - 1, &local_transition_area);
                const auto feature_obstacle_prefix =
                    [&feature_obstacles, &layers](size_t last_layer) {
                        std::vector<const PrintedToolpathObstacle *> prefix;
                        const double maximum_z =
                            layers[last_layer]->print_z + EPSILON;
                        prefix.reserve(feature_obstacles.size());
                        for (const PrintedToolpathObstacle &obstacle :
                             feature_obstacles) {
                            if (obstacle.top_z > maximum_z)
                                break;
                            prefix.emplace_back(&obstacle);
                        }
                        return prefix;
                    };
                std::vector<std::optional<std::vector<
                    const PrintedToolpathObstacle *>>> feature_obstacle_cache(
                        layers.size());
                bool clearance_valid = projection_succeeded;
                for (size_t course_idx = 0; course_idx < courses.size(); ++course_idx) {
                    DeferredCourse &course = courses[course_idx];
                    report_progress(
                        double(run_end), NonplanarProgressStage::FeatureClearance,
                        course_idx + 1, courses.size());
                    if (!feature_obstacle_cache[course.destination_layer]) {
                        const size_t last_layer =
                            accepted.before_current_layer &&
                                    course.destination_layer > 0 ?
                                course.destination_layer - 1 :
                                course.destination_layer;
                        feature_obstacle_cache[course.destination_layer] =
                            feature_obstacle_prefix(last_layer);
                    }
                    const auto &obstacles =
                        *feature_obstacle_cache[course.destination_layer];
                    for (ExtrusionEntity *entity : course.perimeters)
                        clearance_valid &= validate_deferred_feature(
                            *entity, layers[course.destination_layer]->print_z,
                            gantry, fallback, obstacles);
                    for (ExtrusionEntity *entity : course.fills)
                        clearance_valid &= validate_deferred_feature(
                            *entity, layers[course.destination_layer]->print_z,
                            gantry, fallback, obstacles);
                }

                if (clearance_valid) {
                    // At an attached boss/chimney the accepted roof already
                    // owns the outer opening wall for each sloped surface
                    // course. Generating another deferred inset-0 ring makes
                    // two closed walls cross at the same XYZ. Adopt that
                    // existing surface wall as the feature's outer wall and
                    // retain only the additional feature insets. This is the
                    // graph-level splice: it preserves the supporting roof
                    // course instead of clipping it away, while avoiding a
                    // disconnected duplicate loop.
                    for (size_t course_idx = 0;
                         course_idx < surface_course_count; ++course_idx) {
                        DeferredCourse &course = courses[course_idx];
                        ExtrusionEntityCollection &surface_walls =
                            layers[course.destination_layer]
                                ->regions()[region_idx]
                                ->perimeters;
                        const double owner_z =
                            layers[course.destination_layer]->print_z;
                        const auto same_feature_boundary =
                            [owner_z](const ExtrusionPath &surface,
                                      const ExtrusionPath &candidate) {
                                if (paths_have_proper_same_z_crossing(
                                        surface, owner_z, candidate, owner_z,
                                        0.05))
                                    return true;
                                if (surface.polyline.points.size() < 2 ||
                                    candidate.polyline.points.size() < 2)
                                    return false;
                                const double contact_distance =
                                    0.5 * (double(surface.width) +
                                           double(candidate.width)) + 0.05;
                                const double contact_distance_squared =
                                    contact_distance * contact_distance;
                                size_t matching_points = 0;
                                for (const Point3 &point :
                                     candidate.polyline.points) {
                                    const Vec2d sample =
                                        unscale(point.to_point()).cast<double>();
                                    bool matches = false;
                                    for (size_t idx = 1;
                                         idx < surface.polyline.points.size();
                                         ++idx) {
                                        const Point3 &first =
                                            surface.polyline.points[idx - 1];
                                        const Point3 &last =
                                            surface.polyline.points[idx];
                                        const Vec2d a = unscale(
                                            first.to_point()).cast<double>();
                                        const Vec2d b = unscale(
                                            last.to_point()).cast<double>();
                                        const Vec2d delta = b - a;
                                        const double squared =
                                            delta.squaredNorm();
                                        const double ratio =
                                            squared <= EPSILON ? 0. :
                                            std::clamp(
                                                (sample - a).dot(delta) /
                                                    squared,
                                                0., 1.);
                                        if ((sample - (a + ratio * delta))
                                                .squaredNorm() >
                                            contact_distance_squared)
                                            continue;
                                        const double surface_z = owner_z +
                                            unscale_(first.z()) + ratio *
                                                unscale_(last.z() - first.z());
                                        const double candidate_z = owner_z +
                                            unscale_(point.z());
                                        if (std::abs(surface_z - candidate_z) <=
                                            0.08) {
                                            matches = true;
                                            break;
                                        }
                                    }
                                    matching_points += matches;
                                }
                                // Endpoint contact is not ownership. Adopt
                                // only when a meaningful fraction of the
                                // deferred contour follows the already
                                // generated surface boundary.
                                return matching_points >= 2 &&
                                    matching_points * 2 >=
                                        candidate.polyline.points.size();
                            };
                        const auto adopt_surface_wall =
                            [&](auto &&self, ExtrusionEntity &entity,
                                const ExtrusionPath &candidate) -> bool {
                                if (auto *path =
                                        dynamic_cast<ExtrusionPath *>(&entity)) {
                                    if (!path->nonplanar_surface ||
                                        path->nonplanar_feature_transition ||
                                        !is_perimeter(path->role()) ||
                                        path->inset_idx != candidate.inset_idx ||
                                        !same_feature_boundary(
                                            *path, candidate))
                                        return false;
                                    path->nonplanar_feature_transition = true;
                                    path->nonplanar_feature_course =
                                        candidate.nonplanar_feature_course;
                                    separate_wall_from_top_skin(*path, owner_z);
                                    for (Point3 &point : path->polyline.points)
                                        point.z() -= scale_(0.01);
                                    path->polyline.fitting_result.clear();
                                    normalize_fixed_feature_transition(*path);
                                    return true;
                                }
                                if (auto *multipath =
                                        dynamic_cast<ExtrusionMultiPath *>(&entity))
                                    for (ExtrusionPath &path : multipath->paths)
                                        if (self(self, path, candidate))
                                            return true;
                                if (auto *loop =
                                        dynamic_cast<ExtrusionLoop *>(&entity))
                                    for (ExtrusionPath &path : loop->paths)
                                        if (self(self, path, candidate))
                                            return true;
                                if (auto *collection =
                                        dynamic_cast<ExtrusionEntityCollection *>(&entity))
                                    for (ExtrusionEntity *child : collection->entities)
                                        if (self(self, *child, candidate))
                                            return true;
                                return false;
                            };
                        const auto merge_outer =
                            [&](auto &&self, ExtrusionEntity &entity) -> void {
                                if (auto *path =
                                        dynamic_cast<ExtrusionPath *>(&entity)) {
                                    if (path->inset_idx != 0 ||
                                        !path->nonplanar_feature_transition)
                                        return;
                                    for (ExtrusionEntity *surface_entity :
                                         surface_walls.entities)
                                        if (adopt_surface_wall(
                                                adopt_surface_wall,
                                                *surface_entity, *path)) {
                                            path->polyline.points.clear();
                                            return;
                                        }
                                } else if (auto *multipath =
                                               dynamic_cast<ExtrusionMultiPath *>(&entity)) {
                                    for (ExtrusionPath &path : multipath->paths)
                                        self(self, path);
                                } else if (auto *loop =
                                               dynamic_cast<ExtrusionLoop *>(&entity)) {
                                    for (ExtrusionPath &path : loop->paths)
                                        self(self, path);
                                } else if (auto *collection =
                                               dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
                                    for (ExtrusionEntity *child : collection->entities)
                                        self(self, *child);
                                }
                            };
                        for (ExtrusionEntity *entity : course.perimeters)
                            merge_outer(merge_outer, *entity);
                        for (auto it = course.perimeters.begin();
                             it != course.perimeters.end();) {
                            if (prune_empty_extrusion_paths(**it)) {
                                ++it;
                            } else {
                                delete *it;
                                it = course.perimeters.erase(it);
                            }
                        }
                    }

                    // A feature wall may meet a generated raster between two
                    // nominal source planes, so the physical-Z cleanup below
                    // cannot see that handoff. In the wall's destination
                    // owner, terminate accepted top paths at the validated
                    // wall bead. Perimeters and lower support courses are not
                    // touched here.
                    // These masks describe conventional fill at its owner Z,
                    // not the varying height of an accepted surface chain.
                    // Cutting that chain by XY alone removes unrelated height
                    // samples and leaves isolated variable-flow fragments.
                    // Preserve accepted courses; feature-to-skin clearance
                    // must be resolved by the preceding 3D validation.
                    for (const DeferredCourse &course : courses) {
                        Polygons course_wall_material;
                        for (const ExtrusionEntity *entity : course.perimeters)
                            collect_deposited_coverage(
                                *entity, course_wall_material);
                        if (course_wall_material.empty())
                            continue;
                        const ExPolygons course_wall_coverage = intersection_ex(
                            offset_ex(
                                union_ex(course_wall_material), scale_(0.02)),
                            offset_ex(
                                accepted.feature_footprint, scale_(0.05)));
                        if (!course_wall_coverage.empty())
                            remove_nonplanar_top_overlap(
                                layers[course.destination_layer]
                                    ->regions()[region_idx]
                                    ->fills,
                                course_wall_coverage);
                    }

                    // Dependency deferral happens after the main replacement
                    // shell has already clipped its conventional height band.
                    // A transformed feature wall may nevertheless descend
                    // through fill owned by an earlier nominal layer. Carve
                    // only the centerline fragments covered by that wall's
                    // real bead envelope at the fill's physical Z. This is
                    // deliberately delayed until every transformed course
                    // passes projection and gantry clearance so rollback
                    // remains atomic.
                    constexpr double feature_fill_z_tolerance = 0.05;
                    for (size_t source_idx = first_source;
                         source_idx <= accepted.surface_owner; ++source_idx) {
                        Polygons wall_material;
                        const double source_z =
                            layers[source_idx]->print_z;
                        for (const DeferredCourse &course : courses) {
                            const double owner_z =
                                layers[course.destination_layer]->print_z;
                            for (const ExtrusionEntity *entity :
                                 course.perimeters)
                                collect_generated_material_at_z(
                                    *entity, owner_z, source_z,
                                    feature_fill_z_tolerance, false,
                                    wall_material);
                        }
                        // Clipper subtraction treats a centerline exactly on
                        // a bead-envelope boundary as outside. With curved
                        // Arachne rings this left tiny conventional perimeter
                        // fragments which crossed the replacement feature
                        // wall at the same physical Z. Apply the same small
                        // endpoint/join tolerance used by topology validation,
                        // but only to this accepted wall handoff; unrelated
                        // model material and global wall/fill spacing remain
                        // unchanged.
                        const ExPolygons wall_coverage = offset_ex(
                            union_ex(wall_material), scale_(0.05));
                        if (!wall_coverage.empty()) {
                            remove_nonplanar_top_overlap(
                                layers[source_idx]
                                    ->regions()[region_idx]
                                    ->fills,
                                wall_coverage);
                            // The initial marking/roof pass may leave a
                            // conventional stair-wall fragment in the same
                            // physical band even when the complete feature
                            // loop was extracted successfully. Once every
                            // transformed course has passed projection and
                            // clearance, remove only the centerline pieces
                            // covered by its real bead envelope. This avoids
                            // emitting a second planar chimney ring beside or
                            // through the deferred wall while preserving all
                            // perimeter material outside the local handoff.
                        }
                        // A raised brick wall does not occupy source_z. Query
                        // replacement ownership at its actual course height,
                        // separately from nominal walls/fill; using the nominal
                        // band left duplicate chimney arcs at a 50% offset.
                        std::map<coord_t, ExPolygons> coverage_by_offset{{0, wall_coverage}};
                        const std::function<ExPolygons(const ExtrusionPath &)> coverage_at_path =
                            [&](const ExtrusionPath &path) {
                                coord_t offset = 0;
                                if (path.z_contoured && !path.polyline.points.empty() &&
                                    std::all_of(path.polyline.points.begin(), path.polyline.points.end(),
                                        [&](const Point3 &point) {
                                            return point.z() == path.polyline.points.front().z();
                                        }))
                                    offset = path.polyline.points.front().z();
                                auto found = coverage_by_offset.find(offset);
                                if (found == coverage_by_offset.end()) {
                                    Polygons physical_material;
                                    for (const DeferredCourse &course : courses)
                                        for (const ExtrusionEntity *entity : course.perimeters)
                                            collect_generated_material_at_z(*entity,
                                                layers[course.destination_layer]->print_z,
                                                source_z + unscale_(offset),
                                                feature_fill_z_tolerance, false, physical_material);
                                    found = coverage_by_offset.emplace(offset,
                                        offset_ex(union_ex(physical_material), scale_(0.05))).first;
                                }
                                return found->second;
                            };
                        remove_nonplanar_top_overlap(
                            layers[source_idx]->regions()[region_idx]->perimeters,
                            wall_coverage, false, true, &coverage_at_path);
                    }
                }

                if (clearance_valid) {
                    for (DeferredCourse &course : superseded_courses) {
                        destroy_entities(course.perimeters);
                        destroy_entities(course.fills);
                    }
                } else {
                    // Roll the dependency back atomically. None of the
                    // conventional feature courses may be lost when either
                    // projection or the generic gantry check rejects the
                    // transformed schedule.
                    for (DeferredCourse &course : superseded_courses) {
                        LayerRegion &source =
                            *layers[course.source_layer]->regions()[region_idx];
                        append_ordered_entities(
                            source.perimeters, std::move(course.perimeters), 0);
                        append_ordered_entities(
                            source.fills, std::move(course.fills), 0);
                    }
                }

                for (DeferredCourse &course : courses) {
                    if (!clearance_valid && course.synthetic) {
                        destroy_entities(course.perimeters);
                        destroy_entities(course.fills);
                        continue;
                    }
                    const size_t owner = clearance_valid ?
                        course.destination_layer : course.source_layer;
                    if (!clearance_valid) {
                        for (ExtrusionEntity *entity : course.perimeters)
                            restore_deferred_feature(*entity);
                        for (ExtrusionEntity *entity : course.fills)
                            restore_deferred_feature(*entity);
                    } else {
                        // Boundary adoption and overlap cleanup above are the
                        // last operations allowed to alter Z. They may flatten
                        // a locally horizontal bore/boss ring after the first
                        // normalization pass. Normalize again at commit time
                        // so fixed-Z feature walls retain dependency ordering
                        // without being mislabeled as non-planar transitions.
                        for (ExtrusionEntity *entity : course.perimeters)
                            normalize_fixed_feature_transition(*entity);
                        for (ExtrusionEntity *entity : course.fills)
                            normalize_fixed_feature_transition(*entity);
                    }
                    LayerRegion &region = *layers[owner]->regions()[region_idx];
                    append_ordered_entities(
                        region.perimeters, std::move(course.perimeters), 0);
                    append_ordered_entities(
                        region.fills, std::move(course.fills), 0);
                    if (clearance_valid) {
                        // An adopted boundary lives in the destination
                        // region, not necessarily in course.perimeters.
                        // Normalize the committed graph as well so every
                        // fixed-Z feature ring is classified consistently
                        // regardless of which side owned the splice.
                        for (ExtrusionEntity *entity : region.perimeters.entities)
                            normalize_fixed_feature_transition(*entity);
                        for (ExtrusionEntity *entity : region.fills.entities)
                            normalize_fixed_feature_transition(*entity);
                        last_operation_layer = std::max(
                            last_operation_layer, course.destination_layer);
                    }
                }
                BOOST_LOG_TRIVIAL(debug) << "Nonplanar patch "
                                         << accepted.patch_index << " deferred "
                                         << warped_count
                                         << " continuing-feature paths across "
                                         << surface_course_count
                                         << " surface courses and "
                                         << leveling_course_count
                                         << " leveling courses"
                                         << (clearance_valid ? "" :
                                             projection_succeeded ?
                                                " (collision check rejected deferral)" :
                                             " (surface projection rejected deferral)");
            }
            for (auto &patch_snapshots : feature_course_snapshots)
                for (ExtrusionEntitiesPtr &snapshot : patch_snapshots)
                    destroy_entities(snapshot);

            // Feature deferral appends its transformed wall courses after the
            // replacement shell is committed. Restore the configured wall
            // sequence once more over the final graph so those late courses
            // cannot fall back to their extraction (inside-out) order.
            for (size_t owner = run_begin;
                 owner <= last_operation_layer; ++owner) {
                LayerRegion &region = *layers[owner]->regions()[region_idx];
                WallSequence wall_sequence = region_config.wall_sequence.value;
                if (region_config.wall_loops.value >= 3 &&
                    region_config.inner_walls_flow_ratio.get_abs_value(1.) >
                        1.0 + EPSILON && owner > 0)
                    wall_sequence = WallSequence::InnerOuterInner;
                restore_nonplanar_wall_sequence(
                    region.perimeters, wall_sequence,
                    region_config.wall_loops.value, owner > 0);
            }

            layer_idx = last_operation_layer + 1;
            report_progress(
                double(layer_idx), NonplanarProgressStage::Commit,
                std::min(layer_idx, layers.size()), layers.size());
        }
    }
    // Continuing-feature deferral is committed after each surface patch, so
    // two paths which originate in different patches never coexist in either
    // patch's local clearance graph. Validate the completed perimeter graph
    // at absolute Z before restoring final wall order. Keep the earlier
    // course as support; for two fragments on the same owner, keep the larger
    // wall piece. A previous fallback merely relabeled the rejected path as
    // planar, leaving its identical geometry in the final graph and therefore
    // converting one detected collision into an ordinary-vs-nonplanar
    // collision. Mark the redundant extrusion empty and prune it atomically
    // after pointer-based validation instead.
    struct FinalNonplanarWall {
        ExtrusionPath *path;
        size_t owner;
        BoundingBox bounds;
        double minimum_z;
        double maximum_z;
        ExtrusionEntity *continuous_wall = nullptr;
    };
    for (size_t region_idx = 0; region_idx < region_count; ++region_idx) {
        std::vector<FinalNonplanarWall> walls;
        std::vector<FinalNonplanarWall> conventional_walls;
        std::vector<FinalNonplanarWall> top_paths;
        size_t final_crossing_count = 0;
        const auto collect_walls = [&](auto &&self, ExtrusionEntity &entity,
                                       size_t owner,
                                       bool scheduled_nonplanar,
                                       ExtrusionEntity *continuous_wall) -> void {
            if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
                if (!is_perimeter(path->role()) ||
                    path->polyline.points.size() < 2)
                    return;
                // Native Arachne paths may keep their inset only on the
                // parent. All final support/collision comparisons below use
                // the child field, including comparisons with unwarped walls.
                if (continuous_wall != nullptr && continuous_wall->inset_idx >= 0)
                    path->inset_idx = continuous_wall->inset_idx;
                double minimum_z = std::numeric_limits<double>::max();
                double maximum_z = std::numeric_limits<double>::lowest();
                for (const Point3 &point : path->polyline.points) {
                    const double z = layers[owner]->print_z + unscale_(point.z());
                    minimum_z = std::min(minimum_z, z);
                    maximum_z = std::max(maximum_z, z);
                }
                // Recheck this invariant on the final graph.  A candidate may
                // have had a real Z span when it was selected, then collapse
                // onto one owner height while blending from planar support or
                // separating dependency courses.  Keeping the provisional
                // non-planar flag in that state schedules an otherwise
                // ordinary closed wall after the layer, with a retract and a
                // clearance travel.  That is the detached fixed-Z loop seen
                // around the Benchy stern.  A perimeter with no emitted Z
                // variation belongs to its normal owner-layer wall graph.
                if (path->nonplanar_surface &&
                    maximum_z - minimum_z < 0.001)
                    normalize_fixed_feature_transition(*path);
                const BoundingBox bounds =
                    path->polyline.to_polyline().bounding_box();
                if (path->nonplanar_surface || scheduled_nonplanar)
                    walls.push_back({path, owner, bounds, minimum_z, maximum_z,
                                     continuous_wall ? continuous_wall : path});
                if (!path->nonplanar_surface)
                    conventional_walls.push_back(
                        {path, owner, bounds, minimum_z, maximum_z});
            } else if (auto *multipath =
                           dynamic_cast<ExtrusionMultiPath *>(&entity)) {
                const bool emitted_nonplanar = scheduled_nonplanar ||
                    std::any_of(
                        multipath->paths.begin(), multipath->paths.end(),
                        [](const ExtrusionPath &path) {
                            return path.nonplanar_surface;
                        });
                for (ExtrusionPath &path : multipath->paths)
                    self(self, path, owner, emitted_nonplanar, &entity);
            } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
                const bool emitted_nonplanar = scheduled_nonplanar ||
                    std::any_of(
                        loop->paths.begin(), loop->paths.end(),
                        [](const ExtrusionPath &path) {
                            return path.nonplanar_surface;
                        });
                for (ExtrusionPath &path : loop->paths)
                    self(self, path, owner, emitted_nonplanar, &entity);
            } else if (auto *collection =
                           dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
                for (ExtrusionEntity *child : collection->entities)
                    self(self, *child, owner, scheduled_nonplanar, nullptr);
            }
        };
        for (size_t owner = 0; owner < layers.size(); ++owner) {
            collect_walls(
                collect_walls,
                layers[owner]->regions()[region_idx]->perimeters, owner, false, nullptr);
            // Arachne outlines generated with the finishing skin are stored
            // in fills while retaining perimeter roles. They participate in
            // the same physical wall graph and must be visible to duplicate
            // and crossing arbitration.
            collect_walls(
                collect_walls,
                layers[owner]->regions()[region_idx]->fills, owner, false, nullptr);
        }
        // Arachne changes width by splitting one continuous wall into paths.
        // Clearance translations must move that wall atomically: lifting just
        // one child creates a vertical travel at its formerly shared endpoint.
        const auto shift_wall = [&](FinalNonplanarWall &wall, double distance) {
            const coord_t delta = coord_t(std::ceil(scale_(distance)));
            offset_entity_z(*wall.continuous_wall, delta);
            normalize_fixed_feature_transition(*wall.continuous_wall);
            for (FinalNonplanarWall &part : walls)
                if (part.continuous_wall == wall.continuous_wall) {
                    part.path->polyline.fitting_result.clear();
                    part.minimum_z += unscale_(delta);
                    part.maximum_z += unscale_(delta);
                }
        };
        const auto collect_top_paths = [&](auto &&self,
                                           ExtrusionEntity &entity,
                                           size_t owner) -> void {
            if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
                if (!path->nonplanar_surface ||
                    path->role() != erTopSolidInfill ||
                    path->polyline.points.size() < 2)
                    return;
                double minimum_z = std::numeric_limits<double>::max();
                double maximum_z = std::numeric_limits<double>::lowest();
                for (const Point3 &point : path->polyline.points) {
                    const double z = layers[owner]->print_z +
                                     unscale_(point.z());
                    minimum_z = std::min(minimum_z, z);
                    maximum_z = std::max(maximum_z, z);
                }
                top_paths.push_back({
                    path, owner, path->polyline.to_polyline().bounding_box(),
                    minimum_z, maximum_z});
            } else if (auto *multipath =
                           dynamic_cast<ExtrusionMultiPath *>(&entity)) {
                for (ExtrusionPath &path : multipath->paths)
                    self(self, path, owner);
            } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
                for (ExtrusionPath &path : loop->paths)
                    self(self, path, owner);
            } else if (auto *collection =
                           dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
                for (ExtrusionEntity *child : collection->entities)
                    self(self, *child, owner);
            }
        };
        for (size_t owner = 0; owner < layers.size(); ++owner)
            collect_top_paths(
                collect_top_paths,
                layers[owner]->regions()[region_idx]->fills, owner);

        // Skin-clearance correction deliberately moves feature walls only
        // downward.  Re-establish course ordering afterwards, including an
        // outer wall adopted from the surface graph (which is not present in
        // the provisional DeferredCourse array).  Keep the later course on
        // its mesh-derived height and lower the preceding course where their
        // bead neighborhoods overlap.  This preserves the exposed geometry,
        // improves support for the lowered course, and prevents two locally
        // shifted chimney/boss contours from crossing at the same physical Z.
        // Keep enough centerline separation for the two deposited beads, not
        // merely enough numerical distance to make the crossing predicate
        // pass.  The old fixed 0.0601 mm clearance let nominal 0.25 mm
        // chimney courses become a tightly packed stack which the nozzle
        // would drag through.  Derive the limit from the actual pair of path
        // heights so fine and coarse profiles use the same physical rule.
        constexpr double minimum_feature_course_clearance = 0.0601;
        const auto enforce_feature_course_clearance = [&]() {
            std::vector<size_t> feature_wall_order(walls.size());
            std::iota(
                feature_wall_order.begin(), feature_wall_order.end(), size_t(0));
            std::stable_sort(
                feature_wall_order.begin(), feature_wall_order.end(),
                [&walls](size_t left, size_t right) {
                    return walls[left].path->nonplanar_feature_course >
                           walls[right].path->nonplanar_feature_course;
                });
            for (const size_t lower_idx : feature_wall_order) {
            FinalNonplanarWall &lower = walls[lower_idx];
            if (!lower.path->nonplanar_feature_transition ||
                lower.path->nonplanar_feature_course == 0)
                continue;
            size_t next_course = std::numeric_limits<size_t>::max();
            for (const FinalNonplanarWall &candidate : walls)
                if (candidate.path->nonplanar_feature_transition &&
                    candidate.path->inset_idx == lower.path->inset_idx &&
                    candidate.path->nonplanar_feature_course >
                        lower.path->nonplanar_feature_course)
                    next_course = std::min(
                        next_course,
                        candidate.path->nonplanar_feature_course);
            if (next_course == std::numeric_limits<size_t>::max())
                continue;
            const double lower_owner_z = layers[lower.owner]->print_z;
            bool changed = false;
            for (Point3 &point : lower.path->polyline.points) {
                const Vec2d sample = unscale(point.to_point()).cast<double>();
                double limiting_z = std::numeric_limits<double>::infinity();
                for (const FinalNonplanarWall &higher : walls) {
                    if (!higher.path->nonplanar_feature_transition ||
                        higher.path->nonplanar_feature_course != next_course ||
                        higher.path->inset_idx != lower.path->inset_idx)
                        continue;
                    const double contact_distance =
                        0.5 * (double(lower.path->width) +
                               double(higher.path->width)) + 0.05;
                    const double contact_distance_squared =
                        contact_distance * contact_distance;
                    const double higher_owner_z = layers[higher.owner]->print_z;
                    const Points3 &higher_points = higher.path->polyline.points;
                    for (size_t point_idx = 1;
                         point_idx < higher_points.size(); ++point_idx) {
                        const Point3 &first = higher_points[point_idx - 1];
                        const Point3 &last = higher_points[point_idx];
                        const Vec2d a = unscale(first.to_point()).cast<double>();
                        const Vec2d b = unscale(last.to_point()).cast<double>();
                        const Vec2d delta = b - a;
                        const double squared = delta.squaredNorm();
                        const double ratio = squared <= EPSILON ? 0. :
                            std::clamp((sample - a).dot(delta) / squared,
                                       0., 1.);
                        if ((sample - (a + ratio * delta)).squaredNorm() >
                            contact_distance_squared)
                            continue;
                        const double required_course_clearance = std::max(
                            minimum_feature_course_clearance,
                            0.5 * std::max(double(lower.path->height),
                                           double(higher.path->height)) +
                                0.001);
                        limiting_z = std::min(
                            limiting_z,
                            higher_owner_z + unscale_(first.z()) + ratio *
                                unscale_(last.z() - first.z()) -
                                required_course_clearance);
                    }
                }
                const double current_z = lower_owner_z + unscale_(point.z());
                if (std::isfinite(limiting_z) &&
                    current_z > limiting_z + EPSILON) {
                    point.z() = scale_(limiting_z - lower_owner_z);
                    changed = true;
                }
            }
            if (changed) {
                lower.path->polyline.fitting_result.clear();
                lower.minimum_z = std::numeric_limits<double>::max();
                lower.maximum_z = std::numeric_limits<double>::lowest();
                for (const Point3 &point : lower.path->polyline.points) {
                    const double z = lower_owner_z + unscale_(point.z());
                    lower.minimum_z = std::min(lower.minimum_z, z);
                    lower.maximum_z = std::max(lower.maximum_z, z);
                }
            }
            }
        };
        enforce_feature_course_clearance();

        // The roof boundary immediately below an attached feature is course
        // zero of the same physical stack, even though it intentionally does
        // not carry feature-course metadata until full-loop adoption is
        // possible.  Partial boundaries cannot replace a complete chimney
        // ring, so retain both and order the roof segment just below the first
        // nearby deferred feature wall instead of letting near-coincident
        // centerlines survive until G-code rounding.
#if 0
        for (FinalNonplanarWall &surface : walls) {
            if (!surface.path->nonplanar_surface ||
                surface.path->nonplanar_transition ||
                surface.path->nonplanar_feature_transition)
                continue;
            const double owner_z = layers[surface.owner]->print_z;
            bool changed = false;
            for (Point3 &point : surface.path->polyline.points) {
                const Vec2d sample = unscale(point.to_point()).cast<double>();
                const double current_z = owner_z + unscale_(point.z());
                double nearest_upper_z = std::numeric_limits<double>::infinity();
                double nearest_z_distance =
                    std::numeric_limits<double>::infinity();
                for (const FinalNonplanarWall &feature : walls) {
                    if (!feature.path->nonplanar_feature_transition ||
                        feature.path->inset_idx != surface.path->inset_idx)
                        continue;
                    const double contact_distance =
                        0.5 * (double(surface.path->width) +
                               double(feature.path->width)) + 0.05;
                    const double contact_distance_squared =
                        contact_distance * contact_distance;
                    const double feature_owner_z =
                        layers[feature.owner]->print_z;
                    const Points3 &feature_points =
                        feature.path->polyline.points;
                    for (size_t idx = 1; idx < feature_points.size(); ++idx) {
                        const Point3 &first = feature_points[idx - 1];
                        const Point3 &last = feature_points[idx];
                        const Vec2d a = unscale(first.to_point()).cast<double>();
                        const Vec2d b = unscale(last.to_point()).cast<double>();
                        const Vec2d delta = b - a;
                        const double squared = delta.squaredNorm();
                        const double ratio = squared <= EPSILON ? 0. :
                            std::clamp((sample - a).dot(delta) / squared,
                                       0., 1.);
                        if ((sample - (a + ratio * delta)).squaredNorm() >
                            contact_distance_squared)
                            continue;
                        const double feature_z = feature_owner_z +
                            unscale_(first.z()) + ratio *
                                unscale_(last.z() - first.z());
                        const double z_distance =
                            std::abs(feature_z - current_z);
                        if (z_distance < nearest_z_distance) {
                            nearest_z_distance = z_distance;
                            nearest_upper_z = feature_z;
                        }
                    }
                }
                const double maximum_z =
                    nearest_upper_z - minimum_feature_course_clearance;
                if (nearest_z_distance <= 0.08 &&
                    current_z > maximum_z + EPSILON) {
                    point.z() = scale_(maximum_z - owner_z);
                    changed = true;
                }
            }
            if (changed) {
                // A local ordering correction must become a ramp along the
                // existing wall, not a vertical notch at the edge of the
                // feature contact band. Propagate only downward so neither
                // the mesh envelope nor the already validated upper course
                // is exceeded.
                constexpr double maximum_transition_gradient = 0.8;
                for (size_t pass = 0;
                     pass < surface.path->polyline.points.size(); ++pass) {
                    bool smoothed = false;
                    for (size_t idx = 1;
                         idx < surface.path->polyline.points.size(); ++idx) {
                        Point3 &first = surface.path->polyline.points[idx - 1];
                        Point3 &last = surface.path->polyline.points[idx];
                        const double xy = unscale_(
                            (last.to_point() - first.to_point())
                                .cast<double>().norm());
                        if (xy <= EPSILON)
                            continue;
                        const double first_z = owner_z + unscale_(first.z());
                        const double last_z = owner_z + unscale_(last.z());
                        const double maximum_delta =
                            maximum_transition_gradient * xy;
                        if (first_z > last_z + maximum_delta + EPSILON) {
                            first.z() = scale_(last_z + maximum_delta - owner_z);
                            smoothed = true;
                        } else if (last_z > first_z + maximum_delta + EPSILON) {
                            last.z() = scale_(first_z + maximum_delta - owner_z);
                            smoothed = true;
                        }
                    }
                    if (!smoothed)
                        break;
                }
                surface.path->polyline.fitting_result.clear();
                surface.minimum_z = std::numeric_limits<double>::max();
                surface.maximum_z = std::numeric_limits<double>::lowest();
                for (const Point3 &point : surface.path->polyline.points) {
                    const double z = owner_z + unscale_(point.z());
                    surface.minimum_z = std::min(surface.minimum_z, z);
                    surface.maximum_z = std::max(surface.maximum_z, z);
                }
            }
        }
#endif

        // Surface-boundary adoption and deferred feature extraction can
        // produce numerically distinct copies of the same course. Their
        // centerlines may be separated by only a micron, so they do not
        // necessarily form a proper intersection until coordinates are
        // rounded for G-code. Deduplicate by sustained 3D boundary contact,
        // not by intersection alone, and retain the more complete contour.
        std::vector<uint8_t> unique_feature_wall(walls.size(), uint8_t(1));
        const auto coincident_feature_boundary =
            [&layers](const FinalNonplanarWall &reference,
                      const FinalNonplanarWall &probe) {
                const double reference_owner_z =
                    layers[reference.owner]->print_z;
                const double probe_owner_z = layers[probe.owner]->print_z;
                const double contact_distance =
                    0.5 * (double(reference.path->width) +
                           double(probe.path->width)) + 0.05;
                if (!reference.bounds.inflated(scale_(contact_distance))
                         .overlap(probe.bounds) ||
                    reference.minimum_z > probe.maximum_z + 0.08 ||
                    probe.minimum_z > reference.maximum_z + 0.08)
                    return false;
                const double contact_distance_squared =
                    contact_distance * contact_distance;
                size_t matching_points = 0;
                for (const Point3 &point : probe.path->polyline.points) {
                    const Vec2d sample =
                        unscale(point.to_point()).cast<double>();
                    bool matches = false;
                    for (size_t idx = 1;
                         idx < reference.path->polyline.points.size(); ++idx) {
                        const Point3 &first =
                            reference.path->polyline.points[idx - 1];
                        const Point3 &last =
                            reference.path->polyline.points[idx];
                        const Vec2d a =
                            unscale(first.to_point()).cast<double>();
                        const Vec2d b =
                            unscale(last.to_point()).cast<double>();
                        const Vec2d delta = b - a;
                        const double squared = delta.squaredNorm();
                        const double ratio = squared <= EPSILON ? 0. :
                            std::clamp((sample - a).dot(delta) / squared,
                                       0., 1.);
                        if ((sample - (a + ratio * delta)).squaredNorm() >
                            contact_distance_squared)
                            continue;
                        const double reference_z = reference_owner_z +
                            unscale_(first.z()) + ratio *
                                unscale_(last.z() - first.z());
                        const double probe_z = probe_owner_z +
                            unscale_(point.z());
                        if (std::abs(reference_z - probe_z) <= 0.08) {
                            matches = true;
                            break;
                        }
                    }
                    matching_points += matches;
                }
                return matching_points >= 2 &&
                    matching_points * 4 >=
                        3 * probe.path->polyline.points.size();
            };
        const auto covered_by_feature_boundaries =
            [&walls, &layers](const FinalNonplanarWall &conventional) {
                size_t covered_points = 0;
                // This is whole-fragment deduplication, not spatial clipping.
                // A feature touching one end of the dimensional outer wall
                // must not delete the remaining exposed roof edge. Inner and
                // interlocking fragments are deliberately more aggressive:
                // they are redundant when a transformed feature wall owns a
                // meaningful part of them, and retaining them creates a
                // same-Z crossing. Lower outer-wall coverage is left for the
                // physical-Z clipping pass below.
                const size_t point_count =
                    conventional.path->polyline.points.size();
                const size_t required_points = std::max<size_t>(
                    2, conventional.path->inset_idx == 0 ?
                        (3 * point_count + 3) / 4 :
                        (point_count + 3) / 4);
                const double conventional_owner_z =
                    layers[conventional.owner]->print_z;
                for (const Point3 &point : conventional.path->polyline.points) {
                    const Vec2d sample =
                        unscale(point.to_point()).cast<double>();
                    const double conventional_z = conventional_owner_z +
                        unscale_(point.z());
                    bool covered = false;
                    for (const FinalNonplanarWall &feature : walls) {
                        if (!feature.path->nonplanar_feature_transition ||
                            feature.path->polyline.points.size() < 2 ||
                            feature.path->inset_idx !=
                                conventional.path->inset_idx ||
                            !feature.bounds.overlap(conventional.bounds) ||
                            feature.minimum_z > conventional.maximum_z + 0.08 ||
                            conventional.minimum_z > feature.maximum_z + 0.08)
                            continue;
                        const double feature_owner_z =
                            layers[feature.owner]->print_z;
                        const double contact_distance =
                            0.5 * (double(feature.path->width) +
                                   double(conventional.path->width)) + 0.05;
                        for (size_t idx = 1;
                             idx < feature.path->polyline.points.size(); ++idx) {
                            const Point3 &first =
                                feature.path->polyline.points[idx - 1];
                            const Point3 &last =
                                feature.path->polyline.points[idx];
                            const Vec2d a =
                                unscale(first.to_point()).cast<double>();
                            const Vec2d b =
                                unscale(last.to_point()).cast<double>();
                            const Vec2d delta = b - a;
                            const double squared = delta.squaredNorm();
                            const double ratio = squared <= EPSILON ? 0. :
                                std::clamp((sample - a).dot(delta) / squared,
                                           0., 1.);
                            if ((sample - (a + ratio * delta)).norm() >
                                contact_distance)
                                continue;
                            const double feature_z = feature_owner_z +
                                unscale_(first.z()) + ratio *
                                    unscale_(last.z() - first.z());
                            if (std::abs(feature_z - conventional_z) <= 0.08) {
                                covered = true;
                                break;
                            }
                        }
                        if (covered)
                            break;
                    }
                    covered_points += covered;
                    if (covered_points >= required_points)
                        return true;
                }
                return false;
            };
        for (size_t left = 0; left < walls.size(); ++left) {
            if (!unique_feature_wall[left] ||
                !walls[left].path->nonplanar_feature_transition)
                continue;
            for (size_t right = left + 1; right < walls.size(); ++right) {
                if (!unique_feature_wall[right] ||
                    !walls[right].path->nonplanar_feature_transition ||
                    walls[left].path->nonplanar_feature_course !=
                        walls[right].path->nonplanar_feature_course ||
                    walls[left].path->inset_idx !=
                        walls[right].path->inset_idx)
                    continue;
                const size_t shorter =
                    walls[left].path->length() <= walls[right].path->length() ?
                        left : right;
                const size_t longer = shorter == left ? right : left;
                if (!coincident_feature_boundary(
                        walls[longer], walls[shorter]))
                    continue;
                walls[shorter].path->polyline.points.clear();
                unique_feature_wall[shorter] = 0;
                if (shorter == left)
                    break;
            }
        }
        // A successfully deferred feature course replaces the conventional
        // ring at the same physical height. Both can survive when the normal
        // continuation is owned by a later nominal layer than the surface
        // patch's local extraction range. Keeping both produces coincident
        // chimney/boss walls even though the transformed course already
        // supplies the complete dimensional boundary. Remove only a
        // sustained 3D match; partial contacts remain available to connect
        // neighboring planar geometry.
        for (FinalNonplanarWall &conventional : conventional_walls) {
            if (conventional.path->polyline.points.size() < 2)
                continue;
            for (const FinalNonplanarWall &feature : walls) {
                if (!feature.path->nonplanar_feature_transition ||
                    feature.path->polyline.points.size() < 2 ||
                    feature.path->inset_idx != conventional.path->inset_idx ||
                    !feature.bounds.overlap(conventional.bounds) ||
                    !coincident_feature_boundary(feature, conventional))
                    continue;
                conventional.path->polyline.points.clear();
                break;
            }
        }
        walls.erase(std::remove_if(
            walls.begin(), walls.end(), [](const FinalNonplanarWall &wall) {
                return wall.path->polyline.points.size() < 2;
            }), walls.end());

#if 0
        // Adjacent projected courses can approach within the G-code
        // collision tolerance along a short curved section even though their
        // complete loops are distinct and required for support. Resolve the
        // contact locally. The exposed surface remains the mesh authority, so
        // move its buried neighbor downward; between two buried courses move
        // the physically later course upward. The reference course is always
        // within the XY contact band, which keeps the adjusted point directly
        // supported instead of translating or deleting a complete loop.
        const auto average_wall_z = [&layers](const FinalNonplanarWall &wall) {
            const double owner_z = layers[wall.owner]->print_z;
            double sum = 0.;
            for (const Point3 &point : wall.path->polyline.points)
                sum += owner_z + unscale_(point.z());
            return sum / double(wall.path->polyline.points.size());
        };
        constexpr double local_course_xy_tolerance = 0.03;
        constexpr double local_course_clearance = 0.0601;
        for (size_t pass = 0; pass < 4; ++pass) {
            bool any_changed = false;
            for (size_t left = 0; left < walls.size(); ++left) {
                for (size_t right = left + 1; right < walls.size(); ++right) {
                    FinalNonplanarWall &first = walls[left];
                    FinalNonplanarWall &second = walls[right];
                    if (first.path->inset_idx != second.path->inset_idx ||
                        !first.bounds.overlap(second.bounds) ||
                        first.minimum_z > second.maximum_z + 0.08 ||
                        second.minimum_z > first.maximum_z + 0.08)
                        continue;

                    const bool first_exposed =
                        !first.path->nonplanar_transition &&
                        !first.path->nonplanar_feature_transition;
                    const bool second_exposed =
                        !second.path->nonplanar_transition &&
                        !second.path->nonplanar_feature_transition;
                    const bool keep_first_fixed = first_exposed != second_exposed ?
                        first_exposed : average_wall_z(first) <= average_wall_z(second);
                    FinalNonplanarWall &fixed =
                        keep_first_fixed ? first : second;
                    FinalNonplanarWall &adjusted =
                        keep_first_fixed ? second : first;
                    const bool adjusted_is_upper =
                        fixed.path->nonplanar_transition ||
                        fixed.path->nonplanar_feature_transition;
                    const double fixed_owner_z = layers[fixed.owner]->print_z;
                    const double adjusted_owner_z =
                        layers[adjusted.owner]->print_z;
                    bool changed = false;
                    for (Point3 &point : adjusted.path->polyline.points) {
                        const Vec2d sample =
                            unscale(point.to_point()).cast<double>();
                        double nearest_squared =
                            std::numeric_limits<double>::infinity();
                        double reference_z = 0.;
                        for (size_t idx = 1;
                             idx < fixed.path->polyline.points.size(); ++idx) {
                            const Point3 &a3 = fixed.path->polyline.points[idx - 1];
                            const Point3 &b3 = fixed.path->polyline.points[idx];
                            const Vec2d a = unscale(a3.to_point()).cast<double>();
                            const Vec2d delta =
                                unscale(b3.to_point()).cast<double>() - a;
                            const double squared = delta.squaredNorm();
                            const double ratio = squared <= EPSILON ? 0. :
                                std::clamp((sample - a).dot(delta) / squared,
                                           0., 1.);
                            const double distance_squared =
                                (sample - (a + ratio * delta)).squaredNorm();
                            if (distance_squared >= nearest_squared)
                                continue;
                            nearest_squared = distance_squared;
                            reference_z = fixed_owner_z + unscale_(a3.z()) +
                                ratio * unscale_(b3.z() - a3.z());
                        }
                        if (nearest_squared > local_course_xy_tolerance *
                                                  local_course_xy_tolerance)
                            continue;
                        const double current_z = adjusted_owner_z +
                            unscale_(point.z());
                        const double target_z = adjusted_is_upper ?
                            reference_z + local_course_clearance :
                            reference_z - local_course_clearance;
                        if ((adjusted_is_upper && current_z >= target_z - EPSILON) ||
                            (!adjusted_is_upper && current_z <= target_z + EPSILON))
                            continue;
                        point.z() = scale_(target_z - adjusted_owner_z);
                        changed = true;
                    }
                    if (!changed)
                        continue;
                    adjusted.path->polyline.fitting_result.clear();
                    adjusted.minimum_z = std::numeric_limits<double>::max();
                    adjusted.maximum_z = std::numeric_limits<double>::lowest();
                    for (const Point3 &point : adjusted.path->polyline.points) {
                        const double z = adjusted_owner_z + unscale_(point.z());
                        adjusted.minimum_z = std::min(adjusted.minimum_z, z);
                        adjusted.maximum_z = std::max(adjusted.maximum_z, z);
                    }
                    any_changed = true;
                }
            }
            if (!any_changed)
                break;
        }
#endif

#if 0
        // Two successive courses may follow virtually the same XY contour
        // without forming a mathematical intersection at full precision.
        // Millimetre G-code formatting can then turn the pair into coincident
        // or crossing extrusion. Order sustained near-coincident centerlines
        // with one uniform Z translation of the lower course. Uniform motion
        // preserves the already validated slope and endpoint blend; local
        // point edits would create a Z spike at the contact-band boundary.
        const auto wall_average_z = [&layers](const FinalNonplanarWall &wall) {
            const double owner_z = layers[wall.owner]->print_z;
            double total = 0.;
            for (const Point3 &point : wall.path->polyline.points)
                total += owner_z + unscale_(point.z());
            return total / double(wall.path->polyline.points.size());
        };
        std::vector<size_t> physical_wall_order(walls.size());
        std::iota(physical_wall_order.begin(), physical_wall_order.end(), size_t(0));
        std::stable_sort(
            physical_wall_order.begin(), physical_wall_order.end(),
            [&walls, &wall_average_z](size_t left, size_t right) {
                return wall_average_z(walls[left]) >
                       wall_average_z(walls[right]);
            });
        constexpr double coincident_xy_tolerance = 0.03;
        constexpr double rounded_course_clearance = 0.0601;
        for (size_t high_order = 0;
             high_order < physical_wall_order.size(); ++high_order) {
            FinalNonplanarWall &higher =
                walls[physical_wall_order[high_order]];
            const double higher_owner_z = layers[higher.owner]->print_z;
            for (size_t low_order = high_order + 1;
                 low_order < physical_wall_order.size(); ++low_order) {
                FinalNonplanarWall &lower =
                    walls[physical_wall_order[low_order]];
                if (higher.path->inset_idx != lower.path->inset_idx)
                    continue;
                bool has_nearer_course = false;
                for (size_t middle_order = high_order + 1;
                     middle_order < low_order; ++middle_order) {
                    const FinalNonplanarWall &middle =
                        walls[physical_wall_order[middle_order]];
                    if (middle.path->inset_idx == lower.path->inset_idx &&
                        middle.bounds.overlap(lower.bounds)) {
                        has_nearer_course = true;
                        break;
                    }
                }
                if (has_nearer_course)
                    continue;
                const FinalNonplanarWall &probe =
                    higher.path->polyline.points.size() <=
                            lower.path->polyline.points.size() ?
                        higher : lower;
                const FinalNonplanarWall &reference =
                    &probe == &higher ? lower : higher;
                size_t coincident_points = 0;
                for (const Point3 &point : probe.path->polyline.points) {
                    const Vec2d sample =
                        unscale(point.to_point()).cast<double>();
                    double nearest_squared =
                        std::numeric_limits<double>::infinity();
                    for (size_t idx = 1;
                         idx < reference.path->polyline.points.size(); ++idx) {
                        const Vec2d a = unscale(
                            reference.path->polyline.points[idx - 1].to_point())
                                                .cast<double>();
                        const Vec2d b = unscale(
                            reference.path->polyline.points[idx].to_point())
                                                .cast<double>();
                        const Vec2d delta = b - a;
                        const double squared = delta.squaredNorm();
                        const double ratio = squared <= EPSILON ? 0. :
                            std::clamp((sample - a).dot(delta) / squared,
                                       0., 1.);
                        nearest_squared = std::min(
                            nearest_squared,
                            (sample - (a + ratio * delta)).squaredNorm());
                    }
                    coincident_points += nearest_squared <=
                        coincident_xy_tolerance * coincident_xy_tolerance;
                }
                if (coincident_points < 2 ||
                    coincident_points * 2 < probe.path->polyline.points.size())
                    continue;

                const double lower_owner_z = layers[lower.owner]->print_z;
                double required_upward_shift = 0.;
                for (const Point3 &point : lower.path->polyline.points) {
                    const Vec2d sample =
                        unscale(point.to_point()).cast<double>();
                    double nearest_squared =
                        std::numeric_limits<double>::infinity();
                    double nearest_z = 0.;
                    for (size_t idx = 1;
                         idx < higher.path->polyline.points.size(); ++idx) {
                        const Point3 &first = higher.path->polyline.points[idx - 1];
                        const Point3 &last = higher.path->polyline.points[idx];
                        const Vec2d a = unscale(first.to_point()).cast<double>();
                        const Vec2d b = unscale(last.to_point()).cast<double>();
                        const Vec2d delta = b - a;
                        const double squared = delta.squaredNorm();
                        const double ratio = squared <= EPSILON ? 0. :
                            std::clamp((sample - a).dot(delta) / squared,
                                       0., 1.);
                        const double distance_squared =
                            (sample - (a + ratio * delta)).squaredNorm();
                        if (distance_squared >= nearest_squared)
                            continue;
                        nearest_squared = distance_squared;
                        nearest_z = higher_owner_z + unscale_(first.z()) +
                            ratio * unscale_(last.z() - first.z());
                    }
                    if (nearest_squared >
                        coincident_xy_tolerance * coincident_xy_tolerance)
                        continue;
                    const double lower_z = lower_owner_z + unscale_(point.z());
                    required_upward_shift = std::max(
                        required_upward_shift,
                        lower_z + rounded_course_clearance - nearest_z);
                }
                if (required_upward_shift <= EPSILON)
                    continue;
                // Keep the supporting course fixed and move the later course
                // just beyond the rounded contact tolerance. Lowering the
                // support created a larger gap to the course beneath it;
                // raising the already adjacent upper course preserves its
                // deposited support and changes Z by at most one contact
                // tolerance.
                for (Point3 &point : higher.path->polyline.points)
                    point.z() += scale_(required_upward_shift);
                higher.path->polyline.fitting_result.clear();
                higher.minimum_z += required_upward_shift;
                higher.maximum_z += required_upward_shift;
            }
        }

#endif
        // The finishing skin is the model-surface authority. A deferred
        // chimney/boss course which still crosses that skin must move upward
        // in dependency order; lowering the exposed roof to make room for it
        // collapses the preceding transition stack. Repair one crossing at a
        // time and let the wall-stack pass below propagate the new lower bound
        // through subsequent feature courses.
        const size_t maximum_feature_skin_repairs =
            std::max<size_t>(16, 2 * walls.size());
        for (size_t repair = 0;
             repair < maximum_feature_skin_repairs; ++repair) {
            bool repaired = false;
            for (FinalNonplanarWall &feature : walls) {
                if (!feature.path->nonplanar_feature_transition ||
                    feature.path->polyline.points.size() < 2)
                    continue;
                for (const FinalNonplanarWall &top : top_paths) {
                    if (!feature.bounds.overlap(top.bounds) ||
                        feature.minimum_z > top.maximum_z + 0.05 ||
                        top.minimum_z > feature.maximum_z + 0.05 ||
                        !paths_have_proper_same_z_crossing(
                            *feature.path, layers[feature.owner]->print_z,
                            *top.path, layers[top.owner]->print_z, 0.05))
                        continue;
                    constexpr double feature_skin_clearance = 0.061;
                    shift_wall(feature, feature_skin_clearance);
                    repaired = true;
                    break;
                }
                if (repaired)
                    break;
            }
            if (!repaired)
                break;
        }

        // A contoured conventional stair boundary is not part of the model's
        // exposed edge. If it passes through an accepted finishing raster at
        // the same physical Z, the raster owns that interior surface and the
        // stair fragment must not be emitted. This final graph check handles
        // mixed loops whose nominal owner differs from their contoured Z;
        // nominal-layer polygon clipping cannot see that collision.
        const coord_t top_bucket_width = scale_(2.);
        std::map<int64_t, std::vector<const FinalNonplanarWall *>>
            top_paths_by_x;
        for (const FinalNonplanarWall &top : top_paths) {
            const int64_t first_bucket =
                top.bounds.min.x() / top_bucket_width;
            const int64_t last_bucket =
                top.bounds.max.x() / top_bucket_width;
            for (int64_t bucket = first_bucket;
                 bucket <= last_bucket; ++bucket)
                top_paths_by_x[bucket].push_back(&top);
        }
        for (FinalNonplanarWall &wall : conventional_walls) {
            std::set<const FinalNonplanarWall *> candidates;
            const int64_t first_bucket =
                wall.bounds.min.x() / top_bucket_width;
            const int64_t last_bucket =
                wall.bounds.max.x() / top_bucket_width;
            for (int64_t bucket = first_bucket;
                 bucket <= last_bucket; ++bucket) {
                const auto found = top_paths_by_x.find(bucket);
                if (found != top_paths_by_x.end())
                    candidates.insert(found->second.begin(),
                                      found->second.end());
            }
            for (const FinalNonplanarWall *top : candidates) {
                if (!wall.bounds.overlap(top->bounds) ||
                    wall.minimum_z > top->maximum_z + 0.05 ||
                    top->minimum_z > wall.maximum_z + 0.05 ||
                    !paths_have_proper_same_z_crossing(
                        *wall.path, layers[wall.owner]->print_z,
                        *top->path, layers[top->owner]->print_z, 0.05))
                    continue;
                wall.path->polyline.points.clear();
                break;
            }
        }
        walls.erase(std::remove_if(
            walls.begin(), walls.end(), [](const FinalNonplanarWall &wall) {
                return wall.path->polyline.points.size() < 2;
            }), walls.end());
        // Correct one contact and restart the scan. A wall may intersect more
        // than one adjacent course, and excluding it after the first repair
        // leaves the second contact in the emitted G-code. Deferred feature
        // courses form a directed stack, so later courses are moved only
        // upward; same-course duplicate fragments are removed. Rechecking
        // after every edit catches contacts with the next course. The course
        // count bounds convergence without turning malformed input into an
        // unbounded repair loop.
        const size_t maximum_repairs = std::max<size_t>(32, 4 * walls.size());
        // Repairs below only change Z or clear a path; XY bounds and vector
        // indices remain fixed. Cache the broad-phase pairs once, in exactly
        // the original scan order, instead of repeating the same XY tests
        // after every contact. Z and exact crossings must still be rechecked.
        std::vector<std::vector<size_t>> wall_candidates(walls.size());
        for (size_t left = 0; left < walls.size(); ++left)
            for (size_t right = left + 1; right < walls.size(); ++right)
                if (bounds_overlap_with_tolerance(
                        walls[left].bounds, walls[right].bounds, scale_(0.05)))
                    wall_candidates[left].push_back(right);
        for (size_t repair = 0; repair < maximum_repairs; ++repair) {
            bool repaired = false;
            for (size_t left = 0; left < walls.size() && !repaired; ++left) {
                if (walls[left].path->polyline.points.size() < 2)
                    continue;
                for (size_t right : wall_candidates[left]) {
                    if (walls[right].path->polyline.points.size() < 2 ||
                        walls[left].minimum_z > walls[right].maximum_z + 0.05 ||
                        walls[right].minimum_z > walls[left].maximum_z + 0.05 ||
                        !paths_have_proper_same_z_crossing(
                            *walls[left].path, layers[walls[left].owner]->print_z,
                            *walls[right].path, layers[walls[right].owner]->print_z,
                            0.05))
                        continue;
                    const bool duplicate_course =
                        walls[left].path->nonplanar_feature_course != 0 &&
                        walls[left].path->nonplanar_feature_course ==
                            walls[right].path->nonplanar_feature_course &&
                        walls[left].path->inset_idx ==
                            walls[right].path->inset_idx;
                    size_t victim = right;
                    if (duplicate_course) {
                        victim = walls[left].path->length() <=
                                         walls[right].path->length() ?
                            left : right;
                        walls[victim].path->polyline.points.clear();
                    } else if (walls[left].path->inset_idx ==
                                   walls[right].path->inset_idx &&
                               (walls[left].path->nonplanar_feature_course != 0 ||
                                walls[right].path->nonplanar_feature_course != 0) &&
                               walls[left].path->nonplanar_feature_course !=
                                   walls[right].path->nonplanar_feature_course) {
                        // Deferred courses are a directed stack. Repairing a
                        // contact by lowering an arbitrary average-Z path can
                        // make course N collide with N-2 and oscillate. Move
                        // only the later course upward; repeated passes then
                        // converge monotonically while retaining every real
                        // feature wall.
                        victim = walls[left].path->nonplanar_feature_course >
                                         walls[right].path->nonplanar_feature_course ?
                            left : right;
                        constexpr double repair_clearance = 0.061;
                        shift_wall(walls[victim], repair_clearance);
                    } else {
                        // Two unrelated paths are not an ordered feature
                        // stack. A proper centerline crossing is a duplicate
                        // fragment, so retain the more complete wall instead
                        // of displacing valid model geometry.
                        victim = walls[left].path->length() <=
                                         walls[right].path->length() ?
                            left : right;
                        walls[victim].path->polyline.points.clear();
                    }
                    ++final_crossing_count;
                    repaired = true;
                    break;
                }
            }
            if (!repaired)
                break;
        }
        // The crossing repair above may lift a course and leave too little
        // room below the following one. At this late stage the finishing skin
        // is immutable, so propagate the required separation upward through
        // later feature courses rather than lowering a course back through
        // the roof. Process courses in dependency order; a shift made here is
        // therefore included when the next course is evaluated.
        std::vector<size_t> ascending_feature_order(walls.size());
        std::iota(
            ascending_feature_order.begin(), ascending_feature_order.end(),
            size_t(0));
        std::stable_sort(
            ascending_feature_order.begin(), ascending_feature_order.end(),
            [&walls](size_t left, size_t right) {
                return walls[left].path->nonplanar_feature_course <
                       walls[right].path->nonplanar_feature_course;
            });
        for (const size_t higher_idx : ascending_feature_order) {
            FinalNonplanarWall &higher = walls[higher_idx];
            if (!higher.path->nonplanar_feature_transition ||
                higher.path->nonplanar_feature_course == 0 ||
                higher.path->polyline.points.size() < 2)
                continue;
            size_t previous_course = 0;
            for (const FinalNonplanarWall &candidate : walls)
                if (candidate.path->nonplanar_feature_transition &&
                    candidate.path->inset_idx == higher.path->inset_idx &&
                    candidate.path->nonplanar_feature_course <
                        higher.path->nonplanar_feature_course)
                    previous_course = std::max(
                        previous_course,
                        candidate.path->nonplanar_feature_course);
            if (previous_course == 0)
                continue;
            const double higher_owner_z = layers[higher.owner]->print_z;
            double required_shift = 0.;
            const auto evaluate_sample = [&](const Vec2d &sample,
                                             double higher_z) {
                for (const FinalNonplanarWall &lower : walls) {
                    if (!lower.path->nonplanar_feature_transition ||
                        lower.path->nonplanar_feature_course != previous_course ||
                        lower.path->inset_idx != higher.path->inset_idx ||
                        lower.path->polyline.points.size() < 2)
                        continue;
                    const double contact_distance =
                        0.5 * (double(lower.path->width) +
                               double(higher.path->width)) + 0.05;
                    const double contact_distance_squared =
                        contact_distance * contact_distance;
                    const double lower_owner_z = layers[lower.owner]->print_z;
                    const Points3 &lower_points = lower.path->polyline.points;
                    for (size_t point_idx = 1;
                         point_idx < lower_points.size(); ++point_idx) {
                        const Point3 &first = lower_points[point_idx - 1];
                        const Point3 &last = lower_points[point_idx];
                        const Vec2d a = unscale(first.to_point()).cast<double>();
                        const Vec2d b = unscale(last.to_point()).cast<double>();
                        const Vec2d delta = b - a;
                        const double squared = delta.squaredNorm();
                        const double ratio = squared <= EPSILON ? 0. :
                            std::clamp((sample - a).dot(delta) / squared,
                                       0., 1.);
                        if ((sample - (a + ratio * delta)).squaredNorm() >
                            contact_distance_squared)
                            continue;
                        const double lower_z = lower_owner_z +
                            unscale_(first.z()) + ratio *
                                unscale_(last.z() - first.z());
                        const double clearance = std::max(
                            minimum_feature_course_clearance,
                            0.5 * std::max(double(lower.path->height),
                                           double(higher.path->height)) +
                                0.001);
                        required_shift = std::max(
                            required_shift, lower_z + clearance - higher_z);
                    }
                }
            };
            const Points3 &higher_points = higher.path->polyline.points;
            for (const Point3 &point : higher_points)
                evaluate_sample(
                    unscale(point.to_point()).cast<double>(),
                    higher_owner_z + unscale_(point.z()));
            // Curved rings are resampled differently from one nominal course
            // to the next. Vertex-only comparison can therefore miss the
            // closest approach in the middle of a long segment. Include the
            // same segment midpoints used by emitted-path validation.
            for (size_t point_idx = 1;
                 point_idx < higher_points.size(); ++point_idx) {
                const Point3 &first = higher_points[point_idx - 1];
                const Point3 &last = higher_points[point_idx];
                evaluate_sample(
                    0.5 * (unscale(first.to_point()).cast<double>() +
                           unscale(last.to_point()).cast<double>()),
                    higher_owner_z +
                        0.5 * (unscale_(first.z()) + unscale_(last.z())));
            }
            if (required_shift <= EPSILON)
                continue;
            shift_wall(higher, required_shift);
        }
        // The dependency propagation above enforces bead support against the
        // immediately preceding course, but that upward translation can make
        // the settled course cross a non-adjacent course on a curved feature.
        // Re-run same-inset arbitration after propagation. Always move the
        // later dependency course upward, preserving the already deposited
        // support and converging monotonically instead of deleting a real
        // chimney/boss wall.
        const size_t maximum_settled_repairs =
            std::max<size_t>(32, 4 * walls.size());
        for (size_t repair = 0; repair < maximum_settled_repairs; ++repair) {
            bool repaired = false;
            for (size_t left = 0; left < walls.size() && !repaired; ++left) {
                if (walls[left].path->polyline.points.size() < 2)
                    continue;
                for (size_t right = left + 1; right < walls.size(); ++right) {
                    if (walls[right].path->polyline.points.size() < 2 ||
                        walls[left].path->inset_idx !=
                            walls[right].path->inset_idx ||
                        !bounds_overlap_with_tolerance(
                            walls[left].bounds, walls[right].bounds,
                            scale_(0.05)) ||
                        walls[left].minimum_z > walls[right].maximum_z + 0.05 ||
                        walls[right].minimum_z > walls[left].maximum_z + 0.05 ||
                        !paths_have_proper_same_z_crossing(
                            *walls[left].path, layers[walls[left].owner]->print_z,
                            *walls[right].path, layers[walls[right].owner]->print_z,
                            0.05))
                        continue;
                    FinalNonplanarWall &later =
                        walls[left].path->nonplanar_feature_transition !=
                                walls[right].path->nonplanar_feature_transition ?
                            (walls[left].path->nonplanar_feature_transition ?
                                walls[left] : walls[right]) :
                            (walls[left].path->nonplanar_feature_course >
                                    walls[right].path->nonplanar_feature_course ?
                                walls[left] : walls[right]);
                    constexpr double settled_repair_clearance = 0.061;
                    shift_wall(later, settled_repair_clearance);
                    ++final_crossing_count;
                    repaired = true;
                    break;
                }
            }
            if (!repaired)
                break;
        }
        // Brick/interlocking metadata can give neighboring insets different
        // height profiles. The upward dependency propagation above may make
        // two such fragments cross even though each inset's own course stack
        // is valid. They are not an ordered vertical pair, so retain the more
        // complete dimensional boundary and discard only the shorter
        // crossing fragment. This is intentionally after height settlement.
        const size_t maximum_cross_inset_repairs =
            std::max<size_t>(16, 2 * walls.size());
        for (size_t repair = 0;
             repair < maximum_cross_inset_repairs; ++repair) {
            bool repaired = false;
            for (size_t left = 0; left < walls.size() && !repaired; ++left) {
                if (walls[left].path->polyline.points.size() < 2)
                    continue;
                for (size_t right = left + 1; right < walls.size(); ++right) {
                    if (walls[right].path->polyline.points.size() < 2 ||
                        walls[left].path->inset_idx ==
                            walls[right].path->inset_idx ||
                        !walls[left].bounds.overlap(walls[right].bounds) ||
                        walls[left].minimum_z > walls[right].maximum_z + 0.05 ||
                        walls[right].minimum_z > walls[left].maximum_z + 0.05 ||
                        !paths_have_proper_same_z_crossing(
                            *walls[left].path, layers[walls[left].owner]->print_z,
                            *walls[right].path, layers[walls[right].owner]->print_z,
                            0.05))
                        continue;
                    const size_t victim =
                        walls[left].path->length() <=
                                walls[right].path->length() ?
                            left : right;
                    walls[victim].path->polyline.points.clear();
                    ++final_crossing_count;
                    repaired = true;
                    break;
                }
            }
            if (!repaired)
                break;
        }
        // Collision arbitration may move a deferred feature course upward
        // until it becomes the physical replacement for a conventional ring
        // that was not coincident during the first handoff pass. Re-run the
        // same sustained-boundary test on the settled graph so the ordinary
        // ring is not emitted a few hundredths of a millimetre above the
        // transformed wall.
        for (FinalNonplanarWall &conventional : conventional_walls) {
            if (conventional.path->polyline.points.size() < 2)
                continue;
            if (covered_by_feature_boundaries(conventional)) {
                conventional.path->polyline.points.clear();
                continue;
            }
            for (const FinalNonplanarWall &feature : walls) {
                if (!feature.path->nonplanar_feature_transition ||
                    feature.path->polyline.points.size() < 2 ||
                    feature.path->inset_idx != conventional.path->inset_idx ||
                    !feature.bounds.overlap(conventional.bounds) ||
                    !coincident_feature_boundary(feature, conventional))
                    continue;
                conventional.path->polyline.points.clear();
                break;
            }
        }
        // Z corrections may bring a transformed feature boundary into a
        // conventional stair fragment after the patch-local handoff clipping
        // has run. Reapply that same physical-Z material subtraction to the
        // settled graph. Accepted non-planar entities are preserved by the
        // clipping helper; only conventional centerline portions occupying a
        // transformed bead envelope are removed, so a partial planar wall is
        // joined to the feature instead of deleting its complete loop.
        // Build every mask before clipping: `walls` points into collections
        // which clipping replaces. Brick courses occupy their raised physical
        // plane, not the nominal owner plane used by ordinary perimeters.
        std::vector<std::map<coord_t, ExPolygons>> settled_feature_coverage(layers.size());
        // A 2D mask at one Z cannot represent an interlocking wall. Subtract
        // its actual contact intervals instead, retaining XYZ and extrusion
        // metadata and joining consecutive untouched pieces back together.
        std::map<const ExtrusionPath *, ExtrusionPaths> variable_wall_replacements;
        std::vector<bool> variable_wall_changed(layers.size(), false);
        for (const FinalNonplanarWall &wall : conventional_walls) {
            const ExtrusionPath &path = *wall.path;
            if (!path.z_contoured || path.polyline.points.size() < 2 ||
                std::all_of(path.polyline.points.begin(), path.polyline.points.end(),
                    [&](const Point3 &point) { return point.z() == path.polyline.points.front().z(); }))
                continue;
            ExtrusionPaths fragments;
            bool changed = false;
            for (size_t idx = 1; idx < path.polyline.points.size(); ++idx) {
                const Point3 &first = path.polyline.points[idx - 1];
                const Point3 &last = path.polyline.points[idx];
                Vec3d a = first.cast<double>() * SCALING_FACTOR;
                Vec3d b = last.cast<double>() * SCALING_FACTOR;
                a.z() += layers[wall.owner]->print_z;
                b.z() += layers[wall.owner]->print_z;
                std::vector<std::pair<double, double>> covered;
                for (const FinalNonplanarWall &feature : walls) {
                    const bool coincident = feature.path->nonplanar_feature_transition ||
                        (feature.path->nonplanar_surface &&
                         (path.staggered_perimeter || path.staggered_transition || path.nonplanar_interlocking_wall));
                    const bool underpass = feature.path->nonplanar_surface &&
                        (feature.owner > wall.owner ||
                         (feature.owner == wall.owner && !feature.path->nonplanar_before_current_layer));
                    const double maximum_underpass = layers[wall.owner]->height - 0.02;
                    if ((!coincident && !underpass) ||
                        feature.maximum_z < wall.minimum_z -
                            (underpass ? std::max(0.05, maximum_underpass) : 0.05) ||
                        feature.minimum_z > wall.maximum_z + 0.05 ||
                        !feature.bounds.inflated(scale_(0.5 * feature.path->width + 0.05)).overlap(wall.bounds))
                        continue;
                    const Points3 &points = feature.path->polyline.points;
                    for (size_t j = 1; j < points.size(); ++j) {
                        Vec3d c = points[j - 1].cast<double>() * SCALING_FACTOR;
                        Vec3d d = points[j].cast<double>() * SCALING_FACTOR;
                        c.z() += layers[feature.owner]->print_z;
                        d.z() += layers[feature.owner]->print_z;
                        const double xy_tolerance = 0.5 * feature.path->width + 0.05;
                        if (coincident) {
                            const auto intervals = nonplanar_segment_overlap(a, b, c, d,
                                xy_tolerance, 0.05);
                            covered.insert(covered.end(), intervals.begin(), intervals.end());
                        }
                        if (underpass && maximum_underpass > 0.05) {
                            Vec3d lower_a = a, lower_b = b;
                            lower_a.z() -= 0.5 * (maximum_underpass + 0.05);
                            lower_b.z() -= 0.5 * (maximum_underpass + 0.05);
                            const auto intervals = nonplanar_segment_overlap(lower_a, lower_b, c, d,
                                xy_tolerance, 0.5 * (maximum_underpass - 0.05));
                            covered.insert(covered.end(), intervals.begin(), intervals.end());
                        }
                    }
                }
                std::sort(covered.begin(), covered.end());
                const auto append_interval = [&](double begin, double end) {
                    if (end <= begin + 1e-12)
                        return;
                    const auto point_at = [&](double t) -> Point3 {
                        const Vec3d point = first.cast<double>() + t * (last - first).cast<double>();
                        return Point3(coord_t(std::llround(point.x())),
                                      coord_t(std::llround(point.y())),
                                      coord_t(std::llround(point.z())));
                    };
                    const Point3 start = point_at(begin), finish = point_at(end);
                    if (start == finish)
                        return;
                    if (!fragments.empty() && fragments.back().polyline.points.back() == start)
                        fragments.back().polyline.points.push_back(finish);
                    else {
                        fragments.emplace_back(path);
                        fragments.back().polyline.points = {start, finish};
                        fragments.back().polyline.fitting_result.clear();
                    }
                };
                double cursor = 0.;
                for (const auto &[begin, end] : covered) {
                    append_interval(cursor, begin);
                    cursor = std::max(cursor, end);
                    changed = true;
                }
                append_interval(cursor, 1.);
            }
            if (changed) {
                // A cut away from the old seam leaves one surviving arc
                // across that seam, not two disconnected extrusion paths.
                if (fragments.size() > 1 &&
                    fragments.back().polyline.points.back() == fragments.front().polyline.points.front()) {
                    Points3 joined = std::move(fragments.back().polyline.points);
                    const Points3 &head = fragments.front().polyline.points;
                    joined.insert(joined.end(), std::next(head.begin()), head.end());
                    fragments.front().polyline.points = std::move(joined);
                    fragments.pop_back();
                }
                for (ExtrusionPath &fragment : fragments)
                    fragment.nonplanar_replacement_remainder = true;
                variable_wall_replacements.emplace(&path, std::move(fragments));
                variable_wall_changed[wall.owner] = true;
            }
        }
        const auto course_offset = [](const ExtrusionPath &path) -> coord_t {
            if (!path.z_contoured || path.polyline.points.empty())
                return 0;
            const coord_t offset = path.polyline.points.front().z();
            return std::all_of(path.polyline.points.begin(), path.polyline.points.end(),
                [offset](const Point3 &point) { return point.z() == offset; }) ? offset : 0;
        };
        for (const FinalNonplanarWall &wall : conventional_walls)
            settled_feature_coverage[wall.owner].try_emplace(course_offset(*wall.path));
        conventional_walls.clear();
        for (size_t owner = 0; owner < layers.size(); ++owner) {
            for (auto &[offset, coverage] : settled_feature_coverage[owner]) {
                Polygons material;
                const double source_z = layers[owner]->print_z + unscale_(offset);
                for (const FinalNonplanarWall &feature : walls) {
                    if ((!feature.path->nonplanar_feature_transition &&
                         !(offset != 0 && feature.path->nonplanar_surface)) ||
                        feature.path->polyline.points.size() < 2)
                        continue;
                    collect_generated_material_at_z(
                        *feature.path, layers[feature.owner]->print_z, source_z,
                        0.05, false, material);
                }
                if (!material.empty())
                    coverage = offset_ex(union_ex(material), scale_(0.05));
                // Test the settled graph, after feature extraction and wall
                // adoption. A conventional stair wall printed before a lower
                // draped wall must not remain inside its bead footprint. Do
                // not apply this subtraction to a normal supporting course
                // below the draped path or to a later, correctly stacked wall.
                Polygons underpass_material;
                const double maximum_underpass = layers[owner]->height - 0.02;
                if (maximum_underpass > 0.05)
                    for (const FinalNonplanarWall &surface : walls) {
                        if (!surface.path->nonplanar_surface ||
                            surface.owner < owner ||
                            (surface.owner == owner && surface.path->nonplanar_before_current_layer))
                            continue;
                        // Match the finite-endpoint radius of the XYZ overlap
                        // test. Butt caps leave stair-wall tips just beyond a
                        // draped segment's endpoint inside the nozzle footprint.
                        collect_generated_material_at_z(*surface.path,
                            layers[surface.owner]->print_z,
                            source_z - 0.5 * (maximum_underpass + 0.05),
                            0.5 * (maximum_underpass - 0.05), false,
                            underpass_material, true);
                    }
                if (!underpass_material.empty()) {
                    append(coverage, offset_ex(union_ex(underpass_material), scale_(0.05)));
                    coverage = union_ex(coverage);
                }
            }
        }
        for (size_t owner = 0; owner < layers.size(); ++owner) {
            const std::function<ExPolygons(const ExtrusionPath &)> coverage_at_path =
                [&](const ExtrusionPath &path) {
                    if (path.z_contoured && !path.polyline.points.empty() &&
                        std::any_of(path.polyline.points.begin(), path.polyline.points.end(),
                            [&](const Point3 &point) { return point.z() != path.polyline.points.front().z(); }))
                        return ExPolygons{}; // Already checked against the physical 3D segments.
                    const auto found = settled_feature_coverage[owner].find(course_offset(path));
                    return found == settled_feature_coverage[owner].end() ? ExPolygons{} : found->second;
                };
            if (variable_wall_changed[owner] || std::any_of(settled_feature_coverage[owner].begin(),
                            settled_feature_coverage[owner].end(),
                            [](const auto &entry) { return !entry.second.empty(); }))
                remove_nonplanar_top_overlap(
                    layers[owner]->regions()[region_idx]->perimeters,
                    {}, false, true, &coverage_at_path, &variable_wall_replacements);
        }

        for (size_t owner = 0; owner < layers.size(); ++owner) {
            prune_empty_extrusion_paths(
                layers[owner]->regions()[region_idx]->perimeters);
            prune_empty_extrusion_paths(
                layers[owner]->regions()[region_idx]->fills);
        }
        (void) final_crossing_count;
    }
    // Runs and continuing features may span different nominal owner ranges.
    // Normalize the completed perimeter graph once after all insertions so a
    // late transformed course cannot escape the selected wall order.
    for (size_t owner = 0; owner < layers.size(); ++owner)
        for (LayerRegion *region : layers[owner]->regions()) {
            const PrintRegionConfig &config = region->region().config();
            WallSequence wall_sequence = config.wall_sequence.value;
            if (config.wall_loops.value >= 3 &&
                config.inner_walls_flow_ratio.get_abs_value(1.) >
                    1.0 + EPSILON && owner > 0)
                wall_sequence = WallSequence::InnerOuterInner;
            restore_nonplanar_wall_sequence(
                region->perimeters, wall_sequence,
                config.wall_loops.value, owner > 0);
            // Wall-order reconstruction is the final graph mutation and may
            // split/copy a mixed loop after the collision pass normalized it.
            // Reassert the fixed-Z feature invariant on the reconstructed
            // entities before G-code scheduling and preview consume them.
            for (ExtrusionEntity *entity : region->perimeters.entities)
                normalize_fixed_feature_transition(*entity);
        }
    // Wall ordering is intentionally the last structural graph operation.
    // Revalidate feature-to-surface contacts afterwards so no late grouping
    // can reintroduce an emitted-precision collision between an attached
    // feature course and the retained opening boundary.
    for (size_t region_idx = 0; region_idx < region_count; ++region_idx) {
        struct OwnedWallPath {
            ExtrusionPath *path;
            size_t owner;
            BoundingBox bounds;
            ExtrusionEntity *continuous_wall;
        };
        std::vector<OwnedWallPath> feature_walls;
        std::vector<OwnedWallPath> surface_walls;
        std::vector<OwnedWallPath> foundation_walls;
        const auto collect = [&](auto &&self, ExtrusionEntity &entity,
                                 size_t owner, ExtrusionEntity *continuous_wall) -> void {
            if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
                // A horizontal portion of a deferred feature is displayed as
                // an ordinary wall, but still belongs to the same physical
                // course stack and must participate in clearance validation.
                if (!is_perimeter(path->role()) ||
                    path->polyline.points.size() < 2)
                    return;
                OwnedWallPath owned{
                    path, owner, path->polyline.to_polyline().bounding_box(),
                    continuous_wall ? continuous_wall : path};
                if (path->nonplanar_feature_transition)
                    feature_walls.emplace_back(std::move(owned));
                else if (path->nonplanar_surface)
                    surface_walls.emplace_back(std::move(owned));
                else if (!path->nonplanar_schedule_owned)
                    foundation_walls.emplace_back(std::move(owned));
            } else if (auto *multipath =
                           dynamic_cast<ExtrusionMultiPath *>(&entity)) {
                for (ExtrusionPath &path : multipath->paths)
                    self(self, path, owner, &entity);
            } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
                for (ExtrusionPath &path : loop->paths)
                    self(self, path, owner, &entity);
            } else if (auto *collection =
                           dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
                for (ExtrusionEntity *child : collection->entities)
                    self(self, *child, owner, nullptr);
            }
        };
        for (size_t owner = 0; owner < layers.size(); ++owner) {
            LayerRegion &region = *layers[owner]->regions()[region_idx];
            collect(collect, region.perimeters, owner, nullptr);
            collect(collect, region.fills, owner, nullptr);
        }
        for (OwnedWallPath &feature : feature_walls) {
            for (size_t repair = 0; repair < 8; ++repair) {
                bool collision = false;
                for (const OwnedWallPath &surface : surface_walls) {
                    if (feature.path->inset_idx != surface.path->inset_idx ||
                        !bounds_overlap_with_tolerance(
                            feature.bounds, surface.bounds, scale_(0.05)) ||
                        !paths_have_proper_same_z_crossing(
                            *feature.path, layers[feature.owner]->print_z,
                            *surface.path, layers[surface.owner]->print_z,
                            0.05))
                        continue;
                    collision = true;
                    break;
                }
                if (!collision)
                    break;
                offset_entity_z(*feature.continuous_wall, scale_(0.061));
                normalize_fixed_feature_transition(*feature.continuous_wall);
            }
        }
        std::stable_sort(feature_walls.begin(), feature_walls.end(),
            [](const OwnedWallPath &left, const OwnedWallPath &right) {
                return left.path->nonplanar_feature_course <
                       right.path->nonplanar_feature_course;
            });
        for (OwnedWallPath &higher : feature_walls) {
            size_t previous_course = 0;
            for (const OwnedWallPath &candidate : feature_walls)
                if (candidate.path->inset_idx == higher.path->inset_idx &&
                    candidate.path->nonplanar_feature_course <
                        higher.path->nonplanar_feature_course)
                    previous_course = std::max(
                        previous_course,
                        candidate.path->nonplanar_feature_course);
            if (previous_course == 0)
                continue;
            const double higher_owner_z = layers[higher.owner]->print_z;
            double required_shift = 0.;
            for (size_t point_idx = 1;
                 point_idx < higher.path->polyline.points.size(); ++point_idx) {
                const Point3 &first = higher.path->polyline.points[point_idx - 1];
                const Point3 &last = higher.path->polyline.points[point_idx];
                const Vec2d sample = 0.5 * (
                    unscale(first.to_point()).cast<double>() +
                    unscale(last.to_point()).cast<double>());
                const double higher_z = higher_owner_z + 0.5 * (
                    unscale_(first.z()) + unscale_(last.z()));
                for (const OwnedWallPath &lower : feature_walls) {
                    if (lower.path->inset_idx != higher.path->inset_idx ||
                        lower.path->nonplanar_feature_course != previous_course)
                        continue;
                    const double lower_owner_z = layers[lower.owner]->print_z;
                    for (size_t lower_idx = 1;
                         lower_idx < lower.path->polyline.points.size();
                         ++lower_idx) {
                        const Point3 &a3 =
                            lower.path->polyline.points[lower_idx - 1];
                        const Point3 &b3 =
                            lower.path->polyline.points[lower_idx];
                        const Vec2d a =
                            unscale(a3.to_point()).cast<double>();
                        const Vec2d delta =
                            unscale(b3.to_point()).cast<double>() - a;
                        const double squared = delta.squaredNorm();
                        const double ratio = squared <= EPSILON ? 0. :
                            std::clamp((sample - a).dot(delta) / squared,
                                       0., 1.);
                        if ((sample - (a + ratio * delta)).norm() > 0.05)
                            continue;
                        const double lower_z = lower_owner_z +
                            unscale_(a3.z()) + ratio *
                                unscale_(b3.z() - a3.z());
                        const double clearance = 0.5 * std::max(
                            double(lower.path->height),
                            double(higher.path->height));
                        required_shift = std::max(
                            required_shift,
                            lower_z + clearance - higher_z);
                    }
                }
            }
            if (required_shift <= EPSILON)
                continue;
            offset_entity_z(*higher.continuous_wall, coord_t(std::ceil(scale_(required_shift))));
            normalize_fixed_feature_transition(*higher.continuous_wall);
        }
        // Clearance edits may leave a supported feature wall with a larger
        // local bead height than its original planar flow describes. Recover
        // only small, measured gaps directly over preceding feature walls.
        // Do not bridge lateral gaps, move the geometry, or exceed the nozzle
        // and profile height limits. Build all replacements before modifying
        // containers so support pointers remain valid throughout measurement.
        double maximum_height = 0.8 * nozzle_diameter;
        for (double diameter : print_config.nozzle_diameter.values)
            if (diameter > 0.)
                maximum_height = std::min(maximum_height, 0.8 * diameter);
        for (double height : print_config.max_layer_height.values)
            if (height > 0.)
                maximum_height = std::min(maximum_height, height);
        std::map<const ExtrusionPath *, ExtrusionPaths> height_replacements;
        for (const OwnedWallPath &wall : feature_walls) {
            const Points3 &points = wall.path->polyline.points;
            std::vector<const OwnedWallPath *> supports;
            for (const OwnedWallPath &lower : feature_walls)
                if (lower.path->nonplanar_feature_course < wall.path->nonplanar_feature_course &&
                    lower.path->inset_idx == wall.path->inset_idx &&
                    bounds_overlap_with_tolerance(lower.bounds, wall.bounds, scale_(0.05)))
                    supports.push_back(&lower);
            // The first deformed course rests on ordinary sliced walls, not
            // another feature course. Only already printed owners may support
            // it; a future planar wall must never justify extra material here.
            if (wall.path->nonplanar_feature_course == 1)
                for (const OwnedWallPath &lower : foundation_walls)
                    // Ordinary wall generators store the inset on the parent
                    // loop/multipath; its child path may still have inset 0.
                    if (lower.owner < wall.owner &&
                        lower.continuous_wall->inset_idx == wall.path->inset_idx &&
                        bounds_overlap_with_tolerance(lower.bounds, wall.bounds, scale_(0.05)))
                        supports.push_back(&lower);
            std::vector<double> heights(points.size() - 1, wall.path->height);
            bool changed = false;
            for (size_t idx = 1; idx < points.size(); ++idx) {
                double minimum_gap = std::numeric_limits<double>::infinity();
                double maximum_gap = 0.;
                for (double t : {0., 0.5, 1.}) {
                    const Vec3d position = ((1. - t) * points[idx - 1].cast<double>() +
                                           t * points[idx].cast<double>()) * SCALING_FACTOR;
                    const double z = layers[wall.owner]->print_z + position.z();
                    double support_z = -std::numeric_limits<double>::infinity();
                    for (const OwnedWallPath *support : supports) {
                        const OwnedWallPath &lower = *support;
                        const Points3 &below = lower.path->polyline.points;
                        for (size_t j = 1; j < below.size(); ++j) {
                            const Vec2d start = unscale(below[j - 1].to_point());
                            const Vec2d delta = unscale(below[j].to_point()) - start;
                            const double ratio = delta.squaredNorm() <= EPSILON ? 0. :
                                std::clamp((position.head<2>() - start).dot(delta) / delta.squaredNorm(), 0., 1.);
                            if ((position.head<2>() - start - ratio * delta).norm() > 0.05)
                                continue;
                            const double candidate_z = layers[lower.owner]->print_z + unscale_(below[j - 1].z()) +
                                ratio * unscale_(below[j].z() - below[j - 1].z());
                            if (candidate_z < z - EPSILON)
                                support_z = std::max(support_z, candidate_z);
                        }
                    }
                    const double gap = z - support_z;
                    minimum_gap = std::min(minimum_gap, gap);
                    maximum_gap = std::max(maximum_gap, gap);
                }
                if (minimum_gap > double(wall.path->height) + 0.01 &&
                    maximum_gap <= std::min(maximum_height, double(wall.path->width)) &&
                    maximum_gap - minimum_gap <= 0.01) {
                    heights[idx - 1] = 0.5 * (minimum_gap + maximum_gap);
                    changed = true;
                }
            }
            if (changed)
                height_replacements.emplace(wall.path, split_nonplanar_path_by_height(*wall.path, heights));
        }
        const auto replace_paths = [&](ExtrusionPaths &paths) {
            ExtrusionPaths replaced;
            for (ExtrusionPath &path : paths) {
                auto found = height_replacements.find(&path);
                if (found == height_replacements.end())
                    replaced.emplace_back(std::move(path));
                else
                    append(replaced, std::move(found->second));
            }
            paths = std::move(replaced);
        };
        const auto replace_heights = [&](auto &&self, ExtrusionEntity *&entity) -> void {
            if (auto *path = dynamic_cast<ExtrusionPath *>(entity)) {
                auto found = height_replacements.find(path);
                if (found == height_replacements.end())
                    return;
                auto *replacement = new ExtrusionMultiPath(*path);
                replacement->inset_idx = path->inset_idx;
                replacement->paths = std::move(found->second);
                delete entity;
                entity = replacement;
            } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(entity))
                replace_paths(loop->paths);
            else if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(entity))
                replace_paths(multipath->paths);
            else if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(entity))
                for (ExtrusionEntity *&child : collection->entities)
                    self(self, child);
        };
        if (!height_replacements.empty())
            for (Layer *layer : layers) {
                LayerRegion &region = *layer->regions()[region_idx];
                for (ExtrusionEntity *&entity : region.perimeters.entities)
                    replace_heights(replace_heights, entity);
                for (ExtrusionEntity *&entity : region.fills.entities)
                    replace_heights(replace_heights, entity);
            }
    }
    if (progress)
        progress({1., NonplanarProgressStage::Commit, 1, 1});
}

std::vector<uint8_t> select_nonplanar_surface_facets(
    const sla::IndexedMesh &mesh, double maximum_angle_degrees,
    double minimum_vertical_span, double minimum_projected_area)
{
    const indexed_triangle_set &its = *mesh.get_triangle_mesh();
    std::vector<uint8_t> candidates(its.indices.size(), 0);
    const double minimum_normal_z = std::cos(Geometry::deg2rad(
        std::clamp(maximum_angle_degrees, 0., 90.)));
    // Horizontal faces already slice perfectly and must remain ordinary top
    // surfaces. Including them connects unrelated sloped patches through a
    // deck or floor and was the source of the false Benchy selection.
    // Treat tiny modeled/coplanarity slopes as flat too. They slice cleanly as
    // ordinary top surfaces and otherwise let curved rim facets leak onto a
    // nominally flat deck. Two degrees remains well below the user-facing
    // maximum surface angle and the smallest useful non-planar slope.
    constexpr double maximum_flat_normal_z = 0.9993908270190958; // cos(2 degrees)
    std::vector<uint8_t> upper_flat_faces(its.indices.size(), 0);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, its.indices.size()),
        [&](const tbb::blocked_range<size_t> &range) {
            for (size_t face = range.begin(); face < range.end(); ++face) {
                const double normal_z = mesh.normal_by_face_id(int(face)).z();
                candidates[face] = normal_z + 1e-6 >= minimum_normal_z &&
                                   normal_z < maximum_flat_normal_z;
                upper_flat_faces[face] = normal_z >= maximum_flat_normal_z;
                if (!candidates[face] && !upper_flat_faces[face])
                    continue;

                // A finishing skin is valid only on the uppermost model
                // surface at its XY position. Each worker writes a distinct
                // byte and IndexedMesh ray queries are read-only.
                const Vec3i32 &indices = its.indices[face];
                const Vec3d centroid = (
                    its.vertices[size_t(indices[0])].cast<double>() +
                    its.vertices[size_t(indices[1])].cast<double>() +
                    its.vertices[size_t(indices[2])].cast<double>()) / 3.;
                const sla::IndexedMesh::hit_result above = mesh.query_ray_hit(
                    centroid + Vec3d(0., 0., 0.01), Vec3d::UnitZ());
                if (above.is_hit() && above.position().z() > centroid.z() + 0.02) {
                    candidates[face] = 0;
                    upper_flat_faces[face] = 0;
                }
            }
        });

    const std::vector<Vec3i32> neighbors = its_face_neighbors_par(its);

    // A shallow roof may contain a narrow horizontal crown between its two
    // sloped sides. Leaving that crown as ordinary top fill splits every
    // finishing raster into unrelated pieces and leaves a visible seam (the
    // Benchy cabin roof is a common example). Admit only bounded flat
    // connector components that touch eligible shallow faces on at least two
    // sides. Broad decks and other genuinely planar tops remain conventional.
    constexpr double maximum_flat_connector_width = 5.;
    std::vector<uint8_t> visited_flat(its.indices.size(), 0);
    for (size_t seed = 0; seed < upper_flat_faces.size(); ++seed) {
        if (!upper_flat_faces[seed] || visited_flat[seed])
            continue;
        std::vector<size_t> flat_component;
        std::vector<size_t> adjacent_slopes;
        std::queue<size_t> flat_pending;
        BoundingBoxf bbox;
        flat_pending.push(seed);
        visited_flat[seed] = 1;
        while (!flat_pending.empty()) {
            const size_t face = flat_pending.front();
            flat_pending.pop();
            flat_component.push_back(face);
            const Vec3i32 &indices = its.indices[face];
            for (int vertex = 0; vertex < 3; ++vertex)
                bbox.merge(its.vertices[size_t(indices[vertex])].head<2>().cast<double>());
            for (int edge = 0; edge < 3; ++edge) {
                const int neighbor = neighbors[face][edge];
                if (neighbor < 0)
                    continue;
                if (candidates[size_t(neighbor)])
                    adjacent_slopes.push_back(size_t(neighbor));
                else if (upper_flat_faces[size_t(neighbor)] &&
                         !visited_flat[size_t(neighbor)]) {
                    visited_flat[size_t(neighbor)] = 1;
                    flat_pending.push(size_t(neighbor));
                }
            }
        }
        sort_remove_duplicates(adjacent_slopes);
        const Vec2d size = bbox.size();
        if (adjacent_slopes.size() >= 2 &&
            std::min(size.x(), size.y()) <= maximum_flat_connector_width + EPSILON)
            for (size_t face : flat_component)
                candidates[face] = 1;
    }

    std::vector<uint8_t> selected(its.indices.size(), 0);
    std::vector<uint8_t> visited(its.indices.size(), 0);
    std::vector<size_t> component;
    std::queue<size_t> pending;
    for (size_t seed = 0; seed < candidates.size(); ++seed) {
        if (!candidates[seed] || visited[seed])
            continue;
        component.clear();
        pending.push(seed);
        visited[seed] = 1;
        double projected_area = 0.;
        double minimum_z = std::numeric_limits<double>::max();
        double maximum_z = std::numeric_limits<double>::lowest();
        while (!pending.empty()) {
            const size_t face = pending.front();
            pending.pop();
            component.push_back(face);
            const Vec3i32 &indices = its.indices[face];
            const Vec3d a = its.vertices[size_t(indices[0])].cast<double>();
            const Vec3d b = its.vertices[size_t(indices[1])].cast<double>();
            const Vec3d c = its.vertices[size_t(indices[2])].cast<double>();
            projected_area += 0.5 * std::abs(
                (b.x() - a.x()) * (c.y() - a.y()) -
                (b.y() - a.y()) * (c.x() - a.x()));
            minimum_z = std::min({minimum_z, a.z(), b.z(), c.z()});
            maximum_z = std::max({maximum_z, a.z(), b.z(), c.z()});
            for (int edge = 0; edge < 3; ++edge) {
                const int neighbor = neighbors[face][edge];
                if (neighbor >= 0 && candidates[size_t(neighbor)] && !visited[size_t(neighbor)]) {
                    visited[size_t(neighbor)] = 1;
                    pending.push(size_t(neighbor));
                }
            }
        }
        // Eligibility is ultimately gated by exposed top-solid toolpaths, so
        // height within the model is not a safe proxy for whether a surface is
        // a top. A sloped deck may legitimately be low while a bottom or side
        // facet cannot create eligible top-solid coverage.
        if (projected_area + EPSILON < minimum_projected_area ||
            maximum_z - minimum_z + EPSILON < minimum_vertical_span)
            continue;
        for (size_t face : component)
            selected[face] = 1;
    }
    // An empty vector is the fast-path sentinel used by every layer. Keeping
    // a mesh-sized all-zero mask made the slicer ray-test every top path on
    // every layer even after geometry eligibility rejected the whole model.
    if (std::none_of(selected.begin(), selected.end(), [](uint8_t value) {
            return value != 0;
        }))
        selected.clear();
    return selected;
}

std::vector<std::vector<ExPolygons>> project_nonplanar_surface_facets(
    const sla::IndexedMesh &mesh,
    const std::vector<std::vector<uint8_t>> &selected_facets,
    const std::vector<std::pair<double, double>> &layer_z_ranges)
{
    std::vector<std::vector<ExPolygons>> projections(
        selected_facets.size(), std::vector<ExPolygons>(layer_z_ranges.size()));
    if (layer_z_ranges.empty())
        return projections;

    const indexed_triangle_set &its = *mesh.get_triangle_mesh();
    for (size_t region_idx = 0; region_idx < selected_facets.size(); ++region_idx) {
        if (selected_facets[region_idx].empty())
            continue;
        std::vector<Polygons> layer_triangles(layer_z_ranges.size());
        const size_t face_count = std::min(
            selected_facets[region_idx].size(), its.indices.size());
        for (size_t face_idx = 0; face_idx < face_count; ++face_idx) {
            if (!selected_facets[region_idx][face_idx])
                continue;
            const Vec3i32 &indices = its.indices[face_idx];
            const Vec3f &a = its.vertices[size_t(indices[0])];
            const Vec3f &b = its.vertices[size_t(indices[1])];
            const Vec3f &c = its.vertices[size_t(indices[2])];
            const double minimum_z = std::min({double(a.z()), double(b.z()), double(c.z())});
            const double maximum_z = std::max({double(a.z()), double(b.z()), double(c.z())});
            auto range_it = std::lower_bound(
                layer_z_ranges.begin(), layer_z_ranges.end(), minimum_z,
                [](const std::pair<double, double> &range, double z) {
                    return range.second + EPSILON < z;
                });
            if (range_it == layer_z_ranges.end())
                continue;
            Polygon triangle({
                Point(scale_(a.x()), scale_(a.y())),
                Point(scale_(b.x()), scale_(b.y())),
                Point(scale_(c.x()), scale_(c.y()))});
            if (triangle.area() < 0)
                triangle.reverse();
            if (!triangle.is_valid())
                continue;
            for (size_t layer_idx = size_t(std::distance(
                     layer_z_ranges.begin(), range_it));
                 layer_idx < layer_z_ranges.size() &&
                     layer_z_ranges[layer_idx].first <= maximum_z + EPSILON;
                 ++layer_idx)
                layer_triangles[layer_idx].emplace_back(triangle);
        }
        for (size_t layer_idx = 0; layer_idx < layer_triangles.size(); ++layer_idx) {
            if (layer_triangles[layer_idx].empty())
                continue;
            // The source perimeter lies on the interface between an upper
            // surface and its side wall, while the selected facet projection
            // ends exactly at that interface. A one-millimetre query halo is
            // deliberately only a performance filter: the facet mask still
            // decides eligibility, so it cannot select a neighboring surface.
            projections[region_idx][layer_idx] = offset_ex(
                union_ex(layer_triangles[layer_idx]), scale_(1.));
        }
    }
    return projections;
}

namespace {

bool project_support_contact_path(
    ExtrusionPath &path, const SupportLayer &support_layer,
    const sla::IndexedMesh &mesh, double contact_gap,
    double maximum_angle_degrees, double resolution)
{
    if ((path.role() != erSupportMaterialInterface && path.role() != erIroning) ||
        path.polyline.points.size() < 2)
        return false;

    const double owner_z = support_layer.print_z;
    // Start below the deposited support surface. With a soluble zero contact
    // gap, starting above owner_z would put the ray inside the model and skip
    // the underside in favor of a later upward-facing exit triangle.
    const double source_z = mesh.ground_level() + owner_z - 0.01;
    const double max_angle = maximum_angle_degrees * M_PI / 180.;
    const double max_rise = std::max(double(support_layer.height),
                                     double(path.height)) + 0.02;
    const Points3 original = path.polyline.points;
    Points3 projected;
    std::vector<double> distance;
    std::vector<bool> has_underside;

    auto underside_at = [&](const Vec2d &xy) -> std::optional<double> {
        const auto hits = mesh.query_ray_hits(
            Vec3d(xy.x(), xy.y(), source_z), Vec3d::UnitZ());
        for (const sla::IndexedMesh::hit_result &hit : hits) {
            if (!hit.is_hit() || hit.normal().z() >= -EPSILON)
                continue;
            const double angle = std::acos(std::clamp(-hit.normal().z(), 0., 1.));
            if (angle > max_angle + EPSILON)
                return std::nullopt;
            const double target_z = hit.position().z() - mesh.ground_level() - contact_gap;
            const double offset = target_z - owner_z;
            // The interface supplies the inclined foundation for its ironing
            // pass. A larger rise would leave both paths suspended above the
            // preceding support course, so reject the complete path.
            if (offset < -0.02 || offset > max_rise)
                return std::nullopt;
            return offset;
        }
        return std::nullopt;
    };

    resolution = std::max(0.02, resolution);
    for (size_t segment_idx = 0; segment_idx + 1 < original.size(); ++segment_idx) {
        const Vec2d from = unscale(original[segment_idx].to_point());
        const Vec2d to = unscale(original[segment_idx + 1].to_point());
        const Vec2d delta = to - from;
        const int segments = std::max(1, int(std::ceil(delta.norm() / resolution)));
        for (int idx = segment_idx == 0 ? 0 : 1; idx <= segments; ++idx) {
            const Vec2d xy = from + double(idx) / double(segments) * delta;
            const std::optional<double> offset = underside_at(xy);
            projected.emplace_back(scale_(xy.x()), scale_(xy.y()),
                                   scale_(offset.value_or(0.)));
            has_underside.emplace_back(offset.has_value());
            if (distance.empty())
                distance.emplace_back(0.);
            else
                distance.emplace_back(distance.back() +
                    (unscale(projected.back().to_point()) -
                     unscale(projected[projected.size() - 2].to_point())).norm());
        }
    }
    if (projected.size() < 2 || distance.back() <= EPSILON)
        return false;

    // Support-interface lines may extend beyond the model footprint. Keep
    // those spans planar and contour each contiguous span that actually has
    // an underside above it. Entry and exit use the same direction-independent
    // smoothstep, so planar-to-surface is exactly the reverse of
    // surface-to-planar without splitting the original extrusion path.
    for (size_t begin = 0; begin < projected.size();) {
        while (begin < projected.size() && !has_underside[begin])
            ++begin;
        if (begin == projected.size())
            break;
        size_t end = begin;
        while (end + 1 < projected.size() && has_underside[end + 1])
            ++end;
        const double span_length = distance[end] - distance[begin];
        const double transition_length = std::min(
            0.5 * span_length, std::max(2., 4. * double(path.width)));
        for (size_t idx = begin; idx <= end; ++idx) {
            const double boundary = std::min(distance[idx] - distance[begin],
                                             distance[end] - distance[idx]);
            const double t = transition_length <= EPSILON ? 0. :
                std::clamp(boundary / transition_length, 0., 1.);
            const double smooth = t * t * (3. - 2. * t);
            projected[idx].z() = coord_t(std::lround(
                double(projected[idx].z()) * smooth));
        }
        begin = end + 1;
    }

    bool changed = false;
    for (size_t idx = 1; idx < projected.size(); ++idx) {
        const Vec3d from(unscale_(projected[idx - 1].x()),
                         unscale_(projected[idx - 1].y()),
                         unscale_(projected[idx - 1].z()));
        const Vec3d to(unscale_(projected[idx].x()),
                       unscale_(projected[idx].y()),
                       unscale_(projected[idx].z()));
        const double xy = (to.head<2>() - from.head<2>()).norm();
        if (xy > EPSILON && std::atan(std::abs(to.z() - from.z()) / xy) >
                max_angle + EPSILON)
            return false;
        changed |= std::abs(to.z()) > EPSILON;
    }
    if (!changed)
        return false;

    Polyline3 polyline;
    polyline.points = std::move(projected);
    polyline.simplify(scale_(std::min(0.005, resolution * 0.025)));
    if (polyline.points.size() < 2)
        return false;
    path.polyline = std::move(polyline);
    path.polyline.fitting_result.clear();
    path.z_contoured = true;
    path.nonplanar_surface = true;
    path.nonplanar_schedule_owned = false;
    path.nonplanar_clearance_validated = false;
    return true;
}

size_t project_support_contact_entity(
    ExtrusionEntity &entity, const SupportLayer &support_layer,
    ExtrusionRole role,
    const sla::IndexedMesh &mesh, double contact_gap,
    double maximum_angle_degrees, double resolution,
    const SequentialGantryGeometry &gantry,
    const ToolheadClearanceFallback &fallback,
    const std::vector<PrintedToolpathObstacle> &obstacles)
{
    auto project_path = [&](ExtrusionPath &path) -> size_t {
        if (path.role() != role)
            return 0;
        ExtrusionPath candidate(path);
        if (!project_support_contact_path(candidate, support_layer, mesh,
                                          contact_gap, maximum_angle_degrees,
                                          resolution))
            return 0;
        if (!generated_path_is_clear(candidate, support_layer.print_z,
                                     gantry, fallback, obstacles))
            return 0;
        path = std::move(candidate);
        return 1;
    };
    if (auto *path = dynamic_cast<ExtrusionPath *>(&entity))
        return project_path(*path);
    size_t count = 0;
    if (auto *multipath = dynamic_cast<ExtrusionMultiPath *>(&entity)) {
        for (ExtrusionPath &path : multipath->paths)
            count += project_path(path);
    } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
        for (ExtrusionPath &path : loop->paths)
            count += project_path(path);
    } else if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
        for (ExtrusionEntity *child : collection->entities)
            count += project_support_contact_entity(
                *child, support_layer, role, mesh, contact_gap,
                maximum_angle_degrees, resolution, gantry, fallback, obstacles);
    }
    return count;
}

} // namespace

size_t project_nonplanar_support_interface(
    const std::vector<Layer *> &model_layers,
    const std::vector<SupportLayer *> &support_layers,
    const sla::IndexedMesh &mesh, double contact_gap,
    double maximum_angle_degrees, double resolution,
    const std::function<void(size_t, size_t)> &progress)
{
    if (model_layers.empty() || support_layers.empty())
        return 0;
    const PrintConfig &print_config = model_layers.front()->object()->print()->config();
    SequentialGantryGeometry gantry = load_sequential_gantry_geometry(print_config);
    const double nozzle_diameter = print_config.nozzle_diameter.values.empty() ? 0.4 :
        *std::max_element(print_config.nozzle_diameter.values.begin(),
                          print_config.nozzle_diameter.values.end());
    gantry.set_nozzle_tip_diameter(nozzle_diameter);
    const ToolheadClearanceFallback fallback = make_nonplanar_clearance_cone(print_config);

    size_t projected = 0;
    size_t completed = 0;
    for (SupportLayer *support_layer : support_layers) {
        if (support_layer != nullptr && support_layer->has_top_interface_contact) {
            size_t printed_model_layer_count = 0;
            while (printed_model_layer_count < model_layers.size() &&
                   model_layers[printed_model_layer_count]->print_z <=
                       support_layer->print_z + EPSILON)
                ++printed_model_layer_count;
            const std::vector<PrintedToolpathObstacle> obstacles =
                printed_model_layer_count == 0 ?
                    std::vector<PrintedToolpathObstacle>{} :
                    printed_model_obstacles(model_layers,
                                             printed_model_layer_count - 1);

            size_t projected_interface = 0;
            for (ExtrusionEntity *entity : support_layer->support_fills.entities)
                projected_interface += project_support_contact_entity(
                    *entity, *support_layer, erSupportMaterialInterface,
                    mesh, contact_gap, maximum_angle_degrees, resolution,
                    gantry, fallback, obstacles);
            projected += projected_interface;

            // Ironing is a finishing pass, not a replacement for its support
            // foundation. Only contour it after at least one interface path
            // on this contact layer survived the same surface and toolhead
            // validation. This prevents an independently projected ironing
            // line from being deposited in mid-air when the interface had to
            // fall back to its conventional planar geometry.
            if (projected_interface > 0)
                for (ExtrusionEntity *entity : support_layer->support_fills.entities)
                    projected += project_support_contact_entity(
                        *entity, *support_layer, erIroning, mesh, contact_gap,
                        maximum_angle_degrees, resolution, gantry, fallback,
                        obstacles);
        }
        if (progress)
            progress(++completed, support_layers.size());
    }
    return projected;
}

} // namespace Slic3r
