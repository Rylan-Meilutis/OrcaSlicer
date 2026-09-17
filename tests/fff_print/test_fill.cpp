#include <catch2/catch_all.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/AABBTreeLines.hpp"
#include "libslic3r/Fill/Fill.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/Geometry.hpp"
#include "libslic3r/Geometry/ConvexHull.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/MeshBoolean.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/SVG.hpp"
#include "libslic3r/Tesselate.hpp"
#include "libslic3r/libslic3r.h"

#include "test_helpers.hpp"
#include "test_utils.hpp"

using namespace Slic3r;

bool test_if_solid_surface_filled(const ExPolygon& expolygon, double flow_spacing, double angle = 0, double density = 1.0);

static bool print_has_arc_overhang(const Print &print)
{
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions())
                for (const ExtrusionEntity *entity : region->fills.flatten().entities)
                    if (is_arc_fill(entity->role()))
                        return true;
    return false;
}

static std::optional<ExtrusionPath> first_arc_overhang_path(const Print &print)
{
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions())
                for (const ExtrusionEntity *entity : region->fills.flatten().entities)
                    if (const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
                        path != nullptr && is_arc_fill(path->role()) && path->mm3_per_mm > 0.)
                        return *path;
    return std::nullopt;
}

static std::vector<ExtrusionPath> nonplanar_top_paths(const Print &print)
{
    std::vector<ExtrusionPath> paths;
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions())
                for (const ExtrusionEntity *entity : region->fills.flatten().entities)
                    if (const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
                        path != nullptr && path->nonplanar_surface && path->role() == erTopSolidInfill)
                        paths.push_back(*path);
    return paths;
}

static std::vector<ExtrusionPath> nonplanar_perimeter_paths(const Print &print)
{
    std::vector<ExtrusionPath> paths;
    const auto collect = [&paths](const ExtrusionEntityCollection &collection) {
        for (const ExtrusionEntity *entity : collection.flatten().entities)
            if (const auto *path = dynamic_cast<const ExtrusionPath *>(entity)) {
                if (path->nonplanar_surface && is_perimeter(path->role()))
                    paths.push_back(*path);
            } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(entity)) {
                for (const ExtrusionPath &path : loop->paths)
                    if (path.nonplanar_surface && is_perimeter(path.role()))
                        paths.push_back(path);
            } else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(entity)) {
                for (const ExtrusionPath &path : multipath->paths)
                    if (path.nonplanar_surface && is_perimeter(path.role()))
                        paths.push_back(path);
            }
    };
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions()) {
                collect(region->perimeters);
                collect(region->fills);
            }
    return paths;
}

static bool has_connected_nonplanar_wall_loop(const Print &print)
{
    const auto connected = [](const ExtrusionPaths &paths) {
        if (paths.empty() || paths.front().first_point3() != paths.back().last_point3())
            return false;
        bool has_nonplanar = false;
        for (size_t idx = 0; idx < paths.size(); ++idx) {
            if (!is_perimeter(paths[idx].role()))
                return false;
            has_nonplanar |= paths[idx].nonplanar_surface;
            if (idx > 0 && paths[idx - 1].last_point3() != paths[idx].first_point3())
                return false;
        }
        return has_nonplanar;
    };
    const auto inspect = [&](auto &&self, const ExtrusionEntity &entity) -> bool {
        if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity))
            return connected(loop->paths);
        if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity))
            return connected(multipath->paths);
        if (const auto *collection =
                dynamic_cast<const ExtrusionEntityCollection *>(&entity))
            return std::any_of(collection->entities.begin(), collection->entities.end(),
                [&](const ExtrusionEntity *child) { return self(self, *child); });
        return false;
    };
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions())
                if (std::any_of(region->perimeters.entities.begin(),
                        region->perimeters.entities.end(),
                        [&](const ExtrusionEntity *entity) {
                            return inspect(inspect, *entity);
                        }))
                    return true;
    return false;
}

static void append_nonplanar_wall_insets_in_order(
    const ExtrusionEntity &entity, std::vector<unsigned int> &insets)
{
    const auto contains_nonplanar = [](auto &&self,
                                       const ExtrusionEntity &candidate) -> bool {
        if (const auto *path = dynamic_cast<const ExtrusionPath *>(&candidate))
            return path->nonplanar_surface && is_perimeter(path->role());
        if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&candidate))
            return std::any_of(loop->paths.begin(), loop->paths.end(),
                [](const ExtrusionPath &path) {
                    return path.nonplanar_surface && is_perimeter(path.role());
                });
        if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&candidate))
            return std::any_of(multipath->paths.begin(), multipath->paths.end(),
                [](const ExtrusionPath &path) {
                    return path.nonplanar_surface && is_perimeter(path.role());
                });
        if (const auto *collection =
                dynamic_cast<const ExtrusionEntityCollection *>(&candidate))
            return std::any_of(collection->entities.begin(), collection->entities.end(),
                [&](const ExtrusionEntity *child) { return self(self, *child); });
        return false;
    };
    if (!contains_nonplanar(contains_nonplanar, entity))
        return;
    if (const auto *collection =
                   dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        for (const ExtrusionEntity *child : collection->entities) {
            if (child->inset_idx >= 0 && is_perimeter(child->role()))
                insets.emplace_back(child->inset_idx);
            else
                append_nonplanar_wall_insets_in_order(*child, insets);
        }
    } else if (entity.inset_idx >= 0 && is_perimeter(entity.role()))
        insets.emplace_back(entity.inset_idx);
}

static std::vector<unsigned int> nonplanar_wall_inset_order(const Print &print)
{
    std::vector<unsigned int> insets;
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions()) {
                for (const ExtrusionEntity *entity : region->perimeters.entities)
                    append_nonplanar_wall_insets_in_order(*entity, insets);
                for (const ExtrusionEntity *entity : region->fills.entities)
                    append_nonplanar_wall_insets_in_order(*entity, insets);
            }
    insets.erase(std::remove_if(insets.begin(), insets.end(),
        [](unsigned int inset) { return inset > 2; }), insets.end());
    // Clipping a wall at the replacement boundary may split one inset into
    // several consecutive path fragments. They are still one scheduled wall,
    // so compare changes in inset rather than mistaking those fragments for
    // three separately ordered walls.
    insets.erase(std::unique(insets.begin(), insets.end()), insets.end());
    // The first complete three-wall course is sufficient to verify the
    // configured order; subsequent courses repeat that sequence.
    if (insets.size() > 3)
        insets.resize(3);
    return insets;
}

static size_t perimeter_path_count(const Print &print, ExtrusionRole role)
{
    size_t count = 0;
    auto count_paths = [role, &count](const ExtrusionEntity &entity) {
        if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity)) {
            count += path->role() == role;
        } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
            count += std::count_if(loop->paths.begin(), loop->paths.end(),
                [role](const ExtrusionPath &path) { return path.role() == role; });
        } else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
            count += std::count_if(multipath->paths.begin(), multipath->paths.end(),
                [role](const ExtrusionPath &path) { return path.role() == role; });
        }
    };
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions())
                for (const ExtrusionEntity *entity : region->perimeters.flatten().entities)
                    count_paths(*entity);
    return count;
}

static size_t staggered_perimeter_path_count(const Print &print)
{
    size_t count = 0;
    auto count_paths = [&count](const ExtrusionEntity &entity) {
        if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity)) {
            count += path->staggered_perimeter;
        } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
            count += std::count_if(loop->paths.begin(), loop->paths.end(),
                [](const ExtrusionPath &path) { return path.staggered_perimeter; });
        } else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
            count += std::count_if(multipath->paths.begin(), multipath->paths.end(),
                [](const ExtrusionPath &path) { return path.staggered_perimeter; });
        }
    };
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions())
                for (const ExtrusionEntity *entity : region->perimeters.flatten().entities)
                    count_paths(*entity);
    return count;
}

static size_t interlocking_perimeter_path_count(const Print &print)
{
    size_t count = 0;
    const auto count_paths = [&count](auto &&self,
                                      const ExtrusionEntity &entity) -> void {
        if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity)) {
            count += path->nonplanar_interlocking_wall;
        } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
            for (const ExtrusionPath &path : loop->paths)
                count += path.nonplanar_interlocking_wall;
        } else if (const auto *multipath =
                       dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
            for (const ExtrusionPath &path : multipath->paths)
                count += path.nonplanar_interlocking_wall;
        } else if (const auto *collection =
                       dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
            for (const ExtrusionEntity *child : collection->entities)
                self(self, *child);
        }
    };
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions())
                for (const ExtrusionEntity *entity : region->perimeters.entities)
                    count_paths(count_paths, *entity);
    return count;
}

static bool has_nonplanar_path_with_role(const Print &print, ExtrusionRole role)
{
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions())
                for (const ExtrusionEntity *entity : region->fills.flatten().entities)
                    if (const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
                        path != nullptr && path->nonplanar_surface && path->role() == role)
                        return true;
    return false;
}

static bool has_planar_path_inside_nonplanar_top_coverage(const Print &print)
{
    for (const PrintObject *object : print.objects()) {
        for (const Layer *layer : object->layers()) {
            Polygons coverage;
            for (const LayerRegion *region : layer->regions()) {
                for (const ExtrusionEntity *entity : region->fills.flatten().entities) {
                    if (const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
                        path != nullptr && path->nonplanar_surface && path->role() == erTopSolidInfill)
                        path->polygons_covered_by_spacing(coverage, float(scale_(0.02)));
                }
            }
            if (coverage.empty())
                continue;
            const ExPolygons interior = offset_ex(union_ex(coverage), -scale_(0.01));
            const auto overlaps = [&interior](const ExtrusionPath &path) {
                if (path.nonplanar_surface || path.role() != erPerimeter ||
                    path.inset_idx == 0 || path.polyline.points.empty())
                    return false;
                coord_t min_x = path.polyline.points.front().x();
                coord_t max_x = min_x;
                coord_t min_y = path.polyline.points.front().y();
                coord_t max_y = min_y;
                for (const Point3 &point : path.polyline.points) {
                    min_x = std::min(min_x, point.x());
                    max_x = std::max(max_x, point.x());
                    min_y = std::min(min_y, point.y());
                    max_y = std::max(max_y, point.y());
                }
                // Ignore legitimate compact features (for example the
                // chimney bore). The regression is the large inner-wall box
                // spanning the accepted roof patch.
                if (unscale_(max_x - min_x) < 12. ||
                    unscale_(max_y - min_y) < 8.)
                    return false;
                return !intersection_pl(
                    Polylines{path.polyline.to_polyline()}, interior).empty();
            };
            for (const LayerRegion *region : layer->regions()) {
                for (const ExtrusionEntityCollection *collection : {&region->fills, &region->perimeters}) {
                    for (const ExtrusionEntity *entity : collection->flatten().entities) {
                        if (const auto *path = dynamic_cast<const ExtrusionPath *>(entity)) {
                            if (overlaps(*path))
                                return true;
                        } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(entity)) {
                            if (std::any_of(loop->paths.begin(), loop->paths.end(), overlaps))
                                return true;
                        } else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(entity)) {
                            if (std::any_of(multipath->paths.begin(), multipath->paths.end(), overlaps))
                                return true;
                        }
                    }
                }
            }
        }
    }
    return false;
}

static TriangleMesh shallow_top_wedge(double angle_degrees)
{
    TriangleMesh cube = Slic3r::Test::cube(20.);
    std::vector<Vec3f> vertices = cube.its.vertices;
    const std::vector<Vec3i32> faces = cube.its.indices;
    const float top_z = std::max_element(vertices.begin(), vertices.end(),
        [](const Vec3f &lhs, const Vec3f &rhs) { return lhs.z() < rhs.z(); })->z();
    const float slope = float(std::tan(angle_degrees * M_PI / 180.));
    for (Vec3f &vertex : vertices)
        if (std::abs(vertex.z() - top_z) < float(EPSILON))
            vertex.z() += slope * vertex.x();
    return TriangleMesh(std::move(vertices), std::move(faces));
}

template<class BottomZ, class TopZ>
static TriangleMesh extrude_sloped_footprint(const ExPolygon &footprint,
                                              BottomZ bottom_z, TopZ top_z)
{
    std::vector<Vec3f> vertices;
    std::vector<Vec3i32> faces;
    const auto append_triangle_soup = [&](const std::vector<Vec2f> &triangles,
                                          const auto &z_at) {
        REQUIRE(triangles.size() % 3 == 0);
        for (size_t idx = 0; idx < triangles.size(); idx += 3) {
            const int32_t first = int32_t(vertices.size());
            for (size_t point_idx = 0; point_idx < 3; ++point_idx) {
                const Vec2f &point = triangles[idx + point_idx];
                vertices.emplace_back(point.x(), point.y(), float(z_at(point.x())));
            }
            faces.emplace_back(first, first + 1, first + 2);
        }
    };
    append_triangle_soup(triangulate_expolygon_2f(footprint, NORMALS_DOWN), bottom_z);
    append_triangle_soup(triangulate_expolygon_2f(footprint, NORMALS_UP), top_z);

    const auto append_sides = [&](const Polygon &polygon) {
        for (size_t idx = 0; idx < polygon.points.size(); ++idx) {
            const Point &a_scaled = polygon.points[idx];
            const Point &b_scaled = polygon.points[(idx + 1) % polygon.points.size()];
            const Vec2f a(float(unscale_(a_scaled.x())), float(unscale_(a_scaled.y())));
            const Vec2f b(float(unscale_(b_scaled.x())), float(unscale_(b_scaled.y())));
            const int32_t first = int32_t(vertices.size());
            vertices.emplace_back(a.x(), a.y(), float(bottom_z(a.x())));
            vertices.emplace_back(b.x(), b.y(), float(bottom_z(b.x())));
            vertices.emplace_back(b.x(), b.y(), float(top_z(b.x())));
            vertices.emplace_back(a.x(), a.y(), float(top_z(a.x())));
            faces.emplace_back(first, first + 1, first + 2);
            faces.emplace_back(first, first + 2, first + 3);
        }
    };
    append_sides(footprint.contour);
    for (const Polygon &hole : footprint.holes)
        append_sides(hole);
    TriangleMesh mesh(std::move(vertices), std::move(faces));
    its_merge_vertices(mesh.its);
    return mesh;
}

// A compact version of the Benchy cabin/roof topology: a shallow sloped roof
// interrupted by an open, hollow chimney. It catches regressions where the
// non-planar patch closes itself with a new inner-wall box, or consumes either
// side of the chimney wall while replacing the roof's planar top courses.
static TriangleMesh benchy_roof_with_chimney(double angle_degrees)
{
    constexpr double width = 24.;
    constexpr double depth = 20.;
    constexpr double roof_height = 12.;
    constexpr double chimney_x = 7.;
    constexpr double chimney_y = 10.;
    // Deliberately thick enough for three requested walls on both sides of the
    // annulus. A thinner cosmetic tube would legitimately collapse inner
    // insets under Arachne and could not detect a missing-wall regression.
    constexpr double chimney_outer_radius = 5.0;
    constexpr double chimney_inner_radius = 1.5;
    constexpr double chimney_top = 22.;
    const double slope = std::tan(angle_degrees * M_PI / 180.);

    Polygon rectangle({Point::new_scale(0., 0.), Point::new_scale(width, 0.),
                       Point::new_scale(width, depth), Point::new_scale(0., depth)});
    ExPolygon roof_footprint(rectangle);
    TriangleMesh model = extrude_sloped_footprint(
        roof_footprint, [](double) { return 0.; },
        [slope](double x) { return roof_height + slope * x; });

    Polygon chimney_outer = make_circle_num_segments(scale_(chimney_outer_radius), 64);
    chimney_outer.translate(Point::new_scale(chimney_x, chimney_y));
    ExPolygon chimney_footprint(chimney_outer);
    TriangleMesh chimney = extrude_sloped_footprint(
        chimney_footprint,
        [slope](double x) { return roof_height - 0.6 + slope * x; },
        [](double) { return chimney_top; });
    MeshBoolean::cgal::plus(model, chimney);

    // Build a single watertight volume. Merely merging overlapping triangle
    // soups makes the roof/chimney intersection non-manifold and can itself
    // make upper chimney contours disappear, masking the slicing regression
    // this fixture is intended to detect.
    Polygon chimney_bore = make_circle_num_segments(scale_(chimney_inner_radius), 48);
    chimney_bore.translate(Point::new_scale(chimney_x, chimney_y));
    TriangleMesh bore = extrude_sloped_footprint(
        ExPolygon(chimney_bore), [](double) { return -1.; },
        [](double) { return chimney_top + 1.; });
    MeshBoolean::cgal::minus(model, bore);
    return model;
}

static bool arc_paths_have_proper_crossing(const Polylines &paths,
                                           std::string *details = nullptr)
{
    struct PathLine {
        Line line;
        size_t path_idx;
        size_t line_idx;
    };
    std::vector<PathLine> lines;
    for (size_t path_idx = 0; path_idx < paths.size(); ++path_idx)
        for (size_t line_idx = 1; line_idx < paths[path_idx].points.size(); ++line_idx)
            lines.push_back({Line(paths[path_idx].points[line_idx - 1],
                                  paths[path_idx].points[line_idx]),
                             path_idx, line_idx - 1});

    for (size_t first_idx = 0; first_idx < lines.size(); ++first_idx) {
        for (size_t second_idx = first_idx + 1; second_idx < lines.size(); ++second_idx) {
            const PathLine &first = lines[first_idx];
            const PathLine &second = lines[second_idx];
            if (first.path_idx == second.path_idx &&
                second.line_idx <= first.line_idx + 1)
                continue;
            Point intersection;
            if (!first.line.intersection(second.line, &intersection))
                continue;
            const Vec2d p = first.line.a.cast<double>();
            const Vec2d q = second.line.a.cast<double>();
            const Vec2d r = (first.line.b - first.line.a).cast<double>();
            const Vec2d s = (second.line.b - second.line.a).cast<double>();
            const auto cross = [](const Vec2d &lhs, const Vec2d &rhs) {
                return lhs.x() * rhs.y() - lhs.y() * rhs.x();
            };
            const double denominator = cross(r, s);
            if (std::abs(denominator) <= std::numeric_limits<double>::epsilon())
                continue;
            const double first_position = cross(q - p, s) / denominator;
            const double second_position      = cross(q - p, r) / denominator;
            const auto path_endpoint_distance = [&paths, &intersection](size_t path_idx) {
                return std::min((intersection - paths[path_idx].first_point()).cast<double>().norm(),
                                (intersection - paths[path_idx].last_point()).cast<double>().norm());
            };
            const double first_endpoint_distance  = path_endpoint_distance(first.path_idx);
            const double second_endpoint_distance = path_endpoint_distance(second.path_idx);
            // Parent/child arcs intentionally meet at a path endpoint. Arc
            // fitting and three-decimal G-code coordinates can move the
            // centerline intersection away from that endpoint, but it remains
            // an anchored T-junction while it lies within half of the standard
            // 0.45 mm deposited bead. An intersection farther into both
            // complete paths is a real crossing.
            const double endpoint_epsilon = scale_(0.225);
            // A child arc is allowed to terminate on its supporting parent.
            // Only an intersection in the interior of both complete paths
            // means extrusion actually continues across an existing path.
            if (first_endpoint_distance > endpoint_epsilon && second_endpoint_distance > endpoint_epsilon) {
                if (details != nullptr) {
                    std::ostringstream message;
                    message << "paths " << first.path_idx << "/" << second.path_idx << ", segments " << first.line_idx << "/"
                            << second.line_idx << ", intersection " << unscale<double>(intersection.x()) << ","
                            << unscale<double>(intersection.y()) << ", endpoint distances " << unscale<double>(first_endpoint_distance)
                            << "/" << unscale<double>(second_endpoint_distance) << " mm";
                    *details = message.str();
                }
                return true;
            }
        }
    }
    return false;
}

static bool arc_paths_cross_obstacles(const Polylines &arcs,
                                      const Polylines &obstacles,
                                      std::string *details = nullptr)
{
    for (size_t arc_idx = 0; arc_idx < arcs.size(); ++arc_idx) {
        const Polyline &arc = arcs[arc_idx];
        for (size_t arc_line_idx = 1;
             arc_line_idx < arc.points.size(); ++arc_line_idx) {
            const Line arc_line(
                arc.points[arc_line_idx - 1], arc.points[arc_line_idx]);
            for (size_t obstacle_idx = 0;
                 obstacle_idx < obstacles.size(); ++obstacle_idx) {
                const Polyline &obstacle = obstacles[obstacle_idx];
                for (size_t obstacle_line_idx = 1;
                     obstacle_line_idx < obstacle.points.size();
                     ++obstacle_line_idx) {
                    const Line obstacle_line(
                        obstacle.points[obstacle_line_idx - 1],
                        obstacle.points[obstacle_line_idx]);
                    Point intersection;
                    if (!arc_line.intersection(
                            obstacle_line, &intersection))
                        continue;
                    const double arc_endpoint_distance = std::min(
                        (intersection - arc.first_point())
                            .cast<double>().norm(),
                        (intersection - arc.last_point())
                            .cast<double>().norm());
                    // An arc may begin or terminate on a retained wall. Any
                    // intersection farther into the arc means the nozzle
                    // extrudes through that already printed perimeter.
                    if (arc_endpoint_distance <= scale_(0.05))
                        continue;
                    if (details != nullptr) {
                        std::ostringstream message;
                        message << "arc " << arc_idx << " segment "
                                << arc_line_idx - 1 << " crosses obstacle "
                                << obstacle_idx << " segment "
                                << obstacle_line_idx - 1 << " at "
                                << unscale<double>(intersection.x()) << ","
                                << unscale<double>(intersection.y())
                                << ", arc endpoint distance "
                                << unscale<double>(arc_endpoint_distance)
                                << " mm";
                        *details = message.str();
                    }
                    return true;
                }
            }
        }
    }
    return false;
}

static bool arc_paths_track_obstacles(const Polylines &arcs,
                                      const Polylines &obstacles,
                                      coord_t minimum_centerline_clearance,
                                      coord_t allowed_endpoint_lead,
                                      std::string *details = nullptr)
{
    const Lines obstacle_lines = to_lines(obstacles);
    const double sample_step = scale_(0.05);
    for (size_t arc_idx = 0; arc_idx < arcs.size(); ++arc_idx) {
        const Polyline &arc = arcs[arc_idx];
        const double total_length = arc.length();
        double path_position = 0.;
        for (size_t line_idx = 1; line_idx < arc.points.size(); ++line_idx) {
            const Point &start = arc.points[line_idx - 1];
            const Point &end = arc.points[line_idx];
            const Vec2d delta = (end - start).cast<double>();
            const double length = delta.norm();
            const size_t samples = std::max<size_t>(
                1, size_t(std::ceil(length / sample_step)));
            for (size_t sample_idx = 0; sample_idx <= samples; ++sample_idx) {
                const double segment_position =
                    double(sample_idx) / double(samples);
                const double position =
                    path_position + segment_position * length;
                if (position <= allowed_endpoint_lead ||
                    total_length - position <= allowed_endpoint_lead)
                    continue;
                const Vec2d sample =
                    start.cast<double>() + segment_position * delta;
                const Point point(coord_t(std::lround(sample.x())),
                                  coord_t(std::lround(sample.y())));
                const double clearance_squared = std::accumulate(
                    obstacle_lines.begin(), obstacle_lines.end(),
                    std::numeric_limits<double>::max(),
                    [&point](double nearest, const Line &obstacle) {
                        return std::min(
                            nearest,
                            obstacle.distance_to_squared(point));
                    });
                if (clearance_squared >=
                    double(minimum_centerline_clearance) *
                        double(minimum_centerline_clearance))
                    continue;
                if (details != nullptr) {
                    std::ostringstream message;
                    message << "arc " << arc_idx << " tracks a retained wall "
                            << unscale<double>(std::sqrt(clearance_squared))
                            << " mm away at "
                            << unscale<double>(point.x()) << ","
                            << unscale<double>(point.y());
                    *details = message.str();
                }
                return true;
            }
            path_position += length;
        }
    }
    return false;
}

static bool arc_paths_have_sustained_endpoint_retrace(
    const Polylines &arcs, coord_t centerline_clearance,
    coord_t minimum_run, bool trailing,
    std::string *details = nullptr)
{
    Lines printed;
    const double sample_step = scale_(0.05);
    for (size_t arc_idx = 0; arc_idx < arcs.size(); ++arc_idx) {
        Polyline inspected = arcs[arc_idx];
        if (trailing)
            inspected.reverse();
        double close_run = 0.;
        bool left_leading_contact = false;
        for (size_t point_idx = 1;
             point_idx < inspected.points.size(); ++point_idx) {
            const Vec2d start =
                inspected.points[point_idx - 1].cast<double>();
            const Vec2d delta =
                (inspected.points[point_idx] -
                 inspected.points[point_idx - 1]).cast<double>();
            const double length = delta.norm();
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
                const double distance_squared = std::accumulate(
                    printed.begin(), printed.end(),
                    std::numeric_limits<double>::max(),
                    [&sample, &delta](double nearest, const Line &line) {
                        const Vec2d printed_delta =
                            (line.b - line.a).cast<double>();
                        const double direction_norm_squared =
                            delta.squaredNorm() *
                            printed_delta.squaredNorm();
                        const double dot =
                            delta.dot(printed_delta);
                        if (direction_norm_squared <= 0. ||
                            dot * dot <
                                0.82 * direction_norm_squared)
                            return nearest;
                        return std::min(
                            nearest,
                            line.distance_to_squared(sample));
                    });
                if (distance_squared <
                    double(centerline_clearance) *
                        double(centerline_clearance)) {
                    close_run += interval;
                    if (close_run >= minimum_run) {
                        if (details != nullptr) {
                            std::ostringstream message;
                            message << "arc " << arc_idx
                                    << " has a sustained "
                                    << (trailing ? "trailing" : "leading")
                                    << " retrace for "
                                    << unscale<double>(close_run) << " mm";
                            *details = message.str();
                        }
                        return true;
                    }
                } else {
                    left_leading_contact = true;
                    break;
                }
            }
            if (left_leading_contact)
                break;
        }
        append(printed, to_lines(arcs[arc_idx]));
    }
    return false;
}

static TriangleMesh make_arc_bridge_validation_model(double width, double length)
{
    constexpr double support_depth = 3.;
    constexpr double support_height = 1.;
    constexpr double roof_height = 0.6;
    TriangleMesh model = make_cube(width, support_depth, support_height);
    TriangleMesh opposite_support =
        make_cube(width, support_depth, support_height);
    opposite_support.translate(0., length - support_depth, 0.);
    model.merge(opposite_support);
    TriangleMesh roof = make_cube(width, length, roof_height);
    roof.translate(0., 0., support_height);
    model.merge(roof);
    return model;
}

TEST_CASE("Arc overhang fill produces curved paths inside its bridge surface", "[Fill][ArcOverhang]")
{
    const ExPolygon expolygon(Points{
        Point::new_scale(0., 0.),
        Point::new_scale(20., 0.),
        Point::new_scale(20., 10.),
        Point::new_scale(0., 10.)
    });
    Surface surface(stBottomBridge, expolygon);
    surface.bridge_angle = 0.;

    std::unique_ptr<Fill> filler(Fill::new_from_type("arc-overhang"));
    REQUIRE(filler != nullptr);
    filler->spacing = 0.45;
    filler->angle = 0.f;
    filler->bounding_box = get_extents(expolygon);

    FillParams params;
    params.density = 1.f;
    params.resolution = 0.05f;

    const Polylines paths = filler->fill_surface(&surface, params);
    REQUIRE_FALSE(paths.empty());
    CHECK(std::any_of(paths.begin(), paths.end(), [](const Polyline &path) {
        return path.points.size() > 3;
    }));
    CHECK(diff_pl(paths, offset(expolygon, float(SCALED_EPSILON * 10))).empty());
}

TEST_CASE("Arc overhangs keep consecutive starts on the current side", "[Fill][ArcOverhang][Travel]")
{
    const ExPolygon expolygon(Points{
        Point::new_scale(0., 0.),
        Point::new_scale(20., 0.),
        Point::new_scale(20., 10.),
        Point::new_scale(0., 10.)
    });
    Surface surface(stBottomBridge, expolygon);
    surface.bridge_angle = 0.;

    std::unique_ptr<Fill> filler(Fill::new_from_type("arc-overhang"));
    REQUIRE(filler != nullptr);
    filler->spacing = 0.45;
    filler->bounding_box = get_extents(expolygon);

    PrintRegionConfig config;
    config.arc_overhang_recursive_fill.value = false;
    FillParams params;
    params.density = 1.f;
    params.resolution = 0.05f;
    params.config = &config;

    const Polylines paths = filler->fill_surface(&surface, params);
    REQUIRE(paths.size() > 1);
    for (size_t idx = 1; idx < paths.size(); ++idx) {
        const Point previous = paths[idx - 1].last_point();
        INFO("path index: " << idx
             << ", previous: " << unscale<double>(previous.x())
             << ", " << unscale<double>(previous.y())
             << ", first: " << unscale<double>(paths[idx].first_point().x())
             << ", " << unscale<double>(paths[idx].first_point().y())
             << ", last: " << unscale<double>(paths[idx].last_point().x())
             << ", " << unscale<double>(paths[idx].last_point().y()));
        CHECK((paths[idx].first_point() - previous).squaredNorm() <=
              (paths[idx].last_point() - previous).squaredNorm());
    }
}

TEST_CASE("Primary arc overhangs cover the bottom and start on its perimeter", "[Fill][ArcOverhang][Coverage]")
{
    const ExPolygon expolygon(Points{
        Point::new_scale(0., 0.),
        Point::new_scale(30., 0.),
        Point::new_scale(30., 20.),
        Point::new_scale(0., 20.)
    });
    Surface surface(stBottomBridge, expolygon);
    surface.bridge_angle = 0.;

    std::unique_ptr<Fill> filler(Fill::new_from_type("arc-overhang"));
    REQUIRE(filler != nullptr);
    filler->spacing = 0.45;
    filler->bounding_box = get_extents(expolygon);

    PrintRegionConfig config;
    config.arc_overhang_recursive_fill.value = false;
    FillParams params;
    params.density = 1.f;
    params.resolution = 0.05f;
    params.config = &config;

    const Polylines paths = filler->fill_surface(&surface, params);
    REQUIRE_FALSE(paths.empty());

    const coord_t spacing = scale_(filler->spacing);
    ExPolygons uncovered = diff_ex(
        ExPolygons{expolygon}, union_ex(offset(paths, float(0.52 * spacing))));
    uncovered = opening_ex(uncovered, float(0.75 * spacing));
    const double uncovered_area = std::accumulate(
        uncovered.begin(), uncovered.end(), 0.,
        [](double total, const ExPolygon &part) { return total + std::abs(part.area()); });
    CHECK(uncovered_area < 0.01 * std::abs(expolygon.area()));

    const Lines perimeter = to_lines(expolygon);
    // Fill is clipped at the perimeter's inner extrusion edge, half a line
    // width from the model contour. Starts within that overlap are supported by
    // the perimeter even though they are not on its centerline.
    const double anchor_tolerance_squared = std::pow(0.6 * spacing, 2);
    for (const Polyline &path : paths) {
        double nearest_perimeter_distance_squared = std::numeric_limits<double>::max();
        for (const Line &line : perimeter)
            nearest_perimeter_distance_squared = std::min(
                nearest_perimeter_distance_squared,
                line.distance_to_squared(path.first_point()));
        CHECK(nearest_perimeter_distance_squared <= anchor_tolerance_squared);
    }
}

TEST_CASE("Primary arc overhang starts at lower-layer support", "[Fill][ArcOverhang][Anchor]")
{
    const ExPolygon expolygon(Points{
        Point::new_scale(0., 0.),
        Point::new_scale(30., 0.),
        Point::new_scale(30., 20.),
        Point::new_scale(0., 20.)
    });
    const ExPolygons supported{
        ExPolygon(Points{
            Point::new_scale(27., 0.),
            Point::new_scale(30., 0.),
            Point::new_scale(30., 20.),
            Point::new_scale(27., 20.)
        })
    };
    Surface surface(stBottomBridge, expolygon);
    surface.bridge_angle = 0.;

    std::unique_ptr<Fill> filler(Fill::new_from_type("arc-overhang"));
    REQUIRE(filler != nullptr);
    filler->spacing = 0.45;
    filler->bounding_box = get_extents(expolygon);

    FillParams params;
    params.density = 1.f;
    params.resolution = 0.05f;
    params.arc_anchor_regions = &supported;

    const Polylines paths = filler->fill_surface(&surface, params);
    REQUIRE_FALSE(paths.empty());
    CHECK(paths.front().first_point().x() >= scale_(26.9));
    CHECK(supported.front().contains(paths.front().first_point()));

    const Polylines supported_portions = intersection_pl(
        Polylines{paths.front()}, supported);
    const auto anchored_portion = std::find_if(
        supported_portions.begin(), supported_portions.end(),
        [&paths](const Polyline &portion) {
            return portion.first_point() == paths.front().first_point() ||
                   portion.last_point() == paths.front().first_point();
        });
    REQUIRE(anchored_portion != supported_portions.end());
    // Merely touching support at the first coordinate is not printable. Keep a
    // meaningful part of the initial extrusion over the retained wall.
    CHECK(anchored_portion->length() >= scale_(0.4));
}

TEST_CASE("Primary arc prefers the preceding layer perimeter footprint",
          "[Fill][ArcOverhang][Anchor]")
{
    const ExPolygon expolygon(Points{
        Point::new_scale(0., 0.),
        Point::new_scale(30., 0.),
        Point::new_scale(30., 20.),
        Point::new_scale(0., 20.)
    });
    const ExPolygons broad_support{expolygon};
    const ExPolygons previous_perimeter{
        ExPolygon(Points{
            Point::new_scale(27., 0.),
            Point::new_scale(30., 0.),
            Point::new_scale(30., 20.),
            Point::new_scale(27., 20.)
        })
    };
    Surface surface(stBottomBridge, expolygon);
    surface.bridge_angle = 0.;

    std::unique_ptr<Fill> filler(Fill::new_from_type("arc-overhang"));
    REQUIRE(filler != nullptr);
    filler->spacing = 0.45;
    filler->bounding_box = get_extents(expolygon);

    FillParams params;
    params.density = 1.f;
    params.resolution = 0.05f;
    params.arc_anchor_regions = &broad_support;
    params.arc_root_anchor_regions = &previous_perimeter;

    const Polylines paths = filler->fill_surface(&surface, params);
    REQUIRE_FALSE(paths.empty());
    CHECK(paths.front().first_point().x() >= scale_(26.9));
    CHECK(previous_perimeter.front().contains(paths.front().first_point()));

    const Polylines supported_portions =
        intersection_pl(Polylines{paths.front()}, previous_perimeter);
    REQUIRE_FALSE(supported_portions.empty());
    CHECK(std::any_of(
        supported_portions.begin(), supported_portions.end(),
        [](const Polyline &portion) {
            return portion.length() >= scale_(0.4);
        }));
}

TEST_CASE("Arc overhangs stop at retained perimeter paths", "[Fill][ArcOverhang][Collision]")
{
    const ExPolygon expolygon(Points{
        Point::new_scale(0., 0.),
        Point::new_scale(30., 0.),
        Point::new_scale(30., 20.),
        Point::new_scale(0., 20.)
    });
    const ExPolygons supported{
        ExPolygon(Points{
            Point::new_scale(0., 0.),
            Point::new_scale(3., 0.),
            Point::new_scale(3., 20.),
            Point::new_scale(0., 20.)
        })
    };
    const Polylines retained_perimeters{
        Polyline(Points{
            Point::new_scale(12., 2.),
            Point::new_scale(12., 18.)
        })
    };
    const ExPolygons retained_perimeter_cores =
        union_ex(offset(retained_perimeters, float(scale_(0.1))));
    Surface surface(stBottomBridge, expolygon);
    surface.bridge_angle = 0.;

    std::unique_ptr<Fill> filler(Fill::new_from_type("arc-overhang"));
    REQUIRE(filler != nullptr);
    filler->spacing = 0.45;
    filler->bounding_box = get_extents(expolygon);

    PrintRegionConfig config;
    config.arc_overhang_recursive_fill.value = false;
    FillParams params;
    params.density = 1.f;
    params.resolution           = 0.05f;
    params.config               = &config;
    params.arc_anchor_regions   = &supported;
    params.arc_obstacle_paths   = &retained_perimeters;
    params.arc_obstacle_regions = &retained_perimeter_cores;

    const Polylines paths = filler->fill_surface(&surface, params);
    REQUIRE_FALSE(paths.empty());
    CHECK(std::all_of(paths.begin(), paths.end(), [](const Polyline &path) {
        return path.length() >= scale_(0.75 * 0.45);
    }));
    Lines deposited                              = to_lines(retained_perimeters);
    const double maximum_anchor_distance_squared =
        std::pow(scale_(0.45 + params.resolution), 2);
    for (const Polyline& path : paths) {
        const bool end_on_lower_support  = std::any_of(supported.begin(), supported.end(),
                                                       [&path](const ExPolygon& region) { return region.contains(path.last_point()); });
        const bool end_on_deposited_path = std::any_of(deposited.begin(), deposited.end(),
                                                       [&path, maximum_anchor_distance_squared](const Line& line) {
                                                           return line.distance_to_squared(path.last_point()) <=
                                                                  maximum_anchor_distance_squared;
                                                       });
        const bool end_is_supported      = end_on_lower_support || end_on_deposited_path;
        CHECK(end_is_supported);
        append(deposited, to_lines(path));
    }
    CHECK_FALSE(arc_paths_have_proper_crossing(paths));
    CHECK(std::any_of(paths.begin(), paths.end(), [](const Polyline& path) {
        return std::any_of(path.points.begin(), path.points.end(), [](const Point& point) { return point.x() > scale_(14.5); });
    }));
    CHECK(std::any_of(paths.begin(), paths.end(), [](const Polyline &path) {
        const bool begins_at_perimeter =
            path.first_point().x() >= scale_(12.09) &&
            path.first_point().x() <= scale_(12.11);
        return begins_at_perimeter &&
               std::any_of(path.points.begin(), path.points.end(), [](const Point &point) {
                   return point.x() > scale_(12.5);
               });
    }));
    for (const Polyline &path : paths) {
        const Polylines collision =
            intersection_pl(Polylines{path}, retained_perimeter_cores);
        CHECK(std::none_of(
            collision.begin(), collision.end(),
            [](const Polyline &portion) {
                return portion.length() > SCALED_EPSILON;
            }));
    }
    Polylines all_paths = paths;
    append(all_paths, retained_perimeters);
    CHECK_FALSE(arc_paths_have_proper_crossing(all_paths));

    // Retained centerlines are authoritative even if polygon offsetting cannot
    // produce a protected core (for example, for a very short wall fragment).
    // The arc may terminate on the wall, but it must not continue through it.
    const Polylines short_retained_perimeter{
        Polyline(Points{
            Point::new_scale(0.225, 0.),
            Point::new_scale(0.225, 20.)
        }),
        Polyline(Points{
            Point::new_scale(12., 9.9),
            Point::new_scale(12., 10.1)
        })
    };
    params.arc_obstacle_paths = &short_retained_perimeter;
    params.arc_obstacle_regions = nullptr;
    const Polylines centerline_clipped_paths =
        filler->fill_surface(&surface, params);
    REQUIRE_FALSE(centerline_clipped_paths.empty());
    CHECK(std::all_of(
        centerline_clipped_paths.begin(),
        centerline_clipped_paths.end(),
        [](const Polyline &path) {
            return path.length() >= scale_(0.75 * 0.45);
        }));
    std::string crossing_details;
    INFO(crossing_details);
    CHECK_FALSE(arc_paths_cross_obstacles(
        centerline_clipped_paths, short_retained_perimeter,
        &crossing_details));
}

TEST_CASE("Narrow bridge uses multiple small supported arc families", "[Fill][ArcOverhang][Narrow]")
{
    const ExPolygon expolygon(Points{
        Point::new_scale(0., 0.),
        Point::new_scale(3., 0.),
        Point::new_scale(3., 20.),
        Point::new_scale(0., 20.)
    });
    const ExPolygons supported{
        ExPolygon(Points{
            Point::new_scale(0., 0.),
            Point::new_scale(0.6, 0.),
            Point::new_scale(0.6, 20.),
            Point::new_scale(0., 20.)
        })
    };
    Surface surface(stBottomBridge, expolygon);
    surface.bridge_angle = 0.;

    std::unique_ptr<Fill> filler(Fill::new_from_type("arc-overhang"));
    REQUIRE(filler != nullptr);
    filler->spacing = 0.45;
    filler->bounding_box = get_extents(expolygon);

    PrintRegionConfig config;
    FillParams params;
    params.density = 1.f;
    params.resolution = 0.05f;
    params.config = &config;
    params.arc_anchor_regions = &supported;

    config.arc_overhang_recursive_fill.value = false;
    const Polylines primary_paths = filler->fill_surface(&surface, params);
    config.arc_overhang_recursive_fill.value = true;
    const Polylines paths = filler->fill_surface(&surface, params);
    REQUIRE_FALSE(primary_paths.empty());
    REQUIRE(paths.size() >= primary_paths.size());
    for (size_t path_idx = 0; path_idx < primary_paths.size(); ++path_idx)
        CHECK(paths[path_idx].points == primary_paths[path_idx].points);
    if (paths.size() > primary_paths.size()) {
        CHECK(std::any_of(
            paths.begin() + primary_paths.size(), paths.end(),
            [](const Polyline &path) {
                return get_extents(path).size().y() <= scale_(8.);
            }));
    }

    const coord_t spacing = scale_(filler->spacing);
    const auto uncovered_area = [&expolygon, spacing](const Polylines &fill_paths) {
        ExPolygons uncovered = diff_ex(
            ExPolygons{expolygon}, union_ex(offset(fill_paths, float(0.52 * spacing))));
        uncovered = opening_ex(uncovered, float(0.75 * spacing));
        return std::accumulate(
            uncovered.begin(), uncovered.end(), 0.,
            [](double total, const ExPolygon &part) {
                return total + std::abs(part.area());
            });
    };
    CHECK(uncovered_area(primary_paths) < 0.05 * std::abs(expolygon.area()));
    CHECK(uncovered_area(paths) <= uncovered_area(primary_paths));

    Lines printed_bridge;
    for (const Polyline &path : primary_paths) {
        bool anchored = std::any_of(
            supported.begin(), supported.end(),
            [&path](const ExPolygon &region) {
                return region.contains(path.first_point());
            });
        if (!anchored) {
            double nearest_printed_distance_squared =
                std::numeric_limits<double>::max();
            for (const Line &line : printed_bridge)
                nearest_printed_distance_squared = std::min(
                    nearest_printed_distance_squared,
                    line.distance_to_squared(path.first_point()));
            anchored = nearest_printed_distance_squared <=
                       double(spacing) * double(spacing);
        }
        CHECK(anchored);

        if (path.points.size() >= 3) {
            const Line chord(path.first_point(), path.last_point());
            if (chord.length() >= 4. * spacing) {
                double max_sagitta_squared = 0.;
                for (const Point &point : path.points)
                    max_sagitta_squared = std::max(
                        max_sagitta_squared,
                        chord.distance_to_squared(point));
                CHECK(max_sagitta_squared >=
                      double(spacing) * double(spacing));
            }
        }
        append(printed_bridge, to_lines(path));
    }
}

TEST_CASE("Four millimeter bridge uses anchored curved arches instead of shallow chords",
          "[Fill][ArcOverhang][Narrow][Anchor]")
{
    // This is the important scale of the narrow slots seen in real printer
    // test models. It is wider than the old eight-pitch cutoff, but expanding
    // circles clipped to this long rectangle still flatten into bridge lines.
    const ExPolygon expolygon(Points{
        Point::new_scale(0., 0.),
        Point::new_scale(4.5, 0.),
        Point::new_scale(4.5, 40.),
        Point::new_scale(0., 40.)
    });
    const ExPolygons supported{
        ExPolygon(Points{
            Point::new_scale(0., 0.),
            Point::new_scale(0.6, 0.),
            Point::new_scale(0.6, 40.),
            Point::new_scale(0., 40.)
        })
    };
    const Polylines retained_perimeters{
        Polyline(Points{
            Point::new_scale(0.225, 0.),
            Point::new_scale(0.225, 40.)
        })
    };
    const ExPolygons retained_perimeter_cores =
        union_ex(offset(retained_perimeters, float(scale_(0.1))));
    Surface surface(stBottomBridge, expolygon);
    surface.bridge_angle = 0.;

    std::unique_ptr<Fill> filler(Fill::new_from_type("arc-overhang"));
    REQUIRE(filler != nullptr);
    filler->spacing = 0.45;
    filler->bounding_box = get_extents(expolygon);

    PrintRegionConfig config;
    config.arc_overhang_recursive_fill.value = false;
    FillParams params;
    params.density = 1.f;
    params.resolution = 0.05f;
    params.config = &config;
    params.arc_anchor_regions = &supported;
    params.arc_obstacle_paths = &retained_perimeters;
    params.arc_obstacle_regions = &retained_perimeter_cores;

    const Polylines paths = filler->fill_surface(&surface, params);
    // Exact collision trimming may merge or discard a redundant terminal
    // fragment, so validate deposited coverage rather than requiring the old
    // implementation's exact path count.
    REQUIRE(paths.size() >= 19);

    const coord_t spacing = scale_(filler->spacing);
    Lines printed;
    const Lines perimeter_lines = to_lines(retained_perimeters);
    for (const Polyline &path : paths) {
        double anchor_distance_squared = std::numeric_limits<double>::max();
        for (const Line &line : perimeter_lines)
            anchor_distance_squared = std::min(
                anchor_distance_squared,
                line.distance_to_squared(path.first_point()));
        for (const Line &line : printed)
            anchor_distance_squared = std::min(
                anchor_distance_squared,
                line.distance_to_squared(path.first_point()));
        CHECK(anchor_distance_squared <=
              double(spacing) * double(spacing));

        if (path.points.size() >= 3) {
            const Line chord(path.first_point(), path.last_point());
            const bool fully_laterally_supported =
                std::all_of(
                    path.points.begin(), path.points.end(),
                    [&perimeter_lines, &printed,
                     spacing](const Point &point) {
                        double distance_squared =
                            std::numeric_limits<double>::max();
                        for (const Line &line :
                             perimeter_lines)
                            distance_squared = std::min(
                                distance_squared,
                                line.distance_to_squared(
                                    point));
                        for (const Line &line : printed)
                            distance_squared = std::min(
                                distance_squared,
                                line.distance_to_squared(
                                    point));
                        return distance_squared <=
                               double(spacing) *
                                   double(spacing);
                    });
            if (chord.length() >= 4. * spacing &&
                !fully_laterally_supported) {
                double max_sagitta_squared = 0.;
                for (const Point &point : path.points)
                    max_sagitta_squared = std::max(
                        max_sagitta_squared,
                        chord.distance_to_squared(point));
                CHECK(max_sagitta_squared >=
                      double(spacing) * double(spacing));
            }
        }
        append(printed, to_lines(path));
    }
    ExPolygons covered =
        union_ex(offset(paths, float(0.52 * spacing)));
    append(covered, supported);
    ExPolygons uncovered =
        diff_ex(ExPolygons{expolygon}, union_ex(covered));
    uncovered = opening_ex(uncovered, float(0.75 * spacing));
    const double uncovered_area = std::accumulate(
        uncovered.begin(), uncovered.end(), 0.,
        [](double total, const ExPolygon &part) {
            return total + std::abs(part.area());
        });
    CHECK(uncovered_area < 0.08 * std::abs(expolygon.area()));
    CHECK_FALSE(arc_paths_have_proper_crossing(paths));
}

TEST_CASE("Narrow bridge curvature and coverage hold across practical nozzle-scale widths",
          "[Fill][ArcOverhang][Narrow][Coverage]")
{
    const double width = GENERATE(2.4, 3.2, 4.5, 5.8, 6.8);
    const bool recursive_fill = GENERATE(false, true);
    DYNAMIC_SECTION("width " << width << " mm, recursive " << recursive_fill)
    {
        constexpr double length = 40.;
        const ExPolygon expolygon(Points{
            Point::new_scale(0., 0.),
            Point::new_scale(width, 0.),
            Point::new_scale(width, length),
            Point::new_scale(0., length)
        });
        const ExPolygons supported{
            ExPolygon(Points{
                Point::new_scale(0., 0.),
                Point::new_scale(0.6, 0.),
                Point::new_scale(0.6, length),
                Point::new_scale(0., length)
            })
        };
        Surface surface(stBottomBridge, expolygon);
        surface.bridge_angle = 0.;

        std::unique_ptr<Fill> filler(Fill::new_from_type("arc-overhang"));
        REQUIRE(filler != nullptr);
        filler->spacing = 0.45;
        filler->bounding_box = get_extents(expolygon);

        PrintRegionConfig config;
        config.arc_overhang_recursive_fill.value = recursive_fill;
        FillParams params;
        params.density = 1.f;
        params.resolution = 0.05f;
        params.config = &config;
        params.arc_anchor_regions = &supported;

        const Polylines paths = filler->fill_surface(&surface, params);
        REQUIRE_FALSE(paths.empty());
        const coord_t spacing = scale_(filler->spacing);
        Lines printed;
        for (size_t path_idx = 0; path_idx < paths.size(); ++path_idx) {
            const Polyline &path = paths[path_idx];
            INFO("path index: " << path_idx
                 << ", start: " << unscale<double>(path.first_point().x())
                 << ", " << unscale<double>(path.first_point().y()));
            bool anchored = std::any_of(
                supported.begin(), supported.end(),
                [&path](const ExPolygon &region) {
                    return region.contains(path.first_point());
                });
            if (!anchored) {
                double nearest_distance_squared =
                    std::numeric_limits<double>::max();
                for (const Line &line : printed)
                    nearest_distance_squared = std::min(
                        nearest_distance_squared,
                        line.distance_to_squared(path.first_point()));
                anchored = nearest_distance_squared <=
                           double(spacing) * double(spacing);
            }
            CHECK(anchored);

            const bool fully_supported = std::any_of(
                supported.begin(), supported.end(),
                [&path](const ExPolygon &region) {
                    return region.contains(path);
                });
            if (!fully_supported) {
                CHECK(path.points.size() >= 3);
            }
            if (!fully_supported && path.points.size() >= 3) {
                const Line chord(path.first_point(), path.last_point());
                if (chord.length() >= 2. * spacing) {
                    double max_sagitta_squared = 0.;
                    for (const Point &point : path.points)
                        max_sagitta_squared = std::max(
                            max_sagitta_squared,
                            chord.distance_to_squared(point));
                    const double minimum_sagitta =
                        std::min(double(spacing), 0.10 * chord.length());
                    CHECK(max_sagitta_squared >=
                          minimum_sagitta * minimum_sagitta);
                }
            }
            append(printed, to_lines(path));
        }

        ExPolygons covered =
            union_ex(offset(paths, float(0.52 * spacing)));
        append(covered, supported);
        ExPolygons uncovered =
            diff_ex(ExPolygons{expolygon}, union_ex(covered));
        uncovered = opening_ex(uncovered, float(0.75 * spacing));
        const double uncovered_area = std::accumulate(
            uncovered.begin(), uncovered.end(), 0.,
            [](double total, const ExPolygon &part) {
                return total + std::abs(part.area());
            });
        CHECK(uncovered_area < 0.08 * std::abs(expolygon.area()));
        CHECK_FALSE(arc_paths_have_proper_crossing(paths));
    }
}

TEST_CASE("One-sided narrow overhang fills from supported arc origins",
          "[Fill][ArcOverhang][Narrow][Anchor][Coverage]")
{
    const ExPolygon expolygon(Points{
        Point::new_scale(0., 0.),
        Point::new_scale(4.5, 0.),
        Point::new_scale(4.5, 20.),
        Point::new_scale(0., 20.)
    });
    const ExPolygons supported{
        ExPolygon(Points{
            Point::new_scale(0., 0.),
            Point::new_scale(0.6, 0.),
            Point::new_scale(0.6, 20.),
            Point::new_scale(0., 20.)
        })
    };
    Surface surface(stBottomBridge, expolygon);
    surface.bridge_angle = 0.;

    std::unique_ptr<Fill> filler(Fill::new_from_type("arc-overhang"));
    REQUIRE(filler != nullptr);
    filler->spacing = 0.45;
    filler->bounding_box = get_extents(expolygon);

    PrintRegionConfig config;
    config.arc_overhang_recursive_fill.value = true;
    FillParams params;
    params.density = 1.f;
    params.resolution = 0.05f;
    params.config = &config;
    params.arc_anchor_regions = &supported;

    const Polylines paths = filler->fill_surface(&surface, params);
    REQUIRE_FALSE(paths.empty());
    const coord_t spacing = scale_(filler->spacing);
    const coord_t arc_pitch = scale_(
        filler->spacing -
        0.45 * config.arc_overhang_overlap.value / 100.);
    const coord_t retrace_clearance =
        std::max<coord_t>(
            coord_t(std::lround(scale_(0.20 * 0.45))),
            coord_t(std::lround(0.55 * double(arc_pitch))));
    const double maximum_distance_squared =
        std::pow(scale_(0.45 + params.resolution), 2);
    Lines deposited;
    const auto point_is_supported =
        [&supported, &deposited,
         maximum_distance_squared](const Point &point) {
            if (std::any_of(
                    supported.begin(), supported.end(),
                    [&point](const ExPolygon &region) {
                        return region.contains(point);
                    }))
                return true;
            return std::any_of(
                deposited.begin(), deposited.end(),
                [&point, maximum_distance_squared](const Line &line) {
                    return line.distance_to_squared(point) <=
                           maximum_distance_squared;
                });
        };
    const auto lead_is_supported =
        [&point_is_supported](const Polyline &path) {
            double remaining = scale_(0.75 * 0.45);
            const double sample_step = scale_(0.20 * 0.45);
            for (size_t point_idx = 1;
                 point_idx < path.points.size() && remaining > 0.;
                 ++point_idx) {
                const Vec2d start =
                    path.points[point_idx - 1].cast<double>();
                const Vec2d delta =
                    (path.points[point_idx] -
                     path.points[point_idx - 1]).cast<double>();
                const double length = delta.norm();
                if (length <= 0.)
                    continue;
                const double checked = std::min(length, remaining);
                const size_t samples = std::max<size_t>(
                    1, size_t(std::ceil(checked / sample_step)));
                for (size_t sample_idx = 0;
                     sample_idx <= samples; ++sample_idx) {
                    const Vec2d position =
                        start + delta *
                                    (checked * double(sample_idx) /
                                     double(samples) / length);
                    if (!point_is_supported(
                            Point(coord_t(std::lround(position.x())),
                                  coord_t(std::lround(position.y())))))
                        return false;
                }
                remaining -= checked;
            }
            return true;
        };

    for (const Polyline &path : paths) {
        // No sub-bead specks, but retain three-quarter-bead recursive
        // fragments: their deposited width closes visible pockets which a
        // full-bead centerline threshold would leave open.
        CHECK(path.length() >= scale_(0.75 * 0.45));
        CHECK(lead_is_supported(path));
        append(deposited, to_lines(path));
    }

    std::string retrace_details;
    const bool has_sustained_leading_retrace =
        arc_paths_have_sustained_endpoint_retrace(
            paths, retrace_clearance,
            scale_(0.65 * 0.45), false, &retrace_details);
    INFO(retrace_details);
    CHECK_FALSE(has_sustained_leading_retrace);

    retrace_details.clear();
    const bool has_sustained_trailing_retrace =
        arc_paths_have_sustained_endpoint_retrace(
            paths, retrace_clearance,
            scale_(0.65 * 0.45), true, &retrace_details);
    INFO(retrace_details);
    CHECK_FALSE(has_sustained_trailing_retrace);

    ExPolygons covered =
        union_ex(offset(paths, float(0.52 * spacing)));
    append(covered, supported);
    ExPolygons uncovered =
        diff_ex(ExPolygons{expolygon}, union_ex(covered));
    uncovered = opening_ex(uncovered, float(0.75 * spacing));
    const double uncovered_area = std::accumulate(
        uncovered.begin(), uncovered.end(), 0.,
        [](double total, const ExPolygon &part) {
            return total + std::abs(part.area());
        });
    CHECK(uncovered_area < 0.08 * std::abs(expolygon.area()));
    CHECK_FALSE(arc_paths_have_proper_crossing(paths));
}

TEST_CASE("Narrow bridge grows tight arcs from both supported ends", "[Fill][ArcOverhang][Narrow]")
{
    const ExPolygon expolygon(Points{
        Point::new_scale(0., 0.),
        Point::new_scale(4.5, 0.),
        Point::new_scale(4.5, 12.),
        Point::new_scale(0., 12.)
    });
    const ExPolygons supported{
        ExPolygon(Points{
            Point::new_scale(0., 0.),
            Point::new_scale(4.5, 0.),
            Point::new_scale(4.5, 0.6),
            Point::new_scale(0., 0.6)
        }),
        ExPolygon(Points{
            Point::new_scale(0., 11.4),
            Point::new_scale(4.5, 11.4),
            Point::new_scale(4.5, 12.),
            Point::new_scale(0., 12.)
        })
    };
    Surface surface(stBottomBridge, expolygon);
    surface.bridge_angle = 0.;

    std::unique_ptr<Fill> filler(Fill::new_from_type("arc-overhang"));
    REQUIRE(filler != nullptr);
    filler->spacing = 0.45;
    filler->bounding_box = get_extents(expolygon);

    PrintRegionConfig config;
    config.arc_overhang_recursive_fill.value = true;
    FillParams params;
    params.density = 1.f;
    params.resolution = 0.05f;
    params.config = &config;
    params.arc_anchor_regions = &supported;

    const Polylines paths = filler->fill_surface(&surface, params);
    REQUIRE_FALSE(paths.empty());
    CHECK(std::any_of(paths.begin(), paths.end(), [](const Polyline &path) {
        return path.first_point().y() <= scale_(0.65);
    }));
    CHECK(std::any_of(paths.begin(), paths.end(), [](const Polyline &path) {
        return path.first_point().y() >= scale_(11.35);
    }));

    const double spacing = scale_(filler->spacing);
    size_t tight_path_count = 0;
    for (const Polyline &path : paths) {
        if (path.points.size() < 3)
            continue;
        const Line chord(path.first_point(), path.last_point());
        if ((path.last_point() - path.first_point()).cast<double>().norm() < 4. * spacing)
            continue;
        double max_sagitta_squared = 0.;
        for (const Point &point : path.points)
            max_sagitta_squared = std::max(max_sagitta_squared, chord.distance_to_squared(point));
        if (max_sagitta_squared >= std::pow(1.2 * spacing, 2))
            ++tight_path_count;
    }
    // The complete large family is the printable backbone. The later family
    // from the opposite support must still contribute tighter arches instead
    // of leaving every narrow-span path as a shallow chord.
    CHECK(tight_path_count > 0);

    // Independently rounding both supported half-families down used to leave a
    // complete missing arch where they met. The phase-aligned opposite family
    // must reach within one deposited bead of the first family at the bridge
    // midpoint; demanding a rectangular footprint here would incorrectly
    // reject the deliberately curved space inside each arch.
    const double midpoint_y = scale_(6.);
    double lower_reach = std::numeric_limits<double>::lowest();
    double upper_reach = std::numeric_limits<double>::max();
    for (const Polyline &path : paths) {
        for (const Point &point : path.points) {
            const double y = point.y();
            if (y <= midpoint_y)
                lower_reach = std::max(lower_reach, y);
            if (y >= midpoint_y)
                upper_reach = std::min(upper_reach, y);
        }
    }
    REQUIRE(lower_reach > std::numeric_limits<double>::lowest());
    REQUIRE(upper_reach < std::numeric_limits<double>::max());
    CHECK(upper_reach - lower_reach <= spacing);
}

TEST_CASE("Long narrow arc overhang retains complete primary coverage", "[Fill][ArcOverhang][Coverage][Performance]")
{
    const ExPolygon expolygon(Points{
        Point::new_scale(0., 0.),
        Point::new_scale(120., 0.),
        Point::new_scale(120., 8.),
        Point::new_scale(0., 8.)
    });
    const ExPolygons supported{
        ExPolygon(Points{
            Point::new_scale(0., 0.),
            Point::new_scale(2., 0.),
            Point::new_scale(2., 8.),
            Point::new_scale(0., 8.)
        }),
        ExPolygon(Points{
            Point::new_scale(118., 0.),
            Point::new_scale(120., 0.),
            Point::new_scale(120., 8.),
            Point::new_scale(118., 8.)
        })
    };
    Surface surface(stBottomBridge, expolygon);
    surface.bridge_angle = 0.;

    std::unique_ptr<Fill> filler(Fill::new_from_type("arc-overhang"));
    REQUIRE(filler != nullptr);
    filler->spacing = 0.45;
    filler->bounding_box = get_extents(expolygon);

    PrintRegionConfig config;
    FillParams params;
    params.density = 1.f;
    params.resolution = 0.05f;
    params.config = &config;
    params.arc_anchor_regions = &supported;

    const coord_t spacing = scale_(filler->spacing);
    const auto check_coverage = [&](bool recursive_fill) {
        INFO("recursive fill: " << recursive_fill);
        config.arc_overhang_recursive_fill.value = recursive_fill;
        const Polylines paths = filler->fill_surface(&surface, params);
        REQUIRE_FALSE(paths.empty());
        CHECK(std::any_of(paths.begin(), paths.end(), [](const Polyline &path) {
            return path.first_point().x() <= scale_(2.05);
        }));
        CHECK(std::any_of(paths.begin(), paths.end(), [](const Polyline &path) {
            return path.first_point().x() >= scale_(117.95);
        }));

        for (size_t path_idx = 0; path_idx < paths.size(); ++path_idx) {
            const Polyline &path = paths[path_idx];
            if (path.points.size() < 3)
                continue;
            const Line chord(path.first_point(), path.last_point());
            if (chord.length() < 4. * spacing)
                continue;
            double max_sagitta_squared = 0.;
            for (const Point &point : path.points)
                max_sagitta_squared = std::max(
                    max_sagitta_squared, chord.distance_to_squared(point));
            // Bounded-radius chained families must remain visibly curved. A
            // long clipped segment with less than one line width of curvature
            // behaves like the sagging straight bridges this pattern replaces.
            INFO("path index: " << path_idx
                 << ", path length: " << unscale<double>(path.length())
                 << ", chord length: " << unscale<double>(chord.length()));
            CHECK(max_sagitta_squared >= double(spacing) * double(spacing));
        }

        ExPolygons covered = union_ex(offset(paths, float(0.52 * spacing)));
        append(covered, supported);
        ExPolygons uncovered = diff_ex(ExPolygons{expolygon}, union_ex(covered));
        uncovered = opening_ex(uncovered, float(0.75 * spacing));
        const double uncovered_area = std::accumulate(
            uncovered.begin(), uncovered.end(), 0.,
            [](double total, const ExPolygon &part) {
                return total + std::abs(part.area());
            });
        CHECK(uncovered_area < 0.02 * std::abs(expolygon.area()));
        CHECK_FALSE(arc_paths_have_proper_crossing(paths));
    };
    check_coverage(false);
    check_coverage(true);
}

TEST_CASE("Arc overlap reduces free-air line pitch", "[Fill][ArcOverhang][Spacing]")
{
    const ExPolygon expolygon(Points{
        Point::new_scale(0., 0.),
        Point::new_scale(40., 0.),
        Point::new_scale(40., 8.),
        Point::new_scale(0., 8.)
    });
    const ExPolygons supported{
        ExPolygon(Points{
            Point::new_scale(0., 0.),
            Point::new_scale(2., 0.),
            Point::new_scale(2., 8.),
            Point::new_scale(0., 8.)
        }),
        ExPolygon(Points{
            Point::new_scale(38., 0.),
            Point::new_scale(40., 0.),
            Point::new_scale(40., 8.),
            Point::new_scale(38., 8.)
        })
    };
    Surface surface(stBottomBridge, expolygon);
    surface.bridge_angle = 0.;

    std::unique_ptr<Fill> filler(Fill::new_from_type("arc-overhang"));
    REQUIRE(filler != nullptr);
    filler->spacing = 0.45;
    filler->bounding_box = get_extents(expolygon);

    PrintRegionConfig config;
    config.arc_overhang_recursive_fill.value = false;
    FillParams params;
    params.density = 1.f;
    params.resolution = 0.05f;
    params.config = &config;
    params.arc_anchor_regions = &supported;

    config.arc_overhang_overlap.value = 0.;
    const Polylines touching_paths = filler->fill_surface(&surface, params);
    config.arc_overhang_overlap.value = 20.;
    const Polylines overlapping_paths = filler->fill_surface(&surface, params);

    REQUIRE_FALSE(touching_paths.empty());
    CHECK(overlapping_paths.size() > touching_paths.size());
}

TEST_CASE("Recursive arc refinement follows the complete large family", "[Fill][ArcOverhang][Recursive][Ordering]")
{
    const ExPolygon expolygon(Points{
        Point::new_scale(0., 0.),
        Point::new_scale(80., 0.),
        Point::new_scale(80., 8.),
        Point::new_scale(0., 8.)
    });
    const ExPolygons supported{
        ExPolygon(Points{
            Point::new_scale(0., 0.),
            Point::new_scale(0.6, 0.),
            Point::new_scale(0.6, 8.),
            Point::new_scale(0., 8.)
        })
    };
    Surface surface(stBottomBridge, expolygon);
    surface.bridge_angle = 0.;

    std::unique_ptr<Fill> filler(Fill::new_from_type("arc-overhang"));
    REQUIRE(filler != nullptr);
    filler->spacing = 0.45;
    filler->bounding_box = get_extents(expolygon);

    PrintRegionConfig config;
    FillParams params;
    params.density = 1.f;
    params.resolution = 0.05f;
    params.config = &config;
    params.arc_anchor_regions = &supported;

    config.arc_overhang_recursive_fill.value = false;
    const Polylines primary_paths = filler->fill_surface(&surface, params);
    config.arc_overhang_recursive_fill.value = true;
    const Polylines recursive_paths = filler->fill_surface(&surface, params);

    REQUIRE_FALSE(primary_paths.empty());
    REQUIRE(recursive_paths.size() >= primary_paths.size());
    // Recursive output is deliberately append-only: the complete large family
    // remains the prefix, and smaller gap-filling arches follow it.
    for (size_t path_idx = 0; path_idx < primary_paths.size(); ++path_idx)
        CHECK(recursive_paths[path_idx].points == primary_paths[path_idx].points);

    const coord_t spacing = scale_(filler->spacing);
    ExPolygons covered = union_ex(offset(primary_paths, float(0.52 * spacing)));
    append(covered, supported);
    ExPolygons uncovered = diff_ex(ExPolygons{expolygon}, union_ex(covered));
    uncovered = opening_ex(uncovered, float(0.75 * spacing));
    const double uncovered_area = std::accumulate(
        uncovered.begin(), uncovered.end(), 0.,
        [](double total, const ExPolygon &part) {
            return total + std::abs(part.area());
        });
    CHECK(uncovered_area < 0.02 * std::abs(expolygon.area()));
}

TEST_CASE("Recursive arc fill branches into uncovered space", "[Fill][ArcOverhang][Recursive]")
{
    ExPolygon expolygon(Points{
        Point::new_scale(0., 0.),
        Point::new_scale(30., 0.),
        Point::new_scale(30., 20.),
        Point::new_scale(0., 20.)
    });
    expolygon.holes.emplace_back(Points{
        Point::new_scale(8., 6.),
        Point::new_scale(8., 14.),
        Point::new_scale(12., 14.),
        Point::new_scale(12., 6.)
    });
    Surface surface(stBottomBridge, expolygon);
    surface.bridge_angle = 0.;

    std::unique_ptr<Fill> filler(Fill::new_from_type("arc-overhang"));
    REQUIRE(filler != nullptr);
    filler->spacing = 0.45;
    filler->bounding_box = get_extents(expolygon);

    PrintRegionConfig config;
    FillParams params;
    params.density = 1.f;
    params.resolution = 0.05f;
    params.config = &config;

    config.arc_overhang_recursive_fill.value = false;
    const Polylines root_paths = filler->fill_surface(&surface, params);
    config.arc_overhang_recursive_fill.value = true;
    const Polylines recursive_paths = filler->fill_surface(&surface, params);

    // The complete root family may leave only sub-bead pockets around the
    // hole. Recursive mode must preserve that family, but should not emit
    // unprintable micro-arcs merely to increase the path count.
    REQUIRE(recursive_paths.size() >= root_paths.size());

    const coord_t spacing = scale_(filler->spacing);
    const auto uncovered_regions = [&expolygon, spacing](const Polylines &paths) {
        ExPolygons uncovered = diff_ex(
            ExPolygons{expolygon}, union_ex(offset(paths, float(0.52 * spacing))));
        uncovered = opening_ex(uncovered, float(0.75 * spacing));
        return uncovered;
    };
    const auto uncovered_area = [&uncovered_regions](const Polylines &paths) {
        const ExPolygons uncovered = uncovered_regions(paths);
        return std::accumulate(uncovered.begin(), uncovered.end(), 0.,
            [](double total, const ExPolygon &part) { return total + std::abs(part.area()); });
    };
    const double root_uncovered_area = uncovered_area(root_paths);
    INFO("root paths: " << root_paths.size() << ", recursive paths: " << recursive_paths.size());
    INFO("root uncovered: " << root_uncovered_area << ", recursive uncovered: " << uncovered_area(recursive_paths));
    const ExPolygons recursive_uncovered = uncovered_regions(recursive_paths);
    const auto largest_recursive_region = std::max_element(
        recursive_uncovered.begin(), recursive_uncovered.end(),
        [](const ExPolygon &lhs, const ExPolygon &rhs) {
            return std::abs(lhs.area()) < std::abs(rhs.area());
        });
    const double largest_recursive_area =
        largest_recursive_region == recursive_uncovered.end() ? 0. : std::abs(largest_recursive_region->area());
    INFO("recursive regions: " << recursive_uncovered.size() << ", largest: "
         << largest_recursive_area);
    CHECK(uncovered_area(recursive_paths) <= root_uncovered_area);
    CHECK(largest_recursive_area < 10. * double(spacing) * double(spacing));
    CHECK_FALSE(arc_paths_have_proper_crossing(root_paths));
    CHECK_FALSE(arc_paths_have_proper_crossing(recursive_paths));

    Lines printed_lines = to_lines(root_paths);
    const double maximum_anchor_distance_squared = std::pow(1.1 * spacing, 2);
    for (size_t path_idx = root_paths.size(); path_idx < recursive_paths.size(); ++path_idx) {
        double nearest_printed_distance_squared = std::numeric_limits<double>::max();
        for (const Line &line : printed_lines)
            nearest_printed_distance_squared = std::min(
                nearest_printed_distance_squared,
                line.distance_to_squared(recursive_paths[path_idx].first_point()));
        CHECK(nearest_printed_distance_squared <= maximum_anchor_distance_squared);
        append(printed_lines, to_lines(recursive_paths[path_idx]));
    }
}

TEST_CASE("Arc overhang keeps a distinct G-code feature role", "[Fill][ArcOverhang][GCode]")
{
    CHECK(ExtrusionEntity::role_to_string(erArcOverhang) == "Arc overhang");
    CHECK(ExtrusionEntity::string_to_role("Arc overhang") == erArcOverhang);
    CHECK(ExtrusionEntity::role_to_string(erArcBridge) == "Arc bridge");
    CHECK(ExtrusionEntity::string_to_role("Arc bridge") == erArcBridge);
}

TEST_CASE("Wide unsupported roofs reach the preview as arc overhangs", "[Fill][ArcOverhang][GCode]")
{
    TriangleMesh model = make_cube(10., 20., 5.);
    TriangleMesh roof  = make_cube(30., 20., 2.);
    roof.translate(-10., 0., 5.);
    model.merge(roof);

    Print print;
    Slic3r::Test::init_and_process_print({model}, print, {
        {"arc_overhang_enabled", "1"},
        {"arc_overhang_bridge_distance", "0"},
        {"arc_overhang_min_overhang_distance", "0"},
        {"enable_support", "0"}
    });

    CHECK(print_has_arc_overhang(print));

    ScopedTemporaryFile gcode_file(".gcode");
    GCodeProcessorResult preview;
    print.export_gcode(gcode_file.string(), &preview, nullptr);
    const auto first_arc = std::find_if(preview.moves.begin(), preview.moves.end(), [](const auto &move) {
        return move.type == EMoveType::Extrude && move.extrusion_role == erArcOverhang;
    });
    REQUIRE(first_arc != preview.moves.end());
    REQUIRE(first_arc != preview.moves.begin());
    CHECK(std::any_of(preview.moves.begin(), first_arc, [](const auto &move) {
        return move.type == EMoveType::Extrude && is_perimeter(move.extrusion_role);
    }));
    CHECK(std::none_of(preview.moves.begin(), preview.moves.end(), [](const auto &move) {
        return move.type == EMoveType::Extrude && move.extrusion_role == erOverhangPerimeter;
    }));

    const Vec2d arc_start = std::prev(first_arc)->position.head<2>().cast<double>();
    double minimum_bead_clearance = std::numeric_limits<double>::max();
    double previous_layer_z = std::numeric_limits<double>::lowest();
    for (auto move = preview.moves.begin(); move != first_arc; ++move)
        if (move->type == EMoveType::Extrude &&
            move->position.z() < first_arc->position.z() - 0.001)
            previous_layer_z =
                std::max(previous_layer_z, double(move->position.z()));
    REQUIRE(previous_layer_z > std::numeric_limits<double>::lowest());
    for (auto move = std::next(preview.moves.begin()); move != first_arc; ++move) {
        if (move->type != EMoveType::Extrude ||
            !is_perimeter(move->extrusion_role) ||
            std::abs(move->position.z() - previous_layer_z) > 0.001)
            continue;
        const Vec2d line_start = std::prev(move)->position.head<2>().cast<double>();
        const Vec2d line_delta =
            move->position.head<2>().cast<double>() - line_start;
        const double length_squared = line_delta.squaredNorm();
        if (length_squared == 0.)
            continue;
        const double position = std::clamp(
            (arc_start - line_start).dot(line_delta) / length_squared, 0., 1.);
        const double centerline_distance =
            (arc_start - line_start - position * line_delta).norm();
        minimum_bead_clearance = std::min(
            minimum_bead_clearance,
            centerline_distance - 0.5 * (move->width + first_arc->width));
    }
    // The first arc bead must overlap the actual perimeter emitted immediately
    // below it. A same-layer wall or broad lower-layer slice is not sufficient
    // to establish the intended shape before the path enters free air.
    INFO("first arc start: " << arc_start.x() << "," << arc_start.y()
         << ", width: " << first_arc->width
         << ", previous-layer bead clearance: "
         << minimum_bead_clearance);
    CHECK(minimum_bead_clearance <= 0.);
}

TEST_CASE("Bridge and overhang ordering can precede walls on only the affected layer", "[Fill][ArcOverhang][GCode][Ordering]")
{
    const bool use_arc_overhangs = GENERATE(false, true);
    CAPTURE(use_arc_overhangs);
    TriangleMesh model = make_cube(10., 20., 5.);
    TriangleMesh roof  = make_cube(30., 20., 2.);
    roof.translate(-10., 0., 5.);
    model.merge(roof);

    Print print;
    Slic3r::Test::init_and_process_print({model}, print, {
        {"arc_overhang_enabled", use_arc_overhangs ? "1" : "0"},
        {"arc_overhang_bridge_distance", "0"},
        {"arc_overhang_min_overhang_distance", "0"},
        {"bridge_overhang_before_walls", "1"},
        {"enable_support", "0"}
    });

    ScopedTemporaryFile gcode_file(".gcode");
    GCodeProcessorResult preview;
    print.export_gcode(gcode_file.string(), &preview, nullptr);
    const auto first_supporting_fill = std::find_if(preview.moves.begin(), preview.moves.end(), [](const auto &move) {
        return move.type == EMoveType::Extrude && is_bridge_or_arc_fill(move.extrusion_role);
    });
    REQUIRE(first_supporting_fill != preview.moves.end());
    if (use_arc_overhangs)
        CHECK(is_arc_fill(first_supporting_fill->extrusion_role));
    else
        CHECK((first_supporting_fill->extrusion_role == erBridgeInfill ||
               first_supporting_fill->extrusion_role == erInternalBridgeInfill));
    const float affected_z = first_supporting_fill->position.z();

    CHECK(std::none_of(preview.moves.begin(), first_supporting_fill, [affected_z](const auto &move) {
        return move.type == EMoveType::Extrude && is_perimeter(move.extrusion_role) &&
               std::abs(move.position.z() - affected_z) < 0.001f;
    }));
    CHECK(std::any_of(std::next(first_supporting_fill), preview.moves.end(), [affected_z](const auto &move) {
        return move.type == EMoveType::Extrude && is_perimeter(move.extrusion_role) &&
               std::abs(move.position.z() - affected_z) < 0.001f;
    }));
}

TEST_CASE("Narrow bridge G-code keeps every arc path anchored and curved",
          "[Fill][ArcOverhang][GCode][Narrow]")
{
    const double width = GENERATE(3.2, 4.5, 6.8);
    const bool recursive_fill = GENERATE(false, true);
    DYNAMIC_SECTION("width " << width << " mm, recursive " << recursive_fill)
    {
        Print print;
        Slic3r::Test::init_and_process_print(
            {make_arc_bridge_validation_model(width, 30.)}, print, {
                {"layer_height", "0.2"},
                {"initial_layer_print_height", "0.2"},
                {"wall_loops", "2"},
                {"arc_overhang_enabled", "1"},
                {"arc_overhang_bridge_distance", "0"},
                {"arc_overhang_min_overhang_distance", "0"},
                {"arc_overhang_recursive_fill", recursive_fill ? "1" : "0"},
                {"enable_arc_fitting", "1"},
                {"enable_support", "0"}
            });
        REQUIRE(print_has_arc_overhang(print));
        for (const PrintObject *object : print.objects()) {
            for (const Layer *layer : object->layers()) {
                Polylines raw_arc_paths;
                for (const LayerRegion *region : layer->regions()) {
                    for (const ExtrusionEntity *entity :
                         region->fills.flatten().entities) {
                        if (const auto *path =
                                dynamic_cast<const ExtrusionPath *>(entity);
                            path != nullptr &&
                            is_arc_fill(path->role())) {
                            raw_arc_paths.emplace_back();
                            for (const Point3 &point :
                                 path->polyline.points)
                                raw_arc_paths.back().append(
                                    point.to_point());
                        }
                    }
                }
                if (!raw_arc_paths.empty()) {
                    std::string crossing_details;
                    const bool has_crossing =
                        arc_paths_have_proper_crossing(
                            raw_arc_paths, &crossing_details);
                    INFO("raw arc-overhang layer " << layer->id());
                    INFO(crossing_details);
                    CHECK_FALSE(has_crossing);
                }
            }
        }

        ScopedTemporaryFile gcode_file(".gcode");
        GCodeProcessorResult preview;
        print.export_gcode(gcode_file.string(), &preview, nullptr);

        std::ifstream stream(gcode_file.string());
        REQUIRE(stream.good());
        const std::string gcode(
            (std::istreambuf_iterator<char>(stream)),
            std::istreambuf_iterator<char>());
        CHECK(gcode.find(";TYPE:Arc bridge") != std::string::npos);
        CHECK(gcode.find(";TYPE:Arc overhang") == std::string::npos);

        struct PrintedSegment {
            Vec2d a;
            Vec2d b;
            double width;
            unsigned int layer_id;
            ExtrusionRole role;
            double print_z;
        };
        std::vector<PrintedSegment> printed_segments;
        size_t arc_path_count = 0;
        size_t arc_segment_count = 0;
        std::map<unsigned int, Polylines> exported_arc_paths;
        std::map<unsigned int, Polylines> exported_perimeter_paths;

        for (size_t move_idx = 1; move_idx < preview.moves.size(); ++move_idx) {
            const auto &previous = preview.moves[move_idx - 1];
            const auto &move = preview.moves[move_idx];
            const bool is_arc =
                move.type == EMoveType::Extrude &&
                is_arc_fill(move.extrusion_role);
            const bool previous_is_arc =
                previous.type == EMoveType::Extrude &&
                is_arc_fill(previous.extrusion_role) &&
                previous.layer_id == move.layer_id;
            const bool is_wall =
                move.type == EMoveType::Extrude &&
                is_perimeter(move.extrusion_role);
            const bool previous_is_wall =
                previous.type == EMoveType::Extrude &&
                is_perimeter(previous.extrusion_role) &&
                previous.layer_id == move.layer_id;

            if (is_arc) {
                ++arc_segment_count;
                if (!previous_is_arc) {
                    ++arc_path_count;
                    const Vec2d start =
                        previous.position.head<2>().cast<double>();
                    exported_arc_paths[move.layer_id].emplace_back();
                    exported_arc_paths[move.layer_id].back().append(
                        Point::new_scale(start.x(), start.y()));
                    double bead_clearance =
                        std::numeric_limits<double>::max();
                    double any_previous_clearance =
                        std::numeric_limits<double>::max();
                    double previous_layer_z =
                        std::numeric_limits<double>::lowest();
                    const double current_z = move.position.z();
                    for (const PrintedSegment &anchor : printed_segments)
                        if (anchor.print_z < current_z - 0.001)
                            previous_layer_z =
                                std::max(previous_layer_z, anchor.print_z);
                    for (const PrintedSegment &anchor : printed_segments) {
                        const Vec2d delta = anchor.b - anchor.a;
                        const double length_squared = delta.squaredNorm();
                        const double position = length_squared == 0. ? 0. :
                            std::clamp((start - anchor.a).dot(delta) /
                                           length_squared,
                                       0., 1.);
                        const double centerline_distance =
                            (start - anchor.a - position * delta).norm();
                        any_previous_clearance = std::min(
                            any_previous_clearance,
                            centerline_distance -
                                0.5 * (move.width + anchor.width));
                        const bool same_layer_anchor =
                            std::abs(anchor.print_z - current_z) <= 0.001 &&
                            (is_perimeter(anchor.role) ||
                             is_arc_fill(anchor.role));
                        const bool lower_layer_anchor =
                            std::abs(anchor.print_z - previous_layer_z) <= 0.001;
                        if (!same_layer_anchor && !lower_layer_anchor)
                            continue;
                        bead_clearance = std::min(
                            bead_clearance,
                            centerline_distance -
                                0.5 * (move.width + anchor.width));
                    }
                    INFO("arc start " << start.x() << "," << start.y()
                                      << " layer " << move.layer_id
                                      << " clearance " << bead_clearance
                                      << " any " << any_previous_clearance);
                    CHECK(bead_clearance <= 0.051);
                }
                const Vec2d end = move.position.head<2>().cast<double>();
                exported_arc_paths[move.layer_id].back().append(
                    Point::new_scale(end.x(), end.y()));
            }

            if (is_wall) {
                if (!previous_is_wall) {
                    const Vec2d start =
                        previous.position.head<2>().cast<double>();
                    exported_perimeter_paths[move.layer_id].emplace_back();
                    exported_perimeter_paths[move.layer_id].back().append(
                        Point::new_scale(start.x(), start.y()));
                }
                const Vec2d end =
                    move.position.head<2>().cast<double>();
                exported_perimeter_paths[move.layer_id].back().append(
                    Point::new_scale(end.x(), end.y()));
            }

            if (move.type == EMoveType::Extrude) {
                printed_segments.push_back({
                    previous.position.head<2>().cast<double>(),
                    move.position.head<2>().cast<double>(),
                    move.width,
                    move.layer_id,
                    move.extrusion_role,
                    move.position.z()
                });
            }
        }
        for (const auto &[layer_id, paths] : exported_arc_paths) {
            INFO("exported arc-overhang layer " << layer_id);
            std::string crossing_details;
            const bool has_crossing =
                arc_paths_have_proper_crossing(paths, &crossing_details);
            INFO(crossing_details);
            CHECK_FALSE(has_crossing);
            const auto walls = exported_perimeter_paths.find(layer_id);
            if (walls != exported_perimeter_paths.end()) {
                crossing_details.clear();
                const bool crosses_wall = arc_paths_cross_obstacles(
                    paths, walls->second, &crossing_details);
                INFO(crossing_details);
                CHECK_FALSE(crosses_wall);
                crossing_details.clear();
                const bool tracks_wall = arc_paths_track_obstacles(
                    paths, walls->second, scale_(0.25), scale_(0.45),
                    &crossing_details);
                INFO(crossing_details);
                CHECK_FALSE(tracks_wall);
            }
            // This validation model requests two walls. The inner one is
            // replaced on the arc layer, exposing the previous layer as the
            // supported runway while the external dimensional shell remains.
            const unsigned int current_layer_id = layer_id;
            CHECK_FALSE(std::any_of(
                preview.moves.begin(), preview.moves.end(),
                [current_layer_id](const auto &move) {
                    return move.layer_id == current_layer_id &&
                           move.type == EMoveType::Extrude &&
                           move.extrusion_role == erPerimeter;
                }));
        }

        bool arc_fill_role = false;
        double x = 0.;
        double y = 0.;
        size_t arc_command_count = 0;
        std::istringstream lines(gcode);
        std::string line;
        while (std::getline(lines, line)) {
            if (line.rfind(";TYPE:", 0) == 0) {
                arc_fill_role = line == ";TYPE:Arc bridge";
                continue;
            }
            std::istringstream words(line);
            std::string command;
            words >> command;
            const bool is_motion =
                command == "G0" || command == "G1" ||
                command == "G2" || command == "G3";
            if (!is_motion)
                continue;

            double next_x = x;
            double next_y = y;
            double i = 0.;
            double j = 0.;
            bool has_i = false;
            bool has_j = false;
            std::string word;
            while (words >> word) {
                if (word.size() < 2)
                    continue;
                try {
                    const double value = std::stod(word.substr(1));
                    switch (word.front()) {
                    case 'X': next_x = value; break;
                    case 'Y': next_y = value; break;
                    case 'I': i = value; has_i = true; break;
                    case 'J': j = value; has_j = true; break;
                    default: break;
                    }
                } catch (const std::exception &) {
                    // The remainder may be a comment rather than a G-code word.
                    break;
                }
            }

            if (arc_fill_role &&
                (command == "G2" || command == "G3") &&
                (has_i || has_j)) {
                ++arc_command_count;
            }
            x = next_x;
            y = next_y;
        }
        CHECK(arc_path_count >= 4);
        CHECK(arc_segment_count > arc_path_count);
        // Arc-overhang geometry is fitted only after its paths have passed the
        // generator's anchoring and collision checks. The parsed preview checks
        // above exercise the emitted curves rather than just their source
        // polylines, so fitted G2/G3 moves must retain those guarantees.
        CHECK(arc_command_count > 0);
    }
}

TEST_CASE("Arc overhang fitting follows the global arc fitting setting",
          "[Fill][ArcOverhang][GCode][ArcFitting]")
{
    const bool enable_arc_fitting = GENERATE(false, true);
    DYNAMIC_SECTION("arc fitting " << enable_arc_fitting)
    {
        Print print;
        Slic3r::Test::init_and_process_print(
            {make_arc_bridge_validation_model(4.5, 30.)}, print, {
                {"layer_height", "0.2"},
                {"initial_layer_print_height", "0.2"},
                {"wall_loops", "2"},
                {"arc_overhang_enabled", "1"},
                {"arc_overhang_bridge_distance", "0"},
                {"arc_overhang_min_overhang_distance", "0"},
                {"arc_overhang_recursive_fill", "1"},
                {"enable_arc_fitting", enable_arc_fitting ? "1" : "0"},
                {"resolution", "0"},
                {"enable_support", "0"}
            });
        REQUIRE(print_has_arc_overhang(print));

        ScopedTemporaryFile gcode_file(".gcode");
        GCodeProcessorResult preview;
        print.export_gcode(gcode_file.string(), &preview, nullptr);

        std::ifstream stream(gcode_file.string());
        REQUIRE(stream.good());
        bool arc_fill_role = false;
        size_t arc_command_count = 0;
        for (std::string line; std::getline(stream, line);) {
            if (line.rfind(";TYPE:", 0) == 0) {
                arc_fill_role = line == ";TYPE:Arc bridge";
                continue;
            }
            if (!arc_fill_role)
                continue;
            std::istringstream words(line);
            std::string command;
            words >> command;
            if (command == "G2" || command == "G3")
                ++arc_command_count;
        }
        if (enable_arc_fitting)
            CHECK(arc_command_count > 0);
        else
            CHECK(arc_command_count == 0);
    }
}

TEST_CASE("Arc overhang keeps perimeter width while scaling volume", "[Fill][ArcOverhang][Flow]")
{
    TriangleMesh model = make_cube(10., 20., 5.);
    TriangleMesh roof  = make_cube(30., 20., 2.);
    roof.translate(-10., 0., 5.);
    model.merge(roof);

    constexpr double layer_height = 0.2;
    constexpr double perimeter_width = 0.45;
    constexpr double arc_flow_ratio = 0.8;
    const bool thick_bridges = GENERATE(false, true);
    Print print;
    Slic3r::Test::init_and_process_print({model}, print, {
        {"layer_height", std::to_string(layer_height)},
        {"inner_wall_line_width", std::to_string(perimeter_width)},
        {"bridge_flow", "1"},
        {"thick_bridges", thick_bridges ? "1" : "0"},
        {"arc_overhang_flow_ratio", std::to_string(arc_flow_ratio * 100.) + "%"},
        {"arc_overhang_enabled", "1"},
        {"arc_overhang_bridge_distance", "0"},
        {"arc_overhang_min_overhang_distance", "0"},
        {"enable_support", "0"}
    });

    const std::optional<ExtrusionPath> path = first_arc_overhang_path(print);
    REQUIRE(path.has_value());
    const double expected_area =
        PI * perimeter_width * perimeter_width * 0.25 * arc_flow_ratio;
    CHECK_THAT(path->mm3_per_mm, Catch::Matchers::WithinRel(expected_area, 1e-5));
    CHECK_THAT(path->width, Catch::Matchers::WithinRel(perimeter_width, 1e-5));
    CHECK(path->width >= path->height);
    CHECK_THAT(path->height, Catch::Matchers::WithinRel(perimeter_width, 1e-5));
    CHECK(path->height > layer_height);
}

TEST_CASE("Extra overhang perimeters do not consume arc overhang roofs", "[Fill][ArcOverhang][Regression]")
{
    TriangleMesh model = make_cube(10., 20., 5.);
    TriangleMesh roof  = make_cube(30., 20., 2.);
    roof.translate(-10., 0., 5.);
    model.merge(roof);

    Print print;
    Slic3r::Test::init_and_process_print({model}, print, {
        {"arc_overhang_enabled", "1"},
        {"arc_overhang_bridge_distance", "5"},
        {"extra_perimeters_on_overhangs", "1"},
        {"wall_loops", "3"},
        {"enable_support", "0"}
    });

    CHECK(print_has_arc_overhang(print));
}

TEST_CASE("Arc bridge and overhang selectors can be disabled independently", "[Fill][ArcOverhang][Config]")
{
    TriangleMesh model = make_cube(5., 20., 5.);
    TriangleMesh right_support = make_cube(5., 20., 5.);
    right_support.translate(25., 0., 0.);
    model.merge(right_support);
    TriangleMesh roof  = make_cube(30., 20., 2.);
    roof.translate(0., 0., 5.);
    model.merge(roof);

    Print print;
    Slic3r::Test::init_and_process_print({model}, print, {
        {"arc_overhang_enabled", "1"},
        {"arc_overhang_bridges", "0"},
        {"arc_overhang_overhangs", "1"},
        {"arc_overhang_bridge_distance", "0"},
        {"enable_support", "0"}
    });

    CHECK_FALSE(print_has_arc_overhang(print));
}

TEST_CASE("Generated support contact keeps supported bottom surfaces out of arc overhang fill",
          "[Fill][ArcOverhang][Support][Regression]")
{
    TriangleMesh post = make_cube(5., 20., 5.);
    TriangleMesh roof = make_cube(30., 20., 2.);
    roof.translate(-12.5, 0., 5.);
    post.merge(roof);

    Print print;
    Slic3r::Test::init_and_process_print({post}, print, {
        {"arc_overhang_enabled", "1"},
        {"arc_overhang_bridge_distance", "0"},
        {"arc_overhang_min_overhang_distance", "0"},
        {"enable_support", "1"},
        {"support_type", "normal(auto)"},
        {"support_threshold_angle", "30"},
        {"support_top_z_distance", "0"}
    });

    CHECK_FALSE(print_has_arc_overhang(print));
}

TEST_CASE("A middle support applies the arc threshold to each unsupported bridge span",
          "[Fill][ArcOverhang][Support][Regression]")
{
    TriangleMesh left_post = make_cube(5., 20., 5.);
    TriangleMesh middle_post = make_cube(5., 20., 5.);
    middle_post.translate(22.5, 0., 0.);
    left_post.merge(middle_post);
    TriangleMesh right_post = make_cube(5., 20., 5.);
    right_post.translate(45., 0., 0.);
    left_post.merge(right_post);
    TriangleMesh roof = make_cube(50., 20., 2.);
    roof.translate(0., 0., 5.);
    left_post.merge(roof);

    Print print;
    Slic3r::Test::init_and_process_print({left_post}, print, {
        {"arc_overhang_enabled", "1"},
        {"arc_overhang_bridge_distance", "25"},
        {"enable_support", "0"}
    });

    CHECK_FALSE(print_has_arc_overhang(print));
}

#if 0
TEST_CASE("Adjusted solid distance", "[Fill]") {
    int surface_width = 250;
    int distance = Slic3r::Flow::solid_spacing(surface_width, 47);
    REQUIRE(distance == Catch::Approx(50));
    REQUIRE(surface_width % distance == 0);
}
#endif

TEST_CASE("Pattern path length", "[Fill]") {
    std::unique_ptr<Slic3r::Fill> filler(Slic3r::Fill::new_from_type("rectilinear"));
    filler->angle = float(-(PI)/2.0);
	FillParams fill_params;
	filler->spacing = 5;
	fill_params.dont_adjust = true;
	//fill_params.endpoints_overlap = false;
	fill_params.density = float(filler->spacing / 50.0);

    auto test = [&filler, &fill_params] (const ExPolygon& poly) -> Slic3r::Polylines {
        Slic3r::Surface surface(stTop, poly);
        return filler->fill_surface(&surface, fill_params);
    };

    SECTION("Square") {
        Slic3r::Points test_set;
        test_set.reserve(4);
        std::vector<Vec2d> points {Vec2d(0,0), Vec2d(100,0), Vec2d(100,100), Vec2d(0,100)};
        for (size_t i = 0; i < 4; ++i) {
            std::transform(points.cbegin()+i, points.cend(),   std::back_inserter(test_set), [] (const Vec2d& a) -> Point { return Point::new_scale(a.x(), a.y()); } ); 
            std::transform(points.cbegin(), points.cbegin()+i, std::back_inserter(test_set), [] (const Vec2d& a) -> Point { return Point::new_scale(a.x(), a.y()); } );
            Slic3r::Polylines paths = test(Slic3r::ExPolygon(test_set));
            REQUIRE(paths.size() == 1); // one continuous path

            // TODO: determine what the "Expected length" should be for rectilinear fill of a 100x100 polygon. 
            // This check only checks that it's above scale(3*100 + 2*50) + scaled_epsilon.
            // ok abs($paths->[0]->length - scale(3*100 + 2*50)) - scaled_epsilon, 'path has expected length';
            REQUIRE(std::abs(paths[0].length() - static_cast<double>(scale_(3*100 + 2*50))) - SCALED_EPSILON > 0); // path has expected length

            test_set.clear();
        }
    }
    SECTION("Diamond with endpoints on grid") {
        std::vector<Vec2d> points {Vec2d(0,0), Vec2d(100,0), Vec2d(150,50), Vec2d(100,100), Vec2d(0,100), Vec2d(-50,50)};
        Slic3r::Points test_set;
        test_set.reserve(6);
        std::transform(points.cbegin(), points.cend(),   std::back_inserter(test_set), [] (const Vec2d& a) -> Point { return Point::new_scale(a.x(), a.y()); } );
        Slic3r::Polylines paths = test(Slic3r::ExPolygon(test_set));
        REQUIRE(paths.size() == 1); // one continuous path
    }

    SECTION("Square with hole") {
        std::vector<Vec2d> square {Vec2d(0,0), Vec2d(100,0), Vec2d(100,100), Vec2d(0,100)};
        std::vector<Vec2d> hole {Vec2d(25,25), Vec2d(75,25), Vec2d(75,75), Vec2d(25,75) };
        std::reverse(hole.begin(), hole.end());

        Slic3r::Points test_hole;
        Slic3r::Points test_square;

        std::transform(square.cbegin(), square.cend(), std::back_inserter(test_square), [] (const Vec2d& a) -> Point { return Point::new_scale(a.x(), a.y()); } );
        std::transform(hole.cbegin(), hole.cend(), std::back_inserter(test_hole), [] (const Vec2d& a) -> Point { return Point::new_scale(a.x(), a.y()); } );

        for (double angle : {-(PI/2.0), -(PI/4.0), -(PI), PI/2.0, PI}) {
            for (double spacing : {25.0, 5.0, 7.5, 8.5}) {
				fill_params.density = float(filler->spacing / spacing);
                filler->angle = float(angle);
                ExPolygon e(test_square, test_hole);
                Slic3r::Polylines paths = test(e);
#if 0
				{
					BoundingBox bbox = get_extents(e);
					SVG svg("c:\\data\\temp\\square_with_holes.svg", bbox);
					svg.draw(e);
					svg.draw(paths);
					svg.Close();
				}
#endif
                REQUIRE((paths.size() >= 1 && paths.size() <= 3));
                // paths don't cross hole
                REQUIRE(diff_pl(paths, offset(e, float(SCALED_EPSILON*10))).size() == 0);
            }
        }
    }
    SECTION("Regression: Missing infill segments in some rare circumstances") {
        filler->angle = float(PI/4.0);
		fill_params.dont_adjust = false;
        filler->spacing = 0.654498;
        //filler->endpoints_overlap = unscale(359974);
		fill_params.density = 1;
        filler->layer_id = 66;
        filler->z = 20.15;

        Slic3r::Points points {Point(25771516,14142125),Point(14142138,25771515),Point(2512749,14142131),Point(14142125,2512749)};
        Slic3r::Polylines paths = test(Slic3r::ExPolygon(points));
        REQUIRE(paths.size() == 1); // one continuous path

        // TODO: determine what the "Expected length" should be for rectilinear fill of a 100x100 polygon. 
        // This check only checks that it's above scale(3*100 + 2*50) + scaled_epsilon.
        // ok abs($paths->[0]->length - scale(3*100 + 2*50)) - scaled_epsilon, 'path has expected length';
        REQUIRE(std::abs(paths[0].length() - static_cast<double>(scale_(3*100 + 2*50))) - SCALED_EPSILON > 0); // path has expected length
    }

    SECTION("Rotated Square") {
        Slic3r::Points square { Point::new_scale(0,0), Point::new_scale(50,0), Point::new_scale(50,50), Point::new_scale(0,50)};
        Slic3r::ExPolygon expolygon(square);
        std::unique_ptr<Slic3r::Fill> filler(Slic3r::Fill::new_from_type("rectilinear"));
		filler->bounding_box = get_extents(expolygon.contour);
        filler->angle = 0;
        
        Surface surface(stTop, expolygon);
        auto flow = Slic3r::Flow(0.69f, 0.4f, 0.50f);

		FillParams fill_params;
		fill_params.density = 1.0;
		filler->spacing = flow.spacing();

        for (auto angle : { 0.0, 45.0}) {
            surface.expolygon.rotate(angle, Point(0,0));
            Polylines paths = filler->fill_surface(&surface, fill_params);
            REQUIRE(paths.size() == 1);
        }
    }

    #if 0   // Disabled temporarily due to precision issues on the Mac VM
    SECTION("Solid surface fill") {
        Slic3r::Points points {
            Point::new_scale(6883102, 9598327.01296997),
            Point::new_scale(6883102, 20327272.01297),
            Point::new_scale(3116896, 20327272.01297),
            Point::new_scale(3116896, 9598327.01296997) 
        };
        Slic3r::ExPolygon expolygon(points);
         
        REQUIRE(test_if_solid_surface_filled(expolygon, 0.55) == true);
        for (size_t i = 0; i <= 20; ++i)
        {
            expolygon.scale(1.05);
            REQUIRE(test_if_solid_surface_filled(expolygon, 0.55) == true);
        }
    }
    #endif

    SECTION("Solid surface fill") {
        Slic3r::Points points {
                Slic3r::Point(59515297,5422499),Slic3r::Point(59531249,5578697),Slic3r::Point(59695801,6123186),
                Slic3r::Point(59965713,6630228),Slic3r::Point(60328214,7070685),Slic3r::Point(60773285,7434379),
                Slic3r::Point(61274561,7702115),Slic3r::Point(61819378,7866770),Slic3r::Point(62390306,7924789),
                Slic3r::Point(62958700,7866744),Slic3r::Point(63503012,7702244),Slic3r::Point(64007365,7434357),
                Slic3r::Point(64449960,7070398),Slic3r::Point(64809327,6634999),Slic3r::Point(65082143,6123325),
                Slic3r::Point(65245005,5584454),Slic3r::Point(65266967,5422499),Slic3r::Point(66267307,5422499),
                Slic3r::Point(66269190,8310081),Slic3r::Point(66275379,17810072),Slic3r::Point(66277259,20697500),
                Slic3r::Point(65267237,20697500),Slic3r::Point(65245004,20533538),Slic3r::Point(65082082,19994444),
                Slic3r::Point(64811462,19488579),Slic3r::Point(64450624,19048208),Slic3r::Point(64012101,18686514),
                Slic3r::Point(63503122,18415781),Slic3r::Point(62959151,18251378),Slic3r::Point(62453416,18198442),
                Slic3r::Point(62390147,18197355),Slic3r::Point(62200087,18200576),Slic3r::Point(61813519,18252990),
                Slic3r::Point(61274433,18415918),Slic3r::Point(60768598,18686517),Slic3r::Point(60327567,19047892),
                Slic3r::Point(59963609,19493297),Slic3r::Point(59695865,19994587),Slic3r::Point(59531222,20539379),
                Slic3r::Point(59515153,20697500),Slic3r::Point(58502480,20697500),Slic3r::Point(58502480,5422499)
        };
        Slic3r::ExPolygon expolygon(points);
         
        REQUIRE(test_if_solid_surface_filled(expolygon, 0.55) == true);
        REQUIRE(test_if_solid_surface_filled(expolygon, 0.55, PI/2.0) == true);
    }
    SECTION("Solid surface fill") {
        Slic3r::Points points {
            Point::new_scale(0,0),Point::new_scale(98,0),Point::new_scale(98,10), Point::new_scale(0,10)
        };
        Slic3r::ExPolygon expolygon(points);
         
        REQUIRE(test_if_solid_surface_filled(expolygon, 0.5, 45.0, 0.99) == true);
    }
}

/*
{
    my $collection = Slic3r::Polyline::Collection->new(
            Slic3r::Polyline->new([0,15], [0,18], [0,20]),
            Slic3r::Polyline->new([0,10], [0,8], [0,5]),
            );
    is_deeply
        [ map $_->[Y], map @$_, @{$collection->chained_path_from(Slic3r::Point->new(0,30), 0)} ],
        [20, 18, 15, 10, 8, 5],
        'chained path';
}

{
    my $collection = Slic3r::Polyline::Collection->new(
            Slic3r::Polyline->new([4,0], [10,0], [15,0]),
            Slic3r::Polyline->new([10,5], [15,5], [20,5]),
            );
    is_deeply
        [ map $_->[X], map @$_, @{$collection->chained_path_from(Slic3r::Point->new(30,0), 0)} ],
        [reverse 4, 10, 15, 10, 15, 20],
        'chained path';
}

{
    my $collection = Slic3r::ExtrusionPath::Collection->new(
            map Slic3r::ExtrusionPath->new(polyline => $_, role => 0, mm3_per_mm => 1),
            Slic3r::Polyline->new([0,15], [0,18], [0,20]),
            Slic3r::Polyline->new([0,10], [0,8], [0,5]),
            );
    is_deeply
        [ map $_->[Y], map @{$_->polyline}, @{$collection->chained_path_from(Slic3r::Point->new(0,30), 0)} ],
        [20, 18, 15, 10, 8, 5],
        'chained path';
}

{
    my $collection = Slic3r::ExtrusionPath::Collection->new(
            map Slic3r::ExtrusionPath->new(polyline => $_, role => 0, mm3_per_mm => 1),
            Slic3r::Polyline->new([15,0], [10,0], [4,0]),
            Slic3r::Polyline->new([10,5], [15,5], [20,5]),
            );
    is_deeply
        [ map $_->[X], map @{$_->polyline}, @{$collection->chained_path_from(Slic3r::Point->new(30,0), 0)} ],
        [reverse 4, 10, 15, 10, 15, 20],
        'chained path';
}

for my $pattern (qw(rectilinear honeycomb hilbertcurve concentric)) {
    my $config = Slic3r::Config->new_from_defaults;
    $config->set('fill_pattern', $pattern);
    $config->set('external_fill_pattern', $pattern);
    $config->set('perimeters', 1);
    $config->set('skirts', 0);
    $config->set('fill_density', 20);
    $config->set('layer_height', 0.05);
    $config->set('perimeter_extruder', 1);
    $config->set('infill_extruder', 2);
    my $print = Slic3r::Test::init_print('20mm_cube', config => $config, scale => 2);
    ok my $gcode = Slic3r::Test::gcode($print), "successful $pattern infill generation";
    my $tool = undef;
    my @perimeter_points = my @infill_points = ();
    Slic3r::GCode::Reader->new->parse($gcode, sub {
            my ($self, $cmd, $args, $info) = @_;

            if ($cmd =~ /^T(\d+)/) {
            $tool = $1;
            } elsif ($cmd eq 'G1' && $info->{extruding} && $info->{dist_XY} > 0) {
            if ($tool == $config->perimeter_extruder-1) {
            push @perimeter_points, Slic3r::Point->new_scale($args->{X}, $args->{Y});
            } elsif ($tool == $config->infill_extruder-1) {
            push @infill_points, Slic3r::Point->new_scale($args->{X}, $args->{Y});
            }
            }
            });
    my $convex_hull = convex_hull(\@perimeter_points);
    ok !(defined first { !$convex_hull->contains_point($_) } @infill_points), "infill does not exceed perimeters ($pattern)";
}

{
    my $config = Slic3r::Config->new_from_defaults;
    $config->set('infill_only_where_needed', 1);
    $config->set('bottom_solid_layers', 0);
    $config->set('infill_extruder', 2);
    $config->set('infill_extrusion_width', 0.5);
    $config->set('fill_density', 40);
    $config->set('cooling', 0);                 # for preventing speeds from being altered
        $config->set('first_layer_speed', '100%');  # for preventing speeds from being altered

        my $test = sub {
            my $print = Slic3r::Test::init_print('pyramid', config => $config);

            my $tool = undef;
            my @infill_extrusions = ();  # array of polylines
                Slic3r::GCode::Reader->new->parse(Slic3r::Test::gcode($print), sub {
                        my ($self, $cmd, $args, $info) = @_;

                        if ($cmd =~ /^T(\d+)/) {
                        $tool = $1;
                        } elsif ($cmd eq 'G1' && $info->{extruding} && $info->{dist_XY} > 0) {
                        if ($tool == $config->infill_extruder-1) {
                        push @infill_extrusions, Slic3r::Line->new_scale(
                                [ $self->X, $self->Y ],
                                [ $info->{new_X}, $info->{new_Y} ],
                                );
                        }
                        }
                        });
            return 0 if !@infill_extrusions;  # prevent calling convex_hull() with no points

                my $convex_hull = convex_hull([ map $_->pp, map @$_, @infill_extrusions ]);
            return unscale unscale sum(map $_->area, @{offset([$convex_hull], scale(+$config->infill_extrusion_width/2))});
        };

    my $tolerance = 5;  # mm^2

        $config->set('solid_infill_below_area', 0);
    ok $test->() < $tolerance,
       'no infill is generated when using infill_only_where_needed on a pyramid';

    $config->set('solid_infill_below_area', 70);
    ok abs($test->() - $config->solid_infill_below_area) < $tolerance,
       'infill is only generated under the forced solid shells';
}

{
    my $config = Slic3r::Config->new_from_defaults;
    $config->set('skirts', 0);
    $config->set('perimeters', 1);
    $config->set('fill_density', 0);
    $config->set('top_solid_layers', 0);
    $config->set('bottom_solid_layers', 0);
    $config->set('solid_infill_below_area', 20000000);
    $config->set('solid_infill_every_layers', 2);
    $config->set('perimeter_speed', 99);
    $config->set('external_perimeter_speed', 99);
    $config->set('cooling', 0);
    $config->set('first_layer_speed', '100%');

    my $print = Slic3r::Test::init_print('20mm_cube', config => $config);
    my %layers_with_extrusion = ();
    Slic3r::GCode::Reader->new->parse(Slic3r::Test::gcode($print), sub {
            my ($self, $cmd, $args, $info) = @_;

            if ($cmd eq 'G1' && $info->{dist_XY} > 0 && $info->{extruding}) {
            if (($args->{F} // $self->F) != $config->perimeter_speed*60) {
            $layers_with_extrusion{$self->Z} = ($args->{F} // $self->F);
            }
            }
            });

    ok !%layers_with_extrusion,
       "solid_infill_below_area and solid_infill_every_layers are ignored when fill_density is 0";
}

{
    my $config = Slic3r::Config->new_from_defaults;
    $config->set('skirts', 0);
    $config->set('perimeters', 3);
    $config->set('fill_density', 0);
    $config->set('layer_height', 0.2);
    $config->set('first_layer_height', 0.2);
    $config->set('nozzle_diameter', [0.35]);
    $config->set('infill_extruder', 2);
    $config->set('solid_infill_extruder', 2);
    $config->set('infill_extrusion_width', 0.52);
    $config->set('solid_infill_extrusion_width', 0.52);
    $config->set('first_layer_extrusion_width', 0);

    my $print = Slic3r::Test::init_print('A', config => $config);
    my %infill = ();  # Z => [ Line, Line ... ]
        my $tool = undef;
    Slic3r::GCode::Reader->new->parse(Slic3r::Test::gcode($print), sub {
            my ($self, $cmd, $args, $info) = @_;

            if ($cmd =~ /^T(\d+)/) {
            $tool = $1;
            } elsif ($cmd eq 'G1' && $info->{extruding} && $info->{dist_XY} > 0) {
            if ($tool == $config->infill_extruder-1) {
            my $z = 1 * $self->Z;
            $infill{$z} ||= [];
            push @{$infill{$z}}, Slic3r::Line->new_scale(
                    [ $self->X, $self->Y ],
                    [ $info->{new_X}, $info->{new_Y} ],
                    );
            }
            }
            });
    my $grow_d = scale($config->infill_extrusion_width)/2;
    my $layer0_infill = union([ map @{$_->grow($grow_d)}, @{ $infill{0.2} } ]);
    my $layer1_infill = union([ map @{$_->grow($grow_d)}, @{ $infill{0.4} } ]);
    my $diff = diff($layer0_infill, $layer1_infill);
    $diff = offset2_ex($diff, -$grow_d, +$grow_d);
    $diff = [ grep { $_->area > 2*(($grow_d*2)**2) } @$diff ];
    is scalar(@$diff), 0, 'no missing parts in solid shell when fill_density is 0';
}

{
    # GH: #2697
    my $config = Slic3r::Config->new_from_defaults;
    $config->set('perimeter_extrusion_width', 0.72);
    $config->set('top_infill_extrusion_width', 0.1);
    $config->set('infill_extruder', 2);         # in order to distinguish infill
        $config->set('solid_infill_extruder', 2);   # in order to distinguish infill

        my $print = Slic3r::Test::init_print('20mm_cube', config => $config);
    my %infill = ();  # Z => [ Line, Line ... ]
        my %other  = ();  # Z => [ Line, Line ... ]
        my $tool = undef;
    Slic3r::GCode::Reader->new->parse(Slic3r::Test::gcode($print), sub {
            my ($self, $cmd, $args, $info) = @_;

            if ($cmd =~ /^T(\d+)/) {
            $tool = $1;
            } elsif ($cmd eq 'G1' && $info->{extruding} && $info->{dist_XY} > 0) {
            my $z = 1 * $self->Z;
            my $line = Slic3r::Line->new_scale(
                    [ $self->X, $self->Y ],
                    [ $info->{new_X}, $info->{new_Y} ],
                    );
            if ($tool == $config->infill_extruder-1) {
            $infill{$z} //= [];
            push @{$infill{$z}}, $line;
            } else {
            $other{$z} //= [];
            push @{$other{$z}}, $line;
            }
            }
            });
    my $top_z = max(keys %infill);
    my $top_infill_grow_d = scale($config->top_infill_extrusion_width)/2;
    my $top_infill = union([ map @{$_->grow($top_infill_grow_d)}, @{ $infill{$top_z} } ]);
    my $perimeters_grow_d = scale($config->perimeter_extrusion_width)/2;
    my $perimeters = union([ map @{$_->grow($perimeters_grow_d)}, @{ $other{$top_z} } ]);
    my $covered = union_ex([ @$top_infill, @$perimeters ]);
    my @holes = map @{$_->holes}, @$covered;
    ok sum(map unscale unscale $_->area*-1, @holes) < 1, 'no gaps between top solid infill and perimeters';
}
*/

bool test_if_solid_surface_filled(const ExPolygon& expolygon, double flow_spacing, double angle, double density)
{
    std::unique_ptr<Slic3r::Fill> filler(Slic3r::Fill::new_from_type("rectilinear"));
	filler->bounding_box = get_extents(expolygon.contour);
    filler->angle = float(angle);

	Flow flow(float(flow_spacing), 0.4f, float(flow_spacing));
	filler->spacing = flow.spacing();

	FillParams fill_params;
	fill_params.density = float(density);
	fill_params.dont_adjust = false;

	Surface surface(stBottom, expolygon);
	Slic3r::Polylines paths = filler->fill_surface(&surface, fill_params);

    // check whether any part was left uncovered
    Polygons grown_paths;
    grown_paths.reserve(paths.size());

    // figure out what is actually going on here re: data types
    float line_offset = float(scale_(filler->spacing / 2.0 + EPSILON));
    std::for_each(paths.begin(), paths.end(), [line_offset, &grown_paths] (const Slic3r::Polyline& p) {
        polygons_append(grown_paths, offset(p, line_offset));
    });

	// Shrink the initial expolygon a bit, this simulates the infill / perimeter overlap that we usually apply.
    ExPolygons uncovered = diff_ex(offset(expolygon, - float(0.2 * scale_(flow_spacing))), grown_paths, ApplySafetyOffset::Yes);

    // ignore very small dots
    const double scaled_flow_spacing = std::pow(scale_(flow_spacing), 2);
    uncovered.erase(std::remove_if(uncovered.begin(), uncovered.end(), [scaled_flow_spacing](const ExPolygon& poly) { return poly.area() < scaled_flow_spacing; }), uncovered.end());

#if 0
	if (! uncovered.empty()) {
		BoundingBox bbox = get_extents(expolygon.contour);
		bbox.merge(get_extents(uncovered));
		bbox.merge(get_extents(grown_paths));
		SVG svg("c:\\data\\temp\\test_if_solid_surface_filled.svg", bbox);
		svg.draw(expolygon);
		svg.draw(uncovered, "red");
		svg.Close();
	}
#endif

    return uncovered.empty(); // solid surface is fully filled
}

// Length-weighted dominant direction of the layer's role_wanted extrusions, whole degrees
// [0, 180), or -1 if it has none. Needs a line pattern such as monotonic or rectilinear.
template<typename RolePred> static int dominant_fill_angle(const Layer &layer, RolePred role_wanted)
{
    std::map<int, double> weight_per_degree;

    auto account = [&weight_per_degree, &role_wanted](const ExtrusionPath &path) {
        if (!role_wanted(path.role()))
            return;
        const Points3 &pts = path.polyline.points;
        for (size_t i = 1; i < pts.size(); ++i) {
            const double dx = double(pts[i].x() - pts[i - 1].x());
            const double dy = double(pts[i].y() - pts[i - 1].y());
            const double len = std::hypot(dx, dy);
            if (len <= 0.)
                continue;
            int deg = int(std::lround(Geometry::rad2deg(std::atan2(dy, dx)))) % 180;
            if (deg < 0)
                deg += 180;
            weight_per_degree[deg] += len;
        }
    };

    for (const LayerRegion *region : layer.regions())
        for (const ExtrusionEntity *entity : region->fills.flatten().entities) {
            if (auto *path = dynamic_cast<const ExtrusionPath *>(entity))
                account(*path);
            else if (auto *multi = dynamic_cast<const ExtrusionMultiPath *>(entity))
                for (const ExtrusionPath &p : multi->paths)
                    account(p);
            else if (auto *loop = dynamic_cast<const ExtrusionLoop *>(entity))
                for (const ExtrusionPath &p : loop->paths)
                    account(p);
        }

    if (weight_per_degree.empty())
        return -1;
    return std::max_element(weight_per_degree.begin(), weight_per_degree.end(),
                            [](const auto &a, const auto &b) { return a.second < b.second; })->first;
}

template<typename RolePred> static std::vector<int> angles_per_layer(const Print &print, RolePred role_wanted)
{
    std::vector<int> angles;
    for (const Layer *layer : print.objects().front()->layers())
        angles.push_back(dominant_fill_angle(*layer, role_wanted));
    return angles;
}

static bool solid_role(ExtrusionRole role) { return is_solid_infill(role) && role != erIroning; }
static bool sparse_role(ExtrusionRole role) { return role == erInternalInfill; }
static bool ironing_role(ExtrusionRole role) { return role == erIroning; }

TEST_CASE("Infill rotation template is unaffected by a raft", "[Fill][Regression]")
{
    // More angles than raft layers, so a raft cannot alias back to the same angle.
    const std::string template_string = GENERATE("+45", "0,25,50,75,100,125,150");
    const int raft_layers = GENERATE(1, 3);
    CAPTURE(template_string, raft_layers);

    auto angles_for = [&template_string](int rafts) {
        Print print;
        // 100% density makes every layer solid, so the template shows on all 100, not just shells.
        Slic3r::Test::init_and_process_print({Slic3r::Test::cube(20)}, print,
                                            {{"solid_infill_rotate_template", template_string},
                                             {"sparse_infill_density", "100%"},
                                             {"internal_solid_infill_pattern", "monotonic"},
                                             {"layer_height", 0.2},
                                             {"raft_layers", rafts}});
        return angles_per_layer(print, solid_role);
    };

    const std::vector<int> without_raft = angles_for(0);
    const std::vector<int> with_raft    = angles_for(raft_layers);

    REQUIRE(without_raft.size() == 100);
    REQUIRE(with_raft.size() == without_raft.size());
    REQUIRE(std::count(without_raft.begin(), without_raft.end(), -1) == 0);
    CHECK(with_raft == without_raft);
}

TEST_CASE("Sparse infill rotation template turns the infill layer by layer", "[Fill]")
{
    const std::vector<int> expected_cycle = {0, 25, 50, 75, 100, 125, 150};

    Print print;
    // No shells, so every layer is sparse infill rather than solid.
    Slic3r::Test::init_and_process_print({Slic3r::Test::cube(10)}, print,
                                        {{"sparse_infill_rotate_template", "0,25,50,75,100,125,150"},
                                         {"sparse_infill_density", "40%"},
                                         {"sparse_infill_pattern", "rectilinear"},
                                         {"top_shell_layers", 0},
                                         {"bottom_shell_layers", 0},
                                         {"layer_height", 0.2}});

    const std::vector<int> angles = angles_per_layer(print, sparse_role);
    REQUIRE(angles.size() == 50);
    REQUIRE(std::count(angles.begin(), angles.end(), -1) == 0);

    std::vector<int> expected;
    for (size_t i = 0; i < angles.size(); ++i)
        expected.push_back(expected_cycle[i % expected_cycle.size()]);
    CHECK(angles == expected);
}

TEST_CASE("Non-planar infill is a bounded wall-anchored Z wave", "[Fill][NonplanarInfill]")
{
    Print print;
    Slic3r::Test::init_and_process_print(
        {Slic3r::Test::cube(10)}, print,
        {{"nonplanar_infill", 1},
         {"nonplanar_infill_amplitude", 0.08},
         {"nonplanar_infill_wavelength", 6.0},
         {"nonplanar_infill_resolution", 0.5},
         {"sparse_infill_density", "40%"},
         {"sparse_infill_pattern", "rectilinear"},
         {"top_shell_layers", 0},
         {"bottom_shell_layers", 0},
         {"layer_height", 0.2}});

    size_t waved_paths = 0;
    bool has_nonzero_z = false;
    for (const Layer *layer : print.objects().front()->layers())
        for (const LayerRegion *region : layer->regions())
            for (const ExtrusionEntity *entity : region->fills.flatten().entities)
                if (const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
                    path != nullptr && path->nonplanar_infill) {
                    ++waved_paths;
                    REQUIRE(path->z_contoured);
                    REQUIRE(path->polyline.points.size() >= 3);
                    CHECK(path->polyline.points.front().z() == 0);
                    CHECK(path->polyline.points.back().z() == 0);
                    for (const Point3 &point : path->polyline.points) {
                        CHECK(std::abs(unscale_(point.z())) <= 0.08 + EPSILON);
                        has_nonzero_z |= point.z() != 0;
                    }
                }
    CHECK(waved_paths > 0);
    CHECK(has_nonzero_z);
    ScopedTemporaryFile gcode_file(".gcode");
    GCodeProcessorResult preview;
    REQUIRE_NOTHROW(print.export_gcode(gcode_file.string(), &preview, nullptr));
    std::ifstream stream(gcode_file.string());
    REQUIRE(stream.good());
    const std::string gcode((std::istreambuf_iterator<char>(stream)),
                            std::istreambuf_iterator<char>());
    CHECK(gcode.find(";TYPE:Non-planar interlocking infill") != std::string::npos);
    CHECK(std::any_of(preview.moves.begin(), preview.moves.end(),
        [](const GCodeProcessorResult::MoveVertex &move) {
            return move.type == EMoveType::Extrude &&
                   move.extrusion_role == erNonplanarInfill;
        }));
}

TEST_CASE("Fuzzy skin can texture top fill without losing line anchors", "[Fill][FuzzySkin]")
{
    Print print;
    Slic3r::Test::init_and_process_print(
        {Slic3r::Test::cube(10)}, print,
        {{"fuzzy_skin", "all"},
         {"fuzzy_skin_top_surface", 1},
         {"fuzzy_skin_thickness", 0.1},
         {"fuzzy_skin_point_distance", 0.3},
         {"top_surface_pattern", "rectilinear"},
         {"layer_height", 0.2}});

    size_t textured_paths = 0;
    for (const Layer *layer : print.objects().front()->layers())
        for (const LayerRegion *region : layer->regions())
            for (const ExtrusionEntity *entity : region->fills.flatten().entities)
                if (const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
                    path != nullptr && path->role() == erTopSolidInfill &&
                    path->polyline.points.size() > 3) {
                    ++textured_paths;
                    CHECK(path->polyline.points.front().z() == 0);
                    CHECK(path->polyline.points.back().z() == 0);
                }
    CHECK(textured_paths > 0);
}

TEST_CASE("Infill rotation template layer count modifier holds each angle for N layers", "[Fill]")
{
    Print print;
    // "+45#2" turns 45 degrees every 2 layers, so equal angles come in pairs.
    Slic3r::Test::init_and_process_print({Slic3r::Test::cube(10)}, print,
                                        {{"solid_infill_rotate_template", "+45#2"},
                                         {"sparse_infill_density", "100%"},
                                         {"internal_solid_infill_pattern", "monotonic"},
                                         {"layer_height", 0.2}});

    const std::vector<int> angles = angles_per_layer(print, solid_role);
    REQUIRE(angles.size() == 50);
    REQUIRE(std::count(angles.begin(), angles.end(), -1) == 0);

    std::vector<int> run_lengths;
    for (size_t i = 0; i < angles.size();) {
        size_t j = i;
        while (j < angles.size() && angles[j] == angles[i])
            ++j;
        run_lengths.push_back(int(j - i));
        i = j;
    }
    // The first and last runs can be clipped by the start and end of the object.
    REQUIRE(run_lengths.size() > 3);
    const std::vector<int> interior(run_lengths.begin() + 1, run_lengths.end() - 1);
    CHECK(std::count(interior.begin(), interior.end(), 2) == int(interior.size()));
}

TEST_CASE("Z anti-aliasing keeps the infill rotation template's step", "[Fill]")
{
    Print print;
    Slic3r::Test::init_and_process_print({Slic3r::Test::cube(10)}, print,
                                        {{"solid_infill_rotate_template", "+45"},
                                         {"sparse_infill_density", "100%"},
                                         {"internal_solid_infill_pattern", "monotonic"},
                                         {"zaa_enabled", 1},
                                         {"zaa_min_z", 0.05},
                                         {"layer_height", 0.2}});

    // Z contouring varies the layer heights, so the layer count is not 10mm / 0.2mm here.
    const std::vector<int> angles = angles_per_layer(print, solid_role);
    REQUIRE(angles.size() > 10);
    REQUIRE(std::count(angles.begin(), angles.end(), -1) == 0);

    // Z contouring may change when the template advances, but each step must still be 45 degrees.
    int steps = 0;
    for (size_t i = 1; i < angles.size(); ++i) {
        const int delta = ((angles[i] - angles[i - 1]) % 180 + 180) % 180;
        CAPTURE(i, angles[i - 1], angles[i]);
        // Split rather than "delta == 0 || delta == 45" so Catch2 can show the operands.
        REQUIRE(delta % 45 == 0);
        REQUIRE(delta <= 45);
        steps += delta == 45;
    }
    CHECK(steps > 0);
}

TEST_CASE("True non-planar top surfaces work without Z contouring", "[Fill][NonplanarSurface]")
{
    Print print;
    Slic3r::Test::init_and_process_print(
        {shallow_top_wedge(8.)}, print,
        {{"nonplanar_top_surface", 1},
         {"nonplanar_top_surface_max_angle", 45},
         {"nonplanar_top_surface_resolution", 0.2},
         {"nonplanar_top_surface_min_height", 0.08},
         {"zaa_enabled", 0},
         {"layer_height", 0.2}});

    const std::vector<ExtrusionPath> paths = nonplanar_top_paths(print);
    REQUIRE_FALSE(paths.empty());
    const std::vector<ExtrusionPath> boundary_paths = nonplanar_perimeter_paths(print);
    // The complete configured wall stack is draped with the skin. Leaving a
    // nominal-Z inner wall at this boundary reintroduces stair steps and can
    // collide with the top raster. The outer wall retains dimensional
    // accuracy while Arachne closes variable-width wall-to-fill gaps.
    CHECK_FALSE(boundary_paths.empty());
    CHECK(has_connected_nonplanar_wall_loop(print));
    // A completely selected wall has no in-layer planar/non-planar boundary.
    // It therefore remains one closed, continuously contoured loop; forcing
    // its seam down to the nominal owner Z would introduce a notch on every
    // transition course. Attachment is supplied by the aligned flat
    // foundation and the preceding wall course, which are validated below.
    CHECK(std::any_of(boundary_paths.begin(), boundary_paths.end(),
        [](const ExtrusionPath &path) {
            if (path.polyline.points.size() < 3 ||
                path.polyline.points.front().head<2>() !=
                    path.polyline.points.back().head<2>() ||
                path.polyline.points.front().z() !=
                    path.polyline.points.back().z())
                return false;
            const auto [minimum, maximum] = std::minmax_element(
                path.polyline.points.begin(), path.polyline.points.end(),
                [](const Point3 &left, const Point3 &right) {
                    return left.z() < right.z();
                });
            return maximum->z() > minimum->z();
        }));

    // The fill domain is already inset from the generated wall stack. Row
    // phase must not add another half-spacing at either edge, otherwise the
    // skin leaves a visible wedge beside an otherwise valid draped wall.
    // Inspect the interior of the outer non-planar wall and require the
    // deposited wall/skin envelopes to form one closed surface. A 0.05 mm
    // opening removes only polygon rounding noise, not a printable gap.
    Points outer_wall_points;
    Polygons deposited_surface;
    for (const ExtrusionPath &path : boundary_paths) {
        if (path.inset_idx == 0)
            for (const Point3 &point : path.polyline.points)
                outer_wall_points.emplace_back(point.to_point());
        path.polygons_covered_by_width(deposited_surface, 0.f);
    }
    for (const ExtrusionPath &path : paths)
        path.polygons_covered_by_width(deposited_surface, 0.f);
    REQUIRE(outer_wall_points.size() >= 3);
    ExPolygons expected_surface{ExPolygon(
        Geometry::convex_hull(std::move(outer_wall_points)))};
    expected_surface = offset_ex(expected_surface, -scale_(0.05));
    ExPolygons uncovered_surface = diff_ex(
        expected_surface, union_ex(deposited_surface), ApplySafetyOffset::No);
    uncovered_surface = opening_ex(uncovered_surface, scale_(0.05));
    const double uncovered_surface_area = std::accumulate(
        uncovered_surface.begin(), uncovered_surface.end(), 0.,
        [](double total, const ExPolygon &polygon) {
            return total + std::abs(polygon.area()) *
                SCALING_FACTOR * SCALING_FACTOR;
        });
    CAPTURE(uncovered_surface_area);
    CHECK(uncovered_surface_area < 0.05);
    CHECK_FALSE(has_nonplanar_path_with_role(print, erSolidInfill));
    bool has_continuous_z_change = false;
    coord_t maximum_path_z_span = 0;
    size_t x_aligned_segments = 0;
    size_t boundary_anchor_segments = 0;
    size_t short_serpentine_connectors = 0;
    const float maximum_path_width = std::max_element(
        paths.begin(), paths.end(), [](const ExtrusionPath &left,
                                       const ExtrusionPath &right) {
            return left.width < right.width;
        })->width;
    for (const ExtrusionPath &path : paths) {
        REQUIRE(path.role() == erTopSolidInfill);
        REQUIRE(path.z_contoured);
        for (size_t point_idx = 1; point_idx < path.polyline.points.size(); ++point_idx) {
            has_continuous_z_change |= path.polyline.points[point_idx - 1].z() !=
                                       path.polyline.points[point_idx].z();
            const Vec3crd delta = path.polyline.points[point_idx] -
                                  path.polyline.points[point_idx - 1];
            if (std::abs(delta.x()) + std::abs(delta.y()) > scale_(0.05)) {
                if (std::abs(delta.y()) <= scale_(0.001)) {
                    ++x_aligned_segments;
                } else if (std::abs(delta.x()) <= scale_(0.001)) {
                    // Arachne may add a long variable-width anchor along the
                    // rectangular wall band.  It closes the bead-envelope
                    // sliver at the wall rather than connecting unrelated
                    // scanlines, so its length is not bounded by raster
                    // spacing.
                    ++boundary_anchor_segments;
                } else {
                    CHECK(delta.head<2>().cast<double>().norm() <= scale_(0.6));
                    ++short_serpentine_connectors;
                }
            }
        }
        const auto [minimum_z, maximum_z] = std::minmax_element(
            path.polyline.points.begin(), path.polyline.points.end(),
            [](const Point3 &lhs, const Point3 &rhs) { return lhs.z() < rhs.z(); });
        maximum_path_z_span = std::max(maximum_path_z_span, maximum_z->z() - minimum_z->z());
    }
    CHECK(has_continuous_z_change);
    CHECK(x_aligned_segments > 0);
    CHECK(boundary_anchor_segments > 0);
    CHECK(boundary_anchor_segments + short_serpentine_connectors > 0);
    // Every point of the finishing course lies on the source plane. Sampling
    // an interior XY by first projecting it to the patch boundary produces a
    // visibly flat plateau even though the emitted XY remains in the
    // interior. The unknown object translation affects only the intercept,
    // so z - slope*x must remain constant across the complete finishing skin.
    const double expected_slope = std::tan(8. * M_PI / 180.);
    double minimum_plane_intercept = std::numeric_limits<double>::max();
    double maximum_plane_intercept = std::numeric_limits<double>::lowest();
    size_t finishing_surface_points = 0;
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions())
                for (const ExtrusionEntity *entity : region->fills.flatten().entities)
                    if (const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
                        path != nullptr && path->nonplanar_surface &&
                        !path->nonplanar_transition &&
                        path->role() == erTopSolidInfill)
                        for (const Point3 &point : path->polyline.points) {
                            const double intercept = layer->print_z +
                                unscale_(point.z()) -
                                expected_slope * unscale_(point.x());
                            minimum_plane_intercept = std::min(
                                minimum_plane_intercept, intercept);
                            maximum_plane_intercept = std::max(
                                maximum_plane_intercept, intercept);
                            ++finishing_surface_points;
                        }
    REQUIRE(finishing_surface_points > 0);
    CHECK(maximum_plane_intercept - minimum_plane_intercept < 0.01);
    // The boundary is the original slicer-owned wall loop, not a second
    // Arachne outline generated from the finishing patch. This keeps its
    // dimensional contour, configured inset count and seam topology.
    CHECK(std::any_of(boundary_paths.begin(), boundary_paths.end(),
        [](const ExtrusionPath &path) {
            return path.nonplanar_clearance_validated;
        }));
    // The final curved skin is generated from one projected surface rather
    // than nominal-layer fragments, so a raster may cross several ordinary
    // layer heights continuously. The stepped planar shell below supports it.
    CHECK(unscale<double>(maximum_path_z_span) > 0.2 + EPSILON);

    ScopedTemporaryFile gcode_file(".gcode");
    GCodeProcessorResult preview;
    REQUIRE_NOTHROW(print.export_gcode(gcode_file.string(), &preview, nullptr));
    std::ifstream stream(gcode_file.string());
    REQUIRE(stream.good());
    const std::string gcode((std::istreambuf_iterator<char>(stream)),
                            std::istreambuf_iterator<char>());
    CHECK_FALSE(gcode.empty());
    CHECK(gcode.find(";TYPE:Non-planar top surface") != std::string::npos);
    CHECK(gcode.find(";TYPE:Non-planar transition") != std::string::npos);
    // Adjacent variable-width pieces are one wall even when Arachne represents
    // them as separate ExtrusionPaths. Their intentional bead contact must not
    // be mistaken for a travel collision. Genuine moves between disconnected
    // rasters may still need clearance, so reject only a raise whose eventual
    // XYZ target was already within the normal geometry join tolerance.
    Vec3d gcode_position = Vec3d::Zero();
    bool have_gcode_position = false;
    std::optional<Vec3d> position_before_clearance_raise;
    size_t spurious_clearance_raises = 0;
    std::istringstream clearance_lines(gcode);
    for (std::string line; std::getline(clearance_lines, line);) {
        if (line.rfind("G0 ", 0) != 0 && line.rfind("G1 ", 0) != 0)
            continue;
        Vec3d next_position = gcode_position;
        bool has_coordinate = false;
        std::istringstream words(line);
        for (std::string word; words >> word;) {
            if (word.size() < 2)
                continue;
            const int axis = word.front() == 'X' ? 0 :
                             word.front() == 'Y' ? 1 :
                             word.front() == 'Z' ? 2 : -1;
            if (axis >= 0) {
                next_position[axis] = std::stod(word.substr(1));
                has_coordinate = true;
            }
        }
        if (line.find("raise for non-planar toolhead clearance") !=
                std::string::npos && have_gcode_position)
            position_before_clearance_raise = gcode_position;
        if (line.find("set Z for variable-Z surface") != std::string::npos &&
            position_before_clearance_raise) {
            const Vec3d delta = next_position - *position_before_clearance_raise;
            spurious_clearance_raises +=
                delta.head<2>().norm() <= 0.05 && std::abs(delta.z()) <= 0.05;
            position_before_clearance_raise.reset();
        }
        if (has_coordinate) {
            gcode_position = next_position;
            have_gcode_position = true;
        }
    }
    CHECK(spurious_clearance_raises == 0);
    // A mixed wall must be written as one extrusion entity. Splitting it at
    // the planar/non-planar metadata boundary produces the characteristic
    // retract, clearance travel and isolated raised contour that this feature
    // is intended to replace.
    bool wall_started = false;
    bool retracted_after_wall = false;
    bool continuous_wall_transition = false;
    bool nonplanar_started = false;
    std::istringstream gcode_lines(gcode);
    for (std::string line; std::getline(gcode_lines, line);) {
        if (line == ";LAYER_CHANGE") {
            nonplanar_started = false;
            wall_started = false;
            retracted_after_wall = false;
        } else if (line == ";TYPE:Outer wall" || line == ";TYPE:Inner wall" ||
                   line == ";TYPE:Brick wall") {
            wall_started = true;
            retracted_after_wall = false;
        } else if (wall_started && line.find("; retract") != std::string::npos) {
            retracted_after_wall = true;
        } else if (line.rfind(";TYPE:Non-planar", 0) == 0) {
            continuous_wall_transition |= !retracted_after_wall;
            nonplanar_started = true;
        } else if (line.rfind(";TYPE:", 0) == 0 && line != ";TYPE:Custom") {
            CAPTURE(line);
            CHECK_FALSE(nonplanar_started);
        }
    }
    CHECK(continuous_wall_transition);

    std::map<unsigned int, std::pair<float, float>> preview_z_ranges;
    for (const GCodeProcessorResult::MoveVertex &move : preview.moves) {
        if (move.type != EMoveType::Extrude || move.extrusion_role != erNonplanarSurface)
            continue;
        auto range = preview_z_ranges.try_emplace(
            move.layer_id, std::numeric_limits<float>::max(), std::numeric_limits<float>::lowest()).first;
        range->second.first  = std::min(range->second.first, move.position.z());
        range->second.second = std::max(range->second.second, move.position.z());
    }
    REQUIRE_FALSE(preview_z_ranges.empty());
    CHECK(std::any_of(preview_z_ranges.begin(), preview_z_ranges.end(), [](const auto &entry) {
        return entry.second.second - entry.second.first > 0.001f;
    }));
}

TEST_CASE("A connected non-planar top surface transitions across upper layer operations", "[Fill][NonplanarSurface]")
{
    Print print;
    Slic3r::Test::init_and_process_print(
        {shallow_top_wedge(8.)}, print,
        {{"nonplanar_top_surface", 1},
         {"nonplanar_top_surface_max_angle", 45},
         {"nonplanar_top_surface_layers", 5},
         {"zaa_enabled", 0},
         {"layer_height", 0.2}});

    size_t owning_layers = 0;
    coord_t deepest_relative_z = 0;
    size_t planar_anchored_paths = 0;
    size_t unsupported_anchor_points = 0;
    size_t planar_wall_anchor_paths = 0;
    size_t unsupported_wall_anchor_points = 0;
    std::map<std::pair<coord_t, coord_t>, std::vector<double>> shell_z_by_xy;
    std::set<const Layer *> nonplanar_top_owner_layers;
    std::set<const Layer *> nonplanar_wall_owner_layers;
    for (const PrintObject *object : print.objects()) {
        for (const Layer *layer : object->layers()) {
            bool owns_nonplanar_top = false;
            for (const LayerRegion *region : layer->regions()) {
                const auto inspect = [&](const ExtrusionEntityCollection &collection) {
                    for (const ExtrusionEntity *entity : collection.flatten().entities) {
                        const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
                        if (path == nullptr || !path->nonplanar_surface)
                            continue;
                        if (path->role() != erTopSolidInfill)
                            continue;
                        owns_nonplanar_top = true;
                        nonplanar_top_owner_layers.insert(layer);
                        for (const Point3 &point : path->polyline.points) {
                            deepest_relative_z = std::min(deepest_relative_z, point.z());
                            shell_z_by_xy[{point.x(), point.y()}].push_back(
                                layer->print_z + unscale_(point.z()));
                        }

                        if (!path->nonplanar_transition ||
                            path->nonplanar_feature_transition)
                            continue;
                        std::vector<const Layer *> support_layers;
                        support_layers.reserve(path->polyline.points.size());
                        bool lies_on_planar_courses = true;
                        for (const Point3 &point : path->polyline.points) {
                            const double absolute_z = layer->print_z + unscale_(point.z());
                            const auto nearest = std::min_element(
                                object->layers().begin(), object->layers().end(),
                                [absolute_z](const Layer *left, const Layer *right) {
                                    return std::abs(left->print_z - absolute_z) <
                                           std::abs(right->print_z - absolute_z);
                                });
                            if (nearest == object->layers().end() ||
                                std::abs((*nearest)->print_z - absolute_z) > 1e-4) {
                                lies_on_planar_courses = false;
                                break;
                            }
                            support_layers.emplace_back(*nearest);
                        }
                        if (!lies_on_planar_courses)
                            continue;
                        ++planar_anchored_paths;
                        for (size_t point_idx = 0;
                             point_idx < path->polyline.points.size(); ++point_idx)
                            unsupported_anchor_points += !Geometry::contains(
                                support_layers[point_idx]->lslices,
                                path->polyline.points[point_idx].to_point());
                    }
                };
                inspect(region->perimeters);
                inspect(region->fills);

                const auto inspect_wall_anchors = [&](auto &&self,
                                                      const ExtrusionEntity &entity) -> void {
                    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity)) {
                        if (!path->nonplanar_surface ||
                            !is_perimeter(path->role()) || path->polyline.points.empty())
                            return;
                        nonplanar_wall_owner_layers.insert(layer);
                        std::vector<const Layer *> support_layers;
                        support_layers.reserve(path->polyline.points.size());
                        bool lies_on_planar_courses = true;
                        for (const Point3 &point : path->polyline.points) {
                                const double absolute_z =
                                    layer->print_z + unscale_(point.z());
                                const auto support = std::min_element(
                                    object->layers().begin(), object->layers().end(),
                                    [absolute_z](const Layer *left, const Layer *right) {
                                        return std::abs(left->print_z - absolute_z) <
                                               std::abs(right->print_z - absolute_z);
                                    });
                                if (support == object->layers().end() ||
                                    std::abs((*support)->print_z - absolute_z) > 1e-4) {
                                    lies_on_planar_courses = false;
                                    break;
                                }
                                support_layers.emplace_back(*support);
                        }
                        if (!lies_on_planar_courses)
                            return;
                        ++planar_wall_anchor_paths;
                        for (size_t point_idx = 0;
                             point_idx < path->polyline.points.size(); ++point_idx) {
                            const ExPolygons support_area = offset_ex(
                                support_layers[point_idx]->lslices,
                                scale_(0.5 * double(path->width) + 0.02));
                            unsupported_wall_anchor_points += !Geometry::contains(
                                support_area,
                                path->polyline.points[point_idx].to_point());
                        }
                    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
                        for (const ExtrusionPath &path : loop->paths)
                            self(self, path);
                    } else if (const auto *multipath =
                                   dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
                        for (const ExtrusionPath &path : multipath->paths)
                            self(self, path);
                    } else if (const auto *collection =
                                   dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
                        for (const ExtrusionEntity *child : collection->entities)
                            self(self, *child);
                    }
                };
                inspect_wall_anchors(inspect_wall_anchors, region->perimeters);
                inspect_wall_anchors(inspect_wall_anchors, region->fills);
            }
            owning_layers += owns_nonplanar_top;
        }
    }

    // The configured count is the minimum transition depth. A broad slope may
    // require additional intermediate courses so its high edge still remains
    // within one nominal layer of the preceding extrusion.
    REQUIRE(owning_layers >= 5);
    // The boundary-reaching wedge retains a complete multi-course perimeter
    // transition, and that transition begins no later than the generated top
    // stack. Perimeters keep their original slicer owner when one connected
    // wall spans multiple absolute Z courses, so owner IDs need not match the
    // raster owners one-for-one.
    REQUIRE_FALSE(nonplanar_top_owner_layers.empty());
    REQUIRE(nonplanar_wall_owner_layers.size() >= 5);
    const auto earlier_layer = [](const Layer *left, const Layer *right) {
        return left->print_z < right->print_z;
    };
    const Layer *first_top_owner = *std::min_element(
        nonplanar_top_owner_layers.begin(), nonplanar_top_owner_layers.end(),
        earlier_layer);
    const Layer *first_wall_owner = *std::min_element(
        nonplanar_wall_owner_layers.begin(), nonplanar_wall_owner_layers.end(),
        earlier_layer);
    CAPTURE(first_top_owner->id(), first_wall_owner->id());
    CHECK(first_wall_owner->print_z <= first_top_owner->print_z + 1e-4);
    // Generated skin courses remain anchored to their actual XY footprint.
    // Retained walls are checked separately below because their original XY
    // contour is intentionally preserved for dimensional accuracy.
    // The ordinary planar foundation is no longer emitted a second time as a
    // generated blend=0 path. Exact-course generated anchors must therefore
    // be absent; support for the first generated course is checked against
    // the preceding retained layer below.
    CHECK(planar_anchored_paths == 0);
    CHECK(unsupported_anchor_points == 0);
    // The retained perimeter graph starts on the same real flat course as the
    // dense foundation.  This is intentional: it gives the first rising wall
    // an attached planar endpoint instead of beginning on the already stepped
    // contour.  Every such anchor still has to lie inside the supporting bead
    // envelope.
    CHECK(planar_wall_anchor_paths > 0);
    CHECK(unsupported_wall_anchor_points == 0);
    // Corresponding points in the blended stack may be owned by different
    // nominal layers, but no adjacent deposited course may be farther apart
    // than the configured layer height. Selecting the support course below
    // the requested transition depth adds a rounding remainder to the stack
    // and makes the visible finishing shell float above its predecessor.
    size_t aligned_stack_points = 0;
    size_t foundation_points = 0;
    size_t unsupported_foundation_points = 0;
    double maximum_shell_gap = 0.;
    double maximum_foundation_gap = 0.;
    double minimum_inferred_anchor = std::numeric_limits<double>::max();
    double maximum_inferred_anchor = std::numeric_limits<double>::lowest();
    REQUIRE(print.objects().size() == 1);
    const auto object_layers = print.objects().front()->layers();
    for (auto &[xy, z_values] : shell_z_by_xy) {
        std::sort(z_values.begin(), z_values.end());
        z_values.erase(std::unique(z_values.begin(), z_values.end(),
            [](double left, double right) {
                return std::abs(left - right) < 1e-4;
            }), z_values.end());
        if (z_values.size() < 5)
            continue;
        ++aligned_stack_points;
        const double shell_step =
            (z_values.back() - z_values.front()) /
            double(z_values.size() - 1);
        const double inferred_anchor = z_values.front() - shell_step;
        minimum_inferred_anchor = std::min(
            minimum_inferred_anchor, inferred_anchor);
        maximum_inferred_anchor = std::max(
            maximum_inferred_anchor, inferred_anchor);
        const double first_generated_z = z_values.front();
        auto foundation = std::lower_bound(
            object_layers.begin(), object_layers.end(),
            first_generated_z - 1e-4,
            [](const Layer *layer, double z) { return layer->print_z < z; });
        if (foundation != object_layers.begin()) {
            --foundation;
            ++foundation_points;
            const double gap = first_generated_z - (*foundation)->print_z;
            maximum_foundation_gap = std::max(maximum_foundation_gap, gap);
            const ExPolygons support_area = offset_ex(
                (*foundation)->lslices, scale_(0.25));
            unsupported_foundation_points +=
                gap > 0.2 + 1e-4 ||
                !Geometry::contains(support_area, Point(xy.first, xy.second));
        }
        for (size_t idx = 1; idx < z_values.size(); ++idx)
            maximum_shell_gap = std::max(
                maximum_shell_gap, z_values[idx] - z_values[idx - 1]);
    }
    CHECK(aligned_stack_points > 0);
    CHECK(foundation_points > 0);
    CHECK(unsupported_foundation_points == 0);
    CHECK(maximum_foundation_gap <= 0.2 + 1e-4);
    CHECK(maximum_shell_gap <= 0.2 + 1e-4);
    // Every XY column must extrapolate to the same horizontal foundation.
    // Anchoring each point to surface_z - shell_depth reproduces a hidden
    // staircase and yields different inferred anchors across the roof.
    CHECK(maximum_inferred_anchor - minimum_inferred_anchor <= 0.02);
    // Walls remain in the original perimeter graph rather than being cloned
    // from a finishing-patch outline. A vertical-sided wedge legitimately has
    // the same XY wall footprint on every course, so connectivity—not forced
    // lateral motion—is the invariant here.
    CHECK_FALSE(nonplanar_perimeter_paths(print).empty());
    // At least one owner contains a path originating more than one nominal
    // layer below it, preserving continuous multi-layer Z motion.
    CHECK(deepest_relative_z < -scale_(0.2));

    // An open perimeter may be surface-following only when a loop or
    // multipath owns its connection to the neighboring planar wall pieces.
    // A bare path child would be scheduled as a separate operation and leave
    // the characteristic retract/travel-isolated raised Arachne fragment.
    size_t isolated_open_nonplanar_walls = 0;
    const auto inspect_collection = [&](auto &&self,
                                        const ExtrusionEntityCollection &collection) -> void {
        for (const ExtrusionEntity *entity : collection.entities) {
            if (const auto *path = dynamic_cast<const ExtrusionPath *>(entity)) {
                isolated_open_nonplanar_walls +=
                    path->nonplanar_surface && is_perimeter(path->role()) &&
                    path->first_point() != path->last_point();
            } else if (const auto *nested =
                           dynamic_cast<const ExtrusionEntityCollection *>(entity)) {
                self(self, *nested);
            }
            // ExtrusionLoop and ExtrusionMultiPath intentionally own open
            // component paths, so their children are connected by contract.
        }
    };
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions()) {
                inspect_collection(inspect_collection, region->perimeters);
                inspect_collection(inspect_collection, region->fills);
            }
    CHECK(isolated_open_nonplanar_walls == 0);
}

TEST_CASE("Non-planar walls respect the configured wall sequence", "[Fill][NonplanarSurface][WallSequence]")
{
    const auto generate_order = [](const char *wall_sequence) {
        Print print;
        Slic3r::Test::init_and_process_print(
            {shallow_top_wedge(8.)}, print,
            {{"nonplanar_top_surface", 1},
             {"nonplanar_top_surface_max_angle", 45},
             {"wall_loops", 3},
             {"wall_sequence", wall_sequence},
             {"zaa_enabled", 0},
             {"layer_height", 0.2}});
        return nonplanar_wall_inset_order(print);
    };

    CHECK(generate_order("inner wall/outer wall") ==
          std::vector<unsigned int>{2, 1, 0});
    CHECK(generate_order("outer wall/inner wall") ==
          std::vector<unsigned int>{0, 1, 2});
    CHECK(generate_order("inner-outer-inner wall") ==
          std::vector<unsigned int>{2, 0, 1});
}

TEST_CASE("True non-planar top surfaces leave paths planar when disabled", "[Fill][NonplanarSurface]")
{
    Print print;
    Slic3r::Test::init_and_process_print(
        {shallow_top_wedge(8.)}, print,
        {{"nonplanar_top_surface", 0}, {"zaa_enabled", 0}, {"layer_height", 0.2}});

    CHECK(nonplanar_top_paths(print).empty());
}

TEST_CASE("Horizontal top surfaces are not classified as non-planar", "[Fill][NonplanarSurface]")
{
    Print print;
    Slic3r::Test::init_and_process_print(
        {Slic3r::Test::cube(10.)}, print,
        {{"nonplanar_top_surface", 1},
         {"nonplanar_top_surface_max_angle", 45},
         {"zaa_enabled", 0},
         {"layer_height", 0.2}});

    CHECK(nonplanar_top_paths(print).empty());
    CHECK(nonplanar_perimeter_paths(print).empty());
}

TEST_CASE("True non-planar top surfaces respect the configured surface angle", "[Fill][NonplanarSurface]")
{
    Print accepted;
    Slic3r::Test::init_and_process_print(
        {shallow_top_wedge(8.)}, accepted,
        {{"nonplanar_top_surface", 1},
         {"nonplanar_top_surface_max_angle", 10},
         {"zaa_enabled", 0},
         {"layer_height", 0.2}});
    CHECK_FALSE(nonplanar_top_paths(accepted).empty());

    Print rejected;
    Slic3r::Test::init_and_process_print(
        {shallow_top_wedge(8.)}, rejected,
        {{"nonplanar_top_surface", 1},
         {"nonplanar_top_surface_max_angle", 5},
         {"zaa_enabled", 0},
         {"layer_height", 0.2}});
    CHECK(nonplanar_top_paths(rejected).empty());
}

TEST_CASE("Doubly curved top surfaces use closed non-planar tracks", "[Fill][NonplanarSurface]")
{
    TriangleMesh dome = make_sphere(10., PI / 24.);
    dome.translate(Vec3f(0.f, 0.f, 10.f));

    Print print;
    Slic3r::Test::init_and_process_print(
        {dome}, print,
        {{"nonplanar_top_surface", 1},
         {"nonplanar_top_surface_max_angle", 20},
         {"nonplanar_top_surface_resolution", 0.2},
         {"nonplanar_top_surface_min_height", 0.05},
         {"top_surface_pattern", "rectilinear"},
         {"zaa_enabled", 0},
         {"layer_height", 0.2}});

    const std::vector<ExtrusionPath> paths = nonplanar_top_paths(print);
    REQUIRE_FALSE(paths.empty());
    // The configured rectilinear pattern remains appropriate for a plane,
    // even one sloped in both X and Y. A sphere's gradient changes across
    // both XY directions, so the adaptive non-planar generator must instead
    // follow it with closed concentric courses.
    CHECK(std::any_of(paths.begin(), paths.end(), [](const ExtrusionPath &path) {
        return path.polyline.points.size() > 3 &&
               path.polyline.points.front() == path.polyline.points.back();
    }));
}

TEST_CASE("Compact wall-only slopes without a solid source course stay planar",
          "[Fill][NonplanarSurface][PerimeterOnly]")
{
    Polygon footprint = make_circle_num_segments(scale_(1.7), 48);
    const double slope = std::tan(8. * M_PI / 180.);
    TriangleMesh model = extrude_sloped_footprint(
        ExPolygon(footprint), [](double) { return 0.; },
        [slope](double x) { return 8. + slope * x; });

    Print print;
    Slic3r::Test::init_and_process_print(
        {std::move(model)}, print,
        {{"nonplanar_top_surface", 1},
         {"nonplanar_top_surface_max_angle", 45},
         {"wall_loops", 3}, {"zaa_enabled", 0}, {"layer_height", 0.2}});

    // A wall-only component has no dense course from which to build a bonded
    // transition. Warping its rings alone creates an isolated recess (the
    // Benchy flag-hole regression), so preserve the conventional slice until
    // a coherent skin and support course exist.
    CHECK(nonplanar_top_paths(print).empty());
    const std::vector<ExtrusionPath> walls = nonplanar_perimeter_paths(print);
    CHECK(walls.empty());
}

TEST_CASE("Narrow supported top caps remain eligible for non-planar finishing",
          "[Fill][NonplanarSurface][Regression]")
{
    // This has the proportions of Benchy's narrow sloped bow cap: its area is
    // below the broad-patch threshold and its short dimension is below 10 mm,
    // but it still contains a real solid source course. Size alone must not
    // force the hybrid mode to Z-contour geometry that passed support and
    // toolhead-clearance validation.
    Polygon footprint({Point::new_scale(0., 0.), Point::new_scale(6., 0.),
                       Point::new_scale(6., 14.), Point::new_scale(0., 14.)});
    const double slope = std::tan(8. * M_PI / 180.);
    TriangleMesh model = extrude_sloped_footprint(
        ExPolygon(footprint), [](double) { return 0.; },
        [slope](double x) { return 8. + slope * x; });

    Print print;
    Slic3r::Test::init_and_process_print(
        {std::move(model)}, print,
        {{"top_surface_z_mode", "nonplanar_with_z_contouring_fallback"},
         {"nonplanar_top_surface", true}, {"zaa_enabled", true},
         {"nonplanar_top_surface_max_angle", 45},
         {"wall_loops", 3}, {"layer_height", 0.2}});

    CHECK_FALSE(nonplanar_top_paths(print).empty());
    CHECK_FALSE(nonplanar_perimeter_paths(print).empty());
    CHECK_FALSE(has_planar_path_inside_nonplanar_top_coverage(print));
}

TEST_CASE("Non-planar skins exclude normal paths from their top footprint", "[Fill][NonplanarSurface]")
{
    Print nonplanar;
    Slic3r::Test::init_and_process_print(
        {shallow_top_wedge(8.)}, nonplanar,
        {{"nonplanar_top_surface", 1}, {"wall_loops", 3},
         {"staggered_perimeters", 1},
         {"top_surface_expansion", 0.}, {"zaa_enabled", 0}, {"layer_height", 0.2}});

    // Clipping closed loops may split one wall into several path entities, so
    // raw entity counts do not measure replacement. Verify the actual overlap
    // below instead.
    CHECK(perimeter_path_count(nonplanar, erExternalPerimeter) > 0);
    // Brick courses remain active in buried wall regions, while wall sections
    // accepted into the non-planar replacement use only the mesh-following Z
    // schedule. This permits both features in one print without combining two
    // independent offsets on the same physical extrusion.
    CHECK(staggered_perimeter_path_count(nonplanar) > 0);
    const std::vector<ExtrusionPath> nonplanar_walls =
        nonplanar_perimeter_paths(nonplanar);
    CHECK_FALSE(std::any_of(nonplanar_walls.begin(), nonplanar_walls.end(),
        [](const ExtrusionPath &path) { return path.staggered_perimeter; }));
    CHECK(std::any_of(nonplanar_walls.begin(), nonplanar_walls.end(),
        [](const ExtrusionPath &path) { return path.inset_idx == 0; }));
    CHECK(std::any_of(nonplanar_walls.begin(), nonplanar_walls.end(),
        [](const ExtrusionPath &path) { return path.inset_idx == 1; }));
    CHECK(std::any_of(nonplanar_walls.begin(), nonplanar_walls.end(),
        [](const ExtrusionPath &path) { return path.inset_idx == 2; }));
    CHECK_FALSE(has_planar_path_inside_nonplanar_top_coverage(nonplanar));
}

TEST_CASE("Non-planar top skins smooth exposed interlocking wall courses",
          "[Fill][NonplanarSurface][InterlockingWalls]")
{
    TriangleMesh interlocking_control = Slic3r::Test::cube(8.);
    interlocking_control.translate(Vec3f(20.f, 0.f, 0.f));
    Print print;
    Slic3r::Test::init_and_process_print(
        {shallow_top_wedge(8.), std::move(interlocking_control)}, print,
        {{"nonplanar_top_surface", 1},
         {"perimeter_layering", "interlocking_walls"},
         {"interlocking_wall_amplitude", 0.08},
         {"interlocking_wall_wavelength", 6.0},
         {"interlocking_wall_resolution", 0.4},
         {"wall_loops", 3}, {"zaa_enabled", 0}, {"layer_height", 0.2}});

    // Interlocking courses remain available on ordinary geometry in the same
    // print. A wall section whose entire staircase is replaced by the angled
    // transition has one authoritative Z schedule and is intentionally not
    // interlocked a second time.
    CHECK(interlocking_perimeter_path_count(print) > 0);
    const std::vector<ExtrusionPath> exposed = nonplanar_perimeter_paths(print);
    REQUIRE_FALSE(exposed.empty());
    // Once accepted into the exposed skin, the mesh-following schedule is the
    // only Z schedule and preview role on that physical wall section.
    CHECK_FALSE(std::any_of(exposed.begin(), exposed.end(),
        [](const ExtrusionPath &path) { return path.nonplanar_interlocking_wall; }));
    // Covered by the emitted wall/top intersection audit below; keeping the
    // older all-layer polygon scan here makes this focused regression
    // needlessly quadratic for a complete coplanar roof.
}

TEST_CASE("Non-planar Benchy roof output preserves chimney walls and connected blends",
          "[Fill][NonplanarSurface][GCode][Regression]")
{
    const std::string perimeter_mode = GENERATE(
        std::string("standard"), std::string("brick"),
        std::string("interlocking_walls"));
    CAPTURE(perimeter_mode);

    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({
        {"top_surface_z_mode", "nonplanar_top_surface"},
        {"nonplanar_top_surface", true},
        {"nonplanar_top_surface_max_angle", 45.},
        {"nonplanar_top_surface_layers", 5},
        {"nonplanar_top_surface_resolution", 0.3},
        {"nonplanar_top_surface_min_height", 0.05},
        {"perimeter_layering", perimeter_mode},
        {"staggered_perimeter_offset", 1},
        {"interlocking_wall_amplitude", 0.08},
        {"interlocking_wall_wavelength", 6.0},
        {"interlocking_wall_resolution", 0.4},
        {"wall_generator", "arachne"},
        {"wall_loops", 3},
        {"top_shell_layers", 4},
        {"sparse_infill_density", "15%"},
        {"layer_height", 0.25},
        {"initial_layer_print_height", 0.25},
        {"skirt_loops", 0},
        {"brim_type", "no_brim"},
        {"zaa_enabled", false}
    });

    Print print;
    Model model;
    std::vector<TriangleMesh> meshes;
    meshes.emplace_back(benchy_roof_with_chimney(8.));
    Slic3r::Test::init_print(std::move(meshes), print, model, config, nullptr,
                             false);
    REQUIRE_NOTHROW(print.process());

    const std::vector<ExtrusionPath> top_paths = nonplanar_top_paths(print);
    const std::vector<ExtrusionPath> wall_paths = nonplanar_perimeter_paths(print);
    REQUIRE_FALSE(top_paths.empty());
    REQUIRE_FALSE(wall_paths.empty());
    // The finishing skin owns the requested wall stack. Feature-transition
    // rings around the chimney are a separate dependency and must not hide a
    // missing roof inset in this assertion: omitting an inner roof wall while
    // still insetting the raster for it leaves an unextruded eave band.
    std::set<unsigned int> roof_wall_insets;
    for (const ExtrusionPath &path : wall_paths)
        if (!path.nonplanar_feature_transition &&
            !path.nonplanar_transition && path.inset_idx <= 2)
            roof_wall_insets.insert(path.inset_idx);
    CAPTURE(roof_wall_insets);
    CHECK(roof_wall_insets.count(0) == 1);
    CHECK(roof_wall_insets.count(1) == 1);
    CHECK(roof_wall_insets.count(2) == 1);
    // A nominal inner wall crossing the accepted roof patch is the spurious
    // rectangular "box" seen in the Benchy preview. The top patch must replace
    // that section instead of adding a second closed wall around itself.
    CHECK_FALSE(has_planar_path_inside_nonplanar_top_coverage(print));

    const std::string generated_gcode = Slic3r::Test::gcode(print);
    std::ofstream("/tmp/nonplanar-benchy-focused.gcode") << generated_gcode;
    REQUIRE_FALSE(generated_gcode.empty());
    CHECK(generated_gcode.find(";TYPE:Non-planar top surface") != std::string::npos);
    CHECK(generated_gcode.find(";TYPE:Non-planar transition") != std::string::npos);
    // Feature deferral changes the wall schedule only. Source-layer infill is
    // the physical foundation for both the draped roof and those walls; moving
    // it down into a later dependency course creates unsupported paths and may
    // cross a bridge that has already been emitted.
    CHECK(generated_gcode.find("surface-following feature fill") ==
          std::string::npos);

    // Validate support in actual print order. Geometry-only containment can
    // pass while the commit step clips away the wall course immediately below
    // the first draped perimeter. That produced a 0.4--0.6 mm vertical gap on
    // the real Benchy even though every generated point remained inside the
    // model slice. A non-planar wall bead must overlap material deposited
    // earlier and no more than one nominal course below it.
    struct DepositedSegment {
        Vec3d start;
        Vec3d end;
    };
    std::vector<DepositedSegment> deposited_segments;
    size_t unsupported_nonplanar_segments = 0;
    std::map<std::string, size_t> unsupported_nonplanar_by_comment;
    std::vector<std::string> unsupported_nonplanar_samples;
    std::string support_type;
    GCodeReader support_reader;
    support_reader.parse_buffer(generated_gcode,
        [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
            const std::string comment(line.comment());
            if (comment.rfind("TYPE:", 0) == 0)
                support_type = comment.substr(5);
            if (!line.extruding(self) || line.dist_XY(self) <= 0.01)
                return;
            const Vec3d start(self.x(), self.y(), self.z());
            const Vec3d end(line.new_X(self), line.new_Y(self), line.new_Z(self));
            if (support_type == "Non-planar transition" ||
                support_type == "Non-planar top surface" ||
                comment == "non-planar perimeter" ||
                comment == "surface-following feature perimeter" ||
                comment == "non-planar top surface") {
                // A deliberately classified overhang perimeter is printable
                // by cantilever/bridge policy and need not have a bead directly
                // below every midpoint. The regression targets ordinary
                // surface courses that incorrectly float while claiming to be
                // supported; bridge paths are validated by bridge tests.
                if (comment.find("(bridge)") != std::string::npos) {
                    deposited_segments.push_back({
                        Vec3d(self.x(), self.y(), self.z()),
                        Vec3d(line.new_X(self), line.new_Y(self), line.new_Z(self))});
                    return;
                }
                const Vec3d midpoint = 0.5 * (start + end);
                bool supported = false;
                double nearest_support_xy = std::numeric_limits<double>::max();
                double nearest_support_gap = std::numeric_limits<double>::max();
                for (auto candidate = deposited_segments.rbegin();
                     candidate != deposited_segments.rend(); ++candidate) {
                    const Vec2d delta =
                        candidate->end.head<2>() - candidate->start.head<2>();
                    const double denominator = delta.squaredNorm();
                    const double t = denominator <= EPSILON ? 0. : std::clamp(
                        (midpoint.head<2>() - candidate->start.head<2>()).dot(delta) /
                            denominator,
                        0., 1.);
                    const Vec2d projected = candidate->start.head<2>() + t * delta;
                    const double xy_distance =
                        (projected - midpoint.head<2>()).norm();
                    const double support_z = candidate->start.z() +
                        t * (candidate->end.z() - candidate->start.z());
                    const double gap = midpoint.z() - support_z;
                    if (xy_distance < nearest_support_xy) {
                        nearest_support_xy = xy_distance;
                        nearest_support_gap = gap;
                    }
                    if (xy_distance > 0.4)
                        continue;
                    if (gap >= 0.04 && gap <=
                            config.opt_float("layer_height") + 0.06) {
                        supported = true;
                        break;
                    }
                }
                if (!supported) {
                    ++unsupported_nonplanar_segments;
                    ++unsupported_nonplanar_by_comment[
                        comment.empty() ? support_type : comment];
                    if (unsupported_nonplanar_samples.size() < 8)
                        unsupported_nonplanar_samples.emplace_back(
                            "(" + std::to_string(midpoint.x()) + "," +
                            std::to_string(midpoint.y()) + "," +
                            std::to_string(midpoint.z()) + ") " +
                            (comment.empty() ? support_type : comment) +
                            " nearest_xy=" + std::to_string(nearest_support_xy) +
                            " gap=" + std::to_string(nearest_support_gap));
                }
            }
            deposited_segments.push_back({start, end});
        });
    std::string unsupported_summary;
    for (const auto &[comment, count] : unsupported_nonplanar_by_comment)
        unsupported_summary += comment + '=' + std::to_string(count) + ' ';
    CAPTURE(unsupported_summary);
    CAPTURE(unsupported_nonplanar_samples);
    CHECK(unsupported_nonplanar_segments == 0);

    struct EmittedSegment {
        Vec3d start;
        Vec3d end;
        double middle_z;
        size_t chain;
        std::string type;
        std::string comment;
    };
    std::vector<EmittedSegment> emitted_wall_segments;
    std::vector<EmittedSegment> emitted_top_segments;
    std::vector<EmittedSegment> emitted_conventional_segments;
    size_t transition_segments = 0;
    size_t transition_z_moves = 0;
    std::map<std::string, std::pair<size_t, size_t>> transition_moves_by_layer;
    struct EmittedFeatureChain {
        Vec3d start;
        Vec3d end;
        size_t segments {0};
    };
    std::map<size_t, EmittedFeatureChain> emitted_feature_chains;
    double maximum_transition_gradient = 0.;
    std::string emitted_type;
    std::string emitted_layer;
    size_t emitted_chain = 0;
    std::optional<Vec3d> previous_extrusion_end;
    GCodeReader geometry_reader;
    geometry_reader.parse_buffer(generated_gcode,
        [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
            const std::string comment(line.comment());
            if (comment.rfind("TYPE:", 0) == 0)
                emitted_type = comment.substr(5);
            else if (comment.rfind("Z:", 0) == 0)
                emitted_layer = comment.substr(2);
            if (!line.extruding(self) || line.dist_XY(self) <= 0.01)
                return;
            const Vec3d start(self.x(), self.y(), self.z());
            const Vec3d end(line.new_X(self), line.new_Y(self), line.new_Z(self));
            if (!previous_extrusion_end ||
                (start - *previous_extrusion_end).norm() > 0.01)
                ++emitted_chain;
            previous_extrusion_end = end;
            if (comment.find("surface-following feature perimeter") !=
                    std::string::npos) {
                EmittedFeatureChain &chain =
                    emitted_feature_chains[emitted_chain];
                if (chain.segments == 0)
                    chain.start = start;
                chain.end = end;
                ++chain.segments;
            }
            if (emitted_type == "Non-planar transition") {
                const double dz = std::abs(end.z() - start.z());
                ++transition_segments;
                transition_z_moves += dz > 0.001;
                auto &layer_moves = transition_moves_by_layer[emitted_layer];
                ++layer_moves.first;
                layer_moves.second += dz > 0.001;
                maximum_transition_gradient = std::max(
                    maximum_transition_gradient, dz / line.dist_XY(self));
            }
            // GCodeReader preserves optional whitespace/tool annotations in
            // inline comments. Classify by the semantic tag so a top-surface
            // move emitted under the shared transition TYPE is not also
            // counted as a wall by the collision audit below.
            const bool inline_nonplanar_top =
                comment.find("non-planar top surface") != std::string::npos;
            const bool inline_nonplanar_wall =
                comment.find("non-planar perimeter") != std::string::npos ||
                comment.find("surface-following feature perimeter") !=
                    std::string::npos;
            const bool wall = inline_nonplanar_wall ||
                              (!inline_nonplanar_top &&
                              (emitted_type == "Outer wall" ||
                              emitted_type == "Inner wall" ||
                              emitted_type == "Overhang wall" ||
                              emitted_type == "Brick wall" ||
                              emitted_type == "Interlocking inner wall" ||
                              emitted_type == "Non-planar transition"));
            const double middle_z = 0.5 * (start.z() + end.z());
            if (wall)
                emitted_wall_segments.push_back(
                    {start, end, middle_z, emitted_chain, emitted_type, comment});
            if (inline_nonplanar_top ||
                emitted_type == "Non-planar top surface")
                emitted_top_segments.push_back(
                    {start, end, middle_z, emitted_chain, emitted_type, comment});
            if (!inline_nonplanar_top && !inline_nonplanar_wall &&
                emitted_type != "Non-planar transition" &&
                emitted_type != "Non-planar top surface")
                emitted_conventional_segments.push_back(
                    {start, end, middle_z, emitted_chain, emitted_type, comment});
        });
    REQUIRE(transition_segments > 0);
    CHECK(transition_z_moves > 0);
    // A horizontal foundation is an ordinary planar wall/top course.  Do not
    // relabel a newly generated fixed-Z closed loop as a transition: that was
    // the first tell-tale sign that the non-planar stack had been appended as
    // a post-process shell instead of deforming the slicer's wall graph. Every
    // layer advertised as a transition must contain an extrusion that truly
    // changes Z; otherwise the original planar course must remain authoritative.
    for (const auto &[layer, moves] : transition_moves_by_layer) {
        CAPTURE(layer, moves.first, moves.second);
        CHECK(moves.second > 0);
    }
    REQUIRE_FALSE(emitted_feature_chains.empty());
    for (const auto &[chain_index, chain] : emitted_feature_chains) {
        CAPTURE(chain_index, chain.segments, chain.start, chain.end);
        CHECK((chain.end - chain.start).norm() <= 0.05);
    }
    // A transition may traverse a long sloped roof in one move, so validate
    // gradient rather than absolute dZ. It must remain below the configured
    // 45-degree surface limit.
    CHECK(maximum_transition_gradient < 1.);

    // Select a stable band below the highest emitted wall rather than assuming
    // an absolute G-code Z origin. Only the hollow chimney exists there.
    REQUIRE_FALSE(emitted_wall_segments.empty());
    const double highest_wall_z = std::max_element(
        emitted_wall_segments.begin(), emitted_wall_segments.end(),
        [](const EmittedSegment &lhs, const EmittedSegment &rhs) {
            return lhs.middle_z < rhs.middle_z;
        })->middle_z;
    std::vector<EmittedSegment> chimney_wall_segments;
    std::copy_if(emitted_wall_segments.begin(), emitted_wall_segments.end(),
        std::back_inserter(chimney_wall_segments),
        [highest_wall_z](const EmittedSegment &segment) {
            return segment.middle_z > highest_wall_z - 5. &&
                   segment.middle_z < highest_wall_z - 2.5;
        });
    CAPTURE(highest_wall_z);
    REQUIRE_FALSE(chimney_wall_segments.empty());
    Vec2d minimum_xy = chimney_wall_segments.front().start.head<2>();
    Vec2d maximum_xy = minimum_xy;
    for (const EmittedSegment &segment : chimney_wall_segments) {
        minimum_xy = minimum_xy.cwiseMin(segment.start.head<2>());
        minimum_xy = minimum_xy.cwiseMin(segment.end.head<2>());
        maximum_xy = maximum_xy.cwiseMax(segment.start.head<2>());
        maximum_xy = maximum_xy.cwiseMax(segment.end.head<2>());
    }
    const Vec2d chimney_center = 0.5 * (minimum_xy + maximum_xy);
    double bore_wall_length = 0.;
    double outside_wall_length = 0.;
    for (const EmittedSegment &segment : chimney_wall_segments) {
        const Vec2d midpoint =
            0.5 * (segment.start.head<2>() + segment.end.head<2>());
        const double length =
            (segment.end.head<2>() - segment.start.head<2>()).norm();
        if ((midpoint - chimney_center).norm() < 3.2)
            bore_wall_length += length;
        else
            outside_wall_length += length;
    }
    CAPTURE(minimum_xy, maximum_xy, bore_wall_length, outside_wall_length);
    CHECK((maximum_xy.x() - minimum_xy.x()) > 8.);
    CHECK((maximum_xy.y() - minimum_xy.y()) > 8.);
    CHECK(bore_wall_length > 10.);
    CHECK(outside_wall_length > 30.);

    // A top raster may end on the innermost perimeter, but it must never pass
    // through a retained wall and continue on the other side at the same Z.
    // The real Benchy failure presented exactly this way: long roof rows
    // crossed a rectangular retained inner wall, creating duplicate material
    // and a nozzle collision in the preview.
    size_t wall_top_crossings = 0;
    std::vector<std::string> wall_top_crossing_samples;
    std::map<int, std::vector<size_t>> top_segments_by_x;
    for (size_t index = 0; index < emitted_top_segments.size(); ++index) {
        const EmittedSegment &top = emitted_top_segments[index];
        const double minimum_x = std::min(top.start.x(), top.end.x());
        const double maximum_x = std::max(top.start.x(), top.end.x());
        for (int bucket = int(std::floor(minimum_x / 2.));
             bucket <= int(std::floor(maximum_x / 2.)); ++bucket)
            top_segments_by_x[bucket].push_back(index);
    }
    for (const EmittedSegment &wall : emitted_wall_segments) {
        const Vec2d wall_start = wall.start.head<2>();
        const Vec2d wall_delta = wall.end.head<2>() - wall_start;
        std::set<size_t> candidate_indices;
        const double wall_minimum_x = std::min(wall.start.x(), wall.end.x());
        const double wall_maximum_x = std::max(wall.start.x(), wall.end.x());
        for (int bucket = int(std::floor(wall_minimum_x / 2.));
             bucket <= int(std::floor(wall_maximum_x / 2.)); ++bucket) {
            const auto found = top_segments_by_x.find(bucket);
            if (found != top_segments_by_x.end())
                candidate_indices.insert(found->second.begin(),
                                         found->second.end());
        }
        for (const size_t top_index : candidate_indices) {
            const EmittedSegment &top = emitted_top_segments[top_index];
            const Vec2d top_start = top.start.head<2>();
            const Vec2d top_delta = top.end.head<2>() - top_start;
            const Vec2d wall_min =
                wall.start.head<2>().cwiseMin(wall.end.head<2>());
            const Vec2d wall_max =
                wall.start.head<2>().cwiseMax(wall.end.head<2>());
            const Vec2d top_min =
                top.start.head<2>().cwiseMin(top.end.head<2>());
            const Vec2d top_max =
                top.start.head<2>().cwiseMax(top.end.head<2>());
            if (wall_max.x() < top_min.x() || top_max.x() < wall_min.x() ||
                wall_max.y() < top_min.y() || top_max.y() < wall_min.y() ||
                std::max(wall.start.z(), wall.end.z()) <
                    std::min(top.start.z(), top.end.z()) - 0.05 ||
                std::max(top.start.z(), top.end.z()) <
                    std::min(wall.start.z(), wall.end.z()) - 0.05)
                continue;
            const double denominator =
                wall_delta.x() * top_delta.y() -
                wall_delta.y() * top_delta.x();
            if (std::abs(denominator) <= EPSILON)
                continue;
            const Vec2d offset = top_start - wall_start;
            const double wall_t =
                (offset.x() * top_delta.y() -
                 offset.y() * top_delta.x()) / denominator;
            const double top_t =
                (offset.x() * wall_delta.y() -
                 offset.y() * wall_delta.x()) / denominator;
            // Endpoint contact is the intended wall-to-skin bond. Only a
            // proper interior crossing is an overlap.
            if (wall_t <= 0.01 || wall_t >= 0.99 ||
                top_t <= 0.01 || top_t >= 0.99)
                continue;
            const double wall_z = wall.start.z() +
                wall_t * (wall.end.z() - wall.start.z());
            const double top_z = top.start.z() +
                top_t * (top.end.z() - top.start.z());
            if (std::abs(wall_z - top_z) < 0.05) {
                ++wall_top_crossings;
                if (wall_top_crossing_samples.size() < 8)
                    wall_top_crossing_samples.emplace_back(
                        "wall=(" + std::to_string(wall.start.x()) + "," +
                        std::to_string(wall.start.y()) + " -> " +
                        std::to_string(wall.end.x()) + "," +
                        std::to_string(wall.end.y()) + ",z=" +
                        std::to_string(wall_z) + ",chain=" +
                        std::to_string(wall.chain) + "," + wall.type + "," +
                        wall.comment +
                        ") top=(" + std::to_string(top.start.x()) + "," +
                        std::to_string(top.start.y()) + " -> " +
                        std::to_string(top.end.x()) + "," +
                        std::to_string(top.end.y()) + ",z=" +
                        std::to_string(top_z) + ")");
            }
        }
    }
    CAPTURE(wall_top_crossing_samples);
    CHECK(wall_top_crossings == 0);

    // Successive roof/cabin wall courses may have slightly different XY
    // outlines. They may cross geometrically only when their physical Z is
    // separated; emitting both centerlines at the same Z deposits duplicate
    // material and makes the nozzle traverse the preceding course. This
    // catches intersections inside sparse Arachne segments, not just at the
    // source vertices used by the transition constraint.
    size_t wall_wall_same_z_crossings = 0;
    std::vector<std::string> wall_wall_crossing_samples;
    struct EmittedSegmentBounds {
        size_t index;
        Vec2d minimum;
        Vec2d maximum;
        double minimum_z;
        double maximum_z;
    };
    std::vector<EmittedSegmentBounds> wall_segment_bounds;
    wall_segment_bounds.reserve(emitted_wall_segments.size());
    for (size_t index = 0; index < emitted_wall_segments.size(); ++index) {
        const EmittedSegment &segment = emitted_wall_segments[index];
        wall_segment_bounds.push_back({
            index,
            segment.start.head<2>().cwiseMin(segment.end.head<2>()),
            segment.start.head<2>().cwiseMax(segment.end.head<2>()),
            std::min(segment.start.z(), segment.end.z()),
            std::max(segment.start.z(), segment.end.z())});
    }
    std::sort(wall_segment_bounds.begin(), wall_segment_bounds.end(),
        [](const EmittedSegmentBounds &left,
           const EmittedSegmentBounds &right) {
            return left.minimum.x() < right.minimum.x();
        });
    for (size_t left = 0; left < wall_segment_bounds.size(); ++left) {
        const EmittedSegmentBounds &first_bounds = wall_segment_bounds[left];
        const EmittedSegment &first =
            emitted_wall_segments[first_bounds.index];
        const Vec2d first_start = first.start.head<2>();
        const Vec2d first_delta = first.end.head<2>() - first_start;
        for (size_t right = left + 1;
             right < wall_segment_bounds.size(); ++right) {
            const EmittedSegmentBounds &second_bounds =
                wall_segment_bounds[right];
            if (second_bounds.minimum.x() > first_bounds.maximum.x())
                break;
            if (second_bounds.maximum.y() < first_bounds.minimum.y() ||
                second_bounds.minimum.y() > first_bounds.maximum.y() ||
                first_bounds.minimum_z > second_bounds.maximum_z + 0.05 ||
                second_bounds.minimum_z > first_bounds.maximum_z + 0.05)
                continue;
            const EmittedSegment &second =
                emitted_wall_segments[second_bounds.index];
            const Vec2d second_start = second.start.head<2>();
            const Vec2d second_delta = second.end.head<2>() - second_start;
            const double denominator =
                first_delta.x() * second_delta.y() -
                first_delta.y() * second_delta.x();
            const double first_length = first_delta.norm();
            const double second_length = second_delta.norm();
            if (first_length <= EPSILON || second_length <= EPSILON)
                continue;
            double first_z = 0.;
            double second_z = 0.;
            bool collision = false;
            if (std::abs(denominator) <=
                    1e-4 * first_length * second_length) {
                const Vec2d offset = second_start - first_start;
                const double offset_cross =
                    first_delta.x() * offset.y() -
                    first_delta.y() * offset.x();
                if (std::abs(offset_cross) > 0.02 * first_length)
                    continue;
                const double second_start_on_first =
                    offset.dot(first_delta) / first_delta.squaredNorm();
                const double second_end_on_first =
                    (second.end.head<2>() - first_start).dot(first_delta) /
                    first_delta.squaredNorm();
                const double overlap_start = std::max(
                    0.01, std::min(second_start_on_first, second_end_on_first));
                const double overlap_end = std::min(
                    0.99, std::max(second_start_on_first, second_end_on_first));
                if ((overlap_end - overlap_start) * first_length <= 0.05)
                    continue;
                const double first_t = 0.5 * (overlap_start + overlap_end);
                const Vec2d overlap_middle = first_start + first_t * first_delta;
                const double second_t =
                    (overlap_middle - second_start).dot(second_delta) /
                    second_delta.squaredNorm();
                if (second_t <= 0.01 || second_t >= 0.99)
                    continue;
                first_z = first.start.z() +
                    first_t * (first.end.z() - first.start.z());
                second_z = second.start.z() +
                    second_t * (second.end.z() - second.start.z());
                collision = std::abs(first_z - second_z) < 0.05;
            } else {
                const Vec2d offset = second_start - first_start;
                const double first_t =
                    (offset.x() * second_delta.y() -
                     offset.y() * second_delta.x()) / denominator;
                const double second_t =
                    (offset.x() * first_delta.y() -
                     offset.y() * first_delta.x()) / denominator;
                if (first_t <= 0.01 || first_t >= 0.99 ||
                    second_t <= 0.01 || second_t >= 0.99)
                    continue;
                first_z = first.start.z() +
                    first_t * (first.end.z() - first.start.z());
                second_z = second.start.z() +
                    second_t * (second.end.z() - second.start.z());
                collision = std::abs(first_z - second_z) < 0.05;
            }
            wall_wall_same_z_crossings += collision;
            if (collision &&
                wall_wall_crossing_samples.size() < 8)
                wall_wall_crossing_samples.emplace_back(
                    "(" + std::to_string(first.start.x()) + "," +
                    std::to_string(first.start.y()) + "," +
                    std::to_string(first.end.x()) + "," +
                    std::to_string(first.end.y()) + "," +
                    std::to_string(first_z) + ",chain=" +
                    std::to_string(first.chain) + "," + first.type + "," +
                    first.comment + ") x (" +
                    std::to_string(second.start.x()) + "," +
                    std::to_string(second.start.y()) + "," +
                    std::to_string(second.end.x()) + "," +
                    std::to_string(second.end.y()) + "," +
                    std::to_string(second_z) + ",chain=" +
                    std::to_string(second.chain) + "," + second.type + "," +
                    second.comment + ")");
        }
    }
    CAPTURE(wall_wall_crossing_samples);
    const auto gcode_context = [&generated_gcode](const std::string &needle) {
        const size_t position = generated_gcode.find(needle);
        if (position == std::string::npos)
            return std::string{};
        const size_t begin = position > 500 ? position - 500 : 0;
        return generated_gcode.substr(begin, 1000);
    };
    CHECK(wall_wall_same_z_crossings == 0);

    // Conventional bridge/infill centerlines must also be removed where the
    // draped skin replaces them. A retained bridge crossing the first roof
    // raster at the same physical Z is a real nozzle collision even though it
    // is not classified as a wall in preview.
    size_t conventional_top_crossings = 0;
    for (const EmittedSegment &conventional : emitted_conventional_segments) {
        const Vec2d conventional_start = conventional.start.head<2>();
        const Vec2d conventional_delta =
            conventional.end.head<2>() - conventional_start;
        for (const EmittedSegment &top : emitted_top_segments) {
            const Vec2d top_start = top.start.head<2>();
            const Vec2d top_delta = top.end.head<2>() - top_start;
            // The exact intersection test below is intentionally strict, but
            // almost every Benchy segment pair is spatially disjoint. Reject
            // those pairs before doing the line solve so this regression test
            // remains practical as non-planar rasters gain resolution.
            if (std::max(conventional.start.x(), conventional.end.x()) <=
                    std::min(top.start.x(), top.end.x()) ||
                std::max(top.start.x(), top.end.x()) <=
                    std::min(conventional.start.x(), conventional.end.x()) ||
                std::max(conventional.start.y(), conventional.end.y()) <=
                    std::min(top.start.y(), top.end.y()) ||
                std::max(top.start.y(), top.end.y()) <=
                    std::min(conventional.start.y(), conventional.end.y()))
                continue;
            const double denominator =
                conventional_delta.x() * top_delta.y() -
                conventional_delta.y() * top_delta.x();
            if (std::abs(denominator) <= EPSILON)
                continue;
            const Vec2d offset = top_start - conventional_start;
            const double conventional_t =
                (offset.x() * top_delta.y() -
                 offset.y() * top_delta.x()) / denominator;
            const double top_t =
                (offset.x() * conventional_delta.y() -
                 offset.y() * conventional_delta.x()) / denominator;
            if (conventional_t <= 0.01 || conventional_t >= 0.99 ||
                top_t <= 0.01 || top_t >= 0.99)
                continue;
            const double conventional_z = conventional.start.z() +
                conventional_t *
                    (conventional.end.z() - conventional.start.z());
            const double top_z = top.start.z() +
                top_t * (top.end.z() - top.start.z());
            conventional_top_crossings +=
                std::abs(conventional_z - top_z) < 0.05;
        }
    }
    CHECK(conventional_top_crossings == 0);

    // Boundary repair may create tiny variable-width tips beside a compact
    // feature. Distinct top-surface chains must not cross at the same Z; that
    // would deposit the same roof volume twice on every transition course.
    size_t top_top_crossings = 0;
    for (size_t left = 0; left < emitted_top_segments.size(); ++left) {
        const EmittedSegment &first = emitted_top_segments[left];
        const Vec2d first_start = first.start.head<2>();
        const Vec2d first_delta = first.end.head<2>() - first_start;
        for (size_t right = left + 1; right < emitted_top_segments.size(); ++right) {
            const EmittedSegment &second = emitted_top_segments[right];
            const Vec2d second_start = second.start.head<2>();
            const Vec2d second_delta = second.end.head<2>() - second_start;
            const double denominator = first_delta.x() * second_delta.y() -
                                       first_delta.y() * second_delta.x();
            if (std::abs(denominator) <= EPSILON)
                continue;
            const Vec2d offset = second_start - first_start;
            const double first_t =
                (offset.x() * second_delta.y() -
                 offset.y() * second_delta.x()) / denominator;
            const double second_t =
                (offset.x() * first_delta.y() -
                 offset.y() * first_delta.x()) / denominator;
            if (first_t <= 0.01 || first_t >= 0.99 ||
                second_t <= 0.01 || second_t >= 0.99)
                continue;
            const double first_z = first.start.z() +
                first_t * (first.end.z() - first.start.z());
            const double second_z = second.start.z() +
                second_t * (second.end.z() - second.start.z());
            top_top_crossings += std::abs(first_z - second_z) < 0.05;
        }
    }
    CHECK(top_top_crossings == 0);

    // Dependency-aware feature deferral is atomic at the wall-course level.
    // Each transformed chimney course must retain the configured three-wall
    // topology, including inset zero. Accepting a course after projecting only
    // one child path is what made the outer wall disappear at the roof base.
    struct FeatureCourseWalls {
        std::set<unsigned int> insets;
        bool has_external {false};
        bool has_detached_open_fragment {false};
    };
    std::map<size_t, FeatureCourseWalls> feature_course_walls;
    const auto collect_feature_courses =
        [&feature_course_walls](auto &&self,
                                const ExtrusionEntity &entity,
                                bool connected_parent) -> void {
            const auto collect_path = [&feature_course_walls](
                                          const ExtrusionPath &path,
                                          bool path_has_connected_parent) {
                if (!path.nonplanar_feature_transition ||
                    !is_perimeter(path.role()))
                    return;
                FeatureCourseWalls &course =
                    feature_course_walls[path.nonplanar_feature_course];
                course.insets.insert(path.inset_idx);
                course.has_external |= path.role() == erExternalPerimeter;
                // Variable-width Arachne rings contain several path
                // fragments. Feature deferral must retain their parent loop
                // or multipath; appending the open children independently
                // lets the sorter insert travels/retracts through the middle
                // of one configured chimney wall.
                course.has_detached_open_fragment |=
                    !path_has_connected_parent &&
                    path.polyline.points.size() >= 2 &&
                    path.polyline.points.front().to_point().distance_to(
                        path.polyline.points.back().to_point()) > scale_(0.05);
            };
            if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity))
                collect_path(*path, connected_parent);
            else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity))
                for (const ExtrusionPath &path : loop->paths)
                    collect_path(path, true);
            else if (const auto *multipath =
                         dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
                for (const ExtrusionPath &path : multipath->paths)
                    collect_path(path, true);
                const auto feature = std::find_if(
                    multipath->paths.begin(), multipath->paths.end(),
                    [](const ExtrusionPath &path) {
                        return path.nonplanar_feature_transition &&
                               is_perimeter(path.role());
                    });
                if (feature != multipath->paths.end() &&
                    !multipath->paths.empty() &&
                    multipath->paths.front().first_point().distance_to(
                        multipath->paths.back().last_point()) > scale_(0.05))
                    feature_course_walls[feature->nonplanar_feature_course]
                        .has_detached_open_fragment = true;
            }
            else if (const auto *collection =
                         dynamic_cast<const ExtrusionEntityCollection *>(&entity))
                for (const ExtrusionEntity *child : collection->entities)
                    self(self, *child, connected_parent);
        };
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions())
                for (const ExtrusionEntityCollection *collection :
                    {&region->perimeters, &region->fills})
                    for (const ExtrusionEntity *entity : collection->entities)
                        collect_feature_courses(
                            collect_feature_courses, *entity, false);
    REQUIRE_FALSE(feature_course_walls.empty());
    for (const auto &[course_index, walls] : feature_course_walls) {
        CAPTURE(course_index, walls.insets);
        CHECK(walls.has_external);
        CHECK(walls.insets.count(0) == 1);
        CHECK(walls.insets.count(1) == 1);
        CHECK(walls.insets.count(2) == 1);
        CHECK_FALSE(walls.has_detached_open_fragment);
    }

    // Successive deferred feature courses must retain a printable physical
    // separation everywhere around the wall. Blending every chimney ring
    // independently toward the same roof surface collapsed the low side to
    // roughly 0.03 mm while retaining nominal 0.25 mm wall flow. In preview
    // those rings crossed one another; on a printer the nozzle would drag
    // through the preceding course. Compare corresponding external-wall
    // centerlines in XY and require at least half a nominal layer of Z
    // separation (the transition may intentionally vary local thickness).
    struct FeatureWallSegment {
        Vec3d start;
        Vec3d end;
    };
    std::map<size_t, std::vector<FeatureWallSegment>> feature_outer_courses;
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions()) {
                const auto collect = [&](auto &&self,
                                         const ExtrusionEntity &entity) -> void {
                    const auto path = [&](const ExtrusionPath &value) {
                        if (!value.nonplanar_feature_transition ||
                            value.role() != erExternalPerimeter)
                            return;
                        const Points3 &points = value.polyline.points;
                        for (size_t point_idx = 1; point_idx < points.size(); ++point_idx) {
                            const auto absolute = [layer](const Point3 &point) {
                                return Vec3d(unscale_(point.x()), unscale_(point.y()),
                                    layer->print_z + unscale_(point.z()));
                            };
                            feature_outer_courses[value.nonplanar_feature_course].push_back(
                                {absolute(points[point_idx - 1]),
                                 absolute(points[point_idx])});
                        }
                    };
                    if (const auto *value = dynamic_cast<const ExtrusionPath *>(&entity))
                        path(*value);
                    else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity))
                        for (const ExtrusionPath &value : loop->paths)
                            path(value);
                    else if (const auto *multipath =
                                 dynamic_cast<const ExtrusionMultiPath *>(&entity))
                        for (const ExtrusionPath &value : multipath->paths)
                            path(value);
                    else if (const auto *collection =
                                 dynamic_cast<const ExtrusionEntityCollection *>(&entity))
                        for (const ExtrusionEntity *child : collection->entities)
                            self(self, *child);
                };
                for (const ExtrusionEntityCollection *collection :
                     {&region->perimeters, &region->fills})
                    for (const ExtrusionEntity *entity : collection->entities)
                        collect(collect, *entity);
            }
    REQUIRE(feature_outer_courses.size() > 1);
    double minimum_feature_course_separation =
        std::numeric_limits<double>::max();
    for (auto current = std::next(feature_outer_courses.begin());
         current != feature_outer_courses.end(); ++current) {
        const auto previous = std::prev(current);
        for (const FeatureWallSegment &segment : current->second) {
            const Vec3d midpoint = 0.5 * (segment.start + segment.end);
            for (const FeatureWallSegment &support : previous->second) {
                const Vec2d delta = support.end.head<2>() - support.start.head<2>();
                const double denominator = delta.squaredNorm();
                const double t = denominator <= EPSILON ? 0. : std::clamp(
                    (midpoint.head<2>() - support.start.head<2>()).dot(delta) /
                        denominator,
                    0., 1.);
                const Vec2d projected = support.start.head<2>() + t * delta;
                if ((projected - midpoint.head<2>()).norm() > 0.05)
                    continue;
                const double support_z = support.start.z() +
                    t * (support.end.z() - support.start.z());
                minimum_feature_course_separation = std::min(
                    minimum_feature_course_separation, midpoint.z() - support_z);
            }
        }
    }
    CAPTURE(minimum_feature_course_separation);
    CHECK(minimum_feature_course_separation >=
          0.5 * config.opt_float("layer_height"));

    // Once the feature has levelled back to normal layers, every complete
    // chimney layer must still contain all three walls. Brick mode changes an
    // inner wall's Z/order metadata; it may not punch intermittent holes in
    // the annulus or suppress its dimensional outer wall.
    size_t complete_chimney_layers = 0;
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers()) {
            if (layer->print_z < 18. || layer->print_z > 20.)
                continue;
            std::set<unsigned int> insets;
            bool has_external = false;
            bool has_brick = false;
            const auto collect = [&](auto &&self,
                                     const ExtrusionEntity &entity) -> void {
                const auto path = [&](const ExtrusionPath &value) {
                    if (!is_perimeter(value.role()))
                        return;
                    insets.insert(value.inset_idx);
                    has_external |= value.role() == erExternalPerimeter;
                    has_brick |= value.staggered_perimeter;
                };
                if (const auto *value =
                        dynamic_cast<const ExtrusionPath *>(&entity))
                    path(*value);
                else if (const auto *loop =
                             dynamic_cast<const ExtrusionLoop *>(&entity))
                    for (const ExtrusionPath &value : loop->paths)
                        path(value);
                else if (const auto *multipath =
                             dynamic_cast<const ExtrusionMultiPath *>(&entity))
                    for (const ExtrusionPath &value : multipath->paths)
                        path(value);
                else if (const auto *collection =
                             dynamic_cast<const ExtrusionEntityCollection *>(&entity))
                    for (const ExtrusionEntity *child : collection->entities)
                        self(self, *child);
            };
            for (const LayerRegion *region : layer->regions())
                for (const ExtrusionEntity *entity : region->perimeters.entities)
                    collect(collect, *entity);
            CAPTURE(layer->print_z, insets);
            CHECK(has_external);
            CHECK(insets.count(0) == 1);
            CHECK(insets.count(1) == 1);
            CHECK(insets.count(2) == 1);
            if (perimeter_mode == "brick")
                CHECK(has_brick);
            ++complete_chimney_layers;
        }
    CHECK(complete_chimney_layers >= 6);

    // Compare emitted 3D extrusion segments, not only the in-memory entities.
    // Re-emitting a retained roof wall underneath its non-planar replacement
    // produces an identical segment (possibly reversed) and is unprintable.
    std::map<std::array<long long, 6>, std::string> emitted_segments;
    size_t duplicate_roof_wall_segments = 0;
    std::string feature_type;
    GCodeReader reader;
    reader.parse_buffer(generated_gcode,
        [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
            const std::string comment(line.comment());
            if (comment.rfind("TYPE:", 0) == 0)
                feature_type = comment.substr(5);
            if (!line.extruding(self) || line.dist_XY(self) <= 0.01)
                return;
            const bool is_wall = feature_type == "Outer wall" ||
                                 feature_type == "Inner wall" ||
                                 feature_type == "Overhang wall" ||
                                 feature_type == "Brick wall" ||
                                 feature_type == "Interlocking inner wall" ||
                                 feature_type == "Non-planar transition" ||
                                 feature_type == "Non-planar top surface";
            if (!is_wall)
                return;
            const Vec3d start(self.x(), self.y(), self.z());
            const Vec3d end(line.new_X(self), line.new_Y(self), line.new_Z(self));
            if (std::max(start.z(), end.z()) < 10.5)
                return;
            const auto quantize = [](double value) {
                return std::llround(value * 1000.);
            };
            std::array<long long, 3> a{quantize(start.x()), quantize(start.y()),
                                       quantize(start.z())};
            std::array<long long, 3> b{quantize(end.x()), quantize(end.y()),
                                       quantize(end.z())};
            if (b < a)
                std::swap(a, b);
            const std::array<long long, 6> key{
                a[0], a[1], a[2], b[0], b[1], b[2]};
            const auto [found, inserted] = emitted_segments.try_emplace(key, feature_type);
            if (!inserted && found->second != feature_type) {
                const bool previous_nonplanar =
                    found->second.rfind("Non-planar", 0) == 0;
                const bool current_nonplanar =
                    feature_type.rfind("Non-planar", 0) == 0;
                duplicate_roof_wall_segments += previous_nonplanar || current_nonplanar;
            }
        });
    CHECK(duplicate_roof_wall_segments == 0);
}

TEST_CASE("Detailed toolhead geometry limits non-planar surface angle", "[Fill][NonplanarSurface]")
{
    Print print;
    Slic3r::Test::init_and_process_print(
        {shallow_top_wedge(8.)}, print,
        {{"nonplanar_top_surface", 1},
         {"nonplanar_top_surface_max_angle", 45},
         {"sequential_print_gantry_geometry",
          R"({"slices":[{"height":"1","type":"convex","polygons":["-10,-10;10,-10;10,10;-10,10"]}]})"},
         {"zaa_enabled", 0},
         {"layer_height", 0.2}});

    // atan(1 / sqrt(10^2 + 10^2)) is about four degrees, so the requested
    // eight-degree skin cannot clear this physical toolhead definition.
    CHECK(nonplanar_top_paths(print).empty());
    // Rejection is a fallback to the original supported top skin, not a
    // destructive removal of its candidate paths. Otherwise internal solid
    // infill becomes the visible surface wherever a gantry vetoes a patch.
    bool has_planar_top_solid = false;
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions())
                for (const ExtrusionEntity *entity : region->fills.flatten().entities)
                    if (const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
                        path != nullptr && !path->nonplanar_surface &&
                        path->role() == erTopSolidInfill)
                        has_planar_top_solid = true;
    CHECK(has_planar_top_solid);
}

TEST_CASE("Hybrid top surfaces Z contour only paths rejected by non-planar clearance",
          "[Fill][NonplanarSurface][ContourZ][Regression]")
{
    Print rejected;
    Slic3r::Test::init_and_process_print(
        {shallow_top_wedge(8.)}, rejected,
        {{"top_surface_z_mode", "nonplanar_with_z_contouring_fallback"},
         {"nonplanar_top_surface", true}, {"zaa_enabled", true},
         {"nonplanar_top_surface_max_angle", 45},
         {"sequential_print_gantry_geometry",
          R"({"slices":[{"height":"1","type":"convex","polygons":["-10,-10;10,-10;10,10;-10,10"]}]})"},
         {"layer_height", 0.2}});

    CHECK(nonplanar_top_paths(rejected).empty());
    size_t fallback_paths = 0;
    for (const PrintObject *object : rejected.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions())
                for (const ExtrusionEntity *entity : region->fills.flatten().entities)
                    if (const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
                        path != nullptr && !path->nonplanar_surface &&
                        path->role() == erTopSolidInfill && path->z_contoured)
                        ++fallback_paths;
    CHECK(fallback_paths > 0);

    Print accepted;
    Slic3r::Test::init_and_process_print(
        {shallow_top_wedge(8.)}, accepted,
        {{"top_surface_z_mode", "nonplanar_with_z_contouring_fallback"},
         {"nonplanar_top_surface", true}, {"zaa_enabled", true},
         {"nonplanar_top_surface_max_angle", 45},
         {"wall_loops", 3}, {"layer_height", 0.2}});
    CHECK_FALSE(nonplanar_top_paths(accepted).empty());
    CHECK_FALSE(has_planar_path_inside_nonplanar_top_coverage(accepted));
}

TEST_CASE("Ironing follows eligible non-planar top surfaces", "[Fill][NonplanarSurface]")
{
    Print print;
    Slic3r::Test::init_and_process_print(
        {shallow_top_wedge(8.)}, print,
        {{"nonplanar_top_surface", 1},
         {"nonplanar_top_surface_max_angle", 45},
         {"nonplanar_top_surface_layers", 5},
         {"ironing_type", "top"},
         {"top_surface_pattern", "rectilinear"},
         {"layer_height", 0.2}});

    CHECK(has_nonplanar_path_with_role(print, erIroning));
}

TEST_CASE("Ironing follows the solid infill rotation template", "[Fill]")
{
    Print print;
    Slic3r::Test::init_and_process_print({Slic3r::Test::cube(10)}, print,
                                        {{"solid_infill_rotate_template", "+45"},
                                         {"internal_solid_infill_pattern", "monotonic"},
                                         {"top_surface_pattern", "monotonic"},
                                         // Every solid surface, so the comparison covers every layer.
                                         {"ironing_type", "solid"},
                                         {"sparse_infill_density", "100%"},
                                         {"ironing_angle", 0},
                                         {"ironing_angle_fixed", 0},
                                         {"layer_height", 0.2}});

    const std::vector<int> ironing = angles_per_layer(print, ironing_role);
    const std::vector<int> solid   = angles_per_layer(print, solid_role);
    REQUIRE(ironing.size() == solid.size());

    // With no fixed angle and no offset, ironing runs along the template's angle for that layer.
    int compared = 0;
    for (size_t i = 0; i < ironing.size(); ++i)
        if (ironing[i] != -1 && solid[i] != -1) {
            CAPTURE(i, ironing[i], solid[i]);
            CHECK(ironing[i] == solid[i]);
            ++compared;
        }
    // Most of the object, not one lucky layer.
    REQUIRE(compared > int(ironing.size()) / 2);
}


namespace {

PrintRegionConfig ironing_config(IroningType type,
                                 int top_surface_filament_id = 1,
                                 int top_shell_layers        = 3,
                                 int bottom_shell_layers     = 1)
{
    PrintRegionConfig cfg;
    cfg.ironing_type.value            = type;
    cfg.top_surface_filament_id.value = top_surface_filament_id;
    cfg.top_shell_layers.value        = top_shell_layers;
    cfg.bottom_shell_layers.value     = bottom_shell_layers;
    cfg.outer_wall_filament_id.value  = 1;
    cfg.wall_loops.value              = 2;
    return cfg;
}

} // namespace

TEST_CASE("Ironing an all-solid region uses the top surface filament on every layer", "[Fill]")
{
    const PrintRegionConfig cfg = ironing_config(IroningType::AllSolid, /*top_surface_filament_id=*/2);
    const bool is_topmost_layer = GENERATE(false, true);
    CAPTURE(is_topmost_layer);
    REQUIRE(Layer::choose_ironing_extruder(cfg, /*spiral_mode=*/false, is_topmost_layer) == 2);
}

TEST_CASE("Ironing top surfaces uses the top surface filament when the region has top shells", "[Fill]")
{
    const PrintRegionConfig cfg = ironing_config(IroningType::TopSurfaces,
                                                 /*top_surface_filament_id=*/3,
                                                 /*top_shell_layers=*/2);
    REQUIRE(Layer::choose_ironing_extruder(cfg, /*spiral_mode=*/false, /*is_topmost_layer=*/false) == 3);
}

TEST_CASE("Ironing top surfaces without top shells needs spiral mode and more than one bottom shell", "[Fill]")
{
    const PrintRegionConfig one_bottom_shell = ironing_config(IroningType::TopSurfaces,
                                                              /*top_surface_filament_id=*/1,
                                                              /*top_shell_layers=*/0,
                                                              /*bottom_shell_layers=*/1);
    const PrintRegionConfig two_bottom_shells = ironing_config(IroningType::TopSurfaces,
                                                               /*top_surface_filament_id=*/1,
                                                               /*top_shell_layers=*/0,
                                                               /*bottom_shell_layers=*/2);

    REQUIRE(Layer::choose_ironing_extruder(two_bottom_shells, /*spiral_mode=*/true, /*is_topmost_layer=*/false) == 1);
    REQUIRE(Layer::choose_ironing_extruder(one_bottom_shell, /*spiral_mode=*/true, /*is_topmost_layer=*/false) == -1);
    REQUIRE(Layer::choose_ironing_extruder(two_bottom_shells, /*spiral_mode=*/false, /*is_topmost_layer=*/false) == -1);
}

TEST_CASE("Ironing the topmost surface only applies to the topmost layer", "[Fill]")
{
    const PrintRegionConfig cfg = ironing_config(IroningType::TopmostOnly, /*top_surface_filament_id=*/4);
    REQUIRE(Layer::choose_ironing_extruder(cfg, /*spiral_mode=*/false, /*is_topmost_layer=*/true) == 4);
    REQUIRE(Layer::choose_ironing_extruder(cfg, /*spiral_mode=*/false, /*is_topmost_layer=*/false) == -1);
}

TEST_CASE("A region with ironing turned off is never ironed", "[Fill]")
{
    const PrintRegionConfig cfg = ironing_config(IroningType::NoIroning);
    const bool spiral_mode      = GENERATE(false, true);
    CAPTURE(spiral_mode);
    REQUIRE(Layer::choose_ironing_extruder(cfg, spiral_mode, /*is_topmost_layer=*/true) == -1);
}

TEST_CASE("Solid infill direction offsets every layer when no template is set", "[Fill]")
{
    auto angles_for = [](int direction) {
        Print print;
        Slic3r::Test::init_and_process_print({Slic3r::Test::cube(10)}, print,
                                            {{"solid_infill_direction", direction},
                                             {"sparse_infill_density", "100%"},
                                             {"internal_solid_infill_pattern", "monotonic"},
                                             {"layer_height", 0.2}});
        return angles_per_layer(print, solid_role);
    };

    const std::vector<int> at_0  = angles_for(0);
    const std::vector<int> at_30 = angles_for(30);
    REQUIRE(at_0.size() == at_30.size());
    REQUIRE(std::count(at_0.begin(), at_0.end(), -1) == 0);

    for (size_t i = 0; i < at_0.size(); ++i) {
        const int delta = ((at_30[i] - at_0[i]) % 180 + 180) % 180;
        CAPTURE(i, at_0[i], at_30[i]);
        CHECK(delta == 30);
    }
}

// Orca: the spiral inset pattern chains the concentric loops into a single continuous path per
// island, so it has to cope with the degenerate loops offsetting leaves behind and it must not join
// loops that only look adjacent.
namespace {

Slic3r::Polylines spiral_inset_fill(const Slic3r::ExPolygon &surface_shape, double spacing)
{
    std::unique_ptr<Slic3r::Fill> filler(Slic3r::Fill::new_from_type("spiralinset"));
    filler->spacing = spacing;
    // Cancel the half-spacing contraction fill_surface() applies, so the filler sees the shape as given.
    filler->overlap = 0.5 * spacing;

    Slic3r::FillParams fill_params;
    fill_params.density     = 1.f;
    fill_params.dont_adjust = true;

    Slic3r::Surface surface(Slic3r::stBottom, surface_shape);
    return filler->fill_surface(&surface, fill_params);
}

Slic3r::ExPolygon rectangle(double x, double y, double w, double h)
{
    return Slic3r::ExPolygon({Slic3r::Point::new_scale(x, y), Slic3r::Point::new_scale(x + w, y),
                              Slic3r::Point::new_scale(x + w, y + h), Slic3r::Point::new_scale(x, y + h)});
}

// Area of the surface the toolpaths fail to cover, and the largest single patch of it, in mm2. Each
// bead is measured at its own width so the variable width walls are not sold short.
std::pair<double, double> uncovered_area(const Slic3r::ExPolygon &surface_shape, const Slic3r::Polygons &covered)
{
    double total = 0, biggest = 0;
    for (const Slic3r::ExPolygon &gap : Slic3r::diff_ex(Slic3r::ExPolygons{surface_shape}, Slic3r::union_(covered))) {
        const double area = unscale<double>(unscale<double>(gap.area()));
        total += area;
        biggest = std::max(biggest, area);
    }
    return {total, biggest};
}

Slic3r::Polygons beads_of(const Slic3r::Polylines &paths, double width)
{
    return Slic3r::offset(paths, float(scale_(0.5 * width)));
}

Slic3r::Polygons beads_of(const Slic3r::ThickPolylines &paths)
{
    Slic3r::Polygons covered;
    for (const Slic3r::ThickPolyline &path : paths)
        for (size_t i = 0; i + 1 < path.points.size(); ++i) {
            Slic3r::Polyline segment;
            segment.points = {path.points[i], path.points[i + 1]};
            Slic3r::append(covered, Slic3r::offset(Slic3r::Polylines{segment},
                                                  float(0.5 * std::max(path.width[2 * i], path.width[2 * i + 1]))));
        }
    return covered;
}

} // namespace

TEST_CASE("Spiral inset fill drops loops shorter than the end clipping", "[Fill][Regression]")
{
    // A sliver whose whole perimeter is shorter than the length clipped off the end of a loop, so the
    // clipping consumes the path entirely. Such a loop carries no extrusion and must be dropped
    // rather than kept as an empty path and read back from.
    const double spacing = 0.45;

    Slic3r::Polylines paths;
    REQUIRE_NOTHROW(paths = spiral_inset_fill(rectangle(0, 0, 0.05, 0.05), spacing));
    for (const Slic3r::Polyline &path : paths)
        CHECK(path.size() >= 2);

    // The same surface at a size the clipping cannot swallow still gets filled.
    REQUIRE_NOTHROW(paths = spiral_inset_fill(rectangle(0, 0, 5, 5), spacing));
    REQUIRE(paths.size() == 1);
    CHECK(paths.front().size() >= 2);
}

TEST_CASE("Spiral inset fill keeps separate islands on separate paths", "[Fill]")
{
    // Two lobes joined by a neck narrower than the loop spacing: the inward offsets break the surface
    // into two islands, which cannot share one spiral, and no path may leave the surface.
    const double spacing = 0.45;
    Slic3r::ExPolygon dumbbell = rectangle(0, 0, 6, 6);
    dumbbell = Slic3r::union_ex(Slic3r::ExPolygons{dumbbell, rectangle(6, 2.9, 4, 0.2), rectangle(10, 0, 6, 6)}).front();

    const Slic3r::Polylines paths = spiral_inset_fill(dumbbell, spacing);
    REQUIRE(paths.size() >= 2);

    // Inflate by a hair so that loops sitting exactly on the outline still count as contained.
    const Slic3r::ExPolygons within = Slic3r::offset_ex(dumbbell, float(SCALED_EPSILON));
    REQUIRE(within.size() == 1);
    for (const Slic3r::Polyline &path : paths) {
        CHECK(path.size() >= 2);
        CHECK(within.front().contains(path));
    }
}


TEST_CASE("Spiral inset fill stays connected across sharp corners", "[Fill][Regression]")
{
    // At a corner of half-angle a, the next ring inward retreats along the bisector by spacing/sin(a),
    // which leaves it several spacings from the end of the ring it continues. Judging the break by
    // distance broke the spiral into loose rings at every spike; nesting is what decides the island.
    const double spacing = 0.45;
    const Slic3r::ExPolygon spike({Slic3r::Point::new_scale(0, 0), Slic3r::Point::new_scale(30, 0),
                                   Slic3r::Point::new_scale(15, 4)});

    const Slic3r::Polylines paths = spiral_inset_fill(spike, spacing);
    CHECK(paths.size() == 1);

    const Slic3r::ExPolygons within = Slic3r::offset_ex(spike, float(SCALED_EPSILON));
    REQUIRE(within.size() == 1);
    for (const Slic3r::Polyline &path : paths)
        CHECK(within.front().contains(path));
}

TEST_CASE("Spiral inset fill starts on a convex corner", "[Fill][Regression]")
{
    // The only right angle on this outline is the reflex one: the two edges meeting at the origin
    // span 90 degrees exactly as a square corner would, but the material lies outside them. The next
    // ring in steps away from a reflex corner along the bisector instead of hugging it, so starting
    // the spiral there sent it across a long diagonal on every single ring.
    const double spacing = 0.45;
    const Slic3r::ExPolygon notched({Slic3r::Point::new_scale(0, 0), Slic3r::Point::new_scale(0, 10),
                                     Slic3r::Point::new_scale(-16, 18), Slic3r::Point::new_scale(-16, -2),
                                     Slic3r::Point::new_scale(-8, -16), Slic3r::Point::new_scale(18, -16),
                                     Slic3r::Point::new_scale(10, 0)});

    const Slic3r::Polylines paths = spiral_inset_fill(notched, spacing);
    REQUIRE(paths.size() >= 1);

    // Every edge of the outline is at least 45 degrees off the bisector of that reflex corner, and
    // so is every ring offset from it. A long segment running along the bisector can therefore only
    // be the spiral striking out across the rings to reach the next one.
    for (const Slic3r::Polyline &path : paths)
        for (const Slic3r::Line &segment : path.lines()) {
            const Vec2d  v         = (segment.b - segment.a).cast<double>();
            const double direction = std::fmod(std::atan2(v.y(), v.x()) * 180.0 / M_PI + 180.0, 180.0);
            if (std::abs(direction - 45.0) > 25.0)
                continue;
            CAPTURE(direction, unscale<double>(segment.length()));
            CHECK(segment.length() <= scale_(1.5 * spacing));
        }
}

TEST_CASE("Spiral inset fill closes the gaps with variable width walls", "[Fill]")
{
    // Fixed width loops cannot fill a region that is not a whole number of lines across and leave the
    // remainder open, which on a ring shows up as a wedge several lines wide. Plain concentric avoids
    // that by building solid surfaces out of Arachne's variable width walls, and so must this pattern.
    const double spacing = 0.45;
    Slic3r::ExPolygon ring = rectangle(0, 0, 24, 24);
    Slic3r::Polygon   hole;
    for (int i = 0; i < 64; ++i) {
        const double angle = -2.0 * PI * i / 64.0; // clockwise, so it reads as a hole
        hole.points.emplace_back(Slic3r::Point::new_scale(12 + 7.3 * std::cos(angle), 12 + 7.3 * std::sin(angle)));
    }
    ring.holes.emplace_back(hole);

    Slic3r::PrintConfig       print_config;
    Slic3r::PrintObjectConfig object_config;
    auto make_filler = [&]() {
        std::unique_ptr<Slic3r::Fill> filler(Slic3r::Fill::new_from_type("spiralinset"));
        filler->spacing             = spacing;
        filler->overlap             = 0.5 * spacing; // cancel the contraction, so both see the same surface
        filler->print_config        = &print_config;
        filler->print_object_config = &object_config;
        return filler;
    };

    Slic3r::FillParams params;
    params.density      = 1.f;
    params.dont_adjust  = false;
    params.layer_height = 0.2;

    const Slic3r::Surface surface(Slic3r::stTop, ring);

    std::unique_ptr<Slic3r::Fill> fixed = make_filler();
    const Slic3r::Polylines fixed_width = fixed->fill_surface(&surface, params);
    REQUIRE(!fixed_width.empty());
    const auto fixed_gaps = uncovered_area(ring, beads_of(fixed_width, fixed->spacing));

    params.use_arachne = true;
    std::unique_ptr<Slic3r::Fill> variable = make_filler();
    const Slic3r::ThickPolylines variable_width = variable->fill_surface_arachne(&surface, params);
    REQUIRE(!variable_width.empty());
    const auto variable_gaps = uncovered_area(ring, beads_of(variable_width));

    CAPTURE(fixed_gaps.first, fixed_gaps.second, variable_gaps.first, variable_gaps.second);
    // The wedges the fixed width loops leave behind are what the variable width walls take up.
    CHECK(variable_gaps.second < 0.5 * fixed_gaps.second);
    CHECK(variable_gaps.first < fixed_gaps.first);

    // And it is still a spiral: far fewer paths than the ring has loops.
    // And the walls are still chained into spirals rather than printed one path per wall. The ring is
    // at its narrowest (12 - 7.3) mm across and is filled from both sides, so it is at least this many
    // walls thick there and thicker elsewhere. Arachne's short thin feature walls cannot join a spiral,
    // so only the substantial paths count towards this.
    const size_t walls_across = size_t(2.0 * (12.0 - 7.3) / spacing);
    size_t       spirals      = 0;
    for (const Slic3r::ThickPolyline &path : variable_width)
        if (path.length() > scale_(10.0 * spacing))
            ++spirals;
    CAPTURE(spirals, walls_across, variable_width.size(), fixed_width.size());
    CHECK(2 * spirals < walls_across);
}

TEST_CASE("Honeycomb infill rounds its cell corners with the smooth factor", "[Fill]")
{
    // A cell whose sides are several times the line width, so that the corners have room to be rounded.
    const double spacing = 0.45;
    const double density = 0.1;
    auto         fill    = [spacing, density](double smooth_factor) {
        std::unique_ptr<Slic3r::Fill> filler(Slic3r::Fill::new_from_type("honeycomb"));
        filler->spacing = spacing;

        FillParams params;
        params.density = float(density);
        params.dont_adjust = true;
        // Keep the fragments apart, so that only the turns of the pattern itself are measured.
        params.anchor_length_max = 0.f;
        params.smooth_factor     = smooth_factor;

        Slic3r::ExPolygon square{ Slic3r::Points{
            Point::new_scale(0., 0.), Point::new_scale(50., 0.), Point::new_scale(50., 50.), Point::new_scale(0., 50.) } };
        Slic3r::Surface surface(stInternal, square);
        return filler->fill_surface(&surface, params);
    };

    // Cosine of the sharpest turn of any of the paths, 1 meaning none of them turns at all.
    auto sharpest_turn_cosine = [](const Slic3r::Polylines &polylines) {
        double sharpest = 1.;
        for (const Polyline &polyline : polylines)
            for (size_t i = 1; i + 1 < polyline.size(); ++i) {
                const Vec2d incoming = (polyline[i] - polyline[i - 1]).cast<double>().normalized();
                const Vec2d outgoing = (polyline[i + 1] - polyline[i]).cast<double>().normalized();
                sharpest = std::min(sharpest, incoming.dot(outgoing));
            }
        return sharpest;
    };
    auto point_count = [](const Slic3r::Polylines &polylines) {
        return std::accumulate(polylines.begin(), polylines.end(), size_t(0),
                               [](size_t count, const Polyline &polyline) { return count + polyline.size(); });
    };

    const Slic3r::Polylines sharp  = fill(0.);
    const Slic3r::Polylines smooth = fill(1.);

    REQUIRE(!sharp.empty());
    REQUIRE(smooth.size() == sharp.size());
    REQUIRE(point_count(smooth) > point_count(sharp));
    // The cell corners turn by 60 degrees; smoothing replaces them by gentle curves.
    REQUIRE(sharpest_turn_cosine(sharp) < 0.6);
    REQUIRE(sharpest_turn_cosine(smooth) > 0.9);
}

// Point count, number of turns sharper than 25 degrees and length of the sparse infill of a print.
// A rounded corner is a run of much gentler turns, so smoothing shows up as fewer sharp ones.
struct SparseInfillShape {
    size_t point_count { 0 };
    size_t sharp_turns { 0 };
    size_t path_count { 0 };
    double length { 0. };
    // Digest of every point in the order it is printed. The counts above all survive the same
    // extrusions being joined into different polylines, so only this tells two such fills apart.
    uint64_t sequence { 14695981039346656037ull };
};

static SparseInfillShape sparse_infill_shape(const Print &print)
{
    SparseInfillShape shape;

    auto account = [&shape](const ExtrusionPath &path) {
        if (!sparse_role(path.role()))
            return;
        const Points3 &pts = path.polyline.points;
        ++shape.path_count;
        shape.point_count += pts.size();
        for (const auto &pt : pts)
            for (const coord_t coordinate : {pt.x(), pt.y(), pt.z()})
                shape.sequence = (shape.sequence ^ uint64_t(coordinate)) * 1099511628211ull;
        for (size_t i = 1; i < pts.size(); ++i)
            shape.length += (pts[i] - pts[i - 1]).head<2>().cast<double>().norm();
        for (size_t i = 1; i + 1 < pts.size(); ++i) {
            const Vec2d incoming = (pts[i] - pts[i - 1]).head<2>().cast<double>();
            const Vec2d outgoing = (pts[i + 1] - pts[i]).head<2>().cast<double>();
            if (incoming.squaredNorm() > 0. && outgoing.squaredNorm() > 0. &&
                incoming.normalized().dot(outgoing.normalized()) < 0.9)
                ++shape.sharp_turns;
        }
    };

    for (const Layer *layer : print.objects().front()->layers())
        for (const LayerRegion *region : layer->regions())
            for (const ExtrusionEntity *entity : region->fills.flatten().entities) {
                if (auto *path = dynamic_cast<const ExtrusionPath *>(entity))
                    account(*path);
                else if (auto *multi = dynamic_cast<const ExtrusionMultiPath *>(entity))
                    for (const ExtrusionPath &p : multi->paths)
                        account(p);
                else if (auto *loop = dynamic_cast<const ExtrusionLoop *>(entity))
                    for (const ExtrusionPath &p : loop->paths)
                        account(p);
            }
    return shape;
}

TEST_CASE("Lightning infill slices the same model the same way twice", "[Fill][Regression]")
{
    // Slicing twice in one process catches a generator that carries state from one slice to the
    // next, or whose result depends on how the parallel layer fill interleaves.
    auto shape = [] {
        Print print;
        Slic3r::Test::init_and_process_print({Slic3r::Test::cube(20)}, print,
                                            {{"sparse_infill_pattern", "lightning"},
                                             {"sparse_infill_density", "50%"},
                                             {"layer_height", 0.2}});
        return sparse_infill_shape(print);
    };

    const SparseInfillShape first  = shape();
    const SparseInfillShape second = shape();

    REQUIRE(first.path_count > 0);
    REQUIRE(second.path_count == first.path_count);
    REQUIRE(second.point_count == first.point_count);
    REQUIRE(second.sharp_turns == first.sharp_turns);
    // No tolerance: the same extrusions in the same order add up to the very same number.
    REQUIRE_THAT(second.length, Catch::Matchers::WithinAbs(first.length, 0.));
    // All of the above agree when the same branches are joined into different polylines, so the
    // point sequence is what actually decides whether the two slices produced the same infill.
    REQUIRE(second.sequence == first.sequence);
}

TEST_CASE("Lightning infill rounds the turns of its branches with the smooth factor", "[Fill]")
{
    auto shape_for = [](const std::string &smooth_factor) {
        Print print;
        Slic3r::Test::init_and_process_print({Slic3r::Test::cube(20)}, print,
                                            {{"sparse_infill_pattern", "lightning"},
                                             {"sparse_infill_density", "15%"},
                                             {"sparse_infill_smooth_factor", smooth_factor},
                                             {"layer_height", 0.2}});
        return sparse_infill_shape(print);
    };

    const SparseInfillShape sharp  = shape_for("0%");
    const SparseInfillShape smooth = shape_for("100%");

    REQUIRE(sharp.point_count > 0);
    // The branch turns are replaced by curves, which cut the corners off and take more points to
    // describe. The turns where two branches are joined into one path stay sharp.
    REQUIRE(smooth.point_count > sharp.point_count);
    REQUIRE(smooth.sharp_turns < sharp.sharp_turns);
    REQUIRE(smooth.length < sharp.length);
}

TEST_CASE("Concentric infill rounds its loops with the smooth factor", "[Fill]")
{
    auto shape_for = [](const std::string &smooth_factor) {
        Print print;
        Slic3r::Test::init_and_process_print({Slic3r::Test::cube(20)}, print,
                                            {{"sparse_infill_pattern", "concentric"},
                                             {"sparse_infill_density", "20%"},
                                             {"sparse_infill_smooth_factor", smooth_factor},
                                             {"layer_height", 0.2}});
        return sparse_infill_shape(print);
    };

    const SparseInfillShape sharp  = shape_for("0%");
    const SparseInfillShape smooth = shape_for("100%");

    REQUIRE(sharp.point_count > 0);
    REQUIRE(smooth.point_count > sharp.point_count);
    REQUIRE(smooth.sharp_turns < sharp.sharp_turns);
    REQUIRE(smooth.length < sharp.length);
}

TEST_CASE("Cross hatch infill rounds its transition layers with the smooth factor", "[Fill]")
{
    auto shape_for = [](const std::string &smooth_factor) {
        Print print;
        Slic3r::Test::init_and_process_print({Slic3r::Test::cube(20)}, print,
                                            {{"sparse_infill_pattern", "crosshatch"},
                                             {"sparse_infill_density", "20%"},
                                             {"sparse_infill_smooth_factor", smooth_factor},
                                             {"layer_height", 0.2}});
        return sparse_infill_shape(print);
    };

    const SparseInfillShape sharp  = shape_for("0%");
    const SparseInfillShape smooth = shape_for("100%");

    REQUIRE(sharp.point_count > 0);
    REQUIRE(smooth.point_count > sharp.point_count);
    REQUIRE(smooth.sharp_turns < sharp.sharp_turns);
    REQUIRE(smooth.length < sharp.length);
}

TEST_CASE("Trapezoidal grid infill rounds its corners only with more than one line", "[Fill]")
{
    auto shape_for = [](int multiline, const std::string &smooth_factor) {
        Print print;
        Slic3r::Test::init_and_process_print({Slic3r::Test::cube(20)}, print,
                                            {{"sparse_infill_pattern", "grid"},
                                             {"sparse_infill_density", "20%"},
                                             {"fill_multiline", multiline},
                                             {"sparse_infill_smooth_factor", smooth_factor},
                                             {"layer_height", 0.2}});
        return sparse_infill_shape(print);
    };

    const SparseInfillShape sharp  = shape_for(2, "0%");
    const SparseInfillShape smooth = shape_for(2, "100%");

    REQUIRE(sharp.point_count > 0);
    REQUIRE(smooth.point_count > sharp.point_count);
    REQUIRE(smooth.sharp_turns < sharp.sharp_turns);
    REQUIRE(smooth.length < sharp.length);

    // A single line per infill wall is the plain crossing line grid, which has no corner of its own.
    const SparseInfillShape single_sharp  = shape_for(1, "0%");
    const SparseInfillShape single_smooth = shape_for(1, "100%");
    REQUIRE(single_sharp.point_count > 0);
    REQUIRE(single_smooth.point_count == single_sharp.point_count);
    REQUIRE(single_smooth.length == single_sharp.length);
}

TEST_CASE("3D honeycomb infill rounds its octahedral waves with the smooth factor", "[Fill]")
{
    auto shape_for = [](const std::string &smooth_factor) {
        Print print;
        Slic3r::Test::init_and_process_print({Slic3r::Test::cube(20)}, print,
                                            {{"sparse_infill_pattern", "3dhoneycomb"},
                                             {"sparse_infill_density", "20%"},
                                             {"sparse_infill_smooth_factor", smooth_factor},
                                             {"layer_height", 0.2}});
        return sparse_infill_shape(print);
    };

    const SparseInfillShape sharp  = shape_for("0%");
    const SparseInfillShape smooth = shape_for("100%");

    REQUIRE(sharp.point_count > 0);
    REQUIRE(smooth.point_count > sharp.point_count);
    REQUIRE(smooth.sharp_turns < sharp.sharp_turns);
    REQUIRE(smooth.length < sharp.length);
}

TEST_CASE("Smoothed concentric infill stays inside the fill region", "[Fill][Regression]")
{
    // The concentric loops are offsets of the fill region and are never clipped to it, so a corner
    // rounded across its boundary ends up in a hole or over a wall. Rounding cuts toward the inside of
    // the turn, which leaves the region at every corner of a hole, and in a region thinner than the
    // curve even at a corner turning inwards.
    const bool  thin_region = GENERATE(false, true);
    ExPolygon   region;
    if (thin_region) {
        // An L of two 1.2mm wide arms: cutting the corner they meet at crosses both of them.
        region = ExPolygon{ Slic3r::Points{
            Point::new_scale(0., 0.), Point::new_scale(20., 0.), Point::new_scale(20., 1.2),
            Point::new_scale(1.2, 1.2), Point::new_scale(1.2, 20.), Point::new_scale(0., 20.) } };
    } else {
        region = ExPolygon{ Slic3r::Points{ Point::new_scale(0., 0.), Point::new_scale(50., 0.),
                                            Point::new_scale(50., 50.), Point::new_scale(0., 50.) },
                            Slic3r::Points{ Point::new_scale(30., 20.), Point::new_scale(30., 30.),
                                            Point::new_scale(20., 30.), Point::new_scale(20., 20.) } };
    }
    CAPTURE(thin_region);

    auto fill = [&region](double smooth_factor) {
        std::unique_ptr<Slic3r::Fill> filler(Slic3r::Fill::new_from_type("concentric"));
        filler->spacing = 0.45;

        FillParams params;
        params.density       = 0.1f;
        params.dont_adjust   = true;
        params.smooth_factor = smooth_factor;

        Slic3r::Surface surface(stInternal, region);
        return filler->fill_surface(&surface, params);
    };
    auto point_count = [](const Slic3r::Polylines &polylines) {
        return std::accumulate(polylines.begin(), polylines.end(), size_t(0),
                               [](size_t count, const Polyline &polyline) { return count + polyline.size(); });
    };

    const Slic3r::Polylines sharp  = fill(0.);
    const Slic3r::Polylines smooth = fill(1.);
    REQUIRE(!sharp.empty());

    // Nothing leaves the fill region, which the unrounded loops already touch from the inside.
    const ExPolygons bounds = offset_ex(region, float(SCALED_EPSILON));
    REQUIRE(diff_pl(sharp, bounds).empty());
    REQUIRE(diff_pl(smooth, bounds).empty());
    // The corners that the region has room for are still rounded.
    if (!thin_region)
        REQUIRE(point_count(smooth) > point_count(sharp));
}

TEST_CASE("Smoothing multiline lightning infill keeps its outlines connected", "[Fill][Regression]")
{
    // With more than one line per infill wall, the branches are printed as outlines drawn around them,
    // and the outlines of branches that run close to each other merge into one. Rounding the branches
    // before those outlines are built moves them apart, which breaks the merged outlines up into
    // separate loops - many more of them, each needing its own travel move.
    auto shape_for = [](const std::string &smooth_factor) {
        Print print;
        Slic3r::Test::init_and_process_print({Slic3r::Test::cube(20)}, print,
                                            {{"sparse_infill_pattern", "lightning"},
                                             {"sparse_infill_density", "50%"},
                                             {"fill_multiline", 2},
                                             {"sparse_infill_smooth_factor", smooth_factor},
                                             {"layer_height", 0.2}});
        return sparse_infill_shape(print);
    };

    const SparseInfillShape sharp  = shape_for("0%");
    const SparseInfillShape smooth = shape_for("100%");

    REQUIRE(sharp.path_count > 0);
    // The loop count varies by a loop or two between platforms and between runs, so this is not an
    // exact comparison. Smoothing should leave it about where it was; uncapping the smoothing
    // reach, the regression this guards against, adds about 10%.
    const size_t allowed_extra = sharp.path_count / 50; // 2%
    REQUIRE(smooth.path_count <= sharp.path_count + allowed_extra);
    // The outlines are still rounded.
    REQUIRE(smooth.point_count > sharp.point_count);
    REQUIRE(smooth.sharp_turns < sharp.sharp_turns);
}

TEST_CASE("Sparse plane-path anchors match the printed infill", "[Fill][InternalBridge][Regression]")
{
    // Orca: Compare generated anchors with actual extrusion across plane-path patterns,
    // smoothing, multiline and rotations; an origin shift must not pass as valid support.
    const std::string pattern = GENERATE("hilbertcurve", "octagramspiral", "archimedeanchords");
    const std::string smoothing = GENERATE("0%", "100%");
    const int multiline = GENERATE(1, 2);
    const bool rotated = GENERATE(false, true);
    const bool separated = GENERATE(false, true);
    CAPTURE(pattern, smoothing, multiline, rotated, separated);

    auto config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({{"sparse_infill_pattern", pattern},
                                   {"sparse_infill_density", "15%"},
                                   {"sparse_infill_smooth_factor", smoothing},
                                   {"fill_multiline", multiline},
                                   {"infill_direction", 45},
                                   {"sparse_infill_rotate_template", rotated ? "0,25,50" : ""},
                                   {"align_infill_direction_to_model", rotated},
                                   {"separated_infills", separated},
                                   {"top_shell_layers", 0},
                                   {"bottom_shell_layers", 0},
                                   {"top_shell_thickness", 0},
                                   {"bottom_shell_thickness", 0},
                                   {"layer_height", 0.2},
                                   {"initial_layer_print_height", 0.2},
                                   {"resolution", 0.012}});
    Print print;
    Model model;
    TriangleMesh mesh = make_cube(30, 24, 1);
    if (separated) {
        // Orca: Two disconnected bodies in one object must each use their own infill origin.
        TriangleMesh second = make_cube(30, 24, 1);
        second.translate(50, 0, 0);
        mesh.merge(second);
    }
    Slic3r::Test::init_print({mesh}, print, model, config, nullptr, false);
    if (rotated) {
        model.objects.front()->instances.front()->set_rotation(Vec3d(0., 0., Geometry::deg2rad(23.)));
        print.apply(model, config);
    }
    print.process();

    const Layer &layer = *print.objects().front()->get_layer(4);
    Polylines printed;
    for (const LayerRegion *region : layer.regions())
        for (const ExtrusionEntity *entity : region->fills.flatten().entities)
            if (entity->role() == erInternalInfill)
                entity->collect_polylines(printed);
    REQUIRE_FALSE(printed.empty());
    const AABBTreeLines::LinesDistancer<Line> printed_tree(to_lines(printed));

    // Orca: Exclude perimeter connections: anchoring and extrusion can trim those differently.
    const Polylines anchors = intersection_pl(layer.generate_sparse_infill_polylines_for_anchoring(nullptr, nullptr, nullptr),
                                              shrink(to_polygons(layer.lslices), scale_(3.)));
    REQUIRE_FALSE(anchors.empty());
    double max_distance = 0.;
    for (const Polyline &path : anchors)
        for (const Point &point : path.equally_spaced_points(scale_(0.25)))
            max_distance = std::max(max_distance, printed_tree.distance_from_lines<false>(point));
    // Orca: Allow only the configured simplification tolerance; infill-scale offsets
    // would hide anchors that no longer coincide with printed lines.
    CHECK(unscale<double>(max_distance) <= config.opt_float("resolution"));
}
