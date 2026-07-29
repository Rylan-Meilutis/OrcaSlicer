#include <catch2/catch_all.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/Fill/Fill.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/Geometry.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/SVG.hpp"
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
                    if (entity->role() == erArcOverhang)
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
                        path != nullptr && path->role() == erArcOverhang && path->mm3_per_mm > 0.)
                        return *path;
    return std::nullopt;
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

        for (const Polyline &path : paths) {
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
                            path->role() == erArcOverhang) {
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
        CHECK(gcode.find(";TYPE:Arc overhang") != std::string::npos);

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
                move.extrusion_role == erArcOverhang;
            const bool previous_is_arc =
                previous.type == EMoveType::Extrude &&
                previous.extrusion_role == erArcOverhang &&
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
                             anchor.role == erArcOverhang);
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

        bool arc_overhang_role = false;
        double x = 0.;
        double y = 0.;
        size_t arc_command_count = 0;
        std::istringstream lines(gcode);
        std::string line;
        while (std::getline(lines, line)) {
            if (line.rfind(";TYPE:", 0) == 0) {
                arc_overhang_role =
                    line.find("Arc overhang") != std::string::npos;
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

            if (arc_overhang_role &&
                (command == "G2" || command == "G3") &&
                (has_i || has_j)) {
                ++arc_command_count;
            }
            x = next_x;
            y = next_y;
        }
        CHECK(arc_path_count >= 4);
        CHECK(arc_segment_count > arc_path_count);
        // Global arc fitting remains enabled, but these topology-validated
        // paths must retain their exact G1 polyline through final export.
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
