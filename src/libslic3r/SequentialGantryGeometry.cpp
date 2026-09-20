#include "SequentialGantryGeometry.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <limits>
#include <mutex>
#include <sstream>
#include <type_traits>
#include <unordered_map>

#include <boost/filesystem.hpp>
#include <boost/nowide/fstream.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/regex.hpp>

#include "Config.hpp"
#include "ClipperUtils.hpp"
#include "Geometry/ConvexHull.hpp"
#include "TriangleMesh.hpp"
#include "libslic3r.h"
#include "Utils.hpp"

namespace fs = boost::filesystem;
namespace pt = boost::property_tree;

namespace Slic3r {

namespace {

Polygon parse_polygon(std::string serialized)
{
    std::replace(serialized.begin(), serialized.end(), ';', ' ');
    std::replace(serialized.begin(), serialized.end(), ',', ' ');
    std::stringstream stream(serialized);
    Polygon polygon;
    for (double x = 0., y = 0.; stream >> x >> y;)
        polygon.points.emplace_back(Point::new_scale(x, y));
    return polygon;
}

fs::path resolve_geometry_path(const std::string &serialized)
{
    fs::path path(serialized);
    if (!path.is_absolute() && !fs::exists(path))
        path = fs::path(resources_dir()) / path;
    return path;
}

void normalize_geometry(SequentialGantryGeometry &geometry)
{
    for (SequentialGantrySlice &slice : geometry.slices) {
        slice.polygons.erase(std::remove_if(slice.polygons.begin(), slice.polygons.end(),
            [](const Polygon &polygon) {
                return polygon.points.size() < 3 || std::abs(polygon.area()) <= 1.;
            }), slice.polygons.end());
        for (Polygon &polygon : slice.polygons)
            if (polygon.area() < 0.)
                polygon.reverse();
    }
    geometry.slices.erase(std::remove_if(geometry.slices.begin(), geometry.slices.end(),
        [](const SequentialGantrySlice &slice) { return slice.polygons.empty(); }),
        geometry.slices.end());
    std::stable_sort(geometry.slices.begin(), geometry.slices.end(),
        [](const SequentialGantrySlice &left, const SequentialGantrySlice &right) {
            return left.height < right.height;
        });
}

SequentialGantryGeometry parse_geometry(const pt::ptree &entry, const fs::path &model_directory)
{
    SequentialGantryGeometry geometry;
    if (const auto filename = entry.get_optional<std::string>("gantry_model_filename"))
        geometry.model_path = (model_directory / *filename).string();

    const auto slices = entry.get_child_optional("slices");
    if (!slices)
        return geometry;

    for (const auto &item : *slices) {
        SequentialGantrySlice slice;
        slice.height = item.second.get<double>("height");
        slice.is_box = item.second.get<std::string>("type", "convex") == "box";
        if (const auto polygons = item.second.get_child_optional("polygons")) {
            for (const auto &polygon : *polygons) {
                Polygon parsed = parse_polygon(polygon.second.data());
                if (parsed.points.size() >= 3)
                    slice.polygons.emplace_back(std::move(parsed));
            }
        }
        if (!slice.polygons.empty())
            geometry.slices.emplace_back(std::move(slice));
    }
    normalize_geometry(geometry);
    return geometry;
}

SequentialGantryGeometry load_custom_geometry(const std::string &custom)
{
    if (custom.empty())
        return {};

    pt::ptree root;
    fs::path model_directory;
    const size_t first = custom.find_first_not_of(" \t\r\n");
    if (first != std::string::npos && custom[first] == '{') {
        std::stringstream stream(custom);
        pt::read_json(stream, root);
    } else {
        fs::path path = resolve_geometry_path(custom);
        boost::nowide::ifstream input(path.string());
        pt::read_json(input, root);
        model_directory = path.parent_path();
    }
    SequentialGantryGeometry geometry = parse_geometry(root, model_directory);
    if (geometry.slices.empty() && !geometry.model_path.empty())
        geometry = process_sequential_gantry_model(geometry.model_path);
    return geometry;
}

} // namespace

Vec2d SequentialGantryGeometry::clearance_reach() const
{
    Vec2d reach = Vec2d::Zero();
    for (const auto &slice : slices) {
        if (slice.is_box)
            continue;
        for (const auto &polygon : slice.polygons)
            for (const Point &point : polygon.points)
                reach = reach.cwiseMax(point.cast<double>().cwiseAbs() * SCALING_FACTOR);
    }
    return reach;
}

Polygon sequential_clearance_hull(const Polygon &footprint, const Vec2d &half_extent)
{
    Points points;
    points.reserve(footprint.points.size() * 4);
    const coord_t x = scale_(half_extent.x()), y = scale_(half_extent.y());
    for (const Point &point : footprint.points)
        for (int sx : {-1, 1})
            for (int sy : {-1, 1})
                points.emplace_back(point.x() + sx * x, point.y() + sy * y);
    return Geometry::convex_hull(points);
}

double SequentialGantryGeometry::conservative_clearance_radius() const
{
    double radius = 0.;
    for (const SequentialGantrySlice &slice : slices) {
        if (slice.is_box)
            continue;
        for (const Polygon &polygon : slice.polygons)
            for (const Point &point : polygon.points)
                radius = std::max(radius, unscale_(point.cast<double>().norm()));
    }
    return radius;
}

double SequentialGantryGeometry::maximum_xy_reach(double maximum_height) const
{
    double radius = 0.;
    for (const SequentialGantrySlice &slice : slices) {
        if (slice.height > maximum_height)
            continue;
        for (const Polygon &polygon : slice.polygons)
            for (const Point &point : polygon.points)
                radius = std::max(radius, unscale_(point.cast<double>().norm()));
    }
    return radius;
}

void SequentialGantryGeometry::set_nozzle_tip_diameter(double diameter)
{
    if (diameter <= 0.)
        return;

    const Polygon nozzle_tip = make_circle(scale_(0.5 * diameter), scale_(0.02));
    bool replaced = false;
    for (SequentialGantrySlice &slice : slices) {
        if (std::abs(slice.height) > EPSILON)
            continue;
        // Bundled sequential-arrangement data intentionally uses a very
        // conservative whole-toolhead footprint at the deposition plane.
        // Non-planar clearance needs the physical nozzle tip there; the
        // original higher slices still model the heater block and gantry.
        slice.is_box = false;
        slice.polygons = {nozzle_tip};
        replaced = true;
    }
    if (!replaced)
        slices.insert(slices.begin(), {0., false, {nozzle_tip}});
}

double SequentialGantryGeometry::first_box_height() const
{
    double height = std::numeric_limits<double>::max();
    for (const SequentialGantrySlice &slice : slices)
        if (slice.is_box)
            height = std::min(height, slice.height);
    return height == std::numeric_limits<double>::max() ? 0. : height;
}

std::string SequentialGantryGeometry::validation_error() const
{
    if (slices.empty())
        return "gantry geometry contains no usable collision slices";
    double previous_height = -EPSILON;
    for (const SequentialGantrySlice &slice : slices) {
        if (!std::isfinite(slice.height) || slice.height < 0.)
            return "gantry slice heights must be finite and non-negative";
        if (slice.height + EPSILON < previous_height)
            return "gantry slices must be ordered by height";
        previous_height = slice.height;
        if (slice.polygons.empty())
            return "gantry slice contains no polygons";
        for (const Polygon &polygon : slice.polygons) {
            if (polygon.points.size() < 3 || std::abs(polygon.area()) <= 1.)
                return "gantry slice contains a degenerate polygon";
        }
    }
    return {};
}

SequentialGantryGeometry process_sequential_gantry_model(
    const std::string &serialized_path, double band_height)
{
    SequentialGantryGeometry geometry;
    const fs::path path = resolve_geometry_path(serialized_path);
    geometry.model_path = path.string();
    if (serialized_path.empty() || !fs::exists(path) || !fs::is_regular_file(path))
        return geometry;

    band_height = std::clamp(band_height, 0.25, 5.);
    struct CachedGeometry {
        std::time_t stamp {};
        double band_height {};
        SequentialGantryGeometry geometry;
    };
    static std::mutex cache_mutex;
    static std::unordered_map<std::string, CachedGeometry> cache;
    const std::time_t stamp = fs::last_write_time(path);
    {
        std::lock_guard<std::mutex> lock(cache_mutex);
        const auto found = cache.find(path.string());
        if (found != cache.end() && found->second.stamp == stamp &&
            std::abs(found->second.band_height - band_height) <= EPSILON)
            return found->second.geometry;
    }

    TriangleMesh mesh;
    if (!mesh.ReadSTLFile(path.string().c_str(), true) || mesh.empty())
        return geometry;
    const BoundingBoxf3 bbox = mesh.bounding_box();
    if (!bbox.defined || !std::isfinite(bbox.min.z()) || !std::isfinite(bbox.max.z()) ||
        bbox.max.z() - bbox.min.z() <= EPSILON)
        return geometry;

    const double minimum_z = bbox.min.z();
    const double height = bbox.max.z() - minimum_z;
    const size_t band_count = size_t(std::ceil(height / band_height));
    if (band_count == 0 || band_count > 4000)
        return geometry;

    std::vector<double> sample_z;
    sample_z.reserve(band_count);
    for (size_t band = 0; band < band_count; ++band)
        sample_z.emplace_back(minimum_z + std::min(
            height - 1e-4, (double(band) + 0.5) * band_height));
    const std::vector<ExPolygons> center_slices = mesh.slice(sample_z);
    std::vector<Polygons> band_polygons(band_count);
    for (size_t band = 0; band < center_slices.size(); ++band)
        band_polygons[band] = to_polygons(center_slices[band]);

    // A centre cross-section supplies the interior. Adding every triangle's
    // XY projection supplies the swept boundary of the entire slab, including
    // sloped or curved portions between sample planes.
    for (const Vec3i32 &indices : mesh.its.indices) {
        Polygon projection;
        projection.points.reserve(3);
        double triangle_min_z = std::numeric_limits<double>::max();
        double triangle_max_z = std::numeric_limits<double>::lowest();
        for (int vertex = 0; vertex < 3; ++vertex) {
            const Vec3f &point = mesh.its.vertices[size_t(indices[vertex])];
            projection.points.emplace_back(scale_(point.x()), scale_(point.y()));
            triangle_min_z = std::min(triangle_min_z, double(point.z()));
            triangle_max_z = std::max(triangle_max_z, double(point.z()));
        }
        if (std::abs(projection.area()) <= 1.)
            continue;
        if (projection.area() < 0.)
            projection.reverse();
        const size_t first_band = std::min(band_count - 1, size_t(std::max(
            0., std::floor((triangle_min_z - minimum_z) / band_height))));
        const size_t last_band = std::min(band_count - 1, size_t(std::max(
            0., std::floor((triangle_max_z - minimum_z) / band_height))));
        for (size_t band = first_band; band <= last_band; ++band)
            band_polygons[band].push_back(projection);
    }

    geometry.slices.reserve(band_count);
    for (size_t band = 0; band < band_count; ++band) {
        ExPolygons merged = union_ex(band_polygons[band]);
        Polygons polygons = to_polygons(merged);
        if (polygons.empty())
            continue;
        const double lower = double(band) * band_height;
        // Keep the first body band just above the physical nozzle slice. It
        // conservatively represents the complete 0..band_height slab.
        const double slice_height = band == 0 ? std::min(0.1, 0.1 * band_height) : lower;
        geometry.slices.push_back({slice_height, false, std::move(polygons)});
    }
    normalize_geometry(geometry);
    // The collision transforms assume the nozzle axis is XY origin. Reject a
    // misplaced model instead of silently producing plausible-looking but
    // unsafe toolpaths. The lowest slab must surround that axis because its
    // lowest Z is defined as the nozzle tip.
    const Point nozzle_axis(0, 0);
    if (!geometry.slices.empty() && std::none_of(
            geometry.slices.front().polygons.begin(),
            geometry.slices.front().polygons.end(),
            [&nozzle_axis](const Polygon &polygon) {
                return polygon.contains(nozzle_axis);
            })) {
        BOOST_LOG_TRIVIAL(warning) << "Gantry model '" << path.string()
                                   << "' is not nozzle-centred at XY origin";
        geometry.slices.clear();
    }
    if (!geometry.validation_error().empty())
        geometry.slices.clear();
    {
        std::lock_guard<std::mutex> lock(cache_mutex);
        cache[path.string()] = {stamp, band_height, geometry};
    }
    return geometry;
}

SequentialGantryGeometry load_sequential_gantry_geometry(const ConfigBase &config)
{
    try {
        if (config.has("sequential_print_gantry_geometry")) {
            const auto *custom_option = config.option<ConfigOptionString>("sequential_print_gantry_geometry");
            const std::string custom = custom_option == nullptr ? std::string() : custom_option->value;
            if (!custom.empty()) {
                SequentialGantryGeometry geometry = load_custom_geometry(custom);
                if (config.has("sequential_print_gantry_model")) {
                    const auto *model_option = config.option<ConfigOptionString>("sequential_print_gantry_model");
                    const std::string model = model_option == nullptr ? std::string() : model_option->value;
                    if (!model.empty())
                        geometry.model_path = model;
                }
                if (geometry.validation_error().empty())
                    return geometry;
                BOOST_LOG_TRIVIAL(warning) << "Ignoring invalid custom gantry geometry: "
                                           << geometry.validation_error();
                return {};
            }
        }

        // An explicitly supplied STL is a user override. Profiles which leave
        // it empty continue to use the curated bundled envelope below.
        if (config.has("sequential_print_gantry_model")) {
            const auto *model_option = config.option<ConfigOptionString>("sequential_print_gantry_model");
            const std::string model = model_option == nullptr ? std::string() : model_option->value;
            if (!model.empty()) {
                if (model.compare(0, 8, "builtin:") == 0) {
                    const fs::path directory = fs::path(resources_dir()) / "data" / "printer_gantries";
                    boost::nowide::ifstream input((directory / "geometries.json").string());
                    pt::ptree root;
                    pt::read_json(input, root);
                    for (const auto &entry : root.get_child("printers"))
                        if (entry.second.get<std::string>("gantry_model_filename", "") == model.substr(8)) {
                            auto geometry = parse_geometry(entry.second, directory);
                            return geometry.validation_error().empty() ? geometry : SequentialGantryGeometry{};
                        }
                    return {};
                }
                SequentialGantryGeometry geometry = process_sequential_gantry_model(model);
                if (geometry.validation_error().empty())
                    return geometry;
                BOOST_LOG_TRIVIAL(warning) << "Unable to process custom gantry model '"
                                           << model << "': " << geometry.validation_error();
                return {};
            }
        }

        const fs::path directory = fs::path(resources_dir()) / "data" / "printer_gantries";
        boost::nowide::ifstream input((directory / "geometries.json").string());
        if (input) {
            pt::ptree root;
            pt::read_json(input, root);
            const auto *notes_option = config.option<ConfigOptionString>("printer_notes");
            const std::string notes = notes_option == nullptr ? std::string() : notes_option->value;
            if (const auto printers = root.get_child_optional("printers")) {
                for (const auto &printer : *printers) {
                    const boost::regex matcher(printer.second.get<std::string>("printer_notes_regex"));
                    if (boost::regex_match(notes, matcher)) {
                        SequentialGantryGeometry geometry = parse_geometry(printer.second, directory);
                        if (geometry.validation_error().empty())
                            return geometry;
                        BOOST_LOG_TRIVIAL(warning) << "Ignoring invalid bundled gantry geometry: "
                                                   << geometry.validation_error();
                        return {};
                    }
                }
            }
        }

    } catch (const std::exception &error) {
        // Invalid optional geometry must not make an otherwise valid machine
        // profile unusable; callers retain the legacy clearance fallback.
        BOOST_LOG_TRIVIAL(warning) << "Unable to load gantry collision geometry: "
                                   << error.what();
    }
    return {};
}

ToolheadClearanceFallback make_nonplanar_clearance_cone(const ConfigBase &config)
{
    ToolheadClearanceFallback fallback;
    if (const auto *diameters = config.option<ConfigOptionFloats>("nozzle_diameter");
        diameters != nullptr && !diameters->values.empty())
        fallback.radius = 0.5 * *std::max_element(diameters->values.begin(), diameters->values.end());
    if (const auto *height = config.option<ConfigOptionFloat>("extruder_clearance_height_to_rod"))
        fallback.height = std::max(0., height->value);
    if (const auto *angle = config.option<ConfigOptionFloat>("nonplanar_toolhead_clearance_angle"))
        fallback.angle_degrees = std::clamp(angle->value, 0., 89.);
    return fallback;
}

namespace {

// Returns the portion of the move for which a horizontal toolhead slice at
// height_offset is strictly below the obstacle top. Touching the deposition
// plane is expected and is therefore not considered a collision.
bool collision_interval(double start_z, double end_z, double height_offset,
                        double obstacle_top_z, double epsilon,
                        double &from, double &to)
{
    from = 0.;
    to   = 1.;
    const double limit_z = obstacle_top_z - epsilon - height_offset;
    const double dz = end_z - start_z;
    if (std::abs(dz) < EPSILON)
        return start_z < limit_z;

    const double crossing = (limit_z - start_z) / dz;
    if (dz > 0.)
        to = std::min(1., crossing);
    else
        from = std::max(0., crossing);
    return from < to && to > 0. && from < 1.;
}

Point segment_point(const Vec3d &start, const Vec3d &end, double t)
{
    return Point::new_scale(start.x() + (end.x() - start.x()) * t,
                            start.y() + (end.y() - start.y()) * t);
}

Polygon swept_polygon(const Polygon &shape, const Point &from, const Point &to)
{
    Points points;
    points.reserve(2 * shape.points.size());
    for (const Point &point : shape.points) {
        points.emplace_back(point + from);
        points.emplace_back(point + to);
    }
    return Geometry::convex_hull(std::move(points));
}

bool overlaps_obstacle(Polygon swept, const PrintedToolpathObstacle &obstacle,
                       double xy_clearance)
{
    BoundingBox swept_bbox = get_extents(swept.points);
    if (xy_clearance > 0.)
        swept_bbox.offset(scale_(xy_clearance));
    if (obstacle.bbox.defined && !swept_bbox.overlap(obstacle.bbox))
        return false;

    Polygons envelope;
    if (xy_clearance > 0.)
        envelope = offset(swept, float(scale_(xy_clearance)), jtRound);
    else
        envelope.emplace_back(std::move(swept));
    return !intersection(envelope, obstacle.footprint()).empty();
}

bool swept_polygon_overlaps_centerline(const Polygon &swept,
                                       const PrintedToolpathObstacle &obstacle,
                                       double xy_clearance)
{
    if (obstacle.centerline.points.size() < 2 || obstacle.line_width <= 0.)
        return false;

    const coord_t clearance = scale_(0.5 * obstacle.line_width + xy_clearance);
    const double clearance_squared = double(clearance) * double(clearance);
    BoundingBox swept_bbox = get_extents(swept.points);
    swept_bbox.offset(clearance);
    const Lines boundary = swept.lines();

    for (const Line &printed : obstacle.centerline.lines()) {
        const BoundingBox printed_bbox(Points{printed.a, printed.b});
        if (!swept_bbox.overlap(printed_bbox))
            continue;
        if (swept.contains(printed.a) || swept.contains(printed.b))
            return true;
        for (const Line &edge : boundary) {
            Point intersection_point;
            if (edge.intersection(printed, &intersection_point) ||
                std::min({edge.distance_to_squared(printed.a),
                          edge.distance_to_squared(printed.b),
                          printed.distance_to_squared(edge.a),
                          printed.distance_to_squared(edge.b)}) <= clearance_squared)
                return true;
        }
    }
    return false;
}

bool swept_bbox_overlaps(const Polygon &shape, const Point &from, const Point &to,
                         const PrintedToolpathObstacle &obstacle, double xy_clearance)
{
    if (!obstacle.bbox.defined)
        return true;
    BoundingBox shape_bbox = get_extents(shape.points);
    BoundingBox swept_bbox;
    swept_bbox.merge(shape_bbox.min + from);
    swept_bbox.merge(shape_bbox.max + from);
    swept_bbox.merge(shape_bbox.min + to);
    swept_bbox.merge(shape_bbox.max + to);
    if (xy_clearance > 0.)
        swept_bbox.offset(scale_(xy_clearance));
    return swept_bbox.overlap(obstacle.bbox);
}

template<class Obstacles> const PrintedToolpathObstacle &obstacle_at(
    const Obstacles &obstacles, size_t index)
{
    if constexpr (std::is_pointer_v<typename Obstacles::value_type>)
        return *obstacles[index];
    else
        return obstacles[index];
}

template<class Obstacles>
ToolheadCollision check_toolhead_clearance_impl(
    const Vec3d                              &start,
    const Vec3d                              &end,
    const SequentialGantryGeometry           &gantry,
    const Obstacles                          &obstacles,
    const ToolheadClearanceFallback          &fallback,
    double                                    vertical_epsilon,
    double                                    xy_clearance)
{
    if (!gantry.empty()) {
        for (size_t obstacle_idx = 0; obstacle_idx < obstacles.size(); ++obstacle_idx) {
            const PrintedToolpathObstacle &obstacle = obstacle_at(obstacles, obstacle_idx);
            if (obstacle.empty())
                continue;
            for (size_t slice_idx = 0; slice_idx < gantry.slices.size(); ++slice_idx) {
                const SequentialGantrySlice &slice = gantry.slices[slice_idx];
                double from = 0., to = 0.;
                if (!collision_interval(start.z(), end.z(), slice.height,
                                        obstacle.top_z, vertical_epsilon, from, to))
                    continue;
                const Point from_xy = segment_point(start, end, std::max(0., from));
                const Point to_xy   = segment_point(start, end, std::min(1., to));
                for (const Polygon &polygon : slice.polygons)
                    if (swept_bbox_overlaps(polygon, from_xy, to_xy, obstacle,
                                            xy_clearance)) {
                        const Polygon swept = swept_polygon(polygon, from_xy, to_xy);
                        const bool overlaps = obstacle.centerline.points.size() >= 2 && obstacle.line_width > 0. ?
                            swept_polygon_overlaps_centerline(swept, obstacle, xy_clearance) :
                            overlaps_obstacle(swept, obstacle, xy_clearance);
                        if (overlaps)
                            return {ToolheadCollisionType::GantryGeometry, obstacle_idx, slice_idx};
                    }
            }
        }
        return {};
    }

    if (fallback.radius <= 0. || fallback.height <= 0.)
        return {};

    const double angle = std::clamp(fallback.angle_degrees, 0., 89.) * PI / 180.;
    for (size_t obstacle_idx = 0; obstacle_idx < obstacles.size(); ++obstacle_idx) {
        const PrintedToolpathObstacle &obstacle = obstacle_at(obstacles, obstacle_idx);
        if (obstacle.empty())
            continue;
        double from = 0., to = 0.;
        if (!collision_interval(start.z(), end.z(), 0., obstacle.top_z,
                                vertical_epsilon, from, to))
            continue;

        from = std::max(0., from);
        to   = std::min(1., to);
        const double z_from = start.z() + (end.z() - start.z()) * from;
        const double z_to   = start.z() + (end.z() - start.z()) * to;
        const double relative_height = std::clamp(
            obstacle.top_z - std::min(z_from, z_to), 0., fallback.height);
        const double radius = fallback.radius + relative_height * std::tan(angle) + xy_clearance;
        BoundingBox centerline_bbox(Points{segment_point(start, end, from),
                                           segment_point(start, end, to)});
        centerline_bbox.offset(scale_(radius));
        if (obstacle.bbox.defined && !centerline_bbox.overlap(obstacle.bbox))
            continue;
        Polyline centerline({segment_point(start, end, from), segment_point(start, end, to)});
        if (obstacle.centerline.points.size() >= 2 && obstacle.line_width > 0.) {
            const Line move(centerline.points.front(), centerline.points.back());
            const double collision_distance = scale_(radius + 0.5 * obstacle.line_width);
            const double collision_distance_squared = collision_distance * collision_distance;
            for (const Line &printed : obstacle.centerline.lines()) {
                Point intersection_point;
                if (move.intersection(printed, &intersection_point) ||
                    std::min({move.distance_to_squared(printed.a),
                              move.distance_to_squared(printed.b),
                              printed.distance_to_squared(move.a),
                              printed.distance_to_squared(move.b)}) < collision_distance_squared)
                    return {ToolheadCollisionType::FallbackEnvelope, obstacle_idx, size_t(-1)};
            }
            continue;
        }
        const Polygons envelope = offset(centerline, float(scale_(radius)), jtRound,
                                         DefaultLineMiterLimit, ClipperLib::etOpenRound);
        if (!intersection(envelope, obstacle.footprint()).empty())
            return {ToolheadCollisionType::FallbackEnvelope, obstacle_idx, size_t(-1)};
    }
    return {};
}

template<class Obstacles>
ToolheadCollision check_toolpath_clearance_impl(
    const std::vector<Vec3d>                  &points,
    const SequentialGantryGeometry            &gantry,
    const Obstacles                            &obstacles,
    const ToolheadClearanceFallback           &fallback,
    double                                     vertical_epsilon,
    double                                     xy_clearance)
{
    if (points.size() < 2)
        return {};
    for (size_t idx = 1; idx < points.size(); ++idx) {
        ToolheadCollision collision = check_toolhead_clearance_impl(
            points[idx - 1], points[idx], gantry, obstacles, fallback,
            vertical_epsilon, xy_clearance);
        if (collision)
            return collision;
    }
    return {};
}

} // namespace

ToolheadCollision check_toolhead_clearance(
    const Vec3d                              &start,
    const Vec3d                              &end,
    const SequentialGantryGeometry           &gantry,
    const std::vector<PrintedToolpathObstacle> &obstacles,
    const ToolheadClearanceFallback          &fallback,
    double                                    vertical_epsilon,
    double                                    xy_clearance)
{
    return check_toolhead_clearance_impl(start, end, gantry, obstacles, fallback,
                                         vertical_epsilon, xy_clearance);
}

ToolheadCollision check_toolhead_clearance(
    const Vec3d                              &start,
    const Vec3d                              &end,
    const SequentialGantryGeometry           &gantry,
    const std::vector<const PrintedToolpathObstacle *> &obstacles,
    const ToolheadClearanceFallback          &fallback,
    double                                    vertical_epsilon,
    double                                    xy_clearance)
{
    return check_toolhead_clearance_impl(start, end, gantry, obstacles, fallback,
                                         vertical_epsilon, xy_clearance);
}

ToolheadCollision check_toolpath_clearance(
    const std::vector<Vec3d>                  &points,
    const SequentialGantryGeometry            &gantry,
    const std::vector<PrintedToolpathObstacle> &obstacles,
    const ToolheadClearanceFallback           &fallback,
    double                                     vertical_epsilon,
    double                                     xy_clearance)
{
    return check_toolpath_clearance_impl(points, gantry, obstacles, fallback,
                                         vertical_epsilon, xy_clearance);
}

ToolheadCollision check_toolpath_clearance(
    const std::vector<Vec3d>              &points,
    const SequentialGantryGeometry        &gantry,
    const std::vector<const PrintedToolpathObstacle *> &obstacles,
    const ToolheadClearanceFallback       &fallback,
    double                                 vertical_epsilon,
    double                                 xy_clearance)
{
    return check_toolpath_clearance_impl(points, gantry, obstacles, fallback,
                                         vertical_epsilon, xy_clearance);
}

PrintedToolpathObstacle make_printed_toolpath_obstacle(
    const std::vector<Vec3d> &points,
    double                    line_width)
{
    PrintedToolpathObstacle obstacle;
    if (points.size() < 2 || line_width <= 0.)
        return obstacle;

    obstacle.centerline.points.reserve(points.size());
    obstacle.top_z = -std::numeric_limits<double>::max();
    for (const Vec3d &point : points) {
        obstacle.centerline.points.emplace_back(Point::new_scale(point.x(), point.y()));
        obstacle.top_z = std::max(obstacle.top_z, point.z());
    }
    obstacle.line_width = line_width;
    obstacle.bbox = get_extents(obstacle.centerline.points);
    obstacle.bbox.offset(scale_(0.5 * line_width));
    return obstacle;
}

const Polygons &PrintedToolpathObstacle::footprint() const
{
    if (polygons.empty() && centerline.points.size() >= 2 && line_width > 0.)
        polygons = offset_polyline_raw(centerline, float(scale_(0.5 * line_width)), jtRound,
                                       DefaultLineMiterLimit, ClipperLib::etOpenRound);
    return polygons;
}

} // namespace Slic3r
