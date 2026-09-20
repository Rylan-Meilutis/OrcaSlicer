#pragma once

#include <string>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#include "BoundingBox.hpp"
#include "Polygon.hpp"
#include "Point.hpp"
#include "Polyline.hpp"

namespace Slic3r {

class ConfigBase;

struct SequentialGantrySlice
{
    double   height = 0.;
    bool     is_box = false;
    Polygons polygons;
};

struct SequentialGantryGeometry
{
    std::string                        model_path;
    std::vector<SequentialGantrySlice> slices;

    bool empty() const { return slices.empty(); }
    double conservative_clearance_radius() const;
    // Axis-wise nozzle-to-head reach. Rod/lid boxes are checked separately.
    Vec2d clearance_reach() const;
    double maximum_xy_reach(double maximum_height = std::numeric_limits<double>::max()) const;
    void set_nozzle_tip_diameter(double diameter);
    double first_box_height() const;
    // Empty means the geometry is suitable for collision planning. The
    // diagnostic is intentionally human readable so profile editors can
    // report why a newly supplied model was rejected.
    std::string validation_error() const;
};

// Expand a convex footprint by a rectangle in machine space.
Polygon sequential_clearance_hull(const Polygon &footprint, const Vec2d &half_extent);

// Converts a nozzle-centred STL into conservative horizontal collision bands.
// The model's lowest Z is treated as the nozzle tip. Each generated footprint
// encloses the complete mesh slab above its height, so using it for non-planar
// checks cannot miss geometry between sampled planes.
SequentialGantryGeometry process_sequential_gantry_model(
    const std::string &model_path, double band_height = 1.);

// Loads a custom machine definition when configured, otherwise matches the
// printer notes against the bundled Prusa-style geometry catalog.
SequentialGantryGeometry load_sequential_gantry_geometry(const ConfigBase &config);

// A conservative 2.5D representation of material which has already been
// printed. Its exact footprint is constructed lazily from the centerline, then
// treated as occupied from the bed through top_z. This is intentionally
// independent of Layer so it may also describe paths already emitted on the
// current non-planar layer.
struct PrintedToolpathObstacle
{
    double   top_z = 0.;
    mutable Polygons polygons;
    BoundingBox bbox;
    Polyline centerline;
    double   line_width = 0.;

    PrintedToolpathObstacle() = default;
    PrintedToolpathObstacle(double top_z, Polygons polygons) :
        top_z(top_z), polygons(std::move(polygons)), bbox(get_extents(this->polygons)) {}

    bool empty() const { return polygons.empty() && centerline.points.size() < 2; }
    const Polygons &footprint() const;
};

struct ToolheadClearanceFallback
{
    // Radius at the nozzle tip and maximum height of the modeled toolhead.
    double radius = 0.;
    double height = 0.;
    // Outward toolhead angle measured from vertical. Zero models a cylinder.
    double angle_degrees = 0.;
};

// Builds the rotationally symmetric non-planar fallback used when a printer
// has no detailed gantry geometry. The cone starts at the physical nozzle tip;
// the much larger sequential-print clearance radius is intentionally not used.
ToolheadClearanceFallback make_nonplanar_clearance_cone(const ConfigBase &config);

enum class ToolheadCollisionType
{
    None,
    GantryGeometry,
    FallbackEnvelope
};

struct ToolheadCollision
{
    ToolheadCollisionType type = ToolheadCollisionType::None;
    size_t obstacle_index = size_t(-1);
    size_t slice_index = size_t(-1);

    explicit operator bool() const { return type != ToolheadCollisionType::None; }
};

// Tests the complete swept toolhead envelope for a linear XYZ move. Material
// exactly level with the nozzle is allowed; material above the nozzle or a
// gantry slice is a collision. When detailed gantry slices are unavailable,
// fallback is used as a conservative conical envelope.
ToolheadCollision check_toolhead_clearance(
    const Vec3d                         &start,
    const Vec3d                         &end,
    const SequentialGantryGeometry      &gantry,
    const std::vector<PrintedToolpathObstacle> &obstacles,
    const ToolheadClearanceFallback     &fallback = {},
    double                               vertical_epsilon = 1e-4,
    double                               xy_clearance = 0.);

ToolheadCollision check_toolhead_clearance(
    const Vec3d                              &start,
    const Vec3d                              &end,
    const SequentialGantryGeometry           &gantry,
    const std::vector<const PrintedToolpathObstacle *> &obstacles,
    const ToolheadClearanceFallback          &fallback = {},
    double                                    vertical_epsilon = 1e-4,
    double                                    xy_clearance = 0.);

// Convenience helpers for validating and then recording an ordered extrusion
// path. Points are absolute machine-space positions. The recorded obstacle is
// deliberately conservative: the whole XY footprint is occupied through the
// highest Z reached by the path.
ToolheadCollision check_toolpath_clearance(
    const std::vector<Vec3d>             &points,
    const SequentialGantryGeometry       &gantry,
    const std::vector<PrintedToolpathObstacle> &obstacles,
    const ToolheadClearanceFallback      &fallback = {},
    double                                vertical_epsilon = 1e-4,
    double                                xy_clearance = 0.);

ToolheadCollision check_toolpath_clearance(
    const std::vector<Vec3d>              &points,
    const SequentialGantryGeometry        &gantry,
    const std::vector<const PrintedToolpathObstacle *> &obstacles,
    const ToolheadClearanceFallback       &fallback = {},
    double                                 vertical_epsilon = 1e-4,
    double                                 xy_clearance = 0.);

PrintedToolpathObstacle make_printed_toolpath_obstacle(
    const std::vector<Vec3d> &points,
    double                    line_width);

} // namespace Slic3r
