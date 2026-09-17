#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <vector>

#include "ExPolygon.hpp"

namespace Slic3r {

namespace sla { class IndexedMesh; }
class Layer;
class LayerRegion;
class SupportLayer;
class ExtrusionPath;
class ExtrusionEntity;

// Filter an aligned course without permitting removal of a dimensional wall.
// On rejection the caller retains ownership of the unchanged complete course.
bool filter_nonplanar_course_paths(std::vector<ExtrusionEntity *> &paths,
                                   const std::vector<uint8_t> &retained);

// Arachne wall fragments have not been joined yet; only discard short fill.
void remove_short_nonplanar_fill_paths(std::vector<ExtrusionEntity *> &paths);

// Native centerlines at a wall inset, retaining parent loop/multipath ownership.
Polylines nonplanar_wall_centerlines(const ExtrusionEntity &entity, int inset);

// Caller owns the returned projected solid-fill entities. Layer heights are
// in print space; the selected mesh may have a different ground level.
std::vector<ExtrusionEntity *> project_nonplanar_source_course(
    LayerRegion &source_region, const sla::IndexedMesh &mesh,
    const ExPolygons &surface_projection, const std::vector<uint8_t> &facets,
    const ExPolygons &course_domain, coordf_t destination_owner_z,
    double maximum_drape_height);

struct NonplanarCourseProfile {
    double plane_z;
    double blend;
};

// Index of the first solid course required by layer count and physical shell
// thickness. Zero uses the whole available stack; any additional required
// depth must be provided by the conventional material below its anchor.
size_t nonplanar_solid_course_begin(
    const std::vector<NonplanarCourseProfile> &courses,
    double minimum_surface_z, double maximum_surface_z,
    size_t minimum_layers, double minimum_thickness);

// Transfer the existing in-path entry/exit blend onto a planned course.
// Mixed walls keep their planar join; wholly selected paths have no such join.
double nonplanar_wall_course_z(double planar_z, double projected_z,
    double surface_z, double course_z, bool fully_selected_path);

// A flat foundation followed by a monotone, smooth adoption of the surface
// slope. Empty means this anchor cannot satisfy the supplied bead heights.
std::vector<NonplanarCourseProfile> plan_nonplanar_courses(
    double anchor_z, double minimum_surface_z, double maximum_surface_z,
    size_t minimum_courses, double minimum_height, double maximum_height);

// Sample the selected surface, extending its tangent across enclosed feature
// openings. Outside the patch's outer boundary, retain the boundary height.
std::optional<double> project_nonplanar_patch_height(
    const sla::IndexedMesh &mesh, const ExPolygons &projection,
    const std::vector<uint8_t> &facets, const Point &point);

// Preserve the complete XYZ chain and wall metadata while assigning one
// validated local bead height to each segment. Adjacent equal heights merge.
std::vector<ExtrusionPath> split_nonplanar_path_by_height(
    const ExtrusionPath &path, const std::vector<double> &segment_heights);

// Parameter intervals on [a,b] whose XY distance to [c,d] and vertical
// separation are both within the supplied limits. Coordinates are physical
// millimetres, including absolute Z; no nominal layer plane is assumed.
std::vector<std::pair<double, double>> nonplanar_segment_overlap(
    const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &d,
    double xy_tolerance, double z_tolerance);

// Validate wall contact with deposited paths in physical Z. Fill/bridge
// support is planned separately; this check does not certify those roles.
bool nonplanar_walls_are_supported(
    const ExtrusionEntity &entity,
    const std::vector<const ExtrusionPath *> &support_paths,
    double entity_owner_z, double support_owner_z,
    bool allow_short_unsupported_runs = false);

// Require deposited contact for every path, including fill. This does not
// classify unsupported fill as a bridge; bridge planning belongs upstream.
bool nonplanar_course_is_supported(
    const ExtrusionEntity &entity,
    const std::vector<const ExtrusionPath *> &support_paths,
    double entity_owner_z, double support_owner_z);

bool native_foundation_supports_course(
    const std::vector<ExtrusionEntity *> &course, const LayerRegion &foundation,
    double owner_z, double foundation_z, double minimum_height, double maximum_height);

enum class NonplanarProgressStage : uint8_t {
    SurfaceTopology,
    SurfaceProjection,
    BoundaryRepair,
    ShellClearance,
    FeatureScheduling,
    FeatureProjection,
    FeatureClearance,
    Commit
};

struct NonplanarProgress {
    double progress {0.};
    NonplanarProgressStage stage {NonplanarProgressStage::SurfaceTopology};
    size_t current {0};
    size_t total {0};
};

// One entry per mesh face. A face is selected only when it belongs to a
// connected, upward-facing, genuinely sloped surface patch large enough to
// print as a coherent skin.
std::vector<uint8_t> select_nonplanar_surface_facets(
    const sla::IndexedMesh &mesh, double maximum_angle_degrees,
    double minimum_vertical_span, double minimum_projected_area);

// Build a cheap XY query domain for each region and reachable layer Z band.
// Candidate marking uses it to avoid mesh raycasts for ordinary paths which
// are nowhere near a selected non-planar surface at that layer.
std::vector<std::vector<ExPolygons>> project_nonplanar_surface_facets(
    const sla::IndexedMesh &mesh,
    const std::vector<std::vector<uint8_t>> &selected_facets,
    const std::vector<std::pair<double, double>> &layer_z_ranges);

// Replace the configured upper shell depth with complete, mutually supported
// draped skins. Intermediate skins are emitted on successive nominal layers;
// the last skin follows the exposed model surface continuously in XYZ.
void consolidate_nonplanar_top_surfaces(
    const std::vector<Layer *> &layers, const sla::IndexedMesh &mesh,
    const std::vector<std::vector<uint8_t>> &selected_facets,
    const std::function<void(const NonplanarProgress &)> &progress = {});

// Project a top contact interface and its ironing pass onto a reachable
// downward-facing model surface. Each path is committed atomically; an
// incomplete projection or failed support / clearance check retains the
// original planar path.
size_t project_nonplanar_support_interface(
    const std::vector<Layer *> &model_layers,
    const std::vector<SupportLayer *> &support_layers,
    const sla::IndexedMesh &mesh, double contact_gap,
    double maximum_angle_degrees, double resolution,
    const std::function<void(size_t, size_t)> &progress = {});

} // namespace Slic3r
