#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

#include "ExPolygon.hpp"

namespace Slic3r {

namespace sla { class IndexedMesh; }
class Layer;
class SupportLayer;

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
