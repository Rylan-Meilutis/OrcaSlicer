#pragma once

#include "Point.hpp"
#include "ThermoelasticLayers.hpp"
#include <cstddef>
#include <functional>
#include <vector>

namespace Slic3r {
class Layer;
namespace Thermoelastic {

// Intended ownership from a completed native planar layer, not from CAD or
// imported G-code. Caller must hold the matching slicing generation immutable.
// The explicit placement must use the same coordinates as the deposition;
// no plate origin, nozzle offset or physical-tool remap is inferred here.
// filament_to_material maps zero-based logical filament IDs to caller-defined
// mechanical material IDs. Role-dependent assignments must all map identically
// within each region. Empty/unsupported/missing mappings reject, never default.
// Only conflicts may consume these masks: they must NOT become filled geometry.
SectionOwnership native_layer_ownership(const Layer &layer, double section_z_m,
    const Point &layer_to_deposition_xy, double layer_to_deposition_z_m,
    const std::vector<std::size_t> &filament_to_material,
    const LayerGeometryLimits &limits = {}, const std::function<void()> &checkpoint = {});

} // namespace Thermoelastic
} // namespace Slic3r
