#pragma once

#include "Config.hpp"
#include "PrintConfig.hpp"

namespace Slic3r {

// Imported parts may explicitly bind every feature to the same material.
// A user's whole-part assignment replaces that uniform binding. Keep partial
// and mixed feature overrides, which express a separate feature-level choice.
// Call only for an explicit assignment, never when loading or dispatching a job.
inline void clear_uniform_feature_filament_overrides(ModelConfig &config)
{
    constexpr const char *keys[] = {"outer_wall_filament_id", "inner_wall_filament_id",
        "sparse_infill_filament_id", "internal_solid_filament_id",
        "top_surface_filament_id", "bottom_surface_filament_id"};
    int material = 0;
    for (const char *key : keys) {
        const auto *value = config.get().option<ConfigOptionInt>(key);
        if (!value || value->value <= 0 || (material != 0 && material != value->value))
            return;
        material = value->value;
    }
    for (const char *key : keys)
        config.erase(key);
}

} // namespace Slic3r
