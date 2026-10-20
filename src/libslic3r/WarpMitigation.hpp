#pragma once

#include "Model.hpp"
#include "Config.hpp"
#include "Point.hpp"
#include "PrintConfig.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <set>
#include <string>
#include <string_view>

namespace Slic3r {

// Screening heuristics, not a physical prediction of adhesion or deformation.
// See docs/HLSD/warp-mitigation.md for limits and material references.
inline bool warp_prone_material(std::string_view type)
{
    const auto suffix = type.find('-');
    const auto matrix = type.substr(0, suffix);
    return matrix == "ABS" || matrix == "ASA" || matrix == "PC" ||
           matrix == "PA" || matrix == "PA6" || matrix == "PA12" ||
           matrix == "PA66" || matrix == "PP" || matrix == "PPS" ||
           matrix == "HIPS" || matrix == "POM" || matrix == "PEEK" ||
           matrix == "PEKK" || matrix == "PSU" || matrix == "PES" || matrix == "PESU";
}

struct WarpRisk {
    std::set<std::string> materials;
    double span_mm = 0.;
    bool elevated() const { return !materials.empty(); }
    // Bounded adhesion allowance. This is not a calibrated warpage model.
    double brim_mm() const { return span_mm >= 100. ? 10. : span_mm >= 40. ? 8. : 6.; }
};

inline WarpRisk assess_warp_risk(const DynamicPrintConfig &config, const Model &model)
{
    WarpRisk risk;
    const auto *types = config.option<ConfigOptionStrings>("filament_type");
    if (!types) return risk;
    for (const ModelObject *object : model.objects) {
        if (!object->printable || std::none_of(object->instances.begin(), object->instances.end(),
                [](const ModelInstance *instance) { return instance->is_printable(); }))
            continue;
        bool prone = false;
        const auto inspect_slot = [&](int slot) {
            if (slot > 0 && size_t(slot) <= types->values.size() &&
                warp_prone_material(types->values[slot - 1])) {
                risk.materials.insert(types->values[slot - 1]);
                prone = true;
            }
        };
        const auto inspect_features = [&](const DynamicPrintConfig &settings) {
            for (const char *key : {"outer_wall_filament_id", "inner_wall_filament_id",
                                    "sparse_infill_filament_id", "internal_solid_filament_id",
                                    "top_surface_filament_id", "bottom_surface_filament_id"}) {
                if (const auto *value = settings.option<ConfigOptionInt>(key))
                    inspect_slot(value->value);
            }
        };
        for (const ModelVolume *volume : object->volumes) {
            if (!volume->is_model_part() && !volume->is_modifier()) continue;
            const auto slots = volume->get_extruders();
            if (slots.empty()) inspect_slot(project_default_filament(model));
            for (int slot : slots) inspect_slot(slot);
            inspect_features(volume->config.get());
        }
        inspect_features(config);
        inspect_features(object->config.get());
        for (const auto &range : object->layer_config_ranges)
            inspect_features(range.second.get());
        if (prone) {
            for (size_t i = 0; i < object->instances.size(); ++i) {
                if (!object->instances[i]->is_printable()) continue;
                const Vec3d size = object->instance_bounding_box(i).size();
                risk.span_mm = std::max(risk.span_mm, std::hypot(size.x(), size.y()));
            }
        }
    }
    return risk;
}

inline void apply_warp_mitigations(DynamicPrintConfig &config, const WarpRisk &risk)
{
    if (!risk.elevated()) return;
    // Do not alter temperature, cooling, geometry, or material strength settings.
    // An opted-in continuous brim replaces sparse ears/painted coverage; the
    // painting itself stays in the model and is restored by selecting that mode.
    const auto *raft = config.option<ConfigOptionInt>("raft_layers");
    const auto *belt = config.option<ConfigOptionBool>("belt_printer");
    const auto *brim = config.option<ConfigOptionEnum<BrimType>>("brim_type");
    if ((!raft || raft->value == 0) && (!belt || !belt->value) && brim) {
        config.set_key_value("brim_type", new ConfigOptionEnum<BrimType>(
            brim->value == btInnerOnly || brim->value == btOuterAndInner ? btOuterAndInner : btOuterOnly));
        const auto *width = config.option<ConfigOptionFloat>("brim_width");
        config.set_key_value("brim_width", new ConfigOptionFloat(std::max(width ? width->value : 0., risk.brim_mm())));
        const auto *gap = config.option<ConfigOptionFloat>("brim_object_gap");
        if (gap && gap->value > 0.1)
            config.set_key_value("brim_object_gap", new ConfigOptionFloat(0.1));
    }
    for (const char *key : {"enable_overhang_speed", "slowdown_for_curled_perimeters"}) {
        if (auto *value = config.option<ConfigOptionBoolsNullable>(key))
            std::fill(value->values.begin(), value->values.end(), true);
    }
}

} // namespace Slic3r
