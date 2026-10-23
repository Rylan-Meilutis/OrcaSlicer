#pragma once

#include "Config.hpp"
#include "Preset.hpp"
#include "PrintConfig.hpp"
#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <map>
#include <set>
#include <string>
#include <utility>

namespace Slic3r {

// Explicit profile edits only: never called by Print::apply or G-code generation.
// Use resolved printer/material presets for fan percentages and layer-time limits
// instead of assuming equal airflow from equal PWM on different printers.
inline bool apply_material_cooling(DynamicPrintConfig &target, const DynamicPrintConfig &source,
                                   const DynamicPrintConfig &printer, size_t variant_index = 0)
{
    const auto *material = target.option<ConfigOptionStrings>("filament_type");
    const auto *source_material = source.option<ConfigOptionStrings>("filament_type");
    if (!material || material->values.size() != 1 || !source_material ||
        material->values != source_material->values) return false;

    DynamicPrintConfig edited = target;
    const auto *variants = target.option<ConfigOptionStrings>("filament_extruder_variant");
    const auto *source_variants = source.option<ConfigOptionStrings>("filament_extruder_variant");
    if (variant_index >= (variants && !variants->values.empty() ? variants->values.size() : 1)) return false;
    for (const char *key : {"fan_min_speed", "fan_max_speed", "reduce_fan_stop_start_freq",
             "fan_cooling_layer_time", "slow_down_layer_time", "slow_down_min_speed",
             "slow_down_for_layer_cooling", "dont_slow_down_outer_wall",
             "close_fan_the_first_x_layers", "initial_layer_fan_speed", "full_fan_speed_layer",
             "enable_overhang_bridge_fan", "overhang_fan_speed", "overhang_fan_threshold",
             "internal_bridge_fan_speed", "support_material_interface_fan_speed", "ironing_fan_speed",
             "additional_cooling_fan_speed", "close_additional_fan_first_x_layers",
             "first_x_layer_fan_speed", "additional_fan_full_speed_layer"}) {
        auto *dest = dynamic_cast<ConfigOptionVectorBase *>(edited.option(key));
        const auto *src = dynamic_cast<const ConfigOptionVectorBase *>(source.option(key));
        if (!dest || !src || src->size() == 0 || dest->size() == 0 || dest->type() != src->type()) return false;
        const bool per_variant = filament_options_with_variant.count(key) != 0;
        if (per_variant && dest->size() == 1 && variants && variants->values.size() > 1)
            dest->resize(variants->values.size());
        if (per_variant && variant_index >= dest->size()) return false;
        for (size_t i = per_variant ? variant_index : 0; i < (per_variant ? variant_index + 1 : dest->size()); ++i) {
            size_t source_index = 0;
            if (per_variant && variants && source_variants && !variants->values.empty()) {
                if (i >= variants->values.size()) return false;
                const auto found = std::find(source_variants->values.begin(), source_variants->values.end(), variants->values[i]);
                if (found == source_variants->values.end()) return false;
                source_index = size_t(found - source_variants->values.begin());
            } else if (src->size() > 1) {
                // Do not copy ambiguous multi-filament vectors into one preset.
                return false;
            }
            if (src->size() == 1) source_index = 0;
            if (source_index >= src->size()) return false;
            dest->set_at(src, i, source_index);
        }
    }
    const auto clear = [&](const char *key) {
        auto *option = dynamic_cast<ConfigOptionVectorBase *>(edited.option(key));
        if (option) {
            const size_t begin = filament_options_with_variant.count(key) ? variant_index : 0;
            const size_t end = filament_options_with_variant.count(key) ? variant_index + 1 : option->size();
            if (auto *floats = dynamic_cast<ConfigOptionFloats *>(option))
                std::fill(floats->values.begin() + begin, floats->values.begin() + end, 0.);
            else if (auto *ints = dynamic_cast<ConfigOptionInts *>(option))
                std::fill(ints->values.begin() + begin, ints->values.begin() + end, 0);
            else if (auto *bools = dynamic_cast<ConfigOptionBools *>(option))
                std::fill(bools->values.begin() + begin, bools->values.begin() + end, false);
        }
    };
    // Conservative strength/warp policy for these unfilled material types.
    // Blends, filled materials and flexibles retain their own preset's policy.
    const std::string &type = material->values.front();
    const bool heat_sensitive = type == "ASA" || type == "ABS" || type == "PC";
    if (heat_sensitive)
        for (const char *key : {"fan_min_speed", "fan_max_speed", "reduce_fan_stop_start_freq"}) clear(key);
    const auto *aux = printer.option<ConfigOptionBool>("auxiliary_fan");
    if (heat_sensitive || !aux || !aux->value) {
        clear("additional_cooling_fan_speed");
        clear("first_x_layer_fan_speed");
    }
    // Avoid an inherited explicit first-layer override bypassing adhesion cooling.
    clear("initial_layer_fan_speed");
    target = std::move(edited);
    return true;
}

inline const Preset *find_cooling_reference(const PresetCollection &filaments, const Preset &edited,
                                           const PresetWithVendorProfile &printer, size_t variant_index = 0)
{
    std::map<std::string, int> ancestors;
    int distance = 0;
    for (const Preset *p = &edited; p && ancestors.emplace(p->name, distance++).second;
         p = filaments.get_preset_parent(*p)) {}
    const Preset *best = nullptr;
    int best_score = -1;
    for (const Preset &candidate : filaments.get_presets()) {
        if (!candidate.is_system || candidate.is_default ||
            !is_compatible_with_printer(filaments.get_preset_with_vendor_profile(candidate), printer)) continue;
        const bool inherited = ancestors.count(candidate.name) != 0;
        const bool same_product = !edited.filament_id.empty() && edited.filament_id == candidate.filament_id;
        const bool generic = candidate.name.rfind("Generic ", 0) == 0;
        if (!inherited && !same_product && !generic) continue;
        DynamicPrintConfig trial = edited.config;
        if (!apply_material_cooling(trial, candidate.config, printer.preset.config, variant_index)) continue;
        const auto *compatible = candidate.config.option<ConfigOptionStrings>("compatible_printers");
        const int score = (inherited ? 1000 - 20 * ancestors.at(candidate.name) : same_product ? 80 : 40) +
                          (compatible && !compatible->values.empty() ? 10 : 0);
        if (!best || score > best_score || (score == best_score && candidate.name < best->name)) {
            best = &candidate;
            best_score = score;
        }
    }
    return best;
}

} // namespace Slic3r
