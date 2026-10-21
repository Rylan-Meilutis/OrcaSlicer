#pragma once

#include "Model.hpp"
#include "Config.hpp"
#include "Point.hpp"
#include "PrintConfig.hpp"
#include "TriangleMeshSlicer.hpp"
#include "ClipperUtils.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <set>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace Slic3r {

// Material data and applicability are deliberately separate from risk. A material
// name alone must never cause a warning or modify a project.
struct WarpRisk {
    std::set<std::string> materials;
    std::set<std::string> unknown_materials;
    std::set<std::string> resin_screened_materials;
    std::set<std::string> family_screened_materials;
    std::set<ObjectID> elevated_objects;
    double span_mm = 0.;
    double free_contraction_mm = 0.;
    double free_bow_mm = 0.;
    double restrained_energy_j_m2 = 0.;
    double crystallization_contraction_mm = 0.;
    double layer_height_mm = 0.2;
    bool estimated_ambient = false;
    bool geometry_unknown = false;
    // Compare deformation with the actual deposition resolution, not a material
    // whitelist or fixed object-size threshold. This is a potential, not a
    // prediction that the plate bond will fail.
    bool elevated() const { return free_bow_mm > layer_height_mm; }
    // Mitigation policy, NOT a physically predicted required adhesive width.
    double brim_mm() const { return elevated() ? std::clamp(2. * free_bow_mm, 3., 10.) : 0.; }
};

struct WarpMaterial {
    double cte_min, cte_max; // 1/K, representative solid-state CTE
    double tg_min, tg_max;  // degrees C
    double modulus_max;    // Pa, room-temperature upper-envelope assumption
    double melting_max = 0.; // >0 for semicrystalline resins; never substitute Tg for Tm
    double amorphous_density = 0., crystalline_density = 0.; // solid phase endpoints, g/cm3
    bool resin_screening = false; // filled grade: matrix scenario, NOT composite properties
    bool family_screening = false; // generic nylon: envelope of the documented nylon grades

    double phase_contraction_bound() const {
        return amorphous_density > 0. && crystalline_density > amorphous_density ?
            1. - std::cbrt(amorphous_density / crystalline_density) : 0.;
    }
};

inline std::optional<WarpMaterial> warp_resin_material(std::string_view type)
{
    // Published solid-state references; sources and applicability in warp-mitigation.md.
    if (type == "ABS") return WarpMaterial{80e-6, 100e-6, 95., 105., 3.0e9};
    if (type == "ASA") return WarpMaterial{85e-6, 105e-6, 95., 105., 2.9e9};
    if (type == "PC")  return WarpMaterial{75e-6, 80e-6, 140., 150., 2.4e9};
    // Representative unfilled grades, not guaranteed bounds for every blend.
    // HIPS combines published resin CTE/modulus with a measured filament Tg.
    if (type == "HIPS") return WarpMaterial{80e-6, 90e-6, 101.1, 101.1, 2.0e9};
    if (type == "PETG") return WarpMaterial{51e-6, 51e-6, 80., 80., 2.1e9};
    if (type == "PMMA") return WarpMaterial{90e-6, 110e-6, 45., 115., 3.3e9};
    // Crystalline solid may carry stress above Tg. Tm is a deliberately limiting
    // no-relaxation reference, not an assumed measured crystallization temperature.
    if (type == "PLA") return WarpMaterial{80e-6, 85e-6, 45., 65., 2.8e9, 160., 1.248, 1.290};
    if (type == "PP") return WarpMaterial{130e-6, 180e-6, -20., 20., 1.8e9, 165., 0.855, 0.936};
    if (type == "PA6" || type == "NYLON6" || type == "NYLON-6")
        return WarpMaterial{80e-6, 90e-6, 45., 80., 2.8e9, 235., 1.08, 1.23};
    if (type == "PA11" || type == "NYLON11" || type == "NYLON-11")
        return WarpMaterial{85e-6, 120e-6, 40., 55., 1.4e9, 190., 1.01, 1.18};
    if (type == "PA12" || type == "NYLON12" || type == "NYLON-12")
        return WarpMaterial{120e-6, 140e-6, 40., 50., 1.4e9, 180., 0.99, 1.085};
    if (type == "PA" || type == "NYLON")
        return WarpMaterial{80e-6, 140e-6, 40., 80., 2.8e9, 235., 1.01, 1.18, false, true};
    if (type == "POM" || type == "ACETAL")
        return WarpMaterial{160e-6, 180e-6, -85., -75., 3.2e9, 190., 1.215, 1.492};
    if (type == "PPS") return WarpMaterial{50e-6, 70e-6, 85., 100., 3.7e9, 290., 1.32, 1.43};
    if (type == "PEEK") return WarpMaterial{50e-6, 70e-6, 145., 155., 3.7e9, 345., 1.265, 1.400};
    return std::nullopt;
}

inline std::optional<WarpMaterial> warp_material(std::string_view type)
{
    // Normalize only spelling, not arbitrary blends or unknown additive grades.
    std::string name;
    for (char c : type) {
        if (c == ' ' || c == '\t') continue;
        name += c >= 'a' && c <= 'z' ? char(c - 'a' + 'A') : (c == '_' ? '-' : c);
    }
    bool filled = false;
    for (std::string_view suffix : {std::string_view("-CF"), std::string_view("-GF"),
                                   std::string_view("CF"), std::string_view("GF")}) {
        if (name.size() > suffix.size() && name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0) {
            name.resize(name.size() - suffix.size());
            filled = true;
            break;
        }
    }
    auto material = warp_resin_material(name);
    if (material) material->resin_screening = filled;
    return material;
}

// Reduced-order thermoelastic strip with a linear through-thickness gradient.
// Uniform cooling gives epsilon=alpha*dT; gradient gives curvature=alpha*dT/h.
// Free edge bow relative to the centre is curvature*L^2/8 (small deflection).
// Real restraint/relaxation can reduce all these values to zero. Energy is stored
// elastic energy per bonded area, NOT adhesion strength or fracture energy.
inline void assess_thermal_strip(WarpRisk &risk, const WarpMaterial &material,
                                 double span, double thickness, double bed, double ambient,
                                 double deposition)
{
    if (!(span > 0. && thickness > 0.) || !std::isfinite(span + thickness + bed + ambient + deposition))
        return;
    const double solidification = std::min(material.melting_max > 0. ? material.melting_max : material.tg_max, deposition);
    const double strain = material.cte_max * std::max(0., solidification - ambient);
    const double gradient = std::max(0., std::min(bed, solidification) - ambient);
    risk.span_mm = std::max(risk.span_mm, span);
    risk.free_contraction_mm = std::max(risk.free_contraction_mm, strain * span);
    risk.free_bow_mm = std::max(risk.free_bow_mm,
        material.cte_max * gradient * span * span / (8. * thickness));
    if (!material.resin_screening)
        risk.restrained_energy_j_m2 = std::max(risk.restrained_energy_j_m2,
            0.5 * material.modulus_max * strain * strain * thickness * 0.001);
    // Isotropic full phase-change bound. Keep separate from thermal contraction:
    // no cooling history/kinetics is known, and uniform phase shrinkage is NOT bow.
    // Never infer crystallinity from the spool density (fillers/voids confound it).
    if (material.melting_max > ambient && deposition > ambient)
        risk.crystallization_contraction_mm = std::max(risk.crystallization_contraction_mm,
            material.phase_contraction_bound() * span);
}

inline WarpRisk assess_warp_risk(const DynamicPrintConfig &config, const Model &model)
{
    WarpRisk risk;
    const auto *types = config.option<ConfigOptionStrings>("filament_type");
    if (!types) return risk;
    for (const ModelObject *object : model.objects) {
        if (!object->printable || std::none_of(object->instances.begin(), object->instances.end(),
                [](const ModelInstance *instance) { return instance->is_printable(); }))
            continue;
        std::set<size_t> used_slots;
        const auto inspect_slot = [&](int slot) {
            if (slot > 0 && size_t(slot) <= types->values.size())
                used_slots.insert(size_t(slot - 1));
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
        DynamicPrintConfig effective = DynamicPrintConfig::full_print_config();
        effective.apply(config);
        effective.apply(object->config.get());
        const double layer = effective.opt_float("layer_height");
        if (!(layer > 0.)) { risk.geometry_unknown = true; continue; }
        // The result stores the greatest ratio independently of other objects'
        // layer heights below, so a fine unrelated object cannot trigger a brim.
        const auto temperature = [&](const std::string &key, size_t slot, double fallback) {
            const auto *values = effective.option<ConfigOptionInts>(key);
            return values && !values->values.empty() ? double(values->get_at(slot)) : fallback;
        };
        const auto *plate = effective.option<ConfigOptionEnum<BedType>>("curr_bed_type");
        std::vector<std::pair<double, double>> sections; // span and solid base thickness, mm
        if (std::any_of(used_slots.begin(), used_slots.end(), [&](size_t slot) {
                return warp_material(types->values[slot]).has_value(); })) {
            for (const ModelInstance *instance : object->instances) {
                if (!instance->is_printable()) continue;
                TriangleMesh mesh = object->raw_mesh();
                instance->transform_mesh(&mesh);
                const auto bounds = mesh.bounding_box();
                const double height = bounds.size().z();
                const double shell = std::max(effective.opt_float("bottom_shell_thickness"),
                    effective.opt_int("bottom_shell_layers") * layer);
                if (height <= 0. || shell <= 0.) { risk.geometry_unknown = true; continue; }
                const auto contacts = union_ex(slice_mesh(mesh.its,
                    float(bounds.min.z() + std::min(layer, height) * 0.5), MeshSlicingParams{}));
                if (contacts.empty()) risk.geometry_unknown = true;
                for (const auto &contact : contacts) {
                    const Point centre = contact.contour.centroid();
                    double radius = 0.;
                    for (const Point &point : contact.contour.points)
                        radius = std::max(radius, unscale<double>((point - centre).cast<double>().norm()));
                    sections.emplace_back(2. * radius, std::min(height, shell));
                }
            }
        }
        for (size_t slot : used_slots) {
            const auto material = warp_material(types->values[slot]);
            if (!material) { risk.unknown_materials.insert(types->values[slot]); continue; }
            risk.materials.insert(types->values[slot]);
            if (material->resin_screening) risk.resin_screened_materials.insert(types->values[slot]);
            if (material->family_screening) risk.family_screened_materials.insert(types->values[slot]);
            double ambient = temperature("chamber_temperature", slot, 0.);
            if (ambient <= 0.) { ambient = 20.; risk.estimated_ambient = true; }
            const double bed = plate ? temperature(get_bed_temp_key(plate->value), slot, -1.) : -1.;
            const double deposition = temperature("nozzle_temperature", slot, -1.);
            if (bed < 0. || deposition <= 0.) { risk.geometry_unknown = true; continue; }
            for (const auto &[span, thickness] : sections) {
                WarpRisk section;
                assess_thermal_strip(section, *material, span, thickness, bed, ambient, deposition);
                if (section.free_bow_mm > layer)
                    risk.elevated_objects.insert(object->id());
                risk.span_mm = std::max(risk.span_mm, section.span_mm);
                risk.free_contraction_mm = std::max(risk.free_contraction_mm, section.free_contraction_mm);
                risk.restrained_energy_j_m2 = std::max(risk.restrained_energy_j_m2, section.restrained_energy_j_m2);
                risk.crystallization_contraction_mm = std::max(risk.crystallization_contraction_mm, section.crystallization_contraction_mm);
                if (section.free_bow_mm / layer > risk.free_bow_mm / risk.layer_height_mm) {
                    risk.free_bow_mm = section.free_bow_mm;
                    risk.layer_height_mm = layer;
                }
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
