#include "ThermoelasticOwnership.hpp"
#include "ClipperUtils.hpp"
#include "ExPolygon.hpp"
#include "Flow.hpp"
#include "Layer.hpp"
#include "Point.hpp"
#include "Polygon.hpp"
#include "Print.hpp"
#include "PrintConfig.hpp"
#include "Surface.hpp"
#include "SurfaceCollection.hpp"
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticMesh.hpp"
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
namespace {
long long native_z(double metres)
{
    if (!std::isfinite(metres) || std::abs(metres) > 10.)
        throw std::invalid_argument("Invalid native ownership coordinate.");
    return std::llround(metres * 1e9);
}
}

SectionOwnership native_layer_ownership(const Layer &layer, double section_z_m,
    const Point &translation, double z_translation,
    const std::vector<std::size_t> &filament_to_material,
    const LayerGeometryLimits &limits, const std::function<void()> &checkpoint)
{
    if (checkpoint) checkpoint();
    if (filament_to_material.empty() || filament_to_material.size() > 256 ||
        limits.max_boundary_points == 0 || limits.max_boundary_points > 2000000 ||
        limits.max_groups == 0 || limits.max_groups > 10000 ||
        std::abs(double(translation.x())) > 1e10 || std::abs(double(translation.y())) > 1e10)
        throw std::invalid_argument("Invalid native ownership mapping or budget.");
    const auto z = native_z(section_z_m);
    native_z(z_translation);
    const auto bottom = native_z(layer.bottom_z() * .001 + z_translation);
    const auto top = native_z(layer.print_z * .001 + z_translation);
    if (top <= bottom)
        throw std::invalid_argument("Invalid native ownership layer interval.");
    const auto *object = layer.object();
    if (!object || !object->is_step_done(posSlice) || !object->is_step_done(posPerimeters) ||
        !object->is_step_done(posInfill))
        throw std::invalid_argument("Native ownership requires completed slicing and toolpaths.");
    SectionOwnership result{section_z_m, {}};
    if (z < bottom || z >= top) return result;
    if (object->print()->config().spiral_mode.value)
        throw std::invalid_argument("Planar ownership does not describe spiral extrusion.");
    const auto &object_config = object->config();
    if (object_config.flush_into_objects.value || object_config.flush_into_infill.value)
        throw std::invalid_argument("Native ownership requires final purge extrusion assignments.");
    if (object_config.thick_bridges.value || object_config.thick_internal_bridges.value)
        throw std::invalid_argument("Planar ownership does not qualify height-spanning bridge extrusion.");
    std::map<std::size_t, ExPolygons> regions;
    std::size_t point_count = 0;
    for (const auto *layer_region : layer.regions()) {
        if (checkpoint) checkpoint();
        if (layer_region->slices.empty()) continue;
        const auto &region = layer_region->region();
        const auto &config = region.config();
        if (config.fuzzy_skin.value != FuzzySkinType::None && config.fuzzy_skin.value != FuzzySkinType::Disabled_fuzzy)
            throw std::invalid_argument("Native ownership does not qualify fuzzy-skin displacement or variable extrusion.");
        if (config.perimeter_layering.value != PerimeterLayeringMode::Standard ||
            config.top_surface_z_mode.value != TopSurfaceZMode::Disabled ||
            config.seam_slope_type.value != SeamScarfType::None || config.nonplanar_infill.value ||
            config.zaa_enabled.value || config.infill_combination.value)
            throw std::invalid_argument("Planar ownership does not describe height-modified toolpaths.");
        const auto mapped = [&](unsigned filament) {
            if (filament == 0 || filament > filament_to_material.size())
                throw std::invalid_argument("Missing logical filament ownership mapping.");
            const auto &print_config = object->print()->config();
            if ((!print_config.filament_is_mixed.values.empty() &&
                 print_config.filament_is_mixed.get_at(filament - 1)) ||
                (!print_config.filament_mixed_gradient.values.empty() &&
                 print_config.filament_mixed_gradient.get_at(filament - 1)))
                throw std::invalid_argument("Planar ownership does not describe mixed extrusion.");
            return filament_to_material[filament - 1];
        };
        const auto material = mapped(region.extruder(frPerimeter));
        for (const auto role : {frExternalPerimeter, frInfill, frSolidInfill, frTopSolidInfill})
            if (mapped(region.extruder(role)) != material)
                throw std::invalid_argument("Region ownership is ambiguous across extrusion roles.");
        if (config.bottom_surface_filament_id.value <= 0 ||
            mapped(static_cast<unsigned>(config.bottom_surface_filament_id.value)) != material)
            throw std::invalid_argument("Region ownership is ambiguous for bottom surfaces.");
        if (regions.count(material) == 0 && regions.size() >= limits.max_groups)
            throw ResourceLimit("Native ownership exceeds the material-group budget.");
        auto &polygons = regions[material];
        for (const auto &surface : layer_region->slices.surfaces) {
            if (checkpoint) checkpoint();
            const auto &source = surface.expolygon;
            const auto count = [&](const Polygon &polygon) {
                if (polygon.points.size() > limits.max_boundary_points - point_count)
                    throw ResourceLimit("Native ownership exceeds the boundary-point budget.");
                point_count += polygon.points.size();
                for (const auto &point : polygon.points)
                    if (std::abs(double(point.x()) + double(translation.x())) > 1e10 ||
                        std::abs(double(point.y()) + double(translation.y())) > 1e10)
                        throw std::invalid_argument("Native ownership exceeds the coordinate range.");
            };
            count(source.contour);
            for (const auto &hole : source.holes) count(hole);
            polygons.push_back(source);
            polygons.back().translate(translation);
        }
    }
    std::size_t output_points = 0;
    for (auto &[material, polygons] : regions) {
        if (checkpoint) checkpoint();
        polygons = union_ex(polygons);
        for (const auto &polygon : polygons) {
            output_points += polygon.contour.points.size();
            for (const auto &hole : polygon.holes) output_points += hole.points.size();
        }
        if (output_points > limits.max_boundary_points)
            throw ResourceLimit("Native ownership exceeds the output boundary-point budget.");
        result.regions.push_back({material, std::move(polygons)});
    }
    return result;
}
} // namespace Slic3r::Thermoelastic
