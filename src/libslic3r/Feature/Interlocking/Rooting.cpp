#include "Rooting.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include <algorithm>
#include <cmath>

namespace Slic3r {

Rooting::Rooting(const PrintObject &object)
{
    if (!object.config().rooting || object.num_printing_regions() < 2 || object.print()->config().spiral_mode)
        return;
    for (const Layer *layer : object.layers()) {
        auto &regions = m_slices.emplace_back();
        for (const LayerRegion *region : layer->regions())
            regions.push_back(to_expolygons(region->slices.surfaces));
    }
}

namespace {
Polygon root_disk(double x, double y, double radius)
{
    Points points;
    for (unsigned i = 0; i < 24; ++i) {
        const double angle = 2. * PI * i / 24.;
        points.emplace_back(scale_(x + radius * std::cos(angle)), scale_(y + radius * std::sin(angle)));
    }
    return Polygon(std::move(points));
}

bool coherent_material(const PrintRegion &region, const PrintConfig &config)
{
    const auto tool = region.extruder(frExternalPerimeter);
    if (config.filament_is_support.get_at(tool - 1) || config.filament_soluble.get_at(tool - 1))
        return false;
    // A root must use its owning material throughout, not a different material
    // in its core because of a feature-specific filament override.
    for (FlowRole role : {frPerimeter, frInfill, frSolidInfill, frTopSolidInfill})
        if (region.extruder(role) != tool)
            return false;
    return region.config().bottom_surface_filament_id.value == int(tool) &&
           region.config().wall_loops.value > 0 && region.config().top_shell_layers.value > 0 &&
           region.config().bottom_shell_layers.value > 0;
}
}

size_t Rooting::apply(PrintObject &object, const std::function<void()> &cancel) const
{
    if (m_slices.empty())
        return 0;
    const auto &config = object.config();
    const double depth = config.rooting_depth.value;
    const double width = config.rooting_width.value;
    const double skin = config.rooting_skin.value;
    if (!std::isfinite(depth) || !std::isfinite(width) || !std::isfinite(skin) ||
        depth < 0.4 || depth > 20. || width < 0.4 || width > 10. || skin < 0.4 || skin > 10. ||
        !std::isfinite(config.rooting_spacing.value) || config.rooting_spacing.value < 1.)
        return 0;

    size_t count = 0;
    const auto &layers = object.layers();
    for (size_t top = 1; top < layers.size(); ++top) {
        cancel();
        const double contact_z = layers[top - 1]->print_z;
        // Keep full layers within the requested depth, plus an untouched floor.
        size_t bottom = top - 1;
        while (bottom > 0 && contact_z - (layers[bottom - 1]->print_z - layers[bottom - 1]->height) <= depth + EPSILON)
            --bottom;
        if (top - bottom < 3 || bottom == 0)
            continue;
        const double bottom_z = layers[bottom]->print_z - layers[bottom]->height;
        size_t floor = bottom;
        while (floor > 0 && bottom_z - (layers[floor]->print_z - layers[floor]->height) < skin - EPSILON)
            --floor;
        if (bottom_z - (layers[floor]->print_z - layers[floor]->height) < skin - EPSILON)
            continue;
        const double taper_height = contact_z - layers[bottom]->print_z;
        // At most 0.5 mm lateral change per mm of height. Branches always
        // overlap the previous course; the center stem connects the whole root.
        const double spread = std::min(width, 0.5 * taper_height);
        if (spread < 0.25 * width)
            continue;
        const double pitch = std::max(config.rooting_spacing.value, 2. * spread + width + skin);

        for (size_t upper = 0; upper < object.num_printing_regions(); ++upper) {
            if (m_slices[top][upper].empty() || !coherent_material(object.printing_region(upper), object.print()->config()))
                continue;
            const auto upper_tool = object.printing_region(upper).extruder(frExternalPerimeter);
            // Only new downward-facing material, not a continuing vertical seam.
            const ExPolygons new_material = diff_ex(m_slices[top][upper], m_slices[top - 1][upper]);
            if (new_material.empty())
                continue;
            for (size_t host = 0; host < object.num_printing_regions(); ++host) {
                const auto host_tool = object.printing_region(host).extruder(frExternalPerimeter);
                if (host == upper || host_tool == upper_tool || m_slices[top - 1][host].empty() ||
                    !coherent_material(object.printing_region(host), object.print()->config()))
                    continue;
                const auto &nozzles = object.print()->config().nozzle_diameter;
                if (width + EPSILON < 2. * std::max(nozzles.get_at(upper_tool - 1), nozzles.get_at(host_tool - 1)))
                    continue;
                const ExPolygons contact = intersection_ex(new_material, m_slices[top - 1][host]);
                if (contact.empty())
                    continue;

                // Entire root footprint stays under the contact, not underneath
                // exposed ledges. Intersect every native AND current host slice
                // to protect cavities, skins, other roots and existing beams.
                ExPolygons safe = offset_ex(contact, -scale_(skin));
                // Beam interlocking may also have changed the receiving layer.
                // A native contact alone does not guarantee a live attachment.
                safe = intersection_ex(safe, offset_ex(to_expolygons(layers[top]->get_region(upper)->slices.surfaces), -scale_(skin)));
                for (size_t z = floor; z < top && !safe.empty(); ++z) {
                    cancel();
                    safe = intersection_ex(safe, offset_ex(m_slices[z][host], -scale_(skin)));
                    safe = intersection_ex(safe, offset_ex(to_expolygons(layers[z]->get_region(host)->slices.surfaces), -scale_(skin)));
                }
                if (safe.empty())
                    continue;
                std::vector<Polygons> additions(top - bottom);
                for (const ExPolygon &island : safe) {
                    const BoundingBox bounds = get_extents(island);
                    const double cx = unscale_(bounds.center().x()), cy = unscale_(bounds.center().y());
                    const int nx = int(unscale_(bounds.size().x()) / (2. * pitch));
                    const int ny = int(unscale_(bounds.size().y()) / (2. * pitch));
                    for (int ix = -nx; ix <= nx; ++ix) {
                        cancel();
                        for (int iy = -ny; iy <= ny; ++iy) {
                            const double x = cx + ix * pitch, y = cy + iy * pitch;
                            // A conservative disk encloses every branch at every
                            // depth. Reject whole roots; never clip off an anchor.
                            const Polygons envelope{root_disk(x, y, (spread + width * 0.5) / std::cos(PI / 24.) + 0.005)};
                            if (!diff_ex(envelope, ExPolygons{island}).empty())
                                continue;
                            ++count;
                            for (size_t z = bottom; z < top; ++z) {
                                const double offset = spread * (contact_z - layers[z]->print_z) / taper_height;
                                Polygons &polys = additions[z - bottom];
                                polys.push_back(root_disk(x, y, width * 0.5));
                                for (const Vec2d &dir : {Vec2d(1., 0.), Vec2d(-1., 0.), Vec2d(0., 1.), Vec2d(0., -1.)})
                                    polys.push_back(root_disk(x + dir.x() * offset, y + dir.y() * offset, width * 0.5));
                            }
                        }
                    }
                }
                for (size_t z = bottom; z < top; ++z) {
                    cancel();
                    if (additions[z - bottom].empty())
                        continue;
                    const ExPolygons roots = union_ex(additions[z - bottom]);
                    auto &from = layers[z]->get_region(host)->slices;
                    auto &to = layers[z]->get_region(upper)->slices;
                    from.set(diff_ex(to_expolygons(from.surfaces), roots), stInternal);
                    to.set(union_ex(to_expolygons(to.surfaces), roots), stInternal);
                }
            }
        }
    }
    return count;
}
}
