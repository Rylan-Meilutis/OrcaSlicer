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

// Rasterize a tapered, swept sphere in horizontal courses. Unlike translating
// disks with Z, this includes the long horizontal portions of a shallow root.
// Adjacent samples overlap; the rounded parent node blends fork junctions.
std::vector<Polygons> root_branch(const Vec3d &a, const Vec3d &b, double radius_a, double radius_b,
    const std::vector<Layer *> &layers, size_t bottom, size_t top, double contact_z,
    const std::function<void()> &cancel)
{
    const int samples = std::max(1, int(std::ceil((b - a).norm() / (0.25 * radius_b))));
    std::vector<Polygons> result(top - bottom);
    for (size_t z = bottom; z < top; ++z) {
        cancel();
        Polygons disks;
        const double descent = contact_z - layers[z]->print_z;
        for (int i = 0; i <= samples; ++i) {
            const double t = double(i) / samples;
            const Vec3d center = a + t * (b - a);
            const double radius = radius_a + t * (radius_b - radius_a);
            const double dz = descent - center.z();
            const double r2 = radius * radius - dz * dz;
            if (r2 > EPSILON * EPSILON)
                disks.push_back(root_disk(center.x(), center.y(), std::sqrt(r2)));
        }
        result[z - bottom] = union_(disks);
    }
    return result;
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
    std::vector<ExPolygons> occupied(layers.size());
    for (size_t top = 1; top < layers.size(); ++top) {
        cancel();
        const double contact_z = layers[top - 1]->print_z;
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

                // Qualify against the native host, not beam-created pockets.
                // Reserve the local envelope for roots; elsewhere beams remain.
                // Exclude third materials explicitly. Requiring the post-beam
                // union to be watertight magnifies its rounding cracks by skin
                // and falsely rejects a broad collar inside a solid native host.
                const auto other_material = [&](size_t z) {
                    ExPolygons other;
                    for (size_t r = 0; r < object.num_printing_regions(); ++r)
                        if (r != host && r != upper)
                            append(other, to_expolygons(layers[z]->get_region(r)->slices.surfaces));
                    return other;
                };
                ExPolygons available = offset_ex(contact, -scale_(skin)), safe;
                size_t bottom = top, floor = top;
                // Walk down through actual host material. Retain the deepest
                // usable depth plus a complete protective floor; a shallow base
                // or cavity limits depth instead of rejecting all shorter roots.
                for (size_t z = top; z-- > 0 && !available.empty();) {
                    cancel();
                    available = intersection_ex(available, offset_ex(m_slices[z][host], -scale_(skin)));
                    available = diff_ex(available, offset_ex(other_material(z), scale_(skin)));
                    if (available.empty()) break;
                    const double floor_z = layers[z]->print_z - layers[z]->height;
                    size_t candidate = z;
                    while (candidate < top &&
                        (layers[candidate]->print_z - layers[candidate]->height < floor_z + skin - EPSILON ||
                         contact_z - (layers[candidate]->print_z - layers[candidate]->height) > depth + EPSILON))
                        ++candidate;
                    if (candidate < bottom && top - candidate >= 3) {
                        bottom = candidate;
                        floor = z;
                        safe = available;
                    }
                    if (contact_z - floor_z >= depth + skin - EPSILON) break;
                }
                if (safe.empty())
                    continue;
                const double taper_height = contact_z - layers[bottom]->print_z;
                if (taper_height < 0.5 * width) continue;
                // The collar stays under the upper part. Buried limbs may extend
                // beyond it, but only where a full native top/bottom skin remains.
                std::vector<ExPolygons> permitted(top - bottom);
                for (size_t z = bottom; z < top; ++z) {
                    cancel();
                    ExPolygons &area = permitted[z - bottom];
                    area = offset_ex(m_slices[z][host], -scale_(skin));
                    size_t lower = z;
                    const double course_bottom = layers[z]->print_z - layers[z]->height;
                    while (lower > 0 && course_bottom - (layers[lower]->print_z - layers[lower]->height) < skin - EPSILON)
                        --lower;
                    for (size_t k = lower; k < z && !area.empty(); ++k)
                        area = intersection_ex(area, offset_ex(m_slices[k][host], -scale_(skin)));
                    for (size_t k = z + 1; k < layers.size() && layers[k - 1]->print_z < layers[z]->print_z + skin - EPSILON && !area.empty(); ++k)
                        area = intersection_ex(area, offset_ex(union_ex(m_slices[k][host], m_slices[k][upper]), -scale_(skin)));
                    area = diff_ex(area, offset_ex(other_material(z), scale_(skin)));
                }
                const auto fits = [&](const std::vector<Polygons> &sections) {
                    for (size_t z = 0; z < sections.size(); ++z)
                        if (!diff_ex(sections[z], permitted[z]).empty()) return false;
                    return true;
                };
                const double pitch = std::max(config.rooting_spacing.value, 4. * width + skin);
                std::vector<Polygons> additions(top - bottom);
                const ExPolygons contact_collar = offset_ex(contact, -scale_(skin));
                ExPolygons collar;
                Points placed_centers;
                size_t placed = 0;
                for (const ExPolygon &island : safe) {
                    const BoundingBox bounds = get_extents(island);
                    const double cx = unscale_(bounds.center().x()), cy = unscale_(bounds.center().y());
                    const int nx = int(unscale_(bounds.size().x()) / (2. * pitch));
                    const int ny = int(unscale_(bounds.size().y()) / (2. * pitch));
                    Points candidates;
                    for (int ix = -nx; ix <= nx; ++ix)
                        for (int iy = -ny; iy <= ny; ++iy)
                            candidates.emplace_back(scale_(cx + ix * pitch), scale_(cy + iy * pitch));
                    // A concave or annular contact need not contain its bounding
                    // box center. Probe a bounded set of interior contour seeds.
                    for (const ExPolygon &interior : offset_ex(ExPolygons{island}, -scale_(width + 0.02))) {
                        const Points &points = interior.contour.points;
                        const size_t stride = std::max(size_t(1), (points.size() + 7) / 8);
                        for (size_t i = 0; i < points.size(); i += stride) candidates.push_back(points[i]);
                    }
                    for (const Point &candidate : candidates) {
                        cancel();
                        if (std::any_of(placed_centers.begin(), placed_centers.end(), [&](const Point &center) {
                            return unscale(Point(candidate - center)).norm() < pitch - EPSILON;
                        })) continue;
                        const double x = unscale_(candidate.x()), y = unscale_(candidate.y());
                        // Four primary limbs, then binary forks: at most 16
                        // tips. Shallower/smaller contacts use fewer levels.
                        for (int levels = 3; levels >= 1; --levels) {
                            const double radius = 0.5 * width * std::sqrt(4. * std::pow(2., levels - 1));
                            if (radius > taper_height) continue;
                            const Vec3d origin(x, y, 0.);
                            auto tree = root_branch(origin, origin, radius, radius, layers, bottom, top, contact_z, cancel);
                            if (!fits(tree)) continue;
                            std::function<bool(const Vec3d &, double, double, int)> grow;
                            grow = [&](const Vec3d &parent, double parent_radius, double heading, int level) {
                                cancel();
                                const double child_radius = parent_radius / (level == 0 ? 2. : std::sqrt(2.));
                                if (child_radius + EPSILON < 0.5 * width) return false;
                                const double descent = std::min(taper_height - child_radius,
                                    skin + child_radius + std::max(0., taper_height - skin - 2. * child_radius) * (level + 1.) / levels);
                                const double reach = taper_height * std::pow(0.72, level);
                                for (int attempt = 0; attempt < 6; ++attempt) {
                                    const double length = reach * std::pow(0.75, attempt);
                                    if (length < child_radius) break;
                                    const Vec3d child(parent.x() + length * std::cos(heading),
                                        parent.y() + length * std::sin(heading), std::max(parent.z(), descent));
                                    auto branch = root_branch(parent, child, parent_radius, child_radius,
                                        layers, bottom, top, contact_z, cancel);
                                    if (!fits(branch)) continue;
                                    for (size_t z = 0; z < tree.size(); ++z) append(tree[z], std::move(branch[z]));
                                    if (level + 1 < levels)
                                        for (double turn : {-PI / 5., PI / 5.})
                                            grow(child, child_radius, heading + turn, level + 1);
                                    return true;
                                }
                                return false;
                            };
                            int limbs = 0;
                            for (int direction = 0; direction < 4; ++direction)
                                limbs += grow(origin, radius, direction * PI / 2., 0);
                            if (limbs < 2) continue;
                            for (size_t z = 0; z < tree.size(); ++z) append(additions[z], std::move(tree[z]));
                            placed_centers.push_back(candidate);
                            for (const ExPolygon &part : contact_collar)
                                if (part.contains(candidate)) collar.push_back(part);
                            ++placed;
                            break;
                        }
                    }
                }
                if (placed == 0) continue;
                count += placed;
                // A conformal attachment follows the upper model (including
                // holes), not a cylinder. It blends into the buried root trees.
                collar = union_ex(collar);
                const double collar_depth = std::min(0.25 * taper_height, 0.5 * width);
                for (size_t z = bottom; z < top; ++z) {
                    const double descent = contact_z - layers[z]->print_z;
                    if (descent <= collar_depth + EPSILON)
                        append(additions[z - bottom], to_polygons(intersection_ex(
                            offset_ex(collar, -scale_(descent)), permitted[z - bottom])));
                    additions[z - bottom] = union_(additions[z - bottom]);
                }
                // Restore a host skin around accepted branches, without erasing
                // other roots. Beam pockets outside this local guard remain.
                for (size_t z = floor; z < top; ++z) {
                    cancel();
                    Polygons nearby;
                    for (size_t k = bottom; k < top; ++k)
                        if (std::abs(layers[k]->print_z - layers[z]->print_z) <= skin + EPSILON)
                            append(nearby, additions[k - bottom]);
                    auto guard = intersection_ex(offset_ex(union_ex(nearby), scale_(skin)), m_slices[z][host]);
                    guard = diff_ex(guard, union_ex(other_material(z), occupied[z]));
                    auto &from = layers[z]->get_region(host)->slices;
                    auto &to = layers[z]->get_region(upper)->slices;
                    from.set(union_ex(to_expolygons(from.surfaces), guard), stInternal);
                    to.set(diff_ex(to_expolygons(to.surfaces), guard), stInternal);
                }
                for (size_t z = top; z < layers.size(); ++z) {
                    cancel();
                    const auto attachment = diff_ex(intersection_ex(collar, m_slices[z][upper]), other_material(z));
                    if (attachment.empty()) break;
                    auto &from = layers[z]->get_region(host)->slices;
                    auto &to = layers[z]->get_region(upper)->slices;
                    from.set(diff_ex(to_expolygons(from.surfaces), attachment), stInternal);
                    to.set(union_ex(to_expolygons(to.surfaces), attachment), stInternal);
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
                    occupied[z] = union_ex(occupied[z], roots);
                }
            }
        }
    }
    return count;
}
}
