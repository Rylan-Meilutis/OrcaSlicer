#include "Rooting.hpp"
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>
#include "libslic3r/Print.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include <algorithm>
#include <cmath>
#include <tbb/parallel_for.h>

namespace Slic3r {

Rooting::Rooting(const PrintObject &object)
{
    if ((!object.config().rooting && !object.config().rooting_xy) || object.num_printing_regions() < 2 || object.print()->config().spiral_mode)
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
    const std::function<void()> &cancel, double vertical_scale = 1.)
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
            const double dz = (descent - center.z()) / vertical_scale;
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

double root_wall_clearance(const PrintObject &object, const PrintRegion &region, const Layer &layer)
{
    const Flow outer = region.flow(object, frExternalPerimeter, layer.height, layer.id() == 0);
    const Flow inner = region.flow(object, frPerimeter, layer.height, layer.id() == 0);
    // Keep the native wall stack plus a printable host web. A fixed rooting
    // skin alone can put another material inside thick perimeter stacks.
    return std::max(object.config().rooting_skin.value,
        double(outer.width()) + std::max(0, region.config().wall_loops.value - 1) * inner.spacing() + inner.spacing());
}
}

size_t Rooting::apply(PrintObject &object, const std::function<void()> &cancel, const Progress &progress) const
{
    if (m_slices.empty())
        return 0;
    const auto &config = object.config();
    const double depth = config.rooting_depth.value;
    const double width = config.rooting_width.value;
    const double base_width = std::max(width, config.rooting_base_width.value);
    const double skin = config.rooting_skin.value;
    if (!std::isfinite(depth) || !std::isfinite(width) || !std::isfinite(skin) ||
        depth < 0.4 || depth > 20. || width < 0.4 || width > 10. || skin < 0.4 || skin > 10. ||
        !std::isfinite(config.rooting_base_width.value) || config.rooting_base_width.value < 0.4 ||
        config.rooting_base_width.value > 10. ||
        !std::isfinite(config.rooting_spacing.value) || config.rooting_spacing.value < 1.)
        return 0;

    size_t count = 0;
    const auto &layers = object.layers();
    // Freeze beam additions before any roots are applied. Later root contacts
    // must not mistake a previously grown root for part of the beam joint.
    std::vector<std::vector<ExPolygons>> beams(layers.size());
    if (config.interlocking_beam)
        for (size_t z = 0; z < layers.size(); ++z) {
            cancel();
            for (size_t r = 0; r < object.num_printing_regions(); ++r)
                beams[z].push_back(diff_ex(to_expolygons(layers[z]->get_region(r)->slices.surfaces), m_slices[z][r]));
        }
    if (!config.rooting)
        return apply_xy(object, cancel, beams, progress);
    std::vector<ExPolygons> occupied(layers.size());
    for (size_t interface_top = 1; interface_top < layers.size(); ++interface_top) {
        cancel();
        if (progress) progress(false, interface_top - 1, layers.size() - 1);
        for (size_t upper = 0; upper < object.num_printing_regions(); ++upper) {
            if (m_slices[interface_top][upper].empty() || !coherent_material(object.printing_region(upper), object.print()->config()))
                continue;
            const auto upper_tool = object.printing_region(upper).extruder(frExternalPerimeter);
            // Only new downward-facing material, not a continuing vertical seam.
            const ExPolygons new_material = diff_ex(m_slices[interface_top][upper], m_slices[interface_top - 1][upper]);
            if (new_material.empty())
                continue;
            for (size_t host = 0; host < object.num_printing_regions(); ++host) {
                size_t top = interface_top;
                const auto host_tool = object.printing_region(host).extruder(frExternalPerimeter);
                if (host == upper || host_tool == upper_tool || m_slices[top - 1][host].empty() ||
                    !coherent_material(object.printing_region(host), object.print()->config()))
                    continue;
                const auto &nozzles = object.print()->config().nozzle_diameter;
                const double minimum_radius = 0.5 * width;
                const double attachment_radius = std::min(minimum_radius,
                    std::max(nozzles.get_at(upper_tool - 1), nozzles.get_at(host_tool - 1)));
                const auto host_clearance = [&](size_t z) {
                    return root_wall_clearance(object, object.printing_region(host), *layers[z]);
                };
                if (width + EPSILON < 2. * std::max(nozzles.get_at(upper_tool - 1), nozzles.get_at(host_tool - 1)))
                    continue;
                ExPolygons contact = intersection_ex(new_material, m_slices[top - 1][host]);
                if (contact.empty())
                    continue;
                bool below_beams = false;
                if (config.interlocking_beam) {
                    ExPolygons bottom_beams;
                    for (size_t z = interface_top; z-- > 0;) {
                        cancel();
                        const auto overlap = intersection_ex(beams[z][upper], contact);
                        if (overlap.empty()) break;
                        top = z;
                        bottom_beams = overlap;
                        below_beams = true;
                    }
                    if (below_beams) {
                        // A broad collar below the last beam course joins its
                        // teeth without replacing the alternating beam band.
                        ExPolygons attached;
                        for (const auto &part : contact)
                            if (!intersection_ex(offset_ex(ExPolygons{part}, -scale_(skin)), bottom_beams).empty())
                                attached.push_back(part);
                        contact = std::move(attached);
                    }
                }
                if (top == 0 || contact.empty()) continue;
                const double contact_z = layers[top - 1]->print_z;

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
                    available = intersection_ex(available, offset_ex(m_slices[z][host], -scale_(host_clearance(z))));
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
                if (taper_height <= EPSILON) continue;
                const double vertical_scale = std::min(1., taper_height / base_width);
                // The collar stays under the upper part. Buried limbs may extend
                // beyond it, but only where a full native top/bottom skin remains.
                std::vector<ExPolygons> permitted(top - bottom);
                tbb::parallel_for(bottom, top, [&](size_t z) {
                    cancel();
                    ExPolygons &area = permitted[z - bottom];
                    area = offset_ex(m_slices[z][host], -scale_(host_clearance(z)));
                    size_t lower = z;
                    const double course_bottom = layers[z]->print_z - layers[z]->height;
                    while (lower > 0 && course_bottom - (layers[lower]->print_z - layers[lower]->height) < skin - EPSILON)
                        --lower;
                    for (size_t k = lower; k < z && !area.empty(); ++k)
                        area = intersection_ex(area, offset_ex(m_slices[k][host], -scale_(host_clearance(k))));
                    for (size_t k = z + 1; k < layers.size() && layers[k - 1]->print_z < layers[z]->print_z + skin - EPSILON && !area.empty(); ++k)
                        area = intersection_ex(area, offset_ex(union_ex(m_slices[k][host], m_slices[k][upper]), -scale_(skin)));
                    area = diff_ex(area, offset_ex(other_material(z), scale_(skin)));
                });
                const auto fits = [&](const std::vector<Polygons> &sections) {
                    for (size_t z = 0; z < sections.size(); ++z)
                        if (!diff_ex(sections[z], permitted[z]).empty()) return false;
                    return true;
                };
                const double pitch = std::max(config.rooting_spacing.value, 4. * width + skin);
                std::vector<Polygons> additions(top - bottom);
                const double collar_depth = std::min(0.25 * taper_height, 0.5 * width);
                const double web_radius = 0.5 * std::max(skin, 2. * nozzles.get_at(host_tool - 1));
                // Preserve the host's existing topology, including a printable
                // web between neighboring limbs. Testing an eroded host also
                // catches connections that are only a hairline wide.
                std::vector<ExPolygons> host_sections(top - bottom);
                std::vector<std::vector<size_t>> host_core_counts(top - bottom);
                tbb::parallel_for(bottom, top, [&](size_t z) {
                    cancel();
                    host_sections[z - bottom] = diff_ex(m_slices[z][host], occupied[z]);
                    for (const auto &part : host_sections[z - bottom])
                        host_core_counts[z - bottom].push_back(offset_ex(ExPolygons{part}, -scale_(web_radius)).size());
                });
                const auto preserves_host_webs = [&](const std::vector<Polygons> &tree,
                                                     const std::vector<Polygons> &branch) {
                    for (size_t z = bottom; z < top; ++z) {
                        cancel();
                        // The model-shaped attachment is intentionally a plate;
                        // enforce webs in the load-bearing root network below it.
                        if (contact_z - layers[z]->print_z <= collar_depth + EPSILON) continue;
                        const size_t i = z - bottom;
                        Polygons proposed = additions[i];
                        append(proposed, tree[i]);
                        append(proposed, branch[i]);
                        const auto roots = union_ex(proposed);
                        for (size_t h = 0; h < host_sections[i].size(); ++h) {
                            const auto &host_part = host_sections[i][h];
                            if (intersection_ex(ExPolygons{host_part}, roots).empty()) continue;
                            const auto remaining = diff_ex(ExPolygons{host_part}, roots);
                            if (remaining.size() > 1 || remaining.empty()) return false;
                            const auto core = offset_ex(remaining, -scale_(web_radius));
                            if (core.empty() || core.size() > host_core_counts[i][h]) return false;
                        }
                    }
                    return true;
                };
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
                    for (double seed_radius : {width, minimum_radius, attachment_radius})
                        for (const ExPolygon &interior : offset_ex(ExPolygons{island}, -scale_(seed_radius + 0.02))) {
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
                        // Fit the trunk independently of branching depth. A
                        // narrow contact must not lose its lateral roots merely
                        // because an area-preserving trunk for every tip is wide.
                        for (int levels = 3; levels >= 1; --levels) {
                            double radius = std::min(0.5 * base_width, taper_height / vertical_scale);
                            if (radius < 0.5 * width) continue;
                            const Vec3d origin(x, y, 0.);
                            auto tree = root_branch(origin, origin, radius, radius, layers, bottom, top, contact_z, cancel, vertical_scale);
                            const std::vector<Polygons> empty_tree(top - bottom);
                            // Annular contacts around screw holes can fit a
                            // printable branch but not the preferred broad
                            // trunk. Retain the conformal collar and shorten
                            // the attachment locally. This must not reduce the
                            // requested diameter of the buried branches.
                            while (radius > attachment_radius + EPSILON &&
                                   (!fits(tree) || !preserves_host_webs(empty_tree, tree))) {
                                radius = std::max(attachment_radius, radius * 0.75);
                                tree = root_branch(origin, origin, radius, radius,
                                    layers, bottom, top, contact_z, cancel, vertical_scale);
                            }
                            if (!fits(tree) || !preserves_host_webs(empty_tree, tree)) continue;
                            std::function<bool(const Vec3d &, double, double, int)> grow;
                            grow = [&](const Vec3d &parent, double parent_radius, double heading, int level) {
                                cancel();
                                const double child_radius = 0.5 * (base_width +
                                    (width - base_width) * double(level + 1) / levels);
                                const double vertical_radius = child_radius * vertical_scale;
                                const double descent = std::min(taper_height - vertical_radius,
                                    skin + vertical_radius + std::max(0., taper_height - skin - 2. * vertical_radius) * (level + 1.) / levels);
                                const double reach = std::max(taper_height, 2. * width) * std::pow(0.72, level);
                                for (int attempt = 0; attempt < 6; ++attempt) {
                                    const double length = reach * std::pow(0.75, attempt);
                                    if (length < child_radius) break;
                                    const Vec3d child(parent.x() + length * std::cos(heading),
                                        parent.y() + length * std::sin(heading), std::max(parent.z(), descent));
                                    auto branch = root_branch(parent, child, parent_radius, child_radius,
                                        layers, bottom, top, contact_z, cancel, vertical_scale);
                                    // Trim only the side that meets a protected
                                    // wall. Rejecting the entire swept branch
                                    // makes a wider root disappear around bores.
                                    // Keep only geometry connected upward to the
                                    // existing tree; clipping must not leave loose
                                    // islands in the host.
                                    bool extends_tree = true;
                                    if (!fits(branch)) {
                                        ExPolygons above = union_ex(tree.back());
                                        for (size_t z = branch.size(); z-- > 0;) {
                                            const auto clipped = intersection_ex(branch[z], permitted[z]);
                                            const auto attachment = union_ex(above, union_ex(tree[z]));
                                            branch[z].clear();
                                            for (const auto &part : clipped)
                                                if (!intersection_ex(ExPolygons{part}, attachment).empty())
                                                    append(branch[z], to_polygons(part));
                                            above = union_ex(union_ex(tree[z]), branch[z]);
                                        }
                                        extends_tree = false;
                                        for (size_t z = 0; z < tree.size() && !extends_tree; ++z)
                                            extends_tree = !diff_ex(branch[z], tree[z]).empty();
                                    }
                                    // A later fork can originate in a trimmed
                                    // portion of its parent, even if the fork
                                    // itself fits. Check untrimmed forks too.
                                    ExPolygons above = union_ex(tree.back());
                                    for (size_t z = tree.size(); z-- > 0 && extends_tree;) {
                                        auto combined = union_ex(union_ex(tree[z]), branch[z]);
                                        for (const auto &part : combined)
                                            if (intersection_ex(ExPolygons{part}, above).empty()) {
                                                extends_tree = false;
                                                break;
                                            }
                                        above = std::move(combined);
                                    }
                                    if (!extends_tree || !preserves_host_webs(tree, branch)) continue;
                                    for (size_t z = 0; z < tree.size(); ++z) append(tree[z], std::move(branch[z]));
                                    if (level + 1 < levels) {
                                        // Side roots emerge along the limb, not
                                        // only at its terminal tip. Keep the main
                                        // limb intact and taper the offshoots to
                                        // the configured printable minimum.
                                        for (double turn : {-PI / 3., PI / 3.}) {
                                            const double t = turn < 0. ? 0.55 : 0.8;
                                            const Vec3d fork = parent + t * (child - parent);
                                            const double fork_radius = parent_radius + t * (child_radius - parent_radius);
                                            grow(fork, fork_radius, heading + turn, level + 1);
                                        }
                                    }
                                    return true;
                                }
                                return false;
                            };
                            int limbs = 0;
                            for (int direction = 0; direction < 4; ++direction)
                                limbs += grow(origin, radius, direction * PI / 2., 0);
                            // A wall may block some directions without invalidating
                            // the printable, connected limbs growing away from it.
                            if (limbs == 0) continue;
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
                    if (config.interlocking_beam)
                        guard = diff_ex(guard, beams[z][upper]);
                    auto &from = layers[z]->get_region(host)->slices;
                    auto &to = layers[z]->get_region(upper)->slices;
                    from.set(union_ex(to_expolygons(from.surfaces), guard), stInternal);
                    to.set(diff_ex(to_expolygons(to.surfaces), guard), stInternal);
                }
                for (size_t z = top; !below_beams && z < layers.size(); ++z) {
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
    if (progress) progress(false, layers.size() - 1, layers.size() - 1);
    return count + apply_xy(object, cancel, beams, progress);
}

size_t Rooting::apply_xy(PrintObject &object, const std::function<void()> &cancel,
                        const std::vector<std::vector<ExPolygons>> &beams, const Progress &progress) const
{
    const auto &config = object.config();
    if (!config.rooting_xy) return 0;
    const auto &layers = object.layers();
    const auto &print_config = object.print()->config();
    const double width = config.rooting_width.value, skin = config.rooting_skin.value;
    const double pitch = std::max(config.rooting_spacing.value, 4. * width + skin);
    const double radius = 0.5 * std::max(width, config.rooting_base_width.value);
    size_t count = 0;
    std::vector<ExPolygons> occupied(layers.size());
    for (size_t donor = 0; donor < object.num_printing_regions(); ++donor) {
        if (!coherent_material(object.printing_region(donor), print_config)) continue;
        const auto donor_tool = object.printing_region(donor).extruder(frExternalPerimeter);
        for (size_t host = 0; host < object.num_printing_regions(); ++host) {
            cancel();
            if (progress) progress(true, donor * object.num_printing_regions() + host,
                                   object.num_printing_regions() * object.num_printing_regions());
            if (host == donor || !coherent_material(object.printing_region(host), print_config)) continue;
            const auto host_tool = object.printing_region(host).extruder(frExternalPerimeter);
            if (host_tool == donor_tool || width + EPSILON < 2. * std::max(
                print_config.nozzle_diameter.get_at(donor_tool - 1),
                print_config.nozzle_diameter.get_at(host_tool - 1))) continue;
            const double web_radius = 0.5 * std::max(skin, 2. * print_config.nozzle_diameter.get_at(host_tool - 1));
            const auto native_envelope = [&](size_t z) {
                const double clearance = std::max(
                    root_wall_clearance(object, object.printing_region(host), *layers[z]),
                    root_wall_clearance(object, object.printing_region(donor), *layers[z]));
                // Independently transformed meshes can leave sub-micron cracks
                // at a shared contact. Close only rounding-scale gaps before
                // erosion so they are not inflated into a full skin-width void.
                // Transfers remain clipped to actual host material below.
                return offset2_ex(union_ex(m_slices[z][donor], m_slices[z][host]),
                    scale_(0.001), -scale_(clearance + 0.001));
            };
            // Neighboring courses reuse the same clearance envelope. Compute it
            // once per layer; each worker owns its output and reads native slices.
            std::vector<ExPolygons> envelopes(layers.size());
            tbb::parallel_for(size_t(0), layers.size(), [&](size_t z) {
                cancel();
                envelopes[z] = native_envelope(z);
            });
            std::vector<ExPolygons> permitted(layers.size());
            tbb::parallel_for(size_t(0), layers.size(), [&](size_t z) {
                cancel();
                const double low = layers[z]->print_z - layers[z]->height - skin;
                const double high = layers[z]->print_z + skin;
                if (low < layers.front()->print_z - layers.front()->height - EPSILON ||
                    high > layers.back()->print_z + EPSILON) return;
                auto area = envelopes[z];
                // Require a complete native skin above and below each course,
                // not merely clearance from the object's global bounding box.
                for (size_t k = z; k-- > 0 && layers[k]->print_z > low + EPSILON && !area.empty();)
                    area = intersection_ex(area, envelopes[k]);
                for (size_t k = z + 1; k < layers.size() && layers[k]->print_z - layers[k]->height < high - EPSILON && !area.empty(); ++k)
                    area = intersection_ex(area, envelopes[k]);
                permitted[z] = std::move(area);
            });
            envelopes.clear(); // Do not retain the clearance cache while growing trees.
            double next_z = layers.front()->print_z - layers.front()->height + skin + radius;
            for (size_t center = 0; center < layers.size(); ++center) {
                cancel();
                if (layers[center]->print_z + EPSILON < next_z || m_slices[center][donor].empty() ||
                    m_slices[center][host].empty()) continue;
                // Embed the smaller adjoining part into the larger base, not
                // the reverse. Equal areas have a stable region-order tie break.
                const double donor_area = area(m_slices[center][donor]), host_area = area(m_slices[center][host]);
                if (donor_area > host_area || (donor_area == host_area && donor > host)) continue;
                const size_t count_before = count;
                const double center_z = layers[center]->print_z;
                size_t bottom = center, top = center + 1;
                while (bottom > 0 && center_z - layers[bottom - 1]->print_z < radius) --bottom;
                while (top < layers.size() && layers[top]->print_z - center_z < radius) ++top;
                ExPolygons beam_band;
                if (config.interlocking_beam)
                    for (size_t z = bottom; z < top; ++z)
                        append(beam_band, intersection_ex(beams[z][donor], m_slices[z][host]));
                beam_band = union_ex(beam_band);
                // Sample the outer ends of the existing beam teeth, not an
                // artificial contact inside their alternating band.
                const auto attachment = offset2_ex(union_ex(m_slices[center][donor], beam_band),
                    scale_(0.001), -scale_(0.001));
                Points placed;
                for (const Polygon &contour : to_polygons(attachment)) {
                    for (const Line &edge : contour.lines()) {
                        const Vec2d a = unscale(edge.a), b = unscale(edge.b);
                        const Vec2d delta = b - a;
                        const double length = delta.norm();
                        // Curved contacts are tessellated into edges much
                        // shorter than a root. Their edge length does not limit
                        // the available attachment width; the sweep fit does.
                        if (length < EPSILON) continue;
                        const Vec2d tangent = delta / length;
                        for (double s = std::min(length * 0.5, pitch * 0.5); s < length; s += pitch) {
                            cancel();
                            const Vec2d p = a + tangent * s;
                            const Point seed(scale_(p.x()), scale_(p.y()));
                            if (std::any_of(placed.begin(), placed.end(), [&](const Point &q) {
                                return unscale(Point(seed - q)).norm() < pitch - EPSILON;
                            })) continue;
                            Vec2d normal(-tangent.y(), tangent.x());
                            const auto contains = [](const ExPolygons &area, const Vec2d &q) {
                                const Point point(scale_(q.x()), scale_(q.y()));
                                return std::any_of(area.begin(), area.end(), [&](const ExPolygon &part) { return part.contains(point); });
                            };
                            if (contains(attachment, p + normal * 0.02)) normal = -normal;
                            if (!contains(attachment, p - normal * 0.02) ||
                                !contains(m_slices[center][host], p + normal * 0.02)) continue;
                            // Beam teeth are diagonal to the native material
                            // contact. Their edge normals point along the teeth,
                            // leaving little room for a buried branching fan.
                            // Grow into the host normal to the original contact.
                            Vec2d branch_tangent = tangent;
                            if (!beam_band.empty()) {
                                double nearest = std::numeric_limits<double>::max();
                                for (const Polygon &native_contour : to_polygons(m_slices[center][donor]))
                                    for (const Line &native_edge : native_contour.lines()) {
                                        const double distance = native_edge.distance_to_squared(seed);
                                        const Vec2d direction = unscale(Point(native_edge.b - native_edge.a));
                                        if (distance < nearest && direction.norm() > EPSILON) {
                                            nearest = distance;
                                            branch_tangent = direction.normalized();
                                        }
                                    }
                                Vec2d inward(-branch_tangent.y(), branch_tangent.x());
                                if (inward.dot(normal) < 0.) inward = -inward;
                                normal = inward;
                            }
                            std::vector<ExPolygons> from, to;
                            std::vector<std::vector<size_t>> core_counts;
                            for (size_t z = bottom; z < top; ++z) {
                                from.push_back(to_expolygons(layers[z]->get_region(host)->slices.surfaces));
                                to.push_back(to_expolygons(layers[z]->get_region(donor)->slices.surfaces));
                                auto &counts = core_counts.emplace_back();
                                for (const auto &part : from.back())
                                    counts.push_back(offset_ex(ExPolygons{part}, -scale_(web_radius)).size());
                            }
                            const auto transferable = [&](const Polygons &section, size_t i) {
                                return diff_ex(intersection_ex(intersection_ex(section, from[i]),
                                    m_slices[bottom + i][host]), beam_band);
                            };
                            const auto fits = [&](const std::vector<Polygons> &sections) {
                                for (size_t i = 0; i < sections.size(); ++i) {
                                    cancel();
                                    if (!diff_ex(sections[i], permitted[bottom + i]).empty() ||
                                        !intersection_ex(sections[i], occupied[bottom + i]).empty()) return false;
                                    // Check exactly the canonical geometry that
                                    // will be transferred, including rounding at
                                    // overlapping branch/beam intersections.
                                    const auto roots = transferable(union_(sections[i]), i);
                                    if (diff_ex(from[i], roots).size() > from[i].size()) return false;
                                    for (size_t h = 0; h < from[i].size(); ++h) {
                                        const auto &part = from[i][h];
                                        if (intersection_ex(ExPolygons{part}, roots).empty()) continue;
                                        const auto remaining = diff_ex(ExPolygons{part}, roots);
                                        if (remaining.size() != 1) return false;
                                        const auto core = offset_ex(remaining, -scale_(web_radius));
                                        if (core.empty() || core.size() > core_counts[i][h])
                                            return false;
                                    }
                                }
                                return true;
                            };
                            bool accepted = false;
                            for (int attempt = 0; attempt < 6 && !accepted; ++attempt) {
                                const double depth = config.rooting_depth.value * std::pow(0.75, attempt);
                                if (depth < 2. * width) break;
                                const Vec2d root_start = p - normal * (0.25 * width);
                                // Reach inside the host before spreading sideways.
                                // Near a corner, a fan starting halfway down the
                                // trunk can hit the wall even on its open side.
                                double trunk_length = depth - radius;
                                // The fan must clear neighboring teeth as well
                                // as the one carrying its attachment. Otherwise
                                // its shoulders cut pockets off between teeth.
                                for (const Polygon &band : to_polygons(beam_band))
                                    for (const Point &point : band.points) {
                                        const Vec2d relative = unscale(point) - p;
                                        if (std::abs(relative.dot(branch_tangent)) < depth + radius)
                                            trunk_length = std::max(trunk_length,
                                                relative.dot(normal) + radius * 0.8 + web_radius);
                                    }
                                if (trunk_length > depth - radius) continue;
                                const Vec2d root_end = p + normal * trunk_length;
                                const Vec3d origin(root_start.x(), root_start.y(), 0.);
                                const Vec3d trunk(root_end.x(), root_end.y(), 0.);
                                // Match the neck to a beam tooth, then broaden
                                // inside the base. A full-width neck severs the
                                // neighboring host teeth from their material.
                                const double neck_radius = beam_band.empty() ? radius :
                                    std::min(radius, std::max(0.5 * config.interlocking_beam_width.value,
                                        print_config.nozzle_diameter.get_at(donor_tool - 1)));
                                auto tree = root_branch(origin, trunk, neck_radius, radius,
                                    layers, bottom, top, center_z, cancel);
                                for (size_t i = 0; i < tree.size(); ++i)
                                    tree[i] = to_polygons(intersection_ex(tree[i], permitted[bottom + i]));
                                if (!fits(tree)) continue;
                                std::function<bool(const Vec3d &, double, double, int)> grow;
                                grow = [&](const Vec3d &parent, double parent_radius, double side, int level) {
                                    const double child_radius = radius + (0.5 * width - radius) * double(level + 1) / 3.;
                                    const double reach = depth * 0.65 * std::pow(0.65, level);
                                    for (int shorten = 0; shorten < 6; ++shorten) {
                                        const double scale = std::pow(0.75, shorten);
                                        const double forward = std::max(0., std::min((depth - radius) * 0.15 * scale,
                                            depth - child_radius - (parent.head<2>() - p).dot(normal)));
                                        const Vec2d xy = parent.head<2>() + normal * forward +
                                            branch_tangent * (side * reach * scale);
                                        const Vec3d child(xy.x(), xy.y(), 0.);
                                        if ((child - parent).norm() < (level == 0 ? width * 1.25 : child_radius)) break;
                                        auto branch = root_branch(parent, child, parent_radius, child_radius,
                                            layers, bottom, top, center_z, cancel);
                                        for (size_t i = 0; i < branch.size(); ++i)
                                            branch[i] = to_polygons(intersection_ex(branch[i], permitted[bottom + i]));
                                        auto proposed = tree;
                                        for (size_t i = 0; i < tree.size(); ++i) append(proposed[i], branch[i]);
                                        if (!fits(proposed)) continue;
                                        tree = std::move(proposed);
                                        if (level < 2) {
                                            const Vec3d fork = parent + 0.65 * (child - parent);
                                            const double fork_radius = parent_radius + 0.65 * (child_radius - parent_radius);
                                            grow(fork, fork_radius, -side, level + 1);
                                            grow(child, child_radius, side, level + 1);
                                        }
                                        return true;
                                    }
                                    return false;
                                };
                                const bool left = grow(trunk, radius, -1., 0);
                                const bool right = grow(trunk, radius, 1., 0);
                                if (!left && !right) continue;
                                std::vector<ExPolygons> additions(tree.size());
                                std::vector<std::vector<bool>> connected(tree.size());
                                size_t pieces = 0;
                                for (size_t i = 0; i < tree.size(); ++i) {
                                    additions[i] = transferable(union_(tree[i]), i);
                                    for (const auto &part : additions[i])
                                        connected[i].push_back(!intersection_ex(offset_ex(ExPolygons{part}, scale_(0.01)), to[i]).empty());
                                    pieces += additions[i].size();
                                }
                                if (pieces == 0) continue;
                                // Every transferred component must reach the donor,
                                // directly or through adjacent courses. This permits
                                // attachment to alternating beam teeth without
                                // leaving floating material islands in other courses.
                                bool changed = true;
                                while (changed) {
                                    cancel();
                                    changed = false;
                                    for (size_t i = 0; i < additions.size(); ++i)
                                        for (size_t j = 0; j < additions[i].size(); ++j) {
                                            if (connected[i][j]) continue;
                                            for (int step : {-1, 1}) {
                                                const int k = int(i) + step;
                                                if (k < 0 || k >= int(additions.size())) continue;
                                                for (size_t n = 0; n < additions[k].size(); ++n)
                                                    if (connected[k][n] && !intersection_ex(
                                                        ExPolygons{additions[i][j]}, ExPolygons{additions[k][n]}).empty()) {
                                                        connected[i][j] = changed = true;
                                                        break;
                                                    }
                                            }
                                        }
                                }
                                if (std::any_of(connected.begin(), connected.end(), [](const auto &course) {
                                    return std::find(course.begin(), course.end(), false) != course.end();
                                })) continue;
                                for (size_t i = 0; i < additions.size(); ++i) {
                                    layers[bottom + i]->get_region(host)->slices.set(diff_ex(from[i], additions[i]), stInternal);
                                    layers[bottom + i]->get_region(donor)->slices.set(union_ex(to[i], additions[i]), stInternal);
                                    // Reserve the attachment as well as the roots so
                                    // later trees cannot cut off an earlier network.
                                    occupied[bottom + i] = union_ex(occupied[bottom + i],
                                        offset_ex(union_ex(tree[i]), scale_(skin)));
                                }
                                placed.push_back(seed);
                                ++count;
                                accepted = true;
                            }
                        }
                    }
                }
                if (count > count_before) next_z = center_z + pitch;
            }
        }
    }
    if (progress) progress(true, object.num_printing_regions() * object.num_printing_regions(),
                           object.num_printing_regions() * object.num_printing_regions());
    return count;
}
}
