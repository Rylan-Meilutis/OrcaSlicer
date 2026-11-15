#include "ThermoelasticSectionTopology.hpp"
#include "BoundingBox.hpp"
#include "ClipperUtils.hpp"
#include "ExPolygon.hpp"
#include "Point.hpp"
#include "Polygon.hpp"
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticExactSections.hpp"
#include "ThermoelasticSectionFootprint.hpp"
#include "ThermoelasticMesh.hpp"
#include "libslic3r.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
namespace {
long long height(double metres)
{
    if (!std::isfinite(metres) || std::abs(metres) > 10.)
        throw std::invalid_argument("Invalid section topology height.");
    return std::llround(metres * 1e9);
}
std::size_t point_count(const ExPolygon &polygon)
{
    auto count = polygon.contour.points.size();
    for (const auto &hole : polygon.holes) count += hole.points.size();
    return count;
}
void validate(const ExPolygon &polygon)
{
    if (!polygon.is_valid() || polygon.area() <= 0.)
        throw std::invalid_argument("Invalid section topology polygon.");
    // A pinched contour can join two bodies through a single vertex. It is
    // not a positive-area load path, even if a polygon library retains it.
    std::set<std::pair<coord_t, coord_t>> seen;
    const auto ring = [&](const Polygon &boundary) {
        for (const auto &point : boundary.points) {
            if (std::abs(double(point.x())) > 1e10 || std::abs(double(point.y())) > 1e10 ||
                !seen.emplace(point.x(), point.y()).second)
                throw std::invalid_argument("Unresolved point contact in section topology.");
        }
    };
    ring(polygon.contour);
    for (const auto &hole : polygon.holes) ring(hole);
}
}

SectionTopology qualify_section_topology(const std::vector<DepositedSlab> &slabs,
    double bed_z_m, const SectionTopologyLimits &limits, const std::function<void()> &checkpoint)
{
    if (slabs.empty() || limits.max_slabs == 0 || limits.max_slabs > 10000 ||
        limits.max_boundary_points == 0 || limits.max_boundary_points > 2000000 ||
        limits.max_component_pairs == 0 || limits.max_component_pairs > 2000000)
        throw std::invalid_argument("Invalid section topology input or budget.");
    if (slabs.size() > limits.max_slabs) throw ResourceLimit("Section topology slab budget exceeded.");
    const auto bed = height(bed_z_m);
    using Interval = std::pair<long long, long long>;
    std::map<Interval, std::map<std::size_t, std::vector<SectionFootprint>>> grouped;
    std::size_t input_points = 0, work = 0;
    const auto charge = [&](std::size_t amount) {
        if (checkpoint) checkpoint();
        if (amount > limits.max_component_pairs-work)
            throw ResourceLimit("Section topology work budget exceeded.");
        work += amount;
    };
    const auto check_work = [&] { charge(1); };
    const auto remaining = [&] {
        if (work == limits.max_component_pairs) throw ResourceLimit("Exact section topology work budget exhausted.");
        return std::min<std::size_t>(4000000,limits.max_component_pairs-work);
    };
    const auto promote = [&](const SectionFootprint& footprint) {
        if (footprint.is_exact()) return footprint.exact();
        auto converted = make_exact_section(footprint.native(),remaining(),checkpoint);
        charge(converted.work);
        return std::move(converted.section);
    };
    const auto exact_boolean = [&](const ExactSection& a,const ExactSection& b,ExactSectionOperation operation) {
        auto combined = boolean_exact_sections(a,b,operation,remaining(),checkpoint);
        charge(combined.work);
        return std::move(combined.section);
    };
    const bool exact_geometry = std::any_of(slabs.begin(),slabs.end(),[](const DepositedSlab& slab) {
        return slab.footprint.is_exact();
    });
    for (const auto &slab : slabs) {
        if (checkpoint) checkpoint();
        if (slab.boundary_uncertainty_nm)
            throw std::invalid_argument("Boundary uncertainty requires local topology contact qualification.");
        const auto bottom = height(slab.bottom_m), top = height(slab.top_m);
        if (bottom < bed || top <= bottom)
            throw std::invalid_argument("Invalid section topology slab interval.");
        auto &region = grouped[{bottom, top}][slab.material];
        if (slab.footprint.is_exact()) {
            const auto count = slab.footprint.vertex_count();
            if (count > limits.max_boundary_points-input_points)
                throw ResourceLimit("Section topology input point budget exceeded.");
            input_points += count;
            if (!limits.exact_interior_components)
                charge(validate_exact_section_topology(slab.footprint.exact(),remaining(),checkpoint));
            region.push_back(slab.footprint);
            continue;
        }
        for (const auto &polygon : slab.footprint.native()) {
            const auto count = point_count(polygon);
            if (count > limits.max_boundary_points - input_points)
                throw ResourceLimit("Section topology input point budget exceeded.");
            input_points += count;
            validate(polygon);
            region.emplace_back(ExPolygons{polygon});
        }
    }
    const auto normalize = [&](const std::vector<SectionFootprint>& parts) {
        std::vector<SectionFootprint> normalized;
        if (exact_geometry) {
            ExactSection joined;
            for (const auto& part : parts)
                joined = exact_boolean(joined,promote(part),ExactSectionOperation::Union);
            // Canonical inputs are certified by face extraction below in the
            // explicit interior policy; the strict default keeps ring checks.
            if (!limits.exact_interior_components)
                charge(validate_exact_section_topology(joined,remaining(),checkpoint));
            auto components = limits.exact_interior_components ?
                exact_section_interior_components(joined,remaining(),checkpoint) :
                exact_section_components(joined,remaining(),checkpoint);
            charge(components.work);
            for (auto& component : components.sections) normalized.emplace_back(std::move(component));
        } else {
            ExPolygons polygons;
            for (const auto& part : parts)
                polygons.insert(polygons.end(),part.native().begin(),part.native().end());
            for (auto& polygon : union_ex(polygons)) {
                validate(polygon);
                normalized.emplace_back(ExPolygons{std::move(polygon)});
            }
        }
        return normalized;
    };
    const auto overlap = [&](const SectionFootprint& a,const SectionFootprint& b) {
        if (exact_geometry)
            return SectionFootprint(exact_boolean(promote(a),promote(b),ExactSectionOperation::Intersection));
        return SectionFootprint(intersection_ex(a.native(),b.native()));
    };
    const auto contained = [&](const SectionFootprint& a,const SectionFootprint& b) {
        if (exact_geometry)
            return exact_boolean(promote(a),promote(b),ExactSectionOperation::Difference).empty();
        return diff_ex(a.native(),b.native()).empty();
    };
    struct Island { SectionFootprint footprint; BoundingBox bounds; std::size_t parent; bool bed; };
    std::vector<Island> islands;
    SectionTopology result;
    result.m_bed_z_m = double(bed) * 1e-9;
    const auto root = [&](std::size_t i) {
        while (islands[i].parent != i) {
            islands[i].parent = islands[islands[i].parent].parent;
            i = islands[i].parent;
        }
        return i;
    };
    const auto unite = [&](std::size_t a, std::size_t b) {
        a = root(a); b = root(b);
        if (a != b) {
            islands[b].parent = a;
            islands[a].bed = islands[a].bed || islands[b].bed;
        }
    };
    std::vector<std::size_t> previous;
    long long previous_top = bed;
    std::size_t island_points = 0, output_points = 0;
    for (auto &[interval, materials] : grouped) {
        if (checkpoint) checkpoint();
        const auto [bottom, top] = interval;
        if (bottom < previous_top)
            throw std::invalid_argument("Section topology requires aligned nonoverlapping Z intervals.");
        std::vector<SectionFootprint> occupied;
        std::size_t normalized_points = 0;
        for (auto &[material, polygons] : materials) {
            check_work();
            polygons = normalize(polygons);
            for (const auto &polygon : polygons) {
                normalized_points += polygon.vertex_count();
            }
            if (normalized_points > limits.max_boundary_points)
                throw ResourceLimit("Section topology normalized point budget exceeded.");
            if (!exact_geometry) {
                ExPolygons prior, current;
                for (const auto& part : occupied)
                    prior.insert(prior.end(),part.native().begin(),part.native().end());
                for (const auto& part : polygons)
                    current.insert(current.end(),part.native().begin(),part.native().end());
                const auto common = intersection_ex(prior,current);
                if (!common.empty()) {
                    double area = 0.;
                    for (const auto& part : common) area += part.area()*1e-18;
                    throw LayerMaterialOverlap(double(bottom)*1e-9,double(top)*1e-9,area,material);
                }
            } else for (const auto& prior : occupied) for (const auto& polygon : polygons) {
                check_work();
                const auto common = overlap(prior,polygon);
                if (common.empty()) continue;
                const auto bounds = common.bounds();
                if (!bounds) throw std::invalid_argument("Missing exact overlap bounds.");
                const auto measure = exact_section_volume(common.exact(),
                    {bounds->min.x(),bounds->min.y(),0},
                    {bounds->max.x()-bounds->min.x(),bounds->max.y()-bounds->min.y(),1000000000},
                    0,1000000000,remaining(),checkpoint);
                charge(measure.work);
                throw LayerMaterialOverlap(double(bottom)*1e-9, double(top)*1e-9, measure.value, material);
            }
            occupied.insert(occupied.end(), polygons.begin(), polygons.end());
        }
        check_work();
        occupied = normalize(occupied);
        std::vector<std::size_t> current;
        for (auto &polygon : occupied) {
            island_points += polygon.vertex_count();
            if (island_points > limits.max_boundary_points)
                throw ResourceLimit("Section topology island point budget exceeded.");
            const auto index = islands.size();
            const auto bounds = polygon.bounds();
            if (!bounds) throw std::invalid_argument("Missing section topology island bounds.");
            islands.push_back({std::move(polygon), *bounds, index, bottom == bed});
            current.push_back(index);
        }
        if (bottom == previous_top) {
            for (auto i : current) for (auto j : previous) {
                check_work();
                if (!islands[i].bounds.overlap(islands[j].bounds)) continue;
                const auto contact = overlap(islands[i].footprint,islands[j].footprint);
                if (contact.is_exact()) {
                    if (!contact.empty()) unite(i,j);
                } else {
                    double area = 0.;
                    for (const auto& polygon : contact.native()) area += polygon.area();
                    if (area > 0.) unite(i,j);
                }
            }
        }
        // Each normalized material polygon belongs wholly to one occupied
        // island. Do not divide it using coarse-grid cells or centroid tests.
        for (auto &[material, polygons] : materials) for (auto &polygon : polygons) {
            bool assigned = false;
            const auto bounds = polygon.bounds();
            if (!bounds) throw std::invalid_argument("Missing section topology material bounds.");
            for (auto i : current) {
                check_work();
                if (!bounds->overlap(islands[i].bounds)) continue;
                if (!contained(polygon,islands[i].footprint)) continue;
                output_points += polygon.vertex_count();
                if (output_points > limits.max_boundary_points)
                    throw ResourceLimit("Section topology output point budget exceeded.");
                result.m_pieces.push_back({{double(bottom)*1e-9, double(top)*1e-9, material,
                    std::move(polygon)}, i});
                assigned = true;
                break;
            }
            if (!assigned) throw std::invalid_argument("Material polygon has unresolved solid connectivity.");
        }
        previous = std::move(current);
        previous_top = top;
    }
    if (islands.empty()) throw std::invalid_argument("Section topology has no occupied solid.");
    std::map<std::size_t, std::size_t> components;
    for (auto &piece : result.m_pieces) {
        if (checkpoint) checkpoint();
        const auto id = root(piece.component);
        if (!islands[id].bed)
            throw std::invalid_argument("Section topology contains an unanchored solid component.");
        const auto inserted = components.emplace(id, components.size());
        piece.component = inserted.first->second;
    }
    result.m_components = components.size();
    return result;
}
} // namespace Slic3r::Thermoelastic
