#include "ThermoelasticLayers.hpp"
#include "ThermoelasticBoundary.hpp"
#include "ThermoelasticExactSections.hpp"
#include "BoundingBox.hpp"
#include "ClipperUtils.hpp"
#include "ExPolygon.hpp"
#include "ThermoelasticMesh.hpp"
#include "Point.hpp"
#include "Polygon.hpp"
#include "Polyline.hpp"
#include "libslic3r.h"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {

std::size_t qualify_complete_occupied_partition(const std::vector<const DepositedSlab*>& slabs,
    std::size_t max_work, const std::function<void()>& checkpoint)
{
    if (slabs.empty() || slabs.size()>256 || !slabs.front() || !max_work || max_work>4000000)
        throw std::invalid_argument("Invalid occupied partition qualification request.");
    std::size_t work=0;
    const auto charge = [&](std::size_t amount) {
        if (checkpoint) checkpoint();
        if (amount>max_work-work) throw ResourceLimit("Occupied partition qualification budget exceeded.");
        work+=amount;
    };
    charge(slabs.size()+1);
    const auto& proof=slabs.front()->occupied_partition;
    const double bottom=slabs.front()->bottom_m,top=slabs.front()->top_m;
    if (!proof || proof->regions.size()!=slabs.size() || !std::isfinite(bottom) || !std::isfinite(top) ||
        top<=bottom || std::abs(bottom)>10. || std::abs(top)>10. ||
        proof->sample_z_m<bottom || proof->sample_z_m>top)
        throw std::invalid_argument("Missing or incomplete occupied partition evidence.");
    std::set<std::size_t> materials;
    for (const auto* slab : slabs) {
        charge(1);
        if (!slab || slab->occupied_partition!=proof || slab->bottom_m!=bottom || slab->top_m!=top ||
            slab->boundary_uncertainty_nm!=proof->material_allowance_nm || !slab->footprint.is_exact() ||
            !materials.insert(slab->material).second)
            throw std::invalid_argument("Occupied partition members or sampling band changed.");
        const auto found=std::find_if(proof->regions.begin(),proof->regions.end(),[&](const auto& region) {
            charge(1);
            return region.material==slab->material;
        });
        if (found==proof->regions.end())
            throw std::invalid_argument("Occupied partition material changed.");
        if (work==max_work) throw ResourceLimit("Occupied partition qualification budget exhausted.");
        const auto changed=boolean_exact_sections(slab->footprint.exact(),found->footprint.exact(),
            ExactSectionOperation::SymmetricDifference,max_work-work,[&] { charge(0); });
        charge(changed.work);
        if (!changed.section.empty())
            throw std::invalid_argument("Occupied partition geometry changed.");
    }
    return work;
}

namespace {
// Keep integer operations well inside Clipper's range, even after offsets.
coord_t coordinate(double metres)
{
    if (!std::isfinite(metres) || std::abs(metres) > 10.)
        throw std::invalid_argument("Layered deposition coordinate is outside the supported range.");
    return coord_t(std::llround(metres * 1e9));
}
coord_t plane_key(double metres, const LayerGeometryLimits& limits)
{
    const auto native = coordinate(metres); // Same finite/range validation.
    return limits.half_native_section_planes ? coord_t(std::llround(metres * 2e9)) : native;
}
std::size_t points(const ExPolygons &polygons)
{
    std::size_t count = 0;
    for (const auto &polygon : polygons) {
        count += polygon.contour.points.size();
        for (const auto &hole : polygon.holes)
            count += hole.points.size();
    }
    return count;
}
struct Group {
    coord_t bottom, top, width;
    std::size_t material;
    Polylines paths;
    ExPolygons footprint;
};
// Cut an explicitly rounded bead whose top Z varies linearly along the path.
// Native-coordinate ellipse arcs bound the central strip. Original endpoint
// half-discs retain endpoint height; no new round caps appear at interior cuts.
ExPolygons rounded_sloped_section(const Point &start, const Point &end,
    coord_t z0, coord_t z1, coord_t height, coord_t width, double z,
    const LayerGeometryLimits &limits, const std::function<void()> &checkpoint)
{
    const Eigen::Vector2d a(double(start.x()), double(start.y()));
    const Eigen::Vector2d d(double(end.x()-start.x()), double(end.y()-start.y()));
    const double length = d.norm(), dz = double(z1-z0), tolerance = limits.arc_tolerance_m * 1e9;
    const Eigen::Vector2d direction = d / length, normal(-direction.y(), direction.x());
    const double flat_radius = .5 * double(width-height), half_height = .5 * double(height);
    const double pi = std::acos(-1.);
    std::size_t generated = 0;
    const auto count_steps = [&](double angle, double radius) {
        const double step = 2. * std::acos(std::clamp(1. - tolerance / radius, -1., 1.));
        if (!(step > 0.)) throw ResourceLimit("Rounded sloped section exceeds arc resolution.");
        const double count = std::max(1., std::ceil(std::abs(angle) / step));
        if (count > double(limits.max_boundary_points) || count > 100000.)
            throw ResourceLimit("Rounded sloped section exceeds arc point budget.");
        return std::size_t(count);
    };
    const auto append = [&](Polygon &polygon, const Eigen::Vector2d &p) {
        if (checkpoint) checkpoint();
        if (++generated > limits.max_boundary_points)
            throw ResourceLimit("Rounded sloped section exceeds boundary point budget.");
        const Point point(coord_t(std::llround(p.x())), coord_t(std::llround(p.y())));
        if (polygon.points.empty() || polygon.points.back() != point) polygon.points.push_back(point);
    };
    Polygons pieces;
    double lo = double(z-z0) / dz, hi = double(z+height-z0) / dz;
    if (lo > hi) std::swap(lo, hi);
    lo = std::max(0., lo);
    hi = std::min(1., hi);
    if (hi > lo) {
        const auto angle_at = [&](double t) {
            return std::asin(std::clamp((double(z0-z) + t*dz - half_height) / half_height, -1., 1.));
        };
        const double first = angle_at(lo), last = angle_at(hi);
        const double major_radius = std::max(half_height, length * half_height / std::abs(dz));
        const auto steps = count_steps(last-first, major_radius);
        Polygon strip;
        for (unsigned side = 0; side < 2; ++side)
            for (std::size_t i = 0; i <= steps; ++i) {
                const double fraction = double(side == 0 ? i : steps-i) / double(steps);
                const double angle = first + fraction * (last-first);
                const double depth = half_height * (1. + std::sin(angle));
                const double t = std::clamp((depth - double(z0-z)) / dz, lo, hi);
                const double radius = flat_radius + half_height * std::cos(angle);
                append(strip, a + t*d + (side == 0 ? radius : -radius)*normal);
            }
        if (strip.points.size() >= 3) {
            strip.make_counter_clockwise();
            pieces.push_back(std::move(strip));
        }
    }
    for (unsigned endpoint = 0; endpoint < 2; ++endpoint) {
        const coord_t top = endpoint ? z1 : z0;
        if (z < top-height || z >= top) continue;
        const double depth = double(top-z);
        const double radius = flat_radius + std::sqrt(std::max(0., depth * (double(height)-depth)));
        if (!(radius > 0.)) continue;
        const auto steps = count_steps(pi, radius);
        const Eigen::Vector2d center = a + (endpoint ? d : Eigen::Vector2d::Zero());
        Polygon cap;
        for (std::size_t i = 0; i <= steps; ++i) {
            const double angle = pi * double(i) / double(steps);
            append(cap, center + radius * (std::cos(angle)*normal +
                (endpoint ? 1. : -1.) * std::sin(angle)*direction));
        }
        if (cap.points.size() >= 3) {
            cap.make_counter_clockwise();
            pieces.push_back(std::move(cap));
        }
    }
    return union_ex(pieces);
}
// Variable design thickness (e.g. the two complementary scarf ramps). The
// rounded side radius is concave: a linear term plus the geometric mean of
// distance to the affine top and bottom planes. Locate its maximum departure
// from each chord, rather than relying on a midpoint-only sampling heuristic.
ExPolygons tapered_section(const Point &start, const Point &end, coord_t z0, coord_t z1,
    coord_t h0, coord_t h1, coord_t width, double z, BeadSection section,
    const LayerGeometryLimits &limits, const std::function<void()> &checkpoint)
{
    const Eigen::Vector2d a(double(start.x()), double(start.y()));
    const Eigen::Vector2d d(double(end.x()-start.x()), double(end.y()-start.y()));
    const Eigen::Vector2d direction = d.normalized(), normal(-direction.y(), direction.x());
    const double dz = double(z1-z0), dh = double(h1-h0), tolerance = limits.arc_tolerance_m*1e9;
    double lo = 0., hi = 1.;
    const auto clip = [&](double value, double slope, bool strict) {
        if (slope > 0.) lo = std::max(lo, -value/slope);
        else if (slope < 0.) hi = std::min(hi, -value/slope);
        else if (strict ? value <= 0. : value < 0.) hi = lo;
    };
    clip(double(z0-z), dz, true);
    clip(double(z-z0+h0), dh-dz, false);
    const auto radius = [&](double t) {
        if (section == BeadSection::RectangularEnvelope) return .5*double(width);
        const double upper = std::max(0., double(z0-z)+t*dz);
        const double lower = std::max(0., double(z-z0+h0)+t*(dh-dz));
        return .5*(double(width-h0)-t*dh) + std::sqrt(upper*lower);
    };
    std::size_t generated = 0;
    const auto append = [&](Polygon &polygon, const Eigen::Vector2d &p) {
        if (checkpoint) checkpoint();
        if (++generated > std::min<std::size_t>(limits.max_boundary_points, 100000))
            throw ResourceLimit("Tapered bead section exceeds the boundary-point budget.");
        const Point point(coord_t(std::llround(p.x())), coord_t(std::llround(p.y())));
        if (polygon.points.empty() || polygon.points.back() != point) polygon.points.push_back(point);
    };
    Polygons pieces;
    if (hi > lo) {
        std::vector<double> knots{lo}, pending{hi};
        while (!pending.empty()) {
            if (checkpoint) checkpoint();
            const double left = knots.back(), right = pending.back();
            double peak = .5*(left+right);
            if (section == BeadSection::RoundedRectangle) {
                const double slope = (radius(right)-radius(left))/(right-left);
                double l = left, r = right;
                for (unsigned i = 0; i < 64; ++i) {
                    peak = .5*(l+r);
                    const double upper = std::max(0., double(z0-z)+peak*dz);
                    const double lower = std::max(0., double(z-z0+h0)+peak*(dh-dz));
                    const double root = std::sqrt(upper*lower);
                    const double derivative = -.5*dh + (root > 0. ?
                        (dz*lower + (dh-dz)*upper)/(2.*root) : 0.);
                    if (derivative > slope) l = peak; else r = peak;
                }
            }
            const double chord = radius(left) + (peak-left)/(right-left)*(radius(right)-radius(left));
            if (radius(peak)-chord > tolerance) {
                if (!(peak > left && peak < right) || knots.size()+pending.size() >=
                    std::min<std::size_t>(limits.max_boundary_points/2, 50000))
                    throw ResourceLimit("Tapered bead section exceeds arc resolution.");
                pending.push_back(peak);
            } else {
                knots.push_back(right);
                pending.pop_back();
            }
        }
        Polygon strip;
        for (double t : knots) append(strip, a+t*d+radius(t)*normal);
        for (auto it = knots.rbegin(); it != knots.rend(); ++it)
            append(strip, a+*it*d-radius(*it)*normal);
        if (strip.points.size() >= 3) {
            strip.make_counter_clockwise();
            pieces.push_back(std::move(strip));
        }
    }
    for (unsigned endpoint = 0; endpoint < 2; ++endpoint) {
        const coord_t top = endpoint ? z1 : z0, height = endpoint ? h1 : h0;
        if (height == 0 || z < top-height || z >= top) continue;
        const double r = radius(double(endpoint));
        if (!(r > 0.)) continue;
        const double pi = std::acos(-1.);
        const double step = 2.*std::acos(std::clamp(1.-tolerance/r, -1., 1.));
        if (!(step > 0.)) throw ResourceLimit("Tapered bead cap exceeds arc resolution.");
        const double count = std::max(2., std::ceil(pi/step));
        if (count > double(std::min<std::size_t>(limits.max_boundary_points, 100000)))
            throw ResourceLimit("Tapered bead cap exceeds the point budget.");
        Polygon cap;
        for (std::size_t i = 0; i <= std::size_t(count); ++i) {
            const double angle = pi*double(i)/count;
            append(cap, a+double(endpoint)*d+r*(std::cos(angle)*normal+
                (endpoint ? 1. : -1.)*std::sin(angle)*direction));
        }
        cap.make_counter_clockwise();
        pieces.push_back(std::move(cap));
    }
    return union_ex(pieces);
}
void validate_limits(const std::vector<DepositedSegment> &segments, double bed_z_m,
    const LayerGeometryLimits &limits)
{
    if (segments.empty() || !std::isfinite(bed_z_m) ||
        !std::isfinite(limits.arc_tolerance_m) || limits.arc_tolerance_m < 1e-9 ||
        limits.arc_tolerance_m > 1e-4 || limits.max_segments == 0 ||
        limits.max_segments > 2000000 || limits.max_groups == 0 || limits.max_groups > 10000 ||
        limits.max_boundary_points == 0 || limits.max_boundary_points > 2000000 ||
        limits.max_slab_group_visits == 0 || limits.max_slab_group_visits > 8000000 || limits.ownership_boundary_grid_units > 4 ||
        (limits.retain_exact_interface_sections && limits.ownership_boundary_grid_units != 4) ||
        (limits.exact_per_command_composition && !limits.retain_exact_interface_sections))
        throw std::invalid_argument("Invalid layered deposition input or budget.");
    if (segments.size() > limits.max_segments)
        throw ResourceLimit("Layered deposition exceeds the segment budget.");
}
void validate_bead(const DepositedSegment &segment, coord_t bed, const LayerGeometryLimits &limits,
    bool allow_rounded = false)
{
    const double end_height = segment.end_height_m.value_or(segment.height_m);
    if (!segment.start_m.allFinite() || !segment.end_m.allFinite() ||
        !std::isfinite(segment.width_m) || segment.width_m <= 0. ||
        !std::isfinite(segment.height_m) || segment.height_m < 0. ||
        !std::isfinite(end_height) || end_height < 0. ||
        (segment.section != BeadSection::RectangularEnvelope &&
         !(allow_rounded && segment.section == BeadSection::RoundedRectangle)))
        throw std::invalid_argument("Layered deposition requires a supported finite bead section.");
    if (segment.section == BeadSection::RoundedRectangle && segment.width_m < std::max(segment.height_m, end_height))
        throw std::invalid_argument("Rounded bead width must be at least its height.");
    const coord_t bottom = std::min(coordinate(segment.start_m.z()) - coordinate(segment.height_m),
        coordinate(segment.end_m.z()) - coordinate(end_height));
    const Point start(coordinate(segment.start_m.x()), coordinate(segment.start_m.y()));
    const Point end(coordinate(segment.end_m.x()), coordinate(segment.end_m.y()));
    coordinate(segment.start_m.z());
    coordinate(segment.end_m.z());
    if (bottom < bed || coordinate(std::max(segment.height_m, end_height)) < 1 || coordinate(segment.width_m) < 2 || start == end)
        throw std::invalid_argument("Layered deposition has an unresolved bead dimension or bed intersection.");
    if (limits.arc_tolerance_m < segment.width_m * 1e-6)
        throw ResourceLimit("Layered deposition cap resolution exceeds the bounded geometry budget.");
}
} // namespace


DepositedLayers layered_deposition(const std::vector<DepositedSegment> &segments,
    double bed_z_m, const LayerGeometryLimits &limits, const std::function<void()> &checkpoint)
{
    validate_limits(segments, bed_z_m, limits);
    const coord_t bed = coordinate(bed_z_m);
    std::map<std::tuple<coord_t, coord_t, coord_t, std::size_t>, Group> groups;
    std::set<coord_t> levels;
    for (const auto &segment : segments) {
        if (checkpoint) checkpoint();
        validate_bead(segment, bed, limits);
        if (coordinate(segment.start_m.z()) != coordinate(segment.end_m.z()) ||
            coordinate(segment.end_height_m.value_or(segment.height_m)) != coordinate(segment.height_m))
            throw std::invalid_argument("Layered deposition requires horizontal rectangular-envelope beads.");
        const coord_t top = coordinate(segment.end_m.z());
        const coord_t bottom = coordinate(segment.end_m.z() - segment.height_m);
        const coord_t width = coordinate(segment.width_m);
        const Point start(coordinate(segment.start_m.x()), coordinate(segment.start_m.y()));
        const Point end(coordinate(segment.end_m.x()), coordinate(segment.end_m.y()));
        if (bottom < bed || top <= bottom || width < 2 || start == end)
            throw std::invalid_argument("Layered deposition has an unresolved bead dimension or bed intersection.");
        const auto key = std::make_tuple(bottom, top, width, segment.material);
        auto found = groups.find(key);
        if (found == groups.end()) {
            if (groups.size() >= limits.max_groups)
                throw ResourceLimit("Layered deposition exceeds the bead-group budget.");
            found = groups.emplace(key, Group{bottom, top, width, segment.material, {}, {}}).first;
        }
        auto &paths = found->second.paths;
        if (!paths.empty() && paths.back().last_point() == start)
            paths.back().append(end);
        else
            paths.emplace_back(start, end);
        levels.insert(bottom);
        levels.insert(top);
    }
    std::size_t boundary_points = 0;
    for (auto &entry : groups) {
        auto &group = entry.second;
        // Bound each Clipper call independently of the length of a print.
        // Splitting paths only repeats nominal round end caps; their union is
        // the same capsule envelope, with no added safety offset or gap fill.
        Polylines batch;
        std::size_t batch_segments = 0;
        const auto flush = [&] {
            if (batch.empty()) return;
            if (checkpoint) checkpoint();
            auto polygons = offset(batch, float(group.width / 2.), jtRound,
                limits.arc_tolerance_m * 1e9, etOpenRound);
            group.footprint = union_ex(group.footprint, polygons);
            if (points(group.footprint) + boundary_points > limits.max_boundary_points)
                throw ResourceLimit("Layered deposition exceeds the boundary-point budget.");
            batch.clear();
            batch_segments = 0;
        };
        for (const auto &path : group.paths) {
            for (std::size_t i = 1; i < path.points.size(); ++i) {
                if (!batch.empty() && batch.back().last_point() == path.points[i - 1])
                    batch.back().append(path.points[i]);
                else
                    batch.emplace_back(path.points[i - 1], path.points[i]);
                if (++batch_segments == 512) flush();
            }
        }
        flush();
        boundary_points += points(group.footprint);
        if (boundary_points > limits.max_boundary_points)
            throw ResourceLimit("Layered deposition exceeds the boundary-point budget.");
        Polylines().swap(group.paths);
    }
    DepositedLayers result;
    std::size_t visits = 0;
    for (auto lower = levels.begin(), upper = std::next(lower); upper != levels.end(); ++lower, ++upper) {
        if (checkpoint) checkpoint();
        std::map<std::size_t, ExPolygons> material_regions;
        for (const auto &[key, group] : groups) {
            if (++visits > limits.max_slab_group_visits)
                throw ResourceLimit("Layered deposition exceeds the slab-work budget.");
            if (group.bottom <= *lower && group.top >= *upper) {
                auto &regions = material_regions[group.material];
                regions.insert(regions.end(), group.footprint.begin(), group.footprint.end());
            }
        }
        ExPolygons occupied;
        for (auto &[material, regions] : material_regions) {
            regions = union_ex(regions);
            const auto overlap = intersection_ex(occupied, regions);
            if (!overlap.empty()) {
                double area_m2 = 0.;
                for (const auto &region : overlap) area_m2 += region.area() * 1e-18;
                throw LayerMaterialOverlap(double(*lower) * 1e-9, double(*upper) * 1e-9, area_m2, material);
            }
            occupied.insert(occupied.end(), regions.begin(), regions.end());
            result.boundary_points += points(regions);
            if (result.boundary_points > limits.max_boundary_points)
                throw ResourceLimit("Layered deposition exceeds the output boundary-point budget.");
            for (const auto &region : regions)
                result.occupied_volume_m3 += region.area() * 1e-18 * double(*upper - *lower) * 1e-9;
            result.slabs.push_back({double(*lower) * 1e-9, double(*upper) * 1e-9, material, std::move(regions)});
        }
    }
    return result;
}

DepositedSection section_deposition(const std::vector<DepositedSegment> &input_segments,
    double z_m, double bed_z_m, const LayerGeometryLimits &input_limits, const std::function<void()> &checkpoint,
    const SectionOwnership *input_ownership, const ExPolygon *input_clip, ExactOverlayDiagnostics* diagnostics)
{
    const bool frozen_deposited_inputs=input_ownership && bool(input_ownership->resolve_deposited_partition);
    std::vector<DepositedSegment> frozen_segments;
    std::vector<std::size_t> frozen_cut_points;
    LayerGeometryLimits frozen_limits;
    std::optional<ExPolygon> frozen_clip;
    std::optional<SectionOwnership> frozen_ownership;
    std::size_t copy_work=0;
    if (frozen_deposited_inputs) {
        // Freeze before ANY checkpoint, geometry, or callable copy. Evidence
        // and composition must see the same inputs even if a wrapper mutates
        // the caller-owned vectors/limits/clip while resolving this request.
        if (input_segments.size()>4096)
            throw ResourceLimit("Deposited sampler input snapshot command budget exceeded.");
        if (!input_ownership->regions.empty() || input_ownership->resolve_conflict ||
            input_ownership->resolve_sampled_conflict || input_ownership->resolve_partition ||
            input_ownership->partition || input_ownership->deposited_partition)
            throw std::invalid_argument("Deposited sampler requires one deferred ownership descriptor.");
        std::size_t clip_points=0;
        if (input_clip) {
            const auto count=[&](std::size_t amount) {
                if (amount>1024-clip_points)
                    throw ResourceLimit("Deposited sampler clipping snapshot point budget exceeded.");
                clip_points+=amount;
            };
            if (input_clip->holes.size()>1024)
                throw ResourceLimit("Deposited sampler clipping snapshot ring budget exceeded.");
            count(input_clip->contour.points.size());
            for (const auto& hole:input_clip->holes) count(hole.points.size());
        }
        frozen_limits=input_limits;
        copy_work=2*input_segments.size()+clip_points+(input_clip ? input_clip->holes.size() : 0)+1;
        if (!frozen_limits.max_slab_group_visits || frozen_limits.max_slab_group_visits>8000000)
            throw std::invalid_argument("Invalid deposited sampler snapshot work budget.");
        if (copy_work>frozen_limits.max_slab_group_visits)
            throw ResourceLimit("Deposited sampler snapshot work budget exceeded.");
        frozen_segments=input_segments;
        frozen_cut_points.assign(input_segments.size(),0);
        if (input_clip) frozen_clip=*input_clip;
        frozen_ownership=*input_ownership;
        if (checkpoint) checkpoint();
    }
    const auto& segments=frozen_deposited_inputs ? frozen_segments : input_segments;
    const auto& limits=frozen_deposited_inputs ? frozen_limits : input_limits;
    const auto* ownership=frozen_deposited_inputs ? &*frozen_ownership : input_ownership;
    const auto* clip=frozen_deposited_inputs ? (frozen_clip ? &*frozen_clip : nullptr) : input_clip;
    validate_limits(segments, bed_z_m, limits);
    if (clip && (!clip->is_valid() || clip->area() <= 0. || points(ExPolygons{*clip}) > limits.max_boundary_points))
        throw std::invalid_argument("Invalid deposition section clipping polygon.");
    const double z = double(plane_key(z_m, limits)) / (limits.half_native_section_planes ? 2. : 1.);
    const coord_t bed = coordinate(bed_z_m);
    if (z < bed) throw std::invalid_argument("Deposition section is below the bed.");
    const bool deposited_policy=ownership && bool(ownership->resolve_deposited_partition);
    if (limits.exact_per_command_composition && ownership)
        throw std::invalid_argument("Exact command composition does not establish material ownership.");
    const bool exact_composition=deposited_policy || limits.exact_per_command_composition;
    if (deposited_policy && !limits.retain_exact_interface_sections)
        throw std::invalid_argument("Deposited-order ownership requires retained exact sections.");
    struct MaterialRegion {
        ExPolygons footprint, pending;
        ExactSection exact;
        std::array<ExactSection,13> exact_runs;
        std::size_t segments = 0, pending_points = 0;
    };
    std::map<std::size_t, MaterialRegion> regions;
    std::size_t total_points = 0, ownership_work = copy_work, exact_commands = 0;
    SectionSamplerWork sampler_work;
    const auto record_work=[&](DepositedSection& section) {
        section.ownership_work=ownership_work;
        section.sampler_work=sampler_work;
        // Each category is a subset of successfully checked existing debits.
        section.sampler_work.other=ownership_work-sampler_work.exact_construction-sampler_work.exact_composition;
    };
    const auto exact_charge=[&](std::size_t amount,std::size_t* category=nullptr) {
        if (checkpoint) checkpoint();
        if (amount>limits.max_slab_group_visits-ownership_work)
            throw ResourceLimit("Exact deposited sampler work budget exceeded.");
        ownership_work+=amount;
        if (category) *category+=amount;
    };
    const auto exact_budget=[&] {
        if (ownership_work==limits.max_slab_group_visits)
            throw ResourceLimit("Exact deposited sampler work budget exhausted.");
        return std::min<std::size_t>(4000000,limits.max_slab_group_visits-ownership_work);
    };
    const auto merge_exact=[&](ExactSection& destination,ExactSection& source) {
        if (source.empty()) return;
        if (destination.empty()) { destination=std::move(source); source={}; return; }
        auto united=boolean_exact_sections(destination,source,ExactSectionOperation::Union,exact_budget(),checkpoint,diagnostics);
        exact_charge(united.work,&sampler_work.exact_composition);
        const auto retained=total_points-destination.vertex_count()-source.vertex_count();
        if (united.section.vertex_count()>std::min<std::size_t>(4096,limits.max_boundary_points)-retained)
            throw ResourceLimit("Exact deposited sampler retained vertex budget exceeded.");
        total_points=retained+united.section.vertex_count();
        destination=std::move(united.section);
        source={};
    };
    const auto collapse_exact=[&] {
        for (auto& entry:regions)
            for (auto& run:entry.second.exact_runs) merge_exact(entry.second.exact,run);
    };
    const auto flush = [&](MaterialRegion &region) {
        if (checkpoint) checkpoint();
        total_points -= points(region.footprint) + region.pending_points;
        region.footprint = union_ex(region.footprint, region.pending);
        region.pending.clear();
        region.segments = region.pending_points = 0;
        total_points += points(region.footprint);
        if (total_points > limits.max_boundary_points)
            throw ResourceLimit("Deposition section exceeds the boundary-point budget.");
    };
    auto cut_limits=limits;
    if (deposited_policy) cut_limits.max_boundary_points=std::min<std::size_t>(limits.max_boundary_points,1024);
    for (std::size_t bead_index=0;bead_index<segments.size();++bead_index) {
        const auto& bead=segments[bead_index];
        if (checkpoint) checkpoint();
        validate_bead(bead, bed, limits, true);
        const coord_t z0 = coordinate(bead.start_m.z()), z1 = coordinate(bead.end_m.z());
        const coord_t height = coordinate(bead.height_m);
        const coord_t end_height = coordinate(bead.end_height_m.value_or(bead.height_m));
        if (z < std::min(z0-height, z1-end_height) || z >= std::max(z0, z1)) continue;
        const Point start(coordinate(bead.start_m.x()), coordinate(bead.start_m.y()));
        const Point end(coordinate(bead.end_m.x()), coordinate(bead.end_m.y()));
        double radius = .5 * coordinate(bead.width_m);
        if (height == end_height && bead.section == BeadSection::RoundedRectangle && z0 == z1) {
            // The nominal section is a rectangle with semicircular sides,
            // not a full-width rectangle at every Z. Keep its corner voids.
            const double depth = double(z0 - z);
            radius = .5 * (coordinate(bead.width_m) - height) +
                std::sqrt(std::max(0., depth * (double(height) - depth)));
            if (radius == 0.) continue; // Tangent section of a circular bead.
        }
        auto footprint = height != end_height ?
            tapered_section(start, end, z0, z1, height, end_height, coordinate(bead.width_m), z,
                bead.section, cut_limits, checkpoint) :
            bead.section == BeadSection::RoundedRectangle && z0 != z1 ?
            rounded_sloped_section(start, end, z0, z1, height, coordinate(bead.width_m), z, cut_limits, checkpoint) :
            union_ex(offset(Polylines{Polyline(start, end)}, float(radius), jtRound,
                limits.arc_tolerance_m * 1e9, etOpenRound));
        if (height == end_height && z0 != z1 && bead.section == BeadSection::RectangularEnvelope) {
            // Within the central strip, bottom <= z < top selects an interval
            // of the path parameter. Endpoint caps retain endpoint height;
            // they are not replaced by round caps at the cut plane.
            double lo = double(z - z0) / double(z1 - z0);
            double hi = double(z + height - z0) / double(z1 - z0);
            if (lo > hi) std::swap(lo, hi);
            const Eigen::Vector2d a(double(start.x()), double(start.y()));
            const Eigen::Vector2d d(double(end.x() - start.x()), double(end.y() - start.y()));
            const double length = d.norm();
            const Eigen::Vector2d n(-d.y() / length, d.x() / length);
            // Extend beyond the original cap only when its endpoint lies in
            // the section; otherwise cut across the straight capsule strip.
            const double margin = radius + 2.;
            lo = z >= z0-height && z < z0 ? -margin / length : std::max(0., lo);
            hi = z >= z1-height && z < z1 ? 1. + margin / length : std::min(1., hi);
            const Eigen::Vector2d left = a + lo * d, right = a + hi * d;
            const auto point = [](const Eigen::Vector2d &p) {
                return Point(coord_t(std::llround(p.x())), coord_t(std::llround(p.y())));
            };
            const Polygon clip{point(left - margin*n), point(right - margin*n),
                point(right + margin*n), point(left + margin*n)};
            footprint = intersection_ex(footprint, Polygons{clip});
        }
        if (clip) footprint = intersection_ex(footprint, ExPolygons{*clip});
        if (footprint.empty()) continue;
        if (regions.find(bead.material) == regions.end() && regions.size() >= limits.max_groups)
            throw ResourceLimit("Deposition section exceeds the material-group budget.");
        auto &region = regions[bead.material];
        const auto count = points(footprint);
        if (exact_composition) {
            if (deposited_policy && count>cut_limits.max_boundary_points)
                throw ResourceLimit("Deposited command boundary point budget exceeded.");
            if (limits.exact_per_command_composition && regions.size()>1)
                throw std::invalid_argument("Exact command composition requires a single material.");
            // Preserve the same quantized individual bead cut and clip, but
            // compose its material union exactly. A cross-bead integer union
            // could otherwise round away intersections before authentication.
            exact_charge(count+1,&sampler_work.exact_construction);
            if (++exact_commands>4096)
                throw ResourceLimit("Exact deposited sampler command budget exceeded.");
            // Match the singleton sampler's final flush as well as its cut:
            // the factory independently reconstructs exactly that singleton.
            const auto normalized=union_ex(ExPolygons{},footprint);
            const auto normalized_count=points(normalized);
            if (deposited_policy) {
                if (normalized_count>cut_limits.max_boundary_points)
                    throw ResourceLimit("Deposited command boundary point budget exceeded.");
                frozen_cut_points[bead_index]=normalized_count;
            }
            // Stable retained payload counts every run (including overlapping
            // aliases), across all materials. Normalize and exact overlay
            // temporaries use their existing bounded operation workspace.
            // Consolidate before constructing/admitting a cut under pressure.
            const auto ceiling=std::min<std::size_t>(4096,limits.max_boundary_points);
            if (normalized_count>ceiling-total_points) collapse_exact();
            auto cut=make_exact_section(normalized,exact_budget(),checkpoint);
            exact_charge(cut.work,&sampler_work.exact_construction);
            // Cuts are consumed one at a time, not retained as an input batch.
            // Bound the current material unions plus this cut; cumulative
            // construction and overlay work still uses the same work budget.
            if (cut.section.vertex_count()>1024 || cut.section.vertex_count()>ceiling-total_points)
                throw ResourceLimit("Exact deposited sampler live input vertex budget exceeded.");
            total_points+=cut.section.vertex_count();
            // Equal-sized command runs merge in a binary carry tree instead
            // of overlaying every cut against a growing material prefix.
            auto carry=std::move(cut.section);
            for (auto& run:region.exact_runs) {
                if (run.empty()) {
                    run=std::move(carry); carry={};
                    break;
                }
                merge_exact(carry,run);
            }
            if (!carry.empty()) throw ResourceLimit("Exact deposited sampler run budget exceeded.");
            continue;
        }
        total_points += count;
        region.pending_points += count;
        region.pending.insert(region.pending.end(),
            std::make_move_iterator(footprint.begin()), std::make_move_iterator(footprint.end()));
        if (total_points > limits.max_boundary_points)
            throw ResourceLimit("Deposition section exceeds the boundary-point budget.");
        if (++region.segments == 512) flush(region);
    }
    if (exact_composition) collapse_exact();
    DepositedSection result;
    if (limits.exact_per_command_composition) {
        for (const auto& [material,region]:regions) {
            if (region.exact.empty()) continue;
            if (const auto bounds=region.exact.bounds()) {
                const std::array<coord_t,3> start{bounds->min.x(),bounds->min.y(),0};
                const std::array<coord_t,3> size{bounds->max.x()-bounds->min.x(),bounds->max.y()-bounds->min.y(),1};
                const auto area=exact_section_horizontal_area(region.exact,start,size,0,exact_budget(),checkpoint);
                exact_charge(area.work);
                result.area_m2+=area.value;
            }
            result.regions.push_back({material,SectionFootprint(region.exact)});
        }
        result.boundary_points=total_points;
        record_work(result);
        return result;
    }
    ExPolygons occupied, conflicts;
    std::map<std::size_t, ExPolygons> authenticated_claims;
    bool authenticated_partition = false;
    const auto charge = [&](const ExPolygons& a, const ExPolygons& b) {
        if (checkpoint) checkpoint();
        const auto count = points(a) + points(b) + 1;
        if (count > limits.max_slab_group_visits - ownership_work)
            throw ResourceLimit("Section ownership exceeds its polygon-work budget.");
        ownership_work += count;
    };
    const auto bounded = [&](ExPolygons value) {
        if (points(value) > limits.max_boundary_points)
            throw ResourceLimit("Section ownership exceeds its boundary-point budget.");
        return value;
    };
    const auto exact_empty = [&](const ExPolygons& a, const ExPolygons& b, ExactSectionOperation operation) {
        const auto query = [&](const ExPolygons& left, const ExPolygons& right) {
            if (ownership_work == limits.max_slab_group_visits)
                throw ResourceLimit("Exact ownership predicate work budget exhausted.");
            try {
                const auto result = exact_section_empty(left, right, operation, 2048,
                    std::min<size_t>(4000000, limits.max_slab_group_visits - ownership_work), checkpoint);
                ownership_work += result.work;
                return result.empty;
            } catch (const ResourceLimit& error) {
                throw ResourceLimit("Exact ownership predicate operands: vertices=" + std::to_string(points(left)) + "," +
                    std::to_string(points(right)) + " components=" + std::to_string(left.size()) + "," +
                    std::to_string(right.size()) + ": " + error.what());
            }
        };
        if (operation != ExactSectionOperation::Intersection) return query(a, b);
        // Intersection emptiness distributes over component pairs. Difference
        // and union equality do NOT, and retain the full-set query above.
        std::vector<BoundingBox> boxes;
        boxes.reserve(b.size());
        for (const auto& polygon : b) boxes.push_back(get_extents(polygon));
        for (const auto& p : a) {
            const auto box = get_extents(p);
            for (std::size_t i = 0; i < b.size(); ++i) {
                charge({}, {});
                const auto& other = boxes[i];
                if (std::min(box.max.x(), other.max.x()) <= std::max(box.min.x(), other.min.x()) ||
                    std::min(box.max.y(), other.max.y()) <= std::max(box.min.y(), other.min.y())) continue;
                // Holes wholly outside the other operand's outer bounds
                // cannot affect their intersection. Remove only those from
                // the predicate copy, not from the deposited geometry.
                const auto relevant_holes = [&](const ExPolygon& polygon,const BoundingBox& bounds) {
                    ExPolygon copy=polygon;
                    copy.holes.erase(std::remove_if(copy.holes.begin(),copy.holes.end(),[&](const Polygon& hole) {
                        const auto h=get_extents(hole);
                        return h.max.x()<bounds.min.x() || h.min.x()>bounds.max.x() ||
                            h.max.y()<bounds.min.y() || h.min.y()>bounds.max.y();
                    }),copy.holes.end());
                    return copy;
                };
                charge(ExPolygons{p},ExPolygons{b[i]});
                const auto left=relevant_holes(p,other),right=relevant_holes(b[i],box);
                if (count_points(left) + count_points(right) > 2048)
                    throw ResourceLimit("Exact ownership interacting components exceed vertex budget: " +
                        std::to_string(count_points(left)) + "," + std::to_string(count_points(right)));
                if (!query(ExPolygons{left}, ExPolygons{right})) return false;
            }
        }
        return true;
    };
    if (!deposited_policy) for (auto &[material, region] : regions) {
        flush(region);
        auto &footprint = region.footprint;
        if (ownership) charge(occupied, footprint);
        const auto overlap = intersection_ex(occupied, footprint);
        if (!overlap.empty()) {
            double area = 0.;
            for (const auto &region : overlap) area += region.area() * 1e-18;
            if (!ownership) {
                LayerMaterialOverlap error(z_m, z_m, area, material);
                error.section_footprint_complete = true;
                std::size_t retained_points = 0;
                for (const auto &polygon : overlap) {
                    std::size_t count = polygon.contour.points.size();
                    for (const auto &hole : polygon.holes) count += hole.points.size();
                    if (count > 4096 - retained_points) { error.section_footprint_complete = false; continue; }
                    retained_points += count;
                    error.section_footprint.push_back(polygon);
                }
                throw error;
            }
            charge(conflicts, overlap);
            conflicts = bounded(union_ex(conflicts, overlap));
        }
        occupied.insert(occupied.end(), footprint.begin(), footprint.end());
    }
    if (ownership) {
        if (plane_key(ownership->z_m, limits) != plane_key(z_m, limits))
            throw std::invalid_argument("Section ownership belongs to a different plane.");
        if (ownership->partition || ownership->deposited_partition)
            throw std::invalid_argument("A section partition must be produced for the live sampler request.");
        if ((ownership->resolve_conflict || ownership->resolve_sampled_conflict || ownership->resolve_partition ||
             ownership->resolve_deposited_partition) && !ownership->regions.empty())
            throw std::invalid_argument("Section ownership cannot combine masks with a deferred resolver.");
        if (unsigned(bool(ownership->resolve_conflict)) + unsigned(bool(ownership->resolve_sampled_conflict)) +
            unsigned(bool(ownership->resolve_partition)) + unsigned(bool(ownership->resolve_deposited_partition)) > 1)
            throw std::invalid_argument("Section ownership cannot combine different deferred resolvers.");
        if (ownership->resolve_deposited_partition) {
            if (!limits.retain_exact_interface_sections)
                throw std::invalid_argument("Deposited-order ownership requires retained exact sections.");
            LayerMaterialOverlap complete(z_m,z_m,0.,0);
            complete.section_footprint=conflicts;
            complete.section_footprint_complete=true;
            std::vector<DepositedSection::Region> original;
            for (const auto& [material,region]:regions) original.push_back({material,SectionFootprint(region.exact)});
            const SampledSectionConflict evidence(complete,segments,bed_z_m,limits,clip,original,
                frozen_deposited_inputs,&frozen_cut_points);
            const auto token=ownership->resolve_deposited_partition(evidence);
            if (!token || token->identity!=evidence.m_identity || token->m_policy!=DepositedOverlapPolicy::FirstDepositionWins)
                throw std::invalid_argument("Deposited partition does not match its live sampler request.");
            if (token->work>limits.max_slab_group_visits-ownership_work)
                throw ResourceLimit("Deposited partition work budget exceeded.");
            ownership_work+=token->work;
            std::set<std::size_t> materials;
            std::size_t retained=0;
            for (const auto& replacement:token->regions) {
                if (!regions.count(replacement.material) || !materials.insert(replacement.material).second ||
                    !replacement.footprint.is_exact())
                    throw std::invalid_argument("Invalid deposited partition material or representation.");
                if (replacement.footprint.vertex_count()>limits.max_boundary_points-retained)
                    throw ResourceLimit("Deposited partition output budget exceeded.");
                retained+=replacement.footprint.vertex_count();
                if (const auto bounds=replacement.footprint.bounds()) {
                    if (ownership_work==limits.max_slab_group_visits)
                        throw ResourceLimit("Deposited partition area work budget exhausted.");
                    const std::array<coord_t,3> start{bounds->min.x(),bounds->min.y(),0};
                    const std::array<coord_t,3> size{bounds->max.x()-bounds->min.x(),bounds->max.y()-bounds->min.y(),1};
                    const auto area=exact_section_horizontal_area(replacement.footprint.exact(),start,size,0,
                        std::min<std::size_t>(4000000,limits.max_slab_group_visits-ownership_work),checkpoint);
                    ownership_work+=area.work;
                    result.area_m2+=area.value;
                }
            }
            // Sum exact nonnegative per-material exclusions. Subtracting
            // rounded total areas could erase a tiny overlap by cancellation;
            // overlap UNION area would undercount triple-material coverage.
            double excluded_nominal_area=0.;
            for (const auto& entry:regions) {
                exact_charge(1);
                ExactSection removed=entry.second.exact;
                for (const auto& replacement:token->regions) {
                    exact_charge(1);
                    if (replacement.material!=entry.first) continue;
                    const auto difference=boolean_exact_sections(removed,replacement.footprint.exact(),
                        ExactSectionOperation::Difference,exact_budget(),checkpoint);
                    exact_charge(difference.work);
                    removed=difference.section;
                    break;
                }
                if (const auto bounds=removed.bounds()) {
                    const std::array<coord_t,3> start{bounds->min.x(),bounds->min.y(),0};
                    const std::array<coord_t,3> size{bounds->max.x()-bounds->min.x(),bounds->max.y()-bounds->min.y(),1};
                    const auto area=exact_section_horizontal_area(removed,start,size,0,exact_budget(),checkpoint);
                    exact_charge(area.work);
                    if (!std::isfinite(area.value) || !(area.value>0.))
                        throw std::invalid_argument("Invalid deposited excluded nominal area.");
                    excluded_nominal_area+=area.value;
                }
            }
            result.regions=token->regions;
            result.boundary_points=retained;
            record_work(result);
            result.resolved_overlap_area_m2=token->m_overlap_area_m2;
            result.deposited_partition=token;
            // Complete occupied geometry evidence only. The separate policy
            // tag is not native-boundary or physical material qualification.
            auto occupied=std::shared_ptr<OccupiedSectionPartition>(new OccupiedSectionPartition);
            occupied->regions=result.regions;
            occupied->sample_z_m=z_m;
            occupied->deposited_policy=token->m_policy;
            occupied->deposited_overlap_area=token->m_overlap_area_m2;
            occupied->deposited_excluded_nominal_area=excluded_nominal_area;
            result.occupied_partition=std::move(occupied);
            return result;
        }
        if (ownership->resolve_partition && !conflicts.empty()) {
            LayerMaterialOverlap complete(z_m, z_m, 0., 0);
            complete.section_footprint = conflicts;
            complete.section_footprint_complete = true;
            std::vector<DepositedSection::Region> original;
            for (const auto& [material, region] : regions) original.push_back({material, region.footprint});
            const SampledSectionConflict evidence(complete, segments, bed_z_m, limits, clip, original);
            const auto token = ownership->resolve_partition(evidence);
            if (!token || token->identity != evidence.m_identity || token->regions.size() != regions.size() ||
                limits.ownership_boundary_grid_units != 4 || token->interface_allowance_nm != 7 ||
                token->representation_allowance_nm != (limits.retain_exact_interface_sections ? 0u : 1u))
                throw std::invalid_argument("Section partition does not match this authenticated sampler request.");
            if (token->work > limits.max_slab_group_visits - ownership_work)
                throw ResourceLimit("Section partition work budget exceeded.");
            ownership_work += token->work;
            std::size_t retained = 0;
            std::set<std::size_t> replacement_materials;
            for (const auto& replacement : token->regions) {
                const auto found = regions.find(replacement.material);
                if (found == regions.end() || !replacement_materials.insert(replacement.material).second)
                    throw std::invalid_argument("Section partition material is absent or duplicated.");
                if (replacement.footprint.is_exact() != limits.retain_exact_interface_sections)
                    throw std::invalid_argument("Section partition representation does not match its sampler.");
                const auto count = replacement.footprint.vertex_count();
                if (count > limits.max_boundary_points - retained)
                    throw ResourceLimit("Section partition output budget exceeded.");
                retained += count;
                if (!limits.retain_exact_interface_sections)
                    found->second.footprint = replacement.footprint.native();
            }
            authenticated_partition = true;
            result.boundary_error.uncertainty_nm = token->interface_allowance_nm;
            for (const auto& polygon : conflicts) result.resolved_overlap_area_m2 += polygon.area() * 1e-18;
            if (limits.retain_exact_interface_sections) {
                // The authenticated exact producer proved disjoint material
                // sets and equality to these complete sampled envelopes.
                // Retain those sets unchanged; no integer union or difference
                // may substitute a rounded footprint during publication.
                for (const auto& replacement : token->regions) {
                    const auto bounds=replacement.footprint.bounds();
                    if (bounds) {
                        if (ownership_work==limits.max_slab_group_visits)
                            throw ResourceLimit("Exact section area work budget exhausted.");
                        const std::array<coord_t,3> start{bounds->min.x(),bounds->min.y(),0};
                        const std::array<coord_t,3> size{bounds->max.x()-bounds->min.x(),bounds->max.y()-bounds->min.y(),1};
                        const auto area=exact_section_horizontal_area(replacement.footprint.exact(),start,size,0,
                            std::min<std::size_t>(4000000,limits.max_slab_group_visits-ownership_work),checkpoint);
                        ownership_work+=area.work;
                        result.area_m2+=area.value;
                    }
                }
                result.regions=token->regions;
                result.boundary_points=retained;
                record_work(result);
                auto occupied=std::shared_ptr<OccupiedSectionPartition>(new OccupiedSectionPartition);
                occupied->regions=result.regions;
                occupied->sample_z_m=z_m;
                occupied->material_allowance_nm=token->interface_allowance_nm;
                occupied->dual_native_boundary_area_m2=token->dual_native_boundary_area_m2;
                result.occupied_partition=std::move(occupied);
                return result;
            }
        }
    }
    if (ownership && !authenticated_partition) {
        std::vector<DepositedSection::Region> deferred;
        if ((ownership->resolve_conflict || ownership->resolve_sampled_conflict) && !conflicts.empty()) {
            double area = 0.;
            for (const auto& polygon : conflicts) area += polygon.area() * 1e-18;
            LayerMaterialOverlap complete(z_m, z_m, area, 0);
            complete.section_footprint = conflicts;
            complete.section_footprint_complete = true;
            if (ownership->resolve_sampled_conflict) {
                const std::vector<DepositedSection::Region> unused;
                const SampledSectionConflict evidence(complete, segments, bed_z_m, limits, clip, unused);
                deferred = ownership->resolve_sampled_conflict(evidence);
            } else deferred = ownership->resolve_conflict(complete);
        }
        const auto& owner_regions = ownership->resolve_conflict || ownership->resolve_sampled_conflict ? deferred : ownership->regions;
        if (owner_regions.size() > limits.max_groups)
            throw ResourceLimit("Section ownership exceeds its material-group budget.");
        std::map<std::size_t, ExPolygons> masks;
        std::size_t mask_points = 0, merged_mask_points = 0;
        for (const auto& owner : owner_regions) {
            if (checkpoint) checkpoint();
            const auto count = points(owner.footprint.native());
            if (count > limits.max_boundary_points - mask_points)
                throw ResourceLimit("Section ownership exceeds its input-point budget.");
            mask_points += count;
            for (const auto& polygon : owner.footprint.native()) {
                const auto validate = [](const Polygon& ring) {
                    if (ring.points.size() < 3) throw std::invalid_argument("Invalid section ownership polygon.");
                    for (const auto& p : ring.points)
                        if (std::abs(double(p.x())) > 1e10 || std::abs(double(p.y())) > 1e10)
                            throw std::invalid_argument("Section ownership coordinate is outside the supported range.");
                };
                validate(polygon.contour);
                for (const auto& hole : polygon.holes) validate(hole);
            }
            auto& mask = masks[owner.material];
            charge(mask, owner.footprint.native());
            merged_mask_points -= points(mask);
            mask = bounded(union_ex(mask, owner.footprint.native()));
            merged_mask_points += points(mask);
            if (merged_mask_points > limits.max_boundary_points)
                throw ResourceLimit("Section ownership exceeds its merged-mask budget.");
        }
        if (!conflicts.empty()) {
            auto& claims = authenticated_claims;
            ExPolygons assigned;
            std::size_t claim_points = 0;
            for (const auto& [material, mask] : masks) {
                charge(conflicts, mask);
                auto claim = bounded(intersection_ex(conflicts, mask));
                if (claim.empty()) continue;
                charge(assigned, claim);
                const auto ambiguous = intersection_ex(assigned, claim);
                if (!ambiguous.empty() &&
                    (!limits.ownership_boundary_grid_units ||
                     !exact_empty(assigned, claim, ExactSectionOperation::Intersection))) {
                    double area_nm2 = 0.;
                    for (const auto& polygon : ambiguous) area_nm2 += polygon.area();
                    throw std::invalid_argument("Section ownership is ambiguous within a material conflict (z_nm=" +
                        std::to_string(z_m * 1e9) + ", material=" + std::to_string(material) +
                        ", vertices=" + std::to_string(points(assigned)) + "," + std::to_string(points(claim)) +
                        ", components=" + std::to_string(assigned.size()) + "," + std::to_string(claim.size()) +
                        ", overlap_nm2=" + std::to_string(area_nm2) + ").");
                }
                const auto found = regions.find(material);
                if (found == regions.end())
                    throw std::invalid_argument("Section owner has no deposited envelope.");
                charge(claim, found->second.footprint);
                const auto missing = diff_ex(claim, found->second.footprint);
                if (!missing.empty() && (!limits.ownership_boundary_grid_units ||
                    !exact_empty(claim, found->second.footprint, ExactSectionOperation::Difference))) {
                    double area = 0.;
                    for (const auto& polygon : missing) area += polygon.area();
                    throw std::invalid_argument("Section owner has no deposited envelope at the conflict (material=" +
                        std::to_string(material) + ", z_nm=" + std::to_string(z_m * 1e9) +
                        ", missing_nm2=" + std::to_string(area) + ").");
                }
                charge(assigned, claim);
                assigned = bounded(union_ex(assigned, claim));
                claim_points += points(claim);
                if (claim_points > limits.max_boundary_points)
                    throw ResourceLimit("Section ownership exceeds its claim-point budget.");
                claims.emplace(material, std::move(claim));
            }
            charge(conflicts, assigned);
            if (!diff_ex(conflicts, assigned).empty() && (!limits.ownership_boundary_grid_units ||
                !exact_empty(conflicts, assigned, ExactSectionOperation::Difference)))
                throw std::invalid_argument("Section ownership does not cover the material conflict.");
            std::size_t retained_points = 0;
            for (auto& [material, region] : regions) {
                charge(region.footprint, conflicts);
                auto retained = bounded(diff_ex(region.footprint, conflicts));
                const auto claim = claims.find(material);
                if (claim != claims.end()) {
                    charge(retained, claim->second);
                    retained = bounded(union_ex(retained, claim->second));
                }
                retained_points += points(retained);
                if (retained_points > limits.max_boundary_points)
                    throw ResourceLimit("Section ownership exceeds its retained-point budget.");
                region.footprint = std::move(retained);
            }
            for (const auto& polygon : conflicts)
                result.resolved_overlap_area_m2 += polygon.area() * 1e-18;
        }
    }
    ExPolygons resolved;
    for (auto& [material, region] : regions) {
        bool exactly_disjoint = false;
        if (!authenticated_partition && ownership && !conflicts.empty() && limits.ownership_boundary_grid_units && !resolved.empty()) {
            // All macroscopic ownership was established above. This ordered
            // subtraction only normalizes residual integer-grid intersections;
            // a complete boundary/topology check rejects any larger change.
            const auto original = region.footprint;
            for (unsigned attempt = 0; attempt < 4; ++attempt) {
                charge(resolved, region.footprint);
                if (intersection_ex(resolved, region.footprint).empty()) break;
                // Rounded diagnostic rings are not proof of actual overlap.
                // Query the original polygons before altering either material.
                if (exact_empty(resolved, region.footprint, ExactSectionOperation::Intersection)) {
                    exactly_disjoint = true;
                    break;
                }
                region.footprint = bounded(diff_ex(region.footprint, resolved));
            }
            charge(original, region.footprint);
            const auto removed = bounded(diff_ex(original, region.footprint));
            const auto added = bounded(diff_ex(region.footprint, original));
            ExPolygons foreign_claims;
            for (const auto& [owner, claim] : authenticated_claims)
                if (owner < material) foreign_claims.insert(foreign_claims.end(), claim.begin(), claim.end());
            charge(removed, foreign_claims);
            // Restoring an already-authenticated foreign claim is intentional
            // partitioning, not an exterior-boundary representation error.
            // The complete deposited union is still checked independently.
            const bool authorized_removal = added.empty() && diff_ex(removed, foreign_claims).empty();
            if ((!removed.empty() || !added.empty()) && !authorized_removal) {
                charge(original, foreign_claims);
                const auto intended = bounded(diff_ex(original, foreign_claims));
                if (ownership_work == limits.max_slab_group_visits)
                    throw ResourceLimit("Section ownership boundary work budget exhausted.");
                SectionBoundaryError error;
                try {
                    error = qualify_section_boundary_error(intended, region.footprint,
                        limits.ownership_boundary_grid_units, limits.max_boundary_points,
                        limits.max_slab_group_visits - ownership_work, checkpoint);
                } catch (const std::invalid_argument& failure) {
                    throw std::invalid_argument("Material partition representation (material=" + std::to_string(material) + "): " + failure.what());
                }
                ownership_work += error.work;
                result.boundary_error.uncertainty_nm = std::max(result.boundary_error.uncertainty_nm, error.uncertainty_nm);
            }
        }
        result.boundary_points += points(region.footprint);
        if (result.boundary_points > limits.max_boundary_points)
            throw ResourceLimit("Section ownership exceeds its output-point budget.");
        for (const auto& polygon : region.footprint) result.area_m2 += polygon.area() * 1e-18;
        if (ownership && (!conflicts.empty() || limits.ownership_boundary_grid_units)) {
            // A rounded empty intersection is not a proof of disjointness:
            // an exact sub-grid overlap may have disappeared from the initial
            // diagnostic. Never send that unresolved material pair downstream.
            if (limits.ownership_boundary_grid_units && !resolved.empty() &&
                !exact_empty(resolved,region.footprint,ExactSectionOperation::Intersection))
                throw std::invalid_argument("Unresolved exact material overlap in the section output at z_m=" +
                    std::to_string(z_m) + ", material=" + std::to_string(material) +
                    ", authenticated_partition=" + (authenticated_partition ? "true" : "false") +
                    ", initial_conflict=" + (conflicts.empty() ? "false" : "true") + ".");
            charge(resolved, region.footprint);
            const auto residual_overlap = intersection_ex(resolved, region.footprint);
            if (!residual_overlap.empty() && !exactly_disjoint &&
                (!limits.ownership_boundary_grid_units || !exact_empty(resolved, region.footprint, ExactSectionOperation::Intersection))) {
                std::string diagnostic;
                const auto& ring = residual_overlap.front().contour.points;
                for (std::size_t i = 0; i < std::min<std::size_t>(ring.size(), 12); ++i)
                    diagnostic += " " + std::to_string(ring[i].x()) + "," + std::to_string(ring[i].y());
                throw std::invalid_argument("Resolved section materials are not strictly disjoint; first native ring:" + diagnostic);
            }
            // Let the final Boolean check union the complete polygon set once,
            // This matches the original occupied representation above.
            resolved.insert(resolved.end(), region.footprint.begin(), region.footprint.end());
        }
        result.regions.push_back({material, std::move(region.footprint)});
    }
    if (ownership && !conflicts.empty()) {
        charge(resolved, occupied);
        const auto added = diff_ex(resolved, occupied), removed = diff_ex(occupied, resolved);
        if (authenticated_partition || ((!added.empty() || !removed.empty()) &&
            (!limits.ownership_boundary_grid_units || !exact_empty(occupied, resolved, ExactSectionOperation::SymmetricDifference)))) {
            if (limits.ownership_boundary_grid_units) {
                if (ownership_work == limits.max_slab_group_visits)
                    throw ResourceLimit("Section ownership boundary work budget exhausted.");
                std::vector<const ExPolygons*> material_regions;
                for (const auto& region : result.regions) material_regions.push_back(&region.footprint.native());
                SectionBoundaryError error;
                try {
                    error = qualify_section_boundary_error(occupied, resolved,
                        limits.ownership_boundary_grid_units, limits.max_boundary_points,
                        limits.max_slab_group_visits - ownership_work, checkpoint, material_regions);
                } catch (const std::invalid_argument& failure) {
                    std::string location;
                    if (clip) {
                        const auto box = get_extents(*clip);
                        location = ", clip_nm=" + std::to_string(box.min.x()) + "," + std::to_string(box.min.y()) +
                            ":" + std::to_string(box.max.x()) + "," + std::to_string(box.max.y());
                    }
                    throw std::invalid_argument("Occupied union representation at z_m=" + std::to_string(z_m) + location + ": " + failure.what());
                }
                ownership_work += error.work;
                result.boundary_error.uncertainty_nm = std::max(result.boundary_error.uncertainty_nm, error.uncertainty_nm);
                result.boundary_error.added_area_bound_m2 = error.added_area_bound_m2;
                result.boundary_error.removed_area_bound_m2 = error.removed_area_bound_m2;
            } else {
                double added_area = 0., removed_area = 0.;
                for (const auto& polygon : added) added_area += polygon.area();
                for (const auto& polygon : removed) removed_area += polygon.area();
                throw std::invalid_argument("Section ownership changed the deposited union (native square-nanometres: added=" +
                    std::to_string(added_area) + ", removed=" + std::to_string(removed_area) + ").");
            }
        }
    }
    record_work(result);
    return result;
}

struct DepositionSectionCursor::Impl
{
    struct Interval
    {
        coord_t bottom, top;
        std::size_t segment;
    };
    const std::vector<DepositedSegment>& segments;
    std::vector<double> z_m;
    double bed_z_m;
    SectionSweepLimits limits;
    std::function<void()> checkpoint;
    std::function<std::optional<SectionOwnership>(std::size_t, double)> ownership;
    const ExPolygon* clip;
    std::vector<coord_t> planes;
    std::vector<Interval> intervals;
    std::set<std::size_t> active;
    std::multimap<coord_t, std::size_t> ending;
    std::vector<DepositedSegment> selected;
    std::size_t next = 0, visits = 0, plane = 0;
    bool busy = false, failed = false;
    Impl(const std::vector<DepositedSegment>& source,
         std::vector<double> schedule,
         double bed,
         const SectionSweepLimits& bounds,
         std::function<void()> check,
         std::function<std::optional<SectionOwnership>(std::size_t, double)> owner,
         const ExPolygon* clipping)
        : segments(source)
        , z_m(std::move(schedule))
        , bed_z_m(bed)
        , limits(bounds)
        , checkpoint(std::move(check))
        , ownership(std::move(owner))
        , clip(clipping)
    {
        validate_limits(segments, bed_z_m, limits.geometry);
        if (z_m.empty() || limits.max_sections == 0 || limits.max_sections > 10000 || limits.max_segment_visits == 0 ||
            limits.max_segment_visits > 20000000)
            throw std::invalid_argument("Invalid deposition section sweep or budget.");
        if (z_m.size() > limits.max_sections)
            throw ResourceLimit("Deposition section sweep exceeds the plane budget.");
        const coord_t bed_native  = coordinate(bed_z_m);
        const coord_t plane_scale = limits.geometry.half_native_section_planes ? 2 : 1;
        planes.reserve(z_m.size());
        for (const double z : z_m) {
            const coord_t plane = plane_key(z, limits.geometry);
            if (plane < plane_scale * bed_native || (!planes.empty() && plane <= planes.back()))
                throw std::invalid_argument("Deposition section planes must increase at the selected plane resolution.");
            planes.push_back(plane);
        }
        intervals.reserve(segments.size());
        for (std::size_t i = 0; i < segments.size(); ++i) {
            if (checkpoint && i % 1024 == 0)
                checkpoint();
            const auto& bead = segments[i];
            validate_bead(bead, bed_native, limits.geometry, true);
            const coord_t start = coordinate(bead.start_m.z()), end = coordinate(bead.end_m.z());
            // Match the single-plane inclusion test exactly, including its integer
            // height rounding. Do not infer regular layers from nozzle Z values.
            intervals.push_back(
                {plane_scale * std::min(start - coordinate(bead.height_m), end - coordinate(bead.end_height_m.value_or(bead.height_m))),
                 plane_scale * std::max(start, end), i});
        }
        std::sort(intervals.begin(), intervals.end(),
                  [](const auto& a, const auto& b) { return std::tie(a.bottom, a.segment) < std::tie(b.bottom, b.segment); });
    }
    void step(const std::function<void(std::size_t, const DepositedSection&)>& visit)
    {
        if (checkpoint)
            checkpoint();
        const coord_t z = planes[plane];
        while (!ending.empty() && ending.begin()->first <= z) {
            active.erase(ending.begin()->second);
            ending.erase(ending.begin());
        }
        while (next < intervals.size() && intervals[next].bottom <= z) {
            if (checkpoint && next % 1024 == 0)
                checkpoint();
            const auto& interval = intervals[next++];
            if (interval.top <= z)
                continue;
            active.insert(interval.segment);
            ending.emplace(interval.top, interval.segment);
        }
        if (active.size() > limits.max_segment_visits - visits)
            throw ResourceLimit("Deposition section sweep exceeds the bead-visit budget.");
        visits += active.size();
        selected.clear();
        selected.reserve(active.size());
        // Preserve input order independently of sweep activation order, so
        // bounded polygon-union batches match the single-plane implementation.
        for (const auto index : active)
            selected.push_back(segments[index]);
        const auto owner   = ownership ? ownership(plane, z_m[plane]) : std::nullopt;
        const auto section = selected.empty() ? DepositedSection{} :
                                                section_deposition(selected, z_m[plane], bed_z_m, limits.geometry, checkpoint,
                                                                   owner ? &*owner : nullptr, clip);
        visit(plane, section);
        ++plane;
    }
};
DepositionSectionCursor::DepositionSectionCursor(const std::vector<DepositedSegment>& segments,
                                                 std::vector<double> planes,
                                                 double bed,
                                                 const SectionSweepLimits& limits,
                                                 std::function<void()> checkpoint,
                                                 std::function<std::optional<SectionOwnership>(std::size_t, double)> ownership,
                                                 const ExPolygon* clip)
    : m_impl(std::make_unique<Impl>(segments, std::move(planes), bed, limits, std::move(checkpoint), std::move(ownership), clip))
{}
DepositionSectionCursor::~DepositionSectionCursor() = default;
bool DepositionSectionCursor::advance(const std::function<void(std::size_t, const DepositedSection&)>& visit)
{
    if (!visit || m_impl->failed || m_impl->busy) {
        m_impl->failed = true;
        throw std::logic_error("Inactive or reentered deposition section cursor.");
    }
    if (m_impl->plane == m_impl->planes.size())
        return false;
    m_impl->busy = true;
    try {
        m_impl->step(visit);
        if (m_impl->failed)
            throw std::logic_error("Reentered deposition section cursor.");
        m_impl->busy = false;
        return true;
    } catch (...) {
        m_impl->busy   = false;
        m_impl->failed = true;
        throw;
    }
}
std::size_t DepositionSectionCursor::completed_visits() const
{
    if (m_impl->failed || m_impl->busy || m_impl->plane != m_impl->planes.size()) {
        m_impl->failed = true;
        throw std::logic_error("Deposition section cursor is incomplete or inactive.");
    }
    return m_impl->visits;
}
std::size_t visit_deposition_sections(const std::vector<DepositedSegment>& segments,
                                      const std::vector<double>& z_m,
                                      double bed_z_m,
                                      const std::function<void(std::size_t, const DepositedSection&)>& visit,
                                      const SectionSweepLimits& limits,
                                      const std::function<void()>& checkpoint,
                                      const std::function<std::optional<SectionOwnership>(std::size_t, double)>& ownership,
                                      const ExPolygon* clip)
{
    if (!visit)
        throw std::invalid_argument("Missing deposition section consumer.");
    DepositionSectionCursor cursor(segments, z_m, bed_z_m, limits, checkpoint, ownership, clip);
    while (cursor.advance(visit)) {}
    return cursor.completed_visits();
}
} // namespace Slic3r::Thermoelastic
