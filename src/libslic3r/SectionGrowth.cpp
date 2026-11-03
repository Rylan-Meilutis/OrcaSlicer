#include "SectionGrowth.hpp"
#include "SectionGrowthStream.hpp"
#include "ThermoelasticLocalSections.hpp"
#include "ThermoelasticLocalTopology.hpp"
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticExactSections.hpp"
#include "ThermoelasticSectionFootprint.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include "ThermoelasticMesh.hpp"
#include "BoundingBox.hpp"
#include "ClipperUtils.hpp"
#include "ExPolygon.hpp"
#include "Polygon.hpp"
#include "libslic3r.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <array>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
namespace {
std::size_t points(const ExPolygons& polygons)
{
    std::size_t result = 0;
    for (const auto& polygon : polygons) {
        result += polygon.contour.points.size();
        for (const auto& hole : polygon.holes)
            result += hole.points.size();
    }
    return result;
}
coord_t height(double value) { return std::llround(value * 1e9); }
} // namespace

// Internal partition arithmetic, shared by retained and streamed growth.
// This does not produce a connectivity certificate.
namespace GrowthDetail {
// Index construction scales with indexed segments, unlike shared section
// planning, where each cursor traverses the same sampling intervals.
std::vector<std::size_t> index_work_quotas(std::size_t budget,const std::vector<std::size_t>& inputs)
{
    if (!budget || budget>500000000 || inputs.empty() || inputs.size()>258)
        throw std::invalid_argument("Invalid reference index allocation.");
    std::uint64_t total=0;
    for (const auto count : inputs) {
        if (!count || count>4000000 || total+count>8000000)
            throw std::invalid_argument("Invalid reference index input counts.");
        total+=count;
    }
    std::vector<std::size_t> result;
    for (const auto count : inputs) {
        // The validated factors fit uint64_t on every supported platform.
        const auto quota=std::uint64_t(budget)*count/total;
        if (!quota) throw ResourceLimit("Reference activation index allocation exhausted.");
        result.push_back(std::size_t(quota));
    }
    return result;
}
void subtract_partition(const std::vector<SectionPiece>& source,
                        const std::vector<SectionPiece>& cover,
                        bool old,
                        std::size_t point_limit,
                        const std::function<void(std::size_t)>& charge,
                        const std::function<std::size_t()>& remaining_work,
                        const std::function<void(std::size_t, std::size_t)>& contact,
                        const std::function<void(const SectionPiece&)>& added)
{
    const auto available = [&] {
        const auto budget = std::min<std::size_t>(remaining_work(), 4000000);
        if (!budget) throw ResourceLimit("Section growth exact work budget exhausted.");
        return budget;
    };
    const auto promote = [&](const SectionFootprint& footprint) {
        if (footprint.is_exact()) return footprint.exact();
        auto result = make_exact_section(footprint.native(), available(), [&] { charge(0); });
        charge(result.work);
        return std::move(result.section);
    };
    const auto combine = [&](const SectionFootprint& a, const SectionFootprint& b, ExactSectionOperation operation) {
        if (!a.is_exact() && !b.is_exact())
            return SectionFootprint(operation == ExactSectionOperation::Intersection ?
                intersection_ex(a.native(), b.native()) : diff_ex(a.native(), b.native()));
        const auto left = promote(a), right = promote(b);
        auto result = boolean_exact_sections(left, right, operation, available(), [&] { charge(0); });
        charge(result.work);
        return SectionFootprint(std::move(result.section));
    };
    // Both inputs are complete, non-overlapping Z partitions. Index each
    // cover band once so fine layer sampling does not cause an all-history
    // scan for every prism in a tile.
    struct Band
    {
        coord_t bottom, top;
        std::vector<const SectionPiece*> pieces;
    };
    std::map<std::pair<coord_t, coord_t>, std::vector<const SectionPiece*>> intervals;
    for (const auto& p : cover) {
        charge(1);
        intervals[{height(p.slab.bottom_m), height(p.slab.top_m)}].push_back(&p);
    }
    std::vector<Band> bands;
    for (auto& [interval, pieces] : intervals) {
        if (!bands.empty() && interval.first < bands.back().top)
            throw std::invalid_argument("Overlapping growth cover intervals.");
        bands.push_back({interval.first, interval.second, std::move(pieces)});
    }
    for (const auto& piece : source) {
        const auto& slab = piece.slab;
        if (slab.footprint.vertex_count() > point_limit)
            throw ResourceLimit("Section growth input point budget exceeded.");
        const coord_t bottom = height(slab.bottom_m), top = height(slab.top_m);
        charge(slab.footprint.vertex_count());
        const auto bounds = slab.footprint.bounds();
        if (slab.footprint.empty()) continue;
        if (!bounds) throw std::invalid_argument("Missing section growth footprint bounds.");
        std::vector<const SectionPiece*> candidates;
        std::vector<coord_t> planes{bottom, top};
        const auto first = std::lower_bound(bands.begin(), bands.end(), bottom, [&](const Band& band, coord_t z) {
            charge(1);
            return band.top <= z;
        });
        for (auto band = first; band != bands.end() && band->bottom < top; ++band)
            for (const auto* candidate : band->pieces) {
                const auto& other = *candidate;
                charge(1);
                const coord_t lo = std::max(bottom, height(other.slab.bottom_m));
                const coord_t hi = std::min(top, height(other.slab.top_m));
                if (lo >= hi)
                    continue;
                charge(other.slab.footprint.vertex_count());
                const auto other_bounds = other.slab.footprint.bounds();
                if (!other_bounds || !bounds->overlap(*other_bounds))
                    continue;
                candidates.push_back(&other);
                planes.push_back(lo);
                planes.push_back(hi);
            }
        std::sort(planes.begin(), planes.end());
        planes.erase(std::unique(planes.begin(), planes.end()), planes.end());
        for (std::size_t z = 1; z < planes.size(); ++z) {
            SectionFootprint remaining = slab.footprint;
            for (const auto* other : candidates) {
                charge(1);
                if (height(other->slab.bottom_m) > planes[z - 1] || height(other->slab.top_m) < planes[z])
                    continue;
                charge(remaining.vertex_count() + other->slab.footprint.vertex_count());
                const auto overlap = combine(remaining, other->slab.footprint, ExactSectionOperation::Intersection);
                if (overlap.vertex_count() > point_limit)
                    throw ResourceLimit("Section growth intersection budget exceeded.");
                if (overlap.empty())
                    continue;
                if (old) {
                    if (slab.material != other->slab.material)
                        throw std::invalid_argument("Existing section material changed during deposition.");
                    contact(piece.component, other->component);
                }
                charge(remaining.vertex_count() + other->slab.footprint.vertex_count());
                remaining = combine(remaining, other->slab.footprint, ExactSectionOperation::Difference);
                if (remaining.vertex_count() > point_limit)
                    throw ResourceLimit("Section growth difference budget exceeded.");
                if (remaining.empty())
                    break;
            }
            if (old && !remaining.empty())
                throw std::invalid_argument("Existing section geometry disappeared during deposition.");
            if (!old && !remaining.empty()) {
                added({{double(planes[z - 1]) * 1e-9, double(planes[z]) * 1e-9, slab.material, std::move(remaining),
                        slab.boundary_uncertainty_nm},
                       piece.component});
            }
        }
    }
}
} // namespace GrowthDetail

SectionGrowth qualify_section_growth(const SectionTopology& previous,
                                     const SectionTopology& next,
                                     const SectionGrowthLimits& limits,
                                     const std::function<void()>& checkpoint)
{
    if (previous.bed_z_m() != next.bed_z_m() || !limits.max_work || limits.max_work > 200000000 || !limits.max_points ||
        limits.max_points > 2000000 || !limits.max_pieces || limits.max_pieces > 10000)
        throw std::invalid_argument("Invalid section growth domain or limits.");
    SectionGrowth result;
    result.m_bed       = next.bed_z_m();
    const auto missing = std::numeric_limits<std::size_t>::max();
    result.m_components.assign(previous.components(), missing);
    std::size_t work = 0, output_points = 0;
    const auto charge = [&](std::size_t count) {
        if (checkpoint)
            checkpoint();
        if (count > limits.max_work - work)
            throw ResourceLimit("Section growth work budget exceeded.");
        work += count;
    };
    const auto contact = [&](std::size_t old, std::size_t next) {
        auto& component = result.m_components.at(old);
        if (component != missing && component != next)
            throw std::invalid_argument("Existing section component split during deposition.");
        component = next;
    };
    const auto added = [&](const SectionPiece& piece) {
        const auto count = piece.slab.footprint.vertex_count();
        if (count > limits.max_points - output_points || result.m_added.size() >= limits.max_pieces)
            throw ResourceLimit("Section growth output budget exceeded.");
        output_points += count;
        result.m_added.push_back(piece);
    };
    const auto remaining = [&] { return limits.max_work - work; };
    GrowthDetail::subtract_partition(previous.pieces(), next.pieces(), true, limits.max_points, charge, remaining, contact, added);
    if (std::find(result.m_components.begin(), result.m_components.end(), missing) != result.m_components.end())
        throw std::invalid_argument("Missing previous section component coverage.");
    GrowthDetail::subtract_partition(next.pieces(), previous.pieces(), false, limits.max_points, charge, remaining, contact, added);
    return result;
}

static double clip_reference_sections(const std::vector<const DepositedSlab*>& reference,
    const std::vector<const DepositedSlab*>& deposited, const SectionGrowthLimits& limits,
    const std::function<void(std::size_t)>& charge, const std::function<std::size_t()>& remaining_work,
    const std::function<void(DepositedSlab)>& emit,
    std::optional<std::pair<coord_t,coord_t>> window = {})
{
    std::size_t output_points = 0, output_pieces = 0;
    bool has_exact=false;
    std::set<coord_t> planes;
    const auto index = [&](const auto& slabs) {
        for (const auto* p : slabs) {
            charge(1 + p->footprint.vertex_count());
            has_exact=has_exact || p->footprint.is_exact();
            if (!window) {
                planes.insert(height(p->bottom_m));
                planes.insert(height(p->top_m));
            }
            if (planes.size() > limits.max_pieces + 1)
                throw ResourceLimit("Reference activation plane budget exceeded.");
        }
    };
    index(reference);
    index(deposited);
    if (window) {
        if (window->first >= window->second) throw std::invalid_argument("Invalid reference clipping band.");
        planes.insert(window->first);
        planes.insert(window->second);
    }
    const std::vector<coord_t> ordered(planes.begin(),planes.end());
    double excluded = 0.;
    for (std::size_t z = 1; z < ordered.size(); ++z) {
        const auto in_band = [&](const DepositedSlab* p) {
            return height(p->bottom_m)<=ordered[z-1] && height(p->top_m)>=ordered[z];
        };
        bool exact_band=false;
        if (has_exact) {
            for (const auto* p : reference) { charge(1); if (in_band(p) && p->footprint.is_exact()) exact_band=true; }
            for (const auto* p : deposited) { charge(1); if (in_band(p) && p->footprint.is_exact()) exact_band=true; }
        }
        if (exact_band) {
            const auto budget = [&] {
                const auto available=std::min<std::size_t>(4000000,remaining_work());
                if (!available) throw ResourceLimit("Exact reference activation budget exhausted.");
                return available;
            };
            std::vector<const DepositedSlab*> reference_band;
            bool preserve_partition=false;
            for (const auto* p : reference) if (in_band(p)) {
                reference_band.push_back(p);
                preserve_partition=preserve_partition || p->boundary_uncertainty_nm;
            }
            std::map<std::shared_ptr<const OccupiedSectionPartition>,std::vector<const DepositedSlab*>> source_proofs;
            for (const auto* p : deposited) if (in_band(p) && p->boundary_uncertainty_nm) {
                preserve_partition=true;
                source_proofs[p->occupied_partition].push_back(p);
            }
            if (preserve_partition) {
                // Only a complete unchanged reference can transport its scoped
                // evidence. A subset must obtain separate contact qualification.
                charge(qualify_complete_occupied_partition(reference_band,budget(),[&] { charge(0); }));
                for (const auto& [proof,members] : source_proofs)
                    charge(qualify_complete_occupied_partition(members,budget(),[&] { charge(0); }));
            }
            std::map<std::size_t,ExactSection> sources,regions;
            std::size_t live=0,input=0;
            const auto check_live = [&](std::size_t n) {
                if (n>limits.max_points-live) throw ResourceLimit("Exact reference activation live point budget exceeded.");
            };
            const auto collect_exact = [&](const auto& slabs,auto& target) {
                for (const auto* p : slabs) {
                    charge(1);
                    if (!in_band(p)) continue;
                    const auto n=p->footprint.vertex_count();
                    if (n>limits.max_points-input) throw ResourceLimit("Exact reference activation input point budget exceeded.");
                    input+=n;
                    charge(n);
                    ExactSection section;
                    if (p->footprint.is_exact()) section=p->footprint.exact();
                    else {
                        const auto promoted=make_exact_section(p->footprint.native(),budget(),[&] { charge(0); });
                        charge(promoted.work);
                        section=promoted.section;
                    }
                    check_live(section.vertex_count());
                    // Input, previous union and new union coexist through the
                    // Boolean result handoff, even when handles share storage.
                    live+=section.vertex_count();
                    auto& accumulated=target[p->material];
                    const auto merged=boolean_exact_sections(accumulated,section,ExactSectionOperation::Union,budget(),[&] { charge(0); });
                    charge(merged.work);
                    check_live(merged.section.vertex_count());
                    live-=accumulated.vertex_count();
                    accumulated=merged.section;
                    live+=accumulated.vertex_count();
                    live-=section.vertex_count();
                }
            };
            collect_exact(reference,regions);
            collect_exact(deposited,sources);
            std::vector<DepositedSlab> active;
            for (const auto& [material,source] : sources) {
                const auto kept=boolean_exact_sections(source,regions[material],ExactSectionOperation::Intersection,
                    budget(),[&] { charge(0); });
                charge(kept.work);
                check_live(kept.section.vertex_count());
                live+=kept.section.vertex_count();
                const auto removed=boolean_exact_sections(source,regions[material],ExactSectionOperation::Difference,
                    budget(),[&] { charge(0); });
                charge(removed.work);
                check_live(removed.section.vertex_count());
                if (!removed.section.empty()) {
                    const auto box=*removed.section.bounds();
                    const std::array<coord_t,3> base{box.min.x(),box.min.y(),ordered[z-1]};
                    const std::array<coord_t,3> size{box.max.x()-box.min.x(),box.max.y()-box.min.y(),ordered[z]-ordered[z-1]};
                    const auto volume=exact_section_volume(removed.section,base,size,ordered[z-1],ordered[z],
                        budget(),[&] { charge(0); });
                    charge(volume.work);
                    if (!(volume.value>0.) || !std::isfinite(volume.value))
                        throw std::invalid_argument("Unresolved exact excluded reference volume.");
                    excluded+=volume.value;
                }
                if (preserve_partition) {
                    live-=kept.section.vertex_count();
                    continue;
                }
                if (kept.section.empty()) continue;
                if (active.size()>=limits.max_pieces) throw ResourceLimit("Exact reference activation pending piece budget exceeded.");
                active.push_back({double(ordered[z-1])*1e-9,double(ordered[z])*1e-9,material,SectionFootprint(kept.section),0});
            }
            if (preserve_partition) {
                for (const auto* original : reference_band) {
                    const auto kept=boolean_exact_sections(regions.at(original->material),sources[original->material],
                        ExactSectionOperation::Intersection,budget(),[&] { charge(0); });
                    charge(kept.work);
                    check_live(kept.section.vertex_count());
                    live+=kept.section.vertex_count();
                    if (active.size()>=limits.max_pieces) throw ResourceLimit("Exact reference activation pending piece budget exceeded.");
                    auto candidate=*original;
                    candidate.bottom_m=double(ordered[z-1])*1e-9;
                    candidate.top_m=double(ordered[z])*1e-9;
                    candidate.footprint=SectionFootprint(kept.section);
                    // Diagnose before the qualification call: failed exact
                    // work must not be retried without a shared-budget debit.
                    const auto changed=boolean_exact_sections(original->footprint.exact(),kept.section,
                        ExactSectionOperation::SymmetricDifference,budget(),[&] { charge(0); });
                    charge(changed.work);
                    check_live(changed.section.vertex_count());
                    if (!changed.section.empty()) {
                        live+=changed.section.vertex_count();
                        const auto area = [&](const ExactSection& section) {
                            if (section.empty()) return 0.;
                            const auto box=*section.bounds();
                            const std::array<coord_t,3> base{box.min.x(),box.min.y(),0};
                            const std::array<coord_t,3> size{box.max.x()-box.min.x(),box.max.y()-box.min.y(),1};
                            const auto measure=exact_section_horizontal_area(section,base,size,0,
                                budget(),[&] { charge(0); });
                            charge(measure.work);
                            return measure.bounds[1];
                        };
                        const double reference_area=area(original->footprint.exact());
                        const double candidate_area=area(kept.section);
                        const double difference_area=area(changed.section);
                        const auto box=*changed.section.bounds();
                        std::ostringstream message;
                        message<<std::setprecision(17)<<"Occupied partition geometry changed during exact activation: material="
                            <<original->material<<", reference_band_nm=["<<height(original->bottom_m)<<','<<height(original->top_m)
                            <<"], candidate_band_nm=["<<ordered[z-1]<<','<<ordered[z]<<"], reference_area_upper_m2="
                            <<reference_area<<", candidate_area_upper_m2="<<candidate_area<<", difference_area_upper_m2="
                            <<difference_area<<", difference_bounds_nm=["<<box.min.x()<<','<<box.min.y()<<','
                            <<box.max.x()<<','<<box.max.y()<<"], source_bands_nm=";
                        std::size_t shown=0;
                        for (const auto* slab : deposited) {
                            charge(1);
                            if (!in_band(slab) || slab->material!=original->material) continue;
                            if (shown++==16) { message<<"..."; break; }
                            message<<'['<<height(slab->bottom_m)<<','<<height(slab->top_m)<<']';
                        }
                        throw std::invalid_argument(message.str());
                    }
                    // Retain empty members too: completeness is part of the proof.
                    active.push_back(std::move(candidate));
                }
                std::vector<const DepositedSlab*> candidates;
                for (const auto& candidate : active) candidates.push_back(&candidate);
                charge(qualify_complete_occupied_partition(candidates,budget(),[&] { charge(0); }));
            }
            for (auto& slab : active) {
                const auto n=slab.footprint.vertex_count();
                if (n>limits.max_points-output_points || output_pieces>=limits.max_pieces)
                    throw ResourceLimit("Reference activation output budget exceeded.");
                output_points+=n;
                ++output_pieces;
                emit(std::move(slab));
            }
            continue;
        }
        std::map<std::size_t,ExPolygons> source_regions, reference_regions;
        std::size_t live_points = 0;
        unsigned uncertainty = 0;
        const auto collect = [&](const auto& slabs, auto& regions) {
            for (const auto* p : slabs) {
                charge(1);
                if (height(p->bottom_m) > ordered[z-1] || height(p->top_m) < ordered[z]) continue;
                uncertainty = std::max(uncertainty,p->boundary_uncertainty_nm);
                const auto n = points(p->footprint.native());
                charge(n);
                if (n > limits.max_points - live_points)
                    throw ResourceLimit("Reference activation live point budget exceeded.");
                live_points += n;
                auto& polygons = regions[p->material];
                polygons.insert(polygons.end(),p->footprint.native().begin(),p->footprint.native().end());
            }
        };
        collect(reference,reference_regions);
        collect(deposited,source_regions);
        std::vector<DepositedSlab> active;
        std::size_t active_points=0;
        for (auto& [material, polygons] : source_regions) {
            charge(points(polygons));
            auto source = polygons.size() == 1 ? polygons : union_ex(polygons);
            if (points(source) > limits.max_points)
                throw ResourceLimit("Reference activation union point budget exceeded.");
            const auto& region = reference_regions[material];
            charge(points(source) + points(region));
            auto footprint = intersection_ex(source,region);
            charge(points(source) + points(region));
            const auto outside = diff_ex(source,region);
            if (points(outside) > limits.max_points)
                throw ResourceLimit("Reference activation difference point budget exceeded.");
            for (const auto& polygon : outside)
                excluded += polygon.area() * 1e-18 * double(ordered[z]-ordered[z-1]) * 1e-9;
            if (footprint.empty()) continue;
            const auto retained=points(footprint);
            if (retained>limits.max_points-active_points || active.size()>=limits.max_pieces)
                throw ResourceLimit("Reference activation pending output budget exceeded.");
            active_points+=retained;
            active.push_back({double(ordered[z-1])*1e-9,double(ordered[z])*1e-9,material,std::move(footprint),uncertainty});
        }
        std::set<std::size_t> ambiguous;
        const auto erased_island = [&](const ExPolygons& polygons, unsigned band) {
            if (!band) return false;
            for (const auto& polygon : polygons) {
                std::size_t input_points = polygon.contour.points.size();
                for (const auto& hole : polygon.holes) input_points += hole.points.size();
                charge(input_points + 1);
                const auto inner = offset_ex(ExPolygons{polygon}, -float(band + 2), jtRound, .25);
                const auto count = points(inner);
                charge(count + 1);
                if (count > limits.max_points)
                    throw ResourceLimit("Reference activation erosion diagnostic point budget exceeded.");
                if (inner.empty()) return true;
            }
            return false;
        };
        const auto weak_boundary = [&](const ExPolygons& polygons) {
            for (const auto& polygon : polygons) {
                std::set<std::pair<coord_t,coord_t>> vertices;
                const auto repeated = [&](const Polygon& ring) {
                    for (const auto& point : ring.points) {
                        charge(1);
                        if (!vertices.emplace(point.x(),point.y()).second) return true;
                    }
                    return false;
                };
                if (repeated(polygon.contour)) return true;
                for (const auto& hole : polygon.holes)
                    if (repeated(hole)) return true;
            }
            return false;
        };
        bool uncertain_island = false;
        for (const auto& slab : active)
            if (slab.boundary_uncertainty_nm && (weak_boundary(slab.footprint.native()) ||
                erased_island(slab.footprint.native(), slab.boundary_uncertainty_nm))) { uncertain_island = true; break; }
        if (!uncertain_island && uncertainty &&
            (active.size()>1 || (!active.empty() && active.front().footprint.native().size()>1))) {
            // Topology also unions materials before assigning occupied islands.
            // A rounded union can expose a point-touch loop or disconnect a
            // sliver even when every material polygon is individually simple.
            ExPolygons occupied;
            charge(active_points + 1);
            for (const auto& slab : active)
                occupied.insert(occupied.end(),slab.footprint.native().begin(),slab.footprint.native().end());
            auto joined = occupied.size()==1 ? std::move(occupied) : union_ex(occupied);
            const auto count = points(joined);
            charge(count + 1);
            if (count>limits.max_points) throw ResourceLimit("Reference activation occupied diagnostic point budget exceeded.");
            uncertain_island = weak_boundary(joined) || erased_island(joined,uncertainty);
        }
        for (std::size_t i=0;!uncertain_island && i<active.size();++i)
            for (std::size_t j=0;j<i;++j) {
                charge(points(active[i].footprint.native())+points(active[j].footprint.native()));
                if (!get_extents(active[i].footprint.native()).overlap(get_extents(active[j].footprint.native()))) continue;
                // An integer intersection may round a positive rational sliver
                // to empty. Only an exact predicate can rule out an ownership
                // conflict before the joint reference-clipping fallback.
                const auto available=std::min<std::size_t>(remaining_work(),2000000);
                if (!available) throw ResourceLimit("Reference activation overlap budget exhausted.");
                const auto overlap=exact_section_empty(active[i].footprint.native(),active[j].footprint.native(),
                    ExactSectionOperation::Intersection,1024,available,[&] { charge(0); },true);
                charge(overlap.work);
                if (!overlap.empty) {
                    ambiguous.insert(active[i].material);
                    ambiguous.insert(active[j].material);
                }
            }
        if (!ambiguous.empty() || uncertain_island) {
            // Independent integer clipping can move just one side of a shared
            // interface or create a tiny disconnected output component.
            // Reconstruct intersections from the original
            // operands and round their shared junctions together. Never trim
            // an overlap by arbitrarily choosing a material priority.
            std::vector<std::size_t> materials;
            for (const auto& entry : source_regions)
                if (!entry.second.empty() && !reference_regions[entry.first].empty()) materials.push_back(entry.first);
            if (materials.empty() || materials.size()>2)
                throw std::invalid_argument("Joint reference activation clipping requires at most two active materials.");
            std::array<ExPolygons,2> regions, sources;
            for (std::size_t i=0; i<materials.size(); ++i) {
                regions[i] = reference_regions.at(materials[i]);
                sources[i] = source_regions.at(materials[i]);
            }
            auto clipped=exact_clip_reference_sections(regions,sources,4000000,[&] { charge(0); },charge);
            active.clear();
            active_points = 0;
            for (std::size_t i=0;i<materials.size();++i) {
                if (clipped.sections[i].empty()) continue;
                const auto count = points(clipped.sections[i]);
                if (count > limits.max_points-active_points || active.size()>=limits.max_pieces)
                    throw ResourceLimit("Exact reference activation output budget exceeded.");
                active_points += count;
                // This field selects the topology contact-search band; it is
                // not a cumulative physical geometry-error certificate.
                const auto band = std::max(uncertainty, clipped.representation_allowance_nm);
                if (erased_island(clipped.sections[i], band))
                    throw std::invalid_argument("Exact activation clipping retains an uncertainty-erased positive component; material=" +
                        std::to_string(materials[i]) + " z_nm=" + std::to_string(ordered[z-1]) + "," +
                        std::to_string(ordered[z]) + " uncertainty_nm=" + std::to_string(band) + ".");
                active.push_back({double(ordered[z-1])*1e-9,double(ordered[z])*1e-9,materials[i],
                    std::move(clipped.sections[i]),band});
            }
        }
        for (auto& slab : active) {
            if (slab.footprint.native().empty()) continue;
            const auto n = points(slab.footprint.native());
            if (n > limits.max_points - output_points || output_pieces >= limits.max_pieces)
                throw ResourceLimit("Reference activation output budget exceeded.");
            output_points += n;
            ++output_pieces;
            emit(std::move(slab));
        }
    }
    return excluded;
}

ReferenceSectionActivation activate_reference_sections(const SectionTopology& reference,
    const std::vector<SectionTopology>& deposited, const SectionGrowthLimits& limits,
    const SectionTopologyLimits& topology_limits, const std::function<void()>& checkpoint)
{
    if (deposited.empty() || deposited.size() > 256 || !limits.max_work || limits.max_work > 200000000 ||
        !limits.max_points || limits.max_points > 2000000 || !limits.max_pieces || limits.max_pieces > 10000)
        throw std::invalid_argument("Invalid fixed-interface activation input or limits.");
    std::size_t work = 0;
    const auto charge = [&](std::size_t n) {
        if (checkpoint) checkpoint();
        if (n > limits.max_work - work) throw ResourceLimit("Reference activation work budget exceeded.");
        work += n;
    };
    std::vector<const DepositedSlab*> ref, source;
    for (const auto& p : reference.pieces()) { charge(1); ref.push_back(&p.slab); }
    for (const auto& topology : deposited) {
        if (topology.bed_z_m() != reference.bed_z_m())
            throw std::invalid_argument("Reference activation bed planes differ.");
        for (const auto& p : topology.pieces()) { charge(1); source.push_back(&p.slab); }
    }
    std::vector<DepositedSlab> active;
    const double excluded = clip_reference_sections(ref,source,limits,charge,[&] { return limits.max_work-work; },
        [&](DepositedSlab slab) { active.push_back(std::move(slab)); });
    return {qualify_section_topology(active,reference.bed_z_m(),topology_limits,checkpoint),excluded};
}

struct ReferenceActivationCursor::Impl {
    SectionGrowthStreamSettings settings;
    std::function<void()> checkpoint;
    LocalSectionPacketVisitor active_packets;
    unsigned hooks = 0;
    std::size_t work = 0, processed_points = 0, output_points = 0, output_pieces = 0;
    std::function<void(std::size_t)> charge;
    std::map<std::size_t,std::vector<DepositedSegment>> materials;
    std::vector<DepositedSegment> complete_reference;
    bool reuse_complete = false;
    std::shared_ptr<LocalSectionSamplingSchedule> sampling_schedule;
    std::unique_ptr<LocalSectionCursor> fixed;
    std::vector<std::unique_ptr<LocalSectionCursor>> sources;
    std::unique_ptr<LocalTopologyBuilder> fixed_topology, active_topology;
    double excluded = 0.;
    bool failed = false, busy = false, finished = false;
    Impl(const std::vector<DepositedSegment>& reference,const std::vector<DepositedSegment>& deposited,
        const SectionGrowthStreamSettings& input,std::function<void(const SectionPiece&)> visitor,
        std::function<std::optional<SectionOwnership>(double)> reference_ownership,
        std::function<void()> check,LocalSectionPacketVisitor packets)
        : settings(input),checkpoint(std::move(check)),active_packets(std::move(packets))
    {
    hooks = bool(active_packets.begin_tile)+bool(active_packets.band)+bool(active_packets.end_tile);
    if (hooks != 0 && hooks != 3) throw std::invalid_argument("Incomplete reference activation packet hooks.");
    if (!visitor || reference.empty() || deposited.empty() || !settings.max_work || settings.max_work > 400000000 ||
        !settings.max_tile_points || settings.max_tile_points > 2000000 || !settings.max_output_points ||
        settings.max_output_points > 128000000 || !settings.max_output_pieces || settings.max_output_pieces > 1000000 ||
        deposited.size() > settings.geometry.geometry.max_segments)
        throw std::invalid_argument("Invalid streamed reference activation inputs or limits.");
    charge = [this](std::size_t n) {
        if (checkpoint) checkpoint();
        if (n > settings.max_work-work) throw ResourceLimit("Streamed reference activation work budget exceeded.");
        work += n;
    };
    std::size_t compared=0;
    if (settings.reuse_complete_deposited_reference) {
        if (!settings.exact_reference_sources || !settings.geometry.geometry.retain_exact_interface_sections)
            throw std::invalid_argument("Complete deposited reference reuse requires exact source composition.");
        if (!settings.geometry.geometry.max_segments || settings.geometry.geometry.max_segments>2000000)
            throw std::invalid_argument("Invalid complete reference snapshot segment limit.");
        // Freeze the full source BEFORE the first caller checkpoint. A later
        // mutation of either input must not turn an equality test into a stale
        // authorization to omit source sampling. Existing segment/work caps
        // bound this synchronous comparison and copy.
        if (reference.size()==deposited.size()) {
            reuse_complete=true;
            for (std::size_t i=0;i<reference.size();++i) {
                if (++compared>settings.max_work)
                    throw ResourceLimit("Complete reference comparison work budget exceeded.");
                const auto& a=reference[i]; const auto& b=deposited[i];
                if ((a.start_m.array()!=b.start_m.array()).any() || (a.end_m.array()!=b.end_m.array()).any() ||
                    a.width_m!=b.width_m || a.height_m!=b.height_m || a.material!=b.material ||
                    a.section!=b.section || a.end_height_m!=b.end_height_m || a.section_source_index!=b.section_source_index) {
                    reuse_complete=false;
                    break;
                }
            }
            if (reuse_complete) {
                if (reference.size()>settings.max_work-compared)
                    throw ResourceLimit("Complete reference snapshot work budget exceeded.");
                complete_reference=reference;
                compared+=reference.size();
            }
        }
    }
    if (settings.reuse_complete_deposited_reference) charge(compared);
    const auto& fixed_input=reuse_complete ? complete_reference : reference;
    if (!reuse_complete) for (const auto& segment : deposited) {
        charge(1);
        materials[segment.material].push_back(segment);
        if (materials.size() > 256) throw ResourceLimit("Reference activation material budget exceeded.");
    }
    auto geometry = settings.geometry;
    if (settings.exact_reference_sources && !geometry.geometry.retain_exact_interface_sections)
        throw std::invalid_argument("Exact reference activation sources require retained exact sections.");
    const bool shared_sampling=geometry.geometry.retain_exact_interface_sections;
    const auto cursor_count = materials.size()+1+std::size_t(shared_sampling);
    std::vector<std::size_t> index_quotas;
    std::size_t index_quota=0;
    if (shared_sampling) {
        std::vector<std::size_t> input_counts{reuse_complete ? 2*fixed_input.size() : reference.size()+deposited.size(),
            fixed_input.size()};
        for (const auto& entry : materials) input_counts.push_back(entry.second.size());
        index_quotas=GrowthDetail::index_work_quotas(geometry.max_index_work,input_counts);
    } else geometry.max_index_work /= cursor_count;
    geometry.max_sections /= cursor_count;
    geometry.max_segment_visits /= cursor_count;
    // Pool actual index occupancy, not equal per-material shares. A large
    // reference and small secondary material can otherwise reject despite
    // fitting the same total memory budget; smaller spatial tiles then become
    // unusable even when they would reduce the live polygon working set.
    const auto reserve_index = [&](const auto& cursor) {
        const auto used=cursor.index_usage();
        geometry.max_tiles-=used.tiles;
        geometry.max_tile_references-=used.references;
    };
    if (shared_sampling) {
        geometry.max_index_work=index_quotas.at(index_quota++);
        std::vector<const std::vector<DepositedSegment>*> schedule_sources{&fixed_input};
        // The old schedule unions reference with all material subsets of an
        // identical deposited vector. A second identical vector preserves that
        // same multiset of endpoints and partition planes, without sources.
        if (reuse_complete) schedule_sources.push_back(&fixed_input);
        for (const auto& entry : materials) schedule_sources.push_back(&entry.second);
        sampling_schedule=std::make_shared<LocalSectionSamplingSchedule>(schedule_sources,settings.origin_m,settings.tile_size_m,
            settings.bed_z_m,settings.maximum_height_m,geometry,checkpoint,settings.interval_subdivisions);
        reserve_index(*sampling_schedule);
    }
    if (!geometry.max_tiles || !geometry.max_tile_references)
        throw ResourceLimit("Reference activation shared index budget exhausted.");
    if (shared_sampling) geometry.max_index_work=index_quotas.at(index_quota++);
    fixed = std::make_unique<LocalSectionCursor>(fixed_input,settings.origin_m,settings.tile_size_m,settings.bed_z_m,
        settings.maximum_height_m,geometry,std::move(reference_ownership),checkpoint,settings.interval_subdivisions,
        nullptr,nullptr,sampling_schedule);
    reserve_index(*fixed);
    geometry.geometry.exact_per_command_composition=settings.exact_reference_sources;
    for (const auto& [material,segments] : materials) {
        if (!geometry.max_tiles || !geometry.max_tile_references)
            throw ResourceLimit("Reference activation shared index budget exhausted.");
        if (shared_sampling) geometry.max_index_work=index_quotas.at(index_quota++);
        sources.push_back(std::make_unique<LocalSectionCursor>(segments,settings.origin_m,settings.tile_size_m,
            settings.bed_z_m,settings.maximum_height_m,geometry,
            std::function<std::optional<SectionOwnership>(double)>{},checkpoint,settings.interval_subdivisions,
            nullptr,nullptr,sampling_schedule));
        reserve_index(*sources.back());
    }
    if (!reuse_complete)
        fixed_topology = std::make_unique<LocalTopologyBuilder>(settings.bed_z_m,settings.topology,[](const SectionPiece&) {},checkpoint,
            settings.connectivity,nullptr,LocalTileOrder::Lexicographic,LocalAttachmentPolicy::AllowFreeComponents);
    active_topology = std::make_unique<LocalTopologyBuilder>(settings.bed_z_m,settings.topology,std::move(visitor),checkpoint,
        settings.connectivity,nullptr,LocalTileOrder::Lexicographic,LocalAttachmentPolicy::AllowFreeComponents);
    }
    std::optional<std::array<int,2>> next() const {
        auto result = fixed->next_tile_index();
        for (const auto& source : sources) {
            const auto candidate = source->next_tile_index();
            if (candidate && (!result || *candidate < *result)) result = candidate;
        }
        return result;
    }
    bool step();
    bool step_complete();
};

bool ReferenceActivationCursor::Impl::step_complete()
{
    const auto next=fixed->next_tile_index();
    if (!next) return false;
    while (fixed->next_tile_index()==next) {
        charge(1);
        fixed->advance_packet({[&](const LocalSectionTile& tile) {
            active_topology->begin_tile(tile);
            if (hooks) active_packets.begin_tile(tile);
        },[&](const LocalSectionTile& tile) {
            if (tile.slabs.empty()) return;
            if (tile.slabs.size()>256)
                throw std::invalid_argument("Complete reference activation has too many partition members.");
            std::vector<const DepositedSlab*> members;
            std::size_t points=0;
            for (const auto& slab:tile.slabs) {
                charge(1);
                if (slab.boundary_uncertainty_nm)
                    throw std::invalid_argument("Complete deposited reference requires zero boundary uncertainty.");
                const auto n=slab.footprint.vertex_count();
                if (n>settings.max_tile_points-points || n>settings.max_output_points-processed_points ||
                    n>settings.max_output_points-output_points || output_pieces>=settings.max_output_pieces)
                    throw ResourceLimit("Complete reference activation geometry budget exceeded.");
                points+=n; processed_points+=n; output_points+=n; ++output_pieces;
                members.push_back(&slab);
            }
            const auto& proof=tile.slabs.front().occupied_partition;
            if (!proof || proof->deposited_overlap_policy()!=DepositedOverlapPolicy::FirstDepositionWins)
                throw std::invalid_argument("Complete reference activation requires deposited partition evidence.");
            if (work==settings.max_work) throw ResourceLimit("Complete reference activation proof budget exhausted.");
            charge(qualify_complete_occupied_partition(members,std::min<std::size_t>(4000000,settings.max_work-work),checkpoint));
            const double area=proof->deposited_excluded_nominal_area_m2();
            if (!std::isfinite(area) || area<0.) throw std::invalid_argument("Invalid complete reference excluded area.");
            const double volume=area*double(height(tile.slabs.front().top_m)-height(tile.slabs.front().bottom_m))*1e-9;
            if (!std::isfinite(volume) || (area>0. && !(volume>0.)))
                throw std::invalid_argument("Unresolved exact excluded reference volume.");
            excluded+=volume;
            active_topology->add_band(tile);
            if (hooks) active_packets.band(tile);
        },[&] {
            active_topology->end_tile();
            if (hooks) active_packets.end_tile();
        }});
    }
    return true;
}

bool ReferenceActivationCursor::Impl::step()
{
    if (reuse_complete) return step_complete();
    charge(sources.size() + 1);
    const auto next = this->next();
    if (!next)
        return false;
    LocalSectionTile header;
    std::vector<LocalSectionCursor*> cursors{fixed.get()};
    for (auto& source : sources)
        cursors.push_back(source.get());
    std::vector<LocalSectionTile> bands(cursors.size());
    std::size_t live_points = 0;
    bool have_header        = false;
    const auto read         = [&](std::size_t i) {
        auto& cursor            = *cursors[i];
        auto& band              = bands[i];
        const bool is_reference = i == 0;
        for (const auto& slab : band.slabs)
            live_points -= slab.footprint.vertex_count();
        band.slabs.clear();
        try {
            while (cursor.next_tile_index() == next && band.slabs.empty())
                cursor.advance_packet(
                    {[&](const LocalSectionTile& tile) {
                         if (have_header && (tile.minimum != header.minimum || tile.maximum != header.maximum))
                             throw std::invalid_argument("Reference activation tile bounds differ.");
                         header      = tile;
                         have_header = true;
                         if (is_reference)
                             fixed_topology->begin_tile(tile);
                     },
                     [&](const LocalSectionTile& tile) {
                         if (is_reference) {
                             try {
                                 fixed_topology->add_band(tile);
                             } catch (const LayerMaterialOverlap& error) {
                                 throw LayerMaterialOverlap(error.bottom_m, error.top_m, error.area_m2, error.material,
                                                            " Fixed reference input: " + std::string(error.what()));
                             } catch (const std::invalid_argument& error) {
                                 throw std::invalid_argument("Fixed reference sampler topology: " + std::string(error.what()));
                             }
                         }
                         for (const auto& slab : tile.slabs) {
                             const auto n = slab.footprint.vertex_count();
                             charge(n + 1);
                             if (n > settings.max_tile_points - live_points)
                                 throw ResourceLimit(
                                     "Reference activation live band point budget exceeded (tile=" + std::to_string(tile.index[0]) + "," +
                                     std::to_string(tile.index[1]) + ", retained=" + std::to_string(live_points) +
                                     ", incoming=" + std::to_string(n) + ", limit=" + std::to_string(settings.max_tile_points) + ").");
                             if (n > settings.max_output_points - processed_points)
                                 throw ResourceLimit("Reference activation cumulative sampled point budget exceeded (processed=" +
                                                     std::to_string(processed_points) + ", incoming=" + std::to_string(n) +
                                                     ", limit=" + std::to_string(settings.max_output_points) + ").");
                             live_points += n;
                             processed_points += n;
                             band.slabs.push_back(slab);
                         }
                     },
                     [&] {
                         if (is_reference)
                             fixed_topology->end_tile();
                     }});
        } catch (const ResourceLimit& error) {
            throw ResourceLimit(std::string(error.what()) + " Reference activation sampler=" +
                                (is_reference ? std::string("fixed") : "source-" + std::to_string(i - 1)) + ".");
        }
    };
    for (std::size_t i = 0; i < cursors.size(); ++i)
        read(i);
    if (!have_header)
        return true;
    active_topology->begin_tile(header);
    if (hooks)
        active_packets.begin_tile(header);
    SectionGrowthLimits clipping;
    clipping.max_points = settings.max_tile_points;
    clipping.max_pieces = std::min<std::size_t>(10000, settings.max_output_pieces);
    auto z              = std::numeric_limits<coord_t>::max();
    for (const auto& band : bands)
        if (!band.slabs.empty())
            z = std::min(z, height(band.slabs.front().bottom_m));
    while (z != std::numeric_limits<coord_t>::max()) {
        std::vector<const DepositedSlab*> ref, source;
        auto top = std::numeric_limits<coord_t>::max();
        for (std::size_t i = 0; i < bands.size(); ++i) {
            const auto& band = bands[i];
            if (band.slabs.empty())
                continue;
            const auto bottom = height(band.slabs.front().bottom_m);
            if (bottom > z) {
                top = std::min(top, bottom);
                continue;
            }
            top = std::min(top, height(band.slabs.front().top_m));
            for (const auto& slab : band.slabs)
                (i == 0 ? ref : source).push_back(&slab);
        }
        if (top == std::numeric_limits<coord_t>::max())
            break;
        LocalSectionTile packet = header;
        const auto contextual_limit = [&](const ResourceLimit& error,const char* stage) {
            return ResourceLimit(std::string(error.what())+" Reference activation stage="+stage+
                ", tile=("+std::to_string(header.index[0])+","+std::to_string(header.index[1])+
                "), band_nm=["+std::to_string(z)+","+std::to_string(top)+"], activation_work="+
                std::to_string(work)+", activation_limit="+std::to_string(settings.max_work)+".");
        };
        try {
            excluded += clip_reference_sections(
                ref, source, clipping, charge, [this] { return settings.max_work-work; },
                [&](DepositedSlab slab) {
                    const auto n = slab.footprint.vertex_count();
                    if (n > settings.max_output_points - output_points || output_pieces >= settings.max_output_pieces)
                        throw ResourceLimit("Reference activation cumulative output budget exceeded.");
                    output_points += n;
                    ++output_pieces;
                    packet.slabs.push_back(std::move(slab));
                },
                std::make_pair(z, top));
        } catch (const ResourceLimit& error) {
            throw contextual_limit(error,"clipping");
        }
        if (!packet.slabs.empty()) {
            try { active_topology->add_band(packet); }
            catch (const ResourceLimit& error) {
                throw contextual_limit(error,"active-topology");
            }
            catch (const LayerMaterialOverlap& error) {
                throw LayerMaterialOverlap(error.bottom_m,error.top_m,error.area_m2,error.material,
                    " Activated reference clipping: " + std::string(error.what()));
            }
            catch (const std::invalid_argument& error) {
                std::string materials;
                for (std::size_t i=0;i<std::min<std::size_t>(packet.slabs.size(),8);++i)
                    materials+=(i ? "," : "")+std::to_string(packet.slabs[i].material);
                throw std::invalid_argument("Activated reference clipping topology: " + std::string(error.what())+
                    " tile=("+std::to_string(header.index[0])+","+std::to_string(header.index[1])+
                    "), band_nm=["+std::to_string(z)+","+std::to_string(top)+"], materials=["+materials+"].");
            }
            if (hooks) {
                try { active_packets.band(packet); }
                catch (const ResourceLimit& error) {
                    throw contextual_limit(error,"active-packet-consumer");
                }
            }
        }
        z = top;
        for (std::size_t i = 0; i < bands.size(); ++i)
            if (!bands[i].slabs.empty() && height(bands[i].slabs.front().top_m) <= z)
                read(i);
    }
    active_topology->end_tile();
    if (hooks)
        active_packets.end_tile();
    return true;
}

ReferenceActivationCursor::ReferenceActivationCursor(const std::vector<DepositedSegment>& reference,
    const std::vector<DepositedSegment>& deposited,const SectionGrowthStreamSettings& settings,
    std::function<void(const SectionPiece&)> visitor,
    std::function<std::optional<SectionOwnership>(double)> ownership,std::function<void()> checkpoint,
    LocalSectionPacketVisitor packets)
    : m_impl(std::make_unique<Impl>(reference,deposited,settings,std::move(visitor),std::move(ownership),
        std::move(checkpoint),std::move(packets))) {}
ReferenceActivationCursor::~ReferenceActivationCursor() = default;
std::optional<std::array<int,2>> ReferenceActivationCursor::next_tile_index() const
{
    if (m_impl->failed || m_impl->busy || m_impl->finished) {
        m_impl->failed = true;
        throw std::logic_error("Reference activation cursor is not available.");
    }
    return m_impl->next();
}
bool ReferenceActivationCursor::advance()
{
    next_tile_index();
    m_impl->busy = true;
    try {
        const bool result = m_impl->step();
        if (m_impl->failed) throw std::logic_error("Reference activation cursor was reentered.");
        m_impl->busy = false;
        return result;
    } catch (...) { m_impl->failed = true; m_impl->busy = false; throw; }
}
std::vector<std::size_t> ReferenceActivationCursor::active_component_mapping()
{
    if (m_impl->failed || m_impl->finished)
        throw std::logic_error("Reference activation cursor is not available.");
    try {
        auto mapping = m_impl->active_topology->provisional_component_mapping();
        if (m_impl->failed) throw std::logic_error("Reference activation cursor was reentered.");
        return mapping;
    } catch (...) { m_impl->failed = true; throw; }
}

ReferenceActivationStreamResult ReferenceActivationCursor::finish()
{
    if (next_tile_index()) {
        m_impl->failed = true;
        throw std::logic_error("Reference activation cursor is incomplete.");
    }
    m_impl->busy = true;
    try {
        auto result = [&]() -> ReferenceActivationStreamResult {
            if (m_impl->reuse_complete) {
                // The authenticated complete path forwards identical geometry:
                // one finished stream proves both roles, including its labels.
                auto active = m_impl->active_topology->finish();
                m_impl->charge(active.component_mapping().size() + 1);
                auto reference = active;
                return {std::move(reference),std::move(active),m_impl->excluded,m_impl->processed_points};
            }
            return {m_impl->fixed_topology->finish(),m_impl->active_topology->finish(),
                m_impl->excluded,m_impl->processed_points};
        }();
        if (m_impl->failed) throw std::logic_error("Reference activation cursor was reentered.");
        m_impl->busy = false;
        m_impl->finished = true;
        return result;
    } catch (...) { m_impl->failed = true; m_impl->busy = false; throw; }
}

ReferenceActivationStreamResult visit_reference_activation(
    const std::vector<DepositedSegment>& reference, const std::vector<DepositedSegment>& deposited,
    const SectionGrowthStreamSettings& settings, const std::function<void(const SectionPiece&)>& visitor,
    const std::function<std::optional<SectionOwnership>(double)>& reference_ownership,
    const std::function<void()>& checkpoint, const LocalSectionPacketVisitor& active_packets)
{
    ReferenceActivationCursor cursor(reference,deposited,settings,visitor,reference_ownership,checkpoint,active_packets);
    while (cursor.advance()) {}
    return cursor.finish();
}

static SectionGrowthStreamResult visit_section_growth_impl(const std::vector<DepositedSegment>& previous,
                                                     const std::vector<DepositedSegment>& next,
                                                     const SectionGrowthStreamSettings& settings,
                                                     const SectionGrowthStreamVisitor& visitor,
                                                     const std::function<std::optional<SectionOwnership>(double)>& previous_ownership,
                                                     const std::function<std::optional<SectionOwnership>(double)>& next_ownership,
                                                     const std::function<void()>& checkpoint,
                                                     const std::vector<DepositedSegment>* previous_partitions,
                                                     const std::vector<DepositedSegment>* next_partitions,
                                                     const std::vector<DepositedSegment>* reference = nullptr,
                                                     std::array<double,2>* excluded = nullptr)
{
    const unsigned hooks = bool(visitor.next_packets.begin_tile) + bool(visitor.next_packets.band) + bool(visitor.next_packets.end_tile);
    if (hooks != 0 && hooks != 3)
        throw std::invalid_argument("Incomplete local growth packet hooks.");
    if (!visitor.previous || !visitor.next || !visitor.added || !settings.max_work || settings.max_work > 400000000 ||
        !settings.max_tile_points || settings.max_tile_points > 2000000 || !settings.max_output_points ||
        settings.max_output_points > 128000000 || !settings.max_output_pieces || settings.max_output_pieces > 1000000 ||
        !settings.max_component_links || settings.max_component_links > 1000000)
        throw std::invalid_argument("Invalid local growth visitor or limits.");
    std::unique_ptr<LocalSectionCursor> old_cursor,new_cursor;
    if (!reference) {
        old_cursor = std::make_unique<LocalSectionCursor>(previous, settings.origin_m, settings.tile_size_m, settings.bed_z_m,
            settings.maximum_height_m,settings.geometry, previous_ownership, checkpoint, settings.interval_subdivisions, nullptr,
            previous_partitions);
        new_cursor = std::make_unique<LocalSectionCursor>(next, settings.origin_m, settings.tile_size_m, settings.bed_z_m,
            settings.maximum_height_m,settings.geometry, next_ownership, checkpoint, settings.interval_subdivisions, nullptr,next_partitions);
    }
    std::vector<SectionPiece> old_tile, new_tile;
    std::size_t live_points = 0, output_points = 0, output_pieces = 0, work = 0;
    const auto charge = [&](std::size_t n) {
        if (checkpoint)
            checkpoint();
        if (n > settings.max_work - work)
            throw ResourceLimit("Local growth work budget exceeded.");
        work += n;
    };
    const auto capture = [&](std::vector<SectionPiece>& tile, const SectionPiece& p, const auto& callback) {
        const auto n = p.slab.footprint.vertex_count();
        charge(n + 1);
        if (n > settings.max_tile_points - live_points)
            throw ResourceLimit("Local growth tile geometry budget exceeded.");
        live_points += n;
        tile.push_back(p);
        callback(p);
    };
    LocalTopologyBuilder old_builder(
        settings.bed_z_m, settings.topology, [&](const SectionPiece& p) { capture(old_tile, p, visitor.previous); }, checkpoint,
        settings.connectivity, nullptr, LocalTileOrder::Lexicographic, LocalAttachmentPolicy::AllowFreeComponents);
    LocalTopologyBuilder new_builder(
        settings.bed_z_m, settings.topology, [&](const SectionPiece& p) { capture(new_tile, p, visitor.next); }, checkpoint,
        settings.connectivity, nullptr, LocalTileOrder::Lexicographic, LocalAttachmentPolicy::AllowFreeComponents);
    const auto packets = [](LocalTopologyBuilder& builder) {
        return LocalSectionPacketVisitor{[&builder](const LocalSectionTile& h) { builder.begin_tile(h); },
                                         [&builder](const LocalSectionTile& p) { builder.add_band(p); },
                                         [&builder] { builder.end_tile(); }};
    };
    const auto old_packets = packets(old_builder);
    const LocalSectionPacketVisitor new_packets{[&](const LocalSectionTile& h) {
                                                    new_builder.begin_tile(h);
                                                    if (hooks)
                                                        visitor.next_packets.begin_tile(h);
                                                },
                                                [&](const LocalSectionTile& p) {
                                                    new_builder.add_band(p);
                                                    if (hooks)
                                                        visitor.next_packets.band(p);
                                                },
                                                [&] {
                                                    new_builder.end_tile();
                                                    if (hooks)
                                                        visitor.next_packets.end_tile();
                                                }};
    std::unique_ptr<ReferenceActivationCursor> old_reference,new_reference;
    if (reference) {
        // Packet consumers below perform the shared chronology qualification.
        // The activation cursors separately certify each immutable reference.
        old_reference = std::make_unique<ReferenceActivationCursor>(*reference,previous,settings,
            [](const SectionPiece&) {},previous_ownership,checkpoint,old_packets);
        new_reference = std::make_unique<ReferenceActivationCursor>(*reference,next,settings,
            [](const SectionPiece&) {},next_ownership,checkpoint,new_packets);
    }
    std::set<std::pair<std::size_t, std::size_t>> links;
    const auto contact = [&](std::size_t a, std::size_t b) {
        charge(1);
        const auto link = std::make_pair(a, b);
        if (!links.count(link) && links.size() >= settings.max_component_links)
            throw ResourceLimit("Local growth component-link budget exceeded.");
        links.insert(link);
    };
    const auto added = [&](const SectionPiece& p) {
        const auto n = p.slab.footprint.vertex_count();
        if (n > settings.max_output_points - output_points || output_pieces >= settings.max_output_pieces)
            throw ResourceLimit("Local growth output budget exceeded.");
        charge(n + 1);
        output_points += n;
        ++output_pieces;
        visitor.added(p);
    };
    for (;;) {
        const auto a = reference ? old_reference->next_tile_index() : old_cursor->next_tile_index();
        const auto b = reference ? new_reference->next_tile_index() : new_cursor->next_tile_index();
        if (!a && !b)
            break;
        const auto index = (!a ? *b : (!b ? *a : std::min(*a, *b)));
        if (a && *a == index) {
            if (reference) old_reference->advance(); else old_cursor->advance(old_packets);
        }
        if (b && *b == index) {
            if (reference) new_reference->advance(); else new_cursor->advance(new_packets);
        }
        const auto remaining = [&] { return settings.max_work - work; };
        GrowthDetail::subtract_partition(old_tile, new_tile, true, settings.max_tile_points, charge, remaining, contact, added);
        GrowthDetail::subtract_partition(new_tile, old_tile, false, settings.max_tile_points, charge, remaining, contact, added);
        old_tile.clear();
        new_tile.clear();
        live_points = 0;
    }
    if (reference) {
        const auto old_proof = old_reference->finish(), new_proof = new_reference->finish();
        *excluded = {old_proof.excluded_nominal_volume_m3,new_proof.excluded_nominal_volume_m3};
    }
    auto old_certificate = old_builder.finish(), new_certificate = new_builder.finish();
    const auto missing = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> mapping(old_certificate.components(), missing);
    for (const auto& [old_label, new_label] : links) {
        charge(1);
        auto& value       = mapping.at(old_certificate.component_mapping().at(old_label));
        const auto target = new_certificate.component_mapping().at(new_label);
        if (value != missing && value != target)
            throw std::invalid_argument("Existing local component split during deposition.");
        value = target;
    }
    if (std::find(mapping.begin(), mapping.end(), missing) != mapping.end())
        throw std::invalid_argument("Missing local growth component coverage.");
    return {std::move(old_certificate), std::move(new_certificate), std::move(mapping),
        reference ? LocalSectionMetrics{} : old_cursor->metrics(), reference ? LocalSectionMetrics{} : new_cursor->metrics()};
}

SectionGrowthStreamResult visit_local_section_growth(const std::vector<DepositedSegment>& previous,
    const std::vector<DepositedSegment>& next,const SectionGrowthStreamSettings& settings,const SectionGrowthStreamVisitor& visitor,
    const std::function<std::optional<SectionOwnership>(double)>& previous_ownership,
    const std::function<std::optional<SectionOwnership>(double)>& next_ownership,const std::function<void()>& checkpoint,
    const std::vector<DepositedSegment>* previous_partitions,const std::vector<DepositedSegment>* next_partitions)
{
    return visit_section_growth_impl(previous,next,settings,visitor,previous_ownership,next_ownership,checkpoint,
        previous_partitions,next_partitions);
}
ReferenceGrowthStreamResult visit_reference_section_growth(const std::vector<DepositedSegment>& reference,
    const std::vector<DepositedSegment>& previous,const std::vector<DepositedSegment>& next,
    const SectionGrowthStreamSettings& settings,const SectionGrowthStreamVisitor& visitor,
    const std::function<std::optional<SectionOwnership>(double)>& ownership,const std::function<void()>& checkpoint)
{
    std::array<double,2> excluded{};
    auto result = visit_section_growth_impl(previous,next,settings,visitor,ownership,ownership,checkpoint,nullptr,nullptr,
        &reference,&excluded);
    return {std::move(result.previous),std::move(result.next),std::move(result.previous_components),excluded[0],excluded[1]};
}
} // namespace Slic3r::Thermoelastic
