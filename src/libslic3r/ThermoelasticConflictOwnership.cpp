#include "ThermoelasticConflictOwnership.hpp"
#include "ThermoelasticSourceOwnership.hpp"
#include "ThermoelasticToolpath.hpp"
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticSectionFootprint.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticExactSections.hpp"
#include "GCode/GCodeProcessor.hpp"
#include "GCode/DepositionProvenance.hpp"
#include "AABBTreeIndirect.hpp"
#include "BoundingBox.hpp"
#include "ClipperUtils.hpp"
#include "ExPolygon.hpp"
#include "Polygon.hpp"
#include "ExtrusionEntity.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace Slic3r::Thermoelastic {
class DepositedOwnershipPartitionFactory {
public:
    static std::shared_ptr<const AuthenticatedDepositedPartition> make(const SampledSectionConflict& sampled,
        const std::vector<DepositedSegment>& ordered, std::size_t max_work,
        const std::function<void()>& checkpoint, std::size_t& consumed)
    {
        if (!sampled.geometry_limits().retain_exact_interface_sections || !max_work || max_work>4000000)
            throw std::invalid_argument("Deposited partition requires bounded exact sampler geometry.");
        std::size_t work=0;
        const auto charge=[&](std::size_t amount) {
            if (checkpoint) checkpoint();
            if (amount>max_work-work) throw ResourceLimit("Deposited partition work budget exceeded.");
            work+=amount;
        };
        const auto available=[&] {
            if (work==max_work) throw ResourceLimit("Deposited partition work budget exhausted.");
            return max_work-work;
        };
        // Only evidence composed from entry-time frozen inputs may supply
        // geometry directly. Borrowed inputs could change in a callback after
        // sampling but before source authentication; those require replay.
        if (sampled.m_frozen_inputs && !ordered.empty() && sampled.original_regions().size()==1) {
            const auto& original=sampled.original_regions().front();
            bool single=original.footprint.is_exact() && !original.footprint.empty();
            if (ordered.size()>4096) throw ResourceLimit("Deposited command count budget exceeded.");
            for (const auto& bead:ordered) {
                charge(1);
                if (bead.material>=256) throw std::invalid_argument("Invalid first-deposition material ID.");
                single=single && bead.material==original.material;
            }
            if (single) {
                const auto vertices=original.footprint.vertex_count();
                // Preserve the primitive's conservative occupied+owner alias
                // bound even though both share the immutable sampled handle.
                if (vertices>2048) throw ResourceLimit("First-deposition retained vertex budget exceeded.");
                charge(vertices+1);
                auto token=std::shared_ptr<AuthenticatedDepositedPartition>(new AuthenticatedDepositedPartition);
                token->regions.push_back(original);
                token->identity=sampled.m_identity;
                token->work=work;
                consumed=work;
                return token;
            }
        }
        // A material which never returns after its chronological run owns the
        // same run union already composed by the frozen sampler. Reuse only
        // those complete runs; A/B/A must still replay individual commands.
        // Keep the existing streamed incoming-cut and live-state ceilings,
        // even though a larger precomposed run could be geometrically valid.
        std::vector<ExactDepositedCommandSection> runs;
        bool reuse_runs=sampled.m_frozen_inputs && !ordered.empty() &&
            sampled.original_regions().size()>1 && sampled.original_regions().size()<=256;
        if (reuse_runs) {
            std::array<const DepositedSection::Region*,256> originals{};
            std::array<bool,256> seen{};
            std::size_t vertices=0;
            for (const auto& region:sampled.original_regions()) {
                charge(1);
                if (region.material>=256 || originals[region.material] || !region.footprint.is_exact() ||
                    region.footprint.empty() || region.footprint.vertex_count()>1024 ||
                    region.footprint.vertex_count()>4096-vertices) {
                    reuse_runs=false;
                    break;
                }
                vertices+=region.footprint.vertex_count();
                originals[region.material]=&region;
            }
            if (ordered.size()>4096) throw ResourceLimit("Deposited command count budget exceeded.");
            if (reuse_runs) for (const auto& bead:ordered) {
                charge(1);
                if (bead.material>=256) throw std::invalid_argument("Invalid first-deposition material ID.");
                if (!runs.empty() && runs.back().material==bead.material) continue;
                if (seen[bead.material] || !originals[bead.material]) {
                    reuse_runs=false;
                    break;
                }
                seen[bead.material]=true;
                runs.push_back({bead.material,originals[bead.material]->footprint.exact()});
            }
            reuse_runs=reuse_runs && runs.size()==sampled.original_regions().size();
            if (!reuse_runs) runs.clear();
        }
        // Reconstruct one authenticated command at a time. Retaining all cuts
        // makes densely segmented paths consume an input-storage budget even
        // when their occupied union remains small.
        const auto command_at=[&](std::size_t index,std::size_t remaining_work) {
            if (checkpoint) checkpoint();
            if (reuse_runs) return ExactDepositedCommandResult{runs.at(index),1};
            const auto& bead=ordered.at(index);
            auto geometry=sampled.geometry_limits();
            geometry.max_boundary_points=std::min<std::size_t>(geometry.max_boundary_points,1024);
            const auto cut=section_deposition({bead},sampled.conflict().bottom_m,sampled.bed_z_m(),
                geometry,checkpoint,nullptr,sampled.clip());
            const auto cut_work=cut.boundary_points+1;
            if (cut_work>=remaining_work)
                throw ResourceLimit("Deposited command construction budget exhausted.");
            if (cut.regions.empty()) return ExactDepositedCommandResult{{bead.material,{}},cut_work};
            if (cut.regions.size()!=1) throw std::invalid_argument("One deposited command produced multiple materials.");
            auto exact=make_exact_section(cut.regions.front().footprint.native(),remaining_work-cut_work,checkpoint);
            return ExactDepositedCommandResult{{bead.material,std::move(exact.section)},cut_work+exact.work};
        };
        auto partition=exact_first_deposition_partition_streamed(reuse_runs ? runs.size() : ordered.size(),
            command_at,available(),checkpoint);
        charge(partition.work);
        ExactSection expected;
        for (const auto& original:sampled.original_regions()) {
            if (!original.footprint.is_exact())
                throw std::invalid_argument("Deposited partition requires exact sampled material unions.");
            charge(original.footprint.vertex_count()+1);
            auto united=boolean_exact_sections(expected,original.footprint.exact(),ExactSectionOperation::Union,available(),checkpoint);
            charge(united.work);
            expected=std::move(united.section);
        }
        const auto difference=boolean_exact_sections(partition.occupied_union,expected,
            ExactSectionOperation::SymmetricDifference,available(),checkpoint);
        charge(difference.work);
        if (!difference.section.empty())
            throw std::invalid_argument("Deposited-order partition changes the complete sampled union.");
        auto token=std::shared_ptr<AuthenticatedDepositedPartition>(new AuthenticatedDepositedPartition);
        for (auto& region:partition.regions)
            token->regions.push_back({region.material,SectionFootprint(std::move(region.envelope))});
        token->identity=sampled.m_identity;
        token->m_overlap_area_m2=partition.cross_material_overlap_area_m2;
        token->work=work;
        consumed=work;
        return token;
    }
};
class NativeOwnershipPartitionFactory {
public:
    static std::shared_ptr<const AuthenticatedSectionPartition> make(const SampledSectionConflict& sampled,
        const SectionOwnership& native, NativeInterfaceVariant variant, std::size_t max_work,
        const std::function<void()>& checkpoint, std::size_t& consumed, std::size_t& native_winding_consumed,
        bool allow_single_sided_native_boundary,
        bool allow_dual_native_boundary,
        const std::function<void(const NativeOwnershipPartitionFailure&)>& observe_failure)
    {
        const auto& original = sampled.original_regions();
        if (original.size() != 2 || sampled.geometry_limits().ownership_boundary_grid_units != 4 ||
            (variant != NativeInterfaceVariant::FavorFirst && variant != NativeInterfaceVariant::FavorSecond))
            throw std::invalid_argument("Interface sensitivity requires two materials and the explicit 7 nm representation scenario.");
        std::array<ExPolygons, 2> envelopes{original[0].footprint.native(), original[1].footprint.native()}, masks;
        const auto bounds = get_extents(sampled.conflict().section_footprint);
        const auto intersects = [&](const Polygon& ring) {
            const auto box = get_extents(ring.points);
            return box.min.x() <= bounds.max.x() && box.max.x() >= bounds.min.x() &&
                box.min.y() <= bounds.max.y() && box.max.y() >= bounds.min.y();
        };
        for (const auto& region : native.regions) {
            const std::size_t slot = region.material == original[0].material ? 0 :
                region.material == original[1].material ? 1 : 2;
            if (slot == 2) throw std::invalid_argument("Interface source contains an unexpected material.");
            for (const auto& polygon : region.footprint.native()) {
                if (checkpoint) checkpoint();
                const bool outer = intersects(polygon.contour);
                auto retained = polygon;
                retained.holes.erase(std::remove_if(retained.holes.begin(), retained.holes.end(),
                    [&](const Polygon& hole) { return !intersects(hole); }), retained.holes.end());
                if (!outer) {
                    if (retained.holes.empty()) continue;
                    retained.contour = std::move(retained.holes.back());
                    retained.holes.pop_back();
                }
                masks[slot].push_back(std::move(retained));
            }
        }
        auto token = std::shared_ptr<AuthenticatedSectionPartition>(new AuthenticatedSectionPartition);
        bool callback_failed=false;
        const auto checked_checkpoint = [&] {
            if (checkpoint) try { checkpoint(); } catch (...) { callback_failed=true; throw; }
        };
        try {
            if (sampled.geometry_limits().retain_exact_interface_sections) {
                auto selected=exact_partition_uncertain_interface(masks[0],masks[1],envelopes,
                    variant==NativeInterfaceVariant::FavorFirst ? ExactInterfaceChoice::FavorFirst : ExactInterfaceChoice::FavorSecond,
                    7,std::min<std::size_t>(max_work,4000000),checked_checkpoint,NativeWitnessFill::NonzeroWinding,
                    allow_single_sided_native_boundary,allow_dual_native_boundary);
                for (std::size_t i=0;i<2;++i)
                    token->regions.push_back({original[i].material,SectionFootprint(std::move(selected.materials[i]))});
                token->work=selected.work;
                token->interface_allowance_nm=selected.interface_allowance_nm;
                token->representation_allowance_nm=selected.representation_allowance_nm;
                token->dual_native_boundary_area_m2=selected.dual_native_boundary_area_m2;
                token->native_winding_work=selected.native_winding_work;
            } else {
                auto variants = partition_uncertain_interface({}, {}, masks[0], masks[1], 7,
                    std::min<std::size_t>(max_work, 4000000), checked_checkpoint, &envelopes,NativeWitnessFill::NonzeroWinding,
                    allow_single_sided_native_boundary);
                auto& selected = variant == NativeInterfaceVariant::FavorFirst ? variants.favor_first : variants.favor_second;
                for (std::size_t i = 0; i < 2; ++i)
                    token->regions.push_back({original[i].material, std::move(selected[i])});
                token->work = variants.work;
                token->interface_allowance_nm = variants.interface_allowance_nm;
                token->representation_allowance_nm = variants.representation_allowance_nm;
            }
        } catch (const std::invalid_argument& error) {
            if (callback_failed) throw;
            if (observe_failure) {
                const auto bounded = [](const ExPolygons& polygons) {
                    if (polygons.size()>2048) return false;
                    std::size_t points=0;
                    for (const auto& polygon:polygons) {
                        if (polygon.holes.size()>2048 || polygon.contour.points.size()<3 || polygon.contour.points.size()>2048-points) return false;
                        points+=polygon.contour.points.size();
                        for (const auto& hole:polygon.holes) {
                            if (hole.points.size()<3 || hole.points.size()>2048-points) return false;
                            points+=hole.points.size();
                        }
                    }
                    return true;
                };
                if (bounded(masks[0]) && bounded(masks[1]) && bounded(envelopes[0]) && bounded(envelopes[1])) {
                    try {
                        const std::string message=error.what();
                        observe_failure({sampled,masks,envelopes,variant,std::min<std::size_t>(max_work,4000000),
                            allow_single_sided_native_boundary,allow_dual_native_boundary,message});
                    } catch (...) {
                        // Investigation cannot change the original failure.
                    }
                }
            }
            throw std::invalid_argument(std::string(error.what()) + " Materials=" + std::to_string(original[0].material) +
                "," + std::to_string(original[1].material) + "; section_z_nm=" +
                std::to_string(sampled.conflict().bottom_m*1e9) + ".");
        }
        token->identity = sampled.m_identity;
        consumed = token->work;
        native_winding_consumed = token->native_winding_work;
        return token;
    }
};
namespace {
struct PreparedOwnershipWork {
    // Successful debits only; failed operations keep their existing internal
    // diagnostics. Native winding is an exact-factory subset of partition work.
    std::size_t source = 0, partition = 0, native_winding = 0;
    std::string context(std::size_t limit, std::size_t remaining) const
    {
        const auto charged=limit-remaining;
        return ", prepared_charged_work="+std::to_string(charged)+
            ", prepared_source_work="+std::to_string(source)+
            ", prepared_partition_work="+std::to_string(partition)+
            ", prepared_native_winding_work="+std::to_string(native_winding)+
            ", prepared_other_work="+std::to_string(charged-source-partition);
    }
};
using OwnershipTree = AABBTreeIndirect::Tree<3, double>;
OwnershipTree::BoundingBox bead_bounds(const DepositedSegment& bead)
{
    if (!bead.start_m.allFinite() || !bead.end_m.allFinite() || bead.start_m.cwiseAbs().maxCoeff() > 10. ||
        bead.end_m.cwiseAbs().maxCoeff() > 10. || !std::isfinite(bead.width_m) || bead.width_m <= 0. || bead.width_m > 10.)
        throw std::invalid_argument("Invalid native ownership bead bounds.");
    // Reject bounds only, never geometry inflation. Cover native XY rounding
    // and the float radius passed to the existing polygon-offset primitive.
    const double margin = 2e-9 + std::numeric_limits<float>::epsilon() * bead.width_m;
    const double radius = .5 * bead.width_m + margin;
    Eigen::Vector3d lo = bead.start_m.cwiseMin(bead.end_m), hi = bead.start_m.cwiseMax(bead.end_m);
    lo.head<2>().array() -= radius;
    hi.head<2>().array() += radius;
    lo.z() = std::min(bead.start_m.z() - bead.height_m,
        bead.end_m.z() - bead.end_height_m.value_or(bead.height_m)) - margin;
    hi.z() += margin;
    return {lo, hi};
}
OwnershipTree::BoundingBox conflict_bounds(const LayerMaterialOverlap& conflict)
{
    const auto box = get_extents(conflict.section_footprint);
    return {Eigen::Vector3d(double(box.min.x()) * 1e-9, double(box.min.y()) * 1e-9, conflict.bottom_m),
        Eigen::Vector3d(double(box.max.x()) * 1e-9, double(box.max.y()) * 1e-9, conflict.bottom_m)};
}
SectionOwnership resolve_overlap(const GCodeProcessorResult& result,
                                              const LayerMaterialOverlap& conflict,
                                              const std::vector<DepositedSegment>& nominal,
                                              const std::vector<std::size_t>& move_indices,
                                              double bed_z_m,
                                              const std::vector<std::size_t>& filament_to_material,
                                              const ConflictOwnershipLimits& limits,
                                              const std::function<void()>& checkpoint,
                                              std::size_t* remaining = nullptr,
                                              bool retain_unique = false,
                                              const std::vector<std::size_t>* candidates = nullptr,
                                              const SampledSectionConflict* sampled = nullptr,
                                              NativeOwnershipSnapshot::ContainmentCache* containment_cache = nullptr,
                                              PreparedOwnershipWork* prepared_work = nullptr)
{
    if (!result.deposition_provenance || !result.native_ownership || !conflict.section_footprint_complete ||
        conflict.section_footprint.empty() || !std::isfinite(conflict.bottom_m) || conflict.bottom_m != conflict.top_m ||
        !std::isfinite(bed_z_m) || std::abs(bed_z_m) > 10. || std::abs(conflict.bottom_m) > 10. || filament_to_material.empty() ||
        filament_to_material.size() > 256 || limits.max_moves == 0 || limits.max_moves > 2000000 || limits.max_segments == 0 ||
        limits.max_segments > 2000000 || limits.max_polygon_checks == 0 || limits.max_polygon_checks > 4096 ||
        limits.max_polygon_work == 0 || limits.max_polygon_work > 20000000 ||
        limits.max_source_work == 0 || limits.max_source_work > 20000000 ||
        limits.source_centerline_boundary_allowance_nm > 7 ||
        ((limits.source_centerline_boundary_allowance_nm || limits.allow_single_sided_native_boundary || limits.allow_dual_native_boundary) &&
            limits.interface_variant == NativeInterfaceVariant::Strict) ||
        (limits.allow_dual_native_boundary && (!sampled || retain_unique || !limits.geometry.retain_exact_interface_sections)))
        throw std::invalid_argument("Native conflict ownership requires complete bounded live source evidence.");
    std::size_t work = 0, checks = 0;
    const auto charge_total = [&](std::size_t n) {
        if (remaining) {
            if (n > *remaining) throw ResourceLimit("Prepared native ownership total work budget exceeded"+
                (prepared_work ? prepared_work->context(limits.max_total_work,*remaining) : std::string{})+".");
            *remaining -= n;
        }
    };
    const auto charge = [&](const ExPolygons& polygons) {
        if (checkpoint)
            checkpoint();
        for (const auto& polygon : polygons) {
            const auto ring = [&](const Polygon& ring) {
                if (ring.points.size() < 3)
                    throw std::invalid_argument("Invalid native conflict polygon.");
                if (ring.points.size() > limits.max_polygon_work - work)
                    throw ResourceLimit("Native conflict ownership polygon work exceeded.");
                work += ring.points.size();
                charge_total(ring.points.size());
                for (const auto& point : ring.points)
                    if (std::abs(double(point.x())) > 1e10 || std::abs(double(point.y())) > 1e10)
                        throw std::invalid_argument("Native conflict coordinate exceeds bounds.");
            };
            ring(polygon.contour);
            for (const auto& hole : polygon.holes)
                ring(hole);
        }
    };
    charge(conflict.section_footprint);
    const auto bounds = conflict_bounds(conflict);
    if (nominal.size() != move_indices.size() || nominal.size() > limits.max_segments)
        throw std::invalid_argument("Native conflict source geometry is incomplete.");
    std::vector<SourceBead> contributors;
    std::map<std::size_t, ExPolygons> material_sections;
    // Authenticated sampler evidence already establishes the actual conflict.
    // Only unsampled requests and unique-material retention consume a rebuilt
    // material union; the sampled partition uses its original envelopes below.
    const bool reconstruct_material_sections = !sampled || retain_unique;
    const double z = conflict.bottom_m;
    const std::size_t count = candidates ? candidates->size() : nominal.size();
    for (std::size_t candidate = 0; candidate < count; ++candidate) {
        const std::size_t i = candidates ? (*candidates)[candidate] : candidate;
        if (i >= nominal.size()) throw std::invalid_argument("Invalid native ownership spatial index.");
        if (checkpoint)
            checkpoint();
        charge_total(1);
        auto bead           = nominal[i];
        const auto filament = bead.material;
        if (filament >= filament_to_material.size())
            throw std::invalid_argument("Missing native conflict material mapping.");
        bead.material       = filament_to_material[filament];
        if (!bead_bounds(bead).intersects(bounds)) continue;
        if (++checks > limits.max_polygon_checks)
            throw ResourceLimit("Native conflict contributor budget exceeded.");
        const auto sliced = section_deposition({bead}, z, bed_z_m, limits.geometry, checkpoint);
        if (sliced.regions.empty())
            continue;
        const auto& footprint = sliced.regions.front().footprint.native();
        charge(footprint);
        charge(conflict.section_footprint);
        auto contact = intersection_ex(footprint, conflict.section_footprint);
        charge(contact);
        if (contact.empty())
            continue;
        if (move_indices[i] >= result.moves.size()) throw std::invalid_argument("Native conflict source move is unavailable.");
        const auto& move = result.moves[move_indices[i]];
        contributors.push_back({bead, move.deposition_source_id, move.extruder_id, move_indices[i]});
        // Reconstruct original geometry before testing the conflict. Repeatedly
        // clipping each bead to C and then intersecting those clipped unions
        // creates extra integer-rounded vertices on C's own boundary.
        if (reconstruct_material_sections) {
            auto& pieces = material_sections[bead.material];
            pieces.insert(pieces.end(), footprint.begin(), footprint.end());
        }
    }
    if (contributors.empty())
        throw std::invalid_argument("No actual deposition contributes to the conflict.");
    // Verify the requested mask is genuinely occupied by different materials,
    // not just an arbitrary native region supplied by the caller.
    ExPolygons occupied, verified_conflicts;
    for (auto& [material, pieces] : material_sections) {
        charge(pieces);
        pieces = union_ex(pieces);
        charge(pieces);
        charge(occupied);
        const auto overlap = intersection_ex(occupied, pieces);
        charge(overlap);
        verified_conflicts.insert(verified_conflicts.end(), overlap.begin(), overlap.end());
        occupied.insert(occupied.end(), pieces.begin(), pieces.end());
    }
    if (reconstruct_material_sections) {
        charge(verified_conflicts);
        charge(conflict.section_footprint);
    }
    if (!retain_unique && !sampled) {
        const auto missing = diff_ex(conflict.section_footprint, verified_conflicts);
        if (!missing.empty()) {
            std::string details;
            const auto& ring = missing.front().contour.points;
            double twice_area = 0.;
            if (ring.size() >= 3) {
                const auto origin = ring.front();
                for (std::size_t i = 1; i + 1 < ring.size(); ++i) {
                    const auto a = (ring[i] - origin).cast<double>().eval();
                    const auto b = (ring[i + 1] - origin).cast<double>().eval();
                    twice_area += a.x() * b.y() - a.y() * b.x();
                }
            }
            for (std::size_t i = 0; i < std::min<std::size_t>(ring.size(), 12); ++i)
                details += std::to_string(ring[i].x()) + "," + std::to_string(ring[i].y()) + ";";
            throw std::invalid_argument("Conflict mask is not covered by actual material intersections (z_nm=" +
                std::to_string(z * 1e9) + ", first_contour_nm2=" + std::to_string(.5 * twice_area) +
                ", holes=" + std::to_string(missing.front().holes.size()) + ", points=" + details + ").");
        }
    }
    if (retain_unique) {
        charge(occupied);
        if (!diff_ex(conflict.section_footprint, occupied).empty())
            throw std::invalid_argument("Expanded material conflict includes original void; no material owner is established.");
    }
    SourceOwnershipLimits source_limits;
    source_limits.centerline_boundary_allowance_nm = limits.source_centerline_boundary_allowance_nm;
    source_limits.max_work = remaining ? std::min(limits.max_source_work,*remaining) : limits.max_source_work;
    if (!source_limits.max_work) throw ResourceLimit("Prepared native ownership total work budget exceeded"+
        (prepared_work && remaining ? prepared_work->context(limits.max_total_work,*remaining) : std::string{})+".");
    std::size_t source_work = 0;
    SectionOwnership ownership{z,{}};
    bool source_callback_failed=false;
    try {
        ownership=result.native_ownership->section(result.deposition_provenance,contributors,z,filament_to_material,
            source_limits,[&] {
                if (checkpoint) try { checkpoint(); } catch (...) { source_callback_failed=true; throw; }
            },&result,&source_work,containment_cache);
    } catch (const ResourceLimit& error) {
        if (source_callback_failed || !prepared_work || !remaining) throw;
        throw ResourceLimit("Native ownership source budget context"+
            prepared_work->context(limits.max_total_work,*remaining)+"; "+error.what());
    }
    charge_total(source_work);
    if (prepared_work) prepared_work->source+=source_work;
    if (limits.observe_masks) limits.observe_masks(z, false, ownership.regions);
    if (limits.interface_variant != NativeInterfaceVariant::Strict) {
        if (!sampled || retain_unique)
            throw std::invalid_argument("Interface sensitivity requires authenticated unchanged sampler geometry.");
        for (const auto& region : ownership.regions) charge(region.footprint.native());
        std::size_t consumed = 0, native_winding_consumed = 0;
        const auto local_remaining=limits.max_polygon_work-work;
        const auto prepared_remaining=remaining ? *remaining : local_remaining;
        const auto available=std::min(local_remaining,prepared_remaining);
        bool callback_failed=false;
        try {
            ownership.partition = NativeOwnershipPartitionFactory::make(*sampled, ownership, limits.interface_variant,
                available,[&] {
                    if (checkpoint) {
                        try { checkpoint(); } catch (...) { callback_failed=true; throw; }
                    }
                },consumed,native_winding_consumed,limits.allow_single_sided_native_boundary,limits.allow_dual_native_boundary,limits.observe_partition_failure);
        } catch (const ResourceLimit& error) {
            if (callback_failed) throw;
            throw ResourceLimit(std::string("Native ownership partition budget context: limiting_pool=")+
                (remaining && prepared_remaining<local_remaining ? "prepared-total" :
                    remaining && prepared_remaining==local_remaining ? "prepared-total+local-polygon" : "local-polygon")+
                ", prepared_total_remaining="+(remaining ? std::to_string(prepared_remaining) : "not-prepared")+
                ", prepared_total_limit="+std::to_string(limits.max_total_work)+
                ", local_polygon_used="+std::to_string(work)+
                ", local_polygon_limit="+std::to_string(limits.max_polygon_work)+
                ", factory_available="+std::to_string(available)+
                (prepared_work && remaining ? prepared_work->context(limits.max_total_work,*remaining) : std::string{})+
                "; "+error.what());
        }
        charge_total(consumed);
        if (prepared_work) {
            prepared_work->partition+=consumed;
            prepared_work->native_winding+=native_winding_consumed;
        }
        ownership.regions.clear();
        return ownership;
    }
    for (auto& region : ownership.regions) {
        charge(region.footprint.native());
        charge(conflict.section_footprint);
        if (retain_unique) {
            region.footprint.native() = intersection_ex(region.footprint.native(), verified_conflicts);
            charge(region.footprint.native());
        }
        region.footprint.native() = intersection_ex(region.footprint.native(), conflict.section_footprint);
        charge(region.footprint.native());
    }
    if (retain_unique) {
        // All source contributors were authenticated above, including sources
        // contributing only unique original material. A foreign expansion may
        // not displace that material, even outside the original multi-material
        // overlap. No originally empty point is assigned here.
        for (const auto& [material, footprint] : material_sections) {
            charge(footprint);
            charge(verified_conflicts);
            auto unique = diff_ex(footprint, verified_conflicts);
            charge(unique);
            charge(conflict.section_footprint);
            unique = intersection_ex(unique, conflict.section_footprint);
            charge(unique);
            if (unique.empty()) continue;
            const auto owner_material = material;
            auto owner = std::find_if(ownership.regions.begin(), ownership.regions.end(),
                [owner_material](const auto& region) { return region.material == owner_material; });
            if (owner == ownership.regions.end()) ownership.regions.push_back({material, std::move(unique)});
            else {
                charge(owner->footprint.native());
                owner->footprint.native() = union_ex(owner->footprint.native(), unique);
                charge(owner->footprint.native());
            }
        }
    }
    if (limits.observe_masks) limits.observe_masks(z, true, ownership.regions);
    return ownership;
}
} // namespace

SectionOwnership native_ownership_for_overlap(const GCodeProcessorResult& result,
    const LayerMaterialOverlap& conflict, BeadSection section, double bed_z,
    const std::vector<std::size_t>& mapping, const ConflictOwnershipLimits& limits,
    const std::function<void()>& checkpoint)
{
    if (limits.deposited_overlap_policy!=DepositedOverlapPolicy::None)
        throw std::invalid_argument("Deposited-order ownership requires a live sampled prefix.");
    auto path = timed_deposition(result, limits.max_moves, limits.max_segments, checkpoint);
    if (path.state == ToolpathState::BudgetExceeded) throw ResourceLimit("Native conflict path budget exceeded.");
    if (path.state != ToolpathState::Complete || !path.unresolved_deposition.empty())
        throw std::invalid_argument("Native conflict ownership requires complete spatial deposition.");
    const auto nominal = resolve_nominal_deposition(result, path, path.total_s, section, {}, checkpoint);
    std::vector<std::size_t> moves;
    moves.reserve(path.deposition.size());
    for (const auto& entry : path.deposition) moves.push_back(entry.move_index);
    return resolve_overlap(result, conflict, nominal.segments, moves, bed_z, mapping, limits, checkpoint);
}

ConflictOwnershipResolver prepare_native_model_ownership(const GCodeProcessorResult& result,
    BeadSection section, double bed_z, const std::vector<std::size_t>& mapping, double end,
    const std::optional<StationaryExtrusionAllowance>& stationary_allowance,
    const ConflictOwnershipLimits& limits, const std::function<void()>& checkpoint, bool retain_unique,
    SampledConflictOwnershipResolver* sampled_resolver, const std::optional<CollapsedEnvelopeAllowance>& collapsed_allowance)
{
    const bool deposited=limits.deposited_overlap_policy==DepositedOverlapPolicy::FirstDepositionWins;
    if ((limits.deposited_overlap_policy!=DepositedOverlapPolicy::None && !deposited) ||
        (deposited && (!sampled_resolver || retain_unique || collapsed_allowance ||
            !limits.geometry.retain_exact_interface_sections || limits.interface_variant!=NativeInterfaceVariant::Strict ||
            limits.source_centerline_boundary_allowance_nm || limits.allow_single_sided_native_boundary || limits.allow_dual_native_boundary)))
        throw std::invalid_argument("Deposited-order ownership requires a separate exact sampled prefix without collapsed envelopes.");
    if (deposited && (!limits.max_polygon_checks || limits.max_polygon_checks>4096 ||
        !limits.max_polygon_work || limits.max_polygon_work>20000000 ||
        !limits.max_source_work || limits.max_source_work>20000000 ||
        std::any_of(mapping.begin(),mapping.end(),[](std::size_t material) { return material>=256; })))
        throw std::invalid_argument("Invalid deposited ownership material or work limits.");
    if (limits.source_centerline_boundary_allowance_nm > 7 ||
        ((limits.source_centerline_boundary_allowance_nm || limits.allow_single_sided_native_boundary || limits.allow_dual_native_boundary) &&
            limits.interface_variant == NativeInterfaceVariant::Strict) ||
        (limits.allow_dual_native_boundary && (!limits.geometry.retain_exact_interface_sections || retain_unique)))
        throw std::invalid_argument("Native boundary sensitivity requires an explicit non-strict analysis scenario.");
    if (!result.deposition_provenance || !result.native_ownership || mapping.empty() || mapping.size() > 256 ||
        !limits.max_total_work || limits.max_total_work > 1000000000)
        throw std::invalid_argument("Native model ownership requires matching live provenance and material mapping.");
    constexpr std::size_t workspace = 512u * 1024u * 1024u - sizeof(NativeOwnershipSnapshot::ContainmentCache);
    const std::size_t possible = std::min(result.moves.size(), limits.max_segments);
    if (possible > workspace / (sizeof(TimedDeposition) + sizeof(DepositedSegment) + sizeof(double) + sizeof(std::size_t)))
        throw ResourceLimit("Native model ownership input exceeds its workspace budget.");
    auto path = timed_deposition(result, limits.max_moves, limits.max_segments, checkpoint);
    if (path.state == ToolpathState::BudgetExceeded) throw ResourceLimit("Native model ownership path budget exceeded.");
    if (path.state != ToolpathState::Complete || !std::isfinite(end) || end < 0. || end > path.total_s)
        throw std::invalid_argument("Native model ownership requires a complete selected timeline.");
    const auto excluded = [](const auto& entry) {
        return entry.role == erWipeTower || entry.role == erCustom || entry.role == erNone;
    };
    path.deposition.erase(std::remove_if(path.deposition.begin(), path.deposition.end(), excluded), path.deposition.end());
    path.unresolved_deposition.erase(std::remove_if(path.unresolved_deposition.begin(), path.unresolved_deposition.end(), excluded),
        path.unresolved_deposition.end());
    for (const auto& entry : path.deposition)
        if (entry.start_s < end && entry.role == erMixed)
            throw std::invalid_argument("Native model ownership cannot resolve mixed extrusion roles.");
    if (deposited)
        for (const auto& entry:path.deposition)
            if (entry.start_s<end && entry.end_s>end)
                throw std::invalid_argument("Deposited-order ownership requires completed command children.");
    struct Prepared {
        struct SourceIdentity {
            std::uint32_t source, command;
            unsigned filament;
            ExtrusionRole role;
            float extrusion;
        };
        std::vector<DepositedSegment> nominal;
        std::vector<std::size_t> move_indices;
        std::vector<SourceIdentity> source_identities;
        OwnershipTree tree;
    };
    auto prepared = std::make_shared<Prepared>();
    prepared->nominal = resolve_nominal_deposition(result, path, end, section, stationary_allowance, checkpoint, collapsed_allowance).segments;
    if (deposited) {
        if (prepared->nominal.size()>workspace/(sizeof(DepositedSegment)+sizeof(std::size_t)+sizeof(Prepared::SourceIdentity)))
            throw ResourceLimit("Prepared deposited source identities exceed the workspace budget.");
        prepared->source_identities.reserve(prepared->nominal.size());
    }
    prepared->move_indices.reserve(prepared->nominal.size());
    for (std::size_t i = 0; i < prepared->nominal.size(); ++i) {
        prepared->move_indices.push_back(path.deposition[i].move_index);
        if (deposited) {
            const auto& move=result.moves.at(path.deposition[i].move_index);
            prepared->source_identities.push_back({move.deposition_source_id,move.deposition_command_section_id,
                move.extruder_id,move.extrusion_role,move.delta_extruder});
        }
    }
    path = {};
    const auto registry = result.deposition_provenance;
    const auto snapshot = result.native_ownership;
    const auto containment_cache = snapshot->make_containment_cache();
    const auto commands = result.deposition_sections;
    if (result.moves.size() > limits.max_total_work || prepared->nominal.size() > limits.max_total_work - result.moves.size())
        throw ResourceLimit("Prepared native ownership input work budget exceeded.");
    auto remaining = std::make_shared<std::size_t>(limits.max_total_work - result.moves.size() - prepared->nominal.size());
    auto prepared_work = std::make_shared<PreparedOwnershipWork>();
    struct SourceBox {
        std::size_t index;
        OwnershipTree::BoundingBox bounds;
        OwnershipTree::VectorType center;
        std::size_t* remaining;
        const std::function<void()>* checkpoint;
        void tick() const {
            if (!*remaining) throw ResourceLimit("Native ownership index construction work budget exceeded.");
            --*remaining;
            if ((*remaining & 1023) == 0 && *checkpoint) (*checkpoint)();
        }
        std::size_t idx() const { return index; }
        const OwnershipTree::BoundingBox& bbox() const { tick(); return bounds; }
        const OwnershipTree::VectorType& centroid() const { tick(); return center; }
    };
    std::size_t power = 1;
    while (power < prepared->nominal.size()) power *= 2;
    const std::size_t bytes = (2 * power - 1) * sizeof(OwnershipTree::Node) +
        prepared->nominal.size() * (sizeof(SourceBox) + sizeof(DepositedSegment) + sizeof(std::size_t) +
            (deposited ? sizeof(Prepared::SourceIdentity) : 0));
    if (bytes > workspace)
        throw ResourceLimit("Prepared native ownership index exceeds its workspace budget.");
    {
        std::vector<SourceBox> boxes;
        boxes.reserve(prepared->nominal.size());
        for (std::size_t i = 0; i < prepared->nominal.size(); ++i) {
            if (checkpoint) checkpoint();
            const auto bounds = bead_bounds(prepared->nominal[i]);
            boxes.push_back({i, bounds, bounds.center(), remaining.get(), &checkpoint});
        }
        prepared->tree.build_modify_input(boxes);
    }
    const auto candidates_for = [prepared, limits, checkpoint, remaining, prepared_work](const OwnershipTree::BoundingBox& query) {
        std::vector<std::size_t> candidates;
        AABBTreeIndirect::traverse(prepared->tree,
            [&](const OwnershipTree::Node& node) {
                if (checkpoint) checkpoint();
                if (!*remaining) throw ResourceLimit("Prepared native ownership total work budget exceeded"+
                    prepared_work->context(limits.max_total_work,*remaining)+".");
                --*remaining;
                return node.bbox.intersects(query);
            }, [&](const OwnershipTree::Node& node) {
                if (candidates.size() == limits.max_polygon_checks)
                    throw ResourceLimit("Native ownership indexed contributor budget exceeded.");
                candidates.push_back(node.idx);
                return true;
            });
        std::sort(candidates.begin(), candidates.end());
        return candidates;
    };
    const auto resolve = [&result, prepared, registry, snapshot, commands, mapping, bed_z, limits, checkpoint,
        remaining, retain_unique, candidates_for, containment_cache, prepared_work](const LayerMaterialOverlap& conflict, const SampledSectionConflict* sampled) {
        if (checkpoint) checkpoint();
        if (result.deposition_provenance != registry || result.native_ownership != snapshot || result.deposition_sections != commands)
            throw std::invalid_argument("Prepared native ownership no longer matches its source generation.");
        const bool deposited=limits.deposited_overlap_policy==DepositedOverlapPolicy::FirstDepositionWins;
        if (deposited && !sampled)
            throw std::invalid_argument("Deposited-order ownership cannot accept arbitrary conflict masks.");
        if (!conflict.section_footprint_complete || (!deposited && conflict.section_footprint.empty()) || !std::isfinite(conflict.bottom_m))
            throw std::invalid_argument("Native ownership index requires a complete finite conflict.");
        if (sampled) {
            if (retain_unique || sampled->bed_z_m() != bed_z ||
                sampled->geometry_limits().arc_tolerance_m != limits.geometry.arc_tolerance_m ||
                sampled->geometry_limits().half_native_section_planes != limits.geometry.half_native_section_planes ||
                sampled->geometry_limits().retain_exact_interface_sections != limits.geometry.retain_exact_interface_sections ||
                sampled->segments().size() > limits.max_polygon_checks)
                throw std::invalid_argument("Sampled conflict does not match its prepared geometry context.");
            std::set<std::size_t> present;
            for (const auto& actual : sampled->segments()) {
                if (checkpoint) checkpoint();
                if (!*remaining) throw ResourceLimit("Sampled ownership authentication work exceeded"+
                    prepared_work->context(limits.max_total_work,*remaining)+".");
                --*remaining;
                if (!actual.section_source_index || *actual.section_source_index >= prepared->nominal.size() ||
                    !present.insert(*actual.section_source_index).second)
                    throw std::invalid_argument("Sampled ownership has missing, duplicate or unknown source identity.");
                const auto& expected = prepared->nominal[*actual.section_source_index];
                if (expected.material >= mapping.size() || actual.material != mapping[expected.material] ||
                    (actual.start_m.array() != expected.start_m.array()).any() ||
                    (actual.end_m.array() != expected.end_m.array()).any() || actual.width_m != expected.width_m ||
                    actual.height_m != expected.height_m || actual.end_height_m != expected.end_height_m || actual.section != expected.section)
                    throw std::invalid_argument("Sampled ownership source geometry or material was altered.");
            }
            // Check the COMPLETE tile/plane, not merely beads already known to
            // overlap C. Omitting a third material cannot manufacture evidence.
            if (prepared->tree.nodes().empty()) throw std::invalid_argument("Empty sampled ownership prefix.");
            auto query = prepared->tree.nodes().front().bbox;
            if (sampled->clip()) {
                const auto bounds = get_extents(ExPolygons{*sampled->clip()});
                for (unsigned axis = 0; axis < 2; ++axis) {
                    query.min()[axis] = double(bounds.min[axis]) * 1e-9 - 2e-9;
                    query.max()[axis] = double(bounds.max[axis]) * 1e-9 + 2e-9;
                }
            }
            query.min().z() = query.max().z() = conflict.bottom_m;
            for (const auto index : candidates_for(query)) {
                if (checkpoint) checkpoint();
                // Exact geometry for present IDs was checked above. Only a
                // missing conservative AABB candidate needs an actual cut to
                // distinguish harmless tangency from omitted deposition.
                if (present.find(index) != present.end()) continue;
                const auto section = section_deposition({prepared->nominal[index]}, conflict.bottom_m, bed_z,
                    limits.geometry, checkpoint, nullptr, sampled->clip());
                if (section.boundary_points > *remaining) throw ResourceLimit("Sampled ownership completeness work exceeded"+
                    prepared_work->context(limits.max_total_work,*remaining)+".");
                *remaining -= section.boundary_points;
                if (!section.regions.empty())
                    throw std::invalid_argument("Sampled ownership omitted deposition from its observation prefix.");
            }
        }
        if (deposited) {
            const auto debit=[&](std::size_t amount) {
                if (checkpoint) checkpoint();
                if (amount>*remaining) throw ResourceLimit("Prepared deposited ownership budget exceeded.");
                *remaining-=amount;
            };
            const auto* admitted=sampled->frozen_cut_boundary_points();
            if (admitted && admitted->size()!=sampled->segments().size())
                throw std::invalid_argument("Frozen sampled cut admission does not match its inputs.");
            std::vector<std::size_t> order;
            for (std::size_t i=0;i<sampled->segments().size();++i) order.push_back(i);
            std::sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b) {
                debit(1);
                return prepared->move_indices[*sampled->segments()[a].section_source_index]<
                    prepared->move_indices[*sampled->segments()[b].section_source_index];
            });
            std::vector<DepositedSegment> ordered;
            std::size_t source_work=0;
            auto command_geometry=limits.geometry;
            command_geometry.max_boundary_points=std::min<std::size_t>(command_geometry.max_boundary_points,1024);
            // A caller may choose a stricter prepared construction cap than
            // the sampler. Such requests still reconstruct under that cap.
            if (std::min<std::size_t>(sampled->geometry_limits().max_boundary_points,1024)>
                command_geometry.max_boundary_points) admitted=nullptr;
            for (const auto ordinal:order) {
                debit(1);
                const auto index=*sampled->segments()[ordinal].section_source_index;
                auto bead=prepared->nominal[index];
                bead.material=mapping.at(bead.material);
                if (admitted) {
                    // Full frozen geometry/context and completeness were
                    // checked above. Reuse only the actual sampler's cut
                    // admission; every nonempty source still authenticates.
                    debit(1);
                    const auto count=(*admitted)[ordinal];
                    if (count>command_geometry.max_boundary_points)
                        throw ResourceLimit("Deposited command boundary point budget exceeded.");
                    if (!count) continue;
                } else {
                    const auto cut=section_deposition({bead},conflict.bottom_m,bed_z,command_geometry,checkpoint,nullptr,sampled->clip());
                    debit(cut.boundary_points);
                    if (cut.regions.empty()) continue;
                }
                const auto move_index=prepared->move_indices[index];
                if (move_index>=result.moves.size()) throw std::invalid_argument("Deposited source move disappeared.");
                const auto& move=result.moves[move_index];
                const auto& identity=prepared->source_identities[index];
                if (move.deposition_source_id!=identity.source || move.deposition_command_section_id!=identity.command ||
                    move.extruder_id!=identity.filament || move.extrusion_role!=identity.role || move.delta_extruder!=identity.extrusion)
                    throw std::invalid_argument("Deposited source identity changed after preparation.");
                SourceBead source{bead,move.deposition_source_id,move.extruder_id,move_index};
                const auto available=std::min(limits.max_source_work-source_work,*remaining);
                if (!available) throw ResourceLimit("Deposited source authentication budget exhausted.");
                std::size_t used=0;
                bool callback_failed=false;
                try {
                    used=authenticate_live_deposited_source(result,source,available,[&] {
                        if (checkpoint) try { checkpoint(); } catch (...) { callback_failed=true; throw; }
                    });
                } catch (const std::invalid_argument& error) {
                    if (callback_failed) throw;
                    const auto point=[](const Eigen::Vector3d& p) {
                        return std::to_string(p.x()*1e9)+","+std::to_string(p.y()*1e9)+","+std::to_string(p.z()*1e9);
                    };
                    throw std::invalid_argument(std::string(error.what())+" Deposited source context: move="+
                        std::to_string(move_index)+", source="+std::to_string(move.deposition_source_id)+
                        ", command="+std::to_string(move.deposition_command_section_id)+
                        ", role="+std::to_string(unsigned(move.extrusion_role))+", move_type="+std::to_string(unsigned(move.type))+
                        ", filament="+std::to_string(move.extruder_id)+", material="+std::to_string(bead.material)+
                        ", delta_E_filament_nm="+std::to_string(double(move.delta_extruder)*1e6)+
                        ", section_z_nm="+std::to_string(conflict.bottom_m*1e9)+
                        ", start_nm=("+point(bead.start_m)+"), end_nm=("+point(bead.end_m)+").");
                }
                debit(used);
                source_work+=used;
                prepared_work->source+=used;
                if (ordered.size()==4096)
                    throw ResourceLimit("Deposited command count budget exceeded.");
                ordered.push_back(std::move(bead));
            }
            const auto available=std::min({std::size_t(4000000),limits.max_polygon_work,*remaining});
            if (!available) throw ResourceLimit("Deposited partition budget exhausted.");
            std::size_t consumed=0;
            SectionOwnership ownership{conflict.bottom_m,{}};
            bool callback_failed=false;
            const auto partition_checkpoint=[&] {
                if (checkpoint) try { checkpoint(); } catch (...) { callback_failed=true; throw; }
            };
            try {
                ownership.deposited_partition=DepositedOwnershipPartitionFactory::make(*sampled,ordered,available,partition_checkpoint,consumed);
            } catch (const ResourceLimit& error) {
                if (callback_failed) throw;
                throw ResourceLimit("Deposited partition budget context: local_limit="+std::to_string(available)+
                    prepared_work->context(limits.max_total_work,*remaining)+"; "+error.what());
            }
            debit(consumed);
            prepared_work->partition+=consumed;
            return ownership;
        }
        std::size_t points = 0;
        for (const auto& polygon : conflict.section_footprint) {
            if (checkpoint) checkpoint();
            const auto count = [&](std::size_t n) {
                if (n < 3) throw std::invalid_argument("Invalid native ownership query polygon.");
                if (n > limits.max_polygon_work - points) throw ResourceLimit("Native ownership query point budget exceeded.");
                points += n;
            };
            count(polygon.contour.points.size());
            for (const auto& hole : polygon.holes) count(hole.points.size());
        }
        const auto query = conflict_bounds(conflict);
        const auto candidates = candidates_for(query);
        return resolve_overlap(result, conflict, prepared->nominal, prepared->move_indices, bed_z, mapping, limits,
            checkpoint, remaining.get(), retain_unique, &candidates, sampled, containment_cache.get(), prepared_work.get());
    };
    if (sampled_resolver)
        *sampled_resolver = [resolve](const SampledSectionConflict& sampled) { return resolve(sampled.conflict(), &sampled); };
    return [resolve](const LayerMaterialOverlap& conflict) { return resolve(conflict, nullptr); };
}
} // namespace Slic3r::Thermoelastic
