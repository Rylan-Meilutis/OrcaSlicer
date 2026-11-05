#include "ThermoelasticSourceOwnership.hpp"
#include "libslic3r.h"
#include "GCode/DepositionProvenance.hpp"
#include "GCode/DepositionSection.hpp"
#include "GCode/GCodeProcessor.hpp"
#include "ClipperUtils.hpp"
#include "Polyline.hpp"
#include "ExtrusionEntity.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include "Flow.hpp"
#include "Layer.hpp"
#include "Point.hpp"
#include "Print.hpp"
#include "PrintConfig.hpp"
#include "Surface.hpp"
#include "SurfaceCollection.hpp"
#include "ThermoelasticOwnership.hpp"
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticToolpath.hpp"
#include "ThermoelasticExactSections.hpp"
#include "ExPolygon.hpp"
#include "Polygon.hpp"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
namespace {
void validate(const SourceOwnershipLimits& limits)
{
    if (limits.max_sources == 0 || limits.max_sources > DepositionProvenance::maximum_sources || limits.max_boundary_points == 0 ||
        limits.max_boundary_points > 2000000 || limits.max_work == 0 || limits.max_work > 20000000 ||
        limits.centerline_boundary_allowance_nm > 7)
        throw std::invalid_argument("Invalid source ownership budget.");
}
void count_points(const ExPolygons& polygons, std::size_t& used, const SourceOwnershipLimits& limits)
{
    const auto add = [&](const Polygon& polygon) {
        if (polygon.points.size() > limits.max_boundary_points - used)
            throw ResourceLimit("Source ownership boundary budget exceeded.");
        used += polygon.points.size();
    };
    for (const auto& polygon : polygons) {
        add(polygon.contour);
        for (const auto& hole : polygon.holes)
            add(hole);
    }
}
void authenticate_live_span(const GCodeProcessorResult& result,const SourceBead& bead,bool require_command,
    const std::function<void()>& tick,bool allow_auxiliary=false)
{
    if (!bead.move_index || *bead.move_index>=result.moves.size())
        throw std::invalid_argument("Missing live deposition move.");
    const auto& move=result.moves[*bead.move_index];
    const auto* command=result.deposition_sections ? result.deposition_sections->section(move.deposition_command_section_id) : nullptr;
    const auto& geometry=bead.geometry;
    const double end_height=geometry.end_height_m.value_or(geometry.height_m);
    if (move.internal_only || move.type!=EMoveType::Extrude || move.deposition_source_id!=bead.source_id ||
        move.extruder_id!=bead.logical_filament || !geometry.start_m.allFinite() || !geometry.end_m.allFinite() ||
        geometry.start_m.cwiseAbs().maxCoeff()>10. || geometry.end_m.cwiseAbs().maxCoeff()>10. ||
        !std::isfinite(geometry.width_m) || geometry.width_m<=0. || geometry.width_m>10. ||
        !std::isfinite(geometry.height_m) || geometry.height_m<0. || geometry.height_m>10. ||
        !std::isfinite(end_height) || end_height<0. || end_height>10. || std::max(geometry.height_m,end_height)<=0.)
        throw std::invalid_argument("Source geometry is not the matching live deposition move.");
    if (command) {
        if (command->purpose!=DepositionCommandSection::Purpose::Model &&
            !(allow_auxiliary && command->purpose==DepositionCommandSection::Purpose::AuxiliaryBrim))
            throw std::invalid_argument("Auxiliary emission cannot authenticate a native model source.");
        if (command->filament!=bead.logical_filament || !move.deposition_design_heights || !move.deposition_design_width ||
            geometry.height_m!=(*move.deposition_design_heights)[0]*.001 ||
            end_height!=(*move.deposition_design_heights)[1]*.001 || geometry.width_m!=double(*move.deposition_design_width)*.001)
            throw std::invalid_argument("Source dimensions do not match the live command.");
    } else if (require_command || move.deposition_command_section_id || move.deposition_design_heights || move.deposition_design_width ||
        geometry.width_m!=double(move.width)*.001 || geometry.height_m!=double(move.height)*.001)
        throw std::invalid_argument("Source dimensions have no matching live command.");
    std::size_t previous=*bead.move_index;
    while (previous>0) {
        tick();
        --previous;
        if (!result.moves[previous].internal_only && result.moves[previous].type!=EMoveType::Seam) break;
    }
    if (previous==*bead.move_index || result.moves[previous].internal_only || result.moves[previous].type==EMoveType::Seam)
        throw std::invalid_argument("Source start does not match the emitted command.");
    Eigen::Vector3d actual_from=result.moves[previous].position.cast<double>();
    Eigen::Vector3d actual_to=move.position.cast<double>();
    if (const auto z=authenticated_deposition_z(result,*bead.move_index,previous)) {
        actual_from.z()=(*z)[0];
        actual_to.z()=(*z)[1];
    }
    if ((geometry.start_m-actual_from*.001).squaredNorm()!=0. || (geometry.end_m-actual_to*.001).squaredNorm()!=0.)
        throw std::invalid_argument("Source geometry is not the matching live deposition span.");
    if (command) {
        const std::array<double,3> start{actual_from.x(),actual_from.y(),actual_from.z()};
        const std::array<double,3> end{actual_to.x(),actual_to.y(),actual_to.z()};
        DepositionSectionChain membership(*command);
        if (start!=command->start && !membership.append(command->start,start,bead.logical_filament))
            throw std::invalid_argument("Source child starts outside its authenticated parent command.");
        const auto heights=membership.append(start,end,bead.logical_filament);
        if (!heights || *heights!=*move.deposition_design_heights || float(command->width)!=*move.deposition_design_width)
            throw std::invalid_argument("Source child dimensions do not match its authenticated parent command.");
    }
}
} // namespace

std::size_t authenticate_live_deposited_source(const GCodeProcessorResult& result,const SourceBead& bead,
    std::size_t max_work,const std::function<void()>& checkpoint)
{
    if (!max_work || max_work>20000000) throw std::invalid_argument("Invalid live source authentication budget.");
    std::size_t work=0;
    const auto tick=[&] {
        if (checkpoint) checkpoint();
        if (work==max_work) throw ResourceLimit("Live source authentication work budget exceeded.");
        ++work;
    };
    tick();
    if (!result.deposition_provenance || !result.deposition_provenance->sealed() ||
        !result.deposition_sections || !result.deposition_sections->sealed() ||
        result.deposition_sections->generation()!=result.deposition_provenance->generation() ||
        !bead.move_index || *bead.move_index>=result.moves.size())
        throw std::invalid_argument("Live source authentication requires matching sealed registries.");
    const auto* source=result.deposition_provenance->source(bead.source_id);
    const auto& move=result.moves[*bead.move_index];
    const auto* command=result.deposition_sections->section(move.deposition_command_section_id);
    const bool auxiliary=command && command->purpose==DepositionCommandSection::Purpose::AuxiliaryBrim;
    if (!std::isfinite(move.delta_extruder) || move.delta_extruder<=0. || move.deposition_vertical_extrusion ||
        unsigned(move.extrusion_role)>=unsigned(erCount) ||
        move.extrusion_role==erNone || move.extrusion_role==erCustom || move.extrusion_role==erMixed ||
        move.extrusion_role==erWipeTower ||
        (bead.geometry.start_m.head<2>()-bead.geometry.end_m.head<2>()).squaredNorm()==0.)
        throw std::invalid_argument("Unsupported or unresolved live deposition source.");
    if (auxiliary) {
        // An auxiliary emitted shape is not a CAD region. Keep its source ID
        // zero and require its own typed, planar, constant-height witness.
        if (bead.source_id || move.deposition_source_id || move.extrusion_role!=erBrim ||
            command->transform!=DepositionCommandSection::Transform::Ordinary || command->arc_direction ||
            command->nominal_lower_z || command->start[2]!=command->end[2] || command->start_height<=0. ||
            command->start_height!=command->end_height || bead.geometry.start_m.z()!=bead.geometry.end_m.z() ||
            bead.geometry.start_m.z()!=command->start[2]*.001)
            throw std::invalid_argument("Auxiliary brim lacks a matching planar linear emission witness.");
    } else if (!source || source->logical_filament!=bead.logical_filament ||
        (source->height!=DepositionSource::Height::NativePlanar && source->height!=DepositionSource::Height::Transformed) ||
        (command && command->purpose!=DepositionCommandSection::Purpose::Model))
        throw std::invalid_argument("Unsupported or unresolved live deposition source.");
    if (!move.deposition_command_section_id &&
        bead.geometry.end_height_m.value_or(bead.geometry.height_m)!=bead.geometry.height_m)
        throw std::invalid_argument("Untyped live deposition cannot authenticate varying height.");
    authenticate_live_span(result,bead,auxiliary || source->height==DepositionSource::Height::Transformed,tick,auxiliary);
    return work;
}

NativeOwnershipSnapshot capture_native_ownership(const Print& print,
                                                 std::shared_ptr<const DepositionProvenance> registry,
                                                 const Point& parsed_xy_offset,
                                                 double parsed_z_offset_m,
                                                 const SourceOwnershipLimits& limits,
                                                 const std::function<void()>& checkpoint)
{
    validate(limits);
    if (!registry || !registry->sealed() || !std::isfinite(parsed_z_offset_m) || std::abs(parsed_z_offset_m) > 10. ||
        print.config().belt_printer.value || std::abs(double(parsed_xy_offset.x())) > 1e10 || std::abs(double(parsed_xy_offset.y())) > 1e10)
        throw std::invalid_argument("Unsupported source ownership registry or coordinate frame.");
    NativeOwnershipSnapshot out;
    out.m_registry   = std::move(registry);
    std::size_t work = 0, points = 0;
    using ReferenceKey = std::array<std::uint64_t, 5>;
    std::map<ReferenceKey, std::shared_ptr<const NativeOwnershipSnapshot::Reference>> references;
    std::vector<std::size_t> identity_materials(256);
    std::iota(identity_materials.begin(), identity_materials.end(), 0);
    bool callback_failed = false;
    const auto tick      = [&] {
        if (checkpoint) {
            try {
                checkpoint();
            } catch (...) {
                callback_failed = true;
                throw;
            }
        }
        if (++work > limits.max_work)
            throw ResourceLimit("Source ownership work budget exceeded.");
    };
    for (std::uint32_t id = 1; const auto* source = out.m_registry->source(id); ++id) {
        tick();
        if (out.m_masks.size() == limits.max_sources)
            throw ResourceLimit("Source ownership source budget exceeded.");
        if (source->logical_filament >= identity_materials.size())
            throw std::invalid_argument("Source ownership logical filament is out of range.");
        const PrintObject* object = nullptr;
        for (const auto* candidate : print.objects()) {
            tick();
            if (candidate->get_id() == source->object) {
                if (object)
                    throw std::invalid_argument("Ambiguous native source object.");
                object = candidate;
            }
        }
        if (!object || source->instance >= object->instances().size())
            throw std::invalid_argument("Missing native source object or instance.");
        const Layer* layer = nullptr;
        for (const auto* candidate : object->layers()) {
            tick();
            if (candidate->id() == source->layer) {
                if (layer)
                    throw std::invalid_argument("Ambiguous native source layer.");
                layer = candidate;
            }
        }
        if (!layer)
            throw std::invalid_argument("Missing native source layer.");
        const LayerRegion* region = nullptr;
        for (const auto* candidate : layer->regions()) {
            tick();
            if (candidate->region().print_region_id() == source->region) {
                if (region)
                    throw std::invalid_argument("Ambiguous native source region.");
                region = candidate;
            }
        }
        if (!region)
            throw std::invalid_argument("Missing native source region.");
        const auto& object_config=object->config();
        const bool purge_enabled=object_config.flush_into_objects.value || object_config.flush_into_infill.value;
        const bool nominal_assignment=source->assignment==DepositionSource::Assignment::NominalRole;
        if (source->height == DepositionSource::Height::Unknown ||
            source->assignment==DepositionSource::Assignment::WipingOverride ||
            (purge_enabled && !nominal_assignment)) {
            out.m_masks.emplace_back(std::nullopt);
            continue;
        }
        const auto& shift = object->instances()[source->instance].shift;
        if (std::abs(double(shift.x()) + double(parsed_xy_offset.x())) > 1e10 ||
            std::abs(double(shift.y()) + double(parsed_xy_offset.y())) > 1e10)
            throw std::invalid_argument("Native instance placement exceeds coordinate bounds.");
        const Point translation        = shift + parsed_xy_offset;
        const auto configured_filament = region->region().extruder(frPerimeter);
        if (configured_filament == 0 || configured_filament - 1 != source->logical_filament) {
            out.m_masks.emplace_back(std::nullopt);
            continue;
        }
        // Reuse all existing planar/role/mixed/purge guards rather than treating
        // a source ID or the NativePlanar enum as a physical certificate.
        LayerGeometryLimits geometry_limits;
        geometry_limits.max_boundary_points = limits.max_boundary_points;
        try {
            if (source->height == DepositionSource::Height::NativePlanar) {
                // Only native_layer_ownership examines all regions. A
                // transformed source uses the selected immutable references
                // below, whose actual copies are separately charged once.
                for (const auto* candidate : layer->regions())
                    for (const auto& surface : candidate->slices.surfaces) {
                        std::size_t count = surface.expolygon.contour.points.size();
                        for (const auto& hole : surface.expolygon.holes)
                            count += hole.points.size();
                        if (count > limits.max_work - work)
                            throw ResourceLimit("Source ownership polygon work budget exceeded.");
                        work += count;
                    }
                native_layer_ownership(*layer, .001 * (layer->print_z - .5 * layer->height) + parsed_z_offset_m, translation,
                                       parsed_z_offset_m, identity_materials, geometry_limits, tick);
            }
        } catch (const std::invalid_argument&) {
            if (callback_failed)
                throw;
            out.m_masks.emplace_back(std::nullopt);
            continue;
        }
        NativeOwnershipSnapshot::Mask mask{source->logical_filament,
                                           .001 * layer->bottom_z() + parsed_z_offset_m,
                                           .001 * layer->print_z + parsed_z_offset_m,
                                           {},
                                           source->height == DepositionSource::Height::NativePlanar};
        const auto& region_config=region->region().config();
        const auto& print_config=print.config();
        const auto filament=source->logical_filament;
        const bool ordinary_settings=!purge_enabled && !region_config.zaa_enabled.value &&
            !region_config.infill_combination.value && region_config.top_surface_z_mode.value==TopSurfaceZMode::Disabled;
        if (!print_config.spiral_mode.value && !region_config.nonplanar_infill.value &&
            (ordinary_settings || nominal_assignment) &&
            (print_config.filament_is_mixed.values.empty() || !print_config.filament_is_mixed.get_at(filament)) &&
            (print_config.filament_mixed_gradient.values.empty() || !print_config.filament_mixed_gradient.get_at(filament))) {
            mask.infill_requires_command=!ordinary_settings;
            for (const auto& [role,flow] : {std::pair{erInternalInfill,frInfill},
                                          std::pair{erSolidInfill,frSolidInfill},
                                          std::pair{erTopSolidInfill,frTopSolidInfill}})
                if (region->region().extruder(flow)==filament+1)
                    mask.planar_infill_roles|=std::uint32_t(1)<<unsigned(role);
        }
        const auto capture_reference = [&](const Layer& reference_layer,const LayerRegion& reference_region) {
            tick();
            const ReferenceKey key{source->object,source->instance,std::uint64_t(reference_layer.id()),source->region,filament};
            const auto found=references.find(key);
            if (found!=references.end()) return found->second;
            if (references.size()==3*limits.max_sources)
                throw ResourceLimit("Source ownership reference budget exceeded.");
            auto reference=std::make_shared<NativeOwnershipSnapshot::Reference>();
            reference->bottom_m=.001*reference_layer.bottom_z()+parsed_z_offset_m;
            reference->top_m=.001*reference_layer.print_z+parsed_z_offset_m;
            if (!std::isfinite(reference->bottom_m) || !std::isfinite(reference->top_m) ||
                std::abs(reference->bottom_m)>10. || std::abs(reference->top_m)>10. || reference->top_m<=reference->bottom_m)
                throw std::invalid_argument("Invalid native source reference interval.");
            reference->identity=references.size();
            for (const auto& surface:reference_region.slices.surfaces) {
                tick();
                const auto before=points;
                count_points(ExPolygons{surface.expolygon},points,limits);
                if (points-before>limits.max_work-work)
                    throw ResourceLimit("Source ownership reference copy work exceeded.");
                work+=points-before;
                reference->polygons.push_back(surface.expolygon);
                reference->polygons.back().translate(translation);
            }
            std::shared_ptr<const NativeOwnershipSnapshot::Reference> immutable=std::move(reference);
            references.emplace(key,immutable);
            return immutable;
        };
        mask.reference=capture_reference(*layer,*region);
        mask.brick_eligible=!mask.native_planar && layer->lower_layer && !mask.reference->polygons.empty() &&
            region_config.perimeter_layering.value==PerimeterLayeringMode::Brick &&
            region->region().extruder(frExternalPerimeter)==filament+1 &&
            !print_config.spiral_mode.value && !region_config.nonplanar_infill.value &&
            (ordinary_settings || nominal_assignment) &&
            (print_config.filament_is_mixed.values.empty() || !print_config.filament_is_mixed.get_at(filament)) &&
            (print_config.filament_mixed_gradient.values.empty() || !print_config.filament_mixed_gradient.get_at(filament));
        if (mask.brick_eligible) {
            const std::array<const Layer*,2> adjacent{layer->lower_layer,layer->upper_layer};
            for (std::size_t side=0;side<adjacent.size();++side) {
                const auto* neighbor=adjacent[side];
                if (!neighbor || neighbor->object()!=object) continue;
                const double contact=side ? neighbor->bottom_z() : neighbor->print_z;
                const double expected=side ? layer->print_z : layer->bottom_z();
                if (!std::isfinite(contact) || !std::isfinite(expected))
                    throw std::invalid_argument("Invalid brick reference adjacency.");
                if (std::llround(contact*1e6)!=std::llround(expected*1e6)) continue;
                const LayerRegion* matching=nullptr;
                for (const auto* candidate:neighbor->regions()) {
                    tick();
                    if (candidate->region().print_region_id()!=source->region) continue;
                    if (matching) throw std::invalid_argument("Ambiguous brick reference region.");
                    matching=candidate;
                }
                if (!matching || matching->slices.surfaces.empty() || matching->region().extruder(frPerimeter)!=filament+1 ||
                    matching->region().extruder(frExternalPerimeter)!=filament+1)
                    continue;
                mask.brick_neighbors[side]=capture_reference(*neighbor,*matching);
            }
        }
        out.m_masks.push_back(std::move(mask));
    }
    return out;
}

SectionOwnership NativeOwnershipSnapshot::section(std::shared_ptr<const DepositionProvenance> registry,
                                                  const std::vector<SourceBead>& active_beads,
                                                  double z_m,
                                                  const std::vector<std::size_t>& filament_to_material,
                                                  const SourceOwnershipLimits& limits,
                                                  const std::function<void()>& checkpoint,
                                                  const GCodeProcessorResult* emitted_evidence, std::size_t* used_work,
                                                  ContainmentCache* containment_cache) const
{
    if (used_work) *used_work = 0;
    validate(limits);
    if (containment_cache && containment_cache->m_identity != m_identity)
        throw std::invalid_argument("Source containment cache belongs to a different snapshot.");
    if (registry != m_registry || !std::isfinite(z_m) || std::abs(z_m) > 10. || filament_to_material.empty() ||
        filament_to_material.size() > 256)
        throw std::invalid_argument("Source ownership registry or section mismatch.");
    if (active_beads.size() > limits.max_work)
        throw ResourceLimit("Source ownership bead budget exceeded.");
    SectionOwnership result{z_m, {}};
    std::set<std::uint32_t> included;
    std::size_t points = 0, work = active_beads.size();
    for (const auto& bead : active_beads) {
        if (checkpoint)
            checkpoint();
        if (emitted_evidence && emitted_evidence->deposition_sections && bead.move_index &&
            *bead.move_index<emitted_evidence->moves.size()) {
            const auto* command=emitted_evidence->deposition_sections->section(
                emitted_evidence->moves[*bead.move_index].deposition_command_section_id);
            if (command && command->purpose!=DepositionCommandSection::Purpose::Model)
                throw std::invalid_argument("Auxiliary emission cannot establish native material ownership.");
        }
        if (bead.source_id == 0 || bead.source_id > m_masks.size())
            throw std::invalid_argument("Unowned deposition cannot establish native material ownership.");
        if (!m_masks[bead.source_id - 1])
            throw std::invalid_argument("Deposition source has unsupported native ownership.");
        const auto& mask = *m_masks[bead.source_id - 1];
        auto reference=mask.reference;
        if (mask.logical_filament >= filament_to_material.size())
            throw std::invalid_argument("Missing logical filament material mapping.");
        const auto material     = filament_to_material[mask.logical_filament];
        const auto* source      = m_registry->source(bead.source_id);
        const auto& geometry    = bead.geometry;
        const double end_height = geometry.end_height_m.value_or(geometry.height_m);
        // Strict nanometre comparison: no tolerance which could hide a scarf or
        // brick displacement. Float-preview rounding may conservatively reject.
        const auto same = [](double a, double b) { return std::llround(a * 1e9) == std::llround(b * 1e9); };
        const auto interval_error = [&](const char* stage) {
            // Fixed scalar fields only: do not dump the whole path or mask on
            // failure. Diagnostics never change interval/ownership decisions.
            const auto point = [](const Eigen::Vector3d& p) {
                return std::to_string(p.x()*1e9)+","+std::to_string(p.y()*1e9)+","+std::to_string(p.z()*1e9);
            };
            std::string detail="Deposition does not match its native material layer interval: stage="+std::string(stage)+
                ", source="+std::to_string(bead.source_id)+", layer="+std::to_string(source->layer)+
                ", move="+(bead.move_index ? std::to_string(*bead.move_index) : std::string("none"))+
                ", material="+std::to_string(geometry.material)+", expected_material="+std::to_string(material)+
                ", filament="+std::to_string(bead.logical_filament)+", expected_filament="+std::to_string(source->logical_filament)+
                ", query_z_nm="+std::to_string(z_m*1e9)+", bottom_nm="+std::to_string(mask.bottom_m*1e9)+
                ", top_nm="+std::to_string(mask.top_m*1e9)+", start_nm=("+point(geometry.start_m)+
                "), end_nm=("+point(geometry.end_m)+"), width_nm="+std::to_string(geometry.width_m*1e9)+
                ", height_nm="+std::to_string(geometry.height_m*1e9)+", end_height_nm="+std::to_string(end_height*1e9)+
                ", native_planar="+std::to_string(mask.native_planar);
            if (emitted_evidence && bead.move_index && *bead.move_index<emitted_evidence->moves.size()) {
                const auto& move=emitted_evidence->moves[*bead.move_index];
                detail+=", role="+std::to_string(unsigned(move.extrusion_role))+
                    ", command="+std::to_string(move.deposition_command_section_id);
            }
            return std::invalid_argument(detail+".");
        };
        if (!geometry.start_m.allFinite() || !geometry.end_m.allFinite() || geometry.start_m.cwiseAbs().maxCoeff() > 10. ||
            geometry.end_m.cwiseAbs().maxCoeff() > 10. || !std::isfinite(geometry.width_m) || geometry.width_m <= 0. ||
            geometry.width_m > 10. || !std::isfinite(geometry.height_m) || geometry.height_m < 0. || geometry.height_m > 10. ||
            !std::isfinite(end_height) || end_height < 0. || end_height > 10. || std::max(geometry.height_m, end_height) <= 0. ||
            geometry.material != material || bead.logical_filament != source->logical_filament)
            throw interval_error("input or query");
        if (mask.native_planar) {
            if (z_m<mask.bottom_m || z_m>=mask.top_m) throw interval_error("input or query");
            if (!same(end_height, geometry.height_m) || !same(geometry.start_m.z(), mask.top_m) || !same(geometry.end_m.z(), mask.top_m) ||
                !same(geometry.height_m, mask.top_m - mask.bottom_m))
                throw interval_error("native planar envelope");
        } else {
            // Region-wide transformed provenance also covers ordinary infill
            // when scarf is enabled on its walls. Those moves need independent
            // live planar evidence; the source enum alone proves neither case.
            if (!emitted_evidence || emitted_evidence->native_ownership.get() != this ||
                emitted_evidence->deposition_provenance != m_registry || !emitted_evidence->deposition_sections ||
                !emitted_evidence->deposition_sections->sealed() ||
                emitted_evidence->deposition_sections->generation() != m_registry->generation() || !bead.move_index ||
                *bead.move_index >= emitted_evidence->moves.size())
                throw std::invalid_argument("Transformed ownership requires matching live scarf evidence.");
            const auto& move    = emitted_evidence->moves[*bead.move_index];
            const auto* command = emitted_evidence->deposition_sections->section(move.deposition_command_section_id);
            const bool brick=command && command->transform==DepositionCommandSection::Transform::BrickCourse;
            if (brick && !mask.brick_eligible)
                throw std::invalid_argument("Brick command has no eligible captured material references.");
            if (brick && (command->arc_direction!=0 || command->nominal_lower_z || command->start_height!=command->end_height))
                throw std::invalid_argument("Brick ownership requires an unmixed constant-height linear command.");
            if (!brick && (z_m<mask.bottom_m || z_m>=mask.top_m)) throw interval_error("input or query");
            const bool native_plane_infill=!brick && unsigned(move.extrusion_role)<32 &&
                (mask.planar_infill_roles & (std::uint32_t(1)<<unsigned(move.extrusion_role))) &&
                same(end_height,geometry.height_m) && same(geometry.start_m.z(),mask.top_m) &&
                same(geometry.end_m.z(),mask.top_m) && same(geometry.height_m,mask.top_m-mask.bottom_m);
            const bool planar_infill=native_plane_infill && !mask.infill_requires_command &&
                move.deposition_command_section_id==0 && !move.deposition_design_heights && !move.deposition_design_width &&
                geometry.width_m==double(move.width)*.001 && geometry.height_m==double(move.height)*.001;
            if (move.internal_only || move.type != EMoveType::Extrude || move.deposition_source_id != bead.source_id ||
                move.extruder_id != bead.logical_filament)
                throw std::invalid_argument("Source geometry is not the matching live deposition move.");
            const bool supported_role=brick ?
                (move.extrusion_role==erStaggeredPerimeter || move.extrusion_role==erPerimeter || move.extrusion_role==erExternalPerimeter) :
                (move.extrusion_role==erPerimeter || move.extrusion_role==erExternalPerimeter || native_plane_infill);
            if (!planar_infill && (!command || command->filament != bead.logical_filament ||
                !supported_role || !move.deposition_design_heights ||
                !move.deposition_design_width || geometry.height_m != double((*move.deposition_design_heights)[0]) * .001 ||
                end_height != double((*move.deposition_design_heights)[1]) * .001 ||
                geometry.width_m != double(*move.deposition_design_width) * .001))
                throw std::invalid_argument("Scarf geometry is not the authenticated emitted command (move=" +
                    std::to_string(*bead.move_index) + ", source=" + std::to_string(bead.source_id) +
                    ", command=" + std::to_string(move.deposition_command_section_id) + ", role=" +
                    std::to_string(unsigned(move.extrusion_role)) + ", planar_roles=" +
                    std::to_string(mask.planar_infill_roles) + ", native_z_nm=" + std::to_string(mask.top_m * 1e9) +
                    ", native_height_nm=" + std::to_string((mask.top_m - mask.bottom_m) * 1e9) +
                    ", start_z_nm=" + std::to_string(geometry.start_m.z() * 1e9) +
                    ", end_z_nm=" + std::to_string(geometry.end_m.z() * 1e9) +
                    ", height_nm=" + std::to_string(geometry.height_m * 1e9) +
                    ", end_height_nm=" + std::to_string(end_height * 1e9) +
                    ", width_nm=" + std::to_string(geometry.width_m * 1e9) +
                    ", move_height_nm=" + std::to_string(double(move.height) * 1e6) +
                    ", move_width_nm=" + std::to_string(double(move.width) * 1e6) +
                    ", design_heights=" + std::to_string(bool(move.deposition_design_heights)) +
                    ", design_width=" + std::to_string(bool(move.deposition_design_width)) + ").");
            authenticate_live_span(*emitted_evidence,bead,!planar_infill,[&] {
                if (checkpoint) checkpoint();
                if (work==limits.max_work) throw ResourceLimit("Scarf ownership work budget exceeded.");
                ++work;
            });
            const auto native = [](double v) -> coord_t { return std::llround(v * 1e9); };
            double envelope_bottom=mask.bottom_m,envelope_top=mask.top_m;
            if (brick) {
                // The command authenticates the emitted transform, not a
                // vertically extended original mask. Select the actual native
                // interval only after all live parent/child checks above.
                if (mask.brick_neighbors[0]) envelope_bottom=mask.brick_neighbors[0]->bottom_m;
                if (mask.brick_neighbors[1]) envelope_top=mask.brick_neighbors[1]->top_m;
                std::shared_ptr<const Reference> selected;
                for (const auto& candidate:{mask.reference,mask.brick_neighbors[0],mask.brick_neighbors[1]}) {
                    if (work==limits.max_work) throw ResourceLimit("Brick reference selection work exceeded.");
                    ++work;
                    if (checkpoint) checkpoint();
                    if (!candidate || z_m<candidate->bottom_m || z_m>=candidate->top_m) continue;
                    if (selected) throw std::invalid_argument("Brick query has ambiguous native intervals.");
                    selected=candidate;
                }
                if (!selected) throw interval_error("brick query reference");
                reference=std::move(selected);
            }
            if (native(geometry.start_m.z()) > native(envelope_top) || native(geometry.end_m.z()) > native(envelope_top) ||
                native(geometry.start_m.z() - geometry.height_m) < native(envelope_bottom) ||
                native(geometry.end_m.z() - end_height) < native(envelope_bottom))
                throw std::invalid_argument("Scarf envelope crosses its native material layer interval.");
            const ContainmentCache::Key containment_key{bead.source_id,std::int64_t(reference->identity),
                native(geometry.start_m.x()), native(geometry.start_m.y()),
                native(geometry.end_m.x()), native(geometry.end_m.y()), limits.centerline_boundary_allowance_nm};
            const auto cache_tick = [&] {
                if (work == limits.max_work) throw ResourceLimit("Source containment cache work budget exceeded.");
                ++work;
                if (checkpoint) checkpoint();
            };
            bool cached = false;
            std::optional<ContainmentCache::Key>* cache_entry = nullptr;
            if (containment_cache) {
                cache_tick();
                cache_entry = &containment_cache->m_covered[ContainmentCache::slot(containment_key)];
                cached = *cache_entry && **cache_entry == containment_key;
            }
            // Only the immutable-mask XY predicate is reusable. All live source,
            // command, span, geometry and query-Z checks above run on every call.
            if (!cached) {
                std::size_t native_points = 0;
                for (const auto& polygon : reference->polygons) {
                    const auto charge = [&](const Polygon& ring) {
                        if (ring.points.size() > limits.max_work - work)
                            throw ResourceLimit("Scarf ownership polygon work exceeded.");
                        work += ring.points.size();
                        native_points += ring.points.size();
                    };
                    charge(polygon.contour);
                    for (const auto& hole : polygon.holes)
                        charge(hole);
                }
                const Polyline centerline(Point(native(geometry.start_m.x()), native(geometry.start_m.y())),
                                          Point(native(geometry.end_m.x()), native(geometry.end_m.y())));
                const auto outside = diff_pl(Polylines{centerline}, reference->polygons);
                bool exactly_covered = outside.empty();
                bool within_boundary_band = false;
                const unsigned band_nm = limits.centerline_boundary_allowance_nm ? limits.centerline_boundary_allowance_nm : 4;
                if (!exactly_covered) {
                    if (work == limits.max_work) throw ResourceLimit("Exact source containment work budget exceeded.");
                    ExactSegmentCoverage exact;
                    try {
                        exact = exact_native_section_covers_segment(reference->polygons,centerline.points.front(),centerline.points.back(),
                            std::min<std::size_t>(4000000,limits.max_work-work),checkpoint,band_nm);
                    } catch (const ResourceLimit& error) {
                        throw ResourceLimit("Native centerline predicate resource limit (move=" + std::to_string(*bead.move_index) +
                            ", source=" + std::to_string(bead.source_id) + ", command=" +
                            std::to_string(move.deposition_command_section_id) + ", native_points=" + std::to_string(native_points) +
                            ", start_nm=" + std::to_string(centerline.points.front().x()) + "," +
                            std::to_string(centerline.points.front().y()) + ", end_nm=" +
                            std::to_string(centerline.points.back().x()) + "," + std::to_string(centerline.points.back().y()) +
                            ", band_nm=" + std::to_string(band_nm) + "): " + error.what());
                    } catch (const std::invalid_argument& error) {
                        throw std::invalid_argument("Native centerline predicate failed (move=" + std::to_string(*bead.move_index) +
                            ", source=" + std::to_string(bead.source_id) + ", command=" +
                            std::to_string(move.deposition_command_section_id) + "): " + error.what());
                    }
                    work += exact.work;
                    exactly_covered = exact.covered;
                    within_boundary_band = exact.within_boundary_band;
                }
                if (!exactly_covered && !(limits.centerline_boundary_allowance_nm && within_boundary_band)) {
                    double length = 0.;
                    for (const auto& path : outside) length += path.length();
                    throw std::invalid_argument("Command centerline leaves its native material region (move=" +
                        std::to_string(*bead.move_index) + ", source=" + std::to_string(bead.source_id) +
                        ", command=" + std::to_string(move.deposition_command_section_id) +
                        ", arc=" + std::to_string(command ? command->arc_direction : 0) + ", outside_nm=" + std::to_string(length) +
                        ", boundary_band_nm=" + std::to_string(band_nm) +
                        ", within_boundary_band=" + (within_boundary_band ? "true" : "false") +
                        ", start_nm=" + std::to_string(centerline.points.front().x()) + "," +
                        std::to_string(centerline.points.front().y()) + ", end_nm=" +
                        std::to_string(centerline.points.back().x()) + "," + std::to_string(centerline.points.back().y()) + ").");
                }
                if (containment_cache) {
                    cache_tick();
                    *cache_entry = containment_key;
                }
            }
        }
        if (included.insert(bead.source_id).second) {
            if (included.size() > limits.max_sources)
                throw ResourceLimit("Source ownership section budget exceeded.");
            count_points(reference->polygons, points, limits);
            result.regions.push_back({material, reference->polygons});
        }
    }
    if (used_work) *used_work = work;
    return result;
}
} // namespace Slic3r::Thermoelastic
