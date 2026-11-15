#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ThermoelasticConflictOwnership.hpp"
#include "libslic3r/SectionGCodeHistory.hpp"
#include "libslic3r/SectionGrowthStream.hpp"
#include "libslic3r/SectionThermalAssembly.hpp"
#include "libslic3r/ThermalHistory.hpp"
#include <Eigen/Core>
#include "libslic3r/ThermoelasticSourceOwnership.hpp"
#include "libslic3r/ThermoelasticToolpath.hpp"
#include "libslic3r/ThermoelasticLayers.hpp"
#include "libslic3r/ThermoelasticLocalSections.hpp"
#include "libslic3r/ThermoelasticLocalTopology.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticExactSections.hpp"
#include "libslic3r/ThermoelasticSectionFootprint.hpp"
#include "libslic3r/ThermoelasticSectionTopology.hpp"
#include "libslic3r/GCode/DepositionProvenance.hpp"
#include "libslic3r/GCode/DepositionSection.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Surface.hpp"
#include "libslic3r/SurfaceCollection.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/libslic3r.h"
#include "test_helpers.hpp"
#include <cstddef>
#include <algorithm>
#include <optional>
#include <array>
#include <cstdint>
#include <cmath>
#include <initializer_list>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;

namespace {
void populate(GCodeProcessorResult& result, bool transformed = false, coord_t native_mask_shift = 0, unsigned filaments = 2)
{
    const auto config = Test::multifilament_config(filaments, {{"layer_height", "0.2"},
                                                       {"initial_layer_print_height", "0.2"},
                                                       {"thick_bridges", "0"},
                                                       {"thick_internal_bridges", "0"},
                                                       {"elefant_foot_compensation", "0"},
                                                       {"brim_type", "no_brim"},
                                                       {"skirt_loops", "0"}});
    Print print;
    Model model;
    Test::init_print({Test::cube(2.)}, print, model, config);
    auto& object = *model.objects.front();
    object.volumes.front()->config.set("extruder", 1);
    auto second = Test::cube(2.);
    second.translate(2., 0., 0.);
    object.add_volume(std::move(second))->config.set("extruder", 2);
    if (filaments==3) {
        auto third=Test::cube(2.);
        third.translate(4.,0.,0.);
        object.add_volume(std::move(third))->config.set("extruder",3);
    }
    print.apply(model, config);
    print.process();
    const auto* printed = print.get_object(0);
    const auto* layer   = printed->layers().front();
    auto registry       = std::make_shared<DepositionProvenance>();
    std::vector<BoundingBox> boxes(filaments);
    std::vector<std::uint32_t> ids(filaments);
    for (const auto* region : layer->regions()) {
        const auto filament = region->region().extruder(frPerimeter) - 1;
        REQUIRE(filament < filaments);
        ids[filament] = registry->add({printed->get_id(), 0, std::uint64_t(layer->id()), std::uint64_t(region->region().print_region_id()),
                                       filament,
                                       transformed ? DepositionSource::Height::Transformed : DepositionSource::Height::NativePlanar});
        ExPolygons polygons;
        for (const auto& surface : region->slices.surfaces) {
            polygons.push_back(surface.expolygon);
            polygons.back().translate(printed->instances()[0].shift);
        }
        boxes[filament] = get_extents(polygons);
    }
    // An explicit reference-mask discrepancy fixture; bead placement below
    // still uses the original boxes. Production capture itself is unchanged.
    if (native_mask_shift)
        for (auto* region:print.get_object(0)->layers().front()->regions()) {
            const auto filament=region->region().extruder(frPerimeter)-1;
            for (auto& surface:region->slices.surfaces)
                surface.expolygon.translate(Point(filament ? native_mask_shift : -native_mask_shift,0));
        }
    registry->seal();
    result.deposition_provenance = registry;
    result.native_ownership = std::make_shared<const NativeOwnershipSnapshot>(capture_native_ownership(print, registry, Point(0, 0), 0.));
    REQUIRE(boxes[0].max.x() == boxes[1].min.x());
    const float x   = float(double(boxes[0].max.x()) * 1e-6);
    const float y   = float(.5 * double(boxes[0].min.y() + boxes[0].max.y()) * 1e-6);
    const auto line = [&](float px, unsigned filament, std::uint32_t id) {
        GCodeProcessorResult::MoveVertex move;
        move.type        = EMoveType::Travel;
        move.position    = {px, y - .5f, float(layer->print_z)};
        move.time[0]     = 1.f;
        move.extruder_id = filament;
        result.moves.push_back(move);
        move.type                 = EMoveType::Extrude;
        move.position.y()         = y + .5f;
        move.width                = .4f;
        move.height               = float(layer->height);
        move.extrusion_role       = id ? erExternalPerimeter : erSkirt;
        move.deposition_source_id = id;
        move.delta_extruder       = .1f;
        result.moves.push_back(move);
    };
    line(x - .1f, 0, ids[0]);
    line(x + .1f, 1, ids[1]);
    if (filaments==3) line(x+.3f,2,ids[2]);
    line(x + 20.f, 0, 0); // Unowned brim, spatially unrelated to the conflict.
    result.print_statistics.modes[0].time = double(result.moves.size());
}

TEST_CASE("Full native partitions bind one live sampler and preserve strict callback behavior", "[ThermoelasticConflictOwnership]")
{
    STATIC_REQUIRE_FALSE(std::is_default_constructible_v<AuthenticatedSectionPartition>);
    GCodeProcessorResult result;
    populate(result);
    const auto path = timed_deposition(result);
    auto nominal = resolve_nominal_deposition(result, path, path.total_s, BeadSection::RectangularEnvelope);
    for (std::size_t i = 0; i < nominal.segments.size(); ++i) nominal.segments[i].section_source_index = i;
    ConflictOwnershipLimits limits;
    limits.geometry.ownership_boundary_grid_units = 4;
    limits.max_polygon_work = 4000000;
    limits.interface_variant = NativeInterfaceVariant::FavorFirst;
    SampledConflictOwnershipResolver sampled;
    prepare_native_model_ownership(result, BeadSection::RectangularEnvelope, 0., {0,1}, path.total_s,
        {}, limits, {}, false, &sampled);
    SectionOwnership owner{.0001, {}};
    std::shared_ptr<const AuthenticatedSectionPartition> captured;
    owner.resolve_partition = [&](const SampledSectionConflict& evidence) {
        captured = sampled(evidence).partition;
        return captured;
    };
    const auto section = section_deposition(nominal.segments, .0001, 0., limits.geometry, {}, &owner);
    REQUIRE(captured);
    REQUIRE(section.regions.size() == 2);
    CHECK(section.boundary_error.uncertainty_nm == 7);
    CHECK(exact_section_empty(section.regions[0].footprint.native(), section.regions[1].footprint.native(),
        ExactSectionOperation::Intersection).empty);
    owner.resolve_partition = [&](const SampledSectionConflict&) { return captured; };
    CHECK_THROWS_AS(section_deposition(nominal.segments, .0001, 0., limits.geometry, {}, &owner), std::invalid_argument);
    owner.resolve_partition = {};
    owner.partition = captured;
    CHECK_THROWS_AS(section_deposition(nominal.segments, .0001, 0., limits.geometry, {}, &owner), std::invalid_argument);
    owner.partition.reset();
    owner.resolve_partition = [&](const SampledSectionConflict& evidence) { return sampled(evidence).partition; };
    owner.resolve_sampled_conflict = [&](const SampledSectionConflict& evidence) { return sampled(evidence).regions; };
    CHECK_THROWS_AS(section_deposition(nominal.segments, .0001, 0., limits.geometry, {}, &owner), std::invalid_argument);
    owner.resolve_sampled_conflict = {};
    result.deposition_provenance.reset();
    CHECK_THROWS_AS(section_deposition(nominal.segments, .0001, 0., limits.geometry, {}, &owner), std::invalid_argument);
}

TEST_CASE("Deposited order partitions bind the complete live prefix without native qualification", "[ThermoelasticConflictOwnership]")
{
    const bool hidden_later_material=GENERATE(false,true);
    const bool reversed_material_ids=GENERATE(false,true);
    const std::vector<std::size_t> mapping=reversed_material_ids ? std::vector<std::size_t>{1,0} : std::vector<std::size_t>{0,1};
    GCodeProcessorResult result;
    populate(result);
    if (hidden_later_material) {
        result.moves[2].position=result.moves[0].position;
        result.moves[3].position=result.moves[1].position;
    }
    auto commands=std::make_shared<DepositionSections>(result.deposition_provenance->generation());
    commands->seal();
    result.deposition_sections=commands;
    const auto path=timed_deposition(result);
    auto nominal=resolve_nominal_deposition(result,path,path.total_s,BeadSection::RectangularEnvelope);
    for (std::size_t i=0;i<nominal.segments.size();++i) {
        nominal.segments[i].section_source_index=i;
        nominal.segments[i].material=mapping[nominal.segments[i].material];
    }
    const auto center=.5*(nominal.segments[0].start_m+nominal.segments[0].end_m);
    const coord_t x=std::llround(.5*(nominal.segments[0].start_m.x()+nominal.segments[1].start_m.x())*1e9);
    const coord_t y=std::llround(center.y()*1e9);
    const ExPolygon clip(Polygon{{x-400000,y-100000},{x+400000,y-100000},
        {x+400000,y+100000},{x-400000,y+100000}});
    ConflictOwnershipLimits limits;
    CHECK(limits.deposited_overlap_policy==DepositedOverlapPolicy::None);
    limits.deposited_overlap_policy=DepositedOverlapPolicy::FirstDepositionWins;
    limits.geometry.ownership_boundary_grid_units=4;
    limits.geometry.retain_exact_interface_sections=true;
    limits.max_polygon_work=4000000;
    limits.geometry.max_slab_group_visits=4000000;
    SampledConflictOwnershipResolver sampled;
    const auto unsampled=prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,mapping,path.total_s,
        {},limits,{},false,&sampled);
    SectionOwnership owner{.0001,{}};
    std::shared_ptr<const AuthenticatedDepositedPartition> captured;
    owner.resolve_deposited_partition=[&](const SampledSectionConflict& evidence) {
        const auto* admitted=evidence.frozen_cut_boundary_points();
        REQUIRE(admitted);
        REQUIRE(admitted->size()==evidence.segments().size());
        for (std::size_t i=0;i<admitted->size();++i) {
            auto singleton_limits=evidence.geometry_limits();
            singleton_limits.max_boundary_points=std::min<std::size_t>(1024,singleton_limits.max_boundary_points);
            const auto cut=section_deposition({evidence.segments()[i]},evidence.conflict().bottom_m,
                evidence.bed_z_m(),singleton_limits,{},nullptr,evidence.clip());
            CHECK((*admitted)[i]==cut.boundary_points);
        }
        captured=sampled(evidence).deposited_partition;
        return captured;
    };
    const auto resolved=section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner,&clip);
    REQUIRE(captured);
    REQUIRE(resolved.deposited_partition);
    REQUIRE(resolved.occupied_partition);
    CHECK(resolved.occupied_partition->deposited_overlap_policy()==DepositedOverlapPolicy::FirstDepositionWins);
    CHECK(resolved.boundary_error.uncertainty_nm==0);
    CHECK(captured->policy()==DepositedOverlapPolicy::FirstDepositionWins);
    CHECK(captured->overlap_area_m2()>0.);
    REQUIRE(resolved.regions.size()==(hidden_later_material ? 1u : 2u));
    CHECK(resolved.regions.front().material==mapping[0]);
    std::vector<DepositedSlab> slabs;
    for (const auto& region:resolved.regions)
        slabs.push_back({0.,.0002,region.material,region.footprint,0,resolved.occupied_partition});
    std::vector<const DepositedSlab*> members;
    for (const auto& slab:slabs) members.push_back(&slab);
    CHECK(qualify_complete_occupied_partition(members)>0);
    slabs.front().boundary_uncertainty_nm=7;
    CHECK_THROWS_AS(qualify_complete_occupied_partition(members),std::invalid_argument);
    slabs.front().boundary_uncertainty_nm=0;
    const auto first=section_deposition({nominal.segments[0]},.0001,0.,limits.geometry,{},nullptr,&clip);
    REQUIRE(first.regions.size()==1);
    const auto expected=make_exact_section(first.regions.front().footprint.native()).section;
    CHECK(boolean_exact_sections(resolved.regions.front().footprint.exact(),expected,
        ExactSectionOperation::SymmetricDifference).section.empty());
    if (!hidden_later_material)
        CHECK(boolean_exact_sections(resolved.regions[0].footprint.exact(),resolved.regions[1].footprint.exact(),
            ExactSectionOperation::Intersection).section.empty());
    std::reverse(nominal.segments.begin(),nominal.segments.end());
    const auto reversed=section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner,&clip);
    REQUIRE(reversed.regions.size()==resolved.regions.size());
    for (std::size_t i=0;i<resolved.regions.size();++i) {
        CHECK(reversed.regions[i].material==resolved.regions[i].material);
        CHECK(boolean_exact_sections(reversed.regions[i].footprint.exact(),resolved.regions[i].footprint.exact(),
            ExactSectionOperation::SymmetricDifference).section.empty());
    }
    owner.resolve_deposited_partition=[&](const SampledSectionConflict&) { return captured; };
    CHECK_THROWS_AS(section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner,&clip),std::invalid_argument);
    owner.resolve_deposited_partition=[&](const SampledSectionConflict& evidence) { return sampled(evidence).deposited_partition; };
    auto omitted=nominal.segments;
    omitted.erase(std::remove_if(omitted.begin(),omitted.end(),[](const auto& bead) { return bead.section_source_index==0; }),omitted.end());
    CHECK_THROWS_AS(section_deposition(omitted,.0001,0.,limits.geometry,{},&owner,&clip),std::invalid_argument);
    auto all_omitted=nominal.segments;
    all_omitted.erase(std::remove_if(all_omitted.begin(),all_omitted.end(),[](const auto& bead) {
        return bead.section_source_index!=2;
    }),all_omitted.end());
    CHECK_THROWS_AS(section_deposition(all_omitted,.0001,0.,limits.geometry,{},&owner,&clip),std::invalid_argument);
    owner.z_m=.001;
    const auto empty=section_deposition(nominal.segments,.001,0.,limits.geometry,{},&owner,&clip);
    CHECK(empty.regions.empty());
    REQUIRE(empty.deposited_partition);
    CHECK(empty.deposited_partition->policy()==DepositedOverlapPolicy::FirstDepositionWins);
    owner.z_m=.0001;
    const auto move=path.deposition.front().move_index;
    const auto original_position=result.moves[move].position;
    result.moves[move].position.x()+=.01f;
    CHECK_THROWS_AS(section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner,&clip),std::invalid_argument);
    result.moves[move].position=original_position;
    LayerMaterialOverlap arbitrary(.0001,.0001,0.,0);
    CHECK_THROWS_AS(unsampled(arbitrary),std::invalid_argument);
    limits.geometry.retain_exact_interface_sections=false;
    CHECK_THROWS_AS(prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,mapping,path.total_s,
        {},limits,{},false,&sampled),std::invalid_argument);
    limits.geometry.retain_exact_interface_sections=true;
    bool canceled=false;
    prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,mapping,path.total_s,
        {},limits,[&] { if (canceled) throw ResourceLimit("deposited cancellation"); },false,&sampled);
    canceled=true;
    CHECK_THROWS_WITH(section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner,&clip),"deposited cancellation");
    limits.max_total_work=1;
    CHECK_THROWS_AS(prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,mapping,path.total_s,
        {},limits,{},false,&sampled),ResourceLimit);
}

TEST_CASE("Deposited sampling composes slanted individual cuts exactly without changing the native sampler", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    result.moves.resize(4);
    result.moves[0].position={0.f,0.f,.2f};
    result.moves[1].position={1.f,.333333f,.2f};
    result.moves[2].position={.1f,.2f,.2f};
    result.moves[3].position={1.f,-.222222f,.2f};
    result.print_statistics.modes[0].time=4.;
    auto commands=std::make_shared<DepositionSections>(result.deposition_provenance->generation());
    commands->seal();
    result.deposition_sections=commands;
    const auto path=timed_deposition(result);
    auto nominal=resolve_nominal_deposition(result,path,path.total_s,BeadSection::RectangularEnvelope);
    for (std::size_t i=0;i<nominal.segments.size();++i) {
        nominal.segments[i].section_source_index=i;
        nominal.segments[i].material=0;
    }
    ConflictOwnershipLimits limits;
    limits.geometry.ownership_boundary_grid_units=4;
    limits.geometry.retain_exact_interface_sections=true;
    limits.geometry.max_slab_group_visits=4000000;
    limits.max_polygon_work=4000000;
    limits.deposited_overlap_policy=DepositedOverlapPolicy::FirstDepositionWins;
    ExactSection exact_union;
    for (const auto& bead:nominal.segments) {
        const auto cut=section_deposition({bead},.0001,0.,limits.geometry);
        REQUIRE(cut.regions.size()==1);
        const auto exact=make_exact_section(cut.regions.front().footprint.native());
        exact_union=boolean_exact_sections(exact_union,exact.section,ExactSectionOperation::Union).section;
    }
    const auto legacy=section_deposition(nominal.segments,.0001,0.,limits.geometry);
    REQUIRE(legacy.regions.size()==1);
    CHECK_FALSE(legacy.regions.front().footprint.is_exact());
    const auto rounded=make_exact_section(legacy.regions.front().footprint.native()).section;
    CHECK_FALSE(boolean_exact_sections(exact_union,rounded,ExactSectionOperation::SymmetricDifference).section.empty());
    SampledConflictOwnershipResolver sampled;
    prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,{0,0},path.total_s,
        {},limits,{},false,&sampled);
    std::size_t callbacks=0;
    SectionOwnership owner{.0001,{}};
    owner.resolve_deposited_partition=[&](const SampledSectionConflict& evidence) {
        ++callbacks;
        REQUIRE(evidence.original_regions().size()==1);
        CHECK(evidence.original_regions().front().footprint.is_exact());
        return sampled(evidence).deposited_partition;
    };
    const auto exact=section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner);
    REQUIRE(exact.regions.size()==1);
    CHECK(boolean_exact_sections(exact_union,exact.regions.front().footprint.exact(),
        ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK(exact.ownership_work>0);
    CHECK_THAT(exact.deposited_partition->overlap_area_m2(),Catch::Matchers::WithinAbs(0.,0.));
    const auto native_again=section_deposition(nominal.segments,.0001,0.,limits.geometry);
    CHECK(native_again.regions.front().footprint.native()==legacy.regions.front().footprint.native());
    const auto completed=callbacks;
    auto bounded=limits.geometry;
    bounded.max_slab_group_visits=1;
    CHECK_THROWS_AS(section_deposition(nominal.segments,.0001,0.,bounded,{},&owner),ResourceLimit);
    CHECK(callbacks==completed);
    CHECK_THROWS_WITH(section_deposition(nominal.segments,.0001,0.,limits.geometry,
        [] { throw ResourceLimit("exact sampler cancellation"); },&owner),"exact sampler cancellation");

    // Restoring caller inputs after sampling must not replace the immutable
    // geometry bound to evidence, even for a single-material shortcut.
    auto forged=nominal.segments;
    for (auto& bead:forged) bead.width_m*=1.5;
    bool mutated=false;
    std::shared_ptr<const AuthenticatedDepositedPartition> forged_token;
    owner.resolve_deposited_partition=[&](const SampledSectionConflict& evidence) {
        REQUIRE(evidence.original_regions().size()==1);
        REQUIRE(evidence.original_regions().front().footprint.is_exact());
        CHECK_FALSE(boolean_exact_sections(evidence.original_regions().front().footprint.exact(),exact_union,
            ExactSectionOperation::SymmetricDifference).section.empty());
        // Replacing the caller vector cannot rewrite the entry-time snapshot.
        forged=nominal.segments;
        mutated=true;
        forged_token=sampled(evidence).deposited_partition;
        return forged_token;
    };
    CHECK_THROWS_WITH(section_deposition(forged,.0001,0.,limits.geometry,{},&owner),
        "Sampled ownership source geometry or material was altered.");
    CHECK(mutated);
    CHECK_FALSE(forged_token);
    forged=nominal.segments;
    for (auto& bead:forged) bead.width_m*=1.5;
    owner.resolve_deposited_partition=[&](const SampledSectionConflict& evidence) {
        return sampled(evidence).deposited_partition;
    };
    // A snapshot taken only in the evidence constructor would be too late:
    // the very first checkpoint can mutate the caller's composition inputs.
    bool restored=false;
    CHECK_THROWS_WITH(section_deposition(forged,.0001,0.,limits.geometry,[&] {
        if (!restored) { forged=nominal.segments; restored=true; }
    },&owner),"Sampled ownership source geometry or material was altered.");
    CHECK(restored);
    auto mutable_segments=nominal.segments;
    auto mutable_limits=limits.geometry;
    ExPolygon mutable_clip(Polygon{{-1000000,-1000000},{2000000,-1000000},
        {2000000,2000000},{-1000000,2000000}});
    bool changed=false;
    const auto frozen=section_deposition(mutable_segments,.0001,0.,mutable_limits,[&] {
        if (changed) return;
        changed=true;
        for (auto& bead:mutable_segments) bead.width_m*=1.5;
        mutable_limits.retain_exact_interface_sections=false;
        mutable_clip.contour.points.clear();
        owner.z_m=.5;
        owner.resolve_deposited_partition={};
    },&owner,&mutable_clip);
    CHECK(changed);
    REQUIRE(frozen.regions.size()==1);
    CHECK(boolean_exact_sections(frozen.regions.front().footprint.exact(),exact_union,
        ExactSectionOperation::SymmetricDifference).section.empty());

    // Cut evidence cannot inherit a more permissive sampler construction cap.
    auto strict_limits=limits;
    strict_limits.geometry.max_boundary_points=4;
    SampledConflictOwnershipResolver strict_sampled;
    prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,{0,0},path.total_s,
        {},strict_limits,{},false,&strict_sampled);
    owner.z_m=.0001;
    owner.resolve_deposited_partition=[&](const SampledSectionConflict& evidence) {
        return strict_sampled(evidence).deposited_partition;
    };
    CHECK_THROWS_AS(section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner),ResourceLimit);
}

TEST_CASE("Frozen chronological material runs match independent command replay", "[ThermoelasticConflictOwnership]")
{
    const bool returning_material=GENERATE(false,true);
    const bool hidden_foreign=GENERATE(false,true);
    GCodeProcessorResult result;
    populate(result);
    const auto first_start=result.moves[0],first_end=result.moves[1];
    auto second_start=result.moves[2],second_end=result.moves[3];
    if (hidden_foreign) {
        second_start.position=first_start.position;
        second_end.position=first_end.position;
    }
    result.moves.clear();
    for (unsigned material=0;material<2;++material) for (unsigned i=0;i<24;++i) {
        result.moves.push_back(material ? second_start : first_start);
        result.moves.push_back(material ? second_end : first_end);
    }
    if (returning_material) {
        auto start=first_start,end=first_end;
        start.position.x()+=.3f; end.position.x()+=.3f;
        result.moves.push_back(start); result.moves.push_back(end);
    }
    result.print_statistics.modes[0].time=double(result.moves.size());
    auto commands=std::make_shared<DepositionSections>(result.deposition_provenance->generation());
    commands->seal(); result.deposition_sections=commands;
    const auto path=timed_deposition(result);
    REQUIRE(path.state==ToolpathState::Complete);
    auto nominal=resolve_nominal_deposition(result,path,path.total_s,BeadSection::RectangularEnvelope);
    const std::vector<std::size_t> mapping{1,0};
    for (std::size_t i=0;i<nominal.segments.size();++i) {
        nominal.segments[i].section_source_index=i;
        nominal.segments[i].material=mapping[nominal.segments[i].material];
    }
    const auto center=.5*(nominal.segments.front().start_m+nominal.segments.front().end_m);
    const coord_t x=std::llround(center.x()*1e9),y=std::llround(center.y()*1e9);
    const ExPolygon clip(Polygon{{x-500000,y-100000},{x+700000,y-100000},
        {x+700000,y+100000},{x-500000,y+100000}});
    ConflictOwnershipLimits limits;
    limits.deposited_overlap_policy=DepositedOverlapPolicy::FirstDepositionWins;
    limits.geometry.ownership_boundary_grid_units=4;
    limits.geometry.retain_exact_interface_sections=true;
    limits.geometry.max_slab_group_visits=4000000;
    limits.max_polygon_work=4000000;
    std::vector<ExactDepositedCommandSection> independent;
    for (const auto& bead:nominal.segments) {
        const auto cut=section_deposition({bead},.0001,0.,limits.geometry,{},nullptr,&clip);
        REQUIRE(cut.regions.size()==1);
        independent.push_back({bead.material,make_exact_section(cut.regions.front().footprint.native()).section});
    }
    // The retained primitive deliberately does not batch material runs.
    const auto replay=exact_first_deposition_partition(independent);
    SampledConflictOwnershipResolver sampled;
    prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,mapping,path.total_s,
        {},limits,{},false,&sampled);
    SectionOwnership owner{.0001,{}};
    owner.resolve_deposited_partition=[&](const SampledSectionConflict& evidence) { return sampled(evidence).deposited_partition; };
    // Caller order must not replace the authenticated physical chronology.
    std::reverse(nominal.segments.begin(),nominal.segments.end());
    const auto actual=section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner,&clip);
    REQUIRE(actual.regions.size()==replay.regions.size());
    REQUIRE(actual.deposited_partition);
    CHECK(actual.regions.front().material==1);
    for (std::size_t i=0;i<actual.regions.size();++i) {
        CHECK(actual.regions[i].material==replay.regions[i].material);
        CHECK(boolean_exact_sections(actual.regions[i].footprint.exact(),replay.regions[i].envelope,
            ExactSectionOperation::SymmetricDifference).section.empty());
    }
    CHECK_THAT(actual.deposited_partition->overlap_area_m2(),
        Catch::Matchers::WithinRel(replay.cross_material_overlap_area_m2,1e-12));
    CHECK_THROWS_WITH(section_deposition(nominal.segments,.0001,0.,limits.geometry,
        [] { throw std::runtime_error("cancel chronological run reuse"); },&owner,&clip),"cancel chronological run reuse");
    auto bounded=limits;
    bounded.max_polygon_work=1;
    SampledConflictOwnershipResolver limited;
    prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,mapping,path.total_s,
        {},bounded,{},false,&limited);
    SectionOwnership bounded_owner{.0001,{}};
    bounded_owner.resolve_deposited_partition=[&](const SampledSectionConflict& evidence) { return limited(evidence).deposited_partition; };
    CHECK_THROWS_AS(section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&bounded_owner,&clip),ResourceLimit);
    result.moves.back().delta_extruder=0.f;
    CHECK_THROWS_AS(section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner,&clip),std::invalid_argument);
}

TEST_CASE("Complete deposited reference reuse preserves every band and triple-overlap exclusions", "[ThermoelasticConflictOwnership]")
{
    const unsigned shape=GENERATE(0u,1u,2u);
    GCodeProcessorResult result;
    populate(result,false,0,3);
    result.moves.resize(6);
    auto commands=std::make_shared<DepositionSections>(result.deposition_provenance->generation());
    for (std::size_t i=0;i<3;++i) {
        auto& start=result.moves[2*i]; auto& end=result.moves[2*i+1];
        start.position={0.f,0.f,.2f};
        end.position={1.f,.333333f,shape ? .4f : .2f};
        DepositionCommandSection command;
        command.start={start.position.x(),start.position.y(),start.position.z()};
        command.end={end.position.x(),end.position.y(),end.position.z()};
        command.width=.4;
        command.start_height=.2;
        command.end_height=shape==2 ? .3 : .2;
        command.filament=unsigned(i);
        end.deposition_command_section_id=commands->add(command);
        end.deposition_design_width=float(command.width);
        // Match the emitter/processor's child metadata, including legacy
        // float-height rounding for an unanchored parent command.
        DepositionSectionChain chain(*commands->section(end.deposition_command_section_id));
        const auto heights=chain.append(command.start,command.end,command.filament);
        REQUIRE(heights);
        REQUIRE(chain.complete());
        end.deposition_design_heights=*heights;
    }
    // Repeating material0 contributes no additional nominal material union.
    const auto repeated_start=result.moves[0];
    auto repeated_end=result.moves[1];
    const auto repeated_command=*commands->section(repeated_end.deposition_command_section_id);
    repeated_end.deposition_command_section_id=commands->add(repeated_command);
    result.moves.push_back(repeated_start); result.moves.push_back(repeated_end);
    commands->seal(); result.deposition_sections=commands;
    result.print_statistics.modes[0].time=double(result.moves.size());
    const auto path=timed_deposition(result);
    REQUIRE(path.state==ToolpathState::Complete);
    auto nominal=resolve_nominal_deposition(result,path,path.total_s,BeadSection::RectangularEnvelope);
    REQUIRE(nominal.segments.size()==4);
    for (std::size_t i=0;i<nominal.segments.size();++i) nominal.segments[i].section_source_index=i;
    ConflictOwnershipLimits limits;
    limits.deposited_overlap_policy=DepositedOverlapPolicy::FirstDepositionWins;
    limits.geometry.ownership_boundary_grid_units=4;
    limits.geometry.retain_exact_interface_sections=true;
    limits.geometry.max_slab_group_visits=4000000;
    limits.max_polygon_work=4000000;
    SectionGrowthStreamSettings settings;
    CHECK_FALSE(settings.reuse_complete_deposited_reference);
    settings.exact_reference_sources=true;
    settings.geometry.geometry=limits.geometry;
    settings.topology.exact_interior_components=true;
    settings.origin_m={-.001,-.001}; settings.tile_size_m={.004,.004};
    settings.maximum_height_m=.00005;
    const auto run=[&](bool reuse,std::vector<LocalSectionTile>& packets,bool mutate=false,bool partial=false) {
        SampledConflictOwnershipResolver sampled;
        prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,{0,1,2},path.total_s,
            {},limits,{},false,&sampled);
        const auto ownership=[&](double z)->std::optional<SectionOwnership> {
            SectionOwnership owner{z,{}};
            owner.resolve_deposited_partition=[&](const SampledSectionConflict& evidence) { return sampled(evidence).deposited_partition; };
            return owner;
        };
        auto reference=nominal.segments,source=nominal.segments;
        if (partial) source.pop_back();
        auto selected=settings; selected.reuse_complete_deposited_reference=reuse;
        bool changed=false;
        const auto value=visit_reference_activation(reference,source,selected,[](const SectionPiece&) {},ownership,[&] {
            if (mutate && !changed) { reference.clear(); source.clear(); changed=true; }
        },{[](const LocalSectionTile&) {},[&](const LocalSectionTile& tile) { packets.push_back(tile); },[] {}});
        if (mutate) CHECK(changed);
        return value;
    };
    std::vector<LocalSectionTile> baseline,optimized,mutated,partial;
    const auto ordinary=run(false,baseline);
    const auto reused=run(true,optimized);
    const auto frozen=run(true,mutated,true);
    const auto prefix=run(true,partial,false,true);
    REQUIRE_FALSE(baseline.empty());
    const auto compare=[&](const std::vector<LocalSectionTile>& packets) {
        REQUIRE(packets.size()==baseline.size());
        for (std::size_t i=0;i<packets.size();++i) {
            REQUIRE(packets[i].slabs.size()==baseline[i].slabs.size());
            CHECK(packets[i].index==baseline[i].index);
            for (std::size_t j=0;j<packets[i].slabs.size();++j) {
                const auto& a=packets[i].slabs[j]; const auto& b=baseline[i].slabs[j];
                CHECK_THAT(a.bottom_m,Catch::Matchers::WithinAbs(b.bottom_m,0.));
                CHECK_THAT(a.top_m,Catch::Matchers::WithinAbs(b.top_m,0.));
                CHECK(a.material==b.material);
                CHECK(boolean_exact_sections(a.footprint.exact(),b.footprint.exact(),ExactSectionOperation::SymmetricDifference).section.empty());
                const auto bottom=coord_t(std::llround(a.bottom_m*1e9)),top=coord_t(std::llround(a.top_m*1e9));
                const auto am=exact_section_volume_moments(a.footprint.exact(),{-1000000,-1000000,0},{4000000,4000000,1000000},bottom,top);
                const auto bm=exact_section_volume_moments(b.footprint.exact(),{-1000000,-1000000,0},{4000000,4000000,1000000},bottom,top);
                for (std::size_t k=0;k<27;++k) CHECK_THAT(am.values[k],Catch::Matchers::WithinAbs(bm.values[k],1e-25));
                if (a.occupied_partition) {
                    CHECK_THAT(a.occupied_partition->deposited_excluded_nominal_area_m2(),
                        Catch::Matchers::WithinRel(2.*am.values[0]/(a.top_m-a.bottom_m),1e-10));
                }
            }
        }
    };
    compare(optimized); compare(mutated); compare(partial);
    CHECK(ordinary.excluded_nominal_volume_m3>0.);
    CHECK_THAT(reused.excluded_nominal_volume_m3,Catch::Matchers::WithinRel(ordinary.excluded_nominal_volume_m3,1e-10));
    CHECK_THAT(frozen.excluded_nominal_volume_m3,Catch::Matchers::WithinRel(ordinary.excluded_nominal_volume_m3,1e-10));
    CHECK_THAT(prefix.excluded_nominal_volume_m3,Catch::Matchers::WithinRel(ordinary.excluded_nominal_volume_m3,1e-10));
    CHECK(reused.active.components()==ordinary.active.components());
    CHECK(reused.reference.source()==reused.active.source());
    CHECK(reused.reference.component_mapping()==reused.active.component_mapping());
    CHECK(reused.reference.components()==ordinary.reference.components());
    CHECK(reused.reference.free_components()==ordinary.reference.free_components());
    CHECK(frozen.reference.source()==frozen.active.source());
    CHECK(ordinary.reference.source()!=ordinary.active.source());
    CHECK(prefix.reference.source()!=prefix.active.source());
    if (shape==0) {
        SampledConflictOwnershipResolver sampled;
        prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,{0,1,2},path.total_s,
            {},limits,{},false,&sampled);
        const auto ownership=[&](double z)->std::optional<SectionOwnership> {
            SectionOwnership owner{z,{}};
            owner.resolve_deposited_partition=[&](const SampledSectionConflict& evidence) { return sampled(evidence).deposited_partition; };
            return owner;
        };
        ThermalSectionStreamSettings thermal;
        thermal.origin_m={-.001,-.001,0.}; thermal.spacing_m={.001,.001,.0002};
        thermal.tile_size_m=settings.tile_size_m; thermal.maximum_height_m=settings.maximum_height_m;
        thermal.geometry=settings.geometry; thermal.topology=settings.topology;
        const std::vector<SectionThermalMaterial> properties{
            {Eigen::Matrix3d::Identity()*.2,2e6,400.},
            {Eigen::Matrix3d::Identity()*.1,1e6,400.},
            {Eigen::Matrix3d::Identity()*.3,3e6,400.}};
        const auto baseline_thermal=assemble_reference_thermal_volume(nominal.segments,nominal.segments,properties,thermal,settings,ownership);
        auto optimized_settings=settings; optimized_settings.reuse_complete_deposited_reference=true;
        const auto optimized_thermal=assemble_reference_thermal_volume(nominal.segments,nominal.segments,properties,thermal,optimized_settings,ownership);
        CHECK(optimized_thermal.assembly.exposure_complete);
        CHECK(optimized_thermal.assembly.groups.size()==baseline_thermal.assembly.groups.size());
        CHECK_THAT(optimized_thermal.assembly.occupied_volume_m3,
            Catch::Matchers::WithinRel(baseline_thermal.assembly.occupied_volume_m3,1e-10));
        CHECK_THAT(optimized_thermal.excluded_nominal_volume_m3,
            Catch::Matchers::WithinRel(baseline_thermal.excluded_nominal_volume_m3,1e-10));
    }
}

TEST_CASE("Streamed deposited ownership bounds live geometry while authenticating every repeated command", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    const auto travel=result.moves[0],extrusion=result.moves[1];
    constexpr std::size_t repeats=1100;
    result.moves.clear();
    result.moves.reserve(2*repeats);
    for (std::size_t i=0;i<repeats;++i) {
        result.moves.push_back(travel);
        result.moves.push_back(extrusion);
    }
    result.print_statistics.modes[0].time=double(result.moves.size());
    auto commands=std::make_shared<DepositionSections>(result.deposition_provenance->generation());
    commands->seal();
    result.deposition_sections=commands;
    const auto path=timed_deposition(result);
    REQUIRE(path.state==ToolpathState::Complete);
    auto nominal=resolve_nominal_deposition(result,path,path.total_s,BeadSection::RectangularEnvelope);
    REQUIRE(nominal.segments.size()==repeats);
    for (std::size_t i=0;i<nominal.segments.size();++i) nominal.segments[i].section_source_index=i;
    const auto center=.5*(nominal.segments.front().start_m+nominal.segments.front().end_m);
    const coord_t x=std::llround(center.x()*1e9),y=std::llround(center.y()*1e9);
    const ExPolygon clip(Polygon{{x-400000,y-100000},{x+400000,y-100000},
        {x+400000,y+100000},{x-400000,y+100000}});
    ConflictOwnershipLimits limits;
    limits.deposited_overlap_policy=DepositedOverlapPolicy::FirstDepositionWins;
    limits.geometry.ownership_boundary_grid_units=4;
    limits.geometry.retain_exact_interface_sections=true;
    limits.geometry.max_slab_group_visits=4000000;
    limits.max_polygon_work=4000000;
    const auto single=section_deposition({nominal.segments.front()},.0001,0.,limits.geometry,{},nullptr,&clip);
    REQUIRE(single.regions.size()==1);
    const auto expected=make_exact_section(single.regions.front().footprint.native()).section;
    CHECK(expected.vertex_count()*repeats>4096);
    SampledConflictOwnershipResolver sampled;
    prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,{0,1},path.total_s,
        {},limits,{},false,&sampled);
    std::size_t callbacks=0;
    SectionOwnership owner{.0001,{}};
    owner.resolve_deposited_partition=[&](const SampledSectionConflict& evidence) {
        ++callbacks;
        const auto* admitted=evidence.frozen_cut_boundary_points();
        REQUIRE(admitted);
        REQUIRE(admitted->size()==repeats);
        CHECK(std::all_of(admitted->begin(),admitted->end(),[&](std::size_t points) {
            return points==single.boundary_points;
        }));
        return sampled(evidence).deposited_partition;
    };
    const auto resolved=section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner,&clip);
    REQUIRE(resolved.regions.size()==1);
    CHECK(resolved.regions.front().material==0);
    CHECK(boolean_exact_sections(resolved.regions.front().footprint.exact(),expected,
        ExactSectionOperation::SymmetricDifference).section.empty());
    REQUIRE(resolved.deposited_partition);
    CHECK_THAT(resolved.deposited_partition->overlap_area_m2(),Catch::Matchers::WithinAbs(0.,0.));
    CHECK(resolved.boundary_points==expected.vertex_count());
    CHECK(resolved.ownership_work<=limits.geometry.max_slab_group_visits);
    auto composition_limits=limits.geometry;
    composition_limits.exact_per_command_composition=true;
    const auto composed=section_deposition(nominal.segments,.0001,0.,composition_limits,{},nullptr,&clip);
    REQUIRE(resolved.ownership_work>=composed.ownership_work);
    // Independent replay necessarily reconstructs every identical cut. The
    // frozen sole-owner proof costs less than those constructions alone.
    const auto construction=make_exact_section(single.regions.front().footprint.native()).work+single.boundary_points+1;
    CHECK(resolved.ownership_work-composed.ownership_work<repeats*construction);
    // A completely covered later command still requires live authentication.
    result.moves.back().delta_extruder=0.f;
    CHECK_THROWS_AS(section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner,&clip),std::invalid_argument);
    const auto prior_callbacks=callbacks;
    const std::vector<DepositedSegment> excessive(4097,nominal.segments.front());
    CHECK_THROWS_WITH(section_deposition(excessive,.0001,0.,limits.geometry,{},&owner,&clip),
        "Deposited sampler input snapshot command budget exceeded.");
    CHECK(callbacks==prior_callbacks);
}

TEST_CASE("Deposited order preserves retained geometry under the existing bounded stationary omission", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    result.filament_diameters={1.75f,1.75f};
    for (auto& move:result.moves)
        if (move.type==EMoveType::Extrude) move.mm3_per_mm=.08f;
    auto commands=std::make_shared<DepositionSections>(result.deposition_provenance->generation());
    commands->seal();
    result.deposition_sections=commands;
    ConflictOwnershipLimits limits;
    limits.deposited_overlap_policy=DepositedOverlapPolicy::FirstDepositionWins;
    limits.geometry.ownership_boundary_grid_units=4;
    limits.geometry.retain_exact_interface_sections=true;
    limits.max_polygon_work=4000000;
    limits.geometry.max_slab_group_visits=4000000;
    StationaryExtrusionSummary summary;
    const auto resolve=[&](const std::optional<StationaryExtrusionAllowance>& allowance) {
        const auto path=timed_deposition(result);
        auto nominal=resolve_nominal_deposition(result,path,path.total_s,BeadSection::RectangularEnvelope,allowance);
        summary=nominal.stationary_extrusion;
        for (std::size_t i=0;i<nominal.segments.size();++i) nominal.segments[i].section_source_index=i;
        const auto center=.5*(nominal.segments[0].start_m+nominal.segments[0].end_m);
        const coord_t x=std::llround(.5*(nominal.segments[0].start_m.x()+nominal.segments[1].start_m.x())*1e9);
        const coord_t y=std::llround(center.y()*1e9);
        const ExPolygon clip(Polygon{{x-400000,y-100000},{x+400000,y-100000},
            {x+400000,y+100000},{x-400000,y+100000}});
        SampledConflictOwnershipResolver sampled;
        prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,{0,1},path.total_s,
            allowance,limits,{},false,&sampled);
        SectionOwnership owner{.0001,{}};
        owner.resolve_deposited_partition=[&](const SampledSectionConflict& evidence) { return sampled(evidence).deposited_partition; };
        return section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner,&clip);
    };
    const auto baseline=resolve({});
    CHECK(summary.count==0);
    auto stationary=result.moves[1];
    stationary.delta_extruder=1e-9f;
    result.moves.insert(result.moves.begin()+2,stationary);
    result.print_statistics.modes[0].time=double(result.moves.size());
    const StationaryExtrusionAllowance allowance{1e-12,1e-6,1e-3};
    const auto retained=resolve(allowance);
    CHECK(summary.count==1);
    CHECK(summary.commanded_volume_m3>0.);
    CHECK(summary.commanded_volume_m3<=allowance.maximum_volume_m3);
    REQUIRE(retained.regions.size()==baseline.regions.size());
    for (std::size_t i=0;i<baseline.regions.size();++i) {
        CHECK(retained.regions[i].material==baseline.regions[i].material);
        CHECK(boolean_exact_sections(retained.regions[i].footprint.exact(),baseline.regions[i].footprint.exact(),
            ExactSectionOperation::SymmetricDifference).section.empty());
    }
    CHECK_THAT(retained.deposited_partition->overlap_area_m2(),
        Catch::Matchers::WithinAbs(baseline.deposited_partition->overlap_area_m2(),0.));
    SampledConflictOwnershipResolver sampled;
    CHECK_THROWS_AS(prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,{0,1},
        double(result.moves.size()),{},limits,{},false,&sampled),std::invalid_argument);
    CHECK_THROWS_AS(prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,{0,1},
        2.5,allowance,limits,{},false,&sampled),std::invalid_argument);
    result.moves[2].position.z()+=.01f;
    CHECK_THROWS_AS(prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,{0,1},
        double(result.moves.size()),allowance,limits,{},false,&sampled),std::invalid_argument);
}

TEST_CASE("Authenticated sampled ownership avoids rebuilding unused material unions", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    const auto path = timed_deposition(result);
    auto nominal = resolve_nominal_deposition(result, path, path.total_s, BeadSection::RectangularEnvelope);
    for (std::size_t i = 0; i < nominal.segments.size(); ++i) nominal.segments[i].section_source_index = i;
    const auto resolve = [&](bool use_sampled, std::size_t polygon_budget) {
        ConflictOwnershipLimits limits;
        limits.max_polygon_work = polygon_budget;
        SampledConflictOwnershipResolver sampled;
        const auto native = prepare_native_model_ownership(result, BeadSection::RectangularEnvelope, 0., {0, 1},
            path.total_s, {}, limits, {}, false, &sampled);
        SectionOwnership owner{.0001, {}};
        if (use_sampled)
            owner.resolve_sampled_conflict = [&](const SampledSectionConflict& evidence) { return sampled(evidence).regions; };
        else
            owner.resolve_conflict = [&](const LayerMaterialOverlap& conflict) { return native(conflict).regions; };
        return section_deposition(nominal.segments, owner.z_m, 0., limits.geometry, {}, &owner);
    };
    // Discover the required work rather than freezing incidental vertex counts.
    const std::size_t ceiling = 200000;
    const auto baseline = resolve(false, ceiling);
    const auto minimum_budget = [&](bool use_sampled) {
        std::size_t low = 1, high = ceiling;
        while (low < high) {
            const auto middle = low + (high - low) / 2;
            try {
                resolve(use_sampled, middle);
                high = middle;
            } catch (const ResourceLimit&) {
                low = middle + 1;
            }
        }
        return low;
    };
    const auto sampled_budget = minimum_budget(true);
    CHECK(sampled_budget < minimum_budget(false));
    REQUIRE(sampled_budget > 1);
    CHECK_THROWS_AS(resolve(true, sampled_budget - 1), ResourceLimit);
    const auto sampled_section = resolve(true, sampled_budget);
    REQUIRE(sampled_section.regions.size() == baseline.regions.size());
    for (std::size_t i = 0; i < baseline.regions.size(); ++i) {
        CHECK(sampled_section.regions[i].material == baseline.regions[i].material);
        CHECK(sampled_section.regions[i].footprint.native() == baseline.regions[i].footprint.native());
    }
}

TEST_CASE("Partition failure observers borrow bounded operands without replacing rejection", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    ConflictOwnershipLimits limits;
    limits.geometry.ownership_boundary_grid_units=4;
    limits.geometry.retain_exact_interface_sections=true;
    limits.geometry.max_slab_group_visits=4000000;
    limits.max_polygon_work=4000000;
    limits.interface_variant=NativeInterfaceVariant::FavorFirst;
    limits.allow_single_sided_native_boundary=true;
    std::size_t observations=0;
    limits.observe_partition_failure=[&](const NativeOwnershipPartitionFailure& failure) {
        ++observations;
        CHECK(failure.variant==limits.interface_variant);
        CHECK_THAT(failure.sampled.conflict().bottom_m,Catch::Matchers::WithinAbs(.0001,0.));
        CHECK(failure.max_work<=4000000);
        CHECK(failure.allow_single_sided_native_boundary);
        CHECK_FALSE(failure.message.empty());
        for (const auto* pair:{&failure.native_masks,&failure.envelopes})
            for (const auto& polygons:*pair) {
                std::size_t points=0;
                for (const auto& polygon:polygons) {
                    points+=polygon.contour.points.size();
                    for (const auto& hole:polygon.holes) points+=hole.points.size();
                }
                CHECK(points<=2048);
            }
        throw std::logic_error("An observer cannot replace the ownership failure.");
    };
    const auto resolve=[&] {
        const auto path=timed_deposition(result);
        auto nominal=resolve_nominal_deposition(result,path,path.total_s,BeadSection::RectangularEnvelope);
        for (std::size_t i=0;i<nominal.segments.size();++i) nominal.segments[i].section_source_index=i;
        SampledConflictOwnershipResolver sampled;
        prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,{0,1},path.total_s,
            {},limits,{},false,&sampled);
        SectionOwnership owner{.0001,{}};
        owner.resolve_partition=[&](const SampledSectionConflict& evidence) { return sampled(evidence).partition; };
        return section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner);
    };
    REQUIRE(resolve().regions.size()==2);
    CHECK(observations==0);
    // Keep both authenticated bead envelopes together, but move their conflict
    // outside the captured native masks. No owner can be inferred there.
    for (std::size_t i=0;i<4;++i) result.moves[i].position.x()+=10.f;
    CHECK_THROWS_AS(resolve(),std::invalid_argument);
    CHECK(observations==1);
}

TEST_CASE("Authenticated occupied partitions retain dual boundary ambiguity area", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    const coord_t shift=1,half_height=200000;
    populate(result,false,shift);
    const auto path=timed_deposition(result);
    auto nominal=resolve_nominal_deposition(result,path,path.total_s,BeadSection::RectangularEnvelope);
    for (std::size_t i=0;i<nominal.segments.size();++i) nominal.segments[i].section_source_index=i;
    const auto center=.5*(nominal.segments[0].start_m+nominal.segments[0].end_m);
    const coord_t x=std::llround(center.x()*1e9),y=std::llround(center.y()*1e9);
    const ExPolygon clip(Polygon{{x-1000000,y-half_height},{x+1000000,y-half_height},
        {x+1000000,y+half_height},{x-1000000,y+half_height}});
    for (const auto variant:{NativeInterfaceVariant::FavorFirst,NativeInterfaceVariant::FavorSecond}) {
        ConflictOwnershipLimits limits;
        limits.geometry.ownership_boundary_grid_units=4;
        limits.geometry.retain_exact_interface_sections=true;
        limits.geometry.max_slab_group_visits=4000000;
        limits.max_polygon_work=4000000;
        limits.interface_variant=variant;
        limits.allow_dual_native_boundary=true;
        SampledConflictOwnershipResolver sampled;
        prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,{0,1},path.total_s,
            {},limits,{},false,&sampled);
        SectionOwnership owner{.0001,{}};
        owner.resolve_partition=[&](const SampledSectionConflict& evidence) { return sampled(evidence).partition; };
        const auto section=section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner,&clip);
        REQUIRE(section.occupied_partition);
        // Two one-nanometre outward shifts leave a2nm reference gap across
        // the complete straight-sided0.4mm clipped deposition interval.
        const double expected_area=double(2*shift)*double(2*half_height)*1e-18;
        CHECK_THAT(section.occupied_partition->dual_native_boundary_ambiguity_area_m2(),
            Catch::Matchers::WithinRel(expected_area,1e-12));
    }
}

TEST_CASE("Native partition budget failures identify local and prepared work pools", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    const auto path=timed_deposition(result);
    auto nominal=resolve_nominal_deposition(result,path,path.total_s,BeadSection::RectangularEnvelope);
    for (std::size_t i=0;i<nominal.segments.size();++i) nominal.segments[i].section_source_index=i;
    bool contextual_failure=false;
    for (std::size_t budget=128;budget<=65536 && !contextual_failure;budget*=2) {
        ConflictOwnershipLimits limits;
        limits.geometry.ownership_boundary_grid_units=4;
        limits.geometry.retain_exact_interface_sections=true;
        limits.geometry.max_slab_group_visits=4000000;
        limits.max_polygon_work=budget;
        limits.interface_variant=NativeInterfaceVariant::FavorFirst;
        try {
            SampledConflictOwnershipResolver sampled;
            prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,{0,1},path.total_s,
                {},limits,{},false,&sampled);
            SectionOwnership owner{.0001,{}};
            owner.resolve_partition=[&](const SampledSectionConflict& evidence) { return sampled(evidence).partition; };
            section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner);
        } catch (const ResourceLimit& error) {
            const std::string message=error.what();
            if (message.find("Native ownership partition budget context:")==std::string::npos) continue;
            contextual_failure=true;
            CHECK(message.find("limiting_pool=local-polygon")!=std::string::npos);
            CHECK(message.find("prepared_total_remaining=")!=std::string::npos);
            CHECK(message.find("prepared_total_limit="+std::to_string(limits.max_total_work))!=std::string::npos);
            CHECK(message.find("local_polygon_used=")!=std::string::npos);
            CHECK(message.find("local_polygon_limit="+std::to_string(budget))!=std::string::npos);
            CHECK(message.find("factory_available=")!=std::string::npos);
            const auto scalar=[&](const std::string& key) {
                const auto offset=message.find(key+"=");
                REQUIRE(offset!=std::string::npos);
                return std::stoull(message.substr(offset+key.size()+1));
            };
            const auto source_work=scalar("prepared_source_work"),partition_work=scalar("prepared_partition_work");
            CHECK(source_work>0);
            CHECK(partition_work==0);
            CHECK(scalar("prepared_native_winding_work")==0);
            CHECK(source_work+partition_work+scalar("prepared_other_work")==scalar("prepared_charged_work"));
            CHECK(message.find("; ")!=std::string::npos);
        }
    }
    CHECK(contextual_failure);
}

TEST_CASE("Prepared ownership work diagnostics retain accepted stage totals and preserve cancellation", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    const auto path=timed_deposition(result);
    auto nominal=resolve_nominal_deposition(result,path,path.total_s,BeadSection::RectangularEnvelope);
    for (std::size_t i=0;i<nominal.segments.size();++i) nominal.segments[i].section_source_index=i;
    ConflictOwnershipLimits limits;
    limits.geometry.ownership_boundary_grid_units=4;
    limits.geometry.retain_exact_interface_sections=true;
    limits.geometry.max_slab_group_visits=4000000;
    limits.max_polygon_work=4000000;
    limits.max_total_work=1000000;
    limits.interface_variant=NativeInterfaceVariant::FavorFirst;
    SampledConflictOwnershipResolver sampled;
    prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,{0,1},path.total_s,{},limits,{},false,&sampled);
    SectionOwnership owner{.0001,{}};
    owner.resolve_partition=[&](const SampledSectionConflict& evidence) { return sampled(evidence).partition; };
    std::size_t accepted=0;
    std::string failure;
    for (std::size_t attempt=0;attempt<256 && failure.empty();++attempt) {
        try {
            section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner);
            ++accepted;
        } catch (const ResourceLimit& error) { failure=error.what(); }
    }
    REQUIRE(accepted>0);
    REQUIRE_FALSE(failure.empty());
    const auto scalar=[&](const std::string& key) {
        const auto offset=failure.find(key+"=");
        REQUIRE(offset!=std::string::npos);
        return std::stoull(failure.substr(offset+key.size()+1));
    };
    const auto source_work=scalar("prepared_source_work"),partition_work=scalar("prepared_partition_work");
    CHECK(source_work>0);
    CHECK(partition_work>0);
    CHECK(scalar("prepared_native_winding_work")>0);
    CHECK(scalar("prepared_native_winding_work")<=partition_work);
    CHECK(source_work+partition_work+scalar("prepared_other_work")==scalar("prepared_charged_work"));
    CHECK(scalar("prepared_charged_work")<=limits.max_total_work);
    bool armed=false;
    std::size_t callbacks=0;
    limits.observe_masks=[&](double,bool,const auto&) { armed=true; };
    prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,{0,1},path.total_s,{},limits,
        [&] { if (armed && ++callbacks==4) throw ResourceLimit("external cancellation marker"); },false,&sampled);
    try {
        section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner);
        FAIL("The cancellation callback must reject ownership.");
    } catch (const ResourceLimit& error) {
        CHECK(std::string(error.what())=="external cancellation marker");
    }
}

TEST_CASE("Opt in exact ownership retains canonical sections without qualifying native consumers", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    const auto path=timed_deposition(result);
    auto nominal=resolve_nominal_deposition(result,path,path.total_s,BeadSection::RectangularEnvelope);
    for (std::size_t i=0;i<nominal.segments.size();++i) nominal.segments[i].section_source_index=i;
    ConflictOwnershipLimits limits;
    limits.geometry.ownership_boundary_grid_units=4;
    limits.geometry.retain_exact_interface_sections=true;
    limits.geometry.max_slab_group_visits=4000000;
    limits.max_polygon_work=4000000;
    limits.interface_variant=NativeInterfaceVariant::FavorFirst;
    SampledConflictOwnershipResolver sampled;
    prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,{0,1},path.total_s,
        {},limits,{},false,&sampled);
    std::shared_ptr<const AuthenticatedSectionPartition> captured;
    const auto ownership = [&](double z) {
        SectionOwnership owner{z,{}};
        owner.resolve_partition=[&](const SampledSectionConflict& evidence) {
            captured=sampled(evidence).partition;
            return captured;
        };
        return owner;
    };
    auto owner=ownership(.0001);
    const auto section=section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner);
    REQUIRE(captured);
    REQUIRE(section.regions.size()==2);
    CHECK(section.boundary_error.uncertainty_nm==7);
    CHECK_THAT(section.boundary_error.added_area_bound_m2,Catch::Matchers::WithinAbs(0.,0.));
    CHECK_THAT(section.boundary_error.removed_area_bound_m2,Catch::Matchers::WithinAbs(0.,0.));
    std::size_t points=0;
    for (const auto& region : section.regions) {
        REQUIRE(region.footprint.is_exact());
        CHECK_THROWS_AS(region.footprint.native(),std::invalid_argument);
        points+=region.footprint.vertex_count();
    }
    CHECK(section.boundary_points==points);
    REQUIRE(section.occupied_partition);
    std::vector<DepositedSlab> partition;
    for (const auto& region : section.regions)
        partition.push_back({0.,.0002,region.material,region.footprint,7,section.occupied_partition});
    const std::vector<const DepositedSlab*> members{&partition[0],&partition[1]};
    const auto qualification_work=qualify_complete_occupied_partition(members);
    CHECK(qualification_work>0);
    CHECK_THROWS_AS(qualify_complete_occupied_partition({members.front()}),std::invalid_argument);
    CHECK_THROWS_AS(qualify_complete_occupied_partition(members,qualification_work-1),ResourceLimit);
    CHECK_THROWS_AS(qualify_complete_occupied_partition(members,4000000,
        [] { throw std::runtime_error("cancel occupied partition"); }),std::runtime_error);
    const auto original_member=partition[1];
    partition[1].occupied_partition.reset();
    CHECK_THROWS_AS(qualify_complete_occupied_partition(members),std::invalid_argument);
    partition[1]=original_member;
    partition[1].material=partition[0].material;
    CHECK_THROWS_AS(qualify_complete_occupied_partition(members),std::invalid_argument);
    partition[1]=original_member;
    partition[1].footprint=SectionFootprint(ExactSection{});
    CHECK_THROWS_AS(qualify_complete_occupied_partition(members),std::invalid_argument);
    partition[1]=original_member;
    partition[1].boundary_uncertainty_nm=0;
    CHECK_THROWS_AS(qualify_complete_occupied_partition(members),std::invalid_argument);
    partition[1]=original_member;
    partition[1].bottom_m=.0003;
    partition[1].top_m=.0004;
    CHECK_THROWS_AS(qualify_complete_occupied_partition(members),std::invalid_argument);
    partition[1]=original_member;
    for (auto& member : partition) { member.bottom_m=.0003; member.top_m=.0004; }
    CHECK_THROWS_AS(qualify_complete_occupied_partition(members),std::invalid_argument);
    partition[0].bottom_m=0.; partition[0].top_m=.0002;
    partition[1]=original_member;
    auto second_owner=ownership(.0001);
    const auto second_section=section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&second_owner);
    partition[1].occupied_partition=second_section.occupied_partition;
    CHECK_THROWS_AS(qualify_complete_occupied_partition(members),std::invalid_argument);
    partition[1]=original_member;
    CHECK(boolean_exact_sections(section.regions[0].footprint.exact(),section.regions[1].footprint.exact(),
        ExactSectionOperation::Intersection).section.empty());
    const auto occupied=boolean_exact_sections(section.regions[0].footprint.exact(),section.regions[1].footprint.exact(),
        ExactSectionOperation::SymmetricDifference);
    ExPolygons original;
    for (const std::size_t material : {std::size_t(0),std::size_t(1)}) {
        std::vector<DepositedSegment> beads;
        for (const auto& bead : nominal.segments) if (bead.material==material) beads.push_back(bead);
        const auto unpartitioned=section_deposition(beads,.0001,0.,limits.geometry);
        for (const auto& region : unpartitioned.regions)
            original.insert(original.end(),region.footprint.native().begin(),region.footprint.native().end());
    }
    const auto source=make_exact_section(original);
    CHECK(boolean_exact_sections(occupied.section,source.section,ExactSectionOperation::SymmetricDifference).section.empty());
    const auto box=source.section.bounds();
    REQUIRE(box.has_value());
    const std::array<coord_t,3> start{box->min.x(),box->min.y(),0};
    const std::array<coord_t,3> size{box->max.x()-box->min.x(),box->max.y()-box->min.y(),1};
    const auto area=exact_section_horizontal_surface_moments(source.section,start,size,0);
    CHECK_THAT(section.area_m2,Catch::Matchers::WithinRel(area.values[0],1e-12));
    const auto saved=captured;
    owner.resolve_partition=[&](const SampledSectionConflict&) { return saved; };
    CHECK_THROWS_AS(section_deposition(nominal.segments,.0001,0.,limits.geometry,{},&owner),std::invalid_argument);
    owner=ownership(.0001);
    auto mismatched=limits.geometry;
    mismatched.retain_exact_interface_sections=false;
    CHECK_THROWS_AS(section_deposition(nominal.segments,.0001,0.,mismatched,{},&owner),std::invalid_argument);

    LocalSectionLimits local_limits;
    local_limits.geometry=limits.geometry;
    const auto local=local_deposition_sections(nominal.segments,Eigen::Vector2d{-1.,-1.},Eigen::Vector2d{2.,2.},
        0.,.00005,local_limits,ownership);
    std::size_t exact_slabs=0;
    for (const auto& tile : local.tiles)
        for (const auto& slab : tile.slabs) {
            REQUIRE(slab.footprint.is_exact());
            REQUIRE(slab.occupied_partition);
            ++exact_slabs;
            CHECK(slab.boundary_uncertainty_nm==7);
            CHECK(slab.top_m-slab.bottom_m<=.00005+1e-12);
        }
    CHECK(exact_slabs>2);
    CHECK(local.maximum_boundary_uncertainty_nm==7);
    CHECK_THAT(local.volume_m3,Catch::Matchers::WithinRel(section.area_m2*.0002,1e-10));
    CHECK_THAT(local.added_representation_volume_bound_m3,Catch::Matchers::WithinAbs(0.,0.));
    CHECK_THAT(local.removed_representation_volume_bound_m3,Catch::Matchers::WithinAbs(0.,0.));
    auto cleared_uncertainty=local.tiles;
    for (auto& tile : cleared_uncertainty)
        for (auto& slab : tile.slabs) {
            REQUIRE(slab.occupied_partition);
            slab.boundary_uncertainty_nm=0;
        }
    // Clearing every member's sensitivity must not bypass the retained
    // authenticated partition metadata by taking the ordinary zero path.
    CHECK_THROWS_AS(qualify_local_section_topology(cleared_uncertainty,0.),std::invalid_argument);
    const auto qualified=qualify_local_section_topology(local.tiles,0.);
    CHECK(qualified.components()>0);
    for (const auto& piece : qualified.pieces()) {
        CHECK(piece.slab.footprint.is_exact());
        CHECK(piece.slab.boundary_uncertainty_nm==7);
        CHECK_FALSE(piece.slab.occupied_partition);
    }
    SectionGrowthStreamSettings activation;
    activation.origin_m={-1.,-1.};
    activation.tile_size_m={2.,2.};
    activation.maximum_height_m=.00005;
    activation.geometry=local_limits;
    std::size_t activated=0;
    const auto complete=visit_reference_activation(nominal.segments,nominal.segments,activation,
        [&](const SectionPiece& piece) {
            ++activated;
            CHECK(piece.slab.footprint.is_exact());
            CHECK(piece.slab.boundary_uncertainty_nm==7);
        },ownership);
    CHECK(activated>0);
    CHECK(complete.active.components()==complete.reference.components());

    ThermalSectionStreamSettings thermal;
    thermal.origin_m={-1.,-1.,0.};
    thermal.spacing_m={.001,.001,.0002};
    thermal.tile_size_m=activation.tile_size_m;
    thermal.maximum_height_m=activation.maximum_height_m;
    thermal.geometry=local_limits;
    thermal.include_exposure=true;
    const std::vector<SectionThermalMaterial> materials{
        {Eigen::Matrix3d::Identity()*.2,2e6,400.},
        {Eigen::Matrix3d::Identity()*.1,1e6,400.}};
    const auto retained=assemble_section_thermal_volume(qualified,materials,thermal.origin_m,thermal.spacing_m);
    REQUIRE(retained.exposure_complete);
    const auto areas=[](const SectionThermalAssembly& assembly) {
        std::array<double,4> totals{};
        for (const auto& face : assembly.horizontal_faces)
            totals[std::size_t(face.kind)]+=face.area_moments[0];
        for (const auto& face : assembly.vertical_faces) totals[3]+=face.area_moments[0];
        return totals;
    };
    const auto expected=areas(retained);
    for (const unsigned bins : {0u,8u}) {
        thermal.vertical_normal_bins=bins;
        const auto streamed=assemble_local_deposition_thermal_volume(nominal.segments,materials,thermal,ownership);
        REQUIRE(streamed.assembly.exposure_complete);
        CHECK_THAT(streamed.assembly.occupied_volume_m3,
            Catch::Matchers::WithinRel(retained.occupied_volume_m3,1e-10));
        const auto actual=areas(streamed.assembly);
        for (std::size_t i=0;i<actual.size();++i)
            CHECK_THAT(actual[i],Catch::Matchers::WithinAbs(expected[i],1e-16));
    }
    thermal.assembly.max_points=1024;
    CHECK_THROWS_AS(assemble_local_deposition_thermal_volume(nominal.segments,materials,thermal,ownership),ResourceLimit);
}

TEST_CASE("Streamed thermal mechanics consumes matching native material ownership", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    const auto path = timed_deposition(result);
    REQUIRE(path.state == ToolpathState::Complete);
    SectionGCodeScenario scenario(BeadSection::RectangularEnvelope);
    scenario.acknowledge_modeled_timing = true;
    scenario.use_native_ownership = true;
    scenario.thermal_materials = {{Eigen::Matrix3d::Identity()*.2,2e6,400.},
                                  {Eigen::Matrix3d::Identity()*.1,1e6,400.}};
    scenario.expansion = {{{{250.,450.,Eigen::Matrix3d::Identity()*1e-5}},Eigen::Matrix3d::Identity(),400.},
                          {{{250.,450.,Eigen::Matrix3d::Identity()*2e-5}},Eigen::Matrix3d::Identity(),400.}};
    scenario.mechanical_materials = {{Eigen::Matrix<double,6,6>::Identity()*1000.},
                                     {Eigen::Matrix<double,6,6>::Identity()*500.}};
    scenario.filament_to_material = {0,1};
    scenario.spacing_m = {.001,.001,.0002};
    scenario.stream_tile_size_m = Eigen::Vector2d{.002,.002};
    scenario.activation_times_s = {path.total_s};
    scenario.cooldown_s = .1;
    scenario.deposition_temperatures_k = [](double) { return std::vector<double>{400.,400.}; };
    // Explicitly insulated numeric fixture, not an airflow/adhesion scenario.
    scenario.exposure = [](double,const SectionThermalAssembly&) { return std::vector<ThermalHistory::Element>{}; };
    SectionGCodeLimits limits;
    // Explicit nanometre-scale interface sensitivity, not a physical priority
    // or permission to silently discard union loss in the strict default.
    limits.ownership.interface_variant = NativeInterfaceVariant::FavorFirst;
    limits.ownership.geometry.ownership_boundary_grid_units = 4;
    limits.streaming_geometry.geometry.ownership_boundary_grid_units = 4;
    // The adapter prepares matching ownership from the sampler settings;
    // callers must not have to duplicate its plane resolution manually.
    limits.streaming_geometry.geometry.half_native_section_planes = true;
    limits.streaming_geometry.geometry.max_slab_group_visits = 4000000;
    limits.ownership.max_polygon_work = 4000000;
    const auto solved = simulate_section_gcode_history(result,scenario,limits);
    REQUIRE(solved.thermal.state == ThermalHistory::State::Converged);
    REQUIRE(solved.mechanical.has_value());
    CHECK(solved.occupied_volume_m3 > 0.);
    CHECK(solved.activated_snapshots == 1);
    scenario.activation_times_s = {path.deposition[0].end_s,path.deposition[1].end_s,path.total_s};
    // Adding the overlapping foreign bead changes native ownership of volume
    // already activated by the first bead. It is not monotone deposition:
    // accepting it would silently discard/reassign existing material and heat.
    CHECK_THROWS_AS(simulate_section_gcode_history(result,scenario,limits),std::invalid_argument);
    scenario.use_fixed_interface_reference=true;
    const auto fixed_history=simulate_section_gcode_history(result,scenario,limits);
    REQUIRE(fixed_history.thermal.state==ThermalHistory::State::Converged);
    REQUIRE(fixed_history.mechanical.has_value());
    CHECK(fixed_history.uses_fixed_interface_reference);
    CHECK(fixed_history.activated_snapshots==3);
    CHECK(fixed_history.excluded_nominal_volume_m3>0.);
    CHECK_THAT(fixed_history.thermal.energy_balance_error_j,Catch::Matchers::WithinAbs(0.,1e-9));
    scenario.use_fixed_interface_reference=false;
    // Once the interface is present, later spatially disjoint deposition is
    // valid and must retain the earlier material/thermal history.
    scenario.activation_times_s = {path.deposition[1].end_s,path.total_s};
    const auto staged = simulate_section_gcode_history(result,scenario,limits);
    REQUIRE(staged.thermal.state == ThermalHistory::State::Converged);
    REQUIRE(staged.mechanical.has_value());
    CHECK(staged.activated_snapshots == 2);
    CHECK_THAT(staged.occupied_volume_m3,Catch::Matchers::WithinAbs(solved.occupied_volume_m3,1e-18));
    CHECK_THAT(staged.thermal.final_energy_j,Catch::Matchers::WithinAbs(solved.thermal.final_energy_j,1e-9));
    CHECK_THAT(staged.thermal.energy_balance_error_j,Catch::Matchers::WithinAbs(0.,1e-9));
    result.deposition_provenance.reset();
    CHECK_THROWS_AS(simulate_section_gcode_history(result,scenario),std::invalid_argument);
}

LayerMaterialOverlap conflict_for(const GCodeProcessorResult& result)
{
    const auto path = timed_deposition(result);
    REQUIRE(path.state == ToolpathState::Complete);
    const auto nominal = resolve_nominal_deposition(result, path, path.total_s, BeadSection::RectangularEnvelope);
    try {
        section_deposition(nominal.segments, .0001, 0.);
    } catch (const LayerMaterialOverlap& conflict) {
        return conflict;
    }
    throw std::runtime_error("Fixture has no material conflict.");
}
} // namespace

TEST_CASE("Authenticated conflict ownership ignores unrelated unowned deposition without expanding masks",
          "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    const auto conflict  = conflict_for(result);
    const auto ownership = native_ownership_for_overlap(result, conflict, BeadSection::RectangularEnvelope, 0., {0, 1});
    REQUIRE_FALSE(ownership.regions.empty());
    for (const auto& region : ownership.regions)
        CHECK(diff_ex(region.footprint.native(), conflict.section_footprint).empty());
    const auto path     = timed_deposition(result);
    const auto nominal  = resolve_nominal_deposition(result, path, path.total_s, BeadSection::RectangularEnvelope);
    const auto resolved = section_deposition(nominal.segments, .0001, 0., {}, {}, &ownership);
    CHECK(resolved.resolved_overlap_area_m2 > 0.);
    auto uniform = nominal.segments;
    for (auto& bead : uniform)
        bead.material = 0;
    const auto before = section_deposition(uniform, .0001, 0.);
    ExPolygons after;
    for (const auto& region : resolved.regions)
        after.insert(after.end(), region.footprint.native().begin(), region.footprint.native().end());
    CHECK(diff_ex(after, before.regions.front().footprint.native()).empty());
    CHECK(diff_ex(before.regions.front().footprint.native(), after).empty());
}

TEST_CASE("Ownership observers receive borrowed stages and propagate cancellation", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    const auto conflict = conflict_for(result);
    std::vector<bool> stages;
    ConflictOwnershipLimits limits;
    limits.observe_masks = [&](double, bool clipped, const auto& regions) {
        REQUIRE_FALSE(regions.empty());
        stages.push_back(clipped);
    };
    const auto observed = native_ownership_for_overlap(result, conflict, BeadSection::RectangularEnvelope, 0., {0, 1}, limits);
    REQUIRE(stages.size() == 2);
    CHECK_FALSE(stages[0]);
    CHECK(stages[1]);
    REQUIRE_FALSE(observed.regions.empty());
    limits.observe_masks = [](double, bool, const auto&) { throw ResourceLimit("diagnostic cancelled"); };
    CHECK_THROWS_AS(native_ownership_for_overlap(result, conflict, BeadSection::RectangularEnvelope, 0., {0, 1}, limits), ResourceLimit);
}

TEST_CASE("Authenticated conflict ownership rejects missing contributor identity and incomplete evidence",
          "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    auto conflict = conflict_for(result);
    SECTION("Unowned contributor")
    {
        result.moves[1].deposition_source_id = 0;
        CHECK_THROWS_AS(native_ownership_for_overlap(result, conflict, BeadSection::RectangularEnvelope, 0., {0, 1}), std::invalid_argument);
    }
    SECTION("Stale registry")
    {
        auto other = std::make_shared<DepositionProvenance>();
        other->add({});
        other->add({0, 0, 0, 0, 1});
        other->seal();
        result.deposition_provenance = other;
        CHECK_THROWS_AS(native_ownership_for_overlap(result, conflict, BeadSection::RectangularEnvelope, 0., {0, 1}), std::invalid_argument);
    }
    SECTION("A third unowned contributor cannot borrow the first two sources")
    {
        const float middle           = .5f * (result.moves[1].position.x() + result.moves[3].position.x());
        result.moves[4].position.x() = middle;
        result.moves[5].position.x() = middle;
        REQUIRE(result.moves[5].deposition_source_id == 0);
        CHECK_THROWS_AS(native_ownership_for_overlap(result, conflict, BeadSection::RectangularEnvelope, 0., {0, 1}), std::invalid_argument);
    }
    SECTION("A mask outside the actual contact is not ownership evidence")
    {
        for (auto& polygon : conflict.section_footprint)
            polygon.translate(Point(100000000, 0));
        CHECK_THROWS_AS(native_ownership_for_overlap(result, conflict, BeadSection::RectangularEnvelope, 0., {0, 1}), std::invalid_argument);
    }
    SECTION("An explicit material collapse removes the claimed different-material conflict")
    {
        CHECK_THROWS_AS(native_ownership_for_overlap(result, conflict, BeadSection::RectangularEnvelope, 0., {0, 0}), std::invalid_argument);
    }
    SECTION("Partial mask")
    {
        conflict.section_footprint_complete = false;
        CHECK_THROWS_AS(native_ownership_for_overlap(result, conflict, BeadSection::RectangularEnvelope, 0., {0, 1}), std::invalid_argument);
    }
    SECTION("Budget and cancellation")
    {
        ConflictOwnershipLimits limits;
        limits.source_centerline_boundary_allowance_nm = 4;
        CHECK_THROWS_AS(native_ownership_for_overlap(result, conflict, BeadSection::RectangularEnvelope, 0., {0, 1}, limits), std::invalid_argument);
        limits.interface_variant = NativeInterfaceVariant::FavorFirst;
        limits.source_centerline_boundary_allowance_nm = 8;
        CHECK_THROWS_AS(native_ownership_for_overlap(result, conflict, BeadSection::RectangularEnvelope, 0., {0, 1}, limits), std::invalid_argument);
        limits = {};
        limits.max_polygon_checks = 1;
        CHECK_THROWS_AS(native_ownership_for_overlap(result, conflict, BeadSection::RectangularEnvelope, 0., {0, 1}, limits), ResourceLimit);
        limits                  = {};
        limits.max_polygon_work = 1;
        CHECK_THROWS_AS(native_ownership_for_overlap(result, conflict, BeadSection::RectangularEnvelope, 0., {0, 1}, limits), ResourceLimit);
        limits           = {};
        limits.max_moves = result.moves.size() - 1;
        CHECK_THROWS_AS(native_ownership_for_overlap(result, conflict, BeadSection::RectangularEnvelope, 0., {0, 1}, limits), ResourceLimit);
        CHECK_THROWS_AS(native_ownership_for_overlap(result, conflict, BeadSection::RectangularEnvelope, 0., {0, 1}, {},
                                                     [] { throw std::runtime_error("cancel"); }),
                        std::runtime_error);
    }
}

TEST_CASE("Ownership for one complete contact cannot authorize a second unexamined contact", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    const auto first_conflict = conflict_for(result);
    // A second disjoint contact has no source ownership. Both original sources
    // still suffice at the first contact, but must not authorize this one.
    const auto original_size = result.moves.size();
    for (std::size_t i = 0; i < 4; ++i) {
        auto move = result.moves[i];
        move.position.y() += 10.f;
        move.deposition_source_id = 0;
        result.moves.push_back(move);
    }
    result.print_statistics.modes[0].time += double(result.moves.size() - original_size);
    const auto ownership = native_ownership_for_overlap(result, first_conflict, BeadSection::RectangularEnvelope, 0., {0, 1});
    REQUIRE_FALSE(ownership.regions.empty());
    for (const auto& region : ownership.regions)
        CHECK(diff_ex(region.footprint.native(), first_conflict.section_footprint).empty());
    const auto path = timed_deposition(result);
    REQUIRE(path.state == ToolpathState::Complete);
    const auto nominal = resolve_nominal_deposition(result, path, path.total_s, BeadSection::RectangularEnvelope);
    CHECK_THROWS_AS(section_deposition(nominal.segments, .0001, 0., {}, {}, &ownership), std::invalid_argument);
}

TEST_CASE("Authenticated conflict ownership does not admit transformed sources", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result, true);
    const auto conflict = conflict_for(result);
    CHECK_THROWS_AS(native_ownership_for_overlap(result, conflict, BeadSection::RectangularEnvelope, 0., {0, 1}), std::invalid_argument);
}

TEST_CASE("Prepared native ownership uses only the selected model deposition prefix", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    const auto conflict = conflict_for(result);
    const float middle = .5f * (result.moves[1].position.x() + result.moves[3].position.x());
    result.moves[4].position.x() = middle;
    result.moves[5].position.x() = middle;
    REQUIRE(result.moves[5].deposition_source_id == 0);
    const auto prefix = prepare_native_model_ownership(result, BeadSection::RectangularEnvelope, 0., {0, 1}, 4.);
    REQUIRE_FALSE(prefix(conflict).regions.empty());
    CHECK_THROWS_AS(native_ownership_for_overlap(result, conflict, BeadSection::RectangularEnvelope, 0., {0, 1}),
        std::invalid_argument);

    // Model-only analysis explicitly excludes tower/custom deposition. It
    // must not reinstate it solely in the ownership evidence path.
    result.moves[5].extrusion_role = erWipeTower;
    const auto model = prepare_native_model_ownership(result, BeadSection::RectangularEnvelope, 0., {0, 1}, 6.);
    REQUIRE_FALSE(model(conflict).regions.empty());
}

TEST_CASE("Deferred ownership sees the complete actual conflict and preserves union geometry", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    const auto path = timed_deposition(result);
    const auto nominal = resolve_nominal_deposition(result, path, path.total_s, BeadSection::RectangularEnvelope);
    const auto resolver = prepare_native_model_ownership(result, BeadSection::RectangularEnvelope, 0., {0, 1}, path.total_s);
    SectionOwnership deferred{.0001, {}};
    unsigned calls = 0;
    deferred.resolve_conflict = [&](const LayerMaterialOverlap& conflict) {
        ++calls;
        REQUIRE(conflict.section_footprint_complete);
        return resolver(conflict).regions;
    };
    const auto resolved = section_deposition(nominal.segments, .0001, 0., {}, {}, &deferred);
    CHECK(calls == 1);
    CHECK(resolved.resolved_overlap_area_m2 > 0.);
    auto one_material = nominal.segments;
    for (auto& bead : one_material) bead.material = 0;
    const auto original = section_deposition(one_material, .0001, 0.);
    ExPolygons after;
    for (const auto& region : resolved.regions) after.insert(after.end(), region.footprint.native().begin(), region.footprint.native().end());
    CHECK(diff_ex(after, original.regions.front().footprint.native()).empty());
    CHECK(diff_ex(original.regions.front().footprint.native(), after).empty());
}

TEST_CASE("Prepared native ownership rejects introduced overlap stale evidence and exhausted work", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    const auto conflict = conflict_for(result);
    const auto resolver = prepare_native_model_ownership(result, BeadSection::RectangularEnvelope, 0., {0, 1}, 6.);
    auto expanded = conflict;
    expanded.section_footprint = offset_ex(expanded.section_footprint, 1.);
    CHECK_THROWS_AS(resolver(expanded), std::invalid_argument);
    auto empty_polygon = conflict;
    empty_polygon.section_footprint.front().contour.points.clear();
    CHECK_THROWS_AS(resolver(empty_polygon), std::invalid_argument);
    const auto snapshot = result.native_ownership;
    result.native_ownership.reset();
    CHECK_THROWS_AS(resolver(conflict), std::invalid_argument);
    result.native_ownership = snapshot;

    ConflictOwnershipLimits limits;
    limits.max_total_work = result.moves.size() + 4;
    CHECK_THROWS_AS(prepare_native_model_ownership(result, BeadSection::RectangularEnvelope, 0., {0, 1}, 6., {}, limits), ResourceLimit);
    CHECK_THROWS_AS(prepare_native_model_ownership(result, BeadSection::RectangularEnvelope, 0., {0, 1}, 6., {}, {},
        [] { throw std::logic_error("cancel"); }), std::logic_error);
}

TEST_CASE("Strict prepared ownership rejects boundary assumptions before resolving a contact", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    ConflictOwnershipLimits limits;
    limits.allow_single_sided_native_boundary = true;
    CHECK_THROWS_AS(prepare_native_model_ownership(result, BeadSection::RectangularEnvelope,
        0., {0, 1}, 6., {}, limits), std::invalid_argument);
    limits.allow_single_sided_native_boundary = false;
    limits.allow_dual_native_boundary = true;
    limits.geometry.retain_exact_interface_sections = true;
    CHECK_THROWS_AS(prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,{0,1},6.,{},limits),std::invalid_argument);
    limits.interface_variant = NativeInterfaceVariant::FavorFirst;
    limits.geometry.retain_exact_interface_sections = false;
    CHECK_THROWS_AS(prepare_native_model_ownership(result,BeadSection::RectangularEnvelope,0.,{0,1},6.,{},limits),std::invalid_argument);
    limits.geometry.retain_exact_interface_sections = true;
    CHECK_THROWS_AS(native_ownership_for_overlap(result,conflict_for(result),BeadSection::RectangularEnvelope,0.,{0,1},limits),std::invalid_argument);
    limits.allow_dual_native_boundary = false;
    limits.interface_variant = NativeInterfaceVariant::Strict;
    limits.source_centerline_boundary_allowance_nm = 4;
    CHECK_THROWS_AS(prepare_native_model_ownership(result, BeadSection::RectangularEnvelope,
        0., {0, 1}, 6., {}, limits), std::invalid_argument);
}

TEST_CASE("Foreign envelope expansion cannot displace authenticated uniquely occupied original material", "[ThermoelasticConflictOwnership]")
{
    GCodeProcessorResult result;
    populate(result);
    const auto path = timed_deposition(result);
    const auto nominal = resolve_nominal_deposition(result, path, path.total_s, BeadSection::RectangularEnvelope);
    const auto center = .5 * (nominal.segments[0].start_m + nominal.segments[0].end_m);
    const coord_t x = std::llround(center.x() * 1e9), y = std::llround(center.y() * 1e9);
    // Exercise ownership on the exact straight-sided interior. Curved caps
    // require integer Boolean rounding and must retain the strict union guard.
    const ExPolygon clip(Polygon({Point(x - 1000000, y - 200000), Point(x + 1000000, y - 200000),
                                  Point(x + 1000000, y + 200000), Point(x - 1000000, y + 200000)}));
    const auto a = section_deposition({nominal.segments[0]}, .0001, 0., {}, {}, nullptr, &clip);
    const auto b = section_deposition({nominal.segments[1]}, .0001, 0., {}, {}, nullptr, &clip);
    const auto unique_a = diff_ex(a.regions.front().footprint.native(), b.regions.front().footprint.native());
    REQUIRE_FALSE(unique_a.empty());
    const auto resolver = prepare_native_model_ownership(result, BeadSection::RectangularEnvelope, 0., {0, 1}, 6., {}, {}, {}, true);
    // A widened foreign nominal envelope is a simple fixture for added
    // conflict. Original authenticated paths/dimensions remain unchanged.
    auto expanded = nominal.segments;
    expanded[1].width_m = .0008;
    SectionOwnership deferred{.0001, {}};
    deferred.resolve_conflict = [&](const LayerMaterialOverlap& conflict) { return resolver(conflict).regions; };
    CHECK_THROWS_AS(section_deposition(expanded, .0001, 0., {}, {}, &deferred), std::invalid_argument);
    const auto resolved = section_deposition(expanded, .0001, 0., {}, {}, &deferred, &clip);
    ExPolygons original_owner, after;
    for (const auto& region : resolved.regions) {
        if (region.material == 0) original_owner = region.footprint.native();
        after.insert(after.end(), region.footprint.native().begin(), region.footprint.native().end());
    }
    CHECK(diff_ex(unique_a, original_owner).empty());
    auto uniform = expanded;
    for (auto& bead : uniform) bead.material = 0;
    const auto union_geometry = section_deposition(uniform, .0001, 0., {}, {}, nullptr, &clip);
    CHECK(diff_ex(union_geometry.regions.front().footprint.native(), after).empty());
    CHECK(diff_ex(after, union_geometry.regions.front().footprint.native()).empty());

    auto empty_space = conflict_for(result);
    for (auto& polygon : empty_space.section_footprint) polygon.translate(Point(100000000, 0));
    CHECK_THROWS_AS(resolver(empty_space), std::invalid_argument);
    result.moves[1].deposition_source_id = 0;
    const auto unowned = prepare_native_model_ownership(result, BeadSection::RectangularEnvelope, 0., {0, 1}, 6., {}, {}, {}, true);
    CHECK_THROWS_AS(unowned(conflict_for(result)), std::invalid_argument);
}

TEST_CASE("Sampler minted conflicts authenticate the complete unchanged observation prefix", "[ThermoelasticConflictOwnership]")
{
    STATIC_REQUIRE_FALSE(std::is_default_constructible_v<SampledSectionConflict>);
    STATIC_REQUIRE_FALSE(std::is_copy_constructible_v<SampledSectionConflict>);
    GCodeProcessorResult result;
    populate(result);
    const auto path = timed_deposition(result);
    auto nominal = resolve_nominal_deposition(result, path, path.total_s, BeadSection::RectangularEnvelope);
    for (std::size_t i = 0; i < nominal.segments.size(); ++i) nominal.segments[i].section_source_index = i;
    SampledConflictOwnershipResolver sampled;
    prepare_native_model_ownership(result, BeadSection::RectangularEnvelope, 0., {0, 1}, path.total_s,
        {}, {}, {}, false, &sampled);
    SectionOwnership owner{.0001, {}};
    owner.resolve_sampled_conflict = [&](const SampledSectionConflict& conflict) { return sampled(conflict).regions; };
    SECTION("Complete geometry includes unrelated unowned brim without assigning its ownership") {
        const auto section = section_deposition(nominal.segments, .0001, 0., {}, {}, &owner);
        CHECK(section.resolved_overlap_area_m2 > 0.);
    }
    SECTION("Forged geometry with a valid index is refused") {
        nominal.segments.front().start_m.x() += 1e-9;
        CHECK_THROWS_AS(section_deposition(nominal.segments, .0001, 0., {}, {}, &owner), std::invalid_argument);
    }
    SECTION("A shorter prefix cannot borrow a longer observation certificate") {
        nominal.segments.pop_back();
        CHECK_THROWS_AS(section_deposition(nominal.segments, .0001, 0., {}, {}, &owner), std::invalid_argument);
    }
    SECTION("A later prefix cannot borrow an earlier observation certificate") {
        prepare_native_model_ownership(result, BeadSection::RectangularEnvelope, 0., {0, 1}, 4.,
            {}, {}, {}, false, &sampled);
        CHECK_THROWS_AS(section_deposition(nominal.segments, .0001, 0., {}, {}, &owner), std::invalid_argument);
    }
    SECTION("A stale source generation cannot use sampler evidence") {
        result.deposition_provenance.reset();
        CHECK_THROWS_AS(section_deposition(nominal.segments, .0001, 0., {}, {}, &owner), std::invalid_argument);
    }
    SECTION("Mismatched section approximation cannot use sampler evidence") {
        LayerGeometryLimits geometry;
        geometry.arc_tolerance_m *= 2.;
        CHECK_THROWS_AS(section_deposition(nominal.segments, .0001, 0., geometry, {}, &owner), std::invalid_argument);
    }
    SECTION("A different query-plane precision cannot borrow sampler evidence") {
        LayerGeometryLimits geometry;
        geometry.half_native_section_planes = true;
        CHECK_THROWS_AS(section_deposition(nominal.segments, .0001, 0., geometry, {}, &owner), std::invalid_argument);
    }
    SECTION("Authenticated half-native source queries preserve plane identity") {
        ConflictOwnershipLimits limits;
        limits.geometry.half_native_section_planes = true;
        prepare_native_model_ownership(result, BeadSection::RectangularEnvelope, 0., {0,1}, path.total_s,
            {}, limits, {}, false, &sampled);
        owner.z_m = .0001000005;
        const auto section = section_deposition(nominal.segments, owner.z_m, 0., limits.geometry, {}, &owner);
        CHECK(section.resolved_overlap_area_m2 > 0.);
        owner.z_m = .0001;
        CHECK_THROWS_AS(section_deposition(nominal.segments, .0001000005, 0., limits.geometry, {}, &owner), std::invalid_argument);
    }
    SECTION("A third unowned contributing bead cannot be omitted") {
        result.moves[4].position = result.moves[0].position;
        result.moves[5].position = result.moves[1].position;
        prepare_native_model_ownership(result, BeadSection::RectangularEnvelope, 0., {0, 1}, path.total_s,
            {}, {}, {}, false, &sampled);
        nominal.segments.pop_back();
        CHECK_THROWS_AS(section_deposition(nominal.segments, .0001, 0., {}, {}, &owner), std::invalid_argument);
    }
}
