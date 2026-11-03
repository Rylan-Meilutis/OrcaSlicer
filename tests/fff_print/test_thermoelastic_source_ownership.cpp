#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <cstddef>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/ThermoelasticLayers.hpp"
#include <Eigen/Core>
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/GCode/DepositionProvenance.hpp"
#include "libslic3r/GCode/DepositionSection.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Surface.hpp"
#include "libslic3r/SurfaceCollection.hpp"
#include "libslic3r/ThermoelasticSourceOwnership.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "test_helpers.hpp"
#include <cstdint>
#include <array>
#include <memory>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;

TEST_CASE("Unowned planar brim requires its own typed live emission witness", "[ThermoelasticSourceOwnership]")
{
    const int mutation=GENERATE(-1,0,1,2,3,4,5,6,7,8,9);
    auto registry=std::make_shared<DepositionProvenance>();
    const auto model_id=registry->add({1,0,0,0,0,DepositionSource::Height::NativePlanar,
        DepositionSource::Assignment::NominalRole});
    registry->seal();
    auto sections=std::make_shared<DepositionSections>(registry->generation());
    GCodeProcessorResult result;
    result.deposition_provenance=registry;
    result.moves.resize(2);
    auto& first=result.moves[0];
    auto& move=result.moves[1];
    first.type=EMoveType::Travel;
    first.position={1.f,2.f,.2f};
    move.type=EMoveType::Extrude;
    move.position={2.f,2.f,.2f};
    move.extrusion_role=erBrim;
    move.delta_extruder=.03618f;
    move.width=.4f;
    move.height=.2f;
    DepositionCommandSection command{{1.,2.,double(.2f)},{3.,2.,double(.2f)},.4f,.2f,.2f,0};
    command.purpose=mutation==0 ? DepositionCommandSection::Purpose::Model : DepositionCommandSection::Purpose::AuxiliaryBrim;
    move.deposition_command_section_id=sections->add(command);
    move.deposition_design_heights=std::array<double,2>{double(.2f),double(.2f)};
    move.deposition_design_width=.4f;
    sections->seal();
    result.deposition_sections=sections;
    SourceBead bead{{first.position.cast<double>()*.001,move.position.cast<double>()*.001,
        double(.4f)*.001,double(.2f)*.001,7},0,0,1};
    switch (mutation) {
    case 1: bead.source_id=move.deposition_source_id=model_id; break;
    case 2: move.extrusion_role=erPerimeter; break;
    case 3:
        move.position.z()=.201f;
        bead.geometry.end_m=move.position.cast<double>()*.001;
        break;
    case 4:
        move.deposition_design_heights=std::array<double,2>{.1,.1};
        bead.geometry.height_m=.0001;
        break;
    case 5: {
        auto stale=std::make_shared<DepositionSections>(registry->generation()+1);
        stale->add(command);
        stale->seal();
        result.deposition_sections=stale;
        break;
    }
    case 6: move.deposition_command_section_id=0; break;
    case 7: move.deposition_design_width=.5f; break;
    case 8: move.delta_extruder=0.f; break;
    case 9: bead.logical_filament=1; break;
    default: break;
    }
    if (mutation>=0) {
        CHECK_THROWS_AS(authenticate_live_deposited_source(result,bead),std::invalid_argument);
        return;
    }
    CHECK(bead.source_id==0);
    CHECK_FALSE(result.native_ownership);
    const auto work=authenticate_live_deposited_source(result,bead);
    CHECK(work>1);
    CHECK(authenticate_live_deposited_source(result,bead,work)==work);
    CHECK_THROWS_AS(authenticate_live_deposited_source(result,bead,work-1),ResourceLimit);
    CHECK_THROWS_AS(authenticate_live_deposited_source(result,bead,100,
        [] { throw std::runtime_error("cancel auxiliary brim"); }),std::runtime_error);
}

TEST_CASE("Live source shape authentication is independent of native masks and rejects stale witnesses", "[ThermoelasticSourceOwnership]")
{
    const bool typed=GENERATE(false,true);
    auto registry=std::make_shared<DepositionProvenance>();
    const auto id=registry->add({1,0,0,0,0,typed ? DepositionSource::Height::Transformed : DepositionSource::Height::NativePlanar,
        DepositionSource::Assignment::NominalRole});
    const auto other_id=registry->add({2,0,0,0,0,DepositionSource::Height::NativePlanar,
        DepositionSource::Assignment::NominalRole});
    registry->seal();
    auto sections=std::make_shared<DepositionSections>(registry->generation());
    GCodeProcessorResult result;
    result.deposition_provenance=registry;
    result.moves.resize(2);
    auto& first=result.moves[0];
    auto& move=result.moves[1];
    first.type=EMoveType::Travel;
    first.position={1.f,2.f,.2f};
    move.type=EMoveType::Extrude;
    move.position={2.f,2.f,typed ? .3f : .2f};
    move.extrusion_role=erPerimeter;
    move.delta_extruder=.1f;
    move.width=.4f;
    move.height=.2f;
    move.deposition_source_id=id;
    if (typed) {
        DepositionCommandSection command{{first.position.x(),first.position.y(),first.position.z()},
            {move.position.x(),move.position.y(),move.position.z()},move.width,move.height,move.height,0};
        command.transform=DepositionCommandSection::Transform::BrickCourse;
        move.deposition_command_section_id=sections->add(command);
        move.deposition_design_heights=std::array<double,2>{move.height,move.height};
        move.deposition_design_width=move.width;
    }
    sections->seal();
    result.deposition_sections=sections;
    SourceBead bead{{first.position.cast<double>()*.001,move.position.cast<double>()*.001,
        double(move.width)*.001,double(move.height)*.001,7},id,0,1};
    // No native snapshot exists; this proves emitted shape, not CAD ownership.
    const auto work=authenticate_live_deposited_source(result,bead);
    CHECK(work>0);
    CHECK(authenticate_live_deposited_source(result,bead,work)==work);
    CHECK_THROWS_AS(authenticate_live_deposited_source(result,bead,work-1),ResourceLimit);
    CHECK_THROWS_AS(authenticate_live_deposited_source(result,bead,100,
        [] { throw std::runtime_error("cancel live shape"); }),std::runtime_error);
    const int mutation=GENERATE(0,1,2,3,4,5,6,7,8,9,10,11);
    switch (mutation) {
    case 0: move.position.x()+=.01f; break;
    case 1: bead.source_id=other_id; break; // Same filament is insufficient.
    case 2: bead.geometry.width_m*=2.; break;
    case 3: move.delta_extruder=0.f; break;
    case 4: move.deposition_command_section_id=999; break;
    case 5: move.extrusion_role=erMixed; break;
    case 6: move.deposition_vertical_extrusion=true; break;
    case 7: move.deposition_design_heights=std::array<double,2>{.1,.1}; break;
    case 8: result.deposition_sections=std::make_shared<DepositionSections>(registry->generation()+1); break;
    case 9: bead.logical_filament=1; break;
    case 10: move.type=EMoveType::Travel; break;
    case 11: bead.geometry.end_m.x()-=.0001; break;
    }
    CHECK_THROWS_AS(authenticate_live_deposited_source(result,bead),std::invalid_argument);
}

TEST_CASE("Planar infill in a transformed region requires live native-plane evidence", "[ThermoelasticSourceOwnership]")
{
    const bool authenticated = GENERATE(false, true);
    const auto assignment=GENERATE(DepositionSource::Assignment::Unknown,DepositionSource::Assignment::NominalRole,
                                   DepositionSource::Assignment::WipingOverride);
    for (const std::string option : {"", "flush_into_infill", "nonplanar_infill", "infill_combination", "zaa_enabled", "top_surface_z_mode"}) {
        DYNAMIC_SECTION("Infill setting " << option) {
            auto config=Test::multifilament_config(1,{{"layer_height","0.2"},{"initial_layer_print_height","0.2"},
                {"seam_slope_type","all"},{"brim_type","no_brim"},{"skirt_loops","0"}});
            if (option=="top_surface_z_mode") config.set_deserialize_strict({{"top_surface_z_mode","z_contouring"}});
            else if (!option.empty()) config.set(option,true);
            Print print;
            Model model;
            Test::init_print({Test::cube(2.)},print,model,config);
            print.process();
            const auto* object=print.get_object(0);
            const auto* layer=object->layers().front();
            const auto* region=layer->regions().front();
            auto registry=std::make_shared<DepositionProvenance>();
            const auto id=registry->add({object->get_id(),0,std::uint64_t(layer->id()),
                std::uint64_t(region->region().print_region_id()),0,DepositionSource::Height::Transformed,assignment});
            registry->seal();
            auto snapshot=std::make_shared<const NativeOwnershipSnapshot>(capture_native_ownership(print,registry,Point(0,0),0.));
            ExPolygons native;
            for (const auto& surface : region->slices.surfaces) {
                native.push_back(surface.expolygon);
                native.back().translate(object->instances()[0].shift);
            }
            const auto box=get_extents(native);
            const float x=float(.5*double(box.min.x()+box.max.x())*1e-6);
            const float y=float(.5*double(box.min.y()+box.max.y())*1e-6);
            GCodeProcessorResult result;
            result.deposition_provenance=registry;
            result.native_ownership=snapshot;
            auto sections=std::make_shared<DepositionSections>(registry->generation());
            GCodeProcessorResult::MoveVertex first{},last{};
            first.position={x-.1f,y,float(layer->print_z)};
            first.type=EMoveType::Travel;
            last.position={x+.1f,y,float(layer->print_z)};
            last.type=EMoveType::Extrude;
            last.extrusion_role=erTopSolidInfill;
            last.deposition_source_id=id;
            last.width=.4f;
            last.height=float(layer->height);
            if (authenticated) {
                last.deposition_command_section_id=sections->add({
                    {first.position.x(),first.position.y(),first.position.z()},
                    {last.position.x(),last.position.y(),last.position.z()},last.width,last.height,last.height,0});
                last.deposition_design_heights=std::array<double,2>{last.height,last.height};
                last.deposition_design_width=last.width;
            }
            sections->seal();
            result.deposition_sections=sections;
            result.moves={first,last};
            SourceBead bead{{first.position.cast<double>()*.001,last.position.cast<double>()*.001,
                double(last.width)*.001,double(last.height)*.001,0},id,0,1};
            const double plane=(layer->print_z-.5*layer->height)*.001;
            const bool eligible=assignment!=DepositionSource::Assignment::WipingOverride &&
                (option.empty() || (assignment==DepositionSource::Assignment::NominalRole && authenticated &&
                                    option!="nonplanar_infill"));
            if (!eligible) {
                CHECK_THROWS_AS(snapshot->section(registry,{bead},plane,{0},{},{},&result),std::invalid_argument);
                continue;
            }
            const auto admitted=snapshot->section(registry,{bead},plane,{0},{},{},&result);
            REQUIRE(admitted.regions.size()==1);
            CHECK(diff_ex(admitted.regions.front().footprint.native(),native).empty());
            CHECK(diff_ex(native,admitted.regions.front().footprint.native()).empty());
            CHECK_THROWS_AS(snapshot->section(registry,{bead},plane,{0}),std::invalid_argument);
            if (authenticated) {
                auto wrong_generation=std::make_shared<DepositionSections>(registry->generation()+1);
                wrong_generation->seal();
                result.deposition_sections=wrong_generation;
                CHECK_THROWS_AS(snapshot->section(registry,{bead},plane,{0},{},{},&result),std::invalid_argument);
                result.deposition_sections=sections;
                result.moves[1].deposition_command_section_id=sections->size()+1;
                CHECK_THROWS_AS(snapshot->section(registry,{bead},plane,{0},{},{},&result),std::invalid_argument);
                result.moves[1]=last;
                // A valid nominal-pass command still cannot certify a
                // combined-height infill bead or a genuinely raised path.
                for (const bool raised : {false,true}) {
                    auto altered=last;
                    auto altered_sections=std::make_shared<DepositionSections>(registry->generation());
                    const float end_height=raised ? last.height+.01f : 2.f*last.height;
                    const float start_height=raised ? last.height : end_height;
                    if (raised) altered.position.z()+=.01f;
                    altered.deposition_command_section_id=altered_sections->add({
                        {first.position.x(),first.position.y(),first.position.z()},
                        {altered.position.x(),altered.position.y(),altered.position.z()},
                        last.width,start_height,end_height,0});
                    altered.deposition_design_heights=std::array<double,2>{start_height,end_height};
                    altered_sections->seal();
                    result.deposition_sections=altered_sections;
                    result.moves[1]=altered;
                    auto altered_bead=bead;
                    altered_bead.geometry.end_m=altered.position.cast<double>()*.001;
                    altered_bead.geometry.height_m=double(start_height)*.001;
                    altered_bead.geometry.end_height_m=double(end_height)*.001;
                    CHECK_THROWS_AS(snapshot->section(registry,{altered_bead},plane,{0},{},{},&result),std::invalid_argument);
                }
                result.deposition_sections=sections;
                result.moves[1]=last;
            }
            auto wrong=bead;
            wrong.geometry.width_m+=.00001;
            CHECK_THROWS_AS(snapshot->section(registry,{wrong},plane,{0},{},{},&result),std::invalid_argument);
            wrong=bead;
            wrong.geometry.start_m.y()+=.00001;
            CHECK_THROWS_AS(snapshot->section(registry,{wrong},plane,{0},{},{},&result),std::invalid_argument);
            wrong=bead;
            wrong.geometry.end_m.y()+=.00001;
            CHECK_THROWS_AS(snapshot->section(registry,{wrong},plane,{0},{},{},&result),std::invalid_argument);
            wrong=bead;
            wrong.geometry.end_height_m=wrong.geometry.height_m*.5;
            CHECK_THROWS_AS(snapshot->section(registry,{wrong},plane,{0},{},{},&result),std::invalid_argument);
            result.moves[1].extrusion_role=erBridgeInfill;
            CHECK_THROWS_AS(snapshot->section(registry,{bead},plane,{0},{},{},&result),std::invalid_argument);
            result.moves[1]=last;
            result.moves[1].position.z()+=.01f;
            bead.geometry.end_m=result.moves[1].position.cast<double>()*.001;
            CHECK_THROWS_AS(snapshot->section(registry,{bead},plane,{0},{},{},&result),std::invalid_argument);
        }
    }
}

TEST_CASE("Transformed source masks require authenticated contained scarf commands", "[ThermoelasticSourceOwnership]")
{
    const auto config = Test::multifilament_config(1, {{"layer_height", "0.2"},
                                                       {"initial_layer_print_height", "0.2"},
                                                       {"flush_into_infill", "1"},
                                                       {"brim_type", "no_brim"},
                                                       {"skirt_loops", "0"}});
    Print print;
    Model model;
    Test::init_print({Test::cube(2.)}, print, model, config);
    print.process();
    const auto* object = print.get_object(0);
    const auto* layer  = object->layers().front();
    const auto* region = layer->regions().front();
    auto registry      = std::make_shared<DepositionProvenance>();
    const auto id = registry->add({object->get_id(), 0, std::uint64_t(layer->id()), std::uint64_t(region->region().print_region_id()), 0,
                                   DepositionSource::Height::Transformed,DepositionSource::Assignment::NominalRole});
    registry->seal();
    auto snapshot = std::make_shared<const NativeOwnershipSnapshot>(capture_native_ownership(print, registry, Point(0, 0), 0.));
    CHECK(snapshot->source_status(id) == NativeOwnershipSnapshot::SourceStatus::Unsupported);
    ExPolygons native;
    for (const auto& surface : region->slices.surfaces) {
        native.push_back(surface.expolygon);
        native.back().translate(object->instances()[0].shift);
    }
    const auto box = get_extents(native);
    const double x = .5 * double(box.min.x() + box.max.x()) * 1e-6;
    const double y = .5 * double(box.min.y() + box.max.y()) * 1e-6;
    GCodeProcessorResult result;
    result.deposition_provenance = registry;
    result.native_ownership      = snapshot;
    auto sections                = std::make_shared<DepositionSections>(registry->generation());
    const auto command           = sections->add({{float(x - .1), float(y), .15f}, {float(x + .1), float(y), .15f}, .4, .05, .1, 0});
    sections->seal();
    result.deposition_sections = sections;
    GCodeProcessorResult::MoveVertex first{}, last{};
    first.position                     = Vec3f(float(x - .1), float(y), .15f);
    first.type                         = EMoveType::Travel;
    last.position                      = Vec3f(float(x + .1), float(y), .15f);
    last.type                          = EMoveType::Extrude;
    last.extrusion_role                = erExternalPerimeter;
    last.deposition_source_id          = id;
    last.deposition_command_section_id = command;
    last.deposition_design_heights     = std::array<double, 2>{.05f, .1f};
    last.deposition_design_width       = .4f;
    result.moves                       = {first, last};
    SourceBead bead{{first.position.cast<double>() * .001, last.position.cast<double>() * .001, double(.4f) * .001, double(.05f) * .001, 0},
                    id,
                    0,
                    1};
    bead.geometry.end_height_m = double(.1f) * .001;
    std::size_t source_work = 0;
    REQUIRE(snapshot->section(registry, {bead}, .00012, {0}, {}, {}, &result, &source_work).regions.size() == 1);
    CHECK(source_work > 0);
    CHECK(source_work <= SourceOwnershipLimits{}.max_work);
    {
        auto auxiliary_sections=std::make_shared<DepositionSections>(registry->generation());
        DepositionCommandSection auxiliary{{first.position.x(),first.position.y(),first.position.z()},
            {last.position.x(),last.position.y(),last.position.z()},.4f,.05f,.05f,0};
        auxiliary.purpose=DepositionCommandSection::Purpose::AuxiliaryBrim;
        const auto auxiliary_id=auxiliary_sections->add(auxiliary);
        auxiliary_sections->seal();
        result.deposition_sections=auxiliary_sections;
        result.moves[1].deposition_command_section_id=auxiliary_id;
        result.moves[1].deposition_design_heights=std::array<double,2>{double(.05f),double(.05f)};
        result.moves[1].extrusion_role=erBrim;
        result.moves[1].delta_extruder=.01f;
        auto auxiliary_bead=bead;
        auxiliary_bead.geometry.end_height_m=auxiliary_bead.geometry.height_m;
        // Even a real same-filament CAD source cannot adopt an auxiliary parent.
        CHECK_THROWS_AS(authenticate_live_deposited_source(result,auxiliary_bead),std::invalid_argument);
        CHECK_THROWS_AS(snapshot->section(registry,{auxiliary_bead},.00012,{0},{},{},&result),std::invalid_argument);
        auxiliary_bead.source_id=0;
        result.moves[1].deposition_source_id=0;
        CHECK(authenticate_live_deposited_source(result,auxiliary_bead)>0);
        CHECK_THROWS_AS(snapshot->section(registry,{auxiliary_bead},.00012,{0},{},{},&result),std::invalid_argument);
        result.deposition_sections=sections;
        result.moves[1]=last;
    }
    {
        const auto cache = snapshot->make_containment_cache();
        std::size_t cold_work = 0, warm_work = 0;
        const auto cold = snapshot->section(registry, {bead}, .00012, {0}, {}, {}, &result, &cold_work, cache.get());
        const auto warm = snapshot->section(registry, {bead}, .00013, {0}, {}, {}, &result, &warm_work, cache.get());
        CHECK(warm_work < cold_work);
        REQUIRE(cold.regions.size() == 1);
        REQUIRE(warm.regions.size() == 1);
        CHECK(warm.regions.front().footprint.native() == cold.regions.front().footprint.native());
        CHECK_THROWS_AS(snapshot->section(registry, {bead}, .001, {0}, {}, {}, &result, nullptr, cache.get()), std::invalid_argument);
        result.moves[1].deposition_command_section_id = command + 1;
        CHECK_THROWS_AS(snapshot->section(registry, {bead}, .00012, {0}, {}, {}, &result, nullptr, cache.get()), std::invalid_argument);
        result.moves[1].deposition_command_section_id = command;
        result.moves[1].extrusion_role = erInternalInfill;
        CHECK_THROWS_AS(snapshot->section(registry, {bead}, .00012, {0}, {}, {}, &result, nullptr, cache.get()), std::invalid_argument);
        result.moves[1].extrusion_role = erExternalPerimeter;
        result.moves[1].position.y() += .1f;
        CHECK_THROWS_AS(snapshot->section(registry, {bead}, .00012, {0}, {}, {}, &result, nullptr, cache.get()), std::invalid_argument);
        result.moves[1].position = last.position;
        auto wrong_generation = std::make_shared<DepositionSections>(registry->generation() + 1);
        wrong_generation->seal();
        result.deposition_sections = wrong_generation;
        CHECK_THROWS_AS(snapshot->section(registry, {bead}, .00012, {0}, {}, {}, &result, nullptr, cache.get()), std::invalid_argument);
        result.deposition_sections = sections;
        SourceOwnershipLimits limited;
        const auto interrupted_cache = snapshot->make_containment_cache();
        limited.max_work = cold_work - 1;
        CHECK_THROWS_AS(snapshot->section(registry, {bead}, .00012, {0}, limited, {}, &result, nullptr,
            interrupted_cache.get()), ResourceLimit);
        std::size_t retried_work = 0;
        snapshot->section(registry, {bead}, .00012, {0}, {}, {}, &result, &retried_work, interrupted_cache.get());
        CHECK(retried_work == cold_work);
        limited.max_work = warm_work - 1;
        CHECK_THROWS_AS(snapshot->section(registry, {bead}, .00012, {0}, limited, {}, &result, nullptr, cache.get()), ResourceLimit);
        CHECK_THROWS_AS(snapshot->section(registry, {bead}, .00012, {0}, {},
            [] { throw std::runtime_error("cancel"); }, &result, nullptr, cache.get()), std::runtime_error);
        const auto other_snapshot = std::make_shared<const NativeOwnershipSnapshot>(capture_native_ownership(print, registry, Point(0, 0), 0.));
        result.native_ownership = other_snapshot;
        REQUIRE(other_snapshot->section(registry, {bead}, .00012, {0}, {}, {}, &result).regions.size() == 1);
        CHECK_THROWS_AS(other_snapshot->section(registry, {bead}, .00012, {0}, {}, {}, &result, nullptr, cache.get()), std::invalid_argument);
        result.native_ownership = snapshot;
    }
    {
        // The analysis uses exact decimal tops, while both display vertices
        // remain float witnesses of that same authenticated leading ramp.
        auto anchored_sections = std::make_shared<DepositionSections>(registry->generation());
        DepositionCommandSection anchored{{float(x - .1), float(y), .1},
            {float(x + .1), float(y), .15}, .4, .1, .15, 0};
        anchored.nominal_lower_z = std::array<double, 2>{0., 0.};
        const auto anchored_id = anchored_sections->add(anchored);
        anchored_sections->seal();
        result.deposition_sections = anchored_sections;
        result.moves[0].position.z() = .1f;
        result.moves[1].deposition_command_section_id = anchored_id;
        result.moves[1].deposition_design_heights = std::array<double, 2>{.1, .15};
        result.deposition_z_generation = registry->generation();
        result.deposition_z_spans.push_back({1, id, anchored_id, {.1, .15}, {.1f, .15f}});
        auto precise_bead = bead;
        precise_bead.geometry.start_m.z() = .1 * .001;
        precise_bead.geometry.end_m.z() = .15 * .001;
        precise_bead.geometry.height_m = .1 * .001;
        precise_bead.geometry.end_height_m = .15 * .001;
        CHECK(snapshot->section(registry, {precise_bead}, .00008, {0}, {}, {}, &result).regions.size() == 1);
        ++result.deposition_z_spans[0].source_id;
        CHECK_THROWS_AS(snapshot->section(registry, {precise_bead}, .00008, {0}, {}, {}, &result), std::invalid_argument);
        result.deposition_z_spans[0].source_id = id;
        precise_bead.geometry.start_m = result.moves[0].position.cast<double>() * .001;
        CHECK_THROWS_AS(snapshot->section(registry, {precise_bead}, .00008, {0}, {}, {}, &result), std::invalid_argument);
        result.deposition_z_spans.clear();
        result.deposition_z_generation = 0;
        result.deposition_sections = sections;
        result.moves = {first, last};
    }
    CHECK_THROWS_AS(snapshot->section(registry, {bead}, .00012, {0}), std::invalid_argument);
    result.moves[1].deposition_command_section_id = command + 1;
    CHECK_THROWS_AS(snapshot->section(registry, {bead}, .00012, {0}, {}, {}, &result), std::invalid_argument);
    result.moves[1].deposition_command_section_id = command;
    result.moves[1].extrusion_role                = erInternalInfill;
    CHECK_THROWS_AS(snapshot->section(registry, {bead}, .00012, {0}, {}, {}, &result), std::invalid_argument);
    result.moves[1].extrusion_role = erExternalPerimeter;
    auto wrong_generation          = std::make_shared<DepositionSections>(registry->generation() + 1);
    wrong_generation->seal();
    result.deposition_sections = wrong_generation;
    CHECK_THROWS_AS(snapshot->section(registry, {bead}, .00012, {0}, {}, {}, &result), std::invalid_argument);
    result.deposition_sections = sections;
    // Geometry can remain inside the native mask while falsely pointing at a
    // different command span. Parent membership must reject this too.
    result.moves[1].position.y() += .1f;
    bead.geometry.end_m = result.moves[1].position.cast<double>() * .001;
    CHECK_THROWS_AS(snapshot->section(registry, {bead}, .00012, {0}, {}, {}, &result), std::invalid_argument);
    result.moves[1].position = last.position;
    bead.geometry.end_m = last.position.cast<double>() * .001;
    // Even matching optional heights cannot move a scarf outside the native XY mask.
    result.moves[0].position.x() += 100.f;
    result.moves[1].position.x() += 100.f;
    bead.geometry.start_m = result.moves[0].position.cast<double>() * .001;
    bead.geometry.end_m   = result.moves[1].position.cast<double>() * .001;
    CHECK_THROWS_AS(snapshot->section(registry, {bead}, .00012, {0}, {}, {}, &result), std::invalid_argument);
    // An explicit numerical source-band study may admit a2nm centerline
    // excursion, but must leave the original native mask exactly unchanged.
    const Point offset(-box.min.x(),0);
    auto boundary_snapshot = std::make_shared<const NativeOwnershipSnapshot>(capture_native_ownership(print,registry,offset,0.));
    result.native_ownership = boundary_snapshot;
    auto boundary_sections = std::make_shared<DepositionSections>(registry->generation());
    const auto boundary_command = boundary_sections->add({{-.000002f,float(y-.1),.15f},
        {-.000002f,float(y+.1),.15f},.4,.05,.1,0});
    boundary_sections->seal();
    result.deposition_sections = boundary_sections;
    result.moves = {first,last};
    result.moves[0].position = Vec3f(-.000002f,float(y-.1),.15f);
    result.moves[1].position = Vec3f(-.000002f,float(y+.1),.15f);
    result.moves[1].deposition_command_section_id = boundary_command;
    bead.geometry.start_m = result.moves[0].position.cast<double>()*.001;
    bead.geometry.end_m = result.moves[1].position.cast<double>()*.001;
    CHECK_THROWS_AS(boundary_snapshot->section(registry,{bead},.00012,{0},{},{},&result),std::invalid_argument);
    SourceOwnershipLimits sensitivity;
    sensitivity.centerline_boundary_allowance_nm = 1;
    CHECK_THROWS_AS(boundary_snapshot->section(registry,{bead},.00012,{0},sensitivity,{},&result),std::invalid_argument);
    sensitivity.centerline_boundary_allowance_nm = 4;
    const auto admitted = boundary_snapshot->section(registry,{bead},.00012,{0},sensitivity,{},&result);
    const auto boundary_cache = boundary_snapshot->make_containment_cache();
    REQUIRE(boundary_snapshot->section(registry, {bead}, .00012, {0}, sensitivity, {}, &result, nullptr,
        boundary_cache.get()).regions.size() == 1);
    CHECK_THROWS_AS(boundary_snapshot->section(registry, {bead}, .00012, {0}, {}, {}, &result, nullptr,
        boundary_cache.get()), std::invalid_argument);
    sensitivity.centerline_boundary_allowance_nm = 1;
    CHECK_THROWS_AS(boundary_snapshot->section(registry, {bead}, .00012, {0}, sensitivity, {}, &result, nullptr,
        boundary_cache.get()), std::invalid_argument);
    REQUIRE(admitted.regions.size() == 1);
    for (auto& polygon : native) polygon.translate(offset);
    CHECK(diff_ex(admitted.regions.front().footprint.native(),native).empty());
    CHECK(diff_ex(native,admitted.regions.front().footprint.native()).empty());
}

TEST_CASE("Authenticated brick commands select actual native intervals without borrowing cached masks", "[ThermoelasticSourceOwnership]")
{
    //0 unchanged,1 nearby different upper mask,2 unsupported upper XY,
    //3 missing upper,4 missing lower. All are immutable after capture.
    const int mode=GENERATE(0,1,2,3,4);
    const auto config=Test::multifilament_config(1,{{"layer_height","0.2"},{"initial_layer_print_height","0.2"},
        {"perimeter_layering","brick"},{"staggered_perimeter_offset","50%"},{"brim_type","no_brim"},{"skirt_loops","0"}});
    Print print;
    Model model;
    Test::init_print({Test::cube(2.)},print,model,config);
    print.process();
    auto* object=print.get_object(0);
    REQUIRE(object->layers().size()>=3);
    auto* layer=object->layers()[1];
    const auto* region=layer->regions().front();
    auto* upper=layer->upper_layer->regions().front();
    auto* lower=layer->lower_layer->regions().front();
    ExPolygons original;
    for (const auto& surface:region->slices.surfaces) {
        original.push_back(surface.expolygon);
        original.back().translate(object->instances()[0].shift);
    }
    const auto bounds=get_extents(original);
    const double x=.5*double(bounds.min.x()+bounds.max.x())*1e-6;
    const double y=.5*double(bounds.min.y()+bounds.max.y())*1e-6;
    if (mode==1 || mode==2)
        for (auto& surface:upper->slices.surfaces) surface.expolygon.translate(Point(mode==1 ? 100 : 10000000,0));
    if (mode==3) upper->slices.surfaces.clear();
    if (mode==4) lower->slices.surfaces.clear();
    auto registry=std::make_shared<DepositionProvenance>();
    const DepositionSource identity{object->get_id(),0,std::uint64_t(layer->id()),std::uint64_t(region->region().print_region_id()),0,
        DepositionSource::Height::Transformed,DepositionSource::Assignment::NominalRole};
    const auto id=registry->add(identity);
    registry->add(identity); // Repeated source scopes must share reference storage.
    registry->seal();
    std::size_t reference_points=0;
    for (const auto* reference:std::array<const LayerRegion*,3>{lower,region,upper})
        for (const auto& surface:reference->slices.surfaces) {
            reference_points+=surface.expolygon.contour.points.size();
            for (const auto& hole:surface.expolygon.holes) reference_points+=hole.points.size();
        }
    SourceOwnershipLimits capture_limits;
    capture_limits.max_boundary_points=reference_points;
    std::size_t capture_checks=0;
    auto snapshot=std::make_shared<const NativeOwnershipSnapshot>(capture_native_ownership(print,registry,Point(0,0),0.,capture_limits,
        [&] { ++capture_checks; }));
    // Transformed capture visits the identity index and copies each selected
    // reference once. It does not run or charge all-layer planar validation.
    capture_limits.max_work=capture_checks+reference_points;
    CHECK_NOTHROW(capture_native_ownership(print,registry,Point(0,0),0.,capture_limits));
    auto too_small=capture_limits;
    too_small.max_boundary_points=reference_points-1;
    CHECK_THROWS_AS(capture_native_ownership(print,registry,Point(0,0),0.,too_small),ResourceLimit);
    GCodeProcessorResult result;
    result.deposition_provenance=registry;
    result.native_ownership=snapshot;
    DepositionCommandSection command{{float(x-.1),float(y),.5f},{float(x+.1),float(y),.45f},.4f,.3f,.3f,0};
    command.transform=DepositionCommandSection::Transform::BrickCourse;
    auto sections=std::make_shared<DepositionSections>(registry->generation());
    const auto command_id=sections->add(command);
    auto ordinary=command;
    ordinary.transform=DepositionCommandSection::Transform::Ordinary;
    const auto ordinary_id=sections->add(ordinary);
    auto arc=command;
    arc.end[2]=arc.start[2];
    arc.arc_direction=1;
    arc.arc_center={.5*(arc.start[0]+arc.end[0]),arc.start[1]};
    const auto arc_id=sections->add(arc);
    auto taper=command;
    taper.end_height=.2f;
    const auto taper_id=sections->add(taper);
    auto anchored=command;
    anchored.nominal_lower_z=std::array<double,2>{anchored.start[2]-anchored.start_height,anchored.end[2]-anchored.end_height};
    const auto anchored_id=sections->add(anchored);
    sections->seal();
    result.deposition_sections=sections;
    GCodeProcessorResult::MoveVertex first{},last{};
    first.type=EMoveType::Travel;
    first.position={float(x-.1),float(y),.5f};
    last.type=EMoveType::Extrude;
    last.position={float(x+.1),float(y),.45f};
    last.extrusion_role=erStaggeredPerimeter;
    last.deposition_source_id=id;
    last.deposition_command_section_id=command_id;
    last.deposition_design_heights=std::array<double,2>{.3f,.3f};
    last.deposition_design_width=.4f;
    result.moves={first,last};
    SourceBead bead{{first.position.cast<double>()*.001,last.position.cast<double>()*.001,double(.4f)*.001,double(.3f)*.001,0},id,0,1};
    bead.geometry.end_height_m=double(.3f)*.001;
    const auto cache=snapshot->make_containment_cache();
    const auto resolve=[&](double z) { return snapshot->section(registry,{bead},z,{0},{},{},&result,nullptr,cache.get()); };
    if (mode>=3) {
        CHECK_THROWS_AS(resolve(.0003),std::invalid_argument);
        CHECK_THROWS_AS(resolve(.000425),std::invalid_argument);
        return;
    }
    const auto current=resolve(.0003);
    REQUIRE(current.regions.size()==1);
    CHECK(current.regions.front().footprint.native()==original);
    if (mode==2) {
        // Same source/XY/allowance, different actual-Z reference: a warm
        // original-layer containment certificate cannot authorize this mask.
        CHECK_THROWS_AS(resolve(.000425),std::invalid_argument);
        return;
    }
    ExPolygons expected_upper;
    for (const auto& surface:upper->slices.surfaces) {
        expected_upper.push_back(surface.expolygon);
        expected_upper.back().translate(object->instances()[0].shift);
    }
    const auto above=resolve(.000425);
    REQUIRE(above.regions.size()==1);
    CHECK(above.regions.front().footprint.native()==expected_upper);
    CHECK(resolve(.000175).regions.size()==1);
    CHECK_THROWS_AS(resolve(.00065),std::invalid_argument);
    for (const auto unsupported:{arc_id,taper_id,anchored_id}) {
        result.moves[1].deposition_command_section_id=unsupported;
        CHECK_THROWS_AS(resolve(.000425),std::invalid_argument);
    }
    result.moves[1]=last;
    result.moves[1].extrusion_role=erPerimeter;
    CHECK(resolve(.000425).regions.size()==1);
    result.moves[1]=last;
    result.moves[1].deposition_command_section_id=ordinary_id;
    CHECK_THROWS_AS(resolve(.000425),std::invalid_argument);
    result.moves[1].deposition_command_section_id=0;
    CHECK_THROWS_AS(resolve(.000425),std::invalid_argument);
    result.moves[1]=last;
    result.moves[1].extrusion_role=erInternalInfill;
    CHECK_THROWS_AS(resolve(.000425),std::invalid_argument);
    result.moves[1]=last;
    result.moves[1].extruder_id=1;
    CHECK_THROWS_AS(resolve(.000425),std::invalid_argument);
    result.moves[1]=last;
    result.moves[1].position.y()+=.1f;
    CHECK_THROWS_AS(resolve(.000425),std::invalid_argument);
    result.moves[1]=last;
    result.moves[1].deposition_design_width=.5f;
    CHECK_THROWS_AS(resolve(.000425),std::invalid_argument);
    result.moves[1]=last;
    auto stale=std::make_shared<DepositionSections>(registry->generation()+1);
    stale->add(command);
    stale->seal();
    result.deposition_sections=stale;
    CHECK_THROWS_AS(resolve(.000425),std::invalid_argument);
    result.deposition_sections=sections;
    SourceOwnershipLimits limited;
    limited.max_work=1;
    CHECK_THROWS_AS(snapshot->section(registry,{bead},.000425,{0},limited,{},&result,nullptr,cache.get()),ResourceLimit);
    CHECK_THROWS_AS(snapshot->section(registry,{bead},.000425,{0},{},[] { throw std::runtime_error("cancel brick"); },
        &result,nullptr,cache.get()),std::runtime_error);
}

TEST_CASE("Source ownership partitions touching material envelopes without filling native masks", "[ThermoelasticSourceOwnership]")
{
    const auto config = Test::multifilament_config(2, {{"layer_height", "0.2"},
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
    print.apply(model, config);
    print.process();
    const auto* printed = print.get_object(0);
    const auto* layer   = printed->layers().front();
    auto registry       = std::make_shared<DepositionProvenance>();
    std::vector<BoundingBox> boxes(2);
    std::vector<std::uint32_t> ids(2);
    for (const auto* region : layer->regions()) {
        const auto filament = region->region().extruder(frPerimeter) - 1;
        REQUIRE(filament < 2);
        ids[filament] = registry->add({printed->get_id(), 0, std::uint64_t(layer->id()), std::uint64_t(region->region().print_region_id()),
                                       filament, DepositionSource::Height::NativePlanar});
        ExPolygons polygons;
        for (const auto& surface : region->slices.surfaces) {
            polygons.push_back(surface.expolygon);
            polygons.back().translate(printed->instances()[0].shift);
        }
        boxes[filament] = get_extents(polygons);
    }
    registry->seal();
    const auto snapshot = capture_native_ownership(print, registry, Point(0, 0), 0.);
    REQUIRE(boxes[0].max.x() == boxes[1].min.x());
    const double x     = double(boxes[0].max.x()) * 1e-9;
    const double y     = .5 * double(boxes[0].min.y() + boxes[0].max.y()) * 1e-9;
    const double top   = layer->print_z * .001;
    const double plane = top - layer->height * .0005;
    const std::vector<SourceBead>
        sources{{{{x - .0001, y - .0005, top}, {x - .0001, y + .0005, top}, .0004, layer->height * .001, 0}, ids[0], 0},
                {{{x + .0001, y - .0005, top}, {x + .0001, y + .0005, top}, .0004, layer->height * .001, 1}, ids[1], 1}};
    const auto ownership = snapshot.section(registry, sources, plane, {0, 1});
    std::vector<DepositedSegment> beads{sources[0].geometry, sources[1].geometry};
    CHECK_THROWS_AS(section_deposition(beads, plane, 0.), LayerMaterialOverlap);
    const auto resolved = section_deposition(beads, plane, 0., {}, {}, &ownership);
    REQUIRE(resolved.regions.size() == 2);
    CHECK(resolved.resolved_overlap_area_m2 > 0.);
    CHECK(intersection_ex(resolved.regions[0].footprint.native(), resolved.regions[1].footprint.native()).empty());
    beads[1].material   = 0;
    const auto original = section_deposition(beads, plane, 0.);
    CHECK_THAT(resolved.area_m2, Catch::Matchers::WithinAbs(original.area_m2, 1e-12));
    ExPolygons combined = resolved.regions[0].footprint.native();
    combined.insert(combined.end(), resolved.regions[1].footprint.native().begin(), resolved.regions[1].footprint.native().end());
    CHECK(diff_ex(combined, original.regions.front().footprint.native()).empty());
    CHECK(diff_ex(original.regions.front().footprint.native(), combined).empty());
}

TEST_CASE("Native source snapshots distinguish placements of repeated instances", "[ThermoelasticSourceOwnership]")
{
    const auto config = Test::multifilament_config(1, {{"layer_height", "0.2"},
                                                       {"initial_layer_print_height", "0.2"},
                                                       {"thick_bridges", "0"},
                                                       {"thick_internal_bridges", "0"},
                                                       {"brim_type", "no_brim"},
                                                       {"skirt_loops", "0"}});
    Print print;
    Model model;
    Test::init_print({Test::cube(1.)}, print, model, config);
    model.objects.front()->add_instance()->set_offset(Eigen::Vector3d(20., 0., 0.));
    print.apply(model, config);
    print.process();
    const auto* object = print.get_object(0);
    REQUIRE(object->instances().size() == 2);
    const auto* layer  = object->layers().front();
    const auto* region = layer->regions().front();
    auto registry      = std::make_shared<DepositionProvenance>();
    for (std::uint64_t instance = 0; instance < 2; ++instance)
        registry->add({object->get_id(), instance, std::uint64_t(layer->id()), std::uint64_t(region->region().print_region_id()), 0,
                       DepositionSource::Height::NativePlanar});
    registry->seal();
    const auto snapshot = capture_native_ownership(print, registry, Point(0, 0), 0.);
    const double top    = layer->print_z * .001;
    const double plane  = top - layer->height * .0005;
    std::vector<ExPolygons> placed;
    for (std::uint32_t instance = 0; instance < 2; ++instance) {
        ExPolygons expected;
        for (const auto& surface : region->slices.surfaces) {
            expected.push_back(surface.expolygon);
            expected.back().translate(object->instances()[instance].shift);
        }
        const auto box = get_extents(expected);
        const double x = .5 * double(box.min.x() + box.max.x()) * 1e-9;
        const double y = .5 * double(box.min.y() + box.max.y()) * 1e-9;
        SourceBead bead{{{x, y, top}, {x + .0001, y, top}, .0004, layer->height * .001, 0}, instance + 1, 0};
        const auto ownership = snapshot.section(registry, {bead}, plane, {0});
        REQUIRE(ownership.regions.size() == 1);
        CHECK(diff_ex(expected, ownership.regions.front().footprint.native()).empty());
        CHECK(diff_ex(ownership.regions.front().footprint.native(), expected).empty());
        placed.push_back(ownership.regions.front().footprint.native());
    }
    CHECK(intersection_ex(placed[0], placed[1]).empty());
}

TEST_CASE("Native source snapshots retain immutable masks and exact registry identity", "[ThermoelasticSourceOwnership]")
{
    const auto config = Test::multifilament_config(1, {{"layer_height", "0.2"},
                                                       {"initial_layer_print_height", "0.2"},
                                                       {"thick_bridges", "0"},
                                                       {"thick_internal_bridges", "0"},
                                                       {"brim_type", "no_brim"},
                                                       {"skirt_loops", "0"}});
    Print print;
    Model model;
    Test::init_print({Test::cube(1.)}, print, model, config);
    print.process();
    const auto* object = print.get_object(0);
    auto* layer        = print.get_object(0)->layers().front();
    auto* region       = layer->regions().front();
    DepositionSource source{object->get_id(),
                            0,
                            std::uint64_t(layer->id()),
                            std::uint64_t(region->region().print_region_id()),
                            0,
                            DepositionSource::Height::NativePlanar};
    auto registry             = std::make_shared<DepositionProvenance>();
    const auto id             = registry->add(source);
    auto transformed_source   = source;
    transformed_source.height = DepositionSource::Height::Transformed;
    const auto transformed_id = registry->add(transformed_source);
    registry->seal();
    const Point offset(100000000, 200000000);
    const auto snapshot = capture_native_ownership(print, registry, offset, .001);
    CHECK(snapshot.source_status(id) == NativeOwnershipSnapshot::SourceStatus::NativePlanar);
    CHECK(snapshot.source_status(transformed_id) == NativeOwnershipSnapshot::SourceStatus::Unsupported);
    const double top = .001 * layer->print_z + .001;
    SourceBead bead{{{.1, .2, top}, {.1001, .2, top}, .0004, layer->height * .001, 7}, id, 0};
    const double plane   = top - .5 * bead.geometry.height_m;
    const auto ownership = snapshot.section(registry, {bead}, plane, {7});
    REQUIRE(ownership.regions.size() == 1);
    CHECK(ownership.regions.front().material == 7);
    ExPolygons expected;
    for (auto& surface : region->slices.surfaces) {
        expected.push_back(surface.expolygon);
        expected.back().translate(object->instances()[0].shift + offset);
        surface.expolygon.translate(Point(5000000, 0));
    }
    const auto unchanged = snapshot.section(registry, {bead}, plane, {7});
    CHECK(diff_ex(expected, unchanged.regions.front().footprint.native()).empty());
    CHECK(diff_ex(unchanged.regions.front().footprint.native(), expected).empty());
    auto other_registry = std::make_shared<DepositionProvenance>();
    other_registry->add(source);
    other_registry->seal();
    CHECK_THROWS_AS(snapshot.section(other_registry, {bead}, plane, {7}), std::invalid_argument);
    auto invalid      = bead;
    invalid.source_id = transformed_id;
    CHECK_THROWS_AS(snapshot.section(registry, {invalid}, plane, {7}), std::invalid_argument);
    invalid.source_id = 0;
    CHECK_THROWS_AS(snapshot.section(registry, {invalid}, plane, {7}), std::invalid_argument);
    invalid = bead;
    invalid.geometry.end_m.z() += .00001;
    CHECK_THROWS_AS(snapshot.section(registry, {invalid}, plane, {7}), std::invalid_argument);
    invalid                  = bead;
    invalid.logical_filament = 1;
    CHECK_THROWS_AS(snapshot.section(registry, {invalid}, plane, {7}), std::invalid_argument);
    SourceOwnershipLimits limits;
    limits.centerline_boundary_allowance_nm = 8;
    CHECK_THROWS_AS(snapshot.section(registry, {bead}, plane, {7}, limits), std::invalid_argument);
    limits.centerline_boundary_allowance_nm = 4;
    CHECK_THROWS_AS(snapshot.section(registry, {invalid}, plane, {7}, limits), std::invalid_argument);
    const auto sensitivity = snapshot.section(registry, {bead}, plane, {7}, limits);
    CHECK(diff_ex(sensitivity.regions.front().footprint.native(),expected).empty());
    CHECK(diff_ex(expected,sensitivity.regions.front().footprint.native()).empty());
    limits = {};
    limits.max_boundary_points = 3;
    CHECK_THROWS_AS(snapshot.section(registry, {bead}, plane, {7}, limits), ResourceLimit);
    CHECK_THROWS_AS(snapshot.section(registry, {bead}, plane, {7}, {}, [] { throw std::runtime_error("cancel"); }), std::runtime_error);
}

TEST_CASE("Native source capture rejects invalid identity transformed sources and purge ambiguity", "[ThermoelasticSourceOwnership]")
{
    const auto config = Test::multifilament_config(1, {{"layer_height", "0.2"},
                                                       {"initial_layer_print_height", "0.2"},
                                                       {"thick_bridges", "0"},
                                                       {"thick_internal_bridges", "0"},
                                                       {"brim_type", "no_brim"},
                                                       {"skirt_loops", "0"}});
    Print print;
    Model model;
    Test::init_print({Test::cube(1.)}, print, model, config);
    print.process();
    const auto* object = print.get_object(0);
    const auto* layer  = object->layers().front();
    const DepositionSource good{object->get_id(),
                                0,
                                std::uint64_t(layer->id()),
                                std::uint64_t(layer->regions().front()->region().print_region_id()),
                                0,
                                DepositionSource::Height::NativePlanar};
    const auto make_registry = [](DepositionSource source) {
        auto registry = std::make_shared<DepositionProvenance>();
        registry->add(source);
        registry->seal();
        return registry;
    };
    for (int field = 0; field < 4; ++field) {
        auto bad = good;
        if (field == 0)
            bad.object += 100;
        if (field == 1)
            bad.instance = 100;
        if (field == 2)
            bad.layer = 1000;
        if (field == 3)
            bad.region = 1000;
        CHECK_THROWS_AS(capture_native_ownership(print, make_registry(bad), Point(0, 0), 0.), std::invalid_argument);
    }
    SourceOwnershipLimits limits;
    limits.max_work = 1;
    CHECK_THROWS_AS(capture_native_ownership(print, make_registry(good), Point(0, 0), 0., limits), ResourceLimit);
    CHECK_THROWS_AS(capture_native_ownership(print, make_registry(good), Point(0, 0), 0., {}, [] { throw std::runtime_error("cancel"); }),
                    std::runtime_error);
    auto purge = config;
    purge.set("flush_into_infill", true);
    print.apply(model, purge);
    print.process();
    const auto* purge_object  = print.get_object(0);
    const auto* purge_layer   = purge_object->layers().front();
    auto purge_source         = good;
    purge_source.object       = purge_object->get_id();
    purge_source.layer        = purge_layer->id();
    purge_source.region       = purge_layer->regions().front()->region().print_region_id();
    const auto purge_registry = make_registry(purge_source);
    const auto unsupported    = capture_native_ownership(print, purge_registry, Point(0, 0), 0.);
    CHECK(unsupported.source_status(1) == NativeOwnershipSnapshot::SourceStatus::Unsupported);
    SourceBead bead{{{0, 0, .0002}, {.001, 0, .0002}, .0004, .0002, 0}, 1, 0};
    CHECK_THROWS_AS(unsupported.section(purge_registry, {bead}, .0001, {0}), std::invalid_argument);
}
