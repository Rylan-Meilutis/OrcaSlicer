#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/SectionGrowthStream.hpp"
#include "libslic3r/SectionGrowth.hpp"
#include "libslic3r/SectionThermalAssembly.hpp"
#include "libslic3r/ThermoelasticLocalSections.hpp"
#include "libslic3r/ThermoelasticLayers.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticSectionTopology.hpp"
#include "libslic3r/ThermoelasticExactSections.hpp"
#include <cstddef>
#include <Eigen/Core>
#include <stdexcept>
#include <string>
#include <vector>
#include <array>

using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;

TEST_CASE("Exact reference activation composes each source prefix from unchanged individual cuts", "[SectionGrowthStream]")
{
    const std::vector<DepositedSegment> paths{
        {{0.,0.,.0002},{.001,double(.333333f)*.001,.0002},.0004,.0002,0},
        {{.0001,.0002,.0002},{.001,double(-.222222f)*.001,.0002},.0004,.0002,0}};
    SectionGrowthStreamSettings settings;
    CHECK_FALSE(settings.exact_reference_sources);
    CHECK_FALSE(settings.geometry.geometry.exact_per_command_composition);
    settings.origin_m={-.001,-.001};
    settings.tile_size_m={.004,.004};
    settings.maximum_height_m=.0002;
    settings.geometry.geometry.ownership_boundary_grid_units=4;
    settings.geometry.geometry.retain_exact_interface_sections=true;
    settings.geometry.geometry.max_slab_group_visits=4000000;
    const auto native=section_deposition(paths,.0001,0.,settings.geometry.geometry);
    REQUIRE(native.regions.size()==1);
    CHECK_FALSE(native.regions.front().footprint.is_exact());
    ExactSection complete;
    for (const auto& bead:paths) {
        const auto cut=section_deposition({bead},.0001,0.,settings.geometry.geometry);
        REQUIRE(cut.regions.size()==1);
        complete=boolean_exact_sections(complete,make_exact_section(cut.regions.front().footprint.native()).section,
            ExactSectionOperation::Union).section;
    }
    CHECK_FALSE(boolean_exact_sections(complete,make_exact_section(native.regions.front().footprint.native()).section,
        ExactSectionOperation::SymmetricDifference).section.empty());
    // This single-material geometry reference has no ownership certificate.
    // Its explicit composition matches the first-policy sampler's operation.
    settings.geometry.geometry.exact_per_command_composition=true;
    settings.exact_reference_sources=true;
    for (std::size_t count=1;count<=paths.size();++count) {
        const std::vector<DepositedSegment> prefix(paths.begin(),paths.begin()+count);
        const auto expected=section_deposition(prefix,.0001,0.,settings.geometry.geometry);
        REQUIRE(expected.regions.size()==1);
        CHECK_FALSE(expected.deposited_partition);
        CHECK_FALSE(expected.occupied_partition);
        ExactSection activated;
        const auto proof=visit_reference_activation(paths,prefix,settings,[&](const SectionPiece& piece) {
            REQUIRE(piece.slab.footprint.is_exact());
            activated=boolean_exact_sections(activated,piece.slab.footprint.exact(),ExactSectionOperation::Union).section;
        });
        CHECK(boolean_exact_sections(activated,expected.regions.front().footprint.exact(),
            ExactSectionOperation::SymmetricDifference).section.empty());
        CHECK(validate_exact_section_topology(activated)>0);
        CHECK_THAT(proof.excluded_nominal_volume_m3,WithinAbs(0.,1e-25));
    }
    auto legacy=settings.geometry.geometry;
    legacy.exact_per_command_composition=false;
    const auto unchanged=section_deposition(paths,.0001,0.,legacy);
    CHECK(unchanged.regions.front().footprint.native()==native.regions.front().footprint.native());
    auto mixed=paths;
    mixed.back().material=1;
    CHECK_THROWS_AS(section_deposition(mixed,.0001,0.,settings.geometry.geometry),std::invalid_argument);
    auto invalid=settings.geometry.geometry;
    invalid.retain_exact_interface_sections=false;
    CHECK_THROWS_AS(section_deposition(paths,.0001,0.,invalid),std::invalid_argument);
    auto bounded=settings.geometry;
    bounded.max_ownership_work=1;
    CHECK_THROWS_AS(visit_local_deposition_sections(paths,settings.origin_m,settings.tile_size_m,0.,.0002,
        [](const LocalSectionTile&) {},bounded),ResourceLimit);
    CHECK_THROWS_AS(visit_reference_activation(paths,paths,settings,[](const SectionPiece&) {},{},
        [] { throw std::runtime_error("cancel exact source composition"); }),std::runtime_error);
}

TEST_CASE("Complete deposited reference reuse refuses unproven packets and poisons traversal", "[SectionGrowthStream]")
{
    const std::vector<DepositedSegment> paths{{{0.,0.,.0002},{.001,0.,.0002},.0004,.0002,0}};
    SectionGrowthStreamSettings settings;
    settings.origin_m={-.001,-.001}; settings.tile_size_m={.004,.004};
    settings.maximum_height_m=.0002;
    settings.geometry.geometry.retain_exact_interface_sections=true;
    settings.geometry.geometry.ownership_boundary_grid_units=4;
    settings.exact_reference_sources=true;
    CHECK_FALSE(settings.reuse_complete_deposited_reference);
    CHECK_NOTHROW(visit_reference_activation(paths,paths,settings,[](const SectionPiece&) {}));
    settings.reuse_complete_deposited_reference=true;
    ReferenceActivationCursor cursor(paths,paths,settings,[](const SectionPiece&) {});
    CHECK_THROWS_WITH(cursor.advance(),"Complete reference activation requires deposited partition evidence.");
    CHECK_THROWS_AS(cursor.advance(),std::logic_error);
    CHECK_THROWS_AS(cursor.finish(),std::logic_error);
    CHECK_THROWS_WITH(visit_reference_activation(paths,paths,settings,[](const SectionPiece&) {},{},
        [] { throw std::runtime_error("cancel complete reference"); }),"cancel complete reference");
    settings.max_work=1;
    CHECK_THROWS_AS(visit_reference_activation(paths,paths,settings,[](const SectionPiece&) {}),ResourceLimit);
}

TEST_CASE("Reference activation resource failures identify the consumer stage and poison traversal", "[SectionGrowthStream]")
{
    const unsigned failure=GENERATE(0u,1u,2u);
    const std::vector<DepositedSegment> paths{
        {{.001,.001,.0002},{.003,.001,.0002},.0004,.0002,0,BeadSection::RectangularEnvelope},
        {{.0011,.001,.0004},{.0031,.001,.0004},.0004,.0002,0,BeadSection::RectangularEnvelope}};
    SectionGrowthStreamSettings settings;
    settings.tile_size_m={.01,.01};
    settings.maximum_height_m=.0002;
    if (failure==0) settings.max_output_pieces=1;
    LocalSectionPacketVisitor packets;
    packets.begin_tile=[](const LocalSectionTile&) {};
    packets.band=[&](const LocalSectionTile&) {
        if (failure==2) throw ResourceLimit("consumer resource marker");
    };
    packets.end_tile=[] {};
    ReferenceActivationCursor cursor(paths,paths,settings,[&](const SectionPiece&) {
        if (failure==1) throw ResourceLimit("consumer resource marker");
    },{},{},packets);
    try {
        while (cursor.advance()) {}
        FAIL("Expected staged activation resource failure");
    } catch (const ResourceLimit& error) {
        const std::string message=error.what();
        const std::array<const char*,3> stages{"stage=clipping","stage=active-topology","stage=active-packet-consumer"};
        CHECK(message.find(stages[failure])!=std::string::npos);
        CHECK(message.find("tile=(0,0)")!=std::string::npos);
        CHECK(message.find("band_nm=[")!=std::string::npos);
        CHECK(message.find("activation_work=")!=std::string::npos);
        CHECK(message.find("activation_limit=")!=std::string::npos);
        CHECK(message.find(failure==0 ? "Reference activation cumulative output budget exceeded." :
            "consumer resource marker")==0);
    }
    CHECK_THROWS_AS(cursor.next_tile_index(),std::logic_error);
    CHECK_THROWS_AS(cursor.advance(),std::logic_error);
}

TEST_CASE("Reference activation reports which sampler exhausted its work budget", "[SectionGrowthStream]")
{
    const std::vector<DepositedSegment> paths{
        {{.001,.001,.0002},{.003,.001,.0002},.0004,.0002,0,BeadSection::RectangularEnvelope}};
    SectionGrowthStreamSettings settings;
    settings.tile_size_m={.01,.01};
    settings.maximum_height_m=.0002;
    settings.exact_reference_sources=true;
    settings.geometry.geometry.retain_exact_interface_sections=true;
    settings.geometry.geometry.ownership_boundary_grid_units=4;
    settings.geometry.max_ownership_work=1;
    ReferenceActivationCursor cursor(paths,paths,settings,[](const SectionPiece&) {});
    try {
        while (cursor.advance()) {}
        FAIL("Expected the exact source sampler work budget to reject");
    } catch (const ResourceLimit& error) {
        const std::string message=error.what();
        CHECK(message.find("Local section ownership total work budget exceeded.")==0);
        CHECK(message.find("Reference activation sampler=source-0.")!=std::string::npos);
    }
    CHECK_THROWS_AS(cursor.advance(),std::logic_error);
}

TEST_CASE("Reference activation bounds live geometry by bands instead of tile height", "[SectionGrowthStream]")
{
    std::vector<DepositedSegment> paths;
    for (unsigned layer=0;layer<120;++layer) {
        const double z=(layer+1)*.0002, x=.001+(layer%2)*.0001;
        paths.push_back({{x,.001,z},{x+.002,.001,z},.0004,.0002,0,BeadSection::RectangularEnvelope});
    }
    SectionGrowthStreamSettings settings;
    settings.tile_size_m={.01,.01};
    settings.maximum_height_m=.0002;
    settings.max_tile_points=128;
    std::size_t pieces=0;
    double volume=0.;
    const auto proof=visit_reference_activation(paths,paths,settings,[&](const SectionPiece& piece) {
        ++pieces;
        for (const auto& p : piece.slab.footprint.native())
            volume+=p.area()*1e-18*(piece.slab.top_m-piece.slab.bottom_m);
    });
    CHECK(pieces==120);
    CHECK(proof.active.components()==1);
    CHECK_THAT(proof.excluded_nominal_volume_m3,WithinAbs(0.,1e-18));
    const auto retained=local_deposition_sections(paths,settings.origin_m,settings.tile_size_m,0.,.0002);
    CHECK(retained.boundary_points>settings.max_tile_points);
    CHECK_THAT(volume,WithinAbs(retained.volume_m3,1e-16));
}

TEST_CASE("Streamed fixed reference activation matches retained geometry and reports excluded nominal material", "[SectionGrowthStream]")
{
    const std::vector<DepositedSegment> reference{
        {{0.,.0005,.001},{.003,.0005,.001},.001,.001,0,BeadSection::RectangularEnvelope},
        {{0.,.0015,.001},{.003,.0015,.001},.001,.001,1,BeadSection::RectangularEnvelope}};
    auto paths = reference;
    for (auto& path : paths) path.width_m = .0014;
    paths.push_back({{.01,.0005,.001},{.011,.0005,.001},.001,.001,0,BeadSection::RectangularEnvelope});
    SectionGrowthStreamSettings settings;
    settings.tile_size_m = {.001,.001};
    settings.maximum_height_m = .0002;
    const auto fixed = local_deposition_sections(reference,settings.origin_m,settings.tile_size_m,0.,.0002);
    const auto topology = qualify_local_section_topology(fixed.tiles,0.);
    const auto first = local_deposition_sections({paths[0],paths[2]},settings.origin_m,settings.tile_size_m,0.,.0002);
    const auto second = local_deposition_sections({paths[1]},settings.origin_m,settings.tile_size_m,0.,.0002);
    const auto expected = activate_reference_sections(topology,
        {qualify_local_section_topology(first.tiles,0.),qualify_local_section_topology(second.tiles,0.)});
    std::array<double,2> volumes{};
    const auto streamed = visit_reference_activation(reference,paths,settings,[&](const SectionPiece& p) {
        REQUIRE(p.slab.material < volumes.size());
        for (const auto& polygon : p.slab.footprint.native())
            volumes[p.slab.material] += polygon.area()*1e-18*(p.slab.top_m-p.slab.bottom_m);
    });
    std::array<double,2> expected_volumes{};
    for (const auto& p : expected.geometry.pieces())
        for (const auto& polygon : p.slab.footprint.native())
            expected_volumes[p.slab.material] += polygon.area()*1e-18*(p.slab.top_m-p.slab.bottom_m);
    for (std::size_t i=0;i<volumes.size();++i)
        CHECK_THAT(volumes[i],WithinAbs(expected_volumes[i],1e-17));
    CHECK_THAT(streamed.excluded_nominal_volume_m3,WithinAbs(expected.excluded_nominal_volume_m3,1e-17));
    CHECK(streamed.excluded_nominal_volume_m3 > 0.);
    CHECK(streamed.active.components() == expected.geometry.components());
    ThermalSectionStreamSettings thermal;
    thermal.spacing_m = {.001,.001,.001};
    thermal.tile_size_m = settings.tile_size_m;
    thermal.maximum_height_m = settings.maximum_height_m;
    const std::vector<SectionThermalMaterial> properties{
        {Eigen::Matrix3d::Identity(),2e6,400.},{Eigen::Matrix3d::Identity(),3e6,500.}};
    const auto assembly = assemble_reference_thermal_volume(reference,paths,properties,thermal);
    const auto retained = assemble_section_thermal_volume(expected.geometry,properties,thermal.origin_m,thermal.spacing_m);
    CHECK(assembly.assembly.exposure_complete);
    CHECK_THAT(assembly.excluded_nominal_volume_m3,WithinAbs(expected.excluded_nominal_volume_m3,1e-17));
    CHECK_THAT(assembly.assembly.occupied_volume_m3,WithinAbs(retained.occupied_volume_m3,1e-17));
    const auto heat = [](const SectionThermalAssembly& field) {
        double value = 0.;
        for (const auto& cell : field.initial.deposited_cells)
            value += cell.volume_m3*cell.volumetric_heat_capacity_j_m3_k*cell.temperature_k;
        return value;
    };
    const auto area = [](const SectionThermalAssembly& field) {
        double value = 0.;
        for (const auto& face : field.horizontal_faces) value += face.area_moments[0];
        for (const auto& face : field.vertical_faces) value += face.area_moments[0];
        return value;
    };
    CHECK_THAT(heat(assembly.assembly),WithinAbs(heat(retained),1e-8));
    CHECK_THAT(area(assembly.assembly),WithinAbs(area(retained),1e-13));
    bool future_material = false;
    visit_reference_activation(reference,{paths[0]},settings,[&](const SectionPiece& p) {
        future_material |= p.slab.material != 0;
    });
    CHECK_FALSE(future_material);
    double added_volume = 0.;
    const auto growth = visit_reference_section_growth(reference,{paths[0]},paths,settings,
        {[](const SectionPiece&) {},[](const SectionPiece&) {},[&](const SectionPiece& p) {
            CHECK(p.slab.material == 1);
            for (const auto& polygon : p.slab.footprint.native())
                added_volume += polygon.area()*1e-18*(p.slab.top_m-p.slab.bottom_m);
        }});
    CHECK_THAT(added_volume,WithinAbs(expected_volumes[1],1e-17));
    CHECK_THAT(growth.next_excluded_nominal_volume_m3,WithinAbs(streamed.excluded_nominal_volume_m3,1e-17));
    REQUIRE(growth.previous_components.size() == 1);
    CHECK(growth.previous_components[0] == 0);
    ReferenceActivationCursor* observed = nullptr;
    std::vector<std::size_t> representatives;
    LocalSectionPacketVisitor packets;
    packets.begin_tile = [](const LocalSectionTile&) {};
    packets.band = [&](const LocalSectionTile&) { representatives = observed->active_component_mapping(); };
    packets.end_tile = [] {};
    ReferenceActivationCursor inspected(reference,paths,settings,[](const SectionPiece&) {},{},{},packets);
    observed = &inspected;
    while (inspected.advance()) {}
    const auto inspected_proof = inspected.finish();
    REQUIRE_FALSE(representatives.empty());
    for (std::size_t i = 0; i < representatives.size(); ++i)
        CHECK(inspected_proof.active.component_mapping().at(i) ==
            inspected_proof.active.component_mapping().at(representatives[i]));
    CHECK_THROWS_AS(inspected.active_component_mapping(),std::logic_error);
    ReferenceActivationCursor unfinished(reference,paths,settings,[](const SectionPiece&) {});
    CHECK_THROWS_AS(unfinished.finish(),std::logic_error);
    CHECK_THROWS_AS(unfinished.advance(),std::logic_error);
    ReferenceActivationCursor* active = nullptr;
    ReferenceActivationCursor reentrant(reference,paths,settings,[&](const SectionPiece&) {
        try { active->next_tile_index(); } catch (const std::logic_error&) {}
    });
    active = &reentrant;
    CHECK_THROWS_AS(([&] { while (reentrant.advance()) {} }()),std::logic_error);
    CHECK_THROWS_AS(reentrant.finish(),std::logic_error);
    auto limited = settings;
    const std::vector<DepositedSegment> source_a{paths[0],paths[2]}, source_b{paths[1]};
    const LocalSectionCursor index_reference(reference,settings.origin_m,settings.tile_size_m,0.,.0002);
    const LocalSectionCursor index_a(source_a,settings.origin_m,settings.tile_size_m,0.,.0002);
    const LocalSectionCursor index_b(source_b,settings.origin_m,settings.tile_size_m,0.,.0002);
    const auto ref_usage=index_reference.index_usage(), a_usage=index_a.index_usage(), b_usage=index_b.index_usage();
    limited.geometry.max_tiles=ref_usage.tiles+a_usage.tiles+b_usage.tiles;
    limited.geometry.max_tile_references=ref_usage.references+a_usage.references+b_usage.references;
    CHECK(a_usage.tiles>limited.geometry.max_tiles/3);
    CHECK_NOTHROW(visit_reference_activation(reference,paths,limited,[](const SectionPiece&) {}));
    --limited.geometry.max_tiles;
    CHECK_THROWS_AS(visit_reference_activation(reference,paths,limited,[](const SectionPiece&) {}),ResourceLimit);
    limited = settings;
    limited.max_tile_points = 3;
    CHECK_THROWS_AS(visit_reference_activation(reference,paths,limited,[](const SectionPiece&) {}),ResourceLimit);
    limited = settings;
    limited.max_work = 1;
    CHECK_THROWS_AS(visit_reference_activation(reference,paths,limited,[](const SectionPiece&) {}),ResourceLimit);
    CHECK_THROWS_AS(visit_reference_activation(reference,paths,settings,[](const SectionPiece&) {},{},
        [] { throw std::runtime_error("cancel reference"); }),std::runtime_error);
}

TEST_CASE("Paired section streams certify cross tile joins and emit only added material", "[SectionGrowthStream]")
{
    const std::vector<DepositedSegment> previous{{{0., .0005, .001}, {.001, .0005, .001}, .001, .001, 0, BeadSection::RectangularEnvelope},
                                                 {{.002, .0005, .001}, {.003, .0005, .001}, .001, .001, 0, BeadSection::RectangularEnvelope}};
    auto next = previous;
    next.push_back({{.001, .0005, .001}, {.002, .0005, .001}, .001, .001, 0, BeadSection::RectangularEnvelope});
    SectionGrowthStreamSettings settings;
    settings.tile_size_m      = {.001, .001};
    settings.maximum_height_m = .0002;
    double added_volume       = 0.;
    std::vector<std::size_t> labels;
    SectionGrowthStreamVisitor visitor{[](const SectionPiece&) {}, [](const SectionPiece&) {},
                                       [&](const SectionPiece& p) {
                                           for (const auto& polygon : p.slab.footprint.native())
                                               added_volume += polygon.area() * 1e-18 * (p.slab.top_m - p.slab.bottom_m);
                                           labels.push_back(p.component);
                                       }};
    const auto result = visit_local_section_growth(previous, next, settings, visitor);
    REQUIRE(result.previous.components() == 2);
    REQUIRE(result.next.components() == 1);
    REQUIRE(result.previous_components.size() == 2);
    CHECK(result.previous_components[0] == 0);
    CHECK(result.previous_components[1] == 0);
    // Nominal swept paths include their XY end caps. Compare against the
    // independently retained native polygon difference, not length*width.
    const auto old_geometry = local_deposition_sections(previous, settings.origin_m, settings.tile_size_m, 0., .0002);
    const auto new_geometry = local_deposition_sections(next, settings.origin_m, settings.tile_size_m, 0., .0002);
    const auto old_topology = qualify_local_section_topology(old_geometry.tiles, 0.);
    const auto new_topology = qualify_local_section_topology(new_geometry.tiles, 0.);
    const auto retained     = qualify_section_growth(old_topology, new_topology);
    double expected         = 0.;
    for (const auto& p : retained.added_pieces())
        for (const auto& polygon : p.slab.footprint.native())
            expected += polygon.area() * 1e-18 * (p.slab.top_m - p.slab.bottom_m);
    CHECK(expected > 0.);
    CHECK_THAT(added_volume, WithinAbs(expected, 1e-18));
    CHECK_THAT(result.next_geometry.volume_m3 - result.previous_geometry.volume_m3, WithinAbs(added_volume, 1e-18));
    REQUIRE_FALSE(labels.empty());
    for (const auto label : labels)
        CHECK(result.next.component_mapping().at(label) == 0);
    ThermalSectionStreamSettings thermal;
    thermal.spacing_m        = {.001, .001, .001};
    thermal.tile_size_m      = settings.tile_size_m;
    thermal.maximum_height_m = .0002;
    const std::vector<SectionThermalMaterial> materials{{Eigen::Matrix3d::Identity(), 2e6, 400.}};
    const auto assembled = assemble_local_thermal_growth(previous, next, materials, thermal);
    CHECK(assembled.next.exposure_complete);
    CHECK_FALSE(assembled.added.exposure_complete);
    CHECK(assembled.previous_components == result.previous_components);
    CHECK_THAT(assembled.added.occupied_volume_m3, WithinAbs(expected, 1e-18));
    double heat = 0.;
    for (const auto& cell : assembled.added.initial.deposited_cells)
        heat += cell.volume_m3 * cell.volumetric_heat_capacity_j_m3_k * cell.temperature_k;
    CHECK_THAT(heat, WithinAbs(expected * 2e6 * 400., 1e-9));
    const auto unchanged = assemble_local_thermal_growth(next, next, materials, thermal);
    CHECK(unchanged.added.nodes.empty());
    CHECK(unchanged.next.exposure_complete);
    added_volume = 0.;
    visit_local_section_growth(next, next, settings, visitor);
    CHECK_THAT(added_volume, WithinAbs(0., 1e-18));
}

TEST_CASE("Paired section growth rejects missing or relabeled geometry and cumulative resource exhaustion", "[SectionGrowthStream]")
{
    const std::vector<DepositedSegment> previous{{{0., .0005, .001}, {.003, .0005, .001}, .001, .001, 0, BeadSection::RectangularEnvelope}};
    SectionGrowthStreamSettings settings;
    settings.tile_size_m      = {.001, .001};
    settings.maximum_height_m = .0002;
    const SectionGrowthStreamVisitor visitor{[](const SectionPiece&) {}, [](const SectionPiece&) {}, [](const SectionPiece&) {}};
    auto changed             = previous;
    changed.front().material = 1;
    CHECK_THROWS_AS(visit_local_section_growth(previous, changed, settings, visitor), std::invalid_argument);
    changed = previous;
    changed.front().start_m.x() += .001;
    changed.front().end_m.x() += .001;
    CHECK_THROWS_AS(visit_local_section_growth(previous, changed, settings, visitor), std::invalid_argument);
    auto restricted     = settings;
    restricted.max_work = 1;
    CHECK_THROWS_AS(visit_local_section_growth(previous, previous, restricted, visitor), ResourceLimit);
    restricted                 = settings;
    restricted.max_tile_points = 3;
    CHECK_THROWS_AS(visit_local_section_growth(previous, previous, restricted, visitor), ResourceLimit);
    CHECK_THROWS_AS(visit_local_section_growth(previous, previous, settings, visitor, {}, {},
                                               [] { throw std::runtime_error("cancel growth"); }),
                    std::runtime_error);
}
