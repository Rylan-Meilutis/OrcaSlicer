#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/ThermoelasticLayers.hpp"
#include "libslic3r/ThermoelasticLocalSections.hpp"
#include "libslic3r/ThermoelasticLocalTopology.hpp"
#include "libslic3r/ThermoelasticExactSections.hpp"
#include "libslic3r/ThermoelasticSectionFootprint.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticMaterial.hpp"
#include "libslic3r/ThermoelasticSections.hpp"
#include "libslic3r/ThermoelasticSectionStudy.hpp"
#include "libslic3r/ThermoelasticSectionTopology.hpp"
#include "libslic3r/SectionIntegration.hpp"
#include "libslic3r/libslic3r.h"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;

TEST_CASE("Local cursor resource failures identify the sampled interval and poison traversal", "[ThermoelasticLocalSections]")
{
    const bool whole_tile=GENERATE(false,true);
    const std::vector<DepositedSegment> beads{{{.001,.001,.001},{.002,.001,.001},.0001,.001,0}};
    LocalSectionPacketVisitor visitor;
    visitor.begin_tile=[](const LocalSectionTile&) {};
    visitor.band=[](const LocalSectionTile&) {};
    visitor.end_tile=[] {};
    LocalSectionCursor cursor(beads,{0.,0.},{.01,.01},0.,.001,{},
        [](double) -> std::optional<SectionOwnership> { throw ResourceLimit("sampling limit marker"); });
    try {
        if (whole_tile) cursor.advance(visitor); else cursor.advance_packet(visitor);
        FAIL("Expected bounded sampling failure");
    } catch (const ResourceLimit& error) {
        const std::string message=error.what();
        CHECK(message.find("sampling limit marker")==0);
        CHECK(message.find("tile=(0,0)")!=std::string::npos);
        CHECK(message.find("Bounds_nm=[0,0,10000000,10000000]")!=std::string::npos);
        CHECK(message.find("Sample_interval_nm=[0,1000000]")!=std::string::npos);
    }
    CHECK_THROWS_AS(cursor.advance_packet(visitor),std::invalid_argument);
    CHECK_THROWS_AS(cursor.next_tile_index(),std::logic_error);
    LocalSectionCursor cancelled(beads,{0.,0.},{.01,.01},0.,.001,{},
        [](double) -> std::optional<SectionOwnership> { throw std::runtime_error("cancel unchanged"); });
    CHECK_THROWS_WITH(cancelled.advance_packet(visitor),"cancel unchanged");
    CHECK_THROWS_AS(cancelled.advance_packet(visitor),std::invalid_argument);
}

namespace {
LocalSectionTile box_tile(int x, double bottom, double top)
{
    const coord_t left = coord_t(x) * 1000000, right = left + 1000000;
    return {{x, 0}, {left, 0}, {right, 1000000},
        {{bottom, top, 0, {ExPolygon(Polygon{Point(left, 0), Point(right, 0),
            Point(right, 1000000), Point(left, 1000000)})}}}};
}
SectionMaterial elastic_material()
{
    ElasticityMatrix d = ElasticityMatrix::Zero();
    d.diagonal().head<3>().setConstant(1000.);
    d.diagonal().tail<3>().setConstant(500.);
    return {d};
}
} // namespace

TEST_CASE("Exact local topology joins rational tile traces in the shared frontier", "[ThermoelasticLocalSections]")
{
    const auto source=make_exact_section({ExPolygon(Polygon{{0,0},{20,0},{0,1}})}).section;
    const auto left_mask=make_exact_section({ExPolygon(Polygon{{0,0},{10,0},{10,10},{0,10}})}).section;
    const auto right_mask=make_exact_section({ExPolygon(Polygon{{10,0},{20,0},{20,10},{10,10}})}).section;
    const auto left=boolean_exact_sections(source,left_mask,ExactSectionOperation::Intersection).section;
    const auto right=boolean_exact_sections(source,right_mask,ExactSectionOperation::Intersection).section;
    const LocalSectionTile a{{0,0},{0,0},{10,10},{{0.,.001,0,SectionFootprint(left)}}};
    const LocalSectionTile b{{1,0},{10,0},{20,10},{{0.,.001,0,SectionFootprint(right)}}};
    for (const auto order : {LocalTileOrder::Arbitrary,LocalTileOrder::Lexicographic}) {
        StreamTopologyStatistics stats;
        std::vector<SectionPiece> pieces;
        LocalTopologyBuilder builder(0.,{},[&](const SectionPiece& piece) { pieces.push_back(piece); },{},
            StreamTopologyLimits{},&stats,order);
        builder.add_tile(a);
        builder.add_tile(b);
        const auto proof=builder.finish();
        CHECK(proof.components()==1);
        REQUIRE(pieces.size()==2);
        for (const auto& piece : pieces) {
            CHECK(piece.slab.footprint.is_exact());
            CHECK(proof.component_mapping().at(piece.component)==0);
        }
        if (order==LocalTileOrder::Lexicographic) CHECK(stats.face_records==0);
    }
    CHECK(qualify_local_section_topology({b,a},0.).components()==1);
    // One native full-height trace also overlaps the rational half-nm trace.
    auto native_a=a;
    native_a.slabs.front().footprint=ExPolygons{ExPolygon(Polygon{{0,0},{10,0},{10,10},{0,10}})};
    CHECK(qualify_local_section_topology({native_a,b},0.).components()==1);
    CHECK(qualify_local_section_topology({b,native_a},0.).components()==1);
    auto point_b=b;
    point_b.slabs.front().footprint=SectionFootprint(make_exact_section({ExPolygon(Polygon{{10,0},{20,0},{20,1}})}).section);
    CHECK(qualify_local_section_topology({a,point_b},0.).components()==2);
}

TEST_CASE("Exact local packets preserve rational vertical contact without welding points", "[ThermoelasticLocalSections]")
{
    const auto full=make_exact_section({ExPolygon(Polygon{{0,0},{10,0},{10,4}})}).section;
    const auto strip=make_exact_section({ExPolygon(Polygon{{0,0},{1,0},{1,10},{0,10}})}).section;
    const auto small=boolean_exact_sections(full,strip,ExactSectionOperation::Intersection).section;
    const LocalSectionTile header{{0,0},{0,0},{10,10},{}};
    const LocalSectionTile lower{{0,0},{0,0},{10,10},{{0.,.001,0,SectionFootprint(full)}}};
    const LocalSectionTile upper{{0,0},{0,0},{10,10},{{.001,.002,0,SectionFootprint(small)}}};
    std::vector<SectionPiece> pieces;
    LocalTopologyBuilder builder(0.,{},[&](const SectionPiece& piece) { pieces.push_back(piece); });
    builder.begin_tile(header);
    builder.add_band(lower);
    builder.add_band(upper);
    builder.end_tile();
    CHECK(builder.finish().components()==1);
    REQUIRE(pieces.size()==2);
    CHECK(boolean_exact_sections(pieces.back().slab.footprint.exact(),small,
        ExactSectionOperation::SymmetricDifference).section.empty());
    auto floating=upper;
    floating.slabs.front().footprint=SectionFootprint(make_exact_section({ExPolygon(Polygon{{0,0},{0,1},{1,1}})}).section);
    auto both=lower;
    both.slabs.push_back(floating.slabs.front());
    CHECK_THROWS_AS(qualify_local_section_topology({both},0.),std::invalid_argument);
}

TEST_CASE("Exact local topology rejects unsupported geometry and bounds its work", "[ThermoelasticLocalSections]")
{
    auto tile=box_tile(0,0.,.001);
    const auto exact=make_exact_section(tile.slabs.front().footprint.native()).section;
    tile.slabs.front().footprint=SectionFootprint(exact);
    auto uncertain=tile;
    uncertain.slabs.front().boundary_uncertainty_nm=7;
    CHECK_THROWS_AS(qualify_local_section_topology({uncertain},0.),std::invalid_argument);
    auto overlap=tile;
    overlap.slabs.push_back({0.,.001,1,SectionFootprint(exact)});
    CHECK_THROWS_AS(qualify_local_section_topology({overlap},0.),LayerMaterialOverlap);
    auto escaped=tile;
    escaped.maximum[0]=500000;
    CHECK_THROWS_AS(qualify_local_section_topology({escaped},0.),std::invalid_argument);
    StreamTopologyLimits limited;
    limited.max_work=1;
    LocalTopologyBuilder work_limited(0.,{},[](const SectionPiece&) {},{},limited);
    CHECK_THROWS_AS(work_limited.add_tile(tile),ResourceLimit);
    CHECK_THROWS_AS(work_limited.finish(),std::invalid_argument);
    limited={}; limited.max_live_points=32;
    LocalTopologyBuilder live_limited(0.,{},[](const SectionPiece&) {},{},limited);
    CHECK_THROWS_AS(live_limited.add_tile(tile),ResourceLimit);
    limited={}; limited.max_trace_intervals=1;
    LocalTopologyBuilder trace_limited(0.,{},[](const SectionPiece&) {},{},limited);
    CHECK_THROWS_AS(trace_limited.add_tile(tile),ResourceLimit);
    CHECK_THROWS_AS(qualify_local_section_topology({tile},0.,{},
        [] { throw std::runtime_error("cancel exact local topology"); }),std::runtime_error);
}

TEST_CASE("Uncertain tile contacts must survive their coordinate margin", "[ThermoelasticLocalSections]")
{
    auto left = box_tile(0, 0., .001), right = box_tile(1, 0., .001);
    left.slabs.front().boundary_uncertainty_nm = 4;
    right.slabs.front().boundary_uncertainty_nm = 4;
    CHECK(qualify_local_section_topology({left, right}, 0.).components() == 1);
    // Both islands remain bed connected, but their lateral contact is only
    // four nanometres long and is not robust to the stated uncertainty.
    auto& points = right.slabs.front().footprint.native().front().contour.points;
    for (auto& point : points) if (point.y() == 0) point.y() = 999996;
    CHECK_THROWS_AS(qualify_local_section_topology({left, right}, 0.), std::invalid_argument);
}

TEST_CASE("Local topology rejects exact material overlap even inside a representation band", "[ThermoelasticLocalSections]")
{
    auto tile=box_tile(0,0.,.001);
    auto second=tile.slabs.front();
    second.material=1;
    second.boundary_uncertainty_nm=4;
    for (auto& point : second.footprint.native().front().contour.points)
        if (point.x()==0) point.x()=999999;
    tile.slabs.push_back(std::move(second));
    // A one-nanometre strip is still genuine overlap, not a small-area waiver.
    CHECK_THROWS_AS(qualify_local_section_topology({tile},0.),LayerMaterialOverlap);
}

TEST_CASE("Later tile contacts anchor previously separated vertical bands", "[ThermoelasticLocalSections]")
{
    auto split = box_tile(0, 0., .001);
    split.slabs.push_back(box_tile(0, .002, .003).slabs.front());
    auto neighbor = box_tile(1, 0., .001);
    neighbor.slabs.push_back(box_tile(1, .001, .002).slabs.front());
    neighbor.slabs.push_back(box_tile(1, .002, .003).slabs.front());
    CHECK_THROWS_AS(qualify_local_section_topology({split}, 0.), std::invalid_argument);
    const auto forward = qualify_local_section_topology({split, neighbor}, 0.);
    const auto reverse = qualify_local_section_topology({neighbor, split}, 0.);
    CHECK(forward.components() == 1);
    CHECK(reverse.components() == 1);
    CHECK(forward.pieces().size() == 5);
    CHECK(reverse.pieces().size() == forward.pieces().size());
    for (const auto& piece : forward.pieces()) {
        CHECK_FALSE(piece.slab.footprint.native().empty());
        CHECK(piece.component == 0);
    }
}

TEST_CASE("Streaming connectivity seals delayed joins without retaining borrowed tiles", "[ThermoelasticLocalSections]")
{
    STATIC_REQUIRE_FALSE(std::is_default_constructible_v<LocalConnectivityCertificate>);
    STATIC_REQUIRE_FALSE(std::is_copy_constructible_v<LocalTopologyBuilder>);
    std::vector<std::size_t> provisional;
    LocalTopologyBuilder builder(0., {}, [&](const SectionPiece& piece) {
        CHECK_FALSE(piece.slab.footprint.native().empty());
        provisional.push_back(piece.component);
    });
    const auto identity = builder.source();
    auto split = box_tile(0, 0., .001);
    split.slabs.push_back(box_tile(0, .002, .003).slabs.front());
    builder.add_tile(split);
    split.slabs.clear();
    auto neighbor = box_tile(1, 0., .003);
    builder.add_tile(neighbor);
    neighbor.slabs.clear();
    const auto certificate = builder.finish();
    CHECK(certificate.source() == identity);
    CHECK(certificate.components() == 1);
    REQUIRE(provisional.size() == 3);
    for (const auto id : provisional) CHECK(certificate.component_mapping().at(id) == 0);
    CHECK_THROWS_AS(builder.finish(), std::invalid_argument);
    CHECK_THROWS_AS(builder.add_tile(box_tile(2, 0., .001)), std::invalid_argument);
    LocalTopologyBuilder other(0., {}, [](const SectionPiece&) {});
    CHECK(other.source() != identity);
}

TEST_CASE("Streaming pieces share proved vertical connectivity but not gaps", "[ThermoelasticLocalSections]")
{
    std::vector<std::size_t> provisional;
    LocalTopologyBuilder builder(0., {}, [&](const SectionPiece& piece) { provisional.push_back(piece.component); });
    auto tile = box_tile(0, 0., .001);
    tile.slabs.push_back(box_tile(0, .001, .002).slabs.front());
    tile.slabs.push_back(box_tile(0, .003, .004).slabs.front());
    builder.add_tile(tile);
    REQUIRE(provisional.size() == 3);
    CHECK(provisional[0] == provisional[1]);
    CHECK(provisional[0] != provisional[2]);
    CHECK_THROWS_AS(builder.finish(), std::invalid_argument);
}

TEST_CASE("Provisional representatives consolidate only already proved online joins", "[ThermoelasticLocalSections]")
{
    LocalTopologyBuilder builder(0., {}, [](const SectionPiece&) {}, {}, {}, nullptr, LocalTileOrder::Lexicographic);
    auto split = box_tile(0, 0., .001);
    split.slabs.push_back(box_tile(0, .002, .003).slabs.front());
    builder.add_tile(split);
    const auto before = builder.provisional_component_mapping();
    REQUIRE(before.size() == 2);
    CHECK(before[0] != before[1]);
    builder.add_tile(box_tile(1, 0., .003));
    const auto after = builder.provisional_component_mapping();
    REQUIRE(after.size() >= before.size());
    CHECK(after[0] == after[1]);
    for (const auto root : after) CHECK(after.at(root) == root);
    const auto certificate = builder.finish();
    for (std::size_t i = 0; i < after.size(); ++i)
        CHECK(certificate.component_mapping().at(i) == certificate.component_mapping().at(after[i]));
    CHECK_THROWS_AS(builder.provisional_component_mapping(),std::invalid_argument);
}

TEST_CASE("Provisional representative snapshots charge work and reject swallowed reentry", "[ThermoelasticLocalSections]")
{
    StreamTopologyStatistics statistics;
    LocalTopologyBuilder* active = nullptr;
    bool reenter = false;
    LocalTopologyBuilder builder(0., {}, [](const SectionPiece&) {}, [&] {
        if (reenter) {
            try { active->provisional_component_mapping(); } catch (const std::invalid_argument&) {}
        }
    }, StreamTopologyLimits{}, &statistics);
    active = &builder;
    builder.add_tile(box_tile(0,0.,.001));
    const auto before = statistics.work;
    REQUIRE(builder.provisional_component_mapping().size() == 1);
    CHECK(statistics.work > before);
    reenter = true;
    CHECK_THROWS_AS(builder.provisional_component_mapping(),std::invalid_argument);
    CHECK_THROWS_AS(builder.finish(),std::invalid_argument);
}

TEST_CASE("Failed streaming connectivity cannot publish a certificate", "[ThermoelasticLocalSections]")
{
    LocalTopologyBuilder cancelled(0., {}, [](const SectionPiece&) { throw ResourceLimit("cancel piece consumer"); });
    CHECK_THROWS_AS(cancelled.add_tile(box_tile(0, 0., .001)), ResourceLimit);
    CHECK_THROWS_AS(cancelled.finish(), std::invalid_argument);
    CHECK_THROWS_AS(cancelled.add_tile(box_tile(1, 0., .001)), std::invalid_argument);

    LocalTopologyBuilder unanchored(0., {}, [](const SectionPiece&) {});
    unanchored.add_tile(box_tile(0, .001, .002));
    CHECK_THROWS_AS(unanchored.finish(), std::invalid_argument);
    CHECK_THROWS_AS(unanchored.add_tile(box_tile(1, 0., .002)), std::invalid_argument);

    SectionTopologyLimits limits;
    limits.max_slabs = 1;
    LocalTopologyBuilder bounded(0., limits, [](const SectionPiece&) {});
    bounded.add_tile(box_tile(0, 0., .001));
    CHECK_THROWS_AS(bounded.add_tile(box_tile(1, 0., .001)), ResourceLimit);
    CHECK_THROWS_AS(bounded.finish(), std::invalid_argument);
}

TEST_CASE("Explicit free-component topology retains gaps and reports unattached solids", "[ThermoelasticLocalSections]")
{
    LocalTopologyBuilder builder(0., {}, [](const SectionPiece&) {}, {}, {}, nullptr,
        LocalTileOrder::Arbitrary, LocalAttachmentPolicy::AllowFreeComponents);
    auto tile = box_tile(0, 0., .001);
    tile.slabs.push_back(box_tile(0, .002, .003).slabs.front());
    builder.add_tile(tile);
    const auto certificate = builder.finish();
    CHECK(certificate.components() == 2);
    CHECK(certificate.free_components() == 1);
    REQUIRE(certificate.component_mapping().size() == 2);
    CHECK(certificate.component_mapping()[0] != certificate.component_mapping()[1]);
}

TEST_CASE("Streaming topology preflights dense island and face trace storage", "[ThermoelasticLocalSections]")
{
    auto tile = box_tile(0, 0., .001);
    tile.slabs.front().footprint.native().clear();
    for (coord_t x : {coord_t(100000), coord_t(400000), coord_t(700000)})
        tile.slabs.front().footprint.native().emplace_back(Polygon{Point(x, 100000), Point(x + 100000, 100000),
            Point(x + 100000, 900000), Point(x, 900000)});
    StreamTopologyLimits limits;
    limits.max_islands = 2;
    StreamTopologyStatistics statistics;
    unsigned emitted = 0;
    LocalTopologyBuilder dense(0., {}, [&](const SectionPiece&) { ++emitted; }, {}, limits, &statistics);
    CHECK_THROWS_AS(dense.add_tile(tile), ResourceLimit);
    CHECK(statistics.processed_islands == 2);
    CHECK(statistics.labels == 0);
    CHECK(statistics.face_records == 0);
    CHECK(emitted == 0);
    CHECK_THROWS_AS(dense.finish(), std::invalid_argument);

    limits.max_islands = 100;
    limits.max_trace_intervals = 1;
    LocalTopologyBuilder faces(0., {}, [&](const SectionPiece&) { ++emitted; }, {}, limits, &statistics);
    CHECK_THROWS_AS(faces.add_tile(box_tile(0, 0., .001)), ResourceLimit);
    CHECK(statistics.trace_intervals <= 1);
    CHECK(statistics.islands == 0);
    CHECK(emitted == 0);
    CHECK_THROWS_AS(faces.finish(), std::invalid_argument);
}

TEST_CASE("Swallowed reentry still invalidates streaming connectivity", "[ThermoelasticLocalSections]")
{
    LocalTopologyBuilder* active = nullptr;
    LocalTopologyBuilder builder(0., {}, [&](const SectionPiece&) {
        REQUIRE(active != nullptr);
        CHECK_THROWS_AS(active->finish(), std::invalid_argument);
    });
    active = &builder;
    CHECK_THROWS_AS(builder.add_tile(box_tile(0, 0., .001)), std::invalid_argument);
    CHECK_THROWS_AS(builder.finish(), std::invalid_argument);
}

TEST_CASE("Streamed mechanics merges late connectivity before counting physical cells", "[ThermoelasticLocalSections]")
{
    const std::vector<LocalSectionTile> tiles{box_tile(0, 0., .001), box_tile(1, 0., .001)};
    auto material = elastic_material();
    material.eigenstrain.head<3>().setConstant(-.001);
    SectionSolveLimits limits;
    limits.max_cells = 1;
    const auto topology = qualify_local_section_topology(tiles, 0.);
    const auto reference = solve_deposited_sections(topology, {material}, {0., 0., 0.}, {.002, .001, .001}, limits);
    const auto streamed = solve_local_section_tiles(tiles, 0., {material}, {0., 0., 0.}, {.002, .001, .001}, limits);
    CHECK_THAT(streamed.occupied_volume_m3, Catch::Matchers::WithinRel(reference.occupied_volume_m3, 1e-12));
    CHECK_THAT(streamed.elastic_energy_j, Catch::Matchers::WithinRel(reference.elastic_energy_j, 1e-8));
    REQUIRE(streamed.nodes_m.size() == reference.nodes_m.size());
    for (std::size_t i = 0; i < streamed.nodes_m.size(); ++i) {
        CHECK_THAT((streamed.nodes_m[i] - reference.nodes_m[i]).norm(), Catch::Matchers::WithinAbs(0., 1e-15));
        CHECK_THAT((streamed.displacement_m[i] - reference.displacement_m[i]).norm(), Catch::Matchers::WithinAbs(0., 1e-12));
        CHECK_THAT((streamed.generalized_bed_reaction_n[i] - reference.generalized_bed_reaction_n[i]).norm(),
            Catch::Matchers::WithinAbs(0., 1e-12));
    }
}

TEST_CASE("Streamed mechanics preserves separate bodies and missing sample coverage", "[ThermoelasticLocalSections]")
{
    const std::vector<LocalSectionTile> tiles{box_tile(0, 0., .001), box_tile(2, 0., .001)};
    auto material = elastic_material();
    material.eigenstrain.head<3>().setConstant(-.001);
    SectionSolveLimits limits;
    limits.max_cells = 1;
    CHECK_THROWS_AS(solve_local_section_tiles(tiles, 0., {material}, {0., 0., 0.}, {.003, .001, .001}, limits), ResourceLimit);
    limits.max_cells = 2;
    const std::vector<SectionSampleRequest> samples{
        {1, {.0005, .0005, .0005}, 0}, {2, {.0025, .0005, .0005}, 0}, {3, {.0015, .0005, .0005}, 0}};
    const auto streamed = solve_local_section_tiles(tiles, 0., {material}, {0., 0., 0.}, {.003, .001, .001}, limits, {}, {}, samples);
    const auto reference = solve_deposited_sections(qualify_local_section_topology(tiles, 0.),
        {material}, {0., 0., 0.}, {.003, .001, .001}, limits, {}, samples);
    CHECK(streamed.nodes_m.size() == 16);
    CHECK(streamed.unsampled_points == 1);
    REQUIRE(streamed.samples.size() == reference.samples.size());
    for (std::size_t i = 0; i < streamed.samples.size(); ++i) {
        CHECK(streamed.samples[i].id == reference.samples[i].id);
        CHECK_THAT((streamed.samples[i].displacement_m - reference.samples[i].displacement_m).norm(),
            Catch::Matchers::WithinAbs(0., 1e-12));
        CHECK_THAT((streamed.samples[i].stress_pa - reference.samples[i].stress_pa).norm(),
            Catch::Matchers::WithinAbs(0., 1e-9));
    }
    SectionTopologyLimits bounded;
    bounded.max_slabs = 1;
    CHECK_THROWS_AS(solve_local_section_tiles(tiles, 0., {material}, {0., 0., 0.}, {.003, .001, .001}, limits, bounded), ResourceLimit);
    CHECK_THROWS_AS(solve_local_section_tiles(tiles, 0., {material}, {0., 0., 0.}, {.003, .001, .001}, limits, {},
        [] { throw ResourceLimit("cancel streaming mechanics"); }), ResourceLimit);
}

TEST_CASE("Uncertain interlayer contacts cannot create an attachment", "[ThermoelasticLocalSections]")
{
    auto tile = box_tile(0, 0., .001);
    tile.slabs.front().boundary_uncertainty_nm = 4;
    auto upper = box_tile(0, .001, .002).slabs.front();
    upper.boundary_uncertainty_nm = 4;
    tile.slabs.push_back(upper);
    CHECK(qualify_local_section_topology({tile}, 0.).components() == 1);
    for (auto& point : tile.slabs.front().footprint.native().front().contour.points)
        if (point.x() == 1000000) point.x() = 500002;
    for (auto& point : tile.slabs.back().footprint.native().front().contour.points)
        if (point.x() == 0) point.x() = 499998;
    try {
        qualify_local_section_topology({tile}, 0.);
        FAIL("An uncertainty-sensitive interlayer contact was admitted.");
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        CHECK(message.find("first_z_nm=0,1000000") != std::string::npos);
        CHECK(message.find("second_z_nm=1000000,2000000") != std::string::npos);
        CHECK(message.find("uncertainty_nm=4,4") != std::string::npos);
        CHECK(message.find("exact_nominal_contact=1") != std::string::npos);
        CHECK(message.find("exact_offset_contact=0") != std::string::npos);
    }
}

TEST_CASE("An unresolved input island is reported rather than discarded by area", "[ThermoelasticLocalSections]")
{
    auto tile = box_tile(0, 0., .001);
    tile.slabs.front().boundary_uncertainty_nm = 7;
    tile.slabs.front().footprint.native() = {ExPolygon(Polygon{{100,100},{104,101},{102,101}})};
    try {
        qualify_local_section_topology({tile}, 0.);
        FAIL("An unresolved positive-area island was discarded.");
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        CHECK(message.find("area_nm2=1") != std::string::npos);
        CHECK(message.find("raw_input_exact_intersection_empty=0") != std::string::npos);
        CHECK(message.find("raw_input_exact_outside_empty=1") != std::string::npos);
        CHECK(message.find("material_union_exact_intersection_empty=0") != std::string::npos);
        CHECK(message.find("material_union_exact_outside_empty=1") != std::string::npos);
    }
}

TEST_CASE("A single certified occupied polygon retains its narrow boundary chain", "[ThermoelasticLocalSections]")
{
    auto tile = box_tile(0, 0., .001);
    const ExPolygon shape(Polygon{{0,0},{100,0},{100,49},{200,49},{200,0},{300,0},
        {300,100},{200,100},{200,51},{100,51},{100,100},{0,100}});
    tile.slabs.front().footprint.native() = {shape};
    const auto topology = qualify_local_section_topology({tile}, 0.);
    REQUIRE(topology.pieces().size() == 1);
    REQUIRE(topology.pieces().front().slab.footprint.native().size() == 1);
    CHECK(topology.pieces().front().slab.footprint.native().front().contour.points == shape.contour.points);
}

TEST_CASE("Connected uncertainty cores preserve original holes and integration moments", "[ThermoelasticLocalSections]")
{
    const bool merge_holes=GENERATE(false,true);
    auto tile=box_tile(0,0.,.001);
    tile.slabs.front().boundary_uncertainty_nm=7;
    auto& shape=tile.slabs.front().footprint.native().front();
    if (merge_holes) {
        shape.holes.emplace_back(Polygon{{200000,200000},{200000,400000},{400000,400000},{400000,200000}});
        shape.holes.emplace_back(Polygon{{400010,200000},{400010,400000},{600000,400000},{600000,200000}});
    } else {
        // The5nm outer web disappears under the9nm contact-search erosion,
        // opening the hole without disconnecting the remaining solid.
        shape.holes.emplace_back(Polygon{{5,200000},{5,400000},{200000,400000},{200000,200000}});
    }
    const auto original=tile.slabs.front().footprint.native();
    const auto core=offset_ex(original,-9.f,jtRound,.25);
    REQUIRE(core.size()==1);
    REQUIRE(core.front().holes.size()<shape.holes.size());
    const auto topology=qualify_local_section_topology({tile},0.);
    REQUIRE(topology.components()==1);
    REQUIRE(topology.pieces().size()==1);
    const auto& retained=topology.pieces().front().slab;
    CHECK(retained.footprint.native()==original);
    CHECK(retained.boundary_uncertainty_nm==7);
    const std::array<coord_t,3> origin{0,0,0},size{1000000,1000000,1000000};
    const auto before=SectionIntegration::volume_moments(original,origin,size,0,1000000);
    const auto after=SectionIntegration::volume_moments(retained.footprint.native(),origin,size,0,1000000);
    for (std::size_t i=0;i<before.size();++i)
        CHECK_THAT(after[i],Catch::Matchers::WithinAbs(before[i],1e-24));
}

TEST_CASE("Uncertainty cores still reject split necks and erased solids", "[ThermoelasticLocalSections]")
{
    auto tile=box_tile(0,0.,.001);
    tile.slabs.front().boundary_uncertainty_nm=7;
    tile.slabs.front().footprint.native()={ExPolygon(Polygon{{100,100},{200,100},{200,149},{300,149},{300,100},{400,100},
        {400,200},{300,200},{300,151},{200,151},{200,200},{100,200}})};
    REQUIRE(offset_ex(tile.slabs.front().footprint.native(),-9.f,jtRound,.25).size()==2);
    CHECK_THROWS_AS(qualify_local_section_topology({tile},0.),std::invalid_argument);
    tile.slabs.front().footprint.native()={ExPolygon(Polygon{{100,100},{104,100},{104,104},{100,104}})};
    REQUIRE(offset_ex(tile.slabs.front().footprint.native(),-9.f,jtRound,.25).empty());
    CHECK_THROWS_AS(qualify_local_section_topology({tile},0.),std::invalid_argument);
}

TEST_CASE("A point touching tiny island is not hidden by its larger contour", "[ThermoelasticLocalSections]")
{
    auto tile = box_tile(0, 0., .001);
    tile.slabs.front().boundary_uncertainty_nm = 7;
    // The repeated corner joins two positive oriented loops. Eroding the
    // unsplit chain retains the large box, but not its genuine tiny island.
    tile.slabs.front().footprint.native() = {ExPolygon(Polygon{{100,100},{1000,100},
        {1002,99},{1004,99},{1000,100},{1000,1000},{100,1000}})};
    REQUIRE_FALSE(offset_ex(tile.slabs.front().footprint.native(),-9.f,jtRound,.25).empty());
    CHECK_THROWS_AS(qualify_local_section_topology({tile},0.),std::invalid_argument);
}

TEST_CASE("Local component assignment preserves bent shared material boundaries", "[ThermoelasticLocalSections]")
{
    const Point a(0,0), q(107087,259824), b(214844,506414), c(129002,312996), d(107018,259666), left(-100000,100000);
    const ExPolygons first{ExPolygon(Polygon{a,b,q})}, second{ExPolygon(Polygon{a,q,c,left})};
    const auto clipped = exact_clip_reference_sections({first,second},{ExPolygons{ExPolygon(Polygon{a,b,d})},second});
    LocalSectionTile tile{{0,0},{-1000000,-1000000},{1000000,1000000},
        {{0.,.001,0,clipped.sections[0]},{0.,.001,1,clipped.sections[1]}}};
    const auto topology = qualify_local_section_topology({tile},0.);
    CHECK(topology.components() == 1);
    std::array<ExPolygons,2> retained;
    for (const auto& piece : topology.pieces()) {
        REQUIRE(piece.slab.material < retained.size());
        auto& material = retained[piece.slab.material];
        material.insert(material.end(),piece.slab.footprint.native().begin(),piece.slab.footprint.native().end());
    }
    for (std::size_t i=0; i<retained.size(); ++i)
        CHECK(exact_section_empty(retained[i],clipped.sections[i],ExactSectionOperation::SymmetricDifference).empty);
}

TEST_CASE("Offline topology work reservations retain independent memory limits", "[ThermoelasticLocalSections]")
{
    StreamTopologyLimits limits;
    CHECK(limits.max_work == 20000000);
    limits.max_work = 4000000000ULL;
    LocalTopologyBuilder admitted(0.,{},[](const SectionPiece&) {},{},limits);
    admitted.add_tile(box_tile(0,0.,.001));
    CHECK(admitted.finish().components() == 1);
    ++limits.max_work;
    CHECK_THROWS_AS(LocalTopologyBuilder(0.,{},[](const SectionPiece&) {},{},limits),std::invalid_argument);
    --limits.max_work;
    limits.max_live_points = 1;
    LocalTopologyBuilder bounded(0.,{},[](const SectionPiece&) {},{},limits);
    CHECK_THROWS_AS(bounded.add_tile(box_tile(0,0.,.001)),ResourceLimit);
}

TEST_CASE("Uncertain interlayer edges never replace independent certified connectivity", "[ThermoelasticLocalSections]")
{
    auto tile = box_tile(0, 0., .001);
    tile.slabs.front().boundary_uncertainty_nm = 4;
    auto upper = box_tile(0, .001, .002).slabs.front();
    upper.boundary_uncertainty_nm = 4;
    tile.slabs.push_back(upper);
    const coord_t upper_minimum = GENERATE(coord_t(499998), coord_t(500006));
    for (auto& point : tile.slabs.front().footprint.native().front().contour.points)
        if (point.y() == 1000000) point.y() = 500002;
    for (auto& point : tile.slabs.back().footprint.native().front().contour.points)
        if (point.y() == 0) point.y() = upper_minimum;
    auto bridge = box_tile(1, 0., .002);
    CHECK_THROWS_AS(qualify_local_section_topology({tile}, 0.), std::invalid_argument);
    CHECK(qualify_local_section_topology({tile, bridge}, 0.).components() == 1);
    CHECK(qualify_local_section_topology({bridge, tile}, 0.).components() == 1);
    StreamTopologyStatistics statistics;
    LocalTopologyBuilder builder(0., {}, [](const SectionPiece&) {}, {}, StreamTopologyLimits{},
        &statistics, LocalTileOrder::Lexicographic);
    builder.add_tile(tile);
    // No uncertain edge has merged the two vertical labels before the
    // independent neighboring tile supplies its certified lateral path.
    const auto provisional = builder.provisional_component_mapping();
    REQUIRE(provisional.size() == 2);
    CHECK(provisional[0] != provisional[1]);
    builder.add_tile(bridge);
    CHECK(builder.finish().components() == 1);
    CHECK(statistics.redundant_contacts == 1);
}

TEST_CASE("Local section tiles preserve the union volume of overlapping beads", "[ThermoelasticLocalSections]")
{
    const std::vector<DepositedSegment> beads{
        {{-.001, 0., .001}, {.001, 0., .001}, .001, .001, 0},
        {{0., 0., .001}, {.002, 0., .001}, .001, .001, 0}};
    const auto reference = layered_deposition(beads, 0.);
    const auto result = local_deposition_sections(beads, {0., 0.}, {.001, .001}, 0., .001);
    CHECK(result.piecewise_constant);
    CHECK(result.tiles.size() > 1);
    CHECK_THAT(result.volume_m3, Catch::Matchers::WithinRel(reference.occupied_volume_m3, 1e-8));
    for (const auto& tile : result.tiles)
        for (const auto& slab : tile.slabs)
            for (const auto& polygon : slab.footprint.native())
                for (const auto& point : polygon.contour.points) {
                    CHECK(point.x() >= tile.minimum[0]);
                    CHECK(point.x() <= tile.maximum[0]);
                    CHECK(point.y() >= tile.minimum[1]);
                    CHECK(point.y() <= tile.maximum[1]);
                }
}

TEST_CASE("Diagonal bead indexing avoids empty bounding rectangle tiles", "[ThermoelasticLocalSections]")
{
    for (const double direction : {-1.,1.}) {
        // The old rectangle covered more than ten thousand tiles. A narrow
        // diagonal occupies only a few hundred, including conservative caps.
        const std::vector<DepositedSegment> beads{
            {{0.,0.,.0002},{.1,direction*.1,.0002},.0004,.0002,0}};
        LocalSectionLimits limits;
        limits.max_tiles=400;
        LocalSectionCursor cursor(beads,{0.,0.},{.001,.001},0.,.0002,limits);
        CHECK(cursor.index_usage().tiles<400);
        CHECK(cursor.index_usage().references==cursor.index_usage().tiles);
        const auto reference=layered_deposition(beads,0.);
        const auto tiled=local_deposition_sections(beads,{0.,0.},{.001,.001},0.,.0002,limits);
        CHECK_THAT(tiled.volume_m3,Catch::Matchers::WithinRel(reference.occupied_volume_m3,1e-7));
        limits.max_index_work=1;
        CHECK_THROWS_AS(LocalSectionCursor(beads,{0.,0.},{.001,.001},0.,.0002,limits),ResourceLimit);
    }
    // Indexing also works close to the supported coordinate limit.
    const std::vector<DepositedSegment> translated{
        {{9.,9.,.0002},{9.1,9.1,.0002},.0004,.0002,0}};
    LocalSectionLimits limits;
    limits.max_tiles=400;
    LocalSectionCursor shifted(translated,{9.,9.},{.001,.001},0.,.0002,limits);
    CHECK(shifted.index_usage().tiles<400);
    // This long segment produces cross products beyond int64, despite using
    // only a few hundred relevant coarse tiles.
    const std::vector<DepositedSegment> long_bead{
        {{-4.5,-4.5,.0002},{4.5,4.5,.0002},.0004,.0002,0}};
    LocalSectionCursor wide(long_bead,{0.,0.},{.1,.1},0.,.0002,limits);
    CHECK(wide.index_usage().tiles<400);
}

TEST_CASE("Segment box indexing retains tangent candidate tiles", "[ThermoelasticLocalSections]")
{
    const std::vector<DepositedSegment> beads{
        {{0.,0.,.0002},{.003,.003,.0002},.001,.0002,0}};
    LocalSectionCursor cursor(beads,{0.,0.},{.001,.001},0.,.0002);
    bool tangent=false;
    LocalSectionPacketVisitor visitor;
    visitor.begin_tile=[&](const LocalSectionTile& tile) {
        // The expanded box touches the centerline at (1.5,1.5)mm.
        if (tile.index==std::array<int,2>{{0,2}}) tangent=true;
    };
    visitor.band=[](const LocalSectionTile&) {};
    visitor.end_tile=[] {};
    while (cursor.advance(visitor)) {}
    CHECK(tangent);
}

TEST_CASE("Fine tile studies retain independent index storage bounds", "[ThermoelasticLocalSections]")
{
    const std::vector<DepositedSegment> beads{
        {{0.,0.,.0002},{.1,0.,.0002},.00002,.0002,0}};
    LocalSectionLimits limits;
    CHECK_THROWS_AS(LocalSectionCursor(beads,{0.,0.},{.00002,.00002},0.,.0002,limits),ResourceLimit);
    limits.max_tiles=16384;
    LocalSectionCursor cursor(beads,{0.,0.},{.00002,.00002},0.,.0002,limits);
    CHECK(cursor.index_usage().tiles>4096);
    CHECK(cursor.index_usage().tiles<=limits.max_tiles);
    limits.max_tile_references=1;
    CHECK_THROWS_AS(LocalSectionCursor(beads,{0.,0.},{.00002,.00002},0.,.0002,limits),ResourceLimit);
    limits.max_tiles=16385;
    CHECK_THROWS_AS(LocalSectionCursor(beads,{0.,0.},{.00002,.00002},0.,.0002,limits),std::invalid_argument);
}

TEST_CASE("Local tile traversal preserves retained geometry and cumulative budgets", "[ThermoelasticLocalSections]")
{
    const std::vector<DepositedSegment> beads{
        {{-.001, 0., .001}, {.003, 0., .001}, .001, .001, 0},
        {{0., 0., .002}, {.001, 0., .002}, .001, .001, 0}};
    const auto retained = local_deposition_sections(beads, {0., 0.}, {.001, .001}, 0., .0005);
    std::size_t visited = 0;
    const auto metrics = visit_local_deposition_sections(beads, {0., 0.}, {.001, .001}, 0., .0005,
        [&](const LocalSectionTile& tile) {
            REQUIRE(visited < retained.tiles.size());
            const auto& expected = retained.tiles[visited++];
            CHECK(tile.index == expected.index);
            CHECK(tile.minimum == expected.minimum);
            CHECK(tile.maximum == expected.maximum);
            REQUIRE(tile.slabs.size() == expected.slabs.size());
            for (std::size_t i = 0; i < tile.slabs.size(); ++i) {
                CHECK(tile.slabs[i].footprint.native() == expected.slabs[i].footprint.native());
                CHECK(tile.slabs[i].material == expected.slabs[i].material);
                CHECK_THAT(tile.slabs[i].bottom_m, Catch::Matchers::WithinAbs(expected.slabs[i].bottom_m, 1e-15));
                CHECK_THAT(tile.slabs[i].top_m, Catch::Matchers::WithinAbs(expected.slabs[i].top_m, 1e-15));
            }
        });
    CHECK(visited == retained.tiles.size());
    CHECK(metrics.sections == retained.sections);
    CHECK(metrics.segment_visits == retained.segment_visits);
    CHECK(metrics.boundary_points == retained.boundary_points);
    CHECK(metrics.piecewise_constant == retained.piecewise_constant);
    CHECK_THAT(metrics.volume_m3, Catch::Matchers::WithinRel(retained.volume_m3, 1e-12));
    LocalSectionLimits bounded;
    bounded.max_boundary_points = retained.boundary_points - 1;
    CHECK_THROWS_AS(visit_local_deposition_sections(beads, {0., 0.}, {.001, .001}, 0., .0005,
        [](const LocalSectionTile&) {}, bounded), ResourceLimit);
    CHECK_THROWS_AS(visit_local_deposition_sections(beads, {0., 0.}, {.001, .001}, 0., .0005,
        [](const LocalSectionTile&) { throw ResourceLimit("cancel tile consumer"); }), ResourceLimit);
    CHECK_THROWS_AS(visit_local_deposition_sections(beads, {0., 0.}, {.001, .001}, 0., .0005, {}), std::invalid_argument);
}

TEST_CASE("Resumable section cursors preserve geometry and cumulative budgets when interleaved", "[ThermoelasticLocalSections]")
{
    const std::vector<DepositedSegment> beads{
        {{.0005,.0005,.0005},{.0035,.0005,.0005},.001,.0005,0,BeadSection::RoundedRectangle}};
    const auto reference=local_deposition_sections(beads,{0.,0.},{.001,.001},0.,.0001);
    LocalSectionCursor first(beads,{0.,0.},{.001,.001},0.,.0001);
    LocalSectionCursor second(beads,{0.,0.},{.001,.001},0.,.0001);
    std::vector<LocalSectionTile> reconstructed;
    LocalSectionPacketVisitor collect;
    collect.begin_tile=[&](const LocalSectionTile& h) { reconstructed.push_back(h); };
    collect.band=[&](const LocalSectionTile& p) {
        auto& slabs=reconstructed.back().slabs;
        slabs.insert(slabs.end(),p.slabs.begin(),p.slabs.end());
    };
    collect.end_tile=[] {};
    LocalSectionPacketVisitor discard{[](const LocalSectionTile&) {},[](const LocalSectionTile&) {},[] {}};
    while (first.next_tile_index()) {
        CHECK(first.next_tile_index()==second.next_tile_index());
        REQUIRE(first.advance(collect));
        REQUIRE(second.advance(discard));
    }
    CHECK_FALSE(first.advance(collect));
    CHECK_FALSE(second.next_tile_index());
    // The index may include tangent-only empty tiles; retained geometry does not.
    reconstructed.erase(std::remove_if(reconstructed.begin(),reconstructed.end(),
        [](const auto& tile) { return tile.slabs.empty(); }),reconstructed.end());
    REQUIRE(reconstructed.size()==reference.tiles.size());
    for (std::size_t i=0;i<reconstructed.size();++i) {
        const auto& a=reconstructed[i]; const auto& b=reference.tiles[i];
        CHECK(a.index==b.index);
        REQUIRE(a.slabs.size()==b.slabs.size());
        for (std::size_t j=0;j<a.slabs.size();++j) {
            CHECK(a.slabs[j].footprint.native()==b.slabs[j].footprint.native());
            CHECK_THAT(a.slabs[j].bottom_m,Catch::Matchers::WithinAbs(b.slabs[j].bottom_m,1e-15));
            CHECK_THAT(a.slabs[j].top_m,Catch::Matchers::WithinAbs(b.slabs[j].top_m,1e-15));
        }
    }
    CHECK(first.metrics().sections==reference.sections);
    CHECK(first.metrics().segment_visits==reference.segment_visits);
    CHECK_THAT(first.metrics().volume_m3,Catch::Matchers::WithinAbs(reference.volume_m3,1e-18));
    LocalSectionLimits limited; limited.max_sections=5;
    LocalSectionCursor bounded(beads,{0.,0.},{.001,.001},0.,.0001,limited);
    CHECK_THROWS_AS(([&] { while (bounded.advance(discard)) {} }()),ResourceLimit);
    CHECK_THROWS_AS(bounded.metrics(),std::logic_error);
    LocalSectionCursor rejected(beads,{0.,0.},{.001,.001},0.,.0001);
    CHECK_THROWS_AS(rejected.advance({}),std::invalid_argument);
    CHECK_THROWS_AS(rejected.next_tile_index(),std::logic_error);
}

TEST_CASE("Band packets preserve complete material partitions and empty gaps", "[ThermoelasticLocalSections]")
{
    const std::vector<DepositedSegment> beads{
        {{.001, .001, .001}, {.002, .001, .001}, .0005, .001, 0},
        {{.004, .001, .001}, {.005, .001, .001}, .0005, .001, 1},
        {{.001, .001, .003}, {.002, .001, .003}, .0005, .001, 0}};
    const auto retained = local_deposition_sections(beads, {0., 0.}, {.01, .01}, 0., .0002);
    REQUIRE(retained.tiles.size() == 1);
    LocalSectionTile reconstructed;
    unsigned begins = 0, bands = 0, ends = 0;
    LocalSectionPacketVisitor visitor;
    visitor.begin_tile = [&](const LocalSectionTile& header) {
        ++begins;
        CHECK(header.slabs.empty());
        reconstructed = header;
    };
    visitor.band = [&](const LocalSectionTile& packet) {
        ++bands;
        CHECK(packet.index == reconstructed.index);
        REQUIRE_FALSE(packet.slabs.empty());
        CHECK(packet.slabs.size() == (bands == 1 ? 2 : 1));
        for (const auto& slab : packet.slabs) {
            CHECK_THAT(slab.bottom_m, Catch::Matchers::WithinAbs(packet.slabs.front().bottom_m, 1e-15));
            CHECK_THAT(slab.top_m, Catch::Matchers::WithinAbs(packet.slabs.front().top_m, 1e-15));
        }
        reconstructed.slabs.insert(reconstructed.slabs.end(), packet.slabs.begin(), packet.slabs.end());
    };
    visitor.end_tile = [&] { ++ends; };
    const auto metrics = visit_local_deposition_packets(beads, {0., 0.}, {.01, .01}, 0., .0002, visitor);
    CHECK(begins == 1);
    CHECK(bands == 2);
    CHECK(ends == 1);
    CHECK(metrics.boundary_points == retained.boundary_points);
    CHECK(metrics.sections == retained.sections);
    begins = bands = ends = 0;
    LocalSectionCursor cursor(beads,{0.,0.},{.01,.01},0.,.0002);
    while (cursor.next_tile_index()) {
        const auto before = bands;
        REQUIRE(cursor.advance_packet(visitor));
        CHECK(bands <= before+1);
    }
    CHECK(begins == 1);
    CHECK(bands == 2);
    CHECK(ends == 1);
    CHECK(cursor.metrics().boundary_points == metrics.boundary_points);
    REQUIRE(reconstructed.slabs.size() == retained.tiles.front().slabs.size());
    for (std::size_t i = 0; i < reconstructed.slabs.size(); ++i) {
        const auto& expected = retained.tiles.front().slabs[i];
        CHECK(reconstructed.slabs[i].footprint.native() == expected.footprint.native());
        CHECK(reconstructed.slabs[i].material == expected.material);
        CHECK_THAT(reconstructed.slabs[i].bottom_m, Catch::Matchers::WithinAbs(expected.bottom_m, 1e-15));
        CHECK_THAT(reconstructed.slabs[i].top_m, Catch::Matchers::WithinAbs(expected.top_m, 1e-15));
    }
    visitor.band = [](const LocalSectionTile&) { throw ResourceLimit("cancel packet"); };
    ends = 0;
    CHECK_THROWS_AS(visit_local_deposition_packets(beads, {0., 0.}, {.01, .01}, 0., .0002, visitor), ResourceLimit);
    CHECK(ends == 0);
    CHECK_THROWS_AS(visit_local_deposition_packets(beads, {0., 0.}, {.01, .01}, 0., .0002, {}), std::invalid_argument);
}

TEST_CASE("Deposition streaming matches retained geometry mechanics and sample coverage", "[ThermoelasticLocalSections]")
{
    const std::vector<DepositedSegment> beads{
        {{0., 0., .001}, {.002, 0., .001}, .001, .001, 0},
        {{0., 0., .002}, {.002, 0., .002}, .001, .001, 1}};
    auto first = elastic_material(), second = elastic_material();
    first.eigenstrain.head<3>().setConstant(-.001);
    second.eigenstrain.head<3>().setConstant(-.002);
    second.elasticity_pa *= .4;
    const std::vector<SectionMaterial> materials{first, second};
    LocalSectionStreamSettings settings;
    settings.origin_m = {-.001, -.001, 0.};
    settings.basis_spacing_m = {.001, .001, .001};
    settings.tile_size_m = {.001, .001};
    settings.maximum_height_m = .0005;
    const std::vector<SectionSampleRequest> samples{
        {1, {.0012, .0001, .0005}, 0}, {2, {.0012, .0001, .0015}, 1}, {3, {.001, .0008, .0005}, 0}};
    const auto geometry = local_deposition_sections(beads, settings.origin_m.head<2>(), settings.tile_size_m,
        0., settings.maximum_height_m);
    const auto reference = solve_deposited_sections(qualify_local_section_topology(geometry.tiles, 0.),
        materials, settings.origin_m, settings.basis_spacing_m, {}, {}, samples);
    LocalSectionStatistics progress;
    const auto streamed = solve_local_deposition_sections(beads, materials, settings, {}, {}, nullptr, samples,
        nullptr, nullptr, &progress);
    CHECK(progress.completed_tiles == progress.indexed_tiles);
    CHECK(progress.completed_tiles >= geometry.tiles.size());
    CHECK(streamed.geometry.sections == geometry.sections);
    CHECK(streamed.geometry.boundary_points == geometry.boundary_points);
    CHECK(streamed.geometry.segment_visits == geometry.segment_visits);
    CHECK_THAT(streamed.geometry.volume_m3, Catch::Matchers::WithinRel(geometry.volume_m3, 1e-12));
    CHECK_THAT(streamed.field.elastic_energy_j, Catch::Matchers::WithinRel(reference.elastic_energy_j, 1e-8));
    CHECK_THAT(streamed.field.bed_area_m2, Catch::Matchers::WithinRel(reference.bed_area_m2, 1e-12));
    CHECK(streamed.field.unsampled_points == 1);
    REQUIRE(streamed.field.samples.size() == reference.samples.size());
    for (std::size_t i = 0; i < streamed.field.samples.size(); ++i) {
        CHECK(streamed.field.samples[i].id == reference.samples[i].id);
        CHECK_THAT((streamed.field.samples[i].stress_pa - reference.samples[i].stress_pa).norm(),
            Catch::Matchers::WithinAbs(0., 1e-8));
    }
    // Tiles bound traversal memory, not the mechanical basis or physical bonds.
    // Changing their seams must not change the solved field away from seams.
    auto other_tiles = settings;
    other_tiles.tile_size_m = {.0025, .0025};
    const auto retiled = solve_local_deposition_sections(beads, materials, other_tiles, {}, {}, nullptr, samples);
    CHECK_THAT(retiled.field.occupied_volume_m3,
        Catch::Matchers::WithinRel(streamed.field.occupied_volume_m3, 1e-10));
    CHECK_THAT(retiled.field.elastic_energy_j,
        Catch::Matchers::WithinRel(streamed.field.elastic_energy_j, 1e-8));
    REQUIRE(retiled.field.samples.size() == streamed.field.samples.size());
    for (std::size_t i = 0; i < streamed.field.samples.size(); ++i) {
        CHECK(retiled.field.samples[i].id == streamed.field.samples[i].id);
        CHECK_THAT((retiled.field.samples[i].stress_pa - streamed.field.samples[i].stress_pa).norm(),
            Catch::Matchers::WithinAbs(0., 1e-8));
    }
    settings.geometry.max_boundary_points = geometry.boundary_points - 1;
    CHECK_THROWS_AS(solve_local_deposition_sections(beads, materials, settings, {}, {}, nullptr, {},
        nullptr, nullptr, &progress), ResourceLimit);
    CHECK(progress.completed_tiles < progress.indexed_tiles);
    settings.geometry.max_boundary_points = geometry.boundary_points;
    CHECK_THROWS_AS(solve_local_deposition_sections(beads, materials, settings, {},
        [] { throw ResourceLimit("cancel complete source traversal"); }), ResourceLimit);
    CHECK_THROWS_AS(solve_local_deposition_sections(beads, materials, settings,
        [](double) -> std::optional<SectionOwnership> { throw ResourceLimit("cancel source ownership"); }), ResourceLimit);
    CHECK_THROWS_AS(solve_local_deposition_sections(beads, materials, LocalSectionStreamSettings{}), std::invalid_argument);
}

TEST_CASE("Larger cumulative geometry budgets require explicit streaming consumption", "[ThermoelasticLocalSections]")
{
    const std::vector<DepositedSegment> beads{
        {{.001, .001, .0002}, {.002, .001, .0002}, .0004, .0002, 0}};
    LocalSectionLimits limits;
    limits.max_boundary_points = 128000000;
    limits.max_index_work = 500000000;
    CHECK_THROWS_AS(local_deposition_sections(beads, {0., 0.}, {.01, .01}, 0., .0002, limits),
        std::invalid_argument);
    std::size_t bands = 0;
    LocalSectionPacketVisitor visitor;
    visitor.begin_tile = [](const LocalSectionTile&) {};
    visitor.band = [&](const LocalSectionTile&) { ++bands; };
    visitor.end_tile = [] {};
    const auto result = visit_local_deposition_packets(beads, {0., 0.}, {.01, .01}, 0., .0002, visitor, limits);
    CHECK(bands == 1);
    CHECK(result.boundary_points > 0);
    ++limits.max_boundary_points;
    CHECK_THROWS_AS(visit_local_deposition_packets(beads, {0., 0.}, {.01, .01}, 0., .0002, visitor, limits),
        std::invalid_argument);
    CHECK(bands == 1);
    limits.max_boundary_points = 128000000;
    ++limits.max_index_work;
    CHECK_THROWS_AS(visit_local_deposition_packets(beads, {0., 0.}, {.01, .01}, 0., .0002, visitor, limits),
        std::invalid_argument);
    CHECK(bands == 1);
}


TEST_CASE("Explicit partition envelopes preserve source beads and refuse invented bed contact", "[ThermoelasticLocalSections]")
{
    const DepositedSegment bead{{.001,.001,.000203},{.002,.001,.000203},.0004,.000203,0};
    const std::vector<DepositedSegment> source{bead};
    const auto partition = make_section_partition_envelopes(source, 0., .00001);
    REQUIRE(partition.segments.size() == 1);
    CHECK(partition.segments[0].start_m.head<2>() == bead.start_m.head<2>());
    CHECK_THAT(partition.maximum_shift_m, Catch::Matchers::WithinAbs(7e-6,1e-15));
    CHECK_THAT(source[0].start_m.z(), Catch::Matchers::WithinAbs(.000203,1e-15));
    LocalSectionLimits limits;
    limits.geometry.half_native_section_planes = true;
    const auto sampled = local_deposition_sections(source,{0.,0.},{.01,.01},0.,.00005,
        limits,{}, {},1,nullptr,&partition.segments);
    REQUIRE_FALSE(sampled.tiles.empty());
    CHECK(sampled.uses_partition_envelopes);
    limits.partition_increment_m=.00001;
    const auto automatic=local_deposition_sections(source,{0.,0.},{.01,.01},0.,.00005,limits);
    CHECK(automatic.uses_partition_envelopes);
    CHECK_FALSE(automatic.piecewise_constant);
    CHECK_THAT(automatic.volume_m3,Catch::Matchers::WithinAbs(sampled.volume_m3,1e-20));
    REQUIRE(automatic.tiles.size()==sampled.tiles.size());
    for (std::size_t i=0;i<automatic.tiles.size();++i) {
        REQUIRE(automatic.tiles[i].slabs.size()==sampled.tiles[i].slabs.size());
        for (std::size_t j=0;j<automatic.tiles[i].slabs.size();++j)
            CHECK(automatic.tiles[i].slabs[j].footprint.native()==sampled.tiles[i].slabs[j].footprint.native());
    }
    CHECK_THROWS_AS(local_deposition_sections(source,{0.,0.},{.01,.01},0.,.00005,
        limits,{}, {},1,nullptr,&partition.segments),std::invalid_argument);
    limits.partition_increment_m=0.;
    for (const auto& tile : sampled.tiles)
        for (const auto& slab : tile.slabs) {
            const auto reference = section_deposition(source,.5*(slab.bottom_m+slab.top_m),0.,limits.geometry);
            REQUIRE(reference.regions.size() == 1);
            CHECK(slab.footprint.native() == reference.regions[0].footprint.native());
        }
    auto floating = bead;
    floating.height_m = .000202; //1µm real gap cannot become bed contact.
    CHECK_THROWS_AS(make_section_partition_envelopes({floating},0.,.00001),std::invalid_argument);
    floating.height_m = .000204;
    CHECK_THROWS_AS(make_section_partition_envelopes({floating},0.,.00001),std::invalid_argument);
    CHECK_THROWS_AS(make_section_partition_envelopes({},0.,.00001),std::invalid_argument);
    CHECK_THROWS_AS(make_section_partition_envelopes(source,0.,1e-9),std::invalid_argument);
    CHECK_THROWS_AS(make_section_partition_envelopes(source,0.,.002),std::invalid_argument);
    CHECK_THROWS_AS(make_section_partition_envelopes(source,0.,.00001,[] {throw ResourceLimit("cancel");}),ResourceLimit);
}

TEST_CASE("Unrelated tile heights do not multiply every local section", "[ThermoelasticLocalSections]")
{
    std::vector<DepositedSegment> beads;
    for (unsigned i = 0; i < 12; ++i) {
        const double x = .01 * i;
        const double height = .001 + .00001 * i;
        beads.push_back({{x + .002, .003, height}, {x + .003, .003, height}, .001, height, 0});
    }
    const auto result = local_deposition_sections(beads, {0., 0.}, {.01, .01}, 0., .01);
    REQUIRE(result.tiles.size() == beads.size());
    CHECK(result.sections == beads.size());
    CHECK(result.segment_visits == beads.size());
}

TEST_CASE("Identical sampled partitions share storage without changing mechanical response", "[ThermoelasticLocalSections]")
{
    const std::vector<DepositedSegment> beads{
        {{.002, .003, .002}, {.004, .003, .002}, .001, .002, 0, BeadSection::RectangularEnvelope}};
    const auto coarse = local_deposition_sections(beads, {0., 0.}, {.01, .01}, 0., .002);
    LocalSectionLimits limits;
    limits.max_boundary_points = coarse.boundary_points;
    LocalSectionStatistics stats;
    const auto fine = local_deposition_sections(beads, {0., 0.}, {.01, .01}, 0., .0001,
        limits, {}, {}, 1, &stats);
    CHECK(fine.sections == 20);
    CHECK(fine.boundary_points == coarse.boundary_points);
    CHECK(stats.retained_slabs == 1);
    CHECK(stats.retained_points == fine.boundary_points);
    REQUIRE(fine.tiles.size() == 1);
    REQUIRE(fine.tiles.front().slabs.size() == 1);
    CHECK_THAT(fine.volume_m3, Catch::Matchers::WithinRel(coarse.volume_m3, 1e-12));
    auto split = fine.tiles;
    const auto prism = split.front().slabs.front();
    split.front().slabs.clear();
    for (unsigned i = 0; i < 20; ++i) {
        auto slab = prism;
        slab.bottom_m = double(i) * .0001;
        slab.top_m = double(i + 1) * .0001;
        split.front().slabs.push_back(std::move(slab));
    }
    auto material = elastic_material();
    material.eigenstrain.head<3>().setConstant(-.001);
    const auto compact = solve_deposited_sections(qualify_local_section_topology(fine.tiles, 0.),
        {material}, {0., 0., 0.}, {.001, .001, .0005});
    const auto expanded = solve_deposited_sections(qualify_local_section_topology(split, 0.),
        {material}, {0., 0., 0.}, {.001, .001, .0005});
    REQUIRE(compact.elastic_energy_j > 0.);
    CHECK_THAT(compact.elastic_energy_j, Catch::Matchers::WithinRel(expanded.elastic_energy_j, 1e-8));
    CHECK_THAT(compact.occupied_volume_m3, Catch::Matchers::WithinRel(expanded.occupied_volume_m3, 1e-12));
    const auto tensile_reaction = [](const SectionMechanicalResult& field) {
        double demand = 0.;
        for (const auto& force : field.generalized_bed_reaction_n)
            demand += std::max(0., -force.z());
        return demand;
    };
    CHECK_THAT(tensile_reaction(compact), Catch::Matchers::WithinAbs(tensile_reaction(expanded), 1e-12));
}

TEST_CASE("Partition compaction preserves empty intervals and material changes", "[ThermoelasticLocalSections]")
{
    const DepositedSegment lower{{.002, .003, .001}, {.004, .003, .001}, .001, .001, 0,
        BeadSection::RectangularEnvelope};
    auto upper = lower;
    upper.start_m.z() = upper.end_m.z() = .003;
    auto gap = local_deposition_sections({lower, upper}, {0., 0.}, {.01, .01}, 0., .0001);
    REQUIRE(gap.tiles.size() == 1);
    REQUIRE(gap.tiles.front().slabs.size() == 2);
    CHECK_THAT(gap.tiles.front().slabs.front().top_m, Catch::Matchers::WithinAbs(.001, 1e-12));
    CHECK_THAT(gap.tiles.front().slabs.back().bottom_m, Catch::Matchers::WithinAbs(.002, 1e-12));
    CHECK_THROWS_AS(qualify_local_section_topology(gap.tiles, 0.), std::invalid_argument);
    upper.start_m.z() = upper.end_m.z() = .002;
    upper.material = 1;
    const auto materials = local_deposition_sections({lower, upper}, {0., 0.}, {.01, .01}, 0., .0001);
    REQUIRE(materials.tiles.size() == 1);
    REQUIRE(materials.tiles.front().slabs.size() == 2);
    CHECK(materials.tiles.front().slabs.front().material == 0);
    CHECK(materials.tiles.front().slabs.back().material == 1);
    CHECK(qualify_local_section_topology(materials.tiles, 0.).components() == 1);
}

TEST_CASE("Compaction retains common intervals when one material changes shape", "[ThermoelasticLocalSections]")
{
    const std::vector<DepositedSegment> beads{
        {{.002, .003, .001}, {.004, .003, .001}, .001, .001, 0, BeadSection::RectangularEnvelope},
        {{.002, .003, .002}, {.004, .003, .002}, .002, .001, 0, BeadSection::RectangularEnvelope},
        {{.006, .003, .002}, {.008, .003, .002}, .001, .002, 1, BeadSection::RectangularEnvelope}};
    const auto result = local_deposition_sections(beads, {0., 0.}, {.01, .01}, 0., .0001);
    REQUIRE(result.tiles.size() == 1);
    CHECK(result.tiles.front().slabs.size() == 4);
    CHECK(qualify_local_section_topology(result.tiles, 0.).components() == 2);
}

TEST_CASE("Local tapered geometry converges without flattening the ramp", "[ThermoelasticLocalSections]")
{
    const double length = .004, width = .001, height = .001;
    const std::vector<DepositedSegment> beads{
        {{0., 0., 0.}, {length, 0., height}, width, 0., 0, BeadSection::RectangularEnvelope, height}};
    const double expected = length * width * height / 2. + std::acos(-1.) * width * width * height / 8.;
    const auto result = local_deposition_sections(beads, {0., 0.}, {.001, .001}, 0., .00005);
    CHECK_FALSE(result.piecewise_constant);
    CHECK_THAT(result.volume_m3, Catch::Matchers::WithinRel(expected, .005));
    double left = 0., right = 0.;
    for (const auto& tile : result.tiles)
        for (const auto& slab : tile.slabs)
            for (const auto& polygon : slab.footprint.native()) {
                const double volume = polygon.area() * (slab.top_m - slab.bottom_m) * 1e-18;
                if (tile.index[0] == 0) left += volume;
                if (tile.index[0] == 3) right += volume;
            }
    CHECK(right > left);
}

TEST_CASE("Local refinement splits native intervals shorter than the maximum step", "[ThermoelasticLocalSections]")
{
    const std::vector<DepositedSegment> beads{
        {{0., 0., .001}, {.002, 0., .00101}, .001, .001, 0, BeadSection::RoundedRectangle}};
    const auto coarse = local_deposition_sections(beads, {0., 0.}, {.01, .01}, 0., .1);
    const auto refined = local_deposition_sections(beads, {0., 0.}, {.01, .01}, 0., .1, {}, {}, {}, 2);
    CHECK(refined.sections == 2 * coarse.sections);
    CHECK(refined.segment_visits > coarse.segment_visits);
}

TEST_CASE("Local section clipping rejects only conflicts inside the tile", "[ThermoelasticLocalSections]")
{
    // Conservative bead AABBs overlap many tiles, but conflict detection must
    // use the actually clipped nominal geometry, not authorize/discard material.
    const std::vector<DepositedSegment> beads{
        {{0., 0., .001}, {.004, 0., .001}, .001, .001, 0},
        {{0., 0., .001}, {.004, 0., .001}, .001, .001, 1}};
    CHECK_THROWS_AS(local_deposition_sections(beads, {0., 0.}, {.001, .001}, 0., .001), LayerMaterialOverlap);
}

TEST_CASE("Local section budgets and cancellation never return partial geometry", "[ThermoelasticLocalSections]")
{
    const std::vector<DepositedSegment> beads{{{0., 0., .001}, {.004, 0., .001}, .001, .001, 0}};
    auto limits = LocalSectionLimits{};
    limits.max_tile_references = 1;
    CHECK_THROWS_AS(local_deposition_sections(beads, {0., 0.}, {.001, .001}, 0., .001, limits), ResourceLimit);
    limits = {};
    limits.max_boundary_points = 1;
    CHECK_THROWS_AS(local_deposition_sections(beads, {0., 0.}, {.001, .001}, 0., .001, limits), ResourceLimit);
    CHECK_THROWS_AS(local_deposition_sections(beads, {0., 0.}, {.001, .001}, 0., .001, {}, {},
        [] { throw std::logic_error("cancel"); }), std::logic_error);
}

TEST_CASE("Local ownership accounting releases unused per-plane reservations", "[ThermoelasticLocalSections]")
{
    const std::vector<DepositedSegment> beads{{{0., 0., .001}, {.002, 0., .001}, .001, .001, 0}};
    LocalSectionLimits limits;
    limits.max_ownership_work = limits.geometry.max_slab_group_visits + 100000;
    unsigned calls = 0;
    const auto result = local_deposition_sections(beads, {0., 0.}, {.001, .001}, 0., .0001, limits,
        [&](double z) { ++calls; return SectionOwnership{z, {}}; });
    CHECK(calls > 1);
    CHECK(result.ownership_work > 0);
    CHECK(result.ownership_work < 100000);
    CHECK(result.sampler_work.exact_construction==0);
    CHECK(result.sampler_work.exact_composition==0);
    CHECK(result.sampler_work.other==result.ownership_work);
}

TEST_CASE("Local exact sampler work diagnostics retain category totals across reservations", "[ThermoelasticLocalSections]")
{
    const DepositedSegment bead{{0.,0.,.0002},{.001,0.,.0002},.0004,.0002,0};
    const std::vector<DepositedSegment> beads{bead,bead};
    LocalSectionLimits limits;
    limits.geometry.retain_exact_interface_sections=true;
    limits.geometry.ownership_boundary_grid_units=4;
    limits.geometry.exact_per_command_composition=true;
    const auto result=local_deposition_sections(beads,{-.001,-.001},{.004,.004},0.,.0001,limits);
    CHECK(result.sampler_work.exact_construction>0);
    CHECK(result.sampler_work.exact_composition>0);
    CHECK(result.sampler_work.other>0);
    CHECK(result.sampler_work.exact_construction+result.sampler_work.exact_composition+
        result.sampler_work.other==result.ownership_work);
    const ExPolygon clip(Polygon{{-1000000,-1000000},{3000000,-1000000},
        {3000000,3000000},{-1000000,3000000}});
    const auto one=section_deposition(beads,.00005,0.,limits.geometry,{},nullptr,&clip);
    CHECK(result.sections==2);
    CHECK(result.sampler_work.exact_construction==2*one.sampler_work.exact_construction);
    CHECK(result.sampler_work.exact_composition==2*one.sampler_work.exact_composition);
    CHECK(result.sampler_work.other==2*one.sampler_work.other);
    limits.max_ownership_work=limits.geometry.max_slab_group_visits+one.ownership_work-1;
    try {
        local_deposition_sections(beads,{-.001,-.001},{.004,.004},0.,.0001,limits);
        FAIL("Expected the second reservation to exceed the unchanged budget");
    } catch (const ResourceLimit& error) {
        const std::string message=error.what();
        CHECK(message.find("construction="+std::to_string(one.sampler_work.exact_construction))!=std::string::npos);
        CHECK(message.find("composition="+std::to_string(one.sampler_work.exact_composition))!=std::string::npos);
        CHECK(message.find("other_including_reservation="+std::to_string(one.sampler_work.other))!=std::string::npos);
        CHECK(message.find("total="+std::to_string(one.ownership_work))!=std::string::npos);
    }
}

TEST_CASE("Staggered local slabs connect through physical lateral faces", "[ThermoelasticLocalSections]")
{
    auto left = box_tile(0, 0., .001);
    auto right = box_tile(1, .0005, .0015);
    const auto topology = qualify_local_section_topology({left, right}, 0.);
    CHECK(topology.components() == 1);
    REQUIRE(topology.pieces().size() == 2);
    auto material = elastic_material();
    material.eigenstrain[2] = -.001;
    const auto field = solve_deposited_sections(topology, {material}, {0., 0., 0.}, {.001, .001, .0005}, {}, {},
        {{1, {.0002, .0002, .00075}, 0}, {2, {.0012, .0002, .00075}, 0}});
    CHECK_THAT(field.occupied_volume_m3, Catch::Matchers::WithinRel(2e-9, 1e-10));
    CHECK(field.samples.size() == 2);
    CHECK_THAT(field.elastic_energy_j, Catch::Matchers::WithinAbs(0., 1e-20));
    for (const auto& sample : field.samples) {
        CHECK_THAT(sample.displacement_m.z(), Catch::Matchers::WithinAbs(-.001 * sample.position_m.z(), 1e-14));
        CHECK_THAT(sample.stress_pa.norm(), Catch::Matchers::WithinAbs(0., 1e-8));
    }
    CHECK_THROWS_AS(qualify_local_section_topology({right}, 0.), std::invalid_argument);
}

TEST_CASE("Tile edge and point contacts do not anchor detached local material", "[ThermoelasticLocalSections]")
{
    const auto left = box_tile(0, 0., .001);
    const auto right = box_tile(1, .001, .002);
    CHECK_THROWS_AS(qualify_local_section_topology({left, right}, 0.), std::invalid_argument);
    auto separate = box_tile(2, 0., .001);
    const auto topology = qualify_local_section_topology({left, separate}, 0.);
    CHECK(topology.components() == 2);
    CHECK_THROWS_AS(qualify_local_section_topology({left, left}, 0.), std::invalid_argument);
}

TEST_CASE("Local quadrature feeds an independently refined mechanical study", "[ThermoelasticLocalSections]")
{
    const std::vector<DepositedSegment> beads{{{0., 0., .001}, {.002, 0., .001}, .001, .001, 0}};
    SectionStudySettings settings;
    settings.max_section_height_m = .001;
    settings.relative_tolerance = .05;
    settings.local_tile_size_m = Eigen::Vector2d(.001, .001);
    const auto result = study_deposition_sections(beads, {elastic_material()}, {-.001, -.001, 0.},
        {.001, .0005, .001}, 0., settings);
    CHECK(result.state == DepositionSectionStudy::State::NominalGeometryAndBasisConverged);
    CHECK(result.piecewise_constant_geometry);
    REQUIRE(result.field);
    CHECK_THAT(result.field->occupied_volume_m3, Catch::Matchers::WithinRel(
        (.002 * .001 + std::acos(-1.) * .0005 * .0005) * .001, .005));
    CHECK_THROWS_AS(study_deposition_sections(beads, {elastic_material()}, {-.001, -.001, 0.},
        {.001, .001, .001}, 0., settings, {}, {}, [] { throw ResourceLimit("cancel"); }), ResourceLimit);
}

TEST_CASE("Loaded local and global prism integration produce the same mechanical response", "[ThermoelasticLocalSections]")
{
    const auto left = box_tile(0, 0., .001), right = box_tile(1, 0., .001);
    auto material = elastic_material();
    material.eigenstrain.head<3>().setConstant(-.001);
    const auto local = qualify_local_section_topology({left, right}, 0.);
    auto slabs = left.slabs;
    slabs.insert(slabs.end(), right.slabs.begin(), right.slabs.end());
    const auto global = qualify_section_topology(slabs, 0.);
    const auto a = solve_deposited_sections(local, {material}, {0., 0., 0.}, {.0005, .0005, .0005});
    const auto b = solve_deposited_sections(global, {material}, {0., 0., 0.}, {.0005, .0005, .0005});
    REQUIRE(a.elastic_energy_j > 0.);
    CHECK_THAT(a.elastic_energy_j, Catch::Matchers::WithinRel(b.elastic_energy_j, 1e-9));
    REQUIRE(a.nodes_m.size() == b.nodes_m.size());
    double reaction_magnitude = 0.;
    for (std::size_t i = 0; i < a.nodes_m.size(); ++i) {
        CHECK_THAT((a.displacement_m[i] - b.displacement_m[i]).norm(), Catch::Matchers::WithinAbs(0., 1e-14));
        CHECK_THAT((a.generalized_bed_reaction_n[i] - b.generalized_bed_reaction_n[i]).norm(),
            Catch::Matchers::WithinAbs(0., 1e-14));
        reaction_magnitude += a.generalized_bed_reaction_n[i].norm();
    }
    CHECK(reaction_magnitude > 0.);
}

TEST_CASE("Loaded material interfaces are invariant under local tile decomposition", "[ThermoelasticLocalSections]")
{
    const std::vector<DepositedSegment> beads{
        {{0., 0., .001}, {.01, 0., .001}, .004, .001, 0},
        {{0., 0., .001}, {.01, 0., .001}, .004, .001, 1}};
    const ExPolygon left(Polygon{Point(-100000000, -100000000), Point(-45000000, -100000000),
        Point(55000000, 100000000), Point(-100000000, 100000000)});
    const ExPolygon right(Polygon{Point(-45000000, -100000000), Point(100000000, -100000000),
        Point(100000000, 100000000), Point(55000000, 100000000)});
    const auto ownership = [&](double z) -> std::optional<SectionOwnership> {
        return SectionOwnership{z, {{0, {left}}, {1, {right}}}};
    };
    auto material = elastic_material();
    material.eigenstrain.head<3>().setConstant(-.001);
    const auto coarse = local_deposition_sections(beads, {0., 0.}, {.005, .005}, 0., .001, {}, ownership);
    const auto fine = local_deposition_sections(beads, {0., 0.}, {.0025, .0025}, 0., .001, {}, ownership);
    const auto coarse_topology = qualify_local_section_topology(coarse.tiles, 0.);
    const auto fine_topology = qualify_local_section_topology(fine.tiles, 0.);
    REQUIRE(coarse_topology.components() == 1);
    REQUIRE(fine_topology.components() == 1);
    const auto a = solve_deposited_sections(coarse_topology, {material, material},
        {-.005, -.005, 0.}, {.0025, .0025, .001});
    const auto b = solve_deposited_sections(fine_topology, {material, material},
        {-.005, -.005, 0.}, {.0025, .0025, .001});
    REQUIRE(a.elastic_energy_j > 0.);
    CHECK_THAT(a.occupied_volume_m3, Catch::Matchers::WithinRel(b.occupied_volume_m3, 1e-9));
    CHECK_THAT(a.elastic_energy_j, Catch::Matchers::WithinRel(b.elastic_energy_j, 1e-8));
}

TEST_CASE("Local topology keeps interior holes empty in mechanical integration", "[ThermoelasticLocalSections]")
{
    auto tile = box_tile(0, 0., .001);
    Polygon hole{Point(200000, 200000), Point(200000, 800000), Point(800000, 800000), Point(800000, 200000)};
    tile.slabs.front().footprint.native().front().holes.push_back(hole);
    const auto topology = qualify_local_section_topology({tile}, 0.);
    const auto field = solve_deposited_sections(topology, {elastic_material()}, {0., 0., 0.}, {.001, .001, .001}, {}, {},
        {{1, {.0001, .0001, .0005}, 0}, {2, {.0005, .0005, .0005}, 0}});
    CHECK_THAT(field.occupied_volume_m3, Catch::Matchers::WithinRel(.64e-9, 1e-10));
    CHECK(field.unsampled_points == 1);
    REQUIRE(field.samples.size() == 1);
    CHECK(field.samples.front().id == 1);
}

TEST_CASE("Explicit Z approximation retains separate qualification and refinement metadata", "[ThermoelasticLocalSections]")
{
    const std::vector<DepositedSegment> beads{{{0., 0., .000999999}, {.002, 0., .000999999}, .001, .000999999, 0}};
    SectionStudySettings settings;
    settings.max_section_height_m = .001;
    settings.relative_tolerance = .05;
    settings.max_total_work = 300000000;
    settings.local_tile_size_m = Eigen::Vector2d(.001, .001);
    settings.initial_z_quantization_m = 1e-6;
    const auto result = study_deposition_sections(beads, {elastic_material()}, {-.001, -.001, 0.},
        {.001, .0005, .001}, 0., settings);
    CHECK(result.state == DepositionSectionStudy::State::ApproximateGeometryAndBasisConverged);
    CHECK_FALSE(result.piecewise_constant_geometry);
    REQUIRE(result.field);
    REQUIRE(result.history.size() == 3);
    CHECK_THAT(result.history.back().z_quantization_m, Catch::Matchers::WithinRel(.25e-6, 1e-12));
    CHECK_THAT(result.history.back().maximum_partition_shift_m, Catch::Matchers::WithinAbs(1e-9, 1e-15));
    CHECK_THAT(beads.front().start_m.z(), Catch::Matchers::WithinAbs(.000999999, 1e-15));

    auto lifted = beads;
    for (auto& bead : lifted) { bead.start_m.z() += 1e-9; bead.end_m.z() += 1e-9; }
    CHECK_THROWS_AS(study_deposition_sections(lifted, {elastic_material()}, {-.001, -.001, 0.},
        {.001, .0005, .001}, 0., settings), std::invalid_argument);

    auto contracting = elastic_material();
    contracting.eigenstrain.head<3>().setConstant(-.001);
    settings.max_geometry_levels = 1;
    const auto loaded = study_deposition_sections(beads, {contracting}, {-.001, -.001, 0.},
        {.001, .0005, .001}, 0., settings);
    REQUIRE_FALSE(loaded.last_basis_history.empty());
    CHECK(loaded.last_basis_history.front().elastic_energy_j > 0.);
    CHECK_FALSE(loaded.field); // A loaded result cannot bypass geometry refinement.
}

TEST_CASE("Approximate section partitions sample original rounded material geometry", "[ThermoelasticLocalSections]")
{
    std::vector<DepositedSegment> original{{{.0008, .001, .0009}, {.001, .001, .0009}, .001, .0009, 0}};
    original.front().width_m = .001;
    original.front().section = BeadSection::RoundedRectangle;
    auto partition = original;
    partition.front().start_m.z() = partition.front().end_m.z() = .001;
    partition.front().height_m = .001;
    const auto local = local_deposition_sections(original, {0., 0.}, {.002, .002}, 0., .001,
        {}, {}, {}, 1, nullptr, &partition);
    REQUIRE(local.tiles.size() == 1);
    REQUIRE(local.tiles.front().slabs.size() == 1);
    const auto& slab = local.tiles.front().slabs.front();
    const double z = .5 * (slab.bottom_m + slab.top_m);
    const auto expected = section_deposition(original, z, 0.);
    const auto altered = section_deposition(partition, z, 0.);
    CHECK(diff_ex(slab.footprint.native(), expected.regions.front().footprint.native()).empty());
    CHECK(diff_ex(expected.regions.front().footprint.native(), slab.footprint.native()).empty());
    CHECK_FALSE(diff_ex(altered.regions.front().footprint.native(), slab.footprint.native()).empty());
    CHECK(slab.material == original.front().material);

    auto narrow = original;
    narrow.front().height_m = 1e-7;
    CHECK_THROWS_AS(local_deposition_sections(narrow, {0., 0.}, {.002, .002}, 0., .001,
        {}, {}, {}, 1, nullptr, &partition), std::invalid_argument); // New bed anchoring is forbidden.
    partition.front().height_m = .0002;
    const auto retained = local_deposition_sections(narrow, {0., 0.}, {.002, .002}, 0., .001,
        {}, {}, {}, 1, nullptr, &partition);
    REQUIRE(retained.tiles.size() == 1);
    REQUIRE(retained.tiles.front().slabs.size() == 1);
    CHECK_THAT(retained.tiles.front().slabs.front().bottom_m, Catch::Matchers::WithinAbs(.0008999, 1e-12));
    CHECK_THAT(retained.tiles.front().slabs.front().top_m, Catch::Matchers::WithinAbs(.0009, 1e-12));
    LocalSectionLimits bounded;
    bounded.max_sections = 1;
    CHECK_THROWS_AS(local_deposition_sections(narrow, {0., 0.}, {.002, .002}, 0., .001,
        bounded, {}, {}, 1, nullptr, &partition), ResourceLimit);
    const auto resolved = local_deposition_sections(narrow, {0., 0.}, {.002, .002}, 0., 1e-7,
        {}, {}, {}, 1, nullptr, &partition);
    CHECK(resolved.volume_m3 > 0.);
}

TEST_CASE("Retained thin features replace adjacent artificial grid cuts", "[ThermoelasticLocalSections]")
{
    const DepositedSegment thin{{.001, .001, 999999e-9}, {.002, .001, 999999e-9},
        .0001, 99e-9, 0, BeadSection::RectangularEnvelope};
    const DepositedSegment carrier{{.004, .001, .002}, {.005, .001, .002},
        .0001, .002, 0, BeadSection::RectangularEnvelope};
    const std::vector<DepositedSegment> original{thin, carrier};
    auto partition = original;
    partition.front().start_m.z() = partition.front().end_m.z() = .001;
    partition.front().height_m = 1e-6;
    const auto local = local_deposition_sections(original, {0., 0.}, {.01, .01}, 0., .002,
        {}, {}, {}, 1, nullptr, &partition);
    const auto exact = layered_deposition(original, 0.);
    CHECK_THAT(local.volume_m3, Catch::Matchers::WithinRel(exact.occupied_volume_m3, 1e-10));
    REQUIRE(local.tiles.size() == 1);
    bool retained = false;
    for (const auto& slab : local.tiles.front().slabs)
        if (std::abs(slab.top_m - 999999e-9) < 1e-15 && std::abs(slab.bottom_m - 999900e-9) < 1e-15) {
            CHECK(slab.footprint.native().size() == 2);
            retained = true;
        }
    CHECK(retained);
    // Without an explicit approximate partition, original boundaries remain
    // mandatory. A true one-nanometre interval must still refuse sampling.
    auto adjacent = thin;
    adjacent.start_m.x() += .006;
    adjacent.end_m.x() += .006;
    adjacent.start_m.z() = adjacent.end_m.z() = .001;
    CHECK_THROWS_AS(local_deposition_sections({thin, adjacent}, {0., 0.}, {.01, .01}, 0., .002), ResourceLimit);
    LocalSectionLimits half_planes;
    half_planes.geometry.half_native_section_planes = true;
    const auto resolved = local_deposition_sections({thin, adjacent}, {0., 0.}, {.01, .01}, 0., .002, half_planes);
    const auto reference = layered_deposition({thin, adjacent}, 0.);
    CHECK_THAT(resolved.volume_m3, Catch::Matchers::WithinRel(reference.occupied_volume_m3, 1e-10));
    REQUIRE(resolved.tiles.size() == 1);
    bool thin_interval = false;
    for (const auto& slab : resolved.tiles.front().slabs)
        if (std::abs(slab.bottom_m - 999999e-9) < 1e-15) {
            CHECK_THAT(slab.top_m, Catch::Matchers::WithinAbs(.001, 1e-15));
            thin_interval = true;
        }
    CHECK(thin_interval);
    auto native_thin = thin;
    native_thin.height_m = 1e-9;
    const auto native_result = local_deposition_sections({native_thin}, {0., 0.}, {.01, .01}, 0., 1e-9, half_planes);
    REQUIRE(native_result.tiles.size() == 1);
    REQUIRE(native_result.tiles.front().slabs.size() == 1);
    CHECK(native_result.volume_m3 > 0.);
}

TEST_CASE("Local partition refinement rechecks displaced midplanes and duplicate bead intervals", "[ThermoelasticLocalSections]")
{
    const DepositedSegment upper{{.001, .001, .001}, {.002, .001, .001}, .001, 100e-9, 0,
        BeadSection::RectangularEnvelope};
    auto middle = upper;
    middle.start_m.z() = middle.end_m.z() = 999520e-9;
    middle.height_m = 40e-9;
    const std::vector<DepositedSegment> original{upper, upper, middle};
    auto envelopes = original;
    for (auto& bead : envelopes) {
        bead.start_m.z() = bead.end_m.z() = .001;
        bead.height_m = 1e-6;
    }
    // The first coarse midpoint hits only the middle bead. Retaining the upper
    // bead moves that midpoint outside the middle bead, requiring a second pass.
    const auto result = local_deposition_sections(original, {0., 0.}, {.004, .004}, 0., .001,
        {}, {}, {}, 1, nullptr, &envelopes);
    CHECK(result.sections == 4);
    REQUIRE(result.tiles.size() == 1);
    REQUIRE(result.tiles.front().slabs.size() == 2);
    CHECK_THAT(result.tiles.front().slabs.front().bottom_m, Catch::Matchers::WithinAbs(999480e-9, 1e-12));
    CHECK_THAT(result.tiles.front().slabs.front().top_m, Catch::Matchers::WithinAbs(999520e-9, 1e-12));
    CHECK_THAT(result.tiles.front().slabs.back().bottom_m, Catch::Matchers::WithinAbs(999900e-9, 1e-12));
    CHECK_THAT(result.tiles.front().slabs.back().top_m, Catch::Matchers::WithinAbs(.001, 1e-12));
}
