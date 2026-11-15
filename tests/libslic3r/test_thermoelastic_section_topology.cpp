#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/ThermoelasticLayers.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticLocalSections.hpp"
#include "libslic3r/ThermoelasticLocalTopology.hpp"
#include "libslic3r/ThermoelasticSections.hpp"
#include "libslic3r/ThermoelasticSectionTopology.hpp"
#include "libslic3r/ThermoelasticExactSections.hpp"
#include "libslic3r/ThermoelasticSectionFootprint.hpp"
#include "libslic3r/libslic3r.h"
#include <stdexcept>
#include <string>
#include <cstddef>
#include <utility>
#include <array>
#include <map>
#include <Eigen/Core>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;
namespace {
ExPolygon box(coord_t x, coord_t y, coord_t size = 1000000)
{
    return ExPolygon(Polygon{{x,y},{x+size,y},{x+size,y+size},{x,y+size}});
}

TEST_CASE("Retained topology preserves rational components and positive vertical contact", "[ThermoelasticSectionTopology]")
{
    const ExPolygons triangle{ExPolygon(Polygon{{0,0},{10,0},{10,4}})};
    const auto full = make_exact_section(triangle).section;
    const auto strip = make_exact_section({ExPolygon(Polygon{{0,0},{1,0},{1,10},{0,10}})}).section;
    const auto small = boolean_exact_sections(full,strip,ExactSectionOperation::Intersection).section;
    const auto topology = qualify_section_topology({{0.,.001,0,triangle},
        {.001,.002,0,SectionFootprint(small)}},0.);
    REQUIRE(topology.components() == 1);
    REQUIRE(topology.pieces().size() == 2);
    for (const auto& piece : topology.pieces()) {
        REQUIRE(piece.slab.footprint.is_exact());
        CHECK(piece.component == 0);
        const auto& expected = piece.slab.bottom_m == 0. ? full : small;
        CHECK(boolean_exact_sections(piece.slab.footprint.exact(),expected,
            ExactSectionOperation::SymmetricDifference).section.empty());
    }
    // Repeated same-material contributions are a union, not duplicate mass.
    const auto duplicate = qualify_section_topology({{0.,.001,0,SectionFootprint(small)},
        {0.,.001,0,SectionFootprint(small)}},0.);
    REQUIRE(duplicate.pieces().size() == 1);
    CHECK(boolean_exact_sections(duplicate.pieces().front().slab.footprint.exact(),small,
        ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK_THROWS_AS(qualify_section_topology({{0.,.001,0,SectionFootprint(full)},
        {0.,.001,1,SectionFootprint(small)}},0.),LayerMaterialOverlap);
    CHECK_THROWS_WITH(qualify_section_topology({{0.,.001,0,SectionFootprint(small),1}},0.),
        "Boundary uncertainty requires local topology contact qualification.");
    SectionTopologyLimits limited;
    limited.max_component_pairs = 1;
    CHECK_THROWS_AS(qualify_section_topology({{0.,.001,0,SectionFootprint(small)}},0.,limited),ResourceLimit);
    limited = {};
    limited.max_boundary_points = 2;
    CHECK_THROWS_AS(qualify_section_topology({{0.,.001,0,SectionFootprint(small)}},0.,limited),ResourceLimit);
    CHECK_THROWS_AS(qualify_section_topology({{0.,.001,0,SectionFootprint(small)}},0.,{},
        [] { throw std::runtime_error("cancel exact topology"); }),std::runtime_error);
}

TEST_CASE("Retained exact topology rejects point-only joins and pinched components", "[ThermoelasticSectionTopology]")
{
    const auto lower = make_exact_section({box(0,0,10)}).section;
    const auto upper = make_exact_section({box(10,10,10)}).section;
    CHECK_THROWS_AS(qualify_section_topology({{0.,.001,0,SectionFootprint(lower)},
        {0.,.001,0,SectionFootprint(upper)}},0.),std::invalid_argument);
    CHECK_THROWS_WITH(qualify_section_topology({{0.,.001,0,SectionFootprint(lower)},
        {.001,.002,0,SectionFootprint(upper)}},0.),
        "Section topology contains an unanchored solid component.");
    const auto edge_neighbor = make_exact_section({box(10,0,10)}).section;
    CHECK_THROWS_WITH(qualify_section_topology({{0.,.001,0,SectionFootprint(lower)},
        {.001,.002,0,SectionFootprint(edge_neighbor)}},0.),
        "Section topology contains an unanchored solid component.");
    const auto touching_hole = make_exact_section({ExPolygon(Polygon{{0,5},{2,4},{2,6}})}).section;
    const auto pinched = boolean_exact_sections(lower,touching_hole,ExactSectionOperation::Difference).section;
    CHECK_THROWS_AS(qualify_section_topology({{0.,.001,0,SectionFootprint(pinched)}},0.),std::invalid_argument);
}

TEST_CASE("Explicit exact interior topology separates point contacts without filling tangent voids", "[ThermoelasticSectionTopology]")
{
    const unsigned shape=GENERATE(0u,1u,2u);
    auto original=make_exact_section({box(0,0,10)}).section;
    if (shape==0)
        original=boolean_exact_sections(original,make_exact_section({box(10,10,10)}).section,
            ExactSectionOperation::Union).section;
    else {
        const auto subtract=[&](const Polygon& hole) {
            original=boolean_exact_sections(original,make_exact_section({ExPolygon(hole)}).section,
                ExactSectionOperation::Difference).section;
        };
        if (shape==1) {
            subtract(Polygon{{1,5},{3,3},{5,5},{3,7}});
            subtract(Polygon{{5,5},{7,3},{9,5},{7,7}});
        } else {
            subtract(Polygon{{0,5},{2,3},{5,5},{2,7}});
            subtract(Polygon{{5,5},{8,3},{10,5},{8,7}});
        }
    }
    const std::vector<DepositedSlab> slabs{{0.,.001,0,SectionFootprint(original)}};
    const LocalSectionTile tile{{0,0},{0,0},{20,20},slabs};
    SectionTopologyLimits limits;
    CHECK_FALSE(limits.exact_interior_components);
    CHECK_THROWS_AS(qualify_section_topology(slabs,0.,limits),std::invalid_argument);
    CHECK_THROWS_AS(qualify_local_section_topology({tile},0.,limits),std::invalid_argument);
    limits.exact_interior_components=true;
    const auto check=[&](const SectionTopology& topology) {
        const std::size_t expected=shape==1 ? 1 : 2;
        REQUIRE(topology.components()==expected);
        REQUIRE(topology.pieces().size()==expected);
        ExactSection reconstructed;
        for (const auto& piece:topology.pieces())
            reconstructed=boolean_exact_sections(reconstructed,piece.slab.footprint.exact(),ExactSectionOperation::Union).section;
        CHECK(boolean_exact_sections(reconstructed,original,ExactSectionOperation::SymmetricDifference).section.empty());
        if (expected==2) CHECK(topology.pieces()[0].component!=topology.pieces()[1].component);
    };
    check(qualify_section_topology(slabs,0.,limits));
    check(qualify_local_section_topology({tile},0.,limits));
    auto bounded=limits;
    bounded.max_component_pairs=1;
    CHECK_THROWS_AS(qualify_local_section_topology({tile},0.,bounded),ResourceLimit);
    CHECK_THROWS_AS(qualify_local_section_topology({tile},0.,limits,
        [] { throw std::runtime_error("cancel interior topology"); }),std::runtime_error);
}

TEST_CASE("Explicit exact interior topology cannot anchor an upper solid through a point", "[ThermoelasticSectionTopology]")
{
    SectionTopologyLimits limits;
    limits.exact_interior_components=true;
    const auto lower=make_exact_section({box(0,0,10)}).section;
    const auto upper=make_exact_section({box(10,10,10)}).section;
    const std::vector<DepositedSlab> slabs{{0.,.001,0,SectionFootprint(lower)},
        {.001,.002,0,SectionFootprint(upper)}};
    CHECK_THROWS_AS(qualify_section_topology(slabs,0.,limits),std::invalid_argument);
    CHECK_THROWS_AS(qualify_local_section_topology({{{0,0},{0,0},{20,20},slabs}},0.,limits),std::invalid_argument);
    const std::vector<LocalSectionTile> neighbors{
        {{0,0},{0,0},{10,20},{{0.,.001,0,SectionFootprint(lower)}}},
        {{1,0},{10,0},{20,20},{{0.,.001,0,SectionFootprint(upper)}}}};
    CHECK(qualify_local_section_topology(neighbors,0.,limits).components()==2);
}

TEST_CASE("Retained exact topology keeps holes and assigns adjacent materials to one solid", "[ThermoelasticSectionTopology]")
{
    ExPolygon holed = box(0,0,10);
    holed.holes.emplace_back(Polygon{{2,2},{2,4},{4,4},{4,2}});
    const auto left = make_exact_section({holed}).section;
    const auto topology = qualify_section_topology({{0.,.001,0,SectionFootprint(left)},
        {0.,.001,1,{box(10,0,10)}}},0.);
    CHECK(topology.components() == 1);
    REQUIRE(topology.pieces().size() == 2);
    for (const auto& piece : topology.pieces()) {
        REQUIRE(piece.slab.footprint.is_exact());
        CHECK(piece.component == 0);
        if (piece.slab.material == 0)
            CHECK(boolean_exact_sections(piece.slab.footprint.exact(),left,
                ExactSectionOperation::SymmetricDifference).section.empty());
    }
}

TEST_CASE("Ordered topology releases absent neighbors and recycles a bounded face frontier", "[ThermoelasticSectionTopology]")
{
    StreamTopologyLimits limits;
    limits.max_islands = 100;
    StreamTopologyStatistics stats;
    LocalTopologyBuilder ordered(0., {}, [](const SectionPiece&) {}, {}, limits, &stats, LocalTileOrder::Lexicographic);
    LocalTopologyBuilder arbitrary(0., {}, [](const SectionPiece&) {}, {}, limits);
    for (int x = -20; x < 30; ++x) {
        const coord_t lo = coord_t(x) * 1000000;
        const LocalSectionTile tile{{x,0},{lo,0},{lo+1000000,1000000},{{0.,.001,0,{box(lo,0)}}}};
        ordered.add_tile(tile);
        arbitrary.add_tile(tile);
        CHECK(stats.face_records <= 1);
    }
    const auto actual = ordered.finish();
    const auto reference = arbitrary.finish();
    CHECK(actual.components() == reference.components());
    CHECK(actual.components() == 1);
    CHECK(stats.peak_face_records <= 2);
    CHECK(stats.face_records == 0);
    CHECK(stats.trace_intervals == 0);
}

TEST_CASE("Translated native slivers retain exact positive topology area", "[ThermoelasticSectionTopology]")
{
    for (const coord_t offset : {coord_t(0),coord_t(1000000000)}) {
        const ExPolygon triangle(Polygon{{offset,offset},{offset+3,offset},{offset,offset+3}});
        const LocalSectionTile tile{{0,0},{offset,offset},{offset+10,offset+10},{{0.,.001,0,{triangle}}}};
        std::size_t emitted=0;
        LocalTopologyBuilder builder(0.,{},[&](const SectionPiece& piece) {
            REQUIRE(piece.slab.footprint.native().size()==1);
            CHECK(piece.slab.footprint.native()[0].contour.points.size()==3);
            ++emitted;
        });
        builder.add_tile(tile);
        CHECK(builder.finish().components()==1);
        CHECK(emitted==1);
    }
}

TEST_CASE("Translated native holes retain their negative winding", "[ThermoelasticSectionTopology]")
{
    const coord_t offset=1000000000;
    auto polygon=box(offset,offset,100);
    polygon.holes.push_back(Polygon{{offset+2,offset+2},{offset+2,offset+5},{offset+5,offset+2}});
    const LocalSectionTile tile{{0,0},{offset,offset},{offset+100,offset+100},{{0.,.001,0,{polygon}}}};
    std::size_t emitted=0;
    LocalTopologyBuilder builder(0.,{},[&](const SectionPiece& piece) {
        REQUIRE(piece.slab.footprint.native().size()==1);
        CHECK(piece.slab.footprint.native()[0].holes.size()==1);
        ++emitted;
    });
    builder.add_tile(tile);
    CHECK(builder.finish().components()==1);
    CHECK(emitted==1);
}

TEST_CASE("Ordered topology preserves two-axis joins missing tiles and vertical gaps", "[ThermoelasticSectionTopology]")
{
    // A floating branch is anchored across Y, then across X. The final island
    // is separated by a missing tile and must remain a distinct component.
    const std::vector<LocalSectionTile> tiles{
        {{0,0},{0,0},{1000000,1000000},{{0.,.002,0,{box(0,0)}}}},
        {{0,1},{0,1000000},{1000000,2000000},{{.001,.002,0,{box(0,1000000)}}}},
        {{1,1},{1000000,1000000},{2000000,2000000},{{.001,.002,0,{box(1000000,1000000)}}}},
        {{1,3},{1000000,3000000},{2000000,4000000},{{0.,.001,0,{box(1000000,3000000)}}}}
    };
    for (const auto order : {LocalTileOrder::Arbitrary, LocalTileOrder::Lexicographic}) {
        LocalTopologyBuilder builder(0., {}, [](const SectionPiece&) {}, {}, {}, nullptr, order);
        for (const auto& tile : tiles) builder.add_tile(tile);
        CHECK(builder.finish().components() == 2);
        auto disconnected = tiles;
        disconnected[2].slabs[0].bottom_m = .002;
        disconnected[2].slabs[0].top_m = .003;
        LocalTopologyBuilder rejected(0., {}, [](const SectionPiece&) {}, {}, {}, nullptr, order);
        for (const auto& tile : disconnected) rejected.add_tile(tile);
        CHECK_THROWS_AS(rejected.finish(), std::invalid_argument);
    }
}

TEST_CASE("Ordered topology rejects backward traversal and poisons its certificate", "[ThermoelasticSectionTopology]")
{
    const LocalSectionTile a{{0,0},{0,0},{1000000,1000000},{{0.,.001,0,{box(0,0)}}}};
    const LocalSectionTile b{{1,0},{1000000,0},{2000000,1000000},{{0.,.001,0,{box(1000000,0)}}}};
    LocalTopologyBuilder builder(0., {}, [](const SectionPiece&) {}, {}, {}, nullptr, LocalTileOrder::Lexicographic);
    builder.add_tile(b);
    CHECK_THROWS_AS(builder.add_tile(a), std::invalid_argument);
    CHECK_THROWS_AS(builder.finish(), std::invalid_argument);
    CHECK_THROWS_AS(builder.add_tile(b), std::invalid_argument);
}

TEST_CASE("Detached component diagnostics aggregate late joins and identify the largest sampled solid", "[ThermoelasticSectionTopology]")
{
    const std::vector<LocalSectionTile> tiles{
        {{0,0},{0,0},{1000000,1000000},{{.002,.003,0,{box(100000,100000,100000)}}}},
        {{1,0},{1000000,0},{2000000,1000000},{{.001,.003,0,{box(1000000,0)}}}},
        {{2,0},{2000000,0},{3000000,1000000},
            {{.001,.002,0,{box(2000000,0)}},{.002,.003,0,{box(2000000,0)}}}}
    };
    for (const auto order : {LocalTileOrder::Arbitrary, LocalTileOrder::Lexicographic}) {
        LocalTopologyBuilder builder(0., {}, [](const SectionPiece&) {}, {}, {}, nullptr, order);
        for (const auto& tile : tiles) builder.add_tile(tile);
        std::string failure;
        try { builder.finish(); }
        catch (const std::invalid_argument& error) { failure = error.what(); }
        CHECK(failure.find("count=2") != std::string::npos);
        CHECK(failure.find("minimum_nm=1000000,0,1000000") != std::string::npos);
        CHECK(failure.find("maximum_nm=3000000,1000000,3000000") != std::string::npos);
        CHECK(failure.find("pieces=3 labels=2 islands=3") != std::string::npos);
        CHECK(failure.find("bed_gap_nm=1000000") != std::string::npos);
        REQUIRE(failure.find("volume_m3=") != std::string::npos);
        CHECK_THAT(std::stod(failure.substr(failure.find("volume_m3=") + 10)),
            Catch::Matchers::WithinRel(4e-9,1e-12));
        CHECK(failure.size() < 512);
        CHECK_THROWS_AS(builder.finish(),std::invalid_argument);
    }
}

TEST_CASE("Explicit topology streaming separates cumulative visits from live polygon payload", "[ThermoelasticSectionTopology]")
{
    LocalSectionTile first{{0,0},{0,0},{1000000,1000000},{{0.,.001,0,{box(0,0)}}}};
    LocalSectionTile second{{1,0},{1000000,0},{2000000,1000000},{{0.,.001,0,{box(1000000,0)}}}};
    const auto reference = qualify_local_section_topology({first,second}, 0.);
    SectionTopologyLimits old;
    old.max_boundary_points = 4;
    LocalTopologyBuilder legacy(0., old, [](const SectionPiece&) {});
    CHECK_THROWS_AS(legacy.add_tile(first), ResourceLimit);
    CHECK_THROWS_AS(legacy.finish(), std::invalid_argument);
    StreamTopologyLimits limits;
    StreamTopologyStatistics stats;
    std::vector<SectionPiece> emitted;
    LocalTopologyBuilder stream(0., old, [&](const SectionPiece& piece) { emitted.push_back(piece); }, {}, limits, &stats);
    stream.add_tile(first);
    const auto first_peak = stats.peak_live_points;
    stream.add_tile(second);
    const auto certificate = stream.finish();
    CHECK(certificate.components() == reference.components());
    REQUIRE(emitted.size() == reference.pieces().size());
    CHECK(stats.peak_live_points == first_peak);
    CHECK(stats.input_points == 8);
    CHECK(stats.processed_points > stats.input_points);
    CHECK(stats.islands == 2);
    CHECK(stats.trace_intervals == 8);
    for (const auto& piece : emitted)
        CHECK(certificate.component_mapping().at(piece.component) == 0);
    limits.max_live_points = first_peak - 1;
    LocalTopologyBuilder live(0., {}, [](const SectionPiece&) {}, {}, limits);
    CHECK_THROWS_AS(live.add_tile(first), ResourceLimit);
    CHECK_THROWS_AS(live.finish(), std::invalid_argument);
    limits = StreamTopologyLimits{};
    limits.max_processed_points = stats.processed_points - 1;
    LocalTopologyBuilder cumulative(0., {}, [](const SectionPiece&) {}, {}, limits);
    CHECK_NOTHROW(cumulative.add_tile(first));
    CHECK_THROWS_AS(cumulative.add_tile(second), ResourceLimit);
    CHECK_THROWS_AS(cumulative.finish(), std::invalid_argument);
    limits = StreamTopologyLimits{};
    limits.max_live_points = 3;
    StreamTopologyStatistics refused;
    LocalTopologyBuilder borrowed(0., {}, [](const SectionPiece&) {}, {}, limits, &refused);
    CHECK_THROWS_AS(borrowed.add_tile(first), ResourceLimit);
    CHECK(refused.input_points == 0);
    CHECK(refused.islands == 0);
    limits = StreamTopologyLimits{};
    limits.max_processed_points = 128000001;
    CHECK_THROWS_AS(LocalTopologyBuilder(0., {}, [](const SectionPiece&) {}, {}, limits), std::invalid_argument);
    limits = StreamTopologyLimits{};
    limits.max_work = 4000000001ULL;
    CHECK_THROWS_AS(LocalTopologyBuilder(0., {}, [](const SectionPiece&) {}, {}, limits), std::invalid_argument);
    limits = StreamTopologyLimits{};
    limits.max_live_points = 2000001;
    CHECK_THROWS_AS(LocalTopologyBuilder(0., {}, [](const SectionPiece&) {}, {}, limits), std::invalid_argument);
    limits = StreamTopologyLimits{};
    limits.max_slabs = 1;
    LocalTopologyBuilder slab_limited(0., {}, [](const SectionPiece&) {}, {}, limits);
    slab_limited.add_tile(first);
    CHECK_THROWS_AS(slab_limited.add_tile(second), ResourceLimit);
    CHECK_THROWS_AS(slab_limited.finish(), std::invalid_argument);
    limits.max_slabs = 250001;
    CHECK_THROWS_AS(LocalTopologyBuilder(0., {}, [](const SectionPiece&) {}, {}, limits), std::invalid_argument);
    limits.max_slabs = 250000;
    limits.max_islands = 250001;
    CHECK_THROWS_AS(LocalTopologyBuilder(0., {}, [](const SectionPiece&) {}, {}, limits), std::invalid_argument);
    limits.max_islands = 250000;
    limits.max_processed_points = 128000000;
    limits.max_work = 400000000;
    LocalTopologyBuilder offline(0., {}, [](const SectionPiece&) {}, {}, limits);
    offline.add_tile(first);
    offline.add_tile(second);
    CHECK(offline.finish().components() == reference.components());
}
}

TEST_CASE("Section topology retains separate bodies and joins only positive-area layer contacts", "[ThermoelasticSectionTopology]")
{
    const DepositedSlab a{0., .001, 0, {box(0,0)}};
    const DepositedSlab b{0., .001, 1, {box(2000000,0)}};
    const auto separate = qualify_section_topology({a,b}, 0.);
    CHECK(separate.components() == 2);
    REQUIRE(separate.pieces().size() == 2);
    CHECK(separate.pieces()[0].component != separate.pieces()[1].component);
    // A real upper bridge joins both bodies; empty lower space remains empty.
    const ExPolygon bridge(Polygon{{0,0},{3000000,0},{3000000,1000000},{0,1000000}});
    const auto joined = qualify_section_topology({a,b,{.001,.002,0,{bridge}}}, 0.);
    CHECK(joined.components() == 1);
    CHECK(joined.pieces().size() == 3);
    double volume = 0.;
    for (const auto &piece : joined.pieces())
        for (const auto &polygon : piece.slab.footprint.native())
            volume += polygon.area() * 1e-18 * (piece.slab.top_m-piece.slab.bottom_m);
    CHECK_THAT(volume, Catch::Matchers::WithinAbs(5e-9, 1e-20));
    CHECK_THROWS_AS(qualify_section_topology({a,{.001,.002,0,{box(1000000,0)}}},0.), std::invalid_argument);
    CHECK_THROWS_AS(qualify_section_topology({a,{.002,.003,0,{box(0,0)}}},0.), std::invalid_argument);
}

TEST_CASE("Band packets preserve vertical contacts and release completed partition payload", "[ThermoelasticSectionTopology]")
{
    LocalSectionTile tile{{0,0},{0,0},{1000000,1000000},{}};
    for (unsigned i = 0; i < 20; ++i)
        tile.slabs.push_back({double(i) * .001, double(i + 1) * .001, 0, {box(0,0)}});
    StreamTopologyLimits limits;
    StreamTopologyStatistics retained_stats, packet_stats;
    LocalTopologyBuilder retained(0., {}, [](const SectionPiece&) {}, {}, limits, &retained_stats);
    retained.add_tile(tile);
    const auto expected = retained.finish();
    std::vector<SectionPiece> emitted;
    LocalTopologyBuilder packets(0., {}, [&](const SectionPiece& piece) { emitted.push_back(piece); }, {}, limits, &packet_stats);
    LocalSectionTile packet{tile.index,tile.minimum,tile.maximum,{}};
    packets.begin_tile(packet);
    for (const auto& slab : tile.slabs) {
        packet.slabs = {slab};
        packets.add_band(packet);
        packet.slabs.clear();
    }
    packets.end_tile();
    const auto actual = packets.finish();
    CHECK(actual.components() == expected.components());
    CHECK(actual.component_mapping() == expected.component_mapping());
    CHECK(emitted.size() == tile.slabs.size());
    CHECK(packet_stats.input_points == retained_stats.input_points);
    CHECK(packet_stats.processed_points == retained_stats.processed_points);
    CHECK(packet_stats.work == retained_stats.work);
    CHECK(packet_stats.peak_live_points < retained_stats.peak_live_points);
    CHECK(packet_stats.processed_islands == 20);
    CHECK(packet_stats.labels == 1);
    CHECK(packet_stats.face_records == 1);
    CHECK(packet_stats.trace_intervals == 4);
    CHECK(actual.component_mapping().size() == 1);
    limits.max_live_points = packet_stats.peak_live_points;
    LocalTopologyBuilder bounded(0., {}, [](const SectionPiece&) {}, {}, limits);
    CHECK_THROWS_AS(bounded.add_tile(tile), ResourceLimit);
    LocalTopologyBuilder gap(0., {}, [](const SectionPiece&) {});
    gap.begin_tile(packet);
    packet.slabs = {tile.slabs.front()};
    gap.add_band(packet);
    packet.slabs = {tile.slabs[2]};
    gap.add_band(packet);
    gap.end_tile();
    CHECK_THROWS_AS(gap.finish(), std::invalid_argument);
}

TEST_CASE("Incomplete or mismatched band packet sequences poison topology streams", "[ThermoelasticSectionTopology]")
{
    const LocalSectionTile header{{0,0},{0,0},{1000000,1000000},{}};
    LocalSectionTile packet = header;
    packet.slabs = {{0.,.001,0,{box(0,0)}}};
    LocalTopologyBuilder missing(0., {}, [](const SectionPiece&) {});
    CHECK_THROWS_AS(missing.add_band(packet), std::invalid_argument);
    CHECK_THROWS_AS(missing.begin_tile(header), std::invalid_argument);
    LocalTopologyBuilder open(0., {}, [](const SectionPiece&) {});
    open.begin_tile(header);
    open.add_band(packet);
    CHECK_THROWS_AS(open.finish(), std::invalid_argument);
    CHECK_THROWS_AS(open.end_tile(), std::invalid_argument);
    LocalTopologyBuilder mismatch(0., {}, [](const SectionPiece&) {});
    mismatch.begin_tile(header);
    packet.index = {1,0};
    CHECK_THROWS_AS(mismatch.add_band(packet), std::invalid_argument);
    CHECK_THROWS_AS(mismatch.finish(), std::invalid_argument);
    packet.index = header.index;
    packet.slabs.push_back({.001,.002,0,{box(0,0)}});
    LocalTopologyBuilder mixed(0., {}, [](const SectionPiece&) {});
    mixed.begin_tile(header);
    CHECK_THROWS_AS(mixed.add_band(packet), std::invalid_argument);
    CHECK_THROWS_AS(mixed.finish(), std::invalid_argument);
    LocalTopologyBuilder reentry(0., {}, [](const SectionPiece&) {});
    reentry.begin_tile(header);
    CHECK_THROWS_AS(reentry.begin_tile(header), std::invalid_argument);
    CHECK_THROWS_AS(reentry.finish(), std::invalid_argument);
}

TEST_CASE("Merged lateral records preserve real Z gaps during delayed anchoring", "[ThermoelasticSectionTopology]")
{
    const auto stripe = [](coord_t x, coord_t y0, coord_t y1) {
        return ExPolygon(Polygon{{x,y0},{x+1000000,y0},{x+1000000,y1},{x,y1}});
    };
    LocalSectionTile left{{0,0},{0,0},{1000000,1000000},
        {{.002,.003,0,{stripe(0,600000,800000)}}}};
    LocalSectionTile middle{{1,0},{1000000,0},{2000000,1000000},
        {{0.,.001,0,{stripe(1000000,0,200000)}},
         {.001,.002,0,{stripe(1000000,0,200000),stripe(1000000,600000,800000)}},
         {.002,.003,0,{stripe(1000000,0,200000)}}}};
    LocalSectionTile right{{2,0},{2000000,0},{3000000,1000000},
        {{0.,.002,0,{stripe(2000000,600000,800000)}}}};
    StreamTopologyStatistics stats;
    LocalTopologyBuilder invalid(0., {}, [](const SectionPiece&) {}, {}, StreamTopologyLimits{}, &stats);
    invalid.add_tile(left);
    invalid.add_tile(middle);
    invalid.add_tile(right);
    // Middle's first record extends0..3mm; the later1..2mm record is expired
    // at left's2mm bottom. Sharing that Z endpoint cannot anchor left.
    CHECK(stats.processed_islands == 6);
    CHECK(stats.face_records == 4);
    CHECK_THROWS_AS(invalid.finish(), std::invalid_argument);
    LocalTopologyBuilder ordered_invalid(0., {}, [](const SectionPiece&) {}, {}, {}, nullptr, LocalTileOrder::Lexicographic);
    ordered_invalid.add_tile(left);
    ordered_invalid.add_tile(middle);
    ordered_invalid.add_tile(right);
    CHECK_THROWS_AS(ordered_invalid.finish(), std::invalid_argument);
    left.slabs.front().bottom_m = .0015;
    LocalTopologyBuilder valid(0., {}, [](const SectionPiece&) {});
    valid.add_tile(left);
    valid.add_tile(middle);
    valid.add_tile(right);
    CHECK(valid.finish().components() == 2);
    LocalTopologyBuilder ordered_valid(0., {}, [](const SectionPiece&) {}, {}, {}, nullptr, LocalTileOrder::Lexicographic);
    ordered_valid.add_tile(left);
    ordered_valid.add_tile(middle);
    ordered_valid.add_tile(right);
    CHECK(ordered_valid.finish().components() == 2);
    LocalTopologyBuilder reversed(0., {}, [](const SectionPiece&) {});
    reversed.add_tile(right);
    reversed.add_tile(middle);
    reversed.add_tile(left);
    CHECK(reversed.finish().components() == 2);
}

TEST_CASE("Interior vertical runs need one label and no historical face traces", "[ThermoelasticSectionTopology]")
{
    LocalSectionTile tile{{0,0},{0,0},{1000000,1000000},{}};
    for (unsigned i = 0; i < 8; ++i)
        tile.slabs.push_back({double(i)*.001,double(i+1)*.001,0,{box(100000,100000,200000)}});
    StreamTopologyStatistics stats;
    StreamTopologyLimits limits;
    limits.max_islands = 1;
    LocalTopologyBuilder builder(0., {}, [](const SectionPiece&) {}, {}, limits, &stats);
    builder.add_tile(tile);
    const auto certificate = builder.finish();
    CHECK(certificate.components() == 1);
    CHECK(certificate.component_mapping().size() == 1);
    CHECK(stats.processed_islands == 8);
    CHECK(stats.labels == 1);
    CHECK(stats.face_records == 0);
    CHECK(stats.trace_intervals == 0);
}

TEST_CASE("Topology refusal distinguishes escaped coordinates from exact point contacts", "[ThermoelasticSectionTopology]")
{
    LocalSectionTile tile{{3,4},{0,0},{2000000,2000000},{{0.,.001,0,{box(0,0)}}}};
    const auto refusal = [&] {
        LocalTopologyBuilder builder(0., {}, [](const SectionPiece&) {});
        try { builder.add_tile(tile); }
        catch (const std::invalid_argument& error) {
            CHECK_THROWS_AS(builder.finish(), std::invalid_argument);
            return std::string(error.what());
        }
        return std::string{};
    };
    tile.slabs[0].footprint.native() = {box(-1,0)};
    const auto escaped = refusal();
    CHECK(escaped.find("escapes tile at input") != std::string::npos);
    CHECK(escaped.find("tile=3,4") != std::string::npos);
    CHECK(escaped.find("Znm=0..1000000") != std::string::npos);
    tile.maximum = {4000000,4000000};
    ExPolygon unsupported = box(0,0,4000000);
    Polygon touching{{500000,500000},{1500000,1500000},{500000,1500000},{1500000,500000},
        {500000,500000},{500000,250000},{250000,500000}};
    unsupported.holes.push_back(touching);
    tile.slabs[0].footprint.native() = {unsupported};
    const auto contact = refusal();
    CHECK(contact.find("Crossing point-contact cycle edges") != std::string::npos);
    CHECK(contact.find("stage=input") != std::string::npos);
    CHECK(contact.find("tile=3,4") != std::string::npos);
    CHECK(contact.find("Znm=0..1000000") != std::string::npos);
    CHECK(contact.find("ring=hole:0") != std::string::npos);
    CHECK(contact.find("endpoint_on_edge=0") != std::string::npos);
    CHECK(contact.size() < 512);
}

TEST_CASE("A captured repeated hole vertex preserves both clockwise void cycles", "[ThermoelasticSectionTopology]")
{
    const Points captured{{68630935,38339728},{68639366,38353519},{68640221,38356205},{68642182,38359683},
        {68654038,38396886},{68654453,38400870},{68655311,38403563},{68850311,39535563},
        {68852902,39573486},{68852982,39574240},{68844978,39819239},{68842634,39836493},
        {68841481,39846649},{68861766,39766763},{68854978,39519715},{68659234,38386321},
        {68639366,38353519},{68629899,38337890}};
    auto polygon = box(68000000,38000000,2000000);
    polygon.holes.emplace_back(captured);
    LocalSectionTile tile{{6,3},{60000000,30000000},{70000000,40000000},{{0.,.001,0,{polygon}}}};
    const auto topology = qualify_local_section_topology({tile},0.);
    REQUIRE(topology.components() == 1);
    REQUIRE(topology.pieces().size() == 1);
    const auto& normalized = topology.pieces()[0].slab.footprint.native()[0];
    REQUIRE(normalized.holes.size() == 2);
    std::size_t hole_points = 0;
    for (const auto& hole : normalized.holes) {
        CHECK(hole.area() < 0.);
        hole_points += hole.points.size();
    }
    CHECK(hole_points == captured.size());
    const double expected_area_nm2 = 4e12 - (21656852147. + 1208702.) / 2.;
    CHECK_THAT(normalized.area(),Catch::Matchers::WithinAbs(expected_area_nm2,1.));
    SectionMaterial material;
    material.elasticity_pa.setZero();
    material.elasticity_pa.diagonal().head<3>().setConstant(1000.);
    material.elasticity_pa.diagonal().tail<3>().setConstant(500.);
    material.eigenstrain[2] = -.001;
    const Eigen::Vector3d tiny_void{double(68630935LL+68639366LL+68629899LL)/3.*1e-9,
        double(38339728LL+38353519LL+38337890LL)/3.*1e-9,.0005};
    const auto field = solve_deposited_sections(topology,{material},{.06,.03,0.},{.01,.01,.001},{},{},
        {{1,{.0695,.0385,.0005},0},{2,tiny_void,0}});
    CHECK_THAT(field.occupied_volume_m3,Catch::Matchers::WithinAbs(expected_area_nm2*1e-21,1e-19));
    CHECK_THAT(field.elastic_energy_j,Catch::Matchers::WithinAbs(0.,1e-20));
    REQUIRE(field.samples.size() == 1);
    CHECK(field.samples.front().id == 1);
    CHECK(field.unsampled_points == 1);
    CHECK_THAT(field.samples.front().displacement_m.z(),Catch::Matchers::WithinAbs(-.0000005,1e-12));
}

TEST_CASE("Touching void forests cannot hide a disconnected solid island or two outer attachments", "[ThermoelasticSectionTopology]")
{
    auto polygon = box(0,0,4000000);
    for (const Point p : Points{{500000,1500000},{1500000,2500000},{2500000,1500000},{1500000,500000}}) {
        auto hole = box(p.x(),p.y()).contour;
        hole.make_clockwise();
        polygon.holes.push_back(std::move(hole));
    }
    LocalSectionTile tile{{0,0},{0,0},{4000000,4000000},{{0.,.001,0,{polygon}}}};
    CHECK_THROWS_AS(qualify_local_section_topology({tile},0.),std::invalid_argument);
    // Removing one void opens the loop: a tree of point-touching voids leaves
    // the solid connected and every negative polygon remains unfilled.
    polygon.holes.pop_back();
    tile.slabs[0].footprint.native() = {polygon};
    CHECK(qualify_local_section_topology({tile},0.).components() == 1);
    polygon.contour = Polygon{{0,0},{4000000,0},{4000000,2000000},{4000000,4000000},{0,4000000},{0,2000000}};
    polygon.holes = {Polygon{{0,2000000},{1000000,1000000},{2000000,2000000},{1000000,3000000}},
        Polygon{{2000000,2000000},{3000000,1000000},{4000000,2000000},{3000000,3000000}}};
    for (auto& hole : polygon.holes) hole.make_clockwise();
    tile.slabs[0].footprint.native() = {polygon};
    CHECK_THROWS_AS(qualify_local_section_topology({tile},0.),std::invalid_argument);
    polygon.holes = {Polygon{{1000000,1000000},{1000000,500000},{500000,1000000},
        {1000000,1000000},{1500000,1000000},{1000000,1500000}}};
    tile.slabs[0].footprint.native() = {polygon};
    CHECK_THROWS_AS(qualify_local_section_topology({tile},0.),std::invalid_argument);
}

TEST_CASE("Exact outer point contacts split into unwelded islands without changing directed edges", "[ThermoelasticSectionTopology]")
{
    ExPolygon touching(Polygon{{0,0},{1000000,0},{1000000,1000000},{2000000,1000000},
        {2000000,2000000},{1000000,2000000},{1000000,1000000},{0,1000000}});
    Polygon hole = box(200000,200000,200000).contour;
    hole.reverse();
    touching.holes.push_back(hole);
    LocalSectionTile tile{{0,0},{0,0},{2000000,2000000},{{0.,.001,0,{touching}}}};
    const auto topology = qualify_local_section_topology({tile},0.);
    REQUIRE(topology.components() == 2);
    REQUIRE(topology.pieces().size() == 2);
    using Edge = std::array<coord_t,4>;
    const auto edges = [](const ExPolygons& polygons) {
        std::map<Edge,std::size_t> result;
        const auto ring = [&](const Polygon& p) {
            for (std::size_t i = 0; i < p.points.size(); ++i) {
                const auto& a = p.points[i]; const auto& b = p.points[(i+1) % p.points.size()];
                ++result[{a.x(),a.y(),b.x(),b.y()}];
            }
        };
        for (const auto& p : polygons) {
            ring(p.contour);
            for (const auto& h : p.holes) ring(h);
        }
        return result;
    };
    ExPolygons result;
    for (const auto& piece : topology.pieces())
        result.insert(result.end(),piece.slab.footprint.native().begin(),piece.slab.footprint.native().end());
    CHECK(edges(result) == edges({touching}));
    CHECK(topology.pieces()[0].component != topology.pieces()[1].component);
    double area = 0.;
    for (const auto& p : result) area += p.area();
    CHECK_THAT(area, Catch::Matchers::WithinAbs(touching.area(),1.));
    SectionTopologyLimits bounded;
    bounded.max_component_pairs = 4;
    CHECK_THROWS_AS(qualify_local_section_topology({tile},0.,bounded), ResourceLimit);
    touching.holes.clear();
    touching.holes.push_back(Polygon{{1000000,1000000},{1000000,1200000},{1200000,1200000}});
    tile.slabs[0].footprint.native() = {touching};
    CHECK_THROWS_AS(qualify_local_section_topology({tile},0.), std::invalid_argument);
    tile.slabs[0].footprint.native() = {ExPolygon(Polygon{{0,0},{2000000,0},{2000000,2000000},{0,2000000},
        {0,0},{1500000,0},{1500000,1000000},{0,1000000}})};
    CHECK_THROWS_AS(qualify_local_section_topology({tile},0.), std::invalid_argument);
}

TEST_CASE("A captured negative boundary cycle remains an unfilled point-attached void", "[ThermoelasticSectionTopology]")
{
    // Native coordinates of the actual six-edge negative cycle. The surrounding
    // rectangle isolates its topology without depending on the private project.
    const Points cycle{{42815718,120000000},{42814214,127735324},{42816357,127751067},
        {42816320,127747043},{42815210,127733443},{42818206,126638839}};
    Points contour = cycle;
    contour.push_back(cycle.front());
    for (const auto& p : Points{{50000000,120000000},{50000000,130000000},
                               {40000000,130000000},{40000000,120000000}}) contour.push_back(p);
    LocalSectionTile tile{{4,12},{40000000,120000000},{50000000,130000000},
        {{0.,.001,0,{ExPolygon(Polygon(std::move(contour)))}}}};
    const auto topology = qualify_local_section_topology({tile},0.);
    REQUIRE(topology.components() == 1);
    REQUIRE(topology.pieces().size() == 1);
    const auto& polygon = topology.pieces()[0].slab.footprint.native().front();
    REQUIRE(polygon.holes.size() == 1);
    CHECK_THAT(polygon.holes.front().area(), Catch::Matchers::WithinAbs(-30338564527. / 2.,1.));
    const double expected_volume = (1e14 - 30338564527. / 2.) * 1e-21;
    SectionMaterial material;
    material.elasticity_pa.setZero();
    material.elasticity_pa.diagonal().head<3>().setConstant(1000.);
    material.elasticity_pa.diagonal().tail<3>().setConstant(500.);
    material.eigenstrain[2] = -.001;
    const std::vector<SectionSampleRequest> samples{{1,{.045,.125,.0005},0},
        {2,{.0428158375,.126267619333,.0005},0}};
    // One XY cell avoids introducing independently rounded cut vertices into
    // this exact area check; basis refinement/clip error is a separate concern.
    const auto field = solve_deposited_sections(topology,{material},{.04,.12,0.},{.01,.01,.001},{},{},samples);
    CHECK_THAT(field.occupied_volume_m3,Catch::Matchers::WithinAbs(expected_volume,1e-19));
    CHECK_THAT(field.elastic_energy_j,Catch::Matchers::WithinAbs(0.,1e-20));
    REQUIRE(field.samples.size() == 1);
    CHECK(field.samples.front().id == 1);
    CHECK(field.unsampled_points == 1);
    CHECK_THAT(field.samples.front().displacement_m.z(),Catch::Matchers::WithinAbs(-.0000005,1e-12));
}

TEST_CASE("Point-attached void certification rejects multiple contacts crossings and exterior holes", "[ThermoelasticSectionTopology]")
{
    ExPolygon polygon(Polygon{{0,0},{1000000,0},{2000000,0},{2000000,2000000},
                              {1000000,2000000},{0,2000000}});
    LocalSectionTile tile{{0,0},{0,0},{3000000,3000000},{}};
    const auto refused = [&](const Polygon& hole) {
        polygon.holes = {hole};
        tile.slabs = {{0.,.001,0,{polygon}}};
        CHECK_THROWS_AS(qualify_local_section_topology({tile},0.), std::invalid_argument);
    };
    // The void connects two distinct outer endpoints: not a single attached
    // puncture and may disconnect the remaining material.
    refused(Polygon{{1000000,0},{1000000,2000000},{1500000,1000000}});
    // Shares one original endpoint but crosses another outer edge.
    refused(Polygon{{1000000,0},{1500000,1500000},{2500000,500000}});
    // Only touches at a corner, but its interior is outside the solid.
    refused(Polygon{{2000000,2000000},{2000000,2500000},{2500000,2000000}});
}

TEST_CASE("A coarse bounding box cannot anchor material over a real hole", "[ThermoelasticSectionTopology]")
{
    auto ring = box(0,0,3000000);
    auto hole = box(500000,500000,2000000).contour;
    hole.reverse();
    ring.holes.push_back(hole);
    const DepositedSlab lower{0.,.001,0,{ring}};
    const auto valid = qualify_section_topology({lower},0.);
    REQUIRE(valid.pieces().size() == 1);
    CHECK(valid.pieces()[0].slab.footprint.native()[0].holes.size() == 1);
    CHECK_THROWS_AS(qualify_section_topology({lower,{.001,.002,0,{box(1000000,1000000)}}},0.), std::invalid_argument);
}

TEST_CASE("Section topology rejects unresolved overlaps and obeys resource limits", "[ThermoelasticSectionTopology]")
{
    const DepositedSlab a{0.,.001,0,{box(0,0)}};
    const auto duplicate = qualify_section_topology({a,a},0.);
    REQUIRE(duplicate.pieces().size() == 1);
    CHECK(duplicate.components() == 1);
    CHECK_THROWS_AS(qualify_section_topology({a,{0.,.001,1,a.footprint.native()}},0.), LayerMaterialOverlap);
    CHECK_THROWS_AS(qualify_section_topology({a,{.0005,.002,0,a.footprint.native()}},0.), std::invalid_argument);
    CHECK_THROWS_AS(qualify_section_topology({a},.0005), std::invalid_argument);
    CHECK_THROWS_AS(qualify_section_topology({{0.,.001,0,{}}},0.), std::invalid_argument);
    SectionTopologyLimits limits;
    limits.max_slabs = 1;
    CHECK_THROWS_AS(qualify_section_topology({a,a},0.,limits), ResourceLimit);
    limits = {};
    limits.max_boundary_points = 3;
    CHECK_THROWS_AS(qualify_section_topology({a},0.,limits), ResourceLimit);
    limits = {};
    limits.max_component_pairs = 1;
    CHECK_THROWS_AS(qualify_section_topology({a},0.,limits), ResourceLimit);
    CHECK_THROWS_AS(qualify_section_topology({a},0.,{},[] {throw std::runtime_error("cancel");}), std::runtime_error);
}
