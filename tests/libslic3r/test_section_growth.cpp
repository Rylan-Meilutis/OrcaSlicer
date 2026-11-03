#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/SectionGrowth.hpp"
#include "libslic3r/SectionThermalAssembly.hpp"
#include "libslic3r/ThermoelasticSectionTopology.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticExactSections.hpp"
#include "libslic3r/ThermoelasticSectionFootprint.hpp"
#include "libslic3r/ThermoelasticLocalSections.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/libslic3r.h"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;
// Internal arithmetic seam: these raw partitions do not certify topology.
namespace Slic3r::Thermoelastic::GrowthDetail {
std::vector<std::size_t> index_work_quotas(std::size_t budget,const std::vector<std::size_t>& inputs);
void subtract_partition(const std::vector<SectionPiece>& source, const std::vector<SectionPiece>& cover,
    bool old, std::size_t point_limit, const std::function<void(std::size_t)>& charge,
    const std::function<std::size_t()>& remaining_work,
    const std::function<void(std::size_t, std::size_t)>& contact,
    const std::function<void(const SectionPiece&)>& added);
}

TEST_CASE("Shared reference index quotas account for the combined planner without increasing the total", "[SectionGrowth]")
{
    using GrowthDetail::index_work_quotas;
    // Planner contains reference plus all deposited inputs. Fixed reference
    // and the two material indices together contain the same total again.
    const auto quotas=index_work_quotas(100000000,{2000000,1000000,600000,400000});
    REQUIRE(quotas.size()==4);
    CHECK(quotas[0]==50000000);
    CHECK(quotas[1]==25000000);
    CHECK(quotas[2]==15000000);
    CHECK(quotas[3]==10000000);
    const auto rounded=index_work_quotas(101,{14,7,4,3});
    std::size_t total=0;
    for (const auto quota : rounded) { CHECK(quota>0); total+=quota; }
    CHECK(total<=101);
    CHECK_THROWS_AS(index_work_quotas(2,{14,7,4,3}),ResourceLimit);
    CHECK_THROWS_AS(index_work_quotas(100,{1,0}),std::invalid_argument);
}
namespace {
ExPolygon rectangle(coord_t left,coord_t right)
{
    return ExPolygon(Polygon{{left,0},{right,0},{right,1000000000},{left,1000000000}});
}
}

TEST_CASE("Certified growth identifies actual joins and integrates only the connecting material", "[SectionGrowth]")
{
    const auto previous = qualify_section_topology({{0.,1.,0,{rectangle(0,400000000),rectangle(600000000,1000000000)}}},0.);
    const auto next = qualify_section_topology({{0.,1.,0,{rectangle(0,1000000000)}}},0.);
    REQUIRE(previous.components() == 2);
    REQUIRE(next.components() == 1);
    const auto growth = qualify_section_growth(previous,next);
    REQUIRE(growth.previous_components().size() == 2);
    CHECK(growth.previous_components()[0] == 0);
    CHECK(growth.previous_components()[1] == 0);
    const auto added = assemble_section_thermal_addition(growth,
        {{Eigen::Matrix3d::Identity(),2.,400.}},Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones());
    CHECK_THAT(added.occupied_volume_m3,WithinAbs(.2,1e-12));
    double heat = 0.;
    for (const auto& node : added.initial.deposited_cells)
        heat += node.volume_m3*node.volumetric_heat_capacity_j_m3_k*node.temperature_k;
    CHECK_THAT(heat,WithinAbs(160.,1e-9));
    CHECK(added.horizontal_faces.empty());
    CHECK(added.vertical_faces.empty());
}

TEST_CASE("Growth retains only new height and does not require the increment itself to touch the bed", "[SectionGrowth]")
{
    const auto previous = qualify_section_topology({{0.,1.,0,{rectangle(0,1000000000)}}},0.);
    const auto next = qualify_section_topology({{0.,2.,0,{rectangle(0,1000000000)}}},0.);
    const auto growth = qualify_section_growth(previous,next);
    REQUIRE(growth.added_pieces().size() == 1);
    CHECK_THAT(growth.added_pieces()[0].slab.bottom_m,WithinAbs(1.,1e-12));
    const auto added = assemble_section_thermal_addition(growth,
        {{Eigen::Matrix3d::Identity(),2.,400.}},Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones());
    CHECK_THAT(added.occupied_volume_m3,WithinAbs(1.,1e-12));
    const auto repartitioned = qualify_section_topology({{0.,1.,0,{rectangle(0,1000000000)}},
        {1.,2.,0,{rectangle(0,1000000000)}}},0.);
    CHECK(qualify_section_growth(next,repartitioned).added_pieces().empty());
}

TEST_CASE("Growth rejects material replacement and missing old geometry even at equal volume", "[SectionGrowth]")
{
    const auto previous = qualify_section_topology({{0.,1.,0,{rectangle(0,1000000000)}}},0.);
    const auto replaced = qualify_section_topology({{0.,1.,1,{rectangle(0,1000000000)}}},0.);
    CHECK_THROWS_AS(qualify_section_growth(previous,replaced),std::invalid_argument);
    const auto moved = qualify_section_topology({{0.,1.,0,{rectangle(100000000,1100000000)}}},0.);
    CHECK_THROWS_AS(qualify_section_growth(previous,moved),std::invalid_argument);
    const auto split = qualify_section_topology({{0.,1.,0,{rectangle(0,400000000),rectangle(600000000,1000000000)}}},0.);
    CHECK_THROWS_AS(qualify_section_growth(previous,split),std::invalid_argument);
    SectionGrowthLimits limits;
    limits.max_work = 1;
    CHECK_THROWS_AS(qualify_section_growth(previous,previous,limits),ResourceLimit);
    CHECK_THROWS_AS(qualify_section_growth(previous,previous,{},[] { throw std::runtime_error("cancel"); }),std::runtime_error);
}

TEST_CASE("Growth partition subtraction preserves rational overlap and exact added geometry", "[SectionGrowth][Regression]")
{
    const ExPolygons triangle{ExPolygon(Polygon{{0,0},{10,0},{10,4}})};
    const auto full = make_exact_section(triangle).section;
    const auto strip = make_exact_section({ExPolygon(Polygon{{0,0},{1,0},{1,10},{0,10}})}).section;
    // The old section has area .2 nm^2 and a rational vertex at (1,.4).
    // Rounding that section onto native vertices would erase it.
    const auto small = boolean_exact_sections(full,strip,ExactSectionOperation::Intersection).section;
    REQUIRE_FALSE(small.empty());
    const std::vector<SectionPiece> previous{{{0.,.001,0,SectionFootprint(small)},0}};
    const std::vector<SectionPiece> next{{{0.,.001,0,triangle},0}};
    // The first native vertex is the only vertex at X=0. Index bounds must
    // retain it, otherwise the positive exact overlap is never examined.
    const auto native_bounds = next.front().slab.footprint.bounds();
    REQUIRE(native_bounds.has_value());
    CHECK(native_bounds->min == Point(0,0));
    CHECK(native_bounds->max == Point(10,4));
    std::size_t budget = 4000000, work = 0, contacts = 0, checkpoints = 0;
    bool cancel = false;
    const auto charge = [&](std::size_t amount) {
        if (cancel && ++checkpoints == 20) throw std::runtime_error("cancel exact growth");
        if (amount > budget-work) throw ResourceLimit("test growth budget exhausted");
        work += amount;
    };
    const auto remaining = [&] { return budget-work; };
    const auto contact = [&](std::size_t old, std::size_t current) {
        CHECK(old == 0);
        CHECK(current == 0);
        ++contacts;
    };
    std::vector<SectionPiece> additions;
    const auto added_piece = [&](const SectionPiece& piece) { additions.push_back(piece); };
    const auto subtract = [&](const auto& source, const auto& cover, bool old, std::size_t point_limit = 4096) {
        GrowthDetail::subtract_partition(source,cover,old,point_limit,charge,remaining,contact,added_piece);
    };
    subtract(previous,next,true);
    CHECK(contacts == 1);
    CHECK(additions.empty());
    subtract(next,previous,false);
    REQUIRE(additions.size() == 1);
    const auto& added = additions.front();
    REQUIRE(added.slab.footprint.is_exact());
    CHECK(added.component == 0);
    CHECK(added.slab.material == 0);
    const auto expected = boolean_exact_sections(full,small,ExactSectionOperation::Difference).section;
    CHECK(boolean_exact_sections(added.slab.footprint.exact(),expected,
        ExactSectionOperation::SymmetricDifference).section.empty());
    additions.clear();
    subtract(previous,previous,false);
    CHECK(additions.empty());
    CHECK_THROWS_WITH(subtract(next,previous,true),
        "Existing section geometry disappeared during deposition.");
    const std::vector<SectionPiece> replaced{{{0.,.001,1,SectionFootprint(small)},0}};
    CHECK_THROWS_WITH(subtract(previous,replaced,true),
        "Existing section material changed during deposition.");
    budget = work + 20;
    CHECK_THROWS_AS(subtract(previous,next,true),ResourceLimit);
    budget = work + 4000000;
    CHECK_THROWS_AS(subtract(previous,next,true,2),ResourceLimit);
    cancel = true;
    CHECK_THROWS_AS(subtract(previous,next,true),std::runtime_error);
    CHECK(checkpoints == 20);
}

TEST_CASE("Exact reference activation retains rational material cuts and excluded volume", "[SectionGrowth]")
{
    const auto triangle=make_exact_section({ExPolygon(Polygon{{0,0},{1000000001,0},{1000000001,400000000}})});
    const auto reference=qualify_section_topology({{0.,1.,1,SectionFootprint(triangle.section)}},0.);
    const ExPolygon strip(Polygon{{0,0},{500000000,0},{500000000,400000000},{0,400000000}});
    const auto source=qualify_section_topology({{0.,1.,1,{strip}}},0.);
    const auto activated=activate_reference_sections(reference,{source,source});
    REQUIRE(activated.geometry.pieces().size()==1);
    const auto& slab=activated.geometry.pieces().front().slab;
    CHECK(slab.footprint.is_exact());
    CHECK(slab.material==1);
    CHECK(slab.boundary_uncertainty_nm==0);
    const auto exact_strip=make_exact_section({strip});
    const auto expected=boolean_exact_sections(triangle.section,exact_strip.section,ExactSectionOperation::Intersection);
    CHECK(boolean_exact_sections(slab.footprint.exact(),expected.section,ExactSectionOperation::SymmetricDifference).section.empty());
    const double occupied=50000000./1000000001.;
    CHECK_THAT(activated.excluded_nominal_volume_m3,WithinAbs(.2-occupied,1e-14));
    // Exact source and native reference take the same dispatch and retain tags.
    const auto reverse=activate_reference_sections(source,{reference});
    REQUIRE(reverse.geometry.pieces().size()==1);
    CHECK(reverse.geometry.pieces().front().slab.footprint.is_exact());
    CHECK(boolean_exact_sections(reverse.geometry.pieces().front().slab.footprint.exact(),expected.section,
        ExactSectionOperation::SymmetricDifference).section.empty());
    SectionGrowthLimits limits;
    limits.max_work=100;
    CHECK_THROWS_AS(activate_reference_sections(reference,{source},limits),ResourceLimit);
    limits={}; limits.max_points=1;
    CHECK_THROWS_AS(activate_reference_sections(reference,{source},limits),ResourceLimit);
    CHECK_THROWS_AS(activate_reference_sections(reference,{source},{},{},[] { throw std::runtime_error("cancel exact activation"); }),std::runtime_error);
}

TEST_CASE("Shared reference sampling aligns sloped material planes without adding source material", "[SectionGrowth]")
{
    const std::vector<DepositedSegment> first{{{.2,.2,.8},{.8,.2,1.},.1,.2,0}};
    const std::vector<DepositedSegment> second{{{.2,.7,.7},{.8,.7,.7},.1,.15,1}};
    const std::vector<const std::vector<DepositedSegment>*> sources{&first,&second};
    LocalSectionLimits limits;
    auto schedule=std::make_shared<LocalSectionSamplingSchedule>(sources,Eigen::Vector2d::Zero(),
        Eigen::Vector2d::Ones(),0.,.2,limits);
    const auto intervals=schedule->intervals({0,0});
    REQUIRE(intervals.size()>=4);
    CHECK(intervals.front().first==550000000);
    CHECK(intervals.back().second==1000000000);
    std::vector<double> first_planes,second_planes;
    std::size_t first_pieces=0;
    const auto sample = [&](const auto& source,auto& planes) {
        LocalSectionCursor cursor(source,Eigen::Vector2d::Zero(),Eigen::Vector2d::Ones(),0.,.2,limits,
            [&](double z) -> std::optional<SectionOwnership> { planes.push_back(z); return {}; },{},1,nullptr,nullptr,schedule);
        LocalSectionPacketVisitor visitor;
        visitor.begin_tile=[](const LocalSectionTile&) {};
        visitor.band=[&](const LocalSectionTile& packet) {
            for (const auto& slab : packet.slabs) {
                CHECK(slab.material==source.front().material);
                if (slab.material==0) {
                    ++first_pieces;
                    CHECK(slab.bottom_m>=.6-1e-12);
                }
            }
        };
        visitor.end_tile=[] {};
        while (cursor.advance_packet(visitor)) {}
    };
    sample(first,first_planes);
    sample(second,second_planes);
    REQUIRE(first_planes.size()==intervals.size());
    REQUIRE(second_planes.size()==first_planes.size());
    for (std::size_t i=0;i<intervals.size();++i) {
        CHECK_THAT(first_planes[i],WithinAbs(second_planes[i],0.));
        CHECK_THAT(first_planes[i],WithinAbs(double(intervals[i].first+(intervals[i].second-intervals[i].first)/2)*1e-9,0.));
    }
    CHECK(first_pieces>0);
    CHECK_THROWS_AS(schedule->intervals({-1,0}),std::invalid_argument);
    limits.max_sections=1;
    auto exhausted=std::make_shared<LocalSectionSamplingSchedule>(sources,Eigen::Vector2d::Zero(),
        Eigen::Vector2d::Ones(),0.,.2,limits);
    CHECK_THROWS_AS(exhausted->intervals({0,0}),ResourceLimit);
    bool cancel=false;
    LocalSectionSamplingSchedule cancelled(sources,Eigen::Vector2d::Zero(),Eigen::Vector2d::Ones(),0.,.2,{},
        [&] { if (cancel) throw std::runtime_error("cancel shared sampling"); });
    cancel=true;
    CHECK_THROWS_AS(cancelled.intervals({0,0}),std::runtime_error);
}

TEST_CASE("Shared coarse sampling retains a positive one nanometre source interval", "[SectionGrowth]")
{
    const std::vector<DepositedSegment> thin{{{.001,.001,999999e-9},{.002,.001,999999e-9},
        .0001,1e-9,0,BeadSection::RectangularEnvelope}};
    const std::vector<DepositedSegment> carrier{{{.004,.001,.002},{.005,.001,.002},
        .0001,.002,1,BeadSection::RectangularEnvelope}};
    LocalSectionLimits limits;
    limits.partition_increment_m=1e-6;
    limits.geometry.half_native_section_planes=true;
    auto schedule=std::make_shared<LocalSectionSamplingSchedule>(
        std::vector<const std::vector<DepositedSegment>*>{&thin,&carrier},Eigen::Vector2d::Zero(),
        Eigen::Vector2d::Constant(.01),0.,.002,limits);
    const auto intervals=schedule->intervals({0,0});
    const std::pair<coord_t,coord_t> retained{999998,999999};
    CHECK(std::find(intervals.begin(),intervals.end(),retained)!=intervals.end());
    bool sampled_interior=false,retained_material=false;
    LocalSectionCursor cursor(thin,Eigen::Vector2d::Zero(),Eigen::Vector2d::Constant(.01),0.,.002,limits,
        [&](double z) -> std::optional<SectionOwnership> {
            if (z>999998e-9 && z<999999e-9) sampled_interior=true;
            return {};
        },{},1,nullptr,nullptr,schedule);
    LocalSectionPacketVisitor visitor;
    visitor.begin_tile=[](const LocalSectionTile&) {};
    visitor.band=[&](const LocalSectionTile& tile) {
        for (const auto& slab : tile.slabs) {
            retained_material=true;
            CHECK(slab.material==0);
            CHECK_THAT(slab.bottom_m,WithinAbs(999998e-9,1e-15));
            CHECK_THAT(slab.top_m,WithinAbs(999999e-9,1e-15));
        }
    };
    visitor.end_tile=[] {};
    while (cursor.advance_packet(visitor)) {}
    CHECK(sampled_interior);
    CHECK(retained_material);
    CHECK(cursor.metrics().volume_m3>0.);
}

TEST_CASE("Shared schedules support interleaved cursors with different tile subsets", "[SectionGrowth]")
{
    const std::vector<DepositedSegment> long_path{{{.2,.2,.8},{2.8,.2,1.},.1,.2,0}};
    const std::vector<DepositedSegment> middle_path{{{1.2,.7,.7},{1.8,.7,.7},.1,.15,1}};
    auto schedule=std::make_shared<LocalSectionSamplingSchedule>(
        std::vector<const std::vector<DepositedSegment>*>{&long_path,&middle_path},
        Eigen::Vector2d::Zero(),Eigen::Vector2d::Ones(),0.,.2);
    LocalSectionCursor first(long_path,Eigen::Vector2d::Zero(),Eigen::Vector2d::Ones(),0.,.2,
        {},{},{},1,nullptr,nullptr,schedule);
    LocalSectionCursor second(middle_path,Eigen::Vector2d::Zero(),Eigen::Vector2d::Ones(),0.,.2,
        {},{},{},1,nullptr,nullptr,schedule);
    std::array<LocalSectionCursor*,2> cursors{&first,&second};
    std::array<std::vector<std::array<int,2>>,2> visited;
    for (;;) {
        std::optional<std::array<int,2>> next;
        for (const auto* cursor : cursors) {
            const auto candidate=cursor->next_tile_index();
            if (candidate && (!next || *candidate<*next)) next=candidate;
        }
        if (!next) break;
        for (std::size_t i=0;i<cursors.size();++i) {
            if (cursors[i]->next_tile_index()!=next) continue;
            LocalSectionPacketVisitor visitor;
            visitor.begin_tile=[&](const LocalSectionTile& tile) { visited[i].push_back(tile.index); };
            visitor.band=[&](const LocalSectionTile& tile) {
                for (const auto& slab : tile.slabs) CHECK(slab.material==i);
            };
            visitor.end_tile=[] {};
            CHECK(cursors[i]->advance(visitor));
        }
    }
    REQUIRE(visited[0].size()==3);
    REQUIRE(visited[1].size()==1);
    CHECK(visited[0][0][0]==0);
    CHECK(visited[0][2][0]==2);
    CHECK(visited[1][0][0]==1);
    CHECK(first.metrics().volume_m3>0.);
    CHECK(second.metrics().volume_m3>0.);
}

TEST_CASE("Fixed reference activation preserves material ownership without activating future paths", "[SectionGrowth]")
{
    const auto reference = qualify_section_topology({
        {0.,1.,0,{rectangle(0,500000000)}}, {0.,1.,1,{rectangle(500000000,1000000000)}}},0.);
    const auto first_path = qualify_section_topology({{0.,1.,0,{rectangle(0,600000000)}}},0.);
    const auto second_path = qualify_section_topology({{0.,1.,1,{rectangle(400000000,1000000000)}}},0.);
    const auto first = activate_reference_sections(reference,{first_path});
    REQUIRE(first.geometry.pieces().size() == 1);
    CHECK(first.geometry.pieces().front().slab.material == 0);
    CHECK_THAT(first.excluded_nominal_volume_m3,WithinAbs(.1,1e-12));
    const auto complete = activate_reference_sections(reference,{first_path,second_path});
    CHECK_THAT(complete.excluded_nominal_volume_m3,WithinAbs(.2,1e-12));
    // Actual same-material geometric coverage, not merely equal total volume.
    CHECK(qualify_section_growth(reference,complete.geometry).added_pieces().empty());
    CHECK(qualify_section_growth(complete.geometry,reference).added_pieces().empty());
    const auto growth = qualify_section_growth(first.geometry,complete.geometry);
    REQUIRE(growth.added_pieces().size() == 1);
    CHECK(growth.added_pieces().front().slab.material == 1);
    const auto added = assemble_section_thermal_addition(growth,
        {{Eigen::Matrix3d::Identity(),2.,400.},{Eigen::Matrix3d::Identity(),3.,500.}},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones());
    CHECK_THAT(added.occupied_volume_m3,WithinAbs(.5,1e-12));
    double heat = 0.;
    for (const auto& node : added.initial.deposited_cells)
        heat += node.volume_m3*node.volumetric_heat_capacity_j_m3_k*node.temperature_k;
    CHECK_THAT(heat,WithinAbs(.5*3.*500.,1e-9));
    const auto repeated = activate_reference_sections(reference,{first_path,second_path,first_path});
    CHECK(qualify_section_growth(complete.geometry,repeated.geometry).added_pieces().empty());
    CHECK_THAT(repeated.excluded_nominal_volume_m3,WithinAbs(.2,1e-12));
    SectionGrowthLimits limited;
    limited.max_work = 1;
    CHECK_THROWS_AS(activate_reference_sections(reference,{first_path},limited),ResourceLimit);
    limited = {};
    limited.max_points = 1;
    CHECK_THROWS_AS(activate_reference_sections(reference,{first_path},limited),ResourceLimit);
    CHECK_THROWS_AS(activate_reference_sections(reference,{first_path},{},{},[] { throw std::runtime_error("cancel"); }),std::runtime_error);
}

TEST_CASE("Reference activation detects overlaps hidden by integer intersection", "[SectionGrowth][Regression]")
{
    // Cutting the upper reference at x=1 rounds (1,.4) to (1,0).
    // Its new lower edge crosses the unchanged lower material in a sliver
    // whose integer intersection collapses, despite its positive exact area.
    // A shared upward rounding phase can represent the junction at (1,1)
    // within 1 nm without collapsing it onto the lower material's bottom edge.
    const ExPolygons lower{ExPolygon(Polygon{{0,0},{1000000,0},{1000000,400000}})};
    const ExPolygons upper{ExPolygon(Polygon{{0,0},{1000000,400000},{1000000,1000000},{0,1000000}})};
    const ExPolygons cut{ExPolygon(Polygon{{1,0},{1000000,0},{1000000,1000000},{1,1000000}})};
    const auto independent_upper=intersection_ex(upper,cut);
    REQUIRE(intersection_ex(lower,independent_upper).empty());
    REQUIRE_FALSE(exact_section_empty(lower,independent_upper,ExactSectionOperation::Intersection).empty);
    const auto reference=qualify_section_topology({{0.,.001,0,lower},{0.,.001,1,upper}},0.);
    const auto lower_path=qualify_section_topology({{0.,.001,0,lower}},0.);
    const auto upper_path=qualify_section_topology({{0.,.001,1,cut}},0.);
    // Joint clipping introduces a representation allowance. The retained API
    // must refuse that result, since it has no local uncertainty contact proof.
    // In particular it must not silently certify the independently clipped pair.
    CHECK_THROWS_WITH(activate_reference_sections(reference,{lower_path,upper_path}),
        "Boundary uncertainty requires local topology contact qualification.");
    const auto clipped=exact_clip_reference_sections({lower,upper},{lower,cut});
    REQUIRE_FALSE(clipped.sections[0].empty());
    REQUIRE_FALSE(clipped.sections[1].empty());
    CHECK(exact_section_empty(clipped.sections[0],clipped.sections[1],ExactSectionOperation::Intersection).empty);
}
