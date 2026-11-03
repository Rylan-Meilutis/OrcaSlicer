#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "libslic3r/ThermoelasticExactSections.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticLocalSections.hpp"
#include "libslic3r/ThermoelasticLocalTopology.hpp"
#include "libslic3r/ThermoelasticSectionContact.hpp"
#include "libslic3r/SectionIntegration.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/libslic3r.h"
#include <Eigen/Core>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <stdexcept>
#include <algorithm>
#include <cstddef>
#include <array>
#include <cmath>
#include <limits>
#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;
namespace {
ExPolygon box(coord_t left, coord_t right)
{
    return ExPolygon(Polygon({Point(left, 0), Point(right, 0), Point(right, 1000), Point(left, 1000)}));
}
}

TEST_CASE("Optional overlay counters mirror charges without changing exact results or failures", "[ThermoelasticExactSections]")
{
    const auto a=make_exact_section({ExPolygon(Polygon{{0,0},{10,0},{10,10},{0,10}})}).section;
    const auto b=make_exact_section({ExPolygon(Polygon{{5,5},{15,5},{15,15},{5,15}})}).section;
    std::size_t baseline_callbacks=0, measured_callbacks=0;
    const auto baseline=boolean_exact_sections(a,b,ExactSectionOperation::Union,4000000,[&] { ++baseline_callbacks; });
    ExactOverlayDiagnostics counters;
    const auto measured=boolean_exact_sections(a,b,ExactSectionOperation::Union,4000000,
        [&] { ++measured_callbacks; },&counters);
    CHECK(measured_callbacks==baseline_callbacks);
    CHECK(measured.work==baseline.work);
    CHECK(boolean_exact_sections(measured.section,baseline.section,ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK(counters.overlays==1);
    CHECK(counters.edges==8);
    CHECK(counters.boxes_projection==17);
    CHECK(counters.sweep==8);
    CHECK(counters.sorting>0);
    CHECK(counters.candidates>0);
    CHECK(counters.active_scans>=counters.candidates);
    CHECK(counters.event_reservation>0);
    const auto sum=[](const ExactOverlayDiagnostics& value) {
        return value.boxes_projection+value.sorting+value.sweep+value.active_scans+value.event_reservation;
    };
    CHECK(sum(counters)<=measured.work);
    CHECK_FALSE(counters.overflow);
    ExactOverlayDiagnostics partial;
    CHECK_THROWS_AS(boolean_exact_sections(a,b,ExactSectionOperation::Union,baseline.work-1,{},&partial),ResourceLimit);
    CHECK(sum(partial)<=baseline.work-1);
    ExactOverlayDiagnostics cancelled;
    CHECK_THROWS_WITH(boolean_exact_sections(a,b,ExactSectionOperation::Union,4000000,
        [] { throw std::runtime_error("overlay diagnostic cancellation"); },&cancelled),"overlay diagnostic cancellation");
    CHECK(cancelled.overlays==0);
    CHECK(sum(cancelled)==0);
    ExactOverlayDiagnostics saturated;
    saturated.active_scans=std::numeric_limits<std::size_t>::max();
    CHECK(boolean_exact_sections(a,b,ExactSectionOperation::Union,4000000,{},&saturated).work==baseline.work);
    CHECK(saturated.overflow);
    CHECK(saturated.active_scans==std::numeric_limits<std::size_t>::max());
    ExactOverlayDiagnostics fast;
    boolean_exact_sections(a,a,ExactSectionOperation::Union,4000000,{},&fast);
    CHECK(fast.overlays==0);
    CHECK(sum(fast)==0);
}

TEST_CASE("Exact work failures identify operations and effective budgets", "[ThermoelasticExactSections]")
{
    using Catch::Matchers::ContainsSubstring;
    const auto section=make_exact_section({ExPolygon(Polygon{{0,0},{5,0},{5,3},{0,3}})}).section;
    const auto boundary=extract_exact_section_boundary(section).boundary;
    const auto separate=make_exact_section({ExPolygon(Polygon{{0,0},{5,0},{5,3},{0,3}})}).section;
    CHECK_THROWS_WITH(boolean_exact_sections(section,separate,ExactSectionOperation::SymmetricDifference,1),
        ContainsSubstring("operation=boolean_exact_sections") && ContainsSubstring("used=") &&
        ContainsSubstring("requested=") && ContainsSubstring("limit=1,") && ContainsSubstring("per_operation_cap=4000000"));
    CHECK_THROWS_WITH(exact_section_vertical_surface_moments(boundary,{0,0,0},{5,3,10},0,10,
        ExactBoundaryCellOwnership::MaterialSide,1),
        ContainsSubstring("operation=exact_section_vertical_surface_moments") && ContainsSubstring("used=") &&
        ContainsSubstring("requested=") && ContainsSubstring("limit=1,") && ContainsSubstring("per_operation_cap=4000000"));
    CHECK_THROWS_WITH(exact_section_horizontal_surface_moments(section,{0,0,0},{5,3,10},0,1),
        ContainsSubstring("operation=exact_section_horizontal_surface_moments") && ContainsSubstring("used=") &&
        ContainsSubstring("requested=") && ContainsSubstring("limit=1,") && ContainsSubstring("per_operation_cap=4000000"));
    CHECK_THROWS_WITH(exact_section_bilinear_extrema(section,{0,0},{5,3},{0.,1.,1.,2.},1),
        ContainsSubstring("operation=exact_section_bilinear_extrema") && ContainsSubstring("used=") &&
        ContainsSubstring("requested=") && ContainsSubstring("limit=1,") && ContainsSubstring("per_operation_cap=4000000"));
    CHECK_THROWS_WITH(validate_exact_section_topology(section,1),
        ContainsSubstring("operation=validate_exact_section_topology") && ContainsSubstring("used=") &&
        ContainsSubstring("requested=") && ContainsSubstring("limit=1,") && ContainsSubstring("per_operation_cap=4000000"));
}

TEST_CASE("First deposition partitions preserve every prefix union and prior material ownership", "[ThermoelasticExactSections]")
{
    using Catch::Matchers::WithinAbs;
    const std::vector<ExactDepositedCommandSection> commands{
        {0,make_exact_section({box(0,1000)}).section},
        {1,make_exact_section({box(500,1500)}).section},
        {2,make_exact_section({box(1000,2000)}).section},
        {0,make_exact_section({box(250,1750)}).section}};
    ExactSection expected;
    ExactDepositedPartition previous;
    for (std::size_t count=1;count<=commands.size();++count) {
        const std::vector<ExactDepositedCommandSection> prefix(commands.begin(),commands.begin()+count);
        const auto result=exact_first_deposition_partition(prefix);
        CHECK(result.policy==DepositedOverlapPolicy::FirstDepositionWins);
        expected=boolean_exact_sections(expected,commands[count-1].envelope,ExactSectionOperation::Union).section;
        CHECK(boolean_exact_sections(result.occupied_union,expected,ExactSectionOperation::SymmetricDifference).section.empty());
        ExactSection assembled;
        for (const auto& region:result.regions) {
            CHECK(boolean_exact_sections(assembled,region.envelope,ExactSectionOperation::Intersection).section.empty());
            assembled=boolean_exact_sections(assembled,region.envelope,ExactSectionOperation::Union).section;
            const auto old=std::find_if(previous.regions.begin(),previous.regions.end(),[&](const auto& entry) {
                return entry.material==region.material;
            });
            if (old!=previous.regions.end())
                CHECK(boolean_exact_sections(old->envelope,region.envelope,ExactSectionOperation::Difference).section.empty());
        }
        CHECK(boolean_exact_sections(assembled,expected,ExactSectionOperation::SymmetricDifference).section.empty());
        previous=result;
    }
    REQUIRE(previous.regions.size()==3);
    CHECK_THAT(previous.cross_material_overlap_area_m2,WithinAbs(1.25e-12,1e-25));
    CHECK(boolean_exact_sections(previous.cross_material_overlap,make_exact_section({box(500,1750)}).section,
        ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK(boolean_exact_sections(previous.regions.front().envelope,commands.front().envelope,
        ExactSectionOperation::SymmetricDifference).section.empty());
    auto reversed=commands;
    std::reverse(reversed.begin(),reversed.end());
    const auto other=exact_first_deposition_partition(reversed);
    CHECK(boolean_exact_sections(other.occupied_union,previous.occupied_union,ExactSectionOperation::SymmetricDifference).section.empty());
    REQUIRE(other.regions.front().material==0);
    CHECK_FALSE(boolean_exact_sections(other.regions.front().envelope,previous.regions.front().envelope,
        ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK_THAT(other.cross_material_overlap_area_m2,WithinAbs(previous.cross_material_overlap_area_m2,1e-25));
}

TEST_CASE("First deposition merges same material commands without artificial barriers", "[ThermoelasticExactSections]")
{
    const auto first=make_exact_section({box(0,1000)}).section;
    const auto second=make_exact_section({box(500,1500)}).section;
    const auto result=exact_first_deposition_partition({{7,first},{7,second},{7,first}});
    REQUIRE(result.regions.size()==1);
    CHECK(result.regions.front().material==7);
    CHECK(result.cross_material_overlap.empty());
    CHECK_THAT(result.cross_material_overlap_area_m2,Catch::Matchers::WithinAbs(0.,1e-30));
    CHECK(boolean_exact_sections(result.regions.front().envelope,make_exact_section({box(0,1500)}).section,
        ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK(exact_section_components(result.regions.front().envelope).sections.size()==1);
}

TEST_CASE("Single owner deposition uses one exact union per command with bounded work", "[ThermoelasticExactSections]")
{
    std::vector<ExactDepositedCommandSection> commands;
    ExactSection expected;
    std::size_t union_work=0;
    for (coord_t i=0;i<8;++i) {
        commands.push_back({7,make_exact_section({box(i*250,i*250+1000)}).section});
        const auto united=boolean_exact_sections(expected,commands.back().envelope,ExactSectionOperation::Union);
        expected=united.section;
        union_work+=united.work;
        const auto prefix=exact_first_deposition_partition(commands);
        REQUIRE(prefix.regions.size()==1);
        CHECK(boolean_exact_sections(prefix.occupied_union,expected,ExactSectionOperation::SymmetricDifference).section.empty());
        CHECK(boolean_exact_sections(prefix.regions.front().envelope,expected,ExactSectionOperation::SymmetricDifference).section.empty());
        CHECK(prefix.cross_material_overlap.empty());
        // Allow linear metadata/validation overhead, not repeated geometry
        // subtraction, union rebuilding and equality overlays for one owner.
        CHECK(prefix.work<=union_work+8*commands.size()+8);
        CHECK(exact_first_deposition_partition(commands,prefix.work).work==prefix.work);
        CHECK_THROWS_AS(exact_first_deposition_partition(commands,prefix.work-1),ResourceLimit);
    }
    std::size_t checkpoints=0;
    CHECK_THROWS_WITH(exact_first_deposition_partition(commands,4000000,[&] {
        if (++checkpoints==commands.size()+4) throw std::runtime_error("cancel single owner union");
    }),"cancel single owner union");
}

TEST_CASE("Single owner growth retains ambiguity from fully occluded foreign deposition", "[ThermoelasticExactSections]")
{
    const auto original=make_exact_section({box(0,1000)}).section;
    const auto hidden=make_exact_section({box(250,750)}).section;
    const auto extension=make_exact_section({box(500,1500)}).section;
    const auto before=exact_first_deposition_partition({{0,original},{1,hidden}});
    REQUIRE(before.regions.size()==1);
    const auto after=exact_first_deposition_partition({{0,original},{1,hidden},{0,extension}});
    REQUIRE(after.regions.size()==1);
    CHECK(after.regions.front().material==0);
    const auto expected=make_exact_section({box(0,1500)}).section;
    CHECK(boolean_exact_sections(after.occupied_union,expected,ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK(boolean_exact_sections(after.regions.front().envelope,expected,ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK(boolean_exact_sections(after.cross_material_overlap,hidden,ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK_THAT(after.cross_material_overlap_area_m2,Catch::Matchers::WithinAbs(before.cross_material_overlap_area_m2,1e-25));
}

TEST_CASE("First deposition retains holes and unique material outside earlier envelopes", "[ThermoelasticExactSections]")
{
    auto ring=box(0,1000);
    ring.holes.push_back(Polygon{{250,250},{250,750},{750,750},{750,250}});
    const auto first=make_exact_section({ring}).section;
    const auto second=make_exact_section({box(0,1500)}).section;
    const auto result=exact_first_deposition_partition({{0,first},{1,second}});
    REQUIRE(result.regions.size()==2);
    CHECK(boolean_exact_sections(result.regions[0].envelope,first,ExactSectionOperation::SymmetricDifference).section.empty());
    const auto exclusive=boolean_exact_sections(second,first,ExactSectionOperation::Difference).section;
    CHECK(boolean_exact_sections(result.regions[1].envelope,exclusive,ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK(exact_closed_contains_point(result.regions[1].envelope,{500,500}).contained);
    CHECK(exact_closed_contains_point(result.regions[1].envelope,{1250,500}).contained);
}

TEST_CASE("First deposition partition work and storage remain bounded and cancellation propagates", "[ThermoelasticExactSections]")
{
    const auto envelope=make_exact_section({box(0,1000)}).section;
    const std::vector<ExactDepositedCommandSection> commands{{0,envelope},{1,envelope}};
    const auto result=exact_first_deposition_partition(commands);
    REQUIRE(result.work>1);
    CHECK(exact_first_deposition_partition(commands,result.work).work==result.work);
    CHECK_THROWS_AS(exact_first_deposition_partition(commands,result.work-1),ResourceLimit);
    CHECK_THROWS_AS(exact_first_deposition_partition(commands,0),std::invalid_argument);
    CHECK_THROWS_AS(exact_first_deposition_partition(commands,4000001),std::invalid_argument);
    CHECK_THROWS_AS(exact_first_deposition_partition({{256,envelope}}),std::invalid_argument);
    CHECK_THROWS_AS(exact_first_deposition_partition(std::vector<ExactDepositedCommandSection>(4097)),ResourceLimit);
    CHECK_THROWS_AS(exact_first_deposition_partition(std::vector<ExactDepositedCommandSection>(1025,{0,envelope})),ResourceLimit);
    CHECK_THROWS_WITH(exact_first_deposition_partition(commands,4000000,
        [] { throw std::runtime_error("cancel first deposition"); }),"cancel first deposition");
    const auto empty=exact_first_deposition_partition({});
    CHECK(empty.occupied_union.empty());
    CHECK(empty.regions.empty());
}

TEST_CASE("Streamed first deposition bounds live cuts without retaining cumulative input geometry", "[ThermoelasticExactSections]")
{
    const auto envelope=make_exact_section({box(0,1000)}).section;
    std::size_t supplied=0;
    const auto provider=[&](std::size_t index,std::size_t remaining) {
        CHECK(index==supplied++);
        CHECK(remaining>0);
        return ExactDepositedCommandResult{{0,envelope},0};
    };
    const auto result=exact_first_deposition_partition_streamed(1025,provider);
    CHECK(supplied==1025);
    REQUIRE(result.regions.size()==1);
    CHECK(result.cross_material_overlap.empty());
    CHECK(boolean_exact_sections(result.occupied_union,envelope,ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK_THROWS_AS(exact_first_deposition_partition(std::vector<ExactDepositedCommandSection>(1025,{0,envelope})),ResourceLimit);
}

TEST_CASE("Streamed and retained first deposition preserve identical ordered material partitions", "[ThermoelasticExactSections]")
{
    const std::vector<ExactDepositedCommandSection> commands{
        {0,make_exact_section({box(0,1000)}).section},
        {1,make_exact_section({box(500,1500)}).section},
        {2,make_exact_section({box(1000,2000)}).section},
        {0,make_exact_section({box(250,1750)}).section}};
    for (std::size_t count=1;count<=commands.size();++count) {
        const std::vector<ExactDepositedCommandSection> prefix(commands.begin(),commands.begin()+count);
        const auto retained=exact_first_deposition_partition(prefix);
        const auto streamed=exact_first_deposition_partition_streamed(count,[&](std::size_t index,std::size_t) {
            return ExactDepositedCommandResult{commands[index],0};
        });
        REQUIRE(streamed.regions.size()==retained.regions.size());
        CHECK(streamed.policy==retained.policy);
        CHECK(boolean_exact_sections(streamed.occupied_union,retained.occupied_union,ExactSectionOperation::SymmetricDifference).section.empty());
        CHECK(boolean_exact_sections(streamed.cross_material_overlap,retained.cross_material_overlap,ExactSectionOperation::SymmetricDifference).section.empty());
        CHECK_THAT(streamed.cross_material_overlap_area_m2,Catch::Matchers::WithinAbs(retained.cross_material_overlap_area_m2,1e-25));
        for (std::size_t i=0;i<retained.regions.size();++i) {
            CHECK(streamed.regions[i].material==retained.regions[i].material);
            CHECK(boolean_exact_sections(streamed.regions[i].envelope,retained.regions[i].envelope,
                ExactSectionOperation::SymmetricDifference).section.empty());
        }
    }
}

TEST_CASE("Contiguous material runs preserve every requested prefix and occluded ambiguity", "[ThermoelasticExactSections]")
{
    const auto cut=[](coord_t left,coord_t right) { return make_exact_section({box(left,right)}).section; };
    const std::vector<ExactDepositedCommandSection> commands{
        {0,cut(0,1200)},{2,{}},{0,cut(200,1400)},
        {1,cut(250,500)},{1,cut(400,1600)},{1,cut(1500,1800)},
        {0,cut(1000,2000)},{2,cut(600,2200)},{2,cut(0,100)},
        {0,{}},{0,cut(1900,2400)}};
    for (std::size_t count=1;count<=commands.size();++count) {
        const std::vector<ExactDepositedCommandSection> prefix(commands.begin(),commands.begin()+count);
        const auto reference=exact_first_deposition_partition(prefix);
        std::size_t calls=0;
        const auto streamed=exact_first_deposition_partition_streamed(count,[&](std::size_t index,std::size_t) {
            CHECK(index==calls++);
            return ExactDepositedCommandResult{commands[index],0};
        });
        CHECK(calls==count);
        REQUIRE(streamed.regions.size()==reference.regions.size());
        CHECK(boolean_exact_sections(streamed.occupied_union,reference.occupied_union,ExactSectionOperation::SymmetricDifference).section.empty());
        CHECK(boolean_exact_sections(streamed.cross_material_overlap,reference.cross_material_overlap,ExactSectionOperation::SymmetricDifference).section.empty());
        CHECK_THAT(streamed.cross_material_overlap_area_m2,Catch::Matchers::WithinAbs(reference.cross_material_overlap_area_m2,1e-25));
        for (std::size_t i=0;i<reference.regions.size();++i) {
            CHECK(streamed.regions[i].material==reference.regions[i].material);
            CHECK(boolean_exact_sections(streamed.regions[i].envelope,reference.regions[i].envelope,
                ExactSectionOperation::SymmetricDifference).section.empty());
        }
        CHECK(exact_first_deposition_partition_streamed(count,[&](std::size_t index,std::size_t) {
            return ExactDepositedCommandResult{commands[index],0};
        },streamed.work).work==streamed.work);
    }
}

TEST_CASE("Contiguous material runs reduce repeated foreign ownership proof work", "[ThermoelasticExactSections]")
{
    std::vector<ExactDepositedCommandSection> commands{{0,make_exact_section({box(0,1000)}).section}};
    for (coord_t i=0;i<32;++i)
        commands.push_back({1,make_exact_section({box(500+i*10,1500+i*10)}).section});
    const auto reference=exact_first_deposition_partition(commands);
    const auto streamed=exact_first_deposition_partition_streamed(commands.size(),[&](std::size_t index,std::size_t) {
        return ExactDepositedCommandResult{commands[index],0};
    });
    CHECK(streamed.work<reference.work);
    REQUIRE(streamed.regions.size()==reference.regions.size());
    for (std::size_t i=0;i<reference.regions.size();++i)
        CHECK(boolean_exact_sections(streamed.regions[i].envelope,reference.regions[i].envelope,
            ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK(boolean_exact_sections(streamed.cross_material_overlap,reference.cross_material_overlap,
        ExactSectionOperation::SymmetricDifference).section.empty());
}

TEST_CASE("Balanced streamed runs retain exact chronological ownership at binary boundaries", "[ThermoelasticExactSections]")
{
    std::vector<ExactDepositedCommandSection> commands;
    std::size_t sequential_union_work=0;
    ExactSection growing;
    for (coord_t i=0;i<64;++i) {
        const auto cut=make_exact_section({ExPolygon(Polygon{{3*i,0},{3*i+1,0},{3*i+1,1},{3*i,1}})}).section;
        commands.push_back({0,cut});
        const auto next=boolean_exact_sections(growing,cut,ExactSectionOperation::Union);
        sequential_union_work+=next.work;
        growing=next.section;
    }
    const auto first_run=exact_first_deposition_partition_streamed(commands.size(),[&](std::size_t i,std::size_t) {
        return ExactDepositedCommandResult{commands[i],0};
    });
    CHECK(first_run.work<sequential_union_work);
    commands.push_back({1,make_exact_section({ExPolygon(Polygon{{0,0},{192,0},{192,2},{0,2}})}).section});
    // Keep the deliberately unbatched foreign-return oracle within its4M
    // work cap; the first64 commands already exercise the full carry tree.
    for (coord_t i=0;i<8;++i)
        commands.push_back({0,make_exact_section({ExPolygon(Polygon{{3*i+1,0},{3*i+2,0},{3*i+2,1},{3*i+1,1}})}).section});
    for (const std::size_t count:{31u,32u,33u,64u,65u,73u}) {
        CAPTURE(count);
        const std::vector<ExactDepositedCommandSection> prefix(commands.begin(),commands.begin()+count);
        ExactDepositedPartition reference;
        REQUIRE_NOTHROW(reference=exact_first_deposition_partition(prefix));
        const auto balanced=exact_first_deposition_partition_streamed(count,[&](std::size_t i,std::size_t) {
            return ExactDepositedCommandResult{commands[i],0};
        });
        REQUIRE(balanced.regions.size()==reference.regions.size());
        CHECK(boolean_exact_sections(balanced.occupied_union,reference.occupied_union,
            ExactSectionOperation::SymmetricDifference).section.empty());
        CHECK(boolean_exact_sections(balanced.cross_material_overlap,reference.cross_material_overlap,
            ExactSectionOperation::SymmetricDifference).section.empty());
        for (std::size_t i=0;i<reference.regions.size();++i) {
            CHECK(balanced.regions[i].material==reference.regions[i].material);
            CHECK(boolean_exact_sections(balanced.regions[i].envelope,reference.regions[i].envelope,
                ExactSectionOperation::SymmetricDifference).section.empty());
        }
        CHECK_THROWS_AS(exact_first_deposition_partition_streamed(count,[&](std::size_t i,std::size_t) {
            return ExactDepositedCommandResult{commands[i],0};
        },balanced.work-1),ResourceLimit);
    }
    std::size_t callbacks=0;
    CHECK_THROWS_WITH(exact_first_deposition_partition_streamed(commands.size(),[&](std::size_t i,std::size_t) {
        return ExactDepositedCommandResult{commands[i],0};
    },4000000,[&] { if (++callbacks==100) throw std::runtime_error("cancel balanced run"); }),"cancel balanced run");
}

TEST_CASE("Streamed run bins consolidate before requesting an unknown incoming cut", "[ThermoelasticExactSections]")
{
    const auto boxes=[](coord_t begin,std::size_t count) {
        ExPolygons polygons;
        for (std::size_t i=0;i<count;++i) {
            const coord_t x=begin+coord_t(i)*3;
            polygons.emplace_back(Polygon{{x,0},{x+1,0},{x+1,1},{x,1}});
        }
        return make_exact_section(polygons).section;
    };
    const auto foreign=boxes(-450,150), repeated=boxes(0,200);
    // Stored foreign owner+occupied=1200 vertices. Seven repeated cuts leave
    // three binary bins (2400 vertices); a further800-vertex cut cannot be
    // requested at that stable boundary without first consolidating the bins.
    const auto reference=exact_first_deposition_partition({{1,foreign},{0,repeated}});
    std::size_t calls=0;
    const auto balanced=exact_first_deposition_partition_streamed(9,[&](std::size_t index,std::size_t) {
        ++calls;
        return ExactDepositedCommandResult{{index==0 ? 1u : 0u,index==0 ? foreign : repeated},0};
    });
    CHECK(calls==9);
    REQUIRE(balanced.regions.size()==reference.regions.size());
    for (std::size_t i=0;i<reference.regions.size();++i)
        CHECK(boolean_exact_sections(balanced.regions[i].envelope,reference.regions[i].envelope,
            ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK(boolean_exact_sections(balanced.occupied_union,reference.occupied_union,
        ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK(balanced.cross_material_overlap.empty());
}

TEST_CASE("Streamed first deposition charges provider work and propagates cancellation", "[ThermoelasticExactSections]")
{
    const auto envelope=make_exact_section({box(0,1000)}).section;
    const auto provider=[&](std::size_t,std::size_t remaining) {
        CHECK(remaining>=17);
        return ExactDepositedCommandResult{{0,envelope},17};
    };
    const auto result=exact_first_deposition_partition_streamed(2,provider);
    const auto free_input=exact_first_deposition_partition_streamed(2,[&](std::size_t,std::size_t) {
        return ExactDepositedCommandResult{{0,envelope},0};
    });
    CHECK(result.work==free_input.work+34);
    CHECK(exact_first_deposition_partition_streamed(2,provider,result.work).work==result.work);
    CHECK_THROWS_AS(exact_first_deposition_partition_streamed(2,provider,result.work-1),ResourceLimit);
    CHECK_THROWS_AS(exact_first_deposition_partition_streamed(1,[&](std::size_t,std::size_t remaining) {
        return ExactDepositedCommandResult{{0,envelope},remaining+1};
    }),ResourceLimit);
    CHECK_THROWS_WITH(exact_first_deposition_partition_streamed(1,[](std::size_t,std::size_t) -> ExactDepositedCommandResult {
        throw std::runtime_error("provider cancellation");
    }),"provider cancellation");
    CHECK_THROWS_WITH(exact_first_deposition_partition_streamed(1,provider,4000000,
        [] { throw std::runtime_error("stream cancellation"); }),"stream cancellation");
    std::size_t calls=0;
    CHECK_THROWS_AS(exact_first_deposition_partition_streamed(4097,[&](std::size_t,std::size_t) {
        ++calls;
        return ExactDepositedCommandResult{};
    }),ResourceLimit);
    CHECK(calls==0);
    CHECK_THROWS_AS(exact_first_deposition_partition_streamed(0,{}),std::invalid_argument);
    CHECK_THROWS_AS(exact_first_deposition_partition_streamed(1,provider,0),std::invalid_argument);
    CHECK_THROWS_AS(exact_first_deposition_partition_streamed(1,provider,4000001),std::invalid_argument);
    CHECK_THROWS_AS(exact_first_deposition_partition_streamed(1,[](std::size_t,std::size_t) {
        return ExactDepositedCommandResult{{256,{}},0};
    }),std::invalid_argument);
}

TEST_CASE("Streamed first deposition rejects oversized incoming and live geometry", "[ThermoelasticExactSections]")
{
    const auto boxes=[](coord_t begin,std::size_t count) {
        ExPolygons polygons;
        for (std::size_t i=0;i<count;++i) {
            const coord_t x=begin+coord_t(i)*3;
            polygons.emplace_back(Polygon{{x,0},{x+1,0},{x+1,1},{x,1}});
        }
        return make_exact_section(polygons).section;
    };
    const auto first=boxes(0,200),second=boxes(600,200),subset=boxes(0,256);
    const auto united=boolean_exact_sections(first,second,ExactSectionOperation::Union).section;
    REQUIRE(first.vertex_count()==800);
    REQUIRE(united.vertex_count()==1600);
    // The old retained API still admits an individual input larger than1024.
    REQUIRE(exact_first_deposition_partition({{0,united}}).regions.size()==1);
    CHECK_THROWS_WITH(exact_first_deposition_partition_streamed(1,[&](std::size_t,std::size_t) {
        return ExactDepositedCommandResult{{0,united},0};
    }),Catch::Matchers::ContainsSubstring("cut vertex budget"));
    const std::vector<ExactDepositedCommandSection> commands{{0,first},{0,second},{1,boxes(0,1)},{0,subset}};
    REQUIRE(exact_first_deposition_partition(commands).regions.size()==1);
    // After the first run,3200 stored vertices (owner + occupied), a pending
    // four-vertex foreign run and its1024-vertex lookahead exceed4096.
    CHECK_THROWS_WITH(exact_first_deposition_partition_streamed(commands.size(),[&](std::size_t index,std::size_t) {
        return ExactDepositedCommandResult{commands[index],0};
    }),Catch::Matchers::ContainsSubstring("live vertex budget"));
}

TEST_CASE("Exact overlay projection avoids dense scans without changing separated strip unions", "[ThermoelasticExactSections]")
{
    const bool rotated=GENERATE(false,true);
    constexpr coord_t strips=128;
    ExPolygons first,second,all;
    for (coord_t i=0;i<2*strips;++i) {
        const coord_t y=-1000+3*i;
        Polygon ring{{-100000,y},{100000,y},{100000,y+1},{-100000,y+1}};
        if (rotated) for (auto& point:ring.points) {
            const auto x=point.x(); point.x()=-point.y(); point.y()=x;
        }
        ExPolygon polygon(std::move(ring));
        (i%2 ? second : first).push_back(polygon);
        all.push_back(std::move(polygon));
    }
    const auto a=make_exact_section(first).section,b=make_exact_section(second).section;
    const auto actual=boolean_exact_sections(a,b,ExactSectionOperation::Union);
    const auto reference=make_exact_section(all).section;
    CHECK(boolean_exact_sections(actual.section,reference,ExactSectionOperation::SymmetricDifference).section.empty());
    // With an X-only sweep in the unrotated case, each of the2*N long
    // horizontal edges from A scans every one from B: at least4*N*N visits,
    // before sorting/event/output charges. The selected sparse projection
    // must beat that bound including ALL selection and Boolean work.
    CHECK(actual.work<std::size_t(4*strips*strips));
    const std::array<coord_t,3> start{-100000,-100000,0},size{200000,200000,10};
    const auto moments=exact_section_volume_moments(actual.section,start,size,0,10);
    const auto expected=exact_section_volume_moments(reference,start,size,0,10);
    for (std::size_t i=0;i<27;++i)
        CHECK_THAT(moments.values[i],Catch::Matchers::WithinRel(expected.values[i],1e-13));
    CHECK_THROWS_AS(boolean_exact_sections(a,b,ExactSectionOperation::Union,actual.work-1),ResourceLimit);
    CHECK_THROWS_WITH(boolean_exact_sections(a,b,ExactSectionOperation::Union,4000000,
        [] { throw std::runtime_error("projection cancelled"); }),"projection cancelled");
}

TEST_CASE("Exact identity metadata gates never substitute bounds for occupied geometry", "[ThermoelasticExactSections]")
{
    const ExPolygon square(Polygon{{0,0},{10,0},{10,10},{0,10}});
    const auto a=make_exact_section({square}).section;
    const auto same=make_exact_section({square}).section;
    const auto diamond=make_exact_section({ExPolygon(Polygon{{0,5},{5,0},{10,5},{5,10}})}).section;
    const auto shifted=make_exact_section({ExPolygon(Polygon{{5,0},{15,0},{15,10},{5,10}})}).section;
    const auto identity=boolean_exact_sections(a,same,ExactSectionOperation::Union);
    const auto different=boolean_exact_sections(a,diamond,ExactSectionOperation::SymmetricDifference);
    CHECK_FALSE(different.section.empty());
    CHECK(identity.work<different.work);
    const auto containing=boolean_exact_sections(a,diamond,ExactSectionOperation::Union).section;
    CHECK(boolean_exact_sections(containing,a,ExactSectionOperation::SymmetricDifference).section.empty());
    const auto expected=make_exact_section({square,ExPolygon(Polygon{{5,0},{15,0},{15,10},{5,10}})}).section;
    const auto united=boolean_exact_sections(a,shifted,ExactSectionOperation::Union).section;
    CHECK(boolean_exact_sections(united,expected,ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK_THROWS_AS(boolean_exact_sections(a,same,ExactSectionOperation::Union,identity.work-1),ResourceLimit);
}

TEST_CASE("Adaptive exact overlays preserve closed touching and holed boundaries", "[ThermoelasticExactSections]")
{
    const bool rotated=GENERATE(false,true);
    ExPolygon outer(Polygon{{-100,-10},{100,-10},{100,10},{-100,10}});
    outer.holes.emplace_back(Polygon{{-80,-5},{-80,5},{80,5},{80,-5}});
    ExPolygon touching(Polygon{{-100,10},{0,10},{100,10},{100,20},{-100,20}});
    if (rotated) {
        // A quarter turn preserves orientation, including holes.
        for (auto* polygon:{&outer,&touching}) {
            const auto rotate=[](Polygon& ring) {
                for (auto& p:ring.points) { const auto x=p.x(); p.x()=-p.y(); p.y()=x; }
            };
            rotate(polygon->contour);
            for (auto& hole:polygon->holes) rotate(hole);
        }
    }
    const auto a=make_exact_section({outer}).section,b=make_exact_section({touching}).section;
    const auto expected=make_exact_section({outer,touching}).section;
    const auto united=boolean_exact_sections(a,b,ExactSectionOperation::Union).section;
    CHECK(boolean_exact_sections(united,expected,ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK(boolean_exact_sections(a,b,ExactSectionOperation::Intersection).section.empty());
    const auto duplicate=make_exact_section({touching}).section;
    CHECK(boolean_exact_sections(b,duplicate,ExactSectionOperation::SymmetricDifference).section.empty());
}

TEST_CASE("Strictly separated exact bounds avoid overlays without changing geometry", "[ThermoelasticExactSections]")
{
    const coord_t offset=GENERATE(coord_t(0),coord_t(-50));
    const bool vertical=GENERATE(false,true);
    ExPolygon holed(Polygon{{offset,offset},{offset+10,offset},{offset+10,offset+10},{offset,offset+10}});
    holed.holes.emplace_back(Polygon{{offset+3,offset+3},{offset+3,offset+7},{offset+7,offset+7},{offset+7,offset+3}});
    const auto left=make_exact_section({holed}).section;
    const coord_t x=offset+(vertical ? 0 : 30), y=offset+(vertical ? 30 : 0);
    const auto right=make_exact_section({ExPolygon(Polygon{{x,y},{x+10,y},{x,y+3}})}).section;
    const auto intersection=boolean_exact_sections(left,right,ExactSectionOperation::Intersection);
    const auto difference=boolean_exact_sections(left,right,ExactSectionOperation::Difference);
    CHECK(intersection.section.empty());
    CHECK(boolean_exact_sections(difference.section,left,ExactSectionOperation::SymmetricDifference).section.empty());
    const std::array<coord_t,3> start{offset,offset,0}, size{10,10,10};
    const auto expected=exact_section_volume_moments(left,start,size,0,10);
    const auto actual=exact_section_volume_moments(difference.section,start,size,0,10);
    for (std::size_t i=0;i<27;++i)
        CHECK_THAT(actual.values[i],Catch::Matchers::WithinRel(expected.values[i],1e-13));
    // Union and symmetric difference remain independent exact-overlay paths.
    const auto united=boolean_exact_sections(left,right,ExactSectionOperation::Union);
    const auto xor_set=boolean_exact_sections(left,right,ExactSectionOperation::SymmetricDifference);
    CHECK(boolean_exact_sections(united.section,xor_set.section,ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK(difference.work<united.work);
    for (const auto operation:{ExactSectionOperation::Intersection,ExactSectionOperation::Difference}) {
        const auto measured=boolean_exact_sections(left,right,operation);
        CHECK_THROWS_AS(boolean_exact_sections(left,right,operation,measured.work-1),ResourceLimit);
        CHECK_THROWS_AS(boolean_exact_sections(left,right,operation,0),std::invalid_argument);
        CHECK_THROWS_AS(boolean_exact_sections(left,right,operation,4000001),std::invalid_argument);
        unsigned calls=0;
        CHECK_THROWS_WITH(boolean_exact_sections(left,right,operation,4000000,[&] {
            if (++calls==2) throw std::runtime_error("cancel disjoint bounds");
        }),"cancel disjoint bounds");
    }
    CHECK_THROWS_AS(boolean_exact_sections(left,right,static_cast<ExactSectionOperation>(255)),std::invalid_argument);
}

TEST_CASE("Touching and overlapping enclosures retain exact Boolean fallback", "[ThermoelasticExactSections]")
{
    const auto left=make_exact_section({ExPolygon(Polygon{{0,0},{10,0},{10,10},{0,10}})}).section;
    for (const coord_t x:{9,10}) {
        const auto right=make_exact_section({ExPolygon(Polygon{{x,0},{x+10,0},{x+10,10},{x,10}})}).section;
        const auto intersection=boolean_exact_sections(left,right,ExactSectionOperation::Intersection);
        CHECK(intersection.work>2);
        CHECK(intersection.section.empty()==(x==10));
        const auto difference=boolean_exact_sections(left,right,ExactSectionOperation::Difference);
        CHECK(difference.work>2);
        CHECK(boolean_exact_sections(difference.section,left,ExactSectionOperation::SymmetricDifference).section.empty()==(x==10));
    }
    const auto hole=make_exact_section({ExPolygon(Polygon{{2,2},{8,2},{8,8},{2,8}})}).section;
    const auto ring=boolean_exact_sections(left,hole,ExactSectionOperation::Difference).section;
    const auto in_hole=make_exact_section({ExPolygon(Polygon{{3,3},{7,3},{7,7},{3,7}})}).section;
    const auto disjoint=boolean_exact_sections(ring,in_hole,ExactSectionOperation::Intersection);
    CHECK(disjoint.section.empty());
    CHECK(disjoint.work>2);
    // Rational x=1/3 extrema round outwards to x=1. A right set starting at1
    // touches the enclosure, so it must not use strict metadata separation.
    const auto triangle=make_exact_section({ExPolygon(Polygon{{0,0},{1,0},{0,3}})}).section;
    const auto strip=make_exact_section({ExPolygon(Polygon{{0,2},{2,2},{2,3},{0,3}})}).section;
    const auto rational=boolean_exact_sections(triangle,strip,ExactSectionOperation::Intersection).section;
    const auto right=make_exact_section({ExPolygon(Polygon{{1,2},{2,2},{2,3},{1,3}})}).section;
    const auto separate=boolean_exact_sections(rational,right,ExactSectionOperation::Intersection);
    CHECK(separate.section.empty());
    CHECK(separate.work>2);
}

TEST_CASE("Disjoint metadata cannot bypass the exact Boolean input edge cap", "[ThermoelasticExactSections]")
{
    const auto block=[](coord_t start) {
        ExPolygons polygons;
        for (coord_t i=0;i<256;++i) {
            const coord_t x=start+3*i;
            polygons.emplace_back(Polygon{{x,0},{x+1,0},{x+1,1},{x,1}});
        }
        return make_exact_section(polygons).section;
    };
    const auto left=boolean_exact_sections(block(0),block(1000),ExactSectionOperation::Union).section;
    auto right=boolean_exact_sections(block(10000),block(11000),ExactSectionOperation::Union).section;
    const auto extra=make_exact_section({ExPolygon(Polygon{{12000,0},{12001,0},{12001,1},{12000,1}})}).section;
    right=boolean_exact_sections(right,extra,ExactSectionOperation::Union).section;
    REQUIRE(left.vertex_count()+right.vertex_count()>4096);
    for (const auto operation:{ExactSectionOperation::Intersection,ExactSectionOperation::Difference})
        CHECK_THROWS_WITH(boolean_exact_sections(left,right,operation),
            Catch::Matchers::ContainsSubstring("overlay edge budget"));
}

TEST_CASE("Exact Boolean aliases preserve results and resource checks", "[ThermoelasticExactSections]")
{
    const ExPolygons geometry{ExPolygon(Polygon{{0,0},{5,0},{5,3},{0,3}})};
    const auto original=make_exact_section(geometry).section;
    const auto alias=original;
    const auto separate=make_exact_section(geometry).section;
    const auto operation=GENERATE(ExactSectionOperation::Union,ExactSectionOperation::Intersection,
        ExactSectionOperation::Difference,ExactSectionOperation::SymmetricDifference);
    const auto shared=boolean_exact_sections(original,alias,operation);
    const auto distinct=boolean_exact_sections(original,separate,operation);
    CHECK(shared.work<distinct.work);
    CHECK(shared.section.vertex_count()==distinct.section.vertex_count());
    const bool empty=operation==ExactSectionOperation::Difference || operation==ExactSectionOperation::SymmetricDifference;
    CHECK(boolean_exact_sections(original,alias,operation,1).section.empty()==empty);
    CHECK(shared.section.empty()==empty);
    CHECK(boolean_exact_sections(shared.section,distinct.section,ExactSectionOperation::SymmetricDifference).section.empty());
    if (!empty) {
        CHECK(boolean_exact_sections(shared.section,original,ExactSectionOperation::SymmetricDifference).section.empty());
        REQUIRE(shared.section.bounds().has_value());
        REQUIRE(original.bounds().has_value());
        CHECK(shared.section.bounds()->min==original.bounds()->min);
        CHECK(shared.section.bounds()->max==original.bounds()->max);
    }
    for (const auto& right:{alias,separate}) {
        CHECK_THROWS_AS(boolean_exact_sections(original,right,operation,0),std::invalid_argument);
        CHECK_THROWS_AS(boolean_exact_sections(original,right,operation,4000001),std::invalid_argument);
        CHECK_THROWS_AS(boolean_exact_sections(original,right,static_cast<ExactSectionOperation>(255)),std::invalid_argument);
        unsigned calls=0;
        CHECK_THROWS_AS(boolean_exact_sections(original,right,operation,4000000,
            [&] { if (++calls==2) throw std::runtime_error("cancel"); }),std::runtime_error);
    }
    CHECK_THROWS_AS(boolean_exact_sections(original,separate,operation,distinct.work-1),ResourceLimit);
    CHECK_THROWS_AS(boolean_exact_sections(original,separate,operation,1),ResourceLimit);
    const auto vacant=boolean_exact_sections(ExactSection{},ExactSection{},operation);
    CHECK(vacant.section.empty());
    CHECK(boolean_exact_sections(ExactSection{},ExactSection{},operation,1).section.empty());
    CHECK_THROWS_AS(boolean_exact_sections(ExactSection{},ExactSection{},operation,0),std::invalid_argument);
    const auto left_empty=boolean_exact_sections({},original,operation,1);
    const auto right_empty=boolean_exact_sections(original,{},operation,1);
    CHECK(left_empty.section.empty()==(operation==ExactSectionOperation::Intersection || operation==ExactSectionOperation::Difference));
    CHECK(right_empty.section.empty()==(operation==ExactSectionOperation::Intersection));
    for (const auto& result:{left_empty.section,right_empty.section})
        if (!result.empty())
            CHECK(boolean_exact_sections(result,original,ExactSectionOperation::SymmetricDifference,1).section.empty());
}

TEST_CASE("Exact Boolean equality proves directed holes and disconnected boundaries", "[ThermoelasticExactSections]")
{
    ExPolygon holed(Polygon{{0,0},{10,0},{10,10},{0,10}});
    holed.holes.emplace_back(Polygon{{2,2},{2,8},{8,8},{8,2}});
    const ExPolygon island(Polygon{{20,0},{25,0},{25,5},{20,5}});
    const auto first=make_exact_section({holed,island}).section;
    std::rotate(holed.contour.points.begin(),holed.contour.points.begin()+1,holed.contour.points.end());
    std::rotate(holed.holes.front().points.begin(),holed.holes.front().points.begin()+2,holed.holes.front().points.end());
    const auto second=make_exact_section({island,holed}).section;
    const auto operation=GENERATE(ExactSectionOperation::Union,ExactSectionOperation::Intersection,
        ExactSectionOperation::Difference,ExactSectionOperation::SymmetricDifference);
    const auto result=boolean_exact_sections(first,second,operation);
    const auto edges=first.vertex_count()+second.vertex_count();
    CHECK(result.work<edges*edges+1); // Below the former all-pairs reservation.
    const bool empty=operation==ExactSectionOperation::Difference || operation==ExactSectionOperation::SymmetricDifference;
    CHECK(result.section.empty()==empty);
    if (!empty) CHECK(boolean_exact_sections(result.section,first,ExactSectionOperation::SymmetricDifference,1).section.empty());
    CHECK_THROWS_AS(boolean_exact_sections(first,second,operation,result.work-1),ResourceLimit);
    unsigned calls=0;
    CHECK_THROWS_AS(boolean_exact_sections(first,second,operation,4000000,
        [&] { if (++calls==12) throw std::runtime_error("cancel"); }),std::runtime_error);
    holed.holes.front()=Polygon{{3,2},{3,8},{9,8},{9,2}};
    const auto changed=make_exact_section({island,holed}).section;
    CHECK_FALSE(boolean_exact_sections(first,changed,ExactSectionOperation::SymmetricDifference).section.empty());
}

TEST_CASE("Exact Boolean equality retains the overlay for subdivided boundaries", "[ThermoelasticExactSections]")
{
    const auto first=make_exact_section({ExPolygon(Polygon{{0,0},{10,0},{10,10},{0,10}})}).section;
    const auto second=make_exact_section({ExPolygon(Polygon{{0,0},{5,0},{10,0},{10,10},{0,10}})}).section;
    const auto result=boolean_exact_sections(first,second,ExactSectionOperation::SymmetricDifference);
    CHECK(result.section.empty());
    // If canonicalization retained the additional collinear vertex, directed
    // edge identity cannot prove equality; the ordinary overlay must do so.
    if (first.vertex_count()!=second.vertex_count()) {
        const auto edges=first.vertex_count()+second.vertex_count();
        CHECK(result.work>=edges*edges+1);
    }
    const auto intersection=boolean_exact_sections(first,second,ExactSectionOperation::Intersection);
    CHECK(boolean_exact_sections(first,intersection.section,ExactSectionOperation::SymmetricDifference).section.empty());
}

TEST_CASE("Exact sparse overlays reserve geometric event work", "[ThermoelasticExactSections]")
{
    ExPolygons first,second;
    for (coord_t i=0;i<32;++i) {
        const coord_t x=20*i;
        first.emplace_back(Polygon{{x,0},{x+4,0},{x+4,4},{x,4}});
        second.emplace_back(Polygon{{x+10,0},{x+14,0},{x+14,4},{x+10,4}});
    }
    const auto a=make_exact_section(first).section,b=make_exact_section(second).section;
    const auto operation=GENERATE(ExactSectionOperation::Union,ExactSectionOperation::Intersection,
        ExactSectionOperation::Difference,ExactSectionOperation::SymmetricDifference);
    const auto result=boolean_exact_sections(a,b,operation);
    const auto edges=a.vertex_count()+b.vertex_count();
    CHECK(result.work<edges*edges);
    ExPolygons expected;
    if (operation!=ExactSectionOperation::Intersection) expected=first;
    if (operation==ExactSectionOperation::Union || operation==ExactSectionOperation::SymmetricDifference)
        expected.insert(expected.end(),second.begin(),second.end());
    const std::array<coord_t,3> start{0,0,0},size{700,10,10};
    const auto moments=exact_section_volume_moments(result.section,start,size,0,10);
    const auto reference=SectionIntegration::volume_moments(expected,start,size,0,10);
    for (std::size_t i=0;i<27;++i)
        CHECK_THAT(moments.values[i],Catch::Matchers::WithinAbs(reference[i],std::max(1e-40,std::abs(reference[i])*1e-12)));
    CHECK_THROWS_AS(boolean_exact_sections(a,b,operation,result.work-1),ResourceLimit);
    unsigned calls=0;
    CHECK_THROWS_AS(boolean_exact_sections(a,b,operation,4000000,
        [&] { if (++calls==100) throw std::runtime_error("cancel"); }),std::runtime_error);
}

TEST_CASE("Exact overlay broadphase preserves contact containment and holes", "[ThermoelasticExactSections]")
{
    ExPolygon outer(Polygon{{0,0},{20,0},{20,20},{0,20}});
    const int layout=GENERATE(0,1,2,3);
    if (layout==3) outer.holes.emplace_back(Polygon{{5,5},{5,15},{15,15},{15,5}});
    const ExPolygons first{outer};
    const ExPolygons second{layout==0 ? ExPolygon(Polygon{{20,0},{25,0},{25,20},{20,20}}) :
        layout==1 ? ExPolygon(Polygon{{10,0},{25,0},{25,20},{10,20}}) :
        ExPolygon(Polygon{{7,7},{13,7},{13,13},{7,13}})};
    const auto operation=GENERATE(ExactSectionOperation::Union,ExactSectionOperation::Intersection,
        ExactSectionOperation::Difference,ExactSectionOperation::SymmetricDifference);
    ExPolygons expected;
    if (operation==ExactSectionOperation::Union) {
        auto joined=first; joined.insert(joined.end(),second.begin(),second.end()); expected=union_ex(joined);
    } else if (operation==ExactSectionOperation::Intersection) expected=intersection_ex(first,second);
    else {
        expected=diff_ex(first,second);
        if (operation==ExactSectionOperation::SymmetricDifference) {
            const auto rest=diff_ex(second,first); expected.insert(expected.end(),rest.begin(),rest.end());
        }
    }
    const auto result=boolean_exact_sections(make_exact_section(first).section,make_exact_section(second).section,operation);
    const std::array<coord_t,3> start{0,0,0},size{30,30,10};
    const auto moments=exact_section_volume_moments(result.section,start,size,0,10);
    const auto reference=SectionIntegration::volume_moments(expected,start,size,0,10);
    for (std::size_t i=0;i<27;++i)
        CHECK_THAT(moments.values[i],Catch::Matchers::WithinAbs(reference[i],std::max(1e-40,std::abs(reference[i])*1e-12)));
}

TEST_CASE("Exact dense overlays exhaust event work before construction", "[ThermoelasticExactSections]")
{
    ExPolygons vertical,horizontal;
    for (coord_t i=0;i<192;++i) {
        const coord_t p=2*i;
        vertical.emplace_back(Polygon{{p,0},{p+1,0},{p+1,400},{p,400}});
        horizontal.emplace_back(Polygon{{0,p},{400,p},{400,p+1},{0,p+1}});
    }
    const auto a=make_exact_section(vertical).section,b=make_exact_section(horizontal).section;
    CHECK_THROWS_WITH(boolean_exact_sections(a,b,ExactSectionOperation::Intersection),
        Catch::Matchers::ContainsSubstring("operation=boolean_exact_sections") &&
        Catch::Matchers::ContainsSubstring("limit=4000000"));
}

TEST_CASE("Exact bilinear extrema include rational edge stationary points", "[ThermoelasticExactSections]")
{
    const auto triangle=make_exact_section({ExPolygon(Polygon{{0,0},{5,0},{0,3}})});
    const auto strip=make_exact_section({ExPolygon(Polygon{{1,0},{5,0},{5,3},{1,3}})});
    const auto cut=boolean_exact_sections(triangle.section,strip.section,ExactSectionOperation::Intersection);
    const double sign=GENERATE(-1.,1.);
    const std::array<double,4> values{0.,0.,0.,sign};
    const auto result=exact_section_bilinear_extrema(cut.section,{0,0},{5,3},values);
    const double expected_min=sign<0 ? -.25 : 0.,expected_max=sign>0 ? .25 : 0.;
    CHECK(result.minimum_bounds[0]<=expected_min);
    CHECK(result.minimum_bounds[1]>=expected_min);
    CHECK(result.maximum_bounds[0]<=expected_max);
    CHECK(result.maximum_bounds[1]>=expected_max);
    CHECK_THAT(result.minimum,Catch::Matchers::WithinAbs(expected_min,1e-15));
    CHECK_THAT(result.maximum,Catch::Matchers::WithinAbs(expected_max,1e-15));
    const auto& local=sign<0 ? result.minimum_local : result.maximum_local;
    CHECK_THAT(local[0],Catch::Matchers::WithinAbs(.5,1e-15));
    CHECK_THAT(local[1],Catch::Matchers::WithinAbs(.5,1e-15));
    CHECK_THROWS_AS(exact_section_bilinear_extrema(cut.section,{0,0},{5,3},values,result.work-1),ResourceLimit);
    unsigned calls=0;
    CHECK_THROWS_AS(exact_section_bilinear_extrema(cut.section,{0,0},{5,3},values,4000000,
        [&] { if (++calls==3) throw std::runtime_error("cancel"); }),std::runtime_error);
    CHECK_THROWS_AS(exact_section_bilinear_extrema(cut.section,{0,0},{4,3},values),std::invalid_argument);
    CHECK_THROWS_AS(exact_section_bilinear_extrema(cut.section,{0,0},{0,3},values),std::invalid_argument);
    CHECK_THROWS_AS(exact_section_bilinear_extrema(cut.section,{0,0},{5,3},values,4000001),std::invalid_argument);
    CHECK_THROWS_AS(exact_section_bilinear_extrema({}, {0,0},{5,3},values),std::invalid_argument);
    CHECK_THROWS_AS(exact_section_bilinear_extrema(cut.section,{0,0},{5,3},
        {0.,0.,0.,std::numeric_limits<double>::infinity()}),std::invalid_argument);
}

TEST_CASE("Exact bilinear extrema preserve holes and disconnected native traces", "[ThermoelasticExactSections]")
{
    ExPolygon holed(Polygon{{0,0},{4,0},{4,10},{0,10}});
    holed.holes.emplace_back(Polygon{{1,2},{1,8},{3,8},{3,2}});
    const ExPolygons native{holed,ExPolygon(Polygon{{6,0},{10,0},{10,10},{6,10}})};
    const auto exact=make_exact_section(native);
    const auto values=GENERATE((std::array<double,4>{{-3.,-3.,-3.,-3.}}),
        (std::array<double,4>{{-2.,3.,4.,-7.}}), (std::array<double,4>{{0.,1.,1.,2.}}));
    const auto result=exact_section_bilinear_extrema(exact.section,{0,0},{10,10},values);
    const auto reference=section_trace_extrema(native,{0,0},{10,10},values);
    CHECK_THAT(result.minimum,Catch::Matchers::WithinAbs(reference.minimum,1e-13));
    CHECK_THAT(result.maximum,Catch::Matchers::WithinAbs(reference.maximum,1e-13));
    CHECK(result.minimum_bounds[0]<=result.minimum);
    CHECK(result.minimum_bounds[1]>=result.minimum);
    CHECK(result.maximum_bounds[0]<=result.maximum);
    CHECK(result.maximum_bounds[1]>=result.maximum);
}

TEST_CASE("Exact volume moments preserve rational intersections in every polynomial term", "[ThermoelasticExactSections]")
{
    const ExPolygons triangle{ExPolygon(Polygon{{0,0},{5,0},{0,1}})};
    const ExPolygons strip{ExPolygon(Polygon{{4,0},{6,0},{6,1},{4,1}})};
    const std::array<coord_t,3> start{0,0,0}, size{10,10,10};
    const auto moments=exact_section_volume_moments(triangle,strip,ExactSectionOperation::Intersection,start,size,0,10);
    const auto face=boolean_exact_sections(make_exact_section(triangle).section,make_exact_section(strip).section,
        ExactSectionOperation::Intersection);
    const auto horizontal=exact_section_horizontal_surface_moments(face.section,start,size,3);
    // The intersection has exact vertices(4,0),(5,0),(4,1/5).
    // Map the unit simplex to x/10=2/5+u/10,y/10=v/50;
    // integral u^i v^q du dv = i!q!/(i+q+2)!.
    constexpr unsigned factorial[]{1,1,2,6,24,120,720};
    constexpr unsigned choose[3][3]{{1,0,0},{1,1,0},{1,2,1}};
    for (unsigned r=0;r<3;++r)
        for (unsigned q=0;q<3;++q)
            for (unsigned p=0;p<3;++p) {
                long double xy=0.;
                for (unsigned i=0;i<=p;++i)
                    xy+=choose[p][i]*std::pow(.4L,int(p-i))*std::pow(.1L,int(i))*std::pow(.02L,int(q))*
                        factorial[i]*factorial[q]/factorial[i+q+2];
                const double expected=double(.2L*xy*10.L/(r+1)*1e-27L);
                const auto index=p+3*q+9*r;
                CHECK_THAT(moments.values[index],Catch::Matchers::WithinRel(expected,1e-13));
                CHECK(moments.bounds[index][0]<=moments.values[index]);
                CHECK(moments.bounds[index][1]>=moments.values[index]);
                const double expected_area=double(.2L*xy*std::pow(.3L,int(r))*1e-18L);
                CHECK_THAT(horizontal.values[index],Catch::Matchers::WithinRel(expected_area,1e-13));
                CHECK(horizontal.bounds[index][0]<=horizontal.values[index]);
                CHECK(horizontal.bounds[index][1]>=horizontal.values[index]);
            }
    const auto remainder=exact_section_volume_moments(triangle,strip,ExactSectionOperation::Difference,start,size,0,10);
    const auto whole=exact_section_volume_moments(triangle,{},ExactSectionOperation::Difference,start,size,0,10);
    for (std::size_t i=0;i<27;++i)
        CHECK_THAT(moments.values[i]+remainder.values[i],Catch::Matchers::WithinRel(whole.values[i],1e-13));
    CHECK(moments.work>0);
    CHECK_THROWS_AS(exact_section_volume_moments(triangle,strip,ExactSectionOperation::Intersection,
        start,size,0,10,moments.work-1),ResourceLimit);
    CHECK_THROWS_AS(exact_section_volume_moments(triangle,strip,ExactSectionOperation::Intersection,
        start,size,0,10,4000000,[] { throw std::runtime_error("cancel moments"); }),std::runtime_error);
}

TEST_CASE("Exact volume moments retain translated holes and clip the cell and height", "[ThermoelasticExactSections]")
{
    constexpr coord_t base=200000000;
    ExPolygon polygon(Polygon{{base,base},{base+10,base},{base+10,base+10},{base,base+10}});
    polygon.holes.emplace_back(Polygon{{base+2,base+2},{base+2,base+4},{base+4,base+4},{base+4,base+2}});
    const ExPolygons source{polygon};
    const std::array<coord_t,3> start{base,base,100}, size{10,10,10};
    const auto exact=exact_section_volume_moments(source,{},ExactSectionOperation::Difference,start,size,95,108);
    const auto integer=SectionIntegration::volume_moments(source,start,size,100,108);
    const auto handle=make_exact_section(source);
    const auto horizontal=exact_section_horizontal_surface_moments(handle.section,start,size,103);
    const auto integer_horizontal=SectionIntegration::horizontal_surface_moments(source,start,size,103);
    for (std::size_t i=0;i<27;++i)
        CHECK_THAT(exact.values[i],Catch::Matchers::WithinRel(integer[i],1e-12));
    for (std::size_t i=0;i<27;++i)
        CHECK_THAT(horizontal.values[i],Catch::Matchers::WithinRel(integer_horizontal[i],1e-12));
    const ExPolygons expanded{ExPolygon(Polygon{{base-5,base-5},{base+15,base-5},
        {base+15,base+15},{base-5,base+15}})};
    const auto cell=exact_section_volume_moments(expanded,{},ExactSectionOperation::Difference,start,size,95,115);
    for (unsigned r=0;r<3;++r)
        for (unsigned q=0;q<3;++q)
            for (unsigned p=0;p<3;++p)
                CHECK_THAT(cell.values[p+3*q+9*r],Catch::Matchers::WithinRel(1e-24/double((p+1)*(q+1)*(r+1)),1e-13));
    const auto empty=exact_section_volume_moments(source,{},ExactSectionOperation::Difference,start,size,80,90);
    for (const auto value : empty.values) CHECK_THAT(value,Catch::Matchers::WithinAbs(0.,0.));
    auto invalid_size=size;
    invalid_size[1]=0;
    CHECK_THROWS_AS(exact_section_volume_moments(source,{},ExactSectionOperation::Difference,start,invalid_size,100,110),std::invalid_argument);
    CHECK_THROWS_AS(exact_section_volume_moments(source,{},ExactSectionOperation::Difference,start,size,110,100),std::invalid_argument);
    const ExPolygons excessive(257,polygon);
    CHECK_THROWS_AS(exact_section_volume_moments(excessive,{},ExactSectionOperation::Difference,start,size,100,110),ResourceLimit);
}

TEST_CASE("Exact volume moments preserve rational hole corners", "[ThermoelasticExactSections]")
{
    const ExPolygons outer{ExPolygon(Polygon{{0,0},{10,0},{10,10},{0,10}})};
    const ExPolygons removed{ExPolygon(Polygon{{2,2},{7,2},{2,3}}),
        ExPolygon(Polygon{{6,2},{8,2},{8,3},{6,3}})};
    const std::array<coord_t,3> start{0,0,0}, size{10,10,10};
    const auto moments=exact_section_volume_moments(outer,removed,ExactSectionOperation::Difference,start,size,0,10);
    const auto face=boolean_exact_sections(make_exact_section(outer).section,make_exact_section(removed).section,
        ExactSectionOperation::Difference);
    const auto horizontal=exact_section_horizontal_surface_moments(face.section,start,size,7);
    // The enclosed removed union has a rational boundary vertex(6,11/5).
    // Integrate its triangle + rectangle - overlapping triangle analytically.
    constexpr unsigned factorial[]{1,1,2,6,24,120,720};
    constexpr unsigned choose[3][3]{{1,0,0},{1,1,0},{1,2,1}};
    const auto triangle = [&](unsigned p,unsigned q,long double x,long double y,long double dx,long double dy) {
        long double moment=0.;
        for (unsigned i=0;i<=p;++i)
            for (unsigned j=0;j<=q;++j)
                moment+=choose[p][i]*choose[q][j]*std::pow(x,int(p-i))*std::pow(y,int(q-j))*
                    std::pow(dx,int(i))*std::pow(dy,int(j))*factorial[i]*factorial[j]/factorial[i+j+2];
        return moment*dx*dy;
    };
    for (unsigned r=0;r<3;++r)
        for (unsigned q=0;q<3;++q)
            for (unsigned p=0;p<3;++p) {
                const long double rectangle=(std::pow(.8L,int(p+1))-std::pow(.6L,int(p+1)))/(p+1)*
                    (std::pow(.3L,int(q+1))-std::pow(.2L,int(q+1)))/(q+1);
                const long double hole=triangle(p,q,.2L,.2L,.5L,.1L)+rectangle-triangle(p,q,.6L,.2L,.1L,.02L);
                const double expected=double((1.L/((p+1)*(q+1))-hole)*1e-24L/(r+1));
                CHECK_THAT(moments.values[p+3*q+9*r],Catch::Matchers::WithinRel(expected,1e-13));
                const double expected_area=double((1.L/((p+1)*(q+1))-hole)*1e-16L*std::pow(.7L,int(r)));
                CHECK_THAT(horizontal.values[p+3*q+9*r],Catch::Matchers::WithinRel(expected_area,1e-13));
            }
}

TEST_CASE("Exact horizontal moments clip faces and preserve bounded execution", "[ThermoelasticExactSections]")
{
    const auto section=make_exact_section({ExPolygon(Polygon{{-5,-5},{15,-5},{15,15},{-5,15}})});
    const std::array<coord_t,3> start{0,0,100},size{10,10,10};
    const auto lower=exact_section_horizontal_surface_moments(section.section,start,size,100);
    const auto upper=exact_section_horizontal_surface_moments(section.section,start,size,110);
    for (unsigned r=0;r<3;++r)
        for (unsigned q=0;q<3;++q)
            for (unsigned p=0;p<3;++p) {
                const auto index=p+3*q+9*r;
                const double area=1e-16/double((p+1)*(q+1));
                CHECK_THAT(upper.values[index],Catch::Matchers::WithinRel(area,1e-13));
                if (r==0) CHECK_THAT(lower.values[index],Catch::Matchers::WithinRel(area,1e-13));
                else CHECK_THAT(lower.values[index],Catch::Matchers::WithinAbs(0.,0.));
            }
    const auto outside=exact_section_horizontal_surface_moments(section.section,start,size,111);
    for (const auto value : outside.values) CHECK_THAT(value,Catch::Matchers::WithinAbs(0.,0.));
    CHECK_THROWS_AS(exact_section_horizontal_surface_moments(section.section,start,size,110,upper.work-1),ResourceLimit);
    CHECK_THROWS_AS(exact_section_horizontal_surface_moments(section.section,start,size,110,4000000,
        [] { throw std::runtime_error("cancel face moments"); }),std::runtime_error);
    auto invalid=size;
    invalid[2]=0;
    CHECK_THROWS_AS(exact_section_horizontal_surface_moments(section.section,start,invalid,110),std::invalid_argument);
    CHECK_THROWS_AS(exact_section_horizontal_surface_moments(section.section,start,size,coord_t(10000000001LL)),std::invalid_argument);
}

TEST_CASE("Exact section handles retain rational corners through composed operations", "[ThermoelasticExactSections]")
{
    const auto triangle=make_exact_section({ExPolygon(Polygon{{0,0},{5,0},{0,1}})});
    const auto strip=make_exact_section({ExPolygon(Polygon{{4,0},{6,0},{6,1},{4,1}})});
    const auto clipped=boolean_exact_sections(triangle.section,strip.section,ExactSectionOperation::Intersection);
    const auto cutter=make_exact_section({ExPolygon(Polygon{{4,0},{5,0},{5,1}})});
    const auto remaining=boolean_exact_sections(clipped.section,cutter.section,ExactSectionOperation::Difference);
    REQUIRE_FALSE(remaining.section.empty());
    CHECK(remaining.section.vertex_count()==3);
    const auto bounds=remaining.section.bounds();
    REQUIRE(bounds.has_value());
    CHECK(bounds->min==Point(4,0));
    CHECK(bounds->max==Point(5,1));
    // Remaining vertices are(4,0),(25/6,1/6),(4,1/5). Any native-nm
    // vertex round trip would erase this entire triangle. Exact area is1/60.
    const std::array<coord_t,3> start{0,0,0},size{10,10,10};
    const auto moments=exact_section_volume_moments(remaining.section,start,size,0,10);
    CHECK_THAT(moments.values[0],Catch::Matchers::WithinRel(1e-27/6.,1e-13));
    const auto components=exact_section_components(remaining.section);
    REQUIRE(components.sections.size()==1);
    const auto extracted=exact_section_volume_moments(components.sections.front(),start,size,0,10);
    for (std::size_t i=0;i<27;++i)
        CHECK_THAT(extracted.values[i],Catch::Matchers::WithinRel(moments.values[i],1e-13));
    const auto removed=boolean_exact_sections(clipped.section,cutter.section,ExactSectionOperation::Intersection);
    const auto removed_moments=exact_section_volume_moments(removed.section,start,size,0,10);
    const auto original=exact_section_volume_moments(clipped.section,start,size,0,10);
    for (std::size_t i=0;i<27;++i)
        CHECK_THAT(moments.values[i]+removed_moments.values[i],Catch::Matchers::WithinRel(original.values[i],1e-13));
    // Derived operations never mutate another shared handle's canonical set.
    const auto alias=clipped.section;
    CHECK(alias.vertex_count()==3);
    CHECK_THAT(original.values[0],Catch::Matchers::WithinRel(1e-27,1e-13));
    CHECK_THROWS_AS(boolean_exact_sections(clipped.section,cutter.section,ExactSectionOperation::Difference,remaining.work-1),ResourceLimit);
    CHECK_THROWS_AS(boolean_exact_sections(clipped.section,cutter.section,ExactSectionOperation::Difference,4000000,
        [] { throw std::runtime_error("cancel Boolean"); }),std::runtime_error);
    CHECK_THROWS_AS(exact_section_components(remaining.section,components.work-1),ResourceLimit);
    CHECK_THROWS_AS(exact_section_components(remaining.section,4000000,
        [] { throw std::runtime_error("cancel components"); }),std::runtime_error);
    CHECK_THROWS_AS(exact_section_volume_moments(remaining.section,start,size,0,10,moments.work-1),ResourceLimit);
    CHECK_THROWS_AS(exact_section_volume_moments(remaining.section,start,size,0,10,4000000,
        [] { throw std::runtime_error("cancel handle moments"); }),std::runtime_error);
}

TEST_CASE("Exact section handles expose bounded components without native rounding", "[ThermoelasticExactSections]")
{
    const ExPolygons input{box(0,100),box(200,300)};
    const auto section=make_exact_section(input);
    const auto components=exact_section_components(section.section);
    REQUIRE(components.sections.size()==2);
    CHECK(section.section.vertex_count()==8);
    CHECK(components.sections[0].vertex_count()+components.sections[1].vertex_count()==8);
    CHECK_THROWS_AS(make_exact_section(input,section.work-1),ResourceLimit);
    CHECK_THROWS_AS(make_exact_section(input,4000000,[] { throw std::runtime_error("cancel construction"); }),std::runtime_error);
    CHECK_THROWS_AS(make_exact_section(ExPolygons(257,box(0,100))),ResourceLimit);
    const ExactSection empty;
    CHECK(empty.empty());
    CHECK(empty.vertex_count()==0);
    CHECK_FALSE(empty.bounds().has_value());
    CHECK(exact_section_components(empty).sections.empty());
    CHECK(boolean_exact_sections(section.section,section.section,ExactSectionOperation::SymmetricDifference).section.empty());
}

TEST_CASE("Exact occupied interiors separate point contacts while retaining tangent holes", "[ThermoelasticExactSections]")
{
    const unsigned shape = GENERATE(0u, 1u, 2u, 3u);
    const auto polygon = [](std::initializer_list<Point> points) {
        return make_exact_section({ExPolygon(Polygon(points))}).section;
    };
    ExactSection original;
    if (shape == 0) {
        const auto first = polygon({{0,0},{2,0},{2,2},{0,2}});
        const auto second = polygon({{2,2},{4,2},{4,4},{2,4}});
        original = boolean_exact_sections(first, second, ExactSectionOperation::Union).section;
    } else {
        original = polygon({{0,0},{10,0},{10,10},{0,10}});
        const auto subtract = [&](const ExactSection& hole) {
            original = boolean_exact_sections(original, hole, ExactSectionOperation::Difference).section;
        };
        if (shape == 1) {
            subtract(polygon({{0,5},{3,2},{6,5},{3,8}})); // Hole tangent to outer boundary.
        } else if (shape == 2) {
            subtract(polygon({{1,5},{3,3},{5,5},{3,7}}));
            subtract(polygon({{5,5},{7,3},{9,5},{7,7}})); // Two tangent holes, connected solid.
        } else {
            subtract(polygon({{0,5},{2,3},{5,5},{2,7}}));
            subtract(polygon({{5,5},{8,3},{10,5},{8,7}})); // Void chain separates upper/lower interiors.
        }
    }
    // The new geometry primitive does not silently relax the existing policy.
    CHECK_THROWS_AS(validate_exact_section_topology(original), std::invalid_argument);
    const auto alias = original;
    const auto components = exact_section_interior_components(original);
    const std::size_t expected = shape == 0 || shape == 3 ? 2 : 1;
    REQUIRE(components.sections.size() == expected);
    if (expected == 1)
        CHECK(components.work <= 3 * original.vertex_count() + 4);
    const std::array<coord_t, 3> start{0,0,0}, size{10,10,10};
    const auto original_moments = exact_section_volume_moments(original, start, size, 0, 10);
    const auto scalar = exact_section_volume(original, start, size, 0, 10);
    CHECK_THAT(scalar.value, Catch::Matchers::WithinRel(original_moments.values[0], 1e-13));
    CHECK(scalar.bounds[0] <= original_moments.bounds[0][1]);
    CHECK(scalar.bounds[1] >= original_moments.bounds[0][0]);
    CHECK(scalar.work < original_moments.work);
    const std::array<double, 4> areas{8.,82.,84.,80.};
    CHECK_THAT(original_moments.values[0], Catch::Matchers::WithinRel(areas[shape] * 10e-27, 1e-13));
    std::array<double, 27> summed{};
    ExactSection reconstructed;
    for (std::size_t i = 0; i < components.sections.size(); ++i) {
        const auto& component = components.sections[i];
        REQUIRE_FALSE(component.empty());
        CHECK(exact_section_interior_components(component).sections.size() == 1);
        for (std::size_t j = 0; j < i; ++j)
            CHECK(boolean_exact_sections(component, components.sections[j], ExactSectionOperation::Intersection).section.empty());
        reconstructed = boolean_exact_sections(reconstructed, component, ExactSectionOperation::Union).section;
        const auto moments = exact_section_volume_moments(component, start, size, 0, 10);
        for (std::size_t j = 0; j < summed.size(); ++j) summed[j] += moments.values[j];
    }
    CHECK(boolean_exact_sections(reconstructed, original, ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK(boolean_exact_sections(alias, original, ExactSectionOperation::SymmetricDifference).section.empty());
    for (std::size_t i = 0; i < summed.size(); ++i)
        CHECK_THAT(summed[i], Catch::Matchers::WithinRel(original_moments.values[i], 1e-13));
    // Directed boundaries retain tangent-hole points without fabricating a
    // bridge or an artificial cut surface in the occupied domain.
    const auto boundary = extract_exact_section_boundary(original).boundary;
    auto remaining_boundary = boundary;
    for (const auto& component : components.sections) {
        const auto next = extract_exact_section_boundary(component).boundary;
        CHECK(boolean_exact_section_boundaries(next, boundary, ExactBoundaryOperation::Difference,
            ExactBoundaryOrientation::Any).boundary.empty());
        CHECK(boolean_exact_section_boundaries(next, boundary, ExactBoundaryOperation::Intersection,
            ExactBoundaryOrientation::Opposed).boundary.empty());
        remaining_boundary = boolean_exact_section_boundaries(remaining_boundary, next,
            ExactBoundaryOperation::Difference, ExactBoundaryOrientation::Any).boundary;
    }
    CHECK(remaining_boundary.empty());
    CHECK_THROWS_AS(exact_section_interior_components(original, components.work - 1), ResourceLimit);
    CHECK_THROWS_AS(exact_section_interior_components(original, 0), std::invalid_argument);
    CHECK_THROWS_AS(exact_section_interior_components(original, 4000001), std::invalid_argument);
    unsigned checkpoints = 0;
    CHECK_THROWS_WITH(exact_section_interior_components(original, 4000000, [&] {
        if (++checkpoints == 8) throw std::runtime_error("cancel interior decomposition");
    }), "cancel interior decomposition");
    CHECK(exact_section_interior_components({}).sections.empty());
}

TEST_CASE("Strict exact interiors reconstruct locally without copying unrelated islands", "[ThermoelasticExactSections]")
{
    ExPolygons polygons;
    for (coord_t i=0;i<64;++i) {
        const coord_t x=4*i;
        ExPolygon polygon(Polygon{{x,0},{x+3,0},{x+3,4},{x,4}});
        polygon.holes.emplace_back(Polygon{{x+1,1},{x+1,2},{x+2,2},{x+2,1}});
        polygons.push_back(std::move(polygon));
    }
    const auto original=make_exact_section(polygons).section;
    // The old full512-vertex copy reservation alone needs262657 work before
    // selecting one island. Local8-vertex reconstruction fits this bound.
    const auto components=exact_section_interior_components(original,25000);
    REQUIRE(components.sections.size()==64);
    CHECK(components.work<25000);
    ExactSection reconstructed;
    std::array<double,27> sums{};
    const std::array<coord_t,3> start{0,0,0}, size{256,4,5};
    for (const auto& component:components.sections) {
        CHECK(component.vertex_count()==8);
        CHECK_NOTHROW(validate_exact_section_topology(component));
        reconstructed=boolean_exact_sections(reconstructed,component,ExactSectionOperation::Union).section;
        const auto moments=exact_section_volume_moments(component,start,size,0,5);
        for (std::size_t i=0;i<27;++i) sums[i]+=moments.values[i];
    }
    CHECK(boolean_exact_sections(reconstructed,original,ExactSectionOperation::SymmetricDifference).section.empty());
    const auto expected=exact_section_volume_moments(original,start,size,0,5);
    for (std::size_t i=0;i<27;++i)
        CHECK_THAT(sums[i],Catch::Matchers::WithinRel(expected.values[i],1e-13));
    CHECK_THROWS_AS(exact_section_interior_components(original,components.work-1),ResourceLimit);
    unsigned calls=0;
    CHECK_THROWS_WITH(exact_section_interior_components(original,25000,[&] {
        if (++calls==100) throw ResourceLimit("cancel local interior reconstruction");
    }),"cancel local interior reconstruction");
}

TEST_CASE("Mixed strict and tangent interiors retain exact fallback geometry", "[ThermoelasticExactSections]")
{
    const auto outer=make_exact_section({ExPolygon(Polygon{{0,0},{10,0},{10,10},{0,10}})}).section;
    const auto hole=make_exact_section({ExPolygon(Polygon{{0,5},{3,2},{6,5},{3,8}})}).section;
    const auto tangent=boolean_exact_sections(outer,hole,ExactSectionOperation::Difference).section;
    const auto detached=make_exact_section({ExPolygon(Polygon{{12,0},{14,0},{14,2},{12,2}})}).section;
    const auto original=boolean_exact_sections(tangent,detached,ExactSectionOperation::Union).section;
    const auto components=exact_section_interior_components(original);
    REQUIRE(components.sections.size()==2);
    ExactSection reconstructed;
    std::array<double,27> sums{};
    const std::array<coord_t,3> start{0,0,0}, size{14,10,5};
    for (const auto& component:components.sections) {
        reconstructed=boolean_exact_sections(reconstructed,component,ExactSectionOperation::Union).section;
        const auto moments=exact_section_volume_moments(component,start,size,0,5);
        for (std::size_t i=0;i<27;++i) sums[i]+=moments.values[i];
    }
    CHECK(boolean_exact_sections(reconstructed,original,ExactSectionOperation::SymmetricDifference).section.empty());
    const auto expected=exact_section_volume_moments(original,start,size,0,5);
    for (std::size_t i=0;i<27;++i)
        CHECK_THAT(sums[i],Catch::Matchers::WithinRel(expected.values[i],1e-13));
    CHECK_THROWS_WITH(exact_section_interior_components(original,components.work-1),
        Catch::Matchers::ContainsSubstring("occupied_faces=2") && Catch::Matchers::ContainsSubstring("strict_pieces=false"));
}

TEST_CASE("Exact polynomial moments match an independent rational simplex oracle", "[ThermoelasticExactSections]")
{
    const coord_t offset = GENERATE(coord_t(0), coord_t(-30));
    const bool holed = GENERATE(false, true);
    const bool contained = GENERATE(false, true);
    ExPolygon polygon(Polygon{{offset,offset},{offset+20,offset},{offset,offset+12}});
    if (holed) polygon.holes.emplace_back(Polygon{{offset+3,offset+1},{offset+3,offset+2},
        {offset+5,offset+2},{offset+5,offset+1}});
    const auto original = make_exact_section({polygon}).section;
    auto section = original;
    const std::array<coord_t,3> start{offset+1,offset,-4}, size{6,12,10};
    if (contained) {
        const auto cell = make_exact_section({ExPolygon(Polygon{{offset+1,offset},{offset+7,offset},
            {offset+7,offset+12},{offset+1,offset+12}})}).section;
        section = boolean_exact_sections(original,cell,ExactSectionOperation::Intersection).section;
    }
    const auto volume = exact_section_volume_moments(section,start,size,-3,3);
    const auto surface = exact_section_horizontal_surface_moments(section,start,size,3);
    if (contained) {
        // The same rational occupied set now fits the original normalization
        // cell. Skipping the redundant overlay must save work, not change any
        // of its 27 moments or the immutable input handle.
        const auto clipped_volume = exact_section_volume_moments(original,start,size,-3,3);
        const auto clipped_surface = exact_section_horizontal_surface_moments(original,start,size,3);
        CHECK(volume.work < clipped_volume.work);
        CHECK(surface.work < clipped_surface.work);
        for (std::size_t i = 0; i < 27; ++i) {
            CHECK_THAT(volume.values[i],Catch::Matchers::WithinRel(clipped_volume.values[i],1e-13));
            CHECK_THAT(surface.values[i],Catch::Matchers::WithinRel(clipped_surface.values[i],1e-13));
            CHECK(volume.bounds[i][0] <= clipped_volume.bounds[i][1]);
            CHECK(volume.bounds[i][1] >= clipped_volume.bounds[i][0]);
        }
    }
    using Rational = boost::multiprecision::cpp_rational;
    using RationalPoint = std::array<Rational,2>;
    const auto point = [](Rational x, Rational y) -> RationalPoint { return {(x-1)/6,y/12}; };
    // Analytic intersection of y=12-3*x/5 with cell faces x=1 and x=7.
    std::vector<std::vector<RationalPoint>> rings{{point(1,0),point(7,0),
        point(7,Rational(39)/5),point(1,Rational(57)/5)}};
    if (holed) rings.push_back({point(3,1),point(3,2),point(5,2),point(5,1)});
    const auto power = [](Rational x, unsigned n) {
        Rational result(1);
        while (n--) result *= x;
        return result;
    };
    constexpr unsigned factorial[]{1,1,2,6,24,120,720};
    constexpr unsigned binomial[3][3]{{1,0,0},{1,1,0},{1,2,1}};
    const Rational volume_scale = Rational(720)/Rational("1000000000000000000000000000");
    const Rational area_scale = Rational(72)/Rational("1000000000000000000");
    const auto compare = [](const ExactSectionMoments& actual, std::size_t index, const Rational& expected) {
        CHECK_THAT(actual.values[index], Catch::Matchers::WithinRel(expected.convert_to<double>(),1e-13));
        CHECK(Rational(actual.bounds[index][0]) <= expected);
        CHECK(Rational(actual.bounds[index][1]) >= expected);
    };
    for (unsigned q = 0; q < 3; ++q)
        for (unsigned p = 0; p < 3; ++p) {
            Rational xy(0);
            // Signed triangles (origin,a,b), not the production Green-edge
            // polynomial expansion. Simplex integrals are exact factorial
            // ratios; clockwise hole triangles subtract without quadrature.
            for (const auto& ring : rings)
                for (std::size_t edge = 0; edge < ring.size(); ++edge) {
                    const auto& a = ring[edge];
                    const auto& b = ring[(edge+1)%ring.size()];
                    const Rational determinant = a[0]*b[1]-a[1]*b[0];
                    for (unsigned i = 0; i <= p; ++i)
                        for (unsigned j = 0; j <= q; ++j)
                            xy += determinant * binomial[p][i] * binomial[q][j] *
                                power(a[0],i)*power(b[0],p-i)*power(a[1],j)*power(b[1],q-j)*
                                factorial[i+j]*factorial[p+q-i-j]/factorial[p+q+2];
                }
            for (unsigned r = 0; r < 3; ++r) {
                const auto index = p+3*q+9*r;
                const Rational z = (power(Rational(7)/10,r+1)-power(Rational(1)/10,r+1))/(r+1);
                compare(volume,index,xy*z*volume_scale);
                compare(surface,index,xy*power(Rational(7)/10,r)*area_scale);
            }
        }
    CHECK_THROWS_AS(exact_section_volume_moments(section,start,size,-3,3,volume.work-1), ResourceLimit);
    unsigned calls = 0;
    CHECK_THROWS_WITH(exact_section_horizontal_surface_moments(section,start,size,3,4000000,[&] {
        if (++calls == 20) throw std::runtime_error("cancel cached polynomial powers");
    }), "cancel cached polynomial powers");
}

TEST_CASE("Exact scalar measures preserve rational clipping and bounded SI conversion", "[ThermoelasticExactSections]")
{
    const bool clipped = GENERATE(false, true);
    const coord_t offset = GENERATE(coord_t(0), coord_t(-20));
    const auto section = make_exact_section({ExPolygon(Polygon{{offset,offset},
        {offset+10,offset},{offset,offset+3}})}).section;
    const std::array<coord_t,3> start{offset+(clipped ? 1 : 0),offset,2};
    const std::array<coord_t,3> size{clipped ? 5 : 10,3,7};
    const auto volume = exact_section_volume(section, start, size, 0, 20);
    const auto area = exact_section_horizontal_area(section, start, size, 2);
    const auto full_volume = exact_section_volume_moments(section, start, size, 0, 20);
    const auto full_area = exact_section_horizontal_surface_moments(section, start, size, 2);
    const auto compare = [](const ExactSectionMeasure& scalar, const ExactSectionMoments& tensor) {
        CHECK_THAT(scalar.value, Catch::Matchers::WithinRel(tensor.values[0], 1e-13));
        CHECK(scalar.bounds[0] <= scalar.value);
        CHECK(scalar.bounds[1] >= scalar.value);
        CHECK(scalar.bounds[0] <= tensor.bounds[0][1]);
        CHECK(scalar.bounds[1] >= tensor.bounds[0][0]);
        CHECK(scalar.work < tensor.work);
    };
    compare(volume, full_volume);
    compare(area, full_area);
    using Decimal = boost::multiprecision::cpp_dec_float_50;
    const Decimal exact_area = (clipped ? Decimal(39)/4 : Decimal(15)) / Decimal("1e18");
    const Decimal exact_volume = exact_area * 7 / Decimal("1e9");
    CHECK(Decimal(area.bounds[0]) <= exact_area);
    CHECK(Decimal(area.bounds[1]) >= exact_area);
    CHECK(Decimal(volume.bounds[0]) <= exact_volume);
    CHECK(Decimal(volume.bounds[1]) >= exact_volume);
    CHECK_THAT(exact_section_horizontal_area(section,start,size,9).value,
        Catch::Matchers::WithinRel(area.value,1e-13));
    CHECK(exact_section_horizontal_area(section,start,size,1).value == 0.);
    CHECK(exact_section_volume(section,start,size,10,20).value == 0.);
    CHECK(exact_section_volume({},start,size,0,20).value == 0.);
    CHECK_THROWS_AS(exact_section_volume(section,start,size,0,20,volume.work-1), ResourceLimit);
    CHECK_THROWS_AS(exact_section_horizontal_area(section,start,size,2,area.work-1), ResourceLimit);
    CHECK_THROWS_AS(exact_section_volume(section,start,size,0,0), std::invalid_argument);
    CHECK_THROWS_AS(exact_section_horizontal_area(section,start,{0,3,7},2), std::invalid_argument);
    CHECK_THROWS_AS(exact_section_volume(section,start,size,0,20,0), std::invalid_argument);
    CHECK_THROWS_AS(exact_section_horizontal_area(section,start,size,2,4000001), std::invalid_argument);
    unsigned calls = 0;
    CHECK_THROWS_WITH(exact_section_volume(section,start,size,0,20,4000000,[&] {
        if (++calls == 3) throw std::runtime_error("cancel scalar measure");
    }), "cancel scalar measure");
}

TEST_CASE("Single exact interiors validate linearly and retain immutable canonical geometry", "[ThermoelasticExactSections]")
{
    const auto original = make_exact_section({box(0,100)}).section;
    const auto components = exact_section_interior_components(original);
    REQUIRE(components.sections.size() == 1);
    CHECK(components.work <= 3 * original.vertex_count() + 4);
    // The previous polygon reconstruction pays a quadratic reservation even
    // for one simple component. The fast path pays only vertex/face/edge visits.
    CHECK(components.work < exact_section_components(original).work);
    const auto saved = components.sections.front();
    const auto cutter = make_exact_section({box(0,50)}).section;
    const auto changed = boolean_exact_sections(saved, cutter, ExactSectionOperation::Difference).section;
    CHECK_FALSE(boolean_exact_sections(changed, original, ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK(boolean_exact_sections(saved, original, ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK_THROWS_AS(exact_section_interior_components(original, components.work - 1), ResourceLimit);
    unsigned calls = 0;
    CHECK_THROWS_WITH(exact_section_interior_components(original, 4000000, [&] {
        if (++calls == 4) throw std::runtime_error("cancel single interior validation");
    }), "cancel single interior validation");
}

TEST_CASE("Exact boundary contacts retain rational tile face junctions", "[ThermoelasticExactSections]")
{
    const auto source=make_exact_section({ExPolygon(Polygon{{0,0},{5,0},{0,1}})});
    const auto left_box=make_exact_section({ExPolygon(Polygon{{0,0},{4,0},{4,1},{0,1}})});
    const auto right_box=make_exact_section({ExPolygon(Polygon{{4,0},{6,0},{6,1},{4,1}})});
    const auto left=boolean_exact_sections(source.section,left_box.section,ExactSectionOperation::Intersection);
    const auto right=boolean_exact_sections(source.section,right_box.section,ExactSectionOperation::Intersection);
    const auto left_boundary=extract_exact_section_boundary(left.section);
    const auto left_trace=exact_boundary_on_native_plane(left_boundary.boundary,0,4);
    const auto right_trace=exact_boundary_on_native_plane(extract_exact_section_boundary(right.section).boundary,0,4);
    REQUIRE(left_trace.boundary.segment_count()==1);
    REQUIRE(right_trace.boundary.segment_count()==1);
    CHECK(left_trace.boundary.bounds()->min==Point(4,0));
    CHECK(left_trace.boundary.bounds()->max==Point(4,1));
    // The shared vertical edge is only1/5nm high. Split its right side at
    // y=1/10nm with an integer-source oblique cutter, preserving a T junction.
    const auto cutter=make_exact_section({ExPolygon(Polygon{{3,0},{13,1},{13,10},{3,10}})});
    const auto upper=boolean_exact_sections(right.section,cutter.section,ExactSectionOperation::Intersection);
    const auto lower=boolean_exact_sections(right.section,cutter.section,ExactSectionOperation::Difference);
    const auto upper_trace=exact_boundary_on_native_plane(extract_exact_section_boundary(upper.section).boundary,0,4);
    const auto lower_trace=exact_boundary_on_native_plane(extract_exact_section_boundary(lower.section).boundary,0,4);
    const auto contact=boolean_exact_section_boundaries(left_trace.boundary,right_trace.boundary,
        ExactBoundaryOperation::Intersection,ExactBoundaryOrientation::Opposed);
    CHECK(contact.boundary.segment_count()==1);
    CHECK_FALSE(boolean_exact_section_boundaries(left_trace.boundary,upper_trace.boundary,
        ExactBoundaryOperation::Intersection,ExactBoundaryOrientation::Opposed).boundary.empty());
    const auto remainder=boolean_exact_section_boundaries(left_trace.boundary,lower_trace.boundary,
        ExactBoundaryOperation::Difference,ExactBoundaryOrientation::Opposed);
    CHECK_FALSE(remainder.boundary.empty());
    CHECK(boolean_exact_section_boundaries(remainder.boundary,upper_trace.boundary,
        ExactBoundaryOperation::Difference,ExactBoundaryOrientation::Opposed).boundary.empty());
    CHECK(boolean_exact_section_boundaries(left_trace.boundary,left_trace.boundary,
        ExactBoundaryOperation::Intersection,ExactBoundaryOrientation::Opposed).boundary.empty());
    CHECK(boolean_exact_section_boundaries(left_trace.boundary,left_trace.boundary,
        ExactBoundaryOperation::Difference,ExactBoundaryOrientation::Any).boundary.empty());
    CHECK_THROWS_AS(extract_exact_section_boundary(left.section,left_boundary.work-1),ResourceLimit);
    CHECK_THROWS_AS(exact_boundary_on_native_plane(left_boundary.boundary,0,4,left_trace.work-1),ResourceLimit);
    CHECK_THROWS_AS(boolean_exact_section_boundaries(left_trace.boundary,right_trace.boundary,
        ExactBoundaryOperation::Intersection,ExactBoundaryOrientation::Opposed,contact.work-1),ResourceLimit);
    CHECK_THROWS_AS(boolean_exact_section_boundaries(left_trace.boundary,right_trace.boundary,
        ExactBoundaryOperation::Intersection,ExactBoundaryOrientation::Opposed,4000000,
        [] { throw std::runtime_error("cancel exact boundary"); }),std::runtime_error);
    CHECK_THROWS_AS(extract_exact_section_boundary(left.section,4000000,
        [] { throw std::runtime_error("cancel boundary extraction"); }),std::runtime_error);
    CHECK_THROWS_AS(exact_boundary_on_native_plane(left_boundary.boundary,0,4,4000000,
        [] { throw std::runtime_error("cancel plane trace"); }),std::runtime_error);
    CHECK_THROWS_AS(exact_boundary_on_native_plane(left_boundary.boundary,2,4),std::invalid_argument);
}

TEST_CASE("Exact traces compare native intervals without rounding rational endpoints", "[ThermoelasticExactSections]")
{
    const auto source=make_exact_section({ExPolygon(Polygon{{0,0},{5,0},{0,1}})});
    const auto clip=make_exact_section({ExPolygon(Polygon{{0,0},{4,0},{4,1},{0,1}})});
    const auto clipped=boolean_exact_sections(source.section,clip.section,ExactSectionOperation::Intersection);
    const auto boundary=extract_exact_section_boundary(clipped.section).boundary;
    // The x=4 face has positive length1/5nm, although integer endpoint
    // rounding would collapse it. Its outward enclosure is not its trace.
    const auto overlap=exact_boundary_overlaps_native_interval(boundary,0,4,0,1);
    CHECK_FALSE(overlap.empty);
    CHECK(exact_boundary_overlaps_native_interval(boundary,0,4,-1,0).empty);
    CHECK(exact_boundary_overlaps_native_interval(boundary,0,4,1,2).empty);
    CHECK(exact_boundary_overlaps_native_interval(boundary,0,4,0,0).empty);
    CHECK(exact_boundary_overlaps_native_interval(boundary,0,3,0,1).empty);
    CHECK_FALSE(exact_boundary_overlaps_native_interval(boundary,1,0,1,3).empty);
    CHECK(exact_boundary_overlaps_native_interval(boundary,1,0,4,5).empty);
    CHECK(exact_boundary_overlaps_native_interval(ExactSectionBoundary{},0,4,0,1).empty);
    CHECK_THROWS_AS(exact_boundary_overlaps_native_interval(boundary,0,4,0,1,overlap.work-1),ResourceLimit);
    CHECK_THROWS_AS(exact_boundary_overlaps_native_interval(boundary,2,4,0,1),std::invalid_argument);
    CHECK_THROWS_AS(exact_boundary_overlaps_native_interval(boundary,0,4,1,0),std::invalid_argument);
    CHECK_THROWS_AS(exact_boundary_overlaps_native_interval(boundary,0,coord_t(10000000001LL),0,1),std::invalid_argument);
    std::size_t checkpoints=0;
    CHECK_THROWS_AS(exact_boundary_overlaps_native_interval(boundary,0,3,0,1,4000000,[&] {
        if (++checkpoints==2) throw std::runtime_error("cancel mixed trace");
    }),std::runtime_error);
    CHECK(checkpoints==2);
}

TEST_CASE("Exact boundaries orient holes and reject point only contacts", "[ThermoelasticExactSections]")
{
    ExPolygon holed(Polygon{{0,0},{10,0},{10,10},{0,10}});
    holed.holes.emplace_back(Polygon{{2,2},{2,4},{4,4},{4,2}});
    const auto surrounding=extract_exact_section_boundary(make_exact_section({holed}).section);
    const auto insert=extract_exact_section_boundary(make_exact_section({ExPolygon(Polygon{{2,2},{4,2},{4,4},{2,4}})}).section);
    CHECK(surrounding.boundary.segment_count()==8);
    const auto contact=boolean_exact_section_boundaries(surrounding.boundary,insert.boundary,
        ExactBoundaryOperation::Intersection,ExactBoundaryOrientation::Opposed);
    CHECK(contact.boundary.segment_count()==4);
    const auto exposed=boolean_exact_section_boundaries(surrounding.boundary,insert.boundary,
        ExactBoundaryOperation::Difference,ExactBoundaryOrientation::Opposed);
    CHECK(exposed.boundary.segment_count()==4);
    const auto touching=extract_exact_section_boundary(make_exact_section({ExPolygon(Polygon{{10,10},{12,10},{12,12},{10,12}})}).section);
    CHECK(boolean_exact_section_boundaries(surrounding.boundary,touching.boundary,
        ExactBoundaryOperation::Intersection,ExactBoundaryOrientation::Any).boundary.empty());
    CHECK(exact_boundary_on_native_plane(surrounding.boundary,0,5).boundary.empty());
    const ExactSectionBoundary empty;
    CHECK(empty.empty());
    CHECK(empty.segment_count()==0);
    CHECK_FALSE(empty.bounds().has_value());
}

TEST_CASE("Exact surface bounds survive cancellation at translated rational endpoints", "[ThermoelasticExactSections]")
{
    const coord_t n=GENERATE(coord_t(100000),coord_t(10000000),coord_t(999999937),coord_t(1000000000));
    const auto triangle=make_exact_section({ExPolygon(Polygon{{0,n},{n,n+1},{0,n+1}})});
    const auto strip=make_exact_section({ExPolygon(Polygon{{n-1,n},{n,n},{n,n+1},{n-1,n+1}})});
    const auto cut=boolean_exact_sections(triangle.section,strip.section,ExactSectionOperation::Intersection);
    const auto boundary=extract_exact_section_boundary(cut.section);
    REQUIRE(boundary.boundary.segment_count()==3);
    const auto split=split_exact_boundary_segments(boundary.boundary);
    REQUIRE(split.segments.size()==3);
    using Oracle=boost::multiprecision::cpp_dec_float_100;
    const Oracle inverse=Oracle(1)/n;
    const Oracle diagonal_length=boost::multiprecision::sqrt(1+inverse*inverse);
    for (const auto& segment:split.segments) {
        CHECK_THAT(segment.outward_normal[0]*segment.outward_normal[0]+
            segment.outward_normal[1]*segment.outward_normal[1],Catch::Matchers::WithinAbs(1.,1e-14));
        const std::array<Oracle,2> normal=segment.outward_normal[0]<-.5 ?
            std::array<Oracle,2>{{-1,0}} : (segment.outward_normal[1]>.5 ?
                std::array<Oracle,2>{{0,1}} : std::array<Oracle,2>{{inverse/diagonal_length,-1/diagonal_length}});
        for (unsigned axis=0;axis<2;++axis) {
            const auto& bounds=segment.outward_normal_bounds[axis];
            CHECK(Oracle(bounds[0])<=normal[axis]);
            CHECK(Oracle(bounds[1])>=normal[axis]);
            CHECK(std::nextafter(bounds[0],bounds[1])>=bounds[1]);
        }
    }
    const auto moments=exact_section_vertical_surface_moments(boundary.boundary,{n-1,n,0},{1,1,10},0,10,
        ExactBoundaryCellOwnership::AllIncident);
    const Oracle area=(1+inverse+boost::multiprecision::sqrt(1+inverse*inverse))*10/
        boost::multiprecision::pow(Oracle(10),18);
    CHECK(Oracle(moments.bounds[0][0])<=area);
    CHECK(Oracle(moments.bounds[0][1])>=area);
    CHECK_THAT(moments.values[0],Catch::Matchers::WithinRel(area.convert_to<double>(),1e-13));
    CHECK_THROWS_AS(exact_section_vertical_surface_moments(boundary.boundary,{n-1,n,0},{1,1,10},0,10,
        ExactBoundaryCellOwnership::AllIncident,moments.work-1),ResourceLimit);
    unsigned checkpoints=0;
    CHECK_THROWS_AS(exact_section_vertical_surface_moments(boundary.boundary,{n-1,n,0},{1,1,10},0,10,
        ExactBoundaryCellOwnership::AllIncident,4000000,[&] { if (++checkpoints==20) throw std::runtime_error("cancel"); }),
        std::runtime_error);
}

TEST_CASE("Exact vertical moments retain tiny polynomial coefficients", "[ThermoelasticExactSections]")
{
    const auto whole=extract_exact_section_boundary(make_exact_section({ExPolygon(Polygon{{0,0},{1,1},{0,1}})}).section);
    const auto axis=exact_boundary_on_native_plane(whole.boundary,0,0);
    const auto horizontal=exact_boundary_on_native_plane(whole.boundary,1,1);
    auto diagonal=boolean_exact_section_boundaries(whole.boundary,axis.boundary,
        ExactBoundaryOperation::Difference,ExactBoundaryOrientation::Any);
    diagonal=boolean_exact_section_boundaries(diagonal.boundary,horizontal.boundary,
        ExactBoundaryOperation::Difference,ExactBoundaryOrientation::Any);
    REQUIRE(diagonal.boundary.segment_count()==1);
    constexpr coord_t extent=10000000000LL;
    const std::array<coord_t,3> start{0,0,0},size{extent,extent,extent};
    const auto moments=exact_section_vertical_surface_moments(diagonal.boundary,start,size,0,1,
        ExactBoundaryCellOwnership::MaterialSide);
    using Oracle=boost::multiprecision::cpp_dec_float_100;
    for (unsigned r=0;r<3;++r)
        for (unsigned q=0;q<3;++q)
            for (unsigned p=0;p<3;++p) {
                const auto i=p+3*q+9*r;
                const Oracle expected=boost::multiprecision::sqrt(Oracle(2))/
                    boost::multiprecision::pow(Oracle(10),18)/
                    boost::multiprecision::pow(Oracle(extent),int(p+q+r))/(p+q+1)/(r+1);
                CHECK(Oracle(moments.bounds[i][0])<=expected);
                CHECK(Oracle(moments.bounds[i][1])>=expected);
                CHECK_THAT(moments.values[i],Catch::Matchers::WithinRel(expected.convert_to<double>(),1e-13));
            }
    CHECK_THROWS_AS(exact_section_vertical_surface_moments(diagonal.boundary,start,size,0,1,
        ExactBoundaryCellOwnership::MaterialSide,moments.work-1),ResourceLimit);
    unsigned checkpoints=0;
    CHECK_THROWS_AS(exact_section_vertical_surface_moments(diagonal.boundary,start,size,0,1,
        ExactBoundaryCellOwnership::MaterialSide,4000000,[&] { if (++checkpoints==20) throw std::runtime_error("cancel"); }),
        std::runtime_error);
}

TEST_CASE("Exact vertical moments bound oblique lengths and every polynomial term", "[ThermoelasticExactSections]")
{
    const std::array<coord_t,3> start{0,0,0},size{10,10,10};
    for (const Point tip : {Point(3,4),Point(1,1)}) {
        const auto whole=extract_exact_section_boundary(make_exact_section({ExPolygon(Polygon{{0,0},tip,{0,tip.y()}})}).section);
        const auto axis=exact_boundary_on_native_plane(whole.boundary,0,0);
        const auto horizontal=exact_boundary_on_native_plane(whole.boundary,1,tip.y());
        auto oblique=boolean_exact_section_boundaries(whole.boundary,axis.boundary,ExactBoundaryOperation::Difference,ExactBoundaryOrientation::Any);
        oblique=boolean_exact_section_boundaries(oblique.boundary,horizontal.boundary,ExactBoundaryOperation::Difference,ExactBoundaryOrientation::Any);
        REQUIRE(oblique.boundary.segment_count()==1);
        const auto moments=exact_section_vertical_surface_moments(oblique.boundary,start,size,2,8,ExactBoundaryCellOwnership::MaterialSide);
        // On macOS ARM long double is only binary64. Use an independent
        // high-precision oracle with integer-built fractions, not rounded
        // .2/.8 literals, when checking the certified one-ULP bounds.
        using Oracle=boost::multiprecision::cpp_dec_float_100;
        const Oracle length=boost::multiprecision::sqrt(Oracle(tip.x()*tip.x()+tip.y()*tip.y()));
        const auto normals=split_exact_boundary_segments(oblique.boundary);
        REQUIRE(normals.segments.size()==1);
        const std::array<Oracle,2> expected_normal{{Oracle(tip.y())/length,-Oracle(tip.x())/length}};
        for (std::size_t axis_index=0;axis_index<2;++axis_index) {
            const auto& normal=normals.segments.front();
            CHECK_THAT(normal.outward_normal[axis_index],Catch::Matchers::WithinAbs(expected_normal[axis_index].convert_to<double>(),1e-14));
            CHECK(Oracle(normal.outward_normal_bounds[axis_index][0])<=expected_normal[axis_index]);
            CHECK(Oracle(normal.outward_normal_bounds[axis_index][1])>=expected_normal[axis_index]);
        }
        for (unsigned r=0;r<3;++r)
            for (unsigned q=0;q<3;++q)
                for (unsigned p=0;p<3;++p) {
                    const auto index=p+3*q+9*r;
                    const Oracle expected=length*10/boost::multiprecision::pow(Oracle(10),18)*
                        boost::multiprecision::pow(Oracle(tip.x())/10,int(p))*
                        boost::multiprecision::pow(Oracle(tip.y())/10,int(q))/(p+q+1)*
                        (boost::multiprecision::pow(Oracle(8)/10,int(r+1))-
                         boost::multiprecision::pow(Oracle(2)/10,int(r+1)))/(r+1);
                    CHECK_THAT(moments.values[index],Catch::Matchers::WithinRel(expected.convert_to<double>(),1e-13));
                    CHECK(Oracle(moments.bounds[index][0])<=expected);
                    CHECK(Oracle(moments.bounds[index][1])>=expected);
                }
        const auto integer=SectionIntegration::vertical_segment_moments(Eigen::Vector2d(0.,double(tip.y())/10.),
            Eigen::Vector2d(0.,0.),Eigen::Vector3d(1e-8,1e-8,1e-8),.2,.8);
        const auto axial=exact_section_vertical_surface_moments(axis.boundary,start,size,2,8,ExactBoundaryCellOwnership::MaterialSide);
        for (std::size_t i=0;i<27;++i)
            CHECK_THAT(axial.values[i],Catch::Matchers::WithinAbs(integer[i],1e-31));
        if (tip==Point(3,4)) {
            // Clip the oblique edge at interior Y=2: its new endpoint is
            // exactly(3/2,2), its length is5/2nm, and normalizedY runs0..1.
            const std::array<coord_t,3> clipped_size{10,2,10};
            const auto clipped=exact_section_vertical_surface_moments(oblique.boundary,start,clipped_size,
                2,8,ExactBoundaryCellOwnership::MaterialSide);
            for (unsigned r=0;r<3;++r)
                for (unsigned q=0;q<3;++q)
                    for (unsigned p=0;p<3;++p) {
                        const long double expected=2.5L*10e-18L*std::pow(.15L,int(p))/(p+q+1)*
                            (std::pow(.8L,int(r+1))-std::pow(.2L,int(r+1)))/(r+1);
                        CHECK_THAT(clipped.values[p+3*q+9*r],Catch::Matchers::WithinRel(double(expected),1e-13));
                    }
        }
        CHECK_THROWS_AS(exact_section_vertical_surface_moments(oblique.boundary,start,size,2,8,
            ExactBoundaryCellOwnership::MaterialSide,moments.work-1),ResourceLimit);
        CHECK_THROWS_AS(exact_section_vertical_surface_moments(oblique.boundary,start,size,2,8,
            ExactBoundaryCellOwnership::MaterialSide,4000000,[] { throw std::runtime_error("cancel vertical moments"); }),std::runtime_error);
    }
}

TEST_CASE("Exact vertical moments assign cell faces to material and clip rational spans", "[ThermoelasticExactSections]")
{
    const auto source=make_exact_section({ExPolygon(Polygon{{0,0},{5,0},{0,1}})});
    const auto left_box=make_exact_section({ExPolygon(Polygon{{0,0},{4,0},{4,1},{0,1}})});
    const auto left=boolean_exact_sections(source.section,left_box.section,ExactSectionOperation::Intersection);
    const auto edge=exact_boundary_on_native_plane(extract_exact_section_boundary(left.section).boundary,0,4);
    const std::array<coord_t,3> a{0,0,0},b{4,0,0},size{4,10,10};
    const auto owner=exact_section_vertical_surface_moments(edge.boundary,a,size,-2,12,ExactBoundaryCellOwnership::MaterialSide);
    const auto neighbor=exact_section_vertical_surface_moments(edge.boundary,b,size,0,10,ExactBoundaryCellOwnership::MaterialSide);
    CHECK_THAT(owner.values[0],Catch::Matchers::WithinRel(2e-18,1e-13));
    for (const auto value : neighbor.values) CHECK_THAT(value,Catch::Matchers::WithinAbs(0.,0.));
    const auto both=exact_section_vertical_surface_moments(edge.boundary,b,size,0,10,ExactBoundaryCellOwnership::AllIncident);
    CHECK_THAT(both.values[0],Catch::Matchers::WithinRel(owner.values[0],1e-13));
    const auto cutter=make_exact_section({ExPolygon(Polygon{{3,0},{13,1},{13,10},{3,10}})});
    const auto upper=boolean_exact_sections(left.section,cutter.section,ExactSectionOperation::Intersection);
    const auto lower=boolean_exact_sections(left.section,cutter.section,ExactSectionOperation::Difference);
    const auto first=exact_boundary_on_native_plane(extract_exact_section_boundary(upper.section).boundary,0,4);
    const auto second=exact_boundary_on_native_plane(extract_exact_section_boundary(lower.section).boundary,0,4);
    const auto m1=exact_section_vertical_surface_moments(first.boundary,a,size,0,10,ExactBoundaryCellOwnership::MaterialSide);
    const auto m2=exact_section_vertical_surface_moments(second.boundary,a,size,0,10,ExactBoundaryCellOwnership::MaterialSide);
    for (std::size_t i=0;i<27;++i)
        CHECK_THAT(m1.values[i]+m2.values[i],Catch::Matchers::WithinRel(owner.values[i],1e-13));
    const auto zero=exact_section_vertical_surface_moments(edge.boundary,a,size,2,2,ExactBoundaryCellOwnership::MaterialSide);
    CHECK_THAT(zero.values[0],Catch::Matchers::WithinAbs(0.,0.));
    CHECK_THROWS_AS(exact_section_vertical_surface_moments(edge.boundary,a,size,3,2,ExactBoundaryCellOwnership::MaterialSide),std::invalid_argument);
    auto invalid=size; invalid[0]=0;
    CHECK_THROWS_AS(exact_section_vertical_surface_moments(edge.boundary,a,invalid,0,10,ExactBoundaryCellOwnership::MaterialSide),std::invalid_argument);
    // A hole boundary points clockwise: the material cell is outside it.
    ExPolygon holed(Polygon{{0,0},{10,0},{10,10},{0,10}});
    holed.holes.emplace_back(Polygon{{2,2},{2,4},{4,4},{4,2}});
    const auto hole_edge=exact_boundary_on_native_plane(extract_exact_section_boundary(make_exact_section({holed}).section).boundary,0,2);
    const auto outside=exact_section_vertical_surface_moments(hole_edge.boundary,{0,0,0},{2,10,10},0,10,ExactBoundaryCellOwnership::MaterialSide);
    const auto inside=exact_section_vertical_surface_moments(hole_edge.boundary,{2,0,0},{2,10,10},0,10,ExactBoundaryCellOwnership::MaterialSide);
    CHECK_THAT(outside.values[0],Catch::Matchers::WithinRel(20e-18,1e-13));
    CHECK_THAT(inside.values[0],Catch::Matchers::WithinAbs(0.,0.));
}

TEST_CASE("Exact boundary normal records preserve rational spans and hole orientation", "[ThermoelasticExactSections]")
{
    const auto triangle=make_exact_section({ExPolygon(Polygon{{0,0},{5,0},{0,1}})});
    const auto strip=make_exact_section({ExPolygon(Polygon{{4,0},{6,0},{6,1},{4,1}})});
    const auto clipped=boolean_exact_sections(triangle.section,strip.section,ExactSectionOperation::Intersection);
    const auto boundary=extract_exact_section_boundary(clipped.section);
    const auto split=split_exact_boundary_segments(boundary.boundary);
    REQUIRE(split.segments.size()==3);
    auto remaining=boundary.boundary;
    const std::array<coord_t,3> start{0,0,0},size{10,10,10};
    const auto total=exact_section_vertical_surface_moments(boundary.boundary,start,size,0,10,ExactBoundaryCellOwnership::AllIncident);
    SectionIntegration::Moments summed{};
    for (const auto& record : split.segments) {
        CHECK(record.span.segment_count()==1);
        CHECK_THAT(record.outward_normal[0]*record.outward_normal[0]+record.outward_normal[1]*record.outward_normal[1],
            Catch::Matchers::WithinAbs(1.,1e-14));
        const auto moments=exact_section_vertical_surface_moments(record.span,start,size,0,10,ExactBoundaryCellOwnership::AllIncident);
        for (std::size_t i=0;i<27;++i) summed[i]+=moments.values[i];
        remaining=boolean_exact_section_boundaries(remaining,record.span,ExactBoundaryOperation::Difference,ExactBoundaryOrientation::Any).boundary;
    }
    CHECK(remaining.empty());
    for (std::size_t i=0;i<27;++i)
        CHECK_THAT(summed[i],Catch::Matchers::WithinRel(total.values[i],1e-13));
    ExPolygon holed(Polygon{{0,0},{10,0},{10,10},{0,10}});
    holed.holes.emplace_back(Polygon{{2,2},{2,4},{4,4},{4,2}});
    const auto all=extract_exact_section_boundary(make_exact_section({holed}).section);
    const auto hole=split_exact_boundary_segments(exact_boundary_on_native_plane(all.boundary,0,2).boundary);
    REQUIRE(hole.segments.size()==1);
    CHECK_THAT(hole.segments.front().outward_normal[0],Catch::Matchers::WithinAbs(1.,0.));
    CHECK_THAT(hole.segments.front().outward_normal[1],Catch::Matchers::WithinAbs(0.,0.));
    const auto outer=split_exact_boundary_segments(exact_boundary_on_native_plane(all.boundary,0,0).boundary);
    REQUIRE(outer.segments.size()==1);
    CHECK_THAT(outer.segments.front().outward_normal[0],Catch::Matchers::WithinAbs(-1.,0.));
    CHECK_THROWS_AS(split_exact_boundary_segments(boundary.boundary,split.work-1),ResourceLimit);
    std::size_t callbacks=0;
    CHECK_THROWS_AS(split_exact_boundary_segments(all.boundary,4000000,
        [&] { if (++callbacks==20) throw std::runtime_error("cancel boundary normals"); }),std::runtime_error);
    CHECK(callbacks==20);
    CHECK(split_exact_boundary_segments(ExactSectionBoundary{}).segments.empty());
}

TEST_CASE("Exact union retains holes and overlapping area without counting shared boundaries twice", "[ThermoelasticExactSections]")
{
    ExPolygon holed(Polygon{{0,0},{10,0},{10,10},{0,10}});
    holed.holes.emplace_back(Polygon{{2,2},{2,4},{4,4},{4,2}});
    const ExPolygons second{ExPolygon(Polygon{{3,2},{12,2},{12,6},{3,6}})};
    const auto a=make_exact_section({holed}),b=make_exact_section(second);
    const auto joined=boolean_exact_sections(a.section,b.section,ExactSectionOperation::Union);
    const auto overlap=boolean_exact_sections(a.section,b.section,ExactSectionOperation::Intersection);
    const std::array<coord_t,3> start{0,0,0},size{20,20,10};
    const auto ma=exact_section_volume_moments(a.section,start,size,0,10),mb=exact_section_volume_moments(b.section,start,size,0,10);
    const auto mu=exact_section_volume_moments(joined.section,start,size,0,10),mi=exact_section_volume_moments(overlap.section,start,size,0,10);
    for (std::size_t i=0;i<27;++i)
        CHECK_THAT(mu.values[i]+mi.values[i],Catch::Matchers::WithinRel(ma.values[i]+mb.values[i],1e-13));
    const auto native=exact_section_volume_moments({holed},second,ExactSectionOperation::Union,start,size,0,10);
    for (std::size_t i=0;i<27;++i) CHECK_THAT(native.values[i],Catch::Matchers::WithinRel(mu.values[i],1e-13));
    CHECK_FALSE(exact_section_empty({holed},{holed},ExactSectionOperation::Union,2048,4000000,{},true).empty);
    CHECK_FALSE(exact_section_empty({holed},{holed},ExactSectionOperation::Union).empty);
    const auto same=boolean_exact_sections(a.section,a.section,ExactSectionOperation::Union);
    CHECK(boolean_exact_sections(same.section,a.section,ExactSectionOperation::SymmetricDifference).section.empty());
    const auto touching=boolean_exact_sections(make_exact_section({box(0,100)}).section,
        make_exact_section({box(100,200)}).section,ExactSectionOperation::Union);
    CHECK(boolean_exact_sections(touching.section,make_exact_section({box(0,200)}).section,ExactSectionOperation::SymmetricDifference).section.empty());
}

TEST_CASE("Exact point membership distinguishes closed boundaries from rational hole interiors", "[ThermoelasticExactSections]")
{
    const auto outer=make_exact_section({ExPolygon(Polygon{{0,0},{100,0},{100,100},{0,100}})});
    const auto cut=make_exact_section({ExPolygon(Polygon{{20,20},{71,20},{20,31}}),
        ExPolygon(Polygon{{60,20},{80,20},{80,30},{60,30}})});
    const auto section=boolean_exact_sections(outer.section,cut.section,ExactSectionOperation::Difference);
    CHECK(validate_exact_section_topology(section.section)>0);
    for (const Point point : {Point(40,25),Point(61,25),Point(-1,50)}) {
        const auto location=exact_closed_contains_point(section.section,point);
        CHECK_FALSE(location.contained);
        CHECK_FALSE(location.boundary);
    }
    const auto interior=exact_closed_contains_point(section.section,Point(40,28));
    CHECK(interior.contained);
    CHECK_FALSE(interior.boundary);
    for (const Point point : {Point(0,50),Point(20,20),Point(71,30)}) {
        const auto location=exact_closed_contains_point(section.section,point);
        CHECK(location.contained);
        CHECK(location.boundary);
    }
    const auto empty=exact_closed_contains_point(ExactSection{},Point(0,0));
    CHECK_FALSE(empty.contained);
    CHECK_FALSE(empty.boundary);
    CHECK_THROWS_AS(exact_closed_contains_point(section.section,Point(40,28),interior.work-1),ResourceLimit);
    CHECK_THROWS_AS(exact_closed_contains_point(section.section,Point(40,28),4000000,
        [] { throw std::runtime_error("cancel exact membership"); }),std::runtime_error);
    CHECK_THROWS_AS(exact_closed_contains_point(section.section,Point(coord_t(10000000001LL),0)),std::invalid_argument);
}

TEST_CASE("Exact topology rejects point contacts within components without joining separate islands", "[ThermoelasticExactSections]")
{
    ExPolygon touching_outer(Polygon{{0,0},{10,0},{10,10},{0,10}});
    touching_outer.holes.emplace_back(Polygon{{0,5},{3,7},{3,3}});
    const auto pinched=make_exact_section({touching_outer});
    CHECK_THROWS_AS(validate_exact_section_topology(pinched.section),std::invalid_argument);
    ExPolygon touching_holes(Polygon{{0,0},{20,0},{20,20},{0,20}});
    touching_holes.holes.emplace_back(Polygon{{2,2},{2,6},{6,6},{6,2}});
    touching_holes.holes.emplace_back(Polygon{{6,6},{6,10},{10,10},{10,6}});
    const auto pinched_holes=make_exact_section({touching_holes});
    CHECK_THROWS_AS(validate_exact_section_topology(pinched_holes.section),std::invalid_argument);
    const auto touching=make_exact_section({ExPolygon(Polygon{{0,0},{10,0},{10,10},{0,10}}),
        ExPolygon(Polygon{{10,10},{20,10},{20,20},{10,20}})});
    CHECK_THROWS_AS(validate_exact_section_topology(touching.section),std::invalid_argument);
    const auto separate=make_exact_section({ExPolygon(Polygon{{0,0},{10,0},{10,10},{0,10}}),
        ExPolygon(Polygon{{11,11},{21,11},{21,21},{11,21}})});
    REQUIRE(exact_section_components(separate.section).sections.size()==2);
    const auto work=validate_exact_section_topology(separate.section);
    CHECK(work>0);
    const auto point=exact_closed_contains_point(separate.section,Point(10,10));
    CHECK(point.contained);
    CHECK(point.boundary);
    CHECK_THROWS_AS(validate_exact_section_topology(separate.section,work-1),ResourceLimit);
    CHECK_THROWS_AS(validate_exact_section_topology(separate.section,4000000,
        [] { throw std::runtime_error("cancel exact topology"); }),std::runtime_error);
}

TEST_CASE("Occupied assignment certificates preserve pieces and bound representation changes", "[ThermoelasticExactSections]")
{
    const ExPolygons source{box(0,500),box(500,1000)};
    const auto original=source;
    const BoundingBox tile(Point(-100,-100),Point(2000,2000));
    const auto certificate=certify_occupied_assignment(source,{box(1,1001)},tile);
    const std::vector<std::size_t> expected{0,0};
    CHECK(certificate.occupied_indices==expected);
    CHECK(certificate.representation_allowance_nm==1);
    CHECK(source==original);
    CHECK_THROWS_AS(certify_occupied_assignment(source,{box(2,1002)},tile),std::invalid_argument);
    CHECK_THROWS_AS(certify_occupied_assignment(source,{box(1,1001)},BoundingBox(Point(0,0),Point(2000,2000))),std::invalid_argument);
    CHECK_THROWS_AS(certify_occupied_assignment({box(0,600),box(500,1000)},{box(0,1000)},tile),std::invalid_argument);
    try {
        certify_occupied_assignment({box(0,600),box(500,1000)},{box(0,1000)},tile);
        FAIL("Overlapping original pieces must retain typed diagnostic ownership.");
    } catch (const OccupiedAssignmentOverlap& overlap) {
        CHECK(overlap.original);
        CHECK(overlap.first==0);
        CHECK(overlap.second==1);
        CHECK(overlap.work>0);
    }
    CHECK_THROWS_AS(certify_occupied_assignment({box(0,100),box(200,300)},{box(0,300)},tile),std::invalid_argument);
    auto perforated=box(0,1000);
    perforated.holes.emplace_back(Polygon{{1,1},{1,2},{2,1}});
    CHECK_THROWS_AS(certify_occupied_assignment(source,{perforated},tile),std::invalid_argument);
    CHECK_THROWS_AS(certify_occupied_assignment(source,{box(0,1000)},tile,1),ResourceLimit);
    CHECK_THROWS_AS(certify_occupied_assignment(source,{box(0,1000)},tile,4000000,
        [] { throw std::runtime_error("cancel"); }),std::runtime_error);
}

TEST_CASE("Exact integer unions preserve native regions without rounding", "[ThermoelasticExactSections]")
{
    const ExPolygons source{box(0,1000),box(500,1500)};
    const auto result=exact_integer_section_union(source);
    REQUIRE(result.polygons.size()==1);
    CHECK(exact_section_empty(result.polygons,source,ExactSectionOperation::SymmetricDifference).empty);
    CHECK_THROWS_AS(exact_integer_section_union(source,1),ResourceLimit);
    CHECK_THROWS_AS(exact_integer_section_union(source,4000000,[] { throw std::runtime_error("cancel"); }),std::runtime_error);
    const ExPolygons fractional{ExPolygon(Polygon{{0,0},{5,0},{0,3}}),
        ExPolygon(Polygon{{2,-1},{4,-1},{4,4},{2,4}})};
    CHECK_THROWS_AS(exact_integer_section_union(fractional),std::invalid_argument);
}

TEST_CASE("Exact section predicates distinguish face contact from positive area", "[ThermoelasticExactSections]")
{
    const ExPolygons a{box(0, 1000)}, touching{box(1000, 2000)}, overlap{box(999, 2000)};
    CHECK(exact_section_empty(a, touching, ExactSectionOperation::Intersection).empty);
    CHECK_FALSE(exact_section_empty(a, overlap, ExactSectionOperation::Intersection).empty);
    CHECK(exact_section_empty(a, a, ExactSectionOperation::SymmetricDifference).empty);
    CHECK_FALSE(exact_section_empty(a, touching, ExactSectionOperation::Difference).empty);
    CHECK(exact_section_empty(a, {box(-1, 1001)}, ExactSectionOperation::Difference).empty);
}

TEST_CASE("Native exact containment preserves boundary touching holes", "[ThermoelasticExactSections]")
{
    auto touching = box(0,1000);
    touching.holes.emplace_back(Points{{0,200},{0,800},{200,800},{200,200}});
    const ExPolygons notch{ExPolygon(Polygon{{0,0},{1000,0},{1000,1000},{0,1000},
        {0,800},{200,800},{200,200},{0,200}})};
    CHECK_THROWS_AS(exact_section_empty({touching},notch,ExactSectionOperation::Difference),std::invalid_argument);
    // Equality of the original ordered chains needs no arrangement repair.
    const auto identity = exact_section_empty({touching},{touching},ExactSectionOperation::Difference,
        1024,17,{},true);
    CHECK(identity.empty);
    CHECK(identity.work == 17);
    CHECK_THROWS_AS(exact_section_empty({touching},{touching},ExactSectionOperation::Difference,
        1024,16,{},true),ResourceLimit);
    CHECK_THROWS_AS(exact_section_empty({touching},{touching},ExactSectionOperation::Difference),std::invalid_argument);
    const auto equal = exact_section_empty({touching},notch,ExactSectionOperation::SymmetricDifference,
        1024,2000000,{},true);
    CHECK(equal.empty);
    CHECK(equal.work <= 2000000);
    CHECK(exact_section_empty({touching},{box(0,1000)},ExactSectionOperation::Difference,
        1024,2000000,{},true).empty);
    CHECK_FALSE(exact_section_empty({box(0,1000)},{touching},ExactSectionOperation::Difference,
        1024,2000000,{},true).empty);
    CHECK_THROWS_AS(exact_section_empty({touching},notch,ExactSectionOperation::Difference,
        1024,1,{},true),ResourceLimit);
    CHECK_THROWS_AS(exact_section_empty({touching},notch,ExactSectionOperation::Difference,
        1024,2000000,[] { throw std::invalid_argument("cancel"); },true),std::invalid_argument);
}

TEST_CASE("Joint activation clipping shares rounded interface junctions", "[ThermoelasticExactSections]")
{
    // The partial upper material meets the diagonal at (501,166.833...).
    // Integer clipping moves that point in only its polygon, into the lower one.
    const std::array<ExPolygons,2> reference{{
        {ExPolygon(Polygon{{0,0},{1000,333},{1000,1000},{0,1000}})},
        {ExPolygon(Polygon{{0,0},{1000,0},{1000,333}})}
    }};
    const std::array<ExPolygons,2> deposited{{{box(0,501)},{box(0,1000)}}};
    const auto independent_first=intersection_ex(reference[0],deposited[0]);
    const auto independent_second=intersection_ex(reference[1],deposited[1]);
    REQUIRE_FALSE(exact_section_empty(independent_first,independent_second,
        ExactSectionOperation::Intersection).empty);
    const auto clipped=exact_clip_reference_sections(reference,deposited);
    REQUIRE(clipped.sections[0].size()==1);
    REQUIRE(clipped.sections[1].size()==1);
    CHECK(exact_section_empty(clipped.sections[0],clipped.sections[1],
        ExactSectionOperation::Intersection).empty);
    CHECK(clipped.representation_allowance_nm==1);
    // A remote unprinted strip is not activated by the shared rounding.
    CHECK(exact_section_empty(clipped.sections[0],{box(502,1000)},
        ExactSectionOperation::Intersection).empty);
    CHECK_FALSE(exact_section_empty(clipped.sections[1],{},
        ExactSectionOperation::SymmetricDifference).empty);
    const std::array<ExPolygons,2> reversed_reference{{reference[1],reference[0]}};
    const std::array<ExPolygons,2> reversed_deposited{{deposited[1],deposited[0]}};
    const auto reversed=exact_clip_reference_sections(reversed_reference,reversed_deposited);
    for (std::size_t i=0;i<2;++i)
        CHECK(exact_section_empty(clipped.sections[i],reversed.sections[1-i],
            ExactSectionOperation::SymmetricDifference).empty);
    CHECK_THROWS_AS(exact_clip_reference_sections(reference,deposited,1),ResourceLimit);
    CHECK_THROWS_AS(exact_clip_reference_sections(reference,deposited,4000000,
        [] { throw std::runtime_error("cancel"); }),std::runtime_error);
    std::size_t charged=0;
    const auto metered=exact_clip_reference_sections(reference,deposited,4000000,{},
        [&](std::size_t work) { charged+=work; });
    CHECK(charged==metered.work);
    CHECK_THROWS_AS(exact_clip_reference_sections(reference,deposited,4000000,{},
        [](std::size_t) { throw ResourceLimit("outer budget"); }),ResourceLimit);
    auto invalid=reference;
    invalid[1]={box(0,1000)};
    CHECK_THROWS_AS(exact_clip_reference_sections(invalid,deposited),std::invalid_argument);
}

TEST_CASE("Joint activation clipping preserves holes and inactive materials", "[ThermoelasticExactSections]")
{
    auto holed=box(0,1000);
    holed.holes.emplace_back(Points{{200,200},{200,800},{800,800},{800,200}});
    const std::array<ExPolygons,2> reference{{{holed},{box(1000,2000)}}};
    const std::array<ExPolygons,2> deposited{{{box(0,1000)}, {}}};
    const auto result=exact_clip_reference_sections(reference,deposited);
    REQUIRE(result.sections[0].size()==1);
    CHECK(result.sections[0].front().holes.size()==1);
    CHECK(result.sections[1].empty());
    CHECK(exact_section_empty(result.sections[0],reference[0],
        ExactSectionOperation::SymmetricDifference).empty);
}

TEST_CASE("A shared bounded rounding phase separates nearly touching activated interfaces", "[ThermoelasticExactSections]")
{
    // Exact parallel interfaces have a 1/5 nm vertical gap. Nearest rounding
    // moves the first cut at y=20.6 upward and the second at y=20.4 downward,
    // crossing the independently represented boundary chains.
    const std::array<ExPolygons,2> reference{
        ExPolygons{ExPolygon(Polygon{{0,0},{500,0},{500,100}})},
        ExPolygons{ExPolygon(Polygon{{4,1},{499,100},{499,200},{4,200}})}};
    const std::array<ExPolygons,2> deposited{{{box(0,103)},{box(101,500)}}};
    // Explicit nearest rounding of the exact rational intersections, not
    // Clipper's independent intersection representation.
    const ExPolygons nearest_a{ExPolygon(Polygon{{0,0},{103,0},{103,21}})};
    const ExPolygons nearest_b{ExPolygon(Polygon{{101,20},{499,100},{499,200},{101,200}})};
    REQUIRE_FALSE(exact_section_empty(nearest_a,nearest_b,ExactSectionOperation::Intersection).empty);
    const auto result = exact_clip_reference_sections(reference,deposited);
    REQUIRE(result.sections[0].size() == 1);
    REQUIRE(result.sections[1].size() == 1);
    CHECK(exact_section_empty(result.sections[0],result.sections[1],ExactSectionOperation::Intersection).empty);
    CHECK(result.representation_allowance_nm == 1);
    CHECK_THAT(result.vertex_rounding_bound_m,Catch::Matchers::WithinAbs(1e-9,1e-20));
    const std::array<ExPolygons,2> swapped_reference{{reference[1],reference[0]}};
    const std::array<ExPolygons,2> swapped_deposited{{deposited[1],deposited[0]}};
    const auto swapped = exact_clip_reference_sections(swapped_reference,swapped_deposited);
    CHECK(exact_section_empty(result.sections[0],swapped.sections[1],ExactSectionOperation::SymmetricDifference).empty);
    CHECK(exact_section_empty(result.sections[1],swapped.sections[0],ExactSectionOperation::SymmetricDifference).empty);
    CHECK_THROWS_AS(exact_clip_reference_sections(reference,deposited,result.work-1),ResourceLimit);
    std::size_t charged = 0;
    try {
        exact_clip_reference_sections(reference,deposited,4000000,{},[&](std::size_t work) {
            charged += work;
            if (charged > result.work/2) throw std::invalid_argument("external rounding cancellation");
        });
        FAIL("External work cancellation was ignored.");
    } catch (const std::invalid_argument& error) {
        CHECK(std::string(error.what()) == "external rounding cancellation");
    }
}

TEST_CASE("Exact single material activation distinguishes a narrow gap from a genuine tiny island", "[ThermoelasticExactSections]")
{
    const std::array<ExPolygons,2> reference{{{ExPolygon(Polygon{{0,0},{1000,0},{1000,333}})}, {}}};
    const std::array<ExPolygons,2> above{{{ExPolygon(Polygon{{0,1},{1000,334},{0,1000}})}, {}}};
    const auto disjoint = exact_clip_reference_sections(reference, above);
    CHECK(disjoint.sections[0].empty());
    CHECK(disjoint.sections[1].empty());
    const ExPolygon tiny(Polygon{{100,100},{104,101},{102,101}});
    const std::array<ExPolygons,2> full{{{box(0,1000)}, {}}};
    const std::array<ExPolygons,2> deposited{{{tiny}, {}}};
    const auto retained = exact_clip_reference_sections(full, deposited);
    REQUIRE(retained.sections[0].size() == 1);
    CHECK(exact_section_empty(retained.sections[0], {tiny}, ExactSectionOperation::SymmetricDifference).empty);
    // Exact positive geometry is retained, even when it cannot be certified
    // under the requested representation band. No area cutoff deletes it.
    const LocalSectionTile tile{{0,0},{0,0},{1000,1000},{{0.,.001,0,retained.sections[0],7}}};
    CHECK_THROWS_AS(qualify_local_section_topology({tile},0.),std::invalid_argument);
}

TEST_CASE("Quarter axis rounding preserves a bent activated shared interface", "[ThermoelasticExactSections]")
{
    // Translated local coordinates from the failing interface. Intersecting
    // A-B-Q with A-B-D introduces P=(107086.2672...,259822.2220...) on A-Q.
    // Nearest and all eighth-nm phases map P to (107086,259822), placing Q
    // across the rounded P-B edge. A quarter X phase preserves its side.
    const Point a(0,0), q(107087,259824), b(214844,506414), c(129002,312996), d(107018,259666), left(-100000,100000);
    const ExPolygons first{ExPolygon(Polygon{a,b,q})};
    const ExPolygons second{ExPolygon(Polygon{a,q,c,left})};
    const std::array<ExPolygons,2> reference{first,second};
    const std::array<ExPolygons,2> deposited{ExPolygons{ExPolygon(Polygon{a,b,d})},second};
    REQUIRE(exact_section_empty(first,second,ExactSectionOperation::Intersection).empty);
    const auto result = exact_clip_reference_sections(reference,deposited);
    REQUIRE(result.sections[0].size() == 1);
    REQUIRE(result.sections[1].size() == 1);
    CHECK(exact_section_empty(result.sections[0],result.sections[1],ExactSectionOperation::Intersection).empty);
    CHECK(result.representation_allowance_nm == 1);
    CHECK_THAT(result.vertex_rounding_bound_m,Catch::Matchers::WithinAbs(1e-9,1e-20));
    const Point admitted(107087,259822);
    CHECK(std::find(result.sections[0].front().contour.points.begin(),result.sections[0].front().contour.points.end(),admitted) !=
        result.sections[0].front().contour.points.end());
    const auto swapped = exact_clip_reference_sections({second,first},{second,deposited[0]});
    CHECK(exact_section_empty(result.sections[0],swapped.sections[1],ExactSectionOperation::SymmetricDifference).empty);
    CHECK(exact_section_empty(result.sections[1],swapped.sections[0],ExactSectionOperation::SymmetricDifference).empty);
    CHECK_THROWS_AS(exact_clip_reference_sections(reference,deposited,result.work-1),ResourceLimit);
}
TEST_CASE("Exact segment coverage preserves boundary contact and rejects narrow gaps", "[ThermoelasticExactSections]")
{
    CHECK(exact_section_covers_segment({box(0,1000)},{0,0},{1000,0}).covered);
    CHECK(exact_section_covers_segment({box(0,1000),box(1000,2000)},{0,500},{2000,500}).covered);
    CHECK_FALSE(exact_section_covers_segment({box(0,1000),box(1001,2000)},{0,500},{2000,500}).covered);
    auto hole = box(0,1000);
    hole.holes.emplace_back(Points{{499,200},{499,800},{500,800},{500,200}});
    CHECK_FALSE(exact_section_covers_segment({hole},{0,500},{1000,500}).covered);
    CHECK(exact_section_covers_segment({hole},{0,200},{1000,200}).covered);
    // A remote invalid ring cannot change coverage on the bottom edge. This
    // local predicate must not claim that the whole polygon is valid.
    hole.holes.emplace_back(Points{{600,300},{800,700},{600,700},{800,300}});
    CHECK(exact_section_covers_segment({hole},{0,0},{1000,0}).covered);
    CHECK_THROWS_AS(exact_section_covers_segment({hole},{0,500},{1000,500}),std::invalid_argument);
    CHECK_THROWS_AS(exact_section_covers_segment({hole},{0,500},{1000,500},1),ResourceLimit);
    CHECK_THROWS_AS(exact_section_covers_segment({hole},{0,500},{1000,500},4000000,
        [] { throw std::runtime_error("cancel"); }),std::runtime_error);
}
TEST_CASE("Boundary band sensitivity never changes strict segment coverage", "[ThermoelasticExactSections]")
{
    const ExPolygons polygons{box(0,1000)};
    for (const bool reverse : {false,true}) {
        const Point a = reverse ? Point(900,-1) : Point(100,-1);
        const Point b = reverse ? Point(100,-1) : Point(900,-1);
        const auto strict = exact_section_covers_segment(polygons,a,b);
        CHECK_FALSE(strict.covered);
        CHECK_FALSE(strict.within_boundary_band);
        const auto band = exact_section_covers_segment(polygons,a,b,4000000,{},1);
        CHECK_FALSE(band.covered);
        CHECK(band.within_boundary_band);
    }
    const auto farther = exact_section_covers_segment(polygons,{100,-2},{900,-2},4000000,{},1);
    CHECK_FALSE(farther.covered);
    CHECK_FALSE(farther.within_boundary_band);
    const auto point = exact_section_covers_segment(polygons,{0,-1},{0,-1},4000000,{},1);
    CHECK_FALSE(point.covered);
    CHECK(point.within_boundary_band);
    const auto gap = exact_section_covers_segment({box(0,499),box(500,1000)},{0,500},{1000,500},4000000,{},1);
    CHECK_FALSE(gap.covered);
    CHECK(gap.within_boundary_band);
    const auto wide_gap = exact_section_covers_segment({box(0,490),box(510,1000)},{0,500},{1000,500},4000000,{},1);
    CHECK_FALSE(wide_gap.within_boundary_band);
    CHECK_THROWS_AS(exact_section_covers_segment(polygons,{100,-1},{900,-1},4000000,{},8),std::invalid_argument);
    CHECK_THROWS_AS(exact_section_covers_segment(polygons,{100,-1},{900,-1},1,{},1),ResourceLimit);
    CHECK_THROWS_AS(exact_section_covers_segment(polygons,{100,-1},{900,-1},4000000,
        [] { throw std::runtime_error("cancel"); },1),std::runtime_error);
}

TEST_CASE("Invalid source boundaries identify their crossing without being repaired", "[ThermoelasticExactSections]")
{
    const ExPolygons crossed{ExPolygon(Polygon({{0,0},{1000,1000},{0,1000},{1000,0}}))};
    try {
        (void)exact_section_covers_segment(crossed,{0,500},{1000,500});
        FAIL("A self-crossing source boundary must not certify ownership.");
    } catch (const std::invalid_argument& error) {
        const std::string message(error.what());
        CHECK(message.find("crossing_edges=0,2") != std::string::npos);
        CHECK(message.find("(0,0) (1000,1000) (0,1000) (1000,0)") != std::string::npos);
    }
}

TEST_CASE("Native nonzero filling preserves crossing lobes and cancels opposite boundary chains", "[ThermoelasticExactSections]")
{
    const ExPolygons crossed{ExPolygon(Polygon({{0,0},{1000,1000},{0,1000},{1000,0}}))};
    CHECK(exact_native_section_covers_segment(crossed,{400,100},{600,100}).covered);
    CHECK(exact_native_section_covers_segment(crossed,{400,900},{600,900}).covered);
    CHECK_FALSE(exact_native_section_covers_segment(crossed,{100,400},{100,600}).covered);
    CHECK(exact_native_section_covers_segment(crossed,{500,500},{500,500}).covered);
    CHECK_THROWS_AS(exact_section_covers_segment(crossed,{400,100},{600,100}),std::invalid_argument);
    const ExPolygons twice{ExPolygon(Polygon({{0,0},{1000,0},{1000,1000},{0,1000},
        {0,0},{1000,0},{1000,1000},{0,1000}}))};
    CHECK(exact_native_section_covers_segment(twice,{100,500},{900,500}).covered);
    const ExPolygons canceled{ExPolygon(Polygon({{0,0},{1000,0},{1000,1000},{0,1000},
        {0,0},{0,1000},{1000,1000},{1000,0}}))};
    CHECK_FALSE(exact_native_section_covers_segment(canceled,{100,500},{900,500}).covered);
    CHECK_FALSE(exact_native_section_covers_segment(canceled,{0,0},{1000,0}).covered);
    const ExPolygons collinear{ExPolygon(Polygon{{0,0},{1000,0},{500,0}})};
    CHECK_FALSE(exact_native_section_covers_segment(collinear,{0,0},{1000,0}).covered);
    const ExPolygons collapsed{ExPolygon(Polygon{{0,0},{0,0},{0,0}})};
    CHECK_FALSE(exact_native_section_covers_segment(collapsed,{0,0},{0,0}).covered);
    auto reversed = box(0,1000);
    reversed.contour.reverse();
    CHECK_FALSE(exact_native_section_covers_segment({box(0,1000),reversed},{100,500},{900,500}).covered);
    const ExPolygons native_crossing{ExPolygon(Polygon{{52959738,122932547},{49592121,114802399},
        {51249672,116713580},{50510568,117019727}})};
    // The last endpoint is only ~0.162 nm across the first edge. Its lobe
    // survives exact winding reconstruction without snapping to that edge.
    CHECK(exact_native_section_covers_segment(native_crossing,{50510568,117019727},{50510568,117019727}).covered);
    auto misplaced_hole = box(2000,3000);
    misplaced_hole.holes.push_back(reversed.contour);
    // The outer chain is remote, but the clockwise chain still has nonzero
    // winding here. Local filtering must not drop it or reverse it.
    CHECK(exact_native_section_covers_segment({misplaced_hole},{100,500},{900,500}).covered);
    CHECK_THROWS_AS(exact_native_section_covers_segment(crossed,{400,100},{600,100},100),ResourceLimit);
    CHECK_THROWS_AS(exact_native_section_covers_segment(crossed,{400,100},{600,100},4000000,
        [] { throw std::invalid_argument("cancel"); }),std::invalid_argument);
}

TEST_CASE("Disjoint native winding components do not share an exceptional triangulation", "[ThermoelasticExactSections]")
{
    ExPolygons polygons;
    for (coord_t i = 0; i < 65; ++i) polygons.push_back(box(i*2000,i*2000+1000));
    polygons.emplace_back(Polygon{{130000,0},{131000,1000},{130000,1000},{131000,0}});
    const auto coverage = exact_native_section_covers_segment(polygons,{500,100},{130500,100});
    CHECK_FALSE(coverage.covered); // Genuine spaces between the independent components.
    CHECK(coverage.work <= 4000000);
    CHECK_THROWS_AS(exact_native_section_covers_segment(polygons,{500,100},{130500,100},100),ResourceLimit);
}

TEST_CASE("Native winding rejects excessive crossing construction before triangulation", "[ThermoelasticExactSections]")
{
    Points chain;
    for (coord_t i = 0; i < 100; ++i) {
        chain.emplace_back(i,0);
        chain.emplace_back(200-i,1000);
    }
    const ExPolygons tangled{ExPolygon(Polygon(std::move(chain)))};
    CHECK_THROWS_AS(exact_native_section_covers_segment(tangled,{0,500},{200,500}),ResourceLimit);
}

TEST_CASE("Exact native segment coverage crops remote chains without inventing boundary support", "[ThermoelasticExactSections]")
{
    Polygon outer;
    outer.points = {{0,0},{10000,0},{10000,1000}};
    for (coord_t i=1; i<=2600; ++i) outer.points.emplace_back(10000-i*10000/2600,1000);
    ExPolygon native(outer);
    native.holes.emplace_back(Points{{400,50},{400,150},{600,150},{600,50}});
    const ExPolygons polygons{native};
    const auto interior = exact_native_section_covers_segment(polygons,{100,10},{900,10});
    CHECK(interior.covered);
    CHECK(interior.work <= 4000000);
    CHECK_FALSE(exact_native_section_covers_segment(polygons,{100,100},{900,100}).covered);
    const auto near = exact_native_section_covers_segment(polygons,{100,-1},{900,-1},4000000,{},1);
    CHECK_FALSE(near.covered);
    CHECK(near.within_boundary_band);
    const auto hole = exact_native_section_covers_segment(polygons,{450,51},{550,51},4000000,{},1);
    CHECK_FALSE(hole.covered);
    CHECK(hole.within_boundary_band);
    const auto exterior = exact_native_section_covers_segment(polygons,{100,-5},{900,-5},4000000,{},1);
    CHECK_FALSE(exterior.covered);
    CHECK_FALSE(exterior.within_boundary_band);
    CHECK_THROWS_AS(exact_native_section_covers_segment(polygons,{100,10},{900,10},interior.work-1),ResourceLimit);
    CHECK_THROWS_AS(exact_native_section_covers_segment(polygons,{100,10},{900,10},4000000,
        [] { throw std::invalid_argument("cancel native segment crop"); }),std::invalid_argument);
    native.contour.points.back().x() = coord_t(10000000001LL);
    CHECK_THROWS_AS(exact_native_section_covers_segment({native},{100,10},{900,10}),std::invalid_argument);
}

TEST_CASE("Exact topology diagnostics count union components and holes without rounding", "[ThermoelasticExactSections]")
{
    auto holed = box(0, 1000);
    holed.holes.emplace_back(Points{{200, 200}, {200, 800}, {800, 800}, {800, 200}});
    const auto separate = diagnose_exact_section_topology({holed, box(2000, 3000)});
    CHECK(separate.components == 2);
    CHECK(separate.holes == 1);
    const auto filled = diagnose_exact_section_topology({holed, box(0, 1000)});
    CHECK(filled.components == 1);
    CHECK(filled.holes == 0);
    CHECK_THROWS_AS(diagnose_exact_section_topology({holed}, 1), ResourceLimit);
    CHECK_THROWS_AS(diagnose_exact_section_topology({holed}, 4000000,
        [] { throw std::runtime_error("cancel"); }), std::runtime_error);
}
TEST_CASE("Exact topology certificates match each component and its holes", "[ThermoelasticExactSections]")
{
    auto holed = box(0, 1000);
    holed.holes.emplace_back(Points{{200, 200}, {200, 800}, {800, 800}, {800, 200}});
    const ExPolygons before{holed, box(2000, 3000)};
    const ExPolygons reordered{box(2000, 3001), holed};
    const auto work = qualify_exact_section_topology(before, reordered, 4000000);
    CHECK(work > 0);
    CHECK_THROWS_AS(qualify_exact_section_topology(before, reordered, work - 1), ResourceLimit);
    CHECK_THROWS_AS(qualify_exact_section_topology(before, {box(0, 1000), box(2000, 3000)}, 4000000), std::invalid_argument);
    auto moved_hole = box(2000, 3000);
    moved_hole.holes.emplace_back(Points{{2200, 200}, {2200, 800}, {2800, 800}, {2800, 200}});
    // Same total components and holes, but the hole changed owner.
    CHECK_THROWS_AS(qualify_exact_section_topology(before, {box(0, 1000), moved_hole}, 4000000), std::invalid_argument);
    CHECK_THROWS_AS(qualify_exact_section_topology(before, reordered, 4000000,
        [] { throw std::runtime_error("cancel"); }), std::runtime_error);
    // A rounded component may touch a different source component at an edge;
    // that is not a second positive-area correspondence.
    CHECK(qualify_exact_section_topology({box(0, 1000), box(1001, 2000)},
        {box(0, 1001), box(1002, 2000)}, 4000000) > 0);
}

TEST_CASE("Exact topology certificates preserve holes hidden by integer Boolean rounding", "[ThermoelasticExactSections]")
{
    // The slanted lip closes the outlet less than one grid unit from the
    // outer edge. Integer union output opens an artificial portal, whereas
    // both exact unions contain the same enclosed hole.
    const ExPolygon open_box(Polygon{{0, 0}, {1000, 0}, {1000, 500}, {700, 500}, {700, 200}, {200, 200},
        {200, 800}, {700, 800}, {700, 501}, {1000, 501}, {1000, 1000}, {0, 1000}});
    const ExPolygon lip(Polygon{{800, 699}, {1001, 499}, {1001, 1000}, {800, 1000}});
    const ExPolygon moved_lip(Polygon{{800, 698}, {1001, 499}, {1001, 1000}, {800, 1000}});
    const ExPolygons before{open_box, lip}, after{open_box, moved_lip};
    const auto exact_before = diagnose_exact_section_topology(before);
    const auto exact_after = diagnose_exact_section_topology(after);
    CHECK(exact_before.components == 1);
    CHECK(exact_before.holes == 1);
    CHECK(exact_after.components == 1);
    CHECK(exact_after.holes == 1);
    const auto rounded = union_ex(before);
    REQUIRE(rounded.size() == 1);
    CHECK(rounded.front().holes.empty());
    CHECK(qualify_exact_section_topology(before, after, 4000000) > 0);
}
TEST_CASE("Exact section predicates preserve holes and original coordinates", "[ThermoelasticExactSections]")
{
    auto holed = box(0, 1000);
    holed.holes.emplace_back(Points{Point(200, 200), Point(200, 800), Point(800, 800), Point(800, 200)});
    const ExPolygons interior{ExPolygon(Polygon({Point(300, 300), Point(700, 300), Point(700, 700), Point(300, 700)}))};
    CHECK(exact_section_empty({holed}, interior, ExactSectionOperation::Intersection).empty);
    CHECK_FALSE(exact_section_empty({box(0, 1000)}, {holed}, ExactSectionOperation::Difference).empty);
    const ExPolygons lower{ExPolygon(Polygon({Point(0, 0), Point(100000, 99999), Point(100000, -1000)}))};
    auto positive = lower;
    positive.front().contour.reverse();
    const ExPolygons upper{ExPolygon(Polygon({Point(0, 1), Point(100000, 100000), Point(0, 200000)}))};
    CHECK(exact_section_empty(positive, upper, ExactSectionOperation::Intersection).empty);
}
TEST_CASE("Exact predicates refuse captured self-crossing diagnostic rings and bounded exhaustion", "[ThermoelasticExactSections]")
{
    const ExPolygons diagnostic{ExPolygon(Polygon({{54171490,108585150},{54163647,108588647},{54163641,108588649},
        {54163646,108588647},{54171916,108584959},{54185678,108578823},{54182230,108580361},
        {54191215,108576355},{54191216,108576354}}))};
    CHECK_THROWS_AS(exact_section_empty(diagnostic, {}, ExactSectionOperation::Intersection), std::invalid_argument);
    const ExPolygons a{box(0, 1000)};
    CHECK_THROWS_AS(exact_section_empty(a, a, ExactSectionOperation::Intersection, 3), ResourceLimit);
    CHECK_THROWS_AS(exact_section_empty(a, a, ExactSectionOperation::Intersection, 1024, 1), ResourceLimit);
    CHECK_THROWS_AS(exact_section_empty(a, a, ExactSectionOperation::Intersection, 1024, 2000000,
        [] { throw std::runtime_error("cancel"); }), std::runtime_error);
}

TEST_CASE("Dense exact predicates preserve the independent work ceiling", "[ThermoelasticExactSections]")
{
    Polygon boundary;
    for (coord_t x=0;x<275;++x) boundary.points.emplace_back(x,0);
    for (coord_t x=274;x>=0;--x) boundary.points.emplace_back(x,1000);
    const ExPolygons region{ExPolygon(boundary)};
    CHECK_THROWS_AS(exact_section_empty(region,region,ExactSectionOperation::Intersection),ResourceLimit);
    const auto result=exact_section_empty(region,region,ExactSectionOperation::Intersection,2048,4000000);
    CHECK_FALSE(result.empty);
    CHECK(result.work<=4000000);
    CHECK_THROWS_AS(exact_section_empty(region,region,ExactSectionOperation::Intersection,2048,result.work-1),ResourceLimit);
    CHECK_THROWS_AS(exact_section_empty(region,region,ExactSectionOperation::Intersection,2049,4000000),std::invalid_argument);
}

TEST_CASE("Exact predicates retain narrow positive overlap between captured ownership claims", "[ThermoelasticExactSections]")
{
    const ExPolygons assigned{ExPolygon(Polygon{{49673308, 113902865}, {49475995, 113984595},
        {49478209, 113983382}, {49479412, 113983000}, {49480402, 113982459},
        {49830400, 113837455}, {49832902, 113836693}})};
    const ExPolygons claim{ExPolygon(Polygon{{49475375, 113984935}, {49476005, 113984589},
        {49674562, 113902345}})};
    // About two native grid units wide is still positive occupied area, not
    // empty merely because the interface is very narrow.
    CHECK_FALSE(exact_section_empty(assigned, claim, ExactSectionOperation::Intersection).empty);
    const auto witness = diagnose_ownership_witness(assigned, claim, assigned, claim);
    REQUIRE(witness.in_first_native);
    REQUIRE(witness.in_second_native);
    CHECK(*witness.in_first_native);
    CHECK(*witness.in_second_native);
    const auto outside = diagnose_ownership_witness(assigned, claim, {box(0, 1000)}, claim);
    REQUIRE(outside.in_first_native);
    CHECK_FALSE(*outside.in_first_native);
    CHECK_THROWS_AS(diagnose_ownership_witness(assigned, claim, assigned, claim, 1), ResourceLimit);
    CHECK_THROWS_AS(diagnose_ownership_witness(assigned, claim, assigned, claim, 20000000,
        [] { throw std::runtime_error("cancel"); }), std::runtime_error);
    auto malformed = assigned;
    malformed.front().contour.reverse();
    CHECK_THROWS_AS(diagnose_ownership_witness(assigned, claim, malformed, claim), std::invalid_argument);
    const auto absent = diagnose_ownership_witness({box(0, 1000)}, {box(2000, 3000)}, assigned, claim);
    CHECK_FALSE(absent.in_first_native);
    CHECK_FALSE(absent.in_second_native);
}

TEST_CASE("Witness diagnostic vertex admission cannot bypass its total work ceiling", "[ThermoelasticExactSections]")
{
    // 4500 native vertices fit the diagnostic vertex ceiling but their
    // conservative validation reservation alone exceeds twenty million.
    ExPolygons native(1125, box(0, 1000));
    const ExPolygons claim{box(0, 1000)};
    CHECK_THROWS_AS(diagnose_ownership_witness(claim, claim, native, {}, 20000000), ResourceLimit);
    CHECK_THROWS_AS(diagnose_ownership_witness(claim, claim, {}, {}, 20000001), std::invalid_argument);
}

TEST_CASE("Diagnostic nonzero fill respects touching hole cycles without repairing them", "[ThermoelasticExactSections]")
{
    auto native = box(0, 1000);
    native.holes.emplace_back(Points{{500,500},{500,200},{200,200},{200,500},
        {500,500},{500,800},{800,800},{800,500}});
    const ExPolygons hole_claim{ExPolygon(Polygon{{300,300},{400,300},{400,400},{300,400}})};
    CHECK_THROWS_AS(diagnose_ownership_witness(hole_claim, hole_claim, {native}, hole_claim), std::invalid_argument);
    const auto hole = diagnose_ownership_witness(hole_claim, hole_claim, {native}, hole_claim,
        20000000, {}, NativeWitnessFill::NonzeroWinding);
    REQUIRE(hole.in_first_native);
    CHECK_FALSE(*hole.in_first_native);
    REQUIRE(hole.in_second_native);
    CHECK(*hole.in_second_native);
    const ExPolygons solid_claim{ExPolygon(Polygon{{600,300},{700,300},{700,400},{600,400}})};
    const auto solid = diagnose_ownership_witness(solid_claim, solid_claim, {native}, solid_claim,
        20000000, {}, NativeWitnessFill::NonzeroWinding);
    REQUIRE(solid.in_first_native);
    CHECK(*solid.in_first_native);
}

TEST_CASE("Whole intersection ownership requires complete unique native coverage", "[ThermoelasticExactSections]")
{
    const ExPolygons a{box(0, 1000)}, b{box(999, 2000)}, owner{box(0, 2000)}, other{box(3000, 4000)};
    const auto unique = exact_intersection_ownership(a, b, owner, other);
    REQUIRE(unique.nonempty);
    CHECK(unique.first_unique);
    CHECK_FALSE(unique.second_unique);
    CHECK(unique.native_partition);
    const auto competing = exact_intersection_ownership(a, b, owner, owner);
    CHECK_FALSE(competing.first_unique);
    CHECK_FALSE(competing.second_unique);
    CHECK_FALSE(competing.native_partition);
    CHECK(competing.native_covers);
    CHECK_FALSE(competing.native_disjoint);
    const auto missing = exact_intersection_ownership(a, b, {box(0, 999)}, other);
    CHECK_FALSE(missing.first_unique);
    CHECK_FALSE(missing.native_partition);
    CHECK_FALSE(missing.native_covers);
    CHECK(missing.native_disjoint);
    const auto split = exact_intersection_ownership(owner, owner, {box(0, 1000)}, {box(1000, 2000)});
    CHECK(split.native_partition);
    CHECK_FALSE(split.first_unique);
    CHECK_FALSE(split.second_unique);
    CHECK_THROWS_AS(exact_intersection_ownership(a, b, owner, other, 1), ResourceLimit);
    CHECK_THROWS_AS(exact_intersection_ownership(a, b, owner, other, 4000000,
        [] { throw std::runtime_error("cancel"); }), std::runtime_error);
}

TEST_CASE("Analysis interface variants preserve the union and holes without selecting a physical priority", "[ThermoelasticExactSections]")
{
    auto a = box(0, 502), native_a = box(0, 501);
    const Polygon hole{{100,200},{100,800},{200,800},{200,200}};
    a.holes.push_back(hole);
    native_a.holes.push_back(hole);
    const ExPolygons b{box(498, 1000)}, native_b{box(499, 1000)};
    const auto variants = partition_uncertain_interface({a}, b, {native_a}, native_b);
    REQUIRE(variants.favor_first[0].size() == 1);
    REQUIRE(variants.favor_second[0].size() == 1);
    CHECK(variants.favor_first[0].front().holes.size() == 1);
    CHECK(variants.favor_second[0].front().holes.size() == 1);
    CHECK(exact_section_empty(variants.favor_first[0], variants.favor_first[1], ExactSectionOperation::Intersection).empty);
    CHECK(exact_section_empty(variants.favor_second[0], variants.favor_second[1], ExactSectionOperation::Intersection).empty);
    CHECK_FALSE(exact_section_empty(variants.favor_first[0], variants.favor_second[0], ExactSectionOperation::SymmetricDifference).empty);
    ExPolygons before{a};
    before.insert(before.end(), b.begin(), b.end());
    for (const auto* pair : {&variants.favor_first, &variants.favor_second}) {
        auto after = (*pair)[0];
        after.insert(after.end(), (*pair)[1].begin(), (*pair)[1].end());
        CHECK(exact_section_empty(before, after, ExactSectionOperation::SymmetricDifference).empty);
    }
    CHECK(variants.interface_allowance_nm == 7);
    CHECK(variants.representation_allowance_nm == 1);
    CHECK_THROWS_AS(partition_uncertain_interface({a}, b, {native_a}, native_b, 7, 1), ResourceLimit);
    CHECK_THROWS_AS(partition_uncertain_interface({a}, b, {native_a}, native_b, 7, 4000000,
        [] { throw std::runtime_error("cancel"); }), std::runtime_error);
    CHECK_THROWS_AS(partition_uncertain_interface({box(0,510)}, {box(490,1000)}, {native_a}, native_b), std::invalid_argument);
    CHECK_THROWS_AS(partition_uncertain_interface({a}, b, {box(0,499)}, {box(501,1000)}), std::invalid_argument);
}

TEST_CASE("Exact interface budgets apply to separate operands with cumulative work", "[ThermoelasticExactSections]")
{
    const ExPolygons a{box(0,502)}, b{box(498,1000)};
    const auto expected=partition_uncertain_interface(a,b,{box(0,501)},{box(499,1000)});
    // Same-material input masks may overlap. Their unions are unchanged, but
    // one operand exceeds1024, while the independent work ceiling still
    // applies to constructing all operands and every later Boolean operation.
    const ExPolygons native_a(320,box(0,501)), native_b(128,box(499,1000));
    const auto actual=partition_uncertain_interface(a,b,native_a,native_b);
    for (std::size_t i=0;i<2;++i) {
        CHECK(exact_section_empty(actual.favor_first[i],expected.favor_first[i],
            ExactSectionOperation::SymmetricDifference).empty);
        CHECK(exact_section_empty(actual.favor_second[i],expected.favor_second[i],
            ExactSectionOperation::SymmetricDifference).empty);
    }
    CHECK(actual.work>expected.work);
    CHECK_THROWS_AS(partition_uncertain_interface(a,b,native_a,native_b,7,actual.work-1),ResourceLimit);
    CHECK_THROWS_AS(partition_uncertain_interface(a,b,ExPolygons(513,box(0,501)),native_b),ResourceLimit);
}

TEST_CASE("Local native reconstruction ignores remote boundary complexity without rounding its interface", "[ThermoelasticExactSections]")
{
    const std::array<ExPolygons,2> envelopes{{{box(0,502)},{box(498,1000)}}};
    const auto expected=partition_uncertain_interface({}, {}, {box(0,501)}, {box(499,1000)},
        7,4000000,{},&envelopes,NativeWitnessFill::NonzeroWinding);
    Polygon remote;
    remote.points={{-10000,0},{501,0},{501,1000}};
    for (coord_t i=1;i<=1200;++i)
        remote.points.emplace_back(501-i*10501/1200,1000);
    // A second positive lobe meets the main contour at one remote point.
    // This is a valid native nonzero chain, but not a simple polygon ring.
    remote.points.emplace_back(-10000,0);
    remote.points.emplace_back(-11000,-1000);
    remote.points.emplace_back(-9000,-1000);
    const auto actual=partition_uncertain_interface({}, {}, {ExPolygon(remote)}, {box(499,1000)},
        7,4000000,{},&envelopes,NativeWitnessFill::NonzeroWinding);
    for (std::size_t i=0;i<2;++i) {
        CHECK(exact_section_empty(actual.favor_first[i],expected.favor_first[i],
            ExactSectionOperation::SymmetricDifference).empty);
        CHECK(exact_section_empty(actual.favor_second[i],expected.favor_second[i],
            ExactSectionOperation::SymmetricDifference).empty);
    }
    CHECK(actual.work<=4000000);
}

TEST_CASE("Rational native cropping agrees with complete simple masks and preserves holes", "[ThermoelasticExactSections]")
{
    ExPolygon native(Polygon{{-10000,-10000},{501,0},{501,1000},{-10000,10000}});
    native.holes.emplace_back(Points{{499,200},{499,800},{500,800},{500,200}});
    const std::array<ExPolygons,2> envelopes{{{box(0,502)},{box(498,1000)}}};
    const auto whole=partition_uncertain_interface({}, {}, {native}, {box(499,1000)},
        7,4000000,{},&envelopes,NativeWitnessFill::ValidatedSimple);
    const auto cropped=partition_uncertain_interface({}, {}, {native}, {box(499,1000)},
        7,4000000,{},&envelopes,NativeWitnessFill::NonzeroWinding);
    for (std::size_t i=0;i<2;++i) {
        CHECK(exact_section_empty(cropped.favor_first[i],whole.favor_first[i],
            ExactSectionOperation::SymmetricDifference).empty);
        CHECK(exact_section_empty(cropped.favor_second[i],whole.favor_second[i],
            ExactSectionOperation::SymmetricDifference).empty);
    }
}

TEST_CASE("Interface operand reuse retains boundary and hole validation", "[ThermoelasticExactSections]")
{
    const ExPolygons a{box(0, 502)}, b{box(498, 1000)}, native_b{box(499, 1000)};
    auto reversed = box(0, 501);
    reversed.contour.reverse();
    auto misplaced_hole = box(0, 501);
    misplaced_hole.holes.emplace_back(Points{{700, 200}, {700, 800}, {800, 800}, {800, 200}});
    const ExPolygons malformed[] = {
        {ExPolygon(Polygon{{0, 0}, {1, 1}})},
        {reversed},
        {misplaced_hole},
        {box(coord_t(10000000001LL), coord_t(10000000101LL))}
    };
    for (const auto& invalid : malformed) {
        CHECK_THROWS_AS(partition_uncertain_interface(a, b, invalid, native_b), std::invalid_argument);
        const std::array<ExPolygons, 2> envelopes{invalid, b};
        CHECK_THROWS_AS(partition_uncertain_interface({}, {}, {box(0, 501)}, native_b,
            7, 4000000, {}, &envelopes), std::invalid_argument);
    }
}

TEST_CASE("Native interface winding does not relax deposited envelope validity", "[ThermoelasticExactSections]")
{
    const ExPolygons native_a{box(0,501)}, native_b{box(499,1000)};
    const std::array<ExPolygons,2> envelopes{ExPolygons{box(0,502)},ExPolygons{box(498,1000)}};
    auto reversed = native_a;
    reversed.front().contour.reverse();
    const auto expected = partition_uncertain_interface({}, {},native_a,native_b,7,4000000,{},&envelopes);
    const auto actual = partition_uncertain_interface({}, {},reversed,native_b,7,4000000,{},&envelopes,
        NativeWitnessFill::NonzeroWinding);
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(exact_section_empty(expected.favor_first[i],actual.favor_first[i],ExactSectionOperation::SymmetricDifference).empty);
        CHECK(exact_section_empty(expected.favor_second[i],actual.favor_second[i],ExactSectionOperation::SymmetricDifference).empty);
    }
    auto invalid_envelopes = envelopes;
    invalid_envelopes[0].front().contour.reverse();
    CHECK_THROWS_AS(partition_uncertain_interface({}, {},native_a,native_b,7,4000000,{},&invalid_envelopes,
        NativeWitnessFill::NonzeroWinding),std::invalid_argument);
    reversed.push_back(native_a.front());
    CHECK_THROWS_AS(partition_uncertain_interface({}, {},reversed,native_b,7,4000000,{},&envelopes,
        NativeWitnessFill::NonzeroWinding),std::invalid_argument);
}

TEST_CASE("Unassigned deposited overlap reports its area without filling native void", "[ThermoelasticExactSections]")
{
    const std::array<ExPolygons,2> envelopes{ExPolygons{box(0,502)},ExPolygons{box(498,1000)}};
    try {
        (void)partition_uncertain_interface({}, {},{box(0,499)},{box(501,1000)},7,4000000,{},&envelopes,
            NativeWitnessFill::NonzeroWinding);
        FAIL("A deposited conflict without native material ownership must remain unresolved.");
    } catch (const std::invalid_argument& error) {
        CHECK(std::string(error.what()).find("area_nm2=2000.000000") != std::string::npos);
        CHECK(std::string(error.what()).find("Entire uncovered region lies in both native boundary bands") != std::string::npos);
    }
}

TEST_CASE("Native void outside the boundary bands remains diagnostically distinct", "[ThermoelasticExactSections]")
{
    const std::array<ExPolygons,2> envelopes{ExPolygons{box(0,510)},ExPolygons{box(490,1000)}};
    try {
        (void)partition_uncertain_interface({}, {},{box(0,490)},{box(510,1000)},7,4000000,{},&envelopes,
            NativeWitnessFill::NonzeroWinding);
        FAIL("Native void must not receive an arbitrary material.");
    } catch (const std::invalid_argument& error) {
        CHECK(std::string(error.what()).find("area_nm2=20000.000000") != std::string::npos);
        CHECK(std::string(error.what()).find("Analysis interface band is unverified") != std::string::npos);
    }
}

TEST_CASE("Thin deposited overlap outside native masks remains unassigned", "[ThermoelasticExactSections]")
{
    const std::array<ExPolygons,2> envelopes{ExPolygons{box(0,502)},ExPolygons{box(498,1000)}};
    try {
        (void)partition_uncertain_interface({}, {},{box(0,400)},{box(600,1000)},7,4000000,{},&envelopes,
            NativeWitnessFill::NonzeroWinding);
        FAIL("A thin bead overlap does not establish its material owner.");
    } catch (const std::invalid_argument& error) {
        const std::string message(error.what());
        CHECK(message.find("Native interface masks do not cover")!=std::string::npos);
        CHECK(message.find("Deposited boundary band 0 certified")!=std::string::npos);
        CHECK(message.find("Deposited boundary band 1 certified")!=std::string::npos);
    }
}

TEST_CASE("Single-sided boundary analysis requires one provably unique nearby native material", "[ThermoelasticExactSections]")
{
    const std::array<ExPolygons,2> envelopes{ExPolygons{box(0,502)},ExPolygons{box(498,1000)}};
    const ExPolygons first{box(0,400)},second{box(501,1000)};
    CHECK_THROWS_AS(partition_uncertain_interface({}, {},first,second,7,4000000,{},&envelopes),std::invalid_argument);
    for (const bool reverse : {false,true}) {
        const std::array<ExPolygons,2> inputs{envelopes[reverse ? 1 : 0],envelopes[reverse ? 0 : 1]};
        const auto result=partition_uncertain_interface({}, {},reverse ? second : first,reverse ? first : second,
            7,4000000,{},&inputs,NativeWitnessFill::ValidatedSimple,true);
        CHECK_THAT(result.single_sided_native_boundary_area_m2,Catch::Matchers::WithinRel(3000e-18,1e-12));
        for (const auto* variant : {&result.favor_first,&result.favor_second}) {
            CHECK(exact_section_empty((*variant)[reverse ? 1 : 0],{box(0,498)},ExactSectionOperation::SymmetricDifference).empty);
            CHECK(exact_section_empty((*variant)[reverse ? 0 : 1],{box(498,1000)},ExactSectionOperation::SymmetricDifference).empty);
        }
    }
    // Near both boundaries is ambiguous; near neither has no supporting mask.
    CHECK_THROWS_AS(partition_uncertain_interface({}, {},{box(0,499)},{box(501,1000)},7,4000000,{},&envelopes,
        NativeWitnessFill::ValidatedSimple,true),std::invalid_argument);
    CHECK_THROWS_AS(partition_uncertain_interface({}, {},{box(0,400)},{box(600,1000)},7,4000000,{},&envelopes,
        NativeWitnessFill::ValidatedSimple,true),std::invalid_argument);
}

TEST_CASE("Explicit dual boundary sensitivity preserves deposited geometry in both material assignments", "[ThermoelasticExactSections]")
{
    const std::array<ExPolygons,2> envelopes{ExPolygons{box(0,502)},ExPolygons{box(498,1000)}};
    const ExPolygons first{box(0,499)},second{box(501,1000)};
    CHECK_THROWS_AS(exact_partition_uncertain_interface(first,second,envelopes,ExactInterfaceChoice::FavorFirst),std::invalid_argument);
    CHECK_THROWS_AS(exact_partition_uncertain_interface(first,second,envelopes,ExactInterfaceChoice::FavorFirst,
        7,4000000,{},NativeWitnessFill::ValidatedSimple,true),std::invalid_argument);
    const auto whole=make_exact_section({box(0,1000)}).section;
    CHECK(exact_section_empty(first,second,ExactSectionOperation::Intersection).empty);
    for (const auto choice:{ExactInterfaceChoice::FavorFirst,ExactInterfaceChoice::FavorSecond}) {
        const auto partition=exact_partition_uncertain_interface(first,second,envelopes,choice,
            7,4000000,{},NativeWitnessFill::ValidatedSimple,false,true);
        const auto split=choice==ExactInterfaceChoice::FavorFirst ? 501 : 499;
        const auto expected_first=make_exact_section({box(0,split)}).section;
        const auto expected_second=make_exact_section({box(split,1000)}).section;
        CHECK(boolean_exact_sections(partition.materials[0],expected_first,ExactSectionOperation::SymmetricDifference).section.empty());
        CHECK(boolean_exact_sections(partition.materials[1],expected_second,ExactSectionOperation::SymmetricDifference).section.empty());
        CHECK(boolean_exact_sections(partition.materials[0],partition.materials[1],ExactSectionOperation::Intersection).section.empty());
        const auto occupied=boolean_exact_sections(partition.materials[0],partition.materials[1],ExactSectionOperation::Union).section;
        CHECK(boolean_exact_sections(occupied,whole,ExactSectionOperation::SymmetricDifference).section.empty());
        CHECK_THAT(partition.dual_native_boundary_area_m2,Catch::Matchers::WithinRel(2000e-18,1e-12));
        CHECK_THAT(partition.single_sided_native_boundary_area_m2,Catch::Matchers::WithinAbs(0.,0.));
        CHECK(partition.representation_allowance_nm==0);
        CHECK(partition.interface_allowance_nm==7);
    }
    CHECK_THROWS_AS(exact_partition_uncertain_interface({box(0,490)},second,envelopes,ExactInterfaceChoice::FavorFirst,
        7,4000000,{},NativeWitnessFill::ValidatedSimple,false,true),std::invalid_argument);
    // Every corner touches both native hole boundaries, but the square's
    // interior is far outside both bands. Vertex-only admission is invalid.
    ExPolygon surrounding(Polygon{{-100,-100},{1100,-100},{1100,1100},{-100,1100}});
    surrounding.holes.push_back(Polygon{{0,0},{0,1000},{1000,1000},{1000,0}});
    const std::array<ExPolygons,2> square{ExPolygons{box(0,1000)},ExPolygons{box(0,1000)}};
    CHECK_THROWS_AS(exact_partition_uncertain_interface({surrounding},{surrounding},square,ExactInterfaceChoice::FavorFirst,
        7,4000000,{},NativeWitnessFill::ValidatedSimple,false,true),std::invalid_argument);
    CHECK_THROWS_AS(exact_partition_uncertain_interface(first,second,envelopes,ExactInterfaceChoice::FavorFirst,
        6,4000000,{},NativeWitnessFill::ValidatedSimple,false,true),std::invalid_argument);
    CHECK_THROWS_AS(exact_partition_uncertain_interface(first,second,envelopes,ExactInterfaceChoice::FavorFirst,
        7,1,{},NativeWitnessFill::ValidatedSimple,false,true),ResourceLimit);
    CHECK_THROWS_AS(exact_partition_uncertain_interface(first,second,envelopes,ExactInterfaceChoice::FavorFirst,
        7,4000000,[] { throw std::runtime_error("cancel dual sensitivity"); },NativeWitnessFill::ValidatedSimple,false,true),std::runtime_error);
}

TEST_CASE("Single sided boundary failures distinguish unsupported and ambiguous fragments", "[ThermoelasticExactSections]")
{
    const std::array<ExPolygons,2> envelopes{ExPolygons{box(0,502)},ExPolygons{box(498,1000)}};
    for (const bool separated:{false,true}) {
        const ExPolygons first{box(0,separated ? 400 : 499)},second{box(separated ? 600 : 501,1000)};
        try {
            exact_partition_uncertain_interface(first,second,envelopes,ExactInterfaceChoice::FavorFirst,
                7,4000000,{},NativeWitnessFill::ValidatedSimple,true);
            FAIL("Unowned boundary fragments must remain rejected.");
        } catch (const std::invalid_argument& error) {
            const std::string message=error.what();
            CHECK(message.find(separated ? "separation=away-from-both" : "separation=not-separated-from-either")!=std::string::npos);
            CHECK(message.find(separated ? "away_first=true" : "away_first=false")!=std::string::npos);
            CHECK(message.find(separated ? "away_second=true" : "away_second=false")!=std::string::npos);
            CHECK(message.find("allowance_nm=7")!=std::string::npos);
            CHECK(message.find("native_fill=validated-simple")!=std::string::npos);
            CHECK(message.find("variant=favor-first")!=std::string::npos);
            CHECK(message.find("fragment_area_nm2=")!=std::string::npos);
            CHECK(message.find("fragment_holes=0")!=std::string::npos);
            CHECK(message.find("fragment_vertices=4")!=std::string::npos);
            CHECK(message.find("fragment_bounds_nm=")!=std::string::npos);
            CHECK(message.find("fragment_first_vertices_nm=")!=std::string::npos);
            CHECK(message.find("Native interface masks do not cover the complete deposited conflict")!=std::string::npos);
            CHECK(message.find("used_work=")!=std::string::npos);
            CHECK(message.find("max_work=4000000")!=std::string::npos);
            CHECK(message.size()<4096);
        }
    }
}

TEST_CASE("Interface rounding certifies long exterior wedges against their exact variant boundary", "[ThermoelasticExactSections]")
{
    // Splitting the slanted lower edge at x=501 creates a rational ordinate
    // 501/502. Rounding it moves a thin wedge reaching far from the interface,
    // but the entire wedge remains within one nm of the original exact edge.
    const ExPolygons a{ExPolygon(Polygon{{0,0},{502,1},{502,1000},{0,1000}})};
    const auto variants = partition_uncertain_interface(a, {box(498,1000)}, {box(0,501)}, {box(499,1000)});
    CHECK(variants.interface_allowance_nm == 7);
    CHECK(variants.representation_allowance_nm == 1);
    CHECK(exact_section_empty(variants.favor_first[0], variants.favor_first[1], ExactSectionOperation::Intersection).empty);
    CHECK(exact_section_empty(variants.favor_second[0], variants.favor_second[1], ExactSectionOperation::Intersection).empty);
}

TEST_CASE("Two-axis interface rounding retains certified disjoint variants", "[ThermoelasticExactSections]")
{
    const ExPolygons a{ExPolygon(Polygon{{0,0},{502,1},{502,4},{0,4}})};
    const ExPolygons b{ExPolygon(Polygon{{498,0},{1000,0},{1000,4},{498,4}})};
    const ExPolygons native_a{ExPolygon(Polygon{{521,-31},{479,32},{179,-168},{221,-231}})};
    const ExPolygons native_b{ExPolygon(Polygon{{479,31},{521,-32},{821,168},{779,231}})};
    const auto variants = partition_uncertain_interface(a, b, native_a, native_b);
    CHECK(variants.representation_allowance_nm == 1);
    CHECK(exact_section_empty(variants.favor_first[0], variants.favor_first[1], ExactSectionOperation::Intersection).empty);
    CHECK(exact_section_empty(variants.favor_second[0], variants.favor_second[1], ExactSectionOperation::Intersection).empty);
}

TEST_CASE("Interface topology matching preserves separated rational components", "[ThermoelasticExactSections]")
{
    ExPolygons a, b, native_a, native_b;
    constexpr std::size_t components = 30;
    for (std::size_t i = 0; i < components; ++i) {
        const coord_t x = coord_t(i) * 2000;
        a.emplace_back(Polygon{{x, 0}, {x + 502, 1}, {x + 502, 1000}, {x, 1000}});
        b.push_back(box(x + 498, x + 1000));
        native_a.push_back(box(x, x + 501));
        native_b.push_back(box(x + 499, x + 1000));
    }
    const auto variants = partition_uncertain_interface(a, b, native_a, native_b);
    for (const auto* pair : {&variants.favor_first, &variants.favor_second}) {
        CHECK((*pair)[0].size() == components);
        CHECK((*pair)[1].size() == components);
        CHECK(exact_section_empty((*pair)[0], (*pair)[1], ExactSectionOperation::Intersection).empty);
    }
}

TEST_CASE("Interface topology witnesses preserve holes during rational rounding", "[ThermoelasticExactSections]")
{
    ExPolygon a(Polygon{{0, 0}, {600, 1}, {600, 1000}, {0, 1000}});
    a.holes.emplace_back(Points{{100, 200}, {100, 800}, {200, 800}, {200, 200}});
    auto native_a = box(0, 501);
    native_a.holes = a.holes;
    const std::array<ExPolygons, 2> envelopes{{{a}, {box(400, 1000)}}};
    const auto variants = partition_uncertain_interface({}, {}, {native_a}, {box(499, 1000)},
        7, 4000000, {}, &envelopes);
    const ExPolygons inside_hole{ExPolygon(Polygon{{120, 300}, {180, 300}, {180, 700}, {120, 700}})};
    for (const auto* pair : {&variants.favor_first, &variants.favor_second}) {
        REQUIRE((*pair)[0].size() == 1);
        CHECK((*pair)[0].front().holes.size() == 1);
        CHECK(exact_section_empty((*pair)[0], inside_hole, ExactSectionOperation::Intersection).empty);
        CHECK(exact_section_empty((*pair)[0], (*pair)[1], ExactSectionOperation::Intersection).empty);
    }
}

TEST_CASE("Interface rounding removes only consecutive collapsed edge images", "[ThermoelasticExactSections]")
{
    // The native cut intersects the short diagonal less than half a grid unit
    // from its endpoint. Both exact vertices therefore have one rounded image.
    const ExPolygons a{ExPolygon(Polygon{{0,0},{502,1},{502,1000},{0,1000}})};
    const ExPolygons b{box(498,1000)};
    const ExPolygons native_a{ExPolygon(Polygon{{0,0},{501,0},{502,2},{502,1000},{0,1000}})};
    const auto variants = partition_uncertain_interface(a, b, native_a, {box(499,1000)});
    for (const auto* pair : {&variants.favor_first, &variants.favor_second}) {
        CHECK(exact_section_empty((*pair)[0], (*pair)[1], ExactSectionOperation::Intersection).empty);
        for (const auto& material : *pair)
            for (const auto& polygon : material) {
                REQUIRE(polygon.contour.points.size() >= 3);
                for (std::size_t i = 0; i < polygon.contour.points.size(); ++i)
                    CHECK(polygon.contour.points[i] != polygon.contour.points[(i + 1) % polygon.contour.points.size()]);
            }
    }
    CHECK(variants.representation_allowance_nm == 1);
}

TEST_CASE("Complete interface construction preserves deposited exclusive regions without reclipping claims", "[ThermoelasticExactSections]")
{
    const std::array<ExPolygons, 2> envelopes{{{box(0, 600)}, {box(400, 1000)}}};
    const auto variants = partition_uncertain_interface({}, {}, {box(0,501)}, {box(499,1000)},
        7, 4000000, {}, &envelopes);
    CHECK(exact_section_empty(variants.favor_first[0],{box(0,501)},ExactSectionOperation::SymmetricDifference).empty);
    CHECK(exact_section_empty(variants.favor_first[1],{box(501,1000)},ExactSectionOperation::SymmetricDifference).empty);
    CHECK(exact_section_empty(variants.favor_second[0],{box(0,499)},ExactSectionOperation::SymmetricDifference).empty);
    CHECK(exact_section_empty(variants.favor_second[1],{box(499,1000)},ExactSectionOperation::SymmetricDifference).empty);
    const auto disjoint=partition_uncertain_interface({}, {}, {box(0,500)}, {box(500,1000)},7,4000000,{},&envelopes);
    CHECK(exact_section_empty(disjoint.favor_first[0],disjoint.favor_second[0],ExactSectionOperation::SymmetricDifference).empty);
    CHECK(exact_section_empty(disjoint.favor_first[1],disjoint.favor_second[1],ExactSectionOperation::SymmetricDifference).empty);
    for (const auto* pair : {&variants.favor_first, &variants.favor_second}) {
        CHECK(exact_section_empty((*pair)[0], (*pair)[1], ExactSectionOperation::Intersection).empty);
        auto combined = (*pair)[0];
        combined.insert(combined.end(), (*pair)[1].begin(), (*pair)[1].end());
        CHECK(exact_section_empty(combined, {box(0,1000)}, ExactSectionOperation::SymmetricDifference).empty);
    }
    CHECK_THROWS_AS(partition_uncertain_interface({}, {}, {box(0,499)}, {box(501,1000)},
        7, 4000000, {}, &envelopes), std::invalid_argument);
    CHECK_THROWS_AS(partition_uncertain_interface({}, {}, {box(0,510)}, {box(490,1000)},
        7, 4000000, {}, &envelopes), std::invalid_argument);
}

TEST_CASE("Selected exact partitions match both legacy integer variants", "[ThermoelasticExactSections]")
{
    const ExPolygons native_a{box(0,501)},native_b{box(499,1000)};
    const std::array<ExPolygons,2> envelopes{ExPolygons{box(0,510)},ExPolygons{box(490,1000)}};
    const auto legacy=partition_uncertain_interface({}, {}, native_a,native_b,7,4000000,{},&envelopes);
    const auto choice=GENERATE(ExactInterfaceChoice::FavorFirst,ExactInterfaceChoice::FavorSecond);
    const auto exact=exact_partition_uncertain_interface(native_a,native_b,envelopes,choice);
    const auto& expected=choice==ExactInterfaceChoice::FavorFirst ? legacy.favor_first : legacy.favor_second;
    for (unsigned material=0;material<2;++material) {
        CHECK(boolean_exact_sections(exact.materials[material],make_exact_section(expected[material]).section,
            ExactSectionOperation::SymmetricDifference).section.empty());
        const auto exclusive=boolean_exact_sections(make_exact_section(envelopes[material]).section,
            make_exact_section(envelopes[1-material]).section,ExactSectionOperation::Difference);
        CHECK(boolean_exact_sections(exclusive.section,exact.materials[material],ExactSectionOperation::Difference).section.empty());
    }
    CHECK(boolean_exact_sections(exact.materials[0],exact.materials[1],ExactSectionOperation::Intersection).section.empty());
    const auto occupied=boolean_exact_sections(exact.materials[0],exact.materials[1],ExactSectionOperation::Union);
    CHECK(boolean_exact_sections(occupied.section,make_exact_section({box(0,1000)}).section,
        ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK_THROWS_AS(exact_partition_uncertain_interface(native_a,native_b,envelopes,choice,7,exact.work-1),ResourceLimit);
    CHECK_THROWS_AS(exact_partition_uncertain_interface({box(0,499)},{box(501,1000)},envelopes,choice),std::invalid_argument);
    CHECK_THROWS_AS(exact_partition_uncertain_interface({box(0,510)},{box(490,1000)},envelopes,choice),std::invalid_argument);
    unsigned calls=0;
    CHECK_THROWS_AS(exact_partition_uncertain_interface(native_a,native_b,envelopes,choice,7,4000000,
        [&] { if (++calls==100) throw std::runtime_error("cancel"); }),std::runtime_error);
}

TEST_CASE("Exact interface partitions preserve rational occupied holes for both material choices", "[ThermoelasticExactSections]")
{
    auto first=box(0,1000), second=first;
    first.holes.emplace_back(Polygon{{400,400},{400,600},{600,400}});
    second.holes.emplace_back(Polygon{{400,400},{400,600},{600,500}});
    const std::array<ExPolygons,2> envelopes{ExPolygons{first},ExPolygons{second}};
    const auto deposited=make_exact_section({first,second});
    const std::array<coord_t,3> start{0,0,0},size{1000,1000,10};
    std::array<ExactSection,2> first_materials;
    for (const auto choice : {ExactInterfaceChoice::FavorFirst,ExactInterfaceChoice::FavorSecond}) {
        const auto partition=exact_partition_uncertain_interface({box(0,501)},{box(499,1000)},envelopes,choice);
        first_materials[choice==ExactInterfaceChoice::FavorFirst ? 0 : 1]=partition.materials[0];
        CHECK(partition.interface_allowance_nm==7);
        CHECK(partition.representation_allowance_nm==0);
        CHECK(partition.materials[0].vertex_count()+partition.materials[1].vertex_count()<=4096);
        CHECK(boolean_exact_sections(partition.materials[0],partition.materials[1],ExactSectionOperation::Intersection).section.empty());
        // The exact occupied hole is a triangle with third corner
        // (1600/3,1400/3), not an integer-grid approximation.
        const auto occupied=boolean_exact_sections(partition.materials[0],partition.materials[1],ExactSectionOperation::SymmetricDifference);
        CHECK(boolean_exact_sections(occupied.section,deposited.section,ExactSectionOperation::SymmetricDifference).section.empty());
        const auto moments=exact_section_volume_moments(occupied.section,start,size,0,10);
        CHECK_THAT(moments.values[0],Catch::Matchers::WithinRel(double((1000000.L-40000.L/3.L)*10.L*1e-27L),1e-13));
        const auto expected=exact_section_volume_moments(deposited.section,start,size,0,10);
        for (std::size_t i=0;i<27;++i)
            CHECK_THAT(moments.values[i],Catch::Matchers::WithinRel(expected.values[i],1e-13));
        CHECK_THROWS_AS(exact_partition_uncertain_interface({box(0,501)},{box(499,1000)},envelopes,choice,7,partition.work-1),ResourceLimit);
    }
    CHECK_FALSE(boolean_exact_sections(first_materials[0],first_materials[1],ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK_THROWS_AS(exact_partition_uncertain_interface({box(0,499)},{box(501,1000)},envelopes,ExactInterfaceChoice::FavorFirst),std::invalid_argument);
    CHECK_THROWS_AS(exact_partition_uncertain_interface({box(0,510)},{box(490,1000)},envelopes,ExactInterfaceChoice::FavorFirst),std::invalid_argument);
    CHECK_THROWS_AS(exact_partition_uncertain_interface({box(0,501)},{box(499,1000)},envelopes,ExactInterfaceChoice::FavorFirst,
        7,4000000,[] { throw std::runtime_error("cancel exact partition"); }),std::runtime_error);
    CHECK_THROWS_AS(exact_partition_uncertain_interface({box(0,501)},{box(499,1000)},envelopes,static_cast<ExactInterfaceChoice>(9)),std::invalid_argument);
}

TEST_CASE("Exact interface partitions retain single sided sensitivity without rounding allowance", "[ThermoelasticExactSections]")
{
    const std::array<ExPolygons,2> envelopes{ExPolygons{box(0,502)},ExPolygons{box(498,1000)}};
    const ExPolygons first{box(0,400)},second{box(501,1000)};
    CHECK_THROWS_AS(exact_partition_uncertain_interface(first,second,envelopes,ExactInterfaceChoice::FavorSecond),std::invalid_argument);
    const auto partition=exact_partition_uncertain_interface(first,second,envelopes,ExactInterfaceChoice::FavorSecond,
        7,4000000,{},NativeWitnessFill::ValidatedSimple,true);
    CHECK(partition.interface_allowance_nm==7);
    CHECK(partition.representation_allowance_nm==0);
    CHECK_THAT(partition.single_sided_native_boundary_area_m2,Catch::Matchers::WithinRel(3000e-18,1e-12));
    CHECK(boolean_exact_sections(partition.materials[0],make_exact_section({box(0,498)}).section,
        ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK(boolean_exact_sections(partition.materials[1],make_exact_section({box(498,1000)}).section,
        ExactSectionOperation::SymmetricDifference).section.empty());
}

TEST_CASE("Interface diagnostics propagate cancellation inside uncovered boundary proofs", "[ThermoelasticExactSections]")
{
    struct Cancel : std::invalid_argument { Cancel() : std::invalid_argument("cancel uncovered band") {} };
    const std::array<ExPolygons,2> envelopes{ExPolygons{box(0,502)},ExPolygons{box(498,1000)}};
    for (const bool distant : {false,true}) {
        const ExPolygons first{box(0,distant ? 400 : 499)},second{box(distant ? 600 : 501,1000)};
        std::size_t visits=0;
        try {
            (void)exact_partition_uncertain_interface(first,second,envelopes,ExactInterfaceChoice::FavorFirst,
                7,4000000,[&] { ++visits; });
            FAIL("Uncovered native material must remain unresolved.");
        } catch (const std::invalid_argument& error) {
            const std::string message(error.what());
            CHECK(message.find(distant ? "Deposited boundary band 1 certified" :
                "Entire uncovered region lies in both native boundary bands")!=std::string::npos);
        }
        REQUIRE(visits>1);
        // The final callback is inside the uncovered-band proof: either the
        // second native band or the final deposited-band diagnostic. Throwing
        // invalid_argument there must not be mistaken for geometry failure.
        std::size_t current=0;
        CHECK_THROWS_AS(exact_partition_uncertain_interface(first,second,envelopes,ExactInterfaceChoice::FavorFirst,
            7,4000000,[&] { if (++current==visits) throw Cancel(); }),Cancel);
        CHECK(current==visits);
    }
}

TEST_CASE("Interface rounding refuses a subgrid gap that would join distinct material islands", "[ThermoelasticExactSections]")
{
    const ExPolygons envelope{ExPolygon(Polygon{{500,-100},{502,-100},{502,100},{500,100}})};
    const std::array<ExPolygons, 2> envelopes{envelope, envelope};
    const ExPolygons native_a{
        ExPolygon(Polygon{{498,-100},{508,-100},{508,0},{498,0}}),
        ExPolygon(Polygon{{498,0},{508,1},{508,100},{498,100}})};
    const ExPolygons native_b{ExPolygon(Polygon{{499,0},{503,0},{503,1},{499,1}})};
    // Inside the actual envelope the two A pieces are separated by .2-.4 nm.
    // Integer output must not silently join them or erase the intervening B.
    CHECK_THROWS_AS(partition_uncertain_interface({}, {}, native_a, native_b,
        7, 4000000, {}, &envelopes), std::invalid_argument);
}
