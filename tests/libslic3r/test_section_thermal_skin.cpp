#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/SectionThermalSkin.hpp"
#include "libslic3r/SectionThermalAssembly.hpp"
#include "libslic3r/SectionThermalExposure.hpp"
#include "libslic3r/SectionIntegration.hpp"
#include "libslic3r/ThermoelasticLocalSections.hpp"
#include "libslic3r/ThermoelasticLocalTopology.hpp"
#include "libslic3r/ThermoelasticSectionTopology.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticLayers.hpp"
#include "libslic3r/ThermoelasticExactSections.hpp"
#include "libslic3r/ThermoelasticSectionFootprint.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/libslic3r.h"
#include <Eigen/Core>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;

TEST_CASE("Exact streamed skin preserves rational tile seams and layer exposure", "[SectionThermalSkin]")
{
    const auto triangle=make_exact_section({ExPolygon(Polygon{{0,0},{5,0},{0,1}})}).section;
    const auto cutter=make_exact_section({ExPolygon(Polygon{{0,0},{4,0},{4,2},{0,2}})}).section;
    const auto left=boolean_exact_sections(triangle,cutter,ExactSectionOperation::Intersection).section;
    const auto right=boolean_exact_sections(triangle,cutter,ExactSectionOperation::Difference).section;
    const std::array<coord_t,3> start{0,0,0},size{10,10,10};
    for (const coord_t right_top:{coord_t(1),coord_t(2)}) {
        CAPTURE(right_top);
        SectionIntegration::Moments horizontal{},vertical{};
        SectionThermalSkin skin(
            [](std::size_t,std::size_t,const ExPolygons&,coord_t,bool) { FAIL("Native face in exact stream"); },
            [](std::size_t,std::size_t,const VerticalSectionPatch&) { FAIL("Native side in exact stream"); },
            {},{},{},
            [&](std::size_t,std::size_t,const ExactSection& face,coord_t z,bool) {
                const auto moment=exact_section_horizontal_surface_moments(face,start,size,z);
                for (std::size_t i=0;i<27;++i) horizontal[i]+=moment.values[i];
            },
            [&](std::size_t,std::size_t,const ExactSectionBoundary& span,coord_t bottom,coord_t top,const std::array<double,2>& normal) {
                CHECK_THAT(normal[0]*normal[0]+normal[1]*normal[1],WithinAbs(1.,1e-14));
                const auto moment=exact_section_vertical_surface_moments(span,start,size,bottom,top,ExactBoundaryCellOwnership::AllIncident);
                for (std::size_t i=0;i<27;++i) vertical[i]+=moment.values[i];
            });
        skin.begin_tile({{0,0},{0,0},{4,2},{}});
        skin.add_band({{{0.,2e-9,0,SectionFootprint(left)},0}});
        skin.end_tile();
        skin.begin_tile({{1,0},{4,0},{8,2},{}});
        skin.add_band({{{0.,double(right_top)*1e-9,1,SectionFootprint(right)},0}});
        skin.end_tile();
        skin.finish();
        CHECK_THAT(horizontal[0],WithinAbs(5e-18,1e-30));
        if (right_top==2) {
            const auto boundary=extract_exact_section_boundary(triangle).boundary;
            const auto expected=exact_section_vertical_surface_moments(boundary,start,size,0,2,ExactBoundaryCellOwnership::AllIncident);
            const auto bottom=exact_section_horizontal_surface_moments(triangle,start,size,0);
            const auto top=exact_section_horizontal_surface_moments(triangle,start,size,2);
            for (std::size_t i=0;i<27;++i) {
                CHECK_THAT(vertical[i],WithinAbs(expected.values[i],1e-29));
                CHECK_THAT(horizontal[i],WithinAbs(bottom.values[i]+top.values[i],1e-29));
            }
        } else {
            const double perimeter=6.+std::sqrt(26.);
            CHECK_THAT(vertical[0],WithinAbs((2.*perimeter-1.-std::sqrt(1.04)+.2)*1e-18,1e-29));
        }
    }
}

TEST_CASE("Exact horizontal skin tries matching materials without losing cross material coverage", "[SectionThermalSkin]")
{
    const auto first=make_exact_section({ExPolygon(Polygon{{0,0},{4,0},{0,4}})}).section;
    const auto second=make_exact_section({ExPolygon(Polygon{{4,0},{4,4},{0,4}})}).section;
    const auto whole=make_exact_section({ExPolygon(Polygon{{0,0},{4,0},{4,4},{0,4}})}).section;
    const std::array<coord_t,3> start{0,0,0},size{4,4,2};
    const auto run=[&](bool swap_materials,std::size_t budget) {
        std::array<SectionIntegration::Moments,2> exposed{};
        SectionSkinLimits limits;
        limits.max_work=budget;
        SectionThermalSkin skin(
            [](std::size_t,std::size_t,const ExPolygons&,coord_t,bool) { FAIL("Native face in exact stream"); },
            [](std::size_t,std::size_t,const VerticalSectionPatch&) { FAIL("Native side in exact stream"); },
            limits,{},{},
            [&](std::size_t,std::size_t,const ExactSection& face,coord_t z,bool) {
                const auto moments=exact_section_horizontal_surface_moments(face,start,size,z);
                for (std::size_t i=0;i<27;++i) exposed[0][i]+=moments.values[i];
            },
            [&](std::size_t,std::size_t,const ExactSectionBoundary& span,coord_t bottom,coord_t top,const std::array<double,2>&) {
                const auto moments=exact_section_vertical_surface_moments(span,start,size,bottom,top,
                    ExactBoundaryCellOwnership::AllIncident);
                for (std::size_t i=0;i<27;++i) exposed[1][i]+=moments.values[i];
            });
        // Both outward boxes overlap. The first opposing piece has the other
        // geometry; equal same-material coverage comes later in the band.
        // Keep the square inside the tile so seam compaction cannot affect work.
        skin.begin_tile({{0,0},{-1,-1},{5,5},{}});
        skin.add_band({{{0.,1e-9,0,SectionFootprint(first)},0},{{0.,1e-9,1,SectionFootprint(second)},0}});
        skin.add_band({{{1e-9,2e-9,swap_materials?0u:1u,SectionFootprint(second)},0},
                       {{1e-9,2e-9,swap_materials?1u:0u,SectionFootprint(first)},0}});
        skin.end_tile();
        skin.finish();
        return exposed;
    };
    constexpr std::size_t ceiling=100000;
    const auto bottom=exact_section_horizontal_surface_moments(whole,start,size,0);
    const auto top=exact_section_horizontal_surface_moments(whole,start,size,2);
    const auto sides=exact_section_vertical_surface_moments(extract_exact_section_boundary(whole).boundary,
        start,size,0,2,ExactBoundaryCellOwnership::AllIncident);
    for (const bool swapped:{false,true}) {
        const auto exposed=run(swapped,ceiling);
        for (std::size_t i=0;i<27;++i) {
            CHECK_THAT(exposed[0][i],WithinAbs(bottom.values[i]+top.values[i],1e-28));
            CHECK_THAT(exposed[1][i],WithinAbs(sides.values[i],1e-28));
        }
    }
    // Find the actual bounded work, without freezing exact-overlay counters.
    const auto minimum_budget=[&](bool swapped) {
        std::size_t low=1,high=ceiling;
        while (low<high) {
            const auto middle=low+(high-low)/2;
            try { run(swapped,middle); high=middle; }
            catch (const ResourceLimit&) { low=middle+1; }
        }
        return low;
    };
    const auto matching_budget=minimum_budget(false);
    CHECK(matching_budget<minimum_budget(true));
    REQUIRE(matching_budget>1);
    CHECK_THROWS_AS(run(false,matching_budget-1),ResourceLimit);
}

TEST_CASE("Exact skin retains native bands while enforcing cancellation and storage caps", "[SectionThermalSkin]")
{
    const auto h=[](std::size_t,std::size_t,const ExPolygons&,coord_t,bool) {};
    const auto v=[](std::size_t,std::size_t,const VerticalSectionPatch&) {};
    const auto eh=[](std::size_t,std::size_t,const ExactSection&,coord_t,bool) {};
    const auto ev=[](std::size_t,std::size_t,const ExactSectionBoundary&,coord_t,coord_t,const std::array<double,2>&) {};
    const std::vector<SectionPiece> band{{{0.,1e-9,0,{ExPolygon(Polygon{{0,0},{2,0},{2,2},{0,2}})}},0}};
    CHECK_THROWS_AS(SectionThermalSkin(h,v,{},{},{},eh,{}),std::invalid_argument);
    SectionSkinLimits limits;
    limits.max_live_points=10;
    SectionThermalSkin bounded(h,v,limits,{},{},eh,ev);
    bounded.begin_tile({{0,0},{0,0},{2,2},{}});
    CHECK_THROWS_AS(bounded.add_band(band),ResourceLimit);
    CHECK_THROWS_AS(bounded.finish(),std::invalid_argument);
    bool cancel=false;
    SectionThermalSkin cancelled(h,v,{},[&] { if (cancel) throw std::runtime_error("cancel exact skin"); },{},eh,ev);
    cancelled.begin_tile({{0,0},{0,0},{2,2},{}});
    cancel=true;
    CHECK_THROWS_AS(cancelled.add_band(band),std::runtime_error);
    SectionThermalSkin promoted(h,v,{},{},{},eh,ev);
    promoted.begin_tile({{0,0},{0,0},{2,2},{}});
    promoted.add_band(band);
    promoted.end_tile();
    CHECK_NOTHROW(promoted.finish());
}

TEST_CASE("Exact skin compacts adjacent seam bands before reaching its seam cap", "[SectionThermalSkin]")
{
    SectionSkinLimits limits;
    limits.max_seam_patches=3;
    const auto shape=make_exact_section({ExPolygon(Polygon{{0,1},{2,1},{2,3},{0,3}})}).section;
    double area=0.;
    SectionThermalSkin skin(
        [](std::size_t,std::size_t,const ExPolygons&,coord_t,bool) {},
        [](std::size_t,std::size_t,const VerticalSectionPatch&) {},limits,{},{},
        [](std::size_t,std::size_t,const ExactSection&,coord_t,bool) {},
        [&](std::size_t,std::size_t,const ExactSectionBoundary& span,coord_t bottom,coord_t top,const std::array<double,2>&) {
            area+=exact_section_vertical_surface_moments(span,{0,0,0},{4,4,4},bottom,top,
                ExactBoundaryCellOwnership::AllIncident).values[0];
        });
    skin.begin_tile({{0,0},{0,0},{2,4},{}});
    for (coord_t z=0;z<3;++z)
        skin.add_band({{{double(z)*1e-9,double(z+1)*1e-9,0,SectionFootprint(shape)},0}});
    skin.end_tile();
    skin.finish();
    CHECK_THAT(area,WithinAbs(24e-18,1e-29));
}

TEST_CASE("Exact skin preserves hole walls and unsupported undersides beneath a closing cap", "[SectionThermalSkin]")
{
    ExPolygon ring(Polygon{{0,0},{4,0},{4,4},{0,4}});
    ring.holes.emplace_back(Polygon{{1,1},{1,3},{3,3},{3,1}});
    const auto lower=make_exact_section({ring}).section;
    const auto upper=make_exact_section({ExPolygon(Polygon{{0,0},{4,0},{4,4},{0,4}})}).section;
    const auto hole=make_exact_section({ExPolygon(Polygon{{1,1},{3,1},{3,3},{1,3}})}).section;
    const std::array<coord_t,3> start{0,0,0},size{4,4,2};
    SectionIntegration::Moments underside{};
    double bed_area=0.,top_area=0.,side_area=0.,hole_wall_area=0.;
    std::size_t hole_walls=0;
    SectionThermalSkin skin(
        [](std::size_t,std::size_t,const ExPolygons&,coord_t,bool) { FAIL("Native face in exact stream"); },
        [](std::size_t,std::size_t,const VerticalSectionPatch&) { FAIL("Native side in exact stream"); },{},{},{},
        [&](std::size_t,std::size_t,const ExactSection& face,coord_t z,bool upward) {
            const auto moments=exact_section_horizontal_surface_moments(face,start,size,z);
            if (!upward && z==0) bed_area+=moments.values[0];
            else if (upward && z==2) top_area+=moments.values[0];
            else {
                CHECK_FALSE(upward);
                CHECK(z==1);
                for (std::size_t i=0;i<27;++i) underside[i]+=moments.values[i];
            }
        },
        [&](std::size_t,std::size_t,const ExactSectionBoundary& span,coord_t bottom,coord_t top,const std::array<double,2>& normal) {
            const auto moments=exact_section_vertical_surface_moments(span,start,size,bottom,top,ExactBoundaryCellOwnership::AllIncident);
            side_area+=moments.values[0];
            const auto box=*span.bounds();
            if (box.min.x()>=1 && box.max.x()<=3 && box.min.y()>=1 && box.max.y()<=3) {
                ++hole_walls;
                hole_wall_area+=moments.values[0];
                CHECK(bottom==0);
                CHECK(top==1);
                if (box.min.x()==box.max.x())
                    CHECK_THAT(normal[0],WithinAbs(box.min.x()==1 ? 1. : -1.,1e-14));
                else CHECK_THAT(normal[1],WithinAbs(box.min.y()==1 ? 1. : -1.,1e-14));
            }
        });
    skin.begin_tile({{0,0},{0,0},{4,4},{}});
    skin.add_band({{{0.,1e-9,0,SectionFootprint(lower)},0}});
    skin.add_band({{{1e-9,2e-9,1,SectionFootprint(upper)},0}});
    skin.end_tile();
    skin.finish();
    CHECK(hole_walls==4);
    CHECK_THAT(hole_wall_area,WithinAbs(8e-18,1e-29));
    CHECK_THAT(side_area,WithinAbs(40e-18,1e-29));
    CHECK_THAT(bed_area,WithinAbs(12e-18,1e-29));
    CHECK_THAT(top_area,WithinAbs(16e-18,1e-29));
    const auto expected=exact_section_horizontal_surface_moments(hole,start,size,1);
    for (std::size_t i=0;i<27;++i) CHECK_THAT(underside[i],WithinAbs(expected.values[i],1e-29));
}

TEST_CASE("Exact skin does not remove tile seam faces meeting only at a point", "[SectionThermalSkin]")
{
    double side_area=0.,seam_area=0.;
    SectionThermalSkin skin(
        [](std::size_t,std::size_t,const ExPolygons&,coord_t,bool) {},
        [](std::size_t,std::size_t,const VerticalSectionPatch&) {},{},{},{},
        [](std::size_t,std::size_t,const ExactSection&,coord_t,bool) {},
        [&](std::size_t,std::size_t,const ExactSectionBoundary& span,coord_t bottom,coord_t top,const std::array<double,2>&) {
            const auto area=exact_section_vertical_surface_moments(span,{0,0,0},{4,2,1},bottom,top,
                ExactBoundaryCellOwnership::AllIncident).values[0];
            side_area+=area;
            const auto box=*span.bounds();
            if (box.min.x()==2 && box.max.x()==2) seam_area+=area;
        });
    const auto left=make_exact_section({ExPolygon(Polygon{{0,0},{2,0},{2,1},{0,1}})}).section;
    const auto right=make_exact_section({ExPolygon(Polygon{{2,1},{4,1},{4,2},{2,2}})}).section;
    skin.begin_tile({{0,0},{0,0},{2,2},{}});
    skin.add_band({{{0.,1e-9,0,SectionFootprint(left)},0}});
    skin.end_tile();
    skin.begin_tile({{1,0},{2,0},{4,2},{}});
    skin.add_band({{{0.,1e-9,0,SectionFootprint(right)},1}});
    skin.end_tile();
    skin.finish();
    CHECK_THAT(side_area,WithinAbs(12e-18,1e-29));
    CHECK_THAT(seam_area,WithinAbs(2e-18,1e-29));
}

TEST_CASE("Exact skin callback failures permanently poison the stream", "[SectionThermalSkin]")
{
    const auto shape=make_exact_section({ExPolygon(Polygon{{1,1},{2,1},{2,2},{1,2}})}).section;
    const std::vector<SectionPiece> band{{{0.,1e-9,0,SectionFootprint(shape)},0}};
    for (const bool horizontal_failure:{false,true}) {
        CAPTURE(horizontal_failure);
        std::size_t failures=0;
        SectionThermalSkin skin(
            [](std::size_t,std::size_t,const ExPolygons&,coord_t,bool) {},
            [](std::size_t,std::size_t,const VerticalSectionPatch&) {},{},{},{},
            [&](std::size_t,std::size_t,const ExactSection&,coord_t,bool) {
                if (horizontal_failure) { ++failures; throw std::runtime_error("cancel horizontal callback"); }
            },
            [&](std::size_t,std::size_t,const ExactSectionBoundary&,coord_t,coord_t,const std::array<double,2>&) {
                ++failures; throw std::runtime_error("cancel vertical callback");
            });
        skin.begin_tile({{0,0},{0,0},{4,4},{}});
        CHECK_THROWS_AS(skin.add_band(band),std::runtime_error);
        CHECK_THROWS_AS(skin.add_band(band),std::invalid_argument);
        CHECK_THROWS_AS(skin.end_tile(),std::invalid_argument);
        CHECK_THROWS_AS(skin.finish(),std::invalid_argument);
        CHECK(failures==1);
    }
}

TEST_CASE("Exact skin distinguishes zero area horizontal contact from shared vertical faces", "[SectionThermalSkin]")
{
    const auto left=make_exact_section({ExPolygon(Polygon{{0,0},{1,0},{1,1},{0,1}})}).section;
    const auto right=make_exact_section({ExPolygon(Polygon{{1,0},{2,0},{2,1},{1,1}})}).section;
    for (const bool same_band:{false,true}) {
        CAPTURE(same_band);
        double horizontal=0.,vertical=0.;
        SectionThermalSkin skin(
            [](std::size_t,std::size_t,const ExPolygons&,coord_t,bool) {},
            [](std::size_t,std::size_t,const VerticalSectionPatch&) {},{},{},{},
            [&](std::size_t,std::size_t,const ExactSection& face,coord_t plane,bool) {
                horizontal+=exact_section_horizontal_surface_moments(face,{0,0,0},{2,2,2},plane).values[0];
            },
            [&](std::size_t,std::size_t,const ExactSectionBoundary& span,coord_t bottom,coord_t top,const std::array<double,2>&) {
                vertical+=exact_section_vertical_surface_moments(span,{0,0,0},{2,2,2},bottom,top,
                    ExactBoundaryCellOwnership::AllIncident).values[0];
            });
        skin.begin_tile({{0,0},{0,0},{2,2},{}});
        if (same_band) skin.add_band({{{0.,1e-9,0,SectionFootprint(left)},0},{{0.,1e-9,1,SectionFootprint(right)},0}});
        else {
            skin.add_band({{{0.,1e-9,0,SectionFootprint(left)},0}});
            skin.add_band({{{1e-9,2e-9,1,SectionFootprint(right)},0}});
        }
        skin.end_tile(); skin.finish();
        CHECK_THAT(horizontal,WithinAbs(4e-18,1e-29));
        CHECK_THAT(vertical,WithinAbs((same_band ? 6. : 8.)*1e-18,1e-29));
    }
}

TEST_CASE("Exact skin subtraction budget errors retain stage totals", "[SectionThermalSkin]")
{
    const auto first=make_exact_section({ExPolygon(Polygon{{0,0},{2,0},{2,2},{0,2}})}).section;
    const auto second=make_exact_section({ExPolygon(Polygon{{0,0},{2,0},{2,2},{0,2}})}).section;
    bool saw_subtraction=false;
    for (std::size_t budget=64;budget<=8192 && !saw_subtraction;budget+=16) {
        SectionSkinLimits limits; limits.max_work=budget;
        SectionThermalSkin skin(
            [](std::size_t,std::size_t,const ExPolygons&,coord_t,bool) {},
            [](std::size_t,std::size_t,const VerticalSectionPatch&) {},limits,{},{},
            [](std::size_t,std::size_t,const ExactSection&,coord_t,bool) {},
            [](std::size_t,std::size_t,const ExactSectionBoundary&,coord_t,coord_t,const std::array<double,2>&) {});
        try {
            skin.begin_tile({{0,0},{0,0},{2,2},{}});
            skin.add_band({{{0.,1e-9,0,SectionFootprint(first)},0}});
            skin.add_band({{{1e-9,2e-9,0,SectionFootprint(second)},0}});
            skin.end_tile(); skin.finish();
        } catch (const ResourceLimit& error) {
            const std::string message=error.what();
            if (message.find("Exact skin section subtraction failed:")==std::string::npos) continue;
            saw_subtraction=true;
            CHECK(message.find("stage=horizontal face")!=std::string::npos);
            CHECK(message.find("seam subtraction=")!=std::string::npos);
            CHECK(message.find("vertical face=")!=std::string::npos);
            CHECK_THROWS_AS(skin.finish(),std::invalid_argument);
        }
    }
    CHECK(saw_subtraction);
}

TEST_CASE("Exact skin preserves multirow frontier coverage within bounded seam storage", "[SectionThermalSkin]")
{
    const std::array<coord_t,3> start{0,0,0},size{4,6,2};
    SectionIntegration::Moments horizontal{},vertical{};
    SectionSkinLimits limits;
    limits.max_seam_patches=16;
    limits.max_work=100000;
    std::size_t checkpoints=0;
    SectionThermalSkin skin(
        [](std::size_t,std::size_t,const ExPolygons&,coord_t,bool) { FAIL("Native face in exact stream"); },
        [](std::size_t,std::size_t,const VerticalSectionPatch&) { FAIL("Native side in exact stream"); },
        limits,[&] { ++checkpoints; },{},
        [&](std::size_t,std::size_t,const ExactSection& face,coord_t plane,bool) {
            const auto moments=exact_section_horizontal_surface_moments(face,start,size,plane);
            for (std::size_t i=0;i<27;++i) horizontal[i]+=moments.values[i];
        },
        [&](std::size_t,std::size_t,const ExactSectionBoundary& span,coord_t bottom,coord_t top,const std::array<double,2>&) {
            const auto moments=exact_section_vertical_surface_moments(span,start,size,bottom,top,
                ExactBoundaryCellOwnership::AllIncident);
            for (std::size_t i=0;i<27;++i) vertical[i]+=moments.values[i];
        });
    // Lexicographic traversal keeps several next-X-row traces pending while
    // neighboring Y tiles arrive. Two bands require owner-preserving merging.
    for (int x=0;x<2;++x) for (int y=0;y<3;++y) {
        const coord_t left=2*x,bottom=2*y;
        const auto shape=make_exact_section({ExPolygon(Polygon{{left,bottom},{left+2,bottom},
            {left+2,bottom+2},{left,bottom+2}})}).section;
        skin.begin_tile({{x,y},{left,bottom},{left+2,bottom+2},{}});
        for (coord_t z=0;z<2;++z)
            skin.add_band({{{double(z)*1e-9,double(z+1)*1e-9,0,SectionFootprint(shape)},0}});
        skin.end_tile();
    }
    skin.finish();
    CHECK(checkpoints>0);
    const auto whole=make_exact_section({ExPolygon(Polygon{{0,0},{4,0},{4,6},{0,6}})}).section;
    const auto lower=exact_section_horizontal_surface_moments(whole,start,size,0);
    const auto upper=exact_section_horizontal_surface_moments(whole,start,size,2);
    const auto sides=exact_section_vertical_surface_moments(extract_exact_section_boundary(whole).boundary,
        start,size,0,2,ExactBoundaryCellOwnership::AllIncident);
    for (std::size_t i=0;i<27;++i) {
        CHECK_THAT(horizontal[i],WithinAbs(lower.values[i]+upper.values[i],1e-28));
        CHECK_THAT(vertical[i],WithinAbs(sides.values[i],1e-28));
    }
}

TEST_CASE("Skin offline work allowance leaves independent storage guards active", "[SectionThermalSkin]")
{
    const auto horizontal = [](std::size_t, std::size_t, const ExPolygons&, coord_t, bool) {};
    const auto vertical = [](std::size_t, std::size_t, const VerticalSectionPatch&) {};
    SectionSkinLimits limits;
    limits.max_work = 2000000000;
    CHECK_NOTHROW(SectionThermalSkin(horizontal,vertical,limits));
    limits.max_work += 1;
    CHECK_THROWS_AS(SectionThermalSkin(horizontal,vertical,limits),std::invalid_argument);
    limits.max_work = 2000000000;
    limits.max_seam_patches = 250001;
    CHECK_THROWS_AS(SectionThermalSkin(horizontal,vertical,limits),std::invalid_argument);
}

TEST_CASE("Surface line indexing preserves diagonal contacts and native air gaps", "[SectionThermalSkin]")
{
    for (const coord_t gap : {coord_t(0), coord_t(1)}) {
        CAPTURE(gap);
        std::vector<SectionPiece> pieces{
            {{0., .001, 0, {ExPolygon{Polygon{{0,0},{2000000,0},{0,2000000}}}}}, 0},
            {{0., .001, 1, {ExPolygon{Polygon{{2000000+gap,gap},{2000000+gap,2000000+gap},
                                              {gap,2000000+gap},{1000000+gap,1000000+gap}}}}}, 1}};
        std::size_t work = 0;
        const auto patches = partition_vertical_section_patches(pieces,pieces.front(),100,
            [&](std::size_t n) { work += n; });
        double area = 0.;
        for (const auto& patch : patches)
            area += (patch.b-patch.a).cast<double>().norm()*(patch.end-patch.begin)*double(patch.top-patch.bottom)*1e-18;
        CHECK_THAT(area,WithinAbs((4.+(gap ? std::sqrt(8.) : 0.))*1e-6,1e-16));
        CHECK(work < 200);
    }
}

TEST_CASE("Band surface visits share a contact index without changing patches", "[SectionThermalSkin]")
{
    std::vector<SectionPiece> pieces;
    for (coord_t i = 0; i < 24; ++i) {
        const coord_t x = i * 3000000;
        pieces.push_back({{0., .001, 0, {ExPolygon{Polygon{{x,0},{x+2000000,0},
            {x+2000000,2000000},{x,2000000}}}}}, 0});
    }
    std::size_t single_work = 0, batch_work = 0, visited = 0;
    visit_vertical_section_patches(pieces,100,
        [&](const SectionPiece& piece, const std::vector<VerticalSectionPatch>& patches) {
            REQUIRE(&piece == &pieces[visited++]);
            const auto expected = partition_vertical_section_patches(pieces,piece,100,
                [&](std::size_t n) { single_work += n; });
            REQUIRE(patches.size() == expected.size());
            for (std::size_t i = 0; i < patches.size(); ++i) {
                CHECK(patches[i].a == expected[i].a);
                CHECK(patches[i].b == expected[i].b);
                CHECK_THAT(patches[i].begin,WithinAbs(expected[i].begin,1e-15));
                CHECK_THAT(patches[i].end,WithinAbs(expected[i].end,1e-15));
                CHECK(patches[i].bottom == expected[i].bottom);
                CHECK(patches[i].top == expected[i].top);
            }
        }, [&](std::size_t n) { batch_work += n; });
    CHECK(visited == pieces.size());
    CHECK(batch_work * 2 < single_work);
    CHECK_THROWS_AS(visit_vertical_section_patches(pieces,100,
        [](const SectionPiece&, const std::vector<VerticalSectionPatch>&) {},
        [](std::size_t) { throw std::runtime_error("cancel"); }),std::runtime_error);
}

TEST_CASE("Separate islands avoid exact contact indexing without losing exposed area", "[SectionThermalSkin]")
{
    std::vector<SectionPiece> pieces;
    for (coord_t i = 0; i < 8; ++i) {
        const coord_t x = i*2000000;
        pieces.push_back({{0.,.001,0,{ExPolygon{Polygon{{x,0},{x+1000000,0},
            {x+1000000,1000000},{x,1000000}}}}},std::size_t(i)});
    }
    std::size_t work = 0;
    double area = 0.;
    visit_vertical_section_patches(pieces,100,
        [&](const SectionPiece&, const std::vector<VerticalSectionPatch>& patches) {
            for (const auto& patch : patches)
                area += (patch.b-patch.a).cast<double>().norm()*(patch.end-patch.begin)*double(patch.top-patch.bottom)*1e-18;
        }, [&](std::size_t n) { work += n; });
    CHECK_THAT(area,WithinAbs(32e-6,1e-16));
    CHECK(work < 200);
}

TEST_CASE("Horizontal skin retains holes and gaps between otherwise identical bands", "[SectionThermalSkin]")
{
    for (const coord_t gap : {coord_t(0), coord_t(1)}) {
        CAPTURE(gap);
        std::array<double, 2> area{};
        SectionThermalSkin skin(
            [&](std::size_t, std::size_t, const ExPolygons& polygons, coord_t, bool upward) {
                for (const auto& polygon : polygons) area[upward ? 1 : 0] += polygon.area()*1e-18;
            }, [](std::size_t, std::size_t, const VerticalSectionPatch&) {});
        skin.begin_tile({{0,0},{0,0},{10000000,10000000},{}});
        ExPolygon ring{Polygon{{1000000,1000000},{3000000,1000000},{3000000,3000000},{1000000,3000000}}};
        ring.holes.push_back(Polygon{{1500000,1500000},{1500000,2500000},{2500000,2500000},{2500000,1500000}});
        const ExPolygon separate{Polygon{{6000000,1000000},{7000000,1000000},{7000000,2000000},{6000000,2000000}}};
        skin.add_band({{{0.,.001,0,{ring}},0},{{0.,.001,0,{separate}},1}});
        const double bottom = .001+double(gap)*1e-9;
        // Reverse piece order and change material: exact geometric contact,
        // not matching indices or material IDs, controls exposed faces.
        skin.add_band({{{bottom,.002,1,{separate}},1},{{bottom,.002,1,{ring}},0}});
        skin.end_tile();
        skin.finish();
        for (const auto value : area)
            CHECK_THAT(value, WithinAbs((gap ? 8. : 4.)*1e-6,1e-16));
    }
}

TEST_CASE("Skin compacts identical adjacent traces before exhausting its live frontier", "[SectionThermalSkin]")
{
    for (const bool gap : {false,true}) {
        CAPTURE(gap);
        double area = 0.;
        SectionSkinLimits limits;
        limits.max_seam_patches = 8;
        SectionThermalSkin skin([](std::size_t, std::size_t, const ExPolygons&, coord_t, bool) {},
            [&](std::size_t, std::size_t, const VerticalSectionPatch& p) {
                area += (p.b-p.a).cast<double>().norm()*(p.end-p.begin)*double(p.top-p.bottom)*1e-18;
            }, limits);
        skin.begin_tile({{0,0},{0,0},{1000000,1000000},{}});
        const auto add = [&](int i) {
            const double bottom = i*(gap ? .002 : .001);
            skin.add_band({{{bottom,bottom+.001,0,
                {ExPolygon{Polygon{{0,0},{1000000,0},{1000000,1000000},{0,1000000}}}}},0}});
        };
        add(0);
        add(1);
        if (gap) {
            CHECK_THROWS_AS(add(2), ResourceLimit);
        } else {
            for (int i = 2; i < 100; ++i) add(i);
            skin.end_tile();
            skin.finish();
            CHECK_THAT(area, WithinAbs(400e-6,1e-15));
        }
    }
}

TEST_CASE("Seam compaction preserves gaps and ownership across collinear subdivisions", "[SectionThermalSkin]")
{
    for (const unsigned variant : {0u, 1u, 2u, 3u}) {
        CAPTURE(variant);
        std::size_t traces = 0;
        double area = 0.;
        SectionIntegration::Moments moments{};
        SectionThermalSkin skin([](std::size_t, std::size_t, const ExPolygons&, coord_t, bool) {},
            [&](std::size_t material, std::size_t component, const VerticalSectionPatch& p) {
                if (p.a.x() != 1000000 || p.b.x() != 1000000) return;
                ++traces;
                area += (p.b-p.a).cast<double>().norm()*(p.end-p.begin)*double(p.top-p.bottom)*1e-18;
                const auto integrated = vertical_patch_cell_moments(p,{0,0,0},{1000000,1000000,1000000});
                REQUIRE(integrated);
                for (std::size_t i = 0; i < moments.size(); ++i) moments[i] += (*integrated)[i];
                if (variant == 2) CHECK(material == (p.a.y() >= 500000 && p.b.y() >= 500000 ? 1 : 0));
                if (variant == 3) CHECK(component == (p.a.y() >= 500000 && p.b.y() >= 500000 ? 1 : 0));
            });
        skin.begin_tile({{0,0},{0,0},{1000000,1000000},{}});
        const coord_t upper = 500000 + (variant == 1 ? 1 : 0);
        skin.add_band({
            {{0.,.001,0,{ExPolygon{Polygon{{0,0},{1000000,0},{1000000,500000},{0,500000}}}}},0},
            {{0.,.001,variant == 2 ? 1u : 0u,{ExPolygon{Polygon{{0,upper},{1000000,upper},
                {1000000,1000000},{0,1000000}}}}},variant == 3 ? 1u : 0u}});
        skin.end_tile();
        skin.finish();
        CHECK(traces == (variant == 0 ? 1 : 2));
        CHECK_THAT(area,WithinAbs(1e-6-(variant == 1 ? 1e-12 : 0.),1e-18));
        for (unsigned x = 0; x < 3; ++x)
            for (unsigned y = 0; y < 3; ++y)
                for (unsigned z = 0; z < 3; ++z) {
                    const double removed = std::pow(double(upper)/1000000.,y+1)-std::pow(.5,y+1);
                    const double expected = 1e-6*(1.-removed)/double((y+1)*(z+1));
                    CHECK_THAT(moments[x+y*3+z*9],WithinAbs(expected,1e-18));
                }
    }
}

TEST_CASE("A full seam frontier reuses proved provisional representatives", "[SectionThermalSkin]")
{
    for (const bool joined : {false,true}) {
        CAPTURE(joined);
        SectionSkinLimits limits;
        limits.max_seam_patches = 8;
        std::vector<std::size_t> mapping;
        std::size_t queries = 0;
        double area = 0.;
        SectionThermalSkin skin([](std::size_t,std::size_t,const ExPolygons&,coord_t,bool) {},
            [&](std::size_t,std::size_t,const VerticalSectionPatch& p) {
                area += (p.b-p.a).cast<double>().norm()*(p.end-p.begin)*double(p.top-p.bottom)*1e-18;
            },limits,{},[&] { ++queries; return mapping; });
        skin.begin_tile({{0,0},{0,0},{1000000,1000000},{}});
        const auto add = [&](std::size_t i) {
            mapping.push_back(joined ? 0 : i);
            skin.add_band({{{double(i)*.001,double(i+1)*.001,0,
                {ExPolygon{Polygon{{0,0},{1000000,0},{1000000,1000000},{0,1000000}}}}},i}});
        };
        add(0);
        add(1);
        if (!joined) {
            CHECK_THROWS_AS(add(2),ResourceLimit);
            CHECK_THROWS_AS(skin.finish(),std::invalid_argument);
        } else {
            for (std::size_t i = 2; i < 32; ++i) add(i);
            skin.end_tile();
            skin.finish();
            CHECK_THAT(area,WithinAbs(128e-6,1e-15));
        }
        CHECK(queries > 0);
    }
}

TEST_CASE("Invalid provisional seam representatives poison the stream", "[SectionThermalSkin]")
{
    for (const auto& mapping : std::vector<std::vector<std::size_t>>{{1,0},{2},{}}) {
        SectionSkinLimits limits;
        limits.max_seam_patches = 4;
        SectionThermalSkin skin([](std::size_t,std::size_t,const ExPolygons&,coord_t,bool) {},
            [](std::size_t,std::size_t,const VerticalSectionPatch&) {},limits,{},[&] { return mapping; });
        skin.begin_tile({{0,0},{0,0},{1000000,1000000},{}});
        const ExPolygons footprint{ExPolygon{Polygon{{0,0},{1000000,0},{1000000,1000000},{0,1000000}}}};
        skin.add_band({{{0.,.001,0,footprint},0}});
        CHECK_THROWS_AS(skin.add_band({{{.001,.002,0,footprint},0}}),std::invalid_argument);
        CHECK_THROWS_AS(skin.finish(),std::invalid_argument);
    }
}

TEST_CASE("Vertical surface streaming preserves repeated bands and changing owners", "[SectionThermalSkin]")
{
    for (const unsigned variant : {0u,1u,2u,3u}) {
        CAPTURE(variant);
        std::size_t callbacks = 0;
        double area = 0.;
        SectionThermalSkin skin([](std::size_t, std::size_t, const ExPolygons&, coord_t, bool) {},
            [&](std::size_t, std::size_t, const VerticalSectionPatch& p) {
                ++callbacks;
                area += (p.b-p.a).cast<double>().norm()*(p.end-p.begin)*double(p.top-p.bottom)*1e-18;
            });
        skin.begin_tile({{0,0},{0,0},{3000000,3000000},{}});
        for (int i = 0; i < 100; ++i) {
            const double bottom = i*(variant == 1 ? .002 : .001);
            const std::size_t material = variant == 2 ? std::size_t(i%2) : 0;
            const std::size_t component = variant == 3 ? std::size_t(i%2) : 0;
            skin.add_band({{{bottom,bottom+.001,material,
                {ExPolygon{Polygon{{1000000,1000000},{2000000,1000000},{2000000,2000000},{1000000,2000000}}}}},component}});
        }
        skin.end_tile();
        skin.finish();
        CHECK_THAT(area,WithinAbs(400e-6,1e-15));
        CHECK(callbacks == 400);
    }
}

TEST_CASE("Per-piece patch limits preserve faces from separate islands", "[SectionThermalSkin]")
{
    double area = 0.;
    SectionSkinLimits limits;
    limits.max_seam_patches = 4;
    SectionThermalSkin skin([](std::size_t, std::size_t, const ExPolygons&, coord_t, bool) {},
        [&](std::size_t, std::size_t, const VerticalSectionPatch& p) {
            area += (p.b-p.a).cast<double>().norm()*(p.end-p.begin)*double(p.top-p.bottom)*1e-18;
        }, limits);
    skin.begin_tile({{0,0},{0,0},{10000000,10000000},{}});
    for (int i = 0; i < 2; ++i) {
        std::vector<SectionPiece> pieces;
        for (const coord_t x : {coord_t(1000000),coord_t(4000000)})
            pieces.push_back({{i*.001,(i+1)*.001,0,
                {ExPolygon{Polygon{{x,1000000},{x+1000000,1000000},{x+1000000,2000000},{x,2000000}}}}},0});
        skin.add_band(std::move(pieces));
    }
    skin.end_tile();
    skin.finish();
    CHECK_THAT(area,WithinAbs(16e-6,1e-16));
}

TEST_CASE("Streamed skin removes tile cuts and partial height material contacts", "[SectionThermalSkin]")
{
    const std::vector<DepositedSegment> beads{{{.0005, .001, .0005}, {.0035, .001, .0005}, .002, .0005, 0, BeadSection::RectangularEnvelope},
                                              {{.0015, .001, .001}, {.0025, .001, .001}, .001, .0005, 1, BeadSection::RectangularEnvelope}};
    const auto geometry = local_deposition_sections(beads, {0., 0.}, {.001, .001}, 0., .0001);
    const auto topology = qualify_local_section_topology(geometry.tiles, 0.);
    const std::vector<SectionThermalMaterial> materials{{Eigen::Matrix3d::Identity(), 1e6, 300.}, {Eigen::Matrix3d::Identity(), 1e6, 300.}};
    const auto reference = assemble_section_thermal_volume(topology, materials, {0., 0., 0.}, {.001, .001, .0005});
    // Per material: upward, downward (including bed), vertical surface area.
    std::array<std::array<double, 3>, 2> area{}, expected{};
    for (const auto& face : reference.horizontal_faces)
        expected[face.material][face.kind == SectionHorizontalFaceKind::Upward ? 0 : 1] += face.area_moments[0];
    for (const auto& face : reference.vertical_faces)
        expected[face.material][2] += face.area_moments[0];
    SectionThermalSkin skin(
        [&](std::size_t m, std::size_t, const ExPolygons& p, coord_t, bool up) {
            for (const auto& polygon : p)
                area[m][up ? 0 : 1] += polygon.area() * 1e-18;
        },
        [&](std::size_t m, std::size_t, const VerticalSectionPatch& p) {
            const auto delta = (p.b - p.a).cast<double>();
            area[m][2] += delta.norm() * (p.end - p.begin) * double(p.top - p.bottom) * 1e-18;
            CHECK_THAT(p.outward_normal.norm(), WithinAbs(1., 1e-12));
        });
    std::vector<SectionPiece> band;
    LocalTopologyBuilder builder(
        0., {}, [&](const SectionPiece& piece) { band.push_back(piece); }, {}, StreamTopologyLimits{}, nullptr,
        LocalTileOrder::Lexicographic);
    LocalSectionPacketVisitor visitor;
    visitor.begin_tile = [&](const LocalSectionTile& h) {
        builder.begin_tile(h);
        skin.begin_tile(h);
    };
    visitor.band = [&](const LocalSectionTile& p) {
        builder.add_band(p);
        skin.add_band(std::move(band));
        band.clear();
    };
    visitor.end_tile = [&] {
        builder.end_tile();
        skin.end_tile();
    };
    visit_local_deposition_packets(beads, {0., 0.}, {.001, .001}, 0., .0001, visitor);
    const auto certificate = builder.finish();
    skin.finish();
    CHECK(certificate.components() == 1);
    for (std::size_t m = 0; m < 2; ++m)
        for (std::size_t f = 0; f < 3; ++f)
            CHECK_THAT(area[m][f], WithinAbs(expected[m][f], 1e-16));
    CHECK_THROWS_AS(skin.finish(), std::invalid_argument);
}

TEST_CASE("Streamed skin exposes missing neighbors and poisons rejected streams", "[SectionThermalSkin]")
{
    double sides          = 0.;
    const auto horizontal = [](std::size_t, std::size_t, const ExPolygons&, coord_t, bool) {};
    const auto vertical   = [&](std::size_t, std::size_t, const VerticalSectionPatch& p) {
        sides += (p.b - p.a).cast<double>().norm() * (p.end - p.begin) * double(p.top - p.bottom) * 1e-18;
    };
    SectionThermalSkin skin(horizontal, vertical);
    const auto slab = [](coord_t x) {
        return SectionPiece{{0., .001, 0, {ExPolygon{Polygon{{x, 0}, {x + 1000000, 0}, {x + 1000000, 1000000}, {x, 1000000}}}}}, 0};
    };
    skin.begin_tile({{0, 0}, {0, 0}, {1000000, 1000000}, {}});
    skin.add_band({slab(0)});
    skin.end_tile();
    skin.begin_tile({{2, 0}, {2000000, 0}, {3000000, 1000000}, {}});
    skin.add_band({slab(2000000)});
    skin.end_tile();
    skin.finish();
    CHECK_THAT(sides, WithinAbs(8e-6, 1e-16));
    SectionSkinLimits limits;
    limits.max_live_points = 3;
    SectionThermalSkin bounded(horizontal, vertical, limits);
    bounded.begin_tile({{0, 0}, {0, 0}, {1000000, 1000000}, {}});
    CHECK_THROWS_AS(bounded.add_band({slab(0)}), ResourceLimit);
    CHECK_THROWS_AS(bounded.end_tile(), std::invalid_argument);
    SectionThermalSkin canceled(horizontal, vertical, {}, [] { throw std::runtime_error("cancel"); });
    CHECK_THROWS_AS(canceled.begin_tile({{0, 0}, {0, 0}, {1000000, 1000000}, {}}), std::runtime_error);
    CHECK_THROWS_AS(canceled.finish(), std::invalid_argument);
}

TEST_CASE("Streamed skin retains hole walls and exposed undersides across empty bands", "[SectionThermalSkin]")
{
    ExPolygon ring{Polygon{{0, 0}, {4000000, 0}, {4000000, 4000000}, {0, 4000000}}};
    ring.holes.push_back(Polygon{{1000000, 1000000}, {1000000, 3000000}, {3000000, 3000000}, {3000000, 1000000}});
    double upward = 0., downward = 0., vertical_area = 0.;
    Eigen::Vector2d normal_integral = Eigen::Vector2d::Zero();
    SectionThermalSkin skin(
        [&](std::size_t, std::size_t, const ExPolygons& p, coord_t, bool up) {
            for (const auto& polygon : p)
                (up ? upward : downward) += polygon.area() * 1e-18;
        },
        [&](std::size_t, std::size_t, const VerticalSectionPatch& p) {
            const double area = (p.b - p.a).cast<double>().norm() * (p.end - p.begin) * double(p.top - p.bottom) * 1e-18;
            vertical_area += area;
            normal_integral += area * p.outward_normal;
        });
    skin.begin_tile({{0, 0}, {0, 0}, {4000000, 4000000}, {}});
    skin.add_band({{{0., .001, 0, {ring}}, 0}});
    // The second slab is intentionally floating. Skin extraction must expose
    // both sides of the intervening air gap, never bridge that empty interval.
    skin.add_band({{{.002, .003, 0, {ring}}, 1}});
    skin.end_tile();
    skin.finish();
    CHECK_THAT(upward, WithinAbs(24e-6, 1e-16));
    CHECK_THAT(downward, WithinAbs(24e-6, 1e-16));
    CHECK_THAT(vertical_area, WithinAbs(48e-6, 1e-16));
    CHECK_THAT(normal_integral.norm(), WithinAbs(0., 1e-16));
}

TEST_CASE("Streamed skin removes only the actual rectangle shared with a neighboring tile", "[SectionThermalSkin]")
{
    double right = 0., left = 0.;
    SectionThermalSkin skin([](std::size_t, std::size_t, const ExPolygons&, coord_t, bool) {},
                            [&](std::size_t m, std::size_t, const VerticalSectionPatch& p) {
                                if (p.a.x() != 1000000 || p.b.x() != 1000000)
                                    return;
                                const double area = (p.b - p.a).cast<double>().norm() * (p.end - p.begin) * double(p.top - p.bottom) *
                                                    1e-18;
                                if (m == 0) {
                                    right += area;
                                    CHECK_THAT(p.outward_normal.x(), WithinAbs(1., 1e-15));
                                } else {
                                    left += area;
                                    CHECK_THAT(p.outward_normal.x(), WithinAbs(-1., 1e-15));
                                }
                            });
    skin.begin_tile({{0, 0}, {0, 0}, {1000000, 4000000}, {}});
    skin.add_band({{{0., .004, 0, {ExPolygon{Polygon{{0, 0}, {1000000, 0}, {1000000, 4000000}, {0, 4000000}}}}}, 0}});
    skin.end_tile();
    skin.begin_tile({{1, 0}, {1000000, 0}, {2000000, 4000000}, {}});
    skin.add_band(
        {{{.001, .003, 1, {ExPolygon{Polygon{{1000000, 1000000}, {2000000, 1000000}, {2000000, 3000000}, {1000000, 3000000}}}}}, 0}});
    skin.end_tile();
    skin.finish();
    CHECK_THAT(right, WithinAbs(12e-6, 1e-16));
    CHECK_THAT(left, WithinAbs(0., 1e-16));
}

TEST_CASE("Seam height indexing preserves partial contacts across separated bands", "[SectionThermalSkin]")
{
    std::array<double, 2> area{};
    SectionThermalSkin skin([](std::size_t, std::size_t, const ExPolygons&, coord_t, bool) {},
        [&](std::size_t material, std::size_t, const VerticalSectionPatch& p) {
            if (p.a.x() == 1000000 && p.b.x() == 1000000)
                area[material] += (p.b-p.a).cast<double>().norm() * (p.end-p.begin) *
                    double(p.top-p.bottom) * 1e-18;
        });
    constexpr int bands = 100;
    for (int side = 0; side < 2; ++side) {
        const coord_t x = side * 1000000;
        skin.begin_tile({{side,0},{x,0},{x+1000000,1000000},{}});
        for (int band = 0; band < bands; ++band) {
            const double bottom = band * .002 + side * .0005;
            skin.add_band({{{bottom,bottom+.001,std::size_t(side),
                {ExPolygon{Polygon{{x,0},{x+1000000,0},{x+1000000,1000000},{x,1000000}}}}},0}});
        }
        skin.end_tile();
    }
    skin.finish();
    for (const double exposed : area)
        CHECK_THAT(exposed,WithinAbs(bands * .5e-6,1e-16));
}

TEST_CASE("Spatial seam indexing subtracts only overlapping edge spans", "[SectionThermalSkin]")
{
    std::array<double,2> area{};
    SectionThermalSkin skin([](std::size_t, std::size_t, const ExPolygons&, coord_t, bool) {},
        [&](std::size_t material, std::size_t, const VerticalSectionPatch& p) {
            if (p.a.x() == 1000000 && p.b.x() == 1000000)
                area[material] += (p.b-p.a).cast<double>().norm()*(p.end-p.begin)*double(p.top-p.bottom)*1e-18;
        });
    constexpr int strips = 100;
    for (int side = 0; side < 2; ++side) {
        const coord_t x = side*1000000;
        skin.begin_tile({{side,0},{x,0},{x+1000000,strips*30000},{}});
        std::vector<SectionPiece> pieces;
        for (int i = 0; i < strips; ++i) {
            const coord_t y = i*30000+side*10000;
            pieces.push_back({{0.,.001,std::size_t(side),
                {ExPolygon{Polygon{{x,y},{x+1000000,y},{x+1000000,y+20000},{x,y+20000}}}}},std::size_t(i)});
        }
        skin.add_band(std::move(pieces));
        skin.end_tile();
    }
    skin.finish();
    for (const double exposed : area)
        CHECK_THAT(exposed,WithinAbs(strips*1e-8,1e-16));
}

TEST_CASE("Tile traces preserve complete and partial contact at the live frontier limit", "[SectionThermalSkin]")
{
    for (const coord_t gap : {coord_t(0),coord_t(1)}) {
        CAPTURE(gap);
        double exposed = 0.;
        SectionSkinLimits limits;
        limits.max_seam_patches = 4;
        SectionThermalSkin skin([](std::size_t,std::size_t,const ExPolygons&,coord_t,bool) {},
            [&](std::size_t,std::size_t,const VerticalSectionPatch& patch) {
                if (patch.a.x() == 1000000 && patch.b.x() == 1000000)
                    exposed += (patch.b-patch.a).cast<double>().norm()*(patch.end-patch.begin)*
                        double(patch.top-patch.bottom)*1e-18;
            },limits);
        for (int side = 0; side < 2; ++side) {
            const coord_t x = side*1000000;
            skin.begin_tile({{side,0},{x,0},{x+1000000,1000000},{}});
            std::vector<SectionPiece> pieces;
            for (coord_t strip = 0; strip < 2; ++strip) {
                const coord_t lo = side ? x : x+100000, hi = side ? x+900000 : x+1000000;
                const coord_t y = 100000+strip*400000+side*gap;
                pieces.push_back({{0.,.001,std::size_t(side),{ExPolygon{Polygon{
                    {lo,y},{hi,y},{hi,y+200000},{lo,y+200000}}}}},std::size_t(strip)});
            }
            skin.add_band(std::move(pieces));
            skin.end_tile();
        }
        skin.finish();
        // Two strips, each exposing one one-nm remnant on either side.
        CHECK_THAT(exposed,WithinAbs(double(gap)*4e-12,1e-20));
    }
}

TEST_CASE("Fragment streaming subtracts many contacts without retaining all exposed remnants", "[SectionThermalSkin]")
{
    constexpr coord_t strips = 100, pitch = 30000, width = 10000, offset = 100000;
    std::array<double,2> area{};
    SectionSkinLimits limits;
    limits.max_seam_patches = 104;
    SectionThermalSkin skin([](std::size_t,std::size_t,const ExPolygons&,coord_t,bool) {},
        [&](std::size_t material,std::size_t,const VerticalSectionPatch& patch) {
            if (patch.a.x() != 1000000 || patch.b.x() != 1000000) return;
            area[material] += (patch.b-patch.a).cast<double>().norm()*(patch.end-patch.begin)*
                double(patch.top-patch.bottom)*1e-18;
        },limits);
    skin.begin_tile({{0,0},{0,0},{1000000,offset*2+strips*pitch},{}});
    skin.add_band({{{0.,.001,0,{ExPolygon{Polygon{{100000,offset},{1000000,offset},
        {1000000,offset+strips*pitch},{100000,offset+strips*pitch}}}}},0}});
    skin.end_tile();
    skin.begin_tile({{1,0},{1000000,0},{2000000,offset*2+strips*pitch},{}});
    std::vector<SectionPiece> pieces;
    for (coord_t i = 0; i < strips; ++i) {
        const coord_t y = offset+i*pitch;
        pieces.push_back({{0.,.001,1,{ExPolygon{Polygon{{1000000,y},{1900000,y},
            {1900000,y+width},{1000000,y+width}}}}},std::size_t(i)});
    }
    skin.add_band(std::move(pieces));
    skin.end_tile();
    skin.finish();
    CHECK_THAT(area[0],WithinAbs(double(strips*(pitch-width))*1e-12,1e-16));
    CHECK_THAT(area[1],WithinAbs(0.,1e-16));
}

TEST_CASE("Fragment stack exhaustion invalidates the provisional surface stream", "[SectionThermalSkin]")
{
    SectionSkinLimits limits;
    limits.max_seam_patches = 4;
    SectionThermalSkin skin([](std::size_t,std::size_t,const ExPolygons&,coord_t,bool) {},
        [](std::size_t,std::size_t,const VerticalSectionPatch&) {},limits);
    skin.begin_tile({{0,0},{0,0},{1000000,2000000},{}});
    skin.add_band({{{0.,.001,0,{ExPolygon{Polygon{{100000,100000},{1000000,100000},
        {1000000,1100000},{100000,1100000}}}}},0}});
    skin.end_tile();
    skin.begin_tile({{1,0},{1000000,0},{2000000,2000000},{}});
    std::vector<SectionPiece> pieces;
    for (coord_t i = 0; i < 3; ++i) {
        const coord_t y = 300000+i*300000;
        pieces.push_back({{.00025,.00075,1,{ExPolygon{Polygon{{1000000,y},{1900000,y},
            {1900000,y+100000},{1000000,y+100000}}}}},std::size_t(i)});
    }
    skin.add_band(std::move(pieces));
    CHECK_THROWS_AS(skin.end_tile(),ResourceLimit);
    CHECK_THROWS_AS(skin.finish(),std::invalid_argument);
}
