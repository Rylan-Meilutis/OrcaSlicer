#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "libslic3r/ThermoelasticBoundary.hpp"
#include "libslic3r/ThermoelasticExactSections.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/libslic3r.h"
#include <stdexcept>
#include <array>
#include <cstddef>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;

namespace {
ExPolygon rectangle(coord_t x, coord_t y, coord_t width, coord_t height)
{
    return ExPolygon(Polygon({Point(x, y), Point(x + width, y),
        Point(x + width, y + height), Point(x, y + height)}));
}
}

TEST_CASE("Boundary representation changes require an explicit narrow band", "[ThermoelasticBoundary]")
{
    const ExPolygons before{rectangle(200000000, 200000000, 1000, 1000)};
    const ExPolygons after{rectangle(200000000, 200000000, 1001, 1000)};
    CHECK_THROWS_AS(qualify_section_boundary_error(before, after, 0, 1000, 100000), std::invalid_argument);
    const auto error = qualify_section_boundary_error(before, after, 2, 1000, 100000);
    CHECK(error.uncertainty_nm == 5);
    CHECK(error.added_area_bound_m2 > 0.);
    CHECK(error.removed_area_bound_m2 == 0.);
    CHECK_THROWS_AS(qualify_section_boundary_error(before,
        {rectangle(200000000, 200000000, 1010, 1000)}, 2, 1000, 100000), std::invalid_argument);
}

TEST_CASE("Boundary uncertainty preserves individual islands and holes", "[ThermoelasticBoundary]")
{
    const ExPolygons before{rectangle(0, 0, 1000, 1000), rectangle(1001, 0, 1, 1)};
    CHECK_THROWS_AS(qualify_section_boundary_error(before, {rectangle(0, 0, 1001, 1000)},
        2, 1000, 100000), std::invalid_argument);
    auto holed = rectangle(0, 0, 1000, 1000);
    auto hole = rectangle(500, 500, 2, 2).contour;
    hole.reverse();
    holed.holes.push_back(hole);
    auto expanded = holed;
    expanded.contour = rectangle(0, 0, 1001, 1000).contour;
    CHECK_THROWS_AS(qualify_section_boundary_error({holed}, {expanded},
        2, 1000, 100000), std::invalid_argument);
}

TEST_CASE("Boundary uncertainty rejects invalid orientation and bounded work exhaustion", "[ThermoelasticBoundary]")
{
    auto polygon = rectangle(200000000, 200000000, 2, 1);
    polygon.contour.reverse();
    CHECK_THROWS_AS(qualify_section_boundary_error({polygon}, {polygon}, 2, 1000, 100000), std::invalid_argument);
    const ExPolygons before{rectangle(0, 0, 1000, 1000)};
    const ExPolygons after{rectangle(0, 0, 1001, 1000)};
    CHECK_THROWS_AS(qualify_section_boundary_error(before, after, 2, 1000, 1), ResourceLimit);
    CHECK_THROWS_AS(qualify_section_boundary_error(before, after, 2, 1000, 100000,
        [] { throw std::runtime_error("cancel"); }), std::runtime_error);
}

TEST_CASE("Exact boundary identity adds no uncertainty to an existing small hole", "[ThermoelasticBoundary]")
{
    auto original=rectangle(0,0,1000,1000);
    auto hole=rectangle(500,500,2,2).contour;
    hole.reverse();
    original.holes.push_back(hole);
    auto represented=original;
    represented.contour=Polygon{{0,0},{500,0},{1000,0},{1000,1000},{0,1000}};
    const auto identity=qualify_section_boundary_error({original},{represented},2,1000,100000);
    CHECK(identity.uncertainty_nm==0);
    CHECK(identity.added_area_bound_m2==0.);
    CHECK(identity.removed_area_bound_m2==0.);
    CHECK(identity.work>0);
    CHECK_THROWS_AS(qualify_section_boundary_error({original},{represented},2,1000,identity.work-1),ResourceLimit);
    CHECK_THROWS_AS(qualify_section_boundary_error({original},{represented},2,1000,100000,
        [] { throw std::runtime_error("cancel identity"); }),std::runtime_error);
    // Actual exterior movement still requires the original hole to survive
    // the declared uncertainty; exact identity is not a morphology bypass.
    represented.contour=rectangle(0,0,1001,1000).contour;
    CHECK_THROWS_AS(qualify_section_boundary_error({original},{represented},2,1000,100000),std::invalid_argument);
    const ExPolygons tiny{rectangle(0,0,1,1)};
    CHECK_THROWS_AS(qualify_section_boundary_error(tiny,tiny,2,1000,100000),std::invalid_argument);
}

TEST_CASE("Boundary identity rejects exact changes hidden by rounded unions", "[ThermoelasticBoundary]")
{
    // The lip closes this outlet at a rational crossing less than one native
    // unit from the outer edge. Integer union opens a portal into the hole.
    const ExPolygon open_box(Polygon{{0,0},{1000,0},{1000,500},{700,500},{700,200},{200,200},
        {200,800},{700,800},{700,501},{1000,501},{1000,1000},{0,1000}});
    const ExPolygon lip(Polygon{{800,699},{1001,499},{1001,1000},{800,1000}});
    const ExPolygons original{open_box,lip};
    const auto represented=union_ex(original);
    REQUIRE(diff_ex(union_ex(original),union_ex(represented)).empty());
    REQUIRE(diff_ex(union_ex(represented),union_ex(original)).empty());
    CHECK_THROWS_AS(qualify_section_boundary_error(original,represented,2,1000,100000),std::invalid_argument);
}

TEST_CASE("Boundary certificates include rational area omitted by rounded differences", "[ThermoelasticBoundary]")
{
    const bool visible_addition=GENERATE(false,true);
    const coord_t extent=GENERATE(coord_t(1000),coord_t(1000000));
    // The exact union retains a .2 nm^2 triangle at x=0..1 below y=.4x.
    // Its rational corner (1,.4) disappears from the integer union. The large
    // remaining solid retains its component and hole topology under offsets.
    // The large variant also spans source-edge capsule junctions along a long
    // thin exterior ribbon, requiring deeper bounded exact subdivision.
    const ExPolygons original{
        ExPolygon(Polygon{{0,0},{extent,0},{extent,2*extent/5}}),
        ExPolygon(Polygon{{1,0},{extent,2*extent/5},{extent,extent},{1,extent}})};
    auto represented=union_ex(original);
    REQUIRE(diff_ex(union_ex(original),represented).empty());
    if (visible_addition)
        for (auto& polygon : represented)
            for (auto& point : polygon.contour.points)
                if (point.x()==extent) point.x()+=1;
    REQUIRE(diff_ex(union_ex(original),union_ex(represented)).empty());
    const auto error=qualify_section_boundary_error(original,represented,2,1000,4000000);
    CHECK(error.uncertainty_nm==5);
    CHECK(error.removed_area_bound_m2>=2e-19);
    CHECK(error.removed_area_bound_m2<1e-18);
    if (visible_addition) {
        REQUIRE_FALSE(diff_ex(union_ex(represented),union_ex(original)).empty());
        CHECK(error.added_area_bound_m2>0.);
    } else {
        CHECK(error.added_area_bound_m2==0.);
    }
    const auto exact=certify_exact_boundary_difference(original,represented,5);
    CHECK(exact.removed_area_bound_m2==error.removed_area_bound_m2);
    CHECK_THROWS_AS(certify_exact_boundary_difference(original,represented,5,exact.work-1),ResourceLimit);
    CHECK_THROWS_AS(certify_exact_boundary_difference(original,represented,5,4000000,
        [] { throw std::runtime_error("cancel exact band"); }),std::runtime_error);
}

TEST_CASE("Exact hole diagnostics distinguish remote changes from changed hole boundaries", "[ThermoelasticBoundary]")
{
    auto original=rectangle(0,0,1000,1000);
    auto hole=rectangle(400,400,100,100).contour;
    hole.reverse();
    original.holes.push_back(hole);
    auto represented=original;
    represented.contour=rectangle(0,0,1001,1000).contour;
    const auto remote=diagnose_exact_hole_change({original},{represented},hole);
    REQUIRE(remote.source_hole_unique);
    REQUIRE(remote.represented_hole_unique);
    REQUIRE(remote.correspondence_unique);
    REQUIRE(remote.hole_unchanged.has_value());
    CHECK(*remote.hole_unchanged);
    REQUIRE(remote.difference_empty.has_value());
    CHECK_FALSE(*remote.difference_empty);
    REQUIRE(remote.change_distance_squared_nm2.has_value());
    CHECK((*remote.change_distance_squared_nm2)[0]<=250000.);
    CHECK((*remote.change_distance_squared_nm2)[1]>=250000.);
    CHECK((*remote.change_distance_squared_nm2)[0]>0.);
    CHECK_THROWS_AS(diagnose_exact_hole_change({original},{represented},hole,remote.work-1),ResourceLimit);
    std::size_t checkpoints=0;
    CHECK_THROWS_WITH(diagnose_exact_hole_change({original},{represented},hole,4000000,[&] {
        if (++checkpoints==10) throw std::invalid_argument("hole diagnostic cancellation");
    }),"hole diagnostic cancellation");

    represented=original;
    represented.holes.front()=rectangle(400,400,101,100).contour;
    represented.holes.front().reverse();
    const auto local=diagnose_exact_hole_change({original},{represented},hole);
    REQUIRE(local.correspondence_unique);
    REQUIRE(local.hole_unchanged.has_value());
    CHECK_FALSE(*local.hole_unchanged);
    REQUIRE(local.change_distance_squared_nm2.has_value());
    CHECK((*local.change_distance_squared_nm2)[0]<=0.);
    CHECK((*local.change_distance_squared_nm2)[1]<=0.);
    const auto identity=diagnose_exact_hole_change({original},{original},hole);
    REQUIRE(identity.hole_unchanged.has_value());
    CHECK(*identity.hole_unchanged);
    REQUIRE(identity.difference_empty.has_value());
    CHECK(*identity.difference_empty);
    CHECK_FALSE(identity.change_distance_squared_nm2.has_value());
}

TEST_CASE("Exact hole diagnostics leave missing or ambiguous correspondence unknown", "[ThermoelasticBoundary]")
{
    auto original=rectangle(0,0,1000,1000);
    for (const coord_t x : {coord_t(200),coord_t(600)}) {
        auto hole=rectangle(x,400,100,100).contour;
        hole.reverse();
        original.holes.push_back(hole);
    }
    auto hint=rectangle(100,300,700,300).contour;
    hint.reverse();
    const auto ambiguous=diagnose_exact_hole_change({original},{original},hint);
    CHECK_FALSE(ambiguous.source_hole_unique);
    CHECK_FALSE(ambiguous.represented_hole_unique);
    CHECK_FALSE(ambiguous.correspondence_unique);
    CHECK_FALSE(ambiguous.hole_unchanged.has_value());
    CHECK_FALSE(ambiguous.change_distance_squared_nm2.has_value());
    auto represented=original;
    represented.holes.clear();
    const auto missing=diagnose_exact_hole_change({original},{represented},original.holes.front());
    CHECK(missing.source_hole_unique);
    CHECK_FALSE(missing.represented_hole_unique);
    CHECK_FALSE(missing.hole_unchanged.has_value());
    CHECK_FALSE(missing.change_distance_squared_nm2.has_value());
}

TEST_CASE("Boundary uncertainty rejects unresolved necks and invented contacts", "[ThermoelasticBoundary]")
{
    const ExPolygons lobes = union_ex(ExPolygons{rectangle(0, 0, 100, 100),
        rectangle(200, 0, 100, 100), rectangle(99, 49, 102, 2)});
    const ExPolygons modified = union_ex(lobes, ExPolygons{rectangle(0, 0, 101, 100)});
    CHECK_THROWS_AS(qualify_section_boundary_error(lobes, modified, 2, 1000, 100000), std::invalid_argument);
    const ExPolygons separate{rectangle(0, 0, 100, 100), rectangle(102, 0, 100, 100)};
    CHECK_THROWS_AS(qualify_section_boundary_error(separate, {rectangle(0, 0, 202, 100)},
        2, 1000, 100000), std::invalid_argument);
    // Equal component totals cannot substitute a nearby vanished island.
    const ExPolygons before{rectangle(0, 0, 100, 100), rectangle(101, 0, 1, 1)};
    const ExPolygons after{rectangle(0, 0, 100, 100), rectangle(101, 2, 1, 1)};
    CHECK_THROWS_AS(qualify_section_boundary_error(before, after, 2, 1000, 100000), std::invalid_argument);
}

TEST_CASE("Boundary bands do not use small changed area as containment evidence", "[ThermoelasticBoundary]")
{
    const ExPolygons before{rectangle(0, 0, 1000, 1000)};
    const ExPolygons thin{ExPolygon(Polygon({Point(0, 0), Point(1000, 0),
        Point(1003, 999), Point(1000, 1000), Point(0, 1000)}))};
    CHECK(qualify_section_boundary_error(before, thin, 4, 1000, 100000).added_area_bound_m2 > 0.);
    auto far = before;
    far.emplace_back(Polygon({Point(1010, 0), Point(1011, 0), Point(1010, 1)}));
    CHECK_THROWS_AS(qualify_section_boundary_error(before, far, 2, 1000, 100000), std::invalid_argument);
    // A triangle with vertices at different boundary edges has an interior
    // far from the boundary. Nearest-edge tests of vertices alone are unsafe.
    auto holed = rectangle(0, 0, 1000, 1000);
    auto hole = Polygon({Point(1, 1), Point(999, 1), Point(1, 999)});
    hole.reverse();
    holed.holes.push_back(hole);
    CHECK_THROWS_AS(qualify_section_boundary_error(before, {holed}, 2, 1000, 100000), std::invalid_argument);
}

TEST_CASE("Boundary coverage follows a narrow bent strip without bridging its interior", "[ThermoelasticBoundary]")
{
    const coord_t x = 200000000, y = 200000000;
    const ExPolygons before{ExPolygon(Polygon({Point(x, y), Point(x + 1000, y),
        Point(x + 1000, y + 1000), Point(x + 500, y + 998), Point(x, y + 1000)}))};
    auto after = before;
    for (auto& p : after.front().contour.points) if (p.y() > y) p.y() += 3;
    const auto result = qualify_section_boundary_error(before, after, 4, 1000, 100000);
    CHECK(result.uncertainty_nm == 7);
    CHECK(result.added_area_bound_m2 > 0.);
}

TEST_CASE("Directed boundary chains retain collinear stitches but reject missing occupied triangles", "[ThermoelasticBoundary]")
{
    const ExPolygon polygon(Polygon({Point(0, 0), Point(50, 0), Point(100, 0), Point(100, 100), Point(0, 100)}));
    std::vector<std::array<std::size_t, 3>> triangles{{0, 2, 3}, {0, 3, 4}, {0, 1, 2}};
    CHECK(Detail::certify_boundary_triangulation(polygon, triangles, 10000) > 0);
    auto missing = triangles;
    missing.erase(missing.begin());
    CHECK_THROWS_AS(Detail::certify_boundary_triangulation(polygon, missing, 10000), std::invalid_argument);
    auto duplicate = triangles;
    duplicate.push_back(triangles.front());
    CHECK_THROWS_AS(Detail::certify_boundary_triangulation(polygon, duplicate, 10000), std::invalid_argument);
    CHECK_THROWS_AS(Detail::certify_boundary_triangulation(polygon, triangles, 1), ResourceLimit);
    CHECK_THROWS_AS(Detail::certify_boundary_triangulation(polygon, triangles, 10000,
        [] { throw std::runtime_error("cancel"); }), std::runtime_error);
}

TEST_CASE("Directed triangle coverage preserves hole boundaries", "[ThermoelasticBoundary]")
{
    auto polygon = rectangle(0, 0, 100, 100);
    polygon.holes.emplace_back(Points{Point(25, 25), Point(25, 75), Point(75, 75), Point(75, 25)});
    std::vector<std::array<std::size_t, 3>> triangles{{0, 1, 7}, {0, 7, 4}, {1, 2, 6}, {1, 6, 7},
        {2, 3, 5}, {2, 5, 6}, {3, 0, 4}, {3, 4, 5}};
    CHECK(Detail::certify_boundary_triangulation(polygon, triangles, 10000) > 0);
    triangles.push_back({4, 7, 6});
    triangles.push_back({4, 6, 5});
    CHECK_THROWS_AS(Detail::certify_boundary_triangulation(polygon, triangles, 10000), std::invalid_argument);
}
