#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Support/TypedSupportPaint.hpp"

using namespace Slic3r;

namespace {
Polygon rectangle(int x0, int y0, int x1, int y1)
{
    return Polygon{Points{Point(x0, y0), Point(x1, y0), Point(x1, y1), Point(x0, y1)}};
}
double contact_area(const ExPolygons &polygons)
{
    double result = 0;
    for (const ExPolygon &polygon : polygons) result += polygon.area();
    return result;
}
}

TEST_CASE("Typed support contacts partition demand with blocker precedence", "[SupportPaint]")
{
    const ExPolygons demand{ExPolygon(rectangle(0, 0, 10000, 10000))};
    const auto result = partition_painted_support_contacts(demand,
        {rectangle(0, 0, 6000, 10000)}, {rectangle(4000, 0, 8000, 10000)},
        {rectangle(0, 0, 1000, 10000)});
    CHECK_THAT(contact_area(result.normal), Catch::Matchers::WithinAbs(50000000., 1.));
    CHECK_THAT(contact_area(result.tree), Catch::Matchers::WithinAbs(20000000., 1.));
    CHECK_THAT(contact_area(result.automatic), Catch::Matchers::WithinAbs(20000000., 1.));
    CHECK(intersection_ex(result.normal, result.tree).empty());
    CHECK(intersection_ex(result.normal, result.automatic).empty());
    CHECK(intersection_ex(result.tree, result.automatic).empty());
}

TEST_CASE("Small painted normal contacts are retained without expanding demand", "[SupportPaint]")
{
    const ExPolygons demand{ExPolygon(rectangle(0, 0, 10000, 10000))};
    const auto result = partition_painted_support_contacts(demand,
        {rectangle(-1000, -1000, 1000, 1000)}, {}, {});
    REQUIRE(result.normal.size() == 1);
    CHECK_THAT(contact_area(result.normal), Catch::Matchers::WithinAbs(1000000., 1.));
    CHECK(diff_ex(result.normal, demand).empty());
    CHECK_THAT(contact_area(result.normal) + contact_area(result.automatic),
               Catch::Matchers::WithinAbs(contact_area(demand), 1.));
}
