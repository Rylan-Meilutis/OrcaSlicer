#include <catch2/catch_test_macros.hpp>
#include "libslic3r/GCode/LocalExtrusionGroups.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Point.hpp"
#include <cstddef>
#include <vector>

using namespace Slic3r;

TEST_CASE("Local extrusion groups keep walls and enclosed fill together", "[LocalExtrusionGroups]")
{
    const std::vector<BoundingBox> boxes = {
        {Point(0, 0), Point(20, 20)}, {Point(100, 100), Point(120, 120)},
        {Point(2, 2), Point(18, 18)}, {Point(102, 102), Point(118, 118)}};
    const auto groups = local_extrusion_groups(boxes);
    REQUIRE(groups.size() == 2);
    CHECK(groups[0] == std::vector<size_t>{0, 2});
    CHECK(groups[1] == std::vector<size_t>{1, 3});
}

TEST_CASE("A spanning extrusion joins local groups without reordering their members", "[LocalExtrusionGroups]")
{
    const std::vector<BoundingBox> boxes = {
        {Point(0, 0), Point(20, 20)}, {Point(100, 100), Point(120, 120)},
        {Point(10, 10), Point(110, 110)}};
    const auto groups = local_extrusion_groups(boxes);
    REQUIRE(groups.size() == 1);
    CHECK(groups.front() == std::vector<size_t>{0, 1, 2});
    CHECK(local_extrusion_groups({}).empty());
}
