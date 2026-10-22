#include <catch2/catch_test_macros.hpp>
#include "libslic3r/GCode/LocalExtrusionGroups.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/ExtrusionEntityCollection.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Point.hpp"
#include <cstddef>
#include <vector>

using namespace Slic3r;

TEST_CASE("Local grouping exposes nested sortable branches but preserves ordered paths", "[LocalExtrusionGroups]")
{
    ExtrusionEntityCollection branches;
    branches.append(ExtrusionPath(erPerimeter));
    branches.append(ExtrusionPath(erExternalPerimeter));
    ExtrusionEntityCollection ordered = branches;
    ordered.no_sort = true;
    ExtrusionEntityCollection root;
    root.append(branches);
    root.append(ordered);
    std::vector<const ExtrusionEntity *> entities;
    append_local_extrusion_entities(root, entities);
    REQUIRE(entities.size() == 3);
    CHECK(entities[0]->role() == erPerimeter);
    CHECK(entities[1]->role() == erExternalPerimeter);
    CHECK(entities[2] == root.entities[1]);
    CHECK_FALSE(entities[2]->can_sort());
}

TEST_CASE("Mixed bridge fill remains local while nested arcs preserve global order", "[LocalExtrusionGroups]")
{
    ExtrusionEntityCollection mixed;
    mixed.no_sort = true;
    mixed.append(ExtrusionPath(erSolidInfill));
    mixed.append(ExtrusionPath(erBridgeInfill));
    CHECK_FALSE(local_extrusion_requires_global_order(mixed));
    CHECK(mixed.no_sort);
    ExtrusionEntityCollection nested;
    nested.append(mixed);
    CHECK_FALSE(local_extrusion_requires_global_order(nested));
    nested.append(ExtrusionPath(erArcOverhang));
    mixed.append(nested);
    CHECK(local_extrusion_requires_global_order(mixed));
    CHECK(local_extrusion_requires_global_order(ExtrusionPath(erArcBridge)));
}

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
