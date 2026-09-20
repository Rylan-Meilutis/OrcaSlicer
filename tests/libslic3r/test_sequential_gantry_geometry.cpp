#include <catch2/catch_all.hpp>

#include <boost/filesystem.hpp>

#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/SequentialGantryGeometry.hpp"
#include "libslic3r/Format/STL.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/Utils.hpp"
#include "test_utils.hpp"

using namespace Slic3r;

namespace {

Polygon rectangle(double min_x, double min_y, double max_x, double max_y)
{
    return Polygon::new_scale({{min_x, min_y}, {max_x, min_y},
                               {max_x, max_y}, {min_x, max_y}});
}

PrintedToolpathObstacle obstacle(double top_z, double min_x, double min_y,
                                 double max_x, double max_y)
{
    return {top_z, {rectangle(min_x, min_y, max_x, max_y)}};
}

SequentialGantryGeometry square_toolhead(double height, double radius)
{
    SequentialGantryGeometry geometry;
    geometry.slices.push_back({height, false,
        {rectangle(-radius, -radius, radius, radius)}});
    return geometry;
}

} // namespace

TEST_CASE("Swept gantry geometry detects collisions between move endpoints",
          "[SequentialGantryGeometry]")
{
    const SequentialGantryGeometry gantry = square_toolhead(2., 1.);
    const std::vector<PrintedToolpathObstacle> obstacles{
        obstacle(4., 4.5, -0.5, 5.5, 0.5)};

    const ToolheadCollision collision = check_toolhead_clearance(
        Vec3d(0., 0., 1.), Vec3d(10., 0., 1.), gantry, obstacles);

    REQUIRE(collision);
    CHECK(collision.type == ToolheadCollisionType::GantryGeometry);
    CHECK(collision.obstacle_index == 0);
    CHECK(collision.slice_index == 0);
}

TEST_CASE("Material on the deposition plane does not collide with the toolhead",
          "[SequentialGantryGeometry]")
{
    const SequentialGantryGeometry gantry = square_toolhead(0., 1.);
    const std::vector<PrintedToolpathObstacle> obstacles{
        obstacle(1., 4.5, -0.5, 5.5, 0.5)};

    CHECK_FALSE(check_toolhead_clearance(
        Vec3d(0., 0., 1.), Vec3d(10., 0., 1.), gantry, obstacles));
}

TEST_CASE("Rising moves only sweep the gantry below an obstacle",
          "[SequentialGantryGeometry]")
{
    const SequentialGantryGeometry gantry = square_toolhead(0., 0.5);
    const std::vector<PrintedToolpathObstacle> low_obstacle{
        obstacle(3., 8., -0.25, 9., 0.25)};
    const std::vector<PrintedToolpathObstacle> high_obstacle{
        obstacle(3., 1., -0.25, 2., 0.25)};

    CHECK_FALSE(check_toolhead_clearance(
        Vec3d(0., 0., 1.), Vec3d(10., 0., 5.), gantry, low_obstacle));
    CHECK(check_toolhead_clearance(
        Vec3d(0., 0., 1.), Vec3d(10., 0., 5.), gantry, high_obstacle));
}

TEST_CASE("Fallback angle expands the toolhead envelope above the nozzle",
          "[SequentialGantryGeometry]")
{
    const std::vector<PrintedToolpathObstacle> obstacles{
        obstacle(8., 3., -0.25, 4., 0.25)};
    const ToolheadClearanceFallback straight{0.5, 10., 0.};
    const ToolheadClearanceFallback angled{0.5, 10., 30.};

    CHECK_FALSE(check_toolhead_clearance(
        Vec3d(0., 0., 1.), Vec3d(1., 0., 1.), {}, obstacles, straight));
    const ToolheadCollision collision = check_toolhead_clearance(
        Vec3d(0., 0., 1.), Vec3d(1., 0., 1.), {}, obstacles, angled);
    REQUIRE(collision);
    CHECK(collision.type == ToolheadCollisionType::FallbackEnvelope);
}

TEST_CASE("Non-planar fallback cone starts at the largest configured nozzle",
          "[SequentialGantryGeometry]")
{
    FullPrintConfig config;
    config.nozzle_diameter.values = {0.4, 0.6};
    config.extruder_clearance_radius.value = 40.;
    config.extruder_clearance_height_to_rod.value = 25.;
    config.nonplanar_toolhead_clearance_angle.value = 35.;

    const ToolheadClearanceFallback cone = make_nonplanar_clearance_cone(config);

    CHECK_THAT(cone.radius, Catch::Matchers::WithinAbs(0.3, EPSILON));
    CHECK_THAT(cone.height, Catch::Matchers::WithinAbs(25., EPSILON));
    CHECK_THAT(cone.angle_degrees, Catch::Matchers::WithinAbs(35., EPSILON));
}

TEST_CASE("Fallback clearance cone is rotationally symmetric",
          "[SequentialGantryGeometry]")
{
    const ToolheadClearanceFallback cone{0.2, 10., 30.};
    const double obstacle_radius = 3.;
    const std::vector<PrintedToolpathObstacle> along_x{
        obstacle(8., obstacle_radius, -0.2, obstacle_radius + 0.4, 0.2)};
    const std::vector<PrintedToolpathObstacle> along_y{
        obstacle(8., -0.2, obstacle_radius, 0.2, obstacle_radius + 0.4)};

    const bool x_collision = bool(check_toolhead_clearance(
        Vec3d(0., 0., 1.), Vec3d(0.5, 0., 1.), {}, along_x, cone));
    const bool y_collision = bool(check_toolhead_clearance(
        Vec3d(0., 0., 1.), Vec3d(0., 0.5, 1.), {}, along_y, cone));

    CHECK(x_collision == y_collision);
    CHECK(x_collision);
}

TEST_CASE("Detailed gantry geometry takes precedence over the fallback envelope",
          "[SequentialGantryGeometry]")
{
    const SequentialGantryGeometry gantry = square_toolhead(0., 0.5);
    const std::vector<PrintedToolpathObstacle> obstacles{
        obstacle(8., 3., -0.25, 4., 0.25)};
    const ToolheadClearanceFallback fallback{0.5, 10., 45.};

    CHECK_FALSE(check_toolhead_clearance(
        Vec3d(0., 0., 1.), Vec3d(1., 0., 1.), gantry, obstacles, fallback));
}

TEST_CASE("Ordered toolpaths collide with previously deposited extrusion",
          "[SequentialGantryGeometry]")
{
    const SequentialGantryGeometry gantry = square_toolhead(0., 0.5);
    const PrintedToolpathObstacle printed = make_printed_toolpath_obstacle(
        {Vec3d(4., -1., 2.), Vec3d(4., 1., 2.)}, 0.4);

    REQUIRE_THAT(printed.top_z, Catch::Matchers::WithinAbs(2., EPSILON));
    REQUIRE_FALSE(printed.footprint().empty());
    CHECK(check_toolpath_clearance(
        {Vec3d(0., 0., 1.5), Vec3d(8., 0., 1.5)}, gantry, {printed}));
    CHECK_FALSE(check_toolpath_clearance(
        {Vec3d(0., 3., 1.5), Vec3d(8., 3., 1.5)}, gantry, {printed}));
}

TEST_CASE("Symmetric clearance cone checks extrusion centerlines without polygon expansion",
          "[SequentialGantryGeometry]")
{
    const SequentialGantryGeometry gantry;
    const ToolheadClearanceFallback fallback{0.2, 10., 0.};
    const PrintedToolpathObstacle printed = make_printed_toolpath_obstacle(
        {Vec3d(4., -1., 2.), Vec3d(4., 1., 2.)}, 0.4);

    CHECK(check_toolhead_clearance(
        Vec3d(0., 0., 1.5), Vec3d(8., 0., 1.5), gantry, {printed}, fallback));
    CHECK_FALSE(check_toolhead_clearance(
        Vec3d(0., 1.5, 1.5), Vec3d(8., 1.5, 1.5), gantry, {printed}, fallback));
}

TEST_CASE("Gantry broad-phase reach only includes relevant lower slices",
          "[SequentialGantryGeometry]")
{
    SequentialGantryGeometry gantry = square_toolhead(0., 1.);
    gantry.slices.push_back({20., false, {rectangle(-10., -10., 10., 10.)}});

    CHECK_THAT(gantry.maximum_xy_reach(5.),
               Catch::Matchers::WithinAbs(std::sqrt(2.), EPSILON));
    CHECK_THAT(gantry.maximum_xy_reach(20.),
               Catch::Matchers::WithinAbs(std::sqrt(200.), EPSILON));
}

TEST_CASE("Non-planar gantry geometry uses the physical nozzle at deposition height",
          "[SequentialGantryGeometry]")
{
    SequentialGantryGeometry gantry = square_toolhead(0., 30.);
    gantry.slices.push_back({2., false, {rectangle(-40., -40., 40., 40.)}});

    gantry.set_nozzle_tip_diameter(0.4);

    CHECK_THAT(gantry.maximum_xy_reach(0.2),
               Catch::Matchers::WithinAbs(0.2, 0.01));
    CHECK_THAT(gantry.maximum_xy_reach(2.),
               Catch::Matchers::WithinAbs(std::sqrt(3200.), EPSILON));
}

TEST_CASE("Bundled CORE One geometry retains the narrow lower hotend envelope",
          "[SequentialGantryGeometry][COREOne]")
{
    struct ResourcesDirGuard {
        std::string original = resources_dir();
        ~ResourcesDirGuard() { set_resources_dir(original); }
    } resources_dir_guard;
    const boost::filesystem::path source_root =
        boost::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
    set_resources_dir((source_root / "resources").string());

    FullPrintConfig config;
    config.printer_notes.value = "PRINTER_MODEL_COREONE";

    SequentialGantryGeometry gantry = load_sequential_gantry_geometry(config);
    REQUIRE_FALSE(gantry.empty());
    gantry.set_nozzle_tip_diameter(0.4);

    CHECK(gantry.maximum_xy_reach(2.01) < 5.0);
    CHECK(gantry.maximum_xy_reach(2.5) > 20.0);
}

TEST_CASE("Every bundled gantry model supplies analyzed non-planar geometry",
          "[SequentialGantryGeometry][Regression]")
{
    struct ResourcesDirGuard {
        std::string original = resources_dir();
        ~ResourcesDirGuard() { set_resources_dir(original); }
    } resources_dir_guard;
    const boost::filesystem::path source_root =
        boost::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
    set_resources_dir((source_root / "resources").string());

    const std::array<const char *, 9> printer_notes{
        "SEQ_ARRANGE_MODEL_COREONE_INDX", "PRINTER_MODEL_MK4", "PRINTER_MODEL_MK4S",
        "PRINTER_MODEL_MK3", "PRINTER_MODEL_MINI", "PRINTER_MODEL_XL",
        "PRINTER_MODEL_COREONE_L", "PRINTER_MODEL_COREONE", "PRINTER_MODEL_HT90"};
    for (const char *notes : printer_notes) {
        CAPTURE(notes);
        FullPrintConfig config;
        config.printer_notes.value = notes;

        SequentialGantryGeometry gantry = load_sequential_gantry_geometry(config);
        REQUIRE_FALSE(gantry.empty());
        // Explicit selection must reuse the curated collision slices, not
        // reinterpret the visualization mesh or require matching printer notes.
        config.printer_notes.value.clear();
        config.sequential_print_gantry_model.value = "builtin:" + boost::filesystem::path(gantry.model_path).filename().string();
        const auto selected = load_sequential_gantry_geometry(config);
        REQUIRE(selected.slices.size() == gantry.slices.size());
        CHECK(selected.model_path == gantry.model_path);
        CHECK_THAT(selected.conservative_clearance_radius(),
                   Catch::Matchers::WithinAbs(gantry.conservative_clearance_radius(), EPSILON));
        CHECK_THAT(selected.clearance_reach().x(), Catch::Matchers::WithinAbs(gantry.clearance_reach().x(), EPSILON));
        CHECK_THAT(selected.clearance_reach().y(), Catch::Matchers::WithinAbs(gantry.clearance_reach().y(), EPSILON));
        gantry.set_nozzle_tip_diameter(0.4);
        CHECK_THAT(gantry.maximum_xy_reach(0.001),
                   Catch::Matchers::WithinAbs(0.2, 0.01));
    }
}

TEST_CASE("A user gantry STL is converted into conservative collision bands",
          "[SequentialGantryGeometry][Regression]")
{
    ScopedTemporaryFile model_file(".stl");
    TriangleMesh model(its_make_cube(10., 20., 6.));
    model.translate(Vec3f(-5.f, -10.f, 0.f));
    REQUIRE(store_stl(model_file.string().c_str(), &model, true));

    SequentialGantryGeometry gantry =
        process_sequential_gantry_model(model_file.string(), 1.);

    REQUIRE_FALSE(gantry.empty());
    CHECK(gantry.validation_error().empty());
    CHECK(gantry.slices.size() >= 6);
    CHECK_THAT(gantry.clearance_reach().x(), Catch::Matchers::WithinAbs(5., 0.02));
    CHECK_THAT(gantry.clearance_reach().y(), Catch::Matchers::WithinAbs(10., 0.02));
    FullPrintConfig imported;
    imported.sequential_print_gantry_model.value = model_file.string();
    imported.extruder_clearance_radius.value = 75.;
    const auto loaded = load_sequential_gantry_geometry(imported);
    CHECK_THAT(loaded.clearance_reach().x(), Catch::Matchers::WithinAbs(5., 0.02));
    CHECK_THAT(loaded.clearance_reach().y(), Catch::Matchers::WithinAbs(10., 0.02));
    CHECK_THAT(gantry.maximum_xy_reach(),
               Catch::Matchers::WithinAbs(std::sqrt(125.), 0.02));

    gantry.set_nozzle_tip_diameter(0.4);
    CHECK_THAT(gantry.maximum_xy_reach(0.001),
               Catch::Matchers::WithinAbs(0.2, 0.01));
    CHECK(gantry.maximum_xy_reach(0.2) > 10.);
}

TEST_CASE("An unrecognized printer profile automatically analyzes its gantry STL",
          "[SequentialGantryGeometry][Regression]")
{
    ScopedTemporaryFile model_file(".stl");
    TriangleMesh model(its_make_cube(8., 8., 4.));
    model.translate(Vec3f(-4.f, -4.f, 0.f));
    REQUIRE(store_stl(model_file.string().c_str(), &model, true));

    FullPrintConfig config;
    config.printer_notes.value = "CUSTOM_UNCATALOGUED_PRINTER";
    config.sequential_print_gantry_geometry.value.clear();
    config.sequential_print_gantry_model.value = model_file.string();
    REQUIRE(config.has("sequential_print_gantry_model"));
    REQUIRE(config.option<ConfigOptionString>("sequential_print_gantry_model") != nullptr);
    CHECK(config.option<ConfigOptionString>("sequential_print_gantry_model")->value ==
          model_file.string());

    const SequentialGantryGeometry gantry =
        load_sequential_gantry_geometry(config);
    REQUIRE_FALSE(gantry.empty());
    CHECK(gantry.validation_error().empty());
    CHECK_THAT(gantry.maximum_xy_reach(),
               Catch::Matchers::WithinAbs(std::sqrt(32.), 0.02));

    // An explicit user model also overrides a recognized bundled family.
    config.printer_notes.value = "PRINTER_MODEL_COREONE";
    const SequentialGantryGeometry overridden =
        load_sequential_gantry_geometry(config);
    REQUIRE_FALSE(overridden.empty());
    CHECK_THAT(overridden.maximum_xy_reach(),
               Catch::Matchers::WithinAbs(std::sqrt(32.), 0.02));
}

TEST_CASE("Invalid gantry slices are rejected with a profile diagnostic",
          "[SequentialGantryGeometry]")
{
    SequentialGantryGeometry gantry = square_toolhead(-1., 1.);
    CHECK_FALSE(gantry.validation_error().empty());
}

TEST_CASE("An STL whose nozzle is not centered is rejected",
          "[SequentialGantryGeometry]")
{
    ScopedTemporaryFile model_file(".stl");
    TriangleMesh model(its_make_cube(4., 4., 4.));
    model.translate(Vec3f(10.f, 10.f, 0.f));
    REQUIRE(store_stl(model_file.string().c_str(), &model, true));

    const SequentialGantryGeometry gantry =
        process_sequential_gantry_model(model_file.string());
    CHECK(gantry.empty());
    CHECK_FALSE(gantry.validation_error().empty());
}

TEST_CASE("Toolhead clearance accepts spatial-index obstacle references",
          "[SequentialGantryGeometry]")
{
    const SequentialGantryGeometry gantry = square_toolhead(2., 1.);
    const PrintedToolpathObstacle printed = obstacle(4., 4.5, -0.5, 5.5, 0.5);
    const std::vector<const PrintedToolpathObstacle *> obstacles{&printed};

    const ToolheadCollision collision = check_toolhead_clearance(
        Vec3d(0., 0., 1.), Vec3d(10., 0., 1.), gantry, obstacles);

    REQUIRE(collision);
    CHECK(collision.type == ToolheadCollisionType::GantryGeometry);
    CHECK(collision.obstacle_index == 0);
}

TEST_CASE("Unmerged extrusion footprints preserve bent-path occupancy",
          "[SequentialGantryGeometry]")
{
    const SequentialGantryGeometry gantry = square_toolhead(0., 0.2);
    const PrintedToolpathObstacle printed = make_printed_toolpath_obstacle(
        {Vec3d(0., 0., 2.), Vec3d(4., 0., 2.), Vec3d(4., 4., 2.)}, 0.4);

    CHECK(check_toolpath_clearance(
        {Vec3d(3., 2., 1.5), Vec3d(5., 2., 1.5)}, gantry, {printed}));
    CHECK_FALSE(check_toolpath_clearance(
        {Vec3d(0., 3., 1.5), Vec3d(2., 3., 1.5)}, gantry, {printed}));
}
