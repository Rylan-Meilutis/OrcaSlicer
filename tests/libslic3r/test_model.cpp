#include <catch2/catch_all.hpp>

#include "libslic3r/Model.hpp"

using namespace Slic3r;

TEST_CASE("Invalid material mappings leave the model unchanged", "[Model][ToolMapping]")
{
    Model model;
    ModelObject *object = model.add_object();
    object->config.set("extruder", 1);
    const int invalid = GENERATE(-1, int(EnforcerBlockerType::ExtruderMax));
    CHECK_THROWS_AS(remap_model_filament_slots(model, {{0, 7}, {1, invalid}}), std::invalid_argument);
    CHECK(object->config.extruder() == 1);
    CHECK_THROWS_AS(remap_model_filament_slots(model, {{invalid, 0}}), std::invalid_argument);
    CHECK(object->config.extruder() == 1);
}

TEST_CASE("Identity material mappings do not invalidate slicing", "[Model][ToolMapping]")
{
    Model model;
    ModelObject *object = model.add_object();
    ModelVolume *volume = object->add_volume(make_cube(20., 20., 20.));
    object->config.set("extruder", 1);
    TriangleSelector selector(volume->mesh());
    selector.set_facet(0, EnforcerBlockerType::Extruder1);
    volume->mmu_segmentation_facets.set(selector);
    const ModelConfig &config = object->config;
    const auto config_timestamp = config.timestamp();
    const auto painting_timestamp = volume->mmu_segmentation_facets.timestamp();

    remap_model_filament_slots(model, {{0, 0}, {1, 1}});

    CHECK(config.timestamp() == config_timestamp);
    CHECK(volume->mmu_segmentation_facets.timestamp() == painting_timestamp);
}

TEST_CASE("Slice model material remapping preserves project painting and role assignments", "[Model][ToolMapping]")
{
    Model project;
    ModelObject *object = project.add_object();
    ModelVolume *volume = object->add_volume(make_cube(20., 20., 20.));
    object->config.set("extruder", 1);
    object->config.set("support_filament", 2);
    volume->config.set("extruder", 2);
    volume->config.set("inner_wall_filament_id", 1);
    volume->config.set("outer_wall_filament_id", 0);
    volume->config.set("sparse_infill_filament_id", -1);
    object->layer_config_ranges[{1., 2.}].set("top_surface_filament_id", 2);
    TriangleSelector selector(volume->mesh());
    selector.set_facet(0, EnforcerBlockerType::Extruder1);
    selector.set_facet(1, EnforcerBlockerType::Extruder2);
    REQUIRE(volume->mmu_segmentation_facets.set(selector));

    Model slice_model(project);
    // Swap logical material 2 onto T0 while material 1 moves onto T7. The
    // relocation must be simultaneous, not cascade through the new slot 1.
    remap_model_filament_slots(slice_model, {{0, 7}, {1, 0}});
    const auto *mapped = slice_model.objects.front();
    const auto *mapped_volume = mapped->volumes.front();
    CHECK(mapped->config.extruder() == 8);
    CHECK(mapped->config.opt_int("support_filament") == 1);
    CHECK(mapped_volume->config.extruder() == 1);
    CHECK(mapped_volume->config.opt_int("inner_wall_filament_id") == 8);
    CHECK(mapped_volume->config.opt_int("outer_wall_filament_id") == 0);
    CHECK(mapped_volume->config.opt_int("sparse_infill_filament_id") == -1);
    CHECK(mapped->layer_config_ranges.begin()->second.opt_int("top_surface_filament_id") == 1);
    CHECK(mapped_volume->mmu_segmentation_facets.has_facets(*mapped_volume, EnforcerBlockerType::Extruder8));
    CHECK(mapped_volume->mmu_segmentation_facets.has_facets(*mapped_volume, EnforcerBlockerType::Extruder1));
    CHECK_FALSE(mapped_volume->mmu_segmentation_facets.has_facets(*mapped_volume, EnforcerBlockerType::Extruder2));
    CHECK(object->config.extruder() == 1);
    CHECK(object->config.opt_int("support_filament") == 2);
    CHECK(volume->config.extruder() == 2);
    CHECK(volume->mmu_segmentation_facets.has_facets(*volume, EnforcerBlockerType::Extruder2));
    CHECK(object->layer_config_ranges.begin()->second.opt_int("top_surface_filament_id") == 2);
}

TEST_CASE("Material remapping includes layer tool changes without rewriting custom G-code", "[Model][ToolMapping]")
{
    Model project;
    for (int plate : {0, 1}) {
        auto &info = project.plates_custom_gcodes[plate];
        info.mode = CustomGCode::MultiExtruder;
        info.gcodes = {
            {1., CustomGCode::ToolChange, 1, "", ""},
            {2., CustomGCode::ColorChange, 2, "#112233", ""},
            {3., CustomGCode::ToolChange, 0, "", ""},
            {4., CustomGCode::ToolChange, 3, "", ""},
            {5., CustomGCode::Custom, 1, "", "T0\nM104 T0 S210"},
            {6., CustomGCode::PausePrint, 2, "", "Check print"}
        };
    }
    Model job(project);
    remap_model_filament_slots(job, {{0, 7}, {1, 0}});
    for (int plate : {0, 1}) {
        const auto &original = project.plates_custom_gcodes.at(plate);
        const auto &mapped = job.plates_custom_gcodes.at(plate);
        CHECK(mapped.mode == original.mode);
        REQUIRE(mapped.gcodes.size() == original.gcodes.size());
        CHECK(mapped.gcodes[0].extruder == 8);
        CHECK(mapped.gcodes[1].extruder == 1);
        CHECK(mapped.gcodes[1].color == "#112233");
        CHECK(original.gcodes[0].extruder == 1);
        CHECK(original.gcodes[1].extruder == 2);
        for (size_t i = 2; i < mapped.gcodes.size(); ++i)
            CHECK(mapped.gcodes[i] == original.gcodes[i]);
    }
}

// convex_hull_2d does not clip geometry below the bed, so these cases avoid
// sinking transforms.
TEST_CASE("A part's 2D convex hull is its footprint projected onto the bed", "[Model]")
{
    Model model;
    ModelObject* object = model.add_object();
    // Keep the cube's raw coordinates ([0,20] on every axis): the default
    // add_volume re-centers the geometry, which would move the footprint.
    object->add_volume(make_cube(20, 20, 20), ModelVolumeType::MODEL_PART, false);

    SECTION("identity transform yields the 20 mm square") {
        const Polygon hull   = object->convex_hull_2d(Geometry::Transformation{}.get_matrix());
        const BoundingBox bb = hull.bounding_box();
        CHECK(hull.size() == 4);
        CHECK(bb.min.x() == scaled(0.));
        CHECK(bb.min.y() == scaled(0.));
        CHECK(bb.max.x() == scaled(20.));
        CHECK(bb.max.y() == scaled(20.));
    }

    SECTION("scaling and offset move and grow the footprint") {
        Geometry::Transformation t;
        t.set_scaling_factor({2, 2, 2}); // cube now spans [0,40]
        t.set_offset({10, 5, 0});        // then shift +10 in X, +5 in Y

        const Polygon hull   = object->convex_hull_2d(t.get_matrix());
        const BoundingBox bb = hull.bounding_box();
        CHECK(hull.size() == 4);
        CHECK(bb.min.x() == scaled(10.));
        CHECK(bb.min.y() == scaled(5.));
        CHECK(bb.max.x() == scaled(50.));
        CHECK(bb.max.y() == scaled(45.));
    }
}
