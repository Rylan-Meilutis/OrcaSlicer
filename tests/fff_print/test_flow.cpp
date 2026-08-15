#include <catch2/catch_all.hpp>

#include <numeric>
#include <limits>
#include <set>
#include <sstream>

#include "test_helpers.hpp" // get access to init_print, etc

#include "libslic3r/Config.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Config.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/libslic3r.h"
#include "test_utils.hpp"

using namespace Slic3r::Test;
using namespace Slic3r;

TEST_CASE("Feature flows use their selected tool nozzle diameters",
          "[Flow][MultiNozzle][FeatureTools]")
{
    Print print;
    Slic3r::Test::init_and_process_print({cube(10.)}, print, {
        {"nozzle_diameter", "0.2,0.6,0.8"},
        {"filament_diameter", "1.75,1.75,1.75"},
        {"layer_height", 0.1},
        {"initial_layer_print_height", 0.1},
        {"outer_wall_filament_id", 1},
        {"inner_wall_filament_id", 2},
        {"sparse_infill_filament_id", 3},
        {"outer_wall_line_width", 0.0},
        {"inner_wall_line_width", 0.0},
        {"sparse_infill_line_width", 0.0},
        {"wall_loops", 3},
        {"sparse_infill_density", "15%"},
    });

    const PrintObject &object = *print.objects().front();
    const PrintRegion &region = object.printing_region(0);
    const Flow outer = region.flow(object, frExternalPerimeter, 0.1, false);
    const Flow inner = region.flow(object, frPerimeter, 0.1, false);
    const Flow infill = region.flow(object, frInfill, 0.1, false);

    CHECK(outer.nozzle_diameter() == Catch::Approx(0.2));
    CHECK(inner.nozzle_diameter() == Catch::Approx(0.6));
    CHECK(infill.nozzle_diameter() == Catch::Approx(0.8));
    CHECK(outer.width() < inner.width());
    CHECK(inner.width() < infill.width());
}

static double total_positive_extrusion(const std::string &gcode)
{
    double total = 0.;
    GCodeReader reader;
    reader.parse_buffer(gcode, [&total](GCodeReader &self, const GCodeReader::GCodeLine &line) {
        if (line.extruding(self))
            total += line.dist_E(self);
    });
    return total;
}

static std::vector<std::pair<const Layer *, const ExtrusionPath *>> staggered_wall_paths(const Print &print)
{
    std::vector<std::pair<const Layer *, const ExtrusionPath *>> result;
    const auto collect = [&result](const Layer *layer, const ExtrusionEntity &entity, const auto &self) -> void {
        if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity)) {
            if (path->staggered_perimeter)
                result.emplace_back(layer, path);
        } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
            for (const ExtrusionPath &path : loop->paths)
                if (path.staggered_perimeter)
                    result.emplace_back(layer, &path);
        } else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
            for (const ExtrusionPath &path : multipath->paths)
                if (path.staggered_perimeter)
                    result.emplace_back(layer, &path);
        } else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
            for (const ExtrusionEntity *child : collection->entities)
                self(layer, *child, self);
        }
    };
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions())
                collect(layer, region->perimeters, collect);
    return result;
}

TEST_CASE("Staggered perimeters raise only supported odd inner walls", "[Flow][StaggeredPerimeters]")
{
    for (const std::string generator : {"classic", "arachne"}) {
        CAPTURE(generator);
        Print print;
        Slic3r::Test::init_and_process_print({cube(10.)}, print, {
            {"layer_height", 0.2},
            {"initial_layer_print_height", 0.2},
            {"wall_generator", generator},
            {"wall_loops", 4},
            {"top_shell_layers", 0},
            {"bottom_shell_layers", 0},
            {"sparse_infill_density", "0%"},
            {"staggered_perimeters", 1},
            {"staggered_perimeter_offset", "50%"}
        });

        const auto paths = staggered_wall_paths(print);
        REQUIRE_FALSE(paths.empty());
        bool found_course_entry = false;
        for (const auto &[layer, path_ptr] : paths) {
            const ExtrusionPath &path = *path_ptr;
            CHECK(layer->lower_layer != nullptr);
            CHECK(path.role() == erPerimeter);
            CHECK(path.inset_idx > 0);
            CHECK(path.inset_idx % 2 == 1);
            CHECK(path.z_contoured);
            CHECK_FALSE(path.nonplanar_surface);
            REQUIRE_FALSE(path.polyline.points.empty());
            // Entry and exit transition beads are still brick-course
            // geometry. Neither may be exposed as a visible top surface.
            REQUIRE(layer->upper_layer != nullptr);
            const Polygons footprint = path.polygons_covered_by_width(
                float(SCALED_EPSILON));
            REQUIRE_FALSE(footprint.empty());
            CHECK(diff_ex(footprint, layer->upper_layer->lslices,
                          ApplySafetyOffset::No).empty());
            const bool raised = std::any_of(
                path.polyline.points.begin(), path.polyline.points.end(),
                [](const Point3 &point) { return point.z() > 0; });
            if (raised) {
                CHECK(std::all_of(path.polyline.points.begin(), path.polyline.points.end(),
                    [](const Point3 &point) {
                        return std::abs(point.z() - coord_t(scale_(0.1))) <= 1;
                    }));
                found_course_entry |= path.height > 0.2f + float(EPSILON);
            } else {
                FAIL("Nominal exit transitions must not use the brick preview role");
            }
        }
        CHECK(found_course_entry);
    }
}

TEST_CASE("Staggered perimeters can include buried outer walls",
          "[Flow][StaggeredPerimeters]")
{
    Print print;
    Slic3r::Test::init_and_process_print({cube(10.)}, print, {
        {"layer_height", 0.2},
        {"initial_layer_print_height", 0.2},
        {"wall_loops", 4},
        {"top_shell_layers", 0},
        {"bottom_shell_layers", 0},
        {"sparse_infill_density", "0%"},
        {"staggered_perimeters", 1},
        {"staggered_perimeters_inner_only", 0},
        {"staggered_perimeter_offset", "50%"}
    });

    const auto paths = staggered_wall_paths(print);
    REQUIRE_FALSE(paths.empty());
    bool found_outer_wall = false;
    bool found_outer_entry_transition = false;
    for (const auto &[layer, path] : paths) {
        if (path->inset_idx != 0)
            continue;
        found_outer_wall = true;
        CHECK(path->role() == erExternalPerimeter);
        CHECK(layer->lower_layer != nullptr);
        CHECK(layer->upper_layer != nullptr);
        const Polygons footprint = path->polygons_covered_by_width(
            float(SCALED_EPSILON));
        REQUIRE_FALSE(footprint.empty());
        const ExPolygons changed = xor_ex(
            layer->lslices, layer->upper_layer->lslices,
            ApplySafetyOffset::No);
        CHECK((changed.empty() || intersection_ex(
            footprint, changed, ApplySafetyOffset::No).empty()));
        found_outer_entry_transition |=
            path->height > layer->height + float(EPSILON);
    }
    CHECK(found_outer_wall);
    CHECK(found_outer_entry_transition);
}

TEST_CASE("Staggered wall entities are emitted from lower to higher Z",
          "[Flow][StaggeredPerimeters][Ordering]")
{
    Print print;
    Slic3r::Test::init_and_process_print({cube(10.)}, print, {
        {"layer_height", 0.2},
        {"initial_layer_print_height", 0.2},
        {"wall_loops", 4},
        {"top_shell_layers", 0},
        {"bottom_shell_layers", 0},
        {"sparse_infill_density", "0%"},
        {"staggered_perimeters", 1},
        {"staggered_perimeters_inner_only", 0},
        {"staggered_perimeter_offset", "50%"}
    });

    const auto maximum_z = [](const ExtrusionEntity &entity, const auto &self) -> coord_t {
        coord_t result = 0;
        const auto inspect_path = [&result](const ExtrusionPath &path) {
            for (const Point3 &point : path.polyline.points)
                result = std::max(result, point.z());
        };
        if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity))
            inspect_path(*path);
        else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity))
            for (const ExtrusionPath &path : loop->paths)
                inspect_path(path);
        else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity))
            for (const ExtrusionPath &path : multipath->paths)
                inspect_path(path);
        else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity))
            for (const ExtrusionEntity *child : collection->entities)
                result = std::max(result, self(*child, self));
        return result;
    };
    bool checked_dependency = false;
    const auto inspect_collection = [&maximum_z, &checked_dependency](
        const ExtrusionEntity &entity, const auto &self) -> void {
        const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity);
        if (collection == nullptr)
            return;
        coord_t previous = 0;
        bool has_raised = false;
        for (const ExtrusionEntity *child : collection->entities) {
            const coord_t current = maximum_z(*child, maximum_z);
            CHECK(current >= previous);
            previous = current;
            has_raised |= current > 0;
            self(*child, self);
        }
        checked_dependency |= has_raised && collection->entities.size() > 1;
    };
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions())
                inspect_collection(region->perimeters, inspect_collection);

    CHECK(checked_dependency);
}

TEST_CASE("Brick transition beads do not activate non-planar travel clearance",
          "[Flow][StaggeredPerimeters][GCode][Regression]")
{
    const std::string gcode = slice({cube(10.)}, {
        {"layer_height", 0.2},
        {"initial_layer_print_height", 0.2},
        {"wall_loops", 4},
        {"top_shell_layers", 0},
        {"bottom_shell_layers", 0},
        {"sparse_infill_density", "15%"},
        {"staggered_perimeters", 1},
        {"staggered_perimeter_offset", "50%"}
    });

    CHECK(gcode.find(";TYPE:Brick wall") != std::string::npos);
    CHECK(gcode.find("raise for non-planar toolhead clearance") == std::string::npos);
}

TEST_CASE("Staggered inner walls continue through shrinking outlines", "[Flow][StaggeredPerimeters][Regression]")
{
    Print print;
    Slic3r::Test::init_and_process_print(
        {TriangleMesh(its_make_frustum(10., 10., PI / 32.))}, print, {
        {"layer_height", 0.2},
        {"initial_layer_print_height", 0.2},
        {"wall_loops", 3},
        {"top_shell_layers", 0},
        {"bottom_shell_layers", 0},
        {"sparse_infill_density", "0%"},
        {"staggered_perimeters", 1},
        {"staggered_perimeter_offset", "50%"}
    });

    const auto paths = staggered_wall_paths(print);
    std::set<size_t> staggered_layers;
    for (const auto &[layer, path] : paths)
        staggered_layers.insert(layer->id());

    // A 10 mm frustum at 0.2 mm has 50 layers. Every non-boundary layer
    // retains an odd inner wall even though the next outline is smaller.
    CHECK(staggered_layers.size() >= 48);
}

TEST_CASE("Staggered perimeters split partially covered wall paths", "[Flow][StaggeredPerimeters][Regression]")
{
    TriangleMesh model = make_cube(20., 20., 2.);
    TriangleMesh upper = make_cube(10., 20., 2.);
    upper.translate(0., 0., 2.);
    model.merge(upper);

    Print print;
    Slic3r::Test::init_and_process_print({model}, print, {
        {"layer_height", 0.2},
        {"initial_layer_print_height", 0.2},
        {"wall_generator", "arachne"},
        {"wall_loops", 4},
        {"top_shell_layers", 0},
        {"bottom_shell_layers", 0},
        {"sparse_infill_density", "0%"},
        {"staggered_perimeters", 1},
        {"staggered_perimeter_offset", "50%"}
    });

    bool found_partially_staggered_path = false;
    bool found_contiguous_transition = false;
    bool found_half_height_ramp = false;
    bool found_low_course_first = false;
    const auto inspect_paths = [&found_partially_staggered_path,
                                &found_contiguous_transition,
                                &found_half_height_ramp,
                                &found_low_course_first](const ExtrusionPaths &paths) {
        const auto is_raised = [](const ExtrusionPath &path) {
            return path.staggered_perimeter && std::any_of(
                path.polyline.points.begin(), path.polyline.points.end(),
                [](const Point3 &point) { return point.z() > 0; });
        };
        const bool has_staggered = std::any_of(paths.begin(), paths.end(), is_raised);
        const bool has_nominal_inner = std::any_of(paths.begin(), paths.end(),
            [&is_raised](const ExtrusionPath &path) {
                return !is_raised(path) && path.role() == erPerimeter &&
                       path.inset_idx > 0 && path.inset_idx % 2 == 1;
            });
        if (!has_staggered || !has_nominal_inner)
            return;
        found_partially_staggered_path = true;
        const auto path_min_z = [](const ExtrusionPath &path) {
            coord_t minimum = std::numeric_limits<coord_t>::max();
            for (const Point3 &point : path.polyline.points)
                minimum = std::min(minimum, point.z());
            return minimum;
        };
        const auto lowest_path = std::min_element(
            paths.begin(), paths.end(),
            [&path_min_z](const ExtrusionPath &left, const ExtrusionPath &right) {
                return path_min_z(left) < path_min_z(right);
            });
        const coord_t course_minimum = path_min_z(*lowest_path);
        found_low_course_first |= path_min_z(paths.front()) == course_minimum;
        for (size_t idx = 1; idx < paths.size(); ++idx) {
            const ExtrusionPath &before = paths[idx - 1];
            const ExtrusionPath &after = paths[idx];
            if (is_raised(before) != is_raised(after) &&
                before.last_point3().to_point() == after.first_point3().to_point()) {
                found_contiguous_transition = true;
                CHECK(before.last_point3().z() == after.first_point3().z());
                const ExtrusionPath &staggered = is_raised(before) ? before : after;
                const auto z_extents = std::minmax_element(
                    staggered.polyline.points.begin(), staggered.polyline.points.end(),
                    [](const Point3 &left, const Point3 &right) {
                        return left.z() < right.z();
                    });
                const coord_t minimum_z = z_extents.first->z();
                const coord_t maximum_z = z_extents.second->z();
                if (minimum_z != maximum_z)
                    CHECK(staggered.staggered_transition);
                found_half_height_ramp |= std::any_of(
                    staggered.polyline.points.begin(), staggered.polyline.points.end(),
                    [minimum_z, maximum_z](const Point3 &point) {
                        return point.z() > minimum_z && point.z() < maximum_z;
                    });
            }
        }
    };
    const auto visit = [&inspect_paths](const ExtrusionEntity &entity, const auto &self) -> void {
        if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity))
            inspect_paths(loop->paths);
        else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity))
            inspect_paths(multipath->paths);
        else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity))
            for (const ExtrusionEntity *child : collection->entities)
                self(*child, self);
    };
    for (const PrintObject *object : print.objects())
        for (const Layer *layer : object->layers())
            for (const LayerRegion *region : layer->regions())
                visit(region->perimeters, visit);

    CHECK(found_partially_staggered_path);
    CHECK(found_contiguous_transition);
    CHECK(found_half_height_ramp);
    CHECK(found_low_course_first);
}

TEST_CASE("Staggered perimeters use the existing intermediate inner wall flow", "[Flow][StaggeredPerimeters]")
{
    const auto extrusion_for = [](bool enabled, double flow) {
        return total_positive_extrusion(slice({cube(10.)}, {
            {"layer_height", 0.2},
            {"initial_layer_print_height", 0.2},
            {"wall_loops", 4},
            {"top_shell_layers", 0},
            {"bottom_shell_layers", 0},
            {"sparse_infill_density", "0%"},
            {"staggered_perimeters", enabled ? "1" : "0"},
            {"inner_walls_flow_ratio", std::to_string(flow * 100.) + "%"}
        }));
    };

    CHECK_THAT(extrusion_for(true, 1.00),
               Catch::Matchers::WithinAbs(extrusion_for(false, 1.00), 1e-4));
    CHECK(extrusion_for(false, 1.10) > extrusion_for(false, 1.00));
    CHECK(extrusion_for(true, 1.10) > extrusion_for(true, 1.00));
}

TEST_CASE("Intermediate inner wall flow changes only prints with at least three walls", "[Flow][Regression]")
{
    const auto extrusion_for = [](int walls, double inner_walls_flow) {
        return total_positive_extrusion(slice({cube(10.)}, {
            {"wall_generator", "classic"},
            {"wall_loops", walls},
            {"top_shell_layers", 0},
            {"bottom_shell_layers", 0},
            {"sparse_infill_density", "0%"},
            {"set_other_flow_ratios", 0},
            {"inner_walls_flow_ratio", std::to_string(inner_walls_flow * 100.) + "%"}
        }));
    };

    CHECK_THAT(extrusion_for(2, 1.5), Catch::Matchers::WithinAbs(extrusion_for(2, 1.0), 1e-5));
    CHECK(extrusion_for(3, 1.5) > extrusion_for(3, 1.0));
}

TEST_CASE("Intermediate inner wall flow preserves first and exposed top layers", "[Flow][Regression]")
{
    const auto extrusion_for = [](double height, double inner_walls_flow) {
        return total_positive_extrusion(slice({make_cube(10., 10., height)}, {
            {"layer_height", 0.2},
            {"initial_layer_print_height", 0.2},
            {"wall_generator", "classic"},
            {"wall_loops", 3},
            {"top_shell_layers", 3},
            {"bottom_shell_layers", 3},
            {"sparse_infill_density", "0%"},
            {"set_other_flow_ratios", 0},
            {"inner_walls_flow_ratio", std::to_string(inner_walls_flow * 100.) + "%"}
        }));
    };

    CHECK_THAT(extrusion_for(0.4, 1.5), Catch::Matchers::WithinAbs(extrusion_for(0.4, 1.0), 1e-5));
    CHECK(extrusion_for(2.0, 1.5) > extrusion_for(2.0, 1.0));
}

TEST_CASE("Intermediate inner wall flow is represented by the G-code preview width", "[Flow][Preview][Regression]")
{
    const auto maximum_preview_width = [](double inner_walls_flow) {
        Print print;
        Slic3r::Test::init_and_process_print({make_cube(10., 10., 2.)}, print, {
            {"layer_height", 0.2},
            {"initial_layer_print_height", 0.2},
            {"wall_generator", "classic"},
            {"wall_loops", 3},
            {"top_shell_layers", 0},
            {"bottom_shell_layers", 0},
            {"sparse_infill_density", "0%"},
            {"set_other_flow_ratios", 0},
            {"inner_walls_flow_ratio", std::to_string(inner_walls_flow * 100.) + "%"}
        });

        ScopedTemporaryFile gcode_file(".gcode");
        GCodeProcessorResult preview;
        print.export_gcode(gcode_file.string(), &preview, nullptr);

        float maximum = 0.f;
        for (const auto &move : preview.moves)
            if (move.type == EMoveType::Extrude && move.extrusion_role == erPerimeter)
                maximum = std::max(maximum, move.width);
        return maximum;
    };

    const float normal_width = maximum_preview_width(1.);
    CHECK(normal_width > 0.f);
    CHECK(maximum_preview_width(1.5) > normal_width);
}

/// Test the expected behavior for auto-width,
/// spacing, etc
SCENARIO("Flow math for non-bridges", "[Flow]") {
    GIVEN("Nozzle Diameter of 0.4, a desired width of 1mm and layer height of 0.5") {
        ConfigOptionFloatOrPercent	width(1.0, false);
        float nozzle_diameter	= 0.4f;
        float layer_height		= 0.4f;

        // Spacing for non-bridges is has some overlap
        THEN("External perimeter flow has spacing fixed to 1.125 * nozzle_diameter") {
            auto flow = Flow::new_from_config_width(frExternalPerimeter, ConfigOptionFloatOrPercent(0, false), nozzle_diameter, layer_height);
            REQUIRE(flow.spacing() == Catch::Approx(1.125 * nozzle_diameter - layer_height * (1.0 - PI / 4.0)));
        }

        THEN("Internal perimeter flow has spacing fixed to 1.125 * nozzle_diameter") {
            auto flow = Flow::new_from_config_width(frPerimeter, ConfigOptionFloatOrPercent(0, false), nozzle_diameter, layer_height);
            REQUIRE(flow.spacing() == Catch::Approx(1.125 *nozzle_diameter - layer_height * (1.0 - PI / 4.0)));
        }
        THEN("Spacing for supplied width is 0.8927f") {
            auto flow = Flow::new_from_config_width(frExternalPerimeter, width, nozzle_diameter, layer_height);
            REQUIRE(flow.spacing() == Catch::Approx(width.value - layer_height * (1.0 - PI / 4.0)));
            flow = Flow::new_from_config_width(frPerimeter, width, nozzle_diameter, layer_height);
            REQUIRE(flow.spacing() == Catch::Approx(width.value - layer_height * (1.0 - PI / 4.0)));
        }
    }
    /// Check the min/max
    GIVEN("Nozzle Diameter of 0.25") {
        float nozzle_diameter	= 0.25f;
        float layer_height		= 0.5f;
        WHEN("layer height is set to 0.2") {
            layer_height = 0.15f;
            THEN("Max width is set.") {
                auto flow = Flow::new_from_config_width(frPerimeter, ConfigOptionFloatOrPercent(0, false), nozzle_diameter, layer_height);
                REQUIRE(flow.width() == Catch::Approx(1.125 * nozzle_diameter));
            }
        }
        WHEN("Layer height is set to 0.25") {
            layer_height = 0.25f;
            THEN("Min width is set.") {
                auto flow = Flow::new_from_config_width(frPerimeter, ConfigOptionFloatOrPercent(0, false), nozzle_diameter, layer_height);
                REQUIRE(flow.width() == Catch::Approx(1.125 * nozzle_diameter));
            }
        }
    }

#if 0
    /// Check for an edge case in the maths where the spacing could be 0; original
    /// math is 0.99. Slic3r issue #4654
    GIVEN("Input spacing of 0.414159 and a total width of 2") {
        double in_spacing = 0.414159;
        double total_width = 2.0;
        auto flow = Flow::new_from_spacing(1.0, 0.4, 0.3);
        WHEN("solid_spacing() is called") {
            double result = flow.solid_spacing(total_width, in_spacing);
            THEN("Yielded spacing is greater than 0") {
                REQUIRE(result > 0);
            }
        }
    }
#endif    

}

/// Spacing, width calculation for bridge extrusions
SCENARIO("Flow math for bridges", "[Flow]") {
    GIVEN("Nozzle Diameter of 0.4, a desired width of 1mm and layer height of 0.5") {
		float nozzle_diameter	= 0.4f;
		float bridge_flow		= 1.0f;
        WHEN("Flow role is frExternalPerimeter") {
            auto flow = Flow::bridging_flow(nozzle_diameter * sqrt(bridge_flow), nozzle_diameter);
            THEN("Bridge width is same as nozzle diameter") {
                REQUIRE(flow.width() == Catch::Approx(nozzle_diameter));
            }
            THEN("Bridge spacing is same as nozzle diameter + BRIDGE_EXTRA_SPACING") {
                REQUIRE(flow.spacing() == Catch::Approx(nozzle_diameter + BRIDGE_EXTRA_SPACING));
            }
        }
    }
}
