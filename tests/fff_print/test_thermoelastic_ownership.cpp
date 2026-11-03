#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Surface.hpp"
#include "libslic3r/SurfaceCollection.hpp"
#include "libslic3r/ThermoelasticLayers.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticOwnership.hpp"
#include "test_helpers.hpp"
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;

TEST_CASE("Native ownership uses final layer masks and explicit deposition placement", "[ThermoelasticOwnership]")
{
    const auto config = Test::multifilament_config(1, {
        {"layer_height", "0.2"}, {"initial_layer_print_height", "0.2"},
        {"thick_bridges", "0"}, {"thick_internal_bridges", "0"},
        {"brim_type", "no_brim"}, {"skirt_loops", "0"}});
    Print print;
    Model model;
    Test::init_print({Test::cube(1.)}, print, model, config);
    print.process();
    REQUIRE_FALSE(print.get_object(0)->layers().empty());
    auto &layer = *print.get_object(0)->layers().front();
    const Point shift(100000000, 200000000);
    const double z_shift = .001;
    const double plane = .001 * (layer.print_z - .5 * layer.height) + z_shift;
    const auto ownership = native_layer_ownership(layer, plane, shift, z_shift, {7});
    REQUIRE(ownership.regions.size() == 1);
    CHECK(ownership.regions.front().material == 7);
    ExPolygons expected;
    for (const auto *region : layer.regions())
        for (const auto &surface : region->slices.surfaces) {
            expected.push_back(surface.expolygon);
            expected.back().translate(shift);
        }
    CHECK(diff_ex(expected, ownership.regions.front().footprint.native()).empty());
    CHECK(diff_ex(ownership.regions.front().footprint.native(), expected).empty());
    CHECK(native_layer_ownership(layer, .001 * layer.print_z + z_shift, shift, z_shift, {7}).regions.empty());
    CHECK_FALSE(native_layer_ownership(layer, .001 * layer.bottom_z() + z_shift, shift, z_shift, {7}).regions.empty());
    CHECK_THROWS_AS(native_layer_ownership(layer, plane, shift, z_shift, {}), std::invalid_argument);
    LayerGeometryLimits limits;
    limits.max_boundary_points = 3;
    CHECK_THROWS_AS(native_layer_ownership(layer, plane, shift, z_shift, {7}, limits), ResourceLimit);
    CHECK_THROWS_AS(native_layer_ownership(layer, plane, shift, z_shift, {7}, {}, [] {
        throw std::runtime_error("cancel");
    }), std::runtime_error);
    // A slicing ownership edit must be reflected, not replaced by original CAD
    // or raw_slices. Production callers hold this generation immutable.
    for (auto *region : layer.regions())
        for (auto &surface : region->slices.surfaces)
            surface.expolygon.translate(Point(5000000, 0));
    const auto changed = native_layer_ownership(layer, plane, shift, z_shift, {7});
    REQUIRE(changed.regions.size() == 1);
    CHECK(intersection_ex(changed.regions.front().footprint.native(), ownership.regions.front().footprint.native()).empty());
}

TEST_CASE("Native ownership rejects differing role materials rather than guessing a region label", "[ThermoelasticOwnership]")
{
    const auto config = Test::multifilament_config(2, {
        {"layer_height", "0.2"}, {"initial_layer_print_height", "0.2"},
        {"thick_bridges", "0"}, {"thick_internal_bridges", "0"},
        {"inner_wall_filament_id", "2"}, {"outer_wall_filament_id", "1"},
        {"brim_type", "no_brim"}, {"skirt_loops", "0"}});
    Print print;
    Model model;
    Test::init_print({Test::cube(1.)}, print, model, config);
    print.process();
    const auto &layer = *print.get_object(0)->layers().front();
    const double plane = .001 * (layer.print_z - .5 * layer.height);
    CHECK_THROWS_AS(native_layer_ownership(layer, plane, Point(0, 0), 0., {0, 1}), std::invalid_argument);
    CHECK_THROWS_AS(native_layer_ownership(layer, plane, Point(0, 0), 0., {0}), std::invalid_argument);
    // Explicit caller-proven mechanical equivalence is different from silently
    // treating every filament of a given family as the same material.
    const auto equivalent = native_layer_ownership(layer, plane, Point(0, 0), 0., {4, 4});
    REQUIRE(equivalent.regions.size() == 1);
    CHECK(equivalent.regions.front().material == 4);
}

TEST_CASE("Native adjacent part ownership resolves bead overlaps without filling the part interior", "[ThermoelasticOwnership]")
{
    const auto config = Test::multifilament_config(2, {
        {"layer_height", "0.2"}, {"initial_layer_print_height", "0.2"},
        {"thick_bridges", "0"}, {"thick_internal_bridges", "0"},
        {"elefant_foot_compensation", "0"}, {"brim_type", "no_brim"}, {"skirt_loops", "0"}});
    Print print;
    Model model;
    Test::init_print({Test::cube(2.)}, print, model, config);
    auto &object = *model.objects.front();
    object.volumes.front()->config.set("extruder", 1);
    auto second = Test::cube(2.);
    second.translate(2., 0., 0.);
    object.add_volume(std::move(second))->config.set("extruder", 2);
    print.apply(model, config);
    print.process();
    const auto &layer = *print.get_object(0)->layers().front();
    const double plane = .001 * (layer.print_z - .5 * layer.height);
    const auto ownership = native_layer_ownership(layer, plane, Point(0, 0), 0., {0, 1});
    REQUIRE(ownership.regions.size() == 2);
    const auto left = get_extents(ownership.regions[0].footprint.native());
    const auto right = get_extents(ownership.regions[1].footprint.native());
    REQUIRE(left.max.x() == right.min.x());
    const double x = double(left.max.x()) * 1e-9;
    const double y = .5 * double(left.min.y() + left.max.y()) * 1e-9;
    const double top = layer.print_z * .001;
    const std::vector<DepositedSegment> beads{
        {{x-.0001,y-.0005,top},{x-.0001,y+.0005,top},.0004,layer.height*.001,0},
        {{x+.0001,y-.0005,top},{x+.0001,y+.0005,top},.0004,layer.height*.001,1}};
    CHECK_THROWS_AS(section_deposition(beads, plane, 0.), LayerMaterialOverlap);
    const auto resolved = section_deposition(beads, plane, 0., {}, {}, &ownership);
    REQUIRE(resolved.regions.size() == 2);
    CHECK(resolved.resolved_overlap_area_m2 > 0.);
    CHECK(intersection_ex(resolved.regions[0].footprint.native(), resolved.regions[1].footprint.native()).empty());
    double mask_area = 0.;
    for (const auto &region : ownership.regions)
        for (const auto &polygon : region.footprint.native()) mask_area += polygon.area() * 1e-18;
    CHECK(resolved.area_m2 < mask_area * .2);
    auto uniform = beads;
    uniform[1].material = 0;
    const auto original_union = section_deposition(uniform, plane, 0.);
    // Boolean clipping uses a 1 nm coordinate grid along millimetre boundaries.
    CHECK_THAT(resolved.area_m2, Catch::Matchers::WithinAbs(original_union.area_m2, 1e-12));
}

TEST_CASE("Native ownership abstains when purge or transformed paths need additional provenance", "[ThermoelasticOwnership]")
{
    const auto key = GENERATE("flush_into_objects", "flush_into_infill", "thick_bridges",
        "thick_internal_bridges", "infill_combination", "fuzzy_skin");
    CAPTURE(key);
    auto config = Test::multifilament_config(1, {
        {"layer_height", "0.2"}, {"initial_layer_print_height", "0.2"},
        {"thick_bridges", "0"}, {"thick_internal_bridges", "0"},
        {"brim_type", "no_brim"}, {"skirt_loops", "0"}});
    config.set_deserialize_strict({{key, std::string(key) == "fuzzy_skin" ? "external" : "1"}});
    Print print;
    Model model;
    Test::init_print({Test::cube(1.)}, print, model, config);
    print.process();
    const auto &layer = *print.get_object(0)->layers().front();
    const double plane = .001 * (layer.print_z - .5 * layer.height);
    CHECK_THROWS_AS(native_layer_ownership(layer, plane, Point(0, 0), 0., {0}), std::invalid_argument);
}
