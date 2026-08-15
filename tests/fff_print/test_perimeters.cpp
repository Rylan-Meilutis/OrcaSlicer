#include <catch2/catch_all.hpp>

#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/ExtrusionEntityCollection.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Print.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "test_helpers.hpp"
#include "test_utils.hpp"

using namespace Slic3r;
using namespace Slic3r::Test;

namespace {

// The layer at this Z is the last one of the base, so its top surface is the ledge.
const double ledge_z = 5.0;

// The first layer, at initial_layer_print_height.
const double first_layer_z = 0.2;

// TestMesh::step scaled 3x in X/Y: a 60x60x5 base carrying a 54x54 column up to z=10, leaving a 3mm
// top ledge around a feature that keeps rising. That is the geometry both only_one_wall_top and the
// top surface expansion act on. The ledge has to stay wider than the wall band plus two top-infill
// lines, or the expansion discards it as a sliver and the tests below assert nothing.
TriangleMesh step_with_ledge()
{
    TriangleMesh m = Slic3r::Test::mesh(TestMesh::step);
    m.scale(Vec3f(3.f, 3.f, 1.f));
    return m;
}

// A broad 40x40 floor becomes four narrow walls above Z=5, modeling the
// solid-floor-to-cavity transition found in storage boxes and trays.
TriangleMesh storage_box_transition()
{
    constexpr double size = 40.;
    constexpr double floor_height = 5.;
    constexpr double wall_height = 10.;
    constexpr double wall_width = 3.;
    TriangleMesh model = make_cube(size, size, floor_height);
    TriangleMesh wall = make_cube(wall_width, size, wall_height);
    wall.translate(0., 0., floor_height);
    model.merge(wall);
    wall = make_cube(wall_width, size, wall_height);
    wall.translate(size - wall_width, 0., floor_height);
    model.merge(wall);
    wall = make_cube(size - 2. * wall_width, wall_width, wall_height);
    wall.translate(wall_width, 0., floor_height);
    model.merge(wall);
    wall = make_cube(size - 2. * wall_width, wall_width, wall_height);
    wall.translate(wall_width, size - wall_width, floor_height);
    model.merge(wall);
    return model;
}

// Every setting the assertions depend on, so none of them rests on a default.
DynamicPrintConfig base_config(const char *wall_generator)
{
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({
        { "wall_generator",             wall_generator },
        { "layer_height",               0.2 },  // puts a layer boundary exactly on ledge_z
        { "initial_layer_print_height", 0.2 },
        { "wall_loops",                 3 },
        { "sparse_infill_density",      "15%" },
        { "top_shell_layers",           3 },
        { "bottom_shell_layers",        3 },
        { "top_surface_density",        "100%" },
        { "top_surface_expansion",      0.0 },
        { "only_one_wall_top",          false },
        { "only_one_wall_first_layer",  false },
        // Do not let the one-wall threshold discard the 3mm ledge before the feature sees it.
        { "min_width_top_surface",      0.0 },
    });
    return config;
}

double collection_length(const ExtrusionEntityCollection &coll)
{
    double len = 0.;
    for (const ExtrusionEntity *entity : coll.flatten().entities)
        if (! entity->is_collection())
            len += entity->length();
    return len;
}

// Extruded length per layer. Two slices are compared through this rather than through their G-code,
// because the G-code carries a config block that differs whenever any setting differs.
struct SliceLengths {
    std::vector<double> perimeters;
    std::vector<double> fills;
};

SliceLengths slice_lengths(const Print &print)
{
    SliceLengths out;
    for (const Layer *layer : print.objects().front()->layers()) {
        double perimeters = 0., fills = 0.;
        for (const LayerRegion *region : layer->regions()) {
            perimeters += collection_length(region->perimeters);
            fills      += collection_length(region->fills);
        }
        out.perimeters.push_back(perimeters);
        out.fills.push_back(fills);
    }
    return out;
}

double perimeter_length_at(const Print &print, double print_z)
{
    for (const Layer *layer : print.objects().front()->layers())
        if (std::abs(layer->print_z - print_z) < 1e-4) {
            double len = 0.;
            for (const LayerRegion *region : layer->regions())
                len += collection_length(region->perimeters);
            return len;
        }
    return 0.;
}

// Largest per-layer difference between two series; a negative result means they are not comparable.
double max_difference(const std::vector<double> &a, const std::vector<double> &b)
{
    if (a.size() != b.size() || a.empty())
        return -1.;
    double worst = 0.;
    for (size_t i = 0; i < a.size(); ++ i)
        worst = std::max(worst, std::abs(a[i] - b[i]));
    return worst;
}

std::vector<ExtrusionPath> smooth_outer_wall_paths(const Print &print)
{
    std::vector<ExtrusionPath> paths;
    for (const Layer *layer : print.objects().front()->layers())
        for (const LayerRegion *region : layer->regions()) {
            const ExtrusionEntityCollection flat = region->perimeters.flatten();
            for (const ExtrusionEntity *entity : flat.entities) {
                if (const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
                    path != nullptr && path->smooth_outer_wall)
                    paths.emplace_back(*path);
                else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(entity))
                    for (const ExtrusionPath &path : loop->paths)
                        if (path.smooth_outer_wall)
                            paths.emplace_back(path);
            }
        }
    return paths;
}

std::vector<ExtrusionPath> interlocking_wall_paths(const Print &print)
{
    std::vector<ExtrusionPath> paths;
    for (const Layer *layer : print.objects().front()->layers())
        for (const LayerRegion *region : layer->regions()) {
            const ExtrusionEntityCollection flat = region->perimeters.flatten();
            for (const ExtrusionEntity *entity : flat.entities) {
                if (const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
                    path != nullptr && path->nonplanar_interlocking_wall)
                    paths.emplace_back(*path);
                else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(entity))
                    for (const ExtrusionPath &path : loop->paths)
                        if (path.nonplanar_interlocking_wall)
                            paths.emplace_back(path);
            }
        }
    return paths;
}

double perimeter_volume(const Print &print)
{
    double volume = 0.;
    for (const Layer *layer : print.objects().front()->layers())
        for (const LayerRegion *region : layer->regions())
            volume += region->perimeters.total_volume();
    return volume;
}

std::vector<int> nominal_perimeter_order(const LayerRegion &region)
{
    std::vector<int> order;
    const ExtrusionEntityCollection flat = region.perimeters.flatten();
    for (const ExtrusionEntity *entity : flat.entities) {
        const auto *loop = dynamic_cast<const ExtrusionLoop *>(entity);
        if (loop == nullptr || loop->paths.empty())
            continue;
        const bool nominal = std::all_of(loop->paths.begin(), loop->paths.end(),
            [](const ExtrusionPath &path) {
                return std::all_of(path.polyline.points.begin(), path.polyline.points.end(),
                    [](const Point3 &point) { return point.z() == 0; });
            });
        if (nominal)
            order.push_back(loop->inset_idx);
    }
    return order;
}

} // namespace

TEST_CASE("Smooth outer-wall layering subdivides supported loops natively",
          "[Perimeters][SmoothOuterWall]")
{
    const char *wall_generator = GENERATE("classic", "arachne");
    const char *wall_sequence = GENERATE("outer wall/inner wall", "inner wall/outer wall");
    const double target_height = GENERATE(0.1, 0.14);
    CAPTURE(wall_generator);
    CAPTURE(wall_sequence);
    CAPTURE(target_height);
    DynamicPrintConfig config = base_config(wall_generator);
    config.set_deserialize_strict({
        {"perimeter_layering", std::string("smooth_outer_wall")},
        {"smooth_outer_wall_layer_height", target_height},
        {"wall_sequence", std::string(wall_sequence)},
    });
    Print print;
    init_and_process_print({make_cube(20., 20., 4.)}, print, config);

    DynamicPrintConfig standard_config = base_config(wall_generator);
    standard_config.set_deserialize_strict({
        {"wall_sequence", std::string(wall_sequence)},
    });
    Print standard_print;
    init_and_process_print({make_cube(20., 20., 4.)}, standard_print,
                           standard_config);

    const std::vector<ExtrusionPath> paths = smooth_outer_wall_paths(print);
    REQUIRE_FALSE(paths.empty());
    CHECK(std::all_of(paths.begin(), paths.end(), [](const ExtrusionPath &path) {
        return path.role() == erExternalPerimeter && path.z_contoured &&
            std::abs(path.height - 0.1f) < 1e-4f;
    }));
    CHECK(std::any_of(paths.begin(), paths.end(), [](const ExtrusionPath &path) {
        return std::any_of(path.polyline.points.begin(), path.polyline.points.end(),
            [](const Point3 &point) { return point.z() < 0; });
    }));
    CHECK(std::any_of(paths.begin(), paths.end(), [](const ExtrusionPath &path) {
        return std::all_of(path.polyline.points.begin(), path.polyline.points.end(),
            [](const Point3 &point) { return point.z() == 0; });
    }));
    // Subdivision changes the wall's Z schedule, not its deposited volume.
    CHECK(perimeter_volume(print) == Catch::Approx(perimeter_volume(standard_print)).epsilon(1e-6));

    // The lower independent course is dependency-ordered first, while the
    // final course remains in the selected nominal wall sequence.
    const Layer *middle = print.objects().front()->get_layer(5);
    REQUIRE(middle != nullptr);
    REQUIRE_FALSE(middle->regions().empty());
    const std::vector<int> nominal_order = nominal_perimeter_order(*middle->regions().front());
    REQUIRE(nominal_order.size() >= 3);
    if (std::string_view(wall_sequence) == "outer wall/inner wall") {
        CHECK(nominal_order[0] == 0);
        CHECK(nominal_order[1] == 1);
        CHECK(nominal_order[2] == 2);
    } else {
        CHECK(nominal_order[0] == 2);
        CHECK(nominal_order[1] == 1);
        CHECK(nominal_order[2] == 0);
    }

    ScopedTemporaryFile gcode_file(".gcode");
    GCodeProcessorResult preview;
    REQUIRE_NOTHROW(print.export_gcode(gcode_file.string(), &preview, nullptr));
    CHECK(std::any_of(preview.moves.begin(), preview.moves.end(),
        [](const GCodeProcessorResult::MoveVertex &move) {
            return move.type == EMoveType::Extrude &&
                move.extrusion_role == erSmoothOuterWall;
        }));
}

TEST_CASE("Non-planar interlocking walls remain buried and preserve outer walls",
          "[Perimeters][InterlockingWalls]")
{
    const char *wall_generator = GENERATE("classic", "arachne");
    CAPTURE(wall_generator);
    DynamicPrintConfig config = base_config(wall_generator);
    config.set_deserialize_strict({
        {"perimeter_layering", std::string("interlocking_walls")},
        {"interlocking_wall_amplitude", 0.08},
        {"interlocking_wall_wavelength", 6.0},
        {"interlocking_wall_resolution", 0.4},
    });
    Print print;
    init_and_process_print({make_cube(20., 20., 4.)}, print, config);

    const std::vector<ExtrusionPath> paths = interlocking_wall_paths(print);
    REQUIRE_FALSE(paths.empty());
    CHECK(std::all_of(paths.begin(), paths.end(), [](const ExtrusionPath &path) {
        return path.role() == erPerimeter && path.inset_idx > 0 &&
               path.z_contoured && path.polyline.points.front().z() == 0 &&
               path.polyline.points.back().z() == 0;
    }));
    CHECK(std::all_of(paths.begin(), paths.end(), [](const ExtrusionPath &path) {
        return std::all_of(path.polyline.points.begin(), path.polyline.points.end(),
            [](const Point3 &point) {
                return point.z() <= 0 && point.z() >= -scale_(0.08) - SCALED_EPSILON;
            });
    }));
    CHECK(std::any_of(paths.begin(), paths.end(), [](const ExtrusionPath &path) {
        return std::any_of(path.polyline.points.begin(), path.polyline.points.end(),
            [](const Point3 &point) { return point.z() < 0; });
    }));

    ScopedTemporaryFile gcode_file(".gcode");
    GCodeProcessorResult preview;
    REQUIRE_NOTHROW(print.export_gcode(gcode_file.string(), &preview, nullptr));
    CHECK(std::any_of(preview.moves.begin(), preview.moves.end(),
        [](const GCodeProcessorResult::MoveVertex &move) {
            return move.type == EMoveType::Extrude &&
                move.extrusion_role == erNonplanarInterlockingWall;
        }));
}

TEST_CASE("Hull-line reinforcement detects a broad storage-box floor transition", "[Perimeters][HullLine]")
{
    const char *wall_generator = GENERATE("classic", "arachne");
    CAPTURE(wall_generator);

    auto transition_perimeters_for = [wall_generator](int extra_perimeters, double expansion) {
        DynamicPrintConfig config = base_config(wall_generator);
        config.set_deserialize_strict({
            { "hull_line_extra_perimeters", extra_perimeters },
            { "hull_line_perimeter_expansion", expansion },
        });
        Print print;
        init_and_process_print({ storage_box_transition() }, print, config);
        REQUIRE_FALSE(print.objects().empty());
        return perimeter_length_at(print, ledge_z);
    };

    const double disabled = transition_perimeters_for(0, 0.0);
    const double reinforced = transition_perimeters_for(2, 0.0);
    const double expanded = transition_perimeters_for(0, 0.2);
    REQUIRE(disabled > 0.);
    CHECK(reinforced > disabled + scale_(10.));
    CHECK(std::abs(expanded - disabled) > scale_(0.5));
}

TEST_CASE("Localized shrinkage strategies expose their intended strength tradeoffs",
          "[Perimeters][LocalizedShrinkage]")
{
    const char *wall_generator = GENERATE("classic", "arachne");
    CAPTURE(wall_generator);

    auto lengths_for = [wall_generator](const char *strategy, int walls = 0,
                                        double expansion = 0., double relief = 0.) {
        DynamicPrintConfig config = base_config(wall_generator);
        config.set_deserialize_strict({
            {"localized_shrinkage_strategy", strategy},
            {"hull_line_extra_perimeters", walls},
            {"hull_line_perimeter_expansion", expansion},
            {"localized_shrinkage_infill_wall_gap", relief},
        });
        Print print;
        init_and_process_print({storage_box_transition()}, print, config);
        REQUIRE_FALSE(print.objects().empty());
        return slice_lengths(print);
    };

    const SliceLengths baseline = lengths_for("custom");
    const SliceLengths disabled = lengths_for("disabled", 2, 0.2, 0.2);
    const SliceLengths reinforced = lengths_for("reinforced_walls");
    const SliceLengths compensated = lengths_for("dimensional_compensation");
    const SliceLengths decoupled = lengths_for("reduced_wall_coupling");
    const SliceLengths balanced = lengths_for("balanced");
    const SliceLengths fully_decoupled = lengths_for("full_wall_decoupling");
    const SliceLengths sectioned = lengths_for("sectioned_solid_infill");
    const SliceLengths perforated = lengths_for("perforated_wall_relief");

    // Disabled overrides stale tuning values, while Custom keeps legacy
    // projects governed by their existing individual values.
    CHECK_THAT(max_difference(disabled.perimeters, baseline.perimeters),
               Catch::Matchers::WithinAbs(0., 1.));
    CHECK_THAT(max_difference(disabled.fills, baseline.fills),
               Catch::Matchers::WithinAbs(0., 1.));
    CHECK(max_difference(reinforced.perimeters, baseline.perimeters) > scale_(10.));
    CHECK(max_difference(compensated.perimeters, baseline.perimeters) > scale_(0.5));
    // Relief changes only the solid/wall connection, intentionally trading
    // some local bonding for less contraction transferred into the shell.
    CHECK(max_difference(decoupled.fills, baseline.fills) > scale_(0.1));
    CHECK(max_difference(balanced.perimeters, baseline.perimeters) > scale_(10.));
    CHECK(max_difference(balanced.fills, baseline.fills) > scale_(0.1));
    CHECK(max_difference(fully_decoupled.fills, decoupled.fills) > scale_(0.1));
    // Sectioning changes buried fill geometry without adding or moving walls.
    CHECK_THAT(max_difference(sectioned.perimeters, baseline.perimeters),
               Catch::Matchers::WithinAbs(0., 1.));
    CHECK(max_difference(sectioned.fills, baseline.fills) > 1.0);
    // Perforation changes only the buried fill/wall interface. Periodic
    // material bridges remain between pockets, unlike full decoupling.
    CHECK_THAT(max_difference(perforated.perimeters, baseline.perimeters),
               Catch::Matchers::WithinAbs(0., 1.));
    CHECK(max_difference(perforated.fills, baseline.fills) > 1.0);
    CHECK(max_difference(perforated.fills, fully_decoupled.fills) > 1.0);
}

TEST_CASE("Localized shrinkage walls carry a dedicated preview marker",
          "[Perimeters][LocalizedShrinkage][Preview]")
{
    const char *wall_generator = GENERATE("classic", "arachne");
    CAPTURE(wall_generator);
    DynamicPrintConfig config = base_config(wall_generator);
    config.set_deserialize_strict({
        {"localized_shrinkage_strategy", "balanced"},
        {"hull_line_extra_perimeters", 1},
        {"hull_line_perimeter_expansion", 0.15},
        {"localized_shrinkage_infill_wall_gap", 0.1},
    });
    Print print;
    init_and_process_print({storage_box_transition()}, print, config);
    REQUIRE_FALSE(print.objects().empty());

    size_t marked_paths = 0;
    const auto count_marked = [&](const auto &self,
                                  const ExtrusionEntity &entity) -> void {
        if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity))
            marked_paths += path->shrinkage_compensation;
        else if (const auto *multipath =
                     dynamic_cast<const ExtrusionMultiPath *>(&entity))
            for (const ExtrusionPath &path : multipath->paths)
                marked_paths += path.shrinkage_compensation;
        else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity))
            for (const ExtrusionPath &path : loop->paths)
                marked_paths += path.shrinkage_compensation;
        else if (const auto *collection =
                     dynamic_cast<const ExtrusionEntityCollection *>(&entity))
            for (const ExtrusionEntity *child : collection->entities)
                self(self, *child);
    };
    for (const Layer *layer : print.objects().front()->layers())
        for (const LayerRegion *region : layer->regions())
            count_marked(count_marked, region->perimeters);
    CHECK(marked_paths > 0);
}

// The expansion only retypes area as top solid infill, so it can do nothing where there is no top
// fill to begin with: zero top shell layers retypes the top surfaces as internal, and a top surface
// density of 0% leaves the top layer with walls only. The last section is the control - the same
// expansion on the same model does change the slice once a top fill exists - without which the two
// equality checks above it would hold for an unrelated reason.
TEST_CASE("Top surface expansion only acts where there is a top fill", "[Perimeters]")
{
    const char *wall_generator = GENERATE("classic", "arachne");
    CAPTURE(wall_generator);

    auto lengths_for = [wall_generator](int top_shell_layers, const char *top_surface_density, double expansion) {
        DynamicPrintConfig config = base_config(wall_generator);
        config.set_deserialize_strict({
            { "top_shell_layers",      top_shell_layers },
            { "top_surface_density",   top_surface_density },
            { "top_surface_expansion", expansion },
        });
        Print print;
        init_and_process_print({ step_with_ledge() }, print, config);
        REQUIRE_FALSE(print.objects().empty());
        return slice_lengths(print);
    };

    SECTION("no top shell layers") {
        const SliceLengths off = lengths_for(0, "100%", 0.0);
        const SliceLengths on  = lengths_for(0, "100%", 2.0);
        REQUIRE(off.perimeters.size() == on.perimeters.size());
        CHECK_THAT(max_difference(off.perimeters, on.perimeters), Catch::Matchers::WithinAbs(0., 1.0));
        CHECK_THAT(max_difference(off.fills,      on.fills),      Catch::Matchers::WithinAbs(0., 1.0));
    }

    SECTION("zero top surface density") {
        const SliceLengths off = lengths_for(3, "0%", 0.0);
        const SliceLengths on  = lengths_for(3, "0%", 2.0);
        REQUIRE(off.perimeters.size() == on.perimeters.size());
        CHECK_THAT(max_difference(off.perimeters, on.perimeters), Catch::Matchers::WithinAbs(0., 1.0));
        CHECK_THAT(max_difference(off.fills,      on.fills),      Catch::Matchers::WithinAbs(0., 1.0));
    }

    SECTION("with a top fill the same expansion does change the slice") {
        const SliceLengths off = lengths_for(3, "100%", 0.0);
        const SliceLengths on  = lengths_for(3, "100%", 2.0);
        REQUIRE(off.fills.size() == on.fills.size());
        CHECK(max_difference(off.fills, on.fills) > scale_(0.5));
    }
}

// With no top shell the top surfaces are retyped as internal, so the top surface density has nothing
// left to control: there is no top fill, and only_one_wall_top - the one route from the density to the
// perimeters - is itself switched off for want of a top surface to act on.
TEST_CASE("Top surface density does not affect a slice without a top shell", "[Perimeters]")
{
    const char *wall_generator = GENERATE("classic", "arachne");
    CAPTURE(wall_generator);

    auto lengths_for = [wall_generator](const char *top_surface_density) {
        DynamicPrintConfig config = base_config(wall_generator);
        config.set_deserialize_strict({
            { "top_shell_layers",    0 },
            { "only_one_wall_top",   true },
            { "top_surface_density", top_surface_density },
        });
        Print print;
        init_and_process_print({ step_with_ledge() }, print, config);
        REQUIRE_FALSE(print.objects().empty());
        return slice_lengths(print);
    };

    const SliceLengths solid = lengths_for("100%");
    const SliceLengths none  = lengths_for("0%");
    REQUIRE(solid.perimeters.size() == none.perimeters.size());
    CHECK_THAT(max_difference(solid.perimeters, none.perimeters), Catch::Matchers::WithinAbs(0., 1.0));
    CHECK_THAT(max_difference(solid.fills,      none.fills),      Catch::Matchers::WithinAbs(0., 1.0));
}

// On the ledge layer the inner walls are given up to the top fill, so that layer loses wall length.
// The handover needs a top fill that reaches the freed space: at a top surface density of 0% there is
// no top fill at all, and without top_surface_expansion the fill never grows over the walls. Either
// way the feature still runs, through the original generation, which keeps the inner walls up to the
// top boundary - putting that layer back between the plain and the one-wall slice.
TEST_CASE("Only one wall on top surfaces drops inner walls only where a top fill replaces them", "[Perimeters]")
{
    const char *wall_generator = GENERATE("classic", "arachne");
    CAPTURE(wall_generator);

    auto ledge_perimeters_for = [wall_generator](bool only_one_wall_top, const char *top_surface_density, double expansion) {
        DynamicPrintConfig config = base_config(wall_generator);
        config.set_deserialize_strict({
            { "only_one_wall_top",     only_one_wall_top },
            { "top_surface_density",   top_surface_density },
            { "top_surface_expansion", expansion },
        });
        Print print;
        init_and_process_print({ step_with_ledge() }, print, config);
        REQUIRE_FALSE(print.objects().empty());
        return perimeter_length_at(print, ledge_z);
    };

    const double plain              = ledge_perimeters_for(false, "100%", 2.0);
    const double one_wall           = ledge_perimeters_for(true,  "100%", 2.0);
    const double one_wall_no_fill   = ledge_perimeters_for(true,  "0%",   2.0);
    const double one_wall_no_expand = ledge_perimeters_for(true,  "100%", 0.0);

    REQUIRE(plain > 0.);
    CHECK(one_wall < plain);
    // Both fall back to the original generation, which cuts the walls back to the top boundary but not past it.
    CHECK(one_wall_no_fill > one_wall);
    CHECK(one_wall_no_fill < plain);
    CHECK(one_wall_no_expand > one_wall);
    CHECK(one_wall_no_expand < plain);
}

// The bottom counterpart: the first layer is thinned to a single wall only where a bottom shell fills the
// space behind it. With no bottom shell layers the bottom surfaces are retyped as internal, so that wall
// would ring sparse infill on the bed - the option is switched off instead, and the GUI hides it in that
// state so a profile that left it enabled cannot act behind a hidden checkbox.
TEST_CASE("Only one wall on the first layer needs a bottom shell", "[Perimeters]")
{
    const char *wall_generator = GENERATE("classic", "arachne");
    CAPTURE(wall_generator);

    auto first_layer_perimeters_for = [wall_generator](bool only_one_wall_first_layer, int bottom_shell_layers) {
        DynamicPrintConfig config = base_config(wall_generator);
        config.set_deserialize_strict({
            { "only_one_wall_first_layer", only_one_wall_first_layer },
            { "bottom_shell_layers",       bottom_shell_layers },
        });
        Print print;
        init_and_process_print({ step_with_ledge() }, print, config);
        REQUIRE_FALSE(print.objects().empty());
        return perimeter_length_at(print, first_layer_z);
    };

    const double plain             = first_layer_perimeters_for(false, 3);
    const double one_wall          = first_layer_perimeters_for(true,  3);
    // Both at zero bottom shell layers, so everything else that setting changes cancels out between them.
    const double plain_no_shell    = first_layer_perimeters_for(false, 0);
    const double one_wall_no_shell = first_layer_perimeters_for(true,  0);

    REQUIRE(plain > 0.);
    CHECK(one_wall < plain);
    // No bottom shell: the option is inert, down to the same walls an unchecked box gives.
    CHECK_THAT(one_wall_no_shell, Catch::Matchers::WithinAbs(plain_no_shell, 1.0));
}
