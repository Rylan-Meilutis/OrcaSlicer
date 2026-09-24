#include <catch2/catch_all.hpp>

#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/ExtrusionEntityCollection.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/TriangleMesh.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <set>
#include <vector>

#include "test_helpers.hpp"
#include "test_utils.hpp"

using namespace Slic3r;
using namespace Slic3r::Test;

TEST_CASE("Extra overhang wall overlap adds inward courses without moving outer walls", "[Perimeters][OverhangOverlap]")
{
    const std::string generator = GENERATE(std::string("classic"), std::string("arachne"));
    const double arc_threshold = GENERATE(-1., 0., 5., 100.);
    double last_arc_length = 0.;
    auto paths_for = [&](double overlap, bool enabled) {
        TriangleMesh model = make_cube(10., 20., 5.);
        // The mixed case has a 10 mm cantilever on one side and a 4 mm
        // cantilever on the other: only the latter should receive extra walls.
        TriangleMesh roof = make_cube(arc_threshold == 5. ? 24. : 30., 20., 2.);
        roof.translate(-10., 0., 5.);
        model.merge(roof);
        Print print;
        init_and_process_print({model}, print, {
            {"wall_generator", generator}, {"overhang_wall_overlap", overlap},
            {"extra_perimeters_on_overhangs", enabled}, {"detect_overhang_wall", true},
            {"arc_overhang_enabled", arc_threshold >= 0.},
            {"arc_overhang_bridges", true}, {"arc_overhang_overhangs", true},
            {"arc_overhang_bridge_distance", std::max(0., arc_threshold)},
            {"arc_overhang_min_overhang_distance", std::max(0., arc_threshold)},
            {"enable_support", false},
            {"wall_loops", 3}, {"layer_height", 0.2}, {"initial_layer_print_height", 0.2}
        });
        double length = 0.;
        last_arc_length = 0.;
        std::vector<Points> outer;
        for (const Layer *layer : print.objects().front()->layers())
            for (const LayerRegion *region : layer->regions()) {
                const auto entities = region->perimeters.flatten();
                for (const ExtrusionEntity *entity : entities.entities) {
                    length += entity->length();
                    if (entity->role() == erExternalPerimeter || entity->inset_idx == 0)
                        for (const Polyline &line : entity->as_polylines())
                            outer.push_back(line.points);
                }
                const auto fills = region->fills.flatten();
                for (const ExtrusionEntity *entity : fills.entities)
                    if (entity->role() == erArcBridge || entity->role() == erArcOverhang)
                        last_arc_length += entity->length();
            }
        return std::make_pair(length, outer);
    };
    const auto normal = paths_for(0., true);
    const double normal_arc_length = last_arc_length;
    const auto overlap = paths_for(20., true);
    if (arc_threshold == 0. || arc_threshold == 5.) {
        CHECK(normal_arc_length > 0.);
        CHECK(last_arc_length > 0.);
    }
    REQUIRE_FALSE(normal.second.empty());
    if (arc_threshold == 0.)
        CHECK_THAT(overlap.first, Catch::Matchers::WithinAbs(normal.first, 1e-6));
    else
        CHECK(overlap.first > normal.first);
    CHECK(overlap.second == normal.second);
    const auto disabled = paths_for(0., false);
    const auto disabled_overlap = paths_for(20., false);
    CHECK_THAT(disabled_overlap.first, Catch::Matchers::WithinAbs(disabled.first, 1e-6));
    CHECK(disabled_overlap.second == disabled.second);
}

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

TEST_CASE("Shrinkage reinforcement keeps the outer bead envelope unchanged",
          "[Perimeters][LocalizedShrinkage][Regression]")
{
    const char *generator = GENERATE("classic", "arachne");
    CAPTURE(generator);
    const auto envelopes = [](const Print &print) {
        std::vector<std::array<double, 4>> result;
        for (const Layer *layer : print.objects().front()->layers()) {
            std::array<double, 4> bounds{1e30, 1e30, -1e30, -1e30};
            const auto add = [&](const ExtrusionPath &path) {
                if (path.role() != erExternalPerimeter)
                    return;
                for (const Point3 &point : path.polyline.points) {
                    const double x = unscale<double>(point.x()), y = unscale<double>(point.y());
                    bounds[0] = std::min(bounds[0], x - path.width / 2.);
                    bounds[1] = std::min(bounds[1], y - path.width / 2.);
                    bounds[2] = std::max(bounds[2], x + path.width / 2.);
                    bounds[3] = std::max(bounds[3], y + path.width / 2.);
                }
            };
            const auto visit = [&](auto &&self, const ExtrusionEntity &entity) -> void {
                if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity))
                    add(*path);
                else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity))
                    for (const auto &part : loop->paths) add(part);
                else if (const auto *multi = dynamic_cast<const ExtrusionMultiPath *>(&entity))
                    for (const auto &part : multi->paths) add(part);
                else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity))
                    for (const auto *child : collection->entities) self(self, *child);
            };
            for (const LayerRegion *region : layer->regions())
                visit(visit, region->perimeters);
            REQUIRE(bounds[0] < bounds[2]);
            result.push_back(bounds);
        }
        return result;
    };
    DynamicPrintConfig config = base_config(generator);
    config.set_deserialize_strict({{"localized_shrinkage_strategy", "disabled"},
        {"hull_line_extra_perimeters", 0}, {"hull_line_perimeter_expansion", 0.},
        {"inner_walls_flow_ratio", "110%"}});
    Print baseline;
    init_and_process_print({storage_box_transition()}, baseline, config);
    const auto expected = envelopes(baseline);
    config.set_deserialize_strict("localized_shrinkage_strategy", "reinforced_walls");
    Print reinforced;
    init_and_process_print({storage_box_transition()}, reinforced, config);
    const auto actual = envelopes(reinforced);
    REQUIRE(actual.size() == expected.size());
    CHECK(perimeter_length_at(reinforced, ledge_z) > perimeter_length_at(baseline, ledge_z));
    for (size_t layer = 0; layer < actual.size(); ++layer)
        for (size_t side = 0; side < 4; ++side) {
            CAPTURE(layer, side);
            CHECK_THAT(actual[layer][side], Catch::Matchers::WithinAbs(expected[layer][side], 0.002));
        }
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
    std::set<size_t> marked_layers;
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
        for (const LayerRegion *region : layer->regions()) {
            const size_t before = marked_paths;
            count_marked(count_marked, region->perimeters);
            if (marked_paths != before)
                marked_layers.insert(layer->id());
        }
    CHECK(marked_paths > 0);
    // Mitigation is a short structural band, not a one-layer stiffness step.
    CHECK(marked_layers.size() >= 3);
}

TEST_CASE("Compensated brick paths retain the brick preview classification",
          "[Perimeters][StaggeredPerimeters][Preview][Regression]")
{
    auto config = base_config(GENERATE("classic", "arachne"));
    config.set_deserialize_strict({{"perimeter_layering", "brick"},
        {"wall_loops", 3}, {"seam_gap", 0.}, {"seam_start_on_inner_wall", false},
        {"enable_arc_fitting", false}, {"gcode_comments", true}});
    Print print;
    init_and_process_print({make_cube(12., 12., 2.)}, print, config);
    double expected = 0.;
    const auto mark = [&](auto &&self, ExtrusionEntity &entity) -> void {
        if (auto *path = dynamic_cast<ExtrusionPath *>(&entity)) {
            if (path->staggered_perimeter) {
                expected += unscale<double>(path->length());
                path->shrinkage_compensation = true;
            }
        } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(&entity)) {
            for (auto &path : loop->paths) self(self, path);
        } else if (auto *multi = dynamic_cast<ExtrusionMultiPath *>(&entity)) {
            for (auto &path : multi->paths) self(self, path);
        } else if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
            for (auto *child : collection->entities) self(self, *child);
        }
    };
    for (Layer *layer : print.objects().front()->layers())
        for (LayerRegion *region : layer->regions()) mark(mark, region->perimeters);
    REQUIRE(expected > 0.);
    double actual = 0.;
    bool brick = false;
    GCodeReader reader;
    reader.parse_buffer(gcode(print), [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
        const std::string comment(line.comment());
        if (comment.find("TYPE:") != std::string::npos)
            brick = comment.find("Brick wall") != std::string::npos;
        if (brick && line.extruding(self)) actual += line.dist_XY(self);
    });
    CHECK_THAT(actual, Catch::Matchers::WithinAbs(expected, 0.05));
}

TEST_CASE("Compensated brick walls finish exposed tops at nominal height", "[Perimeters][LocalizedShrinkage][StaggeredPerimeters][Regression]")
{
    const char *generator = GENERATE("classic", "arachne");
    const bool inner_only = GENERATE(false, true);
    const char *strategy = GENERATE("balanced", "reinforced_walls");
    const char *surface_mode = GENERATE("disabled", "nonplanar_with_z_contouring_fallback");
    auto config = base_config(generator);
    config.set_deserialize_strict({{"localized_shrinkage_strategy", strategy},
        {"top_surface_z_mode", surface_mode},
        {"nonplanar_top_surface", std::string(surface_mode) != "disabled"},
        {"zaa_enabled", std::string(surface_mode) != "disabled"},
        {"perimeter_layering", "brick"}, {"staggered_perimeter_offset", "50%"},
        {"staggered_perimeters_inner_only", inner_only}, {"hull_line_perimeter_expansion", 0.15}});
    Print print;
    TriangleMesh shape = make_cube(40., 40., 5.);
    TriangleMesh boss = make_cube(10., 10., 5.);
    boss.translate(15., 15., 5.);
    shape.merge(boss);
    init_and_process_print({shape}, print, config);
    const auto paths_in = [](const Layer *layer) {
        std::vector<const ExtrusionPath *> paths;
        const auto collect = [&](auto &&self, const ExtrusionEntity &entity) -> void {
            if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity)) paths.push_back(path);
            else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity))
                for (const auto &p : loop->paths) paths.push_back(&p);
            else if (const auto *multi = dynamic_cast<const ExtrusionMultiPath *>(&entity))
                for (const auto &p : multi->paths) paths.push_back(&p);
            else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity))
                for (const auto *child : collection->entities) self(self, *child);
        };
        for (const LayerRegion *region : layer->regions()) collect(collect, region->perimeters);
        return paths;
    };
    size_t compensated = 0, brick = 0, restored = 0;
    for (const Layer *layer : print.objects().front()->layers()) {
        for (const ExtrusionPath *path : paths_in(layer)) {
            CAPTURE(generator, inner_only, strategy, surface_mode, layer->print_z, path->inset_idx);
            // This fixture's outer wall is vertical throughout the base, so
            // intermediate raised courses must not repeatedly get entry flow.
            if (layer->print_z > 1.0 && layer->print_z < 4.0 &&
                path->inset_idx == 0 && path->staggered_perimeter)
                CHECK_THAT(path->height, Catch::Matchers::WithinAbs(0.2, 1e-5));
            compensated += path->shrinkage_compensation;
            brick += path->staggered_perimeter;
            if (layer->upper_layer == nullptr) {
                CHECK_FALSE(path->staggered_perimeter);
                for (const Point3 &point : path->polyline.points) CHECK(point.z() == 0);
            }
            if (std::abs(layer->print_z - 5.0) > 1e-5 || path->staggered_perimeter ||
                !layer->lower_layer || path->polyline.points.size() < 2)
                continue;
            const Point3 &sample = path->polyline.points[path->polyline.points.size() / 2];
            for (const ExtrusionPath *lower : paths_in(layer->lower_layer)) {
                if (!lower->staggered_perimeter || lower->inset_idx != path->inset_idx)
                    continue;
                const auto &points = lower->polyline.points;
                if (std::any_of(points.begin(), points.end(), [&](const Point3 &p) {
                        return std::abs(unscale<double>(p.z()) - 0.1) < 1e-5;
                    })) {
                    ++restored;
                    CHECK_THAT(path->height, Catch::Matchers::WithinAbs(0.1, 1e-5));
                    CHECK(sample.z() == 0);
                    break;
                }
            }
        }
    }
    CHECK(compensated > 0);
    CHECK(brick > 0);
    CHECK(restored > 0);
}

TEST_CASE("Brick wall caps remain coplanar under hybrid surface fallback",
          "[Perimeters][LocalizedShrinkage][StaggeredPerimeters][ContourZ][Regression]")
{
    const char *generator = GENERATE("classic", "arachne");
    const bool brick = GENERATE(false, true);
    auto config = base_config(generator);
    config.set_deserialize_strict({{"localized_shrinkage_strategy", "reinforced_walls"},
        {"perimeter_layering", brick ? "brick" : "standard"},
        {"staggered_perimeters_inner_only", true}, {"staggered_perimeter_offset", "50%"},
        {"top_surface_z_mode", "nonplanar_with_z_contouring_fallback"},
        {"nonplanar_top_surface", true}, {"zaa_enabled", true}});
    // A hollow bucket with a rim between nominal Z planes. Integral-height
    // cubes miss the conflict: fallback lowers ordinary paths but skips the
    // already Z-tagged half-height brick return course.
    TriangleMesh bucket = make_cube(16., 16., 1.);
    for (const auto &wall : std::vector<std::array<double, 4>>{
             {0., 0., 16., 3.}, {0., 13., 16., 3.},
             {0., 3., 3., 10.}, {13., 3., 3., 10.}}) {
        TriangleMesh part = make_cube(wall[2], wall[3], 3.956);
        part.translate(wall[0], wall[1], 1.);
        bucket.merge(part);
    }
    Print print;
    init_and_process_print({bucket}, print, config);
    const Layer *top = print.objects().front()->layers().back();
    size_t paths = 0, returns = 0, adjusted = 0;
    const auto check = [&](auto &&self, const ExtrusionEntity &entity) -> void {
        if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity)) {
            ++paths;
            returns += path->staggered_transition;
            for (const Point3 &p : path->polyline.points) {
                adjusted += p.z() != 0;
                if (brick) CHECK(p.z() == 0);
            }
        } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
            for (const auto &path : loop->paths) self(self, path);
        } else if (const auto *multi = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
            for (const auto &path : multi->paths) self(self, path);
        } else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
            for (const auto *child : collection->entities) self(self, *child);
        }
    };
    CAPTURE(generator, brick, top->print_z);
    for (const LayerRegion *region : top->regions()) {
        check(check, region->perimeters);
        check(check, region->fills);
    }
    REQUIRE(paths > 0);
    if (brick) CHECK(returns > 0);
    else CHECK(adjusted > 0); // Control: ordinary fallback still follows this mesh.

    // Check the emitted coordinates too: seam handling and Z-state restoration
    // must not reintroduce a lowered skin beside the nominal return course.
    const double cap_z = top->print_z;
    const std::string output = gcode(print);
    double layer_z = 0.;
    size_t cap_moves = 0;
    GCodeReader reader;
    reader.apply_config(config);
    reader.parse_buffer(output, [&](GCodeReader &r, const GCodeReader::GCodeLine &line) {
        if (line.raw().rfind(";Z:", 0) == 0)
            layer_z = std::stod(line.raw().substr(3));
        if (std::abs(layer_z - cap_z) < 1e-5 && line.extruding(r) &&
            (line.has_x() || line.has_y())) {
            ++cap_moves;
            if (brick) CHECK_THAT(line.new_Z(r), Catch::Matchers::WithinAbs(cap_z, 1e-5));
        }
    });
    CHECK(cap_moves > 0);
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

namespace {

// The layer that closes the cavity of box_over_cavity(), the first one printed over air.
const double cavity_ceiling_z = 6.2;

// A cone standing on its tip, flaring by 5mm of radius per mm of height: at a layer height of 0.2 every
// wall of a layer lands a full millimetre outside the one below, entirely off the layer below but right
// alongside the walls printed with it.
TriangleMesh flared_cone()
{
    TriangleMesh cone = make_cone(20., 4.);
    cone.mirror(Z);
    cone.translate(0., 0., 4.);
    return cone;
}

// A 30mm box holding a 20mm cavity from z=2 to z=6, with a 4mm hole punched down through the ceiling
// of that cavity. The layer at cavity_ceiling_z bridges the cavity, and the walls of the hole sit in
// the middle of that bridge, 15mm clear of anything the layer below supports.
Print &box_over_cavity(Print &print, Model &model, const DynamicPrintConfig &config)
{
    ModelObject *object = model.add_object();
    object->name = "box_over_cavity.stl";
    object->add_volume(make_cube(30., 30., 8.), ModelVolumeType::MODEL_PART, false);
    TriangleMesh cavity = make_cube(20., 20., 4.);
    cavity.translate(5.f, 5.f, 2.f);
    object->add_volume(std::move(cavity), ModelVolumeType::NEGATIVE_VOLUME, false);
    TriangleMesh hole = make_cube(4., 4., 6.);
    hole.translate(13.f, 13.f, 5.f);
    object->add_volume(std::move(hole), ModelVolumeType::NEGATIVE_VOLUME, false);
    object->add_instance();
    object->ensure_on_bed();

    print.auto_assign_extruders(object);
    print.apply(model, config);
    print.validate();
    print.set_status_silent();
    return print;
}

// Every setting the assertions below depend on, so none of them rests on a default.
DynamicPrintConfig unsupported_walls_config(const char *wall_generator, bool unsupported_wall_last)
{
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({
        { "wall_generator",             wall_generator },
        { "layer_height",               0.2 },
        { "initial_layer_print_height", 0.2 },
        { "wall_loops",                 3 },
        { "detect_overhang_wall",       true },
        // Outer wall first, so an unsupported loop only ends up last if the feature puts it there.
        { "wall_sequence",              "outer wall/inner wall" },
        { "is_infill_first",            false },
        { "sparse_infill_density",      "15%" },
        { "unsupported_wall_last",      unsupported_wall_last },
        { "gcode_comments",             true },
    });
    return config;
}

// A loop extruded entirely in mid air: every one of its paths is an overhang.
bool unsupported_loop(const ExtrusionEntity *entity)
{
    if (! entity->is_loop())
        return false;
    const ExtrusionPaths &paths = static_cast<const ExtrusionLoop *>(entity)->paths;
    return ! paths.empty() && std::all_of(paths.begin(), paths.end(),
                                          [](const ExtrusionPath &path) { return path.role() == erOverhangPerimeter; });
}

// The loops of every wall island of the print, island by island, in extrusion order.
std::vector<std::vector<const ExtrusionLoop*>> wall_islands(const Print &print)
{
    std::vector<std::vector<const ExtrusionLoop*>> islands;
    for (const Layer *layer : print.objects().front()->layers())
        for (const LayerRegion *region : layer->regions())
            for (const ExtrusionEntity *island : region->perimeters.entities) {
                std::vector<const ExtrusionLoop*> loops;
                for (const ExtrusionEntity *entity : static_cast<const ExtrusionEntityCollection*>(island)->entities)
                    if (entity->is_loop())
                        loops.push_back(static_cast<const ExtrusionLoop*>(entity));
                islands.push_back(std::move(loops));
            }
    return islands;
}

// Islands where a loop that is anchored is extruded after one that is not.
int islands_with_a_supported_loop_last(const Print &print)
{
    int count = 0;
    for (const std::vector<const ExtrusionLoop*> &loops : wall_islands(print)) {
        bool seen_unsupported = false;
        for (const ExtrusionLoop *loop : loops) {
            if (unsupported_loop(loop))
                seen_unsupported = true;
            else if (seen_unsupported) {
                ++ count;
                break;
            }
        }
    }
    return count;
}

// The unsupported loops of the print, and those of them held back for the infill.
std::vector<const ExtrusionLoop*> unsupported_loops(const Print &print, double print_z = -1.)
{
    std::vector<const ExtrusionLoop*> loops;
    for (const Layer *layer : print.objects().front()->layers()) {
        if (print_z >= 0. && std::abs(layer->print_z - print_z) > EPSILON)
            continue;
        for (const LayerRegion *region : layer->regions())
            for (const ExtrusionEntity *island : region->perimeters.entities)
                for (const ExtrusionEntity *entity : static_cast<const ExtrusionEntityCollection*>(island)->entities)
                    if (unsupported_loop(entity))
                        loops.push_back(static_cast<const ExtrusionLoop*>(entity));
    }
    return loops;
}

int loops_held_back_for_infill(const std::vector<const ExtrusionLoop*> &loops)
{
    return int(std::count_if(loops.begin(), loops.end(), [](const ExtrusionLoop *loop) { return loop->print_after_infill; }));
}

// The G-code emitted at `print_z`, so the order of one layer can be read on its own.
std::string layer_gcode(const std::string &gcode, double print_z)
{
    std::string out;
    GCodeReader reader;
    reader.parse_buffer(gcode, [&out, print_z](GCodeReader &self, const GCodeReader::GCodeLine &line) {
        if (std::abs(self.z() - print_z) < EPSILON)
            out += line.raw() + "\n";
    });
    return out;
}

} // namespace

// Whatever the wall order asks for, a loop with nothing under it cannot be extruded before the loops it
// leans on. The flared cone gives every layer an outer wall that lands completely off the one below, and
// the outer wall first sequence would otherwise put it down before any of them.
TEST_CASE("Unsupported wall loops are extruded after the walls that anchor them", "[Perimeters]")
{
    const char *wall_generator = GENERATE("classic", "arachne");
    CAPTURE(wall_generator);

    auto slice_cone = [wall_generator](bool unsupported_wall_last, Print &print) {
        init_and_process_print({ flared_cone() }, print, unsupported_walls_config(wall_generator, unsupported_wall_last));
        REQUIRE_FALSE(print.objects().empty());
    };

    Print on;
    slice_cone(true, on);
    // Without unsupported loops to reorder the rest of the test would pass on an empty print.
    REQUIRE(unsupported_loops(on).size() > 0);
    CHECK(islands_with_a_supported_loop_last(on) == 0);

    SECTION("the held back loops run innermost first") {
        for (const std::vector<const ExtrusionLoop*> &loops : wall_islands(on)) {
            int previous_inset = std::numeric_limits<int>::max();
            for (const ExtrusionLoop *loop : loops)
                if (unsupported_loop(loop)) {
                    CHECK(loop->inset_idx <= previous_inset);
                    previous_inset = loop->inset_idx;
                }
        }
    }

    SECTION("switched off, the configured wall order is left alone") {
        Print off;
        slice_cone(false, off);
        REQUIRE(unsupported_loops(off).size() == unsupported_loops(on).size());
        // Outer wall first puts the unsupported outer wall ahead of the walls behind it.
        CHECK(islands_with_a_supported_loop_last(off) > 0);
    }
}

// A loop the walls cannot reach is a different case: only the bridges of its own layer will ever hold it,
// so it has to wait for them - while a loop that runs alongside a wall keeps its place, because the
// bridges anchor on it instead.
TEST_CASE("A wall loop out of reach of the layer below waits for the infill", "[Perimeters]")
{
    const char *wall_generator = GENERATE("classic", "arachne");
    CAPTURE(wall_generator);

    Print print;
    Model model;
    box_over_cavity(print, model, unsupported_walls_config(wall_generator, true));
    print.process();

    const std::vector<const ExtrusionLoop*> hole_loops = unsupported_loops(print, cavity_ceiling_z);
    REQUIRE(hole_loops.size() > 0);
    CHECK(loops_held_back_for_infill(hole_loops) == int(hole_loops.size()));

    SECTION("a loop alongside a supported wall is not held back") {
        Print cone;
        init_and_process_print({ flared_cone() }, cone, unsupported_walls_config(wall_generator, true));
        const std::vector<const ExtrusionLoop*> loops = unsupported_loops(cone);
        REQUIRE(loops.size() > 0);
        CHECK(loops_held_back_for_infill(loops) == 0);
    }

    SECTION("switched off, no loop is held back") {
        Print off;
        Model off_model;
        box_over_cavity(off, off_model, unsupported_walls_config(wall_generator, false));
        off.process();
        const std::vector<const ExtrusionLoop*> loops = unsupported_loops(off, cavity_ceiling_z);
        REQUIRE(loops.size() == hole_loops.size());
        CHECK(loops_held_back_for_infill(loops) == 0);
    }
}

// The held back loops reach the G-code in a second pass, after the infill of their layer: on the layer
// that closes the cavity the walls of the hole are extruded once the bridge is down, so the layer emits
// perimeters, then infill, then the perimeters that were waiting for it.
TEST_CASE("Loops waiting for the infill are extruded after it", "[Perimeters]")
{
    const char *wall_generator = GENERATE("classic", "arachne");
    CAPTURE(wall_generator);

    auto ceiling_roles = [wall_generator](bool unsupported_wall_last) {
        Print print;
        Model model;
        box_over_cavity(print, model, unsupported_walls_config(wall_generator, unsupported_wall_last));
        const std::string layer = layer_gcode(gcode(print), cavity_ceiling_z);
        REQUIRE_FALSE(layer.empty());
        return role_sequence(layer, { "perimeter", "infill" });
    };

    CHECK(ceiling_roles(true)  == std::vector<std::string>{ "perimeter", "infill", "perimeter" });
    CHECK(ceiling_roles(false) == std::vector<std::string>{ "perimeter", "infill" });
}

namespace {

// The rib spans z=[0,5] and the slab z=[5,6], so this is the slab's first layer - the only one whose
// support comes from the rib rather than from the slab below it.
const double slab_first_layer_z = 5.2;

// Rib widths either side of what the wall generators can print. At a 0.4mm nozzle the classic generator
// builds nothing thinner than nozzle/3 = 0.133mm and Arachne drops anything below min_feature_size, 25%
// of the nozzle = 0.1mm. 0.08mm is under both thresholds, 0.3mm over both.
const double unprintable_rib = 0.08;
const double printable_rib   = 0.3;

// A 4x5mm anchor tower carrying a 20x5mm slab at z=[5,6], with a rib `rib_width` wide running the whole
// length of the slab beneath its y=0 edge; a `rib_width` of 0 leaves the rib out. Nothing else is under
// that edge, so whether the wall along it is an overhang rests entirely on the rib. Overhang detection
// grows the lower slices by half the nozzle diameter before it asks, which carries either rib past the
// 0.21mm from the slab edge to that wall - the unprintable one only fails to reach it once it is filtered
// out for being unprintable.
Print &slab_over_rib(Print &print, Model &model, double rib_width, const DynamicPrintConfig &config)
{
    ModelObject *object = model.add_object();
    object->name = "slab_over_rib.stl";
    object->add_volume(make_cube(4., 5., 6.), ModelVolumeType::MODEL_PART, false);
    if (rib_width > 0.) {
        TriangleMesh rib = make_cube(20., rib_width, 5.);
        rib.translate(4.f, 0.f, 0.f);
        object->add_volume(std::move(rib), ModelVolumeType::MODEL_PART, false);
    }
    TriangleMesh slab = make_cube(20., 5., 1.);
    slab.translate(4.f, 0.f, 5.f);
    object->add_volume(std::move(slab), ModelVolumeType::MODEL_PART, false);
    object->add_instance();
    object->ensure_on_bed();

    print.auto_assign_extruders(object);
    print.apply(model, config);
    print.validate();
    print.set_status_silent();
    return print;
}

// Every setting the assertions below depend on, so none of them rests on a default. The wall line widths
// are pinned because the rib widths above are chosen against the distance from the slab edge to its outer
// wall, and min_feature_size because it is one of the two thresholds under test.
DynamicPrintConfig printable_rib_config(const char *wall_generator, bool detect_thin_wall)
{
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({
        { "wall_generator",                wall_generator },
        { "layer_height",                  0.2 },  // puts a layer boundary exactly on the top of the rib
        { "initial_layer_print_height",    0.2 },
        { "nozzle_diameter",               "0.4" },
        { "outer_wall_line_width",         0.42 },
        { "inner_wall_line_width",         0.45 },
        { "wall_loops",                    2 },
        { "detect_overhang_wall",          true },
        { "detect_thin_wall",              detect_thin_wall },
        { "min_feature_size",              "25%" },
        { "raft_layers",                   0 },
        // Anything that adds, drops or reorders walls would move length between the roles being counted.
        { "extra_perimeters_on_overhangs", false },
        { "overhang_reverse",              false },
        { "only_one_wall_top",             false },
        { "only_one_wall_first_layer",     false },
        { "unsupported_wall_last",         false },
        { "sparse_infill_density",         "15%" },
    });
    return config;
}

// Length of every overhang perimeter path on the layer at `print_z`, loops and open extrusions alike.
double overhang_length_at(const Print &print, double print_z)
{
    double len = 0.;
    const auto add_entity = [&len](const ExtrusionEntity *entity, auto &&self) -> void {
        const auto add_paths = [&len](const ExtrusionPaths &paths) {
            for (const ExtrusionPath &path : paths)
                if (path.role() == erOverhangPerimeter)
                    len += path.length();
        };
        if (const auto *coll = dynamic_cast<const ExtrusionEntityCollection*>(entity)) {
            for (const ExtrusionEntity *child : coll->entities)
                self(child, self);
        } else if (const auto *loop = dynamic_cast<const ExtrusionLoop*>(entity)) {
            add_paths(loop->paths);
        } else if (const auto *multi = dynamic_cast<const ExtrusionMultiPath*>(entity)) {
            add_paths(multi->paths);
        } else if (const auto *path = dynamic_cast<const ExtrusionPath*>(entity)) {
            if (path->role() == erOverhangPerimeter)
                len += path->length();
        }
    };

    for (const Layer *layer : print.objects().front()->layers()) {
        if (std::abs(layer->print_z - print_z) > EPSILON)
            continue;
        for (const LayerRegion *region : layer->regions())
            add_entity(&region->perimeters, add_entity);
    }
    return len;
}

} // namespace

// A sliver the wall generator prints nothing for holds nothing up, so it cannot be what decides that the
// wall above it is not an overhang. The rib under the slab is the only thing that edge of the slab could
// rest on: below the threshold of the active generator the slab has to come out exactly as it does with
// no rib at all, and the last check is the control - a rib the generator does print anchors that wall,
// without which the first check would hold for want of any sensitivity to the rib.
TEST_CASE("A lower layer sliver too thin to print does not support the wall above it", "[Perimeters]")
{
    const char *wall_generator   = GENERATE("classic", "arachne");
    const bool  detect_thin_wall = GENERATE(true, false);
    CAPTURE(wall_generator, detect_thin_wall);

    auto overhang_for = [wall_generator, detect_thin_wall](double rib_width) {
        Print print;
        Model model;
        slab_over_rib(print, model, rib_width, printable_rib_config(wall_generator, detect_thin_wall));
        print.process();
        REQUIRE_FALSE(print.objects().empty());
        return overhang_length_at(print, slab_first_layer_z);
    };

    const double no_rib      = overhang_for(0.);
    const double unprintable = overhang_for(unprintable_rib);
    const double printable   = overhang_for(printable_rib);

    // Only where the slab meets the tower is it held up from below, so both of its 20mm walls overhang.
    REQUIRE(no_rib > scale_(30.));
    CHECK_THAT(unprintable, Catch::Matchers::WithinAbs(no_rib, scale_(1.)));
    // A rib that does get printed takes the 20mm outer wall running along it out of the overhangs.
    CHECK(printable < no_rib - scale_(15.));
}
