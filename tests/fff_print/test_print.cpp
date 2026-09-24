#ifdef WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>
#endif

#include <catch2/catch_all.hpp>

#include "libslic3r/libslic3r.h"
#include "libslic3r/Print.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/PresetBundle.hpp"
#include "libslic3r/Format/SlicedGCode.hpp"
#include "libslic3r/miniz_extension.hpp"

#include "test_helpers.hpp"
#include "test_utils.hpp"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <string_view>

using namespace Slic3r;

TEST_CASE("Dispatch material selection excludes unused palette slots and unprinted objects", "[Print][ToolMapping]")
{
    auto config = Test::multifilament_config(5, {
        {"enable_prime_tower", "0"}, {"skirts", "0"}, {"brim_type", "no_brim"}});
    Print print;
    Model model;
    Test::init_print({make_cube(4., 4., 0.4), make_cube(4., 4., 0.4)}, print, model, config);
    model.objects[0]->config.set("extruder", 5);
    model.objects[1]->config.set("extruder", 2);
    for (auto *instance : model.objects[1]->instances)
        instance->printable = false;
    print.apply(model, config);
    CHECK(print.extruders(true) == std::vector<unsigned int>{4});
}

TEST_CASE("Mapped jobs generate physical tool commands and startup used-tool masks", "[Print][ToolMapping]")
{
    const bool use_second_material = GENERATE(false, true);
    const double destination_nozzle = GENERATE(0.4, 0.6);
    PresetBundle bundle;
    auto &printer = bundle.printers.get_edited_preset().config;
    printer.set_key_value("host_type", new ConfigOptionEnum<PrintHostType>(htOctoPrint));
    printer.set_deserialize_strict("single_extruder_multi_material", "0");
    printer.set_num_extruders(8);
    printer.set_deserialize_strict("nozzle_diameter", "0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4");
    printer.option<ConfigOptionFloats>("nozzle_diameter")->values[7] = destination_nozzle;
    printer.set_deserialize_strict("machine_start_gcode",
        ";MAPPED_INITIAL={initial_tool}\n;MAPPED_USED={is_extruder_used[0]},{is_extruder_used[1]},{is_extruder_used[7]}\n"
        ";MAPPED_TEMPS={nozzle_temperature[0]},{nozzle_temperature[7]}\n");
    printer.set_deserialize_strict("machine_end_gcode", "");
    printer.set_deserialize_strict("change_filament_gcode", "T{next_extruder}\n");
    bundle.set_num_filaments(2);
    // Use distinct temperatures so a T-only string rewrite cannot pass.
    auto material = bundle.filaments.get_edited_preset().config;
    material.set_key_value("nozzle_temperature", new ConfigOptionInts{205});
    material.set_key_value("nozzle_temperature_initial_layer", new ConfigOptionInts{215});
    bundle.filaments.load_preset("", "Dispatch PLA A", material, false);
    material.set_key_value("nozzle_temperature", new ConfigOptionInts{225});
    material.set_key_value("nozzle_temperature_initial_layer", new ConfigOptionInts{235});
    bundle.filaments.load_preset("", "Dispatch PLA B", material, false);
    bundle.filament_presets = {"Dispatch PLA A", "Dispatch PLA B"};
    auto &process = bundle.prints.get_edited_preset().config;
    process.set_deserialize_strict("enable_prime_tower", "0");
    process.set_deserialize_strict("skirts", "0");
    process.set_deserialize_strict("brim_type", "no_brim");
    process.set_deserialize_strict("layer_height", "0.2");
    process.set_deserialize_strict("initial_layer_print_height", "0.2");
    Model project;
    Print setup;
    Test::init_print({make_cube(4., 4., 0.4), make_cube(4., 4., 0.4)}, setup, project, bundle.full_config(false));
    project.objects[0]->config.set("extruder", 1);
    project.objects[1]->config.set("extruder", use_second_material ? 2 : 1);
    setup.apply(project, bundle.full_config(false));
    const std::vector<unsigned int> logical_tools = use_second_material ?
        std::vector<unsigned int>{0, 1} : std::vector<unsigned int>{0};
    CHECK(setup.extruders(true) == logical_tools);
    Model job(project);
    const auto config = bundle.tool_mapped_config(job, {7, 0}, {});
    CHECK_THAT(config.opt_float("nozzle_diameter", 7u), Catch::Matchers::WithinAbs(destination_nozzle, 1e-9));
    CHECK_THAT(config.opt_float("nozzle_diameter", 0u), Catch::Matchers::WithinAbs(0.4, 1e-9));
    CHECK(config.opt_int("nozzle_temperature_initial_layer", 0u) == 235);
    CHECK(config.opt_int("nozzle_temperature_initial_layer", 7u) == 215);
    Print print;
    print.set_status_silent();
    print.apply(job, config);
    const std::vector<unsigned int> physical_tools = use_second_material ?
        std::vector<unsigned int>{0, 7} : std::vector<unsigned int>{7};
    CHECK(print.extruders(true) == physical_tools);
    const std::string gcode = Test::gcode(print);
    const auto used_mask = gcode.find(";MAPPED_USED=");
    REQUIRE(used_mask != std::string::npos);
    INFO(gcode.substr(used_mask, 100));
    const std::string mask = use_second_material ? ";MAPPED_USED=true,false,true" : ";MAPPED_USED=false,false,true";
    CHECK(gcode.find(mask) != std::string::npos);
    CHECK(gcode.find(";MAPPED_TEMPS=225,205") != std::string::npos);
    CHECK(gcode.find("\nT7\n") != std::string::npos);
    CHECK((gcode.find("\nT0\n") != std::string::npos) == use_second_material);
    CHECK(gcode.find("\nT1\n") == std::string::npos);
    CHECK(project.objects[0]->config.extruder() == 1);
    CHECK(project.objects[1]->config.extruder() == (use_second_material ? 2 : 1));
    CHECK(bundle.filament_presets.size() == 2);
    ScopedTemporaryFile raw(".gcode"), package(".gcode.3mf");
    GCodeProcessorResult result;
    print.export_gcode(raw.string(), &result);
    REQUIRE(store_sliced_gcode_3mf(package.string(), print, result, raw.string()));
    mz_zip_archive archive{};
    REQUIRE(open_zip_reader(&archive, package.string()));
    struct CloseArchive {
        mz_zip_archive *archive;
        ~CloseArchive() { close_zip_reader(archive); }
    } guard{&archive};
    const auto entry = [&](const char *name) {
        size_t size = 0;
        void *data = mz_zip_reader_extract_file_to_heap(&archive, name, &size, 0);
        REQUIRE(data != nullptr);
        std::string content(static_cast<const char *>(data), size);
        mz_free(data);
        return content;
    };
    CHECK(entry("Metadata/plate_1.gcode").find(mask) != std::string::npos);
    const auto metadata = entry("Metadata/slice_info.config");
    CHECK(metadata.find("<filament id=\"8\"") != std::string::npos);
    CHECK((metadata.find("<filament id=\"1\"") != std::string::npos) == use_second_material);
    CHECK(metadata.find("<filament id=\"2\"") == std::string::npos);
}

TEST_CASE("Same or dissimilar materials bound to different features use their mapped nozzle geometry", "[Print][ToolMapping][ProjectFilamentBindings]")
{
    const bool with_support = GENERATE(false, true);
    PresetBundle bundle;
    auto &printer = bundle.printers.get_edited_preset().config;
    printer.set_num_extruders(3);
    printer.set_deserialize_strict("single_extruder_multi_material", "0");
    printer.set_deserialize_strict("nozzle_diameter", "0.6,0.2,0.4");
    printer.set_deserialize_strict("change_filament_gcode", "T{next_extruder}\n");
    printer.set_deserialize_strict("layer_change_gcode", "G92 E0\n");
    bundle.set_num_filaments(3); // The same profile in all three logical slots.
    if (GENERATE(false, true)) {
        const std::vector<std::string> materials{"PLA", "PETG", "TPU"};
        for (size_t i = 0; i < materials.size(); ++i) {
            auto filament = bundle.filaments.default_preset().config;
            filament.set_key_value("filament_type", new ConfigOptionStrings{materials[i]});
            const auto name = "mapped " + materials[i];
            bundle.filaments.load_preset("", name, filament, false);
            bundle.filament_presets[i] = name;
        }
    }
    auto &process = bundle.prints.get_edited_preset().config;
    process.set_deserialize_strict({{"enable_prime_tower", "0"}, {"skirts", "0"}, {"brim_type", "no_brim"},
        {"layer_height", "0.1"}, {"initial_layer_print_height", "0.1"}, {"wall_generator", "classic"},
        {"wall_loops", "2"}, {"sparse_infill_density", "30%"}, {"top_shell_layers", "2"}, {"bottom_shell_layers", "2"},
        {"line_width", "0"}, {"initial_layer_line_width", "0"}, {"outer_wall_line_width", "0"},
        {"inner_wall_line_width", "0"}, {"sparse_infill_line_width", "0"}, {"support_line_width", "0"}});
    for (const auto &key : project_filament_role_keys())
        process.set_key_value(key, new ConfigOptionInt(-1));
    process.set_key_value("raft_layers", new ConfigOptionInt(with_support ? 3 : 0));
    bundle.project_config.set_key_value("project_filament_bindings", new ConfigOptionInts{3, 2, 1, 2, 3, 3, 1, 1});
    Model project;
    Print setup;
    Test::init_print({make_cube(10., 10., 1.2)}, setup, project, bundle.full_config(false));
    Model job(project);
    const auto config = bundle.tool_mapped_config(job, {1, 2, 0}, {});
    Print print;
    print.set_status_silent();
    print.apply(job, config);
    const auto validation = print.validate();
    INFO(validation.string);
    REQUIRE(validation.string.empty());
    print.process();
    REQUIRE_FALSE(print.objects().empty());
    const auto *object = print.objects().front();
    REQUIRE_FALSE(object->layers().empty());
    REQUIRE_FALSE(object->layers().front()->regions().empty());
    const auto &region = object->layers().front()->regions().front()->region();
    const auto outer = region.flow(*object, frExternalPerimeter, 0.1, false);
    const auto inner = region.flow(*object, frPerimeter, 0.1, false);
    const auto infill = region.flow(*object, frInfill, 0.1, false);
    CHECK(region.extruder(frExternalPerimeter) == 2);
    CHECK(region.extruder(frPerimeter) == 3);
    CHECK(region.extruder(frInfill) == 1);
    CHECK_THAT(outer.nozzle_diameter(), Catch::Matchers::WithinAbs(0.2, 1e-6));
    CHECK_THAT(inner.nozzle_diameter(), Catch::Matchers::WithinAbs(0.4, 1e-6));
    CHECK_THAT(infill.nozzle_diameter(), Catch::Matchers::WithinAbs(0.6, 1e-6));
    CHECK(outer.width() < inner.width());
    CHECK(inner.width() < infill.width());
    CHECK(outer.mm3_per_mm() < inner.mm3_per_mm());
    CHECK(inner.mm3_per_mm() < infill.mm3_per_mm());
    CHECK_THAT(support_material_flow(object, 0.1f).nozzle_diameter(), Catch::Matchers::WithinAbs(0.6, 1e-6));
    CHECK_THAT(support_material_interface_flow(object, 0.1f).nozzle_diameter(), Catch::Matchers::WithinAbs(0.4, 1e-6));

    ScopedTemporaryFile raw(".gcode");
    GCodeProcessorResult result;
    print.export_gcode(raw.string(), &result);
    std::set<ExtrusionRole> found;
    for (const auto &move : result.moves) {
        if (move.type != EMoveType::Extrude) continue;
        if (move.extrusion_role == erExternalPerimeter) {
            CHECK(int(move.extruder_id) == 1);
            CHECK_THAT(move.width, Catch::Matchers::WithinAbs(outer.width(), 0.01));
            found.insert(move.extrusion_role);
        } else if (move.extrusion_role == erPerimeter) {
            CHECK(int(move.extruder_id) == 2);
            CHECK_THAT(move.width, Catch::Matchers::WithinAbs(inner.width(), 0.01));
            found.insert(move.extrusion_role);
        } else if (move.extrusion_role == erInternalInfill) {
            CHECK(int(move.extruder_id) == 0);
            CHECK_THAT(move.width, Catch::Matchers::WithinAbs(infill.width(), 0.01));
            found.insert(move.extrusion_role);
        } else if (move.extrusion_role == erSupportMaterial) {
            CHECK(int(move.extruder_id) == 0);
            found.insert(move.extrusion_role);
        } else if (move.extrusion_role == erSupportMaterialInterface) {
            CHECK(int(move.extruder_id) == 2);
            found.insert(move.extrusion_role);
        }
    }
    CHECK(found.size() == (with_support ? 5 : 3));
    CHECK(bundle.filament_presets.size() == 3);
    CHECK(process.opt_int("outer_wall_filament_id") == -1);
}

TEST_CASE("Per-print area slots reach sliced wall and infill regions",
          "[Print][ProjectFilamentBindings][Regression]")
{
    DynamicPrintConfig config = Test::multifilament_config(3, {
        {"outer_wall_filament_id", -1}, {"inner_wall_filament_id", -1},
        {"sparse_infill_filament_id", -1}, {"top_surface_filament_id", -1},
        {"support_filament", -1}, {"support_interface_filament", -1}});
    config.set_key_value("project_filament_bindings", new ConfigOptionInts{2, 3, 3, 2, 1, 0, 2, 0});
    Print print;
    Test::init_and_process_print({make_cube(5., 5., 1.)}, print, config);
    REQUIRE_FALSE(print.objects().empty());
    const PrintObject *object = print.objects().front();
    CHECK(object->config().support_filament.value == 2);
    CHECK(object->config().support_interface_filament.value == 3);
    REQUIRE_FALSE(object->layers().empty());
    REQUIRE_FALSE(object->layers().front()->regions().empty());
    const auto &region = object->layers().front()->regions().front()->region().config();
    CHECK(region.outer_wall_filament_id.value == 3);
    CHECK(region.inner_wall_filament_id.value == 2);
    CHECK(region.sparse_infill_filament_id.value == 1);
    CHECK(region.top_surface_filament_id.value == 2);
    CHECK(config.opt_int("outer_wall_filament_id") == -1);
}
using namespace Slic3r::Test;

namespace {

std::vector<unsigned int> perimeter_inset_order(const LayerRegion &region)
{
    ExtrusionEntityCollection flattened = region.perimeters.flatten();
    std::vector<unsigned int> order;
    order.reserve(flattened.entities.size());
    for (const ExtrusionEntity *entity : flattened.entities)
        if (entity->is_loop())
            order.push_back(entity->inset_idx);
    return order;
}

} // namespace

TEST_CASE("Inner wall seam preparation preserves the selected wall order", "[Print][Seam]")
{
    const std::string wall_generator = GENERATE("classic", "arachne");
    const std::string wall_sequence = GENERATE("outer wall/inner wall", "inner wall/outer wall");
    CAPTURE(wall_generator);
    CAPTURE(wall_sequence);
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({
        {"wall_generator", wall_generator},
        {"wall_loops", 3},
        {"wall_sequence", wall_sequence},
        {"seam_start_on_inner_wall", true},
        {"layer_height", 0.2},
        {"initial_layer_print_height", 0.2}
    });

    Print print;
    init_and_process_print({cube(20.0)}, print, config);

    const Layer *layer = print.objects().front()->get_layer(1);
    REQUIRE(layer != nullptr);
    REQUIRE_FALSE(layer->regions().empty());
    const std::vector<unsigned int> order = perimeter_inset_order(*layer->regions().front());
    REQUIRE(order.size() >= 3);
    if (wall_sequence == "outer wall/inner wall") {
        CHECK(order[0] == 0);
        CHECK(order[1] == 1);
        CHECK(order[2] == 2);
    } else {
        CHECK(order[0] == 2);
        CHECK(order[1] == 1);
        CHECK(order[2] == 0);
    }
}

TEST_CASE("Disabled inner wall seam start preserves the selected wall order", "[Print][Seam]")
{
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({
        {"wall_generator", "classic"},
        {"wall_loops", 3},
        {"wall_sequence", "outer wall/inner wall"},
        {"seam_start_on_inner_wall", false},
        {"layer_height", 0.2},
        {"initial_layer_print_height", 0.2}
    });

    Print print;
    init_and_process_print({cube(20.0)}, print, config);

    const Layer *layer = print.objects().front()->get_layer(1);
    REQUIRE(layer != nullptr);
    REQUIRE_FALSE(layer->regions().empty());
    const std::vector<unsigned int> order = perimeter_inset_order(*layer->regions().front());
    REQUIRE(order.size() >= 3);
    CHECK(order[0] == 0);
}

TEST_CASE("Internal seam entries and returns preserve the original wall material", "[Print][Seam]")
{
    const std::string wall_generator = GENERATE("classic", "arachne");
    CAPTURE(wall_generator);
    const auto generate_gcode = [&](bool prepare_outer_seam) {
        return slice({make_cylinder(10., 1., PI / 60.)}, {
            {"wall_generator", wall_generator},
            {"wall_loops", 2},
            {"wall_sequence", "outer wall/inner wall"},
            {"seam_start_on_inner_wall", prepare_outer_seam},
            {"seam_gap", 0.0},
            {"gcode_comments", true},
            {"skirt_loops", 0},
            {"brim_type", "no_brim"},
            {"layer_height", 0.2},
            {"initial_layer_print_height", 0.2}
        });
    };
    const std::string generated_gcode = generate_gcode(true);

    size_t transitions = 0;
    size_t prime_segments = 0;
    double enabled_extrusion_length = 0.0;
    GCodeReader reader;
    reader.parse_buffer(generated_gcode, [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
        const std::string comment(line.comment());
        const bool connector = comment.find("outer wall seam transition") != std::string::npos ||
                               comment.find("outer wall seam return") != std::string::npos;
        if (line.extruding(self) && !connector)
            enabled_extrusion_length += line.dist_XY(self);
        if (line.dist_XY(self) > 0.0 &&
            comment.find("outer wall seam transition") != std::string::npos) {
            ++transitions;
            CHECK(line.extruding(self));
        }
        if (line.has_e() && line.dist_XY(self) > 0.0 &&
            comment.find("outer wall seam prime") != std::string::npos) {
            ++prime_segments;
            CHECK(line.extruding(self));
        }
    });
    CHECK(transitions > 0);
    CHECK(prime_segments > 0);

    double baseline_extrusion_length = 0.0;
    GCodeReader baseline_reader;
    baseline_reader.parse_buffer(
        generate_gcode(false),
        [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
            if (line.extruding(self))
                baseline_extrusion_length += line.dist_XY(self);
        });
    CHECK(enabled_extrusion_length == Catch::Approx(baseline_extrusion_length).epsilon(0.005));
}

TEST_CASE("Internal seam starts and finishes inside without retracing the reserved wall",
          "[Print][Seam][Regression]")
{
    const std::string generator = GENERATE("classic", "arachne");
    const std::string sequence = GENERATE("outer wall/inner wall", "inner wall/outer wall", "inner-outer-inner wall");
    const double gap = GENERATE(0., 0.06);
    const std::string shape = GENERATE("cube", "cylinder", "concave hole");
    const bool reverse_inner = GENERATE(false, true);
    const std::string layering = GENERATE("standard", "brick");
    CAPTURE(generator, sequence, gap, shape, layering, reverse_inner);
    TriangleMesh mesh = shape == "concave hole" ? Test::mesh(Test::TestMesh::cube_with_concave_hole) :
        shape == "cylinder" ? make_cylinder(6., 1., PI / 60.) : make_cube(12., 12., 1.);
    mesh.scale(Vec3f(1.f, 1.f, 1.f / mesh.bounding_box().size().z()));
    Print print;
    Test::init_and_process_print({mesh}, print, {
        {"wall_generator", generator}, {"wall_sequence", sequence}, {"wall_loops", 3},
        {"perimeter_layering", layering}, {"staggered_perimeters_inner_only", true},
        {"seam_start_on_inner_wall", true}, {"seam_slope_type", "none"}, {"seam_gap", gap},
        {"enable_arc_fitting", false}, {"gcode_comments", true}, {"skirt_loops", 0},
        {"brim_type", "no_brim"}, {"layer_height", 0.2}, {"initial_layer_print_height", 0.2},
        {"sparse_infill_density", "0%"}, {"top_shell_layers", 0}, {"bottom_shell_layers", 0}});
    if (reverse_inner) {
        std::function<void(ExtrusionEntity &)> reverse = [&](ExtrusionEntity &entity) {
            if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
                for (auto *child : collection->entities) reverse(*child);
            } else if (entity.inset_idx == 1) {
                entity.reverse();
            }
        };
        for (const Layer *layer : print.objects().front()->layers())
            for (LayerRegion *region : layer->regions())
                reverse(region->perimeters);
    }
    ScopedTemporaryFile file(".gcode");
    print.export_gcode(file.string(), nullptr, nullptr);
    struct Segment { Vec2d a, b; double z; std::string comment; };
    std::vector<Segment> segments;
    bool in_seam = false;
    Vec2d inside_start = Vec2d::Zero();
    size_t entries = 0, returns = 0, finishes = 0;
    GCodeReader reader;
    reader.parse_file(file.string(), [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
        if (line.extruding(self) && line.dist_XY(self) > 0.001)
            segments.push_back({Vec2d(self.x(), self.y()), Vec2d(line.new_X(self), line.new_Y(self)), self.z(), std::string(line.comment())});
        // Concave contours may deliberately fall back to an ordinary seam.
        // They must still never cross or duplicate another extrusion.
        if (shape == "concave hole")
            return;
        const std::string comment(line.comment());
        const bool prime = comment.find("outer wall seam prime") != std::string::npos;
        const bool finish = comment.find("outer wall seam finish") != std::string::npos;
        if (prime && !in_seam && line.extruding(self)) {
            inside_start = Vec2d(self.x(), self.y());
            in_seam = true;
        }
        if (in_seam && line.dist_XY(self) > 0.001)
            CHECK(line.extruding(self));
        if (in_seam && std::abs(line.new_Z(self) - self.z()) > 0.001)
            CHECK(line.extruding(self));
        if (in_seam && line.has_e()) CHECK(line.dist_E(self) >= 0.);
        if (comment.find("outer wall seam transition") != std::string::npos) {
            CHECK(line.extruding(self));
            ++entries;
        }
        if (comment.find("outer wall seam return") != std::string::npos) {
            CHECK(line.extruding(self));
            ++returns;
        }
        if (finish && line.extruding(self) &&
            (Vec2d(line.new_X(self), line.new_Y(self)) - inside_start).norm() < 0.002) {
            in_seam = false;
            ++finishes;
        }
    });
    CHECK_FALSE(in_seam);
    if (shape == "cylinder")
        REQUIRE(entries > 0);
    CHECK(returns == entries);
    CHECK(finishes == entries);
    // Include later inner walls, not just the entry itself: reserving a tail
    // is only correct if that wall never extrudes over it again.
    for (size_t i = 0; i < segments.size(); ++i) {
        const auto &a = segments[i];
        const Vec2d direction = (a.b - a.a).normalized();
        const auto cross = [&](const Vec2d &v) { return direction.x() * v.y() - direction.y() * v.x(); };
        for (size_t j = i + 1; j < segments.size(); ++j) {
            const auto &b = segments[j];
            if (std::abs(a.z - b.z) > 0.001)
                continue;
            const Vec2d other = b.b - b.a;
            const double determinant = cross(other);
            if (std::abs(determinant) > 1e-6) {
                const double u = -cross(b.a - a.a) / determinant;
                const double distance = (b.a + u * other - a.a).dot(direction);
                const bool interior_crossing =
                    u * other.norm() > 0.005 && (1. - u) * other.norm() > 0.005 &&
                    distance > 0.005 && distance < (a.b - a.a).norm() - 0.005;
                CAPTURE(i, j, distance, u, a.comment, b.comment, a.a.x(), a.a.y(), a.b.x(), a.b.y(), b.a.x(), b.a.y(), b.b.x(), b.b.y());
                CHECK_FALSE(interior_crossing);
                continue;
            }
            if (std::abs(cross(b.a - a.a)) > 0.001 || std::abs(cross(b.b - a.a)) > 0.001)
                continue;
            const double p = (b.a - a.a).dot(direction), q = (b.b - a.a).dot(direction);
            const double overlap = std::min((a.b - a.a).norm(), std::max(p, q)) - std::max(0., std::min(p, q));
            CAPTURE(i, j, overlap);
            CHECK(overlap < 0.005);
        }
    }
}

TEST_CASE("Opposite winding inner walls keep connected seam entries and returns", "[Print][Seam][ReverseTail]")
{
    const std::string generator = GENERATE("classic", "arachne");
    const std::string layering = GENERATE("standard", "brick");
    Print print;
    Test::init_and_process_print({make_cylinder(6., 1., PI / 60.)}, print, {
        {"wall_generator", generator}, {"wall_sequence", "outer wall/inner wall"},
        {"wall_loops", 3}, {"perimeter_layering", layering},
        {"staggered_perimeters_inner_only", true}, {"seam_start_on_inner_wall", true},
        {"seam_slope_type", "none"}, {"seam_gap", 0.06}, {"seam_position", "back"},
        {"enable_arc_fitting", false}, {"gcode_comments", true}, {"skirt_loops", 0},
        {"brim_type", "no_brim"}, {"layer_height", 0.2}, {"initial_layer_print_height", 0.2},
        {"sparse_infill_density", "0%"}, {"top_shell_layers", 0}, {"bottom_shell_layers", 0}
    });
    // Opposite winding changes the available tail direction, not the wall
    // footprint. The seam planner must qualify both possible deposit orders.
    std::function<void(ExtrusionEntity &)> reverse_inner = [&](ExtrusionEntity &entity) {
        if (auto *collection = dynamic_cast<ExtrusionEntityCollection *>(&entity)) {
            for (auto *child : collection->entities) reverse_inner(*child);
        } else if (entity.inset_idx == 1) {
            entity.reverse();
        }
    };
    for (const Layer *layer : print.objects().front()->layers())
        for (LayerRegion *region : layer->regions())
            reverse_inner(region->perimeters);
    ScopedTemporaryFile file(".gcode");
    print.export_gcode(file.string(), nullptr, nullptr);
    size_t entries = 0, returns = 0;
    GCodeReader reader;
    reader.parse_file(file.string(), [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
        const auto comment = line.comment();
        // Exclude approach travels, whose comments also contain the seam name.
        if (comment == " outer wall seam transition" || comment == " outer wall seam return") {
            CHECK(line.extruding(self));
            if (comment == " outer wall seam transition") ++entries;
            else ++returns;
        }
    });
    CHECK(entries > 0);
    CHECK(returns == entries);
}

TEST_CASE("Inner wall seam preparation follows the configured seam position", "[Print][Seam]")
{
    const std::string seam_position =
        GENERATE("nearest", "aligned", "aligned_back", "back", "random");
    CAPTURE(seam_position);
    const std::string generated_gcode = slice({make_cylinder(10., 1., PI / 60.)}, {
        {"wall_generator", "classic"},
        {"wall_loops", 2},
        {"wall_sequence", "outer wall/inner wall"},
        {"seam_position", seam_position},
        {"seam_start_on_inner_wall", true},
        {"seam_gap", 0.0},
        {"gcode_comments", true},
        {"skirt_loops", 0},
        {"brim_type", "no_brim"},
        {"layer_height", 0.2},
        {"initial_layer_print_height", 0.2}
    });

    bool awaiting_outer_wall = false;
    Vec2f transition_end = Vec2f::Zero();
    size_t aligned_transitions = 0;
    GCodeReader reader;
    reader.parse_buffer(generated_gcode, [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
        const std::string comment(line.comment());
        if (line.dist_XY(self) > 0. &&
            comment.find("outer wall seam transition") != std::string::npos) {
            transition_end = Vec2f(line.new_X(self), line.new_Y(self));
            awaiting_outer_wall = true;
        } else if (awaiting_outer_wall && line.extruding(self)) {
            // The transition target is the outer-loop seam selected by
            // SeamPlacer, so every placement strategy must begin extrusion at
            // exactly that point.
            CHECK(self.x() == Catch::Approx(transition_end.x()).margin(0.0015));
            CHECK(self.y() == Catch::Approx(transition_end.y()).margin(0.0015));
            awaiting_outer_wall = false;
            ++aligned_transitions;
        }
    });
    CHECK(aligned_transitions > 0);
    CHECK_FALSE(awaiting_outer_wall);
}

TEST_CASE("Qualified scarf seams replace internal entries without removing inner wall material", "[Print][Seam][ScarfFallback][Regression]")
{
    const bool smooth = GENERATE(false, true);
    const bool scarf_inner = GENERATE(false, true);
    const std::string sequence = GENERATE("outer wall/inner wall", "inner wall/outer wall");
    const std::string generator = GENERATE("classic", "arachne");
    CAPTURE(smooth, scarf_inner, sequence, generator);
    const auto output = [&](bool internal) {
        return slice({smooth ? make_cylinder(8., 1., PI / 120.) : make_cube(16., 16., 1.)}, {
            {"wall_generator", generator}, {"wall_loops", 3},
            {"wall_sequence", sequence}, {"seam_start_on_inner_wall", internal},
            {"seam_slope_type", "all"}, {"seam_slope_conditional", true},
            {"scarf_angle_threshold", 155.}, {"scarf_overhang_threshold", 40.},
            {"seam_slope_inner_walls", scarf_inner}, {"seam_gap", 0.},
            {"layer_height", 0.2}, {"initial_layer_print_height", 0.2},
            {"enable_arc_fitting", false}, {"gcode_comments", true},
            {"skirt_loops", 0}, {"brim_type", "no_brim"}
        });
    };
    size_t entries = 0;
    size_t first_layer_entries = 0;
    size_t slopes = 0;
    const auto deposited = [&](const std::string &gcode, bool inspect) {
        std::pair<double, double> result{0., 0.};
        GCodeReader reader;
        reader.parse_buffer(gcode, [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
            const std::string comment(line.comment());
            const bool connector = comment.find("outer wall seam transition") != std::string::npos ||
                                   comment.find("outer wall seam return") != std::string::npos;
            if (line.extruding(self) && line.dist_XY(self) > 0. && !connector) {
                result.first += line.dist_XY(self);
                result.second += line.dist_E(self);
            }
            if (inspect && self.z() <= 0.21 &&
                comment.find("outer wall seam transition") != std::string::npos)
                ++first_layer_entries;
            // The first layer deliberately cannot scarf and must still use
            // the internal fallback. Check subsequent layers separately.
            if (inspect && self.z() > 0.21) {
                if (comment.find("outer wall seam transition") != std::string::npos)
                    ++entries;
                if (line.extruding(self) && line.dist_XY(self) > 0. &&
                    std::abs(line.new_Z(self) - self.z()) > 0.0001)
                    ++slopes;
            }
        });
        return result;
    };
    const auto baseline = deposited(output(false), false);
    const auto combined = deposited(output(true), true);
    if (smooth) CHECK(first_layer_entries > 0);
    CHECK_THAT(combined.first, Catch::Matchers::WithinAbs(baseline.first, 0.05));
    CHECK_THAT(combined.second, Catch::Matchers::WithinAbs(baseline.second, 0.01));
    if (smooth) {
        CHECK(slopes > 0);
        CHECK(entries == 0);
    } else {
        CHECK(slopes == 0);
        CHECK(entries == 0); // Sharp corners fall back when connector beads overlap.
    }
}

TEST_CASE("Internal seam preparation preserves wall material with variable height modes", "[Print][Seam][Regression]")
{
    const std::string mode = GENERATE("interlocking_walls", "nonplanar_top_surface", "z_contouring");
    const bool wall_mode = mode == "brick" || mode == "interlocking_walls";
    const std::string sequence = GENERATE("outer wall/inner wall", "inner wall/outer wall");
    CAPTURE(mode, sequence);
    const auto perimeter_length = [&](bool inner_seam) {
        const std::string output = slice({cube(5.0)}, {
            {"wall_loops", 3},
            {"wall_sequence", sequence},
            {"perimeter_layering", wall_mode ? mode : "standard"},
            {"top_surface_z_mode", wall_mode ? "disabled" : mode},
            {"seam_start_on_inner_wall", inner_seam},
            {"seam_gap", 0.0},
            {"seam_slope_type", "none"},
            {"gcode_comments", true},
            {"skirt_loops", 0},
            {"brim_type", "no_brim"},
            {"layer_height", 0.2},
            {"initial_layer_print_height", 0.2}
        });
        double length = 0.;
        GCodeReader reader;
        reader.parse_buffer(output, [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
            const std::string comment(line.comment());
            const bool original_wall = comment.find("perimeter") != std::string::npos ||
                comment.find("outer wall seam prime") != std::string::npos ||
                comment.find("inner wall seam restore") != std::string::npos ||
                comment.find("outer wall seam finish") != std::string::npos;
            if (line.extruding(self) && original_wall)
                length += std::hypot(line.dist_XY(self), line.new_Z(self) - self.z());
        });
        return length;
    };
    const double baseline = perimeter_length(false);
    REQUIRE(baseline > 0.);
    CHECK_THAT(perimeter_length(true), Catch::Matchers::WithinAbs(baseline, 0.05));
}

TEST_CASE("Brick courses preserve wall material when unsafe corner connectors fall back",
          "[Print][Seam][StaggeredPerimeters][Regression]")
{
    const std::string generator = GENERATE("classic", "arachne");
    const std::string sequence = GENERATE("outer wall/inner wall", "inner wall/outer wall");
    const bool inner_only = GENERATE(false, true);
    const std::string seam_position = GENERATE("aligned", "back");
    const std::string surface_mode = GENERATE("disabled", "nonplanar_with_z_contouring_fallback");
    CAPTURE(generator, sequence, inner_only, surface_mode, seam_position);
    const auto output = [&](bool internal) {
        return slice({make_cube(12., 12., 2.)}, {
            {"wall_generator", generator}, {"wall_sequence", sequence}, {"wall_loops", 3},
            {"perimeter_layering", "brick"}, {"staggered_perimeters_inner_only", inner_only},
            {"top_surface_z_mode", surface_mode},
            {"seam_position", seam_position},
            {"staggered_perimeter_offset", "50%"}, {"seam_start_on_inner_wall", internal},
            {"seam_slope_type", "none"}, {"seam_gap", 0.}, {"enable_arc_fitting", false},
            {"gcode_comments", true}, {"skirt_loops", 0}, {"brim_type", "no_brim"},
            {"layer_height", 0.2}, {"initial_layer_print_height", 0.2},
            {"sparse_infill_density", "0%"}, {"top_shell_layers", 0}, {"bottom_shell_layers", 0}});
    };
    size_t raised_entries = 0, returns = 0, entries = 0;
    const auto volume = [&](const std::string &text, bool inspect) {
        double deposited = 0.;
        bool brick = false, inside = false;
        double seam_z = 0., inner_z = 0.;
        GCodeReader reader;
        reader.parse_buffer(text, [&](GCodeReader &self, const GCodeReader::GCodeLine &line) {
            const std::string comment(line.comment());
            if (comment.find("TYPE:") != std::string::npos)
                brick = comment.find("Brick wall") != std::string::npos;
            const bool entry = comment.find("outer wall seam transition") != std::string::npos;
            const bool finish = comment.find("outer wall seam return") != std::string::npos;
            if (line.extruding(self) && line.dist_XY(self) > 0. && !entry && !finish)
                deposited += line.dist_E(self);
            if (!inspect) return;
            if (entry) {
                ++entries;
                if (brick) ++raised_entries;
                CHECK(line.extruding(self));
                inner_z = self.z();
                seam_z = line.new_Z(self);
                inside = true;
            }
            if (inside && (line.has_x() || line.has_y() || line.has_z())) {
                CHECK_THAT(double(line.new_Z(self)), Catch::Matchers::WithinAbs(finish ? inner_z : seam_z, 0.001));
                if (line.dist_XY(self) > 0.001) CHECK(line.extruding(self));
            }
            if (finish) { ++returns; inside = false; }
        });
        CHECK_FALSE(inside);
        return deposited;
    };
    const double expected = volume(output(false), false);
    CHECK_THAT(volume(output(true), true), Catch::Matchers::WithinAbs(expected, 0.015));
    CHECK(entries == 0);
    CHECK(returns == 0);
    CHECK(raised_entries == 0);
}

TEST_CASE("Timelapse g-code is emitted once per layer for Bambu and non-Bambu printers", "[Print][Regression]")
{
    struct PrinterCase {
        std::string name;
        std::string structure;
        bool        is_bbl;
    };
    const PrinterCase printer = GENERATE(from_range(std::vector<PrinterCase>{
        { "non-BBL undefined", "undefine", false },
        { "non-BBL CoreXY",    "corexy",   false },
        { "non-BBL i3",        "i3",       false },
        { "non-BBL H-Bot",     "hbot",     false },
        { "non-BBL Delta",     "delta",    false },
        { "Bambu CoreXY",      "corexy",   true },
        { "Bambu i3",          "i3",       true },
    }));
    INFO("printer: " << printer.name);

    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({
        { "initial_layer_print_height", 0.2 },
        { "layer_change_gcode",          ";TEST_LAYER_CHANGE" },
        { "layer_height",                0.2 },
        { "printer_structure",           printer.structure },
        { "spiral_mode",                 false },
        { "time_lapse_gcode",            "TIMELAPSE_TAKE_FRAME" },
    });
    Print print;
    print.is_BBL_printer() = printer.is_bbl;
    Model model;
    init_print({ cube(20) }, print, model, config);
    const std::string gcode = Slic3r::Test::gcode(print);

    const auto count = [&gcode](std::string_view token) {
        size_t occurrences = 0;
        size_t pos = 0;
        while ((pos = gcode.find(token, pos)) != std::string::npos) {
            ++occurrences;
            pos += token.size();
        }
        return occurrences;
    };

    const size_t layer_changes = count("\n;TEST_LAYER_CHANGE\n");
    REQUIRE(layer_changes > 0);
    CHECK(count("\nTIMELAPSE_TAKE_FRAME\n") == layer_changes);
}

SCENARIO("Changing the number of solid shell layers does not make all surfaces internal", "[Print]") {
    GIVEN("sliced 20mm cube and config with top_shell_layers = 2 and bottom_shell_layers = 1") {
        Slic3r::DynamicPrintConfig config = Slic3r::DynamicPrintConfig::full_print_config();
		config.set_deserialize_strict({
            { "top_shell_layers",           2 },
            { "bottom_shell_layers",        1 },
            { "layer_height",               0.25 }, // get a known number of layers
            { "initial_layer_print_height", 0.25 }
			});
        Slic3r::Print print;
        Slic3r::Model model;
        Slic3r::Test::init_print({cube(20)}, print, model, config);
        // Precondition: Ensure that the model has 2 solid top layers (79, 78)
        // and one solid bottom layer (0).
		auto test_is_solid_infill = [&print](size_t obj_id, size_t layer_id) {
		    const Layer &layer = *(print.objects().at(obj_id)->get_layer((int)layer_id));
		    // iterate over all of the regions in the layer
		    for (const LayerRegion *region : layer.regions()) {
		        // for each region, iterate over the fill surfaces
		        for (const Surface &surface : region->fill_surfaces.surfaces)
		            CHECK(surface.is_solid());
		    }
		};
        print.process();
        test_is_solid_infill(0,  0); // should be solid
        test_is_solid_infill(0, 79); // should be solid
        test_is_solid_infill(0, 78); // should be solid
        WHEN("Model is re-sliced with top_shell_layers == 3") {
			config.set("top_shell_layers", 3);
			print.apply(model, config);
            print.process();
            THEN("Print object does not have 0 solid bottom layers.") {
                test_is_solid_infill(0, 0);
            }
            AND_THEN("Print object has 3 top solid layers") {
                test_is_solid_infill(0, 79);
                test_is_solid_infill(0, 78);
                test_is_solid_infill(0, 77);
            }
        }
    }
}

// ---------------------------------------------------------------------------
TEST_CASE("Imported nonplanar settings are rejected before slicing", "[Print][Nonplanar]")
{
    const auto [key, value] = GENERATE(table<const char *, const char *>({
        {"top_surface_z_mode", "nonplanar_top_surface"},
        {"perimeter_layering", "interlocking_walls"},
        {"nonplanar_infill", "1"},
        {"support_ironing_nonplanar", "1"}
    }));
    CAPTURE(key);
    Model model;
    Print print;
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict(key, value);
    Slic3r::Test::init_print({make_cube(20., 20., 20.)}, print, model, config);
    const auto error = print.validate();
    CHECK(error.opt_key == key);
    CHECK_FALSE(error.string.empty());
}

TEST_CASE("Nonplanar quarantine leaves standard brick and Z contouring available", "[Print][Nonplanar]")
{
    const auto [key, value] = GENERATE(table<const char *, const char *>({
        {"perimeter_layering", "standard"},
        {"perimeter_layering", "brick"},
        {"top_surface_z_mode", "z_contouring"}
    }));
    CAPTURE(key, value);
    Model model;
    Print print;
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict(key, value);
    config.set_key_value("layer_change_gcode", new ConfigOptionString("G92 E0\n"));
    Slic3r::Test::init_print({make_cube(20., 20., 20.)}, print, model, config);
    CHECK(print.validate().string.empty());
}

// Print::validate() warning collection
//
// validate() returns its warnings in a vector. The warning paths deliberately
// differ in how many entries they produce; these tests pin down each behaviour:
//   * independent checks    -> stack (one entry each)
//   * motion-ability        -> coalesce into one (mutually exclusive, gated)
//   * clumping detection    -> one independent warning
//   * layered clearance     -> many collisions concatenated into one entry
//   * null warnings pointer -> no-op, no crash, no blocking error
// ---------------------------------------------------------------------------
namespace {

// Build `n` 20mm cubes (spread apart, or stacked at the origin when `overlap`) into
// `model`/`print` and apply `config`, leaving the print ready to validate(). No slicing needed.
void build_cubes(Slic3r::Model& model, Slic3r::Print& print,
                 DynamicPrintConfig config, int n, bool overlap)
{
    config.set_key_value("layer_change_gcode", new ConfigOptionString("G92 E0\n")); // validate() relative-E reset

    for (int i = 0; i < n; ++i) {
        ModelObject* object = model.add_object();
        object->add_volume(cube(20));
        ModelInstance* inst = object->add_instance();
        inst->set_offset(Vec3d(overlap ? 0.0 : i * 60.0, 0.0, 0.0));
    }
    for (ModelObject* mo : model.objects) {
        mo->ensure_on_bed();
        print.auto_assign_extruders(mo);
    }
    print.apply(model, config);
}

// Build cubes and run validate(), collecting warnings; returns the blocking error.
StringObjectException validate_cubes(const DynamicPrintConfig& config,
                                     std::vector<StringObjectException>& warnings,
                                     int n = 1, bool overlap = false)
{
    Slic3r::Model model;
    Slic3r::Print print;
    build_cubes(model, print, config, n, overlap);
    return print.validate(&warnings);
}

size_t count_opt_key(const std::vector<StringObjectException>& warnings, const std::string& key)
{
    return std::count_if(warnings.begin(), warnings.end(),
        [&](const StringObjectException& w) { return w.opt_key == key; });
}

// Make `default_acceleration` exceed the machine's extruding-acceleration limit.
void trigger_acceleration_warning(DynamicPrintConfig& c)
{
    c.set_key_value("machine_max_acceleration_extruding", new ConfigOptionFloats{ 100. });
    c.set_key_value("default_acceleration", new ConfigOptionFloatsNullable{ 100000. });
}

// Make `default_jerk` exceed the machine's jerk limit (junction deviation off so
// the jerk check is not skipped).
void trigger_jerk_warning(DynamicPrintConfig& c)
{
    c.set_key_value("machine_max_junction_deviation", new ConfigOptionFloats{ 0. });
    c.set_key_value("machine_max_jerk_x", new ConfigOptionFloats{ 1. });
    c.set_key_value("machine_max_jerk_y", new ConfigOptionFloats{ 1. });
    c.set_key_value("default_jerk", new ConfigOptionFloatsNullable{ 9999. });
}

// Precise outer wall is ignored unless the wall sequence is inner-outer.
void trigger_precise_wall_warning(DynamicPrintConfig& c)
{
    c.set_key_value("precise_outer_wall", new ConfigOptionBool(true));
    c.set_key_value("wall_sequence", new ConfigOptionEnum<WallSequence>(WallSequence::OuterInner));
}

} // namespace

// ---------------------------------------------------------------------------
// {first_object_name} filename placeholder
// ---------------------------------------------------------------------------
namespace {

// Add a printable 20mm cube named `name` to `model`; returns it so the caller can tweak it.
ModelObject* add_named_cube(Model& model, const std::string& name)
{
    ModelObject* obj = model.add_object();
    obj->name = name;
    obj->add_volume(make_cube(20.0, 20.0, 20.0));
    obj->add_instance();
    obj->ensure_on_bed();
    return obj;
}

// Resolve `format` to an output file name for a print of `model`. `filename_base`, when set,
// is the saved-project name passed to output_filename().
std::string resolved_output_name(Model& model, const std::string& format, const std::string& filename_base = {})
{
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_key_value("filename_format", new ConfigOptionString(format));

    Print print;
    for (ModelObject* obj : model.objects)
        print.auto_assign_extruders(obj);
    print.apply(model, config);
    return print.output_filename(filename_base);
}

} // namespace

TEST_CASE("Print: {first_object_name} names the first printable object on the plate", "[Print]")
{
    Model model;

    SECTION("uses the object's name") {
        add_named_cube(model, "WidgetPart");
        CHECK(resolved_output_name(model, "{first_object_name}") == "WidgetPart.gcode");
    }

    SECTION("picks the first when several objects are printable") {
        add_named_cube(model, "FirstPart");
        add_named_cube(model, "SecondPart");
        CHECK(resolved_output_name(model, "{first_object_name}") == "FirstPart.gcode");
    }

    SECTION("skips objects outside the print volume (e.g. on another plate)") {
        // First in model order, but not on the current plate, so is_printable() is false.
        add_named_cube(model, "OtherPlatePart")->instances.front()->print_volume_state = ModelInstancePVS_Fully_Outside;
        add_named_cube(model, "OnPlatePart");
        CHECK(resolved_output_name(model, "{first_object_name}") == "OnPlatePart.gcode");
    }

    SECTION("is empty when the object has no name") {
        add_named_cube(model, "");
        CHECK(resolved_output_name(model, "part_{first_object_name}") == "part_.gcode");
    }
}

TEST_CASE("Print: {first_object_name} is not replaced by the saved-project file name", "[Print]")
{
    // Passing a saved-project file name as the filename_base must not change {first_object_name}.
    Model model;
    add_named_cube(model, "WidgetPart");
    CHECK(resolved_output_name(model, "{first_object_name}", "SavedProject") == "WidgetPart.gcode");
}

TEST_CASE("Print::validate stacks independent warnings", "[Print][validate]")
{
    // Two unrelated checks (region precise-wall + machine acceleration) must each
    // contribute their own entry.
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    trigger_precise_wall_warning(config);
    trigger_acceleration_warning(config);

    std::vector<StringObjectException> warnings;
    StringObjectException err = validate_cubes(config, warnings);

    CHECK(err.string.empty());
    CHECK(warnings.size() >= 2);
    CHECK(count_opt_key(warnings, "precise_outer_wall") == 1);  // jump-to key is preserved
    for (const auto& w : warnings)
        CHECK(w.is_warning);                                   // every collected entry is a warning
}

TEST_CASE("Print::validate coalesces motion-ability warnings into one", "[Print][validate]")
{
    // The jerk/junction/acceleration checks are mutually exclusive (gated on a shared
    // key), so adding a second motion trigger must NOT add a second warning.
    DynamicPrintConfig accel_only = DynamicPrintConfig::full_print_config();
    trigger_acceleration_warning(accel_only);
    std::vector<StringObjectException> w_accel;
    CHECK(validate_cubes(accel_only, w_accel).string.empty());

    DynamicPrintConfig accel_and_jerk = DynamicPrintConfig::full_print_config();
    trigger_acceleration_warning(accel_and_jerk);
    trigger_jerk_warning(accel_and_jerk);
    std::vector<StringObjectException> w_both;
    CHECK(validate_cubes(accel_and_jerk, w_both).string.empty());

    CHECK(w_accel.size() >= 1);
    CHECK(w_both.size() == w_accel.size());  // the extra motion trigger collapses into the same warning
}

TEST_CASE("Print::validate reports the clumping-detection warning", "[Print][validate]")
{
    // A distinct single-shot path: clumping/wrapping detection without a prime tower warns
    // (and carries the enable_prime_tower jump-to key). enable_prime_tower must be off, as
    // the warning lives in the no-prime-tower branch.
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_key_value("enable_prime_tower", new ConfigOptionBool(false));
    config.set_key_value("enable_wrapping_detection", new ConfigOptionBool(true));

    std::vector<StringObjectException> warnings;
    StringObjectException err = validate_cubes(config, warnings);

    CHECK(err.string.empty());
    CHECK(count_opt_key(warnings, "enable_prime_tower") == 1);
}

TEST_CASE("Print::validate concatenates layered-clearance collisions into one warning", "[Print][validate]")
{
    // In by-layer mode, layered_print_cleareance_valid folds every too-close pair into a
    // single warning entry (newline-joined), unlike the per-check stacking above. Isolate
    // that entry by type so unrelated default-config warnings don't affect the assertion.
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();

    std::vector<StringObjectException> warnings;
    StringObjectException err = validate_cubes(config, warnings, /*n=*/3, /*overlap=*/true);

    CHECK(err.string.empty());
    auto is_layered = [](const StringObjectException& w) {
        return w.type == STRING_EXCEPT_OBJECT_COLLISION_IN_LAYER_PRINT; };
    REQUIRE(std::count_if(warnings.begin(), warnings.end(), is_layered) == 1);  // 3 objects, 2 collisions, 1 entry
    auto it = std::find_if(warnings.begin(), warnings.end(), is_layered);
    CHECK(it->string.find('\n') != std::string::npos);  // the collisions were concatenated
}

TEST_CASE("Print::validate tolerates a null warnings pointer", "[Print][validate]")
{
    // Callers may pass no warnings sink: a warning-producing config must not crash
    // and must still return without a blocking error.
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    trigger_precise_wall_warning(config);
    trigger_acceleration_warning(config);

    Slic3r::Model model;
    Slic3r::Print print;
    build_cubes(model, print, config, /*n=*/1, /*overlap=*/false);

    StringObjectException err = print.validate();  // warnings == nullptr
    CHECK(err.string.empty());
}

TEST_CASE("A default slice emits perimeter, infill, and skirt", "[Print]")
{
    const std::string gcode = slice({ cube(20) }, {
        { "layer_height",               0.2 },
        { "initial_layer_print_height", 0.2 },
        { "z_hop",                      0 } // keep recorded Z at the printed height
    });
    CHECK(role_passes(gcode, "perimeter") > 0);
    CHECK(role_passes(gcode, "infill")    > 0);
    CHECK(role_passes(gcode, "skirt")     > 0);
    CHECK_THAT(max_z(gcode), Catch::Matchers::WithinAbs(20.0, 1e-4));
}

// The G-code carries a config-comment block describing the resolved settings. The
// per-region width lines are always present; the support and first-layer lines appear
// only when those features are configured.
TEST_CASE("G-code lists the resolved extrusion-width settings", "[Print]")
{
    const std::string gcode = slice({ cube(20) }, { { "initial_layer_line_width", 0 } });
    CHECK(gcode.find("; external perimeters extrusion width") != std::string::npos);
    CHECK(gcode.find("; perimeters extrusion width")          != std::string::npos);
    CHECK(gcode.find("; infill extrusion width")              != std::string::npos);
    CHECK(gcode.find("; solid infill extrusion width")        != std::string::npos);
    CHECK(gcode.find("; top infill extrusion width")          != std::string::npos);
    CHECK(gcode.find("; support material extrusion width")    == std::string::npos);
    CHECK(gcode.find("; first layer extrusion width")         == std::string::npos);
    CHECK(gcode.find("; layer_height")                        != std::string::npos);
    CHECK(gcode.find("; sparse_infill_density")               != std::string::npos);

    const std::string with_support = slice({ cube(20) }, {
        { "initial_layer_line_width", 0 }, { "enable_support", true }, { "raft_layers", 3 },
    });
    CHECK(with_support.find("; support material extrusion width") != std::string::npos);

    const std::string with_first_layer = slice({ cube(20) }, { { "initial_layer_line_width", "0.5" } });
    CHECK(with_first_layer.find("; first layer extrusion width") != std::string::npos);
}

// gcode_skip_config_block suppresses the resolved-settings block while leaving the
// header and executable blocks intact.
TEST_CASE("gcode_skip_config_block omits the resolved-settings comment block", "[Print]")
{
    const std::string gcode = slice({ cube(20) }, {
        { "gcode_skip_config_block", true },
        { "gcode_comments",         true },
    });
    CHECK(gcode.find("; CONFIG_BLOCK_START")     == std::string::npos);
    CHECK(gcode.find("; CONFIG_BLOCK_END")       == std::string::npos);
    CHECK(gcode.find("; layer_height =")         == std::string::npos);
    CHECK(gcode.find("; fill_density =")         == std::string::npos);
    CHECK(gcode.find("; HEADER_BLOCK_START")     != std::string::npos);
    CHECK(gcode.find("; EXECUTABLE_BLOCK_START") != std::string::npos);
}

// Custom G-code templates substitute placeholders during export.
TEST_CASE("Custom G-code placeholders are substituted", "[Print]")
{
    // [current_extruder] in the start G-code.
    CHECK(slice({ cube(20) }, { { "machine_start_gcode", "; Extruder [current_extruder]" } })
              .find("; Extruder 0") != std::string::npos);

    // [layer_num] / [layer_z] in the end G-code (a 20mm cube at 0.1mm is 200 layers).
    const std::string end_gcode = slice({ cube(20) }, {
        { "machine_end_gcode",          "; Layer_num [layer_num]\n; Layer_z [layer_z]" },
        { "layer_height",               0.1 },
        { "initial_layer_print_height", 0.1 },
    });
    CHECK(end_gcode.find("; Layer_num 199") != std::string::npos);
    CHECK(end_gcode.find("; Layer_z 20")    != std::string::npos);

    // printing_by_object_gcode is emitted between sequentially printed objects.
    CHECK(slice_two_cubes_arranged({
                    { "print_sequence",           "by object" },
                    { "printing_by_object_gcode", "; between-object-gcode" },
                })
              .find("; between-object-gcode") != std::string::npos);

    // [layer_num] keeps counting across sequentially printed objects (199 then 399).
    const std::string per_layer = slice_two_cubes_arranged({
        { "print_sequence",             "by object" },
        { "layer_change_gcode",         ";Layer:[layer_num] ([layer_z] mm)" },
        { "layer_height",               0.1 },
        { "initial_layer_print_height", 0.1 },
    });
    CHECK(per_layer.find(";Layer:199 ") != std::string::npos);
    CHECK(per_layer.find(";Layer:399 ") != std::string::npos);
}

TEST_CASE("export_gcode writes G-code without a result pointer", "[Print][export_gcode]")
{
    Print print;
    Model model;
    Slic3r::Test::init_print({cube(20)}, print, model);
    print.process();

    SECTION("non-BBL printer") {}
    SECTION("BBL printer") { print.is_BBL_printer() = true; }

    ScopedTemporaryFile temp(".gcode");
    REQUIRE_NOTHROW(print.export_gcode(temp.string(), nullptr, nullptr));

    std::ifstream in(temp.string());
    const std::string gcode((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

    REQUIRE_FALSE(gcode.empty());
}

TEST_CASE("Sequential printing follows model order", "[Print]")
{
    // Two objects of different heights, taller one added first. Orca prints
    // sequential objects in model order, so the taller one is printed first.
    const std::string gcode = Slic3r::Test::slice({ cube(20), Slic3r::make_cube(20, 20, 10) }, {
        { "print_sequence",             "by object" },
        { "layer_height",               0.2 },
        { "initial_layer_print_height", 0.2 },
        { "z_hop",                      0 }
    });

    // The first object's height is the peak Z reached before Z drops back to the
    // first layer (the object change). With by-object printing only an object
    // change returns Z to the bottom.
    double first_object_peak_z = 0.0;
    double running_peak        = 0.0;
    GCodeReader reader;
    reader.parse_buffer(gcode, [&] (GCodeReader& self, const GCodeReader::GCodeLine& line) {
        if (first_object_peak_z != 0.0 || !line.extruding(self)) return; // ignore travels (e.g. start-gcode Z lift)
        if (running_peak > 1.0 && self.z() < 1.0)
            first_object_peak_z = running_peak;
        else
            running_peak = std::max(running_peak, static_cast<double>(self.z()));
    });

    REQUIRE_THAT(first_object_peak_z, Catch::Matchers::WithinAbs(20.0, 0.3));
}

TEST_CASE("Sequential object spacing accounts for one moving gantry rather than two", "[Print][SequentialGantryGeometry]")
{
    const double gap = GENERATE(6., 9., 15.);
    Print print;
    Model model;
    place_two_cubes_apart(gap, {
        {"print_sequence", "by object"}, {"extruder_clearance_radius", "75"},
        {"extruder_clearance_height_to_rod", "100"}, {"extruder_clearance_height_to_lid", "100"},
        {"skirts", "0"}, {"brim_type", "no_brim"},
        {"sequential_print_gantry_geometry", R"({"slices":[{"height":0,"type":"convex","polygons":["-8,-8;8,-8;8,8;-8,8"]}]})"}
    }, print, model);
    const auto result = Print::sequential_print_clearance_valid(print);
    INFO(result.string);
    CHECK(result.string.empty() == (gap > 8.));
}

TEST_CASE("Modeled sequential clearance distinguishes side gaps from rear gaps", "[Print][SequentialGantryGeometry]")
{
    const bool along_y = GENERATE(false, true);
    const double gap = GENERATE(25., 35., 95.);
    auto config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({
        {"print_sequence", "by object"}, {"extruder_clearance_radius", "75"},
        {"extruder_clearance_height_to_rod", "100"}, {"extruder_clearance_height_to_lid", "100"},
        {"skirts", "0"}, {"brim_type", "no_brim"},
        {"sequential_print_gantry_geometry", R"({"slices":[{"height":0,"type":"convex","polygons":["-27,-90;28,-90;28,21;-27,21"]}]})"}
    });
    auto a = Test::cube(20.);
    auto b = Test::cube(20.);
    a.translate(50., 50., 0.);
    b.translate(along_y ? 50. : 70. + gap, along_y ? 70. + gap : 50., 0.);
    Print print;
    Model model;
    Test::init_print(std::vector<TriangleMesh>{std::move(a), std::move(b)}, print, model, config, nullptr, false);
    const auto result = Print::sequential_print_clearance_valid(print);
    INFO(result.string);
    CHECK(result.string.empty() == (gap > (along_y ? 90. : 28.)));
}

// A sequential (by-object) print must publish the print-level nozzle group result just
// like a by-layer print, so custom g-code can index the per-nozzle placeholder tables
// (e.g. nozzle_diameter_at_nozzle_id[]) instead of failing on an empty vector.
TEST_CASE("Sequential printing publishes the nozzle group result", "[Print][MultiNozzle]")
{
    SECTION("process() publishes the result") {
        Print print;
        Model model;
        place_two_cubes_apart(60.0, { { "print_sequence", "by object" } }, print, model);
        print.process();
        REQUIRE(print.get_layered_nozzle_group_result() != nullptr);
    }

    SECTION("start g-code can index the per-nozzle diameter table") {
        const std::string gcode = slice_two_cubes_arranged({
            { "print_sequence",      "by object" },
            { "machine_start_gcode", "{if nozzle_diameter_at_nozzle_id[0] > 0}; SEQ-ND-OK\n{endif}" },
        });
        CHECK(gcode.find("; SEQ-ND-OK") != std::string::npos);
    }
}
