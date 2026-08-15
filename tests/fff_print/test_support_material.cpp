#include <catch2/catch_all.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <mutex>
#include <string>

#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/Layer.hpp"

#include "test_helpers.hpp" // get access to init_print, etc

using namespace Slic3r::Test;
using namespace Slic3r;

namespace {

size_t count_nonplanar_support_contact_paths(const ExtrusionEntity &entity,
                                             ExtrusionRole role)
{
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity))
        return path->nonplanar_surface && path->role() == role;
    size_t count = 0;
    if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
        for (const ExtrusionPath &path : multipath->paths)
            count += path.nonplanar_surface && path.role() == role;
    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
        for (const ExtrusionPath &path : loop->paths)
            count += path.nonplanar_surface && path.role() == role;
    } else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        for (const ExtrusionEntity *child : collection->entities)
            count += count_nonplanar_support_contact_paths(*child, role);
    }
    return count;
}

size_t count_nonplanar_support_contact_paths(const PrintObject &object,
                                             ExtrusionRole role)
{
    size_t count = 0;
    for (const SupportLayer *layer : object.support_layers())
        for (const ExtrusionEntity *entity : layer->support_fills.entities)
            count += count_nonplanar_support_contact_paths(*entity, role);
    return count;
}

size_t count_support_paths(const ExtrusionEntity &entity, ExtrusionRole role)
{
    if (const auto *path = dynamic_cast<const ExtrusionPath *>(&entity))
        return path->role() == role;
    size_t count = 0;
    if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath *>(&entity)) {
        for (const ExtrusionPath &path : multipath->paths)
            count += path.role() == role;
    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop *>(&entity)) {
        for (const ExtrusionPath &path : loop->paths)
            count += path.role() == role;
    } else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection *>(&entity)) {
        for (const ExtrusionEntity *child : collection->entities)
            count += count_support_paths(*child, role);
    }
    return count;
}

} // namespace

TEST_CASE("Floating extrusion islands warn unless generated support reaches the island", "[SupportMaterial][Warnings]")
{
    TriangleMesh model = make_cube(12., 12., 2.);
    TriangleMesh floating = make_cube(4., 4., 2.);
    floating.translate(Vec3f(4.f, 4.f, 5.f));
    model.merge(floating);

    auto has_midair_warning = [](const Print &print) {
        const auto state = print.objects().front()->step_state_with_warnings(posSupportMaterial);
        return std::any_of(state.warnings.begin(), state.warnings.end(), [](const PrintStateBase::Warning &warning) {
            return warning.current && warning.message.find("starting in midair") != std::string::npos;
        });
    };

    SECTION("unsupported island is reported") {
        Print print;
        init_and_process_print({model}, print, {{"enable_support", false}});
        REQUIRE(has_midair_warning(print));
    }

    SECTION("enabled supports suppress the warning") {
        Print print;
        init_and_process_print({model}, print, {{"enable_support", true}});
        REQUIRE_FALSE(has_midair_warning(print));
    }

    SECTION("enabled supports still warn when no support reaches the island") {
        Print print;
        init_and_process_print({model}, print, {
            {"enable_support", true},
            {"support_type",  std::string("normal(manual)")},
        });
        REQUIRE(has_midair_warning(print));
    }
}

TEST_CASE("Three raft layers are created", "[SupportMaterial]")
{
	Slic3r::Print print;
	Slic3r::Test::init_and_process_print({ cube(20) }, print, {
        { "enable_support", 1 },
        { "raft_layers",    3 }
		});
    REQUIRE(print.objects().front()->support_layers().size() == 3);
}

TEST_CASE("Surface-following support ironing projects the supported interface",
          "[SupportMaterial][Nonplanar]")
{
    TriangleMesh sloped_plate = make_cube(20., 12., 1.2);
    sloped_plate.rotate_y(float(10. * M_PI / 180.));
    sloped_plate.translate(Vec3f(0.f, 0.f, 8.f));

    Print print;
    init_and_process_print({sloped_plate}, print, {
        {"enable_support", true},
        {"support_type", std::string("normal(auto)")},
        {"support_ironing", true},
        {"support_ironing_nonplanar", true},
        {"support_ironing_nonplanar_max_angle", 45.},
        {"support_ironing_nonplanar_resolution", 0.2},
        {"support_interface_top_layers", 3},
        {"support_top_z_distance", 0.2},
        {"dont_support_bridges", false},
        {"layer_height", 0.2}
    });

    const PrintObject &object = *print.objects().front();
    size_t contact_layers = 0;
    size_t interface_paths = 0;
    size_t ironing_paths = 0;
    for (const SupportLayer *layer : object.support_layers()) {
        contact_layers += layer->has_top_interface_contact;
        for (const ExtrusionEntity *entity : layer->support_fills.entities) {
            interface_paths += count_support_paths(*entity, erSupportMaterialInterface);
            ironing_paths += count_support_paths(*entity, erIroning);
        }
    }
    INFO("support layers=" << object.support_layers().size()
         << ", contact layers=" << contact_layers
         << ", interface paths=" << interface_paths
         << ", ironing paths=" << ironing_paths);
    REQUIRE(count_nonplanar_support_contact_paths(
        object, erSupportMaterialInterface) > 0);
    REQUIRE(count_nonplanar_support_contact_paths(object, erIroning) > 0);
}

TEST_CASE("Enforced support layers are generated", "[SupportMaterial]")
{
    // enforce_support_layers forces support on the first N layers even with support off.
    Slic3r::Print baseline;
    Slic3r::Test::init_and_process_print({ TestMesh::overhang }, baseline, {
        { "enable_support",         0 },
        { "enforce_support_layers", 0 }
    });
    REQUIRE(baseline.objects().front()->support_layers().empty());

    Slic3r::Print enforced;
    Slic3r::Test::init_and_process_print({ TestMesh::overhang }, enforced, {
        { "enable_support",         0 },
        { "enforce_support_layers", 100 }
    });
    REQUIRE(enforced.objects().front()->support_layers().size() > 0);
}

TEST_CASE("Support toolpath generation reports its real internal stages",
          "[SupportMaterial][Progress]")
{
    Print print;
    Model model;
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({{"enable_support", "1"}});
    init_print({TestMesh::overhang}, print, model, config);

    std::mutex messages_mutex;
    std::vector<std::string> messages;
    print.set_status_callback([&](const PrintBase::SlicingStatus &status) {
        std::scoped_lock lock(messages_mutex);
        messages.emplace_back(status.text);
    });
    REQUIRE_NOTHROW(print.process());

    const auto contains = [&](const char *needle) {
        return std::any_of(messages.begin(), messages.end(),
            [needle](const std::string &message) {
                return message.find(needle) != std::string::npos;
            });
    };
    CHECK(contains("Generating support region paths: layer"));
    CHECK(contains("Assembling support extrusion layers: layer"));
}

SCENARIO("Support layer Z honors contact distance", "[SupportMaterial]")
{
    // Box h = 20mm, hole bottom at 5mm, hole height 10mm (top edge at 15mm).
    TriangleMesh mesh = Slic3r::Test::mesh(Slic3r::Test::TestMesh::cube_with_hole);
    mesh.rotate_x(float(M_PI / 2));

	auto check = [](Slic3r::Print &print, bool &first_support_layer_height_ok, bool &layer_height_minimum_ok, bool &layer_height_maximum_ok)
	{
        ConstSupportLayerPtrsAdaptor support_layers = print.objects().front()->support_layers();

		first_support_layer_height_ok = support_layers.front()->print_z == print.config().initial_layer_print_height.value;

		layer_height_minimum_ok = true;
		layer_height_maximum_ok = true;
		double min_layer_height = print.config().min_layer_height.values.front();
		double max_layer_height = print.config().nozzle_diameter.values.front();
		if (print.config().max_layer_height.values.front() > EPSILON)
			max_layer_height = std::min(max_layer_height, print.config().max_layer_height.values.front());
		for (size_t i = 1; i < support_layers.size(); ++ i) {
			if (support_layers[i]->print_z - support_layers[i - 1]->print_z < min_layer_height - EPSILON)
				layer_height_minimum_ok = false;
			if (support_layers[i]->print_z - support_layers[i - 1]->print_z > max_layer_height + EPSILON)
				layer_height_maximum_ok = false;
		}
	};

    GIVEN("A print object having one modelObject") {
        WHEN("Layer height = 0.2 and first layer height = 0.4") {
			Slic3r::Print print;
			Slic3r::Test::init_and_process_print({ mesh }, print, {
                { "enable_support",             1 },
                { "layer_height",               0.2 },
                { "initial_layer_print_height", 0.4 },
                { "dont_support_bridges",       false },
			});
			bool first_layer_ok, layer_min_ok, layer_max_ok;
            check(print, first_layer_ok, layer_min_ok, layer_max_ok);
            THEN("First layer height is honored")			{ REQUIRE(first_layer_ok == true); }
            THEN("No null or negative support layers")		{ REQUIRE(layer_min_ok == true); }
            THEN("No layers thicker than nozzle diameter")	{ REQUIRE(layer_max_ok == true); }
        }
        WHEN("Layer height = 0.2 and first layer height = 0.3") {
			Slic3r::Print print;
			Slic3r::Test::init_and_process_print({ mesh }, print, {
                { "enable_support",             1 },
                { "layer_height",               0.2 },
                { "initial_layer_print_height", 0.3 },
                { "dont_support_bridges",       false },
            });
            bool first_layer_ok, layer_min_ok, layer_max_ok;
            check(print, first_layer_ok, layer_min_ok, layer_max_ok);
            THEN("First layer height is honored")			{ REQUIRE(first_layer_ok == true); }
            THEN("No null or negative support layers")		{ REQUIRE(layer_min_ok == true); }
            THEN("No layers thicker than nozzle diameter")	{ REQUIRE(layer_max_ok == true); }
        }
    }
}

// extrude_support once held a `static` lambda capturing `this`, so a second export in the
// same process dereferenced a returned stack frame (ASan: stack-use-after-return).
TEST_CASE("Support G-code emission survives a second slice in the same process", "[SupportMaterial][Regression]")
{
    const std::string first = slice({ TestMesh::overhang }, { { "enable_support", 1 } });
    REQUIRE(! layers_with_role(first, "support").empty());

    const std::string second = slice({ TestMesh::overhang }, { { "enable_support", 1 } });
    REQUIRE(! layers_with_role(second, "support").empty());
}

TEST_CASE("Top support interface uses and restores its temperature override", "[SupportMaterial][Regression]")
{
    constexpr int interface_temperature = 150;
    constexpr int normal_temperature = 215;
    Print print;
    init_and_process_print({TestMesh::overhang}, print, {
        {"enable_support", 1},
        {"support_interface_top_layers", 2},
        {"support_interface_top_temperature", interface_temperature},
        {"nozzle_temperature", normal_temperature}
    });

    const auto support_layers = print.objects().front()->support_layers();
    REQUIRE(std::any_of(support_layers.begin(), support_layers.end(), [](const SupportLayer *layer) {
        return layer->has_top_interface_contact;
    }));

    const std::string output = gcode(print);
    CHECK(output.find("M109 S" + std::to_string(interface_temperature)) != std::string::npos);
    CHECK(output.find("M109 S" + std::to_string(normal_temperature)) != std::string::npos);
}

namespace {

std::map<double, double> first_model_feedrate_by_layer(const std::string &gcode)
{
    std::map<double, double> feedrates;
    GCodeReader reader;
    reader.parse_buffer(gcode, [&feedrates](GCodeReader &self, const GCodeReader::GCodeLine &line) {
        if (!line.extruding(self) || line.dist_XY(self) <= 0. ||
            line.comment().find("support") != std::string_view::npos)
            return;
        feedrates.try_emplace(self.z(), line.new_F(self) / 60.);
    });
    return feedrates;
}

} // namespace

TEST_CASE("Model speed recovers after a dissimilar top support interface",
          "[SupportMaterial][GCode]")
{
    constexpr double slowed_speed = 12.;
    DynamicPrintConfig config = multifilament_config(2, {
        {"enable_support", true},
        {"dont_support_bridges", false},
        {"support_filament", 1},
        {"support_interface_filament", 2},
        {"support_interface_top_layers", 2},
        {"filament_type", "PLA;PETG"},
        {"slow_down_layer_above_dissimilar_support_interface", true},
        {"dissimilar_support_interface_speed", slowed_speed},
        {"dissimilar_support_interface_speed_layers", 2},
        {"outer_wall_speed", 60},
        {"inner_wall_speed", 60},
        {"sparse_infill_speed", 60},
        {"internal_solid_infill_speed", 60},
        {"top_surface_speed", 60},
        {"bridge_speed", 60},
        {"internal_bridge_speed", 60},
        {"filament_max_volumetric_speed", "100,100"},
    });

    Print print;
    init_and_process_print({TestMesh::overhang}, print, config);
    const PrintObject &object = *print.objects().front();
    const Layer *contact_layer = nullptr;
    for (const Layer *layer : object.layers()) {
        const SupportLayer *support = support_layer_below(*layer);
        if (support != nullptr && support->has_top_interface_contact)
            contact_layer = layer;
    }
    REQUIRE(contact_layer != nullptr);

    const std::map<double, double> feedrates = first_model_feedrate_by_layer(gcode(print));
    const double contact_z = contact_layer->print_z;
    const auto speed_at = [&feedrates](double z) {
        const auto it = std::find_if(feedrates.begin(), feedrates.end(), [z](const auto &entry) {
            return std::abs(entry.first - z) < 0.001;
        });
        REQUIRE(it != feedrates.end());
        return it->second;
    };

    CHECK_THAT(speed_at(contact_z), Catch::Matchers::WithinAbs(slowed_speed, 0.01));
    REQUIRE(contact_layer->upper_layer != nullptr);
    CHECK_THAT(speed_at(contact_layer->upper_layer->print_z),
               Catch::Matchers::WithinAbs(0.5 * (slowed_speed + 60.), 0.01));
    REQUIRE(contact_layer->upper_layer->upper_layer != nullptr);
    CHECK_THAT(speed_at(contact_layer->upper_layer->upper_layer->print_z),
               Catch::Matchers::WithinAbs(60., 0.01));
}

TEST_CASE("Same-material support interface does not slow the model layer",
          "[SupportMaterial][GCode]")
{
    DynamicPrintConfig config = multifilament_config(2, {
        {"enable_support", true},
        {"dont_support_bridges", false},
        {"support_filament", 1},
        {"support_interface_filament", 2},
        {"support_interface_top_layers", 2},
        {"filament_type", "PLA;PLA"},
        {"slow_down_layer_above_dissimilar_support_interface", true},
        {"dissimilar_support_interface_speed", 12},
        {"outer_wall_speed", 60},
        {"inner_wall_speed", 60},
        {"bridge_speed", 60},
        {"filament_max_volumetric_speed", "100,100"},
    });

    Print print;
    init_and_process_print({TestMesh::overhang}, print, config);
    const PrintObject &object = *print.objects().front();
    const Layer *contact_layer = nullptr;
    for (const Layer *layer : object.layers()) {
        const SupportLayer *support = support_layer_below(*layer);
        if (support != nullptr && support->has_top_interface_contact)
            contact_layer = layer;
    }
    REQUIRE(contact_layer != nullptr);

    const std::map<double, double> feedrates = first_model_feedrate_by_layer(gcode(print));
    const auto speed = std::find_if(feedrates.begin(), feedrates.end(), [contact_layer](const auto &entry) {
        return std::abs(entry.first - contact_layer->print_z) < 0.001;
    });
    REQUIRE(speed != feedrates.end());
    CHECK_THAT(speed->second, Catch::Matchers::WithinAbs(60., 0.01));
}
