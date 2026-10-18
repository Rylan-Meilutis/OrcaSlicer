#include <catch2/catch_all.hpp>

#include <algorithm>
#include <sstream>
#include <string>

#include "libslic3r/calib.hpp"
#include "libslic3r/SeamCalibration.hpp"
#include "libslic3r/TemperatureTowerPresets.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/PrintConfig.hpp"

using namespace Slic3r;

TEST_CASE("Seam calibration creates closed repeatable coupons and isolated setting ranges", "[Calib][SeamCalibration]")
{
    const auto kind = GENERATE(SeamCalibrationKind::Length, SeamCalibrationKind::StartHeight,
        SeamCalibrationKind::Flow, SeamCalibrationKind::Speed, SeamCalibrationKind::InnerEntry);
    SeamCalibrationParams params;
    params.kind = kind;
    params.count = 2;
    if (kind == SeamCalibrationKind::Flow) { params.start = 0.8; params.step = 0.1; }
    const auto preset = DynamicPrintConfig::full_print_config();
    Model model;
    make_seam_calibration(model, preset, params);
    const bool inner = kind == SeamCalibrationKind::InnerEntry;
    REQUIRE(model.objects.size() == (inner ? 8 : 12));
    for (size_t i = 0; i < model.objects.size(); ++i) {
        const auto &object = *model.objects[i];
        const auto &mesh = object.volumes.front()->mesh();
        CHECK(its_volume(mesh.its) > 0.);
        CHECK_THAT(mesh.bounding_box().size().z(), Catch::Matchers::WithinAbs(inner ? 8. : 24., 1e-5));
        // Every mesh edge must be shared by exactly two triangles, with opposite winding.
        std::map<std::pair<int,int>, std::pair<int,int>> edges;
        for (const auto &triangle : mesh.its.indices)
            for (int e = 0; e < 3; ++e) {
                const int a = triangle[e], b = triangle[(e + 1) % 3];
                auto &edge = edges[std::minmax(a,b)];
                ++edge.first; edge.second += a < b ? 1 : -1;
            }
        for (const auto &[key, edge] : edges) {
            CHECK(edge.first == 2);
            CHECK(edge.second == 0);
        }
        CHECK_FALSE(object.config.get().opt_bool("arc_overhang_enabled"));
        CHECK(object.config.get().opt_bool("seam_start_on_inner_wall") == (inner && i >= 4));
        CHECK(object.layer_config_ranges.size() == (inner ? 0 : 3));
        if (!inner) {
            CHECK(object.layer_config_ranges.begin()->second.opt_serialize("seam_slope_type") == "none");
            for (auto range = std::next(object.layer_config_ranges.begin()); range != object.layer_config_ranges.end(); ++range)
                CHECK(range->second.size() == 2); // swept value plus required layer height
        }
    }
    CHECK(preset == DynamicPrintConfig::full_print_config());
}

TEST_CASE("Seam calibration rejects unsafe or unbounded sweep values", "[Calib][SeamCalibration]")
{
    SeamCalibrationParams params;
    params.count = 100;
    CHECK_THROWS_AS(seam_calibration_values(params), std::invalid_argument);
    params.count = 4;
    params.start = 39.;
    CHECK_THROWS_AS(seam_calibration_values(params), std::invalid_argument);
    params.kind = SeamCalibrationKind::StartHeight;
    params.start = 100.;
    CHECK_THROWS_AS(seam_calibration_values(params), std::invalid_argument);
}

TEST_CASE("Temperature tower presets have unique materials and valid descending steps", "[Calib]")
{
    for (size_t i = 0; i < temperature_tower_presets.size(); ++i) {
        const auto &preset = temperature_tower_presets[i];
        INFO(preset.material);
        CHECK(preset.start <= 500);
        CHECK(preset.end >= 155);
        CHECK(preset.start > preset.end);
        CHECK(preset.start % 5 == 0);
        CHECK(preset.end % 5 == 0);
        for (size_t j = 0; j < i; ++j)
            CHECK(std::string(preset.material) != temperature_tower_presets[j].material);
    }
    CHECK(std::string(temperature_tower_presets[5].material) == "PA-CF");
    CHECK(temperature_tower_presets[5].start == 320);
    CHECK(std::string(temperature_tower_presets[6].material) == "PET-CF");
    CHECK(temperature_tower_presets[6].start == 320);
}

TEST_CASE("Calibration disables inherited and local arc overhangs without changing the preset", "[Calib][Regression]")
{
    DynamicPrintConfig preset;
    preset.set_key_value("arc_overhang_enabled", new ConfigOptionBool(true));
    Model model;
    for (int i = 0; i < 2; ++i) {
        ModelObject *object = model.add_object("cube", "", make_cube(20, 20, 20));
        object->config.set_key_value("wall_loops", new ConfigOptionInt(3));
        object->volumes.front()->config.set_key_value("arc_overhang_enabled", new ConfigOptionBool(true));
        object->layer_config_ranges[{1., 2.}].set_key_value("arc_overhang_enabled", new ConfigOptionBool(true));
    }
    disable_arc_overhangs_for_calibration(model);
    for (const ModelObject *object : model.objects) {
        DynamicPrintConfig resolved = preset;
        resolved.apply(object->config.get());
        CHECK_FALSE(resolved.opt_bool("arc_overhang_enabled"));
        CHECK_FALSE(object->volumes.front()->config.get().opt_bool("arc_overhang_enabled"));
        CHECK_FALSE(object->layer_config_ranges.begin()->second.get().opt_bool("arc_overhang_enabled"));
        CHECK(object->config.opt_int("wall_loops") == 3);
    }
    CHECK(preset.opt_bool("arc_overhang_enabled"));
}

namespace {

// The width-resolution getters are protected; expose them so the resolution can be asserted directly.
struct PaPatternProbe : public CalibPressureAdvancePattern
{
    using CalibPressureAdvancePattern::CalibPressureAdvancePattern;
    using CalibPressureAdvancePattern::line_width;
    using CalibPressureAdvancePattern::line_width_first_layer;
};

} // namespace

TEST_CASE("Zero calibration line width resolves to a positive default", "[Calib][Regression]")
{
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({
        {"line_width", "0"},
        {"initial_layer_line_width", "0"},
    });

    Model model;
    model.add_object("cube", "", make_cube(20, 20, 20))->add_instance();

    Calib_Params params;
    params.mode = CalibMode::Calib_PA_Pattern;

    PaPatternProbe pattern(params, config, /* is_bbl_machine */ true, *model.objects.front(), Vec3d(0, 0, 0));

    REQUIRE(pattern.line_width() > 0.);
    REQUIRE(pattern.line_width_first_layer() > 0.);
}

namespace {

struct EndState { double final_e; double max_e; };

EndState simulate_absolute_e(const std::string &gcode)
{
    double final_e = 0.;
    double max_e   = 0.;

    std::istringstream lines(gcode);
    std::string        line;
    while (std::getline(lines, line)) {
        std::istringstream words(line);
        std::string        op;
        if (!(words >> op))
            continue;
        if (op != "G1" && op != "G0" && op != "G92")
            continue;

        std::string word;
        while (words >> word) {
            if (word.size() >= 2 && word[0] == 'E') {
                final_e = std::stod(word.substr(1));
                max_e   = std::max(max_e, final_e);
                break;
            }
        }
    }

    return {final_e, max_e};
}

} // namespace

TEST_CASE("PA pattern resets the extruder after the final layer in absolute E mode", "[Calib][Regression]")
{
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({
        {"use_relative_e_distances", "0"},
        {"line_width", "0.45"},
        {"initial_layer_line_width", "0.45"},
    });

    Model model;
    model.add_object("cube", "", make_cube(20, 20, 20))->add_instance();

    Calib_Params params;
    params.mode  = CalibMode::Calib_PA_Pattern;
    params.start = 0.;
    params.end   = 0.08;
    params.step  = 0.002;

    CalibPressureAdvancePattern pattern(params, config, /* is_bbl_machine */ false, *model.objects.front(), Vec3d(0, 0, 0));
    const CustomGCode::Info     info = pattern.generate_custom_gcodes(config, /* is_bbl_machine */ false, *model.objects.front(),
                                                                      Vec3d(0, 0, 0));

    std::string gcode;
    for (const CustomGCode::Item &item : info.gcodes)
        gcode += item.extra;

    const EndState state = simulate_absolute_e(gcode);

    REQUIRE(state.max_e > 1.);
    REQUIRE_THAT(state.final_e, Catch::Matchers::WithinAbs(0., 1e-9));
}
