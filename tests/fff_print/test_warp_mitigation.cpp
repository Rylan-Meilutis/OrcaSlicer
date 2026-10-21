#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "test_helpers.hpp"
#include "libslic3r/Config.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/PrintBase.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/WarpMitigation.hpp"
#include <algorithm>
#include <vector>

using namespace Slic3r;
using Catch::Matchers::WithinAbs;

TEST_CASE("Automatic warp mitigation changes effective adhesion without changing the project", "[WarpMitigation]")
{
    auto config = Test::multifilament_config(3, {
        {"filament_type", "ASA;ASA;TPU"}, {"brim_type", "painted"},
        {"brim_width", "0"}, {"brim_object_gap", "0.2"},
        {"automatic_warp_mitigation", "1"}});
    Model model;
    Print print;
    Test::init_print({make_cube(80., 70., 1.)}, print, model, config);
    REQUIRE(assess_warp_risk(config, model).elevated());
    REQUIRE(print.objects().size() == 1);
    CHECK(print.objects().front()->config().brim_type.value == btOuterOnly);
    CHECK_THAT(print.objects().front()->config().brim_width.value,
        WithinAbs(assess_warp_risk(config, model).brim_mm(), 1e-6));
    CHECK(config.option<ConfigOptionEnum<BrimType>>("brim_type")->value == btPainted);
    CHECK_THAT(config.opt_float("brim_width"), WithinAbs(0., 1e-6));
    print.process();
    CHECK(print.has_brim());
    REQUIRE_FALSE(print.skirt_brim_groups().empty());
    CHECK_FALSE(print.skirt_brim_groups().front().brims.empty());
    std::vector<StringObjectException> mitigated_warnings;
    print.validate(&mitigated_warnings);
    CHECK(std::any_of(mitigated_warnings.begin(), mitigated_warnings.end(), [](const auto &warning) {
        return warning.opt_key == "automatic_warp_mitigation" && !warning.is_warning;
    }));

    config.set_key_value("automatic_warp_mitigation", new ConfigOptionBool(false));
    print.apply(model, config);
    CHECK(print.objects().front()->config().brim_type.value == btPainted);
    CHECK_THAT(print.objects().front()->config().brim_width.value, WithinAbs(0., 1e-6));
    std::vector<StringObjectException> warnings;
    print.validate(&warnings);
    CHECK(std::any_of(warnings.begin(), warnings.end(), [](const auto &warning) {
        return warning.opt_key == "automatic_warp_mitigation" && warning.is_warning;
    }));
}

TEST_CASE("Unused susceptible filament does not change PLA adhesion", "[WarpMitigation]")
{
    auto config = Test::multifilament_config(2, {
        {"filament_type", "PLA;ASA"}, {"brim_type", "no_brim"},
        {"curr_bed_type", "Textured PEI Plate"}, {"textured_plate_temp", "60;100"},
        {"chamber_temperature", "20;20"}, {"bottom_shell_thickness", "1"},
        {"brim_width", "0"}, {"automatic_warp_mitigation", "1"}});
    Model model;
    Print print;
    Test::init_print({Test::cube(10.)}, print, model, config);
    CHECK_FALSE(assess_warp_risk(config, model).elevated());
    REQUIRE(print.objects().size() == 1);
    CHECK(print.objects().front()->config().brim_type.value == btNoBrim);
    std::vector<StringObjectException> warnings;
    print.validate(&warnings);
    CHECK(std::none_of(warnings.begin(), warnings.end(), [](const auto &warning) {
        return warning.opt_key == "automatic_warp_mitigation";
    }));
}

TEST_CASE("Warp risk remains with a brim or raft and clears when the gradient disappears", "[WarpMitigation]")
{
    const int raft_layers = GENERATE(0, 3);
    auto config = Test::multifilament_config(1, {
        {"filament_type", "ASA"}, {"brim_type", "outer_only"},
        {"brim_width", "12"}, {"brim_object_gap", "0"},
        {"raft_layers", raft_layers}, {"automatic_warp_mitigation", "0"}});
    Model model;
    Print print;
    Test::init_print({make_cube(80., 70., 1.)}, print, model, config);
    std::vector<StringObjectException> warnings;
    print.validate(&warnings);
    CHECK(std::count_if(warnings.begin(), warnings.end(), [](const auto &warning) {
        return warning.opt_key == "automatic_warp_mitigation";
    }) == 1);

    config.set_deserialize_strict("chamber_temperature", "150");
    print.apply(model, config);
    warnings.clear();
    print.validate(&warnings);
    CHECK(std::none_of(warnings.begin(), warnings.end(), [](const auto &warning) {
        return warning.opt_key == "automatic_warp_mitigation";
    }));
}

TEST_CASE("New material families apply geometry based mitigation without editing saved settings", "[WarpMitigation]")
{
    const auto type = GENERATE("PLA", "PP", "PA", "PA6", "PA11", "PA12", "POM", "PPS", "PEEK", "PPS-CF", "PA-GF");
    auto config = Test::multifilament_config(1, {
        {"filament_type", type}, {"brim_type", "no_brim"}, {"brim_width", "0"},
        {"curr_bed_type", "Textured PEI Plate"}, {"textured_plate_temp", "100"},
        {"chamber_temperature", "20"}, {"nozzle_temperature", "400"},
        {"bottom_shell_thickness", "1"}, {"layer_height", "0.2"},
        {"automatic_warp_mitigation", "1"}});
    Model model;
    Print print;
    Test::init_print({make_cube(80., 70., 1.)}, print, model, config);
    const auto risk = assess_warp_risk(config, model);
    CHECK(risk.unknown_materials.empty());
    CHECK(risk.elevated());
    CHECK(risk.crystallization_contraction_mm > 0.);
    REQUIRE(print.objects().size() == 1);
    CHECK(print.objects().front()->config().brim_type.value == btOuterOnly);
    CHECK(config.option<ConfigOptionEnum<BrimType>>("brim_type")->value == btNoBrim);
    config.set_deserialize_strict("automatic_warp_mitigation", "0");
    print.apply(model, config);
    CHECK(print.objects().front()->config().brim_type.value == btNoBrim);
}
