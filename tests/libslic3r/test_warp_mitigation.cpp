#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/Config.hpp"
#include "libslic3r/libslic3r.h"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Semver.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/WarpMitigation.hpp"
#include "libslic3r/Format/bbs_3mf.hpp"
#include "libslic3r/Preset.hpp"
#include "../test_utils.hpp"
#include <cstdlib>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

using namespace Slic3r;
using Catch::Matchers::WithinAbs;

TEST_CASE("Warp material coverage distinguishes resin screening from composite properties", "[WarpMitigation]")
{
    const auto type = GENERATE("ABS", "ASA", "PC", "HIPS", "PETG", "PMMA", "PLA", "PP", "PA", "PA6", "PA11", "PA12", "POM", "PPS", "PEEK");
    const auto material = warp_material(type);
    REQUIRE(material.has_value());
    CHECK_FALSE(material->resin_screening);
    CHECK(material->cte_min > 0.);
    CHECK(material->cte_max >= material->cte_min);
    CHECK(material->modulus_max > 0.);
    for (const std::string suffix : {"-CF", "-GF"}) {
        const auto filled = warp_material(std::string(type) + suffix);
        REQUIRE(filled.has_value());
        CHECK(filled->resin_screening);
        WarpRisk risk;
        assess_thermal_strip(risk, *filled, 100., 1., 100., 20., 400.);
        CHECK(risk.elevated());
        // Matrix stiffness must never be reported as composite stored energy.
        CHECK_THAT(risk.restrained_energy_j_m2, WithinAbs(0., 1e-12));
    }
    CHECK_FALSE(warp_material("unknown").has_value());
    CHECK_FALSE(warp_material("unknown-CF").has_value());
    CHECK_FALSE(warp_material("PLA-AERO").has_value());
}

TEST_CASE("Additional amorphous materials retain geometry and temperature dependent screening", "[WarpMitigation]")
{
    const auto type = GENERATE("HIPS", "PETG", "PMMA");
    const auto material = *warp_material(type);
    WarpRisk small, large, uniform;
    assess_thermal_strip(small, material, 5., 1., 100., 20., 250.);
    assess_thermal_strip(large, material, 100., 1., 100., 20., 250.);
    assess_thermal_strip(uniform, material, 100., 1., 20., 20., 250.);
    CHECK_FALSE(small.elevated());
    CHECK(large.elevated());
    CHECK_FALSE(uniform.elevated());
    CHECK(uniform.free_contraction_mm > 0.);
    CHECK_THAT(large.free_bow_mm, WithinAbs(400. * small.free_bow_mm, 1e-9));
    CHECK(warp_material(std::string(type) + "-CF")->resin_screening);
    CHECK(warp_material(std::string(type) + "-GF")->resin_screening);
}

TEST_CASE("Semicrystalline shrinkage has separate thermal and phase contributions", "[WarpMitigation]")
{
    const auto type = GENERATE("PLA", "PP", "PA6", "PA11", "PA12", "POM", "PPS", "PEEK");
    const auto material = *warp_material(type);
    REQUIRE(material.melting_max > material.tg_max);
    REQUIRE(material.phase_contraction_bound() > 0.);
    CHECK(material.phase_contraction_bound() < 0.1);
    WarpRisk reference, uniform, wider, hotter;
    assess_thermal_strip(reference, material, 100., 1., 100., 20., 400.);
    assess_thermal_strip(uniform, material, 100., 1., 20., 20., 400.);
    assess_thermal_strip(wider, material, 200., 1., 100., 20., 400.);
    assess_thermal_strip(hotter, material, 100., 1., 100., 20., 500.);
    CHECK(reference.free_bow_mm > 0.);
    CHECK(reference.free_contraction_mm > 0.);
    CHECK_THAT(reference.crystallization_contraction_mm,
        WithinAbs(100. * (1. - std::cbrt(material.amorphous_density / material.crystalline_density)), 1e-9));
    CHECK_THAT(wider.crystallization_contraction_mm, WithinAbs(2. * reference.crystallization_contraction_mm, 1e-9));
    CHECK_THAT(wider.free_bow_mm, WithinAbs(4. * reference.free_bow_mm, 1e-9));
    CHECK_THAT(hotter.free_contraction_mm, WithinAbs(reference.free_contraction_mm, 1e-9));
    CHECK(uniform.crystallization_contraction_mm > 0.);
    CHECK_FALSE(uniform.elevated()); // Uniform phase contraction is not a curvature measurement.
}

TEST_CASE("Material spelling aliases retain explicit grade uncertainty", "[WarpMitigation]")
{
    const auto type = GENERATE("pps-cf", "PPSCF", " PPS_CF ", "pps cf");
    const auto material = warp_material(type);
    REQUIRE(material.has_value());
    CHECK(material->resin_screening);
    CHECK_THAT(material->melting_max, WithinAbs(290., 1e-12));
    CHECK(warp_material("nylon")->family_screening);
    CHECK_FALSE(warp_material("nylon-12")->family_screening);
    CHECK(warp_material("nylon-cf")->family_screening);
    CHECK_FALSE(warp_material("PA-unknown").has_value());
}

TEST_CASE("Invalid thermal strip geometry never creates a numerical risk", "[WarpMitigation]")
{
    const auto value = GENERATE(0., -1., std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN());
    WarpRisk risk;
    assess_thermal_strip(risk, *warp_material("PPS"), 100., value, 100., 20., 320.);
    CHECK_FALSE(risk.elevated());
    CHECK_THAT(risk.free_contraction_mm, WithinAbs(0., 1e-12));
    CHECK_THAT(risk.crystallization_contraction_mm, WithinAbs(0., 1e-12));
}

TEST_CASE("Thermal bow follows span squared and inverse solid thickness", "[WarpMitigation]")
{
    const auto material = *warp_material("ASA");
    WarpRisk reference, wider, thicker, isothermal;
    assess_thermal_strip(reference, material, 20., 1., 100., 50., 250.);
    assess_thermal_strip(wider, material, 40., 1., 100., 50., 250.);
    assess_thermal_strip(thicker, material, 20., 2., 100., 50., 250.);
    assess_thermal_strip(isothermal, material, 20., 1., 50., 50., 250.);
    CHECK_THAT(reference.free_bow_mm, WithinAbs(105e-6 * 50. * 400. / 8., 1e-9));
    CHECK_THAT(wider.free_bow_mm, WithinAbs(4. * reference.free_bow_mm, 1e-9));
    CHECK_THAT(thicker.free_bow_mm, WithinAbs(0.5 * reference.free_bow_mm, 1e-9));
    CHECK_THAT(isothermal.free_bow_mm, WithinAbs(0., 1e-9));
    CHECK(isothermal.free_contraction_mm > 0.);
    CHECK_THAT(reference.restrained_energy_j_m2,
        WithinAbs(0.5 * 2.9e9 * std::pow(105e-6 * 55., 2) * 0.001, 1e-9));
}

TEST_CASE("Small ABS contacts do not acquire risk from empty space between them", "[WarpMitigation]")
{
    auto config = DynamicPrintConfig::full_print_config();
    config.set_key_value("filament_type", new ConfigOptionStrings({"ABS"}));
    config.set_key_value("curr_bed_type", new ConfigOptionEnum<BedType>(btPEI));
    config.set_key_value("hot_plate_temp", new ConfigOptionInts({100}));
    config.set_key_value("chamber_temperature", new ConfigOptionInts({50}));
    config.set_key_value("nozzle_temperature", new ConfigOptionInts({250}));
    config.set_key_value("bottom_shell_layers", new ConfigOptionInt(5));
    config.set_key_value("bottom_shell_thickness", new ConfigOptionFloat(1.));
    config.set_key_value("layer_height", new ConfigOptionFloat(0.2));
    Model model;
    auto *object = model.add_object();
    object->add_volume(make_cube(5., 5., 5.));
    auto second = make_cube(5., 5., 5.);
    second.translate(100., 0., 0.);
    object->add_volume(second);
    object->add_instance();
    const auto small = assess_warp_risk(config, model);
    CHECK_FALSE(small.elevated());
    CHECK(small.span_mm < 8.);
    const auto original = config;
    apply_warp_mitigations(config, small);
    CHECK(config.diff(original).empty());
    object->add_volume(make_cube(100., 20., 1.));
    const auto large = assess_warp_risk(config, model);
    CHECK(large.elevated());
    object->instances.front()->set_rotation(Z, 0.7);
    const auto rotated = assess_warp_risk(config, model);
    CHECK_THAT(rotated.free_bow_mm, WithinAbs(large.free_bow_mm, 1e-5));
    config.set_key_value("chamber_temperature", new ConfigOptionInts({100}));
    const auto warm = assess_warp_risk(config, model);
    CHECK_FALSE(warm.elevated());
    CHECK(warm.free_contraction_mm > 0.);
    config.set_key_value("filament_type", new ConfigOptionStrings({"ABS-CF"}));
    const auto filled = assess_warp_risk(config, model);
    CHECK_FALSE(filled.elevated());
    CHECK(filled.unknown_materials.empty());
    CHECK(filled.resin_screened_materials.count("ABS-CF") == 1);
}

TEST_CASE("Bed lift mitigation replaces painted-only adhesion with a continuous brim", "[WarpMitigation][Regression]")
{
    // ValveBoard_case uses ASA/ASA/TPU with painted brim, zero width and
    // a 0.2 mm gap. Keep this regression independent of the private model.
    auto config = DynamicPrintConfig::full_print_config();
    config.set_key_value("brim_type", new ConfigOptionEnum<BrimType>(btPainted));
    config.set_key_value("brim_width", new ConfigOptionFloat(0.));
    config.set_key_value("brim_object_gap", new ConfigOptionFloat(0.2));
    WarpRisk risk;
    risk.materials.insert("ASA");
    risk.span_mm = 120.;
    risk.free_bow_mm = 6.;
    apply_warp_mitigations(config, risk);
    CHECK(config.option<ConfigOptionEnum<BrimType>>("brim_type")->value == btOuterOnly);
    CHECK_THAT(config.opt_float("brim_width"), WithinAbs(10., 1e-6));
    CHECK_THAT(config.opt_float("brim_object_gap"), WithinAbs(0.1, 1e-6));
    const auto applied = config;
    apply_warp_mitigations(config, risk);
    CHECK(config.diff(applied).empty());
}

TEST_CASE("Warp mitigation preserves stronger adhesion and raft choices", "[WarpMitigation]")
{
    auto config = DynamicPrintConfig::full_print_config();
    config.set_key_value("brim_type", new ConfigOptionEnum<BrimType>(btInnerOnly));
    config.set_key_value("brim_width", new ConfigOptionFloat(15.));
    config.set_key_value("brim_object_gap", new ConfigOptionFloat(0.05));
    WarpRisk risk;
    risk.materials.insert("PC");
    risk.free_bow_mm = 1.;
    apply_warp_mitigations(config, risk);
    CHECK(config.option<ConfigOptionEnum<BrimType>>("brim_type")->value == btOuterAndInner);
    CHECK_THAT(config.opt_float("brim_width"), WithinAbs(15., 1e-6));
    CHECK_THAT(config.opt_float("brim_object_gap"), WithinAbs(0.05, 1e-6));

    const bool belt = GENERATE(false, true);
    config.set_key_value("raft_layers", new ConfigOptionInt(belt ? 0 : 3));
    config.set_key_value("belt_printer", new ConfigOptionBool(belt));
    config.set_key_value("brim_type", new ConfigOptionEnum<BrimType>(btNoBrim));
    apply_warp_mitigations(config, risk);
    CHECK(config.option<ConfigOptionEnum<BrimType>>("brim_type")->value == btNoBrim);
}

TEST_CASE("Absent warp risk leaves settings untouched", "[WarpMitigation]")
{
    auto config = DynamicPrintConfig::full_print_config();
    const auto original = config;
    apply_warp_mitigations(config, {});
    CHECK(config.diff(original).empty());
    CHECK_FALSE(config.opt_bool("automatic_warp_mitigation"));
    Model empty;
    CHECK_FALSE(assess_warp_risk(config, empty).elevated());
}

TEST_CASE("Warp screening follows model and feature material assignments", "[WarpMitigation]")
{
    auto config = DynamicPrintConfig::full_print_config();
    config.set_key_value("filament_type", new ConfigOptionStrings({"PLA", "ASA", "TPU"}));
    Model model;
    auto *object = model.add_object();
    auto *volume = object->add_volume(make_cube(80., 70., 2.));
    object->add_instance();
    volume->config.set_key_value("extruder", new ConfigOptionInt(1));
    const auto initial = assess_warp_risk(config, model);
    CHECK(initial.materials.count("PLA") == 1);
    CHECK(initial.materials.count("ASA") == 0);

    const auto key = GENERATE("extruder", "outer_wall_filament_id",
        "inner_wall_filament_id", "sparse_infill_filament_id",
        "internal_solid_filament_id", "top_surface_filament_id", "bottom_surface_filament_id");
    volume->config.set_key_value(key, new ConfigOptionInt(2));
    const auto risk = assess_warp_risk(config, model);
    CHECK(risk.materials.count("ASA") == 1);
    CHECK(risk.span_mm > 100.);
}

TEST_CASE("A supplied bed-lift project receives an adhesion recommendation", "[WarpMitigation][.LocalProject]")
{
    // Opt-in local regression: do not ship private project geometry in the suite.
    const char *path = std::getenv("ORCA_WARP_REGRESSION_3MF");
    REQUIRE(path != nullptr);
    ScopedSlic3rTemporaryDir backup{"orca_warp_regression"};
    boost::filesystem::create_directories(boost::filesystem::path(backup.string()) / "Metadata");
    Model model;
    model.set_backup_path(backup.string());
    // Import into a sparse config, as the project loader does. Static defaults
    // do not supply all generic enum lookup maps needed by JSON deserialization.
    DynamicPrintConfig config;
    ConfigSubstitutionContext substitutions{ForwardCompatibilitySubstitutionRule::Enable};
    struct Resources {
        PlateDataPtrs plates;
        std::vector<Preset*> presets;
        ~Resources() { release_PlateData_list(plates); for (auto *preset : presets) delete preset; }
    } resources;
    bool bbs = false, orca = false;
    Semver version;
    REQUIRE(load_bbs_3mf(path, &config, &substitutions, &model, &resources.plates,
        &resources.presets, &bbs, &orca, &version, nullptr,
        LoadStrategy::LoadModel | LoadStrategy::LoadConfig));
    const auto risk = assess_warp_risk(config, model);
    REQUIRE(risk.elevated());
    INFO("Recommended brim width: " << risk.brim_mm());
    apply_warp_mitigations(config, risk);
    CHECK(config.option<ConfigOptionEnum<BrimType>>("brim_type")->value == btOuterOnly);
    CHECK(config.opt_float("brim_width") >= risk.brim_mm());
    CHECK(config.opt_float("brim_object_gap") <= 0.1);
}
