#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/Config.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Semver.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/WarpMitigation.hpp"
#include "libslic3r/Format/bbs_3mf.hpp"
#include "libslic3r/Preset.hpp"
#include "../test_utils.hpp"
#include <cstdlib>
#include <vector>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

using namespace Slic3r;
using Catch::Matchers::WithinAbs;

TEST_CASE("Warp screening recognizes susceptible material families", "[WarpMitigation]")
{
    const auto material = GENERATE("ABS", "ASA", "PC", "PA12", "PP", "PPS-CF", "HIPS");
    CHECK(warp_prone_material(material));
    CHECK_FALSE(warp_prone_material("PLA"));
    CHECK_FALSE(warp_prone_material("TPU"));
    CHECK_FALSE(warp_prone_material("unknown"));
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
    CHECK_FALSE(assess_warp_risk(config, model).elevated());

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
