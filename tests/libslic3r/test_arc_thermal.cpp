#include <catch2/catch_all.hpp>
#include "libslic3r/GCode/ArcCooling.hpp"
#include "libslic3r/GCode/FeatureTemperature.hpp"
#include <limits>

using namespace Slic3r;

TEST_CASE("Feature temperature estimates respect opt in and spool limits", "[ArcThermal][FeatureTemperature]")
{
    CHECK(feature_temperature(false, true, "PLA", 220, 190, 190, 240) == 220);
    CHECK(feature_temperature(true, false, "PLA", 220, 0, 190, 240) == 220);
    CHECK(feature_temperature(true, true, "Unknown", 220, 0, 190, 240) == 220);
    CHECK(feature_temperature(true, true, "PLA", 220, 0, 0, 0) == 220);
    CHECK(feature_temperature(true, true, "PLA", 220, 205, 190, 240) == 205);
    CHECK(feature_temperature(true, true, "PLA", 220, 0, 218, 240) == 218);
    for (const auto material : {"PLA", "PETG", "HIPS", "PP", "PPS-CF"}) {
        const int normal = std::string_view(material) == "PPS-CF" ? 320 : 240;
        const int target = feature_temperature(true, true, material, normal, 200, 180, 350);
        CHECK(target == 200);
        const int automatic = feature_temperature(true, true, material, normal, 0, 180, 350);
        CHECK(automatic < normal);
        CHECK(automatic >= normal - 10);
        CHECK(feature_temperature(true, true, material, normal, 0, 180, 350, true) >= automatic);
    }
}

TEST_CASE("Arc thermal estimates respond to heat storage and cooling", "[ArcThermal]")
{
    const auto base = estimate_arc_cooling("PLA", 1.24, 215., 25., 100., 0.45);
    REQUIRE(base.has_value());
    for (const auto &other : {
        estimate_arc_cooling("PLA", 1.4, 215., 25., 100., 0.45),
        estimate_arc_cooling("PLA", 1.24, 230., 25., 100., 0.45),
        estimate_arc_cooling("PLA", 1.24, 215., 40., 100., 0.45),
        estimate_arc_cooling("PLA", 1.24, 215., 25., 50., 0.45),
        estimate_arc_cooling("PLA", 1.24, 215., 25., 100., 0.6)}) {
        REQUIRE(other.has_value());
        CHECK(other->path_time > base->path_time);
        CHECK(other->cover_speed < base->cover_speed);
    }
    const auto fallback_density = estimate_arc_cooling("PLA", 0., 215., 25., 100., 0.45);
    REQUIRE(fallback_density.has_value());
    CHECK_THAT(fallback_density->path_time, Catch::Matchers::WithinAbs(base->path_time, 1e-9));
}

TEST_CASE("Arc thermal estimates decline unsupported materials and conditions", "[ArcThermal]")
{
    for (const char *type : {"PLA-CF", "PET", "TPU", "ABS-GF", "ASA-AERO", "PC-ABS", "HIPS-CF", "HIPS-FR", "PA", "PP-GF", "unknown"})
        CHECK_FALSE(estimate_arc_cooling(type, 1.24, 215., 25., 100., 0.45));
    CHECK_FALSE(estimate_arc_cooling("PLA", 1.24, 215., 70., 100., 0.45));
    CHECK_FALSE(estimate_arc_cooling("PLA", 1.24, 215., 25., 0., 0.45));
    CHECK_FALSE(estimate_arc_cooling("PLA", 1.24, 215., 25., 100., 0.));
    CHECK_FALSE(estimate_arc_cooling("PLA", -1., 215., 25., 100., 0.45));
    CHECK_FALSE(estimate_arc_cooling("PLA", 1.24, std::numeric_limits<double>::quiet_NaN(), 25., 100., 0.45));
    CHECK_FALSE(estimate_arc_cooling("PETG", 1.27, std::numeric_limits<double>::infinity(), 25., 100., 0.45));
}

TEST_CASE("Amorphous material references resolve independently and honor profile density", "[ArcThermal]")
{
    const auto type = GENERATE("ABS", "ASA", "PC", "HIPS", "PMMA", "PSU", "PES", "PESU");
    CAPTURE(type);
    const auto material = arc_thermal_material(type, 260.);
    REQUIRE(material);
    const double temperature = std::max(260., material->glass_transition + 150.);
    const auto automatic = estimate_arc_cooling(type, 0., temperature, 40., 50., 0.45);
    const auto explicit_density = estimate_arc_cooling(type, material->density, temperature, 40., 50., 0.45);
    REQUIRE(automatic);
    REQUIRE(explicit_density);
    CHECK_THAT(automatic->path_time, Catch::Matchers::WithinAbs(explicit_density->path_time, 1e-9));
    const auto denser = estimate_arc_cooling(type, material->density * 1.1, temperature, 40., 50., 0.45);
    REQUIRE(denser);
    CHECK(denser->path_time > automatic->path_time);
    CHECK_FALSE(estimate_arc_cooling(type, 0., 260., material->glass_transition - 5., 50., 0.45));
    CHECK_FALSE(estimate_arc_cooling(type, 0., 260., 40., 0., 0.45));
}

TEST_CASE("HIPS arc estimates respond to fan cooling at printing temperature", "[ArcThermal]")
{
    const auto full_fan = estimate_arc_cooling("HIPS", 1.04, 240., 40., 100., 0.45);
    const auto reduced_fan = estimate_arc_cooling("HIPS", 1.04, 240., 40., 50., 0.45);
    REQUIRE(full_fan);
    REQUIRE(reduced_fan);
    CHECK(std::isfinite(full_fan->path_time));
    CHECK(full_fan->path_time >= 1.);
    CHECK(full_fan->path_time <= 30.);
    CHECK(reduced_fan->path_time > full_fan->path_time);
    CHECK(reduced_fan->cover_speed < full_fan->cover_speed);
}

TEST_CASE("PETG reference heat capacity interpolates only within its source range", "[ArcThermal]")
{
    REQUIRE(arc_thermal_material("PETG", 175.));
    CHECK_THAT(arc_thermal_material("PETG", 175.)->specific_heat, Catch::Matchers::WithinAbs(1925., 1e-9));
    CHECK_THAT(arc_thermal_material("PETG", 300.)->specific_heat, Catch::Matchers::WithinAbs(2050., 1e-9));
    CHECK_THAT(arc_thermal_material("PETG", 50.)->specific_heat, Catch::Matchers::WithinAbs(1760., 1e-9));
    REQUIRE(estimate_arc_cooling("PETG", 1.27, 240., 25., 70., 0.45));
}

TEST_CASE("Semi crystalline arc references account for heat release and retain bounded estimates", "[ArcThermal]")
{
    const auto type = GENERATE("PP", "PPS", "PPS-CF");
    const double nozzle = std::string_view(type) == "PP" ? 240. : 320.;
    const auto reference = arc_thermal_material(type, nozzle);
    REQUIRE(reference);
    CHECK(reference->crystallization_heat > 0.);
    const auto cool = estimate_arc_cooling(type, 0., nozzle, 40., 100., 0.45);
    const auto warm = estimate_arc_cooling(type, 0., nozzle, 60., 100., 0.45);
    REQUIRE(cool);
    REQUIRE(warm);
    CHECK(warm->path_time > cool->path_time);
    CHECK(cool->path_time <= 30.);
    CHECK(cool->cover_speed >= 0.5);
    CHECK_FALSE(estimate_arc_cooling(type, 0., nozzle, reference->cooling_target, 100., 0.45));
    CHECK_FALSE(estimate_arc_cooling(type, 0., nozzle, 40., 0., 0.45));
    CHECK(reference->matrix_proxy == (std::string_view(type) == "PPS-CF"));
}
