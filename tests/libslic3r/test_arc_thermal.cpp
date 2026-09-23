#include <catch2/catch_all.hpp>
#include "libslic3r/GCode/ArcCooling.hpp"
#include <limits>

using namespace Slic3r;

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
    for (const char *type : {"PLA-CF", "PET", "TPU", "ABS-GF", "ASA-AERO", "PC-ABS", "PA", "PP", "unknown"})
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
    const auto type = GENERATE("ABS", "ASA", "PC");
    CAPTURE(type);
    const auto material = arc_thermal_material(type, 260.);
    REQUIRE(material);
    const auto automatic = estimate_arc_cooling(type, 0., 260., 40., 50., 0.45);
    const auto explicit_density = estimate_arc_cooling(type, material->density, 260., 40., 50., 0.45);
    REQUIRE(automatic);
    REQUIRE(explicit_density);
    CHECK_THAT(automatic->path_time, Catch::Matchers::WithinAbs(explicit_density->path_time, 1e-9));
    const auto denser = estimate_arc_cooling(type, material->density * 1.1, 260., 40., 50., 0.45);
    REQUIRE(denser);
    CHECK(denser->path_time > automatic->path_time);
    CHECK_FALSE(estimate_arc_cooling(type, 0., 260., material->glass_transition - 5., 50., 0.45));
    CHECK_FALSE(estimate_arc_cooling(type, 0., 260., 40., 0., 0.45));
}

TEST_CASE("PETG reference heat capacity interpolates only within its source range", "[ArcThermal]")
{
    REQUIRE(arc_thermal_material("PETG", 175.));
    CHECK_THAT(arc_thermal_material("PETG", 175.)->specific_heat, Catch::Matchers::WithinAbs(1925., 1e-9));
    CHECK_THAT(arc_thermal_material("PETG", 300.)->specific_heat, Catch::Matchers::WithinAbs(2050., 1e-9));
    CHECK_THAT(arc_thermal_material("PETG", 50.)->specific_heat, Catch::Matchers::WithinAbs(1760., 1e-9));
    REQUIRE(estimate_arc_cooling("PETG", 1.27, 240., 25., 70., 0.45));
}
