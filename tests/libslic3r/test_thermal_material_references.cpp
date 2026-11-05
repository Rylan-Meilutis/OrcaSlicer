#include "libslic3r/ThermalMaterialReferences.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <limits>
#include <set>
#include <string_view>

using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinRel;

TEST_CASE("Named ASA transport references retain measured domain and SI units", "[ThermalMaterialReferences]")
{
    const auto* material = thermal_material_reference("stratasys-asa-natural-2025");
    REQUIRE(material);
    const auto k = reference_value(*material, ReferenceQuantity::Conductivity, 30.);
    const auto capacity = reference_value(*material, ReferenceQuantity::VolumetricHeatCapacity, 30.);
    REQUIRE(k);
    REQUIRE(capacity);
    CHECK_THAT(*k, WithinRel(.1642));
    CHECK_THAT(*capacity, WithinRel(.1642 / .096e-6));
    CHECK_THAT(*reference_value(*material, ReferenceQuantity::Conductivity, 45.), WithinRel((.1642 + .1622) / 2.));
    CHECK_FALSE(reference_value(*material, ReferenceQuantity::Conductivity, 110.));
    CHECK_FALSE(reference_value(*material, ReferenceQuantity::Conductivity, -1.));
    CHECK_FALSE(reference_value(*material, ReferenceQuantity::Conductivity, std::numeric_limits<double>::quiet_NaN()));
    CHECK_FALSE(reference_value(*material, ReferenceQuantity::YoungModulus, 30.));
}

TEST_CASE("Unspecified temperature domains cannot silently become material laws", "[ThermalMaterialReferences]")
{
    REQUIRE(thermal_material_reference("ultimaker-tpu95a-2017"));
    CHECK_FALSE(reference_value(*thermal_material_reference("ultimaker-tpu95a-2017"), ReferenceQuantity::YoungModulus, 23.));
    CHECK_FALSE(thermal_material_reference("TPU"));
    CHECK_FALSE(thermal_material_reference("ASA"));
    std::set<std::string_view> ids;
    for (const auto& material : thermal_material_references()) {
        CHECK(ids.insert(material.id).second);
        CHECK_FALSE(material.source_url.empty());
        for (const auto& property : material.properties) {
            CHECK_FALSE(property.notes.empty());
            if (property.constant_si) CHECK(property.points.empty());
        }
    }
    const auto* asa = thermal_material_reference("luran-s778t-tensile-2026");
    REQUIRE(asa);
    CHECK_THAT(*reference_value(*asa, ReferenceQuantity::YoungModulus, 23.), WithinRel(2206e6));
    CHECK_FALSE(reference_value(*asa, ReferenceQuantity::YoungModulus, 90.));
}
