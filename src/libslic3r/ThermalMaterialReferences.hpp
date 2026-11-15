#pragma once

#include <optional>
#include <string_view>
#include <vector>

namespace Slic3r::Thermoelastic {
enum class ReferenceQuantity { Conductivity, SpecificHeatCapacity, VolumetricHeatCapacity, Density, YoungModulus, LinearExpansion };
enum class ReferenceEvidence { Measured, DerivedFromMeasurements, PublishedModelAssumption };
struct ReferencePoint {
    double temperature_c;
    double value_si;
};
struct ReferenceProperty {
    ReferenceQuantity quantity;
    ReferenceEvidence evidence;
    std::vector<ReferencePoint> points;
    // A constant with no specified measurement domain is information for an
    // explicit scenario choice, NOT a law evaluable at arbitrary temperature.
    std::optional<double> constant_si;
    std::string_view notes;
};
struct ThermalMaterialReference {
    std::string_view id, name, source_url;
    std::vector<ReferenceProperty> properties;
};
// Named references only. No match by generic filament name, no mixing grades,
// no inferred Poisson ratio, deposition/reference temperature or interface law.
const std::vector<ThermalMaterialReference>& thermal_material_references();
const ThermalMaterialReference* thermal_material_reference(std::string_view id);
// Linear interpolation is an explicit numerical convention between tabulated
// measurements. No extrapolation, including unspecified-domain constants.
std::optional<double> reference_value(const ThermalMaterialReference& reference,
                                      ReferenceQuantity quantity, double temperature_c);
} // namespace Slic3r::Thermoelastic
