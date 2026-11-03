#include "ThermalMaterialReferences.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <string_view>
#include <vector>

namespace Slic3r::Thermoelastic {
const std::vector<ThermalMaterialReference>& thermal_material_references()
{
    using Q = ReferenceQuantity;
    using E = ReferenceEvidence;
    static const std::vector<ThermalMaterialReference> references{
        {"stratasys-asa-natural-2025", "Stratasys ASA natural: transport measurements",
         "https://www.stratasys.com/siteassets/materials/materials-catalog/fdm-materials/asa/mds_fdm_asa_0925a.pdf",
         {{Q::Conductivity, E::Measured, {{0., .1685}, {30., .1642}, {60., .1622}, {90., .1563}}, {},
           "ASTM E1952, natural ASA. Values are typical, not population bounds."},
          {Q::VolumetricHeatCapacity, E::DerivedFromMeasurements,
           {{0., .1685 / (.108e-6)}, {30., .1642 / (.096e-6)},
            {60., .1622 / (.087e-6)}, {90., .1563 / (.077e-6)}}, {},
           "rho*cp = conductivity / diffusivity at the same temperature; no assumed density correction."},
          {Q::Density, E::Measured, {{23., 1080.}}, {}, "ASTM D792 specific gravity, SI density conversion."}}},
        {"luran-s778t-tensile-2026", "BASF Luran S778T ASA: tensile measurements",
         "https://pmc.ncbi.nlm.nih.gov/articles/PMC12897946/",
         {{Q::YoungModulus, E::Measured, {{-10., 2487e6}, {23., 2206e6}, {65., 2102e6}}, {},
           "Injection-molded, strain rate 0.01/s; means, not printed orthotropy or relaxation. Reported SD: 196,103,156 MPa."}}},
        {"ultimaker-tpu95a-2017", "UltiMaker TPU 95A: printed reference properties",
         "https://ase.au.dk/fileadmin/www.ase.au.dk/Filer/Laboratorier_og_vaerksteder/MakerLab/Ultimaker-TDS_TPU_95A_v3.010_ENG.pdf",
         {{Q::YoungModulus, E::Measured, {}, 26e6,
           "ASTM D638 printed bars; this table does not establish a temperature-dependent modulus law."},
          {Q::LinearExpansion, E::Measured, {}, 100e-6,
           "ASTM E693 typical CTE; temperature interval unspecified. Not the TPU soft-segment glass transition."},
          {Q::Density, E::Measured, {}, 1220., "Specific gravity 1.22; temperature interval unspecified."}}},
        {"polyflex95a-delbart-2024", "PolyFlex TPU 95A: published actuator-model assumptions",
         "https://doi.org/10.1016/j.compositesa.2024.108381",
         {{Q::YoungModulus, E::PublishedModelAssumption, {}, 8.4e6,
           "Initial modulus of a nonlinear elastic model, not a complete linear constitutive calibration."},
          {Q::Conductivity, E::PublishedModelAssumption, {}, .15,
           "Table 3: 0.00015 W/(mm K), converted to W/(m K); borrowed constant, no measured temperature domain."},
          {Q::SpecificHeatCapacity, E::PublishedModelAssumption, {}, 1210.,
           "Table 3: 1.21 J/(g K), converted to J/(kg K). Density is not supplied by this reference entry."},
          {Q::LinearExpansion, E::PublishedModelAssumption, {}, 110e-6,
           "Table 3 constant assumption. Study disregards printing residual stresses; not a deposition calibration."}}}
    };
    return references;
}

const ThermalMaterialReference* thermal_material_reference(std::string_view id)
{
    const auto& references = thermal_material_references();
    const auto found = std::find_if(references.begin(), references.end(),
        [id](const auto& reference) { return reference.id == id; });
    return found == references.end() ? nullptr : &*found;
}

std::optional<double> reference_value(const ThermalMaterialReference& reference,
                                      ReferenceQuantity quantity, double temperature_c)
{
    if (!std::isfinite(temperature_c)) return {};
    const auto property = std::find_if(reference.properties.begin(), reference.properties.end(),
        [quantity](const auto& item) { return item.quantity == quantity; });
    if (property == reference.properties.end() || property->points.empty()) return {};
    const auto& points = property->points;
    double previous = -std::numeric_limits<double>::infinity();
    for (const auto& point : points) {
        if (!std::isfinite(point.temperature_c) || !std::isfinite(point.value_si) ||
            point.temperature_c <= previous || point.value_si <= 0.) return {};
        previous = point.temperature_c;
    }
    if (temperature_c < points.front().temperature_c || temperature_c > points.back().temperature_c) return {};
    const auto upper = std::lower_bound(points.begin(), points.end(), temperature_c,
        [](const auto& point, double temperature) { return point.temperature_c < temperature; });
    if (upper == points.begin() || upper->temperature_c == temperature_c) return upper->value_si;
    const auto& lower = *(upper - 1);
    const double fraction = (temperature_c - lower.temperature_c) / (upper->temperature_c - lower.temperature_c);
    return lower.value_si + fraction * (upper->value_si - lower.value_si);
}
} // namespace Slic3r::Thermoelastic
