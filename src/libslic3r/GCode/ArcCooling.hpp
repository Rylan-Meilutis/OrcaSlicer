#pragma once

#include <algorithm>
#include <cmath>
#include <optional>
#include <string_view>

namespace Slic3r {

// Version 5 reference data, not measured properties of the user's spool.
// Sources, validity and model assumptions: docs/arc-material-cooling.md.
struct ArcThermalMaterial {
    double density;       // g/cm^3
    double specific_heat; // J/(kg K), reference value, not a melt measurement
    double conductivity;  // W/(m K)
    double glass_transition; // degrees C; proxy, not a load-bearing criterion
    // Semi-crystalline references use a solidification scheduling target instead
    // of Tg (PP's Tg is below room temperature), and include latent heat.
    double cooling_target = 0.;
    double crystallization_heat = 0.; // J/kg; conservative fully crystalline reference
    double crystallization_temperature = 0.; // C; lumped heat-release proxy
    bool matrix_proxy = false; // Composite formulation/orientation is unknown.
};

inline std::optional<ArcThermalMaterial> arc_thermal_material(std::string_view type, double temperature)
{
    // Exact matches exclude unspecified blends and PET (not PETG). PPS-CF is
    // explicitly identified below as a matrix proxy, not a composite measurement.
    if (type == "PLA")
        return ArcThermalMaterial{1.24, 1800., 0.2, 62.8};
    // NETZSCH generic ranges: representative midpoints, using the rigid phase
    // transition for ABS/ASA. Do not extend these to blends or filled grades.
    if (type == "ABS")
        return ArcThermalMaterial{1.05, 1470., 0.175, 100.};
    if (type == "ASA")
        return ArcThermalMaterial{1.055, 1350., 0.18, 100.};
    if (type == "PC")
        return ArcThermalMaterial{1.22, 1335., 0.20, 145.};
    if (type == "PMMA")
        return ArcThermalMaterial{1.17, 1460., 0.205, 105.};
    if (type == "PSU")
        return ArcThermalMaterial{1.245, 1370., 0.15, 187.};
    // The filament catalog uses PES for polyethersulfone (PESU).
    if (type == "PES" || type == "PESU")
        return ArcThermalMaterial{1.37, 1370., 0.18, 228.};
    if (type == "PP")
        return ArcThermalMaterial{0.905, 1800., 0.21, 0., 100., 208000., 120.};
    if (type == "PPS" || type == "PPS-CF")
        return ArcThermalMaterial{1.35, 1830., 0.30, 90., 90., 150400., 205., type == "PPS-CF"};
    // HIPS-specific reference: FSRI solid-state cp (40 C) and k (45 C),
    // with Impact 3450 density and the unfilled filament's measured matrix Tg.
    // These are fixed reference values, not extrapolated melt measurements.
    if (type == "HIPS")
        return ArcThermalMaterial{1.04, 1351., 0.128, 101.1};
    if (type == "PETG") {
        // Eastar 6763 DSC data at 100/150/200/250 C. Clamp outside the data range.
        constexpr double cp[] = {1760., 1880., 1970., 2050.};
        if (!std::isfinite(temperature))
            return std::nullopt;
        const double x = (std::clamp(temperature, 100., 250.) - 100.) / 50.;
        const int i = std::min(2, int(x));
        return ArcThermalMaterial{1.27, cp[i] + (cp[i + 1] - cp[i]) * (x - i), 0.21, 80.};
    }
    return std::nullopt;
}

struct ArcCoolingEstimate {
    double path_time;
    double cover_speed;
};

inline std::optional<ArcCoolingEstimate> estimate_arc_cooling(
    std::string_view type, double density, double temperature, double chamber_temperature,
    double fan_percent, double bead_diameter)
{
    const auto material = arc_thermal_material(type, temperature);
    if (!material || !std::isfinite(density) || !std::isfinite(temperature) ||
        !std::isfinite(chamber_temperature) || !std::isfinite(fan_percent) ||
        !std::isfinite(bead_diameter) || density < 0. || bead_diameter <= 0. ||
        bead_diameter > 2. || fan_percent <= 0. || fan_percent > 100.)
        return std::nullopt;
    const double ambient = std::max(25., chamber_temperature);
    const double target = material->cooling_target > 0. ? material->cooling_target : material->glass_transition;
    if (ambient >= target - 5. || temperature <= target)
        return std::nullopt;
    const double rho = 1000. * (density > 0. ? density : material->density);
    if (rho < 500. || rho > 2000.)
        return std::nullopt; // Do not extrapolate to foams or heavily filled material.
    const double d = bead_diameter * 0.001;
    // Uncalibrated convection assumption, not a measured fan PWM -> airflow map.
    const double h = 20. + 80. * fan_percent / 100.;
    const double convection = rho * material->specific_heat * d / (4. * h) *
        std::log((temperature - ambient) / (target - ambient));
    const double diffusion = rho * material->specific_heat * d * d / (16. * material->conductivity);
    // Heat released during crystallization must also leave the strand. This is
    // a lumped scheduling allowance, not a prediction of crystallization kinetics.
    const double latent = material->crystallization_heat > 0. && temperature > material->crystallization_temperature ?
        rho * material->crystallization_heat * d /
            (4. * h * (material->crystallization_temperature - ambient)) : 0.;
    const double time = convection + diffusion + latent;
    if (!std::isfinite(time) || time > 30.)
        return std::nullopt; // Out of model range; do not silently claim adequate cooling.
    const double bounded_time = std::max(1., time);
    // A 15 mm reference traverse is a scheduling heuristic, not a material property.
    return ArcCoolingEstimate{bounded_time, std::clamp(15. / bounded_time, 0.5, 15.)};
}

} // namespace Slic3r
