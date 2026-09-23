#pragma once

#include <algorithm>
#include <cmath>
#include <optional>
#include <string_view>

namespace Slic3r {

// Version 2 reference data, not measured properties of the user's spool.
// Sources, validity and model assumptions: docs/arc-material-cooling.md.
struct ArcThermalMaterial {
    double density;       // g/cm^3
    double specific_heat; // J/(kg K), reference value, not a melt measurement
    double conductivity;  // W/(m K)
    double glass_transition; // degrees C; proxy, not a load-bearing criterion
};

inline std::optional<ArcThermalMaterial> arc_thermal_material(std::string_view type, double temperature)
{
    // Exact matches deliberately exclude filled/foamed blends and PET (not PETG).
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
    if (ambient >= material->glass_transition - 5. || temperature <= material->glass_transition)
        return std::nullopt;
    const double rho = 1000. * (density > 0. ? density : material->density);
    if (rho < 500. || rho > 2000.)
        return std::nullopt; // Do not extrapolate to foams or heavily filled material.
    const double d = bead_diameter * 0.001;
    // Uncalibrated convection assumption, not a measured fan PWM -> airflow map.
    const double h = 20. + 80. * fan_percent / 100.;
    const double convection = rho * material->specific_heat * d / (4. * h) *
        std::log((temperature - ambient) / (material->glass_transition - ambient));
    const double diffusion = rho * material->specific_heat * d * d / (16. * material->conductivity);
    const double time = convection + diffusion;
    if (!std::isfinite(time) || time > 30.)
        return std::nullopt; // Out of model range; do not silently claim adequate cooling.
    const double bounded_time = std::max(1., time);
    // A 15 mm reference traverse is a scheduling heuristic, not a material property.
    return ArcCoolingEstimate{bounded_time, std::clamp(15. / bounded_time, 0.5, 15.)};
}

} // namespace Slic3r
