#pragma once

#include "ArcCooling.hpp"

namespace Slic3r {

// A conservative starting point, not a melt-viscosity or hotend-response model.
// The user's spool temperature range takes precedence over generic material data.
inline int feature_temperature(bool fast_hotend, bool automatic, std::string_view material,
                               int normal, int override_temperature, int minimum, int maximum,
                               bool gentle = false)
{
    if (!fast_hotend || normal <= 0)
        return normal;
    if (override_temperature > 0)
        return override_temperature;
    if (!automatic || minimum <= 0 || maximum < minimum || normal < minimum || normal > maximum)
        return normal;
    const auto properties = arc_thermal_material(material, normal);
    if (!properties)
        return normal;
    const double target = properties->cooling_target > 0. ? properties->cooling_target : properties->glass_transition;
    const int reduction = int(std::lround(std::clamp((normal - target) * (gentle ? 0.02 : 0.04), 3., 10.)));
    return std::max(minimum, normal - reduction);
}

} // namespace Slic3r
