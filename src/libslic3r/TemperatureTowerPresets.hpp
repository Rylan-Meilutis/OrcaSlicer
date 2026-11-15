#pragma once

#include <array>

namespace Slic3r {

struct TemperatureTowerPreset {
    const char *material;
    int start;
    int end;
};

// Editable calibration starting points, not universal material or printer limits.
// Keep labels and ranges together so adding a material cannot shift the GUI mapping.
inline constexpr std::array<TemperatureTowerPreset, 12> temperature_tower_presets {{
    {"PLA", 230, 190},
    {"ABS/ASA", 270, 230},
    {"PETG", 250, 230},
    {"PCTG", 280, 240},
    {"TPU", 240, 210},
    {"PA-CF", 320, 280},
    {"PET-CF", 320, 280},
    // Fillamentum HIPS Extrafill: https://fillamentum.com/blogs/news/4-new-hips-colors/
    {"HIPS", 250, 230},
    // Fiberlogy PP: https://fiberlogy.com/product/filament-pp/
    {"PP", 250, 220},
    // Bambu filament guide: https://bambulab.com/it/filament-guide
    {"PC", 290, 260},
    {"PPA-CF", 310, 280},
    // PPS grades vary; this starts around the bundled Generic PPS profile's 320 C.
    // PPS-CF range: https://bambulab.cn/zh-cn/filament/pps-cf
    {"PPS/PPS-CF", 340, 310},
}};

} // namespace Slic3r
