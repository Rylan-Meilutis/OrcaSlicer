#include "SeamCalibration.hpp"
#include "libslic3r/Config.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/libslic3r.h"
#include <algorithm>
#include <utility>
#include <vector>
#include "Model.hpp"
#include "PrintConfig.hpp"
#include "TriangleMesh.hpp"
#include <cmath>
#include <stdexcept>
#include <sstream>

namespace Slic3r {

std::vector<double> seam_calibration_values(const SeamCalibrationParams &p)
{
    if (p.count < 1 || p.count > 8 || !std::isfinite(p.nozzle) || p.nozzle < 0.15 || p.nozzle > 1.2 ||
        !std::isfinite(p.start) || !std::isfinite(p.step) || int(p.kind) < 0 || int(p.kind) > 4)
        throw std::invalid_argument("Invalid seam calibration range or nozzle diameter.");
    std::vector<double> values;
    for (int i = 0; i < p.count; ++i) {
        const double value = p.start + i * p.step;
        double low = 0., high = 100.;
        switch (p.kind) {
        case SeamCalibrationKind::Length: low = 1.; high = 40.; break;
        case SeamCalibrationKind::StartHeight: high = 90.; break;
        case SeamCalibrationKind::Flow: low = 0.5; high = 1.2; break;
        case SeamCalibrationKind::Speed: low = 5.; high = 150.; break;
        case SeamCalibrationKind::InnerEntry: break;
        }
        if (!std::isfinite(value) || value < low || value > high)
            throw std::invalid_argument("Seam calibration values are outside the supported range.");
        values.push_back(value);
    }
    return values;
}

// Closed annular loft with matching outer/inner rings. The triangular sideways
// displacement repeats in every band, so each setting sees identical slopes.
static TriangleMesh seam_coupon(int shape, int bands, double angle, double wall)
{
    std::vector<Vec2d> contour;
    if (shape == 0)
        contour = {{-10,-10}, {10,-10}, {10,10}, {-10,10}};
    else if (shape == 3)
        contour = {{-10,-10}, {10,-10}, {10,0}, {0,0}, {0,10}, {-10,10}};
    else for (int i = 0; i < 64; ++i) {
        const double a = 2. * PI * i / 64.;
        const double exponent = shape == 1 ? 0.35 : 1.;
        contour.emplace_back(10. * std::copysign(std::pow(std::abs(std::cos(a)), exponent), std::cos(a)),
                             10. * std::copysign(std::pow(std::abs(std::sin(a)), exponent), std::sin(a)));
    }
    std::vector<Vec3f> vertices;
    std::vector<Vec3i32> faces;
    const int n = int(contour.size());
    const double inset = 1. - wall / 10.;
    std::vector<Vec2d> inside;
    if (shape == 3)
        inside = {{-10+wall,-10+wall}, {10-wall,-10+wall}, {10-wall,-wall},
                  {-wall,-wall}, {-wall,10-wall}, {-10+wall,10-wall}};
    else for (const Vec2d &v : contour) inside.push_back(v * inset);
    const int rings = bands * 4 + 1;
    for (int r = 0; r < rings; ++r) {
        const int phase = r % 4;
        const double shift = (phase == 2 ? 2. : 0.) * std::tan(angle * PI / 180.);
        for (int side = 0; side < 2; ++side)
            for (const Vec2d &v : (side ? inside : contour))
                vertices.emplace_back(float(v.x() + shift), float(v.y()), float(r * 2.));
    }
    for (int r = 0; r < rings - 1; ++r)
        for (int side = 0; side < 2; ++side)
            for (int i = 0; i < n; ++i) {
                const int a = r * 2 * n + side * n + i;
                const int b = r * 2 * n + side * n + (i + 1) % n;
                const int c = b + 2 * n, d = a + 2 * n;
                if (side) { faces.emplace_back(a,c,b); faces.emplace_back(a,d,c); }
                else { faces.emplace_back(a,b,c); faces.emplace_back(a,c,d); }
            }
    for (int i = 0; i < n; ++i) {
        const int j = (i + 1) % n, top = (rings - 1) * 2 * n;
        faces.emplace_back(i,n+j,j); faces.emplace_back(i,n+i,n+j);
        faces.emplace_back(top+i,top+j,top+n+j); faces.emplace_back(top+i,top+n+j,top+n+i);
    }
    return TriangleMesh(std::move(vertices), std::move(faces));
}

void make_seam_calibration(Model &model, const DynamicPrintConfig &process, const SeamCalibrationParams &p)
{
    const auto values = seam_calibration_values(p);
    const bool inner = p.kind == SeamCalibrationKind::InnerEntry;
    const char *shapes[] = {"Sharp square", "Rounded square", "Circle / hole", "Concave corner"};
    const char *keys[] = {"seam_slope_min_length", "seam_slope_start_height", "scarf_joint_flow_ratio", "scarf_joint_speed"};
    for (int row = 0; row < (inner ? 2 : 3); ++row)
        for (int shape = 0; shape < 4; ++shape) {
            auto *object = model.add_object();
            const double angle = inner ? 0. : row * 30.;
            std::ostringstream name;
            name << shapes[shape] << (inner ? (row ? " - L-entry ON" : " - L-entry OFF") : " - scarf sweep");
            if (!inner) {
                name << " - " << angle << " deg from vertical - control";
                for (double value : values) name << " / " << value;
            }
            object->name = name.str();
            auto mesh = seam_coupon(shape, inner ? 1 : p.count + 1, angle, std::clamp(6. * p.nozzle, 2.4, 4.));
            // Front witness ticks identify bands without putting labels on the
            // rear seam surface. No ticks means the conventional control band.
            for (int band = inner ? 0 : 1; band < (inner ? 1 : p.count + 1); ++band) {
                const int ticks = inner ? row + 1 : band;
                for (int tick = 0; tick < ticks; ++tick) {
                    auto mark = make_cube(0.6, 1.4, 0.8);
                    mark.translate((tick - (ticks - 1) / 2.) * 1.1 - 0.3, -10.4, band * 8. + 0.8);
                    mesh.merge(mark);
                }
            }
            object->add_volume(std::move(mesh));
            object->add_instance()->set_offset(Vec3d(20. + shape * 30., 20. + row * 30., 0.));
            DynamicPrintConfig cfg;
            // Preserve the user's other scarf settings; only the selected value varies.
            for (const char *key : {"seam_slope_start_height", "seam_slope_min_length", "seam_slope_steps",
                                   "scarf_joint_flow_ratio", "scarf_joint_speed"})
                if (const auto *option = process.option(key)) cfg.set_key_value(key, option->clone());
            cfg.set_deserialize_strict({{"extruder", 1}, {"wall_loops", 3}, {"sparse_infill_density", "0%"},
                {"top_shell_layers", 0}, {"bottom_shell_layers", 0}, {"seam_position", "back"},
                {"seam_slope_type", inner ? "none" : "all"}, {"seam_slope_conditional", false},
                {"seam_slope_entire_loop", false}, {"seam_slope_inner_walls", false},
                {"seam_start_on_inner_wall", inner && row == 1}, {"wall_sequence", "inner wall/outer wall"},
                {"perimeter_layering", "standard"}, {"arc_overhang_enabled", false},
                {"enable_support", false}, {"alternate_extra_wall", false}, {"overhang_reverse", false}});
            cfg.set_key_value("layer_height", new ConfigOptionFloat(p.nozzle / 2.));
            object->config.assign_config(std::move(cfg));
            if (!inner) {
                object->layer_config_ranges[{0., 8.}].set_key_value("seam_slope_type",
                    new ConfigOptionEnum<SeamScarfType>(SeamScarfType::None));
                for (int i = 0; i < p.count; ++i) {
                    auto &range = object->layer_config_ranges[{8. * (i + 1), 8. * (i + 2)}];
                    const char *key = keys[int(p.kind)];
                    if (p.kind == SeamCalibrationKind::StartHeight || p.kind == SeamCalibrationKind::Speed)
                        range.set_key_value(key, new ConfigOptionFloatOrPercent(values[i], p.kind == SeamCalibrationKind::StartHeight));
                    else range.set_key_value(key, new ConfigOptionFloat(values[i]));
                }
                // Layer-range consumers require an explicit height even when
                // the range changes only a seam parameter.
                for (auto &[bounds, range] : object->layer_config_ranges)
                    range.set_key_value("layer_height", new ConfigOptionFloat(p.nozzle / 2.));
            }
        }
}
}
