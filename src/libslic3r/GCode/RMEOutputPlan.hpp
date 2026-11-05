#pragma once

#include <algorithm>
#include <cstddef>
#include <ios>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace Slic3r {

// Collected from deposition, not heater commands: probing, idle temperatures,
// purges without XY motion, and custom startup/cleaning moves are excluded.
struct RMEOutputPlan {
    struct Tool {
        int physical = 0;
        std::string material;
        int initial = 0;
        std::set<int> temperatures;
    };
    std::map<int, Tool> tools;
    std::array<double, 2> minimum{INFINITY, INFINITY};
    std::array<double, 2> maximum{-INFINITY, -INFINITY};
    bool valid = true;

    void point(double x, double y, double width)
    {
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(width) || width <= 0.) {
            valid = false;
            return;
        }
        minimum[0] = std::min(minimum[0], x - width / 2.);
        minimum[1] = std::min(minimum[1], y - width / 2.);
        maximum[0] = std::max(maximum[0], x + width / 2.);
        maximum[1] = std::max(maximum[1], y + width / 2.);
    }

    void temperature(int logical, int physical, const std::string &material, double target)
    {
        auto &tool = tools[logical];
        tool.physical = physical;
        tool.material = material;
        if (!std::isfinite(target) || target <= 0. || target > 500.) {
            tool.temperatures.insert(0); // Do not silently omit unknown deposition targets.
            return;
        }
        const int value = int(std::lround(target));
        if (tool.initial == 0) tool.initial = value;
        tool.temperatures.insert(value);
    }

    void arc(double sx, double sy, double ex, double ey, double cx, double cy, bool clockwise, double width)
    {
        constexpr double pi = 3.14159265358979323846;
        const auto positive = [](double angle) {
            angle = std::fmod(angle, 2. * pi);
            return angle < 0. ? angle + 2. * pi : angle;
        };
        const double start = std::atan2(sy - cy, sx - cx);
        const double end = std::atan2(ey - cy, ex - cx);
        const double radius = std::hypot(sx - cx, sy - cy);
        const bool full = std::hypot(ex - sx, ey - sy) < 1e-6;
        const double sweep = full ? 2. * pi : positive(clockwise ? start - end : end - start);
        point(sx, sy, width);
        point(ex, ey, width);
        for (int i = 0; i < 4; ++i) {
            const double angle = i * pi / 2.;
            const double distance = positive(clockwise ? start - angle : angle - start);
            if (distance <= sweep + 1e-9)
                point(cx + radius * std::cos(angle), cy + radius * std::sin(angle), width);
        }
    }

    std::string pa_command() const
    {
        if (tools.empty()) return {};
        std::ostringstream out;
        out.imbue(std::locale::classic());
        out << "M976 A ";
        bool first = true;
        for (const auto &[logical, tool] : tools) {
            if (logical < 0 || tool.physical < 0 || tool.initial <= 0 || tool.temperatures.count(0) || tool.material.empty() ||
                tool.material.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-") != std::string::npos)
                return {};
            if (!first) out << ',';
            first = false;
            out << tool.physical << ':' << logical << ':' << tool.material << ':' << tool.initial;
            std::vector<int> selected{tool.initial};
            const auto add = [&](int value) {
                if (selected.size() < 8 && std::find(selected.begin(), selected.end(), value) == selected.end())
                    selected.push_back(value);
            };
            add(*tool.temperatures.begin());
            add(*tool.temperatures.rbegin());
            for (int value : tool.temperatures) add(value);
            for (size_t i = 1; i < selected.size(); ++i) out << '/' << selected[i];
        }
        return out.str() + "\n";
    }

    std::string mesh_command() const
    {
        if (!valid || !std::isfinite(minimum[0]) || !std::isfinite(minimum[1]) ||
            !std::isfinite(maximum[0]) || !std::isfinite(maximum[1]) ||
            maximum[0] <= minimum[0] || maximum[1] <= minimum[1])
            return {};
        // Round outwards so serializing the bounds never clips the bead.
        const double x = std::floor(minimum[0] * 10000.) / 10000.;
        const double y = std::floor(minimum[1] * 10000.) / 10000.;
        const double right = std::ceil(maximum[0] * 10000.) / 10000.;
        const double top = std::ceil(maximum[1] * 10000.) / 10000.;
        std::ostringstream out;
        out.imbue(std::locale::classic());
        out << std::fixed << std::setprecision(4) << "M555 X" << x << " Y" << y
            << " W" << right - x << " H" << top - y << "\n";
        return out.str();
    }
};

} // namespace Slic3r
