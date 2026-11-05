#include "DepositionSection.hpp"
#include "../Geometry/ArcWelder.hpp"
#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

namespace Slic3r {
DepositionSections::DepositionSections(std::uint64_t generation) : m_generation(generation)
{
    if (!generation)
        throw std::invalid_argument("Missing deposition section generation.");
}

std::uint32_t DepositionSections::add(const DepositionCommandSection& s)
{
    if (m_sealed || m_sections.size() >= maximum_commands)
        throw std::length_error("Deposition section command budget exceeded.");
    if (s.transform != DepositionCommandSection::Transform::Ordinary &&
        s.transform != DepositionCommandSection::Transform::BrickCourse)
        throw std::invalid_argument("Invalid deposition section transform.");
    if (s.purpose != DepositionCommandSection::Purpose::Model &&
        s.purpose != DepositionCommandSection::Purpose::AuxiliaryBrim)
        throw std::invalid_argument("Invalid deposition section purpose.");
    if (s.purpose == DepositionCommandSection::Purpose::AuxiliaryBrim &&
        (s.transform != DepositionCommandSection::Transform::Ordinary || s.start[2] != s.end[2] ||
         s.start_height <= 0. || s.start_height != s.end_height || s.nominal_lower_z))
        throw std::invalid_argument("Invalid auxiliary brim section geometry.");
    double length2 = 0.;
    for (std::size_t i = 0; i < 3; ++i) {
        if (!std::isfinite(s.start[i]) || !std::isfinite(s.end[i]))
            throw std::invalid_argument("Invalid deposition section coordinates.");
        length2 += (s.end[i] - s.start[i]) * (s.end[i] - s.start[i]);
    }
    if (!std::isfinite(s.width) || s.width <= 0.)
        throw std::invalid_argument("Invalid deposition section width.");
    if (!std::isfinite(s.start_height) || !std::isfinite(s.end_height))
        throw std::invalid_argument("Non-finite deposition section height.");
    if (s.start_height < 0. || s.end_height < 0.)
        throw std::invalid_argument("Negative deposition section height.");
    // A positively extruding command can collapse to its certified lower
    // surface at output precision. Keep that explicit unresolved state; it is
    // not a bead and must never receive a fallback or minimum thickness.
    if (std::max(s.start_height, s.end_height) <= 0. && !s.nominal_lower_z)
        throw std::invalid_argument("Both deposition section heights are zero.");
    if (s.nominal_lower_z &&
        (!std::isfinite((*s.nominal_lower_z)[0]) || !std::isfinite((*s.nominal_lower_z)[1]) ||
         s.start_height != s.start[2] - (*s.nominal_lower_z)[0] ||
         s.end_height != s.end[2] - (*s.nominal_lower_z)[1]))
        throw std::invalid_argument("Deposition lower surface and design heights disagree.");
    if (!std::isfinite(length2) || length2 <= 0.)
        throw std::invalid_argument("Stationary or invalid deposition section length.");
    if (s.arc_direction != 0) {
        if ((s.arc_direction != -1 && s.arc_direction != 1) ||
            !std::isfinite(s.arc_center[0]) || !std::isfinite(s.arc_center[1]) || s.start[2] != s.end[2] ||
            s.start_height != s.end_height)
            throw std::invalid_argument("Invalid planar deposition arc.");
        const double radius = std::hypot(s.start[0]-s.arc_center[0],s.start[1]-s.arc_center[1]);
        if (!(radius > 0.) || !std::isfinite(radius) || radius > 1e7)
            throw std::invalid_argument("Invalid deposition arc radius.");
    }
    // Do not let implementation-specific vector growth exceed the hard cap.
    if (m_sections.size() == m_sections.capacity())
        m_sections.reserve(std::min(maximum_commands, std::max(std::size_t(64), m_sections.capacity() * 2)));
    m_sections.push_back(s);
    return static_cast<std::uint32_t>(m_sections.size());
}

const DepositionCommandSection* DepositionSections::section(std::uint32_t id) const noexcept
{ return id && id <= m_sections.size() ? &m_sections[id - 1] : nullptr; }

std::string DepositionSections::marker(std::uint32_t id) const
{ return ";_ORCA_INTERNAL_SECTION:" + std::to_string(m_generation) + ":" + std::to_string(id) + "\n"; }

bool DepositionSectionCapture::consume(std::string_view line)
{
    constexpr std::string_view prefix = ";_ORCA_INTERNAL_SECTION";
    if (line.substr(0, prefix.size()) != prefix)
        return false;
    line.remove_prefix(prefix.size());
    if (!line.empty() && line.back() == '\n')
        line.remove_suffix(1);
    if (!line.empty() && line.back() == '\r')
        line.remove_suffix(1);
    std::uint64_t generation = 0;
    std::uint32_t id         = 0;
    bool parsed              = !line.empty() && line.front() == ':';
    if (parsed) {
        line.remove_prefix(1);
        const auto first = std::from_chars(line.data(), line.data() + line.size(), generation);
        parsed           = first.ec == std::errc{} && first.ptr != line.data() + line.size() && *first.ptr == ':';
        if (parsed) {
            const auto second = std::from_chars(first.ptr + 1, line.data() + line.size(), id);
            parsed            = second.ec == std::errc{} && second.ptr == line.data() + line.size();
        }
    }
    if (!parsed || generation != m_generation || (id && m_active) || (!id && !m_active) ||
        m_events.size() >= 2 * DepositionSections::maximum_commands)
        m_valid = false;
    if (m_valid) {
        m_active = id;
        m_events.push_back(id);
    }
    return true;
}

bool DepositionSectionCapture::finish(const DepositionSections& sections)
{
    if (!m_valid || m_active || !sections.sealed() || sections.generation() != m_generation || m_events.size() != 2 * sections.size())
        return false;
    for (std::size_t i = 0; i < sections.size(); ++i)
        if (m_events[2 * i] != i + 1 || m_events[2 * i + 1] != 0)
            return false;
    return true;
}

std::optional<std::array<double, 2>> DepositionSectionChain::append(const std::array<double, 3>& start,
                                                                   const std::array<double, 3>& end,
                                                                   std::uint32_t filament)
{
    if (!m_valid || filament != m_parent.filament) {
        m_valid = false;
        return {};
    }
    double length2 = 0.;
    for (std::size_t i = 0; i < 3; ++i)
        length2 += std::pow(m_parent.end[i] - m_parent.start[i], 2);
    auto fraction = [&](const std::array<double, 3>& p) -> std::optional<double> {
        if (m_parent.arc_direction) {
            if (m_parent.arc_direction != -1 && m_parent.arc_direction != 1) return {};
            for (const auto coordinate : p)
                if (!std::isfinite(coordinate)) return {};
            // The last endpoint is explicitly commanded, and can differ from
            // the start radius after decimal I/J and XYZ quantization.
            if (p == m_parent.start) return 0.;
            if (p == m_parent.end) return 1.;
            const auto& c = m_parent.arc_center;
            const double x = m_parent.start[0]-c[0], y = m_parent.start[1]-c[1];
            const auto angle = [&](const std::array<double,3>& q) {
                double a = m_parent.arc_direction * std::atan2(x*(q[1]-c[1])-y*(q[0]-c[0]),
                                                               x*(q[0]-c[0])+y*(q[1]-c[1]));
                if (a < 0.) a += 2.*std::acos(-1.);
                return a;
            };
            const double total = angle(m_parent.end), a = angle(p);
            if (!(total > 0.)) return {};
            const double radius = std::hypot(x,y);
            if (!(radius > 0.) || !std::isfinite(radius) || radius > 1e7) return {};
            // GCodeProcessor's legacy arc interpolation uses a truncated
            // rotation for at most 25 steps before exact correction. Bound
            // its accumulated radial drift analytically, not from the child
            // being authenticated. The Marlin exact-trig rotation is tighter.
            const auto steps = Geometry::ArcWelder::arc_discretization_steps(radius,total,.0125);
            const double theta = total / double(steps), square = theta*theta;
            const double rotation_norm = std::hypot(1.-.5*square,theta-theta*square/6.);
            const double drift = radius * std::abs(std::pow(rotation_norm,double(std::min<std::size_t>(25,steps-1)))-1.) +
                50.*std::numeric_limits<float>::epsilon()*radius;
            if (a > total || std::abs(p[2]-m_parent.start[2]) > .0006 ||
                std::abs(std::hypot(p[0]-c[0],p[1]-c[1])-radius) > .0006+drift)
                return {};
            return a / total;
        }
        double dot = 0.;
        for (std::size_t i = 0; i < 3; ++i) {
            if (!std::isfinite(p[i]))
                return {};
            dot += (p[i] - m_parent.start[i]) * (m_parent.end[i] - m_parent.start[i]);
        }
        const double t = dot / length2;
        if (!std::isfinite(t) || t < -1e-7 || t > 1. + 1e-7)
            return {};
        for (std::size_t i = 0; i < 3; ++i)
            if (std::abs(p[i] - (m_parent.start[i] + t * (m_parent.end[i] - m_parent.start[i]))) > 0.0006)
                return {};
        return std::clamp(t, 0., 1.);
    };
    const auto a = fraction(start), b = fraction(end);
    if (!a || !b || std::abs(*a - m_fraction) > 1e-7 || *b <= *a) {
        m_valid = false;
        return {};
    }
    std::array<double, 2> heights;
    if (m_parent.nominal_lower_z) {
        const auto& lower = *m_parent.nominal_lower_z;
        heights = {start[2] - (lower[0] + *a * (lower[1] - lower[0])),
                   end[2] - (lower[0] + *b * (lower[1] - lower[0]))};
        if (heights[0] < 0. || heights[1] < 0.) {
            m_valid = false;
            return {};
        }
    } else {
        // Preserve the legacy flow-derived envelope for unanchored commands.
        heights = {float(m_parent.start_height + *a * (m_parent.end_height - m_parent.start_height)),
                   float(m_parent.start_height + *b * (m_parent.end_height - m_parent.start_height))};
    }
    m_fraction = *b;
    return heights;
}
} // namespace Slic3r
