#pragma once

#include "BoundingBox.hpp"
#include "ExPolygon.hpp"
#include "Point.hpp"
#include "Polygon.hpp"
#include "ThermoelasticExactSections.hpp"
#include <cstddef>
#include <initializer_list>
#include <optional>
#include <stdexcept>
#include <utility>
#include <variant>

namespace Slic3r::Thermoelastic {

// One authoritative representation, never an exact handle plus a rounded
// fallback. Unported consumers must request native() explicitly and reject
// exact geometry, including an empty exact set, rather than silently drop it.
class SectionFootprint {
public:
    SectionFootprint() = default;
    SectionFootprint(ExPolygons polygons) : m_geometry(std::move(polygons)) {}
    SectionFootprint(std::initializer_list<ExPolygon> polygons) : m_geometry(ExPolygons(polygons)) {}
    explicit SectionFootprint(ExactSection section) : m_geometry(std::move(section)) {}

    bool is_exact() const { return std::holds_alternative<ExactSection>(m_geometry); }
    ExPolygons& native()
    {
        if (is_exact()) throw std::invalid_argument("Exact section footprint requires a geometry-aware consumer.");
        return std::get<ExPolygons>(m_geometry);
    }
    const ExPolygons& native() const
    {
        if (is_exact()) throw std::invalid_argument("Exact section footprint requires a geometry-aware consumer.");
        return std::get<ExPolygons>(m_geometry);
    }
    const ExactSection& exact() const
    {
        if (!is_exact()) throw std::invalid_argument("Native section footprint has no exact canonical handle.");
        return std::get<ExactSection>(m_geometry);
    }
    bool empty() const { return is_exact() ? exact().empty() : native().empty(); }
    std::size_t vertex_count() const
    {
        if (is_exact()) return exact().vertex_count();
        std::size_t count=0;
        for (const auto& polygon : native()) {
            count+=polygon.contour.points.size();
            for (const auto& hole : polygon.holes) count+=hole.points.size();
        }
        return count;
    }
    // Conservative indexing metadata only; never a polygon replacement.
    std::optional<BoundingBox> bounds() const
    {
        if (is_exact()) return exact().bounds();
        std::optional<BoundingBox> box;
        const auto include = [&](const Polygon& ring) {
            for (const Point& point : ring.points) {
                if (!box) box.emplace();
                box->merge(point);
            }
        };
        for (const auto& polygon : native()) {
            include(polygon.contour);
            for (const auto& hole : polygon.holes) include(hole);
        }
        return box;
    }
private:
    std::variant<ExPolygons,ExactSection> m_geometry;
};

} // namespace Slic3r::Thermoelastic
