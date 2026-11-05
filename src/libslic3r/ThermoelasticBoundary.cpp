#include "ThermoelasticBoundary.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticExactSections.hpp"
#include "ClipperUtils.hpp"
#include "ExPolygon.hpp"
#include "Polygon.hpp"
#include "Point.hpp"
#include "BoundingBox.hpp"
#include "Triangulation.hpp"
#include "libslic3r.h"
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
namespace {
boost::multiprecision::int256_t ring_twice_area(const Polygon& ring)
{
    using Integer = boost::multiprecision::int256_t;
    Integer twice = 0;
    if (ring.points.size() < 3) throw std::invalid_argument("Invalid boundary-uncertainty ring.");
    const auto origin = ring.points.front();
    for (std::size_t i = 1; i + 1 < ring.points.size(); ++i) {
        const auto& a = ring.points[i];
        const auto& b = ring.points[i + 1];
        twice += (Integer(a.x()) - origin.x()) * (Integer(b.y()) - origin.y()) -
            (Integer(a.y()) - origin.y()) * (Integer(b.x()) - origin.x());
    }
    return twice;
}
double checked_area(const ExPolygons& polygons)
{
    boost::multiprecision::int256_t result = 0;
    for (const auto& polygon : polygons) {
        auto value = ring_twice_area(polygon.contour);
        if (value <= 0) throw std::invalid_argument("Invalid boundary-uncertainty contour orientation.");
        for (const auto& hole : polygon.holes) {
            const auto hole_area = ring_twice_area(hole);
            if (hole_area >= 0) throw std::invalid_argument("Invalid boundary-uncertainty hole orientation.");
            value += hole_area;
        }
        if (value < 0) throw std::invalid_argument("Invalid boundary-uncertainty polygon orientation.");
        result += value;
    }
    return .5 * result.convert_to<double>();
}
double coverage_area_bound_m2(const ExPolygons& polygons, const std::function<void()>& checkpoint)
{
    using Integer = boost::multiprecision::int256_t;
    Integer twice = 0;
    const auto ring = [&](const Polygon& boundary) {
        const auto& origin = boundary.points.front();
        for (std::size_t i = 1; i + 1 < boundary.points.size(); ++i) {
            if (checkpoint) checkpoint();
            const auto& a = boundary.points[i];
            const auto& b = boundary.points[i + 1];
            const Integer cross = (Integer(a.x()) - origin.x()) * (Integer(b.y()) - origin.y()) -
                (Integer(a.y()) - origin.y()) * (Integer(b.x()) - origin.x());
            twice += cross < 0 ? -cross : cross;
        }
    };
    for (const auto& polygon : polygons) {
        ring(polygon.contour);
        for (const auto& hole : polygon.holes) ring(hole);
    }
    if (twice == 0) return 0.;
    const double area = twice.convert_to<double>() * .5e-18;
    // Round outward after integer conversion and unit conversion. The sum of
    // absolute fan areas bounds nonzero-winding coverage, even for crossing
    // diagnostic rings; it is deliberately not a signed area measurement.
    return std::nextafter(std::nextafter(area, std::numeric_limits<double>::infinity()), std::numeric_limits<double>::infinity());
}
} // namespace

std::size_t Detail::certify_boundary_triangulation(const ExPolygon& polygon,
    std::vector<std::array<std::size_t, 3>>& triangles, std::size_t max_work,
    const std::function<void()>& checkpoint)
{
    using Integer = boost::multiprecision::int256_t;
    std::size_t input_points = polygon.contour.points.size();
    if (input_points > 512 || polygon.holes.size() > 128)
        throw ResourceLimit("Boundary triangulation input exceeds its point budget.");
    for (const auto& hole : polygon.holes) {
        if (hole.points.size() > 512 - input_points)
            throw ResourceLimit("Boundary triangulation input exceeds its point budget.");
        input_points += hole.points.size();
    }
    const auto vertices = to_points(polygon);
    if (vertices.size() < 3 || vertices.size() > 512 || triangles.empty() || triangles.size() > 4 * vertices.size() ||
        !max_work || max_work > 20000000)
        throw std::invalid_argument("Invalid boundary triangulation certificate input.");
    std::size_t work = 0;
    const auto tick = [&] {
        if (checkpoint) checkpoint();
        if (++work > max_work) throw ResourceLimit("Boundary triangulation certificate work budget exceeded.");
    };
    std::set<std::pair<coord_t, coord_t>> unique;
    for (const auto& p : vertices) {
        tick();
        if (std::abs(double(p.x())) > 1e10 || std::abs(double(p.y())) > 1e10 || !unique.emplace(p.x(), p.y()).second)
            throw std::invalid_argument("Boundary certificate requires unique bounded contour vertices.");
    }
    checked_area(ExPolygons{polygon});
    const auto orientation = [&](std::size_t a, std::size_t b, std::size_t c) -> Integer {
        return (Integer(vertices[b].x()) - vertices[a].x()) * (Integer(vertices[c].y()) - vertices[a].y()) -
            (Integer(vertices[b].y()) - vertices[a].y()) * (Integer(vertices[c].x()) - vertices[a].x());
    };
    using EdgeMap = std::map<std::pair<std::size_t, std::size_t>, int>;
    EdgeMap actual, expected;
    const auto add_edge = [&](EdgeMap& map, std::size_t a, std::size_t b) {
        std::vector<std::size_t> subdivisions;
        for (std::size_t p = 0; p < vertices.size(); ++p) {
            tick();
            if (vertices[p].x() >= std::min(vertices[a].x(), vertices[b].x()) && vertices[p].x() <= std::max(vertices[a].x(), vertices[b].x()) &&
                vertices[p].y() >= std::min(vertices[a].y(), vertices[b].y()) && vertices[p].y() <= std::max(vertices[a].y(), vertices[b].y()) &&
                orientation(a, b, p) == 0) subdivisions.push_back(p);
        }
        const auto less = [&](std::size_t p, std::size_t q) {
            return vertices[p].x() < vertices[q].x() || (vertices[p].x() == vertices[q].x() && vertices[p].y() < vertices[q].y());
        };
        std::sort(subdivisions.begin(), subdivisions.end(), less);
        if (less(b, a)) std::reverse(subdivisions.begin(), subdivisions.end());
        for (std::size_t i = 1; i < subdivisions.size(); ++i) {
            const auto p = subdivisions[i - 1], q = subdivisions[i];
            const auto key = std::minmax(p, q);
            auto& count = map[key];
            count += p < q ? 1 : -1;
            if (!count) map.erase(key);
        }
    };
    std::size_t offset = 0;
    const auto add_ring = [&](const Polygon& boundary) {
        for (std::size_t i = 0; i < boundary.points.size(); ++i)
            add_edge(expected, offset + i, offset + (i + 1) % boundary.points.size());
        offset += boundary.points.size();
    };
    add_ring(polygon.contour);
    for (const auto& hole : polygon.holes) add_ring(hole);
    for (auto& triangle : triangles) {
        for (const auto index : triangle)
            if (index >= vertices.size()) throw std::invalid_argument("Boundary certificate triangle index is invalid.");
        const Integer turn = orientation(triangle[0], triangle[1], triangle[2]);
        if (turn == 0) continue;
        if (turn < 0) std::swap(triangle[1], triangle[2]);
        for (unsigned i = 0; i < 3; ++i) add_edge(actual, triangle[i], triangle[(i + 1) % 3]);
    }
    if (actual != expected)
        throw std::invalid_argument("Boundary certificate triangulation changes the directed polygon boundary.");
    return work;
}

SectionBoundaryError qualify_section_boundary_error(const ExPolygons& original,
    const ExPolygons& represented, unsigned band, std::size_t max_points,
    std::size_t max_work, const std::function<void()>& checkpoint,
    const std::vector<const ExPolygons*>& material_regions)
{
    if (band > 4 || !max_points || max_points > 2000000 || !max_work || max_work > 20000000 || material_regions.size() > 256)
        throw std::invalid_argument("Invalid explicit section boundary-uncertainty budget.");
    SectionBoundaryError result;
    const auto points = [&](const ExPolygons& polygons) {
        std::size_t count = 0;
        for (const auto& polygon : polygons) {
            const auto ring = [&](const Polygon& boundary) {
                if (boundary.points.size() < 3 || boundary.points.size() > max_points - count)
                    throw ResourceLimit("Section boundary-uncertainty point budget exceeded.");
                count += boundary.points.size();
                for (const auto& point : boundary.points)
                    if (std::abs(double(point.x())) > 1e10 || std::abs(double(point.y())) > 1e10)
                        throw std::invalid_argument("Boundary-uncertainty coordinate outside the supported range.");
            };
            ring(polygon.contour);
            for (const auto& hole : polygon.holes) ring(hole);
        }
        return count;
    };
    bool callback_failed = false;
    const auto charge = [&](const ExPolygons& a, const ExPolygons& b = ExPolygons{}) {
        if (checkpoint) try { checkpoint(); } catch (...) { callback_failed = true; throw; }
        const auto count = points(a) + points(b) + 1;
        if (count > max_work - result.work) throw ResourceLimit("Section boundary-uncertainty work budget exceeded.");
        result.work += count;
    };
    const auto bounded = [&](ExPolygons polygons) { points(polygons); return polygons; };
    charge(original, represented);
    checked_area(original);
    checked_area(represented);
    if (band) {
        // Test input islands before union: integer clipping may discard a tiny
        // ring, which must not turn an unresolved feature into exact equality.
        const auto input_survives = [&](const ExPolygons& polygons) {
            for (const auto& polygon : polygons) {
                charge(ExPolygons{polygon});
                if (offset_ex(ExPolygons{polygon}, -float(band + 5), jtRound, .25).empty())
                    throw std::invalid_argument("An input island is unresolved within boundary uncertainty.");
            }
        };
        input_survives(original);
        input_survives(represented);
    }
    // A different boundary decomposition can round to different integer
    // unions even when the original occupied sets are identical. Such an
    // identity introduces no representation uncertainty, so do not require
    // existing holes to survive an artificial nonzero perturbation. Retain
    // input validation and island-survival guards above, and the existing
    // morphology checks below for every genuinely changed representation.
    bool exactly_changed=false;
    if (points(original)+points(represented)<=1024) {
        if (result.work==max_work) throw ResourceLimit("Exact boundary identity work budget exhausted.");
        const auto identity=exact_section_empty(original,represented,ExactSectionOperation::SymmetricDifference,
            1024,std::min<std::size_t>(4000000,max_work-result.work),checkpoint,true);
        result.work+=identity.work;
        if (identity.empty) return result;
        exactly_changed=true;
    }
    const auto before = bounded(union_ex(original)), after = bounded(union_ex(represented));
    charge(before, after);
    const auto added = bounded(diff_ex(after, before)), removed = bounded(diff_ex(before, after));
    if (added.empty() && removed.empty()) {
        if (!exactly_changed) throw ResourceLimit("Exact boundary identity input budget exceeded.");
    }
    if (!band) throw std::invalid_argument("Section ownership changed the deposited union.");
    if (exactly_changed) {
        if (result.work==max_work) throw ResourceLimit("Exact boundary difference work budget exhausted.");
        const auto exact=certify_exact_boundary_difference(original,represented,band+3,
            std::min<std::size_t>(4000000,max_work-result.work),checkpoint);
        result.work+=exact.work;
        result.added_area_bound_m2=exact.added_area_bound_m2;
        result.removed_area_bound_m2=exact.removed_area_bound_m2;
        // No early acceptance: the material and occupied morphology checks
        // below remain required, even when rounded differences contain no area.
    }
    ExPolygons boundary_band;
    if (!exactly_changed) {
        charge(before);
        const auto outer = bounded(offset_ex(before, float(band), jtRound, .25));
        const auto inner = bounded(offset_ex(before, -float(band), jtRound, .25));
        charge(outer, inner);
        boundary_band = bounded(diff_ex(outer, inner));
    }
    using Integer = boost::multiprecision::int256_t;
    using DyadicPoint = std::array<std::int64_t, 2>;
    constexpr unsigned maximum_depth = 24;
    constexpr std::int64_t scale = std::int64_t(1) << maximum_depth;
    const auto dyadic = [&](const Point& p) { return DyadicPoint{std::int64_t(p.x()) * scale, std::int64_t(p.y()) * scale}; };
    const auto in_segment_capsule = [&](const DyadicPoint& p, const DyadicPoint& a, const DyadicPoint& b) {
        // Coordinates <=1e10, with24 dyadic bits, need fewer than240 bits
        // after squaring cross products. All midpoints remain exact.
        const Integer x = Integer(p[0]) - a[0], y = Integer(p[1]) - a[1];
        const Integer dx = Integer(b[0]) - a[0], dy = Integer(b[1]) - a[1];
        const Integer length2 = dx * dx + dy * dy, dot = x * dx + y * dy;
        // Use the same declared Euclidean allowance as the rounded band path,
        // not an additional tolerance beyond what the caller will receive.
        const Integer radius = Integer(band + 3) * scale, radius2 = radius * radius;
        if (dot <= 0) return x * x + y * y <= radius2;
        if (dot >= length2) return (x - dx) * (x - dx) + (y - dy) * (y - dy) <= radius2;
        const Integer cross = x * dy - y * dx;
        return cross * cross <= radius2 * length2;
    };
    const auto certify_band = [&](const ExPolygons& changed) {
        charge(changed, boundary_band);
        if (diff_ex(changed, boundary_band).empty()) return;
        // Original-index triangulation avoids over-covering concave ribbons.
        // Its directed boundary is certified below before any triangle is used.
        for (const auto& polygon : changed) {
            struct Triangle { std::array<DyadicPoint, 3> points; unsigned depth; };
            std::vector<Triangle> pending;
            pending.reserve(maximum_depth + 1);
            const auto ring = to_points(polygon);
            if (ring.size() > 512) throw ResourceLimit("Boundary certificate triangulation exceeds its point budget.");
            const auto orient_dyadic = [&](const DyadicPoint& a, const DyadicPoint& b, const DyadicPoint& c) -> Integer {
                return (Integer(b[0]) - a[0]) * (Integer(c[1]) - a[1]) -
                    (Integer(b[1]) - a[1]) * (Integer(c[0]) - a[0]);
            };
            const auto outside_changed = [&](const Triangle& triangle) {
                const auto between = [](const DyadicPoint& p, const DyadicPoint& a, const DyadicPoint& b) {
                    return p[0] >= std::min(a[0], b[0]) && p[0] <= std::max(a[0], b[0]) &&
                        p[1] >= std::min(a[1], b[1]) && p[1] <= std::max(a[1], b[1]);
                };
                int winding = 0;
                const auto& sample = triangle.points[0];
                const auto check_ring = [&](const Polygon& boundary) {
                    for (std::size_t i = 0; i < boundary.points.size(); ++i) {
                        charge({});
                        const auto a = dyadic(boundary.points[i]), b = dyadic(boundary.points[(i + 1) % boundary.points.size()]);
                        const Integer u = orient_dyadic(triangle.points[0], triangle.points[1], a);
                        const Integer v = orient_dyadic(triangle.points[1], triangle.points[2], a);
                        const Integer w = orient_dyadic(triangle.points[2], triangle.points[0], a);
                        if ((u >= 0 && v >= 0 && w >= 0) || (u <= 0 && v <= 0 && w <= 0)) return false;
                        for (unsigned edge = 0; edge < 3; ++edge) {
                            const auto& c = triangle.points[edge];
                            const auto& d = triangle.points[(edge + 1) % 3];
                            const Integer p = orient_dyadic(a, b, c), q = orient_dyadic(a, b, d);
                            const Integer r = orient_dyadic(c, d, a), s = orient_dyadic(c, d, b);
                            if ((p == 0 && between(c, a, b)) || (q == 0 && between(d, a, b)) ||
                                (r == 0 && between(a, c, d)) || (s == 0 && between(b, c, d)) ||
                                (((p < 0 && q > 0) || (p > 0 && q < 0)) && ((r < 0 && s > 0) || (r > 0 && s < 0)))) return false;
                        }
                        const Integer side = orient_dyadic(a, b, sample);
                        if (a[1] <= sample[1] && b[1] > sample[1] && side > 0) ++winding;
                        if (a[1] > sample[1] && b[1] <= sample[1] && side < 0) --winding;
                    }
                    return true;
                };
                if (!check_ring(polygon.contour)) return false;
                for (const auto& hole : polygon.holes) if (!check_ring(hole)) return false;
                return winding == 0;
            };
            const auto certify_triangle = [&](std::size_t a, std::size_t b, std::size_t c) {
                pending.push_back({{dyadic(ring[a]), dyadic(ring[b]), dyadic(ring[c])}, 0});
                while (!pending.empty()) {
                    charge({});
                    const Triangle triangle = pending.back();
                    pending.pop_back();
                    auto lower = triangle.points[0], upper = lower;
                    for (const auto& p : triangle.points) for (unsigned axis = 0; axis < 2; ++axis) {
                        lower[axis] = std::min(lower[axis], p[axis]);
                        upper[axis] = std::max(upper[axis], p[axis]);
                    }
                    bool certified = false;
                    for (const auto& source : before) {
                        const auto check_ring = [&](const Polygon& boundary) {
                            for (std::size_t i = 0; i < boundary.points.size() && !certified; ++i) {
                                charge({});
                                const auto a = dyadic(boundary.points[i]);
                                const auto b = dyadic(boundary.points[(i + 1) % boundary.points.size()]);
                                const std::int64_t radius = std::int64_t(band + 3) * scale;
                                if (lower[0] < std::min(a[0], b[0]) - radius || lower[1] < std::min(a[1], b[1]) - radius ||
                                    upper[0] > std::max(a[0], b[0]) + radius || upper[1] > std::max(a[1], b[1]) + radius) continue;
                                certified = true;
                                for (const auto& p : triangle.points) certified = certified && in_segment_capsule(p, a, b);
                            }
                        };
                        check_ring(source.contour);
                        for (const auto& hole : source.holes) check_ring(hole);
                        if (certified) break;
                    }
                    if (certified) continue;
                    // Only prune fan excess after proving no boundary touches
                    // or crosses it and its constant winding is zero. A touch
                    // remains ambiguous and must be covered or subdivided.
                    if (outside_changed(triangle)) continue;
                    if (triangle.depth == maximum_depth) {
                        std::string diagnostic;
                        const auto nearest_distance = [&](const DyadicPoint& p) {
                            double nearest = std::numeric_limits<double>::infinity();
                            const auto distance_ring = [&](const Polygon& boundary) {
                                for (std::size_t i = 0; i < boundary.points.size(); ++i) {
                                    charge({});
                                    const auto& a = boundary.points[i];
                                    const auto& b = boundary.points[(i + 1) % boundary.points.size()];
                                    const double x = double(p[0]) / scale - double(a.x());
                                    const double y = double(p[1]) / scale - double(a.y());
                                    const double dx = double(b.x() - a.x()), dy = double(b.y() - a.y());
                                    const double length2 = dx * dx + dy * dy;
                                    const double t = length2 > 0. ? std::clamp((x * dx + y * dy) / length2, 0., 1.) : 0.;
                                    nearest = std::min(nearest, std::hypot(x - t * dx, y - t * dy));
                                }
                            };
                            for (const auto& source : before) {
                                distance_ring(source.contour);
                                for (const auto& hole : source.holes) distance_ring(hole);
                            }
                            return nearest;
                        };
                        const auto winding_in = [&](const ExPolygons& polygons, const DyadicPoint& p) {
                            int winding = 0;
                            const auto boundary = [&](const Polygon& boundary_ring) {
                                for (std::size_t i = 0; i < boundary_ring.points.size(); ++i) {
                                    charge({});
                                    const auto a = dyadic(boundary_ring.points[i]);
                                    const auto b = dyadic(boundary_ring.points[(i + 1) % boundary_ring.points.size()]);
                                    const Integer side = orient_dyadic(a, b, p);
                                    if (a[1] <= p[1] && b[1] > p[1] && side > 0) ++winding;
                                    if (a[1] > p[1] && b[1] <= p[1] && side < 0) --winding;
                                }
                            };
                            for (const auto& source : polygons) {
                                boundary(source.contour);
                                for (const auto& hole : source.holes) boundary(hole);
                            }
                            return winding;
                        };
                        {
                            const auto& p = triangle.points.front();
                            diagnostic += " [" + std::to_string(double(p[0]) / scale) + "," +
                                std::to_string(double(p[1]) / scale) + "; nearest_nm=" + std::to_string(nearest_distance(p)) +
                                "; winding=" + std::to_string(winding_in(before, p)) + "->" + std::to_string(winding_in(after, p)) + "]";
                        }
                        double vertex_maximum = 0.;
                        for (const auto& p : ring) vertex_maximum = std::max(vertex_maximum, nearest_distance(dyadic(p)));
                        diagnostic += " max_ring_vertex_distance_nm=" + std::to_string(vertex_maximum);
                        diagnostic += " certificate_work=" + std::to_string(result.work);
                        throw std::invalid_argument("Section geometry cannot be certified within its boundary uncertainty band:" + diagnostic);
                    }
                    unsigned edge = 0;
                    Integer longest = -1;
                    for (unsigned i = 0; i < 3; ++i) {
                        const auto& a = triangle.points[i];
                        const auto& b = triangle.points[(i + 1) % 3];
                        const Integer dx = Integer(a[0]) - b[0], dy = Integer(a[1]) - b[1];
                        const Integer length = dx * dx + dy * dy;
                        if (length > longest) { longest = length; edge = i; }
                    }
                    const auto& a = triangle.points[edge];
                    const auto& b = triangle.points[(edge + 1) % 3];
                    const DyadicPoint midpoint{(a[0] + b[0]) / 2, (a[1] + b[1]) / 2};
                    auto left = triangle, right = triangle;
                    left.depth = right.depth = triangle.depth + 1;
                    left.points[edge] = midpoint;
                    right.points[(edge + 1) % 3] = midpoint;
                    pending.push_back(left);
                    pending.push_back(right);
                }
            };
            const auto orientation = [&](std::size_t a, std::size_t b, std::size_t c) -> Integer {
                return (Integer(ring[b].x()) - ring[a].x()) * (Integer(ring[c].y()) - ring[a].y()) -
                    (Integer(ring[b].y()) - ring[a].y()) * (Integer(ring[c].x()) - ring[a].x());
            };
            std::set<std::pair<coord_t, coord_t>> unique;
            for (const auto& p : ring) {
                charge({});
                if (!unique.emplace(p.x(), p.y()).second)
                    throw std::invalid_argument("Boundary certificate requires unique contour vertices.");
            }
            // Reserve quadratic work before the bounded, noninterruptible CDT.
            for (std::size_t i = 0; i < ring.size(); ++i) charge(ExPolygons{polygon});
            const auto generated = Triangulation::triangulate(ExPolygons{polygon}, ring);
            if (generated.size() > 4 * ring.size()) throw std::invalid_argument("Boundary certificate has an invalid triangle count.");
            std::vector<std::array<std::size_t, 3>> triangles;
            triangles.reserve(generated.size());
            for (const auto& triangle : generated)
                triangles.push_back({std::size_t(triangle[0]), std::size_t(triangle[1]), std::size_t(triangle[2])});
            if (result.work == max_work) throw ResourceLimit("Boundary certificate work budget exceeded.");
            bool certified_triangulation = true, callback_failed = false;
            const auto check = [&] {
                if (checkpoint) try { checkpoint(); } catch (...) { callback_failed = true; throw; }
            };
            try {
                result.work += Detail::certify_boundary_triangulation(polygon, triangles, max_work - result.work, check);
            } catch (const std::invalid_argument&) {
                if (callback_failed) throw;
                certified_triangulation = false;
                // The failed certificate still consumed work; reserve its
                // bounded worst-case edge-comparison count before fallback.
                const auto reserve = 16 * ring.size() * ring.size();
                if (reserve > max_work - result.work) throw ResourceLimit("Boundary certificate fallback budget exceeded.");
                result.work += reserve;
            }
            if (certified_triangulation) {
                for (const auto& triangle : triangles)
                    if (orientation(triangle[0], triangle[1], triangle[2]) != 0)
                        certify_triangle(triangle[0], triangle[1], triangle[2]);
            } else {
                for (std::size_t i = 1; i + 1 < polygon.contour.points.size(); ++i)
                    if (orientation(0, i, i + 1) != 0) certify_triangle(0, i, i + 1);
            }
        }
    };
    if (!exactly_changed) {
        certify_band(added);
        certify_band(removed);
    }

    // Match each component AND each hole separately. Equal totals alone would
    // allow a disappearing island to be hidden by a split elsewhere.
    const auto matching = [&](const ExPolygons& a, const ExPolygons& b, bool holes) {
        if (a.size() != b.size()) throw std::invalid_argument("Boundary uncertainty changes connected components.");
        std::vector<bool> used(b.size(), false);
        for (const auto& first : a) {
            std::size_t match = b.size();
            for (std::size_t j = 0; j < b.size(); ++j) {
                charge(ExPolygons{first}, ExPolygons{b[j]});
                if (!intersection_ex(ExPolygons{first}, ExPolygons{b[j]}).empty()) {
                    if (match != b.size() || used[j]) throw std::invalid_argument("Ambiguous boundary component correspondence.");
                    match = j;
                }
            }
            if (match == b.size()) throw std::invalid_argument("A boundary component disappears within the uncertainty band.");
            used[match] = true;
            if (!holes) continue;
            const auto& second = b[match];
            if (first.holes.size() != second.holes.size()) {
                const auto diagnostic = [](const Polygons& rings) {
                    std::vector<std::pair<double, const Polygon*>> areas;
                    for (const auto& ring : rings) areas.emplace_back(std::abs(ring.area()), &ring);
                    std::sort(areas.begin(), areas.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
                    std::string text = std::to_string(rings.size()) + " smallest_nm2=";
                    for (std::size_t i = 0; i < std::min<std::size_t>(areas.size(), 8); ++i) {
                        const auto box = get_extents(areas[i].second->points);
                        text += std::to_string(areas[i].first) + "[" + std::to_string(box.min.x()) + "," +
                            std::to_string(box.min.y()) + ":" + std::to_string(box.max.x()) + "," +
                            std::to_string(box.max.y()) + "],";
                    }
                    return text;
                };
                throw std::invalid_argument("Boundary uncertainty changes holes: before=" + diagnostic(first.holes) +
                    " after=" + diagnostic(second.holes) + ".");
            }
            std::vector<bool> hole_used(second.holes.size(), false);
            for (const auto& source_hole : first.holes) {
                auto p = source_hole;
                p.reverse();
                std::size_t found = second.holes.size();
                for (std::size_t j = 0; j < second.holes.size(); ++j) {
                    auto q = second.holes[j];
                    q.reverse();
                    charge(ExPolygons{ExPolygon(p)}, ExPolygons{ExPolygon(q)});
                    if (!intersection_ex(ExPolygons{ExPolygon(p)}, ExPolygons{ExPolygon(q)}).empty()) {
                        if (found != second.holes.size() || hole_used[j])
                            throw std::invalid_argument("Ambiguous boundary hole correspondence.");
                        found = j;
                    }
                }
                if (found == second.holes.size()) throw std::invalid_argument("A hole disappears within the uncertainty band.");
                hole_used[found] = true;
            }
        }
    };
    if (!exactly_changed) {
        try {
            matching(before, after, true);
        } catch (const std::invalid_argument&) {
            if (callback_failed) throw;
            // Integer Boolean output may close a sub-grid channel even though
            // the original input unions preserve it. Prove correspondence on
            // those exact unions, not just equal aggregate component/hole counts.
            // The boundary band, surviving-feature and material-contact checks
            // still apply; this does not repair or replace any sampled geometry.
            if (result.work == max_work) throw ResourceLimit("Exact boundary topology work budget exhausted.");
            result.work += qualify_exact_section_topology(original, represented, max_work - result.work, checkpoint);
        }
    }
    // Euclidean allowance: offset chord sag <= .25 and vertex rounding <=
    // sqrt(.5), followed by two Boolean vertex roundings <= sqrt(.5) each.
    // Their sum is <3 native units. Two further units protect the morphological
    // test's own offset and intersection rounding; they are not geometric error.
    const unsigned allowance = band + 3, margin = allowance + 2;
    auto topology_sets = material_regions;
    topology_sets.push_back(&before);
    topology_sets.push_back(&after);
    for (std::size_t set_index = 0; set_index < topology_sets.size(); ++set_index) {
        const auto* polygons = topology_sets[set_index];
        if (!polygons) throw std::invalid_argument("Missing material boundary geometry.");
        charge(*polygons);
        const auto eroded = bounded(offset_ex(*polygons, -float(margin), jtRound, .25));
        const auto dilated = bounded(offset_ex(*polygons, float(margin), jtRound, .25));
        const auto check_robustness = [&](const ExPolygons& perturbed, const char* operation) {
            try {
                matching(*polygons, perturbed, true);
            } catch (const std::invalid_argument& error) {
                if (callback_failed) throw;
                // The robustness reference is the original input union, not
                // a rounded union that may already have an artificial portal.
                const auto& reference = set_index == material_regions.size() ? original :
                    set_index == material_regions.size() + 1 ? represented : *polygons;
                const auto checked_checkpoint = [&] {
                    if (checkpoint) try { checkpoint(); } catch (...) { callback_failed=true; throw; }
                };
                std::ostringstream diagnostic;
                if (set_index>=material_regions.size()) {
                    charge(*polygons);
                    const Polygon* smallest=nullptr;
                    for (const auto& polygon : *polygons)
                        for (const auto& hole : polygon.holes)
                            if (!smallest || std::abs(hole.area())<std::abs(smallest->area())) smallest=&hole;
                    if (smallest) {
                        if (result.work==max_work) throw ResourceLimit("Exact hole diagnostic work budget exhausted.");
                        // Compute before the fallback: a failed exact topology
                        // call cannot return its consumed work for another query.
                        // This is evidence only; it never changes admission.
                        const auto detail=diagnose_exact_hole_change(original,represented,*smallest,
                            std::min<std::size_t>(4000000,max_work-result.work),checked_checkpoint);
                        result.work+=detail.work;
                        diagnostic.precision(17);
                        diagnostic << " Smallest nominal occupied-hole diagnostic: source_unique=" << detail.source_hole_unique
                            << ", represented_unique=" << detail.represented_hole_unique
                            << ", correspondence_unique=" << detail.correspondence_unique
                            << ", hole_unchanged=" << (detail.hole_unchanged ? (*detail.hole_unchanged ? "true" : "false") : "unknown")
                            << ", difference_empty=" << (detail.difference_empty ? (*detail.difference_empty ? "true" : "false") : "unknown");
                        if (detail.change_distance_squared_nm2)
                            diagnostic << ", change_distance_squared_nm2=[" << (*detail.change_distance_squared_nm2)[0]
                                << ',' << (*detail.change_distance_squared_nm2)[1] << ']';
                        else diagnostic << ", change_distance_squared_nm2=unknown";
                        diagnostic << ", diagnostic_work=" << detail.work << '.';
                    }
                }
                if (result.work == max_work) throw ResourceLimit("Exact boundary robustness work budget exhausted.");
                try {
                    result.work += qualify_exact_section_topology(reference, perturbed, max_work - result.work, checked_checkpoint);
                } catch (const std::invalid_argument& exact_error) {
                    if (callback_failed) throw;
                    throw std::invalid_argument("Boundary robustness (set=" + std::to_string(set_index) +
                        ", operation=" + operation + ", margin_nm=" + std::to_string(margin) + "): " +
                        error.what() + " Exact check: " + exact_error.what()+diagnostic.str());
                }
            }
        };
        check_robustness(eroded, "erosion");
        check_robustness(dilated, "dilation");
    }
    result.uncertainty_nm = allowance;
    charge(added, removed);
    if (!exactly_changed) {
        result.added_area_bound_m2 = coverage_area_bound_m2(added, checkpoint);
        result.removed_area_bound_m2 = coverage_area_bound_m2(removed, checkpoint);
    }
    return result;
}
} // namespace Slic3r::Thermoelastic
