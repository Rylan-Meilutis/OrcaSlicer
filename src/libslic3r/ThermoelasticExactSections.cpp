#include "ThermoelasticExactSections.hpp"
#include "ThermoelasticMesh.hpp"
#include "ExPolygon.hpp"
#include "Point.hpp"
#include "Polygon.hpp"
#include "BoundingBox.hpp"
#include "libslic3r.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Lazy.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_repair/repair.h>
#include <CGAL/Polygon_repair/Non_zero_rule.h>
#include <CGAL/Boolean_set_operations_2/Gps_polygon_validation.h>
#include <CGAL/enum.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/number_utils.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <list>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
namespace {
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Polygon2 = CGAL::Polygon_2<Kernel>;
using PolygonHoles = CGAL::Polygon_with_holes_2<Kernel>;
using Set = CGAL::Polygon_set_2<Kernel>;

Kernel::FT segment_distance_squared(const Kernel::Point_2& point, const Kernel::Point_2& p, const Kernel::Point_2& q)
{
    const auto dx = q.x()-p.x(), dy = q.y()-p.y();
    const auto x = point.x()-p.x(), y = point.y()-p.y();
    const auto length2 = dx*dx+dy*dy, projection = x*dx+y*dy;
    if (length2 == 0 || projection <= 0) return x*x+y*y;
    if (projection >= length2) {
        const auto qx = point.x()-q.x(), qy = point.y()-q.y();
        return qx*qx+qy*qy;
    }
    const auto cross = dx*y-dy*x;
    return cross*cross/length2;
}

std::vector<PolygonHoles> convert_polygons(const ExPolygons& polygons, const std::function<void()>& checkpoint,
    const char* stage = "predicate")
{
    std::vector<PolygonHoles> result;
    result.reserve(polygons.size());
    std::size_t boundary_index = 0;
    const auto ring = [&](const Polygon& boundary, bool hole) {
        Polygon2 output;
        for (const auto& p : boundary.points) {
            if (checkpoint) checkpoint();
            output.push_back(Kernel::Point_2(double(p.x()), double(p.y())));
        }
        const bool simple = output.is_simple();
        if (!simple || output.orientation() != (hole ? CGAL::CLOCKWISE : CGAL::COUNTERCLOCKWISE)) {
            std::set<std::pair<coord_t, coord_t>> distinct;
            std::size_t duplicates = 0, zero_edges = 0, collinear = 0;
            std::string crossing;
            for (std::size_t i = 0; i < boundary.points.size(); ++i) {
                if (checkpoint) checkpoint();
                const auto& p = boundary.points[i];
                if (!distinct.emplace(p.x(), p.y()).second) ++duplicates;
                if (output[i] == output[(i + 1) % output.size()]) ++zero_edges;
                if (Kernel().orientation_2_object()(output[i], output[(i + 1) % output.size()],
                    output[(i + 2) % output.size()]) == CGAL::COLLINEAR) ++collinear;
            }
            // Only on rejection: identify an exact non-adjacent edge contact.
            // The caller has already reserved quadratic boundary work. Keep
            // this diagnostic bounded and preserve the original geometry.
            for (std::size_t i = 0; !simple && crossing.empty() && i < output.size(); ++i) {
                for (std::size_t j = i + 2; j < output.size(); ++j) {
                    if (checkpoint) checkpoint();
                    if (i == 0 && j + 1 == output.size()) continue;
                    if (!Kernel().do_intersect_2_object()(Kernel::Segment_2(output[i],output[(i+1)%output.size()]),
                        Kernel::Segment_2(output[j],output[(j+1)%output.size()]))) continue;
                    crossing = ", crossing_edges=" + std::to_string(i) + "," + std::to_string(j);
                    for (const auto index : {i,(i+1)%output.size(),j,(j+1)%output.size()}) {
                        const auto& p = boundary.points[index];
                        crossing += " (" + std::to_string(p.x()) + "," + std::to_string(p.y()) + ")";
                    }
                    break;
                }
            }
            throw std::invalid_argument(std::string("Exact section requires simple oriented rings: stage=") + stage +
                ", ring=" + std::to_string(boundary_index) + ", vertices=" + std::to_string(boundary.points.size()) +
                ", simple=" + (simple ? "true" : "false") + ", duplicates=" + std::to_string(duplicates) +
                ", zero_edges=" + std::to_string(zero_edges) + ", collinear_corners=" + std::to_string(collinear) + crossing + ".");
        }
        ++boundary_index;
        return output;
    };
    for (const auto& polygon : polygons) {
        auto outer = ring(polygon.contour, false);
        std::vector<Polygon2> holes;
        holes.reserve(polygon.holes.size());
        for (const auto& hole : polygon.holes) holes.push_back(ring(hole, true));
        PolygonHoles converted(outer, holes.begin(), holes.end());
        if (!CGAL::is_valid_polygon_with_holes(converted, Set::Traits_2{}))
            throw std::invalid_argument("Exact section polygon has invalid hole placement.");
        result.push_back(std::move(converted));
    }
    return result;
}
std::vector<PolygonHoles> repair_native_boundaries(const std::vector<Polygon2>& boundaries,
    std::size_t max_work, const std::function<void()>& checkpoint, std::size_t& work)
{
    const auto charge = [&](std::size_t amount) {
        if (checkpoint) checkpoint();
        if (amount > max_work-work) throw ResourceLimit("Exact native winding construction work exceeded.");
        work += amount;
    };
    std::size_t count=0;
    std::vector<Kernel::Segment_2> edges;
    for (const auto& ring : boundaries) {
        if (ring.size()>1024-count) throw ResourceLimit("Exact native winding input vertex budget exceeded.");
        count+=ring.size();
        for (std::size_t i=0;i<ring.size();++i) {
            if (checkpoint) checkpoint();
            const auto& a=ring[i];
            const auto& b=ring[(i+1)%ring.size()];
            if (a!=b) edges.emplace_back(a,b);
        }
    }
    bool planar=false;
    for (const auto& edge : edges) {
        charge(1);
        const auto orientation=Kernel().orientation_2_object();
        if (orientation(edges.front().source(),edges.front().target(),edge.source())!=CGAL::COLLINEAR ||
            orientation(edges.front().source(),edges.front().target(),edge.target())!=CGAL::COLLINEAR) {
            planar=true;
            break;
        }
    }
    // Closed chains supported only on a line have no occupied area.
    if (!planar) return {};
    std::size_t vertex_bound=count;
    for (std::size_t i=0;i<edges.size();++i)
        for (std::size_t j=i+1;j<edges.size();++j) {
            charge(1);
            const auto& a=edges[i];
            const auto& b=edges[j];
            if (a.source()==b.source() || a.source()==b.target() ||
                a.target()==b.source() || a.target()==b.target()) continue;
            const auto orientation=Kernel().orientation_2_object();
            if (orientation(a.source(),a.target(),b.source())==CGAL::COLLINEAR &&
                orientation(a.source(),a.target(),b.target())==CGAL::COLLINEAR) continue;
            if (Kernel().do_intersect_2_object()(a,b) && ++vertex_bound>4096)
                throw ResourceLimit("Exact native winding intersection vertex budget exceeded.");
        }
    // Bound the noninterruptible exact repair before entering CGAL.
    charge(2*vertex_bound*vertex_bound+vertex_bound+1);
    const PolygonHoles original(boundaries.front(),std::next(boundaries.begin()),boundaries.end());
    const auto filled=CGAL::Polygon_repair::repair(original,CGAL::Polygon_repair::Non_zero_rule{});
    if (checkpoint) checkpoint();
    std::vector<PolygonHoles> converted;
    std::size_t output_points=0;
    for (const auto& polygon : filled.polygons_with_holes()) {
        output_points+=polygon.outer_boundary().size();
        for (const auto& hole : polygon.holes()) output_points+=hole.size();
        if (output_points>4096) throw ResourceLimit("Exact native winding boundary budget exceeded.");
        converted.push_back(polygon);
    }
    charge(output_points*output_points+1);
    return converted;
}

template<class Charge>
void certify_boundary_band(const Set& changed, const std::list<PolygonHoles>& native, unsigned radius_nm,
                           const std::string& stage, const Charge& charge)
{
    charge(changed.arrangement().number_of_vertices()+1);
    std::list<PolygonHoles> changed_pieces;
    changed.polygons_with_holes(std::back_inserter(changed_pieces));
    const auto distance_squared = [&](const Kernel::Point_2& point, const Kernel::Point_2& p, const Kernel::Point_2& q) {
        charge(1);
        return segment_distance_squared(point,p,q);
    };
        for (const auto& piece : changed_pieces) {
            const auto covered = [&](const auto& vertices) {
                bool certified = false;
                for (const auto& reference : native) {
                    const auto ring = [&](const Polygon2& boundary) {
                        for (std::size_t i = 0; i < boundary.size() && !certified; ++i) {
                            bool contains = true;
                            for (const auto& point : vertices)
                                if (distance_squared(point, boundary[i], boundary[(i + 1) % boundary.size()]) >
                                    Kernel::FT(radius_nm * radius_nm)) {
                                    contains = false;
                                    break;
                                }
                            certified = contains;
                        }
                    };
                    ring(reference.outer_boundary());
                    for (auto h = reference.holes_begin(); h != reference.holes_end() && !certified; ++h) ring(*h);
                    if (certified) break;
                }
                return certified;
            };
            bool certified = covered(piece.outer_boundary());
            if (!certified) {
                const auto& outer = piece.outer_boundary();
                charge(outer.size() * outer.size() + 1);
                struct Triangle { std::array<Kernel::Point_2, 3> p; unsigned depth = 0; };
                std::vector<Triangle> pending;
                const auto orientation = Kernel().orientation_2_object();
                std::vector<Polygon2> cycles;
                if (outer.is_simple()) cycles.push_back(outer);
                else {
                    std::map<Kernel::Point_2, std::size_t> positions;
                    std::vector<Kernel::Point_2> path;
                    for (std::size_t i = 0; i <= outer.size(); ++i) {
                        charge(1);
                        const auto& p = outer[i % outer.size()];
                        const auto found = positions.find(p);
                        if (found == positions.end()) {
                            positions.emplace(p, path.size());
                            path.push_back(p);
                        } else {
                            const auto begin = found->second;
                            cycles.emplace_back(path.begin() + std::ptrdiff_t(begin), path.end());
                            for (std::size_t j = begin + 1; j < path.size(); ++j) {
                                charge(1);
                                positions.erase(path[j]);
                            }
                            path.resize(begin + 1);
                        }
                    }
                    if (path.size() != 1) throw std::invalid_argument("Analysis interface cycle decomposition is incomplete.");
                    using Edge = std::pair<Kernel::Point_2, Kernel::Point_2>;
                    std::map<Edge, std::size_t> original_edges, cycle_edges;
                    for (std::size_t i = 0; i < outer.size(); ++i) {
                        charge(1);
                        ++original_edges[{outer[i], outer[(i + 1) % outer.size()]}];
                    }
                    for (const auto& cycle : cycles)
                        for (std::size_t i = 0; i < cycle.size(); ++i) {
                            charge(1);
                            ++cycle_edges[{cycle[i], cycle[(i + 1) % cycle.size()]}];
                        }
                    if (original_edges != cycle_edges) throw std::invalid_argument("Analysis interface cycle edges changed.");
                }
                for (auto& cycle : cycles) {
                    if (cycle.size() < 3) continue; // Closed point/segment has no regularized area.
                    bool collinear = true;
                    for (std::size_t i = 2; i < cycle.size(); ++i) {
                        charge(1);
                        if (orientation(cycle[0], cycle[1], cycle[i]) != CGAL::COLLINEAR) collinear = false;
                    }
                    if (collinear) continue;
                    if (!cycle.is_simple()) throw std::invalid_argument("Analysis interface cycle has a genuine crossing.");
                    // Preserve directed edges above before reversing. Filling
                    // every signed cycle is conservative, including holes.
                    if (cycle.orientation() == CGAL::CLOCKWISE) cycle.reverse_orientation();
                    if (cycle.is_convex()) {
                        for (std::size_t i = 1; i + 1 < cycle.size(); ++i) {
                            charge(1);
                            if (orientation(cycle[0], cycle[i], cycle[i + 1]) != CGAL::COLLINEAR)
                                pending.push_back({{cycle[0], cycle[i], cycle[i + 1]}, 0});
                        }
                        continue;
                    }
                    // Bounded exact ear clipping of the simple outer ring;
                    // filling holes is conservative for this coverage proof.
                    std::vector<Kernel::Point_2> ring(cycle.vertices_begin(), cycle.vertices_end());
                    while (ring.size() > 3) {
                        bool removed = false;
                        for (std::size_t i = 0; i < ring.size() && !removed; ++i) {
                            const auto prev = (i + ring.size() - 1) % ring.size(), next = (i + 1) % ring.size();
                            const auto& a = ring[prev]; const auto& b = ring[i]; const auto& c = ring[next];
                            charge(1);
                            const auto turn = orientation(a, b, c);
                            if (turn == CGAL::RIGHT_TURN) continue;
                            if (turn == CGAL::COLLINEAR) {
                                if ((b.x() - a.x()) * (b.x() - c.x()) + (b.y() - a.y()) * (b.y() - c.y()) > 0) continue;
                            } else {
                                bool contains = false;
                                for (std::size_t j = 0; j < ring.size() && !contains; ++j) {
                                    if (j == prev || j == i || j == next) continue;
                                    charge(1);
                                    contains = orientation(a, b, ring[j]) != CGAL::RIGHT_TURN &&
                                        orientation(b, c, ring[j]) != CGAL::RIGHT_TURN && orientation(c, a, ring[j]) != CGAL::RIGHT_TURN;
                                }
                                if (contains) continue;
                                pending.push_back({{a, b, c}, 0});
                            }
                            charge(ring.size());
                            ring.erase(ring.begin() + std::ptrdiff_t(i));
                            removed = true;
                        }
                        if (!removed) throw std::invalid_argument("Analysis interface exact triangulation is unresolved.");
                    }
                    if (orientation(ring[0], ring[1], ring[2]) != CGAL::COLLINEAR)
                        pending.push_back({{ring[0], ring[1], ring[2]}, 0});
                }
                certified = true;
                while (!pending.empty()) {
                    charge(1);
                    const auto triangle = pending.back();
                    pending.pop_back();
                    if (covered(triangle.p)) continue;
                    // Long thin ribbons may straddle adjacent source-edge
                    // capsules. Longest-edge bisection can require two levels
                    // per spatial halving; exact rational midpoints permit
                    // this depth without reducing the independent work cap.
                    if (triangle.depth == 48) { certified = false; break; }
                    std::size_t edge = 0;
                    Kernel::FT longest = -1;
                    for (std::size_t i = 0; i < 3; ++i) {
                        const auto dx = triangle.p[i].x() - triangle.p[(i + 1) % 3].x();
                        const auto dy = triangle.p[i].y() - triangle.p[(i + 1) % 3].y();
                        const Kernel::FT length2 = dx * dx + dy * dy;
                        if (length2 > longest) { longest = length2; edge = i; }
                    }
                    const auto& a = triangle.p[edge]; const auto& b = triangle.p[(edge + 1) % 3];
                    const auto& c = triangle.p[(edge + 2) % 3];
                    const Kernel::Point_2 midpoint((a.x() + b.x()) / 2, (a.y() + b.y()) / 2);
                    pending.push_back({{a, midpoint, c}, triangle.depth + 1});
                    pending.push_back({{midpoint, b, c}, triangle.depth + 1});
                }
            }
            // A capsule is convex: containing every outer vertex contains the
            // complete component, including a conservative fill of its holes.
            if (!certified) {
                double maximum_distance = 0.;
                for (const auto& point : piece.outer_boundary()) {
                    std::optional<Kernel::FT> nearest;
                    for (const auto& reference : native) {
                        const auto ring = [&](const Polygon2& boundary) {
                            for (std::size_t i = 0; i < boundary.size(); ++i) {
                                const auto distance = distance_squared(point, boundary[i], boundary[(i + 1) % boundary.size()]);
                                if (!nearest || distance < *nearest) nearest = distance;
                            }
                        };
                        ring(reference.outer_boundary());
                        for (auto h = reference.holes_begin(); h != reference.holes_end(); ++h) ring(*h);
                    }
                    if (nearest) maximum_distance = std::max(maximum_distance, std::sqrt(CGAL::to_double(*nearest)));
                }
                throw std::invalid_argument("Analysis interface band is unverified: " + stage +
                    ", component_vertices=" + std::to_string(piece.outer_boundary().size()) +
                    ", convex=" + (piece.outer_boundary().is_convex() ? "true" : "false") +
                    ", diagnostic_max_nearest_vertex_nm=" + std::to_string(maximum_distance) + ".");
            }
        }

}
Set make_topology_union(const ExPolygons& polygons, std::size_t max_work,
    const std::function<void()>& checkpoint, std::size_t& work, bool native_nonzero_fill = false)
{
    if (!max_work || max_work > 4000000) throw std::invalid_argument("Invalid exact topology diagnostic budget.");
    std::size_t count = 0;
    if (polygons.size() > 2048) throw ResourceLimit("Exact topology diagnostic polygon limit exceeded.");
    for (const auto& polygon : polygons) {
        if (polygon.holes.size() > 2048) throw ResourceLimit("Exact topology diagnostic hole limit exceeded.");
        const auto preflight = [&](const Polygon& ring) {
            if (ring.points.size() < 3) throw std::invalid_argument("Invalid exact topology diagnostic boundary.");
            if (ring.points.size() > 2048 - count) throw ResourceLimit("Exact topology diagnostic vertex limit exceeded.");
            count += ring.points.size();
            for (const auto& p : ring.points) {
                if (checkpoint) checkpoint();
                if (std::abs(double(p.x())) > 1e10 || std::abs(double(p.y())) > 1e10)
                    throw std::invalid_argument("Exact topology diagnostic coordinate outside supported range.");
            }
        };
        preflight(polygon.contour);
        for (const auto& hole : polygon.holes) preflight(hole);
    }
    work = count * count + count + 1;
    if (work > max_work) throw ResourceLimit("Exact topology diagnostic work limit exceeded (vertices=" +
        std::to_string(count) + ", requested=" + std::to_string(work) + ", available=" + std::to_string(max_work) + ").");
    std::vector<PolygonHoles> converted;
    bool callback_failed = false;
    const auto validation_checkpoint = [&] {
        if (checkpoint) {
            try { checkpoint(); }
            catch (...) { callback_failed = true; throw; }
        }
    };
    try {
        converted = convert_polygons(polygons,validation_checkpoint);
    } catch (const std::invalid_argument&) {
        if (!native_nonzero_fill || callback_failed) throw;
        // Bound the exceptional exact-intersection triangulation separately.
        // This fill rule defines the native occupied set; it is not a tolerance
        // repair or permission to round rational intersections to the nm grid.
        if (count > 256) {
            // Disjoint bounding-box components cannot cancel one another's
            // winding. Process them independently instead of triangulating a
            // whole mask because one small component has a crossing.
            if (2*count*count+count+1 > max_work)
                throw ResourceLimit("Exact native winding component budget exceeded.");
            work = 2*count*count+count+1;
            std::vector<BoundingBox> boxes;
            for (const auto& polygon : polygons) {
                BoundingBox box(polygon.contour.points);
                for (const auto& hole : polygon.holes) for (const auto& p : hole.points) box.merge(p);
                boxes.push_back(box);
            }
            std::vector<bool> assigned(polygons.size(),false);
            std::vector<ExPolygons> components;
            for (std::size_t seed = 0; seed < polygons.size(); ++seed) {
                if (assigned[seed]) continue;
                assigned[seed] = true;
                std::vector<std::size_t> connected{seed};
                for (std::size_t cursor = 0; cursor < connected.size(); ++cursor) {
                    if (checkpoint) checkpoint();
                    for (std::size_t i = 0; i < polygons.size(); ++i)
                        if (!assigned[i] && boxes[connected[cursor]].overlap(boxes[i])) {
                            assigned[i] = true;
                            connected.push_back(i);
                        }
                }
                auto& component = components.emplace_back();
                for (const auto i : connected) component.push_back(polygons[i]);
            }
            if (components.size() > 1) {
                Set occupied;
                for (const auto& component : components) {
                    if (work == max_work) throw ResourceLimit("Exact native winding component work exceeded.");
                    std::size_t used = 0;
                    const auto part = make_topology_union(component,max_work-work,checkpoint,used,true);
                    work += used;
                    const auto vertices = occupied.arrangement().number_of_vertices()+part.arrangement().number_of_vertices();
                    if (vertices > 4096) throw ResourceLimit("Exact native winding component boundary budget exceeded.");
                    const auto join_work = vertices*vertices+1;
                    if (join_work > max_work-work) throw ResourceLimit("Exact native winding component union work exceeded.");
                    work += join_work;
                    occupied.join(part);
                    if (checkpoint) checkpoint();
                }
                return occupied;
            }
        }
        const auto ring = [&](const Polygon& input) {
            Polygon2 output;
            for (const auto& point : input.points) {
                if (checkpoint) checkpoint();
                output.push_back(Kernel::Point_2(double(point.x()),double(point.y())));
            }
            return output;
        };
        std::vector<Polygon2> boundaries;
        for (std::size_t i = 0; i < polygons.size(); ++i) {
            boundaries.push_back(ring(polygons[i].contour));
            for (const auto& hole : polygons[i].holes) boundaries.push_back(ring(hole));
        }
        // Nonzero reconstruction uses the ORIGINAL orientation of every
        // boundary, including other contours. Packing them as boundary chains
        // in one input makes cancellation global, matching native clipping;
        // no chain is reversed or independently repaired and re-unioned.
        converted=repair_native_boundaries(boundaries,max_work,checkpoint,work);
    }
    Set occupied;
    occupied.join(converted.begin(), converted.end());
    if (checkpoint) checkpoint();
    if (occupied.arrangement().number_of_vertices() > 4096)
        throw ResourceLimit("Exact topology arrangement limit exceeded.");
    return occupied;
}

Set make_clipped_native_union(const ExPolygons& polygons, const std::array<Kernel::FT,4>& box,
    std::size_t max_work, const std::function<void()>& checkpoint, std::size_t& work)
{
    work=0;
    const auto charge = [&](std::size_t amount) {
        if (checkpoint) checkpoint();
        if (amount>max_work-work) throw ResourceLimit("Exact native clipping work budget exceeded.");
        work+=amount;
    };
    std::vector<Polygon2> boundaries;
    std::size_t retained=0;
    const auto clip_ring = [&](const Polygon& ring) {
        charge(ring.points.size()+1);
        if (ring.points.empty()) return;
        std::vector<Kernel::Point_2> chain;
        // Stream the original chain through the first plane: remote input
        // complexity must not allocate an unbounded rational-point array.
        for (unsigned side=0;side<4;++side) {
            if (side && chain.empty()) break;
            const unsigned axis=side/2;
            const bool lower=side%2==0;
            const auto& boundary=box[side];
            const auto coordinate = [&](const Kernel::Point_2& p) { return axis ? p.y() : p.x(); };
            const auto inside = [&](const Kernel::Point_2& p) {
                return lower ? coordinate(p)>=boundary : coordinate(p)<=boundary;
            };
            std::vector<Kernel::Point_2> clipped;
            const auto append = [&](const Kernel::Point_2& p) {
                if (clipped.empty() || clipped.back()!=p) {
                    if (clipped.size()==4096) throw ResourceLimit("Exact native clipping boundary budget exceeded.");
                    clipped.push_back(p);
                }
            };
            auto previous=side ? chain.back() : Kernel::Point_2(double(ring.points.back().x()),double(ring.points.back().y()));
            bool previous_inside=inside(previous);
            const auto consume = [&](const Kernel::Point_2& point) {
                charge(1);
                const bool point_inside=inside(point);
                if (point_inside!=previous_inside) {
                    const Kernel::FT t=(boundary-coordinate(previous))/(coordinate(point)-coordinate(previous));
                    const auto crossing=previous+t*(point-previous);
                    append(crossing);
                }
                if (point_inside) append(point);
                previous=point;
                previous_inside=point_inside;
            };
            if (side) for (const auto& point : chain) consume(point);
            else for (const auto& point : ring.points) consume(Kernel::Point_2(double(point.x()),double(point.y())));
            if (clipped.size()>1 && clipped.front()==clipped.back()) clipped.pop_back();
            chain=std::move(clipped);
        }
        if (chain.size()<3) return;
        if (chain.size()>1024-retained) throw ResourceLimit("Exact cropped native winding vertex budget exceeded.");
        retained+=chain.size();
        boundaries.emplace_back(chain.begin(),chain.end());
    };
    // Sutherland-Hodgman preserves each oriented chain's winding inside the
    // rectangle, even when it produces retraced connectors for concave input.
    // Keep rational crossings and all original orientations, then resolve the
    // combined nonzero winding once. Never round or repair each ring alone.
    for (const auto& polygon : polygons) {
        clip_ring(polygon.contour);
        for (const auto& hole : polygon.holes) clip_ring(hole);
    }
    charge(retained*retained+retained+1);
    const auto converted=repair_native_boundaries(boundaries,max_work,checkpoint,work);
    Set occupied;
    occupied.join(converted.begin(),converted.end());
    if (checkpoint) checkpoint();
    if (occupied.arrangement().number_of_vertices()>4096)
        throw ResourceLimit("Exact cropped native arrangement limit exceeded.");
    return occupied;
}
} // namespace
struct ExactSection::Impl {
    Set geometry;
    std::size_t vertices = 0;
    std::optional<BoundingBox> enclosure;
};
class ExactSectionAccess {
public:
    static const Set& set(const ExactSection& section)
    {
        static const Set empty;
        return section.m_impl ? section.m_impl->geometry : empty;
    }
    static ExactSection make(Set geometry, const std::function<void(std::size_t)>& charge)
    {
        const auto vertices = geometry.arrangement().number_of_vertices();
        if (vertices > 4096) throw ResourceLimit("Exact section retained vertex budget exceeded.");
        charge(vertices + 1);
        if (geometry.is_empty()) return {};
        auto impl = std::make_shared<ExactSection::Impl>();
        auto vertex = geometry.arrangement().vertices_begin();
        auto bounds = vertex->point().bbox();
        for (; vertex != geometry.arrangement().vertices_end(); ++vertex) {
            charge(1);
            bounds = bounds + vertex->point().bbox();
        }
        if (!std::isfinite(bounds.xmin()) || !std::isfinite(bounds.ymin()) ||
            !std::isfinite(bounds.xmax()) || !std::isfinite(bounds.ymax()) ||
            bounds.xmin() < -1e10 || bounds.ymin() < -1e10 || bounds.xmax() > 1e10 || bounds.ymax() > 1e10)
            throw std::invalid_argument("Exact section bounds exceed the coordinate limit.");
        impl->enclosure = BoundingBox(Point(coord_t(std::floor(bounds.xmin())), coord_t(std::floor(bounds.ymin()))),
            Point(coord_t(std::ceil(bounds.xmax())), coord_t(std::ceil(bounds.ymax()))));
        impl->vertices = vertices;
        impl->geometry = std::move(geometry);
        return ExactSection(std::move(impl));
    }
};
bool ExactSection::empty() const { return !m_impl; }
std::size_t ExactSection::vertex_count() const { return m_impl ? m_impl->vertices : 0; }
std::optional<BoundingBox> ExactSection::bounds() const { return m_impl ? m_impl->enclosure : std::nullopt; }

namespace {
class ExactSectionWork {
public:
    ExactSectionWork(std::size_t maximum, const std::function<void()>& checkpoint, const char* operation)
        : m_maximum(maximum), m_checkpoint(checkpoint), m_operation(operation)
    {
        if (!maximum || maximum > 4000000) throw std::invalid_argument("Invalid exact section work budget.");
    }
    // The callback is borrowed for this operation; reject temporary wrappers.
    ExactSectionWork(std::size_t, std::function<void()>&&, const char*) = delete;
    void charge(std::size_t amount)
    {
        if (m_checkpoint) m_checkpoint();
        if (amount > remaining())
            throw ResourceLimit(std::string("Exact section work budget exceeded: operation=") + m_operation +
                ", used=" + std::to_string(used) + ", requested=" + std::to_string(amount) +
                ", limit=" + std::to_string(m_maximum) + ", per_operation_cap=4000000.");
        used += amount;
    }
    std::size_t remaining() const { return m_maximum - used; }
    std::size_t used = 0;
private:
    std::size_t m_maximum;
    const std::function<void()>& m_checkpoint;
    const char* m_operation;
};
std::array<double,2> certified_nonnegative_root(const Kernel::FT& squared, double estimate,
    const std::function<void(std::size_t)>& charge)
{
    using FT=Kernel::FT;
    if (squared<0 || !std::isfinite(estimate) || estimate<0)
        throw std::invalid_argument("Invalid exact square-root estimate.");
    double low=estimate,high=estimate;
    bool low_valid=false,high_valid=false;
    // Floating sqrt only seeds the search. Exact squared comparisons certify
    // both endpoints. Lazy-number cancellation can put the seed many ULPs
    // away, so adjacent-double refinement is only the fast path.
    for (unsigned pass=0;pass<2;++pass) {
        if (pass) {
            charge(1);
            // Refresh from the exact rational, not the lazy interval that
            // supplied the original seed. Under/overflow still uses bisection.
            const double exact_square=CGAL::to_double(CGAL::exact(squared));
            if (!(exact_square>0.) || !std::isfinite(exact_square)) break;
            low=high=std::sqrt(exact_square);
            low_valid=high_valid=false;
        }
        for (unsigned attempt=0;attempt<64;++attempt) {
            charge(1);
            if (FT(low)*FT(low)<=squared) { low_valid=true; break; }
            low=std::nextafter(low,0.);
        }
        for (unsigned attempt=0;attempt<64;++attempt) {
            charge(1);
            if (FT(high)*FT(high)>=squared) { high_valid=true; break; }
            high=std::nextafter(high,std::numeric_limits<double>::infinity());
            if (!std::isfinite(high)) break;
        }
        if (low_valid && high_valid) break;
    }
    if (!low_valid || !high_valid) {
        low=low_valid ? low : 0.;
        high=high_valid ? high : std::numeric_limits<double>::max();
        charge(1);
        if (FT(high)*FT(high)<squared)
            throw std::invalid_argument("Unresolved exact square-root bounds.");
    }
    // Tighten even successful fast-path bounds: a valid bracket can still be
    // tens of ULPs wide when the original lazy seed suffered cancellation.
    // Halving spans the complete exponent range (including subnormals) in the
    // fallback, but takes only a few iterations for a fast-path bracket.
    constexpr int attempts=std::numeric_limits<double>::max_exponent-
        std::numeric_limits<double>::min_exponent+std::numeric_limits<double>::digits+4;
    for (int attempt=0;attempt<attempts;++attempt) {
        charge(1);
        if (std::nextafter(low,high)==high) return {{low,high}};
        const double middle=low+(high-low)/2.;
        const FT square=FT(middle)*FT(middle);
        if (square==squared) return {{middle,middle}};
        if (square<squared) low=middle;
        else high=middle;
    }
    throw std::invalid_argument("Unresolved exact square-root bounds.");
}
void validate_exact_operation(ExactSectionOperation operation)
{
    if (operation != ExactSectionOperation::Intersection && operation != ExactSectionOperation::Difference &&
        operation != ExactSectionOperation::SymmetricDifference && operation != ExactSectionOperation::Union)
        throw std::invalid_argument("Invalid exact section operation.");
}
void reserve_exact_overlay_work(const Set& a, const Set& b,
    const std::function<void(std::size_t)>& charge, ExactOverlayDiagnostics* diagnostics = nullptr)
{
    const auto record=[&](std::size_t ExactOverlayDiagnostics::* field,std::size_t amount) {
        if (!diagnostics) return;
        auto& value=diagnostics->*field;
        if (amount>std::numeric_limits<std::size_t>::max()-value) {
            value=std::numeric_limits<std::size_t>::max();
            diagnostics->overflow=true;
        } else value+=amount;
    };
    const auto debit=[&](std::size_t amount,std::size_t ExactOverlayDiagnostics::* field) {
        charge(amount);
        record(field,amount);
    };
    // Callers enforce their existing arrangement caps before this preflight.
    const auto edges=a.arrangement().number_of_edges()+b.arrangement().number_of_edges();
    const auto vertices=a.arrangement().number_of_vertices()+b.arrangement().number_of_vertices();
    record(&ExactOverlayDiagnostics::overlays,1);
    record(&ExactOverlayDiagnostics::edges,edges);
    struct EdgeBox { CGAL::Bbox_2 box; unsigned operand; };
    std::vector<EdgeBox> boxes;
    boxes.reserve(edges);
    const auto collect = [&](const Set& geometry,unsigned operand) {
        for (auto edge=geometry.arrangement().edges_begin();edge!=geometry.arrangement().edges_end();++edge) {
            debit(1,&ExactOverlayDiagnostics::boxes_projection);
            boxes.push_back({edge->source()->point().bbox()+edge->target()->point().bbox(),operand});
        }
    };
    collect(a,0);
    collect(b,1);
    // Choose only the broad-phase projection, never the candidate predicate.
    // Summed projected edge lengths / total span estimates active-list density
    // in O(E), avoiding a quadratic trial sweep. Outward boxes and strict
    // expiry still retain every closed-box candidate on either chosen axis.
    unsigned axis=0;
    if (!boxes.empty()) {
        auto extent=boxes.front().box;
        std::array<double,2> lengths{};
        for (const auto& edge:boxes) {
            debit(1,&ExactOverlayDiagnostics::boxes_projection);
            extent=extent+edge.box;
            lengths[0]+=edge.box.xmax()-edge.box.xmin();
            lengths[1]+=edge.box.ymax()-edge.box.ymin();
        }
        debit(1,&ExactOverlayDiagnostics::boxes_projection);
        const double xspan=extent.xmax()-extent.xmin(), yspan=extent.ymax()-extent.ymin();
        if (std::isfinite(xspan) && std::isfinite(yspan) && std::isfinite(lengths[0]) &&
            std::isfinite(lengths[1]) && xspan>0. && yspan>0. && lengths[1]/yspan<lengths[0]/xspan) axis=1;
    }
    const auto lower=[&](const CGAL::Bbox_2& box) { return axis ? box.ymin() : box.xmin(); };
    const auto upper=[&](const CGAL::Bbox_2& box) { return axis ? box.ymax() : box.xmax(); };
    std::sort(boxes.begin(),boxes.end(),[&](const EdgeBox& x,const EdgeBox& y) {
        debit(1,&ExactOverlayDiagnostics::sorting);
        return lower(x.box)<lower(y.box);
    });
    std::array<std::vector<std::size_t>,2> active;
    std::size_t candidates=0;
    for (std::size_t i=0;i<boxes.size();++i) {
        debit(1,&ExactOverlayDiagnostics::sweep);
        const auto& current=boxes[i];
        auto& opposite=active[1-current.operand];
        for (std::size_t j=0;j<opposite.size();) {
            debit(1,&ExactOverlayDiagnostics::active_scans); // Includes expired entries and secondary-axis misses.
            const auto& previous=boxes[opposite[j]].box;
            // CLOSED outward boxes retain endpoint touches, vertical edges,
            // and both endpoints of any possible collinear overlap.
            if (upper(previous)<lower(current.box)) {
                opposite[j]=opposite.back();
                opposite.pop_back();
            } else {
                if (previous.xmin()<=current.box.xmax() && current.box.xmin()<=previous.xmax() &&
                    previous.ymin()<=current.box.ymax() && current.box.ymin()<=previous.ymax()) {
                    ++candidates;
                    record(&ExactOverlayDiagnostics::candidates,1);
                }
                ++j;
            }
        }
        active[current.operand].push_back(i);
    }
    // Each operand is already planar. A cross-operand segment pair introduces
    // at most two split endpoints. Reserve sweep/event complexity plus the
    // input copy; output traversal is charged by the caller. These are bounded
    // geometric-complexity units, NOT a CPU-time or rational-bit-cost guarantee.
    // Dense projected scans are charged above even if their other axis misses.
    // Every candidate consumed a scan charge under the callers'4M ceiling;
    // together with their bounded planar inputs this also bounds size_t
    // arithmetic here on32-bit platforms, before the reservation is tested.
    const auto events=edges+2*candidates+1;
    std::size_t levels=1;
    for (auto count=events;count>1;count=(count+1)/2) ++levels;
    debit(events*levels+edges+vertices+1,&ExactOverlayDiagnostics::event_reservation);
}
Set overlay_exact_sets(const Set& a, const Set& b, ExactSectionOperation operation,
    const std::function<void(std::size_t)>& charge, ExactOverlayDiagnostics* diagnostics = nullptr)
{
    validate_exact_operation(operation);
    if (a.arrangement().number_of_edges()+b.arrangement().number_of_edges()>4096)
        throw ResourceLimit("Exact section overlay edge budget exceeded.");
    reserve_exact_overlay_work(a,b,charge,diagnostics);
    // Zero candidates does not exclude containment: always do the exact
    // operation, never replace it with a boundary-only occupancy decision.
    Set output = a;
    switch (operation) {
    case ExactSectionOperation::Intersection: output.intersection(b); break;
    case ExactSectionOperation::Difference: output.difference(b); break;
    case ExactSectionOperation::SymmetricDifference: output.symmetric_difference(b); break;
    case ExactSectionOperation::Union: output.join(b); break;
    }
    charge(0);
    if (output.arrangement().number_of_vertices() > 4096)
        throw ResourceLimit("Exact section overlay vertex budget exceeded.");
    return output;
}
void count_exact_inputs(const ExPolygons& polygons, std::size_t& count,
    const std::function<void(std::size_t)>& charge)
{
    if (polygons.size() > 1024) throw ResourceLimit("Exact section input polygon budget exceeded.");
    for (const auto& polygon : polygons) {
        if (polygon.holes.size() > 1024) throw ResourceLimit("Exact section input ring budget exceeded.");
        const auto scan = [&](const Polygon& ring) {
            if (ring.points.size() > 1024 - count) throw ResourceLimit("Exact section input vertex budget exceeded.");
            count += ring.points.size();
            charge(ring.points.size() + 1);
        };
        scan(polygon.contour);
        for (const auto& hole : polygon.holes) scan(hole);
    }
}
} // namespace

ExactSectionResult make_exact_section(const ExPolygons& polygons, std::size_t max_work,
    const std::function<void()>& checkpoint)
{
    ExactSectionWork work(max_work, checkpoint, __func__);
    const auto charge = [&](std::size_t amount) { work.charge(amount); };
    std::size_t count = 0, used = 0;
    count_exact_inputs(polygons, count, charge);
    if (!work.remaining()) throw ResourceLimit("Exact section construction budget exhausted.");
    auto set = make_topology_union(polygons, work.remaining(), checkpoint, used, true);
    work.charge(used);
    auto section = ExactSectionAccess::make(std::move(set), charge);
    return {std::move(section), work.used};
}
ExactSectionResult boolean_exact_sections(const ExactSection& left, const ExactSection& right,
    ExactSectionOperation operation, std::size_t max_work, const std::function<void()>& checkpoint,
    ExactOverlayDiagnostics* diagnostics)
{
    ExactSectionWork work(max_work, checkpoint, __func__);
    const auto charge = [&](std::size_t amount) { work.charge(amount); };
    validate_exact_operation(operation);
    charge(0);
    const auto& a=ExactSectionAccess::set(left);
    const auto& b=ExactSectionAccess::set(right);
    if (a.arrangement().number_of_edges()+b.arrangement().number_of_edges()>4096)
        throw ResourceLimit("Exact section overlay edge budget exceeded.");
    // Algebraic identities on immutable canonical handles require metadata
    // only: no set copy, coordinate traversal, or Boolean overlay is performed.
    // Distinct nonempty handles require an exact boundary proof or overlay.
    if (&a==&b || left.empty() || right.empty()) {
        charge(1);
        ExactSection result;
        if (&a==&b) {
            if (operation==ExactSectionOperation::Union || operation==ExactSectionOperation::Intersection) result=left;
        } else if (left.empty()) {
            if (operation==ExactSectionOperation::Union || operation==ExactSectionOperation::SymmetricDifference) result=right;
        } else if (operation!=ExactSectionOperation::Intersection) result=left;
        return {std::move(result),work.used};
    }
    if (operation==ExactSectionOperation::Intersection || operation==ExactSectionOperation::Difference) {
        charge(2);
        const auto left_bounds=left.bounds(), right_bounds=right.bounds();
        // Immutable outward integer enclosures contain all rational fill.
        // Strict separation proves disjointness; touching boxes and all
        // containment cases still require the existing exact machinery.
        if (left_bounds && right_bounds &&
            (left_bounds->max.x()<right_bounds->min.x() || right_bounds->max.x()<left_bounds->min.x() ||
             left_bounds->max.y()<right_bounds->min.y() || right_bounds->max.y()<left_bounds->min.y()))
            return {operation==ExactSectionOperation::Difference ? left : ExactSection{},work.used};
    }
    // Independently sampled canonical sets often have identical directed
    // boundaries. Prove that identity before paying for an overlay. This is
    // deliberately stronger than geometric equality: differently subdivided
    // edges fall through, never rounded, merged, or treated as unequal sets.
    bool try_identity=a.arrangement().number_of_edges()==b.arrangement().number_of_edges();
    if (try_identity) {
        charge(2);
        const auto left_bounds=left.bounds(),right_bounds=right.bounds();
        // This only avoids an optional boundary-equality proof. A mismatch
        // still performs the full Boolean operation; matching bounds prove
        // nothing about occupied geometry and still require exact edges.
        try_identity=left_bounds && right_bounds && left_bounds->min==right_bounds->min &&
            left_bounds->max==right_bounds->max;
    }
    if (try_identity) {
        using Edge=std::pair<Kernel::Point_2,Kernel::Point_2>;
        const auto directed_edges = [&](const Set& geometry) {
            const auto count=geometry.arrangement().number_of_edges();
            charge(count+geometry.arrangement().number_of_vertices()+1);
            std::list<PolygonHoles> polygons;
            geometry.polygons_with_holes(std::back_inserter(polygons));
            std::vector<Edge> edges;
            edges.reserve(count);
            const auto ring = [&](const Polygon2& polygon) {
                for (std::size_t i=0;i<polygon.size();++i) {
                    charge(1);
                    if (edges.size()==count) throw ResourceLimit("Exact equality boundary edge limit exceeded.");
                    edges.emplace_back(polygon[i],polygon[(i+1)%polygon.size()]);
                }
            };
            for (const auto& polygon:polygons) {
                ring(polygon.outer_boundary());
                for (const auto& hole:polygon.holes()) ring(hole);
            }
            const auto less=Kernel().less_xy_2_object();
            std::sort(edges.begin(),edges.end(),[&](const Edge& x,const Edge& y) {
                charge(1);
                if (less(x.first,y.first)) return true;
                if (less(y.first,x.first)) return false;
                return less(x.second,y.second);
            });
            return edges;
        };
        const auto first=directed_edges(a),second=directed_edges(b);
        bool identical=first.size()==second.size();
        for (std::size_t i=0;identical && i<first.size();++i) {
            charge(1);
            identical=first[i].first==second[i].first && first[i].second==second[i].second;
        }
        if (identical) {
            charge(1);
            return {operation==ExactSectionOperation::Union || operation==ExactSectionOperation::Intersection ? left : ExactSection{},work.used};
        }
    }
    auto set = overlay_exact_sets(ExactSectionAccess::set(left), ExactSectionAccess::set(right), operation, charge, diagnostics);
    auto section = ExactSectionAccess::make(std::move(set), charge);
    return {std::move(section), work.used};
}

static ExactDepositedPartition first_deposition_partition_impl(std::size_t command_count,
    const std::function<ExactDepositedCommandResult(std::size_t,std::size_t)>& provider,
    ExactSectionWork& work,const std::function<void()>& checkpoint,bool streamed)
{
    ExactDepositedPartition result;
    const auto combine = [&](const ExactSection& a,const ExactSection& b,ExactSectionOperation operation) {
        work.charge(1);
        if (!work.remaining()) throw ResourceLimit("First-deposition partition budget exhausted.");
        auto combined=boolean_exact_sections(a,b,operation,work.remaining(),checkpoint);
        work.charge(combined.work);
        return std::move(combined.section);
    };
    const auto apply=[&](const ExactDepositedCommandSection& command,std::size_t lookahead) {
        if (command.envelope.empty()) return;
        if (result.regions.empty() || (result.regions.size()==1 && result.regions.front().material==command.material)) {
            // With no foreign owner, this material owns the entire prefix
            // union. One exact union establishes both quantities; sharing its
            // immutable handle proves equality without rebuilding overlays.
            // A fully occluded earlier foreign command may still have added
            // ambiguity area: preserve that evidence, not just the owner.
            auto united=combine(result.occupied_union,command.envelope,ExactSectionOperation::Union);
            work.charge(1);
            const auto overlap_vertices=result.cross_material_overlap.vertex_count();
            if (overlap_vertices+lookahead>4096 || united.vertex_count()>(4096-overlap_vertices-lookahead)/2)
                throw ResourceLimit(lookahead ? "Streamed first-deposition live vertex budget exceeded." :
                    "First-deposition retained vertex budget exceeded.");
            if (result.regions.empty()) result.regions.push_back({command.material,united});
            else result.regions.front().envelope=united;
            result.occupied_union=std::move(united);
            return;
        }
        auto found=std::find_if(result.regions.begin(),result.regions.end(),[&](const auto& region) {
            work.charge(1);
            return region.material==command.material;
        });
        const ExactSection previous=found==result.regions.end() ? ExactSection{} : found->envelope;
        const auto foreign=combine(result.occupied_union,previous,ExactSectionOperation::Difference);
        const auto competing=combine(command.envelope,foreign,ExactSectionOperation::Intersection);
        result.cross_material_overlap=combine(result.cross_material_overlap,competing,ExactSectionOperation::Union);
        const auto added=combine(command.envelope,result.occupied_union,ExactSectionOperation::Difference);
        auto owned=combine(previous,added,ExactSectionOperation::Union);
        // Existing regions remain unchanged. Only this material grows, and
        // its complete new region must remain disjoint from all prior others.
        if (!combine(owned,foreign,ExactSectionOperation::Intersection).empty())
            throw std::invalid_argument("First-deposition material partition is not disjoint.");
        if (found!=result.regions.end()) found->envelope=std::move(owned);
        else if (!owned.empty()) result.regions.push_back({command.material,std::move(owned)});
        const auto expected=combine(result.occupied_union,command.envelope,ExactSectionOperation::Union);
        ExactSection assembled;
        std::size_t retained=result.cross_material_overlap.vertex_count()+expected.vertex_count()+lookahead;
        for (const auto& region:result.regions) {
            work.charge(1);
            if (retained>4096 || region.envelope.vertex_count()>4096-retained)
                throw ResourceLimit(lookahead ? "Streamed first-deposition live vertex budget exceeded." :
                    "First-deposition retained vertex budget exceeded.");
            retained+=region.envelope.vertex_count();
            assembled=combine(assembled,region.envelope,ExactSectionOperation::Union);
        }
        if (!combine(assembled,expected,ExactSectionOperation::SymmetricDifference).empty())
            throw std::invalid_argument("First-deposition partition changes a deposited prefix union.");
        result.occupied_union=expected;
    };
    const auto live_capacity=[&] {
        std::size_t retained=result.occupied_union.vertex_count()+result.cross_material_overlap.vertex_count();
        for (const auto& region:result.regions) {
            work.charge(1);
            if (retained>4096 || region.envelope.vertex_count()>4096-retained)
                throw ResourceLimit("Streamed first-deposition live vertex budget exceeded.");
            retained+=region.envelope.vertex_count();
        }
        if (retained>4096)
            throw ResourceLimit("Streamed first-deposition live vertex budget exceeded.");
        return 4096-retained;
    };
    const auto check_live=[&](std::size_t pending) {
        if (pending>live_capacity())
            throw ResourceLimit("Streamed first-deposition live vertex budget exceeded.");
    };
    std::optional<ExactDepositedCommandSection> run;
    std::array<ExactSection,13> run_bins;
    std::size_t pending_vertices=0;
    const auto merge_run=[&](ExactSection& destination,ExactSection& source,std::size_t lookahead) {
        if (source.empty()) return;
        work.charge(1);
        if (destination.empty()) { destination=std::move(source); source={}; return; }
        auto united=combine(destination,source,ExactSectionOperation::Union);
        pending_vertices-=destination.vertex_count()+source.vertex_count();
        pending_vertices+=united.vertex_count();
        check_live(pending_vertices+lookahead);
        destination=std::move(united);
        source={};
    };
    const auto collapse_run=[&](std::size_t lookahead) {
        if (run) for (auto& bin:run_bins) merge_run(run->envelope,bin,lookahead);
    };
    for (std::size_t index=0;index<command_count;++index) {
        work.charge(1);
        if (streamed && run) {
            const auto capacity=live_capacity();
            // The provider's cut size is unknown. Consolidate BEFORE asking
            // for another (at most1024-vertex) cut, never after admitting an
            // oversized live input. Aliases in every bin count separately.
            if (capacity<1024 || pending_vertices>capacity-1024) collapse_run(0);
        }
        if (!work.remaining()) throw ResourceLimit("First-deposition provider budget exhausted.");
        auto provided=provider(index,work.remaining());
        work.charge(provided.work);
        auto& command=provided.command;
        if (command.material>=256) throw std::invalid_argument("Invalid first-deposition material ID.");
        if (!streamed) {
            apply(command,0);
            continue;
        }
        const auto incoming=command.envelope.vertex_count();
        work.charge(incoming+1);
        if (incoming>1024) throw ResourceLimit("Streamed first-deposition cut vertex budget exceeded.");
        check_live(incoming+pending_vertices);
        if (command.envelope.empty()) continue;
        if (run && run->material!=command.material) {
            // The next material's cut stays live while the preceding run is
            // allocated and certified. It cannot borrow that run's budget.
            collapse_run(incoming);
            apply(*run,incoming);
            run.reset();
            pending_vertices=0;
            check_live(incoming);
        }
        if (!run) run=ExactDepositedCommandSection{command.material,{}};
        pending_vertices+=incoming;
        // Prior foreign owners do not change within this uninterrupted run.
        // Equal-sized command groups merge as binary carries; foreign command
        // boundaries still flush the complete run before ownership is applied.
        auto carry=std::move(command.envelope);
        command.envelope={};
        for (auto& bin:run_bins) {
            if (bin.empty()) { bin=std::move(carry); carry={}; break; }
            merge_run(carry,bin,0);
        }
        if (!carry.empty()) throw ResourceLimit("Streamed first-deposition run budget exceeded.");
        check_live(pending_vertices);
    }
    if (run) { collapse_run(0); apply(*run,0); run.reset(); }
    work.charge(result.cross_material_overlap.vertex_count()+1);
    std::list<PolygonHoles> overlap;
    ExactSectionAccess::set(result.cross_material_overlap).polygons_with_holes(std::back_inserter(overlap));
    Kernel::FT area(0);
    for (const auto& polygon:overlap) {
        work.charge(polygon.outer_boundary().size()+1);
        area+=polygon.outer_boundary().area();
        for (const auto& hole:polygon.holes()) {
            work.charge(hole.size()+1);
            area+=hole.area();
        }
    }
    result.cross_material_overlap_area_m2=CGAL::to_double(area)*1e-18;
    result.work=work.used;
    return result;
}

ExactDepositedPartition exact_first_deposition_partition(
    const std::vector<ExactDepositedCommandSection>& commands,std::size_t max_work,
    const std::function<void()>& checkpoint)
{
    ExactSectionWork work(max_work,checkpoint,__func__);
    if (commands.size()>4096) throw ResourceLimit("First-deposition command budget exceeded.");
    std::size_t input_vertices=0;
    for (const auto& command:commands) {
        work.charge(1);
        if (command.material>=256) throw std::invalid_argument("Invalid first-deposition material ID.");
        if (command.envelope.vertex_count()>4096-input_vertices)
            throw ResourceLimit("First-deposition input vertex budget exceeded.");
        input_vertices+=command.envelope.vertex_count();
    }
    return first_deposition_partition_impl(commands.size(),[&](std::size_t index,std::size_t) {
        return ExactDepositedCommandResult{commands[index],0};
    },work,checkpoint,false);
}

ExactDepositedPartition exact_first_deposition_partition_streamed(std::size_t command_count,
    const std::function<ExactDepositedCommandResult(std::size_t,std::size_t)>& provider,
    std::size_t max_work,const std::function<void()>& checkpoint)
{
    ExactSectionWork work(max_work,checkpoint,__func__);
    if (!provider) throw std::invalid_argument("Missing first-deposition command provider.");
    if (command_count>4096) throw ResourceLimit("First-deposition command budget exceeded.");
    work.charge(0);
    return first_deposition_partition_impl(command_count,provider,work,checkpoint,true);
}

ExactSectionComponents exact_section_components(const ExactSection& section, std::size_t max_work,
    const std::function<void()>& checkpoint)
{
    ExactSectionWork work(max_work, checkpoint, __func__);
    const auto charge = [&](std::size_t amount) { work.charge(amount); };
    const auto& set = ExactSectionAccess::set(section);
    charge(section.vertex_count() + 1);
    std::list<PolygonHoles> polygons;
    set.polygons_with_holes(std::back_inserter(polygons));
    ExactSectionComponents result;
    std::size_t retained = 0;
    for (const auto& polygon : polygons) {
        std::size_t vertices = polygon.outer_boundary().size();
        for (const auto& hole : polygon.holes()) vertices += hole.size();
        if (vertices > 4096 - retained) throw ResourceLimit("Exact section component payload budget exceeded.");
        retained += vertices;
        charge(vertices * vertices + vertices + 1);
        result.sections.push_back(ExactSectionAccess::make(Set(polygon), charge));
    }
    result.work = work.used;
    return result;
}

ExactSectionComponents exact_section_interior_components(const ExactSection& section, std::size_t max_work,
    const std::function<void()>& checkpoint)
{
    bool callback_failed = false;
    const std::function<void()> checked_checkpoint = [&] {
        if (checkpoint) try { checkpoint(); }
        catch (...) { callback_failed = true; throw; }
    };
    ExactSectionWork work(max_work, checked_checkpoint, __func__);
    const auto charge = [&](std::size_t amount) { work.charge(amount); };
    charge(section.vertex_count() + 1);
    const auto& original = ExactSectionAccess::set(section);
    const auto validate_boundary = [&](const Set& geometry) {
        Kernel::FT twice_area(0);
        for (auto edge = geometry.arrangement().edges_begin(); edge != geometry.arrangement().edges_end(); ++edge) {
            charge(1);
            Set::Arrangement_2::Halfedge_const_handle boundary = edge;
            if (!boundary->face()->contained()) boundary = boundary->twin();
            if (!boundary->face()->contained() || boundary->twin()->face()->contained())
                throw std::invalid_argument("Invalid exact occupied interior boundary.");
            const auto& a = boundary->source()->point();
            const auto& b = boundary->target()->point();
            twice_area += a.x() * b.y() - a.y() * b.x();
        }
        if (twice_area <= 0) throw std::invalid_argument("Exact occupied interior has no positive area.");
    };
    std::size_t occupied_faces = 0;
    for (auto face = original.arrangement().faces_begin(); face != original.arrangement().faces_end(); ++face) {
        charge(1);
        if (!face->contained()) continue;
        if (face->is_unbounded()) throw std::invalid_argument("Unbounded exact occupied interior.");
        ++occupied_faces;
    }
    ExactSectionComponents result;
    if (occupied_faces == 1) {
        validate_boundary(original);
        // Most cuts already have one occupied interior, including tangent
        // holes. Preserve the immutable canonical handle without a copy or
        // an exact Boolean operation; validation still charges every edge.
        result.sections.push_back(section);
        result.work = work.used;
        return result;
    }
    bool strict_pieces = false;
    try {
        // Canonical polygon output can merge point-touching faces. Reconstruct
        // locally only when every output has strict, disjoint boundary rings and
        // the output count agrees with the independently counted occupied faces.
        std::list<PolygonHoles> polygons;
        charge(original.arrangement().number_of_edges() + 1);
        original.polygons_with_holes(std::back_inserter(polygons));
        strict_pieces = polygons.size() == occupied_faces;
        for (const auto& polygon : polygons) {
            if (!strict_pieces) break;
            std::set<Kernel::Point_2> seen;
            const auto strict_ring = [&](const Polygon2& ring, bool hole) {
                charge(ring.size() + 1);
                if (ring.size() < 3 || (hole ? ring.area() >= 0 : ring.area() <= 0)) return false;
                for (const auto& point : ring) {
                    std::size_t levels = 1;
                    for (auto count = seen.size(); count; count /= 2) ++levels;
                    charge(levels);
                    if (!seen.insert(point).second) return false;
                }
                return true;
            };
            strict_pieces = !polygon.is_unbounded() && strict_ring(polygon.outer_boundary(),false);
            for (const auto& hole : polygon.holes())
                if (strict_pieces) strict_pieces = strict_ring(hole,true);
        }
        if (strict_pieces) {
            std::size_t retained = 0;
            for (const auto& polygon : polygons) {
                std::size_t vertices = polygon.outer_boundary().size();
                for (const auto& hole : polygon.holes()) vertices += hole.size();
                if (vertices > 4096-retained) throw ResourceLimit("Exact interior component payload budget exceeded.");
                retained += vertices;
                // Set(PWH) uses a no-op default validation policy. Bound and
                // explicitly check its constructor preconditions first.
                charge(vertices*vertices + vertices + 1);
                if (!CGAL::is_valid_polygon_with_holes(polygon,Set::Traits_2{}))
                    throw std::invalid_argument("Invalid strict exact interior polygon.");
                charge(1);
                Set piece(polygon);
                std::size_t faces = 0;
                for (auto face = piece.arrangement().faces_begin(); face != piece.arrangement().faces_end(); ++face) {
                    charge(1);
                    if (!face->contained()) continue;
                    if (face->is_unbounded()) throw std::invalid_argument("Unbounded exact occupied interior.");
                    ++faces;
                }
                if (faces != 1) throw std::invalid_argument("Strict exact polygon has multiple occupied interiors.");
                validate_boundary(piece);
                result.sections.push_back(ExactSectionAccess::make(std::move(piece),charge));
            }
            result.work = work.used;
            return result;
        }
        polygons.clear();
        Set remaining = original;
        std::size_t retained = 0;
        while (!remaining.is_empty()) {
            const auto vertices = remaining.arrangement().number_of_vertices();
            // Bound the arrangement copy and normalization before invoking CGAL.
            charge(vertices * vertices + vertices + 1);
            Set selected = remaining;
            bool found = false;
            for (auto face = selected.arrangement().faces_begin(); face != selected.arrangement().faces_end(); ++face) {
                charge(1);
                if (!face->contained()) continue;
                if (face->is_unbounded()) throw std::invalid_argument("Unbounded exact occupied interior.");
                face->set_contained(!found);
                found = true;
            }
            if (!found) throw std::invalid_argument("Missing exact occupied interior face.");
            // Polygon-set output rings can join distinct faces at a vertex.
            // Keep the face itself instead of rebuilding from those weak rings.
            selected.remove_redundant_edges();
            charge(1);
            validate_boundary(selected);
            const auto count = selected.arrangement().number_of_vertices();
            if (count > 4096 - retained) throw ResourceLimit("Exact interior component payload budget exceeded.");
            retained += count;
            // Removing the selected face avoids depending on iterator order in
            // separate arrangement copies. Regularized difference adds no bonds.
            remaining = overlay_exact_sets(remaining, selected, ExactSectionOperation::Difference, charge);
            result.sections.push_back(ExactSectionAccess::make(std::move(selected), charge));
        }
        result.work = work.used;
        return result;
    } catch (const ResourceLimit& error) {
        if (callback_failed) throw;
        throw ResourceLimit(std::string(error.what()) + " Interior decomposition: original_vertices=" +
            std::to_string(section.vertex_count()) + ", occupied_faces=" + std::to_string(occupied_faces) +
            ", strict_pieces=" + (strict_pieces ? "true." : "false."));
    }
}

std::size_t validate_exact_section_topology(const ExactSection& section, std::size_t max_work,
    const std::function<void()>& checkpoint)
{
    ExactSectionWork work(max_work, checkpoint, __func__);
    work.charge(section.vertex_count()+1);
    std::list<PolygonHoles> polygons;
    ExactSectionAccess::set(section).polygons_with_holes(std::back_inserter(polygons));
    for (const auto& polygon : polygons) {
        // The source is already a canonical planar arrangement: genuine
        // crossings have been split into vertices. Repeated vertices expose
        // point-only contacts in a component's directed boundary chains.
        std::set<Kernel::Point_2> vertices;
        Kernel::FT area=0;
        const auto check = [&](const Polygon2& ring,bool hole) {
            if (ring.size()<3) throw std::invalid_argument("Invalid exact component boundary.");
            work.charge(ring.size()+1);
            const auto signed_area=ring.area();
            if (hole ? signed_area>=0 : signed_area<=0)
                throw std::invalid_argument("Invalid exact component boundary orientation.");
            area+=signed_area;
            for (const auto& point : ring) {
                std::size_t levels=1;
                for (auto count=vertices.size();count;count/=2) ++levels;
                work.charge(levels);
                if (!vertices.insert(point).second)
                    throw std::invalid_argument("Exact component boundary has a nonmanifold point contact.");
            }
        };
        check(polygon.outer_boundary(),false);
        for (const auto& hole : polygon.holes()) check(hole,true);
        if (area<=0) throw std::invalid_argument("Exact component has no positive occupied area.");
    }
    return work.used;
}

ExactSectionPointMembership exact_closed_contains_point(const ExactSection& section, const Point& point,
    std::size_t max_work, const std::function<void()>& checkpoint)
{
    if (std::abs(double(point.x()))>1e10 || std::abs(double(point.y()))>1e10)
        throw std::invalid_argument("Exact point query exceeds the coordinate limit.");
    ExactSectionWork work(max_work, checkpoint, __func__);
    work.charge(section.vertex_count()+1);
    std::list<PolygonHoles> polygons;
    ExactSectionAccess::set(section).polygons_with_holes(std::back_inserter(polygons));
    const Kernel::Point_2 query(double(point.x()),double(point.y()));
    const auto turn=Kernel().orientation_2_object();
    int winding=0;
    bool boundary=false;
    const auto ring = [&](const Polygon2& polygon) {
        for (std::size_t i=0;i<polygon.size();++i) {
            work.charge(1);
            const auto& a=polygon[i];
            const auto& b=polygon[(i+1)%polygon.size()];
            const auto side=turn(a,b,query);
            if (side==CGAL::COLLINEAR &&
                query.x()>=std::min(a.x(),b.x()) && query.x()<=std::max(a.x(),b.x()) &&
                query.y()>=std::min(a.y(),b.y()) && query.y()<=std::max(a.y(),b.y())) boundary=true;
            if (a.y()<=query.y() && b.y()>query.y() && side==CGAL::LEFT_TURN) ++winding;
            if (a.y()>query.y() && b.y()<=query.y() && side==CGAL::RIGHT_TURN) --winding;
        }
    };
    for (const auto& polygon : polygons) {
        ring(polygon.outer_boundary());
        for (const auto& hole : polygon.holes()) ring(hole);
    }
    return {boundary || winding!=0,boundary,work.used};
}

ExactSectionBilinearExtrema exact_section_bilinear_extrema(const ExactSection& section,
    const std::array<coord_t,2>& cell_start, const std::array<coord_t,2>& cell_size,
    const std::array<double,4>& corner_values, std::size_t max_work,
    const std::function<void()>& checkpoint)
{
    using FT=Kernel::FT;
    ExactSectionWork work(max_work, checkpoint, __func__);
    for (unsigned axis=0;axis<2;++axis)
        if (cell_size[axis]<=0 || double(cell_size[axis])>2e10 ||
            std::abs(double(cell_start[axis]))>1e10 ||
            std::abs(double(cell_start[axis])+double(cell_size[axis]))>1e10)
            throw std::invalid_argument("Invalid exact bilinear cell.");
    for (const double value:corner_values)
        if (!std::isfinite(value)) throw std::invalid_argument("Invalid exact bilinear corner value.");
    if (section.empty()) throw std::invalid_argument("Empty exact bilinear section.");
    if (section.vertex_count()>4096) throw ResourceLimit("Exact bilinear vertex limit exceeded.");
    work.charge(section.vertex_count()+1);
    std::list<PolygonHoles> polygons;
    ExactSectionAccess::set(section).polygons_with_holes(std::back_inserter(polygons));
    const FT a(corner_values[0]),b=FT(corner_values[1])-a,c=FT(corner_values[2])-a;
    const FT d=FT(corner_values[3])-FT(corner_values[1])-FT(corner_values[2])+a;
    FT minimum,maximum;
    std::array<FT,2> minimum_local{},maximum_local{};
    bool initialized=false;
    const auto consider = [&](const FT& x,const FT& y) {
        work.charge(1);
        if (x<0 || x>1 || y<0 || y>1)
            throw std::invalid_argument("Exact bilinear section leaves its cell.");
        const FT value=a+b*x+c*y+d*x*y;
        if (!initialized || value<minimum) { minimum=value; minimum_local={{x,y}}; }
        if (!initialized || value>maximum) { maximum=value; maximum_local={{x,y}}; }
        initialized=true;
    };
    std::size_t vertices=0;
    const auto ring = [&](const Polygon2& polygon) {
        if (polygon.size()>4096-vertices) throw ResourceLimit("Exact bilinear vertex limit exceeded.");
        vertices+=polygon.size();
        for (std::size_t i=0;i<polygon.size();++i) {
            work.charge(1);
            const auto& first=polygon[i];
            const auto& second=polygon[(i+1)%polygon.size()];
            const FT x=(first.x()-FT(double(cell_start[0])))/FT(double(cell_size[0]));
            const FT y=(first.y()-FT(double(cell_start[1])))/FT(double(cell_size[1]));
            const FT dx=(second.x()-first.x())/FT(double(cell_size[0]));
            const FT dy=(second.y()-first.y())/FT(double(cell_size[1]));
            consider(x,y);
            const FT quadratic=d*dx*dy;
            if (quadratic!=0) {
                const FT t=-(b*dx+c*dy+d*(x*dy+y*dx))/(FT(2)*quadratic);
                if (t>0 && t<1) consider(x+t*dx,y+t*dy);
            }
        }
    };
    for (const auto& polygon:polygons) {
        ring(polygon.outer_boundary());
        for (const auto& hole:polygon.holes()) ring(hole);
    }
    if (!initialized) throw std::invalid_argument("Empty exact bilinear boundary.");
    ExactSectionBilinearExtrema result;
    const auto convert = [&](const FT& value,const std::array<FT,2>& local,
        double& representative,std::array<double,2>& bounds,std::array<double,2>& position) {
        work.charge(1);
        const auto& exact_value=CGAL::exact(value);
        const auto interval=CGAL::to_interval(exact_value);
        bounds={{interval.first,interval.second}};
        representative=CGAL::to_double(exact_value);
        if (!std::isfinite(bounds[0]) || !std::isfinite(bounds[1]) || !std::isfinite(representative))
            throw std::invalid_argument("Unresolved exact bilinear value conversion.");
        representative=std::clamp(representative,bounds[0],bounds[1]);
        for (unsigned axis=0;axis<2;++axis) {
            position[axis]=CGAL::to_double(local[axis]);
            if (!std::isfinite(position[axis])) throw std::invalid_argument("Unresolved exact bilinear location.");
            position[axis]=std::clamp(position[axis],0.,1.);
        }
    };
    convert(minimum,minimum_local,result.minimum,result.minimum_bounds,result.minimum_local);
    convert(maximum,maximum_local,result.maximum,result.maximum_bounds,result.maximum_local);
    result.work=work.used;
    return result;
}

struct ExactSectionBoundary::Impl {
    std::vector<Kernel::Segment_2> segments;
    std::optional<BoundingBox> enclosure;
};
class ExactSectionBoundaryAccess {
public:
    static const std::vector<Kernel::Segment_2>& segments(const ExactSectionBoundary& boundary)
    {
        static const std::vector<Kernel::Segment_2> empty;
        return boundary.m_impl ? boundary.m_impl->segments : empty;
    }
    static ExactSectionBoundary make(std::vector<Kernel::Segment_2> segments,
        const std::function<void(std::size_t)>& charge)
    {
        if (segments.size() > 4096) throw ResourceLimit("Exact boundary retained segment budget exceeded.");
        charge(segments.size() + 1);
        if (segments.empty()) return {};
        std::optional<CGAL::Bbox_2> box;
        for (const auto& segment : segments) {
            charge(1);
            if (segment.source() == segment.target()) throw std::invalid_argument("Zero-length exact boundary segment.");
            const auto current = segment.source().bbox() + segment.target().bbox();
            box = box ? *box + current : current;
        }
        if (!std::isfinite(box->xmin()) || !std::isfinite(box->ymin()) ||
            !std::isfinite(box->xmax()) || !std::isfinite(box->ymax()) ||
            box->xmin() < -1e10 || box->ymin() < -1e10 || box->xmax() > 1e10 || box->ymax() > 1e10)
            throw std::invalid_argument("Exact boundary bounds exceed the coordinate limit.");
        auto impl = std::make_shared<ExactSectionBoundary::Impl>();
        impl->enclosure = BoundingBox(Point(coord_t(std::floor(box->xmin())),coord_t(std::floor(box->ymin()))),
            Point(coord_t(std::ceil(box->xmax())),coord_t(std::ceil(box->ymax()))));
        impl->segments = std::move(segments);
        return ExactSectionBoundary(std::move(impl));
    }
};
bool ExactSectionBoundary::empty() const { return !m_impl; }
std::size_t ExactSectionBoundary::segment_count() const { return m_impl ? m_impl->segments.size() : 0; }
std::optional<BoundingBox> ExactSectionBoundary::bounds() const { return m_impl ? m_impl->enclosure : std::nullopt; }

ExactSectionBoundaryResult extract_exact_section_boundary(const ExactSection& section, std::size_t max_work,
    const std::function<void()>& checkpoint)
{
    ExactSectionWork work(max_work, checkpoint, __func__);
    const auto charge = [&](std::size_t amount) { work.charge(amount); };
    charge(section.vertex_count() + 1);
    std::list<PolygonHoles> polygons;
    ExactSectionAccess::set(section).polygons_with_holes(std::back_inserter(polygons));
    std::vector<Kernel::Segment_2> segments;
    const auto ring = [&](const Polygon2& polygon) {
        for (std::size_t i = 0; i < polygon.size(); ++i) {
            charge(1);
            const auto& a = polygon[i];
            const auto& b = polygon[(i + 1) % polygon.size()];
            if (a == b) continue;
            if (segments.size() == 4096) throw ResourceLimit("Exact boundary extraction segment budget exceeded.");
            segments.emplace_back(a,b);
        }
    };
    for (const auto& polygon : polygons) {
        ring(polygon.outer_boundary());
        for (const auto& hole : polygon.holes()) ring(hole);
    }
    auto boundary = ExactSectionBoundaryAccess::make(std::move(segments),charge);
    return {std::move(boundary),work.used};
}

ExactSectionBoundaryResult boolean_exact_section_boundaries(const ExactSectionBoundary& left,
    const ExactSectionBoundary& right, ExactBoundaryOperation operation, ExactBoundaryOrientation orientation,
    std::size_t max_work, const std::function<void()>& checkpoint)
{
    if ((operation != ExactBoundaryOperation::Intersection && operation != ExactBoundaryOperation::Difference) ||
        (orientation != ExactBoundaryOrientation::Any && orientation != ExactBoundaryOrientation::Opposed))
        throw std::invalid_argument("Invalid exact boundary operation.");
    ExactSectionWork work(max_work, checkpoint, __func__);
    const auto charge = [&](std::size_t amount) { work.charge(amount); };
    const auto& a = ExactSectionBoundaryAccess::segments(left);
    const auto& b = ExactSectionBoundaryAccess::segments(right);
    charge(a.size() * b.size() + a.size() + b.size() + 1);
    std::vector<Kernel::Segment_2> output;
    using FT = Kernel::FT;
    using Interval = std::pair<FT,FT>;
    for (const auto& segment : a) {
        charge(1);
        const auto& start = segment.source();
        const auto direction = segment.target() - start;
        std::vector<Interval> covered;
        for (const auto& candidate : b) {
            charge(1);
            const auto turn = Kernel().orientation_2_object();
            if (turn(start,segment.target(),candidate.source()) != CGAL::COLLINEAR ||
                turn(start,segment.target(),candidate.target()) != CGAL::COLLINEAR) continue;
            const auto other = candidate.target() - candidate.source();
            if (orientation == ExactBoundaryOrientation::Opposed && direction * other >= FT(0)) continue;
            const auto parameter = [&](const Kernel::Point_2& point) {
                return direction.x() != FT(0) ? (point.x() - start.x()) / direction.x() :
                    (point.y() - start.y()) / direction.y();
            };
            const FT first = parameter(candidate.source()), second = parameter(candidate.target());
            const FT lo = std::max(FT(0),std::min(first,second));
            const FT hi = std::min(FT(1),std::max(first,second));
            if (lo < hi) covered.emplace_back(lo,hi);
        }
        std::size_t levels = 1;
        for (auto count = covered.size(); count > 1; count /= 2) ++levels;
        charge(covered.size() * levels + 1);
        std::sort(covered.begin(),covered.end());
        std::size_t merged = 0;
        for (const auto& interval : covered) {
            charge(1);
            if (merged && interval.first <= covered[merged-1].second)
                covered[merged-1].second = std::max(covered[merged-1].second,interval.second);
            else covered[merged++] = interval;
        }
        covered.resize(merged);
        const auto emit = [&](const FT& lo,const FT& hi) {
            charge(1);
            if (lo >= hi) return;
            if (output.size() == 4096) throw ResourceLimit("Exact boundary result segment budget exceeded.");
            output.emplace_back(start + direction * lo,start + direction * hi);
        };
        FT cursor = 0;
        for (const auto& interval : covered) {
            if (operation == ExactBoundaryOperation::Intersection) emit(interval.first,interval.second);
            else emit(cursor,interval.first);
            cursor = interval.second;
        }
        if (operation == ExactBoundaryOperation::Difference) emit(cursor,FT(1));
    }
    auto boundary = ExactSectionBoundaryAccess::make(std::move(output),charge);
    return {std::move(boundary),work.used};
}

ExactSectionBoundaryResult exact_boundary_on_native_plane(const ExactSectionBoundary& boundary,
    unsigned axis, coord_t coordinate, std::size_t max_work, const std::function<void()>& checkpoint)
{
    if (axis > 1 || std::abs(double(coordinate)) > 1e10)
        throw std::invalid_argument("Invalid exact boundary native plane.");
    ExactSectionWork work(max_work, checkpoint, __func__);
    const auto charge = [&](std::size_t amount) { work.charge(amount); };
    const Kernel::FT plane{double(coordinate)};
    std::vector<Kernel::Segment_2> selected;
    for (const auto& segment : ExactSectionBoundaryAccess::segments(boundary)) {
        charge(1);
        const auto a = axis == 0 ? segment.source().x() : segment.source().y();
        const auto b = axis == 0 ? segment.target().x() : segment.target().y();
        if (a == plane && b == plane) selected.push_back(segment);
    }
    auto result = ExactSectionBoundaryAccess::make(std::move(selected),charge);
    return {std::move(result),work.used};
}

ExactSectionPredicate exact_boundary_overlaps_native_interval(const ExactSectionBoundary& boundary,
    unsigned axis, coord_t plane, coord_t lo, coord_t hi, std::size_t max_work,
    const std::function<void()>& checkpoint)
{
    if (axis > 1 || lo > hi || std::abs(double(plane)) > 1e10 ||
        std::abs(double(lo)) > 1e10 || std::abs(double(hi)) > 1e10)
        throw std::invalid_argument("Invalid exact boundary native interval.");
    ExactSectionWork work(max_work, checkpoint, __func__);
    const Kernel::FT fixed{double(plane)}, lower{double(lo)}, upper{double(hi)};
    work.charge(1);
    if (lo == hi) return {true,work.used};
    for (const auto& segment : ExactSectionBoundaryAccess::segments(boundary)) {
        work.charge(1);
        const auto fixed_a = axis == 0 ? segment.source().x() : segment.source().y();
        const auto fixed_b = axis == 0 ? segment.target().x() : segment.target().y();
        if (fixed_a != fixed || fixed_b != fixed) continue;
        const auto a = axis == 0 ? segment.source().y() : segment.source().x();
        const auto b = axis == 0 ? segment.target().y() : segment.target().x();
        if (std::max(lower,std::min(a,b)) < std::min(upper,std::max(a,b)))
            return {false,work.used};
    }
    return {true,work.used};
}

ExactBoundarySegments split_exact_boundary_segments(const ExactSectionBoundary& boundary,
    std::size_t max_work, const std::function<void()>& checkpoint)
{
    ExactSectionWork work(max_work, checkpoint, __func__);
    const auto charge = [&](std::size_t amount) { work.charge(amount); };
    const auto& source=ExactSectionBoundaryAccess::segments(boundary);
    if (source.size()>4096) throw ResourceLimit("Exact boundary split segment budget exceeded.");
    charge(source.size()+1);
    ExactBoundarySegments result;
    result.segments.reserve(source.size());
    for (const auto& segment : source) {
        charge(1);
        const auto delta=segment.target()-segment.source();
        const auto length_squared=delta.squared_length();
        if (length_squared<=0) throw std::invalid_argument("Exact boundary normal has no positive length.");
        const double length=std::sqrt(CGAL::to_double(length_squared));
        if (!(length>0) || !std::isfinite(length))
            throw std::invalid_argument("Unresolved exact boundary normal length.");
        const std::array<Kernel::FT,2> numerators{{delta.y(),-delta.x()}};
        ExactBoundarySegment record;
        for (std::size_t axis=0;axis<2;++axis) {
            charge(1);
            const auto& numerator=numerators[axis];
            const auto square=numerator*numerator/length_squared;
            const double estimate=CGAL::to_double(numerator)/length;
            const auto bounds=certified_nonnegative_root(square,std::abs(estimate),charge);
            record.outward_normal_bounds[axis]=numerator<0 ?
                std::array<double,2>{{-bounds[1],-bounds[0]}} : bounds;
            record.outward_normal[axis]=std::clamp(estimate,
                record.outward_normal_bounds[axis][0],record.outward_normal_bounds[axis][1]);
        }
        const double norm_squared=record.outward_normal[0]*record.outward_normal[0]+
            record.outward_normal[1]*record.outward_normal[1];
        if (std::abs(norm_squared-1.)>32*std::numeric_limits<double>::epsilon()) {
            std::ostringstream message;
            message.precision(17);
            message << "Exact boundary normal estimate is not unit length: normal=(" << record.outward_normal[0]
                << ',' << record.outward_normal[1] << "), bounds=[" << record.outward_normal_bounds[0][0]
                << ',' << record.outward_normal_bounds[0][1] << "]x[" << record.outward_normal_bounds[1][0]
                << ',' << record.outward_normal_bounds[1][1] << "], squared_norm_error=" << norm_squared-1. << '.';
            throw std::invalid_argument(message.str());
        }
        charge(1); // Reserve the single-span payload before copying endpoints.
        record.span=ExactSectionBoundaryAccess::make({segment},charge);
        result.segments.push_back(std::move(record));
    }
    result.work=work.used;
    return result;
}

static ExactSegmentCoverage section_covers_segment(const ExPolygons& polygons, const Point& start, const Point& end,
    std::size_t max_work, const std::function<void()>& checkpoint, unsigned boundary_allowance_nm, bool native_nonzero_fill)
{
    ExactSegmentCoverage result;
    if (!max_work || max_work > 4000000 || boundary_allowance_nm > 7)
        throw std::invalid_argument("Invalid exact segment coverage budget or boundary allowance.");
    for (const auto& p : {start,end})
        if (std::abs(double(p.x())) > 1e10 || std::abs(double(p.y())) > 1e10)
            throw std::invalid_argument("Exact segment coordinate outside supported range.");
    const auto charge = [&](std::size_t n) {
        if (checkpoint) checkpoint();
        if (n > max_work-result.work) throw ResourceLimit("Exact segment coverage work limit exceeded.");
        result.work += n;
    };
    const coord_t allowance = boundary_allowance_nm;
    const BoundingBox query(Point(std::min(start.x(),end.x())-allowance,std::min(start.y(),end.y())-allowance),
                            Point(std::max(start.x(),end.x())+allowance,std::max(start.y(),end.y())+allowance));
    Set occupied;
    if (native_nonzero_fill) {
        // The query and its boundary capsules are strictly inside this exact
        // crop. Artificial closure edges therefore cannot prove band coverage.
        // Validate/charge the whole input before allocating clipped geometry;
        // retained chains and arrangements keep their existing independent caps.
        const auto validate = [&](const Polygon& ring) {
            charge(ring.points.size()+1);
            for (const auto& p : ring.points)
                if (std::abs(double(p.x())) > 1e10 || std::abs(double(p.y())) > 1e10)
                    throw std::invalid_argument("Exact segment polygon coordinate outside supported range.");
        };
        for (const auto& polygon : polygons) {
            validate(polygon.contour);
            for (const auto& hole : polygon.holes) validate(hole);
        }
        if (result.work == max_work) throw ResourceLimit("Exact segment coverage work limit exceeded.");
        const std::array<Kernel::FT,4> crop{Kernel::FT(double(query.min.x()-2)),Kernel::FT(double(query.max.x()+2)),
            Kernel::FT(double(query.min.y()-2)),Kernel::FT(double(query.max.y()+2))};
        std::size_t construction=0;
        occupied=make_clipped_native_union(polygons,crop,max_work-result.work,checkpoint,construction);
        charge(construction);
    } else {
        std::size_t scanned = 0;
        const auto overlaps = [&](const Polygon& ring) {
            if (ring.points.size() > 2048-scanned) throw ResourceLimit("Exact segment input vertex limit exceeded.");
            scanned += ring.points.size();
            charge(ring.points.size()+1);
            for (const auto& p : ring.points)
                if (std::abs(double(p.x())) > 1e10 || std::abs(double(p.y())) > 1e10)
                    throw std::invalid_argument("Exact segment polygon coordinate outside supported range.");
            return !ring.points.empty() && query.overlap(BoundingBox(ring.points));
        };
        // A polygon or hole with a disjoint bounding box cannot affect any point
        // on this segment. Restrict the query without clipping/rounding boundaries;
        // this certifies segment coverage, not global validity of unrelated rings.
        ExPolygons relevant;
        for (const auto& polygon : polygons) {
            const bool outer = overlaps(polygon.contour);
            ExPolygon local;
            if (outer) local.contour = polygon.contour;
            for (const auto& hole : polygon.holes)
                if (overlaps(hole) && outer) local.holes.push_back(hole);
            if (!local.contour.empty()) relevant.push_back(std::move(local));
        }
        if (result.work == max_work) throw ResourceLimit("Exact segment coverage work limit exceeded.");
        std::size_t construction = 0;
        occupied = make_topology_union(relevant,max_work-result.work,checkpoint,construction,false);
        charge(construction);
    }
    const Kernel::Point_2 a(double(start.x()),double(start.y())), b(double(end.x()),double(end.y()));
    const auto direction = b-a;
    const auto edges = occupied.arrangement().number_of_edges();
    const auto covered = [&](const Kernel::Point_2& p) {
        charge(edges+1);
        return occupied.oriented_side(p) != CGAL::ON_NEGATIVE_SIDE;
    };
    const auto band_contains = [&](const Kernel::Point_2& first,const Kernel::Point_2& last) {
        if (!boundary_allowance_nm) return false;
        const Kernel::FT radius2(boundary_allowance_nm*boundary_allowance_nm);
        for (auto edge = occupied.arrangement().edges_begin(); edge != occupied.arrangement().edges_end(); ++edge) {
            charge(2);
            const auto p = edge->source()->point(), q = edge->target()->point();
            // A segment's radius-r capsule is convex. Both exact endpoints in
            // the SAME capsule certify the whole interval, not just samples.
            if (segment_distance_squared(first,p,q) <= radius2 && segment_distance_squared(last,p,q) <= radius2)
                return true;
        }
        return false;
    };
    result.covered = true;
    result.within_boundary_band = true;
    for (const auto& p : {a,b}) {
        if (covered(p)) continue;
        result.covered = false;
        if (!band_contains(p,p)) { result.within_boundary_band = false; return result; }
    }
    if (a == b) return result;
    const auto cross = [](const Kernel::Vector_2& u,const Kernel::Vector_2& v) { return u.x()*v.y()-u.y()*v.x(); };
    std::vector<Kernel::FT> cuts{0,1};
    for (auto edge = occupied.arrangement().edges_begin(); edge != occupied.arrangement().edges_end(); ++edge) {
        charge(1);
        const auto c = edge->source()->point(), d = edge->target()->point();
        const auto side = d-c, offset = c-a;
        const auto denominator = cross(direction,side);
        if (denominator != 0) {
            const auto t = cross(offset,side)/denominator, u = cross(offset,direction)/denominator;
            if (t > 0 && t < 1 && u >= 0 && u <= 1) cuts.push_back(t);
        } else if (cross(offset,direction) == 0) {
            for (const auto& p : {c,d}) {
                const auto t = direction.x() != 0 ? (p.x()-a.x())/direction.x() : (p.y()-a.y())/direction.y();
                if (t > 0 && t < 1) cuts.push_back(t);
            }
        }
    }
    std::size_t levels = 1;
    for (auto n = cuts.size(); n > 1; n /= 2) ++levels;
    charge(cuts.size()*levels);
    std::sort(cuts.begin(),cuts.end());
    cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    for (std::size_t i = 1; i < cuts.size(); ++i) {
        if (covered(a+direction*((cuts[i-1]+cuts[i])/2))) continue;
        result.covered = false;
        if (!band_contains(a+direction*cuts[i-1],a+direction*cuts[i])) {
            result.within_boundary_band = false;
            return result;
        }
    }
    return result;
}
ExactSegmentCoverage exact_section_covers_segment(const ExPolygons& polygons, const Point& start, const Point& end,
    std::size_t max_work, const std::function<void()>& checkpoint, unsigned boundary_allowance_nm)
{
    return section_covers_segment(polygons,start,end,max_work,checkpoint,boundary_allowance_nm,false);
}
ExactSegmentCoverage exact_native_section_covers_segment(const ExPolygons& polygons, const Point& start, const Point& end,
    std::size_t max_work, const std::function<void()>& checkpoint, unsigned boundary_allowance_nm)
{
    return section_covers_segment(polygons,start,end,max_work,checkpoint,boundary_allowance_nm,true);
}
ExactSectionTopologySummary diagnose_exact_section_topology(const ExPolygons& polygons,
    std::size_t max_work, const std::function<void()>& checkpoint)
{
    ExactSectionTopologySummary result;
    const auto occupied = make_topology_union(polygons, max_work, checkpoint, result.work);
    std::list<PolygonHoles> components;
    occupied.polygons_with_holes(std::back_inserter(components));
    result.components = components.size();
    for (const auto& component : components) {
        if (checkpoint) checkpoint();
        result.holes += component.number_of_holes();
    }
    return result;
}
template<class Charge>
static std::vector<std::size_t> exact_component_correspondence(const Set& left, const Set& right, const Charge& charge)
{
    const auto pieces = [&](const Set& set) {
        charge(set.arrangement().number_of_vertices() + 1);
        std::list<PolygonHoles> list;
        set.polygons_with_holes(std::back_inserter(list));
        return std::vector<PolygonHoles>(list.begin(), list.end());
    };
    const auto a = pieces(left), b = pieces(right);
    const auto match = [&](const std::vector<PolygonHoles>& first, const std::vector<PolygonHoles>& second, const std::string& kind) {
        if (first.size() != second.size()) throw std::invalid_argument("Exact topology changes " + kind + ": " +
            std::to_string(first.size()) + " -> " + std::to_string(second.size()) + ".");
        const auto prepare = [&](const std::vector<PolygonHoles>& polygons) {
            std::vector<Set> sets;
            sets.reserve(polygons.size());
            for (const auto& polygon : polygons) {
                std::size_t count = polygon.outer_boundary().size();
                for (auto hole = polygon.holes_begin(); hole != polygon.holes_end(); ++hole) count += hole->size();
                charge(count * count + count + 1);
                sets.emplace_back(polygon);
            }
            return sets;
        };
        auto first_sets = prepare(first);
        const auto second_sets = prepare(second);
        std::vector<std::size_t> matches;
        std::vector<bool> used(second.size(), false);
        for (std::size_t source = 0; source < first.size(); ++source) {
            const auto& p = first[source];
            std::size_t found = second.size();
            for (std::size_t i = 0; i < second.size(); ++i) {
                charge(1);
                const auto x = p.outer_boundary().bbox(), y = second[i].outer_boundary().bbox();
                if (std::min(x.xmax(), y.xmax()) <= std::max(x.xmin(), y.xmin()) ||
                    std::min(x.ymax(), y.ymax()) <= std::max(x.ymin(), y.ymin())) continue;
                auto& left_set = first_sets[source];
                const auto& right_set = second_sets[i];
                const std::size_t n = left_set.arrangement().number_of_edges(), m = right_set.arrangement().number_of_edges();
                charge(2 * n * m + n + m + left_set.arrangement().number_of_vertices() +
                    right_set.arrangement().number_of_vertices() + 1);
                const bool overlaps = left_set.do_intersect(right_set);
                charge(0);
                if (!overlaps) continue;
                if (found != second.size() || used[i]) throw std::invalid_argument("Exact topology correspondence is ambiguous.");
                found = i;
            }
            if (found == second.size()) throw std::invalid_argument("Exact topology loses a component or hole.");
            used[found] = true;
            matches.push_back(found);
        }
        return matches;
    };
    const auto correspondence = match(a, b, "components");
    const auto holes = [&](const PolygonHoles& polygon) {
        std::vector<PolygonHoles> result;
        for (auto hole = polygon.holes_begin(); hole != polygon.holes_end(); ++hole) {
            charge(hole->size() + 1);
            auto ring = *hole;
            ring.reverse_orientation();
            result.emplace_back(ring);
        }
        return result;
    };
    for (std::size_t i = 0; i < a.size(); ++i)
        match(holes(a[i]), holes(b[correspondence[i]]), "holes in component " + std::to_string(i));

    return correspondence;
}
std::size_t qualify_exact_section_topology(const ExPolygons& before, const ExPolygons& after,
    std::size_t max_work, const std::function<void()>& checkpoint)
{
    if (!max_work || max_work > 20000000) throw std::invalid_argument("Invalid exact topology certificate budget.");
    std::size_t work = 0;
    const auto charge = [&](std::size_t amount) {
        if (checkpoint) checkpoint();
        if (amount > max_work - work) throw ResourceLimit("Exact topology certificate work limit exceeded.");
        work += amount;
    };
    const auto make_union = [&](const ExPolygons& input) {
        if (work == max_work) throw ResourceLimit("Exact topology certificate work limit exceeded.");
        std::size_t consumed = 0;
        auto result = make_topology_union(input, std::min<std::size_t>(4000000, max_work - work), checkpoint, consumed);
        charge(consumed);
        return result;
    };
    const auto left = make_union(before), right = make_union(after);
    exact_component_correspondence(left,right,charge);
    return work;
}
namespace {
void validate_exact_moment_cell(const std::array<coord_t,3>& start, const std::array<coord_t,3>& size)
{
    for (unsigned axis = 0; axis < 3; ++axis)
        if (size[axis] <= 0 || double(size[axis]) > 2e10 || std::abs(double(start[axis])) > 1e10 ||
            std::abs(double(start[axis]) + double(size[axis])) > 1e10)
            throw std::invalid_argument("Invalid exact section moment cell.");
}
std::array<Kernel::FT, 9> exact_xy_moments(const Set& occupied, const std::array<coord_t, 3>& cell_start,
    const std::array<coord_t, 3>& cell_size, const std::function<void(std::size_t)>& charge)
{
    using FT = Kernel::FT;
    const coord_t x0 = cell_start[0], y0 = cell_start[1];
    constexpr unsigned binomial[4][4]{{1,0,0,0},{1,1,0,0},{1,2,1,0},{1,3,3,1}};
    std::array<FT, 9> xy{};
    const auto integrate_ring = [&](const Polygon2& ring) {
        for (std::size_t edge = 0; edge < ring.size(); ++edge) {
            const auto& first = ring[edge];
            const auto& second = ring[(edge + 1) % ring.size()];
            const FT x = (first.x() - FT(double(x0))) / FT(double(cell_size[0]));
            const FT y = (first.y() - FT(double(y0))) / FT(double(cell_size[1]));
            const FT dx = (second.x() - first.x()) / FT(double(cell_size[0]));
            const FT dy = (second.y() - first.y()) / FT(double(cell_size[1]));
            // Reuse the exact power DAG across all nine XY monomials rather
            // than rebuilding lazy products in every expansion term. Existing
            // per-term work/checkpoints remain a conservative bounded charge.
            std::array<FT,4> xp{FT(1),x}, dxp{FT(1),dx};
            for (unsigned i = 2; i < 4; ++i) {
                xp[i] = xp[i-1] * x;
                dxp[i] = dxp[i-1] * dx;
            }
            const std::array<FT,3> yp{FT(1),y,y*y}, dyp{FT(1),dy,dy*dy};
            // Green's theorem: integral x^p y^q dA =
            // integral_boundary x^(p+1) y^q dy / (p+1). Expand each
            // linear edge exactly; directed hole rings subtract naturally.
            for (unsigned q = 0; q < 3; ++q)
                for (unsigned p = 0; p < 3; ++p)
                    for (unsigned i = 0; i <= p + 1; ++i)
                        for (unsigned j = 0; j <= q; ++j) {
                            charge(1);
                            xy[p + 3*q] += FT(int(binomial[p+1][i] * binomial[q][j])) *
                                xp[p+1-i] * dxp[i] * yp[q-j] * dyp[j] * dy /
                                FT(int((p+1) * (i+j+1)));
                        }
        }
    };
    charge(occupied.arrangement().number_of_vertices() + 1);
    std::list<PolygonHoles> polygons;
    occupied.polygons_with_holes(std::back_inserter(polygons));
    for (const auto& polygon : polygons) {
        integrate_ring(polygon.outer_boundary());
        for (const auto& hole : polygon.holes()) integrate_ring(hole);
    }
    return xy;
}
// Cell coordinates are validated by the caller. The returned reference borrows
// either the immutable input or caller-owned scratch, never a local set copy.
const Set& exact_cell_occupancy(const ExactSection& section,
    const std::array<coord_t,3>& cell_start, const std::array<coord_t,3>& cell_size,
    std::optional<Set>& scratch, const std::function<std::size_t()>& remaining,
    const std::function<void(std::size_t)>& charge, const std::function<void()>& checkpoint)
{
    const auto bounds = section.bounds();
    const auto x1 = cell_start[0] + cell_size[0], y1 = cell_start[1] + cell_size[1];
    // The outward integer enclosure contains every rational vertex, so its
    // containment in the convex cell rectangle proves containment of all fill.
    const bool enclosed = bounds && bounds->min.x() >= cell_start[0] && bounds->min.y() >= cell_start[1] &&
        bounds->max.x() <= x1 && bounds->max.y() <= y1;
    const auto& original = ExactSectionAccess::set(section);
    if (section.empty() || enclosed) return original;
    const auto available = remaining();
    if (!available) throw ResourceLimit("Exact section cell clipping budget exhausted.");
    const auto cell = make_exact_section({ExPolygon(Polygon{{cell_start[0],cell_start[1]},
        {x1,cell_start[1]},{x1,y1},{cell_start[0],y1}})}, available, checkpoint);
    charge(cell.work);
    scratch = overlay_exact_sets(original, ExactSectionAccess::set(cell.section), ExactSectionOperation::Intersection, charge);
    return *scratch;
}

ExactSectionMeasure integrate_exact_measure(const ExactSection& section,
    const std::array<coord_t, 3>& cell_start, const std::array<coord_t, 3>& cell_size,
    coord_t bottom, coord_t top, bool horizontal, std::size_t max_work,
    const std::function<void()>& checkpoint)
{
    if ((horizontal ? bottom != top : bottom >= top) ||
        std::abs(double(bottom)) > 1e10 || std::abs(double(top)) > 1e10)
        throw std::invalid_argument("Invalid exact section measure limits.");
    validate_exact_moment_cell(cell_start, cell_size);
    ExactSectionWork work(max_work, checkpoint, horizontal ? "exact_section_horizontal_area" : "exact_section_volume");
    const auto charge = [&](std::size_t amount) { work.charge(amount); };
    charge(section.vertex_count() + 1);
    ExactSectionMeasure result;
    const coord_t lo = std::max(bottom, cell_start[2]);
    const coord_t hi = std::min(top, cell_start[2] + cell_size[2]);
    const bool outside = horizontal ? bottom < cell_start[2] || bottom > cell_start[2] + cell_size[2] : lo >= hi;
    if (outside || section.empty()) {
        result.work = work.used;
        return result;
    }
    std::optional<Set> clipped;
    const auto& occupied = exact_cell_occupancy(section, cell_start, cell_size, clipped,
        [&] { return work.remaining(); }, charge, checkpoint);
    Kernel::FT twice_area(0);
    for (auto edge = occupied.arrangement().edges_begin(); edge != occupied.arrangement().edges_end(); ++edge) {
        charge(1);
        Set::Arrangement_2::Halfedge_const_handle boundary = edge;
        if (!boundary->face()->contained()) boundary = boundary->twin();
        if (!boundary->face()->contained() || boundary->twin()->face()->contained())
            throw std::invalid_argument("Invalid exact section measure boundary.");
        const auto& a = boundary->source()->point();
        const auto& b = boundary->target()->point();
        twice_area += a.x() * b.y() - a.y() * b.x();
    }
    using FT = Kernel::FT;
    const FT metre(1000000000);
    const FT measure = horizontal ? twice_area / (FT(2) * metre * metre) :
        twice_area * FT(double(hi - lo)) / (FT(2) * metre * metre * metre);
    charge(1);
    if (measure < 0) throw std::invalid_argument("Negative exact occupied section measure.");
    const auto interval = CGAL::to_interval(measure);
    result.value = CGAL::to_double(measure);
    if (!std::isfinite(interval.first) || !std::isfinite(interval.second) || !std::isfinite(result.value) ||
        result.value < interval.first || result.value > interval.second)
        throw std::invalid_argument("Unresolved exact occupied section measure conversion.");
    result.bounds = {interval.first, interval.second};
    result.work = work.used;
    return result;
}

ExactSectionMoments integrate_exact_section(const ExactSection& section, const std::array<coord_t, 3>& cell_start,
    const std::array<coord_t, 3>& cell_size, coord_t bottom, coord_t top,
    bool horizontal, std::size_t max_work, const std::function<void()>& checkpoint)
{
    if (!max_work || max_work > 4000000 || (horizontal ? bottom != top : bottom >= top) ||
        std::abs(double(bottom)) > 1e10 || std::abs(double(top)) > 1e10)
        throw std::invalid_argument("Invalid exact section moment limits.");
    validate_exact_moment_cell(cell_start,cell_size);
    ExactSectionMoments result;
    const auto charge = [&](std::size_t amount) {
        if (checkpoint) checkpoint();
        if (amount > max_work - result.work)
            throw ResourceLimit(std::string("Exact section moment work budget exceeded: operation=") +
                (horizontal ? "exact_section_horizontal_surface_moments" : "exact_section_volume_moments") +
                ", used=" + std::to_string(result.work) + ", requested=" + std::to_string(amount) +
                ", limit=" + std::to_string(max_work) + ", per_operation_cap=4000000.");
        result.work += amount;
    };
    charge(section.vertex_count() + 1);
    std::optional<Set> clipped;
    const auto& occupied = exact_cell_occupancy(section, cell_start, cell_size, clipped,
        [&] { return max_work - result.work; }, charge, checkpoint);
    const coord_t lo = std::max(bottom, cell_start[2]);
    const coord_t hi = std::min(top, cell_start[2] + cell_size[2]);
    const bool outside = horizontal ? bottom < cell_start[2] || bottom > cell_start[2] + cell_size[2] : lo >= hi;
    if (outside || occupied.is_empty()) return result;

    using FT = Kernel::FT;
    const auto xy = exact_xy_moments(occupied, cell_start, cell_size, charge);
    const FT z0 = FT(double(lo - cell_start[2])) / FT(double(cell_size[2]));
    const FT z1 = FT(double(hi - cell_start[2])) / FT(double(cell_size[2]));
    std::array<FT,4> z0p{FT(1),z0}, z1p{FT(1),z1};
    for (unsigned i = 2; i < 4; ++i) {
        z0p[i] = z0p[i-1] * z0;
        z1p[i] = z1p[i-1] * z1;
    }
    std::array<FT,3> z_moments;
    for (unsigned r = 0; r < 3; ++r)
        z_moments[r] = horizontal ? z0p[r] : (z1p[r+1]-z0p[r+1])/FT(int(r+1));
    const FT area_scale = FT(double(cell_size[0])) * FT(double(cell_size[1])) /
        (FT(1000000000) * FT(1000000000));
    const FT scale = horizontal ? area_scale : area_scale * FT(double(cell_size[2])) / FT(1000000000);
    for (unsigned r = 0; r < 3; ++r)
        for (unsigned q = 0; q < 3; ++q)
            for (unsigned p = 0; p < 3; ++p) {
                charge(1);
                const FT moment = xy[p + 3*q] * z_moments[r] * scale;
                if (moment < 0) throw std::invalid_argument("Negative exact occupied section moment.");
                const auto interval = CGAL::to_interval(moment);
                const double value = CGAL::to_double(moment);
                if (!std::isfinite(interval.first) || !std::isfinite(interval.second) || !std::isfinite(value) ||
                    value < interval.first || value > interval.second)
                    throw std::invalid_argument("Unresolved exact occupied section moment conversion.");
                const auto index = p + 3*q + 9*r;
                result.values[index] = value;
                result.bounds[index] = {{interval.first, interval.second}};
            }
    return result;
}
} // namespace

ExactSectionMeasure exact_section_volume(const ExactSection& section,
    const std::array<coord_t, 3>& cell_start, const std::array<coord_t, 3>& cell_size,
    coord_t bottom, coord_t top, std::size_t max_work, const std::function<void()>& checkpoint)
{
    return integrate_exact_measure(section, cell_start, cell_size, bottom, top, false, max_work, checkpoint);
}
ExactSectionMeasure exact_section_horizontal_area(const ExactSection& section,
    const std::array<coord_t, 3>& cell_start, const std::array<coord_t, 3>& cell_size,
    coord_t z, std::size_t max_work, const std::function<void()>& checkpoint)
{
    return integrate_exact_measure(section, cell_start, cell_size, z, z, true, max_work, checkpoint);
}

ExactSectionMoments exact_section_volume_moments(const ExactSection& section,
    const std::array<coord_t, 3>& cell_start, const std::array<coord_t, 3>& cell_size,
    coord_t bottom, coord_t top, std::size_t max_work, const std::function<void()>& checkpoint)
{
    return integrate_exact_section(section, cell_start, cell_size, bottom, top, false, max_work, checkpoint);
}
ExactSectionMoments exact_section_horizontal_surface_moments(const ExactSection& section,
    const std::array<coord_t, 3>& cell_start, const std::array<coord_t, 3>& cell_size,
    coord_t z, std::size_t max_work, const std::function<void()>& checkpoint)
{
    return integrate_exact_section(section, cell_start, cell_size, z, z, true, max_work, checkpoint);
}

ExactSectionMoments exact_section_vertical_surface_moments(const ExactSectionBoundary& boundary,
    const std::array<coord_t,3>& cell_start, const std::array<coord_t,3>& cell_size,
    coord_t bottom, coord_t top, ExactBoundaryCellOwnership ownership,
    std::size_t max_work, const std::function<void()>& checkpoint)
{
    validate_exact_moment_cell(cell_start,cell_size);
    if (bottom > top || std::abs(double(bottom)) > 1e10 || std::abs(double(top)) > 1e10 ||
        (ownership != ExactBoundaryCellOwnership::MaterialSide && ownership != ExactBoundaryCellOwnership::AllIncident))
        throw std::invalid_argument("Invalid exact vertical surface request.");
    ExactSectionWork work(max_work, checkpoint, __func__);
    const auto charge = [&](std::size_t amount) { work.charge(amount); };
    charge(boundary.segment_count()+1);
    ExactSectionMoments result;
    const coord_t bottom_clipped=std::max(bottom,cell_start[2]);
    const coord_t top_clipped=std::min(top,cell_start[2]+cell_size[2]);
    if (bottom_clipped>=top_clipped) { result.work=work.used; return result; }
    using FT=Kernel::FT;
    const auto power = [](FT value,unsigned exponent) {
        FT product(1);
        for (unsigned i=0;i<exponent;++i) product*=value;
        return product;
    };
    const FT z0=FT(double(bottom_clipped-cell_start[2]))/FT(double(cell_size[2]));
    const FT z1=FT(double(top_clipped-cell_start[2]))/FT(double(cell_size[2]));
    std::array<FT,3> z_moments;
    for (unsigned r=0;r<3;++r)
        z_moments[r]=(power(z1,r+1)-power(z0,r+1))/FT(int(r+1))*FT(double(cell_size[2]))/
            (FT(1000000000)*FT(1000000000));
    std::array<FT,27> lower{},upper{},representative{};
    constexpr unsigned choose[3][3]{{1,0,0},{1,1,0},{1,2,1}};
    for (const auto& segment : ExactSectionBoundaryAccess::segments(boundary)) {
        charge(1);
        const auto& origin=segment.source();
        const auto direction=segment.target()-origin;
        FT lo=0,hi=1;
        bool outside=false;
        for (unsigned axis=0;axis<2;++axis) {
            charge(1);
            const FT position=axis==0 ? origin.x() : origin.y();
            const FT delta=axis==0 ? direction.x() : direction.y();
            const FT minimum{double(cell_start[axis])},maximum{double(cell_start[axis]+cell_size[axis])};
            if (delta==0) {
                if (position<minimum || position>maximum) outside=true;
            } else {
                const FT first=(minimum-position)/delta,second=(maximum-position)/delta;
                lo=std::max(lo,std::min(first,second));
                hi=std::min(hi,std::max(first,second));
            }
        }
        if (outside || lo>=hi) continue;
        const auto a=origin+direction*lo,b=origin+direction*hi;
        const auto delta=b-a;
        if (ownership==ExactBoundaryCellOwnership::MaterialSide) {
            const FT xmin{double(cell_start[0])},xmax{double(cell_start[0]+cell_size[0])};
            const FT ymin{double(cell_start[1])},ymax{double(cell_start[1]+cell_size[1])};
            // Directed boundary material is LEFT: reject an incident cell
            // lying entirely on the other side of a cell-aligned segment.
            if ((a.x()==xmin && b.x()==xmin && delta.y()>=0) ||
                (a.x()==xmax && b.x()==xmax && delta.y()<=0) ||
                (a.y()==ymin && b.y()==ymin && delta.x()<=0) ||
                (a.y()==ymax && b.y()==ymax && delta.x()>=0)) continue;
        }
        const FT x=(a.x()-FT(double(cell_start[0])))/FT(double(cell_size[0]));
        const FT y=(a.y()-FT(double(cell_start[1])))/FT(double(cell_size[1]));
        const FT dx=delta.x()/FT(double(cell_size[0])),dy=delta.y()/FT(double(cell_size[1]));
        const FT length_squared=delta.squared_length();
        const double length_estimate=std::sqrt(CGAL::to_double(length_squared));
        const auto length_bounds=certified_nonnegative_root(length_squared,length_estimate,charge);
        const FT length_lower(length_bounds[0]),length_upper(length_bounds[1]);
        const FT length_representative(std::clamp(length_estimate,length_bounds[0],length_bounds[1]));
        std::array<FT,9> xy{};
        for (unsigned q=0;q<3;++q)
            for (unsigned p=0;p<3;++p)
                for (unsigned i=0;i<=p;++i)
                    for (unsigned j=0;j<=q;++j) {
                        charge(1);
                        xy[p+3*q]+=FT(int(choose[p][i]*choose[q][j]))*power(x,p-i)*power(dx,i)*
                            power(y,q-j)*power(dy,j)/FT(int(i+j+1));
                    }
        for (unsigned r=0;r<3;++r)
            for (unsigned q=0;q<3;++q)
                for (unsigned p=0;p<3;++p) {
                    charge(1);
                    const auto index=p+3*q+9*r;
                    const FT coefficient=xy[p+3*q]*z_moments[r];
                    if (coefficient<0) throw std::invalid_argument("Negative exact vertical surface moment.");
                    if (coefficient==0) continue;
                    // All moments share this segment's length. Nonnegative
                    // rational scaling preserves its certified bounds without
                    // forming and recovering 27 separate square roots.
                    lower[index]+=coefficient*length_lower;
                    upper[index]+=coefficient*length_upper;
                    representative[index]+=coefficient*length_representative;
                }
    }
    for (std::size_t i=0;i<27;++i) {
        charge(1);
        const double lo=CGAL::to_interval(lower[i]).first,hi=CGAL::to_interval(upper[i]).second;
        const double value=CGAL::to_double(representative[i]);
        if (!std::isfinite(lo) || !std::isfinite(hi) || !std::isfinite(value) || value<lo || value>hi)
            throw std::invalid_argument("Unresolved exact vertical surface moment sum.");
        result.values[i]=value;
        result.bounds[i]={{lo,hi}};
    }
    result.work=work.used;
    return result;
}

ExactSectionMoments exact_section_volume_moments(const ExPolygons& left, const ExPolygons& right,
    ExactSectionOperation operation, const std::array<coord_t, 3>& cell_start,
    const std::array<coord_t, 3>& cell_size, coord_t bottom, coord_t top,
    std::size_t max_work, const std::function<void()>& checkpoint)
{
    ExactSectionWork work(max_work, checkpoint, __func__);
    const auto charge = [&](std::size_t amount) { work.charge(amount); };
    std::size_t count = 0;
    count_exact_inputs(left, count, charge);
    count_exact_inputs(right, count, charge);
    const auto remaining = [&] {
        if (!work.remaining()) throw ResourceLimit("Exact volume moment operation budget exhausted.");
        return work.remaining();
    };
    const auto a = make_exact_section(left, remaining(), checkpoint);
    charge(a.work);
    const auto b = make_exact_section(right, remaining(), checkpoint);
    charge(b.work);
    const auto combined = boolean_exact_sections(a.section, b.section, operation, remaining(), checkpoint);
    charge(combined.work);
    auto result = exact_section_volume_moments(combined.section, cell_start, cell_size, bottom, top, remaining(), checkpoint);
    charge(result.work);
    result.work = work.used;
    return result;
}

ExactBoundaryDifference certify_exact_boundary_difference(const ExPolygons& original,
    const ExPolygons& represented, unsigned allowance_nm, std::size_t max_work,
    const std::function<void()>& checkpoint)
{
    if (!allowance_nm || allowance_nm>7 || !max_work || max_work>4000000)
        throw std::invalid_argument("Invalid exact boundary difference limits.");
    ExactBoundaryDifference result;
    const auto charge = [&](std::size_t amount) {
        if (checkpoint) checkpoint();
        if (amount>max_work-result.work) throw ResourceLimit("Exact boundary difference work budget exceeded.");
        result.work+=amount;
    };
    std::size_t input_points=0;
    for (const auto* polygons : {&original,&represented}) {
        if (polygons->size()>1024) throw ResourceLimit("Exact boundary difference polygon budget exceeded.");
        for (const auto& polygon : *polygons) {
            if (polygon.holes.size()>1024) throw ResourceLimit("Exact boundary difference hole budget exceeded.");
            const auto scan = [&](const Polygon& ring) {
                if (ring.points.size()>1024-input_points)
                    throw ResourceLimit("Exact boundary difference input budget exceeded.");
                input_points+=ring.points.size();
                charge(ring.points.size()+1);
            };
            scan(polygon.contour);
            for (const auto& hole : polygon.holes) scan(hole);
        }
    }
    const auto bound = [](const Set& set) {
        if (set.arrangement().number_of_vertices()>4096)
            throw ResourceLimit("Exact boundary difference arrangement budget exceeded.");
    };
    const auto make_set = [&](const ExPolygons& polygons) {
        if (result.work==max_work) throw ResourceLimit("Exact boundary difference construction budget exhausted.");
        std::size_t used=0;
        auto set=make_topology_union(polygons,max_work-result.work,checkpoint,used,true);
        charge(used);
        bound(set);
        return set;
    };
    const auto before=make_set(original), after=make_set(represented);
    exact_component_correspondence(before,after,charge);
    const auto difference = [&](const Set& a,const Set& b) {
        const auto edges=a.arrangement().number_of_edges()+b.arrangement().number_of_edges();
        if (edges>4096) throw ResourceLimit("Exact boundary difference overlay budget exceeded.");
        charge(edges*edges+1);
        Set output=a;
        output.difference(b);
        charge(0);
        bound(output);
        return output;
    };
    const auto pieces = [&](const Set& set) {
        charge(set.arrangement().number_of_vertices()+1);
        std::list<PolygonHoles> polygons;
        set.polygons_with_holes(std::back_inserter(polygons));
        charge(0);
        return polygons;
    };
    const auto native=pieces(before);
    const auto area_bound = [&](const Set& changed) {
        certify_boundary_band(changed,native,allowance_nm,"exact boundary difference",charge);
        Kernel::FT area=0;
        for (const auto& polygon : pieces(changed)) {
            charge(polygon.outer_boundary().size()+1);
            area+=polygon.outer_boundary().area();
            for (const auto& hole : polygon.holes()) {
                charge(hole.size()+1);
                area+=hole.area();
            }
        }
        if (area<0) throw std::invalid_argument("Invalid exact boundary difference area.");
        // 1e18 is an exactly representable integer. Convert only after the
        // exact unit conversion, retaining the outward interval endpoint.
        return CGAL::to_interval(area/Kernel::FT(1e18)).second;
    };
    result.added_area_bound_m2=area_bound(difference(after,before));
    result.removed_area_bound_m2=area_bound(difference(before,after));
    return result;
}

ExactHoleChangeDiagnostic diagnose_exact_hole_change(const ExPolygons& original,
    const ExPolygons& represented, const Polygon& hole_hint, std::size_t max_work,
    const std::function<void()>& checkpoint)
{
    if (!max_work || max_work>4000000)
        throw std::invalid_argument("Invalid exact hole diagnostic budget.");
    ExactHoleChangeDiagnostic result;
    const auto charge = [&](std::size_t amount) {
        if (checkpoint) checkpoint();
        if (amount>max_work-result.work) throw ResourceLimit("Exact hole diagnostic work budget exceeded.");
        result.work+=amount;
    };
    std::size_t count=0;
    const auto scan = [&](const Polygon& ring) {
        if (ring.points.size()>1024-count) throw ResourceLimit("Exact hole diagnostic input budget exceeded.");
        count+=ring.points.size();
        charge(ring.points.size()+1);
    };
    scan(hole_hint);
    for (const auto* polygons : {&original,&represented}) {
        if (polygons->size()>1024) throw ResourceLimit("Exact hole diagnostic polygon budget exceeded.");
        for (const auto& polygon : *polygons) {
            if (polygon.holes.size()>1024) throw ResourceLimit("Exact hole diagnostic ring budget exceeded.");
            scan(polygon.contour);
            for (const auto& hole : polygon.holes) scan(hole);
        }
    }
    const auto bound = [](const Set& set) {
        if (set.arrangement().number_of_vertices()>4096)
            throw ResourceLimit("Exact hole diagnostic arrangement budget exceeded.");
    };
    const auto make_set = [&](const ExPolygons& polygons) {
        if (result.work==max_work) throw ResourceLimit("Exact hole diagnostic construction budget exhausted.");
        std::size_t used=0;
        auto set=make_topology_union(polygons,max_work-result.work,checkpoint,used,true);
        charge(used);
        bound(set);
        return set;
    };
    const auto before=make_set(original), after=make_set(represented);
    auto hint=hole_hint;
    hint.reverse();
    const auto hint_set=make_set({ExPolygon(std::move(hint))});
    const auto overlay = [&](const Set& a,const Set& b,bool intersection) {
        const auto edges=a.arrangement().number_of_edges()+b.arrangement().number_of_edges();
        if (edges>4096) throw ResourceLimit("Exact hole diagnostic overlay budget exceeded.");
        charge(edges*edges+1);
        Set output=a;
        if (intersection) output.intersection(b);
        else output.symmetric_difference(b);
        charge(0);
        bound(output);
        return output;
    };
    const auto holes = [&](const Set& set) {
        charge(set.arrangement().number_of_vertices()+1);
        std::list<PolygonHoles> polygons;
        set.polygons_with_holes(std::back_inserter(polygons));
        std::vector<Set> regions;
        std::size_t retained=0;
        for (const auto& polygon : polygons)
            for (const auto& hole : polygon.holes()) {
                if (hole.size()>4096-retained) throw ResourceLimit("Exact hole diagnostic retained boundary budget exceeded.");
                retained+=hole.size();
                charge(hole.size()*hole.size()+hole.size()+1);
                auto contour=hole;
                contour.reverse_orientation();
                regions.emplace_back(contour);
            }
        return regions;
    };
    const auto source_holes=holes(before), represented_holes=holes(after);
    const auto unique_match = [&](const Set& candidate,const std::vector<Set>& regions) {
        std::optional<std::size_t> found;
        for (std::size_t i=0;i<regions.size();++i)
            if (!overlay(candidate,regions[i],true).is_empty()) {
                if (found) return std::optional<std::size_t>{};
                found=i;
            }
        return found;
    };
    const auto source=unique_match(hint_set,source_holes), represented_index=unique_match(hint_set,represented_holes);
    result.source_hole_unique=source.has_value();
    result.represented_hole_unique=represented_index.has_value();
    if (!source || !represented_index) return result;
    const auto& source_hole=source_holes[*source];
    const auto& represented_hole=represented_holes[*represented_index];
    if (unique_match(source_hole,represented_holes)!=represented_index ||
        unique_match(represented_hole,source_holes)!=source) return result;
    result.correspondence_unique=true;
    result.hole_unchanged=overlay(source_hole,represented_hole,false).is_empty();
    const auto changed=overlay(before,after,false);
    result.difference_empty=changed.is_empty();
    if (*result.difference_empty) return result;
    Kernel::FT distance=0;
    if (overlay(source_hole,changed,true).is_empty()) {
        std::optional<Kernel::FT> nearest;
        for (auto a=source_hole.arrangement().edges_begin();a!=source_hole.arrangement().edges_end();++a)
            for (auto b=changed.arrangement().edges_begin();b!=changed.arrangement().edges_end();++b) {
                charge(5);
                const auto p=a->source()->point(), q=a->target()->point();
                const auto r=b->source()->point(), s=b->target()->point();
                Kernel::FT candidate=0;
                if (!Kernel().do_intersect_2_object()(Kernel::Segment_2(p,q),Kernel::Segment_2(r,s)))
                    candidate=std::min({segment_distance_squared(p,r,s),segment_distance_squared(q,r,s),
                        segment_distance_squared(r,p,q),segment_distance_squared(s,p,q)});
                if (!nearest || candidate<*nearest) nearest=candidate;
            }
        if (!nearest) return result;
        distance=*nearest;
    }
    const auto interval=CGAL::to_interval(distance);
    result.change_distance_squared_nm2=std::array<double,2>{{interval.first,interval.second}};
    return result;
}

ExactOccupiedAssignment certify_occupied_assignment(const ExPolygons& input,
    const ExPolygons& occupied, const BoundingBox& tile, std::size_t max_work,
    const std::function<void()>& checkpoint)
{
    if (!max_work || max_work>4000000 || tile.min.x()>=tile.max.x() || tile.min.y()>=tile.max.y())
        throw std::invalid_argument("Invalid occupied assignment certificate limits.");
    ExactOccupiedAssignment result;
    const auto charge = [&](std::size_t amount) {
        if (checkpoint) checkpoint();
        if (amount>max_work-result.work) throw ResourceLimit("Occupied assignment work budget exceeded.");
        result.work+=amount;
    };
    std::size_t input_points=0;
    for (const auto* polygons : {&input,&occupied}) {
        if (polygons->size()>1024) throw ResourceLimit("Occupied assignment polygon budget exceeded.");
        for (const auto& polygon : *polygons) {
            if (polygon.holes.size()>1024) throw ResourceLimit("Occupied assignment hole budget exceeded.");
            const auto scan = [&](const Polygon& ring) {
                if (ring.points.size()>1024-input_points) throw ResourceLimit("Occupied assignment input budget exceeded.");
                input_points+=ring.points.size();
                for (const auto& p : ring.points) {
                    charge(1);
                    if (p.x()<tile.min.x() || p.x()>tile.max.x() || p.y()<tile.min.y() || p.y()>tile.max.y())
                        throw std::invalid_argument("Occupied assignment geometry leaves its tile.");
                }
            };
            scan(polygon.contour);
            for (const auto& hole : polygon.holes) scan(hole);
        }
    }
    const auto make_set = [&](const ExPolygons& polygons) {
        if (result.work==max_work) throw ResourceLimit("Occupied assignment construction budget exceeded.");
        std::size_t used=0;
        auto set=make_topology_union(polygons,max_work-result.work,checkpoint,used,true);
        charge(used);
        return set;
    };
    const auto components = [&](const Set& set) {
        charge(set.arrangement().number_of_vertices()+1);
        std::list<PolygonHoles> pieces;
        set.polygons_with_holes(std::back_inserter(pieces));
        return pieces;
    };
    const auto overlay = [&](const Set& a,const Set& b,ExactSectionOperation operation) {
        const auto edges=a.arrangement().number_of_edges()+b.arrangement().number_of_edges();
        if (edges>4096) throw ResourceLimit("Occupied assignment overlay boundary budget exceeded.");
        charge(edges*edges+1);
        Set output=a;
        switch (operation) {
        case ExactSectionOperation::Difference: output.difference(b); break;
        case ExactSectionOperation::Intersection: output.intersection(b); break;
        case ExactSectionOperation::SymmetricDifference: output.symmetric_difference(b); break;
        case ExactSectionOperation::Union: output.join(b); break;
        }
        if (output.arrangement().number_of_vertices()>4096)
            throw ResourceLimit("Occupied assignment overlay output budget exceeded.");
        if (checkpoint) checkpoint();
        return output;
    };
    const auto build_pieces = [&](const ExPolygons& polygons,bool single_component) {
        std::vector<Set> sets;
        std::size_t retained_vertices=0;
        for (const auto& polygon : polygons) {
            auto piece=make_set({polygon});
            retained_vertices+=piece.arrangement().number_of_vertices();
            if (retained_vertices>4096) throw ResourceLimit("Occupied assignment retained boundary budget exceeded.");
            const auto parts=components(piece);
            if (parts.empty() || (single_component && parts.size()!=1))
                throw std::invalid_argument("Occupied assignment has empty or split representation.");
            for (std::size_t previous=0;previous<sets.size();++previous) {
                const auto overlap=overlay(piece,sets[previous],ExactSectionOperation::Intersection);
                if (!overlap.is_empty()) {
                    Kernel::FT area=0;
                    for (const auto& part : components(overlap)) {
                        charge(part.outer_boundary().size()+1);
                        area+=part.outer_boundary().area();
                        for (auto hole=part.holes_begin();hole!=part.holes_end();++hole) {
                            charge(hole->size()+1);
                            area+=hole->area();
                        }
                    }
                    auto vertex=overlap.arrangement().vertices_begin();
                    auto bounds=vertex->point().bbox();
                    for (;vertex!=overlap.arrangement().vertices_end();++vertex) {
                        charge(1);
                        bounds=bounds+vertex->point().bbox();
                    }
                    std::ostringstream message;
                    message.precision(17);
                    message << "Occupied assignment pieces overlap; role=" << (single_component ? "represented" : "original")
                            << " pieces=" << previous << ',' << sets.size() << " area_nm2=" << CGAL::to_double(area)
                            << " bbox_nm=" << bounds.xmin() << ',' << bounds.ymin() << ',' << bounds.xmax() << ',' << bounds.ymax();
                    throw OccupiedAssignmentOverlap(message.str(),!single_component,previous,sets.size(),result.work);
                }
            }
            sets.push_back(std::move(piece));
        }
        return sets;
    };
    const auto originals=build_pieces(input,false), represented=build_pieces(occupied,true);
    const auto source_union=make_set(input), represented_union=make_set(occupied);
    exact_component_correspondence(source_union,represented_union,charge);
    const auto source_components=components(source_union);
    // The full regularized symmetric difference, including its closure, must
    // avoid every tile face. A2D distance band alone cannot bound face traces
    // at a shallow boundary angle.
    const auto changed=overlay(source_union,represented_union,ExactSectionOperation::SymmetricDifference);
    for (auto point=changed.arrangement().vertices_begin();point!=changed.arrangement().vertices_end();++point) {
        charge(1);
        if (point->point().x()<=Kernel::FT(double(tile.min.x())) || point->point().x()>=Kernel::FT(double(tile.max.x())) ||
            point->point().y()<=Kernel::FT(double(tile.min.y())) || point->point().y()>=Kernel::FT(double(tile.max.y())))
            throw std::invalid_argument("Occupied assignment changes a tile-face trace.");
    }
    certify_boundary_band(changed,source_components,1,"occupied representation",charge);
    std::vector<Set> exact_components;
    for (const auto& piece : source_components) {
        std::size_t count=piece.outer_boundary().size();
        for (auto h=piece.holes_begin();h!=piece.holes_end();++h) count+=h->size();
        charge(count*count+count+1);
        exact_components.emplace_back(piece);
    }
    std::vector<std::size_t> target;
    for (const auto& component : exact_components) {
        std::size_t found=represented.size();
        for (std::size_t i=0;i<represented.size();++i)
            if (!overlay(component,represented[i],ExactSectionOperation::Intersection).is_empty()) {
                if (found!=represented.size()) throw std::invalid_argument("Occupied component correspondence is ambiguous.");
                found=i;
            }
        if (found==represented.size()) throw std::invalid_argument("Occupied component correspondence is missing.");
        target.push_back(found);
    }
    for (const auto& original : originals) {
        std::size_t found=exact_components.size();
        for (std::size_t i=0;i<exact_components.size();++i)
            if (overlay(original,exact_components[i],ExactSectionOperation::Difference).is_empty()) {
                if (found!=exact_components.size()) throw std::invalid_argument("Original piece component assignment is ambiguous.");
                found=i;
            }
        if (found==exact_components.size()) throw std::invalid_argument("Original piece spans occupied components.");
        result.occupied_indices.push_back(target[found]);
    }
    return result;
}

ExactIntegerSectionUnion exact_integer_section_union(const ExPolygons& polygons,
    std::size_t max_work, const std::function<void()>& checkpoint)
{
    if (!max_work || max_work > 4000000) throw std::invalid_argument("Invalid exact integer union work budget.");
    std::size_t count = 0;
    if (polygons.size() > 1024) throw ResourceLimit("Exact integer union polygon budget exceeded.");
    for (const auto& polygon : polygons) {
        if (polygon.holes.size() > 1024) throw ResourceLimit("Exact integer union hole budget exceeded.");
        const auto scan = [&](const Polygon& ring) {
            if (checkpoint) checkpoint();
            if (ring.points.size() > 1024-count) throw ResourceLimit("Exact integer union vertex budget exceeded.");
            count += ring.points.size();
        };
        scan(polygon.contour);
        for (const auto& hole : polygon.holes) scan(hole);
    }
    ExactIntegerSectionUnion result;
    auto occupied = make_topology_union(polygons,max_work,checkpoint,result.work,true);
    const auto output_vertices = occupied.arrangement().number_of_vertices();
    if (output_vertices > 4096 || output_vertices > max_work-result.work)
        throw ResourceLimit("Exact integer union output budget exceeded.");
    std::list<PolygonHoles> pieces;
    occupied.polygons_with_holes(std::back_inserter(pieces));
    std::size_t emitted_vertices=0;
    const auto coordinate = [&](const Kernel::FT& value) {
        const double candidate = std::round(CGAL::to_double(value));
        if (!std::isfinite(candidate) || std::abs(candidate)>1e10 || Kernel::FT(candidate)!=value)
            throw NonIntegralSectionUnion(result.work);
        return coord_t(candidate);
    };
    const auto ring = [&](const Polygon2& boundary) {
        Polygon result_ring;
        for (const auto& point : boundary) {
            if (checkpoint) checkpoint();
            // Point-touching components can repeat an arrangement vertex in
            // several output rings. Bound payload, not just unique vertices.
            if (emitted_vertices==4096 || result.work==max_work)
                throw ResourceLimit("Exact integer union emitted vertex budget exceeded.");
            ++emitted_vertices;
            ++result.work;
            result_ring.points.emplace_back(coordinate(point.x()),coordinate(point.y()));
        }
        return result_ring;
    };
    for (const auto& piece : pieces) {
        if (piece.is_unbounded()) throw std::invalid_argument("Unbounded exact occupied union.");
        ExPolygon output(ring(piece.outer_boundary()));
        for (auto hole=piece.holes_begin();hole!=piece.holes_end();++hole)
            output.holes.push_back(ring(*hole));
        result.polygons.push_back(std::move(output));
    }
    return result;
}
ExactSectionPredicate exact_section_empty(const ExPolygons& a, const ExPolygons& b,
    ExactSectionOperation operation, std::size_t max_vertices, std::size_t max_work,
    const std::function<void()>& checkpoint, bool native_nonzero_fill)
{
    if (!max_vertices || max_vertices > 2048 || !max_work || max_work > 4000000 ||
        (operation != ExactSectionOperation::Intersection && operation != ExactSectionOperation::Difference &&
         operation != ExactSectionOperation::SymmetricDifference && operation != ExactSectionOperation::Union))
        throw std::invalid_argument("Invalid exact section predicate limits.");
    std::size_t vertices = 0;
    const auto preflight = [&](const ExPolygons& polygons) {
        if (polygons.size() > max_vertices) throw ResourceLimit("Exact section polygon budget exceeded.");
        for (const auto& polygon : polygons) {
            if (polygon.holes.size() > max_vertices) throw ResourceLimit("Exact section hole budget exceeded.");
            const auto ring = [&](const Polygon& boundary) {
                if (checkpoint) checkpoint();
                if (boundary.points.size() < 3) throw std::invalid_argument("Invalid exact section boundary.");
                if (boundary.points.size() > max_vertices - vertices) throw ResourceLimit("Exact section vertex budget exceeded.");
                vertices += boundary.points.size();
                for (const auto& p : boundary.points)
                    if (std::abs(double(p.x())) > 1e10 || std::abs(double(p.y())) > 1e10)
                        throw std::invalid_argument("Exact section coordinate exceeds the supported range.");
            };
            ring(polygon.contour);
            for (const auto& hole : polygon.holes) ring(hole);
        }
    };
    preflight(a);
    preflight(b);
    if (native_nonzero_fill) {
        if (vertices > 1024) throw ResourceLimit("Exact native predicate vertex budget exceeded.");
        std::size_t work = vertices + 1;
        if (work > max_work) throw ResourceLimit("Exact native predicate work budget exhausted.");
        // Identical ordered chains define the same winding everywhere, even
        // when their touching boundaries need reconstruction for other queries.
        if ((operation == ExactSectionOperation::Difference || operation == ExactSectionOperation::SymmetricDifference) && a == b)
            return {true, work};
        const auto make_set = [&](const ExPolygons& polygons) {
            if (work >= max_work) throw ResourceLimit("Exact native predicate work budget exhausted.");
            std::size_t used = 0;
            auto set = make_topology_union(polygons, max_work-work, checkpoint, used, true);
            work += used;
            return set;
        };
        auto x = make_set(a), y = make_set(b);
        const std::size_t edges = x.arrangement().number_of_edges() + y.arrangement().number_of_edges();
        if (edges > 4096) throw ResourceLimit("Exact native predicate boundary budget exceeded.");
        const auto overlay_work = edges*edges+1;
        if (overlay_work > max_work-work) throw ResourceLimit("Exact native predicate overlay work budget exceeded.");
        work += overlay_work;
        if (checkpoint) checkpoint();
        switch (operation) {
        case ExactSectionOperation::Intersection: x.intersection(y); break;
        case ExactSectionOperation::Difference: x.difference(y); break;
        case ExactSectionOperation::SymmetricDifference: x.symmetric_difference(y); break;
        case ExactSectionOperation::Union: x.join(y); break;
        }
        if (checkpoint) checkpoint();
        return {x.is_empty(), work};
    }
    const std::size_t work = vertices * vertices + vertices + 1;
    if (work > max_work) throw ResourceLimit("Exact section arrangement work budget exceeded.");
    const auto left = convert_polygons(a, checkpoint), right = convert_polygons(b, checkpoint);
    if (checkpoint) checkpoint();
    if (operation == ExactSectionOperation::Intersection) {
        bool candidate = false;
        for (const auto& p : left) for (const auto& q : right) {
            if (checkpoint) checkpoint();
            const auto x = p.outer_boundary().bbox(), y = q.outer_boundary().bbox();
            if (std::min(x.xmax(), y.xmax()) > std::max(x.xmin(), y.xmin()) &&
                std::min(x.ymax(), y.ymax()) > std::max(x.ymin(), y.ymin())) candidate = true;
        }
        if (!candidate) return {true, work};
    }
    Set x, y;
    x.join(left.begin(), left.end());
    if (checkpoint) checkpoint();
    y.join(right.begin(), right.end());
    if (checkpoint) checkpoint();
    switch (operation) {
    case ExactSectionOperation::Intersection: x.intersection(y); break;
    case ExactSectionOperation::Difference: x.difference(y); break;
    case ExactSectionOperation::SymmetricDifference: x.symmetric_difference(y); break;
    case ExactSectionOperation::Union: x.join(y); break;
    }
    if (checkpoint) checkpoint();
    return {x.is_empty(), work};
}
ExactIntersectionOwnership exact_intersection_ownership(const ExPolygons& first_claim, const ExPolygons& second_claim,
    const ExPolygons& first_native, const ExPolygons& second_native, std::size_t max_work,
    const std::function<void()>& checkpoint)
{
    if (!max_work || max_work > 4000000) throw std::invalid_argument("Invalid exact ownership budget.");
    ExPolygons all;
    std::size_t count = 0;
    for (const auto* polygons : {&first_claim, &second_claim, &first_native, &second_native}) {
        if (checkpoint) checkpoint();
        if (polygons->size() > 1024) throw ResourceLimit("Exact ownership polygon budget exceeded.");
        for (const auto& polygon : *polygons)
            if (polygon.holes.size() > 1024) throw ResourceLimit("Exact ownership hole budget exceeded.");
        const auto n = count_points(*polygons);
        if (n > 1024 - count) throw ResourceLimit("Exact ownership vertex budget exceeded.");
        count += n;
        all.insert(all.end(), polygons->begin(), polygons->end());
    }
    if (max_work < 13) throw ResourceLimit("Exact ownership work budget exceeded.");
    // Reserve four operand unions and nine bounded Boolean operations before
    // constructing any. The empty
    // intersection validates original rings and hole placement without output
    // conversion, and supplies the same conservative quadratic reservation.
    const auto validation = exact_section_empty(all, {}, ExactSectionOperation::Intersection, 1024, max_work / 13, checkpoint);
    ExactIntersectionOwnership result;
    result.work = 13 * validation.work;
    const auto make_set = [&](const ExPolygons& polygons) {
        const auto converted = convert_polygons(polygons, checkpoint);
        Set set;
        set.join(converted.begin(), converted.end());
        if (checkpoint) checkpoint();
        return set;
    };
    Set overlap = make_set(first_claim);
    overlap.intersection(make_set(second_claim));
    if (checkpoint) checkpoint();
    if (overlap.is_empty()) return result;
    result.nonempty = true;
    const auto a = make_set(first_native), b = make_set(second_native);
    const auto unique = [&](const Set& owner, const Set& other) {
        Set missing = overlap;
        missing.difference(owner);
        if (checkpoint) checkpoint();
        if (!missing.is_empty()) return false;
        Set competing = overlap;
        competing.intersection(other);
        if (checkpoint) checkpoint();
        return competing.is_empty();
    };
    result.first_unique = unique(a, b);
    result.second_unique = unique(b, a);
    Set covered = a;
    covered.join(b);
    if (checkpoint) checkpoint();
    Set missing = overlap;
    missing.difference(covered);
    if (checkpoint) checkpoint();
    Set competing = a;
    competing.intersection(b);
    if (checkpoint) checkpoint();
    competing.intersection(overlap);
    if (checkpoint) checkpoint();
    result.native_covers = missing.is_empty();
    result.native_disjoint = competing.is_empty();
    result.native_partition = result.native_covers && result.native_disjoint;
    return result;
}

static UnqualifiedInterfaceVariants partition_interface_impl(const ExPolygons& first_claim, const ExPolygons& second_claim,
    const ExPolygons& first_native, const ExPolygons& second_native, unsigned allowance_nm,
    std::size_t max_work, const std::function<void()>& input_checkpoint,
    const std::array<ExPolygons, 2>* complete_envelopes, NativeWitnessFill native_fill,
    bool allow_single_sided_native_boundary, bool activation,
    const std::function<void(std::size_t)>& charge_work,
    ExactInterfacePartition* exact_output = nullptr, ExactInterfaceChoice exact_choice = ExactInterfaceChoice::FavorFirst,
    bool allow_dual_native_boundary = false)
{
    bool external_failure = false;
    const std::function<void()> checkpoint = [&] {
        if (input_checkpoint) try { input_checkpoint(); } catch (...) { external_failure = true; throw; }
    };
    if (allowance_nm == 0 || allowance_nm > 7 || max_work == 0 || max_work > 4000000 ||
        (native_fill != NativeWitnessFill::ValidatedSimple && native_fill != NativeWitnessFill::NonzeroWinding))
        throw std::invalid_argument("Invalid analysis interface uncertainty limits.");
    if (allow_single_sided_native_boundary && !complete_envelopes)
        throw std::invalid_argument("Single-sided boundary analysis requires complete deposited envelopes.");
    if (allow_dual_native_boundary && (!exact_output || !complete_envelopes || activation || allowance_nm!=7))
        throw std::invalid_argument("Dual-boundary sensitivity requires complete exact envelopes and the explicit7nm scenario.");
    if (exact_output && (activation || !complete_envelopes ||
        (exact_choice != ExactInterfaceChoice::FavorFirst && exact_choice != ExactInterfaceChoice::FavorSecond)))
        throw std::invalid_argument("Invalid exact interface partition request.");
    UnqualifiedInterfaceVariants output;
    double dual_native_boundary_area_m2=0.;
    std::size_t native_winding_work=0;
    std::string phase = "input";
    const auto charge = [&](std::size_t work, const char* stage = "interface proof") {
        if (checkpoint) checkpoint();
        if (work > max_work - output.work)
            throw ResourceLimit(std::string("Analysis interface work budget exceeded: ") + stage +
                ", used=" + std::to_string(output.work) + ", requested=" + std::to_string(work) +
                ", limit=" + std::to_string(max_work) + ", phase=" + phase + ".");
        if (charge_work) try { charge_work(work); } catch (...) { external_failure = true; throw; }
        output.work += work;
    };
    std::vector<const ExPolygons*> inputs{&first_claim, &second_claim, &first_native, &second_native};
    if (complete_envelopes) {
        if (!first_claim.empty() || !second_claim.empty())
            throw std::invalid_argument("Complete interface envelopes cannot accompany precomputed claims.");
        inputs.push_back(&(*complete_envelopes)[0]);
        inputs.push_back(&(*complete_envelopes)[1]);
    }
    // Dense local masks may exceed1024 vertices without exceeding the
    // independent4M operation budget. Construction still reserves its
    // quadratic work before CGAL; this is not an unbounded overlay allowance.
    constexpr std::size_t operand_limit = 2048;
    for (const auto* polygons : inputs) {
        if (polygons->size() > operand_limit) throw ResourceLimit("Analysis interface polygon budget exceeded.");
        for (const auto& polygon : *polygons)
            if (polygon.holes.size() > operand_limit) throw ResourceLimit("Analysis interface hole budget exceeded.");
        const auto count = count_points(*polygons);
        if (count > operand_limit)
            throw ResourceLimit(std::string(activation ? "Reference activation" : "Analysis interface") +
                " input vertex budget exceeded (operand=" + std::to_string(count) +
                ", limit=" + std::to_string(operand_limit) + ").");
        // Check cheap input bounds here. Ring and hole validity is checked by
        // convert_polygons during make_set below, once per operand. An empty
        // intersection solely for validation would convert every ring twice.
        charge(count + 1, "input validation");
        for (const auto& polygon : *polygons) {
            const auto validate_ring = [&](const Polygon& ring) {
                if (ring.points.size() < 3) throw std::invalid_argument("Invalid exact section boundary.");
                for (const auto& point : ring.points) {
                    if (checkpoint) checkpoint();
                    if (std::abs(double(point.x())) > 1e10 || std::abs(double(point.y())) > 1e10)
                        throw std::invalid_argument("Exact section coordinate exceeds the supported range.");
                }
            };
            validate_ring(polygon.contour);
            for (const auto& hole : polygon.holes) validate_ring(hole);
        }
    }
    const auto bound = [&](const Set& set) {
        const auto count = set.arrangement().number_of_vertices();
        if (count > 4096) throw ResourceLimit("Analysis interface arrangement vertex budget exceeded.");
        return std::size_t(count);
    };
    std::optional<std::array<Kernel::FT,4>> native_clip;
    if (complete_envelopes && native_fill==NativeWitnessFill::NonzeroWinding) {
        std::optional<std::array<coord_t,4>> bounds;
        for (const auto& envelope : *complete_envelopes)
            for (const auto& polygon : envelope) {
                const auto include = [&](const Polygon& ring) {
                    for (const auto& p : ring.points) {
                        charge(1,"native crop bounds");
                        if (!bounds) bounds=std::array<coord_t,4>{{p.x(),p.x(),p.y(),p.y()}};
                        else {
                            (*bounds)[0]=std::min((*bounds)[0],p.x());
                            (*bounds)[1]=std::max((*bounds)[1],p.x());
                            (*bounds)[2]=std::min((*bounds)[2],p.y());
                            (*bounds)[3]=std::max((*bounds)[3],p.y());
                        }
                    }
                };
                include(polygon.contour);
                for (const auto& hole : polygon.holes) include(hole);
            }
        if (bounds) {
            // All claims and uncovered-conflict proofs lie in the complete
            // envelopes. Artificial crop edges stay strictly farther away
            // than the admitted interface band and cannot supply band proof.
            const coord_t padding=coord_t(allowance_nm)+2;
            native_clip=std::array<Kernel::FT,4>{{Kernel::FT(double((*bounds)[0]-padding)),
                Kernel::FT(double((*bounds)[1]+padding)),Kernel::FT(double((*bounds)[2]-padding)),
                Kernel::FT(double((*bounds)[3]+padding))}};
        }
    }
    const auto make_set = [&](const ExPolygons& polygons, bool native = false) {
        if (native && native_fill == NativeWitnessFill::NonzeroWinding) {
            if (output.work == max_work) throw ResourceLimit("Analysis native winding work budget exceeded.");
            std::size_t construction = 0;
            auto set = native_clip ? make_clipped_native_union(polygons,*native_clip,max_work-output.work,checkpoint,construction) :
                make_topology_union(polygons,max_work-output.work,checkpoint,construction,true);
            charge(construction,"native winding construction");
            native_winding_work+=construction;
            bound(set);
            return set;
        }
        const auto count = count_points(polygons);
        charge(count * count + count + 1,"operand construction");
        const auto converted = convert_polygons(polygons, checkpoint);
        Set set;
        set.join(converted.begin(), converted.end());
        if (checkpoint) checkpoint();
        bound(set);
        return set;
    };
    enum class Operation { Union, Intersection, Difference, SymmetricDifference };
    const auto combine = [&](const Set& a, const Set& b, Operation operation) {
        // Empty-set identities need no arrangement overlay. Charge the copy
        // that remains, rather than a quadratic operation we do not perform.
        if (a.is_empty() || b.is_empty()) {
            if (operation == Operation::Intersection || (operation == Operation::Difference && a.is_empty())) {
                charge(1);
                return Set{};
            }
            const auto& nonempty = a.is_empty() ? b : a;
            charge(bound(nonempty)+1);
            return Set(nonempty);
        }
        const auto count = bound(a) + bound(b);
        if (exact_output && complete_envelopes) {
            reserve_exact_overlay_work(a,b,[&](std::size_t amount) { charge(amount,"Boolean event preflight"); });
        } else {
            // Both operands are already planar arrangements. Only edges from
            // different operands can create new intersections. Charge the cheap
            // broad-phase tests too, then reserve both endpoints of each possible
            // collinear overlap. CGAL's outward-rounded boxes only reject pairs
            // that provably cannot intersect; the actual overlay remains exact.
            const std::size_t edges_a = a.arrangement().number_of_edges();
            const std::size_t edges_b = b.arrangement().number_of_edges();
            charge(edges_a * edges_b + edges_a + edges_b + count + 1, "Boolean broad phase");
            std::size_t possible_intersections = 0;
            for (auto i=a.arrangement().edges_begin();i!=a.arrangement().edges_end();++i) {
                if (checkpoint) checkpoint();
                const auto box_a=i->source()->point().bbox()+i->target()->point().bbox();
                for (auto j=b.arrangement().edges_begin();j!=b.arrangement().edges_end();++j) {
                    const auto box_b=j->source()->point().bbox()+j->target()->point().bbox();
                    if (box_a.xmin()<=box_b.xmax() && box_b.xmin()<=box_a.xmax() &&
                        box_a.ymin()<=box_b.ymax() && box_b.ymin()<=box_a.ymax())
                        ++possible_intersections;
                }
            }
            charge(2 * possible_intersections, "Boolean overlay");
        }
        Set result = a;
        switch (operation) {
        case Operation::Union: result.join(b); break;
        case Operation::Intersection: result.intersection(b); break;
        case Operation::Difference: result.difference(b); break;
        case Operation::SymmetricDifference: result.symmetric_difference(b); break;
        }
        if (checkpoint) checkpoint();
        bound(result);
        return result;
    };
    const auto polygons = [&](const Set& set) {
        charge(bound(set) + 1);
        std::list<PolygonHoles> result;
        set.polygons_with_holes(std::back_inserter(result));
        if (checkpoint) checkpoint();
        return result;
    };
    phase = "native coverage";
    const Set na = make_set(first_native,true), nb = make_set(second_native,true);
    std::array<std::array<Set, 2>, 2> variants;
    Set overlap;
    if (activation) {
        phase = "reference activation";
        if (!combine(na,nb,Operation::Intersection).is_empty())
            throw std::invalid_argument("Exact reference activation materials overlap.");
        variants[0] = {combine(make_set(first_claim),na,Operation::Intersection),
                       combine(make_set(second_claim),nb,Operation::Intersection)};
        // Exact subsets of disjoint reference materials are disjoint. The
        // empty overlap also selects one shared rounding/certification pass.
    } else {
    std::array<Set, 2> envelopes;
    Set conflict;
    Set uncovered_native;
    std::string coverage_error;
    if (complete_envelopes) {
        envelopes = {make_set((*complete_envelopes)[0]), make_set((*complete_envelopes)[1])};
        conflict = combine(envelopes[0], envelopes[1], Operation::Intersection);
        // C \ (A union B) == (C \ A) \ B. Work on the local
        // conflict rather than constructing the complete native union.
        uncovered_native = combine(combine(conflict, na, Operation::Difference), nb, Operation::Difference);
        if (!uncovered_native.is_empty()) {
            const auto pieces = polygons(uncovered_native);
            Kernel::FT area(0);
            std::optional<CGAL::Bbox_2> bounds;
            for (const auto& piece : pieces) {
                area += piece.outer_boundary().area();
                for (const auto& hole : piece.holes()) area += hole.area();
                const auto box = piece.outer_boundary().bbox();
                bounds = bounds ? *bounds+box : box;
            }
            coverage_error = "Native interface masks do not cover the complete deposited conflict: area_nm2=" +
                std::to_string(CGAL::to_double(area)) + ", pieces=" + std::to_string(pieces.size()) +
                (bounds ? ", bounds_nm=" + std::to_string(bounds->xmin()) + "," + std::to_string(bounds->ymin()) +
                    " -> " + std::to_string(bounds->xmax()) + "," + std::to_string(bounds->ymax()) : "") + ".";
            if (!pieces.empty()) {
                const auto& piece=pieces.front();
                coverage_error += " First piece area_nm2=" + std::to_string(CGAL::to_double(piece.outer_boundary().area())) +
                    ", holes=" + std::to_string(piece.number_of_holes()) + ", vertices_nm=";
                for (std::size_t i=0;i<std::min<std::size_t>(piece.outer_boundary().size(),8);++i)
                    coverage_error += "(" + std::to_string(CGAL::to_double(piece.outer_boundary()[i].x())) + "," +
                        std::to_string(CGAL::to_double(piece.outer_boundary()[i].y())) + ")";
            }
        }
    }
    Set a = complete_envelopes ? combine(conflict, na, Operation::Intersection) : make_set(first_claim);
    Set b = complete_envelopes ? combine(conflict, nb, Operation::Intersection) : make_set(second_claim);
    overlap = combine(a, b, Operation::Intersection);
    if (overlap.is_empty() && !complete_envelopes) throw std::invalid_argument("Analysis interface requires an actual claim overlap.");
    // With complete envelopes, a and b were intersected with the native
    // masks above. Their overlap is therefore already covered by both;
    // repeating a native-union overlay proves nothing additional.
    if (!complete_envelopes &&
        !combine(combine(overlap, na, Operation::Difference), nb, Operation::Difference).is_empty())
        throw std::invalid_argument("Analysis interface uncertainty cannot assign native void.");
    // These are boundaries of the exact material UNION, not internal edges of
    // overlapping same-material input pieces.
    const auto native_a = polygons(na), native_b = polygons(nb);
    const auto distance_squared = [&](const Kernel::Point_2& point, const Kernel::Point_2& p, const Kernel::Point_2& q) -> Kernel::FT {
        charge(1);
        return segment_distance_squared(point,p,q);
    };
    const auto certify_band = [&](const Set& changed, const std::list<PolygonHoles>& native, unsigned radius_nm,
                                  const std::string& stage) {
        certify_boundary_band(changed,native,radius_nm,stage,charge);
    };
    if (!uncovered_native.is_empty() && (allow_single_sided_native_boundary || allow_dual_native_boundary)) {
        phase="single-sided native boundary";
        const Kernel::FT radius2=Kernel::FT(int(allowance_nm))*Kernel::FT(int(allowance_nm));
        // Failure of the containing-band proof is NOT proof of separation.
        // Independently exclude every segment capsule of the other material.
        const auto separated=[&](const PolygonHoles& piece,const std::list<PolygonHoles>& native) {
            for (const auto& reference : native) {
                const auto ring_separated=[&](const Polygon2& boundary) {
                    for (std::size_t j=0;j<boundary.size();++j) {
                        const auto& p=boundary[j];
                        const auto& q=boundary[(j+1)%boundary.size()];
                        charge(piece.outer_boundary().size()+1);
                        // Conservative also when this point lies in a hole.
                        if (piece.outer_boundary().bounded_side(p)!=CGAL::ON_UNBOUNDED_SIDE) return false;
                        const auto edges_separated=[&](const Polygon2& ring) {
                            for (std::size_t i=0;i<ring.size();++i) {
                                charge(5);
                                const auto& u=ring[i];
                                const auto& v=ring[(i+1)%ring.size()];
                                if (Kernel().do_intersect_2_object()(Kernel::Segment_2(u,v),Kernel::Segment_2(p,q)) ||
                                    segment_distance_squared(u,p,q)<=radius2 || segment_distance_squared(v,p,q)<=radius2 ||
                                    segment_distance_squared(p,u,v)<=radius2 || segment_distance_squared(q,u,v)<=radius2)
                                    return false;
                            }
                            return true;
                        };
                        if (!edges_separated(piece.outer_boundary())) return false;
                        for (auto h=piece.holes_begin();h!=piece.holes_end();++h)
                            if (!edges_separated(*h)) return false;
                    }
                    return true;
                };
                if (!ring_separated(reference.outer_boundary())) return false;
                for (auto h=reference.holes_begin();h!=reference.holes_end();++h)
                    if (!ring_separated(*h)) return false;
            }
            return true;
        };
        for (const auto& piece : polygons(uncovered_native)) {
            std::size_t count=piece.outer_boundary().size();
            Kernel::FT area=piece.outer_boundary().area();
            for (auto h=piece.holes_begin();h!=piece.holes_end();++h) { count+=h->size(); area+=h->area(); }
            charge(count*count+count+1);
            Set fragment;
            fragment.insert(piece);
            bound(fragment);
            const bool away_first=separated(piece,native_a),away_second=separated(piece,native_b);
            if (!away_first && !away_second && allow_dual_native_boundary) {
                phase="dual native boundary sensitivity";
                // This fragment is already in BOTH deposited envelopes and
                // outside BOTH native references. Its entire exact area must
                // fit both existing bands; separation failures alone prove no
                // coverage. The two priority scenarios assign it differently.
                certify_band(fragment,native_a,allowance_nm,"dual/first native");
                certify_band(fragment,native_b,allowance_nm,"dual/second native");
                a=combine(a,fragment,Operation::Union);
                b=combine(b,fragment,Operation::Union);
                overlap=combine(overlap,fragment,Operation::Union);
                dual_native_boundary_area_m2+=CGAL::to_double(area)*1e-18;
                continue;
            }
            if (away_first==away_second) {
                // Report the already constructed failing fragment, not another
                // overlay or an inferred owner. False means separation was not
                // proved; it does not certify containment in a boundary band.
                const auto bounds=piece.outer_boundary().bbox();
                std::string details="Single-sided native boundary ownership is ambiguous or unsupported: separation="+
                    std::string(away_first ? "away-from-both" : "not-separated-from-either")+
                    ", away_first="+(away_first ? "true" : "false")+
                    ", away_second="+(away_second ? "true" : "false")+
                    ", native_fill="+(native_fill==NativeWitnessFill::NonzeroWinding ? "nonzero-winding" : "validated-simple")+
                    ", variant="+(!exact_output ? "both" : exact_choice==ExactInterfaceChoice::FavorFirst ? "favor-first" : "favor-second")+
                    ", allowance_nm="+std::to_string(allowance_nm)+
                    ", used_work="+std::to_string(output.work)+", max_work="+std::to_string(max_work)+
                    ", fragment_area_nm2="+std::to_string(CGAL::to_double(area))+
                    ", fragment_holes="+std::to_string(piece.number_of_holes())+
                    ", fragment_vertices="+std::to_string(count)+
                    ", fragment_bounds_nm="+std::to_string(bounds.xmin())+","+std::to_string(bounds.ymin())+
                    " -> "+std::to_string(bounds.xmax())+","+std::to_string(bounds.ymax())+
                    ", fragment_first_vertices_nm=";
                for (std::size_t i=0;i<std::min<std::size_t>(piece.outer_boundary().size(),8);++i)
                    details+="("+std::to_string(CGAL::to_double(piece.outer_boundary()[i].x()))+","+
                        std::to_string(CGAL::to_double(piece.outer_boundary()[i].y()))+")";
                throw std::invalid_argument(details+". "+coverage_error);
            }
            if (!allow_single_sided_native_boundary)
                throw std::invalid_argument("Dual-boundary sensitivity does not authorize single-sided native ownership.");
            certify_band(fragment,away_first ? native_b : native_a,allowance_nm,"single-sided native ownership");
            if (away_first) b=combine(b,fragment,Operation::Union);
            else a=combine(a,fragment,Operation::Union);
            output.single_sided_native_boundary_area_m2+=CGAL::to_double(area)*1e-18;
        }
        // Single-sided additions leave a intersect b unchanged. Dual additions
        // explicitly extend it above so even an initially empty native overlap
        // retains both sensitivity choices. No deposited exclusive is altered.
        uncovered_native=Set{};
    }
    if (!uncovered_native.is_empty()) {
        // Diagnostic only: even complete coverage by both uncertainty bands
        // does not authorize filling native void or selecting its material.
        // Prove the whole region, not merely its vertices or small area.
        phase = "uncovered native band";
        try {
            certify_band(uncovered_native,native_a,allowance_nm,"uncovered/first native");
            certify_band(uncovered_native,native_b,allowance_nm,"uncovered/second native");
        } catch (const std::invalid_argument& error) {
            if (external_failure) throw;
            std::string details=coverage_error + " " + error.what();
            try {
                certify_band(uncovered_native,native_b,allowance_nm,"uncovered/second native");
                details += " Second native boundary band certified.";
            } catch (const std::invalid_argument& boundary_error) {
                if (external_failure) throw;
                details += " " + std::string(boundary_error.what());
            }
            // Distinguish native-mask disagreement from a sliver at the
            // deposited envelope boundary. This is diagnostic evidence only:
            // even a certified thin overlap still has no native owner.
            for (unsigned side=0;side<2;++side) {
                phase="uncovered deposited band diagnostic";
                try {
                    certify_band(uncovered_native,polygons(envelopes[side]),allowance_nm,
                        side==0 ? "uncovered/first deposited" : "uncovered/second deposited");
                    details += " Deposited boundary band " + std::to_string(side) + " certified.";
                } catch (const std::invalid_argument& boundary_error) {
                    if (external_failure) throw;
                    details += " " + std::string(boundary_error.what());
                }
            }
            throw std::invalid_argument(details);
        }
        throw std::invalid_argument(coverage_error + " Entire uncovered region lies in both native boundary bands; ownership remains unresolved.");
    }
    phase = "interface band";
    certify_band(overlap, native_a, allowance_nm, "overlap/first native");
    certify_band(overlap, native_b, allowance_nm, "overlap/second native");
    phase = "variant construction";
    const bool selected_only=exact_output && complete_envelopes;
    if (selected_only) {
        // Preserve all source/band admission above, but construct only the
        // requested partition. Both full variants remain required by the
        // legacy API; an exact caller publishes only this selected pair.
        if (overlap.is_empty()) {
            charge(bound(a)+bound(b)+1);
            variants[0]={a,b};
        } else if (exact_choice==ExactInterfaceChoice::FavorFirst) {
            charge(bound(a)+1);
            variants[0]={a,combine(b,a,Operation::Difference)};
        } else {
            charge(bound(b)+1);
            variants[0]={combine(a,b,Operation::Difference),b};
        }
    } else {
        // The already computed exact intersection also proves these differences
        // are identities when empty (single-sided additions cannot change it).
        if (overlap.is_empty()) charge(bound(a)+bound(b)+1);
        const auto exclusive_a = overlap.is_empty() ? a : combine(a, b, Operation::Difference);
        const auto exclusive_b = overlap.is_empty() ? b : combine(b, a, Operation::Difference);
        if (complete_envelopes) {
            // With complete envelopes, priority assignments reduce to
            // (a,b\a) and (a\b,b). Exact-only dual sensitivity may additionally
            // include certified fragments outside both native references.
            charge(bound(a)+bound(b)+bound(exclusive_a)+bound(exclusive_b)+1);
            variants={{{a,exclusive_b},{exclusive_a,b}}};
        } else {
            const auto owned_a = combine(overlap, na, Operation::Intersection), owned_b = combine(overlap, nb, Operation::Intersection);
            variants={{{combine(exclusive_a, owned_a, Operation::Union),
                        combine(exclusive_b, combine(overlap, owned_a, Operation::Difference), Operation::Union)},
                       {combine(exclusive_a, combine(overlap, owned_b, Operation::Difference), Operation::Union),
                        combine(exclusive_b, owned_b, Operation::Union)}}};
        }
    }
    if (complete_envelopes) {
        const std::array<Set,2> deposited_exclusive{
            combine(envelopes[0],conflict,Operation::Difference),
            combine(envelopes[1],conflict,Operation::Difference)};
        for (std::size_t v=0;v<(selected_only || overlap.is_empty() ? 1u : 2u);++v)
            for (std::size_t material = 0; material < 2; ++material)
                variants[v][material] = combine(deposited_exclusive[material],variants[v][material],Operation::Union);
        // The second identical variant is copied only after certification.
        if (exact_output) {
            phase = "exact partition publication";
            auto& pair = variants[0];
            if (!combine(pair[0],pair[1],Operation::Intersection).is_empty())
                throw std::invalid_argument("Exact interface materials are not disjoint.");
            const auto occupied = combine(pair[0],pair[1],Operation::Union);
            const auto deposited = combine(envelopes[0],envelopes[1],Operation::Union);
            if (!combine(occupied,deposited,Operation::SymmetricDifference).is_empty())
                throw std::invalid_argument("Exact interface partition changes the deposited union.");
            const auto retained = bound(pair[0]) + bound(pair[1]);
            if (retained > 4096) throw ResourceLimit("Exact interface partition retained vertex budget exceeded.");
            ExactInterfacePartition certified;
            const auto publication_charge = [&](std::size_t amount) { charge(amount,"exact partition publication"); };
            certified.materials[0] = ExactSectionAccess::make(std::move(pair[0]),publication_charge);
            certified.materials[1] = ExactSectionAccess::make(std::move(pair[1]),publication_charge);
            certified.interface_allowance_nm = allowance_nm;
            // No integer rounding has occurred. Keep this distinct from the
            // inherited native-interface sensitivity scenario.
            certified.representation_allowance_nm = 0;
            certified.single_sided_native_boundary_area_m2 = output.single_sided_native_boundary_area_m2;
            certified.dual_native_boundary_area_m2 = dual_native_boundary_area_m2;
            certified.native_winding_work = native_winding_work;
            certified.work = output.work;
            *exact_output = std::move(certified);
            return output;
        }
    }
    }
    std::map<Kernel::Point_2, Point> shared_vertices;
    std::array<double, 2> rounding_bias{0., 0.};
    const auto round_coordinate = [](const Kernel::FT& original, double bias) {
        const Kernel::FT value = original + Kernel::FT(bias);
        coord_t rounded = coord_t(std::floor(CGAL::to_double(value) + .5));
        // Correct the floating estimate with exact comparisons, including ties
        // and negative coordinates. Coordinates stay inside the input bounds.
        if (value < Kernel::FT(double(rounded)) - Kernel::FT(.5)) --rounded;
        if (value >= Kernel::FT(double(rounded)) + Kernel::FT(.5)) ++rounded;
        if (value < Kernel::FT(double(rounded)) - Kernel::FT(.5) || value >= Kernel::FT(double(rounded)) + Kernel::FT(.5))
            throw std::invalid_argument("Analysis interface coordinate rounding is unresolved.");
        return rounded;
    };
    const auto round_set = [&](const Set& set, const std::vector<Kernel::Point_2>& interface_vertices, bool& unchanged) {
        unchanged = true;
        ExPolygons result;
        const auto ring = [&](const Polygon2& boundary) {
            Polygon rounded;
            std::vector<Kernel::Point_2> split_boundary;
            for (std::size_t i=0;i<boundary.size();++i) {
                const auto& start=boundary[i];
                const auto& end=boundary[(i+1)%boundary.size()];
                const Kernel::Segment_2 edge(start,end);
                std::vector<Kernel::Point_2> interior;
                // A shared straight edge may have a T-junction vertex in
                // only one material's arrangement. Split it on BOTH sides
                // before rounding, so both use the same piecewise chain.
                // The shared set is lexicographically sorted: restrict tests
                // to its exact x interval rather than scanning all vertices.
                std::size_t search_work=1;
                for (auto n=interface_vertices.size();n;n/=2) ++search_work;
                charge(2*search_work);
                const auto begin=std::lower_bound(interface_vertices.begin(),interface_vertices.end(),
                    std::min(start.x(),end.x()),[](const auto& point,const auto& x) { return point.x()<x; });
                const auto finish=std::upper_bound(begin,interface_vertices.end(),
                    std::max(start.x(),end.x()),[](const auto& x,const auto& point) { return x<point.x(); });
                charge(std::size_t(std::distance(begin,finish))+1);
                const auto box=start.bbox()+end.bbox();
                for (auto vertex=begin;vertex!=finish;++vertex) {
                    const auto& point=*vertex;
                    const auto candidate=point.bbox();
                    if (candidate.xmin()>box.xmax() || candidate.xmax()<box.xmin() ||
                        candidate.ymin()>box.ymax() || candidate.ymax()<box.ymin()) continue;
                    charge(1);
                    if (point!=start && point!=end && edge.has_on(point)) interior.push_back(point);
                }
                charge(interior.size()*interior.size()+1);
                std::sort(interior.begin(),interior.end(),[&](const auto& a,const auto& b) {
                    return (a-start).squared_length()<(b-start).squared_length();
                });
                split_boundary.push_back(start);
                split_boundary.insert(split_boundary.end(),interior.begin(),interior.end());
                if (split_boundary.size()>4096) throw ResourceLimit("Analysis interface split boundary exceeds its vertex budget.");
            }
            for (const auto& p : split_boundary) {
                charge(1);
                auto found = shared_vertices.find(p);
                if (found == shared_vertices.end()) {
                    if (shared_vertices.size() == 4096) throw ResourceLimit("Analysis interface output vertex budget exceeded.");
                    found = shared_vertices.emplace(p, Point(round_coordinate(p.x(), rounding_bias[0]),
                        round_coordinate(p.y(), rounding_bias[1]))).first;
                }
                const auto dx = p.x() - Kernel::FT(double(found->second.x()));
                const auto dy = p.y() - Kernel::FT(double(found->second.y()));
                unchanged = unchanged && dx == 0 && dy == 0;
                if (dx * dx + dy * dy > Kernel::FT(activation ? 1. : .5))
                    throw std::invalid_argument("Analysis interface rounding exceeds its representation bound.");
                // Preserve the ordered edge correspondence, including holes.
                // During linear vertex homotopy each swept edge quadrilateral
                // lies in the original edge's convex 1 nm capsule: both moved
                // endpoints are within 1 nm of their originals. Winding
                // membership can change only where a boundary sweeps, proving
                // the whole represented-set change is inside those capsules.
                // Removing a consecutive duplicate only removes a zero mapped
                // edge; it does not alter that boundary chain or its winding.
                if (rounded.points.empty() || rounded.points.back() != found->second) rounded.points.push_back(found->second);
            }
            if (rounded.points.size() > 1 && rounded.points.front() == rounded.points.back()) rounded.points.pop_back();
            if (rounded.points.size() < 3) throw std::invalid_argument("Analysis interface rounding collapses a boundary.");
            return rounded;
        };
        for (const auto& piece : polygons(set)) {
            result.emplace_back(ring(piece.outer_boundary()));
            for (auto h = piece.holes_begin(); h != piece.holes_end(); ++h) result.back().holes.push_back(ring(*h));
        }
        return result;
    };
    const auto piece_set = [&](const PolygonHoles& piece) {
        std::size_t count = piece.outer_boundary().size();
        for (auto hole = piece.holes_begin(); hole != piece.holes_end(); ++hole) count += hole->size();
        charge(count * count + count + 1);
        Set set;
        set.insert(piece);
        if (checkpoint) checkpoint();
        bound(set);
        return set;
    };
    const auto match_pieces = [&](const std::vector<Set>& before, const std::vector<Set>& after) {
        if (before.size() != after.size())
            throw std::invalid_argument("Analysis interface rounding changes component or hole topology.");
        std::vector<std::size_t> matches;
        std::vector<bool> used(after.size(), false);
        const auto bounds = [&](const std::vector<Set>& pieces) {
            std::vector<CGAL::Bbox_2> result;
            result.reserve(pieces.size());
            for (const auto& piece : pieces) {
                std::optional<CGAL::Bbox_2> box;
                for (auto vertex = piece.arrangement().vertices_begin(); vertex != piece.arrangement().vertices_end(); ++vertex) {
                    charge(1);
                    const auto point_box = vertex->point().bbox();
                    box = box ? *box + point_box : point_box;
                }
                if (!box) throw std::invalid_argument("Analysis interface has an empty topology component.");
                result.push_back(*box);
            }
            return result;
        };
        const auto before_bounds = bounds(before), after_bounds = bounds(after);
        const auto intersects = [&](const Set& a, const Set& b) {
            // A point strictly interior to both exact sets certifies positive
            // area overlap without constructing the intersection arrangement.
            // Try bounded inward edge offsets and rational centroids, then fall back to the
            // full Boolean predicate. Failure to find a witness proves nothing.
            std::array<Kernel::Point_2, 3> points;
            std::size_t visited = 0;
            const std::size_t point_work = bound(a) + bound(b) + 1;
            const auto shared_interior = [&](const Kernel::Point_2& point) {
                charge(point_work);
                return a.oriented_side(point) == CGAL::ON_POSITIVE_SIDE &&
                       b.oriented_side(point) == CGAL::ON_POSITIVE_SIDE;
            };
            // Component centroids can lie in a hole or outside a thin curved
            // region. Small exact offsets from its oriented outer boundary
            // also sample the material, without assuming any minimum width.
            for (const auto& piece : polygons(a)) {
                const auto& ring = piece.outer_boundary();
                for (std::size_t i = 0; i < std::min<std::size_t>(ring.size(), 2); ++i) {
                    const auto& p = ring[i];
                    const auto& q = ring[(i + 1) % ring.size()];
                    const auto dx = q.x() - p.x(), dy = q.y() - p.y();
                    const auto norm = CGAL::abs(dx) + CGAL::abs(dy);
                    if (norm == 0) continue;
                    for (const int quarter_nm : {4, 64}) {
                        const Kernel::FT step = Kernel::FT(quarter_nm) / 4;
                        if (shared_interior(Kernel::Point_2((p.x() + q.x()) / 2 - dy * step / norm,
                                                           (p.y() + q.y()) / 2 + dx * step / norm)))
                            return true;
                    }
                }
            }
            for (auto vertex = a.arrangement().vertices_begin();
                 vertex != a.arrangement().vertices_end() && visited < 4; ++vertex) {
                points[visited % 3] = vertex->point();
                if (++visited < 3) continue;
                const Kernel::Point_2 candidate((points[0].x() + points[1].x() + points[2].x()) / 3,
                                               (points[0].y() + points[1].y() + points[2].y()) / 3);
                if (shared_interior(candidate)) return true;
            }
            return !combine(a, b, Operation::Intersection).is_empty();
        };
        for (std::size_t source_index = 0; source_index < before.size(); ++source_index) {
            const auto& source = before[source_index];
            std::size_t match = after.size();
            for (std::size_t i = 0; i < after.size(); ++i) {
                charge(1);
                const auto& a = before_bounds[source_index];
                const auto& b = after_bounds[i];
                // CGAL's outward-rounded boxes contain the exact rational
                // coordinates. Disjoint boxes cannot have positive-area
                // overlap; touching boxes also need no polygon overlay.
                if (std::min(a.xmax(), b.xmax()) <= std::max(a.xmin(), b.xmin()) ||
                    std::min(a.ymax(), b.ymax()) <= std::max(a.ymin(), b.ymin()))
                    continue;
                phase = "rounded topology match " + std::to_string(source_index) + "/" + std::to_string(before.size()) +
                    " vs " + std::to_string(i) + "/" + std::to_string(after.size()) +
                    " vertices=" + std::to_string(bound(source)) + "/" + std::to_string(bound(after[i]));
                if (intersects(source, after[i])) {
                    if (match != after.size() || used[i])
                        throw std::invalid_argument("Analysis interface rounding has ambiguous component correspondence.");
                    match = i;
                }
            }
            if (match == after.size())
                throw std::invalid_argument("Analysis interface rounding loses a component or hole.");
            used[match] = true;
            matches.push_back(match);
        }
        return matches;
    };
    const auto qualify_topology = [&](const Set& before, const Set& after) {
        const auto a_list = polygons(before), b_list = polygons(after);
        const std::vector<PolygonHoles> a(a_list.begin(), a_list.end()), b(b_list.begin(), b_list.end());
        std::vector<Set> a_sets, b_sets;
        for (const auto& piece : a) a_sets.push_back(piece_set(piece));
        for (const auto& piece : b) b_sets.push_back(piece_set(piece));
        const auto matches = match_pieces(a_sets, b_sets);
        const auto holes = [&](const PolygonHoles& piece) {
            std::vector<Set> result;
            for (auto h = piece.holes_begin(); h != piece.holes_end(); ++h) {
                auto contour = *h;
                contour.reverse_orientation();
                result.push_back(piece_set(PolygonHoles(contour)));
            }
            return result;
        };
        for (std::size_t i = 0; i < a.size(); ++i) match_pieces(holes(a[i]), holes(b[matches[i]]));
    };
    // Global biases preserve integer vertices and use one shared mapping for
    // both materials. After the 1/8 nm family, sparse one-axis quarter phases
    // can resolve a bent shared chain without simultaneously biasing both
    // coordinates further. Every candidate retains the exact 1 nm norm check.
    constexpr std::array<double, 3> biases{0., -.125, .125};
    constexpr std::array<std::array<double, 2>, 4> quarter_biases{{
        {{.25,0.}}, {{-.25,0.}}, {{0.,.25}}, {{0.,-.25}}}};
    const std::size_t attempts = activation ? 13 : 1;
    bool certified = false;
    std::string last_geometric_failure;
    for (std::size_t trial = 0; trial < attempts; ++trial) {
        rounding_bias = trial < 9 ? std::array<double, 2>{biases[trial / 3], biases[trial % 3]} : quarter_biases[trial - 9];
        shared_vertices.clear();
        output.favor_first = {};
        output.favor_second = {};
        try {
            for (std::size_t v = 0; v < 2; ++v) {
                if (v==1 && overlap.is_empty()) {
                    // With disjoint exact native claims there is no priority choice:
                    // both expressions above are the same sets. Reuse the first
                    // certified rounding/topology result instead of redoing overlays.
                    charge(count_points(output.favor_first[0])+count_points(output.favor_first[1])+1);
                    output.favor_second=output.favor_first;
                    break;
                }
                phase = "rounded topology";
                auto& destination = v == 0 ? output.favor_first : output.favor_second;
                std::set<Kernel::Point_2> unique_vertices;
                for (const auto& material : variants[v]) {
                    charge(bound(material)+1);
                    for (auto i=material.arrangement().vertices_begin();i!=material.arrangement().vertices_end();++i)
                        unique_vertices.insert(i->point());
                }
                if (unique_vertices.size()>4096) throw ResourceLimit("Analysis interface shared vertices exceed their budget.");
                const std::vector<Kernel::Point_2> interface_vertices(unique_vertices.begin(),unique_vertices.end());
                std::array<Set, 2> rounded;
                bool both_unchanged = true;
                for (std::size_t material = 0; material < 2; ++material) {
                    bool unchanged = false;
                    destination[material] = round_set(variants[v][material], interface_vertices, unchanged);
                    both_unchanged = both_unchanged && unchanged;
                    if (unchanged) {
                        // Exact integer vertices leave every oriented edge and hole
                        // unchanged. Reuse the already valid arrangement: rebuilding
                        // it and matching all its components would prove identity
                        // again, at quadratic cost. Rational changes still require
                        // the full topology check below.
                        charge(bound(variants[v][material]) + 1);
                        rounded[material] = variants[v][material];
                    } else {
                        rounded[material] = make_set(destination[material]);
                        qualify_topology(variants[v][material], rounded[material]);
                    }
                }
                phase = "rounded disjointness";
                // The exact priority partitions are disjoint by construction. When
                // neither boundary moved, another overlay proves the same identity.
                if (!both_unchanged) {
                    const auto rounded_overlap = combine(rounded[0], rounded[1], Operation::Intersection);
                    if (!rounded_overlap.is_empty()) {
                        if (trial + 1 != attempts) throw std::invalid_argument("Analysis interface rounded materials overlap.");
                        std::ostringstream detail;
                        detail.precision(17);
                        auto vertex = rounded_overlap.arrangement().vertices_begin();
                        auto bounds = vertex->point().bbox();
                        for (; vertex != rounded_overlap.arrangement().vertices_end(); ++vertex) {
                            charge(1);
                            bounds = bounds + vertex->point().bbox();
                        }
                        detail << "Analysis interface rounded materials overlap. bias_nm=" << rounding_bias[0] << ',' << rounding_bias[1]
                               << " overlap_bbox_nm=" << bounds.xmin() << ',' << bounds.ymin() << ',' << bounds.xmax() << ',' << bounds.ymax()
                               << " overlap_vertices=" << rounded_overlap.arrangement().number_of_vertices();
                        for (std::size_t material = 0; material < 2; ++material) {
                            std::size_t emitted = 0;
                            for (auto edge = variants[v][material].arrangement().edges_begin();
                                 edge != variants[v][material].arrangement().edges_end(); ++edge) {
                                charge(1);
                                const auto a = edge->source()->point(), b = edge->target()->point();
                                const auto box = a.bbox() + b.bbox();
                                if (box.xmax() < bounds.xmin()-2 || box.xmin() > bounds.xmax()+2 ||
                                    box.ymax() < bounds.ymin()-2 || box.ymin() > bounds.ymax()+2) continue;
                                if (emitted++ == 8) break;
                                const auto ra = shared_vertices.find(a), rb = shared_vertices.find(b);
                                detail << " material=" << material << " exact_edge_nm=" << CGAL::to_double(a.x()) << ',' << CGAL::to_double(a.y())
                                       << ':' << CGAL::to_double(b.x()) << ',' << CGAL::to_double(b.y());
                                if (ra != shared_vertices.end() && rb != shared_vertices.end())
                                    detail << " rounded_edge_nm=" << ra->second.x() << ',' << ra->second.y() << ':' << rb->second.x() << ',' << rb->second.y();
                            }
                        }
                        throw std::invalid_argument(detail.str());
                    }
                }
            }
            certified = true;
            break;
        } catch (const ResourceLimit&) {
            throw;
        } catch (const std::invalid_argument& error) {
            if (!activation || external_failure) throw;
            last_geometric_failure = error.what();
        }
    }
    if (!certified)
        throw std::invalid_argument("Reference activation rounding failed all " + std::to_string(attempts) +
            " bounded phases: " + last_geometric_failure);
    output.interface_allowance_nm = activation ? 0 : allowance_nm;
    output.representation_allowance_nm = 1;
    output.vertex_rounding_bound_m = 1e-9; // Exact per-vertex norm certificate; unchanged public bound.
    return output;
}

UnqualifiedInterfaceVariants partition_uncertain_interface(const ExPolygons& first_claim, const ExPolygons& second_claim,
    const ExPolygons& first_native, const ExPolygons& second_native, unsigned allowance_nm,
    std::size_t max_work, const std::function<void()>& checkpoint,
    const std::array<ExPolygons, 2>* complete_envelopes, NativeWitnessFill native_fill,
    bool allow_single_sided_native_boundary)
{
    return partition_interface_impl(first_claim,second_claim,first_native,second_native,allowance_nm,max_work,
        checkpoint,complete_envelopes,native_fill,allow_single_sided_native_boundary,false,{});
}

ExactInterfacePartition exact_partition_uncertain_interface(const ExPolygons& first_native,
    const ExPolygons& second_native, const std::array<ExPolygons, 2>& complete_envelopes,
    ExactInterfaceChoice choice, unsigned allowance_nm, std::size_t max_work,
    const std::function<void()>& checkpoint, NativeWitnessFill native_fill,
    bool allow_single_sided_native_boundary, bool allow_dual_native_boundary)
{
    ExactInterfacePartition result;
    partition_interface_impl({}, {}, first_native, second_native, allowance_nm, max_work, checkpoint,
        &complete_envelopes, native_fill, allow_single_sided_native_boundary, false, {}, &result, choice, allow_dual_native_boundary);
    return result;
}

ExactReferenceClipping exact_clip_reference_sections(const std::array<ExPolygons, 2>& reference,
    const std::array<ExPolygons, 2>& deposited, std::size_t max_work,
    const std::function<void()>& checkpoint, const std::function<void(std::size_t)>& charge_work)
{
    auto result=partition_interface_impl(deposited[0],deposited[1],reference[0],reference[1],1,max_work,
        checkpoint,nullptr,NativeWitnessFill::ValidatedSimple,false,true,charge_work);
    return {std::move(result.favor_first),result.representation_allowance_nm,result.vertex_rounding_bound_m,result.work};
}

ExactOwnershipWitness diagnose_ownership_witness(const ExPolygons& first_claim, const ExPolygons& second_claim,
    const ExPolygons& first_native, const ExPolygons& second_native, std::size_t max_work,
    const std::function<void()>& checkpoint, NativeWitnessFill native_fill)
{
    if (!max_work || max_work > 20000000 || (native_fill != NativeWitnessFill::ValidatedSimple &&
        native_fill != NativeWitnessFill::NonzeroWinding)) throw std::invalid_argument("Invalid ownership witness budget.");
    ExactOwnershipWitness result;
    const auto charge = [&](std::size_t n) {
        if (checkpoint) checkpoint();
        if (n > max_work - result.work) throw ResourceLimit("Ownership witness work budget exceeded.");
        result.work += n;
    };
    const auto validate = [&](const ExPolygons& polygons, std::size_t& count, std::size_t cap) {
        if (polygons.size() > cap) throw ResourceLimit("Ownership witness polygon budget exceeded.");
        for (const auto& polygon : polygons) {
            if (polygon.holes.size() > cap) throw ResourceLimit("Ownership witness hole budget exceeded.");
            const auto ring = [&](const Polygon& boundary) {
                if (boundary.points.size() < 3) throw std::invalid_argument("Invalid ownership witness ring.");
                if (boundary.points.size() > cap - count) throw ResourceLimit("Ownership witness vertex budget exceeded.");
                count += boundary.points.size();
                charge(boundary.points.size());
                for (const auto& p : boundary.points)
                    if (std::abs(double(p.x())) > 1e10 || std::abs(double(p.y())) > 1e10)
                        throw std::invalid_argument("Ownership witness coordinate exceeds bounds.");
            };
            ring(polygon.contour);
            for (const auto& hole : polygon.holes) ring(hole);
        }
    };
    std::size_t claims = 0, native = 0;
    validate(first_claim, claims, 1024);
    validate(second_claim, claims, 1024);
    validate(first_native, native, 4500);
    validate(second_native, native, 4500);
    charge(claims * claims + native * native + 1);
    const auto a = convert_polygons(first_claim, checkpoint, "first claim"), b = convert_polygons(second_claim, checkpoint, "second claim");
    const auto na = native_fill == NativeWitnessFill::ValidatedSimple ?
        convert_polygons(first_native, checkpoint, "first native") : std::vector<PolygonHoles>{};
    const auto nb = native_fill == NativeWitnessFill::ValidatedSimple ?
        convert_polygons(second_native, checkpoint, "second native") : std::vector<PolygonHoles>{};
    Set intersection, other;
    intersection.join(a.begin(), a.end());
    if (checkpoint) checkpoint();
    other.join(b.begin(), b.end());
    if (checkpoint) checkpoint();
    intersection.intersection(other);
    if (checkpoint) checkpoint();
    std::list<PolygonHoles> pieces;
    intersection.polygons_with_holes(std::back_inserter(pieces));
    const auto contains = [&](const PolygonHoles& polygon, const Kernel::Point_2& point) -> std::optional<bool> {
        charge(polygon.outer_boundary().size());
        const auto outer = polygon.outer_boundary().bounded_side(point);
        if (outer == CGAL::ON_BOUNDARY) return std::nullopt;
        if (outer == CGAL::ON_UNBOUNDED_SIDE) return false;
        for (auto h = polygon.holes_begin(); h != polygon.holes_end(); ++h) {
            charge(h->size());
            const auto side = h->bounded_side(point);
            if (side == CGAL::ON_BOUNDARY) return std::nullopt;
            if (side == CGAL::ON_BOUNDED_SIDE) return false;
        }
        return true;
    };
    const auto classify = [&](const std::vector<PolygonHoles>& polygons, const Kernel::Point_2& point) -> std::optional<bool> {
        bool boundary = false;
        for (const auto& polygon : polygons) {
            const auto inside = contains(polygon, point);
            if (inside && *inside) return true;
            boundary = boundary || !inside;
        }
        return boundary ? std::nullopt : std::optional<bool>(false);
    };
    const auto nonzero_fill = [&](const ExPolygons& polygons, const Kernel::Point_2& point) -> std::optional<bool> {
        bool on_boundary = false;
        for (const auto& polygon : polygons) {
            int winding = 0;
            bool boundary = false;
            const auto ring = [&](const Polygon& contour) {
                for (std::size_t i = 0; i < contour.points.size(); ++i) {
                    charge(1);
                    const auto& p = contour.points[i];
                    const auto& q = contour.points[(i + 1) % contour.points.size()];
                    const Kernel::FT px(double(p.x())), py(double(p.y())), qx(double(q.x())), qy(double(q.y()));
                    const Kernel::FT cross = (qx - px) * (point.y() - py) - (qy - py) * (point.x() - px);
                    if (cross == 0 && point.x() >= std::min(px, qx) && point.x() <= std::max(px, qx) &&
                        point.y() >= std::min(py, qy) && point.y() <= std::max(py, qy)) boundary = true;
                    if (py <= point.y() && point.y() < qy && cross > 0) ++winding;
                    if (qy <= point.y() && point.y() < py && cross < 0) --winding;
                }
            };
            ring(polygon.contour);
            for (const auto& hole : polygon.holes) ring(hole);
            if (!boundary && winding != 0) return true;
            on_boundary = on_boundary || boundary;
        }
        return on_boundary ? std::nullopt : std::optional<bool>(false);
    };
    for (const auto& piece : pieces) {
        const auto& ring = piece.outer_boundary();
        for (std::size_t i = 0; i < ring.size(); ++i) {
            charge(1);
            const auto& p = ring[i];
            const auto& q = ring[(i + 1) % ring.size()];
            const auto& r = ring[(i + 2) % ring.size()];
            const Kernel::Point_2 witness((p.x() + q.x() + r.x()) / 3, (p.y() + q.y() + r.y()) / 3);
            if (contains(piece, witness) != std::optional<bool>(true)) continue;
            result.in_first_native = native_fill == NativeWitnessFill::ValidatedSimple ? classify(na, witness) : nonzero_fill(first_native, witness);
            result.in_second_native = native_fill == NativeWitnessFill::ValidatedSimple ? classify(nb, witness) : nonzero_fill(second_native, witness);
            return result;
        }
    }
    return result;
}
} // namespace Slic3r::Thermoelastic
