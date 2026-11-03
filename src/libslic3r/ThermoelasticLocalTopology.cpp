#include "ThermoelasticLocalSections.hpp"
#include "ThermoelasticLocalTopology.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include "BoundingBox.hpp"
#include "ClipperUtils.hpp"
#include "ExPolygon.hpp"
#include "Polygon.hpp"
#include "Point.hpp"
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticExactSections.hpp"
#include "ThermoelasticSectionFootprint.hpp"
#include "ThermoelasticMesh.hpp"
#include "libslic3r.h"
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
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
using ExactInteger = boost::multiprecision::int128_t;
ExactInteger turn(const Point& a, const Point& b, const Point& c)
{
    return ExactInteger(b.x() - a.x()) * (c.y() - a.y()) - ExactInteger(b.y() - a.y()) * (c.x() - a.x());
}
bool on_segment(const Point& p, const Point& a, const Point& b)
{
    return turn(a,b,p) == 0 && p.x() >= std::min(a.x(),b.x()) && p.x() <= std::max(a.x(),b.x()) &&
        p.y() >= std::min(a.y(),b.y()) && p.y() <= std::max(a.y(),b.y());
}
bool intersect_segments(const Point& a, const Point& b, const Point& c, const Point& d)
{
    if (std::max(a.x(),b.x()) < std::min(c.x(),d.x()) || std::max(c.x(),d.x()) < std::min(a.x(),b.x()) ||
        std::max(a.y(),b.y()) < std::min(c.y(),d.y()) || std::max(c.y(),d.y()) < std::min(a.y(),b.y())) return false;
    const auto ac = turn(a,b,c), ad = turn(a,b,d), ca = turn(c,d,a), cb = turn(c,d,b);
    if (((ac > 0 && ad < 0) || (ac < 0 && ad > 0)) && ((ca > 0 && cb < 0) || (ca < 0 && cb > 0))) return true;
    return on_segment(c,a,b) || on_segment(d,a,b) || on_segment(a,c,d) || on_segment(b,c,d);
}
// -1 boundary, 0 exterior, 1 nonzero-winding interior; no division/rounding.
int exact_location(const Point& p, const Polygon& ring)
{
    int winding = 0;
    for (std::size_t i = 0; i < ring.points.size(); ++i) {
        const auto& a = ring.points[i];
        const auto& b = ring.points[(i+1) % ring.points.size()];
        if (on_segment(p,a,b)) return -1;
        if (a.y() <= p.y() && b.y() > p.y() && turn(a,b,p) > 0) ++winding;
        if (a.y() > p.y() && b.y() <= p.y() && turn(a,b,p) < 0) --winding;
    }
    return winding != 0;
}
coord_t native_height(double z)
{
    if (!std::isfinite(z) || std::abs(z) > 10.)
        throw std::invalid_argument("Invalid local topology height.");
    return std::llround(z * 1e9);
}
std::size_t validate_local_polygon(const ExPolygon& polygon, const LocalSectionTile& tile,
    const char* stage, coord_t bottom, coord_t top)
{
    const auto twice_area=[](const Polygon& ring) {
        ExactInteger sum=0;
        for (std::size_t i=1;i+1<ring.points.size();++i)
            sum+=turn(ring.points.front(),ring.points[i],ring.points[i+1]);
        return sum;
    };
    // Native-coordinate slivers can have positive exact area while the global
    // double shoelace sum rounds to zero. Validate winding in the same integer
    // arithmetic used by the contact normalizer, without changing geometry.
    ExactInteger area=twice_area(polygon.contour);
    bool valid=polygon.contour.points.size()>=3 && area>0;
    for (const auto& hole : polygon.holes) {
        const auto hole_area=twice_area(hole);
        valid=valid && hole.points.size()>=3 && hole_area<0;
        area+=hole_area;
    }
    if (!valid || area<=0) {
        std::ostringstream message;
        message << "Invalid local topology polygon at " << stage << " tile=" << tile.index[0] << ',' << tile.index[1]
                << " Znm=" << bottom << ".." << top << " twice_area_nm2=" << area;
        throw std::invalid_argument(message.str());
    }
    const auto ring = [&](const Polygon& boundary) {
        for (std::size_t i = 0; i < boundary.points.size(); ++i) {
            const auto& p = boundary.points[i];
            if (p.x() < tile.minimum[0] || p.x() > tile.maximum[0] ||
                p.y() < tile.minimum[1] || p.y() > tile.maximum[1]) {
                std::ostringstream message;
                message << "Local topology polygon escapes tile at " << stage << " tile=" << tile.index[0] << ','
                        << tile.index[1] << " Znm=" << bottom << ".." << top << " point=" << p.x() << ',' << p.y()
                        << " bounds=" << tile.minimum[0] << ',' << tile.minimum[1] << ".."
                        << tile.maximum[0] << ',' << tile.maximum[1];
                throw std::invalid_argument(message.str());
            }
            // Every caller first runs the private exact contact normalizer.
            // A certified outer/hole shared endpoint is deliberately retained.
        }
    };
    ring(polygon.contour);
    for (const auto& hole : polygon.holes) ring(hole);
    std::size_t count = polygon.contour.points.size();
    for (const auto& hole : polygon.holes) count += hole.points.size();
    return count;
}
using NativeTrace = std::vector<std::pair<coord_t, coord_t>>;
std::size_t polygon_points(const ExPolygon& polygon)
{
    std::size_t n = polygon.contour.points.size();
    for (const auto& hole : polygon.holes) n += hole.points.size();
    return n;
}
std::size_t polygon_points(const ExPolygons& polygons)
{
    std::size_t n = 0;
    for (const auto& polygon : polygons) n += polygon_points(polygon);
    return n;
}
NativeTrace boundary_trace(const ExPolygon& polygon, unsigned axis, coord_t coordinate,
    std::size_t maximum = std::numeric_limits<std::size_t>::max())
{
    NativeTrace trace;
    const auto ring = [&](const Polygon& boundary) {
        for (std::size_t i = 0; i < boundary.points.size(); ++i) {
            const auto& a = boundary.points[i];
            const auto& b = boundary.points[(i + 1) % boundary.points.size()];
            if (a[axis] == coordinate && b[axis] == coordinate && a[1 - axis] != b[1 - axis]) {
                if (trace.size() == maximum) throw ResourceLimit("Local topology trace budget exceeded.");
                trace.emplace_back(std::min(a[1 - axis], b[1 - axis]), std::max(a[1 - axis], b[1 - axis]));
            }
        }
    };
    ring(polygon.contour);
    for (const auto& hole : polygon.holes) ring(hole);
    std::sort(trace.begin(), trace.end());
    std::size_t merged = 0;
    for (const auto interval : trace) {
        if (!merged || interval.first > trace[merged - 1].second) trace[merged++] = interval;
        else trace[merged - 1].second = std::max(trace[merged - 1].second, interval.second);
    }
    trace.resize(merged);
    return trace;
}
bool face_contact(const NativeTrace& a, const NativeTrace& b, coord_t a_margin = 0, coord_t b_margin = 0)
{
    std::size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (std::min(a[i].second - a_margin, b[j].second - b_margin) >
            std::max(a[i].first + a_margin, b[j].first + b_margin)) return true;
        if (a[i].second < b[j].second) ++i;
        else ++j;
    }
    return false;
}
} // namespace

struct LocalTopologyBuilder::Impl {
    using Index    = std::array<int, 2>;
    using Interval = std::pair<coord_t, coord_t>;
    struct Trace {
        NativeTrace native;
        std::optional<ExactSectionBoundary> exact;
        Trace() = default;
        Trace(NativeTrace value) : native(std::move(value)) {}
        explicit Trace(ExactSectionBoundary value) : exact(std::move(value)) {}
        std::size_t size() const { return exact ? exact->segment_count() : native.size(); }
        bool empty() const { return size()==0; }
        // Exact records are not coalesced by a rounded or identity comparison.
        bool operator==(const Trace& other) const { return !exact && !other.exact && native==other.native; }
    };
    struct Island
    {
        ExPolygon footprint;
        BoundingBox bounds;
        coord_t bottom, top;
        std::size_t label = std::numeric_limits<std::size_t>::max();
        std::array<Trace, 4> traces;
        unsigned uncertainty = 0;
        std::size_t record = std::numeric_limits<std::size_t>::max();
        std::optional<ExactSection> exact;
    };
    struct Label {
        std::size_t parent;
        bool bed;
        // Diagnostic payload only: merged alongside the actual union-find,
        // never used to authorize contact or stabilize detached material.
        std::array<coord_t, 3> minimum{}, maximum{};
        double volume_m3 = 0.;
        std::size_t pieces = 0, members = 1, geometric_islands = 0;
    };
    struct ContactRecord {
        coord_t bottom, top;
        std::size_t label;
        std::array<Trace,4> traces;
        unsigned uncertainty;
    };
    std::vector<Label> labels;
    struct UncertainContact {
        std::size_t first, second;
        Index tile;
        coord_t height;
    };
    // Stable label IDs, never unioned by these edges. One bounded diagnostic
    // describes the first observation; geometry is not retained per edge.
    std::vector<UncertainContact> uncertain_contacts;
    std::string first_uncertain_contact;
    std::vector<ContactRecord> records;
    std::map<Index, std::vector<std::size_t>> tile_islands;
    LocalTileOrder order;
    LocalAttachmentPolicy attachment;
    std::size_t free_components = 0;
    using Face = std::pair<std::size_t, unsigned>;
    std::map<Index, std::vector<Face>> pending;
    std::array<std::vector<std::size_t>, 2> incoming;
    std::array<std::size_t, 2> incoming_first{};
    std::vector<std::size_t> free_records;
    std::optional<Index> last_tile;
    std::size_t active_records = 0;
    coord_t bed;
    SectionTopologyLimits limits;
    std::optional<StreamTopologyLimits> streaming;
    StreamTopologyStatistics* statistics;
    std::size_t processed = 0, live_points = 0, trace_intervals = 0;
    std::function<void(const SectionPiece&)> piece_callback;
    std::function<void()> checkpoint;
    std::array<coord_t, 2> size{}, origin{};
    std::size_t input_points = 0, retained_points = 0, work = 0, slabs = 0, emitted_pieces = 0, components = 0;
    std::size_t containment_work = 0, containment_queries = 0, overlap_work = 0, diagnostic_work = 0;
    bool tile_open = false;
    LocalSectionTile active_tile;
    std::vector<Island> previous;
    coord_t previous_top = 0;
    enum class State { Ready, Active, Finished, Failed } state = State::Ready;

    Impl(double bed_z,
         const SectionTopologyLimits& input_limits,
         std::function<void(const SectionPiece&)> piece,
         std::function<void()> check, std::optional<StreamTopologyLimits> stream_limits,
         StreamTopologyStatistics* stats, LocalTileOrder tile_order, LocalAttachmentPolicy attachment_policy)
        : order(tile_order), attachment(attachment_policy), bed(native_height(bed_z)), limits(input_limits), streaming(stream_limits), statistics(stats),
          piece_callback(std::move(piece)), checkpoint(std::move(check))
    {
        if (order != LocalTileOrder::Arbitrary && order != LocalTileOrder::Lexicographic)
            throw std::invalid_argument("Invalid local tile order.");
        if (attachment != LocalAttachmentPolicy::RequireBed && attachment != LocalAttachmentPolicy::AllowFreeComponents)
            throw std::invalid_argument("Invalid local attachment policy.");
        if (statistics) *statistics = {};
        if (streaming && (!streaming->max_processed_points || streaming->max_processed_points > 128000000 ||
            !streaming->max_work || streaming->max_work > 4000000000ULL || !streaming->max_live_points ||
            streaming->max_live_points > 2000000 || !streaming->max_islands || streaming->max_islands > 250000 ||
            !streaming->max_trace_intervals || streaming->max_trace_intervals > 2000000 ||
            !streaming->max_slabs || streaming->max_slabs > 250000))
            throw std::invalid_argument("Invalid streaming topology limits.");
        if (!piece_callback || !limits.max_slabs || limits.max_slabs > 10000 || !limits.max_boundary_points ||
            limits.max_boundary_points > 2000000 || !limits.max_component_pairs || limits.max_component_pairs > 2000000)
            throw std::invalid_argument("Invalid local topology input or limits.");
    }
    void initialize_grid(const LocalSectionTile& tile)
    {
        for (unsigned axis = 0; axis < 2; ++axis) {
            if (std::abs(double(tile.minimum[axis])) > 1e10 || std::abs(double(tile.maximum[axis])) > 1e10)
                throw std::invalid_argument("Local topology tile coordinate exceeds the supported range.");
            size[axis] = tile.maximum[axis] - tile.minimum[axis];
            if (size[axis] <= 0 || size[axis] > coord_t(1e10))
                throw std::invalid_argument("Invalid local topology tile dimensions.");
            const long double proposed_origin = static_cast<long double>(tile.minimum[axis]) -
                                                static_cast<long double>(tile.index[axis]) * size[axis];
            if (std::abs(proposed_origin) > 1e10L)
                throw std::invalid_argument("Local topology grid origin exceeds the supported range.");
            origin[axis] = coord_t(proposed_origin);
        }
    }
    void tick(std::size_t amount = 1)
    {
        if (checkpoint)
            checkpoint();
        if (amount > (streaming ? streaming->max_work : limits.max_component_pairs) - work)
            throw ResourceLimit("Local topology work budget exceeded.");
        work += amount;
        if (statistics) statistics->work = work;
    }
    void charge(std::size_t& value, std::size_t n)
    {
        if (n > (streaming ? streaming->max_processed_points - processed : limits.max_boundary_points - value))
            throw ResourceLimit("Local topology point budget exceeded.");
        value += n;
        processed += n;
        if (statistics) {
            statistics->processed_points = processed;
            statistics->input_points = input_points;
        }
    }
    void reserve_live(std::size_t n)
    {
        if (streaming && n > streaming->max_live_points - live_points)
            throw ResourceLimit("Local topology live polygon point budget exceeded.");
        live_points += n;
        if (statistics) statistics->peak_live_points = std::max(statistics->peak_live_points, live_points);
    }
    struct TemporaryPoints {
        Impl& owner;
        std::size_t count;
        TemporaryPoints(Impl& o, std::size_t n) : owner(o), count(n) { owner.reserve_live(n); }
        ~TemporaryPoints() { owner.live_points -= count; }
    };
    std::size_t exact_budget() const
    {
        const auto available=(streaming ? streaming->max_work : limits.max_component_pairs)-work;
        if (!available) throw ResourceLimit("Exact local topology work budget exhausted.");
        return std::min<std::size_t>(available,4000000);
    }
    ExactSection canonical(const SectionFootprint& footprint)
    {
        if (footprint.is_exact()) return footprint.exact();
        TemporaryPoints scratch(*this,4096);
        auto result=make_exact_section(footprint.native(),exact_budget(),[&] { tick(0); });
        tick(result.work);
        return std::move(result.section);
    }
    ExactSection combine(const ExactSection& a,const ExactSection& b,ExactSectionOperation operation)
    {
        TemporaryPoints scratch(*this,4096+a.vertex_count()+b.vertex_count());
        auto result=boolean_exact_sections(a,b,operation,exact_budget(),[&] { tick(0); });
        tick(result.work);
        return std::move(result.section);
    }
    void validate_canonical(const ExactSection& section)
    {
        // Immutable canonical inputs are certified when material and occupied
        // unions are decomposed below. Do not reject their boundary pinches
        // before the explicit interior policy can separate unwelded faces.
        if (limits.exact_interior_components) return;
        // The bounded validator retains extracted rings and their exact vertex
        // set concurrently, in addition to the caller's canonical handle.
        TemporaryPoints scratch(*this,8192);
        tick(validate_exact_section_topology(section,exact_budget(),[&] { tick(0); }));
    }
    std::size_t island_points(const Island& island) const
    {
        return island.exact ? island.exact->vertex_count() : polygon_points(island.footprint);
    }
    void discard_trace(Trace& trace)
    {
        if (trace.exact) live_points-=2*trace.size();
        trace={};
    }
    bool contact_traces(const Trace& p,const Trace& q,unsigned axis,unsigned p_uncertainty,unsigned q_uncertainty)
    {
        if (!p.exact && !q.exact) {
            const bool contact=face_contact(p.native,q.native);
            if (p_uncertainty || q_uncertainty) {
                const coord_t sign=contact ? 1 : -1;
                if (contact!=face_contact(p.native,q.native,sign*p_uncertainty,sign*q_uncertainty))
                    throw std::invalid_argument("Tile-face contact is unresolved within boundary uncertainty.");
            }
            return contact;
        }
        if (p_uncertainty || q_uncertainty)
            throw std::invalid_argument("Exact tile contact requires zero occupied boundary uncertainty.");
        if (p.empty() || q.empty()) return false;
        if (p.exact && q.exact) {
            TemporaryPoints scratch(*this,8192);
            const auto common=boolean_exact_section_boundaries(*p.exact,*q.exact,
                ExactBoundaryOperation::Intersection,ExactBoundaryOrientation::Opposed,exact_budget(),[&] { tick(0); });
            tick(common.work);
            return !common.boundary.empty();
        }
        const auto& exact=p.exact ? *p.exact : *q.exact;
        const auto& native=p.exact ? q.native : p.native;
        const auto plane=exact.bounds()->min[axis];
        for (const auto& interval : native) {
            const auto common=exact_boundary_overlaps_native_interval(exact,axis,plane,interval.first,interval.second,
                exact_budget(),[&] { tick(0); });
            tick(common.work);
            if (!common.empty) return true;
        }
        return false;
    }
    std::optional<ExPolygons> split_point_contacts(const ExPolygon& polygon, const LocalSectionTile& tile,
        const char* stage, coord_t bottom, coord_t top)
    {
        std::map<std::pair<coord_t,coord_t>, std::size_t> positions;
        bool repeated = false, bounded = true;
        const auto inspect = [&](const Polygon& ring) {
            for (const auto& p : ring.points) {
                if (checkpoint) checkpoint();
                if (p.x() < tile.minimum[0] || p.x() > tile.maximum[0] || p.y() < tile.minimum[1] || p.y() > tile.maximum[1])
                    bounded = false;
                repeated = !positions.emplace(std::make_pair(p.x(),p.y()), positions.size()).second || repeated;
            }
        };
        // Bound both predicate arithmetic and detection storage before building
        // the exact cycle graph. Ordinary no-contact rings take no quadratic path.
        if (polygon_points(polygon) > (streaming ? streaming->max_live_points : limits.max_boundary_points))
            throw ResourceLimit("Local topology contact input point budget exceeded.");
        inspect(polygon.contour);
        for (const auto& hole : polygon.holes) inspect(hole);
        if (!bounded || !repeated) return {};
        const auto count = polygon_points(polygon);
        tick(count);
        reserve_live(count); // Output has exactly the original directed edges.
        const auto decompose = [&](const Polygon& ring) {
            if (ring.points.size() < 3) throw std::invalid_argument("Invalid point-contact ring.");
            TemporaryPoints stack_storage(*this, ring.points.size());
            positions.clear();
            Points stack;
            Polygons cycles;
            std::size_t cycle_points = 0;
            for (std::size_t i = 0; i <= ring.points.size(); ++i) {
                tick();
                const auto& p = ring.points[i % ring.points.size()];
                const auto key = std::make_pair(p.x(),p.y());
                const auto found = positions.find(key);
                if (found == positions.end()) {
                    positions.emplace(key, stack.size());
                    stack.push_back(p);
                    continue;
                }
                const auto first = found->second;
                if (stack.size() - first < 3)
                    throw std::invalid_argument("Unsupported zero-edge point-contact contour.");
                Points cycle(stack.begin() + first, stack.end());
                cycle_points += cycle.size();
                cycles.emplace_back(std::move(cycle));
                for (std::size_t j = first + 1; j < stack.size(); ++j)
                    positions.erase({stack[j].x(),stack[j].y()});
                stack.resize(first + 1);
            }
            if (stack.size() != 1 || cycle_points != ring.points.size())
                throw std::invalid_argument("Incomplete point-contact edge decomposition.");
            return cycles;
        };
        const auto validate_ring = [&](const Polygon& ring, const char* kind, std::size_t ring_id) {
            if (ring.points.size() < 3) throw std::invalid_argument("Invalid point-contact cycle.");
            for (const auto& p : ring.points)
                if (p.x() < tile.minimum[0] || p.x() > tile.maximum[0] || p.y() < tile.minimum[1] || p.y() > tile.maximum[1])
                    throw std::invalid_argument("Invalid point-contact cycle coordinates.");
            ExactInteger area = 0;
            for (std::size_t i = 0; i < ring.points.size(); ++i) {
                const auto& a = ring.points[i];
                const auto& b = ring.points[(i+1) % ring.points.size()];
                if (a.x() < tile.minimum[0] || a.x() > tile.maximum[0] || a.y() < tile.minimum[1] || a.y() > tile.maximum[1] || a == b)
                    throw std::invalid_argument("Invalid point-contact cycle coordinates.");
                area += ExactInteger(a.x()) * b.y() - ExactInteger(a.y()) * b.x();
                for (std::size_t j = 0; j < i; ++j) {
                    tick();
                    const auto& c = ring.points[j];
                    const auto& d = ring.points[(j+1) % ring.points.size()];
                    if (j+1 == i || (j == 0 && i+1 == ring.points.size())) {
                        // Adjacent edges may meet only at their common endpoint.
                        const Point common = j+1 == i ? a : b;
                        const Point other_a = a == common ? b : a, other_b = c == common ? d : c;
                        if (on_segment(other_a,c,d) || on_segment(other_b,a,b))
                            throw std::invalid_argument("Overlapping point-contact cycle edges.");
                    } else if (intersect_segments(a,b,c,d)) {
                        std::ostringstream message;
                        message << "Crossing point-contact cycle edges: stage=" << stage
                                << " tile=" << tile.index[0] << ',' << tile.index[1]
                                << " Znm=" << bottom << ".." << top << " ring=" << kind << ':' << ring_id
                                << " points=" << ring.points.size() << " edges=" << i << ',' << j
                                << " A=" << a.x() << ',' << a.y() << " B=" << b.x() << ',' << b.y()
                                << " C=" << c.x() << ',' << c.y() << " D=" << d.x() << ',' << d.y()
                                << " endpoint_on_edge=" << (on_segment(a,c,d) || on_segment(b,c,d) ||
                                    on_segment(c,a,b) || on_segment(d,a,b));
                        throw std::invalid_argument(message.str());
                    }
                }
            }
            if (area == 0) throw std::invalid_argument("Zero-area point-contact cycle.");
            return area;
        };
        ExPolygons result;
        Polygons holes;
        auto outer_cycles = decompose(polygon.contour);
        for (std::size_t i = 0; i < outer_cycles.size(); ++i) {
            auto& ring = outer_cycles[i];
            if (validate_ring(ring, "cycle", i) > 0) result.emplace_back(std::move(ring));
            else holes.push_back(std::move(ring));
        }
        for (std::size_t i = 0; i < polygon.holes.size(); ++i) {
            auto cycles = decompose(polygon.holes[i]);
            for (auto& hole : cycles) {
                if (validate_ring(hole, "hole", i) >= 0)
                    throw std::invalid_argument("Unsupported positive island in point-contact hole.");
                holes.push_back(std::move(hole));
            }
        }
        if (result.empty()) throw std::invalid_argument("Point-contact contour has no positive interior.");
        const auto disjoint_boundaries = [&](const Polygon& a, const Polygon& b, bool endpoint_touch) {
            std::set<std::pair<coord_t,coord_t>> contacts;
            for (std::size_t i = 0; i < a.points.size(); ++i)
                for (std::size_t j = 0; j < b.points.size(); ++j) {
                    tick();
                    const auto& p = a.points[i]; const auto& q = a.points[(i+1) % a.points.size()];
                    const auto& r = b.points[j]; const auto& s = b.points[(j+1) % b.points.size()];
                    if (!intersect_segments(p,q,r,s)) continue;
                    const bool shared = p == r || p == s || q == r || q == s;
                    if (!endpoint_touch || !shared || (p != r && p != s && on_segment(p,r,s)) ||
                        (q != r && q != s && on_segment(q,r,s)) || (r != p && r != q && on_segment(r,p,q)) ||
                        (s != p && s != q && on_segment(s,p,q)) || (p == r && q == s) || (p == s && q == r))
                        throw std::invalid_argument("Unsupported crossing or touching point-contact cycles.");
                    if (p == r || p == s) contacts.emplace(p.x(),p.y());
                    if (q == r || q == s) contacts.emplace(q.x(),q.y());
                }
            return contacts.size();
        };
        for (std::size_t i = 0; i < result.size(); ++i)
            for (std::size_t j = 0; j < i; ++j) {
                const auto& a = result[i].contour; const auto& b = result[j].contour;
                disjoint_boundaries(a,b,true);
                for (const auto& p : a.points) {
                    tick(b.points.size());
                    if (exact_location(p,b) > 0) throw std::invalid_argument("Overlapping point-contact cycle interiors.");
                }
                for (const auto& p : b.points) {
                    tick(a.points.size());
                    if (exact_location(p,a) > 0) throw std::invalid_argument("Overlapping point-contact cycle interiors.");
                }
            }
        std::vector<std::size_t> hole_owners, void_parent(holes.size()), outer_contacts(holes.size(), 0);
        for (std::size_t i = 0; i < holes.size(); ++i) void_parent[i] = i;
        const auto void_root = [&](std::size_t i) {
            while (i != void_parent[i]) { tick(); i = void_parent[i]; }
            return i;
        };
        for (std::size_t i = 0; i < holes.size(); ++i) {
            const auto& hole = holes[i];
            std::optional<std::size_t> owner;
            for (std::size_t j = 0; j < result.size(); ++j) {
                const auto contacts = disjoint_boundaries(hole,result[j].contour,true);
                if (contacts > 1) throw std::invalid_argument("Multiple boundary contacts of a point-contact void.");
                int location = -1;
                for (const auto& p : hole.points) {
                    tick(result[j].contour.points.size());
                    location = exact_location(p,result[j].contour);
                    if (location != -1) break;
                }
                if (location > 0) {
                    if (owner) throw std::invalid_argument("Ambiguous point-contact hole ownership.");
                    owner = j;
                    outer_contacts[i] = contacts;
                } else if (contacts || location == -1)
                    throw std::invalid_argument("Uncontained point-contact void.");
            }
            if (!owner) throw std::invalid_argument("Unowned point-contact hole.");
            for (std::size_t j = 0; j < i; ++j) {
                const auto contacts = disjoint_boundaries(hole,holes[j],true);
                const auto nonboundary_location = [&](const Polygon& a, const Polygon& b) {
                    for (const auto& p : a.points) {
                        tick(b.points.size());
                        const int location = exact_location(p,b);
                        if (location != -1) return location;
                    }
                    return -1;
                };
                if (nonboundary_location(hole,holes[j]) != 0 || nonboundary_location(holes[j],hole) != 0)
                    throw std::invalid_argument("Nested point-contact holes.");
                if (contacts) {
                    const auto a = void_root(i), b = void_root(j);
                    // A void chain joining two outer points or closing a loop
                    // can disconnect solid. Do not hide that behind one island.
                    if (contacts != 1 || *owner != hole_owners[j] || a == b || outer_contacts[a] + outer_contacts[b] > 1)
                        throw std::invalid_argument("Point-contact void graph may disconnect solid.");
                    void_parent[b] = a;
                    outer_contacts[a] += outer_contacts[b];
                }
            }
            hole_owners.push_back(*owner);
        }
        for (std::size_t i = 0; i < holes.size(); ++i)
            result[hole_owners[i]].holes.push_back(std::move(holes[i]));
        return result;
    }
    void split_point_contacts(ExPolygons& polygons, const LocalSectionTile& tile,
        const char* stage, coord_t bottom, coord_t top)
    {
        ExPolygons normalized;
        for (auto& polygon : polygons) {
            auto split = split_point_contacts(polygon, tile, stage, bottom, top);
            if (!split) normalized.push_back(std::move(polygon));
            else {
                live_points -= polygon_points(polygon);
                for (auto& part : *split) normalized.push_back(std::move(part));
            }
        }
        polygons = std::move(normalized);
    }
    std::size_t root(std::size_t i)
    {
        while (labels[i].parent != i) {
            labels[i].parent = labels[labels[i].parent].parent;
            i                 = labels[i].parent;
        }
        return i;
    }
    void unite(std::size_t a, std::size_t b)
    {
        a = root(a);
        b = root(b);
        if (a != b) {
            labels[b].parent = a;
            labels[a].bed    = labels[a].bed || labels[b].bed;
            if (labels[b].geometric_islands) {
                for (unsigned axis = 0; axis < 3; ++axis) {
                    labels[a].minimum[axis] = labels[a].geometric_islands ?
                        std::min(labels[a].minimum[axis], labels[b].minimum[axis]) : labels[b].minimum[axis];
                    labels[a].maximum[axis] = labels[a].geometric_islands ?
                        std::max(labels[a].maximum[axis], labels[b].maximum[axis]) : labels[b].maximum[axis];
                }
            }
            labels[a].volume_m3 += labels[b].volume_m3;
            labels[a].pieces += labels[b].pieces;
            labels[a].members += labels[b].members;
            labels[a].geometric_islands += labels[b].geometric_islands;
        }
    }
    void record_geometry(const Island& island)
    {
        auto& label = labels[root(island.label)];
        const std::array<coord_t, 3> lo{island.bounds.min.x(), island.bounds.min.y(), island.bottom};
        const std::array<coord_t, 3> hi{island.bounds.max.x(), island.bounds.max.y(), island.top};
        for (unsigned axis = 0; axis < 3; ++axis) {
            label.minimum[axis] = label.geometric_islands ? std::min(label.minimum[axis], lo[axis]) : lo[axis];
            label.maximum[axis] = label.geometric_islands ? std::max(label.maximum[axis], hi[axis]) : hi[axis];
        }
        if (island.exact) {
            TemporaryPoints scratch(*this,4096);
            const auto moments=exact_section_volume(*island.exact,lo,
                {hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]},island.bottom,island.top,
                exact_budget(),[&] { tick(0); });
            tick(moments.work);
            label.volume_m3+=moments.value;
            ++label.geometric_islands;
            return;
        }
        ExactInteger twice_area = 0;
        const auto area = [&](const Polygon& ring) {
            tick(ring.points.size());
            const auto& origin = ring.points.front();
            for (std::size_t i = 0; i < ring.points.size(); ++i) {
                if (checkpoint && i % 256 == 0) checkpoint();
                twice_area += turn(origin, ring.points[i], ring.points[(i + 1) % ring.points.size()]);
            }
        };
        area(island.footprint.contour);
        for (const auto& hole : island.footprint.holes) area(hole);
        label.volume_m3 += .5 * twice_area.convert_to<double>() * double(island.top - island.bottom) * 1e-27;
        ++label.geometric_islands;
    }
    void record_statistics()
    {
        if (statistics) {
            statistics->islands = statistics->face_records = active_records;
            statistics->peak_face_records = std::max(statistics->peak_face_records, active_records);
            statistics->trace_intervals = trace_intervals;
        }
    }
    void release_face(std::size_t id, unsigned axis)
    {
        tick();
        auto& record = records[id];
        auto& trace = record.traces[2 * axis + 1];
        trace_intervals -= trace.size();
        discard_trace(trace);
        if (record.traces[1].empty() && record.traces[3].empty()) {
            free_records.push_back(id);
            --active_records;
        }
        record_statistics();
    }
    void join_face(const ContactRecord& record, Island& island, unsigned axis)
    {
        tick();
        if (record.top <= island.bottom || island.top <= record.bottom) return;
        const auto& p = record.traces[2 * axis + 1];
        const auto& q = island.traces[2 * axis];
        tick(p.size()+q.size());
        const bool contact = contact_traces(p,q,axis,record.uncertainty,island.uncertainty);
        if (contact) unite(island.label, record.label);
    }
    void ensure_label(Island& island)
    {
        if (island.label != std::numeric_limits<std::size_t>::max()) return;
        if (streaming && labels.size() == streaming->max_islands)
            throw ResourceLimit("Local topology label budget exceeded.");
        island.label = labels.size();
        labels.push_back({island.label, island.bottom == bed});
        if (statistics) statistics->labels = labels.size();
    }
    void begin_tile(const LocalSectionTile& tile)
    {
        if (tile_open || !tile.slabs.empty())
            throw std::invalid_argument("Local topology requires a closed tile and empty packet header.");
        if (order == LocalTileOrder::Lexicographic && last_tile && !( *last_tile < tile.index))
            throw std::invalid_argument("Local topology tiles are not in strict lexicographic order.");
        if (tile_islands.size() == 4096)
            throw ResourceLimit("Local topology tile budget exceeded.");
        if (tile_islands.empty())
            initialize_grid(tile);
        tick();
        for (unsigned axis = 0; axis < 2; ++axis)
            if (std::abs(double(tile.minimum[axis])) > 1e10 || std::abs(double(tile.maximum[axis])) > 1e10 ||
                tile.maximum[axis] - tile.minimum[axis] != size[axis] ||
                static_cast<long double>(tile.minimum[axis]) - static_cast<long double>(tile.index[axis]) * size[axis] != origin[axis])
                throw std::invalid_argument("Local topology tiles must share one nonoverlapping grid.");
        const auto inserted = tile_islands.emplace(tile.index, std::vector<std::size_t>{});
        if (!inserted.second)
            throw std::invalid_argument("Duplicate local topology tile.");
        if (order == LocalTileOrder::Lexicographic) {
            while (!pending.empty() && pending.begin()->first <= tile.index) {
                const auto it = pending.begin();
                for (const auto& [id, axis] : it->second) {
                    if (it->first == tile.index) incoming[axis].push_back(id);
                    else release_face(id, axis); // Missing neighbor proved by ordered traversal.
                }
                pending.erase(it);
            }
            last_tile = tile.index;
        }
        active_tile = tile;
        previous_top = bed;
        tile_open = true;
    }
    void finish_band_frontier(const LocalSectionTile& tile,coord_t bottom,coord_t top,unsigned uncertainty,
        std::vector<Island>& current,std::size_t& temporary_traces)
    {
        auto& tile_ids=tile_islands.at(tile.index);
        // Allocate labels only for genuinely new vertical components.
        // Geometry lives for one further band; historical contact records
        // contain only boundary traces and stable source-label references.
        std::set<std::size_t> extended_records;
        for (auto& island : current) {
            tick();
            ensure_label(island);
            island.label = root(island.label);
            if (order == LocalTileOrder::Lexicographic) {
                for (unsigned axis = 0; axis < 2; ++axis) {
                    auto& first = incoming_first[axis];
                    const auto& ids = incoming[axis];
                    while (first < ids.size() && records[ids[first]].top <= bottom) { tick(); ++first; }
                    for (std::size_t j = first; j < ids.size(); ++j) {
                        if (records[ids[j]].bottom >= top) break;
                        join_face(records[ids[j]], island, axis);
                    }
                    temporary_traces -= island.traces[2 * axis].size();
                    discard_trace(island.traces[2 * axis]);
                    if (tile.index[axis] == std::numeric_limits<int>::max()) {
                        temporary_traces -= island.traces[2 * axis + 1].size();
                        discard_trace(island.traces[2 * axis + 1]);
                    }
                }
                island.label = root(island.label);
            }
            record_geometry(island);
            std::size_t count = 0;
            for (const auto& trace : island.traces) count += trace.size();
            if (!count) continue;
            for (const auto& old : previous) {
                tick();
                if (old.record == std::numeric_limits<std::size_t>::max() || extended_records.count(old.record)) continue;
                auto& record = records[old.record];
                if (record.top != bottom || record.uncertainty != uncertainty || root(record.label) != island.label) continue;
                std::size_t comparisons = count;
                for (const auto& trace : record.traces) comparisons += trace.size();
                tick(comparisons);
                if (record.traces == island.traces) {
                    record.top = top;
                    island.record = old.record;
                    extended_records.insert(old.record);
                    break;
                }
            }
            if (island.record == std::numeric_limits<std::size_t>::max()) {
                if (streaming && active_records == streaming->max_islands)
                    throw ResourceLimit("Local topology face-record budget exceeded.");
                ContactRecord record{bottom,top,island.label,std::move(island.traces),uncertainty};
                if (free_records.empty()) {
                    island.record = records.size();
                    records.push_back(std::move(record));
                } else {
                    island.record = free_records.back();
                    free_records.pop_back();
                    records[island.record] = std::move(record);
                }
                ++active_records;
                if (order == LocalTileOrder::Arbitrary) tile_ids.push_back(island.record);
                else for (unsigned axis = 0; axis < 2; ++axis) {
                    if (records[island.record].traces[2 * axis + 1].empty()) continue;
                    Index target = tile.index;
                    ++target[axis]; // An INT_MAX face was removed above.
                    pending[target].emplace_back(island.record, axis);
                }
                trace_intervals += count;
            }
            temporary_traces -= count;
            record_statistics();
            // A merged record owns the same traces already; don't retain
            // a second vector in the next band's geometry frontier.
            island.traces = {};
        }
    }
    void process_exact_bands(const LocalSectionTile& tile,bool single_band)
    {
        if (!tile_open || tile.index!=active_tile.index || tile.minimum!=active_tile.minimum ||
            tile.maximum!=active_tile.maximum || (single_band && tile.slabs.empty()))
            throw std::invalid_argument("Local topology packet does not match its open tile.");
        std::map<Interval,std::map<std::size_t,std::vector<const DepositedSlab*>>> groups;
        std::size_t borrowed_points=0;
        std::optional<Interval> packet_interval;
        for (const auto& slab : tile.slabs) {
            tick();
            if (slab.boundary_uncertainty_nm>7)
                throw std::invalid_argument("Invalid local boundary uncertainty.");
            const Interval interval{native_height(slab.bottom_m),native_height(slab.top_m)};
            if (interval.first<bed || interval.second<=interval.first ||
                (single_band && packet_interval && *packet_interval!=interval))
                throw std::invalid_argument("Invalid exact local topology interval.");
            packet_interval=interval;
            if (++slabs>(streaming ? streaming->max_slabs : limits.max_slabs))
                throw ResourceLimit("Local topology slab budget exceeded.");
            if (statistics) statistics->slabs=slabs;
            const auto count=slab.footprint.vertex_count();
            charge(input_points,count);
            borrowed_points+=count;
            groups[interval][slab.material].push_back(&slab);
        }
        TemporaryPoints borrowed(*this,borrowed_points);
        const auto rectangle=canonical(SectionFootprint(ExPolygons{ExPolygon(Polygon{
            {tile.minimum[0],tile.minimum[1]},{tile.maximum[0],tile.minimum[1]},
            {tile.maximum[0],tile.maximum[1]},{tile.minimum[0],tile.maximum[1]}})}));
        TemporaryPoints rectangle_storage(*this,rectangle.vertex_count());
        for (const auto& [interval,sources] : groups) {
            const auto [bottom,top]=interval;
            if (bottom<previous_top) throw std::invalid_argument("Overlapping local topology Z intervals.");
            unsigned material_uncertainty=0;
            bool has_partition_evidence=false;
            std::vector<const DepositedSlab*> partition;
            for (const auto& entry : sources) for (const auto* slab : entry.second) {
                tick();
                material_uncertainty=std::max(material_uncertainty,slab->boundary_uncertainty_nm);
                has_partition_evidence=has_partition_evidence || bool(slab->occupied_partition);
                partition.push_back(slab);
            }
            if (material_uncertainty || has_partition_evidence) {
                // Only the COMPLETE unchanged authenticated partition proves
                // its occupied union fixed. Its material sensitivity remains
                // on emitted pieces, but does not perturb occupied contacts.
                TemporaryPoints scratch(*this,4096);
                tick(qualify_complete_occupied_partition(partition,exact_budget(),[&] { tick(0); }));
            }
            std::map<std::size_t,std::vector<ExactSection>> materials;
            std::size_t material_points=0;
            ExactSection occupied;
            TemporaryPoints occupied_storage(*this,0);
            for (const auto& [material,parts] : sources) {
                ExactSection joined;
                TemporaryPoints joined_storage(*this,0);
                for (const auto* slab : parts) {
                    const auto section=canonical(slab->footprint);
                    TemporaryPoints section_storage(*this,section.vertex_count());
                    validate_canonical(section);
                    if (!combine(section,rectangle,ExactSectionOperation::Difference).empty())
                        throw std::invalid_argument("Exact local topology footprint escapes tile.");
                    auto merged=combine(joined,section,ExactSectionOperation::Union);
                    reserve_live(merged.vertex_count());
                    live_points-=joined_storage.count;
                    joined_storage.count=merged.vertex_count();
                    joined=std::move(merged);
                }
                validate_canonical(joined);
                const auto common=combine(occupied,joined,ExactSectionOperation::Intersection);
                TemporaryPoints common_storage(*this,common.vertex_count());
                if (!common.empty()) {
                    const auto box=*common.bounds();
                    TemporaryPoints scratch(*this,4096);
                    const auto measure=exact_section_volume(common,{box.min.x(),box.min.y(),0},
                        {box.max.x()-box.min.x(),box.max.y()-box.min.y(),1000000000},0,1000000000,
                        exact_budget(),[&] { tick(0); });
                    tick(measure.work);
                    throw LayerMaterialOverlap(double(bottom)*1e-9,double(top)*1e-9,measure.value,material);
                }
                auto merged=combine(occupied,joined,ExactSectionOperation::Union);
                reserve_live(merged.vertex_count());
                live_points-=occupied_storage.count;
                occupied_storage.count=merged.vertex_count();
                occupied=std::move(merged);
                TemporaryPoints scratch(*this,limits.exact_interior_components ? 16384 : 4096);
                auto components=limits.exact_interior_components ?
                    exact_section_interior_components(joined,exact_budget(),[&] { tick(0); }) :
                    exact_section_components(joined,exact_budget(),[&] { tick(0); });
                tick(components.work);
                for (const auto& component : components.sections) {
                    charge(retained_points,component.vertex_count());
                    reserve_live(component.vertex_count());
                    material_points+=component.vertex_count();
                }
                materials.emplace(material,std::move(components.sections));
            }
            validate_canonical(occupied);
            ExactSectionComponents components;
            {
                TemporaryPoints scratch(*this,limits.exact_interior_components ? 16384 : 4096);
                components=limits.exact_interior_components ?
                    exact_section_interior_components(occupied,exact_budget(),[&] { tick(0); }) :
                    exact_section_components(occupied,exact_budget(),[&] { tick(0); });
                tick(components.work);
                for (const auto& component : components.sections) reserve_live(component.vertex_count());
            }
            std::vector<Island> current;
            std::size_t temporary_traces=0;
            for (auto& component : components.sections) {
                if (streaming && current.size()==streaming->max_islands)
                    throw ResourceLimit("Local topology island budget exceeded.");
                charge(retained_points,component.vertex_count());
                Island island{{},*component.bounds(),bottom,top};
                island.exact=std::move(component);
                TemporaryPoints scratch(*this,16384);
                const auto boundary=extract_exact_section_boundary(*island.exact,exact_budget(),[&] { tick(0); });
                tick(boundary.work);
                for (unsigned axis=0;axis<2;++axis) for (unsigned side=0;side<2;++side) {
                    const auto trace=exact_boundary_on_native_plane(boundary.boundary,axis,
                        side ? tile.maximum[axis] : tile.minimum[axis],exact_budget(),[&] { tick(0); });
                    tick(trace.work);
                    const auto count=trace.boundary.segment_count();
                    const auto ceiling=streaming ? streaming->max_trace_intervals : limits.max_boundary_points;
                    if (count>ceiling-trace_intervals-temporary_traces)
                        throw ResourceLimit("Local topology trace budget exceeded.");
                    reserve_live(2*count);
                    island.traces[2*axis+side]=Trace(trace.boundary);
                    temporary_traces+=count;
                }
                current.push_back(std::move(island));
                if (statistics) ++statistics->processed_islands;
            }
            if (bottom==previous_top) for (const auto& first : previous) for (auto& second : current) {
                tick();
                if (!first.bounds.overlap(second.bounds)) continue;
                if (first.uncertainty)
                    throw std::invalid_argument("Exact vertical contact requires zero occupied boundary uncertainty.");
                const auto first_section=first.exact ? *first.exact : canonical(SectionFootprint(ExPolygons{first.footprint}));
                TemporaryPoints first_storage(*this,first_section.vertex_count());
                if (combine(first_section,*second.exact,ExactSectionOperation::Intersection).empty()) continue;
                if (second.label==std::numeric_limits<std::size_t>::max()) second.label=root(first.label);
                else unite(first.label,second.label);
            }
            finish_band_frontier(tile,bottom,top,0,current,temporary_traces);
            for (auto& [material,parts] : materials) for (auto& part : parts) {
                std::optional<std::size_t> owner;
                for (std::size_t i=0;i<current.size();++i) {
                    tick();
                    if (!part.bounds()->overlap(current[i].bounds)) continue;
                    if (combine(part,*current[i].exact,ExactSectionOperation::Difference).empty()) {
                        if (owner) throw std::invalid_argument("Ambiguous exact local material connectivity.");
                        owner=i;
                    }
                }
                if (!owner) throw std::invalid_argument("Unresolved exact local material connectivity.");
                const auto label=root(current[*owner].label);
                // Component splitting invalidates the complete-partition
                // evidence: deliberately do not copy its token to this piece.
                const SectionPiece piece{{double(bottom)*1e-9,double(top)*1e-9,material,
                    SectionFootprint(std::move(part)),material_uncertainty},label};
                ++emitted_pieces;
                ++labels[label].pieces;
                piece_callback(piece);
                if (state!=State::Active) throw std::invalid_argument("Local topology callback reentered its stream.");
            }
            live_points-=material_points;
            for (const auto& island : previous) live_points-=island_points(island);
            previous=std::move(current);
            previous_top=top;
        }
    }
    void process_bands(const LocalSectionTile& tile, bool single_band)
    {
        if (std::any_of(tile.slabs.begin(),tile.slabs.end(),[](const DepositedSlab& slab) { return slab.footprint.is_exact(); }) ||
            std::any_of(previous.begin(),previous.end(),[](const Island& island) { return bool(island.exact); })) {
            process_exact_bands(tile,single_band);
            return;
        }
        if (!tile_open || tile.index != active_tile.index || tile.minimum != active_tile.minimum ||
            tile.maximum != active_tile.maximum || (single_band && tile.slabs.empty()))
            throw std::invalid_argument("Local topology packet does not match its open tile.");
        std::size_t borrowed_points = 0;
        std::optional<Interval> band;
        for (const auto& slab : tile.slabs) {
            if (checkpoint) checkpoint();
            const Interval interval{native_height(slab.bottom_m), native_height(slab.top_m)};
            if (single_band && band && interval != *band)
                throw std::invalid_argument("Local topology packet must contain one complete Z interval.");
            band = interval;
            const auto n = polygon_points(slab.footprint.native());
            if (streaming && n > streaming->max_live_points - borrowed_points)
                throw ResourceLimit("Local topology borrowed tile point budget exceeded.");
            borrowed_points += n;
        }
        TemporaryPoints borrowed(*this, borrowed_points);
        std::map<Interval, std::map<std::size_t, ExPolygons>> groups;
        std::map<Interval, unsigned> uncertainties;
        for (const auto& slab : tile.slabs) {
            tick();
            if (++slabs > (streaming ? streaming->max_slabs : limits.max_slabs))
                throw ResourceLimit("Local topology slab budget exceeded.");
            if (statistics) statistics->slabs = slabs;
            const coord_t bottom = native_height(slab.bottom_m), top = native_height(slab.top_m);
            if (bottom < bed || top <= bottom)
                throw std::invalid_argument("Invalid local topology interval.");
            if (slab.boundary_uncertainty_nm > 7)
                throw std::invalid_argument("Invalid local boundary uncertainty.");
            auto& uncertainty = uncertainties[{bottom, top}];
            uncertainty       = std::max(uncertainty, slab.boundary_uncertainty_nm);
            auto& polygons    = groups[{bottom, top}][slab.material];
            for (const auto& polygon : slab.footprint.native()) {
                auto split = split_point_contacts(polygon, tile, "input", bottom, top);
                if (split) {
                    for (auto& part : *split) {
                        charge(input_points, validate_local_polygon(part, tile, "split input", bottom, top));
                        polygons.push_back(std::move(part));
                    }
                } else {
                    reserve_live(polygon_points(polygon));
                    charge(input_points, validate_local_polygon(polygon, tile, "input", bottom, top));
                    polygons.push_back(polygon);
                }
            }
        }
        for (auto& [interval, materials] : groups) {
            tick();
            const auto [bottom, top]   = interval;
            const unsigned uncertainty = uncertainties.at(interval);
            if (bottom < previous_top)
                throw std::invalid_argument("Overlapping local topology Z intervals.");
            ExPolygons occupied;
            for (auto& [material, polygons] : materials) {
                tick();
                // A single validated simple polygon already is its union.
                // Avoid re-clipping its certified interface boundary chain.
                if (polygons.size()>1) {
                    auto normalized = union_ex(polygons);
                    reserve_live(polygon_points(normalized));
                    live_points -= polygon_points(polygons);
                    polygons = std::move(normalized);
                }
                split_point_contacts(polygons, tile, "material union", bottom, top);
                const auto overlap = intersection_ex(occupied, polygons);
                TemporaryPoints overlap_storage(*this, polygon_points(overlap));
                if (!overlap.empty()) {
                    // Integer clipping can round a boundary-only contact into
                    // a positive-area sliver. Reject actual overlap, but do
                    // not infer it from a rounded diagnostic polygon alone.
                    const auto available=(streaming ? streaming->max_work : limits.max_component_pairs)-work;
                    if (!available) throw ResourceLimit("Local topology exact overlap budget exhausted.");
                    const auto exact=exact_section_empty(occupied,polygons,ExactSectionOperation::Intersection,
                        1024,std::min<std::size_t>(available,2000000),checkpoint);
                    tick(exact.work);
                    overlap_work += exact.work;
                    if (!exact.empty) {
                        std::size_t diagnostic_points=0;
                        for (const auto& slab : tile.slabs)
                            if (native_height(slab.bottom_m)==bottom && native_height(slab.top_m)==top && slab.material<=material)
                                diagnostic_points+=polygon_points(slab.footprint.native());
                        TemporaryPoints input_storage(*this,diagnostic_points);
                        ExPolygons input_previous,input_current;
                        for (const auto& slab : tile.slabs) {
                            if (native_height(slab.bottom_m)!=bottom || native_height(slab.top_m)!=top ||
                                slab.material>material) continue;
                            auto& target=slab.material==material ? input_current : input_previous;
                            target.insert(target.end(),slab.footprint.native().begin(),slab.footprint.native().end());
                        }
                        const auto remaining=(streaming ? streaming->max_work : limits.max_component_pairs)-work;
                        if (!remaining) throw ResourceLimit("Local topology input diagnostic budget exhausted.");
                        const auto input_exact=exact_section_empty(input_previous,input_current,ExactSectionOperation::Intersection,
                            1024,std::min<std::size_t>(remaining,2000000),checkpoint);
                        tick(input_exact.work);
                        diagnostic_work += input_exact.work;
                        double area = 0.;
                        for (const auto& polygon : overlap)
                            area += polygon.area() * 1e-18;
                        throw LayerMaterialOverlap(double(bottom) * 1e-9, double(top) * 1e-9, area, material,
                            " Exact local topology overlap after material union; tile=" +
                            std::to_string(tile.index[0]) + "," + std::to_string(tile.index[1]) +
                            "; polygons=" + std::to_string(polygons.size()) +
                            "; input_disjoint=" + (input_exact.empty ? "true" : "false") + ".");
                    }
                }
                reserve_live(polygon_points(polygons));
                occupied.insert(occupied.end(), polygons.begin(), polygons.end());
            }
            if (occupied.size() > 1) {
                auto normalized = union_ex(occupied);
                reserve_live(polygon_points(normalized));
                live_points -= polygon_points(occupied);
                occupied = std::move(normalized);
            }
            split_point_contacts(occupied, tile, "occupied union", bottom, top);
            // Establish material ownership before publishing labels, contact
            // records or pieces. A rounded union may lose a narrow input edge;
            // retry from the original material chains, never by assigning a
            // partially covered material to the nearest island.
            std::size_t material_pieces=0;
            for (const auto& item : materials) material_pieces+=item.second.size();
            TemporaryPoints assignment_storage(*this,material_pieces);
            std::vector<std::size_t> material_islands;
            material_islands.reserve(material_pieces);
            std::string assignment_failure;
            const auto assign_materials = [&, lower=bottom, upper=top, &stage_materials=materials] {
                material_islands.clear();
                for (auto& [material, polygons] : stage_materials)
                    for (auto& polygon : polygons) {
                        charge(retained_points, validate_local_polygon(polygon, tile, "material union", lower, upper));
                        bool assigned = false;
                        const auto material_bounds = get_extents(polygon);
                        std::ostringstream candidates;
                        std::size_t reported = 0;
                        bool reported_boundaries = false;
                        for (std::size_t island_index=0; island_index<occupied.size(); ++island_index) {
                            const auto& island_footprint=occupied[island_index];
                            const auto island_bounds=get_extents(island_footprint);
                            tick();
                            if (!material_bounds.overlap(island_bounds)) continue;
                            tick(polygon_points(polygon) + polygon_points(island_footprint));
                            bool contained = polygon == island_footprint;
                            if (!contained) {
                                TemporaryPoints input_copy(*this, polygon_points(polygon) + polygon_points(island_footprint));
                                const auto outside = diff_ex(ExPolygons{polygon}, ExPolygons{island_footprint});
                                TemporaryPoints output(*this, polygon_points(outside));
                                contained = outside.empty();
                            }
                            if (!contained) {
                                const auto ceiling = streaming ? streaming->max_work : limits.max_component_pairs;
                                if (work >= ceiling) throw ResourceLimit("Local material containment proof budget exhausted.");
                                // A rounded difference fragment is not evidence of
                                // positive-area escape. Check the original integer
                                // boundaries exactly; never assign by proximity or
                                // by ignoring a genuinely uncovered material piece.
                                // Both operands were validated as native oriented
                                // chains; a hole touching its outer ring retains
                                // that exact fill, not strict simple-hole semantics.
                                const auto available = std::min<std::size_t>(ceiling-work,2000000);
                                ExactSectionPredicate exact;
                                ++containment_queries;
                                try {
                                    exact = exact_section_empty({polygon}, {island_footprint}, ExactSectionOperation::Difference,
                                        1024, available, checkpoint, true);
                                } catch (const ResourceLimit& error) {
                                    throw ResourceLimit(std::string(error.what()) + " Local material containment; tile=" +
                                        std::to_string(tile.index[0]) + "," + std::to_string(tile.index[1]) +
                                        " z_nm=" + std::to_string(lower) + "," + std::to_string(upper) +
                                        " material=" + std::to_string(material) + " operand_vertices=" +
                                        std::to_string(polygon_points(polygon)) + "," + std::to_string(polygon_points(island_footprint)) +
                                        " local_available=" + std::to_string(available) + " cumulative_work=" +
                                        std::to_string(work) + " cumulative_limit=" + std::to_string(ceiling) +
                                        " containment_queries=" + std::to_string(containment_queries) +
                                        " containment_work=" + std::to_string(containment_work) +
                                        " overlap_work=" + std::to_string(overlap_work) +
                                        " diagnostic_work=" + std::to_string(diagnostic_work) +
                                        " other_work=" + std::to_string(work-containment_work-overlap_work-diagnostic_work));
                                }
                                tick(exact.work);
                                containment_work += exact.work;
                                contained = exact.empty;
                                if (!contained && reported++ < 8) {
                                    candidates << " island_bbox_nm=" << island_bounds.min.x() << ',' << island_bounds.min.y() << ','
                                               << island_bounds.max.x() << ',' << island_bounds.max.y()
                                               << " island_vertices=" << polygon_points(island_footprint) << " exact_outside=positive";
                                    // One bounded exact-input witness makes a lost
                                    // union boundary reproducible without dumping a
                                    // model or treating its rounded difference as truth.
                                    if (!reported_boundaries && polygon_points(polygon)+polygon_points(island_footprint)<=128 &&
                                        island_bounds.contains(material_bounds)) {
                                        reported_boundaries = true;
                                        const auto boundary = [&](const ExPolygon& value, const char* name) {
                                            candidates << ' ' << name << '=';
                                            const auto ring = [&](const Polygon& shape) {
                                                candidates << '[';
                                                for (const auto& p : shape.points) {
                                                    tick();
                                                    candidates << '(' << p.x() << ',' << p.y() << ')';
                                                }
                                                candidates << ']';
                                            };
                                            ring(value.contour);
                                            for (const auto& hole : value.holes) ring(hole);
                                        };
                                        boundary(polygon,"material_boundaries_nm");
                                        boundary(island_footprint,"occupied_boundaries_nm");
                                    }
                                }
                            }
                            if (contained) {
                                material_islands.push_back(island_index);
                                assigned = true;
                                break;
                            }
                        }
                        if (!assigned) {
                            assignment_failure = "Unresolved local material connectivity; tile=" + std::to_string(tile.index[0]) +
                                "," + std::to_string(tile.index[1]) + " z_nm=" + std::to_string(lower) + "," + std::to_string(upper) +
                                " material=" + std::to_string(material) + " material_bbox_nm=" + std::to_string(material_bounds.min.x()) +
                                "," + std::to_string(material_bounds.min.y()) + "," + std::to_string(material_bounds.max.x()) + "," +
                                std::to_string(material_bounds.max.y()) + " material_vertices=" + std::to_string(polygon_points(polygon)) +
                                " occupied_islands=" + std::to_string(occupied.size()) + candidates.str();
                            return false;
                        }
                    }

                return true;
            };
            if (!assign_materials()) {
                const auto ceiling = streaming ? streaming->max_work : limits.max_component_pairs;
                if (work >= ceiling) throw ResourceLimit("Exact occupied union work budget exhausted.");
                std::size_t source_points=0;
                for (const auto& item : materials) {
                    const auto n=polygon_points(item.second);
                    if (n>1024-source_points) throw ResourceLimit("Exact occupied union input budget exceeded.");
                    source_points+=n;
                }
                TemporaryPoints source_storage(*this,source_points);
                ExPolygons source;
                for (const auto& item : materials)
                    source.insert(source.end(),item.second.begin(),item.second.end());
                bool callback_failed=false;
                const auto checked_checkpoint = [&] {
                    if (checkpoint) try { checkpoint(); }
                    catch (...) { callback_failed=true; throw; }
                };
                ExactIntegerSectionUnion exact;
                const auto available=std::min<std::size_t>(ceiling-work,4000000);
                try {
                    exact=exact_integer_section_union(source,available,checked_checkpoint);
                } catch (const NonIntegralSectionUnion& error) {
                    if (callback_failed) throw;
                    tick(error.work);
                    if (error.work>=available) throw ResourceLimit("Occupied assignment diagnostic budget exhausted.");
                    // Diagnostic gate only: no represented component mapping
                    // is admitted until its contact uncertainty is accounted.
                    ExactOccupiedAssignment certificate;
                    try {
                        try {
                            certificate=certify_occupied_assignment(source,occupied,
                                BoundingBox(Point(tile.minimum[0],tile.minimum[1]),Point(tile.maximum[0],tile.maximum[1])),
                                available-error.work,checked_checkpoint);
                        } catch (const OccupiedAssignmentOverlap& overlap) {
                            if (callback_failed) throw;
                            tick(overlap.work);
                            diagnostic_work+=overlap.work;
                            std::ostringstream raw_detail;
                            if (overlap.original) {
                                const auto owner = [&stage_materials=materials](std::size_t index) -> std::optional<std::size_t> {
                                    for (const auto& item : stage_materials) {
                                        if (index<item.second.size()) return item.first;
                                        index-=item.second.size();
                                    }
                                    return {};
                                };
                                const auto first=owner(overlap.first),second=owner(overlap.second);
                                std::size_t points=0;
                                auto remaining=available-error.work-overlap.work;
                                if (first && second && first!=second && tile.slabs.size()<remaining/2) {
                                    tick(2*tile.slabs.size());
                                    remaining-=2*tile.slabs.size();
                                    for (const auto& slab : tile.slabs)
                                        if (native_height(slab.bottom_m)==bottom && native_height(slab.top_m)==top &&
                                            (slab.material==*first || slab.material==*second))
                                            points+=polygon_points(slab.footprint.native());
                                    if (points<=1024 && remaining) {
                                        TemporaryPoints raw_storage(*this,points);
                                        ExPolygons a,b;
                                        for (const auto& slab : tile.slabs) {
                                            if (native_height(slab.bottom_m)!=bottom || native_height(slab.top_m)!=top) continue;
                                            if (slab.material!=*first && slab.material!=*second) continue;
                                            auto& target=slab.material==*first ? a : b;
                                            target.insert(target.end(),slab.footprint.native().begin(),slab.footprint.native().end());
                                        }
                                        const auto raw=exact_section_empty(a,b,ExactSectionOperation::Intersection,1024,
                                            std::min<std::size_t>(remaining,2000000),checked_checkpoint,true);
                                        tick(raw.work);
                                        diagnostic_work+=raw.work;
                                        raw_detail << " raw_materials=" << *first << ',' << *second
                                                   << " raw_input_disjoint=" << (raw.empty ? "true" : "false");
                                    } else raw_detail << " raw_input_comparison=budget_limited";
                                } else raw_detail << " raw_input_comparison=same_material_or_budget_limited";
                            }
                            throw std::invalid_argument(std::string(overlap.what())+raw_detail.str());
                        }
                    } catch (const std::invalid_argument& failure) {
                        if (callback_failed) throw;
                        std::ostringstream detail;
                        detail << failure.what() << " tile=" << tile.index[0] << ',' << tile.index[1]
                               << " z_nm=" << bottom << ',' << top;
                        std::size_t first=0, reported=0;
                        for (const auto& item : materials) {
                            if (reported++==8) break;
                            detail << " material=" << item.first << " piece_range=" << first << ':' << first+item.second.size();
                            first+=item.second.size();
                        }
                        throw std::invalid_argument(detail.str());
                    }
                    tick(certificate.work);
                    throw std::invalid_argument("Exact component mapping certified; mapped admission not enabled.");
                }
                tick(exact.work);
                reserve_live(polygon_points(exact.polygons));
                live_points-=polygon_points(occupied);
                occupied=std::move(exact.polygons);
                split_point_contacts(occupied,tile,"exact occupied union",bottom,top);
                if (!assign_materials()) throw std::invalid_argument(assignment_failure);
            }
            std::vector<Island> current;
            std::size_t temporary_traces = 0;
            for (auto& polygon : occupied) {
                if (streaming && current.size() == streaming->max_islands)
                    throw ResourceLimit("Local topology island budget exceeded.");
                charge(retained_points, validate_local_polygon(polygon, tile, "occupied union", bottom, top));
                Island island{std::move(polygon), {}, bottom, top};
                island.uncertainty = uncertainty;
                if (uncertainty) {
                    tick();
                    TemporaryPoints input_copy(*this, polygon_points(island.footprint));
                    const auto inner = offset_ex(ExPolygons{island.footprint}, -float(uncertainty + 2), jtRound, .25);
                    TemporaryPoints inner_storage(*this, polygon_points(inner));
                    // This downstream gate certifies a surviving connected
                    // solid, not void topology or local stiffness. Erosion may
                    // merge/open holes without splitting the solid. Keep the
                    // ORIGINAL holes for every piece, trace and integration;
                    // source/representation hole certificates remain separate.
                    if (inner.size() != 1) {
                        const auto bounds = get_extents(island.footprint);
                        std::size_t holes = 0;
                        double inner_area = 0.;
                        for (const auto& polygon : inner) {
                            holes += polygon.holes.size();
                            inner_area += polygon.area();
                        }
                        std::ostringstream detail;
                        detail.precision(17);
                        detail << "Local occupied island is unresolved within its boundary uncertainty. tile="
                               << tile.index[0] << ',' << tile.index[1]
                               << " z_nm=" << bottom << ',' << top
                               << " uncertainty_nm=" << uncertainty << " erosion_nm=" << uncertainty + 2
                               << " bbox_nm=" << bounds.min.x() << ',' << bounds.min.y() << ','
                               << bounds.max.x() << ',' << bounds.max.y()
                               << " tile_bbox_nm=" << tile.minimum[0] << ',' << tile.minimum[1] << ','
                               << tile.maximum[0] << ',' << tile.maximum[1]
                               << " area_nm2=" << island.footprint.area()
                               << " vertices=" << polygon_points(island.footprint)
                               << " holes=" << island.footprint.holes.size()
                               << " eroded_components=" << inner.size() << " eroded_holes=" << holes
                               << " eroded_vertices=" << polygon_points(inner) << " eroded_area_nm2=" << inner_area;
                        if (island.footprint.contour.points.size() <= 16)
                            for (const auto& p : island.footprint.contour.points)
                                detail << " vertex_nm=" << p.x() << ',' << p.y();
                        const auto diagnose_stage = [&, lower = bottom, upper = top, &stage_materials = materials](bool original) {
                            const char* stage = original ? "raw_input" : "material_union";
                            ExPolygons candidates;
                            TemporaryPoints candidate_storage(*this, 0);
                            std::size_t count = 0;
                            bool oversized = false;
                            const auto append = [&](const ExPolygon& polygon) {
                                const auto vertices = polygon_points(polygon);
                                tick(vertices + 1);
                                if (!get_extents(polygon).overlap(bounds)) return;
                                if (vertices > 1024 - std::min<std::size_t>(1024, count)) oversized = true;
                                count += vertices;
                                if (!oversized) {
                                    reserve_live(vertices);
                                    candidate_storage.count += vertices;
                                    candidates.push_back(polygon);
                                }
                            };
                            if (original) {
                                for (const auto& slab : tile.slabs)
                                    if (native_height(slab.bottom_m) == lower && native_height(slab.top_m) == upper)
                                        for (const auto& polygon : slab.footprint.native()) append(polygon);
                            } else {
                                for (const auto& material : stage_materials)
                                    for (const auto& polygon : material.second) append(polygon);
                            }
                            detail << ' ' << stage << "_candidate_vertices=" << count;
                            if (oversized || count + polygon_points(island.footprint) > 1024) {
                                detail << ' ' << stage << "_exact=not_evaluated_vertex_budget";
                                return;
                            }
                            bool cancelled = false;
                            const auto check = [&] {
                                if (checkpoint) try { checkpoint(); } catch (...) { cancelled = true; throw; }
                            };
                            try {
                                for (const auto operation : {ExactSectionOperation::Intersection, ExactSectionOperation::Difference}) {
                                    const auto ceiling = streaming ? streaming->max_work : limits.max_component_pairs;
                                    if (work >= ceiling) throw ResourceLimit("Island witness diagnostic budget exhausted.");
                                    const auto exact = exact_section_empty({island.footprint}, candidates, operation, 1024,
                                        std::min<std::size_t>(ceiling - work, 2000000), check);
                                    tick(exact.work);
                                    diagnostic_work += exact.work;
                                    detail << ' ' << stage << (operation == ExactSectionOperation::Intersection ?
                                        "_exact_intersection_empty=" : "_exact_outside_empty=") << exact.empty;
                                }
                            } catch (const std::invalid_argument& error) {
                                if (cancelled) throw;
                                detail << ' ' << stage << "_exact_unsupported=" << std::string(error.what()).substr(0, 160);
                            }
                        };
                        diagnose_stage(true);
                        diagnose_stage(false);
                        throw std::invalid_argument(detail.str());
                    }
                }
                island.bounds = get_extents(island.footprint);
                for (unsigned axis = 0; axis < 2; ++axis) {
                    for (unsigned side = 0; side < 2; ++side) {
                        auto& trace = island.traces[2 * axis + side];
                        trace = boundary_trace(island.footprint, axis, side ? tile.maximum[axis] : tile.minimum[axis],
                            streaming ? streaming->max_trace_intervals - trace_intervals - temporary_traces :
                                        std::numeric_limits<std::size_t>::max());
                        temporary_traces += trace.size();
                    }
                }
                current.push_back(std::move(island));
                if (statistics) ++statistics->processed_islands;
            }
            if (bottom == previous_top)
                for (const auto& first : previous)
                    for (auto& second : current) {
                        tick();
                        bool contact = false;
                        if (first.bounds.overlap(second.bounds)) {
                            TemporaryPoints input_copy(*this, polygon_points(first.footprint) + polygon_points(second.footprint));
                            const auto overlap = intersection_ex(ExPolygons{first.footprint}, ExPolygons{second.footprint});
                            TemporaryPoints output(*this, polygon_points(overlap));
                            contact = !overlap.empty();
                        }
                        if (first.uncertainty || second.uncertainty) {
                            tick();
                            const float sign = contact ? -1.f : 1.f;
                            TemporaryPoints input_copy(*this, polygon_points(first.footprint) + polygon_points(second.footprint));
                            const auto p     = offset_ex(ExPolygons{first.footprint}, sign * (first.uncertainty + 2), jtRound, .25);
                            TemporaryPoints p_storage(*this, polygon_points(p));
                            const auto q     = offset_ex(ExPolygons{second.footprint}, sign * (second.uncertainty + 2), jtRound, .25);
                            TemporaryPoints q_storage(*this, polygon_points(q));
                            const auto count_points = [&](const ExPolygons& polygons) {
                                for (const auto& polygon : polygons) {
                                    charge(retained_points, polygon.contour.points.size());
                                    for (const auto& hole : polygon.holes)
                                        charge(retained_points, hole.points.size());
                                }
                            };
                            count_points(p);
                            count_points(q);
                            const auto expanded_overlap = intersection_ex(p, q);
                            TemporaryPoints overlap_storage(*this, polygon_points(expanded_overlap));
                            if (contact == expanded_overlap.empty()) {
                                const auto edge_limit = std::min<std::size_t>(100000,
                                    streaming ? streaming->max_islands : limits.max_component_pairs);
                                if (uncertain_contacts.size() >= edge_limit)
                                    throw ResourceLimit("Local topology uncertain contact budget exceeded.");
                                tick();
                                ensure_label(second);
                                // Diagnose the exact nominal/offset relation without
                                // substituting rounded slivers for a contact proof.
                                // This is first-observation, bounded work; it never admits
                                // a contact that the uncertainty test cannot certify.
                                if (uncertain_contacts.empty()) {
                                    std::ostringstream detail;
                                    detail.precision(17);
                                    detail << "Adjacent section contact is unresolved within boundary uncertainty. tile="
                                           << tile.index[0] << ',' << tile.index[1]
                                           << " first_z_nm=" << first.bottom << ',' << first.top
                                           << " second_z_nm=" << second.bottom << ',' << second.top
                                           << " uncertainty_nm=" << first.uncertainty << ',' << second.uncertainty
                                           << " nominal_contact=" << contact
                                           << " tested_offset_nm=" << sign * (first.uncertainty + 2) << ',' << sign * (second.uncertainty + 2)
                                           << " offset_contact=" << !expanded_overlap.empty()
                                           << " first_bbox_nm=" << first.bounds.min.x() << ',' << first.bounds.min.y() << ','
                                           << first.bounds.max.x() << ',' << first.bounds.max.y()
                                           << " second_bbox_nm=" << second.bounds.min.x() << ',' << second.bounds.min.y() << ','
                                           << second.bounds.max.x() << ',' << second.bounds.max.y()
                                           << " tile_bbox_nm=" << tile.minimum[0] << ',' << tile.minimum[1] << ','
                                           << tile.maximum[0] << ',' << tile.maximum[1]
                                           << " points=" << polygon_points(first.footprint) << ',' << polygon_points(second.footprint);
                                    const auto exact_contact = [&](const ExPolygons& a, const ExPolygons& b, const char* label) {
                                        const auto ceiling = streaming ? streaming->max_work : limits.max_component_pairs;
                                        const auto remaining = ceiling > work ? ceiling - work : 0;
                                        if (polygon_points(a) + polygon_points(b) > 1024 || remaining < 1024) {
                                            detail << ' ' << label << "=not_evaluated_budget";
                                            return;
                                        }
                                        const auto exact = exact_section_empty(a, b, ExactSectionOperation::Intersection,
                                            1024, std::min<std::size_t>(remaining, 2000000), checkpoint);
                                        tick(exact.work);
                                        diagnostic_work += exact.work;
                                        detail << ' ' << label << '=' << !exact.empty;
                                    };
                                    exact_contact({first.footprint}, {second.footprint}, "exact_nominal_contact");
                                    exact_contact(p, q, "exact_offset_contact");
                                    double offset_area = 0.;
                                    for (const auto& polygon : expanded_overlap) offset_area += polygon.area();
                                    detail << " offset_overlap_area_nm2=" << offset_area;
                                    first_uncertain_contact = detail.str();
                                }
                                uncertain_contacts.push_back({first.label, second.label, tile.index, bottom});
                                if (statistics) statistics->uncertain_contacts = uncertain_contacts.size();
                                // A possible edge must never establish connectivity.
                                // Finish verifies whether independent certified paths
                                // make this edge irrelevant to the component partition.
                                contact = false;
                            }
                        }
                        if (contact) {
                            if (second.label == std::numeric_limits<std::size_t>::max()) second.label = root(first.label);
                            else unite(first.label, second.label);
                        }
                    }
            finish_band_frontier(tile,bottom,top,uncertainty,current,temporary_traces);
            std::size_t material_index=0;
            for (auto& [material, polygons] : materials)
                for (auto& polygon : polygons) {
                    const auto& island=current.at(material_islands.at(material_index++));
                    ExPolygons footprint;
                    footprint.push_back(std::move(polygon));
                    const SectionPiece piece{{double(bottom)*1e-9,double(top)*1e-9,material,std::move(footprint),uncertainty},
                                             root(island.label)};
                    ++emitted_pieces;
                    ++labels[root(island.label)].pieces;
                    piece_callback(piece);
                    live_points-=polygon_points(piece.slab.footprint.native());
                    if (state!=State::Active)
                        throw std::invalid_argument("Local topology callback reentered its stream.");
                }
            // Full island polygons are needed only for contact with the next
            // vertical band. Later tile joins use the retained native face
            // traces, while material pieces already own their integration
            // geometry. Do not retain a second copy for every historical band.
            for (const auto& island : previous) live_points -= island_points(island);
            previous     = std::move(current);
            previous_top = top;
        }
    }
    void end_tile()
    {
        if (!tile_open)
            throw std::invalid_argument("Local topology has no open tile.");
        for (const auto& island : previous) live_points -= island_points(island);
        previous.clear();
        if (order == LocalTileOrder::Lexicographic)
            for (unsigned axis = 0; axis < 2; ++axis) {
                for (const auto id : incoming[axis]) release_face(id, axis);
                incoming[axis].clear();
                incoming_first[axis] = 0;
            }
        tile_open = false;
    }
    void add_tile(const LocalSectionTile& tile)
    {
        begin_tile(LocalSectionTile{tile.index, tile.minimum, tile.maximum, {}});
        process_bands(tile, false);
        end_tile();
    }
    template<class Operation> void accept(Operation&& operation)
    {
        if (state != State::Ready) {
            state = State::Failed;
            throw std::invalid_argument("Local topology stream is not accepting input.");
        }
        state = State::Active;
        try {
            operation();
            if (state != State::Active)
                throw std::invalid_argument("Local topology callback reentered its stream.");
            state = State::Ready;
        } catch (...) {
            state = State::Failed;
            throw;
        }
    }
    std::vector<std::size_t> finish()
    {
        if (tile_open)
            throw std::invalid_argument("Local topology cannot finish an incomplete tile.");
        if (order == LocalTileOrder::Lexicographic) {
            for (const auto& entry : pending)
                for (const auto& [id, axis] : entry.second) release_face(id, axis);
            pending.clear();
        }
        for (const auto& [index, ids] : tile_islands)
            for (unsigned axis = 0; axis < 2; ++axis) {
                Index neighbor = index;
                if (neighbor[axis] == std::numeric_limits<int>::max())
                    continue;
                ++neighbor[axis];
                const auto found = tile_islands.find(neighbor);
                if (found == tile_islands.end())
                    continue;
                std::size_t first = 0;
                for (const auto a : ids) {
                    tick();
                    // Record bottoms remain sorted. Extended runs have
                    // nonmonotone tops, so also check every candidate's top.
                    while (first < found->second.size() && records[found->second[first]].top <= records[a].bottom) {
                        tick();
                        ++first;
                    }
                    for (std::size_t j = first; j < found->second.size(); ++j) {
                        const auto b = found->second[j];
                        if (records[b].bottom >= records[a].top)
                            break;
                        tick();
                        if (records[b].top <= records[a].bottom) continue;
                        const auto& p      = records[a].traces[2 * axis + 1];
                        const auto& q      = records[b].traces[2 * axis];
                        const bool contact = contact_traces(p,q,axis,records[a].uncertainty,records[b].uncertainty);
                        if (contact)
                            unite(records[a].label, records[b].label);
                    }
                }
            }
        for (const auto& edge : uncertain_contacts) {
            tick();
            if (root(edge.first) != root(edge.second))
                throw std::invalid_argument("Unresolved component-sensitive adjacent contact; tile=" +
                    std::to_string(edge.tile[0]) + "," + std::to_string(edge.tile[1]) +
                    " z_nm=" + std::to_string(edge.height) + ". First ambiguous observation: " + first_uncertain_contact);
            if (statistics) ++statistics->redundant_contacts;
        }
        if (!emitted_pieces)
            throw std::invalid_argument("Local topology contains no occupied solid.");
        // Preserve the retained path's final per-piece work reservation even though
        // the streaming certificate needs only one mapping entry per label.
        tick(emitted_pieces);
        std::size_t unanchored = 0, largest = labels.size();
        for (std::size_t i = 0; i < labels.size(); ++i) {
            tick();
            if (labels[i].parent != i || labels[i].bed) continue;
            ++unanchored;
            if (largest == labels.size() || labels[i].volume_m3 > labels[largest].volume_m3) largest = i;
        }
        free_components = unanchored;
        if (unanchored && attachment == LocalAttachmentPolicy::RequireBed) {
            const auto& label = labels[largest];
            std::ostringstream message;
            message << "Local topology contains unanchored sampled solid components: count=" << unanchored
                    << " largest_root=" << largest << " minimum_nm=" << label.minimum[0] << ','
                    << label.minimum[1] << ',' << label.minimum[2] << " maximum_nm=" << label.maximum[0] << ','
                    << label.maximum[1] << ',' << label.maximum[2] << " volume_m3=" << label.volume_m3
                    << " pieces=" << label.pieces << " labels=" << label.members
                    << " islands=" << label.geometric_islands << " bed_gap_nm=" << label.minimum[2] - bed;
            throw std::invalid_argument(message.str());
        }
        std::map<std::size_t, std::size_t> component_ids;
        std::vector<std::size_t> mapping;
        mapping.reserve(labels.size());
        for (std::size_t i = 0; i < labels.size(); ++i) {
            if (checkpoint)
                checkpoint();
            const auto id = root(i);
            mapping.push_back(component_ids.emplace(id, component_ids.size()).first->second);
        }
        components = component_ids.size();
        return mapping;
    }
};

LocalTopologyBuilder::LocalTopologyBuilder(double bed,
                                           const SectionTopologyLimits& limits,
                                           std::function<void(const SectionPiece&)> piece,
                                           std::function<void()> checkpoint, std::optional<StreamTopologyLimits> stream_limits,
                                           StreamTopologyStatistics* statistics, LocalTileOrder order, LocalAttachmentPolicy attachment)
    : m_impl(std::make_unique<Impl>(bed, limits, std::move(piece), std::move(checkpoint), stream_limits, statistics, order, attachment)),
      m_source(new LocalTopologySource)
{}
LocalTopologyBuilder::~LocalTopologyBuilder() = default;

void LocalTopologyBuilder::add_tile(const LocalSectionTile& tile)
{
    m_impl->accept([&] { m_impl->add_tile(tile); });
}
void LocalTopologyBuilder::begin_tile(const LocalSectionTile& header)
{
    m_impl->accept([&] { m_impl->begin_tile(header); });
}
void LocalTopologyBuilder::add_band(const LocalSectionTile& packet)
{
    m_impl->accept([&] { m_impl->process_bands(packet, true); });
}
void LocalTopologyBuilder::end_tile()
{
    m_impl->accept([&] { m_impl->end_tile(); });
}

std::vector<std::size_t> LocalTopologyBuilder::provisional_component_mapping()
{
    std::vector<std::size_t> mapping;
    m_impl->accept([&] {
        mapping.reserve(m_impl->labels.size());
        for (std::size_t i = 0; i < m_impl->labels.size(); ++i) {
            auto representative = i;
            m_impl->tick();
            while (m_impl->labels[representative].parent != representative) {
                m_impl->tick();
                representative = m_impl->labels[representative].parent;
            }
            mapping.push_back(representative);
        }
    });
    return mapping;
}

LocalConnectivityCertificate LocalTopologyBuilder::finish()
{
    if (m_impl->state != Impl::State::Ready) {
        m_impl->state = Impl::State::Failed;
        throw std::invalid_argument("Local topology stream cannot be finished.");
    }
    m_impl->state = Impl::State::Active;
    try {
        LocalConnectivityCertificate certificate;
        certificate.m_mapping = m_impl->finish();
        if (m_impl->state != Impl::State::Active)
            throw std::invalid_argument("Local topology callback reentered its stream.");
        certificate.m_components = m_impl->components;
        certificate.m_free_components = m_impl->free_components;
        certificate.m_bed_z_m    = double(m_impl->bed) * 1e-9;
        certificate.m_source     = m_source;
        m_impl->state            = Impl::State::Finished;
        return certificate;
    } catch (...) {
        m_impl->state = Impl::State::Failed;
        throw;
    }
}

SectionTopology qualify_local_section_topology(const std::vector<LocalSectionTile>& tiles,
                                               double bed_z,
                                               const SectionTopologyLimits& limits,
                                               const std::function<void()>& checkpoint)
{
    if (tiles.empty() || tiles.size() > 4096)
        throw std::invalid_argument("Invalid local topology input or limits.");
    SectionTopology result;
    LocalTopologyBuilder builder(bed_z, limits, [&](const SectionPiece& piece) { result.m_pieces.push_back(piece); }, checkpoint);
    for (const auto& tile : tiles)
        builder.add_tile(tile);
    const auto certificate = builder.finish();
    result.m_bed_z_m       = certificate.bed_z_m();
    result.m_components    = certificate.components();
    for (auto& piece : result.m_pieces)
        piece.component = certificate.component_mapping().at(piece.component);
    return result;
}
} // namespace Slic3r::Thermoelastic
