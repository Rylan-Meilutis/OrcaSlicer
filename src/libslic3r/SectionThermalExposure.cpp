#include "SectionThermalExposure.hpp"
#include "SectionIntegration.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include "ThermoelasticMesh.hpp"
#include "ExPolygon.hpp"
#include "BoundingBox.hpp"
#include "Polygon.hpp"
#include "Point.hpp"
#include "libslic3r.h"
#include <Eigen/Core>
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
namespace {
using Exact = boost::multiprecision::int128_t;
std::size_t point_count(const ExPolygons& polygons)
{
    std::size_t count = 0;
    for (const auto& p : polygons) {
        count += p.contour.points.size();
        for (const auto& h : p.holes) count += h.points.size();
    }
    return count;
}
using LineKey = std::tuple<coord_t, coord_t, Exact>;
LineKey line_key(const Point& a, const Point& b)
{
    coord_t dx = b.x()-a.x(), dy = b.y()-a.y();
    const auto divisor = std::gcd(dx,dy);
    if (!divisor) throw std::invalid_argument("Degenerate section boundary edge.");
    dx /= divisor;
    dy /= divisor;
    if (dx < 0 || (dx == 0 && dy < 0)) { dx = -dx; dy = -dy; }
    return {dx,dy,Exact(dx)*a.y()-Exact(dy)*a.x()};
}
template<class F> void edges(const ExPolygons& polygons,const F& visit)
{
    for (const auto& polygon : polygons) {
        const auto ring = [&](const Polygon& boundary,bool hole) {
            const bool forward = (boundary.area() > 0.) != hole;
            for (std::size_t i = 0; i < boundary.points.size(); ++i) {
                const auto& a = boundary.points[i];
                const auto& b = boundary.points[(i+1)%boundary.points.size()];
                if (forward) visit(a,b); else visit(b,a);
            }
        };
        ring(polygon.contour,false);
        for (const auto& hole : polygon.holes) ring(hole,true);
    }
}
}

std::vector<VerticalSectionPatch> exposed_vertical_section_patches(const SectionTopology& topology,
    const SectionPiece& piece,std::size_t max_patches,const std::function<void(std::size_t)>& charge)
{
    return partition_vertical_section_patches(topology.pieces(),piece,max_patches,charge);
}

namespace {
struct LineContact {
    std::size_t owner;
    coord_t low,high,bottom,top;
};
using ContactIndex = std::map<LineKey,std::vector<LineContact>>;
ContactIndex partition_contacts(const std::vector<SectionPiece>& partition,
    const std::function<void(std::size_t)>& charge)
{
    std::size_t count = 0;
    for (const auto& piece : partition) {
        const auto n = point_count(piece.slab.footprint.native());
        charge(n+1);
        if (n > 250000-count) throw ResourceLimit("Section band contact index budget exceeded.");
        count += n;
    }
    if (partition.size() < 2) return {};
    // Small bands commonly contain separate islands. Their edges cannot
    // contact, so avoid building an exact line index for those islands.
    // Bound this pairwise broad phase; larger partitions use the index directly.
    std::vector<bool> candidate(partition.size(),true);
    if (partition.size() <= 16) {
        std::fill(candidate.begin(),candidate.end(),false);
        std::vector<BoundingBox> bounds;
        for (const auto& piece : partition) bounds.push_back(get_extents(piece.slab.footprint.native()));
        for (std::size_t i = 0; i < partition.size(); ++i)
            for (std::size_t j = i+1; j < partition.size(); ++j) {
                charge(1);
                if (bounds[i].overlap(bounds[j])) candidate[i] = candidate[j] = true;
            }
    }
    std::size_t lookup_work = 1;
    for (auto n = count; n > 1; n /= 2) ++lookup_work;
    ContactIndex contacts;
    for (std::size_t i = 0; i < partition.size(); ++i) {
        if (!candidate[i]) continue;
        const auto& slab = partition[i].slab;
        const coord_t bottom = std::llround(slab.bottom_m*1e9), top = std::llround(slab.top_m*1e9);
        edges(slab.footprint.native(),[&](const Point& a,const Point& b) {
            charge(lookup_work);
            const unsigned axis = std::abs(b.x()-a.x()) >= std::abs(b.y()-a.y()) ? 0 : 1;
            contacts[line_key(a,b)].push_back({i,std::min(a[axis],b[axis]),std::max(a[axis],b[axis]),bottom,top});
        });
    }
    return contacts;
}

std::vector<VerticalSectionPatch> indexed_vertical_patches(const std::vector<SectionPiece>& partition,
    std::size_t index, const ContactIndex& contacts,
    std::size_t max_patches,const std::function<void(std::size_t)>& charge)
{
    const auto& piece = partition[index];
    const coord_t bottom = std::llround(piece.slab.bottom_m*1e9), top = std::llround(piece.slab.top_m*1e9);
    struct Cover { coord_t low,high,bottom,top; };
    // The whole band's exact supporting-line index is built once, rather than
    // rescanning neighboring polygons for every piece. No snapping hides gaps.
    std::size_t lookup_work = 1;
    for (auto n = contacts.size(); n > 1; n /= 2) ++lookup_work;
    std::size_t contact_count = 0;
    std::vector<VerticalSectionPatch> result;
    edges(piece.slab.footprint.native(),[&](const Point& a,const Point& b) {
        charge(1);
        const coord_t dx = b.x()-a.x(), dy = b.y()-a.y();
        if (!dx && !dy) throw std::invalid_argument("Degenerate section boundary edge.");
        const unsigned axis = std::abs(dx) >= std::abs(dy) ? 0 : 1;
        const coord_t low = std::min(a[axis],b[axis]), high = std::max(a[axis],b[axis]);
        std::vector<Cover> covers;
        std::vector<coord_t> planes{bottom,top};
        if (!contacts.empty()) {
            charge(lookup_work);
            const auto found = contacts.find(line_key(a,b));
            if (found != contacts.end()) for (const auto& contact : found->second) {
                charge(1);
                if (contact.owner == index) continue;
                const coord_t lo = std::max(low,contact.low), hi = std::min(high,contact.high);
                const coord_t lo_z = std::max(bottom,contact.bottom), hi_z = std::min(top,contact.top);
                if (hi <= lo || hi_z <= lo_z) continue;
                if (contact_count++ == max_patches) throw ResourceLimit("Section surface contact budget exceeded.");
                covers.push_back({lo,hi,lo_z,hi_z});
                planes.push_back(lo_z); planes.push_back(hi_z);
            }
        }
        std::sort(planes.begin(),planes.end());
        planes.erase(std::unique(planes.begin(),planes.end()),planes.end());
        const double length = std::hypot(double(dx),double(dy));
        const Eigen::Vector2d normal(double(dy)/length,-double(dx)/length);
        for (std::size_t z = 1; z < planes.size(); ++z) {
            std::vector<std::pair<coord_t,coord_t>> occupied;
            for (const auto& cover : covers) {
                charge(1);
                if (cover.bottom <= planes[z-1] && cover.top >= planes[z]) occupied.emplace_back(cover.low,cover.high);
            }
            std::sort(occupied.begin(),occupied.end());
            const auto emit = [&](coord_t lo,coord_t hi) {
                if (hi <= lo) return;
                charge(1);
                if (result.size() >= max_patches) throw ResourceLimit("Section vertical patch budget exceeded.");
                const double t0 = double(lo-a[axis])/double(b[axis]-a[axis]);
                const double t1 = double(hi-a[axis])/double(b[axis]-a[axis]);
                result.push_back({a,b,std::min(t0,t1),std::max(t0,t1),planes[z-1],planes[z],normal});
            };
            coord_t cursor = low;
            for (const auto& interval : occupied) {
                emit(cursor,interval.first);
                cursor = std::max(cursor,interval.second);
            }
            emit(cursor,high);
        }
    });
    return result;
}
}

std::vector<VerticalSectionPatch> partition_vertical_section_patches(const std::vector<SectionPiece>& partition,
    const SectionPiece& piece,std::size_t max_patches,const std::function<void(std::size_t)>& charge)
{
    if (!charge || !max_patches) throw std::invalid_argument("Missing vertical surface budget.");
    for (std::size_t i = 0; i < partition.size(); ++i) {
        charge(1);
        if (&partition[i] == &piece)
            return indexed_vertical_patches(partition,i,partition_contacts(partition,charge),max_patches,charge);
    }
    throw std::invalid_argument("Section surface piece does not belong to its topology.");
}

void visit_vertical_section_patches(const std::vector<SectionPiece>& partition, std::size_t max_patches,
    const std::function<void(const SectionPiece&, const std::vector<VerticalSectionPatch>&)>& visitor,
    const std::function<void(std::size_t)>& charge)
{
    if (!charge || !max_patches || !visitor) throw std::invalid_argument("Missing vertical surface visitor or budget.");
    const auto contacts = partition_contacts(partition,charge);
    for (std::size_t i = 0; i < partition.size(); ++i)
        visitor(partition[i],indexed_vertical_patches(partition,i,contacts,max_patches,charge));
}

std::optional<SectionIntegration::Moments> vertical_patch_cell_moments(const VerticalSectionPatch& patch,
    const std::array<coord_t,3>& origin,const std::array<coord_t,3>& size)
{
    const coord_t bottom = std::max(patch.bottom,origin[2]), top = std::min(patch.top,origin[2]+size[2]);
    if (top <= bottom) return {};
    const auto xy = vertical_patch_cell_xy_moments(patch,origin,size);
    if (!xy) return {};
    auto moments = SectionIntegration::extrude_line_moments(*xy,double(size[2])*1e-9,
        double(bottom-origin[2])/double(size[2]),double(top-origin[2])/double(size[2]));
    if (!std::isfinite(moments[0]) || moments[0] <= 0.) throw std::invalid_argument("Invalid section vertical face area.");
    return moments;
}

std::optional<SectionIntegration::LineMoments> vertical_patch_cell_xy_moments(const VerticalSectionPatch& patch,
    const std::array<coord_t,3>& origin,const std::array<coord_t,3>& size)
{
    const std::array<coord_t,2> delta{patch.b.x()-patch.a.x(),patch.b.y()-patch.a.y()};
    const std::array<coord_t,2> inward{-delta[1],delta[0]};
    double lo = patch.begin, hi = patch.end;
    for (unsigned axis = 0; axis < 2; ++axis) {
        const coord_t relative = patch.a[axis]-origin[axis];
        if (!delta[axis]) {
            if (relative < 0 || relative > size[axis] || (relative == 0 && inward[axis] < 0) ||
                (relative == size[axis] && inward[axis] > 0)) return {};
        } else {
            const double first = -double(relative)/double(delta[axis]);
            const double last = double(size[axis]-relative)/double(delta[axis]);
            lo = std::max(lo,std::min(first,last));
            hi = std::min(hi,std::max(first,last));
        }
    }
    if (hi <= lo) return {};
    const auto endpoint = [&](double t) {
        Eigen::Vector2d normalized;
        for (unsigned axis = 0; axis < 2; ++axis) {
            const double relative = double(patch.a[axis]-origin[axis])+t*double(delta[axis]);
            // Only floating arithmetic residue, at most 0.001 native nm.
            // A larger mismatch is not repaired by moving the boundary.
            if (relative < -.001 || relative > double(size[axis])+.001)
                throw std::invalid_argument("Unresolved parametric section face clipping.");
            normalized[axis] = std::clamp(relative,0.,double(size[axis]))/double(size[axis]);
        }
        return normalized;
    };
    auto moments = SectionIntegration::vertical_segment_xy_moments(endpoint(lo),endpoint(hi),
        {double(size[0])*1e-9,double(size[1])*1e-9});
    if (!std::isfinite(moments[0]) || moments[0] <= 0.) throw std::invalid_argument("Invalid section vertical face length.");
    return moments;
}
} // namespace Slic3r::Thermoelastic
