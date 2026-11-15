#include "SectionThermalSkin.hpp"
#include "SectionThermalExposure.hpp"
#include "ThermoelasticExactSections.hpp"
#include "ThermoelasticSectionFootprint.hpp"
#include "ThermoelasticLocalSections.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include "ClipperUtils.hpp"
#include "AABBTreeIndirect.hpp"
#include "BoundingBox.hpp"
#include "ExPolygon.hpp"
#include "Point.hpp"
#include "libslic3r.h"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
namespace {
std::size_t points(const ExPolygons& polygons)
{
    std::size_t n = 0;
    for (const auto& p : polygons) {
        n += p.contour.points.size();
        for (const auto& h : p.holes)
            n += h.points.size();
    }
    return n;
}
coord_t native(double z)
{
    if (!std::isfinite(z) || std::abs(z) > 10. || std::abs(z * 1e9 - std::round(z * 1e9)) > 1e-5)
        throw std::invalid_argument("Invalid native skin plane.");
    return std::llround(z * 1e9);
}
struct Seam
{
    std::size_t material, component;
    unsigned axis;
    coord_t plane, low, high, bottom, top;
    bool positive;
};
} // namespace
struct SectionThermalSkin::Impl
{
    using Index = std::array<int, 2>;
    Horizontal horizontal;
    Vertical vertical;
    ExactHorizontal exact_horizontal;
    ExactVertical exact_vertical;
    struct ExactSeam {
        Seam bounds;
        ExactSectionBoundary span;
        std::array<double,2> normal;
    };
    std::map<Index,std::vector<ExactSeam>> exact_pending;
    std::vector<ExactSeam> exact_incoming, exact_negative;
    std::size_t exact_live_points=0;
    SectionSkinLimits limits;
    std::function<void()> checkpoint;
    ComponentMapping component_mapping;
    std::size_t work = 0, seams = 0;
    enum class WorkStage { Stream, Compaction, Indexing, Subtraction, Horizontal, Vertical };
    WorkStage work_stage = WorkStage::Stream;
    std::array<std::size_t, 6> stage_work{};
    std::map<Index, std::vector<Seam>> pending;
    std::vector<Seam> incoming, negative;
    std::vector<SectionPiece> previous;
    std::vector<BoundingBox> previous_bounds;
    std::optional<Index> last;
    LocalSectionTile tile;
    bool open = false, failed = false, finished = false, busy = false;

    Impl(Horizontal h, Vertical v, SectionSkinLimits l, std::function<void()> c, ComponentMapping mapping,
        ExactHorizontal eh, ExactVertical ev)
        : horizontal(std::move(h)), vertical(std::move(v)), limits(l), checkpoint(std::move(c)),
          component_mapping(std::move(mapping))
    {
        if (!horizontal || !vertical || !limits.max_live_points || limits.max_live_points > 2000000 || !limits.max_seam_patches ||
            limits.max_seam_patches > 250000 || !limits.max_work || limits.max_work > 2000000000)
            throw std::invalid_argument("Invalid streamed skin callbacks or limits.");
        if (bool(eh)!=bool(ev)) throw std::invalid_argument("Incomplete exact skin callbacks.");
        exact_horizontal=std::move(eh);
        exact_vertical=std::move(ev);
    }
    template<class F> void invoke(F&& f)
    {
        if (failed || finished || busy) {
            failed = true;
            throw std::invalid_argument("Inactive or reentered skin stream.");
        }
        busy = true;
        try {
            f();
            if (failed)
                throw std::invalid_argument("Reentered skin stream.");
            busy = false;
        } catch (...) {
            failed = true;
            busy   = false;
            throw;
        }
    }
    std::string work_details() const
    {
        static constexpr std::array<const char*,6> names{
            "stream", "seam compaction", "seam indexing", "seam subtraction", "horizontal face", "vertical face"};
        std::string details=std::string("stage=")+names[std::size_t(work_stage)]+", used="+std::to_string(work)+
            ", limit="+std::to_string(limits.max_work);
        for (std::size_t i=0;i<names.size();++i)
            details+=", "+std::string(names[i])+"="+std::to_string(stage_work[i]);
        return details;
    }
    void charge(std::size_t n)
    {
        if (checkpoint)
            checkpoint();
        if (n > limits.max_work - work) {
            throw ResourceLimit("Streamed skin work budget exceeded: "+work_details()+", requested="+std::to_string(n)+".");
        }
        work += n;
        stage_work[std::size_t(work_stage)] += n;
    }
    std::size_t exact_budget() const
    {
        if (work==limits.max_work) throw ResourceLimit("Exact skin work budget exhausted.");
        return std::min<std::size_t>(4000000,limits.max_work-work);
    }
    void retain(std::size_t n)
    {
        if (n>limits.max_live_points-exact_live_points)
            throw ResourceLimit("Exact skin live geometry exceeded.");
        exact_live_points+=n;
    }
    struct Scratch {
        Impl& owner;
        std::size_t count;
        Scratch(Impl& o,std::size_t n):owner(o),count(n) { owner.retain(n); }
        ~Scratch() { owner.exact_live_points-=count; }
    };
    ExactSection exact_geometry(const SectionFootprint& footprint)
    {
        if (footprint.is_exact()) return footprint.exact();
        Scratch scratch(*this,4096+footprint.vertex_count());
        auto result=make_exact_section(footprint.native(),exact_budget(),[&] { charge(0); });
        charge(result.work);
        return std::move(result.section);
    }
    ExactSection subtract(const ExactSection& a,const ExactSection& b)
    {
        Scratch scratch(*this,4096+a.vertex_count()+b.vertex_count());
        bool callback_failed=false;
        ExactSectionResult result;
        try {
            result=boolean_exact_sections(a,b,ExactSectionOperation::Difference,exact_budget(),[&] {
                try { charge(0); } catch (...) { callback_failed=true; throw; }
            });
        } catch (const ResourceLimit& error) {
            if (callback_failed) throw;
            throw ResourceLimit("Exact skin section subtraction failed: "+work_details()+"; "+error.what());
        }
        charge(result.work);
        return std::move(result.section);
    }
    ExactSectionBoundary subtract(const ExactSectionBoundary& a,const ExactSectionBoundary& b)
    {
        Scratch scratch(*this,8192+2*(a.segment_count()+b.segment_count()));
        bool callback_failed=false;
        ExactSectionBoundaryResult result;
        try {
            result=boolean_exact_section_boundaries(a,b,ExactBoundaryOperation::Difference,
                ExactBoundaryOrientation::Opposed,exact_budget(),[&] {
                    try { charge(0); } catch (...) { callback_failed=true; throw; }
                });
        } catch (const ResourceLimit& error) {
            if (callback_failed) throw;
            throw ResourceLimit("Exact skin boundary subtraction failed: "+work_details()+"; "+error.what());
        }
        charge(result.work);
        return std::move(result.boundary);
    }
    void exact_emit(const ExactSeam& s)
    {
        charge(1);
        exact_vertical(s.bounds.material,s.bounds.component,s.span,s.bounds.bottom,s.bounds.top,s.normal);
    }
    void exact_append(std::vector<ExactSeam>& target,ExactSeam s)
    {
        if (seams==limits.max_seam_patches) {
            exact_compact(exact_incoming);
            exact_compact(exact_negative);
            for (auto& entry:exact_pending) exact_compact(entry.second);
        }
        if (seams==limits.max_seam_patches) throw ResourceLimit("Exact skin seam budget exceeded.");
        retain(2*s.span.segment_count());
        target.push_back(std::move(s));
        ++seams;
    }
    void exact_compact(std::vector<ExactSeam>& patches)
    {
        work_stage=WorkStage::Compaction;
        const auto key=[](const ExactSeam& s) {
            return std::tie(s.bounds.axis,s.bounds.plane,s.bounds.low,s.bounds.high,
                s.bounds.material,s.bounds.component,s.bounds.positive);
        };
        std::size_t levels=1;
        for (auto n=patches.size();n>1;n/=2) ++levels;
        charge(patches.size()*levels);
        std::sort(patches.begin(),patches.end(),[&](const auto& a,const auto& b) {
            return key(a)<key(b) || (key(a)==key(b) && a.bounds.bottom<b.bounds.bottom);
        });
        std::size_t kept=0;
        for (std::size_t i=0;i<patches.size();++i) {
            bool equal=false;
            const auto& patch=patches[i];
            if (kept && key(patches[kept-1])==key(patch) && patches[kept-1].bounds.top==patch.bounds.bottom) {
                Scratch scratch(*this,8192+2*(patch.span.segment_count()+patches[kept-1].span.segment_count()));
                const auto covered=[&](const auto& a,const auto& b) {
                    const auto remainder=boolean_exact_section_boundaries(a,b,ExactBoundaryOperation::Difference,
                        ExactBoundaryOrientation::Any,exact_budget(),[&] { charge(0); });
                    charge(remainder.work);
                    return remainder.boundary.empty();
                };
                equal=covered(patch.span,patches[kept-1].span) && covered(patches[kept-1].span,patch.span);
            }
            if (equal) {
                patches[kept-1].bounds.top=patch.bounds.top;
                exact_live_points-=2*patch.span.segment_count();
                patches[i].span={};
                --seams;
            } else {
                if (kept!=i) patches[kept]=std::move(patches[i]);
                ++kept;
            }
        }
        patches.resize(kept);
    }
    void exact_clear(std::vector<ExactSeam>& patches)
    {
        for (const auto& s:patches) exact_live_points-=2*s.span.segment_count();
        seams-=patches.size();
        patches.clear();
    }
    void exact_face(const SectionPiece& p,const std::vector<SectionPiece>& opposing,bool upward)
    {
        work_stage=WorkStage::Horizontal;
        auto exposed=p.slab.footprint.exact();
        Scratch held(*this,exposed.vertex_count());
        const auto plane=native(upward?p.slab.top_m:p.slab.bottom_m);
        // Repeated bands often have an identical same-material face. Try it
        // before unrelated overlapping bounds, but retain every other material
        // in the second pass: only exact subtraction proves face coverage.
        for (const bool same_material:{true,false}) {
            if (exposed.empty()) break;
            for (const auto& other:opposing) {
                charge(1);
                if (exposed.empty()) break;
                if ((other.slab.material==p.slab.material)!=same_material) continue;
                if (native(upward?other.slab.bottom_m:other.slab.top_m)!=plane) continue;
                const auto a=*exposed.bounds(),b=*other.slab.footprint.bounds();
                // Horizontal faces carry area: an outward-box line/point contact
                // cannot cover positive area. Vertical shared edges deliberately
                // retain the separate closed-box test below.
                if (std::max(a.min.x(),b.min.x())>=std::min(a.max.x(),b.max.x()) ||
                    std::max(a.min.y(),b.min.y())>=std::min(a.max.y(),b.max.y())) continue;
                auto cut=subtract(exposed,other.slab.footprint.exact());
                retain(cut.vertex_count());
                exact_live_points-=held.count;
                held.count=cut.vertex_count();
                exposed=std::move(cut);
            }
        }
        if (!exposed.empty()) exact_horizontal(p.slab.material,p.component,exposed,plane,upward);
    }
    void exact_band(std::vector<SectionPiece> current)
    {
        charge(1);
        if (!open || current.empty()) throw std::invalid_argument("Missing exact skin tile or band.");
        const auto bottom=native(current.front().slab.bottom_m),top=native(current.front().slab.top_m);
        if (top<=bottom || (!previous.empty() && native(previous.front().slab.top_m)>bottom))
            throw std::invalid_argument("Overlapping exact skin bands.");
        std::size_t input_points=0;
        for (const auto& p:current) {
            charge(1);
            if (p.slab.footprint.vertex_count()>limits.max_live_points-input_points)
                throw ResourceLimit("Exact skin input geometry exceeded.");
            input_points+=p.slab.footprint.vertex_count();
        }
        Scratch borrowed(*this,input_points);
        for (auto& p:current) {
            charge(1);
            if (native(p.slab.bottom_m)!=bottom || native(p.slab.top_m)!=top)
                throw std::invalid_argument("Exact skin partition is not one complete band.");
            auto section=exact_geometry(p.slab.footprint);
            if (section.empty()) throw std::invalid_argument("Empty exact skin piece.");
            retain(section.vertex_count());
            p.slab.footprint=SectionFootprint(std::move(section));
        }
        for (const auto& p:previous) exact_face(p,current,true);
        std::vector<ExactSectionBoundary> boundaries;
        std::size_t boundary_points=0;
        for (const auto& p:current) {
            Scratch scratch(*this,8192);
            auto boundary=extract_exact_section_boundary(p.slab.footprint.exact(),exact_budget(),[&] { charge(0); });
            charge(boundary.work);
            const auto count=2*boundary.boundary.segment_count();
            retain(count); boundary_points+=count;
            boundaries.push_back(std::move(boundary.boundary));
        }
        for (std::size_t i=0;i<current.size();++i) {
            const auto& p=current[i];
            exact_face(p,previous,false);
            work_stage=WorkStage::Vertical;
            auto exposed=boundaries[i];
            Scratch held(*this,2*exposed.segment_count());
            for (std::size_t j=0;j<current.size();++j) {
                charge(1);
                if (i==j || !p.slab.footprint.bounds()->overlap(*current[j].slab.footprint.bounds())) continue;
                auto cut=subtract(exposed,boundaries[j]);
                retain(2*cut.segment_count());
                exact_live_points-=held.count; held.count=2*cut.segment_count();
                exposed=std::move(cut);
            }
            Scratch scratch(*this,8192);
            auto split=split_exact_boundary_segments(exposed,exact_budget(),[&] { charge(0); });
            charge(split.work);
            for (auto& patch:split.segments) {
                bool seam=false;
                for (unsigned axis=0;axis<2 && !seam;++axis) for (unsigned side=0;side<2 && !seam;++side) {
                    const auto plane=side?tile.maximum[axis]:tile.minimum[axis];
                    const auto selected=exact_boundary_on_native_plane(patch.span,axis,plane,exact_budget(),[&] { charge(0); });
                    charge(selected.work);
                    if (selected.boundary.empty()) continue;
                    const auto box=*patch.span.bounds();
                    ExactSeam s{{p.slab.material,p.component,axis,plane,box.min[1-axis],box.max[1-axis],bottom,top,bool(side)},
                        patch.span,patch.outward_normal};
                    if (side && tile.index[axis]!=std::numeric_limits<int>::max()) {
                        auto target=tile.index; ++target[axis]; exact_append(exact_pending[target],std::move(s));
                    } else if (side) exact_emit(s);
                    else exact_append(exact_negative,std::move(s));
                    seam=true;
                }
                if (!seam) exact_vertical(p.slab.material,p.component,patch.span,bottom,top,patch.outward_normal);
            }
        }
        exact_live_points-=boundary_points;
        for (const auto& p:previous) exact_live_points-=p.slab.footprint.vertex_count();
        previous=std::move(current);
    }
    void emit(const Seam& s)
    {
        charge(1);
        Point a, b;
        a[s.axis] = b[s.axis] = s.plane;
        a[1 - s.axis]         = s.low;
        b[1 - s.axis]         = s.high;
        // Material-left edge orientation: right normal is outward.
        if (s.positive != (s.axis == 0))
            std::swap(a, b);
        Eigen::Vector2d normal = Eigen::Vector2d::Zero();
        normal[s.axis]         = s.positive ? 1. : -1.;
        vertical(s.material, s.component, {a, b, 0., 1., s.bottom, s.top, normal});
    }
    void append(std::vector<Seam>& target, const Seam& s)
    {
        // Reclaim redundant adjacent-band traces before exhausting the live
        // frontier. No callbacks are emitted and no owner/span is changed.
        if (seams == limits.max_seam_patches) {
            const auto stage = work_stage;
            if (component_mapping) {
                const auto mapping = component_mapping();
                if (mapping.size() > 250000) throw ResourceLimit("Skin component mapping budget exceeded.");
                charge(mapping.size());
                for (std::size_t i = 0; i < mapping.size(); ++i)
                    if (mapping[i] >= mapping.size() || mapping[mapping[i]] != mapping[i])
                        throw std::invalid_argument("Invalid provisional skin component mapping.");
                const auto relabel = [&](std::vector<Seam>& patches) {
                    charge(patches.size());
                    for (auto& patch : patches) {
                        if (patch.component >= mapping.size())
                            throw std::invalid_argument("Missing provisional skin component.");
                        patch.component = mapping[patch.component];
                    }
                };
                relabel(incoming);
                relabel(negative);
                for (auto& entry : pending) relabel(entry.second);
            }
            compact(incoming);
            compact(negative);
            for (auto& entry : pending) compact(entry.second);
            work_stage = stage;
        }
        if (seams == limits.max_seam_patches)
            throw ResourceLimit("Streamed skin seam budget exceeded.");
        target.push_back(s);
        ++seams;
    }
    void compact(std::vector<Seam>& patches)
    {
        work_stage = WorkStage::Compaction;
        // Adjacent quadrature bands often carry the identical seam trace.
        // Merge only identical owners and spans; do not merge across gaps or
        // material/component changes. This bounds later contact subtraction
        // by actual trace changes instead of the sampling plane count.
        std::size_t levels = 1;
        for (auto n = patches.size(); n > 1; n /= 2)
            ++levels;
        charge(patches.size() * levels);
        const auto key = [](const Seam& s) { return std::tie(s.axis, s.plane, s.low, s.high, s.material, s.component, s.positive); };
        std::sort(patches.begin(), patches.end(),
                  [&](const Seam& a, const Seam& b) { return key(a) < key(b) || (key(a) == key(b) && a.bottom < b.bottom); });
        std::size_t kept = 0;
        for (const auto& patch : patches) {
            if (kept && key(patches[kept - 1]) == key(patch) && patches[kept - 1].top == patch.bottom)
                patches[kept - 1].top = patch.top;
            else
                patches[kept++] = patch;
        }
        seams -= patches.size() - kept;
        patches.resize(kept);
        // Collinear subdivisions at one height are equally redundant. Keep
        // both ownership labels and the exact Z interval: edge contact alone
        // must not turn separate mechanical components into one component.
        levels = 1;
        for (auto n = patches.size(); n > 1; n /= 2)
            ++levels;
        charge(patches.size() * levels);
        const auto horizontal_key = [](const Seam& s) {
            return std::tie(s.axis, s.plane, s.bottom, s.top, s.material, s.component, s.positive);
        };
        std::sort(patches.begin(), patches.end(), [&](const Seam& a, const Seam& b) {
            return horizontal_key(a) < horizontal_key(b) ||
                (horizontal_key(a) == horizontal_key(b) && a.low < b.low);
        });
        kept = 0;
        for (const auto& patch : patches) {
            if (kept && horizontal_key(patches[kept - 1]) == horizontal_key(patch) && patches[kept - 1].high == patch.low)
                patches[kept - 1].high = patch.high;
            else
                patches[kept++] = patch;
        }
        seams -= patches.size() - kept;
        patches.resize(kept);
    }
    // All native integer coordinates are exactly representable as doubles.
    // Plane and axis are separate dimensions, excluding unrelated tile faces.
    using SeamIndex = AABBTreeIndirect::Tree<4,double>;
    static SeamIndex::BoundingBox seam_bounds(const Seam& s)
    {
        for (const auto value : {s.low,s.high,s.bottom,s.top,s.plane})
            if (std::abs(double(value)) > 1e10) throw std::invalid_argument("Seam index coordinate outside supported range.");
        return {SeamIndex::VectorType(double(s.low),double(s.bottom),double(s.plane),double(s.axis)),
                SeamIndex::VectorType(double(s.high),double(s.top),double(s.plane),double(s.axis))};
    }
    SeamIndex index_seams(const std::vector<Seam>& patches)
    {
        work_stage = WorkStage::Indexing;
        struct Box {
            std::size_t index;
            SeamIndex::BoundingBox bounds;
            SeamIndex::VectorType center;
            Impl* owner;
            std::size_t idx() const { return index; }
            const SeamIndex::BoundingBox& bbox() const { owner->charge(1); return bounds; }
            const SeamIndex::VectorType& centroid() const { owner->charge(1); return center; }
        };
        std::size_t power = 1;
        while (power < patches.size()) power *= 2;
        if ((2*power-1)*sizeof(SeamIndex::Node)+patches.size()*sizeof(Box) > 64*1024*1024)
            throw ResourceLimit("Streamed seam index memory budget exceeded.");
        std::vector<Box> boxes;
        boxes.reserve(patches.size());
        for (std::size_t i = 0; i < patches.size(); ++i) {
            charge(1);
            const auto bounds = seam_bounds(patches[i]);
            boxes.push_back({i,bounds,bounds.center(),this});
        }
        SeamIndex index;
        index.build(std::move(boxes));
        return index;
    }
    // Rectangle subtraction in native (edge-coordinate,Z), preserving the
    // owner of each exposed remnant. Used in both directions at a tile seam.
    void difference(const Seam& source, const std::vector<Seam>& covers, const SeamIndex& index)
    {
        work_stage = WorkStage::Subtraction;
        if (covers.empty()) {
            emit(source);
            return;
        }
        std::size_t levels = 1;
        for (auto n = covers.size(); n > 1; n /= 2)
            ++levels;
        charge(levels);
        const auto query = seam_bounds(source);
        std::vector<std::size_t> candidates;
        std::vector<std::size_t> pending_nodes{0};
        while (!pending_nodes.empty()) {
            charge(1);
            const auto id = pending_nodes.back();
            pending_nodes.pop_back();
            const auto& node = index.node(id);
            if (!node.is_valid() || !node.bbox.intersects(query)) continue;
            if (node.is_inner()) {
                pending_nodes.push_back(SeamIndex::left_child_idx(id));
                pending_nodes.push_back(SeamIndex::right_child_idx(id));
                continue;
            }
            const auto i = node.idx;
            const auto& cover = covers[i];
            // A fully covered trace needs no scratch fragments. This also
            // permits cancellation at a full live frontier without enlarging
            // its storage budget or dropping a partially exposed surface.
            if (source.axis == cover.axis && source.plane == cover.plane &&
                cover.bottom <= source.bottom && cover.top >= source.top &&
                cover.low <= source.low && cover.high >= source.high)
                return;
            if (source.axis == cover.axis && source.plane == cover.plane && source.bottom < cover.top && source.top > cover.bottom &&
                source.low < cover.high && source.high > cover.low)
                candidates.push_back(i);
        }
        if (candidates.empty()) {
            emit(source);
            return;
        }
        // Retain the original subtraction order and resulting patch order.
        charge(candidates.size() * levels);
        std::sort(candidates.begin(), candidates.end());
        // Depth-first subtraction emits a finished fragment immediately.
        // Scratch is fixed-size and independent of the retained tile frontier;
        // a full frontier must not prevent a single partial trace from being
        // processed. Each fragment resumes AFTER the cover that produced it.
        struct Fragment { Seam seam; std::size_t next; };
        const auto scratch_limit = std::min<std::size_t>(4096,limits.max_seam_patches);
        std::vector<Fragment> scratch;
        scratch.reserve(scratch_limit);
        scratch.push_back({source,0});
        while (!scratch.empty()) {
            const auto fragment = scratch.back();
            scratch.pop_back();
            const auto& s = fragment.seam;
            bool split = false;
            for (std::size_t cursor = fragment.next; cursor < candidates.size(); ++cursor) {
                charge(1);
                const auto& cover = covers[candidates[cursor]];
                const auto lo = std::max(s.low, cover.low), hi = std::min(s.high, cover.high);
                const auto bottom = std::max(s.bottom, cover.bottom), top = std::min(s.top, cover.top);
                if (lo >= hi || bottom >= top) continue;
                const auto add = [&](coord_t left,coord_t right,coord_t low_z,coord_t high_z) {
                    if (left >= right || low_z >= high_z) return;
                    if (scratch.size() == scratch_limit)
                        throw ResourceLimit("Streamed skin subtraction depth exceeded.");
                    auto piece = s;
                    piece.low = left;
                    piece.high = right;
                    piece.bottom = low_z;
                    piece.top = high_z;
                    scratch.push_back({piece,cursor+1});
                };
                // Reverse push order preserves the old fragment traversal.
                add(lo,hi,top,s.top);
                add(lo,hi,s.bottom,bottom);
                add(hi,s.high,s.bottom,s.top);
                add(s.low,lo,s.bottom,s.top);
                split = true;
                break;
            }
            if (!split) emit(s);
        }
    }
    void exact_difference(const ExactSeam& source,const std::vector<ExactSeam>& covers,const SeamIndex& index)
    {
        work_stage=WorkStage::Subtraction;
        if (covers.empty()) { exact_emit(source); return; }
        std::vector<std::size_t> candidates, nodes{0};
        std::vector<coord_t> planes{source.bounds.bottom,source.bounds.top};
        const auto query=seam_bounds(source.bounds);
        while (!nodes.empty()) {
            charge(1);
            const auto id=nodes.back(); nodes.pop_back();
            const auto& node=index.node(id);
            if (!node.is_valid() || !node.bbox.intersects(query)) continue;
            if (node.is_inner()) {
                nodes.push_back(SeamIndex::left_child_idx(id));
                nodes.push_back(SeamIndex::right_child_idx(id));
                continue;
            }
            const auto& cover=covers[node.idx].bounds;
            if (cover.axis!=source.bounds.axis || cover.plane!=source.bounds.plane ||
                cover.bottom>=source.bounds.top || cover.top<=source.bounds.bottom) continue;
            candidates.push_back(node.idx);
            planes.push_back(std::max(cover.bottom,source.bounds.bottom));
            planes.push_back(std::min(cover.top,source.bounds.top));
        }
        std::size_t levels=1;
        for (auto n=planes.size();n>1;n/=2) ++levels;
        charge(planes.size()*levels);
        std::sort(planes.begin(),planes.end());
        planes.erase(std::unique(planes.begin(),planes.end()),planes.end());
        for (std::size_t z=1;z<planes.size();++z) {
            auto exposed=source.span;
            Scratch held(*this,2*exposed.segment_count());
            for (const auto id:candidates) {
                charge(1);
                const auto& cover=covers[id];
                if (cover.bounds.bottom>=planes[z] || cover.bounds.top<=planes[z-1]) continue;
                auto cut=subtract(exposed,cover.span);
                retain(2*cut.segment_count());
                exact_live_points-=held.count; held.count=2*cut.segment_count();
                exposed=std::move(cut);
                if (exposed.empty()) break;
            }
            if (!exposed.empty()) {
                auto remnant=source;
                remnant.span=std::move(exposed);
                remnant.bounds.bottom=planes[z-1]; remnant.bounds.top=planes[z];
                exact_emit(remnant);
            }
        }
    }
    SeamIndex exact_index(const std::vector<ExactSeam>& patches)
    {
        std::vector<Seam> boxes;
        boxes.reserve(patches.size());
        for (const auto& patch:patches) { charge(1); boxes.push_back(patch.bounds); }
        return index_seams(boxes);
    }
    void exact_end()
    {
        for (const auto& p:previous) exact_face(p,{},true);
        for (const auto& p:previous) exact_live_points-=p.slab.footprint.vertex_count();
        previous.clear();
        exact_compact(exact_incoming);
        exact_compact(exact_negative);
        // Other queued neighbors are unchanged by this tile. Re-sorting the
        // entire next-row frontier on every tile adds no compaction benefit.
        // Full-frontier reclamation remains in exact_append at the seam cap.
        for (unsigned axis=0;axis<2;++axis) {
            if (tile.index[axis]==std::numeric_limits<int>::max()) continue;
            auto target=tile.index; ++target[axis];
            const auto found=exact_pending.find(target);
            if (found!=exact_pending.end()) exact_compact(found->second);
        }
        const auto negative_index=exact_index(exact_negative), incoming_index=exact_index(exact_incoming);
        for (const auto& s:exact_incoming) exact_difference(s,exact_negative,negative_index);
        for (const auto& s:exact_negative) exact_difference(s,exact_incoming,incoming_index);
        exact_clear(exact_incoming); exact_clear(exact_negative);
        open=false;
    }
    void flush_missing(const Index& index)
    {
        while (!pending.empty() && pending.begin()->first < index) {
            for (const auto& s : pending.begin()->second)
                emit(s);
            seams -= pending.begin()->second.size();
            pending.erase(pending.begin());
        }
    }
    void begin(const LocalSectionTile& header)
    {
        charge(1);
        if (open || (last && header.index <= *last) || !header.slabs.empty())
            throw std::invalid_argument("Invalid skin tile order or header.");
        for (unsigned a = 0; a < 2; ++a)
            if (header.minimum[a] >= header.maximum[a])
                throw std::invalid_argument("Invalid skin tile bounds.");
        if (exact_horizontal) {
            while (!exact_pending.empty() && exact_pending.begin()->first<header.index) {
                for (const auto& s:exact_pending.begin()->second) exact_emit(s);
                exact_clear(exact_pending.begin()->second);
                exact_pending.erase(exact_pending.begin());
            }
        } else flush_missing(header.index);
        tile       = header;
        last       = header.index;
        open       = true;
        if (exact_horizontal) {
            const auto found=exact_pending.find(header.index);
            if (found!=exact_pending.end()) {
                exact_incoming=std::move(found->second);
                exact_pending.erase(found);
            }
            return;
        }
        auto found = pending.find(header.index);
        if (found != pending.end()) {
            incoming = std::move(found->second);
            pending.erase(found);
        }
    }
    void face(const SectionPiece& p, const BoundingBox& bounds, const std::vector<SectionPiece>& opposing,
              const std::vector<BoundingBox>& opposing_bounds, bool upward)
    {
        work_stage = WorkStage::Horizontal;
        ExPolygons exposed  = p.slab.footprint.native();
        const coord_t plane = native(upward ? p.slab.top_m : p.slab.bottom_m);
        for (std::size_t i = 0; i < opposing.size(); ++i) {
            const auto& other = opposing[i];
            charge(1);
            if (native(upward ? other.slab.bottom_m : other.slab.top_m) != plane || !bounds.overlap(opposing_bounds[i]))
                continue;
            charge(points(exposed) + points(other.slab.footprint.native()));
            // Identical adjacent footprints cover the whole face, irrespective
            // of material. Preserve exact vertices; no tolerance hides gaps.
            if (exposed == other.slab.footprint.native())
                return;
            exposed = diff_ex(exposed, other.slab.footprint.native());
            if (points(exposed) > limits.max_live_points)
                throw ResourceLimit("Streamed skin face budget exceeded.");
            if (exposed.empty())
                break;
        }
        if (!exposed.empty())
            horizontal(p.slab.material, p.component, exposed, plane, upward);
    }
    void band(std::vector<SectionPiece> current)
    {
        if (exact_horizontal) { exact_band(std::move(current)); return; }
        work_stage = WorkStage::Stream;
        charge(1);
        if (!open || current.empty())
            throw std::invalid_argument("Missing skin tile or band.");
        const coord_t bottom = native(current.front().slab.bottom_m), top = native(current.front().slab.top_m);
        if (top <= bottom || (!previous.empty() && native(previous.front().slab.top_m) > bottom))
            throw std::invalid_argument("Overlapping skin bands.");
        std::size_t n = 0;
        for (const auto& p : previous)
            n += points(p.slab.footprint.native());
        for (const auto& p : current) {
            n += points(p.slab.footprint.native());
            if (n > limits.max_live_points)
                throw ResourceLimit("Streamed skin live geometry exceeded.");
            if (native(p.slab.bottom_m) != bottom || native(p.slab.top_m) != top)
                throw std::invalid_argument("Skin partition is not one complete band.");
        }
        std::vector<BoundingBox> current_bounds;
        current_bounds.reserve(current.size());
        for (const auto& p : current) {
            charge(points(p.slab.footprint.native())+1);
            current_bounds.push_back(get_extents(p.slab.footprint.native()));
        }
        for (std::size_t i = 0; i < previous.size(); ++i)
            face(previous[i], previous_bounds[i], current, current_bounds, true);
        work_stage = WorkStage::Vertical;
        visit_vertical_section_patches(current, limits.max_seam_patches,
            [&](const SectionPiece& p, const std::vector<VerticalSectionPatch>& patches) {
            face(p, current_bounds[&p-current.data()], previous, previous_bounds, false);
            work_stage = WorkStage::Vertical;
            for (const auto& patch : patches) {
                bool seam = false;
                for (unsigned axis = 0; axis < 2; ++axis) {
                    if (patch.a[axis] != patch.b[axis])
                        continue;
                    const auto plane = patch.a[axis];
                    if (plane != tile.minimum[axis] && plane != tile.maximum[axis])
                        continue;
                    const auto coordinate = [&](double t) {
                        return std::llround(double(patch.a[1 - axis]) + t * double(patch.b[1 - axis] - patch.a[1 - axis]));
                    };
                    const auto a = coordinate(patch.begin), b = coordinate(patch.end);
                    Seam s{p.slab.material, p.component,    axis,
                           plane,           std::min(a, b), std::max(a, b),
                           patch.bottom,    patch.top,      plane == tile.maximum[axis]};
                    if (s.positive && tile.index[axis] != std::numeric_limits<int>::max()) {
                        auto target = tile.index;
                        ++target[axis];
                        append(pending[target], s);
                    } else if (s.positive)
                        emit(s);
                    else
                        append(negative, s);
                    seam = true;
                    break;
                }
                if (!seam)
                    vertical(p.slab.material,p.component,patch);
            }
        }, [&](std::size_t count) { charge(count); });
        previous = std::move(current);
        previous_bounds = std::move(current_bounds);
    }
    void end()
    {
        charge(1);
        if (!open)
            throw std::invalid_argument("No active skin tile.");
        if (exact_horizontal) { exact_end(); return; }
        for (std::size_t i = 0; i < previous.size(); ++i)
            face(previous[i], previous_bounds[i], {}, {}, true);
        previous.clear();
        previous_bounds.clear();
        compact(incoming);
        compact(negative);
        // Pending neighboring tile traces share the same live budget as
        // subtraction scratch. Reclaim their adjacent-band duplicates before
        // allocating fragments, not only when appending the next trace.
        for (unsigned axis = 0; axis < 2; ++axis) {
            if (tile.index[axis] == std::numeric_limits<int>::max())
                continue;
            auto target = tile.index;
            ++target[axis];
            const auto found = pending.find(target);
            if (found != pending.end())
                compact(found->second);
        }
        const auto negative_index = index_seams(negative);
        const auto incoming_index = index_seams(incoming);
        for (const auto& s : incoming)
            difference(s, negative, negative_index);
        for (const auto& s : negative)
            difference(s, incoming, incoming_index);
        seams -= incoming.size() + negative.size();
        incoming.clear();
        negative.clear();
        open = false;
    }
    void finish()
    {
        charge(1);
        if (open || !last)
            throw std::invalid_argument("Incomplete skin stream.");
        if (exact_horizontal) {
            for (auto& entry:exact_pending) {
                for (const auto& s:entry.second) exact_emit(s);
                exact_clear(entry.second);
            }
            exact_pending.clear();
            finished=true;
            return;
        }
        for (const auto& [index, patches] : pending)
            for (const auto& s : patches)
                emit(s);
        pending.clear();
        seams    = 0;
        finished = true;
    }
};
SectionThermalSkin::SectionThermalSkin(Horizontal h, Vertical v, const SectionSkinLimits& l, std::function<void()> c,
    ComponentMapping mapping, ExactHorizontal eh, ExactVertical ev)
    : m_impl(std::make_unique<Impl>(std::move(h), std::move(v), l, std::move(c), std::move(mapping),std::move(eh),std::move(ev)))
{}
SectionThermalSkin::~SectionThermalSkin() = default;
void SectionThermalSkin::begin_tile(const LocalSectionTile& h)
{
    m_impl->invoke([&] { m_impl->begin(h); });
}
void SectionThermalSkin::add_band(std::vector<SectionPiece> p)
{
    m_impl->invoke([&] { m_impl->band(std::move(p)); });
}
void SectionThermalSkin::end_tile()
{
    m_impl->invoke([&] { m_impl->end(); });
}
void SectionThermalSkin::finish()
{
    m_impl->invoke([&] { m_impl->finish(); });
}
} // namespace Slic3r::Thermoelastic
