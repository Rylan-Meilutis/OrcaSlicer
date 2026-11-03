#include "ThermoelasticLocalSections.hpp"
#include "ExPolygon.hpp"
#include "Point.hpp"
#include "Polygon.hpp"
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticMesh.hpp"
#include "libslic3r.h"
#include <Eigen/Core>
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
namespace {
coord_t native(double x)
{
    if (!std::isfinite(x) || std::abs(x) > 10.)
        throw std::invalid_argument("Local deposition coordinates exceed the supported range.");
    return std::llround(x * 1e9);
}
coord_t floor_div(coord_t a, coord_t b)
{
    return a / b - (a % b < 0 ? 1 : 0);
}
void charge(std::size_t& used, std::size_t amount, std::size_t maximum, const char* message)
{
    if (amount > maximum - used)
        throw ResourceLimit(std::string(message) + " Used=" + std::to_string(used) +
            ", requested=" + std::to_string(amount) + ", limit=" + std::to_string(maximum) + ".");
    used += amount;
}
std::size_t points(const ExPolygons& polygons)
{
    std::size_t count = 0;
    for (const auto& polygon : polygons) {
        count += polygon.contour.points.size();
        for (const auto& hole : polygon.holes)
            count += hole.points.size();
    }
    return count;
}
bool separated_from_segment(const Point& a, const Point& b,
    const std::array<coord_t,2>& minimum, const std::array<coord_t,2>& maximum)
{
    // The box axes were tested by the candidate AABB enumeration. The only
    // remaining separating axis is perpendicular to the segment. Strict signs
    // retain tangent contact. Coordinates are bounded to tens of metres, but
    // the cross product can exceed int64; use the portable exact integer type.
    using Integer = boost::multiprecision::int128_t;
    const Integer dx = Integer(b.x())-a.x(), dy = Integer(b.y())-a.y();
    bool positive = true, negative = true;
    for (const coord_t x : {minimum[0],maximum[0]})
        for (const coord_t y : {minimum[1],maximum[1]}) {
            const Integer cross = dx*(Integer(y)-a.y())-dy*(Integer(x)-a.x());
            positive = positive && cross>0;
            negative = negative && cross<0;
        }
    return positive || negative;
}
} // namespace

struct LocalSectionCursor::Impl {
    const std::vector<DepositedSegment>& segments;
    double bed_z;
    LocalSectionLimits limits;
    std::function<std::optional<SectionOwnership>(double)> ownership;
    std::function<void()> checkpoint;
    unsigned interval_subdivisions;
    LocalSectionStatistics* statistics;
    std::optional<SectionPartitionEnvelopes> owned_partitions;
    const std::vector<DepositedSegment>* partition_envelopes;
    std::shared_ptr<LocalSectionSamplingSchedule> sampling_schedule;
    std::array<coord_t,2> o, s;
    coord_t bed, height, minimum_interval;
    using Index = std::array<int,2>;
    std::map<Index,std::vector<std::size_t>> buckets;
    decltype(buckets)::const_iterator cursor;
    std::size_t references = 0, index_work = 0;
    LocalSectionMetrics result;
    bool failed = false, active = false;
    struct TileState {
        LocalSectionTile tile;
        ExPolygon clip;
        std::vector<DepositedSegment> selected;
        std::vector<std::pair<coord_t,coord_t>> intervals;
        std::unique_ptr<DepositionSectionCursor> sweep;
        std::size_t reserved_ownership = 0, previous_count = 0;
        std::size_t sampled_plane = 0;
        bool sampling = false;
        coord_t previous_top = 0;
    };
    std::unique_ptr<TileState> current;

    Impl(const std::vector<DepositedSegment>& source,
    const Eigen::Vector2d& origin, const Eigen::Vector2d& size,
    double bed_plane, double maximum_height, const LocalSectionLimits& input_limits,
    std::function<std::optional<SectionOwnership>(double)> owner,
    std::function<void()> check, unsigned subdivisions, LocalSectionStatistics* stats,
    const std::vector<DepositedSegment>* partitions)
    : segments(source),bed_z(bed_plane),limits(input_limits),ownership(std::move(owner)),checkpoint(std::move(check)),
      interval_subdivisions(subdivisions),statistics(stats),partition_envelopes(partitions)
{
    if (statistics) *statistics = {};
    if (!interval_subdivisions || interval_subdivisions > 8 ||
        (interval_subdivisions & (interval_subdivisions - 1)) ||
        segments.empty() || segments.size() > 2000000 ||
        (partition_envelopes && partition_envelopes->size() != segments.size()) || !origin.allFinite() || !size.allFinite() ||
        (size.array() <= 0.).any() || !limits.max_tiles || limits.max_tiles > 16384 ||
        !limits.max_tile_references || limits.max_tile_references > 8000000 ||
        !limits.max_sections || limits.max_sections > 200000 || !limits.max_segment_visits ||
        limits.max_segment_visits > 80000000 || !limits.max_boundary_points ||
        limits.max_boundary_points > 128000000u ||
        !limits.max_index_work || limits.max_index_work > 500000000u ||
        !limits.max_ownership_work || limits.max_ownership_work > 1000000000)
        throw std::invalid_argument("Invalid local deposition section inputs or limits.");
    o = {native(origin.x()), native(origin.y())};
    s = {native(size.x()), native(size.y())};
    bed = native(bed_z); height = native(maximum_height);
    minimum_interval = limits.geometry.half_native_section_planes ? 1 : 2;
    if (s[0] < 2 || s[1] < 2 || height < minimum_interval)
        throw std::invalid_argument("Local deposition sampling must exceed native coordinate resolution.");
    if (limits.partition_increment_m != 0.) {
        if (partition_envelopes)
            throw std::invalid_argument("Automatic and explicit section partition envelopes cannot be combined.");
        owned_partitions=make_section_partition_envelopes(segments,bed_z,limits.partition_increment_m,checkpoint);
        partition_envelopes=&owned_partitions->segments;
    }
    result.uses_partition_envelopes=partition_envelopes != nullptr;
    bool exact = !ownership && !partition_envelopes;
    for (std::size_t i = 0; i < segments.size(); ++i) {
        if (checkpoint) checkpoint();
        const auto& bead = segments[i];
        const coord_t width = native(bead.width_m), h0 = native(bead.height_m);
        const coord_t h1 = native(bead.end_height_m.value_or(bead.height_m));
        if (!bead.start_m.allFinite() || !bead.end_m.allFinite() || width <= 0 || h0 < 0 || h1 < 0 ||
            std::max(h0, h1) < minimum_interval || (bead.end_m.head<2>() - bead.start_m.head<2>()).squaredNorm() == 0. ||
            (bead.section == BeadSection::RoundedRectangle && width < std::max(h0, h1)))
            throw std::invalid_argument("Invalid local deposition bead dimensions.");
        const coord_t z0 = native(bead.start_m.z()), z1 = native(bead.end_m.z());
        if (std::min(z0 - h0, z1 - h1) < bed)
            throw std::invalid_argument("Local deposition bead intersects the bed.");
        if (partition_envelopes) {
            const auto& partition = partition_envelopes->at(i);
            const auto p0 = native(partition.start_m.z()), p1 = native(partition.end_m.z());
            const auto ph0 = native(partition.height_m), ph1 = native(partition.end_height_m.value_or(partition.height_m));
            if (partition.start_m.head<2>() != bead.start_m.head<2>() || partition.end_m.head<2>() != bead.end_m.head<2>() ||
                partition.width_m != bead.width_m || partition.material != bead.material || partition.section != bead.section ||
                p0 < z0 || p1 < z1 || ph0 < 0 || ph1 < 0 || p0 - ph0 > z0 - h0 || p1 - ph1 > z1 - h1 ||
                p0 - ph0 < bed || p1 - ph1 < bed ||
                (z0 - h0 > bed && p0 - ph0 == bed) || (z1 - h1 > bed && p1 - ph1 == bed))
                throw std::invalid_argument("Invalid or newly bed-anchored section partition envelope.");
        }
        exact = exact && z0 == z1 && h0 == h1 && bead.section == BeadSection::RectangularEnvelope;
        const Point start(native(bead.start_m.x()),native(bead.start_m.y()));
        const Point end(native(bead.end_m.x()),native(bead.end_m.y()));
        const coord_t radius = (width + 1) / 2;
        const coord_t filter_radius = std::max(radius,coord_t(std::ceil(double(float(.5*double(width))))))+2;
        Index lower, upper;
        for (unsigned axis = 0; axis < 2; ++axis) {
            const coord_t a = native(bead.start_m[axis]), b = native(bead.end_m[axis]);
            // Round outwards, retaining tangent buckets. Empty clipped sections
            // are discarded later; an AABB cannot omit a rounded endpoint.
            const coord_t lo = floor_div(std::min(a, b) - radius - o[axis], s[axis]);
            const coord_t hi = floor_div(std::max(a, b) + radius - o[axis], s[axis]);
            if (lo < std::numeric_limits<int>::min() || hi > std::numeric_limits<int>::max())
                throw ResourceLimit("Local deposition tile index exceeds its supported range.");
            lower[axis] = int(lo);
            upper[axis] = int(hi);
        }
        const std::size_t nx = std::size_t(std::int64_t(upper[0]) - lower[0] + 1);
        const std::size_t ny = std::size_t(std::int64_t(upper[1]) - lower[1] + 1);
        if (nx > limits.max_index_work / ny)
            throw ResourceLimit("Local deposition tile-index work budget exceeded.");
        charge(index_work, nx * ny, limits.max_index_work, "Local deposition tile-index work budget exceeded.");
        for (std::int64_t x = lower[0]; x <= upper[0]; ++x)
            for (std::int64_t y = lower[1]; y <= upper[1]; ++y) {
                if (checkpoint) checkpoint();
                charge(index_work,4,limits.max_index_work,"Local deposition segment-box work budget exceeded.");
                std::array<coord_t,2> minimum, maximum;
                const std::array<coord_t,2> position{{coord_t(x),coord_t(y)}};
                for (unsigned axis=0;axis<2;++axis) {
                    const coord_t edge=o[axis]+position[axis]*s[axis];
                    // A width-radius square contains every circular/rounded
                    // section offset, including its float radius conversion.
                    // Two extra native units cover endpoint and section-coordinate
                    // rounding; this only rejects
                    // provably irrelevant buckets, never clips bead geometry.
                    minimum[axis]=edge-filter_radius;
                    maximum[axis]=edge+s[axis]+filter_radius;
                }
                if (separated_from_segment(start,end,minimum,maximum)) continue;
                const Index index{int(x), int(y)};
                auto found = buckets.find(index);
                if (found == buckets.end()) {
                    if (buckets.size() == limits.max_tiles)
                        throw ResourceLimit("Local deposition tile budget exceeded.");
                    found = buckets.emplace(index, std::vector<std::size_t>{}).first;
                }
                charge(references, 1, limits.max_tile_references, "Local deposition reference budget exceeded.");
                found->second.push_back(i);
                if (statistics) statistics->indexed_tiles = buckets.size();
            }
    }
    result.piecewise_constant = exact;
    cursor = buckets.begin();
}
    void prepare(const LocalSectionPacketVisitor& visitor,bool schedule_only=false) {
        const auto* packets = &visitor;
        const auto& [index, ids] = *cursor;
        if (checkpoint) checkpoint();
        current = std::make_unique<TileState>();
        current->previous_top = bed;
        auto& tile = current->tile;
        tile.index = index;
        for (unsigned axis = 0; axis < 2; ++axis) {
            tile.minimum[axis] = o[axis] + coord_t(index[axis]) * s[axis];
            tile.maximum[axis] = tile.minimum[axis] + s[axis];
        }
        current->clip  = ExPolygon(Polygon{Point(tile.minimum[0], tile.minimum[1]), Point(tile.maximum[0], tile.minimum[1]),
                                           Point(tile.maximum[0], tile.maximum[1]), Point(tile.minimum[0], tile.maximum[1])});
        auto& selected = current->selected;
        selected.reserve(ids.size());
        std::set<coord_t> boundaries;
        for (const auto id : ids) {
            if (checkpoint)
                checkpoint();
            selected.push_back(segments[id]);
            const auto& bead = partition_envelopes ? partition_envelopes->at(id) : segments[id];
            const coord_t a = native(bead.start_m.z()), b = native(bead.end_m.z());
            for (const coord_t z : {a, b, a - native(bead.height_m), b - native(bead.end_height_m.value_or(bead.height_m))})
                boundaries.insert(z);
            if (boundaries.size() > limits.max_sections - result.sections + 1)
                throw ResourceLimit("Local deposition plane budget exceeded: local endpoints=" + std::to_string(boundaries.size()) +
                                    ", previous planes=" + std::to_string(result.sections) +
                                    ", indexed tiles=" + std::to_string(buckets.size()) + ".");
        }
        auto& intervals = current->intervals;
        std::vector<double> planes;
        std::set<coord_t> retained_endpoints;
        if (sampling_schedule) {
            intervals=sampling_schedule->intervals(index);
            for (const auto& interval : intervals) {
                if (checkpoint) checkpoint();
                const auto lo=interval.first,hi=interval.second;
                planes.push_back((limits.geometry.half_native_section_planes ? double(lo)+.5*double(hi-lo) :
                    double(lo+(hi-lo)/2))*1e-9);
            }
        } else for (;;) {
            intervals.clear();
            planes.clear();
            std::size_t planned = 0;
            for (auto it = boundaries.begin(), next = std::next(it); next != boundaries.end(); ++it, ++next) {
                const coord_t span  = *next - *it;
                const coord_t count = ((span + height - 1) / height) * interval_subdivisions;
                charge(planned, std::size_t(count), limits.max_sections - result.sections, "Local deposition plane budget exceeded.");
                charge(index_work, std::size_t(count), limits.max_index_work, "Local partition construction work budget exceeded.");
                if (statistics)
                    statistics->planned_sections = result.sections + planned;
                for (coord_t j = 0; j < count; ++j) {
                    if (checkpoint)
                        checkpoint();
                    const coord_t lo = *it + span * j / count, hi = *it + span * (j + 1) / count;
                    // Never discard an interval below native sampling resolution.
                    if (hi - lo < minimum_interval)
                        throw ResourceLimit("Local deposition interval " + std::to_string(lo) + ".." + std::to_string(hi) +
                                            " nm in tile (" + std::to_string(index[0]) + "," + std::to_string(index[1]) +
                                            ") is below native sampling resolution; indexed tiles=" + std::to_string(buckets.size()) +
                                            ", references=" + std::to_string(references) +
                                            ", previous planes=" + std::to_string(result.sections) +
                                            ", retained endpoints=" + std::to_string(retained_endpoints.count(lo)) + "/" +
                                            std::to_string(retained_endpoints.count(hi)) + ".");
                    intervals.emplace_back(lo, hi);
                    planes.push_back(
                        (limits.geometry.half_native_section_planes ? double(lo) + .5 * double(hi - lo) : double(lo + (hi - lo) / 2)) *
                        1e-9);
                }
            }
            if (!partition_envelopes)
                break;
            bool added = false, missed = false;
            // Keep exact endpoints locally when the coarser partition misses a
            // positive bead. Recheck all beads: inserting a cut changes nearby
            // midplanes. Every retry adds boundaries and consumes bounded work.
            for (const auto& bead : selected) {
                if (checkpoint)
                    checkpoint();
                charge(index_work, 1, limits.max_index_work, "Local partition coverage work budget exceeded.");
                const auto lo    = std::min(native(bead.start_m.z()) - native(bead.height_m),
                                            native(bead.end_m.z()) - native(bead.end_height_m.value_or(bead.height_m)));
                const auto hi    = std::max(native(bead.start_m.z()), native(bead.end_m.z()));
                const auto first = std::upper_bound(planes.begin(), planes.end(), double(lo) * 1e-9, [&](double value, double plane) {
                    charge(index_work, 1, limits.max_index_work, "Local partition coverage work budget exceeded.");
                    return value < plane;
                });
                if (first != planes.end() && *first < double(hi) * 1e-9)
                    continue;
                missed = true;
                retained_endpoints.insert(lo);
                retained_endpoints.insert(hi);
                const bool new_lo = boundaries.insert(lo).second, new_hi = boundaries.insert(hi).second;
                added = added || new_lo || new_hi;
                if (boundaries.size() > limits.max_sections - result.sections + 1)
                    throw ResourceLimit("Local bead coverage exceeds the section budget.");
            }
            if (!missed)
                break;
            if (!added)
                throw ResourceLimit("A positive deposited bead remains unresolved at native partition resolution.");
            // Quantized envelope cuts are quadrature choices, not deposited
            // boundaries. An exact retained endpoint can lie one native unit
            // from such a cut, leaving no integer midplane. Remove only that
            // redundant internal grid cut, never an exact retained endpoint or
            // the enclosing extremes. The next pass rechecks every bead.
            if (!limits.geometry.half_native_section_planes)
                for (const coord_t endpoint : retained_endpoints) {
                    if (checkpoint)
                        checkpoint();
                    charge(index_work, 2, limits.max_index_work, "Local partition coverage work budget exceeded.");
                    for (const coord_t neighbor : {endpoint - 1, endpoint + 1}) {
                        const auto cut = boundaries.find(neighbor);
                        if (cut != boundaries.end() && cut != boundaries.begin() && std::next(cut) != boundaries.end() &&
                            retained_endpoints.count(neighbor) == 0)
                            boundaries.erase(cut);
                    }
                }
        }
        charge(result.sections, planes.size(), limits.max_sections, "Local deposition plane budget exceeded.");
        if (schedule_only) return;
        if (planes.empty()) {
            current.reset();
            ++cursor;
            return;
        }
        packets->begin_tile(tile);
        // The sweep retains active-bead order and does not rescan the whole
        // project for each plane. Clip before material conflict reconciliation.
        SectionSweepLimits sweep;
        sweep.geometry           = limits.geometry;
        sweep.max_sections       = std::min<std::size_t>(10000, limits.max_sections);
        sweep.max_segment_visits = std::min<std::size_t>(20000000, limits.max_segment_visits - result.segment_visits);
        const auto resolve       = [this](std::size_t, double z) {
            current->reserved_ownership = limits.geometry.max_slab_group_visits;
            charge(result.ownership_work, current->reserved_ownership, limits.max_ownership_work,
                   "Local section ownership total work budget exceeded.");
            result.sampler_work.other+=current->reserved_ownership;
            return ownership ? ownership(z) : std::optional<SectionOwnership>{};
        };
        current->sweep = std::make_unique<DepositionSectionCursor>(selected, std::move(planes), bed_z, sweep, checkpoint,
                                                                   (ownership || limits.geometry.exact_per_command_composition) ?
                                                                       std::function<std::optional<SectionOwnership>(std::size_t, double)>(
                                                                           resolve) :
                                                                       std::function<std::optional<SectionOwnership>(std::size_t, double)>{},
                                                                   &current->clip);
    }
    bool packet(const LocalSectionPacketVisitor& visitor)
    {
        if (!current)
            prepare(visitor);
        if (!current)
            return true;
        auto& tile               = current->tile;
        auto& previous_count     = current->previous_count;
        auto& previous_top       = current->previous_top;
        auto& reserved_ownership = current->reserved_ownership;
        bool emitted             = false;
        const auto flush_packet  = [&] {
            if (!tile.slabs.empty()) {
                visitor.band(tile);
                tile.slabs.clear();
                emitted = true;
            }
        };
        current->sampling=true;
        while (current->sweep->advance([&](std::size_t plane, const DepositedSection& section) {
            current->sampled_plane=plane;
            if (checkpoint)
                checkpoint();
            if (ownership || limits.geometry.exact_per_command_composition) {
                if (section.ownership_work > reserved_ownership)
                    throw ResourceLimit("Local section ownership exceeded its reserved work.");
                result.ownership_work -= reserved_ownership - section.ownership_work;
                // All additions are bounded by the successful reservation
                // above; category sums never exceed max_ownership_work.
                result.sampler_work.other-=reserved_ownership-section.sampler_work.other;
                result.sampler_work.exact_construction+=section.sampler_work.exact_construction;
                result.sampler_work.exact_composition+=section.sampler_work.exact_composition;
                reserved_ownership = 0;
            }
            const auto [lo, hi]    = current->intervals[plane];
            const double thickness = double(hi - lo) * 1e-9;
            result.maximum_boundary_uncertainty_nm = std::max(result.maximum_boundary_uncertainty_nm, section.boundary_error.uncertainty_nm);
            const auto accumulate_bound = [&](double& total, double area_bound) {
                if (area_bound == 0.)
                    return;
                const double infinity = std::numeric_limits<double>::infinity();
                const double product  = std::nextafter(area_bound * std::nextafter(thickness, infinity), infinity);
                total                 = std::nextafter(total + product, infinity);
            };
            accumulate_bound(result.added_representation_volume_bound_m3, section.boundary_error.added_area_bound_m2);
            accumulate_bound(result.removed_representation_volume_bound_m3, section.boundary_error.removed_area_bound_m2);
            // Merge only an entire identical material partition. Merging
            // one material while another changes would create overlapping
            // Z intervals for the topology qualifier. Every sampling plane
            // is still evaluated; this only compacts the sampled solid.
            bool same = previous_count != 0 && previous_top == lo && previous_count == section.regions.size();
            for (std::size_t i = 0; same && i < previous_count; ++i) {
                if (checkpoint)
                    checkpoint();
                const auto& previous = tile.slabs[i];
                const auto& region   = section.regions[i];
                charge(index_work, region.footprint.vertex_count() + 1, limits.max_index_work, "Local section comparison work budget exceeded.");
                // Exact handles are retained, not compared by pointer identity
                // or replaced by native polygons to coalesce sampled bands.
                same = !previous.footprint.is_exact() && !region.footprint.is_exact() &&
                       previous.material == region.material && previous.boundary_uncertainty_nm == section.boundary_error.uncertainty_nm &&
                       previous.footprint.native() == region.footprint.native();
            }
            if (same) {
                for (std::size_t i = 0; i < previous_count; ++i)
                    tile.slabs[i].top_m = double(hi) * 1e-9;
            } else {
                flush_packet();
            }
            const bool has_exact=std::any_of(section.regions.begin(),section.regions.end(),
                [](const auto& region) { return region.footprint.is_exact(); });
            for (const auto& region : section.regions) {
                if (!same) {
                    charge(result.boundary_points, region.footprint.vertex_count(), limits.max_boundary_points,
                           "Local deposition cumulative geometry exceeds its point budget.");
                    tile.slabs.push_back(
                        {double(lo) * 1e-9, double(hi) * 1e-9, region.material, region.footprint,
                         section.boundary_error.uncertainty_nm, section.occupied_partition});
                    if (statistics) {
                        ++statistics->retained_slabs;
                        statistics->retained_points = result.boundary_points;
                    }
                }
                if (!has_exact)
                    for (const auto& polygon : region.footprint.native())
                        result.volume_m3 += polygon.area() * double(hi - lo) * 1e-27;
            }
            // Exact section area was integrated and charged by the sampler.
            // This metric never substitutes a rounded geometry for a slab.
            if (has_exact) result.volume_m3+=section.area_m2*thickness;
            previous_count = section.regions.size();
            previous_top   = hi;
            current->sampled_plane=plane+1;
        })) {
            if (emitted)
                return true;
        }
        current->sampling=false;
        result.segment_visits += current->sweep->completed_visits();
        flush_packet();
        visitor.end_tile();
        if (statistics)
            ++statistics->completed_tiles;
        current.reset();
        ++cursor;
        return true;
    }
};

struct LocalSectionSamplingSchedule::Impl {
    std::vector<const std::vector<DepositedSegment>*> sources;
    std::vector<DepositedSegment> combined;
    std::unique_ptr<LocalSectionCursor> planner;
    std::optional<std::array<int,2>> last;
    std::vector<std::pair<coord_t,coord_t>> cached;
    bool failed=false,busy=false;
};
LocalSectionSamplingSchedule::LocalSectionSamplingSchedule(
    const std::vector<const std::vector<DepositedSegment>*>& sources,const Eigen::Vector2d& origin,
    const Eigen::Vector2d& size,double bed,double height,const LocalSectionLimits& limits,
    const std::function<void()>& checkpoint,unsigned subdivisions) : m_impl(std::make_unique<Impl>())
{
    if (sources.empty() || sources.size()>257) throw std::invalid_argument("Invalid shared sampling sources.");
    const auto maximum=std::min<std::size_t>(2000000,limits.geometry.max_segments);
    for (const auto* source : sources) {
        if (checkpoint) checkpoint();
        if (!source || source->empty()) throw std::invalid_argument("Missing shared sampling source.");
        if (source->size()>maximum-m_impl->combined.size()) throw ResourceLimit("Shared sampling source budget exceeded.");
        for (const auto& bead : *source) {
            if (checkpoint) checkpoint();
            m_impl->combined.push_back(bead);
        }
        m_impl->sources.push_back(source);
    }
    m_impl->planner=std::make_unique<LocalSectionCursor>(m_impl->combined,origin,size,bed,height,limits,
        std::function<std::optional<SectionOwnership>(double)>{},checkpoint,subdivisions);
    auto& planner=*m_impl->planner->m_impl;
    charge(planner.index_work,m_impl->combined.size(),limits.max_index_work,"Shared sampling copy work exceeded.");
}
LocalSectionSamplingSchedule::~LocalSectionSamplingSchedule() = default;
LocalSectionCursor::IndexUsage LocalSectionSamplingSchedule::index_usage() const
{
    if (m_impl->failed || m_impl->busy) throw std::logic_error("Inactive shared sampling schedule.");
    return m_impl->planner->index_usage();
}
void LocalSectionSamplingSchedule::validate(const std::vector<DepositedSegment>& source,const Eigen::Vector2d& origin,
    const Eigen::Vector2d& size,double bed,double height,const LocalSectionLimits& limits,unsigned subdivisions) const
{
    const auto& planner=*m_impl->planner->m_impl;
    if (m_impl->failed || m_impl->busy || std::find(m_impl->sources.begin(),m_impl->sources.end(),&source)==m_impl->sources.end() ||
        planner.o!=std::array<coord_t,2>{native(origin.x()),native(origin.y())} ||
        planner.s!=std::array<coord_t,2>{native(size.x()),native(size.y())} || planner.bed!=native(bed) ||
        planner.height!=native(height) || planner.interval_subdivisions!=subdivisions ||
        planner.limits.partition_increment_m!=limits.partition_increment_m ||
        planner.limits.geometry.half_native_section_planes!=limits.geometry.half_native_section_planes)
        throw std::invalid_argument("Shared sampling source or quadrature settings changed.");
}
const std::vector<std::pair<coord_t,coord_t>>& LocalSectionSamplingSchedule::intervals(const std::array<int,2>& tile)
{
    if (m_impl->failed || m_impl->busy) {
        m_impl->failed=true;
        throw std::logic_error("Inactive or reentered shared sampling schedule.");
    }
    m_impl->busy=true;
    try {
        auto& planner=*m_impl->planner->m_impl;
        if (planner.checkpoint) planner.checkpoint();
        charge(planner.index_work,1,planner.limits.max_index_work,"Shared sampling lookup budget exceeded.");
        if (m_impl->last && tile<*m_impl->last) throw std::invalid_argument("Shared sampling tile order reversed.");
        if (!m_impl->last || tile!=*m_impl->last) {
            planner.cursor=planner.buckets.find(tile);
            if (planner.cursor==planner.buckets.end()) throw std::invalid_argument("Shared sampling tile is missing.");
            planner.prepare({},true);
            m_impl->cached=std::move(planner.current->intervals);
            planner.current.reset();
            m_impl->last=tile;
        }
        if (m_impl->failed) throw std::logic_error("Reentered shared sampling schedule.");
        m_impl->busy=false;
        return m_impl->cached;
    } catch (...) {
        m_impl->busy=false;
        m_impl->failed=true;
        throw;
    }
}

LocalSectionCursor::LocalSectionCursor(const std::vector<DepositedSegment>& source,
                                       const Eigen::Vector2d& origin,
                                       const Eigen::Vector2d& size,
                                       double bed,
                                       double height,
                                       const LocalSectionLimits& limits,
                                       std::function<std::optional<SectionOwnership>(double)> ownership,
                                       std::function<void()> checkpoint,
                                       unsigned subdivisions,
                                       LocalSectionStatistics* statistics,
                                       const std::vector<DepositedSegment>* partitions,
                                       std::shared_ptr<LocalSectionSamplingSchedule> schedule)
    : m_impl(std::make_unique<Impl>(
          source, origin, size, bed, height, limits, std::move(ownership), std::move(checkpoint), subdivisions, statistics, partitions))
{
    if (schedule) {
        if (partitions) throw std::invalid_argument("Shared section schedules do not accept independent partition envelopes.");
        schedule->validate(source,origin,size,bed,height,limits,subdivisions);
        m_impl->sampling_schedule=std::move(schedule);
    }
}
LocalSectionCursor::~LocalSectionCursor() = default;
std::optional<std::array<int,2>> LocalSectionCursor::next_tile_index() const
{
    if (m_impl->failed || m_impl->active) { m_impl->failed = true; throw std::logic_error("Inactive or reentered section cursor."); }
    if (m_impl->cursor == m_impl->buckets.end()) return {};
    return m_impl->cursor->first;
}
const LocalSectionMetrics& LocalSectionCursor::metrics() const
{
    if (m_impl->failed || m_impl->active) { m_impl->failed = true; throw std::logic_error("Inactive or reentered section cursor."); }
    return m_impl->result;
}
LocalSectionCursor::IndexUsage LocalSectionCursor::index_usage() const
{
    if (m_impl->failed || m_impl->active) { m_impl->failed = true; throw std::logic_error("Inactive or reentered section cursor."); }
    return {m_impl->buckets.size(),m_impl->references};
}
bool LocalSectionCursor::advance(const LocalSectionPacketVisitor& visitor)
{
    if (m_impl->current) {
        m_impl->failed = true;
        throw std::logic_error("Whole-tile advance cannot interrupt a packet traversal.");
    }
    const auto index = next_tile_index();
    if (!advance_packet(visitor)) return false;
    while (next_tile_index() == index) advance_packet(visitor);
    return true;
}
bool LocalSectionCursor::advance_packet(const LocalSectionPacketVisitor& visitor)
{
    if (m_impl->failed || m_impl->active || !visitor.begin_tile || !visitor.band || !visitor.end_tile) {
        m_impl->failed = true;
        throw std::invalid_argument("Invalid or reentered section cursor advance.");
    }
    if (m_impl->cursor == m_impl->buckets.end()) return false;
    m_impl->active = true;
    try {
        m_impl->packet(visitor);
        if (m_impl->failed) throw std::logic_error("Reentered section cursor.");
        m_impl->active = false;
        return true;
    } catch (const ResourceLimit& error) {
        m_impl->active=false;
        m_impl->failed=true;
        std::string context;
        if (m_impl->cursor!=m_impl->buckets.end()) {
            const auto& index=m_impl->cursor->first;
            context=" Local section cursor tile=("+std::to_string(index[0])+","+std::to_string(index[1])+").";
        }
        if (m_impl->current) {
            const auto& state=*m_impl->current;
            const auto& tile=state.tile;
            context+=" Bounds_nm=["+std::to_string(tile.minimum[0])+","+std::to_string(tile.minimum[1])+","+
                std::to_string(tile.maximum[0])+","+std::to_string(tile.maximum[1])+"].";
            if (state.sampling && state.sampled_plane<state.intervals.size()) {
                const auto interval=state.intervals[state.sampled_plane];
                context+=" Sample_interval_nm=["+std::to_string(interval.first)+","+std::to_string(interval.second)+"].";
            }
        }
        const auto& categories=m_impl->result.sampler_work;
        context+=" Sampler_work={construction="+std::to_string(categories.exact_construction)+
            ", composition="+std::to_string(categories.exact_composition)+
            ", other_including_reservation="+std::to_string(categories.other)+
            ", total="+std::to_string(m_impl->result.ownership_work)+"}.";
        throw ResourceLimit(std::string(error.what())+context);
    } catch (...) { m_impl->active = false; m_impl->failed = true; throw; }
}

static LocalSections build_local_sections(const std::vector<DepositedSegment>& segments,
    const Eigen::Vector2d& origin, const Eigen::Vector2d& size,
    double bed_z, double maximum_height, const LocalSectionLimits& limits,
    const std::function<std::optional<SectionOwnership>(double)>& ownership,
    const std::function<void()>& checkpoint, unsigned interval_subdivisions, LocalSectionStatistics* statistics,
    const std::vector<DepositedSegment>* partition_envelopes,
    const std::function<void(const LocalSectionTile&)>& visit,
    const LocalSectionPacketVisitor* packets)
{
    if (!visit && !packets && (limits.max_boundary_points > 8000000 || limits.max_index_work > 20000000))
        throw std::invalid_argument("Invalid retained local section limits.");
    LocalSectionCursor cursor(segments,origin,size,bed_z,maximum_height,limits,ownership,checkpoint,
        interval_subdivisions,statistics,partition_envelopes);
    LocalSections result;
    LocalSectionTile retained;
    LocalSectionPacketVisitor collector;
    collector.begin_tile = [&](const LocalSectionTile& header) { retained = header; };
    collector.band = [&](const LocalSectionTile& band) {
        retained.slabs.insert(retained.slabs.end(),band.slabs.begin(),band.slabs.end());
    };
    collector.end_tile = [&] {
        if (!retained.slabs.empty()) {
            if (visit) visit(retained);
            else result.tiles.push_back(std::move(retained));
        }
    };
    while (cursor.advance(packets ? *packets : collector)) {}
    static_cast<LocalSectionMetrics&>(result) = cursor.metrics();
    return result;
}

LocalSections local_deposition_sections(const std::vector<DepositedSegment>& segments,
    const Eigen::Vector2d& origin, const Eigen::Vector2d& size,
    double bed_z, double maximum_height, const LocalSectionLimits& limits,
    const std::function<std::optional<SectionOwnership>(double)>& ownership,
    const std::function<void()>& checkpoint, unsigned interval_subdivisions, LocalSectionStatistics* statistics,
    const std::vector<DepositedSegment>* partition_envelopes)
{
    return build_local_sections(segments, origin, size, bed_z, maximum_height, limits,
        ownership, checkpoint, interval_subdivisions, statistics, partition_envelopes, {}, nullptr);
}

LocalSectionMetrics visit_local_deposition_sections(const std::vector<DepositedSegment>& segments,
    const Eigen::Vector2d& origin, const Eigen::Vector2d& size, double bed_z, double maximum_height,
    const std::function<void(const LocalSectionTile&)>& visit, const LocalSectionLimits& limits,
    const std::function<std::optional<SectionOwnership>(double)>& ownership,
    const std::function<void()>& checkpoint, unsigned interval_subdivisions, LocalSectionStatistics* statistics,
    const std::vector<DepositedSegment>* partition_envelopes)
{
    if (!visit) throw std::invalid_argument("Local section traversal requires a tile consumer.");
    return build_local_sections(segments, origin, size, bed_z, maximum_height, limits,
        ownership, checkpoint, interval_subdivisions, statistics, partition_envelopes, visit, nullptr);
}

LocalSectionMetrics visit_local_deposition_packets(const std::vector<DepositedSegment>& segments,
    const Eigen::Vector2d& origin, const Eigen::Vector2d& size, double bed_z, double maximum_height,
    const LocalSectionPacketVisitor& visit, const LocalSectionLimits& limits,
    const std::function<std::optional<SectionOwnership>(double)>& ownership,
    const std::function<void()>& checkpoint, unsigned interval_subdivisions, LocalSectionStatistics* statistics,
    const std::vector<DepositedSegment>* partition_envelopes)
{
    if (!visit.begin_tile || !visit.band || !visit.end_tile)
        throw std::invalid_argument("Local section packet traversal requires complete callbacks.");
    return build_local_sections(segments, origin, size, bed_z, maximum_height, limits,
        ownership, checkpoint, interval_subdivisions, statistics, partition_envelopes, {}, &visit);
}
} // namespace Slic3r::Thermoelastic
