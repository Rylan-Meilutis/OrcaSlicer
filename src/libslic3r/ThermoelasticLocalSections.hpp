#pragma once

#include "ThermoelasticLayers.hpp"
#include "ThermoelasticLocalTopology.hpp"
#include "ThermoelasticSectionStudy.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include "ThermoelasticSections.hpp"
#include "libslic3r.h"
#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {

struct LocalSectionTile {
    std::array<int, 2> index;
    // Native nanometre coordinates. Slabs are clipped to this closed box;
    // neighboring boxes share only a face, never positive volume.
    std::array<coord_t, 2> minimum, maximum;
    std::vector<DepositedSlab> slabs;
};

struct LocalSectionLimits {
    LayerGeometryLimits geometry;
    // Explicit analysis-only quadrature coarsening (zero disables). Original
    // beads are still sampled; independently refine before using results.
    // Mutually exclusive with caller-supplied partition envelopes.
    double partition_increment_m = 0.;
    // Offline fine-tile studies may explicitly request up to16384 buckets.
    // Independent reference/point/work limits still bound their storage.
    std::size_t max_tiles = 4096, max_tile_references = 8000000;
    std::size_t max_sections = 200000, max_segment_visits = 80000000;
    std::size_t max_boundary_points = 8000000, max_index_work = 20000000;
    std::size_t max_ownership_work = 1000000000;
};

struct LocalSectionMetrics {
    std::size_t sections = 0, segment_visits = 0, boundary_points = 0;
    std::size_t ownership_work = 0;
    // Sum equals ownership_work. During a failed/in-progress sample, other
    // includes its unrefunded work reservation, not measured partial work.
    SectionSamplerWork sampler_work;
    double volume_m3 = 0.;
    unsigned maximum_boundary_uncertainty_nm = 0;
    double added_representation_volume_bound_m3 = 0., removed_representation_volume_bound_m3 = 0.;
    bool piecewise_constant = false;
    bool uses_partition_envelopes = false;
};
struct LocalSections : LocalSectionMetrics {
    std::vector<LocalSectionTile> tiles;
};
struct LocalSectionStatistics {
    // Admission diagnostics only; partial counters never certify geometry.
    std::size_t indexed_tiles = 0, completed_tiles = 0, retained_slabs = 0, planned_sections = 0, retained_points = 0;
};

// Analysis quadrature boundaries only. These envelopes must never replace the
// original beads passed to the sampler or be treated as deposited geometry.
struct SectionPartitionEnvelopes {
    std::vector<DepositedSegment> segments;
    double maximum_shift_m = 0.;
    std::size_t collapsed_intervals = 0, expanded_zero_endpoints = 0;
};
SectionPartitionEnvelopes make_section_partition_envelopes(const std::vector<DepositedSegment>& source,
    double bed_z_m, double increment_m, const std::function<void()>& checkpoint = {});

// Tile-local midplane quadrature of the actual nominal bead sections. Every
// local bead endpoint top/bottom is retained, including scarf ramps. A
// sampling interval never crosses a local discontinuity or exceeds maximum_height_m.
// Adjacent intervals with exactly identical full material partitions and
// uncertainty are stored as one prism without changing the sampled solid.
// Opt-in exact canonical partitions retain their individual sampled intervals;
// they are not coalesced through native conversion or handle identity.
// Empty intervals, changed materials and changed footprints are never bridged.
// The result is an UNQUALIFIED quadrature geometry: callers must independently refine it
// and certify physical connectivity before using it for mechanical analysis.
// Material ownership, if supplied, is consulted at the actual sampled plane.
// Bounds are enforced across all tiles, not restarted for each tile.
// Optional partition envelopes supply ONLY Z boundaries; they must enclose the
// original endpoint intervals without new bed contact. The section sweep always
// evaluates the original beads. Missed positive beads add their original local
// endpoints to the partition. Internal approximate cuts one native unit from
// these retained endpoints are removed; retained endpoints and the enclosing
// extremes are never removed. Every bead is rechecked after repartitioning.
// With explicit half-native section planes, one-nanometre intervals retain
// their exact interior half-nanometre plane and no adjacent cut is removed.
// Budget or native-resolution exhaustion rejects
// the study instead of silently dropping material.
LocalSections local_deposition_sections(const std::vector<DepositedSegment>& segments,
    const Eigen::Vector2d& origin_m, const Eigen::Vector2d& tile_size_m,
    double bed_z_m, double maximum_height_m, const LocalSectionLimits& limits = {},
    const std::function<std::optional<SectionOwnership>(double)>& ownership = {},
    const std::function<void()>& checkpoint = {}, unsigned interval_subdivisions = 1,
    LocalSectionStatistics* statistics = nullptr,
    const std::vector<DepositedSegment>* partition_envelopes = nullptr);

// Synchronous bounded alternative: each complete nonempty tile is borrowed by
// the consumer, then released before the next tile. Metrics are returned only
// after the entire traversal succeeds; callbacks/cancellation propagate.
// The consumer must not retain references to the tile. The same cumulative
// geometry/work defaults apply, and emitted tiles are not a topology certificate.
// Offline callers may explicitly allow up to 128 million cumulative boundary
// points for either visitor form; the retained form stays capped at eight million.
// Consumers must separately bound their live payload and accumulated metadata.
LocalSectionMetrics visit_local_deposition_sections(const std::vector<DepositedSegment>& segments,
    const Eigen::Vector2d& origin_m, const Eigen::Vector2d& tile_size_m,
    double bed_z_m, double maximum_height_m,
    const std::function<void(const LocalSectionTile&)>& visit,
    const LocalSectionLimits& limits = {},
    const std::function<std::optional<SectionOwnership>(double)>& ownership = {},
    const std::function<void()>& checkpoint = {}, unsigned interval_subdivisions = 1,
    LocalSectionStatistics* statistics = nullptr,
    const std::vector<DepositedSegment>* partition_envelopes = nullptr);

struct LocalSectionPacketVisitor {
    // Header has no slabs. Each band borrows exactly one complete coalesced
    // material partition; bands increase in Z and may be separated by gaps.
    std::function<void(const LocalSectionTile&)> begin_tile, band;
    std::function<void()> end_tile;
};
// Resumable lexicographic sampler. Indexes the source once and advances one
// tile at a time, allowing two snapshots to be compared without retaining
// either whole polygon history or rebuilding an index for every tile.
// Borrows segments, partition_envelopes and statistics: they must outlive the
// cursor and remain unchanged. Budgets accumulate across all advances.
class LocalSectionSamplingSchedule;
class LocalSectionCursor {
public:
    struct IndexUsage { std::size_t tiles = 0, references = 0; };
    LocalSectionCursor(const std::vector<DepositedSegment>& segments,
        const Eigen::Vector2d& origin_m, const Eigen::Vector2d& tile_size_m,
        double bed_z_m, double maximum_height_m, const LocalSectionLimits& limits = {},
        std::function<std::optional<SectionOwnership>(double)> ownership = {},
        std::function<void()> checkpoint = {}, unsigned interval_subdivisions = 1,
        LocalSectionStatistics* statistics = nullptr,
        const std::vector<DepositedSegment>* partition_envelopes = nullptr,
        std::shared_ptr<LocalSectionSamplingSchedule> sampling_schedule = {});
    ~LocalSectionCursor();
    LocalSectionCursor(const LocalSectionCursor&) = delete;
    LocalSectionCursor& operator=(const LocalSectionCursor&) = delete;
    std::optional<std::array<int,2>> next_tile_index() const;
    bool advance(const LocalSectionPacketVisitor& visitor);
    // Advances at most one coalesced height band, retaining only its pending
    // successor. May also emit begin/end events; keep the same consumer for
    // the entire tile. Mixing with advance is allowed at a tile boundary only.
    bool advance_packet(const LocalSectionPacketVisitor& visitor);
    // Final only after exhaustion; partial metrics are admission diagnostics.
    const LocalSectionMetrics& metrics() const;
    IndexUsage index_usage() const;
private:
    friend class LocalSectionSamplingSchedule;
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};
// Shared quadrature only, never shared material geometry. All registered source
// vectors must remain unchanged and outlive the schedule and its cursors. One
// bounded index covers their union; only the latest lexicographic tile's native
// intervals are cached. Every participating cursor still samples its own beads.
// Scope one synchronized activation traversal: independently nested history
// traversals need separate schedules, even when their source vectors overlap.
// Separate schedules do not certify geometric nesting across activation times.
class LocalSectionSamplingSchedule {
public:
    LocalSectionSamplingSchedule(const std::vector<const std::vector<DepositedSegment>*>& sources,
        const Eigen::Vector2d& origin_m,const Eigen::Vector2d& tile_size_m,double bed_z_m,double maximum_height_m,
        const LocalSectionLimits& limits = {},const std::function<void()>& checkpoint = {},unsigned subdivisions = 1);
    ~LocalSectionSamplingSchedule();
    LocalSectionCursor::IndexUsage index_usage() const;
    const std::vector<std::pair<coord_t,coord_t>>& intervals(const std::array<int,2>& tile);
private:
    friend class LocalSectionCursor;
    void validate(const std::vector<DepositedSegment>& source,const Eigen::Vector2d& origin,
        const Eigen::Vector2d& size,double bed,double height,const LocalSectionLimits& limits,unsigned subdivisions) const;
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};
// Band-packet form of the same sampler. Only the current coalesced partition
// is retained, not an entire tile's history. All callbacks are required; failure
// propagates without an end event or successful metrics for a partial tile.
LocalSectionMetrics visit_local_deposition_packets(const std::vector<DepositedSegment>& segments,
    const Eigen::Vector2d& origin_m, const Eigen::Vector2d& tile_size_m,
    double bed_z_m, double maximum_height_m, const LocalSectionPacketVisitor& visit,
    const LocalSectionLimits& limits = {},
    const std::function<std::optional<SectionOwnership>(double)>& ownership = {},
    const std::function<void()>& checkpoint = {}, unsigned interval_subdivisions = 1,
    LocalSectionStatistics* statistics = nullptr,
    const std::vector<DepositedSegment>* partition_envelopes = nullptr);

SectionTopology qualify_local_section_topology(const std::vector<LocalSectionTile>& tiles,
    double bed_z_m, const SectionTopologyLimits& limits = {}, const std::function<void()>& checkpoint = {});

struct StreamIntegrationLimits {
    // Explicit offline streaming only; defaults stay modest. Hard cumulative
    // Without sparse-factor opt-in, cumulative ceilings are 600M work,32M
    // points,250k pieces/labels. Sparse-factor mode allows 3.2B work,64M
    // points,500k pieces; labels and live point storage retain their bounds.
    std::size_t max_work = 20000000, max_total_points = 2000000;
    std::size_t max_live_points = 200000, max_groups = 31250, max_memberships = 2000000;
    std::size_t max_pieces = 100000, max_provisional_components = 100000;
};
// Admission diagnostics, including partial traversal on failure; never a field
// or a geometry certificate. Counters reset at entry.
struct StreamIntegrationStatistics {
    std::size_t processed_points = 0, integration_work = 0, assembly_work = 0;
    std::size_t peak_groups = 0, peak_memberships = 0, peak_input_points = 0, peak_clip_points = 0;
    std::size_t piece_count = 0, provisional_mapping_size = 0;
};
struct LocalSectionStreamSettings {
    Eigen::Vector3d origin_m = Eigen::Vector3d::Zero();
    Eigen::Vector3d basis_spacing_m = Eigen::Vector3d::Zero();
    // Optional explicit Q1 Z planes (metres), including bed and entire solid.
    // Empty keeps the regular grid; otherwise replaces only basis_spacing_m.z().
    // Strictly increasing native-nanometre coordinates, 2..128 planes.
    std::vector<double> basis_z_planes_m;
    std::vector<SectionBedReleaseRegion> released_bed_regions;
    std::optional<SectionCellContact> cell_contact;
    SectionLoadField load;
    Eigen::Vector2d tile_size_m = Eigen::Vector2d::Zero();
    double bed_z_m = 0., maximum_height_m = 0.;
    unsigned interval_subdivisions = 1;
    LocalSectionLimits geometry;
    SectionSolveLimits mechanics;
    SectionTopologyLimits topology;
    // Explicit streaming-only cumulative budget. Live geometry/moment storage
    // and sparse solve bounds remain independent; absent preserves legacy caps.
    std::optional<StreamIntegrationLimits> integration;
    std::optional<StreamTopologyLimits> streaming_topology;
};
struct LocalSectionMechanicalResult {
    LocalSectionMetrics geometry;
    SectionMechanicalResult field;
};
// Closed complete sampler -> connectivity -> moments -> sparse solve pipeline.
// Returns a SINGLE fixed-basis/quadrature result, still UNQUALIFIED: independent
// refinement and physical material/interface laws remain necessary. No partial
// metrics or field escape failed traversal, ownership, connectivity or solve.
// Existing cumulative geometry, topology and mechanics budgets remain separate
// and unchanged; streaming does not itself authorize larger inputs.
LocalSectionMechanicalResult solve_local_deposition_sections(const std::vector<DepositedSegment>& segments,
    const std::vector<SectionMaterial>& materials, const LocalSectionStreamSettings& settings,
    const std::function<std::optional<SectionOwnership>(double)>& ownership = {},
    const std::function<void()>& checkpoint = {}, const std::vector<DepositedSegment>* partition_envelopes = nullptr,
    const std::vector<SectionSampleRequest>& samples = {}, StreamIntegrationStatistics* statistics = nullptr,
    StreamTopologyStatistics* topology_statistics = nullptr, LocalSectionStatistics* geometry_statistics = nullptr);

// Four streamed solves on identical input/partition: reference, then X/Y/Z
// halvings. Uses the same selector as the retained topology probe. Only scalar
// diagnostics survive each run; no field or convergence verdict is returned.
// All directions must finish before a recommendation is made. Up to four of
// each configured per-solve geometry/topology/integration/mechanics budget are
// consumed. Callbacks and ownership must be deterministic across the four runs.
SectionDirectionalProbe probe_local_deposition_refinement(const std::vector<DepositedSegment>& segments,
    const std::vector<SectionMaterial>& materials, const LocalSectionStreamSettings& settings,
    const std::function<std::optional<SectionOwnership>(double)>& ownership = {},
    const std::function<void(unsigned)>& progress = {}, const std::function<void()>& checkpoint = {},
    const std::vector<DepositedSegment>* partition_envelopes = nullptr);

// Actual mechanical study over the local quadrature. Requires independently
// converged Q1 bases and two successive stable geometry comparisons, except
// exactly piecewise-constant rectangular geometry. No physical validation or
// bed-adhesion prediction is implied by numerical convergence.
DepositionSectionStudy study_local_deposition_sections(const std::vector<DepositedSegment>& segments,
    const std::vector<SectionMaterial>& materials, const Eigen::Vector3d& origin_m,
    const Eigen::Vector3d& basis_spacing_m, const Eigen::Vector2d& tile_size_m,
    double bed_z_m, const SectionStudySettings& settings,
    const std::function<std::optional<SectionOwnership>(double)>& ownership = {},
    const std::function<void(unsigned, unsigned)>& progress = {},
    const std::function<void()>& checkpoint = {}, const std::vector<SectionSampleRequest>& samples = {});

} // namespace Slic3r::Thermoelastic
