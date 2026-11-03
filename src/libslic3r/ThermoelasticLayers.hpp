#pragma once

#include "ExPolygon.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticBoundary.hpp"
#include "ThermoelasticSectionFootprint.hpp"
#include "ThermoelasticExactSections.hpp"
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace Slic3r::Thermoelastic {

class OccupiedSectionPartition;
class AuthenticatedDepositedPartition;

class LayerMaterialOverlap : public std::invalid_argument {
public:
    LayerMaterialOverlap(double bottom, double top, double area, std::size_t slot, const std::string& detail = {})
        : std::invalid_argument("Different materials overlap in the nominal deposited envelopes: bottom_nm=" +
              std::to_string(bottom*1e9) + ", top_nm=" + std::to_string(top*1e9) +
              ", area_nm2=" + std::to_string(area*1e18) + ", material=" + std::to_string(slot) + "." + detail),
          bottom_m(bottom), top_m(top), area_m2(area), material(slot) {}
    double bottom_m, top_m, area_m2;
    std::size_t material;
    // Optional bounded diagnostic mask at bottom_m == top_m. Never an
    // ownership decision. Retains complete polygons up to 4096 points.
    ExPolygons section_footprint;
    bool section_footprint_complete = false;
};

struct DepositedSlab {
    double bottom_m = 0., top_m = 0.;
    std::size_t material = 0;
    // Canonical native or exact geometry, including holes. No safety offset.
    SectionFootprint footprint;
    unsigned boundary_uncertainty_nm = 0;
    // Valid only for the complete, unchanged sampled material partition.
    // Subsets and material-specific activation must not inherit this proof.
    std::shared_ptr<const OccupiedSectionPartition> occupied_partition;
};
struct DepositedLayers {
    std::vector<DepositedSlab> slabs;
    double occupied_volume_m3 = 0.;
    std::size_t boundary_points = 0;
};
struct LayerGeometryLimits {
    std::size_t max_segments = 2000000, max_groups = 10000;
    std::size_t max_boundary_points = 2000000, max_slab_group_visits = 2000000;
    // Explicit group-work ceiling may reach8M for interface + topology studies;
    // the default remains2M and vertex/storage limits are independent.
    double arc_tolerance_m = 1e-6;
    // Opt-in Boolean representation band in native 1nm units; zero is strict.
    // Maximum4. Affects analysis only, never print paths or physical laws.
    unsigned ownership_boundary_grid_units = 0;
    // Opt-in query-plane precision only: native bead endpoints/heights and XY
    // remain on the 1 nm grid. Enables interior sampling of a 1 nm Z interval.
    bool half_native_section_planes = false;
    // Opt-in canonical storage for authenticated interface sensitivity only.
    // Preserves rational boundaries; does not reduce inherited uncertainty or
    // qualify contacts. Unported geometry consumers reject exact footprints.
    bool retain_exact_interface_sections = false;
    // Geometry-only opt-in for a single material: exact union of the same
    // normalized quantized command cuts. No ownership/qualification token.
    // Requires retained exact sections; default/native composition is unchanged.
    bool exact_per_command_composition = false;
};

// Compact snapshot of nominal horizontal rectangular-envelope beads. Overlap
// within a material is a union; overlap between materials is unresolved and
// rejected. Z intervals are split at every bead boundary, not rounded to a
// layer number. Coordinate quantization is 1 nm; circular caps are polygonal
// at the explicit arc tolerance. This is geometry, not a thermal history or a
// mechanically qualified coarse mesh. Sloped/rounded-section beads reject.
DepositedLayers layered_deposition(const std::vector<DepositedSegment> &segments,
    double bed_z_m, const LayerGeometryLimits &limits = {},
    const std::function<void()> &checkpoint = {});

// Disjoint categories of the existing sampler ownership-work charge, not
// timings or additional budgets. Native construction denotes exactification
// and its input visit charge; composition denotes material-union overlays.
struct SectionSamplerWork {
    std::size_t exact_construction = 0, exact_composition = 0, other = 0;
};
struct DepositedSection {
    struct Region {
        std::size_t material;
        SectionFootprint footprint;
    };
    std::vector<Region> regions;
    double area_m2 = 0.;
    std::size_t boundary_points = 0;
    double resolved_overlap_area_m2 = 0.;
    // Polygon-work charges consumed by optional ownership reconciliation.
    std::size_t ownership_work = 0;
    SectionSamplerWork sampler_work;
    SectionBoundaryError boundary_error;
    std::shared_ptr<const OccupiedSectionPartition> occupied_partition;
    std::shared_ptr<const AuthenticatedDepositedPartition> deposited_partition;
};
class SampledSectionConflict;
class AuthenticatedSectionPartition;
class NativeOwnershipPartitionFactory;
class DepositedOwnershipPartitionFactory;
class AuthenticatedDepositedPartition;
struct SectionOwnership {
    double z_m;
    // Intended native material ownership, not an adhesion/bonding model.
    std::vector<DepositedSection::Region> regions;
    // Optional deferred evidence: receives the complete actual (possibly
    // tile-clipped) conflict, never just a requested Z plane. Cannot accompany
    // precomputed regions. Its returned masks undergo the same union and
    // unique-owner validation; this does not authorize invented deposition.
    std::function<std::vector<DepositedSection::Region>(const LayerMaterialOverlap&)> resolve_conflict = {};
    std::function<std::vector<DepositedSection::Region>(const SampledSectionConflict&)> resolve_sampled_conflict = {};
    // Analysis-only full replacement, minted only by the generation-verified
    // native resolver. Its live sampler identity cannot be reused elsewhere.
    std::function<std::shared_ptr<const AuthenticatedSectionPartition>(const SampledSectionConflict&)> resolve_partition = {};
    std::shared_ptr<const AuthenticatedSectionPartition> partition;
    // Separate chronological material assumption, never native ownership.
    std::function<std::shared_ptr<const AuthenticatedDepositedPartition>(const SampledSectionConflict&)> resolve_deposited_partition = {};
    std::shared_ptr<const AuthenticatedDepositedPartition> deposited_partition;
};
// Horizontal section through nominal beads with linearly varying top Z
// and optional linearly varying design thickness. End caps retain endpoint Z
// and thickness; interior cuts never acquire new round caps.
// Unlike layered_deposition, sloped paths are retained rather than flattened.
// The bottom is included and the top excluded. This defines a nominal shape,
// not a reconstruction of variable flow from display metadata.
// Query planes use the native 1 nm grid by default. Explicit half-native mode
// evaluates 0.5 nm planes without changing the native bead or XY coordinates;
// direct sampling, sweep inclusion and ownership plane identity share the mode.
// Explicit RoundedRectangle beads retain the semicircular-side corner voids.
// Sloped rounded beads use linear top Z, vertical height and original endpoint
// half-caps. Tapered sides are ruled normal cross-sections, not the Minkowski
// union of variable-radius disks. This geometry does not infer extrusion flow.
// A section may be empty. Without ownership, material overlap rejects. Explicit
// ownership partitions conflicts only: each conflict point must have one owner
// with a deposited envelope there. Unambiguous beads outside all masks remain.
// The deposited union/holes must remain unchanged; this is intended ownership,
// not measured redistribution of extrusion volume or proof of interface bonding.
// Optional diagnostics observe raw exact material-union preflight only, not
// ownership reconciliation. The caller owns the sink and serializes access.
DepositedSection section_deposition(const std::vector<DepositedSegment> &segments,
    double z_m, double bed_z_m, const LayerGeometryLimits &limits = {},
    const std::function<void()> &checkpoint = {}, const SectionOwnership *ownership = nullptr,
    const ExPolygon *clip = nullptr, ExactOverlayDiagnostics* diagnostics = nullptr);

// Synchronous, noncopyable evidence minted by the actual section operation.
// The receiver must still bind/verify the full source prefix and generation;
// an arbitrary externally supplied LayerMaterialOverlap is not this evidence.
class SampledSectionConflict {
public:
    SampledSectionConflict(const SampledSectionConflict&) = delete;
    SampledSectionConflict& operator=(const SampledSectionConflict&) = delete;
    const LayerMaterialOverlap& conflict() const { return m_conflict; }
    const std::vector<DepositedSegment>& segments() const { return m_segments; }
    const LayerGeometryLimits& geometry_limits() const { return m_limits; }
    double bed_z_m() const { return m_bed; }
    const ExPolygon* clip() const { return m_clip; }
    const std::vector<DepositedSection::Region>& original_regions() const { return m_regions; }
    // Entry-frozen deposited inputs only, in segments() order. Zero means an
    // empty clipped cut; nonzero is its admitted normalized native point count.
    // No source authentication is implied. Other evidence requires recutting.
    const std::vector<std::size_t>* frozen_cut_boundary_points() const {
        return m_frozen_inputs ? m_cut_boundary_points : nullptr;
    }
private:
    friend class NativeOwnershipPartitionFactory;
    friend class DepositedOwnershipPartitionFactory;
    friend DepositedSection section_deposition(const std::vector<DepositedSegment>&, double, double,
        const LayerGeometryLimits&, const std::function<void()>&, const SectionOwnership*, const ExPolygon*, ExactOverlayDiagnostics*);
    SampledSectionConflict(const LayerMaterialOverlap& conflict, const std::vector<DepositedSegment>& segments,
        double bed, const LayerGeometryLimits& limits, const ExPolygon* clip,
        const std::vector<DepositedSection::Region>& regions, bool frozen_inputs = false,
        const std::vector<std::size_t>* cut_boundary_points = nullptr)
        : m_conflict(conflict), m_segments(segments), m_limits(limits), m_bed(bed), m_clip(clip), m_regions(regions),
          m_frozen_inputs(frozen_inputs), m_cut_boundary_points(cut_boundary_points) {}
    const LayerMaterialOverlap& m_conflict;
    const std::vector<DepositedSegment>& m_segments;
    const LayerGeometryLimits& m_limits;
    double m_bed;
    const ExPolygon* m_clip;
    const std::vector<DepositedSection::Region>& m_regions;
    // True only for deposited composition using bounded owned input copies
    // taken before its first checkpoint. Other sampler evidence stays strict.
    const bool m_frozen_inputs;
    const std::vector<std::size_t>* m_cut_boundary_points;
    const std::shared_ptr<const unsigned char> m_identity = std::make_shared<const unsigned char>(0);
};

class AuthenticatedSectionPartition {
    friend class NativeOwnershipPartitionFactory;
    friend DepositedSection section_deposition(const std::vector<DepositedSegment>&, double, double,
        const LayerGeometryLimits&, const std::function<void()>&, const SectionOwnership*, const ExPolygon*, ExactOverlayDiagnostics*);
    AuthenticatedSectionPartition() = default;
    std::shared_ptr<const unsigned char> identity;
    std::vector<DepositedSection::Region> regions;
    std::size_t work = 0;
    unsigned interface_allowance_nm = 0, representation_allowance_nm = 0;
    double dual_native_boundary_area_m2 = 0.;
    std::size_t native_winding_work = 0;
};

// A live-sampler-bound deposited-order partition, not a native or physical
// qualification. Fully hidden later materials may have no output region.
class AuthenticatedDepositedPartition {
public:
    DepositedOverlapPolicy policy() const { return m_policy; }
    double overlap_area_m2() const { return m_overlap_area_m2; }
private:
    friend class DepositedOwnershipPartitionFactory;
    friend DepositedSection section_deposition(const std::vector<DepositedSegment>&, double, double,
        const LayerGeometryLimits&, const std::function<void()>&, const SectionOwnership*, const ExPolygon*, ExactOverlayDiagnostics*);
    AuthenticatedDepositedPartition() = default;
    std::shared_ptr<const unsigned char> identity;
    std::vector<DepositedSection::Region> regions;
    std::size_t work = 0;
    DepositedOverlapPolicy m_policy = DepositedOverlapPolicy::FirstDepositionWins;
    double m_overlap_area_m2 = 0.;
};

// Proves no additional occupied-union representation change relative to the
// admitted sampled envelopes. It does not remove material-allocation
// sensitivity, bead-model error or quadrature uncertainty.
class OccupiedSectionPartition {
public:
    // Area at this sampling plane, not a cumulative area or a physical bound.
    double dual_native_boundary_ambiguity_area_m2() const { return dual_native_boundary_area_m2; }
    DepositedOverlapPolicy deposited_overlap_policy() const { return deposited_policy; }
    double deposited_overlap_area_m2() const { return deposited_overlap_area; }
    // Sum of per-material nominal area excluded by this complete partition.
    // Unlike overlap union area, a triple overlap contributes twice its area.
    double deposited_excluded_nominal_area_m2() const { return deposited_excluded_nominal_area; }
private:
    friend DepositedSection section_deposition(const std::vector<DepositedSegment>&, double, double,
        const LayerGeometryLimits&, const std::function<void()>&, const SectionOwnership*, const ExPolygon*, ExactOverlayDiagnostics*);
    friend std::size_t qualify_complete_occupied_partition(const std::vector<const DepositedSlab*>&,
        std::size_t, const std::function<void()>&);
    OccupiedSectionPartition() = default;
    std::vector<DepositedSection::Region> regions;
    double sample_z_m = 0.;
    unsigned material_allowance_nm = 0;
    double dual_native_boundary_area_m2 = 0.;
    DepositedOverlapPolicy deposited_policy = DepositedOverlapPolicy::None;
    double deposited_overlap_area = 0.;
    double deposited_excluded_nominal_area = 0.;
};

// Requires every member of one unchanged authenticated partition in one Z
// band containing its sampling plane. Throws for missing/mixed proof, members,
// altered geometry or metadata. Returned work debits the caller's budget.
std::size_t qualify_complete_occupied_partition(const std::vector<const DepositedSlab*>& slabs,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {});

struct SectionSweepLimits {
    LayerGeometryLimits geometry;
    std::size_t max_sections = 10000, max_segment_visits = 20000000;
};
// Resumable form of the same sweep, advancing one sampled plane at a time.
// Borrows immutable segments and clip; owns the plane schedule. All callbacks
// are provisional until exhausted. Failure or reentry permanently poisons it.
class DepositionSectionCursor {
public:
    DepositionSectionCursor(const std::vector<DepositedSegment>& segments,
        std::vector<double> z_m, double bed_z_m, const SectionSweepLimits& limits = {},
        std::function<void()> checkpoint = {},
        std::function<std::optional<SectionOwnership>(std::size_t,double)> ownership = {},
        const ExPolygon* clip = nullptr);
    ~DepositionSectionCursor();
    DepositionSectionCursor(const DepositionSectionCursor&) = delete;
    DepositionSectionCursor& operator=(const DepositionSectionCursor&) = delete;
    bool advance(const std::function<void(std::size_t,const DepositedSection&)>& visit);
    // Only available after complete traversal; no partial success counters.
    std::size_t completed_visits() const;
private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};
// Stream strictly increasing sections without retaining a whole-print mesh or
// rescanning every bead for every plane. The callback owns a section only for
// the duration of its call. Returns the number of intersecting bead visits.
// Geometry is identical to section_deposition, including holes and material
// overlap rejection. Exceptions invalidate the sweep; earlier callbacks are
// partial geometry, never a completed simulation or a qualification result.
std::size_t visit_deposition_sections(const std::vector<DepositedSegment> &segments,
    const std::vector<double> &z_m, double bed_z_m,
    const std::function<void(std::size_t, const DepositedSection &)> &visit,
    const SectionSweepLimits &limits = {}, const std::function<void()> &checkpoint = {},
    const std::function<std::optional<SectionOwnership>(std::size_t, double)> &ownership = {},
    const ExPolygon *clip = nullptr);

} // namespace Slic3r::Thermoelastic
