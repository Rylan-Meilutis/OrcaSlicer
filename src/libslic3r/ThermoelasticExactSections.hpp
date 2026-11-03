#pragma once
#include "ExPolygon.hpp"
#include "Point.hpp"
#include "BoundingBox.hpp"
#include "SectionIntegration.hpp"
#include <cstddef>
#include <array>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
enum class ExactSectionOperation { Intersection, Difference, SymmetricDifference, Union };
// Immutable canonical rational geometry. Copies share bounded storage; no
// integer polygon conversion is provided. Default construction is empty.
// Consumers must bound the total retained vertex counts across live handles.
class ExactSection {
public:
    ExactSection() = default;
    bool empty() const;
    std::size_t vertex_count() const;
    // Outward integer enclosure for indexing only, never exact geometry.
    std::optional<BoundingBox> bounds() const;
private:
    friend class ExactSectionAccess;
    struct Impl;
    explicit ExactSection(std::shared_ptr<const Impl> impl) : m_impl(std::move(impl)) {}
    std::shared_ptr<const Impl> m_impl;
};
struct ExactSectionResult {
    ExactSection section;
    std::size_t work = 0;
};
// Optional caller-owned diagnostic sink, never a budget or admission input.
// Categories count successful existing preflight debits only; fast paths and
// other Boolean work are excluded. Failure retains the charged prefix.
// Saturation marks overflow rather than changing the operation's outcome.
struct ExactOverlayDiagnostics {
    std::size_t overlays=0, edges=0, candidates=0;
    std::size_t boxes_projection=0, sorting=0, sweep=0, active_scans=0, event_reservation=0;
    bool overflow=false;
};
struct ExactSectionComponents {
    std::vector<ExactSection> sections;
    std::size_t work = 0;
};
// Canonical exact nonzero union of native boundaries, without integer rounding.
// At most1024 input vertices,4096 arrangement vertices and4M work.
ExactSectionResult make_exact_section(const ExPolygons& polygons, std::size_t max_work = 4000000,
    const std::function<void()>& checkpoint = {});
// Boolean composition of immutable exact sets. Combined input edges and
// output vertices are capped at4096; work is capped at4M. Algebraic identities
// may return shared immutable storage without an overlay.
ExactSectionResult boolean_exact_sections(const ExactSection& left, const ExactSection& right,
    ExactSectionOperation operation, std::size_t max_work = 4000000,
    const std::function<void()>& checkpoint = {}, ExactOverlayDiagnostics* diagnostics = nullptr);
enum class DepositedOverlapPolicy { None, FirstDepositionWins };
struct ExactDepositedCommandSection {
    std::size_t material = 0;
    ExactSection envelope;
};
struct ExactDepositedPartition {
    DepositedOverlapPolicy policy = DepositedOverlapPolicy::FirstDepositionWins;
    // One canonical union per material, with no same-material command barriers.
    std::vector<ExactDepositedCommandSection> regions;
    ExactSection occupied_union, cross_material_overlap;
    double cross_material_overlap_area_m2 = 0.;
    std::size_t work = 0;
};
// Untrusted geometry primitive: vector order is the caller's chronological
// assumption, NOT authenticated command provenance or physical displacement.
// Each command claims only its envelope outside all earlier envelopes. Every
// prefix preserves the exact occupied union and disjoint material interiors.
// Overlap is the spatial union of competing-material envelopes, counted once.
// No native-boundary or mechanical qualification token is created. Callers
// must independently authenticate order before using this in a history model.
// At most4096 commands/input vertices,256 material IDs (0..255),4096 combined
// output vertices, and4M cumulative work; standard exact overlay limits apply.
ExactDepositedPartition exact_first_deposition_partition(
    const std::vector<ExactDepositedCommandSection>& commands, std::size_t max_work = 4000000,
    const std::function<void()>& checkpoint = {});
struct ExactDepositedCommandResult {
    ExactDepositedCommandSection command;
    std::size_t work = 0;
};
// Streaming variant of the same untrusted ordered-geometry primitive. The
// provider is called once per index, in order, with the remaining shared work
// budget. It must bound construction to that budget and report consumed work;
// provider work, input visits and partition proofs together cannot exceed4M.
// At most4096 commands,1024 canonical vertices per incoming cut, and4096
// stored state plus pending same-material run union and incoming vertices at
// stable streaming boundaries (aliases conservatively counted). While a run
// is certified, old operands and replacement geometry use bounded operation
// scratch; this is not a simultaneous-allocation limit for that transition.
// Only contiguous nonempty runs are combined; each
// requested command prefix ends and certifies its own final run.
// Existing exact-operation scratch limits also apply. No cumulative input
// geometry is retained. Exceptions expose no partial partition.
ExactDepositedPartition exact_first_deposition_partition_streamed(std::size_t command_count,
    const std::function<ExactDepositedCommandResult(std::size_t, std::size_t)>& provider,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {});
ExactSectionComponents exact_section_components(const ExactSection& section,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {});
// Components of the occupied open interior, not connected closed boundaries.
// Selects canonical arrangement faces: point-touching solids stay separate,
// while tangent holes retain their exact directed boundaries within one face.
// No cuts, rounding, or contact/mechanical qualification are introduced.
// At most4096 combined retained vertices and4M work; exceptions publish nothing.
// This does not change the stricter policy of validate_exact_section_topology.
ExactSectionComponents exact_section_interior_components(const ExactSection& section,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {});
// Reject repeated rational boundary vertices within/across rings of each
// canonical occupied component. A canonical point-touch pinch is rejected,
// never arbitrarily split or treated as a positive-length/area load path.
// Returns consumed work, never changes the canonical geometry.
std::size_t validate_exact_section_topology(const ExactSection& section,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {});
struct ExactSectionPointMembership {
    bool contained = false, boundary = false;
    std::size_t work = 0;
};
// CLOSED-set membership of an exact native-nm query point. Outer and hole
// boundaries are contained, hole interiors are not. Strict interior requires
// contained && !boundary. This predicate does not certify component identity.
ExactSectionPointMembership exact_closed_contains_point(const ExactSection& section, const Point& point,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {});
struct ExactSectionBilinearExtrema {
    double minimum = 0., maximum = 0.;
    std::array<double,2> minimum_bounds{}, maximum_bounds{};
    std::array<double,2> minimum_local{}, maximum_local{};
    std::size_t work = 0;
};
// Extrema on a nonempty closed section wholly inside the native-nm cell.
// Corner order is(0,0),(1,0),(0,1),(1,1); finite double values are interpreted
// exactly. Rational contour/hole vertices and edge stationary points determine
// the extrema: a bilinear field has no strict interior maximum or minimum.
// Value bounds are outward certified; local XY coordinates are representative
// doubles only, never geometry witnesses. Empty/outside sections are rejected.
// At most4096 boundary vertices and4M work; caller accounts extraction scratch.
ExactSectionBilinearExtrema exact_section_bilinear_extrema(const ExactSection& section,
    const std::array<coord_t,2>& cell_start, const std::array<coord_t,2>& cell_size,
    const std::array<double,4>& corner_values, std::size_t max_work = 4000000,
    const std::function<void()>& checkpoint = {});
// Directed rational boundary segments, with material on the left. No public
// endpoint constructor/conversion: only canonical sections and exact boundary
// operations can create these handles. Copies share immutable storage.
class ExactSectionBoundary {
public:
    ExactSectionBoundary() = default;
    bool empty() const;
    std::size_t segment_count() const;
    std::optional<BoundingBox> bounds() const;
private:
    friend class ExactSectionBoundaryAccess;
    struct Impl;
    explicit ExactSectionBoundary(std::shared_ptr<const Impl> impl) : m_impl(std::move(impl)) {}
    std::shared_ptr<const Impl> m_impl;
};
struct ExactSectionBoundaryResult {
    ExactSectionBoundary boundary;
    std::size_t work = 0;
};
enum class ExactBoundaryOperation { Intersection, Difference };
enum class ExactBoundaryOrientation { Any, Opposed };
// Positive-length collinear spans only: isolated crossings/touching points are
// never contact. Results retain the left operand's direction. Opposed admits
// only oppositely directed right spans; Any ignores the right orientation.
// Each operation caps output at4096 segments and work at4M. Callers also debit
// total retained segment counts across live handles (two rational endpoints each).
ExactSectionBoundaryResult extract_exact_section_boundary(const ExactSection& section,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {});
ExactSectionBoundaryResult boolean_exact_section_boundaries(const ExactSectionBoundary& left,
    const ExactSectionBoundary& right, ExactBoundaryOperation operation, ExactBoundaryOrientation orientation,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {});
// Select existing boundary segments exactly on native X(axis0) or Y(axis1)
// plane. This is a tile-FACE trace only after exact tile clipping; it does not
// compute an interior plane section or snap crossing segments onto the plane.
ExactSectionBoundaryResult exact_boundary_on_native_plane(const ExactSectionBoundary& boundary,
    unsigned axis, coord_t coordinate, std::size_t max_work = 4000000,
    const std::function<void()>& checkpoint = {});
struct ExactBoundarySegment {
    ExactSectionBoundary span;
    std::array<double,2> outward_normal{};
    std::array<std::array<double,2>,2> outward_normal_bounds{};
};
struct ExactBoundarySegments {
    std::vector<ExactBoundarySegment> segments;
    std::size_t work = 0;
};
// One immutable directed span per record, with finite unit XY outward normal
// (right of the material-left direction). Irrational normalization is bounded
// by exact comparisons; native endpoint rounding is never used. Splitting
// does not certify exposure or contacts. At most4096 spans and4M work; callers
// must account all retained spans across live results.
ExactBoundarySegments split_exact_boundary_segments(const ExactSectionBoundary& boundary,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {});
struct ExactSectionMoments {
    SectionIntegration::Moments values{};
    // Outward bounds after exact polynomial integration and SI conversion.
    std::array<std::array<double, 2>, 27> bounds{};
    std::size_t work = 0;
};
struct ExactSectionMeasure {
    double value = 0.;
    std::array<double, 2> bounds{};
    std::size_t work = 0;
};
// Scalar equivalents of the zeroth volume/horizontal-area moment below.
// Same exact cell clipping, native-nm domain, SI units and closed horizontal
// Z interval. Directed rational boundaries (including tangent holes) are
// integrated before one representative/outward-interval conversion. An exact
// enclosing-bounds proof avoids redundant XY intersection. No contact or
// component qualification; standard4096 retained vertices/4M work limits apply.
ExactSectionMeasure exact_section_volume(const ExactSection& section,
    const std::array<coord_t, 3>& cell_start, const std::array<coord_t, 3>& cell_size,
    coord_t bottom, coord_t top, std::size_t max_work = 4000000,
    const std::function<void()>& checkpoint = {});
ExactSectionMeasure exact_section_horizontal_area(const ExactSection& section,
    const std::array<coord_t, 3>& cell_start, const std::array<coord_t, 3>& cell_size,
    coord_t z, std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {});
// Native nonzero unions, then the requested Boolean and exact cell clipping.
// Difference with an empty right operand integrates the left union. Cell and
// Z coordinates are native nm; moments use normalized cell coordinates and
// SI volume. Rational vertices are never rounded before integration. This is
// numerical integration only, not component/contact or source qualification.
// At most1024 combined input vertices,4096 arrangement vertices and4M work.
ExactSectionMoments exact_section_volume_moments(const ExPolygons& left, const ExPolygons& right,
    ExactSectionOperation operation, const std::array<coord_t, 3>& cell_start,
    const std::array<coord_t, 3>& cell_size, coord_t bottom, coord_t top,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {});
ExactSectionMoments exact_section_volume_moments(const ExactSection& section,
    const std::array<coord_t, 3>& cell_start, const std::array<coord_t, 3>& cell_size,
    coord_t bottom, coord_t top, std::size_t max_work = 4000000,
    const std::function<void()>& checkpoint = {});
// Exact cell-clipped horizontal face at fixed native z. The27 moments carry
// SI AREA and normalized cell powers, not volume. A plane outside the closed
// cell Z interval contributes zero; either cell boundary plane is included.
// Caller must subtract all covered/contact regions before supplying the face.
// Assign each physical face to one cell at a shared Z plane to avoid counting
// it twice. No surface normal or boundary-condition orientation is inferred.
ExactSectionMoments exact_section_horizontal_surface_moments(const ExactSection& section,
    const std::array<coord_t, 3>& cell_start, const std::array<coord_t, 3>& cell_size,
    coord_t z, std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {});
enum class ExactBoundaryCellOwnership { MaterialSide, AllIncident };
// SI-area moments of the supplied directed boundary extruded over[bottom,top],
// exactly clipped to the cell. MaterialSide assigns a cell-aligned edge only
// to its material-left cell; AllIncident includes either incident cell and
// requires external ownership to avoid double counting. The caller must first
// subtract internal/contact spans; this primitive does not establish exposure
// or infer a single normal for differently oriented segments. Irrational edge
// lengths receive exact-comparison-certified outward floating-point bounds.
ExactSectionMoments exact_section_vertical_surface_moments(const ExactSectionBoundary& boundary,
    const std::array<coord_t, 3>& cell_start, const std::array<coord_t, 3>& cell_size,
    coord_t bottom, coord_t top, ExactBoundaryCellOwnership ownership,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {});
struct ExactSectionPredicate {
    bool empty = false;
    std::size_t work = 0;
};
// Positive-length overlap with a native interval on an axis-aligned plane.
// axis selects the fixed coordinate; [lo,hi] uses the other coordinate.
// Orientation is ignored, isolated endpoint contacts are empty, and no
// uncertainty margin or endpoint rounding is applied. Requires lo <= hi.
ExactSectionPredicate exact_boundary_overlaps_native_interval(const ExactSectionBoundary& boundary,
    unsigned axis, coord_t plane, coord_t lo, coord_t hi,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {});
struct ExactSegmentCoverage {
    bool covered = false;
    // Separate numerical sensitivity result; never changes strict coverage.
    bool within_boundary_band = false;
    std::size_t work = 0;
};
// Closed-set coverage of the entire native-coordinate segment. Exact boundary
// crossings partition the segment; no rounded clipping fragments or tolerance.
// Optional0..7nm band reports a SEPARATE analysis sensitivity, not coverage or
// material assignment. It does not offset polygons or repair source geometry.
ExactSegmentCoverage exact_section_covers_segment(const ExPolygons& polygons, const Point& start, const Point& end,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {}, unsigned boundary_allowance_nm = 0);
// Native slicer clipping uses nonzero winding across the original oriented
// boundary chains. Crossing vertices remain exact rational
// intersections; no rounding, offset or modification of the input masks.
// This establishes coverage under that fill rule, not simple-ring validity.
ExactSegmentCoverage exact_native_section_covers_segment(const ExPolygons& polygons, const Point& start, const Point& end,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {}, unsigned boundary_allowance_nm = 0);
struct ExactSectionTopologySummary {
    std::size_t components = 0, holes = 0, work = 0;
};
struct ExactIntegerSectionUnion {
    ExPolygons polygons;
    std::size_t work = 0;
};
class NonIntegralSectionUnion : public std::invalid_argument {
public:
    explicit NonIntegralSectionUnion(std::size_t consumed)
        : std::invalid_argument("Exact occupied union requires a noninteger boundary vertex."), work(consumed) {}
    std::size_t work;
};
struct ExactOccupiedAssignment {
    std::vector<std::size_t> occupied_indices;
    unsigned representation_allowance_nm = 1;
    std::size_t work = 0;
};
struct ExactBoundaryDifference {
    double added_area_bound_m2 = 0., removed_area_bound_m2 = 0.;
    std::size_t work = 0;
};
// Whole exact native symmetric difference inside the original boundary's
// Euclidean capsules, with exact component/hole correspondence. Returns
// outward-rounded area bounds without rounding the difference geometry.
// This does not replace downstream material and morphology qualifications.
// At most1024 combined input vertices,4096 arrangement vertices and4M work.
ExactBoundaryDifference certify_exact_boundary_difference(const ExPolygons& original,
    const ExPolygons& represented, unsigned allowance_nm, std::size_t max_work = 4000000,
    const std::function<void()>& checkpoint = {});
struct ExactHoleChangeDiagnostic {
    bool source_hole_unique = false, represented_hole_unique = false, correspondence_unique = false;
    std::optional<bool> hole_unchanged, difference_empty;
    // Outward interval in squared native nanometres. Missing means unknown
    // correspondence or an empty difference, not zero separation.
    std::optional<std::array<double,2>> change_distance_squared_nm2;
    std::size_t work = 0;
};
// Diagnostic only: locate a hole by an oriented native hole-ring hint, then
// compare its exact source/represented regions and their separation from the
// complete exact occupied-set change. The filled hole-ring region includes its
// boundary and any separately occupied islands nested inside the hole.
// No material-uncertainty localization or admission follows from this result.
ExactHoleChangeDiagnostic diagnose_exact_hole_change(const ExPolygons& original,
    const ExPolygons& represented, const Polygon& hole_hint, std::size_t max_work = 4000000,
    const std::function<void()>& checkpoint = {});
class OccupiedAssignmentOverlap : public std::invalid_argument {
public:
    OccupiedAssignmentOverlap(const std::string& message, bool source, std::size_t a, std::size_t b, std::size_t consumed)
        : std::invalid_argument(message), original(source), first(a), second(b), work(consumed) {}
    bool original;
    std::size_t first, second, work;
};
// Diagnostic certificate for an existing occupied representation. Original
// pieces remain unchanged. Requires disjoint pieces, exact component/hole
// correspondence, full component containment and a fixed1nm difference band.
// Differences touching a tile face reject, preserving exact face traces.
ExactOccupiedAssignment certify_occupied_assignment(const ExPolygons& pieces,
    const ExPolygons& occupied, const BoundingBox& tile, std::size_t max_work = 4000000,
    const std::function<void()>& checkpoint = {});
// Exact native nonzero union, with no rounded output. Rejects any rational
// boundary vertex not exactly on the native integer grid. Bounded to1024 input
// and4096 output vertices under the independent4M work ceiling.
ExactIntegerSectionUnion exact_integer_section_union(const ExPolygons& polygons,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {});
// Diagnostic only: counts alone do not establish component correspondence.
// Uses the exact union, without rounding Boolean output to the native grid.
ExactSectionTopologySummary diagnose_exact_section_topology(const ExPolygons& polygons,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {});
// Exact one-to-one correspondence of components AND each component's holes.
// Throws on changed or ambiguous topology; returns charged work on success.
std::size_t qualify_exact_section_topology(const ExPolygons& before, const ExPolygons& after,
    std::size_t max_work, const std::function<void()>& checkpoint = {});
// Local exact predicate on ORIGINAL native-coordinate polygons. No rounded
// diagnostic difference or rational output mesh is substituted for the input.
// Default1024 total vertices; explicit dense queries may request2048 while
// retaining the independent4M quadratic-work ceiling. Work is reserved before
// the noninterruptible CGAL operation, with cancellation before and after.
// Explicit native fill interprets already validated slicer boundary chains by
// nonzero winding, including touching holes; it retains a hard1024 input limit.
// The default still requires valid simple polygons with properly placed holes.
ExactSectionPredicate exact_section_empty(const ExPolygons& a, const ExPolygons& b,
    ExactSectionOperation operation, std::size_t max_vertices = 1024,
    std::size_t max_work = 2000000, const std::function<void()>& checkpoint = {}, bool native_nonzero_fill = false);
struct ExactOwnershipWitness {
    std::optional<bool> in_first_native, in_second_native;
    std::size_t work = 0;
};
enum class NativeWitnessFill { ValidatedSimple, NonzeroWinding };
struct ExactReferenceClipping {
    std::array<ExPolygons, 2> sections;
    unsigned representation_allowance_nm = 0;
    double vertex_rounding_bound_m = 0.;
    std::size_t work = 0;
};
// Intersect each deposited region with its fixed reference using exact
// arithmetic, then jointly round their shared boundary vertices. No material
// priority or activation of whole reference regions is used. The exact result
// is a subset of deposited material; the returned integer representation has
// the reported boundary rounding allowance. Components/holes and disjointness
// are certified after rounding. Rejects overlapping reference materials.
// charge_work is called BEFORE each bounded operation, for a caller's shared
// cumulative budget; callback exceptions propagate unchanged.
// Each operand is capped at2048 vertices; all constructions, overlays and
// certification still share the independent4M work ceiling.
ExactReferenceClipping exact_clip_reference_sections(const std::array<ExPolygons, 2>& reference,
    const std::array<ExPolygons, 2>& deposited, std::size_t max_work = 4000000,
    const std::function<void()>& checkpoint = {},
    const std::function<void(std::size_t)>& charge_work = {});
struct UnqualifiedInterfaceVariants {
    // Full claim-pair outputs, not patches to subtract/rejoin independently.
    std::array<ExPolygons, 2> favor_first, favor_second;
    unsigned interface_allowance_nm = 0;
    unsigned representation_allowance_nm = 0;
    double vertex_rounding_bound_m = 0.;
    double single_sided_native_boundary_area_m2 = 0.;
    std::size_t work = 0;
};
enum class ExactInterfaceChoice { FavorFirst, FavorSecond };
struct ExactInterfacePartition {
    std::array<ExactSection, 2> materials;
    unsigned interface_allowance_nm = 0, representation_allowance_nm = 0;
    double single_sided_native_boundary_area_m2 = 0.;
    // Already deposited conflict outside both native references, certified in
    // both7nm bands and assigned only by an explicit material sensitivity.
    double dual_native_boundary_area_m2 = 0.;
    // Subset of work spent constructing the two native nonzero-winding sets.
    std::size_t native_winding_work = 0;
    std::size_t work = 0;
};
// Same native ownership/coverage proof as the analysis variants below, but
// retains the selected full material sets before integer vertex rounding.
// Exact disjointness and occupied-union conservation are proved against the
// complete envelopes. Zero representation allowance does NOT remove inherited
// interface sensitivity. This is not downstream contact/solver qualification.
// Both retained arrangements together are capped at4096 vertices; all stages
// share the4M work ceiling. No source geometry or printer toolpath is changed.
// Optional dual-native-boundary sensitivity assigns already deposited conflict
// outside both references only after whole-region certification in BOTH7nm
// bands. FavorFirst and FavorSecond must be assessed separately; neither is a
// physical ownership certificate or an exhaustive uncertainty bound. Default
// false; this option is deliberately unavailable to rounded legacy variants.
ExactInterfacePartition exact_partition_uncertain_interface(const ExPolygons& first_native,
    const ExPolygons& second_native, const std::array<ExPolygons, 2>& complete_envelopes,
    ExactInterfaceChoice choice, unsigned allowance_nm = 7, std::size_t max_work = 4000000,
    const std::function<void()>& checkpoint = {}, NativeWitnessFill native_fill = NativeWitnessFill::ValidatedSimple,
    bool allow_single_sided_native_boundary = false, bool allow_dual_native_boundary = false);
// Explicit ANALYSIS ONLY prototype. Opposite assignments of ambiguous native
// interface material, never a physical priority or an exhaustive force bound.
// Requires complete native coverage and exact dual-side interface-band proof.
// Caller MUST still validate full occupied-union representation, individual
// components/holes and contacts, then independently converge BOTH studies.
// With complete_envelopes, claims must be empty: the exact arrangement derives
// claims from native masks and the complete deposited intersection, and returns
// full deposited regions in one rounding pass (not independently clipped masks).
// Input validation is capped at2048 vertices per separately constructed
// operand, not their sum; validation and Boolean work share one cumulative
// budget. Each later arrangement operation reserves its own work before CGAL.
// Explicit NonzeroWinding interprets native masks only; deposited claims and
// envelopes still require simple oriented boundaries. No input is rewritten.
// Explicit single-sided mode additionally allows uncovered DEPOSITED overlap
// entirely within one native boundary band and provably outside the other's.
// It reports the reassigned area, never fills deposited void, and leaves
// ambiguous two-sided gaps unsupported. Default requires exact native coverage.
UnqualifiedInterfaceVariants partition_uncertain_interface(const ExPolygons& first_claim, const ExPolygons& second_claim,
    const ExPolygons& first_native, const ExPolygons& second_native, unsigned allowance_nm = 7,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {},
    const std::array<ExPolygons, 2>* complete_envelopes = nullptr,
    NativeWitnessFill native_fill = NativeWitnessFill::ValidatedSimple,
    bool allow_single_sided_native_boundary = false);
struct ExactIntersectionOwnership {
    bool nonempty = false, first_unique = false, second_unique = false;
    bool native_partition = false;
    bool native_covers = false, native_disjoint = false;
    std::size_t work = 0;
};
// Whole regularized claim intersection, constructed without rounded output.
// A unique owner covers all of it and the competing native set covers none.
ExactIntersectionOwnership exact_intersection_ownership(const ExPolygons& first_claim, const ExPolygons& second_claim,
    const ExPolygons& first_native, const ExPolygons& second_native,
    std::size_t max_work = 4000000, const std::function<void()>& checkpoint = {});
// Diagnostic only: a strictly interior exact rational witness of the small
// claim intersection is classified against independently validated native
// polygons. Unknown is returned if bounded centroid search finds no witness.
// This does not certify ownership or permit a material partition repair.
// Explicit NonzeroWinding mode diagnoses native fill semantics even for weakly
// simple rings; it does NOT validate their topology or authorize those inputs.
ExactOwnershipWitness diagnose_ownership_witness(const ExPolygons& first_claim, const ExPolygons& second_claim,
    const ExPolygons& first_native, const ExPolygons& second_native,
    std::size_t max_work = 20000000, const std::function<void()>& checkpoint = {},
    NativeWitnessFill native_fill = NativeWitnessFill::ValidatedSimple);
} // namespace Slic3r::Thermoelastic
