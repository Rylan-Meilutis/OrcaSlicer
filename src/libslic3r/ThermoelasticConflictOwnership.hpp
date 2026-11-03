#pragma once
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticToolpath.hpp"
#include "ExPolygon.hpp"
#include <cstddef>
#include <array>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace Slic3r {
struct GCodeProcessorResult;
namespace Thermoelastic {
enum class NativeInterfaceVariant { Strict, FavorFirst, FavorSecond };
struct NativeOwnershipPartitionFailure
{
    const SampledSectionConflict& sampled;
    const std::array<ExPolygons, 2>& native_masks;
    const std::array<ExPolygons, 2>& envelopes;
    NativeInterfaceVariant variant;
    std::size_t max_work;
    bool allow_single_sided_native_boundary;
    bool allow_dual_native_boundary;
    const std::string& message;
};
struct ConflictOwnershipLimits
{
    std::size_t max_moves = 2000000, max_segments = 2000000;
    std::size_t max_polygon_checks = 4096, max_polygon_work = 2000000;
    std::size_t max_source_work = 1000000;
    std::size_t max_total_work = 80000000; // Prepared resolver: all calls combined.
    LayerGeometryLimits geometry;
    // Explicit numerical sensitivity scenario, never a physical material priority.
    NativeInterfaceVariant interface_variant = NativeInterfaceVariant::Strict;
    // Exact-only chronological assumption, independent of native boundaries.
    DepositedOverlapPolicy deposited_overlap_policy = DepositedOverlapPolicy::None;
    // Optional0..7nm source-boundary sensitivity, only with a non-strict
    // interface scenario. Neither source masks nor deposited claims expand.
    unsigned source_centerline_boundary_allowance_nm = 0;
    // Explicit analysis-only representation assumption. An uncovered deposited
    // conflict must fit one native 7nm boundary band and be strictly separated
    // from every other native boundary band. Strict ownership never enables it.
    bool allow_single_sided_native_boundary = false;
    // Exact-only conditional material-location sensitivity. Already deposited
    // conflict outside both native references must fit BOTH7nm bands. Both
    // priority scenarios need separate assessment; neither is physical truth.
    bool allow_dual_native_boundary = false;
    // Optional read-only investigation hook. References are borrowed only for
    // the synchronous call; true denotes conflict-clipped output, false native
    // snapshot input. No callback is installed by production scenario setup.
    std::function<void(double, bool, const std::vector<DepositedSection::Region>&)> observe_masks;
    // Failure-only borrowed operands, at most2048 vertices each. Diagnostic
    // observer exceptions are ignored so they cannot replace the rejection.
    std::function<void(const NativeOwnershipPartitionFailure&)> observe_partition_failure;
};

// Full modeled timeline only. Reconstructs nominal paths from the matching result,
// never accepts caller-invented bead/source associations. All geometric contributors
// to the supplied complete conflict are checked, including brim/support/tower.
// Unrelated unowned deposition is not a veto. Unresolved spatial deposition remains
// unsupported. Returned masks are clipped to this conflict; they cannot authorize
// other contacts. No geometry/mass is changed and no physical interface is proven.
SectionOwnership native_ownership_for_overlap(const GCodeProcessorResult& result,
                                              const LayerMaterialOverlap& conflict,
                                              BeadSection section,
                                              double bed_z_m,
                                              const std::vector<std::size_t>& filament_to_material,
                                              const ConflictOwnershipLimits& limits   = {},
                                              const std::function<void()>& checkpoint = {});

// Prepare once for exactly the selected model-deposition prefix. Uses the same
// role exclusions and explicit stationary approximation as section studies.
// The returned single-threaded callback borrows result, which must remain immutable and alive.
// It verifies conflicts against ORIGINAL authenticated geometry. The explicit
// approximation option retains uniquely occupied original material where a
// foreign envelope has expanded into it; original void remains unsupported.
// The separate FirstDepositionWins policy requires the sampled callback and
// exact complete move children. The existing explicitly bounded stationary
// omission is allowed; its physical E is not part of the retained envelope
// union and its summary must be reported separately. Collapsed-envelope
// omission remains unsupported. It authenticates live command
// geometry/order without native boundary admission and returns only the
// deposited_partition member; wire it to resolve_deposited_partition.
using ConflictOwnershipResolver = std::function<SectionOwnership(const LayerMaterialOverlap&)>;
using SampledConflictOwnershipResolver = std::function<SectionOwnership(const SampledSectionConflict&)>;
ConflictOwnershipResolver prepare_native_model_ownership(const GCodeProcessorResult& result,
    BeadSection section, double bed_z_m, const std::vector<std::size_t>& filament_to_material,
    double observation_time_s, const std::optional<StationaryExtrusionAllowance>& stationary_allowance = {},
    const ConflictOwnershipLimits& limits = {}, const std::function<void()>& checkpoint = {},
    bool retain_original_unique_material = false, SampledConflictOwnershipResolver* sampled_resolver = nullptr,
    const std::optional<CollapsedEnvelopeAllowance>& collapsed_allowance = {});
} // namespace Thermoelastic
} // namespace Slic3r
