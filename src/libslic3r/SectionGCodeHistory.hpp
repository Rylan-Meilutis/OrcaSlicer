#pragma once
#include "ThermalHistory.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include "SectionGrowth.hpp"
#include "SectionIntegration.hpp"
#include "SectionThermalAssembly.hpp"
#include "SectionThermalLoad.hpp"
#include "SectionThermalSequence.hpp"
#include "SectionGrowthStream.hpp"
#include "SectionThermalSkin.hpp"
#include "ThermoelasticLocalSections.hpp"
#include "ThermoelasticLocalTopology.hpp"
#include "ThermalDeposition.hpp"
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticSections.hpp"
#include "ThermoelasticToolpath.hpp"
#include "ThermoelasticConflictOwnership.hpp"
#include "ThermoelasticExactSections.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

namespace Slic3r {
struct GCodeProcessorResult;
namespace Thermoelastic {
struct SectionGCodeScenario
{
    explicit SectionGCodeScenario(BeadSection chosen_section) : section(chosen_section) {}
    BeadSection section;
    bool acknowledge_modeled_timing = false;
    std::vector<SectionThermalMaterial> thermal_materials;
    std::vector<SectionExpansionMaterial> expansion;
    std::vector<SectionMaterial> mechanical_materials;
    std::vector<std::size_t> filament_to_material;
    Eigen::Vector3d origin_m = Eigen::Vector3d::Zero(), spacing_m = Eigen::Vector3d::Zero();
    double bed_z_m = 0., cooldown_s = 0.;
    // Opt-in bounded tile-streamed thermal history and final mechanics. Uses
    // the same regular origin/spacing above; no retained polygon snapshots.
    std::optional<Eigen::Vector2d> stream_tile_size_m;
    double stream_section_height_m = .0001;
    // Thermal surface angular quadrature; zero is exact. See stream settings.
    unsigned stream_vertical_normal_bins = 0;
    // Native ownership for streamed prefixes, queried by processor time and Z.
    // IDs use filament_to_material's target domain. Captured source data must
    // remain immutable throughout the run and return fresh sampler evidence.
    std::function<std::optional<SectionOwnership>(double, double)> stream_ownership;
    // Requires a matching live generation snapshot; never inferred from an
    // imported file's filament names. Mutually exclusive with stream_ownership.
    bool use_native_ownership = false;
    // Explicit first-deposition overlap assumption, using authenticated emitted
    // envelopes rather than native material boundaries. Not physical bonding.
    // Mutually exclusive with native and custom ownership.
    bool use_deposited_ownership = false;
    // Explicit analysis approximation: activate cumulative paths inside the
    // final fixed material regions. Reports excluded nominal volume; does not
    // simulate molten-bead displacement or conserve full commanded volume.
    bool use_fixed_interface_reference = false;
    // Strictly increasing processor-time endpoints; last equals total_s.
    // New material within each interval is activated at its endpoint, NOT
    // continuously along its path. Refine this schedule independently of dt.
    // Exposure is refreshed at these endpoints and at exposure_times_s below;
    // it is constant between those events, including during cooldown.
    std::vector<double> activation_times_s;
    // Optional boundary-law changes without new geometry/mass. Strictly
    // increasing in [0, print duration + cooldown]; coincident activations
    // already refresh exposure. Laws must include their full prior state.
    std::vector<double> exposure_times_s;
    std::function<std::vector<double>(double)> deposition_temperatures_k;
    std::function<std::vector<ThermalHistory::Element>(double, const SectionThermalAssembly&)> exposure;
    std::optional<StationaryExtrusionAllowance> stationary_allowance;
    std::optional<CollapsedEnvelopeAllowance> collapsed_envelope_allowance;
};
struct SectionGCodeLimits
{
    std::size_t max_moves = 2000000, max_segments = 250000, max_prefix_visits = 2000000;
    LayerGeometryLimits geometry;
    SectionTopologyLimits topology;
    SectionGrowthLimits growth;
    SectionThermalAssemblyLimits assembly;
    ThermalHistory::HistoryLimits history;
    SectionSolveLimits mechanics;
    LocalSectionLimits streaming_geometry;
    StreamTopologyLimits streaming_topology;
    SectionGrowthStreamLimits streaming_growth;
    StreamIntegrationLimits streaming_integration;
    SectionSkinLimits streaming_skin;
    ConflictOwnershipLimits ownership;
    std::size_t streaming_thermal_points = 8000000, streaming_thermal_work = 20000000, streaming_thermal_groups = 100000;
};
struct SectionGCodeResult
{
    ThermalHistory::HistoryResult thermal;
    std::optional<SectionMechanicalResult> mechanical;
    StationaryExtrusionSummary stationary_extrusion;
    CollapsedEnvelopeSummary collapsed_envelope;
    double modeled_print_time_s = 0., occupied_volume_m3 = 0.;
    std::size_t excluded_segments = 0, activated_snapshots = 0;
    std::size_t exposure_changes = 0;
    bool uses_modeled_timing = true, uses_batched_activation = true;
    bool uses_fixed_interface_reference = false;
    bool allows_single_sided_native_boundary = false;
    bool allows_dual_native_boundary = false;
    NativeInterfaceVariant native_interface_variant = NativeInterfaceVariant::Strict;
    DepositedOverlapPolicy deposited_overlap_policy = DepositedOverlapPolicy::None;
    unsigned vertical_normal_bins = 0;
    double section_partition_increment_m = 0.;
    double excluded_nominal_volume_m3 = 0.;
    // Volume conduction quadrature; boundary operators come from the explicit
    // exposure callback and must be assessed separately.
    SectionIntegration::ThermalQuadrature thermal_quadrature = SectionIntegration::ThermalQuadrature::Consistent;
};
// Read-only bounded prefix adapter. The retained reference uses planar
// rectangular beads; optional streaming also supports rounded/scarf sampling.
// Material overlap without resolved ownership, missing source data,
// non-monotone geometry and exhausted budgets reject; no gap filling or nozzle
// target substitution. No result here qualifies discretization or adhesion.
SectionGCodeResult simulate_section_gcode_history(const GCodeProcessorResult& source,
                                                  const SectionGCodeScenario& scenario,
                                                  const SectionGCodeLimits& limits        = {},
                                                  const std::function<void()>& checkpoint = {});
} // namespace Thermoelastic
} // namespace Slic3r
