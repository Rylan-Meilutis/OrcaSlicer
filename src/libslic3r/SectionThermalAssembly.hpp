#pragma once

#include "ThermalDeposition.hpp"
#include "SectionIntegration.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include "SectionGrowth.hpp"
#include "ThermoelasticLocalSections.hpp"
#include "ThermoelasticLocalTopology.hpp"
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticMesh.hpp"
#include "SectionThermalSkin.hpp"
#include "SectionGrowthStream.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <array>
#include <functional>
#include <optional>
#include <vector>

namespace Slic3r::Thermoelastic {

struct SectionThermalMaterial {
    Eigen::Matrix3d conductivity_w_m_k = Eigen::Matrix3d::Zero();
    double volumetric_capacity_j_m3_k = 0.;
    double initial_temperature_k = 0.;
};
struct SectionThermalAssemblyLimits {
    std::size_t max_cells = 20000, max_nodes = 100000, max_groups = 100000;
    std::size_t max_points = 2000000, max_work = 20000000;
    SectionIntegration::ThermalQuadrature quadrature = SectionIntegration::ThermalQuadrature::Consistent;
};
struct SectionThermalNode {
    Eigen::Vector3d position_m;
    std::size_t component;
};
struct SectionThermalGroup {
    std::array<std::size_t,8> nodes;
    std::size_t material, component;
    std::array<int,3> cell;
    Eigen::Vector3d minimum_m, maximum_m;
    double occupied_volume_m3;
};
enum class SectionHorizontalFaceKind { Bed, Upward, Downward };
struct SectionHorizontalThermalFace {
    std::array<std::size_t,8> nodes;
    std::size_t material;
    SectionHorizontalFaceKind kind;
    SectionIntegration::Moments area_moments;
};
struct SectionVerticalThermalFace {
    std::array<std::size_t,8> nodes;
    std::size_t material;
    Eigen::Vector3d outward_normal;
    SectionIntegration::Moments area_moments;
};
struct SectionThermalAssembly {
    // All supplied geometry is initialized at this event's time (zero).
    // This is a thermal initial condition, NOT inferred print activation.
    ThermalHistory::DepositionEvent initial;
    std::vector<SectionThermalNode> nodes;
    // Same order as initial.added_elements, for material-specific field
    // recovery on the identical qualified geometry and background grid.
    std::vector<SectionThermalGroup> groups;
    // External skin of the supplied piecewise-prismatic geometry. Internal
    // layer/material contacts are excluded; holes and exposed undersides remain.
    std::vector<SectionHorizontalThermalFace> horizontal_faces;
    std::vector<SectionVerticalThermalFace> vertical_faces;
    bool exposure_complete = false;
    unsigned vertical_normal_bins = 0;
    double occupied_volume_m3 = 0.;
    std::size_t occupied_cells = 0;
    SectionIntegration::ThermalQuadrature quadrature = SectionIntegration::ThermalQuadrature::Consistent;
};

// Geometry-certified fixed-domain assembly on a regular native-nanometre Q1
// grid. Integrates only occupied polygons, retains constitutive interfaces and
// separates physically disconnected components even in the same grid cell.
// Shared nodes imply ideal thermal contact within a certified component.
// The caller must independently refine the cut-cell basis and geometry.
// No convection is invented: all surfaces are insulated until explicit
// boundary operators are supplied. Not a G-code deposition-time adapter.
SectionThermalAssembly assemble_section_thermal_volume(const SectionTopology& topology,
    const std::vector<SectionThermalMaterial>& materials, const Eigen::Vector3d& origin_m,
    const Eigen::Vector3d& spacing_m, const SectionThermalAssemblyLimits& limits = {},
    const std::function<void()>& checkpoint = {});

// Integrates ONLY certified new material, using the next snapshot's component
// identities. No exposed skin is returned: skin comes from the complete next
// geometry, including coverage of previously exposed faces. Its node indices
// are local and must be mapped into the active history before adding mass.
SectionThermalAssembly assemble_section_thermal_addition(const SectionGrowth& growth,
    const std::vector<SectionThermalMaterial>& materials, const Eigen::Vector3d& origin_m,
    const Eigen::Vector3d& spacing_m, const SectionThermalAssemblyLimits& limits = {},
    const std::function<void()>& checkpoint = {});

struct ThermalSectionStreamSettings {
    Eigen::Vector3d origin_m = Eigen::Vector3d::Zero(), spacing_m = Eigen::Vector3d::Zero();
    Eigen::Vector2d tile_size_m = Eigen::Vector2d::Zero();
    double bed_z_m = 0., maximum_height_m = 0.;
    unsigned interval_subdivisions = 1;
    LocalSectionLimits geometry;
    SectionTopologyLimits topology;
    StreamTopologyLimits streaming_topology;
    SectionThermalAssemblyLimits assembly;
    bool include_exposure = false;
    // Zero preserves exact normals. Otherwise round side normals to this
    // many uniform angular directions (multiple of four, 8..4096). Area
    // moments remain exact; refine independently for directional heat transfer.
    unsigned vertical_normal_bins = 0;
    SectionSkinLimits skin;
    // Explicit streaming ceilings, separate from fixed-domain retained limits.
    // Offline callers may request up to2B cumulative CPU work; live storage
    // remains bounded independently by the assembly/group limits.
    std::size_t max_processed_points = 8000000, max_work = 20000000;
    std::size_t max_provisional_groups = 100000;
};
struct ThermalSectionStreamResult {
    LocalSectionMetrics geometry;
    StreamTopologyStatistics topology;
    SectionThermalAssembly assembly;
};
// Closed sampler -> connectivity -> thermal volume assembly. Borrows each
// packet and releases its polygons; provisional components merge only after
// complete traversal and source-matched certification. No partial result.
// include_exposure also streams complete horizontal/vertical skin, removing
// tile seams and material contacts before integrating its area moments.
// Otherwise exposure_complete is false: empty face lists are not permission
// to silently assume whole-part insulation. Boundary coefficients remain
// explicit inputs; this does not infer deposition timing or fan heat transfer.
ThermalSectionStreamResult assemble_local_deposition_thermal_volume(const std::vector<DepositedSegment>& segments,
    const std::vector<SectionThermalMaterial>& materials, const ThermalSectionStreamSettings& settings,
    const std::function<std::optional<SectionOwnership>(double)>& ownership = {},
    const std::function<void()>& checkpoint = {}, const std::vector<DepositedSegment>* partition_envelopes = nullptr);

struct SectionThermalGrowthAssembly {
    SectionThermalAssembly next, added;
    std::vector<std::size_t> previous_components;
    double previous_excluded_nominal_volume_m3 = 0., next_excluded_nominal_volume_m3 = 0.;
};
struct ReferenceThermalAssembly {
    SectionThermalAssembly assembly;
    double excluded_nominal_volume_m3 = 0.;
};
// Explicit fixed-interface reference geometry. Complete skin is always
// assembled. Excluded nominal volume is reported, never deposited as heat.
// This is a snapshot; chronological increments still require growth proof.
ReferenceThermalAssembly assemble_reference_thermal_volume(
    const std::vector<DepositedSegment>& reference,const std::vector<DepositedSegment>& deposited,
    const std::vector<SectionThermalMaterial>& materials,const ThermalSectionStreamSettings& settings,
    const SectionGrowthStreamLimits& activation_limits = {},
    const std::function<std::optional<SectionOwnership>(double)>& reference_ownership = {},
    const std::function<void()>& checkpoint = {});
SectionThermalGrowthAssembly assemble_reference_thermal_growth(
    const std::vector<DepositedSegment>& reference,const std::vector<DepositedSegment>& previous,
    const std::vector<DepositedSegment>& next,const std::vector<SectionThermalMaterial>& materials,
    const ThermalSectionStreamSettings& settings,const SectionGrowthStreamLimits& growth_limits = {},
    const std::function<std::optional<SectionOwnership>(double)>& reference_ownership = {},
    const std::function<void()>& checkpoint = {});
// Complete next skin/volume and ONLY the certified increment, assembled in a
// single paired traversal. All piece labels come from those same topology
// streams. Added material gets the supplied initial temperatures; this routine
// does not replace or evolve temperatures already present in a history.
SectionThermalGrowthAssembly assemble_local_thermal_growth(
    const std::vector<DepositedSegment>& previous,const std::vector<DepositedSegment>& next,
    const std::vector<SectionThermalMaterial>& materials,const ThermalSectionStreamSettings& settings,
    const SectionGrowthStreamLimits& growth_limits = {},
    const std::function<std::optional<SectionOwnership>(double)>& previous_ownership = {},
    const std::function<std::optional<SectionOwnership>(double)>& next_ownership = {},
    const std::function<void()>& checkpoint = {},
    const std::vector<DepositedSegment>* previous_partitions = nullptr,
    const std::vector<DepositedSegment>* next_partitions = nullptr);

} // namespace Slic3r::Thermoelastic
