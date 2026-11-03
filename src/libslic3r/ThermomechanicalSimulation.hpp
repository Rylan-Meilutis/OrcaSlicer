#pragma once
#include "ThermalDeposition.hpp"
#include "ThermalHistory.hpp"
#include "ThermoelasticMaterial.hpp"
#include "ThermoelasticReduction.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticPeeling.hpp"
#include "ThermoelasticSolver.hpp"
#include "ThermoelasticThermalStrain.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

namespace Slic3r::Thermoelastic {
struct ThermalMechanicalMaterial {
    double conductivity_w_m_k;
    double volumetric_heat_capacity_j_m3_k;
    double stress_free_temperature_k;
    std::vector<ExpansionInterval> expansion;
    Eigen::Matrix3d local_to_world = Eigen::Matrix3d::Identity();
    // Already in WORLD coordinates (unlike local expansion tensors). Call
    // rotate_elasticity explicitly when supplying material-axis stiffness.
    std::optional<ElasticityMatrix> elasticity_pa;
};
struct CellActivation {
    double time_s;
    double temperature_k;
};
struct BedReleaseEvent {
    double time_s;
    // Stable index into the supplied deposition grid, not an append-order
    // thermal ID or an index into a refined mesh. Whole bottom face only.
    std::size_t voxel;
    // Explicit effective exchange across the released gap. Zero is insulated;
    // neither a gap coefficient nor reservoir temperature is inferred from PWM.
    double coefficient_w_m2_k;
    double temperature_k;
};
struct SurfaceExchange {
    double coefficient_w_m2_k; // Face-surface h; center-to-face conduction is added.
    double temperature_k;
    ThermalHistory::BoundaryKind kind = ThermalHistory::BoundaryKind::Ambient;
    // Bits 0..5 select -X,+X,-Y,+Y,-Z,+Z outward faces. Bed faces
    // remain governed exclusively by the separate bed contact condition.
    unsigned face_mask = 0x3f;
    struct Bounds {
        Eigen::Vector3d minimum_m, maximum_m;
    };
    // Inclusive world-space face-center selection, not fractional face coverage.
    std::optional<Bounds> bounds;
};
struct BoundaryScenario {
    double time_s;
    // Applied to matching currently exposed non-bed faces. These are explicit
    // effective convection coefficients, not fan PWM or inferred airflow.
    std::vector<SurfaceExchange> exposed;
    double bed_coefficient_w_m2_k; // Face-surface contact h, zero means insulated.
    double bed_temperature_k;
};
struct SimulationLimits {
    ThermalHistory::HistoryLimits thermal;
    std::size_t max_schedule_entries = 1000000;
    std::size_t max_mechanical_nodes = 12000;
    struct Reduction {
        unsigned axis = 0;
        // Nested integer strides, each dividing its predecessor. At most five.
        std::vector<unsigned> strides{4, 2, 1};
        double displacement_relative_tolerance = 0.;
        double energy_relative_tolerance = 0.;
        ReductionLimits solve;
    };
    // Explicit opt-in with caller-specified accuracy targets; no Auto adoption.
    std::optional<Reduction> reduction;
};
struct ReductionStudy {
    enum class State { Unreconciled, Stabilized, FineReference, BudgetExceeded };
    State state = State::Unreconciled;
    std::optional<ReducedResult> latest;
    std::vector<unsigned> solved_strides;
    double displacement_relative_change = 0.;
    double energy_relative_change = 0.;
};
struct SimulationResult {
    ThermalHistory::HistoryResult thermal;
    // Maps append-order thermal IDs to the original deposition-grid indices.
    std::vector<std::size_t> active_voxel_indices;
    Mesh mesh;
    std::optional<BedResult> mechanical;
    // Present instead of mechanical for prescribed contact-loss histories.
    // Do not qualify this as a fully bonded field or discard its lift/contact.
    std::optional<ContactBedResult> partial_mechanical;
    // Never converted to local bed traction. Stabilization concerns this fixed
    // fine geometry/basis only, not thermal/spatial or physical convergence.
    std::optional<ReductionStudy> reduction;
};

// Occupied grid cells are counted once; duplicate indices are rejected. An
// activation time and initial temperature are required for each occupied cell.
// Only face-adjacent ACTIVE cells conduct. Inter-material thermal resistance
// is explicit (zero means ideal thermal contact). Mesh mechanics currently
// assumes perfect mechanical bonding at shared material faces.
// Returns no mechanical field on thermal failure. Invalid input, unsupported
// expansion domains or singular mechanics throw; there is no safe fallback.
// Mechanical output is a final elastic snapshot on the original grid, NOT
// stress-free layer activation onto a deformed body or viscoelastic evolution.
// Optional prescribed releases replace bed heat exchange from their event
// times onward and remove adhesion in the final unilateral contact solve.
// Events are irreversible, ordered, unique bottom-cell releases after activation.
// This is one-way prescribed damage, NOT prediction of when peeling begins or
// gap-dependent heat-transfer feedback. Reduction with releases is unsupported.
SimulationResult simulate_deposition_grid(const DepositionGrid &grid,
    const std::vector<CellActivation> &activation,
    const std::vector<ThermalMechanicalMaterial> &materials,
    const std::vector<BoundaryScenario> &boundary_scenarios,
    double intermaterial_resistance_m2_k_w, double end_time_s,
    const SimulationLimits &limits = {}, const std::function<void()> &checkpoint = {},
    const std::vector<BedReleaseEvent> &bed_releases = {},
    const BedContactLimits &contact_limits = {});

struct GridReleasePatch {
    std::vector<std::size_t> voxels;
    double fracture_energy_j_m2;
    double gap_coefficient_w_m2_k;
    double gap_temperature_k;
};
struct PeelingHistoryLimits {
    SimulationLimits simulation;
    BedReleaseEvolutionLimits release;
    unsigned max_observations = 16;
    unsigned max_equilibrium_solves = 256;
};
struct PeelingHistorySample {
    double time_s = 0., released_area_m2 = 0., maximum_lift_m = 0.;
    double elastic_energy_j = 0.;
    std::size_t new_releases = 0;
    double tensile_force_n = 0.;
    // Demand immediately before the observation's release update must not
    // disappear from sampled peaks merely because damage relieves stress.
    double pre_release_force_n = 0., pre_release_energy_j = 0.;
};
struct PeelingHistoryResult {
    // Final thermal state and mesh only. Mechanical fields live exclusively
    // in final_release, including damage accepted at the final observation.
    SimulationResult final_snapshot;
    BedReleaseEvolution final_release;
    std::vector<BedReleaseEvent> releases;
    std::vector<PeelingHistorySample> samples;
    unsigned equilibrium_solves = 0;
    double maximum_accepted_thermal_interval_s = 0.;
};
// Explicit sampled thermal/peeling feedback, not a converged adhesion forecast.
// At each strictly increasing observation, evaluate patch release under the
// current thermal state. Accepted releases change subsequent bed heat exchange
// and adhesion. Replaying the bounded thermal history preserves deposited mass
// and energy without resetting temperatures at observation times. At most 16
// runs are made; equilibrium trials share one global budget. Only the final
// fields are retained. Failures/cancellation throw without partial publication.
// Patches are disjoint bottom-cell sets, all active by the first observation.
// Thermal integration, observation spacing, grid and patch convergence remain
// caller obligations. No stress-free deformed activation, relaxation, inferred
// interface data, or gap-dependent re-contact heat transfer is supplied.
PeelingHistoryResult simulate_peeling_history(const DepositionGrid &grid,
    const std::vector<CellActivation> &activation,
    const std::vector<ThermalMechanicalMaterial> &materials,
    const std::vector<BoundaryScenario> &boundary_scenarios,
    double intermaterial_resistance_m2_k_w, const std::vector<double> &observation_times_s,
    const std::vector<GridReleasePatch> &patches, const PeelingHistoryLimits &limits = {},
    const std::function<void()> &checkpoint = {});

struct PeelingTimeStudyLimits {
    PeelingHistoryLimits run;
    unsigned max_history_runs = 6, max_total_observations = 64, max_total_equilibrium_solves = 256;
    double relative_tolerance = .05;
    double temperature_tolerance_k = .1, release_time_tolerance_s = .1;
    double lift_tolerance_m = 1e-6, area_tolerance_m2 = 1e-12, energy_tolerance_j = 1e-10;
    double force_tolerance_n = 1e-8;
};
struct PeelingTimeStudy {
    unsigned completed_runs = 0, completed_observations = 0, equilibrium_solves = 0;
    // Only present after independent thermal refinement at each of three
    // schedules and two successive midpoint observation refinements agree.
    // This qualifies time discretization on ONE fixed grid/patch/material
    // scenario, not spatial accuracy, actual adhesion or Auto suitability.
    std::optional<PeelingHistoryResult> temporally_qualified;
};
// Six bounded history runs. Compare temperatures, common-time mechanical
// samples, sampled maxima, released cell identities and per-cell release times.
// Refinement must not disguise a different release location with equal area.
// Requires adaptive thermal error control; exceptions and cancellation
// propagate. Exhausted global budgets throw, never return a partial field.
PeelingTimeStudy study_peeling_time(const DepositionGrid &grid,
    const std::vector<CellActivation> &activation,
    const std::vector<ThermalMechanicalMaterial> &materials,
    const std::vector<BoundaryScenario> &boundary_scenarios,
    double intermaterial_resistance_m2_k_w, const std::vector<double> &observation_times_s,
    const std::vector<GridReleasePatch> &patches, const PeelingTimeStudyLimits &limits,
    const std::function<void()> &checkpoint = {});

struct PeelingSpaceStudyLimits {
    PeelingTimeStudyLimits time;
    unsigned max_history_runs = 18, max_total_observations = 192, max_total_equilibrium_solves = 768;
};
struct QualifiedPeelingGrid {
    DepositionGrid grid;
    PeelingHistoryResult history;
    unsigned subdivisions = 1;
};
struct PeelingSpaceStudy {
    unsigned completed_levels = 0, completed_runs = 0, completed_observations = 0, equilibrium_solves = 0;
    std::optional<QualifiedPeelingGrid> discretization_qualified;
};
// Three grids of the SAME occupied volume (1x,2x,4x subdivision on each axis).
// Empty cells remain empty, material/activation is inherited, and each release
// patch retains its physical footprint and fracture work per area. This is
// refinement of a prescribed voxel model, not of the G-code rasterization or
// nominal deposition timing. Each grid must independently qualify in time.
// Both spatial comparisons must agree in sampled mechanics, heat, parent-cell
// average temperature and patch release identity/time. It does not qualify
// local stress singularities, alternative cracks or real material/bed data.
PeelingSpaceStudy study_peeling_space(const DepositionGrid &grid,
    const std::vector<CellActivation> &activation,
    const std::vector<ThermalMechanicalMaterial> &materials,
    const std::vector<BoundaryScenario> &boundary_scenarios,
    double intermaterial_resistance_m2_k_w, const std::vector<double> &observation_times_s,
    const std::vector<GridReleasePatch> &patches, const PeelingSpaceStudyLimits &limits,
    const std::function<void()> &checkpoint = {});
} // namespace Slic3r::Thermoelastic
