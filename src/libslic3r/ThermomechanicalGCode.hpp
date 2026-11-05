#pragma once
#include "ThermomechanicalSimulation.hpp"
#include "ThermalHistory.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticToolpath.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

namespace Slic3r {
struct GCodeProcessorResult;
namespace Thermoelastic {
enum class TimingScenario { RejectModeledIntervals, AcknowledgedModeledIntervals };
struct GCodeThermalScenario {
    TimingScenario timing = TimingScenario::RejectModeledIntervals;
    std::vector<Voxel> mechanical_materials;
    std::vector<ThermalMechanicalMaterial> thermal_materials;
    // Polymer temperature at deposition, explicitly supplied, not nozzle target.
    std::vector<double> deposited_temperature_k;
    std::vector<BoundaryScenario> boundaries;
    Eigen::Vector3d spacing_m = Eigen::Vector3d::Zero();
    double bed_z_m = 0., intermaterial_resistance_m2_k_w = 0., cooldown_s = 0.;
    // Optional snapshot within the modeled print/cooldown timeline. Future
    // deposition is not included in the active thermal or mechanical body.
    std::optional<double> observation_time_s;
    // Explicit scenario assumption. Imported display widths do not establish
    // the physical section; mixed bridge/scarf paths require separate handling.
    BeadSection section = BeadSection::RectangularEnvelope;
    // Explicit approximation only: bounded stationary extrusion immediately
    // at the preceding same-material bead endpoint may leave the nominal
    // geometry/thermal mass unchanged. Missing volume is always reported.
    // Absent by default; vertical extrusion and new isolated deposits reject.
    std::optional<StationaryExtrusionAllowance> stationary_extrusion_allowance;
};
struct GCodeSimulationLimits {
    SimulationLimits simulation;
    std::size_t max_moves = 2000000, max_segments = 250000, max_cell_visits = 2000000;
};
struct GCodeSimulationResult {
    SimulationResult simulation;
    double modeled_print_time_s = 0.;
    double observation_time_s = 0.;
    bool uses_modeled_intervals = true;
    bool uses_nominal_first_hit_cell_activation = true;
    bool physical_deposition_timing_resolved = false;
    // Separate diagnostics: never rescale density or fill gaps to match E.
    // Commanded volume uses available raw E areas, without firmware overrides.
    double commanded_path_volume_m3 = 0.;
    std::size_t segments_without_commanded_volume = 0;
    // Section area times XY length excludes end caps and counts overlaps;
    // occupied volume includes nominal caps but counts overlapping cells once.
    double nominal_section_path_volume_m3 = 0., occupied_volume_m3 = 0.;
    StationaryExtrusionSummary stationary_extrusion;
};
// Bounded scenario evaluation, not a measured temperature/adhesion prediction.
// No implicit material properties, nozzle-temperature substitution or timing
// repair. Exceptions propagate and no partial mechanical field is published.
// An in-progress nominal capsule is clipped at the modeled observation point;
// its end cap is the same nominal footprint assumption, not a measured melt
// front. Cell fractions are transformed back to the original move interval.
GCodeSimulationResult simulate_gcode_deposition(const GCodeProcessorResult &source,
    const GCodeThermalScenario &scenario, const GCodeSimulationLimits &limits = {},
    const std::function<void()> &checkpoint = {});

struct GCodeReleaseRegion {
    // Half-open world XY bounds select bottom-cell centers on the initial
    // raster. Refinements retain that selected physical voxel footprint.
    Eigen::Vector2d minimum_m, maximum_m;
    std::optional<std::size_t> material;
    double fracture_energy_j_m2 = 0., gap_coefficient_w_m2_k = 0., gap_temperature_k = 0.;
};
struct GCodePeelingLimits {
    PeelingSpaceStudyLimits solve;
    std::size_t max_moves = 2000000, max_segments = 250000, max_cell_visits = 2000000;
};
struct GCodePeelingStudy {
    PeelingSpaceStudy solve;
    StationaryExtrusionSummary stationary_extrusion;
};
// Parses/rasterizes once through the same strict preparation as the snapshot
// solver, then evaluates coupled peeling with time and fixed-volume spatial
// refinement. No scalar result is a real adhesion safety verdict. Regions must
// select disjoint, nonempty bottom-cell sets deposited by the first observation.
// scenario.observation_time_s must be absent; the supplied times are authoritative.
// Unresolved extrusion, unsupported geometry/materials and budgets still reject
// the study. This does not qualify the initial bead rasterization itself.
GCodePeelingStudy study_gcode_peeling(const GCodeProcessorResult &source,
    const GCodeThermalScenario &scenario, const std::vector<double> &observation_times_s,
    const std::vector<GCodeReleaseRegion> &regions, const GCodePeelingLimits &limits,
    const std::function<void()> &checkpoint = {});

enum class GCodeStudyState { Unconverged, Converged, BudgetExceeded, SolverFailure };
struct GCodeStudyLimits {
    GCodeSimulationLimits run;
    unsigned spatial_levels = 3; // Initial grid and at least two halvings.
    unsigned max_runs = 8; // Shared-history callers may impose a smaller budget.
    double relative_tolerance = .05;
    double force_absolute_tolerance_n = 1e-8;
    double elastic_energy_absolute_tolerance_j = 1e-12;
    double heat_absolute_tolerance_j = 1e-8;
    double contact_area_absolute_tolerance_m2 = 1e-12;
};
struct GCodeStudyResult {
    GCodeStudyState state = GCodeStudyState::Unconverged;
    unsigned completed_levels = 0, completed_runs = 0;
    std::optional<ThermalHistory::State> failed_thermal_state;
    // Published only after temporal comparison at each grid and two successive
    // spatial comparisons. No partial field masquerades as qualified output.
    std::optional<GCodeSimulationResult> qualified;
};
// Fixed explicit material/boundary/observation scenario. At each grid, compare
// time-step/error-tolerance refinement independently; then halve spatial grid
// spacing and compare integrated mechanical demand, heat and bed area. This
// qualifies discretization at ONE observation, not material data, deposition
// timing, peak print-history demand, local crack traction or adhesion capacity.
// Requires positive local thermal error control. At most eight bounded runs;
// only one simulation field is retained between runs. Cancellation propagates.
GCodeStudyResult study_gcode_deposition(const GCodeProcessorResult &source,
    const GCodeThermalScenario &scenario, const GCodeStudyLimits &limits,
    const std::function<void(unsigned, bool)> &progress = {},
    const std::function<void()> &checkpoint = {});

struct GCodeHistorySample {
    double time_s = 0., tensile_force_n = 0., elastic_energy_j = 0.;
};
struct GCodeHistoryStudy {
    GCodeStudyState state = GCodeStudyState::Unconverged;
    unsigned completed_observations = 0, completed_runs = 0;
    std::optional<std::size_t> failed_observation;
    // Populated only if EVERY observation qualifies. Scalar history only: no
    // accumulation of meshes or temperature fields between observations.
    std::vector<GCodeHistorySample> qualified_samples;
    std::optional<std::size_t> peak_force_sample, peak_energy_sample;
};
// Explicit strictly increasing modeled times, at most 16 observations and 128
// simulation runs in total. Each observation independently qualifies thermal
// integration and spatial resolution. Sampled maxima are NOT bounds on the
// unsampled history or evidence of adhesion safety. scenario.observation_time_s
// must be absent; the supplied times are the sole observation schedule.
GCodeHistoryStudy study_gcode_history(const GCodeProcessorResult &source,
    const GCodeThermalScenario &scenario, const std::vector<double> &times_s,
    const GCodeStudyLimits &limits, unsigned max_total_runs = 64,
    const std::function<void(std::size_t, unsigned, bool)> &progress = {},
    const std::function<void()> &checkpoint = {});
} // namespace Thermoelastic
} // namespace Slic3r
