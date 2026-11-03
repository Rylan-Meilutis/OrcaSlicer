#pragma once
#include "ThermoelasticMesh.hpp"
#include "ExtrusionEntity.hpp"
#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <vector>

namespace Slic3r {
struct GCodeProcessorResult;
namespace Thermoelastic {
// Optional precise model-frame Z for authenticated anchored commands. Rejects
// stale/mismatched evidence; display vertices remain the XY and projection
// witnesses. The caller supplies the preceding physical (not seam) move.
std::optional<std::array<double, 2>> authenticated_deposition_z(
    const GCodeProcessorResult &source, std::size_t move, std::size_t previous);
struct TimedDeposition {
    DepositedSegment geometry;
    double start_s = 0., end_s = 0.;
    std::size_t move_index = 0;
    unsigned gcode_id = 0, temperature_targets_id = 0;
    ExtrusionRole role = erNone;
    double fan_percent = 0., nozzle_target_c = 0.;
    // Raw commanded E volume per path length, not firmware M221-adjusted flow
    // or a measured bead area. Zero/missing processor data remains unavailable.
    std::optional<double> commanded_area_m2;
};
enum class ToolpathState { Complete, MissingTiming, InvalidInput, BudgetExceeded };
// Positive extrusion without a resolved deposited volume (missing XY footprint
// or an authenticated collapsed envelope). Preserve its mass indicator and
// timing, but do not invent a line, thickness, dot radius or deposited mesh.
// Authenticated pure-Z BrickCourse joins remain here even though their legacy
// preview/timing classification is Travel. Unmarked travels remain unchanged.
struct UnresolvedSpatialDeposition {
    enum class Reason { UnresolvedXYFootprint, CollapsedEnvelope };
    Eigen::Vector3d start_m, end_m;
    double start_s = 0., end_s = 0., extruded_filament_mm = 0.;
    std::size_t material = 0, move_index = 0;
    unsigned gcode_id = 0;
    ExtrusionRole role = erNone;
    Reason reason = Reason::UnresolvedXYFootprint;
    std::uint32_t source_id = 0, command_id = 0;
};
enum class ToolpathFailure { None, MoveBudget, MissingTotal, InternalVertexDuration,
    InvalidDurationOrPosition, TimeOverflow, MissingExtrusionDuration,
    InvalidBeadDimensions, InvalidFanOrTemperature, UnresolvedXYFootprint,
    SegmentBudget, TotalMismatch, InvalidCommandedArea, InvalidVerticalDeposition };
struct TimedToolpath {
    ToolpathState state = ToolpathState::InvalidInput;
    ToolpathFailure failure = ToolpathFailure::None;
    std::optional<std::size_t> failing_move_index, previous_physical_move_index;
    std::vector<TimedDeposition> deposition;
    std::vector<UnresolvedSpatialDeposition> unresolved_deposition;
    bool spatial_deposition_resolved = true;
    double total_s = 0.;
    // Processor intervals can include attributed heater/tool-change allowances.
    // They are modeled intervals, not measured motion/deposition timestamps.
    bool deposition_phase_timing_resolved = false;
};
struct StationaryExtrusionAllowance {
    double maximum_volume_m3 = 0.;
    double maximum_total_fraction = 0.;
    double maximum_preceding_command_fraction = 0.;
};
struct StationaryExtrusionSummary {
    std::size_t count = 0;
    double commanded_volume_m3 = 0.;
};
// Explicit numerical omission assumptions, not stress-error or bonding bounds.
// Defaults reject collapsed envelopes; callers must opt in within hard ceilings.
struct CollapsedEnvelopeAllowance {
    double maximum_span_volume_m3 = 0., maximum_total_volume_m3 = 0.;
    double maximum_total_fraction = 0., maximum_neighbor_fraction = 0.;
};
struct CollapsedEnvelopeSummary {
    std::size_t count = 0;
    double commanded_volume_m3 = 0.;
};
struct NominalDeposition {
    std::vector<DepositedSegment> segments;
    std::vector<double> retained_fractions;
    double commanded_path_volume_m3 = 0., nominal_section_path_volume_m3 = 0.;
    std::size_t segments_without_commanded_volume = 0;
    StationaryExtrusionSummary stationary_extrusion;
    CollapsedEnvelopeSummary collapsed_envelope;
};
// Shared geometry admission for static and transient studies. Input must be a
// complete timed_deposition result (or ordered role-filtered subset) from the
// same processor source. Retains the
// nominal part of each path up to observation_time_s; later cooldown is allowed.
// Stationary extrusion is rejected unless the explicit, tightly bounded
// same-material preceding-endpoint approximation passes. Missing mass is
// reported, never converted to invented line/dot geometry or density changes.
NominalDeposition resolve_nominal_deposition(const GCodeProcessorResult &source,
    const TimedToolpath &path, double observation_time_s, BeadSection section,
    const std::optional<StationaryExtrusionAllowance> &stationary_allowance = {},
    const std::function<void()> &checkpoint = {},
    const std::optional<CollapsedEnvelopeAllowance> &collapsed_allowance = {});

TimedToolpath timed_deposition(const GCodeProcessorResult &source,
    std::size_t max_moves = 2000000, std::size_t max_segments = 250000,
    const std::function<void()> &checkpoint = {});
} // namespace Thermoelastic
} // namespace Slic3r
