#pragma once

#include "ThermoelasticMesh.hpp"
#include "ThermalAirflow.hpp"
#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <functional>
#include <iosfwd>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace Slic3r {
struct GCodeProcessorResult;
class DynamicPrintConfig;
namespace Thermoelastic {

enum class AnalysisState { Converged, UnsupportedMaterial, InvalidThermalInputs,
    BudgetExceeded, Unconverged, InvalidGeometry };

struct AnalysisLimits {
    std::size_t max_nodes = 12000;
    std::size_t max_segments = 250000;
    unsigned max_subdivisions = 4;
    double relative_tolerance = .05;
};

struct ScenarioDemand {
    double tensile_force_n = 0.;
    double elastic_energy_j = 0.;
};

// The caller supplies the matching source stream and known firmware capability.
// Consumed from its current position to EOF; never rewound or reopened here.
struct AirflowTraceInput {
    std::istream &source;
    std::map<unsigned, ThermalHistory::FanRole> channels;
    ThermalHistory::FirmwareFanPolicy policy;
    std::size_t max_bytes = 256 * 1024 * 1024;
    std::size_t max_events = 100000;
};

enum class AirflowTraceState { NotProvided, Complete, BudgetExceeded, ReadError };
struct AirflowTraceSummary {
    AirflowTraceState state = AirflowTraceState::NotProvided;
    ThermalHistory::AirflowState final_state;
    // Conservative: movement commands carrying E may deposit under absolute or
    // relative extrusion. Unknown control is not zero convection.
    std::array<bool, static_cast<unsigned>(ThermalHistory::FanRole::Count)> unresolved_during_motion {};
    std::size_t source_lines = 0, events = 0;
    // Retained only for a complete stream, bounded by AirflowTraceInput.
    // Line IDs are relative to the supplied stream's starting position. The
    // caller must establish a matching complete G-code source before timing
    // these events; retaining them does not couple convection to mechanics.
    std::vector<ThermalHistory::AirflowEvent> decoded_events;
    // Current mechanics uses prescribed static temperatures, not this history.
    bool incorporated_in_thermal_scenarios = false;
};

struct ConfiguredAnalysis {
    AnalysisState state = AnalysisState::InvalidGeometry;
    std::vector<ScenarioDemand> scenarios;
    // IDs are zero-based filament indices from processor extrusion moves.
    std::vector<std::size_t> unsupported_filaments;
    // Exact material/temperature/provenance identity; geometry is excluded so
    // candidate toolpaths can be compared under the same prescribed scenarios.
    std::vector<double> scenario_parameters;
    std::vector<std::string> scenario_sources;
    AirflowTraceSummary airflow;
};

// Explicit small-strain, isotropic study inputs, not inferred filament data.
// Slots correspond to processor filament IDs and must match the configured type.
// Provenance identifies the supplied data/assumption; it is not a certification.
struct MaterialScenario {
    std::string filament_type;
    std::string provenance;
    double young_pa = 0.;
    double poisson = 0.;
    double cte_min = 0., cte_max = 0.;
    double stress_free_temperature_c = 0.;
    double maximum_chamber_temperature_c = 0.;
};

// Published-property, fully bonded, no-relaxation scenarios. A convergence
// result is NOT an adhesion safety margin. Work is bounded in memory and never
// writes a mesh or temporary G-code file. Cancellation exceptions propagate.
ConfiguredAnalysis analyze_configured_gcode(GCodeProcessorResult &gcode,
    const DynamicPrintConfig &config, const AnalysisLimits &limits = {},
    const std::function<void(unsigned, unsigned)> &progress = {},
    const std::function<void()> &checkpoint = {},
    const std::vector<std::optional<MaterialScenario>> &material_scenarios = {},
    const AirflowTraceInput *airflow_input = nullptr);

// Auto may rank a candidate only when every matched scenario improves both
// measures by more than the discretization tolerance. Missing analysis is not
// a zero-risk baseline. This does not authorize a reduction of bed adhesion.
bool mechanically_improves(const ConfiguredAnalysis &baseline,
    const ConfiguredAnalysis &candidate, double relative_tolerance);

// Inputs are explicit scenarios, not properties inferred from a filament name.
// Call only while owning the result (before publishing it to the GUI). Invalid
// or cancelled analyses clear old samples and never publish a partial field.
RefinementResult analyze_gcode(GCodeProcessorResult &gcode, const std::vector<Voxel> &materials,
    const Eigen::Vector3d &spacing_m, double bed_z_m, std::size_t max_nodes,
    unsigned max_subdivisions, double relative_tolerance,
    const std::function<void(unsigned)> &progress = {},
    const std::function<void()> &checkpoint = {});

} // namespace Thermoelastic
} // namespace Slic3r
