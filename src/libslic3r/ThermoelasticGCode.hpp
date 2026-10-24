#pragma once

#include "ThermoelasticMesh.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <functional>
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

struct ConfiguredAnalysis {
    AnalysisState state = AnalysisState::InvalidGeometry;
    std::vector<ScenarioDemand> scenarios;
    // IDs are zero-based filament indices from processor extrusion moves.
    std::vector<std::size_t> unsupported_filaments;
};

// Published-property, fully bonded, no-relaxation scenarios. A convergence
// result is NOT an adhesion safety margin. Work is bounded in memory and never
// writes a mesh or temporary G-code file. Cancellation exceptions propagate.
ConfiguredAnalysis analyze_configured_gcode(GCodeProcessorResult &gcode,
    const DynamicPrintConfig &config, const AnalysisLimits &limits = {},
    const std::function<void(unsigned, unsigned)> &progress = {},
    const std::function<void()> &checkpoint = {});

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
