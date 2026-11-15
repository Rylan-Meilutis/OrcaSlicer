#pragma once

#include "ThermomechanicalGCode.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace Slic3r::Thermoelastic {
struct LayerMechanicalSample
{
    std::size_t move_index;
    Eigen::Vector3d position_m, displacement_m;
    Eigen::Matrix<double, 6, 1> stress_pa;
    double maximum_principal_pa, von_mises_pa;
};
struct LayerBedSample
{
    Eigen::Vector3d position_m, reaction_n;
    double tributary_area_m2, tensile_traction_pa;
};
struct LayerPreviewLimits
{
    GCodeStudyLimits study;
    std::size_t max_samples = 250000, max_bed_samples = 12000;
};
struct LayerMechanicalPreview
{
    enum class State { Ready, MissingLayer, MissingTiming, UnsupportedInput, BudgetExceeded, Unconverged, SolverFailure };
    State state               = State::UnsupportedInput;
    unsigned layer_id         = 0;
    double observation_time_s = 0.;
    std::vector<LayerMechanicalSample> samples;
    std::vector<LayerBedSample> bed_samples;
    std::size_t unsampled_moves = 0;
    unsigned completed_runs     = 0;
    // Internal bounded diagnostic, not a translated user-facing sentence.
    std::string detail;
    // Qualification applies to integrated measures at this modeled observation.
    // Midpoint samples are not path maxima; local fields/edge traction have NOT
    // received an independent pointwise-convergence or physical validation.
    bool integrated_discretization_qualified = false;
};

// Explicit material/thermal/boundary scenario only. End of the selected
// processor layer's last MODEL deposition command defines the observation;
// later layers/cooldown do not contribute. Layer IDs are processor IDs, not Z.
// Returns stress/displacement and a separate area-normalized bed reaction map,
// never converts stress into a fictitious curl force or adhesion safety score.
// Unsupported geometry, absent data or incomplete qualification publish no map.
// The fine-grid transient backend rejects unsupported geometry and oversized
// cases; the section solver is not substituted as thermal history.
LayerMechanicalPreview preview_mechanical_layer(const GCodeProcessorResult& source,
                                                unsigned layer_id,
                                                const GCodeThermalScenario& scenario,
                                                const LayerPreviewLimits& limits                    = {},
                                                const std::function<void(unsigned, bool)>& progress = {},
                                                const std::function<void()>& checkpoint             = {});
} // namespace Slic3r::Thermoelastic
