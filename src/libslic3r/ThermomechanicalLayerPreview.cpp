#include "ThermomechanicalLayerPreview.hpp"
#include "ThermomechanicalGCode.hpp"
#include "ThermoelasticToolpath.hpp"
#include "ThermoelasticField.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermomechanicalSimulation.hpp"
#include "ThermoelasticSolver.hpp"
#include "GCode/GCodeProcessor.hpp"
#include "ExtrusionEntity.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>

namespace Slic3r::Thermoelastic {
namespace {
bool model_role(ExtrusionRole role) { return role != erNone && role != erCustom && role != erWipeTower; }
} // namespace

LayerMechanicalPreview preview_mechanical_layer(const GCodeProcessorResult& source,
                                                unsigned layer_id,
                                                const GCodeThermalScenario& scenario,
                                                const LayerPreviewLimits& limits,
                                                const std::function<void(unsigned, bool)>& progress,
                                                const std::function<void()>& checkpoint)
{
    LayerMechanicalPreview output;
    output.layer_id = layer_id;
    if (scenario.observation_time_s || limits.max_samples == 0 || limits.max_samples > 2000000 || limits.max_bed_samples == 0 ||
        limits.max_bed_samples > 60000)
        throw std::invalid_argument("Invalid selected-layer mechanical preview limits or observation.");
    std::exception_ptr callback_failure;
    const auto check = [&] {
        if (checkpoint)
            try {
                checkpoint();
            } catch (...) {
                callback_failure = std::current_exception();
                throw;
            }
    };
    const auto report = [&](unsigned level, bool finer) {
        if (progress)
            try {
                progress(level, finer);
            } catch (...) {
                callback_failure = std::current_exception();
                throw;
            }
    };
    try {
        check();
        const auto path = timed_deposition(source, limits.study.run.max_moves, limits.study.run.max_segments, check);
        if (path.state != ToolpathState::Complete) {
            output.state  = path.state == ToolpathState::BudgetExceeded ? LayerMechanicalPreview::State::BudgetExceeded :
                                                                          LayerMechanicalPreview::State::MissingTiming;
            output.detail = "Complete modeled deposition timing is unavailable.";
            return output;
        }
        bool found = false;
        for (const auto& entry : path.deposition) {
            check();
            if (model_role(entry.role) && source.moves.at(entry.move_index).layer_id == layer_id) {
                found                     = true;
                output.observation_time_s = std::max(output.observation_time_s, entry.end_s);
            }
        }
        if (!found) {
            output.state  = LayerMechanicalPreview::State::MissingLayer;
            output.detail = "Selected processor layer has no resolved model deposition.";
            return output;
        }
        std::size_t sample_count = 0;
        for (const auto& entry : path.deposition)
            if (model_role(entry.role) && entry.end_s <= output.observation_time_s)
                ++sample_count;
        if (sample_count > limits.max_samples)
            throw ResourceLimit("Selected-layer sample budget exceeded.");
        auto observation               = scenario;
        observation.observation_time_s = output.observation_time_s;
        auto study                     = study_gcode_deposition(source, observation, limits.study, report, check);
        if (callback_failure)
            std::rethrow_exception(callback_failure);
        output.completed_runs = study.completed_runs;
        if (study.state != GCodeStudyState::Converged || !study.qualified || !study.qualified->simulation.mechanical) {
            output.state  = study.state == GCodeStudyState::BudgetExceeded ? LayerMechanicalPreview::State::BudgetExceeded :
                            study.state == GCodeStudyState::Unconverged    ? LayerMechanicalPreview::State::Unconverged :
                                                                             LayerMechanicalPreview::State::SolverFailure;
            output.detail = "Selected-layer thermal and spatial discretization did not qualify.";
            return output;
        }
        auto& simulation = study.qualified->simulation;
        auto& bed        = *simulation.mechanical;
        if (bed.contact_area_m2.size() != simulation.mesh.nodes_m.size() ||
            bed.tensile_traction_pa.size() != simulation.mesh.nodes_m.size() ||
            bed.equilibrium.reaction_n.size() != simulation.mesh.nodes_m.size())
            throw std::runtime_error("Mechanical bed field dimensions are inconsistent.");
        LayerMechanicalPreview ready;
        ready.layer_id           = layer_id;
        ready.observation_time_s = output.observation_time_s;
        ready.completed_runs     = output.completed_runs;
        for (std::size_t i = 0; i < simulation.mesh.nodes_m.size(); ++i) {
            check();
            if (bed.contact_area_m2[i] <= 0.)
                continue;
            if (ready.bed_samples.size() == limits.max_bed_samples)
                throw ResourceLimit("Selected-layer bed map budget exceeded.");
            ready.bed_samples.push_back(
                {simulation.mesh.nodes_m[i], bed.equilibrium.reaction_n[i], bed.contact_area_m2[i], bed.tensile_traction_pa[i]});
        }
        const Field field(std::move(simulation.mesh.nodes_m), std::move(simulation.mesh.elements), std::move(bed.equilibrium));
        ready.samples.reserve(sample_count);
        for (const auto& entry : path.deposition) {
            check();
            if (!model_role(entry.role) || entry.end_s > ready.observation_time_s)
                continue;
            const auto& bead       = entry.geometry;
            Eigen::Vector3d center = .5 * (bead.start_m + bead.end_m);
            center.z() -= .25 * (bead.height_m + bead.end_height_m.value_or(bead.height_m));
            const auto sample = field.sample(center);
            if (!sample) {
                ++ready.unsampled_moves;
                continue;
            }
            if (!sample->stress_pa.allFinite() || !sample->displacement_m.allFinite() || !std::isfinite(sample->maximum_principal_pa) ||
                !std::isfinite(sample->von_mises_pa))
                throw std::runtime_error("Selected-layer mechanical sample is not finite.");
            ready.samples.push_back(
                {entry.move_index, center, sample->displacement_m, sample->stress_pa, sample->maximum_principal_pa, sample->von_mises_pa});
        }
        ready.state                               = LayerMechanicalPreview::State::Ready;
        ready.integrated_discretization_qualified = true;
        return ready;
    } catch (const ResourceLimit& error) {
        if (callback_failure)
            std::rethrow_exception(callback_failure);
        output.state  = LayerMechanicalPreview::State::BudgetExceeded;
        output.detail = std::string(error.what()).substr(0, 512);
    } catch (const std::invalid_argument& error) {
        if (callback_failure)
            std::rethrow_exception(callback_failure);
        output.state  = LayerMechanicalPreview::State::UnsupportedInput;
        output.detail = std::string(error.what()).substr(0, 512);
    } catch (const std::runtime_error& error) {
        if (callback_failure)
            std::rethrow_exception(callback_failure);
        output.state  = LayerMechanicalPreview::State::SolverFailure;
        output.detail = std::string(error.what()).substr(0, 512);
    }
    return output;
}
} // namespace Slic3r::Thermoelastic
