#include "ThermoelasticSectionPreview.hpp"
#include "ThermoelasticGCodeSections.hpp"
#include "ThermoelasticToolpath.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticSectionStudy.hpp"
#include "ThermoelasticSections.hpp"
#include "ThermoelasticLayers.hpp"
#include "GCode/GCodeProcessor.hpp"
#include "ExtrusionEntity.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cstddef>
#include <exception>
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
SectionLayerPreview preview_section_layer(const GCodeProcessorResult& source,
    unsigned layer_id, const GCodeSectionScenario& scenario,
    std::size_t max_moves, std::size_t max_segments, std::size_t max_samples,
    const std::function<std::optional<SectionOwnership>(double)>& ownership,
    const std::function<void(unsigned, unsigned)>& progress, const std::function<void()>& checkpoint)
{
    if (scenario.observation_time_s || max_samples == 0 || max_samples > 200000)
        throw std::invalid_argument("Invalid static layer preview selection or sample budget.");
    SectionLayerPreview result;
    result.layer_id = layer_id;
    std::exception_ptr callback_error;
    const auto check = [&] {
        try { if (checkpoint) checkpoint(); }
        catch (...) { callback_error = std::current_exception(); throw; }
    };
    const auto report = [&](unsigned a, unsigned b) {
        try { if (progress) progress(a, b); }
        catch (...) { callback_error = std::current_exception(); throw; }
    };
    const auto masks = [&](double z) -> std::optional<SectionOwnership> {
        try {
            auto result = ownership ? ownership(z) : std::nullopt;
            const auto preserve_error = [&](auto deferred) {
                return [&, deferred = std::move(deferred)](const auto& conflict) {
                    try { return deferred(conflict); }
                    catch (...) { callback_error = std::current_exception(); throw; }
                };
            };
            if (result && result->resolve_conflict)
                result->resolve_conflict = preserve_error(std::move(result->resolve_conflict));
            if (result && result->resolve_sampled_conflict)
                result->resolve_sampled_conflict = preserve_error(std::move(result->resolve_sampled_conflict));
            if (result && result->resolve_partition)
                result->resolve_partition = preserve_error(std::move(result->resolve_partition));
            return result;
        }
        catch (...) { callback_error = std::current_exception(); throw; }
    };
    try {
        check();
        auto path = timed_deposition(source, max_moves, max_segments, check);
        if (path.state == ToolpathState::BudgetExceeded) throw ResourceLimit("Static layer path budget exceeded.");
        if (path.state != ToolpathState::Complete) throw std::invalid_argument("Complete modeled deposition timing is unavailable.");
        std::vector<SectionSampleRequest> samples;
        for (const auto& entry : path.deposition) {
            check();
            if (source.moves.at(entry.move_index).layer_id != layer_id || entry.role == erWipeTower ||
                entry.role == erCustom || entry.role == erNone) continue;
            if (entry.geometry.material >= scenario.filament_to_material.size())
                throw std::invalid_argument("Selected layer has an unmapped filament.");
            if (samples.size() == max_samples) throw ResourceLimit("Static layer sample budget exceeded.");
            const auto& bead = entry.geometry;
            Eigen::Vector3d center = .5 * (bead.start_m + bead.end_m);
            center.z() -= .25 * (bead.height_m + bead.end_height_m.value_or(bead.height_m));
            samples.push_back({entry.move_index, center, scenario.filament_to_material[bead.material]});
            result.observation_time_s = std::max(result.observation_time_s, entry.end_s);
        }
        if (samples.empty()) { result.state = SectionLayerPreview::State::MissingLayer; return result; }
        // Do not retain two complete timed-path copies during section assembly.
        path = {};
        auto selected = scenario;
        selected.observation_time_s = result.observation_time_s;
        auto analysis = study_gcode_sections(source, selected, max_moves, max_segments,
            ownership ? std::function<std::optional<SectionOwnership>(double)>(masks) :
                        std::function<std::optional<SectionOwnership>(double)>(), report, check, samples);
        if (callback_error) std::rethrow_exception(callback_error);
        result.basis_history = std::move(analysis.study.last_basis_history);
        result.geometry_history = std::move(analysis.study.history);
        result.geometry_tiles = analysis.study.last_geometry_tiles;
        result.geometry_slabs = analysis.study.last_geometry_slabs;
        result.section_count = analysis.study.section_count;
        result.segment_visits = analysis.study.segment_visits;
        result.boundary_points = analysis.study.boundary_points;
        result.approximate_geometry = scenario.study.initial_z_quantization_m.has_value() ||
            scenario.study.sections.geometry.ownership_boundary_grid_units != 0;
        if (!analysis.study.field || (analysis.study.state != DepositionSectionStudy::State::NominalGeometryAndBasisConverged &&
            analysis.study.state != DepositionSectionStudy::State::ApproximateGeometryAndBasisConverged)) {
            result.state = analysis.study.state == DepositionSectionStudy::State::ResourceLimited ?
                SectionLayerPreview::State::ResourceLimited : SectionLayerPreview::State::Unconverged;
            result.detail = analysis.study.failure_reason.empty() ?
                "Prescribed-strain layer geometry or basis did not converge within the bounded study." : analysis.study.failure_reason;
            return result;
        }
        auto& field = *analysis.study.field;
        result.samples = std::move(field.samples);
        result.unsampled_moves = field.unsampled_points;
        result.elastic_energy_j = field.elastic_energy_j;
        for (const auto& reaction : field.generalized_bed_reaction_n)
            result.tensile_generalized_bed_force_n += std::max(0., -reaction.z());
        result.state = SectionLayerPreview::State::Ready;
    } catch (const ResourceLimit& error) {
        if (callback_error) std::rethrow_exception(callback_error);
        result.state = SectionLayerPreview::State::ResourceLimited;
        result.detail = std::string(error.what()).substr(0, 512);
    } catch (const std::invalid_argument& error) {
        if (callback_error) std::rethrow_exception(callback_error);
        result.state = SectionLayerPreview::State::UnsupportedInput;
        result.detail = std::string(error.what()).substr(0, 512);
    } catch (const std::runtime_error& error) {
        if (callback_error) std::rethrow_exception(callback_error);
        result.state = SectionLayerPreview::State::SolverFailure;
        result.detail = std::string(error.what()).substr(0, 512);
    }
    return result;
}
} // namespace Slic3r::Thermoelastic
