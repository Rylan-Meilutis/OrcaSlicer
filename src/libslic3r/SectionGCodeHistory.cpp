#include "SectionGCodeHistory.hpp"
#include "SectionGrowthStream.hpp"
#include "SectionThermalSequence.hpp"
#include "SectionThermalAssembly.hpp"
#include "SectionThermalLoad.hpp"
#include "ThermalDeposition.hpp"
#include "ThermalHistory.hpp"
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticSections.hpp"
#include "ThermoelasticLocalSections.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include "ThermoelasticToolpath.hpp"
#include "ThermoelasticConflictOwnership.hpp"
#include "ThermoelasticExactSections.hpp"
#include "ExtrusionEntity.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
SectionGCodeResult simulate_section_gcode_history(const GCodeProcessorResult& source,
                                                  const SectionGCodeScenario& scenario,
                                                  const SectionGCodeLimits& limits,
                                                  const std::function<void()>& checkpoint)
{
    const auto count = scenario.thermal_materials.size();
    if (limits.streaming_geometry.partition_increment_m != 0. && !scenario.stream_tile_size_m)
        throw std::invalid_argument("Section partition approximation requires streamed geometry.");
    if (scenario.stream_vertical_normal_bins && (!scenario.stream_tile_size_m ||
        scenario.stream_vertical_normal_bins < 8 || scenario.stream_vertical_normal_bins > 4096 ||
        scenario.stream_vertical_normal_bins % 4))
        throw std::invalid_argument("Invalid streamed thermal surface angular resolution.");
    if ((scenario.stream_ownership || scenario.use_native_ownership || scenario.use_deposited_ownership || scenario.use_fixed_interface_reference) && !scenario.stream_tile_size_m)
        throw std::invalid_argument("Native prefix ownership requires streamed geometry.");
    if (scenario.stream_ownership && scenario.use_native_ownership)
        throw std::invalid_argument("Native evidence and custom prefix ownership cannot be combined.");
    if ((scenario.use_deposited_ownership && (scenario.use_native_ownership || scenario.stream_ownership ||
            limits.ownership.deposited_overlap_policy != DepositedOverlapPolicy::FirstDepositionWins)) ||
        (!scenario.use_deposited_ownership && limits.ownership.deposited_overlap_policy != DepositedOverlapPolicy::None))
        throw std::invalid_argument("Deposited overlap ownership requires its separate explicit scenario.");
    if (!scenario.acknowledge_modeled_timing || (!scenario.stream_tile_size_m && scenario.section != BeadSection::RectangularEnvelope) ||
        (scenario.section != BeadSection::RectangularEnvelope && scenario.section != BeadSection::RoundedRectangle) || !count ||
        count > 256 || scenario.expansion.size() != count || scenario.mechanical_materials.size() != count ||
        scenario.filament_to_material.empty() || scenario.filament_to_material.size() > 256 || !scenario.deposition_temperatures_k ||
        !scenario.exposure || scenario.activation_times_s.empty() || scenario.activation_times_s.size() > 10000 ||
        scenario.exposure_times_s.size() > 10000 - scenario.activation_times_s.size() ||
        !std::isfinite(scenario.cooldown_s) || scenario.cooldown_s < 0. || scenario.cooldown_s > 86400. || !limits.max_prefix_visits ||
        limits.max_prefix_visits > 200000000)
        throw std::invalid_argument("Invalid section G-code thermal scenario or timing acknowledgement.");
    for (auto material : scenario.filament_to_material)
        if (material >= count)
            throw std::invalid_argument("Missing section G-code material mapping.");
    auto path = timed_deposition(source, limits.max_moves, limits.max_segments, checkpoint);
    if (path.state == ToolpathState::BudgetExceeded)
        throw ResourceLimit("Section G-code source budget exceeded.");
    if (path.state != ToolpathState::Complete)
        throw std::invalid_argument("Incomplete section G-code source timing.");
    double previous = 0.;
    for (double time : scenario.activation_times_s) {
        if (!std::isfinite(time) || time <= previous || time > path.total_s)
            throw std::invalid_argument("Invalid section G-code activation schedule.");
        previous = time;
    }
    if (previous != path.total_s)
        throw std::invalid_argument("Section activation schedule omits the end of the print.");
    previous = -1.;
    for (double time : scenario.exposure_times_s) {
        if (!std::isfinite(time) || time < 0. || time <= previous || time > path.total_s + scenario.cooldown_s)
            throw std::invalid_argument("Invalid section G-code exposure schedule.");
        previous = time;
    }
    const auto excluded = [](const auto& entry) { return entry.role == erWipeTower || entry.role == erCustom || entry.role == erNone; };
    SectionGCodeResult result;
    result.thermal_quadrature=limits.assembly.quadrature;
    result.uses_fixed_interface_reference=scenario.use_fixed_interface_reference;
    result.allows_single_sided_native_boundary=scenario.use_native_ownership && limits.ownership.allow_single_sided_native_boundary;
    result.allows_dual_native_boundary=scenario.use_native_ownership && limits.ownership.allow_dual_native_boundary;
    result.native_interface_variant=scenario.use_native_ownership ? limits.ownership.interface_variant : NativeInterfaceVariant::Strict;
    result.deposited_overlap_policy=scenario.use_deposited_ownership ? limits.ownership.deposited_overlap_policy : DepositedOverlapPolicy::None;
    result.vertical_normal_bins=scenario.stream_vertical_normal_bins;
    result.section_partition_increment_m=limits.streaming_geometry.partition_increment_m;
    result.modeled_print_time_s = path.total_s;
    const auto before           = path.deposition.size();
    path.deposition.erase(std::remove_if(path.deposition.begin(), path.deposition.end(), excluded), path.deposition.end());
    result.excluded_segments = before - path.deposition.size();
    path.unresolved_deposition.erase(std::remove_if(path.unresolved_deposition.begin(), path.unresolved_deposition.end(), excluded),
                                     path.unresolved_deposition.end());
    for (const auto& entry : path.deposition)
        if (entry.role == erMixed)
            throw std::invalid_argument("Mixed-role G-code has no resolved section ownership.");
    const auto ownership_at = [&](double time) {
        SectionThermalStreamSequence::Ownership ownership;
        if (scenario.stream_ownership)
            ownership = [provider=scenario.stream_ownership,time](double z) { return provider(time,z); };
        if (scenario.use_native_ownership || scenario.use_deposited_ownership) {
            SampledConflictOwnershipResolver resolver;
            auto ownership_limits = limits.ownership;
            ownership_limits.max_moves = limits.max_moves;
            ownership_limits.max_segments = limits.max_segments;
            // Authentication must describe the actual sampler's geometry,
            // including half-native planes used by one-nanometre intervals.
            ownership_limits.geometry.arc_tolerance_m = limits.streaming_geometry.geometry.arc_tolerance_m;
            ownership_limits.geometry.half_native_section_planes = limits.streaming_geometry.geometry.half_native_section_planes;
            prepare_native_model_ownership(source,scenario.section,scenario.bed_z_m,scenario.filament_to_material,
                time,scenario.stationary_allowance,ownership_limits,checkpoint,false,&resolver,scenario.collapsed_envelope_allowance);
            const bool partition=ownership_limits.interface_variant!=NativeInterfaceVariant::Strict;
            const bool deposited=scenario.use_deposited_ownership;
            ownership=[resolver=std::move(resolver),partition,deposited](double z) -> std::optional<SectionOwnership> {
                SectionOwnership request{z,{}};
                if (deposited)
                    request.resolve_deposited_partition=[resolver](const SampledSectionConflict& conflict) { return resolver(conflict).deposited_partition; };
                else if (partition)
                    request.resolve_partition=[resolver](const SampledSectionConflict& conflict) { return resolver(conflict).partition; };
                else
                    request.resolve_sampled_conflict=[resolver](const SampledSectionConflict& conflict) { return resolver(conflict).regions; };
                return request;
            };
        }
        return ownership;
    };
    const auto map_materials = [&](auto& segments) {
        for (std::size_t index=0;index<segments.size();++index) {
            auto& segment=segments[index];
            if (segment.material>=scenario.filament_to_material.size())
                throw std::invalid_argument("Missing section G-code filament mapping.");
            segment.material=scenario.filament_to_material[segment.material];
            segment.section_source_index=index;
        }
    };
    std::vector<DepositedSegment> reference;
    SectionThermalStreamSequence::Ownership reference_ownership;
    std::size_t reference_visits=0;
    if (scenario.use_fixed_interface_reference) {
        reference_visits=path.deposition.size()+path.unresolved_deposition.size();
        if (reference_visits>limits.max_prefix_visits) throw ResourceLimit("Fixed reference source budget exceeded.");
        reference=resolve_nominal_deposition(source,path,path.total_s,scenario.section,scenario.stationary_allowance,checkpoint,
            scenario.collapsed_envelope_allowance).segments;
        map_materials(reference);
        reference_ownership=ownership_at(path.total_s);
        if (reference.empty()) throw std::invalid_argument("No part geometry for fixed reference history.");
    }
    SectionThermalSequence sequence(scenario.thermal_materials, scenario.origin_m, scenario.spacing_m, limits.assembly, limits.growth,
                                    scenario.activation_times_s.size());
    std::unique_ptr<SectionThermalStreamSequence> streamed;
    auto reference_growth=limits.streaming_growth;
    reference_growth.exact_reference_sources=scenario.use_fixed_interface_reference && scenario.use_deposited_ownership;
    reference_growth.reuse_complete_deposited_reference=reference_growth.exact_reference_sources;
    auto reference_topology=limits.topology;
    reference_topology.exact_interior_components=reference_growth.exact_reference_sources;
    if (scenario.stream_tile_size_m) {
        ThermalSectionStreamSettings settings;
        settings.origin_m               = scenario.origin_m;
        settings.spacing_m              = scenario.spacing_m;
        settings.tile_size_m            = *scenario.stream_tile_size_m;
        settings.bed_z_m                = scenario.bed_z_m;
        settings.maximum_height_m       = scenario.stream_section_height_m;
        settings.geometry               = limits.streaming_geometry;
        settings.topology               = reference_topology;
        settings.streaming_topology     = limits.streaming_topology;
        settings.assembly               = limits.assembly;
        settings.skin                   = limits.streaming_skin;
        settings.vertical_normal_bins   = scenario.stream_vertical_normal_bins;
        settings.max_processed_points   = limits.streaming_thermal_points;
        settings.max_work               = limits.streaming_thermal_work;
        settings.max_provisional_groups = limits.streaming_thermal_groups;
        streamed = std::make_unique<SectionThermalStreamSequence>(scenario.thermal_materials, settings, reference_growth,
                                                                  scenario.activation_times_s.size(),reference,reference_ownership);
    }
    std::vector<DepositedSegment> final_segments;
    SectionThermalStreamSequence::Ownership final_ownership;
    std::optional<SectionTopology> final_geometry;
    std::size_t cursor = 0, exposure_cursor = 0, visits = reference_visits;
    result.thermal = ThermalHistory::evolve_generated(
        [&](ThermalHistory::DepositionEvent& event) {
            while (cursor < scenario.activation_times_s.size() || exposure_cursor < scenario.exposure_times_s.size()) {
                if (checkpoint)
                    checkpoint();
                if (exposure_cursor < scenario.exposure_times_s.size() &&
                    (cursor == scenario.activation_times_s.size() ||
                     scenario.exposure_times_s[exposure_cursor] < scenario.activation_times_s[cursor])) {
                    const double time = scenario.exposure_times_s[exposure_cursor++];
                    if (!result.activated_snapshots)
                        continue; // No material yet; activation evaluates the full current law.
                    const auto exposure = [&](const SectionThermalAssembly& assembly) { return scenario.exposure(time, assembly); };
                    event = streamed ? streamed->replace_exposure(time, exposure, checkpoint) : sequence.replace_exposure(time, exposure, checkpoint);
                    ++result.exposure_changes;
                    return true;
                }
                const auto scan = path.deposition.size() + path.unresolved_deposition.size();
                if (scan > limits.max_prefix_visits - visits)
                    throw ResourceLimit("Section G-code prefix traversal budget exceeded.");
                visits += scan;
                const double time = scenario.activation_times_s[cursor++];
                if (exposure_cursor < scenario.exposure_times_s.size() && scenario.exposure_times_s[exposure_cursor] == time)
                    ++exposure_cursor;
                auto prefix = resolve_nominal_deposition(source, path, time, scenario.section, scenario.stationary_allowance, checkpoint,
                    scenario.collapsed_envelope_allowance);
                result.stationary_extrusion = prefix.stationary_extrusion;
                result.collapsed_envelope = prefix.collapsed_envelope;
                if (prefix.segments.empty())
                    continue;
                map_materials(prefix.segments);
                const auto exposure = [&](const SectionThermalAssembly& assembly) { return scenario.exposure(time, assembly); };
                if (streamed) {
                    auto ownership = scenario.use_fixed_interface_reference ? SectionThermalStreamSequence::Ownership{} : ownership_at(time);
                    event = streamed->append(prefix.segments, time, scenario.deposition_temperatures_k(time), exposure, ownership, nullptr,
                                             checkpoint);
                    if (cursor == scenario.activation_times_s.size()) {
                        final_segments = std::move(prefix.segments);
                        final_ownership = std::move(ownership);
                    }
                } else {
                    const auto layers = layered_deposition(prefix.segments, scenario.bed_z_m, limits.geometry, checkpoint);
                    auto geometry     = qualify_section_topology(layers.slabs, scenario.bed_z_m, limits.topology, checkpoint);
                    event             = sequence.append(geometry, time, scenario.deposition_temperatures_k(time), exposure, checkpoint);
                    final_geometry    = std::move(geometry);
                }
                ++result.activated_snapshots;
                return true;
            }
            return false;
        },
        path.total_s + scenario.cooldown_s, limits.history, checkpoint);
    if (result.thermal.state != ThermalHistory::State::Converged)
        return result;
    if ((!streamed && !final_geometry) || (streamed && final_segments.empty()))
        throw std::invalid_argument("No part geometry in section G-code history.");
    const auto& assembly      = streamed ? streamed->assembly() : sequence.assembly();
    result.occupied_volume_m3 = assembly.occupied_volume_m3;
    if (scenario.use_fixed_interface_reference)
        result.excluded_nominal_volume_m3=streamed->excluded_nominal_volume_m3();
    const auto load           = section_thermal_load(assembly, result.thermal, scenario.expansion, checkpoint);
    if (streamed) {
        LocalSectionStreamSettings settings;
        settings.origin_m           = scenario.origin_m;
        settings.basis_spacing_m    = scenario.spacing_m;
        settings.tile_size_m        = *scenario.stream_tile_size_m;
        settings.bed_z_m            = scenario.bed_z_m;
        settings.maximum_height_m   = scenario.stream_section_height_m;
        settings.geometry           = limits.streaming_geometry;
        settings.topology           = reference_topology;
        settings.streaming_topology = limits.streaming_topology;
        settings.mechanics          = limits.mechanics;
        settings.integration        = limits.streaming_integration;
        settings.load               = load;
        if (scenario.use_fixed_interface_reference)
            result.mechanical=solve_reference_deposition_sections(reference,final_segments,scenario.mechanical_materials,
                settings,reference_growth,reference_ownership,checkpoint).field;
        else
            result.mechanical = solve_local_deposition_sections(final_segments, scenario.mechanical_materials, settings, final_ownership, checkpoint).field;
    } else
        result.mechanical = solve_deposited_sections(*final_geometry, scenario.mechanical_materials, scenario.origin_m, scenario.spacing_m,
                                                     limits.mechanics, checkpoint, {}, {}, {}, load);
    return result;
}
} // namespace Slic3r::Thermoelastic
