#include "ThermoelasticLocalSections.hpp"
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticMaterial.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticSections.hpp"
#include "ThermoelasticSectionStudy.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
SectionPartitionEnvelopes make_section_partition_envelopes(const std::vector<DepositedSegment>& source, double bed, double increment,
    const std::function<void()>& checkpoint)
{
    if (source.empty() || source.size() > 2000000 || !std::isfinite(bed) || std::abs(bed) > 10. ||
        !std::isfinite(increment) || increment < 2e-9 || increment > .001)
        throw std::invalid_argument("Invalid explicit Z geometry quantization.");
    const auto grid = std::llround(increment * 1e9), bed_native = std::llround(bed * 1e9);
    SectionPartitionEnvelopes result;
    result.segments.reserve(source.size());
    std::array<std::set<long long>, 2> original, changed;
    const auto coordinate = [](double z) {
        if (!std::isfinite(z) || std::abs(z) > 10.) throw std::invalid_argument("Invalid bead quantization coordinate.");
        return std::llround(z * 1e9);
    };
    const auto snap = [&](long long native, bool top) {
        const auto delta = native - bed_native;
        auto cell = delta / grid - (delta % grid < 0 ? 1 : 0);
        if (top && delta % grid != 0) ++cell;
        const auto mapped = bed_native + cell * grid;
        original[top].insert(native);
        changed[top].insert(mapped);
        if (original[0].size() + original[1].size() > 200000 || changed[0].size() + changed[1].size() > 200000)
            throw ResourceLimit("Explicit Z geometry endpoint budget exceeded.");
        result.maximum_shift_m = std::max(result.maximum_shift_m, std::abs(double(mapped - native)) * 1e-9);
        return mapped;
    };
    for (const auto& bead : source) {
        if (checkpoint) checkpoint();
        if (!bead.start_m.allFinite() || !bead.end_m.allFinite() || !std::isfinite(bead.height_m) || bead.height_m < 0. ||
            !std::isfinite(bead.end_height_m.value_or(bead.height_m)) || bead.end_height_m.value_or(bead.height_m) < 0.)
            throw std::invalid_argument("Invalid bead dimensions for explicit Z approximation.");
        const auto original_a = coordinate(bead.start_m.z()), original_b = coordinate(bead.end_m.z());
        const auto h0 = coordinate(bead.height_m), h1 = coordinate(bead.end_height_m.value_or(bead.height_m));
        if (original_a - h0 < bed_native || original_b - h1 < bed_native)
            throw std::invalid_argument("Section partition source intersects the bed.");
        if (h0 == 0 && h1 == 0) throw std::invalid_argument("Z approximation cannot create an originally empty bead.");
        const auto a = snap(original_a, true), b = snap(original_b, true);
        const auto lo_a = snap(original_a - h0, false), lo_b = snap(original_b - h1, false);
        if ((original_a - h0 > bed_native && lo_a == bed_native) ||
            (original_b - h1 > bed_native && lo_b == bed_native))
            throw std::invalid_argument("Z approximation would introduce a new bed contact.");
        if (a < lo_a || b < lo_b || (a == lo_a && b == lo_b))
            throw std::invalid_argument("Explicit Z approximation would collapse a complete deposited bead.");
        auto moved = bead;
        moved.start_m.z() = double(a) * 1e-9;
        moved.end_m.z() = double(b) * 1e-9;
        moved.height_m = double(a - lo_a) * 1e-9;
        moved.end_height_m = double(b - lo_b) * 1e-9;
        result.expanded_zero_endpoints += (h0 == 0 && a > lo_a ? 1 : 0) + (h1 == 0 && b > lo_b ? 1 : 0);
        result.segments.push_back(std::move(moved));
    }
    result.collapsed_intervals = original[0].size() - changed[0].size() + original[1].size() - changed[1].size();
    return result;
}
DepositionSectionStudy study_local_deposition_sections(const std::vector<DepositedSegment>& segments,
    const std::vector<SectionMaterial>& materials, const Eigen::Vector3d& origin,
    const Eigen::Vector3d& spacing, const Eigen::Vector2d& tile_size,
    double bed_z, const SectionStudySettings& settings,
    const std::function<std::optional<SectionOwnership>(double)>& ownership,
    const std::function<void(unsigned, unsigned)>& progress,
    const std::function<void()>& checkpoint, const std::vector<SectionSampleRequest>& samples)
{
    if (segments.empty() || segments.size() > 2000000 || materials.empty() || materials.size() > 256 ||
        !settings.max_geometry_levels || settings.max_geometry_levels > 4 ||
        !settings.max_basis_levels || settings.max_basis_levels > 5 || settings.max_initial_basis_probes > 3 ||
        !std::isfinite(settings.max_section_height_m) || settings.max_section_height_m < 2e-9 ||
        !std::isfinite(settings.relative_tolerance) || settings.relative_tolerance <= 0. || settings.relative_tolerance > 1. ||
        !settings.max_total_sections || settings.max_total_sections > 40000 ||
        !settings.max_total_segment_visits || settings.max_total_segment_visits > 80000000 ||
        !settings.max_total_boundary_points || settings.max_total_boundary_points > 8000000 ||
        !settings.max_total_work || settings.max_total_work > 1000000000)
        throw std::invalid_argument("Invalid local section study settings.");
    double maximum_energy_density = 0., maximum_eigenstress = 0.;
    for (const auto& material : materials) {
        if (!material.eigenstrain.allFinite() || !material.eigenstrain_gradient_per_m.allFinite())
            throw std::invalid_argument("Invalid local study eigenstrain.");
        const auto stress = (validated_elasticity(material.elasticity_pa) * material.eigenstrain).eval();
        maximum_energy_density = std::max(maximum_energy_density, .5 * material.eigenstrain.dot(stress));
        maximum_eigenstress = std::max(maximum_eigenstress, stress.norm());
    }
    for (const auto& bead : segments)
        if (bead.material >= materials.size()) throw std::invalid_argument("Unmapped local deposition material.");
    DepositionSectionStudy result;
    bool callback_failed = false;
    const auto check = [&] {
        try { if (checkpoint) checkpoint(); }
        catch (...) { callback_failed = true; throw; }
    };
    const auto charge_work = [&](std::size_t n) {
        if (n > settings.max_total_work - result.reserved_work)
            throw ResourceLimit("Local section study work budget exceeded.");
        result.reserved_work += n;
    };
    const auto stable = [&](double a, double b, double scale) {
        return std::abs(a - b) <= settings.relative_tolerance * std::max(std::abs(a), std::abs(b)) +
            128. * std::numeric_limits<double>::epsilon() * scale;
    };
    unsigned consecutive = 0;
    LocalSectionStatistics current_geometry;
    std::size_t previous_sections = 0, previous_points = 0;
    try {
        for (unsigned level = 0; level < settings.max_geometry_levels; ++level) {
            check();
            current_geometry = {};
            previous_sections = result.section_count;
            previous_points = result.boundary_points;
            LocalSectionLimits limits;
            limits.geometry = settings.sections.geometry;
            limits.max_sections = settings.max_total_sections - result.section_count;
            limits.max_segment_visits = settings.max_total_segment_visits - result.segment_visits;
            limits.max_boundary_points = settings.max_total_boundary_points - result.boundary_points;
            if (!limits.max_sections || !limits.max_segment_visits || !limits.max_boundary_points)
                throw ResourceLimit("Local section study total geometry budget exhausted.");
            const double height = settings.max_section_height_m / double(1u << level);
            if (height < 2e-9) throw ResourceLimit("Local section refinement exceeds native resolution.");
            charge_work(limits.max_index_work + limits.max_segment_visits + settings.topology.max_component_pairs);
            limits.max_ownership_work = settings.max_total_work - result.reserved_work;
            if (!limits.max_ownership_work) throw ResourceLimit("Local section study ownership budget exhausted.");
            const auto owner = [&](double z) -> std::optional<SectionOwnership> {
                try {
                    auto owner = ownership ? ownership(z) : std::nullopt;
                    if (owner && owner->resolve_conflict) {
                        const auto deferred = std::move(owner->resolve_conflict);
                        owner->resolve_conflict = [&, deferred](const LayerMaterialOverlap& conflict) {
                            try { return deferred(conflict); }
                            catch (...) { callback_failed = true; throw; }
                        };
                    }
                    if (owner && owner->resolve_sampled_conflict) {
                        const auto deferred = std::move(owner->resolve_sampled_conflict);
                        owner->resolve_sampled_conflict = [&, deferred](const SampledSectionConflict& conflict) {
                            try { return deferred(conflict); }
                            catch (...) { callback_failed = true; throw; }
                        };
                    }
                    if (owner && owner->resolve_partition) {
                        const auto deferred = std::move(owner->resolve_partition);
                        owner->resolve_partition = [&, deferred](const SampledSectionConflict& conflict) {
                            try { return deferred(conflict); }
                            catch (...) { callback_failed = true; throw; }
                        };
                    }
                    return owner;
                }
                catch (...) { callback_failed = true; throw; }
            };
            std::optional<SectionPartitionEnvelopes> quantized;
            const double quantization = settings.initial_z_quantization_m ?
                *settings.initial_z_quantization_m / double(1u << level) : 0.;
            if (settings.initial_z_quantization_m)
                quantized = make_section_partition_envelopes(segments, bed_z, quantization, check);
            auto geometry = local_deposition_sections(segments,
                origin.head<2>(), tile_size, bed_z, settings.max_section_height_m,
                limits, ownership ? std::function<std::optional<SectionOwnership>(double)>(owner) :
                    std::function<std::optional<SectionOwnership>(double)>{}, check, 1u << level, &current_geometry,
                quantized ? &quantized->segments : nullptr);
            result.section_count += geometry.sections;
            result.segment_visits += geometry.segment_visits;
            result.boundary_points += geometry.boundary_points;
            charge_work(geometry.ownership_work);
            result.last_geometry_tiles = geometry.tiles.size();
            result.last_geometry_slabs = 0;
            for (const auto& tile : geometry.tiles) result.last_geometry_slabs += tile.slabs.size();
            result.piecewise_constant_geometry = geometry.piecewise_constant && !quantized;
            const double maximum_shift = quantized ? quantized->maximum_shift_m : 0.;
            const std::size_t collapsed = quantized ? quantized->collapsed_intervals : 0;
            const std::size_t expanded_zero = quantized ? quantized->expanded_zero_endpoints : 0;
            quantized.reset();
            const std::size_t section_count = geometry.sections;
            const auto topology = qualify_local_section_topology(geometry.tiles, bed_z, settings.topology, check);
            std::vector<LocalSectionTile>().swap(geometry.tiles);
            unsigned run = 0;
            const auto report = [&](unsigned) {
                charge_work(settings.mechanics.max_work);
                try { if (progress) progress(level, run++); }
                catch (...) { callback_failed = true; throw; }
            };
            Eigen::Vector3d selected_spacing = spacing;
            result.last_initial_basis_probes.clear();
            for (unsigned i = 0; i < settings.max_initial_basis_probes; ++i) {
                auto probe = probe_section_refinement(topology, materials, origin, selected_spacing,
                    settings.mechanics, report, check);
                const auto state = probe.state;
                const auto recommendation = probe.recommended_spacing_m;
                result.last_initial_basis_probes.push_back(std::move(probe));
                if (state == SectionDirectionalProbe::State::ResourceLimited)
                    throw ResourceLimit("Local section study grid probe budget exhausted.");
                if (!recommendation) break;
                selected_spacing = *recommendation;
            }
            auto basis = refine_deposited_sections(topology, materials, origin, selected_spacing,
                settings.max_basis_levels, settings.relative_tolerance, settings.mechanics, report, check, samples);
            result.last_basis_history = std::move(basis.history);
            if (basis.state == SectionBasisStudy::State::ResourceLimited)
                throw ResourceLimit("Local section study mechanical budget exhausted.");
            if (basis.state != SectionBasisStudy::State::BasisConverged || !basis.field) return result;
            auto& field = *basis.field;
            double force = 0.;
            for (const auto& reaction : field.generalized_bed_reaction_n) force += std::max(0., -reaction.z());
            if (!result.history.empty()) {
                const auto& previous = result.history.back();
                const double volume_scale = std::max(previous.volume_m3, field.occupied_volume_m3);
                consecutive = (section_count > previous.sections || (quantization > 0. && quantization < previous.z_quantization_m)) &&
                    previous.components == topology.components() &&
                    stable(previous.volume_m3, field.occupied_volume_m3, volume_scale) &&
                    stable(previous.elastic_energy_j, field.elastic_energy_j,
                        std::max(maximum_energy_density * volume_scale,field.prescribed_strain_energy_j)) &&
                    stable(previous.tensile_generalized_bed_force_n, force,
                        std::max(maximum_eigenstress,field.prescribed_stress_bound_pa) * field.bed_area_m2) ?
                    consecutive + 1 : 0;
            }
            result.history.push_back({height, field.occupied_volume_m3, field.elastic_energy_j, force, section_count,
                topology.components(), basis.state, quantization, maximum_shift, collapsed, expanded_zero,
                double(geometry.maximum_boundary_uncertainty_nm) * 1e-9,
                geometry.added_representation_volume_bound_m3, geometry.removed_representation_volume_bound_m3});
            if (result.piecewise_constant_geometry || consecutive >= 2) {
                result.state = settings.initial_z_quantization_m || settings.sections.geometry.ownership_boundary_grid_units ?
                    DepositionSectionStudy::State::ApproximateGeometryAndBasisConverged :
                    DepositionSectionStudy::State::NominalGeometryAndBasisConverged;
                result.field = std::move(field);
                return result;
            }
        }
    } catch (const ResourceLimit& error) {
        if (callback_failed) throw;
        result.state = DepositionSectionStudy::State::ResourceLimited;
        result.failure_reason = std::string(error.what()).substr(0, 512);
        result.last_geometry_tiles = current_geometry.indexed_tiles;
        result.last_geometry_slabs = current_geometry.retained_slabs;
        result.section_count = previous_sections + current_geometry.planned_sections;
        result.boundary_points = previous_points + current_geometry.retained_points;
    }
    return result;
}
} // namespace Slic3r::Thermoelastic
