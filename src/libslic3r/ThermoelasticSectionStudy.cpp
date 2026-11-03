#include "ThermoelasticSectionStudy.hpp"
#include "ThermoelasticLocalSections.hpp"
#include "ExPolygon.hpp"
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticMaterial.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticSections.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include "libslic3r.h"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iterator>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
namespace {
coord_t native_z(double value)
{
    if (!std::isfinite(value) || std::abs(value) > 10.)
        throw std::invalid_argument("Invalid deposition study height.");
    return std::llround(value * 1e9);
}
} // namespace
DepositionSectionStudy study_deposition_sections(const std::vector<DepositedSegment>& segments,
                                                 const std::vector<SectionMaterial>& materials,
                                                 const Eigen::Vector3d& origin,
                                                 const Eigen::Vector3d& basis_spacing,
                                                 double bed_z,
                                                 const SectionStudySettings& settings,
                                                 const std::function<std::optional<SectionOwnership>(double)>& ownership,
                                                 const std::function<void(unsigned, unsigned)>& progress,
                                                 const std::function<void()>& checkpoint,
                                                 const std::vector<SectionSampleRequest>& samples)
{
    if (settings.local_tile_size_m)
        return study_local_deposition_sections(segments, materials, origin, basis_spacing, *settings.local_tile_size_m,
            bed_z, settings, ownership, progress, checkpoint, samples);
    if (settings.initial_z_quantization_m || settings.sections.geometry.ownership_boundary_grid_units)
        throw std::invalid_argument("Explicit Z geometry quantization requires local section quadrature.");
    if (samples.size() > 200000 || segments.empty() || segments.size() > 2000000 || materials.empty() || materials.size() > 256 || settings.max_geometry_levels == 0 ||
        settings.max_geometry_levels > 4 || settings.max_basis_levels == 0 || settings.max_basis_levels > 5 ||
        settings.max_initial_basis_probes > 3 ||
        !std::isfinite(settings.relative_tolerance) || settings.relative_tolerance <= 0. || settings.relative_tolerance > 1. ||
        !settings.max_total_sections || settings.max_total_sections > 40000 || !settings.max_total_segment_visits ||
        settings.max_total_segment_visits > 80000000 || !settings.max_total_boundary_points ||
        settings.max_total_boundary_points > 8000000 || !settings.max_total_work || settings.max_total_work > 1000000000)
        throw std::invalid_argument("Invalid deposition section study or total budgets.");
    const coord_t bed = native_z(bed_z), maximum_height = native_z(settings.max_section_height_m);
    if (maximum_height < 2 || native_z(origin.z()) != bed)
        throw std::invalid_argument("Invalid deposition study section spacing or bed plane.");
    bool callback_failed = false;
    const auto check     = [&] {
        try {
            if (checkpoint)
                checkpoint();
        } catch (...) {
            callback_failed = true;
            throw;
        }
    };
    DepositionSectionStudy study;
    std::set<coord_t> boundaries;
    bool exact                    = !ownership;
    double maximum_energy_density = 0., maximum_eigenstress = 0.;
    for (const auto& material : materials) {
        if (!material.eigenstrain.allFinite() || !material.eigenstrain_gradient_per_m.allFinite())
            throw std::invalid_argument("Invalid study eigenstrain.");
        const auto stress      = (validated_elasticity(material.elasticity_pa) * material.eigenstrain).eval();
        maximum_energy_density = std::max(maximum_energy_density, .5 * material.eigenstrain.dot(stress));
        maximum_eigenstress    = std::max(maximum_eigenstress, stress.norm());
    }
    for (const auto& segment : segments) {
        check();
        const double end_height = segment.end_height_m.value_or(segment.height_m);
        if (segment.material >= materials.size() || !segment.start_m.allFinite() || !segment.end_m.allFinite() ||
            !std::isfinite(segment.width_m) || segment.width_m <= 0. || !std::isfinite(segment.height_m) || segment.height_m < 0. ||
            !std::isfinite(end_height) || end_height < 0.)
            throw std::invalid_argument("Invalid or unmapped deposition study bead.");
        const auto a = native_z(segment.start_m.z()), b = native_z(segment.end_m.z()), h = native_z(segment.height_m);
        const auto eh = native_z(end_height);
        if (std::max(h, eh) < 2 || std::min(a-h, b-eh) < bed)
            throw std::invalid_argument("Unresolved deposition study bead or bed intersection.");
        for (coord_t z : {a - h, b - eh, a, b})
            boundaries.insert(z);
        exact = exact && a == b && h == eh && segment.section == BeadSection::RectangularEnvelope;
        if (boundaries.size() > settings.max_total_sections + 1) {
            study.state = DepositionSectionStudy::State::ResourceLimited;
            study.failure_reason = "Distinct nominal bead Z boundaries exceed the section budget.";
            return study;
        }
    }
    study.piecewise_constant_geometry = exact;
    const auto charge                 = [](std::size_t& used, std::size_t amount, std::size_t maximum) {
        if (amount > maximum - used)
            throw ResourceLimit("Deposition study total budget exceeded.");
        used += amount;
    };
    const auto stable = [&](double a, double b, double scale) {
        return std::abs(a - b) <=
               settings.relative_tolerance * std::max(std::abs(a), std::abs(b)) + 128. * std::numeric_limits<double>::epsilon() * scale;
    };
    unsigned stable_geometry = 0;
    try {
        for (unsigned level = 0; level < (exact ? 1u : settings.max_geometry_levels); ++level) {
            check();
            const coord_t height = maximum_height / (coord_t(1) << level);
            if (height < 2)
                throw ResourceLimit("Deposition section refinement exceeds native height resolution.");
            std::vector<coord_t> cuts;
            for (auto lower = boundaries.begin(), upper = std::next(lower); upper != boundaries.end(); ++lower, ++upper) {
                const coord_t length = *upper - *lower, count = (length + height - 1) / height;
                if (count > length / 2)
                    throw ResourceLimit("A deposition section has no interior native sample plane.");
                if (std::size_t(count) > settings.max_total_sections - study.section_count - cuts.size())
                    throw ResourceLimit("Deposition study total section budget exceeded.");
                for (coord_t i = 0; i < count; ++i)
                    cuts.push_back(*lower + length * i / count);
            }
            cuts.push_back(*boundaries.rbegin());
            std::vector<double> planes;
            planes.reserve(cuts.size() - 1);
            for (std::size_t i = 1; i < cuts.size(); ++i)
                planes.push_back(double(cuts[i - 1] + (cuts[i] - cuts[i - 1]) / 2) * 1e-9);
            charge(study.section_count, planes.size(), settings.max_total_sections);
            // Reserve primitive maxima before entering non-interruptible work.
            charge(study.reserved_work, settings.topology.max_component_pairs, settings.max_total_work);
            std::vector<DepositedSlab> slabs;
            auto sweep               = settings.sections;
            sweep.max_segment_visits = std::min(sweep.max_segment_visits, settings.max_total_segment_visits - study.segment_visits);
            if (sweep.max_segment_visits == 0)
                throw ResourceLimit("Deposition study bead visit budget exhausted.");
            charge(study.reserved_work, sweep.max_segment_visits, settings.max_total_work);
            const auto visits = visit_deposition_sections(
                segments, planes, bed_z,
                [&](std::size_t i, const DepositedSection& section) {
                    // The provider reserves one plane before its polygon work.
                    // Keep consumed charges, releasing unused capacity for the
                    // next plane rather than reserving every plane's maximum.
                    study.reserved_work -= sweep.geometry.max_slab_group_visits - section.ownership_work;
                    charge(study.boundary_points, section.boundary_points, settings.max_total_boundary_points);
                    for (const auto& region : section.regions)
                        if (!region.footprint.native().empty()) {
                            if (slabs.size() >= settings.topology.max_slabs)
                                throw ResourceLimit("Deposition study slab budget exceeded.");
                            slabs.push_back({double(cuts[i]) * 1e-9, double(cuts[i + 1]) * 1e-9, region.material, region.footprint.native()});
                        }
                },
                sweep, check,
                [&](std::size_t, double z) -> std::optional<SectionOwnership> {
                    charge(study.reserved_work, sweep.geometry.max_slab_group_visits, settings.max_total_work);
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
                    } catch (...) {
                        callback_failed = true;
                        throw;
                    }
                });
            charge(study.segment_visits, visits, settings.max_total_segment_visits);
            const auto topology = qualify_section_topology(slabs, bed_z, settings.topology, check);
            // The factory owns its normalized copy; release sampled polygons
            // before allocating any mechanical matrices or output fields.
            std::vector<DepositedSlab>().swap(slabs);
            unsigned basis_run = 0;
            const auto basis_progress = [&](unsigned) {
                charge(study.reserved_work, settings.mechanics.max_work, settings.max_total_work);
                try {
                    if (progress)
                        progress(level, basis_run);
                    ++basis_run;
                } catch (...) {
                    callback_failed = true;
                    throw;
                }
            };
            Eigen::Vector3d selected_spacing = basis_spacing;
            study.last_initial_basis_probes.clear();
            for (unsigned i = 0; i < settings.max_initial_basis_probes; ++i) {
                auto probe = probe_section_refinement(topology, materials, origin, selected_spacing,
                                                       settings.mechanics, basis_progress, check);
                const auto state = probe.state;
                const auto recommendation = probe.recommended_spacing_m;
                study.last_initial_basis_probes.push_back(std::move(probe));
                if (state == SectionDirectionalProbe::State::ResourceLimited)
                    throw ResourceLimit("Deposition study initial-grid probe budget exhausted.");
                if (!recommendation) break;
                selected_spacing = *recommendation;
            }
            auto basis = refine_deposited_sections(topology, materials, origin, selected_spacing, settings.max_basis_levels,
                                                   settings.relative_tolerance, settings.mechanics, basis_progress, check, samples);
            study.last_basis_history = std::move(basis.history);
            if (basis.state == SectionBasisStudy::State::ResourceLimited)
                throw ResourceLimit("Deposition study basis budget exhausted.");
            if (basis.state != SectionBasisStudy::State::BasisConverged || !basis.field)
                return study;
            auto& field  = *basis.field;
            double force = 0.;
            for (const auto& reaction : field.generalized_bed_reaction_n)
                force += std::max(0., -reaction.z());
            if (!study.history.empty()) {
                const auto& previous      = study.history.back();
                const double volume_scale = std::max(previous.volume_m3, field.occupied_volume_m3);
                const double force_scale  = std::max(maximum_eigenstress,field.prescribed_stress_bound_pa) * field.bed_area_m2;
                stable_geometry = planes.size() > previous.sections && previous.components == topology.components() &&
                                          stable(previous.volume_m3, field.occupied_volume_m3, volume_scale) &&
                                          stable(previous.elastic_energy_j, field.elastic_energy_j,
                                              std::max(maximum_energy_density * volume_scale,field.prescribed_strain_energy_j)) &&
                                          stable(previous.tensile_generalized_bed_force_n, force, force_scale) ?
                                      stable_geometry + 1 :
                                      0;
            }
            study.history.push_back({double(height) * 1e-9, field.occupied_volume_m3, field.elastic_energy_j, force, planes.size(),
                                     topology.components(), basis.state});
            if (exact || stable_geometry >= 2) {
                study.field = std::move(field);
                study.state = DepositionSectionStudy::State::NominalGeometryAndBasisConverged;
                return study;
            }
        }
    } catch (const ResourceLimit& error) {
        if (callback_failed)
            throw;
        study.state = DepositionSectionStudy::State::ResourceLimited;
        study.failure_reason = std::string(error.what()).substr(0, 512);
    }
    return study;
}
} // namespace Slic3r::Thermoelastic
