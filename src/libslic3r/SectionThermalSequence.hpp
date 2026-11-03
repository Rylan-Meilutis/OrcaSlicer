#pragma once
#include "SectionGrowth.hpp"
#include "SectionGrowthStream.hpp"
#include "ThermoelasticLayers.hpp"
#include "SectionThermalAssembly.hpp"
#include "ThermalDeposition.hpp"
#include "ThermalHistory.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

namespace Slic3r::Thermoelastic {
// Bounded retained-snapshot adapter. Caller supplies qualified prefix geometry,
// physical activation times, deposited-polymer temperatures and an explicit
// exposed-surface law. This does not infer any of them from final geometry.
class SectionThermalSequence
{
public:
    using Exposure = std::function<std::vector<ThermalHistory::Element>(const SectionThermalAssembly&)>;
    SectionThermalSequence(std::vector<SectionThermalMaterial> materials,
                           Eigen::Vector3d origin_m,
                           Eigen::Vector3d spacing_m,
                           SectionThermalAssemblyLimits assembly_limits = {},
                           SectionGrowthLimits growth_limits            = {},
                           std::size_t max_snapshots                    = 1000);
    ThermalHistory::DepositionEvent append(const SectionTopology& geometry,
                                           double time_s,
                                           const std::vector<double>& deposited_temperatures_k,
                                           const Exposure& exposure,
                                           const std::function<void()>& checkpoint = {});
    // Latest full geometry with IDs matching the emitted history. Its initial
    // temperatures are NOT the evolved field. Use only a successful history's
    // temperatures for mechanics; discard provisional results after failure.
    const SectionThermalAssembly& assembly() const;
    // Replace only surface operators on the current node numbering; no mass,
    // geometry reconstruction or snapshot-budget consumption.
    ThermalHistory::DepositionEvent replace_exposure(double time_s, const Exposure& exposure,
        const std::function<void()>& checkpoint = {});

private:
    std::vector<SectionThermalMaterial> m_materials;
    Eigen::Vector3d m_origin, m_spacing;
    SectionThermalAssemblyLimits m_assembly_limits;
    SectionGrowthLimits m_growth_limits;
    std::size_t m_max_snapshots, m_snapshots = 0;
    double m_time = 0.;
    std::size_t m_exposure_changes = 0;
    std::optional<SectionTopology> m_geometry;
    SectionThermalAssembly m_assembly;
};

// Same conservative event/state semantics with streamed polygon geometry.
// Retains path inputs and the assembled field, not complete polygon snapshots.
// Ownership providers must capture immutable data valid for the sequence's
// lifetime: they are replayed when certifying the next material increment.
class SectionThermalStreamSequence
{
public:
    using Exposure  = SectionThermalSequence::Exposure;
    using Ownership = std::function<std::optional<SectionOwnership>(double)>;
    SectionThermalStreamSequence(std::vector<SectionThermalMaterial> materials,
                                 ThermalSectionStreamSettings settings,
                                 SectionGrowthStreamLimits growth_limits = {},
                                 std::size_t max_snapshots               = 1000,
                                 std::vector<DepositedSegment> fixed_reference = {},
                                 Ownership fixed_reference_ownership = {});
    ThermalHistory::DepositionEvent append(const std::vector<DepositedSegment>& segments,
                                           double time_s,
                                           const std::vector<double>& deposited_temperatures_k,
                                           const Exposure& exposure,
                                           const Ownership& ownership                               = {},
                                           const std::vector<DepositedSegment>* partition_envelopes = nullptr,
                                           const std::function<void()>& checkpoint                  = {});
    const SectionThermalAssembly& assembly() const;
    ThermalHistory::DepositionEvent replace_exposure(double time_s, const Exposure& exposure,
        const std::function<void()>& checkpoint = {});
    double excluded_nominal_volume_m3() const;

private:
    std::vector<SectionThermalMaterial> m_materials;
    ThermalSectionStreamSettings m_settings;
    SectionGrowthStreamLimits m_growth_limits;
    std::size_t m_max_snapshots, m_snapshots = 0;
    double m_time = 0.;
    std::size_t m_exposure_changes = 0;
    std::vector<DepositedSegment> m_segments;
    std::optional<std::vector<DepositedSegment>> m_partitions;
    Ownership m_ownership;
    std::vector<DepositedSegment> m_reference;
    Ownership m_reference_ownership;
    double m_excluded_volume = 0.;
    SectionThermalAssembly m_assembly;
};
} // namespace Slic3r::Thermoelastic
