#include "SectionThermalSequence.hpp"
#include "SectionGrowth.hpp"
#include "SectionGrowthStream.hpp"
#include "ThermoelasticLayers.hpp"
#include "SectionThermalAssembly.hpp"
#include "ThermalDeposition.hpp"
#include "ThermalHistory.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include <Eigen/Core>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
namespace {
using NodeKey = std::pair<std::array<long long, 3>, std::size_t>;
NodeKey key(const SectionThermalNode& node, std::size_t component)
{
    return {{{std::llround(node.position_m.x() * 1e9), std::llround(node.position_m.y() * 1e9), std::llround(node.position_m.z() * 1e9)}},
            component};
}
void remap(std::array<std::size_t, 8>& nodes, const std::vector<std::size_t>& indices)
{
    for (auto& node : nodes)
        node = indices.at(node);
}
} // namespace

SectionThermalSequence::SectionThermalSequence(std::vector<SectionThermalMaterial> materials,
                                               Eigen::Vector3d origin,
                                               Eigen::Vector3d spacing,
                                               SectionThermalAssemblyLimits assembly_limits,
                                               SectionGrowthLimits growth_limits,
                                               std::size_t max_snapshots)
    : m_materials(std::move(materials))
    , m_origin(origin)
    , m_spacing(spacing)
    , m_assembly_limits(assembly_limits)
    , m_growth_limits(growth_limits)
    , m_max_snapshots(max_snapshots)
{
    if (!max_snapshots || max_snapshots > 10000)
        throw std::invalid_argument("Invalid thermal snapshot limit.");
}

const SectionThermalAssembly& SectionThermalSequence::assembly() const
{
    if (!m_geometry)
        throw std::logic_error("No thermal geometry has been activated.");
    return m_assembly;
}

static void validate_exposure(const SectionThermalAssembly& next,
    const std::vector<ThermalHistory::Element>& boundaries, std::size_t boundary_limit,
    const std::function<void()>& checkpoint)
{
    if (!next.exposure_complete)
        throw std::invalid_argument("Incomplete section thermal exposure.");
    if (boundaries.size() > boundary_limit)
        throw ResourceLimit("Section exposure operator budget exceeded.");
    for (const auto& boundary : boundaries) {
        if (checkpoint) checkpoint();
        if (!boundary.boundary || !ThermalHistory::valid_element(boundary, next.nodes.size()))
            throw std::invalid_argument("Invalid section exposure operator.");
    }
}

static ThermalHistory::DepositionEvent exposure_event(const SectionThermalAssembly& assembly,
    double time, double previous, std::size_t count, const SectionThermalSequence::Exposure& exposure,
    std::size_t boundary_limit, const std::function<void()>& checkpoint)
{
    if (!std::isfinite(time) || time < previous || !exposure)
        throw std::invalid_argument("Invalid section exposure change time or provider.");
    if (count >= 10000)
        throw ResourceLimit("Section exposure change budget exceeded.");
    if (checkpoint) checkpoint();
    ThermalHistory::DepositionEvent event;
    event.time_s = time;
    event.exposed_elements = exposure(assembly);
    validate_exposure(assembly, event.exposed_elements, boundary_limit, checkpoint);
    return event;
}

ThermalHistory::DepositionEvent SectionThermalSequence::replace_exposure(double time,
    const Exposure& exposure, const std::function<void()>& checkpoint)
{
    auto event = exposure_event(assembly(), time, m_time, m_exposure_changes, exposure, m_assembly_limits.max_groups, checkpoint);
    m_time = time;
    ++m_exposure_changes;
    return event;
}

ThermalHistory::DepositionEvent SectionThermalStreamSequence::replace_exposure(double time,
    const Exposure& exposure, const std::function<void()>& checkpoint)
{
    auto event = exposure_event(assembly(), time, m_time, m_exposure_changes, exposure, m_settings.assembly.max_groups, checkpoint);
    m_time = time;
    ++m_exposure_changes;
    return event;
}

static ThermalHistory::DepositionEvent make_event(const SectionThermalAssembly* previous,
                                                  SectionThermalAssembly& next,
                                                  const SectionThermalAssembly& added,
                                                  const std::vector<std::size_t>& components,
                                                  double time,
                                                  std::vector<ThermalHistory::Element> boundaries,
                                                  std::size_t boundary_limit,
                                                  const std::function<void()>& checkpoint)
{
    validate_exposure(next, boundaries, boundary_limit, checkpoint);
    ThermalHistory::DepositionEvent event;
    if (!previous) {
        event = next.initial;
    } else {
        std::map<NodeKey, std::size_t> next_nodes;
        for (std::size_t i = 0; i < next.nodes.size(); ++i) {
            if (checkpoint)
                checkpoint();
            if (!next_nodes.emplace(key(next.nodes[i], next.nodes[i].component), i).second)
                throw std::invalid_argument("Duplicate thermal geometry node.");
        }
        const auto missing = std::numeric_limits<std::size_t>::max();
        std::vector<std::size_t> order(next.nodes.size(), missing), old_map(previous->nodes.size());
        std::size_t existing = 0;
        for (std::size_t i = 0; i < previous->nodes.size(); ++i) {
            if (checkpoint)
                checkpoint();
            const auto& node = previous->nodes[i];
            const auto id    = next_nodes.at(key(node, components.at(node.component)));
            if (order[id] == missing)
                order[id] = existing++;
            old_map[i] = order[id];
        }
        std::size_t count = existing;
        for (auto& id : order)
            if (id == missing)
                id = count++;
        event.coalesced_nodes = std::move(old_map);
        std::vector<std::size_t> added_map(added.nodes.size());
        std::vector<ThermalHistory::Cell> new_cells(count - existing);
        std::vector<bool> assigned(count - existing, false);
        for (std::size_t i = 0; i < added.nodes.size(); ++i) {
            if (checkpoint)
                checkpoint();
            const auto& node = added.nodes[i];
            const auto id    = order[next_nodes.at(key(node, node.component))];
            added_map[i]     = id;
            const auto& cell = added.initial.deposited_cells[i];
            if (id < existing)
                event.added_mass.push_back({id, cell.volume_m3, cell.volumetric_heat_capacity_j_m3_k, cell.temperature_k});
            else {
                if (assigned[id - existing])
                    throw std::invalid_argument("Duplicate new thermal node.");
                assigned[id - existing]  = true;
                new_cells[id - existing] = cell;
            }
        }
        for (bool present : assigned)
            if (!present)
                throw std::invalid_argument("New thermal node has no deposited mass.");
        event.deposited_cells = std::move(new_cells);
        event.added_elements  = added.initial.added_elements;
        for (auto& element : event.added_elements)
            remap(element.nodes, added_map);
        for (auto& boundary : boundaries)
            remap(boundary.nodes, order);
        auto nodes = next.nodes;
        auto cells = next.initial.deposited_cells;
        for (std::size_t i = 0; i < order.size(); ++i) {
            next.nodes[order[i]]                   = nodes[i];
            next.initial.deposited_cells[order[i]] = cells[i];
        }
        for (auto& group : next.groups)
            remap(group.nodes, order);
        for (auto& face : next.horizontal_faces)
            remap(face.nodes, order);
        for (auto& face : next.vertical_faces)
            remap(face.nodes, order);
        for (auto& element : next.initial.added_elements)
            remap(element.nodes, order);
    }
    event.time_s           = time;
    event.exposed_elements = std::move(boundaries);
    return event;
}

ThermalHistory::DepositionEvent SectionThermalSequence::append(const SectionTopology& geometry,
                                                               double time,
                                                               const std::vector<double>& temperatures,
                                                               const Exposure& exposure,
                                                               const std::function<void()>& checkpoint)
{
    if (!std::isfinite(time) || time < m_time || !exposure || temperatures.size() != m_materials.size())
        throw std::invalid_argument("Invalid section activation time, temperatures or exposure provider.");
    if (m_snapshots >= m_max_snapshots)
        throw ResourceLimit("Section thermal snapshot budget exceeded.");
    auto materials = m_materials;
    for (std::size_t i = 0; i < temperatures.size(); ++i)
        materials[i].initial_temperature_k = temperatures[i];
    auto next = assemble_section_thermal_volume(geometry, materials, m_origin, m_spacing, m_assembly_limits, checkpoint);
    SectionThermalAssembly added;
    std::vector<std::size_t> components;
    if (m_geometry) {
        const auto growth = qualify_section_growth(*m_geometry, geometry, m_growth_limits, checkpoint);
        added             = assemble_section_thermal_addition(growth, materials, m_origin, m_spacing, m_assembly_limits, checkpoint);
        components        = growth.previous_components();
    }
    auto event = make_event(m_geometry ? &m_assembly : nullptr, next, added, components, time, exposure(next), m_assembly_limits.max_groups,
                            checkpoint);
    // Commit only after all geometry, material, node and boundary checks pass.
    auto retained_geometry = geometry;
    m_geometry             = std::move(retained_geometry);
    m_assembly             = std::move(next);
    m_time                 = time;
    ++m_snapshots;
    return event;
}

SectionThermalStreamSequence::SectionThermalStreamSequence(std::vector<SectionThermalMaterial> materials,
                                                           ThermalSectionStreamSettings settings,
                                                           SectionGrowthStreamLimits growth_limits,
                                                           std::size_t max_snapshots,
                                                           std::vector<DepositedSegment> reference,Ownership reference_ownership)
    : m_materials(std::move(materials)), m_settings(std::move(settings)), m_growth_limits(growth_limits), m_max_snapshots(max_snapshots),
      m_reference(std::move(reference)),m_reference_ownership(std::move(reference_ownership))
{
    if (!max_snapshots || max_snapshots > 10000)
        throw std::invalid_argument("Invalid streamed thermal snapshot limit.");
    if (m_reference.empty() && m_reference_ownership)
        throw std::invalid_argument("Reference ownership requires fixed reference geometry.");
    m_settings.include_exposure = true;
}
const SectionThermalAssembly& SectionThermalStreamSequence::assembly() const
{
    if (!m_snapshots)
        throw std::logic_error("No streamed thermal geometry has been activated.");
    return m_assembly;
}
double SectionThermalStreamSequence::excluded_nominal_volume_m3() const
{
    assembly();
    return m_excluded_volume;
}
ThermalHistory::DepositionEvent SectionThermalStreamSequence::append(const std::vector<DepositedSegment>& segments,
                                                                     double time,
                                                                     const std::vector<double>& temperatures,
                                                                     const Exposure& exposure,
                                                                     const Ownership& ownership,
                                                                     const std::vector<DepositedSegment>* partitions,
                                                                     const std::function<void()>& checkpoint)
{
    if (!std::isfinite(time) || time < m_time || !exposure || temperatures.size() != m_materials.size())
        throw std::invalid_argument("Invalid streamed activation time, temperatures or exposure provider.");
    if (m_snapshots >= m_max_snapshots)
        throw ResourceLimit("Streamed thermal snapshot budget exceeded.");
    if (!m_reference.empty() && (ownership || partitions))
        throw std::invalid_argument("Fixed reference history cannot replace its ownership or sampling envelopes.");
    auto materials = m_materials;
    for (std::size_t i = 0; i < temperatures.size(); ++i)
        materials[i].initial_temperature_k = temperatures[i];
    SectionThermalGrowthAssembly grown;
    double excluded = 0.;
    if (!m_reference.empty()) {
        if (!m_snapshots) {
            auto initial=assemble_reference_thermal_volume(m_reference,segments,materials,m_settings,m_growth_limits,
                m_reference_ownership,checkpoint);
            grown.next=std::move(initial.assembly);
            excluded=initial.excluded_nominal_volume_m3;
        } else {
            grown=assemble_reference_thermal_growth(m_reference,m_segments,segments,materials,m_settings,m_growth_limits,
                m_reference_ownership,checkpoint);
            excluded=grown.next_excluded_nominal_volume_m3;
        }
    } else if (!m_snapshots)
        grown.next = assemble_local_deposition_thermal_volume(segments, materials, m_settings, ownership, checkpoint, partitions).assembly;
    else
        grown = assemble_local_thermal_growth(m_segments, segments, materials, m_settings, m_growth_limits, m_ownership, ownership,
                                              checkpoint, m_partitions ? &*m_partitions : nullptr, partitions);
    auto event = make_event(m_snapshots ? &m_assembly : nullptr, grown.next, grown.added, grown.previous_components, time,
                            exposure(grown.next), m_settings.assembly.max_groups, checkpoint);
    // Keep only path/index inputs between snapshots, never the polygon history.
    // Prepare all potentially throwing copies before committing the new state.
    auto saved_segments = segments;
    std::optional<std::vector<DepositedSegment>> saved_partitions;
    if (partitions)
        saved_partitions = *partitions;
    auto saved_ownership = ownership;
    m_segments           = std::move(saved_segments);
    m_partitions         = std::move(saved_partitions);
    m_ownership          = std::move(saved_ownership);
    m_assembly           = std::move(grown.next);
    m_time               = time;
    m_excluded_volume    = excluded;
    ++m_snapshots;
    return event;
}
} // namespace Slic3r::Thermoelastic
