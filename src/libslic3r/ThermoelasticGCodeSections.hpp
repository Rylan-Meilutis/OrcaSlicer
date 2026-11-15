#pragma once
#include "ThermoelasticSectionStudy.hpp"
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticSections.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticToolpath.hpp"
#include "ExtrusionEntity.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

namespace Slic3r {
struct GCodeProcessorResult;
namespace Thermoelastic {
struct GCodeSectionScenario {
    explicit GCodeSectionScenario(BeadSection chosen_section) : section(chosen_section) {}
    std::vector<SectionMaterial> materials;
    // Explicit processor-ID mapping, not an inferred physical tool mapping or
    // an assumption that same-named filaments have identical properties.
    std::vector<std::size_t> filament_to_material;
    BeadSection section;
    Eigen::Vector3d origin_m = Eigen::Vector3d::Zero();
    Eigen::Vector3d initial_basis_spacing_m = Eigen::Vector3d::Zero();
    double bed_z_m = 0.;
    SectionStudySettings study;
    std::optional<double> observation_time_s;
    std::optional<StationaryExtrusionAllowance> stationary_extrusion_allowance;
    // Live generation evidence only. Imported G-code without matching source
    // registry/snapshot remains unsupported, never assigned synthetic masks.
    bool use_native_ownership = false;
};
struct GCodeSectionResult {
    DepositionSectionStudy study;
    StationaryExtrusionSummary stationary_extrusion;
    double commanded_path_volume_m3 = 0., nominal_section_path_volume_m3 = 0.;
    std::size_t segments_without_commanded_volume = 0, excluded_segments = 0;
};
struct SectionOverlapContributor {
    std::size_t move_index, filament;
    unsigned gcode_id;
    ExtrusionRole role;
    double overlap_area_m2;
};
struct SectionOverlapDiagnostic {
    enum class State { Complete, BudgetExceeded, MaskUnavailable, PartialMask };
    State state = State::MaskUnavailable;
    std::vector<SectionOverlapContributor> contributors;
    std::size_t scanned_segments = 0, polygon_checks = 0;
};
// Read-only evidence for a rejected section, not a resolution of its material
// interface. Input is the same timed source used by the geometry audit; retain
// final filament IDs. Explicit shape and bed plane must match that audit.
// Reports at most four contributors per filament so early curves do not hide
// the opposite side of the contact; this truncation is explicitly labeled.
// Incomplete diagnostics are labeled; callbacks propagate. No G-code is written.
SectionOverlapDiagnostic diagnose_section_overlap(const TimedToolpath& path,
    const LayerMaterialOverlap& overlap, BeadSection section, double bed_z_m,
    std::size_t max_polygon_checks = 256, std::size_t max_reports = 32,
    const std::function<void()>& checkpoint = {});
// Cancellable read-only adapter. No implicit material or thermal properties,
// file writes, G-code mutation, or preview publication. Uses actual processor
// filament IDs and sloped endpoints, excluding tower/custom/unclassified paths.
// Observation timing is modeled processor timing, not measured deposition time.
// The chosen nominal section is explicit; parsed display dimensions alone do
// not establish real bead shape. Unsupported geometry and incomplete input
// throw; resource-limited studies retain no mechanical field.
GCodeSectionResult study_gcode_sections(const GCodeProcessorResult &source,
    const GCodeSectionScenario &scenario, std::size_t max_moves = 2000000,
    std::size_t max_segments = 250000,
    const std::function<std::optional<SectionOwnership>(double)> &ownership = {},
    const std::function<void(unsigned, unsigned)> &progress = {},
    const std::function<void()> &checkpoint = {},
    const std::vector<SectionSampleRequest>& samples = {});
} // namespace Thermoelastic
} // namespace Slic3r
