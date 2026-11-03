#include "ThermoelasticGCodeSections.hpp"
#include "ThermoelasticConflictOwnership.hpp"
#include "BoundingBox.hpp"
#include "ClipperUtils.hpp"
#include "ExPolygon.hpp"
#include "ExtrusionEntity.hpp"
#include "GCode/GCodeProcessor.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticSectionStudy.hpp"
#include "ThermoelasticSections.hpp"
#include "ThermoelasticToolpath.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <optional>
#include <stdexcept>
#include <vector>

namespace Slic3r::Thermoelastic {
SectionOverlapDiagnostic diagnose_section_overlap(const TimedToolpath& path,
    const LayerMaterialOverlap& overlap, BeadSection section, double bed_z_m,
    std::size_t max_checks, std::size_t max_reports, const std::function<void()>& checkpoint)
{
    if (path.state != ToolpathState::Complete || path.deposition.size() > 2000000 ||
        max_checks == 0 || max_checks > 4096 || max_reports == 0 || max_reports > 256 ||
        !std::isfinite(bed_z_m) || !std::isfinite(overlap.bottom_m) ||
        overlap.bottom_m != overlap.top_m ||
        (section != BeadSection::RectangularEnvelope && section != BeadSection::RoundedRectangle))
        throw std::invalid_argument("Invalid overlap diagnostic inputs.");
    SectionOverlapDiagnostic result;
    if (overlap.section_footprint.empty()) return result;
    std::size_t points = 0;
    for (const auto &polygon : overlap.section_footprint) {
        points += polygon.contour.points.size();
        for (const auto &hole : polygon.holes) points += hole.points.size();
    }
    if (points > 4096) throw ResourceLimit("Overlap diagnostic mask budget exceeded.");
    const auto box = get_extents(overlap.section_footprint);
    const double z = overlap.bottom_m;
    std::map<std::size_t, std::size_t> reports_by_filament;
    bool truncated = false;
    for (const auto &entry : path.deposition) {
        if (checkpoint) checkpoint();
        ++result.scanned_segments;
        auto bead = entry.geometry;
        bead.section = section;
        const double end_height = bead.end_height_m.value_or(bead.height_m);
        if (!bead.start_m.allFinite() || !bead.end_m.allFinite() ||
            !std::isfinite(bead.width_m) || bead.width_m <= 0. ||
            !std::isfinite(bead.height_m) || bead.height_m < 0. ||
            !std::isfinite(end_height) || end_height < 0. || std::max(bead.height_m, end_height) <= 0.)
            throw std::invalid_argument("Invalid diagnostic bead geometry.");
        if (z < std::min(bead.start_m.z()-bead.height_m, bead.end_m.z()-end_height) ||
            z >= std::max(bead.start_m.z(), bead.end_m.z())) continue;
        const double radius = bead.width_m * .5;
        if (std::max(bead.start_m.x(), bead.end_m.x()) + radius < double(box.min.x()) * 1e-9 ||
            std::min(bead.start_m.x(), bead.end_m.x()) - radius > double(box.max.x()) * 1e-9 ||
            std::max(bead.start_m.y(), bead.end_m.y()) + radius < double(box.min.y()) * 1e-9 ||
            std::min(bead.start_m.y(), bead.end_m.y()) - radius > double(box.max.y()) * 1e-9) continue;
        // Avoid filling the diagnostic with one early filament's tessellated
        // curve before reaching the other side of the material contact.
        const auto found = reports_by_filament.find(bead.material);
        if (found != reports_by_filament.end() && found->second >= std::min<std::size_t>(4, max_reports)) {
            truncated = true;
            continue;
        }
        if (result.polygon_checks >= max_checks || result.contributors.size() >= max_reports) {
            result.state = SectionOverlapDiagnostic::State::BudgetExceeded;
            return result;
        }
        ++result.polygon_checks;
        LayerGeometryLimits limits;
        limits.max_segments = 1;
        limits.max_boundary_points = 4096;
        const auto cut = section_deposition({bead}, z, bed_z_m, limits, checkpoint);
        double area = 0.;
        for (const auto &region : cut.regions)
            for (const auto &polygon : intersection_ex(region.footprint.native(), overlap.section_footprint))
                area += polygon.area() * 1e-18;
        if (area > 0.) {
            result.contributors.push_back({entry.move_index, bead.material, entry.gcode_id, entry.role, area});
            ++reports_by_filament[bead.material];
        }
    }
    result.state = truncated ? SectionOverlapDiagnostic::State::BudgetExceeded :
                   overlap.section_footprint_complete ? SectionOverlapDiagnostic::State::Complete :
                                                       SectionOverlapDiagnostic::State::PartialMask;
    return result;
}

GCodeSectionResult study_gcode_sections(const GCodeProcessorResult &source,
    const GCodeSectionScenario &scenario, std::size_t max_moves, std::size_t max_segments,
    const std::function<std::optional<SectionOwnership>(double)> &ownership,
    const std::function<void(unsigned, unsigned)> &progress, const std::function<void()> &checkpoint,
    const std::vector<SectionSampleRequest>& samples)
{
    if (scenario.materials.empty() || scenario.materials.size() > 256 ||
        scenario.filament_to_material.empty() || scenario.filament_to_material.size() > 256 ||
        max_moves == 0 || max_moves > 2000000 || max_segments == 0 || max_segments > 2000000)
        throw std::invalid_argument("Invalid G-code section study inputs.");
    for (auto material : scenario.filament_to_material)
        if (material >= scenario.materials.size())
            throw std::invalid_argument("Missing G-code section material mapping.");
    auto path = timed_deposition(source, max_moves, max_segments, checkpoint);
    if (path.state == ToolpathState::BudgetExceeded) throw ResourceLimit("G-code section path budget exceeded.");
    if (path.state != ToolpathState::Complete)
        throw std::invalid_argument("G-code section study requires complete processor timing and geometry.");
    const double end = scenario.observation_time_s.value_or(path.total_s);
    if (!std::isfinite(end) || end < 0. || end > path.total_s)
        throw std::invalid_argument("G-code section observation is outside the modeled timeline.");
    const auto excluded = [](const auto &entry) {
        return entry.role == erWipeTower || entry.role == erCustom || entry.role == erNone;
    };
    GCodeSectionResult result;
    const auto initial_count = path.deposition.size();
    path.deposition.erase(std::remove_if(path.deposition.begin(), path.deposition.end(), excluded), path.deposition.end());
    result.excluded_segments = initial_count - path.deposition.size();
    path.unresolved_deposition.erase(std::remove_if(path.unresolved_deposition.begin(), path.unresolved_deposition.end(), excluded),
        path.unresolved_deposition.end());
    for (const auto &entry : path.deposition)
        if (entry.start_s < end && entry.role == erMixed)
            throw std::invalid_argument("Mixed extrusion roles have no resolved mechanical material.");
    auto nominal = resolve_nominal_deposition(source, path, end, scenario.section,
        scenario.stationary_extrusion_allowance, checkpoint);
    result.stationary_extrusion = nominal.stationary_extrusion;
    result.commanded_path_volume_m3 = nominal.commanded_path_volume_m3;
    result.nominal_section_path_volume_m3 = nominal.nominal_section_path_volume_m3;
    result.segments_without_commanded_volume = nominal.segments_without_commanded_volume;
    for (std::size_t i = 0; i < nominal.segments.size(); ++i) {
        auto& segment = nominal.segments[i];
        if (checkpoint) checkpoint();
        if (segment.material >= scenario.filament_to_material.size())
            throw std::invalid_argument("Missing processor filament in mechanical material mapping.");
        segment.material = scenario.filament_to_material[segment.material];
        segment.section_source_index = i;
    }
    SampledConflictOwnershipResolver resolver;
    std::function<std::optional<SectionOwnership>(double)> selected_ownership = ownership;
    if (scenario.use_native_ownership) {
        if (ownership) throw std::invalid_argument("Native evidence and custom ownership cannot be combined.");
        // The prepared resolver owns its compact original-prefix geometry.
        // Release this first timed copy before preparing it; no source G-code
        // file or complete processor result is copied.
        path = {};
        ConflictOwnershipLimits limits;
        limits.max_moves = max_moves;
        limits.max_segments = max_segments;
        limits.geometry = scenario.study.sections.geometry;
        prepare_native_model_ownership(source, scenario.section, scenario.bed_z_m, scenario.filament_to_material,
            end, scenario.stationary_extrusion_allowance, limits, checkpoint, false, &resolver);
        selected_ownership = [&](double z) -> std::optional<SectionOwnership> {
            SectionOwnership request{z, {}};
            request.resolve_sampled_conflict = [&](const SampledSectionConflict& conflict) { return resolver(conflict).regions; };
            return request;
        };
    }
    result.study = study_deposition_sections(nominal.segments, scenario.materials,
        scenario.origin_m, scenario.initial_basis_spacing_m, scenario.bed_z_m,
        scenario.study, selected_ownership, progress, checkpoint, samples);
    return result;
}
} // namespace Slic3r::Thermoelastic
