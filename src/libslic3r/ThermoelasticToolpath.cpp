#include "ThermoelasticToolpath.hpp"
#include "ThermoelasticMesh.hpp"
#include "GCode/GCodeProcessor.hpp"
#include "GCode/DepositionSection.hpp"
#include "GCode/DepositionProvenance.hpp"
#include "ExtrusionEntity.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>

namespace Slic3r::Thermoelastic {
namespace {
bool authenticated_vertical_deposition(const GCodeProcessorResult& source,std::size_t index,std::size_t previous)
{
    if (previous>=index || index>=source.moves.size() || !source.deposition_provenance || !source.deposition_provenance->sealed() ||
        !source.deposition_sections || !source.deposition_sections->sealed() ||
        source.deposition_sections->generation()!=source.deposition_provenance->generation()) return false;
    const auto& move=source.moves[index];
    const auto& start=source.moves[previous];
    const auto* owner=source.deposition_provenance->source(move.deposition_source_id);
    const auto* parent=source.deposition_sections->section(move.deposition_command_section_id);
    if (!move.deposition_vertical_extrusion || move.type!=EMoveType::Travel || move.internal_only ||
        start.internal_only || start.type==EMoveType::Seam || !owner || owner->height!=DepositionSource::Height::Transformed ||
        owner->logical_filament!=move.extruder_id || !parent || parent->filament!=move.extruder_id ||
        parent->transform!=DepositionCommandSection::Transform::BrickCourse || parent->arc_direction || parent->nominal_lower_z ||
        parent->start_height!=parent->end_height || parent->start_height<=0. ||
        (move.extrusion_role!=erPerimeter && move.extrusion_role!=erExternalPerimeter && move.extrusion_role!=erStaggeredPerimeter) ||
        !std::isfinite(move.delta_extruder) || move.delta_extruder<=0. ||
        !move.deposition_design_heights || !move.deposition_design_width ||
        !std::isfinite(*move.deposition_design_width) || *move.deposition_design_width<=0. ||
        float(parent->width)!=*move.deposition_design_width ||
        parent->start[0]!=parent->end[0] || parent->start[1]!=parent->end[1] || parent->start[2]==parent->end[2] ||
        start.position.x()!=move.position.x() || start.position.y()!=move.position.y() || start.position.z()==move.position.z() ||
        double(start.position.x())!=parent->start[0] || double(start.position.y())!=parent->start[1]) return false;
    const std::array<double,3> from{start.position.x(),start.position.y(),start.position.z()};
    const std::array<double,3> to{move.position.x(),move.position.y(),move.position.z()};
    DepositionSectionChain chain(*parent);
    if (from!=parent->start && !chain.append(parent->start,from,move.extruder_id)) return false;
    const auto heights=chain.append(from,to,move.extruder_id);
    return heights && *heights==*move.deposition_design_heights &&
        std::isfinite((*heights)[0]) && std::isfinite((*heights)[1]) && (*heights)[0]>0. && (*heights)[1]>0.;
}
} // namespace

std::optional<std::array<double, 2>> authenticated_deposition_z(
    const GCodeProcessorResult &source, std::size_t index, std::size_t previous)
{
    const auto *span = source.deposition_z_for_move(index);
    const auto invalid = [] { throw std::invalid_argument("Precise deposition Z does not match its authenticated command."); };
    if (!span) {
        if (index < source.moves.size() && source.deposition_sections) {
            const auto *command = source.deposition_sections->section(source.moves[index].deposition_command_section_id);
            if (command && command->nominal_lower_z) invalid();
        }
        return {};
    }
    if (index >= source.moves.size() || previous >= index ||
        !source.deposition_sections || !source.deposition_sections->sealed() ||
        !source.deposition_provenance || !source.deposition_provenance->sealed() ||
        source.deposition_z_generation != source.deposition_sections->generation() ||
        source.deposition_z_generation != source.deposition_provenance->generation()) invalid();
    const auto &move = source.moves[index], &start = source.moves[previous];
    const auto *command = source.deposition_sections->section(move.deposition_command_section_id);
    const auto *owner = source.deposition_provenance->source(move.deposition_source_id);
    if (span->move_index != index || span->source_id != move.deposition_source_id ||
        span->command_id != move.deposition_command_section_id || !owner ||
        owner->logical_filament != move.extruder_id || !command || !command->nominal_lower_z ||
        command->filament != move.extruder_id || move.type != EMoveType::Extrude || move.internal_only ||
        start.internal_only || start.type == EMoveType::Seam ||
        !move.deposition_design_heights || !move.deposition_design_width ||
        !std::isfinite(span->z_mm[0]) || !std::isfinite(span->z_mm[1]) ||
        span->preview_z_mm[0] != start.position.z() || span->preview_z_mm[1] != move.position.z()) invalid();
    const std::array<double, 3> from{start.position.x(), start.position.y(), span->z_mm[0]};
    const std::array<double, 3> to{move.position.x(), move.position.y(), span->z_mm[1]};
    DepositionSectionChain chain(*command);
    if (from != command->start && !chain.append(command->start, from, move.extruder_id)) invalid();
    const auto heights = chain.append(from, to, move.extruder_id);
    if (!heights || *heights != *move.deposition_design_heights ||
        float(command->width) != *move.deposition_design_width) invalid();
    return span->z_mm;
}

NominalDeposition resolve_nominal_deposition(const GCodeProcessorResult &source,
    const TimedToolpath &path, double end, BeadSection section,
    const std::optional<StationaryExtrusionAllowance> &stationary_allowance,
    const std::function<void()> &checkpoint,
    const std::optional<CollapsedEnvelopeAllowance> &collapsed_allowance)
{
    if (path.state != ToolpathState::Complete || path.deposition.empty() ||
        path.deposition.size() > 2000000 || path.unresolved_deposition.size() > 2000000 ||
        !std::isfinite(end) || end < 0. ||
        (section != BeadSection::RectangularEnvelope && section != BeadSection::RoundedRectangle))
        throw std::invalid_argument("Invalid nominal deposition input.");
    if (stationary_allowance) {
        const auto &allow = *stationary_allowance;
        if (!std::isfinite(allow.maximum_volume_m3) || allow.maximum_volume_m3<=0. || allow.maximum_volume_m3>1e-12 ||
            !std::isfinite(allow.maximum_total_fraction) || allow.maximum_total_fraction<=0. || allow.maximum_total_fraction>1e-6 ||
            !std::isfinite(allow.maximum_preceding_command_fraction) || allow.maximum_preceding_command_fraction<=0. ||
            allow.maximum_preceding_command_fraction>1e-3)
            throw std::invalid_argument("Stationary extrusion approximation exceeds its admission limits.");
    }
    if (collapsed_allowance) {
        const auto &a = *collapsed_allowance;
        if (!std::isfinite(a.maximum_span_volume_m3) || a.maximum_span_volume_m3 <= 0. || a.maximum_span_volume_m3 > 1e-12 ||
            !std::isfinite(a.maximum_total_volume_m3) || a.maximum_total_volume_m3 <= 0. || a.maximum_total_volume_m3 > 1e-11 ||
            !std::isfinite(a.maximum_total_fraction) || a.maximum_total_fraction <= 0. || a.maximum_total_fraction > 1e-6 ||
            !std::isfinite(a.maximum_neighbor_fraction) || a.maximum_neighbor_fraction <= 0. || a.maximum_neighbor_fraction > .1)
            throw std::invalid_argument("Collapsed envelope omission exceeds its admission limits.");
        // A caller may supply a role-filtered timeline, but cannot move an
        // unresolved command or its proposed neighbor earlier in time.
        const auto ordered = [](const auto &entries) {
            for (std::size_t i = 1; i < entries.size(); ++i)
                if (entries[i - 1].move_index >= entries[i].move_index) return false;
            return true;
        };
        if (source.moves.size() > 2000000 || !ordered(path.deposition) || !ordered(path.unresolved_deposition))
            throw std::invalid_argument("Invalid collapsed envelope witness timeline.");
        double time = 0.;
        std::size_t bead = 0, unresolved = 0;
        for (std::size_t i = 0; i < source.moves.size(); ++i) {
            if (checkpoint && i % 1024 == 0) checkpoint();
            const double duration = source.moves[i].time[0], begin = time;
            if (!std::isfinite(duration) || duration < 0.) throw std::invalid_argument("Invalid collapsed envelope timing.");
            time += duration;
            if (bead < path.deposition.size() && path.deposition[bead].move_index == i) {
                const auto &entry = path.deposition[bead++];
                if (entry.start_s != begin || entry.end_s != time) throw std::invalid_argument("Invalid collapsed neighbor timing.");
            }
            if (unresolved < path.unresolved_deposition.size() && path.unresolved_deposition[unresolved].move_index == i) {
                const auto &entry = path.unresolved_deposition[unresolved++];
                if (entry.start_s != begin || entry.end_s != time) throw std::invalid_argument("Invalid collapsed envelope timing.");
            }
        }
        if (bead != path.deposition.size() || unresolved != path.unresolved_deposition.size())
            throw std::invalid_argument("Missing collapsed envelope timeline move.");
    }
    NominalDeposition result;
    std::map<std::size_t, double> omitted_at_neighbor;
    std::size_t witness_work = 0;
    const auto preceding_move = [&](std::size_t index) {
        if (index >= source.moves.size()) throw std::invalid_argument("Invalid collapsed envelope move.");
        while (index > 0) {
            if (++witness_work > 2000000) throw std::invalid_argument("Collapsed envelope witness budget exceeded.");
            --index;
            if (!source.moves[index].internal_only && source.moves[index].type != EMoveType::Seam) return index;
        }
        throw std::invalid_argument("Missing collapsed envelope origin.");
    };
    const auto commanded_volume = [&](std::size_t index, std::size_t material) {
        if (index >= source.moves.size() || material >= source.filament_diameters.size())
            throw std::invalid_argument("Missing collapsed envelope extrusion reference.");
        const double diameter = source.filament_diameters[material] * .001;
        const double e = source.moves[index].delta_extruder * .001;
        const double volume = std::acos(-1.) * .25 * diameter * diameter * e;
        if (!std::isfinite(diameter) || diameter <= 0. || !std::isfinite(e) || e <= 0. || !std::isfinite(volume) || volume <= 0.)
            throw std::invalid_argument("Invalid collapsed envelope extrusion reference.");
        return volume;
    };
    auto &segments = result.segments;
    // Process preceding commands once, including all tessellated arc pieces.
    // Never reconstruct an XY line or silently erase the stationary E amount.
    std::size_t preceding = 0;
    unsigned preceding_command = 0;
    double preceding_volume = 0., admitted_at_command = 0.;
    double previous_stationary_end = -1.;
    bool preceding_volume_known = true;
    for (const auto &point : path.unresolved_deposition) {
        if (point.start_s > end) break;
        if (checkpoint) checkpoint();
        if (point.reason == UnresolvedSpatialDeposition::Reason::CollapsedEnvelope) {
            if (!collapsed_allowance || point.end_s > end)
                throw std::invalid_argument("Unresolved collapsed envelope exists at the observation time.");
            const auto previous = preceding_move(point.move_index);
            const auto z = authenticated_deposition_z(source, point.move_index, previous);
            const auto &move = source.moves[point.move_index];
            Eigen::Vector3d from = source.moves[previous].position.cast<double>() * .001;
            Eigen::Vector3d to = move.position.cast<double>() * .001;
            if (!z || !move.deposition_design_heights || (*move.deposition_design_heights)[0] != 0. ||
                (*move.deposition_design_heights)[1] != 0. || point.source_id != move.deposition_source_id ||
                point.command_id != move.deposition_command_section_id || point.material != move.extruder_id ||
                point.extruded_filament_mm != double(move.delta_extruder))
                throw std::invalid_argument("Collapsed envelope does not match its authenticated extrusion.");
            from.z() = (*z)[0] * .001;
            to.z() = (*z)[1] * .001;
            if ((from - point.start_m).squaredNorm() != 0. || (to - point.end_m).squaredNorm() != 0.)
                throw std::invalid_argument("Collapsed envelope endpoints do not match its authenticated extrusion.");
            const double volume = commanded_volume(point.move_index, point.material);
            const auto &a = *collapsed_allowance;
            if (volume > a.maximum_span_volume_m3 || result.collapsed_envelope.commanded_volume_m3 + volume > a.maximum_total_volume_m3)
                throw std::invalid_argument("Collapsed envelope exceeds its absolute volume allowance.");
            const auto next = std::lower_bound(path.deposition.begin(), path.deposition.end(), point.move_index,
                [](const auto &bead, std::size_t index) { return bead.move_index < index; });
            const TimedDeposition *candidates[2] = {next == path.deposition.end() ? nullptr : &*next,
                next == path.deposition.begin() ? nullptr : &*(next - 1)};
            bool admitted = false;
            for (const auto *neighbor : candidates) {
                if (!neighbor || neighbor->end_s > end || neighbor->move_index >= source.moves.size()) continue;
                const auto &nm = source.moves[neighbor->move_index];
                const auto &g = neighbor->geometry;
                const bool touches = neighbor->move_index < point.move_index ?
                    neighbor->end_s == point.start_s && (g.end_m - point.start_m).squaredNorm() == 0. :
                    neighbor->start_s == point.end_s && (g.start_m - point.end_m).squaredNorm() == 0.;
                if (!touches || nm.deposition_source_id != point.source_id || g.material != point.material ||
                    nm.extruder_id != point.material || nm.extrusion_role != point.role ||
                    !nm.deposition_design_heights || !nm.deposition_design_width ||
                    g.height_m != (*nm.deposition_design_heights)[0] * .001 ||
                    g.end_height_m.value_or(g.height_m) != (*nm.deposition_design_heights)[1] * .001 ||
                    std::max(g.height_m, g.end_height_m.value_or(g.height_m)) <= 0. ||
                    g.width_m != double(*nm.deposition_design_width) * .001 ||
                    (point.end_m - point.start_m).norm() > g.width_m * .5) continue;
                const auto np = preceding_move(neighbor->move_index);
                const auto nz = authenticated_deposition_z(source, neighbor->move_index, np);
                if (!nz) continue;
                Eigen::Vector3d nf = source.moves[np].position.cast<double>() * .001, nt = nm.position.cast<double>() * .001;
                nf.z() = (*nz)[0] * .001;
                nt.z() = (*nz)[1] * .001;
                if ((nf - g.start_m).squaredNorm() != 0. || (nt - g.end_m).squaredNorm() != 0.) continue;
                const double reference = commanded_volume(neighbor->move_index, point.material);
                if (omitted_at_neighbor[neighbor->move_index] + volume > a.maximum_neighbor_fraction * reference) continue;
                omitted_at_neighbor[neighbor->move_index] += volume;
                admitted = true;
                break;
            }
            if (!admitted) throw std::invalid_argument("Collapsed envelope lacks a bounded active same-source neighbor.");
            ++result.collapsed_envelope.count;
            result.collapsed_envelope.commanded_volume_m3 += volume;
            continue;
        }
        if (point.reason != UnresolvedSpatialDeposition::Reason::UnresolvedXYFootprint ||
            !stationary_allowance || point.end_s>end || (point.end_m-point.start_m).squaredNorm()!=0.)
            throw std::invalid_argument("Unresolved spatial deposition exists at the observation time.");
        while (preceding<path.deposition.size() && path.deposition[preceding].end_s<=point.start_s) {
            const auto &entry = path.deposition[preceding++];
            if (entry.gcode_id != preceding_command) {
                preceding_command = entry.gcode_id;
                preceding_volume = 0.;
                admitted_at_command = 0.;
                preceding_volume_known = true;
            }
            if (entry.commanded_area_m2)
                preceding_volume += *entry.commanded_area_m2*(entry.geometry.end_m-entry.geometry.start_m).norm();
            else preceding_volume_known = false;
        }
        if (preceding==0 || point.material>=source.filament_diameters.size())
            throw std::invalid_argument("Stationary extrusion has no preceding bead or filament diameter.");
        const auto &neighbor = path.deposition[preceding-1];
        if (neighbor.geometry.material!=point.material || neighbor.role!=point.role ||
            (neighbor.geometry.end_m-point.start_m).squaredNorm()!=0. ||
            (neighbor.end_s!=point.start_s && previous_stationary_end!=point.start_s))
            throw std::invalid_argument("Stationary extrusion is not at the preceding same-material bead endpoint.");
        const double diameter = source.filament_diameters[point.material]*.001;
        const double volume = std::acos(-1.)*.25*diameter*diameter*point.extruded_filament_mm*.001;
        if (!std::isfinite(diameter) || diameter<=0. || !std::isfinite(volume) || volume<=0. ||
            !preceding_volume_known || !std::isfinite(preceding_volume) || preceding_volume<=0.)
            throw std::invalid_argument("Invalid stationary extrusion volume reference.");
        admitted_at_command += volume;
        result.stationary_extrusion.commanded_volume_m3 += volume;
        ++result.stationary_extrusion.count;
        previous_stationary_end = point.end_s;
        const auto &allow = *stationary_allowance;
        if (result.stationary_extrusion.commanded_volume_m3>allow.maximum_volume_m3 ||
            admitted_at_command>allow.maximum_preceding_command_fraction*preceding_volume)
            throw std::invalid_argument("Stationary extrusion exceeds its absolute or local volume allowance.");
    }
    auto &retained_fractions = result.retained_fractions;
    segments.reserve(path.deposition.size());
    for (const auto &entry : path.deposition) {
        if (checkpoint) checkpoint();
        if (entry.start_s >= end) break;
        const double fraction = std::min(1., (end - entry.start_s) / (entry.end_s - entry.start_s));
        auto geometry = entry.geometry;
        geometry.section = section;
        // A complete child already carries the authenticated endpoint. Even
        // start + 1*(end-start) can lose an ULP through cancellation; only a
        // genuinely partial observation needs reconstructed geometry.
        if (fraction < 1.) {
            geometry.end_m = geometry.start_m + fraction * (geometry.end_m - geometry.start_m);
            if (geometry.end_height_m)
                geometry.end_height_m = geometry.height_m + fraction * (*geometry.end_height_m - geometry.height_m);
        }
        const Eigen::Vector3d delta = geometry.end_m - geometry.start_m;
        if (entry.commanded_area_m2)
            result.commanded_path_volume_m3 += *entry.commanded_area_m2 * delta.norm();
        else
            ++result.segments_without_commanded_volume;
        const double h0 = geometry.height_m, h1 = geometry.end_height_m.value_or(h0);
        const double area = geometry.width_m * .5*(h0+h1) -
            (geometry.section == BeadSection::RoundedRectangle ?
                (1. - std::acos(-1.) / 4.) * (h0*h0+h0*h1+h1*h1)/3. : 0.);
        result.nominal_section_path_volume_m3 += area * delta.head<2>().norm();
        segments.push_back(geometry);
        retained_fractions.push_back(fraction);
    }
    if (segments.empty()) throw std::invalid_argument("No deposition exists at the observation time.");
    if (stationary_allowance && result.stationary_extrusion.commanded_volume_m3>
        stationary_allowance->maximum_total_fraction*result.commanded_path_volume_m3)
        throw std::invalid_argument("Stationary extrusion exceeds its aggregate volume allowance.");
    if (result.collapsed_envelope.count && (!collapsed_allowance || result.segments_without_commanded_volume ||
        result.collapsed_envelope.commanded_volume_m3 > collapsed_allowance->maximum_total_fraction * result.commanded_path_volume_m3))
        throw std::invalid_argument("Collapsed envelope exceeds its aggregate volume fraction allowance.");
    return result;
}

TimedToolpath timed_deposition(const GCodeProcessorResult &source,
    std::size_t max_moves, std::size_t max_segments, const std::function<void()> &checkpoint)
{
    TimedToolpath result;
    const auto fail = [&](ToolpathState state, ToolpathFailure failure) {
        result.state = state;
        result.failure = failure;
        result.deposition.clear();
        result.unresolved_deposition.clear();
        return result;
    };
    if (max_moves == 0 || max_moves > 2000000 || max_segments > 2000000 || source.moves.size() > max_moves)
        return fail(ToolpathState::BudgetExceeded, ToolpathFailure::MoveBudget);
    if (source.deposition_z_spans.size() > GCodeProcessorResult::maximum_deposition_z_spans)
        return fail(ToolpathState::BudgetExceeded, ToolpathFailure::MoveBudget);
    for (std::size_t i = 0; i < source.deposition_z_spans.size(); ++i) {
        if (checkpoint && i % 1024 == 0) checkpoint();
        const auto index = source.deposition_z_spans[i].move_index;
        if (index >= source.moves.size() || (i && source.deposition_z_spans[i - 1].move_index >= index) ||
            source.moves[index].internal_only || source.moves[index].type != EMoveType::Extrude)
            return fail(ToolpathState::InvalidInput, ToolpathFailure::InvalidBeadDimensions);
    }
    constexpr std::size_t mode = static_cast<std::size_t>(PrintEstimatedStatistics::ETimeMode::Normal);
    const double expected = source.print_statistics.modes[mode].time;
    if (!std::isfinite(expected) || expected <= 0.) return fail(ToolpathState::MissingTiming, ToolpathFailure::MissingTotal);
    std::size_t previous = 0;
    bool have_previous = false;
    for (std::size_t i = 0; i < source.moves.size(); ++i) {
        if (checkpoint && i % 1024 == 0) checkpoint();
        const auto &move = source.moves[i];
        result.failing_move_index = i;
        if (have_previous) result.previous_physical_move_index = previous;
        // Display-only acceleration subdivisions have zero duration and must
        // not shorten the original deposition path or invent activation times.
        if (move.internal_only) {
            if (move.deposition_vertical_extrusion)
                return fail(ToolpathState::InvalidInput,ToolpathFailure::InvalidVerticalDeposition);
            if (move.time[mode] != 0.) return fail(ToolpathState::InvalidInput, ToolpathFailure::InternalVertexDuration);
            continue;
        }
        const double duration = move.time[mode];
        if (!std::isfinite(duration) || duration < 0. || !move.position.allFinite())
            return fail(ToolpathState::InvalidInput, ToolpathFailure::InvalidDurationOrPosition);
        const double start = result.total_s;
        result.total_s += duration;
        if (!std::isfinite(result.total_s)) return fail(ToolpathState::InvalidInput, ToolpathFailure::TimeOverflow);
        // Seam markers can sit slightly off the physical endpoint. Retain any
        // attributed legacy timing, but never use their display position as a
        // deposition origin.
        if (move.deposition_vertical_extrusion) {
            if (!have_previous || !authenticated_vertical_deposition(source,i,previous))
                return fail(ToolpathState::InvalidInput,ToolpathFailure::InvalidVerticalDeposition);
            if (duration<=0.) return fail(ToolpathState::MissingTiming,ToolpathFailure::MissingExtrusionDuration);
            if (result.deposition.size()+result.unresolved_deposition.size()>=max_segments)
                return fail(ToolpathState::BudgetExceeded,ToolpathFailure::SegmentBudget);
            // The emitter authenticated positive E and the actual pure-Z span,
            // but supplied no XY volume model. Retain commanded filament/time;
            // preview Travel dimensions and mm3_per_mm cannot invent a bead.
            result.unresolved_deposition.push_back({source.moves[previous].position.cast<double>()*.001,
                move.position.cast<double>()*.001,start,result.total_s,move.delta_extruder,move.extruder_id,i,
                move.gcode_id,move.extrusion_role,UnresolvedSpatialDeposition::Reason::UnresolvedXYFootprint,
                move.deposition_source_id,move.deposition_command_section_id});
            result.spatial_deposition_resolved=false;
            previous=i;
            have_previous=true;
            continue;
        }
        if (move.type == EMoveType::Seam) continue;
        if (move.type == EMoveType::Extrude) {
            if (!have_previous || duration <= 0.) return fail(ToolpathState::MissingTiming, ToolpathFailure::MissingExtrusionDuration);
            if (!std::isfinite(move.width) || move.width <= 0. || !std::isfinite(move.height) || move.height <= 0.)
                return fail(ToolpathState::InvalidInput, ToolpathFailure::InvalidBeadDimensions);
            if (!std::isfinite(move.fan_speed) || move.fan_speed < 0. || move.fan_speed > 100. || !std::isfinite(move.temperature))
                return fail(ToolpathState::InvalidInput, ToolpathFailure::InvalidFanOrTemperature);
            if (!std::isfinite(move.mm3_per_mm) || move.mm3_per_mm < 0.)
                return fail(ToolpathState::InvalidInput, ToolpathFailure::InvalidCommandedArea);
            Eigen::Vector3d from = source.moves[previous].position.cast<double>() * .001;
            Eigen::Vector3d to = move.position.cast<double>() * .001;
            bool precise_anchored = false;
            try {
                if (const auto z = authenticated_deposition_z(source, i, previous)) {
                    from.z() = (*z)[0] * .001;
                    to.z() = (*z)[1] * .001;
                    precise_anchored = true;
                }
            } catch (const std::invalid_argument &) {
                return fail(ToolpathState::InvalidInput, ToolpathFailure::InvalidBeadDimensions);
            }
            // Nominal sloped beads are supported; purely vertical extrusion has
            // no supported XY footprint in the deposition rasterizer.
            if (result.deposition.size() + result.unresolved_deposition.size() >= max_segments)
                return fail(ToolpathState::BudgetExceeded, ToolpathFailure::SegmentBudget);
            const bool collapsed = precise_anchored && move.deposition_design_heights &&
                (*move.deposition_design_heights)[0] == 0. && (*move.deposition_design_heights)[1] == 0.;
            if (collapsed || (to - from).head<2>().squaredNorm() == 0.) {
                if (!std::isfinite(move.delta_extruder) || move.delta_extruder <= 0.)
                    return fail(ToolpathState::InvalidInput, ToolpathFailure::UnresolvedXYFootprint);
                result.unresolved_deposition.push_back({from, to, start, result.total_s,
                    move.delta_extruder, move.extruder_id, i, move.gcode_id, move.extrusion_role,
                    collapsed ? UnresolvedSpatialDeposition::Reason::CollapsedEnvelope :
                        UnresolvedSpatialDeposition::Reason::UnresolvedXYFootprint,
                    move.deposition_source_id, move.deposition_command_section_id});
                result.spatial_deposition_resolved = false;
                previous = i;
                have_previous = true;
                continue;
            }
            result.deposition.push_back({{from, to, move.width * .001, move.height * .001, move.extruder_id},
                start, result.total_s, i, move.gcode_id, move.temperature_targets_id,
                move.extrusion_role, move.fan_speed, move.temperature});
            if (move.deposition_design_heights) {
                const auto &heights = *move.deposition_design_heights;
                if (!source.deposition_sections || !std::isfinite(heights[0]) || !std::isfinite(heights[1]) ||
                    heights[0] < 0.f || heights[1] < 0.f || std::max(heights[0], heights[1]) <= 0.f)
                    return fail(ToolpathState::InvalidInput, ToolpathFailure::InvalidBeadDimensions);
                result.deposition.back().geometry.height_m = double(heights[0]) * .001;
                result.deposition.back().geometry.end_height_m = double(heights[1]) * .001;
                if (!move.deposition_design_width || !std::isfinite(*move.deposition_design_width) || *move.deposition_design_width <= 0.f)
                    return fail(ToolpathState::InvalidInput, ToolpathFailure::InvalidBeadDimensions);
                result.deposition.back().geometry.width_m = double(*move.deposition_design_width) * .001;
            }
            if (move.mm3_per_mm > 0.) result.deposition.back().commanded_area_m2 = move.mm3_per_mm * 1e-6;
        }
        previous = i;
        have_previous = true;
    }
    // Float move durations versus double processor total: allow accumulation
    // rounding only, never fill a missing interval from distance/feedrate.
    const double tolerance = std::max(1e-5, expected * 8. * std::numeric_limits<float>::epsilon());
    result.failing_move_index.reset();
    result.previous_physical_move_index.reset();
    if (std::abs(result.total_s - expected) > tolerance) return fail(ToolpathState::MissingTiming, ToolpathFailure::TotalMismatch);
    result.state = ToolpathState::Complete;
    return result;
}
} // namespace Slic3r::Thermoelastic
