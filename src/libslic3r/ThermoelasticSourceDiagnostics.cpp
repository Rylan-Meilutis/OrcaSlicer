#include "ThermoelasticSourceDiagnostics.hpp"
#include "GCode/DepositionProvenance.hpp"
#include "GCode/GCodeProcessor.hpp"
#include "ThermoelasticMesh.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <stdexcept>

namespace Slic3r::Thermoelastic {
SourceGeometryDiagnostics source_geometry_diagnostics(const GCodeProcessorResult& result,
                                                      std::size_t max_moves,
                                                      const std::function<void()>& checkpoint)
{
    if (max_moves == 0 || max_moves > 2000000)
        throw std::invalid_argument("Invalid source geometry diagnostic budget.");
    if (result.moves.size() > max_moves)
        throw ResourceLimit("Source geometry diagnostic move budget exceeded.");
    SourceGeometryDiagnostics out;
    if (!result.deposition_provenance)
        return out;
    if (!result.deposition_provenance->sealed())
        throw std::invalid_argument("Source geometry diagnostics require a sealed registry.");
    std::size_t source_count = 0;
    while (result.deposition_provenance->source(static_cast<std::uint32_t>(source_count + 1)))
        ++source_count;
    out.sources.resize(source_count + 1);
    std::optional<std::size_t> previous;
    const auto extend = [](std::optional<SourceGeometryDiagnostics::Range>& range, double value) {
        if (!std::isfinite(value))
            return;
        if (!range)
            range = SourceGeometryDiagnostics::Range{value, value};
        else {
            range->minimum = std::min(range->minimum, value);
            range->maximum = std::max(range->maximum, value);
        }
    };
    for (std::size_t i = 0; i < result.moves.size(); ++i) {
        if (checkpoint)
            checkpoint();
        ++out.visited_moves;
        const auto& move = result.moves[i];
        if (move.internal_only || move.type == EMoveType::Seam)
            continue;
        const auto id = move.deposition_source_id;
        if (id >= out.sources.size())
            throw std::invalid_argument("Unknown source geometry provenance ID.");
        if (id != 0 && result.deposition_provenance->source(id)->logical_filament != move.extruder_id)
            throw std::invalid_argument("Source geometry filament does not match provenance.");
        auto& source = out.sources[id];
        if (move.type == EMoveType::Extrude) {
            ++source.extrusion_moves;
            extend(source.top_z_m, double(move.position.z()) * .001);
            if (previous)
                extend(source.top_z_m, double(result.moves[*previous].position.z()) * .001);
            extend(source.height_m, double(move.height) * .001);
            const bool dimensions_valid = std::isfinite(move.height) && move.height > 0.f && std::isfinite(move.width) && move.width > 0.f;
            if (!previous || !dimensions_valid || !move.position.allFinite() || !result.moves[*previous].position.allFinite()) {
                ++source.stationary_or_unresolved_moves;
            } else {
                const auto& start = result.moves[*previous].position;
                if (start.x() == move.position.x() && start.y() == move.position.y())
                    ++source.stationary_or_unresolved_moves; // Includes vertical-only extrusion.
                else if (start.z() == move.position.z())
                    ++source.horizontal_moves;
                else
                    ++source.sloped_moves;
            }
        } else if (move.delta_extruder > 0.f) {
            // Positive E on travel/unretract/custom records has no certified
            // spatial bead section. Inventory it; never invent a deposition path.
            ++source.stationary_or_unresolved_moves;
        }
        previous = i;
    }
    out.state = SourceGeometryDiagnostics::State::Complete;
    return out;
}
} // namespace Slic3r::Thermoelastic
