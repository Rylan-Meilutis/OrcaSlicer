#pragma once

#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

namespace Slic3r {
struct GCodeProcessorResult;
namespace Thermoelastic {

struct SourceGeometryDiagnostics
{
    struct Range
    { double minimum, maximum; };
    struct Source
    {
        std::size_t extrusion_moves                = 0;
        std::size_t horizontal_moves               = 0;
        std::size_t sloped_moves                   = 0;
        std::size_t stationary_or_unresolved_moves = 0;
        std::optional<Range> top_z_m;
        std::optional<Range> height_m;
        bool height_varies() const noexcept { return height_m && height_m->minimum != height_m->maximum; }
    };
    enum class State { MissingRegistry, Complete };
    State state = State::MissingRegistry;
    // Index zero is unowned; remaining indices equal verified source IDs.
    std::vector<Source> sources;
    std::size_t visited_moves = 0;
};

// Diagnostics of parsed nominal paths only. No native-mask admission, physical
// thickness reconstruction, or transformed-source qualification is performed.
// Exact stored float comparisons intentionally expose even tiny Z/height changes.
// Internal display vertices and seam markers do not replace the physical start.
SourceGeometryDiagnostics source_geometry_diagnostics(const GCodeProcessorResult& result,
                                                      std::size_t max_moves                   = 2000000,
                                                      const std::function<void()>& checkpoint = {});

} // namespace Thermoelastic
} // namespace Slic3r
