#pragma once

#include "ExtrusionEntity.hpp"
#include "libslic3r.h"

namespace Slic3r::detail {
// Internal course geometry planner, exposed for focused ramp regressions.
void smooth_staggered_transitions(ExtrusionPaths &paths, coord_t z_offset, bool closed);
}
