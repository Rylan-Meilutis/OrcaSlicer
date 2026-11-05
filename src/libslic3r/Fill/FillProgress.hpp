#pragma once

#include <cstddef>
#include <functional>

namespace Slic3r {

enum class FillProgressStage { GroupSurfaces, GenerateArcs, RecursiveArcs, ValidateArcs, RefineArcs, Complete };
using FillProgressCallback = std::function<void(FillProgressStage, size_t, size_t)>;

} // namespace Slic3r
