#pragma once
#include "ThermoelasticGCodeSections.hpp"
#include "ThermoelasticSections.hpp"
#include "ThermoelasticSectionStudy.hpp"
#include "ThermoelasticLayers.hpp"
#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace Slic3r {
struct GCodeProcessorResult;
namespace Thermoelastic {
struct SectionLayerPreview {
    enum class State { Ready, MissingLayer, UnsupportedInput, ResourceLimited, Unconverged, SolverFailure };
    State state = State::UnsupportedInput;
    unsigned layer_id = 0;
    double observation_time_s = 0., tensile_generalized_bed_force_n = 0., elastic_energy_j = 0.;
    std::vector<SectionSpatialSample> samples; // id is the processor move index.
    std::size_t unsampled_moves = 0;
    std::vector<SectionBasisSample> basis_history;
    std::vector<SectionGeometrySample> geometry_history;
    bool approximate_geometry = false;
    std::size_t geometry_tiles = 0, geometry_slabs = 0, section_count = 0, segment_visits = 0, boundary_points = 0;
    std::string detail;
};
// Static prescribed-eigenstrain field on geometry present at the selected
// modeled layer end. NOT thermal history, predicted peeling or local traction.
// Only selected-layer moves are sampled, independently of preview visibility.
// Results publish only after geometry and basis criteria pass; opt-in integration
// partition and boundary-representation approximations retain their refinement
// history. These integrated criteria do not independently qualify peak local
// stresses, and do not alter the original deposition endpoints.
SectionLayerPreview preview_section_layer(const GCodeProcessorResult& source,
    unsigned layer_id, const GCodeSectionScenario& scenario,
    std::size_t max_moves = 2000000, std::size_t max_segments = 2000000,
    std::size_t max_samples = 200000,
    const std::function<std::optional<SectionOwnership>(double)>& ownership = {},
    const std::function<void(unsigned, unsigned)>& progress = {},
    const std::function<void()>& checkpoint = {});
} // namespace Thermoelastic
} // namespace Slic3r
