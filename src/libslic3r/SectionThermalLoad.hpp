#pragma once
#include "SectionThermalAssembly.hpp"
#include "ThermalDeposition.hpp"
#include "ThermoelasticSections.hpp"
#include "ThermoelasticThermalStrain.hpp"
#include <Eigen/Core>
#include <functional>
#include <vector>

namespace Slic3r::Thermoelastic {
struct SectionExpansionMaterial {
    std::vector<ExpansionInterval> intervals;
    Eigen::Matrix3d local_to_world = Eigen::Matrix3d::Identity();
    double stress_free_temperature_k = 0.;
};
// Owns the converted field; does not borrow the history or assembly. Only a
// successful complete history is accepted. Consumer MUST use the same qualified
// topology and grid: component IDs are snapshot-local, not persistent IDs.
// Query bounds/volume are cross-checked but are not a geometry fingerprint.
// Rejects expansion-law breakpoints within a cell instead of fitting them.
SectionLoadField section_thermal_load(const SectionThermalAssembly& assembly,
    const ThermalHistory::HistoryResult& history, const std::vector<SectionExpansionMaterial>& materials,
    const std::function<void()>& checkpoint = {});
} // namespace Slic3r::Thermoelastic
