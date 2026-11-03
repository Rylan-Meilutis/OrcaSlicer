#include "SectionThermalLoad.hpp"
#include "SectionThermalAssembly.hpp"
#include "ThermalDeposition.hpp"
#include "ThermalHistory.hpp"
#include "ThermoelasticSections.hpp"
#include "ThermoelasticThermalStrain.hpp"
#include <Eigen/Core>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
SectionLoadField section_thermal_load(const SectionThermalAssembly& assembly,
    const ThermalHistory::HistoryResult& history, const std::vector<SectionExpansionMaterial>& materials,
    const std::function<void()>& checkpoint)
{
    if (history.state != ThermalHistory::State::Converged || assembly.groups.empty() ||
        assembly.groups.size() > 100000 || assembly.nodes.size() > 200000 ||
        history.temperatures_k.size() != assembly.nodes.size() || materials.empty() || materials.size() > 256)
        throw std::invalid_argument("Incomplete section thermal history or invalid field limits.");
    using Key = std::tuple<std::array<int,3>,std::size_t,std::size_t>;
    struct Entry {
        Eigen::Vector3d minimum, maximum;
        double volume;
        TrilinearThermalStrain strain;
    };
    auto field = std::make_shared<std::map<Key,Entry>>();
    for (const auto& group : assembly.groups) {
        if (checkpoint) checkpoint();
        if (group.material >= materials.size() || !group.minimum_m.allFinite() || !group.maximum_m.allFinite() ||
            (group.maximum_m.array() <= group.minimum_m.array()).any() ||
            !std::isfinite(group.occupied_volume_m3) || group.occupied_volume_m3 <= 0.)
            throw std::invalid_argument("Invalid section thermal field group.");
        std::array<double,8> temperatures;
        for (unsigned corner = 0; corner < 8; ++corner) {
            const auto id = group.nodes[corner];
            if (id >= assembly.nodes.size() || assembly.nodes[id].component != group.component ||
                !assembly.nodes[id].position_m.allFinite())
                throw std::invalid_argument("Missing section thermal node coverage.");
            Eigen::Vector3d expected;
            for (unsigned axis = 0; axis < 3; ++axis)
                expected[axis] = corner & (1u << axis) ? group.maximum_m[axis] : group.minimum_m[axis];
            if ((expected-assembly.nodes[id].position_m).cwiseAbs().maxCoeff() > 1e-12)
                throw std::invalid_argument("Section thermal node does not match the field grid.");
            temperatures[corner] = history.temperatures_k[id];
        }
        const auto& material = materials[group.material];
        Entry entry{group.minimum_m,group.maximum_m,group.occupied_volume_m3,
            trilinear_thermal_eigenstrain(material.intervals,material.local_to_world,
                material.stress_free_temperature_k,temperatures)};
        if (!field->emplace(Key{group.cell,group.component,group.material},std::move(entry)).second)
            throw std::invalid_argument("Duplicate section thermal field group.");
    }
    // Const shared ownership keeps callback copies inexpensive and immutable.
    const std::shared_ptr<const std::map<Key,Entry>> owned = field;
    return [owned](const SectionLoadQuery& query) -> SectionSpatialStrain {
        const auto found = owned->find({query.cell,query.component,query.material});
        if (found == owned->end()) throw std::out_of_range("Missing section thermal field coverage.");
        const auto& entry = found->second;
        if (!query.minimum_m.allFinite() || !query.maximum_m.allFinite() ||
            !std::isfinite(query.occupied_volume_m3) ||
            (entry.minimum-query.minimum_m).cwiseAbs().maxCoeff() > 1e-12 ||
            (entry.maximum-query.maximum_m).cwiseAbs().maxCoeff() > 1e-12 ||
            std::abs(entry.volume-query.occupied_volume_m3) > 1e-10*entry.volume)
            throw std::invalid_argument("Thermal and mechanical section domains differ.");
        return entry.strain;
    };
}
} // namespace Slic3r::Thermoelastic
