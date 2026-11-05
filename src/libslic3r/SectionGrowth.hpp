#pragma once
#include "ThermoelasticSectionTopology.hpp"
#include <cstddef>
#include <functional>
#include <vector>

namespace Slic3r::Thermoelastic {
struct SectionGrowthLimits {
    std::size_t max_work = 20000000, max_points = 2000000, max_pieces = 10000;
};
// Certificate of monotone material-preserving growth between two qualified
// snapshots. Components may join, but existing material cannot move, disappear,
// change ownership or split. Added pieces retain the NEXT snapshot's component
// IDs; unlike a complete topology, the added subset need not touch the bed.
class SectionGrowth {
public:
    const std::vector<SectionPiece>& added_pieces() const { return m_added; }
    const std::vector<std::size_t>& previous_components() const { return m_components; }
    double bed_z_m() const { return m_bed; }
private:
    SectionGrowth() = default;
    std::vector<SectionPiece> m_added;
    std::vector<std::size_t> m_components;
    double m_bed = 0.;
    friend SectionGrowth qualify_section_growth(const SectionTopology&, const SectionTopology&,
        const SectionGrowthLimits&, const std::function<void()>&);
};
SectionGrowth qualify_section_growth(const SectionTopology& previous, const SectionTopology& next,
    const SectionGrowthLimits& limits = {}, const std::function<void()>& checkpoint = {});

struct ReferenceSectionActivation {
    SectionTopology geometry;
    // Nominal source volume outside its fixed material region. This is an
    // explicit geometric approximation, not vanished physical mass or heat.
    double excluded_nominal_volume_m3 = 0.;
};
// Fixed-interface reference model, NOT a replacement for chronological bead
// deformation. Source topologies may overlap each other, but each is certified
// independently. Activate only their intersection with SAME-material reference
// regions. Never activate future paths just because the final region exists.
// Reuse the identical reference across events and still qualify actual growth.
// A final snapshot must separately prove coverage of the full reference.
ReferenceSectionActivation activate_reference_sections(const SectionTopology& reference,
    const std::vector<SectionTopology>& deposited, const SectionGrowthLimits& limits = {},
    const SectionTopologyLimits& topology_limits = {}, const std::function<void()>& checkpoint = {});
} // namespace Slic3r::Thermoelastic
