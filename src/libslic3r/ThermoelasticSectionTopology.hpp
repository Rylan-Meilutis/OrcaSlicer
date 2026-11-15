#pragma once

#include "ThermoelasticLayers.hpp"
#include <cstddef>
#include <functional>
#include <vector>

namespace Slic3r::Thermoelastic {

struct LocalSectionTile;

struct SectionPiece {
    DepositedSlab slab;
    std::size_t component;
};
struct SectionTopologyLimits {
    std::size_t max_slabs = 10000, max_boundary_points = 2000000;
    std::size_t max_component_pairs = 2000000;
    // Exact geometry only: certify occupied interior faces separately rather
    // than refusing every boundary pinch. Point contacts never create bonds.
    // Default/native topology and ownership qualification are unchanged.
    bool exact_interior_components = false;
};

// Only the factory can create a certificate. It applies to this normalized
// piecewise-constant geometry, NOT the unsampled physical extrusion volume.
class SectionTopology {
public:
    const std::vector<SectionPiece> &pieces() const { return m_pieces; }
    std::size_t components() const { return m_components; }
    double bed_z_m() const { return m_bed_z_m; }
private:
    SectionTopology() = default;
    std::vector<SectionPiece> m_pieces;
    std::size_t m_components = 0;
    double m_bed_z_m = 0.;
    friend SectionTopology qualify_section_topology(const std::vector<DepositedSlab> &,
        double, const SectionTopologyLimits &, const std::function<void()> &);
    friend SectionTopology qualify_local_section_topology(const std::vector<LocalSectionTile> &,
        double, const SectionTopologyLimits &, const std::function<void()> &);
};

// Requires disjoint Z intervals (identical intervals for different materials
// are allowed). Within an interval, unions same-material polygons and rejects
// cross-material overlap. Different intervals connect only through a common
// Z face with strictly positive occupied overlap area, never edge/point contact.
// Every resulting solid component must have positive bed contact; no pins,
// gap filling, or bonds across voids are introduced. Native-only input retains
// integer clipping. If any footprint is canonical exact, all XY unions,
// components and contacts retain rational geometry; Z remains on the native
// nanometre grid. Nonmanifold boundary contacts reject by default; explicit
// interior-face decomposition separates point-touching solids without bonds.
// Boundary uncertainty still requires local topology qualification. Exact
// operation work shares max_component_pairs; retained vertices share the point
// limits. Returned pieces are disjoint and retain material IDs.
SectionTopology qualify_section_topology(const std::vector<DepositedSlab> &slabs,
    double bed_z_m, const SectionTopologyLimits &limits = {},
    const std::function<void()> &checkpoint = {});

} // namespace Slic3r::Thermoelastic
