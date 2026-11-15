#pragma once

#include "SectionThermalExposure.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include "ThermoelasticExactSections.hpp"
#include "ExPolygon.hpp"
#include "libslic3r.h"
#include <cstddef>
#include <array>
#include <functional>
#include <memory>
#include <vector>

namespace Slic3r::Thermoelastic {
struct LocalSectionTile;
// Explicit offline callers may allow up to2B cumulative CPU work; this does
// not enlarge the separate live-point or seam-frontier storage limits.
struct SectionSkinLimits
{ std::size_t max_live_points = 200000, max_seam_patches = 100000, max_work = 20000000; };
// Numerical exposed-skin visitor for a lexicographic, certified topology
// stream. Callbacks are provisional until that SAME topology stream finishes;
// labels must then be mapped through its connectivity certificate. Holds only
// adjacent Z bands and tile-seam traces, not historical polygons. Tile/material
// interfaces are removed, but true holes and unsupported undersides remain.
class SectionThermalSkin
{
public:
    using Horizontal = std::function<void(std::size_t material, std::size_t component, const ExPolygons&, coord_t plane, bool upward)>;
    using Vertical   = std::function<void(std::size_t material, std::size_t component, const VerticalSectionPatch&)>;
    // Optional original-label representatives from the SAME topology stream.
    // Queried only to reclaim a full frontier, never to infer new connections.
    using ComponentMapping = std::function<std::vector<std::size_t>()>;
    using ExactHorizontal = std::function<void(std::size_t, std::size_t, const ExactSection&, coord_t, bool)>;
    using ExactVertical = std::function<void(std::size_t, std::size_t, const ExactSectionBoundary&,
        coord_t, coord_t, const std::array<double,2>&)>;
    // Supplying both exact callbacks selects exact processing for the entire
    // stream, including native bands. This keeps all queued seam geometry in
    // one representation; the default native path remains unchanged.
    SectionThermalSkin(Horizontal horizontal,
                       Vertical vertical,
                       const SectionSkinLimits& limits  = {},
                       std::function<void()> checkpoint = {},
                       ComponentMapping component_mapping = {},
                       ExactHorizontal exact_horizontal = {}, ExactVertical exact_vertical = {});
    ~SectionThermalSkin();
    SectionThermalSkin(const SectionThermalSkin&)            = delete;
    SectionThermalSkin& operator=(const SectionThermalSkin&) = delete;
    void begin_tile(const LocalSectionTile& header);
    // Complete disjoint material partition for ONE common increasing Z band,
    // carrying the provisional labels emitted by LocalTopologyBuilder.
    void add_band(std::vector<SectionPiece> pieces);
    void end_tile();
    void finish();

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};
} // namespace Slic3r::Thermoelastic
