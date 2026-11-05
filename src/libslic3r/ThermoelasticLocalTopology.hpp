#pragma once

#include "ThermoelasticSectionTopology.hpp"
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <vector>

namespace Slic3r::Thermoelastic {
struct LocalSectionTile;
class LocalTopologyBuilder;
enum class LocalTileOrder { Arbitrary, Lexicographic };
enum class LocalAttachmentPolicy { RequireBed, AllowFreeComponents };
struct StreamTopologyLimits {
    // Explicit offline consumers may reserve up to4 billion work units. This
    // bounds CPU operations only; all geometry and live-memory caps still apply.
    std::size_t max_processed_points = 8000000, max_work = 20000000;
    // max_islands independently bounds current-band islands, persistent labels
    // and retained face records; it does not count already released islands.
    std::size_t max_live_points = 200000, max_islands = 100000, max_trace_intervals = 2000000;
    std::size_t max_slabs = 10000;
};
// Partial admission counters only. Live points conservatively count borrowed
// and owned polygon payload, not Clipper workspace or allocator/RSS overhead.
struct StreamTopologyStatistics {
    std::size_t input_points = 0, processed_points = 0, work = 0, slabs = 0;
    std::size_t islands = 0, trace_intervals = 0, peak_live_points = 0;
    // Diagnostic counts, not stable geometry IDs: islands aliases face_records.
    std::size_t labels = 0, face_records = 0, processed_islands = 0;
    std::size_t peak_face_records = 0;
    // Ambiguous vertical edges may be redundant to independently certified
    // connectivity; this does not certify their local contact area.
    std::size_t uncertain_contacts = 0, redundant_contacts = 0;
};

// Identity of this topology stream, not proof of printer/G-code provenance.
class LocalTopologySource {
    friend class LocalTopologyBuilder;
    LocalTopologySource() = default;
};
class LocalConnectivityCertificate {
public:
    const std::vector<std::size_t>& component_mapping() const { return m_mapping; }
    std::size_t components() const { return m_components; }
    std::size_t free_components() const { return m_free_components; }
    double bed_z_m() const { return m_bed_z_m; }
    const std::shared_ptr<const LocalTopologySource>& source() const { return m_source; }
private:
    friend class LocalTopologyBuilder;
    LocalConnectivityCertificate() = default;
    std::vector<std::size_t> m_mapping;
    std::size_t m_components = 0;
    std::size_t m_free_components = 0;
    double m_bed_z_m = 0.;
    std::shared_ptr<const LocalTopologySource> m_source;
};

// Incremental common engine for retained and streamed tile topology. A piece
// callback borrows its polygon only during that synchronous call; its component
// is provisional until finish(). Earlier callbacks are partial, not qualified
// geometry. Consumers must match the stream source and final mapping before
// using accumulated data. No mechanical result is certified here.
// All existing cumulative topology caps apply across every tile. Any exception
// permanently poisons the builder; additions/reentry after failure or finish
// are rejected. No historical integration polygons are retained by this engine.
// Canonical exact footprints use the same labels and frontier, with rational
// unions, components, containment, positive-area Z contacts and tile-face spans.
// Exact nonmanifold point pinches reject rather than welding or rounding them.
// Exact geometry requires zero occupied uncertainty, either explicitly or from
// a complete unchanged authenticated occupied-partition proof. The latter keeps
// material-interface sensitivity on emitted pieces, without copying its group
// token to split components. Mixed contacts with uncertain native geometry reject.
// Exact frontier endpoints and
// bounded operation scratch count against live storage; no rounded proxy is kept.
// Lexicographic mode validates strictly increasing tile indices and joins faces
// online. Only positive faces targeting future tiles remain in a recyclable
// frontier; arbitrary-order mode retains faces until finish().
// Ambiguous vertical edges never join labels. Finish accepts them only when
// independent certified paths already connect their endpoints. Pending edges
// are capped at min(100000,max_islands) in streaming mode (component-pair limit
// otherwise). This resolves component identity, not local contact geometry.
class LocalTopologyBuilder {
public:
    LocalTopologyBuilder(double bed_z_m, const SectionTopologyLimits& limits,
        std::function<void(const SectionPiece&)> piece,
        std::function<void()> checkpoint = {}, std::optional<StreamTopologyLimits> stream_limits = {},
        StreamTopologyStatistics* statistics = nullptr, LocalTileOrder order = LocalTileOrder::Arbitrary,
        LocalAttachmentPolicy attachment = LocalAttachmentPolicy::RequireBed);
    ~LocalTopologyBuilder();
    LocalTopologyBuilder(const LocalTopologyBuilder&) = delete;
    LocalTopologyBuilder& operator=(const LocalTopologyBuilder&) = delete;
    void add_tile(const LocalSectionTile& tile);
    // Packet mode borrows one complete material partition per increasing Z
    // interval. An open tile must be ended before another tile or finish().
    void begin_tile(const LocalSectionTile& empty_header);
    void add_band(const LocalSectionTile& packet);
    void end_tile();
    // Current union-find representatives in the ORIGINAL label namespace.
    // Only between packet operations; not a completed connectivity certificate.
    // Consumers may consolidate provisional data but must still finish this
    // same source and apply its final certificate before using any result.
    std::vector<std::size_t> provisional_component_mapping();
    LocalConnectivityCertificate finish();
    const std::shared_ptr<const LocalTopologySource>& source() const { return m_source; }
private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
    std::shared_ptr<const LocalTopologySource> m_source;
};
} // namespace Slic3r::Thermoelastic
