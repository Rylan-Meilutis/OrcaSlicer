#pragma once
#include "ThermoelasticLocalSections.hpp"
#include "ThermoelasticLocalTopology.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include "ThermoelasticSections.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>
#include <memory>
#include <array>

namespace Slic3r::Thermoelastic {
struct SectionGrowthStreamLimits
{
    std::size_t max_work = 20000000, max_tile_points = 200000;
    std::size_t max_output_points = 8000000, max_output_pieces = 100000, max_component_links = 100000;
    // Explicit representation choice for fixed first-deposition references,
    // not ownership evidence or a resource-limit increase. Only activation's
    // per-material source cursors use this; the reference keeps its resolver.
    bool exact_reference_sources = false;
    // Explicit complete first-deposition reference only. Every forwarded
    // packet still requires its bound complete ownership proof; no native or
    // partial-prefix shortcut and no point-contact/physical qualification.
    bool reuse_complete_deposited_reference = false;
};
struct SectionGrowthStreamSettings : SectionGrowthStreamLimits
{
    Eigen::Vector2d origin_m = Eigen::Vector2d::Zero(), tile_size_m = Eigen::Vector2d::Zero();
    double bed_z_m = 0., maximum_height_m = 0.;
    unsigned interval_subdivisions = 1;
    LocalSectionLimits geometry;
    SectionTopologyLimits topology;
    StreamTopologyLimits connectivity;
};
struct SectionGrowthStreamVisitor
{
    // Borrowed provisional pieces; only successful return validates the
    // complete traversal. Added pieces use NEXT provisional component labels.
    std::function<void(const SectionPiece&)> previous, next, added;
    // Optional complete packet hooks for the NEXT stream, after its topology
    // callbacks. A consumer can finish surface bands using labeled pieces.
    LocalSectionPacketVisitor next_packets;
};
struct SectionGrowthStreamResult
{
    LocalConnectivityCertificate previous, next;
    std::vector<std::size_t> previous_components;
    LocalSectionMetrics previous_geometry, next_geometry;
};
struct ReferenceActivationStreamResult {
    LocalConnectivityCertificate reference, active;
    double excluded_nominal_volume_m3 = 0.;
    std::size_t processed_points = 0;
};
// Resumable form for pairing chronological snapshots. Borrows the immutable
// reference; deposited paths are owned by its per-material samplers. Budgets
// persist across advances. Exceptions/reentry poison the cursor; only finish()
// after exhaustion publishes connectivity certificates.
class ReferenceActivationCursor {
public:
    ReferenceActivationCursor(const std::vector<DepositedSegment>& reference,
        const std::vector<DepositedSegment>& deposited,const SectionGrowthStreamSettings& settings,
        std::function<void(const SectionPiece&)> visitor,
        std::function<std::optional<SectionOwnership>(double)> reference_ownership = {},
        std::function<void()> checkpoint = {},LocalSectionPacketVisitor active_packets = {});
    ~ReferenceActivationCursor();
    ReferenceActivationCursor(const ReferenceActivationCursor&) = delete;
    ReferenceActivationCursor& operator=(const ReferenceActivationCursor&) = delete;
    std::optional<std::array<int,2>> next_tile_index() const;
    bool advance();
    // Available in active packet hooks after the topology operation completes.
    // Representatives remain provisional and in the original label namespace.
    std::vector<std::size_t> active_component_mapping();
    ReferenceActivationStreamResult finish();
private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};
// Explicit fixed-interface approximation. Samples the immutable reference and
// each deposited material separately, retaining only their current tile. No
// source needs artificial bed attachment. Emitted pieces remain provisional
// until both topology streams finish. Keep the same reference across events;
// actual monotone growth and final reference coverage still require proof.
ReferenceActivationStreamResult visit_reference_activation(
    const std::vector<DepositedSegment>& reference, const std::vector<DepositedSegment>& deposited,
    const SectionGrowthStreamSettings& settings, const std::function<void(const SectionPiece&)>& visitor,
    const std::function<std::optional<SectionOwnership>(double)>& reference_ownership = {},
    const std::function<void()>& checkpoint = {}, const LocalSectionPacketVisitor& active_packets = {});
struct ReferenceGrowthStreamResult {
    LocalConnectivityCertificate previous, next;
    std::vector<std::size_t> previous_components;
    double previous_excluded_nominal_volume_m3 = 0., next_excluded_nominal_volume_m3 = 0.;
};
struct ReferenceSectionMechanicalResult {
    SectionMechanicalResult field;
    double excluded_nominal_volume_m3 = 0.;
};
ReferenceSectionMechanicalResult solve_reference_deposition_sections(
    const std::vector<DepositedSegment>& reference,const std::vector<DepositedSegment>& deposited,
    const std::vector<SectionMaterial>& materials,const LocalSectionStreamSettings& settings,
    const SectionGrowthStreamLimits& activation_limits = {},
    const std::function<std::optional<SectionOwnership>(double)>& reference_ownership = {},
    const std::function<void()>& checkpoint = {});
ReferenceGrowthStreamResult visit_reference_section_growth(
    const std::vector<DepositedSegment>& reference,const std::vector<DepositedSegment>& previous,
    const std::vector<DepositedSegment>& next,const SectionGrowthStreamSettings& settings,
    const SectionGrowthStreamVisitor& visitor,
    const std::function<std::optional<SectionOwnership>(double)>& reference_ownership = {},
    const std::function<void()>& checkpoint = {});
// Closed paired sampling/topology/growth traversal. Retains at most the two
// current tiles, not either complete snapshot. Checks polygon coverage by
// material in every Z interval; equal aggregate volume is not a certificate.
// No emitted callback is usable if this traversal later fails. Returned source
// certificates identify each side's labels; component mapping joins old nodes
// only at matching spatial positions, never merely by component membership.
SectionGrowthStreamResult visit_local_section_growth(const std::vector<DepositedSegment>& previous,
                                                     const std::vector<DepositedSegment>& next,
                                                     const SectionGrowthStreamSettings& settings,
                                                     const SectionGrowthStreamVisitor& visitor,
                                                     const std::function<std::optional<SectionOwnership>(double)>& previous_ownership = {},
                                                     const std::function<std::optional<SectionOwnership>(double)>& next_ownership     = {},
                                                     const std::function<void()>& checkpoint                                          = {},
                                                     const std::vector<DepositedSegment>* previous_partitions = nullptr,
                                                     const std::vector<DepositedSegment>* next_partitions     = nullptr);
} // namespace Slic3r::Thermoelastic
