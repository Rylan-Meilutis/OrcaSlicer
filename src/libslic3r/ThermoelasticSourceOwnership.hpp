#pragma once

#include "GCode/DepositionProvenance.hpp"
#include "Point.hpp"
#include "ExPolygon.hpp"
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticMesh.hpp"
#include <cstddef>
#include <cstdint>
#include <array>
#include <functional>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace Slic3r {
class Print;
struct GCodeProcessorResult;
namespace Thermoelastic {

struct SourceOwnershipLimits
{
    std::size_t max_sources         = 4096;
    std::size_t max_boundary_points = 2000000;
    std::size_t max_work            = 1000000;
    // Explicit analysis-only source sensitivity. Does not alter native masks
    // or authenticate another command; strict coverage is the default.
    unsigned centerline_boundary_allowance_nm = 0;
};

struct SourceBead
{
    DepositedSegment geometry;
    std::uint32_t source_id        = 0;
    std::uint32_t logical_filament = 0;
    std::optional<std::size_t> move_index;
};

// Authenticate one complete, spatially resolved live move, independently of
// native material masks. Mapped geometry.material is the caller's responsibility.
// Transformed sources require typed command evidence; ordinary NativePlanar
// moves may use unchanged processor dimensions. Source0 is admitted only for
// typed AuxiliaryBrim planar LINE commands; it remains explicitly non-CAD.
// No partial-command truncation,
// unresolved extrusion, native coverage, or chronological-order claim is made.
// Returns work to debit from the caller's budget; cancellation propagates.
std::size_t authenticate_live_deposited_source(const GCodeProcessorResult& result, const SourceBead& bead,
    std::size_t max_work = 1000000, const std::function<void()>& checkpoint = {});

class NativeOwnershipSnapshot
{
public:
    // Opaque successful XY predicates, never source authentication or sections.
    // A cache belongs to one snapshot and one serial prepared resolver.
    class ContainmentCache {
        friend class NativeOwnershipSnapshot;
        explicit ContainmentCache(std::shared_ptr<const unsigned char> identity) : m_identity(std::move(identity)) {}
        std::shared_ptr<const unsigned char> m_identity;
        using Key = std::array<std::int64_t, 7>;
        // Direct-mapped slots bound storage exactly; collisions merely evict a
        // certificate. Full key equality, never the hash, establishes a hit.
        std::array<std::optional<Key>, 16384> m_covered;
        static std::size_t slot(const Key& key)
        {
            std::uint64_t hash = 14695981039346656037ull;
            for (const auto value : key) {
                hash ^= static_cast<std::uint64_t>(value);
                hash *= 1099511628211ull;
            }
            return (hash ^ (hash >> 32)) % 16384;
        }
    };
    std::shared_ptr<ContainmentCache> make_containment_cache() const
    {
        return std::shared_ptr<ContainmentCache>(new ContainmentCache(m_identity));
    }
    enum class SourceStatus { Missing, Unsupported, NativePlanar };
    SourceStatus source_status(std::uint32_t id) const noexcept
    {
        return id == 0 || id > m_masks.size()                    ? SourceStatus::Missing :
               m_masks[id - 1] && m_masks[id - 1]->native_planar ? SourceStatus::NativePlanar :
                                                                   SourceStatus::Unsupported;
    }
    // Supply every bead contributing to the conflict being resolved. Unrelated
    // unowned deposition need not be included; contributors must never be omitted.
    SectionOwnership section(std::shared_ptr<const DepositionProvenance> registry,
                             const std::vector<SourceBead>& active_beads,
                             double z_m,
                             const std::vector<std::size_t>& filament_to_material,
                             const SourceOwnershipLimits& limits          = {},
                             const std::function<void()>& checkpoint      = {},
                             const GCodeProcessorResult* emitted_evidence = nullptr,
                             std::size_t* used_work = nullptr,
                             ContainmentCache* containment_cache = nullptr) const;

private:
    struct Reference {
        double bottom_m, top_m;
        ExPolygons polygons;
        std::size_t identity;
    };
    struct Mask
    {
        std::uint32_t logical_filament;
        double bottom_m, top_m;
        std::shared_ptr<const Reference> reference;
        // Transformed masks require authenticated local scarf geometry or a
        // live eligible planar-infill move. Neither extends them into adjacent Z.
        bool native_planar = false;
        // Nominal infill roles with unmixed material. Settings that could alter
        // other paths require authenticated nominal-pass command evidence;
        // each live move must still prove the exact native plane/thickness.
        std::uint32_t planar_infill_roles = 0;
        bool infill_requires_command = false;
        // Only typed, authenticated BrickCourse commands may select one of
        // these actual adjacent same-material intervals, never extend a mask.
        bool brick_eligible = false;
        std::array<std::shared_ptr<const Reference>, 2> brick_neighbors;
    };
    NativeOwnershipSnapshot() = default;
    std::shared_ptr<const unsigned char> m_identity = std::make_shared<const unsigned char>(0);
    std::shared_ptr<const DepositionProvenance> m_registry;
    std::vector<std::optional<Mask>> m_masks; // Empty means unsupported; IDs index at id-1.
    friend NativeOwnershipSnapshot capture_native_ownership(const Print&,
                                                            std::shared_ptr<const DepositionProvenance>,
                                                            const Point&,
                                                            double,
                                                            const SourceOwnershipLimits&,
                                                            const std::function<void()>&);
};

// MUST run inside the matching live export while Print is immutable, after the
// registry's control transcript and move filaments have been verified. Numeric
// IDs alone cannot authenticate an unrelated Print. Lifecycle integration must
// attach this snapshot and registry together; this factory does not establish it.
// Offsets explicitly map native instance coordinates to parsed deposition space.
// Only ordinary Cartesian output is supported; callers must exclude IMEX/custom
// coordinate transforms. No CAD solid or deposited envelope is generated here.
NativeOwnershipSnapshot capture_native_ownership(const Print& print,
                                                 std::shared_ptr<const DepositionProvenance> registry,
                                                 const Point& parsed_xy_offset,
                                                 double parsed_z_offset_m,
                                                 const SourceOwnershipLimits& limits     = {},
                                                 const std::function<void()>& checkpoint = {});

} // namespace Thermoelastic
} // namespace Slic3r
