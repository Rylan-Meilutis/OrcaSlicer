#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace Slic3r {

// Generation-local source identity, not a certificate of physical material ownership.
// All IDs refer to the immutable native slicing snapshot, never to object pointers.
struct DepositionSource
{
    enum class Height : std::uint8_t { Unknown, NativePlanar, Transformed };
    enum class Assignment : std::uint8_t { Unknown, NominalRole, WipingOverride };
    std::uint64_t object           = 0;
    std::uint64_t instance         = 0;
    std::uint64_t layer            = 0;
    std::uint64_t region           = 0;
    std::uint32_t logical_filament = 0;
    Height height                  = Height::Unknown;
    // Authenticated emission pass after per-instance wiping overrides. This
    // identifies the nominal role assignment, not measured material purity.
    Assignment assignment          = Assignment::Unknown;
};

class DepositionProvenance
{
public:
    static constexpr std::size_t maximum_sources = 65536;
    explicit DepositionProvenance(std::size_t capacity = maximum_sources);
    DepositionProvenance(const DepositionProvenance&)            = delete;
    DepositionProvenance& operator=(const DepositionProvenance&) = delete;
    // One-based IDs; zero means explicitly unowned. No deduplication by geometry.
    std::uint32_t add(DepositionSource source);
    void seal() noexcept { m_sealed = true; }
    bool sealed() const noexcept { return m_sealed; }
    std::string begin(std::uint32_t id) const;
    std::string clear() const;
    std::uint64_t generation() const noexcept { return m_generation; }
    const DepositionSource* source(std::uint32_t id) const noexcept;

private:
    std::uint64_t m_generation;
    std::size_t m_capacity;
    bool m_sealed = false;
    std::vector<DepositionSource> m_sources;
};

// Consumes complete internal lines AFTER string processors and BEFORE writing or
// parsing public G-code. Recognized control lines must never reach the output file.
// Imported comments cannot create a registry. Integration must prove marker
// preservation through every enabled processor. An independently captured expected
// transcript detects changed control-event sequences, but not ordinary lines moving
// across scopes or exchanging two occurrences with identical source IDs.
class DepositionProvenanceDecoder
{
public:
    explicit DepositionProvenanceDecoder(std::shared_ptr<const DepositionProvenance> registry);
    static constexpr std::size_t maximum_events = 1048576;
    // Copies and seals the expected source-ID sequence before decoding. Zero is
    // an explicit clear. Capture this from emission, NOT from filtered output.
    DepositionProvenanceDecoder(std::shared_ptr<const DepositionProvenance> registry, const std::vector<std::uint32_t>& expected_events);
    bool consume(std::string_view line);
    std::optional<std::uint32_t> source_id(std::uint32_t final_logical_filament);
    void finish() noexcept;
    bool valid() const noexcept { return m_valid; }
    // Existing no-transcript mode never certifies completeness. Published move
    // provenance must remain provisional until this returns true after finish().
    bool transcript_verified() const noexcept { return m_finished && m_valid && m_expected.has_value(); }

private:
    std::shared_ptr<const DepositionProvenance> m_registry;
    std::uint32_t m_active = 0;
    bool m_valid           = true;
    bool m_finished        = false;
    std::optional<std::vector<std::uint32_t>> m_expected;
    std::size_t m_event_index = 0;
};

// Output-thread capture while the producer's registry is still mutable. IDs are
// provisional and MUST NOT authorize ownership until finish succeeds and every
// attributed move's final logical filament matches the sealed source. No registry
// is read while capturing; reserved lines are stripped even after invalidation.
class DepositionProvenanceCapture
{
public:
    explicit DepositionProvenanceCapture(std::uint64_t generation, std::size_t event_capacity = DepositionProvenanceDecoder::maximum_events);
    bool consume(std::string_view line);
    std::uint32_t source_id() const noexcept { return m_valid && !m_finished ? m_active : 0; }
    bool valid() const noexcept { return m_valid; }
    bool finish(std::shared_ptr<const DepositionProvenance> registry, const std::vector<std::uint32_t>& expected_events);

private:
    std::uint64_t m_generation;
    std::size_t m_capacity;
    std::vector<std::uint32_t> m_events;
    std::uint32_t m_active = 0;
    bool m_valid           = true;
    bool m_finished        = false;
};

} // namespace Slic3r
