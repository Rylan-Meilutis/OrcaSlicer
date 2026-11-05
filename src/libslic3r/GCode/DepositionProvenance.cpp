#include "DepositionProvenance.hpp"

#include <atomic>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace Slic3r {
namespace {
constexpr std::string_view prefix = ";_ORCA_INTERNAL_SOURCE:";
std::atomic<std::uint64_t> next_generation{1};

std::uint64_t new_generation()
{
    auto value = next_generation.load();
    do {
        if (value == std::numeric_limits<std::uint64_t>::max())
            throw std::overflow_error("Deposition provenance generation exhausted.");
    } while (!next_generation.compare_exchange_weak(value, value + 1));
    return value;
}

template<class Integer> bool parse(std::string_view text, Integer& value)
{
    if (text.empty() || text.front() < '0' || text.front() > '9')
        return false;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    return parsed.ec == std::errc{} && parsed.ptr == text.data() + text.size();
}

struct ControlLine
{
    bool reserved            = false;
    bool valid               = false;
    std::uint64_t generation = 0;
    std::uint32_t id         = 0;
};

ControlLine control_line(std::string_view line)
{
    ControlLine out;
    const auto reserved = prefix.substr(0, prefix.size() - 1);
    out.reserved        = line.substr(0, reserved.size()) == reserved;
    if (!out.reserved || line.substr(0, prefix.size()) != prefix)
        return out;
    line.remove_prefix(prefix.size());
    if (!line.empty() && line.back() == '\n')
        line.remove_suffix(1);
    if (!line.empty() && line.back() == '\r')
        line.remove_suffix(1);
    const auto separator = line.find(':');
    out.valid            = separator != std::string_view::npos && parse(line.substr(0, separator), out.generation) &&
                           parse(line.substr(separator + 1), out.id);
    return out;
}
} // namespace

DepositionProvenance::DepositionProvenance(std::size_t capacity) : m_generation(new_generation()), m_capacity(capacity)
{
    if (capacity == 0 || capacity > maximum_sources)
        throw std::invalid_argument("Invalid deposition provenance capacity.");
}

std::uint32_t DepositionProvenance::add(DepositionSource source)
{
    if (m_sealed)
        throw std::logic_error("Deposition provenance is sealed.");
    if (m_sources.size() == m_capacity)
        throw std::length_error("Deposition provenance source budget exceeded.");
    if (source.height != DepositionSource::Height::Unknown && source.height != DepositionSource::Height::NativePlanar &&
        source.height != DepositionSource::Height::Transformed)
        throw std::invalid_argument("Invalid deposition height provenance.");
    if (source.assignment != DepositionSource::Assignment::Unknown &&
        source.assignment != DepositionSource::Assignment::NominalRole &&
        source.assignment != DepositionSource::Assignment::WipingOverride)
        throw std::invalid_argument("Invalid deposition material assignment provenance.");
    m_sources.push_back(source);
    return static_cast<std::uint32_t>(m_sources.size());
}

const DepositionSource* DepositionProvenance::source(std::uint32_t id) const noexcept
{ return id != 0 && id <= m_sources.size() ? &m_sources[id - 1] : nullptr; }

std::string DepositionProvenance::begin(std::uint32_t id) const
{
    if (!source(id))
        throw std::invalid_argument("Unknown deposition source.");
    return std::string(prefix) + std::to_string(m_generation) + ":" + std::to_string(id) + "\n";
}

std::string DepositionProvenance::clear() const { return std::string(prefix) + std::to_string(m_generation) + ":0\n"; }

DepositionProvenanceDecoder::DepositionProvenanceDecoder(std::shared_ptr<const DepositionProvenance> registry)
    : m_registry(std::move(registry))
{
    if (!m_registry || !m_registry->sealed())
        throw std::invalid_argument("A sealed local deposition registry is required.");
}

DepositionProvenanceDecoder::DepositionProvenanceDecoder(std::shared_ptr<const DepositionProvenance> registry,
                                                         const std::vector<std::uint32_t>& expected_events)
    : DepositionProvenanceDecoder(std::move(registry))
{
    if (expected_events.size() > maximum_events)
        throw std::length_error("Deposition provenance event budget exceeded.");
    std::uint32_t active = 0;
    for (const auto id : expected_events) {
        if (id != 0 && (!m_registry->source(id) || active != 0))
            throw std::invalid_argument("Invalid expected deposition provenance scope.");
        active = id;
    }
    if (active != 0)
        throw std::invalid_argument("Unclosed expected deposition provenance scope.");
    m_expected = expected_events;
}

bool DepositionProvenanceDecoder::consume(std::string_view line)
{
    if (m_finished) {
        m_valid  = false;
        m_active = 0;
    }
    const auto control = control_line(line);
    if (!control.reserved)
        return false;
    const auto id = control.id;
    if (!control.valid || control.generation != m_registry->generation() || (id != 0 && (!m_registry->source(id) || m_active != 0))) {
        m_valid  = false;
        m_active = 0;
    } else if (m_valid) {
        if (m_expected && (m_event_index >= m_expected->size() || (*m_expected)[m_event_index] != id)) {
            m_valid  = false;
            m_active = 0;
        } else {
            if (m_expected)
                ++m_event_index;
            m_active = id;
        }
    }
    return true;
}

std::optional<std::uint32_t> DepositionProvenanceDecoder::source_id(std::uint32_t final_logical_filament)
{
    if (!m_valid || m_active == 0)
        return std::nullopt;
    if (m_registry->source(m_active)->logical_filament != final_logical_filament) {
        m_valid  = false;
        m_active = 0;
        return std::nullopt;
    }
    return m_active;
}

void DepositionProvenanceDecoder::finish() noexcept
{
    if (m_active != 0 || (m_expected && m_event_index != m_expected->size()))
        m_valid = false;
    m_active   = 0;
    m_finished = true;
}

DepositionProvenanceCapture::DepositionProvenanceCapture(std::uint64_t generation, std::size_t event_capacity)
    : m_generation(generation), m_capacity(event_capacity)
{
    if (generation == 0 || event_capacity == 0 || event_capacity > DepositionProvenanceDecoder::maximum_events)
        throw std::invalid_argument("Invalid deposition provenance capture budget or generation.");
}

bool DepositionProvenanceCapture::consume(std::string_view line)
{
    if (m_finished)
        m_valid = false;
    const auto control = control_line(line);
    if (!control.reserved)
        return false;
    if (!control.valid || control.generation != m_generation || (control.id != 0 && m_active != 0) || m_events.size() == m_capacity) {
        m_valid = false;
    }
    if (m_valid) {
        m_events.push_back(control.id);
        m_active = control.id;
    } else {
        m_active = 0;
    }
    return true;
}

bool DepositionProvenanceCapture::finish(std::shared_ptr<const DepositionProvenance> registry,
                                         const std::vector<std::uint32_t>& expected_events)
{
    if (m_finished || m_active != 0 || !registry || !registry->sealed() || registry->generation() != m_generation)
        m_valid = false;
    m_finished = true;
    m_active   = 0;
    if (!m_valid)
        return false;
    // Unknown provisional IDs are rejected before constructing token strings.
    for (const auto id : m_events)
        if (id != 0 && !registry->source(id))
            return m_valid = false;
    m_valid = false;
    DepositionProvenanceDecoder decoder(registry, expected_events);
    for (const auto id : m_events)
        decoder.consume(id == 0 ? registry->clear() : registry->begin(id));
    decoder.finish();
    m_valid = decoder.transcript_verified();
    return m_valid;
}
} // namespace Slic3r
