#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace Slic3r {
// Nominal design envelope, NOT a measurement or a reconstruction from E.
// Coordinates and dimensions are in mm in the emitted Cartesian command frame.
struct DepositionCommandSection
{
    enum class Transform : std::uint8_t { Ordinary, BrickCourse };
    enum class Purpose : std::uint8_t { Model, AuxiliaryBrim };
    std::array<double, 3> start, end;
    double width = 0., start_height = 0., end_height = 0.;
    std::uint32_t filament = 0;
    std::array<double, 2> arc_center{};
    // Zero denotes a line; +1 is counterclockwise and -1 is clockwise XY.
    int arc_direction = 0;
    // Explicit intended lower surface in machine mm, separate from the flow
    // ratio. Present only where the generator knows that geometric surface.
    // With both heights zero this certifies a collapsed, unresolved envelope,
    // not a volumetric bead; consumers must preserve that distinction.
    std::optional<std::array<double, 2>> nominal_lower_z;
    // Generator-authenticated geometry, independent of the preview role.
    Transform transform = Transform::Ordinary;
    // Auxiliary command geometry is not a fabricated CAD/object source.
    Purpose purpose = Purpose::Model;
};

class DepositionSections
{
public:
    // At most 80 MiB of nominal parent data per registry on supported ABIs;
    // final validation temporarily retains both machine and parsed frames.
    static constexpr std::size_t maximum_commands = 80 * 1024 * 1024 / sizeof(DepositionCommandSection);
    static_assert(sizeof(DepositionCommandSection) * maximum_commands <= 80 * 1024 * 1024,
                  "Deposition parent storage exceeds its bounded workspace");
    explicit DepositionSections(std::uint64_t generation);
    std::uint32_t add(const DepositionCommandSection& section);
    void seal() noexcept { m_sealed = true; }
    bool sealed() const noexcept { return m_sealed; }
    std::uint64_t generation() const noexcept { return m_generation; }
    std::size_t size() const noexcept { return m_sections.size(); }
    const DepositionCommandSection* section(std::uint32_t id) const noexcept;
    std::string marker(std::uint32_t id) const;

private:
    std::uint64_t m_generation;
    bool m_sealed = false;
    std::vector<DepositionCommandSection> m_sections;
};

class DepositionSectionCapture
{
public:
    explicit DepositionSectionCapture(std::uint64_t generation) : m_generation(generation) {}
    bool consume(std::string_view line);
    std::uint32_t active() const noexcept { return m_valid ? m_active : 0; }
    bool finish(const DepositionSections& sections);

private:
    std::uint64_t m_generation;
    std::uint32_t m_active = 0;
    bool m_valid           = true;
    std::vector<std::uint32_t> m_events;
};

// Validate an ordered child chain against one quantized parent command. Children
// may be split by cooling/pressure processors or arc tessellation, but cannot move off its line or circle,
// reverse, skip or repeat a span. Tolerance covers textual XYZ quantization only.
class DepositionSectionChain
{
public:
    explicit DepositionSectionChain(const DepositionCommandSection& parent) : m_parent(parent) {}
    std::optional<std::array<double, 2>> append(const std::array<double, 3>& start, const std::array<double, 3>& end, std::uint32_t filament);
    bool complete() const noexcept { return m_valid && m_fraction >= 1. - 1e-7; }

private:
    DepositionCommandSection m_parent;
    double m_fraction = 0.;
    bool m_valid      = true;
};
} // namespace Slic3r
