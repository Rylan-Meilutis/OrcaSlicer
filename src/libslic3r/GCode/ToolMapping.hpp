#pragma once

#include "../PrintConfig.hpp"
#include "../Config.hpp"
#include "../PlaceholderParser.hpp"
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <string>
#include <vector>

namespace Slic3r {

struct GCodeTemplateEnvironment {
    DynamicConfig parser_config, external_config;
};

// A self-contained expansion snapshot. It owns every parser input; no pointer
// into a Print, GCode writer or preset survives slicing. Only custom sections
// are replayed: the surrounding generated motion is copied verbatim.
struct GCodeTemplateSnapshot {
    std::string name, text, rendered;
    unsigned current_tool = 0;
    std::shared_ptr<const GCodeTemplateEnvironment> environment;
    DynamicConfig parser_changes, overrides, outputs_before, outputs_after;
    std::vector<std::string> removed_parser_keys;
    DynamicConfig globals_before, globals_after;
    std::mt19937 rng_before, rng_after;
    std::set<std::string> inactive_before, inactive_after;
    std::optional<std::size_t> legacy_machine_mode_index;

    static GCodeTemplateSnapshot capture(const std::string &name, const std::string &text,
        unsigned current_tool, const PlaceholderParser &parser, const DynamicConfig *overrides,
        const DynamicConfig &outputs, const PlaceholderParser::ContextData &context,
        std::shared_ptr<const GCodeTemplateEnvironment> environment = {});
    void finish(const std::string &output, const DynamicConfig &outputs,
                const PlaceholderParser::ContextData &context);
};

struct GCodeDispatchData {
    std::vector<GCodeTemplateSnapshot> templates;
    std::size_t memory_bytes = 0;
    bool valid = true;
    bool sealed = false;
    std::uint64_t source_bytes = 0;
    std::uint64_t source_fingerprint = 0;
};

// Seal only the generator-owned, final postprocessed stream. This also protects
// unmarked native motion from edits made before dispatch preflight begins.
void seal_gcode_dispatch(std::istream &source, GCodeDispatchData &data,
                        const std::function<void()> &checkpoint = {});

struct GCodeDispatchEdit {
    std::uint64_t offset = 0;
    std::string expected, replacement;
};
struct GCodeDispatchPlan {
    std::vector<GCodeDispatchEdit> edits;
    std::uint64_t source_bytes = 0;
    std::uint64_t source_fingerprint = 14695981039346656037ULL;
};

// Preflight is read-only and bounded: preserve all motion, replay only captured
// macro spans, relocate explicitly recognised native tool operands, and replace
// metadata. Unknown native commands, modified macro spans or missing captures
// require normal regeneration. No output file is allocated during preflight.
GCodeDispatchPlan prepare_tool_mapped_gcode(std::istream &source, const GCodeDispatchData &data,
    const DynamicPrintConfig &target, const std::vector<int> &tools, const std::string &full_config,
    const std::function<void(std::uint64_t)> &progress = {});
void apply_tool_mapped_gcode(std::istream &source, std::ostream &target, const GCodeDispatchPlan &plan,
    const std::function<void(std::uint64_t)> &progress = {});

// Replays the original first to certify that all dependencies were captured,
// then evaluates with relocated runtime tool IDs and physical-slot arrays.
// Rejects state changes which would require regenerating surrounding motion.
std::string replay_tool_mapped_template(const GCodeTemplateSnapshot &snapshot,
    const DynamicPrintConfig &target, const std::vector<int> &tools);

// A binding is emitted by the generator, never inferred by searching for "T".
// It can refer to a zero-based firmware tool or a one-based metadata slot.
// Offsets refer to the final source file, after timing/cooling post-processing.
struct GCodeToolBinding {
    std::uint64_t offset;
    std::size_t logical_tool;
    unsigned index_base = 0;
};

// Copies the motion stream and changes only certified numeric bindings.
// Rejects stale, overlapping, unsorted or out-of-range bindings. The caller
// writes a private temporary file and publishes it only after success. Progress
// reports actual source bytes consumed; exceptions from the callback cancel.
void remap_gcode_bindings(std::istream &source, std::ostream &target,
                         const std::vector<GCodeToolBinding> &bindings,
                         const std::vector<int> &tools,
                         const std::function<void(std::uint64_t)> &progress = {});

// Reuse of a sliced job is safe for an unchanged assignment only when the
// complete effective configuration also matches. Comparing nozzle diameters
// alone misses custom scripts, offsets, cooling, flow, purge and feature roles.
// Callers supply configurations applied to the same plate/model snapshot.
bool can_reuse_identity_tool_mapping(const DynamicPrintConfig &source,
                                    const DynamicPrintConfig &target,
                                    const std::vector<int> &tools);

// Compare the used material rows, not padding for unloaded physical slots.
// All unrecognised differences fail closed. A successful comparison establishes
// configuration equivalence only; custom scripts and generated tool references
// still need certified replay before a changed assignment can reuse motion.
bool equivalent_tool_mapping_config(const DynamicPrintConfig &source,
                                    const DynamicPrintConfig &target,
                                    const std::vector<int> &tools);

} // namespace Slic3r
