#pragma once

#include <cstddef>
#include <functional>
#include <istream>
#include <optional>
#include <string>
#include <vector>

namespace Slic3r {

struct LearnedTimingFeatures {
    // Caller includes host, firmware, macro and timing-relevant machine/material
    // configuration identity. Unknown identities must not be merged.
    std::string compatibility_key;
    std::string gcode_hash;
    double baseline_seconds = 0.;
    double motion_seconds = 0.; // Zero when unavailable, never substituted with baseline.
    unsigned int tool_changes = 0;
    unsigned int heater_waits = 0;
    unsigned int cleaner_calls = 0;
    unsigned int calibration_calls = 0;
};

// Context identifies the host and firmware/setup revision. A user-declared
// setup revision is valid provenance, but is not automatically detected firmware.
// Empty/unknown context must abstain or use exact-file-only identity.
std::optional<LearnedTimingFeatures> extract_learned_timing_features(
    std::istream &input, const std::string &compatibility_context, const std::string &gcode_hash,
    const std::function<void()> &checkpoint = {});

struct LearnedTimingObservation {
    std::string job_id; // Distinct completed execution, not just the file identity.
    LearnedTimingFeatures features;
    double active_seconds = 0.; // Includes heating; excludes user pause duration.
    bool completed = false;
    bool pause_excluded = false;
};

struct LearnedTimingEstimate {
    double seconds = 0.;
    // Observed residual range, NOT a statistical confidence/prediction interval.
    double lower_seconds = 0.;
    double upper_seconds = 0.;
    std::size_t sample_count = 0;
};

class LearnedTimingModel {
public:
    static constexpr std::size_t maximum_observations = 128;
    static constexpr std::size_t minimum_samples = 3;
    bool add_observation(const LearnedTimingObservation &observation);
    std::optional<LearnedTimingEstimate> estimate(const LearnedTimingFeatures &features) const;
    const std::vector<LearnedTimingObservation> &observations() const { return m_observations; }
private:
    std::vector<LearnedTimingObservation> m_observations;
};

} // namespace Slic3r
