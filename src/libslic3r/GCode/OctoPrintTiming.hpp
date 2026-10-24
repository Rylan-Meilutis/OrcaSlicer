#pragma once

#include "LearnedTiming.hpp"
#include <functional>
#include <cstddef>
#include <optional>
#include <string>

namespace Slic3r::OctoPrintTiming {

enum class Status { Success, Disabled, NetworkError, InvalidResponse, NoCompletedJob, HashMismatch,
                    Cancelled, InvalidSource, StorageError };
struct CompletedPrint {
    std::string sha1;
    std::string job_id;
    double active_seconds = 0.;
};
struct RefreshResult {
    Status status = Status::InvalidSource;
    std::optional<LearnedTimingEstimate> estimate;
    std::size_t accepted_samples = 0;
};

// Pure bounded metadata parsing; no printer access or implicit learning.
std::string selected_local_path(const std::string &job_json);
Status completed_print(const std::string &file_json, const std::string &expected_sha1,
                       const std::string &server_key, CompletedPrint &result);
std::string file_sha1(const std::string &path, const std::function<void()> &checkpoint = {});
std::string server_key(const std::string &host);

// At most 128 observations / 1 MiB. Credentials and printer URLs are never
// persisted. Writes use a private sibling and atomic rename; no large files.
LearnedTimingModel load(const std::string &path);
void save(const std::string &path, const LearnedTimingModel &model);

} // namespace Slic3r::OctoPrintTiming
