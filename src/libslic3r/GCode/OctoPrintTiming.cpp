#include "OctoPrintTiming.hpp"
#include "LearnedTiming.hpp"
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/nowide/fstream.hpp>
#include <boost/uuid/detail/sha1.hpp>
#include <boost/system/error_code.hpp>
#include <nlohmann/json.hpp>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <ios>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>

namespace Slic3r::OctoPrintTiming {
namespace {
constexpr std::size_t maximum_bytes = 1024 * 1024;
std::mutex store_mutex;
bool valid_hash(const std::string &value)
{
    if (value.size() != 40) return false;
    for (const char c : value) if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) return false;
    return true;
}
std::string digest(boost::uuids::detail::sha1 &hash)
{
    boost::uuids::detail::sha1::digest_type words;
    hash.get_digest(words);
    std::ostringstream out;
    for (const auto word : words) out << std::hex << std::setfill('0') << std::setw(8) << word;
    return out.str();
}
nlohmann::json parse(const std::string &body)
{
    if (body.size() > maximum_bytes) throw std::runtime_error("Timing metadata exceeds its size limit.");
    return nlohmann::json::parse(body, [](int depth, nlohmann::json::parse_event_t, nlohmann::json &) {
        if (depth > 16) throw std::runtime_error("Timing metadata exceeds its nesting limit.");
        return true;
    });
}
}

std::string server_key(const std::string &host)
{
    boost::uuids::detail::sha1 hash;
    hash.process_bytes(host.data(), host.size());
    return digest(hash);
}

std::string selected_local_path(const std::string &body)
{
    const auto file = parse(body).at("job").at("file");
    if (file.value("origin", "") != "local") return {};
    const auto path = file.value("path", "");
    if (path.empty() || path.size() > 1024 || path.front() == '/' || path.find('\\') != std::string::npos) return {};
    std::istringstream parts(path);
    std::string part;
    while (std::getline(parts, part, '/'))
        if (part.empty() || part == "." || part == "..") return {};
    for (const unsigned char c : path) if (c < 32 || c == 127) return {};
    return path;
}

Status completed_print(const std::string &body, const std::string &expected_sha1,
                       const std::string &server, CompletedPrint &result)
{
    result = {};
    try {
        if (!valid_hash(expected_sha1) || !valid_hash(server)) return Status::InvalidResponse;
        const auto file = parse(body);
        if (file.at("hash").get<std::string>() != expected_sha1) return Status::HashMismatch;
        if (!file.contains("prints") || !file.at("prints").contains("last")) return Status::NoCompletedJob;
        const auto &last = file.at("prints").at("last");
        if (!last.at("success").is_boolean() || !last.at("success").get<bool>()) return Status::NoCompletedJob;
        const double seconds = last.at("printTime").get<double>();
        const double timestamp = last.at("date").get<double>();
        const double uploaded = file.at("date").get<double>();
        if (!std::isfinite(seconds) || seconds <= 0. || seconds > 365. * 86400. ||
            !std::isfinite(timestamp) || timestamp <= 0. || timestamp > 1e12 ||
            !std::isfinite(uploaded) || uploaded <= 0. || uploaded > 1e12) return Status::InvalidResponse;
        if (timestamp < uploaded) return Status::NoCompletedJob;
        result.sha1 = expected_sha1;
        result.active_seconds = seconds;
        std::ostringstream identity;
        identity << server << ':' << expected_sha1 << ':' << std::setprecision(17) << timestamp;
        result.job_id = identity.str();
        return Status::Success;
    } catch (const nlohmann::json::exception &) { return Status::InvalidResponse; }
      catch (const std::runtime_error &) { return Status::InvalidResponse; }
}

std::string file_sha1(const std::string &path, const std::function<void()> &checkpoint)
{
    boost::nowide::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Could not read timing source file.");
    boost::uuids::detail::sha1 hash;
    std::array<char, 65536> buffer;
    std::size_t total = 0;
    while (input) {
        if (checkpoint) checkpoint();
        input.read(buffer.data(), buffer.size());
        const auto count = static_cast<std::size_t>(input.gcount());
        if (count > 256 * 1024 * 1024 - total) throw std::runtime_error("Timing source exceeds its size limit.");
        total += count;
        hash.process_bytes(buffer.data(), count);
    }
    if (!input.eof()) throw std::runtime_error("Could not read timing source file.");
    return digest(hash);
}

LearnedTimingModel load(const std::string &path)
{
    const std::lock_guard<std::mutex> lock(store_mutex);
    LearnedTimingModel model;
    if (!boost::filesystem::exists(path)) return model;
    if (boost::filesystem::file_size(path) > maximum_bytes) throw std::runtime_error("Timing history exceeds its size limit.");
    boost::nowide::ifstream input(path, std::ios::binary);
    std::string body(maximum_bytes + 1, '\0');
    input.read(body.data(), body.size());
    body.resize(static_cast<std::size_t>(input.gcount()));
    const auto document = parse(body);
    const auto &records = document.at("observations");
    if (document.at("version") != 1 || !records.is_array() || records.size() > LearnedTimingModel::maximum_observations)
        throw std::runtime_error("Invalid timing history.");
    for (const auto &record : records) {
        LearnedTimingObservation observation;
        observation.job_id = record.at("job_id");
        observation.active_seconds = record.at("active_seconds");
        observation.completed = record.at("completed");
        observation.pause_excluded = record.at("pause_excluded");
        const auto &f = record.at("features");
        observation.features = {f.at("compatibility_key"), f.at("gcode_hash"), f.at("baseline_seconds"),
            f.at("motion_seconds"), f.at("tool_changes"), f.at("heater_waits"), f.at("cleaner_calls"),
            f.at("calibration_calls")};
        if (!model.add_observation(observation)) throw std::runtime_error("Invalid timing history record.");
    }
    return model;
}

void save(const std::string &path, const LearnedTimingModel &model)
{
    const std::lock_guard<std::mutex> lock(store_mutex);
    nlohmann::json records = nlohmann::json::array();
    for (const auto &o : model.observations()) {
        const auto &f = o.features;
        records.push_back({{"job_id", o.job_id}, {"active_seconds", o.active_seconds},
            {"completed", o.completed}, {"pause_excluded", o.pause_excluded},
            {"features", {{"compatibility_key", f.compatibility_key}, {"gcode_hash", f.gcode_hash},
                {"baseline_seconds", f.baseline_seconds}, {"motion_seconds", f.motion_seconds},
                {"tool_changes", f.tool_changes}, {"heater_waits", f.heater_waits}, {"cleaner_calls", f.cleaner_calls},
                {"calibration_calls", f.calibration_calls}}}});
    }
    const auto body = nlohmann::json{{"version", 1}, {"observations", records}}.dump();
    if (body.size() > maximum_bytes) throw std::runtime_error("Timing history exceeds its size limit.");
    // General config helpers may remove or overwrite the destination after
    // rename failure. Learning history must preserve the previous complete
    // record set instead, so use a single strict same-directory rename.
    const auto temporary = boost::filesystem::unique_path(path + ".%%%%-%%%%.tmp").string();
    try {
        boost::nowide::ofstream output(temporary, std::ios::binary | std::ios::trunc);
        output.write(body.data(), body.size());
        output.close();
        if (!output) throw std::runtime_error("Could not save timing history.");
        boost::system::error_code error;
        boost::filesystem::rename(temporary, path, error);
        if (error) throw std::runtime_error("Could not replace timing history.");
    } catch (...) {
        boost::system::error_code ignored;
        boost::filesystem::remove(temporary, ignored);
        throw;
    }
}
} // namespace Slic3r::OctoPrintTiming
