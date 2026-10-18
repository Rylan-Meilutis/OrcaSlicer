#pragma once

#include "GCode/SpoolManagerMetadata.hpp"
#include <boost/uuid/detail/sha1.hpp>
#include <iomanip>
#include <sstream>

namespace Slic3r::OctoPrintSpoolCache {

// Scope snapshots to the server (including its reverse-proxy path) and provider.
// Do not persist connection credentials in the cache key.
inline std::string key(const std::string &host, const std::string &endpoint)
{
    const std::string identity = nlohmann::json::array({host, endpoint}).dump();
    boost::uuids::detail::sha1 hash;
    hash.process_bytes(identity.data(), identity.size());
    boost::uuids::detail::sha1::digest_type digest;
    hash.get_digest(digest);
    std::ostringstream result;
    for (auto word : digest)
        result << std::hex << std::setfill('0') << std::setw(8) << word;
    return result.str();
}

inline std::string encode(const std::vector<SpoolManagerMetadata::Filament> &slots,
                          const std::string &updated_at)
{
    return nlohmann::json{{"version", 1}, {"updated_at", updated_at}, {"slots", slots}}.dump();
}

inline bool decode(const std::string &data, std::vector<SpoolManagerMetadata::Filament> &slots,
                   std::string &updated_at)
{
    slots.clear();
    updated_at.clear();
    try {
        const auto json = nlohmann::json::parse(data);
        if (json.at("version").get<int>() != 1 || !json.at("slots").is_array())
            return false;
        // from_json permits absent fields for live providers; a cached slot must
        // be an object, never a null or primitive silently interpreted as empty.
        for (const auto &slot : json.at("slots"))
            if (!slot.is_object())
                return false;
        auto decoded = json.at("slots").get<std::vector<SpoolManagerMetadata::Filament>>();
        auto timestamp = json.at("updated_at").get<std::string>();
        slots = std::move(decoded);
        updated_at = std::move(timestamp);
        return true;
    } catch (const nlohmann::json::exception &) {
        return false;
    }
}

} // namespace Slic3r::OctoPrintSpoolCache
