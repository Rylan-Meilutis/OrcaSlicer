#pragma once

#include "Semver.hpp"
#include <regex>

namespace Slic3r {

// Release tags use SemVer, including dotted RME revisions. Never extract a
// version substring: a partial match could advertise an unrelated build.
inline boost::optional<Semver> parse_app_release_version(std::string tag)
{
    if (!tag.empty() && tag.front() == 'v')
        tag.erase(0, 1);
    static const std::regex syntax(
        R"(^[0-9]+\.[0-9]+\.[0-9]+(-[0-9A-Za-z-]+(\.[0-9A-Za-z-]+)*)?(\+[0-9A-Za-z-]+(\.[0-9A-Za-z-]+)*)?$)");
    if (!std::regex_match(tag, syntax))
        return boost::none;
    return Semver::parse(tag);
}

inline bool app_release_is_skipped(const std::string &candidate, const std::string &skipped)
{
    const auto version = parse_app_release_version(candidate);
    const auto skip = parse_app_release_version(skipped);
    return version && skip && *version <= *skip;
}

} // namespace Slic3r
