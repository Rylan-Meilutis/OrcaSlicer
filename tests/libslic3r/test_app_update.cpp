#include <catch2/catch_all.hpp>
#include "libslic3r/AppUpdate.hpp"

using namespace Slic3r;

TEST_CASE("Release versions preserve dotted RME revisions", "[AppUpdate]")
{
    const auto older = parse_app_release_version("v2.5.0-rme.2");
    const auto newer = parse_app_release_version("v2.5.0-rme.10");
    REQUIRE(older.has_value());
    REQUIRE(newer.has_value());
    CHECK(*older < *newer);
    CHECK(newer->to_string_sf() == "2.5.0-rme.10");
    CHECK_FALSE(app_release_is_skipped("2.5.0-rme.10", "2.5.0-rme.2"));
    CHECK(app_release_is_skipped("2.5.0-rme.2", "2.5.0-rme.2"));
    CHECK_FALSE(app_release_is_skipped("2.5.0-rme.2", "invalid"));
}

TEST_CASE("Release version parsing rejects partial tags", "[AppUpdate]")
{
    for (const char *tag : {"", "release-2.5.0", "2.5.0 garbage", "2.5", "2.5.0-rme."}) {
        CAPTURE(tag);
        CHECK_FALSE(parse_app_release_version(tag).has_value());
    }
    REQUIRE(parse_app_release_version("2.5.0-rme.3+build.42").has_value());
    CHECK(parse_app_release_version("2.5.0").has_value());
}
