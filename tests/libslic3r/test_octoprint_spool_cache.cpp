#include <catch2/catch_all.hpp>
#include "libslic3r/OctoPrintSpoolCache.hpp"
#include "libslic3r/AppConfig.hpp"
#include "libslic3r/Thread.hpp"
#include "test_utils.hpp"

using namespace Slic3r;

TEST_CASE("OctoPrint loaded spools survive an application configuration reload", "[OctoPrintSpoolCache]")
{
    ScopedTemporaryDir directory("orca-spool-cache");
    struct RestoreDataDir {
        std::string previous = data_dir();
        ~RestoreDataDir() { set_data_dir(previous); }
    } restore;
    set_data_dir(directory.string());
    save_main_thread_id();
    const auto key = OctoPrintSpoolCache::key("https://printer", "");
    std::vector<SpoolManagerMetadata::Filament> slots(3);
    slots[2].spool_id = "123";
    {
        AppConfig config;
        config.set("octoprint_loaded_spools", key, OctoPrintSpoolCache::encode(slots, "saved"));
        config.save();
    }
    AppConfig reloaded;
    REQUIRE(reloaded.load().empty());
    std::vector<SpoolManagerMetadata::Filament> restored;
    std::string timestamp;
    REQUIRE(OctoPrintSpoolCache::decode(reloaded.get("octoprint_loaded_spools", key), restored, timestamp));
    REQUIRE(restored.size() == 3);
    CHECK(restored[0].spool_id.empty());
    CHECK(restored[2].spool_id == "123");
    CHECK(timestamp == "saved");
}

TEST_CASE("OctoPrint refresh distinguishes unloading every tool from a missing assignment response", "[OctoPrintSpoolCache]")
{
    std::vector<SpoolManagerMetadata::Filament> slots;
    std::string error;
    REQUIRE(SpoolManagerMetadata::parse_selected_spools(
        R"({"selectedSpools":[null,null]})", slots, error, {}, true));
    REQUIRE(slots.size() == 2);
    CHECK(slots[0].spool_id.empty());
    CHECK(slots[1].spool_id.empty());
    CHECK_FALSE(SpoolManagerMetadata::parse_selected_spools(
        R"({"selectedSpools":[null,null]})", slots, error));
    CHECK_FALSE(SpoolManagerMetadata::parse_selected_spools(
        R"({"allSpools":[]})", slots, error, {}, true));
}

TEST_CASE("OctoPrint snapshots preserve loaded and unloaded tool indices", "[OctoPrintSpoolCache]")
{
    std::vector<SpoolManagerMetadata::Filament> slots(4);
    slots[2].material = "HIPS";
    slots[2].spool_id = "42";
    slots[2].color = "#ABCDEF";
    slots[2].vendor = "Example";
    slots[2].provider = "Spoolman";
    std::vector<SpoolManagerMetadata::Filament> restored;
    std::string timestamp;
    REQUIRE(OctoPrintSpoolCache::decode(OctoPrintSpoolCache::encode(slots, "2026-10-17 12:00:00"), restored, timestamp));
    REQUIRE(restored.size() == 4);
    CHECK(restored[0].material.empty());
    CHECK(restored[1].spool_id.empty());
    CHECK(restored[3].material.empty());
    CHECK(nlohmann::json(restored) == nlohmann::json(slots));
    CHECK(timestamp == "2026-10-17 12:00:00");
    REQUIRE(OctoPrintSpoolCache::decode(OctoPrintSpoolCache::encode({}, "newer"), restored, timestamp));
    CHECK(restored.empty());
    CHECK(timestamp == "newer");
}

TEST_CASE("OctoPrint snapshots reject corrupt and incompatible data without retaining slots", "[OctoPrintSpoolCache]")
{
    const auto data = GENERATE("", "{", "[]", "null",
        R"({"version":2,"updated_at":"now","slots":[]})",
        R"({"version":1,"updated_at":"now","slots":[null]})",
        R"({"version":1,"updated_at":"now","slots":[{"material":42}]})",
        R"({"version":1,"slots":[{}]})");
    std::vector<SpoolManagerMetadata::Filament> slots(8);
    std::string timestamp = "old";
    CHECK_FALSE(OctoPrintSpoolCache::decode(data, slots, timestamp));
    CHECK(slots.empty());
    CHECK(timestamp.empty());
}

TEST_CASE("OctoPrint snapshot keys isolate servers and provider endpoints", "[OctoPrintSpoolCache]")
{
    const auto key = OctoPrintSpoolCache::key("https://printer/octoprint", "plugin/provider");
    CHECK(key == OctoPrintSpoolCache::key("https://printer/octoprint", "plugin/provider"));
    CHECK(key != OctoPrintSpoolCache::key("https://other/octoprint", "plugin/provider"));
    CHECK(key != OctoPrintSpoolCache::key("https://printer/other", "plugin/provider"));
    CHECK(key != OctoPrintSpoolCache::key("https://printer/octoprint", "plugin/other"));
    CHECK(key.find("printer") == std::string::npos);
}
