#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include "libslic3r/GCode/OctoPrintTiming.hpp"
#include "libslic3r/GCode/LearnedTiming.hpp"
#include "libslic3r/Utils.hpp"
#include "test_utils.hpp"
#include <nlohmann/json.hpp>
#include <string>
#include <stdexcept>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

using namespace Slic3r;

TEST_CASE("OctoPrint timing accepts only exact-file successful last attempts", "[OctoPrintTiming]")
{
    const std::string hash(40, 'a'), server(40, 'b');
    nlohmann::json metadata = {{"hash", hash}, {"date", 1699999999}, {"prints", {{"last", {
        {"success", true}, {"date", 1700000000}, {"printTime", 7200}}}}}};
    OctoPrintTiming::CompletedPrint result;
    CHECK(OctoPrintTiming::completed_print(metadata.dump(), hash, server, result) == OctoPrintTiming::Status::Success);
    CHECK_THAT(result.active_seconds, Catch::Matchers::WithinAbs(7200., 1e-9));
    CHECK(result.sha1 == hash);
    const auto job = result.job_id;
    CHECK(OctoPrintTiming::completed_print(metadata.dump(), std::string(40, 'c'), server, result) == OctoPrintTiming::Status::HashMismatch);
    CHECK(result.job_id.empty());
    metadata["prints"]["last"]["success"] = false;
    CHECK(OctoPrintTiming::completed_print(metadata.dump(), hash, server, result) == OctoPrintTiming::Status::NoCompletedJob);
    metadata["prints"]["last"]["success"] = true;
    metadata["prints"]["last"]["date"] = 1700000001;
    CHECK(OctoPrintTiming::completed_print(metadata.dump(), hash, server, result) == OctoPrintTiming::Status::Success);
    CHECK(result.job_id != job);
    metadata["date"] = 1700000002;
    CHECK(OctoPrintTiming::completed_print(metadata.dump(), hash, server, result) == OctoPrintTiming::Status::NoCompletedJob);
    metadata["date"] = 1699999999;
    metadata["prints"]["last"]["printTime"] = -1;
    CHECK(OctoPrintTiming::completed_print(metadata.dump(), hash, server, result) == OctoPrintTiming::Status::InvalidResponse);
    CHECK(OctoPrintTiming::completed_print(std::string(1024 * 1024 + 1, ' '), hash, server, result) == OctoPrintTiming::Status::InvalidResponse);
}

TEST_CASE("OctoPrint timing paths cannot escape local file metadata", "[OctoPrintTiming]")
{
    nlohmann::json job = {{"job", {{"file", {{"origin", "local"}, {"path", "folder/case.gcode"}}}}}};
    CHECK(OctoPrintTiming::selected_local_path(job.dump()) == "folder/case.gcode");
    for (const auto *path : {"../case", "/case", "folder/../case", "folder\\case", "folder/./case"}) {
        job["job"]["file"]["path"] = path;
        CHECK(OctoPrintTiming::selected_local_path(job.dump()).empty());
    }
    job["job"]["file"] = {{"origin", "sdcard"}, {"path", "case.gcode"}};
    CHECK(OctoPrintTiming::selected_local_path(job.dump()).empty());
}

TEST_CASE("Learned timing history persists bounded observations without credentials", "[OctoPrintTiming]")
{
    ScopedTemporaryFile file(".json");
    LearnedTimingModel model;
    LearnedTimingFeatures features{"hashed-context", std::string(40, 'a'), 100., 0., 2, 3, 4, 5};
    REQUIRE(model.add_observation({"completed-job", features, 120., true, true}));
    OctoPrintTiming::save(file.string(), model);
    const auto restored = OctoPrintTiming::load(file.string());
    REQUIRE(restored.observations().size() == 1);
    CHECK(restored.observations().front().features.calibration_calls == 5);
    REQUIRE_FALSE(write_file_atomically(file.string(), "abc"));
    CHECK(OctoPrintTiming::file_sha1(file.string()) == "a9993e364706816aba3e25717850c26c9cd0d89d");
    CHECK_THROWS_AS(OctoPrintTiming::file_sha1(file.string(), [] { throw std::runtime_error("cancelled"); }), std::runtime_error);
    CHECK_THROWS(OctoPrintTiming::load(file.string()));
}

TEST_CASE("Failed timing history replacement preserves the existing destination", "[OctoPrintTiming]")
{
    ScopedTemporaryDir destination;
    const auto existing = (destination.path() / "existing-record").string();
    REQUIRE_FALSE(write_file_atomically(existing, "abc"));
    LearnedTimingModel model;
    CHECK_THROWS(OctoPrintTiming::save(destination.string(), model));
    CHECK(boost::filesystem::is_directory(destination.path()));
    CHECK(OctoPrintTiming::file_sha1(existing) == "a9993e364706816aba3e25717850c26c9cd0d89d");
}
