#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <algorithm>
#include <boost/filesystem/path.hpp>
#include <iterator>
#include "libslic3r/Config.hpp"
#include <string>
#include <boost/filesystem.hpp>
#include <fstream>
#include <sstream>
#include <nlohmann/json.hpp>
#include "libslic3r/PlaceholderParser.hpp"
#include "libslic3r/GCode/RMEOutputPlan.hpp"

using namespace Slic3r;

TEST_CASE("INDX startup homes before vent motion and shutdown parks before closing it", "[IndxProfile][INDX]")
{
    const auto root = boost::filesystem::path(PROFILES_DIR);
    std::ifstream input((root / "Prusa" / "machine" / "Prusa CORE One INDX RME 0.4 nozzle.json").string());
    REQUIRE(input.good());
    const auto profile = nlohmann::json::parse(input);
    const std::string start = profile.at("machine_start_gcode");
    const auto home = start.find("\nG28 XY\n");
    const auto vent = start.find("\nM870 ");
    const auto pickup = start.find("T0 S1 L2 D0");
    REQUIRE(home != std::string::npos);
    REQUIRE(vent != std::string::npos);
    REQUIRE(pickup != std::string::npos);
    CHECK(home < vent);
    CHECK(vent < pickup);
    CHECK(start.find("\nG28 XY\n", home + 1) == std::string::npos);
    const std::string end = profile.at("machine_end_gcode");
    const auto park = end.find("P0 S1");
    const auto close = end.find("M870 C");
    REQUIRE(park != std::string::npos);
    REQUIRE(close != std::string::npos);
    CHECK(park < close);
    CHECK(end.find("P0 S1", park + 1) == std::string::npos);
    CHECK(end.find("M400", park) < close);
    std::ifstream end_file((root.parent_path() / "gcode" / "indx_rme_end.gcode").string());
    REQUIRE(end_file.good());
    CHECK(std::string(std::istreambuf_iterator<char>(end_file), {}) == end);
}

TEST_CASE("RME PA plans preserve initial temperature and both deposition extremes", "[IndxProfile][RMEOutput]")
{
    RMEOutputPlan plan;
    plan.temperature(0, 3, "HIPS", 250);
    plan.temperature(0, 3, "HIPS", 240);
    plan.temperature(0, 3, "HIPS", 250);
    CHECK(plan.pa_command() == "M976 A 3:0:HIPS:250/240\n");
    plan.temperature(1, 5, "ASA", 260);
    CHECK(plan.pa_command() == "M976 A 3:0:HIPS:250/240,5:1:ASA:260\n");
    for (int temp = 230; temp <= 260; ++temp) plan.temperature(0, 3, "HIPS", temp);
    const auto command = plan.pa_command();
    CHECK(command.find("3:0:HIPS:250/230/260") != std::string::npos);
    CHECK(std::count(command.begin(), command.end(), '/') == 7);
    plan.temperature(0, 3, "HIPS", 0);
    CHECK(plan.pa_command().empty());
}

TEST_CASE("RME mesh bounds include exact arc extrema and bead width", "[IndxProfile][RMEOutput]")
{
    RMEOutputPlan plan;
    CHECK(plan.mesh_command().empty());
    plan.arc(10., 0., -10., 0., 0., 0., false, 0.4);
    CHECK_THAT(plan.minimum[0], Catch::Matchers::WithinAbs(-10.2, 1e-6));
    CHECK_THAT(plan.maximum[0], Catch::Matchers::WithinAbs(10.2, 1e-6));
    CHECK_THAT(plan.maximum[1], Catch::Matchers::WithinAbs(10.2, 1e-6));
    CHECK_THAT(plan.minimum[1], Catch::Matchers::WithinAbs(-0.2, 1e-6));
    plan.arc(10., 0., 10., 0., 0., 0., true, 0.4);
    CHECK_THAT(plan.minimum[1], Catch::Matchers::WithinAbs(-10.2, 1e-6));
    plan.point(100., 50., 1.);
    CHECK_THAT(plan.maximum[0], Catch::Matchers::WithinAbs(100.5, 1e-6));
    CHECK_THAT(plan.maximum[1], Catch::Matchers::WithinAbs(50.5, 1e-6));
}

TEST_CASE("Bundled INDX mesh bounds stay local for small corner prints", "[IndxProfile]")
{
    const auto path = boost::filesystem::path(PROFILES_DIR) / "Prusa" / "machine" /
                      "Prusa CORE One INDX RME 0.4 nozzle.json";
    std::ifstream input(path.string());
    REQUIRE(input.good());
    const auto profile = nlohmann::json::parse(input);
    CHECK_FALSE(profile.contains("printhost_apikey"));
    CHECK_FALSE(profile.contains("print_host"));
    CHECK(profile.at("nozzle_diameter").size() == 8);
    CHECK(profile.at("bed_mesh_min") == "0,0");
    CHECK(profile.at("bed_mesh_max") == "248,205");
    const std::string start = profile.at("machine_start_gcode");
    std::ifstream template_file((boost::filesystem::path(PROFILES_DIR).parent_path() / "gcode" / "indx_rme_start.gcode").string());
    REQUIRE(template_file.good());
    const std::string bundled_template((std::istreambuf_iterator<char>(template_file)), std::istreambuf_iterator<char>());
    CHECK(start == bundled_template);
    CHECK(start.find("; RME_DEPOSITION_BOUNDS=1") != std::string::npos);
    CHECK(start.find("; RME_TEMPERATURE_ADAPTIVE_PA=0") != std::string::npos);
    // Plain P1 allows the RME host's file-derived bounds to take precedence.
    CHECK(start.find("G29 P1 X") == std::string::npos);
    const auto probe = start.find("\nG29 P1 ");
    REQUIRE(probe != std::string::npos);
    const auto bounds = start.rfind("\nM555 ", probe);
    REQUIRE(bounds != std::string::npos);
    CHECK(bounds > start.find("\nG427 "));
    const auto command = start.substr(bounds + 1, probe - bounds - 1);
    // These are the clamped, margin-expanded slicer bounds for a 20mm corner part.
    const double x = GENERATE(0., 227.);
    const double y = GENERATE(0., 184.);
    PlaceholderParser parser;
    parser.set("adaptive_bed_mesh_min", new ConfigOptionFloats({x, y}));
    parser.set("adaptive_bed_mesh_max", new ConfigOptionFloats({x + 21., y + 21.}));
    std::ostringstream expected;
    expected << "M555 X" << x << " Y" << y << " W21 H21";
    CHECK(parser.process(command) == expected.str());
}
