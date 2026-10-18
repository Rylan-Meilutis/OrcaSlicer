#include <catch2/catch_all.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include <sstream>
#include <nlohmann/json.hpp>
#include "libslic3r/PlaceholderParser.hpp"

using namespace Slic3r;

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
