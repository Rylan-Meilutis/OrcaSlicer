#include <catch2/catch_test_macros.hpp>

#include "slic3r/Utils/ProfileSourceManager.hpp"

#include <boost/filesystem.hpp>
#include <boost/nowide/fstream.hpp>
#include <nlohmann/json.hpp>

namespace fs = boost::filesystem;
using json = nlohmann::json;

namespace {

struct TemporaryDirectory
{
    TemporaryDirectory()
        : path(fs::temp_directory_path() / fs::unique_path("orca-profile-source-%%%%-%%%%"))
    {
        fs::create_directories(path);
    }

    ~TemporaryDirectory()
    {
        boost::system::error_code error;
        fs::remove_all(path, error);
    }

    fs::path path;
};

void write_file(const fs::path &path, const std::string &contents)
{
    fs::create_directories(path.parent_path());
    boost::nowide::ofstream stream(path.string());
    stream << contents;
}

json read_json(const fs::path &path)
{
    boost::nowide::ifstream stream(path.string());
    return json::parse(stream);
}

} // namespace

TEST_CASE("Prusa profile sources convert bundles and standalone presets", "[ProfileSources][Regression]")
{
    TemporaryDirectory temporary;
    const fs::path input = temporary.path / "input";
    const fs::path output = temporary.path / "output";

    write_file(input / "vendor" / "1.2.3.ini",
        "[print:*base*]\n"
        "perimeters = 3\n"
        "fill_pattern = rectilinear\n"
        "perimeter_speed = 60\n"
        "\n"
        "[print:Converted Process]\n"
        "inherits = *base*\n"
        "external_perimeter_speed = 50%\n"
        "complete_objects = 1\n"
        "\n"
        "[filament:Converted PET]\n"
        "filament_type = PET\n"
        "filament_max_volumetric_speed = 0\n"
        "bed_temperature = 80\n");
    write_file(input / "printer" / "Standalone Printer.ini",
        "nozzle_diameter = 0.4\n"
        "bed_shape = 0x0,250x0,250x210,0x210\n"
        "gcode_flavor = reprap\n"
        "max_print_height = 220\n");

    const Slic3r::ProfileSourceSyncResult result =
        Slic3r::ProfileSourceManager::convert_prusa_profiles(input.string(), output.string());

    REQUIRE(result.success());
    REQUIRE(result.processes == 1);
    REQUIRE(result.filaments == 1);
    REQUIRE(result.printers == 1);

    const json process = read_json(output / "process" / "Converted_Process.json");
    CHECK(process.at("wall_loops") == "3");
    CHECK(process.at("sparse_infill_pattern") == "zig-zag");
    CHECK(process.at("outer_wall_speed") == "30");
    CHECK(process.at("print_sequence") == "by object");

    const json filament = read_json(output / "filament" / "Converted_PET.json");
    CHECK(filament.at("filament_type") == std::vector<std::string>{"PETG"});
    CHECK(filament.at("filament_max_volumetric_speed") == std::vector<std::string>{"10"});
    CHECK(filament.at("hot_plate_temp") == std::vector<std::string>{"80"});
    CHECK(filament.at("textured_plate_temp") == std::vector<std::string>{"80"});

    const json printer = read_json(output / "machine" / "Standalone_Printer.json");
    CHECK(printer.at("gcode_flavor") == "reprapfirmware");
    CHECK(printer.at("printable_height") == "220");
    CHECK(printer.at("printable_area").size() == 4);
}
