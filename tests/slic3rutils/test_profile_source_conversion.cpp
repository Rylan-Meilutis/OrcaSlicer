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

TEST_CASE("Prusa profile sources replace malformed UTF-8 without terminating synchronization",
          "[ProfileSources][Regression]")
{
    TemporaryDirectory temporary;
    const fs::path input = temporary.path / "input";
    const fs::path output = temporary.path / "output";

    std::string contents =
        "[filament:Imported filament]\n"
        "filament_type = PLA\n"
        "notes = Bad ";
    contents.push_back(static_cast<char>(0xc3));
    contents += "_description\n";
    write_file(input / "malformed.ini", contents);

    const Slic3r::ProfileSourceSyncResult result =
        Slic3r::ProfileSourceManager::convert_prusa_profiles(input.string(), output.string());

    REQUIRE(result.success());
    REQUIRE(result.filaments == 1);
    const json filament = read_json(output / "filament" / "Imported_filament.json");
    CHECK(filament.at("notes").get<std::string>().find("\xef\xbf\xbd") != std::string::npos);
}

TEST_CASE("Imported profile names are isolated by source", "[ProfileSources][Regression]")
{
    TemporaryDirectory temporary;
    const fs::path input = temporary.path / "input";
    const fs::path first_output = temporary.path / "first";
    const fs::path second_output = temporary.path / "second";
    write_file(input / "Vendor" / "profiles.ini",
        "[printer:Shared printer]\n"
        "nozzle_diameter = 0.4\n"
        "bed_shape = 0x0,200x0,200x200,0x200\n"
        "max_print_height = 200\n"
        "default_print_profile = Shared process\n"
        "\n"
        "[print:Shared process]\n"
        "layer_height = 0.2\n"
        "compatible_printers = Shared printer\n");

    REQUIRE(Slic3r::ProfileSourceManager::convert_prusa_profiles(
        input.string(), first_output.string(), "First source").success());
    REQUIRE(Slic3r::ProfileSourceManager::convert_prusa_profiles(
        input.string(), second_output.string(), "Second source").success());

    const json first = read_json(first_output / "machine" / "Shared_printer_First_source_Vendor.json");
    const json second = read_json(second_output / "machine" / "Shared_printer_Second_source_Vendor.json");
    CHECK(first.at("name") == "Shared printer [First source / Vendor]");
    CHECK(second.at("name") == "Shared printer [Second source / Vendor]");
    CHECK(first.at("name") != second.at("name"));
    CHECK(first.at("default_print_profile") == "Shared process [First source / Vendor]");
    const json process = read_json(first_output / "process" / "Shared_process_First_source_Vendor.json");
    CHECK(process.at("compatible_printers") ==
          std::vector<std::string>{"Shared printer [First source / Vendor]"});
}

TEST_CASE("Native profile sources preserve namespaced dependencies", "[ProfileSources][Regression]")
{
    TemporaryDirectory temporary;
    const fs::path input = temporary.path / "input";
    const fs::path output = temporary.path / "output";
    const fs::path vendor = input / "resources" / "profiles" / "Example";

    write_file(vendor / "machine" / "base.json",
        R"({"name":"Shared machine","type":"machine","nozzle_diameter":["0.4"]})");
    write_file(vendor / "machine" / "derived.json",
        R"({"name":"Derived machine","type":"machine","inherits":"Shared machine","default_print_profile":"Shared process"})");
    write_file(vendor / "process" / "process.json",
        R"({"name":"Shared process","type":"process","compatible_printers":["Derived machine"]})");
    write_file(vendor / "filament" / "filament.json",
        R"({"name":"Shared filament","type":"filament","compatible_printers":["Derived machine"]})");

    const Slic3r::ProfileSourceSyncResult result =
        Slic3r::ProfileSourceManager::convert_orca_profiles(
            input.string(), output.string(), "example-source", "Example source");

    REQUIRE(result.success());
    REQUIRE(result.printers == 2);
    REQUIRE(result.processes == 1);
    REQUIRE(result.filaments == 1);

    const json machine = read_json(output / "machine" / "Example__Derived_machine.json");
    CHECK(machine.at("name") == "Derived machine [Example source / Example]");
    CHECK(machine.at("inherits") == "Shared machine [Example source / Example]");
    CHECK(machine.at("default_print_profile") == "Shared process [Example source / Example]");

    const json process = read_json(output / "process" / "Example__Shared_process.json");
    CHECK(process.at("compatible_printers") ==
          std::vector<std::string>{"Derived machine [Example source / Example]"});
}
