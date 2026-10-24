#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include "libslic3r/Config.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/ThermoelasticGCode.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "test_helpers.hpp"
#include <Eigen/Core>
#include <boost/nowide/fstream.hpp>
#include <cstddef>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace Slic3r;
using Catch::Matchers::WithinAbs;

// Optional user-file coverage audit, not a physical validation or an adhesion
// test. No output files are produced and the mesh resource limits stay fixed.
TEST_CASE("Mixed material case coverage remains unavailable around the observed print height", "[.][ThermoelasticGCodeProject]")
{
    const char *path = std::getenv("ORCA_MECHANICAL_GCODE");
    if (!path) SKIP("Set ORCA_MECHANICAL_GCODE to the mixed ASA/TPU case export.");
    auto config = DynamicPrintConfig::full_print_config();
    // Read only the scenario inputs, avoiding unrelated GUI preset migration
    // in the complete configuration loader. Require every consumed field.
    std::set<std::string> required{"filament_type", "chamber_temperature", "nozzle_temperature",
        "curr_bed_type", "hot_plate_temp"};
    boost::nowide::ifstream input(path);
    REQUIRE(input.good());
    bool in_config = false;
    std::size_t source_line = 0, layer_start = 0, layer_end = 0;
    std::string line;
    while (std::getline(input, line)) {
        ++source_line;
        if (line.rfind(";Z:", 0) == 0) {
            const double layer_z = std::stod(line.substr(3));
            if (layer_start == 0 && std::abs(layer_z - 11.4) < 1e-4) layer_start = source_line;
            else if (layer_start != 0 && layer_end == 0 && layer_z > 11.4 + 1e-4) layer_end = source_line;
        }
        if (line == "; CONFIG_BLOCK_START") { in_config = true; continue; }
        if (line == "; CONFIG_BLOCK_END") break;
        if (!in_config || line.rfind("; ", 0) != 0) continue;
        const auto equals = line.find(" = ");
        if (equals == std::string::npos) continue;
        const auto key = line.substr(2, equals - 2);
        if (required.erase(key)) config.set_deserialize_strict(key, line.substr(equals + 3));
    }
    REQUIRE(required.empty());
    REQUIRE(layer_start > 0);
    REQUIRE(layer_end > layer_start);
    GCodeProcessor processor;
    processor.process_file(path);
    auto &gcode = processor.result();
    const auto deposited = [](const GCodeProcessorResult::MoveVertex &move) {
        return move.type == EMoveType::Extrude && move.extrusion_role != erWipeTower &&
            move.extrusion_role != erCustom && move.extrusion_role != erNone;
    };
    std::size_t before = gcode.moves.size(), through = before;
    for (std::size_t i = 1; i < gcode.moves.size(); ++i) {
        if (before == gcode.moves.size() && gcode.moves[i].gcode_id >= layer_start) before = i;
        if (gcode.moves[i].gcode_id >= layer_end) { through = i; break; }
    }
    REQUIRE(before < through);
    REQUIRE(through < gcode.moves.size());
    const auto audit = [&](const char *label) {
        std::map<unsigned, std::size_t> material_counts;
        std::size_t count = 0, sloped = 0, first_sloped_line = 0;
        for (std::size_t i = 1; i < gcode.moves.size(); ++i) {
            const auto &move = gcode.moves[i];
            if (!deposited(move)) continue;
            ++count;
            ++material_counts[move.extruder_id];
            if (std::abs(double(move.position.z()) - double(gcode.moves[i-1].position.z())) > 1e-6) {
                if (sloped++ == 0) first_sloped_line = move.gcode_id;
            }
        }
        REQUIRE(count > Thermoelastic::AnalysisLimits{}.max_segments);
        const auto bounded = Thermoelastic::analyze_configured_gcode(gcode, config);
        CHECK(bounded.state == Thermoelastic::AnalysisState::BudgetExceeded);
        CHECK(bounded.scenarios.empty());
        CHECK(gcode.thermal_stress_mpa.empty());
        // Admit the already parsed segment count only to expose material
        // coverage. This does not increase nodes, raster work or refinement.
        Thermoelastic::AnalysisLimits diagnostic_limits;
        diagnostic_limits.max_segments = count;
        const auto diagnostic = Thermoelastic::analyze_configured_gcode(gcode, config, diagnostic_limits);
        CHECK(diagnostic.state == Thermoelastic::AnalysisState::UnsupportedMaterial);
        CHECK(diagnostic.scenarios.empty());
        CHECK(gcode.thermal_stress_mpa.empty());
        REQUIRE_FALSE(diagnostic.unsupported_filaments.empty());
        std::cout << "Mechanical coverage " << label << ": deposited=" << count
                  << " nonplanar=" << sloped << " first_nonplanar_gcode_id=" << first_sloped_line
                  << " default=BudgetExceeded admitted=UnsupportedMaterial material_counts=";
        for (const auto &[id, n] : material_counts) std::cout << ' ' << id << ':' << n;
        std::cout << " unsupported=";
        for (const auto id : diagnostic.unsupported_filaments) std::cout << ' ' << id;
        std::cout << " (zero-based IDs; no force or detachment prediction)" << std::endl;
    };
    audit("complete export");
    // Chronological prefixes retain every predecessor/travel/event vertex.
    // They bracket the observed layer; its actual execution point is unknown.
    gcode.moves.resize(through);
    audit("through deposition at Z11.4mm");
    gcode.moves.resize(before);
    audit("before deposition at Z11.4mm");
}

TEST_CASE("Configured mechanical scenarios never substitute properties for unsupported filaments", "[ThermoelasticGCode]")
{
    auto config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({{"filament_type", "TPU"}, {"chamber_temperature", "50"},
        {"nozzle_temperature", "265"}, {"curr_bed_type", "Textured PEI Plate"},
        {"textured_plate_temp", "110"}});
    GCodeProcessorResult gcode;
    gcode.moves.resize(2);
    gcode.moves[0].position = Eigen::Vector3f(0.f, 0.f, .2f);
    gcode.moves[1] = gcode.moves[0];
    gcode.moves[1].position.x() = 1.f;
    gcode.moves[1].type = EMoveType::Extrude;
    gcode.moves[1].extrusion_role = erPerimeter;
    gcode.moves[1].width = .4f;
    gcode.moves[1].height = .2f;
    gcode.thermal_stress_mpa = {123.f, 123.f};
    auto result = Thermoelastic::analyze_configured_gcode(gcode, config);
    CHECK(result.state == Thermoelastic::AnalysisState::UnsupportedMaterial);
    REQUIRE(result.unsupported_filaments.size() == 1);
    CHECK(result.unsupported_filaments.front() == 0);
    CHECK(gcode.thermal_stress_mpa.empty());
    config.set_deserialize_strict({{"filament_type", "ASA"}, {"chamber_temperature", "0"}});
    result = Thermoelastic::analyze_configured_gcode(gcode, config);
    CHECK(result.state == Thermoelastic::AnalysisState::InvalidThermalInputs);
    config.set_deserialize_strict({{"chamber_temperature", "50"}});
    Thermoelastic::AnalysisLimits limits;
    limits.max_nodes = 4;
    result = Thermoelastic::analyze_configured_gcode(gcode, config, limits);
    CHECK(result.state == Thermoelastic::AnalysisState::BudgetExceeded);
    CHECK(result.scenarios.empty());
    CHECK(gcode.thermal_stress_mpa.empty());
    // A glassy modulus above Tg is outside the model, not a zero-risk result.
    config.set_deserialize_strict({{"chamber_temperature", "110"}});
    result = Thermoelastic::analyze_configured_gcode(gcode, config);
    CHECK(result.state == Thermoelastic::AnalysisState::InvalidThermalInputs);
    CHECK(gcode.thermal_stress_mpa.empty());
    CHECK_THROWS_AS(Thermoelastic::analyze_configured_gcode(gcode, config, {}, {},
        [] { throw std::runtime_error("cancelled"); }), std::runtime_error);
}

TEST_CASE("Equivalent ASA colors share cells only with identical thermal scenario inputs", "[ThermoelasticGCode]")
{
    auto config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({{"filament_type", "ASA;ASA"}, {"chamber_temperature", "50,50"},
        {"nozzle_temperature", "265,265"}, {"curr_bed_type", "Textured PEI Plate"},
        {"textured_plate_temp", "110,110"}});
    GCodeProcessorResult gcode;
    gcode.moves.resize(4);
    gcode.moves[0].position = Eigen::Vector3f(0.f, 0.f, .2f);
    gcode.moves[1] = gcode.moves[0];
    gcode.moves[1].position.x() = 1.f;
    gcode.moves[1].type = EMoveType::Extrude;
    gcode.moves[1].extrusion_role = erPerimeter;
    gcode.moves[1].width = .4f;
    gcode.moves[1].height = .2f;
    gcode.moves[2] = gcode.moves[0];
    gcode.moves[2].type = EMoveType::Travel;
    gcode.moves[3] = gcode.moves[1];
    gcode.moves[3].extruder_id = 1;
    Thermoelastic::AnalysisLimits limits;
    limits.max_subdivisions = 1;
    auto result = Thermoelastic::analyze_configured_gcode(gcode, config, limits);
    CHECK(result.state == Thermoelastic::AnalysisState::Unconverged);
    config.set_deserialize_strict({{"chamber_temperature", "50,55"}});
    result = Thermoelastic::analyze_configured_gcode(gcode, config, limits);
    CHECK(result.state == Thermoelastic::AnalysisState::InvalidGeometry);
}

TEST_CASE("Mechanical candidate ranking requires improvement beyond uncertainty in every scenario", "[ThermoelasticGCode]")
{
    Thermoelastic::ConfiguredAnalysis baseline, candidate;
    baseline.state = candidate.state = Thermoelastic::AnalysisState::Converged;
    baseline.scenarios = {{10., 2.}, {20., 4.}};
    candidate.scenarios = {{8., 1.5}, {15., 3.}};
    CHECK(Thermoelastic::mechanically_improves(baseline, candidate, .05));
    candidate.scenarios[1].tensile_force_n = 19.9;
    CHECK_FALSE(Thermoelastic::mechanically_improves(baseline, candidate, .05));
    candidate.scenarios[1].tensile_force_n = 15.;
    candidate.state = Thermoelastic::AnalysisState::Unconverged;
    CHECK_FALSE(Thermoelastic::mechanically_improves(baseline, candidate, .05));
    candidate.state = Thermoelastic::AnalysisState::Converged;
    candidate.scenarios.pop_back();
    CHECK_FALSE(Thermoelastic::mechanically_improves(baseline, candidate, .05));
    CHECK_FALSE(Thermoelastic::mechanically_improves(baseline, baseline, .05));
}

TEST_CASE("Mechanical preview status follows result assignment and is cleared on reset", "[ThermoelasticGCode]")
{
    using Status = GCodeProcessorResult::ThermalAnalysis::Status;
    GCodeProcessorResult source;
    source.thermal_stress_mpa = {1.f, 2.f};
    source.thermal_analysis = {Status::Converged, 3., 4., 2};
    GCodeProcessorResult copied;
    copied = source;
    CHECK(copied.thermal_analysis.status == Status::Converged);
    CHECK_THAT(copied.thermal_analysis.tensile_force_n, WithinAbs(3., 1e-12));
    CHECK_THAT(copied.thermal_analysis.elastic_energy_j, WithinAbs(4., 1e-12));
    CHECK(copied.thermal_analysis.scenario_count == 2);
    CHECK(copied.thermal_stress_mpa == source.thermal_stress_mpa);
    GCodeProcessorResult moved;
    moved = std::move(copied);
    CHECK(moved.thermal_analysis.status == Status::Converged);
    REQUIRE(moved.thermal_stress_mpa.size() == 2);
    moved.reset();
    CHECK(moved.thermal_analysis.status == Status::NotRun);
    CHECK(moved.thermal_analysis.scenario_count == 0);
    CHECK(moved.thermal_stress_mpa.empty());
}

TEST_CASE("Mechanical G-code analysis excludes travel and prime tower and publishes only converged samples", "[ThermoelasticGCode]")
{
    GCodeProcessorResult gcode;
    gcode.moves.resize(4);
    gcode.moves[0].position = Eigen::Vector3f(0.f, 0.f, .2f);
    auto &line = gcode.moves[1];
    line.position = Eigen::Vector3f(1.f, 0.f, .2f);
    line.type = EMoveType::Extrude;
    line.extrusion_role = erPerimeter;
    line.width = .4f;
    line.height = .2f;
    gcode.moves[2].position = Eigen::Vector3f(10.f, 10.f, .2f);
    gcode.moves[2].type = EMoveType::Travel;
    gcode.moves[3] = line;
    gcode.moves[3].position = Eigen::Vector3f(11.f, 10.f, .2f);
    gcode.moves[3].extrusion_role = erWipeTower;
    // Zero eigenstrain is an exact equilibrium control, not an assumed ASA
    // load case. Geometry still comes from the actual processor move data.
    const std::vector<Thermoelastic::Voxel> materials{{{0,0,0}, 2.9e9, .35}};
    const Eigen::Vector3d spacing(.0002,.0002,.0001);
    const auto analysis = Thermoelastic::analyze_gcode(gcode, materials, spacing, 0., 50000, 4, .05);
    REQUIRE(analysis.converged);
    REQUIRE(gcode.thermal_stress_mpa.size() == gcode.moves.size());
    CHECK_THAT(gcode.thermal_stress_mpa[1], WithinAbs(0., 1e-8));
    CHECK(gcode.thermal_stress_mpa[0] < 0.f);
    CHECK(gcode.thermal_stress_mpa[2] < 0.f);
    CHECK(gcode.thermal_stress_mpa[3] < 0.f);
    const auto coarse = Thermoelastic::analyze_gcode(gcode, materials, spacing, 0., 50000, 1, .05);
    CHECK_FALSE(coarse.converged);
    CHECK(gcode.thermal_stress_mpa.empty());
    auto contracting = materials;
    contracting[0].eigenstrain.head<3>().setConstant(-.005);
    const auto loaded = Thermoelastic::analyze_gcode(gcode, contracting, spacing, 0., 50000, 1, .05);
    CHECK(loaded.bed.equilibrium.elastic_energy_j > 0.);
    CHECK(loaded.bed.tensile_force_n > 0.);
    contracting[0].young_pa *= .5;
    const auto softer = Thermoelastic::analyze_gcode(gcode, contracting, spacing, 0., 50000, 1, .05);
    CHECK_THAT(softer.bed.equilibrium.elastic_energy_j,
        WithinAbs(loaded.bed.equilibrium.elastic_energy_j * .5, 1e-14));
    CHECK_THAT(softer.bed.tensile_force_n, WithinAbs(loaded.bed.tensile_force_n * .5, 1e-10));
    gcode.thermal_stress_mpa.assign(4, 123.f);
    CHECK_THROWS_AS(Thermoelastic::analyze_gcode(gcode, materials, spacing, 0., 50000, 4, .05,
        {}, [] { throw std::runtime_error("cancelled"); }), std::runtime_error);
    CHECK(gcode.thermal_stress_mpa.empty());
}

TEST_CASE("Sliced extrusion geometry can be assembled into a bonded mechanical mesh", "[ThermoelasticGCode]")
{
    auto config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({{"layer_height", "0.2"}, {"initial_layer_print_height", "0.2"},
        {"nozzle_diameter", "0.4"}, {"brim_type", "no_brim"}, {"skirt_loops", "0"}});
    Print print;
    Model model;
    Test::init_print({Test::cube(1.)}, print, model, config);
    GCodeProcessorResult gcode;
    REQUIRE_FALSE(Test::gcode(print, &gcode).empty());
    const std::vector<Thermoelastic::Voxel> materials{{{0,0,0}, 2.9e9, .35}};
    const auto analysis = Thermoelastic::analyze_gcode(gcode, materials, {.0001,.0001,.0001}, 0., 50000, 1, .05);
    CHECK_FALSE(analysis.mesh.elements.empty());
    CHECK(analysis.bed.contact_area_total_m2 > 0.);
    CHECK_FALSE(analysis.converged);
    CHECK(gcode.thermal_stress_mpa.empty());
}
