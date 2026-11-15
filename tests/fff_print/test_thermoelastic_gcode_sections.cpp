#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ThermoelasticGCodeSections.hpp"
#include "libslic3r/ThermoelasticToolpath.hpp"
#include "libslic3r/ThermoelasticSectionStudy.hpp"
#include "libslic3r/ThermoelasticSections.hpp"
#include "libslic3r/ThermoelasticMaterial.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticLayers.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Utils.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "test_utils.hpp"
#include <boost/nowide/fstream.hpp>
#include <stdexcept>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;

TEST_CASE("Rejected material contacts identify source moves without assigning ownership", "[ThermoelasticGCodeSections]")
{
    TimedToolpath path;
    path.state = ToolpathState::Complete;
    TimedDeposition first;
    first.geometry = {{0., 0., .001}, {.002, 0., .001}, .001, .001, 0};
    first.move_index = 7;
    first.gcode_id = 99;
    first.role = erExternalPerimeter;
    auto second = first;
    second.geometry.start_m.y() = second.geometry.end_m.y() = .00075;
    second.geometry.material = 2;
    second.move_index = 9;
    second.gcode_id = 105;
    second.role = erSolidInfill;
    auto unrelated = first;
    unrelated.geometry.start_m.z() = unrelated.geometry.end_m.z() = .003;
    path.deposition = {first, unrelated, second};
    bool rejected = false;
    try {
        section_deposition({first.geometry, second.geometry}, .0005, 0.);
    } catch (const LayerMaterialOverlap &overlap) {
        rejected = true;
        REQUIRE_FALSE(overlap.section_footprint.empty());
        const auto report = diagnose_section_overlap(path, overlap, BeadSection::RectangularEnvelope, 0.);
        CHECK(report.state == SectionOverlapDiagnostic::State::Complete);
        REQUIRE(report.contributors.size() == 2);
        CHECK(report.polygon_checks == 2);
        CHECK(report.scanned_segments == 3);
        CHECK(report.contributors[0].move_index == first.move_index);
        CHECK(report.contributors[0].gcode_id == first.gcode_id);
        CHECK(report.contributors[0].role == first.role);
        CHECK(report.contributors[1].filament == 2);
        CHECK(report.contributors[1].role == second.role);
        CHECK(report.contributors[1].overlap_area_m2 > 0.);
        auto crowded = path;
        crowded.deposition.insert(crowded.deposition.begin() + 1, 5, first);
        const auto representative = diagnose_section_overlap(crowded, overlap, BeadSection::RectangularEnvelope, 0.);
        CHECK(representative.state == SectionOverlapDiagnostic::State::BudgetExceeded);
        REQUIRE(representative.contributors.size() == 5);
        CHECK(representative.contributors.back().filament == 2);
        const auto bounded = diagnose_section_overlap(path, overlap, BeadSection::RectangularEnvelope, 0., 1);
        CHECK(bounded.state == SectionOverlapDiagnostic::State::BudgetExceeded);
        CHECK(bounded.contributors.size() == 1);
        auto partial = overlap;
        partial.section_footprint_complete = false;
        CHECK(diagnose_section_overlap(path, partial, BeadSection::RectangularEnvelope, 0.).state ==
              SectionOverlapDiagnostic::State::PartialMask);
        CHECK_THROWS_AS(diagnose_section_overlap(path, overlap, BeadSection::RectangularEnvelope, 0., 256, 32,
            [] { throw std::runtime_error("cancel"); }), std::runtime_error);
        auto absent = overlap;
        absent.section_footprint.clear();
        CHECK(diagnose_section_overlap(path, absent, BeadSection::RectangularEnvelope, 0.).state ==
              SectionOverlapDiagnostic::State::MaskUnavailable);
    }
    REQUIRE(rejected);
    CHECK(path.deposition[2].geometry.material == 2);
}

TEST_CASE("G-code section studies retain explicit materials and partial deposition", "[ThermoelasticGCodeSections]")
{
    const bool was_bbl = GCodeProcessor::s_IsBBLPrinter;
    const ScopeGuard restore([was_bbl] { GCodeProcessor::s_IsBBLPrinter = was_bbl; });
    GCodeProcessor::s_IsBBLPrinter = false;
    ScopedTemporaryFile input(".gcode");
    {
        boost::nowide::ofstream stream(input.string());
        stream << "G90\nM83\n;TYPE:Outer wall\n;WIDTH:1\n;HEIGHT:1\n"
                  "G1 X1 Y1 Z1 F600\nG1 X3 Y1 E0.4 F600\n";
        REQUIRE(stream.good());
    }
    GCodeProcessor processor;
    processor.apply_config(FullPrintConfig{});
    processor.process_file(input.string());
    const auto path = timed_deposition(processor.get_result());
    REQUIRE(path.state == ToolpathState::Complete);
    REQUIRE(path.deposition.size() == 1);
    GCodeSectionScenario scenario(BeadSection::RectangularEnvelope);
    ElasticityMatrix d = ElasticityMatrix::Zero();
    d.diagonal().head<3>().setConstant(1000.);
    d.diagonal().tail<3>().setConstant(500.);
    scenario.materials = {{d}};
    scenario.filament_to_material = {0};
    scenario.initial_basis_spacing_m = {.001, .0005, .001};
    scenario.study.max_section_height_m = .001;
    scenario.study.relative_tolerance = .05;
    const auto full = study_gcode_sections(processor.get_result(), scenario);
    REQUIRE(full.study.field.has_value());
    CHECK(full.study.state == DepositionSectionStudy::State::NominalGeometryAndBasisConverged);
    CHECK(full.commanded_path_volume_m3 > 0.);
    const auto &bead = path.deposition.front();
    scenario.observation_time_s = (bead.start_s + bead.end_s) / 2.;
    const auto partial = study_gcode_sections(processor.get_result(), scenario);
    REQUIRE(partial.study.field.has_value());
    CHECK(partial.study.field->occupied_volume_m3 < full.study.field->occupied_volume_m3);
    CHECK_THAT(partial.commanded_path_volume_m3,
               Catch::Matchers::WithinRel(full.commanded_path_volume_m3 / 2., 1e-8));
    scenario.observation_time_s.reset();
    GCodeProcessorResult remapped;
    remapped.moves = processor.get_result().moves;
    remapped.print_statistics = processor.get_result().print_statistics;
    remapped.filament_diameters = processor.get_result().filament_diameters;
    for (auto &move : remapped.moves)
        if (move.type == EMoveType::Extrude) move.extruder_id = 1;
    CHECK_THROWS_AS(study_gcode_sections(remapped, scenario), std::invalid_argument);
    scenario.filament_to_material = {0, 0};
    const auto mapped = study_gcode_sections(remapped, scenario);
    REQUIRE(mapped.study.field.has_value());
    CHECK_THAT(mapped.study.field->occupied_volume_m3,
               Catch::Matchers::WithinRel(full.study.field->occupied_volume_m3, 1e-8));
    CHECK(remapped.moves.back().extruder_id == 1);
    for (auto &move : remapped.moves)
        if (move.type == EMoveType::Extrude) move.extrusion_role = erMixed;
    CHECK_THROWS_AS(study_gcode_sections(remapped, scenario), std::invalid_argument);
    for (auto &move : remapped.moves)
        if (move.type == EMoveType::Extrude) move.extrusion_role = erWipeTower;
    // A tower-only input cannot masquerade as a qualified part field.
    CHECK_THROWS_AS(study_gcode_sections(remapped, scenario), std::invalid_argument);
    scenario.filament_to_material = {1};
    CHECK_THROWS_AS(study_gcode_sections(processor.get_result(), scenario), std::invalid_argument);
    scenario.filament_to_material = {0};
    scenario.study.max_total_work = 1;
    const auto bounded = study_gcode_sections(processor.get_result(), scenario);
    CHECK(bounded.study.state == DepositionSectionStudy::State::ResourceLimited);
    CHECK_FALSE(bounded.study.field.has_value());
    CHECK_THROWS_AS(study_gcode_sections(processor.get_result(), scenario, 1), ResourceLimit);
    CHECK_THROWS_AS(study_gcode_sections(processor.get_result(), scenario, 2000000, 250000,
        {}, {}, [] { throw std::runtime_error("cancel"); }), std::runtime_error);
}
