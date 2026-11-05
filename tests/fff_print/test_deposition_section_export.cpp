#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <cstddef>
#include <cstdint>
#include <array>
#include "libslic3r/GCode.hpp"
#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/GCodeWriter.hpp"
#include "libslic3r/GCode/DepositionSection.hpp"
#include "libslic3r/GCode/DepositionProvenance.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/ThermoelasticToolpath.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/Config.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/ExtrusionEntityCollection.hpp"
#include "libslic3r/Polyline.hpp"
#include "libslic3r/libslic3r.h"
#include "test_helpers.hpp"
#include "test_utils.hpp"
#include <boost/nowide/cstdio.hpp>
#include <boost/nowide/fstream.hpp>
#include <iterator>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <utility>

using namespace Slic3r;
namespace {
class SectionStreamAccess : public GCode
{
public:
    using GCode::GCodeOutputStream;
    using GCode::emit_deposition_section;
    using GCode::m_deposition_sections;
};
} // namespace

TEST_CASE("Precise arc evidence uses the exact parsed decimal endpoints", "[DepositionSectionExport]")
{
    const auto purpose = GENERATE(DepositionCommandSection::Purpose::Model, DepositionCommandSection::Purpose::AuxiliaryBrim);
    SectionStreamAccess generator;
    PrintConfig config;
    config.gcode_flavor.value = gcfKlipper;
    config.use_relative_e_distances.value = true;
    generator.apply_print_config(config);
    generator.writer().set_extruders({0});
    generator.writer().toolchange(0, 0);
    generator.writer().set_position(Vec3d(2.037, 6.694, .4));
    generator.m_deposition_sections = std::make_shared<DepositionSections>(42);
    std::string command;
    generator.emit_deposition_section(command, Vec3d(1.269, 6.880, .4), .1, "", false,
        .4, .2, .2, true, Vec2d(-2.152, -7.208), true, {}, DepositionCommandSection::Transform::Ordinary, purpose);
    REQUIRE(generator.m_deposition_sections->size() == 1);
    const auto& parent = *generator.m_deposition_sections->section(1);
    CHECK(parent.purpose == purpose);
    GCodeReader reader;
    bool found = false;
    reader.parse_buffer(command, [&](GCodeReader&, const GCodeReader::GCodeLine& line) {
        if (!line.cmd_is("G3")) return;
        found = true;
        // Exact equality is the contract: reciprocal multiplication used to
        // shift the decimal 1.269 by one ULP and miss the commanded endpoint.
        CHECK((parent.end[0] == line.precise_value(X)));
        CHECK((parent.end[1] == line.precise_value(Y)));
        DepositionSectionChain chain(parent);
        REQUIRE(chain.append(parent.start, {line.precise_value(X), line.precise_value(Y), parent.end[2]}, 0));
        CHECK(chain.complete());
    });
    CHECK(found);
}

TEST_CASE("Only authenticated vertical brick deposition retains material without changing travel timing", "[DepositionSectionExport]")
{
    // Missing source/section, ordinary transform, damaged span, wrong filament,
    // nonpositive E, non-wall role and untransformed source cannot authenticate.
    const int invalid = GENERATE(0, 1, 2, 3, 4, 5, 6, 7, 8, 9);
    CAPTURE(invalid);
    auto sources = std::make_shared<DepositionProvenance>();
    DepositionSource source;
    source.height = invalid == 9 ? DepositionSource::Height::NativePlanar : DepositionSource::Height::Transformed;
    source.logical_filament = invalid == 5 ? 1 : 0;
    const auto source_id = sources->add(source);
    auto sections = std::make_shared<DepositionSections>(sources->generation());
    DepositionCommandSection command{{1., 2., .4}, {1., 2., .5}, .45, .3, .3, 0};
    command.transform = invalid == 3 ? DepositionCommandSection::Transform::Ordinary :
        DepositionCommandSection::Transform::BrickCourse;
    const auto id = sections->add(command);
    GCodeProcessor processor, baseline;
    FullPrintConfig config;
    config.gcode_flavor.value = gcfKlipper;
    processor.apply_config(config);
    baseline.apply_config(config);
    processor.initialize_result_moves();
    baseline.initialize_result_moves();
    ScopedTemporaryFile output(".gcode");
    SectionStreamAccess::GCodeOutputStream stream(boost::nowide::fopen(output.string().c_str(), "wb"), processor);
    REQUIRE(stream.is_open());
    stream.enable_provenance(sources->generation());
    const std::string header = "G90\nM83\n;" + GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Role) +
        (invalid == 8 ? "Sparse infill" : "Brick wall") + "\nG1 X1 Y2 Z0.4 F600\n";
    const std::string vertical = std::string("G1 X1 Y2 Z") + (invalid == 4 ? "0.45" : "0.5") +
        " E" + (invalid == 6 ? "0" : invalid == 7 ? "-0.01" : "0.01") + "\n";
    // Neither unmarked vertical extrusion nor stationary unretraction acquires
    // the authenticated flag, even while the same semantic role is active.
    const std::string tail = "G1 X1 Y2 Z0.6 E0.01\nG1 E0.01\n";
    stream.write(header);
    if (invalid != 1) stream.write(sources->begin(source_id));
    if (invalid != 2) stream.write(sections->marker(id));
    stream.write(vertical);
    if (invalid != 2) stream.write(sections->marker(0));
    if (invalid != 1) stream.write(sources->clear());
    stream.write(tail);
    sources->seal();
    sections->seal();
    REQUIRE(stream.finish_provenance(sources, invalid == 1 ? std::vector<std::uint32_t>{} :
        std::vector<std::uint32_t>{source_id, 0}));
    const bool complete_transcript = stream.finish_sections(*sections);
    CHECK(complete_transcript == (invalid != 2));
    processor.finish_deposition_sections(complete_transcript ? sections : nullptr);
    for (const auto& move : processor.result().moves) CHECK_FALSE(move.deposition_vertical_extrusion);
    processor.finalize(false);
    processor.finish_deposition_provenance(sources);
    stream.close();
    baseline.process_buffer(header + vertical + tail);
    baseline.finalize(false);
    const auto& actual = processor.result();
    const auto& expected = baseline.result();
    REQUIRE(actual.moves.size() == expected.moves.size());
    std::size_t flagged = 0;
    for (std::size_t i = 0; i < actual.moves.size(); ++i) {
        const auto& move = actual.moves[i];
        CHECK(move.type == expected.moves[i].type);
        CHECK(move.time == expected.moves[i].time);
        CHECK_THAT(move.delta_extruder, Catch::Matchers::WithinAbs(expected.moves[i].delta_extruder, 1e-12));
        if (!move.deposition_vertical_extrusion) continue;
        ++flagged;
        CHECK(move.type == EMoveType::Travel);
        REQUIRE(move.deposition_design_width);
        REQUIRE(move.deposition_design_heights);
        CHECK_THAT(*move.deposition_design_width, Catch::Matchers::WithinAbs(.45, 1e-6));
        CHECK_THAT((*move.deposition_design_heights)[0], Catch::Matchers::WithinAbs(.3, 1e-6));
        CHECK_THAT((*move.deposition_design_heights)[1], Catch::Matchers::WithinAbs(.3, 1e-6));
        CHECK_THAT(move.delta_extruder, Catch::Matchers::WithinAbs(.01, 1e-6));
    }
    CHECK(flagged == (invalid == 0 ? 1 : 0));
    boost::nowide::ifstream input(output.string());
    const std::string bytes((std::istreambuf_iterator<char>(input)), {});
    CHECK(bytes == header + vertical + tail);
    processor.finish_deposition_provenance(nullptr);
    for (const auto& move : processor.result().moves) CHECK_FALSE(move.deposition_vertical_extrusion);
}

TEST_CASE("Perimeter and ordinary infill evidence preserves every public straight and circular command", "[DepositionSectionExport]")
{
    const bool transformed = GENERATE(false, true);
    const bool arcs = GENERATE(false, true);
    const auto flavor = GENERATE("marlin", "marlin2", "klipper");
    CAPTURE(transformed, arcs, flavor);
    const auto config = Test::multifilament_config(1, {{"layer_height", "0.2"},
        {"gcode_flavor", flavor},
        {"initial_layer_print_height", "0.2"}, {"brim_type", "no_brim"}, {"seam_slope_type", "none"}, {"enable_arc_fitting", arcs ? "1" : "0"},
        {"enable_prime_tower", "0"}, {"time_lapse_gcode", ""}, {"spiral_mode", "0"},
        {"max_volumetric_extrusion_rate_slope", "0.5"}, {"thick_bridges", transformed ? "1" : "0"},
        {"thick_internal_bridges", "0"}, {"top_shell_layers", "3"}, {"bottom_shell_layers", "3"},
        {"top_shell_thickness", "0"}, {"bottom_shell_thickness", "0"}, {"sparse_infill_density", "15%"}});
    Print print;
    Model model;
    Test::init_print({arcs ? make_cylinder(8., 3.) : make_cube(8., 8., 3.)}, print, model, config);
    print.set_status_silent();
    print.process();
    ScopedTemporaryFile plain(".gcode"), tracked(".gcode");
    GCode normal, authenticated;
    GCodeProcessorResult result;
    normal.do_export(&print, plain.string().c_str());
    authenticated.enable_deposition_provenance();
    authenticated.do_export(&print, tracked.string().c_str(), &result);
    INFO(result.deposition_sections_failure);
    REQUIRE(result.deposition_sections);
    REQUIRE(result.deposition_provenance);
    CHECK(result.deposition_sections->sealed());
    CHECK(result.deposition_sections->generation() == result.deposition_provenance->generation());
    CHECK(result.deposition_sections->size() <= DepositionSections::maximum_commands);
    std::size_t constant_perimeters = 0;
    std::array<std::size_t, 4> infill_counts{};
    std::size_t bottom_moves = 0;
    for (const auto& move : result.moves) {
        if (!move.internal_only && move.type == EMoveType::Extrude && move.extrusion_role == erBottomSurface) {
            ++bottom_moves;
            REQUIRE(move.deposition_source_id != 0);
            const auto* source = result.deposition_provenance->source(move.deposition_source_id);
            REQUIRE(source);
            CHECK((source->height == DepositionSource::Height::Transformed) == transformed);
            CHECK(bool(move.deposition_command_section_id) == transformed);
            CHECK(bool(move.deposition_design_heights) == transformed);
        }
        if (!move.deposition_design_heights) continue;
        const bool perimeter = move.extrusion_role == erPerimeter || move.extrusion_role == erExternalPerimeter;
        const bool infill = move.extrusion_role == erInternalInfill || move.extrusion_role == erSolidInfill ||
                            move.extrusion_role == erTopSolidInfill || move.extrusion_role == erBottomSurface ||
                            move.extrusion_role == erBridgeInfill || move.extrusion_role == erInternalBridgeInfill ||
                            move.extrusion_role == erGapFill;
        CHECK((perimeter || infill));
        CHECK_THAT((*move.deposition_design_heights)[0], Catch::Matchers::WithinAbs((*move.deposition_design_heights)[1], 1e-8));
        REQUIRE(move.deposition_design_width);
        const auto* parent = result.deposition_sections->section(move.deposition_command_section_id);
        REQUIRE(parent);
        CHECK(parent->filament == move.extruder_id);
        CHECK_THAT(*move.deposition_design_width, Catch::Matchers::WithinAbs(parent->width, 1e-6));
        CHECK_THAT((*move.deposition_design_heights)[0], Catch::Matchers::WithinAbs(parent->start_height, 1e-6));
        CHECK_THAT((*move.deposition_design_heights)[1], Catch::Matchers::WithinAbs(parent->end_height, 1e-6));
        constant_perimeters += perimeter;
        if (move.extrusion_role == erInternalInfill) ++infill_counts[0];
        if (move.extrusion_role == erSolidInfill) ++infill_counts[1];
        if (move.extrusion_role == erTopSolidInfill) ++infill_counts[2];
        if (move.extrusion_role == erBottomSurface) ++infill_counts[3];
    }
    REQUIRE(bottom_moves > 0);
    CHECK((constant_perimeters > 0) == transformed);
    for (const auto count : infill_counts) CHECK((count > 0) == transformed);
    std::size_t arc_parents = 0;
    for (std::uint32_t id = 1; id <= result.deposition_sections->size(); ++id)
        arc_parents += result.deposition_sections->section(id)->arc_direction != 0;
    CHECK((arc_parents > 0) == (transformed && arcs));
    // Export timestamps and estimate annotations are not printer commands.
    // Every non-comment byte, including whitespace and inline comments, stays
    // unchanged by the opt-in evidence collection.
    const auto commands = [](const std::string& filename) {
        boost::nowide::ifstream input(filename);
        std::string line, output;
        while (std::getline(input, line)) {
            CHECK(line.find("_ORCA_INTERNAL_SECTION") == std::string::npos);
            if (!line.empty() && line.front() != ';') output += line + '\n';
        }
        return output;
    };
    const bool identical_commands = commands(plain.string()) == commands(tracked.string());
    CHECK(identical_commands);
}

TEST_CASE("Actual ContourZ paths retain anchored envelopes without changing public commands", "[DepositionSectionExport]")
{
    const bool compensated = GENERATE(false, true);
    const bool inner_seam = GENERATE(false, true);
    CAPTURE(compensated, inner_seam);
    const auto config = Test::multifilament_config(1, {{"layer_height", "0.2"},
        {"initial_layer_print_height", "0.2"}, {"perimeter_layering", "standard"},
        {"top_surface_z_mode", "z_contouring"}, {"zaa_enabled", "1"}, {"zaa_min_z", "0.05"},
        {"nonplanar_top_surface", "0"}, {"wall_loops", "3"}, {"wall_generator", "classic"},
        {"seam_start_on_inner_wall", inner_seam ? "1" : "0"}, {"wall_sequence", "inner wall/outer wall"},
        {"brim_type", "no_brim"}, {"skirt_loops", "0"}, {"seam_slope_type", "none"},
        {"ironing_type", "no ironing"}, {"enable_arc_fitting", "0"},
        {"enable_prime_tower", "0"}, {"time_lapse_gcode", ""}, {"spiral_mode", "0"},
        {"max_volumetric_extrusion_rate_slope", "0"}});
    auto mesh = make_cube(8., 8., 2.).its;
    for (auto& vertex : mesh.vertices)
        if (vertex.z() > 0.) vertex.z() += .25f * vertex.x();
    Print print;
    Model model;
    Test::init_print({TriangleMesh(mesh)}, print, model, config);
    print.set_status_silent();
    print.process();
    // Isolate the orthogonal preview modifier from the geometric witness.
    // PerimeterGenerator's shrinkage pass sets this same bool after its XY
    // edits, before ContourZ; it changes neither scheduled height nor Z.
    const auto mark_compensated = [&](auto&& self, ExtrusionEntity& entity) -> void {
        if (auto* path = dynamic_cast<ExtrusionPath*>(&entity)) {
            if (path->contour_z_generated) path->shrinkage_compensation = compensated;
        } else if (auto* loop = dynamic_cast<ExtrusionLoop*>(&entity)) {
            for (auto& path : loop->paths) self(self, path);
        } else if (auto* multi = dynamic_cast<ExtrusionMultiPath*>(&entity)) {
            for (auto& path : multi->paths) self(self, path);
        } else if (auto* collection = dynamic_cast<ExtrusionEntityCollection*>(&entity)) {
            for (auto* child : collection->entities) self(self, *child);
        }
    };
    for (auto* object : print.objects())
        for (auto* layer : object->layers())
            for (auto* region : layer->regions()) mark_compensated(mark_compensated, region->perimeters);
    ScopedTemporaryFile plain(".gcode"), tracked(".gcode"), flag_only(".gcode");
    GCode normal, authenticated;
    GCodeProcessorResult result;
    normal.do_export(&print, plain.string().c_str());
    authenticated.enable_deposition_provenance();
    authenticated.do_export(&print, tracked.string().c_str(), &result);
    INFO(result.deposition_sections_failure);
    REQUIRE(result.deposition_sections);
    REQUIRE(result.deposition_provenance);
    bool sloped = false, flat = false, compensated_envelope = false;
    for (const auto& move : result.moves) {
        if (move.internal_only || !move.deposition_command_section_id) continue;
        const auto* parent = result.deposition_sections->section(move.deposition_command_section_id);
        REQUIRE(parent);
        // Standard-mode seam preparation must not certify its connectors as
        // brick courses when genuine ContourZ paths coexist in this layer.
        CHECK(parent->transform == DepositionCommandSection::Transform::Ordinary);
        if (!parent->nominal_lower_z) continue;
        compensated_envelope |= move.extrusion_role == erShrinkageCompensation;
        const auto* source = result.deposition_provenance->source(move.deposition_source_id);
        REQUIRE(source);
        const auto* layer = print.objects().front()->get_layer(source->layer);
        REQUIRE(layer);
        const double base = std::round(layer->bottom_z() * 1e6) / 1e6;
        CHECK(parent->purpose == DepositionCommandSection::Purpose::Model);
        CHECK(parent->transform == DepositionCommandSection::Transform::Ordinary);
        CHECK_THAT((*parent->nominal_lower_z)[0], Catch::Matchers::WithinAbs(base, 1e-12));
        CHECK_THAT((*parent->nominal_lower_z)[1], Catch::Matchers::WithinAbs(base, 1e-12));
        CHECK_THAT(parent->start_height, Catch::Matchers::WithinAbs(parent->start[2] - base, 1e-12));
        CHECK_THAT(parent->end_height, Catch::Matchers::WithinAbs(parent->end[2] - base, 1e-12));
        REQUIRE(move.deposition_design_heights);
        REQUIRE(move.deposition_design_width);
        sloped |= std::abs(parent->start[2] - parent->end[2]) > 1e-6;
        flat |= std::abs(parent->start[2] - parent->end[2]) < 1e-12;
    }
    REQUIRE(sloped);
    REQUIRE(flat);
    CHECK(compensated_envelope == compensated);

    // Preserve exactly the generated geometry and shared z_contoured flag.
    // Without ContourZ provenance those same paths must not acquire an anchor.
    std::size_t tagged = 0;
    const auto clear_tag = [&](auto&& self, ExtrusionEntity& entity) -> void {
        if (auto* path = dynamic_cast<ExtrusionPath*>(&entity)) {
            if (path->contour_z_generated) {
                REQUIRE(path->z_contoured);
                // The emission pipeline splits and moves paths: all value
                // operations must preserve the geometric provenance tag.
                ExtrusionPath copied(*path), moved(std::move(copied));
                ExtrusionPath replaced(path->polyline, moved);
                ExtrusionPath replaced_move(Polyline3(path->polyline), moved);
                copied = replaced;
                replaced = std::move(replaced_move);
                CHECK(moved.contour_z_generated);
                CHECK(copied.contour_z_generated);
                CHECK(replaced.contour_z_generated);
                path->contour_z_generated = false;
                ++tagged;
            }
        } else if (auto* loop = dynamic_cast<ExtrusionLoop*>(&entity)) {
            for (auto& path : loop->paths) self(self, path);
        } else if (auto* multi = dynamic_cast<ExtrusionMultiPath*>(&entity)) {
            for (auto& path : multi->paths) self(self, path);
        } else if (auto* collection = dynamic_cast<ExtrusionEntityCollection*>(&entity)) {
            for (auto* child : collection->entities) self(self, *child);
        }
    };
    for (auto* object : print.objects())
        for (auto* layer : object->layers())
            for (auto* region : layer->regions()) {
                clear_tag(clear_tag, region->perimeters);
                clear_tag(clear_tag, region->fills);
            }
    REQUIRE(tagged > 0);
    GCode untagged;
    GCodeProcessorResult untagged_result;
    untagged.enable_deposition_provenance();
    untagged.do_export(&print, flag_only.string().c_str(), &untagged_result);
    INFO(untagged_result.deposition_sections_failure);
    REQUIRE(untagged_result.deposition_sections);
    for (const auto& move : untagged_result.moves) {
        if (move.internal_only || !move.deposition_command_section_id) continue;
        const auto* parent = untagged_result.deposition_sections->section(move.deposition_command_section_id);
        REQUIRE(parent);
        CHECK_FALSE(parent->nominal_lower_z);
    }
    const auto commands = [](const std::string& filename) {
        boost::nowide::ifstream input(filename);
        std::string line, output;
        while (std::getline(input, line)) {
            CHECK(line.find("_ORCA_INTERNAL_SECTION") == std::string::npos);
            if (!line.empty() && line.front() != ';') output += line + '\n';
        }
        return output;
    };
    CHECK(commands(plain.string()) == commands(tracked.string()));
    CHECK(commands(plain.string()) == commands(flag_only.string()));
}

TEST_CASE("Planar overhang command evidence retains scheduled bridge dimensions and public commands", "[DepositionSectionExport]")
{
    const bool thick = GENERATE(false, true);
    const bool overhang_speed = GENERATE(false, true);
    CAPTURE(thick, overhang_speed);
    const auto config = Test::multifilament_config(1, {{"layer_height", "0.2"},
        {"initial_layer_print_height", "0.2"}, {"nozzle_diameter", "0.4"}, {"bridge_line_width", "0.4"},
        {"bridge_flow", "1"}, {"thick_bridges", thick ? "1" : "0"},
        // Keep the region explicitly transformed even in the thin-bridge case.
        {"thick_internal_bridges", "1"}, {"wall_generator", "classic"}, {"wall_loops", "2"},
        {"detect_overhang_wall", "1"}, {"arc_overhang_enabled", "0"},
        {"enable_overhang_speed", overhang_speed ? "1" : "0"}, {"slowdown_for_curled_perimeters", "0"},
        {"brim_type", "no_brim"}, {"skirt_loops", "0"}, {"seam_slope_type", "none"},
        {"enable_arc_fitting", "0"}, {"enable_prime_tower", "0"}, {"time_lapse_gcode", ""},
        {"spiral_mode", "0"}, {"max_volumetric_extrusion_rate_slope", "0.5"}});
    auto mesh = make_cube(4., 4., 2.);
    auto roof = make_cube(8., 8., 1.);
    roof.translate(-2., -2., 2.);
    mesh.merge(roof);
    Print print;
    Model model;
    Test::init_print({mesh}, print, model, config);
    print.set_status_silent();
    print.process();
    ScopedTemporaryFile plain(".gcode"), tracked(".gcode");
    GCode normal, authenticated;
    GCodeProcessorResult result;
    normal.do_export(&print, plain.string().c_str());
    authenticated.enable_deposition_provenance();
    authenticated.do_export(&print, tracked.string().c_str(), &result);
    INFO(result.deposition_sections_failure);
    REQUIRE(result.deposition_sections);
    REQUIRE(result.deposition_provenance);
    std::size_t overhang_moves = 0, bridge_moves = 0;
    for (const auto& move : result.moves) {
        if (move.internal_only || move.type != EMoveType::Extrude ||
            (move.extrusion_role != erOverhangPerimeter && move.extrusion_role != erBridgeInfill)) continue;
        overhang_moves += move.extrusion_role == erOverhangPerimeter;
        bridge_moves += move.extrusion_role == erBridgeInfill;
        const auto* source = result.deposition_provenance->source(move.deposition_source_id);
        REQUIRE(source);
        CHECK(source->height == DepositionSource::Height::Transformed);
        const auto* parent = result.deposition_sections->section(move.deposition_command_section_id);
        REQUIRE(parent);
        REQUIRE(move.deposition_design_heights);
        REQUIRE(move.deposition_design_width);
        CHECK(parent->purpose == DepositionCommandSection::Purpose::Model);
        CHECK(parent->transform == DepositionCommandSection::Transform::Ordinary);
        CHECK_THAT(parent->start[2], Catch::Matchers::WithinAbs(parent->end[2], 1e-12));
        CHECK_THAT(parent->width, Catch::Matchers::WithinAbs(.4, 1e-6));
        // Thick bridge height is its configured thread diameter, not the
        // owning 0.2mm layer height. Thin bridges retain the layer height.
        const double height = thick ? .4 : .2;
        CHECK_THAT(parent->start_height, Catch::Matchers::WithinAbs(height, 1e-6));
        CHECK_THAT(parent->end_height, Catch::Matchers::WithinAbs(height, 1e-6));
        CHECK_THAT((*move.deposition_design_heights)[0], Catch::Matchers::WithinAbs(height, 1e-6));
        CHECK_THAT((*move.deposition_design_heights)[1], Catch::Matchers::WithinAbs(height, 1e-6));
    }
    REQUIRE(overhang_moves > 0);
    REQUIRE(bridge_moves > 0);
    const auto commands = [](const std::string& filename) {
        boost::nowide::ifstream input(filename);
        std::string line, output;
        while (std::getline(input, line)) {
            CHECK(line.find("_ORCA_INTERNAL_SECTION") == std::string::npos);
            if (!line.empty() && line.front() != ';') output += line + '\n';
        }
        return output;
    };
    CHECK(commands(plain.string()) == commands(tracked.string()));
}

TEST_CASE("Internal bridges and gap fill retain generated dimensions without changing public commands", "[DepositionSectionExport]")
{
    const auto role = GENERATE(erInternalBridgeInfill, erGapFill);
    CAPTURE(role);
    const auto config = Test::multifilament_config(1, {{"layer_height", "0.2"},
        {"initial_layer_print_height", "0.2"}, {"thick_bridges", "1"}, {"thick_internal_bridges", "1"},
        {"bridge_line_width", "0.4"}, {"bridge_flow", "1"}, {"internal_bridge_flow", "0.8"},
        {"sparse_infill_pattern", "rectilinear"}, {"sparse_infill_density", "15%"},
        {"top_shell_layers", "3"}, {"bottom_shell_layers", "2"},
        {"top_shell_thickness", "0"}, {"bottom_shell_thickness", "0"},
        {"top_surface_pattern", "monotonicline"}, {"gap_fill_target", "topbottom"}, {"filter_out_gap_fill", "0"},
        {"brim_type", "no_brim"}, {"skirt_loops", "0"}, {"seam_slope_type", "none"},
        {"enable_arc_fitting", "0"}, {"enable_prime_tower", "0"}, {"time_lapse_gcode", ""},
        {"spiral_mode", "0"}, {"max_volumetric_extrusion_rate_slope", "0"}});
    auto mesh = make_cube(30., 24., 3.);
    auto tower = make_cube(14., 10., 1.);
    tower.translate(8., 7., 3.);
    mesh.merge(tower);
    if (role == erGapFill) mesh = Test::mesh(Test::TestMesh::A);
    Print print;
    Model model;
    Test::init_print({mesh}, print, model, config);
    print.set_status_silent();
    print.process();
    std::vector<std::array<double, 2>> scheduled;
    const auto collect = [&](auto&& self, const ExtrusionEntity& entity) -> void {
        if (const auto* path = dynamic_cast<const ExtrusionPath*>(&entity)) {
            if (path->role() == role) scheduled.push_back({path->width, path->height});
        } else if (const auto* loop = dynamic_cast<const ExtrusionLoop*>(&entity)) {
            for (const auto& path : loop->paths) self(self, path);
        } else if (const auto* multi = dynamic_cast<const ExtrusionMultiPath*>(&entity)) {
            for (const auto& path : multi->paths) self(self, path);
        } else if (const auto* collection = dynamic_cast<const ExtrusionEntityCollection*>(&entity)) {
            for (const auto* child : collection->entities) self(self, *child);
        }
    };
    for (const auto* object : print.objects())
        for (const auto* layer : object->layers())
            for (const auto* region : layer->regions()) collect(collect, region->fills);
    REQUIRE_FALSE(scheduled.empty());
    ScopedTemporaryFile plain(".gcode"), tracked(".gcode");
    GCode normal, authenticated;
    GCodeProcessorResult result;
    normal.do_export(&print, plain.string().c_str());
    authenticated.enable_deposition_provenance();
    authenticated.do_export(&print, tracked.string().c_str(), &result);
    INFO(result.deposition_sections_failure);
    REQUIRE(result.deposition_sections);
    REQUIRE(result.deposition_provenance);
    std::size_t captured = 0;
    for (const auto& move : result.moves) {
        if (move.internal_only || move.type != EMoveType::Extrude || move.extrusion_role != role) continue;
        ++captured;
        const auto* source = result.deposition_provenance->source(move.deposition_source_id);
        REQUIRE(source);
        CHECK(source->height == DepositionSource::Height::Transformed);
        const auto* parent = result.deposition_sections->section(move.deposition_command_section_id);
        REQUIRE(parent);
        CHECK(parent->purpose == DepositionCommandSection::Purpose::Model);
        CHECK(parent->transform == DepositionCommandSection::Transform::Ordinary);
        CHECK_FALSE(parent->nominal_lower_z);
        CHECK_THAT(parent->start[2], Catch::Matchers::WithinAbs(parent->end[2], 1e-12));
        CHECK_THAT(parent->start_height, Catch::Matchers::WithinAbs(parent->end_height, 1e-12));
        // Compare with the actual generator schedule, not preview dimensions
        // or E (internal_bridge_flow deliberately scales E independently).
        bool matches_schedule = false;
        for (const auto& dimensions : scheduled)
            matches_schedule |= std::abs(parent->width - dimensions[0]) < 1e-6 &&
                                std::abs(parent->start_height - dimensions[1]) < 1e-6;
        CHECK(matches_schedule);
        REQUIRE(move.deposition_design_heights);
        REQUIRE(move.deposition_design_width);
        CHECK_THAT(*move.deposition_design_width, Catch::Matchers::WithinAbs(parent->width, 1e-6));
        CHECK_THAT((*move.deposition_design_heights)[0], Catch::Matchers::WithinAbs(parent->start_height, 1e-6));
        CHECK_THAT((*move.deposition_design_heights)[1], Catch::Matchers::WithinAbs(parent->end_height, 1e-6));
    }
    REQUIRE(captured > 0);
    const auto commands = [](const std::string& filename) {
        boost::nowide::ifstream input(filename);
        std::string line, output;
        while (std::getline(input, line)) {
            CHECK(line.find("_ORCA_INTERNAL_SECTION") == std::string::npos);
            if (!line.empty() && line.front() != ';') output += line + '\n';
        }
        return output;
    };
    CHECK(commands(plain.string()) == commands(tracked.string()));
}

TEST_CASE("Exported brims retain typed auxiliary command evidence without fabricated model sources", "[DepositionSectionExport]")
{
    const bool arcs = GENERATE(false, true);
    const auto config = Test::multifilament_config(1, {{"layer_height", "0.2"},
        {"initial_layer_print_height", "0.2"}, {"brim_type", "outer_only"}, {"brim_width", "2"},
        {"brim_object_gap", "0.1"}, {"skirt_loops", "0"}, {"seam_slope_type", "none"},
        {"enable_arc_fitting", arcs ? "1" : "0"}, {"enable_prime_tower", "0"},
        {"time_lapse_gcode", ""}, {"spiral_mode", "0"}, {"max_volumetric_extrusion_rate_slope", "0.5"}});
    Print print;
    Model model;
    Test::init_print({make_cylinder(4., 1.)}, print, model, config);
    print.set_status_silent();
    print.process();
    ScopedTemporaryFile plain(".gcode"), tracked(".gcode");
    GCode normal, authenticated;
    GCodeProcessorResult result;
    normal.do_export(&print, plain.string().c_str());
    authenticated.enable_deposition_provenance();
    authenticated.do_export(&print, tracked.string().c_str(), &result);
    INFO(result.deposition_sections_failure);
    REQUIRE(result.deposition_sections);
    REQUIRE(result.deposition_provenance);
    CHECK(result.deposition_sections->generation() == result.deposition_provenance->generation());
    CHECK(result.deposition_sections->sealed());
    std::size_t brim_moves = 0;
    for (const auto& move : result.moves) {
        if (move.internal_only || move.type != EMoveType::Extrude) continue;
        const auto* parent = result.deposition_sections->section(move.deposition_command_section_id);
        if (move.extrusion_role != erBrim) {
            if (parent) CHECK(parent->purpose == DepositionCommandSection::Purpose::Model);
            continue;
        }
        ++brim_moves;
        REQUIRE(parent);
        CHECK(move.deposition_source_id == 0);
        CHECK(parent->purpose == DepositionCommandSection::Purpose::AuxiliaryBrim);
        CHECK(parent->transform == DepositionCommandSection::Transform::Ordinary);
        CHECK(parent->filament == move.extruder_id);
        CHECK_THAT(parent->start[2], Catch::Matchers::WithinAbs(parent->end[2], 1e-12));
        CHECK_THAT(parent->start_height, Catch::Matchers::WithinAbs(.2, 1e-6));
        CHECK_THAT(parent->end_height, Catch::Matchers::WithinAbs(.2, 1e-6));
        REQUIRE(move.deposition_design_heights);
        REQUIRE(move.deposition_design_width);
        CHECK_THAT(*move.deposition_design_width, Catch::Matchers::WithinAbs(parent->width, 1e-6));
        CHECK_THAT((*move.deposition_design_heights)[0], Catch::Matchers::WithinAbs(parent->start_height, 1e-6));
        CHECK_THAT((*move.deposition_design_heights)[1], Catch::Matchers::WithinAbs(parent->end_height, 1e-6));
    }
    CHECK(brim_moves > 0);
    const auto commands = [](const std::string& filename) {
        boost::nowide::ifstream input(filename);
        std::string line, output;
        while (std::getline(input, line)) {
            CHECK(line.find("_ORCA_INTERNAL_SECTION") == std::string::npos);
            if (!line.empty() && line.front() != ';') output += line + '\n';
        }
        return output;
    };
    CHECK(commands(plain.string()) == commands(tracked.string()));
}

TEST_CASE("Brick courses publish scheduled entry ordinary and exit envelopes without changing public commands", "[DepositionSectionExport]")
{
    const bool overhang_speed = GENERATE(false, true);
    const bool sloping_roof = GENERATE(false, true);
    CAPTURE(overhang_speed, sloping_roof);
    const auto config = Test::multifilament_config(1, {{"layer_height", "0.2"},
        {"initial_layer_print_height", "0.2"}, {"perimeter_layering", "brick"},
        {"staggered_perimeter_offset", "50%"}, {"staggered_perimeters_inner_only", "1"},
        {"wall_loops", "3"}, {"seam_slope_type", "none"}, {"enable_arc_fitting", "0"},
        {"enable_overhang_speed", overhang_speed ? "1" : "0"}, {"slowdown_for_curled_perimeters", "0"},
        {"enable_prime_tower", "0"}, {"time_lapse_gcode", ""}, {"spiral_mode", "0"},
        {"max_volumetric_extrusion_rate_slope", "0"}});
    Print print;
    Model model;
    // A horizontal roof has whole, constant-Z brick courses. A sloping roof
    // cuts across each inner loop, producing joined raised/unraised fragments
    // and the generator's actual endpoint ramps at their coverage boundary.
    auto mesh = make_cube(8., 8., 2.).its;
    if (sloping_roof)
        for (auto& vertex : mesh.vertices)
            if (vertex.z() > 0.) vertex.z() += .25f * vertex.x();
    Test::init_print({TriangleMesh(mesh)}, print, model, config);
    print.set_status_silent();
    print.process();
    ScopedTemporaryFile plain(".gcode"), tracked(".gcode");
    GCode normal, authenticated;
    GCodeProcessorResult result;
    normal.do_export(&print, plain.string().c_str());
    authenticated.enable_deposition_provenance();
    authenticated.do_export(&print, tracked.string().c_str(), &result);
    INFO(result.deposition_sections_failure);
    REQUIRE(result.deposition_sections);
    REQUIRE(result.deposition_provenance);
    bool entry = false, ordinary = false, exit = false, ramp = false;
    for (const auto& move : result.moves) {
        if (move.internal_only || !move.deposition_command_section_id) continue;
        const auto* command = result.deposition_sections->section(move.deposition_command_section_id);
        REQUIRE(command);
        if (command->transform != DepositionCommandSection::Transform::BrickCourse) continue;
        REQUIRE(move.deposition_source_id != 0);
        REQUIRE(move.deposition_design_heights);
        CHECK_THAT(command->start_height, Catch::Matchers::WithinAbs(command->end_height, 1e-12));
        CHECK_FALSE(command->nominal_lower_z);
        if (command->start_height > .2 + 1e-6) {
            entry = true;
            CHECK_THAT(command->start_height, Catch::Matchers::WithinAbs(.3, 1e-6));
        }
        ordinary |= std::abs(command->start_height - .2) < 1e-6;
        if (command->start_height < .2 - 1e-6) {
            exit = true;
            CHECK_THAT(command->start_height, Catch::Matchers::WithinAbs(.1, 1e-6));
            CHECK(move.extrusion_role != erStaggeredPerimeter);
        }
        ramp |= command->start[2] != command->end[2];
    }
    CHECK(entry);
    CHECK(ordinary);
    CHECK(exit);
    CHECK(ramp == sloping_roof);
    const auto commands = [](const std::string& filename) {
        boost::nowide::ifstream input(filename);
        std::string line, output;
        while (std::getline(input, line)) {
            CHECK(line.find("_ORCA_INTERNAL_SECTION") == std::string::npos);
            if (!line.empty() && line.front() != ';') output += line + '\n';
        }
        return output;
    };
    CHECK(commands(plain.string()) == commands(tracked.string()));
}

TEST_CASE("Live scarf export publishes varying design heights with unchanged public commands", "[DepositionSectionExport]")
{
    const auto config = Test::multifilament_config(1, {{"layer_height", "0.2"},
                                                       {"initial_layer_print_height", "0.2"},
                                                       {"brim_type", "no_brim"},
                                                       {"seam_slope_type", "all"},
                                                       {"seam_slope_conditional", "0"},
                                                       {"enable_prime_tower", "0"},
                                                       {"time_lapse_gcode", ""},
                                                       {"spiral_mode", "0"},
                                                       {"max_volumetric_extrusion_rate_slope", "0.5"}});
    Print print;
    Model model;
    Test::init_print({make_cube(8., 8., 1.)}, print, model, config);
    print.set_status_silent();
    print.process();
    ScopedTemporaryFile output(".gcode");
    GCode generator;
    generator.enable_deposition_provenance();
    GCodeProcessorResult result;
    generator.do_export(&print, output.string().c_str(), &result);
    INFO(result.deposition_sections_failure);
    REQUIRE(result.deposition_sections);
    REQUIRE(result.deposition_sections->size() > 0);
    REQUIRE(result.deposition_provenance);
    bool varying = false;
    for (const auto& move : result.moves)
        if (move.deposition_design_heights)
            varying |= (*move.deposition_design_heights)[0] != (*move.deposition_design_heights)[1];
    CHECK(varying);
    const auto timed = Thermoelastic::timed_deposition(result);
    REQUIRE(timed.state == Thermoelastic::ToolpathState::Complete);
    bool variable_geometry = false;
    std::size_t leading_children=0, tail_children=0;
    for (const auto& entry : timed.deposition) {
        const auto& move = result.moves.at(entry.move_index);
        if (!move.deposition_design_heights) continue;
        REQUIRE(entry.geometry.end_height_m);
        REQUIRE(move.deposition_design_width);
        CHECK_THAT(entry.geometry.height_m, Catch::Matchers::WithinAbs(double((*move.deposition_design_heights)[0]) * .001, 1e-12));
        CHECK_THAT(*entry.geometry.end_height_m, Catch::Matchers::WithinAbs(double((*move.deposition_design_heights)[1]) * .001, 1e-12));
        CHECK_THAT(entry.geometry.width_m, Catch::Matchers::WithinAbs(double(*move.deposition_design_width) * .001, 1e-12));
        variable_geometry |= entry.geometry.height_m != *entry.geometry.end_height_m;
        REQUIRE(move.deposition_source_id!=0);
        const auto* source=result.deposition_provenance->source(move.deposition_source_id);
        REQUIRE(source);
        const auto* layer=print.get_object(0)->get_layer(source->layer);
        REQUIRE(layer);
        const auto native=[](double metres) { return std::llround(metres*1e9); };
        const auto bottom=native(layer->bottom_z()*.001), top=native(layer->print_z*.001);
        INFO("command=" << move.deposition_command_section_id << " source=" << move.deposition_source_id);
        CHECK(native(entry.geometry.start_m.z()-entry.geometry.height_m)>=bottom);
        CHECK(native(entry.geometry.end_m.z()-*entry.geometry.end_height_m)>=bottom);
        CHECK(native(entry.geometry.start_m.z())<=top);
        CHECK(native(entry.geometry.end_m.z())<=top);
        const auto* command=result.deposition_sections->section(move.deposition_command_section_id);
        REQUIRE(command);
        if (command->end[2]>command->start[2]) {
            ++leading_children;
            REQUIRE(command->nominal_lower_z);
            CHECK(native((*command->nominal_lower_z)[0]*.001)==bottom);
            CHECK(native((*command->nominal_lower_z)[1]*.001)==bottom);
            // Leading scarf commands rise from one constant native base.
            // Check every pressure-processed child, not just parent endpoints:
            // independently rounded child Z must not create under-layer beads.
            CHECK(native(entry.geometry.start_m.z()-entry.geometry.height_m)==bottom);
            CHECK(native(entry.geometry.end_m.z()-*entry.geometry.end_height_m)==bottom);
        } else if (command->end[2]==command->start[2] && command->start_height!=command->end_height) {
            ++tail_children;
            // The return pass fills above a rising scarf base, not the layer
            // bottom. It must not receive the leading-ramp base correction.
            CHECK_FALSE(command->nominal_lower_z);
        }
    }
    CHECK(variable_geometry);
    CHECK(leading_children>0);
    CHECK(tail_children>0);
    boost::nowide::ifstream input(output.string());
    const std::string bytes((std::istreambuf_iterator<char>(input)), {});
    CHECK(bytes.find("_ORCA_INTERNAL_SECTION") == std::string::npos);
}

TEST_CASE("Private section scopes publish nominal child heights without leaking to exported Gcode", "[DepositionSectionExport]")
{
    const bool offset = GENERATE(false, true);
    auto sections = std::make_shared<DepositionSections>(42);
    const auto id = sections->add({{.123, .321, .2}, {10.123, .321, .2}, .4, .2, 0., 0});
    GCodeProcessor processor;
    processor.apply_config(FullPrintConfig{});
    if (offset) processor.set_xy_offset(100., 200.);
    ScopedTemporaryFile output(".gcode");
    SectionStreamAccess::GCodeOutputStream stream(boost::nowide::fopen(output.string().c_str(), "wb"), processor);
    REQUIRE(stream.is_open());
    stream.enable_provenance(42);
    const std::string header = "G90\nM83\n;" + GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Role) +
        "Outer wall\nG1 X0.123 Y0.321 Z0.2 F600\n";
    stream.write(header);
    // A quantized stationary scarf command remains in the file and processor
    // inventory, but has no invented spatial section or scope.
    stream.write("G1 E0.00001\n");
    stream.write(sections->marker(id));
    stream.write("G1 X4.123 E0.04\nM400\nG1 X10.123 E0.06\n");
    stream.write(sections->marker(0));
    sections->seal();
    REQUIRE(stream.finish_sections(*sections));
    processor.finish_deposition_sections(sections);
    stream.close();
    REQUIRE(processor.result().deposition_sections);
    const auto* published = processor.result().deposition_sections->section(id);
    REQUIRE(published);
    CHECK_THAT(published->start[0], Catch::Matchers::WithinAbs(float((offset ? 100. : 0.) + .123), 1e-7));
    CHECK_THAT(published->end[1], Catch::Matchers::WithinAbs(float((offset ? 200. : 0.) + .321), 1e-7));
    std::size_t count = 0;
    for (const auto& move : processor.result().moves) {
        if (move.internal_only || !move.deposition_design_heights)
            continue;
        if (count == 0)
            CHECK_THAT((*move.deposition_design_heights)[1], Catch::Matchers::WithinAbs(.12, 1e-7));
        ++count;
    }
    CHECK(count == 2);
    boost::nowide::ifstream input(output.string());
    const std::string bytes((std::istreambuf_iterator<char>(input)), {});
    CHECK(bytes == header + "G1 E0.00001\nG1 X4.123 E0.04\nM400\nG1 X10.123 E0.06\n");
}

TEST_CASE("A straight section following an arc uses the same parser coordinate frame", "[DepositionSectionExport]")
{
    const auto flavor = GENERATE(gcfMarlinFirmware, gcfKlipper);
    const bool displaced = GENERATE(false, true);
    const auto coordinate = [&](double value) { return flavor == gcfKlipper ? value : double(float(value)); };
    auto sections = std::make_shared<DepositionSections>(42);
    const auto id = sections->add({{coordinate(126.854), coordinate(46.654), coordinate(.2)},
        {coordinate(126.844), coordinate(46.936), coordinate(.2)}, .4, .2, .2, 0});
    GCodeProcessor processor;
    FullPrintConfig config;
    config.gcode_flavor.value = flavor;
    processor.apply_config(config);
    ScopedTemporaryFile output(".gcode");
    SectionStreamAccess::GCodeOutputStream stream(boost::nowide::fopen(output.string().c_str(), "wb"), processor);
    REQUIRE(stream.is_open());
    stream.enable_provenance(42);
    stream.write("G90\nM83\n;" + GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Role) +
        "Outer wall\nG1 X125.854 Y46.654 Z0.2 F600\nG2 X126.854 Y46.654 I0.5 J0 E0.1\n");
    stream.write(sections->marker(id));
    stream.write(displaced ? "G1 X126.9 Y46.936 E0.01\n" : "G1 X126.844 Y46.936 E0.01\n");
    stream.write(sections->marker(0));
    sections->seal();
    REQUIRE(stream.finish_sections(*sections));
    processor.finish_deposition_sections(sections);
    stream.close();
    INFO(processor.result().deposition_sections_failure);
    if (displaced) {
        CHECK_FALSE(processor.result().deposition_sections);
        return;
    }
    REQUIRE(processor.result().deposition_sections);
    REQUIRE_FALSE(processor.result().moves.empty());
    CHECK(processor.result().moves.back().deposition_command_section_id == id);
}

TEST_CASE("Anchored Z evidence preserves zero-height starts beyond float preview precision", "[DepositionSectionExport]")
{
    const double base = GENERATE(1.4, 20.2);
    const bool relative = GENERATE(false, true);
    const double offset = GENERATE(0., .125);
    const bool collapsed = GENERATE(false, true);
    const double rise = collapsed ? 0. : .001;
    CAPTURE(base, relative, offset, collapsed);
    DepositionCommandSection parent{{2., 0., base}, {3., 0., base + rise}, .4, 0., rise, 0};
    parent.nominal_lower_z = std::array<double, 2>{base, base};
    parent.end_height = parent.end[2] - base;
    auto sections = std::make_shared<DepositionSections>(42);
    const auto id = sections->add(parent);
    GCodeProcessor processor;
    // G92 synchronizes pending timing blocks, whose first real move needs the
    // same initial sentinel vertex that live export initializes.
    processor.initialize_result_moves();
    FullPrintConfig config;
    config.gcode_flavor.value = gcfMarlinFirmware;
    config.z_offset.value = offset;
    processor.apply_config(config);
    ScopedTemporaryFile output(".gcode");
    SectionStreamAccess::GCodeOutputStream stream(boost::nowide::fopen(output.string().c_str(), "wb"), processor);
    REQUIRE(stream.is_open());
    stream.enable_provenance(42);
    std::string header = "G90\nM83\n;" + GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Role) +
        "Outer wall\nG1 X1 Y0 Z" + std::to_string(base) + " F600\nG1 X2 Y0\n";
    if (relative) header += "G92 Z0\nG91\n";
    stream.write(header);
    stream.write(sections->marker(id));
    const std::string emitted = relative ? "G1 X1 Z" + std::to_string(rise) + " E0.01\n" :
        "G1 X3 Z" + std::to_string(base + rise) + " E0.01\n";
    stream.write(emitted);
    stream.write(sections->marker(0));
    sections->seal();
    REQUIRE(stream.finish_sections(*sections));
    processor.finish_deposition_sections(sections);
    stream.close();
    const auto& result = processor.result();
    INFO(result.deposition_sections_failure);
    REQUIRE(result.deposition_sections);
    REQUIRE(result.deposition_z_spans.size() == 1);
    CHECK(result.deposition_z_generation == 42);
    const auto& span = result.deposition_z_spans.front();
    CHECK(result.deposition_z_for_move(span.move_index) == &span);
    CHECK(span.command_id == id);
    CHECK_THAT(span.z_mm[0], Catch::Matchers::WithinAbs(base - offset, 1e-12));
    CHECK_THAT(span.z_mm[1], Catch::Matchers::WithinAbs(base + rise - offset, 1e-12));
    CHECK((span.preview_z_mm[0] == float(base - offset)));
    const auto& move = result.moves.at(span.move_index);
    REQUIRE(move.deposition_design_heights);
    CHECK_THAT((*move.deposition_design_heights)[0], Catch::Matchers::WithinAbs(0., 1e-12));
    CHECK_THAT(span.z_mm[1] - (*move.deposition_design_heights)[1], Catch::Matchers::WithinAbs(base - offset, 1e-12));
    if (collapsed) CHECK(((*move.deposition_design_heights)[1] == 0.));
    boost::nowide::ifstream input(output.string());
    const std::string bytes((std::istreambuf_iterator<char>(input)), {});
    CHECK(bytes == header + emitted);
}
