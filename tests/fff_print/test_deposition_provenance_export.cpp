#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/GCode/DepositionProvenance.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Config.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/GCode.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/ThermoelasticSourceOwnership.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Surface.hpp"
#include "libslic3r/SurfaceCollection.hpp"
#include "libslic3r/libslic3r.h"
#include <Eigen/Core>
#include "libslic3r/TriangleMesh.hpp"
#include "test_helpers.hpp"
#include "test_utils.hpp"
#include <boost/nowide/fstream.hpp>
#include <boost/nowide/cstdio.hpp>
#include <cstddef>
#include <memory>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <utility>
#include <cstdint>

using namespace Slic3r;

namespace {
class ProvenanceStreamAccess : public GCode {
public:
    using GCode::GCodeOutputStream;
};

std::string tag(GCodeProcessor::ETags type, const std::string &value)
{
    return ";" + GCodeProcessor::reserved_tag(type) + value + "\n";
}

std::string planar_header()
{
    return "G90\nM83\n" + tag(GCodeProcessor::ETags::Role, "Outer wall") +
        tag(GCodeProcessor::ETags::Width, "0.4") + tag(GCodeProcessor::ETags::Height, "0.2") +
        "G1 X1 Y1 Z0.2 F600\n";
}
}

TEST_CASE("Export stream strips fragmented source controls without changing public line numbering", "[DepositionProvenanceExport]")
{
    auto registry = std::make_shared<DepositionProvenance>();
    const auto id = registry->add({});
    GCodeProcessor processor;
    processor.apply_config(FullPrintConfig{});
    ScopedTemporaryFile output(".gcode");
    ProvenanceStreamAccess::GCodeOutputStream stream(boost::nowide::fopen(output.string().c_str(),"wb"),processor);
    REQUIRE(stream.is_open());
    stream.enable_provenance(registry->generation());
    const std::string header = planar_header();
    stream.write(header);
    const auto begin = registry->begin(id), end = registry->clear();
    stream.write(begin.substr(0,5));
    stream.write(begin.substr(5));
    stream.write("G1 X2");
    stream.write(" Y1 E0.04\n");
    stream.write(end.substr(0,3));
    stream.write(end.substr(3));
    stream.write("; unterminated ordinary tail");
    registry->seal();
    REQUIRE(stream.finish_provenance(registry,{id,0}));
    processor.finish_deposition_provenance(registry);
    stream.close();
    boost::nowide::ifstream file(output.string());
    const std::string actual{std::istreambuf_iterator<char>(file),std::istreambuf_iterator<char>()};
    const std::string expected = header+"G1 X2 Y1 E0.04\n; unterminated ordinary tail";
    CHECK(actual == expected);
    GCodeProcessor baseline;
    baseline.apply_config(FullPrintConfig{});
    baseline.process_buffer(expected);
    REQUIRE(processor.get_result().moves.size() == baseline.get_result().moves.size());
    bool attributed = false;
    for (std::size_t i = 0; i < processor.get_result().moves.size(); ++i) {
        const auto &move = processor.get_result().moves[i];
        INFO("line=" << move.gcode_id << " role=" << int(move.extrusion_role)
             << " filament=" << int(move.extruder_id) << " source=" << move.deposition_source_id);
        CHECK(move.gcode_id == baseline.get_result().moves[i].gcode_id);
        if (move.type == EMoveType::Extrude) {
            CHECK(move.deposition_source_id == id);
            attributed = true;
        }
    }
    CHECK(attributed);
}

TEST_CASE("Processor source IDs require a sealed matching logical filament registry", "[DepositionProvenanceExport]")
{
    const int invalid = GENERATE(0,1,2,3,4);
    auto registry = std::make_shared<DepositionProvenance>();
    DepositionSource source;
    source.logical_filament = invalid == 2 ? 1 : 0;
    const auto id = registry->add(source);
    if (invalid != 3) registry->seal();
    GCodeProcessor processor;
    processor.apply_config(FullPrintConfig{});
    processor.process_buffer(planar_header() + ";_ORCA_INTERNAL_SOURCE:1:1\nG1 X2 Y1 E0.04\n");
    // Imported comments cannot activate tracking.
    for (const auto &move : processor.get_result().moves) CHECK(move.deposition_source_id == 0);
    processor.set_deposition_source(invalid == 4 ? id+1 : id);
    processor.process_buffer("G1 X3 Y1 E0.04\nG1 X4 Y1\n" +
        tag(GCodeProcessor::ETags::Role, "Custom") + "G1 X5 Y1 E0.04\n");
    std::size_t attributed = 0;
    for (const auto &move : processor.get_result().moves) {
        if (move.deposition_source_id != 0) {
            ++attributed;
            CHECK(move.type == EMoveType::Extrude);
            CHECK(move.extrusion_role == erExternalPerimeter);
        }
    }
    REQUIRE(attributed == 1);
    CHECK_FALSE(processor.get_result().deposition_provenance);
    processor.finish_deposition_provenance(invalid == 1 ? nullptr : registry);
    CHECK(bool(processor.get_result().deposition_provenance) == (invalid == 0));
    std::size_t published = 0;
    for (const auto &move : processor.get_result().moves) published += move.deposition_source_id != 0;
    CHECK(published == (invalid == 0 ? 1 : 0));
    GCodeProcessorResult copied, moved;
    copied = processor.get_result();
    CHECK(copied.deposition_provenance == processor.get_result().deposition_provenance);
    moved = std::move(copied);
    CHECK(moved.deposition_provenance == processor.get_result().deposition_provenance);
    CHECK(moved.moves.size() == processor.get_result().moves.size());
    processor.reset();
    CHECK_FALSE(processor.get_result().deposition_provenance);
}

TEST_CASE("Live exports retain native source identity without private file markers", "[DepositionProvenanceExport]")
{
    const bool offsets = GENERATE(false, true);
    auto config = Test::multifilament_config(2, {
        {"layer_height","0.2"}, {"initial_layer_print_height","0.2"},
        {"enable_prime_tower","0"}, {"time_lapse_gcode",""},
        {"spiral_mode","0"}, {"gcode_label_objects","0"},
        {"thick_bridges","0"}, {"thick_internal_bridges","0"},
        {"max_volumetric_extrusion_rate_slope","0.5"}});
    if (offsets)
        config.set_key_value("extruder_offset", new ConfigOptionPoints{Vec2d(5., -3.)});
    const std::vector<std::vector<ConfigBase::SetDeserializeItem>> overrides{
        {{"extruder","1"}},{{"extruder","2"}}};
    const auto run = [&](bool tracked) {
        Print print;
        Model model;
        std::vector<TriangleMesh> meshes{make_cube(4.,4.,1.),make_cube(4.,4.,1.)};
        Test::init_print(std::move(meshes),print,model,config,&overrides);
        print.set_status_silent();
        print.process();
        ScopedTemporaryFile output(".gcode");
        GCode generator;
        if (offsets) generator.set_gcode_offset(100., 200.);
        generator.enable_deposition_provenance(tracked);
        GCodeProcessorResult result;
        generator.do_export(&print,output.string().c_str(),&result);
        boost::nowide::ifstream file(output.string());
        const std::string text{std::istreambuf_iterator<char>(file),std::istreambuf_iterator<char>()};
        REQUIRE_FALSE(text.empty());
        CHECK(text.find("_ORCA_INTERNAL_SOURCE") == std::string::npos);
        CHECK(bool(result.deposition_provenance) == tracked);
        CHECK(bool(result.native_ownership) == tracked);
        GCodeProcessorResult delivered;
        delivered = result;
        CHECK(delivered.native_ownership == result.native_ownership);
        GCodeProcessorResult moved;
        moved = std::move(delivered);
        CHECK(moved.native_ownership == result.native_ownership);
        moved.reset();
        CHECK_FALSE(moved.native_ownership);
        CHECK_FALSE(moved.deposition_provenance);
        std::set<unsigned> filaments;
        std::set<std::uint64_t> objects;
        std::size_t sourced = 0;
        for (const auto &move : result.moves) {
            if (!tracked) CHECK(move.deposition_source_id == 0);
            if (tracked && move.type == EMoveType::Extrude &&
                (move.extrusion_role == erPerimeter || move.extrusion_role == erExternalPerimeter ||
                 move.extrusion_role == erInternalInfill || move.extrusion_role == erSolidInfill ||
                 move.extrusion_role == erTopSolidInfill))
                CHECK(move.deposition_source_id != 0);
            if (move.deposition_source_id == 0) continue;
            ++sourced;
            REQUIRE(result.deposition_provenance);
            const auto *source = result.deposition_provenance->source(move.deposition_source_id);
            REQUIRE(source);
            REQUIRE(result.native_ownership);
            CHECK(result.native_ownership->source_status(move.deposition_source_id) ==
                  Thermoelastic::NativeOwnershipSnapshot::SourceStatus::NativePlanar);
            CHECK(source->logical_filament == move.extruder_id);
            CHECK(source->height == DepositionSource::Height::NativePlanar);
            CHECK(source->assignment == DepositionSource::Assignment::NominalRole);
            CHECK(source->instance == 0);
            CHECK(source->layer < print.get_object(0)->layer_count());
            const PrintObject *owner = nullptr;
            for (const auto *object : print.objects())
                if (object->get_id() == source->object) owner = object;
            REQUIRE(owner);
            BoundingBox native_bounds;
            for (const auto *region : owner->get_layer(source->layer)->regions())
                if (region->region().print_region_id() == source->region)
                    for (const auto &surface : region->slices.surfaces)
                        native_bounds.merge(get_extents(surface.expolygon));
            native_bounds.translate(owner->instances()[source->instance].shift);
            CHECK(native_bounds.contains(Point(scaled(move.position.x()), scaled(move.position.y()))));
            filaments.insert(source->logical_filament);
            objects.insert(source->object);
        }
        if (tracked) {
            CHECK(sourced > 0);
            CHECK(filaments == std::set<unsigned>{0,1});
            CHECK(objects.size() == 2);
        }
        // Timestamps and model IDs may differ; every actual machine command
        // must remain byte-identical, including all feedrates and tool changes.
        std::istringstream lines(text);
        std::string commands,line;
        while (std::getline(lines,line)) {
            line.resize(line.find(';') == std::string::npos ? line.size() : line.find(';'));
            if (line.find_first_not_of(" \t\r") != std::string::npos) commands += line+'\n';
        }
        return commands;
    };
    const auto plain = run(false);
    CHECK(run(true) == plain);
}
