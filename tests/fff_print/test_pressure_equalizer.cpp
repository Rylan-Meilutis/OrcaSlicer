#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/GCode.hpp"
#include "libslic3r/GCode/PressureEqualizer.hpp"
#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <algorithm>
#include <cstddef>
#include <map>

using namespace Slic3r;

TEST_CASE("Unsmoothed arcs retain their input modal feed after linear smoothing", "[PressureEqualizer]")
{
    const bool full_circle = GENERATE(false, true);
    GCodeConfig config;
    config.filament_diameter.values = {1.75};
    config.use_relative_e_distances.value = true;
    config.max_volumetric_extrusion_rate_slope.value = .5;
    config.max_volumetric_extrusion_rate_slope_segment_length.value = 1.;
    config.extrusion_rate_smoothing_external_perimeter_only.value = false;
    std::string input = "G1 X0 Y0 Z0.2 F600\n;_EXTRUSION_ROLE:" + std::to_string(int(erInternalInfill)) + "\n";
    for (int i = 1; i <= 80; ++i) {
        input += i == 80 ? ";_EXTRUDE_SET_SPEED\n" :
            "G1 F" + std::to_string(i % 3 == 0 ? 600 : 6000) + " ;_EXTRUDE_SET_SPEED\n";
        input += "G1 X" + std::to_string(i * 4) + " Y0 E0.4\n;_EXTRUDE_END\n";
    }
    const std::string arc = full_circle ? "G2 X320 Y0 I0 J10 E0.5" : "G2 X330 Y10 I0 J10 E0.5";
    input += arc + '\n';
    PressureEqualizer filter(config);
    CHECK(filter.process_layer(LayerResult{input, 0}).nop_layer_result);
    const auto output = filter.process_layer(LayerResult::make_nop_layer_result()).gcode;
    REQUIRE(output.find(arc + '\n') != std::string::npos);
    GCodeReader reader;
    double feed = 0.;
    bool found_arc = false, reduced_feed = false;
    reader.parse_buffer(output, [&](GCodeReader&, const GCodeReader::GCodeLine& line) {
        if (line.has_f()) feed = line.f();
        if (line.cmd_is("G1") && line.has_e() && feed > 600. && feed < 6000.) reduced_feed = true;
        if (line.cmd_is("G2")) {
            found_arc = true;
            CHECK_THAT(feed, Catch::Matchers::WithinAbs(6000., 1e-6));
        }
    });
    CHECK(found_arc);
    CHECK(reduced_feed);
}

TEST_CASE("Pressure smoothing starts linear subdivisions at the actual arc endpoint", "[PressureEqualizer]")
{
    const bool relative_e = GENERATE(false, true);
    const std::string arc_command = GENERATE("G2", "G3");
    GCodeConfig config;
    config.filament_diameter.values = {1.75};
    config.use_relative_e_distances.value = relative_e;
    config.max_volumetric_extrusion_rate_slope.value = .5;
    config.max_volumetric_extrusion_rate_slope_segment_length.value = 1.;
    config.extrusion_rate_smoothing_external_perimeter_only.value = false;
    const std::string arc = arc_command + " X10 Y10 Z0.4 I0 J10 E0.5 F6000";
    std::string input = "G1 X0 Y0 Z0.2 F600\n;_EXTRUSION_ROLE:" + std::to_string(int(erInternalInfill)) +
        "\n" + arc + "\n";
    constexpr int lines = 80;
    for (int i = 1; i <= lines; ++i) {
        // The first line inherits the arc's feed; later lines force genuine
        // smoothing and subdivision rather than a pass-through-only test.
        input += i == 1 ? ";_EXTRUDE_SET_SPEED\n" :
            "G1 F" + std::to_string(i % 3 == 0 ? 600 : 6000) + " ;_EXTRUDE_SET_SPEED\n";
        input += "G1 X" + std::to_string(10 + i * 4) + " Y10 E" +
            std::to_string(relative_e ? .4 : .5 + i * .4) + "\n;_EXTRUDE_END\n";
    }
    PressureEqualizer filter(config);
    CHECK(filter.process_layer(LayerResult{input, 0}).nop_layer_result);
    const auto output = filter.process_layer(LayerResult::make_nop_layer_result()).gcode;
    REQUIRE(output.find(arc + '\n') != std::string::npos);
    GCodeReader reader;
    bool after_arc = false;
    double y = 0., z = 0., extrusion = 0.;
    std::size_t subdivisions = 0;
    reader.parse_buffer(output, [&](GCodeReader&, const GCodeReader::GCodeLine& line) {
        if (line.has_y()) y = line.y();
        if (line.has_z()) z = line.z();
        if (line.has_e()) extrusion = relative_e ? extrusion + line.e() : line.e();
        if (line.cmd_is("G2") || line.cmd_is("G3")) { after_arc = true; return; }
        if (!after_arc || !line.cmd_is("G1") || !line.has_e()) return;
        ++subdivisions;
        CHECK_THAT(y, Catch::Matchers::WithinAbs(10., .0006));
        CHECK_THAT(z, Catch::Matchers::WithinAbs(.4, 1e-6));
        if (line.has_x()) {
            CHECK(line.x() >= 10.);
            CHECK(line.x() <= 10. + lines * 4.);
        }
    });
    CHECK(subdivisions > lines);
    CHECK_THAT(extrusion, Catch::Matchers::WithinAbs(.5 + lines * .4, .001));
}

TEST_CASE("Pressure smoothing preserves private boundaries without changing toolpaths", "[PressureEqualizer]")
{
    GCodeConfig config;
    config.filament_diameter.values = {1.75};
    config.use_relative_e_distances.value = true;
    config.max_volumetric_extrusion_rate_slope.value = .5;
    config.max_volumetric_extrusion_rate_slope_segment_length.value = 1.;
    config.extrusion_rate_smoothing_external_perimeter_only.value = false;
    const std::string prefix = GENERATE(";_ORCA_INTERNAL_SOURCE:1:", ";_ORCA_INTERNAL_SECTION:1:");
    const std::string begin = prefix + "1\n";
    const std::string end = prefix + "0\n";
    std::vector<std::string> layers;
    layers.push_back("G1 X0 Y0 Z0.2 F600\n;_EXTRUSION_ROLE:" + std::to_string(int(erExternalPerimeter)) + "\n");
    for (int i = 1; i <= 180; ++i) {
        layers.back() += "G1 F" + std::to_string(i%3 == 0 ? 600 : 6000) + " ;_EXTRUDE_SET_SPEED\n";
        layers.back() += "G1 X" + std::to_string(i*4) + " Y0 E0.4\n;_EXTRUDE_END\n";
    }
    layers.push_back(""); // A marker-only layer still has to flush its tail.
    layers.push_back("G1 X724 Y0 E0.4 F600\n");
    std::map<unsigned,double> source_endpoints;
    const auto run = [&](bool marked) {
        PressureEqualizer filter(config);
        std::string output;
        unsigned source = 2;
        for (std::size_t i = 0; i < layers.size(); ++i) {
            std::string input;
            if (marked) input += begin + end;
            std::istringstream stream(layers[i]);
            std::string line;
            while (std::getline(stream,line)) {
                if (marked) {
                    input += prefix + std::to_string(source) + '\n';
                    if (line.find("G1 X") == 0 && line.find('E') != std::string::npos)
                        source_endpoints[source] = std::stod(line.substr(4));
                    ++source;
                }
                input += line + '\n';
                if (marked) input += end;
            }
            if (marked) input += begin + end;
            output += filter.process_layer(LayerResult{input,i}).gcode;
        }
        output += filter.process_layer(LayerResult::make_nop_layer_result()).gcode;
        return output;
    };
    const auto plain = run(false), marked = run(true);
    REQUIRE_FALSE(plain.empty());
    std::istringstream plain_stream(plain);
    std::string plain_line;
    std::size_t extrusion_lines = 0;
    while (std::getline(plain_stream,plain_line))
        if (plain_line.find("G1 X") == 0 && plain_line.find('E') != std::string::npos) ++extrusion_lines;
    // Exercise actual line subdivisions, not only a pass-through filter.
    REQUIRE(extrusion_lines > 181);
    std::istringstream stream(marked);
    std::string stripped, line;
    unsigned active = 0;
    std::size_t scopes = 0;
    while (std::getline(stream,line)) {
        if (line+'\n' == end) { CHECK(active != 0); active = 0; }
        else if (line.find(prefix) == 0) {
            CHECK(active == 0);
            active = unsigned(std::stoul(line.substr(line.rfind(':')+1)));
            ++scopes;
        } else {
            if (line.find("G1 X") == 0 && line.find('E') != std::string::npos) {
                const auto owner = source_endpoints.find(active);
                REQUIRE(owner != source_endpoints.end());
                const double x = std::stod(line.substr(4));
                CHECK(x > owner->second-4.);
                CHECK(x <= owner->second);
            }
            stripped += line+'\n';
        }
    }
    CHECK(active == 0);
    std::size_t expected_scopes = 2*layers.size();
    for (const auto &layer : layers)
        expected_scopes += std::count(layer.begin(),layer.end(),'\n');
    CHECK(scopes == expected_scopes);
    CHECK(stripped == plain);

    // A source scope can cross the filter's delayed layer-output boundary.
    const auto spanning = [&](bool marked) {
        PressureEqualizer filter(config);
        const auto first = filter.process_layer(LayerResult{
            (marked ? begin : "") + std::string("G1 X0 Y0 Z0.2 F600\nG1 X4 Y0 E0.4\n"),0});
        CHECK(first.nop_layer_result);
        std::string output = filter.process_layer(LayerResult{
            std::string("G1 X8 Y0 E0.4\n") + (marked ? end : ""),1}).gcode;
        output += filter.process_layer(LayerResult::make_nop_layer_result()).gcode;
        return output;
    };
    const auto spanning_plain = spanning(false), spanning_marked = spanning(true);
    REQUIRE(spanning_marked.find(begin) == 0);
    REQUIRE(spanning_marked.size() >= begin.size()+end.size());
    CHECK(spanning_marked.substr(spanning_marked.size()-end.size()) == end);
    CHECK(spanning_marked.substr(begin.size(),spanning_marked.size()-begin.size()-end.size()) == spanning_plain);
}
