#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "libslic3r/GCode.hpp"
#include "libslic3r/GCode/AdaptivePAProcessor.hpp"
#include "libslic3r/PrintConfig.hpp"
#include <sstream>
#include <string>
#include <utility>

using namespace Slic3r;

TEST_CASE("Adaptive pressure advance ignores private boundaries in feedrate lookahead", "[AdaptivePAProcessor]")
{
    PrintConfig config;
    config.enable_pressure_advance.values = {true};
    config.adaptive_pressure_advance.values = {true};
    config.gcode_comments.value = true;
    config.adaptive_pressure_advance_model.values = {"0.02,1,1000\n0.04,10,1000\n0.03,1,2000\n0.05,10,2000"};
    const std::string tag = "; PA_CHANGE:T0 MM3MM:0.08 ACCEL:1000 BR:0 RC:1 OV:0\n";
    const std::string prefix = GENERATE(";_ORCA_INTERNAL_SOURCE:1:", ";_ORCA_INTERNAL_SECTION:1:");
    const std::string begin = prefix + "1\n", end = prefix + "0\n";
    const auto run = [&](bool marked) {
        GCode generator;
        generator.apply_print_config(config);
        AdaptivePAProcessor filter(generator);
        std::string input = "G1 F6000\n" + tag;
        if (marked) input += begin;
        input += "G1 F600\nG1 X1 Y1 E0.1\nG1 X2 Y1 E0.1\n";
        if (marked) input += end;
        input += "G1 X3 Y1\n";
        return filter.process_layer(std::move(input));
    };
    const auto plain = run(false), marked = run(true);
    REQUIRE(plain.find("; APA Current Speed: 10.000000") != std::string::npos);
    std::string stripped,line;
    std::istringstream stream(marked);
    unsigned boundaries = 0;
    while (std::getline(stream,line)) {
        if (line+'\n' == begin || line+'\n' == end) ++boundaries;
        else stripped += line+'\n';
    }
    CHECK(boundaries == 2);
    CHECK(stripped == plain);
}
