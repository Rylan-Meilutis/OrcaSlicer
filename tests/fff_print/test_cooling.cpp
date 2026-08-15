#include <catch2/catch_all.hpp>

#include "test_helpers.hpp"
#include "libslic3r/GCode/CoolingBuffer.hpp"

#include <string>

using namespace Slic3r;
using namespace Slic3r::Test;

// The fan is held off for the first close_fan_the_first_x_layers layers, so an explicit
// fan-off command is emitted.
TEST_CASE("Fan is held off for the initial layers", "[Cooling]")
{
    const std::string gcode = slice({ cube(20) }, {
        { "cooling",                      true },
        { "close_fan_the_first_x_layers", 5 },
    });
    CHECK(gcode.find("M106 S0") != std::string::npos);
}

// The cooling pass resolves and strips its internal speed placeholders; none leak into
// the final G-code.
TEST_CASE("Cooling consumes its internal speed markers", "[Cooling]")
{
    const std::string gcode = slice({ cube(20) }, { { "layer_height", 0.2 } });
    CHECK(gcode.find(";_EXTRUDE_SET_SPEED") == std::string::npos);
}

TEST_CASE("Hull-line mitigation limits abrupt layer-time reductions", "[Cooling][HullLine]")
{
    CHECK_THAT(CoolingBuffer::hull_line_target_time(40.f, 20.f, 25.f),
               Catch::Matchers::WithinAbs(32.f, EPSILON));
    CHECK_THAT(CoolingBuffer::hull_line_target_time(40.f, 35.f, 25.f),
               Catch::Matchers::WithinAbs(35.f, EPSILON));
    CHECK_THAT(CoolingBuffer::hull_line_target_time(0.f, 20.f, 25.f),
               Catch::Matchers::WithinAbs(20.f, EPSILON));
}
