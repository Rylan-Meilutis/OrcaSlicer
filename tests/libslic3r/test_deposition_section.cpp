#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include "libslic3r/GCode/DepositionSection.hpp"
#include <array>
#include <cstddef>
#include <stdexcept>

using namespace Slic3r;

TEST_CASE("Deposition sections retain authenticated transforms and reject unknown values", "[DepositionSection]")
{
    DepositionSections registry(42);
    DepositionCommandSection command{{0., 0., .3}, {1., 0., .3}, .4, .2, .2, 0};
    CHECK(command.transform == DepositionCommandSection::Transform::Ordinary);
    command.transform = DepositionCommandSection::Transform::BrickCourse;
    const auto id = registry.add(command);
    REQUIRE(registry.section(id));
    CHECK(registry.section(id)->transform == DepositionCommandSection::Transform::BrickCourse);
    command.transform = static_cast<DepositionCommandSection::Transform>(255);
    CHECK_THROWS_AS(registry.add(command), std::invalid_argument);
    CHECK(registry.size() == 1);
}

TEST_CASE("Auxiliary brim sections require an explicit valid purpose and conventional envelope", "[DepositionSection]")
{
    DepositionSections registry(42);
    DepositionCommandSection command{{0., 0., .2}, {1., 0., .2}, .4, .2, .2, 0};
    CHECK(command.purpose == DepositionCommandSection::Purpose::Model);
    command.purpose = DepositionCommandSection::Purpose::AuxiliaryBrim;
    const auto id = registry.add(command);
    REQUIRE(registry.section(id));
    CHECK(registry.section(id)->purpose == DepositionCommandSection::Purpose::AuxiliaryBrim);
    auto invalid = command;
    invalid.purpose = static_cast<DepositionCommandSection::Purpose>(255);
    CHECK_THROWS_AS(registry.add(invalid), std::invalid_argument);
    invalid = command;
    invalid.transform = DepositionCommandSection::Transform::BrickCourse;
    CHECK_THROWS_AS(registry.add(invalid), std::invalid_argument);
    invalid = command;
    invalid.end[2] += .1;
    CHECK_THROWS_AS(registry.add(invalid), std::invalid_argument);
    invalid = command;
    invalid.end_height = .1;
    CHECK_THROWS_AS(registry.add(invalid), std::invalid_argument);
    invalid = command;
    invalid.start_height = invalid.end_height = 0.;
    CHECK_THROWS_AS(registry.add(invalid), std::invalid_argument);
    invalid = command;
    invalid.nominal_lower_z = std::array<double, 2>{0., 0.};
    CHECK_THROWS_AS(registry.add(invalid), std::invalid_argument);
    CHECK(registry.size() == 1);
}

TEST_CASE("Explicit scarf bases survive rounded child tops without inventing thickness", "[DepositionSection]")
{
    DepositionCommandSection parent{{0., 0., 3.816}, {10., 0., 3.820}, .4, 0., 0., 0};
    parent.nominal_lower_z = std::array<double, 2>{3.8, 3.8};
    parent.start_height = parent.start[2] - 3.8;
    parent.end_height = parent.end[2] - 3.8;
    DepositionSections registry(42);
    REQUIRE(registry.add(parent) == 1);
    DepositionSectionChain chain(parent);
    const std::array<double, 3> rounded_child{5., 0., 3.8179};
    const auto first = chain.append(parent.start, rounded_child, 0);
    REQUIRE(first);
    CHECK_THAT(rounded_child[2] - (*first)[1], Catch::Matchers::WithinAbs(3.8, 1e-12));
    const auto second = chain.append(rounded_child, parent.end, 0);
    REQUIRE(second);
    CHECK_THAT(rounded_child[2] - (*second)[0], Catch::Matchers::WithinAbs(3.8, 1e-12));
    CHECK(chain.complete());
    auto inconsistent = parent;
    inconsistent.start_height += .001;
    CHECK_THROWS_AS(registry.add(inconsistent), std::invalid_argument);
    auto collapsed = parent;
    collapsed.nominal_lower_z = std::array<double, 2>{parent.start[2], parent.end[2]};
    collapsed.start_height = collapsed.end_height = 0.;
    REQUIRE(registry.add(collapsed) == 2);
    DepositionSectionChain collapsed_chain(collapsed);
    const auto empty_envelope = collapsed_chain.append(collapsed.start, collapsed.end, 0);
    REQUIRE(empty_envelope);
    CHECK(((*empty_envelope)[0] == 0.));
    CHECK(((*empty_envelope)[1] == 0.));
    CHECK(collapsed_chain.complete());
    collapsed.nominal_lower_z.reset();
    CHECK_THROWS_AS(registry.add(collapsed), std::invalid_argument);
    auto negative = parent;
    negative.nominal_lower_z = std::array<double, 2>{parent.start[2] + .001, 3.8};
    negative.start_height = -.001;
    CHECK_THROWS_AS(registry.add(negative), std::invalid_argument);
}

TEST_CASE("Circular command evidence checks direction radius and ordered children", "[DepositionSection]")
{
    for (int direction : {-1,1}) {
        DepositionCommandSection parent{{1.,0.,.2},{-1.,0.,.2},.4,.2,.2,0};
        parent.arc_center = {0.,0.};
        parent.arc_direction = direction;
        DepositionSections registry(1);
        REQUIRE(registry.add(parent) == 1);
        const std::array<double,3> midpoint{0.,double(direction),.2};
        DepositionSectionChain valid(parent), opposite(parent), radial(parent), skipped(parent);
        REQUIRE(valid.append(parent.start,midpoint,0));
        CHECK_FALSE(valid.complete());
        REQUIRE(valid.append(midpoint,parent.end,0));
        CHECK(valid.complete());
        CHECK_FALSE(opposite.append(parent.start,{0.,-double(direction),.2},0));
        CHECK_FALSE(radial.append(parent.start,{0.,1.01*direction,.2},0));
        CHECK_FALSE(skipped.append(midpoint,parent.end,0));
    }
}

TEST_CASE("Deposition command storage rejects work beyond its finite parent budget", "[DepositionSection]")
{
    DepositionSections registry(1);
    const DepositionCommandSection command{{0., 0., .2}, {1., 0., .2}, .4, .2, .2, 0};
    for (std::size_t i = 0; i < DepositionSections::maximum_commands; ++i)
        registry.add(command);
    CHECK(registry.size() == DepositionSections::maximum_commands);
    CHECK_THROWS_AS(registry.add(command), std::length_error);
    CHECK(registry.size() == DepositionSections::maximum_commands);
}

TEST_CASE("Stationary commands cannot acquire a nominal scarf section", "[DepositionSection]")
{
    DepositionSections registry(1);
    CHECK_THROWS_AS(registry.add({{1., 2., .2}, {1., 2., .2}, .4, .1, .2, 0}), std::invalid_argument);
    CHECK(registry.size() == 0);
    CHECK(registry.add({{1., 2., .2}, {2., 2., .2}, .4, .1, .2, 0}) == 1);
}

TEST_CASE("Scarf command sections preserve tapered design heights across ordered child splits", "[DepositionSection]")
{
    DepositionCommandSection parent{{0., 0., .1}, {10., 0., .2}, .4, .1, .2, 2};
    DepositionSectionChain chain(parent);
    const auto first = chain.append({0., 0., .1}, {4., 0., .14}, 2);
    REQUIRE(first);
    CHECK_THAT((*first)[0], Catch::Matchers::WithinAbs(.1, 1e-7));
    CHECK_THAT((*first)[1], Catch::Matchers::WithinAbs(.14, 1e-7));
    CHECK_FALSE(chain.complete());
    const auto second = chain.append({4., 0., .14}, {10., 0., .2}, 2);
    REQUIRE(second);
    CHECK_THAT((*second)[1], Catch::Matchers::WithinAbs(.2, 1e-7));
    CHECK(chain.complete());
    DepositionSectionChain wrong_material(parent), skipped(parent), displaced(parent);
    CHECK_FALSE(wrong_material.append(parent.start, parent.end, 1));
    CHECK_FALSE(skipped.append({1., 0., .11}, parent.end, 2));
    CHECK_FALSE(displaced.append(parent.start, {10., .01, .2}, 2));
}

TEST_CASE("Section transcripts reject lost duplicated reordered and stale command scopes", "[DepositionSection]")
{
    DepositionSections registry(42);
    const auto a = registry.add({{0., 0., .2}, {1., 0., .2}, .4, .2, 0., 0});
    const auto b = registry.add({{1., 0., .2}, {2., 0., .2}, .4, 0., .2, 0});
    registry.seal();
    DepositionSectionCapture valid(42), missing(42), duplicate(42), reversed(42), stale(43);
    for (auto id : {a, 0u, b, 0u}) {
        valid.consume(registry.marker(id));
        stale.consume(registry.marker(id));
    }
    for (auto id : {a, 0u})
        missing.consume(registry.marker(id));
    for (auto id : {a, 0u, a, 0u})
        duplicate.consume(registry.marker(id));
    for (auto id : {b, 0u, a, 0u})
        reversed.consume(registry.marker(id));
    CHECK(valid.finish(registry));
    CHECK_FALSE(missing.finish(registry));
    CHECK_FALSE(duplicate.finish(registry));
    CHECK_FALSE(reversed.finish(registry));
    CHECK_FALSE(stale.finish(registry));
    CHECK_THROWS_AS(registry.add({}), std::length_error);
}

TEST_CASE("Section chains require endpoints in one consistent coordinate representation", "[DepositionSection]")
{
    const DepositionCommandSection parent{{double(float(147.819)), double(float(44.993)), double(float(.2))},
        {double(float(147.985)), double(float(45.424)), double(float(.2))}, .4, .2, .2, 0};
    auto start = parent.start;
    start[0] -= .244141e-6;
    start[1] -= .030518e-6;
    DepositionSectionChain mixed_representation(parent);
    CHECK_FALSE(mixed_representation.append(start, parent.end, 0));
    for (auto& coordinate : start) coordinate = double(float(coordinate));
    DepositionSectionChain equivalent(parent);
    CHECK(equivalent.append(start, parent.end, 0));
    CHECK(equivalent.complete());
    start[0] -= 2e-6;
    DepositionSectionChain different(parent);
    CHECK_FALSE(different.append(start, parent.end, 0));
}
