#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ThermalAirflowTiming.hpp"
#include "libslic3r/ThermalAirflow.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/GCodeReader.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Utils.hpp"
#include "test_utils.hpp"
#include <boost/nowide/fstream.hpp>
#include <boost/nowide/cstdlib.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::ThermalHistory;
using Catch::Matchers::WithinAbs;

namespace {
constexpr auto normal = static_cast<std::size_t>(PrintEstimatedStatistics::ETimeMode::Normal);

void attributed_moves(GCodeProcessorResult& result)
{
    result.moves.clear();
    const auto add = [&](unsigned line, float duration, EMoveType type) {
        GCodeProcessorResult::MoveVertex move;
        move.gcode_id = line;
        move.time[normal] = duration;
        move.type = type;
        result.moves.push_back(move);
    };
    add(5, 2.f, EMoveType::Travel);
    add(5, 3.f, EMoveType::Travel); // Multiple physical children on one line.
    add(8, 4.f, EMoveType::Travel); // Includes attributed dwell/tool-change time.
    add(8, 1.f, EMoveType::Seam); // Legacy duration retained, position ignored.
    add(1, 0.f, EMoveType::Travel); // Display-only order is irrelevant.
    result.moves.back().internal_only = true;
    result.print_statistics.modes[normal].time = 10.;
}
} // namespace

TEST_CASE("Airflow line boundaries retain all modeled durations and same-line state order", "[ThermalAirflowTiming]")
{
    GCodeProcessorResult source;
    attributed_moves(source);
    AirflowState first, second;
    first.pwm[static_cast<unsigned>(FanRole::Part)] = 128.;
    second.dock_automatic = true;
    second.unsupported_firmware_command = true;
    const std::vector<AirflowEvent> events{{1, {}}, {5, first}, {5, second}, {6, first}, {8, second}, {100, {}}};
    const auto trace = time_airflow_events(source, events);
    REQUIRE(trace.events.size() == events.size());
    CHECK(trace.uses_modeled_timing);
    CHECK_THAT(trace.total_s, WithinAbs(10., 1e-12));
    const std::vector<double> expected{0., 5., 5., 5., 10., 10.};
    for (std::size_t i = 0; i < events.size(); ++i) {
        CHECK(trace.events[i].line_id == events[i].line_id);
        CHECK(trace.events[i].state == events[i].state);
        CHECK_THAT(trace.events[i].time_s, WithinAbs(expected[i], 1e-12));
    }
    CHECK_FALSE(trace.events[2].state.pwm[static_cast<unsigned>(FanRole::Dock)]);
    CHECK(trace.events[2].state.dock_automatic.value());
    CHECK(time_airflow_events(source, {}).events.empty());
}

TEST_CASE("Airflow timing maps a matching parsed stream with arc children and a dwell", "[ThermalAirflowTiming]")
{
    const bool was_bbl = GCodeProcessor::s_IsBBLPrinter;
    const ScopeGuard restore([was_bbl] { GCodeProcessor::s_IsBBLPrinter = was_bbl; });
    GCodeProcessor::s_IsBBLPrinter = false;
    const std::vector<std::string> lines{
        "M106 S0", "G90", "M83", ";TYPE:Outer wall", ";HEIGHT:1",
        "G1 X1 Y1 Z1 F600", "G1 X3 Y1 E1.2 F600", "M106 S128", "G4 S2",
        "G1 X4 Y1 F600", "M106 P6 S255", "G2 X6 Y1 I1 J0 E1.2 F600", "M107"};
    AirflowTimeline airflow({{0, FanRole::Part}, {6, FanRole::Dock}}, FirmwareFanPolicy::RmeLegacy);
    GCodeReader reader;
    ScopedTemporaryFile input(".gcode");
    {
        boost::nowide::ofstream stream(input.string());
        for (std::size_t i = 0; i < lines.size(); ++i) {
            stream << lines[i] << '\n';
            reader.parse_line(lines[i], [&](auto&, const auto& line) { airflow.consume(line, i + 1); });
        }
        REQUIRE(stream.good());
    }
    GCodeProcessor processor;
    processor.apply_config(FullPrintConfig{});
    processor.process_file(input.string());
    const auto& source = processor.get_result();
    const auto trace = time_airflow_events(source, airflow.events());
    REQUIRE(trace.events.size() == airflow.events().size());
    CHECK_THAT(trace.events.front().time_s, WithinAbs(0., 1e-12));
    CHECK_THAT(trace.events.back().time_s, WithinAbs(trace.total_s, 1e-12));
    std::size_t arc_children = 0;
    double arc_end = 0.;
    for (const auto& move : source.moves) {
        if (move.internal_only) continue;
        if (move.gcode_id <= 12) arc_end += move.time[normal];
        if (move.gcode_id == 12 && move.type == EMoveType::Extrude) ++arc_children;
    }
    CHECK(arc_children > 1);
    const auto arc_event = time_airflow_events(source, {{12, {}}});
    CHECK_THAT(arc_event.events.front().time_s, WithinAbs(arc_end, 1e-12));
    // A second parse of the identical motion without the dwell demonstrates
    // that attributed waiting time is retained, not just extrusion duration.
    {
        boost::nowide::ofstream stream(input.string());
        for (std::size_t i = 0; i < lines.size(); ++i)
            stream << (i == 8 ? "G4 S0" : lines[i]) << '\n';
        REQUIRE(stream.good());
    }
    GCodeProcessor without_dwell;
    without_dwell.apply_config(FullPrintConfig{});
    without_dwell.process_file(input.string());
    const auto faster = time_airflow_events(without_dwell.get_result(), {});
    CHECK_THAT(trace.total_s - faster.total_s, WithinAbs(2., 1e-4));
}

TEST_CASE("Airflow timing rejects missing durations unordered identities and invalid states", "[ThermalAirflowTiming]")
{
    GCodeProcessorResult source;
    attributed_moves(source);
    source.print_statistics.modes[normal].time = 0.;
    CHECK_THROWS_AS(time_airflow_events(source, {{1, {}}}), std::invalid_argument);
    attributed_moves(source);
    source.moves[1].time[normal] = -1.f;
    CHECK_THROWS_AS(time_airflow_events(source, {}), std::invalid_argument);
    attributed_moves(source);
    source.moves[1].gcode_id = 4;
    CHECK_THROWS_AS(time_airflow_events(source, {}), std::invalid_argument);
    attributed_moves(source);
    CHECK_THROWS_AS(time_airflow_events(source, {{5, {}}, {4, {}}}), std::invalid_argument);
    CHECK_THROWS_AS(time_airflow_events(source, {{0, {}}}), std::invalid_argument);
    AirflowState invalid;
    invalid.pwm[0] = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(time_airflow_events(source, {{1, invalid}}), std::invalid_argument);
    source.moves.back().time[normal] = 1.f;
    source.print_statistics.modes[normal].time = 11.;
    CHECK_THROWS_AS(time_airflow_events(source, {}), std::invalid_argument);
}

TEST_CASE("Airflow timing bounds work and propagates cancellation", "[ThermalAirflowTiming]")
{
    GCodeProcessorResult source;
    attributed_moves(source);
    AirflowTimingLimits limits;
    limits.max_events = 1;
    CHECK_THROWS_AS(time_airflow_events(source, {{1, {}}, {2, {}}}, limits), std::invalid_argument);
    limits = {};
    limits.max_moves = 1;
    CHECK_THROWS_AS(time_airflow_events(source, {}, limits), std::invalid_argument);
    limits = {};
    limits.max_segments = 0;
    CHECK_THROWS_AS(time_airflow_events(source, {}, limits), std::invalid_argument);
    CHECK_THROWS_AS(time_airflow_events(source, {{1, {}}}, {}, [] { throw std::runtime_error("cancel"); }), std::runtime_error);
}

// Hidden: opt-in read-only admission of an external final G-code stream. This
// does not infer firmware-specific channels or validate physical cooling.
TEST_CASE("An external final G-code stream retains ordered modeled airflow timing", "[.ThermalAirflowTimingProject]")
{
    const char* selected = boost::nowide::getenv("ORCA_AIRFLOW_GCODE");
    REQUIRE(selected != nullptr);
    const boost::filesystem::path path(selected);
    REQUIRE(boost::filesystem::is_regular_file(path));
    constexpr std::size_t max_bytes = 256u * 1024u * 1024u;
    constexpr std::size_t max_moves = 2000000;
    const auto bytes = boost::filesystem::file_size(path);
    const auto modified = boost::filesystem::last_write_time(path);
    REQUIRE(bytes > 0);
    REQUIRE(bytes <= max_bytes);
    const auto started = std::chrono::steady_clock::now();
    GCodeProcessor processor;
    const auto checkpoint = [&] {
        if (std::chrono::steady_clock::now() - started > std::chrono::seconds(120))
            throw std::runtime_error("External airflow timing exceeded its time budget.");
        if (processor.get_result().moves.size() > max_moves)
            throw std::runtime_error("External airflow timing exceeded its move budget.");
        if (boost::filesystem::file_size(path) != bytes || boost::filesystem::last_write_time(path) != modified)
            throw std::runtime_error("External G-code changed during airflow timing admission.");
    };
    const bool was_bbl = GCodeProcessor::s_IsBBLPrinter;
    const ScopeGuard restore([was_bbl] { GCodeProcessor::s_IsBBLPrinter = was_bbl; });
    GCodeProcessor::s_IsBBLPrinter = false;
    processor.apply_config(FullPrintConfig{});
    processor.process_file(path.string(), checkpoint);
    checkpoint();
    AirflowTimeline airflow({{0, FanRole::Part}}, FirmwareFanPolicy::GenericM106, 100000);
    GCodeReader reader;
    std::size_t line_id = 0;
    // Both consumers read this same file. Size/time guards detect ordinary
    // concurrent edits, not malicious content replacement or source identity.
    REQUIRE(reader.parse_file(path.string(), [&](auto&, const auto& line) {
        if (line_id % 1024 == 0) checkpoint();
        airflow.consume(line, ++line_id);
    }));
    checkpoint();
    const auto trace = time_airflow_events(processor.get_result(), airflow.events(),
                                          {max_moves, max_moves, 100000}, checkpoint);
    checkpoint();
    REQUIRE(trace.events.size() == airflow.events().size());
    CHECK(trace.uses_modeled_timing);
    double previous = 0.;
    for (std::size_t i = 0; i < trace.events.size(); ++i) {
        if (i % 1024 == 0) checkpoint();
        REQUIRE(trace.events[i].time_s >= previous);
        REQUIRE(trace.events[i].time_s <= trace.total_s);
        REQUIRE(trace.events[i].line_id == airflow.events()[i].line_id);
        REQUIRE(trace.events[i].state == airflow.events()[i].state);
        previous = trace.events[i].time_s;
    }
    std::cout << "Modeled airflow timing: lines=" << line_id
              << " moves=" << processor.get_result().moves.size()
              << " events=" << trace.events.size() << " total_s=" << trace.total_s
              << " last_event_s=" << previous << '\n';
}
