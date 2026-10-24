#include <catch2/catch_test_macros.hpp>
#include "libslic3r/GCode/ToolMapping.hpp"
#include "libslic3r/PlaceholderParser.hpp"
#include "libslic3r/Config.hpp"
#include "libslic3r/PrintConfig.hpp"
#include <cstdint>
#include <cstddef>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Slic3r;

TEST_CASE("Identity dispatch reuse requires the entire effective configuration to match", "[ToolMapping]")
{
    DynamicPrintConfig source;
    source.set_key_value("nozzle_diameter", new ConfigOptionFloats{0.4, 0.4});
    source.set_key_value("machine_start_gcode", new ConfigOptionString("T0\nM151 T500\n"));
    auto target = source;
    CHECK(can_reuse_identity_tool_mapping(source, target, {0, 1}));
    CHECK_FALSE(can_reuse_identity_tool_mapping(source, target, {1, 0}));
    CHECK_FALSE(can_reuse_identity_tool_mapping(source, target, {-1, 1}));
    CHECK_FALSE(can_reuse_identity_tool_mapping(source, target, {}));
    target.set_key_value("nozzle_diameter", new ConfigOptionFloats{0.4, 0.6});
    CHECK_FALSE(can_reuse_identity_tool_mapping(source, target, {0, 1}));
    target = source;
    target.set_key_value("machine_start_gcode", new ConfigOptionString("T1\nM151 T500\n"));
    CHECK_FALSE(can_reuse_identity_tool_mapping(source, target, {0, 1}));
    target = source;
    target.set_key_value("filament_flow_ratio", new ConfigOptionFloats{1., 0.98});
    CHECK_FALSE(can_reuse_identity_tool_mapping(source, target, {0, 1}));
    CHECK_FALSE(can_reuse_identity_tool_mapping(target, source, {0, 1}));
}

TEST_CASE("Bound tool substitution leaves motion and literal printer commands untouched", "[ToolMapping]")
{
    const std::string original = "T0\nM151 T500\nM104 T1 S265\n; filament 2\nT0 ; literal custom command\nG1 X10 Y20 E.2\n";
    const std::vector<GCodeToolBinding> bindings{
        {1, 0}, {original.find("T1") + 1, 1}, {original.find("filament 2") + 9, 1, 1}};
    std::istringstream source(original);
    std::ostringstream target;
    std::vector<std::uint64_t> progress;
    remap_gcode_bindings(source, target, bindings, {12, 7},
        [&](std::uint64_t bytes) { progress.push_back(bytes); });
    CHECK(target.str() == "T12\nM151 T500\nM104 T7 S265\n; filament 8\nT0 ; literal custom command\nG1 X10 Y20 E.2\n");
    REQUIRE_FALSE(progress.empty());
    CHECK(progress.back() == original.size());
    for (std::size_t i = 1; i < progress.size(); ++i) CHECK(progress[i] >= progress[i-1]);
}

TEST_CASE("Dispatch equivalence checks used material rows and destination nozzle properties", "[ToolMapping]")
{
    DynamicPrintConfig source;
    source.set_key_value("filament_map", new ConfigOptionInts{1, 2});
    source.set_key_value("filament_type", new ConfigOptionStrings{"ASA", "TPU"});
    source.set_key_value("filament_flow_ratio", new ConfigOptionFloats{1., 0.98});
    source.set_key_value("nozzle_diameter", new ConfigOptionFloats{0.4, 0.4, 0.4});
    source.set_key_value("sparse_infill_filament_id", new ConfigOptionInt(1));
    auto target = source;
    CHECK_FALSE(equivalent_tool_mapping_config(source, target, {1, 0}));
    target.set_key_value("filament_map", new ConfigOptionInts{1, 2, 3});
    target.set_key_value("filament_type", new ConfigOptionStrings{"ASA", "TPU", "ASA"});
    target.set_key_value("filament_flow_ratio", new ConfigOptionFloats{1., 0.98, 1.});
    CHECK(can_reuse_identity_tool_mapping(source, target, {0, 1}));

    target.set_key_value("filament_type", new ConfigOptionStrings{"TPU", "ASA", "ASA"});
    target.set_key_value("filament_flow_ratio", new ConfigOptionFloats{0.98, 1., 1.});
    target.set_key_value("sparse_infill_filament_id", new ConfigOptionInt(3));
    CHECK(equivalent_tool_mapping_config(source, target, {2, 0}));
    CHECK_FALSE(can_reuse_identity_tool_mapping(source, target, {2, 0}));
    source.set_key_value("filament_settings_id", new ConfigOptionStrings{"ASA preset", "TPU preset"});
    target.set_key_value("filament_settings_id", new ConfigOptionStrings{"TPU preset", "ASA preset", "ASA preset"});
    source.set_key_value("filament_ids", new ConfigOptionStrings{"ASA id", "TPU id"});
    target.set_key_value("filament_ids", new ConfigOptionStrings{"TPU id", "ASA id", "ASA id"});
    source.set_key_value("filament_self_index", new ConfigOptionInts{1, 2});
    target.set_key_value("filament_self_index", new ConfigOptionInts{1, 2, 3});
    CHECK(equivalent_tool_mapping_config(source, target, {2, 0}));
    target.set_key_value("filament_self_index", new ConfigOptionInts{1, 1, 2});
    CHECK_FALSE(equivalent_tool_mapping_config(source, target, {2, 0}));
    target.set_key_value("filament_self_index", new ConfigOptionInts{1, 2, 3});
    target.set_key_value("nozzle_diameter", new ConfigOptionFloats{0.4, 0.4, 0.6});
    CHECK_FALSE(equivalent_tool_mapping_config(source, target, {2, 0}));
    target.set_key_value("nozzle_diameter", new ConfigOptionFloats{0.4, 0.4, 0.4});
    target.set_key_value("sparse_infill_filament_id", new ConfigOptionInt(1));
    CHECK_FALSE(equivalent_tool_mapping_config(source, target, {2, 0}));
    target.set_key_value("sparse_infill_filament_id", new ConfigOptionInt(3));
    target.set_key_value("filament_flow_ratio", new ConfigOptionFloats{0.98, 1., 0.9});
    CHECK_FALSE(equivalent_tool_mapping_config(source, target, {2, 0}));
    target.set_key_value("filament_flow_ratio", new ConfigOptionFloats{0.98, 1., 1.});
    // Logical materials may share one nozzle-property row during preview,
    // without being a single-extruder multi-material printer.
    source.set_key_value("filament_map", new ConfigOptionInts{1, 1});
    CHECK(equivalent_tool_mapping_config(source, target, {2, 0}));
    source.set_key_value("nozzle_diameter", new ConfigOptionFloats{0.6, 0.4, 0.4});
    CHECK_FALSE(equivalent_tool_mapping_config(source, target, {2, 0}));
}

TEST_CASE("Custom tool templates replay expressions without rewriting literal command parameters", "[ToolMapping]")
{
    PlaceholderParser parser;
    parser.set("initial_tool", 0);
    parser.set("current_filament_id", 1);
    parser.set("current_hotend", 0);
    parser.set("nozzle_temperature", new ConfigOptionInts{265, 220});
    parser.set("is_extruder_used", new ConfigOptionBools{true, true});
    DynamicConfig outputs;
    outputs.set_key_value("e_retracted", new ConfigOptionFloats{0., 0.});
    PlaceholderParser::ContextData context;
    context.global_config = std::make_unique<DynamicConfig>();
    const std::string text =
        "M151 W255 D0 T500\n"
        "T{initial_tool} S1 L2 D0\n"
        "M104 T{initial_tool} S{nozzle_temperature[initial_tool]}\nM104 T{current_hotend} S220\n"
        "{if is_extruder_used[0]}M976 A 0:0:{nozzle_temperature[0]}\n{endif}"
        "{if is_extruder_used[2]}M976 A 2:2:{nozzle_temperature[2]}\n{endif}";
    // Physical-slot loops already see the printer's full capacity when slicing.
    parser.set("nozzle_temperature", new ConfigOptionInts{265, 220, 0});
    parser.set("is_extruder_used", new ConfigOptionBools{true, true, false});
    outputs.set_key_value("e_retracted", new ConfigOptionFloats{0., 0., 0.});
    auto snapshot = GCodeTemplateSnapshot::capture("machine_start_gcode", text, 0, parser, nullptr, outputs, context);
    const auto rendered = parser.process(text, 0, nullptr, &outputs, &context);
    snapshot.finish(rendered, outputs, context);
    DynamicPrintConfig target;
    target.set_key_value("filament_type", new ConfigOptionStrings{"TPU", "ASA", "ASA"});
    target.set_key_value("nozzle_temperature", new ConfigOptionInts{220, 265, 265});
    const auto mapped = replay_tool_mapped_template(snapshot, target, {2, 0});
    CHECK(mapped.find("M151 W255 D0 T500") != std::string::npos);
    CHECK(mapped.find("T2 S1 L2 D0") != std::string::npos);
    CHECK(mapped.find("M104 T2 S265") != std::string::npos);
    CHECK(mapped.find("M104 T0 S220") != std::string::npos);
    CHECK(mapped.find("M976 A 0:0:220") != std::string::npos);
    CHECK(mapped.find("M976 A 2:2:265") != std::string::npos);
    CHECK_THROWS_AS(replay_tool_mapped_template(snapshot, target, {0, 0}), std::invalid_argument);
    snapshot.rendered += "changed";
    CHECK_THROWS_AS(replay_tool_mapped_template(snapshot, target, {2, 0}), std::runtime_error);
}

TEST_CASE("Dispatch accepts only unresolved or static physical nozzle maps", "[ToolMapping]")
{
    DynamicPrintConfig source;
    source.set_key_value("filament_map", new ConfigOptionInts{1, 2, 3});
    source.set_key_value("filament_nozzle_map", new ConfigOptionInts{0, 0, 0});
    source.set_key_value("nozzle_diameter", new ConfigOptionFloats(std::vector<double>(8, 0.4)));
    source.set_key_value("enable_filament_dynamic_map", new ConfigOptionBool(false));
    auto target = source;
    target.set_key_value("filament_map", new ConfigOptionInts{1, 2, 3, 4, 5, 6, 7, 8});
    target.set_key_value("filament_nozzle_map", new ConfigOptionInts(std::vector<int>(8, 0)));
    const std::vector<int> tools{1, 2, 4};
    CHECK(equivalent_tool_mapping_config(source, target, tools));
    // Source after processing versus target before processing is the GUI path.
    source.set_key_value("filament_nozzle_map", new ConfigOptionInts{0, 1, 2});
    CHECK(equivalent_tool_mapping_config(source, target, tools));
    target.set_key_value("filament_nozzle_map", new ConfigOptionInts{0, 1, 2, 3, 4, 5, 6, 7});
    CHECK(equivalent_tool_mapping_config(source, target, tools));
    source.set_key_value("filament_nozzle_map", new ConfigOptionInts{0, 0, 0});
    CHECK(equivalent_tool_mapping_config(source, target, tools));
    source.set_key_value("filament_nozzle_map", new ConfigOptionInts{0, 2, 1});
    CHECK_FALSE(equivalent_tool_mapping_config(source, target, tools));
    source.set_key_value("filament_nozzle_map", new ConfigOptionInts{0, 1, 2});
    source.set_key_value("enable_filament_dynamic_map", new ConfigOptionBool(true));
    CHECK_FALSE(equivalent_tool_mapping_config(source, target, tools));
    source.set_key_value("enable_filament_dynamic_map", new ConfigOptionBool(false));
    target.set_key_value("enable_filament_dynamic_map", new ConfigOptionBool(true));
    CHECK_FALSE(equivalent_tool_mapping_config(source, target, tools));
}

TEST_CASE("Tool replay relocates first-layer lists and effective retract aliases", "[ToolMapping]")
{
    PlaceholderParser parser;
    parser.set("first_tools", new ConfigOptionInts{1, 0});
    parser.set("first_filaments", new ConfigOptionInts{0, 1});
    parser.set("initial_tool", 0);
    parser.set("retract_length", new ConfigOptionFloats{0.7, 1.2});
    DynamicConfig outputs;
    PlaceholderParser::ContextData context;
    context.global_config = std::make_unique<DynamicConfig>();
    const std::string text = "T{first_tools[0]}\nT{first_filaments[0]}\nG1 E-{retract_length[initial_tool]}\n";
    auto snapshot = GCodeTemplateSnapshot::capture("machine_start_gcode", text, 0, parser, nullptr, outputs, context);
    snapshot.finish(parser.process(text, 0, nullptr, &outputs, &context), outputs, context);
    CHECK(replay_tool_mapped_template(snapshot, DynamicPrintConfig{}, {2, 4}) == "T4\nT2\nG1 E-0.7\n");
}

TEST_CASE("Custom tool replay rejects changes to surrounding extrusion state", "[ToolMapping]")
{
    PlaceholderParser parser;
    parser.set("initial_tool", 0);
    DynamicConfig outputs;
    outputs.set_key_value("e_retracted", new ConfigOptionFloats{0., 0.});
    PlaceholderParser::ContextData context;
    context.global_config = std::make_unique<DynamicConfig>();
    const std::string text = "{e_retracted[initial_tool] = initial_tool + 1}T{initial_tool}\n";
    auto snapshot = GCodeTemplateSnapshot::capture("change_filament_gcode", text, 0, parser, nullptr, outputs, context);
    snapshot.finish(parser.process(text, 0, nullptr, &outputs, &context), outputs, context);
    DynamicPrintConfig target;
    target.set_key_value("filament_type", new ConfigOptionStrings{"ASA", "ASA"});
    CHECK_THROWS_AS(replay_tool_mapped_template(snapshot, target, {1, 0}), std::runtime_error);
}

TEST_CASE("Fast dispatch preserves motion while replaying macros and remapping native tool operands", "[ToolMapping]")
{
    PlaceholderParser parser;
    parser.set("initial_tool", 0);
    DynamicConfig outputs;
    PlaceholderParser::ContextData context;
    context.global_config = std::make_unique<DynamicConfig>();
    const std::string text = "M151 W255 D0 T500\nT{initial_tool} S1 L2 D0\nM106 S200\nM555 X0 Y0 W20 H20\n";
    auto snapshot = GCodeTemplateSnapshot::capture("machine_start_gcode", text, 0, parser, nullptr, outputs, context);
    snapshot.finish(parser.process(text, 0, nullptr, &outputs, &context), outputs, context);
    GCodeDispatchData data;
    data.templates.push_back(snapshot);
    auto processed_macro = snapshot.rendered;
    processed_macro.insert(processed_macro.find("M106 S200"), "M106 S255 ; enable fan\n");
    processed_macro.insert(processed_macro.find("T0 S1"), "M73 P4 R12\n");
    processed_macro += "M73 P5 R12\n";
    processed_macro.replace(processed_macro.find("M555 X0 Y0 W20 H20"), std::string("M555 X0 Y0 W20 H20").size(),
                            "M555 X1.0000 Y2.0000 W22.0000 H23.0000");
    const std::string motion = "G1 X10.125 Y20.250 Z.2 E.424 F1200 ; external perimeter\n";
    const std::string original = "; ORCA_DISPATCH_BEGIN 0\n" + processed_macro + "; ORCA_DISPATCH_END 0\n"
        "T0\nM104 T1 S220\nM204 P7000 R2500 T7000\nM205 S0 T0\n" + motion +
        "; filament: 1,2\n; CONFIG_BLOCK_START\n; filament_type = ASA;TPU\n; CONFIG_BLOCK_END\n"
        "; filament used [mm] = 12.5, 3.25\n";
    DynamicPrintConfig target;
    target.set_key_value("gcode_flavor", new ConfigOptionEnum<GCodeFlavor>(gcfMarlinFirmware));
    target.set_key_value("filament_type", new ConfigOptionStrings{"TPU", "ASA", "ASA"});
    target.set_key_value("filament_map", new ConfigOptionInts{1, 2, 3});
    std::istringstream generated(original);
    seal_gcode_dispatch(generated, data);
    std::istringstream input(original);
    const auto plan = prepare_tool_mapped_gcode(input, data, target, {2, 0},
        "; filament_type = TPU;ASA;ASA\n; filament_nozzle_map = 0,0,0\n");
    CHECK(plan.source_bytes == original.size());
    std::istringstream source(original);
    std::ostringstream mapped;
    apply_tool_mapped_gcode(source, mapped, plan);
    CHECK(mapped.str().find(motion) != std::string::npos);
    CHECK(mapped.str().find("M151 W255 D0 T500\nM73 P4 R12\nT2 S1 L2 D0\n") != std::string::npos);
    CHECK(mapped.str().find("M73 P5 R12\n; ORCA_DISPATCH_END") != std::string::npos);
    CHECK(mapped.str().find("T2 S1 L2 D0\nM106 S255 ; enable fan\nM106 S200\n") != std::string::npos);
    CHECK(mapped.str().find("T2\nM104 T0 S220\n") != std::string::npos);
    CHECK(mapped.str().find("M204 P7000 R2500 T7000\nM205 S0 T0\n") != std::string::npos);
    CHECK(mapped.str().find("M555 X1.0000 Y2.0000 W22.0000 H23.0000\n") != std::string::npos);
    CHECK(mapped.str().find("; filament: 3,1\n") != std::string::npos);
    CHECK(mapped.str().find("; filament used [mm] = 3.25, 0, 12.5\n") != std::string::npos);
    CHECK(mapped.str().find("; filament_type = TPU;ASA;ASA\n") != std::string::npos);
    CHECK(mapped.str().find("; filament_nozzle_map = 0,1,2\n") != std::string::npos);

    auto rewritten_macro = original;
    rewritten_macro.replace(rewritten_macro.find("T0 S1"), std::string("T0 S1").size(), "T9 S1");
    auto rewritten_data = data;
    std::istringstream rewritten_seal(rewritten_macro);
    seal_gcode_dispatch(rewritten_seal, rewritten_data);
    std::istringstream rewritten_input(rewritten_macro);
    CHECK_THROWS_AS(prepare_tool_mapped_gcode(rewritten_input, rewritten_data, target, {2, 0},
        "; filament_type = TPU;ASA;ASA\n"), std::runtime_error);

    auto changed = original;
    changed.replace(changed.find("X10.125"), 7, "X90.125");
    std::istringstream stale(changed);
    std::ostringstream discarded;
    CHECK_THROWS_AS(apply_tool_mapped_gcode(stale, discarded, plan), std::runtime_error);
    std::istringstream edited_before_preflight(changed);
    CHECK_THROWS_AS(prepare_tool_mapped_gcode(edited_before_preflight, data, target, {2, 0},
        "; filament_type = TPU;ASA;ASA\n"), std::runtime_error);
    std::istringstream cancel_source(original);
    std::ostringstream cancel_target;
    CHECK_THROWS_AS(apply_tool_mapped_gcode(cancel_source, cancel_target, plan,
        [](std::uint64_t) { throw std::runtime_error("cancelled"); }), std::runtime_error);
}

TEST_CASE("Template edits retain relocated cooling and split motions without cascading substitutions", "[ToolMapping]")
{
    PlaceholderParser parser;
    parser.set("initial_tool", 0);
    parser.set("next_extruder", 1);
    DynamicConfig outputs;
    PlaceholderParser::ContextData context;
    context.global_config = std::make_unique<DynamicConfig>();
    const std::string text = "T{initial_tool}\nG1 X10 E1 F1200\nM106 S200\nT{next_extruder}\n";
    auto snapshot = GCodeTemplateSnapshot::capture("change_filament_gcode", text, 0, parser, nullptr, outputs, context);
    snapshot.finish(parser.process(text, 0, nullptr, &outputs, &context), outputs, context);
    const std::string split_motion = "G1 X5 E.5 F1200\nM106 S200\nG1 X10 E.5 F1200\n";
    const std::string processed = "T0\n" + split_motion + "T1\n";
    DynamicPrintConfig target;
    target.set_key_value("gcode_flavor", new ConfigOptionEnum<GCodeFlavor>(gcfMarlinFirmware));
    const auto remap = [&](const std::string &body) {
        GCodeDispatchData data;
        data.templates.push_back(snapshot);
        const auto original = "; ORCA_DISPATCH_BEGIN 0\n" + body + "; ORCA_DISPATCH_END 0\n";
        std::istringstream seal(original);
        seal_gcode_dispatch(seal, data);
        std::istringstream preflight(original);
        const auto plan = prepare_tool_mapped_gcode(preflight, data, target, {1, 0}, "");
        std::istringstream input(original);
        std::ostringstream output;
        apply_tool_mapped_gcode(input, output, plan);
        return output.str();
    };
    // Changing T0 to T1 must not subsequently remap that newly-written T1.
    CHECK(remap(processed) == "; ORCA_DISPATCH_BEGIN 0\nT1\n" + split_motion + "T0\n; ORCA_DISPATCH_END 0\n");
    CHECK_THROWS_AS(remap(processed + "T0\n"), std::runtime_error);
    CHECK_THROWS_AS(remap(split_motion + "T1\n"), std::runtime_error);
    auto rewritten = processed;
    rewritten.replace(0, 3, "T9\n");
    CHECK_THROWS_AS(remap(rewritten), std::runtime_error);
}

TEST_CASE("Custom tool replay rejects changed RME cleaning parking and driver control", "[ToolMapping]")
{
    for (const std::string command : {"G12 S", "G27 Z", "G4 P", "M906 P", "M569 S", "M593 F"}) {
        PlaceholderParser parser;
        parser.set("initial_tool", 0);
        DynamicConfig outputs;
        PlaceholderParser::ContextData context;
        context.global_config = std::make_unique<DynamicConfig>();
        const auto text = command + "{initial_tool}\n";
        auto snapshot = GCodeTemplateSnapshot::capture("change_filament_gcode", text, 0, parser, nullptr, outputs, context);
        snapshot.finish(parser.process(text, 0, nullptr, &outputs, &context), outputs, context);
        CHECK_THROWS_AS(replay_tool_mapped_template(snapshot, DynamicPrintConfig{}, {1}), std::runtime_error);
    }
}

TEST_CASE("Tool replay retains effective material retraction overrides", "[ToolMapping]")
{
    PlaceholderParser parser;
    parser.set("initial_tool", 0);
    parser.set("retraction_length", new ConfigOptionFloats{0.7});
    DynamicConfig outputs;
    outputs.set_key_value("e_retracted", new ConfigOptionFloats{0.});
    PlaceholderParser::ContextData context;
    context.global_config = std::make_unique<DynamicConfig>();
    const std::string text = "T{initial_tool}\n{e_retracted[initial_tool] = retraction_length[initial_tool]}";
    auto snapshot = GCodeTemplateSnapshot::capture("change_filament_gcode", text, 0, parser, nullptr, outputs, context);
    snapshot.finish(parser.process(text, 0, nullptr, &outputs, &context), outputs, context);
    DynamicPrintConfig target;
    target.set_key_value("filament_map", new ConfigOptionInts{1, 2, 3});
    target.set_key_value("retraction_length", new ConfigOptionFloats{1., 1., 1.});
    target.set_key_value("filament_retraction_length", new ConfigOptionFloatsNullable{1., 1., 0.7});
    CHECK(replay_tool_mapped_template(snapshot, target, {2}) == "T2\n");
    target.set_key_value("filament_retraction_length", new ConfigOptionFloatsNullable{1., 1., 0.9});
    CHECK_THROWS_AS(replay_tool_mapped_template(snapshot, target, {2}), std::runtime_error);
}

TEST_CASE("Tool replay retains implicit machine modes without changing explicit mode indices", "[ToolMapping]")
{
    PlaceholderParser parser;
    parser.set("machine_max_acceleration_e", new ConfigOptionFloats{5000., 2500.});
    parser.set("initial_tool", 0);
    const std::string text = "T{initial_tool}\nM201 E[machine_max_acceleration_e]\nM201 E[machine_max_acceleration_e_1]\n";
    // Normal generation remains backward-compatible: its implicit legacy
    // index is still the supplied tool. Only certified replay pins that index.
    CHECK(parser.process("[machine_max_acceleration_e]", 1) == "2500");
    DynamicConfig outputs;
    PlaceholderParser::ContextData context;
    context.global_config = std::make_unique<DynamicConfig>();
    auto snapshot = GCodeTemplateSnapshot::capture("change_filament_gcode", text, 0, parser, nullptr, outputs, context);
    snapshot.finish(parser.process(text, 0, nullptr, &outputs, &context), outputs, context);
    DynamicPrintConfig target;
    target.set_key_value("machine_max_acceleration_e", new ConfigOptionFloats{5000., 2500.});
    CHECK(replay_tool_mapped_template(snapshot, target, {1}) == "T1\nM201 E5000\nM201 E2500\n");

    context.legacy_machine_mode_index = 1;
    snapshot = GCodeTemplateSnapshot::capture("change_filament_gcode", text, 0, parser, nullptr, outputs, context);
    snapshot.finish(parser.process(text, 0, nullptr, &outputs, &context), outputs, context);
    CHECK(replay_tool_mapped_template(snapshot, target, {1}) == "T1\nM201 E2500\nM201 E2500\n");
}

TEST_CASE("Fast dispatch refuses unproven commands and transformed macro spans", "[ToolMapping]")
{
    DynamicPrintConfig target;
    target.set_key_value("gcode_flavor", new ConfigOptionEnum<GCodeFlavor>(gcfMarlinFirmware));
    GCodeDispatchData data;
    std::istringstream seal_source("M151 T500\nG1 X0\n");
    seal_gcode_dispatch(seal_source, data);
    std::istringstream unknown("M151 T500\nG1 X0\n");
    CHECK_THROWS_AS(prepare_tool_mapped_gcode(unknown, data, target, {1}, ""), std::runtime_error);
    std::istringstream malformed("T0.5\n");
    CHECK_THROWS_AS(prepare_tool_mapped_gcode(malformed, data, target, {1}, ""), std::runtime_error);
    std::istringstream missing("; ORCA_DISPATCH_BEGIN 0\nT0\n; ORCA_DISPATCH_END 0\n");
    CHECK_THROWS_AS(prepare_tool_mapped_gcode(missing, data, target, {1}, ""), std::runtime_error);
    target.set_key_value("post_process", new ConfigOptionStrings{"script.sh"});
    std::istringstream postprocessed("T0\n");
    CHECK_THROWS_AS(prepare_tool_mapped_gcode(postprocessed, data, target, {1}, ""), std::runtime_error);
}

TEST_CASE("Bound tool substitution rejects stale or ambiguous bindings and cancels copying", "[ToolMapping]")
{
    const auto rewrite = [](const std::string &text, const std::vector<GCodeToolBinding> &bindings,
                            const std::vector<int> &tools) {
        std::istringstream source(text);
        std::ostringstream target;
        remap_gcode_bindings(source, target, bindings, tools);
        return target.str();
    };
    CHECK_THROWS_AS(rewrite("T0\n", {{1,0}}, {}), std::invalid_argument);
    CHECK_THROWS_AS(rewrite("T0\n", {{1,0}}, {1,1}), std::invalid_argument);
    CHECK_THROWS_AS(rewrite("T0\n", {{1,0}}, {-1}), std::invalid_argument);
    CHECK_THROWS_AS(rewrite("T0\n", {{1,1}}, {1}), std::invalid_argument);
    CHECK_THROWS_AS(rewrite("T0\n", {{1,0}, {1,0}}, {1}), std::invalid_argument);
    CHECK_THROWS_AS(rewrite("T1\n", {{1,0}}, {1}), std::runtime_error);
    CHECK_THROWS_AS(rewrite("T0\n", {{100,0}}, {1}), std::runtime_error);
    CHECK_THROWS_AS(rewrite("T01\n", {{1,0}}, {1}), std::runtime_error);
    CHECK_THROWS_AS(rewrite("T10\n", {{2,0}}, {1}), std::runtime_error);
    CHECK_THROWS_AS(rewrite("T-0\n", {{2,0}}, {1}), std::runtime_error);
    CHECK_THROWS_AS(rewrite("T+0\n", {{2,0}}, {1}), std::runtime_error);
    CHECK_THROWS_AS(rewrite("T0.5\n", {{1,0}}, {1}), std::runtime_error);
    CHECK_THROWS_AS(rewrite("T.0\n", {{2,0}}, {1}), std::runtime_error);
    std::istringstream source(std::string(200000, 'x'));
    std::ostringstream target;
    CHECK_THROWS_AS(remap_gcode_bindings(source, target, {}, {0}, [](std::uint64_t bytes) {
        if (bytes >= 65536) throw std::runtime_error("cancelled");
    }), std::runtime_error);
    CHECK(target.str().size() == 65536);
}
