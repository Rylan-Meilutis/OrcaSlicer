#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "libslic3r/CoolingAssistant.hpp"
#include "libslic3r/Config.hpp"
#include "libslic3r/Preset.hpp"
#include "libslic3r/PresetBundle.hpp"
#include "libslic3r/PrintConfig.hpp"
#include <string>

using namespace Slic3r;

TEST_CASE("Bundled material cooling is available on custom RME printers", "[CoolingAssistant]")
{
    const std::string path = std::string(TEST_DATA_DIR) + "/../../resources/profiles";
    PresetBundle library, bundle;
    library.load_vendor_configs_from_json(path, "OrcaFilamentLibrary", PresetBundle::LoadSystem,
                                         ForwardCompatibilitySubstitutionRule::EnableSilent, nullptr, false);
    bundle.load_vendor_configs_from_json(path, "Prusa", PresetBundle::LoadSystem,
                                        ForwardCompatibilitySubstitutionRule::EnableSilent, &library, false);
    const std::string material = GENERATE("TPU", "PLA", "PETG", "ABS", "ASA", "PC");
    INFO(material);
    const auto *base = library.filaments.find_preset("Generic " + material + " @System", false);
    const auto *machine = bundle.printers.find_preset("Prusa CORE One INDX RME 0.4 nozzle", false);
    REQUIRE(base);
    REQUIRE(machine);
    Preset edited = *base;
    edited.name = "Custom " + material;
    edited.is_system = false;
    edited.inherits() = base->name;
    Preset printer = *machine;
    printer.name = "Custom RME";
    printer.is_system = false;
    printer.inherits() = machine->name;
    const auto active = PresetWithVendorProfile(printer, machine->vendor);
    if (material == "TPU") {
        // The printer-specific reference only defines Standard; unused target
        // High Flow and Bowden variants must not disqualify it.
        const auto *tpu = bundle.filaments.find_preset("Generic TPU @Prusa CORE One", false);
        REQUIRE(tpu);
        auto trial = edited.config;
        REQUIRE(apply_material_cooling(trial, tpu->config, printer.config));
    }
    // Include the library fallback present in the app's merged collection.
    auto &fallback = bundle.filaments.load_preset("", base->name, base->config, false);
    fallback.is_system = true;
    fallback.vendor = base->vendor;
    fallback.filament_id = base->filament_id;
    const auto *candidate = find_cooling_reference(bundle.filaments, edited, active);
    REQUIRE(candidate);
    CHECK(is_compatible_with_printer(bundle.filaments.get_preset_with_vendor_profile(*candidate), active));
    CHECK(apply_material_cooling(edited.config, candidate->config, printer.config));
    CHECK(find_cooling_reference(bundle.filaments, edited, active) != nullptr);
}

TEST_CASE("Cooling assistance follows material references without changing temperatures", "[CoolingAssistant]")
{
    const std::string type = GENERATE("PLA", "PETG", "TPU", "ASA", "ABS", "PC", "PA", "PP", "PVA", "PETG-CF");
    const bool auxiliary = GENERATE(false, true);
    auto config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict("filament_type", type);
    auto source = config;
    source.set_deserialize_strict("fan_min_speed", "30");
    source.set_deserialize_strict("fan_max_speed", "60");
    source.set_deserialize_strict("additional_cooling_fan_speed", "25");
    source.set_deserialize_strict("overhang_fan_speed", "80");
    source.set_deserialize_strict("nozzle_temperature", "300");
    auto printer = DynamicPrintConfig::full_print_config();
    printer.set_key_value("auxiliary_fan", new ConfigOptionBool(auxiliary));
    const auto original = config;
    REQUIRE(apply_material_cooling(config, source, printer));
    const bool thermal = type == "ASA" || type == "ABS" || type == "PC";
    CHECK(config.opt_serialize("fan_min_speed") == (thermal ? "0" : "30"));
    CHECK(config.opt_serialize("fan_max_speed") == (thermal ? "0" : "60"));
    CHECK(config.opt_serialize("additional_cooling_fan_speed") == (thermal || !auxiliary ? "0" : "25"));
    CHECK(config.opt_serialize("initial_layer_fan_speed") == "0");
    CHECK(config.opt_serialize("overhang_fan_speed") == "80");
    for (const char *key : {"nozzle_temperature", "hot_plate_temp", "chamber_temperature",
             "activate_air_filtration", "during_print_exhaust_fan_speed", "filament_flow_ratio"})
        CHECK(config.opt_serialize(key) == original.opt_serialize(key));
    const auto applied = config;
    REQUIRE(apply_material_cooling(config, source, printer));
    CHECK(config.diff(applied).empty());
}

TEST_CASE("Cooling assistance maps nozzle variants by name and rejects missing matches atomically", "[CoolingAssistant]")
{
    auto config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict("filament_type", "PLA");
    config.set_deserialize_strict("filament_extruder_variant", "Direct Drive High Flow;Direct Drive Standard");
    config.set_deserialize_strict("fan_min_speed", "1,2");
    auto source = config;
    source.set_deserialize_strict("filament_extruder_variant", "Direct Drive Standard;Direct Drive High Flow");
    source.set_deserialize_strict("fan_min_speed", "40,60");
    REQUIRE(apply_material_cooling(config, source, config));
    CHECK(config.opt_serialize("fan_min_speed") == "60,2");
    REQUIRE(apply_material_cooling(config, source, config, 1));
    CHECK(config.opt_serialize("fan_min_speed") == "60,40");
    const auto original = config;
    CHECK_FALSE(apply_material_cooling(config, source, config, 2));
    CHECK(config.diff(original).empty());
    source.set_deserialize_strict("filament_extruder_variant", "Bowden Standard");
    CHECK_FALSE(apply_material_cooling(config, source, config));
    CHECK(config.diff(original).empty());
    source.set_deserialize_strict("filament_type", "PETG");
    CHECK_FALSE(apply_material_cooling(config, source, config));
    CHECK(config.diff(original).empty());
}

TEST_CASE("Cooling references exclude other printers and materials", "[CoolingAssistant]")
{
    PresetCollection filaments(Preset::TYPE_FILAMENT, Preset::filament_options(),
                              static_cast<const PrintRegionConfig &>(FullPrintConfig::defaults()));
    Preset printer(Preset::TYPE_PRINTER, "Test printer");
    printer.config = DynamicPrintConfig::full_print_config();
    auto config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict("filament_type", "PETG");
    config.set_deserialize_strict("compatible_printers", "Test printer");
    filaments.load_preset("", "Generic PETG", config, false).is_system = true;
    config.set_deserialize_strict("compatible_printers", "Other printer");
    filaments.load_preset("", "Generic PETG other", config, false).is_system = true;
    Preset edited(Preset::TYPE_FILAMENT, "My PETG");
    edited.config = config;
    const PresetWithVendorProfile active(printer, nullptr);
    const auto *reference = find_cooling_reference(filaments, edited, active);
    REQUIRE(reference);
    CHECK(reference->name == "Generic PETG");
    edited.config.set_deserialize_strict("filament_type", "Unknown");
    CHECK(find_cooling_reference(filaments, edited, active) == nullptr);
}

TEST_CASE("Thermal cooling policy leaves unselected variants unchanged", "[CoolingAssistant]")
{
    auto config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict("filament_type", "ASA");
    config.set_deserialize_strict("filament_extruder_variant", "Direct Drive Standard;Direct Drive High Flow");
    config.set_deserialize_strict("fan_min_speed", "30,45");
    config.set_deserialize_strict("fan_max_speed", "60,75");
    config.set_deserialize_strict("additional_cooling_fan_speed", "20,35");
    const auto source = config;
    REQUIRE(apply_material_cooling(config, source, config));
    CHECK(config.opt_serialize("fan_min_speed") == "0,45");
    CHECK(config.opt_serialize("fan_max_speed") == "0,75");
    CHECK(config.opt_serialize("additional_cooling_fan_speed") == "0,35");
}
