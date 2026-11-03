#include "SlicedGCode.hpp"
#include "libslic3r/Config.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/GCode/ToolMapping.hpp"
#include <boost/filesystem/operations.hpp>
#include <cstddef>
#include <string>
#include <stdexcept>
#include <vector>
#include "bbs_3mf.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/ProjectTask.hpp"
#include <boost/filesystem.hpp>

namespace Slic3r {
namespace {
bool store_package(const std::string &path, const Print &print, const Print &statistics_print,
                   GCodeProcessorResult &result, const std::string &gcode_path,
                   const std::vector<int> *tools)
{
    Model model(print.model());
    // The 3MF writer stages its configuration through a model backup folder.
    // Give this temporary model its own folder, independent of GUI data_dir.
    model.set_backup_path((boost::filesystem::temp_directory_path() /
        boost::filesystem::unique_path("orca-dispatch-%%%%-%%%%-%%%%")).string());
    DynamicPrintConfig config(print.full_print_config());
    PlateData plate;
    plate.plate_index = print.get_plate_index();
    plate.gcode_file = gcode_path;
    plate.is_sliced_valid = true;
    plate.is_label_object_enabled = result.label_object_enabled;
    plate.is_support_used = statistics_print.is_support_used();
    plate.filament_maps = config.option<ConfigOptionInts>("filament_map")->values;
    plate.parse_filament_info(&result);
    if (tools) {
        // Equivalence admits only static independent nozzles. A logical
        // layered grouping must not override the mapped physical identities.
        plate.nozzle_group_result.reset();
        for (auto &filament : plate.slice_filaments_info) {
            filament.id = tools->at(static_cast<std::size_t>(filament.id));
            filament.group_id.clear();
            filament.nozzle_diameter = 0.;
            filament.nozzle_volume_type.clear();
        }
    }
    for (auto &filament : plate.slice_filaments_info) {
        filament.type = config.opt_string("filament_type", unsigned(filament.id));
        filament.color = config.opt_string("filament_colour", unsigned(filament.id));
        filament.filament_id = config.opt_string("filament_ids", unsigned(filament.id));
    }
    plate.gcode_prediction = std::to_string(int(result.print_statistics.modes[0].time));
    plate.gcode_weight = std::to_string(statistics_print.print_statistics().total_weight);
    plate.printer_model_id = config.opt_string("printer_model");
    plate.nozzle_diameters = config.option("nozzle_diameter")->serialize();
    for (size_t object = 0; object < model.objects.size(); ++object)
        for (size_t instance = 0; instance < model.objects[object]->instances.size(); ++instance)
            if (model.objects[object]->instances[instance]->is_printable())
                plate.objects_and_instances.emplace_back(int(object), int(instance));
    StoreParams params;
    params.path = path;
    params.model = &model;
    params.config = &config;
    params.plate_data_list = {&plate};
    params.export_plate_idx = plate.plate_index;
    params.strategy = SaveStrategy::Silence | SaveStrategy::SkipModel | SaveStrategy::WithGcode | SaveStrategy::SkipAuxiliary;
    return store_bbs_3mf(params);
}
} // namespace
bool store_sliced_gcode_3mf(const std::string &path, const Print &print,
                          GCodeProcessorResult &result, const std::string &gcode_path)
{
    return store_package(path, print, print, result, gcode_path, nullptr);
}

void validate_tool_mapped_sliced_gcode_3mf(const Print &source_print, const Print &mapped_print,
    const GCodeProcessorResult &source_result, const std::vector<int> &tools)
{
    if (!equivalent_tool_mapping_config(source_print.full_print_config(), mapped_print.full_print_config(), tools))
        throw std::invalid_argument("Sliced package settings require normal regeneration.");
    if (!source_result.warnings.empty())
        throw std::invalid_argument("Sliced package warning metadata requires normal regeneration.");
    for (const char *key : {"filament_type", "filament_colour", "filament_ids"}) {
        const auto *values = mapped_print.full_print_config().option<ConfigOptionStrings>(key);
        if (!values || values->empty())
            throw std::invalid_argument("Sliced package material metadata is incomplete.");
    }
    for (const auto &[material, volume] : source_result.print_statistics.total_volumes_per_extruder) {
        (void)volume;
        if (material >= tools.size() || material >= source_result.filament_diameters.size() ||
            material >= source_result.filament_densities.size())
            throw std::invalid_argument("Sliced package material statistics are incomplete.");
    }
}

bool store_tool_mapped_sliced_gcode_3mf(const std::string &path, const Print &source_print,
    const Print &mapped_print, GCodeProcessorResult &source_result,
    const std::vector<int> &tools, const std::string &mapped_gcode_path)
{
    validate_tool_mapped_sliced_gcode_3mf(source_print, mapped_print, source_result, tools);
    return store_package(path, mapped_print, source_print, source_result, mapped_gcode_path, &tools);
}
}
