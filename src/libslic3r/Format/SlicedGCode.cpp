#include "SlicedGCode.hpp"
#include "bbs_3mf.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/ProjectTask.hpp"
#include <boost/filesystem.hpp>

namespace Slic3r {
bool store_sliced_gcode_3mf(const std::string &path, const Print &print,
                          GCodeProcessorResult &result, const std::string &gcode_path)
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
    plate.is_support_used = print.is_support_used();
    plate.filament_maps = config.option<ConfigOptionInts>("filament_map")->values;
    plate.parse_filament_info(&result);
    for (auto &filament : plate.slice_filaments_info) {
        filament.type = config.opt_string("filament_type", unsigned(filament.id));
        filament.color = config.opt_string("filament_colour", unsigned(filament.id));
        filament.filament_id = config.opt_string("filament_ids", unsigned(filament.id));
    }
    plate.gcode_prediction = std::to_string(int(result.print_statistics.modes[0].time));
    plate.gcode_weight = std::to_string(print.print_statistics().total_weight);
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
}
