#include "SequentialGantryGeometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/nowide/fstream.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/regex.hpp>

#include "Config.hpp"
#include "libslic3r.h"
#include "Utils.hpp"

namespace fs = boost::filesystem;
namespace pt = boost::property_tree;

namespace Slic3r {

namespace {

Polygon parse_polygon(std::string serialized)
{
    std::replace(serialized.begin(), serialized.end(), ';', ' ');
    std::replace(serialized.begin(), serialized.end(), ',', ' ');
    std::stringstream stream(serialized);
    Polygon polygon;
    for (double x = 0., y = 0.; stream >> x >> y;)
        polygon.points.emplace_back(Point::new_scale(x, y));
    return polygon;
}

SequentialGantryGeometry parse_geometry(const pt::ptree &entry, const fs::path &model_directory)
{
    SequentialGantryGeometry geometry;
    if (const auto filename = entry.get_optional<std::string>("gantry_model_filename"))
        geometry.model_path = (model_directory / *filename).string();

    const auto slices = entry.get_child_optional("slices");
    if (!slices)
        return geometry;

    for (const auto &item : *slices) {
        SequentialGantrySlice slice;
        slice.height = item.second.get<double>("height");
        slice.is_box = item.second.get<std::string>("type", "convex") == "box";
        if (const auto polygons = item.second.get_child_optional("polygons")) {
            for (const auto &polygon : *polygons) {
                Polygon parsed = parse_polygon(polygon.second.data());
                if (parsed.points.size() >= 3)
                    slice.polygons.emplace_back(std::move(parsed));
            }
        }
        if (!slice.polygons.empty())
            geometry.slices.emplace_back(std::move(slice));
    }
    return geometry;
}

SequentialGantryGeometry load_custom_geometry(const std::string &custom)
{
    if (custom.empty())
        return {};

    pt::ptree root;
    fs::path model_directory;
    const size_t first = custom.find_first_not_of(" \t\r\n");
    if (first != std::string::npos && custom[first] == '{') {
        std::stringstream stream(custom);
        pt::read_json(stream, root);
    } else {
        fs::path path(custom);
        if (!path.is_absolute() && !fs::exists(path))
            path = fs::path(resources_dir()) / path;
        boost::nowide::ifstream input(path.string());
        pt::read_json(input, root);
        model_directory = path.parent_path();
    }
    return parse_geometry(root, model_directory);
}

} // namespace

double SequentialGantryGeometry::conservative_clearance_radius() const
{
    double radius = 0.;
    for (const SequentialGantrySlice &slice : slices) {
        if (slice.is_box)
            continue;
        for (const Polygon &polygon : slice.polygons)
            for (const Point &point : polygon.points)
                radius = std::max(radius, unscale_(point.cast<double>().norm()));
    }
    return radius;
}

double SequentialGantryGeometry::first_box_height() const
{
    double height = std::numeric_limits<double>::max();
    for (const SequentialGantrySlice &slice : slices)
        if (slice.is_box)
            height = std::min(height, slice.height);
    return height == std::numeric_limits<double>::max() ? 0. : height;
}

SequentialGantryGeometry load_sequential_gantry_geometry(const ConfigBase &config)
{
    try {
        if (config.has("sequential_print_gantry_geometry")) {
            const auto *custom_option = config.option<ConfigOptionString>("sequential_print_gantry_geometry");
            const std::string custom = custom_option == nullptr ? std::string() : custom_option->value;
            if (!custom.empty()) {
                SequentialGantryGeometry geometry = load_custom_geometry(custom);
                if (config.has("sequential_print_gantry_model")) {
                    const auto *model_option = config.option<ConfigOptionString>("sequential_print_gantry_model");
                    const std::string model = model_option == nullptr ? std::string() : model_option->value;
                    if (!model.empty())
                        geometry.model_path = model;
                }
                return geometry;
            }
        }

        const fs::path directory = fs::path(resources_dir()) / "data" / "printer_gantries";
        boost::nowide::ifstream input((directory / "geometries.json").string());
        if (!input)
            return {};
        pt::ptree root;
        pt::read_json(input, root);
        const auto *notes_option = config.option<ConfigOptionString>("printer_notes");
        const std::string notes = notes_option == nullptr ? std::string() : notes_option->value;
        for (const auto &printer : root.get_child("printers")) {
            const boost::regex matcher(printer.second.get<std::string>("printer_notes_regex"));
            if (boost::regex_match(notes, matcher))
                return parse_geometry(printer.second, directory);
        }
    } catch (...) {
        // Invalid optional geometry must not make an otherwise valid machine
        // profile unusable; callers retain the legacy clearance fallback.
    }
    return {};
}

} // namespace Slic3r
