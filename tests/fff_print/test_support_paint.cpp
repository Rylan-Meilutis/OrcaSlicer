#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "test_helpers.hpp"
#include "test_utils.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Config.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/ExtrusionEntityCollection.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/PrintBase.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Preset.hpp"
#include "libslic3r/PresetBundle.hpp"
#include "libslic3r/Semver.hpp"
#include "libslic3r/Utils.hpp"
#include "libslic3r/Format/bbs_3mf.hpp"
#include "libslic3r/Format/ProjectStrategy.hpp"
#include "libslic3r/Support/TreeSupport.hpp"
#include "libslic3r/Support/TypedSupportPaint.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/TriangleSelector.hpp"
#include "libslic3r/libslic3r.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

using namespace Slic3r;

namespace {
TriangleMesh supported_roof(float width, float depth)
{
    TriangleMesh mesh = make_cube(2, depth, 8.5);
    mesh.translate(width / 2 - 1, 0, 0);
    TriangleMesh roof = make_cube(width, depth, 2);
    roof.translate(0, 0, 8);
    mesh.merge(roof);
    return mesh;
}

void paint_roof(ModelVolume &volume, EnforcerBlockerType first = EnforcerBlockerType::NORMAL_SUPPORT,
                EnforcerBlockerType second = EnforcerBlockerType::TREE_SUPPORT)
{
    TriangleSelector selector(volume.mesh());
    const indexed_triangle_set &mesh = volume.mesh().its;
    const float contact_z = float(volume.mesh().bounding_box().max.z()) - 2.f;
    int painted = 0;
    for (int facet = 0; facet < int(mesh.indices.size()); ++facet) {
        const auto &triangle = mesh.indices[facet];
        if (its_face_normal(mesh, facet).z() < -0.99f &&
            std::abs(mesh.vertices[triangle[0]].z() - contact_z) < 0.001f) {
            selector.set_facet(facet, painted++ == 0 ? first : second);
        }
    }
    REQUIRE(painted == 2);
    volume.supported_facets.set(selector);
}

DynamicPrintConfig support_config(const char *type)
{
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    config.set_deserialize_strict({
        {"enable_support", 1}, {"support_type", type}, {"support_style", "tree_hybrid"},
        {"support_on_build_plate_only", 0}, {"support_threshold_angle", 30},
        {"raft_layers", 0}, {"layer_height", 0.2}, {"initial_layer_print_height", 0.2},
        {"max_bridge_length", 0}, {"brim_type", "no_brim"}, {"skirt_loops", 0}
    });
    return config;
}
}

TEST_CASE("One painted object prints normal columns and tree branches together", "[SupportPaint]")
{
    const char *type = GENERATE("tree(auto)", "tree(manual)");
    const char *style = GENERATE("tree_hybrid", "tree_slim", "tree_strong");
    // Each triangular half is below/above the automatic Hybrid 100 mm2 cutoff.
    const bool large = GENERATE(false, true);
    Print print;
    Model model;
    auto config = support_config(type);
    config.set_deserialize_strict("support_style", style);
    Test::init_print(std::vector<TriangleMesh>{supported_roof(large ? 32.f : 12.f, large ? 12.f : 4.f)},
                     print, model, config);
    paint_roof(*model.objects.front()->volumes.front());
    REQUIRE(has_typed_support_paint(*model.objects.front()));
    print.apply(model, config);
    // Keep only this test's reference to the generator's transient node data.
    // Production deliberately releases it when the support step completes.
    std::shared_ptr<TreeSupportData> data;
    print.set_status_callback([&](const PrintBase::SlicingStatus &status) {
        if (status.text == "Generating tree support toolpaths")
            data = const_cast<PrintObject *>(print.objects().front())->alloc_tree_support_preview_cache();
    });
    print.process();
    REQUIRE(data);
    const std::string gcode = Test::gcode(print);
    REQUIRE(print.objects().size() == 1);
    const PrintObject *object = print.objects().front();
    REQUIRE_FALSE(object->support_layers().empty());
    Points support_points;
    for (const SupportLayer *layer : object->support_layers())
        layer->support_fills.collect_points(support_points);
    CHECK_FALSE(support_points.empty());
    CHECK(gcode.find(";TYPE:Support") != std::string::npos);

    std::size_t normal_contacts = 0, tree_contacts = 0;
    std::size_t normal_printed_nodes = 0, tree_printed_nodes = 0;
    bool normal_contact_has_extrusion = false, tree_contact_has_extrusion = false;
    for (const auto &node : data->contact_nodes) {
        if (node->parent == SupportNode::NO_PARENT) {
            if (node->type == ePolygon) ++normal_contacts;
            if (node->type == eCircle) ++tree_contacts;
            const bool has_extrusion = std::any_of(support_points.begin(), support_points.end(),
                [&](const Point &point) { return node->overhang.contains(point); });
            if (node->type == ePolygon) normal_contact_has_extrusion |= has_extrusion;
            if (node->type == eCircle) tree_contact_has_extrusion |= has_extrusion;
        }
        if (node->valid && node->distance_to_top >= 0 && node->parent != SupportNode::NO_PARENT) {
            if (node->type == ePolygon) ++normal_printed_nodes;
            if (node->type == eCircle) ++tree_printed_nodes;
        }
    }
    CHECK(normal_contacts > 0);
    CHECK(tree_contacts > 0);
    CHECK(normal_printed_nodes > 0);
    CHECK(tree_printed_nodes > 0);
    CHECK(normal_contact_has_extrusion);
    CHECK(tree_contact_has_extrusion);
}

TEST_CASE("Normal paint does not silently replace an Organic default", "[SupportPaint]")
{
    const char *style = GENERATE("organic", "default");
    Print print;
    Model model;
    auto config = support_config("tree(auto)");
    config.set_deserialize_strict("support_style", style);
    Test::init_print(std::vector<TriangleMesh>{supported_roof(12, 4)}, print, model, config);
    paint_roof(*model.objects.front()->volumes.front());
    print.apply(model, config);
    CHECK(print.validate().opt_key == "support_style");
    config.set_deserialize_strict("enable_support", "0");
    print.apply(model, config);
    CHECK(print.validate().opt_key != "support_style");
}

TEST_CASE("Tree paint changes only marked contacts under a normal default", "[SupportPaint]")
{
    const bool automatic = GENERATE(false, true);
    Print print;
    Model model;
    auto config = support_config(automatic ? "normal(auto)" : "normal(manual)");
    config.set_deserialize_strict("support_style", "snug");
    Test::init_print(std::vector<TriangleMesh>{supported_roof(32, 12)}, print, model, config);
    paint_roof(*model.objects.front()->volumes.front(), EnforcerBlockerType::TREE_SUPPORT,
               EnforcerBlockerType::NONE);
    print.apply(model, config);
    CHECK_FALSE(is_tree(print.objects().front()->config().support_type));
    CHECK(print.objects().front()->config().support_style.value == smsSnug);
    CHECK(print.validate().opt_key != "support_style");
    std::shared_ptr<TreeSupportData> data;
    print.set_status_callback([&](const PrintBase::SlicingStatus &status) {
        if (status.text == "Generating tree support toolpaths")
            data = const_cast<PrintObject *>(print.objects().front())->alloc_tree_support_preview_cache();
    });
    print.process();
    REQUIRE(data);
    std::size_t normal_contacts = 0, tree_contacts = 0;
    std::vector<Polygons> painted_layers;
    print.objects().front()->project_and_append_custom_facets(false, EnforcerBlockerType::TREE_SUPPORT, painted_layers);
    for (const auto &node : data->contact_nodes)
        if (node->parent == SupportNode::NO_PARENT) {
            const size_t source_layer = size_t(node->obj_layer_nr + 1);
            REQUIRE(source_layer < painted_layers.size());
            const ExPolygons tree_interior = offset_ex(painted_layers[source_layer], -scale_(0.1));
            const Layer *layer = print.objects().front()->get_layer(source_layer);
            REQUIRE(layer->lower_layer);
            const double supported_reach = layer->lower_layer->height /
                std::tan((config.opt_int("support_threshold_angle") + 1.) * PI / 180.);
            const ExPolygons automatic_demand = diff_ex(layer->lslices_extrudable,
                offset_ex(layer->lower_layer->lslices_extrudable, scale_(supported_reach)));
            const ExPolygons unpainted_interior = diff_ex(automatic_demand,
                offset_ex(painted_layers[source_layer], scale_(0.1)));
            if (node->type == ePolygon) {
                ++normal_contacts;
                CHECK(intersection_ex(ExPolygons{node->overhang}, tree_interior).empty());
            }
            if (node->type == eCircle) {
                ++tree_contacts;
                const auto beyond_margin = diff_ex(ExPolygons{node->overhang},
                    offset_ex(painted_layers[source_layer], scale_(config.opt_float("tree_support_tip_diameter") + 0.01)));
                CHECK(beyond_margin.empty());
                if (automatic)
                    CHECK(intersection_ex(ExPolygons{node->overhang}, unpainted_interior).empty());
            }
        }
    CHECK(tree_contacts > 0);
    if (automatic) CHECK(normal_contacts > 0);
    else CHECK(normal_contacts == 0);
    CHECK_FALSE(is_tree(print.objects().front()->config().support_type));
    CHECK(print.objects().front()->config().support_style.value == smsSnug);
    CHECK(Test::gcode(print).find(";TYPE:Support") != std::string::npos);
}

TEST_CASE("Tree paint preserves the Organic default and manual detection scope", "[SupportPaint]")
{
    Print print;
    Model model;
    auto config = support_config("tree(manual)");
    config.set_deserialize_strict("support_style", "organic");
    Test::init_print(std::vector<TriangleMesh>{supported_roof(12, 4)}, print, model, config);
    paint_roof(*model.objects.front()->volumes.front(), EnforcerBlockerType::TREE_SUPPORT,
               EnforcerBlockerType::NONE);
    print.apply(model, config);
    CHECK(print.validate().opt_key != "support_style");
    CHECK(Test::gcode(print).find(";TYPE:Support") != std::string::npos);
    CHECK(print.objects().front()->config().support_type.value == stTree);
    CHECK(print.objects().front()->config().support_style.value == smsTreeOrganic);
}

TEST_CASE("Exclusively vertical typed support paint reports its contact limitation", "[SupportPaint]")
{
    Print print;
    Model model;
    const auto config = support_config("tree(manual)");
    Test::init_print(std::vector<TriangleMesh>{Test::cube(4)}, print, model, config);
    ModelVolume &volume = *model.objects.front()->volumes.front();
    TriangleSelector selector(volume.mesh());
    for (int facet = 0; facet < int(volume.mesh().its.indices.size()); ++facet)
        if (its_face_normal(volume.mesh().its, facet).x() > 0.99f)
            selector.set_facet(facet, EnforcerBlockerType::NORMAL_SUPPORT);
    volume.supported_facets.set(selector);
    print.apply(model, config);
    print.process();
    const auto warnings = print.objects().front()->step_state_with_warnings(posSupportMaterial).warnings;
    CHECK(std::count_if(warnings.begin(), warnings.end(), [](const auto &warning) {
        return warning.message.find("no downward-facing contact area") != std::string::npos;
    }) == 1);
}

TEST_CASE("An imported support project uses project defaults after explicit reset", "[SupportPaint][.SupportPaintProject]")
{
    const char *path = std::getenv("ORCA_SUPPORT_PAINT_PROJECT");
    if (!path) SKIP("Set ORCA_SUPPORT_PAINT_PROJECT to a local project with typed support paint");
    const auto previous_resources = resources_dir();
    set_resources_dir(std::string(TEST_DATA_DIR) + "/../../resources");
    const ScopeGuard restore_resources([&] { set_resources_dir(previous_resources); });
    ScopedSlic3rTemporaryDir temporary("support_paint_project");
    DynamicPrintConfig imported;
    ConfigSubstitutionContext substitutions{ForwardCompatibilitySubstitutionRule::Enable};
    PlateDataPtrs plates;
    std::vector<Preset *> presets;
    const ScopeGuard cleanup([&] {
        release_PlateData_list(plates);
        for (Preset *preset : presets) delete preset;
    });
    bool is_bbl = false;
    Semver version;
    auto model = Model::read_from_file(path, &imported, &substitutions,
        LoadStrategy::LoadModel | LoadStrategy::LoadConfig | LoadStrategy::AddDefaultInstances,
        &plates, &presets, &is_bbl, &version);
    Preset::normalize(imported);
    PresetBundle bundle;
    bundle.load_config_model(path, imported, version);
    auto config = bundle.full_config(false);
    if (!plates.empty()) config.apply(plates.front()->config);
    // This audit deliberately isolates the reported fuselage and leaves the
    // user's file untouched. Only the three explicit support overrides are reset.
    for (size_t i = model.objects.size(); i-- > 0;)
        if (model.objects[i]->name.find("Fuselage") == std::string::npos)
            model.delete_object(i);
    REQUIRE(model.objects.size() == 1);
    Print print;
    print.set_status_silent();
    print.apply(model, config);
    WARN("Saved effective support: " << print.objects().front()->config().support_type.serialize()
         << ", " << print.objects().front()->config().support_style.serialize());
    for (const char *key : {"enable_support", "support_type", "support_style"})
        model.objects.front()->config.erase(key);
    print.apply(model, config);
    PrintObject &object = *const_cast<PrintObject *>(print.objects().front());
    REQUIRE(object.config().support_type.value == config.opt_enum<SupportType>("support_type"));
    REQUIRE(object.config().support_style.value == config.opt_enum<SupportMaterialStyle>("support_style"));
    REQUIRE(has_typed_support_paint(*object.model_object()));
    std::shared_ptr<TreeSupportData> data;
    print.set_status_callback([&](const PrintBase::SlicingStatus &status) {
        if (status.text == "Generating tree support toolpaths")
            data = object.alloc_tree_support_preview_cache();
    });
    object.generate_support_preview();
    REQUIRE(data);
    size_t normal = 0, tree = 0;
    for (const auto &node : data->contact_nodes)
        if (node->parent == SupportNode::NO_PARENT) {
            if (node->type == ePolygon) ++normal;
            if (node->type == eCircle) ++tree;
        }
    WARN("Reset effective support: " << object.config().support_type.serialize() << ", "
         << object.config().support_style.serialize() << "; normal contacts=" << normal
         << "; tree contacts=" << tree << "; support layers=" << object.support_layers().size());
    CHECK(normal > 0);
    CHECK(tree > 0);
    CHECK_FALSE(object.support_layers().empty());
}
