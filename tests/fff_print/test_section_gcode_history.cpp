#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/SectionGCodeHistory.hpp"
#include "libslic3r/ThermoelasticConflictOwnership.hpp"
#include "libslic3r/SectionThermalAssembly.hpp"
#include "libslic3r/SectionIntegration.hpp"
#include "libslic3r/ThermalHistory.hpp"
#include "libslic3r/ThermoelasticToolpath.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/GCode/DepositionProvenance.hpp"
#include "libslic3r/GCode/DepositionSection.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Utils.hpp"
#include "test_utils.hpp"
#include <Eigen/Core>
#include <boost/nowide/fstream.hpp>
#include <stdexcept>
#include <cstddef>
#include <optional>
#include <array>
#include <memory>
#include "libslic3r/ThermoelasticLayers.hpp"
#include "libslic3r/ThermoelasticSourceOwnership.hpp"
#include "libslic3r/ThermoelasticExactSections.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Config.hpp"
#include "test_helpers.hpp"
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;

TEST_CASE("Deposited ownership requires a matching explicit scenario without competing ownership", "[SectionGCodeHistory][DepositedOwnership]")
{
    const int invalid = GENERATE(0, 1, 2, 3, 4);
    GCodeProcessorResult source;
    SectionGCodeScenario scenario(BeadSection::RectangularEnvelope);
    SectionGCodeLimits limits;
    scenario.stream_tile_size_m = Eigen::Vector2d{.002, .002};
    scenario.use_deposited_ownership = invalid != 1;
    limits.ownership.deposited_overlap_policy = invalid == 0 ? DepositedOverlapPolicy::None :
        invalid == 4 ? static_cast<DepositedOverlapPolicy>(255) : DepositedOverlapPolicy::FirstDepositionWins;
    scenario.use_native_ownership = invalid == 2;
    if (invalid == 3)
        scenario.stream_ownership = [](double, double) { return std::optional<SectionOwnership>{}; };
    // Exact error matching prevents absent synthetic timing from making this
    // pass merely because a later, unrelated precondition also rejects.
    CHECK_THROWS_WITH(simulate_section_gcode_history(source, scenario, limits),
        "Deposited overlap ownership requires its separate explicit scenario.");
}

TEST_CASE("Streamed first deposition ownership preserves occupied prefixes and reports its assumption", "[SectionGCodeHistory][DepositedOwnership]")
{
    const bool fixed_reference = GENERATE(false, true);
    const auto config = Test::multifilament_config(2, {{"layer_height", "0.2"},
        {"initial_layer_print_height", "0.2"}, {"brim_type", "no_brim"}, {"skirt_loops", "0"}});
    Print print;
    Model model;
    Test::init_print({Test::cube(2.)}, print, model, config);
    auto& object = *model.objects.front();
    object.volumes.front()->config.set("extruder", 1);
    auto second = Test::cube(2.);
    second.translate(2., 0., 0.);
    object.add_volume(std::move(second))->config.set("extruder", 2);
    print.apply(model, config);
    print.process();
    const auto* printed = print.get_object(0);
    const auto* layer = printed->layers().front();
    auto provenance = std::make_shared<DepositionProvenance>();
    std::array<std::uint32_t, 2> owners{};
    for (const auto* region : layer->regions()) {
        const unsigned filament = region->region().extruder(frPerimeter) - 1;
        REQUIRE(filament < owners.size());
        owners[filament] = provenance->add({printed->get_id(), 0, std::uint64_t(layer->id()),
            std::uint64_t(region->region().print_region_id()), filament, DepositionSource::Height::NativePlanar});
    }
    REQUIRE(owners[0] != 0);
    REQUIRE(owners[1] != 0);
    provenance->seal();
    GCodeProcessorResult source;
    source.deposition_provenance = provenance;
    source.native_ownership = std::make_shared<const NativeOwnershipSnapshot>(capture_native_ownership(print, provenance, Point(0, 0), 0.));
    auto commands = std::make_shared<DepositionSections>(provenance->generation());
    commands->seal();
    source.deposition_sections = commands;
    // Synthetic coincident physical commands with real generation-local
    // source records. This isolates overlap policy, not native XY ownership.
    for (unsigned filament = 0; filament < 2; ++filament) {
        GCodeProcessorResult::MoveVertex move;
        move.type = EMoveType::Travel;
        move.position = {0.f, 0.f, .2f};
        move.time[0] = 1.;
        move.extruder_id = filament;
        source.moves.push_back(move);
        move.type = EMoveType::Extrude;
        move.position.x() = 1.f;
        move.width = .4f;
        move.height = .2f;
        move.mm3_per_mm = .08f;
        move.delta_extruder = .1f;
        move.extrusion_role = erPerimeter;
        move.deposition_source_id = owners[filament];
        source.moves.push_back(move);
    }
    source.print_statistics.modes[0].time = 4.;
    const auto path = timed_deposition(source);
    REQUIRE(path.state == ToolpathState::Complete);
    REQUIRE(path.deposition.size() == 2);
    const auto first_prefix = layered_deposition({path.deposition.front().geometry}, 0.);
    SectionGCodeScenario scenario(BeadSection::RectangularEnvelope);
    scenario.acknowledge_modeled_timing = true;
    scenario.use_deposited_ownership = true;
    scenario.use_fixed_interface_reference = fixed_reference;
    scenario.thermal_materials = {{Eigen::Matrix3d::Identity() * .2, 2e6, 400.},
                                  {Eigen::Matrix3d::Identity() * .2, 2e6, 400.}};
    scenario.expansion = {{{{250., 450., Eigen::Matrix3d::Identity() * 1e-5}}, Eigen::Matrix3d::Identity(), 400.},
                          {{{250., 450., Eigen::Matrix3d::Identity() * 1e-5}}, Eigen::Matrix3d::Identity(), 400.}};
    scenario.mechanical_materials = {{Eigen::Matrix<double, 6, 6>::Identity() * 1000.},
                                     {Eigen::Matrix<double, 6, 6>::Identity() * 1000.}};
    // Physical first filament maps to the HIGHER material ID. Sorting by
    // material rather than move chronology would therefore fail this test.
    scenario.filament_to_material = {1, 0};
    scenario.origin_m = {-.001, -.001, 0.};
    scenario.spacing_m = {.001, .001, .0002};
    scenario.stream_tile_size_m = Eigen::Vector2d{.002, .002};
    scenario.activation_times_s = {path.deposition.front().end_s, path.total_s};
    scenario.deposition_temperatures_k = [](double) { return std::vector<double>{400., 400.}; };
    std::vector<double> prefix_volumes;
    scenario.exposure = [&](double, const SectionThermalAssembly& assembly) {
        prefix_volumes.push_back(assembly.occupied_volume_m3);
        for (const auto& group : assembly.groups) CHECK(group.material == 1);
        return std::vector<ThermalHistory::Element>{};
    };
    SectionGCodeLimits limits;
    limits.ownership.deposited_overlap_policy = DepositedOverlapPolicy::FirstDepositionWins;
    limits.ownership.geometry.retain_exact_interface_sections = true;
    limits.streaming_geometry.geometry.retain_exact_interface_sections = true;
    limits.ownership.geometry.ownership_boundary_grid_units = 4;
    limits.streaming_geometry.geometry.ownership_boundary_grid_units = 4;
    limits.ownership.max_polygon_work = 4000000;
    limits.streaming_geometry.geometry.max_slab_group_visits = 4000000;
    const auto result = simulate_section_gcode_history(source, scenario, limits);
    REQUIRE(result.thermal.state == ThermalHistory::State::Converged);
    REQUIRE(result.mechanical);
    CHECK(result.deposited_overlap_policy == DepositedOverlapPolicy::FirstDepositionWins);
    CHECK(result.native_interface_variant == NativeInterfaceVariant::Strict);
    CHECK_FALSE(result.allows_single_sided_native_boundary);
    CHECK_FALSE(result.allows_dual_native_boundary);
    CHECK(result.uses_fixed_interface_reference == fixed_reference);
    if (fixed_reference) CHECK(result.excluded_nominal_volume_m3 > 0.);
    else CHECK_THAT(result.excluded_nominal_volume_m3, WithinAbs(0., 1e-24));
    CHECK(result.activated_snapshots == 2);
    REQUIRE(prefix_volumes.size() == 2);
    CHECK(prefix_volumes.front() > 0.);
    CHECK_THAT(prefix_volumes.back(), WithinAbs(prefix_volumes.front(), 1e-20));
    CHECK_THAT(result.occupied_volume_m3, WithinAbs(prefix_volumes.front(), 1e-20));
    // A final first-owner reference must retain the complete earlier owner,
    // not carve that prefix according to the later material's lower ID.
    CHECK_THAT(result.occupied_volume_m3, WithinAbs(first_prefix.occupied_volume_m3, 1e-20));
    CHECK_THAT(result.thermal.final_energy_j, WithinAbs(result.occupied_volume_m3 * 2e6 * 400., 1e-8));
    CHECK_THAT(result.thermal.energy_balance_error_j, WithinAbs(0., 1e-8));
    CHECK(source.moves.size() == 4);
    CHECK(source.deposition_provenance == provenance);
    CHECK(source.deposition_sections == commands);
}

TEST_CASE("Section histories reject collapsed envelopes unless their bounded omission is explicit", "[SectionGCodeHistory]")
{
    GCodeProcessorResult source;
    auto provenance = std::make_shared<DepositionProvenance>();
    const auto owner = provenance->add({1, 0, 0, 0, 0, DepositionSource::Height::Transformed,
                                      DepositionSource::Assignment::NominalRole});
    provenance->seal();
    auto sections = std::make_shared<DepositionSections>(provenance->generation());
    DepositionCommandSection bead{{0., 0., .2}, {1., 0., .2}, .4, .2, .2, 0};
    bead.nominal_lower_z = std::array<double, 2>{0., 0.};
    const auto bead_id = sections->add(bead);
    DepositionCommandSection collapsed{{1., 0., .2}, {double(1.05f), 0., .2}, .4, 0., 0., 0};
    collapsed.nominal_lower_z = std::array<double, 2>{.2, .2};
    const auto collapsed_id = sections->add(collapsed);
    sections->seal();
    source.deposition_provenance = provenance;
    source.deposition_sections = sections;
    source.deposition_z_generation = provenance->generation();
    source.filament_diameters = {1.75f};
    source.moves.resize(3);
    for (std::size_t i = 0; i < source.moves.size(); ++i) {
        auto& move = source.moves[i];
        move.position = {i == 0 ? 0.f : i == 1 ? 1.f : 1.05f, 0.f, .2f};
        move.type = i == 0 ? EMoveType::Travel : EMoveType::Extrude;
        move.time[0] = 1.;
        move.width = .4f;
        move.height = .2f;
        move.extrusion_role = erPerimeter;
        move.mm3_per_mm = .08;
        if (i == 0) continue;
        move.delta_extruder = i == 1 ? .01f : 1e-9f;
        move.deposition_source_id = owner;
        move.deposition_command_section_id = i == 1 ? bead_id : collapsed_id;
        move.deposition_design_heights = i == 1 ? std::array<double, 2>{.2, .2} : std::array<double, 2>{0., 0.};
        move.deposition_design_width = .4f;
        source.deposition_z_spans.push_back({i, owner, move.deposition_command_section_id, {.2, .2}, {.2f, .2f}});
    }
    source.print_statistics.modes[0].time = 3.;
    const auto path = timed_deposition(source);
    REQUIRE(path.state == ToolpathState::Complete);
    REQUIRE(path.unresolved_deposition.size() == 1);
    SectionGCodeScenario scenario(BeadSection::RectangularEnvelope);
    scenario.acknowledge_modeled_timing = true;
    scenario.thermal_materials = {{Eigen::Matrix3d::Identity() * .2, 2e6, 400.}};
    scenario.expansion = {{{{250., 450., Eigen::Matrix3d::Identity() * 1e-5}}, Eigen::Matrix3d::Identity(), 400.}};
    scenario.mechanical_materials = {{Eigen::Matrix<double, 6, 6>::Identity() * 1000.}};
    scenario.filament_to_material = {0};
    scenario.spacing_m = {.001, .001, .001};
    scenario.activation_times_s = {path.total_s};
    scenario.deposition_temperatures_k = [](double) { return std::vector<double>{400.}; };
    scenario.exposure = [](double, const SectionThermalAssembly&) { return std::vector<ThermalHistory::Element>{}; };
    CHECK_THROWS_AS(simulate_section_gcode_history(source, scenario), std::invalid_argument);
    scenario.collapsed_envelope_allowance = CollapsedEnvelopeAllowance{1e-12, 1e-11, 1e-6, .1};
    const auto nominal = resolve_nominal_deposition(source, path, path.total_s, scenario.section, {}, {},
        scenario.collapsed_envelope_allowance);
    const auto history = simulate_section_gcode_history(source, scenario);
    REQUIRE(history.thermal.state == ThermalHistory::State::Converged);
    REQUIRE(history.mechanical);
    CHECK(history.collapsed_envelope.count == 1);
    CHECK(history.collapsed_envelope.commanded_volume_m3 > 0.);
    CHECK_THAT(history.collapsed_envelope.commanded_volume_m3,
        WithinAbs(nominal.collapsed_envelope.commanded_volume_m3, 1e-24));
    CHECK(history.stationary_extrusion.count == 0);
}

TEST_CASE("G-code prefixes drive section heat history and mechanical recovery without rewriting the source", "[SectionGCodeHistory]")
{
    const bool was_bbl = GCodeProcessor::s_IsBBLPrinter;
    const ScopeGuard restore([was_bbl] { GCodeProcessor::s_IsBBLPrinter = was_bbl; });
    GCodeProcessor::s_IsBBLPrinter = false;
    ScopedTemporaryFile input(".gcode");
    {
        boost::nowide::ofstream stream(input.string());
        stream << "G90\nM83\n;TYPE:Outer wall\n;WIDTH:2\n;HEIGHT:1\n"
                  "G1 X1 Y1 Z1 F600\nG1 X3 Y1 E1.2 F600\nG1 X5 Y1 E1.2 F600\n";
        REQUIRE(stream.good());
    }
    GCodeProcessor processor;
    processor.apply_config(FullPrintConfig{});
    processor.process_file(input.string());
    const auto& source = processor.get_result();
    const auto path = timed_deposition(source);
    REQUIRE(path.state == ToolpathState::Complete);
    REQUIRE(path.deposition.size() == 2);
    // Extrusion volume must support a bead at least as wide as its height
    // when the same parsed source is tested with a rounded cross-section.
    for (const auto& entry : path.deposition) REQUIRE(entry.geometry.width_m>=entry.geometry.height_m);
    const auto move_count = source.moves.size();
    SectionGCodeScenario scenario(BeadSection::RectangularEnvelope);
    scenario.acknowledge_modeled_timing = true;
    scenario.thermal_materials = {{Eigen::Matrix3d::Identity()*.2,2e6,400.}};
    scenario.expansion = {{{{250.,450.,Eigen::Matrix3d::Identity()*1e-5}},Eigen::Matrix3d::Identity(),400.}};
    scenario.mechanical_materials = {{Eigen::Matrix<double,6,6>::Identity()*1000.}};
    scenario.filament_to_material = {0};
    scenario.spacing_m = {.001,.001,.001};
    scenario.activation_times_s = {path.deposition.front().end_s,path.total_s};
    scenario.deposition_temperatures_k = [](double) { return std::vector<double>{400.}; };
    scenario.exposure = [](double,const SectionThermalAssembly&) { return std::vector<ThermalHistory::Element>{}; };
    const auto result = simulate_section_gcode_history(source,scenario);
    REQUIRE(result.thermal.state == ThermalHistory::State::Converged);
    REQUIRE(result.mechanical.has_value());
    CHECK(result.activated_snapshots == 2);
    CHECK(result.uses_batched_activation);
    CHECK(result.uses_modeled_timing);
    CHECK(result.deposited_overlap_policy == DepositedOverlapPolicy::None);
    CHECK_THAT(result.thermal.final_energy_j,WithinAbs(result.occupied_volume_m3*2e6*400.,1e-8));
    CHECK_THAT(result.thermal.energy_balance_error_j,WithinAbs(0.,1e-8));
    CHECK_THAT(result.mechanical->elastic_energy_j,WithinAbs(0.,1e-12));
    CHECK(source.moves.size() == move_count);
    for (double temperature : result.thermal.temperatures_k) CHECK_THAT(temperature,WithinAbs(400.,1e-7));
    // Boundary-only events do not rescan prefixes or deposit material. Initial
    // and coincident changes are folded into the next activation's full law.
    scenario.cooldown_s = .1;
    scenario.exposure_times_s = {0., scenario.activation_times_s.front(),
        (scenario.activation_times_s.front()+path.total_s)*.5, path.total_s, path.total_s+.05};
    std::vector<double> exposure_calls;
    scenario.exposure = [&](double time,const SectionThermalAssembly&) {
        exposure_calls.push_back(time);
        return std::vector<ThermalHistory::Element>{};
    };
    SectionGCodeLimits boundary_limits;
    boundary_limits.max_prefix_visits = 2*(path.deposition.size()+path.unresolved_deposition.size());
    const auto changed = simulate_section_gcode_history(source,scenario,boundary_limits);
    REQUIRE(changed.thermal.state == ThermalHistory::State::Converged);
    CHECK(changed.activated_snapshots == 2);
    CHECK(changed.exposure_changes == 2);
    CHECK(exposure_calls.size() == 4);
    CHECK_THAT(changed.thermal.deposited_energy_j,WithinAbs(result.thermal.deposited_energy_j,1e-8));
    scenario.exposure_times_s = {.1,.1};
    CHECK_THROWS_AS(simulate_section_gcode_history(source,scenario),std::invalid_argument);
    scenario.exposure_times_s = {path.total_s+scenario.cooldown_s+1.};
    CHECK_THROWS_AS(simulate_section_gcode_history(source,scenario),std::invalid_argument);
    scenario.exposure_times_s.clear();
    // Explicit synthetic airflow, insulated bed. This checks source-to-heat-to-
    // stress coupling, not a measured fan coefficient or real adhesion margin.
    scenario.cooldown_s = .1;
    scenario.exposure = [](double,const SectionThermalAssembly& assembly) {
        std::vector<ThermalHistory::Element> boundaries;
        const auto add = [&](const auto& face) {
            const auto robin = SectionIntegration::thermal_boundary_element(face.area_moments,20.,300.,assembly.quadrature);
            ThermalHistory::Element boundary;
            boundary.nodes = face.nodes;
            boundary.conductance = robin.conductance;
            boundary.load = robin.load;
            boundary.boundary = true;
            boundary.kind = ThermalHistory::BoundaryKind::PartFan;
            boundaries.push_back(boundary);
        };
        for (const auto& face : assembly.horizontal_faces)
            if (face.kind != SectionHorizontalFaceKind::Bed) add(face);
        for (const auto& face : assembly.vertical_faces) add(face);
        return boundaries;
    };
    const auto cooled = simulate_section_gcode_history(source,scenario);
    REQUIRE(cooled.thermal.state == ThermalHistory::State::Converged);
    REQUIRE(cooled.mechanical.has_value());
    SectionGCodeLimits lumped_limits;
    lumped_limits.assembly.quadrature=SectionIntegration::ThermalQuadrature::LumpedDiagonal;
    const auto lumped=simulate_section_gcode_history(source,scenario,lumped_limits);
    REQUIRE(lumped.thermal.state==ThermalHistory::State::Converged);
    REQUIRE(lumped.mechanical.has_value());
    CHECK(lumped.thermal_quadrature==SectionIntegration::ThermalQuadrature::LumpedDiagonal);
    CHECK(cooled.thermal_quadrature==SectionIntegration::ThermalQuadrature::Consistent);
    for (const double temperature : lumped.thermal.temperatures_k) {
        CHECK(temperature>=300.-1e-8);
        CHECK(temperature<=400.+1e-8);
    }
    CHECK(cooled.thermal.boundary_energy_j < 0.);
    CHECK(cooled.mechanical->elastic_energy_j > 0.);
    CHECK_THAT(cooled.occupied_volume_m3,WithinAbs(result.occupied_volume_m3,1e-18));
    CHECK_THAT(cooled.thermal.energy_balance_error_j,WithinAbs(0.,1e-8));
    scenario.stream_tile_size_m=Eigen::Vector2d{.002,.002};
    scenario.stream_section_height_m=.0002;
    const auto streamed=simulate_section_gcode_history(source,scenario);
    REQUIRE(streamed.thermal.state==ThermalHistory::State::Converged);
    REQUIRE(streamed.mechanical.has_value());
    CHECK_THAT(streamed.occupied_volume_m3,WithinAbs(cooled.occupied_volume_m3,1e-17));
    CHECK_THAT(streamed.thermal.final_energy_j,WithinAbs(cooled.thermal.final_energy_j,1e-8));
    CHECK(streamed.mechanical->elastic_energy_j>0.);
    scenario.stream_vertical_normal_bins=32;
    const auto angular=simulate_section_gcode_history(source,scenario);
    CHECK(angular.vertical_normal_bins==32);
    REQUIRE(angular.thermal.state==ThermalHistory::State::Converged);
    CHECK_THAT(angular.thermal.final_energy_j,WithinAbs(streamed.thermal.final_energy_j,1e-8));
    scenario.stream_vertical_normal_bins=7;
    CHECK_THROWS_AS(simulate_section_gcode_history(source,scenario),std::invalid_argument);
    scenario.stream_vertical_normal_bins=0;
    // Exercise fixed-reference streamed assembly with actual parsed paths.
    // A partial timeline must not activate the final path merely because its
    // geometry is already available as the immutable reference.
    const auto nominal_final=resolve_nominal_deposition(source,path,path.total_s,BeadSection::RectangularEnvelope);
    const auto nominal_first=resolve_nominal_deposition(source,path,path.deposition.front().end_s,BeadSection::RectangularEnvelope);
    ThermalSectionStreamSettings fixed_settings;
    fixed_settings.spacing_m=scenario.spacing_m;
    fixed_settings.tile_size_m=*scenario.stream_tile_size_m;
    fixed_settings.maximum_height_m=scenario.stream_section_height_m;
    const auto fixed_first=assemble_reference_thermal_volume(nominal_final.segments,nominal_first.segments,
        scenario.thermal_materials,fixed_settings);
    const auto fixed_final=assemble_reference_thermal_volume(nominal_final.segments,nominal_final.segments,
        scenario.thermal_materials,fixed_settings);
    CHECK(fixed_first.assembly.exposure_complete);
    CHECK(fixed_first.assembly.occupied_volume_m3<fixed_final.assembly.occupied_volume_m3);
    CHECK_THAT(fixed_final.assembly.occupied_volume_m3,WithinAbs(streamed.occupied_volume_m3,1e-17));
    CHECK_THAT(fixed_final.excluded_nominal_volume_m3,WithinAbs(0.,1e-17));
    scenario.use_fixed_interface_reference=true;
    const auto fixed_history=simulate_section_gcode_history(source,scenario);
    REQUIRE(fixed_history.thermal.state==ThermalHistory::State::Converged);
    REQUIRE(fixed_history.mechanical.has_value());
    CHECK(fixed_history.uses_fixed_interface_reference);
    CHECK(fixed_history.mechanical->elastic_energy_j>0.);
    CHECK_THAT(fixed_history.thermal.final_energy_j,WithinAbs(streamed.thermal.final_energy_j,1e-8));
    CHECK_THAT(fixed_history.occupied_volume_m3,WithinAbs(streamed.occupied_volume_m3,1e-17));
    CHECK_THAT(fixed_history.excluded_nominal_volume_m3,WithinAbs(0.,1e-17));
    SectionGCodeLimits partitioned_limits;
    partitioned_limits.streaming_geometry.partition_increment_m=.00001;
    const auto partitioned=simulate_section_gcode_history(source,scenario,partitioned_limits);
    REQUIRE(partitioned.thermal.state==ThermalHistory::State::Converged);
    REQUIRE(partitioned.mechanical.has_value());
    CHECK_THAT(partitioned.section_partition_increment_m,WithinAbs(.00001,1e-15));
    CHECK_THAT(partitioned.occupied_volume_m3,WithinAbs(fixed_history.occupied_volume_m3,1e-17));
    CHECK_THAT(partitioned.thermal.final_energy_j,WithinAbs(fixed_history.thermal.final_energy_j,1e-8));
    CHECK_THAT(partitioned.mechanical->elastic_energy_j,WithinAbs(fixed_history.mechanical->elastic_energy_j,1e-12));
    scenario.use_fixed_interface_reference=false;
    // Tile size controls memory partitioning, not the physical discretization.
    // Changing it must preserve both the thermal history and recovered stress.
    scenario.stream_tile_size_m=Eigen::Vector2d{.003,.002};
    const auto retiled=simulate_section_gcode_history(source,scenario);
    REQUIRE(retiled.thermal.state==ThermalHistory::State::Converged);
    REQUIRE(retiled.mechanical.has_value());
    CHECK_THAT(retiled.occupied_volume_m3,WithinAbs(streamed.occupied_volume_m3,1e-17));
    CHECK_THAT(retiled.thermal.final_energy_j,WithinAbs(streamed.thermal.final_energy_j,1e-8));
    CHECK_THAT(retiled.mechanical->elastic_energy_j,WithinAbs(streamed.mechanical->elastic_energy_j,1e-14));
    // Exercise adaptive time integration through the full streamed adapter,
    // independently of the geometry tile and activation schedule.
    SectionGCodeLimits adaptive_limits;
    adaptive_limits.history.max_step_temperature_error_k=1e-4;
    const auto adaptive=simulate_section_gcode_history(source,scenario,adaptive_limits);
    REQUIRE(adaptive.thermal.state==ThermalHistory::State::Converged);
    REQUIRE(adaptive.mechanical.has_value());
    CHECK(adaptive.thermal.temporal_error_controlled);
    CHECK(adaptive.thermal.maximum_step_temperature_error_k<=1e-4);
    CHECK_THAT(adaptive.thermal.energy_balance_error_j,WithinAbs(0.,1e-8));
    scenario.stream_tile_size_m=Eigen::Vector2d{.002,.002};
    scenario.section=BeadSection::RoundedRectangle;
    std::size_t ownership_queries=0;
    scenario.stream_ownership=[&](double,double) -> std::optional<SectionOwnership> { ++ownership_queries; return {}; };
    const auto rounded=simulate_section_gcode_history(source,scenario);
    REQUIRE(rounded.thermal.state==ThermalHistory::State::Converged);
    REQUIRE(rounded.mechanical.has_value());
    CHECK(rounded.occupied_volume_m3<streamed.occupied_volume_m3);
    CHECK(rounded.thermal.boundary_energy_j<0.);
    CHECK(rounded.mechanical->elastic_energy_j>0.);
    CHECK(source.moves.size()==move_count);
    CHECK(ownership_queries>0);
    scenario.use_native_ownership=true;
    CHECK_THROWS_AS(simulate_section_gcode_history(source,scenario),std::invalid_argument);
    scenario.stream_ownership={};
    // A parsed file has no authenticated native generation snapshot. Never
    // substitute filament labels or synthetic masks for that missing evidence.
    CHECK_THROWS_AS(simulate_section_gcode_history(source,scenario),std::invalid_argument);
    scenario.use_native_ownership=false;
    scenario.stream_tile_size_m.reset();
    scenario.section=BeadSection::RectangularEnvelope;
    scenario.acknowledge_modeled_timing = false;
    CHECK_THROWS_AS(simulate_section_gcode_history(source,scenario),std::invalid_argument);
    scenario.acknowledge_modeled_timing = true;
    SectionGCodeLimits limited;
    limited.max_prefix_visits = 1;
    CHECK_THROWS_AS(simulate_section_gcode_history(source,scenario,limited),ResourceLimit);
    limited = {};
    limited.history.max_cell_steps = 1;
    const auto failed = simulate_section_gcode_history(source,scenario,limited);
    CHECK(failed.thermal.state == ThermalHistory::State::BudgetExceeded);
    CHECK_FALSE(failed.mechanical.has_value());
    scenario.activation_times_s.pop_back();
    CHECK_THROWS_AS(simulate_section_gcode_history(source,scenario),std::invalid_argument);
}
