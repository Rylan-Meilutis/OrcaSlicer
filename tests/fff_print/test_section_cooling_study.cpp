#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "libslic3r/SectionCoolingStudy.hpp"
#include "libslic3r/SectionGCodeHistory.hpp"
#include "libslic3r/SectionThermalAssembly.hpp"
#include "libslic3r/ThermalHistory.hpp"
#include "libslic3r/SectionIntegration.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/ThermoelasticToolpath.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticConflictOwnership.hpp"
#include "libslic3r/Utils.hpp"
#include "test_utils.hpp"
#include <boost/nowide/fstream.hpp>
#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;

TEST_CASE("Paired cooling studies isolate boundary laws and retain material accounting", "[SectionCoolingStudy]")
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
    SectionGCodeScenario common(BeadSection::RectangularEnvelope);
    common.acknowledge_modeled_timing = true;
    // Synthetic material/law constants test coupling, not ASA or adhesion.
    common.thermal_materials = {{Eigen::Matrix3d::Identity()*.2,2e6,400.}};
    common.expansion = {{{{250.,450.,Eigen::Matrix3d::Identity()*1e-5}},Eigen::Matrix3d::Identity(),400.}};
    common.mechanical_materials = {{Eigen::Matrix<double,6,6>::Identity()*1000.}};
    common.filament_to_material = {0};
    common.spacing_m = {.001,.001,.001};
    common.activation_times_s = {path.total_s};
    common.cooldown_s = 1.;
    common.deposition_temperatures_k = [](double) { return std::vector<double>{400.}; };
    const auto insulated = [](double,const SectionThermalAssembly&) { return std::vector<ThermalHistory::Element>{}; };
    std::array<SectionCoolingCase,2> cases{{{"baseline",{},insulated},{"candidate",{},insulated}}};
    const auto same = compare_section_cooling(source,common,cases);
    CHECK(same[0].label == "baseline");
    CHECK(same[1].label == "candidate");
    for (const auto& metrics:same) {
        CHECK_FALSE(metrics.uses_fixed_interface_reference);
        CHECK_FALSE(metrics.allows_single_sided_native_boundary);
        CHECK_FALSE(metrics.allows_dual_native_boundary);
        CHECK(metrics.native_interface_variant==NativeInterfaceVariant::Strict);
    }
    SectionGCodeLimits invalid_policy;
    invalid_policy.ownership.allow_dual_native_boundary=true;
    CHECK_THROWS_WITH(compare_section_cooling(source,common,cases,invalid_policy),
        Catch::Matchers::ContainsSubstring("requires an explicit non-strict ownership policy"));
    invalid_policy.ownership.allow_dual_native_boundary=false;
    invalid_policy.ownership.allow_single_sided_native_boundary=true;
    CHECK_THROWS_WITH(compare_section_cooling(source,common,cases,invalid_policy),
        Catch::Matchers::ContainsSubstring("requires an explicit non-strict ownership policy"));
    CHECK_THAT(same[0].final_energy_j,WithinAbs(same[1].final_energy_j,1e-12));
    CHECK_THAT(same[0].elastic_energy_j,WithinAbs(same[1].elastic_energy_j,1e-12));
    cases[1].exposure_times_s = {path.total_s+.5};
    cases[1].exposure = [start=path.total_s+.5](double time,const SectionThermalAssembly& assembly) {
        std::vector<ThermalHistory::Element> output;
        if (time < start) return output;
        const auto add = [&](const auto& face) {
            const auto law = SectionIntegration::thermal_boundary_element(face.area_moments,20.,300.,assembly.quadrature);
            ThermalHistory::Element boundary;
            boundary.nodes=face.nodes; boundary.conductance=law.conductance; boundary.load=law.load; boundary.boundary=true;
            output.push_back(boundary);
        };
        for (const auto& face : assembly.horizontal_faces) add(face);
        for (const auto& face : assembly.vertical_faces) add(face);
        return output;
    };
    const auto changed = compare_section_cooling(source,common,cases);
    CHECK(changed[1].final_energy_j < changed[0].final_energy_j);
    CHECK(changed[1].elastic_energy_j > changed[0].elastic_energy_j);
    CHECK(changed[1].exposure_changes == 1);
    CHECK_THAT(changed[1].occupied_volume_m3,WithinAbs(changed[0].occupied_volume_m3,1e-20));
    CHECK_THAT(changed[1].deposited_energy_j,WithinAbs(changed[0].deposited_energy_j,1e-12));
    CHECK_THAT(changed[1].energy_balance_error_j,WithinAbs(0.,1e-8));
    CHECK(changed[1].collapsed_envelope.count == 0);
    CHECK(changed[1].stationary_extrusion.count == 0);
    CHECK_THAT(changed[1].modeled_observation_time_s,WithinAbs(path.total_s+common.cooldown_s,1e-12));
    CHECK_THROWS_AS(compare_section_cooling(source,common,cases,{},[] { throw std::runtime_error("cancel"); }),std::runtime_error);
    auto limited = SectionGCodeLimits{};
    limited.max_prefix_visits = 1;
    CHECK_THROWS_AS(compare_section_cooling(source,common,cases,limited),ResourceLimit);
    cases[1].exposure_times_s = {0.,0.};
    CHECK_THROWS_AS(compare_section_cooling(source,common,cases),std::invalid_argument);
    cases[1].exposure_times_s.clear();
    cases[1].label = cases[0].label;
    CHECK_THROWS_AS(compare_section_cooling(source,common,cases),std::invalid_argument);
    cases[1] = {"candidate",{},insulated};
    std::size_t temperature_queries = 0;
    common.deposition_temperatures_k = [&](double) {
        return std::vector<double>{++temperature_queries == 1 ? 400. : 390.};
    };
    // Stateful common inputs may not disguise changed deposition heat as cooling.
    CHECK_THROWS_AS(compare_section_cooling(source,common,cases),std::invalid_argument);
}
