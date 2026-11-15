#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/SectionThermalSequence.hpp"
#include "libslic3r/SectionGrowth.hpp"
#include "libslic3r/SectionThermalAssembly.hpp"
#include "libslic3r/SectionThermalLoad.hpp"
#include "libslic3r/SectionIntegration.hpp"
#include "libslic3r/ThermalDeposition.hpp"
#include "libslic3r/ThermalHistory.hpp"
#include "libslic3r/ThermoelasticSectionTopology.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticLayers.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/libslic3r.h"
#include <Eigen/Core>
#include <cstddef>
#include <stdexcept>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;
namespace {
ExPolygon rectangle(coord_t left,coord_t right)
{
    return ExPolygon(Polygon{{left,0},{right,0},{right,1000000000},{left,1000000000}});
}
std::vector<ThermalHistory::Element> insulated(const SectionThermalAssembly&) { return {}; }
}

TEST_CASE("Streamed activation preserves old heat through component joins and unchanged snapshots", "[SectionThermalSequence]")
{
    const std::vector<DepositedSegment> split{
        {{0.,.0005,.001},{.001,.0005,.001},.001,.001,0,BeadSection::RectangularEnvelope},
        {{.002,.0005,.001},{.003,.0005,.001},.001,.001,0,BeadSection::RectangularEnvelope}};
    auto joined=split;
    joined.push_back({{.001,.0005,.001},{.002,.0005,.001},.001,.001,0,BeadSection::RectangularEnvelope});
    ThermalSectionStreamSettings settings;
    settings.spacing_m={.001,.001,.001}; settings.tile_size_m={.001,.001}; settings.maximum_height_m=.0002;
    SectionThermalStreamSequence sequence({{Eigen::Matrix3d::Identity()*.2,2e6,400.}},settings);
    const auto first=sequence.append(split,0.,{400.},insulated);
    const double first_volume=sequence.assembly().occupied_volume_m3;
    const auto second=sequence.append(joined,1.,{500.},insulated);
    const double volume=sequence.assembly().occupied_volume_m3;
    REQUIRE(second.coalesced_nodes.has_value());
    const auto unchanged=sequence.append(joined,2.,{999.},insulated);
    CHECK(unchanged.added_mass.empty());
    CHECK(unchanged.deposited_cells.empty());
    CHECK(unchanged.added_elements.empty());
    auto moved=joined;
    for (auto& p : moved) { p.start_m.x()+=.01; p.end_m.x()+=.01; }
    CHECK_THROWS_AS(sequence.append(moved,3.,{500.},insulated),std::invalid_argument);
    // Failed append must not advance time or replace the prior source snapshot.
    const auto recovered=sequence.append(joined,2.5,{999.},insulated);
    const auto history=ThermalHistory::evolve({first,second,unchanged,recovered},3.);
    REQUIRE(history.state==ThermalHistory::State::Converged);
    const double expected=first_volume*2e6*400.+(volume-first_volume)*2e6*500.;
    CHECK_THAT(history.final_energy_j,WithinAbs(expected,1e-8));
    CHECK_THAT(history.deposited_energy_j,WithinAbs(expected,1e-8));
    CHECK_THAT(history.energy_balance_error_j,WithinAbs(0.,1e-8));
    REQUIRE(history.temperatures_k.size()==sequence.assembly().nodes.size());
    for (const auto t : history.temperatures_k) { CHECK(t>=399.99999); CHECK(t<=500.00001); }
}

TEST_CASE("Geometry driven thermal activation joins components without redepositing old material", "[SectionThermalSequence]")
{
    const auto split = qualify_section_topology({{0.,1.,0,{rectangle(0,400000000),rectangle(600000000,1000000000)}}},0.);
    const auto joined = qualify_section_topology({{0.,1.,0,{rectangle(0,1000000000)}}},0.);
    SectionThermalSequence sequence({{Eigen::Matrix3d::Identity(),1.,400.}},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones());
    const auto first = sequence.append(split,0.,{400.},insulated);
    REQUIRE(first.deposited_cells.size() == 16);
    const auto second = sequence.append(joined,1.,{300.},insulated);
    REQUIRE(second.coalesced_nodes.has_value());
    CHECK(second.deposited_cells.empty());
    const auto result = ThermalHistory::evolve({first,second},1.);
    REQUIRE(result.state == ThermalHistory::State::Converged);
    REQUIRE(result.temperatures_k.size() == sequence.assembly().nodes.size());
    CHECK(result.temperatures_k.size() == 8);
    CHECK_THAT(result.deposited_energy_j,WithinAbs(380.,1e-9));
    CHECK_THAT(result.final_energy_j,WithinAbs(380.,1e-9));
    CHECK_THAT(result.energy_balance_error_j,WithinAbs(0.,1e-9));
    for (double t : result.temperatures_k) CHECK_THAT(t,WithinAbs(380.,1e-9));
    CHECK_NOTHROW(section_thermal_load(sequence.assembly(),result,
        {{{{250.,450.,Eigen::Matrix3d::Identity()*1e-5}},Eigen::Matrix3d::Identity(),400.}}));
}

TEST_CASE("Exposure changes preserve active geometry and deposited heat without consuming snapshots", "[SectionThermalSequence]")
{
    const auto geometry = qualify_section_topology({{0.,1.,0,{rectangle(0,1000000000)}}},0.);
    const auto cooling = [](const SectionThermalAssembly& assembly) {
        std::vector<ThermalHistory::Element> result;
        for (const auto& face : assembly.horizontal_faces) {
            const auto integrated = SectionIntegration::thermal_boundary_element(face.area_moments,1.,300.);
            ThermalHistory::Element boundary;
            boundary.nodes = face.nodes;
            boundary.conductance = integrated.conductance;
            boundary.load = integrated.load;
            boundary.boundary = true;
            result.push_back(boundary);
        }
        return result;
    };
    const auto verify = [&](auto& sequence, const auto& input) {
        CHECK_THROWS_AS(sequence.replace_exposure(0.,cooling),std::logic_error);
        const auto first = sequence.append(input,0.,{400.},insulated);
        const auto* nodes = sequence.assembly().nodes.data();
        const auto count = sequence.assembly().nodes.size();
        const auto volume = sequence.assembly().occupied_volume_m3;
        const auto on = sequence.replace_exposure(1.,cooling);
        CHECK(on.deposited_cells.empty());
        CHECK(on.added_mass.empty());
        CHECK(on.added_elements.empty());
        CHECK(on.added_connections.empty());
        CHECK_FALSE(on.coalesced_nodes);
        REQUIRE_FALSE(on.exposed_elements.empty());
        CHECK_THROWS_AS(sequence.replace_exposure(.5,cooling),std::invalid_argument);
        CHECK_THROWS_AS(sequence.replace_exposure(3.,[](const SectionThermalAssembly&) {
            return std::vector<ThermalHistory::Element>(1); // Not a boundary operator.
        }),std::invalid_argument);
        const auto off = sequence.replace_exposure(2.,insulated);
        CHECK(off.exposed_elements.empty());
        CHECK(sequence.assembly().nodes.data() == nodes);
        CHECK(sequence.assembly().nodes.size() == count);
        CHECK_THAT(sequence.assembly().occupied_volume_m3,WithinAbs(volume,1e-15));
        const auto history = ThermalHistory::evolve({first,on,off},3.);
        REQUIRE(history.state == ThermalHistory::State::Converged);
        CHECK(history.final_energy_j < history.deposited_energy_j);
        CHECK_THAT(history.deposited_energy_j,WithinAbs(volume*400.,1e-7));
        CHECK_THAT(history.energy_balance_error_j,WithinAbs(0.,1e-7));
    };
    SectionThermalSequence retained({{Eigen::Matrix3d::Identity(),1.,400.}},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),{},{},1);
    verify(retained,geometry);
    ThermalSectionStreamSettings settings;
    settings.spacing_m=Eigen::Vector3d::Ones(); settings.tile_size_m={1.,1.}; settings.maximum_height_m=1.;
    SectionThermalStreamSequence streamed({{Eigen::Matrix3d::Identity(),1.,400.}},settings,{},1);
    const std::vector<DepositedSegment> segments{{{0.,.5,1.},{1.,.5,1.},1.,1.,0,BeadSection::RectangularEnvelope}};
    verify(streamed,segments);
}

TEST_CASE("Fixed interface activation adds foreign material without rewriting deposited heat", "[SectionThermalSequence]")
{
    const auto reference = qualify_section_topology({{0.,1.,0,{rectangle(0,500000000)}},
        {0.,1.,1,{rectangle(500000000,1000000000)}}},0.);
    const auto a = qualify_section_topology({{0.,1.,0,{rectangle(0,600000000)}}},0.);
    const auto b = qualify_section_topology({{0.,1.,1,{rectangle(400000000,1000000000)}}},0.);
    const auto first_geometry = activate_reference_sections(reference,{a});
    const auto final_geometry = activate_reference_sections(reference,{a,b});
    SectionThermalSequence sequence({{Eigen::Matrix3d::Identity(),2.,400.},
        {Eigen::Matrix3d::Identity(),3.,500.}},Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones());
    const auto first = sequence.append(first_geometry.geometry,0.,{400.,500.},insulated);
    const auto second = sequence.append(final_geometry.geometry,1.,{999.,500.},insulated);
    const auto repeat = sequence.append(final_geometry.geometry,2.,{999.,999.},insulated);
    const auto history = ThermalHistory::evolve({first,second,repeat},3.);
    REQUIRE(history.state == ThermalHistory::State::Converged);
    const double expected = .5*2.*400.+.5*3.*500.;
    CHECK_THAT(history.final_energy_j,WithinAbs(expected,1e-8));
    CHECK_THAT(history.deposited_energy_j,WithinAbs(expected,1e-8));
    CHECK_THAT(history.energy_balance_error_j,WithinAbs(0.,1e-8));
    CHECK(repeat.added_mass.empty());
    CHECK(repeat.deposited_cells.empty());
}

TEST_CASE("Streamed fixed interfaces preserve chronological heat when foreign paths overlap", "[SectionThermalSequence]")
{
    const std::vector<DepositedSegment> reference{
        {{0.,.0005,.001},{.003,.0005,.001},.001,.001,0,BeadSection::RectangularEnvelope},
        {{0.,.0015,.001},{.003,.0015,.001},.001,.001,1,BeadSection::RectangularEnvelope}};
    auto deposited=reference;
    for (auto& path : deposited) path.width_m=.0014;
    ThermalSectionStreamSettings settings;
    settings.spacing_m={.001,.001,.001}; settings.tile_size_m={.001,.001}; settings.maximum_height_m=.0002;
    SectionThermalStreamSequence sequence({{Eigen::Matrix3d::Identity()*.2,2e6,400.},
        {Eigen::Matrix3d::Identity()*.1,3e6,500.}},settings,{},1000,reference);
    const auto first=sequence.append({deposited[0]},0.,{400.,500.},insulated);
    const double first_volume=sequence.assembly().occupied_volume_m3;
    CHECK(sequence.excluded_nominal_volume_m3()>0.);
    const auto second=sequence.append(deposited,1.,{999.,500.},insulated);
    const double final_volume=sequence.assembly().occupied_volume_m3;
    const double excluded=sequence.excluded_nominal_volume_m3();
    const auto repeat=sequence.append(deposited,2.,{999.,999.},insulated);
    CHECK(repeat.added_mass.empty());
    CHECK(repeat.deposited_cells.empty());
    CHECK_THAT(sequence.excluded_nominal_volume_m3(),WithinAbs(excluded,1e-18));
    CHECK_THROWS_AS(sequence.append({deposited[0]},3.,{400.,500.},insulated),std::invalid_argument);
    const auto recovered=sequence.append(deposited,2.5,{999.,999.},insulated);
    const auto history=ThermalHistory::evolve({first,second,repeat,recovered},3.);
    REQUIRE(history.state==ThermalHistory::State::Converged);
    const double expected=first_volume*2e6*400.+(final_volume-first_volume)*3e6*500.;
    CHECK_THAT(history.final_energy_j,WithinAbs(expected,1e-8));
    CHECK_THAT(history.deposited_energy_j,WithinAbs(expected,1e-8));
    CHECK_THAT(history.energy_balance_error_j,WithinAbs(0.,1e-8));
}

TEST_CASE("Geometry driven thermal growth preserves old cooling and replaces covered exposure", "[SectionThermalSequence]")
{
    const auto lower = qualify_section_topology({{0.,1.,0,{rectangle(0,1000000000)}}},0.);
    const auto taller = qualify_section_topology({{0.,2.,0,{rectangle(0,1000000000)}}},0.);
    std::vector<double> areas;
    const auto cooling = [&](const SectionThermalAssembly& assembly) {
        std::vector<ThermalHistory::Element> elements;
        double area = 0.;
        const auto add = [&](const auto& face) {
            area += face.area_moments[0];
            const auto integrated = SectionIntegration::thermal_boundary_element(face.area_moments,1.,300.);
            ThermalHistory::Element boundary;
            boundary.nodes = face.nodes;
            boundary.conductance = integrated.conductance;
            boundary.load = integrated.load;
            boundary.boundary = true;
            elements.push_back(boundary);
        };
        for (const auto& face : assembly.horizontal_faces) add(face);
        for (const auto& face : assembly.vertical_faces) add(face);
        areas.push_back(area);
        return elements;
    };
    SectionThermalSequence sequence({{Eigen::Matrix3d::Identity(),1.,400.}},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones());
    const auto first = sequence.append(lower,0.,{400.},cooling);
    const auto second = sequence.append(taller,1.,{400.},cooling);
    REQUIRE(areas.size() == 2);
    CHECK_THAT(areas[0],WithinAbs(6.,1e-12));
    CHECK_THAT(areas[1],WithinAbs(10.,1e-12));
    const auto result = ThermalHistory::evolve({first,second},1.);
    REQUIRE(result.state == ThermalHistory::State::Converged);
    CHECK_THAT(result.deposited_energy_j,WithinAbs(800.,1e-8));
    CHECK_THAT(result.final_energy_j,WithinAbs(5000./7.,1e-8));
    CHECK_THAT(result.boundary_energy_j,WithinAbs(-600./7.,1e-8));
    CHECK_THAT(result.energy_balance_error_j,WithinAbs(0.,1e-8));
}

TEST_CASE("Rejected thermal geometry transitions leave the prior sequence intact", "[SectionThermalSequence]")
{
    const auto full = qualify_section_topology({{0.,1.,0,{rectangle(0,1000000000)}}},0.);
    const auto smaller = qualify_section_topology({{0.,1.,0,{rectangle(0,500000000)}}},0.);
    SectionThermalSequence sequence({{Eigen::Matrix3d::Identity(),1.,400.}},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),{},{},2);
    CHECK_THROWS_AS(sequence.assembly(),std::logic_error);
    const auto first = sequence.append(full,0.,{400.},insulated);
    CHECK_THROWS_AS(sequence.append(smaller,1.,{400.},insulated),std::invalid_argument);
    CHECK_THROWS_AS(sequence.append(full,-1.,{400.},insulated),std::invalid_argument);
    const auto unchanged = sequence.append(full,1.,{300.},insulated);
    CHECK(unchanged.added_mass.empty());
    CHECK(unchanged.deposited_cells.empty());
    CHECK(unchanged.added_elements.empty());
    const auto result = ThermalHistory::evolve({first,unchanged},1.);
    REQUIRE(result.state == ThermalHistory::State::Converged);
    CHECK_THAT(result.final_energy_j,WithinAbs(400.,1e-9));
    CHECK_THROWS_AS(sequence.append(full,2.,{400.},insulated),ResourceLimit);
}

TEST_CASE("Thermal snapshot activation mixes material capacities without resetting earlier temperatures", "[SectionThermalSequence]")
{
    const auto lower = qualify_section_topology({{0.,1.,0,{rectangle(0,1000000000)}}},0.);
    const auto both = qualify_section_topology({{0.,1.,0,{rectangle(0,1000000000)}},
        {1.,2.,1,{rectangle(0,1000000000)}}},0.);
    SectionThermalSequence sequence({{Eigen::Matrix3d::Identity(),2.,400.},
        {Eigen::Matrix3d::Identity(),4.,300.}},Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones());
    unsigned cursor = 0;
    const auto result = ThermalHistory::evolve_generated([&](ThermalHistory::DepositionEvent& event) {
        if (cursor == 2) return false;
        // The material-0 deposition target changes, but none of that material
        // is added at t=1; its previously deposited mass must NOT be reheated.
        event = cursor == 0 ? sequence.append(lower,0.,{400.,300.},insulated) :
            sequence.append(both,1.,{999.,300.},insulated);
        ++cursor;
        return true;
    },1.);
    REQUIRE(result.state == ThermalHistory::State::Converged);
    CHECK_THAT(result.deposited_energy_j,WithinAbs(2000.,1e-8));
    CHECK_THAT(result.energy_balance_error_j,WithinAbs(0.,1e-8));
    REQUIRE(result.temperatures_k.size() == sequence.assembly().nodes.size());
    for (std::size_t i = 0; i < result.temperatures_k.size(); ++i) {
        const double z = sequence.assembly().nodes[i].position_m.z();
        const double expected = z < .5 ? 400. : z < 1.5 ? 1000./3. : 300.;
        CHECK_THAT(result.temperatures_k[i],WithinAbs(expected,1e-9));
    }
    const SectionExpansionMaterial expansion{{{250.,450.,Eigen::Matrix3d::Identity()*1e-5}},
        Eigen::Matrix3d::Identity(),400.};
    CHECK_NOTHROW(section_thermal_load(sequence.assembly(),result,{expansion,expansion}));
}
