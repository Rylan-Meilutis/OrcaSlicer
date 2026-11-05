#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/SectionThermalAssembly.hpp"
#include "libslic3r/SectionThermalExposure.hpp"
#include "libslic3r/SectionThermalLoad.hpp"
#include "libslic3r/SectionIntegration.hpp"
#include "libslic3r/SectionGrowth.hpp"
#include "libslic3r/ThermoelasticExactSections.hpp"
#include "libslic3r/ThermoelasticSectionFootprint.hpp"
#include "libslic3r/ThermalDeposition.hpp"
#include "libslic3r/ThermalHistory.hpp"
#include "libslic3r/ThermoelasticSectionTopology.hpp"
#include "libslic3r/ThermoelasticLayers.hpp"
#include "libslic3r/ThermoelasticLocalSections.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticSections.hpp"
#include "libslic3r/ThermoelasticThermalStrain.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/libslic3r.h"
#include <Eigen/Core>
#include <cstddef>
#include <cmath>
#include <limits>
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
SectionThermalMaterial material(double capacity,double temperature)
{
    return {Eigen::Matrix3d::Identity(),capacity,temperature};
}
}

TEST_CASE("Section exposure rejects detached geometry instead of treating it as covered", "[SectionThermalAssembly]")
{
    const auto topology = qualify_section_topology({{0.,1.,0,{rectangle(0,1000000000)}}},0.);
    const auto detached = topology.pieces().front();
    const auto charge = [](std::size_t) {};
    CHECK_THROWS_AS(exposed_vertical_section_patches(topology,detached,100,charge),std::invalid_argument);
    CHECK(exposed_vertical_section_patches(topology,topology.pieces().front(),100,charge).size() == 4);
}

TEST_CASE("Certified section geometry assembles material-specific thermal mass and shared nodes", "[SectionThermalAssembly]")
{
    const auto topology = qualify_section_topology({{0.,1.,0,{rectangle(0,1000000000)}},
        {1.,2.,1,{rectangle(0,1000000000)}}},0.);
    const auto assembly = assemble_section_thermal_volume(topology,{material(2.,400.),material(4.,300.)},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones());
    REQUIRE(assembly.nodes.size() == 12);
    REQUIRE(assembly.initial.deposited_cells.size() == 12);
    CHECK(assembly.initial.added_elements.size() == 2);
    CHECK_THAT(assembly.occupied_volume_m3,WithinAbs(2.,1e-12));
    double bed_area = 0., upper_area = 0., lower_area = 0.;
    for (const auto& face : assembly.horizontal_faces) {
        if (face.kind == SectionHorizontalFaceKind::Bed) bed_area += face.area_moments[0];
        else if (face.kind == SectionHorizontalFaceKind::Upward) upper_area += face.area_moments[0];
        else lower_area += face.area_moments[0];
    }
    CHECK_THAT(bed_area,WithinAbs(1.,1e-12));
    CHECK_THAT(upper_area,WithinAbs(1.,1e-12));
    CHECK_THAT(lower_area,WithinAbs(0.,1e-12));
    double heat = 0., capacity = 0.;
    for (const auto& cell : assembly.initial.deposited_cells) {
        capacity += cell.volume_m3*cell.volumetric_heat_capacity_j_m3_k;
        heat += cell.volume_m3*cell.volumetric_heat_capacity_j_m3_k*cell.temperature_k;
    }
    CHECK_THAT(capacity,WithinAbs(6.,1e-12));
    CHECK_THAT(heat,WithinAbs(2000.,1e-9));
    ThermalHistory::HistoryLimits limits;
    limits.solve.relative_tolerance = 1e-12;
    const auto history = ThermalHistory::evolve({assembly.initial},100.,limits);
    REQUIRE(history.state == ThermalHistory::State::Converged);
    CHECK_THAT(history.energy_balance_error_j,WithinAbs(0.,1e-7));
    CHECK_THAT(history.final_energy_j,WithinAbs(2000.,1e-7));
    for (double temperature : history.temperatures_k)
        CHECK_THAT(temperature,WithinAbs(2000./6.,.001));
}

TEST_CASE("Exact growth assembles thermal mass without rounding the deposited cut cell", "[SectionThermalAssembly]")
{
    const auto triangle = make_exact_section({ExPolygon(Polygon{{0,0},{10000001,0},{10000001,4000000}})});
    const auto strip = make_exact_section({ExPolygon(Polygon{{0,0},{1000000,0},{1000000,4000000},{0,4000000}})});
    const auto cut = boolean_exact_sections(triangle.section,strip.section,ExactSectionOperation::Intersection);
    const auto previous = qualify_section_topology({{0.,.001,0,SectionFootprint(cut.section)}},0.);
    const auto next = qualify_section_topology({{0.,.002,0,SectionFootprint(cut.section)}},0.);
    const auto growth = qualify_section_growth(previous,next);
    REQUIRE(growth.added_pieces().size() == 1);
    REQUIRE(growth.added_pieces().front().slab.footprint.is_exact());
    const auto assembly = assemble_section_thermal_addition(growth,{material(2000000.,400.)},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Constant(.001));
    const double expected_volume = .002 / 10000001.;
    CHECK_THAT(assembly.occupied_volume_m3,Catch::Matchers::WithinRel(expected_volume,1e-12));
    double capacity = 0., heat = 0.;
    for (const auto& cell : assembly.initial.deposited_cells) {
        capacity += cell.volume_m3 * cell.volumetric_heat_capacity_j_m3_k;
        heat += cell.volume_m3 * cell.volumetric_heat_capacity_j_m3_k * cell.temperature_k;
    }
    CHECK_THAT(capacity,Catch::Matchers::WithinRel(expected_volume * 2000000.,1e-12));
    CHECK_THAT(heat,Catch::Matchers::WithinRel(expected_volume * 2000000. * 400.,1e-12));
    CHECK(assembly.horizontal_faces.empty());
    CHECK(assembly.vertical_faces.empty());
    CHECK_FALSE(assembly.exposure_complete);
    const auto full=assemble_section_thermal_volume(next,{material(2000000.,400.)},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Constant(.001));
    CHECK(full.exposure_complete);
    CHECK_THAT(full.occupied_volume_m3,Catch::Matchers::WithinRel(2.*expected_volume,1e-12));
    double horizontal=0.,vertical=0.;
    Eigen::Vector3d normal_area=Eigen::Vector3d::Zero();
    for (const auto& face : full.horizontal_faces) horizontal+=face.area_moments[0];
    for (const auto& face : full.vertical_faces) {
        vertical+=face.area_moments[0];
        normal_area+=face.outward_normal*face.area_moments[0];
        CHECK_THAT(face.outward_normal.norm(),WithinAbs(1.,1e-14));
    }
    const double height=.004*1000000./10000001.;
    CHECK_THAT(horizontal,Catch::Matchers::WithinRel(.001*height,1e-12));
    CHECK_THAT(vertical,Catch::Matchers::WithinRel(.002*(.001+height+std::hypot(.001,height)),1e-12));
    CHECK_THAT(normal_area.norm(),WithinAbs(0.,1e-18));
    SectionThermalAssemblyLimits limits;
    limits.max_work = 100;
    CHECK_THROWS_AS(assemble_section_thermal_addition(growth,{material(2000000.,400.)},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Constant(.001),limits),ResourceLimit);
    CHECK_THROWS_AS(assemble_section_thermal_addition(growth,{material(2000000.,400.)},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Constant(.001),{},
        [] { throw std::runtime_error("cancel exact thermal mass"); }),std::runtime_error);
}

TEST_CASE("Exact retained thermal skin excludes partial height mixed material contacts", "[SectionThermalAssembly]")
{
    const auto left=make_exact_section({rectangle(0,1000000000)});
    const auto topology=qualify_section_topology({{0.,.5,0,SectionFootprint(left.section)},
        {.5,1.,0,SectionFootprint(left.section)},
        {0.,.5,1,{rectangle(1000000000,2000000000)}}},0.);
    const auto full=assemble_section_thermal_volume(topology,{material(2.,300.),material(3.,400.)},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Constant(.5));
    CHECK(full.exposure_complete);
    CHECK_THAT(full.occupied_volume_m3,WithinAbs(1.5,1e-12));
    double horizontal=0.,vertical=0.,positive_x=0.;
    for (const auto& face : full.horizontal_faces) horizontal+=face.area_moments[0];
    for (const auto& face : full.vertical_faces) {
        vertical+=face.area_moments[0];
        if (face.outward_normal.x()>.9) positive_x+=face.area_moments[0];
    }
    CHECK_THAT(horizontal,WithinAbs(4.,1e-12));
    CHECK_THAT(vertical,WithinAbs(5.,1e-12));
    CHECK_THAT(positive_x,WithinAbs(1.,1e-12));
    SectionThermalAssemblyLimits limits;
    limits.max_points=8;
    CHECK_THROWS_AS(assemble_section_thermal_volume(topology,{material(2.,300.),material(3.,400.)},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Constant(.5),limits),ResourceLimit);
    CHECK_THROWS_AS(assemble_section_thermal_volume(topology,{material(2.,300.),material(3.,400.)},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Constant(.5),{},
        [] { throw std::runtime_error("cancel exact skin"); }),std::runtime_error);
}

TEST_CASE("Exact retained thermal skin preserves hole walls on cell faces", "[SectionThermalAssembly]")
{
    auto ring=rectangle(0,1000000000);
    ring.holes.push_back(Polygon{{250000000,250000000},{250000000,750000000},
        {750000000,750000000},{750000000,250000000}});
    const auto exact=make_exact_section({ring});
    const auto topology=qualify_section_topology({{0.,1.,0,SectionFootprint(exact.section)}},0.);
    const auto full=assemble_section_thermal_volume(topology,{material(2.,300.)},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Constant(.25));
    double horizontal=0.,vertical=0.;
    for (const auto& face : full.horizontal_faces) horizontal+=face.area_moments[0];
    for (const auto& face : full.vertical_faces) vertical+=face.area_moments[0];
    CHECK(full.exposure_complete);
    CHECK_THAT(full.occupied_volume_m3,WithinAbs(.75,1e-12));
    CHECK_THAT(horizontal,WithinAbs(1.5,1e-12));
    CHECK_THAT(vertical,WithinAbs(6.,1e-12));
}

TEST_CASE("Disconnected sections never share thermal nodes inside one background cell", "[SectionThermalAssembly]")
{
    const auto topology = qualify_section_topology({{0.,1.,0,{rectangle(0,250000000)}},
        {0.,1.,1,{rectangle(750000000,1000000000)}}},0.);
    REQUIRE(topology.components() == 2);
    const auto assembly = assemble_section_thermal_volume(topology,{material(1.,400.),material(1.,300.)},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones());
    REQUIRE(assembly.nodes.size() == 16);
    const auto history = ThermalHistory::evolve({assembly.initial},10.);
    REQUIRE(history.state == ThermalHistory::State::Converged);
    for (std::size_t i = 0; i < assembly.nodes.size(); ++i)
        CHECK_THAT(history.temperatures_k[i],WithinAbs(assembly.initial.deposited_cells[i].temperature_k,1e-8));
    CHECK_THAT(assembly.occupied_volume_m3,WithinAbs(.5,1e-12));
}

TEST_CASE("Section thermal assembly preserves hole volume when the grid is refined", "[SectionThermalAssembly]")
{
    auto polygon = rectangle(0,1000000000);
    Polygon hole{{250000000,250000000},{750000000,250000000},{750000000,750000000},{250000000,750000000}};
    hole.make_clockwise();
    polygon.holes.push_back(hole);
    const auto topology = qualify_section_topology({{0.,1.,0,{polygon}}},0.);
    for (double width : {1.,.5,.25}) {
        const auto assembly = assemble_section_thermal_volume(topology,{material(2.,300.)},
            Eigen::Vector3d::Zero(),{width,width,1.});
        CHECK_THAT(assembly.occupied_volume_m3,WithinAbs(.75,1e-12));
        double capacity = 0.;
        for (const auto& cell : assembly.initial.deposited_cells)
            capacity += cell.volume_m3*cell.volumetric_heat_capacity_j_m3_k;
        CHECK_THAT(capacity,WithinAbs(1.5,1e-11));
        double exposed_area = 0.;
        for (const auto& face : assembly.horizontal_faces) exposed_area += face.area_moments[0];
        CHECK_THAT(exposed_area,WithinAbs(1.5,1e-12));
        double side_area = 0.;
        Eigen::Vector3d side_normal = Eigen::Vector3d::Zero();
        for (const auto& face : assembly.vertical_faces) {
            side_area += face.area_moments[0];
            side_normal += face.area_moments[0]*face.outward_normal;
        }
        CHECK_THAT(side_area,WithinAbs(6.,1e-11));
        CHECK_THAT(side_normal.norm(),WithinAbs(0.,1e-11));
    }
}

TEST_CASE("Certified section geometry drives a bed to air thermal history without cooling internal faces", "[SectionThermalAssembly]")
{
    const auto topology = qualify_section_topology({{0.,1.,0,{rectangle(0,1000000000)}},
        {1.,2.,1,{rectangle(0,1000000000)}}},0.);
    auto assembly = assemble_section_thermal_volume(topology,{material(1.,350.),material(1.,350.)},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones());
    for (const auto& face : assembly.horizontal_faces) {
        REQUIRE(face.kind != SectionHorizontalFaceKind::Downward);
        const bool bed = face.kind == SectionHorizontalFaceKind::Bed;
        const auto integrated = SectionIntegration::thermal_boundary_element(face.area_moments,1.,bed ? 400. : 300.);
        ThermalHistory::Element boundary;
        boundary.nodes = face.nodes;
        boundary.conductance = integrated.conductance;
        boundary.load = integrated.load;
        boundary.boundary = true;
        boundary.kind = bed ? ThermalHistory::BoundaryKind::Bed : ThermalHistory::BoundaryKind::Exhaust;
        assembly.initial.exposed_elements.push_back(boundary);
    }
    ThermalHistory::HistoryLimits limits;
    limits.solve.relative_tolerance = 1e-12;
    const auto result = ThermalHistory::evolve({assembly.initial},100.,limits);
    REQUIRE(result.state == ThermalHistory::State::Converged);
    for (std::size_t i = 0; i < assembly.nodes.size(); ++i)
        CHECK_THAT(result.temperatures_k[i],WithinAbs(375.-25.*assembly.nodes[i].position_m.z(),1e-6));
    CHECK_THAT(result.energy_balance_error_j,WithinAbs(0.,1e-7));
    const std::vector<SectionExpansionMaterial> expansion{
        {{{250.,450.,Eigen::Matrix3d::Identity()*1e-5}},Eigen::Matrix3d::Identity(),400.},
        {{{250.,450.,Eigen::Matrix3d::Identity()*2e-5}},Eigen::Matrix3d::Identity(),400.}};
    const auto load = section_thermal_load(assembly,result,expansion);
    std::vector<SectionMaterial> elastic(2,SectionMaterial{Eigen::Matrix<double,6,6>::Identity()*1000.});
    const std::vector<SectionSampleRequest> samples{{0,{.3,.4,.4},0},{1,{.3,.4,1.4},1}};
    const auto coupled = solve_deposited_sections(topology,elastic,Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones(),{},{},samples,{},{},load);
    for (std::size_t i = 0; i < elastic.size(); ++i) {
        const double alpha = double(i+1)*1e-5;
        elastic[i].eigenstrain.head<3>().setConstant(-25.*alpha);
        elastic[i].eigenstrain_gradient_per_m.col(2).head<3>().setConstant(-25.*alpha);
    }
    const auto expected = solve_deposited_sections(topology,elastic,Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones(),{},{},samples);
    REQUIRE(coupled.samples.size() == expected.samples.size());
    for (std::size_t i = 0; i < coupled.samples.size(); ++i) {
        CHECK_THAT((coupled.samples[i].stress_pa-expected.samples[i].stress_pa).norm(),WithinAbs(0.,1e-7));
        CHECK_THAT((coupled.samples[i].displacement_m-expected.samples[i].displacement_m).norm(),WithinAbs(0.,1e-10));
    }
}

TEST_CASE("Section exposure retains ledges and unsupported lower faces while excluding positive contacts", "[SectionThermalAssembly]")
{
    // Upper slab overhangs half of a lower slab, with a half-square contact.
    const auto topology = qualify_section_topology({{0.,1.,0,{rectangle(0,1000000000)}},
        {1.,2.,0,{rectangle(500000000,1500000000)}}},0.);
    const auto assembly = assemble_section_thermal_volume(topology,{material(1.,300.)},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones());
    double up = 0., down = 0., bed = 0.;
    for (const auto& face : assembly.horizontal_faces) {
        if (face.kind == SectionHorizontalFaceKind::Upward) up += face.area_moments[0];
        else if (face.kind == SectionHorizontalFaceKind::Downward) down += face.area_moments[0];
        else bed += face.area_moments[0];
    }
    CHECK_THAT(up,WithinAbs(1.5,1e-12));
    CHECK_THAT(down,WithinAbs(.5,1e-12));
    CHECK_THAT(bed,WithinAbs(1.,1e-12));
}

TEST_CASE("Section thermal assembly handles negative grid indices without changing mass", "[SectionThermalAssembly]")
{
    const auto topology = qualify_section_topology({{0.,1.,0,{rectangle(-500000000,500000000)}}},0.);
    const auto assembly = assemble_section_thermal_volume(topology,{material(1.,300.)},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones());
    CHECK(assembly.occupied_cells == 2);
    CHECK_THAT(assembly.occupied_volume_m3,WithinAbs(1.,1e-12));
    const auto result = ThermalHistory::evolve({assembly.initial},1.);
    REQUIRE(result.state == ThermalHistory::State::Converged);
    CHECK_THAT(result.final_energy_j,WithinAbs(300.,1e-9));
}

TEST_CASE("Section thermal assembly rejects missing properties and honors resource limits", "[SectionThermalAssembly]")
{
    const auto topology = qualify_section_topology({{0.,1.,0,{rectangle(0,1000000000)}}},0.);
    const Eigen::Vector3d origin = Eigen::Vector3d::Zero(), spacing = Eigen::Vector3d::Ones();
    SectionThermalAssemblyLimits limits;
    limits.max_nodes = 7;
    CHECK_THROWS_AS(assemble_section_thermal_volume(topology,{material(1.,300.)},origin,spacing,limits),ResourceLimit);
    limits = {};
    limits.max_work = 1;
    CHECK_THROWS_AS(assemble_section_thermal_volume(topology,{material(1.,300.)},origin,spacing,limits),ResourceLimit);
    CHECK_THROWS_AS(assemble_section_thermal_volume(topology,{},origin,spacing),std::invalid_argument);
    CHECK_THROWS_AS(assemble_section_thermal_volume(topology,{material(0.,300.)},origin,spacing),std::invalid_argument);
    CHECK_THROWS_AS(assemble_section_thermal_volume(topology,{material(1.,std::numeric_limits<double>::quiet_NaN())},origin,spacing),std::invalid_argument);
    CHECK_THROWS_AS(assemble_section_thermal_volume(topology,{material(1.,300.)},origin,spacing,{},
        [] { throw std::runtime_error("cancelled"); }),std::runtime_error);
}

TEST_CASE("Section thermal refinement resolves serial material resistance instead of trusting a coarse average", "[SectionThermalAssembly]")
{
    const auto topology = qualify_section_topology({{0.,.5,0,{rectangle(0,1000000000)}},
        {.5,1.,1,{rectangle(0,1000000000)}}},0.);
    auto high_conductivity = material(1.,350.);
    high_conductivity.conductivity_w_m_k *= 10.;
    const auto bed_flux = [&](double spacing_z) {
        auto assembly = assemble_section_thermal_volume(topology,{material(1.,350.),high_conductivity},
            Eigen::Vector3d::Zero(),{1.,1.,spacing_z});
        for (const auto& face : assembly.horizontal_faces) {
            const bool bed = face.kind == SectionHorizontalFaceKind::Bed;
            const auto integrated = SectionIntegration::thermal_boundary_element(face.area_moments,100.,bed ? 400. : 300.);
            ThermalHistory::Element boundary;
            boundary.nodes = face.nodes;
            boundary.conductance = integrated.conductance;
            boundary.load = integrated.load;
            boundary.boundary = true;
            boundary.kind = bed ? ThermalHistory::BoundaryKind::Bed : ThermalHistory::BoundaryKind::Ambient;
            assembly.initial.exposed_elements.push_back(boundary);
        }
        ThermalHistory::HistoryLimits limits;
        limits.solve.relative_tolerance = 1e-12;
        const auto result = ThermalHistory::evolve({assembly.initial},100.,limits);
        REQUIRE(result.state == ThermalHistory::State::Converged);
        double flux = 0.;
        for (const auto& face : assembly.initial.exposed_elements) {
            if (face.kind != ThermalHistory::BoundaryKind::Bed) continue;
            Eigen::Matrix<double,8,1> temperature;
            for (unsigned i = 0; i < 8; ++i) temperature[i] = result.temperatures_k[face.nodes[i]];
            flux += (face.load-face.conductance*temperature).sum();
        }
        return flux; // One square metre of bed contact in this analytical case.
    };
    const double exact = 100./(.01+.5/1.+.5/10.+.01);
    CHECK(std::abs(bed_flux(1.)-exact) > .1*exact);
    CHECK_THAT(bed_flux(.5),WithinAbs(exact,1e-5));
    CHECK_THAT(bed_flux(.25),WithinAbs(exact,1e-5));
}

TEST_CASE("Vertical section exposure excludes internal material interfaces on and between grid faces", "[SectionThermalAssembly]")
{
    const auto topology = qualify_section_topology({{0.,1.,0,{rectangle(0,500000000)}},
        {0.,1.,1,{rectangle(500000000,1000000000)}}},0.);
    for (double spacing : {1.,.5,.25}) {
        const auto assembly = assemble_section_thermal_volume(topology,{material(1.,300.),material(1.,300.)},
            Eigen::Vector3d::Zero(),{spacing,spacing,1.});
        double area = 0.;
        for (const auto& face : assembly.vertical_faces) area += face.area_moments[0];
        CHECK_THAT(area,WithinAbs(4.,1e-11));
    }
}

TEST_CASE("Vertical section exposure subtracts only the shared height of a neighboring prism", "[SectionThermalAssembly]")
{
    const std::vector<LocalSectionTile> tiles{
        {{0,0},{0,0},{1000000000,1000000000},{{0.,2.,0,{rectangle(0,1000000000)}}}},
        {{1,0},{1000000000,0},{2000000000,1000000000},{{1.,2.,1,{rectangle(1000000000,2000000000)}}}}};
    const auto topology = qualify_local_section_topology(tiles,0.);
    const auto assembly = assemble_section_thermal_volume(topology,{material(1.,300.),material(1.,300.)},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones());
    double area = 0.;
    for (const auto& face : assembly.vertical_faces) area += face.area_moments[0];
    CHECK_THAT(area,WithinAbs(10.,1e-11));
}

TEST_CASE("Oblique section side area does not gain artificial grid-cut faces", "[SectionThermalAssembly]")
{
    const ExPolygon triangle(Polygon{{0,0},{1000000000,0},{0,700000001}});
    const auto topology = qualify_section_topology({{0.,1.,0,{triangle}}},0.);
    const double height = .700000001;
    const double exact = 1.+height+std::hypot(1.,height);
    for (double spacing : {1.,.3,.2}) {
        const auto assembly = assemble_section_thermal_volume(topology,{material(1.,300.)},
            Eigen::Vector3d::Zero(),{spacing,spacing,1.});
        double area = 0.;
        Eigen::Vector3d normal = Eigen::Vector3d::Zero();
        for (const auto& face : assembly.vertical_faces) {
            area += face.area_moments[0];
            normal += face.area_moments[0]*face.outward_normal;
        }
        CHECK_THAT(area,WithinAbs(exact,1e-9));
        CHECK_THAT(normal.norm(),WithinAbs(0.,1e-9));
    }
}

TEST_CASE("Geometry-derived side cooling participates in the global heat balance", "[SectionThermalAssembly]")
{
    const bool exact = GENERATE(false,true);
    const ExPolygons polygon{rectangle(0,1000000000)};
    const SectionFootprint footprint = exact ? SectionFootprint(make_exact_section(polygon).section) : SectionFootprint(polygon);
    const auto topology = qualify_section_topology({{0.,1.,0,footprint}},0.);
    auto assembly = assemble_section_thermal_volume(topology,{material(1.,400.)},
        Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones());
    for (const auto& face : assembly.vertical_faces) {
        const auto integrated = SectionIntegration::thermal_boundary_element(face.area_moments,1.,300.);
        ThermalHistory::Element boundary;
        boundary.nodes = face.nodes;
        boundary.conductance = integrated.conductance;
        boundary.load = integrated.load;
        boundary.boundary = true;
        boundary.kind = ThermalHistory::BoundaryKind::PartFan;
        assembly.initial.exposed_elements.push_back(boundary);
    }
    const auto result = ThermalHistory::evolve({assembly.initial},1.);
    REQUIRE(result.state == ThermalHistory::State::Converged);
    for (double temperature : result.temperatures_k) CHECK_THAT(temperature,WithinAbs(320.,1e-9));
    CHECK_THAT(result.boundary_energy_j,WithinAbs(-80.,1e-9));
    CHECK_THAT(result.energy_balance_error_j,WithinAbs(0.,1e-9));

    const std::vector<SectionExpansionMaterial> expansion{{
        {{250.,450.,Eigen::Matrix3d::Identity()*1e-5}},Eigen::Matrix3d::Identity(),400.}};
    const auto load = section_thermal_load(assembly,result,expansion);
    SectionMaterial elastic{Eigen::Matrix<double,6,6>::Identity()*1000.};
    const std::vector<SectionSampleRequest> samples{{0,{.3,.4,.6},0}};
    const auto coupled = solve_deposited_sections(topology,{elastic},Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones(),{},{},samples,{},{},load);
    elastic.eigenstrain.head<3>().setConstant(-80e-5);
    const auto expected = solve_deposited_sections(topology,{elastic},Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones(),{},{},samples);
    REQUIRE(coupled.samples.size() == 1);
    REQUIRE(expected.samples.size() == 1);
    CHECK(coupled.elastic_energy_j > 0.);
    CHECK_THAT(coupled.elastic_energy_j,WithinAbs(expected.elastic_energy_j,1e-12));
    CHECK_THAT((coupled.samples[0].stress_pa-expected.samples[0].stress_pa).norm(),WithinAbs(0.,1e-9));
    CHECK_THAT((coupled.samples[0].displacement_m-expected.samples[0].displacement_m).norm(),WithinAbs(0.,1e-12));

    auto failed = result;
    failed.state = ThermalHistory::State::Unconverged;
    CHECK_THROWS_AS(section_thermal_load(assembly,failed,expansion),std::invalid_argument);
    failed = result;
    failed.temperatures_k.pop_back();
    CHECK_THROWS_AS(section_thermal_load(assembly,failed,expansion),std::invalid_argument);
    const auto& group = assembly.groups.front();
    SectionLoadQuery query{group.material,group.component,group.cell,group.minimum_m,group.maximum_m,group.occupied_volume_m3};
    CHECK_NOTHROW(load(query));
    query.maximum_m.z() += .01;
    CHECK_THROWS_AS(load(query),std::invalid_argument);
    query.maximum_m = group.maximum_m;
    query.material = 1;
    CHECK_THROWS_AS(load(query),std::out_of_range);
}
