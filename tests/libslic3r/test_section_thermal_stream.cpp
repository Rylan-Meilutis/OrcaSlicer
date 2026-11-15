#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/SectionThermalAssembly.hpp"
#include "libslic3r/SectionIntegration.hpp"
#include "libslic3r/ThermalDeposition.hpp"
#include "libslic3r/ThermalHistory.hpp"
#include "libslic3r/ThermoelasticLocalSections.hpp"
#include "libslic3r/ThermoelasticLocalTopology.hpp"
#include "libslic3r/ThermoelasticSectionTopology.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Offline thermal work allowance does not relax live assembly limits", "[SectionThermalStream]")
{
    const std::vector<DepositedSegment> beads{
        {{0.,0.,.001},{.002,0.,.001},.001,.001,0,BeadSection::RectangularEnvelope}};
    const std::vector<SectionThermalMaterial> materials{{Eigen::Matrix3d::Identity(),1e6,300.}};
    ThermalSectionStreamSettings settings;
    settings.origin_m = {-.001,-.001,0.};
    settings.spacing_m = {.001,.001,.001};
    settings.tile_size_m = {.001,.001};
    settings.maximum_height_m = .001;
    settings.max_work = 2000000000;
    CHECK_NOTHROW(assemble_local_deposition_thermal_volume(beads,materials,settings));
    ++settings.max_work;
    CHECK_THROWS_AS(assemble_local_deposition_thermal_volume(beads,materials,settings),std::invalid_argument);
    --settings.max_work;
    settings.assembly.max_nodes = 1;
    CHECK_THROWS_AS(assemble_local_deposition_thermal_volume(beads,materials,settings),ResourceLimit);
}

TEST_CASE("Repeated edge integration preserves layered material ownership and grid clipping", "[SectionThermalStream]")
{
    std::vector<DepositedSegment> beads;
    for (unsigned layer=0;layer<12;++layer)
        beads.push_back({{0.,.0005,(layer+1)*.0002},{.002,.0005,(layer+1)*.0002},
            .001,.0002,layer%2,BeadSection::RectangularEnvelope});
    const std::vector<SectionThermalMaterial> materials{{Eigen::Matrix3d::Identity(),1e6,300.},
        {Eigen::Matrix3d::Identity()*.2,2e6,400.}};
    ThermalSectionStreamSettings settings;
    settings.spacing_m={.001,.001,.0004};
    settings.tile_size_m={.002,.002};
    settings.maximum_height_m=.0001;
    settings.include_exposure=true;
    const auto geometry=local_deposition_sections(beads,settings.origin_m.head<2>(),settings.tile_size_m,
        0.,settings.maximum_height_m);
    const auto topology=qualify_local_section_topology(geometry.tiles,0.);
    const auto reference=assemble_section_thermal_volume(topology,materials,settings.origin_m,settings.spacing_m);
    const auto actual=assemble_local_deposition_thermal_volume(beads,materials,settings).assembly;
    using Key=std::pair<std::size_t,std::array<double,3>>;
    const auto moments=[](const SectionThermalAssembly& assembly) {
        std::map<Key,Slic3r::SectionIntegration::Moments> sums;
        for (const auto& face : assembly.vertical_faces) {
            const Key key{face.material,{face.outward_normal.x(),face.outward_normal.y(),face.outward_normal.z()}};
            for (std::size_t i=0;i<face.area_moments.size();++i) sums[key][i]+=face.area_moments[i];
        }
        return sums;
    };
    const auto expected=moments(reference), observed=moments(actual);
    REQUIRE(observed.size()==expected.size());
    for (const auto& [key,values] : expected) {
        REQUIRE(observed.count(key)==1);
        for (std::size_t i=0;i<values.size();++i)
            CHECK_THAT(observed.at(key)[i],WithinAbs(values[i],1e-16));
    }
    CHECK_THAT(actual.occupied_volume_m3,WithinRel(reference.occupied_volume_m3,1e-12));
}

TEST_CASE("Streamed rounded bead thermal volume matches retained cut cells and material heat", "[SectionThermalStream]")
{
    const std::vector<DepositedSegment> beads{
        {{0.,0.,.0005},{.002,0.,.0005},.001,.0005,0,BeadSection::RoundedRectangle},
        {{0.,0.,.001},{.002,0.,.001},.001,.0005,1,BeadSection::RoundedRectangle}};
    const std::vector<SectionThermalMaterial> materials{{Eigen::Matrix3d::Identity()*.2,2e6,400.},
        {Eigen::Matrix3d::Identity()*.1,1e6,300.}};
    ThermalSectionStreamSettings settings;
    settings.origin_m = {-.001,-.001,0.};
    settings.spacing_m = {.001,.001,.0005};
    settings.tile_size_m = {.001,.001};
    settings.maximum_height_m = .0001;
    const auto geometry = local_deposition_sections(beads,settings.origin_m.head<2>(),settings.tile_size_m,
        0.,settings.maximum_height_m);
    const auto topology = qualify_local_section_topology(geometry.tiles,0.);
    const auto reference = assemble_section_thermal_volume(topology,materials,settings.origin_m,settings.spacing_m);
    const auto streamed = assemble_local_deposition_thermal_volume(beads,materials,settings);
    CHECK_FALSE(streamed.assembly.exposure_complete);
    CHECK(streamed.assembly.horizontal_faces.empty());
    CHECK(streamed.assembly.vertical_faces.empty());
    CHECK(reference.exposure_complete);
    CHECK_THAT(streamed.assembly.occupied_volume_m3,WithinRel(reference.occupied_volume_m3,1e-10));
    REQUIRE(streamed.assembly.nodes.size() == reference.nodes.size());
    REQUIRE(streamed.assembly.initial.added_elements.size() == reference.initial.added_elements.size());
    for (std::size_t i = 0; i < reference.nodes.size(); ++i) {
        CHECK_THAT((streamed.assembly.nodes[i].position_m-reference.nodes[i].position_m).norm(),WithinAbs(0.,1e-15));
        const auto& a = streamed.assembly.initial.deposited_cells[i];
        const auto& b = reference.initial.deposited_cells[i];
        CHECK_THAT(a.volume_m3,WithinRel(b.volume_m3,1e-10));
        CHECK_THAT(a.temperature_k,WithinAbs(b.temperature_k,1e-8));
        CHECK_THAT(a.volumetric_heat_capacity_j_m3_k,WithinRel(b.volumetric_heat_capacity_j_m3_k,1e-10));
    }
    for (std::size_t i = 0; i < reference.initial.added_elements.size(); ++i)
        CHECK_THAT((streamed.assembly.initial.added_elements[i].conductance-reference.initial.added_elements[i].conductance).norm(),
            WithinAbs(0.,1e-12));
    // Explicitly insulated numerical comparison; missing skin is NOT automatic
    // insulation for a whole-part cooling scenario.
    const auto a = Slic3r::ThermalHistory::evolve({streamed.assembly.initial},1.);
    const auto b = Slic3r::ThermalHistory::evolve({reference.initial},1.);
    REQUIRE(a.state == Slic3r::ThermalHistory::State::Converged);
    REQUIRE(b.state == Slic3r::ThermalHistory::State::Converged);
    CHECK_THAT(a.final_energy_j,WithinAbs(b.final_energy_j,1e-8));

    settings.include_exposure = true;
    auto exposed = assemble_local_deposition_thermal_volume(beads,materials,settings).assembly;
    REQUIRE(exposed.exposure_complete);
    // Separable XY/Z integration must preserve ownership of matching XY faces
    // at different heights. Compare every normalized moment against retained
    // geometry, separately for each material.
    std::vector<Slic3r::SectionIntegration::Moments> streamed_moments(materials.size()), retained_moments(materials.size());
    for (const auto& face : exposed.vertical_faces)
        for (std::size_t i = 0; i < face.area_moments.size(); ++i)
            streamed_moments.at(face.material)[i] += face.area_moments[i];
    for (const auto& face : reference.vertical_faces)
        for (std::size_t i = 0; i < face.area_moments.size(); ++i)
            retained_moments.at(face.material)[i] += face.area_moments[i];
    for (std::size_t m = 0; m < materials.size(); ++m)
        for (std::size_t i = 0; i < streamed_moments[m].size(); ++i)
            CHECK_THAT(streamed_moments[m][i],WithinAbs(retained_moments[m][i],1e-16));
    auto retained = reference;
    const auto cool = [](SectionThermalAssembly& assembly) {
        const auto add = [&](const auto& face,bool bed) {
            const auto integrated = Slic3r::SectionIntegration::thermal_boundary_element(face.area_moments,bed ? 100. : 20.,bed ? 380. : 290.);
            Slic3r::ThermalHistory::Element element;
            element.nodes = face.nodes; element.conductance = integrated.conductance;
            element.load = integrated.load; element.boundary = true;
            element.kind = bed ? Slic3r::ThermalHistory::BoundaryKind::Bed : Slic3r::ThermalHistory::BoundaryKind::Ambient;
            assembly.initial.exposed_elements.push_back(element);
        };
        for (const auto& face : assembly.horizontal_faces) add(face,face.kind == SectionHorizontalFaceKind::Bed);
        for (const auto& face : assembly.vertical_faces) add(face,false);
        return Slic3r::ThermalHistory::evolve({assembly.initial},.1);
    };
    const auto streamed_cooling = cool(exposed), retained_cooling = cool(retained);
    REQUIRE(streamed_cooling.state == Slic3r::ThermalHistory::State::Converged);
    REQUIRE(retained_cooling.state == Slic3r::ThermalHistory::State::Converged);
    CHECK_THAT(streamed_cooling.final_energy_j,WithinAbs(retained_cooling.final_energy_j,1e-9));
    REQUIRE(streamed_cooling.temperatures_k.size() == retained_cooling.temperatures_k.size());
    for (std::size_t i=0;i<streamed_cooling.temperatures_k.size();++i)
        CHECK_THAT(streamed_cooling.temperatures_k[i],WithinAbs(retained_cooling.temperatures_k[i],1e-6));
    // Angular quadrature changes only the sampled direction of exposure, not
    // occupied material, surface moments, or isotropic heat exchange.
    for (const unsigned bins : {8u,32u}) {
        settings.vertical_normal_bins=bins;
        auto angular=assemble_local_deposition_thermal_volume(beads,materials,settings).assembly;
        CHECK(angular.vertical_normal_bins==bins);
        CHECK(angular.vertical_faces.size()<exposed.vertical_faces.size());
        Slic3r::SectionIntegration::Moments exact_area{}, angular_area{};
        for (const auto& face : exposed.vertical_faces)
            for (std::size_t i=0;i<exact_area.size();++i) exact_area[i]+=face.area_moments[i];
        for (const auto& face : angular.vertical_faces) {
            CHECK_THAT(face.outward_normal.norm(),WithinAbs(1.,1e-12));
            for (std::size_t i=0;i<angular_area.size();++i) angular_area[i]+=face.area_moments[i];
        }
        for (std::size_t i=0;i<exact_area.size();++i)
            CHECK_THAT(angular_area[i],WithinAbs(exact_area[i],1e-16));
        const auto directional_area=[](const SectionThermalAssembly& assembly) {
            const Eigen::Vector3d wind=Eigen::Vector3d{1.,2.,0.}.normalized();
            double area=0.;
            for (const auto& face : assembly.vertical_faces)
                area+=face.area_moments[0]*std::max(0.,face.outward_normal.dot(wind));
            return area;
        };
        // A positive projected airflow law is Lipschitz in the normal. The
        // nearest-bin angular error is at most pi/N; no fitted tolerance.
        CHECK(std::abs(directional_area(angular)-directional_area(exposed))<=
              exact_area[0]*2.*std::sin(std::acos(-1.)/(2.*bins))+1e-16);
        const auto angular_cooling=cool(angular);
        REQUIRE(angular_cooling.state==Slic3r::ThermalHistory::State::Converged);
        REQUIRE(angular_cooling.temperatures_k.size()==streamed_cooling.temperatures_k.size());
        for (std::size_t i=0;i<angular_cooling.temperatures_k.size();++i)
            CHECK_THAT(angular_cooling.temperatures_k[i],WithinAbs(streamed_cooling.temperatures_k[i],1e-6));
    }
    settings.vertical_normal_bins=7;
    CHECK_THROWS_AS(assemble_local_deposition_thermal_volume(beads,materials,settings),std::invalid_argument);
}

TEST_CASE("Bounded surface batches preserve every moment across repeated flushes", "[SectionThermalStream]")
{
    const std::vector<DepositedSegment> beads{
        {{.002,.005,.001},{.008,.005,.001},.002,.001,0,BeadSection::RoundedRectangle}};
    const std::vector<SectionThermalMaterial> materials{{Eigen::Matrix3d::Identity()*.2,2e6,400.}};
    ThermalSectionStreamSettings settings;
    settings.origin_m = {0.,0.,0.};
    settings.spacing_m = {.01,.01,.002};
    settings.tile_size_m = {.01,.01};
    settings.maximum_height_m = .000001;
    settings.include_exposure = true;
    settings.max_work = 50000000;
    settings.skin.max_work = 50000000;
    const auto exact = assemble_local_deposition_thermal_volume(beads,materials,settings).assembly;
    // More than4096 distinct band/normal keys must be flushed, not retained
    // until the full stream completes or silently dropped at capacity.
    for (const unsigned bins : {8u,32u}) {
        settings.vertical_normal_bins = bins;
        const auto batched = assemble_local_deposition_thermal_volume(beads,materials,settings).assembly;
        Slic3r::SectionIntegration::Moments expected{}, actual{};
        for (const auto& face : exact.vertical_faces)
            for (std::size_t i=0;i<expected.size();++i) expected[i]+=face.area_moments[i];
        for (const auto& face : batched.vertical_faces)
            for (std::size_t i=0;i<actual.size();++i) actual[i]+=face.area_moments[i];
        for (std::size_t i=0;i<actual.size();++i)
            CHECK_THAT(actual[i],WithinAbs(expected[i],1e-16));
        CHECK_THAT(batched.occupied_volume_m3,WithinRel(exact.occupied_volume_m3,1e-12));
    }
    settings.max_provisional_groups = 3;
    try {
        (void)assemble_local_deposition_thermal_volume(beads,materials,settings);
        FAIL("Flushed batches must respect the side-face group limit.");
    } catch (const ResourceLimit& error) {
        CHECK(std::string(error.what()).find("side-face budget exceeded") != std::string::npos);
    }
}

TEST_CASE("Surface batches retain separate material and component moments in one cell", "[SectionThermalStream]")
{
    const std::vector<DepositedSegment> beads{
        {{.002,.003,.001},{.008,.003,.001},.002,.001,0,BeadSection::RoundedRectangle},
        {{.002,.007,.001},{.008,.007,.001},.002,.001,1,BeadSection::RoundedRectangle}};
    const std::vector<SectionThermalMaterial> materials{
        {Eigen::Matrix3d::Identity()*.2,2e6,400.},{Eigen::Matrix3d::Identity()*.1,1e6,300.}};
    ThermalSectionStreamSettings settings;
    settings.origin_m = {0.,0.,0.};
    settings.spacing_m = {.01,.01,.002};
    settings.tile_size_m = {.01,.01};
    settings.maximum_height_m = .0001;
    settings.include_exposure = true;
    const auto exact = assemble_local_deposition_thermal_volume(beads,materials,settings).assembly;
    settings.vertical_normal_bins = 32;
    const auto batched = assemble_local_deposition_thermal_volume(beads,materials,settings).assembly;
    using Key = std::pair<std::array<std::size_t,8>,std::size_t>;
    const auto collect = [](const SectionThermalAssembly& field) {
        std::map<Key,Slic3r::SectionIntegration::Moments> totals;
        for (const auto& face : field.vertical_faces) {
            auto& target = totals[{face.nodes,face.material}];
            for (std::size_t i=0;i<target.size();++i) target[i]+=face.area_moments[i];
        }
        return totals;
    };
    const auto expected = collect(exact), actual = collect(batched);
    REQUIRE(expected.size() == 2);
    REQUIRE(actual.size() == expected.size());
    for (const auto& [key,moments] : expected) {
        REQUIRE(actual.count(key) == 1);
        for (std::size_t i=0;i<moments.size();++i)
            CHECK_THAT(actual.at(key)[i],WithinAbs(moments[i],1e-16));
    }
}

TEST_CASE("Streamed scarf thermal volume applies physical cell caps after component completion", "[SectionThermalStream]")
{
    const std::vector<DepositedSegment> beads{
        {{.001,.001,.0005},{.003,.001,.0008},.001,.0005,0,BeadSection::RoundedRectangle,.0008}};
    ThermalSectionStreamSettings settings;
    settings.spacing_m = {.01,.01,.01};
    settings.tile_size_m = {.001,.001};
    settings.maximum_height_m = .0001;
    settings.assembly.max_cells = 1;
    settings.assembly.max_groups = 1;
    const std::vector<SectionThermalMaterial> materials{{Eigen::Matrix3d::Identity(),2e6,400.}};
    const auto result = assemble_local_deposition_thermal_volume(beads,materials,settings);
    CHECK(result.assembly.occupied_cells == 1);
    CHECK_THAT(result.assembly.occupied_volume_m3,WithinRel(result.geometry.volume_m3,1e-10));
    CHECK(result.geometry.sections > 1);
    auto restricted = settings;
    restricted.max_processed_points = 1;
    CHECK_THROWS_AS(assemble_local_deposition_thermal_volume(beads,materials,restricted),ResourceLimit);
    restricted = settings;
    restricted.max_work = 1;
    CHECK_THROWS_AS(assemble_local_deposition_thermal_volume(beads,materials,restricted),ResourceLimit);
    CHECK_THROWS_AS(assemble_local_deposition_thermal_volume(beads,materials,settings,{},
        [] { throw std::runtime_error("cancel stream"); }),std::runtime_error);
}

TEST_CASE("Streamed thermal cells do not conduct across an air gap inside one background cell", "[SectionThermalStream]")
{
    const std::vector<DepositedSegment> beads{
        {{.001,.001,.0005},{.003,.001,.0005},.0005,.0005,0,BeadSection::RoundedRectangle},
        {{.001,.003,.0005},{.003,.003,.0005},.0005,.0005,1,BeadSection::RoundedRectangle}};
    const std::vector<SectionThermalMaterial> materials{
        {Eigen::Matrix3d::Identity()*.2,2e6,400.},
        {Eigen::Matrix3d::Identity()*.1,1e6,300.}};
    ThermalSectionStreamSettings settings;
    settings.spacing_m = {.01,.01,.01};
    settings.tile_size_m = {.001,.001};
    settings.maximum_height_m = .0001;
    const auto result = assemble_local_deposition_thermal_volume(beads,materials,settings);
    REQUIRE(result.assembly.groups.size() == 2);
    const auto& first = result.assembly.groups[0];
    const auto& second = result.assembly.groups[1];
    CHECK(first.cell == second.cell);
    CHECK(first.component != second.component);
    for (const auto a : first.nodes)
        for (const auto b : second.nodes)
            CHECK(a != b);
    const auto history = Slic3r::ThermalHistory::evolve({result.assembly.initial},10.);
    REQUIRE(history.state == Slic3r::ThermalHistory::State::Converged);
    for (const auto& group : result.assembly.groups)
        for (const auto node : group.nodes)
            CHECK_THAT(history.temperatures_k[node],WithinAbs(materials[group.material].initial_temperature_k,1e-7));
    // Physical component limits must still reject disconnected geometry after
    // provisional tile labels have been collapsed.
    settings.assembly.max_groups = 1;
    CHECK_THROWS_AS(assemble_local_deposition_thermal_volume(beads,materials,settings),ResourceLimit);
}

TEST_CASE("Thermal streaming partition envelopes never enlarge deposited material", "[SectionThermalStream]")
{
    const std::vector<DepositedSegment> deposited{
        {{.001,.001,.0005},{.003,.001,.0005},.001,.0005,0,BeadSection::RectangularEnvelope}};
    auto partition = deposited;
    partition.front().start_m.z() = partition.front().end_m.z() = .001;
    partition.front().height_m = .001;
    ThermalSectionStreamSettings settings;
    settings.spacing_m = {.001,.001,.0005};
    settings.tile_size_m = {.001,.001};
    settings.maximum_height_m = .0001;
    // The enclosing prism supplies sampling boundaries, not additional
    // volume, capacity or heat above the original bead.
    const std::vector<SectionThermalMaterial> materials{{Eigen::Matrix3d::Identity()*.2,2e6,400.}};
    const auto reference = assemble_local_deposition_thermal_volume(deposited,materials,settings);
    const auto result = assemble_local_deposition_thermal_volume(deposited,materials,settings,{}, {},&partition);
    CHECK_THAT(result.assembly.occupied_volume_m3,WithinRel(reference.assembly.occupied_volume_m3,1e-10));
    double energy = 0., reference_energy = 0.;
    for (const auto& cell : result.assembly.initial.deposited_cells)
        energy += cell.volume_m3*cell.volumetric_heat_capacity_j_m3_k*cell.temperature_k;
    for (const auto& cell : reference.assembly.initial.deposited_cells)
        reference_energy += cell.volume_m3*cell.volumetric_heat_capacity_j_m3_k*cell.temperature_k;
    CHECK_THAT(energy,WithinRel(reference_energy,1e-10));
    for (const auto& group : result.assembly.groups) CHECK(group.material == 0);
}
