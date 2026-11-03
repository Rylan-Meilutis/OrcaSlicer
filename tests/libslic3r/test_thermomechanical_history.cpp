#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ThermalHistory.hpp"
#include "libslic3r/ThermalDeposition.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticMaterial.hpp"
#include "libslic3r/ThermoelasticSolver.hpp"
#include "libslic3r/ThermoelasticThermalStrain.hpp"
#include "libslic3r/ThermoelasticSections.hpp"
#include "libslic3r/ThermoelasticLocalSections.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Polygon.hpp"
#include <Eigen/Core>
#include <cmath>
#include <cstddef>
#include <vector>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
namespace TH = Slic3r::ThermalHistory;
namespace TE = Slic3r::Thermoelastic;

namespace {
// Synthetic SI benchmark, not calibrated ASA, airflow or adhesion parameters.
constexpr double length = .008, width = .002, height = .002;
constexpr double conductivity = .2, alpha = 1e-5, reference_k = 400.;

std::vector<double> slab_temperatures(int layers, double bed_k, double ambient_k, double h)
{
    const double dz = height / layers, area = length * width;
    std::vector<TH::Cell> cells(layers, {area * dz, 1e6, bed_k, 0.});
    std::vector<TH::Connection> links;
    for (int z = 1; z < layers; ++z)
        links.push_back({std::size_t(z - 1), std::size_t(z), conductivity * area / dz});
    const std::vector<TH::Boundary> boundaries{
        {0, 2. * conductivity * area / dz, bed_k, TH::BoundaryKind::Bed},
        {std::size_t(layers - 1), area / (dz / (2. * conductivity) + 1. / h),
            ambient_k, TH::BoundaryKind::PartFan}};
    // A long implicit step approaches the steady resistance-network solution
    // from a uniform initial field. The boundary includes the top half-cell.
    TH::Limits limits;
    limits.relative_tolerance = 1e-13;
    const auto result = TH::step(cells, links, boundaries, 1e10, limits);
    REQUIRE(result.state == TH::State::Converged);
    const double flux = (bed_k - ambient_k) / (height / conductivity + 1. / h);
    for (int z = 0; z < layers; ++z)
        CHECK_THAT(result.temperatures_k[z],
            WithinAbs(bed_k - flux * (z + .5) * dz / conductivity, 1e-6));
    return result.temperatures_k;
}

std::size_t node_at(const TE::Mesh &mesh, const Eigen::Vector3d &position)
{
    for (std::size_t n = 0; n < mesh.nodes_m.size(); ++n)
        if ((mesh.nodes_m[n] - position).norm() < 1e-12) return n;
    FAIL("Benchmark node is missing");
    return 0;
}

TE::Mesh coupon(int resolution, const std::vector<double> &temperatures)
{
    std::vector<TE::Voxel> voxels;
    for (int x = 0; x < 4 * resolution; ++x)
        for (int y = 0; y < resolution; ++y)
            for (int z = 0; z < resolution; ++z) {
                TE::Voxel voxel{{x, y, z}, 1e6, 0.};
                voxel.eigenstrain.head<3>().setConstant(alpha * (temperatures[z] - reference_k));
                voxels.push_back(voxel);
            }
    return TE::mesh_voxels(voxels, Eigen::Vector3d::Zero(),
        {length / (4 * resolution), width / resolution, height / resolution}, 10000);
}

TE::Result free_equilibrium(const TE::Mesh &mesh)
{
    // Remove only six rigid-body modes: these are not bed adhesion constraints.
    const auto origin = node_at(mesh, {0., 0., 0.});
    const auto end = node_at(mesh, {length, 0., 0.});
    const auto side = node_at(mesh, {0., width, 0.});
    auto result = TE::solve(mesh.nodes_m, mesh.elements,
        {{origin, 0}, {origin, 1}, {origin, 2}, {end, 1}, {end, 2}, {side, 2}});
    CHECK(result.relative_residual < 1e-8);
    Eigen::Vector3d net_force = Eigen::Vector3d::Zero();
    for (const auto &force : result.reaction_n) net_force += force;
    CHECK_THAT(net_force.norm(), WithinAbs(0., 1e-9));
    return result;
}
}

TEST_CASE("Uniform thermal cooling maps to stress-free coupon contraction", "[ThermomechanicalHistory]")
{
    const auto temperatures = slab_temperatures(2, 300., 300., 20.);
    const auto mesh = coupon(2, temperatures);
    const auto result = free_equilibrium(mesh);
    for (std::size_t n = 0; n < mesh.nodes_m.size(); ++n)
        CHECK_THAT((result.displacement_m[n] - alpha * (300. - reference_k) * mesh.nodes_m[n]).norm(),
            WithinAbs(0., 1e-12));
    for (const auto &stress : result.stress_pa)
        CHECK_THAT(stress.norm(), WithinAbs(0., 1e-5));
    CHECK_THAT(result.elastic_energy_j, WithinAbs(0., 1e-15));
}

TEST_CASE("Transient temperatures load individual section cells without replacing geometry", "[ThermomechanicalHistory]")
{
    TE::LocalSectionTile tile{{0,0},{0,0},{1000000,1000000},
        {{0.,.001,0,{Slic3r::ExPolygon(Slic3r::Polygon{{0,0},{1000000,0},{1000000,1000000},{0,1000000}})}}}};
    std::vector<TH::Cell> cells(2,{5e-10,1e6,350.,0.});
    const std::vector<TH::Connection> links{{0,1,4e-4}};
    const std::vector<TH::Boundary> boundaries{{1,1e-4,300.,TH::BoundaryKind::PartFan}};
    Eigen::Matrix3d expansion = Eigen::Matrix3d::Zero();
    expansion(2,2) = 1e-5;
    const std::vector<TE::ExpansionInterval> law{{250.,400.,expansion}};
    TE::SectionMaterial elastic{TE::ElasticityMatrix::Identity()*1e6};
    const std::vector<TE::SectionSampleRequest> samples{{0,{.0003,.0003,.00025},0},{1,{.0003,.0003,.00075},0}};
    TH::HistoryObservations observations;
    observations.times_s = {.25,.5,.75};
    unsigned calls = 0;
    observations.observe = [&](double, const std::vector<TH::Cell>& active) {
        ++calls;
        REQUIRE(active.size() == 2);
        CHECK(active[1].temperature_k < active[0].temperature_k);
        const TE::SectionLoadField load = [&](const TE::SectionLoadQuery& q) {
            REQUIRE(q.cell[2] >= 0);
            REQUIRE(q.cell[2] < 2);
            return TE::linear_thermal_eigenstrain(law,Eigen::Matrix3d::Identity(),360.,
                q.minimum_m,active.at(q.cell[2]).temperature_k,Eigen::Vector3d::Zero(),q.minimum_m,q.maximum_m);
        };
        const auto field = TE::solve_local_section_tiles({tile},0.,{elastic},Eigen::Vector3d::Zero(),
            {.0005,.0005,.0005},{},{},{},samples,{},{},{},load);
        REQUIRE(field.samples.size() == 2);
        const double lower = 1e-5*(active[0].temperature_k-360.);
        const double upper = 1e-5*(active[1].temperature_k-360.);
        CHECK_THAT(field.samples[0].displacement_m.z(),WithinAbs(.00025*lower,1e-14));
        CHECK_THAT(field.samples[1].displacement_m.z(),WithinAbs(.0005*lower+.00025*upper,1e-14));
        for (const auto& sample : field.samples) CHECK_THAT(sample.stress_pa.norm(),WithinAbs(0.,1e-6));
        CHECK_THAT(field.elastic_energy_j,WithinAbs(0.,1e-15));
        CHECK_THAT(field.occupied_volume_m3,WithinAbs(1e-9,1e-20));
    };
    TH::HistoryLimits limits;
    limits.max_step_temperature_error_k = .001;
    const auto history = TH::evolve({{0.,cells,links,boundaries}},.75,limits,{},observations);
    REQUIRE(history.state == TH::State::Converged);
    CHECK(history.completed_observations == 3);
    CHECK(calls == 3);
    CHECK_THAT(history.energy_balance_error_j,WithinAbs(0.,1e-10));
}

TEST_CASE("A resolved thermal slab gradient converges toward free coupon curvature", "[ThermomechanicalHistory]")
{
    double previous_bending = 0.;
    for (const double h : {10., 100.}) {
        const double gradient = -(380. - 300.) / (height + conductivity / h);
        // The compatible continuum field has w_xx = -alpha*dT/dz.
        // This second difference removes the arbitrary rigid-body rotation.
        const double exact = -alpha * gradient * length * length / 4.;
        std::vector<double> errors;
        double bending = 0.;
        for (const int resolution : {4, 8, 12}) {
            const auto mesh = coupon(resolution, slab_temperatures(resolution, 380., 300., h));
            const auto result = free_equilibrium(mesh);
            bending = result.displacement_m[node_at(mesh, {0., 0., 0.})].z()
                - 2. * result.displacement_m[node_at(mesh, {length / 2., 0., 0.})].z()
                + result.displacement_m[node_at(mesh, {length, 0., 0.})].z();
            CAPTURE(h, resolution, bending, exact);
            CHECK(bending > 0.);
            errors.push_back(std::abs(bending - exact));
        }
        // Constant-strain tetrahedra cannot reproduce quadratic bending and
        // are artificially stiff on the coarse grid. Require refinement to
        // reduce the error substantially, and apply the accuracy criterion to
        // the resolved grid, not the deliberately under-resolved control.
        CHECK(errors.back() < .5 * errors.front());
        CHECK(errors[2] < errors[1]);
        CHECK_THAT(bending, WithinAbs(exact, .05 * exact));
        // Only this fixed-property, steady, free slab comparison is claimed;
        // it is not a fan-PWM law or a prediction of bonded-bed corner lift.
        CHECK(bending > previous_bending);
        previous_bending = bending;
    }
}

TEST_CASE("A solved thermal gradient drives section curvature without a prescribed shrinkage fit", "[ThermomechanicalHistory]")
{
    TE::LocalSectionTile tile;
    tile.index = {0,0}; tile.minimum = {0,0}; tile.maximum = {8000000,2000000};
    // Elevation excludes bed attachment: this checks free curvature, not an
    // adhesive strength, a crack initiation law or a fan PWM calibration.
    tile.slabs = {{height,2.*height,0,{Slic3r::ExPolygon(Slic3r::Polygon{
        {0,0},{8000000,0},{8000000,2000000},{0,2000000}})}}};
    const std::vector<TE::ExpansionInterval> expansion{{250.,400.,alpha*Eigen::Matrix3d::Identity()}};
    TE::SectionSolveLimits limits;
    limits.allow_free_components = true;
    double previous_curvature = 0.;
    for (const double h : {10.,100.}) {
        const auto temperatures = slab_temperatures(4,380.,300.,h);
        const double dz = height/4.;
        const double gradient = (temperatures.back()-temperatures.front())/(3.*dz);
        const auto strain = TE::linear_thermal_eigenstrain(expansion,Eigen::Matrix3d::Identity(),reference_k,
            {0.,0.,height+.5*dz},temperatures.front(),{0.,0.,gradient},
            {0.,0.,height},{length,width,2.*height});
        TE::SectionMaterial material{TE::ElasticityMatrix::Identity()*1e6,strain.at_world_zero,strain.gradient_per_m};
        const std::vector<TE::SectionSampleRequest> samples{
            {0,{.001,.001,.003},0},{1,{.004,.001,.003},0},{2,{.007,.001,.003},0}};
        const auto solved = TE::solve_local_section_tiles({tile},0.,{material},Eigen::Vector3d::Zero(),
            {.002,.001,.001},limits,{},{},samples);
        REQUIRE(solved.samples.size() == 3);
        CHECK(solved.quadratic_free_components == 1);
        const double curvature = (solved.samples[0].displacement_m.z()-2.*solved.samples[1].displacement_m.z()
            +solved.samples[2].displacement_m.z())/(.003*.003);
        const double analytical_gradient = -(380.-300.)/(height+conductivity/h);
        CHECK_THAT(curvature,WithinAbs(-alpha*analytical_gradient,1e-8));
        CHECK_THAT(solved.elastic_energy_j,WithinAbs(0.,1e-15));
        CHECK(curvature > previous_curvature);
        previous_curvature = curvature;
    }
}

TEST_CASE("Thermally driven section release follows linear displacement and quadratic energy scaling", "[ThermomechanicalHistory]")
{
    TE::LocalSectionTile tile;
    tile.index = {0,0}; tile.minimum = {0,0}; tile.maximum = {8000000,2000000};
    tile.slabs = {{0.,height,0,{Slic3r::ExPolygon(Slic3r::Polygon{
        {0,0},{8000000,0},{8000000,2000000},{0,2000000}})}}};
    const std::vector<TE::ExpansionInterval> expansion{{250.,400.,alpha*Eigen::Matrix3d::Identity()}};
    double previous_gradient = 0., previous_lift = 0., previous_release = 0.;
    for (const double h : {10.,100.}) {
        const auto temperatures = slab_temperatures(4,380.,300.,h);
        const double dz = height/4.;
        const double gradient = (temperatures.back()-temperatures.front())/(3.*dz);
        // Stress-free bed temperature removes a uniform offset so this
        // linear-heat-flow benchmark scales only the contraction gradient.
        const auto strain = TE::linear_thermal_eigenstrain(expansion,Eigen::Matrix3d::Identity(),380.,
            {0.,0.,.5*dz},temperatures.front(),{0.,0.,gradient},
            {0.,0.,0.},{length,width,height});
        TE::SectionMaterial material{TE::ElasticityMatrix::Identity()*1e6,strain.at_world_zero,strain.gradient_per_m};
        const auto bonded = TE::solve_local_section_tiles({tile},0.,{material},Eigen::Vector3d::Zero(),{.001,.001,.0005});
        const auto released = TE::solve_local_section_tiles({tile},0.,{material},Eigen::Vector3d::Zero(),
            {.001,.001,.0005},{},{},{},{},{},{{{.002,0.},{length,width}}},TE::SectionCellContact{});
        const double release = bonded.elastic_energy_j-released.elastic_energy_j;
        CHECK(release > 0.);
        CHECK(released.maximum_footprint_lift_m > 0.);
        if (previous_gradient != 0.) {
            const double ratio = gradient/previous_gradient;
            CHECK_THAT(released.maximum_footprint_lift_m,WithinRel(previous_lift*ratio,1e-5));
            CHECK_THAT(release,WithinRel(previous_release*ratio*ratio,1e-5));
        }
        previous_gradient = gradient;
        previous_lift = released.maximum_footprint_lift_m;
        previous_release = release;
    }
}
