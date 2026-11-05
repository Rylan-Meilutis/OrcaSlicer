#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/libslic3r.h"
#include "libslic3r/ThermoelasticLayers.hpp"
#include "libslic3r/ThermoelasticLocalSections.hpp"
#include "libslic3r/ThermoelasticLocalTopology.hpp"
#include "libslic3r/ThermoelasticMaterial.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticSections.hpp"
#include "libslic3r/ThermoelasticExactSections.hpp"
#include "libslic3r/ThermoelasticSectionFootprint.hpp"
#include "libslic3r/ThermoelasticThermalStrain.hpp"
#include "libslic3r/ThermoelasticSectionContact.hpp"
#include "libslic3r/ThermoelasticSectionTopology.hpp"
#include "libslic3r/ThermoelasticSolver.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cstddef>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
ExPolygon square() { return ExPolygon(Polygon{{0, 0}, {1000000000, 0}, {1000000000, 1000000000}, {0, 1000000000}}); }
SectionMaterial material()
{
    ElasticityMatrix matrix = ElasticityMatrix::Zero();
    matrix.diagonal().head<3>().setConstant(1000.);
    matrix.diagonal().tail<3>().setConstant(500.);
    return {matrix};
}
} // namespace

TEST_CASE("Exact section mechanics integrates rational support without moving the clipping cell", "[ThermoelasticSections]")
{
    const auto triangle=make_exact_section({ExPolygon(Polygon{{100000000,100000000},
        {1100000001,100000000},{1100000001,500000000}})});
    const auto strip=make_exact_section({ExPolygon(Polygon{{100000000,100000000},
        {600000000,100000000},{600000000,500000000},{100000000,500000000}})});
    const auto cut=boolean_exact_sections(triangle.section,strip.section,ExactSectionOperation::Intersection);
    const auto topology=qualify_section_topology({{0.,1.,0,SectionFootprint(cut.section)}},0.);
    auto elastic=material();
    elastic.eigenstrain[2]=-.001;
    const std::vector<SectionSampleRequest> samples{{1,{.2,.11,.5},0},{2,{.2,.1,.5},0}};
    const auto result=solve_deposited_sections(topology,{elastic},Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones(),{},{},samples);
    const double volume=50000000./1000000001.;
    CHECK_THAT(result.occupied_volume_m3,WithinRel(volume,1e-12));
    CHECK_THAT(result.bed_area_m2,WithinRel(volume,1e-12));
    CHECK_THAT(result.elastic_energy_j,WithinAbs(0.,1e-12));
    REQUIRE(result.samples.size()==1);
    CHECK(result.samples.front().id==1);
    CHECK(result.unsampled_points==1);
    for (std::size_t i=0;i<result.nodes_m.size();++i)
        CHECK_THAT((result.displacement_m[i]-Eigen::Vector3d(0.,0.,-.001*result.nodes_m[i].z())).norm(),WithinAbs(0.,1e-11));
    SectionSolveLimits limits;
    limits.max_work=100;
    CHECK_THROWS_AS(solve_deposited_sections(topology,{elastic},Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones(),limits),ResourceLimit);
    CHECK_THROWS_AS(solve_deposited_sections(topology,{elastic},Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones(),{},[] { throw std::runtime_error("cancel exact mechanics"); }),std::runtime_error);
}

TEST_CASE("Cut bed trace extrema include edge interiors but exclude empty cell corners", "[ThermoelasticSections]")
{
    const ExPolygons triangle{ExPolygon(Polygon{{0,0},{1000,0},{0,1000}})};
    const auto saddle = section_trace_extrema(triangle,{0,0},{1000,1000},{0.,0.,0.,1.});
    CHECK_THAT(saddle.minimum,WithinAbs(0.,1e-15));
    CHECK_THAT(saddle.maximum,WithinAbs(.25,1e-15));
    CHECK_THAT((saddle.maximum_local-Eigen::Vector2d(.5,.5)).norm(),WithinAbs(0.,1e-15));
    // The excluded (1,1) corner is below the plate; the actual triangular
    // footprint is nonpenetrating, with its minimum inside the diagonal edge.
    const auto gap = section_trace_extrema(triangle,{0,0},{1000,1000},{1.,1.,1.,-3.});
    CHECK_THAT(gap.minimum,WithinAbs(0.,1e-15));
    CHECK_THAT(gap.maximum,WithinAbs(1.,1e-15));
    CHECK_THAT((gap.minimum_local-Eigen::Vector2d(.5,.5)).norm(),WithinAbs(0.,1e-15));
    CHECK_THROWS_AS(section_trace_extrema(triangle,{0,0},{999,1000},{0.,0.,0.,1.}),std::invalid_argument);
    CHECK_THROWS_AS(section_trace_extrema({}, {0,0},{1000,1000},{0.,0.,0.,1.}),std::invalid_argument);
    CHECK_THROWS_AS(section_trace_extrema(triangle,{0,0},{1000,1000},{0.,0.,0.,1.},[] {
        throw std::runtime_error("cancel trace");
    }),std::runtime_error);
}

TEST_CASE("Spatial section loads retain geometry and match material-side loading", "[ThermoelasticSections]")
{
    LocalSectionTile tile{{0,0},{0,0},{1000000000,1000000000},{{0.,1.,0,{square()}}}};
    auto lower = material(), upper = material();
    lower.eigenstrain.head<3>().setConstant(-.001);
    upper.eigenstrain.head<3>().setConstant(-.002);
    const std::vector<SectionSampleRequest> samples{{0,{.3,.3,.25},0},{1,{.3,.3,.75},0}};
    std::size_t calls = 0;
    const SectionLoadField load = [&](const SectionLoadQuery& q) {
        ++calls;
        CHECK(q.material == 0);
        CHECK_THAT(q.occupied_volume_m3,WithinAbs(.125,1e-12));
        CHECK_THAT((q.maximum_m-q.minimum_m-Eigen::Vector3d(.5,.5,.5)).norm(),WithinAbs(0.,1e-12));
        return LinearThermalStrain{q.cell[2] == 0 ? lower.eigenstrain : upper.eigenstrain,
            Eigen::Matrix<double,6,3>::Zero()};
    };
    const auto actual = solve_local_section_tiles({tile},0.,{material()},Eigen::Vector3d::Zero(),
        {.5,.5,.5},{},{},{},samples,{},{},{},load);
    CHECK(actual.uses_spatial_load);
    CHECK(calls == 8);
    tile.slabs = {{0.,.5,0,{square()}},{.5,1.,1,{square()}}};
    auto reference_samples = samples;
    reference_samples[1].material = 1;
    const auto reference = solve_local_section_tiles({tile},0.,{lower,upper},Eigen::Vector3d::Zero(),
        {.5,.5,.5},{},{},{},reference_samples);
    CHECK_FALSE(reference.uses_spatial_load);
    CHECK_THAT(actual.occupied_volume_m3,WithinAbs(reference.occupied_volume_m3,1e-12));
    CHECK_THAT(actual.elastic_energy_j,WithinRel(reference.elastic_energy_j,1e-10));
    REQUIRE(actual.samples.size() == reference.samples.size());
    for (std::size_t i = 0; i < actual.samples.size(); ++i) {
        CHECK_THAT((actual.samples[i].displacement_m-reference.samples[i].displacement_m).norm(),WithinAbs(0.,1e-12));
        CHECK_THAT((actual.samples[i].stress_pa-reference.samples[i].stress_pa).norm(),WithinAbs(0.,1e-9));
    }
}

TEST_CASE("Spatial section loading rejects incomplete and nonfinite fields", "[ThermoelasticSections]")
{
    LocalSectionTile tile{{0,0},{0,0},{1000000000,1000000000},{{0.,1.,0,{square()}}}};
    const auto solve = [&](const SectionLoadField& load) {
        return solve_local_section_tiles({tile},0.,{material()},Eigen::Vector3d::Zero(),{.5,.5,.5},
            {},{},{},{},{},{},{},load);
    };
    CHECK_THROWS_AS(solve([](const SectionLoadQuery&) -> LinearThermalStrain {
        throw std::out_of_range("Missing thermal coverage");
    }),std::out_of_range);
    CHECK_THROWS_AS(solve([](const SectionLoadQuery&) {
        return LinearThermalStrain{Eigen::Matrix<double,6,1>::Constant(std::numeric_limits<double>::quiet_NaN()),
            Eigen::Matrix<double,6,3>::Zero()};
    }),std::invalid_argument);
    CHECK_THROWS_AS(solve([](const SectionLoadQuery&) {
        TrilinearThermalStrain strain;
        for (auto& corner : strain.at_corners) corner.setZero();
        strain.at_corners[7][0] = std::numeric_limits<double>::quiet_NaN();
        return strain;
    }),std::invalid_argument);
}

TEST_CASE("Trilinear section loads preserve compatible cross terms through assembly and sampling", "[ThermoelasticSections]")
{
    auto footprint = square();
    footprint.holes.push_back(Polygon{{600000000,600000000},{800000000,600000000},
        {800000000,800000000},{600000000,800000000}});
    footprint.holes.back().make_clockwise();
    const LocalSectionTile tile{{0,0},{0,0},{1000000000,1000000000},{{0.,1.,0,{footprint}}}};
    // Manufactured u=(0,0,a*x*y*z) is Q1 and zero on the fixed bed. Its
    // compatible engineering strain contains xy, xz and yz cross terms.
    const double a = .001;
    const SectionLoadField load = [a](const SectionLoadQuery& q) {
        TrilinearThermalStrain field;
        for (unsigned corner = 0; corner < 8; ++corner) {
            Eigen::Vector3d p;
            for (unsigned axis = 0; axis < 3; ++axis)
                p[axis] = corner & (1u << axis) ? q.maximum_m[axis] : q.minimum_m[axis];
            field.at_corners[corner] << 0.,0.,a*p.x()*p.y(),0.,a*p.x()*p.z(),a*p.y()*p.z();
        }
        return field;
    };
    const Eigen::Vector3d point(.3,.4,.6);
    for (double spacing : {1.,.5}) {
        auto stale = material();
        stale.eigenstrain.setConstant(.1); // Callback replaces, not adds to it.
        const auto result = solve_local_section_tiles({tile},0.,{stale},Eigen::Vector3d::Zero(),
            Eigen::Vector3d::Constant(spacing),{},{},{},{{0,point,0}},{},{},{},load);
        REQUIRE(result.samples.size() == 1);
        CHECK_THAT(result.elastic_energy_j,WithinAbs(0.,1e-12));
        CHECK_THAT(result.samples[0].stress_pa.norm(),WithinAbs(0.,1e-9));
        CHECK_THAT((result.samples[0].displacement_m-Eigen::Vector3d(0.,0.,a*point.prod())).norm(),WithinAbs(0.,1e-12));
    }
    auto floating = tile;
    floating.slabs.front().bottom_m = 1.;
    floating.slabs.front().top_m = 2.;
    SectionSolveLimits limits;
    limits.allow_free_components = true;
    const auto free = solve_local_section_tiles({floating},0.,{material()},Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones(),limits,{},{},{{0,{.3,.4,1.6},0}},{},{},{},load);
    REQUIRE(free.samples.size() == 1);
    CHECK(free.affine_free_components == 0);
    CHECK(free.quadratic_free_components == 0);
    CHECK_THAT(free.samples[0].stress_pa.norm(),WithinAbs(0.,1e-9));
    const coord_t x=180001000,y=110001000,z=4000100,side=1000;
    const ExPolygon tiny(Polygon{{x,y},{x+side,y},{x+side,y+side},{x,y+side}});
    const LocalSectionTile small{{0,0},{0,0},{1000000000,1000000000},
        {{double(z)*1e-9,double(z+side)*1e-9,0,{tiny}}}};
    const Eigen::Vector3d point_small(double(x+side/2)*1e-9,double(y+side/2)*1e-9,double(z+side/2)*1e-9);
    const auto small_field=solve_local_section_tiles({small},0.,{material()},Eigen::Vector3d::Zero(),
        {.005,.005,.002},limits,{},{},{{0,point_small,0}},{},{},{},load);
    REQUIRE(small_field.samples.size()==1);
    CHECK_THAT(small_field.samples.front().stress_pa.norm(),WithinAbs(0.,1e-9));
    CHECK_THAT(small_field.elastic_energy_j,WithinAbs(0.,small_field.occupied_volume_m3*1e-12));
}

TEST_CASE("Free section components use spatial rather than stale material strain", "[ThermoelasticSections]")
{
    LocalSectionTile tile{{0,0},{0,0},{1000000000,1000000000},{{1.,2.,0,{square()}}}};
    auto elastic = material();
    elastic.eigenstrain.head<3>().setConstant(-.001);
    elastic.eigenstrain_gradient_per_m(0,2) = .0002;
    SectionSolveLimits limits;
    limits.allow_free_components = true;
    const std::vector<SectionSampleRequest> samples{{0,{.3,.3,1.3},0}};
    const auto expected = solve_local_section_tiles({tile},0.,{elastic},Eigen::Vector3d::Zero(),
        {.5,.5,.5},limits,{},{},samples);
    const SectionLoadField load = [&](const SectionLoadQuery&) {
        return LinearThermalStrain{elastic.eigenstrain,elastic.eigenstrain_gradient_per_m};
    };
    const auto actual = solve_local_section_tiles({tile},0.,{material()},Eigen::Vector3d::Zero(),
        {.5,.5,.5},limits,{},{},samples,{},{},{},load);
    CHECK(actual.quadratic_free_components == 1);
    REQUIRE(actual.samples.size() == 1);
    REQUIRE(expected.samples.size() == 1);
    CHECK_THAT((actual.samples[0].displacement_m-expected.samples[0].displacement_m).norm(),WithinAbs(0.,1e-14));
    CHECK_THAT(actual.samples[0].stress_pa.norm(),WithinAbs(0.,1e-14));
    REQUIRE(actual.displacement_m.size() == expected.displacement_m.size());
    for (std::size_t i = 0; i < actual.displacement_m.size(); ++i)
        CHECK_THAT((actual.displacement_m[i]-expected.displacement_m[i]).norm(),WithinAbs(0.,1e-14));
}

TEST_CASE("Fixed-grid bed reactions preserve force moment and coarse virtual work", "[ThermoelasticSections]")
{
    SectionMechanicalResult field;
    field.basis_origin_m = {-.5,-.5,0.};
    field.basis_spacing_m = {.5,.5,.5};
    const auto virtual_displacement = [](const Eigen::Vector3d& p) {
        return Eigen::Vector3d(p.x()*p.y(), p.x()+2.*p.y(), 1.+p.x()-p.y());
    };
    Eigen::Vector3d force = Eigen::Vector3d::Zero(), moment = Eigen::Vector3d::Zero();
    double work = 0.;
    for (int x = 0; x < 3; ++x)
        for (int y = 0; y < 3; ++y) {
            const Eigen::Vector3d p(.5*x-.5,.5*y-.5,0.), r(2.*x,-double(y),double(x-y));
            field.nodes_m.push_back(p);
            field.generalized_bed_reaction_n.push_back(r);
            field.bed_connected_nodes.push_back(true);
            force += r;
            moment += p.cross(r);
            work += r.dot(virtual_displacement(p));
        }
    const auto projected = project_section_bed_reactions(field,{0.,0.},{1.,1.});
    Eigen::Vector3d projected_force = Eigen::Vector3d::Zero(), projected_moment = Eigen::Vector3d::Zero();
    double projected_work = 0.;
    for (std::size_t i = 0; i < projected.nodes_m.size(); ++i) {
        const auto& r = projected.generalized_reaction_n[i];
        const auto& p = projected.nodes_m[i];
        projected_force += r;
        projected_moment += p.cross(r);
        projected_work += r.dot(virtual_displacement(p));
    }
    CHECK_THAT((projected_force-force).norm(),WithinAbs(0.,1e-12));
    CHECK_THAT((projected_moment-moment).norm(),WithinAbs(0.,1e-12));
    CHECK_THAT(projected_work,WithinAbs(work,1e-12));
    CHECK_THROWS_AS(project_section_bed_reactions(field,{0.,0.},{.75,1.}),std::invalid_argument);
    CHECK_THROWS_AS(project_section_bed_reactions(field,{.25,0.},{1.,1.}),std::invalid_argument);
    CHECK_THROWS_AS(project_section_bed_reactions(field,{0.,0.},{1.,1.},[] {
        throw ResourceLimit("cancel projection");
    }),ResourceLimit);
    field.nodes_m.back().z() = .5;
    CHECK_THROWS_AS(project_section_bed_reactions(field,{0.,0.},{1.,1.}),std::invalid_argument);
}

TEST_CASE("Fixed-grid reaction measurement keeps unresolved sign cancellation explicit", "[ThermoelasticSections]")
{
    SectionMechanicalResult field;
    field.basis_spacing_m = {.25,.25,.25};
    field.nodes_m = {{.25,0.,0.},{.5,0.,0.},{.75,0.,0.}};
    field.generalized_bed_reaction_n = {{0.,0.,-1.},{0.,0.,2.},{0.,0.,-1.}};
    field.bed_connected_nodes.assign(3,true);
    const auto fine = project_section_bed_reactions(field,{0.,0.},{.25,.25});
    const auto coarse = project_section_bed_reactions(field,{0.,0.},{1.,1.});
    CHECK_THAT(fine.tensile_generalized_force_n,WithinAbs(2.,1e-12));
    CHECK_THAT(coarse.tensile_generalized_force_n,WithinAbs(0.,1e-12));
    // Coarse zero is deliberately NOT evidence of zero local adhesion demand.
}

TEST_CASE("Graded section planes preserve uniform-grid solutions and nested energy bounds", "[ThermoelasticSections]")
{
    LocalSectionTile tile;
    tile.index = {0,0};
    tile.minimum = {0,0};
    tile.maximum = {1000000000,1000000000};
    tile.slabs = {{0.,1.,0,{square()}}};
    auto elastic = material();
    elastic.eigenstrain.head<3>().setConstant(-.001);
    const std::vector<SectionSampleRequest> samples{{1,{.3,.3,.3},0},{2,{.3,.3,2.},0}};
    const auto solve = [&](double dz, const std::vector<double>& planes) {
        return solve_local_section_tiles({tile},0.,{elastic},Eigen::Vector3d::Zero(),
            {.5,.5,dz},{},{},{},samples,planes);
    };
    const auto regular = solve(.5,{}), explicit_uniform = solve(.5,{0.,.5,1.});
    CHECK_THAT(explicit_uniform.elastic_energy_j,WithinRel(regular.elastic_energy_j,1e-12));
    REQUIRE(explicit_uniform.nodes_m.size() == regular.nodes_m.size());
    for (std::size_t i = 0; i < regular.nodes_m.size(); ++i) {
        CHECK_THAT((explicit_uniform.displacement_m[i]-regular.displacement_m[i]).norm(),WithinAbs(0.,1e-12));
        CHECK_THAT((explicit_uniform.generalized_bed_reaction_n[i]-regular.generalized_bed_reaction_n[i]).norm(),WithinAbs(0.,1e-10));
    }
    const auto graded = solve(.5,{0.,.25,.5,1.}), fine = solve(.25,{});
    CHECK(graded.nodes_m.size() < fine.nodes_m.size());
    CHECK(graded.elastic_energy_j <= regular.elastic_energy_j * (1.+1e-10));
    CHECK(fine.elastic_energy_j <= graded.elastic_energy_j * (1.+1e-10));
    CHECK_THAT(graded.occupied_volume_m3,WithinAbs(1.,1e-12));
    CHECK_THAT(graded.bed_area_m2,WithinAbs(1.,1e-12));
    REQUIRE(graded.samples.size() == 1);
    CHECK(graded.unsampled_points == 1);
    CHECK(graded.relative_residual < 1e-8);
    for (const auto& planes : std::vector<std::vector<double>>{{0.},{0.,.5,.5,1.},{0.,1.,.5},
            {.1,.5,1.},{0.,.5},{0.,.5000000001,1.},std::vector<double>(129,0.)})
        CHECK_THROWS_AS(solve(.5,planes),std::invalid_argument);
}

TEST_CASE("Section factor budgets reject fill before returning a field without changing equilibrium", "[ThermoelasticSections]")
{
    LocalSectionTile tile;
    tile.index = {0,0};
    tile.minimum = {0,0};
    tile.maximum = {1000000000,1000000000};
    tile.slabs = {{0.,1.,0,{square()}}};
    auto elastic = material();
    elastic.eigenstrain.head<3>().setConstant(-.001);
    SectionSolveLimits limits;
    const auto solve = [&] {
        return solve_local_section_tiles({tile},0.,{elastic},Eigen::Vector3d::Zero(),{.25,.25,.25},limits);
    };
    const auto reference = solve();
    REQUIRE(reference.factor_entries > 1);
    REQUIRE(reference.solved_dofs > 0);
    CHECK(reference.assembly_triplets >= reference.matrix_entries);
    CHECK(reference.matrix_entries > 0);
    limits.sparse_factor = SectionFactorLimits{};
    limits.sparse_factor->max_entries = reference.factor_entries;
    const auto bounded = solve();
    CHECK(bounded.factor_entries == reference.factor_entries);
    CHECK(bounded.solved_dofs == reference.solved_dofs);
    CHECK(bounded.assembly_triplets == reference.assembly_triplets);
    CHECK(bounded.matrix_entries == reference.matrix_entries);
    CHECK_THAT(bounded.elastic_energy_j,WithinRel(reference.elastic_energy_j,1e-12));
    for (std::size_t i = 0; i < reference.nodes_m.size(); ++i) {
        CHECK_THAT((bounded.displacement_m[i]-reference.displacement_m[i]).norm(),WithinAbs(0.,1e-12));
        CHECK_THAT((bounded.generalized_bed_reaction_n[i]-reference.generalized_bed_reaction_n[i]).norm(),WithinAbs(0.,1e-10));
    }
    --limits.sparse_factor->max_entries;
    CHECK_THROWS_AS(solve(),ResourceLimit);
    limits.sparse_factor = SectionFactorLimits{};
    limits.sparse_factor->max_work = 1;
    CHECK_THROWS_AS(solve(),ResourceLimit);
    limits.sparse_factor->max_work = 256000001;
    CHECK_THROWS_AS(solve(),std::invalid_argument);
    limits.sparse_factor = SectionFactorLimits{};
    limits.max_nodes = 96000;
    CHECK_NOTHROW(solve());
    limits.sparse_factor.reset();
    CHECK_THROWS_AS(solve(),std::invalid_argument);
}

TEST_CASE("Exact released bed traces match native lift contact and reactions around holes", "[ThermoelasticSections]")
{
    for (const bool hole:{false,true}) for (const double contraction:{-.001,.001}) {
        CAPTURE(hole,contraction);
        auto shape=square();
        if (hole) shape.holes.emplace_back(Polygon{{250000000,250000000},{250000000,750000000},
            {750000000,750000000},{750000000,250000000}});
        auto lower=material(),upper=material();
        upper.eigenstrain[0]=contraction;
        const auto retained=make_exact_section({shape}).section;
        SectionSolveLimits limits;
        limits.allow_free_components=true;
        const std::vector<SectionBedReleaseRegion> release{{{.5,0.},{1.,1.}}};
        const auto solve=[&](bool exact,bool released) {
            const auto footprint=exact ? SectionFootprint(retained) : SectionFootprint(ExPolygons{shape});
            LocalSectionTile tile{{0,0},{0,0},{1000000000,1000000000},
                {{0.,.5,0,footprint},{.5,1.,1,footprint}}};
            return solve_local_section_tiles({tile},0.,{lower,upper},Eigen::Vector3d::Zero(),
                {.25,.25,.25},limits,{},{},{},{},released ? release : std::vector<SectionBedReleaseRegion>{},
                released ? std::optional<SectionCellContact>(SectionCellContact{}) : std::nullopt);
        };
        for (const bool released:{false,true}) {
            CAPTURE(released);
            const auto native=solve(false,released),exact=solve(true,released);
            if (released) {
                const auto topology=qualify_section_topology({{0.,.5,0,SectionFootprint(retained)},
                    {.5,1.,1,SectionFootprint(retained)}},0.);
                const auto retained_result=solve_deposited_sections(topology,{lower,upper},Eigen::Vector3d::Zero(),
                    {.25,.25,.25},limits,{},{},release,SectionCellContact{});
                CHECK_THAT(retained_result.elastic_energy_j,WithinAbs(exact.elastic_energy_j,1e-10));
                CHECK_THAT(retained_result.maximum_footprint_lift_m,WithinAbs(exact.maximum_footprint_lift_m,1e-10));
                CHECK_THAT(retained_result.released_bed_area_m2,WithinAbs(exact.released_bed_area_m2,1e-12));
            }
            CHECK_THAT(exact.elastic_energy_j,WithinAbs(native.elastic_energy_j,1e-10));
            CHECK_THAT(exact.bed_area_m2,WithinAbs(native.bed_area_m2,1e-12));
            CHECK_THAT(exact.released_bed_area_m2,WithinAbs(native.released_bed_area_m2,1e-12));
            CHECK_THAT(exact.maximum_footprint_lift_m,WithinAbs(native.maximum_footprint_lift_m,1e-10));
            CHECK_THAT(exact.maximum_footprint_penetration_m,WithinAbs(native.maximum_footprint_penetration_m,1e-10));
            CHECK(exact.conservative_cell_contact==native.conservative_cell_contact);
            REQUIRE(exact.nodes_m.size()==native.nodes_m.size());
            for (std::size_t i=0;i<exact.nodes_m.size();++i) {
                CHECK_THAT((exact.nodes_m[i]-native.nodes_m[i]).norm(),WithinAbs(0.,1e-12));
                CHECK_THAT((exact.generalized_bed_reaction_n[i]-native.generalized_bed_reaction_n[i]).norm(),WithinAbs(0.,1e-9));
            }
        }
        limits.max_bed_trace_points=1;
        CHECK_THROWS_AS(solve(true,true),ResourceLimit);
    }
}

TEST_CASE("Exact refinement scales preserve area and propagate resource limits", "[ThermoelasticSections]")
{
    auto elastic=material();
    elastic.eigenstrain[2]=-.001;
    const auto exact=make_exact_section({square()}).section;
    const auto topology=qualify_section_topology({{0.,1.,0,SectionFootprint(exact)}},0.);
    const auto native=qualify_section_topology({{0.,1.,0,{square()}}},0.);
    const auto solve=[&](const SectionTopology& input,const SectionSolveLimits& limits) {
        return refine_deposited_sections(input,{elastic},Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),3,.01,limits);
    };
    const auto baseline=solve(native,{}),retained=solve(topology,{});
    CHECK(retained.state==baseline.state);
    SectionSolveLimits limits;
    limits.max_work=1;
    CHECK_THROWS_AS(solve(topology,limits),ResourceLimit);
    limits.max_work=1000000;
    limits.max_points=1;
    CHECK_THROWS_AS(solve(topology,limits),ResourceLimit);
    CHECK_THROWS_AS(refine_deposited_sections(topology,{elastic},Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),3,.01,{},
        {},[] { throw std::runtime_error("cancel exact refinement"); }),std::runtime_error);
}

TEST_CASE("Prescribed section bed release lowers stored energy and permits upward curl", "[ThermoelasticSections]")
{
    LocalSectionTile tile;
    tile.index = {0,0}; tile.minimum = {0,0}; tile.maximum = {1000000000,1000000000};
    tile.slabs = {{0.,.5,0,{square()}},{.5,1.,1,{square()}}};
    auto bottom = material(), top = material();
    top.eigenstrain[0] = -.001;
    const std::vector<SectionBedReleaseRegion> release{{{.5,0.},{1.,1.}}};
    SectionSolveLimits limits;
    limits.allow_free_components = true;
    std::optional<SectionCellContact> contact;
    const auto solve = [&](const std::vector<SectionBedReleaseRegion>& regions) {
        return solve_local_section_tiles({tile},0.,{bottom,top},Eigen::Vector3d::Zero(),
            {.25,.25,.25},limits,{},{},{},{},regions,contact);
    };
    const auto bonded = solve({}), lifted = solve(release);
    CHECK_THAT(lifted.bed_area_m2,WithinAbs(1.,1e-12));
    CHECK_THAT(lifted.released_bed_area_m2,WithinAbs(.5,1e-12));
    CHECK_THAT(lifted.occupied_volume_m3,WithinRel(bonded.occupied_volume_m3,1e-12));
    CHECK(lifted.elastic_energy_j < bonded.elastic_energy_j);
    CHECK(lifted.bed_lift_upper_bound_m > 0.);
    CHECK_THAT(lifted.maximum_footprint_lift_m,WithinRel(lifted.bed_lift_upper_bound_m,1e-12));
    CHECK(lifted.free_components == 0);
    for (std::size_t i = 0; i < lifted.nodes_m.size(); ++i) {
        if (lifted.nodes_m[i].z() != 0.) continue;
        if (lifted.nodes_m[i].x() <= .5)
            CHECK_THAT(lifted.displacement_m[i].norm(),WithinAbs(0.,1e-12));
        else
            CHECK_THAT(lifted.generalized_bed_reaction_n[i].norm(),WithinAbs(0.,1e-12));
    }
    const auto duplicate = solve({release.front(),release.front()});
    CHECK_THAT(duplicate.released_bed_area_m2,WithinAbs(.5,1e-12));
    CHECK_THAT(duplicate.elastic_energy_j,WithinRel(lifted.elastic_energy_j,1e-12));
    contact = SectionCellContact{};
    const auto noncontacting = solve(release);
    CHECK(noncontacting.contact_iterations == 1);
    CHECK_FALSE(noncontacting.conservative_cell_contact);
    CHECK_THAT(noncontacting.elastic_energy_j,WithinRel(lifted.elastic_energy_j,1e-12));
    CHECK_THAT(noncontacting.bed_lift_upper_bound_m,WithinRel(lifted.bed_lift_upper_bound_m,1e-12));
    contact.reset();
    top.eigenstrain *= .5;
    const auto half = solve(release);
    CHECK_THAT(half.bed_lift_upper_bound_m,WithinRel(.5*lifted.bed_lift_upper_bound_m,1e-10));
    CHECK_THAT(half.elastic_energy_j,WithinRel(.25*lifted.elastic_energy_j,1e-10));
    top.eigenstrain *= -1.;
    CHECK_THROWS_AS(solve(release),std::runtime_error); // Curvature into the bed needs contact.
    const auto reverse_bonded = solve({});
    contact = SectionCellContact{};
    const auto pressed = solve(release);
    CHECK(pressed.conservative_cell_contact);
    CHECK(pressed.contact_iterations > 1);
    CHECK(pressed.elastic_energy_j < reverse_bonded.elastic_energy_j);
    double normal_balance = 0.;
    bool tangential_motion = false;
    for (std::size_t i = 0; i < pressed.nodes_m.size(); ++i) {
        normal_balance += pressed.generalized_bed_reaction_n[i].z()+pressed.generalized_contact_reaction_n[i];
        CHECK(pressed.generalized_contact_reaction_n[i] >= -1e-10);
        if (pressed.nodes_m[i].z() != 0.) continue;
        CHECK(pressed.displacement_m[i].z() >= -1e-12);
        if (pressed.generalized_contact_reaction_n[i] > 1e-10)
            CHECK_THAT(pressed.displacement_m[i].z(),WithinAbs(0.,1e-12));
        if (pressed.nodes_m[i].x() > .5 && std::abs(pressed.displacement_m[i].x()) > 1e-6)
            tangential_motion = true;
    }
    CHECK(tangential_motion);
    CHECK_THAT(normal_balance,WithinAbs(0.,1e-10));
    contact->max_iterations = 1;
    CHECK_THROWS_AS(solve(release),ResourceLimit);
    contact.reset();
    CHECK_THROWS_AS(solve({{{0.,0.},{1.,1.}}}),std::runtime_error);
    CHECK_THROWS_AS(solve({{{.3,0.},{1.,1.}}}),std::invalid_argument);
    CHECK_THROWS_AS(solve({{{2.,0.},{3.,1.}}}),std::invalid_argument);
    limits.max_bed_trace_points = 1;
    CHECK_THROWS_AS(solve(release),ResourceLimit);
    limits.max_bed_trace_points = 200000;
    top.eigenstrain.setZero();
    const auto unloaded = solve(release);
    CHECK_THAT(unloaded.elastic_energy_j,WithinAbs(0.,1e-15));
    CHECK_THAT(unloaded.bed_lift_upper_bound_m,WithinAbs(0.,1e-15));
    bottom.eigenstrain[2] = top.eigenstrain[2] = -.001;
    contact = SectionCellContact{};
    const auto compatible = solve(release);
    CHECK(compatible.contact_iterations == 1);
    CHECK_THAT(compatible.elastic_energy_j,WithinAbs(0.,1e-15));
    for (std::size_t i = 0; i < compatible.nodes_m.size(); ++i) {
        const Eigen::Vector3d expected(0.,0.,-.001*compatible.nodes_m[i].z());
        CHECK_THAT((compatible.displacement_m[i]-expected).norm(),WithinAbs(0.,1e-12));
        CHECK_THAT(compatible.generalized_contact_reaction_n[i],WithinAbs(0.,1e-12));
    }
}

TEST_CASE("Released section strip approaches the thermal laminate bending solution", "[ThermoelasticSections]")
{
    const ExPolygon footprint(Polygon{{0,0},{1000000000,0},{1000000000,100000000},{0,100000000}});
    LocalSectionTile tile;
    tile.index = {0,0}; tile.minimum = {0,0}; tile.maximum = {1000000000,100000000};
    tile.slabs = {{0.,.01,0,{footprint}},{.01,.02,1,{footprint}}};
    auto bottom = material(), top = material();
    top.eigenstrain[0] = -.001;
    // Equal-modulus layers, nu=0: integral(z*eigenstrain)/integral(z*z)
    // gives natural curvature 3*delta_strain/(2*thickness). Euler-Bernoulli
    // free-end lift is kappa*length^2/2. The finite anchored base here is
    // elastic rather than an ideal clamp, so only its slender-strip limit agrees.
    const double beam_lift = .5 * (3.*.001/(2.*.02)) * .9*.9;
    double previous_energy = 1.;
    double last_lift = 0.;
    SectionSolveLimits limits;
    limits.max_cells = 2000;
    limits.max_free_dofs = 12000;
    for (double dx : {.05,.025,.0125,.00625}) {
        const auto field = solve_local_section_tiles({tile},0.,{bottom,top},Eigen::Vector3d::Zero(),
            {dx,.05,.005},limits,{},{},{},{},{{{.1,0.},{1.,.1}}});
        CAPTURE(dx,field.bed_lift_upper_bound_m,beam_lift,field.elastic_energy_j);
        CHECK(field.elastic_energy_j <= previous_energy*(1.+1e-10));
        CHECK(field.bed_lift_upper_bound_m > last_lift);
        previous_energy = field.elastic_energy_j;
        last_lift = field.bed_lift_upper_bound_m;
    }
    CHECK_THAT(last_lift,WithinRel(beam_lift,.1));
}

TEST_CASE("Linear vertical eigenstrain recovers quadratic displacement and its exact discretization energy", "[ThermoelasticSections]")
{
    LocalSectionTile tile;
    tile.index = {0,0}; tile.minimum = {0,0}; tile.maximum = {1000000000,1000000000};
    tile.slabs = {{0.,1.,0,{square()}}};
    auto elastic = material();
    constexpr double gradient = -.001;
    elastic.eigenstrain_gradient_per_m(2,2) = gradient;
    for (const auto& planes : std::vector<std::vector<double>>{{0.,.5,1.},{0.,.25,.5,.75,1.},{0.,.125,.5,1.}}) {
        std::vector<SectionSampleRequest> samples;
        double sum_cubes = 0.;
        for (std::size_t i = 1; i < planes.size(); ++i) {
            const double dz = planes[i]-planes[i-1];
            sum_cubes += dz*dz*dz;
            samples.push_back({i,{.3,.3,(planes[i]+planes[i-1])*.5},0});
        }
        const auto field = solve_local_section_tiles({tile},0.,{elastic},Eigen::Vector3d::Zero(),
            {.5,.5,.5},{},{},{},samples,planes);
        // nu=0 and epsilon_zz=g*z: the exact displacement is g*z^2/2.
        // Q1 recovers its nodal values; the residual strain is g*(mid-z).
        // Integrating E/2 times its square gives E*g^2*sum(dz^3)/24.
        CHECK_THAT(field.elastic_energy_j,WithinRel(1000.*gradient*gradient*sum_cubes/24.,1e-9));
        CHECK_THAT(field.prescribed_strain_energy_j,WithinRel(1000.*gradient*gradient/6.,1e-12));
        CHECK_THAT(field.prescribed_stress_bound_pa,WithinRel(1.,1e-12));
        for (std::size_t i = 0; i < field.nodes_m.size(); ++i) {
            const double z = field.nodes_m[i].z();
            CHECK_THAT((field.displacement_m[i]-Eigen::Vector3d(0.,0.,gradient*z*z*.5)).norm(),WithinAbs(0.,1e-12));
        }
        REQUIRE(field.samples.size() == samples.size());
        for (const auto& sample : field.samples)
            CHECK_THAT(sample.stress_pa.norm(),WithinAbs(0.,1e-10));
    }
    elastic.eigenstrain_gradient_per_m(2,2) = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(solve_local_section_tiles({tile},0.,{elastic},Eigen::Vector3d::Zero(),{.5,.5,.5}),std::invalid_argument);
}

TEST_CASE("A released strip with linear through-thickness shrinkage approaches its natural curvature", "[ThermoelasticSections]")
{
    const ExPolygon footprint(Polygon{{0,0},{1000000000,0},{1000000000,100000000},{0,100000000}});
    LocalSectionTile tile;
    tile.index = {0,0}; tile.minimum = {0,0}; tile.maximum = {1000000000,100000000};
    tile.slabs = {{0.,.02,0,{footprint}}};
    auto elastic = material();
    elastic.eigenstrain_gradient_per_m(0,2) = -.05;
    // epsilon_xx = -kappa*z: the stress-free slender strip has
    // u_z = kappa*x^2/2. Its remaining free span is 0.9 m.
    const double beam_lift = .05*.9*.9*.5;
    double previous_energy = 1., previous_lift = 0.;
    SectionSolveLimits limits;
    limits.max_cells = 2000;
    limits.max_free_dofs = 12000;
    for (double dx : {.025,.0125,.00625}) {
        const auto field = solve_local_section_tiles({tile},0.,{elastic},Eigen::Vector3d::Zero(),
            {dx,.05,.005},limits,{},{},{},{},{{{.1,0.},{1.,.1}}});
        CAPTURE(dx,field.maximum_footprint_lift_m,beam_lift);
        CHECK(field.elastic_energy_j < previous_energy);
        CHECK(field.maximum_footprint_lift_m > previous_lift);
        previous_energy = field.elastic_energy_j;
        previous_lift = field.maximum_footprint_lift_m;
    }
    CHECK_THAT(previous_lift,WithinRel(beam_lift,.05));
}

TEST_CASE("Linear eigenstrain exactly reproduces a compatible bilinear displacement on cut cells", "[ThermoelasticSections]")
{
    LocalSectionTile tile;
    tile.index = {0,0}; tile.minimum = {0,0}; tile.maximum = {1000000000,1000000000};
    const ExPolygon triangle(Polygon{{0,0},{1000000000,0},{0,1000000000}});
    tile.slabs = {{0.,1.,0,{triangle}}};
    auto elastic = material();
    constexpr double a = -.001, b = -.002;
    // u=(a*x*z,b*y*z,0) is exactly Q1 and vanishes everywhere on the bed.
    // The xz/yz entries are engineering shear, not tensor shear.
    elastic.eigenstrain_gradient_per_m(0,2) = a;
    elastic.eigenstrain_gradient_per_m(1,2) = b;
    elastic.eigenstrain_gradient_per_m(5,0) = a;
    elastic.eigenstrain_gradient_per_m(4,1) = b;
    const auto field = solve_local_section_tiles({tile},0.,{elastic},Eigen::Vector3d::Zero(),
        {.5,.5,.5},{},{},{},{{1,{.2,.3,.4},0},{2,{.6,.1,.7},0}});
    CHECK_THAT(field.occupied_volume_m3,WithinRel(.5,1e-12));
    CHECK_THAT(field.elastic_energy_j,WithinAbs(0.,1e-12));
    for (std::size_t i = 0; i < field.nodes_m.size(); ++i) {
        const auto& p = field.nodes_m[i];
        CHECK_THAT((field.displacement_m[i]-Eigen::Vector3d(a*p.x()*p.z(),b*p.y()*p.z(),0.)).norm(),WithinAbs(0.,1e-12));
    }
    REQUIRE(field.samples.size() == 2);
    for (const auto& sample : field.samples)
        CHECK_THAT(sample.stress_pa.norm(),WithinAbs(0.,1e-10));
}

TEST_CASE("Linear section eigenstrain uses world coordinates without the uniform free-body shortcut", "[ThermoelasticSections]")
{
    LocalSectionTile tile;
    tile.index = {0,0}; tile.minimum = {0,0}; tile.maximum = {1000000000,1000000000};
    tile.slabs = {{1.,2.,0,{square()}}};
    auto elastic = material();
    elastic.eigenstrain_gradient_per_m(2,2) = -.001;
    const auto shifted = solve_local_section_tiles({tile},1.,{elastic},{0.,0.,1.},{.5,.5,.25});
    for (std::size_t i = 0; i < shifted.nodes_m.size(); ++i) {
        const double z = shifted.nodes_m[i].z();
        CHECK_THAT(shifted.displacement_m[i].z(),WithinAbs(-.0005*(z*z-1.),1e-12));
    }
    SectionSolveLimits limits;
    limits.allow_free_components = true;
    const auto free = solve_local_section_tiles({tile},0.,{elastic},Eigen::Vector3d::Zero(),{.5,.5,.25},limits);
    CHECK(free.free_components == 1);
    CHECK(free.affine_free_components == 0);
    CHECK(free.quadratic_free_components == 1);
    CHECK_THAT(free.elastic_energy_j,WithinAbs(0.,1e-15));
    CHECK(shifted.elastic_energy_j > 0.); // The attached body still uses Q1.
}

TEST_CASE("A compatible linear strain on a free component has exact zero stress across elasticity interfaces", "[ThermoelasticSections]")
{
    LocalSectionTile tile;
    tile.index = {0,0}; tile.minimum = {0,0}; tile.maximum = {1000000000,1000000000};
    tile.slabs = {{1.,1.5,0,{square()}},{1.5,2.,1,{square()}}};
    auto first = material();
    first.eigenstrain << -.001,-.002,-.003,.0004,.0005,.0006;
    first.eigenstrain_gradient_per_m << .001,.002,.003, .004,.005,.006,
        .007,.008,.009, .010,.011,.012, .013,.014,.015, .016,.017,.018;
    auto second = first;
    second.elasticity_pa *= 3.;
    const Eigen::Vector3d point(.4,.4,1.6);
    constexpr double step = 1e-4;
    std::vector<SectionSampleRequest> requests;
    for (unsigned axis = 0; axis < 3; ++axis) {
        auto lo = point, hi = point;
        lo[axis] -= step; hi[axis] += step;
        requests.push_back({2*axis,lo,1});
        requests.push_back({2*axis+1,hi,1});
    }
    SectionSolveLimits limits;
    limits.allow_free_components = true;
    const auto solve = [&] {
        return solve_local_section_tiles({tile},0.,{first,second},Eigen::Vector3d::Zero(),{.5,.5,.5},limits,{},{},requests);
    };
    const auto field = solve();
    CHECK(field.quadratic_free_components == 1);
    CHECK(field.reference_dofs == 0);
    CHECK_THAT(field.elastic_energy_j,WithinAbs(0.,1e-15));
    REQUIRE(field.samples.size() == 6);
    Eigen::Matrix3d derivative;
    for (unsigned axis = 0; axis < 3; ++axis)
        derivative.col(axis) = (field.samples[2*axis+1].displacement_m-field.samples[2*axis].displacement_m)/(2.*step);
    Eigen::Matrix<double,6,1> recovered;
    recovered << derivative(0,0),derivative(1,1),derivative(2,2),
        derivative(0,1)+derivative(1,0),derivative(1,2)+derivative(2,1),derivative(0,2)+derivative(2,0);
    CHECK_THAT((recovered-first.eigenstrain-first.eigenstrain_gradient_per_m*point).norm(),WithinAbs(0.,1e-12));
    for (const auto& sample : field.samples) {
        CHECK_FALSE(sample.bed_connected);
        CHECK_THAT(sample.stress_pa.norm(),WithinAbs(0.,1e-15));
    }
    second.eigenstrain_gradient_per_m(0,0) += .01;
    const auto incompatible = solve();
    CHECK(incompatible.quadratic_free_components == 0);
    CHECK(incompatible.reference_dofs == 6);
    CHECK(incompatible.elastic_energy_j > 0.);
}

TEST_CASE("Graded section sampling recovers affine free contraction", "[ThermoelasticSections]")
{
    LocalSectionTile tile;
    tile.index = {0,0};
    tile.minimum = {0,0};
    tile.maximum = {1000000000,1000000000};
    tile.slabs = {{1.,2.,0,{square()}}};
    auto elastic = material();
    elastic.eigenstrain << -.001,-.002,-.003,.0001,.0002,.0003;
    SectionSolveLimits limits;
    limits.allow_free_components = true;
    const auto field = solve_local_section_tiles({tile},0.,{elastic},Eigen::Vector3d::Zero(),
        {.5,.5,.5},limits,{},{},{{1,{.3,.3,1.1},0},{2,{.3,.3,1.7},0}}, {0.,1.,1.25,2.});
    REQUIRE(field.samples.size() == 2);
    CHECK(field.affine_free_components == 1);
    CHECK_THAT(field.elastic_energy_j,WithinAbs(0.,1e-12));
    for (const auto& sample : field.samples) {
        CHECK_FALSE(sample.bed_connected);
        CHECK_THAT(sample.stress_pa.norm(),WithinAbs(0.,1e-10));
    }
}

TEST_CASE("Graded attached sections recover compatible vertical contraction at cell interfaces", "[ThermoelasticSections]")
{
    LocalSectionTile tile;
    tile.index = {0,0};
    tile.minimum = {0,0};
    tile.maximum = {1000000000,1000000000};
    tile.slabs = {{0.,1.,0,{square()}}};
    auto elastic = material();
    elastic.eigenstrain.z() = -.001;
    const std::vector<SectionSampleRequest> samples{
        {1,{.3,.3,.125},0}, {2,{.3,.3,.25},0}, {3,{.3,.3,.5},0}, {4,{.3,.3,.9},0}};
    const auto field = solve_local_section_tiles({tile},0.,{elastic},Eigen::Vector3d::Zero(),
        {.5,.5,.5},{},{},{},samples,{0.,.25,.5,1.});
    REQUIRE(field.samples.size() == samples.size());
    CHECK(field.affine_free_components == 0);
    CHECK_THAT(field.elastic_energy_j,WithinAbs(0.,1e-12));
    for (const auto& sample : field.samples) {
        CHECK(sample.bed_connected);
        CHECK_THAT(sample.stress_pa.norm(),WithinAbs(0.,1e-10));
    }
    for (std::size_t i = 0; i < field.nodes_m.size(); ++i) {
        const Eigen::Vector3d expected(0.,0.,-.001 * field.nodes_m[i].z());
        CHECK_THAT((field.displacement_m[i]-expected).norm(),WithinAbs(0.,1e-12));
        CHECK_THAT(field.generalized_bed_reaction_n[i].norm(),WithinAbs(0.,1e-10));
    }
}

TEST_CASE("Free section references preserve contraction without bed reactions", "[ThermoelasticSections]")
{
    LocalSectionTile tile;
    tile.index = {0,0};
    tile.minimum = {0,0};
    tile.maximum = {1000000000,1000000000};
    tile.slabs = {{1.,2.,0,{square()}}};
    auto elastic = material();
    elastic.eigenstrain.head<3>().setConstant(-.001);
    SectionSolveLimits limits;
    const auto solve = [&] {
        return solve_local_section_tiles({tile},0.,{elastic},Eigen::Vector3d::Zero(),
            {.5,.5,.5},limits,{}, {}, {{1,{.5,.5,1.5},0}});
    };
    CHECK_THROWS(solve());
    limits.allow_free_components = true;
    const auto field = solve();
    CHECK(field.free_components == 1);
    CHECK(field.reference_dofs == 0);
    CHECK(field.affine_free_components == 1);
    CHECK_THAT(field.bed_area_m2,WithinAbs(0.,1e-15));
    CHECK_THAT(field.elastic_energy_j,WithinAbs(0.,1e-12));
    REQUIRE(field.samples.size() == 1);
    CHECK_FALSE(field.samples[0].bed_connected);
    CHECK_THAT(field.samples[0].stress_pa.norm(),WithinAbs(0.,1e-10));
    for (std::size_t i = 0; i < field.nodes_m.size(); ++i) {
        CHECK_FALSE(field.bed_connected_nodes[i]);
        CHECK_THAT(field.generalized_bed_reaction_n[i].norm(),WithinAbs(0.,1e-12));
        const Eigen::Vector3d delta = field.nodes_m[i]-field.nodes_m.front();
        const Eigen::Vector3d displacement = field.displacement_m[i]-field.displacement_m.front();
        CHECK_THAT(delta.dot(displacement),WithinAbs(-.001*delta.squaredNorm(),1e-12));
    }
}

TEST_CASE("Free components leave attached section forces unchanged", "[ThermoelasticSections]")
{
    LocalSectionTile tile;
    tile.index = {0,0};
    tile.minimum = {0,0};
    tile.maximum = {3000000000LL,1000000000};
    tile.slabs = {{0.,1.,0,{square()}}};
    auto elastic = material();
    elastic.eigenstrain.head<3>().setConstant(-.001);
    SectionSolveLimits limits;
    limits.allow_free_components = true;
    const auto solve = [&] {
        return solve_local_section_tiles({tile},0.,{elastic},Eigen::Vector3d::Zero(),
            {4.,2.,1.},limits,{}, {}, {{1,{.5,.5,.5},0},{2,{2.5,.5,1.5},0}});
    };
    const auto attached = solve();
    auto floating = square();
    floating.translate(2000000000,0);
    tile.slabs.push_back({1.,2.,0,{floating}});
    const auto mixed = solve();
    CHECK(mixed.free_components == 1);
    REQUIRE(mixed.samples.size() == 2);
    REQUIRE(attached.samples.size() == 1);
    CHECK(mixed.samples[0].bed_connected);
    CHECK_FALSE(mixed.samples[1].bed_connected);
    CHECK_THAT(mixed.elastic_energy_j,WithinRel(attached.elastic_energy_j,1e-10));
    CHECK_THAT((mixed.samples[0].stress_pa-attached.samples[0].stress_pa).norm(),WithinAbs(0.,1e-10));
    for (std::size_t i = 0; i < attached.nodes_m.size(); ++i) {
        bool found = false;
        for (std::size_t j = 0; j < mixed.nodes_m.size(); ++j) {
            if (!mixed.bed_connected_nodes[j] || (attached.nodes_m[i]-mixed.nodes_m[j]).norm() > 1e-12) continue;
            found = true;
            CHECK_THAT((attached.generalized_bed_reaction_n[i]-mixed.generalized_bed_reaction_n[j]).norm(),WithinAbs(0.,1e-10));
        }
        CHECK(found);
    }
    Eigen::Vector3d first = Eigen::Vector3d::Zero(), second = Eigen::Vector3d::Zero();
    for (const auto& force : attached.generalized_bed_reaction_n) first += force;
    for (const auto& force : mixed.generalized_bed_reaction_n) second += force;
    CHECK_THAT((first-second).norm(),WithinAbs(0.,1e-10));
}

TEST_CASE("Small free sections retain affine shear on a larger background grid", "[ThermoelasticSections]")
{
    LocalSectionTile tile;
    tile.index = {0,0};
    tile.minimum = {180000000,120000000};
    tile.maximum = {190000000,130000000};
    const ExPolygon footprint(Polygon{{185001494,120654239},{185662510,120654239},
        {185662510,121292759},{185001494,121292759}});
    tile.slabs = {{.0098,.0114,0,{footprint}}};
    auto elastic = material();
    elastic.eigenstrain << -.001,-.002,-.003,.0001,.0002,.0003;
    SectionSolveLimits limits;
    limits.allow_free_components = true;
    const auto field = solve_local_section_tiles({tile},0.,{elastic},Eigen::Vector3d::Zero(),
        {.02,.02,.005},limits,{}, {}, {{1,{.1853,.1209,.0105},0}});
    REQUIRE(field.samples.size() == 1);
    CHECK_FALSE(field.samples[0].bed_connected);
    CHECK_THAT(field.samples[0].stress_pa.norm(),WithinAbs(0.,1e-6));
    CHECK_THAT(field.elastic_energy_j,WithinAbs(0.,1e-15));
    for (const auto& force : field.generalized_bed_reaction_n)
        CHECK_THAT(force.norm(),WithinAbs(0.,1e-15));
}

TEST_CASE("Nonuniform free contraction uses a reaction-free reference frame", "[ThermoelasticSections]")
{
    LocalSectionTile tile;
    tile.index = {0,0};
    tile.minimum = {0,0};
    tile.maximum = {1000000000,1000000000};
    tile.slabs = {{1.,1.5,0,{square()}},{1.5,2.,1,{square()}}};
    auto first = material(), second = material();
    first.eigenstrain.head<3>().setConstant(-.001);
    second.eigenstrain.head<3>().setConstant(-.002);
    SectionSolveLimits limits;
    limits.allow_free_components = true;
    const auto field = solve_local_section_tiles({tile},0.,{first,second},Eigen::Vector3d::Zero(),
        {.5,.5,.5},limits);
    CHECK(field.free_components == 1);
    CHECK(field.affine_free_components == 0);
    CHECK(field.reference_dofs == 6);
    CHECK(field.elastic_energy_j > 0.);
    CHECK(field.relative_residual < 1e-8);
    for (const auto& force : field.generalized_bed_reaction_n)
        CHECK_THAT(force.norm(),WithinAbs(0.,1e-12));
    second.eigenstrain = first.eigenstrain;
    second.elasticity_pa *= 1000.;
    const auto uniform = solve_local_section_tiles({tile},0.,{first,second},Eigen::Vector3d::Zero(),
        {.5,.5,.5},limits);
    CHECK(uniform.affine_free_components == 1);
    CHECK(uniform.reference_dofs == 0);
    CHECK_THAT(uniform.elastic_energy_j,WithinAbs(0.,1e-12));
}

TEST_CASE("Small detached components retain contraction energy under support scaling", "[ThermoelasticSections]")
{
    double reference_density=0.;
    Eigen::Matrix<double,6,1> reference_stress;
    for (const coord_t side : {coord_t(1000000),coord_t(1000)}) {
        const coord_t x=180001000,y=110001000,z=4000100;
        const ExPolygon footprint(Polygon{{x,y},{x+side,y},{x+side,y+side},{x,y+side}});
        LocalSectionTile tile;
        tile.index={0,0}; tile.minimum={0,0}; tile.maximum={1000000000,1000000000};
        tile.slabs={{double(z)*1e-9,double(z+side/2)*1e-9,0,{footprint}},
                    {double(z+side/2)*1e-9,double(z+side)*1e-9,1,{footprint}}};
        auto first=material(),second=material();
        first.eigenstrain.head<3>().setConstant(-.001);
        second.eigenstrain.head<3>().setConstant(-.002);
        SectionSolveLimits limits;
        limits.allow_free_components=true;
        const Eigen::Vector3d sample(double(x+side/2)*1e-9,double(y+side/2)*1e-9,double(z+side/4)*1e-9);
        const auto field=solve_local_section_tiles({tile},0.,{first,second},Eigen::Vector3d::Zero(),
            {.005,.005,.002},limits,{}, {},{{1,sample,0}});
        REQUIRE(field.samples.size()==1);
        CHECK(field.free_components==1);
        CHECK(field.reference_dofs==6);
        CHECK(field.elastic_energy_j>0.);
        CHECK(field.relative_residual<1e-8);
        const double density=field.elastic_energy_j/field.occupied_volume_m3;
        if (side==1000000) {
            reference_density=density;
            reference_stress=field.samples.front().stress_pa;
        } else {
            CHECK_THAT(density,WithinRel(reference_density,1e-7));
            CHECK_THAT((field.samples.front().stress_pa-reference_stress).norm(),WithinAbs(0.,1e-7));
        }
        for (const auto& reaction : field.generalized_bed_reaction_n)
            CHECK_THAT(reaction.norm(),WithinAbs(0.,0.));
    }
}

TEST_CASE("Explicit streaming budgets separate integration from bounded live storage and assembly", "[ThermoelasticSections]")
{
    const std::vector<DepositedSegment> beads{
        {{0., 0., .001}, {.002, 0., .001}, .001, .001, 0},
        {{0., 0., .002}, {.002, 0., .002}, .001, .001, 0}};
    LocalSectionStreamSettings settings;
    settings.origin_m = {-.001, -.001, 0.};
    settings.basis_spacing_m = {.001, .001, .001};
    settings.tile_size_m = {.001, .001};
    settings.maximum_height_m = .0005;
    auto elastic = material();
    elastic.eigenstrain.head<3>().setConstant(-.001);
    StreamIntegrationStatistics legacy, measured;
    const auto reference = solve_local_deposition_sections(beads, {elastic}, settings, {}, {}, nullptr, {}, &legacy);
    auto spatial_settings = settings;
    spatial_settings.load = [&](const SectionLoadQuery&) {
        return LinearThermalStrain{elastic.eigenstrain,elastic.eigenstrain_gradient_per_m};
    };
    const auto spatial = solve_local_deposition_sections(beads,{material()},spatial_settings);
    CHECK(spatial.field.uses_spatial_load);
    CHECK_THAT(spatial.field.elastic_energy_j,WithinRel(reference.field.elastic_energy_j,1e-12));
    CHECK_THAT(spatial.field.occupied_volume_m3,WithinRel(reference.field.occupied_volume_m3,1e-12));
    settings.integration = StreamIntegrationLimits{};
    settings.streaming_topology = StreamTopologyLimits{};
    const auto solve = [&] {
        return solve_local_deposition_sections(beads, {elastic}, settings, {}, {}, nullptr, {}, &measured);
    };
    const auto streamed = solve();
    CHECK(measured.integration_work == legacy.integration_work);
    CHECK(measured.processed_points == legacy.processed_points);
    CHECK(measured.peak_groups == legacy.peak_groups);
    CHECK_THAT(streamed.field.elastic_energy_j, WithinRel(reference.field.elastic_energy_j, 1e-12));
    REQUIRE(measured.integration_work > measured.assembly_work);
    settings.mechanics.max_work = measured.assembly_work;
    CHECK_NOTHROW(solve());
    settings.integration.reset();
    CHECK_THROWS_AS(solve(), ResourceLimit);
    settings.integration = StreamIntegrationLimits{};
    CHECK_NOTHROW(solve());
    const auto good = measured;
    REQUIRE(good.piece_count > 1);
    REQUIRE(good.provisional_mapping_size > 1);
    settings.integration->max_pieces = good.piece_count - 1;
    CHECK_THROWS_AS(solve(), ResourceLimit);
    settings.integration = StreamIntegrationLimits{};
    settings.integration->max_provisional_components = good.provisional_mapping_size - 1;
    CHECK_THROWS_AS(solve(), ResourceLimit);
    settings.integration = StreamIntegrationLimits{};
    CHECK_NOTHROW(solve());
    settings.integration->max_work = good.integration_work - 1;
    CHECK_THROWS_AS(solve(), ResourceLimit);
    settings.integration = StreamIntegrationLimits{};
    settings.integration->max_total_points = good.processed_points - 1;
    CHECK_THROWS_AS(solve(), ResourceLimit);
    settings.integration = StreamIntegrationLimits{};
    settings.integration->max_live_points = good.peak_input_points - 1;
    CHECK_THROWS_AS(solve(), ResourceLimit);
    settings.integration = StreamIntegrationLimits{};
    settings.integration->max_groups = good.peak_groups - 1;
    CHECK_THROWS_AS(solve(), ResourceLimit);
    settings.integration = StreamIntegrationLimits{};
    settings.integration->max_work = 600000001;
    CHECK_THROWS_AS(solve(), std::invalid_argument);
    CHECK(measured.processed_points == 0);
    settings.integration = StreamIntegrationLimits{};
    settings.integration->max_total_points = 32000001;
    CHECK_THROWS_AS(solve(), std::invalid_argument);
    settings.integration = StreamIntegrationLimits{};
    settings.integration->max_groups = 31251;
    CHECK_THROWS_AS(solve(), std::invalid_argument);
    settings.integration = StreamIntegrationLimits{};
    settings.integration->max_pieces = 250001;
    CHECK_THROWS_AS(solve(), std::invalid_argument);
    settings.integration = StreamIntegrationLimits{};
    settings.integration->max_provisional_components = 250001;
    CHECK_THROWS_AS(solve(), std::invalid_argument);
    settings.integration->max_provisional_components = 250000;
    settings.integration->max_pieces = 250000;
    settings.integration->max_total_points = 32000000;
    settings.integration->max_work = 600000000;
    CHECK_NOTHROW(solve());
    ++settings.integration->max_work;
    CHECK_THROWS_AS(solve(), std::invalid_argument);
    settings.mechanics.sparse_factor = SectionFactorLimits{};
    settings.integration->max_work = 600000000;
    settings.integration->max_pieces = 500000;
    CHECK_NOTHROW(solve());
    ++settings.integration->max_pieces;
    CHECK_THROWS_AS(solve(), std::invalid_argument);
    --settings.integration->max_pieces;
    settings.integration->max_work = 3200000000u;
    CHECK_NOTHROW(solve());
    ++settings.integration->max_work;
    CHECK_THROWS_AS(solve(), std::invalid_argument);
    --settings.integration->max_work;
    settings.integration->max_total_points = 64000000;
    settings.integration->max_groups = 125000;
    CHECK_NOTHROW(solve());
    ++settings.integration->max_total_points;
    CHECK_THROWS_AS(solve(), std::invalid_argument);
    --settings.integration->max_total_points;
    ++settings.integration->max_groups;
    CHECK_THROWS_AS(solve(), std::invalid_argument);
    settings.mechanics.sparse_factor.reset();
    settings.integration = StreamIntegrationLimits{};
    CHECK_THROWS_AS(solve_local_deposition_sections(beads, {elastic}, settings, {},
        [] { throw ResourceLimit("cancel stream integration"); }, nullptr, {}, &measured), ResourceLimit);
}

TEST_CASE("Mechanical integration borrows normalized geometry from both topology factories", "[ThermoelasticSections]")
{
    auto footprint = square();
    footprint.holes.push_back(Polygon{{200000000,200000000},{200000000,400000000},
                                      {400000000,400000000},{400000000,200000000}});
    const std::vector<DepositedSlab> repeated{{0.,1.,0,{footprint,footprint}}};
    const auto global = qualify_section_topology(repeated, 0.);
    LocalSectionTile tile;
    tile.index = {0,0};
    tile.minimum = {0,0};
    tile.maximum = {1000000000,1000000000};
    tile.slabs = repeated;
    const auto local = qualify_local_section_topology({tile}, 0.);
    REQUIRE(global.pieces().size() == 1);
    REQUIRE(local.pieces().size() == 1);
    auto elastic = material();
    elastic.eigenstrain.head<3>().setConstant(-.002);
    const std::vector<SectionSampleRequest> samples{{1,{.1,.6,.3},0},{2,{.3,.3,.3},0}};
    const auto solve = [&](const SectionTopology& topology) {
        return solve_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(), {.5,.5,.5}, {}, {}, samples);
    };
    const auto a = solve(global), b = solve(local);
    CHECK_THAT(a.occupied_volume_m3, WithinAbs(.96, 1e-12));
    CHECK_THAT(b.occupied_volume_m3, WithinAbs(a.occupied_volume_m3, 1e-12));
    CHECK_THAT(b.elastic_energy_j, WithinRel(a.elastic_energy_j, 1e-12));
    REQUIRE(a.samples.size() == 1);
    REQUIRE(b.samples.size() == 1);
    CHECK(a.unsampled_points == 1);
    CHECK(b.unsampled_points == 1);
    CHECK_THAT((a.samples[0].stress_pa - b.samples[0].stress_pa).norm(), WithinAbs(0., 1e-10));
    REQUIRE(a.nodes_m.size() == b.nodes_m.size());
    for (std::size_t i = 0; i < a.nodes_m.size(); ++i) {
        CHECK_THAT((a.displacement_m[i] - b.displacement_m[i]).norm(), WithinAbs(0., 1e-12));
        CHECK_THAT((a.generalized_bed_reaction_n[i] - b.generalized_bed_reaction_n[i]).norm(), WithinAbs(0., 1e-10));
    }
}

TEST_CASE("Equilibration preserves independently scaled anchored component fields", "[ThermoelasticSections]")
{
    auto right = square();
    right.translate(2000000000,0);
    const DepositedSlab left_slab{0.,1.,0,{square()}}, right_slab{0.,1.,1,{right}};
    auto stiff = material(), soft = material();
    stiff.eigenstrain.head<3>().setConstant(-.001);
    soft.eigenstrain = stiff.eigenstrain;
    soft.elasticity_pa *= 1e-15;
    const std::vector<SectionMaterial> materials{stiff,soft};
    const std::vector<SectionSampleRequest> samples{{1,{.5,.5,.5},0},{2,{2.5,.5,.5},1}};
    const auto solve = [&](const std::vector<DepositedSlab>& slabs) {
        return solve_deposited_sections(qualify_section_topology(slabs,0.),materials,
            Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),{},{},samples);
    };
    const auto combined = solve({left_slab,right_slab});
    const auto first = solve({left_slab}), second = solve({right_slab});
    REQUIRE(combined.samples.size() == 2);
    REQUIRE(first.samples.size() == 1);
    REQUIRE(second.samples.size() == 1);
    CHECK_THAT((combined.samples[0].displacement_m-first.samples[0].displacement_m).norm(),WithinAbs(0.,1e-12));
    CHECK_THAT((combined.samples[1].displacement_m-second.samples[0].displacement_m).norm(),WithinAbs(0.,1e-12));
    CHECK_THAT((combined.samples[1].stress_pa/1e-15-first.samples[0].stress_pa).norm(),WithinAbs(0.,1e-10));
    CHECK_THAT(second.elastic_energy_j/1e-15,WithinRel(first.elastic_energy_j,1e-10));
    CHECK_THAT(combined.elastic_energy_j,WithinRel(first.elastic_energy_j+second.elastic_energy_j,1e-10));
    CHECK(combined.relative_residual < 1e-8);
    for (std::size_t i = 0; i < combined.nodes_m.size(); ++i) {
        const auto& reference = combined.nodes_m[i].x() < 1.5 ? first : second;
        bool matched = false;
        for (std::size_t j = 0; j < reference.nodes_m.size(); ++j)
            if (combined.nodes_m[i] == reference.nodes_m[j]) {
                const double factor = combined.nodes_m[i].x() < 1.5 ? 1. : 1e15;
                CHECK_THAT((factor*(combined.generalized_bed_reaction_n[i]-reference.generalized_bed_reaction_n[j])).norm(),
                    WithinAbs(0.,1e-10));
                matched = true;
                break;
            }
        REQUIRE(matched);
    }
}

TEST_CASE("Cell material moments preserve loaded holed multimaterial fields under slab splitting", "[ThermoelasticSections]")
{
    auto left = square();
    left.holes.push_back(Polygon{{200000000,200000000},{200000000,400000000},
                                 {400000000,400000000},{400000000,200000000}});
    auto right = square();
    right.translate(1000000000, 0);
    std::vector<DepositedSlab> whole{{0.,1.,0,{left}}, {0.,1.,1,{right}}}, split;
    for (unsigned i = 0; i < 8; ++i) {
        split.push_back({i / 8., (i + 1) / 8., 0, {left}});
        split.push_back({i / 8., (i + 1) / 8., 1, {right}});
    }
    auto first = material(), second = material();
    first.eigenstrain.head<3>().setConstant(-.002);
    second.elasticity_pa *= 2.;
    second.eigenstrain.head<3>().setConstant(-.001);
    const std::vector<SectionSampleRequest> samples{
        {1,{.1,.6,.3},0}, {2,{1.4,.6,.3},1}, {3,{.3,.3,.3},0}, {4,{1.4,.6,.3},0}};
    const auto solve = [&](const std::vector<DepositedSlab>& slabs) {
        return solve_deposited_sections(qualify_section_topology(slabs, 0.), {first,second},
            Eigen::Vector3d::Zero(), {.75,.5,.5}, {}, {}, samples);
    };
    const auto a = solve(whole), b = solve(split);
    LocalSectionTile tile;
    tile.index = {0,0};
    tile.minimum = {0,0};
    tile.maximum = {2000000000,1000000000};
    tile.slabs = split;
    const auto streamed = solve_local_section_tiles({tile}, 0., {first,second}, Eigen::Vector3d::Zero(),
        {.75,.5,.5}, {}, {}, {}, samples);
    REQUIRE(a.nodes_m.size() == b.nodes_m.size());
    REQUIRE(a.samples.size() == 2);
    REQUIRE(b.samples.size() == a.samples.size());
    CHECK(a.unsampled_points == b.unsampled_points);
    CHECK_THAT(b.occupied_volume_m3, WithinRel(a.occupied_volume_m3, 1e-12));
    CHECK_THAT(b.bed_area_m2, WithinRel(a.bed_area_m2, 1e-12));
    CHECK_THAT(b.elastic_energy_j, WithinRel(a.elastic_energy_j, 1e-10));
    CHECK_THAT(streamed.elastic_energy_j, WithinRel(a.elastic_energy_j, 1e-10));
    CHECK_THAT(streamed.occupied_volume_m3, WithinRel(a.occupied_volume_m3, 1e-12));
    REQUIRE(streamed.nodes_m.size() == a.nodes_m.size());
    REQUIRE(streamed.samples.size() == a.samples.size());
    CHECK(streamed.unsampled_points == a.unsampled_points);
    for (std::size_t i = 0; i < a.nodes_m.size(); ++i) {
        CHECK_THAT((a.nodes_m[i] - b.nodes_m[i]).norm(), WithinAbs(0., 1e-14));
        CHECK_THAT((a.displacement_m[i] - b.displacement_m[i]).norm(), WithinAbs(0., 1e-12));
        CHECK_THAT((a.generalized_bed_reaction_n[i] - b.generalized_bed_reaction_n[i]).norm(), WithinAbs(0., 1e-9));
        CHECK_THAT((a.displacement_m[i] - streamed.displacement_m[i]).norm(), WithinAbs(0., 1e-12));
        CHECK_THAT((a.generalized_bed_reaction_n[i] - streamed.generalized_bed_reaction_n[i]).norm(), WithinAbs(0., 1e-9));
    }
    for (std::size_t i = 0; i < a.samples.size(); ++i) {
        CHECK(a.samples[i].id == b.samples[i].id);
        CHECK_THAT((a.samples[i].stress_pa - b.samples[i].stress_pa).norm(), WithinAbs(0., 1e-9));
        CHECK_THAT((a.samples[i].displacement_m - b.samples[i].displacement_m).norm(), WithinAbs(0., 1e-12));
        CHECK(streamed.samples[i].id == a.samples[i].id);
        CHECK_THAT((a.samples[i].stress_pa - streamed.samples[i].stress_pa).norm(), WithinAbs(0., 1e-9));
        CHECK_THAT((a.samples[i].displacement_m - streamed.samples[i].displacement_m).norm(), WithinAbs(0., 1e-12));
    }
}

TEST_CASE("Section field samples respect material holes and recover affine strain", "[ThermoelasticSections]")
{
    auto footprint = square();
    footprint.holes.push_back(Polygon{{400000000,400000000},{400000000,600000000},
                                      {600000000,600000000},{600000000,400000000}});
    const auto topology = qualify_section_topology({{0., 1., 0, {footprint}}}, 0.);
    auto elastic = material();
    elastic.eigenstrain[2] = -.001;
    const std::vector<SectionSampleRequest> requests{
        {7, {.2,.3,.4}, 0}, {8, {.5,.5,.4}, 0}, {9, {1.2,.3,.4}, 0}};
    const auto result = solve_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(),
                                                  {.5,.5,.5}, {}, {}, requests);
    REQUIRE(result.samples.size() == 1);
    CHECK(result.unsampled_points == 2);
    CHECK(result.samples.front().id == 7);
    CHECK_THAT(result.samples.front().displacement_m.z(), WithinAbs(-.0004, 1e-12));
    CHECK_THAT(result.samples.front().stress_pa.norm(), WithinAbs(0., 1e-8));
    CHECK_THAT(result.samples.front().von_mises_pa, WithinAbs(0., 1e-8));
    elastic.eigenstrain.setZero();
    elastic.eigenstrain[0] = -.001;
    const auto loaded = solve_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(),
                                                  {.5,.5,.5}, {}, {}, {{10,{.2,.3,0.},0}});
    REQUIRE(loaded.samples.size() == 1);
    // At the bonded plane all tangential derivatives of displacement vanish.
    // With nu=0, sigma_xx is exactly -E*eigenstrain_xx, independently of dz u.
    CHECK_THAT(loaded.samples.front().stress_pa[0], WithinAbs(1., 1e-10));
    CHECK_THAT(loaded.samples.front().displacement_m.norm(), WithinAbs(0., 1e-12));
}

TEST_CASE("Section basis refinement requires two stable comparisons and never publishes incomplete fields", "[ThermoelasticSections]")
{
    const auto topology = qualify_section_topology({{0., 1., 0, {square()}}}, 0.);
    auto elastic        = material();
    for (double axial_strain : {0., -.001}) {
        elastic.eigenstrain[2]  = axial_strain;
        unsigned progress_calls = 0;
        const auto study = refine_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), 3, .01, {},
                                                     [&](unsigned level) {
                                                         CHECK(level == progress_calls);
                                                         ++progress_calls;
                                                     });
        CHECK(study.state == SectionBasisStudy::State::BasisConverged);
        CHECK(study.field.has_value());
        REQUIRE(study.history.size() == 3);
        CHECK(progress_calls == 3);
        CHECK(study.history[0].nodes < study.history[1].nodes);
        CHECK(study.history[1].nodes < study.history[2].nodes);
        const auto incomplete = refine_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), 2, .01);
        CHECK(incomplete.state == SectionBasisStudy::State::Incomplete);
        CHECK_FALSE(incomplete.field.has_value());
    }
    SectionSolveLimits limits;
    limits.max_cells   = 1;
    const auto bounded = refine_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), 3, .01, limits);
    CHECK(bounded.state == SectionBasisStudy::State::ResourceLimited);
    CHECK_FALSE(bounded.field.has_value());
    REQUIRE(bounded.history.size() == 1);
    CHECK_THROWS_AS(refine_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), 6, .01),
                    std::invalid_argument);
    CHECK_THROWS_AS(refine_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), 3, .01, {}, {},
                                              [] { throw std::runtime_error("cancel"); }),
                    std::runtime_error);
    CHECK_THROWS_AS(refine_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), 3, .01, {},
                                              [](unsigned) { throw ResourceLimit("progress cancellation"); }),
                    ResourceLimit);
    bool started = false;
    CHECK_THROWS_AS(refine_deposited_sections(
                        topology, {elastic}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), 3, .01, {}, [&](unsigned) { started = true; },
                        [&] {
                            if (started)
                                throw ResourceLimit("checkpoint cancellation");
                        }),
                    ResourceLimit);
}

TEST_CASE("Basis convergence requires through-thickness refinement rather than extra lateral nodes", "[ThermoelasticSections]")
{
    const auto topology = qualify_section_topology({{0., .01, 0, {square()}}}, 0.);
    // Zero load makes all energies and reactions agree exactly, but the Z
    // basis remains a single affine cell through all three attempted levels.
    const auto unresolved = refine_deposited_sections(topology, {material()}, Eigen::Vector3d::Zero(),
                                                       Eigen::Vector3d::Ones(), 3, .01);
    REQUIRE(unresolved.history.size() == 3);
    CHECK(unresolved.history.back().nodes > unresolved.history.front().nodes);
    CHECK(unresolved.history.back().occupied_axis_cells[2] == 1);
    CHECK(unresolved.state == SectionBasisStudy::State::Incomplete);
    CHECK_FALSE(unresolved.field.has_value());
    const auto resolved = refine_deposited_sections(topology, {material()}, Eigen::Vector3d::Zero(),
                                                     {1., 1., .01}, 3, .01);
    CHECK(resolved.state == SectionBasisStudy::State::BasisConverged);
    CHECK(resolved.field.has_value());
}

TEST_CASE("Section energy roundoff scales with stiffness and bed planes use native coordinates", "[ThermoelasticSections]")
{
    const auto topology = qualify_section_topology({{1.2, 2.2, 0, {square()}}}, 1.2);
    const Eigen::Vector3d origin(0., 0., 1200000000. * 1e-9);
    for (double factor : {1e-12, 1., 1e12}) {
        auto elastic = material();
        elastic.elasticity_pa *= factor;
        elastic.eigenstrain[2] = -.001;
        const auto result      = solve_deposited_sections(topology, {elastic}, origin, Eigen::Vector3d::Ones());
        CHECK(result.elastic_energy_j >= 0.);
        CHECK(result.elastic_energy_j <= 1e-12 * factor);
    }
    CHECK_THROWS_AS(solve_deposited_sections(topology, {material()}, {0., 0., 1.2 + 1e-9}, Eigen::Vector3d::Ones()), std::invalid_argument);
}

TEST_CASE("Section mechanics reproduces a stress-free affine contraction and shear patch", "[ThermoelasticSections]")
{
    const auto topology = qualify_section_topology({{0., 1., 0, {square()}}}, 0.);
    auto elastic        = material();
    const auto zero     = solve_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones());
    CHECK_THAT(zero.occupied_volume_m3, WithinRel(1., 1e-12));
    CHECK_THAT(zero.bed_area_m2, WithinRel(1., 1e-12));
    CHECK_THAT(zero.elastic_energy_j, WithinAbs(0., 1e-15));
    elastic.eigenstrain[2] = -.001;
    elastic.eigenstrain[4] = .002;
    elastic.eigenstrain[5] = -.003;
    for (double size : {1., .5}) {
        const auto result = solve_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Constant(size));
        for (std::size_t i = 0; i < result.nodes_m.size(); ++i) {
            const Eigen::Vector3d expected(-.003 * result.nodes_m[i].z(), .002 * result.nodes_m[i].z(), -.001 * result.nodes_m[i].z());
            CHECK_THAT((result.displacement_m[i] - expected).norm(), WithinAbs(0., 1e-12));
            CHECK_THAT(result.generalized_bed_reaction_n[i].norm(), WithinAbs(0., 1e-10));
        }
        CHECK_THAT(result.elastic_energy_j, WithinAbs(0., 1e-12));
        CHECK(result.relative_residual < 1e-8);
    }
}

TEST_CASE("Section integration retains polygon holes and thin occupied portions", "[ThermoelasticSections]")
{
    auto ring = square();
    ring.holes.push_back(Polygon{{250000000, 250000000}, {250000000, 750000000}, {750000000, 750000000}, {750000000, 250000000}});
    const auto topology    = qualify_section_topology({{0., .25, 0, {ring}}}, 0.);
    auto elastic           = material();
    elastic.eigenstrain[2] = -.001;
    const auto result      = solve_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones());
    CHECK_THAT(result.occupied_volume_m3, WithinRel(.75 * .25, 1e-12));
    CHECK_THAT(result.bed_area_m2, WithinRel(.75, 1e-12));
    CHECK_THAT(result.elastic_energy_j, WithinAbs(0., 1e-12));
    for (std::size_t n = 0; n < result.nodes_m.size(); ++n)
        CHECK_THAT(result.displacement_m[n].z(), WithinAbs(-.001 * result.nodes_m[n].z(), 1e-12));
}

TEST_CASE("Translated millimetre sections retain exact volume and aligned layer topology", "[ThermoelasticSections]")
{
    ExPolygon footprint(Polygon{{0, 0}, {1000000, 0}, {1000000, 1000000}, {0, 1000000}});
    auto elastic             = material();
    elastic.eigenstrain[2]   = -.001;
    const auto base_topology = qualify_section_topology({{0., .0002, 0, {footprint}}, {.0002, .0004, 0, {footprint}}}, 0.);
    const Eigen::Vector3d spacing(.001, .001, .0002);
    const auto base = solve_deposited_sections(base_topology, {elastic}, Eigen::Vector3d::Zero(), spacing);
    REQUIRE(base.nodes_m.size() == 12);
    footprint.translate(Point(200000000, -200000000));
    const auto translated_topology = qualify_section_topology({{.2, .2002, 0, {footprint}}, {.2002, .2004, 0, {footprint}}}, .2);
    const auto translated          = solve_deposited_sections(translated_topology, {elastic}, {.2, -.2, .2}, spacing);
    REQUIRE(translated.nodes_m.size() == 12);
    CHECK_THAT(translated.occupied_volume_m3, WithinRel(4e-10, 1e-12));
    CHECK_THAT(translated.occupied_volume_m3, WithinRel(base.occupied_volume_m3, 1e-12));
    CHECK_THAT(translated.bed_area_m2, WithinRel(1e-6, 1e-12));
    for (std::size_t n = 0; n < base.nodes_m.size(); ++n)
        CHECK_THAT((base.displacement_m[n] - translated.displacement_m[n]).norm(), WithinAbs(0., 1e-15));
}

TEST_CASE("Section mechanics does not couple independent bodies sharing a background cell", "[ThermoelasticSections]")
{
    ExPolygon left(Polygon{{0, 0}, {250000000, 0}, {250000000, 1000000000}, {0, 1000000000}});
    ExPolygon right(Polygon{{750000000, 0}, {1000000000, 0}, {1000000000, 1000000000}, {750000000, 1000000000}});
    const auto topology = qualify_section_topology({{0., 1., 0, {left}}, {0., 1., 1, {right}}}, 0.);
    REQUIRE(topology.components() == 2);
    auto first = material(), second = material();
    first.eigenstrain[2]  = -.001;
    second.eigenstrain[2] = -.002;
    const auto result     = solve_deposited_sections(topology, {first, second}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones());
    REQUIRE(result.nodes_m.size() == 16);
    CHECK_THAT(result.occupied_volume_m3, WithinRel(.5, 1e-12));
    CHECK_THAT(result.elastic_energy_j, WithinAbs(0., 1e-12));
    std::size_t first_top = 0, second_top = 0;
    for (std::size_t n = 0; n < result.nodes_m.size(); ++n)
        if (result.nodes_m[n].z() > .5) {
            if (std::abs(result.displacement_m[n].z() + .001) < 1e-10)
                ++first_top;
            if (std::abs(result.displacement_m[n].z() + .002) < 1e-10)
                ++second_top;
        }
    CHECK(first_top == 4);
    CHECK(second_top == 4);
}

TEST_CASE("Section mechanics rejects unanchored singular and over-budget inputs", "[ThermoelasticSections]")
{
    CHECK_THROWS_AS(qualify_section_topology({{1., 2., 0, {square()}}}, 0.), std::invalid_argument);
    const auto topology = qualify_section_topology({{0., 1., 0, {square()}}}, 0.);
    auto elastic        = material();
    SectionSolveLimits limits;
    limits.max_free_dofs = 1;
    CHECK_THROWS_AS(solve_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), limits), ResourceLimit);
    limits           = {};
    limits.max_cells = 1;
    CHECK_THROWS_AS(solve_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Constant(.5), limits),
                    ResourceLimit);
    limits = {};
    limits.max_cells = 4000;
    limits.max_free_dofs = 12000;
    CHECK_NOTHROW(solve_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), limits));
    limits.max_cells = 4001;
    CHECK_THROWS_AS(solve_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), limits), std::invalid_argument);
    limits.max_cells = 4000;
    limits.max_free_dofs = 12001;
    CHECK_THROWS_AS(solve_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), limits), std::invalid_argument);
    CHECK_THROWS_AS(solve_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), {},
                                             [] { throw std::runtime_error("cancel"); }),
                    std::runtime_error);
    elastic.elasticity_pa(0, 0) = 0.;
    CHECK_THROWS_AS(solve_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones()), std::invalid_argument);
    // A narrow XY pillar has a poorly resolved torsional mode u=(-yz,xz,0):
    // its shear energy scales with the squared transverse radius. A thin Z
    // slab with its whole XY bed fixed does NOT provide that singular control.
    const ExPolygon pillar(Polygon{{0, 0}, {1, 0}, {1, 1}, {0, 1}});
    const auto thin = qualify_section_topology({{0., 1., 0, {pillar}}}, 0.);
    CHECK_THROWS_AS(solve_deposited_sections(thin, {material()}, Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones()), std::runtime_error);
}

TEST_CASE("Loaded section mechanics converges toward an independent tetrahedral contraction reference", "[ThermoelasticSections]")
{
    const double young = 1000., poisson = .3;
    const double mu     = young / (2. * (1. + poisson));
    const double lambda = young * poisson / ((1. + poisson) * (1. - 2. * poisson));
    auto elastic        = material();
    elastic.elasticity_pa.setZero();
    elastic.elasticity_pa.topLeftCorner<3, 3>().setConstant(lambda);
    elastic.elasticity_pa.diagonal().head<3>().array() += 2. * mu;
    elastic.elasticity_pa.diagonal().tail<3>().setConstant(mu);
    elastic.eigenstrain.head<3>().setConstant(-.001);
    const auto tetrahedral_energy = [&](int divisions) {
        std::vector<Voxel> voxels;
        for (int x = 0; x < divisions; ++x)
            for (int y = 0; y < divisions; ++y)
                for (int z = 0; z < divisions; ++z) {
                    Voxel voxel{{x, y, z}, young, poisson};
                    voxel.eigenstrain = elastic.eigenstrain;
                    voxels.push_back(voxel);
                }
        const auto mesh = mesh_voxels(voxels, Eigen::Vector3d::Zero(), Eigen::Vector3d::Constant(1. / divisions), 2000);
        return solve_bonded_bed(mesh.nodes_m, mesh.elements, mesh.bed_faces, 0.).equilibrium.elastic_energy_j;
    };
    const double fine_reference = tetrahedral_energy(8);
    REQUIRE(fine_reference > 0.);
    const auto topology    = qualify_section_topology({{0., 1., 0, {square()}}}, 0.);
    double previous_energy = 1., first_difference = 0., final_difference = 0.;
    for (int divisions : {1, 2, 4}) {
        const auto result               = solve_deposited_sections(topology, {elastic}, Eigen::Vector3d::Zero(),
                                                                   Eigen::Vector3d::Constant(1. / divisions));
        const double matching_reference = tetrahedral_energy(divisions);
        CAPTURE(divisions, result.elastic_energy_j, matching_reference, fine_reference);
        CHECK(result.elastic_energy_j > 0.);
        CHECK(result.elastic_energy_j < previous_energy);
        CHECK(matching_reference >= fine_reference * (1. - 1e-10));
        previous_energy  = result.elastic_energy_j;
        final_difference = std::abs(result.elastic_energy_j - fine_reference) / fine_reference;
        if (divisions == 1)
            first_difference = final_difference;
    }
    CHECK(final_difference < first_difference);
    CHECK(final_difference < .1);
}
