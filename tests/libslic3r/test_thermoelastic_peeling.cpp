#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include "libslic3r/ThermoelasticPeeling.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticSolver.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <stdexcept>
#include <vector>

using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

TEST_CASE("Finite-patch bed release compares elastic energy with explicit fracture work", "[ThermoelasticPeeling][BedReleaseEvolution]")
{
    std::vector<Voxel> voxels;
    for (int x = 0; x < 8; ++x)
        for (int z = 0; z < 2; ++z) {
            Voxel voxel{{x,0,z}, 2.9e9, 0.};
            voxel.eigenstrain[0] = -.005 * (z + .5) / 2.;
            voxels.push_back(voxel);
        }
    const auto mesh = mesh_voxels(voxels, Eigen::Vector3d::Zero(), {.001,.001,.0005}, 1000);
    std::vector<std::size_t> all, broad, narrow, outer, middle;
    for (std::size_t i = 0; i < mesh.bed_faces.size(); ++i) {
        all.push_back(i);
        double x = 0.;
        for (const auto n : mesh.bed_faces[i]) x += mesh.nodes_m[n].x() / 3.;
        if (x < .004) broad.push_back(i); else outer.push_back(i);
        if (x < .001) narrow.push_back(i); else if (x < .004) middle.push_back(i);
    }
    const auto reference = bed_release_demands(mesh.nodes_m, mesh.elements, mesh.bed_faces, {broad}, 0., {}, BedContactLimits{});
    REQUIRE(reference.size() == 1);
    const double critical = reference.front().incremental_energy_release_j_m2;
    REQUIRE(critical > 0.);
    // These resistances deliberately bracket the independently evaluated
    // finite-increment demand. They are not claimed ASA/PEI properties.
    const auto strong = evolve_bed_release(mesh.nodes_m, mesh.elements, mesh.bed_faces, all, {{outer, 2. * critical}}, 0.);
    CHECK(strong.steps.empty());
    CHECK(strong.retained_faces == all);
    CHECK_THAT(strong.equilibrium.maximum_lift_m, WithinAbs(0., 1e-14));
    const auto weak = evolve_bed_release(mesh.nodes_m, mesh.elements, mesh.bed_faces, all, {{outer, .5 * critical}}, 0.);
    REQUIRE(weak.steps.size() == 1);
    CHECK(weak.steps.front().patch == 0);
    CHECK(weak.retained_faces == broad);
    CHECK_THAT(weak.steps.front().released_elastic_energy_j,
        WithinRel(2. * weak.steps.front().fracture_work_j, 1e-8));
    CHECK_THAT(weak.equilibrium.maximum_lift_m, WithinAbs(reference.front().maximum_lift_m, 1e-12));
    CHECK(weak.equilibrium.maximum_lift_m > 0.);
    const std::vector<BedReleasePatch> patches{{outer, 0.}, {middle, 0.}};
    const auto progressive = evolve_bed_release(mesh.nodes_m, mesh.elements, mesh.bed_faces, all, patches, 0.);
    REQUIRE(progressive.steps.size() == 2);
    CHECK(progressive.steps[0].patch != progressive.steps[1].patch);
    CHECK(progressive.retained_faces == narrow);
    CHECK(progressive.equilibrium.nonpenetrating);
    for (const auto &step : progressive.steps) {
        CHECK(step.released_elastic_energy_j > step.fracture_work_j);
        CHECK(step.newly_released_area_m2 > 0.);
    }
    const auto final = solve_partially_bonded_contact(mesh.nodes_m, mesh.elements, mesh.bed_faces, narrow, 0.);
    CHECK_THAT(progressive.equilibrium.bonded.equilibrium.elastic_energy_j,
        WithinAbs(final.bed.bonded.equilibrium.elastic_energy_j, 1e-14));
    auto stiffer_elements = mesh.elements;
    for (auto &element : stiffer_elements) element.young_pa *= 2.;
    const auto original_strength = evolve_bed_release(mesh.nodes_m, mesh.elements, mesh.bed_faces,
        all, {{outer, 1.5 * critical}}, 0.);
    CHECK(original_strength.steps.empty());
    const auto stiffer = evolve_bed_release(mesh.nodes_m, stiffer_elements, mesh.bed_faces,
        all, {{outer, 1.5 * critical}}, 0.);
    CHECK(stiffer.steps.size() == 1);
    auto relaxed_elements = mesh.elements;
    for (auto &element : relaxed_elements) element.eigenstrain.setZero();
    const auto relaxed = evolve_bed_release(mesh.nodes_m, relaxed_elements, mesh.bed_faces, all, patches, 0.);
    CHECK(relaxed.steps.empty());
    // Unloading may close the gap, but must not heal adhesion. Continuing with
    // the original patch IDs must not pay the same fracture work again.
    const auto unloaded = evolve_bed_release(mesh.nodes_m, relaxed_elements, mesh.bed_faces,
        progressive.retained_faces, patches, 0.);
    CHECK(unloaded.steps.empty());
    CHECK(unloaded.retained_faces == narrow);
    CHECK_THAT(unloaded.equilibrium.maximum_lift_m, WithinAbs(0., 1e-14));
    const auto reloaded = evolve_bed_release(mesh.nodes_m, mesh.elements, mesh.bed_faces,
        unloaded.retained_faces, patches, 0.);
    CHECK(reloaded.steps.empty());
    CHECK(reloaded.retained_faces == narrow);
    CHECK_THAT(reloaded.equilibrium.maximum_lift_m, WithinRel(progressive.equilibrium.maximum_lift_m, 1e-8));
    CHECK(reloaded.equilibrium_solves == final.iterations);
    auto partly_bonded = all;
    partly_bonded.erase(partly_bonded.begin());
    CHECK_THROWS_AS(evolve_bed_release(mesh.nodes_m, mesh.elements, mesh.bed_faces,
        partly_bonded, {{all, critical}}, 0.), std::invalid_argument);
    CHECK_THROWS_AS(evolve_bed_release(mesh.nodes_m, mesh.elements, mesh.bed_faces,
        all, {{all, 0.}}, 0.), std::runtime_error);
    BedReleaseEvolutionLimits limited;
    limited.max_equilibrium_solves = 1;
    CHECK_THROWS_AS(evolve_bed_release(mesh.nodes_m, mesh.elements, mesh.bed_faces, all, patches, 0., limited), ResourceLimit);
    limited = {};
    limited.max_steps = 1;
    CHECK_THROWS_AS(evolve_bed_release(mesh.nodes_m, mesh.elements, mesh.bed_faces, all, patches, 0., limited), ResourceLimit);
    CHECK_THROWS_AS(evolve_bed_release(mesh.nodes_m, mesh.elements, mesh.bed_faces, all,
        {{outer, critical}, {outer, critical}}, 0.), std::invalid_argument);
    CHECK_THROWS_AS(evolve_bed_release(mesh.nodes_m, mesh.elements, mesh.bed_faces, all,
        {{outer, -1.}}, 0.), std::invalid_argument);
    CHECK_THROWS_AS(evolve_bed_release(mesh.nodes_m, mesh.elements, mesh.bed_faces, all,
        {}, 0.), std::invalid_argument);
    CHECK_THROWS_AS(evolve_bed_release(mesh.nodes_m, mesh.elements, mesh.bed_faces, all, patches, 0., {},
        [] { throw std::runtime_error("cancelled release evolution"); }), std::runtime_error);
}

TEST_CASE("An attached brim does not prevent released case material from rising into later layers", "[ThermoelasticPeeling][BrimLift]")
{
    // Mechanism experiment, not a calibrated ASA prediction: a 20 x 2 x 2 mm
    // case strip joins a 0.2 mm thick, fully attached brim. Only the case bed
    // contact is released. Keep the case and imposed strain identical when
    // increasing brim width; never count attached brim area as zero lift risk.
    const int brim_mm = GENERATE(4, 8);
    const int refinement = GENERATE(1, 2);
    CAPTURE(brim_mm, refinement);
    const double dx = .001 / refinement, dz = .0002 / refinement;
    const int brim_cells = brim_mm * refinement;
    std::vector<Voxel> voxels;
    for (int x = 0; x < brim_cells + 20 * refinement; ++x)
        for (int y = 0; y < 2 * refinement; ++y)
            for (int z = 0; z < (x < brim_cells ? refinement : 10 * refinement); ++z) {
                Voxel voxel{{x,y,z}, 2.9e9, 0.};
                if (x >= brim_cells)
                    voxel.eigenstrain[0] = -.005 * ((z + .5) * dz) / .002;
                voxels.push_back(voxel);
            }
    auto mesh = mesh_voxels(voxels, Eigen::Vector3d::Zero(), {dx,dx,dz}, 12000);
    std::vector<std::size_t> retained;
    for (std::size_t i = 0; i < mesh.bed_faces.size(); ++i) {
        double center_x = 0.;
        for (const auto n : mesh.bed_faces[i]) center_x += mesh.nodes_m[n].x() / 3.;
        if (center_x < brim_mm * .001) retained.push_back(i);
    }
    const auto lifted = solve_partially_bonded_bed(mesh.nodes_m, mesh.elements, mesh.bed_faces, retained, 0.);
    INFO("Case lift in mm: " << lifted.maximum_lift_m * 1000.);
    CHECK(lifted.nonpenetrating);
    CHECK(lifted.maximum_lift_m > .0002);
    for (std::size_t n = 0; n < mesh.nodes_m.size(); ++n)
        if (lifted.bonded.contact_area_m2[n] > 0.)
            CHECK_THAT(lifted.footprint_gap_m[n], WithinAbs(0., 1e-14));
    // Changing the driving strain, not the attachment area: the linear model
    // must halve lift and quarter stored energy. This is a sensitivity result,
    // not a claim that a fan setting halves actual locked-in strain.
    for (auto &element : mesh.elements) element.eigenstrain *= .5;
    const auto reduced = solve_partially_bonded_bed(mesh.nodes_m, mesh.elements, mesh.bed_faces, retained, 0.);
    CHECK_THAT(reduced.maximum_lift_m, WithinRel(.5 * lifted.maximum_lift_m, 1e-8));
    CHECK_THAT(reduced.bonded.equilibrium.elastic_energy_j,
        WithinRel(.25 * lifted.bonded.equilibrium.elastic_energy_j, 1e-8));
}

TEST_CASE("Prescribed release reports energy demand without assuming adhesion capacity", "[ThermoelasticPeeling]")
{
    std::vector<Voxel> voxels;
    for (int x = 0; x < 8; ++x)
        for (int y = 0; y < 2; ++y)
            for (int z = 0; z < 2; ++z) {
                Voxel voxel{{x,y,z}, 2.9e9, 0.};
                voxel.eigenstrain[0] = -.005 * (z + .5) / 2.;
                voxels.push_back(voxel);
            }
    auto mesh = mesh_voxels(voxels, Eigen::Vector3d::Zero(), {.001,.001,.0005}, 1000);
    std::vector<std::size_t> broad, narrow;
    for (std::size_t i = 0; i < mesh.bed_faces.size(); ++i) {
        double x = 0.;
        for (const auto n : mesh.bed_faces[i]) x += mesh.nodes_m[n].x() / 3.;
        if (x < .004) broad.push_back(i);
        if (x < .001) narrow.push_back(i);
    }
    const auto result = bed_release_demands(mesh.nodes_m, mesh.elements, mesh.bed_faces, {broad,narrow}, 0.);
    REQUIRE(result.size() == 2);
    CHECK(result[1].released_area_m2 > result[0].released_area_m2);
    CHECK(result[0].incremental_energy_release_j_m2 > 0.);
    CHECK(result[1].incremental_energy_release_j_m2 > 0.);
    CHECK(result[1].maximum_lift_m > 0.);
    const double delta_area = result[1].released_area_m2 - result[0].released_area_m2;
    CHECK_THAT(result[1].incremental_energy_release_j_m2 * delta_area,
        WithinRel(result[0].elastic_energy_j - result[1].elastic_energy_j, 1e-10));
    // Linear elasticity: doubling imposed strain quadruples driving energy.
    for (auto &element : mesh.elements) element.eigenstrain *= 2.;
    const auto doubled = bed_release_demands(mesh.nodes_m, mesh.elements, mesh.bed_faces, {broad,narrow}, 0.);
    CHECK_THAT(doubled[1].incremental_energy_release_j_m2,
        WithinRel(4. * result[1].incremental_energy_release_j_m2, 1e-8));
    CHECK_THROWS_AS(bed_release_demands(mesh.nodes_m, mesh.elements, mesh.bed_faces, {narrow,broad}, 0.), std::invalid_argument);
    CHECK_THROWS_AS(bed_release_demands(mesh.nodes_m, mesh.elements, mesh.bed_faces, {broad,broad}, 0.), std::invalid_argument);
    CHECK_THROWS_AS(bed_release_demands(mesh.nodes_m, mesh.elements, mesh.bed_faces, {{}}, 0.), std::invalid_argument);
    for (auto &element : mesh.elements) element.eigenstrain *= -1.;
    CHECK_THROWS_AS(bed_release_demands(mesh.nodes_m, mesh.elements, mesh.bed_faces, {narrow}, 0.), std::runtime_error);
    const auto contacting = bed_release_demands(mesh.nodes_m, mesh.elements, mesh.bed_faces, {broad,narrow},
        0., {}, BedContactLimits{});
    REQUIRE(contacting.size() == 2);
    CHECK(contacting[0].incremental_energy_release_j_m2 >= 0.);
    CHECK(contacting[1].incremental_energy_release_j_m2 >= 0.);
    CHECK(contacting[1].elastic_energy_j <= contacting[0].elastic_energy_j);
    BedContactLimits exhausted;
    exhausted.max_iterations = 1;
    CHECK_THROWS_AS(bed_release_demands(mesh.nodes_m, mesh.elements, mesh.bed_faces, {broad,narrow},
        0., {}, exhausted), std::runtime_error);
    for (auto &element : mesh.elements) element.eigenstrain.setZero();
    const auto zero = bed_release_demands(mesh.nodes_m, mesh.elements, mesh.bed_faces, {broad,narrow}, 0.);
    CHECK_THAT(zero.back().incremental_energy_release_j_m2, WithinAbs(0., 1e-15));
}
