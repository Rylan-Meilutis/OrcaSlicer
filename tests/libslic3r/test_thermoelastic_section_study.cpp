#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/ThermoelasticLayers.hpp"
#include "libslic3r/ThermoelasticLocalSections.hpp"
#include "libslic3r/ThermoelasticMaterial.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticSections.hpp"
#include "libslic3r/ThermoelasticThermalStrain.hpp"
#include "libslic3r/ThermoelasticSectionStudy.hpp"
#include "libslic3r/ThermoelasticSectionTopology.hpp"
#include <algorithm>
#include <cmath>
#include <Eigen/Core>
#include <optional>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::Thermoelastic;
namespace {
SectionMaterial material()
{
    ElasticityMatrix d = ElasticityMatrix::Zero();
    d.diagonal().head<3>().setConstant(1000.);
    d.diagonal().tail<3>().setConstant(500.);
    return {d};
}
SectionStudySettings settings()
{
    SectionStudySettings result;
    result.max_section_height_m                    = .001;
    result.relative_tolerance                      = .05;
    result.sections.geometry.max_slab_group_visits = 100000;
    result.topology.max_component_pairs            = 100000;
    return result;
}
const std::vector<DepositedSegment> horizontal{{{0., 0., .001}, {.002, 0., .001}, .001, .001, 0}};
} // namespace

TEST_CASE("Horizontal nominal envelopes qualify geometry exactly after basis convergence", "[ThermoelasticSectionStudy]")
{
    const auto result = study_deposition_sections(horizontal, {material()}, {-.001, -.001, 0.}, {.001, .0005, .001}, 0., settings());
    CHECK(result.state == DepositionSectionStudy::State::NominalGeometryAndBasisConverged);
    CHECK(result.piecewise_constant_geometry);
    REQUIRE(result.field.has_value());
    REQUIRE(result.history.size() == 1);
    CHECK(result.history.front().sections == 1);
    CHECK_THAT(result.field->occupied_volume_m3, Catch::Matchers::WithinRel((.002 * .001 + 3.141592653589793 * .0005 * .0005) * .001, .01));
    auto insufficient             = settings();
    insufficient.max_basis_levels = 2;
    const auto incomplete = study_deposition_sections(horizontal, {material()}, {-.001, -.001, 0.}, {.001, .001, .001}, 0., insufficient);
    CHECK(incomplete.state == DepositionSectionStudy::State::Incomplete);
    CHECK_FALSE(incomplete.field.has_value());
}

TEST_CASE("Sloped deposition geometry requires independent section refinements", "[ThermoelasticSectionStudy]")
{
    auto beads                  = horizontal;
    beads.front().end_m.z()     = .002;
    auto limits                 = settings();
    limits.max_section_height_m = .0005;
    limits.max_geometry_levels  = 2;
    const auto incomplete       = study_deposition_sections(beads, {material()}, {-.001, -.001, 0.}, {.001, .0005, .001}, 0., limits);
    CHECK(incomplete.state == DepositionSectionStudy::State::Incomplete);
    CHECK_FALSE(incomplete.piecewise_constant_geometry);
    CHECK_FALSE(incomplete.field.has_value());
    REQUIRE(incomplete.history.size() == 2);
    CHECK(incomplete.history[1].sections > incomplete.history[0].sections);
    limits.max_geometry_levels = 3;
    const auto refined         = study_deposition_sections(beads, {material()}, {-.001, -.001, 0.}, {.001, .0005, .001}, 0., limits);
    CHECK(refined.state == DepositionSectionStudy::State::NominalGeometryAndBasisConverged);
    CHECK(refined.field.has_value());
    CHECK(refined.history.size() == 3);
}

TEST_CASE("Tapered deposition retains its variable boundary through section refinement", "[ThermoelasticSectionStudy]")
{
    const double h = .001, w = .001, length = .002;
    const std::vector<DepositedSegment> ramp{
        {{0.,0.,0.},{length,0.,h},w,0.,0,BeadSection::RectangularEnvelope,h}};
    auto limits = settings();
    limits.max_section_height_m = .0005;
    const auto result = study_deposition_sections(ramp,{material()},{-.001,-.001,0.},{.001,.0005,.001},0.,limits);
    CHECK_FALSE(result.piecewise_constant_geometry);
    REQUIRE(result.state == DepositionSectionStudy::State::NominalGeometryAndBasisConverged);
    REQUIRE(result.field);
    CHECK(result.history.size() >= 3);
    const double volume = length*w*h/2. + std::acos(-1.)*w*w*h/8.;
    CHECK_THAT(result.field->occupied_volume_m3,Catch::Matchers::WithinRel(volume,.01));
}

TEST_CASE("Unconverged contraction studies do not qualify from exact geometry alone", "[ThermoelasticSectionStudy]")
{
    auto loaded = material();
    loaded.eigenstrain.head<3>().setConstant(-.01);
    auto limits = settings();
    limits.relative_tolerance = .2;
    limits.max_basis_levels = 4;
    limits.mechanics.max_cells = 2000;
    const auto result = study_deposition_sections(horizontal, {loaded}, {-.001, -.001, 0.},
                                                   {.001, .001, .001}, 0., limits);
    // This supported bead has non-affine, bed-constrained contraction. Its
    // nominal geometry is exact, but four basis levels do not stabilize the
    // energy and tensile reaction at this tolerance. Do not accept a zero-load
    // geometry shortcut as evidence of loaded mechanical convergence.
    CHECK(result.piecewise_constant_geometry);
    CHECK(result.state == DepositionSectionStudy::State::Incomplete);
    CHECK_FALSE(result.field.has_value());
    REQUIRE(result.last_basis_history.size() == 4);
    for (const auto &sample : result.last_basis_history) {
        INFO("spacing=" << sample.spacing_m.transpose() << " energy=" << sample.elastic_energy_j
             << " tensile=" << sample.tensile_generalized_bed_force_n << " nodes=" << sample.nodes);
        CHECK(sample.relative_residual < 1e-8);
    }
    CHECK(result.last_basis_history.back().elastic_energy_j > 0.);
    CHECK(result.last_basis_history.back().nodes > result.last_basis_history.front().nodes);
}

TEST_CASE("Aspect-aware starting grids resolve loaded contraction within the existing limits", "[ThermoelasticSectionStudy]")
{
    auto loaded = material();
    loaded.eigenstrain.head<3>().setConstant(-.01);
    auto limits = settings();
    limits.relative_tolerance = .2;
    limits.mechanics.max_cells = 2000;
    const auto result = study_deposition_sections(horizontal, {loaded}, {-.001, -.001, 0.},
                                                   {.0005, .00025, .001}, 0., limits);
    for (const auto &sample : result.last_basis_history) {
        INFO("spacing=" << sample.spacing_m.transpose() << " energy=" << sample.elastic_energy_j
             << " tensile=" << sample.tensile_generalized_bed_force_n << " nodes=" << sample.nodes);
        CHECK(sample.relative_residual < 1e-8);
    }
    REQUIRE(result.state == DepositionSectionStudy::State::NominalGeometryAndBasisConverged);
    REQUIRE(result.field.has_value());
    const auto topology = qualify_section_topology(layered_deposition(horizontal, 0.).slabs, 0.);
    const auto reference = solve_deposited_sections(topology, {loaded}, {-.001, -.001, 0.},
                                                   {.000125, .000125, .000125}, limits.mechanics);
    double reference_force = 0.;
    for (const auto &reaction : reference.generalized_bed_reaction_n)
        reference_force += std::max(0., -reaction.z());
    CHECK_THAT(result.field->elastic_energy_j, Catch::Matchers::WithinRel(reference.elastic_energy_j, .05));
    CHECK_THAT(result.history.back().tensile_generalized_bed_force_n, Catch::Matchers::WithinRel(reference_force, .05));
}

TEST_CASE("Diagnostic seed selection still requires full-axis basis convergence", "[ThermoelasticSectionStudy]")
{
    auto loaded = material();
    loaded.eigenstrain.head<3>().setConstant(-.01);
    auto limits = settings();
    limits.relative_tolerance = .2;
    limits.mechanics.max_cells = 2000;
    limits.max_initial_basis_probes = 3;
    unsigned calls = 0;
    const auto result = study_deposition_sections(horizontal, {loaded}, {-.001, -.001, 0.},
        {.001, .001, .001}, 0., limits, {}, [&](unsigned geometry, unsigned run) {
            CHECK(geometry == 0);
            CHECK(run == calls);
            ++calls;
        });
    REQUIRE(result.last_initial_basis_probes.size() == 3);
    for (const auto &probe : result.last_initial_basis_probes) {
        REQUIRE(probe.recommended_spacing_m.has_value());
        INFO("selected spacing=" << probe.recommended_spacing_m->transpose());
        CHECK(probe.state == SectionDirectionalProbe::State::Complete);
    }
    REQUIRE(result.state == DepositionSectionStudy::State::NominalGeometryAndBasisConverged);
    REQUIRE(result.field.has_value());
    CHECK(result.last_basis_history.size() == 3);
    CHECK(calls == 15);
    CHECK(result.reserved_work <= limits.max_total_work);
    limits.max_basis_levels = 2;
    const auto incomplete = study_deposition_sections(horizontal, {loaded}, {-.001, -.001, 0.},
                                                       {.001, .001, .001}, 0., limits);
    CHECK(incomplete.state == DepositionSectionStudy::State::Incomplete);
    CHECK_FALSE(incomplete.field.has_value());
    limits.max_basis_levels = 3;
    limits.max_total_work = 35000000;
    const auto bounded = study_deposition_sections(horizontal, {loaded}, {-.001, -.001, 0.},
                                                    {.001, .001, .001}, 0., limits);
    CHECK(bounded.state == DepositionSectionStudy::State::ResourceLimited);
    CHECK_FALSE(bounded.field.has_value());
    CHECK_THROWS_AS(study_deposition_sections(horizontal, {loaded}, {-.001, -.001, 0.},
        {.001, .001, .001}, 0., limits, {}, [](unsigned, unsigned) { throw ResourceLimit("cancel probe"); }), ResourceLimit);
    limits.max_initial_basis_probes = 4;
    CHECK_THROWS_AS(study_deposition_sections(horizontal, {loaded}, {-.001, -.001, 0.},
        {.001, .001, .001}, 0., limits), std::invalid_argument);
}

TEST_CASE("Directional basis refinement reduces constrained contraction energy", "[ThermoelasticSectionStudy]")
{
    const auto layers = layered_deposition(horizontal, 0.);
    const auto topology = qualify_section_topology(layers.slabs, 0.);
    auto loaded = material();
    loaded.eigenstrain.head<3>().setConstant(-.01);
    SectionSolveLimits limits;
    limits.max_cells = 2000;
    for (double xy : {.00025, .000125}) {
        double previous_energy = 1.;
        for (double z : {.001, .0005, .00025, .000125}) {
            const auto result = solve_deposited_sections(topology, {loaded}, {-.001, -.001, 0.}, {xy, xy, z}, limits);
            double tensile = 0.;
            for (const auto &reaction : result.generalized_bed_reaction_n)
                tensile += std::max(0., -reaction.z());
            INFO("xy=" << xy << " z=" << z << " energy=" << result.elastic_energy_j
                 << " tensile=" << tensile << " nodes=" << result.nodes_m.size());
            CHECK(result.elastic_energy_j <= previous_energy * (1. + 1e-10));
            CHECK(result.relative_residual < 1e-8);
            previous_energy = result.elastic_energy_j;
        }
    }
}

TEST_CASE("Deposition studies preserve ownership policy and reject missing materials", "[ThermoelasticSectionStudy]")
{
    auto beads = horizontal;
    beads.push_back(horizontal.front());
    beads.back().material = 1;
    CHECK_THROWS_AS(study_deposition_sections(beads, {material(), material()}, {-.001, -.001, 0.}, {.001, .001, .001}, 0., settings()),
                    LayerMaterialOverlap);
    CHECK_THROWS_AS(study_deposition_sections(beads, {material()}, {-.001, -.001, 0.}, {.001, .001, .001}, 0., settings()),
                    std::invalid_argument);
    auto limits                 = settings();
    limits.max_section_height_m = .0005;
    const auto owner            = [](double z) -> std::optional<SectionOwnership> {
        return SectionOwnership{z,
                                {{0,
                                  {ExPolygon(Polygon{{-1000000, -1000000}, {3000000, -1000000}, {3000000, 1000000}, {-1000000, 1000000}})}}}};
    };
    const auto resolved = study_deposition_sections(beads, {material(), material()}, {-.001, -.001, 0.}, {.001, .0005, .001}, 0., limits,
                                                    owner);
    CHECK_FALSE(resolved.piecewise_constant_geometry);
    CHECK(resolved.state == DepositionSectionStudy::State::NominalGeometryAndBasisConverged);
    CHECK(resolved.field.has_value());
}

TEST_CASE("Deferred partition cancellation propagates through global and local studies", "[ThermoelasticSectionStudy]")
{
    auto beads = horizontal;
    beads.push_back(horizontal.front());
    beads.back().material = 1;
    for (const bool local : {false, true}) {
        auto limits = settings();
        if (local) limits.local_tile_size_m = Eigen::Vector2d(.01, .01);
        unsigned calls = 0;
        const auto owner = [&](double z) -> std::optional<SectionOwnership> {
            SectionOwnership result{z, {}};
            result.resolve_partition = [&](const SampledSectionConflict&) -> std::shared_ptr<const AuthenticatedSectionPartition> {
                ++calls;
                throw ResourceLimit("cancel deferred partition");
            };
            return result;
        };
        CHECK_THROWS_AS(study_deposition_sections(beads, {material(), material()},
            {-.001, -.001, 0.}, {.001, .001, .001}, 0., limits, owner), ResourceLimit);
        CHECK(calls == 1);
    }
}

TEST_CASE("Directional probes identify refinement demand without publishing a field", "[ThermoelasticSectionStudy]")
{
    const auto topology = qualify_section_topology(layered_deposition(horizontal, 0.).slabs, 0.);
    auto loaded = material();
    loaded.eigenstrain.head<3>().setConstant(-.01);
    const Eigen::Vector3d spacing(.00025, .00025, .0005);
    unsigned calls = 0;
    const auto probe = probe_section_refinement(topology, {loaded}, {-.001, -.001, 0.}, spacing, {},
        [&](unsigned run) { CHECK(run == calls); ++calls; });
    REQUIRE(probe.state == SectionDirectionalProbe::State::Complete);
    REQUIRE(probe.reference.has_value());
    REQUIRE(probe.recommended_spacing_m.has_value());
    CHECK(calls == 4);
    CHECK_THAT((*probe.recommended_spacing_m - probe.directions[1]->spacing_m).norm(),
               Catch::Matchers::WithinAbs(0., 1e-15));
    for (unsigned axis = 0; axis < 3; ++axis) {
        REQUIRE(probe.directions[axis].has_value());
        CHECK(probe.directions[axis]->nodes > probe.reference->nodes);
        CHECK(probe.energy_change[axis] >= 0.);
        CHECK(probe.force_change[axis] >= 0.);
        INFO("axis=" << axis << " energy change=" << probe.energy_change[axis]
             << " force change=" << probe.force_change[axis]);
        CHECK(probe.directions[axis]->relative_residual < 1e-8);
    }
    auto limited = SectionSolveLimits{};
    limited.max_cells = 1;
    const auto rejected = probe_section_refinement(topology, {loaded}, {-.001, -.001, 0.}, spacing, limited);
    CHECK(rejected.state == SectionDirectionalProbe::State::ResourceLimited);
    CHECK_FALSE(rejected.recommended_spacing_m.has_value());
    limited = SectionSolveLimits{};
    limited.max_nodes = probe.reference->nodes;
    const auto partial = probe_section_refinement(topology, {loaded}, {-.001, -.001, 0.}, spacing, limited);
    CHECK(partial.reference.has_value());
    CHECK(partial.state == SectionDirectionalProbe::State::ResourceLimited);
    CHECK_FALSE(partial.recommended_spacing_m.has_value());
    const auto unloaded = probe_section_refinement(topology, {material()}, {-.001, -.001, 0.}, spacing);
    CHECK(unloaded.state == SectionDirectionalProbe::State::Complete);
    CHECK_FALSE(unloaded.recommended_spacing_m.has_value());
    CHECK_THROWS_AS(probe_section_refinement(topology, {loaded}, {-.001, -.001, 0.}, spacing, {},
        [](unsigned) { throw ResourceLimit("cancel progress"); }), ResourceLimit);
    CHECK_THROWS_AS(probe_section_refinement(topology, {loaded}, {-.001, -.001, 0.}, spacing, {}, {},
        [] { throw ResourceLimit("cancel checkpoint"); }), ResourceLimit);
}

TEST_CASE("Streamed directional selection matches retained geometry and propagates cancellation", "[ThermoelasticSectionStudy]")
{
    LocalSectionStreamSettings stream;
    stream.origin_m = {-.001,-.001,0.};
    stream.basis_spacing_m = {.00025,.00025,.0005};
    stream.tile_size_m = {.004,.004};
    stream.maximum_height_m = .001;
    auto loaded = material();
    loaded.eigenstrain.head<3>().setConstant(-.01);
    const auto topology = qualify_section_topology(layered_deposition(horizontal,0.).slabs,0.);
    const auto retained = probe_section_refinement(topology,{loaded},stream.origin_m,stream.basis_spacing_m);
    unsigned calls = 0;
    const auto probe = probe_local_deposition_refinement(horizontal,{loaded},stream,{},
        [&](unsigned run) { CHECK(run == calls); ++calls; });
    REQUIRE(probe.state == SectionDirectionalProbe::State::Complete);
    REQUIRE(probe.reference.has_value());
    REQUIRE(probe.recommended_spacing_m.has_value());
    REQUIRE(retained.recommended_spacing_m.has_value());
    CHECK(calls == 4);
    CHECK_THAT((*probe.recommended_spacing_m-*retained.recommended_spacing_m).norm(),Catch::Matchers::WithinAbs(0.,1e-15));
    for (unsigned axis = 0; axis < 3; ++axis) {
        REQUIRE(probe.directions[axis].has_value());
        CHECK_THAT(probe.directions[axis]->elastic_energy_j,
            Catch::Matchers::WithinRel(retained.directions[axis]->elastic_energy_j,1e-10));
        CHECK(probe.directions[axis]->elastic_energy_j <= probe.reference->elastic_energy_j * (1.+1e-9));
        CHECK_THAT(probe.directions[axis]->tensile_generalized_bed_force_n,
            Catch::Matchers::WithinRel(retained.directions[axis]->tensile_generalized_bed_force_n,1e-10));
    }
    auto limited = stream;
    limited.mechanics.max_nodes = probe.reference->nodes;
    const auto partial = probe_local_deposition_refinement(horizontal,{loaded},limited);
    CHECK(partial.state == SectionDirectionalProbe::State::ResourceLimited);
    CHECK(partial.reference.has_value());
    CHECK_FALSE(partial.recommended_spacing_m.has_value());
    CHECK_THROWS_AS(probe_local_deposition_refinement(horizontal,{loaded},stream,{},
        [](unsigned) { throw ResourceLimit("cancel progress"); }),ResourceLimit);
    CHECK_THROWS_AS(probe_local_deposition_refinement(horizontal,{loaded},stream,{}, {},
        [] { throw ResourceLimit("cancel checkpoint"); }),ResourceLimit);
    CHECK_THROWS_AS(probe_local_deposition_refinement(horizontal,{loaded},stream,
        [](double) -> std::optional<SectionOwnership> { throw ResourceLimit("cancel owner"); }),ResourceLimit);
    auto overlap = horizontal;
    overlap.push_back(horizontal.front());
    overlap.back().material = 1;
    const auto owner = [](double z) -> std::optional<SectionOwnership> {
        SectionOwnership ownership{z,{}};
        ownership.resolve_partition = [](const SampledSectionConflict&) -> std::shared_ptr<const AuthenticatedSectionPartition> {
            throw ResourceLimit("cancel deferred owner");
        };
        return ownership;
    };
    CHECK_THROWS_AS(probe_local_deposition_refinement(overlap,{loaded,loaded},stream,owner),ResourceLimit);
    stream.load = [](const SectionLoadQuery&) -> LinearThermalStrain { throw ResourceLimit("cancel thermal load"); };
    CHECK_THROWS_AS(probe_local_deposition_refinement(horizontal,{loaded},stream),ResourceLimit);
}

TEST_CASE("Streaming graded planes preserve geometry without pretending to refine a regular grid", "[ThermoelasticSectionStudy]")
{
    LocalSectionStreamSettings stream;
    stream.origin_m = {-.001,-.001,0.};
    stream.basis_spacing_m = {.001,.001,.001};
    stream.tile_size_m = {.004,.004};
    stream.maximum_height_m = .001;
    stream.basis_z_planes_m = {0.,.00025,.0005,.001};
    auto loaded = material();
    loaded.eigenstrain.head<3>().setConstant(-.01);
    const auto graded = solve_local_deposition_sections(horizontal,{loaded},stream);
    CHECK_THROWS_AS(probe_local_deposition_refinement(horizontal,{loaded},stream),std::invalid_argument);
    stream.basis_z_planes_m.clear();
    stream.basis_spacing_m.z() = .00025;
    const auto fine = solve_local_deposition_sections(horizontal,{loaded},stream);
    CHECK_THAT(graded.field.occupied_volume_m3,Catch::Matchers::WithinRel(fine.field.occupied_volume_m3,1e-12));
    CHECK(graded.field.nodes_m.size() < fine.field.nodes_m.size());
    CHECK(fine.field.elastic_energy_j <= graded.field.elastic_energy_j * (1.+1e-9));
}

TEST_CASE("Section studies release unused per-plane ownership reservations", "[ThermoelasticSectionStudy]")
{
    auto limits = settings();
    limits.max_section_height_m = .00001;
    limits.sections.geometry.max_slab_group_visits = 2000000;
    limits.max_total_work = 30000000;
    const auto result = study_deposition_sections(horizontal, {material()}, {-.001, -.001, 0.},
                                                   {.001, .0005, .001}, 0., limits);
    REQUIRE(result.state == DepositionSectionStudy::State::NominalGeometryAndBasisConverged);
    CHECK(result.section_count == 100);
    CHECK(result.reserved_work < limits.max_total_work);
    CHECK(result.field.has_value());
}

TEST_CASE("Deposition study budgets and callback cancellation never publish partial fields", "[ThermoelasticSectionStudy]")
{
    auto limits           = settings();
    limits.max_total_work = 1;
    const auto bounded    = study_deposition_sections(horizontal, {material()}, {-.001, -.001, 0.}, {.001, .001, .001}, 0., limits);
    CHECK(bounded.state == DepositionSectionStudy::State::ResourceLimited);
    CHECK_FALSE(bounded.field.has_value());
    limits                           = settings();
    limits.max_total_boundary_points = 1;
    const auto points = study_deposition_sections(horizontal, {material()}, {-.001, -.001, 0.}, {.001, .001, .001}, 0., limits);
    CHECK(points.state == DepositionSectionStudy::State::ResourceLimited);
    CHECK_FALSE(points.field.has_value());
    CHECK_THROWS_AS(study_deposition_sections(horizontal, {material()}, {-.001, -.001, 0.}, {.001, .001, .001}, 0., settings(), {},
                                              [](unsigned, unsigned) { throw ResourceLimit("cancel progress"); }),
                    ResourceLimit);
    CHECK_THROWS_AS(study_deposition_sections(horizontal, {material()}, {-.001, -.001, 0.}, {.001, .001, .001}, 0., settings(), {}, {},
                                              [] { throw ResourceLimit("cancel checkpoint"); }),
                    ResourceLimit);
}
