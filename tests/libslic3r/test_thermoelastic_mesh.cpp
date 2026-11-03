#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ThermoelasticMesh.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <vector>

using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinAbs;

TEST_CASE("Tapered bead rasters retain changing thickness and nominal volume", "[ThermoelasticMesh]")
{
    const std::vector<Voxel> materials{{{0,0,0}, 1e6, .3}};
    for (const auto section : {BeadSection::RectangularEnvelope, BeadSection::RoundedRectangle}) {
        DepositedSegment bead{{0.,0.,0.}, {.002,0.,.0002}, .0004, 0., 0, section};
        bead.end_height_m = .0002;
        const Eigen::Vector3d spacing(.0001,.00000625,.000003125);
        const auto grid = voxelize_deposition({bead}, materials, spacing, 0., 100000);
        REQUIRE_FALSE(grid.voxels.empty());
        double interior_volume = 0.;
        for (const auto &voxel : grid.voxels) {
            const Eigen::Vector3d center = grid.origin_m + spacing.cwiseProduct(
                Eigen::Vector3d(voxel.index[0] + .5, voxel.index[1] + .5, voxel.index[2] + .5));
            CHECK(center.z() >= 0.);
            CHECK(center.z() < std::clamp(center.x() / .002, 0., 1.) * .0002);
            if (center.x() >= 0. && center.x() < .002) interior_volume += spacing.prod();
        }
        const double nominal = .002 * (.0004 * .0001 -
            (section == BeadSection::RoundedRectangle ? (1. - std::acos(-1.) / 4.) * .0002 * .0002 / 3. : 0.));
        CHECK_THAT(interior_volume, Catch::Matchers::WithinRel(nominal, .015));
        std::swap(bead.start_m, bead.end_m);
        bead.height_m = .0002;
        bead.end_height_m = 0.;
        const auto reverse = voxelize_deposition({bead}, materials, spacing, 0., 100000);
        CHECK_THAT(double(reverse.voxels.size()), Catch::Matchers::WithinRel(double(grid.voxels.size()), .001));
        bead.end_height_m = -.0001;
        CHECK_THROWS_AS(voxelize_deposition({bead}, materials, spacing, 0., 100000), std::invalid_argument);
    }
}

TEST_CASE("Rounded bead sections converge to their nominal area without filling corner voids", "[ThermoelasticMesh]")
{
    const std::vector<Voxel> materials{{{0,0,0}, 1e6, .3}};
    DepositedSegment bead{{0.,0.,.0002}, {.002,0.,.0002}, .0004, .0002, 0,
        BeadSection::RoundedRectangle};
    const double expected = bead.width_m * bead.height_m -
        (1. - std::acos(-1.) / 4.) * bead.height_m * bead.height_m;
    const auto area = [&](unsigned resolution) {
        const Eigen::Vector3d spacing(bead.width_m / 2., bead.height_m / resolution, bead.height_m / resolution);
        const auto grid = voxelize_deposition({bead}, materials, spacing, 0., 100000);
        // Interior section excludes the nominal end caps and any overlap.
        const int column = int(std::floor((.001 - grid.origin_m.x()) / spacing.x()));
        const auto count = std::count_if(grid.voxels.begin(), grid.voxels.end(),
            [&](const Voxel &voxel) { return voxel.index[0] == column; });
        return double(count) * spacing.y() * spacing.z();
    };
    const double coarse = area(8), fine = area(32);
    CHECK(std::abs(fine - expected) < std::abs(coarse - expected));
    CHECK_THAT(fine, Catch::Matchers::WithinRel(expected, .01));
    CHECK(fine < bead.width_m * bead.height_m);
    std::swap(bead.start_m, bead.end_m);
    CHECK_THAT(area(32), WithinAbs(fine, 1e-20));
    bead.width_m = bead.height_m;
    CHECK_THAT(area(32), Catch::Matchers::WithinRel(
        std::acos(-1.) * bead.height_m * bead.height_m / 4., .01));
    bead.width_m = bead.height_m / 2.;
    CHECK_THROWS_AS(voxelize_deposition({bead}, materials, {.000025,.000025,.000025}, 0., 100000),
        std::invalid_argument);
}

TEST_CASE("Explicit rectangular bead sections preserve legacy occupancy and provenance", "[ThermoelasticMesh]")
{
    const std::vector<Voxel> materials{{{0,0,0}, 1e6, .3}};
    DepositedSegment bead{{0.,0.,.0002}, {.001,0.,.0002}, .0004, .0002, 0};
    const Eigen::Vector3d spacing(.0001,.000025,.000025);
    const auto legacy = voxelize_deposition({bead}, materials, spacing, 0., 10000);
    bead.section = BeadSection::RectangularEnvelope;
    const auto explicit_rectangle = voxelize_deposition({bead}, materials, spacing, 0., 10000);
    REQUIRE(explicit_rectangle.voxels.size() == legacy.voxels.size());
    for (std::size_t i = 0; i < legacy.voxels.size(); ++i) {
        CHECK(explicit_rectangle.voxels[i].index == legacy.voxels[i].index);
        CHECK_THAT(explicit_rectangle.deposition[i].path_fraction,
            WithinAbs(legacy.deposition[i].path_fraction, 1e-15));
    }
    bead.section = BeadSection::RoundedRectangle;
    const auto rounded = voxelize_deposition({bead}, materials, spacing, 0., 10000);
    CHECK(rounded.voxels.size() < legacy.voxels.size());
    for (const auto &voxel : rounded.voxels)
        CHECK(std::any_of(legacy.voxels.begin(), legacy.voxels.end(),
            [&](const Voxel &other) { return voxel.index == other.index; }));
    bead.section = static_cast<BeadSection>(99);
    CHECK_THROWS_AS(voxelize_deposition({bead}, materials, spacing, 0., 10000), std::invalid_argument);
}

TEST_CASE("Sloped mechanical beads retain their changing elevation and empty space", "[ThermoelasticMesh]")
{
    const std::vector<Voxel> materials{{{0,0,0}, 2.9e9, .35}};
    DepositedSegment bead{{0.,0.,.0002}, {.002,0.,.0006}, .0004, .0002, 0};
    const Eigen::Vector3d spacing(.0001,.0001,.00005);
    const auto grid = voxelize_deposition({bead}, materials, spacing, 0., 10000);
    REQUIRE_FALSE(grid.voxels.empty());
    for (const auto &voxel : grid.voxels) {
        const auto center = (grid.origin_m + spacing.cwiseProduct(
            Eigen::Vector3d(voxel.index[0] + .5, voxel.index[1] + .5, voxel.index[2] + .5))).eval();
        const double top = .0002 + std::clamp(center.x() / .002, 0., 1.) * .0004;
        CHECK(center.z() >= top - bead.height_m);
        CHECK(center.z() < top);
    }
    std::swap(bead.start_m, bead.end_m);
    const auto reverse = voxelize_deposition({bead}, materials, spacing, 0., 10000);
    REQUIRE(reverse.voxels.size() == grid.voxels.size());
    for (std::size_t i = 0; i < grid.voxels.size(); ++i)
        CHECK(reverse.voxels[i].index == grid.voxels[i].index);
    CHECK_THROWS_AS(voxelize_deposition({bead}, materials, spacing, 0., 10000, {}, 10), ResourceLimit);
    bead.start_m = {0.,0.,.0001};
    CHECK_THROWS_AS(voxelize_deposition({bead}, materials, spacing, 0., 10000), std::invalid_argument);
    bead.start_m = {0.,0.,.0004};
    bead.end_m = {0.,0.,.0008};
    CHECK_THROWS_AS(voxelize_deposition({bead}, materials, spacing, 0., 10000), std::invalid_argument);
}

TEST_CASE("Deposition provenance preserves first coverage and path direction", "[ThermoelasticMesh]")
{
    const std::vector<Voxel> materials{{{0,0,0}, 1e6, .3}};
    const DepositedSegment bead{{0.,0.,.0002}, {.001,0.,.0002}, .0004, .0002, 0};
    auto reversed = bead;
    std::swap(reversed.start_m, reversed.end_m);
    const Eigen::Vector3d spacing(.0001,.0001,.0001);
    const auto forward = voxelize_deposition({bead, reversed}, materials, spacing, 0., 1000);
    const auto backward = voxelize_deposition({reversed}, materials, spacing, 0., 1000);
    REQUIRE(forward.voxels.size() == backward.voxels.size());
    REQUIRE(forward.deposition.size() == forward.voxels.size());
    REQUIRE(backward.deposition.size() == backward.voxels.size());
    for (std::size_t i = 0; i < forward.voxels.size(); ++i) {
        CHECK(forward.voxels[i].index == backward.voxels[i].index);
        CHECK(forward.deposition[i].segment_index == 0);
        CHECK(forward.deposition[i].path_fraction >= 0.);
        CHECK(forward.deposition[i].path_fraction <= 1.);
        CHECK_THAT(forward.deposition[i].path_fraction + backward.deposition[i].path_fraction,
            WithinAbs(1., 1e-12));
    }
}
