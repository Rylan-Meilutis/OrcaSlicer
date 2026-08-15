#include <catch2/catch_all.hpp>

#include <algorithm>

#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/ExtrusionEntityCollection.hpp"
#include "libslic3r/NonplanarSurface.hpp"
#include "libslic3r/SLA/IndexedMesh.hpp"
#include "libslic3r/TriangleMesh.hpp"

using namespace Slic3r;

namespace {

void append_sloped_plane(indexed_triangle_set &mesh, float z)
{
    const int first = int(mesh.vertices.size());
    mesh.vertices.emplace_back(0.f, 0.f, z);
    mesh.vertices.emplace_back(10.f, 0.f, z + 1.f);
    mesh.vertices.emplace_back(10.f, 10.f, z + 1.f);
    mesh.vertices.emplace_back(0.f, 10.f, z);
    mesh.indices.emplace_back(first, first + 1, first + 2);
    mesh.indices.emplace_back(first, first + 2, first + 3);
}

indexed_triangle_set crowned_surface(float flat_begin, float flat_end)
{
    indexed_triangle_set mesh;
    const std::array<float, 4> x{0.f, flat_begin, flat_end, flat_end + flat_begin};
    const std::array<float, 4> z{0.f, 1.f, 1.f, 0.f};
    for (size_t column = 0; column < x.size(); ++column) {
        mesh.vertices.emplace_back(x[column], 0.f, z[column]);
        mesh.vertices.emplace_back(x[column], 10.f, z[column]);
    }
    for (int column = 0; column < 3; ++column) {
        const int first = 2 * column;
        mesh.indices.emplace_back(first, first + 2, first + 3);
        mesh.indices.emplace_back(first, first + 3, first + 1);
    }
    return mesh;
}

} // namespace

TEST_CASE("Non-planar finishing selects only the uppermost sloped surface",
          "[NonplanarSurface][Regression]")
{
    indexed_triangle_set triangles;
    append_sloped_plane(triangles, 0.f);
    append_sloped_plane(triangles, 5.f);
    const TriangleMesh mesh(triangles);
    const sla::IndexedMesh indexed(mesh);

    const std::vector<uint8_t> selected = select_nonplanar_surface_facets(
        indexed, 45., 0.2, 1.);

    REQUIRE(selected.size() == 4);
    CHECK(selected[0] == 0);
    CHECK(selected[1] == 0);
    CHECK(selected[2] == 1);
    CHECK(selected[3] == 1);
}

TEST_CASE("Non-planar finishing joins a narrow flat roof crown but not a broad deck",
          "[NonplanarSurface][Regression]")
{
    for (const auto &[flat_end, expected_flat] :
         std::array<std::pair<float, uint8_t>, 2>{{{4.f, 1}, {20.f, 0}}}) {
        CAPTURE(flat_end);
        const TriangleMesh mesh(crowned_surface(1.f, flat_end));
        const sla::IndexedMesh indexed(mesh);

        const std::vector<uint8_t> selected = select_nonplanar_surface_facets(
            indexed, 45., 0.2, 1.);

        REQUIRE(selected.size() == 6);
        CHECK(selected[0] == 1);
        CHECK(selected[1] == 1);
        CHECK(selected[2] == expected_flat);
        CHECK(selected[3] == expected_flat);
        CHECK(selected[4] == 1);
        CHECK(selected[5] == 1);
    }
}

TEST_CASE("Non-planar candidate projections are limited to reachable layer bands",
          "[NonplanarSurface][Regression]")
{
    indexed_triangle_set triangles;
    append_sloped_plane(triangles, 5.f);
    const TriangleMesh mesh(triangles);
    const sla::IndexedMesh indexed(mesh);
    const std::vector<uint8_t> selected = select_nonplanar_surface_facets(
        indexed, 45., 0.2, 1.);

    const std::vector<std::vector<ExPolygons>> projections =
        project_nonplanar_surface_facets(
            indexed, {selected}, {{0., 0.4}, {4.9, 5.4}, {5.8, 6.2}});

    REQUIRE(projections.size() == 1);
    REQUIRE(projections.front().size() == 3);
    CHECK(projections.front()[0].empty());
    CHECK_FALSE(projections.front()[1].empty());
    CHECK_FALSE(projections.front()[2].empty());
}

TEST_CASE("Brick walls round trip through the preview role tag",
          "[NonplanarSurface][GCodePreview]")
{
    CHECK(ExtrusionEntity::role_to_string(erStaggeredPerimeter) == "Brick wall");
    CHECK(ExtrusionEntity::string_to_role("Brick wall") == erStaggeredPerimeter);
    // Preserve preview compatibility with G-code produced by older builds.
    CHECK(ExtrusionEntity::string_to_role("Brick inner wall") == erStaggeredPerimeter);
    CHECK(ExtrusionEntity::role_to_string(erSmoothOuterWall) == "Smooth outer wall");
    CHECK(ExtrusionEntity::string_to_role("Smooth outer wall") == erSmoothOuterWall);
    CHECK(ExtrusionEntity::role_to_string(erNonplanarInterlockingWall) ==
          "Interlocking inner wall");
    CHECK(ExtrusionEntity::string_to_role("Interlocking inner wall") ==
          erNonplanarInterlockingWall);
}

TEST_CASE("Non-planar infill round trips through its preview role tag",
          "[NonplanarSurface][GCodePreview]")
{
    CHECK(ExtrusionEntity::role_to_string(erNonplanarInfill) ==
          "Non-planar interlocking infill");
    CHECK(ExtrusionEntity::string_to_role("Non-planar interlocking infill") ==
          erNonplanarInfill);
}

TEST_CASE("Non-planar transition metadata survives path copies",
          "[NonplanarSurface][GCodePreview]")
{
    ExtrusionPath original(erTopSolidInfill);
    original.nonplanar_surface = true;
    original.nonplanar_transition = true;
    original.nonplanar_clearance_validated = true;
    original.nonplanar_before_current_layer = true;
    original.nonplanar_feature_transition = true;
    original.nonplanar_feature_course = 7;
    original.nonplanar_leveling_transition = true;

    const ExtrusionPath copy(original);
    CHECK(copy.nonplanar_surface);
    CHECK(copy.nonplanar_transition);
    CHECK(copy.nonplanar_clearance_validated);
    CHECK(copy.nonplanar_before_current_layer);
    CHECK(copy.nonplanar_feature_transition);
    CHECK(copy.nonplanar_feature_course == 7);
    CHECK(copy.nonplanar_leveling_transition);
}

TEST_CASE("Clipping a contoured extrusion preserves interpolated Z",
          "[NonplanarSurface][Regression]")
{
    ExtrusionPath source(erExternalPerimeter);
    source.z_contoured = true;
    source.nonplanar_surface = true;
    source.polyline.points = {
        Point3(scale_(0.), scale_(0.), scale_(0.)),
        Point3(scale_(10.), scale_(0.), scale_(1.))};

    ExPolygon clip;
    clip.contour.points = {
        Point(scale_(2.), scale_(-1.)), Point(scale_(8.), scale_(-1.)),
        Point(scale_(8.), scale_(1.)), Point(scale_(2.), scale_(1.))};

    ExtrusionEntityCollection intersection;
    source.intersect_expolygons({clip}, &intersection);
    REQUIRE(intersection.entities.size() == 1);
    const auto *inside =
        dynamic_cast<const ExtrusionPath *>(intersection.entities.front());
    REQUIRE(inside != nullptr);
    REQUIRE(inside->polyline.points.size() == 2);
    const double first_z = unscale_(inside->polyline.points.front().z());
    const double last_z = unscale_(inside->polyline.points.back().z());
    CHECK_THAT(std::min(first_z, last_z),
               Catch::Matchers::WithinAbs(0.2, 1e-6));
    CHECK_THAT(std::max(first_z, last_z),
               Catch::Matchers::WithinAbs(0.8, 1e-6));
    CHECK(inside->nonplanar_surface);

    ExtrusionEntityCollection difference;
    source.subtract_expolygons({clip}, &difference);
    REQUIRE(difference.entities.size() == 2);
    for (const ExtrusionEntity *entity : difference.entities) {
        const auto *outside = dynamic_cast<const ExtrusionPath *>(entity);
        REQUIRE(outside != nullptr);
        REQUIRE(outside->polyline.points.size() == 2);
        CHECK(outside->polyline.points.front().z() !=
              outside->polyline.points.back().z());
        CHECK(outside->nonplanar_surface);
    }
}
