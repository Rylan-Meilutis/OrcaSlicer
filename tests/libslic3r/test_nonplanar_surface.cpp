#include <catch2/catch_all.hpp>

#include <algorithm>
#include <array>
#include <memory>

#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/ExtrusionEntityCollection.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/NonplanarSurface.hpp"
#include "libslic3r/SLA/IndexedMesh.hpp"
#include "libslic3r/TriangleMesh.hpp"

using namespace Slic3r;

TEST_CASE("Non-planar wall anchors retain native inset ownership and variable widths",
          "[NonplanarSurface][NativeWallAnchors][Regression]")
{
    const bool multipath = GENERATE(false, true);
    ExtrusionEntityCollection native;
    for (int inset = 0; inset < 3; ++inset) {
        ExtrusionPaths parts;
        for (const double width : {0.15, 0.6}) {
            ExtrusionPath path(inset == 0 ? erExternalPerimeter : erPerimeter,
                               0.08, width, 0.2);
            path.inset_idx = 0; // The native parent, not the child, owns the inset.
            path.polyline.points = {
                Point3(scale_(width), scale_(double(inset)), scale_(0.1)),
                Point3(scale_(width + 1.), scale_(double(inset)), scale_(0.1))};
            parts.emplace_back(std::move(path));
        }
        ExtrusionEntity *owner = multipath ?
            static_cast<ExtrusionEntity *>(new ExtrusionMultiPath(std::move(parts))) :
            static_cast<ExtrusionEntity *>(new ExtrusionLoop(std::move(parts)));
        owner->inset_idx = inset;
        native.entities.push_back(owner);
    }
    for (int inset = 0; inset < 3; ++inset) {
        const Polylines lines = nonplanar_wall_centerlines(native, inset);
        REQUIRE(lines.size() == 2);
        for (const Polyline &line : lines)
            for (const Point &point : line.points)
                CHECK(point.y() == scale_(double(inset)));
    }
    CHECK(nonplanar_wall_centerlines(native, 3).empty());
}

TEST_CASE("A native foundation must contact fill as well as walls",
          "[NonplanarSurface][NativeFoundation][Regression]")
{
    const ExtrusionRole role = GENERATE(erExternalPerimeter, erPerimeter,
        erInternalInfill, erSolidInfill, erTopSolidInfill, erBridgeInfill);
    const double owner = GENERATE(0., 12.);
    ExtrusionPath path(role, 0.032, 0.4, 0.08);
    path.polyline.points = {Point3(scale_(0.), scale_(0.), scale_(0.08 - owner)),
                           Point3(scale_(10.), scale_(0.), scale_(0.08 - owner))};
    ExtrusionPath below(erBridgeInfill, 0.08, 0.4, 0.2);
    below.polyline.points = {Point3(scale_(0.), scale_(0.), scale_(0.)),
                            Point3(scale_(10.), scale_(0.), scale_(0.))};
    std::vector<const ExtrusionPath *> support{&below};
    CHECK(nonplanar_course_is_supported(path, support, owner, 0.));

    // A wall on either side does not certify fill spanning the empty middle.
    below.polyline.points.back().x() = scale_(1.);
    CHECK_FALSE(nonplanar_course_is_supported(path, support, owner, 0.));
    if (!is_perimeter(role))
        CHECK(nonplanar_walls_are_supported(path, support, owner, 0.));
    below.polyline.points.back().x() = scale_(10.);
    for (Point3 &point : below.polyline.points)
        point.z() = scale_(-0.2);
    CHECK_FALSE(nonplanar_course_is_supported(path, support, owner, 0.));
}

TEST_CASE("Dense non-planar courses cannot use omitted sparse rows as support",
          "[NonplanarSurface][CourseOccupancy][Regression]")
{
    const double owner = GENERATE(0., 20.);
    const double slope = GENERATE(0., 0.1);
    const auto row = [&](double y, double height) {
        ExtrusionPath path(erSolidInfill, 0.04, 0.4, 0.1);
        path.polyline.points = {
            Point3(scale_(0.), scale_(y), scale_(height - owner)),
            Point3(scale_(10.), scale_(y), scale_(height + slope * 10. - owner))};
        return path;
    };
    const ExtrusionPath left = row(0., 1.);
    const ExtrusionPath omitted = row(1., 1.);
    const ExtrusionPath right = row(2., 1.);
    const ExtrusionPath dense = row(1., 1.1);
    const std::vector<const ExtrusionPath *> prototype{&left, &omitted, &right};
    const std::vector<const ExtrusionPath *> deposited{&left, &right};
    CHECK(nonplanar_course_is_supported(dense, prototype, owner, owner));
    CHECK_FALSE(nonplanar_course_is_supported(dense, deposited, owner, owner));

    // Routing must also validate connectors: both rows are supported, but a
    // join across their ends would traverse the missing middle row.
    const ExtrusionPath first = row(0., 1.1);
    const ExtrusionPath last = row(2., 1.1);
    CHECK(nonplanar_course_is_supported(first, deposited, owner, owner));
    CHECK(nonplanar_course_is_supported(last, deposited, owner, owner));
    ExtrusionPath connector(first);
    connector.polyline.points = {first.last_point3(), last.last_point3()};
    CHECK_FALSE(nonplanar_course_is_supported(connector, deposited, owner, owner));

    // A provisional nominal height can falsely certify a later, thinner
    // course. Validation must use the same effective height as emission.
    ExtrusionPath thin = row(0., 1.2);
    thin.height = 0.25;
    CHECK(nonplanar_course_is_supported(thin, deposited, owner, owner));
    thin.height = 0.08;
    CHECK_FALSE(nonplanar_course_is_supported(thin, deposited, owner, owner));
}

TEST_CASE("Short Arachne wall fragments survive non-planar fill cleanup",
          "[NonplanarSurface][WallSelection][Regression]")
{
    const ExtrusionRole role = GENERATE(erExternalPerimeter, erPerimeter,
                                        erOverhangPerimeter, erTopSolidInfill);
    ExtrusionPath path(role, 0.08, 0.4, 0.2);
    path.polyline.points = {Point3(scale_(0.), scale_(0.), scale_(0.)),
                           Point3(scale_(0.2), scale_(0.), scale_(0.02))};
    ExtrusionEntityCollection course;
    course.entities.push_back(path.clone());
    remove_short_nonplanar_fill_paths(course.entities);
    if (is_perimeter(role)) {
        REQUIRE(course.entities.size() == 1);
        const auto *wall = dynamic_cast<const ExtrusionPath *>(course.entities.front());
        REQUIRE(wall != nullptr);
        CHECK(wall->polyline.points == path.polyline.points);
        CHECK(wall->role() == role);
    } else {
        CHECK(course.empty());
    }
}

TEST_CASE("Non-planar conflict filtering cannot delete a required perimeter",
          "[NonplanarSurface][WallSelection][Regression]")
{
    const size_t rejected = GENERATE(size_t(0), size_t(1), size_t(2), size_t(3));
    ExtrusionEntityCollection course;
    for (int inset = 0; inset < 4; ++inset) {
        ExtrusionPath path(inset == 0 ? erExternalPerimeter :
                           inset < 3 ? erPerimeter : erTopSolidInfill,
                           0.08, 0.4, 0.2);
        path.inset_idx = inset;
        path.polyline.points = {Point3(scale_(0.), scale_(double(inset)), scale_(0.)),
                               Point3(scale_(10.), scale_(double(inset)), scale_(1.))};
        course.entities.push_back(path.clone());
    }
    const ExtrusionEntitiesPtr original = course.entities;
    std::vector<uint8_t> retained(4, 1);
    retained[rejected] = 0;
    retained[3] = 0; // Rejecting a wall must not partly delete the fill either.
    const bool accepted = filter_nonplanar_course_paths(course.entities, retained);
    if (rejected < 3) {
        CHECK_FALSE(accepted);
        CHECK(course.entities == original);
    } else {
        CHECK(accepted);
        REQUIRE(course.entities.size() == 3);
        for (size_t i = 0; i < course.entities.size(); ++i)
            CHECK(course.entities[i] == original[i]);
    }
}

TEST_CASE("Non-planar solid courses honor physical thickness as well as layer count",
          "[NonplanarSurface][SolidCourseDepth][Regression]")
{
    const double translation = GENERATE(0., -24., 100.);
    std::vector<NonplanarCourseProfile> courses;
    for (size_t i = 0; i < 6; ++i)
        courses.push_back({translation + 1. + 0.2 * double(i), 0.2 * double(i)});
    const double low = translation + 1.5;
    const double high = translation + 2.;
    // Low-side courses are only 0.1 mm thick despite 0.2 mm nominal spacing.
    CHECK(nonplanar_solid_course_begin(courses, low, high, 2, 0.) == 4);
    CHECK(nonplanar_solid_course_begin(courses, low, high, 2, 0.3) == 3);
    CHECK(nonplanar_solid_course_begin(courses, low, high, 2, 0.35) == 2);
    CHECK(nonplanar_solid_course_begin(courses, low, high, 4, 0.1) == 2);
    CHECK(nonplanar_solid_course_begin(courses, low, high, 20, 0.1) == 0);
    CHECK(nonplanar_solid_course_begin(courses, low, high, 2, 2.) == 0);
    CHECK(nonplanar_solid_course_begin({}, low, high, 2, 0.3) == 0);
    // A flat stack already meets the requested thickness with two courses.
    CHECK(nonplanar_solid_course_begin(courses, high, high, 2, 0.35) == 4);
}

TEST_CASE("Non-planar wall support covers the whole strand at its actual bead height",
          "[NonplanarSurface][WallSupport][Regression]")
{
    const double owner_z = GENERATE(0., 12.);
    const double slope = GENERATE(0., 0.15);
    const bool reverse = GENERATE(false, true);
    ExtrusionPath wall(erExternalPerimeter, 0.032, 0.4, 0.08);
    wall.polyline.points = {
        Point3(scale_(0.), scale_(0.), scale_(0.08 - owner_z)),
        Point3(scale_(10.), scale_(0.), scale_(0.08 + 10. * slope - owner_z))};
    ExtrusionPath below(erExternalPerimeter, 0.08, 0.4, 0.2);
    below.polyline.points = {
        Point3(0, 0, 0), Point3(scale_(10.), scale_(0.), scale_(10. * slope))};
    if (reverse) {
        wall.reverse();
        below.reverse();
    }
    std::vector<const ExtrusionPath *> supports{&below};
    CHECK(nonplanar_walls_are_supported(wall, supports, owner_z, 0.));

    // A thin course cannot borrow the old hard-coded 0.26 mm allowance.
    CHECK_FALSE(nonplanar_walls_are_supported(wall, supports, owner_z + 0.17, 0.));
    CHECK_FALSE(nonplanar_walls_are_supported(wall, supports, owner_z - 0.08, 0.));

    // Only the midpoint is supported. Both long ends are still floating.
    below.polyline.points = {
        Point3(scale_(4.8), scale_(0.), scale_(4.8 * slope)),
        Point3(scale_(5.2), scale_(0.), scale_(5.2 * slope))};
    CHECK_FALSE(nonplanar_walls_are_supported(wall, supports, owner_z, 0.));
    CHECK_FALSE(nonplanar_walls_are_supported(wall, supports, owner_z, 0., true));

    // Adjacent deposited paths may jointly provide the complete foundation.
    ExtrusionPath other = below;
    below.polyline.points = {Point3(0, 0, 0), Point3(scale_(5.), scale_(0.), scale_(5. * slope))};
    other.polyline.points = {
        Point3(scale_(5.), scale_(0.), scale_(5. * slope)),
        Point3(scale_(10.), scale_(0.), scale_(10. * slope))};
    supports.push_back(&other);
    CHECK(nonplanar_walls_are_supported(wall, supports, owner_z, 0.));
}

TEST_CASE("Splitting a non-planar wall does not reset its unsupported span",
          "[NonplanarSurface][WallSupport][Regression]")
{
    ExtrusionPath below(erExternalPerimeter, 0.08, 0.4, 0.2);
    below.polyline.points = {Point3(0, 0, 0), Point3(scale_(2.), scale_(0.), scale_(0.))};
    const std::vector<const ExtrusionPath *> supports{&below};
    ExtrusionMultiPath wall;
    for (size_t i = 0; i < 3; ++i) {
        ExtrusionPath part(erExternalPerimeter, 0.032, 0.4, 0.08);
        part.polyline.points = {
            Point3(scale_(0.4 * i), scale_(0.), scale_(0.25)),
            Point3(scale_(0.4 * (i + 1)), scale_(0.), scale_(0.25))};
        // Every individual short span fits the cantilever allowance; their
        // connected 1.2 mm wall does not. Arachne may split it by bead width.
        CHECK(nonplanar_walls_are_supported(part, supports, 0., 0., true));
        wall.paths.emplace_back(std::move(part));
    }
    CHECK_FALSE(nonplanar_walls_are_supported(wall, supports, 0., 0., true));
    CHECK(nonplanar_walls_are_supported(wall, supports, -0.17, 0., true));
}

TEST_CASE("Mixed non-planar walls preserve planar joins in both directions",
          "[NonplanarSurface][WallCourseBlend][Regression]")
{
    const double planar_z = 10.;
    const double surface_z = GENERATE(8., 12.);
    const double course_z = GENERATE(9., 11.);
    std::array<double, 21> entry;
    for (size_t i = 0; i < entry.size(); ++i) {
        const double t = double(i) / double(entry.size() - 1);
        const double influence = t * t * (3. - 2. * t);
        const double projected_z = planar_z + influence * (surface_z - planar_z);
        entry[i] = nonplanar_wall_course_z(
            planar_z, projected_z, surface_z, course_z, false);
        CHECK_THAT(entry[i], Catch::Matchers::WithinAbs(
            planar_z + influence * (course_z - planar_z), 1e-9));
        const double reverse_t = 1. - t;
        const double reverse_influence = reverse_t * reverse_t * (3. - 2. * reverse_t);
        CHECK_THAT(nonplanar_wall_course_z(planar_z,
            planar_z + reverse_influence * (surface_z - planar_z),
            surface_z, course_z, false), Catch::Matchers::WithinAbs(
                planar_z + reverse_influence * (course_z - planar_z), 1e-9));
    }
    CHECK_THAT(entry.front(), Catch::Matchers::WithinAbs(planar_z, 1e-9));
    CHECK_THAT(entry.back(), Catch::Matchers::WithinAbs(course_z, 1e-9));
    CHECK_THAT(nonplanar_wall_course_z(planar_z, planar_z, planar_z, course_z, false),
        Catch::Matchers::WithinAbs(planar_z, 1e-9));
    CHECK_THAT(nonplanar_wall_course_z(planar_z, planar_z, surface_z, course_z, true),
        Catch::Matchers::WithinAbs(course_z, 1e-9));
    CHECK_THAT(nonplanar_wall_course_z(planar_z, planar_z, planar_z, course_z, true),
        Catch::Matchers::WithinAbs(course_z, 1e-9));
}

TEST_CASE("Non-planar courses rise monotonically while adopting the surface slope",
          "[NonplanarSurface][NonplanarCoursePlan][Regression]")
{
    const double span = GENERATE(0., 0.4, 2.2265);
    const double minimum_z = 38.036;
    const double maximum_z = minimum_z + span;
    const double anchor = 36.6;
    const double minimum_height = 0.08, maximum_height = 0.3;
    const auto courses = plan_nonplanar_courses(anchor, minimum_z, maximum_z,
        5, minimum_height, maximum_height);
    REQUIRE(courses.size() >= 6);
    CHECK_THAT(courses.front().plane_z, Catch::Matchers::WithinAbs(anchor, 1e-9));
    CHECK_THAT(courses.front().blend, Catch::Matchers::WithinAbs(0., 1e-9));
    CHECK_THAT(courses.back().plane_z, Catch::Matchers::WithinAbs(maximum_z, 1e-9));
    CHECK_THAT(courses.back().blend, Catch::Matchers::WithinAbs(1., 1e-9));
    for (int sample = 0; sample <= 20; ++sample) {
        const double surface = minimum_z + span * sample / 20.;
        double previous_z = anchor;
        for (size_t idx = 1; idx < courses.size(); ++idx) {
            const auto &course = courses[idx];
            const double z = course.plane_z + course.blend * (surface - maximum_z);
            CAPTURE(span, sample, idx, z, previous_z);
            CHECK(z - previous_z >= minimum_height - 1e-9);
            CHECK(z - previous_z <= maximum_height + 1e-9);
            CHECK(course.blend >= courses[idx - 1].blend);
            previous_z = z;
        }
        CHECK_THAT(previous_z, Catch::Matchers::WithinAbs(surface, 1e-9));
    }
    // The previous Benchy anchor did not leave enough low-side depth for
    // printable courses. It must request an earlier foundation, not reverse
    // or squeeze courses when entering the solid top-shell portion.
    CHECK(plan_nonplanar_courses(37.2, 38.036, 40.2625, 5, 0.08, 0.3).empty());
    CHECK(plan_nonplanar_courses(anchor, minimum_z, maximum_z, 5, 0., 0.3).empty());
    CHECK(plan_nonplanar_courses(minimum_z, minimum_z, maximum_z, 5, 0.08, 0.3).empty());
}

TEST_CASE("Non-planar course planning respects a nominal layer height ceiling",
          "[NonplanarSurface][NonplanarCoursePlan][Regression]")
{
    const double minimum_z = 38.036, maximum_z = 40.2625;
    const double nominal_height = 0.2, minimum_height = 0.08;
    // A shallower foundation cannot fit all courses at both ends of the slope.
    CHECK(plan_nonplanar_courses(36.6, minimum_z, maximum_z,
        5, minimum_height, nominal_height).empty());
    const auto courses = plan_nonplanar_courses(35.6, minimum_z, maximum_z,
        5, minimum_height, nominal_height);
    REQUIRE_FALSE(courses.empty());
    for (const double surface : {minimum_z, (minimum_z + maximum_z) / 2., maximum_z}) {
        double previous = courses.front().plane_z;
        for (size_t idx = 1; idx < courses.size(); ++idx) {
            const auto &course = courses[idx];
            const double current = course.plane_z + course.blend * (surface - maximum_z);
            CAPTURE(surface, idx);
            CHECK(current - previous >= minimum_height - 1e-9);
            CHECK(current - previous <= nominal_height + 1e-9);
            previous = current;
        }
        CHECK_THAT(previous, Catch::Matchers::WithinAbs(surface, 1e-9));
    }
}

TEST_CASE("Continuing feature courses extend the roof slope across an opening",
          "[NonplanarSurface][NonplanarProjection][Regression]")
{
    const double slope_x = GENERATE(-0.2, 0., 0.2);
    const double slope_y = GENERATE(-0.1, 0.1);
    const auto height = [&](double x, double y) { return 10. + slope_x * x + slope_y * y; };
    indexed_triangle_set triangles;
    const auto rectangle = [&](float x0, float y0, float x1, float y1, double lift = 0.) {
        const int base = int(triangles.vertices.size());
        for (const Vec2f &xy : std::array<Vec2f, 4>{
                 Vec2f(x0, y0), Vec2f(x1, y0), Vec2f(x1, y1), Vec2f(x0, y1)})
            triangles.vertices.emplace_back(xy.x(), xy.y(), float(height(xy.x(), xy.y()) + lift));
        triangles.indices.emplace_back(base, base + 1, base + 2);
        triangles.indices.emplace_back(base, base + 2, base + 3);
    };
    // Actual mesh opening, not a masked but otherwise continuous plane.
    rectangle(0, 0, 10, 4);
    rectangle(0, 6, 10, 10);
    rectangle(0, 4, 4, 6);
    rectangle(6, 4, 10, 6);
    const size_t roof_faces = triangles.indices.size();
    rectangle(0, 0, 10, 10, 5.); // Unselected feature surface must not supply Z.
    const TriangleMesh mesh(triangles);
    sla::IndexedMesh indexed(mesh);
    indexed.ground_level_offset(GENERATE(-100., 0., 100.));
    std::vector<uint8_t> facets(triangles.indices.size(), 0);
    std::fill_n(facets.begin(), roof_faces, uint8_t(1));
    const auto polygon = [](double low, double high) {
        return Polygon(Points{Point(scale_(low), scale_(low)), Point(scale_(high), scale_(low)),
                              Point(scale_(high), scale_(high)), Point(scale_(low), scale_(high))});
    };
    ExPolygon domain(polygon(0., 10.));
    Polygon hole = polygon(4., 6.);
    hole.reverse();
    domain.holes.push_back(hole);
    const ExPolygons projection{domain};
    // Cross the opening in both directions, including its boundary. Every
    // sample must remain on the same plane, regardless of the nearest edge.
    for (int axis = 0; axis < 2; ++axis)
        for (int step = 0; step <= 16; ++step) {
            const double x = axis == 0 ? 3. + 0.25 * step : 5.;
            const double y = axis == 1 ? 3. + 0.25 * step : 5.;
            CAPTURE(slope_x, slope_y, x, y);
            const auto z = project_nonplanar_patch_height(indexed, projection, facets,
                Point(scale_(x), scale_(y)));
            REQUIRE(z.has_value());
            CHECK_THAT(*z, Catch::Matchers::WithinAbs(height(x, y), 2e-5));
        }
    const auto exterior = project_nonplanar_patch_height(indexed, projection, facets,
        Point(scale_(-1.), scale_(5.)));
    REQUIRE(exterior.has_value());
    CHECK_THAT(*exterior, Catch::Matchers::WithinAbs(height(0., 5.), 0.002));
    CHECK_FALSE(project_nonplanar_patch_height(indexed, {}, facets, Point(0, 0)).has_value());
    CHECK_FALSE(project_nonplanar_patch_height(indexed, projection,
        std::vector<uint8_t>(facets.size(), 0), Point(scale_(5.), scale_(5.))).has_value());
}

TEST_CASE("Local non-planar wall heights preserve the connected geometry and flow multiplier",
          "[NonplanarSurface][NonplanarLocalHeight][Regression]")
{
    const double multiplier = 1.1;
    ExtrusionPath source(erPerimeter, Flow(0.45f, 0.2f, 0.f).mm3_per_mm() * multiplier, 0.45, 0.2);
    source.inset_idx = 2;
    source.nonplanar_feature_transition = true;
    source.nonplanar_surface = source.z_contoured = source.nonplanar_transition = true;
    source.nonplanar_feature_course = 3;
    source.polyline.points = {
        Point3(0, 0, 0), Point3(scale_(1.), scale_(0.), scale_(0.1)),
        Point3(scale_(2.), scale_(0.), scale_(0.1)), Point3(scale_(3.), scale_(0.), scale_(0.1)),
        Point3(0, 0, 0)};
    const auto parts = split_nonplanar_path_by_height(source, {0.2, 0.29, 0.29, 0.2});
    REQUIRE(parts.size() == 3);
    CHECK(parts[1].polyline.points.size() == 3);
    CHECK(parts.front().first_point3() == source.first_point3());
    CHECK(parts.back().last_point3() == source.last_point3());
    Points3 joined = parts.front().polyline.points;
    for (size_t idx = 0; idx < parts.size(); ++idx) {
        const auto &part = parts[idx];
        if (idx > 0) {
            CHECK(parts[idx - 1].last_point3() == part.first_point3());
            joined.insert(joined.end(), part.polyline.points.begin() + 1, part.polyline.points.end());
        }
        CHECK(part.inset_idx == 2);
        CHECK(part.nonplanar_feature_transition);
        CHECK(part.nonplanar_feature_course == 3);
        CHECK(part.role() == source.role());
        CHECK_THAT(part.width, Catch::Matchers::WithinAbs(source.width, 1e-6));
        CHECK_THAT(part.mm3_per_mm,
            Catch::Matchers::WithinAbs(Flow(part.width, part.height, 0.f).mm3_per_mm() * multiplier, 1e-7));
    }
    CHECK_THAT(parts[1].height, Catch::Matchers::WithinAbs(0.29, 1e-6));
    CHECK(parts[1].mm3_per_mm > source.mm3_per_mm);
    CHECK_FALSE(parts[1].nonplanar_transition);
    CHECK(parts[1].nonplanar_schedule_owned);
    CHECK(joined == source.polyline.points);
}

TEST_CASE("Connected variable-width walls preserve their inset when copied",
          "[NonplanarSurface][Regression]")
{
    ExtrusionPath first(erPerimeter, 0.08, 0.4, 0.2);
    first.inset_idx = 2;
    first.nonplanar_surface = first.z_contoured = true;
    first.polyline.points = {Point3(0, 0, 0), Point3(scale_(1.), scale_(0.), scale_(0.1))};
    ExtrusionPath second(first);
    second.width = 0.5;
    second.mm3_per_mm = 0.1;
    second.polyline.points = {first.last_point3(), Point3(scale_(2.), scale_(0.), scale_(0.2))};
    ExtrusionMultiPath source(ExtrusionPaths{first, second});
    source.inset_idx = 2;
    source.set_reverse();
    ExtrusionMultiPath copy(source);
    ExtrusionMultiPath moved(std::move(copy));
    ExtrusionMultiPath assigned;
    assigned = source;
    ExtrusionMultiPath move_assigned;
    move_assigned = std::move(assigned);
    std::unique_ptr<ExtrusionEntity> cloned(source.clone());
    ExtrusionMultiPath move_source(source);
    std::unique_ptr<ExtrusionEntity> move_cloned(move_source.clone_move());
    for (const ExtrusionEntity *entity : std::array<const ExtrusionEntity *, 4>{
             &moved, &move_assigned, cloned.get(), move_cloned.get()}) {
        const auto *wall = dynamic_cast<const ExtrusionMultiPath *>(entity);
        REQUIRE(wall != nullptr);
        CHECK(wall->inset_idx == 2);
        CHECK_FALSE(wall->can_reverse());
        REQUIRE(wall->paths.size() == 2);
        CHECK(wall->paths.front().last_point3() == wall->paths.back().first_point3());
        CHECK_THAT(wall->paths.front().width, Catch::Matchers::WithinAbs(0.4, 1e-6));
        CHECK_THAT(wall->paths.back().width, Catch::Matchers::WithinAbs(0.5, 1e-6));
        CHECK_THAT(wall->paths.back().mm3_per_mm, Catch::Matchers::WithinAbs(0.1, 1e-6));
    }
}

TEST_CASE("Non-planar replacement clips the physical Z band of a sloping segment",
          "[NonplanarSurface][NonplanarOverlap]")
{
    const bool reversed = GENERATE(false, true);
    Vec3d a(0., 0., 0.), b(10., 0., 1.);
    if (reversed)
        std::swap(a, b);
    const auto intervals = nonplanar_segment_overlap(a, b,
        Vec3d(0., 0., 0.3), Vec3d(10., 0., 0.3), 0.2, 0.05);
    REQUIRE(intervals.size() == 1);
    CHECK_THAT(intervals.front().first, Catch::Matchers::WithinAbs(reversed ? 0.65 : 0.25, 1e-9));
    CHECK_THAT(intervals.front().second, Catch::Matchers::WithinAbs(reversed ? 0.75 : 0.35, 1e-9));
    const auto touching = nonplanar_segment_overlap(a, b,
        Vec3d(0., 0.1, 0.3), Vec3d(10., 0.1, 0.3), 0.2, 0.05);
    REQUIRE(touching.size() == 1);
    CHECK_THAT(touching.front().first, Catch::Matchers::WithinAbs(intervals.front().first, 1e-9));
    CHECK_THAT(touching.front().second, Catch::Matchers::WithinAbs(intervals.front().second, 1e-9));
    CHECK(nonplanar_segment_overlap(a, b, Vec3d(0., 1., 0.3),
        Vec3d(10., 1., 0.3), 0.2, 0.05).empty());
    CHECK(nonplanar_segment_overlap(a, b, Vec3d(0., 0., 2.),
        Vec3d(10., 0., 2.), 0.2, 0.05).empty());
}

TEST_CASE("Non-planar replacement respects finite segment endpoints",
          "[NonplanarSurface][NonplanarOverlap]")
{
    const auto intervals = nonplanar_segment_overlap(
        Vec3d(0., 0., 1.), Vec3d(10., 0., 1.),
        Vec3d(4., 0., 1.), Vec3d(6., 0., 1.), 0.2, 0.05);
    REQUIRE_FALSE(intervals.empty());
    CHECK_THAT(intervals.front().first, Catch::Matchers::WithinAbs(0.38, 1e-9));
    CHECK_THAT(intervals.back().second, Catch::Matchers::WithinAbs(0.62, 1e-9));
    double covered_length = 0.;
    for (const auto &[begin, end] : intervals)
        covered_length += end - begin;
    CHECK_THAT(covered_length, Catch::Matchers::WithinAbs(0.24, 1e-9));
    const auto crossing = nonplanar_segment_overlap(
        Vec3d(0., 0., 1.), Vec3d(10., 0., 1.),
        Vec3d(5., -1., 1.), Vec3d(5., 1., 1.), 0.2, 0.05);
    REQUIRE(crossing.size() == 1);
    CHECK_THAT(crossing.front().first, Catch::Matchers::WithinAbs(0.48, 1e-9));
    CHECK_THAT(crossing.front().second, Catch::Matchers::WithinAbs(0.52, 1e-9));
    const auto vertical = nonplanar_segment_overlap(
        Vec3d(0., 0., 1.), Vec3d(10., 0., 1.),
        Vec3d(5., 0., 0.5), Vec3d(5., 0., 1.5), 0.2, 0.05);
    REQUIRE(vertical.size() == 1);
    CHECK_THAT(vertical.front().first, Catch::Matchers::WithinAbs(0.48, 1e-9));
    CHECK_THAT(vertical.front().second, Catch::Matchers::WithinAbs(0.52, 1e-9));
}

TEST_CASE("Non-planar clearance detects a raised wall near a roof endpoint",
          "[NonplanarSurface][NonplanarOverlap][Regression]")
{
    // Reduced from the emitted Benchy roof: the centerlines do not cross,
    // but the retained brick endpoint is only 0.055 mm away in XY and is
    // above the draped wall. A crossing-only test misses this collision.
    const bool reversed = GENERATE(false, true);
    Vec3d a(137.404, 101.022, 38.6), b(137.413, 101.036, 38.6);
    const Vec3d c(137.424, 101.09, 38.527), d(137.485, 101.176, 38.529);
    if (reversed)
        std::swap(a, b);
    const Polyline retained({Point::new_scale(a.x(), a.y()), Point::new_scale(b.x(), b.y())});
    const Polyline roof({Point::new_scale(c.x(), c.y()), Point::new_scale(d.x(), d.y())});
    const float half_width = float(scale_(0.38292 / 2.));
    const Polygons blunt_coverage = offset(offset(roof, half_width), float(scale_(0.05)));
    const Polygons rounded_coverage = offset(offset(roof, half_width,
        ClipperLib::jtRound, scale_(0.005), ClipperLib::etOpenRound), float(scale_(0.05)));
    CHECK(intersection_pl(retained, blunt_coverage).empty());
    CHECK_FALSE(intersection_pl(retained, rounded_coverage).empty());
    const auto overlap = [&](Vec3d start, Vec3d end) {
        // The clipping band for a 0.2 mm course is 0.05--0.18 mm
        // above the non-planar path; a course below it is support.
        start.z() -= 0.115;
        end.z() -= 0.115;
        return nonplanar_segment_overlap(start, end, c, d, 0.06, 0.065);
    };
    REQUIRE_FALSE(overlap(a, b).empty());
    a.z() -= 0.2;
    b.z() -= 0.2;
    CHECK(overlap(a, b).empty());
}

TEST_CASE("Sloped bead contact satisfies XY and Z limits at the same support point",
          "[NonplanarSurface][ContactIntervals][Regression]")
{
    // The XY-nearest point is too high. A point 0.15 mm along the
    // supporting slope meets the Z band and remains within the bead radius.
    const Vec3d a(0.4, 0., 0.2), b(0.6, 0., 0.2);
    for (bool reverse : {false, true}) {
        const Vec3d c = reverse ? Vec3d(1., 0., 1.) : Vec3d(0., 0., 0.);
        const Vec3d d = reverse ? Vec3d(0., 0., 0.) : Vec3d(1., 0., 1.);
        const auto contact = nonplanar_segment_overlap(a, b, c, d, 0.2, 0.05);
        REQUIRE(contact.size() == 1);
        CHECK_THAT(contact.front().first, Catch::Matchers::WithinAbs(0., 1e-9));
        CHECK_THAT(contact.front().second, Catch::Matchers::WithinAbs(0.25, 1e-9));
        CHECK(nonplanar_segment_overlap(a, b, c, d, 0.1, 0.05).empty());
    }
}

TEST_CASE("Bead contact intervals agree with independent support sampling",
          "[NonplanarSurface][ContactIntervals][Regression]")
{
    const std::array<std::pair<Vec3d, Vec3d>, 4> supports{{
        {Vec3d(0., 0.1, 0.), Vec3d(1., -0.1, 1.)},
        {Vec3d(1., -0.1, 1.), Vec3d(0., 0.1, 0.)},
        {Vec3d(0.5, 0., 0.), Vec3d(0.5, 0., 1.)},
        {Vec3d(0., 0., 0.4), Vec3d(1., 0., 0.4)}}};
    size_t comparisons = 0;
    for (const auto &[c, d] : supports) {
        const Vec3d a(0., 0., 0.3), b(1., 0.25, 0.65);
        constexpr double xy_limit = 0.2, z_limit = 0.05;
        const auto intervals = nonplanar_segment_overlap(a, b, c, d, xy_limit, z_limit);
        constexpr size_t samples = 1024;
        // Lipschitz bound on the grid's possible missed minimum. Only test
        // definite hits/misses, rather than rounding uncertain border cases.
        const double sampling_error = std::max((d - c).head<2>().norm() / xy_limit,
                                               std::abs(d.z() - c.z()) / z_limit) / samples;
        for (size_t i = 0; i < 100; ++i) {
            const double t = (double(i) + 0.5) / 100.;
            const Vec3d point = a + t * (b - a);
            double margin = std::numeric_limits<double>::max();
            for (size_t j = 0; j <= samples; ++j) {
                const Vec3d delta = point - c - (double(j) / samples) * (d - c);
                margin = std::min(margin, std::max(delta.head<2>().norm() / xy_limit,
                                                  std::abs(delta.z()) / z_limit) - 1.);
            }
            const bool reported = std::any_of(intervals.begin(), intervals.end(),
                [t](const auto &interval) { return interval.first <= t && t <= interval.second; });
            CAPTURE(t, margin, c, d);
            if (margin < -1e-8) {
                CHECK(reported);
                ++comparisons;
            } else if (margin > sampling_error + 1e-8) {
                CHECK_FALSE(reported);
                ++comparisons;
            }
        }
    }
    CHECK(comparisons > 350);
}

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
