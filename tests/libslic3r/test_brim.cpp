#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/BrimEarsPoint.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/libslic3r.h"
#include "libslic3r/Brim.hpp"
#include "libslic3r/ClipperUtils.hpp"

using namespace Slic3r;

TEST_CASE("Brim brushes paint square corners and erase an editable mask", "[Brim]")
{
    BrimPoint paint(Vec3f::Zero(), 5.f);
    paint.shape = BrimPoint::Shape::Square;
    const auto square = painted_brim_area({paint}, Transform3d::Identity());
    REQUIRE(square.size() == 1);
    CHECK(square.front().contains(Point::new_scale(4.5, 4.5)));
    CHECK_THAT(area(square) * SCALING_FACTOR * SCALING_FACTOR, Catch::Matchers::WithinAbs(100., 1e-6));
    BrimPoint erase(Vec3f::Zero(), 2.f);
    erase.shape = BrimPoint::Shape::EraseCircle;
    const auto cut = painted_brim_area({paint, erase}, Transform3d::Identity());
    REQUIRE(cut.size() == 1);
    CHECK(cut.front().holes.size() == 1);
    CHECK_FALSE(cut.front().contains(Point(0, 0)));
    CHECK(cut.front().contains(Point::new_scale(4.5, 4.5)));
    const auto filled_again = painted_brim_area({paint, erase, paint}, Transform3d::Identity());
    CHECK(diff_ex(square, filled_again).empty());
    CHECK(diff_ex(filled_again, square).empty());
}

TEST_CASE("Legacy brim ears remain round and brush geometry follows instance translation", "[Brim]")
{
    BrimPoint ear(Vec3f::Zero(), 5.f);
    CHECK(ear.shape == BrimPoint::Shape::Ear);
    const auto round = painted_brim_area({ear}, Transform3d::Identity());
    REQUIRE(round.size() == 1);
    CHECK_FALSE(round.front().contains(Point::new_scale(4.5, 4.5)));
    Transform3d transform = Transform3d::Identity();
    transform.translate(Vec3d(30., 20., 0.));
    const auto moved = painted_brim_area({ear}, transform);
    REQUIRE(moved.size() == 1);
    CHECK(moved.front().contains(Point::new_scale(30., 20.)));
    CHECK_FALSE(moved.front().contains(Point(0, 0)));
}

TEST_CASE("Erased brim fragments must remain anchored to the first layer", "[Brim]")
{
    BrimPoint paint(Vec3f::Zero(), 5.f);
    paint.shape = BrimPoint::Shape::Square;
    BrimPoint erase(Vec3f::Zero(), 1.f);
    erase.shape = BrimPoint::Shape::EraseSquare;
    BrimPoints operations{paint};
    for (int y = -5; y <= 5; ++y) {
        erase.pos.y() = float(y);
        operations.push_back(erase);
    }
    const auto mask = painted_brim_area(operations, Transform3d::Identity());
    REQUIRE(mask.size() == 2);
    const ExPolygons model{ExPolygon(Polygon::new_scale({Vec2d(-6., -6.), Vec2d(-2., -6.), Vec2d(-2., 6.), Vec2d(-6., 6.)}))};
    const auto attached = connected_brim_paint(mask, model);
    REQUIRE(attached.size() == 1);
    CHECK(attached.front().contains(Point::new_scale(-3., 0.)));
    CHECK_FALSE(attached.front().contains(Point::new_scale(3., 0.)));
}
