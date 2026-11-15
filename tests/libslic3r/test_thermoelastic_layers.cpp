#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "libslic3r/ThermoelasticLayers.hpp"
#include "libslic3r/ThermoelasticSectionFootprint.hpp"
#include "libslic3r/ThermoelasticExactSections.hpp"
#include "libslic3r/ThermoelasticSectionTopology.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Line.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/libslic3r.h"
#include <cmath>
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <optional>
#include <memory>
#include <utility>
#include <vector>
#include <array>
#include <chrono>
#include <iostream>
#include <limits>

TEST_CASE("Grouped planar wall strokes report representation and cost differences",
    "[ThermoelasticLayers][.GroupedWallComparison]")
{
    using namespace Slic3r;
    using namespace Slic3r::Thermoelastic;
    const bool closed = GENERATE(false, true);
    const auto started = std::chrono::steady_clock::now();
    std::size_t checkpoints = 0;
    const auto check = [&] {
        ++checkpoints;
        if (std::chrono::steady_clock::now() - started > std::chrono::seconds(30))
            throw std::runtime_error("Grouped wall diagnostic time budget exceeded.");
    };
    constexpr unsigned count = 128;
    constexpr double pi = 3.14159265358979323846;
    std::vector<DepositedSegment> beads;
    const auto point = [&](unsigned i) {
        if (closed && i == count) i = 0;
        const double angle = (closed ? 2. : 1.5) * pi * i / count;
        return Eigen::Vector3d(.01 * std::cos(angle), .01 * std::sin(angle), .0002);
    };
    for (unsigned i = 0; i < count; ++i)
        beads.push_back({point(i), point(i + 1), .0004, .0002, 0});
    LayerGeometryLimits limits;
    limits.max_segments = count;
    limits.max_boundary_points = 20000;
    limits.arc_tolerance_m = 1e-6;
    auto exact_limits = limits;
    exact_limits.retain_exact_interface_sections = true;
    exact_limits.exact_per_command_composition = true;
    exact_limits.ownership_boundary_grid_units = 4;
    exact_limits.max_slab_group_visits = 4000000;
    const auto baseline_start = std::chrono::steady_clock::now();
    const auto baseline = section_deposition(beads, .0001, 0., exact_limits, check);
    const double baseline_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - baseline_start).count();
    const auto baseline_checks = checkpoints;
    const auto grouped_start = std::chrono::steady_clock::now();
    const auto grouped = layered_deposition(beads, 0., limits, check);
    const double grouped_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - grouped_start).count();
    const auto grouped_checks = checkpoints - baseline_checks;
    REQUIRE(baseline.regions.size() == 1);
    REQUIRE(grouped.slabs.size() == 1);
    // Native output is a separate diagnostic proxy for boundary sampling only;
    // it never replaces the exact baseline in the XOR comparison or timing.
    const auto native_proxy = section_deposition(beads, .0001, 0., limits, check);
    REQUIRE(native_proxy.regions.size() == 1);
    const auto& a = native_proxy.regions.front().footprint.native();
    const auto& b = grouped.slabs.front().footprint.native();
    const auto& exact_a = baseline.regions.front().footprint.exact();
    const auto exact_b = make_exact_section(b, 4000000, check);
    const auto difference = boolean_exact_sections(exact_a, exact_b.section,
        ExactSectionOperation::SymmetricDifference, 4000000, check);
    const auto area = exact_section_horizontal_area(difference.section,
        {-20000000,-20000000,0}, {40000000,40000000,200000}, 100000, 4000000, check);
    const auto boundary = [](const ExPolygons& polygons) {
        Linesf result;
        const auto append = [&](const Polygon& polygon) {
            for (std::size_t i = 0; i < polygon.points.size(); ++i)
                result.emplace_back(polygon.points[i].cast<double>() * 1e-9,
                    polygon.points[(i + 1) % polygon.points.size()].cast<double>() * 1e-9);
        };
        for (const auto& polygon : polygons) {
            append(polygon.contour);
            for (const auto& hole : polygon.holes) append(hole);
        }
        return result;
    };
    const auto edges_a = boundary(a), edges_b = boundary(b);
    REQUIRE_FALSE(edges_a.empty());
    REQUIRE_FALSE(edges_b.empty());
    std::size_t distance_work = 0;
    const auto directed_bound = [&](const Linesf& from, const Linesf& to) {
        std::array<double, 2> bound{};
        for (const auto& edge : from) {
            const double length = (edge.b - edge.a).norm();
            const auto divisions = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(length / .00005)));
            double local_max = 0.;
            for (std::size_t i = 0; i <= divisions; ++i) {
                check();
                const Vec2d query = edge.a + (edge.b - edge.a) * (double(i) / divisions);
                double nearest = std::numeric_limits<double>::infinity();
                for (const auto& candidate : to) {
                    if (++distance_work > 20000000) throw std::runtime_error("Grouped wall distance budget exceeded.");
                    nearest = std::min(nearest, line_alg::distance_to_squared(candidate, query));
                }
                local_max = std::max(local_max, std::sqrt(nearest));
            }
            bound[0] = std::max(bound[0], local_max);
            // Distance to a closed set is 1-Lipschitz. Every edge point lies
            // within half a sample interval of an evaluated boundary point.
            bound[1] = std::max(bound[1], local_max + length / (2. * divisions));
        }
        return bound;
    };
    const auto ab = directed_bound(edges_a, edges_b), ba = directed_bound(edges_b, edges_a);
    const auto holes = [](const ExPolygons& polygons) {
        std::size_t count_holes = 0;
        for (const auto& polygon : polygons) count_holes += polygon.holes.size();
        return count_holes;
    };
    std::cout << "Grouped wall approximation diagnostic baseline=exact_per_command closed=" << closed
              << " baseline_s=" << baseline_seconds << " grouped_s=" << grouped_seconds
              << " baseline_checkpoints=" << baseline_checks << " grouped_checkpoints=" << grouped_checks
              << " (checkpoints are not comparable Boolean operation counts)"
              << " exact_xor_area_m2=" << area.value << " exact_xor_bounds=" << area.bounds[0] << ',' << area.bounds[1]
              << " baseline_sampler_work=" << baseline.ownership_work
              << " exact_comparison_work=" << exact_b.work + difference.work + area.work
              << " native_proxy_boundary_distance_sample_lower_m=" << std::max(ab[0],ba[0])
              << " native_proxy_boundary_distance_lipschitz_upper_m=" << std::max(ab[1],ba[1])
              << " distance_pairs=" << distance_work << " native_proxy_components=" << a.size() << ',' << b.size()
              << " native_proxy_holes=" << holes(a) << ',' << holes(b) << '\n';
    CHECK(a.size() == b.size());
    CHECK(holes(a) == holes(b));
    CHECK(std::isfinite(area.value));
    auto capped = limits;
    capped.max_segments = count - 1;
    CHECK_THROWS_AS(layered_deposition(beads, 0., capped), ResourceLimit);
    CHECK_THROWS_WITH(layered_deposition(beads, 0., limits,
        [] { throw std::runtime_error("cancel grouped wall"); }), "cancel grouped wall");
    check();
}

using namespace Slic3r::Thermoelastic;
using Catch::Matchers::WithinRel;

TEST_CASE("Exact sampler work categories partition existing charges", "[ThermoelasticLayers]")
{
    const DepositedSegment bead{{0.,0.,.0002},{.001,0.,.0002},.0004,.0002,0};
    LayerGeometryLimits limits;
    const auto native=section_deposition({bead},.0001,0.,limits);
    REQUIRE(native.regions.size()==1);
    CHECK(native.sampler_work.exact_construction==0);
    CHECK(native.sampler_work.exact_composition==0);
    CHECK(native.sampler_work.other==native.ownership_work);
    const auto constructed=make_exact_section(native.regions.front().footprint.native());
    limits.retain_exact_interface_sections=true;
    limits.ownership_boundary_grid_units=4;
    limits.exact_per_command_composition=true;
    for (const unsigned count:{1u,2u}) {
        const auto exact=section_deposition(std::vector<DepositedSegment>(count,bead),.0001,0.,limits);
        CHECK(exact.sampler_work.exact_construction==count*(native.boundary_points+1+constructed.work));
        CHECK(exact.sampler_work.exact_construction+exact.sampler_work.exact_composition+
            exact.sampler_work.other==exact.ownership_work);
        CHECK(exact.sampler_work.other>0);
        if (count==2) CHECK(exact.sampler_work.exact_composition>0);
        auto bounded=limits;
        bounded.max_slab_group_visits=exact.ownership_work;
        CHECK(section_deposition(std::vector<DepositedSegment>(count,bead),.0001,0.,bounded).ownership_work==exact.ownership_work);
        --bounded.max_slab_group_visits;
        CHECK_THROWS_AS(section_deposition(std::vector<DepositedSegment>(count,bead),.0001,0.,bounded),ResourceLimit);
    }
    CHECK_THROWS_WITH(section_deposition({bead},.0001,0.,limits,
        [] { throw std::runtime_error("cancel sampler work categories"); }),"cancel sampler work categories");
}

TEST_CASE("Sampler overlay diagnostics preserve composition and admission", "[ThermoelasticLayers]")
{
    const std::vector<DepositedSegment> beads{
        {{0.,0.,.0002},{.001,0.,.0002},.0004,.0002,0},
        {{0.,.0001,.0002},{.001,.0001,.0002},.0004,.0002,0}};
    LayerGeometryLimits limits;
    limits.retain_exact_interface_sections=true;
    limits.exact_per_command_composition=true;
    limits.ownership_boundary_grid_units=4;
    const auto baseline=section_deposition(beads,.0001,0.,limits);
    ExactOverlayDiagnostics counters;
    const auto measured=section_deposition(beads,.0001,0.,limits,{},nullptr,nullptr,&counters);
    REQUIRE(measured.regions.size()==1);
    REQUIRE(baseline.regions.size()==1);
    CHECK(measured.ownership_work==baseline.ownership_work);
    CHECK(measured.sampler_work.exact_composition==baseline.sampler_work.exact_composition);
    CHECK(counters.overlays>0);
    CHECK(counters.boxes_projection+counters.sorting+counters.sweep+counters.active_scans+
        counters.event_reservation<=measured.sampler_work.exact_composition);
    CHECK(boolean_exact_sections(measured.regions.front().footprint.exact(),baseline.regions.front().footprint.exact(),
        ExactSectionOperation::SymmetricDifference).section.empty());
}

TEST_CASE("Balanced exact command composition preserves prefixes with less growing-union work", "[ThermoelasticLayers]")
{
    LayerGeometryLimits limits;
    limits.retain_exact_interface_sections=true;
    limits.exact_per_command_composition=true;
    limits.ownership_boundary_grid_units=4;
    limits.max_slab_group_visits=8000000;
    std::vector<DepositedSegment> beads;
    ExactSection left_fold;
    std::size_t left_work=0;
    for (std::size_t i=0;i<64;++i) {
        const double y=double(i)*.00001;
        beads.push_back({{0.,y,.0002},{.001,y+.000333333,.0002},.0004,.0002,0});
        auto native_limits=limits;
        native_limits.exact_per_command_composition=false;
        const auto singleton=section_deposition({beads.back()},.0001,0.,native_limits);
        REQUIRE(singleton.regions.size()==1);
        const auto& polygon=singleton.regions.front().footprint.native();
        auto cut=make_exact_section(polygon);
        left_work+=singleton.boundary_points+1+cut.work;
        auto joined=boolean_exact_sections(left_fold,cut.section,ExactSectionOperation::Union);
        left_work+=joined.work;
        left_fold=std::move(joined.section);
        const auto composed=section_deposition(beads,.0001,0.,limits);
        REQUIRE(composed.regions.size()==1);
        CHECK(boolean_exact_sections(composed.regions.front().footprint.exact(),left_fold,
            ExactSectionOperation::SymmetricDifference).section.empty());
        CHECK_FALSE(composed.deposited_partition);
        CHECK_FALSE(composed.occupied_partition);
        if (i==63) CHECK(composed.ownership_work<left_work);
    }
    const auto complete=section_deposition(beads,.0001,0.,limits);
    auto bounded=limits;
    bounded.max_slab_group_visits=complete.ownership_work-1;
    CHECK_THROWS_AS(section_deposition(beads,.0001,0.,bounded),ResourceLimit);
    bounded=limits;
    bounded.max_boundary_points=4;
    CHECK_THROWS_AS(section_deposition(beads,.0001,0.,bounded),ResourceLimit);
    CHECK_THROWS_AS(section_deposition(beads,.0001,0.,limits,
        [] { throw std::runtime_error("cancel balanced union"); }),std::runtime_error);
}

TEST_CASE("Balanced exact composition collapses overlapping runs before admitting a cut under pressure", "[ThermoelasticLayers]")
{
    const DepositedSegment bead{{0.,0.,.0002},{.001,.000333333,.0002},.0004,.0002,0};
    LayerGeometryLimits limits;
    limits.retain_exact_interface_sections=true;
    limits.exact_per_command_composition=true;
    limits.ownership_boundary_grid_units=4;
    const auto singleton=section_deposition({bead},.0001,0.,limits);
    REQUIRE(singleton.regions.size()==1);
    const auto& expected=singleton.regions.front().footprint.exact();
    const auto count=expected.vertex_count();
    REQUIRE(count>0);
    auto native_limits=limits;
    native_limits.exact_per_command_composition=false;
    REQUIRE(section_deposition({bead},.0001,0.,native_limits).boundary_points==count);
    limits.max_boundary_points=2*count;
    REQUIRE(limits.max_boundary_points<3*count);
    // After three identical commands, binary levels 0 and 1 each retain this
    // same count (aliases count separately). The fourth cut would need 3*count:
    // existing runs must collapse BEFORE it can enter the two-cut ceiling.
    const auto result=section_deposition(std::vector<DepositedSegment>(4,bead),.0001,0.,limits);
    REQUIRE(result.regions.size()==1);
    CHECK(result.boundary_points==count);
    CHECK(boolean_exact_sections(result.regions.front().footprint.exact(),expected,
        ExactSectionOperation::SymmetricDifference).section.empty());
    CHECK_FALSE(result.deposited_partition);
    CHECK_FALSE(result.occupied_partition);
    // A ceiling smaller than two cuts cannot admit even the second command;
    // consolidation does not waive the incoming-plus-retained payload limit.
    limits.max_boundary_points=2*count-1;
    CHECK_THROWS_AS(section_deposition(std::vector<DepositedSegment>(2,bead),.0001,0.,limits),ResourceLimit);
}

TEST_CASE("Frozen cut admission preserves construction limits before clipping simplifies a tapered bead", "[ThermoelasticLayers]")
{
    DepositedSegment bead{{0.,0.,.0002},{.001,0.,.0002},.001,.0002,0};
    bead.end_height_m=.0001;
    LayerGeometryLimits limits;
    limits.arc_tolerance_m=bead.width_m*1e-6;
    limits.retain_exact_interface_sections=true;
    limits.ownership_boundary_grid_units=4;
    const Slic3r::ExPolygon clip(Slic3r::Polygon{{450000,-100},{550000,-100},{550000,100},{450000,100}});
    // The unbounded native construction simplifies to four clipped points,
    // but its intermediate tapered caps require more than1,024 points.
    const auto native=section_deposition({bead},.0001,0.,limits,{},nullptr,&clip);
    REQUIRE(native.regions.size()==1);
    REQUIRE(native.boundary_points==4);
    bool resolved=false;
    SectionOwnership owner{.0001,{}};
    owner.resolve_deposited_partition=[&](const SampledSectionConflict&) {
        resolved=true;
        return std::shared_ptr<const AuthenticatedDepositedPartition>{};
    };
    CHECK_THROWS_WITH(section_deposition({bead},.0001,0.,limits,{},&owner,&clip),
        "Tapered bead section exceeds the boundary-point budget.");
    CHECK_FALSE(resolved);
}

TEST_CASE("Canonical section footprints retain exact handles and reject native fallback", "[ThermoelasticLayers]")
{
    const Slic3r::ExPolygons polygons{Slic3r::ExPolygon(Slic3r::Polygon{{0,0},{5,0},{0,1}})};
    const auto source=make_exact_section(polygons);
    const auto clip=make_exact_section({Slic3r::ExPolygon(Slic3r::Polygon{{4,0},{6,0},{6,1},{4,1}})});
    const auto rational=boolean_exact_sections(source.section,clip.section,ExactSectionOperation::Intersection);
    SectionFootprint footprint(rational.section);
    REQUIRE(footprint.is_exact());
    CHECK_FALSE(footprint.empty());
    CHECK(footprint.vertex_count()==3);
    REQUIRE(footprint.bounds().has_value());
    CHECK(footprint.bounds()->min==Slic3r::Point(4,0));
    CHECK(footprint.bounds()->max==Slic3r::Point(5,1));
    CHECK_THROWS_AS(footprint.native(),std::invalid_argument);
    const SectionFootprint& immutable=footprint;
    CHECK_THROWS_AS(immutable.native(),std::invalid_argument);
    const DepositedSection::Region region{2,footprint};
    const DepositedSlab slab{0.,.001,region.material,region.footprint,0};
    footprint=SectionFootprint{};
    CHECK_FALSE(footprint.is_exact());
    CHECK(footprint.native().empty());
    CHECK_THROWS_AS(footprint.exact(),std::invalid_argument);
    REQUIRE(slab.footprint.is_exact());
    const std::array<coord_t,3> origin{0,0,0},size{10,10,10};
    const auto copied=exact_section_volume_moments(slab.footprint.exact(),origin,size,0,10);
    CHECK_THAT(copied.values[0],WithinRel(1e-27,1e-13));
    const auto topology=qualify_section_topology({slab},0.);
    REQUIRE(topology.pieces().size()==1);
    CHECK(topology.pieces().front().slab.footprint.is_exact());

    const SectionFootprint empty_exact(ExactSection{});
    REQUIRE(empty_exact.is_exact());
    CHECK(empty_exact.empty());
    CHECK(empty_exact.vertex_count()==0);
    CHECK_FALSE(empty_exact.bounds().has_value());
    CHECK_THROWS_AS(empty_exact.native(),std::invalid_argument);
    CHECK_THROWS_WITH(qualify_section_topology({{0.,.001,2,empty_exact,0}},0.),
        "Section topology has no occupied solid.");
    const SectionFootprint native(polygons);
    CHECK_FALSE(native.is_exact());
    CHECK(native.native()==polygons);
    CHECK(native.vertex_count()==3);
}


TEST_CASE("Half-native query planes sample a one-nanometre bead without changing native geometry", "[ThermoelasticLayers]")
{
    const DepositedSegment bead{{0.,0.,1e-9},{.001,0.,1e-9},.0004,1e-9,0};
    LayerGeometryLimits half;
    half.half_native_section_planes = true;
    CHECK(section_deposition({bead}, .5e-9, 0.).regions.empty());
    const auto interior = section_deposition({bead}, .5e-9, 0., half);
    REQUIRE(interior.area_m2 > 0.);
    CHECK_THAT(interior.area_m2, WithinRel(section_deposition({bead}, 0., 0.).area_m2, 1e-12));
    CHECK(section_deposition({bead}, 1e-9, 0., half).regions.empty());
    SectionSweepLimits sweep;
    sweep.geometry = half;
    const std::vector<double> planes{0., .5e-9, 1e-9};
    std::size_t seen = 0;
    CHECK(visit_deposition_sections({bead}, planes, 0., [&](std::size_t i, const auto& section) {
        ++seen;
        CHECK_THAT(section.area_m2, Catch::Matchers::WithinAbs(section_deposition({bead}, planes[i], 0., half).area_m2, 1e-30));
    }, sweep) == 2);
    CHECK(seen == planes.size());
    CHECK_THROWS_AS(visit_deposition_sections({bead}, planes, 0., [](std::size_t, const auto&) {}), std::invalid_argument);
}

TEST_CASE("Half-native sloped and tapered sections agree with their sweep and retain integer-plane results", "[ThermoelasticLayers]")
{
    LayerGeometryLimits half;
    half.half_native_section_planes = true;
    half.arc_tolerance_m = 1e-9;
    auto integral = half;
    integral.half_native_section_planes = false;
    for (const auto shape : {BeadSection::RectangularEnvelope, BeadSection::RoundedRectangle})
        for (const bool tapered : {false, true}) {
            DepositedSegment bead{{0.,0.,2e-9},{1e-6,0.,4e-9},4e-7,2e-9,0,shape};
            if (tapered) bead.end_height_m = 0.;
            SectionSweepLimits sweep;
            sweep.geometry = half;
            const std::vector<double> planes{.5e-9,1e-9,1.5e-9,2e-9,2.5e-9,3e-9,3.5e-9,4e-9};
            visit_deposition_sections({bead}, planes, 0., [&](std::size_t i, const auto& section) {
                CHECK_THAT(section.area_m2, Catch::Matchers::WithinAbs(section_deposition({bead}, planes[i], 0., half).area_m2, 1e-30));
                if (i % 2 == 1)
                    CHECK_THAT(section.area_m2, Catch::Matchers::WithinAbs(section_deposition({bead}, planes[i], 0., integral).area_m2, 1e-30));
            }, sweep);
        }
}

namespace {
Slic3r::ExPolygons ownership_box(coord_t left, coord_t right)
{
    return {Slic3r::ExPolygon(Slic3r::Polygon{{left,-1000000},{right,-1000000},
        {right,11000000},{left,11000000}})};
}
}

TEST_CASE("Native ownership partitions conflicts without filling holes or clipping unambiguous beads", "[ThermoelasticLayers]")
{
    std::vector<DepositedSegment> ring{
        {{0.,0.,.0002},{.01,0.,.0002},.0004,.0002,0},
        {{.01,0.,.0002},{.01,.01,.0002},.0004,.0002,0},
        {{.01,.01,.0002},{0.,.01,.0002},.0004,.0002,0},
        {{0.,.01,.0002},{0.,0.,.0002},.0004,.0002,0}};
    const auto reference = section_deposition(ring,.0001,0.);
    const auto base = ring;
    for (auto bead : base) { bead.material = 1; ring.push_back(bead); }
    CHECK_THROWS_AS(section_deposition(ring,.0001,0.),LayerMaterialOverlap);
    SectionOwnership ownership{.0001,{{0,ownership_box(-1000000,5000000)},
        {1,ownership_box(5000000,11000000)}}};
    const auto resolved = section_deposition(ring,.0001,0.,{},{},&ownership);
    CHECK_THAT(resolved.resolved_overlap_area_m2,WithinRel(reference.area_m2,1e-10));
    CHECK_THAT(resolved.area_m2,WithinRel(reference.area_m2,1e-10));
    Slic3r::ExPolygons combined;
    for (const auto& region : resolved.regions)
        combined.insert(combined.end(),region.footprint.native().begin(),region.footprint.native().end());
    combined = Slic3r::union_ex(combined);
    REQUIRE(combined.size() == 1);
    CHECK(combined.front().holes.size() == 1);
    CHECK(Slic3r::diff_ex(combined,reference.regions.front().footprint.native()).empty());
    CHECK(Slic3r::intersection_ex(resolved.regions[0].footprint.native(),resolved.regions[1].footprint.native()).empty());
    ring.push_back({{.02,0.,.0002},{.03,0.,.0002},.0004,.0002,0});
    const auto extra = section_deposition({ring.back()},.0001,0.);
    const auto retained = section_deposition(ring,.0001,0.,{},{},&ownership);
    CHECK_THAT(retained.area_m2,WithinRel(reference.area_m2 + extra.area_m2,1e-10));
}

TEST_CASE("Section ownership rejects ambiguous missing and unsupported assignments", "[ThermoelasticLayers]")
{
    const std::vector<DepositedSegment> beads{
        {{0.,0.,.0002},{.01,0.,.0002},.0004,.0002,0},
        {{0.,0.,.0002},{.01,0.,.0002},.0004,.0002,1}};
    SectionOwnership owner{.0001,{{0,ownership_box(-1000000,11000000)}}};
    CHECK_NOTHROW(section_deposition(beads,.0001,0.,{},{},&owner));
    owner.regions.push_back({1,ownership_box(-1000000,11000000)});
    CHECK_THROWS_AS(section_deposition(beads,.0001,0.,{},{},&owner),std::invalid_argument);
    LayerGeometryLimits exact_limits;
    exact_limits.ownership_boundary_grid_units = 4;
    CHECK_THROWS_AS(section_deposition(beads,.0001,0.,exact_limits,{},&owner),std::invalid_argument);
    owner.regions.resize(1);
    owner.regions[0].material = 2;
    CHECK_THROWS_AS(section_deposition(beads,.0001,0.,{},{},&owner),std::invalid_argument);
    CHECK_THROWS_AS(section_deposition(beads,.0001,0.,exact_limits,{},&owner),std::invalid_argument);
    auto unrelated = beads;
    unrelated.push_back({{.02,0.,.0002},{.03,0.,.0002},.0004,.0002,2});
    CHECK_THROWS_AS(section_deposition(unrelated,.0001,0.,{},{},&owner),std::invalid_argument);
    CHECK_THROWS_AS(section_deposition(unrelated,.0001,0.,exact_limits,{},&owner),std::invalid_argument);
    owner.regions[0].material = 0;
    owner.regions[0].footprint.native() = ownership_box(-1000000,5000000);
    CHECK_THROWS_AS(section_deposition(beads,.0001,0.,{},{},&owner),std::invalid_argument);
    CHECK_THROWS_AS(section_deposition(beads,.0001,0.,exact_limits,{},&owner),std::invalid_argument);
    owner.regions[0].footprint.native() = ownership_box(-1000000,11000000);
    owner.z_m = .00015;
    CHECK_THROWS_AS(section_deposition(beads,.0001,0.,{},{},&owner),std::invalid_argument);
    owner.z_m = .0001;
    LayerGeometryLimits limits;
    limits.max_slab_group_visits = 1;
    CHECK_THROWS_AS(section_deposition(beads,.0001,0.,limits,{},&owner),ResourceLimit);
    CHECK_THROWS_AS(section_deposition(beads,.0001,0.,{},[]{throw std::runtime_error("cancel");},&owner),std::runtime_error);
    std::size_t callbacks = 0;
    visit_deposition_sections(beads,{.0001},0.,[&](std::size_t,const auto& section) {
        CHECK(section.resolved_overlap_area_m2 > 0.);
        ++callbacks;
    },{},{},[&](std::size_t,double) { return std::optional<SectionOwnership>{owner}; });
    CHECK(callbacks == 1);
}

TEST_CASE("Authenticated interior material ownership is not exterior boundary error", "[ThermoelasticLayers]")
{
    const std::vector<DepositedSegment> beads{
        {{0., 0., .0002}, {.01, 0., .0002}, .004, .0002, 0},
        {{0., 0., .0002}, {.01, 0., .0002}, .004, .0002, 1}};
    const Slic3r::ExPolygons inner{Slic3r::ExPolygon(Slic3r::Polygon{
        {4000000, -1000000}, {6000000, -1000000}, {6000000, 1000000}, {4000000, 1000000}})};
    const auto outer = section_deposition({beads.front()}, .0001, 0.).regions.front().footprint.native();
    SectionOwnership owner{.0001, {{0, inner}, {1, Slic3r::diff_ex(outer, inner)}}};
    LayerGeometryLimits limits;
    limits.ownership_boundary_grid_units = 4;
    const auto result = section_deposition(beads, .0001, 0., limits, {}, &owner);
    REQUIRE(result.regions.size() == 2);
    CHECK(Slic3r::diff_ex(inner, result.regions[0].footprint.native()).empty());
    CHECK(Slic3r::intersection_ex(inner, result.regions[1].footprint.native()).empty());
    CHECK(Slic3r::diff_ex(outer, Slic3r::union_ex(result.regions[0].footprint.native(), result.regions[1].footprint.native())).empty());
    CHECK(result.boundary_error.uncertainty_nm == 0);
}

TEST_CASE("Streamed deposition sections match independent sections without whole-print rescans", "[ThermoelasticLayers]")
{
    // Deliberately unsorted, overlapping height ranges; separate materials do
    // not touch. Include a sloped bead and planes outside the occupied body.
    const std::vector<DepositedSegment> beads{
        {{0.,0.,.0006},{.01,0.,.0006},.0004,.0002,0},
        {{0.,.002,.0002},{.01,.002,.0006},.0004,.0002,1},
        {{0.,0.,.0002},{.01,0.,.0002},.0004,.0002,0}};
    const std::vector<double> planes{0.,.0001,.0002,.0003,.0004,.0005,.0006,.0007};
    std::size_t callbacks = 0;
    const auto visits = visit_deposition_sections(beads, planes, 0., [&](std::size_t i, const auto &section) {
        const auto reference = section_deposition(beads, planes[i], 0.);
        REQUIRE(section.regions.size() == reference.regions.size());
        CHECK(section.boundary_points == reference.boundary_points);
        CHECK_THAT(section.area_m2, Catch::Matchers::WithinAbs(reference.area_m2, 1e-18));
        for (std::size_t j = 0; j < section.regions.size(); ++j) {
            CHECK(section.regions[j].material == reference.regions[j].material);
            CHECK(section.regions[j].footprint.native() == reference.regions[j].footprint.native());
        }
        ++callbacks;
    });
    CHECK(callbacks == planes.size());
    CHECK(visits == 10);
    CHECK(visits < beads.size() * planes.size());
}

TEST_CASE("Section sweeps preserve rounded bead voids and reject partial work as completion", "[ThermoelasticLayers]")
{
    const DepositedSegment bead{{0.,0.,.0002},{.01,0.,.0002},.0004,.0002,0,
        BeadSection::RoundedRectangle};
    const std::vector<double> planes{.00001,.0001,.00019};
    const auto callback = [&](std::size_t i, const auto &section) {
        const auto reference = section_deposition({bead}, planes[i], 0.);
        CHECK_THAT(section.area_m2, WithinRel(reference.area_m2, 1e-12));
    };
    CHECK(visit_deposition_sections({bead}, planes, 0., callback) == 3);
    SectionSweepLimits extended;
    extended.geometry.max_slab_group_visits = 8000000;
    CHECK(visit_deposition_sections({bead}, planes, 0., callback, extended) == 3);
    extended.geometry.max_slab_group_visits = 8000001;
    CHECK_THROWS_AS(visit_deposition_sections({bead}, planes, 0., callback, extended), std::invalid_argument);
    SectionSweepLimits limits;
    limits.max_segment_visits = 2;
    CHECK_THROWS_AS(visit_deposition_sections({bead}, planes, 0., callback, limits), ResourceLimit);
    limits = {};
    limits.max_sections = 2;
    CHECK_THROWS_AS(visit_deposition_sections({bead}, planes, 0., callback, limits), ResourceLimit);
    CHECK_THROWS_AS(visit_deposition_sections({bead}, {.0001,.0001}, 0., callback), std::invalid_argument);
    CHECK_THROWS_AS(visit_deposition_sections({bead}, {.0001,.00005}, 0., callback), std::invalid_argument);
    CHECK_THROWS_AS(visit_deposition_sections({bead}, planes, 0., {}), std::invalid_argument);
    CHECK_THROWS_AS(visit_deposition_sections({bead}, planes, 0., callback, {}, [] {
        throw std::runtime_error("cancel");
    }), std::runtime_error);
    CHECK_THROWS_AS(visit_deposition_sections({bead}, planes, 0., [](auto, const auto &) {
        throw std::runtime_error("consumer cancelled");
    }), std::runtime_error);
    auto other = bead;
    other.material = 1;
    CHECK_THROWS_AS(visit_deposition_sections({bead, other}, planes, 0., callback), LayerMaterialOverlap);
    // Invalid off-plane geometry is not hidden by the spatial index.
    other.start_m.z() = other.end_m.z() = .1;
    other.width_m = -1.;
    CHECK_THROWS_AS(visit_deposition_sections({bead, other}, planes, 0., callback), std::invalid_argument);
}

TEST_CASE("Resumable section sweeps preserve plane order and poison incomplete or reentered work", "[ThermoelasticLayers]")
{
    const std::vector<DepositedSegment> beads{{{0.,0.,.0002},{.01,0.,.0002},.0004,.0002,0,BeadSection::RoundedRectangle}};
    const std::vector<double> planes{.00001,.0001,.00019};
    std::size_t received=0;
    const auto consume=[&](std::size_t i,const DepositedSection& section) {
        CHECK(i==received++);
        const auto expected=section_deposition(beads,planes[i],0.);
        REQUIRE(section.regions.size()==expected.regions.size());
        for (std::size_t r=0;r<section.regions.size();++r)
            CHECK(section.regions[r].footprint.native()==expected.regions[r].footprint.native());
    };
    DepositionSectionCursor cursor(beads,planes,0.);
    for (std::size_t i=0;i<planes.size();++i) {
        CHECK(cursor.advance(consume));
        CHECK(received==i+1);
    }
    CHECK_FALSE(cursor.advance(consume));
    CHECK(cursor.completed_visits()==3);
    DepositionSectionCursor incomplete(beads,planes,0.);
    CHECK_THROWS_AS(incomplete.completed_visits(),std::logic_error);
    CHECK_THROWS_AS(incomplete.advance(consume),std::logic_error);
    DepositionSectionCursor canceled(beads,planes,0.);
    CHECK_THROWS_AS(canceled.advance([](auto,const auto&) { throw std::runtime_error("cancel"); }),std::runtime_error);
    CHECK_THROWS_AS(canceled.completed_visits(),std::logic_error);
    DepositionSectionCursor reentered(beads,planes,0.);
    CHECK_THROWS_AS(reentered.advance([&](auto,const auto&) {
        // Even a swallowed nested failure invalidates the outer callback.
        try { reentered.advance(consume); } catch (const std::logic_error&) {}
    }),std::logic_error);
    CHECK_THROWS_AS(reentered.completed_visits(),std::logic_error);
}

TEST_CASE("Layered deposited envelopes preserve gaps and union repeated paths", "[ThermoelasticLayers]")
{
    DepositedSegment first{{0.,0.,.0002},{.01,0.,.0002},.0004,.0002,0};
    auto second = first;
    second.start_m.y() = second.end_m.y() = .001;
    const auto result = layered_deposition({first, second, first}, 0.);
    REQUIRE(result.slabs.size() == 1);
    CHECK(result.slabs.front().footprint.native().size() == 2);
    const double capsule_area = .01 * .0004 + std::acos(-1.) * .0002 * .0002;
    CHECK_THAT(result.occupied_volume_m3, WithinRel(2. * capsule_area * .0002, .001));
    CHECK(result.boundary_points < 200);
}

TEST_CASE("Layered deposition retains holes and splits unequal height intervals", "[ThermoelasticLayers]")
{
    std::vector<DepositedSegment> ring{
        {{0.,0.,.0002},{.01,0.,.0002},.0004,.0002,0},
        {{.01,0.,.0002},{.01,.01,.0002},.0004,.0002,0},
        {{.01,.01,.0002},{0.,.01,.0002},.0004,.0002,0},
        {{0.,.01,.0002},{0.,0.,.0002},.0004,.0002,0}};
    const auto baseline = layered_deposition(ring, 0.);
    REQUIRE(baseline.slabs.size() == 1);
    REQUIRE(baseline.slabs.front().footprint.native().size() == 1);
    CHECK(baseline.slabs.front().footprint.native().front().holes.size() == 1);
    ring.push_back({{.02,0.,.0003},{.03,0.,.0003},.0004,.0002,1});
    const auto split = layered_deposition(ring, 0.);
    CHECK(split.slabs.size() == 4);
    const double added = (.01 * .0004 + std::acos(-1.) * .0002 * .0002) * .0002;
    CHECK_THAT(split.occupied_volume_m3, WithinRel(baseline.occupied_volume_m3 + added, .001));
}

TEST_CASE("Layered deposition rejects ambiguous geometry and bounded work", "[ThermoelasticLayers]")
{
    DepositedSegment first{{0.,0.,.0002},{.01,0.,.0002},.0004,.0002,0};
    auto other = first;
    other.material = 1;
    CHECK_THROWS_AS(layered_deposition({first, other}, 0.), std::invalid_argument);
    other = first;
    other.end_m.z() += .0001;
    CHECK_THROWS_AS(layered_deposition({other}, 0.), std::invalid_argument);
    other = first;
    other.section = BeadSection::RoundedRectangle;
    CHECK_THROWS_AS(layered_deposition({other}, 0.), std::invalid_argument);
    CHECK_THROWS_AS(layered_deposition({first}, .0001), std::invalid_argument);
    LayerGeometryLimits limits;
    limits.max_segments = 1;
    CHECK_THROWS_AS(layered_deposition({first, first}, 0., limits), ResourceLimit);
    limits.max_boundary_points = 3;
    CHECK_THROWS_AS(layered_deposition({first}, 0., limits), ResourceLimit);
    CHECK_THROWS_AS(layered_deposition({first}, 0., {}, [] { throw std::runtime_error("cancel"); }),
        std::runtime_error);
}

TEST_CASE("Layered cap geometry converges independently of path direction and translation", "[ThermoelasticLayers]")
{
    DepositedSegment bead{{-.003,-.005,.0002},{-.002,-.005,.0002},.0004,.0002,0};
    const double expected = (.001 * .0004 + std::acos(-1.) * .0002 * .0002) * .0002;
    LayerGeometryLimits limits;
    limits.arc_tolerance_m = 1e-5;
    const double coarse = layered_deposition({bead}, 0., limits).occupied_volume_m3;
    limits.arc_tolerance_m = 1e-7;
    const double fine = layered_deposition({bead}, 0., limits).occupied_volume_m3;
    CHECK(std::abs(fine - expected) < std::abs(coarse - expected));
    CHECK_THAT(fine, WithinRel(expected, .0005));
    std::swap(bead.start_m, bead.end_m);
    CHECK_THAT(layered_deposition({bead}, 0., limits).occupied_volume_m3, WithinRel(fine, 1e-9));
    bead.start_m.x() += .01;
    bead.end_m.x() += .01;
    CHECK_THAT(layered_deposition({bead}, 0., limits).occupied_volume_m3, WithinRel(fine, 1e-9));
}

TEST_CASE("Bounded polygon batches retain the union of a long subdivided bead", "[ThermoelasticLayers]")
{
    const DepositedSegment bead{{0.,0.,.0002},{.01,0.,.0002},.0004,.0002,0};
    std::vector<DepositedSegment> pieces;
    for (std::size_t i = 0; i < 1200; ++i) {
        auto piece = bead;
        piece.start_m.x() = .01 * double(i) / 1200.;
        piece.end_m.x() = .01 * double(i + 1) / 1200.;
        pieces.push_back(piece);
    }
    const auto single = layered_deposition({bead}, 0.);
    const auto split = layered_deposition(pieces, 0.);
    REQUIRE(split.slabs.size() == 1);
    CHECK(split.slabs.front().footprint.native().size() == 1);
    CHECK(split.slabs.front().footprint.native().front().holes.empty());
    CHECK_THAT(split.occupied_volume_m3, WithinRel(single.occupied_volume_m3, 1e-5));
}

TEST_CASE("Sloped bead sections retain straight cuts and original endpoint caps", "[ThermoelasticLayers]")
{
    DepositedSegment bead{{0.,0.,.0002},{.01,0.,.0006},.0004,.0002,0};
    const double half_cap = .5 * std::acos(-1.) * .0002 * .0002;
    const double end_area = .0025 * .0004 + half_cap;
    CHECK_THAT(section_deposition({bead}, .0001, 0.).area_m2, WithinRel(end_area, .001));
    CHECK_THAT(section_deposition({bead}, .0003, 0.).area_m2, WithinRel(.005 * .0004, 1e-6));
    CHECK_THAT(section_deposition({bead}, .0005, 0.).area_m2, WithinRel(end_area, .001));
    // The start cap's top is excluded; a new round cap at this cut would
    // incorrectly add material where the nozzle was already above the plane.
    CHECK_THAT(section_deposition({bead}, .0002, 0.).area_m2, WithinRel(.005 * .0004, 1e-6));
    CHECK(section_deposition({bead}, .0006, 0.).regions.empty());
    std::swap(bead.start_m, bead.end_m);
    CHECK_THAT(section_deposition({bead}, .0001, 0.).area_m2, WithinRel(end_area, .001));
    CHECK_THAT(section_deposition({bead}, .0003, 0.).area_m2, WithinRel(.005 * .0004, 1e-6));
    CHECK_THAT(section_deposition({bead}, .0002, 0.).area_m2, WithinRel(.005 * .0004, 1e-6));
}

TEST_CASE("Integrating sloped bead sections conserves nominal vertical-section volume", "[ThermoelasticLayers]")
{
    const DepositedSegment bead{{-.01,-.01,.0002},{0.,0.,.0006},.0004,.0002,0};
    double volume = 0.;
    constexpr double dz = .00001;
    for (int i = 0; i < 60; ++i)
        volume += section_deposition({bead}, (i + .5) * dz, 0.).area_m2 * dz;
    const double expected = (.01 * std::sqrt(2.) * .0004 + std::acos(-1.) * .0002 * .0002) * .0002;
    CHECK_THAT(volume, WithinRel(expected, .001));
}

TEST_CASE("Horizontal sections agree with slab footprints and enforce budgets and materials", "[ThermoelasticLayers]")
{
    DepositedSegment bead{{0.,0.,.0002},{.01,0.,.0002},.0004,.0002,0};
    const auto layers = layered_deposition({bead}, 0.);
    CHECK_THAT(section_deposition({bead}, .0001, 0.).area_m2 * .0002,
        WithinRel(layers.occupied_volume_m3, 1e-9));
    CHECK(section_deposition({bead}, .0003, 0.).regions.empty());
    auto second = bead;
    second.material = 1;
    CHECK_THROWS_AS(section_deposition({bead, second}, .0001, 0.), LayerMaterialOverlap);
    second.start_m.y() = second.end_m.y() = .001;
    CHECK(section_deposition({bead, second}, .0001, 0.).regions.size() == 2);
    LayerGeometryLimits limits;
    limits.max_boundary_points = 3;
    CHECK_THROWS_AS(section_deposition({bead}, .0001, 0., limits), ResourceLimit);
    CHECK_THROWS_AS(section_deposition({bead}, -.0001, 0.), std::invalid_argument);
    CHECK_THROWS_AS(section_deposition({bead}, .0001, 0., {}, [] { throw std::runtime_error("cancel"); }),
        std::runtime_error);
}

TEST_CASE("Rounded horizontal sections preserve corner voids and their integrated volume", "[ThermoelasticLayers]")
{
    const double w = .00045, h = .0002, length = .01, pi = std::acos(-1.);
    DepositedSegment bead{{0., 0., h}, {length, 0., h}, w, h, 0, BeadSection::RoundedRectangle};
    LayerGeometryLimits limits;
    limits.arc_tolerance_m = 1e-8;
    const auto capsule_area = [&](double width) { return length * width + pi * width * width / 4.; };
    CHECK_THAT(section_deposition({bead}, 0., 0., limits).area_m2, WithinRel(capsule_area(w - h), 1e-5));
    CHECK_THAT(section_deposition({bead}, h / 2., 0., limits).area_m2, WithinRel(capsule_area(w), 1e-5));
    CHECK(section_deposition({bead}, h, 0., limits).regions.empty());
    constexpr int samples = 1000;
    double volume = 0.;
    for (int i = 0; i < samples; ++i)
        volume += section_deposition({bead}, (i + .5) * h / samples, 0., limits).area_m2 * h / samples;
    const double flat = w - h;
    const double section_area = flat * h + pi * h * h / 4.;
    // Integrate the nominal capsule caps separately from the straight path.
    const double cap_volume = pi / 4. * (flat * flat * h + flat * pi * h * h / 2. + 2. * h * h * h / 3.);
    CHECK_THAT(volume, WithinRel(length * section_area + cap_volume, 2e-5));
    bead.width_m = h / 2.;
    CHECK_THROWS_AS(section_deposition({bead}, h / 2., 0., limits), std::invalid_argument);
}

TEST_CASE("Rounded sloped sections conserve volume under reversal and rotation", "[ThermoelasticLayers]")
{
    const double w = .00045, h = .0002, length = .01, rise = .0004, pi = std::acos(-1.);
    DepositedSegment bead{{0.,0.,h},{length,0.,h+rise},w,h,0,BeadSection::RoundedRectangle};
    LayerGeometryLimits limits;
    limits.arc_tolerance_m = 1e-8;
    const double flat = w-h;
    const double expected = length*(flat*h+pi*h*h/4.) +
        pi/4.*(flat*flat*h+flat*pi*h*h/2.+2.*h*h*h/3.);
    constexpr int samples = 2000;
    const double dz = (h+rise)/samples;
    double volume = 0.;
    for (int i = 0; i < samples; ++i)
        volume += section_deposition({bead},(i+.5)*dz,0.,limits).area_m2*dz;
    CHECK_THAT(volume, WithinRel(expected, 1e-4));
    const std::vector<double> planes{h/2.,h,h+rise/2.,rise+h/2.};
    std::size_t visited = 0;
    CHECK(visit_deposition_sections({bead},planes,0.,[&](std::size_t index,const auto &section) {
        ++visited;
        CHECK_THAT(section.area_m2,WithinRel(section_deposition({bead},planes[index],0.).area_m2,1e-10));
    }) == planes.size());
    CHECK(visited == planes.size());
    for (double z : {h/2., h, h+rise/2., rise+h/2.}) {
        const double area = section_deposition({bead},z,0.,limits).area_m2;
        auto reversed = bead;
        std::swap(reversed.start_m,reversed.end_m);
        CHECK_THAT(section_deposition({reversed},z,0.,limits).area_m2, WithinRel(area,1e-8));
        auto rotated = bead;
        rotated.end_m.x() = 0.;
        rotated.end_m.y() = length;
        CHECK_THAT(section_deposition({rotated},z,0.,limits).area_m2, WithinRel(area,1e-8));
        auto envelope = bead;
        envelope.section = BeadSection::RectangularEnvelope;
        CHECK(area <= section_deposition({envelope},z,0.,limits).area_m2);
    }
    CHECK(section_deposition({bead},h+rise,0.,limits).regions.empty());
    limits.max_boundary_points = 3;
    CHECK_THROWS_AS(section_deposition({bead},h,0.,limits),ResourceLimit);
    CHECK_THROWS_AS(section_deposition({bead},h,0.,{},[] { throw std::runtime_error("cancel"); }),std::runtime_error);
}

TEST_CASE("Rounded sloped sections handle circular beads and nearly horizontal paths", "[ThermoelasticLayers]")
{
    const double h = .0002, length = .002, pi = std::acos(-1.);
    DepositedSegment bead{{0.,0.,h},{length,0.,3.*h},h,h,0,BeadSection::RoundedRectangle};
    LayerGeometryLimits limits;
    limits.arc_tolerance_m = 1e-8;
    // With width equal to height the side section is circular and has no flat.
    double volume = 0.;
    constexpr int samples = 2000;
    const double dz = 3.*h/samples;
    for (int i = 0; i < samples; ++i)
        volume += section_deposition({bead},(i+.5)*dz,0.,limits).area_m2*dz;
    CHECK_THAT(volume,WithinRel(length*pi*h*h/4.+pi*h*h*h/6.,1e-4));
    bead.width_m = 2.*h;
    bead.end_m.z() = h;
    const auto horizontal = section_deposition({bead},h/2.,0.,limits).area_m2;
    bead.end_m.z() += 1e-9;
    CHECK_THAT(section_deposition({bead},h/2.,0.,limits).area_m2,WithinRel(horizontal,1e-5));
}

TEST_CASE("Tapered scarf sections retain their design volume and complementary interface", "[ThermoelasticLayers]")
{
    const double w = .00045, h = .0002, length = .01, pi = std::acos(-1.);
    LayerGeometryLimits limits;
    limits.arc_tolerance_m = 1e-8;
    for (const auto shape : {BeadSection::RectangularEnvelope, BeadSection::RoundedRectangle}) {
        DepositedSegment ramp{{0.,0.,0.},{length,0.,h},w,0.,0,shape,h};
        auto complement = ramp;
        complement.start_m.z() = h;
        complement.height_m = h;
        complement.end_height_m = 0.;
        complement.material = 1;
        double volume = 0.;
        constexpr unsigned samples = 1000;
        const double flat = w-h;
        const double cap = shape == BeadSection::RectangularEnvelope ? pi*w*w*h/8. :
            pi/8.*(flat*flat*h + flat*pi*h*h/2. + 2.*h*h*h/3.);
        const double expected = length*(w*h/2.-(shape == BeadSection::RoundedRectangle ?
            (1.-pi/4.)*h*h/3. : 0.)) + cap;
        for (unsigned i = 0; i < samples; ++i) {
            const double z = (double(i)+.5)*h/samples;
            const auto section = section_deposition({ramp},z,0.,limits);
            volume += section.area_m2*h/samples;
            if (i % 100 != 0) continue;
            auto reversed = ramp;
            std::swap(reversed.start_m,reversed.end_m);
            reversed.height_m = h;
            reversed.end_height_m = 0.;
            CHECK_THAT(section_deposition({reversed},z,0.,limits).area_m2,WithinRel(section.area_m2,1e-7));
            auto rotated = ramp;
            rotated.end_m.x() = 0.;
            rotated.end_m.y() = length;
            CHECK_THAT(section_deposition({rotated},z,0.,limits).area_m2,WithinRel(section.area_m2,1e-7));
            // The complementary ramps meet at a plane; they must not overlap
            // as the former path-wide constant-height preview beads did.
            CHECK_NOTHROW(section_deposition({ramp,complement},z,0.,limits));
        }
        CHECK_THAT(volume,WithinRel(expected,2e-4));
        const std::vector<double> planes{h*.1,h*.5,h*.9};
        CHECK(visit_deposition_sections({ramp},planes,0.,[&](std::size_t i,const auto &section) {
            CHECK_THAT(section.area_m2,WithinRel(section_deposition({ramp},planes[i],0.).area_m2,1e-10));
        }) == planes.size());
        CHECK_THROWS_AS(layered_deposition({ramp},0.),std::invalid_argument);
        CHECK_THROWS_AS(section_deposition({ramp},h/2.,0.,{},[] { throw std::runtime_error("cancel"); }),std::runtime_error);
        auto small = limits;
        small.max_boundary_points = 3;
        CHECK_THROWS_AS(section_deposition({ramp},h/2.,0.,small),ResourceLimit);
        ramp.end_height_m = 0.;
        CHECK_THROWS_AS(section_deposition({ramp},h/2.,0.),std::invalid_argument);
        ramp.end_height_m = -h;
        CHECK_THROWS_AS(section_deposition({ramp},h/2.,0.),std::invalid_argument);
    }
}
