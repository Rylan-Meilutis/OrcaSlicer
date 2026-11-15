#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/MeshBoolean.hpp"
#include "libslic3r/GCode/DepositionGeometry.hpp"
#include "libslic3r/GCode/DepositionMesh.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "test_utils.hpp"
#include <STEPControl_Reader.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS.hxx>
#include <BRep_Tool.hxx>
#include <Poly_Triangulation.hxx>
#include <TopLoc_Location.hxx>
#include <gp_Pnt.hxx>
#include <boost/nowide/fstream.hpp>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ostream>
#include <sstream>
#include <locale>
#include <cstdlib>
#include <iostream>
#include <streambuf>
#include <stdexcept>
#include <string>
#include <chrono>
#include <cmath>
#include <limits>

using namespace Slic3r;
namespace {
DepositionMesh::Snapshot beads()
{
    return {{{{0,0,.2f},{2,0,.2f},.45f,.2f},
             {{0,1,.2f},{2,1,.8f},.45f,.2f}},0};
}
}

TEST_CASE("Sliced STEP preserves separate valid solids and millimetres", "[DepositionGeometry]")
{
    ScopedTemporaryFile file(".step");
    boost::nowide::ofstream out(file.string());
    DepositionMesh::GeometryOptions options;
    options.angular_samples = GENERATE(8u, 16u, 32u);
    const auto report = DepositionMesh::export_geometry(beads(), out, DepositionMesh::GeometryFormat::Step, options);
    out.close();
    REQUIRE(report.triangle_count > 0);
    STEPControl_Reader reader;
    REQUIRE(reader.ReadFile(file.string().c_str()) == IFSelect_RetDone);
    REQUIRE(reader.TransferRoots() == 2);
    const auto shape = reader.OneShape();
    REQUIRE_FALSE(shape.IsNull());
    std::size_t count = 0;
    for (TopExp_Explorer it(shape, TopAbs_SOLID); it.More(); it.Next()) {
        ++count;
        CHECK(BRepCheck_Analyzer(it.Current()).IsValid());
        GProp_GProps properties;
        BRepGProp::VolumeProperties(it.Current(), properties);
        CHECK(properties.Mass() > .1);
        CHECK(properties.Mass() < .3);
    }
    CHECK(count == 2);
    std::size_t face_count = 0;
    for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next()) ++face_count;
    CHECK(face_count < report.triangle_count);
    // The CAD face consolidation must not alter the triangulated envelope.
    std::ostringstream obj;
    DepositionMesh::export_geometry(beads(), obj, DepositionMesh::GeometryFormat::Obj, options);
    indexed_triangle_set reference;
    std::istringstream lines(obj.str());
    std::string line;
    while (std::getline(lines, line)) {
        if (line.rfind("v ", 0) == 0) {
            std::istringstream fields(line.substr(2));
            Vec3f p; REQUIRE(bool(fields >> p.x() >> p.y() >> p.z()));
            reference.vertices.push_back(p);
        } else if (line.rfind("f ", 0) == 0) {
            std::istringstream fields(line.substr(2));
            Vec3i32 f; REQUIRE(bool(fields >> f[0] >> f[1] >> f[2]));
            for (unsigned k = 0; k < 3; ++k) {
                REQUIRE(f[k] < 0);
                f[k] += static_cast<int>(reference.vertices.size());
                REQUIRE(f[k] >= 0);
                REQUIRE(static_cast<std::size_t>(f[k]) < reference.vertices.size());
            }
            reference.indices.push_back(f);
        }
    }
    GProp_GProps properties;
    BRepGProp::VolumeProperties(shape, properties);
    CHECK_THAT(properties.Mass(), Catch::Matchers::WithinAbs(its_volume(reference), 1e-6));
}

TEST_CASE("Sliced geometry snapshots retain processor source identities", "[DepositionGeometry]")
{
    GCodeProcessorResult result;
    GCodeProcessorResult::MoveVertex start, end;
    start.position = {0,0,.2f};
    end.position = {2,0,.2f};
    end.type = EMoveType::Extrude;
    end.width = .45f; end.height = .2f;
    end.gcode_id = 1234; end.extruder_id = 2; end.cp_color_id = 7;
    result.moves = {start,end};
    const auto input = DepositionMesh::snapshot(result);
    REQUIRE(input.segments.size() == 1);
    CHECK(input.segments.front().gcode_id == 1234);
    CHECK(input.segments.front().extruder_id == 2);
    CHECK(input.segments.front().color_id == 7);
    std::ostringstream out;
    DepositionMesh::export_geometry(input, out, DepositionMesh::GeometryFormat::Obj);
    CHECK(out.str().find("o extrusion_1_gcode_1234_extruder_2_color_7\n") != std::string::npos);
    result.moves.back().height = .08f;
    CHECK_THROWS(DepositionMesh::snapshot(result)); // 0.1 mm voxels are too coarse.
    CHECK_NOTHROW(DepositionMesh::snapshot(result, {}, {}, false));
}

TEST_CASE("Collinear path compaction removes internal caps without changing the envelope", "[DepositionGeometry]")
{
    const auto format = GENERATE(DepositionMesh::GeometryFormat::Stl, DepositionMesh::GeometryFormat::Obj,
                                DepositionMesh::GeometryFormat::Step, DepositionMesh::GeometryFormat::StepTessellated);
    auto split = beads();
    split.segments.resize(1);
    split.segments.front().end.x() = 1;
    auto continuation = split.segments.front();
    continuation.start.x() = 1;
    continuation.end.x() = 2;
    split.segments.push_back(continuation);
    auto whole = split;
    whole.segments.resize(1);
    whole.segments.front().end.x() = 2;
    std::ostringstream compact, reference, separate;
    const auto report = DepositionMesh::export_geometry(split, compact, format);
    DepositionMesh::export_geometry(whole, reference, format);
    CHECK(compact.str() == reference.str());
    DepositionMesh::GeometryOptions options;
    options.merge_collinear = false;
    const auto original = DepositionMesh::export_geometry(split, separate, format, options);
    CHECK(report.triangle_count < original.triangle_count);
    CHECK(report.output_bytes < original.output_bytes);
}

TEST_CASE("Path compaction preserves gaps bends layers and material boundaries", "[DepositionGeometry]")
{
    const unsigned change = GENERATE(0u, 1u, 2u, 3u, 4u, 5u, 6u);
    auto input = beads();
    input.segments.resize(1);
    auto next = input.segments.front();
    next.start = next.end;
    next.end.x() += 1;
    switch (change) {
    case 0: next.start.x() += .001f; break;
    case 1: next.end.y() += .001f; break;
    case 2: next.start.z() += .001f; next.end.z() += .001f; break;
    case 3: ++next.extruder_id; break;
    case 4: ++next.color_id; break;
    case 5: next.width += .001f; break;
    case 6: next.end.x() = next.start.x() - 1; break;
    }
    input.segments.push_back(next);
    std::ostringstream actual, reference;
    DepositionMesh::GeometryOptions options;
    options.merge_collinear = false;
    DepositionMesh::export_geometry(input, actual, DepositionMesh::GeometryFormat::Obj);
    DepositionMesh::export_geometry(input, reference, DepositionMesh::GeometryFormat::Obj, options);
    CHECK(actual.str() == reference.str());
}

TEST_CASE("Compact sliced STEP retains separate tessellated products", "[DepositionGeometry]")
{
    ScopedTemporaryFile file(".step");
    boost::nowide::ofstream out(file.string());
    const auto report = DepositionMesh::export_geometry(beads(), out, DepositionMesh::GeometryFormat::StepTessellated);
    out.close();
    STEPControl_Reader reader;
    REQUIRE(reader.ReadFile(file.string().c_str()) == IFSelect_RetDone);
    REQUIRE(reader.TransferRoots() == 2);
    std::size_t triangles = 0;
    std::size_t faces = 0;
    for (TopExp_Explorer it(reader.OneShape(), TopAbs_FACE); it.More(); it.Next()) {
        ++faces;
        TopLoc_Location location;
        const auto mesh = BRep_Tool::Triangulation(TopoDS::Face(it.Current()), location);
        REQUIRE_FALSE(mesh.IsNull());
        triangles += mesh->NbTriangles();
        for (int n = 1; n <= mesh->NbNodes(); ++n) {
            const auto p = mesh->Node(n).Transformed(location.Transformation());
            CHECK(p.X() > -.23);
            CHECK(p.X() < 2.23);
            CHECK((p.Y() < .23 || p.Y() > .77));
            CHECK(p.Z() > -1e-6);
            CHECK(p.Z() < .800001);
        }
    }
    CHECK(faces == 2);
    CHECK(triangles == report.triangle_count);
    CHECK(report.output_bytes < 40000);
}

TEST_CASE("Curved path compaction bounds displacement in both directions", "[DepositionGeometry]")
{
    const double tolerance = GENERATE(.001, .005);
    DepositionMesh::Snapshot input;
    const unsigned count = 512;
    const auto point = [&](unsigned i) {
        const double angle = 6.283185307179586 * i / count;
        return Vec3f(float(10 * std::cos(angle)), float(10 * std::sin(angle)), .2f);
    };
    for (unsigned i = 0; i < count; ++i)
        input.segments.push_back({point(i), point(i + 1), .45f, .2f, i, 1, 2});
    input.segments.back().end = input.segments.front().start;
    const auto compact = DepositionMesh::compact_paths(input, tolerance);
    REQUIRE(compact.segments.size() < input.segments.size());
    REQUIRE(compact.segments.size() > 3);
    const auto distance = [](const Vec3d& p, const DepositionMesh::Snapshot& target) {
        double minimum = std::numeric_limits<double>::max();
        for (const auto& s : target.segments) {
            const Vec3d a = s.start.cast<double>(), d = s.end.cast<double>() - a;
            const double t = std::clamp((p - a).dot(d) / d.squaredNorm(), 0., 1.);
            minimum = std::min(minimum, (p - a - t * d).norm());
        }
        return minimum;
    };
    for (const auto& s : input.segments)
        CHECK(distance(s.start.cast<double>(), compact) <= tolerance);
    for (const auto& s : compact.segments) {
        CHECK(s.extruder_id == 1);
        CHECK(s.color_id == 2);
        for (unsigned i = 0; i <= 16; ++i)
            CHECK(distance(s.start.cast<double>() + (s.end.cast<double>() - s.start.cast<double>()) * (double(i) / 16), input) <= tolerance);
    }
}

TEST_CASE("Curved path compaction preserves material changes and rejects invalid hidden segments", "[DepositionGeometry]")
{
    auto input = beads();
    input.segments[1] = input.segments.front();
    input.segments[1].start = input.segments.front().end;
    input.segments[1].end.x() += 1;
    input.segments[1].color_id = 1;
    const auto compact = DepositionMesh::compact_paths(input, .005);
    REQUIRE(compact.segments.size() == 2);
    CHECK(compact.segments[0].color_id == 0);
    CHECK(compact.segments[1].color_id == 1);
    input.segments[1].color_id = 0;
    input.segments[1].end = input.segments[1].start;
    CHECK_THROWS(DepositionMesh::compact_paths(input, .005));
    CHECK_THROWS(DepositionMesh::compact_paths(beads(), .005, [] { throw std::runtime_error("Canceled"); }));
}

TEST_CASE("Joined curved paths retain closed outward surfaces without internal caps", "[DepositionGeometry]")
{
    const float width = GENERATE(.2f, .45f);
    const unsigned samples = GENERATE(8u, 16u, 32u);
    DepositionMesh::Snapshot input;
    const auto point = [](unsigned i) { return Vec3f(float(10 * std::cos(i * .01)), float(10 * std::sin(i * .01)), .2f); };
    for (unsigned i = 0; i < 64; ++i) input.segments.push_back({point(i), point(i + 1), width, .2f});
    DepositionMesh::GeometryOptions options;
    options.angular_samples = samples;
    options.path_tolerance_mm = .005;
    options.join_paths = true;
    std::ostringstream output;
    const auto report = DepositionMesh::export_geometry(input, output, DepositionMesh::GeometryFormat::Obj, options);
    indexed_triangle_set mesh;
    std::istringstream lines(output.str());
    std::string line;
    while (std::getline(lines, line)) {
        std::istringstream fields(line.substr(2));
        if (line.rfind("v ", 0) == 0) {
            Vec3f p; REQUIRE(bool(fields >> p.x() >> p.y() >> p.z())); mesh.vertices.push_back(p);
        } else if (line.rfind("f ", 0) == 0) {
            Vec3i32 f; REQUIRE(bool(fields >> f[0] >> f[1] >> f[2]));
            for (unsigned k = 0; k < 3; ++k) {
                REQUIRE(f[k] < 0);
                f[k] += static_cast<int>(mesh.vertices.size());
                REQUIRE(f[k] >= 0);
                REQUIRE(static_cast<std::size_t>(f[k]) < mesh.vertices.size());
            }
            mesh.indices.push_back(f);
        }
    }
    CHECK(its_num_open_edges(mesh) == 0);
    CHECK(its_volume(mesh) > 0.);
    CHECK_FALSE(MeshBoolean::cgal::does_self_intersect(TriangleMesh(mesh)));
    for (const auto& p : mesh.vertices) {
        CHECK(p.z() >= -1e-7f);
        CHECK(p.z() <= .200001f);
        CHECK(p.head<2>().norm() >= 10 - width / 2 - .011);
        CHECK(p.head<2>().norm() <= 10 + width / 2 + .011);
    }
    options.join_paths = false;
    std::ostringstream separate;
    const auto original = DepositionMesh::export_geometry(input, separate, DepositionMesh::GeometryFormat::Obj, options);
    CHECK(report.triangle_count < original.triangle_count / 2);
}

TEST_CASE("Joined path STEP imports as valid CAD solids with material boundaries", "[DepositionGeometry]")
{
    DepositionMesh::Snapshot input;
    const auto point = [](unsigned i) { return Vec3f(float(10 * std::cos(i * .01)), float(10 * std::sin(i * .01)), .2f); };
    for (unsigned i = 0; i < 64; ++i) {
        DepositionMesh::Segment segment{point(i), point(i + 1), .45f, .2f};
        segment.extruder_id = i < 32 ? 0 : 1;
        input.segments.push_back(segment);
    }
    DepositionMesh::GeometryOptions options;
    options.path_tolerance_mm = .005;
    options.join_paths = true;
    ScopedTemporaryFile file(".step");
    boost::nowide::ofstream out(file.string());
    DepositionMesh::export_geometry(input, out, DepositionMesh::GeometryFormat::Step, options);
    out.close();
    STEPControl_Reader reader;
    REQUIRE(reader.ReadFile(file.string().c_str()) == IFSelect_RetDone);
    REQUIRE(reader.TransferRoots() == 2);
    const auto shape = reader.OneShape();
    REQUIRE_FALSE(shape.IsNull());
    std::size_t count = 0;
    for (TopExp_Explorer it(shape, TopAbs_SOLID); it.More(); it.Next()) {
        ++count;
        CHECK(BRepCheck_Analyzer(it.Current()).IsValid());
        GProp_GProps properties;
        BRepGProp::VolumeProperties(it.Current(), properties);
        CHECK(properties.Mass() > .2);
        CHECK(properties.Mass() < .4);
    }
    CHECK(count == 2);
}

TEST_CASE("Sliced geometry writes locale independent data and propagates stream failure", "[DepositionGeometry]")
{
    struct CommaDecimal : std::numpunct<char> {
        char do_decimal_point() const override { return ','; }
    };
    std::ostringstream ordinary, localized;
    localized.imbue(std::locale(std::locale::classic(), new CommaDecimal));
    const auto format = GENERATE(DepositionMesh::GeometryFormat::Obj, DepositionMesh::GeometryFormat::Step,
                                 DepositionMesh::GeometryFormat::StepTessellated, DepositionMesh::GeometryFormat::Stl);
    DepositionMesh::export_geometry(beads(), ordinary, format);
    DepositionMesh::export_geometry(beads(), localized, format);
    CHECK(ordinary.str() == localized.str());
    std::ostream failed(nullptr);
    CHECK_THROWS(DepositionMesh::export_geometry(beads(), failed, format));
}

TEST_CASE("Sliced OBJ retains separate objects without filling the intervening gap", "[DepositionGeometry]")
{
    std::ostringstream out;
    const auto report = DepositionMesh::export_geometry(beads(), out, DepositionMesh::GeometryFormat::Obj);
    CHECK(report.output_bytes == out.str().size());
    CHECK(out.str().find("o extrusion_1_gcode_0_extruder_0_color_0\n") != std::string::npos);
    CHECK(out.str().find("o extrusion_2_gcode_0_extruder_0_color_0\n") != std::string::npos);
    std::istringstream lines(out.str());
    std::string line;
    while (std::getline(lines, line)) {
        if (line.rfind("v ", 0) != 0) continue;
        std::istringstream vertex(line.substr(2));
        double x, y, z;
        REQUIRE(bool(vertex >> x >> y >> z));
        CHECK((y < .23 || y > .77));
    }
}

TEST_CASE("Sliced geometry enforces output limits and cancellation", "[DepositionGeometry]")
{
    const auto format = GENERATE(DepositionMesh::GeometryFormat::Step, DepositionMesh::GeometryFormat::StepTessellated,
                                 DepositionMesh::GeometryFormat::Obj, DepositionMesh::GeometryFormat::Stl);
    DepositionMesh::GeometryOptions options;
    options.max_output_bytes = 1024;
    std::ostringstream out;
    CHECK_THROWS(DepositionMesh::export_geometry(beads(), out, format, options));
    CHECK(out.str().size() <= options.max_output_bytes);
    CHECK_THROWS(DepositionMesh::export_geometry(beads(), out, format, {}, {},
        [] { throw std::runtime_error("Canceled"); }));
    auto invalid = beads();
    invalid.segments.front().height = -1;
    CHECK_THROWS(DepositionMesh::export_geometry(invalid, out, format));
}

TEST_CASE("Swept beads remain closed and retain their convex envelope", "[DepositionGeometry]")
{
    const unsigned samples = GENERATE(8u,16u,32u);
    for (unsigned i=0;i<24;++i) {
        const float angle = float(i)*.37f;
        DepositionMesh::Snapshot input;
        input.segments.push_back({{0,0,1},{3.f*std::cos(angle),3.f*std::sin(angle),1.f+float(int(i%3)-1)},
                                  i%2 ? .4f : .2f,.2f});
        DepositionMesh::GeometryOptions options;
        options.angular_samples=samples;
        std::ostringstream output;
        DepositionMesh::export_geometry(input,output,DepositionMesh::GeometryFormat::Obj,options);
        indexed_triangle_set mesh;
        std::istringstream lines(output.str());
        std::string line;
        while(std::getline(lines,line)) {
            std::istringstream fields(line.substr(2));
            if(line.rfind("v ",0)==0) {
                Vec3f p; REQUIRE(bool(fields>>p.x()>>p.y()>>p.z())); mesh.vertices.push_back(p);
            } else if(line.rfind("f ",0)==0) {
                Vec3i32 f; REQUIRE(bool(fields>>f[0]>>f[1]>>f[2]));
                for (unsigned k = 0; k < 3; ++k) {
                    REQUIRE(f[k] < 0);
                    f[k] += static_cast<int>(mesh.vertices.size());
                    REQUIRE(f[k] >= 0);
                    REQUIRE(static_cast<std::size_t>(f[k]) < mesh.vertices.size());
                }
                mesh.indices.push_back(f);
            }
        }
        CHECK(its_num_open_edges(mesh)==0);
        CHECK_THAT(its_volume(mesh),Catch::Matchers::WithinRel(its_volume(its_convex_hull(mesh.vertices)),1e-5f));
        std::ostringstream binary;
        const auto report=DepositionMesh::export_geometry(input,binary,DepositionMesh::GeometryFormat::Stl,options);
        CHECK(report.output_bytes==84+50*report.triangle_count);
        CHECK(report.output_bytes==binary.str().size());
        const auto bytes=binary.str();
        const auto count=std::uint32_t(static_cast<unsigned char>(bytes[80])) |
            (std::uint32_t(static_cast<unsigned char>(bytes[81]))<<8) |
            (std::uint32_t(static_cast<unsigned char>(bytes[82]))<<16) |
            (std::uint32_t(static_cast<unsigned char>(bytes[83]))<<24);
        CHECK(count==report.triangle_count);
        REQUIRE(mesh.indices.size() == count);
        // Identical indexed triangles and orientation as the independent STL
        // writer; OBJ introduces only its documented 1 nm coordinate grid.
        for (std::size_t triangle = 0; triangle < mesh.indices.size(); ++triangle)
            for (unsigned vertex = 0; vertex < 3; ++vertex)
                for (unsigned axis = 0; axis < 3; ++axis) {
                    const auto offset = 84 + triangle * 50 + 12 + vertex * 12 + axis * 4;
                    std::uint32_t bits = 0;
                    for (unsigned byte = 0; byte < 4; ++byte)
                        bits |= std::uint32_t(static_cast<unsigned char>(bytes[offset + byte])) << (8 * byte);
                    float coordinate;
                    std::memcpy(&coordinate, &bits, sizeof(coordinate));
                    CHECK_THAT(mesh.vertices[mesh.indices[triangle][vertex]][axis],
                        Catch::Matchers::WithinAbs(coordinate, 1e-6));
                }
    }
}

TEST_CASE("Sliced geometry exports a distributed real G-code sample without disk output", "[.][DepositionGeometryGCode]")
{
    const char *path = std::getenv("ORCA_DEPOSITION_GEOMETRY_GCODE");
    if (!path) SKIP("Set ORCA_DEPOSITION_GEOMETRY_GCODE for this integration check.");
    GCodeProcessor processor;
    processor.process_file(path);
    DepositionMesh::Options options;
    auto input = DepositionMesh::snapshot(processor.get_result(), options, {}, false);
    REQUIRE_FALSE(input.segments.empty());
    DepositionMesh::Snapshot sample;
    const std::size_t stride = std::getenv("ORCA_DEPOSITION_GEOMETRY_FULL") ? 1 : std::max<std::size_t>(1, input.segments.size() / 1000);
    for (std::size_t i = 0; i < input.segments.size(); i += stride)
        sample.segments.push_back(input.segments[i]);
    struct CountingBuffer : std::streambuf {
        std::uint64_t bytes = 0, position = 0;
        std::streamsize xsputn(const char *, std::streamsize n) override {
            position += n; bytes=std::max(bytes,position); return n;
        }
        pos_type seekoff(off_type offset,std::ios_base::seekdir way,std::ios_base::openmode) override {
            const auto base=way==std::ios_base::beg?0:way==std::ios_base::cur?off_type(position):off_type(bytes);
            if(base+offset<0) return pos_type(off_type(-1));
            position=std::uint64_t(base+offset); return pos_type(position);
        }
        pos_type seekpos(pos_type target,std::ios_base::openmode mode) override {
            return seekoff(off_type(target),std::ios_base::beg,mode);
        }
    } buffer;
    std::ostream out(&buffer);
    DepositionMesh::GeometryOptions geometry;
    geometry.max_output_bytes=32ULL*1024*1024*1024;
    if (const char* tolerance = std::getenv("ORCA_DEPOSITION_PATH_TOLERANCE"))
        geometry.path_tolerance_mm = std::strtod(tolerance, nullptr);
    geometry.join_paths = std::getenv("ORCA_DEPOSITION_JOIN_PATHS") != nullptr;
    const auto start=std::chrono::steady_clock::now();
    const auto report = DepositionMesh::export_geometry(sample, out, std::getenv("ORCA_DEPOSITION_GEOMETRY_STL") ?
        DepositionMesh::GeometryFormat::Stl : std::getenv("ORCA_DEPOSITION_GEOMETRY_CAD") ?
        DepositionMesh::GeometryFormat::Step : DepositionMesh::GeometryFormat::StepTessellated,geometry);
    CHECK(report.output_bytes == buffer.bytes);
    std::cout << "Geometry real G-code sample: seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
              << " total_segments=" << input.segments.size()
              << " sampled_segments=" << sample.segments.size() << " triangles=" << report.triangle_count
              << " bytes=" << report.output_bytes << '\n';
}
