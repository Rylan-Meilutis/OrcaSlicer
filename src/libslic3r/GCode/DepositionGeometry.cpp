#include "DepositionGeometry.hpp"
#include "DepositionMesh.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/libslic3r.h"
#include "libslic3r/I18N.hpp"
#include "libslic3r/LocalesUtils.hpp"
#include <Eigen/Geometry>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <limits>
#include <map>
#include <utility>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <ios>
#include <iosfwd>
#include <locale>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdio>

namespace Slic3r::DepositionMesh {
namespace {
// Same 13 significant digits as scientific/precision(12), without padded
// zeroes. STEP requires a decimal point in real literals, including integers.
class GeometryNumPut final : public std::num_put<char> {
    iter_type do_put(iter_type out, std::ios_base&, char_type, double value) const override
    {
        char text[64];
        const int length = std::snprintf(text, sizeof(text), "%.13g", value);
        if (length < 0 || std::size_t(length) >= sizeof(text))
            throw std::runtime_error(_u8L("Could not write sliced geometry."));
        const char* begin = text;
        const char* end = text + length;
        const char* exponent = std::find(begin, end, 'e');
        const bool decimal = std::find(begin, exponent, '.') != exponent;
        out = std::copy(begin, exponent, out);
        if (!decimal) *out++ = '.';
        if (exponent != end) {
            *out++ = 'E';
            out = std::copy(exponent + 1, end, out);
        }
        return out;
    }
};

// OBJ alone uses a 1 nm coordinate grid (units mm): at most 0.5 nm per
// coordinate, independent of the explicitly selected bead angular faceting.
void obj_coordinate(std::ostream& out, double value)
{
    char text[64];
    const int length = std::snprintf(text, sizeof(text), "%.6f", value);
    if (length < 0 || std::size_t(length) >= sizeof(text))
        throw std::runtime_error(_u8L("Could not write sliced geometry."));
    char* end = text + length;
    while (end > text && end[-1] == '0') --end;
    if (end > text && end[-1] == '.') --end;
    if (end - text == 2 && text[0] == '-' && text[1] == '0') {
        out << '0';
        return;
    }
    out.write(text, end - text);
}

bool continues_straight_path(const Segment& a, const Segment& b)
{
    if (a.end != b.start || a.width != b.width || a.height != b.height ||
        a.extruder_id != b.extruder_id || a.color_id != b.color_id ||
        a.start.z() != a.end.z() || b.start.z() != b.end.z() ||
        !b.end.allFinite() || b.end.cwiseAbs().maxCoeff() > 1e6f)
        return false;
    // Double differences preserve the source float coordinates. No angular
    // tolerance: tiny real bends, gaps and width changes must not disappear.
    const Vec3d u = a.end.cast<double>() - a.start.cast<double>();
    const Vec3d v = b.end.cast<double>() - b.start.cast<double>();
    return v.head<2>().squaredNorm() >= 1e-12 && u.cross(v).squaredNorm() == 0. && u.dot(v) > 0.;
}

struct BeadTemplate {
    indexed_triangle_set mesh;
    struct Edge { int a, b, first, second; };
    std::vector<Edge> edges;
};

bool same_planar_path(const Segment& a, const Segment& b)
{
    return a.end == b.start && a.width == b.width && a.height == b.height &&
        a.extruder_id == b.extruder_id && a.color_id == b.color_id &&
        a.start.z() == a.end.z() && a.end.z() == b.end.z();
}

// A mitered join differs from the rounded corner by at most
// half_width * (sec(turn/2)-1). Also keep the join far from either end
// of a short segment, where neighbouring miters could cross.
bool gentle_join(const Segment& a, const Segment& b, double tolerance)
{
    if (!same_planar_path(a, b)) return false;
    const Vec2d u = (a.end.cast<double>() - a.start.cast<double>()).head<2>();
    const Vec2d v = (b.end.cast<double>() - b.start.cast<double>()).head<2>();
    const double lu = u.norm(), lv = v.norm();
    if (lu < 1e-6 || lv < 1e-6) return false;
    const double cosine = std::clamp(u.dot(v) / (lu * lv), -1., 1.);
    const double half_cosine = std::sqrt(.5 * (1. + cosine));
    if (half_cosine <= 0.) return false;
    const double radius = .5 * a.width;
    return radius * (1. / half_cosine - 1.) <= tolerance &&
        radius * std::sqrt(.5 * (1. - cosine)) / half_cosine < .25 * std::min(lu, lv);
}

bool separated_path(const std::vector<Segment>& segments, std::size_t begin, std::size_t end, double tolerance)
{
    const double width = segments[begin].width;
    std::vector<double> length{0.};
    for (std::size_t i = begin; i < end; ++i)
        length.push_back(length.back() + (segments[i].end.cast<double>() - segments[i].start.cast<double>()).norm());
    const auto cross = [](const Vec2d& a, const Vec2d& b) { return a.x() * b.y() - a.y() * b.x(); };
    const auto distance2 = [](const Vec2d& p, const Vec2d& a, const Vec2d& b) {
        const Vec2d d = b - a;
        return (p - a - std::clamp((p - a).dot(d) / d.squaredNorm(), 0., 1.) * d).squaredNorm();
    };
    const double clearance2 = (width + 2 * tolerance) * (width + 2 * tolerance);
    for (std::size_t i = begin; i < end; ++i) for (std::size_t j = i + 2; j < end; ++j) {
        // Immediate neighbours along a smooth tube naturally overlap as
        // capsules. Only distant portions can cause a global self-contact.
        if (length[j - begin] - length[i + 1 - begin] <= 2 * width) continue;
        const Vec2d a = segments[i].start.head<2>().cast<double>(), b = segments[i].end.head<2>().cast<double>();
        const Vec2d c = segments[j].start.head<2>().cast<double>(), d = segments[j].end.head<2>().cast<double>();
        if ((a.cwiseMin(b) - c.cwiseMax(d)).maxCoeff() > width + 2 * tolerance ||
            (c.cwiseMin(d) - a.cwiseMax(b)).maxCoeff() > width + 2 * tolerance) continue;
        const double denominator = cross(b - a, d - c);
        if (denominator != 0.) {
            const double t = cross(c - a, d - c) / denominator, u = cross(c - a, b - a) / denominator;
            if (t >= 0. && t <= 1. && u >= 0. && u <= 1.) return false;
        }
        if (std::min({distance2(a,c,d), distance2(b,c,d), distance2(c,a,b), distance2(d,a,b)}) < clearance2)
            return false;
    }
    return true;
}

indexed_triangle_set joined_path(const std::vector<Segment>& segments, std::size_t begin, std::size_t end, unsigned samples)
{
    const auto& first = segments[begin];
    const double radius = .5 * first.height, flat = .5 * (first.width - first.height);
    const unsigned half = samples / 2;
    std::vector<Vec2d> profile;
    std::vector<int> positive(half + 1), negative(half + 1);
    const auto section = [&](unsigned latitude) {
        const double angle = -PI / 2 + PI * latitude / half;
        return Vec2d(flat + radius * (latitude == 0 || latitude == half ? 0. : std::cos(angle)),
                     radius * (latitude == 0 ? -1. : latitude == half ? 1. : std::sin(angle)));
    };
    for (unsigned j = 0; j <= half; ++j) { positive[j] = int(profile.size()); profile.push_back(section(j)); }
    const int top = flat == 0. ? positive[half] : int(profile.size());
    if (flat != 0.) profile.emplace_back(0., radius);
    for (unsigned j = half + 1; j-- > 0;) {
        if (flat == 0. && (j == 0 || j == half)) { negative[j] = positive[j]; continue; }
        negative[j] = int(profile.size());
        const Vec2d p = section(j); profile.emplace_back(-p.x(), p.y());
    }
    const int bottom = flat == 0. ? positive[0] : int(profile.size());
    if (flat != 0.) profile.emplace_back(0., -radius);
    const std::size_t rings = end - begin + 1, count = profile.size();
    indexed_triangle_set mesh;
    mesh.vertices.reserve(rings * count + 2 * (half + 1) * (half + 1));
    const Vec3d origin = first.start.cast<double>();
    const auto direction = [&](std::size_t i) -> Vec2d {
        return (segments[i].end.cast<double>() - segments[i].start.cast<double>()).head<2>().normalized();
    };
    for (std::size_t i = 0; i < rings; ++i) {
        const Vec3d center = (i == 0 ? first.start : segments[begin + i - 1].end).cast<double>() - origin;
        Vec2d tangent = direction(i == 0 ? begin : begin + i - 1);
        double scale = 1.;
        if (i > 0 && i + 1 < rings) {
            const Vec2d next = direction(begin + i);
            const Vec2d average = (tangent + next).normalized();
            scale = 1. / average.dot(tangent); tangent = average;
        }
        const Vec2d normal(-tangent.y() * scale, tangent.x() * scale);
        for (const auto& p : profile)
            mesh.vertices.emplace_back(float(center.x() + normal.x() * p.x()), float(center.y() + normal.y() * p.x()),
                                       float(center.z() - radius + p.y()));
    }
    const auto face = [&](int a, int b, int c, bool reverse = false) {
        if (a != b && b != c && c != a) mesh.indices.emplace_back(a, reverse ? c : b, reverse ? b : c);
    };
    for (std::size_t i = 0; i + 1 < rings; ++i) for (std::size_t j = 0; j < count; ++j) {
        const int a = int(i * count + j), b = int(i * count + (j + 1) % count);
        face(a, b, b + int(count)); face(a, b + int(count), a + int(count));
    }
    const auto cap = [&](bool start) {
        const int base = start ? 0 : int((rings - 1) * count);
        const Vec3d center = (start ? first.start : segments[end - 1].end).cast<double>() - origin;
        const Vec2d tangent = direction(start ? begin : end - 1), normal(-tangent.y(), tangent.x());
        std::vector<int> grid((half + 1) * (half + 1));
        for (unsigned j = 0; j <= half; ++j) for (unsigned k = 0; k <= half; ++k) {
            int index;
            if (k == 0) index = base + negative[j];
            else if (k == half) index = base + positive[j];
            else if (flat == 0. && (j == 0 || j == half)) index = base + positive[j];
            else {
                const Vec2d p = section(j);
                const double angle = -PI / 2 + PI * k / half;
                const double x = p.x() * std::cos(angle) * (start ? -1. : 1.), y = p.x() * std::sin(angle);
                index = int(mesh.vertices.size());
                mesh.vertices.emplace_back(float(center.x() + tangent.x() * x + normal.x() * y),
                                           float(center.y() + tangent.y() * x + normal.y() * y), float(center.z() - radius + p.y()));
            }
            grid[j * (half + 1) + k] = index;
        }
        for (unsigned j = 0; j < half; ++j) for (unsigned k = 0; k < half; ++k) {
            const int a = grid[j * (half + 1) + k], b = grid[j * (half + 1) + k + 1];
            const int c = grid[(j + 1) * (half + 1) + k + 1], d = grid[(j + 1) * (half + 1) + k];
            face(a,b,c,start); face(a,c,d,start);
        }
        if (flat != 0.) for (unsigned k = 0; k < half; ++k) {
            face(base + top, grid[half * (half + 1) + k], grid[half * (half + 1) + k + 1], start);
            face(base + bottom, grid[k + 1], grid[k], start);
        }
    };
    cap(true); cap(false);
    return mesh;
}

BeadTemplate bead_template(float width, float height, unsigned samples)
{
    const double radius = .5 * height, flat = .5 * (width - height);
    std::vector<Vec3f> points;
    points.reserve(samples * (samples / 2 + 1));
    for (unsigned latitude = 0; latitude <= samples / 2; ++latitude) {
        const double angle = PI * latitude / (samples / 2);
        const double xy = flat + radius * std::sin(angle);
        const double z = radius * (std::cos(angle) - 1.);
        for (unsigned longitude = 0; longitude < samples; ++longitude) {
            const double azimuth = 2. * PI * longitude / samples;
            points.emplace_back(float(xy * std::cos(azimuth)), float(xy * std::sin(azimuth)), float(z));
        }
    }
    BeadTemplate result;
    result.mesh = its_convex_hull(points);
    if (result.mesh.indices.empty() || its_num_open_edges(result.mesh) != 0)
        throw std::runtime_error(_u8L("Could not construct a closed extrusion bead."));
    std::map<std::pair<int,int>, std::size_t> edges;
    for (std::size_t i = 0; i < result.mesh.indices.size(); ++i) {
        const auto& face = result.mesh.indices[i];
        for (unsigned j = 0; j < 3; ++j) {
            const int a = face[j], b = face[(j+1)%3];
            const auto key = std::make_pair(std::min(a,b),std::max(a,b));
            const auto [it, inserted] = edges.emplace(key,result.edges.size());
            if (inserted) result.edges.push_back({a,b,int(i),-1});
            else result.edges[it->second].second = int(i);
        }
    }
    return result;
}

indexed_triangle_set bead(const Segment &s, const BeadTemplate& shape)
{
    if (!s.start.allFinite() || !s.end.allFinite() || !std::isfinite(s.width) ||
        !std::isfinite(s.height) || s.height <= 0 || s.width < s.height ||
        !(s.end - s.start).allFinite() || s.start.cwiseAbs().maxCoeff() > 1e6f ||
        s.end.cwiseAbs().maxCoeff() > 1e6f || s.width > 1000.f ||
        (s.end - s.start).head<2>().squaredNorm() < 1e-12f)
        throw std::invalid_argument(_u8L("Invalid spatial extrusion bead."));
    // Sweep the same sampled disc+sphere polytope along the line. Only the
    // silhouette needs new faces: a per-segment 3D convex hull is unnecessary.
    // Zero-facing facets stay at the start; their forward silhouette extends
    // them in their own plane, without duplicate overlapping caps.
    const Vec3f delta = s.end - s.start;
    const auto& base = shape.mesh;
    const int count = int(base.vertices.size());
    indexed_triangle_set mesh;
    mesh.vertices = base.vertices;
    mesh.vertices.reserve(2*count);
    for (const auto& p : base.vertices) mesh.vertices.push_back(p+delta);
    std::vector<bool> forward;
    forward.reserve(base.indices.size());
    mesh.indices.reserve(base.indices.size()+2*shape.edges.size());
    for (const auto& face : base.indices) {
        const Vec3d a = base.vertices[face[0]].cast<double>();
        const Vec3d normal = (base.vertices[face[1]].cast<double>()-a).cross(base.vertices[face[2]].cast<double>()-a);
        const bool end = normal.dot(delta.cast<double>()) > 0.;
        forward.push_back(end);
        mesh.indices.emplace_back(face[0]+(end?count:0),face[1]+(end?count:0),face[2]+(end?count:0));
    }
    for (const auto& edge : shape.edges) {
        if (edge.second < 0) throw std::runtime_error(_u8L("Could not construct a closed extrusion bead."));
        if (forward[edge.first] == forward[edge.second]) continue;
        const int a = forward[edge.first] ? edge.b : edge.a;
        const int b = forward[edge.first] ? edge.a : edge.b;
        mesh.indices.emplace_back(b,a,a+count);
        mesh.indices.emplace_back(b,a+count,b+count);
    }
    its_compactify_vertices(mesh, false);
    return mesh;
}

void write_u32(std::ostream& out, std::uint32_t value)
{
    const char bytes[]{char(value),char(value>>8),char(value>>16),char(value>>24)};
    out.write(bytes,4);
}
void stl_bead(std::ostream& out, const indexed_triangle_set& mesh, const Vec3f& origin)
{
    std::string records(50 * mesh.indices.size(), '\0');
    std::size_t offset = 0;
    const auto vector = [&](const Vec3f& p) {
        for (unsigned axis=0;axis<3;++axis) {
            std::uint32_t bits;
            const float value=p[axis];
            std::memcpy(&bits,&value,sizeof(bits));
            for (unsigned byte = 0; byte < 4; ++byte)
                records[offset++] = char(bits >> (8 * byte));
        }
    };
    for (const auto& face : mesh.indices) {
        const Vec3f a=mesh.vertices[face[0]]+origin,b=mesh.vertices[face[1]]+origin,c=mesh.vertices[face[2]]+origin;
        const Vec3f cross=(b-a).cross(c-a);
        if (!cross.allFinite() || cross.squaredNorm()==0.f)
            throw std::runtime_error(_u8L("Degenerate extrusion bead face."));
        vector(cross.normalized()); vector(a); vector(b); vector(c);
        offset += 2; // Zero attribute-byte count.
    }
    out.write(records.data(), static_cast<std::streamsize>(records.size()));
}

void step_header(std::ostream &out, bool tessellated)
{
    out << "ISO-10303-21;\nHEADER;\n"
           "FILE_DESCRIPTION(('Nominal faceted G-code beads; separate overlapping solids'),'2;1');\n"
           "FILE_NAME('sliced_geometry','',(''),(''),'OrcaSlicer','OrcaSlicer','');\n"
        << (tessellated ? "FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF'));\n" :
                          "FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\n")
        << "ENDSEC;\nDATA;\n"
        << (tessellated ?
           "#1=APPLICATION_CONTEXT('managed model based 3d engineering');\n"
           "#2=APPLICATION_PROTOCOL_DEFINITION('international standard','ap242_managed_model_based_3d_engineering',2014,#1);\n" :
           "#1=APPLICATION_CONTEXT('automotive design');\n"
           "#2=APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#1);\n")
        << "#3=PRODUCT_CONTEXT('',#1,'mechanical');\n"
           "#4=PRODUCT_DEFINITION_CONTEXT('part definition',#1,'design');\n"
           "#5=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));\n"
           "#6=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));\n"
           "#7=(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT());\n"
           "#8=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#5,#6,#7))"
           "REPRESENTATION_CONTEXT('',''));\n";
}

void step_tessellated_bead(std::ostream &out, const indexed_triangle_set &mesh, const Vec3f &origin,
                          const std::string &name, std::uint64_t &next)
{
    const auto base = next;
    out << '#' << next++ << "=COORDINATES_LIST(''," << mesh.vertices.size() << ",(";
    for (std::size_t i = 0; i < mesh.vertices.size(); ++i) {
        const Vec3d p = mesh.vertices[i].cast<double>() + origin.cast<double>();
        out << (i ? ",(" : "(") << p.x() << ',' << p.y() << ',' << p.z() << ')';
    }
    out << "));\n#" << next++ << "=TRIANGULATED_FACE('',#" << base << ',' << mesh.vertices.size()
        << ",(),$,(),(";
    for (std::size_t i = 0; i < mesh.indices.size(); ++i) {
        const auto &t = mesh.indices[i];
        out << (i ? ",(" : "(") << t[0] + 1 << ',' << t[1] + 1 << ',' << t[2] + 1 << ')';
    }
    out << "));\n#" << next++ << "=TESSELLATED_SOLID('" << name << "',(#" << base + 1 << "),$);\n";
    out << '#' << next++ << "=TESSELLATED_SHAPE_REPRESENTATION('',(#" << base + 2 << "),#8);\n";
    out << '#' << next++ << "=PRODUCT('" << name << "','" << name << "','',(#3));\n";
    out << '#' << next++ << "=PRODUCT_DEFINITION_FORMATION('','',#" << base + 4 << ");\n";
    out << '#' << next++ << "=PRODUCT_DEFINITION('design','',#" << base + 5 << ",#4);\n";
    out << '#' << next++ << "=PRODUCT_DEFINITION_SHAPE('','',#" << base + 6 << ");\n";
    out << '#' << next++ << "=SHAPE_DEFINITION_REPRESENTATION(#" << base + 7 << ",#" << base + 3 << ");\n";
}

// CAD STEP needs a plane and placement for every face. Keep connected coplanar
// triangles as one polygon instead of repeating that topology per triangle.
// Only exact coplanarity is accepted; ambiguous boundaries retain triangles.
std::vector<std::vector<int>> planar_faces(const indexed_triangle_set& mesh)
{
    std::vector<std::size_t> parent(mesh.indices.size());
    for (std::size_t i = 0; i < parent.size(); ++i) parent[i] = i;
    const auto root = [&](std::size_t i) {
        while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; }
        return i;
    };
    const auto normal = [&](std::size_t i) -> Vec3d {
        const auto& t = mesh.indices[i];
        return (mesh.vertices[t[1]].cast<double>() - mesh.vertices[t[0]].cast<double>()).cross(
            mesh.vertices[t[2]].cast<double>() - mesh.vertices[t[0]].cast<double>());
    };
    std::map<std::pair<int,int>, std::size_t> edges;
    for (std::size_t i = 0; i < mesh.indices.size(); ++i) {
        const auto& t = mesh.indices[i];
        for (unsigned j = 0; j < 3; ++j) {
            const auto key = std::minmax(t[j], t[(j + 1) % 3]);
            const auto [it, inserted] = edges.emplace(key, i);
            if (inserted) continue;
            const auto other = it->second;
            const Vec3d n = normal(i);
            bool coplanar = n.dot(normal(other)) > 0.;
            for (int v : mesh.indices[other])
                coplanar = coplanar && n.dot(mesh.vertices[v].cast<double>() - mesh.vertices[t[0]].cast<double>()) == 0.;
            if (coplanar) parent[root(other)] = root(i);
        }
    }
    std::map<std::size_t, std::vector<std::size_t>> groups;
    for (std::size_t i = 0; i < parent.size(); ++i) groups[root(i)].push_back(i);
    std::vector<std::vector<int>> result;
    for (const auto& entry : groups) {
        const auto& group = entry.second;
        const auto& first = mesh.indices[group.front()];
        const Vec3d n = normal(group.front()), origin = mesh.vertices[first[0]].cast<double>();
        std::map<std::pair<int,int>, bool> boundary;
        bool valid = true;
        for (const auto i : group) {
            const auto& t = mesh.indices[i];
            for (unsigned j = 0; j < 3; ++j) {
                valid = valid && n.dot(mesh.vertices[t[j]].cast<double>() - origin) == 0.;
                const auto reversed = boundary.find({t[(j + 1) % 3], t[j]});
                if (reversed != boundary.end()) boundary.erase(reversed);
                else if (!boundary.emplace(std::make_pair(t[j], t[(j + 1) % 3]), true).second) valid = false;
            }
        }
        std::map<int,int> next;
        for (const auto& edge : boundary)
            if (!next.emplace(edge.first.first, edge.first.second).second) valid = false;
        std::vector<int> polygon;
        if (valid && !next.empty()) {
            const int start = next.begin()->first;
            int vertex = start;
            do {
                polygon.push_back(vertex);
                const auto it = next.find(vertex);
                if (it == next.end()) { valid = false; break; }
                vertex = it->second;
                next.erase(it);
            } while (vertex != start);
            valid = valid && next.empty() && polygon.size() >= 3;
        } else valid = false;
        if (valid) result.push_back(std::move(polygon));
        else for (const auto i : group) {
            const auto& t = mesh.indices[i];
            result.push_back({t[0], t[1], t[2]});
        }
    }
    return result;
}

void step_bead(std::ostream &out, const indexed_triangle_set &mesh, const Vec3f &origin,
               const std::string &name, std::uint64_t &next)
{
    const auto vertices = next;
    for (const auto &local : mesh.vertices) {
        const Vec3d p = local.cast<double>() + origin.cast<double>();
        out << '#' << next++ << "=CARTESIAN_POINT('',(" << p.x() << ',' << p.y() << ',' << p.z() << "));\n";
    }
    std::vector<std::uint64_t> faces;
    faces.reserve(mesh.indices.size());
    for (const auto &tri : planar_faces(mesh)) {
        const auto base = next;
        const Vec3d a = mesh.vertices[tri[0]].cast<double>();
        const Vec3d b = mesh.vertices[tri[1]].cast<double>();
        Vec3d normal = Vec3d::Zero();
        // Boundary edges can be collinear after adjoining planar triangles.
        for (std::size_t i = 2; i < tri.size() && normal.squaredNorm() == 0.; ++i)
            normal = (b - a).cross(mesh.vertices[tri[i]].cast<double>() - a);
        normal.normalize();
        const Vec3d tangent = (b - a).normalized();
        if (!normal.allFinite() || !tangent.allFinite() || normal.squaredNorm() == 0. || tangent.squaredNorm() == 0.)
            throw std::runtime_error(_u8L("Degenerate extrusion bead face."));
        out << '#' << next++ << "=POLY_LOOP('',(";
        for (std::size_t i = 0; i < tri.size(); ++i) out << (i ? ",#" : "#") << vertices + tri[i];
        out << "));\n";
        out << '#' << next++ << "=FACE_OUTER_BOUND('',#" << base << ",.T.);\n";
        out << '#' << next++ << "=DIRECTION('',(" << normal.x() << ',' << normal.y() << ',' << normal.z() << "));\n";
        out << '#' << next++ << "=DIRECTION('',(" << tangent.x() << ',' << tangent.y() << ',' << tangent.z() << "));\n";
        out << '#' << next++ << "=AXIS2_PLACEMENT_3D('',#" << vertices + tri[0] << ",#" << base + 2 << ",#" << base + 3 << ");\n";
        out << '#' << next++ << "=PLANE('',#" << base + 4 << ");\n";
        faces.push_back(next);
        out << '#' << next++ << "=FACE_SURFACE('',(#" << base + 1 << "),#" << base + 5 << ",.T.);\n";
    }
    const auto shell = next;
    out << '#' << next++ << "=CLOSED_SHELL('',(";
    for (std::size_t i = 0; i < faces.size(); ++i) out << (i ? ",#" : "#") << faces[i];
    out << "));\n";
    out << '#' << next++ << "=FACETED_BREP('" << name << "',#" << shell << ");\n";
    out << '#' << next++ << "=FACETED_BREP_SHAPE_REPRESENTATION('',(#" << shell + 1 << "),#8);\n";
    out << '#' << next++ << "=PRODUCT('" << name << "','" << name << "','',(#3));\n";
    out << '#' << next++ << "=PRODUCT_DEFINITION_FORMATION('','',#" << shell + 3 << ");\n";
    out << '#' << next++ << "=PRODUCT_DEFINITION('design','',#" << shell + 4 << ",#4);\n";
    out << '#' << next++ << "=PRODUCT_DEFINITION_SHAPE('','',#" << shell + 5 << ");\n";
    out << '#' << next++ << "=SHAPE_DEFINITION_REPRESENTATION(#" << shell + 6 << ",#" << shell + 2 << ");\n";
}
} // namespace

Snapshot compact_paths(const Snapshot& input, double tolerance, const std::function<void()>& check)
{
    if (!std::isfinite(tolerance) || tolerance <= 0. || tolerance > .05 || input.segments.size() > 2000000)
        throw std::invalid_argument(_u8L("Invalid extrusion path tolerance."));
    Snapshot result;
    result.unrepresented_extruder_only_moves = input.unrepresented_extruder_only_moves;
    result.unrepresented_motion_extrusion_moves = input.unrepresented_motion_extrusion_moves;
    result.unrepresented_stationary_extrusion_moves = input.unrepresented_stationary_extrusion_moves;
    result.segments.reserve(std::min<std::size_t>(input.segments.size(), 65536));
    // Bound both temporary storage and the quadratic worst case of RDP. Chunk
    // endpoints are retained, so a boundary cannot add approximation error.
    constexpr std::size_t max_run = 4096;
    for (std::size_t begin = 0; begin < input.segments.size();) {
        if (check) check();
        const auto& first = input.segments[begin];
        std::size_t end = begin + 1;
        while (end < input.segments.size() && end - begin < max_run) {
            const auto& previous = input.segments[end - 1];
            const auto& s = input.segments[end];
            if (previous.end != s.start || first.width != s.width || first.height != s.height ||
                first.extruder_id != s.extruder_id || first.color_id != s.color_id ||
                first.start.z() != first.end.z() || s.start.z() != first.start.z() || s.end.z() != first.start.z())
                break;
            ++end;
        }
        for (std::size_t i = begin; i < end; ++i) {
            const auto& s = input.segments[i];
            if (!s.start.allFinite() || !s.end.allFinite() || !std::isfinite(s.width) || !std::isfinite(s.height) ||
                s.height <= 0.f || s.width < s.height || s.width > 1000.f ||
                s.start.cwiseAbs().maxCoeff() > 1e6f || s.end.cwiseAbs().maxCoeff() > 1e6f ||
                (s.end - s.start).head<2>().squaredNorm() < 1e-12f)
                throw std::invalid_argument(_u8L("Invalid spatial extrusion bead."));
        }
        const auto point = [&](std::size_t i) -> Vec3d {
            return (i == begin ? first.start : input.segments[i - 1].end).cast<double>();
        };
        std::vector<std::pair<std::size_t, std::size_t>> pending{{begin, end}};
        while (!pending.empty()) {
            if (check) check();
            const auto [a, b] = pending.back();
            pending.pop_back();
            const Vec3d start = point(a), delta = point(b) - start;
            const double length2 = delta.squaredNorm();
            std::size_t split = a;
            double worst = -1.;
            for (std::size_t i = a + 1; i < b; ++i) {
                if (check && (i - a) % 256 == 0) check();
                const Vec3d relative = point(i) - start;
                const double t = length2 > 0. ? std::clamp(relative.dot(delta) / length2, 0., 1.) : 0.;
                const double distance2 = (relative - t * delta).squaredNorm();
                if (distance2 > worst) { worst = distance2; split = i; }
            }
            // Leave a rounding margin when deciding whether to replace a path.
            // Keeping an endpoint never introduces a new coordinate.
            const double allowed = std::max(0., tolerance - 1e-9);
            if (b > a + 1 && (worst > allowed * allowed || delta.head<2>().squaredNorm() < 1e-12)) {
                pending.emplace_back(split, b);
                pending.emplace_back(a, split);
            } else {
                Segment s = input.segments[a];
                s.end = input.segments[b - 1].end;
                result.segments.push_back(s);
            }
        }
        begin = end;
    }
    return result;
}

Report export_geometry(const Snapshot &source, std::ostream &out, GeometryFormat format,
                       const GeometryOptions &options, const std::function<void(int)> &progress,
                       const std::function<void()> &checkpoint)
{
    const CNumericLocalesSetter numeric_locale;
    if (!std::isfinite(options.path_tolerance_mm) || options.path_tolerance_mm < 0. || options.path_tolerance_mm > .05 ||
        (options.join_paths && options.path_tolerance_mm == 0.))
        throw std::invalid_argument(_u8L("Invalid extrusion path tolerance."));
    Snapshot compact;
    if (options.path_tolerance_mm > 0.) compact = compact_paths(source, options.path_tolerance_mm, checkpoint);
    const Snapshot& input = options.path_tolerance_mm > 0. ? compact : source;
    if (input.segments.empty() || input.segments.size() > 2000000 ||
        options.angular_samples < 8 || options.angular_samples > 32 || options.angular_samples % 4 != 0 ||
        options.max_output_bytes < 1024 || options.max_output_bytes > 64ULL * 1024 * 1024 * 1024 ||
        (format != GeometryFormat::Step && format != GeometryFormat::StepTessellated && format != GeometryFormat::Obj &&
         format != GeometryFormat::Stl))
        throw std::invalid_argument(_u8L("Invalid sliced geometry export options."));
    Report report;
    const auto start = format == GeometryFormat::Stl ? out.tellp() : std::streampos(0);
    if (format == GeometryFormat::Stl && start == std::streampos(-1))
        throw std::invalid_argument(_u8L("Deposition STL requires a seekable private output stream."));
    // Bounded per-segment staging also lets us enforce bytes before writing, on
    // non-seekable streams and even when the caller uses a non-classic locale.
    auto write = [&](const std::string &data) {
        if (data.size() > options.max_output_bytes - report.output_bytes)
            throw std::runtime_error(_u8L("Sliced geometry exceeds the output size limit."));
        out.write(data.data(), static_cast<std::streamsize>(data.size()));
        if (!out) throw std::runtime_error(_u8L("Could not write sliced geometry."));
        report.output_bytes += data.size();
    };
    std::ostringstream chunk;
    chunk.imbue(std::locale(std::locale::classic(), new GeometryNumPut));
    if (checkpoint) checkpoint();
    if (format == GeometryFormat::Stl) {
        std::string header="OrcaSlicer nominal separate beads; millimetres; overlaps are not fused";
        header.resize(80,'\0');
        chunk.write(header.data(),header.size()); write_u32(chunk,0);
    }
    else if (format != GeometryFormat::Obj) step_header(chunk, format == GeometryFormat::StepTessellated);
    else chunk << "# OrcaSlicer nominal separate extrusion beads; units millimetres; overlaps are not fused\n";
    write(chunk.str());
    std::uint64_t next = 9;
    // Bounded exact-dimension cache. No width/height quantization or reduction
    // of angular detail; unusual varying-flow paths simply rebuild templates.
    std::map<std::pair<float,float>,BeadTemplate> shapes;
    for (std::size_t i = 0; i < input.segments.size(); ++i) {
        if (checkpoint) checkpoint();
        const std::size_t first_segment = i;
        auto segment = input.segments[i];
        if (options.merge_collinear) {
            while (i + 1 < input.segments.size() && continues_straight_path(segment, input.segments[i + 1])) {
                if (checkpoint) checkpoint();
                segment.end = input.segments[++i].end;
            }
        }
        if (!std::isfinite(segment.width) || !std::isfinite(segment.height) || segment.height <= 0.f ||
            segment.width < segment.height || segment.width > 1000.f)
            throw std::invalid_argument(_u8L("Invalid spatial extrusion bead."));
        const auto key = std::make_pair(segment.width,segment.height);
        auto found = shapes.find(key);
        if (found == shapes.end()) {
            if (shapes.size() == 64) shapes.clear();
            found = shapes.emplace(key,bead_template(segment.width,segment.height,options.angular_samples)).first;
        }
        std::size_t run_end = i + 1;
        if (options.join_paths) {
            while (run_end < input.segments.size() && run_end - first_segment < 128 &&
                   gentle_join(input.segments[run_end - 1], input.segments[run_end], options.path_tolerance_mm))
                ++run_end;
        }
        const bool joined = options.join_paths && run_end > first_segment + 1 && run_end - first_segment <= 128 &&
            separated_path(input.segments, first_segment, run_end, options.path_tolerance_mm);
        const auto mesh = joined ? joined_path(input.segments, first_segment, run_end, options.angular_samples) : bead(segment, found->second);
        if (joined) i = run_end - 1;
        const std::string name = "extrusion_" + std::to_string(first_segment + 1) + "_gcode_" + std::to_string(segment.gcode_id) +
            "_extruder_" + std::to_string(segment.extruder_id) + "_color_" + std::to_string(segment.color_id);
        chunk.str({}); chunk.clear();
        if (format == GeometryFormat::Stl) {
            if (mesh.indices.size() > std::numeric_limits<std::uint32_t>::max()-report.triangle_count)
                throw std::runtime_error(_u8L("Sliced geometry exceeds the output size limit."));
            stl_bead(chunk,mesh,segment.start);
        }
        else if (format == GeometryFormat::Step) step_bead(chunk, mesh, segment.start, name, next);
        else if (format == GeometryFormat::StepTessellated) step_tessellated_bead(chunk, mesh, segment.start, name, next);
        else {
            chunk << "o " << name << '\n';
            for (const auto &local : mesh.vertices) {
                const Vec3d p = local.cast<double>() + segment.start.cast<double>();
                chunk << "v "; obj_coordinate(chunk, p.x()); chunk << ' ';
                obj_coordinate(chunk, p.y()); chunk << ' '; obj_coordinate(chunk, p.z()); chunk << '\n';
            }
            // Relative indices refer to this path's just-written vertices;
            // they preserve connectivity without growing with print length.
            const auto count = static_cast<std::int64_t>(mesh.vertices.size());
            for (const auto &face : mesh.indices)
                chunk << "f " << std::int64_t(face[0]) - count << ' ' << std::int64_t(face[1]) - count << ' '
                      << std::int64_t(face[2]) - count << '\n';
        }
        write(chunk.str());
        report.triangle_count += mesh.indices.size();
        if (progress) progress(static_cast<int>(100 * (i + 1) / input.segments.size()));
    }
    if (checkpoint) checkpoint();
    if (format == GeometryFormat::Stl) {
        const auto end=out.tellp();
        out.seekp(start+std::streamoff(80)); write_u32(out,std::uint32_t(report.triangle_count)); out.seekp(end);
        if (!out) throw std::runtime_error(_u8L("Could not write sliced geometry."));
    }
    else if (format != GeometryFormat::Obj) write("ENDSEC;\nEND-ISO-10303-21;\n");
    return report;
}
} // namespace Slic3r::DepositionMesh
