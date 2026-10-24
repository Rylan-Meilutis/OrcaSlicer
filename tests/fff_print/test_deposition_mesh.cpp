#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "libslic3r/GCode/DepositionMesh.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Utils.hpp"
#include "libslic3r/TextureBake/TextureBakeDecimate.hpp"
#include "libslic3r/TextureBake/TextureBakeIndex.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <openvdb/math/Vec3.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <map>
#include <iostream>
#include <iomanip>
#include <ios>
#include <limits>
#include <ostream>
#include <streambuf>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace Slic3r;
using Catch::Matchers::WithinAbs;

namespace {
DepositionMesh::Snapshot line(float y = 0.f)
{
    return {{{Eigen::Vector3f(0.f,y,.4f),Eigen::Vector3f(4.f,y,.4f),.4f,.4f}},0};
}
std::uint32_t read_u32(const std::string &s, std::size_t offset)
{
    std::uint32_t value = 0;
    for (unsigned i=0;i<4;++i) value |= std::uint32_t(static_cast<unsigned char>(s.at(offset+i))) << (8*i);
    return value;
}
float read_f32(const std::string &s, std::size_t offset)
{
    const auto bits = read_u32(s,offset);
    float value;
    std::memcpy(&value,&bits,sizeof(value));
    return value;
}
double volume(const std::string &s)
{
    double result = 0.;
    for(std::uint32_t i=0;i<read_u32(s,80);++i) {
        Eigen::Vector3d p[3];
        for(unsigned j=0;j<3;++j) for(unsigned k=0;k<3;++k)
            p[j][k]=read_f32(s,84+50*std::size_t(i)+12+12*j+4*k);
        result += p[0].dot(p[1].cross(p[2]))/6.;
    }
    return std::abs(result);
}
void require_closed_surface(const std::string &s)
{
    using Point=std::array<float,3>;
    std::map<std::pair<Point,Point>,unsigned> incidence;
    for(std::uint32_t i=0;i<read_u32(s,80);++i) {
        Point points[3];
        for(unsigned j=0;j<3;++j) for(unsigned k=0;k<3;++k)
            points[j][k]=read_f32(s,84+50*std::size_t(i)+12+12*j+4*k);
        for(unsigned j=0;j<3;++j) {
            auto a=points[j], b=points[(j+1)%3];
            if(b<a) std::swap(a,b);
            ++incidence[{a,b}];
        }
    }
    for(const auto &edge:incidence) REQUIRE(edge.second==2);
}
std::string mesh(const DepositionMesh::Snapshot &s,unsigned tile_size=64,bool reduce=false)
{
    std::ostringstream out(std::ios::binary);
    auto options=DepositionMesh::Options{};
    options.voxel_size_mm=.05;
    options.tile_size_cells=tile_size;
    options.reduce_horizontal_faces=reduce;
    const auto report=DepositionMesh::export_stl(s,out,options);
    REQUIRE(report.triangle_count > 0);
    REQUIRE(report.output_bytes == out.str().size());
    require_closed_surface(out.str());
    return out.str();
}
std::vector<std::array<float,9>> canonical_triangles(const std::string &s)
{
    std::vector<std::array<float,9>> result;
    for(std::uint32_t i=0;i<read_u32(s,80);++i) {
        std::array<std::array<float,3>,3> vertices;
        for(unsigned j=0;j<3;++j) for(unsigned k=0;k<3;++k)
            vertices[j][k]=read_f32(s,84+50*std::size_t(i)+12+12*j+4*k);
        const auto first=std::min_element(vertices.begin(),vertices.end())-vertices.begin();
        std::array<float,9> triangle;
        for(unsigned j=0;j<3;++j) for(unsigned k=0;k<3;++k) triangle[3*j+k]=vertices[(first+j)%3][k];
        result.push_back(triangle);
    }
    std::sort(result.begin(),result.end());
    return result;
}
}

TEST_CASE("Tiled deposition preserves whole-grid triangles across lattice faces edges and corners", "[DepositionMesh]")
{
    DepositionMesh::Snapshot s;
    s.segments.push_back({Eigen::Vector3f(-1.6f,-1.6f,.4f),Eigen::Vector3f(1.6f,1.6f,1.2f),.4f,.4f});
    s.segments.push_back({Eigen::Vector3f(-1.6f,.8f,.8f),Eigen::Vector3f(1.6f,.8f,.8f),.6f,.4f});
    s.segments.push_back({Eigen::Vector3f(-.8f,-1.6f,.8f),Eigen::Vector3f(-.8f,1.6f,.8f),.6f,.4f});
    const auto reference=canonical_triangles(mesh(s,0));
    const unsigned size=GENERATE(16u,32u,64u);
    CHECK(canonical_triangles(mesh(s,size))==reference);
}

TEST_CASE("Binary deposition records match scalar little-endian serialization", "[DepositionMesh]")
{
    const auto data=mesh(line());
    bool equivalent=true;
    for(std::uint32_t i=0;i<read_u32(data,80);++i) {
        openvdb::math::Vec3<float> p[3];
        const std::size_t begin=84+50*std::size_t(i);
        for(unsigned j=0;j<3;++j) for(unsigned k=0;k<3;++k) p[j][k]=read_f32(data,begin+12+12*j+4*k);
        auto normal=(openvdb::math::Vec3<double>(p[1])-openvdb::math::Vec3<double>(p[0])).cross(
                    openvdb::math::Vec3<double>(p[2])-openvdb::math::Vec3<double>(p[0]));
        normal/=normal.length();
        std::ostringstream reference(std::ios::binary);
        for(const auto &v:{openvdb::math::Vec3<float>(normal),p[0],p[1],p[2]})
            for(unsigned k=0;k<3;++k) {
                const float value=v[k];
                std::uint32_t bits;
                std::memcpy(&bits,&value,sizeof(bits));
                const char bytes[4]={char(bits),char(bits>>8),char(bits>>16),char(bits>>24)};
                reference.write(bytes,4);
            }
        reference.put(0); reference.put(0);
        equivalent &= reference.str()==data.substr(begin,50);
    }
    CHECK(equivalent);
}

TEST_CASE("Collinear sampled face stitches preserve a closed bead across tile boundaries", "[DepositionMesh]")
{
    const bool translated=GENERATE(false,true);
    const Eigen::Vector3f offset=translated?Eigen::Vector3f(1.5f,.85f,.6f):Eigen::Vector3f::Zero();
    DepositionMesh::Snapshot s;
    s.segments.push_back({Eigen::Vector3f(49.5019989f,104.912003f,3.f)+offset,
                          Eigen::Vector3f(50.0180511f,104.493271f,3.f)+offset,
                          .449999988f,.400000006f});
    const auto reference=canonical_triangles(mesh(s,0));
    for(unsigned size:{16u,32u,64u}) {
        const bool same_surface=canonical_triangles(mesh(s,size))==reference;
        CAPTURE(translated,size);
        CHECK(same_surface);
    }
}

TEST_CASE("Exact horizontal reduction preserves curved faces holes and closed tile seams", "[DepositionMesh]")
{
    DepositionMesh::Snapshot s;
    for(float z:{.4f,.6f}) {
        s.segments.push_back({Eigen::Vector3f(0,0,z),Eigen::Vector3f(4,0,z),.8f,.2f});
        s.segments.push_back({Eigen::Vector3f(4,0,z),Eigen::Vector3f(4,4,z),.8f,.2f});
        s.segments.push_back({Eigen::Vector3f(4,4,z),Eigen::Vector3f(0,4,z),.8f,.2f});
        s.segments.push_back({Eigen::Vector3f(0,4,z),Eigen::Vector3f(0,0,z),.8f,.2f});
    }
    const auto raw=mesh(s,32,false),reduced=mesh(s,32,true);
    CHECK(read_u32(reduced,80)<read_u32(raw,80));
    CHECK_THAT(volume(reduced),WithinAbs(volume(raw),1e-5));
    const auto curved=[](const std::string &data) {
        auto triangles=canonical_triangles(data);
        triangles.erase(std::remove_if(triangles.begin(),triangles.end(),[](const auto &t) {
            return t[2]==t[5] && t[2]==t[8];
        }),triangles.end());
        return triangles;
    };
    const bool same_curved=curved(raw)==curved(reduced);
    CHECK(same_curved);
    const auto area=[](const std::string &data) {
        double total=0.;
        for(const auto &t:canonical_triangles(data)) {
            const Eigen::Vector3d a(t[0],t[1],t[2]),b(t[3],t[4],t[5]),c(t[6],t[7],t[8]);
            total+=(b-a).cross(c-a).norm()*.5;
        }
        return total;
    };
    CHECK_THAT(area(reduced),WithinAbs(area(raw),1e-5));
    for(const auto &t:canonical_triangles(reduced)) {
        const double x=(t[0]+t[3]+t[6])/3.,y=(t[1]+t[4]+t[7])/3.;
        const bool fills_hole=x>1. && x<3. && y>1. && y<3.;
        CHECK_FALSE(fills_hole);
    }
}

TEST_CASE("Exact planar reduction pins disconnected vertex fans", "[DepositionMesh]")
{
    // Synthetic bowtie: two independently triangulated planar sheets share
    // only one vertex identity. Their fans must not be collapsed together.
    TextureBake::TriSoup input;
    const Eigen::Vector3f center(.013f,.007f,.5f);
    for(float angle:{0.f,.37f}) {
        const auto point=[&](int x,int y) {
            return Eigen::Vector3f(center.x()+x*std::cos(angle)-y*std::sin(angle),
                                   center.y()+x*std::sin(angle)+y*std::cos(angle),center.z());
        };
        for(int y=-3;y<3;++y) for(int x=-3;x<3;++x) {
            for(const auto &p:{point(x,y),point(x+1,y),point(x+1,y+1),
                              point(x,y),point(x+1,y+1),point(x,y+1)}) input.pos.push_back(p);
        }
    }
    const auto fans=[&](const TextureBake::TriSoup &s) {
        using Point=std::array<float,3>;
        std::vector<std::array<Point,2>> links;
        for(std::size_t i=0;i<s.pos.size();i+=3) for(unsigned k=0;k<3;++k)
            if(s.pos[i+k]==center) {
                const auto &a=s.pos[i+(k+1)%3],&b=s.pos[i+(k+2)%3];
                links.push_back({Point{a.x(),a.y(),a.z()},Point{b.x(),b.y(),b.z()}});
            }
        std::vector<bool> visited(links.size(),false);
        unsigned components=0;
        for(std::size_t seed=0;seed<links.size();++seed) if(!visited[seed]) {
            ++components; visited[seed]=true;
            bool changed=true;
            while(changed) {
                changed=false;
                for(std::size_t i=0;i<links.size();++i) if(visited[i])
                    for(std::size_t j=0;j<links.size();++j) if(!visited[j])
                        if(links[i][0]==links[j][0] || links[i][0]==links[j][1] ||
                           links[i][1]==links[j][0] || links[i][1]==links[j][1]) {
                            visited[j]=true; changed=true;
                        }
            }
        }
        return components;
    };
    REQUIRE(fans(input)==2);
    const auto output=TextureBake::decimate(input,input.pos.size()/3,true,0.,{},{},{},true);
    CHECK(fans(output.geometry)==2);
}

TEST_CASE("Partially overlapping beads have the same union as one continuous bead", "[DepositionMesh]")
{
    auto overlap=line();
    auto second=overlap.segments.front();
    second.start.x()=2.f; second.end.x()=6.f;
    overlap.segments.push_back(second);
    auto continuous=line();
    continuous.segments.front().end.x()=6.f;
    const bool same_union=mesh(overlap)==mesh(continuous);
    CHECK(same_union);
}

TEST_CASE("Touching flat bead layers export a closed union without internal interface skins", "[DepositionMesh]")
{
    auto s=line();
    s.segments.front().width=.8f;
    auto top=s.segments.front();
    top.start.z()+=.4f; top.end.z()+=.4f;
    s.segments.push_back(top);
    const auto result=mesh(s);
    // The shared flat contact is strictly inside the union. A triangle there
    // would be an unwanted interior skin, not an external printed surface.
    for(std::uint32_t i=0;i<read_u32(result,80);++i) {
        bool internal=true;
        for(unsigned j=0;j<3;++j) {
            const auto offset=84+50*std::size_t(i)+12+12*j;
            const float x=read_f32(result,offset),y=read_f32(result,offset+4),z=read_f32(result,offset+8);
            internal &= x>.5f && x<3.5f && std::abs(y)<.1f && std::abs(z-.4f)<1e-5f;
        }
        CHECK_FALSE(internal);
    }
}

TEST_CASE("Translated touching layers do not gain float-roundoff interface skins", "[DepositionMesh]")
{
    const double contact=GENERATE(2.,4.2,20.2,200.4,4.225,20.225,200.425);
    auto s=line();
    s.segments.front().width=.8f;
    s.segments.front().height=.2f;
    s.segments.front().start.z()=s.segments.front().end.z()=float(contact);
    auto top=s.segments.front();
    top.start.z()=top.end.z()=float(contact+.2);
    s.segments.push_back(top);
    const auto result=mesh(s);
    unsigned internal_triangles=0;
    for(std::uint32_t i=0;i<read_u32(result,80);++i) {
        bool internal=true;
        for(unsigned j=0;j<3;++j) {
            const auto offset=84+50*std::size_t(i)+12+12*j;
            const float x=read_f32(result,offset),y=read_f32(result,offset+4),z=read_f32(result,offset+8);
            internal &= x>.5f && x<3.5f && std::abs(y)<.1f && std::abs(double(z)-contact)<2e-5;
        }
        if(internal) ++internal_triangles;
    }
    CAPTURE(contact);
    CHECK(internal_triangles==0);
}

TEST_CASE("Flat-layer gaps larger than float uncertainty retain their interface surfaces", "[DepositionMesh]")
{
    const double contact=GENERATE(4.2,200.4);
    auto s=line();
    s.segments.front().width=.8f;
    s.segments.front().height=.2f;
    s.segments.front().start.z()=s.segments.front().end.z()=float(contact-.0005);
    auto top=s.segments.front();
    top.start.z()=top.end.z()=float(contact+.2005);
    s.segments.push_back(top);
    const auto result=mesh(s);
    unsigned interface_triangles=0,crossing_triangles=0;
    for(std::uint32_t i=0;i<read_u32(result,80);++i) {
        bool interface=true;
        double low=std::numeric_limits<double>::infinity(),high=-low;
        for(unsigned j=0;j<3;++j) {
            const auto offset=84+50*std::size_t(i)+12+12*j;
            const float x=read_f32(result,offset),y=read_f32(result,offset+4),z=read_f32(result,offset+8);
            interface &= x>.5f && x<3.5f && std::abs(y)<.1f && std::abs(double(z)-contact)<.002;
            low=std::min(low,double(z)); high=std::max(high,double(z));
        }
        if(interface) ++interface_triangles;
        if(low<contact && high>contact) ++crossing_triangles;
    }
    CAPTURE(contact);
    CHECK(interface_triangles>0);
    CHECK(crossing_triangles==0);
}

TEST_CASE("Deposition union contains no duplicate overlapping bead skins", "[DepositionMesh]")
{
    auto s=line();
    const auto first=mesh(s);
    s.segments.push_back(s.segments.front());
    const auto duplicate=mesh(s);
    const bool same_union=duplicate==first;
    CHECK(same_union);
    const double radius=.2, length=4.;
    const double expected=std::acos(-1.)*radius*radius*length+4./3.*std::acos(-1.)*radius*radius*radius;
    CHECK_THAT(volume(first),WithinAbs(expected,expected*.08));
}

TEST_CASE("Disconnected deposited lines retain their air gap in one STL", "[DepositionMesh]")
{
    auto s=line();
    const auto first=mesh(s);
    s.segments.push_back(line(1.f).segments.front());
    const auto both=mesh(s);
    CHECK_THAT(volume(both),WithinAbs(2.*volume(first),.001));
    for(std::uint32_t i=0;i<read_u32(both,80);++i) {
        float low=10.f, high=-10.f;
        for(unsigned j=0;j<3;++j) {
            const float y=read_f32(both,84+50*std::size_t(i)+12+12*j+4);
            low=std::min(low,y); high=std::max(high,y);
        }
        const bool spans_gap = low < .3f && high > .7f;
        CHECK_FALSE(spans_gap);
    }
}

TEST_CASE("Deposition reconstruction rejects cancellation and budgets before output", "[DepositionMesh]")
{
    auto s=line();
    std::ostringstream out;
    DepositionMesh::Options options;
    options.max_output_bytes=64ULL*1024*1024*1024+1;
    CHECK_THROWS_AS(DepositionMesh::export_stl(s,out,options),std::invalid_argument);
    CHECK(out.str().empty());
    options={};
    options.max_sample_visits=50000000001ULL;
    CHECK_THROWS_AS(DepositionMesh::export_stl(s,out,options),std::invalid_argument);
    CHECK(out.str().empty());
    options={};
    options.max_sample_visits=1;
    CHECK_THROWS_AS(DepositionMesh::export_stl(s,out,options),std::runtime_error);
    CHECK(out.str().empty());
    options={}; options.max_output_bytes=84;
    CHECK_THROWS_AS(DepositionMesh::export_stl(s,out,options),std::runtime_error);
    CHECK(out.str().empty());
    CHECK_THROWS_AS(DepositionMesh::export_stl(s,out,{}, {}, []{throw std::runtime_error("cancel");}),std::runtime_error);
    CHECK(out.str().empty());
    s=line();
    s.segments[0].end=s.segments[0].start;
    s.segments[0].end.z()+=.000005f;
    CHECK_THROWS_AS(DepositionMesh::export_stl(s,out),std::runtime_error);
    s=line(); s.segments[0].width=std::numeric_limits<float>::max();
    CHECK_THROWS_AS(DepositionMesh::export_stl(s,out),std::runtime_error);
}

TEST_CASE("Sloped deposited beads retain their three dimensional path and closed surface", "[DepositionMesh]")
{
    auto s=line();
    s.segments.front().end.z()+=1.f;
    const auto result=mesh(s);
    const double radius=.2, length=std::sqrt(17.);
    const double expected=std::acos(-1.)*radius*radius*length+4./3.*std::acos(-1.)*radius*radius*radius;
    CHECK_THAT(volume(result),WithinAbs(expected,expected*.08));
    float highest=-1.f, lowest=10.f;
    for(std::uint32_t i=0;i<read_u32(result,80);++i)
        for(unsigned j=0;j<3;++j) {
            const float z=read_f32(result,84+50*std::size_t(i)+12+12*j+8);
            highest=std::max(highest,z); lowest=std::min(lowest,z);
        }
    CHECK_THAT(highest,WithinAbs(1.4,.05));
    CHECK_THAT(lowest,WithinAbs(0.,.05));
}

TEST_CASE("Deposition snapshots include tower and support independent of preview visibility", "[DepositionMesh]")
{
    GCodeProcessorResult result;
    result.moves.resize(3);
    result.moves[0].position=Eigen::Vector3f(0,0,.2f);
    for(std::size_t i=1;i<3;++i) {
        auto &m=result.moves[i];
        m.position=Eigen::Vector3f(float(i),0,.2f);
        m.type=EMoveType::Extrude; m.width=.4f; m.height=.2f;
        m.extrusion_role=i==1?erWipeTower:erSupportMaterial;
    }
    const auto s=DepositionMesh::snapshot(result);
    REQUIRE(s.segments.size()==2);
    CHECK_THAT(s.segments.back().end.x(),WithinAbs(2.,1e-6));
}

TEST_CASE("Parsed arc extrusions are reconstructed as interpolated segments", "[DepositionMesh]")
{
    const bool was_bbl=GCodeProcessor::s_IsBBLPrinter;
    const ScopeGuard restore([was_bbl]{GCodeProcessor::s_IsBBLPrinter=was_bbl;});
    GCodeProcessor::s_IsBBLPrinter=false;
    GCodeProcessor processor;
    processor.apply_config(FullPrintConfig{});
    processor.process_buffer("G90\nM83\n;TYPE:Outer wall\n;WIDTH:0.4\n;HEIGHT:0.2\nG1 X10 Y10 Z0.2\nG2 X12 Y10 I1 J0 E0.1 F600\n");
    const auto s=DepositionMesh::snapshot(processor.get_result());
    REQUIRE(s.segments.size()>2);
    bool curved=false;
    for(const auto &segment:s.segments) curved |= std::abs(segment.end.y()-10.f)>.1f;
    CHECK(curved);
}

TEST_CASE("Unrepresented positive extrusion is inventoried rather than silently treated as travel", "[DepositionMesh]")
{
    GCodeProcessorResult result;
    result.moves.resize(5);
    result.moves[0].position=Eigen::Vector3f(0,0,.2f);
    auto &bead=result.moves[1];
    bead.position=Eigen::Vector3f(1,0,.2f); bead.type=EMoveType::Extrude;
    bead.width=.4f; bead.height=.2f; bead.delta_extruder=.1f;
    result.moves[2]=bead; result.moves[2].type=EMoveType::Travel;
    result.moves[3]=bead; result.moves[3].type=EMoveType::Wipe;
    result.moves[4]=bead; result.moves[4].type=EMoveType::Unretract;
    const auto s=DepositionMesh::snapshot(result);
    CHECK(s.segments.size()==1);
    CHECK(s.unrepresented_motion_extrusion_moves==2);
    CHECK(s.unrepresented_extruder_only_moves==1);
}

TEST_CASE("Rounded arc endpoint extrusion is inventoried without dropping valid neighboring beads", "[DepositionMesh]")
{
    GCodeProcessorResult result;
    result.moves.resize(4);
    result.moves[0].position=Eigen::Vector3f(0,0,.2f);
    auto &bead=result.moves[1];
    bead.position=Eigen::Vector3f(1,0,.2f); bead.type=EMoveType::Extrude;
    bead.width=.4f; bead.height=.2f; bead.delta_extruder=.1f;
    result.moves[2]=bead; result.moves[2].delta_extruder=.00001f;
    result.moves[3]=bead; result.moves[3].position.x()=2.f;
    const auto s=DepositionMesh::snapshot(result);
    REQUIRE(s.segments.size()==2);
    CHECK(s.unrepresented_stationary_extrusion_moves==1);
    CHECK_THAT(s.segments[1].start.x(),WithinAbs(1.,1e-6));
    CHECK_THAT(s.segments[1].end.x(),WithinAbs(2.,1e-6));
}

TEST_CASE("A saved G-code project can be benchmarked without writing an STL", "[.][DepositionMeshProject]")
{
    const char *path=std::getenv("ORCA_DEPOSITION_GCODE");
    if(!path || !*path) SKIP("Set ORCA_DEPOSITION_GCODE to run this read-only project benchmark.");
    class CountingBuffer : public std::streambuf {
    public:
        std::uint64_t bytes=0;
        std::uint64_t position=0;
        std::streamsize xsputn(const char *, std::streamsize count) override {
            position+=std::uint64_t(count); bytes=std::max(bytes,position); return count;
        }
        int_type overflow(int_type ch) override {
            if(!traits_type::eq_int_type(ch,traits_type::eof())) { ++position; bytes=std::max(bytes,position); }
            return traits_type::not_eof(ch);
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
    std::ostream sink(&buffer);
    DepositionMesh::Options options;
    if(const char *resolution=std::getenv("ORCA_DEPOSITION_VOXEL")) options.voxel_size_mm=std::stod(resolution);
    if(const char *work=std::getenv("ORCA_DEPOSITION_MAX_VISITS")) options.max_sample_visits=std::stoull(work);
    if(const char *bytes=std::getenv("ORCA_DEPOSITION_MAX_BYTES")) options.max_output_bytes=std::stoull(bytes);
    if(const char *bytes=std::getenv("ORCA_DEPOSITION_MEMORY")) options.max_memory_bytes=std::stoull(bytes);
    const auto start=std::chrono::steady_clock::now();
    GCodeProcessor processor;
    processor.process_file(path);
    const auto &moves=processor.get_result().moves;
    std::size_t slopes=0;
    float minimum_width=std::numeric_limits<float>::infinity(),minimum_height=minimum_width;
    std::size_t degenerate=0,undersized=0;
    for(std::size_t i=1;i<moves.size();++i) if(moves[i].type==EMoveType::Extrude) {
        const auto &m=moves[i];
        minimum_width=std::min(minimum_width,m.width); minimum_height=std::min(minimum_height,m.height);
        const bool stationary=(m.position-moves[i-1].position).head<2>().squaredNorm()<1e-12f;
        const bool small=options.voxel_size_mm>.5*std::min(double(m.width),double(m.height))+1e-7;
        if((stationary && degenerate++<3) || (small && undersized++<3))
            std::cout << "Deposition unresolved move: index=" << i << " gcode_id=" << m.gcode_id
                      << " previous_id=" << moves[i-1].gcode_id << " role=" << int(m.extrusion_role)
                      << " start=" << moves[i-1].position.transpose() << " end=" << m.position.transpose()
                      << " width=" << m.width << " height=" << m.height << " delta_e=" << m.delta_extruder
                      << " internal=" << m.internal_only << std::endl;
    }
    std::cout << "Deposition minimum width=" << minimum_width << " height=" << minimum_height
              << " stationary=" << degenerate << " unresolved_size=" << undersized << std::endl;
    for(std::size_t i=1;i<moves.size();++i)
        if(moves[i].type==EMoveType::Extrude && std::abs(moves[i].position.z()-moves[i-1].position.z())>1e-5) {
            if(slopes++<5)
                std::cout << "Deposition sloped move: index=" << i << " gcode_id=" << moves[i].gcode_id
                          << " start=" << moves[i-1].position.transpose() << " end=" << moves[i].position.transpose()
                          << " width=" << moves[i].width << " height=" << moves[i].height << std::endl;
        }
    std::cout << "Deposition sloped move count=" << slopes << std::endl;
    try {
        auto s=DepositionMesh::snapshot(processor.get_result(),options);
        if(const char *region=std::getenv("ORCA_DEPOSITION_REGION")) {
            std::string values(region);
            std::replace(values.begin(),values.end(),',',' ');
            std::istringstream input(values);
            std::array<int,3> tile{};
            if(!(input>>tile[0]>>tile[1]>>tile[2])) throw std::invalid_argument("Expected diagnostic tile x,y,z");
            const double h=options.voxel_size_mm,band=3*h;
            s.segments.erase(std::remove_if(s.segments.begin(),s.segments.end(),[&](const auto &segment) {
                const double r=segment.width*.5+band;
                for(unsigned axis=0;axis<3;++axis) {
                    const double low=std::min(segment.start[axis],segment.end[axis])-(axis==2?segment.height+band:r);
                    const double high=std::max(segment.start[axis],segment.end[axis])+(axis==2?band:r);
                    if(high<(tile[axis]*int(options.tile_size_cells)-3)*h ||
                       low>((tile[axis]+1)*int(options.tile_size_cells)+3)*h) return true;
                }
                return false;
            }),s.segments.end());
        }
        if(std::getenv("ORCA_DEPOSITION_MINIMIZE")) {
            std::array<std::array<float,3>,3> target{};
            bool have_target=false;
            const auto reproduces=[&](const DepositionMesh::Snapshot &candidate) {
                CountingBuffer discard;
                std::ostream stream(&discard);
                try { DepositionMesh::export_stl(candidate,stream,options); }
                catch(const DepositionMesh::DegenerateTriangleError &error) {
                    if(!have_target) { target=error.vertices; have_target=true; }
                    return error.vertices==target && !error.planar_reduced;
                } catch(const std::exception &) {}
                return false;
            };
            REQUIRE(reproduces(s));
            unsigned attempts=0;
            for(std::size_t chunk=s.segments.size()/2;chunk>0 && attempts<256;chunk/=2) {
                for(std::size_t begin=0;begin<s.segments.size() && attempts<256;) {
                    auto candidate=s;
                    candidate.segments.erase(candidate.segments.begin()+std::ptrdiff_t(begin),
                        candidate.segments.begin()+std::ptrdiff_t(std::min(begin+chunk,candidate.segments.size())));
                    ++attempts;
                    if(!candidate.segments.empty() && reproduces(candidate)) s=std::move(candidate);
                    else begin+=chunk;
                }
            }
            std::cout << "Minimized segments=" << s.segments.size() << " attempts=" << attempts << std::endl;
            if(s.segments.size()<=32) for(const auto &segment:s.segments)
                std::cout << std::setprecision(9) << "Segment " << segment.start.transpose() << " -> "
                          << segment.end.transpose() << " width=" << segment.width << " height=" << segment.height << std::endl;
        }
        std::cout << "Deposition benchmark: segments=" << s.segments.size()
                  << " extruder_only=" << s.unrepresented_extruder_only_moves
                  << " unrepresented_motion_extrusion=" << s.unrepresented_motion_extrusion_moves
                  << " unrepresented_stationary_extrusion=" << s.unrepresented_stationary_extrusion_moves
                  << " resolution_mm=" << options.voxel_size_mm << std::endl;
        const auto report=DepositionMesh::export_stl(s,sink,options,[last=-1,&buffer,&start](int progress) mutable {
            if(progress%5==0 && progress!=last) {
                std::cout << "Deposition benchmark progress=" << progress << " bytes=" << buffer.bytes
                          << " seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
                          << std::endl;
                last=progress;
            }
        });
        std::cout << "Deposition benchmark: triangles=" << report.triangle_count
                  << " bytes=" << report.output_bytes << " sample_visits=" << report.sample_visits << std::endl;
        REQUIRE(buffer.bytes==report.output_bytes);
    } catch(const std::exception &error) {
        if(const auto *triangle=dynamic_cast<const DepositionMesh::DegenerateTriangleError *>(&error)) {
            Eigen::Vector3d p[3];
            std::cout << std::setprecision(17) << "Degenerate triangle reduced=" << triangle->planar_reduced;
            for(unsigned i=0;i<3;++i) {
                for(unsigned k=0;k<3;++k) p[i][k]=triangle->vertices[i][k];
                std::cout << " p" << i << "=" << p[i].transpose();
            }
            std::cout << " double_area_squared=" << (p[1]-p[0]).cross(p[2]-p[0]).squaredNorm() << std::endl;
        }
        std::cout << "Deposition benchmark stopped: " << error.what()
                  << " bytes_before_stop=" << buffer.bytes << std::endl;
        FAIL_CHECK(error.what());
    }
    std::cout << "Deposition benchmark seconds="
              << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count() << std::endl;
}
