#include "DepositionMesh.hpp"
#include "GCodeProcessor.hpp"
#include "libslic3r/I18N.hpp"
#include "libslic3r/AABBTreeIndirect.hpp"
#include "libslic3r/TextureBake/TextureBakeDecimate.hpp"
#include "libslic3r/TextureBake/TextureBakeIndex.hpp"
#include <openvdb/openvdb.h>
#include <openvdb/Types.h>
#include <openvdb/math/Transform.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <ios>
#include <iosfwd>
#include <vector>
#include <limits>
#include <map>
#include <ostream>
#include <stdexcept>
#include <utility>

namespace Slic3r::DepositionMesh {
DegenerateTriangleError::DegenerateTriangleError(const std::array<std::array<float,3>,3> &points,bool reduced)
    : std::runtime_error(_u8L("Could not construct a valid deposition surface at the selected resolution.")),
      vertices(points),planar_reduced(reduced) {}
namespace {
void validate(const Options &o)
{
    if (!std::isfinite(o.voxel_size_mm) || o.voxel_size_mm < .005 || o.voxel_size_mm > 1. ||
        o.max_memory_bytes < 1024 * 1024 || o.max_memory_bytes > 1024ULL * 1024 * 1024 ||
        o.max_output_bytes < 84 || o.max_output_bytes > 64ULL * 1024 * 1024 * 1024 ||
        o.max_segments == 0 || o.max_segments > 2000000 || o.max_sample_visits == 0 ||
        o.max_sample_visits > 50000000000ULL || o.max_tiles==0 || o.max_tiles>262144 ||
        (o.tile_size_cells!=0 && (o.tile_size_cells<16 || o.tile_size_cells>128)))
        throw std::invalid_argument(_u8L("Invalid deposition mesh resolution or resource limits."));
}
void validate(const Segment &s, double resolution)
{
    if (!s.start.allFinite() || !s.end.allFinite() || !std::isfinite(s.width) || !std::isfinite(s.height) ||
        s.height <= 0.f || s.width < s.height)
        throw std::runtime_error(_u8L("G-code contains undefined bead dimensions or a bead narrower than its height."));
    if ((s.end - s.start).head<2>().squaredNorm() < 1e-12f)
        throw std::runtime_error(_u8L("Stationary extrusion has no reconstructable spatial bead path."));
    if (resolution > .5 * std::min(double(s.width), double(s.height)) + 1e-7)
        throw std::runtime_error(_u8L("Choose a voxel size no larger than half the smallest extrusion width or layer height."));
    if (s.width / resolution > 100000000. || s.height / resolution > 100000000. ||
        s.start.cwiseAbs().maxCoeff() / resolution > 100000000. ||
        s.end.cwiseAbs().maxCoeff() / resolution > 100000000.)
        throw std::runtime_error(_u8L("G-code coordinates exceed the deposition mesh grid range."));
}
void u32(std::ostream &out, std::uint32_t v)
{
    const char b[4] = {char(v), char(v >> 8), char(v >> 16), char(v >> 24)};
    out.write(b, 4);
}
double area_squared(const openvdb::Vec3s &a, const openvdb::Vec3s &b, const openvdb::Vec3s &c)
{
    const auto normal=(openvdb::Vec3d(b)-openvdb::Vec3d(a)).cross(openvdb::Vec3d(c)-openvdb::Vec3d(a));
    return normal.lengthSqr();
}
bool alternate_diagonal(const openvdb::Vec3s &a, const openvdb::Vec3s &b,
                        const openvdb::Vec3s &c, const openvdb::Vec3s &d)
{
    // Dual contouring can put three vertices of a valid quad on one line.
    // Keep all four boundary edges: discarding the zero-area triangle would
    // create a T-junction; splitting the other diagonal preserves that vertex.
    const double ac=std::min(area_squared(a,b,c),area_squared(a,c,d));
    const double bd=std::min(area_squared(a,b,d),area_squared(b,c,d));
    return bd>ac;
}
struct QuadSplit {
    std::array<unsigned,6> indices{};
    unsigned triangles=0;
};
using SurfacePoint=std::array<float,3>;
using SurfaceEdge=std::array<SurfacePoint,2>;
using SurfaceTriangle=std::array<openvdb::Vec3s,3>;
SurfacePoint point_key(const openvdb::Vec3s &p) { return {p.x(),p.y(),p.z()}; }
SurfaceEdge edge_key(SurfacePoint a,SurfacePoint b)
{
    if(b<a) std::swap(a,b);
    return {a,b};
}
bool strictly_on_edge(const SurfacePoint &p,const SurfacePoint &a,const SurfacePoint &b)
{
    if(p==a || p==b) return false;
    for(unsigned k=0;k<3;++k) if(p[k]<std::min(a[k],b[k]) || p[k]>std::max(a[k],b[k])) return false;
    const openvdb::Vec3d pa(double(p[0])-a[0],double(p[1])-a[1],double(p[2])-a[2]);
    const openvdb::Vec3d ba(double(b[0])-a[0],double(b[1])-a[1],double(b[2])-a[2]);
    return pa.cross(ba).lengthSqr()==0.;
}
QuadSplit split_quad(const openvdb::Vec3s &a, const openvdb::Vec3s &b,
                     const openvdb::Vec3s &c, const openvdb::Vec3s &d)
{
    const std::array<openvdb::Vec3s,4> p{a,b,c,d};
    std::array<unsigned,4> boundary{};
    unsigned n=0;
    for(unsigned i=0;i<4;++i)
        if(n==0 || p[i]!=p[boundary[n-1]]) boundary[n++]=i;
    if(n>1 && p[boundary[0]]==p[boundary[n-1]]) --n;
    // Exact float-coordinate coincidences can occur at a sampled bead pole.
    // Collapse only zero-length edges; never epsilon-weld distinct geometry.
    if(n<3) return {};
    if(n==3) return {{boundary[0],boundary[1],boundary[2],0,0,0},1};
    if(alternate_diagonal(a,b,c,d)) return {{0,1,3,1,2,3},2};
    return {{0,1,2,0,2,3},2};
}
void triangle(std::ostream &out, const openvdb::Vec3s &a, const openvdb::Vec3s &b, const openvdb::Vec3s &c,bool reduced=false)
{
    auto normal = (openvdb::Vec3d(b)-openvdb::Vec3d(a)).cross(openvdb::Vec3d(c)-openvdb::Vec3d(a));
    const double length = normal.length();
    if (!(length > 0.) || !std::isfinite(length))
        throw DegenerateTriangleError({{{a.x(),a.y(),a.z()},{b.x(),b.y(),b.z()},{c.x(),c.y(),c.z()}}},reduced);
    normal /= length;
    std::array<char,50> record{};
    std::size_t offset=0;
    for (const auto &v : {openvdb::Vec3s(normal), a, b, c})
        for (int axis = 0; axis < 3; ++axis) {
            const float value=v[axis];
            std::uint32_t bits;
            static_assert(sizeof(bits)==sizeof(value));
            std::memcpy(&bits,&value,sizeof(bits));
            for(unsigned byte=0;byte<4;++byte) record[offset++]=char(bits>>(8*byte));
        }
    out.write(record.data(),std::streamsize(record.size()));
}
} // namespace

Snapshot snapshot(const GCodeProcessorResult &gcode, const Options &o, const std::function<void()> &check)
{
    validate(o);
    std::size_t count = 0;
    Snapshot result;
    for (std::size_t i = 1; i < gcode.moves.size(); ++i) {
        if (check && i % 4096 == 0) check();
        const auto &m = gcode.moves[i];
        if (m.type == EMoveType::Extrude) {
            if (m.position.allFinite() && gcode.moves[i-1].position.allFinite() &&
                (m.position-gcode.moves[i-1].position).head<2>().squaredNorm()<1e-12f) {
                if(m.delta_extruder>0.f) {
                    if(std::abs(double(m.position.z())-gcode.moves[i-1].position.z())<1e-5)
                        ++result.unrepresented_stationary_extrusion_moves;
                    else ++result.unrepresented_motion_extrusion_moves;
                }
                continue;
            }
            validate(Segment{gcode.moves[i-1].position, m.position, m.width, m.height}, o.voxel_size_mm);
            if (++count > o.max_segments || count > o.max_memory_bytes / (2 * sizeof(Segment)))
                throw std::runtime_error(_u8L("Deposition snapshot exceeds the memory or segment budget."));
        } else if (m.type == EMoveType::Unretract && m.delta_extruder > 0.f)
            ++result.unrepresented_extruder_only_moves;
        else if ((m.type == EMoveType::Travel || m.type == EMoveType::Wipe) && m.delta_extruder > 0.f)
            ++result.unrepresented_motion_extrusion_moves;
    }
    if (!count) throw std::runtime_error(_u8L("There are no spatial extrusion paths to reconstruct."));
    result.segments.reserve(count);
    for (std::size_t i = 1; i < gcode.moves.size(); ++i) {
        if (check && i % 4096 == 0) check();
        const auto &m = gcode.moves[i];
        if (m.type == EMoveType::Extrude &&
            (m.position-gcode.moves[i-1].position).head<2>().squaredNorm()>=1e-12f)
            result.segments.push_back({gcode.moves[i-1].position, m.position, m.width, m.height});
    }
    if (check) check();
    return result;
}

namespace {
struct Tile {
    std::array<int,3> index, low, high;
};
void write_header(std::ostream &out)
{
    std::array<char,80> header{};
    const char label[]="OrcaSlicer nominal G-code bead union; millimetres; no material data";
    std::copy(label,label+sizeof(label)-1,header.begin());
    out.write(header.data(),header.size()); u32(out,0);
}

void export_tile(const Snapshot &input, const std::vector<std::size_t> *selected, const Tile *tile,
                 std::size_t persistent_bytes, std::ostream &out, const Options &o,
                 std::uint64_t &visits, std::uint64_t &triangles_written,
                 const std::function<void(int)> &progress, const std::function<void()> &check)
{
    if(check) check();
    if(persistent_bytes>o.max_memory_bytes-65536)
        throw std::runtime_error(_u8L("Deposition distance field exceeds the memory budget."));
    const std::size_t count_segments=selected?selected->size():input.segments.size();
    const double h = o.voxel_size_mm, band = 3. * h;
    openvdb::initialize();
    auto grid = openvdb::FloatGrid::create(float(band));
    grid->setTransform(openvdb::math::Transform::createLinearTransform(h));
    grid->setGridClass(openvdb::GRID_LEVEL_SET);
    auto accessor = grid->getAccessor();
    using Tree = openvdb::FloatTree;
    using Upper = Tree::RootNodeType::ChildNodeType;
    using Lower = Upper::ChildNodeType;
    using Leaf = Tree::LeafNodeType;
    // Account for tree allocations before inserting a new path, including
    // conservative allocator/root-map overhead. Meshing gets a separate bound.
    std::size_t tree_budget = 65536;
    int last_progress = -1;
    const auto tick = [&] {
        if (++visits > o.max_sample_visits) throw std::runtime_error(_u8L("Deposition reconstruction exceeds the sampling work budget."));
        if (check && visits % 4096 == 0) check();
    };
    for (std::size_t si = 0; si < count_segments; ++si) {
        const auto &s = input.segments[selected?(*selected)[si]:si];
        const double dx = double(s.end.x())-s.start.x(), dy = double(s.end.y())-s.start.y();
        const double dz = double(s.end.z())-s.start.z();
        const bool along_x = std::abs(dx) >= std::abs(dy);
        const int major = along_x ? 0 : 1, minor = 1-major;
        const double delta_major = double(s.end[major])-s.start[major];
        const double delta_minor = double(s.end[minor])-s.start[minor];
        const double radius = double(s.width)*.5 + band;
        int begin = int(std::floor((std::min(s.start[major],s.end[major])-radius)/h));
        int end = int(std::ceil((std::max(s.start[major],s.end[major])+radius)/h));
        if(tile) { begin=std::max(begin,tile->low[major]); end=std::min(end,tile->high[major]); }
        const double length2 = dx*dx+dy*dy, r = double(s.height)*.5;
        const double length3_squared=length2+dz*dz;
        const double parameter_tolerance=h*1e-4/std::sqrt(length3_squared);
        const auto half_ulp=[](float value) {
            return .5*std::max(double(std::nextafter(value,std::numeric_limits<float>::infinity()))-value,
                              double(value)-std::nextafter(value,-std::numeric_limits<float>::infinity()));
        };
        const double flat_z_uncertainty=half_ulp(s.start.z())+half_ulp(s.height)+
            16*std::numeric_limits<double>::epsilon()*(std::abs(double(s.start.z()))+s.height+h);
        for (int a = begin; a <= end; ++a) {
            tick();
            const double t = std::clamp((a*h-s.start[major])/delta_major,0.,1.);
            const double center = s.start[minor]+t*delta_minor;
            // A diagonal strip, not its potentially enormous enclosing box.
            const double strip = radius * (1.+std::abs(delta_minor/delta_major));
            int b0 = int(std::floor((center-strip)/h)), b1 = int(std::ceil((center+strip)/h));
            if(tile) { b0=std::max(b0,tile->low[minor]); b1=std::min(b1,tile->high[minor]); }
            // Only the segment portion within this XY strip can contribute.
            // Bound its Z extent as well, rather than scanning the full ramp's
            // enclosing box for every sample along a sloping/scarf extrusion.
            const double t_radius=radius/std::abs(delta_major);
            const double za=s.start.z()+std::clamp(t-t_radius,0.,1.)*dz;
            const double zb=s.start.z()+std::clamp(t+t_radius,0.,1.)*dz;
            int z0=int(std::floor((std::min(za,zb)-s.height-band)/h));
            int z1=int(std::ceil((std::max(za,zb)+band)/h));
            if(tile) { z0=std::max(z0,tile->low[2]); z1=std::min(z1,tile->high[2]); }
            for (int b = b0; b <= b1; ++b) for (int z = z0; z <= z1; ++z) {
                tick();
                const double x = (along_x?a:b)*h, y = (along_x?b:a)*h;
                double u = std::clamp(((x-s.start.x())*dx+(y-s.start.y())*dy)/length2,0.,1.);
                const double disc_radius=(double(s.width)-s.height)*.5;
                if (dz!=0.) {
                    const double qx=x-s.start.x(),qy=y-s.start.y(),qz=z*h-(double(s.start.z())-r);
                    const double nearest=std::clamp((qx*dx+qy*dy+qz*dz)/length3_squared,0.,1.);
                    const double px=qx-nearest*dx,py=qy-nearest*dy,pz=qz-nearest*dz;
                    // Distance to the segment minus the disc radius is a
                    // lower bound. Reject distant samples before minimization.
                    if(px*px+py*py+pz*pz>radius*radius) continue;
                    // Distance to segment + horizontal disc + sphere. The
                    // squared distance to the swept disc is convex in t. Its
                    // monotone derivative permits bounded bisection, while
                    // horizontal paths retain their exact projection fastpath.
                    const auto derivative=[&](double parameter) {
                        const double qx=x-s.start.x()-parameter*dx, qy=y-s.start.y()-parameter*dy;
                        const double radial=std::hypot(qx,qy);
                        const double qz=z*h-(double(s.start.z())-r)-parameter*dz;
                        return -dz*qz-(radial>disc_radius ? (1.-disc_radius/radial)*(dx*qx+dy*qy) : 0.);
                    };
                    if (derivative(0.)>=0.) u=0.;
                    else if (derivative(1.)<=0.) u=1.;
                    else {
                        double lo=0.,hi=1.;
                        for(unsigned iteration=0;iteration<32 && hi-lo>parameter_tolerance;++iteration) {
                            tick();
                            const double mid=.5*(lo+hi);
                            if (derivative(mid)>0.) hi=mid; else lo=mid;
                        }
                        u=.5*(lo+hi);
                    }
                }
                const double planar = std::hypot(x-s.start.x()-u*dx,y-s.start.y()-u*dy);
                const double side = std::max(0.,planar-disc_radius);
                double distance = std::hypot(side,z*h-(double(s.start.z())+u*dz-r))-r;
                if (distance >= band) continue;
                // Include exact boundary samples in the solid; adjacent flat
                // layers share an interior plane, not two internal STL skins.
                if (std::abs(distance) < h*1e-7) distance = -h*1e-7;
                // Parsed coordinates and heights are floats. Independently
                // rounded touching layers can straddle a lattice plane and
                // otherwise create false internal skins. Snap only the flat
                // horizontal face within its input rounding uncertainty;
                // curved surfaces and resolvable gaps are not expanded.
                if (dz==0. && disc_radius>0. && side==0. && std::abs(distance)<=flat_z_uncertainty)
                    distance=-h*1e-7;
                const openvdb::Coord p(along_x?a:b,along_x?b:a,z);
                const float value = float(std::max(-band,distance));
                if (value >= accessor.getValue(p)) continue;
                if (!grid->tree().probeLeaf(p)) {
                    std::size_t growth = 2*(sizeof(Leaf)+Leaf::SIZE*sizeof(float)+128);
                    if (!grid->tree().probeNode<Lower>(p)) growth += 2*sizeof(Lower);
                    if (!grid->tree().probeNode<Upper>(p)) growth += 2*sizeof(Upper);
                    if (growth > o.max_memory_bytes-persistent_bytes ||
                        tree_budget > o.max_memory_bytes-persistent_bytes-growth)
                        throw std::runtime_error(_u8L("Deposition distance field exceeds the memory budget."));
                    tree_budget += growth;
                }
                if (distance <= -band) accessor.setValueOff(p,value);
                else accessor.setValue(p,value);
            }
        }
        const int percent = int(60*(si+1)/count_segments);
        if (progress && percent != last_progress) { progress(percent); last_progress=percent; }
    }
    // Every intersected cell's lower corner lies within the 3-voxel active
    // band. Count these before non-interruptible VolumeToMesh allocations.
    std::uint64_t surface_cells = 0;
    for (auto it = grid->cbeginValueOn(); it; ++it) {
        tick();
        const auto p = it.getCoord();
        bool inside = false, outside = false;
        for (int z=0;z<2;++z) for(int y=0;y<2;++y) for(int x=0;x<2;++x) {
            const bool sign = accessor.getValue(p.offsetBy(x,y,z)) < 0.f;
            inside |= sign; outside |= !sign;
        }
        if (inside && outside) ++surface_cells;
    }
    // Clipping a tile's field creates additional artificial surfaces at its
    // halo boundary. Budget those too; ownership excludes them from the STL.
    if(tile) surface_cells+=6*std::uint64_t(o.tile_size_cells+8)*(o.tile_size_cells+8);
    if (!surface_cells) return;
    // Nonadaptive dual contouring has at most four vertices/cell and three
    // edge quads/cell. Use a larger bound for auxiliary topology and pools.
    const std::uint64_t remaining = o.max_memory_bytes-persistent_bytes;
    if (tree_budget > remaining/8 || surface_cells > (remaining-8*tree_budget)/256)
        throw std::runtime_error(_u8L("Deposition surface meshing exceeds the conservative memory budget."));
    if (surface_cells > (o.max_output_bytes-84)/600)
        throw std::runtime_error(_u8L("Deposition STL exceeds the conservative output-size budget."));
    if (check) check();
    if (progress) progress(65);
    openvdb::tools::VolumeToMesh mesher(0.,0.,false);
    mesher(*grid);
    if (check) check();
    if (progress) progress(85);
    std::uint64_t count = 0;
    const auto &points = mesher.pointList();
    const auto owns=[&](const auto &indices,unsigned n) {
        if(!tile) return true;
        for(unsigned axis=0;axis<3;++axis) {
            double center=0.;
            for(unsigned k=0;k<n;++k) center+=double(points[indices[k]][axis]);
            if(std::floor(center/(n*o.tile_size_cells*h))!=tile->index[axis]) return false;
        }
        return true;
    };
    const auto visit_raw=[&](const auto &visit) {
        for(std::size_t i=0;i<mesher.polygonPoolListSize();++i) {
            if(check) check();
            const auto &pool=mesher.polygonPoolList()[i];
            for(std::size_t j=0;j<pool.numTriangles();++j) {
                tick();
                const auto &t=pool.triangle(j);
                visit(SurfaceTriangle{points[t[0]],points[t[1]],points[t[2]]},owns(t,3));
            }
            for(std::size_t j=0;j<pool.numQuads();++j) {
                tick();
                const auto &q=pool.quad(j);
                const bool owned=owns(q,4);
                const auto split=split_quad(points[q[0]],points[q[1]],points[q[2]],points[q[3]]);
                for(unsigned t=0;t<split.triangles;++t)
                    visit(SurfaceTriangle{points[q[split.indices[3*t]]],points[q[split.indices[3*t+1]]],
                                          points[q[split.indices[3*t+2]]]},owned);
            }
        }
    };
    // A collinear VDB face can stitch AC to AB+BC. Removing it alone creates
    // a T-junction. Insert its exact vertices in adjacent faces across the
    // full halo, retaining original polygon ownership through every split.
    std::size_t zero_faces=0;
    visit_raw([&](const SurfaceTriangle &t,bool) { if(area_squared(t[0],t[1],t[2])==0.) ++zero_faces; });
    const std::uint64_t repair_available=remaining-8*tree_budget-256*surface_cells;
    constexpr std::size_t face_storage=2*(sizeof(SurfaceTriangle)+3*sizeof(SurfacePoint));
    if(zero_faces>repair_available/face_storage)
        throw std::runtime_error(_u8L("Deposition surface repair exceeds the memory budget."));
    std::uint64_t repair_bytes=zero_faces*face_storage;
    std::vector<SurfaceTriangle> stitches;
    std::vector<SurfacePoint> cuts;
    stitches.reserve(zero_faces); cuts.reserve(zero_faces*3);
    visit_raw([&](const SurfaceTriangle &t,bool) {
        if(area_squared(t[0],t[1],t[2])==0.) {
            stitches.push_back(t);
            for(const auto &p:t) cuts.push_back(point_key(p));
        }
    });
    std::sort(cuts.begin(),cuts.end());
    cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    const auto cut_range=[&](const SurfacePoint &a,const SurfacePoint &b) {
        const auto key=edge_key(a,b);
        const auto compare=[&](const SurfacePoint &lhs,const SurfacePoint &rhs) { tick(); return lhs<rhs; };
        const auto begin=std::upper_bound(cuts.begin(),cuts.end(),key[0],compare);
        return std::make_pair(begin,std::lower_bound(begin,cuts.end(),key[1],compare));
    };
    struct Incidence { unsigned count=0; int orientation=0; };
    std::map<SurfaceEdge,Incidence> affected;
    const auto watch=[&](const SurfacePoint &a,const SurfacePoint &b) {
        if(a==b) return;
        const auto key=edge_key(a,b);
        if(affected.find(key)==affected.end()) {
            if(repair_bytes>repair_available || 256>repair_available-repair_bytes)
                throw std::runtime_error(_u8L("Deposition surface repair exceeds the memory budget."));
            repair_bytes+=256;
            affected.emplace(key,Incidence{});
        }
    };
    for(const auto &t:stitches) for(unsigned edge=0;edge<3;++edge) {
        const auto key=edge_key(point_key(t[edge]),point_key(t[(edge+1)%3]));
        auto previous=key[0];
        const auto range=cut_range(key[0],key[1]);
        for(auto it=range.first;it!=range.second;++it) {
            tick();
            const auto &p=*it;
            if(strictly_on_edge(p,key[0],key[1])) { watch(previous,p); previous=p; }
        }
        watch(previous,key[1]);
    }
    const auto split_edges=[&](const auto &self,const SurfaceTriangle &t,unsigned depth,const auto &emit)->void {
        tick();
        for(unsigned edge=0;edge<3;++edge) {
            const auto a=point_key(t[edge]),b=point_key(t[(edge+1)%3]);
            const auto range=cut_range(a,b);
            for(auto it=range.first;it!=range.second;++it) {
                tick();
                const auto &p=*it;
                if(!strictly_on_edge(p,a,b)) continue;
                if(depth>=64) throw std::runtime_error(_u8L("Could not resolve deposition surface topology within its work limits."));
                const openvdb::Vec3s middle(p[0],p[1],p[2]);
                self(self,SurfaceTriangle{t[edge],middle,t[(edge+2)%3]},depth+1,emit);
                self(self,SurfaceTriangle{middle,t[(edge+1)%3],t[(edge+2)%3]},depth+1,emit);
                return;
            }
        }
        emit(t);
    };
    visit_raw([&](const SurfaceTriangle &t,bool owned) {
        if(area_squared(t[0],t[1],t[2])==0.) return;
        split_edges(split_edges,t,0,[&](const SurfaceTriangle &part) {
            if(owned) ++count;
            for(unsigned e=0;e<3;++e) {
                const auto a=point_key(part[e]),b=point_key(part[(e+1)%3]);
                const auto it=affected.find(edge_key(a,b));
                if(it!=affected.end()) { ++it->second.count; it->second.orientation+=a<b?1:-1; }
            }
        });
    });
    for(const auto &edge:affected) if(edge.second.count!=2 || edge.second.orientation!=0)
        throw std::runtime_error(_u8L("Could not resolve deposition surface topology at the selected resolution."));
    if(!count) return;
    const auto check_output=[&]() {
        if (count > std::numeric_limits<std::uint32_t>::max()-triangles_written ||
            count > (o.max_output_bytes-84)/50-triangles_written)
            throw std::runtime_error(_u8L("Deposition STL exceeds its triangle or output-size budget."));
    };
    TextureBake::TriSoup raw;
    if(o.reduce_horizontal_faces) {
        // Includes input/output soups, exact indexing, incidence, quadrics and
        // the strictly capped 3F candidate heap, alongside the live VDB mesh.
        if(count>(repair_available-repair_bytes)/2000)
            throw std::runtime_error(_u8L("Deposition planar reduction exceeds the memory budget."));
        raw.pos.reserve(std::size_t(count)*3);
    } else {
        check_output();
        if(triangles_written==0) write_header(out);
    }
    const auto emit=[&](const auto &a,const auto &b,const auto &c) {
        if(o.reduce_horizontal_faces) {
            if(area_squared(a,b,c)==0.)
                throw DegenerateTriangleError({{{a.x(),a.y(),a.z()},{b.x(),b.y(),b.z()},{c.x(),c.y(),c.z()}}},false);
            for(const auto &p:{a,b,c}) raw.pos.emplace_back(p.x(),p.y(),p.z());
        } else triangle(out,a,b,c);
    };
    visit_raw([&](const SurfaceTriangle &t,bool owned) {
        if(!owned || area_squared(t[0],t[1],t[2])==0.) return;
        split_edges(split_edges,t,0,[&](const SurfaceTriangle &part) {
            emit(part[0],part[1],part[2]);
            if(!out) throw std::runtime_error(_u8L("Failed writing deposition STL; check free disk space."));
        });
    });
    if(o.reduce_horizontal_faces) {
        const auto reduced=TextureBake::decimate(raw,std::size_t(count),true,0.,{},
            [&](double) { if(check) check(); return true; },{},true);
        count=reduced.geometry.pos.size()/3;
        check_output();
        if(triangles_written==0) write_header(out);
        for(std::size_t i=0;i<reduced.geometry.pos.size();i+=3) {
            if(check && i%4095==0) check();
            const auto point=[&](std::size_t j) {
                const auto &p=reduced.geometry.pos[j]; return openvdb::Vec3s(p.x(),p.y(),p.z());
            };
            triangle(out,point(i),point(i+1),point(i+2),true);
        }
        if(!out) throw std::runtime_error(_u8L("Failed writing deposition STL; check free disk space."));
    }
    triangles_written+=count;
}
} // namespace

Report export_stl(const Snapshot &input, std::ostream &out, const Options &o,
                  const std::function<void(int)> &progress, const std::function<void()> &check)
{
    validate(o);
    if(check) check();
    if(input.segments.empty() || input.segments.size()>o.max_segments ||
       input.segments.size()>o.max_memory_bytes/(2*sizeof(Segment)))
        throw std::runtime_error(_u8L("Deposition snapshot is empty or exceeds its resource budget."));
    for(std::size_t i=0;i<input.segments.size();++i) {
        if(check && i%4096==0) check();
        validate(input.segments[i],o.voxel_size_mm);
    }
    const auto start=out.tellp();
    if(start==std::streampos(-1)) throw std::runtime_error(_u8L("Deposition STL requires a seekable private output stream."));
    const std::size_t snapshot_bytes=input.segments.capacity()*sizeof(Segment);
    if(snapshot_bytes>=o.max_memory_bytes) throw std::runtime_error(_u8L("Deposition snapshot exceeds the memory budget."));
    std::uint64_t visits=0,triangles=0;
    if(o.tile_size_cells==0) {
        // Whole-grid reference used to verify tile-boundary equivalence.
        export_tile(input,nullptr,nullptr,snapshot_bytes,out,o,visits,triangles,progress,check);
    } else {
        using BVH=AABBTreeIndirect::Tree<3,float>;
        struct Source {
            std::size_t id;
            BVH::BoundingBox box;
            std::size_t idx() const { return id; }
            const BVH::BoundingBox &bbox() const { return box; }
            BVH::VectorType centroid() const { return box.center(); }
        };
        std::size_t power=1;
        while(power<input.segments.size()) power*=2;
        const std::size_t index_bound=(2*power-1)*sizeof(BVH::Node);
        if(index_bound+input.segments.size()*sizeof(Source)>o.max_memory_bytes-snapshot_bytes)
            throw std::runtime_error(_u8L("Deposition spatial index exceeds the memory budget."));
        BVH tree;
        {
            std::vector<Source> source;
            source.reserve(input.segments.size());
            const float band=float(3*o.voxel_size_mm);
            for(std::size_t i=0;i<input.segments.size();++i) {
                if(check && i%4096==0) check();
                const auto &s=input.segments[i];
                const float r=s.width*.5f+band;
                BVH::VectorType low=s.start.cwiseMin(s.end)-Eigen::Vector3f(r,r,s.height+band);
                BVH::VectorType high=s.start.cwiseMax(s.end)+Eigen::Vector3f(r,r,band);
                for(unsigned a=0;a<3;++a) {
                    low[a]=std::nextafter(low[a],-std::numeric_limits<float>::infinity());
                    high[a]=std::nextafter(high[a],std::numeric_limits<float>::infinity());
                }
                source.push_back({i,BVH::BoundingBox(low,high)});
            }
            tree.build_modify_input(source);
        }
        if(check) check();
        const std::size_t persistent=snapshot_bytes+tree.nodes().capacity()*sizeof(BVH::Node);
        const auto &bounds=tree.node(0).bbox;
        const double tile_mm=o.tile_size_cells*o.voxel_size_mm;
        std::array<int,3> first{},last{};
        std::size_t tile_count=1;
        for(unsigned a=0;a<3;++a) {
            first[a]=int(std::floor(bounds.min()[a]/tile_mm));
            last[a]=int(std::floor(bounds.max()[a]/tile_mm));
            const std::size_t n=std::size_t(last[a]-first[a])+1;
            if(n>o.max_tiles/tile_count) throw std::runtime_error(_u8L("Deposition reconstruction exceeds the spatial tile budget."));
            tile_count*=n;
        }
        std::vector<std::size_t> selected;
        std::size_t done=0;
        int last_percent=-1;
        for(int z=first[2];z<=last[2];++z) for(int y=first[1];y<=last[1];++y) for(int x=first[0];x<=last[0];++x) {
            if(check) check();
            Tile tile{{x,y,z},{},{}};
            BVH::VectorType low,high;
            for(unsigned a=0;a<3;++a) {
                tile.low[a]=tile.index[a]*int(o.tile_size_cells)-3;
                tile.high[a]=(tile.index[a]+1)*int(o.tile_size_cells)+3;
                low[a]=std::nextafter(float(tile.low[a]*o.voxel_size_mm),-std::numeric_limits<float>::infinity());
                high[a]=std::nextafter(float(tile.high[a]*o.voxel_size_mm),std::numeric_limits<float>::infinity());
            }
            const BVH::BoundingBox query(low,high);
            selected.clear();
            const auto collect=[&](const auto &self,std::size_t index)->void {
                if(++visits>o.max_sample_visits) throw std::runtime_error(_u8L("Deposition reconstruction exceeds the sampling work budget."));
                if(check && visits%4096==0) check();
                const auto &node=tree.node(index);
                if(!node.is_valid() || !node.bbox.intersects(query)) return;
                if(node.is_leaf()) {
                    if(selected.size()==selected.capacity()) {
                        const std::size_t capacity=std::max(std::size_t(64),selected.capacity()*2);
                        if((capacity+selected.capacity())*sizeof(std::size_t)>o.max_memory_bytes-persistent)
                            throw std::runtime_error(_u8L("Deposition tile selection exceeds the memory budget."));
                        selected.reserve(capacity);
                    }
                    selected.push_back(node.idx);
                } else { self(self,BVH::left_child_idx(index)); self(self,BVH::right_child_idx(index)); }
            };
            collect(collect,0);
            if(!selected.empty()) export_tile(input,&selected,&tile,persistent+selected.capacity()*sizeof(std::size_t),
                                              out,o,visits,triangles,{},check);
            ++done;
            const int percent=int(100*done/tile_count);
            if(progress && percent!=last_percent) { progress(percent); last_percent=percent; }
        }
    }
    if(!triangles) throw std::runtime_error(_u8L("No resolved deposition surface was found; choose a finer resolution."));
    if(check) check();
    const auto end=out.tellp();
    out.seekp(start+std::streamoff(80)); u32(out,std::uint32_t(triangles)); out.seekp(end);
    if(!out) throw std::runtime_error(_u8L("Failed writing deposition STL; check free disk space."));
    return {triangles,84+50*triangles,o.voxel_size_mm,visits};
}
} // namespace Slic3r::DepositionMesh
