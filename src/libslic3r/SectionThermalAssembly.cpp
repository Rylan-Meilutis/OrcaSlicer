#include "SectionThermalAssembly.hpp"
#include "SectionGrowth.hpp"
#include "SectionGrowthStream.hpp"
#include "SectionIntegration.hpp"
#include "SectionThermalExposure.hpp"
#include "SectionThermalSkin.hpp"
#include "ThermalHistory.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include "ThermoelasticLocalSections.hpp"
#include "ThermoelasticLocalTopology.hpp"
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticExactSections.hpp"
#include "BoundingBox.hpp"
#include "ClipperUtils.hpp"
#include "ExPolygon.hpp"
#include "Point.hpp"
#include "Polygon.hpp"
#include "libslic3r.h"
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <boost/container_hash/hash.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <unordered_map>
#include <vector>

namespace Slic3r::Thermoelastic {
namespace {
using Index = std::array<coord_t,4>;
std::size_t points(const ExPolygons& polygons)
{
    std::size_t count = 0;
    for (const auto& p : polygons) {
        count += p.contour.points.size();
        for (const auto& h : p.holes) count += h.points.size();
    }
    return count;
}
}

static SectionThermalAssembly assemble_volume(const std::vector<SectionPiece>& pieces, double bed,
    const SectionTopology* topology,
    const std::vector<SectionThermalMaterial>& materials, const Eigen::Vector3d& origin,
    const Eigen::Vector3d& spacing, const SectionThermalAssemblyLimits& limits,
    const std::function<void()>& checkpoint,
    const std::function<std::vector<std::size_t>(const std::function<void(const SectionPiece&)>&,
        const SectionThermalSkin::Horizontal&,const SectionThermalSkin::Vertical&,
        const SectionThermalSkin::ExactHorizontal&,const SectionThermalSkin::ExactVertical&)>& produce = {},
    const ThermalSectionStreamSettings* streaming = nullptr, bool allow_empty_increment = false)
{
    if ((pieces.empty() && !produce) || !origin.allFinite() || !spacing.allFinite() || origin.cwiseAbs().maxCoeff() > 10. ||
        spacing.maxCoeff() > 10. || spacing.minCoeff() < 1e-9 || materials.empty() || materials.size() > 256 ||
        !limits.max_cells || limits.max_cells > 100000 || !limits.max_nodes || limits.max_nodes > 200000 ||
        !limits.max_groups || limits.max_groups > 100000 || !limits.max_points || limits.max_points > 2000000 ||
        !limits.max_work || limits.max_work > 200000000)
        throw std::invalid_argument("Invalid section thermal grid or limits.");
    if (streaming && (!produce || topology || !streaming->max_processed_points || streaming->max_processed_points > 128000000 ||
        !streaming->max_work || streaming->max_work > 2000000000 || !streaming->max_provisional_groups ||
        streaming->max_provisional_groups > 400000))
        throw std::invalid_argument("Invalid streamed thermal integration limits.");
    const unsigned normal_bins = streaming ? streaming->vertical_normal_bins : 0;
    if (normal_bins && (normal_bins < 8 || normal_bins > 4096 || normal_bins % 4))
        throw std::invalid_argument("Invalid thermal surface angular resolution.");
    const std::size_t point_limit = streaming ? streaming->max_processed_points : limits.max_points;
    const std::size_t group_limit = streaming ? streaming->max_provisional_groups : limits.max_groups;
    const std::size_t cell_limit = streaming ? streaming->max_provisional_groups : limits.max_cells;
    const std::size_t work_limit = streaming ? streaming->max_work : limits.max_work;
    for (const auto& material : materials)
        if (!material.conductivity_w_m_k.allFinite() ||
            !material.conductivity_w_m_k.isApprox(material.conductivity_w_m_k.transpose(),1e-12) ||
            material.conductivity_w_m_k.llt().info() != Eigen::Success ||
            !std::isfinite(material.volumetric_capacity_j_m3_k) || material.volumetric_capacity_j_m3_k <= 0. ||
            !std::isfinite(material.initial_temperature_k) || material.initial_temperature_k <= 0.)
            throw std::invalid_argument("Invalid section thermal material properties.");
    std::array<coord_t,3> start, size;
    for (unsigned axis = 0; axis < 3; ++axis) {
        for (double value : {origin[axis],spacing[axis]})
            if (std::abs(value*1e9-std::round(value*1e9)) > 1e-5)
                throw std::invalid_argument("Section thermal grid must use native nanometres.");
        start[axis] = std::llround(origin[axis]*1e9);
        size[axis] = std::llround(spacing[axis]*1e9);
    }
    if (start[2] != std::llround(bed*1e9))
        throw std::invalid_argument("Section thermal grid uses a different bed plane.");
    Eigen::Vector3d native_spacing;
    for (unsigned axis = 0; axis < 3; ++axis) native_spacing[axis] = double(size[axis])*1e-9;
    const auto index = [&](coord_t coordinate,unsigned axis) {
        const coord_t delta = coordinate-start[axis];
        const coord_t value = delta/size[axis]-(delta%size[axis] < 0 ? 1 : 0);
        if (value < -10000000 || value > 10000000) throw ResourceLimit("Section thermal grid range exceeded.");
        return value;
    };
    std::size_t work = 0, input_points = 0;
    std::size_t batched_side_samples = 0, expanded_side_batches = 0;
    enum class WorkStage { Volume, Horizontal, Vertical, Assembly };
    WorkStage stage = WorkStage::Assembly;
    std::array<std::size_t,4> stage_work{};
    const auto charge = [&](std::size_t count) {
        if (checkpoint) checkpoint();
        if (count > work_limit-work)
            throw ResourceLimit("Section thermal integration work exceeded: used=" + std::to_string(work) +
                ", requested=" + std::to_string(count) + ", input_points=" + std::to_string(input_points) +
                ", volume=" + std::to_string(stage_work[0]) + ", horizontal=" + std::to_string(stage_work[1]) +
                ", vertical=" + std::to_string(stage_work[2]) + ", assembly=" + std::to_string(stage_work[3]) +
                ", batched_side_samples=" + std::to_string(batched_side_samples) +
                ", expanded_side_batches=" + std::to_string(expanded_side_batches) + ".");
        work += count;
        stage_work[std::size_t(stage)] += count;
    };
    std::map<std::pair<Index,std::size_t>,SectionIntegration::Moments> groups;
    using GroupKey = std::pair<Index,std::size_t>;
    std::map<std::pair<GroupKey,SectionHorizontalFaceKind>,SectionIntegration::Moments> faces;
    struct SideFace { GroupKey group; Eigen::Vector2d normal; SectionIntegration::Moments area; };
    std::vector<SideFace> sides;
    // Streaming bands split the same cell surface many times. Preserve all
    // area moments but combine identical material/component/cell/normal groups
    // rather than retaining one operator per quadrature band.
    using SideKey = std::pair<GroupKey,std::array<double,2>>;
    std::map<SideKey,SectionIntegration::Moments> streamed_sides;
    std::map<coord_t,std::vector<std::size_t>> bottoms, tops;
    const bool exact_exposure = topology && std::any_of(pieces.begin(),pieces.end(),
        [](const SectionPiece& p) { return p.slab.footprint.is_exact(); });
    const auto exact_budget = [&] {
        const auto remaining=work_limit-work;
        if (!remaining) throw ResourceLimit("Exact thermal geometry budget exhausted.");
        return std::min<std::size_t>(4000000,remaining);
    };
    std::vector<ExactSection> exact_sections;
    std::vector<ExactSectionBoundary> exact_boundaries;
    std::size_t retained_exact_points=0;
    const auto retain_exact = [&](std::size_t n) {
        if (n > limits.max_points-retained_exact_points)
            throw ResourceLimit("Exact thermal live geometry budget exceeded.");
        retained_exact_points+=n;
    };
    if (exact_exposure) {
        for (const auto& piece : pieces) {
            charge(1);
            ExactSection section;
            if (piece.slab.footprint.is_exact()) section=piece.slab.footprint.exact();
            else {
                const auto promoted=make_exact_section(piece.slab.footprint.native(),exact_budget(),[&] { charge(0); });
                charge(promoted.work);
                section=promoted.section;
            }
            retain_exact(section.vertex_count());
            const auto boundary=extract_exact_section_boundary(section,exact_budget(),[&] { charge(0); });
            charge(boundary.work);
            retain_exact(2*boundary.boundary.segment_count());
            exact_sections.push_back(std::move(section));
            exact_boundaries.push_back(boundary.boundary);
        }
    }
    for (std::size_t i = 0; topology && i < pieces.size(); ++i) {
        charge(1);
        const auto& slab = pieces[i].slab;
        bottoms[std::llround(slab.bottom_m*1e9)].push_back(i);
        tops[std::llround(slab.top_m*1e9)].push_back(i);
    }
    std::set<Index> cells;
    const auto accept = [&](const SectionPiece& piece) {
        stage = WorkStage::Volume;
        const auto& slab = piece.slab;
        if (slab.material >= materials.size()) throw std::invalid_argument("Missing section thermal material.");
        const bool exact = exact_exposure || slab.footprint.is_exact();
        const auto owner = exact_exposure ? std::size_t(&piece-pieces.data()) : 0;
        const ExactSection* exact_geometry = exact_exposure ? &exact_sections.at(owner) :
            (exact ? &slab.footprint.exact() : nullptr);
        const auto count = slab.footprint.vertex_count();
        if (count > point_limit-input_points) throw ResourceLimit("Section thermal input points exceeded.");
        input_points += count;
        const auto bounds = slab.footprint.bounds();
        if (slab.footprint.empty()) return;
        if (!bounds) throw std::invalid_argument("Missing section thermal footprint bounds.");
        const auto& box = *bounds;
        const coord_t lo_z = std::llround(slab.bottom_m*1e9), hi_z = std::llround(slab.top_m*1e9);
        const auto exposed = [&](const auto& opposing,coord_t plane) {
            if (!topology || exact_exposure) return ExPolygons{};
            ExPolygons face = slab.footprint.native();
            const auto found = opposing.find(plane);
            if (found != opposing.end()) for (const auto neighbor : found->second) {
                charge(1);
                const auto& cover = pieces[neighbor].slab.footprint.native();
                if (!box.overlap(get_extents(cover))) continue;
                charge(points(face)+points(cover));
                face = diff_ex(face,cover);
                if (points(face) > limits.max_points) throw ResourceLimit("Section thermal exposed-face points exceeded.");
                if (face.empty()) break;
            }
            return face;
        };
        const auto lower_face = exposed(tops,lo_z), upper_face = exposed(bottoms,hi_z);
        const auto vertical_patches = topology && !exact_exposure ? exposed_vertical_section_patches(*topology,piece,limits.max_groups,charge) :
            std::vector<VerticalSectionPatch>{};
        ExactSection exact_lower,exact_upper;
        struct ExactPatch { ExactSectionBoundary span; coord_t bottom,top; Eigen::Vector2d normal; };
        std::vector<ExactPatch> exact_patches;
        const auto saved_retained_points=retained_exact_points;
        if (exact_exposure) {
            const auto face = [&](const auto& opposing,coord_t plane) {
                auto result=*exact_geometry;
                const auto found=opposing.find(plane);
                if (found!=opposing.end()) for (const auto neighbor : found->second) {
                    charge(1);
                    const auto old_points=result.vertex_count();
                    if (old_points>limits.max_points-retained_exact_points)
                        throw ResourceLimit("Exact thermal face geometry budget exceeded.");
                    const auto cut=boolean_exact_sections(result,exact_sections[neighbor],ExactSectionOperation::Difference,
                        exact_budget(),[&] { charge(0); });
                    charge(cut.work);
                    if (cut.section.vertex_count()>limits.max_points-retained_exact_points-old_points)
                        throw ResourceLimit("Exact thermal face transition budget exceeded.");
                    result=cut.section;
                    if (result.vertex_count()>limits.max_points-retained_exact_points)
                        throw ResourceLimit("Exact thermal face geometry budget exceeded.");
                    if (result.empty()) break;
                }
                retain_exact(result.vertex_count());
                return result;
            };
            exact_lower=face(tops,lo_z);
            exact_upper=face(bottoms,hi_z);
            std::vector<std::size_t> neighbors;
            std::vector<coord_t> planes{lo_z,hi_z};
            for (std::size_t j=0;j<pieces.size();++j) {
                charge(1);
                if (j==owner || exact_sections[j].empty() || !box.overlap(*exact_sections[j].bounds())) continue;
                const auto lo=std::max(lo_z,coord_t(std::llround(pieces[j].slab.bottom_m*1e9)));
                const auto hi=std::min(hi_z,coord_t(std::llround(pieces[j].slab.top_m*1e9)));
                if (hi<=lo) continue;
                if (neighbors.size()>=limits.max_groups) throw ResourceLimit("Exact thermal contact budget exceeded.");
                neighbors.push_back(j);
                planes.push_back(lo); planes.push_back(hi);
            }
            std::sort(planes.begin(),planes.end());
            planes.erase(std::unique(planes.begin(),planes.end()),planes.end());
            for (std::size_t z=1;z<planes.size();++z) {
                auto boundary=exact_boundaries[owner];
                for (const auto j : neighbors) {
                    charge(1);
                    if (std::llround(pieces[j].slab.bottom_m*1e9)>planes[z-1] ||
                        std::llround(pieces[j].slab.top_m*1e9)<planes[z]) continue;
                    const auto old_points=2*boundary.segment_count();
                    if (old_points>limits.max_points-retained_exact_points)
                        throw ResourceLimit("Exact thermal boundary geometry budget exceeded.");
                    const auto cut=boolean_exact_section_boundaries(boundary,exact_boundaries[j],
                        ExactBoundaryOperation::Difference,ExactBoundaryOrientation::Opposed,
                        exact_budget(),[&] { charge(0); });
                    charge(cut.work);
                    if (2*cut.boundary.segment_count()>limits.max_points-retained_exact_points-old_points)
                        throw ResourceLimit("Exact thermal boundary transition budget exceeded.");
                    boundary=cut.boundary;
                    if (2*boundary.segment_count()>limits.max_points-retained_exact_points)
                        throw ResourceLimit("Exact thermal boundary geometry budget exceeded.");
                }
                // The unsplit boundary and immutable per-span outputs coexist
                // until this interval has been appended to the patch list.
                if (4*boundary.segment_count()>limits.max_points-retained_exact_points)
                    throw ResourceLimit("Exact thermal split geometry budget exceeded.");
                const auto split=split_exact_boundary_segments(boundary,exact_budget(),[&] { charge(0); });
                charge(split.work);
                for (const auto& span : split.segments) {
                    if (exact_patches.size()>=limits.max_groups) throw ResourceLimit("Exact thermal patch budget exceeded.");
                    retain_exact(2);
                    exact_patches.push_back({span.span,planes[z-1],planes[z],
                        Eigen::Vector2d(span.outward_normal[0],span.outward_normal[1])});
                }
            }
        }
        for (coord_t x = index(box.min.x(),0); x <= index(box.max.x(),0); ++x)
            for (coord_t y = index(box.min.y(),1); y <= index(box.max.y(),1); ++y) {
                charge(count+1);
                const coord_t left = start[0]+x*size[0], bottom = start[1]+y*size[1];
                const auto clipped = exact ? ExPolygons{} : intersection_ex(slab.footprint.native(),Polygons{Polygon{
                    {left,bottom},{left+size[0],bottom},{left+size[0],bottom+size[1]},{left,bottom+size[1]}}});
                if (!exact && clipped.empty()) continue;
                const auto clipped_points = exact ? count : points(clipped);
                if (clipped_points > limits.max_points) throw ResourceLimit("Section thermal clip points exceeded.");
                for (coord_t z = index(lo_z,2); z <= index(hi_z-1,2); ++z) {
                    charge(64+9*clipped_points);
                    const std::array<coord_t,3> base{left,bottom,start[2]+z*size[2]};
                    SectionIntegration::Moments moments;
                    if (exact) {
                        const auto available = std::min<std::size_t>(4000000,work_limit-work);
                        if (!available) throw ResourceLimit("Exact section thermal moment budget exhausted.");
                        const auto integrated = exact_section_volume_moments(*exact_geometry,base,size,
                            lo_z,hi_z,available,[&] { charge(0); });
                        charge(integrated.work);
                        moments = integrated.values;
                        if (integrated.bounds[0][1] == 0.) continue;
                    } else {
                        moments = SectionIntegration::volume_moments(clipped,base,size,
                            std::max(lo_z,base[2]),std::min(hi_z,base[2]+size[2]));
                    }
                    if (!std::isfinite(moments[0]) || moments[0] <= 0.)
                        throw std::invalid_argument("Invalid section thermal occupied volume.");
                    const Index cell{x,y,z,coord_t(piece.component)};
                    if (!cells.count(cell) && cells.size() >= cell_limit)
                        throw ResourceLimit("Section thermal cell budget exceeded.");
                    cells.insert(cell);
                    const auto key = std::make_pair(cell,slab.material);
                    if (!groups.count(key) && groups.size() >= group_limit)
                        throw ResourceLimit("Section thermal material-group budget exceeded.");
                    auto& group = groups[key];
                    for (std::size_t i = 0; i < group.size(); ++i) group[i] += moments[i];
                    const auto add_face = [&](const ExPolygons& exposed_face,coord_t plane,SectionHorizontalFaceKind kind) {
                        if (exposed_face.empty()) return;
                        charge(points(exposed_face)+clipped_points);
                        const auto face = intersection_ex(exposed_face,clipped);
                        if (face.empty()) return;
                        const auto n = points(face);
                        if (n > limits.max_points) throw ResourceLimit("Section thermal face clip points exceeded.");
                        charge(64+9*n);
                        const auto face_key = std::make_pair(key,kind);
                        if (!faces.count(face_key) && faces.size() >= limits.max_groups)
                            throw ResourceLimit("Section thermal face-group budget exceeded.");
                        auto& target = faces[face_key];
                        const auto area = SectionIntegration::horizontal_surface_moments(face,base,size,plane);
                        for (std::size_t i = 0; i < target.size(); ++i) target[i] += area[i];
                    };
                    if (lo_z >= base[2]) add_face(lower_face,lo_z,lo_z == start[2] ?
                        SectionHorizontalFaceKind::Bed : SectionHorizontalFaceKind::Downward);
                    if (hi_z <= base[2]+size[2]) add_face(upper_face,hi_z,SectionHorizontalFaceKind::Upward);
                    const auto add_exact_face = [&](const ExactSection& face,coord_t plane,SectionHorizontalFaceKind kind) {
                        if (face.empty()) return;
                        const auto area=exact_section_horizontal_surface_moments(face,base,size,plane,exact_budget(),[&] { charge(0); });
                        charge(area.work);
                        if (area.bounds[0][1]==0.) return;
                        if (!(area.values[0]>0.)) throw std::invalid_argument("Unresolved exact thermal face area.");
                        const auto face_key=std::make_pair(key,kind);
                        if (!faces.count(face_key) && faces.size()>=limits.max_groups)
                            throw ResourceLimit("Section thermal face-group budget exceeded.");
                        auto& target=faces[face_key];
                        for (std::size_t i=0;i<target.size();++i) target[i]+=area.values[i];
                    };
                    if (exact_exposure) {
                        if (lo_z>=base[2]) add_exact_face(exact_lower,lo_z,lo_z==start[2] ? SectionHorizontalFaceKind::Bed : SectionHorizontalFaceKind::Downward);
                        if (hi_z<=base[2]+size[2]) add_exact_face(exact_upper,hi_z,SectionHorizontalFaceKind::Upward);
                        for (const auto& patch : exact_patches) {
                            const auto area=exact_section_vertical_surface_moments(patch.span,base,size,patch.bottom,patch.top,
                                ExactBoundaryCellOwnership::MaterialSide,exact_budget(),[&] { charge(0); });
                            charge(area.work);
                            if (area.bounds[0][1]==0.) continue;
                            if (!(area.values[0]>0.)) throw std::invalid_argument("Unresolved exact thermal side area.");
                            if (sides.size()>=limits.max_groups) throw ResourceLimit("Section thermal side-face budget exceeded.");
                            sides.push_back({key,patch.normal,area.values});
                        }
                    }
                    for (const auto& patch : vertical_patches) {
                        charge(1);
                        const auto area = vertical_patch_cell_moments(patch,base,size);
                        if (!area) continue;
                        charge(64);
                        if (sides.size() >= limits.max_groups) throw ResourceLimit("Section thermal side-face budget exceeded.");
                        sides.push_back({key,patch.outward_normal,*area});
                    }
                }
            }
        retained_exact_points=saved_retained_points;
    };
    if (produce) {
        const double angular_step = normal_bins ? 2.*std::acos(-1.)/normal_bins : 0.;
        std::vector<std::array<double,2>> normal_directions;
        normal_directions.reserve(normal_bins);
        for (unsigned direction=0;direction<normal_bins;++direction) {
            charge(1);
            normal_directions.push_back({std::cos(direction*angular_step),std::sin(direction*angular_step)});
        }
        // Extrusion boundaries in one band share their Z integral. Sum the
        // nine XY moments first, then form the 27 surface moments once per
        // cell/owner/normal/Z interval. This is linear integration, not a
        // geometric merge; voids and material/component labels are untouched.
        using BandSideKey = std::pair<SideKey,std::array<coord_t,2>>;
        std::unordered_map<BandSideKey,SectionIntegration::LineMoments,boost::hash<BandSideKey>> band_sides;
        if (normal_bins) band_sides.reserve(std::min<std::size_t>(4096,group_limit));
        const auto add_side = [&](const SideKey& key,const SectionIntegration::Moments& area) {
            auto it = streamed_sides.find(key);
            if (it == streamed_sides.end()) {
                if (streamed_sides.size() == group_limit)
                    throw ResourceLimit("Streamed thermal side-face budget exceeded.");
                it = streamed_sides.try_emplace(key).first;
            }
            for (std::size_t i=0;i<area.size();++i) it->second[i]+=area[i];
        };
        const auto flush_sides = [&] {
            stage = WorkStage::Vertical;
            for (const auto& [key,xy] : band_sides) {
                charge(27);
                ++expanded_side_batches;
                const coord_t base_z = start[2]+key.first.first.first[2]*size[2];
                add_side(key.first,SectionIntegration::extrude_line_moments(xy,double(size[2])*1e-9,
                    double(key.second[0]-base_z)/double(size[2]),double(key.second[1]-base_z)/double(size[2])));
            }
            band_sides.clear();
        };
        const auto horizontal = [&](std::size_t material,std::size_t component,const ExPolygons& polygons,coord_t plane,bool upward) {
            stage = WorkStage::Horizontal;
            const auto box=get_extents(polygons);
            const coord_t z=index(upward ? plane-1 : plane,2);
            const auto kind=upward ? SectionHorizontalFaceKind::Upward :
                (plane==start[2] ? SectionHorizontalFaceKind::Bed : SectionHorizontalFaceKind::Downward);
            for (coord_t x=index(box.min.x(),0);x<=index(box.max.x(),0);++x)
                for (coord_t y=index(box.min.y(),1);y<=index(box.max.y(),1);++y) {
                    charge(points(polygons)+1);
                    const std::array<coord_t,3> base{start[0]+x*size[0],start[1]+y*size[1],start[2]+z*size[2]};
                    const auto clipped=intersection_ex(polygons,Polygons{Polygon{{base[0],base[1]},
                        {base[0]+size[0],base[1]},{base[0]+size[0],base[1]+size[1]},{base[0],base[1]+size[1]}}});
                    if (clipped.empty()) continue;
                    if (points(clipped)>limits.max_points) throw ResourceLimit("Streamed thermal face points exceeded.");
                    charge(64+9*points(clipped));
                    const auto key=std::make_pair(GroupKey{Index{x,y,z,coord_t(component)},material},kind);
                    if (!faces.count(key) && faces.size()>=group_limit) throw ResourceLimit("Streamed thermal face budget exceeded.");
                    auto& target=faces[key];
                    const auto moments=SectionIntegration::horizontal_surface_moments(clipped,base,size,plane);
                    for (std::size_t i=0;i<target.size();++i) target[i]+=moments[i];
                }
        };
        const auto vertical = [&](std::size_t material,std::size_t component,const VerticalSectionPatch& patch) {
            stage = WorkStage::Vertical;
            std::array<double,2> normal{patch.outward_normal.x(),patch.outward_normal.y()};
            if (normal_bins) {
                const auto direction = (std::llround(std::atan2(normal[1],normal[0])/angular_step)+normal_bins)%normal_bins;
                normal = normal_directions[direction];
            }
            const auto first_cell = [&](coord_t coordinate,unsigned axis) {
                const coord_t cell = index(coordinate,axis);
                // Only a boundary exactly on a grid plane can belong to the
                // preceding cell; material-side clipping decides which one.
                return cell-(coordinate == start[axis]+cell*size[axis] ? 1 : 0);
            };
            for (coord_t x=first_cell(std::min(patch.a.x(),patch.b.x()),0);x<=index(std::max(patch.a.x(),patch.b.x()),0);++x)
                for (coord_t y=first_cell(std::min(patch.a.y(),patch.b.y()),1);y<=index(std::max(patch.a.y(),patch.b.y()),1);++y) {
                    charge(1);
                    const auto xy = vertical_patch_cell_xy_moments(patch,{start[0]+x*size[0],start[1]+y*size[1],0},size);
                    if (!xy) continue;
                    charge(37);
                    for (coord_t z=index(patch.bottom,2);z<=index(patch.top-1,2);++z) {
                        const std::array<coord_t,3> base{start[0]+x*size[0],start[1]+y*size[1],start[2]+z*size[2]};
                        const auto bottom = std::max(patch.bottom,base[2]), top = std::min(patch.top,base[2]+size[2]);
                        if (top <= bottom) continue;
                        const SideKey key{{Index{x,y,z,coord_t(component)},material},normal};
                        if (normal_bins) {
                            charge(9);
                            ++batched_side_samples;
                            const BandSideKey band_key{key,{bottom,top}};
                            auto found=band_sides.find(band_key);
                            if (found==band_sides.end()) {
                                if (band_sides.size()==std::min<std::size_t>(4096,group_limit)) flush_sides();
                                found=band_sides.try_emplace(band_key).first;
                            }
                            for (std::size_t i=0;i<xy->size();++i) found->second[i]+=(*xy)[i];
                        } else {
                            charge(27);
                            add_side(key,SectionIntegration::extrude_line_moments(*xy,double(size[2])*1e-9,
                                double(bottom-base[2])/double(size[2]),double(top-base[2])/double(size[2])));
                        }
                    }
                }
        };
        const SectionThermalSkin::ExactHorizontal exact_horizontal = [&](std::size_t material,std::size_t component,
            const ExactSection& face,coord_t plane,bool upward) {
            stage=WorkStage::Horizontal;
            // Cell clipping and extracted rational rings coexist with the
            // borrowed face. Reserve the backend's bounded arrangement scratch.
            if (material>=materials.size() || limits.max_points<8192 ||
                face.vertex_count()>limits.max_points-8192)
                throw ResourceLimit("Exact streamed thermal face limits exceeded.");
            if (face.empty()) return;
            const auto box=*face.bounds();
            const coord_t z=index(upward ? plane-1 : plane,2);
            const auto kind=upward ? SectionHorizontalFaceKind::Upward :
                (plane==start[2] ? SectionHorizontalFaceKind::Bed : SectionHorizontalFaceKind::Downward);
            for (coord_t x=index(box.min.x(),0);x<=index(box.max.x(),0);++x)
                for (coord_t y=index(box.min.y(),1);y<=index(box.max.y(),1);++y) {
                    charge(face.vertex_count()+1);
                    const std::array<coord_t,3> base{start[0]+x*size[0],start[1]+y*size[1],start[2]+z*size[2]};
                    const auto area=exact_section_horizontal_surface_moments(face,base,size,plane,exact_budget(),[&] { charge(0); });
                    charge(area.work);
                    if (area.bounds[0][1]==0.) continue;
                    if (!(area.values[0]>0.) || !std::isfinite(area.values[0]))
                        throw std::invalid_argument("Unresolved exact streamed horizontal area.");
                    const auto key=std::make_pair(GroupKey{Index{x,y,z,coord_t(component)},material},kind);
                    if (!faces.count(key) && faces.size()>=group_limit)
                        throw ResourceLimit("Streamed thermal face budget exceeded.");
                    auto& target=faces[key];
                    for (std::size_t i=0;i<target.size();++i) target[i]+=area.values[i];
                }
        };
        const SectionThermalSkin::ExactVertical exact_vertical = [&](std::size_t material,std::size_t component,
            const ExactSectionBoundary& span,coord_t bottom,coord_t top,const std::array<double,2>& outward) {
            stage=WorkStage::Vertical;
            if (material>=materials.size() || limits.max_points<8192 ||
                2*span.segment_count()>limits.max_points-8192)
                throw ResourceLimit("Exact streamed thermal side limits exceeded.");
            if (top<=bottom || !std::isfinite(outward[0]) || !std::isfinite(outward[1]))
                throw std::invalid_argument("Invalid exact streamed thermal side.");
            if (span.empty()) return;
            auto normal=outward;
            if (normal_bins) {
                const auto direction=(std::llround(std::atan2(normal[1],normal[0])/angular_step)+normal_bins)%normal_bins;
                normal=normal_directions[direction];
            }
            const auto box=*span.bounds();
            const auto first_cell = [&](coord_t coordinate,unsigned axis) {
                const auto cell=index(coordinate,axis);
                return cell-(coordinate==start[axis]+cell*size[axis] ? 1 : 0);
            };
            for (coord_t x=first_cell(box.min.x(),0);x<=index(box.max.x(),0);++x)
                for (coord_t y=first_cell(box.min.y(),1);y<=index(box.max.y(),1);++y)
                    for (coord_t z=index(bottom,2);z<=index(top-1,2);++z) {
                        charge(1);
                        const std::array<coord_t,3> base{start[0]+x*size[0],start[1]+y*size[1],start[2]+z*size[2]};
                        const auto area=exact_section_vertical_surface_moments(span,base,size,bottom,top,
                            ExactBoundaryCellOwnership::MaterialSide,exact_budget(),[&] { charge(0); });
                        charge(area.work);
                        if (area.bounds[0][1]==0.) continue;
                        if (!(area.values[0]>0.) || !std::isfinite(area.values[0]))
                            throw std::invalid_argument("Unresolved exact streamed vertical area.");
                        add_side({{Index{x,y,z,coord_t(component)},material},normal},area.values);
                    }
        };
        const auto mapping = produce(accept,horizontal,vertical,exact_horizontal,exact_vertical);
        flush_sides();
        stage = WorkStage::Assembly;
        if (mapping.empty() || mapping.size() > 1000000 || (groups.empty() && !allow_empty_increment))
            throw std::invalid_argument("Incomplete thermal topology mapping.");
        if (groups.empty()) return {};
        decltype(groups) final_groups;
        cells.clear();
        for (const auto& [key,moments] : groups) {
            charge(27);
            if (key.first[3] < 0 || std::size_t(key.first[3]) >= mapping.size() || mapping[key.first[3]] >= 10000)
                throw std::invalid_argument("Invalid completed thermal component mapping.");
            auto mapped = key;
            mapped.first[3] = coord_t(mapping[key.first[3]]);
            if (!final_groups.count(mapped) && final_groups.size() >= limits.max_groups)
                throw ResourceLimit("Section thermal final material-group budget exceeded.");
            if (!cells.count(mapped.first) && cells.size() >= limits.max_cells)
                throw ResourceLimit("Section thermal final cell budget exceeded.");
            cells.insert(mapped.first);
            auto& final = final_groups[mapped];
            for (std::size_t i = 0; i < final.size(); ++i) final[i] += moments[i];
        }
        groups.swap(final_groups);
        const auto remap = [&](GroupKey key) {
            charge(1);
            if (key.first[3]<0 || std::size_t(key.first[3])>=mapping.size())
                throw std::invalid_argument("Invalid thermal skin component mapping.");
            key.first[3]=coord_t(mapping[key.first[3]]);
            if (!groups.count(key)) throw std::invalid_argument("Thermal skin does not belong to occupied material.");
            return key;
        };
        decltype(faces) final_faces;
        for (const auto& [key,moments] : faces) {
            charge(27);
            const auto mapped=std::make_pair(remap(key.first),key.second);
            if (!final_faces.count(mapped) && final_faces.size()>=limits.max_groups)
                throw ResourceLimit("Streamed thermal final face budget exceeded.");
            auto& target=final_faces[mapped];
            for (std::size_t i=0;i<target.size();++i) target[i]+=moments[i];
        }
        faces.swap(final_faces);
        decltype(streamed_sides) final_sides;
        for (const auto& [key,moments] : streamed_sides) {
            charge(27);
            const SideKey mapped{remap(key.first),key.second};
            if (!final_sides.count(mapped) && final_sides.size()>=limits.max_groups)
                throw ResourceLimit("Streamed thermal final side-face budget exceeded.");
            auto& target=final_sides[mapped];
            for (std::size_t i=0;i<target.size();++i) target[i]+=moments[i];
        }
        streamed_sides.clear();
        for (const auto& [key,moments] : final_sides)
            sides.push_back({key.first,{key.second[0],key.second[1]},moments});
    } else {
        for (const auto& piece : pieces) accept(piece);
    }
    SectionThermalAssembly result;
    result.quadrature=limits.quadrature;
    stage = WorkStage::Assembly;
    result.vertical_normal_bins = normal_bins;
    result.exposure_complete = topology != nullptr || (streaming && streaming->include_exposure);
    result.occupied_cells = cells.size();
    std::map<Index,std::size_t> nodes;
    struct Mass { double volume = 0., capacity = 0., heat = 0.; };
    std::vector<Mass> mass;
    for (const auto& [key,moments] : groups) {
        charge(64);
        const auto& material = materials[key.second];
        const auto integrated = SectionIntegration::thermal_element(moments,native_spacing,
            material.volumetric_capacity_j_m3_k,material.conductivity_w_m_k,limits.quadrature);
        ThermalHistory::Element element;
        element.conductance = integrated.conduction;
        for (unsigned corner = 0; corner < 8; ++corner) {
            auto node = key.first;
            Eigen::Vector3d position;
            for (unsigned axis = 0; axis < 3; ++axis) {
                if (corner & (1u << axis)) ++node[axis];
                position[axis] = double(start[axis]+node[axis]*size[axis])*1e-9;
            }
            if (!nodes.count(node) && nodes.size() >= limits.max_nodes)
                throw ResourceLimit("Section thermal node budget exceeded.");
            const auto [entry,inserted] = nodes.emplace(node,nodes.size());
            const std::size_t id = entry->second;
            if (inserted) { result.nodes.push_back({position,std::size_t(node[3])}); mass.emplace_back(); }
            element.nodes[corner] = id;
            const double capacity = integrated.capacity.row(corner).sum();
            if (!std::isfinite(capacity) || capacity <= 0.)
                throw std::invalid_argument("Unresolved section thermal nodal mass.");
            mass[id].volume += capacity/material.volumetric_capacity_j_m3_k;
            mass[id].capacity += capacity;
            mass[id].heat += capacity*material.initial_temperature_k;
        }
        result.occupied_volume_m3 += moments[0];
        const Eigen::Vector3d minimum = result.nodes[element.nodes[0]].position_m;
        result.groups.push_back({element.nodes,key.second,std::size_t(key.first[3]),
            {int(key.first[0]),int(key.first[1]),int(key.first[2])},minimum,minimum+native_spacing,moments[0]});
        result.initial.added_elements.push_back(std::move(element));
    }
    for (const auto& [key,area] : faces) {
        charge(8);
        SectionHorizontalThermalFace face;
        face.material = key.first.second;
        face.kind = key.second;
        face.area_moments = area;
        for (unsigned corner = 0; corner < 8; ++corner) {
            auto node = key.first.first;
            for (unsigned axis = 0; axis < 3; ++axis) if (corner & (1u << axis)) ++node[axis];
            face.nodes[corner] = nodes.at(node);
        }
        result.horizontal_faces.push_back(std::move(face));
    }
    for (const auto& side : sides) {
        charge(8);
        SectionVerticalThermalFace face;
        face.material = side.group.second;
        face.outward_normal = {side.normal.x(),side.normal.y(),0.};
        face.area_moments = side.area;
        for (unsigned corner = 0; corner < 8; ++corner) {
            auto node = side.group.first;
            for (unsigned axis = 0; axis < 3; ++axis) if (corner & (1u << axis)) ++node[axis];
            face.nodes[corner] = nodes.at(node);
        }
        result.vertical_faces.push_back(std::move(face));
    }
    for (const auto& node : mass) {
        const double density = node.capacity/node.volume, temperature = node.heat/node.capacity;
        if (!std::isfinite(node.volume) || node.volume <= 0. || !std::isfinite(density) || density <= 0. ||
            !std::isfinite(temperature) || temperature <= 0.)
            throw std::invalid_argument("Section thermal nodal mass overflow.");
        result.initial.deposited_cells.push_back({node.volume,density,temperature,0.});
    }
    return result;
}

SectionThermalAssembly assemble_section_thermal_volume(const SectionTopology& topology,
    const std::vector<SectionThermalMaterial>& materials, const Eigen::Vector3d& origin,
    const Eigen::Vector3d& spacing, const SectionThermalAssemblyLimits& limits, const std::function<void()>& checkpoint)
{
    return assemble_volume(topology.pieces(),topology.bed_z_m(),&topology,materials,origin,spacing,limits,checkpoint);
}

SectionThermalAssembly assemble_section_thermal_addition(const SectionGrowth& growth,
    const std::vector<SectionThermalMaterial>& materials, const Eigen::Vector3d& origin,
    const Eigen::Vector3d& spacing, const SectionThermalAssemblyLimits& limits, const std::function<void()>& checkpoint)
{
    if (checkpoint) checkpoint();
    if (growth.added_pieces().empty()) return {};
    return assemble_volume(growth.added_pieces(),growth.bed_z_m(),nullptr,materials,origin,spacing,limits,checkpoint);
}

ThermalSectionStreamResult assemble_local_deposition_thermal_volume(const std::vector<DepositedSegment>& segments,
    const std::vector<SectionThermalMaterial>& materials,const ThermalSectionStreamSettings& settings,
    const std::function<std::optional<SectionOwnership>(double)>& ownership,
    const std::function<void()>& checkpoint,const std::vector<DepositedSegment>* partition_envelopes)
{
    ThermalSectionStreamResult result;
    result.assembly = assemble_volume({},settings.bed_z_m,nullptr,materials,settings.origin_m,settings.spacing_m,
        settings.assembly,checkpoint,[&](const auto& accept,const auto& horizontal,const auto& vertical,const auto& exact_horizontal,const auto& exact_vertical) {
            std::unique_ptr<SectionThermalSkin> skin;
            std::vector<SectionPiece> band;
            std::size_t band_points=0;
            LocalTopologyBuilder builder(settings.bed_z_m,settings.topology,[&](const SectionPiece& p) {
                accept(p);
                if (skin) {
                    const auto n=p.slab.footprint.vertex_count();
                    if (n>settings.skin.max_live_points-band_points) throw ResourceLimit("Thermal skin band points exceeded.");
                    band_points+=n; band.push_back(p);
                }
            },checkpoint,
                settings.streaming_topology,&result.topology,LocalTileOrder::Lexicographic,
                LocalAttachmentPolicy::AllowFreeComponents);
            if (settings.include_exposure)
                skin=std::make_unique<SectionThermalSkin>(horizontal,vertical,settings.skin,checkpoint,
                    [&] { return builder.provisional_component_mapping(); },
                    settings.geometry.geometry.retain_exact_interface_sections ? exact_horizontal : SectionThermalSkin::ExactHorizontal{},
                    settings.geometry.geometry.retain_exact_interface_sections ? exact_vertical : SectionThermalSkin::ExactVertical{});
            const auto source = builder.source();
            LocalSectionPacketVisitor visitor;
            visitor.begin_tile = [&](const LocalSectionTile& header) { builder.begin_tile(header); if (skin) skin->begin_tile(header); };
            visitor.band = [&](const LocalSectionTile& packet) {
                builder.add_band(packet);
                if (skin && !band.empty()) { skin->add_band(std::move(band)); band.clear(); band_points=0; }
            };
            visitor.end_tile = [&] { builder.end_tile(); if (skin) skin->end_tile(); };
            result.geometry = visit_local_deposition_packets(segments,settings.origin_m.head<2>(),settings.tile_size_m,
                settings.bed_z_m,settings.maximum_height_m,visitor,settings.geometry,ownership,checkpoint,
                settings.interval_subdivisions,nullptr,partition_envelopes);
            const auto certificate = builder.finish();
            if (skin) skin->finish();
            if (certificate.source() != source || certificate.components() > 10000)
                throw std::invalid_argument("Thermal topology source identity mismatch.");
            return certificate.component_mapping();
        },&settings);
    return result;
}

ReferenceThermalAssembly assemble_reference_thermal_volume(
    const std::vector<DepositedSegment>& reference,const std::vector<DepositedSegment>& deposited,
    const std::vector<SectionThermalMaterial>& materials,const ThermalSectionStreamSettings& settings,
    const SectionGrowthStreamLimits& activation_limits,
    const std::function<std::optional<SectionOwnership>(double)>& reference_ownership,
    const std::function<void()>& checkpoint)
{
    auto full_settings = settings;
    full_settings.include_exposure = true;
    SectionGrowthStreamSettings activation;
    static_cast<SectionGrowthStreamLimits&>(activation) = activation_limits;
    activation.origin_m = settings.origin_m.head<2>();
    activation.tile_size_m = settings.tile_size_m;
    activation.bed_z_m = settings.bed_z_m;
    activation.maximum_height_m = settings.maximum_height_m;
    activation.interval_subdivisions = settings.interval_subdivisions;
    activation.geometry = settings.geometry;
    activation.topology = settings.topology;
    activation.connectivity = settings.streaming_topology;
    ReferenceThermalAssembly result;
    result.assembly = assemble_volume({},settings.bed_z_m,nullptr,materials,settings.origin_m,settings.spacing_m,
        settings.assembly,checkpoint,[&](const auto& accept,const auto& horizontal,const auto& vertical,const auto& exact_horizontal,const auto& exact_vertical) {
            std::unique_ptr<SectionThermalSkin> skin;
            std::vector<SectionPiece> band;
            std::size_t band_points = 0;
            const auto piece = [&](const SectionPiece& p) {
                accept(p);
                const auto n = p.slab.footprint.vertex_count();
                if (n > settings.skin.max_live_points-band_points)
                    throw ResourceLimit("Reference thermal skin band budget exceeded.");
                band_points += n;
                band.push_back(p);
            };
            LocalSectionPacketVisitor packets;
            packets.begin_tile = [&](const LocalSectionTile& h) { skin->begin_tile(h); };
            packets.band = [&](const LocalSectionTile&) {
                if (!band.empty()) { skin->add_band(std::move(band)); band.clear(); band_points=0; }
            };
            packets.end_tile = [&] { skin->end_tile(); };
            ReferenceActivationCursor cursor(reference,deposited,activation,piece,reference_ownership,checkpoint,packets);
            skin=std::make_unique<SectionThermalSkin>(horizontal,vertical,settings.skin,checkpoint,
                [&] { return cursor.active_component_mapping(); },
                settings.geometry.geometry.retain_exact_interface_sections ? exact_horizontal : SectionThermalSkin::ExactHorizontal{},
                settings.geometry.geometry.retain_exact_interface_sections ? exact_vertical : SectionThermalSkin::ExactVertical{});
            while (cursor.advance()) {}
            const auto proof = cursor.finish();
            skin->finish();
            result.excluded_nominal_volume_m3 = proof.excluded_nominal_volume_m3;
            return proof.active.component_mapping();
        },&full_settings);
    return result;
}

static SectionThermalGrowthAssembly assemble_thermal_growth_impl(
    const std::vector<DepositedSegment>& previous,const std::vector<DepositedSegment>& next,
    const std::vector<SectionThermalMaterial>& materials,const ThermalSectionStreamSettings& settings,
    const SectionGrowthStreamLimits& growth_limits,
    const std::function<std::optional<SectionOwnership>(double)>& previous_ownership,
    const std::function<std::optional<SectionOwnership>(double)>& next_ownership,
    const std::function<void()>& checkpoint,const std::vector<DepositedSegment>* previous_partitions,
    const std::vector<DepositedSegment>* next_partitions,const std::vector<DepositedSegment>* reference = nullptr)
{
    auto full_settings=settings, added_settings=settings;
    full_settings.include_exposure=true;
    added_settings.include_exposure=false;
    SectionGrowthStreamSettings growth;
    static_cast<SectionGrowthStreamLimits&>(growth)=growth_limits;
    growth.origin_m=settings.origin_m.head<2>(); growth.tile_size_m=settings.tile_size_m;
    growth.bed_z_m=settings.bed_z_m; growth.maximum_height_m=settings.maximum_height_m;
    growth.interval_subdivisions=settings.interval_subdivisions;
    growth.geometry=settings.geometry; growth.topology=settings.topology; growth.connectivity=settings.streaming_topology;
    SectionThermalGrowthAssembly result;
    result.next=assemble_volume({},settings.bed_z_m,nullptr,materials,settings.origin_m,settings.spacing_m,settings.assembly,
        checkpoint,[&](const auto& accept_next,const auto& horizontal,const auto& vertical,const auto& exact_horizontal,const auto& exact_vertical) {
            SectionThermalSkin skin(horizontal,vertical,settings.skin,checkpoint,{},
                settings.geometry.geometry.retain_exact_interface_sections ? exact_horizontal : SectionThermalSkin::ExactHorizontal{},
                settings.geometry.geometry.retain_exact_interface_sections ? exact_vertical : SectionThermalSkin::ExactVertical{});
            std::vector<SectionPiece> band;
            std::size_t band_points=0;
            std::vector<std::size_t> next_mapping;
            result.added=assemble_volume({},settings.bed_z_m,nullptr,materials,settings.origin_m,settings.spacing_m,settings.assembly,
                checkpoint,[&](const auto& accept_added,const auto&,const auto&,const auto&,const auto&) {
                    SectionGrowthStreamVisitor visitor;
                    visitor.previous=[](const SectionPiece&) {};
                    visitor.next=[&](const SectionPiece& piece) {
                        accept_next(piece);
                        const auto n=piece.slab.footprint.vertex_count();
                        if (n>settings.skin.max_live_points-band_points) throw ResourceLimit("Thermal growth skin band exceeded.");
                        band_points+=n; band.push_back(piece);
                    };
                    visitor.added=accept_added;
                    visitor.next_packets.begin_tile=[&](const LocalSectionTile& h) { skin.begin_tile(h); };
                    visitor.next_packets.band=[&](const LocalSectionTile&) {
                        if (!band.empty()) { skin.add_band(std::move(band)); band.clear(); band_points=0; }
                    };
                    visitor.next_packets.end_tile=[&] { skin.end_tile(); };
                    if (reference) {
                        auto proof=visit_reference_section_growth(*reference,previous,next,growth,visitor,previous_ownership,checkpoint);
                        result.previous_components=std::move(proof.previous_components);
                        result.previous_excluded_nominal_volume_m3=proof.previous_excluded_nominal_volume_m3;
                        result.next_excluded_nominal_volume_m3=proof.next_excluded_nominal_volume_m3;
                        next_mapping=proof.next.component_mapping();
                    } else {
                        auto proof=visit_local_section_growth(previous,next,growth,visitor,previous_ownership,next_ownership,
                            checkpoint,previous_partitions,next_partitions);
                        result.previous_components=std::move(proof.previous_components);
                        next_mapping=proof.next.component_mapping();
                    }
                    skin.finish();
                    return next_mapping;
                },&added_settings,true);
            return next_mapping;
        },&full_settings);
    return result;
}
SectionThermalGrowthAssembly assemble_local_thermal_growth(
    const std::vector<DepositedSegment>& previous,const std::vector<DepositedSegment>& next,
    const std::vector<SectionThermalMaterial>& materials,const ThermalSectionStreamSettings& settings,
    const SectionGrowthStreamLimits& growth_limits,
    const std::function<std::optional<SectionOwnership>(double)>& previous_ownership,
    const std::function<std::optional<SectionOwnership>(double)>& next_ownership,
    const std::function<void()>& checkpoint,const std::vector<DepositedSegment>* previous_partitions,
    const std::vector<DepositedSegment>* next_partitions)
{
    return assemble_thermal_growth_impl(previous,next,materials,settings,growth_limits,previous_ownership,next_ownership,
        checkpoint,previous_partitions,next_partitions);
}
SectionThermalGrowthAssembly assemble_reference_thermal_growth(
    const std::vector<DepositedSegment>& reference,const std::vector<DepositedSegment>& previous,
    const std::vector<DepositedSegment>& next,const std::vector<SectionThermalMaterial>& materials,
    const ThermalSectionStreamSettings& settings,const SectionGrowthStreamLimits& growth_limits,
    const std::function<std::optional<SectionOwnership>(double)>& ownership,const std::function<void()>& checkpoint)
{
    return assemble_thermal_growth_impl(previous,next,materials,settings,growth_limits,ownership,ownership,
        checkpoint,nullptr,nullptr,&reference);
}
} // namespace Slic3r::Thermoelastic
