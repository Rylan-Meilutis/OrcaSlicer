#include "ThermoelasticMesh.hpp"
#include "ThermoelasticSolver.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {

std::vector<double> component_top_z(const DepositionGrid &grid, const std::function<void()> &checkpoint)
{
    if (!grid.origin_m.allFinite() || !grid.spacing_m.allFinite() || grid.spacing_m.minCoeff() <= 0.)
        throw std::invalid_argument("Invalid component-height grid.");
    std::map<std::array<int, 3>, std::size_t> ids;
    for (std::size_t i = 0; i < grid.voxels.size(); ++i) {
        for (const auto index : grid.voxels[i].index)
            if (index < 0 || index >= std::numeric_limits<int>::max() - 1)
                throw std::invalid_argument("Invalid component-height index.");
        if (!ids.emplace(grid.voxels[i].index, i).second)
            throw std::invalid_argument("Duplicate component-height voxel.");
    }
    std::vector<double> tops(grid.voxels.size(), 0.);
    std::vector<bool> seen(grid.voxels.size(), false);
    for (std::size_t start = 0; start < grid.voxels.size(); ++start) {
        if (seen[start]) continue;
        std::vector<std::size_t> component{start};
        seen[start] = true;
        int top = 0;
        for (std::size_t cursor = 0; cursor < component.size(); ++cursor) {
            if (checkpoint && cursor % 128 == 0) checkpoint();
            const auto index = grid.voxels[component[cursor]].index;
            top = std::max(top, index[2] + 1);
            for (unsigned axis = 0; axis < 3; ++axis)
                for (int sign : {-1, 1}) {
                    auto neighbor = index;
                    neighbor[axis] += sign;
                    const auto found = ids.find(neighbor);
                    if (found != ids.end() && !seen[found->second]) {
                        seen[found->second] = true;
                        component.push_back(found->second);
                    }
                }
        }
        for (const auto id : component) tops[id] = grid.origin_m.z() + top * grid.spacing_m.z();
    }
    return tops;
}

DepositionGrid voxelize_deposition(const std::vector<DepositedSegment> &segments,
    const std::vector<Voxel> &materials, const Eigen::Vector3d &spacing_m,
    double bed_z_m, std::size_t max_cells, const std::function<void()> &checkpoint,
    std::size_t max_cell_visits)
{
    if (segments.empty() || materials.empty() || !spacing_m.allFinite() ||
        spacing_m.minCoeff() <= 0. || !std::isfinite(bed_z_m) || max_cells == 0 || max_cell_visits == 0)
        throw std::invalid_argument("Invalid mechanical deposition grid.");
    Eigen::Vector3d origin = segments.front().start_m;
    for (const auto &segment : segments) {
        if (checkpoint) checkpoint();
        if (!segment.start_m.allFinite() || !segment.end_m.allFinite() ||
            !std::isfinite(segment.width_m) || !std::isfinite(segment.height_m) ||
            segment.width_m <= 0. || segment.height_m <= 0. || segment.material >= materials.size())
            throw std::invalid_argument("Invalid mechanical deposition segment.");
        if (std::abs(segment.start_m.z() - segment.end_m.z()) > 1e-9)
            throw std::invalid_argument("Nonplanar deposition requires a different mechanical discretization.");
        if (segment.end_m.z() - segment.height_m < bed_z_m - 1e-9)
            throw std::invalid_argument("Mechanical deposition extends below the bed.");
        if (spacing_m.head<2>().maxCoeff() > segment.width_m / 2. || spacing_m.z() > segment.height_m / 2.)
            throw std::invalid_argument("Mechanical grid does not resolve the deposited bead thickness.");
        origin.head<2>() = origin.head<2>().cwiseMin(
            (segment.start_m.head<2>().cwiseMin(segment.end_m.head<2>()).array() - segment.width_m / 2.).matrix());
    }
    origin.z() = bed_z_m;
    std::map<std::array<int, 3>, std::size_t> occupied;
    std::size_t visits = 0;
    const auto visit = [&] {
        if (visits >= max_cell_visits) throw ResourceLimit("Mechanical raster work budget exceeded.");
        if (++visits % 1024 == 0 && checkpoint) checkpoint();
    };
    const auto grid_index = [&](double coordinate, unsigned axis) {
        const double index = std::floor((coordinate - origin[axis]) / spacing_m[axis]);
        if (!std::isfinite(index) || index > double(std::numeric_limits<int>::max() - 2))
            throw std::invalid_argument("Mechanical deposition grid coordinates overflow.");
        return int(std::max(0., index));
    };
    for (const auto &segment : segments) {
        const Eigen::Vector2d start = segment.start_m.head<2>(), delta = segment.end_m.head<2>() - start;
        const unsigned major = std::abs(delta.x()) >= std::abs(delta.y()) ? 0 : 1, minor = 1 - major;
        const double radius = segment.width_m / 2., length_squared = delta.squaredNorm();
        const int first = grid_index(std::min(start[major], segment.end_m[major]) - radius, major);
        const int last = grid_index(std::max(start[major], segment.end_m[major]) + radius, major);
        const int bottom = grid_index(segment.end_m.z() - segment.height_m, 2);
        const int top = grid_index(segment.end_m.z(), 2);
        const double slope = delta[major] == 0. ? 0. : delta[minor] / delta[major];
        const double margin = radius * std::sqrt(1. + slope * slope);
        for (int i = first; i <= last; ++i) {
            visit();
            if (checkpoint) checkpoint();
            Eigen::Vector2d point;
            point[major] = origin[major] + (double(i) + .5) * spacing_m[major];
            const double t = delta[major] == 0. ? 0. : std::clamp((point[major] - start[major]) / delta[major], 0., 1.);
            const double center = start[minor] + t * delta[minor];
            const int lo = grid_index(center - margin, minor), hi = grid_index(center + margin, minor);
            for (int j = lo; j <= hi; ++j) {
                visit();
                point[minor] = origin[minor] + (double(j) + .5) * spacing_m[minor];
                const double projection = length_squared == 0. ? 0. : std::clamp((point - start).dot(delta) / length_squared, 0., 1.);
                if ((point - start - projection * delta).squaredNorm() > radius * radius) continue;
                for (int k = bottom; k <= top; ++k) {
                    visit();
                    const double z = origin.z() + (double(k) + .5) * spacing_m.z();
                    if (z < segment.end_m.z() - segment.height_m || z >= segment.end_m.z()) continue;
                    std::array<int, 3> index;
                    index[major] = i; index[minor] = j; index[2] = k;
                    const auto found = occupied.find(index);
                    if (found != occupied.end()) {
                        if (found->second != segment.material)
                            throw std::runtime_error("Mechanical grid cannot resolve a material interface.");
                    } else {
                        if (occupied.size() >= max_cells)
                            throw ResourceLimit("Mechanical deposition exceeds the cell budget; analysis is unresolved.");
                        occupied.emplace(index, segment.material);
                    }
                }
            }
        }
    }
    DepositionGrid grid{origin, spacing_m, {}};
    grid.voxels.reserve(occupied.size());
    for (const auto &[index, material] : occupied) {
        Voxel voxel = materials[material];
        voxel.index = index;
        voxel.material = material;
        grid.voxels.push_back(voxel);
    }
    return grid;
}

Mesh mesh_voxels(const std::vector<Voxel> &voxels, const Eigen::Vector3d &origin_m,
                 const Eigen::Vector3d &spacing_m, std::size_t max_nodes,
                 const std::function<void()> &checkpoint)
{
    if (!origin_m.allFinite() || !spacing_m.allFinite() || spacing_m.minCoeff() <= 0. || max_nodes < 4)
        throw std::invalid_argument("Invalid thermoelastic voxel grid.");
    using Index = std::array<int, 3>;
    std::map<Index, std::size_t> occupied;
    for (std::size_t i = 0; i < voxels.size(); ++i) {
        if (checkpoint && i % 128 == 0) checkpoint();
        for (int v : voxels[i].index)
            if (v < 0 || v >= std::numeric_limits<int>::max() - 1)
                throw std::invalid_argument("Invalid thermoelastic voxel index.");
        if (!occupied.emplace(voxels[i].index, i).second)
            throw std::invalid_argument("Duplicate thermoelastic voxel.");
    }
    // Merge corners through actual shared faces, not merely by position and
    // connected component. Two sides of a slit may belong to the same body
    // elsewhere without being bonded across that slit.
    if (voxels.size() > std::numeric_limits<std::size_t>::max() / 8)
        throw std::invalid_argument("Thermoelastic voxel count overflows.");
    std::vector<std::size_t> parent(8 * voxels.size());
    std::iota(parent.begin(), parent.end(), 0);
    const auto root = [&parent](std::size_t node) {
        while (parent[node] != node) {
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        return node;
    };
    for (std::size_t i = 0; i < voxels.size(); ++i) {
        if (checkpoint && i % 128 == 0) checkpoint();
        for (unsigned axis = 0; axis < 3; ++axis) {
            Index neighbor = voxels[i].index;
            ++neighbor[axis];
            const auto found = occupied.find(neighbor);
            if (found == occupied.end()) continue;
            const unsigned bit = 1u << axis;
            for (unsigned corner = 0; corner < 8; ++corner)
                if (corner & bit)
                    parent[root(8 * i + corner)] = root(8 * found->second + (corner ^ bit));
        }
    }
    Mesh mesh;
    std::map<std::size_t, std::size_t> node_ids;
    constexpr std::array<std::array<unsigned, 4>, 6> tetrahedra{{
        {{0,1,3,7}}, {{0,3,2,7}}, {{0,2,6,7}}, {{0,6,4,7}}, {{0,4,5,7}}, {{0,5,1,7}}}};
    for (std::size_t i = 0; i < voxels.size(); ++i) {
        if (checkpoint && i % 128 == 0) checkpoint();
        const auto &voxel = voxels[i];
        std::array<std::size_t, 8> nodes;
        for (unsigned corner = 0; corner < 8; ++corner) {
            Index index = voxel.index;
            for (unsigned axis = 0; axis < 3; ++axis) index[axis] += (corner >> axis) & 1;
            const auto key = root(8 * i + corner);
            auto found = node_ids.find(key);
            if (found == node_ids.end()) {
                if (mesh.nodes_m.size() >= max_nodes)
                    throw ResourceLimit("Mechanical mesh exceeds the node budget; analysis is unresolved.");
                const Eigen::Vector3d point = origin_m + spacing_m.cwiseProduct(
                    Eigen::Vector3d(double(index[0]), double(index[1]), double(index[2])));
                if (!point.allFinite()) throw std::invalid_argument("Thermoelastic voxel coordinates overflow.");
                found = node_ids.emplace(key, mesh.nodes_m.size()).first;
                mesh.nodes_m.push_back(point);
            }
            nodes[corner] = found->second;
        }
        for (const auto &tet : tetrahedra) {
            Element element{{nodes[tet[0]], nodes[tet[1]], nodes[tet[2]], nodes[tet[3]]},
                voxel.young_pa, voxel.poisson};
            element.eigenstrain = voxel.eigenstrain;
            mesh.elements.push_back(element);
        }
        if (voxel.index[2] == 0) {
            mesh.bed_faces.push_back({nodes[0], nodes[1], nodes[3]});
            mesh.bed_faces.push_back({nodes[0], nodes[3], nodes[2]});
        }
    }
    return mesh;
}

RefinementResult solve_refined_voxels(const std::vector<Voxel> &voxels,
    const Eigen::Vector3d &origin_m, const Eigen::Vector3d &spacing_m,
    std::size_t max_nodes, unsigned max_subdivisions, double relative_tolerance,
    const std::function<void(unsigned)> &progress, const std::function<void()> &checkpoint)
{
    if (voxels.empty() || max_subdivisions == 0 || max_subdivisions > 16 ||
        !std::isfinite(relative_tolerance) || relative_tolerance <= 0. || relative_tolerance >= 1.)
        throw std::invalid_argument("Invalid mechanical refinement parameters.");
    RefinementResult result;
    unsigned stable_refinements = 0;
    const auto change = [](double a, double b) {
        const double magnitude = std::max(std::abs(a), std::abs(b));
        return magnitude == 0. ? 0. : std::abs(a - b) / magnitude;
    };
    for (unsigned subdivisions = 1; subdivisions <= max_subdivisions; subdivisions *= 2) {
        if (checkpoint) checkpoint();
        if (progress) progress(subdivisions);
        const std::size_t per_voxel = std::size_t(subdivisions) * subdivisions * subdivisions;
        // Each occupied cell needs its own lower corner, even in a fully dense
        // grid. Bound allocation before expanding the grid or assembling K.
        if (voxels.size() > max_nodes / per_voxel) {
            result.resource_limited = true;
            break;
        }
        std::vector<Voxel> refined;
        refined.reserve(voxels.size() * per_voxel);
        for (std::size_t i = 0; i < voxels.size(); ++i) {
            if (checkpoint && i % 128 == 0) checkpoint();
            for (int coordinate : voxels[i].index)
                if (coordinate < 0 || coordinate > (std::numeric_limits<int>::max() - int(subdivisions)) / int(subdivisions))
                    throw std::invalid_argument("Refined mechanical coordinates overflow.");
            for (unsigned z = 0; z < subdivisions; ++z)
                for (unsigned y = 0; y < subdivisions; ++y)
                    for (unsigned x = 0; x < subdivisions; ++x) {
                        Voxel voxel = voxels[i];
                        voxel.index = {voxel.index[0] * int(subdivisions) + int(x),
                                       voxel.index[1] * int(subdivisions) + int(y),
                                       voxel.index[2] * int(subdivisions) + int(z)};
                        refined.push_back(voxel);
                    }
        }
        // Sparse/slit meshes can have eight nodes per voxel; mesh_voxels
        // checks the exact count as it constructs the mesh.
        Mesh mesh = mesh_voxels(refined, origin_m, spacing_m / double(subdivisions), max_nodes, checkpoint);
        BedResult bed = solve_bonded_bed(mesh.nodes_m, mesh.elements, mesh.bed_faces, origin_m.z(), checkpoint);
        if (result.subdivisions != 0) {
            result.force_relative_change = change(bed.tensile_force_n, result.bed.tensile_force_n);
            result.energy_relative_change = change(bed.equilibrium.elastic_energy_j, result.bed.equilibrium.elastic_energy_j);
            const bool stable = result.force_relative_change <= relative_tolerance &&
                                result.energy_relative_change <= relative_tolerance;
            stable_refinements = stable ? stable_refinements + 1 : 0;
        }
        result.mesh = std::move(mesh);
        result.bed = std::move(bed);
        result.subdivisions = subdivisions;
        result.converged = stable_refinements >= 2;
        if (result.converged) break;
    }
    return result;
}

} // namespace Slic3r::Thermoelastic
