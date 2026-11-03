#include "PeriodicTetrahedralGeometry.hpp"
#include <Eigen/Core>
#include <Eigen/LU>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace Slic3r::Thermoelastic {
PeriodicTetrahedralGeometry validate_periodic_tetrahedral_geometry(const std::vector<Eigen::Vector3d>& nodes_m,
                                                                   const std::vector<std::array<std::size_t, 4>>& elements,
                                                                   const Eigen::Vector3d& minimum_m,
                                                                   const Eigen::Vector3d& maximum_m,
                                                                   std::size_t max_nodes,
                                                                   const std::function<void()>& checkpoint)
{
    const Eigen::Vector3d extent = maximum_m - minimum_m;
    if (!minimum_m.allFinite() || !maximum_m.allFinite() || !extent.allFinite() || extent.minCoeff() <= 0. || max_nodes == 0 ||
        max_nodes > 2048 || nodes_m.empty() || nodes_m.size() > max_nodes || elements.empty() || elements.size() > 32768)
        throw std::invalid_argument("Invalid or oversized periodic RVE.");
    const double cell_volume = extent.prod();
    if (!std::isfinite(cell_volume) || cell_volume <= 0.)
        throw std::invalid_argument("Invalid periodic RVE volume.");
    using Point = std::array<double, 3>;
    std::set<Point> distinct;
    std::map<Point, std::size_t> classes;
    std::vector<std::size_t> node_classes;
    node_classes.reserve(nodes_m.size());
    std::vector<unsigned> class_members(nodes_m.size(), 0), expected_members(nodes_m.size(), 1);
    for (std::size_t i = 0; i < nodes_m.size(); ++i) {
        if (checkpoint && i % 128 == 0)
            checkpoint();
        const auto& node = nodes_m[i];
        if (!node.allFinite() || (node.array() < minimum_m.array()).any() || (node.array() > maximum_m.array()).any())
            throw std::invalid_argument("Periodic RVE node lies outside its cell.");
        Point point{node.x(), node.y(), node.z()};
        if (!distinct.insert(point).second)
            throw std::invalid_argument("Periodic RVE contains coincident nodes.");
        unsigned expected = 1;
        for (unsigned axis = 0; axis < 3; ++axis) {
            if (node[axis] == minimum_m[axis] || node[axis] == maximum_m[axis])
                expected *= 2;
            if (node[axis] == maximum_m[axis])
                point[axis] = minimum_m[axis];
        }
        const auto entry = classes.emplace(point, classes.size()).first;
        const auto id    = entry->second;
        node_classes.push_back(id);
        ++class_members[id];
        expected_members[id] = expected;
    }
    for (std::size_t i = 0; i < classes.size(); ++i)
        if (class_members[i] != expected_members[i])
            throw std::invalid_argument("Periodic RVE has unmatched opposite boundary nodes.");

    using Face = std::array<std::size_t, 3>;
    std::map<Face, unsigned> faces;
    std::set<std::array<std::size_t, 4>> unique_elements;
    std::vector<double> volumes;
    std::vector<bool> referenced(nodes_m.size(), false);
    volumes.reserve(elements.size());
    double solid_volume = 0.;
    for (std::size_t i = 0; i < elements.size(); ++i) {
        if (checkpoint && i % 128 == 0)
            checkpoint();
        const auto& element = elements[i];
        auto sorted         = element;
        std::sort(sorted.begin(), sorted.end());
        if (sorted.back() >= nodes_m.size() || !unique_elements.insert(sorted).second)
            throw std::invalid_argument("Invalid or duplicate periodic tetrahedron.");
        for (const auto node : element)
            referenced[node] = true;
        Eigen::Matrix3d jacobian;
        for (unsigned axis = 0; axis < 3; ++axis)
            jacobian.col(axis) = nodes_m[element[axis + 1]] - nodes_m[element[0]];
        const double volume = std::abs(jacobian.determinant()) / 6.;
        if (!std::isfinite(volume) || volume <= 0.)
            throw std::invalid_argument("Degenerate periodic tetrahedron.");
        volumes.push_back(volume);
        solid_volume += volume;
        for (unsigned omit = 0; omit < 4; ++omit) {
            Face face;
            unsigned cursor = 0;
            for (unsigned k = 0; k < 4; ++k)
                if (k != omit)
                    face[cursor++] = element[k];
            std::sort(face.begin(), face.end());
            if (++faces[face] > 2)
                throw std::invalid_argument("Non-manifold periodic RVE.");
        }
    }
    if (std::find(referenced.begin(), referenced.end(), false) != referenced.end())
        throw std::invalid_argument("Periodic RVE contains unused nodes.");
    if (!std::isfinite(solid_volume) || solid_volume > cell_volume * (1. + 1e-10))
        throw std::invalid_argument("Periodic solid volume exceeds the cell volume.");
    // Matching nodes alone is insufficient: differently triangulated opposing
    // faces would not have identical piecewise-affine displacement traces.
    using FaceGeometry = std::array<Point, 3>;
    std::array<std::map<FaceGeometry, unsigned>, 3> paired_faces;
    for (const auto& [face, incidence] : faces) {
        if (checkpoint)
            checkpoint();
        if (incidence != 1)
            continue;
        for (unsigned axis = 0; axis < 3; ++axis) {
            const double value = nodes_m[face[0]][axis];
            if (value != minimum_m[axis] && value != maximum_m[axis])
                continue;
            if (nodes_m[face[1]][axis] != value || nodes_m[face[2]][axis] != value)
                continue;
            FaceGeometry geometry;
            for (unsigned k = 0; k < 3; ++k) {
                const auto& node  = nodes_m[face[k]];
                geometry[k]       = {node.x(), node.y(), node.z()};
                geometry[k][axis] = minimum_m[axis];
            }
            std::sort(geometry.begin(), geometry.end());
            const unsigned side = value == minimum_m[axis] ? 1 : 2;
            auto& mask          = paired_faces[axis][geometry];
            if (mask & side)
                throw std::invalid_argument("Duplicate periodic boundary face.");
            mask |= side;
        }
    }
    for (const auto& axis : paired_faces)
        for (const auto& [geometry, mask] : axis)
            if (mask != 3)
                throw std::invalid_argument("Periodic RVE boundary triangulations do not match.");

    return {node_classes, volumes, classes.size(), cell_volume, solid_volume};
}
} // namespace Slic3r::Thermoelastic
