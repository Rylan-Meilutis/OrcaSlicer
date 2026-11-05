#include "ThermoelasticReduction.hpp"
#include "ThermoelasticSolver.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Slic3r::Thermoelastic {
LineBasis occupied_line_basis(const std::vector<Eigen::Vector3d>& nodes,
                              const std::vector<Element>& elements,
                              const std::vector<std::size_t>& pinned,
                              unsigned axis,
                              unsigned stride,
                              const ReductionLimits& limits,
                              const std::function<void()>& checkpoint)
{
    if (nodes.empty() || limits.max_fine_nodes == 0 || limits.max_fine_nodes > 100000 || elements.empty() || limits.max_elements == 0 ||
        limits.max_elements > 200000 || limits.max_coarse_nodes == 0 || limits.max_coarse_nodes > 12000 || limits.max_assembly_work == 0 ||
        limits.max_assembly_work > 20000000 || pinned.size() > nodes.size() || axis > 2 || stride == 0 || stride > 16)
        throw std::invalid_argument("Invalid occupied-line reduction inputs.");
    if (nodes.size() > limits.max_fine_nodes || elements.size() > limits.max_elements)
        throw ReductionBudget("Occupied-line reduction input budget exceeded.");
    const auto none = std::numeric_limits<std::size_t>::max();
    std::vector<std::array<std::size_t, 2>> neighbors(nodes.size(), {none, none});
    std::vector<bool> fixed(nodes.size(), false), visited(nodes.size(), false);
    for (const auto& node : nodes)
        if (!node.allFinite())
            throw std::invalid_argument("Nonfinite reduction node.");
    for (auto node : pinned) {
        if (node >= nodes.size())
            throw std::invalid_argument("Invalid pinned reduction node.");
        fixed[node] = true;
    }
    std::size_t work = 0;
    for (const auto& element : elements) {
        if (++work > limits.max_assembly_work)
            throw ReductionBudget("Reduction topology work budget exceeded.");
        if (checkpoint && work % 128 == 0)
            checkpoint();
        for (auto node : element.nodes)
            if (node >= nodes.size())
                throw std::invalid_argument("Invalid reduction element node.");
        for (unsigned a = 0; a < 4; ++a)
            for (unsigned b = a + 1; b < 4; ++b) {
                auto i = element.nodes[a], j = element.nodes[b];
                const Eigen::Vector3d delta = nodes[j] - nodes[i];
                if (delta[axis] == 0. || delta[(axis + 1) % 3] != 0. || delta[(axis + 2) % 3] != 0.)
                    continue;
                if (delta[axis] < 0.)
                    std::swap(i, j);
                if ((neighbors[i][1] != none && neighbors[i][1] != j) || (neighbors[j][0] != none && neighbors[j][0] != i))
                    throw std::invalid_argument("Ambiguous occupied-grid axis topology.");
                neighbors[i][1] = j;
                neighbors[j][0] = i;
            }
    }
    LineBasis result;
    result.rows.resize(nodes.size());
    const auto anchor = [&](std::size_t node) {
        if (result.rows[node].empty()) {
            if (result.nodes_m.size() >= limits.max_coarse_nodes || result.nodes_m.size() >= 12000)
                throw ReductionBudget("Reduction basis node budget exceeded.");
            result.rows[node] = {{result.nodes_m.size(), 1.}};
            result.nodes_m.push_back(nodes[node]);
        }
        return result.rows[node].front().node;
    };
    for (std::size_t start = 0; start < nodes.size(); ++start) {
        if (neighbors[start][0] != none)
            continue;
        if (checkpoint)
            checkpoint();
        std::vector<std::size_t> line;
        for (auto node = start; node != none; node = neighbors[node][1]) {
            if (visited[node])
                throw std::invalid_argument("Invalid occupied-grid line cycle.");
            visited[node] = true;
            line.push_back(node);
        }
        std::size_t left = 0;
        anchor(line[0]);
        for (std::size_t right = 1; right < line.size(); ++right) {
            if (right - left < stride && !fixed[line[right]] && right + 1 != line.size())
                continue;
            const auto a = anchor(line[left]), b = anchor(line[right]);
            const double length = nodes[line[right]][axis] - nodes[line[left]][axis];
            for (auto k = left + 1; k < right; ++k) {
                const double t       = (nodes[line[k]][axis] - nodes[line[left]][axis]) / length;
                result.rows[line[k]] = {{a, 1. - t}, {b, t}};
            }
            left = right;
        }
    }
    if (std::find(visited.begin(), visited.end(), false) != visited.end())
        throw std::invalid_argument("Incomplete occupied-grid line topology.");
    return result;
}
} // namespace Slic3r::Thermoelastic
