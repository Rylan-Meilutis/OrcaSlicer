#include "ThermoelasticField.hpp"
#include "ThermoelasticSolver.hpp"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
namespace {
constexpr double position_tolerance_m = 1e-10;
struct TreeCell {
    std::size_t index;
    Eigen::AlignedBox3d bounds;
    std::size_t idx() const { return index; }
    const Eigen::AlignedBox3d &bbox() const { return bounds; }
    Eigen::Vector3d centroid() const { return bounds.center(); }
};
} // namespace

Field::Field(std::vector<Eigen::Vector3d> nodes_m, std::vector<Element> elements, Result result)
    : m_nodes(std::move(nodes_m)), m_elements(std::move(elements)), m_result(std::move(result))
{
    if (m_nodes.size() != m_result.displacement_m.size() || m_elements.size() != m_result.stress_pa.size())
        throw std::invalid_argument("Thermoelastic field does not match its mesh.");
    for (std::size_t n = 0; n < m_nodes.size(); ++n)
        if (!m_nodes[n].allFinite() || !m_result.displacement_m[n].allFinite())
            throw std::invalid_argument("Non-finite thermoelastic field node.");
    std::vector<TreeCell> boxes;
    boxes.reserve(m_elements.size());
    m_cells.reserve(m_elements.size());
    for (std::size_t i = 0; i < m_elements.size(); ++i) {
        const auto &element = m_elements[i];
        TreeCell box{i, {}};
        for (std::size_t n : element.nodes) {
            if (n >= m_nodes.size()) throw std::invalid_argument("Invalid thermoelastic field node index.");
            box.bounds.extend(m_nodes[n]);
        }
        Eigen::Matrix3d jacobian;
        for (int j = 0; j < 3; ++j)
            jacobian.col(j) = m_nodes[element.nodes[j + 1]] - m_nodes[element.nodes[0]];
        const Eigen::FullPivLU<Eigen::Matrix3d> lu(jacobian);
        if (!lu.isInvertible()) throw std::invalid_argument("Degenerate thermoelastic field element.");
        const auto &s = m_result.stress_pa[i];
        if (!s.allFinite()) throw std::invalid_argument("Non-finite thermoelastic field stress.");
        Eigen::Matrix3d tensor;
        tensor << s[0], s[3], s[5], s[3], s[1], s[4], s[5], s[4], s[2];
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> principal(tensor, Eigen::EigenvaluesOnly);
        if (principal.info() != Eigen::Success)
            throw std::runtime_error("Thermoelastic principal stress calculation failed.");
        const double vm = std::sqrt(0.5 * (std::pow(s[0]-s[1], 2) + std::pow(s[1]-s[2], 2) +
            std::pow(s[2]-s[0], 2)) + 3. * s.tail<3>().squaredNorm());
        if (!std::isfinite(vm)) throw std::invalid_argument("Thermoelastic stress magnitude overflow.");
        m_cells.push_back({lu.inverse(), vm, principal.eigenvalues().maxCoeff()});
        box.bounds.min().array() -= position_tolerance_m;
        box.bounds.max().array() += position_tolerance_m;
        boxes.push_back(std::move(box));
    }
    m_tree.build(std::move(boxes));
}

std::optional<FieldSample> Field::sample(const Eigen::Vector3d &position_m) const
{
    std::optional<FieldSample> result;
    if (position_m.allFinite() && !m_tree.empty()) sample_node(0, position_m, result);
    return result;
}

void Field::sample_node(std::size_t tree_node, const Eigen::Vector3d &position_m,
                       std::optional<FieldSample> &result) const
{
    const auto &node = m_tree.node(tree_node);
    if (!node.is_valid() || !node.bbox.contains(position_m)) return;
    if (node.is_inner()) {
        sample_node(m_tree.left_child_idx(tree_node), position_m, result);
        sample_node(m_tree.right_child_idx(tree_node), position_m, result);
        return;
    }
    const auto &element = m_elements[node.idx];
    const auto &cell = m_cells[node.idx];
    const Eigen::Vector3d b = cell.inverse_jacobian * (position_m - m_nodes[element.nodes[0]]);
    const double b0 = 1. - b.sum();
    // Small coordinate tolerance only; never bridge a hole by sampling the
    // nearest element. On an element boundary retain the larger tensile stress
    // instead of smoothing across a material interface.
    const double tolerance = position_tolerance_m * cell.inverse_jacobian.cwiseAbs().rowwise().sum().maxCoeff();
    if (b.minCoeff() < -tolerance || b0 < -tolerance) return;
    if (result && result->maximum_principal_pa >= cell.maximum_principal_pa) return;
    Eigen::Vector3d displacement = b0 * m_result.displacement_m[element.nodes[0]];
    for (int i = 0; i < 3; ++i) displacement += b[i] * m_result.displacement_m[element.nodes[i + 1]];
    result = FieldSample{displacement, m_result.stress_pa[node.idx], cell.von_mises_pa, cell.maximum_principal_pa};
}

} // namespace Slic3r::Thermoelastic
