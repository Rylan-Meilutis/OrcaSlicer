#pragma once

#include "ThermoelasticSolver.hpp"
#include "AABBTreeIndirect.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <optional>
#include <vector>

namespace Slic3r::Thermoelastic {

struct FieldSample {
    Eigen::Vector3d displacement_m;
    Eigen::Matrix<double, 6, 1> stress_pa;
    double von_mises_pa;
    double maximum_principal_pa;
};

// Spatial lookup of a solved, element-wise field. No nearest-surface
// extrapolation: gaps, unsupported islands and points outside the supplied
// volume return no value, rather than a safe-looking zero stress.
class Field {
public:
    Field(std::vector<Eigen::Vector3d> nodes_m, std::vector<Element> elements, Result result);
    std::optional<FieldSample> sample(const Eigen::Vector3d &position_m) const;

private:
    struct Cell {
        Eigen::Matrix3d inverse_jacobian;
        double von_mises_pa;
        double maximum_principal_pa;
    };
    void sample_node(std::size_t tree_node, const Eigen::Vector3d &position_m,
                     std::optional<FieldSample> &result) const;
    std::vector<Eigen::Vector3d> m_nodes;
    std::vector<Element> m_elements;
    Result m_result;
    std::vector<Cell> m_cells;
    AABBTreeIndirect::Tree<3, double> m_tree;
};

} // namespace Slic3r::Thermoelastic
