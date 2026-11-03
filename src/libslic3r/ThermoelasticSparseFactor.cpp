#include "ThermoelasticSparseFactor.hpp"
#include "ThermoelasticMesh.hpp"
#include <Eigen/SparseCore>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

namespace Slic3r::Thermoelastic {
std::size_t count_ldlt_factor_entries(const Eigen::SparseMatrix<double>& matrix,
    std::size_t max_entries, std::size_t max_work, const std::function<void()>& checkpoint)
{
    if (matrix.rows() != matrix.cols() || matrix.rows() > 192000 ||
        !max_entries || max_entries > 64000000 || !max_work || max_work > 256000000)
        throw std::invalid_argument("Invalid section factor preflight limits.");
    const int size = int(matrix.rows());
    std::vector<int> parent(size,-1), visited(size,-1);
    std::size_t entries = 0, work = 0;
    const auto charge = [&] {
        if (work == max_work) throw ResourceLimit("Section factor symbolic work budget exceeded.");
        ++work;
        if (checkpoint && (work & 1023) == 1) checkpoint();
    };
    for (int column = 0; column < size; ++column) {
        charge();
        visited[column] = column;
        for (Eigen::SparseMatrix<double>::InnerIterator entry(matrix,column); entry; ++entry) {
            charge();
            int row = int(entry.row());
            if (row >= column) continue;
            // Each newly reached ancestor is one structural L(column,row).
            // A parent is first established by its earliest later neighbor.
            while (visited[row] != column) {
                charge();
                if (entries == max_entries)
                    throw ResourceLimit("Section factor storage budget exceeded.");
                ++entries;
                visited[row] = column;
                if (parent[row] == -1) parent[row] = column;
                row = parent[row];
            }
        }
    }
    if (checkpoint) checkpoint();
    return entries;
}
}
