#pragma once

#include <Eigen/SparseCore>
#include <cstddef>
#include <functional>

namespace Slic3r::Thermoelastic {
// Count the strict lower factor pattern for unpivoted LDLT in the supplied
// ordering. Only the upper triangle is read; explicit zeros remain structural.
// An elimination-tree reach needs O(rows) auxiliary storage, not factor storage.
// This is an allocation preflight, NOT a test of positive definiteness. The
// numerical factorization must use this same ordering and still check pivots.
std::size_t count_ldlt_factor_entries(const Eigen::SparseMatrix<double>& matrix,
    std::size_t max_entries, std::size_t max_work,
    const std::function<void()>& checkpoint = {});
}
