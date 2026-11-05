#include <catch2/catch_test_macros.hpp>
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticSparseFactor.hpp"
#include <Eigen/Core>
#include <Eigen/OrderingMethods>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <cstddef>
#include <stdexcept>
#include <vector>

using namespace Slic3r::Thermoelastic;

TEST_CASE("Symbolic section factor counts match numerical factors including structural zeros", "[ThermoelasticSparseFactor]")
{
    for (int pattern = 0; pattern < 4; ++pattern) {
        DYNAMIC_SECTION("pattern " << pattern) {
            constexpr int size = 24;
            std::vector<Eigen::Triplet<double>> entries;
            for (int i = 0; i < size; ++i) {
                entries.emplace_back(i,i,100.);
                for (int j = 0; j < i; ++j)
                    if (pattern == 0 || (pattern == 1 && i-j == 1) ||
                        (pattern == 2 && ((i*7+j*11)%13 < 3)) || (pattern == 3 && j == 0)) {
                        const double value = pattern == 3 ? 0. : -.1;
                        entries.emplace_back(i,j,value);
                        entries.emplace_back(j,i,value);
                    }
            }
            Eigen::SparseMatrix<double> matrix(size,size);
            matrix.setFromTriplets(entries.begin(),entries.end());
            const auto count = count_ldlt_factor_entries(matrix,10000,100000);
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>,Eigen::Lower,Eigen::NaturalOrdering<int>> solver;
            solver.compute(matrix);
            REQUIRE(solver.info() == Eigen::Success);
            CHECK(count == std::size_t(solver.matrixL().nestedExpression().nonZeros()));
            REQUIRE(count > 1);
            CHECK(count_ldlt_factor_entries(matrix,count,100000) == count);
            CHECK_THROWS_AS(count_ldlt_factor_entries(matrix,count-1,100000),ResourceLimit);
            CHECK_THROWS_AS(count_ldlt_factor_entries(matrix,10000,1),ResourceLimit);
            CHECK_THROWS_AS(count_ldlt_factor_entries(matrix,10000,100000,[] {
                throw std::runtime_error("cancel symbolic count");
            }),std::runtime_error);
        }
    }
}

TEST_CASE("Symbolic factor preflight validates dimensions and budgets", "[ThermoelasticSparseFactor]")
{
    Eigen::SparseMatrix<double> diagonal(3,3);
    diagonal.setIdentity();
    CHECK(count_ldlt_factor_entries(diagonal,1,100) == 0);
    CHECK_THROWS_AS(count_ldlt_factor_entries(diagonal,0,100),std::invalid_argument);
    CHECK_THROWS_AS(count_ldlt_factor_entries(diagonal,64000001,100),std::invalid_argument);
    CHECK_THROWS_AS(count_ldlt_factor_entries(diagonal,1,0),std::invalid_argument);
    CHECK_THROWS_AS(count_ldlt_factor_entries(diagonal,1,256000001),std::invalid_argument);
    Eigen::SparseMatrix<double> nonsquare(2,3), oversized(192001,192001);
    CHECK_THROWS_AS(count_ldlt_factor_entries(nonsquare,1,100),std::invalid_argument);
    CHECK_THROWS_AS(count_ldlt_factor_entries(oversized,1,100),std::invalid_argument);
}

TEST_CASE("Large diagonal factor patterns require only linear symbolic work", "[ThermoelasticSparseFactor]")
{
    constexpr int size = 192000;
    Eigen::SparseMatrix<double> matrix(size,size);
    matrix.setIdentity();
    CHECK(count_ldlt_factor_entries(matrix,1,2*size) == 0);
    CHECK_THROWS_AS(count_ldlt_factor_entries(matrix,1,2*size-1),ResourceLimit);
}
