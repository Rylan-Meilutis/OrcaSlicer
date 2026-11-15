#include "ThermoelasticMaterial.hpp"
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <cmath>
#include <stdexcept>

namespace Slic3r::Thermoelastic {
ElasticityMatrix validated_elasticity(const ElasticityMatrix &c)
{
    if (!c.allFinite()) throw std::invalid_argument("Non-finite elasticity tensor.");
    const double scale = c.cwiseAbs().maxCoeff();
    if (!(scale > 0.) || (c / scale - c.transpose() / scale).cwiseAbs().maxCoeff() > 1e-12)
        throw std::invalid_argument("Elasticity tensor must be symmetric positive definite.");
    const ElasticityMatrix normalized = .5 * (c / scale) + .5 * (c.transpose() / scale);
    const Eigen::LLT<ElasticityMatrix> factor(normalized);
    if (factor.info() != Eigen::Success)
        throw std::invalid_argument("Elasticity tensor must be positive definite.");
    return .5 * c + .5 * c.transpose();
}

ElasticityMatrix rotate_elasticity(const ElasticityMatrix &input, const Eigen::Matrix3d &q)
{
    const ElasticityMatrix c = validated_elasticity(input);
    if (!q.allFinite() || (q.transpose() * q - Eigen::Matrix3d::Identity()).cwiseAbs().maxCoeff() > 1e-12 ||
        std::abs(q.determinant() - 1.) > 1e-12)
        throw std::invalid_argument("Material orientation must be a proper orthogonal rotation.");
    ElasticityMatrix rotated;
    for (int column = 0; column < 6; ++column) {
        Eigen::Matrix3d strain = Eigen::Matrix3d::Zero();
        if (column < 3) strain(column, column) = 1.;
        else {
            const int first = column == 4 ? 1 : 0;
            const int second = column == 3 ? 1 : 2;
            strain(first, second) = strain(second, first) = .5;
        }
        const Eigen::Matrix3d local = q.transpose() * strain * q;
        Eigen::Matrix<double, 6, 1> e;
        e << local(0,0), local(1,1), local(2,2), 2.*local(0,1), 2.*local(1,2), 2.*local(0,2);
        const Eigen::Matrix<double, 6, 1> s = c * e;
        Eigen::Matrix3d stress;
        stress << s[0],s[3],s[5], s[3],s[1],s[4], s[5],s[4],s[2];
        const Eigen::Matrix3d world = q * stress * q.transpose();
        rotated.col(column) << world(0,0),world(1,1),world(2,2),world(0,1),world(1,2),world(0,2);
    }
    return validated_elasticity(rotated);
}
} // namespace Slic3r::Thermoelastic
