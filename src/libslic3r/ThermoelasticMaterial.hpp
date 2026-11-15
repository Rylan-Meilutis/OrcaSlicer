#pragma once
#include <Eigen/Core>

namespace Slic3r::Thermoelastic {
using ElasticityMatrix = Eigen::Matrix<double, 6, 6>;

// Engineering strain / physical stress order: xx, yy, zz, xy, yz, xz.
// Reject non-finite, nonsymmetric or non-positive-definite matrices. Only
// relative roundoff asymmetry (1e-12) is symmetrized, not invalid material data.
ElasticityMatrix validated_elasticity(const ElasticityMatrix &elasticity_pa);

// Columns of local_to_world are the local material axes in world coordinates.
// A proper orthogonal rotation is required. Eigenstrain is NOT transformed by
// this function: callers must prescribe it in the element's world frame.
ElasticityMatrix rotate_elasticity(const ElasticityMatrix &elasticity_pa,
                                  const Eigen::Matrix3d &local_to_world);
} // namespace Slic3r::Thermoelastic
