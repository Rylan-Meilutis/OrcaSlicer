#pragma once
#include <Eigen/Core>
#include <array>
#include <vector>

namespace Slic3r::Thermoelastic {
// Piecewise-constant local expansion tensor over a stated temperature interval.
// Symmetric tensor shear entries are tensor components, not engineering shear.
// No grade, stress-free temperature, or extrapolation is inferred.
struct ExpansionInterval {
    double minimum_k;
    double maximum_k;
    Eigen::Matrix3d expansion_per_k;
};

// Integrate alpha(T) dT from explicitly supplied stress-free temperature to the
// element temperature. Intervals must be ordered, contiguous and nonoverlapping.
// local_to_world columns are material axes. The returned vector is world-frame
// engineering strain (xx,yy,zz,xy,yz,xz), suitable for Element::eigenstrain.
Eigen::Matrix<double, 6, 1> thermal_eigenstrain(
    const std::vector<ExpansionInterval> &intervals, const Eigen::Matrix3d &local_to_world,
    double stress_free_temperature_k, double temperature_k);

struct LinearThermalStrain {
    Eigen::Matrix<double, 6, 1> at_world_zero;
    Eigen::Matrix<double, 6, 3> gradient_per_m;
};
// Q1 strain on the queried cell. Corner bits encode high X/Y/Z. This keeps
// cross terms of the thermal field instead of fitting an affine gradient.
struct TrilinearThermalStrain {
    std::array<Eigen::Matrix<double,6,1>,8> at_corners;
};
// Exact for a constant expansion tensor over the spatial temperature range.
// A varying expansion law requires subdivision, not interpolation of strain
// through a constitutive breakpoint. The reference may cross breakpoints.
TrilinearThermalStrain trilinear_thermal_eigenstrain(
    const std::vector<ExpansionInterval>& intervals, const Eigen::Matrix3d& local_to_world,
    double stress_free_temperature_k, const std::array<double,8>& temperatures_k);
// Exact conversion of an explicitly supplied affine temperature field on a
// world-space box. This does not solve heat transfer or infer fan coefficients.
// The stress-free temperature may cross expansion intervals, but the spatial
// temperature range must have constant expansion: otherwise strain is not
// affine and the caller must subdivide, not replace it by a fitted gradient.
LinearThermalStrain linear_thermal_eigenstrain(
    const std::vector<ExpansionInterval>& intervals, const Eigen::Matrix3d& local_to_world,
    double stress_free_temperature_k, const Eigen::Vector3d& anchor_m, double anchor_temperature_k,
    const Eigen::Vector3d& temperature_gradient_k_per_m,
    const Eigen::Vector3d& minimum_m, const Eigen::Vector3d& maximum_m);
} // namespace Slic3r::Thermoelastic
