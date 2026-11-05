#include "ThermoelasticThermalStrain.hpp"
#include <Eigen/Core>
#include <Eigen/LU>
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace Slic3r::Thermoelastic {
TrilinearThermalStrain trilinear_thermal_eigenstrain(const std::vector<ExpansionInterval>& intervals,
    const Eigen::Matrix3d& orientation, double reference, const std::array<double,8>& temperatures)
{
    for (double temperature : temperatures)
        if (!std::isfinite(temperature)) throw std::invalid_argument("Invalid trilinear thermal temperature.");
    const auto range = std::minmax_element(temperatures.begin(),temperatures.end());
    // Q1 interpolation is a convex combination of corner values throughout
    // the cell. Reuse the affine converter's expansion-law range validation.
    linear_thermal_eigenstrain(intervals,orientation,reference,Eigen::Vector3d::Zero(),*range.first,
        Eigen::Vector3d(*range.second-*range.first,0.,0.),Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones());
    TrilinearThermalStrain result;
    for (unsigned i = 0; i < 8; ++i)
        result.at_corners[i] = thermal_eigenstrain(intervals,orientation,reference,temperatures[i]);
    return result;
}

Eigen::Matrix<double, 6, 1> thermal_eigenstrain(const std::vector<ExpansionInterval> &intervals,
    const Eigen::Matrix3d &q, double reference, double temperature)
{
    if (intervals.empty() || intervals.size() > 64 || !std::isfinite(reference) || !std::isfinite(temperature) ||
        !q.allFinite() || (q.transpose()*q - Eigen::Matrix3d::Identity()).cwiseAbs().maxCoeff() > 1e-12 ||
        std::abs(q.determinant() - 1.) > 1e-12)
        throw std::invalid_argument("Invalid thermal expansion inputs or orientation.");
    double previous = intervals.front().minimum_k;
    for (const auto &interval : intervals) {
        if (!std::isfinite(interval.minimum_k) || !std::isfinite(interval.maximum_k) || interval.minimum_k <= 0. ||
            interval.maximum_k <= interval.minimum_k || interval.minimum_k != previous || !interval.expansion_per_k.allFinite())
            throw std::invalid_argument("Invalid thermal expansion temperature intervals.");
        const double scale = interval.expansion_per_k.cwiseAbs().maxCoeff();
        if (scale > 0. && (interval.expansion_per_k / scale - interval.expansion_per_k.transpose() / scale).cwiseAbs().maxCoeff() > 1e-12)
            throw std::invalid_argument("Thermal expansion tensor must be symmetric.");
        previous = interval.maximum_k;
    }
    const double low = std::min(reference, temperature), high = std::max(reference, temperature);
    if (low < intervals.front().minimum_k || high > intervals.back().maximum_k)
        throw std::invalid_argument("Temperature lies outside the measured thermal expansion domain.");
    Eigen::Matrix3d local = Eigen::Matrix3d::Zero();
    for (const auto &interval : intervals) {
        const double span = std::max(0., std::min(high, interval.maximum_k) - std::max(low, interval.minimum_k));
        local += span * (.5 * interval.expansion_per_k + .5 * interval.expansion_per_k.transpose());
    }
    if (temperature < reference) local = -local;
    const Eigen::Matrix3d world = q * local * q.transpose();
    Eigen::Matrix<double,6,1> result;
    result << world(0,0),world(1,1),world(2,2),2.*world(0,1),2.*world(1,2),2.*world(0,2);
    if (!result.allFinite()) throw std::overflow_error("Thermal strain overflow.");
    return result;
}

LinearThermalStrain linear_thermal_eigenstrain(const std::vector<ExpansionInterval>& intervals,
    const Eigen::Matrix3d& q, double reference, const Eigen::Vector3d& anchor, double temperature,
    const Eigen::Vector3d& gradient, const Eigen::Vector3d& minimum, const Eigen::Vector3d& maximum)
{
    if (!anchor.allFinite() || !gradient.allFinite() || !minimum.allFinite() || !maximum.allFinite() ||
        (minimum.array() > maximum.array()).any() || (anchor.array() < minimum.array()).any() ||
        (anchor.array() > maximum.array()).any())
        throw std::invalid_argument("Invalid affine thermal field domain.");
    const auto at_anchor = thermal_eigenstrain(intervals,q,reference,temperature);
    double low = temperature, high = temperature;
    for (unsigned axis = 0; axis < 3; ++axis) {
        const double a = gradient[axis]*(minimum[axis]-anchor[axis]);
        const double b = gradient[axis]*(maximum[axis]-anchor[axis]);
        low += std::min(a,b);
        high += std::max(a,b);
    }
    // Reuse the measured-domain and finite-value checks, including the
    // stress-free reference even when the entire body is currently isothermal.
    thermal_eigenstrain(intervals,q,reference,low);
    thermal_eigenstrain(intervals,q,reference,high);
    LinearThermalStrain result{at_anchor,Eigen::Matrix<double,6,3>::Zero()};
    if (low != high) {
        const Eigen::Matrix3d* expansion = nullptr;
        for (const auto& interval : intervals) {
            if (interval.maximum_k <= low || interval.minimum_k >= high) continue;
            if (expansion && (*expansion-interval.expansion_per_k).cwiseAbs().maxCoeff() != 0.)
                throw std::invalid_argument("Spatial thermal strain crosses an expansion-law change; subdivide the domain.");
            expansion = &interval.expansion_per_k;
        }
        const Eigen::Matrix3d world = q*(.5*(*expansion)+.5*expansion->transpose())*q.transpose();
        Eigen::Matrix<double,6,1> derivative;
        derivative << world(0,0),world(1,1),world(2,2),2.*world(0,1),2.*world(1,2),2.*world(0,2);
        result.gradient_per_m = derivative*gradient.transpose();
        result.at_world_zero -= result.gradient_per_m*anchor;
    }
    if (!result.at_world_zero.allFinite() || !result.gradient_per_m.allFinite())
        throw std::overflow_error("Affine thermal strain overflow.");
    return result;
}
} // namespace Slic3r::Thermoelastic
