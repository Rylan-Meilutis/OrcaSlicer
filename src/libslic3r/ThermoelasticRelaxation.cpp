#include "ThermoelasticRelaxation.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace Slic3r::Thermoelastic {

RelaxationHistory::RelaxationHistory(double modulus, double temperature, const std::vector<PronyTerm> &terms)
    : m_modulus_pa(modulus), m_temperature_c(temperature), m_equilibrium_fraction(1.)
{
    if (!std::isfinite(modulus) || modulus <= 0. || !std::isfinite(temperature) || temperature <= -273.15 || terms.size() > 64)
        throw std::invalid_argument("Invalid scalar relaxation material.");
    double sum = 0.;
    for (const auto &term : terms) {
        if (!std::isfinite(term.fraction) || term.fraction < 0. ||
            !std::isfinite(term.relaxation_seconds) || term.relaxation_seconds <= 0.)
            throw std::invalid_argument("Invalid Prony spectrum.");
        sum += term.fraction;
    }
    if (sum > 1.) throw std::invalid_argument("Prony spectrum is not passive.");
    m_equilibrium_fraction = 1. - sum;
    m_terms = terms;
    m_branch_stress_pa.assign(terms.size(), 0.);
}

double RelaxationHistory::advance(double strain, double dt, double temperature)
{
    if (!std::isfinite(strain) || !std::isfinite(dt) || dt < 0. || temperature != m_temperature_c)
        throw std::invalid_argument("Relaxation update requires finite strain, nonnegative time and the reference temperature.");
    const double delta = strain - m_strain;
    if (!std::isfinite(delta)) throw std::overflow_error("Relaxation strain increment overflow.");
    std::array<double, 64> next{};
    double stress = m_modulus_pa * m_equilibrium_fraction * strain;
    for (std::size_t i = 0; i < m_terms.size(); ++i) {
        const double x = dt / m_terms[i].relaxation_seconds;
        const double decay = std::exp(-x);
        // expm1 avoids cancellation for short steps. An underflowed x is
        // indistinguishable from an instantaneous step at double precision.
        const double ramp = x == 0. ? 1. : -std::expm1(-x) / x;
        next[i] = decay * m_branch_stress_pa[i] + (m_modulus_pa * m_terms[i].fraction) * ramp * delta;
        if (!std::isfinite(next[i])) throw std::overflow_error("Relaxation branch stress overflow.");
        stress += next[i];
    }
    if (!std::isfinite(stress)) throw std::overflow_error("Relaxation stress overflow.");
    for (std::size_t i = 0; i < m_terms.size(); ++i) m_branch_stress_pa[i] = next[i];
    m_strain = strain;
    m_stress_pa = stress;
    return stress;
}

} // namespace Slic3r::Thermoelastic
