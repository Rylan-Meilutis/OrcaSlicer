#pragma once

#include <vector>

namespace Slic3r::Thermoelastic {

struct PronyTerm {
    double fraction; // Fraction of the instantaneous scalar modulus.
    double relaxation_seconds;
};

// Isothermal scalar generalized-Maxwell material-point building block, not a
// coupled FE material. Elastic strain excludes thermal/free strain. Parameters
// are explicit caller data; no grade selection or temperature shift is inferred.
class RelaxationHistory {
public:
    RelaxationHistory(double instantaneous_modulus_pa, double reference_temperature_c,
                      const std::vector<PronyTerm> &terms);

    // Exact update for a linear elastic-strain ramp over elapsed_seconds.
    // A zero duration denotes an instantaneous step. Other temperatures are
    // rejected: using a reference-temperature spectrum there requires shift data.
    double advance(double elastic_strain, double elapsed_seconds, double temperature_c);
    double stress_pa() const { return m_stress_pa; }

private:
    double m_modulus_pa;
    double m_temperature_c;
    double m_equilibrium_fraction;
    double m_strain = 0.;
    double m_stress_pa = 0.;
    std::vector<PronyTerm> m_terms;
    std::vector<double> m_branch_stress_pa;
};

} // namespace Slic3r::Thermoelastic
