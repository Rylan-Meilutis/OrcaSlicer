#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include "libslic3r/ThermoelasticRelaxation.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using Slic3r::Thermoelastic::PronyTerm;
using Slic3r::Thermoelastic::RelaxationHistory;
using Catch::Matchers::WithinAbs;

TEST_CASE("Scalar relaxation reproduces step and hold responses", "[ThermoelasticRelaxation]")
{
    RelaxationHistory history(1000., 20., {{.6, 2.}});
    CHECK_THAT(history.advance(.01, 0., 20.), WithinAbs(10., 1e-12));
    CHECK_THAT(history.advance(.01, 2., 20.), WithinAbs(4. + 6. * std::exp(-1.), 1e-12));
    CHECK_THAT(history.advance(.01, 2000., 20.), WithinAbs(4., 1e-12));
    RelaxationHistory elastic(1000., 20., {});
    CHECK_THAT(elastic.advance(-.02, 100., 20.), WithinAbs(-20., 1e-12));
    RelaxationHistory fluid(1000., 20., {{1., 2.}});
    fluid.advance(.01, 0., 20.);
    CHECK_THAT(fluid.advance(.01, 2000., 20.), WithinAbs(0., 1e-12));
}

TEST_CASE("Linear strain ramps are independent of time partition", "[ThermoelasticRelaxation]")
{
    RelaxationHistory whole(1000., 20., {{.6, 2.}});
    RelaxationHistory split = whole;
    const double expected = 4. + 6. * (1. - std::exp(-1.));
    CHECK_THAT(whole.advance(.01, 2., 20.), WithinAbs(expected, 1e-12));
    for (int i = 1; i <= 100; ++i) split.advance(.01 * i / 100., .02, 20.);
    CHECK_THAT(split.stress_pa(), WithinAbs(expected, 1e-12));
    // The same semigroup property applies on unloading, with signed branch stress.
    whole.advance(-.01, 2., 20.);
    for (int i = 1; i <= 100; ++i) split.advance(.01 - .02 * i / 100., .02, 20.);
    CHECK_THAT(split.stress_pa(), WithinAbs(whole.stress_pa(), 1e-12));
}

TEST_CASE("Relaxation updates remain stable at very short and long times", "[ThermoelasticRelaxation]")
{
    RelaxationHistory short_step(1000., 20., {{1., 2.}});
    CHECK_THAT(short_step.advance(.01, 1e-20, 20.), WithinAbs(10., 1e-12));
    RelaxationHistory long_step(1000., 20., {{1., 1e-200}});
    CHECK_THAT(long_step.advance(.01, 1e200, 20.), WithinAbs(0., 1e-12));
}

TEST_CASE("Invalid relaxation parameters and temperature shifts are rejected", "[ThermoelasticRelaxation]")
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(RelaxationHistory(0., 20., {}), std::invalid_argument);
    CHECK_THROWS_AS(RelaxationHistory(1000., nan, {}), std::invalid_argument);
    CHECK_THROWS_AS(RelaxationHistory(1000., -273.15, {}), std::invalid_argument);
    CHECK_THROWS_AS(RelaxationHistory(1000., 20., {{-.1, 2.}}), std::invalid_argument);
    CHECK_THROWS_AS(RelaxationHistory(1000., 20., {{1.1, 2.}}), std::invalid_argument);
    CHECK_THROWS_AS(RelaxationHistory(1000., 20., {{.6, 2.}, {.6, 3.}}), std::invalid_argument);
    CHECK_THROWS_AS(RelaxationHistory(1000., 20., {{.5, 0.}}), std::invalid_argument);
    CHECK_THROWS_AS(RelaxationHistory(1000., 20., {{.5, nan}}), std::invalid_argument);
    CHECK_THROWS_AS(RelaxationHistory(1000., 20., std::vector<PronyTerm>(65, {.001, 1.})), std::invalid_argument);
    RelaxationHistory history(1000., 20., {{.5, 2.}});
    history.advance(.01, 0., 20.);
    CHECK_THROWS_AS(history.advance(.1, -1., 20.), std::invalid_argument);
    CHECK_THROWS_AS(history.advance(.1, 1., 55.), std::invalid_argument);
    CHECK_THROWS_AS(history.advance(nan, 1., 20.), std::invalid_argument);
    CHECK_THAT(history.stress_pa(), WithinAbs(10., 1e-12));
    CHECK_THROWS_AS(history.advance(std::numeric_limits<double>::max(), 0., 20.), std::overflow_error);
    CHECK_THAT(history.advance(.01, 2., 20.), WithinAbs(5. + 5. * std::exp(-1.), 1e-12));
}

TEST_CASE("Published room temperature TPU spectrum reproduces its normalized relaxation", "[ThermoelasticRelaxation]")
{
    // Filaflex 82A at 20 C, Polymers 2021 13(14)2365 Table 7.
    // This checks the integrator, not applicability to another grade or temperature.
    const std::vector<PronyTerm> terms{{.14292,.054063}, {.021163,.14614}, {.029420,.39504},
        {.035212,1.0678}, {.013225,2.8866}, {.036500,7.8028}, {.00080567,21.092}, {.063713,57.015}};
    RelaxationHistory history(1., 20., terms);
    history.advance(1., 0., 20.);
    double expected = 1.;
    for (const auto &term : terms) expected -= term.fraction * (1. - std::exp(-10. / term.relaxation_seconds));
    CHECK_THAT(history.advance(1., 10., 20.), WithinAbs(expected, 1e-13));
    // Paper reports E(0)/E(10)=1.390 +/- .018 across specimens.
    CHECK_THAT(1. / history.stress_pa(), WithinAbs(1.390, .018));
}
