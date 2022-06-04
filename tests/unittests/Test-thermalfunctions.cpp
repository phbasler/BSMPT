// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Approx = Catch::Approx;
#include <BSMPT/ThermalFunctions/thermalcoefficientcalculator.h>
#include <cmath>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_zeta.h>

TEST_CASE("Check FermionInterpolatedLowCoefficentCalculator", "[thermal]")
{
  BSMPT::ThermalFunctions::ThermalCoefficientCalculator
  FermionInterpolatedLowCoefficientCalculator(
      [](std::size_t l) -> double {
        if (l < 2)
        {
          return 0;
        }
        else
        {
          return gsl_sf_doublefact(2 * l - 3) * gsl_sf_zeta(2 * l - 1) /
                 (gsl_sf_doublefact(2 * l) * (l + 1)) *
                 (std::pow(2, 2 * l - 1) - 1);
        }
      },
      5);

  std::map<std::size_t, double> expected{{0, 0},
                                         {1, 0},
                                         {2, .3505999301},
                                         {3, .5022618813},
                                         {4, 1.000471548},
                                         {5, 2.333453140}};

  for (const auto &el : expected)
  {
    REQUIRE(FermionInterpolatedLowCoefficientCalculator.GetCoefficentAtOrder(
                el.first) == Approx(el.second).margin(1e-4));
  }
}

TEST_CASE("Check BosonInterpolatedLowCoefficientCalculator", "[thermal]")
{
  BSMPT::ThermalFunctions::ThermalCoefficientCalculator
  BosonInterpolatedLowCoefficientCalculator(
      [](std::size_t l) -> double {
        if (l < 2)
        {
          return 0;
        }
        else
        {
          return gsl_sf_doublefact(2 * l - 3) * gsl_sf_zeta(2 * l - 1) /
                 (gsl_sf_doublefact(2 * l) * (l + 1));
        }
      },
      5);

  std::map<std::size_t, double> expected{{0, 0},
                                         {1, 0},
                                         {2, 0.5008570430e-1},
                                         {3, 0.1620199617e-1},
                                         {4, 0.7877728727e-2},
                                         {5, 0.4566444500e-2}};

  for (const auto &el : expected)
  {
    REQUIRE(BosonInterpolatedLowCoefficientCalculator.GetCoefficentAtOrder(
                el.first) == Approx(el.second).margin(1e-4));
  }
}

TEST_CASE("Check JInterpolatedHighCoefficientCalculator", "[thermal]")
{
  BSMPT::ThermalFunctions::ThermalCoefficientCalculator
  JInterpolatedHighCoefficientCalculator(
      [](std::size_t l) -> double {
        return 1 / (std::pow(2, l) * gsl_sf_fact(l)) * gsl_sf_gamma(2.5 + l) /
               gsl_sf_gamma(2.5 - l);
      },
      5);

  std::map<std::size_t, double> expected{{0, 1.},
                                         {1, 1.875000000},
                                         {2, .8203125000},
                                         {3, -.3076171875},
                                         {4, .3172302246},
                                         {5, -.5154991150}};

  for (const auto &el : expected)
  {
    REQUIRE(JInterpolatedHighCoefficientCalculator.GetCoefficentAtOrder(
                el.first) == Approx(el.second).margin(1e-4));
  }
}
