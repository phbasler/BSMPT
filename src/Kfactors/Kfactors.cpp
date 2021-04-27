// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/Kfactors/Kfactors.h>
#include <cmath>
#include <iostream>

#include <gsl/gsl_integration.h>

/**
 * @file
 */

namespace BSMPT
{
namespace Kfactors
{

double distribution_f0(double E0, int s, double Temp, int diff)
{
  double beta = 1.0 / Temp;
  double res  = 0;
  double expv = std::exp(beta * E0);
  if (diff == 0)
  {
    res = 1.0 / (expv + s);
  }
  else if (diff == 1)
  {
    res = -beta * expv / std::pow(expv + s, 2);
  }
  else if (diff == 2)
  {
    res = beta * beta * expv * (expv - s) / std::pow(expv + s, 3);
  }
  return res;
}

double K_integrand(const std::vector<double> &p,
                   double masssquared,
                   int switchvalue,
                   int s,
                   double Temp)
{
  double res = 0;
  double px = p.at(0), py = p.at(1), pz = p.at(2);
  double E0z = std::sqrt(pz * pz + masssquared);
  double E0  = std::sqrt(px * px + py * py + pz * pz + masssquared);
  switch (switchvalue)
  {
  case 1: res = -pz * pz / E0 * distribution_f0(E0, s, Temp, 2); break;
  case 2: res = distribution_f0(E0, s, Temp, 2) / (2.0 * E0); break;
  case 3: res = distribution_f0(E0, s, Temp, 1) / (2.0 * E0); break;
  case 4: res = (pz * pz) / (E0 * E0) * distribution_f0(E0, s, Temp, 1); break;
  case 5: res = (pz * pz) / E0 * distribution_f0(E0, s, Temp, 1); break;
  case 6:
    res = (E0 * E0 - pz * pz) / (2.0 * std::pow(E0, 3)) *
          distribution_f0(E0, s, Temp, 1);
    break;
  case 7:
    res = std::abs(pz) / (2 * E0 * E0 * E0z) *
          (distribution_f0(E0, s, Temp, 1) / E0 -
           distribution_f0(E0, s, Temp, 2));
    break;
  case 8:
    res =
        std::abs(pz) * distribution_f0(E0, s, Temp, 1) / (2.0 * E0 * E0 * E0z);
    break;
  case 9:
    res = std::abs(pz) / (4 * std::pow(E0, 3) * E0z) *
          (distribution_f0(E0, s, Temp, 1) / E0 -
           distribution_f0(E0, s, Temp, 2));
    break;
  case 10:
    res = std::abs(pz) * distribution_f0(E0, s, Temp, 0) /
          (2.0 * std::pow(E0, 3) * E0z);
    break;
  default: std::cerr << "Wrong call for " << __func__ << std::endl; break;
  }

  return res;
}

double K_integrand_gsl(double *x, std::size_t dim, void *p)
{
  std::vector<double> momentum(dim);
  for (std::size_t i = 0; i < dim; i++)
    momentum[i] = x[i];
  struct GSL_integration *params =
      static_cast<GSL_integration *>(p); //(struct GSL_integration *) p;
  double res = K_integrand(momentum,
                           params->masssquared,
                           params->switchval,
                           params->s,
                           params->Temp);
  return res;
}

double K_integration(double masssquared, double Temp, int switchvalue, int s)
{

  struct GSL_integration p;
  p.Temp        = Temp;
  p.s           = s;
  p.switchval   = switchvalue;
  p.masssquared = masssquared;

  const std::size_t dim = 3;

  double numint = 1e3;

  double xl[dim] = {0, 0, 0}; // We integrate from 0 to + infinity and multiply
                              // by 2^3 using the symmetry p -> -p
  double xu[dim] = {numint, numint, numint};

  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_monte_function G;
  G.f      = &K_integrand_gsl;
  G.dim    = dim;
  G.params = &p;

  // std::size_t calls = 500000;
  std::size_t calls = 5e5;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  gsl_monte_vegas_state *state = gsl_monte_vegas_alloc(dim);

  double res, err;

  gsl_monte_vegas_integrate(&G, xl, xu, dim, 10000, r, state, &res, &err);

  do
  {
    gsl_monte_vegas_integrate(&G, xl, xu, 3, calls / 5, r, state, &res, &err);
  } while (std::fabs(gsl_monte_vegas_chisq(state) - 1) > 0.5);

  gsl_monte_vegas_free(state);
  gsl_rng_free(r);

  res *= 8;
  return res;
}

double K_functions(double masssquared, double Temp, int switchvalue, int s)
{
  double numerator = K_integration(masssquared, Temp, switchvalue, s);
  double norm      = 0;
  if (switchvalue > 10)
  {
    std::cerr << "Wrong switch value for K_functions ! " << std::endl;
    exit(EXIT_FAILURE);
  }
  else if (switchvalue == 5 or switchvalue == 6 or switchvalue == 10)
  {
    norm = Ktilde_normalization(Temp, s, masssquared);
  }
  else
  {
    norm = -std::pow(M_PI, 3) * std::pow(Temp, 2) * 2.0 / 3.0;
  }
  return numerator / norm;
}

void display_results(std::string title, double result, double error)
{
  std::cout << std::scientific;
  std::cout << title << "==================\n";
  // printf ("%s ==================\n", title);
  std::cout << "result = " << result << std::endl;
  std::cout << "sigma = " << error << std::endl;
}

double Ktilde_normalization_func(double x, void *p)
{
  struct GSL_integration *params =
      static_cast<GSL_integration *>(p); //(struct GSL_integration * )p;
  double a   = params->masssquared * std::pow(1.0 / params->Temp, 2);
  int s      = params->s;
  double res = std::pow(x, 2) / (std::exp(std::sqrt(a + std::pow(x, 2))) + s);
  return res;
}

double Ktilde_normalization(double Temp, int s, double masssquared)
{
  struct GSL_integration p;
  p.Temp        = Temp;
  p.s           = s;
  p.masssquared = masssquared;

  std::size_t workspace_size = 1000;
  gsl_integration_workspace *w =
      gsl_integration_workspace_alloc(workspace_size);
  double result, error;
  gsl_function F;
  F.function = &Ktilde_normalization_func;
  F.params   = &p;

  double xmax = 30; // around x = 30 the integrand becomes < 10^-10 for m^2 = 0.
                    // For m^2 > 0 this happens earlierer.
  gsl_integration_qags(
      &F, 0, xmax, 0, 1e-9, workspace_size, w, &result, &error);
  gsl_integration_workspace_free(w);

  result *= 4 * M_PI * std::pow(Temp, 3);
  return result;
}

} // namespace Kfactors
} // namespace BSMPT
