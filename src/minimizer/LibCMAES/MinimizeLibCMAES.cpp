// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#include <BSMPT/minimizer/LibCMAES/MinimizeLibCMAES.h>
#include <BSMPT/minimizer/MinimizePlane.h>
#include <BSMPT/models/ClassPotentialOrigin.h>

#include <libcmaes/candidate.h>       // for Candidate
#include <libcmaes/cmaes.h>           // for cmaes
#include <libcmaes/cmaparameters.h>   // for CMAParameters
#include <libcmaes/cmasolutions.h>    // for CMASolutions
#include <libcmaes/esoptimizer.h>     // for aCMAES
#include <libcmaes/esostrategy.h>     // for FitFunc
#include <libcmaes/genopheno.h>       // for GenoPheno
#include <libcmaes/noboundstrategy.h> // for libcmaes

namespace BSMPT
{
namespace Minimizer
{
namespace LibCMAES
{

using namespace libcmaes;

LibCMAESReturn min_cmaes_gen_all(const Class_Potential_Origin &model,
                                 const double &Temp,
                                 const std::vector<double> &Start)
{

  const auto dim = model.get_nVEV();

  std::vector<double> x0 = Start;

  double sigma;
  sigma = 5;

  // sigma *= 0.5;//0.5;
  double ftol = 1e-5;

  // CMAParameters<> cmaparams(dim,&x0.front(),sigma);
  CMAParameters<> cmaparams(x0, sigma, 1);

  // cmaparams.set_mt_feval(true);
  cmaparams.set_algo(aCMAES);
  // cmaparams.set_elitism(1);
  // cmaparams.set_noisy();
  cmaparams.set_ftolerance(ftol);

  FitFunc cmafunc = [&](const double *v, const int &N)
  {
    (void)N;
    std::vector<double> vev;
    for (std::size_t i{0}; i < dim; ++i)
      vev.push_back(v[i]);
    auto minOrdVEV = model.MinimizeOrderVEV(vev);
    auto VeffVal   = model.VEff(minOrdVEV, Temp);
    return VeffVal;
  };

  CMASolutions cmasols = cmaes<>(cmafunc, cmaparams);

  Candidate bcand = cmasols.best_candidate();

  std::vector<double> xsol = bcand.get_x();

  std::vector<double> sol;

  for (std::size_t i = 0; i < dim; i++)
  {
    sol.push_back(xsol.at(i));
  }

  LibCMAESReturn res;
  res.CMAESStatus = cmasols.run_status();
  res.result      = sol;

  return res;
}

LibCMAESReturn
CMAES_Minimize_Plane_gen_all(const struct PointerContainerMinPlane &params,
                             const std::vector<double> &Start)
{

  std::shared_ptr<Class_Potential_Origin> modelPointer;

  modelPointer = params.modelPointer;

  auto dim = params.nVEV - 1;
  std::vector<double> x0(dim);

  double sigma;
  sigma = 5;

  for (std::size_t i = 0; i < dim; i++)
    x0[i] = Start.at(i);

  double ftol = 1e-5;

  CMAParameters<> cmaparams(x0, sigma);

  cmaparams.set_algo(aCMAES);
  cmaparams.set_ftolerance(ftol);

  FitFunc CMAES_VEFF_Minimize_Plane = [=](const double *v, const int &N)
  {
    (void)N;
    std::size_t nVEVs = params.modelPointer->get_nVEV();

    std::vector<double> vMinTilde;
    vMinTilde.resize(nVEVs - 1);
    for (std::size_t i = 0; i < nVEVs - 1; i++)
      vMinTilde[i] = v[i];
    auto vev = TransformCoordinates(vMinTilde, params);

    double res = params.modelPointer->VEff(
        params.modelPointer->MinimizeOrderVEV(vev), params.Temp);
    return res;
  };

  CMASolutions cmasols = cmaes<>(CMAES_VEFF_Minimize_Plane, cmaparams);

  Candidate bcand = cmasols.best_candidate();

  std::vector<double> xsol = bcand.get_x();

  std::vector<double> sol;

  for (std::size_t i = 0; i < dim; i++)
  {
    sol.push_back(xsol.at(i));
  }

  LibCMAESReturn res;
  res.result      = TransformCoordinates(sol, params);
  res.CMAESStatus = cmasols.run_status();

  return res;
}

} // namespace LibCMAES
} // namespace Minimizer
} // namespace BSMPT
