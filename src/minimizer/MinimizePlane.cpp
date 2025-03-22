// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#include <BSMPT/minimizer/MinimizePlane.h>
#include <BSMPT/utility/Logger.h>
#include <gsl/gsl_errno.h>    // for gsl_set_error_handler...
#include <gsl/gsl_multimin.h> // for gsl_multimin_fminimizer

#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>     // for FChoose
#include <iomanip>                             // for setprecision
#include <iostream>                            // for operator<<, endl, ost...
#include <limits>                              // for numeric_limits, numer...
#include <random>                              // for default_random_engine
#include <stdexcept>                           // for runtime_error
#include <stdio.h>                             // for printf
#include <string>                              // for string, operator+
#include <time.h>                              // for std::size_t, time, NULL

#include <BSMPT/config.h>

#ifdef cmaes_FOUND
#include <BSMPT/minimizer/LibCMAES/MinimizeLibCMAES.h>
#endif

#ifdef NLopt_FOUND
#include <BSMPT/minimizer/LibNLOPT/MinimizeNLOPT.h>
#endif

namespace BSMPT
{
namespace Minimizer
{

std::vector<double>
TransformCoordinates(const std::vector<double> &vMinTilde,
                     const struct PointerContainerMinPlane &params)
{

  if (vMinTilde.size() != params.nVEV - 1)
  {
    std::string throwstring = "No valid input in ";
    throwstring += __func__;
    throwstring += ". Expected " + std::to_string(params.nVEV - 1);
    throwstring += " but found " + std::to_string(vMinTilde.size()) + ".";
    throw std::runtime_error(throwstring.c_str());
  }

  std::vector<double> vMin;

  std::size_t nVEVs = params.modelPointer->get_nVEV();
  vMin.resize(nVEVs);
  for (std::size_t i = 0; i < params.Index; i++)
  {
    vMin[i] = vMinTilde.at(i);
  }
  vMin[params.Index] = 0;
  for (std::size_t i = 1 + params.Index; i < nVEVs; i++)
  {
    vMin[i] = vMinTilde.at(i - 1);
  }

  double plane_const = 0, plane_coord = 0;
  for (std::size_t i = 0; i < nVEVs; i++)
  {
    plane_const += (params.normalvector.at(i)) * (params.Point.at(i));
    plane_coord += (params.normalvector.at(i)) * (vMin.at(i));
  }
  vMin[params.Index] = 1.0 / ((params.normalvector.at(params.Index))) *
                       (plane_const - plane_coord);

  return vMin;
}

MinimizePlaneReturn MinimizePlane(const std::vector<double> &basepoint,
                                  const std::vector<double> &VEVSymmetric,
                                  const std::vector<double> &VEVBroken,
                                  const ModelID::ModelIDs &Model,
                                  const std::vector<double> &par,
                                  const std::vector<double> &parCT,
                                  const double &Temp,
                                  const int &WhichMinimizer)
{
  return MinimizePlane(basepoint,
                       VEVSymmetric,
                       VEVBroken,
                       Model,
                       par,
                       parCT,
                       GetSMConstants(),
                       Temp,
                       WhichMinimizer);
}

MinimizePlaneReturn MinimizePlane(const std::vector<double> &basepoint,
                                  const std::vector<double> &VEVSymmetric,
                                  const std::vector<double> &VEVBroken,
                                  const ModelID::ModelIDs &Model,
                                  const std::vector<double> &par,
                                  const std::vector<double> &parCT,
                                  const ISMConstants &SMConstant,
                                  const double &Temp,
                                  const int &WhichMinimizer)
{
  std::shared_ptr<Class_Potential_Origin> modelPointer =
      ModelID::FChoose(Model, SMConstant);
  modelPointer->set_All(par, parCT);
  return MinimizePlane(
      basepoint, VEVSymmetric, VEVBroken, modelPointer, Temp, WhichMinimizer);
}

MinimizePlaneReturn
MinimizePlane(const std::vector<double> &basepoint,
              const std::vector<double> &VEVSymmetric,
              const std::vector<double> &VEVBroken,
              const std::shared_ptr<Class_Potential_Origin> &modelPointer,
              const double &Temp,
              const int &WhichMinimizer)
{

  const auto UseMinimizer = GetMinimizers(WhichMinimizer);

  std::vector<double> PotValues;
  std::vector<std::vector<double>> Minima;

  std::vector<double> normalvector;
  for (std::size_t i = 0; i < VEVBroken.size(); i++)
    normalvector.push_back(VEVBroken.at(i) - VEVSymmetric.at(i));

  struct PointerContainerMinPlane params;
  params.Temp         = Temp;
  params.nVEV         = modelPointer->get_nVEV();
  params.modelPointer = modelPointer;
  params.normalvector = normalvector;
  params.Point        = basepoint;
  params.VEVSymmetric = VEVSymmetric;
  params.VEVBroken    = VEVBroken;

  int IndexNorm = -1;
  for (std::size_t i = 0; i < normalvector.size(); i++)
  {
    if (std::abs(normalvector.at(i)) >= 1e-3)
    {
      IndexNorm = static_cast<int>(i);
      break;
    }
  }
  if (IndexNorm == -1)
  {
    throw std::runtime_error(
        "The normalvector in the minimization of the plane is zero. \
                                 This means the broken and symmetric minimum are identical and\
                                  no tunnelpath exists");
  }

  params.Index = static_cast<size_t>(IndexNorm);

  auto dimensionnames = modelPointer->addLegendTemp();

  if (UseMinimizer.UseGSL)
  {
    // Find the minimum provided by GSL
    auto GSLResult = GSL_Minimize_Plane_gen_all(params, 3, 50);
    PotValues.push_back(GSLResult.PotVal);
    Minima.push_back(GSLResult.Minimum);
  }

#ifdef cmaes_FOUND
  if (UseMinimizer.UseCMAES and modelPointer->get_nVEV() >= 3)
  {
    std::vector<double> startCMAES(params.nVEV - 1);
    for (std::size_t i = 0; i < params.Index; i++)
    {
      startCMAES.at(i) = params.VEVBroken.at(i);
    }
    for (std::size_t i = static_cast<size_t>(params.Index); i < params.nVEV - 1;
         i++)
    {
      startCMAES.at(i) = params.VEVBroken.at(i + 1);
    }
    auto LibCMAESResult =
        LibCMAES::CMAES_Minimize_Plane_gen_all(params, startCMAES);
    PotValues.push_back(modelPointer->VEff(
        modelPointer->MinimizeOrderVEV(LibCMAESResult.result), params.Temp));
    Minima.push_back(LibCMAESResult.result);
  }
#endif

#ifdef NLopt_FOUND
  if (UseMinimizer.UseNLopt)
  {
    auto NLOPTResult = LibNLOPT::MinimizePlaneUsingNLOPT(params);
    PotValues.push_back(NLOPTResult.PotVal);
    Minima.push_back(NLOPTResult.Minimum);
  }
#endif

  std::size_t MinIndex = 0;
  for (std::size_t i = 1; i < Minima.size(); i++)
  {
    if (PotValues.at(i) < PotValues.at(MinIndex))
    {
      MinIndex = i;
    }
  }

  MinimizePlaneReturn res;
  res.PotVal  = PotValues.at(MinIndex);
  res.Minimum = Minima.at(MinIndex);

  return res;
}

double GSL_VEFF_Minimize_Plane(const gsl_vector *v, void *p)
{
  struct PointerContainerMinPlane *params =
      static_cast<PointerContainerMinPlane *>(p);
  std::vector<double> vMinTilde;
  std::size_t nVEVs = params->nVEV;
  for (std::size_t i = 0; i < nVEVs - 1; i++)
  {
    vMinTilde.push_back(gsl_vector_get(v, i));
  }
  auto vMin  = TransformCoordinates(vMinTilde, *params);
  auto vIn   = params->modelPointer->MinimizeOrderVEV(vMin);
  double res = params->modelPointer->VEff(vIn, params->Temp, 0);
  return res;
}

int GSL_Minimize_Plane_From_S_gen_all(
    const struct PointerContainerMinPlane &params,
    std::vector<double> &sol,
    const std::vector<double> &start)
{
  const double GSL_Tolerance = std::pow(10, -6);
  gsl_set_error_handler_off();

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s            = nullptr;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  const double ftol   = GSL_Tolerance;
  std::size_t MinIter = 600;

  std::size_t iter = 0;
  int status;
  double size;

  std::size_t dim = params.nVEV - 1;

  /* Starting point */
  x = gsl_vector_alloc(dim);
  for (std::size_t k = 0; k < dim; k++)
    gsl_vector_set(x, k, start.at(k));
  ss = gsl_vector_alloc(dim);
  gsl_vector_set_all(ss, 1.0);

  struct PointerContainerMinPlane params_nonconst = params;
  /* Initialize method and iterate */
  minex_func.n      = dim;
  minex_func.f      = &GSL_VEFF_Minimize_Plane;
  minex_func.params = &params_nonconst;
  s                 = gsl_multimin_fminimizer_alloc(T, dim);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  do
  {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    if (status) break;

    size   = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, ftol);
  } while (status == GSL_CONTINUE && iter < MinIter);

  if (status == GSL_SUCCESS)
  {
    for (std::size_t k = 0; k < dim; k++)
      sol.push_back(gsl_vector_get(s->x, k));
  }
  else
  {
    for (std::size_t k = 0; k < dim; k++)
      sol.push_back(0);
  }

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);

  return status;
}

GSLPlaneReturn
GSL_Minimize_Plane_gen_all(const struct PointerContainerMinPlane &params,
                           std::size_t seed,
                           std::size_t MinSol)
{
  std::vector<std::vector<double>> saveAllMinima;
  auto result = GSL_Minimize_Plane_gen_all(params, seed, saveAllMinima, MinSol);
  return result;
}

GSLPlaneReturn
GSL_Minimize_Plane_gen_all(const struct PointerContainerMinPlane &params,
                           std::size_t seed)
{
  std::vector<std::vector<double>> saveAllMinima;
  std::size_t MinSol = 20;
  auto result = GSL_Minimize_Plane_gen_all(params, seed, saveAllMinima, MinSol);
  return result;
}

GSLPlaneReturn
GSL_Minimize_Plane_gen_all(const struct PointerContainerMinPlane &params,
                           std::size_t seed,
                           std::vector<std::vector<double>> &saveAllMinima,
                           std::size_t MinSol)
{
  GSLPlaneReturn res;

  std::shared_ptr<Class_Potential_Origin> modelPointer;
  modelPointer = params.modelPointer;

  std::size_t dim = modelPointer->get_nVEV() - 1;

  std::default_random_engine randGen(seed);
  double RNDMin        = 500;
  std::size_t MinTries = 600;
  std::size_t tries    = 0;
  std::size_t numOfSol = 0;
  std::size_t nCol     = dim + 2 + 1;
  std::vector<double> start, vPot, soltilde;
  do
  {
    start.resize(dim);
    for (std::size_t i = 0; i < dim; i++)
      start[i] =
          RNDMin *
          (-1 +
           2 * std::generate_canonical<double,
                                       std::numeric_limits<double>::digits>(
                   randGen));

    GSL_Minimize_Plane_From_S_gen_all(params, soltilde, start);

    vPot.resize(modelPointer->get_NHiggs());
    auto sol = TransformCoordinates(soltilde, params);
    vPot     = modelPointer->MinimizeOrderVEV(sol);

    std::vector<double> row(nCol);
    for (std::size_t i = 0; i < dim + 1; i++)
      row.at(i) = sol.at(i);
    row.at(dim + 1) = modelPointer->EWSBVEV(vPot);
    row.at(dim + 2) = modelPointer->VEff(vPot, params.Temp, 0);

    saveAllMinima.push_back(row);

    numOfSol++;
    if (numOfSol == MinSol) break;

    start.clear();
    sol.clear();
    soltilde.clear();
    tries++;
  } while (tries <= MinTries);
  if (numOfSol == 0)
  {
    Logger::Write(LoggingLevel::MinimizerDetailed,
                  "No solutions found during the GSL minimization at T = " +
                      std::to_string(params.Temp) + " GeV ");
    res.StatusFlag = false;
    return res;
  }

  std::size_t minIndex = 0;
  double VMin          = saveAllMinima[0][dim + 2];
  for (std::size_t k = 1; k < numOfSol; k++)
  {
    if (saveAllMinima[k][dim + 2] <= VMin)
    {
      VMin     = saveAllMinima[k][dim + 2];
      minIndex = k;
    }
  }
  res.StatusFlag = true;
  res.PotVal     = saveAllMinima[minIndex][dim + 2];
  res.vc         = saveAllMinima[minIndex][dim + 1];
  std::vector<double> Minimum;
  for (std::size_t i = 0; i < dim + 1; i++)
    Minimum.push_back(saveAllMinima[minIndex][i]);
  res.Minimum = Minimum;
  return res;
}

} // namespace Minimizer
} // namespace BSMPT
