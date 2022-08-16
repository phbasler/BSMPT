// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#include <BSMPT/WallThickness/WallThicknessCommon.h>
#include <BSMPT/WallThickness/WallThicknessLib.h>
#include <BSMPT/minimizer/MinimizePlane.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/utility/utility.h>

#include <gsl/gsl_min.h>

#include <BSMPT/models/IncludeAllModels.h>
#include <boost/math/tools/minima.hpp>
#include <fstream>
#include <random>

#include <atomic>
#include <mutex>
#include <queue>
#include <thread>

namespace BSMPT
{
namespace Wall
{

namespace
{
const double GSL_Tolerance = std::pow(10, -4);
}

double GSL_VEFF_gen_all_maximum_line(double t, void *p)
{
  struct GSL_params *params = static_cast<GSL_params *>(p);
  std::vector<double> vIn, vMin;
  std::size_t nVEVs = params->modelPointer->get_nVEV();

  for (std::size_t i = 0; i < nVEVs; i++)
  {
    vMin.push_back(t * params->VevMinimum.at(i) +
                   (1 - t) * params->VeVSymmetric.at(i));
  }

  vIn.resize(params->modelPointer->get_NHiggs());
  vIn = params->modelPointer->MinimizeOrderVEV(vMin);

  double res =
      -1 * params->modelPointer->VEff(
               vIn, params->Temp, 0); // -1 makes the maximum the minimum
  return res;
}

double calculate_wall_thickness_1D(
    const std::shared_ptr<Class_Potential_Origin> &modelPointer,
    const double &Temp,
    const std::vector<double> &vcritical,
    const std::vector<double> &vevsymmetric)
{
  std::vector<double> vbarrier;
  GSL_Find_Maximum_line(modelPointer, Temp, vcritical, vevsymmetric, vbarrier);
  double LW1D = 0, Vb1D = 0;
  Vb1D    = vbarrier.at(vcritical.size());
  auto vc = modelPointer->EWSBVEV(modelPointer->MinimizeOrderVEV(vcritical));
  if (Vb1D != 0) LW1D = vc / std::sqrt(8 * Vb1D);
  return LW1D;
}

bool GSL_Find_Maximum_line(
    const std::shared_ptr<Class_Potential_Origin> &modelPointer,
    const double &Temp,
    const std::vector<double> &vcritical,
    const std::vector<double> &vevsymmetric,
    std::vector<double> &solV)
{

  struct GSL_params params;
  params.Temp         = Temp;
  params.nVEV         = modelPointer->get_nVEV();
  params.modelPointer = modelPointer;
  params.VevMinimum   = vcritical;
  params.VeVSymmetric = vevsymmetric;

  double tmax = 0;

  std::vector<std::vector<double>> solutions;
  std::size_t nSolutions = 20;
  int ntry               = 0;
  int seed               = 5;
  std::default_random_engine randGen(seed);
  do
  {
    double initial_guess =
        std::generate_canonical<double, std::numeric_limits<double>::digits>(
            randGen);
    GSL_Maximize_From_S_gen_line(params, solutions, initial_guess);
    ntry++;
  } while (solutions.size() != nSolutions);

  tmax        = solutions.at(0).at(0);
  double fmax = solutions.at(0).at(1);
  for (std::size_t i = 1; i < solutions.size(); i++)
  {
    if (solutions.at(i).at(1) > fmax)
    {
      fmax = solutions.at(i).at(1);
      tmax = solutions.at(i).at(0);
    }
  }

  if (solV.size() != 0 and solV.size() != vcritical.size())
  {
    std::cerr << "The length of solV in " << __func__
              << " does not match with vcritical and is not zero." << std::endl;
  }
  else if (solV.size() != 0)
  {
    for (std::size_t i = 0; i < solV.size(); i++)
    {
      solV.at(i) = tmax * vcritical.at(i);
    }
  }
  else if (solV.size() == 0)
  {
    for (std::size_t i = 0; i < vcritical.size(); i++)
    {
      solV.push_back(tmax * vcritical.at(i));
    }
  }
  std::vector<double> vIn;

  vIn.resize(params.modelPointer->get_NHiggs());
  vIn = params.modelPointer->MinimizeOrderVEV(solV);

  double MaxVal = params.modelPointer->VEff(vIn, params.Temp, 0);

  vIn            = params.modelPointer->MinimizeOrderVEV(vcritical);
  double MinValC = params.modelPointer->VEff(vIn, params.Temp, 0);
  for (std::size_t i = 0; i < vIn.size(); i++)
    vIn.at(i) = 0;
  double MinVal0 = params.modelPointer->VEff(vIn, params.Temp, 0);

  double MinVal =
      std::min(MinVal0, MinValC); // They can differ through the numerical error
                                  // of the minimization

  double Vb = MaxVal - MinVal;

  solV.push_back(Vb);
  return true;
}

bool GSL_Maximize_From_S_gen_line(struct GSL_params &params,
                                  std::vector<std::vector<double>> &solution,
                                  double initial_guess)
{
  gsl_set_error_handler_off();

  gsl_min_fminimizer *s;
  const gsl_min_fminimizer_type *T;

  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc(T);

  std::size_t MaxIter = 600;

  std::size_t iter = 0;
  int status;

  gsl_function F;
  F.function = &GSL_VEFF_gen_all_maximum_line;
  F.params   = &params;

  gsl_min_fminimizer_set(s, &F, initial_guess, 0, 1);

  double xsol, xlow, xup, fval;

  do
  {
    iter++;
    status = gsl_min_fminimizer_iterate(s);
    if (status) break;
    xlow   = gsl_min_fminimizer_x_lower(s);
    xup    = gsl_min_fminimizer_x_upper(s);
    xsol   = gsl_min_fminimizer_x_minimum(s);
    fval   = gsl_min_fminimizer_f_minimum(s);
    status = gsl_min_test_interval(xlow, xup, GSL_Tolerance, 0);
  } while (status == GSL_CONTINUE && iter < MaxIter);

  std::vector<double> row;

  if (status == GSL_SUCCESS)
  {
    row.push_back(xsol);
    row.push_back(fval);
    solution.push_back(row);
  }

  gsl_min_fminimizer_free(s);

  return status == GSL_SUCCESS;
}

} // namespace Wall
} // namespace BSMPT
