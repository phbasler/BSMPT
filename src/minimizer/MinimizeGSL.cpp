/*
 * MinimizeGSL.cpp
 *
 *  Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file
 * Using the Nelder-Mead Simplex algorithm, implemented in gsl, to find multiple
 * local minima of the model and compare them to find a candidate for the global
 * minimum.
 */

#include <BSMPT/minimizer/MinimizeGSL.h>

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <algorithm>                           // for copy, max
#include <gsl/gsl_errno.h>                     // for gsl_set_error_handler...
#include <gsl/gsl_multimin.h>                  // for gsl_multimin_fminimizer
#include <gsl/gsl_vector_double.h>             // for gsl_vector_get, gsl_v...
#include <iostream>                            // for operator<<, endl, bas...
#include <limits>                              // for numeric_limits, numer...
#include <memory>                              // for shared_ptr, __shared_...
#include <random>                              // for default_random_engine
#include <stdio.h>                             // for printf
#include <time.h>                              // for time, NULL, std::size_t
#include <vector>                              // for vector

#include <atomic>
#include <mutex>
#include <queue>
#include <thread>

namespace BSMPT
{
namespace Minimizer
{

double GSL_VEFF_gen_all(const gsl_vector *v, void *p)
{

  struct GSL_params *params = static_cast<GSL_params *>(p);

  std::vector<double> vMin;
  auto nVEVs = params->modelPointer->get_nVEV();

  for (std::size_t i = 0; i < nVEVs; i++)
  {
    vMin.push_back(gsl_vector_get(v, i));
  }

  auto vIn = params->modelPointer->MinimizeOrderVEV(vMin);

  double res = params->modelPointer->VEff(vIn, params->Temp, 0);

  return res;
}

int GSL_Minimize_From_S_gen_all(struct GSL_params &params,
                                std::vector<double> &sol,
                                const std::vector<double> &start)
{
  gsl_set_error_handler_off();

  //	struct GSL_params * params = (struct GSL_params *) p;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s            = nullptr;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  double ftol         = GSL_Tolerance;
  std::size_t MaxIter = 600;

  std::size_t iter = 0;
  int status;
  double size;

  std::size_t dim = params.modelPointer->get_nVEV();

  /* Starting point */
  x = gsl_vector_alloc(dim);
  for (std::size_t k = 0; k < dim; k++)
    gsl_vector_set(x, k, start.at(k));
  ss = gsl_vector_alloc(dim);
  gsl_vector_set_all(ss, 1.0);

  /* Initialize method and iterate */
  minex_func.n      = dim;
  minex_func.f      = &GSL_VEFF_gen_all;
  minex_func.params = &params;
  s                 = gsl_multimin_fminimizer_alloc(T, dim);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  do
  {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    if (status) break;

    size   = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, ftol);

  } while (status == GSL_CONTINUE && iter < MaxIter);

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

std::pair<std::vector<double>, bool> GSL_Minimize_gen_all(
    const std::shared_ptr<Class_Potential_Origin> &modelPointer,
    const double &Temp,
    const int &seed,
    const std::size_t &MaxSol)
{
  std::vector<std::vector<double>> saveAllMinima;
  auto result =
      GSL_Minimize_gen_all(modelPointer, Temp, seed, saveAllMinima, MaxSol);
  return result;
}

std::pair<std::vector<double>, bool> GSL_Minimize_gen_all(
    const std::shared_ptr<Class_Potential_Origin> &modelPointer,
    const double &Temp,
    const int &seed)
{
  std::vector<std::vector<double>> saveAllMinima;
  std::size_t MaxSol = 20;
  auto result =
      GSL_Minimize_gen_all(modelPointer, Temp, seed, saveAllMinima, MaxSol);
  return result;
}

std::pair<std::vector<double>, bool> GSL_Minimize_gen_all(
    const std::shared_ptr<Class_Potential_Origin> &modelPointer,
    const double &Temp,
    const int &seed,
    std::vector<std::vector<double>> &saveAllMinima,
    const std::size_t &MaxSol)
{
  struct GSL_params params(modelPointer, Temp);

  std::size_t dim = modelPointer->get_nVEV();

  std::default_random_engine randGen(seed);
  double RNDMax        = 500;
  std::size_t MaxTries = 600; // 600;
  std::atomic<std::size_t> tries{0};
  std::size_t nCol = dim + 2;

  std::queue<std::vector<double>> StartingPoints;
  for (std::size_t i{0}; i < MaxTries; ++i)
  {
    std::vector<double> start(dim);
    for (std::size_t j = 0; j < dim; ++j)
    {
      start.at(j) =
          RNDMax *
          (-1 +
           2 * std::generate_canonical<double,
                                       std::numeric_limits<double>::digits>(
                   randGen));
    }
    StartingPoints.push(start);
  }

  struct MinimizeParams
  {
    std::vector<double> start;
    std::vector<double> sol;
    int GSL_status;
  };

  std::atomic<std::size_t> FoundSolutions{0};
  std::vector<std::thread> MinThreads;
  std::mutex WriteResultLock;
  std::queue<std::vector<double>> Results;

  auto thread_Job = [](std::atomic<std::size_t> &mFoundSolutions,
                       std::size_t mMaxSol,
                       std::queue<std::vector<double>> &mStartingPoints,
                       GSL_params &mparams,
                       std::mutex &mWriteResultLock,
                       std::queue<std::vector<double>> &mResults) {
    while (mFoundSolutions < mMaxSol and not mStartingPoints.empty())
    {
      std::vector<double> start;
      {
        std::unique_lock<std::mutex> lock(mWriteResultLock);
        start = mStartingPoints.front();
        mStartingPoints.pop();
      }

      std::vector<double> sol;
      auto status = GSL_Minimize_From_S_gen_all(mparams, sol, start);
      if (status == GSL_SUCCESS)
      {
        std::unique_lock<std::mutex> lock(mWriteResultLock);
        ++mFoundSolutions;
        mResults.push(sol);
      }
    }
  };

  for (std::size_t i = 0; i < Num_threads; ++i)
  {
    MinThreads.push_back(std::thread([&]() {
      thread_Job(FoundSolutions,
                 MaxSol,
                 StartingPoints,
                 params,
                 WriteResultLock,
                 Results);
    }));
  }

  for (auto &thr : MinThreads)
  {
    if (thr.joinable()) thr.join();
  }

  while (not Results.empty())
  {
    auto res = Results.front();
    Results.pop();
    auto vpot = modelPointer->MinimizeOrderVEV(res);
    std::vector<double> row(nCol);
    for (std::size_t i = 0; i < dim; ++i)
      row.at(i) = res.at(i);
    row.at(dim)     = modelPointer->EWSBVEV(vpot);
    row.at(dim + 1) = modelPointer->VEff(vpot, Temp, 0);
    saveAllMinima.push_back(row);
  }

  if (saveAllMinima.size() == 0)
  {
    std::cerr << "No solutions found during the GSL minimization at T = "
              << Temp << " GeV "
              << "\n";
    return std::make_pair(std::vector<double>{}, false);
  }

  if (saveAllMinima.size() < MaxSol)
  {
    std::cerr << "Found " << saveAllMinima.size() << " of  " << MaxSol
              << " solutions at T = " << Temp << std::endl;
  }

  std::size_t minIndex = 0;
  double VMin          = saveAllMinima[0][dim + 1];
  for (std::size_t k = 1; k < saveAllMinima.size(); k++)
  {
    if (saveAllMinima[k][dim + 1] <= VMin)
    {
      VMin     = saveAllMinima[k][dim + 1];
      minIndex = k;
    }
  }

  std::vector<double> solV;
  for (std::size_t k = 0; k < dim; k++)
    solV.push_back(saveAllMinima[minIndex][k]);
  return std::make_pair(solV, true);
}

} // namespace Minimizer
} // namespace BSMPT
