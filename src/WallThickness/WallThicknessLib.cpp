/*
 * WallThicknessLib.cpp
 *
 *  Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller

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
 */

#include <BSMPT/WallThickness/WallThicknessLib.h>
#include <BSMPT/minimizer/MinimizePlane.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/utility.h>

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
const double GSL_Tolerance    = std::pow(10, -4);
const std::size_t Num_threads = std::thread::hardware_concurrency();

/**
 * struct containing the required Parameters of the model for the gsl interface
 */
struct GSL_params
{
  /**
   * @brief nVEV number of VEVs
   */
  std::size_t nVEV;
  /**
   * @brief VevMinimum electroweak minimum in the broken phase
   */
  std::vector<double> VevMinimum;
  /**
   * @brief VeVSymmetric electroweak minimum in the symmetric phase
   */
  std::vector<double> VeVSymmetric;
  /**
   * @brief modelPointer shared_ptr for the parameter point
   */
  std::shared_ptr<Class_Potential_Origin> modelPointer;
  /**
   * @brief Temp temperature at which to evaluate the parameter point
   */
  double Temp;
  /**
   * @brief spline cubic spline used to find the potential barrier
   */
  boost_cubic_b_spline<double> spline;
  /**
   * @brief UseSpline Decides if the spline is to be used or not
   */
  bool UseSpline = false;
};

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

double calculate_wall_thickness_plane(
    const std::shared_ptr<Class_Potential_Origin> &modelPointer,
    const double &Temp,
    const std::vector<double> &vcritical,
    const std::vector<double> &vevsymmetric,
    const int &WhichMinimizer)
{

  if (vcritical.size() != vevsymmetric.size())
  {
    throw std::runtime_error(
        "vcritical and vevsymmetric in calculate_wall_thickness_plane do not "
        "have the same size ");
  }
  double LW       = 0;
  int maxstep     = 10;
  double Stepsize = 1.0 / (maxstep);
  std::vector<double> Data_min_negative(maxstep + 1);
  std::vector<std::thread> Data_min_threads(Num_threads);

  std::atomic<std::size_t> DataIndex{0};
  std::mutex WriteLock;

  std::vector<std::pair<int, std::vector<double>>> BasePoints;
  for (int ncounter = 0; ncounter <= maxstep; ++ncounter)
  {
    double line_parameter = Stepsize * ncounter;
    std::vector<double> basepoint;
    for (std::size_t i = 0; i < vcritical.size(); i++)
    {
      basepoint.push_back(vevsymmetric.at(i) * (1 - line_parameter) +
                          vcritical.at(i) * line_parameter);
    }
    BasePoints.push_back(std::make_pair(ncounter, basepoint));
  }

  auto thread_Job =
      [](std::atomic<std::size_t> &mDataIndex,
         std::vector<std::pair<int, std::vector<double>>> &a_BasePoints,
         std::vector<double> &Data_min,
         std::mutex &mWriteResultLock,
         const std::vector<double> &a_vevsymmetric,
         const std::vector<double> &a_vcritical,
         std::shared_ptr<Class_Potential_Origin> a_modelPointer,
         const double &a_Temp,
         const int &a_WhichMinimizer) {
        while (mDataIndex < a_BasePoints.size())
        {
          const auto data         = a_BasePoints.at(mDataIndex++);
          auto MinimumPlaneResult = Minimizer::MinimizePlane(data.second,
                                                             a_vevsymmetric,
                                                             a_vcritical,
                                                             a_modelPointer,
                                                             a_Temp,
                                                             a_WhichMinimizer);
          {
            std::unique_lock<std::mutex> lock(mWriteResultLock);
            Data_min.at(data.first) = -MinimumPlaneResult.PotVal;
          }
        }
      };

  for (auto &thr : Data_min_threads)
  {
    thr = std::thread([&]() {
      thread_Job(DataIndex,
                 BasePoints,
                 Data_min_negative,
                 WriteLock,
                 vevsymmetric,
                 vcritical,
                 modelPointer,
                 Temp,
                 WhichMinimizer);
    });
  }

  for (auto &thr : Data_min_threads)
  {
    if (thr.joinable())
    {
      thr.join();
    }
  }

  struct GSL_params spline;
  boost_cubic_b_spline<double> splinef(
      Data_min_negative.data(), Data_min_negative.size(), 0, Stepsize);
  spline.spline    = splinef;
  spline.UseSpline = true;

  int precision = std::numeric_limits<double>::digits;
  std::pair<double, double> MaxPair =
      boost::math::tools::brent_find_minima(splinef, 0., 1., precision);

  //	double Max_pos = MaxPair.first;
  double VMax     = -MaxPair.second;
  double Vmin     = std::min(-splinef(0), -splinef(1));
  double Vbarrier = VMax - Vmin;

  auto vc = modelPointer->EWSBVEV(modelPointer->MinimizeOrderVEV(vcritical));

  LW = vc / std::sqrt(8 * Vbarrier);

  return LW;
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
  std::cout << std::scientific;

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
