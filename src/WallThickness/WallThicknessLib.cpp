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

const std::size_t Num_threads = std::thread::hardware_concurrency();

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

} // namespace Wall
} // namespace BSMPT
