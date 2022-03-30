#ifndef WALLTHICKNESSCOMMON_H
#define WALLTHICKNESSCOMMON_H

// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#include <BSMPT/utility/utility.h>
#include <memory>
#include <vector>

namespace BSMPT
{
class Class_Potential_Origin;
namespace Wall
{
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
} // namespace Wall

} // namespace BSMPT

#endif // WALLTHICKNESSCOMMON_H
