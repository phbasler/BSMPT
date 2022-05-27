// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef MINIMIZELIBCMAES_H
#define MINIMIZELIBCMAES_H

#include <memory>
#include <vector>
/**
 * @file
 * Minimization using the libcmaes algorithm
 */

namespace BSMPT
{
class Class_Potential_Origin;
namespace Minimizer
{
struct PointerContainerMinPlane;
namespace LibCMAES
{

/**
 * @brief The LibCMAESReturn struct which is returned by the routine
 */
struct LibCMAESReturn
{
  std::vector<double> result;
  int CMAESStatus;
};

/**
 * Calculating the global minimum with libcmaes in the 2HDM for the parameter
 * point par and counterterms parCT and write the solution in sol. The initial
 * guess is given in start.
 * @return the libcmaes run_status of the system
 */
LibCMAESReturn min_cmaes_gen_all(const Class_Potential_Origin &model,
                                 const double &Temp,
                                 const std::vector<double> &VevMinimum);

/**
 * Finds a candidate for the local minimum using the CMAES algorithm.
 * @param params PointerContainerMinPlane with the model information and
 * potential
 * @param Start starting point for the algorithm
 * @return A LibCMAESReturn struct with the CMAES run_status and the candidate
 * for the global minimum
 */
LibCMAESReturn
CMAES_Minimize_Plane_gen_all(const struct PointerContainerMinPlane &params,
                             const std::vector<double> &Start);

} // namespace LibCMAES
} // namespace Minimizer
} // namespace BSMPT

#endif // MINIMIZELIBCMAES_H
