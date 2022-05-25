// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#ifndef MINIMIZEGSL_H_
#define MINIMIZEGSL_H_

#include <cmath>
#include <gsl/gsl_vector_double.h> // for gsl_vector
#include <memory>
#include <vector>

namespace BSMPT
{
class Class_Potential_Origin;
namespace Minimizer
{

/**
 * @brief GSL_Tolerance Tolerance used in the GSL routines
 */
const double GSL_Tolerance = std::pow(10, -4);

/**
 * struct containing the required Parameters of the model for the gsl interface
 */
struct GSL_params
{
  std::shared_ptr<Class_Potential_Origin> modelPointer;
  double Temp;
  GSL_params(const std::shared_ptr<Class_Potential_Origin> &model,
             const double &temperature)
      : modelPointer{model}
      , Temp{temperature} {};
};

/**
 * Calculates the value of the effective potential at the vev v and temperature
 * p->Temp for the gsl interface
 */
double GSL_VEFF_gen_all(const gsl_vector *v, void *p);

/**
 * Calculates the next local minimum in the model from the point start
 * @returns The final status of the gsl minimization process.
 */
int GSL_Minimize_From_S_gen_all(struct GSL_params &p,
                                std::vector<double> &sol,
                                const std::vector<double> &start);

/**
 * Minimize the Potential from different random starting points and choose the
 * local minimum with the deepest potential value as the candidate for the
 * global minimum
 * @param model model reference
 * @param Temp Temperature at which to minimise the parameter point
 * @param seed seed used to find the random starting points for the local
 * optimisations
 * @param UseMultiThreading Decides if the algorithm should use multithreading
 * or not
 * @return first: vector with candidate for the global minimum, second: True if
 * a candidate for the global minimum is found and false otherwise
 */
std::pair<std::vector<double>, bool>
GSL_Minimize_gen_all(const Class_Potential_Origin &model,
                     const double &Temp,
                     const int &seed,
                     bool UseMultiThreading = true);

/**
 * Minimize the Potential from different random starting points and choose the
 * local minimum with the deepest potential value as the candidate for the
 * global minimum
 * @param model model reference
 * @param Temp Temperature at which to minimise the parameter point
 * @param seed seed used to find the random starting points for the local
 * optimisations
 * @param MaxSol numbers of local minima to find
 * @param UseMultiThreading Decides if the algorithm should use multithreading
 * or not
 * @return first: vector with candidate for the global minimum, second: True if
 * a candidate for the global minimum is found and false otherwise
 */
std::pair<std::vector<double>, bool>
GSL_Minimize_gen_all(const Class_Potential_Origin &model,
                     const double &Temp,
                     const int &seed,
                     const std::size_t &MaxSol,
                     bool UseMultiThreading = true);

/**
 * Minimize the Potential from different random starting points and choose the
 * local minimum with the deepest potential value as the candidate for the
 * global minimum
 * @param model model reference
 * @param Temp Temperature at which to minimise the parameter point
 * @param seed seed used to find the random starting points for the local
 * optimisations
 * @param saveAllMinima List of all local minima
 * @param MaxSol numbers of local minima to find
 * @param UseMultiThreading Decides if the algorithm should use multithreading
 * or not
 * @return first: vector with the solution, second: True if a candidate for the
 * global minimum is found and false otherwise
 */
std::pair<std::vector<double>, bool>
GSL_Minimize_gen_all(const Class_Potential_Origin &model,
                     const double &Temp,
                     const int &seed,
                     std::vector<std::vector<double>> &saveAllMinima,
                     const std::size_t &MaxSol,
                     bool UseMultiThreading = true);

} // namespace Minimizer
} // namespace BSMPT

#endif // MINIMIZEGSL_H_
