/*
 * WallThicknessLib.h
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

#ifndef WALLTHICKNESS_H
#define WALLTHICKNESS_H

/**
 * @file
 */

#include <memory>
#include <vector>

#include <BSMPT/minimizer/Minimizer.h>

namespace BSMPT
{
class Class_Potential_Origin;
namespace Wall
{

/**
 * @brief Calculate the Wall thickness using the minimization along the normal
 * plane
 * @param modelPointer_input Pointer to parameter point
 * @param Temp temperature at which the wall thickness is to be calculated
 * @param vcritical the electroweak VEV in the broken minimum
 * @param vevsymmetric the electroweak VEV in the symmetric minimum
 * @return The thickness of the wall in 1/GeV
 */
double
calculate_wall_thickness_plane(
    const std::shared_ptr<Class_Potential_Origin> &modelPointer_input,
    const double &Temp,
    const std::vector<double> &vcritical,
    const std::vector<double> &vevsymmetric,
    const int &WhichMinimizer = Minimizer::WhichMinimizerDefault);

/**
 * @brief calculate_wall_thickness_1D calculates calculates the wall thickness
 * across the straight line from the symmetrical minimum to the broken minimum
 * @param modelPointer Pointer to parameter point
 * @param Temp temperature at which the wall thickness is to be calculated
 * @param vcritical the electroweak VEV in the broken minimum
 * @param vevsymmetric the electroweak VEV in the symmetric minimum
 * @return
 */
double
calculate_wall_thickness_1D(
    const std::shared_ptr<Class_Potential_Origin> &modelPointer,
    const double &Temp,
    const std::vector<double> &vcritical,
    const std::vector<double> &vevsymmetric);
/**
 * @brief GSL_Find_Maximum_line finds the maximum of the function along a line
 * using GSL for the maximisation
 * @param modelPointer Pointer to parameter point
 * @param Temp temperature at which the wall thickness is to be calculated
 * @param vcritical the electroweak VEV in the broken minimum
 * @param vevsymmetric the electroweak VEV in the symmetric minimum
 * @param solV vector to store the solution
 * @return true if successfull, false if not
 */
bool
GSL_Find_Maximum_line(
    const std::shared_ptr<Class_Potential_Origin> &modelPointer,
    const double &Temp,
    const std::vector<double> &vcritical,
    const std::vector<double> &vevsymmetric,
    std::vector<double> &solV);

/**
 * @brief GSL_Maximize_From_S_gen_line finds the next local maximum from the
 * initial guess
 * @param params params struct containing the parameters needed for the
 * maximisation
 * @param solution List of solutions which will be expanded with the new
 * solution
 * @param initial_guess starting point for the local optimisation
 * @return true if successfull, false if not
 */
bool
GSL_Maximize_From_S_gen_line(struct GSL_params &params,
                             std::vector<std::vector<double>> &solution,
                             double initial_guess);

/**
 * @brief GSL_VEFF_gen_all_maximum_line Maximization along the line
 * @param t parametrizes the VEV along (1-t)*vevsymmetric + t*vevcritical
 * @param p pointer to the params struct storing the parameters
 * @return -1*VEff at the parameter point
 */
double
GSL_VEFF_gen_all_maximum_line(
    double t,
    void *p); // 1D Maximization along the line from 0 to omega

} // namespace Wall
} // namespace BSMPT

#endif /* WALLTHICKNESS_H */
