/*
 * gen_calc.h
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

#ifndef GEN_CALC_H
#define GEN_CALC_H

/**
 * @file
 */

#include <BSMPT/baryo_calculation/transport_equations.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <boost/numeric/odeint.hpp>
#include <iostream>

#include <boost/any.hpp>

namespace BSMPT
{
namespace Baryo
{
/**
 * @brief set_up_nL_grid Calculation of the grid points needed for the cubic
 * spline.
 * @param n_step Number of grid points used for the cubic spline grid.
 * @param container Container with all needed parameters.
 * @param classpointer Class reference to the transport equation class
 * (top,bot,tau)
 * @return The pair with (z , nL(z) )
 */
std::pair<std::vector<double>, std::vector<double>>
set_up_nL_grid(std::size_t n_step,
               GSL_integration_mubl &container,
               boost::any const &classpointer);

} // namespace Baryo
} // namespace BSMPT

#endif
