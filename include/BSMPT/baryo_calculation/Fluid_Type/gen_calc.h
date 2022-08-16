// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

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
