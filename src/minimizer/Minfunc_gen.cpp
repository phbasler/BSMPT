// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <algorithm>                           // for copy, max
#include <iostream>                            // for operator<<, endl, cout
#include <libcmaes/candidate.h>                // for Candidate
#include <libcmaes/cmaes.h>                    // for cmaes
#include <libcmaes/cmaparameters.h>            // for CMAParameters
#include <libcmaes/cmasolutions.h>             // for CMASolutions
#include <libcmaes/esoptimizer.h>              // for aCMAES
#include <libcmaes/esostrategy.h>              // for FitFunc
#include <libcmaes/genopheno.h>                // for GenoPheno
#include <libcmaes/noboundstrategy.h>          // for libcmaes
#include <memory>                              // for shared_ptr, __shared_...
#include <vector>                              // for vector

/**
 * @file
 * Finding a candidate for the global minimum using the CMA-ES algorithm
 */

namespace BSMPT
{
namespace Minimizer
{

}
} // namespace BSMPT
