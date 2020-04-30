/*
 * Minfunc_gen.cpp
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

#include <libcmaes/candidate.h>                 // for Candidate
#include <libcmaes/cmaes.h>                     // for cmaes
#include <libcmaes/cmaparameters.h>             // for CMAParameters
#include <libcmaes/cmasolutions.h>              // for CMASolutions
#include <libcmaes/esoptimizer.h>               // for aCMAES
#include <libcmaes/esostrategy.h>               // for FitFunc
#include <libcmaes/genopheno.h>                 // for GenoPheno
#include <libcmaes/noboundstrategy.h>           // for libcmaes
#include <algorithm>                            // for copy, max
#include <iostream>                             // for operator<<, endl, cout
#include <memory>                               // for shared_ptr, __shared_...
#include <vector>                               // for vector
#include <BSMPT/models/ClassPotentialOrigin.h>  // for Class_Potential_Origin
#include <BSMPT/minimizer/Minimizer.h>

/**
 * @file
 * Finding a candidate for the global minimum using the CMA-ES algorithm
 */

namespace BSMPT {
namespace Minimizer{






}
}
