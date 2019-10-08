/*
 * Minimizer.h
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

#ifndef MINIMIZER_H_
#define MINIMIZER_H_

#include "../models/IncludeAllModels.h"
#include "libcmaes/cmaes.h"
#include "libcmaes/cmaparameters.h"
#include "libcmaes/esoptimizer.h"
#include "libcmaes/cmastrategy.h"
#include "libcmaes/llogging.h"
#include <random>
#include <memory>

// Minimizer.cpp
void Minimize_gen_all(const std::shared_ptr<Class_Potential_Origin>& modelPointer,
		double Temp, std::vector<double>& sol, std::vector<double>& Check,const std::vector<double>& start, int WhichMinimizer=3);
void PTFinder_gen_all(const std::shared_ptr<Class_Potential_Origin>& modelPointer,
		double TempStart, double TempEnde, std::vector<double>& sol, int WhichMinimizer=3);

// Tree level minimization
void Minimize_gen_all_tree_level(int Model, const std::vector<double>& par, const std::vector<double>& parCT,
		std::vector<double>& sol, std::vector<double>& Check,const std::vector<double>& start, int WhichMinimizer=3);

// MinimizeGSL.cpp
double GSL_VEFF_gen_all(const gsl_vector *v, void *p);
int GSL_Minimize_From_S_gen_all(struct GSL_params& p, std::vector<double>& sol,const std::vector<double>& start);
bool GSL_Minimize_gen_all(const std::shared_ptr<Class_Potential_Origin>& modelPointer,
		double Temp, std::vector<double>& solV, int seed);
bool GSL_Minimize_gen_all(const std::shared_ptr<Class_Potential_Origin>& modelPointer,
		double Temp, std::vector<double>& solV, int seed,int MaxSol);
bool GSL_Minimize_gen_all(const std::shared_ptr<Class_Potential_Origin>& modelPointer,
		double Temp, std::vector<double>& solV, int seed, std::vector<std::vector<double>>& saveAllMinima, int MaxSol);


// Minfunc_gen.cpp
double min_cmaes_gen_all(const std::shared_ptr<Class_Potential_Origin>& modelPointer,
		double Temp, std::vector<double>& sol, const std::vector<double>& VevMinimum);




struct PointerContainer{
//  Class_Potential_Origin * modelPointer;
	std::shared_ptr<Class_Potential_Origin> modelPointer;
  double Temp;
};

#endif /* MINIMIZER_H_ */
