/*
 * Minimizer.h
 *
 *  Created on: Mar 8, 2017
 *      Author: basler
 */

#ifndef MINIMIZER_H_
#define MINIMIZER_H_

#include "../models/IncludeAllModels.h"
#include "cmaes.h"
#include <cmaparameters.h>
#include "esoptimizer.h"
#include "cmastrategy.h"
#include "llogging.h"
#include <random>
#include <memory>

// Minimizer.cpp
void Minimize_gen_all(int Model, const std::vector<double>& par, const std::vector<double>&  parCT,
		double Temp, std::vector<double>& sol, std::vector<double>& Check,const std::vector<double>& start, int WhichMinimizer=3);
void PTFinder_gen_all(int Model, const std::vector<double>& par, const std::vector<double>&  parCT,
		double TempStart, double TempEnde, std::vector<double>& sol, int WhichMinimizer=3);

// MinimizeGSL.cpp
double GSL_VEFF_gen_all(const gsl_vector *v, void *p);
int GSL_Minimize_From_S_gen_all(const std::vector<double>& p, std::vector<double>& sol,const std::vector<double>& start);
bool GSL_Minimize_gen_all(int Model, const std::vector<double>& par, const std::vector<double>&  parCT,
		double Temp, std::vector<double>& solV, int seed);
bool GSL_Minimize_gen_all(int Model, const std::vector<double>& par, const std::vector<double>&  parCT,
		double Temp, std::vector<double>& solV, int seed,int MaxSol);
bool GSL_Minimize_gen_all(int Model, const std::vector<double>& par, const std::vector<double>&  parCT,
		double Temp, std::vector<double>& solV, int seed, std::vector<std::vector<double>>& saveAllMinima, int MaxSol);


// Minfunc_gen.cpp
double min_cmaes_gen_all(int Model, const std::vector<double>& par, const std::vector<double>&  parCT,
		double Temp, std::vector<double>& sol, const std::vector<double>& Start);




struct PointerContainer{
//  Class_Potential_Origin * modelPointer;
	std::shared_ptr<Class_Potential_Origin> modelPointer;
  double Temp;
};

#endif /* MINIMIZER_H_ */

