/*
 * MinimizeGSL.cpp
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

/**
 * @file
 * Using the Nelder-Mead Simplex algorithm, implemented in gsl, to find multiple local minima of the model
 * and compare them to find a candidate for the global minimum.
 */


#include <BSMPT/minimizer/MinimizeGSL.h>

#include <ext/alloc_traits.h>                   // for __alloc_traits<>::val...
#include <gsl/gsl_errno.h>                      // for gsl_set_error_handler...
#include <gsl/gsl_multimin.h>                   // for gsl_multimin_fminimizer
#include <gsl/gsl_vector_double.h>              // for gsl_vector_get, gsl_v...
#include <stdio.h>                              // for printf
#include <time.h>                               // for time, NULL, std::size_t
#include <algorithm>                            // for copy, max
#include <iostream>                             // for operator<<, endl, bas...
#include <limits>                               // for numeric_limits, numer...
#include <memory>                               // for shared_ptr, __shared_...
#include <random>                               // for default_random_engine
#include <vector>                               // for vector
#include <BSMPT/models/ClassPotentialOrigin.h>  // for Class_Potential_Origin
#include <BSMPT/minimizer/Minimizer.h>

namespace BSMPT {
namespace Minimizer{



double GSL_VEFF_gen_all(const gsl_vector *v, void *p)
{

    struct GSL_params * params =  static_cast<GSL_params*>(p);

    std::vector<double> vMin;
    auto nVEVs = params->modelPointer->get_nVEV();

    for(size_t i=0;i<nVEVs;i++)
    {
        vMin.push_back(gsl_vector_get(v,i));
    }

    auto vIn=params->modelPointer->MinimizeOrderVEV(vMin);

    double res = params->modelPointer->VEff(vIn,params->Temp,0);

    return res;

}



int GSL_Minimize_From_S_gen_all(struct GSL_params& params, std::vector<double>& sol,const std::vector<double>& start)
{
    gsl_set_error_handler_off();

    //	struct GSL_params * params = (struct GSL_params *) p;
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;

    double ftol = GSL_Tolerance;
    std::size_t MaxIter = 600;

    std::size_t iter = 0;
    int status;
    double size;

    std::size_t dim = params.nVEV;

    /* Starting point */
    x = gsl_vector_alloc (dim);
    for(size_t k=0;k<dim;k++) gsl_vector_set(x,k,start.at(k));
    ss = gsl_vector_alloc (dim);
    gsl_vector_set_all (ss, 1.0);



    /* Initialize method and iterate */
    minex_func.n = dim;
    minex_func.f = &GSL_VEFF_gen_all;
    minex_func.params = &params;
    s = gsl_multimin_fminimizer_alloc (T, dim);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status) break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, ftol);

    }
    while (status == GSL_CONTINUE && iter < MaxIter);

    if(status == GSL_SUCCESS)
    {
        for(size_t k=0;k<dim;k++) sol.push_back(gsl_vector_get(s->x,k));
    }
    else{
        for(size_t k=0;k<dim;k++) sol.push_back(0);
    }

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);

    return status;

}



std::pair<std::vector<double>,bool> GSL_Minimize_gen_all(
        const std::shared_ptr<Class_Potential_Origin>& modelPointer,
        const double& Temp,
        const int& seed,
        const std::size_t& MaxSol){
    std::vector<std::vector<double>> saveAllMinima;
    auto result = GSL_Minimize_gen_all(modelPointer,Temp, seed ,saveAllMinima, MaxSol);
    return result;
}


std::pair<std::vector<double>,bool> GSL_Minimize_gen_all(
        const std::shared_ptr<Class_Potential_Origin>& modelPointer,
        const double& Temp,
        const int& seed){
    std::vector<std::vector<double>> saveAllMinima;
    std::size_t MaxSol = 20;
    auto result = GSL_Minimize_gen_all(modelPointer,Temp, seed ,saveAllMinima, MaxSol);
    return result;
}

std::pair<std::vector<double>,bool> GSL_Minimize_gen_all(
        const std::shared_ptr<Class_Potential_Origin>& modelPointer,
        const double& Temp,
        const int& seed ,
        std::vector<std::vector<double>>& saveAllMinima,
        const std::size_t& MaxSol)
{
    struct GSL_params params;
    params.Temp = Temp;
    params.nVEV = modelPointer->get_nVEV();
    params.modelPointer = modelPointer;

    std::size_t dim = modelPointer->get_nVEV();

    std::default_random_engine randGen(seed);
    double RNDMax= 500;
    std::size_t MaxTries = 600;
    std::size_t tries = 0;
    std::size_t numOfSol = 0;
    //	int MaxSol = 20;
    std::size_t nCol = dim+2;
    std::vector<double> start,sol,vPot;
    do{
        start.resize(dim);
        for(size_t i=0;i<dim;i++) start[i] = RNDMax*(-1+2*std::generate_canonical<double,std::numeric_limits<double>::digits>(randGen));


        auto status = GSL_Minimize_From_S_gen_all(params,sol,start);

        if(status == GSL_SUCCESS){
            vPot=modelPointer->MinimizeOrderVEV(sol);

            std::vector<double> row(nCol);
            for(size_t i=0;i<dim;i++) row.at(i) = sol.at(i);
            row.at(dim) = modelPointer->EWSBVEV(vPot);
            row.at(dim+1) = modelPointer->VEff(vPot,Temp,0);

            saveAllMinima.push_back(row);
            numOfSol ++;
            if(numOfSol == MaxSol) break;
        }


        start.clear();
        sol.clear();
        tries++;
    }while(tries <= MaxTries);

    if(numOfSol == 0)
    {
//        std::cerr << "No solutions found during the GSL minimization at T = " << Temp << " GeV " << "\n";
        return std::make_pair(std::vector<double>{}, false);
    }

    std::size_t minIndex = 0;
    double VMin = saveAllMinima[0][dim+1];
    for(size_t k=1;k<numOfSol;k++)
    {
        if(saveAllMinima[k][dim+1] <= VMin)
        {
            VMin = saveAllMinima[k][dim+1];
            minIndex = k;
        }
    }

    std::vector<double> solV;
    for(size_t k=0;k<dim;k++) solV.push_back(saveAllMinima[minIndex][k]);
    return std::make_pair(solV,true);

}





}
}
