
/*
 * Minimizer.cpp
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


#include <vector>
#include <memory>
#include <math.h>                               // for abs, log10
#include <time.h>                               // for time, NULL
#include <algorithm>                            // for copy, max
#include <iostream>                             // for operator<<, cout, endl
#include <BSMPT/models/ClassPotentialOrigin.h>  // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>      // for FChoose
#include <BSMPT/minimizer/Minimizer.h>

#include <BSMPT/config.h>

#include <BSMPT/minimizer/MinimizeGSL.h>

#ifdef CMAES_FOUND
#include <BSMPT/minimizer/LibCMAES/MinimizeLibCMAES.h>
#endif

#ifdef NLopt_FOUND
#include <BSMPT/minimizer/LibNLOPT/MinimizeNLOPT.h>
#endif


/**
 * @file
 * @param WhichMinimizer = 1/2 for single minimizer, sum of those for multiple
 *
 * WhichMinimizer = 1 -> Use CMA-ES Minimizer
 *
 * WhichMinimizer = 2 -> Use the local Minimization of GSL
 *
 *
 *
 *
 *
 *
 * @param Model Which model is considered?
 * @param par Array with parameter point used for set_gen
 * @param parCT Array with Counterterm values used for set_CT_pot_par
 * @param TempStart Starting value of the Temperature for the bisection Method
 * @param TempEnde End value of the Temperature for the bisection Method
 * @param sol std::vector<double> in which the solution of the PTFinder will be stored in the format (T_C,v_C,-+1, single vevs) where +1 means sucessfull find of v_c/T_c > 1 and -1 hit one of the abortion criteria
 *
 */


namespace BSMPT {
namespace Minimizer {

/**
 * @brief Minimization of the Model
* Minimizes the given Potential with parameters par and CT-parameters parCT at a given Temperature Temp and writes the solution in the std::vector sol.
* The Minimization Debugging Options are written in the std::vector Check.
* The std::vector Start gives the start value for the CMA-ES Minimization.
*/
std::vector<double> Minimize_gen_all(
        const std::shared_ptr<Class_Potential_Origin>& modelPointer,
        const double& Temp,
        std::vector<double>& Check,
        const std::vector<double>& start,
        const int& WhichMinimizer)
{
    std::vector<double> PotValues;
    std::vector<std::vector<double>> Minima;

    bool UseCMAES = false;
    bool UseGSLLocal = false;
    bool UseNLOPT = false;


    int PGSL,PCMAES,PNLOPT;
    int WMx = WhichMinimizer;
    PCMAES = WMx%2;
    WMx = WMx/2;
    PGSL = WMx%2;
    WMx = WMx/2;
    PNLOPT = WMx%2;

    UseNLOPT = (PNLOPT == 1);
    UseCMAES = (PCMAES == 1);
    UseGSLLocal = (PGSL == 1);

#ifndef CMAES_FOUND
    UseCMAES = false;
    (void) start;
#endif

#ifndef NLopt_FOUND
    UseNLOPT = false;
#endif


    bool CheckZero = true; // Check if zero is the global minimum explicitly
    if(CheckZero){
        PotValues.push_back(modelPointer->VEff(std::vector<double>(modelPointer->get_NHiggs(),0),Temp));
        Minima.push_back(std::vector<double>(modelPointer->get_nVEV(),0));
    }

    if(modelPointer->get_nVEV() <=2)
      {
				UseCMAES = false;
				UseGSLLocal = true;
      }

    std::vector<double> solGSLMin,solGSLMinPot;

    bool gslMinSuc = false;
    if(UseGSLLocal) {
        if(UseCMAES or UseNLOPT) {
            std::tie(solGSLMin,gslMinSuc) = GSL_Minimize_gen_all(modelPointer, Temp, 5); // If additionally CMAES is minimising GSL does not need as much solutions
        }
        else {
            std::size_t MaxSol = 50;
            std::tie(solGSLMin,gslMinSuc) = GSL_Minimize_gen_all(modelPointer, Temp, 5, MaxSol);
        }

        if(gslMinSuc){
            solGSLMinPot=modelPointer->MinimizeOrderVEV(solGSLMin);
            PotValues.push_back(modelPointer->VEff(solGSLMinPot,Temp));
            Minima.push_back(solGSLMin);
        }


    }
#ifdef CMAES_FOUND
    if(UseCMAES) {
        auto LibCMAES = LibCMAES::min_cmaes_gen_all(modelPointer,Temp,start);
        auto errC = LibCMAES.CMAESStatus;
        auto solCMAES = LibCMAES.result;
        auto solCMAESPotIn=modelPointer->MinimizeOrderVEV(solCMAES);
        PotValues.push_back(modelPointer->VEff(solCMAESPotIn,Temp));
        Minima.push_back(solCMAES);
        Check.push_back(errC);
    }
#endif

#ifdef NLopt_FOUND
    if(UseNLOPT)
    {
//        std::cout<<"NLO opt called"<<std::endl;
        auto NLOPTResult = LibNLOPT::MinimizeUsingNLOPT(modelPointer,Temp);
        PotValues.push_back(NLOPTResult.PotVal);
        Minima.push_back(NLOPTResult.Minimum);
    }
#endif


    std::size_t minIndex=0;
    for(size_t i=1;i<PotValues.size();i++){
        if(PotValues.at(i) < PotValues.at(minIndex)) minIndex = i;
    }

    auto sol = Minima.at(minIndex);
    auto EWVEV = modelPointer->EWSBVEV(modelPointer->MinimizeOrderVEV(sol));
    if(EWVEV <= 0.5) sol = std::vector<double>(modelPointer->get_nVEV(),0);

    solGSLMin.clear();
    if(UseGSLLocal and  gslMinSuc) Check.push_back(1);
    else Check.push_back(-1);

    return sol;

}




EWPTReturnType PTFinder_gen_all(
        const std::shared_ptr<Class_Potential_Origin>& modelPointer,
        const double& TempStart,
        const double& TempEnd,
        const int& WhichMinimizer)
{

    EWPTReturnType result;

    std::size_t dim = modelPointer->get_nVEV();

    double vStart,vEnd,vMitte;
    std::vector<double> solStart,solMitte,solEnd;
    std::vector<double> solStartPot,solMittePot,solEndPot;
    std::vector<double> checkStart,checkMitte,checkEnde;
    std::vector<double> startStart,startMitte,startEnde;
    double Distance=1e-2;
    double TA = TempStart;
    double TM;
    double TE = TempEnd;

    std::vector<double> pSol(dim+1);

    for(size_t k=0;k<dim;k++) startEnde.push_back(modelPointer->get_vevTreeMin(k));
    solEnd = Minimize_gen_all(modelPointer,TempEnd,checkEnde,startEnde,WhichMinimizer);
    solEndPot=modelPointer->MinimizeOrderVEV(solEnd);
    vEnd = modelPointer->EWSBVEV(solEndPot);



    if(vEnd > C_threshold)
    {
        result.Tc = TempEnd;
        result.vc = vEnd;
        result.StatusFlag = MinimizerStatus::NOTVANISHINGATFINALTEMP;
        result.EWMinimum = solEnd;
        return result;
    }

    for(size_t k=0;k<dim;k++) startStart.push_back(modelPointer->get_vevTreeMin(k));
    solStart=Minimize_gen_all(modelPointer,TempStart,checkStart,startStart,WhichMinimizer);
    solStartPot=modelPointer->MinimizeOrderVEV(solStart);
    vStart = modelPointer->EWSBVEV(solStartPot);



    bool SurviveNLO = modelPointer->CheckNLOVEV(solStart);

    if(not SurviveNLO and TempStart == 0 )
    {
        result.Tc = TempEnd;
        result.vc = vStart;
        result.StatusFlag = MinimizerStatus::NOTNLOSTABLE;
        result.EWMinimum = solStart;
        return result;

    }

    if(vStart <= C_threshold or vStart >= 255.0)
    {
        result.Tc = TempEnd;
        result.vc = 0;
        result.StatusFlag = MinimizerStatus::NLOVEVZEROORINF;
        result.EWMinimum = std::vector<double>(dim,0);
        return result;
    }

    for(size_t k=0;k<dim;k++) startMitte.push_back(modelPointer->get_vevTreeMin(k));
    do{
	    TM = 0.5*(TA+TE);
	    solMitte.clear();
	    checkMitte.clear();
	    solMittePot.clear();
        solMitte=Minimize_gen_all(modelPointer,TM,checkMitte,startMitte,WhichMinimizer);
        solMittePot=modelPointer->MinimizeOrderVEV(solMitte);
	    vMitte = modelPointer->EWSBVEV(solMittePot);


	     if(vMitte >= 255.0)
	    {
             result.Tc = TM;
             result.vc = vMitte;
             result.StatusFlag = MinimizerStatus::NUMERICALLYUNSTABLE;
             result.EWMinimum = solMitte;
             return result;
	    }
	    else if(vMitte != 0 and vMitte/TM < C_PT)
	    {
             result.Tc = TM;
             result.vc = vMitte;
             result.StatusFlag = MinimizerStatus::BELOWTHRESHOLD;
             result.EWMinimum = solMitte;
             return result;

	    }
	    if(vMitte >= Distance)
	    {
		    TA = TM;
		    startMitte.clear();
            for(size_t k=0;k<dim;k++) pSol[k+1] = solMitte.at(k);
		    pSol[0] = vMitte;
            for(size_t k=0;k<dim;k++) startMitte.push_back(pSol[k+1]);
		    vStart = vMitte;
	    }
	    else{
		    TE = TM;
	    }

    }while(std::abs(TE-TA) > Distance);

    result.Tc = TA;
    result.vc = pSol[0];
    result.StatusFlag = MinimizerStatus::SUCCESS;
    std::vector<double> SolMin;
    result.EWMinimum.clear();
    for(size_t k=1;k<=dim;k++)  result.EWMinimum.push_back(pSol[k]);
    return result;
}


std::vector<double> Minimize_gen_all_tree_level(
        const ModelID::ModelIDs& Model,
        const std::vector<double>& par,
        const std::vector<double>& parCT,
        std::vector<double>& Check,
        const std::vector<double>& start,
        int WhichMinimizer)
{
    std::shared_ptr<Class_Potential_Origin> modelPointer = ModelID::FChoose(Model);
    modelPointer->set_All(par,parCT);
	modelPointer->SetUseTreeLevel(true);
    auto sol=Minimize_gen_all(modelPointer,0,Check,start,WhichMinimizer);
    modelPointer->SetUseTreeLevel(false);
    return sol;
}

}
}
