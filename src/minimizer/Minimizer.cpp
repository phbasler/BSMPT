
/*
 * Minimizer.cpp
 *
 *  Created on: Jun 15, 2016
 *      Author: basler
 */


#include "Minimizer.h"

//#include "EqSolGenericInput.cpp"
//#include "Minfunc_gen.cpp"

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




/**
 * @brief Minimization of the Model
* Minimizes the given Potential with parameters par and CT-parameters parCT at a given Temperature Temp and writes the solution in the std::vector sol.
* The Minimization Debugging Options are written in the std::vector Check.
* The std::vector Start gives the start value for the CMA-ES Minimization.
*/
void Minimize_gen_all(int Model, const std::vector<double>& par, const std::vector<double>&  parCT, double Temp,
		std::vector<double>& sol, std::vector<double>& Check,const std::vector<double>& start, int WhichMinimizer)
{
    bool Debug = !C_TurnOffDebug;
    Debug = false;

    if(Debug) std::cout << "Debug turned on in " << __func__ << std::endl;

    bool UseCMAES = false;
    bool UseGSLLocal = false;


    int PGSL,PCMAES;
    int WMx = WhichMinimizer;
    PCMAES = WMx%2;
    WMx = WMx/2;
    PGSL = WMx%2;


    if(Debug) std::cout << WhichMinimizer << "\t" << PGSL << "\t" << PCMAES << std::endl;


    UseCMAES = (PCMAES == 1);
    UseGSLLocal = (PGSL == 1);


//    Class_Potential_Origin * modelPointer;
//    Fchoose(modelPointer,Model);
    std::unique_ptr<Class_Potential_Origin> modelPointer = FChoose(Model);
    int npar = modelPointer->nPar;
//    double *p = (double *)par;
//    double *pct = (double *)parCT;
    modelPointer->set_All(par,parCT);
    int dim = modelPointer->nVEV;
    int NHiggs = modelPointer->NHiggs;

    if(Debug) std::cout << "dim = " << dim << std::endl;

    double timeStart=0,timeEnd=0;

    if(dim <3)
      {
	UseCMAES = false;
	UseGSLLocal = true;
      }



    if(Debug) std::cout << "\n\n\nBeginne Minimization at T = " << Temp << std::endl;
    std::vector<double> solGSLMin,solGSLMinPot;
    bool gslMinSuc = false;
    if(UseGSLLocal) {
		if(Debug) std::cout << "Start with GSL Minimization " << std::endl;
		timeStart = time(NULL);
		gslMinSuc = GSL_Minimize_gen_all(Model, par, parCT, Temp,solGSLMin, 5);
		timeEnd = time(NULL);
		if(Debug) std::cout << "GSL Runtime = " << timeEnd-timeStart << "seconds" <<  std::endl;
		modelPointer->MinimizeOrderVEV(solGSLMin,solGSLMinPot);
		if(Debug) std::cout << "End of GSL Minimization " << std::endl;
    }
    if(Debug)
    	{
			std::cout << "GSL Minimizer done and ";
			if(gslMinSuc) {
				std::cout << "was sucessful \n";

				std::cout << "GSL Solution at \n";
				for(int i=0;i<dim;i++) std::cout << solGSLMin.at(i) << "\t";
				std::cout << std::endl;
				std::cout << "V_{GSL} = " << solGSLMin.at(dim+1) << std::endl;
			}
			else std::cout << "was unsucessful \n";
    	}











    std::vector<double> solCMAESt1,solCMAES,solCMAESPotIn;

    double stA[dim];
    if(Debug) std::cout << "Setze initial guess for CMAES " << std::endl;
    for(int k=0;k<dim;k++) stA[k] = start.at(k);
    if(Debug) std::cout << "Starting at : ";
    if(Debug) for(int k=0;k<dim;k++) std::cout << stA[k] << "\t";
    if(Debug) std::cout << std::endl;
    int errC = 0;
    if(UseCMAES) {
		if(Debug) std::cout << "Start CMAES " << std::endl;
		errC = min_cmaes_gen_all(Model,par,parCT,Temp,solCMAES,start);
		if(Debug) std::cout << "Set Solution " << std::endl;
		modelPointer->MinimizeOrderVEV(solCMAES,solCMAESPotIn);
		if(Debug ) std::cout << "Done " << std::endl;
		if(Debug) {
			std::cout << "CMAES Solution at \n";
			for(int i=0;i<dim;i++) std::cout << solCMAES.at(i) << "\t";
			std::cout << std::endl;
			std::cout << "V_{CMAES} = " << modelPointer->VEff(solCMAESPotIn,Temp,0) << std::endl;
		}
    }






    bool GSLCheck=false;




    std::vector<double> solTmp;

    bool MinCMAES=true;

    int MinMethod=-1;
    if(Debug) std::cout << "Begin comparission" << std::endl;

    bool cond1 = gslMinSuc and (not UseCMAES);
    bool cond2 = UseCMAES and (not gslMinSuc);
    bool cond3 = UseCMAES and gslMinSuc;

    if(cond1) MinMethod = 0;

    else if(cond2) MinMethod = 1;


    double ValGSL,ValCMAES;
    if(gslMinSuc) ValGSL = solGSLMin.at(dim+1);
    if(UseCMAES) ValCMAES = modelPointer->VEff(solCMAESPotIn,Temp,0);



    double mintmp=0;

    if(cond3){
    	MinMethod = 0;
    	if(ValCMAES < ValGSL) MinMethod = 1;
    }






    if(MinMethod == 0)
    {
	    for(int k=0;k<dim;k++) solTmp.push_back(solGSLMin.at(k));
    }
    else if(MinMethod == 1)
    {
	    for(int k=0;k<dim;k++) solTmp.push_back(solCMAES.at(k));
    }


    if(Debug) std::cout << "Global Minimum given by (0 GSL Min : 1 CMAES ) " << MinMethod << std::endl;
    bool CheckZero = true;
    double zeroDiffval;
    bool MinZero = false;
    std::vector<double> solTmpPot;
    modelPointer->MinimizeOrderVEV(solTmp,solTmpPot);
    if(CheckZero)
    {
      std::vector<double> vecZero;
      for(int i=0;i<NHiggs;i++) vecZero.push_back(0);
      double VZero = modelPointer->VEff(vecZero,Temp,0);

      double VMin = modelPointer->VEff(solTmpPot,Temp,0);
      double DifferenceBoarder = 1;
      double OrderVMin = std::log10(std::abs(VMin));
      MinZero = (VZero <= VMin);
      if(!MinZero) MinZero = (VZero <= 0 && VMin <= 0 && VZero - VMin <= std::pow(10,-OrderVMin)*VMin);
      MinZero = (std::abs(modelPointer->VEff(vecZero,Temp,0) - modelPointer->VEff(solTmpPot,Temp,0)) <= DifferenceBoarder);
      if(Debug) std::cout << "CheckZero : =" <<  modelPointer->VEff(vecZero,Temp,0)  - modelPointer->VEff(solTmpPot,Temp,0) << "\n";
      if(Debug) std::cout << "VZero = " << VZero << "\t VMin = " << VMin << "\n";
    }

    double vevsolTMp = 0;
    vevsolTMp = modelPointer->EWSBVEV(solTmpPot);
    if(MinZero or vevsolTMp <= 0.5)
    {
	    for(int k=0;k<dim;k++) sol.push_back(modelPointer->ModifiedVEVVectorDim[k]);
    }
    else{
	    for(int k=0;k<dim;k++) sol.push_back(solTmp.at(k));
    }


    if(Debug)
      {
	for(int i=0;i<NHiggs;i++) std::cout << modelPointer->vevTree[i] << "\t";
      }






//    delete modelPointer;
    solCMAES.clear();
    solTmp.clear();
    Check.push_back(errC);
    solGSLMin.clear();
    if(UseGSLLocal and  gslMinSuc) Check.push_back(1);
        else Check.push_back(-1);
    if(MinZero) Check.push_back(1);
    else Check.push_back(0);
    return ;

}



/**
* Uses a bisection method between the Temperature TempStart and TempEnde to find the phase transition in the Model and writes the solution in the vector sol.
*/
void PTFinder_gen_all(int Model, const std::vector<double>& par, const std::vector<double>&  parCT, double TempStart, double TempEnde, std::vector<double>& sol, int WhichMinimizer)
{
	bool Debug = !C_TurnOffDebug;


    if(Debug) std::cout << "Starte PTFinder from Temperature " << TempStart << " until " << TempEnde <<  "\n";

//    Class_Potential_Origin * modelPointer;
//    Fchoose(modelPointer,Model);
    std::unique_ptr<Class_Potential_Origin> modelPointer = FChoose(Model);
    int npar = modelPointer->nPar;
//    double *p = (double *)par;
//    double *pct = (double *)parCT;
    modelPointer->set_All(par,parCT);
    int dim = modelPointer->nVEV;
    int NHiggs = modelPointer->NHiggs;



    double vStart,vEnd,vMitte;
    std::vector<double> solStart,solMitte,solEnd;
    std::vector<double> solStartPot,solMittePot,solEndPot;
    std::vector<double> checkStart,checkMitte,checkEnde;
    std::vector<double> startStart,startMitte,startEnde;
    double Distance=1e-2;
    double TA = TempStart;
    double TM;
    double TE = TempEnde;

    double pSol[dim+1];






    for(int k=0;k<dim;k++) startEnde.push_back(modelPointer->vevTreeMin.at(k));
    Minimize_gen_all(Model,par,parCT,TempEnde,solEnd,checkEnde,startEnde,WhichMinimizer);
    modelPointer->MinimizeOrderVEV(solEnd,solEndPot);
    vEnd = modelPointer->EWSBVEV(solEndPot);



    if(vEnd > C_threshold)
    {
	    if(Debug) std::cout << "v != 0 at  T = " << TempEnde << "\n";
	    sol.push_back(TempEnde);
	    sol.push_back(vEnd);
	    sol.push_back(-1);
	    for(int k=0;k<dim;k++) sol.push_back(solEnd.at(k));

	    if(Debug) std::cout << "vev is given through : ";
	    if(Debug) for(int k=0;k<dim;k++) std::cout << solEnd.at(k) << "\t";
	    if(Debug) std::cout << "\n";

//	    delete modelPointer;
	    return ;
    }

    for(int k=0;k<dim;k++) startStart.push_back(modelPointer->vevTreeMin.at(k));
    Minimize_gen_all(Model,par,parCT,TempStart,solStart,checkStart,startStart,WhichMinimizer);
    modelPointer->MinimizeOrderVEV(solStart,solStartPot);
    vStart = modelPointer->EWSBVEV(solStartPot);

    if(vStart <= C_threshold or vStart >= 255.0)
    {
    if(Debug) std::cout << "v(T=0) = 0" << std::endl;
	    sol.push_back(TempEnde);
	    sol.push_back(0);
	    sol.push_back(-1);
	    for(int k=0;k<dim;k++) sol.push_back(0);
//	    delete modelPointer;
	    return ;
    }

    bool SurviveNLO = false;
	SurviveNLO = modelPointer->CheckNLOVEV(solStart);

    if(!SurviveNLO and TempStart == 0 )
    {
	    if(Debug) std::cout << "v_Tree is not the global Minimum at T=0 \n";
	    sol.push_back(TempEnde);
	    sol.push_back(0);
	    sol.push_back(-1);
	    for(int k=0;k<dim;k++) sol.push_back(0);
//	    delete modelPointer;
	    return ;

    }

    for(int k=0;k<dim;k++) startMitte.push_back(modelPointer->vevTreeMin.at(k));
    do{
	    TM = 0.5*(TA+TE);
	    solMitte.clear();
	    checkMitte.clear();
	    solMittePot.clear();
	    Minimize_gen_all(Model,par,parCT,TM,solMitte,checkMitte,startMitte,WhichMinimizer);
	    modelPointer->MinimizeOrderVEV(solMitte,solMittePot);
	    vMitte = modelPointer->EWSBVEV(solMittePot);


	     if(vMitte >= 255.0)
	    {
	      sol.push_back(TM);
	      sol.push_back(vMitte);
	      sol.push_back(-1);
	      for(int k=0;k<dim;k++) sol.push_back(solMitte.at(k));
//	      delete modelPointer;
	      return;
	    }
	    else if(vMitte != 0 and vMitte/TM < C_PT)
	    {
		    sol.push_back(TM);
		    sol.push_back(vMitte);
		    sol.push_back(-1);
		    for(int k=0;k<dim;k++) sol.push_back(solMitte.at(k));
//		    delete modelPointer;
		    return;
	    }
	    if(vMitte >= Distance)
	    {
		    TA = TM;
		    startMitte.clear();
		    for(int k=0;k<dim;k++) pSol[k+1] = solMitte.at(k);
		    pSol[0] = vMitte;
		    for(int k=0;k<dim;k++) startMitte.push_back(pSol[k+1]);
//			std::cout << "Array : ";
//			for(int k=0;k<=3;k++) std::cout << pSol[k] << "\t";
//			std::cout << "\n";
		    vStart = vMitte;
	    }
	    else{
		    TE = TM;
	    }

    }while(std::abs(TE-TA) > Distance);

    sol.push_back(TA);
    sol.push_back(pSol[0]);
    sol.push_back(1);
    for(int k=1;k<=dim;k++) sol.push_back(pSol[k]);
//    delete modelPointer;
    return ;

}



