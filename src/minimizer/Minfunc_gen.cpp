/*
 * Minfunc_gen.cpp
 *
 *  Created on: Jun 15, 2016
 *      Author: basler
 */

#include "Minimizer.h"

/**
 * @file
 * Finding a candidate for the global minimum using the CMA-ES algorithm
 */


double Tempcom;


PointerContainer PCon;
using namespace libcmaes;


/**
 *
 *
 * Interface for libcmaes to calculate the value of the effective Potential at vev = v and Temp = Tempcom
 */
FitFunc fsphere_gen_all = [](const double *v, const int N)
{

    double *p = (double *)v;



    int nVEVs = PCon.modelPointer->nVEV;



    std::vector<double> vIn,vMin;
    vMin.resize(nVEVs);
    for(int i=0;i<nVEVs;i++) vMin[i] = p[i];
    PCon.modelPointer->MinimizeOrderVEV(vMin,vIn);







    double res = PCon.modelPointer->VEff(vIn,Tempcom,0);
    return res;

};


/**
 * Calculating the global minimum with libcmaes in the 2HDM for the parameter point par and counterterms parCT
 * and write the solution in sol. The initial guess is given in start.
 * @return the libcmaes run_status of the system
 */
double min_cmaes_gen_all(int Model, const std::vector<double>& par, const std::vector<double>&  parCT,
		double Temp, std::vector<double>& sol, const std::vector<double>& Start){

    bool Debug = !C_TurnOffDebug;
    Debug = false;
    if(Debug) std::cout << "Start Debuging " << __func__ << std::endl;
//	double *p = (double *)par;
//    double *pCT =(double *)parCT;
    //double *Start = (double *)st;

    int nPar,nParCT;
//    ModelCom = Model;
    Tempcom = Temp;

//    Class_Potential_Origin * modelPointer;
//    Fchoose(modelPointer,Model);

    std::shared_ptr<Class_Potential_Origin> modelPointer = FChoose(Model);

    nPar = modelPointer->nPar;
    nParCT = modelPointer->nParCT;



    modelPointer->set_All(par,parCT);

    PCon.modelPointer = modelPointer;

    int dim = modelPointer->nVEV;
    double startval[dim];
    for(int i=0;i<dim;i++) startval[i] = Start.at(i);
    std::vector<double> x0(dim);


    double sigma;
    sigma = 5;

    for(int i=0;i<dim;i++) x0[i]=startval[i];





    if(Debug) std::cout << "Setup done \n";


    //sigma *= 0.5;//0.5;
    double ftol=1e-5;

    if(Debug)
    {
    	std::cout << "dim = " << dim << std::endl;
    	std::cout << "x0 = ";
    	for(int i=0;i<dim;i++) std::cout << x0[i] << "\t";
    	std::cout << std::endl;
    }

    //CMAParameters<> cmaparams(dim,&x0.front(),sigma);
    CMAParameters<> cmaparams(x0,sigma);
    if(Debug) std::cout << "Declaring cmaparams done " << std::endl;

    //cmaparams.set_mt_feval(true);
    cmaparams.set_algo(aCMAES);
    //cmaparams.set_elitism(1);
	//cmaparams.set_noisy();
    cmaparams.set_ftolerance(ftol);


    //CMASolutions cmasols = cmaes<>(fsphere,cmaparams,CMAStrategy<CovarianceUpdate>::_defaultPFunc, grad_fsphere);
    CMASolutions cmasols = cmaes<>(fsphere_gen_all,cmaparams);

    Candidate bcand = cmasols.best_candidate();

    std::vector<double> xsol = bcand.get_x();





    bool checkZero = false;// (Temp >= model.TMIN());
    bool MinZero = false;

    sol.clear();

    for(int i = 0;i<dim;i++)
        	{
        	    	sol.push_back(xsol.at(i));
        	}










//    delete modelPointer;

    return cmasols.run_status();


}


