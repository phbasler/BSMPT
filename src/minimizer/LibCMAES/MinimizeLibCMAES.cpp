/*
 * bot_source.h
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

/**
 * @file
 */

#include <BSMPT/minimizer/LibCMAES/MinimizeLibCMAES.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/minimizer/MinimizePlane.h>

#include <libcmaes/candidate.h>                 // for Candidate
#include <libcmaes/cmaes.h>                     // for cmaes
#include <libcmaes/cmaparameters.h>             // for CMAParameters
#include <libcmaes/cmasolutions.h>              // for CMASolutions
#include <libcmaes/esoptimizer.h>               // for aCMAES
#include <libcmaes/esostrategy.h>               // for FitFunc
#include <libcmaes/genopheno.h>                 // for GenoPheno
#include <libcmaes/noboundstrategy.h>           // for libcmaes

namespace BSMPT {
namespace Minimizer{
PointerContainerMinPlane PCon_Plane;
namespace LibCMAES {



PointerContainer PCon;
using namespace libcmaes;


/**
 * @brief Interface for libcmaes to calculate the value of the effective Potential at vev = v and Temp = PCon_Planes.Temp
 */
FitFunc CMAES_VEFF_Minimize_Plane = [](const double *v, const int N)
{

//    double *p = static_cast<double*>(v);//(double *)v;
    (void) N;

    size_t nVEVs = PCon_Plane.modelPointer->get_nVEV();

    std::vector<double> vIn,vMinTilde;
    vMinTilde.resize(nVEVs-1);
    for(size_t i=0;i<nVEVs-1;i++) vMinTilde[i] = v[i];
    auto vMin = TransformCoordinates(vMinTilde,PCon_Plane);
    vIn=PCon_Plane.modelPointer->MinimizeOrderVEV(vMin);


    double res = PCon_Plane.modelPointer->VEff(vIn,PCon_Plane.Temp,0);
    return res;

};


/**
 *
 *
 * Interface for libcmaes to calculate the value of the effective Potential at vev = v and Temp = Tempcom
 */
FitFunc fsphere_gen_all = [](const double *v, const int N)
{
    double *p = (double *)v;
    (void) N;

    int nVEVs = PCon.modelPointer->get_nVEV();

    std::vector<double> vIn,vMin;
    vMin.resize(nVEVs);
    for(int i=0;i<nVEVs;i++) vMin[i] = p[i];
    vIn=PCon.modelPointer->MinimizeOrderVEV(vMin);

    double res = PCon.modelPointer->VEff(vIn,PCon.Temp,0);
    return res;

};



LibCMAESReturn min_cmaes_gen_all(
        const std::shared_ptr<Class_Potential_Origin>& modelPointer,
        const double& Temp,
        const std::vector<double>& Start){

    PCon.Temp = Temp;



    PCon.modelPointer = modelPointer;

    int dim = modelPointer->get_nVEV();

    std::vector<double> x0(dim);


    double sigma;
    sigma = 5;

    for(int i=0;i<dim;i++) x0[i]=Start.at(i);


    //sigma *= 0.5;//0.5;
    double ftol=1e-5;


    //CMAParameters<> cmaparams(dim,&x0.front(),sigma);
    CMAParameters<> cmaparams(x0,sigma);


    //cmaparams.set_mt_feval(true);
    cmaparams.set_algo(aCMAES);
    //cmaparams.set_elitism(1);
    //cmaparams.set_noisy();
    cmaparams.set_ftolerance(ftol);


    //CMASolutions cmasols = cmaes<>(fsphere,cmaparams,CMAStrategy<CovarianceUpdate>::_defaultPFunc, grad_fsphere);
    CMASolutions cmasols = cmaes<>(fsphere_gen_all,cmaparams);

    Candidate bcand = cmasols.best_candidate();

    std::vector<double> xsol = bcand.get_x();


    std::vector<double> sol;

    for(int i = 0;i<dim;i++)
    {
        sol.push_back(xsol.at(i));
    }

    LibCMAESReturn res;
    res.CMAESStatus = cmasols.run_status();
    res.result = sol;

    return res;


}


LibCMAESReturn CMAES_Minimize_Plane_gen_all(
        const struct PointerContainerMinPlane& params,
        const std::vector<double>& Start){


    PCon_Plane = params;

    std::shared_ptr<Class_Potential_Origin> modelPointer;

    modelPointer=params.modelPointer;

    int dim = params.nVEV - 1;
    std::vector<double> x0(dim);


    double sigma;
    sigma = 5;

    for(int i=0;i<dim;i++) x0[i]=Start.at(i);

    double ftol=1e-5;


    CMAParameters<> cmaparams(x0,sigma);

    cmaparams.set_algo(aCMAES);
    cmaparams.set_ftolerance(ftol);


    CMASolutions cmasols = cmaes<>(CMAES_VEFF_Minimize_Plane,cmaparams);

    Candidate bcand = cmasols.best_candidate();

    std::vector<double> xsol = bcand.get_x();

    std::vector<double> sol;

    for(int i = 0;i<dim;i++)
    {
        sol.push_back(xsol.at(i));
    }

    LibCMAESReturn res;
    res.result = sol;
    res.CMAESStatus = cmasols.run_status();

    return res;


}

}
}
}
