#include <BSMPT/minimizer/LibNLOPT/MinimizeNLOPT.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/minimizer/MinimizePlane.h>

/**
 *@file
 */
namespace BSMPT {
namespace Minimizer {
namespace LibNLOPT {



double NLOPTVEff(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    auto settings = *static_cast<ShareInformationNLOPT*>(data);
    (void) grad;
    auto PotVEV = settings.model->MinimizeOrderVEV(x);
    return settings.model->VEff(PotVEV,settings.Temp);
}

NLOPTReturnType MinimizeUsingNLOPT(
        const std::shared_ptr<Class_Potential_Origin>& model,
        const double& Temp)
{
    ShareInformationNLOPT settings;
    settings.Temp = Temp;
    settings.model = model;
    std::vector<double> VEV(model->get_nVEV());

    nlopt::opt opt(nlopt::LN_COBYLA,static_cast<unsigned int>(model->get_nVEV()));

    opt.set_min_objective(NLOPTVEff,&settings);
    opt.set_xtol_rel(1e-4);
    double minf;
    auto result = opt.optimize(VEV,minf);
    NLOPTReturnType res;
    res.PotVal = minf;
    res.Minimum = VEV;
    res.NLOPResult = result;
    return  res;
}

double NLOPTVEffPlane(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    (void) grad;
    auto params = *static_cast<PointerContainerMinPlane*>(data);
    std::vector<double> vMinTilde(std::begin(x),std::end(x));
    auto vev = TransformCoordinates(vMinTilde,params);
    return params.modelPointer->VEff(params.modelPointer->MinimizeOrderVEV(vev),params.Temp);
}

NLOPTReturnType MinimizePlaneUsingNLOPT(
        const struct PointerContainerMinPlane& params)
{
    nlopt::opt opt(nlopt::LN_COBYLA,static_cast<unsigned int>(params.modelPointer->get_nVEV()-1));
    std::vector<double> VEV(params.modelPointer->get_nVEV()-1);
    auto copy = params;
    opt.set_min_objective(NLOPTVEffPlane,&copy);
    opt.set_xtol_rel(1e-4);
    double minf;
    auto result = opt.optimize(VEV,minf);
    NLOPTReturnType res;
    res.PotVal = minf;
    res.Minimum = TransformCoordinates(VEV,params);
    res.NLOPResult = result;
    return  res;
}

}
}
}
