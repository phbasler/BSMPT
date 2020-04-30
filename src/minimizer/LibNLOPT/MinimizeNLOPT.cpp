#include <BSMPT/minimizer/LibNLOPT/MinimizeNLOPT.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/IncludeAllModels.h>

/**
 *@file
 */
namespace BSMPT {
namespace Minimizer {
namespace LibNLOPT {



double NLOPTVEff(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    struct ShareInformationNLOPT * settings =  static_cast<ShareInformationNLOPT*>(data);
    (void) grad;
    auto PotVEV = settings->model->MinimizeOrderVEV(x);
    return settings->model->VEff(PotVEV,settings->Temp);
}

NLOPTReturnType MinimizeUsingNLOPT(
        const std::shared_ptr<Class_Potential_Origin>& model,
        const double& Temp)
{
    ShareInformationNLOPT settings;
    settings.Temp = Temp;
    settings.model = model;
    std::vector<double> VEV(model->get_nVEV());

    nlopt::opt opt(nlopt::LN_COBYLA,static_cast<uint>(model->get_nVEV()));

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
}
}
}
