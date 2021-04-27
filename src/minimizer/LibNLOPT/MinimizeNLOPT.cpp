// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/minimizer/LibNLOPT/MinimizeNLOPT.h>
#include <BSMPT/minimizer/MinimizePlane.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility.h>

/**
 *@file
 */
namespace BSMPT
{
namespace Minimizer
{
namespace LibNLOPT
{

double
NLOPTVEff(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
  auto settings = *static_cast<ShareInformationNLOPT *>(data);
  (void)grad;
  auto PotVEV = settings.model->MinimizeOrderVEV(x);
  return settings.model->VEff(PotVEV, settings.Temp);
}

NLOPTReturnType
MinimizeUsingNLOPT(const std::shared_ptr<Class_Potential_Origin> &model,
                   const double &Temp)
{
  ShareInformationNLOPT settings(model, Temp);
  std::vector<double> VEV(model->get_nVEV());

  nlopt::opt opt(nlopt::GN_ORIG_DIRECT_L,
                 static_cast<unsigned int>(model->get_nVEV()));
  std::vector<double> LowerBound(model->get_nVEV(), -300),
      UpperBound(model->get_nVEV(), 300);
  for (std::size_t i{0}; i < model->get_nVEV(); ++i)
  {
    if (std::abs(model->get_vevTreeMin(i)) > UpperBound.at(i))
    {
      UpperBound.at(i) = 1.5 * model->get_vevTreeMin(i);
    }
    if (-std::abs(model->get_vevTreeMin(i)) < LowerBound.at(i))
    {
      LowerBound.at(i) = -1.5 * std::abs(model->get_vevTreeMin(i));
    }
  }
  opt.set_lower_bounds(LowerBound);
  opt.set_upper_bounds(UpperBound);

  opt.set_min_objective(NLOPTVEff, &settings);
  opt.set_xtol_rel(1e-4);
  double minf;
  auto result  = opt.optimize(VEV, minf);
  bool Success = (result == nlopt::SUCCESS) or
                 (result == nlopt::FTOL_REACHED) or
                 (result == nlopt::XTOL_REACHED);
  NLOPTReturnType res(VEV, minf, result, Success);
  return res;
}

double NLOPTVEffPlane(const std::vector<double> &x,
                      std::vector<double> &grad,
                      void *data)
{
  (void)grad;
  auto params = *static_cast<PointerContainerMinPlane *>(data);
  std::vector<double> vMinTilde(std::begin(x), std::end(x));
  auto vev = TransformCoordinates(vMinTilde, params);
  return params.modelPointer->VEff(params.modelPointer->MinimizeOrderVEV(vev),
                                   params.Temp);
}

NLOPTReturnType
MinimizePlaneUsingNLOPT(const struct PointerContainerMinPlane &params)
{
  auto dim = params.modelPointer->get_nVEV() - 1;
  nlopt::opt opt(nlopt::GN_ORIG_DIRECT_L, static_cast<unsigned int>(dim));
  std::vector<double> LowerBound(dim, -300), UpperBound(dim, 300);
  for (std::size_t i{0}; i < params.modelPointer->get_nVEV(); ++i)
  {
    if (i < params.Index)
    {
      if (std::abs(params.modelPointer->get_vevTreeMin(i)) > UpperBound.at(i))
      {
        UpperBound.at(i) = 1.5 * params.modelPointer->get_vevTreeMin(i);
      }
      if (-std::abs(params.modelPointer->get_vevTreeMin(i)) < LowerBound.at(i))
      {
        LowerBound.at(i) =
            -1.5 * std::abs(params.modelPointer->get_vevTreeMin(i));
      }
    }
    else if (i > params.Index)
    {
      if (std::abs(params.modelPointer->get_vevTreeMin(i)) >
          UpperBound.at(i - 1))
      {
        UpperBound.at(i - 1) = 1.5 * params.modelPointer->get_vevTreeMin(i);
      }
      if (-std::abs(params.modelPointer->get_vevTreeMin(i)) <
          LowerBound.at(i - 1))
      {
        LowerBound.at(i - 1) =
            -1.5 * std::abs(params.modelPointer->get_vevTreeMin(i));
      }
    }
  }
  opt.set_lower_bounds(LowerBound);
  opt.set_upper_bounds(UpperBound);

  std::vector<double> VEV(params.modelPointer->get_nVEV() - 1);
  auto copy = params;
  opt.set_min_objective(NLOPTVEffPlane, &copy);
  opt.set_xtol_rel(1e-4);
  double minf;
  auto result  = opt.optimize(VEV, minf);
  auto minimum = TransformCoordinates(VEV, params);
  bool Success = (result == nlopt::SUCCESS) or
                 (result == nlopt::FTOL_REACHED) or
                 (result == nlopt::XTOL_REACHED);
  NLOPTReturnType res(minimum, minf, result, Success);
  return res;
}

NLOPTReturnType
FindLocalMinimum(const std::shared_ptr<Class_Potential_Origin> &model,
                 const std::vector<double> &Start,
                 const double &Temp)
{
  nlopt::opt opt(nlopt::LN_COBYLA,
                 static_cast<unsigned int>(model->get_nVEV()));
  ShareInformationNLOPT settings(model, Temp);
  auto VEV = Start;
  opt.set_min_objective(NLOPTVEff, &settings);
  opt.set_xtol_rel(1e-4);

  double minf;
  auto result  = opt.optimize(VEV, minf);
  bool Success = (result == nlopt::SUCCESS) or
                 (result == nlopt::FTOL_REACHED) or
                 (result == nlopt::XTOL_REACHED);
  NLOPTReturnType res(VEV, minf, result, Success);

  return res;
}

} // namespace LibNLOPT
} // namespace Minimizer
} // namespace BSMPT
