// SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/bounce_solution/action_calculation.h>
#include <BSMPT/gravitational_waves/gw.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/minimum_tracer/minimum_tracer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/modeltests/ModelTestfunctions.h>
#include <BSMPT/transition_tracer/transition_tracer.h>
#include <BSMPT/utility/Logger.h> // for Logger Class
#include <fstream>

// C++ version of the python function ", ".join()
std::string join(const std::vector<double> &sequence,
                 const std::string &separator)
{
  std::string result;
  for (size_t i = 0; i < sequence.size(); ++i)
    result += std::to_string(sequence[i]) +
              ((i != sequence.size() - 1) ? separator : "");
  return result;
}

using namespace BSMPT;
using namespace Minimizer;

int main()
{
  // Set the logger
  SetLogger({"--logginglevel::complete=true"});

  const std::vector<double> example_point_CXSM{/* v = */ 245.34120667410863,
                                               /* vs = */ 0,
                                               /* va = */ 0,
                                               /* msq = */ -15650,
                                               /* lambda = */ 0.52,
                                               /* delta2 = */ 0.55,
                                               /* b2 = */ -8859,
                                               /* d2 = */ 0.5,
                                               /* Reb1 = */ 0,
                                               /* Imb1 = */ 0,
                                               /* Rea1 = */ 0,
                                               /* Ima1 = */ 0};

  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM, SMConstants);
  modelPointer->initModel(example_point_CXSM);

  user_input input{
      modelPointer,                       /*modelPointer*/
      0,                                  /*templow*/
      400,                                /*temphigh*/
      0.95,                               /*UserDefined_vwall*/
      .71,                                /*perc_prbl*/
      .01,                                /*compl_prbl*/
      0.1,                                /*UserDefined_epsturb*/
      7,                                  /*MaxPathIntegrations*/
      -1,                                 /*UseMultiStepPTMode*/
      10,                                 /*num_check_pts*/
      0,                                  /*CheckEWSymmetryRestoration*/
      1,                                  /*CheckNLOStability*/
      WhichMinimizerDefault,              /*WhichMinimizer*/
      false,                              /*GW calculation*/
      true,                               /*gw_calculation*/
      TransitionTemperature::Percolation, /*WhichTransitionTemperature*/
      1};                                 /*userDefined_PNLO_scaling*/

  TransitionTracer trans(input);
  auto output = trans.output_store;

  for (auto &bounce : trans.ListBounceSolution)
  {
    // Calculate percolation temperature
    bounce.CalculatePercolationTemp();
    if (bounce.GetPercolationTemp() == -1) continue;
    double errorTtoTp = 1e100;
    BounceActionInt const *ClosestBounceActionInt{nullptr};
    std::cout << "Found a transitions with Tp =\t"
              << bounce.GetPercolationTemp() << " GeV.\n";
    for (const auto &BAInt : bounce.SolutionList)
    {
      if (abs(BAInt.T - bounce.GetPercolationTemp()) < errorTtoTp)
      {
        errorTtoTp             = abs(BAInt.T - bounce.GetPercolationTemp());
        ClosestBounceActionInt = &BAInt;
      }
    }

    if (ClosestBounceActionInt == nullptr)
    {
      throw std::runtime_error("ClosestBounceActionInt was not set");
    }
    std::cout << "The closest solution is at a distance of " << errorTtoTp
              << " GeV from Tp.\n The tunnenling path is\n";
    ClosestBounceActionInt->Spline.print_path();
    std::cout << "\n ------ The tunneling profile in Mathematica format ------ "
                 "\n {rho, l, point}\n\n";

    // Save the path and VEV profile into, e.g, Mathematica
    for (size_t i = 0; i < ClosestBounceActionInt->rho_sol.size(); i++)
      std::cout << "{" << ClosestBounceActionInt->rho_sol[i] << ", "
                << ClosestBounceActionInt->l_sol[i] << ", "
                << join(ClosestBounceActionInt->Spline(
                            ClosestBounceActionInt->l_sol[i]),
                        ", ")
                << "},\n";
    std::cout << "\n ------ The tunneling profile in python format ------ "
                 "\n {rho, l, point}\n\n";
    // Save the path and VEV profile into, e.g, python
    for (size_t i = 0; i < ClosestBounceActionInt->rho_sol.size(); i++)
      std::cout << "[" << ClosestBounceActionInt->rho_sol[i] << ", "
                << ClosestBounceActionInt->l_sol[i] << ", "
                << join(ClosestBounceActionInt->Spline(
                            ClosestBounceActionInt->l_sol[i]),
                        ", ")
                << "],\n";
  }
  return 0;
}
