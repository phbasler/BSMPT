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

using namespace BSMPT;
using namespace Minimizer;

class Class_Potential_OriginDerived : public Class_Potential_Origin
{
public:
  ////////////////////////////////////////////////////////////
  void ReadAndSet(const std::string &linestr, std::vector<double> &par) override
  {
    (void)linestr;
    (void)par;
    return;
  };
  std::vector<std::string> addLegendCT() const override { return {}; };
  std::vector<std::string> addLegendTemp() const override { return {}; };
  std::vector<std::string> addLegendTripleCouplings() const override
  {
    return {};
  };

  void set_gen(const std::vector<double> &par) override
  {
    (void)par;
    return;
  };
  void set_CT_Pot_Par(const std::vector<double> &par) override
  {
    (void)par;
    return;
  };
  void write() const override { return; };
  void SetCurvatureArrays() override { return; };
  bool CalculateDebyeSimplified() override { return true; };
  bool CalculateDebyeGaugeSimplified() override { return true; };
  double VTreeSimplified(const std::vector<double> &v) const override
  {
    (void)v;
    return 0;
  };
  double VCounterSimplified(const std::vector<double> &v) const override
  {
    (void)v;
    return 0;
  };
  void AdjustRotationMatrix() override { return; };
  void TripleHiggsCouplings() override { return; };
  std::vector<double> calc_CT() const override { return {0}; };
  void Debugging(const std::vector<double> &input,
                 std::vector<double> &output) const override
  {
    (void)input;
    (void)output;
    return;
  };
  ////////////////////////////////////////////////////////////
  Class_Potential_OriginDerived() : Class_Potential_Origin(GetSMConstants())
  {
    NNeutralHiggs = 4; // number of neutral Higgs bosons at T = 0, w/ goldstones
    NChargedHiggs = 0; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

    NHiggs = NNeutralHiggs + NChargedHiggs;

    nVEV = 1; // number of VEVs to minimize the potential

    // Treelevel minimum
    vevTreeMin = {246.22};

    // This is the mapping from the "manifold vacuum fields" into the full
    // fields.
    VevOrder.resize(nVEV);
    VevOrder[0] = 0; // omegaCB
  }

  // Labels for the VEVS
  std::vector<std::string> addLegendVEV() const override { return {"v"}; };

  /**
   * This function calculates the EW breaking VEV from all contributing field
   * configurations.
   */
  double EWSBVEV(const std::vector<double> &v) const override
  {
    (void)v;
    return 0;
  }

  // Set all EW breaking directions to zero
  void SetEWVEVZero(std::vector<double> &sol) const override
  {
    sol = std::vector<double>(nVEV, 0);
  }

  double VEff(const std::vector<double> &v,
              double Temp = 0,
              int diff    = 0,
              int Order   = 1) const override
  {
    (void)diff;
    (void)Order;
    double r = 0.152808 * (pow(Temp, 2) - pow(160, 2)) * pow(v[0], 2) +
               0.0322634 * pow(v[0], 4);
    return r;
  }
};

int main()
{
  // Set the logger
  SetLogger({"--logginglevel::complete=true"});

  std::shared_ptr<Class_Potential_OriginDerived> modelPointerDerived;
  modelPointerDerived = std::make_shared<Class_Potential_OriginDerived>();

  std::shared_ptr<Class_Potential_Origin> modelPointer =
      std::dynamic_pointer_cast<Class_Potential_OriginDerived>(
          modelPointerDerived);

  user_input input{
      modelPointer,                       /* modelPointer */
      0,                                  /* templow */
      400,                                /* temphigh */
      0.95,                               /* UserDefined_vwall */
      .71,                                /* perc_prbl */
      .01,                                /* compl_prbl */
      0.1,                                /* UserDefined_epsturb */
      7,                                  /* MaxPathIntegrations */
      -1,                                 /* UseMultiStepPTMode */
      10,                                 /* num_check_pts  */
      0,                                  /* CheckEWSymmetryRestoration*/
      0,                                  /* CheckNLOStability*/
      WhichMinimizerDefault,              /* WhichMinimizer*/
      false,                              /* use multithreading */
      true,                               /* gw calculation */
      TransitionTemperature::Percolation, /* WhichTransitionTemperature */
      1};                                 /* UserDefined_PNLO_scaling */

  TransitionTracer trans(input);

  auto output = trans.output_store; // Results are here. Check the other classes
                                    // to learn how to access it.

  return 0;
}
