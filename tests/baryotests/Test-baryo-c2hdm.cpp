// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch.hpp"

#include <BSMPT/baryo_calculation/CalculateEtaInterface.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/ModelTestfunctions.h>

#include "GenerateTestCompares/C2HDM.h"

const std::vector<double> example_point_C2HDM{/* lambda_1 = */ 3.29771,
                                              /* lambda_2 = */ 0.274365,
                                              /* lambda_3 = */ 4.71019,
                                              /* lambda_4 = */ -2.23056,
                                              /* Re(lambda_5) = */ -2.43487,
                                              /* Im(lambda_5) = */ 0.124948,
                                              /* Re(m_{12}^2) = */ 2706.86,
                                              /* tan(beta) = */ 4.64487,
                                              /* Yukawa Type = */ 1};

const Compare_C2HDM Expected;

TEST_CASE("Checking EWBG for C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);

  const auto WhichMin = Minimizer::WhichMinimizerDefault;

  const auto EWPT = Expected.EWPTPerSetting.at(WhichMin);

  std::vector<double> vevsymmetricSolution, checksym, startpoint;
  for (const auto &el : EWPT.EWMinimum)
    startpoint.push_back(0.5 * el);
  vevsymmetricSolution = Minimizer::Minimize_gen_all(
      modelPointer, EWPT.Tc + 1, checksym, startpoint, WhichMin, true);
  auto min_expected      = Expected.vevSymmetricPerSetting.at(WhichMin);
  const double threshold = 1e-4;
  REQUIRE(vevsymmetricSolution.size() == min_expected.size());
  for (std::size_t i{0}; i < vevsymmetricSolution.size(); ++i)
  {
    auto res      = std::abs(vevsymmetricSolution.at(i));
    auto expected = std::abs(min_expected.at(i));
    if (expected > threshold)
    {
      UNSCOPED_INFO("Current Option for Minimizer:\t" << WhichMin);
      UNSCOPED_INFO("This ist the position:"
                    << i << "\tFound solution =" << vevsymmetricSolution.at(i)
                    << "\tExpected solution = " << min_expected.at(i));
      REQUIRE(res == Approx(expected).epsilon(1e-4));
    }
    else
    {
      REQUIRE(res <= threshold);
    }
  }

  // Call: Calculation of eta in the different implemented approaches

  auto config =
      std::pair<std::vector<bool>, int>{std::vector<bool>(5, true), 1};
  Baryo::CalculateEtaInterface EtaInterface(config);
  const double testVW = Expected.testVW;

  if (EWPT.vc / EWPT.Tc > 1)
  {
    auto eta      = EtaInterface.CalcEta(testVW,
                                    EWPT.EWMinimum,
                                    vevsymmetricSolution,
                                    EWPT.Tc,
                                    modelPointer,
                                    WhichMin);
    const auto LW = EtaInterface.getLW();
    REQUIRE(LW == Approx(Expected.LWPerSetting.at(WhichMin)).epsilon(1e-4));

    auto expectedEta = Expected.etaPerSetting.at(WhichMin);
    REQUIRE(eta.size() == expectedEta.size());
    const double etaThreshold = 1e-15;
    for (std::size_t i{0}; i < vevsymmetricSolution.size(); ++i)
    {
      auto res      = std::abs(eta.at(i));
      auto expected = std::abs(expectedEta.at(i));
      if (expected > etaThreshold)
      {
        UNSCOPED_INFO("Current Option for Minimizer:\t" << WhichMin);
        UNSCOPED_INFO("This ist the position:"
                      << i << "\tFound solution =" << vevsymmetricSolution.at(i)
                      << "\tExpected solution = " << min_expected.at(i));
        REQUIRE(res == Approx(expected).epsilon(1e-4));
      }
      else
      {
        REQUIRE(res <= threshold);
      }
    }
  }
}
