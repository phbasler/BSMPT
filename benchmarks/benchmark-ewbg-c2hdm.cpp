// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"
#include <BSMPT/baryo_calculation/CalculateEtaInterface.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/ModelTestfunctions.h>

#include "GenerateTestCompares/C2HDM.h"

namespace
{
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

template <typename T>
void CompareValues(T expected, T result, double epsilon, double threshold)
{
  if (std::abs(expected) >= threshold)
  {
    REQUIRE(result == Approx(expected).epsilon(epsilon));
  }
  else
  {
    REQUIRE(std::abs(result) <= threshold);
  }
}
} // namespace

TEST_CASE("Benchmark EWBG for C2HDM", "[c2hdm]")
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
  auto min_expected = Expected.vevSymmetricPerSetting.at(WhichMin);
  REQUIRE(vevsymmetricSolution.size() == min_expected.size());
  for (std::size_t i{0}; i < vevsymmetricSolution.size(); ++i)
  {
    auto res      = std::abs(vevsymmetricSolution.at(i));
    auto expected = std::abs(min_expected.at(i));
    UNSCOPED_INFO("Current Option for Minimizer:\t" << WhichMin);
    UNSCOPED_INFO("This is the position:"
                  << i << "\tFound solution =" << vevsymmetricSolution.at(i)
                  << "\tExpected solution = " << min_expected.at(i));
    CompareValues(expected, res, 1e-4, 1e-4);
  }

  // Call: Calculation of eta in the different implemented approaches

  auto config =
      std::pair<std::vector<bool>, int>{std::vector<bool>(5, true), 1};
  const double testVW = Expected.testVW;

  BENCHMARK("Eta C2HDM Calc Benchmark")
  {
    Baryo::CalculateEtaInterface benchmarkEtaInterface(config);
    return benchmarkEtaInterface.CalcEta(testVW,
                                         EWPT.EWMinimum,
                                         vevsymmetricSolution,
                                         EWPT.Tc,
                                         modelPointer,
                                         WhichMin);
  };
}
