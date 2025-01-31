// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <benchmark/benchmark.h>

#include <BSMPT/baryo_calculation/CalculateEtaInterface.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/modeltests/ModelTestfunctions.h>

#include "C2HDM.h"

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

static void BM_EWBG(benchmark::State &state)
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM, SMConstants);
  modelPointer->initModel(example_point_C2HDM);

  const auto WhichMin = Minimizer::WhichMinimizerDefault;

  const auto EWPT = Expected.EWPTPerSetting.at(WhichMin);

  std::vector<double> vevsymmetricSolution, checksym, startpoint;
  for (const auto &el : EWPT.EWMinimum)
    startpoint.push_back(0.5 * el);
  vevsymmetricSolution = Minimizer::Minimize_gen_all(
      modelPointer, EWPT.Tc + 1, checksym, startpoint, WhichMin, true);
  auto min_expected = Expected.vevSymmetricPerSetting.at(WhichMin);

  auto config =
      std::pair<std::vector<bool>, int>{std::vector<bool>(5, true), 1};
  const double testVW = Expected.testVW;

  for (auto _ : state)
  {
    Baryo::CalculateEtaInterface benchmarkEtaInterface(config, SMConstants);
    auto result = benchmarkEtaInterface.CalcEta(testVW,
                                                EWPT.EWMinimum,
                                                vevsymmetricSolution,
                                                EWPT.Tc,
                                                modelPointer,
                                                WhichMin);
    (void)result;
  }
}

BENCHMARK(BM_EWBG)->Repetitions(5);
