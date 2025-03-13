// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <benchmark/benchmark.h>

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
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
}

static void BM_NLOVEV(benchmark::State &state)
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM, SMConstants);
  modelPointer->initModel(example_point_C2HDM);
  for (auto _ : state)
  {
    std::vector<double> Check;
    auto result = Minimizer::Minimize_gen_all(modelPointer,
                                              0,
                                              Check,
                                              modelPointer->get_vevTreeMin(),
                                              Minimizer::WhichMinimizerDefault);
    (void)result;
  }
}

static void BM_EWPT(benchmark::State &state)
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM, SMConstants);
  modelPointer->initModel(example_point_C2HDM);
  for (auto _ : state)
  {
    std::vector<double> Check;
    auto result = Minimizer::PTFinder_gen_all(
        modelPointer, 0, 300, Minimizer::WhichMinimizerDefault);
    (void)result;
  }
}

BENCHMARK(BM_NLOVEV);
BENCHMARK(BM_EWPT)->Repetitions(5);
BENCHMARK_MAIN();
