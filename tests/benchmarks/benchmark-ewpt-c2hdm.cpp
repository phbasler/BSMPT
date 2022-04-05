// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch.hpp"
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
}

TEST_CASE("Benchmark NLOVEV for C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  BENCHMARK("NLOVEV")
  {
    std::vector<double> Check;
    return Minimizer::Minimize_gen_all(modelPointer,
                                       0,
                                       Check,
                                       modelPointer->get_vevTreeMin(),
                                       Minimizer::WhichMinimizerDefault);
  };
}

TEST_CASE("Benchmark EWPT for C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  BENCHMARK("EWPT")
  {
    std::vector<double> Check;
    return Minimizer::PTFinder_gen_all(
        modelPointer, 0, 300, Minimizer::WhichMinimizerDefault);
  };
}
