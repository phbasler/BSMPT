// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch.hpp"

#include "GenerateTestCompares/C2HDM.h"
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/IncludeAllModels.h>

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

TEST_CASE("Check f_{abcd}", "[origin]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  modelPointer->resetScale(200);

  {
    double expected = -0.30575009e-3;
    double result   = modelPointer->fbaseFour(10, 20, 30, 40);
    REQUIRE(result == Approx(expected).margin(1e-4));
  }
  {
    double expected = -3.518370054e-7;
    double result   = modelPointer->fbaseFour(1e4, 200, 300, 5);
    REQUIRE(result == Approx(expected).margin(1e-4));
  }

  {
    double expected = -0.23248145e-3;
    double result   = modelPointer->fbaseFour(20, 20, 30, 40);
    REQUIRE(result == Approx(expected).margin(1e-4));
  }

  {
    double expected = -0.284264097e-3;
    double result   = modelPointer->fbaseFour(20, 20, 20, 40);
    REQUIRE(result == Approx(expected).margin(1e-4));
  }

  {
    double expected = -0.4166666667e-3;
    double result   = modelPointer->fbaseFour(20, 20, 20, 20);
    REQUIRE(result == Approx(expected).margin(1e-4));
  }

  {
    double expected = -3.696784963e-9;
    double result   = modelPointer->fbaseFour(
        0, std::pow(100, 2), std::pow(200, 2), std::pow(50, 2));
    REQUIRE(result == Approx(expected).margin(1e-4));
  }
}
