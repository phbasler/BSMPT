// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Approx = Catch::Approx;

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/SMparam.h>
#include <BSMPT/models/modeltests/ModelTestfunctions.h>
#include <BSMPT/utility/utility.h>

TEST_CASE("Checking CKM Unitarity", "[general]")
{
  using namespace BSMPT;
  REQUIRE(ModelTests::CheckCKMUnitarity(GetSMConstants()) ==
          ModelTests::TestResults::Pass);
}

TEST_CASE("Check get model", "[general]")
{
  using namespace BSMPT;
  auto result                = ModelID::getModel("c2hdm");
  ModelID::ModelIDs expected = ModelID::ModelIDs::C2HDM;
  REQUIRE(expected == result);
}

TEST_CASE("Check no model ids are set twice", "[general]")
{
  REQUIRE_NOTHROW(BSMPT::InvertMap(BSMPT::ModelID::ModelNames, "double names"));
}

TEST_CASE("Check calculating of minimizer selection", "[general]")
{
  std::vector<BSMPT::Minimizer::MinimizersToUse> minimizers;
  ;
  for (const auto &gsl : {true, false})
  {
    for (const auto &cmaes : {true, false})
    {
      for (const auto &nlopt : {true, false})
      {
        auto whichMin = BSMPT::Minimizer::CalcWhichMinimizer(gsl, cmaes, nlopt);
        auto minToUse = BSMPT::Minimizer::GetMinimizers(whichMin);
        REQUIRE(minToUse.UseCMAES == cmaes);
        REQUIRE(minToUse.UseGSL == gsl);
        REQUIRE(minToUse.UseNLopt == nlopt);
      }
    }
  }
}
