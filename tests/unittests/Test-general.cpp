// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Approx = Catch::Approx;

#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/ModelTestfunctions.h>

TEST_CASE("Checking CKM Unitarity", "[general]")
{
  using namespace BSMPT;
  REQUIRE(ModelTests::CheckCKMUnitarity() == ModelTests::TestResults::Pass);
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
  std::string foundException;
  try
  {
    auto ids = BSMPT::ModelID::InvertModelNames();
  }
  catch (std::exception &e)
  {
    foundException = e.what();
  }
  REQUIRE(foundException == std::string());
}
