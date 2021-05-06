// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch.hpp"
#include <BSMPT/models/ModelTestfunctions.h>

TEST_CASE("Checking CKM Unitarity", "[general]")
{
  using namespace BSMPT;
  REQUIRE(ModelTests::CheckCKMUnitarity() == ModelTests::TestResults::Pass);
}
