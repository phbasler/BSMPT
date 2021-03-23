#include "catch.hpp"
#include <BSMPT/models/ModelTestfunctions.h>

TEST_CASE("Checking CKM Unitarity", "[general]")
{
  using namespace BSMPT;
  REQUIRE(ModelTests::CheckCKMUnitarity() == ModelTests::TestResults::Pass);
}
