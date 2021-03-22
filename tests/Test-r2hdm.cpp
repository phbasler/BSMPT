#include "catch.hpp"

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/ModelTestfunctions.h>
#include <BSMPT/utility.h>

const std::vector<double> example_point_R2HDM{/* lambda_1 = */ 2.740595,
                                              /* lambda_2 = */ 0.242356,
                                              /* lambda_3 = */ 5.534491,
                                              /* lambda_4 = */ -2.585467,
                                              /* lambda_5 = */ -2.225991,
                                              /* m_{12}^2 = */ 7738.56,
                                              /* tan(beta) = */ 4.63286,
                                              /* Yukawa Type = */ 1};

TEST_CASE("Checking NLOVEV for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM);
  modelPointer->initModel(example_point_R2HDM);
  std::vector<double> Check;
  auto sol = Minimizer::Minimize_gen_all(modelPointer,
                                         0,
                                         Check,
                                         modelPointer->get_vevTreeMin(),
                                         Minimizer::WhichMinimizerDefault);
  for (std::size_t i{0}; i < sol.size(); ++i)
  {
    auto expected = std::abs(modelPointer->get_vevTreeMin(i));
    auto res      = std::abs(sol.at(i));
    REQUIRE(std::abs(res - expected) <= 1e-4);
  }
}

TEST_CASE("Checking EWPT for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM);
  modelPointer->initModel(example_point_R2HDM);
  std::vector<double> Check;
  auto EWPT = Minimizer::PTFinder_gen_all(
      modelPointer, 0, 300, Minimizer::WhichMinimizerDefault);
  const double omega_c_expected = 190.7376740105783;
  const double Tc_expected      = 155.2825927734375;
  const std::vector<double> min_expected{
      0, -54.246714903745641, -182.86102430293158, 0};
  REQUIRE(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS);
  REQUIRE(std::abs(omega_c_expected - EWPT.vc) / omega_c_expected <= 1e-4);
  REQUIRE(std::abs(Tc_expected - EWPT.Tc) / Tc_expected <= 1e-4);
  for (std::size_t i{0}; i < EWPT.EWMinimum.size(); ++i)
  {
    auto res      = std::abs(EWPT.EWMinimum.at(i));
    auto expected = std::abs(min_expected.at(i));
    if (expected != 0)
    {
      REQUIRE(std::abs(res - expected) / expected <= 1e-4);
    }
    else
    {
      REQUIRE(res <= 1e-4);
    }
  }
}

TEST_CASE("Checking number of CT parameters for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckNumberOfCTParameters(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking number of VEV labels for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckNumberOfVEVLabels(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE(
    "Checking number of labels for temperature dependend results for R2HDM",
    "[r2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckLegendTemp(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking number of triple Higgs couplings for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckNumberOfTripleCouplings(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking Gauge Boson masses for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckGaugeBosonMasses(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking fermion and quark masses masses for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckFermionicMasses(*modelPointer);
  REQUIRE(result.first == ModelTests::TestResults::Pass);
  REQUIRE(result.second == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking tree level minimum for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckTreeLevelMin(*modelPointer,
                                              Minimizer::WhichMinimizerDefault);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking tree level tadpoles for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckTadpoleRelations(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking NLO masses matching tree level masses for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckNLOMasses(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking VTreeSimplified for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM);
  if (modelPointer->UseVTreeSimplified)
  {
    modelPointer->initModel(example_point_R2HDM);
    auto result = ModelTests::CheckVTreeSimplified(*modelPointer);
    REQUIRE(result == ModelTests::TestResults::Pass);
  }
  else
  {
    REQUIRE(true);
  }
}

TEST_CASE("Checking VCounterSimplified for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM);
  if (modelPointer->UseVCounterSimplified)
  {
    modelPointer->initModel(example_point_R2HDM);
    auto result = ModelTests::CheckVCounterSimplified(*modelPointer);
    REQUIRE(result == ModelTests::TestResults::Pass);
  }
  else
  {
    REQUIRE(true);
  }
}

TEST_CASE("Checking CKM Unitarity for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  REQUIRE(ModelTests::CheckCKMUnitarity() == ModelTests::TestResults::Pass);
}
