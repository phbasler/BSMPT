// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Approx = Catch::Approx;

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/modeltests/ModelTestfunctions.h>
#include <BSMPT/utility/utility.h>

#include "N2HDM.h"

const std::vector<double> example_point_N2HDM{/* lambda_1 = */ 0.300812,
                                              /* lambda_2 = */ 0.321809,
                                              /* lambda_3 = */ -0.133425,
                                              /* lambda_4 = */ 4.11105,
                                              /* lambda_5 = */ -3.84178,
                                              /* lambda_6 = */ 9.46329,
                                              /* lambda_7 = */ -0.750455,
                                              /* lambda_8 = */ 0.743982,
                                              /* tan(beta) = */ 5.91129,
                                              /* v_s = */ 293.035,
                                              /* m_{12}^2 = */ 4842.28,
                                              /* Yukawa Type = */ 1};

constexpr auto Model = BSMPT::ModelID::ModelIDs::N2HDM;
const Compare_N2HDM Expected;

TEST_CASE("Checking NLOVEV for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(Model, SMConstants);
  modelPointer->initModel(example_point_N2HDM);
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
    REQUIRE(res == Approx(expected).margin(1e-4));
  }
}

TEST_CASE("Checking EWPT for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(Model, SMConstants);
  modelPointer->initModel(example_point_N2HDM);
  std::vector<double> Check;
  auto EWPT = Minimizer::PTFinder_gen_all(
      modelPointer, 0, 300, Minimizer::WhichMinimizerDefault);
  const double omega_c_expected =
      Expected.EWPTPerSetting.at(Minimizer::WhichMinimizerDefault).vc;
  const double Tc_expected =
      Expected.EWPTPerSetting.at(Minimizer::WhichMinimizerDefault).Tc;
  const std::vector<double> min_expected =
      Expected.EWPTPerSetting.at(Minimizer::WhichMinimizerDefault).EWMinimum;
  REQUIRE(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS);
  REQUIRE(std::abs(EWPT.vc) == Approx(omega_c_expected).epsilon(1e-2));
  REQUIRE(EWPT.Tc == Approx(Tc_expected).epsilon(1e-2));
  const double threshold = 1e-2;
  for (std::size_t i{0}; i < EWPT.EWMinimum.size(); ++i)
  {
    auto res      = std::abs(EWPT.EWMinimum.at(i));
    auto expected = std::abs(min_expected.at(i));
    if (expected > threshold)
    {
      REQUIRE(res == Approx(expected).epsilon(1e-2));
    }
    else
    {
      REQUIRE(res <= threshold);
    }
  }
}

TEST_CASE("Checking number of CT parameters for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);
  auto result = ModelTests::CheckNumberOfCTParameters(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking number of VEV labels for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);
  auto result = ModelTests::CheckNumberOfVEVLabels(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE(
    "Checking number of labels for temperature dependend results for N2HDM",
    "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);
  auto result = ModelTests::CheckLegendTemp(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking number of triple Higgs couplings for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);
  auto result = ModelTests::CheckNumberOfTripleCouplings(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking Gauge Boson masses for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);
  auto result = ModelTests::CheckGaugeBosonMasses(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking fermion and quark masses masses for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);
  auto result = ModelTests::CheckFermionicMasses(*modelPointer);
  REQUIRE(result.first == ModelTests::TestResults::Pass);
  REQUIRE(result.second == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking tree level minimum for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);
  auto result = ModelTests::CheckTreeLevelMin(*modelPointer,
                                              Minimizer::WhichMinimizerDefault);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking tree level tadpoles for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);
  auto result = ModelTests::CheckTadpoleRelations(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking NLO masses matching tree level masses for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);
  auto result = ModelTests::CheckNLOMasses(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking VTreeSimplified for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  if (modelPointer->UseVTreeSimplified)
  {
    modelPointer->initModel(example_point_N2HDM);
    auto result = ModelTests::CheckVTreeSimplified(*modelPointer);
    REQUIRE(result == ModelTests::TestResults::Pass);
  }
  else
  {
    REQUIRE(true);
  }
}

TEST_CASE("Checking VCounterSimplified for N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  if (modelPointer->UseVCounterSimplified)
  {
    modelPointer->initModel(example_point_N2HDM);
    auto result = ModelTests::CheckVCounterSimplified(*modelPointer);
    REQUIRE(result == ModelTests::TestResults::Pass);
  }
  else
  {
    REQUIRE(true);
  }
}

TEST_CASE("Checking first derivative of the sum of CT and CW in the N2HDM",
          "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);
  auto result = ModelTests::CheckCTConditionsFirstDerivative(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking second derivative of the sum of CT and CW in the N2HDM",
          "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);
  auto result = ModelTests::CheckCTConditionsSecondDerivative(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking triple higgs NLO couplings in the N2HDM", "[n2hdm]")
{

  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);
  modelPointer->Prepare_Triple();
  modelPointer->TripleHiggsCouplings();

  auto Check = [](auto result, auto expected)
  {
    if (std::abs(expected) > 1e-4)
    {
      REQUIRE(result == Approx(expected).epsilon(1e-4));
    }
    else
    {
      REQUIRE(std::abs(result) < 1e-4);
    }
  };

  auto NHiggs = modelPointer->get_NHiggs();
  for (std::size_t i{0}; i < NHiggs; ++i)
  {
    for (std::size_t j{0}; j < NHiggs; ++j)
    {
      for (std::size_t k{0}; k < NHiggs; ++k)
      {
        Check(modelPointer->get_TripleHiggsCorrectionsTreePhysical(i, j, k),
              Expected.CheckTripleTree.at(i).at(j).at(k));
        Check(modelPointer->get_TripleHiggsCorrectionsCTPhysical(i, j, k),
              Expected.CheckTripleCT.at(i).at(j).at(k));
        Check(modelPointer->get_TripleHiggsCorrectionsCWPhysical(i, j, k),
              Expected.CheckTripleCW.at(i).at(j).at(k));
      }
    }
  }
}

TEST_CASE("Check number of calculated CT parameters in the N2HDM", "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);
  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckCTNumber(*modelPointer));
}

TEST_CASE("Check symmetric properties of the scalar tensor Lij in the N2HDM",
          "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorScalarSecond(
              modelPointer->Get_Curvature_Higgs_L2()));
}

TEST_CASE("Check symmetric properties of the scalar tensor Lijk in the N2HDM",
          "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorScalarThird(
              modelPointer->Get_Curvature_Higgs_L3()));
}

TEST_CASE("Check symmetric properties of the scalar tensor Lijkl in the N2HDM",
          "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorScalarFourth(
              modelPointer->Get_Curvature_Higgs_L4()));
}

TEST_CASE("Check symmetric properties of the gauge tensor in the N2HDM",
          "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorGauge(
              modelPointer->Get_Curvature_Gauge_G2H2()));
}

TEST_CASE("Check symmetric properties of the Lepton tensor in the N2HDM",
          "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorLeptonsThird(
              modelPointer->Get_Curvature_Lepton_F2H1()));
}

TEST_CASE("Check symmetric properties of the mass Lepton tensor in the N2HDM",
          "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorLeptons(
              modelPointer->Get_Curvature_Lepton_F2()));
}

TEST_CASE("Check symmetric properties of the mass quark tensor in the N2HDM",
          "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorQuarks(
              modelPointer->Get_Curvature_Quark_F2()));
}

TEST_CASE("Check symmetric properties of the quark tensor in the N2HDM",
          "[n2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  modelPointer->initModel(example_point_N2HDM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorQuarksThird(
              modelPointer->Get_Curvature_Quark_F2H1()));
}
