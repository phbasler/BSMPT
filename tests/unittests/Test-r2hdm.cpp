// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Approx = Catch::Approx;

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/ClassPotentialR2HDM.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/modeltests/ModelTestfunctions.h>

#include "R2HDM.h"

const Compare_R2HDM Expected;

const std::vector<double> example_point_R2HDM{/* lambda_1 = */ 2.740595,
                                              /* lambda_2 = */ 0.242356,
                                              /* lambda_3 = */ 5.534491,
                                              /* lambda_4 = */ -2.585467,
                                              /* lambda_5 = */ -2.225991,
                                              /* m_{12}^2 = */ 7738.56,
                                              /* tan(beta) = */ 4.63286,
                                              /* Yukawa Type = */ 1};

const std::vector<double> example_point_R2HDM_negTanBeta{
    /* lambda_1 = */ 2.740595,
    /* lambda_2 = */ 0.242356,
    /* lambda_3 = */ 5.534491,
    /* lambda_4 = */ -2.585467,
    /* lambda_5 = */ -2.225991,
    /* m_{12}^2 = */ -7738.56,
    /* tan(beta) = */ -4.63286,
    /* Yukawa Type = */ 1};

TEST_CASE("Checking NLOVEV for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
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
    REQUIRE(res == Approx(expected).margin(1e-4));
  }
}

TEST_CASE("Checking EWPT for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);
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

  REQUIRE(std::abs(EWPT.vc) == Approx(omega_c_expected).epsilon(1e-4));
  REQUIRE(EWPT.Tc == Approx(Tc_expected).epsilon(1e-4));
  const double threshold = 1e-4;
  for (std::size_t i{0}; i < EWPT.EWMinimum.size(); ++i)
  {
    auto res      = std::abs(EWPT.EWMinimum.at(i));
    auto expected = std::abs(min_expected.at(i));
    if (expected > threshold)
    {
      REQUIRE(res == Approx(expected).epsilon(1e-4));
    }
    else
    {
      REQUIRE(res <= threshold);
    }
  }
}

TEST_CASE("Checking sign of SinBeta for pos. TanBeta", "[r2hdm]")
{
  using namespace BSMPT;
  using namespace Models;
  const auto SMConstants = GetSMConstants();
  Class_Potential_R2HDM point(SMConstants);
  point.set_gen(example_point_R2HDM);
  REQUIRE(point.C_SinBeta >= 0);
}

TEST_CASE("Checking sign of SinBeta for neg. TanBeta", "[r2hdm]")
{
  using namespace BSMPT;
  using namespace Models;
  const auto SMConstants = GetSMConstants();
  Class_Potential_R2HDM point(SMConstants);
  point.set_gen(example_point_R2HDM_negTanBeta);
  REQUIRE(point.C_SinBeta <= 0);
}

TEST_CASE("Checking sign of CosBeta for pos. TanBeta", "[r2hdm]")
{
  using namespace BSMPT;
  using namespace Models;
  const auto SMConstants = GetSMConstants();
  Class_Potential_R2HDM point(SMConstants);
  point.set_gen(example_point_R2HDM);
  REQUIRE(point.C_CosBeta >= 0);
}

TEST_CASE("Checking sign of CosBeta for neg. TanBeta", "[r2hdm]")
{
  using namespace BSMPT;
  using namespace Models;
  const auto SMConstants = GetSMConstants();
  Class_Potential_R2HDM point(SMConstants);
  point.set_gen(example_point_R2HDM_negTanBeta);
  REQUIRE(point.C_CosBeta >= 0);
}

TEST_CASE("Checking number of CT parameters for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckNumberOfCTParameters(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking number of VEV labels for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckNumberOfVEVLabels(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE(
    "Checking number of labels for temperature dependend results for R2HDM",
    "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckLegendTemp(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking number of triple Higgs couplings for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckNumberOfTripleCouplings(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking Gauge Boson masses for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckGaugeBosonMasses(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking fermion and quark masses masses for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckFermionicMasses(*modelPointer);
  REQUIRE(result.first == ModelTests::TestResults::Pass);
  REQUIRE(result.second == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking tree level minimum for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckTreeLevelMin(*modelPointer,
                                              Minimizer::WhichMinimizerDefault);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking tree level tadpoles for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckTadpoleRelations(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking NLO masses matching tree level masses for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckNLOMasses(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking VTreeSimplified for R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
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
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
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

TEST_CASE("Checking first derivative of the sum of CT and CW in the R2HDM",
          "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckCTConditionsFirstDerivative(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking second derivative of the sum of CT and CW in the R2HDM",
          "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);
  auto result = ModelTests::CheckCTConditionsSecondDerivative(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking triple higgs NLO couplings in the R2HDM", "[r2hdm]")
{

  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);
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

TEST_CASE("Check number of calculated CT parameters in the R2HDM", "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);
  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckCTNumber(*modelPointer));
}

TEST_CASE("Check symmetric properties of the scalar tensor Lij in the R2HDM",
          "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorScalarSecond(
              modelPointer->Get_Curvature_Higgs_L2()));
}

TEST_CASE("Check symmetric properties of the scalar tensor Lijk in the R2HDM",
          "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorScalarThird(
              modelPointer->Get_Curvature_Higgs_L3()));
}

TEST_CASE("Check symmetric properties of the scalar tensor Lijkl in the R2HDM",
          "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorScalarFourth(
              modelPointer->Get_Curvature_Higgs_L4()));
}

TEST_CASE("Check symmetric properties of the gauge tensor in the R2HDM",
          "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorGauge(
              modelPointer->Get_Curvature_Gauge_G2H2()));
}

TEST_CASE("Check symmetric properties of the Lepton tensor in the R2HDM",
          "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorLeptonsThird(
              modelPointer->Get_Curvature_Lepton_F2H1()));
}

TEST_CASE("Check symmetric properties of the mass Lepton tensor in the R2HDM",
          "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorLeptons(
              modelPointer->Get_Curvature_Lepton_F2()));
}

TEST_CASE("Check symmetric properties of the mass quark tensor in the R2HDM",
          "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorQuarks(
              modelPointer->Get_Curvature_Quark_F2()));
}

TEST_CASE("Check symmetric properties of the quark tensor in the R2HDM",
          "[r2hdm]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorQuarksThird(
              modelPointer->Get_Curvature_Quark_F2H1()));
}
