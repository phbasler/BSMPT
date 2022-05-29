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
#include <BSMPT/models/ModelTestfunctions.h>
#include <BSMPT/utility/Logger.h>

#include "CXSM.h"
#include <fstream>

const std::vector<double> example_point_CXSM{/* vh = */ 246.219651,
                                             /* vs = */ 540.51152,
                                             /* va = */ 0,
                                             /* ms = */ -10201.707997,
                                             /* lambda = */ 0.516782,
                                             /* delta2 = */ -0.037398,
                                             /* b2 = */ -370585.40704,
                                             /* d2 = */ 2.570175,
                                             /* Reb1 = */ -3722.817741,
                                             /* Imb1 = */ 0,
                                             /* Rea1 = */ 0,
                                             /* Ima1 = */ 0};

const Compare_CXSM Expected;

TEST_CASE("Checking NLOVEV for CXSM", "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);
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

    INFO(i << " (res/expected) = (" << res << "/" << expected << ")");

    REQUIRE(res == Approx(expected).margin(1e-4));
  }
}

TEST_CASE("Checking EWPT for CXSM", "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);
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
      UNSCOPED_INFO("Current Option for Minimizer:\t"
                    << Minimizer::WhichMinimizerDefault);
      UNSCOPED_INFO("This ist the position:"
                    << i << "\tFound solution =" << EWPT.EWMinimum.at(i)
                    << "\tExpected solution = " << min_expected.at(i));
      REQUIRE(res == Approx(expected).epsilon(1e-4));
    }
    else
    {
      REQUIRE(res <= threshold);
    }
  }
}

TEST_CASE("Checking number of CT parameters for CXSM", "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);
  auto result = ModelTests::CheckNumberOfCTParameters(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking number of VEV labels for CXSM", "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);
  auto result = ModelTests::CheckNumberOfVEVLabels(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE(
    "Checking number of labels for temperature dependend results for CXSM",
    "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);
  auto result = ModelTests::CheckLegendTemp(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking number of triple Higgs couplings for CXSM", "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);
  auto result = ModelTests::CheckNumberOfTripleCouplings(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking Gauge Boson masses for CXSM", "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);
  auto result = ModelTests::CheckGaugeBosonMasses(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking fermion and quark masses masses for CXSM", "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);
  auto result = ModelTests::CheckFermionicMasses(*modelPointer);
  REQUIRE(result.first == ModelTests::TestResults::Pass);
  REQUIRE(result.second == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking tree level minimum for CXSM", "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);
  auto result = ModelTests::CheckTreeLevelMin(*modelPointer,
                                              Minimizer::WhichMinimizerDefault);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking tree level tadpoles for CXSM", "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);
  auto result = ModelTests::CheckTadpoleRelations(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking NLO masses matching tree level masses for CXSM", "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);
  auto result = ModelTests::CheckNLOMasses(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking VTreeSimplified for CXSM", "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  if (modelPointer->UseVTreeSimplified)
  {
    modelPointer->initModel(example_point_CXSM);
    auto result = ModelTests::CheckVTreeSimplified(*modelPointer);
    REQUIRE(result == ModelTests::TestResults::Pass);
  }
  else
  {
    REQUIRE(true);
  }
}

TEST_CASE("Checking VCounterSimplified for CXSM", "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  if (modelPointer->UseVCounterSimplified)
  {
    modelPointer->initModel(example_point_CXSM);
    auto result = ModelTests::CheckVCounterSimplified(*modelPointer);
    REQUIRE(result == ModelTests::TestResults::Pass);
  }
  else
  {
    REQUIRE(true);
  }
}

TEST_CASE("Checking first derivative of the sum of CT and CW in the CXSM",
          "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);
  auto result = ModelTests::CheckCTConditionsFirstDerivative(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking second derivative of the sum of CT and CW in the CXSM",
          "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);
  auto result = ModelTests::CheckCTConditionsSecondDerivative(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking the identities required to vanish for the CT in the CXSM",
          "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);
  auto result = ModelTests::CheckCTIdentities(*modelPointer);
  REQUIRE(result == ModelTests::TestResults::Pass);
}

TEST_CASE("Checking triple higgs NLO couplings in the CXSM", "[CXSM]")
{

  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);
  modelPointer->Prepare_Triple();
  modelPointer->TripleHiggsCouplings();

  auto Check = [](auto result, auto expected) {
    if (expected != 0)
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
        INFO("Checking TreePhysical");
        Check(modelPointer->get_TripleHiggsCorrectionsTreePhysical(i, j, k),
              Expected.CheckTripleTree.at(i).at(j).at(k));
        INFO("Checking CTPhysical");
        Check(modelPointer->get_TripleHiggsCorrectionsCTPhysical(i, j, k),
              Expected.CheckTripleCT.at(i).at(j).at(k));
        INFO("Checking CWPhysical");
        Check(modelPointer->get_TripleHiggsCorrectionsCWPhysical(i, j, k),
              Expected.CheckTripleCW.at(i).at(j).at(k));
      }
    }
  }
}

TEST_CASE("Check number of calculated CT parameters in the CXSM", "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);
  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckCTNumber(*modelPointer));
}

TEST_CASE("Check symmetric properties of the scalar tensor Lij in the CXSM",
          "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorScalarSecond(
              modelPointer->Get_Curvature_Higgs_L2()));
}

TEST_CASE("Check symmetric properties of the scalar tensor Lijk in the CXSM",
          "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorScalarThird(
              modelPointer->Get_Curvature_Higgs_L3()));
}

TEST_CASE("Check symmetric properties of the scalar tensor Lijkl in the CXSM",
          "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorScalarFourth(
              modelPointer->Get_Curvature_Higgs_L4()));
}

TEST_CASE("Check symmetric properties of the gauge tensor in the CXSM",
          "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorGauge(
              modelPointer->Get_Curvature_Gauge_G2H2()));
}

TEST_CASE("Check symmetric properties of the Lepton tensor in the CXSM",
          "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorLeptonsThird(
              modelPointer->Get_Curvature_Lepton_F2H1()));
}

TEST_CASE("Check symmetric properties of the mass Lepton tensor in the CXSM",
          "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorLeptons(
              modelPointer->Get_Curvature_Lepton_F2()));
}

TEST_CASE("Check symmetric properties of the mass quark tensor in the CXSM",
          "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorQuarks(
              modelPointer->Get_Curvature_Quark_F2()));
}

TEST_CASE("Check symmetric properties of the quark tensor in the CXSM",
          "[CXSM]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM);
  modelPointer->initModel(example_point_CXSM);

  REQUIRE(ModelTests::TestResults::Pass ==
          ModelTests::CheckSymmetricTensorQuarksThird(
              modelPointer->Get_Curvature_Quark_F2H1()));
}
