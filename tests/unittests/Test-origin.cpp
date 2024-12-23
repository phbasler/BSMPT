// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <gmock/gmock.h>

using Approx = Catch::Approx;

#include "C2HDM.h"
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/modeltests/ModelTestfunctions.h>

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

} // namespace

/**
 * @test Check if the automatic Debye corrections match the SM one in the SM
 * case. This should be y_t^2/4.
 */
TEST_CASE("Test Calculate Debye", "[origin]")
{

  using namespace BSMPT;
  ISMConstants SM;
  SM.C_MassTop = 172;
  SM.C_vev0    = 246;

  // This is not a legal point but just a dummy point to only have the top
  // coupling in the Debye Contributions
  const std::vector<double> example_point_CXSM{/* vh = */ SM.C_vev0,
                                               /* vs = */ 0,
                                               /* va = */ 0,
                                               /* ms = */ 41.67,
                                               /* lambda = */ 0,
                                               /* delta2 = */ 0,
                                               /* b2 = */ 0,
                                               /* d2 = */ 0,
                                               /* Reb1 = */ 0,
                                               /* Imb1 = */ 0,
                                               /* Rea1 = */ 0,
                                               /* Ima1 = */ 0};

  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM, SM);
  modelPointer->initModel(example_point_CXSM);
  modelPointer->CalculateDebye(true);
  auto debye = modelPointer->get_DebyeHiggs();

  double topCoupling = SM.C_MassTop * std::sqrt(2) / SM.C_vev0;
  double expected    = std::pow(topCoupling, 2) / 4.0;

  auto calculated = debye.at(3).at(3);

  REQUIRE(calculated != 0);

  REQUIRE(calculated == Approx(expected).margin(1e-4));
}

TEST_CASE("Check f_{abcd}", "[origin]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM, SMConstants);
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

TEST_CASE("Check CheckImplementation", "[origin]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM, SMConstants);
  modelPointer->initModel(example_point_C2HDM);

  ModelTests::CheckImplementation(*modelPointer,
                                  Minimizer::WhichMinimizerDefault);
}

class Mock_Class_Potential_Origin : public BSMPT::Class_Potential_Origin
{
public:
  Mock_Class_Potential_Origin()
      : BSMPT::Class_Potential_Origin(BSMPT::GetSMConstants()) {};
  MOCK_METHOD(double,
              VEff,
              (const std::vector<double> &v, double Temp, int diff, int Order),
              (override, const));
  MOCK_METHOD(double,
              EWSBVEV,
              (const std::vector<double> &v),
              (override, const));
  MOCK_METHOD(void,
              SetEWVEVZero,
              (std::vector<double> & sol),
              (override, const));
  MOCK_METHOD(void,
              ReadAndSet,
              (const std::string &linestr, std::vector<double> &par),
              (override));
  MOCK_METHOD(std::vector<std::string>, addLegendCT, (), (override, const));
  MOCK_METHOD(std::vector<std::string>, addLegendTemp, (), (override, const));
  MOCK_METHOD(std::vector<std::string>,
              addLegendTripleCouplings,
              (),
              (override, const));
  MOCK_METHOD(std::vector<std::string>, addLegendVEV, (), (override, const));
  MOCK_METHOD(void, set_gen, (const std::vector<double> &par), (override));
  MOCK_METHOD(void,
              set_CT_Pot_Par,
              (const std::vector<double> &par),
              (override));
  MOCK_METHOD(void, write, (), (override, const));
  MOCK_METHOD(void, SetCurvatureArrays, (), (override));
  MOCK_METHOD(bool, CalculateDebyeSimplified, (), (override));
  MOCK_METHOD(bool, CalculateDebyeGaugeSimplified, (), (override));
  MOCK_METHOD(double,
              VTreeSimplified,
              (const std::vector<double> &v),
              (override, const));
  MOCK_METHOD(double,
              VCounterSimplified,
              (const std::vector<double> &v),
              (override, const));
  MOCK_METHOD(void, TripleHiggsCouplings, (), (override));
  MOCK_METHOD(std::vector<double>, calc_CT, (), (override, const));
  MOCK_METHOD(void,
              Debugging,
              (const std::vector<double> &input, std::vector<double> &output),
              (override, const));
  MOCK_METHOD(std::vector<double>, GetCTIdentities, (), (override, const));
};

TEST_CASE("Check CheckImplementation", "[origin]")
{
  using namespace BSMPT;
  using ::testing::AtLeast;
  using ::testing::Return;
  SetLogger({"--logginglevel::complete=true"});
  Mock_Class_Potential_Origin mock_point;

  EXPECT_CALL(mock_point, addLegendCT()).Times(AtLeast(0));
  EXPECT_CALL(mock_point, addLegendVEV()).Times(1);
  EXPECT_CALL(mock_point, addLegendTemp()).Times(1);
  EXPECT_CALL(mock_point, addLegendTripleCouplings()).Times(1);

  ModelTests::CheckImplementation(mock_point, Minimizer::WhichMinimizerDefault);
}
