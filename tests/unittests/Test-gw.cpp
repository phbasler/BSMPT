// SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Approx = Catch::Approx;

#include <BSMPT/bounce_solution/action_calculation.h>
#include <BSMPT/gravitational_waves/gw.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/minimum_tracer/minimum_tracer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/ModelTestfunctions.h>
#include <BSMPT/transition_tracer/transition_tracer.h>
#include <BSMPT/utility/Logger.h> // for Logger Class
#include <fstream>

TEST_CASE("Solve bounce equation with analytical derivative", "[gw]")
{
  // Tests bounce solver with analytical derivative
  using namespace BSMPT;
  std::function<double(std::vector<double>)> V = [&](std::vector<double> x)
  {
    double c  = 5;
    double fx = 0;
    double fy = 80;

    double r1 = x[0] * x[0] + c * x[1] * x[1];
    double r2 = c * pow(x[0] - 1, 2) + pow(x[1] - 1, 2);
    double r3 = fx * (0.25 * pow(x[0], 4) - pow(x[0], 3) / 3.);
    r3 += fy * (0.25 * pow(x[1], 4) - pow(x[1], 3) / 3.);

    return (r1 * r2 + r3);
  };

  std::function<std::vector<double>(std::vector<double>)> dV =
      [&](std::vector<double> l0)
  {
    int dim = 2;
    std::vector<double> result(dim);
    result = {2 * l0[0] * (5 * pow(-1 + l0[0], 2) + pow(-1 + l0[1], 2)) +
                  10 * (-1 + l0[0]) * (pow(l0[0], 2) + 5 * pow(l0[1], 2)),
              10 * (5 * pow(-1 + l0[0], 2) + pow(-1 + l0[1], 2)) * l0[1] +
                  2 * (-1 + l0[1]) * (pow(l0[0], 2) + 5 * pow(l0[1], 2)) +
                  80 * (-1. * pow(l0[1], 2) + 1. * pow(l0[1], 3))};
    return result;
  };

  std::vector<double> FalseVacuum = {0, 0};
  std::vector<double> TrueVacuum  = {1, 1};

  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, dV, 0, 6);
  bc.CalculateAction();

  REQUIRE(bc.action == Approx(4.5011952256).epsilon(5e-2));
}

TEST_CASE("Solve bounce equation with numerical derivative", "[gw]")
{
  // Tests bounce solver with numerical derivative
  using namespace BSMPT;

  std::function<double(std::vector<double>)> V = [&](std::vector<double> x)
  {
    double c  = 5;
    double fx = 0;
    double fy = 80;

    double r1 = x[0] * x[0] + c * x[1] * x[1];
    double r2 = c * pow(x[0] - 1, 2) + pow(x[1] - 1, 2);
    double r3 = fx * (0.25 * pow(x[0], 4) - pow(x[0], 3) / 3.);
    r3 += fy * (0.25 * pow(x[1], 4) - pow(x[1], 3) / 3.);

    return (r1 * r2 + r3);
  };

  std::vector<double> FalseVacuum = {0, 0};
  std::vector<double> TrueVacuum  = {1, 1};

  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);
  bc.CalculateAction();

  REQUIRE(bc.action == Approx(4.5011952256).epsilon(5e-2));
}

TEST_CASE("Solve bounce equation with numerical derivative and displaced "
          "potential in VEV space",
          "[gw]")
{
  // Tests bounce solver with numerical derivative and displaced potential in
  // the VEV space
  using namespace BSMPT;

  double disp_0 = 1;
  double disp_1 = 1;

  std::function<double(std::vector<double>)> V = [&](std::vector<double> x)
  {
    x[0] += disp_0;
    x[1] += disp_1;

    double c  = 5;
    double fx = 0;
    double fy = 80;

    double r1 = x[0] * x[0] + c * x[1] * x[1];
    double r2 = c * pow(x[0] - 1, 2) + pow(x[1] - 1, 2);
    double r3 = fx * (0.25 * pow(x[0], 4) - pow(x[0], 3) / 3.);
    r3 += fy * (0.25 * pow(x[1], 4) - pow(x[1], 3) / 3.);

    return (r1 * r2 + r3);
  };

  std::vector<double> FalseVacuum = {0, 0};
  std::vector<double> TrueVacuum  = {1, 1};

  FalseVacuum[0] -= disp_0;
  FalseVacuum[1] -= disp_1;

  TrueVacuum[0] -= disp_0;
  TrueVacuum[1] -= disp_1;

  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);

  bc.CalculateAction();

  REQUIRE(bc.action == Approx(4.5011952256).epsilon(5e-2));
}

TEST_CASE("Solve bounce equation with numerical derivative and displaced "
          "potential in energy",
          "[gw]")
{
  // Tests bounce solver with numerical derivative and displaced potential in
  // the VEV space
  using namespace BSMPT;

  double disp_V = 1;

  std::function<double(std::vector<double>)> V = [&](std::vector<double> x)
  {
    double c  = 5;
    double fx = 0;
    double fy = 80;

    double r1 = x[0] * x[0] + c * x[1] * x[1];
    double r2 = c * pow(x[0] - 1, 2) + pow(x[1] - 1, 2);
    double r3 = fx * (0.25 * pow(x[0], 4) - pow(x[0], 3) / 3.);
    r3 += fy * (0.25 * pow(x[1], 4) - pow(x[1], 3) / 3.);

    return ((r1 * r2 + r3) - disp_V);
  };

  std::vector<double> FalseVacuum = {0, 0};
  std::vector<double> TrueVacuum  = {1, 1};

  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);

  bc.CalculateAction();

  REQUIRE(bc.action == Approx(4.5011952256).epsilon(5e-2));
}

TEST_CASE("Solve bounce equation with numerical derivative and displaced "
          "potential in VEV space and in energy",
          "[gw]")
{
  // Tests bounce solver with numerical derivative and displaced potential in
  // the VEV space
  using namespace BSMPT;

  double disp_0 = 1;
  double disp_1 = 1;
  double disp_V = 1;

  std::function<double(std::vector<double>)> V = [&](std::vector<double> x)
  {
    x[0] += disp_0;
    x[1] += disp_1;

    double c  = 5;
    double fx = 0;
    double fy = 80;

    double r1 = x[0] * x[0] + c * x[1] * x[1];
    double r2 = c * pow(x[0] - 1, 2) + pow(x[1] - 1, 2);
    double r3 = fx * (0.25 * pow(x[0], 4) - pow(x[0], 3) / 3.);
    r3 += fy * (0.25 * pow(x[1], 4) - pow(x[1], 3) / 3.);

    return (r1 * r2 + r3) - disp_V;
  };

  std::vector<double> FalseVacuum = {0, 0};
  std::vector<double> TrueVacuum  = {1, 1};

  FalseVacuum[0] -= disp_0;
  FalseVacuum[1] -= disp_1;

  TrueVacuum[0] -= disp_0;
  TrueVacuum[1] -= disp_1;

  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);

  bc.CalculateAction();

  REQUIRE(bc.action == Approx(4.5011952256).epsilon(5e-2));
}

////////////////////

TEST_CASE("Solve bounce equation with numerical derivative thin walled", "[gw]")
{
  // Tests bounce solver with numerical derivative
  using namespace BSMPT;

  std::function<double(std::vector<double>)> V = [&](std::vector<double> x)
  {
    double c  = 5;
    double fx = 0;
    double fy = 2.;

    double r1 = x[0] * x[0] + c * x[1] * x[1];
    double r2 = c * pow(x[0] - 1, 2) + pow(x[1] - 1, 2);
    double r3 = fx * (0.25 * pow(x[0], 4) - pow(x[0], 3) / 3.);
    r3 += fy * (0.25 * pow(x[1], 4) - pow(x[1], 3) / 3.);

    return (r1 * r2 + r3);
  };

  std::vector<double> FalseVacuum = {0, 0};
  std::vector<double> TrueVacuum  = {1, 1};

  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);

  bc.CalculateAction();

  REQUIRE(bc.action == Approx(1946.3823079011).epsilon(5e-2));
}

TEST_CASE("Checking phase tracking for SM", "[gw]")
{
  const std::vector<double> example_point_SM{
      /* muSq = */ -7823.7540500000005,
      /* lambda = */ 0.12905349405143487};

  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::SM, SMConstants);
  modelPointer->initModel(example_point_SM);

  std::shared_ptr<MinimumTracer> MinTracer(
      new MinimumTracer(modelPointer, Minimizer::WhichMinimizerDefault, false));
  Vacuum vac(0, 300, MinTracer, modelPointer, -1, 10, true);

  REQUIRE(vac.PhasesList.size() == 2);
}

TEST_CASE("Checking phase tracking for BP1", "[gw]")
{
  const std::vector<double> example_point_R2HDM{
      /* lambda_1 = */ 6.9309437685026,
      /* lambda_2 = */ 0.26305141403285998,
      /* lambda_3 = */ 1.2865950045595,
      /* lambda_4 = */ 4.7721306931875001,
      /* lambda_5 = */ 4.7275722046239004,
      /* m_{12}^2 = */ 18933.440789693999,
      /* tan(beta) = */ 16.577896825227999,
      /* Yukawa Type = */ 1};

  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);

  std::shared_ptr<MinimumTracer> MinTracer(
      new MinimumTracer(modelPointer, Minimizer::WhichMinimizerDefault, false));
  Vacuum vac(0, 300, MinTracer, modelPointer, -1, 10, true);

  REQUIRE(vac.PhasesList.size() == 2);
}

TEST_CASE("Checking phase tracking for BP2", "[gw]")
{
  const std::vector<double> example_point_R2HDM{
      /* lambda_1 = */ 6.8467197321288999,
      /* lambda_2 = */ 0.25889890874393001,
      /* lambda_3 = */ 1.4661775278406,
      /* lambda_4 = */ 4.4975594646125998,
      /* lambda_5 = */ 4.4503516057569996,
      /* m_{12}^2 = */ 6629.9728323804002,
      /* tan(beta) = */ 45.319927369307997,
      /* Yukawa Type = */ 1};

  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);

  std::shared_ptr<MinimumTracer> MinTracer(
      new MinimumTracer(modelPointer, Minimizer::WhichMinimizerDefault, false));
  Vacuum vac(0, 300, MinTracer, modelPointer, -1, 10, true);

  REQUIRE(vac.PhasesList.size() == 2);
}

TEST_CASE("Checking phase tracking for BP3", "[gw]")
{
  const std::vector<double> example_point_CXSM{/* v = */ 245.34120667410863,
                                               /* vs = */ 0,
                                               /* va = */ 0,
                                               /* msq = */ -15650,
                                               /* lambda = */ 0.52,
                                               /* delta2 = */ 0.55,
                                               /* b2 = */ -8859,
                                               /* d2 = */ 0.5,
                                               /* Reb1 = */ 0,
                                               /* Imb1 = */ 0,
                                               /* Rea1 = */ 0,
                                               /* Ima1 = */ 0};

  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM, SMConstants);
  modelPointer->initModel(example_point_CXSM);

  std::shared_ptr<MinimumTracer> MinTracer(
      new MinimumTracer(modelPointer, Minimizer::WhichMinimizerDefault, false));
  Vacuum vac(0, 300, MinTracer, modelPointer, -1, 10, true);

  REQUIRE(vac.PhasesList.size() == 2);
}
