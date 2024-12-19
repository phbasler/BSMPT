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
#include <gsl/gsl_sf_expint.h>

TEST_CASE("Test GetStatusEWSR", "[gw]")
{
  // Tests bounce solver with analytical derivative
  using namespace BSMPT;

  MinimumTracer MinTracer;
  REQUIRE(MinTracer.GetStatusEWSR(3) == StatusEWSR::EWSymRes);
  REQUIRE(MinTracer.GetStatusEWSR(2) == StatusEWSR::EWSymNonRes);
  REQUIRE(MinTracer.GetStatusEWSR(1) == StatusEWSR::FlatRegion);
  REQUIRE(MinTracer.GetStatusEWSR(0) == StatusEWSR::Failure);
  REQUIRE(MinTracer.GetStatusEWSR(-1) == StatusEWSR::NotBFB);
}

TEST_CASE("Test Create1DimGrid (point, k, low_value, high_value, nsteps)",
          "[gw]")
{
  using namespace BSMPT;

  auto grid                               = Create1DimGrid({1, 3}, 0, 0, 2, 4);
  std::vector<std::vector<double>> result = {
      {0, 3}, {0.5, 3}, {1, 3}, {1.5, 3}, {2.0, 3.0}};

  REQUIRE(grid == result);
}

TEST_CASE("Test Create1DimGrid (min_start,min_end,npoints)", "[gw] ")
{
  using namespace BSMPT;

  auto grid                               = Create1DimGrid({0, 3}, {2, 3}, 4);
  std::vector<std::vector<double>> result = {
      {0, 3}, {0.5, 3}, {1, 3}, {1.5, 3}, {2.0, 3.0}};

  REQUIRE(grid == result);
}

TEST_CASE("Test almost_the_same", "[gw]")
{
  using namespace BSMPT;

  REQUIRE(almost_the_same({0, 1}, {0, 0.991}, false, 0.01, 1e-5));
  REQUIRE(not almost_the_same({0, 1}, {0, 1.02}, false, 0.01, 1e-5));
  REQUIRE(not almost_the_same({0, 1}, {0, 0.991}, false, 0.01, 0));
}

TEST_CASE("Test I_alpha", "[gw]")
{
  // Tests bounce solver with analytical derivative
  using namespace BSMPT;

  BounceActionInt BACalc;

  REQUIRE(BACalc.BesselI(3, 1) == Approx(0.0221684249).epsilon(1e-8));
  REQUIRE(BACalc.BesselI(1, 3) == Approx(3.953370217).epsilon(1e-8));
  REQUIRE(BACalc.BesselI(1, 1.5) == Approx(0.9816664285779).epsilon(1e-8));
}

TEST_CASE("Test J_1", "[gw]")
{
  // Tests bounce solver with analytical derivative
  using namespace BSMPT;

  BounceActionInt BACalc;

  REQUIRE(BACalc.BesselJ(1) == Approx(0.440050585744).epsilon(1e-8));
  REQUIRE(BACalc.BesselJ(3) == Approx(0.3390589585259).epsilon(1e-8));
  REQUIRE(BACalc.BesselJ(1.5) == Approx(0.5579365079).epsilon(1e-8));
}

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

  REQUIRE(bc.Action == Approx(4.5011952256).epsilon(5e-2));
}

TEST_CASE(
    "Solve bounce equation with analytical derivative and Alpha = 3 (T = 0)",
    "[gw]")
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

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);
  bc.Alpha = 3;
  bc.CalculateAction();

  REQUIRE(bc.Action == Approx(13.3129767888).epsilon(5e-2));
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

  REQUIRE(bc.Action == Approx(4.5011952256).epsilon(5e-2));
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

  REQUIRE(bc.Action == Approx(4.5011952256).epsilon(5e-2));
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

  REQUIRE(bc.Action == Approx(4.5011952256).epsilon(5e-2));
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

  REQUIRE(bc.Action == Approx(4.5011952256).epsilon(5e-2));
}

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

  REQUIRE(bc.Action == Approx(1946.3823079011).epsilon(5e-2));
}

TEST_CASE("Catch if path if backwards on the bounce solver", "[gw]")
{
  // Tests bounce solver with numerical derivative
  using namespace BSMPT;
  std::function<double(std::vector<double>)> V = [&](std::vector<double> x)
  {
    double c  = 5;
    double fx = 0;
    double fy = 80.;

    double r1 = x[0] * x[0] + c * x[1] * x[1];
    double r2 = c * pow(x[0] - 1, 2) + pow(x[1] - 1, 2);
    double r3 = fx * (0.25 * pow(x[0], 4) - pow(x[0], 3) / 3.);
    r3 += fy * (0.25 * pow(x[1], 4) - pow(x[1], 3) / 3.);

    return (r1 * r2 + r3);
  };

  std::vector<double> TrueVacuum  = {0, 0};
  std::vector<double> FalseVacuum = {1, 1};

  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);

  bc.CalculateAction();

  REQUIRE(BSMPT::BounceActionInt::ActionStatus::BackwardsPropagationFailed ==
          bc.StateOfBounceActionInt);
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

TEST_CASE("Checking phase tracking for BP1 - Mode auto", "[gw]")
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

TEST_CASE("Checking phase tracking for BP1 - Mode 0", "[gw]")
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
  Vacuum vac(0, 300, MinTracer, modelPointer, 0, 10, true);

  REQUIRE(vac.PhasesList.size() == 2);
}

TEST_CASE("Checking phase tracking for BP1 - Mode 1", "[gw]")
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
  Vacuum vac(0, 300, MinTracer, modelPointer, 1, 10, true);

  REQUIRE(vac.PhasesList.size() == 2);
}
TEST_CASE("Checking phase tracking for BP1 - Mode 2", "[gw]")
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
  Vacuum vac(0, 300, MinTracer, modelPointer, 2, 10, true);

  REQUIRE(vac.PhasesList.size() == 2);
}

TEST_CASE("Checking phase tracking for BP2 - Mode auto", "[gw]")
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

TEST_CASE("Checking phase tracking for BP2 - Mode 0", "[gw]")
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
  Vacuum vac(0, 300, MinTracer, modelPointer, 0, 10, true);

  REQUIRE(vac.PhasesList.size() == 2);
}

TEST_CASE("Checking phase tracking for BP2 - Mode 1", "[gw]")
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
  Vacuum vac(0, 300, MinTracer, modelPointer, 0, 10, true);

  REQUIRE(vac.PhasesList.size() == 2);
}

TEST_CASE("Checking phase tracking for BP2 - Mode 2", "[gw]")
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
  Vacuum vac(0, 300, MinTracer, modelPointer, 0, 10, true);

  REQUIRE(vac.PhasesList.size() == 2);
}

TEST_CASE("Checking phase tracking for BP3 with Mode 0", "[gw]")
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
  Vacuum vac(0, 300, MinTracer, modelPointer, 0, 10, true);

  REQUIRE(vac.PhasesList.size() == 2);
}

TEST_CASE("Checking phase tracking for SM with Mode 1", "[gw]")
{
  const std::vector<double> example_point_SM{
      /* muSq = */ -7823.7540500000005,
      /* lambda = */ 0.12905349405143487};

  using namespace BSMPT;
  // Check the ASCIIPlotter
  SetLogger({"--logginglevel::mintracerdetailed=true"});
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::SM, SMConstants);
  modelPointer->initModel(example_point_SM);

  std::shared_ptr<MinimumTracer> MinTracer(
      new MinimumTracer(modelPointer, Minimizer::WhichMinimizerDefault, false));
  Vacuum vac(0, 300, MinTracer, modelPointer, 1, 10, true);
  SetLogger({"--logginglevel::mintracerdetailed=false"});

  REQUIRE(vac.PhasesList.size() == 2);
}

TEST_CASE("Checking phase tracking for SM with Mode 2", "[gw]")
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
  Vacuum vac(0, 300, MinTracer, modelPointer, 2, 10, true);

  REQUIRE(vac.PhasesList.size() == 2);
}

TEST_CASE("Check maximal thermal mass squared over temperature ratio")
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

  user_input input;
  input.modelPointer   = modelPointer;
  input.gw_calculation = true;
  TransitionTracer trans(input);
  trans.ListBounceSolution.at(0).CalculatePercolationTemp();

  auto output = trans.output_store;

  REQUIRE(0.781639 == Approx(trans.CheckMassRatio(
                                 input,
                                 output.vec_trans_data.at(0).crit_false_vev,
                                 output.vec_trans_data.at(0).crit_temp.value()))
                          .epsilon(1e-2));
  REQUIRE(0.781639 ==
          Approx(trans.CheckMassRatio(
                     input,
                     output.vec_trans_data.at(0).nucl_approx_false_vev,
                     output.vec_trans_data.at(0).nucl_approx_temp.value()))
              .epsilon(1e-2));
  REQUIRE(0.781639 == Approx(trans.CheckMassRatio(
                                 input,
                                 output.vec_trans_data.at(0).nucl_false_vev,
                                 output.vec_trans_data.at(0).nucl_temp.value()))
                          .epsilon(1e-2));
  REQUIRE(0.781639 == Approx(trans.CheckMassRatio(
                                 input,
                                 output.vec_trans_data.at(0).perc_false_vev,
                                 output.vec_trans_data.at(0).perc_temp.value()))
                          .epsilon(1e-2));
  REQUIRE(0.781639 ==
          Approx(trans.CheckMassRatio(
                     input,
                     output.vec_trans_data.at(0).compl_false_vev,
                     output.vec_trans_data.at(0).compl_temp.value()))
              .epsilon(1e-2));
}

TEST_CASE("Checking phase tracking and GW for BP3", "[gw]")
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

  user_input input;
  input.modelPointer   = modelPointer;
  input.gw_calculation = true;
  TransitionTracer trans(input);
  trans.ListBounceSolution.at(0).CalculatePercolationTemp();

  auto output = trans.output_store;

  REQUIRE(126.0223716 ==
          Approx(output.vec_trans_data.at(0).crit_temp.value()).epsilon(1e-2));
  REQUIRE(121.0869527 ==
          Approx(output.vec_trans_data.at(0).nucl_approx_temp.value())
              .epsilon(1e-2));
  REQUIRE(121.212833 ==
          Approx(output.vec_trans_data.at(0).nucl_temp.value()).epsilon(1e-2));
  REQUIRE(120.7670659 ==
          Approx(output.vec_trans_data.at(0).perc_temp.value()).epsilon(1e-2));
  REQUIRE(120.7267244 ==
          Approx(output.vec_trans_data.at(0).compl_temp.value()).epsilon(1e-2));

  REQUIRE(0.00537281 ==
          Approx(output.vec_gw_data.at(0).alpha.value()).epsilon(1e-2));
  REQUIRE(7658.8931 ==
          Approx(output.vec_gw_data.at(0).beta_over_H.value()).epsilon(1e-2));
  REQUIRE(4.40964e-05 ==
          Approx(output.vec_gw_data.at(0).K_sw.value()).epsilon(1e-2));
  REQUIRE(4.40964e-06 ==
          Approx(output.vec_gw_data.at(0).K_turb.value()).epsilon(1e-2));
  REQUIRE(0.0884755 ==
          Approx(output.vec_gw_data.at(0).fpeak_sw.value()).epsilon(1e-2));
  REQUIRE(0.269136 ==
          Approx(output.vec_gw_data.at(0).fpeak_turb.value()).epsilon(1e-2));
  REQUIRE(1.70812e-20 ==
          Approx(output.vec_gw_data.at(0).h2Omega_sw.value()).epsilon(1e-2));
  REQUIRE(3.69052e-16 ==
          Approx(output.vec_gw_data.at(0).h2Omega_turb.value()).epsilon(1e-2));
  REQUIRE(1.23742e-09 ==
          Approx(output.vec_gw_data.at(0).SNR_sw.value()).epsilon(5e-2));
  REQUIRE(1.28789e-20 ==
          Approx(output.vec_gw_data.at(0).SNR_turb.value()).epsilon(5e-2));
  REQUIRE(1.23742e-09 ==
          Approx(output.vec_gw_data.at(0).SNR.value()).epsilon(5e-2));

  // Check different vwalls
  trans.ListBounceSolution.at(0).UserDefined_vwall = -1;
  trans.ListBounceSolution.at(0).CalculatePTStrength();
  REQUIRE(0.374931042806113 ==
          Approx(trans.ListBounceSolution.at(0).vwall).epsilon(1e-2));

  trans.ListBounceSolution.at(0).UserDefined_vwall = -2;
  trans.ListBounceSolution.at(0).CalculatePTStrength();
  REQUIRE(0.303761086384691 ==
          Approx(trans.ListBounceSolution.at(0).vwall).epsilon(1e-2));
}

TEST_CASE("Checking phase tracking and GW for BP3 (low sample)", "[gw]")
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

  user_input input;
  input.modelPointer                        = modelPointer;
  input.gw_calculation                      = true;
  input.number_of_initial_scan_temperatures = 6;
  TransitionTracer trans(input);
  trans.ListBounceSolution.at(0).CalculatePercolationTemp();

  auto output = trans.output_store;

  REQUIRE(126.0223716 ==
          Approx(output.vec_trans_data.at(0).crit_temp.value()).epsilon(1e-2));
  REQUIRE(121.0869527 ==
          Approx(output.vec_trans_data.at(0).nucl_approx_temp.value())
              .epsilon(1e-2));
  REQUIRE(121.212833 ==
          Approx(output.vec_trans_data.at(0).nucl_temp.value()).epsilon(1e-2));
  REQUIRE(120.7670659 ==
          Approx(output.vec_trans_data.at(0).perc_temp.value()).epsilon(1e-2));
  REQUIRE(120.7267244 ==
          Approx(output.vec_trans_data.at(0).compl_temp.value()).epsilon(1e-2));

  REQUIRE(0.00537281 ==
          Approx(output.vec_gw_data.at(0).alpha.value()).epsilon(1e-2));
  REQUIRE(7658.8931 ==
          Approx(output.vec_gw_data.at(0).beta_over_H.value()).epsilon(1e-2));
  REQUIRE(4.40964e-05 ==
          Approx(output.vec_gw_data.at(0).K_sw.value()).epsilon(1e-2));
  REQUIRE(4.40964e-06 ==
          Approx(output.vec_gw_data.at(0).K_turb.value()).epsilon(1e-2));
  REQUIRE(0.0884755 ==
          Approx(output.vec_gw_data.at(0).fpeak_sw.value()).epsilon(1e-2));
  REQUIRE(0.269136 ==
          Approx(output.vec_gw_data.at(0).fpeak_turb.value()).epsilon(1e-2));
  REQUIRE(1.70812e-20 ==
          Approx(output.vec_gw_data.at(0).h2Omega_sw.value()).epsilon(1e-2));
  REQUIRE(3.69052e-16 ==
          Approx(output.vec_gw_data.at(0).h2Omega_turb.value()).epsilon(1e-2));
  REQUIRE(1.23742e-09 ==
          Approx(output.vec_gw_data.at(0).SNR_sw.value()).epsilon(5e-2));
  REQUIRE(1.28789e-20 ==
          Approx(output.vec_gw_data.at(0).SNR_turb.value()).epsilon(5e-2));
  REQUIRE(1.23742e-09 ==
          Approx(output.vec_gw_data.at(0).SNR.value()).epsilon(5e-2));

  // Check different vwalls
  trans.ListBounceSolution.at(0).UserDefined_vwall = -1;
  trans.ListBounceSolution.at(0).CalculatePTStrength();
  REQUIRE(0.374931042806113 ==
          Approx(trans.ListBounceSolution.at(0).vwall).epsilon(1e-2));

  trans.ListBounceSolution.at(0).UserDefined_vwall = -2;
  trans.ListBounceSolution.at(0).CalculatePTStrength();
  REQUIRE(0.303761086384691 ==
          Approx(trans.ListBounceSolution.at(0).vwall).epsilon(1e-2));
}

TEST_CASE(
    "Checking phase tracking and GW for BP3 (low sample) (supposed to fail)",
    "[gw]")
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

  user_input input;
  input.modelPointer   = modelPointer;
  input.gw_calculation = true;

  input.maxpathintegrations = 1;
  TransitionTracer trans(input);
  auto output = trans.output_store;
  REQUIRE(trans.ListBounceSolution.at(0).status_bounce_sol ==
          StatusGW::Failure);
  REQUIRE(126.0223716 ==
          Approx(output.vec_trans_data.at(0).crit_temp.value()).epsilon(1e-2));
  REQUIRE(not output.vec_trans_data.at(0).nucl_approx_temp.has_value());
  REQUIRE(not output.vec_trans_data.at(0).nucl_temp.has_value());
  REQUIRE(not output.vec_trans_data.at(0).perc_temp.has_value());
  REQUIRE(not output.vec_trans_data.at(0).compl_temp.has_value());
  REQUIRE(not output.vec_gw_data.at(0).alpha.has_value());
  REQUIRE(not output.vec_gw_data.at(0).beta_over_H.has_value());
  REQUIRE(not output.vec_gw_data.at(0).K_sw.has_value());
  REQUIRE(not output.vec_gw_data.at(0).K_turb.has_value());
  REQUIRE(not output.vec_gw_data.at(0).fpeak_sw.has_value());
  REQUIRE(not output.vec_gw_data.at(0).fpeak_turb.has_value());
  REQUIRE(not output.vec_gw_data.at(0).h2Omega_sw.has_value());
  REQUIRE(not output.vec_gw_data.at(0).h2Omega_turb.has_value());
  REQUIRE(not output.vec_gw_data.at(0).SNR_sw.has_value());
  REQUIRE(not output.vec_gw_data.at(0).SNR_turb.has_value());
  REQUIRE(not output.vec_gw_data.at(0).SNR.has_value());
}

TEST_CASE("Test for SO(3)", "[gw]")
{
  const std::vector<double> example_point_CXSM{/* v = */ 10,
                                               /* vs = */ 10,
                                               /* va = */ 10,
                                               /* msq = */ -100,
                                               /* lambda = */ 0,
                                               /* delta2 = */ 0,
                                               /* b2 = */ -100,
                                               /* d2 = */ 0,
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
  MinTracer->FindFlatDirections();
  REQUIRE(MinTracer->flat_3D_dirs.size() == 1);
}

TEST_CASE("Test for EW symmetry restoration BP1", "[gw]")
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
  REQUIRE(MinTracer->IsThereEWSymmetryRestoration() == -1);
}

TEST_CASE("Test for EW symmetry restoration BP2", "[gw]")
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
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);

  std::shared_ptr<MinimumTracer> MinTracer(
      new MinimumTracer(modelPointer, Minimizer::WhichMinimizerDefault, false));
  REQUIRE(MinTracer->IsThereEWSymmetryRestoration() == -1);
}

TEST_CASE("Test for EW symmetry restoration BP3", "[gw]")
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
  REQUIRE(MinTracer->IsThereEWSymmetryRestoration() == 3);
}

TEST_CASE("Check percolation temperature for point of issue #173", "[gw]")
{
  const std::vector<double> example_point_R2HDM_machine_dep_temp_failure{
      /* lambda_1 = */ 1.939019085953685,
      /* lambda_2 = */ 2.368448996569993,
      /* lambda_3 = */ 7.6135777364421875,
      /* lambda_4 = */ -5.450539233086667,
      /* lambda_5 = */ -3.9705820766209374,
      /* m_{12}^2 = */ 162500.89873320065,
      /* tan(beta) = */ 1.2197198842183743,
      /* Yukawa Type = */ 1};

  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM_machine_dep_temp_failure);

  std::shared_ptr<MinimumTracer> MinTracer(
      new MinimumTracer(modelPointer, Minimizer::WhichMinimizerDefault, false));

  user_input input;
  input.modelPointer   = modelPointer;
  input.gw_calculation = true;
  input.T_high         = 1000;
  TransitionTracer trans(input);
  trans.ListBounceSolution.at(0).CalculatePercolationTemp();

  auto output = trans.output_store;

  REQUIRE(314.566 ==
          Approx(output.vec_trans_data.at(1).perc_temp.value()).epsilon(1e-2));
}

TEST_CASE("Test string conversion of enums", "[gw]")
{
  using namespace BSMPT;
  std::stringstream ss;
  REQUIRE_NOTHROW(
      ss << StatusNLOStability::NotSet << StatusNLOStability::Off
         << StatusNLOStability::NoNLOStability << StatusNLOStability::Success
         << StatusEWSR::EWSymNonRes << StatusEWSR::EWSymRes
         << StatusEWSR::Failure << StatusEWSR::FlatRegion << StatusEWSR::NotBFB
         << StatusEWSR::NotSet << StatusEWSR::Off << StatusTracing::Failure
         << StatusTracing::NoCoverage << StatusTracing::NoGlobMinCoverage
         << StatusTracing::NoGlobMinCoverage
         << StatusTracing::NoMinsAtBoundaries << StatusTracing::NotSet
         << StatusTracing::Success << StatusCoexPair::NoCoexPairs
         << StatusCoexPair::NotSet << StatusCoexPair::Success
         << StatusCrit::Failure << StatusCrit::FalseLower << StatusCrit::NotSet
         << StatusCrit::Success << StatusCrit::TrueLower << StatusGW::Failure
         << StatusGW::NotSet << StatusGW::Success << StatusTemperature::NaN
         << StatusTemperature::NotMet << StatusTemperature::NotSet
         << StatusTemperature::Success);
}

TEST_CASE(
    "Espinosa-Konstandin - Example A: Polynomial Vt - Phi0 = 0.99 (thin wall)",
    "[gw]")
{
  // Espinosa-Konstandin examples from arXiv:2312.12360
  using namespace BSMPT;
  double phi0                                  = 0.99;
  std::function<double(std::vector<double>)> V = [&](std::vector<double> x)
  {
    if (x[0] == 0) return 0.;
    if (x[0] == 1) return -1.;
    return x[0] * x[0] *
           (2 * x[0] - 3 +
            pow(1 - x[0], 2) * log((pow(1. - x[0], 2) * phi0 * phi0) /
                                   (pow(1. - phi0, 2) * x[0] * x[0])));
  };

  std::vector<double> FalseVacuum = {0};
  std::vector<double> TrueVacuum  = {1};

  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);
  bc.Alpha = 3;
  bc.CalculateAction();

  REQUIRE(bc.Action ==
          Approx(-pow(M_PI, 2) / 3. * (phi0 + Li2(phi0 / (phi0 - 1))))
              .epsilon(1e-3));
}

TEST_CASE(
    "Espinosa-Konstandin - Example A: Polynomial Vt - Phi0 = 0.5 (thick wall)",
    "[gw]")
{
  // Espinosa-Konstandin examples from arXiv:2312.12360
  using namespace BSMPT;
  double phi0                                  = 0.5;
  std::function<double(std::vector<double>)> V = [&](std::vector<double> x)
  {
    if (x[0] == 0) return 0.;
    if (x[0] == 1) return -1.;
    return x[0] * x[0] *
           (2 * x[0] - 3 +
            pow(1 - x[0], 2) * log((pow(1. - x[0], 2) * phi0 * phi0) /
                                   (pow(1. - phi0, 2) * x[0] * x[0])));
  };

  std::vector<double> FalseVacuum = {0};
  std::vector<double> TrueVacuum  = {1 - 0.17};

  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);
  bc.Alpha = 3;
  bc.CalculateAction();

  REQUIRE(bc.Action ==
          Approx(-pow(M_PI, 2) / 3. * (phi0 + Li2(phi0 / (phi0 - 1))))
              .epsilon(1e-3));
}

TEST_CASE("Espinosa-Konstandin - Example B: Trigonometric tunneling potential "
          "- Phi0 = 1.4 (thick wall)",
          "[gw]")
{
  // Espinosa-Konstandin examples from arXiv:2312.12360
  using namespace BSMPT;
  double phi0                                  = 1.4;
  std::function<double(std::vector<double>)> V = [&](std::vector<double> x)
  {
    if (x[0] == 0) return 0.;
    x[0] = abs(x[0]);
    return pow(sin(x[0]), 2) * (-1. + 2. / 3. * pow(cos(x[0]), 2) *
                                          log(tan(phi0) / abs(tan(x[0]))));
  };

  std::vector<double> FalseVacuum = {0};
  std::vector<double> TrueVacuum  = {1.54751};

  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);
  bc.Alpha = 3;
  bc.CalculateAction();

  REQUIRE(bc.Action ==
          Approx(pow(M_PI, 2) / 4. *
                 (pow(M_PI, 2) / 2. + 6 * pow(log(1. / tan(phi0)), 2) +
                  3 * Li2(-pow(1. / tan(phi0), 2))))
              .epsilon(1e-3));
}

TEST_CASE("Espinosa-Konstandin - Example B: Trigonometric tunneling potential "
          "- Phi0 = 1.54 (thin wall)",
          "[gw]")
{
  // Espinosa-Konstandin examples from arXiv:2312.12360
  using namespace BSMPT;
  double phi0                                  = 1.54;
  std::function<double(std::vector<double>)> V = [&](std::vector<double> x)
  {
    if (x[0] == 0) return 0.;
    x[0] = abs(x[0]);
    return pow(sin(x[0]), 2) * (-1. + 2. / 3. * pow(cos(x[0]), 2) *
                                          log(tan(phi0) / abs(tan(x[0]))));
  };

  std::vector<double> FalseVacuum = {0};
  std::vector<double> TrueVacuum  = {1.5666275};

  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);
  bc.Alpha = 3;
  bc.CalculateAction();

  REQUIRE(bc.Action ==
          Approx(pow(M_PI, 2) / 4. *
                 (pow(M_PI, 2) / 2. + 6 * pow(log(1. / tan(phi0)), 2) +
                  3 * Li2(-pow(1. / tan(phi0), 2))))
              .epsilon(1e-3));
}

TEST_CASE("Espinosa-Konstandin - Example C: Finite mass "
          "- Phi0 = 0.8",
          "[gw]")
{
  // Espinosa-Konstandin examples from arXiv:2312.12360
  using namespace BSMPT;
  double phi0    = 0.8;
  double cut_off = 1.3;

  std::function<double(std::vector<double>)> V = [&](std::vector<double> x)
  {
    if (x[0] < 0) x[0] = -x[0];
    if (x[0] > cut_off) x[0] = 2 * cut_off - x[0];
    if (x[0] == 0) x[0] = 1e-100;

    return pow(x[0], 2) / (-0.5 + log(x[0])) +
           (8. * pow(x[0], 2) * pow(-1 + log(x[0]), 2) *
            (2. * pow(log(x[0]), 2) - 2. * pow(log(phi0), 2) +
             log((-1 + log(x[0])) / (-1 + log(phi0))))) /
               (3. * pow(1. - 2. * log(x[0]), 4));
  };

  std::vector<double> FalseVacuum = {0};
  std::vector<double> TrueVacuum  = {cut_off};

  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);
  bc.Alpha = 3;
  bc.CalculateAction();

  REQUIRE(bc.Action == Approx(3.236736977787636).epsilon(1e-3));
}

TEST_CASE(
    "Espinosa-Konstandin - Example D: Derivative of the tunneling potential "
    "- Phi0 = 0.8",
    "[gw]")
{
  // Espinosa-Konstandin examples from arXiv:2312.12360
  using namespace BSMPT;
  double phi0    = 0.8;
  double cut_off = 0.95;

  std::function<double(std::vector<double>)> V = [&](std::vector<double> x)
  {
    if (x[0] < 0) x[0] = -x[0];
    if (x[0] > cut_off) x[0] = 2 * cut_off - x[0];

    return gsl_sf_expint_Ei(2 * log(x[0] + 1e-100)) +
           1. / 6. * pow(x[0], 2) * (1. - pow(log(phi0) / log(x[0]), 2));
  };

  std::vector<double> FalseVacuum = {0};
  std::vector<double> TrueVacuum  = {cut_off};

  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);
  bc.Alpha = 3;
  bc.CalculateAction();

  REQUIRE(bc.Action == Approx(5.1968779132).epsilon(1e-3));
}

TEST_CASE("Espinosa-Konstandin - Example E: Finite mass with thin-wall limit "
          "- Phi0 = 0.8",
          "[gw]")
{
  // Espinosa-Konstandin examples from arXiv:2312.12360
  using namespace BSMPT;
  double phi0    = 0.8;
  double cut_off = 0.9;

  std::function<double(std::vector<double>)> V = [&](std::vector<double> x)
  {
    if (x[0] < 0) x[0] = -x[0];
    if (x[0] > cut_off) x[0] = 2 * cut_off - x[0];
    if (x[0] == 0) x[0] = 1e-100;

    return exp(2) * gsl_sf_expint_Ei(-2 + 2 * log(x[0])) -
           exp(3) * gsl_sf_expint_Ei(-3 + 3 * log(x[0])) +
           (pow(-1 + x[0], 2) * pow(x[0], 2) *
            (pow(log(x[0]), 2) - pow(log(phi0), 2) +
             2 * log(((-1 + x[0]) * phi0) / (x[0] * (-1 + phi0))) +
             2 * Li2(1 - x[0]) - 2 * Li2(1 - phi0))) /
               (6. * pow(-1 + log(x[0]), 2));
  };

  std::vector<double> FalseVacuum = {0};
  std::vector<double> TrueVacuum  = {cut_off};

  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);
  bc.Alpha = 3;
  bc.CalculateAction();

  REQUIRE(bc.Action == Approx(37.878999).epsilon(1e-3));
}

TEST_CASE("Espinosa-Konstandin - Two-Field Examples", "[gw]")
{
  // Espinosa-Konstandin examples from arXiv:2312.12360
  using namespace BSMPT;
  double phi0  = 0.999;
  double alpha = 1 / 2.;

  std::function<double(std::vector<double>)> V = [&](std::vector<double> x)
  {
    x[0] = x[0] + 1e-100;
    x[1] = x[1] + 1e-100;
    return -4.5 * (-2. + x[0]) * pow(1. + pow(x[0], 2) - 1. * pow(x[1], 2), 2) +
           pow(EllipIntSecond(asinh(x[0])), 2) *
               (-3 + 2 * EllipIntSecond(asinh(x[0])) +
                pow(-1 + EllipIntSecond(asinh(x[0])), 2) *
                    log((pow(phi0, 2) *
                         pow(-1 + EllipIntSecond(asinh(x[0])), 2)) /
                        (pow(-1 + phi0, 2) *
                         pow(EllipIntSecond(asinh(x[0])), 2)))) +
           ((-sqrt(1 + pow(x[0], 2)) + x[1]) *
            (-1 + EllipIntSecond(asinh(x[0]))) * EllipIntSecond(asinh(x[0])) *
            pow(1. / cosh(x[0] / alpha), 3) *
            (4 * pow(EllipIntSecond(asinh(x[0])), 2) *
                 log((pow(phi0, 2) * pow(-1 + EllipIntSecond(asinh(x[0])), 2)) /
                     (pow(-1 + phi0, 2) *
                      pow(EllipIntSecond(asinh(x[0])), 2))) -
             4 * alpha * pow(cosh(x[0] / alpha), 2) *
                 (-4 + log((pow(phi0, 2) *
                            pow(-1 + EllipIntSecond(asinh(x[0])), 2)) /
                           (pow(-1 + phi0, 2) *
                            pow(EllipIntSecond(asinh(x[0])), 2)))) *
                 sinh(x[0] / alpha) +
             2 * EllipIntSecond(asinh(x[0])) *
                 log((pow(phi0, 2) * pow(-1 + EllipIntSecond(asinh(x[0])), 2)) /
                     (pow(-1 + phi0, 2) *
                      pow(EllipIntSecond(asinh(x[0])), 2))) *
                 (-2 + alpha * sinh(x[0] / alpha) +
                  alpha * sinh((3 * x[0]) / alpha)))) /
               (2. * alpha);
  };

  std::vector<double> FalseVacuum = {0., 1.};
  std::vector<double> TrueVacuum  = {0.9181398979435043, 1.3575643168526};

  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);
  bc.Alpha = 3;
  bc.CalculateAction();

  REQUIRE(bc.Action == Approx(80.59).epsilon(2e-2));
}

TEST_CASE("Espinosa-Konstandin - Three-Field Examples - rho = 1/6", "[gw]")
{
  // Espinosa-Konstandin examples from arXiv:2312.12360
  using namespace BSMPT;
  double phi0  = 0.999;
  double alpha = 1 / 2.;
  double rho   = 1 / 6.;

  std::function<double(std::vector<double>)> V = [&](std::vector<double> x)
  {
    x[0] = x[0] + 1e-100;
    x[1] = x[1] + 1e-100;
    x[2] = x[2] + 1e-100;
    if (x[0] == 0.5) x[0] -= 0.0000001;
    return 25 * pow(x[2] - rho +
                        rho * cos((sqrt(1 - pow(alpha, 2)) * x[0]) /
                                  (alpha * rho)),
                    2) +
           (pow(x[0], 2) *
            (alpha * (-3 * alpha + 2 * x[0]) +
             pow(alpha - x[0], 2) * log((pow(phi0, 2) * pow(alpha - x[0], 2)) /
                                        (pow(-1 + phi0, 2) * pow(x[0], 2))))) /
               pow(alpha, 4) +
           25 * pow(x[1] - rho * sin((sqrt(1 - pow(alpha, 2)) * x[0]) /
                                     (alpha * rho)),
                    2) +
           (2 * (alpha - x[0]) * x[0] *
            (x[1] -
             rho * sin((sqrt(1 - pow(alpha, 2)) * x[0]) / (alpha * rho))) *
            (alpha * sqrt(1 - pow(alpha, 2)) * rho *
                 cos((sqrt(1 - pow(alpha, 2)) * x[0]) / (alpha * rho)) *
                 (-4 * alpha + (alpha - 2 * x[0]) *
                                   log((pow(phi0, 2) * pow(alpha - x[0], 2)) /
                                       (pow(-1 + phi0, 2) * pow(x[0], 2)))) +
             (-1 + pow(alpha, 2)) * (alpha - x[0]) * x[0] *
                 log((pow(phi0, 2) * pow(alpha - x[0], 2)) /
                     (pow(-1 + phi0, 2) * pow(x[0], 2))) *
                 sin((sqrt(1 - pow(alpha, 2)) * x[0]) / (alpha * rho)))) /
               (pow(alpha, 4) * rho) -
           (2 * x[0] * (-alpha + x[0]) *
            (x[2] - rho +
             rho * cos((sqrt(1 - pow(alpha, 2)) * x[0]) / (alpha * rho))) *
            (-((-1 + pow(alpha, 2)) * (alpha - x[0]) * x[0] *
               cos((sqrt(1 - pow(alpha, 2)) * x[0]) / (alpha * rho)) *
               log((pow(phi0, 2) * pow(alpha - x[0], 2)) /
                   (pow(-1 + phi0, 2) * pow(x[0], 2)))) +
             alpha * sqrt(1 - pow(alpha, 2)) * rho *
                 (-4 * alpha + (alpha - 2 * x[0]) *
                                   log((pow(phi0, 2) * pow(alpha - x[0], 2)) /
                                       (pow(-1 + phi0, 2) * pow(x[0], 2)))) *
                 sin((sqrt(1 - pow(alpha, 2)) * x[0]) / (alpha * rho)))) /
               (pow(alpha, 4) * rho);
  };

  std::vector<double> FalseVacuum = {0., 0., 0.};
  std::vector<double> TrueVacuum  = {0.500068, -0.147487, 0.0890438};

  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);
  bc.Alpha = 3;
  bc.CalculateAction();

  REQUIRE(bc.Action == Approx(80.59).epsilon(2e-2));
}

TEST_CASE("Espinosa-Konstandin - Three-Field Examples - rho = 1/8", "[gw]")
{
  // Espinosa-Konstandin examples from arXiv:2312.12360
  using namespace BSMPT;
  double phi0  = 0.999;
  double alpha = 1 / 2.;
  double rho   = 1 / 8.;

  std::function<double(std::vector<double>)> V = [&](std::vector<double> x)
  {
    x[0] = x[0] + 1e-100;
    x[1] = x[1] + 1e-100;
    x[2] = x[2] + 1e-100;
    if (x[0] == 0.5) x[0] -= 0.0000001;
    return 25 * pow(x[2] - rho +
                        rho * cos((sqrt(1 - pow(alpha, 2)) * x[0]) /
                                  (alpha * rho)),
                    2) +
           (pow(x[0], 2) *
            (alpha * (-3 * alpha + 2 * x[0]) +
             pow(alpha - x[0], 2) * log((pow(phi0, 2) * pow(alpha - x[0], 2)) /
                                        (pow(-1 + phi0, 2) * pow(x[0], 2))))) /
               pow(alpha, 4) +
           25 * pow(x[1] - rho * sin((sqrt(1 - pow(alpha, 2)) * x[0]) /
                                     (alpha * rho)),
                    2) +
           (2 * (alpha - x[0]) * x[0] *
            (x[1] -
             rho * sin((sqrt(1 - pow(alpha, 2)) * x[0]) / (alpha * rho))) *
            (alpha * sqrt(1 - pow(alpha, 2)) * rho *
                 cos((sqrt(1 - pow(alpha, 2)) * x[0]) / (alpha * rho)) *
                 (-4 * alpha + (alpha - 2 * x[0]) *
                                   log((pow(phi0, 2) * pow(alpha - x[0], 2)) /
                                       (pow(-1 + phi0, 2) * pow(x[0], 2)))) +
             (-1 + pow(alpha, 2)) * (alpha - x[0]) * x[0] *
                 log((pow(phi0, 2) * pow(alpha - x[0], 2)) /
                     (pow(-1 + phi0, 2) * pow(x[0], 2))) *
                 sin((sqrt(1 - pow(alpha, 2)) * x[0]) / (alpha * rho)))) /
               (pow(alpha, 4) * rho) -
           (2 * x[0] * (-alpha + x[0]) *
            (x[2] - rho +
             rho * cos((sqrt(1 - pow(alpha, 2)) * x[0]) / (alpha * rho))) *
            (-((-1 + pow(alpha, 2)) * (alpha - x[0]) * x[0] *
               cos((sqrt(1 - pow(alpha, 2)) * x[0]) / (alpha * rho)) *
               log((pow(phi0, 2) * pow(alpha - x[0], 2)) /
                   (pow(-1 + phi0, 2) * pow(x[0], 2)))) +
             alpha * sqrt(1 - pow(alpha, 2)) * rho *
                 (-4 * alpha + (alpha - 2 * x[0]) *
                                   log((pow(phi0, 2) * pow(alpha - x[0], 2)) /
                                       (pow(-1 + phi0, 2) * pow(x[0], 2)))) *
                 sin((sqrt(1 - pow(alpha, 2)) * x[0]) / (alpha * rho)))) /
               (pow(alpha, 4) * rho);
  };

  std::vector<double> FalseVacuum = {0., 0., 0.};
  std::vector<double> TrueVacuum  = {0.500068, 0.0752454, 0.0251845};

  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);
  bc.Alpha = 3;
  bc.CalculateAction();

  REQUIRE(bc.Action == Approx(55.6).epsilon(2e-2));
}
