// SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using Approx = Catch::Approx;

#include <BSMPT/bounce_solution/action_calculation.h>
#include <BSMPT/gravitational_waves/gw.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/minimum_tracer/minimum_tracer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/modeltests/ModelTestfunctions.h>
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

TEST_CASE("Check calculation of Chapman-Jouget velocity", "[gw]")
{
  const std::vector<double> example_point_SM{
      /* muSq = */ -7823.7540500000005,
      /* lambda = */ 0.12905349405143487};

  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::SM, SMConstants);
  modelPointer->initModel(example_point_SM);

  user_input input;
  input.modelPointer   = modelPointer;
  input.gw_calculation = true;
  TransitionTracer trans(input);
  trans.ListBounceSolution.at(0).SetAndCalculateGWParameters(
      TransitionTemperature::Percolation);

  REQUIRE(0.583037 ==
          Approx(trans.ListBounceSolution.at(0).GetChapmanJougetVelocity())
              .epsilon(1e-2));
}

TEST_CASE("Check calculation of reheating temperature", "[gw]")
{
  const std::vector<double> example_point_SM{
      /* muSq = */ -7823.7540500000005,
      /* lambda = */ 0.12905349405143487};

  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::SM, SMConstants);
  modelPointer->initModel(example_point_SM);

  user_input input;
  input.modelPointer   = modelPointer;
  input.gw_calculation = true;
  TransitionTracer trans(input);
  trans.ListBounceSolution.at(0).SetAndCalculateGWParameters(
      TransitionTemperature::Percolation);

  REQUIRE(
      159.07847220175202 ==
      Approx(trans.ListBounceSolution.at(0).GetReheatingTemp()).epsilon(1e-2));
}

TEST_CASE("Check maximal thermal mass squared over temperature ratio", "[gw]")
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

TEST_CASE("Check gstar implementation", "[gw]")
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

  BounceSolution BASolution(modelPointer);
  BASolution.InitializeGstarProfile();
  REQUIRE(108.75 == Approx(BASolution.CalcGstarPureRad()).epsilon(1e-2));
  REQUIRE(108.75 == Approx(BASolution.GetGstar()).epsilon(1e-2));
  REQUIRE(108.75 == Approx(BASolution.GetGstar(1000)).epsilon(1e-2));
  REQUIRE(104.65 == Approx(BASolution.GetGstar(100)).epsilon(1e-2));
  REQUIRE(86.72 == Approx(BASolution.GetGstar(10)).epsilon(1e-2));
  REQUIRE(76.52 == Approx(BASolution.GetGstar(1)).epsilon(1e-2));
  REQUIRE(10.63 == Approx(BASolution.GetGstar(1e-3)).epsilon(1e-2));
  REQUIRE(3.366 == Approx(BASolution.GetGstar(1e-6)).epsilon(1e-2));
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
  trans.ListBounceSolution.at(0).SetAndCalculateGWParameters(
      TransitionTemperature::Percolation);

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

  REQUIRE(0.0056735067 ==
          Approx(output.vec_gw_data.at(0).alpha.value()).epsilon(1e-2));
  REQUIRE(7658.8931 ==
          Approx(output.vec_gw_data.at(0).beta_over_H.value()).epsilon(1e-2));
  REQUIRE(0.0084807773 ==
          Approx(output.vec_gw_data.at(0).kappa_sw.value()).epsilon(1e-2));

  REQUIRE(0.01699469177 ==
          Approx(output.vec_gw_data.at(0).fb_col.value()).epsilon(1e-2));
  REQUIRE(0. ==
          Approx(output.vec_gw_data.at(0).omegab_col.value()).epsilon(1e-2));

  REQUIRE(0.007441950626 ==
          Approx(output.vec_gw_data.at(0).f1_sw.value()).epsilon(1e-2));
  REQUIRE(0.04661408636 ==
          Approx(output.vec_gw_data.at(0).f2_sw.value()).epsilon(1e-2));
  REQUIRE(9.481654892e-20 ==
          Approx(output.vec_gw_data.at(0).omega_2_sw.value()).epsilon(1e-2));

  REQUIRE(2.729908945e-05 ==
          Approx(output.vec_gw_data.at(0).f1_turb.value()).epsilon(1e-2));
  REQUIRE(0.08186145689 ==
          Approx(output.vec_gw_data.at(0).f2_turb.value()).epsilon(1e-2));
  REQUIRE(3.71764571e-25 ==
          Approx(output.vec_gw_data.at(0).omega_2_turb.value()).epsilon(1e-2));

  REQUIRE(0 == Approx(output.vec_gw_data.at(0).SNR_col.value()).epsilon(5e-2));
  REQUIRE(9.681903249e-08 ==
          Approx(output.vec_gw_data.at(0).SNR_sw.value()).epsilon(5e-2));
  REQUIRE(1.210536738e-12 ==
          Approx(output.vec_gw_data.at(0).SNR_turb.value()).epsilon(5e-2));
  REQUIRE(9.682005317e-08 ==
          Approx(output.vec_gw_data.at(0).SNR.value()).epsilon(5e-2));

  // Check different vwalls
  trans.ListBounceSolution.at(0).UserDefined_vwall = -1;
  trans.ListBounceSolution.at(0).SetAndCalculateGWParameters(
      TransitionTemperature::Percolation);
  REQUIRE(0.374931042806113 ==
          Approx(trans.ListBounceSolution.at(0).vwall).epsilon(1e-2));

  trans.ListBounceSolution.at(0).UserDefined_vwall = -2;
  trans.ListBounceSolution.at(0).SetAndCalculateGWParameters(
      TransitionTemperature::Percolation);
  REQUIRE(0.6539310662 ==
          Approx(trans.ListBounceSolution.at(0).vwall).epsilon(1e-2));
}

TEST_CASE("Checking phase tracking and GW for BP3 (low sample) and not "
          "transition temperature not set",
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
  input.modelPointer                        = modelPointer;
  input.gw_calculation                      = true;
  input.number_of_initial_scan_temperatures = 6;
  TransitionTracer trans(input);

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

  REQUIRE(0.0056735067 ==
          Approx(output.vec_gw_data.at(0).alpha.value()).epsilon(1e-2));
  REQUIRE(7658.8931 ==
          Approx(output.vec_gw_data.at(0).beta_over_H.value()).epsilon(1e-2));
  REQUIRE(0.0084807773 ==
          Approx(output.vec_gw_data.at(0).kappa_sw.value()).epsilon(1e-2));

  REQUIRE(0.01699469177 ==
          Approx(output.vec_gw_data.at(0).fb_col.value()).epsilon(1e-2));
  REQUIRE(0. ==
          Approx(output.vec_gw_data.at(0).omegab_col.value()).epsilon(1e-2));

  REQUIRE(0.007441950626 ==
          Approx(output.vec_gw_data.at(0).f1_sw.value()).epsilon(1e-2));
  REQUIRE(0.04661408636 ==
          Approx(output.vec_gw_data.at(0).f2_sw.value()).epsilon(1e-2));
  REQUIRE(9.481654892e-20 ==
          Approx(output.vec_gw_data.at(0).omega_2_sw.value()).epsilon(1e-2));

  REQUIRE(2.729908945e-05 ==
          Approx(output.vec_gw_data.at(0).f1_turb.value()).epsilon(1e-2));
  REQUIRE(0.08186145689 ==
          Approx(output.vec_gw_data.at(0).f2_turb.value()).epsilon(1e-2));
  REQUIRE(3.71764571e-25 ==
          Approx(output.vec_gw_data.at(0).omega_2_turb.value()).epsilon(1e-2));

  REQUIRE(0 == Approx(output.vec_gw_data.at(0).SNR_col.value()).epsilon(5e-2));
  REQUIRE(9.681903249e-08 ==
          Approx(output.vec_gw_data.at(0).SNR_sw.value()).epsilon(5e-2));
  REQUIRE(1.210536738e-12 ==
          Approx(output.vec_gw_data.at(0).SNR_turb.value()).epsilon(5e-2));
  REQUIRE(9.682005317e-08 ==
          Approx(output.vec_gw_data.at(0).SNR.value()).epsilon(5e-2));

  // Check different vwalls
  trans.ListBounceSolution.at(0).UserDefined_vwall = -1;
  trans.ListBounceSolution.at(0).SetAndCalculateGWParameters(
      TransitionTemperature::Percolation);
  REQUIRE(0.374931042806113 ==
          Approx(trans.ListBounceSolution.at(0).vwall).epsilon(1e-2));

  trans.ListBounceSolution.at(0).UserDefined_vwall = -2;
  trans.ListBounceSolution.at(0).SetAndCalculateGWParameters(
      TransitionTemperature::Percolation);
  REQUIRE(0.6539310662 ==
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
  REQUIRE(not output.vec_gw_data.at(0).kappa_col.has_value());
  REQUIRE(not output.vec_gw_data.at(0).kappa_sw.has_value());
  REQUIRE(not output.vec_gw_data.at(0).fb_col.has_value());
  REQUIRE(not output.vec_gw_data.at(0).omegab_col.has_value());
  REQUIRE(not output.vec_gw_data.at(0).f1_sw.has_value());
  REQUIRE(not output.vec_gw_data.at(0).f2_sw.has_value());
  REQUIRE(not output.vec_gw_data.at(0).omega_2_sw.has_value());
  REQUIRE(not output.vec_gw_data.at(0).f1_turb.has_value());
  REQUIRE(not output.vec_gw_data.at(0).f2_turb.has_value());
  REQUIRE(not output.vec_gw_data.at(0).omega_2_turb.has_value());
  REQUIRE(not output.vec_gw_data.at(0).SNR_col.has_value());
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
  modelPointer->initModel(example_point_CXSM, false);
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
  modelPointer->initModel(example_point_R2HDM, false);
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

TEST_CASE("Test dydv", "[gw]")
{
  using namespace BSMPT::kappa;
  double v   = 0.1;
  double y[] = {0.6, 1};
  double dydv[2];
  double cs2 = 0.57735;
  dfdv(v, y, dydv, &cs2);

  REQUIRE(dydv[0] == Approx(-1.4525696202948652).epsilon(1e-10));
  REQUIRE(dydv[1] == Approx(1.4678979234569798).epsilon(1e-10));

  v = 0.7;
  dfdv(v, y, dydv, &cs2);
  REQUIRE(dydv[0] == Approx(-0.4623000345532002).epsilon(1e-10));
  REQUIRE(dydv[1] == Approx(-0.9236144743536612).epsilon(1e-10));

  y[0] = 0.3;
  y[1] = 0.4;
  v    = 0.2;
  dfdv(v, y, dydv, &cs2);
  REQUIRE(dydv[0] == Approx(-0.7199796242092907).epsilon(1e-10));
  REQUIRE(dydv[1] == Approx(0.12110157868520084).epsilon(1e-10));
}

TEST_CASE("Test getKandWow", "[gw]")
{
  using namespace BSMPT::kappa;
  std::vector<std::vector<double>> vprofile;
  auto [K, wow] = getKandWow(0.6, 0.4, 1 / sqrt(3.), vprofile);

  REQUIRE(K == Approx(0.1943013335557658).epsilon(1e-4));
  REQUIRE(wow == Approx(1.938501886826841).epsilon(1e-4));
}

TEST_CASE("Test kappa_sw", "[gw]")
{
  auto [cs2b, cs2s, al, vw, expected] = GENERATE(
      std::tuple{0.37735026918962583,
                 0.37735026918962583,
                 0.0001,
                 0.05,
                 7.3179378666420855e-06},
      std::tuple{0.37735026918962583,
                 0.37735026918962583,
                 0.0001,
                 0.35,
                 0.0002008455736394683},
      std::tuple{0.37735026918962583,
                 0.37735026918962583,
                 0.0001,
                 0.65,
                 0.001985541586762711},
      std::tuple{0.37735026918962583,
                 0.37735026918962583,
                 0.0001,
                 0.95,
                 0.00021420942862449336},
      std::tuple{0.37735026918962583,
                 0.37735026918962583,
                 0.0031622776601683794,
                 0.05,
                 0.00023266326269007358},
      std::tuple{0.37735026918962583,
                 0.37735026918962583,
                 0.0031622776601683794,
                 0.35,
                 0.0062915254891997945},
      std::tuple{0.37735026918962583,
                 0.37735026918962583,
                 0.0031622776601683794,
                 0.65,
                 0.11046950045060759},
      std::tuple{0.37735026918962583,
                 0.37735026918962583,
                 0.0031622776601683794,
                 0.95,
                 0.006728525079006433},
      std::tuple{0.37735026918962583,
                 0.37735026918962583,
                 0.1,
                 0.05,
                 0.007326900011342696},
      std::tuple{0.37735026918962583,
                 0.37735026918962583,
                 0.1,
                 0.35,
                 0.1583633860049138},
      std::tuple{0.37735026918962583,
                 0.37735026918962583,
                 0.1,
                 0.65,
                 0.4832524716454113},
      std::tuple{0.37735026918962583,
                 0.37735026918962583,
                 0.1,
                 0.95,
                 0.17634385707782466},
      std::tuple{0.37735026918962583,
                 0.37735026918962583,
                 3.1622776601683795,
                 0.05,
                 0},
      std::tuple{0.37735026918962583,
                 0.37735026918962583,
                 3.1622776601683795,
                 0.35,
                 0},
      std::tuple{0.37735026918962583,
                 0.37735026918962583,
                 3.1622776601683795,
                 0.65,
                 0},
      std::tuple{0.37735026918962583,
                 0.37735026918962583,
                 3.1622776601683795,
                 0.95,
                 0.9727036578377455},
      std::tuple{0.37735026918962583,
                 0.5106836025229592,
                 0.0001,
                 0.05,
                 7.540003466141407e-06},
      std::tuple{0.37735026918962583,
                 0.5106836025229592,
                 0.0001,
                 0.35,
                 0.00018923147050079953},
      std::tuple{0.37735026918962583,
                 0.5106836025229592,
                 0.0001,
                 0.65,
                 0.001985541586762711},
      std::tuple{0.37735026918962583,
                 0.5106836025229592,
                 0.0001,
                 0.95,
                 0.00021420942862449336},
      std::tuple{0.37735026918962583,
                 0.5106836025229592,
                 0.0031622776601683794,
                 0.05,
                 0.00023817340171740378},
      std::tuple{0.37735026918962583,
                 0.5106836025229592,
                 0.0031622776601683794,
                 0.35,
                 0.005945714385861376},
      std::tuple{0.37735026918962583,
                 0.5106836025229592,
                 0.0031622776601683794,
                 0.65,
                 0.03607648814742543},
      std::tuple{0.37735026918962583,
                 0.5106836025229592,
                 0.0031622776601683794,
                 0.95,
                 0.006728525079006433},
      std::tuple{0.37735026918962583,
                 0.5106836025229592,
                 0.1,
                 0.05,
                 0.007501379181425136},
      std::tuple{0.37735026918962583,
                 0.5106836025229592,
                 0.1,
                 0.35,
                 0.15655928946830946},
      std::tuple{0.37735026918962583,
                 0.5106836025229592,
                 0.1,
                 0.65,
                 0.3954665998695123},
      std::tuple{0.37735026918962583,
                 0.5106836025229592,
                 0.1,
                 0.95,
                 0.17634385707782466},
      std::tuple{
          0.37735026918962583, 0.5106836025229592, 3.1622776601683795, 0.05, 0},
      std::tuple{
          0.37735026918962583, 0.5106836025229592, 3.1622776601683795, 0.35, 0},
      std::tuple{
          0.37735026918962583, 0.5106836025229592, 3.1622776601683795, 0.65, 0},
      std::tuple{0.37735026918962583,
                 0.5106836025229592,
                 3.1622776601683795,
                 0.95,
                 0.970642029722613},
      std::tuple{0.37735026918962583,
                 0.6440169358562926,
                 0.0001,
                 0.05,
                 7.692929307907293e-06},
      std::tuple{0.37735026918962583,
                 0.6440169358562926,
                 0.0001,
                 0.35,
                 0.00018545963266591788},
      std::tuple{0.37735026918962583,
                 0.6440169358562926,
                 0.0001,
                 0.65,
                 0.001985541586762711},
      std::tuple{0.37735026918962583,
                 0.6440169358562926,
                 0.0001,
                 0.95,
                 0.00021420942862449336},
      std::tuple{0.37735026918962583,
                 0.6440169358562926,
                 0.0031622776601683794,
                 0.05,
                 0.00024211617232036147},
      std::tuple{0.37735026918962583,
                 0.6440169358562926,
                 0.0031622776601683794,
                 0.35,
                 0.005833736573734788},
      std::tuple{0.37735026918962583,
                 0.6440169358562926,
                 0.0031622776601683794,
                 0.65,
                 0.03316056910212684},
      std::tuple{0.37735026918962583,
                 0.6440169358562926,
                 0.0031622776601683794,
                 0.95,
                 0.006728525079006433},
      std::tuple{0.37735026918962583,
                 0.6440169358562926,
                 0.1,
                 0.05,
                 0.007624504332890382},
      std::tuple{0.37735026918962583,
                 0.6440169358562926,
                 0.1,
                 0.35,
                 0.15694210231703978},
      std::tuple{0.37735026918962583,
                 0.6440169358562926,
                 0.1,
                 0.65,
                 0.3436082030718621},
      std::tuple{0.37735026918962583,
                 0.6440169358562926,
                 0.1,
                 0.95,
                 0.17634385707782466},
      std::tuple{
          0.37735026918962583, 0.6440169358562926, 3.1622776601683795, 0.05, 0},
      std::tuple{
          0.37735026918962583, 0.6440169358562926, 3.1622776601683795, 0.35, 0},
      std::tuple{
          0.37735026918962583, 0.6440169358562926, 3.1622776601683795, 0.65, 0},
      std::tuple{0.37735026918962583,
                 0.6440169358562926,
                 3.1622776601683795,
                 0.95,
                 0.9687133749239248},
      std::tuple{0.37735026918962583,
                 0.7773502691896259,
                 0.0001,
                 0.05,
                 7.707876457599111e-06},
      std::tuple{0.37735026918962583,
                 0.7773502691896259,
                 0.0001,
                 0.35,
                 0.00018455369243278248},
      std::tuple{0.37735026918962583,
                 0.7773502691896259,
                 0.0001,
                 0.65,
                 0.001985541586762711},
      std::tuple{0.37735026918962583,
                 0.7773502691896259,
                 0.0001,
                 0.95,
                 0.00021420942862449336},
      std::tuple{0.37735026918962583,
                 0.7773502691896259,
                 0.0031622776601683794,
                 0.05,
                 0.00024499815032913154},
      std::tuple{0.37735026918962583,
                 0.7773502691896259,
                 0.0031622776601683794,
                 0.35,
                 0.005802184299445743},
      std::tuple{0.37735026918962583,
                 0.7773502691896259,
                 0.0031622776601683794,
                 0.65,
                 0.03240315684256059},
      std::tuple{0.37735026918962583,
                 0.7773502691896259,
                 0.0031622776601683794,
                 0.95,
                 0.006728525079006433},
      std::tuple{0.37735026918962583,
                 0.7773502691896259,
                 0.1,
                 0.05,
                 0.007717618473164975},
      std::tuple{0.37735026918962583,
                 0.7773502691896259,
                 0.1,
                 0.35,
                 0.1580415225935992},
      std::tuple{0.37735026918962583,
                 0.7773502691896259,
                 0.1,
                 0.65,
                 0.3150514998815157},
      std::tuple{0.37735026918962583,
                 0.7773502691896259,
                 0.1,
                 0.95,
                 0.17634385707782466},
      std::tuple{
          0.37735026918962583, 0.7773502691896259, 3.1622776601683795, 0.05, 0},
      std::tuple{
          0.37735026918962583, 0.7773502691896259, 3.1622776601683795, 0.35, 0},
      std::tuple{
          0.37735026918962583, 0.7773502691896259, 3.1622776601683795, 0.65, 0},
      std::tuple{0.37735026918962583,
                 0.7773502691896259,
                 3.1622776601683795,
                 0.95,
                 0.9668918501603003},
      std::tuple{0.5106836025229592,
                 0.37735026918962583,
                 0.0001,
                 0.05,
                 7.411159383997479e-06},
      std::tuple{0.5106836025229592,
                 0.37735026918962583,
                 0.0001,
                 0.35,
                 0.00021535515537291937},
      std::tuple{0.5106836025229592,
                 0.37735026918962583,
                 0.0001,
                 0.65,
                 0.21355709597248396},
      std::tuple{0.5106836025229592,
                 0.37735026918962583,
                 0.0001,
                 0.95,
                 0.00048331431171071257},
      std::tuple{0.5106836025229592,
                 0.37735026918962583,
                 0.0031622776601683794,
                 0.05,
                 0.0002334094314291312},
      std::tuple{0.5106836025229592,
                 0.37735026918962583,
                 0.0031622776601683794,
                 0.35,
                 0.006751762441268846},
      std::tuple{0.5106836025229592,
                 0.37735026918962583,
                 0.0031622776601683794,
                 0.65,
                 0.25693735830010916},
      std::tuple{0.5106836025229592,
                 0.37735026918962583,
                 0.0031622776601683794,
                 0.95,
                 0.01509800636192268},
      std::tuple{0.5106836025229592,
                 0.37735026918962583,
                 0.1,
                 0.05,
                 0.0073523723090345204},
      std::tuple{0.5106836025229592,
                 0.37735026918962583,
                 0.1,
                 0.35,
                 0.17225172320523233},
      std::tuple{0.5106836025229592,
                 0.37735026918962583,
                 0.1,
                 0.65,
                 0.5731438341546534},
      std::tuple{0.5106836025229592,
                 0.37735026918962583,
                 0.1,
                 0.95,
                 0.349203204340969},
      std::tuple{
          0.5106836025229592, 0.37735026918962583, 3.1622776601683795, 0.05, 0},
      std::tuple{
          0.5106836025229592, 0.37735026918962583, 3.1622776601683795, 0.35, 0},
      std::tuple{
          0.5106836025229592, 0.37735026918962583, 3.1622776601683795, 0.65, 0},
      std::tuple{0.5106836025229592,
                 0.37735026918962583,
                 3.1622776601683795,
                 0.95,
                 1.2516692504124665},
      std::tuple{0.5106836025229592,
                 0.5106836025229592,
                 0.0001,
                 0.05,
                 7.540003466141407e-06},
      std::tuple{0.5106836025229592,
                 0.5106836025229592,
                 0.0001,
                 0.35,
                 0.00020412867981207998},
      std::tuple{0.5106836025229592,
                 0.5106836025229592,
                 0.0001,
                 0.65,
                 0.0016405613872688411},
      std::tuple{0.5106836025229592,
                 0.5106836025229592,
                 0.0001,
                 0.95,
                 0.00048331431171071257},
      std::tuple{0.5106836025229592,
                 0.5106836025229592,
                 0.0031622776601683794,
                 0.05,
                 0.00023893967203783564},
      std::tuple{0.5106836025229592,
                 0.5106836025229592,
                 0.0031622776601683794,
                 0.35,
                 0.006412043572008043},
      std::tuple{0.5106836025229592,
                 0.5106836025229592,
                 0.0031622776601683794,
                 0.65,
                 0.04575694072054808},
      std::tuple{0.5106836025229592,
                 0.5106836025229592,
                 0.0031622776601683794,
                 0.95,
                 0.01509800636192268},
      std::tuple{0.5106836025229592,
                 0.5106836025229592,
                 0.1,
                 0.05,
                 0.0075278738492524895},
      std::tuple{0.5106836025229592,
                 0.5106836025229592,
                 0.1,
                 0.35,
                 0.1708924639957514},
      std::tuple{0.5106836025229592,
                 0.5106836025229592,
                 0.1,
                 0.65,
                 0.4755894145419821},
      std::tuple{
          0.5106836025229592, 0.5106836025229592, 0.1, 0.95, 0.349203204340969},
      std::tuple{
          0.5106836025229592, 0.5106836025229592, 3.1622776601683795, 0.05, 0},
      std::tuple{
          0.5106836025229592, 0.5106836025229592, 3.1622776601683795, 0.35, 0},
      std::tuple{
          0.5106836025229592, 0.5106836025229592, 3.1622776601683795, 0.65, 0},
      std::tuple{0.5106836025229592,
                 0.5106836025229592,
                 3.1622776601683795,
                 0.95,
                 1.2491837681467124},
      std::tuple{0.5106836025229592,
                 0.6440169358562926,
                 0.0001,
                 0.05,
                 7.692929307907293e-06},
      std::tuple{0.5106836025229592,
                 0.6440169358562926,
                 0.0001,
                 0.35,
                 0.00020084139139559613},
      std::tuple{0.5106836025229592,
                 0.6440169358562926,
                 0.0001,
                 0.65,
                 0.0008542268840703054},
      std::tuple{0.5106836025229592,
                 0.6440169358562926,
                 0.0001,
                 0.95,
                 0.00048331431171071257},
      std::tuple{0.5106836025229592,
                 0.6440169358562926,
                 0.0031622776601683794,
                 0.05,
                 0.00024279891079803805},
      std::tuple{0.5106836025229592,
                 0.6440169358562926,
                 0.0031622776601683794,
                 0.35,
                 0.006314226722618959},
      std::tuple{0.5106836025229592,
                 0.6440169358562926,
                 0.0031622776601683794,
                 0.65,
                 0.025928478833256597},
      std::tuple{0.5106836025229592,
                 0.6440169358562926,
                 0.0031622776601683794,
                 0.95,
                 0.01509800636192268},
      std::tuple{0.5106836025229592,
                 0.6440169358562926,
                 0.1,
                 0.05,
                 0.007651743217987089},
      std::tuple{0.5106836025229592,
                 0.6440169358562926,
                 0.1,
                 0.35,
                 0.17178262559292462},
      std::tuple{0.5106836025229592,
                 0.6440169358562926,
                 0.1,
                 0.65,
                 0.41723240786672683},
      std::tuple{
          0.5106836025229592, 0.6440169358562926, 0.1, 0.95, 0.349203204340969},
      std::tuple{
          0.5106836025229592, 0.6440169358562926, 3.1622776601683795, 0.05, 0},
      std::tuple{
          0.5106836025229592, 0.6440169358562926, 3.1622776601683795, 0.35, 0},
      std::tuple{
          0.5106836025229592, 0.6440169358562926, 3.1622776601683795, 0.65, 0},
      std::tuple{0.5106836025229592,
                 0.6440169358562926,
                 3.1622776601683795,
                 0.95,
                 1.2466103418640628},
      std::tuple{0.5106836025229592,
                 0.7773502691896259,
                 0.0001,
                 0.05,
                 7.806693241205771e-06},
      std::tuple{0.5106836025229592,
                 0.7773502691896259,
                 0.0001,
                 0.35,
                 0.0002001547078822659},
      std::tuple{0.5106836025229592,
                 0.7773502691896259,
                 0.0001,
                 0.65,
                 0.0006609907513608744},
      std::tuple{0.5106836025229592,
                 0.7773502691896259,
                 0.0001,
                 0.95,
                 0.00048331431171071257},
      std::tuple{0.5106836025229592,
                 0.7773502691896259,
                 0.0031622776601683794,
                 0.05,
                 0.0002457888939963252},
      std::tuple{0.5106836025229592,
                 0.7773502691896259,
                 0.0031622776601683794,
                 0.35,
                 0.006298444498776185},
      std::tuple{0.5106836025229592,
                 0.7773502691896259,
                 0.0031622776601683794,
                 0.65,
                 0.020373559368388087},
      std::tuple{0.5106836025229592,
                 0.7773502691896259,
                 0.0031622776601683794,
                 0.95,
                 0.01509800636192268},
      std::tuple{0.5106836025229592,
                 0.7773502691896259,
                 0.1,
                 0.05,
                 0.00774540086331589},
      std::tuple{0.5106836025229592,
                 0.7773502691896259,
                 0.1,
                 0.35,
                 0.17337307239197647},
      std::tuple{0.5106836025229592,
                 0.7773502691896259,
                 0.1,
                 0.65,
                 0.3849851141570427},
      std::tuple{
          0.5106836025229592, 0.7773502691896259, 0.1, 0.95, 0.349203204340969},
      std::tuple{
          0.5106836025229592, 0.7773502691896259, 3.1622776601683795, 0.05, 0},
      std::tuple{
          0.5106836025229592, 0.7773502691896259, 3.1622776601683795, 0.35, 0},
      std::tuple{
          0.5106836025229592, 0.7773502691896259, 3.1622776601683795, 0.65, 0},
      std::tuple{0.5106836025229592,
                 0.7773502691896259,
                 3.1622776601683795,
                 0.95,
                 1.2440939368345787},
      std::tuple{0.6440169358562926,
                 0.37735026918962583,
                 0.0001,
                 0.05,
                 7.411159383997479e-06},
      std::tuple{0.6440169358562926,
                 0.37735026918962583,
                 0.0001,
                 0.35,
                 0.00022461569610901424},
      std::tuple{0.6440169358562926,
                 0.37735026918962583,
                 0.0001,
                 0.65,
                 0.2264581276563679},
      std::tuple{0.6440169358562926,
                 0.37735026918962583,
                 0.0001,
                 0.95,
                 0.0010832003949572397},
      std::tuple{0.6440169358562926,
                 0.37735026918962583,
                 0.0031622776601683794,
                 0.05,
                 0.00023378296386764998},
      std::tuple{0.6440169358562926,
                 0.37735026918962583,
                 0.0031622776601683794,
                 0.35,
                 0.007045529152416038},
      std::tuple{0.6440169358562926,
                 0.37735026918962583,
                 0.0031622776601683794,
                 0.65,
                 0.27524400608876515},
      std::tuple{0.6440169358562926,
                 0.37735026918962583,
                 0.0031622776601683794,
                 0.95,
                 0.033497529479630764},
      std::tuple{0.6440169358562926,
                 0.37735026918962583,
                 0.1,
                 0.05,
                 0.007367376780085396},
      std::tuple{0.6440169358562926,
                 0.37735026918962583,
                 0.1,
                 0.35,
                 0.18140381964529517},
      std::tuple{0.6440169358562926,
                 0.37735026918962583,
                 0.1,
                 0.65,
                 0.6328843746021567},
      std::tuple{0.6440169358562926,
                 0.37735026918962583,
                 0.1,
                 0.95,
                 0.6516798361681999},
      std::tuple{
          0.6440169358562926, 0.37735026918962583, 3.1622776601683795, 0.05, 0},
      std::tuple{
          0.6440169358562926, 0.37735026918962583, 3.1622776601683795, 0.35, 0},
      std::tuple{
          0.6440169358562926, 0.37735026918962583, 3.1622776601683795, 0.65, 0},
      std::tuple{0.6440169358562926,
                 0.37735026918962583,
                 3.1622776601683795,
                 0.95,
                 1.4744151097706777},
      std::tuple{0.6440169358562926,
                 0.5106836025229592,
                 0.0001,
                 0.05,
                 7.540003466141407e-06},
      std::tuple{0.6440169358562926,
                 0.5106836025229592,
                 0.0001,
                 0.35,
                 0.00021351032531940882},
      std::tuple{0.6440169358562926,
                 0.5106836025229592,
                 0.0001,
                 0.65,
                 0.0018062037169213605},
      std::tuple{0.6440169358562926,
                 0.5106836025229592,
                 0.0001,
                 0.95,
                 0.0010832003949572397},
      std::tuple{0.6440169358562926,
                 0.5106836025229592,
                 0.0031622776601683794,
                 0.05,
                 0.00023932330850725875},
      std::tuple{0.6440169358562926,
                 0.5106836025229592,
                 0.0031622776601683794,
                 0.35,
                 0.006711124762791142},
      std::tuple{0.6440169358562926,
                 0.5106836025229592,
                 0.0031622776601683794,
                 0.65,
                 0.05034578173466377},
      std::tuple{0.6440169358562926,
                 0.5106836025229592,
                 0.0031622776601683794,
                 0.95,
                 0.033497529479630764},
      std::tuple{0.6440169358562926,
                 0.5106836025229592,
                 0.1,
                 0.05,
                 0.007543541096817625},
      std::tuple{0.6440169358562926,
                 0.5106836025229592,
                 0.1,
                 0.35,
                 0.1803659060719546},
      std::tuple{0.6440169358562926,
                 0.5106836025229592,
                 0.1,
                 0.65,
                 0.5287300815907237},
      std::tuple{0.6440169358562926,
                 0.5106836025229592,
                 0.1,
                 0.95,
                 0.6516798361681999},
      std::tuple{
          0.6440169358562926, 0.5106836025229592, 3.1622776601683795, 0.05, 0},
      std::tuple{
          0.6440169358562926, 0.5106836025229592, 3.1622776601683795, 0.35, 0},
      std::tuple{
          0.6440169358562926, 0.5106836025229592, 3.1622776601683795, 0.65, 0},
      std::tuple{0.6440169358562926,
                 0.5106836025229592,
                 3.1622776601683795,
                 0.95,
                 1.4709593598880144},
      std::tuple{0.6440169358562926,
                 0.6440169358562926,
                 0.0001,
                 0.05,
                 7.692929307907293e-06},
      std::tuple{0.6440169358562926,
                 0.6440169358562926,
                 0.0001,
                 0.35,
                 0.00021061167167705277},
      std::tuple{0.6440169358562926,
                 0.6440169358562926,
                 0.0001,
                 0.65,
                 0.0009456342249578272},
      std::tuple{0.6440169358562926,
                 0.6440169358562926,
                 0.0001,
                 0.95,
                 0.0010832003949572397},
      std::tuple{0.6440169358562926,
                 0.6440169358562926,
                 0.0031622776601683794,
                 0.05,
                 0.00024328716972870837},
      std::tuple{0.6440169358562926,
                 0.6440169358562926,
                 0.0031622776601683794,
                 0.35,
                 0.00662410929054976},
      std::tuple{0.6440169358562926,
                 0.6440169358562926,
                 0.0031622776601683794,
                 0.65,
                 0.0287227378668273},
      std::tuple{0.6440169358562926,
                 0.6440169358562926,
                 0.0031622776601683794,
                 0.95,
                 0.033497529479630764},
      std::tuple{0.6440169358562926,
                 0.6440169358562926,
                 0.1,
                 0.05,
                 0.007667778068094371},
      std::tuple{0.6440169358562926,
                 0.6440169358562926,
                 0.1,
                 0.35,
                 0.18161578458489833},
      std::tuple{0.6440169358562926,
                 0.6440169358562926,
                 0.1,
                 0.65,
                 0.4661174138210306},
      std::tuple{0.6440169358562926,
                 0.6440169358562926,
                 0.1,
                 0.95,
                 0.6516798361681999},
      std::tuple{
          0.6440169358562926, 0.6440169358562926, 3.1622776601683795, 0.05, 0},
      std::tuple{
          0.6440169358562926, 0.6440169358562926, 3.1622776601683795, 0.35, 0},
      std::tuple{
          0.6440169358562926, 0.6440169358562926, 3.1622776601683795, 0.65, 0},
      std::tuple{0.6440169358562926,
                 0.6440169358562926,
                 3.1622776601683795,
                 0.95,
                 1.4681186958771253},
      std::tuple{0.6440169358562926,
                 0.7773502691896259,
                 0.0001,
                 0.05,
                 7.806693241205771e-06},
      std::tuple{0.6440169358562926,
                 0.7773502691896259,
                 0.0001,
                 0.35,
                 0.00021014521843789488},
      std::tuple{0.6440169358562926,
                 0.7773502691896259,
                 0.0001,
                 0.65,
                 0.0007352512140041711},
      std::tuple{0.6440169358562926,
                 0.7773502691896259,
                 0.0001,
                 0.95,
                 0.0010832003949572397},
      std::tuple{0.6440169358562926,
                 0.7773502691896259,
                 0.0031622776601683794,
                 0.05,
                 0.0002461847436252291},
      std::tuple{0.6440169358562926,
                 0.7773502691896259,
                 0.0031622776601683794,
                 0.35,
                 0.006619503935970726},
      std::tuple{0.6440169358562926,
                 0.7773502691896259,
                 0.0031622776601683794,
                 0.65,
                 0.022673459914223072},
      std::tuple{0.6440169358562926,
                 0.7773502691896259,
                 0.0031622776601683794,
                 0.95,
                 0.033497529479630764},
      std::tuple{0.6440169358562926,
                 0.7773502691896259,
                 0.1,
                 0.05,
                 0.0077618374741581795},
      std::tuple{0.6440169358562926,
                 0.7773502691896259,
                 0.1,
                 0.35,
                 0.18355303095597814},
      std::tuple{0.6440169358562926,
                 0.7773502691896259,
                 0.1,
                 0.65,
                 0.4315289844972973},
      std::tuple{0.6440169358562926,
                 0.7773502691896259,
                 0.1,
                 0.95,
                 0.6516798361681999},
      std::tuple{
          0.6440169358562926, 0.7773502691896259, 3.1622776601683795, 0.05, 0},
      std::tuple{
          0.6440169358562926, 0.7773502691896259, 3.1622776601683795, 0.35, 0},
      std::tuple{
          0.6440169358562926, 0.7773502691896259, 3.1622776601683795, 0.65, 0},
      std::tuple{0.6440169358562926,
                 0.7773502691896259,
                 3.1622776601683795,
                 0.95,
                 1.4652553264278354},
      std::tuple{0.7773502691896259,
                 0.37735026918962583,
                 0.0001,
                 0.05,
                 7.411159383997479e-06},
      std::tuple{0.7773502691896259,
                 0.37735026918962583,
                 0.0001,
                 0.35,
                 0.00023116136272582555},
      std::tuple{0.7773502691896259,
                 0.37735026918962583,
                 0.0001,
                 0.65,
                 0.23563606395424863},
      std::tuple{0.7773502691896259,
                 0.37735026918962583,
                 0.0001,
                 0.95,
                 0.0030546091145512624},
      std::tuple{0.7773502691896259,
                 0.37735026918962583,
                 0.0031622776601683794,
                 0.05,
                 0.00023406330922685052},
      std::tuple{0.7773502691896259,
                 0.37735026918962583,
                 0.0031622776601683794,
                 0.35,
                 0.007249139201052769},
      std::tuple{0.7773502691896259,
                 0.37735026918962583,
                 0.0031622776601683794,
                 0.65,
                 0.288806813335828},
      std::tuple{0.7773502691896259,
                 0.37735026918962583,
                 0.0031622776601683794,
                 0.95,
                 0.09256075757974097},
      std::tuple{0.7773502691896259,
                 0.37735026918962583,
                 0.1,
                 0.05,
                 0.007377231824157576},
      std::tuple{0.7773502691896259,
                 0.37735026918962583,
                 0.1,
                 0.35,
                 0.18788862212722682},
      std::tuple{0.7773502691896259,
                 0.37735026918962583,
                 0.1,
                 0.65,
                 0.6791707222540323},
      std::tuple{0.7773502691896259, 0.37735026918962583, 0.1, 0.95, 0},
      std::tuple{
          0.7773502691896259, 0.37735026918962583, 3.1622776601683795, 0.05, 0},
      std::tuple{
          0.7773502691896259, 0.37735026918962583, 3.1622776601683795, 0.35, 0},
      std::tuple{
          0.7773502691896259, 0.37735026918962583, 3.1622776601683795, 0.65, 0},
      std::tuple{0.7773502691896259,
                 0.37735026918962583,
                 3.1622776601683795,
                 0.95,
                 1.6545035829889259},
      std::tuple{0.7773502691896259,
                 0.5106836025229592,
                 0.0001,
                 0.05,
                 7.540003466141407e-06},
      std::tuple{0.7773502691896259,
                 0.5106836025229592,
                 0.0001,
                 0.35,
                 0.00022017376031957072},
      std::tuple{0.7773502691896259,
                 0.5106836025229592,
                 0.0001,
                 0.65,
                 0.0019291854113340373},
      std::tuple{0.7773502691896259,
                 0.5106836025229592,
                 0.0001,
                 0.95,
                 0.0030546091145512624},
      std::tuple{0.7773502691896259,
                 0.5106836025229592,
                 0.0031622776601683794,
                 0.05,
                 0.00023961116825033191},
      std::tuple{0.7773502691896259,
                 0.5106836025229592,
                 0.0031622776601683794,
                 0.35,
                 0.006919498291816282},
      std::tuple{0.7773502691896259,
                 0.5106836025229592,
                 0.0031622776601683794,
                 0.65,
                 0.053762638707409126},
      std::tuple{0.7773502691896259,
                 0.5106836025229592,
                 0.0031622776601683794,
                 0.95,
                 0.09256075757974097},
      std::tuple{0.7773502691896259,
                 0.5106836025229592,
                 0.1,
                 0.05,
                 0.007553835443188259},
      std::tuple{0.7773502691896259,
                 0.5106836025229592,
                 0.1,
                 0.35,
                 0.1870907093282112},
      std::tuple{0.7773502691896259,
                 0.5106836025229592,
                 0.1,
                 0.65,
                 0.5699414855437102},
      std::tuple{0.7773502691896259,
                 0.5106836025229592,
                 0.1,
                 0.95,
                 1.3456667575818533},
      std::tuple{
          0.7773502691896259, 0.5106836025229592, 3.1622776601683795, 0.05, 0},
      std::tuple{
          0.7773502691896259, 0.5106836025229592, 3.1622776601683795, 0.35, 0},
      std::tuple{
          0.7773502691896259, 0.5106836025229592, 3.1622776601683795, 0.65, 0},
      std::tuple{0.7773502691896259,
                 0.5106836025229592,
                 3.1622776601683795,
                 0.95,
                 1.6513526733164734},
      std::tuple{0.7773502691896259,
                 0.6440169358562926,
                 0.0001,
                 0.05,
                 7.692929307907293e-06},
      std::tuple{0.7773502691896259,
                 0.6440169358562926,
                 0.0001,
                 0.35,
                 0.00021746340965321482},
      std::tuple{0.7773502691896259,
                 0.6440169358562926,
                 0.0001,
                 0.65,
                 0.0010141027513213725},
      std::tuple{0.7773502691896259,
                 0.6440169358562926,
                 0.0001,
                 0.95,
                 0.0030546091145512624},
      std::tuple{0.7773502691896259,
                 0.6440169358562926,
                 0.0031622776601683794,
                 0.05,
                 0.00024358036052285205},
      std::tuple{0.7773502691896259,
                 0.6440169358562926,
                 0.0031622776601683794,
                 0.35,
                 0.006840181608917257},
      std::tuple{0.7773502691896259,
                 0.6440169358562926,
                 0.0031622776601683794,
                 0.65,
                 0.030816935402392814},
      std::tuple{0.7773502691896259,
                 0.6440169358562926,
                 0.0031622776601683794,
                 0.95,
                 0.09256075757974097},
      std::tuple{0.7773502691896259,
                 0.6440169358562926,
                 0.1,
                 0.05,
                 0.0076783467954213485},
      std::tuple{0.7773502691896259,
                 0.6440169358562926,
                 0.1,
                 0.35,
                 0.1886076171694935},
      std::tuple{0.7773502691896259,
                 0.6440169358562926,
                 0.1,
                 0.65,
                 0.5040849323535525},
      std::tuple{0.7773502691896259,
                 0.6440169358562926,
                 0.1,
                 0.95,
                 1.2811313744400128},
      std::tuple{
          0.7773502691896259, 0.6440169358562926, 3.1622776601683795, 0.05, 0},
      std::tuple{
          0.7773502691896259, 0.6440169358562926, 3.1622776601683795, 0.35, 0},
      std::tuple{
          0.7773502691896259, 0.6440169358562926, 3.1622776601683795, 0.65, 0},
      std::tuple{0.7773502691896259,
                 0.6440169358562926,
                 3.1622776601683795,
                 0.95,
                 1.6483214411861051},
      std::tuple{0.7773502691896259,
                 0.7773502691896259,
                 0.0001,
                 0.05,
                 7.806693241205771e-06},
      std::tuple{0.7773502691896259,
                 0.7773502691896259,
                 0.0001,
                 0.35,
                 0.0002173828693434697},
      std::tuple{0.7773502691896259,
                 0.7773502691896259,
                 0.0001,
                 0.65,
                 0.0007909289023483763},
      std::tuple{0.7773502691896259,
                 0.7773502691896259,
                 0.0001,
                 0.95,
                 0.0030546091145512624},
      std::tuple{0.7773502691896259,
                 0.7773502691896259,
                 0.0031622776601683794,
                 0.05,
                 0.00024648183988229506},
      std::tuple{0.7773502691896259,
                 0.7773502691896259,
                 0.0031622776601683794,
                 0.35,
                 0.006844294417901413},
      std::tuple{0.7773502691896259,
                 0.7773502691896259,
                 0.0031622776601683794,
                 0.65,
                 0.024405928427880138},
      std::tuple{0.7773502691896259,
                 0.7773502691896259,
                 0.0031622776601683794,
                 0.95,
                 0.09256075757974097},
      std::tuple{0.7773502691896259,
                 0.7773502691896259,
                 0.1,
                 0.05,
                 0.0077725405473192104},
      std::tuple{0.7773502691896259,
                 0.7773502691896259,
                 0.1,
                 0.35,
                 0.1908016693914296},
      std::tuple{0.7773502691896259,
                 0.7773502691896259,
                 0.1,
                 0.65,
                 0.4677255472039136},
      std::tuple{0.7773502691896259,
                 0.7773502691896259,
                 0.1,
                 0.95,
                 1.1733319381342693},
      std::tuple{
          0.7773502691896259, 0.7773502691896259, 3.1622776601683795, 0.05, 0},
      std::tuple{
          0.7773502691896259, 0.7773502691896259, 3.1622776601683795, 0.35, 0},
      std::tuple{
          0.7773502691896259, 0.7773502691896259, 3.1622776601683795, 0.65, 0},
      std::tuple{0.7773502691896259,
                 0.7773502691896259,
                 3.1622776601683795,
                 0.95,
                 1.6452148494846968});
  REQUIRE(BSMPT::kappa::kappaNuMuModel(cs2b, cs2s, al, vw) ==
          Approx(expected).epsilon(1e-3));
}