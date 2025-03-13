// SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later

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

using namespace BSMPT;

int main()
{
  // Set the logger
  SetLogger({"--logginglevel::complete=true"});

  // Define the potential
  // BounceActionInt only takes std::vector functions as arguments.
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

  // Define the potential's gradient (optional)
  std::function<std::vector<double>(std::vector<double>)> dV =
      [&](std::vector<double> l0)
  {
    std::vector<double> result = {
        2 * l0[0] * (5 * pow(-1 + l0[0], 2) + pow(-1 + l0[1], 2)) +
            10 * (-1 + l0[0]) * (pow(l0[0], 2) + 5 * pow(l0[1], 2)),
        10 * (5 * pow(-1 + l0[0], 2) + pow(-1 + l0[1], 2)) * l0[1] +
            2 * (-1 + l0[1]) * (pow(l0[0], 2) + 5 * pow(l0[1], 2)) +
            80 * (-1. * pow(l0[1], 2) + 1. * pow(l0[1], 3))};
    return result;
  };

  // Define the false and true vacuum
  std::vector<double> FalseVacuum = {0, 0};
  std::vector<double> TrueVacuum  = {1, 1};

  // Your best guess for the path
  std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

  // Calculate the action if NO gradient was provided
  BounceActionInt bc(path, TrueVacuum, FalseVacuum, V, 0, 6);
  bc.CalculateAction();

  // Calculate the action if gradient was provided
  BounceActionInt bcdV(path, TrueVacuum, FalseVacuum, V, dV, 0, 6);
  bcdV.CalculateAction();

  std::cout << "Action calculated using numerical derivatives is " << bc.Action
            << "\n";
  std::cout << "Action calculated using analytical derivatives is "
            << bcdV.Action << "\n";

  return 0;
}