// SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#include <BSMPT/bounce_solution/action_calculation.h>
#include <BSMPT/utility/NumericalDerivatives.h>

namespace BSMPT
{

BounceActionInt::BounceActionInt()
{
}

BounceActionInt::BounceActionInt(
    std::vector<std::vector<double>> InitPath_In,
    std::vector<double> TrueVacuum_In,
    std::vector<double> FalseVacuum_In,
    std::function<double(std::vector<double>)> &V_In,
    std::function<std::vector<double>(std::vector<double>)> &dV_In,
    double T_In,
    int MaxPathIntegrations_In)
{
  // Initialization of the class when the derivative is provided
  this->dim    = InitPath_In.at(0).size();
  this->Vfalse = V_In(FalseVacuum_In);
  this->V  = [&](std::vector<double> vev) { return V_In(vev) - this->Vfalse; };
  this->dV = dV_In;
  this->Hessian = [=](auto const &arg)
  { return HessianNumerical(arg, V_In, this->eps); };
  this->TrueVacuum          = TrueVacuum_In;
  this->FalseVacuum         = FalseVacuum_In;
  this->InitPath            = InitPath_In;
  this->T                   = T_In;
  this->MaxPathIntegrations = MaxPathIntegrations_In;
  // Set Spline path
  SetPath(InitPath_In);
}

BounceActionInt::BounceActionInt(
    std::vector<std::vector<double>> InitPath_In,
    std::vector<double> TrueVacuum_In,
    std::vector<double> FalseVacuum_In,
    std::function<double(std::vector<double>)> &V_In,
    double T_In,
    int MaxPathIntegrations_In)
{
  // Initialization of the class when the derivative is not provided
  this->dim    = InitPath_In.at(0).size();
  this->Vfalse = V_In(FalseVacuum_In);
  this->V = [&](std::vector<double> vev) { return V_In(vev) - this->Vfalse; };
  // Use numerical derivative
  this->dV = [=](auto const &arg)
  { return NablaNumerical(arg, this->V, this->eps); };
  this->Hessian = [=](auto const &arg)
  { return HessianNumerical(arg, V_In, this->eps); };
  this->TrueVacuum          = TrueVacuum_In;
  this->FalseVacuum         = FalseVacuum_In;
  this->InitPath            = InitPath_In;
  this->T                   = T_In;
  this->MaxPathIntegrations = MaxPathIntegrations_In;
  // Set Spline path
  SetPath(InitPath_In);
}

void BounceActionInt::SetPath(std::vector<std::vector<double>> InitPath_In)
{
  // Method to be called when the path is changed manually
  this->Path = InitPath_In;
  // Find minimums near the initial and last position

  this->Spline = cvspline(this->Path);
  this->Path   = Spline.phipath;

  if (V(InitPath_In.front()) > V(InitPath_In.back()))
  {
    std::stringstream ss;
    ss << "-----------------------------------------------\n";
    ss << "Error with new path with length\t" << Spline.L << "\n";
    ss << "V(TrueVacuum) = " << V(InitPath_In.front()) << "\n";
    ss << "V(FalseVacuum) = " << V(InitPath_In.back()) << "\n";
    ss << "V(TrueVacuum) > V(FalseVacuum) < ----  This cannot be! Path might "
          "be backwards\n";
    ss << "-----------------------------------------------\n";
    BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
  }

  // Rasterize dVdl
  RasterizedVdl();
}

void BounceActionInt::RasterizedVdl(double l_start)
{
  // This method is used to calculate all dVdl before hand. Otherwise the code
  // will probably slow down
  std::vector<double> l_temp, dVdl_temp;

  for (int it = 0; it <= 1000; it++)
  {
    l_temp.push_back(l_start + it / 1000.0 * (Spline.L - l_start));
    dVdl_temp.push_back(Calc_dVdl(l_temp.back()));
  }
  // Set the not-a-knot boundary conditions
  RasterizeddVdl.set_boundary(
      tk::spline::not_a_knot, 0.0, tk::spline::not_a_knot, 0.0);
  RasterizeddVdl.set_points(l_temp, dVdl_temp);
}

void BounceActionInt::PrintVector(std::vector<double> vec)
{
  std::stringstream ss;
  ss << std::setprecision(15) << "[";
  for (double i : vec)
    ss << i << " ";
  ss << "]";
  BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
}

std::vector<double>
BounceActionInt::NormalForce(const double &l,
                             const double &dldrho,
                             const std::vector<double> &gradient)
{
  std::vector<double> r;                        // Result
  std::vector<double> d2phidl2 = Spline.d2l(l); // d2Phi/dl2
  std::vector<double> dphidl =
      Spline.dl(l); // dPhi/dl which norm is 1 (needed to calculate the
                    // perpendicular component)
  r = std::pow(dldrho, 2) * d2phidl2 -
      (gradient - (gradient * dphidl) * dphidl);

  return r;
}

std::vector<double>
BounceActionInt::NormalForceBernstein(const double &dldrho,
                                      const std::vector<double> &gradient,
                                      const std::vector<double> &dphidl,
                                      const std::vector<double> &d2phidl2)
{
  std::vector<double> r; // Result
  r = std::pow(dldrho, 2) * d2phidl2 -
      (gradient - (gradient * dphidl) * dphidl);

  return r;
}

double BounceActionInt::d2ldrho2(double l, double rho, double dldrho)
{
  if (dldrho == 0) // This is stupid, we should take a manual Runge-Kutta step
  {
    return (RasterizeddVdl(l));
  }
  return (RasterizeddVdl(l) - Alpha * dldrho / rho);
}

void BounceActionInt::AuxFunctionDev(const double &rho,
                                     const std::vector<double> &dvs,
                                     std::vector<double> &aks)
{
  // Array of the derivatives
  // Auxiliary function for the Runge Kutta 5 order method
  aks = {dvs[1], d2ldrho2(dvs[0], rho, dvs[1])};
  return;
}

void BounceActionInt::RK5_step(const std::vector<double> &y,
                               const std::vector<double> &dydx,
                               int n,
                               float rho,
                               float h,
                               std::vector<double> &yout,
                               std::vector<double> &yerr)
{
  int i;

  double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875, b21 = 0.2,
         b31 = 3.0 / 40.0, b32 = 9.0 / 40.0, b41 = 0.3, b42 = -0.9, b43 = 1.2,
         b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0, b54 = 35.0 / 27.0,
         b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0, b63 = 575.0 / 13824.0,
         b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0, c1 = 37.0 / 378.0,
         c3 = 250.0 / 621.0, c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0,
         dc5 = -277.00 / 14336.0;
  double dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0,
         dc4 = c4 - 13525.0 / 55296.0, dc6 = c6 - 0.25;

  std::vector<double> ak2(2);
  std::vector<double> ak3(2);
  std::vector<double> ak4(2);
  std::vector<double> ak5(2);
  std::vector<double> ak6(2);
  std::vector<double> ytemp(2);
  ytemp = std::vector<double>(2);

  // First step.
  for (i = 0; i < n; i++)
    ytemp[i] = y[i] + b21 * h * dydx[i];
  this->AuxFunctionDev(rho + a2 * h, ytemp, ak2);
  // Second step.
  for (i = 0; i < n; i++)
    ytemp[i] = y[i] + h * (b31 * dydx[i] + b32 * ak2[i]);
  this->AuxFunctionDev(rho + a3 * h, ytemp, ak3);
  // Third step.
  for (i = 0; i < n; i++)
    ytemp[i] = y[i] + h * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
  this->AuxFunctionDev(rho + a4 * h, ytemp, ak4);
  // Fourth step.
  for (i = 0; i < n; i++)
    ytemp[i] =
        y[i] + h * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);
  this->AuxFunctionDev(rho + a5 * h, ytemp, ak5);
  // Fifth step.
  for (i = 0; i < n; i++)
    ytemp[i] = y[i] + h * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] +
                           b64 * ak4[i] + b65 * ak5[i]);
  this->AuxFunctionDev(rho + a6 * h, ytemp, ak6);
  // Sixth step.
  for (i = 0; i < n; i++) // Accumulate increments with proper weights.
    yout[i] =
        y[i] + h * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
  for (i = 0; i < n; i++)
    yerr[i] = h * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] +
                   dc6 * ak6[i]);
  return;
}

double BounceActionInt::BesselI(double alpha, double x, int terms)
{
  // This implementation seems to converge quite quicly
  // https://en.wikipedia.org/wiki/Bessel_function#:~:text=Modified%20Bessel%20functions%3A%20I%CE%B1%2C%20K%CE%B1%5B,first%20and%20second%20kind%20and%20are%20defined%20as%5B19%5D
  double r0 = 1e100;
  double r  = 0;
  int m     = 0;
  while ((m < terms) && (abs((r - r0) / r) > 1e-15))
  {
    r0 = r; // Save step
    r += 1 / (tgamma(m + alpha + 1) * tgamma(m + 1)) *
         std::pow(x / 2.0, 2.0 * m + alpha); // One more term
    m++; // Update later to not mess up summation
  }
  return r;
}

double BounceActionInt::BesselJ(double x, int terms)
{
  // This implementation seems to converge quite quicly
  // https://en.wikipedia.org/wiki/Bessel_function#:~:text=Modified%20Bessel%20functions%3A%20I%CE%B1%2C%20K%CE%B1%5B,first%20and%20second%20kind%20and%20are%20defined%20as%5B19%5D
  double r0 = 1e100;
  double r  = 0;
  int m     = 0;
  while ((m < terms) && (abs((r - r0) / r) > 1e-15))
  {
    r0 = r; // Save step
    r += 1 / (tgamma(m + 2.) * tgamma(m + 1)) * std::pow(-1, m) *
         std::pow(x / 2.0, 2.0 * m + 1.); // One more term
    m++; // Update later to not mess up summation
  }
  return r;
}

std::vector<double> BounceActionInt::ExactSolutionFromMinimum(double l)
{
  double rho_down   = 1e-100;
  double rho_up     = 1;
  double rho_middle = 0;
  std::function<double(double)> LinearSolution, LinearSolutionDerivative;

  TrueVacuumHessian = Calc_d2Vdl2(Initial_lmin); // Sanity check
  if (Alpha == 2)
  {
    LinearSolution = [=](const double rho_in)
    {
      return l - (Initial_lmin + l0_minus_lmin *
                                     sinh(sqrt(TrueVacuumHessian) * rho_in) /
                                     (sqrt(TrueVacuumHessian) * rho_in));
    };
    LinearSolutionDerivative = [=](const double rho_in)
    {
      return ((l0_minus_lmin * cosh(rho_in * std::sqrt(TrueVacuumHessian))) /
              rho_in) -
             (l0_minus_lmin * sinh(rho_in * std::sqrt(TrueVacuumHessian))) /
                 (pow(rho_in, 2) * std::sqrt(TrueVacuumHessian));
    };
  }

  // T = 0
  if (Alpha == 3)
  {
    LinearSolution = [=](const double rho_in)
    {
      return l -
             (Initial_lmin + 2 * l0_minus_lmin *
                                 BesselI(1, sqrt(TrueVacuumHessian) * rho_in) /
                                 (sqrt(TrueVacuumHessian) * rho_in));
    };
    LinearSolutionDerivative = [=](const double rho_in)
    {
      return (
          (2 * l0_minus_lmin * BesselI(2, rho_in * sqrt(TrueVacuumHessian))) /
          rho_in);
    };
  }

  // Check if lower limit is viable
  assert(LinearSolution(rho_down) > 0);

  // Check if upper limit is viable
  while (LinearSolution(rho_up) > 0 and rho_up < 100)
  {
    // This is always possible
    rho_up += 1;
  }

  // Do binary search
  int cc = 0; // For safety
  while (rho_up - rho_down > 1e-10 and cc < 150)
  {
    rho_middle = (rho_up + rho_down) / 2;
    if (LinearSolution(rho_middle) > 0)
    {
      rho_down = rho_middle;
    }
    else
    {
      rho_up = rho_middle;
    }
    cc++;
  }
  return {rho_middle, l, LinearSolutionDerivative(rho_middle)};
}

std::vector<double> BounceActionInt::ExactSolutionLin(double l0,
                                                      double l,
                                                      double dVdl,
                                                      double d2Vdl2)
{
  // Numerical solution of equation of motion with V'(phi) = dVdl + (phi -
  // phi0) * d2Vdl2
  std::stringstream ss;
  double nu          = (Alpha - 1.0) / 2.0; // To keep things tidy
  double Abs_d2Vdl2  = abs(d2Vdl2);
  double Sign_d2Vdl2 = (d2Vdl2 > 0) - (d2Vdl2 < 0); // Sign(d2Vdl2)
  double rho_down    = 1e-100;
  double rho_middle  = 0;
  double rho_up      = 1;
  std::function<double(double)> LinearSolution, LinearSolutionDerivative;

  // Set maximum values that rho can have
  if (Alpha == 2 and d2Vdl2 > 0)
  {
    rho_up = 1;
  }
  else if (Alpha == 2 and d2Vdl2 < 0)
  {
    // Maximum value in oscillatory behaviour 4.493409457909 is the solution of
    // tan(x) = x
    rho_up = 4.493409457909 / std::sqrt(Abs_d2Vdl2);
  }
  else if (Alpha == 3 and d2Vdl2 > 0)
  {
    rho_up = 1;
  }
  else if (Alpha == 3 and d2Vdl2 < 0)
  {
    // Maximum value in oscillatory behaviour 5.13562230184068 is the solution
    // for BesselJ[0, x] - (2 BesselJ[1, x])/x - BesselJ[2, x]
    rho_up = 5.13562230184068 / std::sqrt(Abs_d2Vdl2);
  }
  // Write difference of l(rho) - l0
  // Goal is to find zero of this function
  if (Alpha == 2 and d2Vdl2 > 0)
  {
    LinearSolution = [=](const double rho_in)
    {
      return l - (l0 - dVdl / d2Vdl2 +
                  Sign_d2Vdl2 * dVdl * sinh(sqrt(Abs_d2Vdl2) * rho_in) /
                      (pow(Abs_d2Vdl2, 1.5) * rho_in));
    };
    LinearSolutionDerivative = [=](const double rho_in)
    {
      return (dVdl * cosh(sqrt(d2Vdl2) * rho_in)) / (d2Vdl2 * rho_in) -
             (dVdl * sinh(sqrt(d2Vdl2) * rho_in)) /
                 (pow(d2Vdl2, 1.5) * std::pow(rho_in, 2));
    };
  }
  else if (Alpha == 2 and d2Vdl2 < 0)
  {
    LinearSolution = [=](const double rho_in)
    {
      return l - (l0 - dVdl / d2Vdl2 +
                  Sign_d2Vdl2 * dVdl * sin(sqrt(Abs_d2Vdl2) * rho_in) /
                      (pow(Abs_d2Vdl2, 1.5) * rho_in));
    };

    LinearSolutionDerivative = [=](const double rho_in)
    {
      return -((dVdl * cos(sqrt(Abs_d2Vdl2) * rho_in)) /
               (Abs_d2Vdl2 * rho_in)) +
             (dVdl * sin(sqrt(Abs_d2Vdl2) * rho_in)) /
                 (pow(Abs_d2Vdl2, 1.5) * std::pow(rho_in, 2));
    };
  }
  else if (Alpha == 3 and d2Vdl2 > 0)
  {
    LinearSolution = [=](const double rho_in)
    {
      return l - (l0 - dVdl / d2Vdl2 +
                  2 * dVdl * BesselI(nu, std::sqrt(d2Vdl2) * rho_in) /
                      (pow(d2Vdl2, 1.5) * rho_in));
    };
    LinearSolutionDerivative = [=](const double rho_in)
    {
      return (-2 * dVdl * BesselI(1, std::sqrt(d2Vdl2) * rho_in)) /
                 (pow(d2Vdl2, 1.5) * std::pow(rho_in, 2)) +
             (dVdl * (BesselI(0, std::sqrt(d2Vdl2) * rho_in) +
                      BesselI(2, std::sqrt(d2Vdl2) * rho_in))) /
                 (d2Vdl2 * rho_in);
    };
  }
  else if (Alpha == 3 and d2Vdl2 < 0)
  {
    LinearSolution = [=](const double rho_in)
    {
      return l - (l0 - dVdl / d2Vdl2 -
                  2 * dVdl * BesselJ(sqrt(Abs_d2Vdl2) * rho_in) /
                      (pow(Abs_d2Vdl2, 1.5) * rho_in));
    };
    LinearSolutionDerivative = [=](const double rho_in)
    {
      // We used numerical derivative here because Bessel function J are not
      // completely implemented
      return (
          -(LinearSolution(rho_in + 0.001) - LinearSolution(rho_in - 0.001)) /
          (0.002));
    };
  }
  else
  {
    ss << "An error occured, Alpha must be 2 or 3 (Alpha = D - 1), instead "
          "Alpha is \t"
       << Alpha << "\n";
    BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
    ss.str(std::string());
    double small_step = 1e-5;
    return {
        small_step, l0 + small_step * small_step * dVdl / 2, small_step * dVdl};
  }

  // Check if lower limit is viable
  if (LinearSolution(rho_down) < 0)
  {
    ss << "Error in\t(LinearSolution(rho_down) < 0)\n";
    BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
    ss.str(std::string());
    double small_step = 1e-5;
    return {
        small_step, l0 + small_step * small_step * dVdl / 2, small_step * dVdl};
  }

  // Check if upper limit is viable
  // If d2Vdl2 then there is a possibility that there is no solution
  if (d2Vdl2 < 0)
  {
    if (LinearSolution(rho_up) > 0)
    {
      // TODO
      ss << "Error in\t(LinearSolution(rho_up) < 0). Call function again with "
            "smaller argument\n";
      ss << rho_up << "\t" << l << "\t" << dVdl << "\t" << d2Vdl2 << "\n";
      BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
      ss.str(std::string());
      if (abs((l - l0) / Spline.L) < 1e-10)
      {
        // Maximum numerical precision reached.
        StateOfBounceActionInt = ActionStatus::Integration1DFailed;
        // Abort calculation
        return {0, 0, 0};
      }
      return (ExactSolutionLin(l0, l0 + (l - l0) / 10., dVdl, d2Vdl2));
    }
  }
  else
  {
    // In this case we can solve the problem of the limit
    int counter = 0;
    while (LinearSolution(rho_up) > 0 and counter < 100)
    {
      counter++;
      // This is always possible
      rho_up += 1;
    }
  }
  // Do binary search
  int cc = 0; // For safety
  while (rho_up - rho_down > 1e-10 and cc < 150)
  {
    rho_middle = (rho_up + rho_down) / 2;
    if (LinearSolution(rho_middle) > 0)
    {
      rho_down = rho_middle;
    }
    else
    {
      rho_up = rho_middle;
    }
    cc++;
  }
  return {rho_middle, l, LinearSolutionDerivative(rho_middle)};
}

double BounceActionInt::Calc_dVdl(double l)
{
  // Function to calculate dV/dl
  return dV(Spline(l)) * Spline.dl(l);
}

double BounceActionInt::Calc_d2Vdl2(double l)
{
  // Function to calculate d2V/dl2
  return (dV(Spline(l)) * Spline.d2l(l)) +
         ((Hessian(Spline(l)) * Spline.dl(l)) * Spline.dl(l));
}

double BounceActionInt::LogisticFunction(const double &x)
{
  if (x >= 10) return 1;
  if (x <= -10) return 0;
  return 1 / (1 + exp(-x));
}

void BounceActionInt::CalculateExactSolutionThreshold(double MinError)
{
  double NumberOfSteps = 1000;
  double Error, l0, l, inital_exponent, final_exponent;

  std::vector<double> MinSol, LinSol;
  inital_exponent = log((Spline.L - Initial_lmin) / 100.) -
                    10; // Search from 4e-7% spline path
  final_exponent =
      log((Spline.L - Initial_lmin) / 100.); // Search until 10% spline path

  for (double exponent = inital_exponent; exponent <= final_exponent;
       exponent += (final_exponent - inital_exponent) / NumberOfSteps)
  {
    l0_minus_lmin = exp(exponent);
    l0            = Initial_lmin + l0_minus_lmin;
    l             = l0 + Spline.L * FractionOfThePathExact;
    MinSol        = ExactSolutionFromMinimum(l);
    LinSol        = ExactSolutionLin(l0, l, Calc_dVdl(l0), Calc_d2Vdl2(l0));
    Error         = 0.5 *
            abs((LinSol.at(0) - MinSol.at(0)) / (LinSol.at(0) + MinSol.at(0)));
    if (Error < MinError)
    {
      // Found a better ExactSolutionThreshold
      MinError               = Error;
      ExactSolutionThreshold = l0_minus_lmin;
    }
  }

  if (MinError > 1e-2) // Error not small enough
  {
    if (FractionOfThePathExact <= 1e-4) return;
    FractionOfThePathExact /= 10.;
    CalculateExactSolutionThreshold(MinError);
  }
}
std::vector<double> BounceActionInt::ExactSolution(double l0)
{
  // FractionOfThePathExact = 1e-5;
  // Solving FractionOfThePathExact of the path
  double l = l0 + Spline.L * FractionOfThePathExact;
  // Computing the second derivative of the potential beforehand
  double d2Vdl2 = Calc_d2Vdl2(l0);
  // dVdl
  double dVdl = Calc_dVdl(l0);

  std::stringstream ss;

  // In the case Backwards propagation failed we use only the Linear Solution
  if (not ExactSolutionThreshold.has_value())
    return ExactSolutionLin(l0, l, dVdl, d2Vdl2);

  if (dVdl <= 0)
  {
    // Assume negative grad is numerical error if close to the true minimum
    if (l0_minus_lmin / Spline.L < 1e-2) return ExactSolutionFromMinimum(l);
    // If we are not close to the minimum then probably l0 in not a solution as
    // it will roll backwards
    ss << " \n l = " << l0 << std::endl;
    ss << " dVdl = " << dVdl << std::endl;
    ss << " d2Vd2l = " << d2Vdl2 << std::endl;
    BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
    StateOfBounceActionInt = ActionStatus::UndershootOvershootNegativeGrad;
    return {};
  }

  std::vector<double> MinSol = ExactSolutionFromMinimum(l);
  std::vector<double> LinSol = ExactSolutionLin(l0, l, dVdl, d2Vdl2);

  // Calculate the LogisticExponent which indicates the contribution from both
  // branches
  double LogisticExponent = 100 *
                            (l0_minus_lmin - ExactSolutionThreshold.value()) /
                            (ExactSolutionThreshold.value() - Initial_lmin);

  // If the contribution is too close to 0 or 1 we assume only a single type of
  // solution. The protects the code againts NANs from the LinSol when
  // the gradient is negative due to numerical instabilities.
  if (LogisticFunction(-LogisticExponent) == 1) return MinSol;
  if (LogisticFunction(LogisticExponent) == 1) return LinSol;

  // Interpolate between both solution at ExactSolutionThreshold)
  std::vector<double> MinLinInterpolated = LinSol;

  MinLinInterpolated.at(0) = MinSol.at(0) * LogisticFunction(-LogisticExponent);
  MinLinInterpolated.at(0) += LinSol.at(0) * LogisticFunction(LogisticExponent);
  MinLinInterpolated.at(2) = MinSol.at(2) * LogisticFunction(-LogisticExponent);
  MinLinInterpolated.at(2) += LinSol.at(2) * LogisticFunction(LogisticExponent);

  return MinLinInterpolated;
}

void BounceActionInt::IntegrateBounce(double l0,
                                      UndershootOvershootStatus &conv,
                                      std::vector<double> &rho,
                                      std::vector<double> &l,
                                      std::vector<double> &dl_drho,
                                      std::vector<double> &d2l_drho2,
                                      int maxiter,
                                      double error,
                                      double eps_abs,
                                      double max_step)
{
  std::stringstream ss;
  // Integrate the bounce equation with x(rho) = "x0" until dx/drho < 0
  // (undershoot) or x(rho) > L (overshoot) or
  // converges
  double L = Spline.L;
  double step; // Integration step
  std::vector<double> ExactSol = ExactSolution(l0);

  if (StateOfBounceActionInt != ActionStatus::NotCalculated) return;

  rho       = {0, ExactSol.at(0)};  // Initial integration value for abcissas
  l         = {l0, ExactSol.at(1)}; // Guess for the bounce solution
  dl_drho   = {0, ExactSol.at(2)};  // Initial derivative = 0
  d2l_drho2 = {d2ldrho2(l0, 0, 0)};
  d2l_drho2.push_back(d2ldrho2(l.back(), rho.back(), dl_drho.back()));
  step = rho.back() / 100;

  std::vector<double> next_l_dldrho(
      2); // Save "l" and "dldrho" from the Runge-Kutta 5th order step.
  std::vector<double> err(
      2); // Save error from the "l" and "dldrho" Runge-Kutta 5th order
          // step, used to upgrade the step size.

  double delta0; // Wanted precision
  double delta1; // Step precision
  int it;        // Counter

  for (it = 0;
       (it < maxiter) &&
       (((dl_drho.back() > error) && ((l.back() - L) / L < error)) || it < 5);
       it++) // Take at least 3 steps (due to dldrho < 0 due to numerical
             // errors)
  {
    RK5_step({l.back(), dl_drho.back()},
             {dl_drho.back(), d2l_drho2.back()},
             2,
             rho.back(),
             step,
             next_l_dldrho,
             err);

    delta1 = std::max(abs(err[0]), abs(err[1]));
    delta0 = eps_abs * std::max(abs(l.back() + next_l_dldrho[0]),
                                abs(dl_drho.back() + next_l_dldrho[1]));

    // Update step list
    rho.push_back(rho.back() + step);
    l.push_back(next_l_dldrho[0]);
    dl_drho.push_back(next_l_dldrho[1]);
    d2l_drho2.push_back(d2ldrho2(l.back(), rho.back(), dl_drho.back()));

    if (max_step > 0)
    {
      step = std::min((step * std::pow(delta0 / delta1, 0.2)), max_step);
    }
    else
    {
      step = step * std::pow(delta0 / delta1, 0.2) / 2;
    }
  }
  if ((abs(dl_drho.back()) <= error) && (abs(l.back() - L) / L <= error))
  {
    ss << "Converged\t" << it << "\t" << l0 << "\t" << rho.back() << "\t"
       << l.back() << "\t" << dl_drho.back();
    conv = UndershootOvershootStatus::Converged;
  }
  else if (dl_drho.back() <= error)
  {
    ss << "Undershoot\t" << it << "\t" << l0 << "\t" << rho.back() << "\t"
       << l.back() << "\t" << dl_drho.back();
    rho.pop_back();
    l.pop_back();
    dl_drho.pop_back();
    d2l_drho2.pop_back();
    conv          = UndershootOvershootStatus::Undershoot;
    UndershotOnce = true;
  }
  else if ((l.back() - L) / L >= error)
  {
    ss << "Overshoot\t" << it << "\t" << l0 << "\t" << rho.back() << "\t"
       << l.back() << "\t" << dl_drho.back();
    conv         = UndershootOvershootStatus::Overshoot;
    OvershotOnce = true;
  }
  else
  {
    // This shouldnt happen
    conv = UndershootOvershootStatus::Overshoot;
  }
  BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
}

void BounceActionInt::BackwardsPropagation()
{
  std::stringstream ss;
  // Backwards propagation starting point finder
  double l0  = 0;
  double l00 = 1e100;
  for (int i = 0; i < 100; i++)
  {
    l00 = l0;
    l0 -= Calc_dVdl(l0) / Calc_d2Vdl2(l0);
    if (abs((l0 - l00) / Spline.L) < 1e-8 and Calc_d2Vdl2(l0) > 0 and
        l0 <= Spline.L / 100)
    {
      // Calculate the threshold between linear solution and solution from
      // minimum
      CalculateExactSolutionThreshold();
      Initial_lmin = l0;
      return;
    }
  }
  ss << "Backwards propagation did not work...\t" << l0
     << "\t using minus gradient method instead\n";
  BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
  ss.str(std::string());

  // Restart the loop without the second derivative.
  l0 = 0;
  for (int i = 0; i < 1000; i++)
  {
    l00 = l0;
    l0 -= Calc_dVdl(l0) / 100;
    if (abs((l0 - l00) / Spline.L) < 1e-8 and Calc_d2Vdl2(l0) > 0 and
        l0 <= Spline.L / 100)
    {
      // Calculate the threshold between linear solution and solution from
      // minimum
      CalculateExactSolutionThreshold();
      Initial_lmin = l0;
      return;
    }
  }
  ss << "Backwards propagation not converging\t" << l0
     << "\t using minus 0.1% Spline length as backwards propagation\n";
  BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
  ExactSolutionThreshold
      .reset(); // Use always the linear solution in exact solution
  Initial_lmin = -1 * Spline.L / 1000.0;
  return;
}

void BounceActionInt::Solve1DBounce(
    std::vector<double> &rho,
    std::vector<double> &l,
    std::vector<double> &dl_drho,
    std::vector<double> &d2l_drho2,
    double error,
    int maxiter) // Alpha = 2 at T > 0 and Alpha = 3 at T = 0
{

  std::stringstream ss;
  // Method to solve the bounce ODE
  double lmin, l0, lmax, L;
  UndershootOvershootStatus conv; // Converged?
  L = Spline.L;

  BackwardsPropagation();
  RasterizedVdl(Initial_lmin); // Update rasterized dVdl
  lmin = this->Initial_lmin;   // Lower interval
  lmax = L;                    // Uppter interval
  ss << "Backwards propagation : \t" << lmin << "\t" << Calc_dVdl(lmin) << "\t"
     << Calc_d2Vdl2(lmin) << "\n";
  if (ExactSolutionThreshold.has_value())
  {
    ss << "l_threshold =\t" << ExactSolutionThreshold.value() << "\n";
  }
  else
  {
    ss << "l_threshold was not been calculated\n";
  }

  if (V(Spline(lmin)) > V(Spline(lmax)))
  {
    ss << "Backwards propagation produced V(lmin) > V(L). Abort.";
    BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
    StateOfBounceActionInt = ActionStatus::BackwardsPropagationFailed;
    return;
  }
  TrueVacuumHessian = Calc_d2Vdl2(lmin);
  int j             = 0;
  int resolution    = 1000;
  while (j < resolution)
  {
    // This procedure is done such that the allowed search intervals goes
    // from true VEV up to the point where potential becomes lower than
    // Vfalse. This method does not focus on speed. This method allows for
    // only one solution to the bounce equation in this interval.
    j++;
    if (V(Spline(lmin + (L - lmin) * double(j) / double(resolution))) > 0)
    {
      lmax = lmin + (L - lmin) * double(j) / double(resolution);
      break;
    }
    if (Calc_dVdl((L - lmin) * double(j) / double(resolution)) < 0)
    {
      lmax = lmin + (L - lmin) * double(j) / double(resolution);
      break;
    }
  }

  ss << "Upper limit : l = \t" << lmax
     << "\t | V(l = 0) - V(TrueVacuum) = " << V(Spline(lmax)) - V(TrueVacuum)
     << "\n";
  BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
  ss.str(std::string());

  StateOf1DIntegration = Integration1DStatus::NotConverged;

  UndershotOnce = false;
  OvershotOnce  = false;

  l0 = (lmax + lmin) / 2.0; // Perform binary search

  // mu ~= log(l0 - lmin)
  double mu_min    = -200;
  double mu_max    = log(lmax - lmin);
  double mu_middle = (mu_min + mu_max) / 2;

  int mode = 0; // Binary search. 0 = linear, 1 = log
  for (int i = 0; i < maxiter; i++)
  {
    if (mode == 0)
    {
      l0            = (lmax + lmin) / 2.0; // Perform binary search
      l0_minus_lmin = l0 - Initial_lmin;
      IntegrateBounce(
          l0, conv, rho, l, dl_drho, d2l_drho2, 100000, error, error * 0.0015);
      if (StateOfBounceActionInt != ActionStatus::NotCalculated) return;
      if (rho.size() <= 7)
      {
        ss << "rho\t = ";
        BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
        ss.str(std::string());
        PrintVector(rho);
        ss << "l\t = ";
        BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
        ss.str(std::string());
        PrintVector(l);
        ss << "dl_drho\t = ";
        BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
        ss.str(std::string());
        PrintVector(dl_drho);
        ss << "d2l_drho2\t = ";
        BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
        ss.str(std::string());
        PrintVector(d2l_drho2);
        BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed,
                             "dVdl\t = " + std::to_string(Calc_dVdl(l0)));
        BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed,
                             "d2Vdl2\t = " + std::to_string(Calc_d2Vdl2(l0)));
        ss << "\n Overshoot/Undershoot method failed!\t";
        StateOf1DIntegration   = Integration1DStatus::NotConverged;
        StateOfBounceActionInt = ActionStatus::Integration1DFailed;
        break;
      }
      if (conv == UndershootOvershootStatus::Converged) // Solved!
      {
        ss << "\nFound Solution!\t" << l0 << " in\t" << i << "\titerations.";
        StateOf1DIntegration = Integration1DStatus::Converged;
        break;
      }
      if ((lmax - lmin) / L < error * 0.0000001)
      {
        // A solution was found
        if (OvershotOnce == true)
        {
          ss << "\nConverged due to proximity!\t" << l0 << "\t"
             << " in\t" << i << "\titerations.\t"
             << "\t" << (abs(dl_drho.back())) << "\t"
             << "\t" << (abs(l.back() - L) / L);
          StateOf1DIntegration = Integration1DStatus::Converged;
          break;
        }
        // Method never overshot. Switch to log scale
        mode   = 1;
        mu_max = log(lmax - lmin) + 2; // Give some margin for the binary search
      }
      if (conv == UndershootOvershootStatus::Undershoot) // Undershoot!
      {
        lmax = double(l0);
      }
      if (conv == UndershootOvershootStatus::Overshoot) // Overshoot!
      {
        lmin = double(l0);
      }
    }
    if (mode == 1)
    {
      mu_middle     = (mu_min + mu_max) / 2;
      l0_minus_lmin = exp(mu_middle);
      l0 = Initial_lmin + l0_minus_lmin; // Perform binary search in log space

      IntegrateBounce(
          l0, conv, rho, l, dl_drho, d2l_drho2, 100000, error, error * 0.0015);
      if (StateOfBounceActionInt != ActionStatus::NotCalculated) return;
      if (rho.size() <= 7)
      {
        ss << "rho\t = ";
        BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
        ss.str(std::string());
        PrintVector(rho);
        ss << "l\t = ";
        BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
        ss.str(std::string());
        PrintVector(l);
        ss << "dl_drho\t = ";
        BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
        ss.str(std::string());
        PrintVector(dl_drho);
        ss << "d2l_drho2\t = ";
        BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
        ss.str(std::string());
        PrintVector(d2l_drho2);
        ss << "\n Overshoot/Undershoot method failed!\t";

        StateOf1DIntegration   = Integration1DStatus::NotConverged;
        StateOfBounceActionInt = ActionStatus::Integration1DFailed;
        break;
      }
      if (conv == UndershootOvershootStatus::Converged) // Solved!
      {
        ss << "\nFound Solution!\t" << l0 << " in\t" << i << "\titerations.";
        StateOf1DIntegration = Integration1DStatus::Converged;
        break;
      }
      if (abs(mu_max - mu_min) < 0.0000001)
      {
        ss << "\nConverged due to proximity!\t" << l0 << "\t"
           << " in\t" << i << "\titerations.\t"
           << "\t" << (abs(dl_drho.back())) << "\t"
           << "\t" << (abs(l.back() - L) / L);
        StateOf1DIntegration = Integration1DStatus::Converged;
        break;
      }
      if (conv == UndershootOvershootStatus::Undershoot) // Undershoot!
      {
        mu_max = mu_middle;
      }
      if (conv == UndershootOvershootStatus::Overshoot) // Overshoot!
      {
        mu_min = mu_middle;
      }
    }
  }
  BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
}

double BounceActionInt::Bernstein(int n, int nu, double x)
{
  std::stringstream ss;
  // Implementation of Bernstein polynomials
  // https://en.wikipedia.org/wiki/Bernstein_polynomial#:~:text=In%20the%20mathematical%20field%20of,named%20after%20Sergei%20Natanovich%20Bernstein.
  if (nu < 0 || nu > n)
  {
    return 0;
  }
  if (x < 0 || x > 1.0001)
  {
    ss << "Incorrect argument in Bernstein polynomial ->\t" << x;
    BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
  }
  return (nChoosek(n, nu) * pow(x, nu) * pow(1 - x, n - nu));
}

double BounceActionInt::ReductorCalculator(const double &MaximumGradient)
{
  {
    return MaximumGradient / (Spline.L);
  };
}

bool BounceActionInt::PathDeformationCheck(std::vector<double> &l,
                                           tk::spline &rho_l_spl)
{
  std::stringstream ss;
  // First try at path deformation a
  // Calculate the 1D bounce and then deform the knots until normal force
  // vanishes Problems:
  // -> Cubic splines are too unstable so some smoothing algorithm has to
  // be used

  std::vector<double> gradient(dim, 0);
  std::vector<double> force(dim, 0);

  double delta = (l.back() - l.front()) /
                 (10 * NumberPathKnots);    // Difference between to knots
  double np                    = l.front(); // Start of the new path
  double MaximumForce          = 0;         // Save maximum force
  double PerpendicularGradient = 0; // Save maximum perpendicular gradient
  double MaximumGradient       = 0; // Save maximum gradient
  double MaximumRelativeError  = 0; // Save maximum force relative to gradient
  double Maximum_dldrho        = 0; // Save maximum dl/drho
  std::vector<double> phi;          // Temporary variables for calculating force

  // Creates new list of knots for the new Spline, that then are going to
  // be moved with a force
  for (np = l.front() + delta; np <= l.back() - delta / 10.0; np += delta)
  {
    phi      = Spline(np); // New knot on the splind
    gradient = dV(phi);    // Grandient of knot
    force    = NormalForce(np,
                        1 / rho_l_spl.deriv(1, np),
                        gradient); // Calculate force in the knot

    Maximum_dldrho = std::max(Maximum_dldrho, 1 / rho_l_spl.deriv(1, np));

    PerpendicularGradient = std::max(
        PerpendicularGradient, L2NormVector(NormalForce(np, 0, gradient)));

    MaximumGradient =
        std::max(MaximumGradient,
                 L2NormVector(gradient)); // Calculate maximum gradient
    MaximumForce         = std::max(MaximumForce,
                            L2NormVector(force)); // Calculates maximum force
    MaximumRelativeError = std::max(
        MaximumRelativeError,
        L2NormVector(force) /
            L2NormVector(gradient)); // Calculate maximum force relative
                                     // to gradient on that point
  }

  ss << "----------------\t Path deformation check\t----------------\n";

  double reductor = ReductorCalculator(MaximumGradient);

  ss << "\nMaximmum dl/drho\t" << Maximum_dldrho << "\n";
  ss << "Maximmum gradient\t" << MaximumGradient << "\n";
  ss << "Maximmum perpendicular gradient\t" << PerpendicularGradient << "\n";
  ss << "Maximmum force\t" << MaximumForce << "\n";
  ss << "Maximmum relative error\t" << MaximumRelativeError << "\n";
  ss << "Reductor is\t" << reductor << "\n";
  ss << "Spline length is\t" << Spline.L << "\n";
  // ss << "{" << Maximum_dldrho << ", " << PerpendicularGradient << ", "
  //    << MaximumGradient << ", " << Spline.L << ", " << reductor << "}\n";

  BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
  ss.str(std::string());

  // Check for convergence!

  if (MaximumRelativeError < 0.05)
  {
    // It converged!
    StateOfPathDeformation = PathDeformationStatus::Converged;
    ss << "Everything went well. Path deformation converged!\n";
    BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
    ss.str(std::string());
    return true;
  }
  return false;
}

void BounceActionInt::SinglePathDeformation(
    double &stepsize,
    double &reductor,
    std::vector<double> &l,
    tk::spline &rho_l_spl,
    std::vector<double> &l_fornextpath,
    std::vector<std::vector<double>> &best_path,
    std::vector<std::vector<double>> &next_path,
    double &MaximumGradient,
    double &MaximumForce,
    double &MaximumRelativeError,
    double &Maximum_dldrho,
    double &PerpendicularGradient,
    MatrixXd &inverseK,
    std::vector<std::vector<double>> &forces)
{
  double stepIncrease = 1.5;
  double stepDecrease = 5.;
  double reverseCheck = .15;
  double maxstep      = .1;
  double minstep      = 1e-4;

  std::vector<double> temp_phi(dim, 0);
  std::vector<double> temp_dphi(dim, 0);
  std::vector<double> temp_d2phi(dim, 0);
  std::vector<double> gradient(dim, 0);
  std::vector<double> force(dim, 0);
  std::vector<Eigen::VectorXd> BernsteinCoefficients(dim);

  double l_to_Bernstein, np;
  // Save initial and final parameterization
  double l0 = l.front();
  double lf = l.back();
  // Difference between to knots to calculate Bernstein kernel
  // Transposes the path of differences
  double delta = (lf - l0) / 300;
  std::vector<std::vector<double>> transposed_next_path =
      BSMPT::Transpose(next_path);
  std::vector<std::vector<double>> last_forces = forces;

  for (int d = 0; d < dim; d++)
  {
    // Calculation of the Bernstein Spline coefficient
    tk::spline next_path_spline(l_fornextpath, transposed_next_path[d]);
    std::vector<double> IntegralVector(BernsteinDegree, 0);
    for (int b_it = 0; b_it < BernsteinDegree; b_it++)
    {
      for (np = l0; np <= lf - delta / 10.0; np += delta)
      {
        IntegralVector[b_it] =
            IntegralVector[b_it] +
            Bernstein(BernsteinDegree, b_it, (np - l0) / (lf - l0)) *
                next_path_spline(np);
        IntegralVector[b_it] =
            IntegralVector[b_it] +
            4 *
                Bernstein(
                    BernsteinDegree, b_it, (np + delta / 2 - l0) / (lf - l0)) *
                next_path_spline(np + delta / 2);
        IntegralVector[b_it] =
            IntegralVector[b_it] +
            Bernstein(BernsteinDegree, b_it, (np + delta - l0) / (lf - l0)) *
                next_path_spline(np + delta);
      }
    }

    // Normalization in Simpson Integration
    IntegralVector = (delta / 6) * IntegralVector;

    // Converts into eigenvector
    Eigen::VectorXd EigenIntegralVector =
        Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(IntegralVector.data(),
                                                      IntegralVector.size());

    // Solves K s = b
    BernsteinCoefficients[d] = inverseK * EigenIntegralVector;

    for (int it_path = 0; it_path < NumberPathKnots; it_path++)
    {
      double temp = 0;
      for (int b_it = 0; b_it < BernsteinDegree; b_it++)
      {
        temp += Bernstein(BernsteinDegree,
                          b_it,
                          (l_fornextpath[it_path] - l0) / (lf - l0)) *
                BernsteinCoefficients[d](b_it);
      }
      next_path[it_path][d] = temp;
    }
  }

  double oldMaximumGradient       = MaximumGradient;
  double oldMaximumForce          = MaximumForce;
  double oldMaximumRelativeError  = MaximumRelativeError;
  double oldMaximum_dldrho        = Maximum_dldrho;
  double oldPerpendicularGradient = PerpendicularGradient;

  MaximumGradient       = 0;
  MaximumForce          = 0;
  MaximumRelativeError  = 0;
  Maximum_dldrho        = 0;
  PerpendicularGradient = 0;

  forces.clear();

  for (int it_path = 0; it_path < NumberPathKnots; it_path++)
  {
    l_to_Bernstein = (l_fornextpath[it_path] - l0) / (lf - l0);

    std::fill(temp_phi.begin(), temp_phi.end(), 0);
    std::fill(temp_dphi.begin(), temp_dphi.end(), 0);
    std::fill(temp_d2phi.begin(), temp_d2phi.end(), 0);

    for (int d = 0; d < dim; d++)
    {
      for (int b_it = 0; b_it < BernsteinDegree; b_it++)
      {
        temp_phi[d] += Bernstein(BernsteinDegree, b_it, l_to_Bernstein) *
                       BernsteinCoefficients[d](b_it);
        temp_dphi[d] +=
            BernsteinDegree *
            (Bernstein(BernsteinDegree - 1, b_it - 1, l_to_Bernstein) -
             Bernstein(BernsteinDegree - 1, b_it, l_to_Bernstein)) *
            BernsteinCoefficients[d](b_it) / (lf - l0);
        temp_d2phi[d] +=
            BernsteinDegree * (BernsteinDegree - 1) *
            (Bernstein(BernsteinDegree - 2, b_it - 2, l_to_Bernstein) -
             2 * Bernstein(BernsteinDegree - 2, b_it - 1, l_to_Bernstein) +
             Bernstein(BernsteinDegree - 2, b_it, l_to_Bernstein)) *
            BernsteinCoefficients[d](b_it) / std::pow(lf - l0, 2);
      }
    }

    gradient = dV(FalseVacuum + temp_phi);
    force = NormalForceBernstein(1 / rho_l_spl.deriv(1, l_fornextpath[it_path]),
                                 gradient,
                                 temp_dphi,
                                 temp_d2phi);

    forces.push_back(force);

    next_path[it_path] = next_path[it_path] + force / reductor;

    Maximum_dldrho = std::max(Maximum_dldrho,
                              1 / rho_l_spl.deriv(1, l_fornextpath[it_path]));

    PerpendicularGradient = std::max(
        PerpendicularGradient,
        L2NormVector(NormalForceBernstein(0, gradient, temp_dphi, temp_d2phi)));

    MaximumGradient =
        std::max(MaximumGradient,
                 L2NormVector(gradient)); // Calculate maximum gradient
    MaximumForce = std::max(MaximumForce,
                            L2NormVector(force)); // Calculates maximum force
    /*MaximumRelativeError = std::max(
        MaximumRelativeError,
        L2NormVector(force) /
            L2NormVector(gradient)); // Calculate maximum force relative
    // to gradient on that point*/
    MaximumRelativeError = MaximumForce / MaximumGradient;
  }

  // Convergence check

  std::vector<double> reverser(forces.size(), 0);

  BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed,
                       "Path deformation error (before "
                       "integrating): " +
                           std::to_string(MaximumRelativeError));

  if (oldMaximumRelativeError > MaximumRelativeError)
  {
    best_path = next_path;
    BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed,
                         "Next best path found with error: " +
                             std::to_string(MaximumRelativeError));
  }
  else
  {
    MaximumGradient       = oldMaximumGradient;
    MaximumForce          = oldMaximumForce;
    MaximumRelativeError  = oldMaximumRelativeError;
    Maximum_dldrho        = oldMaximum_dldrho;
    PerpendicularGradient = oldPerpendicularGradient;
  }

  //  Update stepsize
  //  Calculates fraction of forces that switches direction betweeen iteration
  //  It this fraction is too big, the stepsize gets reduced , otherwise it
  //  gets increased
  if (last_forces.size() > 0)
  {
    std::transform(forces.begin(),
                   forces.end(),
                   last_forces.begin(),
                   reverser.begin(),
                   [](std::vector<double> &a, std::vector<double> &b)
                   { return a * b < 0; });
    if (std::accumulate(reverser.begin(), reverser.end(), 0.0) >
        forces.size() * reverseCheck)
    {
      next_path = best_path;
      stepsize /= stepDecrease;
    }
    else
    {
      stepsize *= stepIncrease;
    }
    stepsize = std::min(stepsize, maxstep);
    stepsize = std::max(stepsize, minstep);
  }

  last_forces = forces;
}

void BounceActionInt::PathDeformation(std::vector<double> &l,
                                      tk::spline &rho_l_spl)
{
  // First try at path deformation a
  // Calculate the 1D bounce and then deform the knots until normal force
  // vanishes Problems:
  // -> Cubic splines are too unstable so some smoothing algorithm has to
  // be used
  // -> Smoothing algorithm always spoils, in some way, the solution

  double reductor; // Computed factor to reduce the force
  int NoBestPathCounter = 0;

  double stepsize          = 2e-5;
  double SatisfactoryError = 0.05; // Maximum relative error which we consider a
                                   // success (we must still integrate again)

  double delta =
      (l.back() - l.front()) / NumberPathKnots; // Difference between to knots

  double oldMaximumRelativeError;
  double MaximumForce          = 0; // Save maximum force
  double MaximumGradient       = 0; // Save maximum gradient
  double PerpendicularGradient = 0; // Save maximum perpendicular gradient
  double Maximum_dldrho        = 0; // Save maximum dl/drho
  double MaximumRelativeError =
      1e100; // Save maximum force relative to gradient
  std::vector<double> phi, l_fornextpath;
  std::vector<std::vector<double>> next_path, old_path, transposed_next_path,
      best_path, last_forces, forces; // Next iteration path
  std::vector<double> gradient(dim, 0);
  std::vector<double> force(dim, 0);

  // Converting into Berenstein Basis!
  // Initialize K matrix
  // K_ij = int_0^1 Bi(x)Bj(x) dx
  MatrixXd K = MatrixXd::Zero(BernsteinDegree, BernsteinDegree);
  for (int i = 0; i < BernsteinDegree; i++)
  {
    for (int j = 0; j < BernsteinDegree; j++)
    {
      K(i, j) =
          (l.back() - l.front()) * double(nChoosek(BernsteinDegree, i)) *
          nChoosek(BernsteinDegree, j) /
          (nChoosek(2 * BernsteinDegree, i + j) * (2 * BernsteinDegree + 1));
    }
  }
  MatrixXd inverseK = K.inverse();

  // Creates new list of knots for the new Spline, that then are going to
  // be moved with a force
  for (double np = l.front(); np <= l.back() - delta / 10.0; np += delta)
  {
    phi = Spline(np);                       // New knot on the splind
    best_path.push_back(phi - FalseVacuum); // Add point to list
    l_fornextpath.push_back(np);            // Save the parameter for each point
    gradient = dV(phi);                     // Grandient of knot
    force    = NormalForce(np,
                        1 / rho_l_spl.deriv(1, np),
                        gradient); // Calculate force in the knot

    Maximum_dldrho = std::max(Maximum_dldrho, 1 / rho_l_spl.deriv(1, np));

    PerpendicularGradient = std::max(
        PerpendicularGradient, L2NormVector(NormalForce(np, 0, gradient)));

    MaximumGradient =
        std::max(MaximumGradient,
                 L2NormVector(gradient)); // Calculate maximum gradient
    MaximumForce         = std::max(MaximumForce,
                            L2NormVector(force)); // Calculates maximum force
    MaximumRelativeError = std::max(
        MaximumRelativeError,
        L2NormVector(force) /
            L2NormVector(gradient)); // Calculate maximum force relative
                                     // to gradient on that point
  }

  phi = FalseVacuum;                      // Last point
  best_path.push_back(phi - FalseVacuum); // Add last point to list
  l_fornextpath.push_back(Spline.L);      // Save the parameter for each point
  old_path  = best_path; // Save path in the case "oh no" happens
  next_path = best_path; // Starting path if the last iteration
  BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed,
                       "----------------\tPath deformation\t----------------");
  for (int it_maxpath = 0; it_maxpath < MaxSinglePathDeformations; it_maxpath++)
  {
    NoBestPathCounter++;
    reductor = ReductorCalculator(MaximumGradient) / stepsize;

    oldMaximumRelativeError = MaximumRelativeError;

    SinglePathDeformation(stepsize,
                          reductor,
                          l,
                          rho_l_spl,
                          l_fornextpath,
                          best_path,
                          next_path,
                          MaximumGradient,
                          MaximumForce,
                          MaximumRelativeError,
                          Maximum_dldrho,
                          PerpendicularGradient,
                          inverseK,
                          forces);

    if (MaximumRelativeError != oldMaximumRelativeError) NoBestPathCounter = 0;

    if (MaximumRelativeError < SatisfactoryError)
    {
      BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed,
                           "Enough convergence!\t" +
                               std::to_string(MaximumRelativeError));

      PathDeformationConvergedWithout1D = true;
      break;
    }
    if (MaximumRelativeError > 5) break; // Things went very wrong
    if (NoBestPathCounter > 20)
      break; // Algorithm does not seem to find new best paths
  }

  // Reshift all point to their correct locations
  if (best_path.size() == 0)
  {
    StateOfBounceActionInt = ActionStatus::PathDeformationCrashed;
    BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed,
                         "Path deformation exploded");
    return;
  }
  else if (MaximumRelativeError > SatisfactoryError)
  {
    BSMPT::Logger::Write(
        BSMPT::LoggingLevel::BounceDetailed,
        "Maximum iterations reached without error increasing. Integrate "
        "again!\n");
  }

  for (int it_path = 0; it_path <= NumberPathKnots; it_path++)
  {
    best_path[it_path] = FalseVacuum + best_path[it_path];
  }
  // Change the class path into the new path
  SetPath(best_path);

  return;
}

unsigned BounceActionInt::nChoosek(unsigned n, unsigned k)
{
  // Auxiliary function to calculate n choose k combination
  if (k > n) return 0;
  if (k * 2 > n) k = n - k;
  if (k == 0) return 1;

  int result = n;
  for (std::size_t i = 2; i <= k; ++i)
  {
    result *= (n - i + 1);
    result /= i;
  }
  return result;
}

double
BounceActionInt::CalculateKineticTermAction(const std::vector<double> &rho,
                                            const tk::spline &dl_drho_spl)
{
  double integral  = 0;
  double int_delta = rho[rho.size() - 1] / 2000;
  if (Alpha == 2)
  {
    for (double r = 0.0; r <= rho[rho.size() - 1]; r += int_delta)
    {
      // Simpson Integration (1 + 4 + 1)/ 6 * step
      integral += r * r * (0.5 * std::pow(dl_drho_spl(r), 2));
      integral +=
          4 * r * r * (0.5 * std::pow(dl_drho_spl(r + int_delta / 2.0), 2));
      integral += r * r * (0.5 * std::pow(dl_drho_spl(r + int_delta), 2));
    }
    integral = integral * 4 * M_PI * int_delta /
               6.0; // Angular integration and Simpson step
    return integral;
  }
  else if (Alpha == 3)
  {
    for (double r = 0.0; r <= rho[rho.size() - 1]; r += int_delta)
    {
      // Simpson Integration (1 + 4 + 1)/ 6 * step
      integral += r * r * r * (0.5 * std::pow(dl_drho_spl(r), 2));
      integral +=
          4 * r * r * r * (0.5 * std::pow(dl_drho_spl(r + int_delta / 2.0), 2));
      integral += r * r * r * (0.5 * std::pow(dl_drho_spl(r + int_delta), 2));
    }
    integral = integral * 2 * M_PI * M_PI * int_delta /
               6.0; // Angular integration and Simpson step
    return integral;
  }
  return -1;
}

double
BounceActionInt::CalculatePotentialTermAction(const std::vector<double> &rho,
                                              const tk::spline &l_rho_spl)
{
  double integral  = 0;
  double int_delta = rho[rho.size() - 1] / 2000;
  if (Alpha == 2)
  {
    for (double r = 0.0; r <= rho[rho.size() - 1]; r += int_delta)
    {
      // Simpson Integration (1 + 4 + 1)/ 6 * step
      integral += r * r * (V(Spline(l_rho_spl(r))));
      integral += 4 * r * r * (V(Spline(l_rho_spl(r + int_delta / 2.0))));
      integral += r * r * (V(Spline(l_rho_spl(r + int_delta))));
    }
    integral = integral * 4 * M_PI * int_delta /
               6.0; // Angular integration and Simpson step
    return integral;
  }
  else if (Alpha == 3)
  {
    for (double r = 0.0; r <= rho[rho.size() - 1]; r += int_delta)
    {
      // Simpson Integration (1 + 4 + 1)/ 6 * step
      integral += r * r * r * (V(Spline(l_rho_spl(r))));
      integral += 4 * r * r * r * (V(Spline(l_rho_spl(r + int_delta / 2.0))));
      integral += r * r * r * (V(Spline(l_rho_spl(r + int_delta))));
    }
    integral = integral * 2 * M_PI * M_PI * int_delta /
               6.0; // Angular integration and Simpson step
    return integral;
  }
  return -1;
}

void BounceActionInt::CalculateAction(
    double error) // Alpha = 2 at T > 0 and Alpha = 3 at T = 0
{
  std::stringstream ss;
  if (Calc_d2Vdl2(Spline.L) < 0)
  {
    StateOfBounceActionInt = ActionStatus::FalseVacuumNotMinimum;
    ss << "False vacuum is not a minimum!\n";
    BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
    ss.str(std::string());
    return;
  }

  std::vector<double> rho;       // List of rho
  std::vector<double> l;         // List of l(rho)
  std::vector<double> dl_drho;   // List of d^x/drho(rho)
  std::vector<double> d2l_drho2; // List of d2^x/drho2(rho)

  tk::spline rho_l_spl;     // Spline to find rho as a function of l
  tk::spline l_rho_spl;     // Spline to find l as a function of rho
  tk::spline dl_drho_l_spl; // Spline to find dldrho as a function of l
  tk::spline dl_drho_spl;   // Spline to find dldrho as a function of rho
  tk::spline d2l_drho2_spl; // Spline to find d2ldrho2 as a function of rho

  std::vector<double> force, originalforce;

  StateOf1DIntegration =
      Integration1DStatus ::NotConverged; // Know if a solution was found
  StateOfPathDeformation =
      PathDeformationStatus::NotConverged; // To record if path deformation
                                           // converged or not

  if (dim > 1) // If dim = 0 then we only need the solution for the 1D
               // bounce equation
  {
    BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed,
                         "--------------------------------\t1\t----------------"
                         "-----------------");
    // Solves 1D bounce equation
    Solve1DBounce(rho, l, dl_drho, d2l_drho2, error);

    if (StateOfBounceActionInt != ActionStatus::NotCalculated) return;

    if (rho.size() < 4)
    {
      if (rho.size() < 4)
        StateOfBounceActionInt = ActionStatus::NotEnoughPointsForSpline;
      // actions is less than 1 in case of an error. Abort
      // calculation
      Spline.print_path();
      return;
    }
    // Checks if convergence was met
    rho_l_spl.set_points(l, rho);
    l_rho_spl.set_points(rho, l);
    dl_drho_spl.set_points(rho, dl_drho);
    PathDeformationCheck(l, rho_l_spl);
    // Save solution
    rho_sol    = rho;
    l_sol      = l;
    dldrho_sol = dl_drho;
    // Checks if convergence was met

    for (int i = 2;
         i <= MaxPathIntegrations and
         StateOfPathDeformation == PathDeformationStatus::NotConverged;
         i++)
    {
      BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
      ss.str(std::string());
      if (UndershotOnce == false and Action >= -1)
      {
        BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed,
                             "Method never undershot. Terrible news!");
        StateOfBounceActionInt = ActionStatus::NeverUndershootOvershoot;
      }
      if (OvershotOnce == false and Action >= -1)
      {
        BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed,
                             "Method never overshot. Terrible news!");
        StateOfBounceActionInt = ActionStatus::NeverUndershootOvershoot;
      }

      if (StateOfBounceActionInt != ActionStatus::NotCalculated)
      {
        // actions is less than 1 in case of an error. Abort
        // calculation
        Spline.print_path();
        return;
      }

      if (rho.size() < 4)
      {
        // Not enough point to populate the Spline
        StateOfBounceActionInt = ActionStatus::NotEnoughPointsForSpline;
        Spline.print_path();
        return;
      }
      rho_l_spl.set_points(l, rho);
      l_rho_spl.set_points(rho, l);
      dl_drho_spl.set_points(rho, dl_drho);
      d2l_drho2_spl.set_points(rho, d2l_drho2);
      dl_drho_l_spl.set_points(l, dl_drho);

      BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed,
                           "--------------------------------\t" +
                               std::to_string(i) +
                               "\t---------------------------------\n");
      // Deform path
      PathDeformation(l, rho_l_spl);

      // Solves 1D bounce equation
      Solve1DBounce(rho,
                    l,
                    dl_drho,
                    d2l_drho2,
                    error); // Solves bounce equation

      if (UndershotOnce == false and Action >= -1)
      {
        BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed,
                             "Method never undershot. Terrible news!");
        StateOfBounceActionInt = ActionStatus::NeverUndershootOvershoot;
      }
      if (OvershotOnce == false and Action >= -1)
      {
        BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed,
                             "Method never overshot. Terrible news!");
        StateOfBounceActionInt = ActionStatus::NeverUndershootOvershoot;
      }
      if (StateOfBounceActionInt != ActionStatus::NotCalculated)
      {
        // actions is less than 1 in case of an error. Abort
        // calculation
        Spline.print_path();
        return;
      }

      // Path deformation converged
      if (PathDeformationConvergedWithout1D)
      {
        // Save solution
        StateOfPathDeformation = PathDeformationStatus::Converged;
        rho_sol                = rho;
        l_sol                  = l;
        dldrho_sol             = dl_drho;
        break;
      }

      // Checks if convergence was met
      if (PathDeformationCheck(l, rho_l_spl))
      {
        // Save solution
        rho_sol    = rho;
        l_sol      = l;
        dldrho_sol = dl_drho;
        break;
      }

      if (StateOfPathDeformation == PathDeformationStatus::Converged)
      {
        // PathDeformation Converged!
        break;
      }

      if (StateOf1DIntegration == Integration1DStatus ::NotConverged)
      {
        StateOfBounceActionInt =
            ActionStatus::Integration1DFailed; // Undershoot/overshoot failed.
        BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
      }

      if (StateOfBounceActionInt != ActionStatus::NotCalculated)
      {
        // actions is less than 1 in case of an error. Abort
        // calculation
        Spline.print_path();
        return;
      }
    }
    if (StateOfPathDeformation == PathDeformationStatus::NotConverged)
    {
      StateOfBounceActionInt =
          ActionStatus::PathDeformationNotConverged; // Path deformation did not
                                                     // converged in time.
      BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
      return;
    }
  }
  else
  {
    // Path deformation in not necessary
    Solve1DBounce(rho,
                  l,
                  dl_drho,
                  d2l_drho2,
                  error); // Solves bounce equation once
  }

  if (StateOf1DIntegration == Integration1DStatus::NotConverged)
  {
    StateOfBounceActionInt =
        ActionStatus::Integration1DFailed; // No solution was found.
    BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
    return;
  }

  rho_l_spl.set_points(l, rho);
  l_rho_spl.set_points(rho, l);
  dl_drho_spl.set_points(rho, dl_drho);

  double KineticPart = CalculateKineticTermAction(rho, dl_drho_spl);
  double KineticAction =
      2 * KineticPart /
      (1 +
       Alpha); // Calculate the Action using only the kinetical contributions
  double PotentialPart = CalculatePotentialTermAction(rho, l_rho_spl);
  double PotentialAction =
      2 * PotentialPart /
      (1 -
       Alpha); // Calculate the Action using only the potential contributions

  Action = KineticPart + PotentialPart;

  // Print warning the actions dffer by 10%
  if (abs(Action / KineticAction - 1) > 0.1)
  {
    ss << "Warning! Mismatch between Action and Action calculated using "
          "only kinetic term : "
       << Action << " and " << KineticAction
       << ". Relative error = " << abs(Action / KineticAction - 1) << "\n";
  }
  // Print warning the actions dffer by 10%
  if (abs(Action / PotentialAction - 1) > 0.1)
  {
    ss << "Warning! Mismatch between Action and Action calculated using "
          "only potential term : "
       << Action << "  and " << PotentialAction
       << ". Relative error = " << abs(Action / PotentialAction - 1) << "\n";
  }

  BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
  BSMPT::Logger::Write(
      BSMPT::LoggingLevel::BounceDetailed,
      "Distance from true vacuum to spline(l0) = " +
          std::to_string(L2NormVector(TrueVacuum - Spline(l.at(0)))));
  BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed,
                       "\nAction =\t" + std::to_string(Action) + "\tat T =\t" +
                           std::to_string(T) + "\n");
  return;
}
} // namespace BSMPT