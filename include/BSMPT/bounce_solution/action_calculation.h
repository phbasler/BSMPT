// Copyright (C) 2024 Lisa Biermann, Margarete Mühlleitner, Rui Santos, João
// Viana
//
// SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

/**
 * @file
 */

// Remove warning from <experimental/filesystem>
#define _LIBCPP_NO_EXPERIMENTAL_DEPRECATION_WARNING_FILESYSTEM

// M_PI in Windows
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#include <BSMPT/utility/const_velocity_spline.h>
#include <BSMPT/utility/utility.h>
#include <Eigen/Dense>
#include <numeric>
#include <sys/stat.h>
#include <sys/types.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * @brief This classes calculates the Bounce action of the potential with a set
 * temperature
 */

namespace BSMPT
{

class BounceActionInt
{
private:
  // These two variables check if the method undershot and overshot at least
  // once. If only one type of convergence was achieved it is an indication that
  // no solution has been found.
  /**
   * @brief Records if overshoot/undershoot method undershots at least once
   */
  bool undershot_once = false;
  /**
   * @brief Records if overshoot/undershoot method overshots at least once
   */
  bool overshot_once = false;

  /**
   * @brief Value of d2Vdl^2 near the true vacuum
   *
   */
  double TrueVacuumHessian;

  /**
   * @brief Store the value of backwards propagation
   *
   */
  double Initial_lmin;

  /**
   * @brief Potential at the false vacuum in the unshifted potential (in the
   * shifted Vfalse = 0)
   *
   */
  double Vfalse;

  /**
   * @brief l0 - Initial_lmin for solutions starting very near the true vacuum
   *
   */
  double l0_Initial_lmin;

  /**
   * @brief True if path deformation reached the desired results without solving
   * the 1D equation one more time.
   *
   */
  bool PathDeformationConvergedWithout1D = false;

public:
  /**
   * @brief Dimension of the VEV space
   *
   */
  int dim = -1;

  /**
   * @brief Action calculated by the class
   *
   * * = -1 if not yet calculated.
   * * = -2 if 1D integration did not work.
   * * = -3 if the path deformation did not converged in time.
   * * = -4 if the path deformation crashed.
   * * = -5 if the false vacuum in not a minimum.
   * * = -6 path deformation created a path with V(True) > V(False) (even with
   * backwards propagation).
   * * = -7 1D integration never overshot/undershot
   * * = -8 overshooting/undershooting found a negative gradient point
   * * = -9 if the solution does not have enough points to populate the \f$
   * l(\rho) \f$ spline
   */
  enum action_status
  {
    not_calculated                     = -1,
    integration_1d_failed              = -2,
    path_deformation_not_converged     = -3,
    path_deformation_crashed           = -4,
    false_vacuum_not_a_minimum         = -5,
    backwards_propagation_failed       = -6,
    never_undershoot_overshoot         = -7,
    undershoot_overshoot_negative_grad = -8,
    not_enough_points_for_spline       = -9
  };

  /**
   * @brief either returns a action_status or the value of the action
   *
   */
  double action = not_calculated;

  /**
   * @brief Factor produced by the spherical symmetry of the potential.
   * * = 2 if \f$ T > 0\f$ (\f$O(3)\f$ symmetry).
   * * = 3 if \f$ T = 0\f$ (\f$O(4)\f$ symmetry).
   *
   * Default value is \f$ 2 \f$ since most of our calculation are done at finite
   * temperature.
   */
  double alpha = 2; // alpha = 2 if T > 0 | alpha = 3 if T = 0

  /**
   * @brief Temperature of the potential. Irrelevant but helpful
   *
   */
  double T = -1;

  /**
   * @brief Number of integration of the bounce
   *
   */
  int MaxPathIntegrations;

  /**
   * @brief Number of path deformations before integrating again
   *
   */
  int MaxSinglePathDeformations = 200;

  /**
   * @brief list of \f$ \rho \f$ of the solution
   */
  std::vector<double> rho_sol;

  /**
   * @brief list of \f$ l(\rho) \f$ of the solution
   *
   */
  std::vector<double> l_sol;

  /**
   * @brief list of \f$ \frac{l}{\rho} \f$ of the solution
   *
   */
  std::vector<double> dldrho_sol;

  double eps = 0.01; // Step for the numerical derivative
  /**
   * @brief Records if 1D profile converged or not
   *
   */
  int convS3 = -1; // Did 1D profile converged or not?
  /**
   * @brief Records if path deformation converged or not
   *
   */
  int convPathDeformation = -1;
  /**
   * @brief  Number of basis function that are used + 1
   *
   */
  int BernsteinDegree = 10;

  /**
   * @brief // Number of knots in the new path
   *
   */
  double NumberPathKnots = 50;
  /**
   * @brief True vacuum candidate
   *
   */
  std::vector<double> TrueVacuum;
  /**
   * @brief False vacuum. Should be the same as the last path knot
   *
   */
  std::vector<double> FalseVacuum;
  /**
   * @brief Potential of the class
   *
   */

  std::function<double(std::vector<double>)> V; // Potential
  /**
   * @brief Potential gradient of the class, can be either numerical or
   * analytical
   *
   */
  std::function<std::vector<double>(std::vector<double>)>
      dV; // Potential gradient
  /**
   * @brief Potential hessian of the class, completly numerical
   *
   * TODO: Calculate hessian from analytical gradient
   */
  std::function<std::vector<std::vector<double>>(std::vector<double>)> Hessian;
  /**
   * @brief First path given to class
   *
   */
  std::vector<std::vector<double>> init_path; // Initial path given to class
  /**
   * @brief Current class path, can be changed by @ref PathDeformation
   *
   */
  std::vector<std::vector<double>> path; // Initial path given to class
  /**
   * @brief We describe the tunneling path using a cubid spline.
   *  The parameterization is the length along the spline.
   *
   */
  cvspline::cvspline spline; // Constant velocity spline object
  /**
   * @brief Spline used to save \f$ \frac{dV}{dl} \f$
   *
   */
  tk::spline rasterizedVdl; // RasterizeddVdl

  /**
   * @brief Construct a new Bounce Action Int object
   *
   * @param init_path is the initial path guess
   * @param TrueVacuumIn is the true vacuum candidate of the potential
   * @param FalseVacuumIn is the false vacuum
   * @param V is the class potential
   * @param dV is the class gradient
   */
  BounceActionInt(std::vector<std::vector<double>> init_path_In,
                  std::vector<double> TrueVacuum_In,
                  std::vector<double> FalseVacuum_In,
                  std::function<double(std::vector<double>)> V_In,
                  std::function<std::vector<double>(std::vector<double>)> dV_In,
                  double T_In,
                  int MaxPathIntegrations_in);

  /**
   * @brief Construct a new Bounce Action Int object
   *
   * @param init_path is the initial path guess
   * @param TrueVacuumIn is the true vacuum candidate of the potential
   * @param FalseVacuumIn is the false vacuum
   * @param V is the class potential
   */
  BounceActionInt(std::vector<std::vector<double>> init_path_In,
                  std::vector<double> TrueVacuum_In,
                  std::vector<double> FalseVacuum_In,
                  std::function<double(std::vector<double>)> &V_In,
                  double T_In,
                  int MaxPathIntegrations_in);

  /**
   * @brief Used set the path of the class.
   *
   * @param init_path Knots that are used to describe the path
   */
  void setPath(std::vector<std::vector<double>> init_path_In);

  /**
   * @brief Precalculates dVdl and creates a spline with the result.
   * This is done to increase the runtime in large dimensional models
   *
   *  * @param l_start is the starting position produced in the
   * backwardspropagation part
   */
  void rasterize_dVdl(double l_start = 0);

  /**
   * @brief Prints a vector
   *
   * @param vec Vector to be printed
   */
  static void PrintVector(std::vector<double> vec);

  /**
   * @brief Numerical method to calculate the potential's (or other functions's)
   * gradient using finite differences method.
   *
   * This method is used while BSMPT is not able to
   * calculate the potential derivative analytically. We used the 4th order
   * method
   *
   * \f$\frac{\partial V}{\partial \phi_i} = \frac{1}{12
   * \epsilon}\left(-V(\dots ,\vec{\phi}_i + 2  \epsilon ) + 8 V(\dots
   * ,\vec{\phi}_i + \epsilon )- 8 V(\dots ,\vec{\phi}_i - \epsilon ) +
   * V(\dots ,\vec{\phi}_i - 2  \epsilon )\right)\f$
   *
   * where \f$ \epsilon \f$ is a small step.
   *
   * @param phi Where we want to calculate the gradient
   * @param V Potential (or other function)
   * @param eps Size of finite differences step
   * @param dim Dimensions of the VEV space (or dimensions of V argument)
   * @return std::vector<double> The \f$ dim \times 1 \f$ gradient of V taken at
   * phi
   */
  static std::vector<double>
  NablaNumerical(std::vector<double> phi0,
                 const std::function<double(std::vector<double>)> &V,
                 const double &eps,
                 const int &dim);
  /**
   * @brief Numerical method to calculate the potential's (or other functions's)
   * hessian matrix using finite differences method.
   *
   * \f$\frac{\partial^2 V}{\partial \phi_i \phi_j} = \frac{1}{4
   * \epsilon^2}\left(V(\dots, \vec{\phi}_i + \epsilon , \vec{\phi}_j +
   * \epsilon) - V(\dots, \vec{\phi}_i - \epsilon , \vec{\phi}_j +
   * \epsilon) - V(\dots, \vec{\phi}_i + \epsilon , \vec{\phi}_j -
   * \epsilon) + V(\dots, \vec{\phi}_i - \epsilon , \vec{\phi}_j -
   * \epsilon) \right)\f$
   *
   * where \f$ \epsilon \f$ is a small step.
   *
   * @param phi Where we want to calculate the Hessian matrix
   * @param V Potential (or other function)
   * @param eps Size of finite differences step
   * @param dim Dimensions of the VEV space (or dimensions of V argument)
   * @return std::vector<std::vector<double>> The \f$ dim \times \dim \f$
   *  hessian matrix of V taken at phi
   */
  static std::vector<std::vector<double>>
  HessianNumerical(std::vector<double> phi,
                   const std::function<double(std::vector<double>)> &V,
                   double eps,
                   const int &dim);

  /**
   * @brief Calculates \f$ \frac{dV}{dl} \f$ using the spline and potential
   * derivatives
   *
   * \f$ \frac{dV}{dl}  = \frac{\partial V}{\partial \vec{\phi}} \cdot  \frac{d
   * \vec{\phi}}{dl} = \nabla V(\vec{\phi}) \cdot  \frac{d
   * \vec{\phi}}{dl} \f$
   *
   * @param l Spline parameterization point where \f$ \frac{dV}{dl} \f$ is
   * calculated
   * @return double Returns \f$ \frac{dV}{dl} \f$
   */
  double Calc_dVdl(double l);

  /**
   * @brief Calculates \f$ \frac{d^2V}{dl^2} \f$ using the spline and potential
   * derivatives
   *
   * \f$ \frac{d^2V}{dl^2}  = \nabla V(\vec{\phi}) \cdot  \frac{d^2
   * \vec{\phi}}{dl^2} +  \left(\frac{d\vec{\phi}}{dl}\right)^T H(\vec{\phi})
   * \frac{d\vec{\phi}}{dl} \f$
   *
   * @return double Returns \f$ \frac{d^2V}{dl^2}  \f$
   */
  double Calc_d2Vdl2(double l);

  //
  /**
   * @brief Finds stationary points of a function (not only minimas).
   *
   * @param guess is the initial guess for the minimum
   * @param dV gradient of the function to be minimized
   * @param Hessian hessian of the function
   * @param error Maximum size of \f$ | \vec{dV} | \f$ that is considered a
   * minimum
   * @param const_multiplier If \f$ \det{Hessian} = 0\f$ this method does not
   * work. In that case we move the guess as \f$ \vec{p} \rightarrow \vec{p} -
   * const\_multiplier * \vec{dV}\f$
   * @param maxiter Maximum iteration exiting function
   * @return std::vector<double>
   */
  static std::vector<double> LocateMinimum(
      std::vector<double> guess,
      std::function<std::vector<double>(std::vector<double>)> &dV,
      std::function<std::vector<std::vector<double>>(std::vector<double>)>
          &Hessian,
      double error            = 1e-4,
      double const_multiplier = 0.01,
      int maxiter             = 1000);

  /**
   * @brief This functions uses @ref LocateMinimum to locate a minimum in the
   * potential
   *
   * @param guess of the minimum location.
   * @return std::vector<double>
   */
  std::vector<double> LocateMinimumPotential(std::vector<double> guess);

  // param x
  // return int

  /**
   * @brief Auxiliary function to locate the maximum absolute value of a vector.
   *
   * @param x is the vector.
   * @return int is the index of the maximum absolute value.
   */
  static int IndexMaximumAbsolute(std::vector<double> x);

  /**
   * @brief Calculates the classic BSMPT tunneling path
   *
   * Start with a straight line between @ref TrueVacuum and @ref FalseVacuum and
   * divide it into \f$ n \f$ base points \f$ \vec{b}_i \f$.
   *
   * At each basepoint \f$ \vec{b}_i \f$ draw a plane perpendicular to the
   * straight line path, i.e. the normal is TrueVacuum - FalseVacuum, and find
   * the potential minimum on that plane.
   *
   * @return std::vector<std::vector<double>> the classic BSMPT tunneling path
   */
  std::vector<std::vector<double>> BSMPTApproximateTunnelingPath();

  /**
   * @brief Calculates the classic BSMPT tunneling path using @ref
   * BSMPTApproximateTunnelingPath and sets the spline @ref path to it.
   *
   */
  void setBSMPTApproximateTunnelingPath();

  /**
   * @brief Auxiliary finction to compute the transpose of a matrix
   *
   * @param A is the \f$ n \times m \f$ matrix to be tranposed.
   * @return std::vector<std::vector<double>> \f$ m \times n \f$ transposed
   * matrix.
   */
  std::vector<std::vector<double>>
  transpose(std::vector<std::vector<double>> &A);

  /**
   * @brief Calculated the normal force \f$ \vec{N} \f$ on a @ref spline point.
   *
   * @param l is the @ref spline parameter where the force is calculated.
   * @param dldrho is \f$ \frac{dl}{d\rho}\f$.
   * @param gradient is the gradient evaluated at spline parameter l.
   * @return std::vector<double> is the \f$ \vec{N} \f$ at spline parameter l.
   */
  std::vector<double>
  NormalForce(double l, double dldrho, std::vector<double> gradient);

  /**
   * @brief Auxiliary function used in the Runge-Kutta 5th order
   * adaptative step @ref RK5_step.
   *
   * @param rho is the integration variable \f$ \rho \f$.
   * @param dvs are the functions values.
   * @param aks are used to return the function derivatives.
   */
  void
  AuxFunctionDev(double rho, std::vector<double> dvs, std::vector<double> &aks);

  /**
   * @brief Runge-Kutta 5th order step.
   *
   * Although the Runge-Kutta methods are valid for \f$ y' =
   * f(t,y) \f$ integration one can generalize the method for higher order ODE
   * using auxiliary functions. One can write \f$ y'' = f(t,y, y') \f$ as
   *
   * * \f$ y' = m \f$
   * * \f$ m' = f(t, y, m) \f$
   *
   * linearizing the ODE system allowing the application of Runge-Kutta method
   * to each one.
   *
   * This method takes a 5th order Runge-Kutta step that can be embedded into a
   * 4th order Runge-Kutta step. By comparing the result of the 4th and 5th
   * order we can control the error within our calculation by adjusting the step
   * size accordingly.
   *
   * @param y function values \f$ \{y(\rho_i), m(\rho_i) \} \f$.
   * @param dydx function derivatives \f$ \{y'(\rho_i), m'(\rho_i) \} \f$.
   * @param n is the number of linearizations, i.e. 2.
   * @param rho is the integration variable \f$ \rho \f$.
   * @param h is the step size.
   * @param yout is the 5th order Runge-Kutta integration result.
   * @param yerr is the difference between the 4th order and 5th order
   * Runge-Kutta result.
   */
  void RK5_step(std::vector<double> y,
                std::vector<double> dydx,
                int n,
                float rho,
                float h,
                std::vector<double> &yout,
                std::vector<double> &yerr);

  /**
   * @brief Modified Bessel function \f$I_\alpha (x) \f$ of the first kind
   *
   *  \f$I_\alpha (x) = \sum_{m=0}^\infty \frac{1}{m! \Gamma(m + \alpha + 1)}
   * \left(\frac{x}{2}\right)^{2m + \alpha}  \f$
   *
   * The convergence should be quite fast. We use as many terms to get the last
   * term to impact the result in \f$ 10^{-15}\f$ (relatively). The default
   * maximum number of terms is 50.
   *
   * @param alpha is an complex number.
   * @param x where to calculate the Bessel function.
   * @param terms are the maximum number of terms allowed in the calculation.
   * @return double returns \f$I_\alpha (x) \f$
   */
  static double BesselI(double alpha, double x, int terms = 100);

  /**
   * @brief Modified Bessel function i \f$J_1 (i x) \f$ of the first kind
   *
   *  \f$i J_1 (i x) = i \sum_{m=0}^\infty (-1)^{m + 1} \frac{1}{m! \Gamma(m +
   * \alpha + 1)} \left(\frac{x}{2}\right)^{2m + \alpha}  \f$
   *
   * The convergence should be quite fast. We use as many terms to get the last
   * term to impact the result in \f$ 10^{-15}\f$ (relatively). The default
   * maximum number of terms is 50.
   *
   * @param alpha is an complex number.
   * @param x where to calculate the Bessel function.
   * @param terms are the maximum number of terms allowed in the calculation.
   * @return double returns \f$I_\alpha (x) \f$
   */
  static double BesselJ(double x, int terms = 100);

  /**
   * @brief Integrates the 1D profile assuming \f$ \frac{dV}{dl} \f$ is a
   * constant. The solution is
   *
   * \f$ l(\rho) = l_0 + \left(\frac{1}{2(1 + \alpha)}
   * \frac{dV}{dl}\Big|_{l_0}\right) \rho^2 \f$
   *
   * @param l0 is the integration starting point
   * @param l  is the integration final point
   * @param dVdl \f$ \frac{dV}{dl} \f$ at \f$ l_0 \f$
   * @return std::vector<double>
   */
  std::vector<double> ExactSolutionCons(double l0, double l, double dVdl);

  /**
   * @brief Integrates the 1D profile assuming \f$ \frac{dV}{dl} \f$ is a
   * linear in l, i.e. \f$ \frac{dV}{dl} \approx dV + H (l - l_0) \f$. The
   * solution is for \f$ \alpha = 2 \f$ is
   *
   * \f$ l(\rho) = l_0 - \frac{dV}{H} +
   * \frac{dV}{H^{3/2}}\frac{\sinh(\rho\sqrt{H})}{\rho} \f$
   *
   * and the solution for  \f$ \alpha = 3 \f$ is
   *
   * \f$ l(\rho) = l_0 - \frac{dV}{H} + \frac{2 dV}{H^{3/2}} \frac{I_\alpha(\rho
   * \sqrt{H})}{\rho}\f$
   *
   * where \f$ I_\alpha(\rho) \f$ is the modified Bessel function of the first
   * kind.
   *
   * @param l0 is the integration starting point.
   * @param l is the integration final point.
   * @param dVdl \f$ \frac{dV}{dl} \f$ at \f$ l_0 \f$.
   * @param d2Vdl2 \f$ \frac{d^2V}{dl^2} \f$ at \f$ l_0 \f$.
   * @param maxiter Maximum iteration when calculating \f$ \rho(l) \f$.
   * @return std::vector<double>
   */
  std::vector<double>
  ExactSolutionLin(double l0, double l, double dVdl, double d2Vdl2);

  /**
   * @brief Integrates the 1D profile assuming \f$ \frac{dV}{dl} \f$ is a
   * linear in l, i.e. \f$ \frac{dV}{dl} \approx H (l - l_{min}) \f$. This
   * correspondes to a purely qudratic potential.
   *
   * @param l0 is the integration starting point.
   * @param l is the integration final point.
   * @return std::vector<double>
   */
  std::vector<double> ExactSolutionFromMinimum(double l0, double l);

  /**
   * @brief Calculates the 1D solution by comparing the @ref
   * ExactSolutionCons and @ref ExactSolutionLin so that the analytical step is
   * appropriate.
   *
   * @param l0 starting point
   * @return std::vector<double>
   */
  std::vector<double> ExactSolution(double l0);

  /**
   * @brief This is done to make sure that we can still find solution after path
   * deformation. This propagates the spline into negative values (by
   * extrapolating).
   *
   * @return double negative or zero value.
   */
  double BackwardsPropagation();

  /**
   * @brief Integrates 1D bounce equation once
   *
   * @param l0 is the starting value.
   * @param conv checks type of convergence. 0 = converged. -1 = undershot. +1 =
   * overshot.
   * @param rho vector of integration variable \f$ \rho \f$ steps.
   * @param l vector of variable \f$ l \f$ steps.
   * @param dl_drho vector of variable \f$ \frac{dl}{d\rho} \f$ steps.
   * @param d2l_drho2 vector of variable \f$ \frac{d^2l}{d\rho^2} \f$ steps.
   * @param maxiter is the maximum integration steps.
   * @param error is the acceptance of undershoot/overshoot.
   * @param eps_abs is used to control the step size error (RK4 vs RK5).
   * @param max_step in the case you want to set a maximum step size in \f$ \rho
   * \f$.
   */
  void IntegrateBounce(double l0,
                       int &conv,
                       std::vector<double> &rho,
                       std::vector<double> &l,
                       std::vector<double> &dl_drho,
                       std::vector<double> &d2l_drho2,
                       int maxiter,
                       double error,
                       double eps_abs,
                       double max_step = 0);
  /**
   * @brief Performs a binary search using the overshooting/undershooting method
   * to find the solution to the 1D bounce equation
   *
   * @param rho is the vector of \f$ \rho \f$ values from the solution
   * @param l is the vector of \f$ l \f$ values from the solution
   * @param dl_drho is the vector of \f$ \frac{dl}{d\rho} \f$ values from the
   * solution
   * @param d2l_drho2is the vector of \f$ \frac{d^2l}{d\rho^2} \f$ values from
   * the solution
   * @param error is the @ref IntegrateBounce error (acceptance of
   * undershoot/overshoot)
   * @param maxiter is the maximum number of binary searches
   */
  void Solve1DBounce(std::vector<double> &rho,
                     std::vector<double> &l,
                     std::vector<double> &dl_drho,
                     std::vector<double> &d2l_drho2,
                     double error = 1e-7,
                     int maxiter  = 100);

  /**
   * @brief Calculates the normalization of the force \f$ \vec{\phi} \rightarrow
   * \vec{\phi} + \vec{N}/reductor \f$. We have that \f$ reductor = \varepsilon
   * \max{\nabla V}/L \f$, where \f$ 10^{-4} \le \varepsilon \le 10^{-1} \f$ is
   * a small parameter
   *
   * @param MaximumGradient
   * @return double
   */
  double reductorCalculator(const double &MaximumGradient);

  /**
   * @brief Check if the force in each point is sufficient small compared to the
   * gradient in each point
   *
   * @param l list of spline parameter \f$ l \f$ of the knots
   * @param rho_l_spl list of \f$ \frac{dl}{d\rho}\f$ of the knots
   * @return true if converged
   * @return false if not converged
   */
  bool PathDeformationCheck(std::vector<double> &l, tk::spline &rho_l_spl);

  /**
   * @brief Takes a single path deformation step
   *
   * @param stepsize \f$ \varepsilon \f$
   * @param reductor is the reductor
   * @param l list of \f$ l \f$ at the knots of the old solution
   * @param rho_l_spl list of \f$ \frac{dl}{d\rho} \f$ at the knots of the old
   * solution
   * @param l_fornextpath list of new \f$ l \f$ at the new path iteration
   * @param best_path save the best path
   * @param saves the current iteration on the fly
   * @param MaximumGradient maximum \f$ \nabla V \f$
   * @param MaximumForce maximum \f$ \vec{N} \f$
   * @param MaximumRelativeError maximum \f$ \frac{|\vec{N}|}{|\nabla V|} \f$
   * @param Maximum_dldrho maximum \f$ \frac{dl}{d\rho} \f$
   * @param PerpendicularGradient maximum \f$ \nabla_\perp V \f$
   * @param inverseK inverse of the kernel matrix
   * @param forces list of forces to check if path is converging or not
   */
  void SinglePathDeformation(double &stepsize,
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
                             std::vector<std::vector<double>> &forces);

  /**
   * @brief Deforms the path minimizing the force \f$ \vec{N} \f$ without
   * solving
   *
   * @param l is the vector of \f$ \rho \f$ values from the solution
   * @param rho_l_spl tk::spline of \f$ \rho(l) \f$
   */
  void PathDeformation(std::vector<double> &l, tk::spline &rho_l_spl);
  /**
   * @brief Number of combinations of chosing k in n
   *
   * \f$ {n \choose k}  = \frac{n!}{k!(n-k)!}\f$
   *
   * @param n
   * @param k
   * @return unsigned
   */
  static unsigned nChoosek(unsigned n, unsigned k);

  /**
   * @brief Calculates the \f$ k^{th}\f$ Bernstein polynomial of degree \f$ n
   * \f$ at \f$ x \f$.
   *
   * \f$ B_{\nu,n}(x) = {n \choose \nu} x^\nu (1-x)^{n-\nu} \f$
   *
   * @param n Bernstein degree.
   * @param nu Bernstein basis function, \f$ 0 \leq \nu \leq n \f$.
   * @param x is the independent parameter.
   * @return double \f$ B_{\nu,n}(x) \f$.
   */
  static double Bernstein(int n, int nu, double x);

  /**
   * @brief Calculates the force vector \f$ \vec{N} \f$ of multiple path knots
   * at the same time. Use in the path deformation algorithm.
   *
   * @param dldrho is the vector of \f$ \frac{dl}{d\rho} \f$ values from the
   * solution.
   * @param gradient is the vector of \f$ \nabla(\vec{\phi}(l)) \f$ values from
   * the solution.
   * @param dphidl is the vector of \f$ \frac{d\vec{\phi}}{dl} \f$ values from
   * the solution.
   * @param d2phidl2 is the vector of \f$ \frac{d^2\vec{\phi}}{dl^2} \f$ values
   * from the solution.
   * @return std::vector<double> list of force vectors \f$ \vec{N} \f$ applied
   * to each path knot.
   */
  static std::vector<double>
  NormalForceBernstein(double dldrho,
                       std::vector<double> &gradient,
                       std::vector<double> &dphidl,
                       std::vector<double> &d2phidl2);

  static void ClearFolder(const char *path_in);

  /**
   * @brief Calculates \f$ \frac{d^2l}{d\rho^2} \f$
   *
   * \f$ \frac{d^2l}{d\rho^2} = \frac{dV}{dl} -
   * \frac{\alpha}{\rho}\frac{dl}{d\rho}\f$
   *
   * @param l \f$ = l(\rho) \f$
   * @param rho where we want to calculate.
   * @param dldrho \f$ = \frac{dl(\rho)}{d\rho} \f$
   * @return double \f$ \frac{d^2l}{d\rho^2} \f$ at \f$ \rho \f$
   */
  double d2ldrho2(double l, double rho, double dldrho);

  /**
   * @brief Calculate kinect term of the action
   *
   * @param rho list of rho coming from the integration
   * @param dl_drho_spl Spline of \f$\frac{dl(\rho)}{d\rho}\f$
   * @return double kinetic part of the action
   */
  double CalculateKineticTermAction(std::vector<double> &rho,
                                    tk::spline &dl_drho_spl);

  /**
   * @brief Calculate potential term of the action
   *
   * @param rho list of rho coming from the integration
   * @param l_rho_spl Spline of \f$l(\rho)\f$
   * @return double potential part of the action
   */
  double CalculatePotentialTermAction(std::vector<double> &rho,
                                      tk::spline &l_rho_spl);

  /**
   * @brief Calculates the action of the bounce equation by deforming the
   * given
   * @ref path and minimizing the normal force \f$ \vec{N} \f$ until it gets
   * sufficiently small.
   *
   * @param error is the acceptance of undershoot/overshoot method.
   */
  void CalculateAction(double error = 1e-6);
};
} // namespace BSMPT