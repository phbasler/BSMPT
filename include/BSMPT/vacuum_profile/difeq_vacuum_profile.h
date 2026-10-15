// https://github.com/NJdevPro/Numerical-Recipes

#pragma once

#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/relaxation/difeq.h>
#include <BSMPT/utility/utility.h>
#include <fstream>
#include <gsl/gsl_spline.h>

using BSMPT::Delta;

namespace BSMPT
{
namespace VacuumProfileNS
{
/**
 * @brief Modes of solving the field equation
 *
 */
enum class ProfileSolverMode
{
  Field,
  Deriv
};

struct Difeq_VacuumProfile : Difeq
{
  /**
   * @brief How to solve the equations
   *
   */
  const ProfileSolverMode &mode;
  /**
   * @brief Dimension of VEV space
   *
   */
  const size_t &dim;
  /**
   * @brief List of the z grid of points
   *
   */
  std::vector<double> &z;
  /**
   * @brief True and False Vacuum
   *
   */
  const std::vector<double> &TrueVacuum, &FalseVacuum;
  /**
   * @brief Potential
   *
   */
  const std::function<double(std::vector<double>)> &V;
  /**
   * @brief Gradient
   *
   */
  const std::function<std::vector<double>(std::vector<double>)> &dV;
  /**
   * @brief Hessian
   *
   */
  const std::function<std::vector<std::vector<double>>(std::vector<double>)>
      &Hessian;
  /**
   * @brief Friction coefficient
   * - 0 if V(True) = V(False) -> Domain Wall
   * - > 0 if V(True) < V(False) -> Bubble Wall
   */
  double mu = 0;
  /**
   * @brief V(False) - V(True)
   *
   */
  double DeltaV;

  Difeq_VacuumProfile(
      // Solver mode
      const ProfileSolverMode &mode_In,
      // Dimension of VEV space
      const size_t &dim_In,
      // z list
      std::vector<double> &z_In,
      // True
      const std::vector<double> &TrueVacuum_In,
      // False
      const std::vector<double> &FalseVacuum_In,
      // Potential
      const std::function<double(std::vector<double>)> &V_In,
      // Gradient
      const std::function<std::vector<double>(std::vector<double>)> &dV_In,
      // Hessian
      const std::function<std::vector<std::vector<double>>(std::vector<double>)>
          &Hessian_In)
      : mode(mode_In)
      , dim(dim_In)
      , z(z_In)
      , TrueVacuum(TrueVacuum_In)
      , FalseVacuum(FalseVacuum_In)
      , V(V_In)
      , dV(dV_In)
      , Hessian(Hessian_In)
      , DeltaV(V(FalseVacuum) - V(TrueVacuum))
  {
  }

  void calc_mu(MatDoub &y)
  {
    const size_t n = y.cols();

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline    = gsl_spline_alloc(gsl_interp_cspline, n);

    std::vector<double> zz, phigrad;
    for (size_t k = 0; k < y.cols(); k++)
    {
      zz.push_back(z[k]);
      double t = 0;
      for (size_t j = 0; j < dim; j++)
        t += pow(y[j][k], 2);
      phigrad.push_back(t);
    }

    gsl_spline_init(spline, zz.data(), phigrad.data(), n);
    double EnergyDissipated =
        gsl_spline_eval_integ(spline, zz.front(), zz.back(), acc);

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    mu = DeltaV / EnergyDissipated;
  }

  void smatrix(const int k,
               const int k1,
               const int k2,
               const int jsf,
               VecInt &indexv,
               MatDoub &s,
               MatDoub &y)
  {
    s.zero(); // Set matrix s = 0
    if (k == k1)
    {
      // Calculate mu before anything else
      calc_mu(y);
      if (mode == ProfileSolverMode::Deriv)
      {
        // Boundary conditions dv/dz = 0 on first boundary
        for (size_t field = 0; field < dim; field++)
        {
          // Sn at the first boundary
          s[dim + field][2 * dim + field] = 1.0;
          // B0
          s[dim + field][jsf] = y[field][0];
        }
      }
      else if (mode == ProfileSolverMode::Field)
      {
        for (size_t field = 0; field < dim; field++)
        {
          // Sn at the first boundary
          s[dim + field][3 * dim + field] = 1.0;
          // B0
          s[dim + field][jsf] = y[dim + field][0] - TrueVacuum[field];
        }
      }
    }
    else if (k > k2 - 1)
    {
      if (mode == ProfileSolverMode::Deriv)
      {
        // Boundary conditionsd dv/dz = 0 on second boundary
        for (size_t field = 0; field < dim; field++)
        {
          // Sn at the last boundary
          s[field][2 * dim + field] = 1.0;
          // C0
          s[field][jsf] = y[field][k2 - 1];
        }
      }
      else if (mode == ProfileSolverMode::Field)
      {
        // Fix the fields at |z| -> Infinity
        for (size_t field = 0; field < dim; field++)
        {
          // Sn at the last boundary
          s[field][3 * dim + field] = 1.0;
          // C0
          s[field][jsf] = y[dim + field][k2 - 1] - FalseVacuum[field];
        }
      }
    }
    else
    {
      const double dz = z[k] - z[k - 1];
      std::vector<double> point; // middle point (phi(k) + phi(k-1)) / 2
      for (size_t j = 0; j < dim; j++)
        point.push_back((y[dim + j][k] + y[dim + j][k - 1]) / 2);
      const std::vector<double> dv                   = dV(point);
      const std::vector<std::vector<double>> hessian = Hessian(point);

      // 0 <= j < dim -> phi'(z)
      for (size_t j = 0; j < dim; j++)
      {
        for (size_t n = 0; n < dim; n++)
        {
          s[j][0 * dim + n] = -Delta(j, n) * (1 - dz * mu / 2); // 1
          s[j][1 * dim + n] = -1. / 2. * dz * hessian[j][n];    // 3
          s[j][2 * dim + n] = Delta(j, n) * (1 + dz * mu / 2);  // 5
          s[j][3 * dim + n] = -1. / 2. * dz * hessian[j][n];    // 7
        }
        //  Equations for E(k,k-1)
        s[j][jsf] = y[j][k] - y[j][k - 1] -
                    dz * (dv[j] - mu * (y[j][k] + y[j][k - 1]) / 2);
      }

      // dim <= j < 2 * dim -> phi(z)
      for (size_t j = 0; j < dim; j++)
      {
        for (size_t n = 0; n < dim; n++)
        {
          s[dim + j][0 * dim + n] = -dz / 2. * Delta(j, n); // 2
          s[dim + j][1 * dim + n] = -Delta(j, n);           // 4
          s[dim + j][2 * dim + n] = -dz / 2. * Delta(j, n); // 6
          s[dim + j][3 * dim + n] = Delta(j, n);            // 8
        }
        //  Equations for E(k,k-1)
        s[dim + j][jsf] = y[dim + j][k] - y[dim + j][k - 1] -
                          dz * (y[j][k] + y[j][k - 1]) / 2;
      }
    }
    // Reorder columns
    MatDoub ss = s;
    for (size_t j = 0; j < 2 * dim; j++)
    {
      for (size_t n = 0; n < 2 * dim; n++)
      {
        s[j][n]           = ss[j][indexv[n]];
        s[j][2 * dim + n] = ss[j][2 * dim + indexv[n]];
      }
    }
  }
};
} // namespace VacuumProfileNS
} // namespace BSMPT