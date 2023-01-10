// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/baryo_calculation/Fluid_Type/top_source.h>

/**
 * @file
 */

namespace BSMPT
{
namespace Baryo
{

typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

double top_source::Calc_nL(double z_start, double z_end) const
{
  /*
  omega[0] -> q
  omega[1] -> t
  omega[2] -> h1
  omega[3] -> h2
  omega[4] -> q_prime
  omega[5] -> t_prime
  omega[6] -> h1_prime
  omega[7] -> h2_prime
  */
  state_type mu(8);
  mu                    = {0, 0, 0, 0, 0, 0, 0, 0};
  const double C_AbsErr = 1e-9;
  const double C_RelErr = 1e-5;
  double stepsize_initial{0};
  if (z_start < z_end)
    stepsize_initial = 1e-8;
  else
    stepsize_initial = -1e-8;
  double abs_err = C_AbsErr;
  double rel_err = C_RelErr;
  integrate_adaptive(make_controlled(abs_err, rel_err, error_stepper_type()),
                     *this,
                     mu,
                     z_start,
                     z_end,
                     stepsize_initial);

  return 3 * mu[0] +
         2 * mu[1]; // as defined in 1811.11104; used q1=-2b and b = -(q+t)
}

void top_source::operator()(const state_type &omega,
                            state_type &domega,
                            const double z)
{
  /*
      Definition of all transport coefficients
  */
  std::vector<double> quark_mass;
  std::vector<double> quark_mass_prime;
  top_func(z, quark_mass, quark_mass_prime);
  double mt        = quark_mass[0]; // top quark mass
  double sym_phase = gen_fluid::symmetric_CP_violating_phase;
  double brk_phase = gen_fluid::broken_CP_violating_phase;
  auto theta_vec   = Calc_theta(z, sym_phase, brk_phase);
  // double theta    =   theta_vec[0];//phase factor of the top quark
  double theta_prime =
      theta_vec[1]; // derivative of the phase factor of the top quark
  // Calculation of kappa_t
  Calc_kappa_obj.set_class(Temp, mt);
  double num_int  = NIntegrate_kappa(Calc_kappa_obj);
  double kappa_q  = kappa_QL_0 * num_int; // left-handed kappa
  double kappa_tR = kappa_QR_0 * num_int; // right-handed kappa

  // Rescaled chemical potential like in 1811.11104
  double mu_M = (omega[1] / kappa_tR - omega[0] / kappa_q);
  double bR   = -(omega[1] + omega[0]); // local baryon conservation used to
                                        // expresse the right-handed b density
  double mu_SS = 2 * omega[0] / kappa_q - omega[1] / kappa_tR - bR / kappa_QR_0;
  double mu_Y  = omega[1] / kappa_tR - omega[0] / kappa_q -
                omega[2] / kappa_H_0 - omega[3] / kappa_H_0;

  Calc_Gam_obj.set_class(Temp, vw, mt, msqrt_thermal_top, dmsqrt_thermal_top);
  double Gam_M = Nintegrate_GamM(Calc_Gam_obj);
  Calc_Scp_obj.set_class(
      Temp, vw, mt, theta_prime, msqrt_thermal_top, dmsqrt_thermal_top);
  double Scp = Nintegrate_Scp(Calc_Scp_obj);

  /*
      omega[0] -> q
      omega[1] -> t
      omega[2] -> h1
      omega[3] -> h2
      omega[4] -> q_prime
      omega[5] -> t_prime
      omega[6] -> h1_prime
      omega[7] -> h2_prime
      */
  domega[0] = omega[4];
  domega[1] = omega[5];
  domega[2] = omega[6];
  domega[3] = omega[7];
  domega[4] =
      (Scp - Gam_M * mu_M + 2 * Gam_SS * mu_SS - Gam_Y * mu_Y + vw * omega[4]) /
      Dq;
  domega[5] = (-Scp + Gam_M * mu_M + Gam_Y * mu_Y + vw * omega[5]) / Dt;
  domega[6] = (-Gam_Y * mu_Y + vw * omega[6]) / Dh;
  domega[7] = (-Gam_Y * mu_Y + vw * omega[7]) / Dh;
}

} // namespace Baryo
} // namespace BSMPT
