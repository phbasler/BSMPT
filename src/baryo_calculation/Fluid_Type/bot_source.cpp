// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/baryo_calculation/Fluid_Type/bot_source.h>
#include <BSMPT/utility/Logger.h>

/**
 * @file
 */

namespace BSMPT
{
namespace Baryo
{

typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

void bot_source::operator()(const state_type &omega,
                            state_type &domega,
                            const double z)
{
  /*
      omega[0] -> q
      omega[1] -> t
      omega[2] -> b
      omega[3] -> h1
      omega[4] -> h2
      omega[5] -> u
      omega[6]  -> q_prime
      omega[7]  -> t_prime
      omega[8]  -> b_prime
      omega[9]  -> h1_prime
      omega[10] -> h2_prime
      omega[11] -> u_prime
  */

  /*
      Definition of all transport coefficients
  */
  std::vector<double> quark_mass;
  std::vector<double> quark_mass_prime;
  // TOP and BOT quark mass calculation
  top_func(z, quark_mass, quark_mass_prime);
  double mt = quark_mass[0];
  double mb{0};
  // BOTTOM quark mass can be set to zero to have a crosscheck
  if (bot_mass_flag == 1)
    mb = quark_mass[1];
  else if (bot_mass_flag == 2)
    mb = 0;
  else
    throw std::runtime_error("No valid option for the bottom mass @ ()");

  // Phase Calculation
  // top
  auto theta_vec_top = Calc_theta(z,
                                  gen_fluid::TOP_symmetric_CP_violating_phase,
                                  gen_fluid::TOP_broken_CP_violating_phase);
  // double theta_top        =   theta_vec_top[0];
  double theta_prime_top = theta_vec_top[1];
  // bot
  auto theta_vec_bot = Calc_theta(z,
                                  gen_fluid::BOT_symmetric_CP_violating_phase,
                                  BOT_broken_CP_violating_phase);
  // double theta_bot        =   theta_vec_bot[0];
  double theta_prime_bot = theta_vec_bot[1];
  // TOP statistical factor
  Calc_kappa_obj.set_class(Temp, mt);
  double num_int  = NIntegrate_kappa(Calc_kappa_obj);
  double kappa_tL = kappa_QL_0 * num_int;
  double kappa_tR = kappa_QR_0 * num_int;
  // BOT statistical factor
  Calc_kappa_obj.set_class(Temp, mb);
  num_int         = NIntegrate_kappa(Calc_kappa_obj);
  double kappa_bL = kappa_QL_0 * num_int;
  double kappa_bR = kappa_QR_0 * num_int;
  // Effective statistical factor
  double kappa_q = kappa_tL * kappa_bL / (kappa_tL + kappa_bL);

  // Rescaled chemical potential like in 1811.11104
  double mu_M_t = omega[1] / kappa_tR - omega[0] / kappa_q;
  double mu_M_b = omega[2] / kappa_bR - omega[0] / kappa_q;
  double mu_Y_t = omega[1] / kappa_tR - omega[0] / kappa_q -
                  omega[3] / kappa_H_0 - omega[4] / kappa_H_0;
  double mu_Y_b = omega[2] / kappa_bR - omega[0] / kappa_q -
                  omega[3] / kappa_H_0 - omega[4] / kappa_H_0;
  double mu_SS = -4 * omega[5] * (2 / kappa_QL_0 + 1 / kappa_QR_0) +
                 2 * omega[0] / kappa_q - omega[1] / kappa_tR -
                 omega[2] / kappa_bR;
  Calc_Gam_obj.set_class(Temp, vw, mt, msqrt_thermal_top, dmsqrt_thermal_top);
  double Gam_M_t = Nintegrate_GamM(Calc_Gam_obj);
  Calc_Gam_obj.set_class(Temp, vw, mb, msqrt_thermal_bot, dmsqrt_thermal_bot);
  double Gam_M_b = Nintegrate_GamM(Calc_Gam_obj);

  // Gam_Y_b has to be zero if the bottom mass vanishes
  if (bot_mass_flag == 2) Gam_Y_b = 0;

  Calc_Scp_obj.set_class(
      Temp, vw, mt, theta_prime_top, msqrt_thermal_top, dmsqrt_thermal_top);
  double Scp_t = Nintegrate_Scp(Calc_Scp_obj);
  double Scp_b = 0;
  if (bot_mass_flag == 1)
  {
    Calc_Scp_obj.set_class(
        Temp, vw, mb, theta_prime_bot, msqrt_thermal_bot, dmsqrt_thermal_bot);
    Scp_b = Nintegrate_Scp(Calc_Scp_obj);
  }
  if (bot_mass_flag == 2)
  {
    Scp_b = 0;
  }
  if ((bot_mass_flag != 1) and (bot_mass_flag != 2))
  {
    Logger::Write(LoggingLevel::EWBGDetailed,
                  "bot_mass_flag = " + std::to_string(bot_mass_flag));
    throw std::runtime_error("No valid bot_mass_flag @ operator()");
  }

  domega[0] = omega[6];
  domega[1] = omega[7];
  domega[2] = omega[8];
  domega[3] = omega[9];
  domega[4] = omega[10];
  domega[5] = omega[11];

  /*
      dmu qmu = vw qprime - Dq q_drpime = C
      --> q_dprime = ( vw qprime - C)/Dq
  */

  domega[6] = (vw * omega[6] -
               (Gam_M_t * mu_M_t + Gam_M_b * mu_M_b + Gam_Y_t * mu_Y_t +
                Gam_Y_b * mu_Y_b - 2 * Gam_SS * mu_SS - Scp_t - Scp_b)) /
              Dq;
  domega[7] = (vw * omega[7] - (-Gam_M_t * mu_M_t - Gam_Y_t * mu_Y_t +
                                Gam_SS * mu_SS + Scp_t)) /
              Dq;
  domega[8] = (vw * omega[8] - (-Gam_M_b * mu_M_b - Gam_Y_b * mu_Y_b +
                                Gam_SS * mu_SS + Scp_b)) /
              Dq;
  domega[9]  = (vw * omega[9] - (Gam_Y_t * mu_Y_t - Gam_Y_b * mu_Y_b)) / Dh;
  domega[10] = (vw * omega[10] - (Gam_Y_t * mu_Y_t - Gam_Y_b * mu_Y_b)) / Dh;
  domega[11] = (vw * omega[11] - Gam_SS * mu_SS) / Dq;
}

double bot_source::Calc_nL(double z_start, double z_end) const
{
  /*
      omega[0] -> q
      omega[1] -> t
      omega[2] -> b
      omega[3] -> h1
      omega[4] -> h2
      omega[5] -> u
      omega[6]  -> q_prime
      omega[7]  -> t_prime
      omega[8]  -> b_prime
      omega[9]  -> h1_prime
      omega[10] -> h2_prime
      omega[11] -> u_prime
  */
  state_type mu(12);
  mu                    = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  const double C_AbsErr = 1e-9;
  const double C_RelErr = 1e-5;
  double stepsize_initial{0};
  if (z_start < z_end) stepsize_initial = 1e-8;
  if (z_start > z_end) stepsize_initial = -1e-8;
  double abs_err = C_AbsErr;
  double rel_err = C_RelErr;
  integrate_adaptive(make_controlled(abs_err, rel_err, error_stepper_type()),
                     *this,
                     mu,
                     z_start,
                     z_end,
                     stepsize_initial);

  return mu[0] - 2 * mu[5]; // Additional left-handed up-type quarks;  q1 = -2 u
}

} // namespace Baryo
} // namespace BSMPT
