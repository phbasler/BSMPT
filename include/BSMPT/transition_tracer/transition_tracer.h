// Copyright (C) 2024 Lisa Biermann, Margarete Mühlleitner, Rui Santos, João
// Viana
//
// SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

/**
 * @file transition history evaluator
 */

#include "BSMPT/bounce_solution/bounce_solution.h" // BounceSolution
#include "BSMPT/gravitational_waves/gw.h"          // GravitationalWaves
#include "BSMPT/minimum_tracer/minimum_tracer.h"   // MinimumTracer

namespace BSMPT
{

/**
 * @brief user_input struct to store user input and distribute to the classes
 * @param modelPointer model pointer
 * @param T_low lowest temperature, default: 0 GeV
 * @param T_high highest temperature, default: 300 GeV
 * @param vwall wall velocity, default: 0.95
 * @param perc_prbl false vacuum fraction at percolation temperature, default:
 * 71%
 * @param compl_prbl false vacuum fraction at completion temperature, default:
 * 1%
 * @param epsturb epsilon value of turbulence contribution, default: 0.1
 * @param maxpathintegrations maximal number of path integrations, default: 7
 * @param multistepmode choose multi-step PT modes: default (= -1), 0, 1, 2,
 * auto (= 3)
 * @param num_points number of equally-spaced intermediate points to check,
 * default: 10
 * @param ewsr_check check of electroweak symmetry restoration, default: off (=
 * 0)
 * @param nlo_check check of nlo stability, default: on (= 1)
 * @param use_gsl whether GSL minimizer is used
 * @param use_cmaes whether CMAES minimizer is used
 * @param use_nlopt whether NLopt minimizer is used
 * @param which_minimizer which minimizers are used
 * @param use_multithreading whether multi-threading is used
 * @param gw_calculation bool to turn GW parameter calculation on/off
 * @param which_transition_temp which transition temperature is chosen: 1 =
 * nucl_approx, 2 = nucl, 3 = perc (default), 4 = compl
 */
struct user_input
{
  std::shared_ptr<Class_Potential_Origin> modelPointer;
  double T_low            = 0;
  double T_high           = 300;
  double vwall            = 0.95;
  double perc_prbl        = 0.71;
  double compl_prbl       = 0.01;
  double epsturb          = 0.1;
  int maxpathintegrations = 7;
  int multistepmode       = -1;
  int num_points          = 10;
  int ewsr_check          = 0;
  int nlo_check           = 1;

  int which_minimizer     = Minimizer::WhichMinimizerDefault;
  bool use_multithreading = false;

  bool gw_calculation       = false;
  int which_transition_temp = 3;
};

/**
 * @brief status codes struct
 */
struct status_codes
{
  std::string status_nlo_stability = "not_set";
  std::string status_ewsr          = "not_set";
  std::string status_tracing       = "not_set";
  std::string status_coex_pairs    = "not_set";
  // index of vectors is coex_phase_id
  std::vector<std::string> status_crit;
  std::vector<std::string> status_bounce_sol;
  std::vector<std::string> status_nucl_approx;
  std::vector<std::string> status_nucl;
  std::vector<std::string> status_perc;
  std::vector<std::string> status_compl;
};

/**
 * @brief transition data struct
 */
struct transition_data
{
  double crit_temp        = NAN;
  double nucl_approx_temp = NAN;
  double nucl_temp        = NAN;
  double perc_temp        = NAN;
  double compl_temp       = NAN;

  std::vector<double> crit_true_vev;
  std::vector<double> crit_false_vev;
  std::vector<double> nucl_approx_true_vev;
  std::vector<double> nucl_approx_false_vev;
  std::vector<double> nucl_true_vev;
  std::vector<double> nucl_false_vev;
  std::vector<double> perc_true_vev;
  std::vector<double> perc_false_vev;
  std::vector<double> compl_true_vev;
  std::vector<double> compl_false_vev;
};

/**
 * @brief gravitational wave data struct
 */
struct gw_data
{
  double vwall = NAN;

  double alpha       = NAN;
  double beta_over_H = NAN;

  double K_sw   = NAN;
  double K_turb = NAN;

  double fpeak_sw     = NAN;
  double fpeak_turb   = NAN;
  double h2Omega_sw   = NAN;
  double h2Omega_turb = NAN;

  double SNR_sw   = NAN;
  double SNR_turb = NAN;
  double SNR      = NAN;

  std::string status_gw = "not_set";
  double trans_temp     = NAN;
};

struct output
{
  std::vector<std::string> legend;
  status_codes status;
  std::vector<transition_data> vec_trans_data;
  std::vector<gw_data> vec_gw_data;
  std::size_t num_coex_phase_pairs = 0;

  // stores string-identifier of phase id's throughout the transition history of
  // the universe
  std::string transition_history = "not_set";
};

class TransitionTracer
{
protected:
  /**
   * @brief modelPointer for the used parameter point
   */
  std::shared_ptr<Class_Potential_Origin> modelPointer;

private:
  /**
   * @brief number of VEVs of model
   */
  std::size_t num_vev;

  /**
   * @brief vector of all found coexisting phase regions
   */
  std::vector<CoexPhases> vec_coex;

public:
  /**
   * @brief TransitionTracer constructor
   * @param input user input
   */
  TransitionTracer(user_input &input);
  ~TransitionTracer();

  /**
   * @brief Store the list of bounce solutions
   *
   */
  std::vector<BounceSolution> ListBounceSolution;

  /**
   * @brief output data storage
   */
  output output_store;
};

} // namespace BSMPT