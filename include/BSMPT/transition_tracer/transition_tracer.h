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
 * @param which_minimizer which minimizers are used
 * @param use_multithreading whether multi-threading is used
 * @param gw_calculation bool to turn GW parameter calculation on/off
 * @param which_transition_temp which transition temperature is chosen
 * @param PNLO_scaling pressure scaling at NLO, 1 -> N processes at bubble
 * wall
 * @param number_of_initial_scan_temperatures number of temperature steps in the
 * initial scan of the bounce solver
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

  bool gw_calculation = false;
  TransitionTemperature which_transition_temp =
      TransitionTemperature::Percolation;
  int PNLO_scaling                           = 1;
  size_t number_of_initial_scan_temperatures = 25;
};

/**
 * @brief status codes struct
 */
struct status_codes
{
  StatusNLOStability status_nlo_stability = StatusNLOStability::NotSet;
  StatusEWSR status_ewsr                  = StatusEWSR::NotSet;
  StatusTracing status_tracing            = StatusTracing::NotSet;
  StatusCoexPair status_coex_pairs        = StatusCoexPair::NotSet;
  // index of vectors is coex_phase_id
  std::vector<StatusCrit> status_crit;
  std::vector<StatusGW> status_bounce_sol;
  std::vector<StatusTemperature> status_nucl_approx;
  std::vector<StatusTemperature> status_nucl;
  std::vector<StatusTemperature> status_perc;
  std::vector<StatusTemperature> status_compl;
};

/**
 * @brief transition data struct
 */
struct transition_data
{
  std::optional<double> crit_temp;
  std::optional<double> nucl_approx_temp;
  std::optional<double> nucl_temp;
  std::optional<double> perc_temp;
  std::optional<double> compl_temp;

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
  std::optional<double> vwall;

  std::optional<double> alpha;
  std::optional<double> beta_over_H;

  std::optional<double> kappa_col;
  std::optional<double> kappa_sw;
  std::optional<double> Epsilon_Turb;
  std::optional<double> cs_f;
  std::optional<double> cs_t;

  std::optional<double> fb_col;
  std::optional<double> omegab_col;

  std::optional<double> f1_sw;
  std::optional<double> f2_sw;
  std::optional<double> omega_2_sw;

  std::optional<double> f1_turb;
  std::optional<double> f2_turb;
  std::optional<double> omega_2_turb;

  std::optional<double> SNR_col;
  std::optional<double> SNR_sw;
  std::optional<double> SNR_turb;
  std::optional<double> SNR;

  StatusGW status_gw = StatusGW::NotSet;
  std::optional<double> trans_temp;
  std::optional<double> reh_temp;
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

public:
  /**
   * @brief vector of all found coexisting phase regions
   */
  std::vector<CoexPhases> vec_coex;

  /**
   * @brief TransitionTracer constructor
   * @param input user input
   */
  TransitionTracer(user_input &input);
  ~TransitionTracer();

  /**
   * @brief Store the list of bounce solutions
   */
  std::vector<BounceSolution> ListBounceSolution;

  /**
   * @brief output data storage
   */
  output output_store;

  /**
   * @brief CheckMassRatio Prints the scalar and gauge boson squard masses over
   * temperature squared ratios at the given point and temperature and returns
   * the maximal value
   * @param input
   * @param vev
   * @param temp
   * @return maximal ratio
   */
  double CheckMassRatio(const user_input &input,
                        const std::vector<double> &vec,
                        const double &temp) const;
};

} // namespace BSMPT
