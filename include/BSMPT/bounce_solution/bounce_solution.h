// Copyright (C) 2024 Lisa Biermann, Margarete Mühlleitner, Rui Santos, João
// Viana
//
// SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

/**
 * @file bounce solution calculation
 */

#include <BSMPT/bounce_solution/action_calculation.h>
#include <BSMPT/minimum_tracer/minimum_tracer.h> // MinimumTracer
#include <BSMPT/models/SMparam.h>
#include <BSMPT/utility/spline/spline.h>
#include <Eigen/Dense>
#include <algorithm>             // std::max
#include <gsl/gsl_deriv.h>       // numerical derivative
#include <gsl/gsl_integration.h> // numerical integration

namespace BSMPT
{

/**
 * struct to store result and error
 */
struct resultErrorPair
{
  double result;
  double error;
};

/**
 * @brief BounceSolution class that handles the calculation of the bounce
 * solution as well as the calculation of the charateristic temperature scales
 */
class BounceSolution
{
public:
  /**
   * @brief modelPointer for the used parameter point
   */
  std::shared_ptr<Class_Potential_Origin> modelPointer;

  /**
   * @brief MinTracer object
   */
  std::shared_ptr<MinimumTracer> MinTracer;

  /**
   * @brief epsilon of turbulence efficiency factor
   */
  double epsturb = 0.1;

  /**
   * @brief wall velocity
   */
  double vwall = 0.95;

  /**
   * @brief number of temperature steps in the initial scan of the bounce solver
   *
   */
  size_t NumberOfInitialScanTemperatures;

  /**
   * @brief set to true if nucleation temperature is set
   */
  bool nucleation_temp_set = false;

  /**
   * @brief set to true if percolation temperature is set
   */
  bool percolation_temp_set = false;

  /**
   * @brief set to true if completion temperature is set
   */
  bool completion_temp_set = false;

  /**
   * @brief critical temperature/highest temperature when transition can occur
   */
  double Tc = -1;

  /**
   * @brief lowest temperature when a transition can occur
   */
  double Tm = -1;

  /**
   * @brief nucleation temperature
   */
  double Tnucl = -1;

  /**
   * @brief approximate nucleation temperature
   */
  double Tnucl_approx = -1;

  /**
   * @brief percolation temperature
   */
  double Tperc = -1;

  /**
   * @brief completion temperature
   */
  double Tcompl = -1;

  /**
   * @brief stored temperature
   */
  double store_Temp;

  /**
   * @brief PT strength
   */
  double alpha = -1;

  /**
   * @brief Inverse time scale \f$ \frac{\beta}{H} \f$
   *
   */
  double betaH = -1;

  /**
   * @brief number of effective degrees of freedom
   */
  double gstar;

  /**
   * @brief index of the true vacuum phase candidate in the coex list
   *
   */
  int indexTrueCandidatePhase;

  /**
   * @brief spline used to interpolate the action as a function of the
   * temperature
   *
   */
  tk::spline S3ofT_spline;

  /**
   * @brief Set of BounceActionInt objects with valid solutions.
   *
   */
  std::vector<BounceActionInt> SolutionList;

  /**
   * @brief List of group elements allowed by the potential
   *
   */
  std::vector<Eigen::MatrixXd> GroupElements;

  /**
   * @brief Store symmetry that produces the best tunneling rate
   *
   */
  Eigen::MatrixXd OptimalDiscreteSymmetry;

  /**
   * @brief Transforms the vev to the optimal vev.
   *
   * @param vev to be converted
   * @return std::vector<double> transformed vev
   */
  std::vector<double>
  TransformIntoOptimalDiscreteSymmetry(const std::vector<double> &vev);

  /**
   * @brief Calculates which is the optimal symmetry from the group of
   * symmetries
   *
   */
  void CalculateOptimalDiscreteSymmetry();

  /**
   * @brief Storage of the tunneling rate per volume of the transition from
   * false to true vacuum
   * @param Temp temperature
   */
  double TunnelingRate(const double &Temp);
  /**
   * @brief Storage of the temperature-dependent Hubble rate
   * @param Temp temperature
   */
  double HubbleRate(const double &Temp);

  /**
   * @brief inner_integrand friend to define inner integrand of percolation
   * temperature integral
   */
  friend double inner_integrand(double var, void *params);

  /**
   * @brief outer_integrand friend to define outer integrand of percolation
   * temperature integral
   */
  friend double outer_integrand(double var, void *params);

  /**
   * @brief Calculate euclidian action at temperature T
   */
  double GetBounceSol(const double &Temp) const;

  /**
   * @brief action_ratio friend to define input of numerical derivative in
   * calculation of inverse time scale
   */
  friend double action_ratio(double var, void *params);

public:
  /**
   * @brief AbsErr absolute error for numerical integration
   */
  const double AbsErr = 0;

  /**
   * @brief RelErr relative error for numerical integration
   */
  const double RelErr = 1e-6;

  /**
   * @brief pair of coexisiting phases
   */
  CoexPhases phase_pair;

  /**
   * @brief status of bounce solver
   *
   */
  StatusGW status_bounce_sol = StatusGW::NotSet;

  /**
   * @brief status of approximate nucleation temperature calculation
   *
   */
  BSMPT::StatusTemperature status_nucl_approx =
      BSMPT::StatusTemperature::NotSet;

  /**
   * @brief status of nucleation temperature calculation
   *
   */
  BSMPT::StatusTemperature status_nucl = BSMPT::StatusTemperature::NotSet;

  /**
   * @brief status of percolation temperature calculation
   *
   */
  BSMPT::StatusTemperature status_perc = BSMPT::StatusTemperature::NotSet;

  /**
   * @brief status of completion temperature calculation
   *
   */
  BSMPT::StatusTemperature status_compl = BSMPT::StatusTemperature::NotSet;

  /**
   * @brief \f$ v_{\text{wall}}\f$ defined by the user as an input parameter.
   *  If \f$ v_{\text{wall}}\f = -1$ then we use the approximation coming from
   * https://arxiv.org/abs/2210.16305
   * If \f$ v_{\text{wall}}\f = -2$ then we use the upper bound from
   * https://arxiv.org/abs/2305.02357
   *
   */
  double UserDefined_vwall = 0.95;

  /**
   * @brief Number of integration of the bounce
   *
   */
  int MaxPathIntegrations = 7;

  /**
   * @brief Construct a new Bounce Sol Calc object. Used for testing
   * @param pointer_in model pointer
   */
  BounceSolution(const std::shared_ptr<Class_Potential_Origin> &pointer_in);

  /**
   * @brief Construct a new Bounce Sol Calc object. This class takes as input a
   * pair of coexisting phases and delegates to constructor with provided
   * symmetry group.
   *
   * @param pointer_in model pointer
   * @param MinTracerIn minimum tracer pointer
   * @param phase_pair_in pair of coexisting phases
   * @param UserDefined_vwall_in is the input value for v_wall. If = -1$ then we
   * use the approximation coming from https://arxiv.org/abs/2210.16305. If =
   * -2$ then we use the upper bound from https://arxiv.org/abs/2305.02357
   * @param UserDefined_epsturb_in is the input value for epsturb. If [0..1] set
   * to value, for -1 we use the upper bound from
   * https://arxiv.org/abs/1704.05871
   * @param MaxPathIntegrations_in max number of path integrations
   * @param NumberOfInitialScanTemperatures_in number of temperature steps in
   * the initial scan of the bounce solver
   */
  BounceSolution(const std::shared_ptr<Class_Potential_Origin> &pointer_in,
                 const std::shared_ptr<MinimumTracer> &MinTracer_in,
                 const CoexPhases &phase_pair_in,
                 const double &UserDefined_vwall_in,
                 const double &UserDefined_epsturb_in,
                 const int &MaxPathIntegrations_in,
                 const size_t &NumberOfInitialScanTemperatures_in);

  /**
   * @brief Construct a new Bounce Sol Calc object. This class takes as input a
   * pair of coexisting phases.
   *
   * @param pointer_in model pointer
   * @param MinTracerIn minimum tracer pointer
   * @param phase_pair_in pair of coexisting phases
   * @param UserDefined_vwall_in is the input value for v_wall. If = -1$ then we
   * use the approximation coming from https://arxiv.org/abs/2210.16305. If =
   * -2$ then we use the upper bound from https://arxiv.org/abs/2305.02357
   * @param UserDefined_epsturb_in is the input value for epsturb. If [0..1] set
   * to value, for -1 we use the upper bound from
   * https://arxiv.org/abs/1704.05871
   * @param MaxPathIntegrations_in max number of path integrations
   * @param GroupElements_In List of allowed potential symmetries
   * @param NumberOfInitialScanTemperatures_in number of temperature steps in
   * the initial scan of the bounce solver
   */
  BounceSolution(const std::shared_ptr<Class_Potential_Origin> &pointer_in,
                 const std::shared_ptr<MinimumTracer> &MinTracer_in,
                 const CoexPhases &phase_pair_in,
                 const double &UserDefined_vwall_in,
                 const double &UserDefined_epsturb_in,
                 const int &MaxPathIntegrations_in,
                 const size_t &NumberOfInitialScanTemperatures_in,
                 std::vector<Eigen::MatrixXd> GroupElements_in);

  /**
   * @brief Initially we have no idea where the transition can occur, therefore
   * we scan the complete temperature range
   */
  void GWInitialScan();

  /**
   * @brief Calculate the euclidian action of the transition from false to true
   * phase of phase pair.
   *
   * @param T temperature
   * @param smart
   */

  void CalculateActionAt(double T, bool smart = true);

  /**
   * @brief If solution were found by the GWInitialScan() then we scan
   * temperature range in the vicinity such that we are get a enough sample to
   * then do the extrapolation.
   *
   */
  void GWSecondaryScan();

  /**
   * @brief Do linear extrapolations to calculate action at higher temperatures
   *
   */
  void GWScanTowardsHighAction();

  /**
   * @brief Do linear extrapolations to calculate action at lower temperatures
   *
   */
  void GWScanTowardsLowAction();

  /**
   * @brief Set the Bounce Sol object
   *
   */
  void SetBounceSol();

  /**
   * @brief Get the bubble wall velocity
   * @return vb
   */
  double GetWallVelocity() const;

  /**
   * @brief Get epsturb
   * @return epsturb
   */
  double GetEpsTurb() const;

  /**
   * @brief SetGstar Set gstar
   */
  void SetGstar(const double &gstar_in);

  /**
   * @brief GetGstar Get gstar
   */
  double GetGstar() const;

  /**
   * @brief SetCriticalTemp Set critical temperature
   */
  void SetCriticalTemp(const double &T_in);

  /**
   * @brief GetCriticalTemp Get critical temperature
   */
  double GetCriticalTemp() const;

  /**
   * @brief SetStoredTemp Set stored temperature
   */
  void SetStoredTemp(const double &T_in);
  /**
   * @brief GetStoredTemp Get stored temperature
   */
  double GetStoredTemp() const;

  /**
   * @brief GetNucleationTemp Get nucleation temperature via exact method
   */
  double GetNucleationTemp() const;

  /**
   * @brief GetNucleationTempApprox Get nucleation temperature via approximate
   * method
   */
  double GetNucleationTempApprox() const;

  /**
   * @brief GetPercolationTemp Get percolation temperature
   */
  double GetPercolationTemp() const;

  /**
   * @brief GetCompletionTemp Get percolation temperature
   */
  double GetCompletionTemp() const;

  /**
   * @brief CalcTransitionTemp Get transition temperature from int
   */
  double CalcTransitionTemp(const int &which_transition_temp);

  /**
   * @brief GetPTStrength Get PT strength alpha
   */
  double GetPTStrength() const;

  /**
   * @brief CalcGstarPureRad Calculate the number of effective degrees of
   * freedom assuming a purely radiative universe
   */
  void CalcGstarPureRad();

  /**
   * @brief Calculation of nucleation temperature
   */
  void CalculateNucleationTemp();

  /**
   * @brief Approximate calculation of nucleation temperature
   */
  void CalculateNucleationTempApprox();

  /**
   * @brief CalcTempAtFalseVacFraction calculates the temperature at which the
   * false vacuum fraction drops below val
   * @param false_vac_frac desired false vacuum fraction value
   * @return temperature at which false vacuum fraction drops below val
   */
  double CalcTempAtFalseVacFraction(const double &false_vac_frac);

  /**
   * @brief CalcFalseVacFraction calculates false vacuum fraction as function of
   * temperature
   * @param temp temperature
   * @return false vacuum fraction
   */
  double CalcFalseVacFraction(const double &temp);

  /**
   * @brief CalculatePercolationTemp calculation of the temperature when the
   * false vacuum fraction drops below 71 % (default)
   * @param false_vac_frac false vacuum fraction at percolation temperature, by
   * default set to 71 %
   */
  void CalculatePercolationTemp(const double &false_vac_frac = 0.71);

  /**
   * @brief CalculateCompletionTemp calculation of the temperature when the
   * false vacuum fraction drops below 1 % (default)
   * @param false_vac_frac false vacuum fraction at completion temperature, by
   * default set to 1 %
   */
  void CalculateCompletionTemp(const double &false_vac_frac = 0.01);

  /**
   * @brief Calculate phase transition strength alpha at percolation temperature
   */
  void CalculatePTStrength();

  /**
   * @brief Calculate wall velocity
   * @param false_min initial, false minimum
   * @param true_min final, true minimum
   */
  void CalculateWallVelocity(const Minimum &false_min, const Minimum &true_min);

  /**
   * @brief Calculate inverse time scale of phase transition
   */
  void CalculateInvTimeScale();

  /**
   * @brief Get inverse time scale of phase transition
   */
  double GetInvTimeScale();
};

/**
 * @brief Nintegrate_Inner Numerical integration of inner integral over inverse
 * Hubble rate for the percolation temperature calculation
 * @param obj Class reference to pass all needed parameters
 * @param T upper integration boundary
 * @return Numerical value of integral and absolute error
 */
struct resultErrorPair Nintegrate_Inner(BounceSolution &obj,
                                        const double &Tprime);

/**
 * @brief Nintegrate_Outer Numerical integration of outer integral for the
 * percolation temperature calculation
 * @param obj Class reference to pass all needed parameters
 * @return Numerical value of integral and absolute error
 */
struct resultErrorPair Nintegrate_Outer(BounceSolution &obj);

/**
 * @brief Nderive_BounceRatio Numerical derivative for the inverse time scale
 * calculation
 * @param obj Class reference to pass all needed parameters
 * @return Numerical value of derivative and absolute error
 */
struct resultErrorPair Nderive_BounceRatio(BounceSolution &obj);

} // namespace BSMPT