// Copyright (C) 2024 Lisa Biermann, Margarete Mühlleitner, Rui Santos, João
// Viana
//
// SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

/**
 * @file gravitational wave calculation
 */

#include <BSMPT/bounce_solution/bounce_solution.h> // BounceSolution
#include <BSMPT/models/SMparam.h>
#include <BSMPT/utility/Logger.h>
#include <gsl/gsl_math.h>

namespace BSMPT
{

/**
 * @brief struct to store all calculated GW data
 */
struct GravitationalWaveData
{
  bool swON             = true;  // enable sound wave contribution by default
  bool turbON           = true;  // enable turbulence contribution by default
  double transitionTemp = false; // transition temperature
  double PTStrength     = false; // strength of EW phase transition
  double InvTimeScale   = false; // inverse time scale beta/H
  double vb             = false; // bubble wall velocity

  double kappa_sw = false;   // efficiency factor for sound waves
  double K_sw     = false;   // kinetic energy fraction in fluid of total bubble
                             // energy for sound waves
  double HR         = false; // time scale times max. velocity for sound waves
  double kappa_turb = false; // efficiency factor for turbulence
  double K_turb     = false; // kinetic energy fraction in fluid of total bubble
                             // energy for turbulence
  double gstar = false;      // number of eff. d.o.f.
  double Hstar = false;      // Hubble rate at percolation temperature

  double fPeakSoundWave        = false; // peak frequency for sound wave
  double h2OmegaPeakSoundWave  = false; // peak amplitude for sound wave
  double fPeakTurbulence       = false; // peak frequency for turbulence
  double h2OmegaPeakTurbulence = false; // peak amplitude for turbulence

  StatusGW status = StatusGW::NotSet; // gw calculation status
};

class GravitationalWave
{
private:
public:
  GravitationalWaveData data;
  GravitationalWave(BounceSolution &BACalc,
                    const int &which_transition_temp = 3);
  ~GravitationalWave();

  /**
   * @brief CalcEpsTurb calculate epsilon for turbulence contribution
   * @param epsturb_in is the input value for epsturb. If [0..1] set
   * to value, for -1 we use the upper bound from
   * https://arxiv.org/abs/1704.05871
   * @return value for epsturb
   */
  double CalcEpsTurb(double epsturb_in);

  /**
   * @brief AbsErr absolute error for numerical integration
   */
  const double AbsErr = 0;

  /**
   * @brief RelErr relative error for numerical integration
   */
  const double RelErr = 1e-6;

  /**
   * @brief reduced Hubble constant
   *
   */
  double h = 0.674;

  /**
   * @brief Calculate peak frequency of GW signal for sound waves
   */
  void CalcPeakFrequencySoundWave();

  /**
   * @brief Calculate peak amplitude of GW signal for sound waves
   */
  void CalcPeakAmplitudeSoundWave();

  /**
   * @brief Calculate peak frequency of GW signal from turbulence
   */
  void CalcPeakFrequencyTurbulence();

  /**
   * @brief Calculate peak amplitude of GW signal from turbulence
   */
  void CalcPeakAmplitudeTurbulence();

  /**
   * @brief Amplitude of GW signal as a function of
   * @param f frequency
   * @param swON true = contribution from sound waves switched on,
   * false = switched off
   * @param turbON true = contribution from turbulence switched on,
   * false = switched off
   * @return h2OmegaGW
   */
  double CalcGWAmplitude(double f, bool swON, bool turbON);

  /**
   * @brief GetSNR
   * @param fmin minimal frequency
   * @param fmax maximal frequency
   * @param T duration of exp. data acquisition, default value: 3 years
   * @return signal-to-noise (SNR) ratio at LISA
   */
  double GetSNR(const double fmin, const double fmax, const double T = 3);

  /**
   * @brief snr_integrand friend to define inner integrand of SNR integral
   */
  friend double snr_integrand(double freq, void *params);
};

/**
 * @brief SIfunc
 * @param f frequency
 * @return value of SI-function from LISA mission performance requirement for
 * power spectral density
 */
double SIfunc(const double f);

/**
 * @brief Rfunc
 * @param f frequency
 * @return value of R-function from LISA mission performance requirement for
 * power spectral density
 */
double Rfunc(const double f);

/**
 * @brief powspec_density
 * @param f frequency
 * @return value of power spectral density-function from LISA mission
 * performance requirement
 */
double powspec_density(const double f);

/**
 * @brief return the value of LISA mission nominal sensitivity
 *
 * @param f frequency
 * @return value of LISA's nominal sensitivity for a given frequency
 */
double h2OmSens(const double f);

/**
 * @brief Nintegrate_SNR Numerical integration of SNR integral
 * @param obj Class reference to pass all needed parameters
 * @return Numerical value of integral and absolute error
 */
struct resultErrorPair
Nintegrate_SNR(GravitationalWave &obj, const double fmin, const double fmax);

/**
 * @brief Get efficiency factor kappa_sw for sound waves
 *
 * @param alpha strength of the phase transition
 * @param vwall bubble wall velocity
 * @param Csound speed of sound
 * @return double
 */
double
Getkappa_sw(const double &alpha, const double &vwall, const double &Csound);

/**
 * @brief Get K for sound waves
 *
 * @param alpha strength of the phase transition
 * @param vwall bubble wall velocity
 * @param Csound speed of sound
 * @return double
 */
double GetK_sw(const double &alpha, const double &vwall, const double &Csound);

/**
 * @brief Get HR for sound waves
 *
 * @param invTimeScale beta/H, inverse time scale
 * @param vwall bubble wall velocity
 * @param Csound speed of sound
 * @return double
 */
double
GetHR(const double &invTimeScale, const double &vwall, const double &Csound);

/**
 * @brief Get K for turbulence
 *
 * @param alpha strength of the phase transition
 * @param kappa turbulence efficiency factor
 * @return double
 */
double GetK_turb(const double &alpha, const double &kappa);

/**
 * @brief Determine fluid turnover time regime
 *
 * @param HR mean bubble separation
 * @param K fraction of the kinetic energy in the fluid
 * @return true if H*tauSH approx. 1, false if smaller than 1
 */
bool IsFluidTurnoverApproxOne(const double &HR, const double &K);

} // namespace BSMPT