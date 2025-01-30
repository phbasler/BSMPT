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
#include <algorithm>
#include <cmath>
#include <functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <numeric>
#include <vector>
namespace BSMPT
{

struct BPLParameters
{
  std::optional<double> Omega_b;
  std::optional<double> f_b;
  std::optional<double> n1;
  std::optional<double> n2;
  std::optional<double> a1;

  inline bool IsDefined() const
  {
    return Omega_b.has_value() && f_b.has_value() && n1.has_value() &&
           n2.has_value() && a1.has_value();
  }
};

struct DBPLParameters
{
  std::optional<double> Omega_2;
  std::optional<double> f_1;
  std::optional<double> f_2;
  std::optional<double> n1;
  std::optional<double> n2;
  std::optional<double> n3;
  std::optional<double> a1;
  std::optional<double> a2;

  inline bool IsDefined() const
  {
    return Omega_2.has_value() && f_1.has_value() && f_2.has_value() &&
           n1.has_value() && n2.has_value() && n3.has_value() &&
           a1.has_value() && a2.has_value();
  }
};

/**
 * @brief struct to store all calculated GW data
 */
struct GravitationalWaveData
{
  bool swON             = true;  // enable sound wave contribution by default
  bool turbON           = true;  // enable turbulence contribution by default
  bool collisionON      = true;  // enable collision contribution by default
  double transitionTemp = false; // transition temperature
  double PTStrength     = false; // strength of EW phase transition
  double betaH          = false; // inverse time scale beta/H
  double vw             = false; // bubble wall velocity
  double Csound_false   = false; // speed sound false vacuum
  double Csound_true    = false; // speed sound true vacuum
  double kappa = false; // kinetic energy fraction of a single expanding bubble
  double K     = false; // kinetic energy fraction
  double HR    = false; // time scale times max. velocity for sound waves
  double Epsilon_Turb = false; //  fraction of overall kinetic energy in bulk
                               //  motion that is converted to MHD
  double gstar = false;        // number of eff. d.o.f.
  double Hstar = false;        // Hubble rate at transition temperature
  double Hstar0 =
      false; // Reduced Hubble rate at transition temperature redshift today
  double FGW0 = false; // Redshift factor for the fractional energy density

  BPLParameters CollisionParameter;
  DBPLParameters SoundWaveParameter;
  DBPLParameters TurbulanceParameter;

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
   * @brief Calculate peak amplitude and frequency for GW signal from collision
   */
  void CalcPeakCollision();

  /**
   * @brief Calculate peak amplitude and frequency for GW signal from sound
   * waves
   */
  void CalcPeakSoundWave();

  /**
   * @brief Calculate peak amplitude and frequency for GW signal from turbulence
   */
  void CalcPeakTurbulence();

  /**
   * @brief Calculate sound speed \f$ c_s^2 = \frac{1}{T}
   * \frac{V'_{\text{eff}}(T)}{V''_{\text{eff}}(T)} \f$
   *
   * @param Tstar transition temperature.
   * @param phase phase where we calculate the sound speed.
   * @param modelPointer model pointer.
   * @return double sound speed in that phase at that temperature.
   */
  double CalculateSoundSpeed(
      const double &Tstar,
      Phase &phase,
      const std::shared_ptr<Class_Potential_Origin> &modelPointer);

  /**
   * @brief Broken power law spectrum \f$
   * \Omega_{\mathrm{GW}}^{\mathrm{BPL}}\left(f,
   * \vec{\theta}_{\mathrm{Cosmo}}\right)=\Omega_b\left(\frac{f}{f_b}\right)^{n_1}\left[\frac{1}{2}+\frac{1}{2}\left(\frac{f}{f_b}\right)^{a_1}\right]^{\frac{n_2-n_1}{a_1}}
   * \f$
   *
   * @param f frequency
   * @param par spectrum parameters
   * @return double amplitude at that frequency
   */
  double BPL(const double &f, const BPLParameters &par) const;

  /**
   * @brief Double broken power law spectrum \f$
   * 2^\frac{n_2-n_3}{a_2}\left(1+\left(\frac{f_2}{f_1}\right)^{a_1}\right)^{\frac{n_1-n_2}{a_1}}\left(\frac{f}{f_2}\right)^{n_1}\left(1+\left(\frac{f}{f_1}\right)^{a_1}\right)^{\frac{n_2-n_1}{a_1}}\left(1+\left(\frac{f}{f_2}\right)^{a_2}\right)^{\frac{n_3-n_2}{a_2}}
   * \f$
   *
   * @param f frequency
   * @param par spectrum parameters
   * @return double amplitude at that frequency
   */
  double DBPL(const double &f, const DBPLParameters &par) const;

  /**
   * @brief Amplitude of GW signal as a function of
   * @param f frequency
   * @return h2OmegaGW
   */
  double CalcGWAmplitude(double f);

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
 * @brief Get the kinetic energy fraction \f$ K \f$
 *
 * @param alpha strength of the phase transition
 * @param kappa efficiency factor
 * @return double
 */
double GetK(const double &alpha, const double &kappa);

/**
 * @brief Get HR
 *
 * @param betaH beta/H, inverse time scale
 * @param vwall bubble wall velocity
 * @param Csound speed of sound
 * @return double
 */
double GetHR(const double &betaH, const double &vwall, const double &Csound);

/**
 * @brief Calculate the Hubble rate at transition time refshifted to today
 *
 * @param temp transition temperature
 * @param gstar effective d.o.f. at transition time
 */
double GetHstar0(const double &temp, const double &gstar);

/**
 * @brief Get \f$ \tilde{K} = \frac{\alpha}{1+\alpha}\f$
 *
 * @param alpha strength of the phase transition
 * @return double
 */
double GetKtilde(const double &alpha);

namespace kappa
{
// Compute kappa https://arxiv.org/abs/2010.09744
double mu(double a, double b);
double getwow(double a, double b);
void custom_error_handler(const char *reason,
                          const char *file,
                          int line,
                          int gsl_errno);
std::pair<double, int> getvm(double al, double vw, double cs2b);
int dfdv(double v, const double y[], double dydv[], void *params);
std::vector<std::vector<double>> solve_ode(double vw, double v0, double cs2);
double integrate(const std::vector<double> &y, const std::vector<double> &x);
std::pair<double, double> getKandWow(double vw, double v0, double cs2);
double alN(double al, double wow, double cs2b, double cs2s);
double kappaNuMuModel(double cs2b, double cs2s, double al, double vw);
} // namespace kappa
} // namespace BSMPT