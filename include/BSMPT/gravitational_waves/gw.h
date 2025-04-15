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
  double reheatingTemp  = false; // reheating temperature
  double PTStrength     = false; // strength of EW phase transition
  double betaH          = false; // inverse time scale beta/H
  double vw             = false; // bubble wall velocity
  double vCJ            = false; // Chapman-Jouguet velocity
  double XiShock        = false; // Shock speed
  double Csound_false   = false; // speed sound false vacuum
  double Csound_true    = false; // speed sound true vacuum
  double kappa_col      = false; // efficiency factor for collision
  double kappa_sw       = false; // efficiency factor for sw
  double K_sw           = false; // kinetic energy fraction for sound waves
  double HR             = false; // time scale x max. velocity for sound waves
  int pnlo_scaling      = false; // pressure scaling at NLO, 1 -> N processes
  double Epsilon_Turb   = false; //  fraction of overall kinetic energy in bulk
                                 //  motion that is converted to MHD
  double gstar = false;          // number of eff. d.o.f.
  double FGW0  = false; // Redshift factor for the fractional energy density
  double Hstar0 =
      false; // Reduced Hubble rate at transition temperature redshift today

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
                    const TransitionTemperature &which_transition_temp =
                        TransitionTemperature::Percolation);
  ~GravitationalWave();

  /**
   * @brief CalcEpsTurb calculate epsilon for turbulence contribution
   * @param epsturb_in is the input value for epsturb. If [0..1] set
   * to value, for -1 we use the upper bound sqrt(1 - Upsilon)
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
   * @brief Calcualte the fluid shell thickness \f$ \xi_\text{front} -
   * \xi_\text{rear} \f$
   *
   * @return double
   */
  double CalculateXiShell();

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
 * @param kappa_sw efficiency factor
 * @return double
 */
double GetK_sw(const double &alpha, const double &kappa_sw);

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

/**
 * @brief Calculate \f$ H_*\tau_\text{sh} = H_* R_* / \sqrt{\bar v_f^2} \f$
 *
 * @param HR
 * @param K_sw
 * @return double
 */
double GetHtauSH(const double HR, const double K_sw);

/**
 * @brief Calculate \f$ H_*\tau_\text{sw} = \min(1,H_*\tau_\text{sh}) \f$
 *
 * @param HR
 * @param K_sw
 * @return double
 */
double GetHtauSW(const double HR, const double K_sw);

/**
 * @brief Calculate \f$
 * \Upsilon=1-\frac{1}{\sqrt{1+2H_{\ast}\tau_{\mathrm{sw}}}} \f$ from
 * https://arxiv.org/abs/1903.09642
 * @param HR
 * @param K_sw
 * @return double
 */
double GetYpsilon(const double HR, const double K_sw);

namespace kappa
{
// Compute kappa_sw https://arxiv.org/abs/2010.09744

/**
 * @brief Expansion modes of the fluid shell
 *
 */
enum class ExpansionMode
{
  Deflagration,
  Hybrid,
  Detonation
};

/**
 * @brief Lorentz boost
 *
 * @param a
 * @param b
 * @return double
 */
double mu(double a, double b);
/**
 * @brief encode the (special relativistic) relative velocity and
the ratio of the enthalpies across the bubble wall.
 *
 * @param a
 * @param b
 * @return double
 */
double getwow(double a, double b);
void custom_error_handler(const char *reason,
                          const char *file,
                          int line,
                          int gsl_errno);
/**
 * @brief returns the fluid velocity behind the wall, v−, and the expansion mode
 *
 * @param al \f$ \alpha \f$
 * @param vw \f$ v_w \f$
 * @param cs2b \f$ c_{s,b}^2 $ sound speed squared in the broken phase
 * @return std::pair<double, ExpansionMode>
 */
std::pair<double, ExpansionMode> getvm(double al, double vw, double cs2b);
/**
 * @brief encodes the differential equation solved in the shock/rarefaction wave
 * and returns (dξ/dv, dw/dv).
 *
 * @param v
 * @param y
 * @param dydv
 * @param params
 * @return int
 */
int dfdv(double v, const double y[], double dydv[], void *params);
/**
 * @brief solve the hydrodynamics equation
 *
 * @param vw
 * @param v0
 * @param cs2
 * @return std::vector<std::vector<double>>
 */
std::vector<std::vector<double>> solve_ode(double vw, double v0, double cs2);
double integrate(const std::vector<double> &y, const std::vector<double> &x);
/**
 * @brief Get the Kand Wow object returns the enthalpy-weighted kinetic energy
 * in the shock/rarefaction wave and the ratio between the enthalpy density at
 * the start of the shock/rarefaction compared to its end (for the shocks, the
 * end is in the phase in front of the shock; for the rarefaction wave, the
 * enthalpy density is normalized to 1 behind the wall and has to be rescaled in
 * the other part of the code).
 *
 * @param vw
 * @param v0
 * @param cs2
 * @param vprofile velocity profile
 * @return std::pair<double, double>
 */
std::pair<double, double>
getKandWow(double vw,
           double v0,
           double cs2,
           std::vector<std::vector<double>> &vprofile);
/**
 * @brief returns α¯θn in the nucleation phase (in front of the shock) for a
 * given α¯θ+ value at the wall.
 *
 * @param al
 * @param wow
 * @param cs2b
 * @param cs2s
 * @return double
 */
double alN(double al, double wow, double cs2b, double cs2s);
/**
 * @brief Calculate the \f$ \kappa_{sw} \f$
 *
 * @param cs2b sound speed in the true vacuum \f$ c_{s,b}^2 =
 * \frac{1}{T}\frac{\frac{dV(\phi_t)}{dT}}{\frac{d^2V(\phi_t)}{dT^2}} \f$
 * @param cs2s sound speed in the false vacuum \f$ c_{s,b}^2 =
 * \frac{1}{T}\frac{\frac{dV(\phi_f)}{dT}}{\frac{d^2V(\phi_f)}{dT^2}} \f$
 * @param al \f$ \alpha \f$
 * @param vw \f$ v_w \f$
 * @param vprofile velocity profile
 * @return double effiency factor for sound waves
 */
double kappaNuMuModel(double cs2b,
                      double cs2s,
                      double al,
                      double vw,
                      std::vector<std::vector<double>> &vprofile);
/**
 * @brief Calculate the \f$ \kappa_{sw} \f$
 *
 * @param cs2b sound speed in the true vacuum \f$ c_{s,b}^2 =
 * \frac{1}{T}\frac{\frac{dV(\phi_t)}{dT}}{\frac{d^2V(\phi_t)}{dT^2}} \f$
 * @param cs2s sound speed in the false vacuum \f$ c_{s,b}^2 =
 * \frac{1}{T}\frac{\frac{dV(\phi_f)}{dT}}{\frac{d^2V(\phi_f)}{dT^2}} \f$
 * @param al \f$ \alpha \f$
 * @param vw \f$ v_w \f$
 * @return double effiency factor for sound waves
 */
double kappaNuMuModel(double cs2b, double cs2s, double al, double vw);

/**
 * @brief Calculate \f$ \kappa_{col} \f$ the efficiency factor for collisions
 * https://arxiv.org/pdf/1903.09642
 * https://arxiv.org/pdf/2208.11697
 * https://arxiv.org/pdf/2210.07075
 *
 * @param Tstar transition temperature
 * @param pnlo_scaling scaling of pressure at NLO, 1 -> N processes
 * @param HR Hubble rate x Mean bubble size
 * @param BACalc BACalc object
 * @return double return the collision effiency factor
 */
double Getkappa_col(const double &Tstar,
                    const int &pnlo_scaling,
                    const double &HR,
                    BounceSolution &BACalc);
} // namespace kappa
} // namespace BSMPT