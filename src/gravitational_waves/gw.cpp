// SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file gravitational wave class
 */

#include <BSMPT/gravitational_waves/gw.h>

namespace BSMPT
{

GravitationalWave::GravitationalWave(BounceSolution &BACalc,
                                     const int &which_transition_temp)
{
  data.transitionTemp = BACalc.CalcTransitionTemp(which_transition_temp);
  data.PTStrength     = BACalc.GetPTStrength();
  data.betaH          = BACalc.GetInvTimeScale();
  data.vw             = BACalc.GetWallVelocity();
  data.kappa          = Getkappa(
      data.PTStrength, data.vw, BACalc.modelPointer->SMConstants.Csound);
  data.K =
      GetK(data.PTStrength, data.vw, BACalc.modelPointer->SMConstants.Csound);
  data.HR = GetHR(data.betaH, data.vw, BACalc.modelPointer->SMConstants.Csound);
  data.Hstar0 = GetHstar0(BACalc.CalcTransitionTemp(which_transition_temp),
                          BACalc.GetGstar());
  data.Epsilon_Turb = BACalc.GetEpsTurb();
  data.kappa_turb   = CalcEpsTurb(BACalc.GetEpsTurb()) * data.kappa;
  data.gstar        = BACalc.GetGstar();
  data.Hstar        = BACalc.HubbleRate(data.transitionTemp);
  data.FGW0         = 1.64 * 1.e-5 * pow(100. / BACalc.GetGstar(), 1 / 3.);
  data.Csound       = BACalc.modelPointer->SMConstants.Csound;

  if (data.betaH < 1)
  {
    data.status = StatusGW::Failure;
    Logger::Write(
        LoggingLevel::GWDetailed,
        "beta/H < 1 detected, with beta/H = " + std::to_string(data.betaH) +
            ". Transition is assumed to happen, but no GWs calculated.");
  }
  else if (data.transitionTemp == -1)
  {
    data.status = StatusGW::Failure;
    Logger::Write(LoggingLevel::GWDetailed,
                  "Requested transition temperature could not be calculated.");
  }
}

GravitationalWave::~GravitationalWave()
{
}

double GravitationalWave::CalcEpsTurb(double epsturb_in)
{
  if (epsturb_in == -1)
  {
    double HtauSW = 2. / std::sqrt(3) * data.HR / std::sqrt(data.K);
    return std::pow((1 - std::min(HtauSW, 1.)), 2. / 3.);
  }
  else
  {
    return epsturb_in;
  }
}

void GravitationalWave::CalcPeakCollision()
{
  const double n1 = 2.4;
  const double n2 = -2.4;
  const double a1 = 1.2;

  data.CollisionParameter.n1 = n1;
  data.CollisionParameter.n2 = n2;
  data.CollisionParameter.a1 = a1;

  // Calculate characteristic frequency

  const double f_p            = 0.11 * data.Hstar0 / data.betaH;
  data.CollisionParameter.f_b = f_p * pow(-n1 / n2, -1 / a1);

  // Calculate amplitude for collisions
  const double A_str  = 0.05;
  const double Ktilde = GetKtilde(data.PTStrength);
  const double Omega_p =
      data.FGW0 * A_str * pow(Ktilde, 2) / pow(data.betaH, 2);

  data.CollisionParameter.Omega_b =
      Omega_p * pow(1. / 2. * pow(-n2 / n1, n1 / (n1 - n2)) +
                        1. / 2. * pow(-n1 / n2, -n2 / (n1 - n2)),
                    (n1 - n2) / a1);
}

void GravitationalWave::CalcPeakSoundWave()
{
  data.SoundWaveParameter.n1 = 3.;
  data.SoundWaveParameter.n2 = 1.;
  data.SoundWaveParameter.n3 = -3.;
  data.SoundWaveParameter.a1 = 2.;
  data.SoundWaveParameter.a2 = 4.;

  // Characteristic frequencies

  const double xi_shell = abs(data.vw - data.Csound);
  const double delta_w  = xi_shell / max(data.vw, data.Csound);

  const double f_1 = 0.2 * data.Hstar0 / data.HR;
  const double f_2 = 0.5 * data.Hstar0 / (data.HR * delta_w);

  // Sound wave amplitude

  const double Asw = 0.11; // Numerical simulation
  const double HtauSW =
      std::min(1., 2. / std::sqrt(3) * data.HR / std::sqrt(data.K));

  const double Omega_int = data.FGW0 * Asw * pow(data.K, 2) * HtauSW * data.HR;

  // Convert Omega_int to Omega_2

  data.SoundWaveParameter.f_1 = f_1;
  data.SoundWaveParameter.f_2 = f_2;
  data.SoundWaveParameter.Omega_2 =
      Omega_int * (sqrt(2) + (2 * f_2 / f_1) / (1 + pow(f_2 / f_1, 2))) / M_PI;
}

void GravitationalWave::CalcPeakTurbulence()
{
  data.TurbulanceParameter.n1 = 3.;
  data.TurbulanceParameter.n2 = 1.;
  data.TurbulanceParameter.n3 = -8. / 3.;
  data.TurbulanceParameter.a1 = 4.;
  data.TurbulanceParameter.a2 = 2.15;

  // Characteristic frequencies
  const double A       = 0.085;
  const double N       = 2.;
  const double Omega_s = data.Epsilon_Turb * data.K;

  data.TurbulanceParameter.f_1 =
      sqrt(3 * Omega_s * data.Hstar0) / (2. * N * data.HR);
  data.TurbulanceParameter.f_2 = 2.2 * data.Hstar0 / (data.HR);

  // Calculate Omega_2 for turbulence

  const double A_MHD = 3 * 2.2 * A / (4 * pow(M_PI, 2)) *
                       pow(2, -11 / (3 * data.TurbulanceParameter.a2.value()));

  data.TurbulanceParameter.Omega_2 =
      data.FGW0 * A_MHD * pow(Omega_s, 2) * pow(data.HR, 2);
}

double GravitationalWave::BPL(const double &f, const BPLParameters &par) const
{
  const double Omega_b = par.Omega_b.value();
  const double f_b     = par.f_b.value();
  const double n1      = par.n1.value();
  const double n2      = par.n2.value();
  const double a1      = par.a1.value();

  return Omega_b * pow(f / f_b, n1) *
         pow(0.5 + 0.5 * pow(f / f_b, a1), (n2 - n1) / a1);
}

double GravitationalWave::DBPL(const double &f, const DBPLParameters &par) const
{
  const double Omega_2 = par.Omega_2.value();
  const double f_1     = par.f_1.value();
  const double f_2     = par.f_2.value();
  const double n1      = par.n1.value();
  const double n2      = par.n2.value();
  const double n3      = par.n3.value();
  const double a1      = par.a1.value();
  const double a2      = par.a2.value();

  const double Sf = pow(f / f_1, n1) *
                    pow(1 + pow(f / f_1, a1), (-n1 + n2) / a1) *
                    pow(1 + pow(f / f_2, a2), (-n2 + n3) / a2);
  const double S2 = pow(f_2 / f_1, n1) *
                    pow(1 + pow(f_2 / f_1, a1), (-n1 + n2) / a1) *
                    pow(1 + pow(f_2 / f_2, a2), (-n2 + n3) / a2);
  return Omega_2 * Sf / S2;
}

double GravitationalWave::CalcGWAmplitude(double f)
{
  double res = 0;
  if (data.swON)
  {
    if (!data.SoundWaveParameter.IsDefined()) this->CalcPeakSoundWave();
    res += DBPL(f, data.SoundWaveParameter);
  }
  if (data.turbON)
  {
    if (!data.TurbulanceParameter.IsDefined()) this->CalcPeakTurbulence();
    res += DBPL(f, data.TurbulanceParameter);
  }
  if (data.collisionON)
  {
    if (!data.CollisionParameter.IsDefined()) this->CalcPeakTurbulence();
    res += BPL(f, data.CollisionParameter);
  }
  return h * h * res; // Reduced hubble factor \f$ h^2 \f$
}

double
GravitationalWave::GetSNR(const double fmin, const double fmax, const double T)
{
  auto integral     = Nintegrate_SNR(*this, fmin, fmax);
  double res        = std::sqrt(86400 * 365.25 * T * integral.result);
  this->data.status = StatusGW::Success;
  return res;
}

double SIfunc(const double f)
{
  double f1 = 0.4e-3;
  return 5.76e-48 * (1 + std::pow(f1 / f, 2));
}

double Rfunc(const double f)
{
  double f2 = 25e-3;
  return 1. + std::pow(f / f2, 2);
}

double powspec_density(const double f)
{
  double SIIfunc = 3.6e-41;
  return 1. / 2 * 20. / 3 * (SIfunc(f) / std::pow(2 * M_PI * f, 4) + SIIfunc) *
         Rfunc(f);
}

double h2OmSens(const double f)
{
  double H0 = 100 / 3.09e19;
  return (4 * std::pow(M_PI, 2)) / (3 * std::pow(H0, 2)) * std::pow(f, 3) *
         powspec_density(f);
}

double snr_integrand(double f, void *params)
{
  class GravitationalWave &obj = *static_cast<GravitationalWave *>(params);

  double func = std::pow(obj.CalcGWAmplitude(f) / (h2OmSens(f)), 2);
  return func;
}

struct resultErrorPair
Nintegrate_SNR(GravitationalWave &obj, const double fmin, const double fmax)
{
  double abs_err = obj.AbsErr;
  double rel_err = obj.RelErr;

  std::size_t workspace_size = 1000;
  gsl_integration_workspace *w =
      gsl_integration_workspace_alloc(workspace_size);
  gsl_function F;
  F.function = &snr_integrand;
  F.params   = static_cast<void *>(&obj);

  struct resultErrorPair res;

  gsl_integration_qags(&F,
                       fmin,
                       fmax,
                       abs_err,
                       rel_err,
                       workspace_size,
                       w,
                       &res.result,
                       &res.error);

  gsl_integration_workspace_free(w);

  return res;
}

double Getkappa(const double &alpha, const double &vwall, const double &Csound)
{
  double kappa;
  double kappaA = std::pow(vwall, 6.0 / 5.0) * 6.9 * alpha /
                  (1.36 - 0.037 * std::sqrt(alpha) + alpha);
  double kappaB =
      std::pow(alpha, 2.0 / 5.0) / (0.017 + std::pow(0.997 + alpha, 2.0 / 5.0));
  double kappaC = std::sqrt(alpha) / (0.135 + std::sqrt(0.98 + alpha));
  double kappaD = alpha / (0.73 + 0.083 * std::sqrt(alpha) + alpha);
  double xiJ =
      (sqrt((2.0 / 3.0) * alpha + alpha * alpha) + std::sqrt(1.0 / 3.0)) /
      (1 + alpha);
  double deltaK = -0.9 * log((sqrt(alpha) / (1 + std::sqrt(alpha))));

  if (vwall < Csound)
    kappa = std::pow(Csound, 11.0 / 5.0) * kappaA * kappaB /
            ((pow(Csound, 11.0 / 5.0) - std::pow(vwall, 11.0 / 5.0)) * kappaB +
             vwall * std::pow(Csound, 6.0 / 5.0) * kappaA);
  else if (vwall > xiJ)
    kappa = std::pow(xiJ - 1, 3.0) * std::pow(xiJ, 5.0 / 2.0) *
            std::pow(vwall, -5.0 / 2.0) * kappaC * kappaD /
            ((pow(xiJ - 1, 3.0) - std::pow(vwall - 1, 3.0)) *
                 std::pow(xiJ, 5.0 / 2.0) * kappaC +
             std::pow(vwall - 1, 3.0) * kappaD);
  else
    kappa = kappaB + (vwall - Csound) * deltaK +
            (pow(vwall - Csound, 3.0) / std::pow(xiJ - Csound, 3.0)) *
                (kappaC - kappaB - (xiJ - Csound) * deltaK);

  return kappa;
}

double GetK(const double &alpha, const double &vwall, const double &Csound)
{
  double kappa = Getkappa(alpha, vwall, Csound);
  return 0.6 * kappa * alpha / (1. + alpha);
}

double GetHR(const double &betaH, const double &vwall, const double &Csound)
{
  double max_velo = std::max(vwall, Csound);
  return 1. / betaH * std::pow(8 * M_PI, 1. / 3) * max_velo;
}

double GetHstar0(const double &temp, const double &gstar)
{
  return 1.65 * 1e-5 * pow(gstar / 100, 1 / 6.) * (temp / 100);
}

double GetKtilde(const double &alpha)
{
  return alpha / (1. + alpha);
}

} // namespace BSMPT