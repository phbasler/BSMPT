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
  data.InvTimeScale   = BACalc.GetInvTimeScale();
  data.vb             = BACalc.GetWallVelocity();
  data.kappa_sw       = Getkappa_sw(
      data.PTStrength, data.vb, BACalc.modelPointer->SMConstants.Csound);
  data.K_sw = GetK_sw(
      data.PTStrength, data.vb, BACalc.modelPointer->SMConstants.Csound);
  data.HR = GetHR(
      data.InvTimeScale, data.vb, BACalc.modelPointer->SMConstants.Csound);
  data.kappa_turb = CalcEpsTurb(BACalc.GetEpsTurb()) * data.kappa_sw;
  data.K_turb     = GetK_turb(data.PTStrength, data.kappa_turb);
  data.gstar      = BACalc.GetGstar();
  data.Hstar      = BACalc.HubbleRate(data.transitionTemp);

  if (data.InvTimeScale < 1)
  {
    data.status = StatusGW::Failure;
    Logger::Write(
        LoggingLevel::GWDetailed,
        "beta/H < 1 detected, with beta/H = " +
            std::to_string(data.InvTimeScale) +
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
    double HtauSW = 2. / std::sqrt(3) * data.HR / std::sqrt(data.K_sw);
    return std::pow((1 - std::min(HtauSW, 1.)), 2. / 3.);
  }
  else
  {
    return epsturb_in;
  }
}

void GravitationalWave::CalcPeakFrequencySoundWave()
{
  double res = 26e-6 * (1. / this->data.HR) *
               (this->data.transitionTemp / 100) *
               std::pow(this->data.gstar / 100, 1. / 6);
  this->data.fPeakSoundWave = res;
}

void GravitationalWave::CalcPeakAmplitudeSoundWave()
{
  double ratio = 0;
  if (IsFluidTurnoverApproxOne(this->data.HR, this->data.K_sw))
  {
    ratio = this->data.HR * std::pow(this->data.K_sw, 2.);
  }
  else
  {
    ratio = 2. / std::sqrt(3) * std::pow(this->data.HR, 2) *
            std::pow(this->data.K_sw, 3. / 2);
  }

  //  Taken from erratum of https://arxiv.org/pdf/1704.05871.pdf
  this->data.h2OmegaPeakSoundWave = h * h * 2.061 * 1.2e-2 * 3.57e-5 *
                                    std::pow(100. / this->data.gstar, 1. / 3.) *
                                    ratio;
}

void GravitationalWave::CalcPeakFrequencyTurbulence()
{
  double res = 7.909e-5 * (1. / this->data.HR) *
               (this->data.transitionTemp / 100) *
               std::pow(this->data.gstar / 100., 1. / 6.);
  this->data.fPeakTurbulence = res;
}

void GravitationalWave::CalcPeakAmplitudeTurbulence()
{
  double res = 1.144e-4 * std::pow(100. / this->data.gstar, 1. / 3.) *
               this->data.HR * std::pow(this->data.K_turb, 3. / 2.);
  this->data.h2OmegaPeakTurbulence = res;
}

double GravitationalWave::CalcGWAmplitude(double f, bool swON, bool turbON)
{
  double res = 0;
  if (swON)
  {
    if (!this->data.h2OmegaPeakSoundWave)
    {
      this->CalcPeakAmplitudeSoundWave();
    }
    else if (!this->data.fPeakSoundWave)
    {
      this->CalcPeakFrequencySoundWave();
    }

    res += this->data.h2OmegaPeakSoundWave * std::pow(4. / 7, -7. / 2) *
           std::pow(f / this->data.fPeakSoundWave, 3) *
           std::pow((1 + 3. / 4 * std::pow(f / this->data.fPeakSoundWave, 2)),
                    -7. / 2);
  }
  if (turbON)
  {
    if (!this->data.h2OmegaPeakTurbulence)
    {
      this->CalcPeakAmplitudeTurbulence();
    }
    else if (!this->data.fPeakTurbulence)
    {
      this->CalcPeakFrequencyTurbulence();
    }

    res += this->data.h2OmegaPeakTurbulence *
           std::pow(f / this->data.fPeakTurbulence, 3) /
           std::pow(1 + f / this->data.fPeakTurbulence, 11 / 3) /
           (1 + 8 * M_PI * f / this->data.Hstar);
  }

  return res;
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
  return (2 * std::pow(M_PI, 2)) / (3 * std::pow(H0, 2)) * std::pow(f, 3) *
         powspec_density(f);
}

double snr_integrand(double f, void *params)
{
  class GravitationalWave &obj = *static_cast<GravitationalWave *>(params);

  double func = std::pow(
      obj.CalcGWAmplitude(f, obj.data.swON, obj.data.turbON) / (h2OmSens(f)),
      2);
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

double
Getkappa_sw(const double &alpha, const double &vwall, const double &Csound)
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

double GetK_sw(const double &alpha, const double &vwall, const double &Csound)
{
  double kappa = Getkappa_sw(alpha, vwall, Csound);
  return kappa * alpha / (1. + alpha);
}

double
GetHR(const double &invTimeScale, const double &vwall, const double &Csound)
{
  double max_velo = std::max(vwall, Csound);
  return 1. / invTimeScale * std::pow(8 * M_PI, 1. / 3) * max_velo;
}

double GetK_turb(const double &alpha, const double &kappa)
{
  return kappa * alpha / (1. + alpha);
}

bool IsFluidTurnoverApproxOne(const double &HR, const double &K)
{
  double ratio = 2 * HR / std::sqrt(3 * K);
  return (ratio < 1) ? false : true;
}

} // namespace BSMPT