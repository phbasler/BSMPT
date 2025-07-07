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

GravitationalWave::GravitationalWave(
    BounceSolution &BACalc,
    const TransitionTemperature &which_transition_temp)
{
  BACalc.SetAndCalculateGWParameters(which_transition_temp);
  data.transitionTemp = BACalc.GetTransitionTemp();
  data.reheatingTemp  = BACalc.GetReheatingTemp();
  data.PTStrength     = BACalc.GetPTStrength();
  data.betaH          = BACalc.GetInvTimeScale();
  data.vw             = BACalc.GetWallVelocity();
  Logger::Write(LoggingLevel::GWDetailed,
                "\n--------------- Transition parameters ---------------\n");
  Logger::Write(LoggingLevel::GWDetailed,
                "T* = " + std::to_string(data.transitionTemp));

  Logger::Write(LoggingLevel::GWDetailed,
                "Th = " + std::to_string(data.reheatingTemp));
  Logger::Write(LoggingLevel::GWDetailed,
                "alpha = " + std::to_string(data.PTStrength));
  Logger::Write(LoggingLevel::GWDetailed,
                "beta/H = " + std::to_string(data.betaH));
  Logger::Write(LoggingLevel::GWDetailed,
                "vw = " + std::to_string(data.vw) + "\n");
  // Sound speeds
  data.Csound_false = BACalc.GetSoundSpeedFalse();
  data.Csound_true  = BACalc.GetSoundSpeedTrue();

  Logger::Write(LoggingLevel::GWDetailed,
                "Csound_false = " + std::to_string(data.Csound_false));
  Logger::Write(LoggingLevel::GWDetailed,
                "Csound_true = " + std::to_string(data.Csound_true));
  // General purpose GW parameters
  data.HR    = BACalc.GetRstar() * BACalc.HubbleRate(data.transitionTemp);
  data.gstar = BACalc.GetGstar(data.transitionTemp);
  data.FGW0  = 1.64 / pow(h, 2) * 1.e-5 *
              pow(100. / BACalc.GetGstar(data.transitionTemp), 1 / 3.);
  data.Hstar0 =
      GetHstar0(data.reheatingTemp, BACalc.GetGstar(data.transitionTemp));
  // Collisions
  data.pnlo_scaling = BACalc.pnlo_scaling;
  data.kappa_col    = kappa::Getkappa_col(
      data.transitionTemp, data.pnlo_scaling, data.HR, BACalc);
  // Turbulence
  data.Epsilon_Turb = BACalc.GetEpsTurb();
  // Sound waves
  const double alpha_eff =
      (1 - data.kappa_col) *
      data.PTStrength; // remove energy that goes into collisions
  data.vCJ = BACalc.GetChapmanJougetVelocity();
  std::vector<std::vector<double>> vprofile; // Fluid profile
  data.kappa_sw = (alpha_eff / data.PTStrength) *
                  kappa::kappaNuMuModel(pow(data.Csound_true, 2),
                                        pow(data.Csound_false, 2),
                                        alpha_eff,
                                        BACalc.vwall,
                                        vprofile);
  Logger::Write(LoggingLevel::GWDetailed, "HR = " + std::to_string(data.HR));
  Logger::Write(LoggingLevel::GWDetailed,
                "gstar = " + std::to_string(data.gstar));
  Logger::Write(LoggingLevel::GWDetailed, "T* = " + std::to_string(data.FGW0));
  Logger::Write(LoggingLevel::GWDetailed,
                "Hstar0 = " + std::to_string(data.Hstar0));
  Logger::Write(LoggingLevel::GWDetailed,
                "kappa_col = " + std::to_string(data.kappa_col));
  Logger::Write(LoggingLevel::GWDetailed,
                "Epsilon_Turb = " + std::to_string(data.Epsilon_Turb));
  Logger::Write(LoggingLevel::GWDetailed,
                "alpha_eff = " + std::to_string(alpha_eff));
  Logger::Write(LoggingLevel::GWDetailed, "vCJ = " + std::to_string(data.vCJ));
  Logger::Write(LoggingLevel::GWDetailed,
                "kappa_sw = " + std::to_string(data.kappa_sw));
  // XiShock = fastest fluid shell
  if (data.kappa_sw > 0)
  {
    data.XiShock = vprofile.back().front();
  }
  else
  {
    data.XiShock = 1.; // Does not matter as efficiency if zero
  }

  data.K_sw = GetK_sw(data.PTStrength, data.kappa_sw);
  Logger::Write(LoggingLevel::GWDetailed,
                "XiShock= " + std::to_string(data.XiShock));
  Logger::Write(LoggingLevel::GWDetailed,
                "K_sw = " + std::to_string(data.K_sw));
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

double GetHtauSH(const double HR, const double K_sw)
{
  return 2. / std::sqrt(3) * HR / std::sqrt(K_sw);
}

double GetHtauSW(const double HR, const double K_sw)
{
  return min(GetHtauSH(HR, K_sw), 1.);
}

double GetYpsilon(const double HR, const double K_sw)
{
  return 1. - 1. / sqrt(1. + 2. * GetHtauSW(HR, K_sw));
}

double GravitationalWave::CalcEpsTurb(double epsturb_in)
{
  if (epsturb_in == -1)
  {
    return std::sqrt(1 - GetYpsilon(data.HR, data.K_sw));
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

  const double f_p            = 0.11 * data.Hstar0 * data.betaH;
  data.CollisionParameter.f_b = f_p * pow(-n1 / n2, -1 / a1);

  // Calculate amplitude for collisions
  const double A_str  = 0.05;
  const double Ktilde = data.kappa_col * GetKtilde(data.PTStrength);
  const double Omega_p =
      data.FGW0 * A_str * pow(Ktilde, 2) / pow(data.betaH, 2);

  data.CollisionParameter.Omega_b =
      Omega_p * pow(1. / 2. * pow(-n2 / n1, n1 / (n1 - n2)) +
                        1. / 2. * pow(-n1 / n2, -n2 / (n1 - n2)),
                    (n1 - n2) / a1);

  Logger::Write(LoggingLevel::GWDetailed,
                "\n--------------- Collision ---------------\n");
  Logger::Write(LoggingLevel::GWDetailed,
                "f_b = " + std::to_string(data.CollisionParameter.f_b.value()));
  Logger::Write(LoggingLevel::GWDetailed,
                "Omega_b = " +
                    std::to_string(data.CollisionParameter.Omega_b.value()));
}

double GravitationalWave::CalculateXiShell()
{
  double XiFront, XiRear;
  if (data.vw >= data.vCJ)
    XiFront = data.vw;
  else
    XiFront = data.XiShock;

  if (data.vw < data.Csound_true)
    XiRear = data.vw;
  else
    XiRear = data.Csound_true;

  return XiFront - XiRear;
}

void GravitationalWave::CalcPeakSoundWave()
{
  data.SoundWaveParameter.n1 = 3.;
  data.SoundWaveParameter.n2 = 1.;
  data.SoundWaveParameter.n3 = -3.;
  data.SoundWaveParameter.a1 = 2.;
  data.SoundWaveParameter.a2 = 4.;

  // Characteristic frequencies

  const double xi_shell = CalculateXiShell();
  const double delta_w  = xi_shell / max(data.vw, data.Csound_false);

  const double f_1 = 0.2 * data.Hstar0 / data.HR;
  const double f_2 = 0.5 * data.Hstar0 / (data.HR * delta_w);

  // Sound wave amplitude

  const double Asw = 0.11; // Numerical simulation
  const double Omega_int =
      data.FGW0 * Asw * pow(data.K_sw, 2) *
      GetYpsilon(data.HR, data.K_sw) /* arxiv:2007.08537 */ * data.HR;

  // Convert Omega_int to Omega_2

  data.SoundWaveParameter.f_1 = f_1;
  data.SoundWaveParameter.f_2 = f_2;
  data.SoundWaveParameter.Omega_2 =
      Omega_int * (sqrt(2) + (2 * f_2 / f_1) / (1 + pow(f_2 / f_1, 2))) / M_PI;
  Logger::Write(LoggingLevel::GWDetailed,
                "\n--------------- Sound Wave ---------------\n");
  Logger::Write(LoggingLevel::GWDetailed,
                "f_1 = " + std::to_string(data.SoundWaveParameter.f_1.value()));
  Logger::Write(LoggingLevel::GWDetailed,
                "f_2 = " + std::to_string(data.SoundWaveParameter.f_2.value()));
  Logger::Write(LoggingLevel::GWDetailed,
                "Omega_2 = " +
                    std::to_string(data.SoundWaveParameter.Omega_2.value()));
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
  const double Omega_s = data.Epsilon_Turb * data.K_sw;

  data.TurbulanceParameter.f_1 =
      sqrt(3 * Omega_s) * data.Hstar0 / (2. * N * data.HR);
  data.TurbulanceParameter.f_2 = 2.2 * data.Hstar0 / data.HR;

  // Calculate Omega_2 for turbulence

  const double A_MHD = 3 * 2.2 * A / (4 * pow(M_PI, 2)) *
                       pow(2, -11 / (3 * data.TurbulanceParameter.a2.value()));

  data.TurbulanceParameter.Omega_2 =
      data.FGW0 * A_MHD * pow(Omega_s, 2) * pow(data.HR, 2);

  Logger::Write(LoggingLevel::GWDetailed,
                "\n--------------- Turbulence ---------------\n");
  Logger::Write(LoggingLevel::GWDetailed,
                "f_1 = " +
                    std::to_string(data.TurbulanceParameter.f_1.value()));
  Logger::Write(LoggingLevel::GWDetailed,
                "f_2 = " +
                    std::to_string(data.TurbulanceParameter.f_2.value()));
  Logger::Write(LoggingLevel::GWDetailed,
                "Omega_2 = " +
                    std::to_string(data.TurbulanceParameter.Omega_2.value()));
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
  double abs_err             = obj.AbsErr;
  double rel_err             = obj.RelErr;
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

double GetK_sw(const double &alpha, const double &kappa_sw)
{
  return 0.6 * kappa_sw * alpha / (1. + alpha);
}

double GetHstar0(const double &temp, const double &gstar)
{
  return 1.65 * 1e-5 * pow(gstar / 100, 1 / 6.) * (temp / 100);
}

double GetKtilde(const double &alpha)
{
  return alpha / (1. + alpha);
}

namespace kappa
{
// Compute kappa https://arxiv.org/abs/2010.09744

double mu(double a, double b)
{
  return (a - b) / (1.0 - a * b);
}

double getwow(double a, double b)
{
  return a / (1.0 - a * a) / b * (1.0 - b * b);
}

std::pair<double, ExpansionMode> getvm(double al, double vw, double cs2b)
{
  if (vw * vw < cs2b) // Deflagration
  {
    return {vw, ExpansionMode::Deflagration};
  }
  double cc   = 1.0 - 3.0 * al + vw * vw * (1.0 / cs2b + 3.0 * al);
  double disc = -4.0 * vw * vw / cs2b + cc * cc;
  if (disc < 0.0 || cc < 0.0)
  {
    return {std::sqrt(cs2b), ExpansionMode::Hybrid}; // Hybrid
  }
  return {(cc + std::sqrt(disc)) / 2.0 * cs2b / vw,
          ExpansionMode::Detonation}; // Detonation
}

// Differential equation system for dfdv
int dfdv(double v, const double y[], double dydv[], void *params)
{
  double cs2 = *(double *)params;
  double xi  = y[0];
  double w   = y[1];

  dydv[0] = (std::pow(mu(xi, v), 2) / cs2 - 1.0) * (1.0 - v * xi) * xi /
            (2.0 * v * (1.0 - v * v));
  dydv[1] = (1.0 + 1.0 / cs2) * mu(xi, v) * w / (1.0 - v * v);

  return GSL_SUCCESS;
}

// Solve ODE using GSL
std::vector<std::vector<double>> solve_ode(double vw, double v0, double cs2)
{
  const size_t dim          = 2;
  gsl_odeiv2_system sys     = {dfdv, nullptr, dim, &cs2};
  gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
      &sys, gsl_odeiv2_step_rkf45, -1e-10, 1e-10, 0.0);

  std::vector<std::vector<double>> results;

  double v    = v0;
  double y[2] = {vw, 1.}; // Initial conditions
  results.push_back({v, y[0], y[1]});

  const size_t steps = 64 * 1024;
  for (size_t i = 1; i < steps; ++i)
  {
    double v_next = v0 * (1.0 - static_cast<double>(i) / (steps - 1));
    if (i == steps - 1) v_next = v / 100.; // Avoid 1/0
    gsl_odeiv2_driver_apply(driver, &v, v_next, y);
    results.push_back({v_next, y[0], y[1]});
  }
  gsl_odeiv2_driver_free(driver);
  return results;
}

double integrate(const std::vector<double> &y, const std::vector<double> &x)
{
  if (y.size() != x.size() || y.size() < 2)
  {

    throw std::invalid_argument(
        "Vectors x and y must have the same size and "
        "contain at least two elements. The x and y sizes are " +
        std::to_string(x.size()) + " and " + std::to_string(y.size()));
  }

  // Check if x is in ascending or descending order
  bool is_descending = x.front() > x.back();

  // Create local copies if x is descending
  std::vector<double> x_sorted = x, y_sorted = y;
  if (is_descending)
  {
    std::reverse(x_sorted.begin(), x_sorted.end());
    std::reverse(y_sorted.begin(), y_sorted.end());
  }

  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(2000);
  gsl_function F;
  double result, error;

  // Interpolation function
  F.function = [](double xi, void *params) -> double
  {
    auto data = static_cast<
        std::pair<const std::vector<double> *, const std::vector<double> *> *>(
        params);
    const auto &Y = *data->first;
    const auto &X = *data->second;

    if (xi < X.front() || xi > X.back())
    {
      throw std::out_of_range("xi is out of bounds of the x vector.");
    }

    // Find the interval
    size_t idx = std::lower_bound(X.begin(), X.end(), xi) - X.begin();
    if (idx == 0) return Y[0];
    if (idx >= X.size()) return Y.back();

    // Linear interpolation
    double x1 = X[idx - 1], x2 = X[idx];
    double y1 = Y[idx - 1], y2 = Y[idx];
    return y1 + (xi - x1) * (y2 - y1) / (x2 - x1);
  };

  // Pass sorted data as parameters
  auto params = std::make_pair(&y_sorted, &x_sorted);
  F.params    = &params;

  // Integrate
  gsl_set_error_handler_off(); // Disable custom error handler if not defined
  gsl_integration_qags(&F,
                       x_sorted.front(),
                       x_sorted.back(),
                       1e-10,
                       1e-10,
                       2000,
                       workspace,
                       &result,
                       &error);

  // Clean up
  gsl_integration_workspace_free(workspace);

  // If x was descending, negate the result because integration limits are
  // reversed
  if (is_descending)
  {
    result = -result;
  }

  return result;
}

// Get K and Wow
std::pair<double, double> getKandWow(double vw,
                                     double v0,
                                     double cs2,
                                     std::vector<std::vector<double>> &vprofile)
{
  if (v0 == 0)
  {
    return {0, 1};
  }
  auto results = solve_ode(vw, v0, cs2);
  std::vector<double> vs, xis, wows, yis;

  for (const auto &r : results)
  {
    vs.push_back(r[0]);
    xis.push_back(r[1]);
    wows.push_back(r[2]);
  }

  if (mu(vw, v0) * vw <= cs2)
  {
    size_t ll = 1;
    for (size_t i = 0; i < xis.size(); ++i)
    {
      if (mu(xis[i], vs[i]) * xis[i] <= cs2)
      {
        ll = std::max(static_cast<size_t>(1), i + 1);
      }
    }

    vs.resize(std::min(ll + 1, vs.size()));
    xis.resize(std::min(ll + 1, vs.size()));
    wows.resize(std::min(ll + 1, vs.size()));

    for (size_t i = 0; i < wows.size(); ++i)
      wows[i] *= getwow(xis.back(), mu(xis.back(), vs.back())) / wows.back();
  }

  for (size_t i = 0; i < wows.size(); ++i)
    yis.push_back(wows[i] * pow(xis[i] * vs[i], 2) /
                  (1 - pow(vs[i], 2))); // w (\xi v)/gam^2

  vprofile = Transpose({xis, vs});

  if (xis.size() == 1) return {0, wows[0]};

  double Kint = integrate(yis, xis);
  return {Kint * 4.0 / std::pow(vw, 3), wows[0]};
}

// Remaining functions
double alN(double al, double wow, double cs2b, double cs2s)
{
  double da = (1.0 / cs2b - 1.0 / cs2s) / (1.0 / cs2s + 1.0) / 3.0;
  return (al + da) * wow - da;
}

std::pair<double, double>
getalNwow(double vp, double vm, double vw, double cs2b, double cs2s)
{
  std::vector<std::vector<double>> vprofile;
  auto [Ksh, wow] = getKandWow(vw, mu(vw, vp), cs2s, vprofile);
  double al = (vp / vm - 1.0) * (vp * vm / cs2b - 1.0) / (1.0 - vp * vp) / 3.0;
  return {alN(al, wow, cs2b, cs2s), wow};
}

void custom_error_handler(const char *reason,
                          const char *file,
                          int line,
                          int gsl_errno)
{
  std::cerr << "GSL Warning: " << reason << " at " << file << ":" << line
            << std::endl;
  (void)gsl_errno;
}

double kappaNuMuModel(double cs2b, double cs2s, double al, double vw)
{
  std::vector<std::vector<double>> vprofile;
  return kappaNuMuModel(cs2b, cs2s, al, vw, vprofile);
}

double kappaNuMuModel(double cs2b,
                      double cs2s,
                      double al,
                      double vw,
                      std::vector<std::vector<double>> &vprofile)
{
  if (vw == 1) vw = 0.999; // 1/0 if vw = 1, 0.999 is a decent approximation
  gsl_set_error_handler(&custom_error_handler);
  auto [vm, mode] = getvm(al, vw, cs2b); // Calculate vm = v-
  double Ksh = 0, wow = 1, vp = vw;
  double vpm;

  std::vector<std::vector<double>> vprofile_false, vprofile_true;

  if (mode == ExpansionMode::Deflagration or
      mode == ExpansionMode::Hybrid) // Deflagration or hybrid
  {
    auto [almax, wow2] = getalNwow(0, vm, vw, cs2b, cs2s);
    if (almax < al) return 0;

    vp                 = std::min(cs2s / vw, vw);
    auto [almin, wow1] = getalNwow(vp, vm, vw, cs2b, cs2s);

    if (almin > al) return 0;

    std::vector<std::vector<double>> iv = {{vp, almin}, {0, almax}};
    while (abs(iv[1][0] - iv[0][0]) > 1.e-7)
    {
      vpm        = (iv[1][0] + iv[0][0]) / 2.0;
      double alm = getalNwow(vpm, vm, vw, cs2b, cs2s).first;
      if (alm > al)
      {
        iv = {iv[0], {vpm, alm}};
      }
      else
      {
        iv = {{vpm, alm}, iv[1]};
      }
    }
    vp                 = (iv[1][0] + iv[0][0]) / 2.0;
    std::tie(Ksh, wow) = getKandWow(vw, mu(vw, vp), cs2s, vprofile_false);
  }
  double Krf = 0;
  if (mode == ExpansionMode::Hybrid or
      mode == ExpansionMode::Detonation) // Not deflagration
  {
    auto [Krf_val, wow3] = getKandWow(vw, mu(vw, vm), cs2b, vprofile_true);
    Krf                  = -wow * getwow(vp, vm) * Krf_val;
  }
  vprofile = vprofile_true;
  std::reverse(vprofile.begin(), vprofile.end());
  vprofile.insert(vprofile.end(), vprofile_false.begin(), vprofile_false.end());
  return (Ksh + Krf) / al;
}

double Getkappa_col(const double &Tstar,
                    const int &pnlo_scaling,
                    const double &HR,
                    BounceSolution &BACalc)
{
  // Calculate the false and true vacuum
  const std::vector<double> FalseVacuum = BACalc.modelPointer->MinimizeOrderVEV(
      BACalc.phase_pair.false_phase.Get(Tstar).point);
  const std::vector<double> TrueVacuum = BACalc.modelPointer->MinimizeOrderVEV(
      BACalc.phase_pair.true_phase.Get(Tstar).point);
  // Potential difference
  const double dV = BACalc.modelPointer->VEff(FalseVacuum, Tstar) -
                    BACalc.modelPointer->VEff(TrueVacuum, Tstar);
  // Potential energy contribution to the energy of the initial bubble
  const double E0V = BACalc.GetBounceSol(Tstar) * Tstar / 2.;
  // Initial bubble radius
  double R0 = pow(3 * E0V / (4 * M_PI * dV), 1. / 3.);

  // Mean bubble seperation
  const double Rstar = HR / BACalc.HubbleRate(Tstar);
  // Energy densitity radiation dominated Universe
  const double RhoGamma = BACalc.CalculateRhoGamma(Tstar);

  // Calculate masses squared in false and true vacuum
  std::vector<double> HiggsSq_False =
      BACalc.modelPointer->HiggsMassesSquared(FalseVacuum, Tstar);
  std::vector<double> LeptonSq_False =
      BACalc.modelPointer->LeptonMassesSquared(FalseVacuum);
  std::vector<double> QuarkSq_False =
      BACalc.modelPointer->QuarkMassesSquared(FalseVacuum);
  std::vector<double> GaugeSq_False =
      BACalc.modelPointer->GaugeMassesSquared(FalseVacuum, Tstar);

  std::vector<double> HiggsSq_True =
      BACalc.modelPointer->HiggsMassesSquared(TrueVacuum, Tstar);
  std::vector<double> LeptonSq_True =
      BACalc.modelPointer->LeptonMassesSquared(TrueVacuum);
  std::vector<double> QuarkSq_True =
      BACalc.modelPointer->QuarkMassesSquared(TrueVacuum);
  std::vector<double> GaugeSq_True =
      BACalc.modelPointer->GaugeMassesSquared(TrueVacuum, Tstar);

  double P_LO = 0; // Pressure at LO. 1 -> 1

  for (const auto &massSq : HiggsSq_True)
    P_LO += massSq;
  for (const auto &massSq : GaugeSq_True)
    P_LO += 3 * massSq; // 3 from polarization
  for (const auto &massSq : LeptonSq_True)
    P_LO += massSq / 2.;
  for (const auto &massSq : QuarkSq_True)
    P_LO += 3 * massSq / 2.; // 3 from colour
  for (const auto &massSq : HiggsSq_False)
    P_LO -= massSq;
  for (const auto &massSq : GaugeSq_False)
    P_LO -= 3 * massSq; // 3 from polarization
  for (const auto &massSq : LeptonSq_False)
    P_LO -= massSq / 2.;
  for (const auto &massSq : QuarkSq_False)
    P_LO -= 3 * massSq / 2.; // 3 from colour

  P_LO *= pow(Tstar, 2) / 24.; // pressure LO normalization

  // Pressure at NLO. 1 -> N
  // This only works for the SU(2)_L x U(1)_Y gauge group. But can be easily
  // generalized.
  const double W_coupling = BACalc.modelPointer->SMConstants.C_g / sqrt(2);
  const double Z_coupling =
      BACalc.modelPointer->SMConstants.C_g /
      sqrt(1. - BACalc.modelPointer->SMConstants.C_sinsquaredWeinberg);
  const double A_coupling =
      BACalc.modelPointer->SMConstants.C_g *
      sqrt(BACalc.modelPointer->SMConstants.C_sinsquaredWeinberg);
  double P_NLO = 0;
  if (pnlo_scaling == 1)
  {
    // Mass matrix is ordered
    P_NLO += sqrt(GaugeSq_True.at(0)) * pow(A_coupling, 2);
    P_NLO += 2 * sqrt(GaugeSq_True.at(1)) * pow(W_coupling, 2);
    P_NLO += sqrt(GaugeSq_True.at(3)) * pow(Z_coupling, 2);

    P_NLO -= sqrt(GaugeSq_False.at(0)) * pow(A_coupling, 2);
    P_NLO -= 2 * sqrt(GaugeSq_False.at(1)) * pow(W_coupling, 2);
    P_NLO -= sqrt(GaugeSq_False.at(3)) * pow(Z_coupling, 2);

    P_NLO *= pow(Tstar, 3);
  }
  else if (pnlo_scaling == 2)
  {
    P_NLO = 2 * pow(W_coupling, 2) + pow(Z_coupling, 2) + pow(A_coupling, 2);
    P_NLO *= pow(Tstar, 4); // pressure NLO normalization
  }
  const double alpha_infty    = P_LO / RhoGamma;
  const double gamma_run_away = Rstar / (3. * R0);
  if (BACalc.GetPTStrength() - alpha_infty < 0)
    return 0.; // Not enough pressure to drive the bubble wall.
  const double gamma_eq = pow((dV - P_LO) / P_NLO, 1. / pnlo_scaling);
  if (gamma_eq < 1) return 0.; // unphysical: dV - P_LO < P_NLO
  const double R_eq       = 3. * R0 * gamma_eq;
  const double gamma_star = min(gamma_eq, gamma_run_away);
  const double kappa_col  = (1 - alpha_infty / BACalc.GetPTStrength()) *
                           (1 - 1 / pow(gamma_eq, pnlo_scaling)) *
                           (R_eq / Rstar) * (gamma_star / gamma_eq);

  if (kappa_col < 0) return 0.; // If something unphysical happened
  return kappa_col;
}
} // namespace kappa
} // namespace BSMPT