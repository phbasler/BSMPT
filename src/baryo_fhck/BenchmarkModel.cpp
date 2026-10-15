#include <BSMPT/baryo_fhck/BenchmarkModel.h>
#include <cmath>

namespace BSMPT
{
namespace Baryo
{
namespace FHCK
{

BenchmarkModel::BenchmarkModel(const double Tn_in,
                               const double vw_in,
                               const TruncationScheme &truncationscheme_in,
                               const double &truncationR_in)
    : TransportModel(nullptr,
                     {},
                     {},
                     vw_in,
                     Tn_in,
                     VevProfileMode::Kink,
                     truncationscheme_in,
                     truncationR_in)
    , vn(Tn_in)
    , wn(2. * Tn_in)
    , Tn(Tn_in)
    , LAM(1000.)
{
  Lw = 5. / Tn_in;
  Ls = Lw;
}

BenchmarkModel::BenchmarkModel(const double vn_in,
                               const double wn_in,
                               const double Tn_in,
                               const double LAM_in,
                               const double Lw_in,
                               const double Ls_in,
                               const double vw_in,
                               const TruncationScheme &truncationscheme_in,
                               const double &truncationR_in)
    : TransportModel(nullptr,
                     {},
                     {},
                     vw_in,
                     Tn_in,
                     VevProfileMode::Kink,
                     truncationscheme_in,
                     truncationR_in)
    , vn(vn_in)
    , wn(wn_in)
    , Tn(Tn_in)
    , LAM(LAM_in)
    , Ls(Ls_in)
{
  Lw = Lw_in;
}

void BenchmarkModel::Initialize()
{
  status = TransportModelStatus::Success;
  Logger::Write(LoggingLevel::FHCK, "Initialized");
}

void BenchmarkModel::GenerateFermionMass(const std::vector<double> &zList,
                                         const bool &MakeTopMassPlot)
{
  (void)zList;
  (void)MakeTopMassPlot;
  Logger::Write(LoggingLevel::FHCK, "FermionsMassesGenerated");
}

double BenchmarkModel::hvev(const double &z, const int &deriv)
{
  double res = 0.;
  if (deriv == 0)
    res += vn / 2. * (1 - tanh(z / Lw));
  else if (deriv == 1)
    res += -vn / (2. * Lw) / pow(cosh(z / Lw), 2);
  return res;
}

double BenchmarkModel::svev(const double &z, const int &deriv)
{
  double res = 0.;
  double u   = z / Ls;
  if (deriv == 0)
    res += wn / 2. * (1 + tanh(u));
  else if (deriv == 1)
    res += wn / (2. * Ls) / pow(cosh(u), 2);
  else if (deriv == 2)
    res += -wn / (Ls * Ls) * tanh(u) / pow(cosh(u), 2);
  return res;
}

void BenchmarkModel::GetFermionRatio(const double &z,
                                    const size_t &fermion,
                                    double &x2,
                                    double &x2prime,
                                    double &thetaprime,
                                    double &theta2prime)
{
  // bottom quark has no mass in this model
  if (fermion == 2)
  {
    x2          = 0.;
    x2prime     = 0.;
    thetaprime  = 0.;
    theta2prime = 0.;
  }
  else
  {
    double h    = hvev(z, 0);
    double dh   = hvev(z, 1);
    double s    = svev(z, 0) / LAM;
    double ds   = svev(z, 1) / LAM;
    double d2s  = svev(z, 2) / LAM;
    double temp = sqrt(1. + s * s);
    x2          = yt * h * temp;
    x2prime     = yt * ((1. + s * s) * dh + h * s * ds) / temp;
    x2prime     = 2. * x2 * x2prime / (Tstar * Tstar);
    x2 *= x2 / (Tstar * Tstar);
    thetaprime  = ds / (1 + s * s);
    theta2prime = ((1 + s * s) * d2s - 2. * s * ds * ds) / pow(1. + s * s, 2);
  }
}

double BenchmarkModel::GetWRatio(const double &z)
{
  double h = hvev(z, 0);
  return h * g / 2. / Tn;
}

double BenchmarkModel::EWSBVEV(const double &z)
{
  return hvev(z, 0);
}

} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT