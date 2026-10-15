#define _USE_MATH_DEFINES
#include <BSMPT/baryo_fhck/Kernels.h>
#include <cmath>

namespace BSMPT
{
namespace Baryo
{
namespace FHCK
{

static const double KERNELPRECISION = 1e-6;

double ipow(const double x, const int n)
{
  double res = 1.;
  for (int i = 0; i < std::abs(n); i++)
    res *= x;
  return (n > 0 ? res : 1 / res);
}
double f0(const double x, const double s, const int diff)
{
  if (x > 50) return ipow(-1, diff) * std::exp(-x);

  const double expw = std::exp(-x / 2.);
  if (diff == 0)
    return expw / (1. / expw + s * expw);
  else if (diff == 1)
    return -1. / ipow(1. / expw + s * expw, 2);
  else if (diff == 2)
    return (1. / expw - s * expw) / ipow(1. / expw + s * expw, 3);
  return 0.;
}

void KernelInty::set_all(const double u)
{
  w   = x + (1 - u) / u;
  pwt = std::sqrt(w * w - x * x);
  pre = f0(w, s, k);
}

double KernelInty::operator()(const double y)
{
  const double pzt = gamw * (y * pwt - w * vw);
  const double Et  = gamw * (w - vw * y * pwt);
  double V         = pwt;
  if (structure == 0)
    V = pwt;
  else if (structure == 1)
    V = pwt / std::sqrt(1. + x * x / (pzt * pzt));
  else if (structure == 2)
    V = w * pzt * pzt / (pzt * pzt + x * x);
  return ipow(pzt, n) / ipow(Et, m - 1) * V;
}

double KernelIntw::operator()(const double u)
{
  Integrand.set_all(u);
  if (std::abs(Integrand.pre) == 0.) return 0.;
  return Integrand.pre *
         adap_gauss_kronrod_15(Integrand, -1, 1, KERNELPRECISION) / (u * u);
}

void Q9KernelInty::set_all(const double u)
{
  w    = x + (1 - u) / u;
  pwt  = std::sqrt(w * w - x * x);
  pre1 = pwt * f0(w, s, 1);
  pre2 = pwt * gamw * f0(w, s, 2);
}

double Q9KernelInty::operator()(const double y)
{
  const double pzt = gamw * (y * pwt - w * vw);
  const double Et  = gamw * (w - vw * y * pwt);
  double V         = 1;
  if (structure >= 1) V *= 1. / std::sqrt(1. + x * x / (pzt * pzt));
  if (structure >= 2) V = V * V / std::sqrt(1. - x * x / (w * w));
  return ipow(pzt, l - 2) / ipow(Et, l + 1) * V * (pre1 - pre2 * Et);
}

double Q9KernelIntw::operator()(const double u)
{
  Integrand.set_all(u);
  if ((std::abs(Integrand.pre1) == 0.) && (std::abs(Integrand.pre2) == 0.))
    return 0.;
  return adap_gauss_kronrod_15(Integrand, -1, 1, KERNELPRECISION) / (u * u);
}

double N0Int::operator()(const double u)
{
  const double w = x + (1 - u) / u;
  return std::sqrt(w * w - x * x) * w * f0(w, s, 0) / (u * u);
}

double RbarInt::operator()(const double u)
{
  const double w   = x + (1 - u) / u;
  const double pwt = std::sqrt(w * w - x * x);
  return M_PI / (gamw * gamw) * log(std::abs((pwt - vw * w) / (pwt + vw * w))) *
         f0(w, s, 0) / (u * u);
}

double K4Int::operator()(const double u)
{
  const double w   = x + (1 - u) / u;
  const double pwt = std::sqrt(w * w - x * x);
  return -2 / (M_PI * M_PI) * pwt * pwt * pwt / w * f0(w, s, 1) / (u * u);
}

double Kernel::operator()(const KernelType Kern,
                          const ParticleType P,
                          const double x,
                          const double vw)
{
  const double statistic = (P == ParticleType::Boson ? -1. : 1.);
  double gamw            = 1. / std::sqrt(1. - vw * vw);
  double res             = -3 / (M_PI * M_PI * gamw);

  switch (Kern)
  {
  case KernelType::K:
  {
    KernelIntw integrand(0, 0, l, l, vw, gamw, statistic, x);
    res *= -adap_gauss_kronrod_15(integrand, 0., 1., KERNELPRECISION);
  }
  break;
  case KernelType::D:
  {
    KernelIntw integrand(1, 0, l, l, vw, gamw, statistic, x);
    res *= adap_gauss_kronrod_15(integrand, 0., 1., KERNELPRECISION);
  }
  break;
  case KernelType::Q:
  {
    KernelIntw integrand(2, 0, l - 1, l, vw, gamw, statistic, x);
    res *= adap_gauss_kronrod_15(integrand, 0., 1., KERNELPRECISION) / 2.;
  }
  break;
  case KernelType::Qe:
  {
    KernelIntw integrand(1, 0, l - 1, l, vw, gamw, statistic, x);
    res *= adap_gauss_kronrod_15(integrand, 0., 1., KERNELPRECISION) / 2.;
  }
  break;
  case KernelType::Q8o:
  {
    KernelIntw integrand(1, structure, l - 2, l, vw, gamw, statistic, x);
    res *= adap_gauss_kronrod_15(integrand, 0., 1., KERNELPRECISION) / 2.;
  }
  break;
  case KernelType::Q9o:
  {
    Q9KernelIntw integrand(structure, l, vw, gamw, statistic, x);
    res *= adap_gauss_kronrod_15(integrand, 0., 1., KERNELPRECISION) / 4.;
  }
  break;
  case KernelType::Rb:
  {
    N0Int integrand1(statistic, x);
    const double res1 =
        4. * M_PI * adap_gauss_kronrod_15(integrand1, 0., 1., KERNELPRECISION);
    RbarInt integrand2(vw, gamw, statistic, x);
    const double res2 =
        adap_gauss_kronrod_15(integrand2, 0., 1., KERNELPRECISION);
    return res2 / res1;
  }
  break;
  case KernelType::K4FH:
  {
    K4Int integrand(statistic, x);
    const double est = kronrod_61(integrand, 0., 1.);
    return h_adap_gauss_kronrod_15(integrand, 0., 1., est, 1e-4);
  }
  break;

  default: std::cout << "This kernel type is not supported\n"; break;
  }

  return res;
}

} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT