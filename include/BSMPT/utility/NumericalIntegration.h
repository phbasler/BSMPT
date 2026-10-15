#pragma once

#include "abscissa.h"
#include <iostream>

template <class FUNC> double kronrod_61(FUNC &f, const double l, const double r)
{
  double m = 0.5 * (r + l);
  double h = 0.5 * (r - l);
  double res{0};

  for (size_t i = 0; i < 30; i++)
  {
    double dx = h * kronx_61[i];
    res += wkron_61[i] * (f(m + dx) + f(m - dx));
  }
  res += wkron_61[30] * f(m);
  return h * res;
}

template <class FUNC>
double adap_gauss_kronrod_15(FUNC &f,
                             const double l,
                             const double r,
                             const double err,
                             const size_t depth = 0)
{
  double I1, I2, y[15];
  double h = (r - l) / 2;
  for (int i = 0; i < 15; i++)
  {
    y[i] = f((kronx_15[i] + 1) * h + l);
  }
  I1 = h * (0.022935322010529224963732008058970 * (y[0] + y[14]) +
            0.063092092629978553290700663189204 * (y[1] + y[13]) +
            0.104790010322250183839876322541518 * (y[2] + y[12]) +
            0.140653259715525918745189590510238 * (y[3] + y[11]) +
            0.169004726639267902826583426598550 * (y[4] + y[10]) +
            0.190350578064785409913256402421014 * (y[5] + y[9]) +
            0.204432940075298892414161999234649 * (y[6] + y[8]) +
            0.209482141084727828012999174891714 * y[7]);
  if ((I1 == 0) || (depth > 20))
  {
    return I1;
  }
  I2 = h * (0.129484966168869693270611432679082 * (y[1] + y[13]) +
            0.279705391489276667901467771423780 * (y[3] + y[11]) +
            0.381830050505118944950369775488975 * (y[5] + y[9]) +
            0.417959183673469387755102040816327 * y[7]);

  if (std::abs((I2 / I1 - 1)) < err)
  {
    return I1;
  }
  double m = (l + r) / 2;
  return adap_gauss_kronrod_15(f, l, m, err, depth + 1) +
         adap_gauss_kronrod_15(f, m, r, err, depth + 1);
}

template <class FUNC>
double h_adap_gauss_kronrod_15(FUNC &f,
                               const double l,
                               const double r,
                               const double est,
                               const double err,
                               const size_t depth = 0)
{
  double I1, I2, y[15];
  double h = (r - l) / 2;
  for (int i = 0; i < 15; i++)
  {
    y[i] = f((kronx_15[i] + 1) * h + l);
  }
  I1 = h * (0.022935322010529224963732008058970 * (y[0] + y[14]) +
            0.063092092629978553290700663189204 * (y[1] + y[13]) +
            0.104790010322250183839876322541518 * (y[2] + y[12]) +
            0.140653259715525918745189590510238 * (y[3] + y[11]) +
            0.169004726639267902826583426598550 * (y[4] + y[10]) +
            0.190350578064785409913256402421014 * (y[5] + y[9]) +
            0.204432940075298892414161999234649 * (y[6] + y[8]) +
            0.209482141084727828012999174891714 * y[7]);
  if ((I1 == 0) || (depth > 20))
  {
    return I1;
  }
  I2 = h * (0.129484966168869693270611432679082 * (y[1] + y[13]) +
            0.279705391489276667901467771423780 * (y[3] + y[11]) +
            0.381830050505118944950369775488975 * (y[5] + y[9]) +
            0.417959183673469387755102040816327 * y[7]);

  if (std::abs((I1 - I2)) < err * std::abs(est))
  {
    return I1;
  }
  double m = (l + r) / 2;
  return h_adap_gauss_kronrod_15(f, l, m, est, err, depth + 1) +
         h_adap_gauss_kronrod_15(f, m, r, est, err, depth + 1);
}

template <class FUNC>
double adap_simpson38(FUNC &f,
                      const double l,
                      const double r,
                      double *f0,
                      const double err,
                      const size_t depth = 0)
{
  double I1, I2, f1[4];
  double m  = (r + l) / 2.;
  double h  = (r - l) / 8.;
  double Ia = h * (f0[0] + 3 * f0[1] + 3 * f0[2] + f0[3]);
  f1[0]     = f(m);
  f1[1]     = f0[2];
  f1[2]     = f((l + 5 * r) / 6);
  f1[3]     = f0[3];
  f0[3]     = f1[0];
  f0[2]     = f0[1];
  f0[1]     = f((5 * l + r) / 6);
  I1        = h / 2 * (f0[0] + 3 * f0[1] + 3 * f0[2] + f0[3]);
  I2        = h / 2 * (f1[0] + 3 * f1[1] + 3 * f1[2] + f1[3]);
  double Ib = I1 + I2;

  if ((std::abs(Ia / Ib - 1) < err) || (depth > 16))
  {
    return Ib;
  }
  return adap_simpson38(f, l, m, f0, err, depth + 1) +
         adap_simpson38(f, m, r, f1, err, depth + 1);
}

template <class FUNC>
double h_adap_simpson38(FUNC &f,
                        const double l,
                        const double r,
                        double *f0,
                        const double est,
                        const double err,
                        const size_t depth = 0)
{
  double I1, I2, f1[4];
  double m  = (r + l) / 2.;
  double h  = (r - l) / 8.;
  double Ia = h * (f0[0] + 3 * f0[1] + 3 * f0[2] + f0[3]);
  f1[0]     = f(m);
  f1[1]     = f0[2];
  f1[2]     = f((l + 5 * r) / 6);
  f1[3]     = f0[3];
  f0[3]     = f1[0];
  f0[2]     = f0[1];
  f0[1]     = f((5 * l + r) / 6);
  I1        = h / 2 * (f0[0] + 3 * f0[1] + 3 * f0[2] + f0[3]);
  I2        = h / 2 * (f1[0] + 3 * f1[1] + 3 * f1[2] + f1[3]);
  double Ib = I1 + I2;

  if ((std::abs(Ia - Ib) < err * std::abs(est)) || (depth > 16))
  {
    return Ib;
  }
  return h_adap_simpson38(f, l, m, f0, est, err, depth + 1) +
         h_adap_simpson38(f, m, r, f1, est, err, depth + 1);
}