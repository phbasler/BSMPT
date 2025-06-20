// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/utility.h>
#include <complex>
#include <ostream>
#include <sstream>
#include <string>

/**
 * @file
 */
namespace BSMPT
{

std::vector<std::string> split(const std::string &str, char delimiter)
{
  // Using str in a string stream
  std::stringstream ss(str);
  std::vector<std::string> res;
  std::string token;
  while (getline(ss, token, delimiter))
  {
    res.push_back(token);
  }
  return res;
}

bool StringStartsWith(const std::string &str, const std::string &prefix)
{
  return str.size() >= prefix.size() and str.find(prefix) == 0;
}

int factorial(const int &a)
{
  return (a == 1 || a == 0) ? 1 : factorial(a - 1) * a;
}

double L2NormVector(const std::vector<double> &vec)
{
  double r = 0.0;
  int dim  = vec.size();
  for (int i = 0; i < dim; i++)
  {
    r += vec[i] * vec[i];
  }
  return std::sqrt(r);
}

std::vector<std::vector<double>>
Transpose(const std::vector<std::vector<double>> &A)
{
  int rows = A.size();
  if (rows == 0) return {{}};
  int cols = A[0].size();
  std::vector<std::vector<double>> r(cols, std::vector<double>(rows));
  for (int i = 0; i < rows; ++i)
  {
    for (int j = 0; j < cols; ++j)
    {
      r[j][i] = A[i][j];
    }
  }
  return r;
}

double Li2(const double &x)
{
  if (x == 0) return 0;
  if (x == 1) return pow(M_PI, 2) / 6.;
  if (x < -1) return -pow(M_PI, 2) / 6. - pow(log(-x), 2) / 2. - Li2(1. / x);
  if (x < 0) return 1 / 2. * Li2(-x * -x) - Li2(-x);
  if (x > 0.5) return pow(M_PI, 2) / 6. - log(x) * log(1 - x) - Li2(1 - x);
  double sum = x;                // k = 1 element of the sum
  for (int k = 2; k <= 1e5; k++) // Sum starts at k = 2
  {
    if (abs((pow(x, k) / pow(k, 2.)) / sum) < 1e-10)
    {
      sum += pow(x, k) / pow(k, 2.);
      return sum;
    }
    sum += pow(x, k) / pow(k, 2.);
  }
  return sum;
}

bool StringEndsWith(const std::string &str, const std::string &suffix)
{
  return str.size() >= suffix.size() and
         str.substr(str.size() - suffix.size(), str.size()) == suffix;
}

double EllipIntSecond(const double &x)
{
  std::function<double(double)> integrand = [&](double x_int)
  { return sqrt(1 + 2 * pow(sinh(x_int), 2)); };

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

  double result, error;

  gsl_function F = {[](double d, void *vf) -> double
                    {
                      auto &f =
                          *static_cast<std::function<double(double)> *>(vf);
                      return f(d);
                    },
                    &integrand};

  gsl_integration_qags(&F, 0, x, 0, 1e-7, 1000, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result;
}

bool almost_the_same(const double &a,
                     const double &b,
                     const double &rel_precision,
                     const double &num_zero)
{
  if (std::abs(a) < num_zero and std::abs(b) < num_zero)
  {
    return true;
  }
  return std::abs(a - b) < std::abs(a + b) / 2 * rel_precision;
}

bool almost_the_same(const std::complex<double> &a,
                     const std::complex<double> &b,
                     const double &rel_precision,
                     const double &num_zero)
{
  bool real_part = almost_the_same(a.real(), b.real(), rel_precision, num_zero);
  bool imag_part = almost_the_same(a.imag(), b.imag(), rel_precision, num_zero);
  return (real_part and imag_part);
}

bool almost_the_same(const std::vector<double> &a,
                     const std::vector<double> &b,
                     const bool &allow_for_sign_flip,
                     const double &rel_precision,
                     const double &num_zero)
{
  if (a.size() != b.size())
  {
    throw std::runtime_error("Error. Vectors must have the same size.");
  }
  int count_true = 0;
  for (std::size_t i = 0; i < a.size(); i++)
  {
    if (allow_for_sign_flip)
    {
      count_true +=
          int(almost_the_same(a.at(i), b.at(i), rel_precision, num_zero));
    }
    else
    {
      count_true += int(almost_the_same(
          std::abs(a.at(i)), std::abs(b.at(i)), rel_precision, num_zero));
    }
  }
  if (std::size_t(count_true) == a.size())
  {
    return true;
  }
  else
  {
    return false;
  }
}

} // namespace BSMPT
