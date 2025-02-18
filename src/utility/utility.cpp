// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/utility.h>
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

} // namespace BSMPT
