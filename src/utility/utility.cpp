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

std::ostream &operator<<(std::ostream &os, const ModelID::ModelIDs &Model)
{
  static auto IMN = BSMPT::ModelID::InvertModelNames();
  os << IMN.at(Model);
  return os;
}

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

std::string ModelIDToString(const ModelID::ModelIDs &Model)
{
  std::stringstream ss;
  ss << Model;
  return ss.str();
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

bool StringEndsWith(const std::string &str, const std::string &suffix)
{
  return str.size() >= suffix.size() and
         str.substr(str.size() - suffix.size(), str.size()) == suffix;
}

} // namespace BSMPT
