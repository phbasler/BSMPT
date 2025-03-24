// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <BSMPT/config.h>
#include <algorithm>
#include <functional>
#include <gsl/gsl_integration.h>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#ifdef Boost_FOUND
#include <boost/version.hpp>
#if BOOST_VERSION >= 107200
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#else
#include <boost/math/interpolators/cubic_b_spline.hpp>
#endif
#endif

/**
 * @file
 */
namespace BSMPT
{

/**
 * @brief Inverts a map
 * @param originalMap Map to be inverted
 * @param errorOnDuplicateValue Error message on duplicate value
 * @throw std::runtime_error if the new map would have duplicate keys
 */
template <typename key, typename value>
std::unordered_map<value, key>
InvertMap(const std::unordered_map<key, value> &originalMap,
          const std::string &errorOnDuplicateValue)
{
  std::unordered_map<value, key> result;
  for (const auto &[orig_key, orig_value] : originalMap)
  {
    auto success = result.emplace(orig_value, orig_key);
    if (not success.second)
    {
      throw std::runtime_error(errorOnDuplicateValue);
    }
  }

  return result;
}

/**
 * @brief StringStartsWith checks if str starts with prefix
 */
bool StringStartsWith(const std::string &str, const std::string &prefix);

/**
 * @brief StringEndsWith tests if str ends with suffix
 * @param str
 * @param suffix
 * @return
 */
bool StringEndsWith(const std::string &str, const std::string &suffix);

/**
 * @brief seperator used to write into output files
 */
const std::string sep = "\t";

/**
 * @brief factorial function
 */
int factorial(const int &a);

/**
 * @brief push back vector into vector
 */
template <typename T>
std::vector<T> push_back(std::vector<T> &a, const std::vector<T> &b)
{
  return a.insert(a.end(), b.begin(), b.end());
}

/**
 * @brief vector to_string
 */
template <typename T> std::string vec_to_string(const std::vector<T> &vec)
{
  std::string res;
  bool first = true;
  for (const auto &el : vec)
  {
    if (not first)
    {
      res += sep + std::to_string(el);
    }
    else
    {
      res   = std::to_string(el);
      first = false;
    }
  }
  return res;
}

/**
 * @brief split string separated by delimiter into substrings
 */
std::vector<std::string> split(const std::string &str, char delimiter);

/**
 * @brief Overload to print out vectors with the << operator
 */
template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec)
{
  bool first = true;
  for (const auto &el : vec)
  {
    if (not first)
    {
      os << sep;
    }
    else
    {
      first = false;
    }
    os << el;
  }
  return os;
}

/**
 * @brief vector addition
 */
template <typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b)
{
  if (a.size() != b.size())
    throw std::runtime_error(
        "Vector cannot be added. Must have the same size.");
  std::vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(),
                 a.end(),
                 b.begin(),
                 std::back_inserter(result),
                 std::plus<T>());
  return result;
}

/**
 * @brief vector subtraction
 */
template <typename T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b)
{
  if (a.size() != b.size())
    throw("Vector cannot be subtracted. Must have the same size.");

  std::vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(),
                 a.end(),
                 b.begin(),
                 std::back_inserter(result),
                 std::minus<T>());
  return result;
}

/**
 * @brief multiplication of vector with scalar
 */
template <typename T, typename T2>
std::vector<T> operator*(const T2 &a, const std::vector<T> &b)
{
  std::vector<T> result;
  result.reserve(b.size());

  std::transform(b.begin(),
                 b.end(),
                 std::back_inserter(result),
                 [&a](T i) { return a * i; });
  return result;
}

/**
 * @brief division of vector by scalar
 */
template <typename T, typename T2>
std::vector<T> operator/(const std::vector<T> &a, const T2 &b)
{
  std::vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(),
                 a.end(),
                 std::back_inserter(result),
                 [&b](T i) { return i / b; });
  return result;
}

/**
 * @brief dot product of two vectors
 */
template <typename T>
T operator*(const std::vector<T> &a, const std::vector<T> &b)
{
  if (a.size() != b.size())
    throw(
        "Dot product between vectors cannot be done. Must have the same size.");

  std::vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(),
                 a.end(),
                 b.begin(),
                 std::back_inserter(result),
                 [](T i, T j) { return (i * j); });

  T result1 = std::accumulate(result.begin(), result.end(), 0.0);

  return result1;
}

/**
 * @brief multiplication of matrix with vector
 */
template <typename T>
std::vector<T> operator*(const std::vector<std::vector<T>> &a,
                         const std::vector<T> &b)
{
  if (a.size() != b.size())
    throw("Multiplication of matrix with vector cannot be done. Must have the "
          "same size.");

  std::vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(),
                 a.end(),
                 std::back_inserter(result),
                 [&](std::vector<T> i) { return (i * b); });

  return result;
}

/**
 * @brief flatten matrix
 */
template <typename T>
std::vector<T> flatten(std::vector<std::vector<T>> const &vec)
{
  std::vector<T> flattened;
  for (auto const &v : vec)
  {
    flattened.insert(flattened.end(), v.begin(), v.end());
  }
  return flattened;
}

/**
 * @brief L2NormVector
 * @param vec vector
 * @return L2 norm of vector
 */
double L2NormVector(const std::vector<double> &vec);

/**
 * @brief Calculates the tranpose of a matrix
 * @param A matrix to be transposed
 * @return std::vector<std::vector<double>> transposed matrix
 */
std::vector<std::vector<double>>
Transpose(const std::vector<std::vector<double>> &A);

/**
 * @brief Dilogarithm of x
 *
 * https://en.wikipedia.org/wiki/Dilogarithm
 *
 * @param x real argument of from \f$ (-\infty, 1)\f$
 * @return double
 */
double Li2(const double &x);

/**
 * @brief Incomplete elliptic integral of the second kind of x with a different
 * parameterization and k^2 = -2
 *
 * \f$ \text{EllipIntSecond}(x) = -i E(i \phi, 2) = \int_0^\phi
 * \sqrt{1+2\sinh^2{\theta}}\,d\theta \f$
 *
 * https://en.wikipedia.org/wiki</Elliptic_integral#Incomplete_elliptic_integral_of_the_second_kind
 * https://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html
 *
 * @param x real argument
 * @return double
 */
double EllipIntSecond(const double &x);

#ifdef Boost_FOUND
#if BOOST_VERSION >= 107200
template <typename T>
using boost_cubic_b_spline =
    boost::math::interpolators::cardinal_cubic_b_spline<T>;
#else
template <typename T>
using boost_cubic_b_spline = boost::math::cubic_b_spline<T>;
#endif
#endif

} // namespace BSMPT
