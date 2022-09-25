// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <BSMPT/config.h>
#include <iostream>
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
 * Overload to print out vectors with the << operator
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
 * @brief operator << overload for the model parameter
 */
namespace ModelID
{
enum class ModelIDs;
}
std::ostream &operator<<(std::ostream &os, const ModelID::ModelIDs &Model);
std::string ModelIDToString(const ModelID::ModelIDs &Model);

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
