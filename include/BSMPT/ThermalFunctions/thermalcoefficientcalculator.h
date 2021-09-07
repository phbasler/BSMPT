// Copyright (C) 2021  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <functional>
#include <map>
namespace BSMPT
{
namespace ThermalFunctions
{
/**
 * @brief The ThermalCoefficientCalculator class is a thread-safe wrapper around
 * the calculation of the coefficients for the thermal expansions used in
 * ThermalFunctions
 */
class ThermalCoefficientCalculator
{
public:
  /**
   * @brief ThermalCoefficientCalculator
   * @param func Describes the function to calculate the coefficents
   * @param maxOrderToPrecalc defines until which order the coefficents are
   * precalculated
   */
  ThermalCoefficientCalculator(std::function<double(int)> func,
                               int maxOrderToPrecalc);
  /**
   * @brief GetCoefficentAtOrder
   * @param n describes the order at which the coefficient should be calculated
   * @return The coefficent at the given order. If it was precalculated the
   * result stored in the map will be returned, otherwise the result will be
   * calculated.
   */
  [[nodiscard]] double GetCoefficentAtOrder(int n) const;

private:
  /**
   * @brief Calculater stores the function to calculate the coefficient
   */
  std::function<double(int)> Calculater;
  /**
   * @brief MaxOrderToSave defines the highest order until which the coefficents
   * are precalculated
   */
  const int MaxOrderToSave;

  /**
   * @brief PreCalculatedCoefficents Stores the precalculated coefficents. This
   * will only be changed during the constructor to remain thread-safe!
   */
  std::map<int, double> PreCalculatedCoefficents;
};

} // namespace ThermalFunctions
} // namespace BSMPT
