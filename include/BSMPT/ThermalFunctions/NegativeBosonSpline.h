// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#define _USE_MATH_DEFINES
#include <BSMPT/utility/spline/spline.h>
#include <cmath>
#include <vector>

#ifndef INCLUDE_BSMPT_THERMALFUNCTIONS_NEGATIVEBOSONSPLINE_H_
#define INCLUDE_BSMPT_THERMALFUNCTIONS_NEGATIVEBOSONSPLINE_H_

/**
 * @brief C_NegLine Number of data points used for the interpolation of J_(m^2 <
 * 0)
 */
extern const int C_NegLine;

/**
 * @brief NegLinearInt 2D Array containing the pairs (m^2 , J_(m^2)) for m^2 < 0
 */
extern const double NegLinearInt[3001][2];

/**
 * @brief NegLinearInt Transpose of NegLinearInt
 */
extern const std::vector<std::vector<double>> NegLinearIntTransposed;

/**
 * @brief Cubic spline with the \f$ J_b(\frac{m}{T}) \f$ for \f$ \frac{m}{T} < 0
 * \f$
 *
 */
extern const tk::spline JbosonNegativeSpline;

#endif /* INCLUDE_BSMPT_THERMALFUNCTIONS_NEGATIVEBOSONSPLINE_H_ */
