// Copyright (C) 2025 Lisa Biermann, Margarete Mühlleitner, Rui Santos, João
// Viana
//
// SPDX-FileCopyrightText: 2025 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file spline to calculate gstar
 */

#pragma once
#include <vector>

namespace BSMPT
{

/**
 * @brief List of temperatures below T = 214 GeV for the gstar spline
 * construction
 *
 */
extern const std::vector<double> TGstarLowT;
/**
 * @brief  List of temperatures above T = 214 GeV for the gstar spline
 * construction
 */
extern const std::vector<double> TGstarHighT;
/**
 * @brief List of gstar below T = 214 GeV for the gstar spline
 * construction
 */
extern const std::vector<double> GstarLowT;
/**
 * @brief List of gstar above T = 214 GeV for the gstar spline
 * construction
 */
extern const std::vector<double> GstarHighT;
} // namespace BSMPT
