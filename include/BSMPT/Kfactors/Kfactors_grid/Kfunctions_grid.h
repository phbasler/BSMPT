/*
* ../include/BSMPT/Kfactors/Kfactors_grid/Kfunctions_grid.h
*
*  Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller

       This program is free software: you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published by
       the Free Software Foundation, either version 3 of the License, or
       (at your option) any later version.

       This program is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @file
 * Declaration of the vectors containing the data points used for the Kfactor
 * interpolations
 */
#pragma once
#include <array>
namespace BSMPT
{
namespace Kfactors
{
namespace Data
{
/**
 * @brief msg_size number of data points sampled in msg
 */
const std::size_t msg_size = 12999;
/**
 * @brief Tg_size number of data points sampled in Tg
 */
const std::size_t Tg_size = 121;
/**
 * @brief msg Data points used to sample in the m^2 direction of the Kfunctions
 * starting from 0 to 5 in steps of 10^-3 and then continuing to 200^2 in steps
 * of 5
 */
extern const std::array<double, msg_size> msg;
/**
 * @brief Tg Data points used to sample in the Temperature direction of the
 * Kfunctions starting from 10 to 250 in steps of 2
 */
extern const std::array<double, Tg_size> Tg;
/**
 * @brief K1p Data set used for the interpolation of the K1 functions for
 * fermions
 */
extern const std::array<std::array<double, Tg_size>, msg_size> K1p;
/**
 * @brief K1m Data set used for the interpolation of the K1 functions for bosons
 */
extern const std::array<std::array<double, Tg_size>, msg_size> K1m;
/**
 * @brief K2p Data set used for the interpolation of the K2 functions for
 * fermions
 */
extern const std::array<std::array<double, Tg_size>, msg_size> K2p;
/**
 * @brief K4p Data set used for the interpolation of the K4 functions for
 * fermions
 */
extern const std::array<std::array<double, Tg_size>, msg_size> K4p;
/**
 * @brief K4m Data set used for the interpolation of the K4 functions for bosons
 */
extern const std::array<std::array<double, Tg_size>, msg_size> K4m;
/**
 * @brief K5p Data set used for the interpolation of the K5 functions for
 * fermions
 */
extern const std::array<std::array<double, Tg_size>, msg_size> K5p;
/**
 * @brief K5m Data set used for the interpolation of the K5 functions for bosons
 */
extern const std::array<std::array<double, Tg_size>, msg_size> K5m;
/**
 * @brief K6p Data set used for the interpolation of the K6 functions for
 * fermions
 */
extern const std::array<std::array<double, Tg_size>, msg_size> K6p;
/**
 * @brief K8p Data set used for the interpolation of the K8 functions for
 * fermions
 */
extern const std::array<std::array<double, Tg_size>, msg_size> K8p;
/**
 * @brief K9p Data set used for the interpolation of the K9 functions for
 * fermions
 */
extern const std::array<std::array<double, Tg_size>, msg_size> K9p;
} // namespace Data
} // namespace Kfactors
} // namespace BSMPT
