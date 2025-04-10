// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

// #ifndef INCLUDEALLMODELS_H_
// #define INCLUDEALLMODELS_H_

#pragma once

#include <BSMPT/utility/ModelIDs.h>
#include <algorithm> // for max
#include <memory>
#include <string> // for string
#include <unordered_map>
#include <vector>

/**
 * @file
 */
namespace BSMPT
{
struct ISMConstants;
class Class_Potential_Origin;
namespace ModelID
{

/**
 * @param choice ModelIDs for the Model under investigation
 * @return smart pointer to the instance of the class matching the ModelIDs
 * choice. If choice == NotSet the function throws an runtime error
 * @throw Runtime error if an invalid model was given into choice
 */
[[deprecated(
    "Will call FChoose with GetSMConstants(). Please use the detailed overload "
    "to ensure consistent SM constants through all routines.")]] std::
    unique_ptr<Class_Potential_Origin>
    FChoose(ModelIDs choice);

/**
 * @param choice ModelIDs for the Model under investigation
 * @return smart pointer to the instance of the class matching the ModelIDs
 * choice. If choice == NotSet the function throws an runtime error
 * @param smConstants The SM Constants to use for the parameter
 * @throw Runtime error if an invalid model was given into choice
 */

std::unique_ptr<Class_Potential_Origin>
FChoose(ModelIDs choice, const ISMConstants &smConstants);

/**
 *
 * @param s The input string, which is turned to lower case and then compared to
 * the entries of the ModelNames map
 * @return If a match in ModelNames is found, return the corresponding ModelIDs
 * entry, otherwise return ModelIDs::NoSet
 */
ModelIDs getModel(const std::string &s);
} // namespace ModelID

} // namespace BSMPT

// #endif /* INCLUDEALLMODELS_H_ */
