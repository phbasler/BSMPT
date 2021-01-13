/*
 * IncludeAllModels.h
 *
 *  Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner

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

//#ifndef INCLUDEALLMODELS_H_
//#define INCLUDEALLMODELS_H_

#pragma once

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
class Class_Potential_Origin;
namespace ModelID
{

/**
 * @brief The ModelIDs enum containing all IDs for identifying the Models
 */
enum class ModelIDs
{
  NotSet,
  C2HDM,
  R2HDM,
  RN2HDM,
  CXSM,

  // Here you start adding your models
  TEMPLATE,

  // DO NOT EDIT the part below
  stop
};

/**
 * @brief Mapping between the model name which is given as the first argument to
 * the binary and the ModelIDs element
 */
const std::unordered_map<std::string, ModelIDs> ModelNames{
    {"c2hdm", ModelIDs::C2HDM},
    {"r2hdm", ModelIDs::R2HDM},
    {"n2hdm", ModelIDs::RN2HDM},
    {"cxsm", ModelIDs::CXSM},
    {"template", ModelIDs::TEMPLATE},
};

/**
 * @brief InvertModelNames
 * @return The switched map to ModelNames
 */
std::unordered_map<ModelIDs, std::string> InvertModelNames();

/**
 * @param choice ModelIDs for the Model under investigation
 * @return smart pointer to the instance of the class matching the ModelIDs
 * choice. If choice == NotSet the function throws an runtime error
 * @throw Runtime error if an invalid model was given into choice
 */

std::unique_ptr<Class_Potential_Origin> FChoose(ModelIDs choice);

/**
 *
 * @param s The input string, which is turned to lower case and then compared to
 * the entries of the ModelNames map
 * @return If a match in ModelNames is found, return the corresponding ModelIDs
 * entry, otherwise return ModelIDs::NoSet
 */
ModelIDs getModel(const std::string &s);
} // namespace ModelID

/**
 * @brief ShowInputError shows all the available models in the terminal
 */
void ShowInputError();

} // namespace BSMPT

//#endif /* INCLUDEALLMODELS_H_ */
