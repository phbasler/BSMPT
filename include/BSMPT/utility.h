/*
 * utility.h
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

#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <iostream>
#include <vector>
#include <random>

/**
 * @file
 */
namespace BSMPT {

/**
 * @brief StringStartsWith checks if str starts with prefix
 */
bool StringStartsWith(const std::string& str, const std::string& prefix);

/**
 * @brief seperator used to write into output files
 */
const std::string sep = "\t";

/**
 * Overload to print out vectors with the << operator
 */
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
    bool first = true;
    for (const auto& el : vec)
    {
        if(not first) {
            os << sep;
        }
        else{
            first = false;
        }
        os << el;
    }
    return os;
}


/**
 * @brief operator << overload for the model parameter
 */
namespace ModelID {
  enum class ModelIDs;
}
std::ostream& operator<<(std::ostream& os, const ModelID::ModelIDs& Model);


}

#endif // UTILITY_H
