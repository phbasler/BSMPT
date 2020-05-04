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


#include <BSMPT/utility.h>
#include <string>
#include <BSMPT/models/IncludeAllModels.h>

/**
 * @file
 */
namespace BSMPT {

std::ostream& operator<<(std::ostream& os, const ModelID::ModelIDs& Model)
{
    static auto IMN = BSMPT::ModelID::InvertModelNames();
    os << IMN.at(Model);
    return os;
}
bool StringStartsWith(const std::string& str, const std::string& prefix)
{
    return str.size() >= prefix.size() and str.find(prefix) == 0;
}

}
