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


#include "ClassPotentialOrigin.h"

#include <memory>

/**
 * Here are the Models of the already implemented models
 */

const int C_ModelC2HDM = 0;
const int C_ModelR2HDM = 1;
const int C_ModelRN2HDM = 2;
const int C_ModelTemplate=5;



std::unique_ptr<Class_Potential_Origin> FChoose(int choice);

bool is_number(const std::string& s);
int getModel(const std::string& s);
void ShowInputError();


//#endif /* INCLUDEALLMODELS_H_ */
