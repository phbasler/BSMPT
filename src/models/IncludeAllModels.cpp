/*
 * IncludeModels.cpp
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

#include "IncludeAllModels.h"
#include "ClassPotentialR2HDM.h"
#include "ClassPotentialC2HDM.h"
#include "ClassPotentialRN2HDM.h"


/**
 * This method calls the pointer to the specific model.
 * If you add a model you have to add this case here as well.
 */

std::unique_ptr<Class_Potential_Origin> FChoose(int choice){
	if(choice == C_ModelR2HDM)
	{
		return std::unique_ptr<Class_Potential_Origin> { new Class_Potential_R2HDM};
	}
	else if(choice == C_ModelC2HDM)
	{
		return std::unique_ptr<Class_Potential_Origin> { new Class_Potential_C2HDM};
	}
	else if(choice == C_ModelRN2HDM){
		return std::unique_ptr<Class_Potential_Origin> { new Class_Potential_RN2HDM};
	}
//	else if(choice == C_ModelTemplate)
//	{
//		return std::unique_ptr<Class_Potential_Origin> { new Class_Template };
//	}

	throw std::runtime_error("Invalid model");
}

