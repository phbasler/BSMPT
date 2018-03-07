/*
 * IncludeModels.cpp
 *
 *  Created on: Mar 6, 2018
 *      Author: basler
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

