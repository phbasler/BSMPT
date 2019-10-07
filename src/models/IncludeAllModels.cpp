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


/**
 * Checks if a string is only made up of numbers or not.
 */
bool is_number(const std::string& s)
{
    return !s.empty() && std::find_if(s.begin(),
        s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
}

/**
 * This function offers the possibility to call your model with a string. For that the the Model Input is converted to lower
 * case and then compared with the Model Name, the value of Model is then set to the C_Model Parameter.
 */
int getModel(const std::string& s){
	bool Debug=false;
	if(Debug) std::cout << "Start Debugging in " << __func__ << std::endl;
	int Model = -1;
	if(Debug) std::cout << "s =  " << s << std::endl;
	if(!is_number(s)){
		std::string ModelInput=s;
		if(Debug) std::cout << "ModelInput = " << ModelInput << std::endl;
		std::transform(ModelInput.begin(),ModelInput.end(),ModelInput.begin(),::tolower);

		if(ModelInput.compare("c2hdm")==0){
			// std::cout << "Found C2HDM as an Input. Setting Model = " << C_ModelC2HDM << "." << std::endl;
			Model=C_ModelC2HDM;
		}
		else if(ModelInput.compare("r2hdm")==0){
//			std::cout << "Found R2HDM as an Input. Setting Model = " << C_ModelR2HDM <<  "." << std::endl;
			Model=C_ModelR2HDM;
		}
		else if(ModelInput.compare("n2hdm")==0){
//			std::cout << "Found N2HDM as an Input. Setting Model = " << C_ModelRN2HDM << "." << std::endl;
			Model=C_ModelRN2HDM;
		}
	}
	else{
		Model = std::stoi(s);
	}

  if(Debug) std::cout << "Return " << Model << std::endl;
	if(Debug) std::cout << "End Debugging in " << __func__ << std::endl;

	return Model;
}

void ShowInputError(){
	std::cerr << "The chosen Method for the thermal mass corrections is ";
	if(C_UseParwani) std::cerr << "Parwani ";
	else std::cerr << "Arnold Espinosa\n";
	std::cerr << "The implemented models are \n"
			<< C_ModelC2HDM <<" : C2HDM\n"
			<< C_ModelR2HDM <<  " : R2HDM\n"
			<< C_ModelRN2HDM << " : N2HDM\n"
			<< std::endl;
}

