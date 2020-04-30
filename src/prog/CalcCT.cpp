/*
 * CalcCT.cpp
 *
 *
 *      Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller

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
 * Calculates the Counterterms for mulitple points and adds them at the end of the line. The sequence is the
 * same as used in Potential::set_CT_Pot_Par
 */

#include <bits/exception.h>                     // for exception
#include <stdlib.h>                             // for atoi, EXIT_FAILURE
#include <memory>                               // for unique_ptr
#include <string>                               // for operator<<, getline
#include <utility>                              // for pair
#include <vector>                               // for vector
#include <BSMPT/models/ClassPotentialOrigin.h>  // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility.h>
#include <fstream>
#include <iostream>
using namespace std;
using namespace BSMPT;
//#include "Minimizer.h"
int main(int argc, char *argv[]) try{

	if(!( argc == 6) )
	{
		std::cerr << "./CalcCT Model Inputfile Outputfile LineStart LineEnd \n";
		ShowInputError();
		return EXIT_FAILURE;
	}

	double LineStart,LineEnd;
	char* in_file;char* out_file;

    auto Model=ModelID::getModel(argv[1]);
    if(Model==ModelID::ModelIDs::NotSet) {
        std::cerr << "Your Model parameter does not match with the implemented Models." << std::endl;
        ShowInputError();
        return EXIT_FAILURE;
    }


	in_file = argv[2];
	out_file = argv[3];




	LineStart = atoi(argv[4]);
	LineEnd = atoi(argv[5]);

	if(LineStart < 1)
	{
		std::cerr << "Start line counting with 1" << std::endl;
		return EXIT_FAILURE;
	}
	if(LineStart > LineEnd)
	{
		std::cerr << "LineEnd is smaller then LineStart " << std::endl;
		return EXIT_FAILURE;
	}


	int linecounter = 1;
	std::ifstream infile(in_file);
	if(!infile.good()) {
			std::cout << "Input file not found " << std::endl;
			return EXIT_FAILURE;
	}
	std::ofstream outfile(out_file);
	if(!outfile.good())
	{
		std::cout << "Can not create file " << out_file << std::endl;
		return EXIT_FAILURE;
	}
	std::string linestr;

    std::unique_ptr<Class_Potential_Origin> modelPointer = ModelID::FChoose(Model);
    size_t nPar,nParCT;
    nPar = modelPointer->get_nPar();
    nParCT = modelPointer->get_nParCT();
	std::vector<double> par(nPar);
	std::vector<double> parCT(nParCT);




	while(getline(infile,linestr))
	{


		if(linecounter > LineEnd) break;
		if(linecounter == 1)
		  {
			modelPointer->setUseIndexCol(linestr);
            outfile << linestr;
            for(auto x: modelPointer->addLegendCT()) outfile << sep << x;
		    outfile << std::endl;
		  }

		if(linecounter >= LineStart and linecounter <= LineEnd and linecounter != 1)
		{
			std::pair<std::vector<double>,std::vector<double>> parameters = modelPointer->initModel(linestr);
			par=parameters.first;
			parCT = parameters.second;

            outfile << linestr << sep << parCT << std::endl;
        }
		linecounter++;
		if(infile.eof()) break;
	}
	if(LineStart == LineEnd) {
		std::vector<double> DebugIn,DebugOut;
		modelPointer->Debugging(DebugIn,DebugOut);
	}

	outfile.close();


	return EXIT_SUCCESS;
}
catch(exception& e){
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
}
