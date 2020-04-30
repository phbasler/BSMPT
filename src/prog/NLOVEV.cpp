/*
 * NLOVEV.cpp
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
 * Calculates the VEV at T = 0 at NLO for a given Inputfile for a given subset of lines in the file and
 *  adds it at the end of the line . One parameter point per line.
 *
 */

#include <bits/exception.h>                     // for exception
#include <stdlib.h>                             // for atoi, EXIT_FAILURE
#include <algorithm>                            // for copy, max
#include <memory>                               // for shared_ptr, __shared_...
#include <string>                               // for string, operator<<
#include <utility>                              // for pair
#include <vector>                               // for vector
#include <BSMPT/models/ClassPotentialOrigin.h>  // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/utility.h>
#include <iostream>
#include <fstream>
using namespace std;
using namespace BSMPT;



//#include "Minimizer.h"

int main(int argc, char *argv[]) try{

	if(!( argc == 6) )
	{
		std::cerr << "./NLOVEV Model Inputfile Outputfile  LineStart LineEnd \n";
		ShowInputError();
		return EXIT_FAILURE;
	}


    auto Model=ModelID::getModel(argv[1]);
    if(Model==ModelID::ModelIDs::NotSet) {
		std::cerr << "Your Model parameter does not match with the implemented Models." << std::endl;
		ShowInputError();
		return EXIT_FAILURE;
	}

	double LineStart,LineEnd;
	char* in_file;char* out_file;


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
			std::cerr << "Input file not found " << std::endl;
			return EXIT_FAILURE;
	}
	std::ofstream outfile(out_file);
	if(!outfile.good())
	{
		std::cerr << "Can not create file " << out_file << std::endl;
		return EXIT_FAILURE;
	}
	std::string linestr;


    std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer = ModelID::FChoose(Model);


    size_t nPar,nParCT;
    nPar = modelPointer->get_nPar();
    nParCT = modelPointer->get_nParCT();

    size_t ndim = modelPointer->get_nVEV();
	std::vector<double> par(nPar);
	std::vector<double> parCT(nParCT);




	std::vector<double> sol,Check,Start,start;


	std::vector<double> Weinberg;
	while(true)
	{

		getline(infile,linestr);
		if(linecounter > LineEnd) break;
		if(linecounter == 1)
		  {
			modelPointer->setUseIndexCol(linestr);
            outfile << linestr;
            auto legendCT = modelPointer->addLegendCT();
            for(auto x: legendCT) outfile << sep << x;
            auto legendVEV = modelPointer->addLegendVEV();
            for(auto x: legendVEV) outfile << sep << x;
            outfile << sep << "v_NLO";
		    outfile << std::endl;
		  }
		if(linecounter >= LineStart and linecounter <= LineEnd and linecounter != 1)
		{
			std::pair<std::vector<double>,std::vector<double>> parameters = modelPointer->initModel(linestr);
			par=parameters.first;
			parCT = parameters.second;

			if(LineStart == LineEnd ) modelPointer->write();
			sol.clear();
			Check.clear();
			start.clear();
            for(size_t i=0;i<ndim;i++) start.push_back(modelPointer->get_vevTreeMin(i));
            sol = Minimizer::Minimize_gen_all(modelPointer,0,Check,start);


			std::vector<double> solPot,solSym;
            solPot=modelPointer->MinimizeOrderVEV(sol);
			double vev = modelPointer->EWSBVEV(solPot);

			outfile << linestr;
            for(size_t i=0;i<nParCT;i++) outfile << sep << parCT[i];
            for(size_t i=0;i<ndim;i++) outfile << sep << sol.at(i);
            outfile << sep << vev;
			outfile << std::endl;

			if(LineStart==LineEnd){
                auto dimensionnames = modelPointer->addLegendVEV();
                for(size_t i=0;i<ndim;i++){
					std::cout << dimensionnames.at(i) << " = " << sol.at(i) << " GeV" << std::endl;
				}
			}
		}
		linecounter++;
		if(infile.eof()) break;
	}

	outfile.close();

	return EXIT_SUCCESS;
}
catch(exception& e){
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
}



