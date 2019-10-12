/*
 * BSMPT.cpp
 *
 *
 *      Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner

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
 * Calculates the electroweak phase transition for a given Inputfile for a given subset of lines in the file
 *  and adds it at the end of the line in the format
 * T_c v_c all single vevs. One parameter point per line.
 *
 */

#include "../models/IncludeAllModels.h"
#include "../minimizer/Minimizer.h"
#include <iostream>
using namespace std;



int main(int argc, char *argv[]) try{
	bool PrintErrorLines=true;
	int Model=-1;

	if(!( argc == 4) )
	{
		std::cerr << "./Test Model Inputfile Line \n";
		ShowInputError();
		return EXIT_FAILURE;
	}


	Model=getModel(argv[1]);
	// std::cout << "Model parameter in BSMPT = " << Model << std::endl;
	if(Model==-1) {
		std::cerr << "Your Model parameter does not match with the implemented Models." << std::endl;
		ShowInputError();
		return EXIT_FAILURE;
	}
	int  Line;
	char* in_file;

	in_file = argv[2];

	Line = atoi(argv[3]);



	if(Line < 1)
	{
		std::cout << "Start line counting with 1" << std::endl;
		return EXIT_FAILURE;
	}


	int linecounter = 1;
	std::ifstream infile(in_file);
	if(!infile.good()) {
		std::cout << "Input file not found " << std::endl;
		return EXIT_FAILURE;
	}


	std::string linestr;

//	Class_Potential_Origin * modelPointer;
//	Fchoose(modelPointer,Model);

	std::unique_ptr<Class_Potential_Origin> modelPointer = FChoose(Model);


	int Type;
	double tmp;

	int nPar,nParCT;
	nPar = modelPointer->nPar;
	nParCT = modelPointer->nParCT;

	int ndim = modelPointer->nVEV;


	std::vector<double> par(nPar);
	std::vector<double> parCT(nParCT);




	std::vector<double> sol;



	while(getline(infile,linestr))
	{
		if(linecounter > Line) break;

		if(linecounter == 1)
		  {

		    modelPointer->setUseIndexCol(linestr);

		  }
		if(linecounter == Line and linecounter != 1)
		{
			std::pair<std::vector<double>,std::vector<double>> parameters = modelPointer->initModel(linestr);
			par=parameters.first;
			parCT = parameters.second;

			modelPointer->write();
			std::vector<double> dummy;
			modelPointer->Debugging(dummy,dummy);
			modelPointer->CheckImplementation(par,parCT);




		}
		linecounter++;
		if(infile.eof()) break;
	}

//	delete modelPointer;
	return EXIT_SUCCESS;





}

catch(exception& e){
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
}
