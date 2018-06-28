/*
 * CalcCT.cpp
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
 * Calculates the Counterterms for mulitple points and adds them at the end of the line. The sequence is the
 * same as used in Potential::set_CT_Pot_Par
 */

#include "../models/IncludeAllModels.h"

#include <iostream>
using namespace std;
//#include "Minimizer.h"
int main(int argc, char *argv[]) try{

	if(!( argc == 6) )
	{
		std::cout << "./CalcCT Model Inputfile Outputfile LineStart LineEnd \n";
		return EXIT_FAILURE;
	}

	int Model;
	double LineStart,LineEnd;
	char* in_file;char* out_file;


	in_file = argv[2];
	out_file = argv[3];



	Model = atoi(argv[1]);
	LineStart = atoi(argv[4]);
	LineEnd = atoi(argv[5]);

	if(LineStart < 1)
	{
		std::cout << "Start line counting with 1" << std::endl;
		return EXIT_FAILURE;
	}
	if(LineStart > LineEnd)
	{
		std::cout << "LineEnd is smaller then LineStart " << std::endl;
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
	int Type;
	double tmp;

	std::unique_ptr<Class_Potential_Origin> modelPointer = FChoose(Model);
	int nPar,nParCT;
	nPar = modelPointer->nPar;
	nParCT = modelPointer->nParCT;
	std::vector<double> par(nPar);
	std::vector<double> parCT(nParCT);




	while(getline(infile,linestr))
	{


		if(linecounter > LineEnd) break;
		if(linecounter == 1)
		  {
		    outfile << linestr << "\t" << modelPointer->addLegendCT();
		    outfile << std::endl;
		  }

		if(linecounter >= LineStart and linecounter <= LineEnd and linecounter != 1)
		{
			modelPointer->resetbools();
			modelPointer->ReadAndSet(linestr,par);
//			std::cout << linestr << std::endl;
//			modelPointer->write();
			modelPointer->calc_CT(parCT);
			outfile << linestr;
			for(int i=0;i<nParCT;i++) outfile << "\t" << parCT[i];
			outfile << std::endl;
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
