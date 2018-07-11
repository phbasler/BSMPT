/*
 * CheckNLOVEV.cpp
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



/**
 * @file
 * Calculates the VEV at T = 0 at NLO for a given Inputfile for a given subset of lines in the file and
 *  adds it at the end of the line . One parameter point per line.
 *
 */

#include "../models/IncludeAllModels.h"
#include "../minimizer/Minimizer.h"
#include <iostream>
using namespace std;



//#include "Minimizer.h"

int main(int argc, char *argv[]) try{

	if(!( argc == 6) )
	{
		std::cout << "./NLOVEV Model Inputfile Outputfile  LineStart LineEnd \n";
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


	std::unique_ptr<Class_Potential_Origin> modelPointer = FChoose(Model);

	int Type;
	double tmp;

	int nPar,nParCT;
	nPar = modelPointer->nPar;
	nParCT = modelPointer->nParCT;

	int ndim = modelPointer->nVEV;
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
		    outfile << linestr << "\t" << modelPointer->addLegendCT() << "\t";
		    outfile << modelPointer->addLegendVEV();
		    outfile << "\t" << "v_NLO";
		    outfile << std::endl;
		  }
		if(linecounter >= LineStart and linecounter <= LineEnd and linecounter != 1)
		{
			modelPointer->resetbools();
			modelPointer->ReadAndSet(linestr,par);
			modelPointer->calc_CT(parCT);
			modelPointer->set_CT_Pot_Par(parCT);
			if(LineStart == LineEnd ) modelPointer->write();
			sol.clear();
			Check.clear();
			start.clear();
			for(int i=0;i<ndim;i++) start.push_back(modelPointer->vevTreeMin.at(i));
			Minimize_gen_all(Model,par,parCT,0,sol,Check,start);


			std::vector<double> solPot,solSym;
			modelPointer->MinimizeOrderVEV(sol,solPot);
			double vev = modelPointer->EWSBVEV(solPot);

//			outfile << linecounter << "\t";
			outfile << linestr;
			for(int i=0;i<nParCT;i++) outfile << "\t" << parCT[i];
			for(int i=0;i<ndim;i++) outfile << "\t" << sol.at(i);
			outfile << "\t" << vev;
			outfile << std::endl;

			if(LineStart==LineEnd){
				std::string labels=modelPointer->addLegendVEV();
				std::string delimiter = "\t";
				std::vector<std::string> dimensionnames;
				size_t pos = 0;
				while((pos = labels.find(delimiter)) != std::string::npos){
					dimensionnames.push_back(labels.substr(0,pos));
					labels.erase(0,pos+delimiter.length());
				}
				dimensionnames.push_back(labels);
				for(int i=0;i<ndim;i++){
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



