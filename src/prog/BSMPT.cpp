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






//#include "Minimizer.h"

int main(int argc, char *argv[]) try{

	if(!( argc == 6 or argc == 7) )
	{
		std::cout << "./BSMPT Model Inputfile Outputfile  LineStart LineEnd \n";
		std::cout << "The chosen Method is ";
		if(C_UseParwani) std::cout << "Parwani ";
		else std::cout << "Arnold Espinosa\n";
		std::cout << "The implemented models are \n"
				<< "0 : C2HDM\n"
				<< "1 : R2HDM\n"
				<< "2 : N2HDM"
				<< std::endl;
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

	bool TerminalOutput = false;
	if(argc == 7) {
		std::string s7 = argv[6];
		std::cout << s7 << std::endl;
		TerminalOutput = ("y" == s7);

	}

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


	std::vector<double> Weinberg;

	while(getline(infile,linestr))
	{
		if(linecounter > LineEnd) break;

		if(linecounter == 1)
		  {
		    outfile << linestr << "\t" << modelPointer->addLegendCT() << "\t";
		    outfile << modelPointer->addLegendTemp();
		    outfile << std::endl;
		  }
		if(linecounter >= LineStart and linecounter <= LineEnd and linecounter != 1)
		{
			if(TerminalOutput)
			{
				std::cout << "Currently at line " << linecounter << std::endl;
			}
			modelPointer->resetbools();
			modelPointer->ReadAndSet(linestr,par);

//
			modelPointer->calc_CT(parCT);

			modelPointer->set_CT_Pot_Par(parCT);
			if(LineStart == LineEnd ) modelPointer->write();

			/*std::vector<double> res;
			modelPointer->HiggsMassesSquared(res,modelPointer->vevTree,0,0);
			for(int i=0;i<modelPointer->NHiggs;i++) std::cout << std::sqrt(res.at(i)) << std::endl;
            */


			sol.clear();
			PTFinder_gen_all(Model,par,parCT,0,300,sol,3);
			if(LineStart == LineEnd) {
				std::string labels=modelPointer->addLegendTemp();
				std::string delimiter = "\t";
				std::vector<std::string> dimensionnames;
				size_t pos = 0;
				while((pos = labels.find(delimiter)) != std::string::npos){
					dimensionnames.push_back(labels.substr(0,pos));
					labels.erase(0,pos+delimiter.length());
				}
				dimensionnames.push_back(labels);
				if(dimensionnames.size() != ndim +3){
					std::cout << "The number of names in the function addLegendTemp does not match the number of vevs, going to default naming."
							<< "You should fix this as this will result in errors in your output file." << std::endl;
					std::cout << "Succeded ? " << sol.at(2) << "\t (1 = Success , -1 = v/T reached a value below " << C_PT << " during the calculation) \n";
					std::cout << "omega_c = " << sol.at(1) << " GeV\n";
					std::cout << "T_c = " << sol.at(0) << " GeV\n";
					std::cout << "xi_c = omega_c/T_c =  " << sol.at(1)/sol.at(0) << std::endl;
					for(int i=3;i<ndim+3 ;i++) {
						std::cout << "omega_" << i-2 << " = " << sol.at(i) << " GeV\n";}
				}
				else{
					std::cout << "Succeded ? " << sol.at(2) << "\t (1 = Success , -1 = v/T reached a value below " << C_PT << " during the calculation) \n";
					std::cout << dimensionnames.at(1) << " = " << sol.at(1) << " GeV\n";
					std::cout << dimensionnames.at(0) << " = " << sol.at(0) << " GeV\n";
					std::cout << "xi_c = " << dimensionnames.at(ndim+2)  << " = " << sol.at(1)/sol.at(0) << std::endl;
					for(int i=3;i<ndim + 3; i++){
						std::cout << dimensionnames.at(i-1) << " = " << sol.at(i) << " GeV\n";
					}
				}
			}
			if(sol.at(2) == 1)
			{
				if(C_PT*sol.at(0) < sol.at(1))
				{
					outfile << linestr;
					for(int i=0;i<nParCT;i++) outfile << "\t" << parCT[i];
					outfile << "\t" << sol.at(0) << "\t" << sol.at(1);
					for(int i=0;i<ndim;i++) outfile << "\t" << sol.at(i+3);
					outfile << "\t" << sol.at(1) / sol.at(0);
					outfile << std::endl;
				}
			}


		}
		linecounter++;
		if(infile.eof()) break;
	}
	if(TerminalOutput) std::cout << std::endl;
	outfile.close();

//	delete modelPointer;
	return EXIT_SUCCESS;





}

catch(exception& e){
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
}
