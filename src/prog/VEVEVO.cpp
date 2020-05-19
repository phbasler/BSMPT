/*
 * VEVEVO.cpp
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
 * This program calculates the development of the VeVs with the Temperature.
 */


#include <bits/exception.h>                     // for exception
#include <math.h>                               // for sqrt, abs
#include <stdlib.h>                             // for atof, EXIT_FAILURE, atoi
#include <algorithm>                            // for copy, max
#include <iomanip>                              // for operator<<, setprecision
#include <memory>                               // for shared_ptr, __shared_...
#include <string>                               // for getline, operator<<
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

int main(int argc, char *argv[]) try{
	if(!(argc == 8))
	{
		std::cerr << "./VEVEVO Model Inputfile Outputfile Line Tempstart Tempstep Tempend \n";
		ShowInputError();
		return EXIT_FAILURE;
	}
	char* in_file; char* out_file;
	in_file = argv[2];
	out_file = argv[3];
	double LineNumb;
	double TempStart,TempEnd,TempStep;

    auto Model=ModelID::getModel(argv[1]);
    if(Model==ModelID::ModelIDs::NotSet) {
		std::cerr << "Your Model parameter does not match with the implemented Models." << std::endl;
		ShowInputError();
		return EXIT_FAILURE;
	}



	LineNumb = atoi(argv[4]);
	TempStart = atof(argv[5]);
	TempStep = atof(argv[6]);
	TempEnd = atof(argv[7]);

	if(LineNumb < 1)
	{
		std::cerr << "Start line counting with 1" << std::endl;
		return EXIT_FAILURE;
	}


	if(TempStart < 0)
	{
		std::cout << "The starting value of your Temperature was negative. This was corrected to TempStart = 0"
				<<std::endl;
		TempStart = 0;
	}
	if(TempEnd < TempStart)
	{
		std::cout << "The value of Tempend was lower then the value of Tempstart. This was corrected by swapping them"
				<<std::endl;
		double tmp=TempStart;
		TempStart = TempEnd;
		TempEnd=tmp;

	}

    if(TempStep == 0){
        std::cout << "The given stepsize is zero. This will cause an infinite loop. Therefore the stepsize has been set to 1." << std::endl;
        TempStep = 1;
    }



	std::vector<double> sol,start,solPot;
	std::vector<double> Weinberg,parCTVec;


    std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer = ModelID::FChoose(Model);

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
	int linecounter = 1;

    std::size_t nPar,nParCT;
    nPar = modelPointer->get_nPar();
    nParCT = modelPointer->get_nParCT();

    std::size_t dim = modelPointer->get_nVEV();
	std::vector<double> par(nPar);
	std::vector<double> parCT(nParCT);

	bool found=false;
	while(true)
	{
	   if(infile.eof()) break;
	   std::getline(infile,linestr);
	   if(linecounter == 1){
		   modelPointer->setUseIndexCol(linestr);
	   }
	   else if(linecounter == LineNumb)
	   {
		   std::pair<std::vector<double>,std::vector<double>> parameters = modelPointer->initModel(linestr);
		   par=parameters.first;
		   parCT = parameters.second;
           found=true;
	   }

	   else if(linecounter > LineNumb) break;
	   linecounter++;
	   if(infile.eof()) break;
	}
	infile.close();
	if(!found) {std::cout << "Line not found !\n"; return -1;}

	std::vector<double> vTree;


    for(size_t k=0;k<dim;k++) vTree.push_back(modelPointer->get_vevTreeMin(k));

    std::vector<double> Check;
	double vev=0;
	std::vector<double> NTempMass;

    std::vector<double> zero(dim,0);

	std::vector<double> sol1D;


	std::cout << std::scientific;
	std::cout << std::setprecision(16);
	outfile << std::setprecision(16);


    outfile << "T" << sep << "v";
    for(auto x: modelPointer->addLegendVEV()) outfile << sep << x;
    outfile << sep << "Veff(v,T)"
    		<< std::endl;

   for(double Temp = TempStart; Temp<=TempEnd; Temp+=TempStep)
   {

	   start.clear();
	   if(Temp==TempStart)
	   {
           for(size_t k=0;k<dim;k++) start.push_back(vTree.at(k));
	   }
	   else{
           for(size_t k=0;k<dim;k++) start.push_back(sol.at(k));
	   }
	   sol.clear();
	   Check.clear();
	   solPot.clear();
       sol = Minimizer::Minimize_gen_all(modelPointer,Temp,Check,start);
       solPot=modelPointer->MinimizeOrderVEV(sol);
	   vev = modelPointer->EWSBVEV(solPot);

       outfile << Temp << sep;
	   outfile << vev;
       for(auto x: sol) outfile << sep << x;
       outfile << sep << modelPointer->VEff(solPot,Temp,0);
	   outfile << std::endl;


   }
   outfile.close();

	return EXIT_SUCCESS;
}
catch(exception& e){
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
}
