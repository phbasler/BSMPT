/*
 * CreateMuGrid.cpp
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
 * Calculates the electroweak phase transition for a given Inputfile for a given subset of lines in the file
 *  and adds it at the end of the line in the format
 * T_c v_c all single vevs. One parameter point per line.
 *
 */

#include <bits/exception.h>                               // for exception
#include <stdlib.h>                                       // for EXIT_FAILURE
#include <algorithm>                                      // for copy, max
#include <memory>                                         // for shared_ptr
#include <string>                                         // for getline
#include <utility>                                        // for pair
#include <vector>                                         // for vector
#include "BSMPT/models/ClassPotentialOrigin.h"            // for Class_Poten...
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/baryo_calculation/transport_equations.h>
#include <BSMPT/utility.h>
#include <iostream>
#include <fstream>
using namespace std;
using namespace BSMPT;


int main(int argc, char *argv[]) try{

    const double vw = 0.1;

	if(!(argc == 5))
	{
	std::cerr << "./CreateMuGrid Model Inputfile Outputfile Line \n";
	ShowInputError();
	return EXIT_FAILURE;
	}
	char* in_file; char* out_file;
	in_file = argv[2];
	out_file = argv[3];
    int LineNumb;

    auto Model=ModelID::getModel(argv[1]);
    if(Model==ModelID::ModelIDs::NotSet) {
        std::cerr << "Your Model parameter does not match with the implemented Models." << std::endl;
        ShowInputError();
        return EXIT_FAILURE;
    }



	LineNumb = atoi(argv[4]);

	if(LineNumb < 1)
	{
	std::cerr << "Start line counting with 1" << std::endl;
	return EXIT_FAILURE;
	}





	std::vector<double> sol,start,solPot;
	std::vector<double> Weinberg,parCTVec;


    std::shared_ptr<Class_Potential_Origin> modelPointer = ModelID::FChoose(Model);

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








	std::vector<double> parStart,parEnd;
    parStart = std::vector<double>(modelPointer->get_NHiggs(),0);

    auto EWPT = Minimizer::PTFinder_gen_all(modelPointer,0,300);

	// find the minimum in the symmetric phase. For this minimise at T = Tc + 1
    std::vector<double> vevsymmetricSolution,checksym, startpoint;
    for(size_t i=0;i<modelPointer->get_nVEV();i++) startpoint.push_back(0.5*EWPT.EWMinimum.at(i));
    vevsymmetricSolution=Minimizer::Minimize_gen_all(modelPointer,EWPT.Tc+1,checksym,startpoint);


    double absvevsymmetricSolution = 0;
    for(auto x:vevsymmetricSolution) absvevsymmetricSolution+=std::pow(x,2);

    if(absvevsymmetricSolution != 0){
		std::cout << "Check for sign of non vanishing components of symmetric vacuum : " << std::endl;
	}



    struct Baryo::GSL_integration_mubl p;
    p.init(vw,EWPT.EWMinimum,vevsymmetricSolution,EWPT.Tc,modelPointer);




	std::cout << "vw = " << vw << std::endl;
    std::cout << "LW = " << p.getLW()*EWPT.Tc << "/TC" << std::endl;
    std::cout << "T_C = " << EWPT.Tc << std::endl;
    for(size_t i=0;i<modelPointer->get_NHiggs();i++) std::cout << "v_" << i << " = " << EWPT.EWMinimum.at(i) << std::endl;

	// double res = K1_fermion_interp(2.5,3.7);
	// std::cout << "res = " << res << std::endl;


    std::size_t nstep = 1000;
    double zmax = p.getZMAX();
	double zmin = 0;
	double stepsize = (zmax-zmin)/nstep;
	outfile << "z\tmu_{B_L}" << std::endl;
    for(size_t i=0;i<=nstep;i++){
		double z= zmin + stepsize*i;
        outfile << z << sep << Baryo::mubl_func(z,&p) << std::endl;

	}




	outfile.close();



	return EXIT_SUCCESS;





}

catch(exception& e){
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
}
