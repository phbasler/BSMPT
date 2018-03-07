/*
 * CheckNLOVEV.cpp
 *
 *  Created on: Mar 17, 2017
 *      Author: basler
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


			outfile << linestr;
			for(int i=0;i<nParCT;i++) outfile << "\t" << parCT[i];
			for(int i=0;i<ndim;i++) outfile << "\t" << sol.at(i);
			outfile << "\t" << vev;
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



