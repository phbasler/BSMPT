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
#include <iomanip>
#include <fstream>
using namespace std;
using namespace BSMPT;


auto getCLIArguments(int argc, char *argv[])
{
    struct ReturnType{
        BSMPT::ModelID::ModelIDs Model{};
        int Line{};
        std::string InputFile, OutputFile;
        double vw{0.1};
    };

    std::vector<std::string> args;
    for(int i{1};i<argc;++i) args.push_back(argv[i]);

    if(argc < 5 or args.at(0) == "--help")
    {
        int SizeOfFirstColumn = std::string("--TerminalOutput=           ").size();
        std::cout << "CreateMuGrid calculates the mu_{BL} potential in front of the bubble wall" << std::endl
                  << "It is called either by " << std::endl
                  << "./CreateMuGrid model input output Line" << std::endl
                  << "or with the following arguments" << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<< "--help"
                  << "Shows this menu" << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left << "--model="
                  << "The model you want to investigate"<<std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--input="
                  << "The input file in tsv format" << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--output="
                  << "The output file in tsv format" << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--Line="
                  <<"The line in the input file to calculate the EWPT. Expects line 1 to be a legend." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left << "--vw="
                  << "Wall velocity for the EWBG calculation. Default value of 0.1." << std::endl;
        ShowInputError();
    }

    if(args.size() > 0 and args.at(0)=="--help")
    {
        throw int{0};
    }
    else if(argc < 5)
    {
        throw std::runtime_error("Too few arguments.");
    }


    ReturnType res;
    std::string prefix{"--"};
    bool UsePrefix = StringStartsWith(args.at(0),prefix);
    if(UsePrefix)
    {
        for(const auto& arg: args)
        {
            auto el = arg;
            std::transform(el.begin(), el.end(), el.begin(), ::tolower);
            if(StringStartsWith(el,"--model="))
            {
                res.Model = BSMPT::ModelID::getModel(el.substr(std::string("--model=").size()));
            }
            else if(StringStartsWith(el,"--input="))
            {
                res.InputFile = arg.substr(std::string("--input=").size());
            }
            else if(StringStartsWith(el,"--output="))
            {
                res.OutputFile = arg.substr(std::string("--output=").size());
            }
            else if(StringStartsWith(el,"--firstline="))
            {
                res.Line = std::stoi(el.substr(std::string("--line=").size()));
            }
            else if(StringStartsWith(el,"--vw="))
            {
                res.vw = std::stod(el.substr(std::string("--vw=").size()));
            }
        }
    }
    else{
        res.Model = ModelID::getModel(args.at(0));
        res.InputFile = args.at(1);
        res.OutputFile = args.at(2);
        res.Line = std::stoi(args.at(3));
    }

    return res;
}


int main(int argc, char *argv[]) try{
    const auto args = getCLIArguments(argc,argv);
    if(args.Model==ModelID::ModelIDs::NotSet) {
        std::cerr << "Your Model parameter does not match with the implemented Models." << std::endl;
        ShowInputError();
        return EXIT_FAILURE;
    }

    if(args.Line < 1)
	{
	std::cerr << "Start line counting with 1" << std::endl;
	return EXIT_FAILURE;
	}





	std::vector<double> sol,start,solPot;
	std::vector<double> Weinberg,parCTVec;


    std::shared_ptr<Class_Potential_Origin> modelPointer = ModelID::FChoose(args.Model);

    std::ifstream infile(args.InputFile);
	if(!infile.good()) {
		std::cout << "Input file not found " << std::endl;
		return EXIT_FAILURE;
	}


    std::ofstream outfile(args.OutputFile);
	if(!outfile.good())
	{
    std::cout << "Can not create file " << args.OutputFile << std::endl;
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
     else if(linecounter == args.Line)
	 {
		 std::pair<std::vector<double>,std::vector<double>> parameters = modelPointer->initModel(linestr);
		 par=parameters.first;
		 parCT = parameters.second;
		 found=true;
	 }

     else if(linecounter > args.Line) break;
	 linecounter++;
	 if(infile.eof()) break;
	}
	infile.close();
	if(!found) {std::cout << "Line not found !\n"; return -1;}




	std::vector<double> vTree;


    for(std::size_t k=0;k<dim;k++) vTree.push_back(modelPointer->get_vevTreeMin(k));

	std::vector<double> parStart,parEnd;
    parStart = std::vector<double>(modelPointer->get_NHiggs(),0);
    auto EWPT = Minimizer::PTFinder_gen_all(modelPointer,0,300);


	// find the minimum in the symmetric phase. For this minimise at T = Tc + 1
    std::vector<double> vevsymmetricSolution,checksym, startpoint;
    for(std::size_t i=0;i<modelPointer->get_nVEV();i++) startpoint.push_back(0.5*EWPT.EWMinimum.at(i));
    vevsymmetricSolution=Minimizer::Minimize_gen_all(modelPointer,EWPT.Tc+1,checksym,startpoint);


    double absvevsymmetricSolution = 0;
    for(auto x:vevsymmetricSolution) absvevsymmetricSolution+=std::pow(x,2);

    if(absvevsymmetricSolution != 0){
		std::cout << "Check for sign of non vanishing components of symmetric vacuum : " << std::endl;
	}



    struct Baryo::GSL_integration_mubl p;
    p.init(args.vw,EWPT.EWMinimum,vevsymmetricSolution,EWPT.Tc,modelPointer);

    std::cout << "vw = " << args.vw << std::endl;
    std::cout << "LW = " << p.getLW()*EWPT.Tc << "/TC" << std::endl;
    std::cout << "T_C = " << EWPT.Tc << std::endl;
    for(std::size_t i=0;i<modelPointer->get_NHiggs();i++) std::cout << "v_" << i << " = " << EWPT.EWMinimum.at(i) << std::endl;

	// double res = K1_fermion_interp(2.5,3.7);
	// std::cout << "res = " << res << std::endl;


    std::size_t nstep = 1000;
    double zmin = 0;
    double stepsize = (p.getZMAX()-zmin)/nstep;
	outfile << "z\tmu_{B_L}" << std::endl;
    for(std::size_t i=0;i<=nstep;i++){
		double z= zmin + stepsize*i;
        outfile << z << sep << Baryo::mubl_func(z,&p) << std::endl;
	}




	outfile.close();



	return EXIT_SUCCESS;





}
catch(int)
{
    return EXIT_SUCCESS;
}
catch(exception& e){
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
}
