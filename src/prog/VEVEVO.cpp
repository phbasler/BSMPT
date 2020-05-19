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


auto getCLIArguments(int argc, char *argv[])
{
    struct ReturnType{
        BSMPT::ModelID::ModelIDs Model{};
        int Line{};
        std::string InputFile, OutputFile;
        double TemperatureStart{}, TemperatureStep{}, TemperatureEnd{};
    };

    std::vector<std::string> args;
    for(int i{1};i<argc;++i) args.push_back(argv[i]);

    if(argc < 8 or args.at(0) == "--help")
    {
        int SizeOfFirstColumn = std::string("--TemperatureStart=           ").size();
        std::cout << "VEVEVO calculates the evolution of the global minimum with rising temperature for a given parameter point" << std::endl
                  << "It is called either by " << std::endl
                  << "./VEVEVO Model Inputfile Outputfile Line TemperatureStart TemperatureStep TemperatureEnd" << std::endl
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
                  <<"The line in the input file with the given parameter point. Expects line 1 to be a legend." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--TemperatureStart="
                  <<"The starting temperature to calculate the global minimum." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--TemperatureStep="
                  <<"The stepsize for the temperature." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--TemperatureEnd="
                  <<"The last temperature to calculate the global minimum." << std::endl;
        ShowInputError();
    }

    if(args.size() > 0 and args.at(0)=="--help")
    {
        throw int{0};
    }
    else if(argc < 8)
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
            else if(StringStartsWith(el,"--line="))
            {
                res.Line = std::stoi(el.substr(std::string("--line=").size()));
            }
            else if(StringStartsWith(el,"--temperaturestart="))
            {
                res.TemperatureStart = std::stod(el.substr(std::string("--temperaturestart=").size()));
            }
            else if(StringStartsWith(el,"--temperaturestep="))
            {
                res.TemperatureStep = std::stod(el.substr(std::string("--temperaturestep=").size()));
            }
            else if(StringStartsWith(el,"--temperatureend="))
            {
                res.TemperatureEnd = std::stod(el.substr(std::string("--temperatureend=").size()));
            }
        }
    }
    else{
        res.Model = ModelID::getModel(args.at(0));
        res.InputFile = args.at(1);
        res.OutputFile = args.at(2);
        res.Line = std::stoi(args.at(3));
        res.TemperatureStart = std::stod(args.at(4));
        res.TemperatureStep = std::stod(args.at(5));
        res.TemperatureEnd = std::stod(args.at(6));
    }


    if(res.TemperatureStart < 0)
    {
        std::cout << "The starting value of your Temperature was negative. This was corrected to TemperatureStart = 0."
                <<std::endl;
        res.TemperatureStart = 0;
    }
    if(res.TemperatureEnd < res.TemperatureStart)
    {
        std::cout << "The value of Tempend was lower then the value of Tempstart. This was corrected by swapping them."
                <<std::endl;
        double tmp{res.TemperatureEnd};
        res.TemperatureEnd = res.TemperatureStart;
        res.TemperatureStart = tmp;

    }

    if(res.TemperatureStep == 0){
        std::cout << "The given stepsize is zero. This will cause an infinite loop. Therefore the stepsize has been set to 1." << std::endl;
        res.TemperatureStep = 1;
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


    std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer = ModelID::FChoose(args.Model);

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

    std::vector<double> Check;
    double vev{0.0};

	std::cout << std::scientific;
	std::cout << std::setprecision(16);
	outfile << std::setprecision(16);


    outfile << "T" << sep << "v";
    for(auto x: modelPointer->addLegendVEV()) outfile << sep << x;
    outfile << sep << "Veff(v,T)"
    		<< std::endl;

   for(double Temp = args.TemperatureStart; Temp<=args.TemperatureEnd; Temp+=args.TemperatureStep)
   {

	   start.clear();
       if(Temp==args.TemperatureStart)
	   {
           start = modelPointer->get_vevTreeMin();
	   }
	   else{
           start = sol;
	   }
	   sol.clear();
	   Check.clear();
	   solPot.clear();
       sol = Minimizer::Minimize_gen_all(modelPointer,Temp,Check,start);
       solPot=modelPointer->MinimizeOrderVEV(sol);
	   vev = modelPointer->EWSBVEV(solPot);

       outfile << Temp << sep;
       outfile << vev << sep;
       outfile << sol;
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
