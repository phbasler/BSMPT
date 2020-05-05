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
#include <iomanip>
using namespace std;
using namespace BSMPT;


auto getCLIArguments(int argc, char *argv[])
{
    struct ReturnType{
        BSMPT::ModelID::ModelIDs Model{};
        int FirstLine{}, LastLine{};
        std::string InputFile, OutputFile;
    };

    std::vector<std::string> args;
    for(int i{1};i<argc;++i) args.push_back(argv[i]);

    if(argc < 6 or args.at(0) == "--help")
    {
        int SizeOfFirstColumn = std::string("--TerminalOutput=           ").size();
        std::cout << "NLOVEV calculates the EW VEV at NLO" << std::endl
                  << "It is called either by " << std::endl
                  << "./NLOVEV model input output FirstLine LastLine" << std::endl
                  << "or with the following arguments" << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<< "--help"
                  << "Shows this menu" << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left << "--model="
                  << "The model you want to investigate"<<std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--input="
                  << "The input file in tsv format" << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--output="
                  << "The output file in tsv format" << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--FirstLine="
                  <<"The first line in the input file to calculate the NLO EW VEV. Expects line 1 to be a legend." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--LastLine="
                  <<"The last line in the input file to calculate the NLO EW VEV." << std::endl;
        ShowInputError();
    }

    if(args.size() > 0 and args.at(0)=="--help")
    {
        throw int{0};
    }
    else if(argc < 6)
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
                res.FirstLine = std::stoi(el.substr(std::string("--firstline=").size()));
            }
            else if(StringStartsWith(el,"--lastline="))
            {
                res.LastLine = std::stoi(el.substr(std::string("--lastline=").size()));
            }
        }
    }
    else{
        res.Model = ModelID::getModel(args.at(0));
        res.InputFile = args.at(1);
        res.OutputFile = args.at(2);
        res.FirstLine = std::stoi(args.at(3));
        res.LastLine = std::stoi(args.at(4));
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

    if(args.FirstLine < 1)
	{
		std::cerr << "Start line counting with 1" << std::endl;
		return EXIT_FAILURE;
	}
    if(args.FirstLine > args.LastLine)
	{
		std::cerr << "LineEnd is smaller then LineStart " << std::endl;
		return EXIT_FAILURE;
	}

	int linecounter = 1;
    std::ifstream infile(args.OutputFile);
	if(!infile.good()) {
			std::cerr << "Input file not found " << std::endl;
			return EXIT_FAILURE;
	}
    std::ofstream outfile(args.OutputFile);
	if(!outfile.good())
	{
        std::cerr << "Can not create file " << args.OutputFile << std::endl;
		return EXIT_FAILURE;
	}
	std::string linestr;


    std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer = ModelID::FChoose(args.Model);


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
        if(linecounter > args.LastLine) break;
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
        if(linecounter >= args.FirstLine and linecounter <= args.LastLine and linecounter != 1)
		{
			std::pair<std::vector<double>,std::vector<double>> parameters = modelPointer->initModel(linestr);
			par=parameters.first;
			parCT = parameters.second;

            if(args.FirstLine == args.LastLine ) modelPointer->write();
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

            if(args.FirstLine==args.LastLine){
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
catch(int)
{
    return EXIT_SUCCESS;
}
catch(exception& e){
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
}



