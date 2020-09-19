/*
 * Test.cpp
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

#include <BSMPT/models/IncludeAllModels.h>
#include <iostream>
#include <iomanip>
#include <bits/exception.h>                     // for exception
#include <stdlib.h>                             // for EXIT_FAILURE, atoi
#include <algorithm>                            // for copy
#include <memory>                               // for unique_ptr
#include <string>                               // for getline, string
#include <utility>                              // for pair
#include <vector>                               // for vector
#include <BSMPT/models/ClassPotentialOrigin.h>  // for Class_Potential_Origin
#include <BSMPT/utility.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <fstream>
using namespace std;
using namespace BSMPT;

struct CLIOptions{
    BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
    int Line{};
    std::string InputFile;
    bool UseGSL { Minimizer::UseGSLDefault};
    bool UseCMAES {Minimizer::UseLibCMAESDefault};
    bool UseNLopt{Minimizer::UseNLoptDefault};
    int WhichMinimizer{Minimizer::WhichMinimizerDefault};

    CLIOptions(int argc, char *argv[]);
    bool good() const;
};
int main(int argc, char *argv[]) try{

    const CLIOptions args(argc,argv);
    if(not args.good())
    {
        return EXIT_FAILURE;
    }

	int linecounter = 1;
    std::ifstream infile(args.InputFile);
	if(!infile.good()) {
		std::cout << "Input file not found " << std::endl;
		return EXIT_FAILURE;
	}


	std::string linestr;
    std::unique_ptr<Class_Potential_Origin> modelPointer = ModelID::FChoose(args.Model);


	while(getline(infile,linestr))
	{
        if(linecounter > args.Line) break;

		if(linecounter == 1)
		  {

		    modelPointer->setUseIndexCol(linestr);

		  }
        if(linecounter == args.Line and linecounter != 1)
		{
            modelPointer->initModel(linestr);
			modelPointer->write();
			std::vector<double> dummy;
			modelPointer->Debugging(dummy,dummy);
            modelPointer->CheckImplementation(args.WhichMinimizer);
		}
		linecounter++;
		if(infile.eof()) break;
	}
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

CLIOptions::CLIOptions(int argc, char *argv[])
{


    std::vector<std::string> args;
    for(int i{1};i<argc;++i) args.push_back(argv[i]);

    if(argc < 4 or args.at(0) == "--help")
    {
        int SizeOfFirstColumn = std::string("--TerminalOutput=           ").size();
        std::cout << "Test performs a serious of tests on the given model. Intended for testing new models." << std::endl
                  << "It is called either by " << std::endl
                  << "./Test model input Line" << std::endl
                  << "or with the following arguments" << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<< "--help"
                  << "Shows this menu" << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left << "--model="
                  << "The model you want to test"<<std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--input="
                  << "The input file in tsv format" << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--Line="
                  <<"The line in the input file with the parameter point used to check the model." << std::endl;
        std::string GSLhelp{"--UseGSL="};
        GSLhelp += Minimizer::UseGSLDefault?"true":"false";
        std::cout << std::setw(SizeOfFirstColumn) << std::left <<GSLhelp
                  << "Use the GSL library to minimize the effective potential" << std::endl;
        std::string CMAEShelp{"--UseCMAES="};
        CMAEShelp += Minimizer::UseLibCMAESDefault?"true":"false";
        std::cout << std::setw(SizeOfFirstColumn) << std::left <<CMAEShelp
                  << "Use the CMAES library to minimize the effective potential" << std::endl;
        std::string NLoptHelp{"--UseNLopt="};
        NLoptHelp += Minimizer::UseNLoptDefault?"true":"false";
        std::cout << std::setw(SizeOfFirstColumn) << std::left <<NLoptHelp
                  << "Use the NLopt library to minimize the effective potential" << std::endl;
        ShowInputError();
    }

    if(args.size() > 0 and args.at(0)=="--help")
    {
        throw int{0};
    }
    else if(argc < 4)
    {
        throw std::runtime_error("Too few arguments.");
    }

    const std::string prefix{"--"};
    bool UsePrefix = StringStartsWith(args.at(0),prefix);
    if(UsePrefix)
    {
        for(const auto& arg: args)
        {
            auto el = arg;
            std::transform(el.begin(), el.end(), el.begin(), ::tolower);
            if(StringStartsWith(el,"--model="))
            {
                Model = BSMPT::ModelID::getModel(el.substr(std::string("--model=").size()));
            }
            else if(StringStartsWith(el,"--input="))
            {
                InputFile = arg.substr(std::string("--input=").size());
            }
            else if(StringStartsWith(el,"--line="))
            {
                Line = std::stoi(el.substr(std::string("--line=").size()));
            }
            else if(StringStartsWith(el,"--usegsl="))
            {
                UseGSL = el.substr(std::string("--usegsl=").size()) == "true";
            }
            else if(StringStartsWith(el,"--usecmaes="))
            {
                UseCMAES = el.substr(std::string("--usecmaes=").size()) == "true";
            }
            else if(StringStartsWith(el,"--usenlopt="))
            {
                UseNLopt = el.substr(std::string("--usenlopt=").size()) == "true";
            }
        }
        WhichMinimizer = Minimizer::CalcWhichMinimizer(UseGSL,UseCMAES,UseNLopt);
    }
    else{
        Model = ModelID::getModel(args.at(0));
        InputFile = args.at(1);
        Line = std::stoi(args.at(2));
    }
}

bool CLIOptions::good() const
{
    if(UseGSL and not Minimizer::UseGSLDefault)
    {
        throw std::runtime_error("You set --UseGSL=true but GSL was not found during compilation.");
    }
    if(UseCMAES and not Minimizer::UseLibCMAESDefault)
    {
        throw std::runtime_error("You set --UseCMAES=true but CMAES was not found during compilation.");
    }
    if(UseNLopt and not Minimizer::UseNLoptDefault)
    {
        throw std::runtime_error("You set --UseNLopt=true but NLopt was not found during compilation.");
    }
    if(WhichMinimizer == 0)
    {
        throw std::runtime_error("You disabled all minimizers. You need at least one.");
    }
    if(Model==ModelID::ModelIDs::NotSet) {
        std::cerr << "Your Model parameter does not match with the implemented Models." << std::endl;
        ShowInputError();
        return false;
    }
    if(Line < 1)
    {
        std::cerr << "Start line counting with 1" << std::endl;
        return false;
    }
    return true;
}
