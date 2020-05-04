/*
 * BSMPT.cpp
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
#include <iomanip>
#include <fstream>

using namespace std;
using namespace BSMPT;

auto getCLIArguments(int argc, char *argv[])
{
    struct ReturnType{
        BSMPT::ModelID::ModelIDs Model{};
        int FirstLine{}, LastLine{};
        std::string InputFile, OutputFile;
        bool TerminalOutput{false};
    };

    std::vector<std::string> args;
    for(int i{1};i<argc;++i) args.push_back(argv[i]);

    if(argc < 6 or args.at(0) == "--help")
    {
        int SizeOfFirstColumn = std::string("--TerminalOutput=           ").size();
        std::cout << "BSMPT calculates the strength of the electroweak phase transition" << std::endl
                  << "It is called either by " << std::endl
                  << "./BSMPT model input output FirstLine LastLine" << std::endl
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
                  <<"The first line in the input file to calculate the EWPT. Expects line 1 to be a legend." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--LastLine="
                  <<"The last line in the input file to calculate the EWPT." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--TerminalOutput="
                  <<"y/n Turns on additional information in the terminal during the calculation." << std::endl;
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
            else if(StringStartsWith(el,"--terminaloutput="))
            {
                res.TerminalOutput = el.substr(std::string("--lastline=").size()) == "y";
            }
        }
    }
    else{
        res.Model = ModelID::getModel(args.at(0));
        res.InputFile = args.at(1);
        res.OutputFile = args.at(2);
        res.FirstLine = std::stoi(args.at(3));
        res.LastLine = std::stoi(args.at(4));
        if(argc == 7) {
            std::string s7 = argv[6];
            res.TerminalOutput = ("y" == s7);
        }
    }


    return res;
}

int main(int argc, char *argv[]) try{

	/**
	 * PrintErrorLines decides if parameter points with no valid EWPT (no NLO stability or T=300 vanishing VEV)
	 * are printed in the output file
	 */
	bool PrintErrorLines=true;

    const auto args = getCLIArguments(argc,argv);
    if(args.Model==ModelID::ModelIDs::NotSet) {

        std::cerr << "Your Model parameter does not match with the implemented Models." << std::endl;
        ShowInputError();
        return EXIT_FAILURE;
    }

    if(args.FirstLine < 1)
	{
		std::cout << "Start line counting with 1" << std::endl;
		return EXIT_FAILURE;
	}
    if(args.FirstLine > args.LastLine)
	{
		std::cout << "LineEnd is smaller then LineStart " << std::endl;
		return EXIT_FAILURE;
	}


	int linecounter = 1;
    std::ifstream infile(args.InputFile);
	if(!infile.good()) {
        std::cout << "Input file " << args.InputFile << " not found " << std::endl;
		return EXIT_FAILURE;
	}

    std::ofstream outfile(args.OutputFile);
	if(!outfile.good())
	{
        std::cout << "Can not create file " << args.OutputFile << std::endl;
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


	while(getline(infile,linestr))
	{
        if(linecounter > args.LastLine) break;

		if(linecounter == 1)
		  {
            outfile << linestr << sep << modelPointer->addLegendCT()
                    << sep << modelPointer->addLegendTemp() << std::endl;

		    modelPointer->setUseIndexCol(linestr);
		  }
        if(linecounter >= args.FirstLine and linecounter <= args.LastLine and linecounter != 1)
		{
            if(args.TerminalOutput)
			{
				std::cout << "Currently at line " << linecounter << std::endl;
			}
			std::pair<std::vector<double>,std::vector<double>> parameters = modelPointer->initModel(linestr);
			par=parameters.first;
			parCT = parameters.second;
            if(args.FirstLine == args.LastLine ) {
                 modelPointer->write();
			}

            auto EWPT = Minimizer::PTFinder_gen_all(modelPointer,0,300);
            std::vector<double> vevsymmetricSolution,checksym, startpoint;
            for(size_t i=0;i<modelPointer->get_nVEV();i++) startpoint.push_back(0.5*EWPT.EWMinimum.at(i));
            auto VEVsym = Minimizer::Minimize_gen_all(modelPointer,EWPT.Tc+1,checksym,startpoint);


            if(args.FirstLine == args.LastLine) {
                auto dimensionnames = modelPointer->addLegendTemp();
                std::cout << "Success ? " << EWPT.StatusFlag << sep << " (1 = Yes , -1 = No, v/T reached a value below " << C_PT << " during the calculation) \n";
                if(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS){
                    std::cout << dimensionnames.at(1) << " = " << EWPT.vc << " GeV\n";
                    std::cout << dimensionnames.at(0) << " = " << EWPT.Tc << " GeV\n";
                    std::cout << "xi_c = " << dimensionnames.at(2)  << " = " << EWPT.vc/EWPT.Tc << std::endl;
                    for(size_t i=0;i<ndim; i++){
                        std::cout << dimensionnames.at(i+3) << " = " << EWPT.EWMinimum.at(i) << " GeV\n";
                    }
                    std::cout<< "Symmetric VEV config"<<std::endl;
                    for(size_t i=0;i<ndim; i++){
                        std::cout << dimensionnames.at(i+3) << " = " << VEVsym.at(i) << " GeV\n";
                    }


                }
                else if(EWPT.Tc == 300){
                    std::cout << dimensionnames.at(1) << " != 0 GeV at T = 300 GeV." << std::endl;
                }
                else if(EWPT.Tc == 0){
                    std::cout << "This point is not vacuum stable." << std::endl;
                }
			}
			if(PrintErrorLines){
				outfile << linestr;
                for(auto x: parCT) outfile << sep << x;
                outfile << sep << EWPT.Tc << sep << EWPT.vc;
                if(EWPT.vc>C_PT*EWPT.Tc and EWPT.StatusFlag==Minimizer::MinimizerStatus::SUCCESS) outfile << sep << EWPT.vc/EWPT.Tc;
                else outfile << sep << static_cast<int>(EWPT.StatusFlag);
                for(auto x: EWPT.EWMinimum) outfile << sep << x;
				outfile << std::endl;
			}
            else if(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS)
			{
                if(C_PT* EWPT.Tc < EWPT.vc)
				{
                    outfile << linestr << sep << parCT;
                    outfile << sep << EWPT.Tc << sep << EWPT.vc;
                    outfile << sep << EWPT.vc / EWPT.Tc;
                    for(auto x: EWPT.EWMinimum) outfile << sep << x;
					outfile << std::endl;
				}
			}


		}
		linecounter++;
		if(infile.eof()) break;
	}
    if(args.TerminalOutput) std::cout << std::endl;
	outfile.close();

//	delete modelPointer;
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



