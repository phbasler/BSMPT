/*
 * WallThickness.cpp
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

#include <bits/exception.h>                        // for exception
#include <stdlib.h>                                // for atoi, EXIT_FAILURE
#include <algorithm>                               // for copy, max
#include <memory>                                  // for shared_ptr, __shar...
#include <string>                                  // for string, operator<<
#include <utility>                                 // for pair
#include <vector>                                  // for vector
#include <BSMPT/models/ClassPotentialOrigin.h>     // for Class_Potential_Or...
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/minimizer/MinimizePlane.h>
#include <BSMPT/WallThickness/WallThicknessLib.h>
#include <BSMPT/utility.h>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace BSMPT;

struct CLIOptions{
    BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
    int FirstLine{}, LastLine{};
    std::string InputFile, OutputFile;
    bool TerminalOutput{false};
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

    std::ofstream outfile(args.OutputFile);
	if(!outfile.good())
	{
        std::cout << "Can not create file " << args.OutputFile << std::endl;
		return EXIT_FAILURE;
	}
	std::string linestr;


    std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer = ModelID::FChoose(args.Model);
	std::vector<double> Weinberg;

	while(getline(infile,linestr))
	{
        if(linecounter > args.LastLine) break;

		if(linecounter == 1)
		  {
			modelPointer->setUseIndexCol(linestr);
            outfile << linestr;
            for(auto x: modelPointer->addLegendCT()) outfile << sep << x;
            for(auto x: modelPointer->addLegendTemp()) outfile << sep << x;
            outfile << sep << "L_W (1D)";
            outfile << sep << "L_W (Plane)";
		    outfile << std::endl;
		  }
        if(linecounter >= args.FirstLine and linecounter <= args.LastLine and linecounter != 1)
		{
            if(args.TerminalOutput)
			{
				std::cout << "Currently at line " << linecounter << std::endl;
			}
			std::pair<std::vector<double>,std::vector<double>> parameters = modelPointer->initModel(linestr);

            if(args.FirstLine == args.LastLine ) modelPointer->write();

            auto EWPT = Minimizer::PTFinder_gen_all(modelPointer,0,300,args.WhichMinimizer);
            auto EWPTFlag = EWPT.StatusFlag;

			double Vb1D=0, LW1D=0;
			double LWplane=0;



            if(EWPTFlag == Minimizer::MinimizerStatus::SUCCESS)
			{
                if(C_PT*EWPT.Tc < EWPT.vc)
				{
					std::vector<double> vcritical,vbarrier;
                    vcritical = EWPT.EWMinimum;

                    std::vector<double> VEVSymmetric(modelPointer->get_nVEV());


                    std::vector<double> basepoint;
                    for(std::size_t i=0;i<modelPointer->get_nVEV();i++) basepoint.push_back(VEVSymmetric.at(i)
							+ 0.5*(vcritical.at(i) - VEVSymmetric.at(i)));

                    auto SolMinPlane = Minimizer::MinimizePlane(basepoint,VEVSymmetric,vcritical,modelPointer,EWPT.Tc,args.WhichMinimizer);
                    auto MaximumPlane = SolMinPlane.Minimum;


                    LW1D = Wall::calculate_wall_thickness_1D(modelPointer,EWPT.Tc,vcritical,VEVSymmetric);



                    LWplane = Wall::calculate_wall_thickness_plane(modelPointer,EWPT.Tc,vcritical,VEVSymmetric,args.WhichMinimizer);



					outfile << linestr;
                    outfile << sep << parameters.second;
                    outfile << sep << EWPT.Tc << sep << EWPT.vc;
                    outfile << sep << EWPT.vc / EWPT.Tc;
                    outfile << sep << EWPT.EWMinimum;
                    outfile << sep << LW1D;
                    outfile << sep << LWplane;
					outfile << std::endl;
				}
			}

            if(args.FirstLine == args.LastLine) {
                auto dimensionnames = modelPointer->addLegendTemp();
                std::cout << "Succeded ? " << static_cast<int>(EWPT.StatusFlag)
                          << sep <<" (1 = Success , -1 = v/T reached a value below " << C_PT << " during the calculation) \n";
                if(EWPT.StatusFlag==Minimizer::MinimizerStatus::SUCCESS)
                {
                    std::cout << std::scientific;
                    std::cout << dimensionnames.at(1) << " = " << EWPT.vc << " GeV\n";
                    std::cout << dimensionnames.at(0) << " = " << EWPT.Tc << " GeV\n";
                    std::cout << "xi_c = " << dimensionnames.at(2)  << " = " << EWPT.vc/EWPT.Tc << std::endl;
                    for(std::size_t i=0;i<modelPointer->get_nVEV(); i++){
                        std::cout << dimensionnames.at(i+3) << " = " << EWPT.EWMinimum.at(i) << " GeV\n";
                    }
                    std::cout << "The potential barrier in the 1D approximation is given by Vb = " << Vb1D << " GeV^4" << std::endl;
                    std::cout << "The wall thickness is then given as L_W = " << LW1D << " GeV^-1 " << std::endl;
                    std::cout << "The thick wall parameter is given by T_c L_W = " << EWPT.Tc * LW1D << std::endl;


                    std::cout << "Results of all dimensions : \n";
                    std::cout << "The wall thickness is then given as L_W = " << LWplane << " GeV^-1 "
                              << sep <<" = " << LWplane * 1.239842e-15 << " m "
                              << std::endl;
                    std::cout << "The thick wall parameter is given by T_c L_W = " << EWPT.Tc * LWplane << std::endl;

                }
			}
		}
		linecounter++;
		if(infile.eof()) break;
	}
    if(args.TerminalOutput) std::cout << std::endl;
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

CLIOptions::CLIOptions(int argc, char *argv[])
{
    std::vector<std::string> args;
    for(int i{1};i<argc;++i) args.push_back(argv[i]);

    if(argc < 6 or args.at(0) == "--help")
    {
        int SizeOfFirstColumn = std::string("--TerminalOutput=           ").size();
        std::cout << "WallThickness calculates the wall thickness" << std::endl
                  << "It is called either by " << std::endl
                  << "./WallThickness Model Inputfile Outputfile  LineStart LineEnd" << std::endl
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
    else if(argc < 6)
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
            else if(StringStartsWith(el,"--output="))
            {
                OutputFile = arg.substr(std::string("--output=").size());
            }
            else if(StringStartsWith(el,"--firstline="))
            {
                FirstLine = std::stoi(el.substr(std::string("--firstline=").size()));
            }
            else if(StringStartsWith(el,"--lastline="))
            {
                LastLine = std::stoi(el.substr(std::string("--lastline=").size()));
            }
            else if(StringStartsWith(el,"--terminaloutput="))
            {
                TerminalOutput = el.substr(std::string("--lastline=").size()) == "y";
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
        OutputFile = args.at(2);
        FirstLine = std::stoi(args.at(3));
        LastLine = std::stoi(args.at(4));
        if(argc == 7) {
            TerminalOutput = ("y" == std::string(argv[5]));
        }
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
    if(FirstLine < 1)
    {
        std::cout << "Start line counting with 1" << std::endl;
        return false;
    }
    if(FirstLine > LastLine)
    {
        std::cout << "Firstline is smaller then LastLine " << std::endl;
        return false;
    }
    return true;
}
