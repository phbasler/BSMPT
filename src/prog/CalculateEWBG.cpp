/*
 * CalculateEWBG.cpp
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
 * Calculates the electroweak baryogenesis for a given Inputfile for a given subset of lines in the file
 *  and adds it at the end of the line in the format
 * T_c v_c all single vevs. One parameter point per line.
 *
 */

#include <stdlib.h>                                         // for atoi, std::size_t
#include <algorithm>                                        // for max, copy
#include <memory>                                           // for shared_ptr
#include <string>                                           // for string
#include <vector>                                           // for vector
#include <BSMPT/baryo_calculation/transport_equations.h>    // for GSL_integ...
#include <BSMPT/models/ClassPotentialOrigin.h>              // for Class_Pot...
#include <iostream>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/baryo_calculation/CalculateEtaInterface.h>
#include <BSMPT/utility.h>
#include <fstream>

using namespace std;
using namespace BSMPT;

struct CLIOptions{
    BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
    int FirstLine{}, LastLine{};
    std::string InputFile, OutputFile,ConfigFile;
    bool TerminalOutput{false};
    double vw{0.1};
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


//Init: Interface Class for the different transport methods 
    Baryo::CalculateEtaInterface EtaInterface(args.ConfigFile);


	int linecounter = 1;
    std::ifstream infile(args.InputFile);
	if(!infile.good()) 
	{
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
    std::shared_ptr<Class_Potential_Origin> modelPointer = ModelID::FChoose(args.Model);// Declare the model pointer with the necessary parameters
    std::vector<std::string> etaLegend = EtaInterface.legend();// Declare the vector for the PTFinder algorithm
//Begin: Input Read
	while(getline(infile,linestr))
	{
        if(linecounter > args.LastLine) break;
		else if(linecounter == 1)
		{
			// Write legend
			modelPointer->setUseIndexCol(linestr);
            outfile << linestr;
            for(const auto& x: modelPointer->addLegendCT()) outfile << sep << x+"_EWBG";
            for(const auto& x: modelPointer->addLegendTemp()) outfile << sep << x+"_EWBG";
            outfile << sep << "vw";
            outfile << sep << "L_W";
            outfile << sep << "top_sym_phase";
            outfile << sep << "top_brk_phase";
            outfile << sep << "bot_sym_phase";
            outfile << sep << "bot_brk_phase";
            outfile << sep << "tau_sym_phase";
            outfile << sep << "tau_brk_phase";
            outfile << sep << etaLegend;
		    outfile << std::endl;
		}
        else if(linecounter >= args.FirstLine and linecounter <= args.LastLine and linecounter != 1)
		{
            if(args.TerminalOutput)
			{
				std::cout << "Currently at line " << linecounter << std::endl;
			}
//Begin: Parameter Set Up for BSMPT
            auto parameters = modelPointer->initModel(linestr);
			modelPointer->FindSignSymmetries();
            if(args.FirstLine == args.LastLine ) {
                modelPointer->write();
                std::cout << "vw = " << args.vw << std::endl;
            }
            if(args.TerminalOutput) std::cout<<"Calling PTFinder"<<std::endl;

//Call: BSMPT 
            auto EWPT = Minimizer::PTFinder_gen_all(modelPointer,0,300,args.WhichMinimizer);
//Define parameters for eta
			std::vector<double> eta;
            if(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS and C_PT*EWPT.Tc < EWPT.vc){
                if(args.TerminalOutput) std::cout<<"SFOEWPT found..."<<std::endl;
			//Find the minimum in the symmetric phase. For this minimise at T = Tc + 1
				std::vector<double> vevsymmetricSolution,checksym, startpoint;
                for(const auto& el: EWPT.EWMinimum) startpoint.push_back(0.5*el);
                vevsymmetricSolution=Minimizer::Minimize_gen_all(modelPointer,EWPT.Tc+1,checksym,startpoint,args.WhichMinimizer);
//Call: Calculation of eta in the different implemented approaches 
                if(args.TerminalOutput) std::cout<<"Calling CalcEta..."<<std::endl;
                eta = EtaInterface.CalcEta(args.vw,EWPT.EWMinimum,vevsymmetricSolution,EWPT.Tc,modelPointer,args.WhichMinimizer);
//Outfile
				outfile << linestr;
                outfile << sep << parameters.second;
                outfile << sep << EWPT.Tc << sep << EWPT.vc;
                outfile << sep << EWPT.vc / EWPT.Tc;
                outfile << sep << EWPT.EWMinimum;
                outfile << sep << args.vw;
                outfile << sep << EtaInterface.getLW();
                outfile << sep << EtaInterface.getSymmetricCPViolatingPhase_top();
                outfile << sep << EtaInterface.getBrokenCPViolatingPhase_top();
                outfile << sep << EtaInterface.getSymmetricCPViolatingPhase_bot();
                outfile << sep << EtaInterface.getBrokenCPViolatingPhase_bot();
                outfile << sep << EtaInterface.getSymmetricCPViolatingPhase_tau();
                outfile << sep << EtaInterface.getBrokenCPViolatingPhase_tau();
                outfile << sep << eta;
				outfile << std::endl;
			}//END: SFOEWPT found
			else{ // No SFOEWPT provided
				outfile << linestr;
                outfile << sep <<parameters.second;
                outfile << sep << EWPT.Tc << sep << EWPT.vc;
                outfile << sep << EWPT.EWMinimum;
                outfile << sep << EWPT.vc / EWPT.Tc;
                outfile << sep << args.vw;
                outfile << sep << -1; //LW
                outfile << sep << -50;//top sym CP phase
                outfile << sep << -50;//top brk CP phase
                outfile << sep << -50;//bot sym CP phase
                outfile << sep << -50;//bot brk CP phase
                outfile << sep << -50;//tau sym CP phase
                outfile << sep << -50;//tau brk CP phase
                for(std::size_t i=0;i<etaLegend.size();i++) outfile << sep << 0;
				outfile << std::endl;
			}//END: No SFOEWPT

            if(args.FirstLine == args.LastLine) {
                auto dimensionnames = modelPointer->addLegendTemp();
                std::cout << "Succeded ? " << static_cast<int>(EWPT.StatusFlag)
                          << sep << " (1 = Success , -1 = v/T reached a value below " << C_PT << " during the calculation) \n";
                if(EWPT.StatusFlag==Minimizer::MinimizerStatus::SUCCESS)
                {
                    std::cout << std::scientific;
                    std::cout << dimensionnames.at(1) << " = " << EWPT.vc << " GeV\n";
                    std::cout << dimensionnames.at(0) << " = " << EWPT.Tc << " GeV\n";
                    std::cout << "xi_c = " << dimensionnames.at(2)  << " = " << EWPT.vc/EWPT.Tc << std::endl;
                    for(std::size_t i=0;i<modelPointer->get_nVEV(); i++){
                        std::cout << dimensionnames.at(i+3) << " = " << EWPT.EWMinimum.at(i) << " GeV\n";
                    }
                    std::cout << "The Wall thickness is given by L_W  = " << EtaInterface.getLW()
                              << "GeV^-2\n"
                              << "L_W * T = " << EtaInterface.getLW() * EWPT.Tc << "\n";
                    for(std::size_t i=0;i<etaLegend.size();i++) std::cout << etaLegend.at(i) << " = " << eta.at(i) << std::endl;
                }

			}//END: LineStart == LineEnd
		}//END: Valid Line
		linecounter++;
		if(infile.eof()) break;
	}//END: Input Read
    if(args.TerminalOutput) std::cout << std::endl;
//Closing & Free
	outfile.close();
	return EXIT_SUCCESS;
}//END: Try
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

    if(argc < 7 or args.at(0) == "--help")
    {
        int SizeOfFirstColumn = std::string("--TerminalOutput=           ").size();
        std::cout << "CalculateEWBG calculates the strength of the electroweak baryogenesis" << std::endl
                  << "It is called either by " << std::endl
                  << "./CalculateEWBG model input output FirstLine LastLine ConfigFile" << std::endl
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
                  << std::setw(SizeOfFirstColumn) << std::left << "--config="
                  << "The EWBG config file." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--TerminalOutput="
                  <<"y/n Turns on additional information in the terminal during the calculation." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left << "--vw="
                  << "Wall velocity for the EWBG calculation. Default value of 0.1." << std::endl;
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
    else if(argc < 7)
    {
        throw std::runtime_error("Too few arguments.");
    }

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
                TerminalOutput = el.substr(std::string("--terminaloutput=").size()) == "y";
            }
            else if(StringStartsWith(el,"--vw="))
            {
                vw = std::stod(el.substr(std::string("--vw=").size()));
            }
            else if(StringStartsWith(el,"--config="))
            {
                ConfigFile = arg.substr(std::string("--config=").size());
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
        ConfigFile = args.at(5);
        if(argc == 8) {
            std::string s7 = argv[6];
            TerminalOutput = ("y" == std::string(args.at(6)));
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
    if(vw <= 0 or vw > 1)
    {
        throw std::runtime_error("The wall velocity has to be between 0 and 1.");
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
