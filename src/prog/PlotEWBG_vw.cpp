/*
 * PlotEWBG_vw.cpp
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
 * This program calculates the EWBG eta as a function of vw and varies vw over a given array.
 */


#include <bits/exception.h>                                 // for exception
#include <stdlib.h>                                         // for atof, EXI...
#include <algorithm>                                        // for copy, max
#include <memory>                                           // for shared_ptr
#include <string>                                           // for operator<<
#include <vector>                                           // for vector
#include <BSMPT/models/ClassPotentialOrigin.h>              // for Class_Pot...
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/baryo_calculation/CalculateEtaInterface.h>
#include <BSMPT/utility.h>
#include <iostream>
#include <fstream>
using namespace std;
using namespace BSMPT;

struct CLIoptions{
    BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
    int Line{};
    std::string InputFile, OutputFile,ConfigFile;
    bool TerminalOutput{false};
    double vw_min{},vw_max{},vw_Stepsize{};
    bool UseGSL { Minimizer::UseGSLDefault};
    bool UseCMAES {Minimizer::UseLibCMAESDefault};
    bool UseNLopt{Minimizer::UseNLoptDefault};
    int WhichMinimizer{Minimizer::WhichMinimizerDefault};
};

CLIoptions getCLIArguments(int argc, char *argv[]);

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
//Set up of BSMPT/Baryo Classes
    Baryo::CalculateEtaInterface EtaInterface(argv[8] /* = Config File */);
    std::shared_ptr<Class_Potential_Origin> modelPointer = ModelID::FChoose(args.Model);

    std::vector<double> start,solPot;

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
	bool found=false;


	while(true)
	{
	   if(infile.eof()) break;
	   std::getline(infile,linestr);
	   if(linecounter == 1){
		   modelPointer->setUseIndexCol(linestr);
		//outfile: LEGEND
            outfile<<linestr<<sep;
            outfile << "T_c_var"<<sep<<"omega_c_var"<<sep<<"vw_var"<<sep<<"LW_var";
			auto legend = EtaInterface.legend();
            for(const auto& x: legend) outfile<<sep<< x+"_var";
			outfile << std::endl;
	   }
       else if(linecounter == args.Line)
	   {
			modelPointer->initModel(linestr);
			modelPointer->FindSignSymmetries();
			found = true;
			break;
	   }
       else if(linecounter > args.Line) break;
	   linecounter++;
	   if(infile.eof()) break;
	}
	infile.close();
	if(!found) {std::cout << "Line not found !\n"; return -1;}

    if(args.TerminalOutput) modelPointer->write();
//CALL: BSMPT-->Phasetransition
    if(args.TerminalOutput) std::cout<<"PTFinder called..."<<std::endl;
    auto EWPT = Minimizer::PTFinder_gen_all(modelPointer,0,300);
//SFOEWPT FOUND 
    if(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS and C_PT*EWPT.Tc < EWPT.vc){
        if(args.TerminalOutput) std::cout<<"SFOEWPT found..."<<std::endl;
		std::vector<double> vcritical,vbarrier;
        vcritical = EWPT.EWMinimum;
        double TC = EWPT.Tc;
        double vc = EWPT.vc;
		std::vector<double> MinimumPlane;
//Find the minimum in the symmetric phase. For this minimise at T = Tc + 1
		std::vector<double> vevsymmetricSolution,checksym, startpoint;
        for(std::size_t i=0;i<modelPointer->get_nVEV();i++) startpoint.push_back(0.5*vcritical.at(i));
        vevsymmetricSolution=Minimizer::Minimize_gen_all(modelPointer,TC+1,checksym,startpoint);
		double vw = 0;
        if(args.TerminalOutput) std::cout<<"Currently calculating vw:"<<std::endl;
        for(vw=args.vw_min;vw<=args.vw_max;vw+=args.vw_Stepsize)
		{
            std::cout<<"\rvw = "<<vw<<"\n";
            auto eta = EtaInterface.CalcEta(vw,vcritical,vevsymmetricSolution,TC,modelPointer,args.WhichMinimizer);
            outfile<<linestr<<sep;
            outfile << TC << sep << vc << sep << vw << sep << EtaInterface.getLW();
            for(auto x:eta) outfile << sep << x;
			outfile << std::endl;
		}//END: vw loop
	}//END: SFOEWPT FOUND
	else{
		outfile<<-1<<-1<<-1<<-1<<-1<<std::endl;
	}//NO SFOEWPT
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


CLIoptions getCLIArguments(int argc, char *argv[])
{


    std::vector<std::string> args;
    for(int i{1};i<argc;++i) args.push_back(argv[i]);

    if(argc < 9 or args.at(0) == "--help")
    {
        int SizeOfFirstColumn = std::string("--TerminalOutput=           ").size();
        std::cout << "PlotEWBG_vw calculates the EWBG for varying wall velocity for a given parameter point." << std::endl
                  << "It is called either by " << std::endl
                  << "./PlotEWBG_vw Model Inputfile Outputfile Line vwMin vwStepsize vwMax EWBGConfigFile TerminalOutput(y/n)" << std::endl
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
                  <<"The line with the given parameter point. Expects line 1 to be a legend." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left << "--config="
                  << "The EWBG config file." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--TerminalOutput="
                  <<"y/n Turns on additional information in the terminal during the calculation." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left << "--vw_min="
                  << "The minimum wall velocity." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left << "--vw_max="
                  << "The maximum wall velocity." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left << "--vw_Stepsize="
                  << "The stepsize to increase the wall velocity." << std::endl;
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
    else if(argc < 9)
    {
        throw std::runtime_error("Too few arguments.");
    }


    CLIoptions res;
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
            else if(StringStartsWith(el,"--terminaloutput="))
            {
                res.TerminalOutput = el.substr(std::string("--lastline=").size()) == "y";
            }
            else if(StringStartsWith(el,"--vw_max="))
            {
                res.vw_min = std::stod(el.substr(std::string("--vw_min=").size()));
            }
            else if(StringStartsWith(el,"--vw_max="))
            {
                res.vw_max = std::stod(el.substr(std::string("--vw_max=").size()));
            }
            else if(StringStartsWith(el,"--vw_stepsize="))
            {
                res.vw_Stepsize = std::stod(el.substr(std::string("--vw_stepsize=").size()));
            }
            else if(StringStartsWith(el,"--config="))
            {
                res.ConfigFile = arg.substr(std::string("--config").size());
            }
            else if(StringStartsWith(el,"--usegsl="))
            {
                res.UseGSL = arg.substr(std::string("--usegsl=").size()) == "true";
                if(res.UseGSL and not Minimizer::UseGSLDefault)
                {
                    throw std::runtime_error("You set --UseGSL=true but GSL was not found during compilation.");
                }
            }
            else if(StringStartsWith(el,"--usecmaes="))
            {
                res.UseCMAES = arg.substr(std::string("--usecmaes=").size()) == "true";
                if(res.UseCMAES and not Minimizer::UseLibCMAESDefault)
                {
                    throw std::runtime_error("You set --UseCMAES=true but CMAES was not found during compilation.");
                }
            }
            else if(StringStartsWith(el,"--usenlopt="))
            {
                res.UseNLopt = arg.substr(std::string("--usenlopt=").size()) == "true";
                if(res.UseNLopt and not Minimizer::UseNLoptDefault)
                {
                    throw std::runtime_error("You set --UseNLopt=true but NLopt was not found during compilation.");
                }
            }
        }
        res.WhichMinimizer = Minimizer::CalcWhichMinimizer(res.UseGSL,res.UseCMAES,res.UseNLopt);
    }
    else{
        res.Model = ModelID::getModel(args.at(0));
        res.InputFile = args.at(1);
        res.OutputFile = args.at(2);
        res.Line = std::stoi(args.at(3));
        res.vw_min = std::stod(args.at(4));
        res.vw_Stepsize = std::stod(args.at(5));
        res.vw_max = std::stod(args.at(6));
        res.ConfigFile = args.at(7);
        if(argc == 10) {
            res.TerminalOutput = ("y" == std::string(argv[8]));
        }
    }

    if(res.WhichMinimizer == 0)
    {
        throw std::runtime_error("You disabled all minimizers. You need at least one.");
    }


    return res;
}
