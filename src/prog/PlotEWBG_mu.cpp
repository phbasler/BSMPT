/*
 * PlotEWBG_mu.cpp
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
 * Calculates eta as a function of the renormalised scale mu. The renormalisation scale mu is varied from 1/2 to 1.5 C_vev0 in NumberOfStep steps.
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
#include <BSMPT/baryo_calculation/transport_equations.h>
#include <BSMPT/baryo_calculation/CalculateEtaInterface.h>
#include <BSMPT/utility.h>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace BSMPT;

struct CLIOptions{
    BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
    int Line{}, NumberOfSteps{};
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
    Baryo::CalculateEtaInterface EtaInterface(argv[6] /* = Config file */);
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
    std::vector<std::string> etaLegend = EtaInterface.legend();// Declare the vector for the PTFinder algorithm
    while(getline(infile,linestr))
    {
        if(linecounter > args.Line) break;

        if(linecounter == 1)
          {
            outfile << linestr<<sep << "mu_factor"<<sep<<"mu";
            for(const auto& x: modelPointer->addLegendTemp()) outfile << sep <<x+"_mu";
            outfile << sep << "BSMPT_StatusFlag";
            outfile << sep << "vw";
            outfile << sep << "L_W";
            outfile << sep << "top_sym_phase";
            outfile << sep << "top_brk_phase";
            outfile << sep << "bot_sym_phase";
            outfile << sep << "bot_brk_phase";
            outfile << sep << "tau_sym_phase";
            outfile << sep << "tau_brk_phase";
            for(const auto& x: etaLegend) outfile << sep << x+"_muvar";
            outfile << std::endl;

            modelPointer->setUseIndexCol(linestr);
          }
        if(args.Line==linecounter)
        {
            std::pair<std::vector<double>,std::vector<double>> parameters = modelPointer->initModel(linestr);
            if(args.TerminalOutput) modelPointer->write();
            if(args.TerminalOutput)
            {
                std::cout<<"Calculating EWPT in default settings with:\n mu = "<<modelPointer->get_scale()
                        <<std::endl;
            }
            if(args.TerminalOutput)std::cout<<"Start of mu variation"<<std::endl;
            for(int step=0;step<args.NumberOfSteps;step++){
                double mu_factor = 1/2. + (step/static_cast<double>(args.NumberOfSteps));
                if(args.TerminalOutput) std::cout<<"\r currently mu_factor = "<<mu_factor<<std::endl;
                auto VEVnames = modelPointer->addLegendTemp();
                auto CT_mu=modelPointer->resetScale(C_vev0*mu_factor);
                auto EWPT_mu = Minimizer::PTFinder_gen_all(modelPointer,0,300,args.WhichMinimizer);
                std::vector<double> startpoint;
                for(const auto& x : EWPT_mu.EWMinimum) startpoint.push_back(x/2.);
                if(EWPT_mu.StatusFlag== Minimizer::MinimizerStatus::SUCCESS ){
                    std::vector<double>checkmu;
                    auto VEV_mu_sym = Minimizer::Minimize_gen_all(modelPointer,EWPT_mu.Tc+1,checkmu,startpoint,args.WhichMinimizer);
                    auto VEV_mu_brk = EWPT_mu.EWMinimum;
                    auto eta_mu = EtaInterface.CalcEta(args.vw,EWPT_mu.EWMinimum,VEV_mu_sym,EWPT_mu.Tc,modelPointer,args.WhichMinimizer);

                    outfile << linestr;
                    outfile << sep << mu_factor <<sep<< mu_factor*C_vev0;
                    outfile << sep << EWPT_mu.Tc<<sep<<EWPT_mu.vc<<sep<<EWPT_mu.vc/EWPT_mu.Tc<<sep<<EWPT_mu.EWMinimum;
                    outfile << sep << static_cast<int>(EWPT_mu.StatusFlag);
                    outfile << sep << args.vw;
                    outfile << sep << EtaInterface.getLW();
                    outfile << sep << EtaInterface.getSymmetricCPViolatingPhase_top();
                    outfile << sep << EtaInterface.getBrokenCPViolatingPhase_top();
                    outfile << sep << EtaInterface.getSymmetricCPViolatingPhase_bot();
                    outfile << sep << EtaInterface.getBrokenCPViolatingPhase_bot();
                    outfile << sep << EtaInterface.getSymmetricCPViolatingPhase_tau();
                    outfile << sep << EtaInterface.getBrokenCPViolatingPhase_tau();
                    outfile << sep << eta_mu;

                    outfile << std::endl;


                }//EWPT@mu found
                else {
                    if(args.TerminalOutput) std::cout<<"\tNo SFOEWPT found for given scale"<<std::endl;
                    continue;
                }
        }//END: Mu Factor
        }//END:LineCounter
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

    if(argc < 7 or args.at(0) == "--help")
    {
        int SizeOfFirstColumn = std::string("--TerminalOutput=           ").size();
        std::cout << "EWBGRenormScale calculates the strength of the EWBG while varying the MSBar renormalisation scale" << std::endl
                  << "It is called either by " << std::endl
                  << "./EWBGRenormScale Model Inputfile Outputfile  Line NumberOfSteps Configfile" << std::endl
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
                  <<"The line in the input file with the parameter point. Expects line 1 to be a legend." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left << "--config="
                  << "The EWBG config file." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--TerminalOutput="
                  <<"y/n Turns on additional information in the terminal during the calculation." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left << "--vw="
                  << "Wall velocity for the EWBG calculation. Default value of 0.1." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left << "--NumberOfSteps="
                  << "Number of Steps to vary the scale between 0.5 and 1.5 times the original scale." << std::endl;
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
            else if(StringStartsWith(el,"--line="))
            {
                Line = std::stoi(el.substr(std::string("--line=").size()));
            }
            else if(StringStartsWith(el,"--numberofsteps="))
            {
                NumberOfSteps = std::stoi(el.substr(std::string("--numberofsteps=").size()));
            }
            else if(StringStartsWith(el,"--terminaloutput="))
            {
                TerminalOutput = el.substr(std::string("--lastline=").size()) == "y";
            }
            else if(StringStartsWith(el,"--vw="))
            {
                vw = std::stod(el.substr(std::string("--vw=").size()));
            }
            else if(StringStartsWith(el,"--config="))
            {
                ConfigFile = arg.substr(std::string("--config").size());
            }
            else if(StringStartsWith(el,"--usegsl="))
            {
                UseGSL = arg.substr(std::string("--usegsl=").size()) == "true";
            }
            else if(StringStartsWith(el,"--usecmaes="))
            {
                UseCMAES = arg.substr(std::string("--usecmaes=").size()) == "true";
            }
            else if(StringStartsWith(el,"--usenlopt="))
            {
                UseNLopt = arg.substr(std::string("--usenlopt=").size()) == "true";
            }
        }
        WhichMinimizer = Minimizer::CalcWhichMinimizer(UseGSL,UseCMAES,UseNLopt);
    }
    else{
        Model = ModelID::getModel(args.at(0));
        InputFile = args.at(1);
        OutputFile = args.at(2);
        Line = std::stoi(args.at(3));
        NumberOfSteps = std::stoi(args.at(4));
        ConfigFile = args.at(5);
        if(argc == 8) {
            TerminalOutput = ("y" == std::string(argv[6]));
        }
    }
}

bool CLIOptions::good() const
{
    if(NumberOfSteps == 0)
    {
        throw std::runtime_error("You have set the number of steps to zero.");
    }
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
    if(vw < 0 or vw > 1)
    {
        throw std::runtime_error("The wall velocity has to be between 0 and 1.");
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
