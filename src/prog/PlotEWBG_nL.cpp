/*
 * PlotEWBG_nL.cpp
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
 * This program calculates the left-handed density in front of the bubble wall as a function of the distance z in both approaches.
 */


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
#include <BSMPT/baryo_calculation/Fluid_Type/tau_source.h>
#include <iostream>
#include <fstream>
using namespace std;
using namespace BSMPT;
using namespace Baryo;

struct CLIOptions{
    BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
    int Line{};
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

//Set up of BSMPT/Baryo Classes
    Baryo::CalculateEtaInterface EtaInterface(args.ConfigFile);
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
            outfile << linestr<<sep;
            outfile << "T_c_nLvar"<<sep<<"omega_c_nLvar"<<sep<<"vw_nLvar"<<sep<<"LW_nLvar" <<sep;
            outfile << "z" << sep <<"z/LW" << sep << "nL_VIA" << sep << "muL_FH" << sep << "nL_FH";
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
    auto EWPT = Minimizer::PTFinder_gen_all(modelPointer,0,300,args.WhichMinimizer);
//SFOEWPT FOUND
    if(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS and C_PT*EWPT.Tc < EWPT.vc){
        if(args.TerminalOutput) std::cout<<"SFOEWPT found..."<<std::endl;
        std::vector<double> vcritical,vbarrier;
        vcritical = EWPT.EWMinimum;
        double TC = EWPT.Tc;
        std::vector<double> MinimumPlane;
//Find the minimum in the symmetric phase. For this minimise at T = Tc + 1
        std::vector<double> vevsymmetricSolution,checksym, startpoint;
        for(std::size_t i=0;i<modelPointer->get_nVEV();i++) startpoint.push_back(0.5*vcritical.at(i));
        vevsymmetricSolution=Minimizer::Minimize_gen_all(modelPointer,TC+1,checksym,startpoint,args.WhichMinimizer);

/////////////////////////////////////////////////////////////////////////////////
        std::size_t nstep = 100;

        if(args.TerminalOutput) std::cout<<"Set up the numerics "<<std::endl;
        EtaInterface.setNumerics(args.vw , EWPT.EWMinimum , vevsymmetricSolution , TC , modelPointer,args.WhichMinimizer);//Set up parameter container for Baryo Calculation-->Calls container.init
        if(args.TerminalOutput) std::cout<<"Starting setting the class instances"<<std::endl;
        BSMPT::Baryo::tau_source C_tau;
        EtaInterface.set_transport_method(TransportMethod::tau);//setting to tau class
        bool botflag = true;
        auto class_GamM = EtaInterface.get_class_CalcGamM();
        auto class_ScP  = EtaInterface.get_class_Scp();
        auto class_kappa = EtaInterface.get_class_kappa();
        auto Integration_mubl{EtaInterface.getGSL_integration_mubl_container()};
        C_tau.set_class(botflag, Integration_mubl,class_GamM,class_ScP,class_kappa);

        auto tau_arr_nL = set_up_nL_grid(nstep,Integration_mubl,C_tau);
        auto FH_GSL = generate_mubl_spline(Integration_mubl, static_cast<int>(nstep));

///////
/// Outfile
///////
        for(std::size_t i=0;i<nstep;i++){
            outfile<<linestr<<sep;
            outfile<<TC<<sep<<EWPT.vc<<sep<<args.vw<<sep<<EtaInterface.getLW()<<sep;
            outfile<<tau_arr_nL.first.at(i)<<sep;
            outfile<<tau_arr_nL.first.at(i)/EtaInterface.getLW()<<sep;
            outfile<<tau_arr_nL.second.at(i)<<sep;
            outfile<<FH_GSL.spline(tau_arr_nL.first.at(i))<<sep;
            outfile<<FH_GSL.spline(tau_arr_nL.first.at(i))*std::pow(TC,2);
            outfile<<std::endl;
        }
/////////////////////////////////////////////////////////////////////////////////

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

CLIOptions::CLIOptions(int argc, char *argv[])
{


    std::vector<std::string> args;
    for(int i{1};i<argc;++i) args.push_back(argv[i]);

    if(argc < 7 or args.at(0) == "--help")
    {
        int SizeOfFirstColumn = std::string("--TerminalOutput=           ").size();
        std::cout << "Calculation of the left-handed chemical potentials or particle densities triggering the EW"
                  << " sphaleron transitions as a function of the wall distance z ." << std::endl
                  << "It is called either by " << std::endl
                  << "./PlotEWBG_nL Model Inputfile Outputfile Line vw EWBGConfigFile TerminalOutput(y/n)" << std::endl
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
                  << std::setw(SizeOfFirstColumn) << std::left << "--vw="
                  << "The wall velocity. Default value of 0.1." << std::endl;
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
        Line = std::stoi(args.at(3));
        vw = std::stod(args.at(4));
        ConfigFile = args.at(5);
        if(argc == 8) {
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
    if(Line < 1)
    {
        std::cerr << "Start line counting with 1" << std::endl;
        return false;
    }
    return true;
}

