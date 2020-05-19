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

#include <bits/exception.h>                                 // for exception
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



int main(int argc, char *argv[]) try{


	// Setting the default value of 0.1 for the wall velocity
	//TODO: vw --> As additional input read parameter 
	const double vw = 0.1;


	if((argc != 7) and (argc != 8 ))
	{
		std::cerr << "./CalculateEWBG Model Inputfile Outputfile  LineStart LineEnd ConfigFile TerminalOutput\n";
		ShowInputError();
		return EXIT_FAILURE;
	}


    auto Model=ModelID::getModel(argv[1]);
    if(Model==ModelID::ModelIDs::NotSet) {
        std::cerr << "Your Model parameter does not match with the implemented Models." << std::endl;
        ShowInputError();
        return EXIT_FAILURE;
    }


//Init: Interface Class for the different transport methods 
    Baryo::CalculateEtaInterface EtaInterface(argv[6] /* = Config file */);
//Begin: Input Read parameters
	double LineStart,LineEnd;
	char* in_file;char* out_file;

	in_file = argv[2];
	out_file = argv[3];

	LineStart = atoi(argv[4]);
	LineEnd = atoi(argv[5]);

	bool TerminalOutput = false;
	if(argc == 8) {
		std::string s7 = argv[7];
		std::cout << s7 << std::endl;
		TerminalOutput = ("y" == s7);

	}

	if(LineStart < 1)
	{
		std::cout << "Start line counting with 1" << std::endl;
		return EXIT_FAILURE;
	}
	if(LineStart > LineEnd)
	{
		std::cout << "LineEnd is smaller then LineStart " << std::endl;
		return EXIT_FAILURE;
	}

	int linecounter = 1;
	std::ifstream infile(in_file);
	if(!infile.good()) 
	{
		std::cout << "Input file not found " << std::endl;
		return EXIT_FAILURE;
	}

	std::ofstream outfile(out_file);
	if(!outfile.good())
	{
		std::cout << "Can not create file " << out_file << std::endl;
		return EXIT_FAILURE;
	}
	std::string linestr;
    std::shared_ptr<Class_Potential_Origin> modelPointer = ModelID::FChoose(Model);// Declare the model pointer with the necessary parameters
    std::vector<std::string> etaLegend = EtaInterface.legend();// Declare the vector for the PTFinder algorithm
//Begin: Input Read
	while(getline(infile,linestr))
	{
		if(linecounter > LineEnd) break;
		else if(linecounter == 1)
		{
			// Write legend
			modelPointer->setUseIndexCol(linestr);
            outfile << linestr;
            for(auto x: modelPointer->addLegendCT()) outfile << sep << x+"_EWBG";
            for(auto x: modelPointer->addLegendTemp()) outfile << sep << x+"_EWBG";
            outfile << sep << "vw";
            outfile << sep << "L_W";
            outfile << sep << "top_sym_phase";
            outfile << sep << "top_brk_phase";
            outfile << sep << "bot_sym_phase";
            outfile << sep << "bot_brk_phase";
            outfile << sep << "tau_sym_phase";
            outfile << sep << "tau_brk_phase";
            for(auto x: etaLegend) outfile << sep << x;
		    outfile << std::endl;
		}
		else if(linecounter >= LineStart and linecounter <= LineEnd and linecounter != 1)
		{
			if(TerminalOutput)
			{
				std::cout << "Currently at line " << linecounter << std::endl;
			}
//Begin: Parameter Set Up for BSMPT
			modelPointer->initModel(linestr);
			modelPointer->FindSignSymmetries();
			if(LineStart == LineEnd ) modelPointer->write();
            if(TerminalOutput) std::cout<<"Calling PTFinder"<<std::endl;

//Call: BSMPT 
            auto EWPT = Minimizer::PTFinder_gen_all(modelPointer,0,300);
//Define parameters for eta
			std::vector<double> eta;
            if(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS and C_PT*EWPT.Tc < EWPT.vc){
                if(TerminalOutput) std::cout<<"SFOEWPT found..."<<std::endl;
				std::vector<double> vcritical,vbarrier;
                vcritical = EWPT.EWMinimum;
                double TC = EWPT.Tc;
				std::vector<double> MinimumPlane;
			//Find the minimum in the symmetric phase. For this minimise at T = Tc + 1
				std::vector<double> vevsymmetricSolution,checksym, startpoint;
                for(size_t i=0;i<modelPointer->get_nVEV();i++) startpoint.push_back(0.5*vcritical.at(i));
                vevsymmetricSolution=Minimizer::Minimize_gen_all(modelPointer,TC+1,checksym,startpoint);
//Call: Calculation of eta in the different implemented approaches 
                if(TerminalOutput) std::cout<<"Calling CalcEta..."<<std::endl;
				eta = EtaInterface.CalcEta(vw,vcritical,vevsymmetricSolution,TC,modelPointer);
//Outfile
				outfile << linestr;
                for(auto x:modelPointer->get_parCTStored()) outfile << sep << x;
                outfile << sep << EWPT.Tc << sep << EWPT.vc;
                outfile << sep << EWPT.vc / EWPT.Tc;
                for(auto x: EWPT.EWMinimum) outfile << sep << x;
                outfile << sep << vw;
                outfile << sep << EtaInterface.getLW();
                outfile << sep << EtaInterface.GSL_integration_mubl_container.getSymmetricCPViolatingPhase_top();
                outfile << sep << EtaInterface.GSL_integration_mubl_container.getBrokenCPViolatingPhase_top();
                outfile << sep << EtaInterface.GSL_integration_mubl_container.getSymmetricCPViolatingPhase_bot();
                outfile << sep << EtaInterface.GSL_integration_mubl_container.getBrokenCPViolatingPhase_bot();
                outfile << sep << EtaInterface.GSL_integration_mubl_container.getSymmetricCPViolatingPhase_tau();
                outfile << sep << EtaInterface.GSL_integration_mubl_container.getBrokenCPViolatingPhase_tau();
                for(auto x: eta) outfile << sep << x;
				outfile << std::endl;
			}//END: SFOEWPT found
			else{ // No SFOEWPT provided
				outfile << linestr;
                for(auto x:modelPointer->get_parCTStored()) outfile << sep << x;
                outfile << sep << EWPT.Tc << sep << EWPT.vc;
                for(auto x: EWPT.EWMinimum) outfile << sep << x;
                outfile << sep << EWPT.vc / EWPT.Tc;
                outfile << sep << vw;
                outfile << sep << -1; //LW
                outfile << sep << -50;//top sym CP phase
                outfile << sep << -50;//top brk CP phase
                outfile << sep << -50;//bot sym CP phase
                outfile << sep << -50;//bot brk CP phase
                outfile << sep << -50;//tau sym CP phase
                outfile << sep << -50;//tau brk CP phase
                for(size_t i=0;i<etaLegend.size();i++) outfile << sep << 0;
				outfile << std::endl;
			}//END: No SFOEWPT

			if(LineStart == LineEnd) {
                auto dimensionnames = modelPointer->addLegendTemp();
                std::cout << "Succeded ? " << static_cast<int>(EWPT.StatusFlag)
                          << sep << " (1 = Success , -1 = v/T reached a value below " << C_PT << " during the calculation) \n";
                if(EWPT.StatusFlag==Minimizer::MinimizerStatus::SUCCESS)
                {
                    std::cout << std::scientific;
                    std::cout << dimensionnames.at(1) << " = " << EWPT.vc << " GeV\n";
                    std::cout << dimensionnames.at(0) << " = " << EWPT.Tc << " GeV\n";
                    std::cout << "xi_c = " << dimensionnames.at(2)  << " = " << EWPT.vc/EWPT.Tc << std::endl;
                    for(size_t i=0;i<modelPointer->get_nVEV(); i++){
                        std::cout << dimensionnames.at(i+3) << " = " << EWPT.EWMinimum.at(i) << " GeV\n";
                    }
                    std::cout << "The Wall thickness is given by L_W  = " << EtaInterface.getLW()
                              << "GeV^-2\n"
                              << "L_W * T = " << EtaInterface.getLW() * EWPT.Tc << "\n";
                    for(size_t i=0;i<etaLegend.size();i++) std::cout << etaLegend.at(i) << " = " << eta.at(i) << std::endl;
                }

			}//END: LineStart == LineEnd
		}//END: Valid Line
		linecounter++;
		if(infile.eof()) break;
	}//END: Input Read
	if(TerminalOutput) std::cout << std::endl;
//Closing & Free
	outfile.close();
	return EXIT_SUCCESS;
}//END: Try
catch(exception& e){
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
}
