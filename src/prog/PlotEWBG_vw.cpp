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


int main(int argc, char *argv[]) try{
    if((argc != 9) and (argc !=10))
	{
        std::cerr << "./PlotEWBG_vw Model Inputfile Outputfile Line vwMin vwStepsize vwMax EWBGConfigFile TerminalOutput(y/n)\n";
		ShowInputError();
		return EXIT_FAILURE;
	}
	char* in_file; char* out_file;
	in_file = argv[2];
	out_file = argv[3];
	double LineNumb;

    auto Model=ModelID::getModel(argv[1]);
    if(Model==ModelID::ModelIDs::NotSet) {
        std::cerr << "Your Model parameter does not match with the implemented Models." << std::endl;
        ShowInputError();
        return EXIT_FAILURE;
    }

	LineNumb = atoi(argv[4]);
	double vwMin, vwMax, vwStepsize;
	vwMin = atof(argv[5]);
	vwStepsize = atof(argv[6]);
	vwMax = atof(argv[7]);
	if(LineNumb < 1)
	{
		std::cerr << "Start line counting with 1" << std::endl;
		return EXIT_FAILURE;
	}
    bool TerminalOutput = false;
    if(argc ==10){
        std::string s7 = argv[9];
        std::cout<<"Terminal Output:"<<s7<<std::endl;
        TerminalOutput = ("y" == s7);
    }
//Set up of BSMPT/Baryo Classes
    Baryo::CalculateEtaInterface EtaInterface(argv[8] /* = Config File */);
    std::shared_ptr<Class_Potential_Origin> modelPointer = ModelID::FChoose(Model);

    std::vector<double> start,solPot;

	std::ifstream infile(in_file);
	if(!infile.good()) {
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
	   else if(linecounter == LineNumb)
	   {
			modelPointer->initModel(linestr);
			modelPointer->FindSignSymmetries();
			found = true;
			break;
	   }
	   else if(linecounter > LineNumb) break;
	   linecounter++;
	   if(infile.eof()) break;
	}
	infile.close();
	if(!found) {std::cout << "Line not found !\n"; return -1;}

    if(TerminalOutput) modelPointer->write();
//CALL: BSMPT-->Phasetransition
    if(TerminalOutput) std::cout<<"PTFinder called..."<<std::endl;
////////////////////////////////
    int WhichMinimizer = 3;/////
////////////////////////////////
    auto EWPT = Minimizer::PTFinder_gen_all(modelPointer,0,300,WhichMinimizer);//Change to NLOpt
//SFOEWPT FOUND 
    if(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS and C_PT*EWPT.Tc < EWPT.vc){
        if(TerminalOutput) std::cout<<"SFOEWPT found..."<<std::endl;
		std::vector<double> vcritical,vbarrier;
        vcritical = EWPT.EWMinimum;
        double TC = EWPT.Tc;
        double vc = EWPT.vc;
		std::vector<double> MinimumPlane;
//Find the minimum in the symmetric phase. For this minimise at T = Tc + 1
		std::vector<double> vevsymmetricSolution,checksym, startpoint;
        for(size_t i=0;i<modelPointer->get_nVEV();i++) startpoint.push_back(0.5*vcritical.at(i));
        vevsymmetricSolution=Minimizer::Minimize_gen_all(modelPointer,TC+1,checksym,startpoint,WhichMinimizer);//TODO: Change to NLopt
		double vw = 0;
        if(TerminalOutput) std::cout<<"Currently calculating vw:"<<std::endl;
		for(vw=vwMin;vw<=vwMax;vw+=vwStepsize)
		{
            std::cout<<"\rvw = "<<vw<<"\n";
			auto eta = EtaInterface.CalcEta(vw,vcritical,vevsymmetricSolution,TC,modelPointer);
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
catch(exception& e){
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
}
