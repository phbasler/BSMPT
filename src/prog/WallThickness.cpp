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
#include <fstream>

using namespace std;
using namespace BSMPT;





//#include "Minimizer.h"

int main(int argc, char *argv[]) try{

	if(!( argc == 6 or argc == 7) )
	{
		std::cerr << "./WallThickness Model Inputfile Outputfile  LineStart LineEnd \n";
		ShowInputError();
		return EXIT_FAILURE;
	}


    auto Model=ModelID::getModel(argv[1]);
    if(Model==ModelID::ModelIDs::NotSet) {
        std::cerr << "Your Model parameter does not match with the implemented Models." << std::endl;
        ShowInputError();
        return EXIT_FAILURE;
    }




	double LineStart,LineEnd;
	char* in_file;char* out_file;

	in_file = argv[2];
	out_file = argv[3];


	LineStart = atoi(argv[4]);
	LineEnd = atoi(argv[5]);

	bool TerminalOutput = false;
	if(argc == 7) {
		std::string s7 = argv[6];
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


    std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer = ModelID::FChoose(Model);

    size_t nPar,nParCT;
    nPar = modelPointer->get_nPar();
    nParCT = modelPointer->get_nParCT();

    size_t ndim = modelPointer->get_nVEV();


	std::vector<double> par(nPar);
	std::vector<double> parCT(nParCT);


	std::vector<double> Weinberg;

	while(getline(infile,linestr))
	{
		if(linecounter > LineEnd) break;

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
		if(linecounter >= LineStart and linecounter <= LineEnd and linecounter != 1)
		{
			if(TerminalOutput)
			{
				std::cout << "Currently at line " << linecounter << std::endl;
			}
			std::pair<std::vector<double>,std::vector<double>> parameters = modelPointer->initModel(linestr);
			par=parameters.first;
			parCT = parameters.second;

			if(LineStart == LineEnd ) modelPointer->write();

            auto EWPT = Minimizer::PTFinder_gen_all(modelPointer,0,300);
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
                    for(size_t i=0;i<modelPointer->get_nVEV();i++) basepoint.push_back(VEVSymmetric.at(i)
							+ 0.5*(vcritical.at(i) - VEVSymmetric.at(i)));

                    auto SolMinPlane = Minimizer::MinimizePlane(basepoint,VEVSymmetric,vcritical,Model,par,parCT,EWPT.Tc);
                    auto MaximumPlane = SolMinPlane.Minimum;


                    LW1D = Wall::calculate_wall_thickness_1D(modelPointer,EWPT.Tc,vcritical,VEVSymmetric);



                    LWplane = Wall::calculate_wall_thickness_plane(modelPointer,EWPT.Tc,vcritical,VEVSymmetric);



					outfile << linestr;
                    for(size_t i=0;i<nParCT;i++) {
                        outfile << sep << parCT[i];
                    }
                    outfile << sep << EWPT.Tc << sep << EWPT.vc;
                    outfile << sep << EWPT.vc / EWPT.Tc;
                    for(size_t i=0;i<ndim;i++) outfile << sep << EWPT.EWMinimum.at(i);
                    outfile << sep << LW1D;
                    outfile << sep << LWplane;
					outfile << std::endl;
				}
			}

			if(LineStart == LineEnd) {
                auto dimensionnames = modelPointer->addLegendTemp();
                std::cout << "Succeded ? " << EWPT.StatusFlag << sep <<" (1 = Success , -1 = v/T reached a value below " << C_PT << " during the calculation) \n";
                if(EWPT.StatusFlag==1)
                {
                    std::cout << std::scientific;
                    std::cout << dimensionnames.at(1) << " = " << EWPT.vc << " GeV\n";
                    std::cout << dimensionnames.at(0) << " = " << EWPT.Tc << " GeV\n";
                    std::cout << "xi_c = " << dimensionnames.at(2)  << " = " << EWPT.vc/EWPT.Tc << std::endl;
                    for(size_t i=0;i<ndim; i++){
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
	if(TerminalOutput) std::cout << std::endl;
	outfile.close();

//	delete modelPointer;
	return EXIT_SUCCESS;





}

catch(exception& e){
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
}
