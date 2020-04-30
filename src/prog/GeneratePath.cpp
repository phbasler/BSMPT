/*
 * GeneratePath.cpp
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
#include <ext/alloc_traits.h>                   // for __alloc_traits<>::val...
#include <math.h>                               // for sqrt
#include <stdlib.h>                             // for size_t, atoi, EXIT_FA...
#include <algorithm>                            // for copy, max
#include <memory>                               // for shared_ptr, __shared_...
#include <string>                               // for string, operator<<
#include <utility>                              // for pair
#include <vector>                               // for vector
#include <BSMPT/models/ClassPotentialOrigin.h>  // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/minimizer/MinimizePlane.h>
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

//	Class_Potential_Origin * modelPointer;
//	Fchoose(modelPointer,Model);

    std::shared_ptr<Class_Potential_Origin> modelPointer = ModelID::FChoose(Model);


    size_t nPar,nParCT;
    nPar = modelPointer->get_nPar();
    nParCT = modelPointer->get_nParCT();

    size_t ndim = modelPointer->get_nVEV();


	std::vector<double> par(nPar);
	std::vector<double> parCT(nParCT);


	std::vector<double> Weinberg;

    auto dimensionnames = modelPointer->addLegendTemp();

	while(getline(infile,linestr))
	{
		if(linecounter > LineEnd) break;

		if(linecounter == 1)
		  {
			modelPointer->setUseIndexCol(linestr);
            for(auto x: modelPointer->addLegendCT()) outfile << sep << x;
            for(auto x: modelPointer->addLegendTemp()) outfile << sep << x;
            outfile << sep << "line_parameter";
            for(size_t i=3;i<dimensionnames.size();i++) outfile << sep << dimensionnames.at(i) << "_line";
		    outfile << "\tV_eff_line";
            for(size_t i=3;i<dimensionnames.size();i++) outfile << sep << dimensionnames.at(i) << "_min";
		    outfile << "\tV_eff_min";
		    outfile << "\tDistance_line_min";
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



            if(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS)
			{
                if(C_PT*EWPT.Tc < EWPT.vc)
				{
					std::vector<double> vcritical,vbarrier;
                    vcritical = EWPT.EWMinimum;

                    std::vector<double> VEVSymmetric(modelPointer->get_nVEV());
					for(size_t i=0;i<VEVSymmetric.size();i++) VEVSymmetric.at(i) = 0;


                    std::vector<double> basepoint,basepointPot;



					double DeltaT = 1e-2;
					for(int tcounter=0;tcounter<=100;tcounter++){
						basepoint.clear();
						basepointPot.clear();
						double lineparam = tcounter*DeltaT;
                        for(size_t i=0;i<modelPointer->get_nVEV();i++) {
							basepoint.push_back(VEVSymmetric.at(i)+ lineparam*(vcritical.at(i) - VEVSymmetric.at(i)));

						}

                        double Temp = EWPT.Tc;
                        auto MinPlaneResult = Minimizer::MinimizePlane(basepoint,VEVSymmetric,vcritical,Model,par,parCT,Temp);
                        double Vmin = MinPlaneResult.PotVal;
                        auto MinimumPlane = MinPlaneResult.Minimum;
                        basepointPot=modelPointer->MinimizeOrderVEV(basepoint);
						double Vline = modelPointer->VEff(basepointPot,Temp);
						double DistanceLineMinimum = 0;
                        for(size_t i=0;i<modelPointer->get_nVEV();i++){
							DistanceLineMinimum += std::pow(basepoint.at(i)-MinimumPlane.at(i),2);
						}
						DistanceLineMinimum = std::sqrt(DistanceLineMinimum);





						outfile <<  linestr;
                        for(auto x: parCT) outfile << sep << x;
                        outfile << sep << EWPT.Tc << sep << EWPT.vc;
                        outfile << sep << EWPT.vc / EWPT.Tc;
                        for(auto x: EWPT.EWMinimum) outfile << sep << x;
                        outfile << sep <<lineparam;
                        for(size_t i=0;i<basepoint.size();i++) outfile << sep << basepoint.at(i);
                        outfile << sep << Vline;
                        for(size_t i=0;i<MinimumPlane.size();i++) outfile << sep << MinimumPlane.at(i);
                        outfile << sep << Vmin;
                        outfile << sep << DistanceLineMinimum;
						outfile << std::endl;
					}







				}
			}

			if(LineStart == LineEnd) {

				if(dimensionnames.size() != ndim +3){
					std::cout << "The number of names in the function addLegendTemp does not match the number of vevs, going to default naming."
							<< "You should fix this as this will result in errors in your output file." << std::endl;
                    std::cout << "Succeded ? " << EWPT.StatusFlag << sep << " (1 = Success , -1 = v/T reached a value below " << C_PT << " during the calculation) \n";
                    if(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS){std::cout << "omega_c = " << EWPT.vc << " GeV\n";
                    std::cout << "T_c = " << EWPT.Tc << " GeV\n";
                    std::cout << "xi_c = omega_c/T_c =  " << EWPT.vc/EWPT.Tc << std::endl;
                    for(size_t i=0;i<ndim ;i++) {
                        std::cout << "omega_" << i+1 << " = " << EWPT.EWMinimum.at(i) << " GeV\n";}
					}
				}
				else{
                    std::cout << "Succeded ? " << EWPT.StatusFlag << sep <<" (1 = Success , -1 = v/T reached a value below " << C_PT << " during the calculation) \n";
                    if(EWPT.StatusFlag==Minimizer::MinimizerStatus::SUCCESS)
					{
						std::cout << std::scientific;
                        std::cout << dimensionnames.at(1) << " = " << EWPT.vc << " GeV\n";
                        std::cout << dimensionnames.at(0) << " = " << EWPT.Tc << " GeV\n";
                        std::cout << "xi_c = " << dimensionnames.at(2)  << " = " << EWPT.vc/EWPT.Tc << std::endl;
                        for(size_t i=0;i<ndim; i++){
                            std::cout << dimensionnames.at(i+3) << " = " << EWPT.EWMinimum.at(i) << " GeV\n";
						}

					}
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
