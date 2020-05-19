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
#include <stdlib.h>                             // for std::size_t, atoi, EXIT_FA...
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
#include <iomanip>
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
        std::cout << "Generate Path calculates the tunnel path from the broken to the symmetric minimum" << std::endl
                  << "It is called either by " << std::endl
                  << "./GeneratePath model input output FirstLine LastLine" << std::endl
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
                  <<"The first line in the input file to calculate the Wallthickness. Expects line 1 to be a legend." << std::endl
                  << std::setw(SizeOfFirstColumn) << std::left<<"--LastLine="
                  <<"The last line in the input file to calculate the Wallthickness." << std::endl
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


    std::shared_ptr<Class_Potential_Origin> modelPointer = ModelID::FChoose(args.Model);

    std::size_t ndim = modelPointer->get_nVEV();


	std::vector<double> Weinberg;

    auto dimensionnames = modelPointer->addLegendTemp();

	while(getline(infile,linestr))
	{
        if(linecounter > args.LastLine) break;

		if(linecounter == 1)
		  {
			modelPointer->setUseIndexCol(linestr);
            for(auto x: modelPointer->addLegendCT()) outfile << sep << x;
            for(auto x: modelPointer->addLegendTemp()) outfile << sep << x;
            outfile << sep << "line_parameter";
            for(std::size_t i=3;i<dimensionnames.size();i++) outfile << sep << dimensionnames.at(i) << "_line";
		    outfile << "\tV_eff_line";
            for(std::size_t i=3;i<dimensionnames.size();i++) outfile << sep << dimensionnames.at(i) << "_min";
		    outfile << "\tV_eff_min";
		    outfile << "\tDistance_line_min";
		    outfile << std::endl;
		  }
        if(linecounter >= args.FirstLine and linecounter <= args.LastLine and linecounter != 1)
		{
            if(args.TerminalOutput)
			{
				std::cout << "Currently at line " << linecounter << std::endl;
			}
            auto parameters = modelPointer->initModel(linestr);


            if(args.FirstLine == args.LastLine ) modelPointer->write();

            auto EWPT = Minimizer::PTFinder_gen_all(modelPointer,0,300);



            if(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS)
			{
                if(C_PT*EWPT.Tc < EWPT.vc)
				{
					std::vector<double> vcritical,vbarrier;
                    vcritical = EWPT.EWMinimum;

                    std::vector<double> VEVSymmetric(modelPointer->get_nVEV());
					for(std::size_t i=0;i<VEVSymmetric.size();i++) VEVSymmetric.at(i) = 0;


                    std::vector<double> basepoint,basepointPot;



					double DeltaT = 1e-2;
					for(int tcounter=0;tcounter<=100;tcounter++){
						basepoint.clear();
						basepointPot.clear();
						double lineparam = tcounter*DeltaT;
                        for(std::size_t i=0;i<modelPointer->get_nVEV();i++) {
							basepoint.push_back(VEVSymmetric.at(i)+ lineparam*(vcritical.at(i) - VEVSymmetric.at(i)));

						}

                        double Temp = EWPT.Tc;
                        auto MinPlaneResult = Minimizer::MinimizePlane(basepoint,VEVSymmetric,vcritical,modelPointer,Temp);
                        double Vmin = MinPlaneResult.PotVal;
                        auto MinimumPlane = MinPlaneResult.Minimum;
                        basepointPot=modelPointer->MinimizeOrderVEV(basepoint);
						double Vline = modelPointer->VEff(basepointPot,Temp);
						double DistanceLineMinimum = 0;
                        for(std::size_t i=0;i<modelPointer->get_nVEV();i++){
							DistanceLineMinimum += std::pow(basepoint.at(i)-MinimumPlane.at(i),2);
						}
						DistanceLineMinimum = std::sqrt(DistanceLineMinimum);





						outfile <<  linestr;
                        for(auto x: parameters.second) outfile << sep << x;
                        outfile << sep << EWPT.Tc << sep << EWPT.vc;
                        outfile << sep << EWPT.vc / EWPT.Tc;
                        for(auto x: EWPT.EWMinimum) outfile << sep << x;
                        outfile << sep <<lineparam;
                        for(std::size_t i=0;i<basepoint.size();i++) outfile << sep << basepoint.at(i);
                        outfile << sep << Vline;
                        for(std::size_t i=0;i<MinimumPlane.size();i++) outfile << sep << MinimumPlane.at(i);
                        outfile << sep << Vmin;
                        outfile << sep << DistanceLineMinimum;
						outfile << std::endl;
					}







				}
			}

            if(args.FirstLine == args.LastLine) {

				if(dimensionnames.size() != ndim +3){
					std::cout << "The number of names in the function addLegendTemp does not match the number of vevs, going to default naming."
							<< "You should fix this as this will result in errors in your output file." << std::endl;
                    std::cout << "Succeded ? " << EWPT.StatusFlag << sep << " (1 = Success , -1 = v/T reached a value below " << C_PT << " during the calculation) \n";
                    if(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS){std::cout << "omega_c = " << EWPT.vc << " GeV\n";
                    std::cout << "T_c = " << EWPT.Tc << " GeV\n";
                    std::cout << "xi_c = omega_c/T_c =  " << EWPT.vc/EWPT.Tc << std::endl;
                    for(std::size_t i=0;i<ndim ;i++) {
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
                        for(std::size_t i=0;i<ndim; i++){
                            std::cout << dimensionnames.at(i+3) << " = " << EWPT.EWMinimum.at(i) << " GeV\n";
						}

					}
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
