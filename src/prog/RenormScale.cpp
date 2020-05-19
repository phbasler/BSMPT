/*
 * RenormScale.cpp
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
 * Calculates xi_C as a function of the renormalised scale mu. The renormalisation scale mu is varied from 1/2 to 1.5 C_vev0 in NumberOfStep steps.
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
#include <BSMPT/utility.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace BSMPT;

int main(int argc, char *argv[]) try{
    /**
     * PrintErrorLines decides if parameter points with no valid EWPT (no NLO stability or T=300 vanishing VEV)
     * are printed in the output file
     */
    bool PrintErrorLines=true;

    if(!( argc == 6 or argc == 7) )
    {
        std::cerr << "./RenormScale Model Inputfile Outputfile  Line NumberOfSteps TerminalOutput(y)\n";
        ShowInputError();
        return EXIT_FAILURE;
    }


    auto Model=ModelID::getModel(argv[1]);
    if(Model==ModelID::ModelIDs::NotSet) {
        std::cerr << "Your Model parameter does not match with the implemented Models." << std::endl;
        ShowInputError();
        return EXIT_FAILURE;
    }

    double LineNumber,NumberOfStep;
    char* in_file;char* out_file;

    in_file = argv[2];
    out_file = argv[3];

    LineNumber = atoi(argv[4]);
    NumberOfStep = atoi(argv[5]);


    bool TerminalOutput = false;
    if(argc == 7) {
        std::string s7 = argv[6];
        std::cout << s7 << std::endl;
        TerminalOutput = ("y" == s7);
        if(s7!="y") throw std::runtime_error("To many arguments? Do you want terminal output?");
    }
    if(LineNumber < 1)
    {
        std::cout << "Start line counting with 1" << std::endl;
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
    std::size_t nPar,nParCT;
    nPar = modelPointer->get_nPar();
    nParCT = modelPointer->get_nParCT();
    std::vector<double> par(nPar);
    std::vector<double> parCT(nParCT);

    while(getline(infile,linestr))
    {
        if(linecounter > LineNumber) break;

        if(linecounter == 1)
          {

            outfile << linestr<<sep << "mu_factor"<<sep<<"mu"<<sep<<modelPointer->addLegendTemp()<<sep<<modelPointer->addLegendVEV()<<sep<<modelPointer->addLegendCT()<<sep<<"BSMPT_StatusFlag"<<std::endl;

            modelPointer->setUseIndexCol(linestr);
          }
        if(LineNumber==linecounter)
        {
            std::pair<std::vector<double>,std::vector<double>> parameters = modelPointer->initModel(linestr);
            par=parameters.first;
            parCT = parameters.second;
            if(TerminalOutput) modelPointer->write();

            if(TerminalOutput){
                std::cout<<"Calculating EWPT in default settings with:\n mu = "
                                      <<modelPointer->get_scale()<<std::endl;
            }
            auto EWPTdefault = Minimizer::PTFinder_gen_all(modelPointer,0,300);
            std::vector<double> vevsymmetricSolution,checksym, startpoint;
            for(std::size_t i=0;i<modelPointer->get_nVEV();i++) startpoint.push_back(0.5*EWPTdefault.EWMinimum.at(i));
            auto VEVsym = Minimizer::Minimize_gen_all(modelPointer,EWPTdefault.Tc+1,checksym,startpoint);
            if(TerminalOutput)std::cout<<"Start of mu variation"<<std::endl;
            double nstep = NumberOfStep;
            for(int step=0;step<nstep;step++){
                double mu_factor = 1/2. + (step/nstep);
                auto VEVnames = modelPointer->addLegendTemp();
                auto CT_mu=modelPointer->resetScale(C_vev0*mu_factor);
                auto EWPT_mu = Minimizer::PTFinder_gen_all(modelPointer,0,300);
                std::vector<double>checkmu;
                auto VEV_NLO_mu = Minimizer::Minimize_gen_all(modelPointer,0,checkmu,startpoint);
                if(TerminalOutput){
                    std::cout<<"\tStatusFlag = "<<EWPT_mu.StatusFlag<<" @ mu = "<<modelPointer->get_scale()<<" = "<<mu_factor<<"*C_vev0"<<std::endl;
                    for(std::size_t i=0;i<EWPT_mu.EWMinimum.size();i++)
                    {
                        std::cout<<sep<<sep<< VEVnames.at(i+3)<<sep<<EWPT_mu.EWMinimum.at(i)
                                <<" = " <<EWPT_mu.EWMinimum.at(i)/EWPTdefault.EWMinimum.at(i)
                               <<"/"<<VEVnames.at(i+3)<<std::endl;
                    }
                }
                if(PrintErrorLines){
                    outfile << linestr;
                    outfile << sep << mu_factor <<sep<< mu_factor*C_vev0;
                    outfile << sep << EWPT_mu.Tc<<sep<<EWPT_mu.vc<<sep<<EWPT_mu.vc/EWPT_mu.Tc<<sep<<EWPT_mu.EWMinimum;
                    outfile << sep << VEV_NLO_mu;
                    outfile << sep << CT_mu;
                    outfile << sep << EWPT_mu.StatusFlag;
                    outfile << std::endl;
                }
                else if(EWPTdefault.StatusFlag == Minimizer::MinimizerStatus::SUCCESS)
                {
                    if(C_PT* EWPTdefault.Tc < EWPTdefault.vc)
                    {
                        outfile << linestr;
                        outfile << sep << mu_factor <<sep<< mu_factor*C_vev0;
                        outfile << sep << EWPT_mu.Tc<<sep<<EWPT_mu.vc<<sep<<EWPT_mu.vc/EWPT_mu.Tc<<sep<<EWPT_mu.EWMinimum;
                        outfile << sep << VEV_NLO_mu;
                        outfile << sep << CT_mu;
                        outfile << sep << EWPT_mu.StatusFlag;
                        outfile << std::endl;
                    }
            }

        }//END: Mu Factor
        }//END:LineCounter
        linecounter++;
        if(infile.eof()) break;
    }
    if(TerminalOutput) std::cout << std::endl;
    outfile.close();

    return EXIT_SUCCESS;

}

catch(exception& e){
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
}
