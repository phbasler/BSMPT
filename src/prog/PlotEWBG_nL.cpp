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
using namespace Baryo;

int main(int argc, char *argv[]) try{
    if((argc != 7) and (argc !=8))
    {
        std::cerr << "./PlotEWBG_vw Model Inputfile Outputfile Line vw EWBGConfigFile TerminalOutput(y/n)\n";
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
    double vwIn;
    vwIn = atof(argv[5]);

    if(LineNumb < 1)
    {
        std::cerr << "Start line counting with 1" << std::endl;
        return EXIT_FAILURE;
    }
    bool TerminalOutput = false;
    if(argc ==8){
        std::string s7 = argv[7];
        std::cout<<"Terminal Output:"<<s7<<std::endl;
        TerminalOutput = ("y" == s7);
    }
//Set up of BSMPT/Baryo Classes
    Baryo::CalculateEtaInterface EtaInterface(argv[6] /* = Config File */);
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
            outfile << linestr<<sep;
            outfile << "T_c_nLvar"<<sep<<"omega_c_nLvar"<<sep<<"vw_nLvar"<<sep<<"LW_nLvar" <<sep;
            outfile << "z" << sep <<"z/LW" << sep << "nL_VIA" << sep << "muL_FH" << sep << "nL_FH";
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
    int WhichMinimizer = 1;/////
////////////////////////////////
    auto EWPT = Minimizer::PTFinder_gen_all(modelPointer,0,300,WhichMinimizer);//Change to NLOpt
//SFOEWPT FOUND
    if(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS and C_PT*EWPT.Tc < EWPT.vc){
        if(TerminalOutput) std::cout<<"SFOEWPT found..."<<std::endl;
        std::vector<double> vcritical,vbarrier;
        vcritical = EWPT.EWMinimum;
        double TC = EWPT.Tc;
        std::vector<double> MinimumPlane;
//Find the minimum in the symmetric phase. For this minimise at T = Tc + 1
        std::vector<double> vevsymmetricSolution,checksym, startpoint;
        for(size_t i=0;i<modelPointer->get_nVEV();i++) startpoint.push_back(0.5*vcritical.at(i));
        vevsymmetricSolution=Minimizer::Minimize_gen_all(modelPointer,TC+1,checksym,startpoint,WhichMinimizer);//TODO: Change to NLopt

/////////////////////////////////////////////////////////////////////////////////
        size_t nstep = 100;

        if(TerminalOutput) std::cout<<"Set up the numerics "<<std::endl;
        EtaInterface.setNumerics(vwIn , EWPT.EWMinimum , vevsymmetricSolution , TC , modelPointer);//Set up parameter container for Baryo Calculation-->Calls container.init
        if(TerminalOutput) std::cout<<"Starting setting the class instances"<<std::endl;
        BSMPT::Baryo::tau_source C_tau;
        EtaInterface.GSL_integration_mubl_container.set_transport_method(3);//setting to tau class
        bool botflag = true;
        auto class_GamM = EtaInterface.get_class_CalcGamM();
        auto class_ScP  = EtaInterface.get_class_Scp();
        auto class_kappa = EtaInterface.get_class_kappa();
        C_tau.set_class(botflag,EtaInterface.GSL_integration_mubl_container,class_GamM,class_ScP,class_kappa);

        auto tau_arr_nL = set_up_nL_grid(nstep,EtaInterface.GSL_integration_mubl_container,C_tau);
        auto FH_GSL = generate_mubl_spline(EtaInterface.GSL_integration_mubl_container , static_cast<int>(nstep));

///////
/// Outfile
///////
        for(size_t i=0;i<nstep;i++){
            outfile<<linestr<<sep;
            outfile<<TC<<sep<<EWPT.vc<<sep<<vwIn<<sep<<EtaInterface.getLW()<<sep;
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
catch(exception& e){
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
}
