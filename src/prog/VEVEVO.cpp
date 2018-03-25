/*
 * Plot_gen.cpp
 *
 *  Created on: Mar 3, 2017
 *      Author: basler
 */




/**
 * @file
 * This program calculates the development of the VeVs with the Temperature.
 */


#include "../models/IncludeAllModels.h"
#include "../minimizer/Minimizer.h"
#include <iostream>
using namespace std;

int main(int argc, char *argv[]) try{


    bool Debug = false;
	if(!(argc == 8))
	{
		std::cout << "./VEVEVO Model Inputfile Outputfile Line Tempstart Tempstep Tempend \n";
		std::cout << "\t0: C2HDM \n\t1: R2HDM\n\t2: RN2HDM \n\t3: DarkN2HDM \n";
		return -1;
	}
	char* in_file; char* out_file;
	in_file = argv[2];
	out_file = argv[3];
	double LineNumb,TempStartIn,TempStepIn,TempEndIn;
	double TempStart,TempEnd,TempStep;
	double Model;



	Model = atoi(argv[1]);
	LineNumb = atoi(argv[4]);
	TempStart = atof(argv[5]);
	TempStep = atof(argv[6]);
	TempEnd = atof(argv[7]);

	if(LineNumb < 1)
	{
		std::cout << "Start line counting with 1" << std::endl;
		return EXIT_FAILURE;
	}


	if(TempStart < 0)
	{
		std::cout << "The starting value of your Temperature was negative. This was corrected to TempStart = 0"
				<<std::endl;
		TempStart = 0;
	}
	if(TempEnd < TempStart)
	{
		std::cout << "The value of Tempend was lower then the value of Tempstart. This was corrected by swapping them"
				<<std::endl;
		double tmp=TempStart;
		TempStart = TempEnd;
		TempEnd=tmp;

	}



	std::vector<double> sol,start,solPot;
	std::vector<double> Weinberg,parCTVec;


	std::unique_ptr<Class_Potential_Origin> modelPointer = FChoose(Model);

	if(Debug) std::cout << "Set model pointer " << std::endl;

	std::ifstream infile(in_file);
	if(!infile.good()) {
			std::cout << "Input file not found " << std::endl;
			return EXIT_FAILURE;
	}

	if(Debug) std::cout << "found file " << std::endl;

	std::ofstream outfile(out_file);
	if(!outfile.good())
	{
		std::cout << "Can not create file " << out_file << std::endl;
		return EXIT_FAILURE;
	}

	std::string linestr;
	int linecounter = 1;

	double tmp;
	int nPar,nParCT;
	nPar = modelPointer->nPar;
	nParCT = modelPointer->nParCT;

	int dim = modelPointer->nVEV;
	std::vector<double> par(nPar);
	std::vector<double> parCT(nParCT);

	bool found=false;

	if(Debug) std::cout << "Read start " << std::endl;

	while(true)
	{
	   if(infile.eof()) break;
	   std::getline(infile,linestr);
	   if(linecounter == LineNumb)
	   {
            modelPointer->ReadAndSet(linestr,par);
            found=true;
	   }

	   else if(linecounter > LineNumb) break;
	   linecounter++;
	   if(infile.eof()) break;
	}
	infile.close();
	if(!found) {std::cout << "Line not found !\n"; return -1;}

	if(Debug) std::cout << "Read done " << std::endl;

	if(Debug)std::cout<<"Calculating counterterms"<<std::endl;
	modelPointer->calc_CT(parCT);





	if(Debug)std::cout<<"set function call"<<std::endl;
	modelPointer->set_All(par,parCT);

	if(Debug) modelPointer->write();



	std::vector<double> vTree;


	for(int k=0;k<dim;k++) vTree.push_back(modelPointer->vevTreeMin.at(k));

	if(Debug)
	{
		std::vector<double> TreeMasses;
		modelPointer->HiggsMassesSquared(TreeMasses,modelPointer->vevTree,0,0);
		for(int i=0;i<modelPointer->NHiggs;i++)
		{
			std::cout << TreeMasses[i] << "\t" << std::sqrt(std::abs(TreeMasses[i])) << std::endl;
		}
	}








	std::vector<double> Check;
	double vev=0;
	double v1T,v2T,v3T;
	std::vector<double> NTempMass;

	std::vector<double> zero;
	for(int k=0;k<dim;k++) zero.push_back(0);

	std::vector<double> sol1D;


	double Temp = TempStart;
	std::cout << std::scientific;
	std::cout << std::setprecision(16);
	outfile << std::setprecision(16);


    outfile << "T" << "\t" << "v" << "\t";
    outfile << modelPointer->addLegendVEV()
    		<< "\t" << "Veff(v,T)"
    		<< std::endl;

   for(double Temp = TempStart; Temp<=TempEnd; Temp+=TempStep)
   {

	   start.clear();
	   if(Temp==TempStart)
	   {
		   for(int k=0;k<dim;k++) start.push_back(vTree.at(k));
	   }
	   else{
		   for(int k=0;k<dim;k++) start.push_back(sol.at(k));
	   }
	   sol.clear();
	   Check.clear();
	   solPot.clear();
	   if(Debug)std::cout<<"Minimization start"<<std::endl;
	   Minimize_gen_all(Model,par,parCT,Temp,sol,Check,start,3);
	   if(Debug)std::cout<<"Minimization end"<<std::endl;
	   if(Debug){
		   for(int i=0;i<dim;i++) std::cout << sol[i] << std::endl;
	   }
	   modelPointer->MinimizeOrderVEV(sol,solPot);
	   vev = modelPointer->EWSBVEV(solPot);




	   outfile << Temp << "\t";
	   outfile << vev;
	   for(int k=0;k<dim;k++) outfile << "\t" << sol.at(k);
	   outfile << "\t" << modelPointer->VEff(solPot,Temp,0);
	   outfile << std::endl;


   }
   outfile.close();

	return EXIT_SUCCESS;
}
catch(exception& e){
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
}
