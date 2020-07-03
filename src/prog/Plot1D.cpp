/*
 * Plot1D.cpp
 *
 *  Created on: Jul 6, 2016
 *      Author: basler
 */

/**
 * @file
 * This program calculates the Minimum for a given parameterpoint at a given Temperature and plots
 * the potential as a function V(t*v,Temp) where v is the vector showing to the minimum coming from 0.
 * t is in the range of [FacStart,FacEnd] with FacStep number of Steps.
 */

#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <iostream>
using namespace std;

int main(int argc, char *argv[]) {


    bool Debug = true;
    if(!(argc == 9))
    	{
    		std::cout << "./Plot1D InputType Infile Outfile LinNum Temp FacStart FacEnd FacStep \n";
    		return -1;
    	}
	char* in_file; char* out_file;
    	double InputType;
    	sscanf(argv[1],"%lf",&InputType);
    	in_file = argv[2];
    	out_file = argv[3];
    	double LineNumb,Temp,FacStart,FacEnd,FacStep;
    	sscanf(argv[4],"%lf",&LineNumb);
    	sscanf(argv[5],"%lf",&Temp);
    	FacStart = atof(argv[6]);
    	FacEnd = atof(argv[7]);
    	FacStep = atof(argv[8]);



	std::vector<double> sol,start;
	std::vector<double> Weinberg,parCTVec;

	Class_Potential_Origin * modelPointer;
	Fchoose(modelPointer,InputType);

	std::ifstream infile(in_file);
	if(!infile.good()) return -1;


	std::string linestr;
	int linecounter = 1;

	double tmp;
	int nPar,nParCT;
	nPar = modelPointer->nPar;
	nParCT = modelPointer->nParCT;

	int dim = modelPointer->nVEV;
	double par[nPar];
	double parCT[nParCT];

	bool found=false;


	while(true)
	{
	   if(infile.eof()) break;
	   std::getline(infile,linestr);
	   if(linecounter == 1){
		   modelPointer->setUseIndexCol(linestr);
	   }
	   else if(linecounter == LineNumb)
	   {
            modelPointer->ReadAndSet(linestr,par);
            found=true;
	   }

	   else if(linecounter > LineNumb) break;
	   linecounter++;
	   if(infile.eof()) break;
	}
	infile.close();
	if(!found) {std::cout << "Line not found!" << std::endl; return -1;}
	std::pair<std::vector<double>,std::vector<double>> parameters = modelPointer->initModel(linestr);
	par=parameters.first;
	parCT = parameters.second;

	modelPointer->write();
	std::ofstream outfile(out_file);
	std::vector<double> vTree;
	for(int k=0;k<dim;k++) vTree.push_back(modelPointer->vevTreeMin.at(k));

	std::vector<double> Check;
	double vev=0;
	double v1T,v2T,v3T;
	std::vector<double> NTempMass;
	std::vector<double> zero;
	for(int k=0;k<dim;k++) zero.push_back(0);
	std::vector<double> point,pointPot;
	if(Debug) std::cout << "Start of Minimization" << std::endl;
	Minimize_gen_all(InputType,par,parCT,Temp,sol,Check,vTree);
	if(Debug) std::cout << "End of Minimization" << std::endl;
	double VR1=0,VR2=0;
	std::vector<double> solPot;
	modelPointer->MinimizeOrderVEV(sol,solPot);
	vev = modelPointer->EWSBVeV(solPot);
	if(vev == 0)
	  {
	    std::cout << "Set to Tree-Level " << std::endl;
	    for(int i=0;i<dim;i++) sol[i] = modelPointer->vevTreeMin.at(i);
	  }

	for(int i=0;i<dim;i++) outfile << sol.at(i) << sep;
	outfile << std::endl;
	for(int nStep=0;nStep <= FacStep;nStep++)
		{
		        double tStep = FacStart + (FacEnd-FacStart)/FacStep * nStep;
		        point.clear();
		        pointPot.clear();
		        for(int i=0;i<dim;i++) point.push_back(tStep*sol.at(i));
		        modelPointer->MinimizeOrderVEV(point,pointPot);
		        VR1 = modelPointer->VEff(pointPot,Temp,0);
		        point.clear();
			outfile << tStep;
			outfile <<  sep << VR1;
			outfile << std::endl;
		}

   outfile.close();
   delete modelPointer;



	return 0;
}


