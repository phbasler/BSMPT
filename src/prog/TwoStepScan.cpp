/*
*Two Step Scan for the N2HDM
*Created on: 3.4.2918
*   Author: Jonas
*/



/**
 * @file
 * This program scans for possible candidates for a two step PT
 */

#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <iostream>
using namespace std;

int main(int argc, char *argv[]) try{
	bool Debug = false;
	double steps = 10	;
	if(!(argc == 5))
	{
		std::cout << "./twsep input_file output_file linestart lineend \n";
		return -1;
	}
	char* in_file; char* out_file;
	in_file=argv[1];
	out_file=argv[2];
	double Model=2;
	double LineNumb,LineStart,LineEnd;
	LineStart=atoi(argv[3]);
	LineEnd=atoi(argv[4]);

	if(LineStart < 1)
	{
		std::cout << "Start line counting with 1" << std::endl;
		return EXIT_FAILURE;
	}


	std::vector<double> sol,start,solPot,Check;
	std::vector<double> Weinberg,parCTVec;
	std::unique_ptr<Class_Potential_Origin> modelPointer = FChoose(Model);
	if(Debug)std::cout<<"Setting model parameter"<<std::endl;
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

	int nPar,nParCT;
	nPar = modelPointer->nPar;
	nParCT = modelPointer->nParCT;

	int dim = modelPointer->nVEV;
	std::vector<double> par(nPar);
	std::vector<double> parCT(nParCT);
	bool found=false;

	while(true)
	{
	   if(linecounter>LineEnd)break;
	   std::getline(infile,linestr);
	   if(linecounter == 1)
		 {
		   outfile << linestr << sep;
		   for(int i=0;i<steps;i++){

			   outfile<<"vCB_T"<<100*i<<sep<<"vCP_T"<<100*i<<sep<<"v1_T"<<100*i<<sep<<"v2_T"<<100*i<<sep<<"vS_T"<<100*i<<sep;
		   }
		   outfile << std::endl;
		 }
		 std::cout << std::scientific;
		 std::cout << std::setprecision(16);
		 outfile << std::setprecision(16);
		 std::vector<std::vector<double> > result;
	   if(linecounter>=LineStart and linecounter<=LineEnd and linecounter!=1)
	   {
		   	modelPointer->resetbools();
            modelPointer->ReadAndSet(linestr,par);
	   		modelPointer->calc_CT(parCT);
			modelPointer->set_All(par,parCT);

			if(Debug and LineStart==LineEnd) modelPointer->write();
	//Calculating stuff begins :)
			std::vector<double> vTree;
			for(int k=0;k<dim;k++) vTree.push_back(modelPointer->vevTreeMin.at(k));
			double temp;
			double vev;
			if(Debug)std::cout<<"Loop over temperature steps starts"<<std::endl;
			for(int n=0;n<steps;n++){
				temp=100*n;
				start.clear();
				for(int i =0;i<dim;i++)start.push_back(vTree.at(i));
				Check.clear();
				solPot.clear();
				sol.clear();
				Minimize_gen_all(Model,par,parCT,temp,sol,Check,start);
				result.push_back(sol);
			}
			if(Debug)std::cout<<"Loop over temperature steps ends"<<std::endl;

//Outfile Action
			outfile << linestr;
			for(int j=0;j<steps;j++){
				for(int i=0;i<dim;i++)outfile<<sep<<result[j][i];
				}
			outfile<<std::endl;
			result.clear();
			}
			linecounter++;
			if(infile.eof())break;
		}
		outfile.close();
		return EXIT_SUCCESS;
	}
catch(exception& e){
	std::cerr << e.what() << std::endl;
	return EXIT_FAILURE;
	}
