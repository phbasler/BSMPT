/*
 * ClassTemplate.cpp
 *
 *  Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner

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

#include "ClassTemplate.h"
#include "IncludeAllModels.h"
using namespace Eigen;

/**
 * @file
 * Template for adding a new model class
 */

/**
 * Here you have to adjust NNeutralHiggs, NChargedHiggs, nPar (number of Lagrangian parameters AFTER
 *  using the tadpole conditions),
 * nParCT (number of counterterms) as well as nVEV (number of VEVs for minimization)
 */
Class_Template::Class_Template ()
{
  Model = C_ModelTemplate; // global int constant which will be used to tell the program which model is called
  NNeutralHiggs = 1; // number of neutral Higgs bosons at T = 0
  NChargedHiggs=0; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar = 2; // number of parameters in the tree-Level Lagrangian
  nParCT = 3; // number of parameters in the counterterm potential

  nVEV=1; // number of VEVs to minimize the potential

  NHiggs = NNeutralHiggs+NChargedHiggs;

}

Class_Template::~Class_Template ()
{
}

/**
 * returns a string which tells the user the chronological order of the counterterms. Use this to
 * complement the legend of the given input file
 */
std::string Class_Template::addLegendCT(){
  std::string out = "dT\tdlambda\tdmsquared";
  return out;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and the critical temperature. Use this to
 * complement the legend of the given input file
 */
std::string Class_Template::addLegendTemp(){
  std::string out = "T_c\tv_c";
  //out += "Your VEV order";
  out += "omega";
  return out;
}

/**
 * returns a string which tells the user the chronological order of the Triple Higgs couplings. Use this to
 * complement the legend of the given input file
 *
 */
std::string Class_Template::addLegendTripleCouplings(){
	std::vector<std::string> particles;
	particles.resize(NHiggs);
	//here you have to define the particle names in the vector particles

	particles[0]="H";

	     std::string out = "Tree_";
	     for(int i=0;i<NHiggs;i++)
	       {
	 	for(int j=i;j<NHiggs;j++)
	 	  {
	 	    for(int k=j;k<NHiggs;k++)
	 	      {
	 	    	if(!(i==0 and j==0 and k==0)) out += "\tTree_";
	 		out+=particles.at(i);
	 		out+=particles.at(j);
	 		out+=particles.at(k);
	 		out+="\tCT_";
	 		out+=particles.at(i);
	 		out+=particles.at(j);
	 		out+=particles.at(k);
	 		out+="\tCW_";
	 		out+=particles.at(i);
	 		out+=particles.at(j);
	 		out+=particles.at(k);
	 	      }
	 	  }
	       }

	     return out;
}

/**
 * returns a string which tells the user the chronological order of the VEVs. Use this to
 * complement the legend of the given input file
 */
std::string Class_Template::addLegendVEV(){
	std::string out;
	//out = "Your VEV order";
	out = "omega";
  return out;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_Template::ReadAndSet(const std::string& linestr, std::vector<double>& par )
{
	std::stringstream ss(linestr);
	double tmp;

	/**
	 * Reads first number if the legend started with a tabulator. This assumes that the first column is then an index.
	 */
	if (UseIndexCol){
		ss >> tmp;
	}

	double lms,llambda;


	for(int k=1;k<=2;k++)
	{
	      ss>>tmp;
	      if(k==1) lms = tmp;
	      else if(k==2) llambda = tmp;
	}
	par[0] = lms;
	par[1] = llambda;


	set_gen(par); // This you have to call so that everything will be set
	return ;
}


/**
 * Set Class Object as well as the VEV configuration
 */
void Class_Template::set_gen(const std::vector<double>& par) {


	ms = par[0];
	lambda = par[1];

	g=C_g;

	yt = std::sqrt(2)/C_vev0 * C_MassTop;

	scale = C_vev0;

	vevTreeMin.resize(nVEV);
	vevTree.resize(NHiggs);

	// Here you have to set the vector vevTreeMin. The vector vevTree will then be set by the function MinimizeOrderVEV
	vevTreeMin[0] = C_vev0;

	MinimizeOrderVEV(vevTreeMin,vevTree);
	if(!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the entries of Curvature_Higgs_CT_L1 to
 * Curvature_Higgs_CT_L4.
 */
void Class_Template::set_CT_Pot_Par(const std::vector<double>& par){

	dT = par[0];
	dlambda = par[1];
	dms = par[2];

	Curvature_Higgs_CT_L1[0] = dT;
	Curvature_Higgs_CT_L2[0][0] = dms;
	Curvature_Higgs_CT_L4[0][0][0][0] = dlambda;
}


/**
 * console output of all Parameters
 */
void Class_Template::write() {

	std::cout << "The parameters are : " << std::endl;
	std::cout << "lambda = " << lambda << std::endl
			<< "\tm^2 = " << ms << std::endl;

	std::cout << "The counterterm parameters are : " << std::endl;
	std::cout << "dT = "<< dT << std::endl
			<< "dlambda = " << dlambda << std::endl
			<< "dm^2 = "<< dms << std::endl;

	std::cout << "The scale is given by mu = " << scale << " GeV " << std::endl;

}


/**
 * Calculates the counterterms. Here you need to work out the scheme and implement the formulas.
 */
void Class_Template::calc_CT(std::vector<double>& par){
	bool Debug=false;
	if(Debug) std::cout << "Start" << __func__ << std::endl;

	if(!SetCurvatureDone)SetCurvatureArrays();
	if(!CalcCouplingsdone)CalculatePhysicalCouplings();
	if(Debug) {
	std::cout << "Couplings done " << std::endl;
	}
	std::vector<double> WeinbergNabla,WeinbergHesse;
	WeinbergFirstDerivative(WeinbergNabla);
	WeinbergSecondDerivative(WeinbergHesse);

	if(Debug) std::cout << "Finished Derivatives " << std::endl;

	VectorXd NablaWeinberg(NHiggs);
	MatrixXd HesseWeinberg(NHiggs,NHiggs),HiggsRot(NHiggs,NHiggs);
	for(int i=0;i<NHiggs;i++)
	{
		NablaWeinberg[i] = WeinbergNabla[i];
		for(int j=0;j<NHiggs;j++) HesseWeinberg(i,j) = WeinbergHesse.at(j*NHiggs+i);
	}

	// Here you have to use your formulas for the counterterm scheme
	double t = 0;
	dT = t;
	dlambda = 3.0*t/std::pow(C_vev0,3) + 3.0/std::pow(C_vev0,3) * NablaWeinberg(0) -3.0/std::pow(C_vev0,2) *HesseWeinberg(0,0);
	dms = -3.0/(2*std::pow(C_vev0,2)) *NablaWeinberg(0) + 1.0/2.0 *HesseWeinberg(0,0) -3.0*t/(2*C_vev0);

	par[0] = dT;
	par[1] = dlambda;
	par[2] = dms;

	set_CT_Pot_Par(par);

}




void Class_Template::TripleHiggsCouplings()
{
	bool Debug=false;
	if(Debug) std::cout << "Debug turned on in " << __func__ << std::endl;

	if(!SetCurvatureDone)SetCurvatureArrays();
	if(!CalcCouplingsdone)CalculatePhysicalCouplings();


	std::vector<double> HiggsOrder(NHiggs);
	// Here you have to set the vector HiggsOrder. By telling e.g. HiggsOrder[0] = 5 you always want your 6th lightest
	// particle to be the first particle in the vector (which has the index 5 because they are sorted by mass)

	// example for keeping the mass order
	for(int i=0;i<NHiggs;i++) {
		HiggsOrder[i]=i;
	}

	if(Debug) std::cout << "Calculate Derivative" << std::endl;
	std::vector<double> TripleDeriv;
	WeinbergThirdDerivative(TripleDeriv);
	if(Debug) std::cout << "Finished calculating derivatives " << std::endl;
	double GaugeBasis[NHiggs][NHiggs][NHiggs];
	for(int i=0;i<NHiggs;i++)
	  {
		for(int j=0;j<NHiggs;j++)
		{
		  for(int k=0;k<NHiggs;k++)
			{
			  GaugeBasis[i][j][k] = TripleDeriv.at(i+j*NHiggs+k*NHiggs*NHiggs);
			}
		}
	  }

	MatrixXd HiggsRot(NHiggs,NHiggs);
	for(int i=0;i<NHiggs;i++)
	{
		for(int j=0;j<NHiggs;j++)
		{
			HiggsRot(i,j) = HiggsRotationMatrix[i][j];
		}
	}


	MatrixXd HiggsRotSort(NHiggs,NHiggs);




	for(int i=0;i<NHiggs;i++)
	{
		HiggsRotSort.row(i) = HiggsRot.row(HiggsOrder[i]);
	}

	TripleHiggsCorrectionsCWPhysical.resize(NHiggs);
	TripleHiggsCorrectionsTreePhysical.resize(NHiggs);
	TripleHiggsCorrectionsCTPhysical.resize(NHiggs);
	for(int i=0;i<NHiggs;i++) {
		TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);
		TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);
		TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);
		for(int j=0;j<NHiggs;j++) {
			TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);
			TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);
			TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);
		}
	}

	if(Debug) std::cout << "Setup done " << std::endl;


	for(int i=0;i<NHiggs;i++)
	  {
		for(int j=0;j<NHiggs;j++)
		{
			for(int k=0;k<NHiggs;k++)
			{
			  TripleHiggsCorrectionsCWPhysical[i][j][k] = 0;
			  TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;
			  TripleHiggsCorrectionsCTPhysical[i][j][k] = 0;
			  for(int l=0;l<NHiggs;l++)
			  {
				  for(int m=0;m<NHiggs;m++)
				  {
					  for(int n=0;n<NHiggs;n++)
					  {
						  double RotFac = HiggsRotSort(i,l)*HiggsRotSort(j,m)*HiggsRotSort(k,n);
						  TripleHiggsCorrectionsCWPhysical[i][j][k] += RotFac*GaugeBasis[l][m][n];
						  TripleHiggsCorrectionsTreePhysical[i][j][k] += RotFac*LambdaHiggs_3[l][m][n];
						  TripleHiggsCorrectionsCTPhysical[i][j][k] += RotFac*LambdaHiggs_3_CT[l][m][n];

					  }
				  }
			  }
			}
		}
	  }



}

void Class_Template::SetCurvatureArrays(){
  /*
   *  Here you have to set the vectors
   *  Curvature_Higgs_L1,Curvature_Higgs_L2,Curvature_Higgs_L3,Curvature_Higgs_L4
   *  Curvature_Gauge_G2H2
   *  Curvature_Quark_F2H1, Curvature_Lepton_F2H1
   *  as described in the potential in the paper.
   */

	initVectors();
	SetCurvatureDone=true;
	for(int i=0;i<NHiggs;i++) HiggsVev[i] = vevTree[i];


	Curvature_Higgs_L2[0][0] = ms;
	Curvature_Higgs_L4[0][0][0][0] = lambda;

	Curvature_Gauge_G2H2[0][0][0][0] = 4*std::pow(g,2);

	Curvature_Quark_F2H1[1][0][0] = yt;
	Curvature_Quark_F2H1[0][1][0] = yt;

}


void Class_Template::MinimizeOrderVEV(const std::vector<double>& vevMinimizer, std::vector<double>& vevFunction){
  /*
   * Here you do the conversion from the result vector from the minimizer (which has nVEV entries) to a
   * vector with NHiggs entries which you can use to call the member functions
   */

	VevOrder.resize(nVEV);
	// Here you have to tell which scalar field gets which VEV.
	VevOrder[0] = 0;

	int count=0;
	if(vevFunction.size() != 0)
	{
		for(int i=0;i<NHiggs;i++)
			{
				if(i==VevOrder[count]) {
					vevFunction[i] =  vevMinimizer.at(count);
					count++;
				}
				else vevFunction[i] = 0;

			}
	 }
	 else{
		 for(int i=0;i<NHiggs;i++)
		 {
			 if(i==VevOrder[count]) {
				 vevFunction.push_back(vevMinimizer.at(count));
				 count++;
			 }
			 else vevFunction.push_back(0);
		 }
	}
}
bool Class_Template::CalculateDebyeSimplified(){
  return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass matrix and implement
   * your formula here and return true. The vector is given by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool Class_Template::CalculateDebyeGaugeSimplified()
{
	bool Debug = false;
  if(Debug) std::cout << "Debug turned on in Class_Template :: " << __func__ << std::endl;

  /*
     * Use this function if you calculated the Debye corrections to the gauge mass matrix and implement
     * your formula here and return true. The vector is given by DebyeGauge[NGauge][NGauge]
     */


  return false;
}
double Class_Template::VTreeSimplified(const std::vector<double>& v){
	UseVTreeSimplified = false;
	UseVTreeSimplified = true; // To use the simplified version
	double res = 0;

	double vIn = v[0];
	res = 0.5*ms*std::pow(vIn,2) + 1.0/24.0 * lambda * std::pow(vIn,4);

	return res;
}

double Class_Template::VCounterSimplified(const std::vector<double>& v)
{
	UseVCounterSimplified = false;
	UseVCounterSimplified = true; // To use the simplified version
	if(not UseVCounterSimplified) return 0;
	double res = 0;
	double vIn = v[0];
		res = 0.5*dms*std::pow(vIn,2) + 1.0/24.0 * dlambda * std::pow(vIn,4) + dT * vIn;
	return res;
}

void Class_Template::Debugging(const std::vector<double>& input, std::vector<double>& output){

}
