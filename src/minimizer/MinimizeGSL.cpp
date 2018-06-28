/*
 * MinimizeGSL.cpp
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

/**
 * @file
 * Using the Nelder-Mead Simplex algorithm, implemented in gsl, to find multiple local minima of the model
 * and compare them to find a candidate for the global minimum.
 */


#include "Minimizer.h"

const double GSL_Tolerance=std::pow(10,-4);

bool Minimize1D = false;

//print_state (size_t iter, gsl_multiroot_fdfsolver * s)


/**
 * struct containing the required Parameters of the model for the gsl interface
 */
struct GSL_params {
//int Model;
int nVEV;
//Class_Potential_Origin * modelPointer;
std::shared_ptr<Class_Potential_Origin> modelPointer;
double Temp;
};




/**
 * Calculates the value of the effective potential at the vev v and temperature p->Temp for the gsl interface
 */
double GSL_VEFF_gen_all(const gsl_vector *v, void *p)
{

	struct GSL_params * params = (struct GSL_params *) p;

	bool Debug = false;
	if(Debug) std::cout << "Start Debuging " << __func__ << std::endl;

	std::vector<double> vIn,vMin;
	int nVEVs = params->modelPointer->nVEV;

	for(int i=0;i<nVEVs;i++)
	{
	    vMin.push_back(gsl_vector_get(v,i));
	}

	vIn.resize(params->modelPointer->NHiggs);
	params->modelPointer->MinimizeOrderVEV(vMin,vIn);





	double res = params->modelPointer->VEff(vIn,params->Temp,0);


	return res;

}


/**
 * Calculates the next local minimum in the model from the point start
 * @returns The final status of the gsl minimization process.
 */
int GSL_Minimize_From_S_gen_all(void* p, std::vector<double>& sol,std::vector<double> start)
{

  bool Debug = false;
  if(Debug) std::cout << "Start Debuging " << __func__ << std::endl;
    gsl_set_error_handler_off();

	struct GSL_params * params = (struct GSL_params *) p;
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;

	double ftol = GSL_Tolerance;
	int MaxIter = 600;

	size_t iter = 0;
	int status;
	double size;

	int dim = params->nVEV;

	/* Starting point */
	x = gsl_vector_alloc (dim);
	for(int k=0;k<dim;k++) gsl_vector_set(x,k,start.at(k));
	ss = gsl_vector_alloc (dim);
	gsl_vector_set_all (ss, 1.0);



	/* Initialize method and iterate */
	minex_func.n = dim;
	minex_func.f = &GSL_VEFF_gen_all;
	minex_func.params = params;
	s = gsl_multimin_fminimizer_alloc (T, dim);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
	{
	  iter++;
	  status = gsl_multimin_fminimizer_iterate(s);

	  if (status)
		break;

	  size = gsl_multimin_fminimizer_size (s);
	  status = gsl_multimin_test_size (size, ftol);

	  	 if(Debug){ printf ("%5d %10.3e %10.3e %10.3e %10.3e f() = %7.7f size = %.7f\n",
	  	  	                iter,
	  	  	                gsl_vector_get (s->x, 0),
	  	  	                gsl_vector_get (s->x, 1),
	  				 gsl_vector_get (s->x, 2),
	  				 gsl_vector_get (s->x, 3),
	  	  	                s->fval, size);}


	}
  while (status == GSL_CONTINUE && iter < MaxIter);

  if(status == GSL_SUCCESS)
  {
	  for(int k=0;k<dim;k++) sol.push_back(gsl_vector_get(s->x,k));
  }
  else{
	  for(int k=0;k<dim;k++) sol.push_back(0);
  }

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return status;

}


/**
 * Minimize the Potential from different random starting points and choose the local minimum with the
 * deepest potential value as the candidate for the global minimum
 * @return True if a candidate for the global minimum is found and false otherwise
 */
bool GSL_Minimize_gen_all(int Model, const std::vector<double>& par, const std::vector<double>&  parCT,
		double Temp, std::vector<double>& solV, int seed,int MaxSol){
	std::vector<std::vector<double>> saveAllMinima;
	bool result = GSL_Minimize_gen_all(Model,par, parCT,Temp,solV, seed ,saveAllMinima, MaxSol);

//	std::cout << "rows = " << saveAllMinima.size() << std::endl;
	return result;
}


bool GSL_Minimize_gen_all(int Model, const std::vector<double>& par, const std::vector<double>&  parCT,
		double Temp, std::vector<double>& solV, int seed){
	std::vector<std::vector<double>> saveAllMinima;
	int MaxSol = 20;
	bool result = GSL_Minimize_gen_all(Model,par, parCT,Temp,solV, seed ,saveAllMinima, MaxSol);

//	std::cout << "rows = " << saveAllMinima.size() << std::endl;
	return result;
}

bool GSL_Minimize_gen_all(int Model, const std::vector<double>& par, const std::vector<double>&  parCT, double Temp,
		std::vector<double>& solV, int seed , std::vector<std::vector<double>>& saveAllMinima, int MaxSol)
{
	bool Debug = false;
	if(Debug) std::cout << "Start Debuging " << __func__ << std::endl;
	if(Debug) std::cout << "Searching for " << MaxSol << " solutions " << std::endl;
//	Class_Potential_Origin * modelPointer;
//	Fchoose(modelPointer,Model);
//	int nPar = modelPointer->nPar;
	std::shared_ptr<Class_Potential_Origin> modelPointer = FChoose(Model);
	int nParCT = modelPointer->nParCT;
//	double *p = (double *)par;
//	double *pCT = (double *)parCT;
	struct GSL_params params;
	params.Temp = Temp;
	params.nVEV = modelPointer->nVEV;
	modelPointer->set_All(par,parCT);
	params.modelPointer = modelPointer;





	int dim = modelPointer->nVEV;

	if(Debug) std::cout << "Setup done" << std::endl;




	std::default_random_engine randGen(seed);
	double RNDMax= 500;
	int MaxTries = 600;
	int tries = 0;
	int numOfSol = 0;
//	int MaxSol = 20;
	int nCol = dim+2;
//	double SolAr[MaxSol][nCol];
//	std::vector<std:vector<double>> SolAr;
	std::vector<double> start,sol,vPot;
	if(Debug) std::cout << "Start Loop" << std::endl;
	double timeStart = time(NULL);
	do{
		start.resize(dim);
		for(int i=0;i<dim;i++) start[i] = RNDMax*(-1+2*std::generate_canonical<double,std::numeric_limits<double>::digits>(randGen));


		GSL_Minimize_From_S_gen_all(&params,sol,start);

		vPot.resize(modelPointer->NHiggs);
		modelPointer->MinimizeOrderVEV(sol,vPot);



		if(Debug) std::cout << "Finished from S" << std::endl;
//		for(int i=0;i<dim;i++) SolAr[numOfSol][i] = sol.at(i);
//		SolAr[numOfSol][dim] = modelPointer->EWSBVEV(vPot);
//		SolAr[numOfSol][dim+1] = modelPointer->VEff(vPot,Temp,0);


		std::vector<double> row(nCol);
		for(int i=0;i<dim;i++) row.at(i) = sol.at(i);
		row.at(dim) = modelPointer->EWSBVEV(vPot);
		row.at(dim+1) = modelPointer->VEff(vPot,Temp,0);

		saveAllMinima.push_back(row);


		numOfSol ++;
		if(numOfSol == MaxSol) break;


		start.clear();
		sol.clear();
		tries++;
	}while(tries <= MaxTries);
	double timeEnd = time(NULL);
	if(Debug) std::cout << "Time for GSL = " << timeEnd-timeStart << std::endl;
	if(Debug) std::cout << "End Loop" << std::endl;
	if(Debug) std::cout << "GSL used " << tries << " Starting points to find  " << numOfSol << " Soltuions " << std::endl;

	if(numOfSol == 0)
	{
		std::cout << "No solutions found \n";
//		delete modelPointer;
		return false;
	}

	int minIndex = 0;
	double VMin = saveAllMinima[0][dim+1];
	for(int k=1;k<numOfSol;k++)
	{
		if(saveAllMinima[k][dim+1] <= VMin)
		{
			VMin = saveAllMinima[k][dim+1];
			minIndex = k;
		}
	}


	for(int k=0;k<dim+2;k++) solV.push_back(saveAllMinima[minIndex][k]);

	if(Debug){
	    std::cout << "Solution found by GSL is \n";
	    for(int k=0;k<dim+2;k++) std::cout << saveAllMinima[minIndex][k] << "\t";
	    std::cout << std::endl;
	}


//	delete modelPointer;
	return true;








}





