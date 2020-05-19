/*
 * Kfactors.h
 *
 *  Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller

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


#include <BSMPT/models/IncludeAllModels.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte.h>
#include <vector>

/**
 * @file
 */

namespace BSMPT{
namespace Kfactors{


/**
 * @brief The GSL_integration struct used to pass information to the GSL minimisation routine
 */
struct GSL_integration{
  double Temp;
  int switchval;
  double masssquared;
  int s;
  double vw;
};


/**
 * Helper function which displays the numbers in a more readable way
 */
void display_results (std::string title, double result, double error);

/**
 * Calculates the distribution f0, defined in Eq(13) in 0604159
 * @param E0 Energy at which the distribution is calculated
 * @param s s=1 for fermions and s=-1 for bosons
 * @param Temp Temperature at which the distribution should be evaluated
 * @param diff number of derivatives w.r.t. E0, diff = 0 no derivative and diff = 1 or 2 accordingly
 */
double distribution_f0(double E0, int s, double Temp, int diff=0);

/**
 * Evaluates the integrand in the K functions
 * @param p momentum at which to evaluate the K function
 * @param masssquared the m^2 value at which K should be evaluated
 * @param Temp the temperature at which K should be evaluated
 * @param switchvalue The index which K to choose
 * @param s s=-1 yields the function for a bosonic distribution, s=+1 for a fermionic one
 */
double K_integrand(const std::vector<double>& p, double masssquared, int switchvalue,int s, double Temp);

/**
 * Interface to communicate the value of K_integrand() to the GSL integration routine
 */
double K_integrand_gsl(double *x, std::size_t dim, void *params );

/**
 * Calculates the values of the K functions given in Eq (23) in 0604159 without the normalization
 * @param masssquared the m^2 value at which K should be evaluated
 * @param Temp the temperature at which K should be evaluated
 * @param switchvalue The index which K to choose
 * @param s s=-1 yields the function for a bosonic distribution, s=+1 for a fermionic one
 */
double K_integration(double masssquared, double Temp, int switchvalue, int s);

/**
 * Calculates the normalized K function
 * @param masssquared the m^2 value at which K should be evaluated
 * @param Temp the temperature at which K should be evaluated
 * @param switchvalue The index which K to choose
 * @param s s=-1 yields the function for a bosonic distribution, s=+1 for a fermionic one
 */
double K_functions(double masssquared, double Temp, int switchvalue,int s);
/**
 * Integrand to calculate the normalization for the Ktilde functions. Used in Ktilde_normalization()
 * @param x vector which containts the momentum
 * @param params GSL_integration struct
 */
double Ktilde_normalization_func(double x, void *params);
/**
 * Calculates the normalization for the non tilde K functions
 * @param Temp Temperature
 * @param s +1 for fermions and -1 for bosons
 * @param masssquared m^2 value
 */
double Ktilde_normalization(double Temp, int s, double masssquared);

}
}
