/*
 * Minimizer.h
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
 */

#ifndef MINIMIZER_H_
#define MINIMIZER_H_

#include <vector>                   // for vector
#include <memory>
#include <BSMPT/config.h>
#include <BSMPT/models/IncludeAllModels.h>



namespace BSMPT {
class Class_Potential_Origin;
namespace Minimizer{


// Here you can decide the default setting for the mnimizer

/**
 * @brief UseGSLDefault Use the GSL minimizer in the default settings
 */
const bool UseGSLDefault = true;

/**
 * @brief UseLibCMAESDefault Use the Libcmaes minimizer in the default settings
 */
#ifdef CMAES_FOUND
const bool UseLibCMAESDefault = true;
#else
const bool UseLibCMAESDefault = false;
#endif

/**
 * @brief UseNLoptDefault Use the NLopt minimizer in the default settings
 */
#ifdef NLopt_FOUND
const bool UseNLoptDefault = true;
#else
const bool UseNLoptDefault = false;
#endif

/**
 * @brief CalcWhichMinimizer Calculates the WhichMinimizer value with the given Minimizer options
 * @param UseGSL Should GSL be used?
 * @param UseCMAES Should CMAES be used?
 * @param UseNLopt Should NLopt be used?
 * @return
 */
constexpr int CalcWhichMinimizer(bool UseGSL = UseGSLDefault,
                                 bool UseCMAES =UseLibCMAESDefault,
                                 bool UseNLopt = UseNLoptDefault)
{
    return  static_cast<int>(UseCMAES)
            +2*static_cast<int>(UseGSL)
            + 4*static_cast<int>(UseNLopt);
}

/**
 * @brief WhichMinimizerDefault default value for the Minimizers to use
 */
constexpr int WhichMinimizerDefault = CalcWhichMinimizer();

// Minimizer.cpp
std::vector<double> Minimize_gen_all(
        const std::shared_ptr<Class_Potential_Origin>& modelPointer,
        const double& Temp,
        std::vector<double>& Check,
        const std::vector<double>& start,
        const int& WhichMinimizer=WhichMinimizerDefault);


/**
 * @brief The MinimizerStatus enum for the Statusflags of the minimizer
 */
enum class MinimizerStatus{
    SUCCESS = 1,
    NOTVANISHINGATFINALTEMP=-1,
    NOTNLOSTABLE=-2,
    NUMERICALLYUNSTABLE=-3,
    BELOWTHRESHOLD=-4,
    NLOVEVZEROORINF=-5

};

/**
 * @brief The EWPTReturnType struct
 * Contains the following information
 * @param StatusFlag = SUCCESS : successful
 * @param StatusFlag = NOTVANISHINGATFINALTEMP: v(T=TempEnd) != 0 => vc = v(TempEnd) , Tc = TempEnd
 * @param StatusFlag = NOTNLOSTABLE: Not vacuum stable at T=0 (if TempStart = 0) or v(TempStart) = 0 or v(TempStart) > 255 => Tc = TempEnd, vc = 0
 * @param StatusFlag = NUMERICALLYUNSTABLE: v(T) > 255 during the bisection (this implies the parameter point is numerically unstable) => vc and Tc are the last VEVs encountered
 * @param StatusFlag = BELOWTHRESHOLD: v/T < C_PT during the bisection =>  vc = Last VEV, TC = Last Temp
 * @param vc = critical VEV
 * @param Tc = critical Temperature
 * @param EWMinimum: The broken EW minimum
 */
struct EWPTReturnType{
    MinimizerStatus StatusFlag;
    double vc;
    double Tc;
    std::vector<double> EWMinimum;
};

/**
* Uses a bisection method between the Temperature TempStart and TempEnde to find the phase transition
*  in the Model and writes the solution in the vector sol.
*  @param modelPointer shared_ptr to the parameter point
*  @param TempStart Low temperature for the starting interval of the bisection method
*  @param TempEnd High temperature for the starting interval of the bisection method
*  @param WhichMinimizer Which minimizers should be taken? 1 = CMAES, 2 = GSL, 4 = NLOPT, to use multiple add the numbers
*  @return The information are returned in a EWPTReturnType struct
*/
EWPTReturnType PTFinder_gen_all(
        const std::shared_ptr<Class_Potential_Origin>& modelPointer,
        const double& TempStart,
        const double& TempEnd,
        const int& WhichMinimizer=WhichMinimizerDefault);

/**
 * @brief Minimize_gen_all_tree_level Minimizes the tree-level potential
 * @param Model Which Model to minimize
 * @param par parameters of the point
 * @param parCT counterterm parameters
 * @param Check Vector to safe the error flags during the minimization
 * @param start Starting point for the minimization
 * @param WhichMinimizer Which minimizers should be taken? 1 = CMAES, 2 = GSL, 4 = NLOPT, to use multiple add the numbers
 * @return the global minimum
 */
std::vector<double> Minimize_gen_all_tree_level(
        const ModelID::ModelIDs& Model,
        const std::vector<double>& par,
        const std::vector<double>& parCT,
        std::vector<double>& Check,
        const std::vector<double>& start,
        int WhichMinimizer=WhichMinimizerDefault);


/**
 * @brief FindNextLocalMinima finds the local minima from the given starting point at the given temperature
 * @param model Parameter point to minimize
 * @param StartingPoint Starting point from where to look for the next local minima
 * @param temperature Temperature at which to minimize the potential
 * @param WhichMinimizer Which Minimizer should be used? CMAES is not used here as it is only a global minimizer
 * @return A vector containing the local minima. In case different Minima find different minima all solutions are given.
 */
std::vector<std::vector<double>> FindNextLocalMinima(
        const std::shared_ptr<Class_Potential_Origin>& model,
        const std::vector<double>& StartingPoint,
        const double& temperature,
        int WhichMinimizer = WhichMinimizerDefault);


}
}

#endif /* MINIMIZER_H_ */
