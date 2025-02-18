// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#ifndef MINIMIZER_H_
#define MINIMIZER_H_

#include <BSMPT/config.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <memory>
#include <thread>
#include <vector> // for vector

namespace BSMPT
{
class Class_Potential_Origin;
namespace Minimizer
{

// Here you can decide the default setting for the mnimizer

/**
 * @brief UseGSLDefault Use the GSL minimizer in the default settings
 */
const bool UseGSLDefault = true;

/**
 * @brief UseLibCMAESDefault Use the Libcmaes minimizer in the default settings
 */
#ifdef cmaes_FOUND
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

const std::size_t Num_threads = std::thread::hardware_concurrency();

/**
 * @brief CalcWhichMinimizer Calculates the WhichMinimizer value with the given
 * Minimizer options
 * @param UseGSL Should GSL be used?
 * @param UseCMAES Should CMAES be used?
 * @param UseNLopt Should NLopt be used?
 * @return
 */
constexpr int CalcWhichMinimizer(bool UseGSL   = UseGSLDefault,
                                 bool UseCMAES = UseLibCMAESDefault,
                                 bool UseNLopt = UseNLoptDefault)
{
  return static_cast<int>(UseCMAES) + 2 * static_cast<int>(UseGSL) +
         4 * static_cast<int>(UseNLopt);
}

/**
 * @brief The MinimizersToUse struct used as a return of GetMinimizers
 */
struct MinimizersToUse
{
  bool UseCMAES{UseLibCMAESDefault};
  bool UseGSL{UseGSLDefault};
  bool UseNLopt{UseNLoptDefault};
  MinimizersToUse(bool useCMAES, bool useGSL, bool useNLopt)
      : UseCMAES{useCMAES}
      , UseGSL{useGSL}
      , UseNLopt{useNLopt}
  {
  }
};

/**
 * @brief GetMinimizers returns a struct containing the bools deciding which
 * minimizers should be used or not
 * @param WhichMinimizer
 * @return
 */
MinimizersToUse GetMinimizers(int WhichMinimizer);

/**
 * @brief WhichMinimizerDefault default value for the Minimizers to use
 */
constexpr int WhichMinimizerDefault = CalcWhichMinimizer();

/**
 * @brief Minimization of the Model
 * Minimizes the given Potential with parameters par and CT-parameters parCT at
 * a given Temperature Temp and writes the solution in the std::vector sol. The
 * Minimization Debugging Options are written in the std::vector Check. The
 * std::vector Start gives the start value for the CMA-ES Minimization.
 */
std::vector<double>
Minimize_gen_all(const std::shared_ptr<Class_Potential_Origin> &modelPointer,
                 const double &Temp,
                 std::vector<double> &Check,
                 const std::vector<double> &start,
                 const int &WhichMinimizer = WhichMinimizerDefault,
                 bool UseMultithreading    = true);

/**
 * @brief The MinimizerStatus enum for the Statusflags of the minimizer
 */
enum class MinimizerStatus
{
  SUCCESS                 = 1,
  NOTVANISHINGATFINALTEMP = -1,
  NOTNLOSTABLE            = -2,
  NUMERICALLYUNSTABLE     = -3,
  BELOWTHRESHOLD          = -4,
  NLOVEVZEROORINF         = -5

};

/**
 * @brief The EWPTReturnType struct
 * Contains the following information
 * @param StatusFlag = SUCCESS : successful
 * @param StatusFlag = NOTVANISHINGATFINALTEMP: v(T=TempEnd) != 0 => vc =
 * v(TempEnd) , Tc = TempEnd
 * @param StatusFlag = NOTNLOSTABLE: Not vacuum stable at T=0 (if TempStart = 0)
 * or v(TempStart) = 0 or v(TempStart) > 255 => Tc = TempEnd, vc = 0
 * @param StatusFlag = NUMERICALLYUNSTABLE: v(T) > 255 during the bisection
 * (this implies the parameter point is numerically unstable) => vc and Tc are
 * the last VEVs encountered
 * @param StatusFlag = BELOWTHRESHOLD: v/T < C_PT during the bisection =>  vc =
 * Last VEV, TC = Last Temp
 * @param vc = critical VEV
 * @param Tc = critical Temperature
 * @param EWMinimum: The broken EW minimum
 */
struct EWPTReturnType
{
  MinimizerStatus StatusFlag;
  double vc;
  double Tc;
  std::vector<double> EWMinimum;
};

/**
 * Uses a bisection method between the Temperature TempStart and TempEnde to
 * find the phase transition in the Model and writes the solution in the vector
 * sol.
 *  @param modelPointer shared_ptr to the parameter point
 *  @param TempStart Low temperature for the starting interval of the bisection
 * method
 *  @param TempEnd High temperature for the starting interval of the bisection
 * method
 *  @param WhichMinimizer Which minimizers should be taken? 1 = CMAES, 2 = GSL,
 * 4 = NLOPT, to use multiple add the numbers
 *  @return The information are returned in a EWPTReturnType struct
 */
EWPTReturnType
PTFinder_gen_all(const std::shared_ptr<Class_Potential_Origin> &modelPointer,
                 const double &TempStart,
                 const double &TempEnd,
                 const int &WhichMinimizer = WhichMinimizerDefault,
                 bool UseMultithreading    = true);

/**
 * @brief Minimize_gen_all_tree_level Minimizes the tree-level potential
 * @param Model Which Model to minimize
 * @param par parameters of the point
 * @param parCT counterterm parameters
 * @param Check Vector to safe the error flags during the minimization
 * @param start Starting point for the minimization
 * @param WhichMinimizer Which minimizers should be taken? 1 = CMAES, 2 = GSL, 4
 * = NLOPT, to use multiple add the numbers
 * @return the global minimum
 */
[[deprecated("Will call Minimize_gen_all_tree_level with GetSMConstants(). "
             "Please use the "
             "detailed overload "
             "to ensure consistent SM constants through all "
             "routines.")]] std::vector<double>
Minimize_gen_all_tree_level(const ModelID::ModelIDs &Model,
                            const std::vector<double> &par,
                            const std::vector<double> &parCT,
                            std::vector<double> &Check,
                            const std::vector<double> &start,
                            int WhichMinimizer     = WhichMinimizerDefault,
                            bool UseMultithreading = true);

/**
 * @brief Minimize_gen_all_tree_level Minimizes the tree-level potential
 * @param Model Which Model to minimize
 * @param par parameters of the point
 * @param parCT counterterm parameters
 * @param SMConstants The SM Constants used for the minimisation
 * @param Check Vector to safe the error flags during the minimization
 * @param start Starting point for the minimization
 * @param WhichMinimizer Which minimizers should be taken? 1 = CMAES, 2 = GSL, 4
 * = NLOPT, to use multiple add the numbers
 * @return the global minimum
 */
std::vector<double>
Minimize_gen_all_tree_level(const ModelID::ModelIDs &Model,
                            const std::vector<double> &par,
                            const std::vector<double> &parCT,
                            const ISMConstants &SMConstants,
                            std::vector<double> &Check,
                            const std::vector<double> &start,
                            int WhichMinimizer     = WhichMinimizerDefault,
                            bool UseMultithreading = true);

/**
 * @brief FindNextLocalMinima finds the local minima from the given starting
 * point at the given temperature
 * @param model Parameter point to minimize
 * @param StartingPoint Starting point from where to look for the next local
 * minima
 * @param temperature Temperature at which to minimize the potential
 * @param WhichMinimizer Which Minimizer should be used? CMAES is not used here
 * as it is only a global minimizer
 * @return A vector containing the local minima. In case different Minima find
 * different minima all solutions are given.
 */
std::vector<std::vector<double>>
FindNextLocalMinima(const std::shared_ptr<Class_Potential_Origin> &model,
                    const std::vector<double> &StartingPoint,
                    const double &temperature,
                    int WhichMinimizer = WhichMinimizerDefault);

/**
 * @brief MinimaDevelopmentWithTemperature calculates the temperature
 * development of several local minima
 * @param model The used model and parameter point
 * @param StartingTemperature The first temperature at which to calculate the
 * local minima from the random points
 * @param FinalTemperature The final temperature for the temperature development
 * @param StepsizeTemperature The stepsize for the temperature development
 * @param RNGRanges The ranges for the RNG generated starting point at
 * StartingTemperature
 * @param seed The seed used for the RNG
 * @param NumberOfStartingPoints How many RNG points should be generated
 * @param WhichMinimizer Which minimizers should be used?
 * @return List of Minima development, each member being a list for one point
 * with the entries <Temperature, Minimum>
 */
std::vector<std::vector<std::pair<double, std::vector<double>>>>
MinimaDevelopmentWithTemperature(
    const std::shared_ptr<Class_Potential_Origin> &model,
    const double &StartingTemperature,
    const double &FinalTemperature,
    const double &StepsizeTemperature,
    const std::vector<std::pair<double, double>> &RNGRanges,
    const std::size_t &seed,
    const std::size_t &NumberOfStartingPoints,
    const int &WhichMinimizer);

} // namespace Minimizer
} // namespace BSMPT

#endif /* MINIMIZER_H_ */
