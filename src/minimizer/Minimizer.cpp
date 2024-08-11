// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>     // for FChoose
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
#include <algorithm> // for copy, max
#include <iostream>  // for operator<<, cout, endl
#include <math.h>    // for abs, log10
#include <memory>
#include <random>
#include <thread>
#include <time.h> // for time, NULL
#include <vector>

#include <BSMPT/config.h>

#include <BSMPT/minimizer/MinimizeGSL.h>

#include <exception>

#ifdef cmaes_FOUND
#include <BSMPT/minimizer/LibCMAES/MinimizeLibCMAES.h>
#endif

#ifdef NLopt_FOUND
#include <BSMPT/minimizer/LibNLOPT/MinimizeNLOPT.h>
#endif

/**
 * @file
 * @param WhichMinimizer = 1/2 for single minimizer, sum of those for multiple
 *
 * WhichMinimizer = 1 -> Use CMA-ES Minimizer
 *
 * WhichMinimizer = 2 -> Use the local Minimization of GSL
 *
 *
 *
 *
 *
 *
 * @param Model Which model is considered?
 * @param par Array with parameter point used for set_gen
 * @param parCT Array with Counterterm values used for set_CT_pot_par
 * @param TempStart Starting value of the Temperature for the bisection Method
 * @param TempEnde End value of the Temperature for the bisection Method
 * @param sol std::vector<double> in which the solution of the PTFinder will be
 * stored in the format (T_C,v_C,-+1, single vevs) where +1 means sucessfull
 * find of v_c/T_c > 1 and -1 hit one of the abortion criteria
 *
 */

namespace BSMPT
{
namespace Minimizer
{

MinimizersToUse GetMinimizers(int WhichMinimizer)
{
  bool UseCMAES = (WhichMinimizer % 2 != 0);
  WhichMinimizer /= 2;
  bool UseGSL = (WhichMinimizer % 2 != 0);
  WhichMinimizer /= 2;
  bool UseNLopt = (WhichMinimizer % 2 != 0);

#ifndef cmaes_FOUND
  UseCMAES = false;
#endif

#ifndef NLopt_FOUND
  UseNLopt = false;
#endif

  return MinimizersToUse(UseCMAES, UseGSL, UseNLopt);
}

std::vector<double>
Minimize_gen_all(const std::shared_ptr<Class_Potential_Origin> &modelPointer,
                 const double &Temp,
                 std::vector<double> &Check,
                 const std::vector<double> &start,
                 const int &WhichMinimizer,
                 bool UseMultithreading)
{
  std::vector<double> PotValues;
  std::vector<std::vector<double>> Minima;

  auto UseMinimizer = GetMinimizers(WhichMinimizer);

  bool CheckZero = true; // Check if zero is the global minimum explicitly
  if (CheckZero)
  {
    PotValues.push_back(modelPointer->VEff(
        std::vector<double>(modelPointer->get_NHiggs(), 0), Temp));
    Minima.push_back(std::vector<double>(modelPointer->get_nVEV(), 0));
  }

  if (modelPointer->get_nVEV() <= 2)
  {
    UseMinimizer.UseCMAES = false;
    UseMinimizer.UseGSL   = true;
  }

  std::vector<double> solGSLMin, solGSLMinPot;

  bool gslMinSuc = false;
  std::thread thread_GSL;
  if (UseMinimizer.UseGSL)
  {

    if (UseMinimizer.UseCMAES or UseMinimizer.UseNLopt)
    {
      if (UseMultithreading)
      {
        thread_GSL = std::thread(
            [&solGSLMin, &gslMinSuc, &modelPointer, &Temp]()
            {
              std::tie(solGSLMin, gslMinSuc) = GSL_Minimize_gen_all(
                  *modelPointer,
                  Temp,
                  5); // If additionally CMAES is minimising
                      // GSL does not need as much solutions
            });
      }
      else
      {
        std::tie(solGSLMin, gslMinSuc) =
            GSL_Minimize_gen_all(*modelPointer, Temp, 5, UseMultithreading);
      }
    }
    else
    {
      std::size_t MaxSol = 50;
      if (UseMultithreading)
      {
        thread_GSL = std::thread(
            [&solGSLMin, &gslMinSuc, &modelPointer, &Temp, &MaxSol]()
            {
              std::tie(solGSLMin, gslMinSuc) =
                  GSL_Minimize_gen_all(*modelPointer, Temp, 5, MaxSol);
            });
      }
      else
      {
        std::tie(solGSLMin, gslMinSuc) = GSL_Minimize_gen_all(
            *modelPointer, Temp, 5, MaxSol, UseMultithreading);
      }
    }
  }
#ifdef cmaes_FOUND
  std::thread thread_CMAES;
  LibCMAES::LibCMAESReturn LibCMAES;
  if (UseMinimizer.UseCMAES)
  {
    if (UseMultithreading)
    {
      thread_CMAES = std::thread(
          [&modelPointer, &Temp, &start, &LibCMAES]() {
            LibCMAES = LibCMAES::min_cmaes_gen_all(*modelPointer, Temp, start);
          });
    }
    else
    {
      LibCMAES = LibCMAES::min_cmaes_gen_all(*modelPointer, Temp, start);
    }
  }
#else
  (void)start;
#endif

#ifdef NLopt_FOUND
  LibNLOPT::NLOPTReturnType NLOPTResult;
  std::thread thread_NLopt;
  if (UseMinimizer.UseNLopt)
  {
    if (UseMultithreading)
    {
      thread_NLopt = std::thread(
          [&NLOPTResult, &modelPointer, &Temp]()
          { NLOPTResult = LibNLOPT::MinimizeUsingNLOPT(*modelPointer, Temp); });
    }
    else
    {
      NLOPTResult = LibNLOPT::MinimizeUsingNLOPT(*modelPointer, Temp);
    }
  }
#endif

#ifdef NLopt_FOUND
  if (UseMultithreading and thread_NLopt.joinable())
  {
    Logger::Write(LoggingLevel::MinimizerDetailed, "Waiting for NLopt thread");
    thread_NLopt.join();
  }

  if (NLOPTResult.Success)
  {
    PotValues.push_back(NLOPTResult.PotVal);
    Minima.push_back(NLOPTResult.Minimum);
    std::stringstream ss;
    ss << "NLopt candidate at T = " << Temp << " :  " << NLOPTResult.Minimum
       << " with potential value " << NLOPTResult.PotVal << std::endl;
    Logger::Write(LoggingLevel::MinimizerDetailed, ss.str());
  }
#endif

#ifdef cmaes_FOUND
  if (UseMultithreading and thread_CMAES.joinable())
  {
    Logger::Write(LoggingLevel::MinimizerDetailed, "Waiting for CMAES Thread");
    thread_CMAES.join();
    auto errC          = LibCMAES.CMAESStatus;
    auto solCMAES      = LibCMAES.result;
    auto solCMAESPotIn = modelPointer->MinimizeOrderVEV(solCMAES);
    PotValues.push_back(modelPointer->VEff(solCMAESPotIn, Temp));
    Minima.push_back(solCMAES);
    Check.push_back(errC);

    std::stringstream ss;
    ss << "CMAES candidate at T = " << Temp << " : " << solCMAES
       << " with potential value = " << PotValues.at(PotValues.size() - 1)
       << std::endl;
    Logger::Write(LoggingLevel::MinimizerDetailed, ss.str());
  }
#endif

  if (UseMultithreading and thread_GSL.joinable())
  {
    Logger::Write(LoggingLevel::MinimizerDetailed, "Waiting for GSL thread");
    thread_GSL.join();
  }

  if (gslMinSuc)
  {
    solGSLMinPot = modelPointer->MinimizeOrderVEV(solGSLMin);
    PotValues.push_back(modelPointer->VEff(solGSLMinPot, Temp));
    Minima.push_back(solGSLMin);

    std::stringstream ss;
    ss << "GSL found a minimum at T = " << Temp << ": (" << solGSLMin
       << ") with Potential value = " << PotValues.at(PotValues.size() - 1)
       << std::endl;
    Logger::Write(LoggingLevel::MinimizerDetailed, ss.str());
  }

  std::size_t minIndex = 0;
  for (std::size_t i = 1; i < PotValues.size(); i++)
  {
    if (PotValues.at(i) < PotValues.at(minIndex)) minIndex = i;
  }

  auto sol   = Minima.at(minIndex);
  auto EWVEV = modelPointer->EWSBVEV(modelPointer->MinimizeOrderVEV(sol));
  if (EWVEV <= 0.5) modelPointer->SetEWVEVZero(sol);

  solGSLMin.clear();
  if (UseMinimizer.UseGSL and gslMinSuc)
    Check.push_back(1);
  else
    Check.push_back(-1);

  return sol;
}

EWPTReturnType
PTFinder_gen_all(const std::shared_ptr<Class_Potential_Origin> &modelPointer,
                 const double &TempStart,
                 const double &TempEnd,
                 const int &WhichMinimizer,
                 bool UseMultithreading)
{

  EWPTReturnType result;

  std::size_t dim = modelPointer->get_nVEV();

  double vStart, vEnd, vMitte;
  std::vector<double> solStart, solMitte, solEnd;
  std::vector<double> solStartPot, solMittePot, solEndPot;
  std::vector<double> checkStart, checkMitte, checkEnde;
  std::vector<double> startStart, startMitte, startEnde;
  double Distance = 1e-2;
  double TA       = TempStart;
  double TM;
  double TE = TempEnd;

  std::vector<double> pSol(dim + 1);

  for (std::size_t k = 0; k < dim; k++)
    startEnde.push_back(modelPointer->get_vevTreeMin(k));
  solEnd    = Minimize_gen_all(modelPointer,
                            TempEnd,
                            checkEnde,
                            startEnde,
                            WhichMinimizer,
                            UseMultithreading);
  solEndPot = modelPointer->MinimizeOrderVEV(solEnd);
  vEnd      = modelPointer->EWSBVEV(solEndPot);

  if (vEnd > C_threshold)
  {
    result.Tc         = TempEnd;
    result.vc         = vEnd;
    result.StatusFlag = MinimizerStatus::NOTVANISHINGATFINALTEMP;
    result.EWMinimum  = solEnd;
    return result;
  }

  for (std::size_t k = 0; k < dim; k++)
    startStart.push_back(modelPointer->get_vevTreeMin(k));
  solStart    = Minimize_gen_all(modelPointer,
                              TempStart,
                              checkStart,
                              startStart,
                              WhichMinimizer,
                              UseMultithreading);
  solStartPot = modelPointer->MinimizeOrderVEV(solStart);
  vStart      = modelPointer->EWSBVEV(solStartPot);

  if (vStart <= C_threshold or vStart >= 255.0)
  {
    result.Tc         = TempEnd;
    result.vc         = 0;
    result.StatusFlag = MinimizerStatus::NLOVEVZEROORINF;
    result.EWMinimum  = std::vector<double>(dim, 0);
    return result;
  }

  bool SurviveNLO = modelPointer->CheckNLOVEV(solStart);

  if (not SurviveNLO and TempStart == 0)
  {
    result.Tc         = TempEnd;
    result.vc         = vStart;
    result.StatusFlag = MinimizerStatus::NOTNLOSTABLE;
    result.EWMinimum  = solStart;
    return result;
  }

  for (std::size_t k = 0; k < dim; k++)
    startMitte.push_back(modelPointer->get_vevTreeMin(k));
  do
  {
    TM = 0.5 * (TA + TE);
    solMitte.clear();
    checkMitte.clear();
    solMittePot.clear();
    solMitte    = Minimize_gen_all(modelPointer,
                                TM,
                                checkMitte,
                                startMitte,
                                WhichMinimizer,
                                UseMultithreading);
    solMittePot = modelPointer->MinimizeOrderVEV(solMitte);
    vMitte      = modelPointer->EWSBVEV(solMittePot);

    if (vMitte >= 255.0)
    {
      result.Tc         = TM;
      result.vc         = vMitte;
      result.StatusFlag = MinimizerStatus::NUMERICALLYUNSTABLE;
      result.EWMinimum  = solMitte;
      return result;
    }
    else if (vMitte != 0 and vMitte / TM < C_PT)
    {
      result.Tc         = TM;
      result.vc         = vMitte;
      result.StatusFlag = MinimizerStatus::BELOWTHRESHOLD;
      result.EWMinimum  = solMitte;
      return result;
    }
    if (vMitte >= Distance)
    {
      TA = TM;
      startMitte.clear();
      for (std::size_t k = 0; k < dim; k++)
        pSol[k + 1] = solMitte.at(k);
      pSol[0] = vMitte;
      for (std::size_t k = 0; k < dim; k++)
        startMitte.push_back(pSol[k + 1]);
      vStart = vMitte;
    }
    else
    {
      TE = TM;
    }

  } while (std::abs(TE - TA) > Distance);

  result.Tc         = TA;
  result.vc         = pSol[0];
  result.StatusFlag = MinimizerStatus::SUCCESS;
  std::vector<double> SolMin;
  result.EWMinimum.clear();
  for (std::size_t k = 1; k <= dim; k++)
    result.EWMinimum.push_back(pSol[k]);
  return result;
}

std::vector<double>
Minimize_gen_all_tree_level(const ModelID::ModelIDs &Model,
                            const std::vector<double> &par,
                            const std::vector<double> &parCT,
                            std::vector<double> &Check,
                            const std::vector<double> &start,
                            int WhichMinimizer,
                            bool UseMultithreading)
{
  return Minimize_gen_all_tree_level(Model,
                                     par,
                                     parCT,
                                     GetSMConstants(),
                                     Check,
                                     start,
                                     WhichMinimizer,
                                     UseMultithreading);
}

std::vector<double>
Minimize_gen_all_tree_level(const ModelID::ModelIDs &Model,
                            const std::vector<double> &par,
                            const std::vector<double> &parCT,
                            const ISMConstants &SMConstants,
                            std::vector<double> &Check,
                            const std::vector<double> &start,
                            int WhichMinimizer,
                            bool UseMultithreading)
{
  std::shared_ptr<Class_Potential_Origin> modelPointer =
      ModelID::FChoose(Model, SMConstants);
  modelPointer->set_All(par, parCT);
  modelPointer->SetUseTreeLevel(true);
  auto sol = Minimize_gen_all(
      modelPointer, 0, Check, start, WhichMinimizer, UseMultithreading);
  modelPointer->SetUseTreeLevel(false);
  return sol;
}

std::vector<std::vector<double>>
FindNextLocalMinima(const std::shared_ptr<Class_Potential_Origin> &model,
                    const std::vector<double> &StartingPoint,
                    const double &temperature,
                    int WhichMinimizer)
{
  const auto UseMinimizer = GetMinimizers(WhichMinimizer);

  std::vector<std::vector<double>> Minima;

  if (UseMinimizer.UseGSL)
  {
    std::vector<double> GSLSolution;
    std::size_t tries{0}, MaxTries{600};
    int status;
    GSL_params params(*model, temperature);
    do
    {
      GSLSolution.clear();
      status = GSL_Minimize_From_S_gen_all(params, GSLSolution, StartingPoint);
      tries++;
    } while (status != GSL_SUCCESS and tries < MaxTries);
    if (status == GSL_SUCCESS)
    {
      Minima.push_back(GSLSolution);
    }
  }

#ifdef NLopt_FOUND

  if (UseMinimizer.UseNLopt)
  {
    auto NLOPTres =
        LibNLOPT::FindLocalMinimum(model, StartingPoint, temperature);
    if (NLOPTres.Success)
    {
      Minima.push_back(NLOPTres.Minimum);
    }
  }
#endif

  return Minima;
}

std::vector<std::vector<std::pair<double, std::vector<double>>>>
MinimaDevelopmentWithTemperature(
    const std::shared_ptr<Class_Potential_Origin> &model,
    const double &StartingTemperature,
    const double &FinalTemperature,
    const double &StepsizeTemperature,
    const std::vector<std::pair<double, double>> &RNGRanges,
    const std::size_t &seed,
    const std::size_t &NumberOfStartingPoints,
    const int &WhichMinimizer)
{
  using MinimaDevelopmentType =
      std::vector<std::pair<double, std::vector<double>>>;
  std::vector<MinimaDevelopmentType> res;
  std::default_random_engine randGen(seed);

  std::vector<std::vector<double>> StartingPoints;
  for (std::size_t i{0}; i < NumberOfStartingPoints; ++i)
  {
    std::vector<double> Point;
    for (const auto &el : RNGRanges)
    {
      Point.push_back(
          el.first +
          (el.second - el.first) *
              std::generate_canonical<double,
                                      std::numeric_limits<double>::digits>(
                  randGen));
    }
    StartingPoints.push_back(Point);
  }

  for (const auto &StartingPoint : StartingPoints)
  {
    auto LocalMinima = FindNextLocalMinima(
        model, StartingPoint, StartingTemperature, WhichMinimizer);
    for (const auto &el : LocalMinima)
    {
      auto Min = std::make_pair(StartingTemperature, el);
      res.push_back(MinimaDevelopmentType{Min});
    }
  }

  auto StoppingCriteria = [&](const double &Temp)
  {
    double epsilon = std::abs((FinalTemperature - StartingTemperature) /
                              StepsizeTemperature) *
                     1e-2; // because 1+1 = 2.0000000000000001
    if (StartingTemperature < FinalTemperature)
    {
      return Temp < FinalTemperature + epsilon;
    }
    else
    {
      return Temp > FinalTemperature - epsilon;
    }
  };

  for (double Temp = StartingTemperature + StepsizeTemperature;
       StoppingCriteria(Temp);
       Temp += StepsizeTemperature)
  {
    auto resOld = std::move(res);
    res.clear();
    for (const auto &Point : resOld)
    {
      auto LatestMinima = Point.at(Point.size() - 1);
      auto NextMinima =
          FindNextLocalMinima(model, LatestMinima.second, Temp, WhichMinimizer);
      for (const auto &NM : NextMinima)
      {
        auto base{Point};
        base.push_back(std::make_pair(Temp, NM));
        res.push_back(base);
      }
    }
  }

  return res;
}

} // namespace Minimizer
} // namespace BSMPT
