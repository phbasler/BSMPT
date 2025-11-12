// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/minimizer/LibNLOPT/MinimizeNLOPT.h>
#include <BSMPT/minimizer/MinimizePlane.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>

#include <mutex>
#include <queue>

/**
 *@file
 */
namespace BSMPT
{
namespace Minimizer
{
namespace LibNLOPT
{

double
NLOPTVEff(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
  auto settings = *static_cast<ShareInformationNLOPT *>(data);
  (void)grad;
  auto PotVEV = settings.model.MinimizeOrderVEV(x);
  return settings.model.VEff(PotVEV, settings.Temp);
}

NLOPTReturnType MinimizeUsingNLOPT(const Class_Potential_Origin &model,
                                   const double &Temp)
{
  ShareInformationNLOPT settings(model, Temp);
  std::vector<double> VEV(model.get_nVEV());

  nlopt::opt opt(nlopt::GN_ORIG_DIRECT_L,
                 static_cast<unsigned int>(model.get_nVEV()));
  std::vector<double> LowerBound(model.get_nVEV(), -300),
      UpperBound(model.get_nVEV(), 300);
  for (std::size_t i{0}; i < model.get_nVEV(); ++i)
  {
    if (std::abs(model.get_vevTreeMin(i)) > UpperBound.at(i))
    {
      UpperBound.at(i) = 1.5 * model.get_vevTreeMin(i);
    }
    if (-std::abs(model.get_vevTreeMin(i)) < LowerBound.at(i))
    {
      LowerBound.at(i) = -1.5 * std::abs(model.get_vevTreeMin(i));
    }
  }
  opt.set_lower_bounds(LowerBound);
  opt.set_upper_bounds(UpperBound);

  opt.set_min_objective(NLOPTVEff, &settings);
  opt.set_xtol_rel(1e-4);
  opt.set_maxeval(1000);

  double minf;
  try
  {
    auto result  = opt.optimize(VEV, minf);
    bool Success = (result == nlopt::SUCCESS) or
                   (result == nlopt::FTOL_REACHED) or
                   (result == nlopt::XTOL_REACHED);
    NLOPTReturnType res(VEV, minf, result, Success);
    return res;
  }
  catch (std::exception &e)
  {
    (void)e;
    NLOPTReturnType res(
        std::vector<double>(), 0, nlopt::result::FORCED_STOP, false);
    return res;
  }
}

double NLOPTVEffPlane(const std::vector<double> &x,
                      std::vector<double> &grad,
                      void *data)
{
  (void)grad;
  auto params = *static_cast<PointerContainerMinPlane *>(data);
  std::vector<double> vMinTilde(std::begin(x), std::end(x));
  auto vev = TransformCoordinates(vMinTilde, params);
  return params.modelPointer->VEff(params.modelPointer->MinimizeOrderVEV(vev),
                                   params.Temp);
}

NLOPTReturnType
MinimizePlaneUsingNLOPT(const struct PointerContainerMinPlane &params)
{
  auto dim = params.modelPointer->get_nVEV() - 1;
  nlopt::opt opt(nlopt::GN_ORIG_DIRECT_L, static_cast<unsigned int>(dim));
  std::vector<double> LowerBound(dim, -300), UpperBound(dim, 300);
  for (std::size_t i{0}; i < params.modelPointer->get_nVEV(); ++i)
  {
    if (i < params.Index)
    {
      if (std::abs(params.modelPointer->get_vevTreeMin(i)) > UpperBound.at(i))
      {
        UpperBound.at(i) = 1.5 * params.modelPointer->get_vevTreeMin(i);
      }
      if (-std::abs(params.modelPointer->get_vevTreeMin(i)) < LowerBound.at(i))
      {
        LowerBound.at(i) =
            -1.5 * std::abs(params.modelPointer->get_vevTreeMin(i));
      }
    }
    else if (i > params.Index)
    {
      if (std::abs(params.modelPointer->get_vevTreeMin(i)) >
          UpperBound.at(i - 1))
      {
        UpperBound.at(i - 1) = 1.5 * params.modelPointer->get_vevTreeMin(i);
      }
      if (-std::abs(params.modelPointer->get_vevTreeMin(i)) <
          LowerBound.at(i - 1))
      {
        LowerBound.at(i - 1) =
            -1.5 * std::abs(params.modelPointer->get_vevTreeMin(i));
      }
    }
  }
  opt.set_lower_bounds(LowerBound);
  opt.set_upper_bounds(UpperBound);

  std::vector<double> VEV(params.modelPointer->get_nVEV() - 1);
  auto copy = params;
  opt.set_min_objective(NLOPTVEffPlane, &copy);
  opt.set_xtol_rel(1e-4);
  double minf;
  auto result  = opt.optimize(VEV, minf);
  auto minimum = TransformCoordinates(VEV, params);
  bool Success = (result == nlopt::SUCCESS) or
                 (result == nlopt::FTOL_REACHED) or
                 (result == nlopt::XTOL_REACHED);
  NLOPTReturnType res(minimum, minf, result, Success);
  return res;
}

NLOPTReturnType
FindLocalMinimum(const std::shared_ptr<Class_Potential_Origin> &model,
                 const std::vector<double> &Start,
                 const double &Temp)
{
  nlopt::opt opt(nlopt::LN_COBYLA,
                 static_cast<unsigned int>(model->get_nVEV()));
  ShareInformationNLOPT settings(*model, Temp);
  auto VEV = Start;
  opt.set_min_objective(NLOPTVEff, &settings);
  opt.set_xtol_rel(1e-4);

  double minf;
  auto result  = opt.optimize(VEV, minf);
  bool Success = (result == nlopt::SUCCESS) or
                 (result == nlopt::FTOL_REACHED) or
                 (result == nlopt::XTOL_REACHED);
  NLOPTReturnType res(VEV, minf, result, Success);

  return res;
}

std::pair<std::vector<double>, bool>
NLOPT_SBPLX_Find_Global_Minimum(const Class_Potential_Origin &model,
                                const double &Temp,
                                const int &seed)
{
  std::vector<double> solutions;

  (void)model;
  (void)Temp;
  (void)seed;

  const std::size_t MaxSol{20};

  std::vector<std::vector<double>> saveAllMinima;

  std::size_t dim = model.get_nVEV();

  std::default_random_engine randGen(seed);
  double RNDMax        = 500;
  std::size_t MaxTries = 600; // 600;

  std::queue<std::vector<double>> StartingPoints;
  for (std::size_t i{0}; i < MaxTries; ++i)
  {
    std::vector<double> start(dim);
    for (std::size_t j = 0; j < dim; ++j)
    {
      start.at(j) =
          RNDMax *
          (-1 +
           2 * std::generate_canonical<double,
                                       std::numeric_limits<double>::digits>(
                   randGen));
    }
    StartingPoints.push(start);
  }

  std::atomic<std::size_t> FoundSolutions{0};
  std::vector<std::thread> MinThreads;
  std::mutex WriteResultLock;
  std::vector<std::pair<std::vector<double>, double>> Results;

  auto thread_job =
      [&FoundSolutions, &WriteResultLock](
          std::size_t maxSol,
          std::vector<std::pair<std::vector<double>, double>> &mResults,
          std::queue<std::vector<double>> &mStartingPoints,
          const Class_Potential_Origin &mModel,
          const double &mTemp)
  {
    auto loop_condition = [&]() -> bool
    {
      std::unique_lock<std::mutex> lock(WriteResultLock);
      return FoundSolutions < maxSol and not mStartingPoints.empty();
    };
    while (loop_condition())
    {
      std::vector<double> start;
      {
        std::unique_lock<std::mutex> lock(WriteResultLock);
        start = mStartingPoints.front();
        mStartingPoints.pop();
      }

      auto solution = NLOPT_SBPLX_Find_Global_Minimum(mModel, mTemp, start);
      if (not solution.first.empty())
      {
        std::unique_lock<std::mutex> lock(WriteResultLock);
        ++FoundSolutions;
        mResults.push_back(std::move(solution));
      }
    }
  };

  for (std::size_t i = 0; i < Num_threads; ++i)
  {
    MinThreads.push_back(std::thread(
        [&]() { thread_job(MaxSol, Results, StartingPoints, model, Temp); }));
  }

  for (auto &thr : MinThreads)
  {
    if (thr.joinable())
    {
      std::stringstream ss;
      ss << "Wait for thread" << thr.get_id();
      Logger::Write(LoggingLevel::MinimizerDetailed, ss.str());
      thr.join();
    }
    else
    {
      std::stringstream ss;
      ss << "Thread " << thr.get_id() << " is not joinable";
      Logger::Write(LoggingLevel::MinimizerDetailed, ss.str());
    }
  }

  Logger::Write(LoggingLevel::MinimizerDetailed, "Finished threads");

  if (Results.size() == 0)
  {
    Logger::Write(LoggingLevel::Default,
                  "No solutions found during the GSL minimization at T = " +
                      std::to_string(Temp) + " GeV ");
    return std::make_pair(std::vector<double>{}, false);
  }

  if (Results.size() < MaxSol)
  {
    Logger::Write(LoggingLevel::MinimizerDetailed,
                  "Found " + std::to_string(saveAllMinima.size()) + " of  " +
                      std::to_string(MaxSol) +
                      " solutions at T = " + std::to_string(Temp));
  }

  auto minIter = std::min_element(Results.begin(),
                                  Results.end(),
                                  [dim](const auto &lhs, const auto &rhs)
                                  { return lhs.second <= rhs.second; });

  return {(*minIter).first, true};
}

std::pair<std::vector<double>, double>
NLOPT_SBPLX_Find_Global_Minimum(const Class_Potential_Origin &model,
                                const double &Temp,
                                const std::vector<double> &start)
{
  nlopt::opt opt(nlopt::LN_SBPLX, static_cast<unsigned int>(model.get_nVEV()));
  ShareInformationNLOPT settings(model, Temp);
  auto VEV = start;
  opt.set_min_objective(NLOPTVEff, &settings);
  // opt.set_xtol_rel(1e-4);
  opt.set_ftol_abs(1e-4);

  double minf;
  auto result  = opt.optimize(VEV, minf);
  bool Success = (result == nlopt::SUCCESS) or
                 (result == nlopt::FTOL_REACHED) or
                 (result == nlopt::XTOL_REACHED);

  if (Success)
  {
    return {VEV, minf};
  }
  else
  {
    return {};
  }
  // NLOPTReturnType res(VEV, minf, result, Success);

  // return res;

  // return {};
}

} // namespace LibNLOPT
} // namespace Minimizer
} // namespace BSMPT
