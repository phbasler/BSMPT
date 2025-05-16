// SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file bounce solution calculator class
 */

#include <BSMPT/bounce_solution/bounce_solution.h>

namespace BSMPT
{

BounceSolution::BounceSolution(
    const std::shared_ptr<Class_Potential_Origin> &pointer_in)
{
  modelPointer = pointer_in;
}

BounceSolution::BounceSolution(
    const std::shared_ptr<Class_Potential_Origin> &pointer_in,
    const std::shared_ptr<MinimumTracer> &MinTracer_in,
    const CoexPhases &phase_pair_in,
    const double &UserDefined_vwall_in,
    const double &UserDefined_epsturb_in,
    const int &MaxPathIntegrations_in,
    const size_t &NumberOfInitialScanTemperatures_in,
    const std::vector<Eigen::MatrixXd> &GroupElements_in,
    const int &UserDefined_PNLO_scaling_in)
{
  modelPointer = pointer_in;
  MinTracer    = MinTracer_in;

  phase_pair = phase_pair_in;
  Tc         = phase_pair.crit_temp;
  Tm         = phase_pair.T_low;

  UserDefined_vwall               = UserDefined_vwall_in;
  epsturb                         = UserDefined_epsturb_in;
  pnlo_scaling                    = UserDefined_PNLO_scaling_in;
  MaxPathIntegrations             = MaxPathIntegrations_in;
  NumberOfInitialScanTemperatures = NumberOfInitialScanTemperatures_in;
  InitializeGstarProfile();
  GroupElements = GroupElements_in;

  if (Tc > 0)
  // Calculate which of the VEV has the best change of tunneling
  {
    CalculateOptimalDiscreteSymmetry();
    BounceSolution::GWInitialScan();
  }
}

BounceSolution::BounceSolution(
    const std::shared_ptr<Class_Potential_Origin> &pointer_in,
    const std::shared_ptr<MinimumTracer> &MinTracer_in,
    const CoexPhases &phase_pair_in,
    const double &UserDefined_vwall_in,
    const double &UserDefined_epsturb_in,
    const int &MaxPathIntegrations_in,
    const size_t &NumberOfInitialScanTemperatures_in,
    const int &UserDefined_PNLO_scaling_in)
    : BounceSolution(pointer_in,
                     MinTracer_in,
                     phase_pair_in,
                     UserDefined_vwall_in,
                     UserDefined_epsturb_in,
                     MaxPathIntegrations_in,
                     NumberOfInitialScanTemperatures_in,
                     {Eigen::MatrixXd::Identity(pointer_in->get_nVEV(),
                                                pointer_in->get_nVEV())},
                     UserDefined_PNLO_scaling_in)
{
}

void BounceSolution::CalculateOptimalDiscreteSymmetry()
{
  std::stringstream ss;
  ss << "Calculating optimal symmetry\n";

  Eigen::VectorXd TrueVacuum = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
      phase_pair.true_phase.Get(Tc).point.data(),
      phase_pair.true_phase.Get(Tc).point.size());
  std::vector<double> FalseVacuum = phase_pair.false_phase.Get(Tc).point;

  // Set the optimal discrete symmetro to the identity matrix
  OptimalDiscreteSymmetry =
      Eigen::MatrixXd::Identity(FalseVacuum.size(), FalseVacuum.size());

  double MaximumDistance = 1e100;
  for (const auto &GroupElement : GroupElements)
  {
    Eigen::VectorXd CandidateTrueVacuum = GroupElement * TrueVacuum;
    std::vector<double> DeltaVEV        = FalseVacuum;
    for (std::size_t el = 0; el < FalseVacuum.size(); el++)
      DeltaVEV.push_back(FalseVacuum[el] - TrueVacuum[el]);

    if (L2NormVector(DeltaVEV) < MaximumDistance)
    {
      MaximumDistance         = L2NormVector(DeltaVEV);
      OptimalDiscreteSymmetry = GroupElement;
    }
  }

  ss << "Optimal symmetry is \n\n" << OptimalDiscreteSymmetry << "\n";
  Logger::Write(LoggingLevel::BounceDetailed, ss.str());
}

void BounceSolution::SetAndCalculateGWParameters(
    const TransitionTemperature &which_transition_temp_in)
{
  which_transition_temp = which_transition_temp_in;
  CalcTransitionTemp();
  Logger::Write(LoggingLevel::TransitionDetailed,
                "Calculate PT strength and inverse time scale at the chosen "
                "transition temperature.");
  CalculatePTStrength();
  CalculateInvTimeScale();

  CalculateReheatingTemp();
}

std::vector<double> BounceSolution::TransformIntoOptimalDiscreteSymmetry(
    const std::vector<double> &vev)
{
  std::vector<double> result(vev.size(), 0);
  // Transform the vector
  for (std::size_t i = 0; i < vev.size(); i++)
    for (std::size_t j = 0; j < vev.size(); j++)
      result[i] += OptimalDiscreteSymmetry(i, j) * vev[j];
  return result;
}

void BounceSolution::GWInitialScan()
{
  if (Tc < 0)
  {
    // Transition is never viable
    return;
  }

  double last_action = -1;
  double dT = (Tc - phase_pair.T_low) / NumberOfInitialScanTemperatures;
  std::vector<double> TrueVacuum, FalseVacuum, last_TrueVacuum,
      last_FalseVacuum;
  std::vector<std::vector<double>> last_path, path;

  for (double T = Tc - dT; T >= phase_pair.T_low + dT; T -= dT)
  {
    Logger::Write(LoggingLevel::BounceDetailed, "T = " + std::to_string(T));

    // Check if transition is energetically viable
    if (phase_pair.true_phase.Get(T).potential >=
        phase_pair.false_phase.Get(T).potential)
      continue;

    TrueVacuum = TransformIntoOptimalDiscreteSymmetry(
        phase_pair.true_phase.Get(T).point);
    FalseVacuum = phase_pair.false_phase.Get(T).point;
    std::function<double(std::vector<double>)> V = [&](std::vector<double> vev)
    {
      // Potential wrapper
      return modelPointer->VEff(modelPointer->MinimizeOrderVEV(vev), T);
    };
    if (last_action < 0)
    {
      path = {TrueVacuum, FalseVacuum};
    }
    else
    {
      path = MinTracer->WarpPath(last_path,
                                 last_TrueVacuum,
                                 last_FalseVacuum,
                                 TrueVacuum,
                                 FalseVacuum);
    }
    BounceActionInt bc(
        path, TrueVacuum, FalseVacuum, V, T, MaxPathIntegrations);
    bc.CalculateAction();

    last_path        = bc.Path;
    last_TrueVacuum  = bc.TrueVacuum;
    last_FalseVacuum = bc.FalseVacuum;

    // Comment this is you want dumb paths!!
    last_action = bc.Action;
    if (bc.Action / T > 0)
    {
      SolutionList.insert(std::upper_bound(SolutionList.begin(),
                                           SolutionList.end(),
                                           bc,
                                           [](const BounceActionInt &a,
                                              const BounceActionInt &b)
                                           { return a.T < b.T; }),
                          bc);
    }

    if (bc.Action / T < 40 and bc.Action > 0) break;
  }
  GWSecondaryScan();
}

void BounceSolution::CalculateActionAt(double T, bool smart)
{
  // Action outside allowed range
  if (T < Tm or T > Tc) return;
  Logger::Write(LoggingLevel::BounceDetailed, " T = " + std::to_string(T));
  // Find the closest solution to our goal temperature
  if (SolutionList.size() > 0)
  {
    auto it =
        std::min_element(SolutionList.begin(),
                         SolutionList.end(),
                         [T](const BounceActionInt &a, const BounceActionInt &b)
                         { return std::abs(T - a.T) < std::abs(T - b.T); });
    BounceActionInt Nearest_bc = *it;

    if (abs(Nearest_bc.T - T) < 0.001) return;

    // Check if transition is energetically viable
    if (phase_pair.true_phase.Get(T).potential >=
        phase_pair.false_phase.Get(T).potential)
      return;

    std::vector<double> TrueVacuum = TransformIntoOptimalDiscreteSymmetry(
        phase_pair.true_phase.Get(T).point);
    std::vector<double> FalseVacuum = phase_pair.false_phase.Get(T).point;
    std::function<double(std::vector<double>)> V = [&](std::vector<double> vev)
    {
      // Potential wrapper
      return modelPointer->VEff(modelPointer->MinimizeOrderVEV(vev), T);
    };
    std::vector<std::vector<double>> path;

    if (smart)
      path = MinTracer->WarpPath(Nearest_bc.Path,
                                 Nearest_bc.TrueVacuum,
                                 Nearest_bc.FalseVacuum,
                                 TrueVacuum,
                                 FalseVacuum);
    else
      path = {TrueVacuum, FalseVacuum};

    BounceActionInt bc(
        path, TrueVacuum, FalseVacuum, V, T, MaxPathIntegrations);
    bc.CalculateAction();
    if (bc.Action / T > 0)
    {
      SolutionList.insert(std::upper_bound(SolutionList.begin(),
                                           SolutionList.end(),
                                           bc,
                                           [](const BounceActionInt &a,
                                              const BounceActionInt &b)
                                           { return a.T < b.T; }),
                          bc);
    }
  }
  else
  {
    std::vector<double> TrueVacuum = TransformIntoOptimalDiscreteSymmetry(
        phase_pair.true_phase.Get(T).point);
    std::vector<double> FalseVacuum = phase_pair.false_phase.Get(T).point;
    std::function<double(std::vector<double>)> V = [&](std::vector<double> vev)
    {
      // Potential wrapper
      return modelPointer->VEff(modelPointer->MinimizeOrderVEV(vev), T);
    };
    std::vector<std::vector<double>> path = {TrueVacuum, FalseVacuum};

    BounceActionInt bc(
        path, TrueVacuum, FalseVacuum, V, T, MaxPathIntegrations);
    bc.CalculateAction();
    if (bc.Action / T > 0)
    {
      SolutionList.push_back(bc);
    }
  }
}

void BounceSolution::GWSecondaryScan()
{
  if (SolutionList.size() == 0)
  {
    Logger::Write(LoggingLevel::BounceDetailed,
                  "No solution was found during the initial scan.\n Abort!\n");
    // Failed to find the solution
    status_bounce_sol = StatusGW::Failure;
    return;
  }
  else if (SolutionList.size() == 1)
  {
    Logger::Write(
        LoggingLevel::BounceDetailed,
        "Only one solution was found. Searching near what was found\n");
    CalculateActionAt(SolutionList[0].T -
                      (SolutionList[0].T - phase_pair.T_low) / 10.);
    if (SolutionList.size() == 1)
    {
      // Failed to find the solution
      SetBounceSol();
      return;
    };
  }

  GWScanTowardsLowAction();
  GWScanTowardsHighAction();

  std::size_t NumOfSol = SolutionList.size();

  for (int i = 0; i < 2; i++)
  {
    std::vector<double> nextTList;

    for (std::size_t sol = 0; sol < SolutionList.size() - 1; sol++)
    {
      double t1 = SolutionList[sol].T;
      double t2 = SolutionList[sol + 1].T;
      double s1 = SolutionList[sol].Action / SolutionList[sol].T;
      double s2 = SolutionList[sol + 1].Action / SolutionList[sol + 1].T;

      // Outside our range
      if (s1 > 200) continue;
      if (s2 < 50) continue;

      double NumOfT = ceil((s2 - s1) / 20);

      double dT = (SolutionList[sol + 1].T - SolutionList[sol].T) / NumOfT;

      double newT;
      double ActionProjection;
      for (int it_T = 1; it_T < NumOfT; it_T++)
      {
        newT             = SolutionList[sol].T + it_T * dT;
        ActionProjection = (newT - t1) * (s2 - s1) / (t2 - t1) + s1;
        // Action expected to fall within the wanted range
        if (ActionProjection < 400 and ActionProjection > 0)
          nextTList.push_back(newT);
      }
    }

    for (auto j : nextTList)
    {
      CalculateActionAt(j);
    }

    if (NumOfSol == SolutionList.size()) break;
    NumOfSol = SolutionList.size();
  }
  GWScanTowardsLowAction();
  GWScanTowardsHighAction();
  SetBounceSol();
}

void BounceSolution::GWScanTowardsHighAction()
{
  for (int i = 0;
       i <=
       1. + (200 - SolutionList.back().Action / SolutionList.back().T) / 10.;
       i++)
  {
    if (SolutionList.back().Action / SolutionList.back().T > 200) break;

    // Try to calculate until s3/T = 200
    double t1            = SolutionList[SolutionList.size() - 2].T;
    double t2            = SolutionList[SolutionList.size() - 1].T;
    double s1            = SolutionList[SolutionList.size() - 2].Action / t1;
    double s2            = SolutionList[SolutionList.size() - 1].Action / t2;
    double goal          = s2 + 10;
    std::size_t NumOfSol = SolutionList.size();

    double T = ((s1 - goal) * t2 - (s2 - goal) * t1) / (s1 - s2);

    // Action is not monotonic
    if (T < t2) return;

    if (T > this->Tc) // Action is flat
    {
      if (t2 + (t2 - t1) > this->Tc) return; // Close already
      CalculateActionAt(t2 + (t2 - t1));
      i--;
    }
    else
    {
      CalculateActionAt(T);
      if (NumOfSol == SolutionList.size())
        CalculateActionAt(
            T, false); // No solution was found. Calculate using dumb method
    }
    if (NumOfSol == SolutionList.size()) return; // No solution was found. Abort
  }
}

void BounceSolution::GWScanTowardsLowAction()
{
  for (int i = 0;
       i <=
       1. + (SolutionList.front().Action / SolutionList.front().T - 50) / 10.;
       i++)
  {
    if (SolutionList.front().Action / SolutionList.front().T < 50) break;
    // Try to calculate it at s3/T = 300
    double t1            = SolutionList[0].T;
    double t2            = SolutionList[1].T;
    double s1            = SolutionList[0].Action / t1;
    double s2            = SolutionList[1].Action / t2;
    double goal          = s1 - 10;
    std::size_t NumOfSol = SolutionList.size();

    double T = ((s1 - goal) * t2 - (s2 - goal) * t1) / (s1 - s2);

    // Action is not monotonic
    if (T > t1) return;

    if (T < this->Tm) // Action is flat
    {
      if (t1 - (t2 - t1) < this->Tm) return; // Close already
      CalculateActionAt(t1 - (t2 - t1));
      i--;
    }
    else
    {
      CalculateActionAt(T);
      if (NumOfSol == SolutionList.size())
        CalculateActionAt(
            T, false); // No solution was found. Calculate using dumb method
    }
    if (NumOfSol == SolutionList.size()) return; // No solution was found. Abort
  }
}

void BounceSolution::SetBounceSol()
{
  std::vector<double> list_T, list_S3, list_S3_T, list_140;
  std::stringstream ss;
  ss << "------------ Solution list ------------\n";
  for (auto sol : SolutionList)
  {
    ss << std::setprecision(10) << "{" << sol.T << ",\t" << sol.Action << ",\t"
       << sol.Action / sol.T << "},\n";
    list_T.push_back(sol.T);
    list_S3.push_back(sol.Action);
    list_S3_T.push_back(abs(log(sol.Action / sol.T)));
    list_140.push_back(log(140));
  }
  ss << "---------------------------------------\n";
  Logger::Write(LoggingLevel::BounceDetailed, ss.str());

  if (SolutionList.size() > 1)
  {
    // Plot log of action/T as function of T
    ss.clear();
    AsciiPlotter plotter("S_3/T", 100, 35);
    plotter.addPlot(list_T, list_140, " ", '.');
    plotter.addPlot(list_T, list_S3_T, " ", '*');
    plotter.xlabel("T (GeV)");
    plotter.ylabel("log(S_3/T)");
    plotter.show(ss);
    Logger::Write(LoggingLevel::BounceDetailed, ss.str());
  }

  if (SolutionList.size() < 4)
  {
    status_bounce_sol = StatusGW::Failure;
    Logger::Write(
        LoggingLevel::BounceDetailed,
        "There are not enough points to calculate the path. Abort.\n");
    return; // Not enough points to calculate path
  }

  S3ofT_spline.set_boundary(
      tk::spline::not_a_knot, 0.0, tk::spline::not_a_knot, 0.0);
  S3ofT_spline.set_points(list_T, list_S3);

  InitializedVSpline(); // If there are 4 solutions we can construct a spline

  status_bounce_sol = StatusGW::Success;
}

double BounceSolution::GetWallVelocity() const
{
  return vwall;
}

double BounceSolution::GetChapmanJougetVelocity() const
{
  return vCJ;
}

double BounceSolution::GetSoundSpeedFalse() const
{
  return Csound_false;
}

double BounceSolution::GetSoundSpeedTrue() const
{
  return Csound_true;
}

double BounceSolution::GetEpsTurb() const
{
  return epsturb;
}

void BounceSolution::SetGstar(const double &gstar_in)
{
  gstar = gstar_in;
}

void BounceSolution::InitializeGstarProfile()
{
  this->SetGstar(CalcGstarPureRad());
  GstarProfileLowT.set_boundary(
      tk::spline::not_a_knot, 0.0, tk::spline::not_a_knot, 0.0);
  GstarProfileLowT.set_points(TGstarLowT, GstarLowT);
  GstarProfileHighT.set_boundary(
      tk::spline::not_a_knot, 0.0, tk::spline::not_a_knot, 0.0);
  GstarProfileHighT.set_points(TGstarHighT, GstarHighT);
}

void BounceSolution::ConstructSplineVofT(Phase &phase, tk::spline &spline)
{
  std::vector<double> T_list, V_list;
  // Extract knots from phase
  for (const auto &m : phase.MinimumPhaseVector)
  {
    T_list.push_back(m.temp);
    V_list.push_back(m.potential);
  }
  // Consctruct spline with not-a-know conditions
  spline.set_boundary(tk::spline::not_a_knot, 0.0, tk::spline::not_a_knot, 0.0);
  spline.set_points(T_list, V_list);
}

void BounceSolution::InitializedVSpline()
{
  // Construct V(phi_false,T) spline
  ConstructSplineVofT(phase_pair.false_phase, FalsePhaseVSpline);

  // Construct V(phi_true,T) spline
  ConstructSplineVofT(phase_pair.true_phase, TruePhaseVSpline);
}

double BounceSolution::GetGstar(const double &T) const
{
  const double TinMeV =
      T * 1000.; // Multiplied by 1000 because the fit was done in MeV
  const double TQCD = TGstarLowT.back(); // QCD transition = 214 MeV
  if (TinMeV < TGstarLowT.front())
    return GstarLowT.front(); // Set to \f$ N_\nu \f$
  if (TinMeV > TGstarHighT.back()) return gstar;
  if (TinMeV < TQCD) return GstarProfileLowT(TinMeV);

  return pow(TinMeV / TQCD,
             (log(gstar / GstarHighT.back()) /
              log(TGstarHighT.back() / TQCD))) *
         GstarProfileHighT(TinMeV);
}

double BounceSolution::GetGstar()
{
  return this->CalcGstarPureRad();
}

void BounceSolution::SetCriticalTemp(const double &T_in)
{
  Tc = T_in;
}

double BounceSolution::GetCriticalTemp() const
{
  return Tc;
}

void BounceSolution::SetStoredTemp(const double &T_in)
{
  store_Temp = T_in;
}

double BounceSolution::GetStoredTemp() const
{
  return store_Temp;
}

double BounceSolution::GetNucleationTemp() const
{
  return Tnucl;
}

double BounceSolution::GetNucleationTempApprox() const
{
  return Tnucl_approx;
}

double BounceSolution::GetPercolationTemp() const
{
  return Tperc;
}

double BounceSolution::GetCompletionTemp() const
{
  return Tcompl;
}

double BounceSolution::GetTransitionTemp() const
{
  return Tstar;
}

double BounceSolution::GetReheatingTemp() const
{
  return Treh;
}

void BounceSolution::CalcTransitionTemp()
{
  if (which_transition_temp == TransitionTemperature::NotSet)
  {
    Logger::Write(
        LoggingLevel::TransitionDetailed,
        "'which_transition_temp' not set. Default to percolation temperature");
    which_transition_temp = TransitionTemperature::Percolation;
  }
  // Calculate all temperatures
  CalculateNucleationTempApprox();
  CalculateNucleationTemp();
  CalculatePercolationTemp();
  CalculateCompletionTemp();
  if (status_nucl_approx == BSMPT::StatusTemperature::Success and
      which_transition_temp == TransitionTemperature::ApproxNucleation)
  {
    Tstar = GetNucleationTempApprox();
    Logger::Write(
        LoggingLevel::TransitionDetailed,
        "Approximate nucleation temperature T = " + std::to_string(Tstar) +
            " chosen as transition temperature.");
  }
  else if (status_nucl == BSMPT::StatusTemperature::Success and
           which_transition_temp == TransitionTemperature::Nucleation)
  {
    Tstar = GetNucleationTemp();
    Logger::Write(LoggingLevel::TransitionDetailed,
                  "Nucleation temperature T = " + std::to_string(Tstar) +
                      " chosen as transition temperature.");
  }
  else if (status_perc == BSMPT::StatusTemperature::Success and
           which_transition_temp == TransitionTemperature::Percolation)
  {
    Tstar = GetPercolationTemp();
    Logger::Write(LoggingLevel::TransitionDetailed,
                  "Percolation temperature T = " + std::to_string(Tstar) +
                      " chosen as transition temperature.");
  }
  else if (status_compl == BSMPT::StatusTemperature::Success and
           which_transition_temp == TransitionTemperature::Completion)
  {
    Tstar = GetCompletionTemp();
    Logger::Write(LoggingLevel::TransitionDetailed,
                  "Completion temperature T = " + std::to_string(Tstar) +
                      " chosen as transition temperature.");
  }

  CalculateSoundSpeeds(); // calculate sound speeds in false and true vacuum
}

double BounceSolution::GetPTStrength() const
{
  return alpha;
}

double BounceSolution::TunnelingRate(const double &Temp)
{
  if (Temp < SolutionList.front().T or Temp > SolutionList.back().T)
    return 0; // Never extrapolate
  double Shat3 = GetBounceSol(Temp);
  if (Shat3 < 0)
  {
    Logger::Write(LoggingLevel::BounceDetailed,
                  "Action spline became negative somewhere. T = " +
                      std::to_string(Temp));
    return 1e100; // Spline is unstable or ridiculous action. Return values that
                  // kills the tunneling rate
  }

  double amp =
      std::pow(Temp, 4) * std::pow((Shat3 / (2 * M_PI * Temp)), 3. / 2);
  return amp * std::exp(-Shat3 / Temp);
}

double BounceSolution::HubbleRate(const double &Temp)
{
  const double rhoR = this->GetGstar(Temp) * M_PI * M_PI / 30. *
                      std::pow(Temp, 4); // radiation energy density

  const double DeltaV = FalsePhaseVSpline(Temp) - TruePhaseVSpline(Temp);
  return 1. / (std::sqrt(3) * modelPointer->SMConstants.MPl) *
         std::sqrt(rhoR + DeltaV);
}

double BounceSolution::CalcGstarPureRad()
{
  std::size_t NHiggs = this->modelPointer->get_NHiggs();

  double gb   = 8 * 2 + 4 * 2 + NHiggs;
  double gf   = 6 * 3 * 2 * 2 + 3 * 2 * 2 + 3 * 2;
  double geff = gb + 7. / 8 * gf;

  return geff;
}

void BounceSolution::CalculateNucleationTemp()
{
  if (status_bounce_sol == StatusGW::Success)
  {
    double T_up   = -1;
    double T_down = -1;
    double T_middle;

    for (auto sol = SolutionList.rbegin(); sol != SolutionList.rend(); sol++)
    {
      // Catches the first interval with the nucleation temperature
      if (T_up == -1 and
          TunnelingRate(sol->T) / std::pow(HubbleRate(sol->T), 4) < 1)
        T_up = sol->T;

      if (T_down == -1 and
          TunnelingRate(sol->T) / std::pow(HubbleRate(sol->T), 4) > 1)
        T_down = sol->T;

      if (T_up > 0 and T_down > 0) break;
    }

    if (T_up > 0 and T_down > 0)
    {
      // There is a Tn to be calculated! Use bisection method
      for (int i = 0; i < 50; i++)
      {
        T_middle = (T_up + T_down) / 2;

        if (TunnelingRate(T_middle) / std::pow(HubbleRate(T_middle), 4) < 1)
        {
          T_up = T_middle;
        }
        else
        {
          T_down = T_middle;
        }
        if (std::abs(T_up / T_down - 1) < 1e-10)
        {
          Tnucl               = T_middle;
          nucleation_temp_set = true;
          status_nucl         = BSMPT::StatusTemperature::Success;
          return;
        }
      }
    }
    std::stringstream ss;
    ss << "Nucleation temperature calculation failed\n";

    if (T_up < 0) ss << "Tunneling rate/H never is less than 1.\n";
    if (T_down < 0) ss << "Tunneling rate/H never is more than 1.\n";

    Logger::Write(LoggingLevel::TransitionDetailed, ss.str());

    status_nucl = BSMPT::StatusTemperature::NotMet;
    return;
  }

  return;
}

void BounceSolution::CalculateNucleationTempApprox()
{
  if (status_bounce_sol == StatusGW::Success)
  {
    double T_up   = -1;
    double T_down = -1;
    double T_middle;

    for (auto sol = SolutionList.rbegin(); sol != SolutionList.rend(); sol++)
    {
      // Catches the first interval with the nucleation temperature
      if (T_up == -1 and GetBounceSol(sol->T) / sol->T < 140) T_up = sol->T;

      if (T_down == -1 and GetBounceSol(sol->T) / sol->T > 140) T_down = sol->T;

      if (T_up > 0 and T_down > 0) break;
    }

    if (T_up > 0 and T_down > 0)
    {
      // There is a Tn to be calculated! Use bisection method
      for (int i = 0; i < 50; i++)
      {
        T_middle = (T_up + T_down) / 2;

        if (GetBounceSol(T_middle) / T_middle < 140)
        {
          T_up = T_middle;
        }
        else
        {
          T_down = T_middle;
        }
        if (std::abs(T_up / T_down - 1) < 1e-10)
        {
          Tnucl_approx       = T_middle;
          status_nucl_approx = BSMPT::StatusTemperature::Success;
          return;
        }
      }
    }
    std::stringstream ss;
    ss << "Approximate nucleation temperature calculation failed\n";

    if (T_up < 0) ss << "S3(T)/T never is less than 140.\n";
    if (T_down < 0) ss << "S3(T)/T never is more than 140.\n";

    Logger::Write(LoggingLevel::TransitionDetailed, ss.str());

    status_nucl_approx = BSMPT::StatusTemperature::NotMet;
    return;
  }
  return;
}

double BounceSolution::FalseVacFractionExponent_I(const double &T)
{
  const double prefac = 4. * M_PI / 3. * std::pow(vwall, 3);
  this->SetStoredTemp(T);
  return prefac * Nintegrate_Outer(*this).result;
}

double BounceSolution::CalcFalseVacFraction(const double &temp)
{
  return std::exp(-FalseVacFractionExponent_I(temp));
}

double BounceSolution::CalcTempAtFalseVacFraction(const double &false_vac_frac)
{
  double res_Temp = -1;

  double int_at_false_vac_frac = -std::log(false_vac_frac);
  double T_up                  = -1;
  double T_down                = -1;

  double T_middle;

  for (auto sol = SolutionList.rbegin(); sol != SolutionList.rend(); sol++)
  {
    // catch the first interval containing res_Temp
    double IatT_solT = FalseVacFractionExponent_I(sol->T);

    if (T_up == -1 and IatT_solT < int_at_false_vac_frac) T_up = sol->T;
    if (T_down == -1 and IatT_solT > int_at_false_vac_frac) T_down = sol->T;
    if (T_up > 0 and T_down > 0) break;
  }

  Logger::Write(LoggingLevel::BounceDetailed,
                "T ( Pf = " + std::to_string(std::exp(-int_at_false_vac_frac)) +
                    " ) is in interval [ " + std::to_string(T_down) + ", " +
                    std::to_string(T_up) + " ]");

  if (T_up > 0 and T_down > 0)
  {
    T_middle    = (T_up + T_down) / 2.;
    double IatT = FalseVacFractionExponent_I(T_middle);

    while (not(std::abs(T_up / T_down - 1) <
                   RelativeTemperatureInCalcTempAtFalseVacFraction and
               almost_the_same(int_at_false_vac_frac,
                               IatT,
                               RelativeErrorInCalcTempAtFalseVacFraction)))
    {
      T_middle = (T_up + T_down) / 2.;
      IatT     = FalseVacFractionExponent_I(T_middle);

      Logger::Write(LoggingLevel::BounceDetailed,
                    "Pf ( T = " + std::to_string(T_middle) +
                        " ) = " + std::to_string(std::exp(-IatT)));

      if (IatT < int_at_false_vac_frac)
      {
        T_up = T_middle;
      }
      else
      {
        T_down = T_middle;
      }

      // Condition for success
      if (std::abs(T_up / T_down - 1) <
              RelativeTemperatureInCalcTempAtFalseVacFraction and
          almost_the_same(int_at_false_vac_frac,
                          IatT,
                          RelativeErrorInCalcTempAtFalseVacFraction))
      {
        res_Temp = T_middle;
        break;
      }
    }
  }
  // Not numerically stable
  return res_Temp;
}

void BounceSolution::CalculatePercolationTemp(const double &false_vac_frac)
{
  if (status_bounce_sol == StatusGW::Success)
  {
    Tperc = CalcTempAtFalseVacFraction(false_vac_frac);

    if (Tperc > 0 and percolation_temp_set == false)
    {
      // Try to calculate action at Tp
      // CalculateActionAt(Tperc);
      for (std::size_t i = 0; i < SolutionList.size() - 1; i++)
      {
        if (Tperc > SolutionList[i].T and Tperc < SolutionList[i + 1].T)
          CalculateActionAt((SolutionList[i].T + SolutionList[i + 1].T) / 2.);
      }
      percolation_temp_set = true;
      SetBounceSol();
      Tperc       = CalcTempAtFalseVacFraction(false_vac_frac);
      status_perc = BSMPT::StatusTemperature::Success;
    }
    else if (Tperc < 0)
    {
      Logger::Write(LoggingLevel::TransitionDetailed,
                    "Calculation of the percolation temperature failed.");
      status_perc = BSMPT::StatusTemperature::NotMet;
    }
    return;
  }
  return;
}

void BounceSolution::CalculateCompletionTemp(const double &false_vac_frac)
{
  if (status_bounce_sol == StatusGW::Success)
  {
    Tcompl = CalcTempAtFalseVacFraction(false_vac_frac);

    if (Tcompl > 0 and completion_temp_set == false)
    {
      for (std::size_t i = 0; i < SolutionList.size() - 1; i++)
      {
        if (Tcompl > SolutionList[i].T and Tcompl < SolutionList[i + 1].T)
          CalculateActionAt((SolutionList[i].T + SolutionList[i + 1].T) / 2.);
      }
      completion_temp_set = true;
      Tcompl              = CalcTempAtFalseVacFraction(false_vac_frac);
      status_compl        = BSMPT::StatusTemperature::Success;
    }
    else if (Tcompl < 0)
    {
      Logger::Write(LoggingLevel::TransitionDetailed,
                    "Calculation of the completion temperature failed.");
      status_compl = BSMPT::StatusTemperature::NotMet;
    }
    return;
  }
  return;
}

void BounceSolution::CalculateReheatingTemp()
{
  if (alpha < 1)
  {
    Treh = GetTransitionTemp(); // no supercooling
  }
  else // supercooling
  {
    if (vwall < vCJ)
    {
      std::stringstream ss;
      ss << "\n\033[1;93mWARNING: Detected wall velocity outside detonation "
            "regime with vw = "
         << std::to_string(vwall) << " < vCJ = " << std::to_string(vCJ)
         << " for alpha = " << std::to_string(alpha) << ".\033[0m\n";
      Logger::Write(LoggingLevel::TransitionDetailed, ss.str());
    }
    Treh = GetTransitionTemp() * std::pow((1 + alpha), 1. / 4);
  }
}

double inner_integrand(double Temp, void *params)
{
  class BounceSolution &obj = *static_cast<BounceSolution *>(params);
  double func               = 1. / obj.HubbleRate(Temp);
  return func;
}

double outer_integrand(double Temp, void *params)
{
  class BounceSolution &obj = *static_cast<BounceSolution *>(params);

  double func = obj.TunnelingRate(Temp) /
                (std::pow(Temp, 4) * obj.HubbleRate(Temp)) *
                std::pow(Nintegrate_Inner(obj, Temp).result, 3);
  return func;
}

struct resultErrorPair Nintegrate_Inner(BounceSolution &obj,
                                        const double &Tprime)
{
  double abs_err = obj.AbsErr;
  double rel_err = obj.RelErr;

  std::size_t workspace_size = 1000;
  gsl_integration_workspace *w =
      gsl_integration_workspace_alloc(workspace_size);
  gsl_function F;
  F.function = &inner_integrand;
  F.params   = static_cast<void *>(&obj);

  struct resultErrorPair res;

  gsl_integration_qags(&F,
                       obj.GetStoredTemp(),
                       Tprime,
                       abs_err,
                       rel_err,
                       workspace_size,
                       w,
                       &res.result,
                       &res.error);

  gsl_integration_workspace_free(w);

  return res;
}

struct resultErrorPair Nintegrate_Outer(BounceSolution &obj)
{
  double abs_err = obj.AbsErr;
  double rel_err = obj.RelErr;

  std::size_t workspace_size = 1000;
  gsl_integration_workspace *w =
      gsl_integration_workspace_alloc(workspace_size);
  gsl_function F;
  F.function = &outer_integrand;
  F.params   = static_cast<void *>(&obj);

  struct resultErrorPair res;

  gsl_integration_qags(&F,
                       obj.GetStoredTemp(),
                       obj.GetCriticalTemp(),
                       abs_err,
                       rel_err,
                       workspace_size,
                       w,
                       &res.result,
                       &res.error);

  gsl_integration_workspace_free(w);

  return res;
}

double BounceSolution::CalculateRhoGamma(const double &T) const
{
  return this->GetGstar(T) * std::pow(M_PI, 2) / 30 * std::pow(T, 4);
}

void BounceSolution::CalculatePTStrength()
{
  if (UserDefined_vwall >= 0)
    vwall = UserDefined_vwall;
  else
    vwall = 0.95; // Initial guess

  double old_alpha; // To keep track of the convergence

  for (int c = 0; c < 20; c++)
  {
    // Use recursive method to find solution of
    // alpha = alpha(T_*)
    // T_* =  T_*(alpha, v_wall)
    // v_wall = v_wall(alpha, T_*)
    // Should converge quickly, if fails use default value of v_wall = 0.95
    old_alpha = alpha;
    CalcTransitionTemp(); // Calculation all Ts
    Minimum true_min  = phase_pair.true_phase.Get(GetTransitionTemp());
    Minimum false_min = phase_pair.false_phase.Get(GetTransitionTemp());

    double Vi =
        false_min
            .potential; // potential at false vacuum and percolation temperature
    double Vf =
        true_min
            .potential; // potential at true vacuum and percolation temperature
    double dTVi = this->modelPointer->VEff(
        this->modelPointer->MinimizeOrderVEV(false_min.point),
        GetTransitionTemp(),
        -1); // temperature-derivative at false vacuum
    double dTVf = this->modelPointer->VEff(
        this->modelPointer->MinimizeOrderVEV(true_min.point),
        GetTransitionTemp(),
        -1); // temperature-derivative at true vacuum const

    double rho_gam = CalculateRhoGamma(GetTransitionTemp());
    alpha = 1 / rho_gam * (Vi - Vf - GetTransitionTemp() / 4. * (dTVi - dTVf));
    CalculateWallVelocity(false_min, true_min);
    if (abs(alpha / old_alpha - 1) < 1e-7) return; // Found a solution
  }
  // We could not find the solution for the system. use default value of .95
  // instead
  Logger::Write(LoggingLevel::TransitionDetailed,
                "v_wall could not be calculated. Using 0.95 for v_wall.");
  UserDefined_vwall = 0.95;
  // Call this function recursively
  CalculatePTStrength();
  return;
}

void BounceSolution::CalcChapmanJougetVelocity()
{
  vCJ = (1 + std::sqrt(3 * alpha *
                       (1 - Csound_false * Csound_false +
                        3 * Csound_false * Csound_false * alpha))) /
        (1. / Csound_false + 3 * Csound_false * alpha);
}

void BounceSolution::CalculateWallVelocity(const Minimum &false_min,
                                           const Minimum &true_min)
{
  const double Temp = false_min.temp; // temperature

  // User defined
  if (UserDefined_vwall >= 0) vwall = UserDefined_vwall;
  if (UserDefined_vwall == -1)
  {
    // vwall https://arxiv.org/abs/2210.16305
    double Vi = false_min.potential; // potential at false vacuum
    double Vf = true_min.potential;  // potential at true vacuum

    double rho_gam = CalculateRhoGamma(Temp);
    // Candidate wall velocity
    vwall = std::sqrt((Vi - Vf) / (alpha * rho_gam));
    // If candidate is bigger than chapman jouget velocity, v = 1
    if (vwall > vCJ) vwall = 1;
  }
  if (UserDefined_vwall == -2)
  {
    double dTVi = this->modelPointer->VEff(
        this->modelPointer->MinimizeOrderVEV(false_min.point),
        Temp,
        -1); // temperature-derivative at false vacuum
    double dTVf = this->modelPointer->VEff(
        this->modelPointer->MinimizeOrderVEV(true_min.point),
        Temp,
        -1); // temperature-derivative at true vacuum const

    double psi = dTVf / dTVi;
    double a   = 0.2233;
    double b   = 1.704;
    double p   = -3.433;

    vwall = std::pow(
        pow(abs((3 * alpha + psi - 1) / (2 * (2 - 3 * psi + std::pow(psi, 3)))),
            p / 2) +
            std::pow(abs(vCJ * (1 - a * std::pow(1 - psi, b) / alpha)), p / 2),
        1 / p);
  }
}

double BounceSolution::CalculateSoundSpeed(Phase &phase)
{
  const double eps                         = 0.01;
  std::function<double(Minimum minimum)> V = [&](Minimum minimum)
  {
    // Potential wrapper
    std::vector<double> res = modelPointer->MinimizeOrderVEV(minimum.point);
    return modelPointer->VEff(res, minimum.temp);
  };
  const double V_before = V(phase.Get(Tstar + eps));
  const double V_tstar  = V(phase.Get(Tstar));
  const double V_after  = V(phase.Get(Tstar - eps));
  const double dVdT     = (V_before - V_after) / (2. * eps);
  const double d2VdT2   = (V_before - 2. * V_tstar + V_after) / (eps * eps);
  const double cs       = sqrt(dVdT / (d2VdT2 * Tstar));
  if (isnan(cs))
  {
    stringstream ss;
    ss << "Sound speed calculation failed!" << "\n";
    ss << "V(T-eps) V(T) V(T + eps) " << V_after << " " << V_tstar << " "
       << V_before << "\n";
    ss << "dVdT = \t" << dVdT << "\n";
    ss << "d2VdT2 = \t" << d2VdT2 << "\n";
    ss << "Using cs = 1/sqrt(3) instead.";
    Logger::Write(LoggingLevel::GWDetailed, ss.str());
    return 1. / sqrt(3.);
  }
  return cs;
}

void BounceSolution::CalculateSoundSpeeds()
{
  Csound_false = CalculateSoundSpeed(phase_pair.false_phase);
  Csound_true  = CalculateSoundSpeed(phase_pair.true_phase);

  CalcChapmanJougetVelocity();
}

double BounceSolution::GetBounceSol(const double &Temp) const
{
  return S3ofT_spline(Temp);
}

double action_ratio(double Temp, void *params)
{
  class BounceSolution &obj = *static_cast<BounceSolution *>(params);
  double func               = obj.GetBounceSol(Temp) / Temp;
  return func;
}

struct resultErrorPair Nderive_BounceRatio(BounceSolution &obj)
{
  double step_size = 1e-8;

  gsl_function F;
  F.function = &action_ratio;
  F.params   = static_cast<void *>(&obj);

  struct resultErrorPair res;

  gsl_deriv_central(
      &F, obj.GetTransitionTemp(), step_size, &res.result, &res.error);

  return res;
}

void BounceSolution::CalculateInvTimeScale()
{
  struct resultErrorPair res = Nderive_BounceRatio(*this);
  this->betaH                = this->GetTransitionTemp() * res.result;
}

double BounceSolution::GetInvTimeScale()
{
  if (this->betaH == -1) CalculateInvTimeScale();
  return this->betaH;
}

double BounceSolution::GetRstar()
{
  if (this->Rstar == -1) CalculateRstar();
  return this->Rstar;
}

void BounceSolution::CalculateRstar()
{
  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);

  gsl_function F;
  F.function = [](double T, void *params) -> double
  {
    return static_cast<BounceSolution *>(params)->TunnelingRate(T) *
           exp(-static_cast<BounceSolution *>(params)
                    ->FalseVacFractionExponent_I(T)) /
           (static_cast<BounceSolution *>(params)->HubbleRate(T) * pow(T, 4));
  };
  F.params = this; // Pass `this` pointer as parameters

  double result, error;
  gsl_integration_qags(
      &F, Tstar, Tc, RelErr, AbsErr, 1000, workspace, &result, &error);

  gsl_integration_workspace_free(workspace);

  this->Rstar = pow(pow(Tstar, 3) * result, -1 / 3.);
}
} // namespace BSMPT
