// Copyright (C) 2024 Lisa Biermann, Margarete Mühlleitner, Rui Santos, João
// Viana
//
// SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

/**
 * @file minimum tracer class
 */
#include "Eigen/Eigenvalues"                   // Eigenvalues utility
#include <BSMPT/minimizer/Minimizer.h>         // for Minimizer
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/utility/Logger.h>              // for Logger Class
#include <BSMPT/utility/asciiplotter/asciiplotter.h>
#include <BSMPT/utility/utility.h>
#include <Eigen/Dense> // Eigenvalues matrix
#include <chrono>
#include <cmath>    // std::pow
#include <memory>   // for shared_ptr
#include <optional> // std::optional
#include <stdlib.h> // std::strtol

namespace BSMPT
{

/**
 * @brief Value to be store in the columns without any values
 *
 */
const double EmptyValue = NAN;
/**
 * @brief Possible NLO stability status
 *
 */
enum class StatusNLOStability
{
  NotSet,
  Off,
  Success,
  NoNLOStability
};
/**
 * @brief Map to convert StatusNLOStability to strins
 *
 */
const std::unordered_map<StatusNLOStability, std::string>
    StatusNLOStabilityToString{
        {StatusNLOStability::NotSet, "not_set"},
        {StatusNLOStability::Off, "off"},
        {StatusNLOStability::Success, "success"},
        {StatusNLOStability::NoNLOStability, "no_nlo_stability"}};
/**
 * @brief Possible electroweak symmetry restoration status
 *
 */
enum class StatusEWSR
{
  NotSet,
  Off,
  Failure,
  NotBFB,
  FlatRegion,
  EWSymNonRes,
  EWSymRes
};
/**
 * @brief Map to convert StatusEWSRToString to strins
 *
 */
const std::unordered_map<StatusEWSR, std::string> StatusEWSRToString{
    {StatusEWSR::NotSet, "not_set"},
    {StatusEWSR::Off, "off"},
    {StatusEWSR::Failure, "failure"},
    {StatusEWSR::NotBFB, "non_bfb"},
    {StatusEWSR::FlatRegion, "flat_region"},
    {StatusEWSR::EWSymNonRes, "ew_sym_non_res"},
    {StatusEWSR::EWSymRes, "ew_sym_res"}};
/**
 * @brief Possible tracing results
 *
 */
enum class StatusTracing
{
  NotSet,
  Success,
  NoCoverage,
  NoMinsAtBoundaries,
  NoGlobMinCoverage,
  Failure
};
/**
 * @brief Map to convert StatusTracingToString to strins
 *
 */
const std::unordered_map<StatusTracing, std::string> StatusTracingToString{
    {StatusTracing::NotSet, "not_set"},
    {StatusTracing::Success, "success"},
    {StatusTracing::NoCoverage, "no_coverage"},
    {StatusTracing::NoMinsAtBoundaries, "no_mins_at_boundaries"},
    {StatusTracing::NoGlobMinCoverage, "no_glob_min_coverage"},
    {StatusTracing::Failure, "failure"}};
/**
 * @brief Possible status for the coex phase
 *
 */
enum class StatusCoexPair
{
  NotSet,
  Success,
  NoCoexPairs
};
/**
 * @brief Map to convert StatusCoexPairToString to strings
 *
 */
const std::unordered_map<StatusCoexPair, std::string> StatusCoexPairToString{
    {StatusCoexPair::NotSet, "not_set"},
    {StatusCoexPair::Success, "success"},
    {StatusCoexPair::NoCoexPairs, "no_coex_pair"}};
/**
 * @brief Possible status for the critical temperature
 *
 */
enum class StatusCrit
{
  NotSet,
  Success,
  FalseLower,
  TrueLower,
  Failure
};
/**
 * @brief Map to convert StatusCritToString to strings
 *
 */
const std::unordered_map<StatusCrit, std::string> StatusCritToString{
    {StatusCrit::NotSet, "not_set"},
    {StatusCrit::Success, "success"},
    {StatusCrit::FalseLower, "false_lower"},
    {StatusCrit::TrueLower, "true_lower"},
    {StatusCrit::Failure, "failure"}};
/**
 * @brief Possible status for the approximated nucleation, exact nucleation,
 * percolation and completion temperature
 *
 */
enum class StatusTemperature
{
  NotSet,
  Success,
  NotMet,
  NaN
};
/**
 * @brief Map to convert StatusTemperature to strings
 *
 */
const std::unordered_map<StatusTemperature, std::string>
    StatusTemperatureToString{{StatusTemperature::NotSet, "not_set"},
                              {StatusTemperature::Success, "success"},
                              {StatusTemperature::NotMet, "not_met"},
                              {StatusTemperature::NaN, "nan"}};
/**
 * @brief Possible results for the GW and bounce_sol class.
 *
 */
enum class StatusGW
{
  NotSet,
  Success,
  Failure
};
/**
 * @brief Map to convert StatusGWToString to strings
 *
 */
const std::unordered_map<StatusGW, std::string> StatusGWToString{
    {StatusGW::NotSet, "not_set"},
    {StatusGW::Success, "success"},
    {StatusGW::Failure, "failure"}};

/**
 * @brief Possible transitions temperatures
 *
 */
enum class TransitionTemperature
{
  NotSet,
  ApproxNucleation,
  Nucleation,
  Percolation,
  Completion
};

/**
 * @brief Override << operator to handle StatusNLOStability
 *
 * @param os ostream buffer
 * @param status status to be printed
 * @return std::ostream& buffer
 */
std::ostream &operator<<(std::ostream &os, const StatusNLOStability &status);
/**
 * @brief Override << operator to handle StatusEWSR
 *
 * @param os ostream buffer
 * @param status status to be printed
 * @return std::ostream& buffer
 */
std::ostream &operator<<(std::ostream &os, const StatusEWSR &status);
/**
 * @brief Override << operator to handle StatusTracing
 *
 * @param os ostream buffer
 * @param status status to be printed
 * @return std::ostream& buffer
 */
std::ostream &operator<<(std::ostream &os, const StatusTracing &status);
/**
 * @brief Override << operator to handle StatusCoexPair
 *
 * @param os ostream buffer
 * @param status status to be printed
 * @return std::ostream& buffer
 */
std::ostream &operator<<(std::ostream &os, const StatusCoexPair &status);
/**
 * @brief Override << operator to handle StatusCrit
 *
 * @param os ostream buffer
 * @param status status to be printed
 * @return std::ostream& buffer
 */
std::ostream &operator<<(std::ostream &os, const StatusCrit &status);
/**
 * @brief Override << operator to handle StatusGW
 *
 * @param os ostream buffer
 * @param status status to be printed
 * @return std::ostream& buffer
 */
std::ostream &operator<<(std::ostream &os, const StatusGW &status);
/**
 * @brief Override << operator to handle StatusTemperature
 *
 * @param os ostream buffer
 * @param status status to be printed
 * @return std::ostream& buffer
 */
std::ostream &operator<<(std::ostream &os, const StatusTemperature &status);

/**
 * @brief struct to store minimum and temperature
 * @param point coordinates in field space
 * @param temp temperature
 * @param potential value of the potential
 * @param is_glob_min true if minimum is global minimum
 * @param EdgeOfPhase 1 = starting minimum | 0 = Middle minimum | -1 = Ending
 * minimum (sum of EdgeOfPhase is the number of coexisting phases)
 */
struct Minimum
{
  std::vector<double> point;
  double temp;
  double potential;
  bool is_glob_min = false;
  int EdgeOfPhase  = 0;

  bool operator<(const Minimum &a) const { return temp < a.temp; }
};

class MinimumTracer
{
private:
  int WhichMinimizer;
  bool UseMultithreading;

protected:
  /**
   * @brief modelPointer for the used parameter point
   */
  std::shared_ptr<Class_Potential_Origin> modelPointer;

public:
  /**
   * @brief Threshold for the acceptable gradient
   *
   */
  double GradientThreshold = 1e-3;

  /**
   * @brief Add a constant to the diagonals of the hessian matrix in the
   * LocateMinimum function. Helps with convergence.
   *
   */
  double HessianDiagonalShift = 1e-3;

  /**
   * @brief Minimum found in IsThereEWSymmetryRestoration()
   *
   */
  std::vector<double> HighTemperatureVEV;

  /**
   * @brief Vector to store minima that appeared from VEV splittings
   *
   */
  std::vector<Minimum> SavedMinimaFromVEVSplitting;

  /**
   * @brief default constructor
   */
  MinimumTracer();

  /**
   * @brief constructor
   * @param pointer_in this->modelPointer for used parameter point
   * @param WhichMinimizer_in which minimizers are used
   * @param UseMultithreading_in whether or not multithreading is used
   */
  MinimumTracer(const std::shared_ptr<Class_Potential_Origin> &pointer_in,
                const int &WhichMinimizer_in,
                const bool &UseMultithreading_in);

  /**
   * @brief Calculates flat field directions
   */
  void FindFlatDirections();

  /**
   * @brief Convert point into minimal non-flat space, reduces dimension in case
   * of flat directions, point has to have VEV dimension
   */
  void ConvertToNonFlatDirections(std::vector<double> &point);

  /**
   * @brief Calculates the list of symmetries that leave V unchanged
   *
   */
  void FindDiscreteSymmetries();

  /**
   * @brief ConvertToVEVDim converts point from full to reduced (VEV) dimension
   * @param point point in full field dimension
   * @return point in reduced VEV dimension
   */
  std::vector<double> ConvertToVEVDim(const std::vector<double> &point);

  /**
   * @brief get global minimum of effective potential
   * @param Temp temperature
   * @param check storage for minimization debugging options
   * @param start start value for CMA-ES minimization
   * @return global minimum at temperature Temp
   */
  std::vector<double> GetGlobalMinimum(const double &Temp,
                                       std::vector<double> &check,
                                       const std::vector<double> &start);

  /**
   * @brief get global minimum of effective potential
   * @param Temp temperature
   * @param start start value for CMA-ES minimization
   * @return global minimum at temperature Temp
   */
  std::vector<double> GetGlobalMinimum(const double &Temp,
                                       const std::vector<double> &start);

  /**
   * @brief get global minimum of effective potential
   * @param Temp temperature
   * @return global minimum at temperature Temp
   */
  std::vector<double> GetGlobalMinimum(const double &Temp);

  /**
   * @brief IsGlobMin checks whether current minimum is the global minimum
   * @param min Minimum to check, sets is_glob_min to true if min is global
   * minimum
   */
  void IsGlobMin(Minimum &min);

  /**
   * @brief GetStatusNLOVEV convert bool output of CheckNLOVEV to
   * status string
   * @param out bool output of CheckNLOVEV
   * @return status string for output
   */
  StatusNLOStability GetStatusNLOVEV(const bool &out);

  /**
   * @brief GetStatusEWSR convert double output of IsThereEWSymmetryRestoration
   * to status string
   * @param out int output of IsThereEWSymmetryRestoration
   * @return status string for output
   */
  StatusEWSR GetStatusEWSR(const int &out);

  /**
   * @brief IsThereEWSymmetryRestoration checks if there is EW symmetry
   * restoration at high temperatures
   * @return status; if status < -1: no minimum if status = 0: calculation
   * failed, if status = 1: flat region with minimum, if status = 2: there is a
   * minimum but no symmetry restoration, if status = 3: symmetry restoration
   */
  int IsThereEWSymmetryRestoration();

  /**
   * @brief SmallestEigenvalue calculate Eigenvalues of Hessian and returns
   * smallest
   * @param point point where to evaluate the Hessian
   * @param Hessian Hessian function
   * @return smallest Eigenvalue of Hessian
   */
  double SmallestEigenvalue(
      const std::vector<double> &point,
      const std::function<std::vector<std::vector<double>>(std::vector<double>)>
          &Hessian);

  /**
   * @brief FindZeroSmallestEigenvalue
   * @param point_1 first point
   * @param T_1 temperature of first point
   * @param point_2 second point
   * @param T_2 temperature of second point
   * @return stationary point found in between point_1 and point_2
   */
  std::vector<double> FindZeroSmallestEigenvalue(std::vector<double> point_1,
                                                 double T_1,
                                                 std::vector<double> point_2,
                                                 double T_2);

  /**
   * @brief TrackPhase with enforced global minimum tracing (= phase is checked
   * if it is still the global minimum until it is no longer, then the current
   * temperature is stored in globMinEndT)
   * @param globMinEndT temperature at which phase is no longer global minimum
   * @param point_In start point for tracking
   * @param currentT_In start point temperature for phase tracking
   * @param finalT end point temperature
   * @param dT_In initial temperature step size
   * @param output if true tracking output is printed on the screen
   * @param unprotected if true we dont check the hessian
   */
  std::vector<Minimum> TrackPhase(double &globMinEndT,
                                  const std::vector<double> &point_In,
                                  const double &currentT_In,
                                  const double &finalT,
                                  const double &dT_In     = 1,
                                  const bool &output      = true,
                                  const bool &unprotected = false);

  /**
   * @brief TrackPhase
   * @param point_In start point for tracking
   * @param currentT_In start point temperature for phase tracking
   * @param finalT end point temperature
   * @param dT_In initial temperature step size
   * @param output if true tracking output is printed on the screen
   * @param unprotected if true we dont check the hessian
   */
  std::vector<Minimum> TrackPhase(const std::vector<double> &point_In,
                                  const double &currentT_In,
                                  const double &finalT,
                                  const double &dT_In     = 1,
                                  const bool &output      = true,
                                  const bool &unprotected = false);

  /**
   * @brief bool to store whether flat directions are found
   */
  bool flat_dirs_found = false;

  /**
   * @brief storage of all non-flat VEV-directions
   */
  std::vector<int> NonFlatDirections;

  /**
   * @brief storage of indices of flat 1D directions in VEV basis
   */
  std::vector<std::size_t> flat_1D_dirs;

  /**
   * @brief storage of indices of flat 2D directions in VEV basis
   */
  std::vector<std::vector<std::size_t>> flat_2D_dirs;

  /**
   * @brief storage of indices of flat 3D directions in VEV basis
   */
  std::vector<std::vector<std::size_t>> flat_3D_dirs;

  /**
   * @brief List of group elements allowed by the potential
   *
   */
  std::vector<Eigen::MatrixXd> GroupElements;

  /**
   * @brief Reduce the VEV into the same principal quadrant.
   *
   * @param point to be rotated
   */
  void ReduceVEV(std::vector<double> &vev);

  /**
   * @brief Reduce the VEV of the minimum into the same principal quadrant.
   *
   * @param min to be rotated
   */
  void ReduceVEV(Minimum &min);

  /**
   * @brief WarpPath
   * @param path
   * @param T1
   * @param F1
   * @param T2
   * @param F2
   * @return wraped path
   */
  const std::vector<std::vector<double>>
  WarpPath(const std::vector<std::vector<double>> &path,
           const std::vector<double> &T1,
           const std::vector<double> &F1,
           const std::vector<double> &T2,
           const std::vector<double> &F2);

  /**
   * @brief Finds stationary points of a function (not only minimas).
   *
   * @param guess_In is the initial guess for the minimum
   * @param df gradient of the function to be minimized
   * @param Hessian hessian of the function
   * @param error Maximum size of \f$ | \vec{df} | \f$ that is considered a
   * minimum
   * @param const_multiplier If \f$ \det{Hessian} = 0\f$ this method does not
   * work. In that case we move the guess as \f$ \vec{p} \rightarrow \vec{p} -
   * const\_multiplier * \vec{df}\f$
   * @param maxiter Maximum iteration exiting function
   * @return std::vector<double>
   */
  std::vector<double> LocateMinimum(
      const std::vector<double> &guess_In,
      std::function<std::vector<double>(std::vector<double>)> &df,
      std::function<std::vector<std::vector<double>>(std::vector<double>)>
          &Hessian,
      const double &error            = 1e-4,
      const double &const_multiplier = 1e-2,
      const int &maxiter             = 100);

  /**
   * @brief GetLegend derive legend
   * @param num_coex_phases number of coexisting phase regions
   * @param do_gw_calc bool that determines whether gw calculation is performed
   * @return vector of column label strings
   */
  std::vector<std::string> GetLegend(const int &num_coex_phases,
                                     const bool &do_gw_calc);
};

/**
 * @brief Create1DimGrid creates a 1-dim grid of given size in index-direction
 * @param point
 * @param k index direction in which to create grid
 * @param low_value
 * @param high_value
 * @param nsteps
 * @return gridsize-dim vector of vectors
 */
std::vector<std::vector<double>>
Create1DimGrid(const std::vector<double> &point,
               const int k,
               const double low_value,
               const double high_value,
               const int nsteps = 100);

/**
 * @brief Create1DimGrid creates a 1-dim grid of given size between two points
 * @param min_start
 * @param min_end
 * @param npoints
 * @return npoints long vector of steps on connecting line between min_start and
 * min_end
 */
std::vector<std::vector<double>>
Create1DimGrid(const std::vector<double> &min_start,
               const std::vector<double> &min_end,
               const int npoints = 100);

/**
 * Returns true if two values are the same given some relative precision
 */
bool almost_the_same(const double &a,
                     const double &b,
                     const double &rel_precision = 0.01,
                     const double &num_zero      = 1e-10);

/**
 * Returns true if two vectors are the element-wise the same given some relative
 * precision
 */
bool almost_the_same(const std::vector<double> &a,
                     const std::vector<double> &b,
                     const bool &allow_for_sign_flip = false,
                     const double &rel_precision     = 0.01,
                     const double &num_zero          = 1e-10);

/**
 * @brief Phase object
 *
 */
struct Phase
{
  /**
   * @brief phase ID
   */
  int id = 0;

  /**
   * @brief Lowest temperature of the phase
   */
  double T_low = 0;

  /**
   * @brief Highest temperature of the phase
   */
  double T_high = 0;

  /**
   * @brief Set of Minimum that compose the phase
   */
  std::vector<Minimum> MinimumPhaseVector;

  /**
   * @brief MinTracer object
   */
  std::shared_ptr<MinimumTracer> MinTracer;

  /**
   * @brief Function that adds min to MinimumPhaseVector
   *
   * @param min Minimum
   */
  void Add(Minimum min);

  /**
   * @brief Calculates the minimum of the phase at temperature T. This function
   * assumes that T is inside the temperature range of the phase.
   *
   * @param T temperature
   * @return Minimum at temperature T for the corresponding phase
   */
  Minimum Get(double T);

  /**
   * @brief empty constructor
   */
  Phase();

  /**
   * @brief Construct a new Phase object with enforced global minimum tracing
   * @param phase_start Initial starting point for the phase
   * @param initialT Initial temperature
   * @param finalT Final temperature
   * @param globMinEndT Temperature where phase is first detected to no longer
   * be global minimum
   * @param MinTracerIn MinTracer pointer
   */
  Phase(const std::vector<double> &phase_start,
        const double &initialT,
        const double &finalT,
        double &globMinEndT,
        std::shared_ptr<MinimumTracer> &MinTracerIn);

  /**
   * @brief Construct a new Phase:: Phase object
   *
   * @param phase_start Initial starting point for the phase
   * @param initialT Temperature of the phase given as input
   * @param LowT Lowest temperature
   * @param HighT Highest temperature
   * @param globMinEndT Temperature where phase is first detected to no longer
   * be global minimum
   * @param MinTracerIn MinTracer pointer
   */
  Phase(const double &initialT,
        const std::vector<double> &phase_start,
        const double &LowT,
        const double &HighT,
        double &globMinEndT,
        std::shared_ptr<MinimumTracer> &MinTracerIn);

  /**
   * @brief Construct a new Phase object
   *
   * @param phase_start Initial starting point for the phase
   * @param initialT Initial temperature
   * @param finalT Final temperature
   * @param globMinEndT Temperature where phase is first detected to no longer
   * be global minimum
   * @param MinTracerIn MinTracer pointer
   */
  Phase(const std::vector<double> &phase_start,
        const double &initialT,
        const double &finalT,
        std::shared_ptr<MinimumTracer> &MinTracerIn);

  /**
   * @brief Construct a new Phase:: Phase object
   *
   * @param phase_start Initial starting point for the phase
   * @param initialT Temperature of the phase given as input
   * @param LowT Lowest temperature
   * @param HighT Highest temperature
   * @param MinTracerIn MinTracer pointer
   */
  Phase(const double &initialT,
        const std::vector<double> &phase_start,
        const double &LowT,
        const double &HighT,
        std::shared_ptr<MinimumTracer> &MinTracerIn);

  /**
   * @brief Construct a new Phase:: Phase object
   *
   * @param initialT Temperature of the phase given as input
   * @param LowT Lowest temperature
   * @param HighT Highest temperature
   * @param MinTracerIn MinTracer pointer
   */
  Phase(const double &initialT,
        const double &LowT,
        const double &HighT,
        std::shared_ptr<MinimumTracer> &MinTracerIn);
};

/**
 * @brief CoexPhases struct to save pair of coexisting phases (false and true
 * phase)
 */
struct CoexPhases
{
  int coex_pair_id;
  double T_high;
  double T_low;
  Phase false_phase;
  Phase true_phase;
  double crit_temp       = -1;
  StatusCrit crit_status = StatusCrit::NotSet;

  /**
   * @brief empty constructor
   */
  CoexPhases();

  /**
   * @brief constructor
   */
  CoexPhases(const int &pair_id,
             const Phase &false_phase,
             const Phase &true_phase,
             const double &Tlow_in,
             const double &Thigh_in);

  /**
   * @brief CalculateTc critical temperature for coexising phase pair
   */
  void CalculateTc();
};

/**
 * @brief Complete vacuum structure of the theory for this parameter point
 *
 */
struct Vacuum
{
  /**
   * @brief Lowest temperature 0 GeV
   *
   */
  double T_low = 0;

  /**
   * @brief Highest temperature, 300 GeV or set in input file
   *
   */
  double T_high = 300;

  /**
   * @brief Lowest temperature at which high-temperature phase is found to
   * exist
   *
   */
  double T_low_highTempPhase = -1;

  /**
   * @brief Highest temperature at which low-temperature phase is found to
   * exist
   *
   */
  double T_high_lowTempPhase = -1;

  /**
   * @brief number of equally-spaced intermediate points to check for new
   * phases
   */
  int num_points = 0;

  /**
   * @brief vacuum status code = success, no_coverage, no_glob_min_coverage
   */
  StatusTracing status_vacuum = StatusTracing::NotSet;

  /**
   * @brief coexisting phases status code = success, no_coex_pairs
   */
  StatusCoexPair status_coex_pairs = StatusCoexPair::NotSet;

  /**
   * @brief MinTracer object
   *
   */
  std::shared_ptr<MinimumTracer> MinTracer;

  /**
   * @brief Model pointer
   *
   */
  std::shared_ptr<Class_Potential_Origin> modelPointer;

  /**
   * @brief List of different phases
   *
   */
  std::vector<Phase> PhasesList;

  /**
   * @brief List of different phase pairs
   *
   */
  std::vector<CoexPhases> CoexPhasesList;

  /**
   * @brief Construct a new Vacuum object
   *
   * @param T_lowIn Lowest temperature, 0 GeV or set in input file
   * @param T_highIn Highest temperature, 300 GeV or set in input file
   * @param MinTracerIn  MinTracer object
   * @param modelPointerIn Model pointer
   * @param UseMultiStepPTModeIn choose multi-step PT modes: default (= -1),
   * 0, 1, 2, auto (= 3)
   * @param num_pointsIn number of equally-spaced intermediate points to check
   * for new phases
   * @param do_only_tracing if true only tracing and no identification of all
   * possible coexisting phase pairs and their critical temperatures is done, if
   * false identification and calculation of Tc is done
   */
  Vacuum(const double &T_lowIn,
         const double &T_highIn,
         std::shared_ptr<MinimumTracer> &MinTracerIn,
         std::shared_ptr<Class_Potential_Origin> &modelPointerIn,
         const int &UseMultiStepPTModeIn,
         const int &num_pointsIn     = 10,
         const bool &do_only_tracing = false);

  /**
   * @brief MultiStepPTTracer traces all phases between T_high and T_low
   * assuming that we start and end in a global minimum (absolute stability)
   * @param Temp temperature at which to evaluate phase
   * @param deltaT temperature step to take if already traced phase is found
   * again
   */
  void MultiStepPTTracer(const double &Temp, const double &deltaT = 0);

  /**
   * @brief print info on phase
   * @param phase Phase object
   */
  void print(const Phase &phase);

  /**
   * @brief setCoexPhases Calculates all coexisting phase pairs irrespective of
   * borders of coexisiting phase regions, this way we take into account the
   * full region of coexistance for coexisting phases and do not split at the
   * border of coexisting regions which are split as soon as the number of
   * coexisting phases changes which excludes transitions that have a critical
   * temperature in the first region but only complete in the second region (in
   * case the critical temperature is close to the region border or the
   * temperature difference is large enough)
   */
  void setCoexPhases();

  /**
   * @brief setCoexRegion Calculates all coexisting phase regions with phase
   * pairs included from the phase vector
   * @param UseMultiStepPTMode int to distinguish multistep PT mode, for all
   * modes except mode 0 we try to patch up holes in tracing
   */
  void setCoexRegion(const int &UseMultiStepPTMode);

  /**
   * @brief Adds a phase to the phase list
   *
   * @param phase
   */
  void addPhase(Phase &phase);

  /**
   * @brief This function checks if the minimum already exists in one of the
   * phases in the phase list
   * @param minimum to be checked
   * @return int = -1 if minimum is not in phase. Otherwise returns index of
   * phase.
   */
  int MinimumFoundAlready(const Minimum &minimum);

  /**
   * @brief MultiStepPTMode0 single-step PT mode
   * @param LowTempPoint tracing starting point at T_low
   * @param HighTempPoint tracing starting point at T_high
   */
  void MultiStepPTMode0(const std::vector<double> &LowTempPoint,
                        const std::vector<double> &HighTempPoint);

  /**
   * @brief MultiStepPTMode1 multi-step coverage PT mode
   * @param LowTempPoint tracing starting point at T_low
   * @param HighTempPoint tracing starting point at T_high
   */
  void MultiStepPTMode1(const std::vector<double> &LowTempPoint,
                        const std::vector<double> &HighTempPoint);

  /**
   * @brief DoPhasesOverlap checks if two phases overlap
   */
  bool DoPhasesOverlap(Phase &new_phase, Phase &old_phase);

  /**
   * @brief DoGlobMinOverlap check for global minimum at left-/rightmost
   * overlap and choose endpoint of previous phase if new phase does not overlap
   */
  bool DoGlobMinOverlap(const Phase &new_phase, const Phase &old_phase);

  /**
   * @brief MultiStepPTMode2 multi-step global minimum coverage PT mode
   * @param LowTempPoint tracing starting point at T_low
   * @param HighTempPoint tracing starting point at T_high
   */
  void MultiStepPTMode2(const std::vector<double> &LowTempPoint,
                        const std::vector<double> &HighTempPoint);

  /**
   * @brief prints the phases from T_low up to T_high on the terminal
   *
   * @param size of diagram in terminal
   */
  void PrintPhasesDiagram(int size = 100);
};

} // namespace BSMPT