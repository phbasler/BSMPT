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
#include <BSMPT/utility/asciiplotter.h>
#include <BSMPT/utility/utility.h>
#include <Eigen/Dense> // Eigenvalues matrix
#include <chrono>
#include <cmath>    // std::pow
#include <memory>   // for shared_ptr
#include <stdlib.h> // std::strtol

namespace BSMPT
{

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
  MinimumTracer(std::shared_ptr<Class_Potential_Origin> &pointer_in,
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
  std::string GetStatusNLOVEV(const bool &out);

  /**
   * @brief GetStatusEWSR convert double output of IsThereEWSymmetryRestoration
   * to status string
   * @param out int output of IsThereEWSymmetryRestoration
   * @return status string for output
   */
  std::string GetStatusEWSR(const int &out);

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
   * @brief This function check if the gradient increses from point_1 to
   * point_2. Not used at the moment.
   *
   * @param T is the temperature.
   * @param point_1 initial point.
   * @param point_2 final point.
   * @param V is the potential.
   * @param dV is the potential gradient.
   * @return true if there is a barrier between the minima.
   * @return false if there is no barrier between the minima.
   */
  bool GradientIncreases(
      const double &T,
      const std::vector<double> &point_1,
      const std::vector<double> &point_2,
      const std::function<double(std::vector<double>)> &V,
      const std::function<std::vector<double>(std::vector<double>)> &dV);

  /**
   * @brief trackPhase with enforced global minimum tracing (= phase is checked
   * if it is still the global minimum until it is no longer, then the current
   * temperature is stored in globMinEndT)
   * @param globMinEndT temperature at which phase is no longer global minimum
   * @param point start point for tracking
   * @param currentT start point temperature for phase tracking
   * @param finalT end point temperature
   * @param dT initial temperature step size
   * @param output if true tracking output is printed on the screen
   * @param unprotected if true we dont check the hessian
   */
  std::vector<Minimum> trackPhase(double &globMinEndT,
                                  std::vector<double> point,
                                  double currentT,
                                  double finalT,
                                  double dT               = 3,
                                  const bool &output      = true,
                                  const bool &unprotected = false);

  /**
   * @brief trackPhase
   * @param point start point for tracking
   * @param currentT start point temperature for phase tracking
   * @param finalT end point temperature
   * @param dT initial temperature step size
   * @param output if true tracking output is printed on the screen
   * @param unprotected if true we dont check the hessian
   */
  std::vector<Minimum> trackPhase(std::vector<double> point,
                                  double currentT,
                                  double finalT,
                                  double dT               = 3,
                                  const bool &output      = true,
                                  const bool &unprotected = false);

  /**
   * @brief Calculates the VEV splitting when Hessian matrix gets a single
   * negative eigenvalue.
   *
   * @param point where grad is zero and hessian is not positive definite
   * @param T is the temperature.
   */
  void CalculateVEVSplittings(const std::vector<double> &point,
                              const double &T);

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
   * @param guess is the initial guess for the minimum
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
      std::vector<double> guess,
      std::function<std::vector<double>(std::vector<double>)> &df,
      std::function<std::vector<std::vector<double>>(std::vector<double>)>
          &Hessian,
      double error            = 1e-4,
      double const_multiplier = 1e-2,
      int maxiter             = 100);

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
 * @brief Numerical method to calculate the
 * gradient of a function f using finite differences method.
 *
 * This method is used while BSMPT is not able to
 * calculate the potential derivative analytically. We used the 4th order
 * method
 *
 * \f$\frac{\partial f}{\partial \phi_i} = \frac{1}{12
 * \epsilon}\left(-f(\dots ,\vec{\phi}_i + 2  \epsilon ) + 8 f(\dots
 * ,\vec{\phi}_i + \epsilon )- 8 f(\dots ,\vec{\phi}_i - \epsilon ) +
 * f(\dots ,\vec{\phi}_i - 2  \epsilon )\right)\f$
 *
 * where \f$ \epsilon \f$ is a small step.
 *
 * @param phi Where we want to calculate the gradient
 * @param f function
 * @param eps Size of finite differences step
 * @param dim Dimensions of the VEV space (or dimensions of V argument)
 * @return std::vector<double> The \f$ dim \times 1 \f$ gradient of V taken at
 * phi
 */
std::vector<double>
NablaNumerical(const std::vector<double> &phi,
               const std::function<double(std::vector<double>)> &f,
               const double &eps,
               const int &dim);

/**
 * @brief Numerical method to calculate the
 * hessian matrix of a function f using finite differences method.
 *
 * \f$\frac{\partial^2 f}{\partial \phi_i \phi_j} = \frac{1}{4
 * \epsilon^2}\left(V(\dots, \vec{\phi}_i + \epsilon , \vec{\phi}_j +
 * \epsilon) - f(\dots, \vec{\phi}_i - \epsilon , \vec{\phi}_j +
 * \epsilon) - f(\dots, \vec{\phi}_i + \epsilon , \vec{\phi}_j -
 * \epsilon) + f(\dots, \vec{\phi}_i - \epsilon , \vec{\phi}_j -
 * \epsilon) \right)\f$
 *
 * where \f$ \epsilon \f$ is a small step.
 *
 * @param phi Where we want to calculate the Hessian matrix
 * @param f Potential (or other function)
 * @param eps Size of finite differences step
 * @param dim Dimensions of the VEV space (or dimensions of f argument)
 * @return std::vector<std::vector<double>> The \f$ dim \times \dim \f$
 *  hessian matrix of f taken at phi
 */
std::vector<std::vector<double>>
HessianNumerical(const std::vector<double> &phi,
                 const std::function<double(std::vector<double>)> &f,
                 const double &eps,
                 const int &dim);

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
 * @param min_false
 * @param min_true
 * @param npoints
 * @return npoints long vector of steps on connecting line between min_false and
 * min_true
 */
std::vector<std::vector<double>>
Create1DimGrid(const std::vector<double> &min_false,
               const std::vector<double> &min_true,
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
 * @brief NoSignFlip project Minimum to positive quadrant in vev space
 * @param a Minimum
 */
void NoSignFlip(Minimum &a);

/**
 * @brief NoSignFlip project Minimum to positive quadrant in vev space
 * @param v std::vector<double>
 */
void NoSignFlip(std::vector<double> &v);

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
  Phase(std::vector<double> phase_start,
        double initialT,
        double finalT,
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
  Phase(double initialT,
        std::vector<double> phase_start,
        double LowT,
        double HighT,
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
  Phase(std::vector<double> phase_start,
        double initialT,
        double finalT,
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
  Phase(double initialT,
        std::vector<double> phase_start,
        double LowT,
        double HighT,
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
  double crit_temp        = -1;
  std::string crit_status = "not_set";

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
  std::string status_vacuum = "not_set";

  /**
   * @brief coexisting phases status code = success, no_coex_pairs
   */
  std::string status_coex_pairs = "not_set";

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