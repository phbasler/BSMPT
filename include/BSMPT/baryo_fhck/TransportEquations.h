#pragma once

#include <BSMPT/baryo_fhck/Kernels.h>
#include <BSMPT/baryo_fhck/LocalCalcEta.h>
#include <BSMPT/baryo_fhck/TransportModel.h>
#include <BSMPT/bounce_solution/action_calculation.h>
#include <BSMPT/bounce_solution/bounce_solution.h>
#include <BSMPT/gravitational_waves/gw.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/utility/NumericalDerivatives.h>
#include <BSMPT/utility/relaxation/solvde.h>
#include <BSMPT/utility/spline/spline.h>
#include <BSMPT/utility/utility.h>
#include <BSMPT/vacuum_profile/vacuum_profile.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>

namespace BSMPT
{
namespace Baryo
{
namespace FHCK
{

/**
 * @brief Status of the FHCK baryo calculation
 *
 */
enum class FHCKStatus
{
  NotSet,
  Success,
  SmallIntegrationRegion,
  UnphysicalBoundary
};

/**
 * @brief Convert FHCK baryo into strings
 *
 */
const std::unordered_map<FHCKStatus, std::string> FHCKStatusToString{
    {FHCKStatus::NotSet, "not_set"},
    {FHCKStatus::Success, "success"},
    {FHCKStatus::SmallIntegrationRegion, "small_integration_region"},
    {FHCKStatus::UnphysicalBoundary, "unphysical_boundary"}};

class TransportEquations : Difeq
{
public:
  /**
   * @brief Transport model under consideration
   */
  std::shared_ptr<TransportModel> transportmodel;

  /**
   * @brief Status of the FHCK baryo calculation
   *
   */
  FHCKStatus Status = FHCKStatus::NotSet;

  /**
   * @brief Transition temperature
   *
   */
  double Tstar;

  /**
   * @brief Number of moments used to solve transport equation
   *
   */
  std::vector<size_t> moments;

  /**
   * @brief Number of times we let the solution relax without getting better.
   *
   */
  const size_t NotBetterThreshold = 2;

  /**
   * @brief Temp var to store \f$ \eta \f$ at moment
   *
   */
  double bau;

  /**
   * @brief \f$ \eta = \frac{n_B}{n_\gamma}\f$
   *
   */
  std::vector<std::optional<double>> BAUeta;

  /**
   * @brief Moment to consider
   *
   */
  size_t moment;

  /**
   * @brief Interpolated kernel functions for different moments
   *
   */
  std::vector<tk::spline> Klf, Klb, Dlf, Dlb, Qlf, Qlb, Q8ol, Q9ol;

  /**
   * @brief Special interpolated kernel functions
   *
   */
  tk::spline Rbarb, Rbarf, K4FHf, K4FHb;

  /**
   * @brief Rel. gamma factor of the wall
   *
   */
  double gamwall;

  /**
   * @brief List of point on the z-axis.
   *
   */
  std::vector<double> zList;

  /**
   * @brief \f$ M \f$ matrix at the z-negative boundary.
   *        \f$ S \f$ vector at the z-positive boundary.
   */
  MatDoub MtildeM, MtildeP;

  /**
   * @brief \f$ S \f$ vector at the z - negative boundary.
   *        \f$ M \f$ matrix at the z - positive boundary.
   */
  VecDoub StildeM, StildeP;

  /**
   * @brief EtaInterface object. To  use BSMPTv2 functions
   *
   */
  std::shared_ptr<LocalCalcEta> localCalcEta;

  /**
   * @brief Phase of the top mass at the false vacuum
   *
   */
  std::optional<double> Theta_False;

  /**
   * @brief Phase of the top mass at the true vacuum
   *
   */
  std::optional<double> Theta_True;

  /**
   * @brief Store the solution from the relaxation method
   *
   */
  MatDoub Solution;

  /**
   * @brief Number of Weil fermions
   *
   */
  size_t nFermions = 3;

  /**
   * @brief Number of bosons
   *
   */
  size_t nBosons = 1;

  /**
   * @brief Number of particles
   *
   */
  size_t nParticles;

  /**
   * @brief Number of equations
   *
   */
  size_t nEqs;

  /**
   * @brief Number of fermionic families
   *
   */
  const double nf = 3;

  /**
   * @brief When \f$ S = 0 \f$ the \f$ \mu = \mu_0 e^{-\lambda z}\f$, by taking
   * the imaginary part we can calculate how rapidly the functions are
   * oscillating. This sets the goal number of points per cycle/wavelength
   *
   */
  double StepsPerCycle;

  /**
   * @brief When \f$ S = 0 \f$ the \f$ \mu = \mu_0 e^{-\lambda z}\f$, by taking
   * the imaginary part we can calculate how rapidly the functions are
   * oscillating. This sets the goal number of points per cycle/wavelength for a
   * low precision calculation.
   *
   */
  const double StepsPerCycleLow = 35;

  /**
   * @brief When \f$ S = 0 \f$ the \f$ \mu = \mu_0 e^{-\lambda z}\f$, by taking
   * the imaginary part we can calculate how rapidly the functions are
   * oscillating. This sets the goal number of points per cycle/wavelength for a
   * high precision calculation.
   *
   */
  const double StepsPerCycleHigh = 45;

  /**
   * @brief Used for quick calculations (should not be that relevant)
   *
   */
  const size_t DefaultNumberOfSteps = 10000;

  /**
   * @brief Number of steps in space
   *
   */
  size_t NumberOfSteps;

  /**
   * @brief Maximum number of steps in \f$ z \f$
   *
   */
  const size_t MaxNumberOfSteps = 1e5;

  /**
   * @brief When \f$ S = 0 \f$ the \f$ \mu = \mu_0 e^{-\lambda z}\f$. We use
   * this so calculate when \f$ e^{-\lambda z} = \text{LwMultiplierCutoff}\f$
   * using the slowest decaying $\lambda$ (highest negative real part)
   *
   */
  const double LwMultiplierCutoff = 1e-10;

  /**
   * @brief Used for quick calculations (should not be that relevant)
   *
   */
  const double DefaultLwMultiplier = 100;

  /**
   * @brief The integration goes from \f$ - LwMultiplier * Lw \f$ up to \f$
   * LwMultiplier * Lw \f$
   *
   */
  double LwMultiplier;

  /**
   * @brief Minimum LwMultiplier used
   *
   */
  const double MinLwMultiplier = 10.;

  /**
   * @brief Threshold for which the length of the S vector must be smaller. If
   * not then the integration region must be inscreased.
   *
   */
  const double STildeThreshold = 1e-10;

  /**
   * @brief Uncertainty threshold for accepting a BAU
   *
   */
  const double UncertaintyThreshold = 1e-2;

  /**
   * @brief Construct a new Transport Equations object
   *
   * @param model_in Transport model
   * @param pointer_in Model pointer
   * @param Tstar_in Transition temperature
   */
  TransportEquations(const std::shared_ptr<TransportModel> &model_in,
                     const double &Tstar_in,
                     const std::vector<size_t> &moments_in = {2});
  /**
   * @brief Create the VEV vectors and initalize **Ki()** and **Kfac()**
   * objects.
   *
   */
  void Initialize();

  /**
   * @brief Calculates the decaying eigenvalues to calculate how long we have to
   * go in \f$ z \f$
   *
   */
  void GenerateIntegrationSpace();

  /**
   * @brief Initialize things for the this specific moment
   *
   * @param moment
   */
  void InitializeMoment(const size_t &moment);

  /**
   * @brief Generates the spline for the <kernel_Name>
   *
   * @param kernel_Name Kernel function name
   * @param is_1D is 1D or 2D
   * @return tk::spline
   */
  tk::spline InterpolateKernel(const std::string &kernel_Name,
                               const bool is_1D);

  /**
   * @brief Load all the Kernel functions and interpolate them
   *
   *
   */
  void BuildKernelInterpolation();

  /**
   * @brief Calculate the collision matrix \f$ \deltaC^{1/2} \f$
   *
   * @param z distance to the bubble wall.
   * @return MatDoub Collision matrix
   */
  MatDoub CalculateCollisionMatrix(const double &xW,
                                   VecDoub &FermionMasses,
                                   VecDoub &BosonMasses);

  /**
   * @brief Function that calculates the Ri vector;
   *
   * @param particle which particle
   * @param k index of ponit
   * @return std::vector<double> \f$ R_i \f$
   */
  std::vector<double> calc_Ri(const size_t &particle = 0, const size_t &k = 0);

  /**
   * @brief Calculate the 2x2 submatrix of Ainv for 1 particle
   *
   * @param x mass/temperature of the particle
   * @param type type type of particle e.g. fermion/boson
   * @param particle which particle
   * @param k index in ydifeq
   * @return MatDoub of Ainv
   */
  MatDoub calc_Ainv(const double &x,
                    const ParticleType &type,
                    const size_t &particle,
                    const size_t &k = 0);

  /**
   * @brief Calculate the 2x2 submatrix of m2'B for 1 particle
   *
   * @param x mass/temperature of the particle
   * @param dx2 derivative of the squared particle mass
   * @param type type of particle e.g. fermion/boson
   * @return MatDoub of m2'B
   */
  MatDoub
  calc_m2B(const double &x, const double &dx2, const ParticleType &type);

  /**
   * @brief Calculate the souce term S for 1 particle
   *
   * @param x mass/temperature of the particle
   * @param dx2 derivative of the squared particle mass
   * @param dth derivative of theta
   * @param d2th second derivative theta
   * @param h helicity
   * @return VecDoub source term
   */
  VecDoub calc_source(const double &x,
                      const double &dx2,
                      const double &dth,
                      const double &d2th,
                      const int &h);

  /**
   * @brief Check if the boundary have enough decaying modes for the solution to
   * exist.
   */

  /**
   * @brief Check if the boundary have enough decaying modes for the solution to
   * exist.
   *
   * @param HighestNegReEigenvalue Highest real part of negative eigenvalues
   * (close to zero from below)
   * @param HighestImEigenvalue Highest magnitude of all eigenvalues with
   * negative real part
   */
  void CheckBoundary(double &HighestNegReEigenvalue,
                     double &HighestImEigenvalue);

  /**
   * @brief Check if the boundary have enough decaying modes for the solution to
   * exist.
   */
  void CheckBoundary();

  /**
   * @brief Insert a sub matrix into a block diagonal matrix
   *
   * @param full \f$ reference to the block diagonal matrix \f$
   * @param sub \f$ sub matrix to be inserted \f$
   * @param Stilde \f$ position at which we want to insert the sub matrix \f$
   */
  void InsertBlockDiagonal(MatDoub &full, MatDoub &sub, const size_t position);

  /**
   * @brief Calculate the equations
   *
   * @param z distance to the bubble wall
   * @param Mtilde \f$ A^{-1} \Gamma \f$
   * @param Stilde \f$ A^{-1}\right( \deltaC - m^2' B \left) \f$
   * @param k index of point
   */
  void Equations(const double &z,
                 MatDoub &Mtilde,
                 VecDoub &Stilde,
                 const size_t &k = 0);

  /**
   * @brief Solve the transport equation using the relaxation method.
   * Recommended.
   *
   */
  void RelaxationMethod();

  /**
   * @brief Distributes the points used in the relaxation method
   *
   * @param amplitude amplitude of points distributed
   * @param npoints number of points
   */
  void MakeDistribution(const double amplitude, const size_t npoints);

  /**
   * @brief Smatrix needed of the relaxation method to solve the ODE
   */
  void smatrix(const int k,
               const int k1,
               const int k2,
               const int jsf,
               VecInt &indexv,
               MatDoub &s,
               MatDoub &y);

  /**
   * @brief Solve the transport equation for the moment \f$ \ell \f$
   *
   * @param ell \f$ \ell \f$
   * @return double BAU
   */
  double SolveTransportEquationEll(const size_t &ell);

  /**
   * @brief Solve the transport equations.
   *
   */
  void SolveTransportEquation();

  /**
   * @brief Interpolated spharelon factor suppression
   *
   * @param z
   * @return double
   */
  double Gws(const double &z);

  /**
   * @brief \f$ \mathrm{exp}\left[-\frac{A
   * n_{f}}{2v_{w}\gamma_{w}}\int_{-\infty}^{z}\Gamma_{w s}(u)d u\right] \f$
   *
   * @param z
   * @return double
   */
  double WashoutFactor(const double &z);

  /**
   * @brief Calculate the baryonic assymetry of the Universe
   *
   */
  void CalculateBAU();

  /**
   * @brief
   *
   * @param size
   * @param Particle
   * @param MuOrU
   */
  void PrintTransportEquation(const int &size,
                              const std::string &Particle,
                              const std::string &MuOrU);
};
} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT
