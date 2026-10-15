#pragma once

#include <BSMPT/baryo_fhck/Kernels.h>
#include <BSMPT/baryo_fhck/LocalCalcEta.h>
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
 * @brief Truncation schemes
 *
 */
enum class TruncationScheme
{
  Const,
  MinusVw,
  Variance
};

/**
 * @brief Transform Truncation schemes to strings
 *
 */
const std::unordered_map<TruncationScheme, std::string>
    TruncationSchemeToString{{TruncationScheme::Const, "R = const"},
                             {TruncationScheme::MinusVw, "R = -vw"},
                             {TruncationScheme::Variance, "Variance"}};

/**
 * @brief Types of vev profiles
 *
 */
enum class VevProfileMode
{
  Unset,
  Kink,
  FieldEquation
};

/**
 * @brief Transform VevProfileMode to strings
 *
 */
const std::unordered_map<VevProfileMode, std::string> VevProfileModeToString{
    {VevProfileMode::Unset, "unset"},
    {VevProfileMode::Kink, "Kink"},
    {VevProfileMode::FieldEquation, "FieldEquation"}};

/**
 * @brief Status of the Transport Model class
 *
 */
enum class TransportModelStatus
{
  Unset,
  Failed,
  Success
};

class TransportModel
{
private:
  /**
   * @brief modelPointer for the used parameter point
   */
  std::shared_ptr<Class_Potential_Origin> modelPointer;

  /**
   * @brief CalculateEtaInterface implementation with local minimizer
   *
   */
  std::shared_ptr<LocalCalcEta> localCalcEta;

  /**
   * @brief False vacuum
   *
   */
  std::vector<double> FalseVacuum;

  /**
   * @brief True vacuum
   *
   */
  std::vector<double> TrueVacuum;

  /**
   * @brief Empty vacuum
   *
   */
  std::vector<double> EmptyVacuum;

  /**
   * @brief Coex phase object
   *
   */
  std::shared_ptr<CoexPhases> CoexPhase;

  /**
   * @brief Splines to store the variation of the quark masses
   *
   */
  std::vector<tk::spline> QuarkMassesRe, QuarkMassesIm;

public:
  /**
   * @brief pointer to the vacuum profile solver
   *
   */
  std::shared_ptr<VacuumProfileNS::VacuumProfile> vacuumprofile;

  /**
   * @brief The standard 4 particle that appear in the transport equations
   *
   */
  std::vector<ParticleType> involvedparticles = {ParticleType::LeftFermion,
                                                 ParticleType::RightFermion,
                                                 ParticleType::LeftFermion,
                                                 ParticleType::Boson};
  /**
   * @brief Wall thickness
   *
   */
  double Lw = -1;

  /**
   * @brief Transition temperature
   *
   */
  double Tstar;

  /**
   * @brief Number of points per \f$ L_w \f$ when constructiong the mass spline.
   *
   */
  double PointsPerLwMassSpline = 100.;

  /**
   * @brief Bubble wall velocity
   *
   */
  double vwall;

  /**
   * @brief Status of the TransportModel class
   *
   */
  TransportModelStatus status = TransportModelStatus::Unset;

  /**
   * @brief Mode of the vev profile.
   * - Kink -> Kink solution
   * - FieldEquation -> Calculate the vev profile
   * using the tunneling solution
   */
  VevProfileMode VevProfile = VevProfileMode::Unset;

  /**
   * @brief Constant \f$ R \f$
   *
   */
  double truncationR;

  /**
   * @brief Truncation scheme used in the transport equations
   *
   */
  TruncationScheme truncationscheme = TruncationScheme::MinusVw;

  /**
   * @brief Default constructor
   */
  TransportModel() {};

  /**
   * @brief Construct a new Transport Model object
   *
   * @param pointer_in  Model pointer
   * @param TrueVacuum_In True Vacum
   * @param FalseVacuum_In False Vacuum
   * @param vwall_in Wall velocity
   * @param Tstar_in  Transition temperature
   * @param VevProfile_In Solver mode. Default: field EOM solution
   * @param truncationscheme_in Truncation scheme to be used
   * @param truncationR_in Truncation constant \f$ R \f$
   */
  TransportModel(
      const std::shared_ptr<Class_Potential_Origin> &pointer_in,
      const std::vector<double> TrueVacuum_In,
      const std::vector<double> FalseVacuum_In,
      const double &vwall_in,
      const double &Tstar_in,
      const VevProfileMode &VevProfile_In = VevProfileMode::FieldEquation,
      const TruncationScheme &truncationscheme_in = TruncationScheme::MinusVw,
      const double &truncationR_in                = 0);

  /**
   * @brief Construct a new Transport Model object
   *
   * @param pointer_in Model pointer
   * @param CoexPhase_in Coexphase pointer
   * @param vwall_in Wall velocity
   * @param Tstar_in Transition temperature
   * @param VevProfile_In Solver mode.
   * @param truncationscheme_in Truncation scheme to be used
   * @param truncationR_in Truncation constant \f$ R \f$
   */
  TransportModel(
      const std::shared_ptr<Class_Potential_Origin> &pointer_in,
      const std::shared_ptr<CoexPhases> &CoexPhase_in,
      const double &vwall_in,
      const double &Tstar_in,
      const VevProfileMode &VevProfile_In = VevProfileMode::FieldEquation,
      const TruncationScheme &truncationscheme_in = TruncationScheme::MinusVw,
      const double &truncationR_in                = 0);

  /**
   * @brief Create the VEV vectors
   *
   */
  virtual void Initialize();

  /**
   * @brief Set EtaInterface obejct to use BSMPTv2 functions.
   *
   */
  void SetEtaInterface();

  /**
   * @brief Calculates the VEV \f$ \vec{v} \f$ and derivatives \f$
   * \drac{d^n\vec{v}}{dz^n} \f$ at the positions z.
   *
   * @param z
   * @param diff 0 = vev, 1 = \f$ \drac{d\vec{v}}{dz} \f$, 2 = \f$
   * \drac{d^2\vec{v}}{dz^2} \f$
   * @return std::vector<double> result
   */
  std::vector<double> Vev(const double &z, const int &diff = 0);

  /**
   * @brief Generate the quark masses splines
   *
   * @param zList z distribution where we interpolate the fermion mass
   */

  /**
   * @brief Generate the quark masses splines
   *
   * @param zList z distribution where we interpolate the fermion mass
   * @param MakeTopMassPlot Make the plots or not
   */
  virtual void GenerateFermionMass(const std::vector<double> &zList,
                                   const bool &MakeTopMassPlot = false);

  /**
   * @brief Get the Fermion Mass object Calculate the fermion mass and its
   * derivatives
   *
   * @param z distance to the bubble wall
   * @param fermion which fermion, 0 = most massive
   * @param x2 \f$ x^2 \f$
   * @param x2prime \f$ x'^2 \f$
   * @param thetaprime \f$ \theta' \f$
   * @param theta2prime \f$ \theta'' \f$
   */
  virtual void GetFermionRatio(const double &z,
                              const size_t &fermion,
                              double &x2,
                              double &x2prime,
                              double &thetaprime,
                              double &theta2prime);
  /**
   * @brief Calculate the W boson mass/temperature. (same code as BSMPTv2)
   *
   * @param vev VEV
   * @return double W boson mass/temperature
   */
  virtual double GetWRatio(const double &z);

  /**
   * @brief This function calculates the EW breaking VEV from all contributing
   * field configurations.
   *
   * @param z distance from the bubble wall
   */
  virtual double EWSBVEV(const double &z);

  virtual ~TransportModel() {};
};

} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT