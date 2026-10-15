#pragma once

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/SMparam.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/relaxation/solvde.h>
#include <BSMPT/utility/spline/spline.h>
#include <BSMPT/utility/utility.h>
#include <BSMPT/vacuum_profile/difeq_vacuum_profile.h>
#include <iostream>
#include <string>
#include <vector>

namespace BSMPT
{
namespace VacuumProfileNS
{

/**
 * @brief Status of the Vacuum Profile calculation
 *
 */
enum class VacuumProfileStatus
{
  Unset,
  Failed,
  Success
};

struct VacuumProfile
{
  /**
   * @brief Method used to solve the field equation
   * - Field - Fixed \f$ \vec\phi(z \to \infty) = \vec\phi_f \f$ and \f$
   * \vec\phi(z \to -\infty) = \vec\phi_t \f$
   * - Deriv - Fixes \frac{d\vec\phi}{dz}(|z| \to \infty) = 0
   *
   */
  ProfileSolverMode mode = ProfileSolverMode::Field;

  /**
   * @brief State of the vacuum profile
   *
   */
  VacuumProfileStatus status = VacuumProfileStatus::Unset;

  /**
   * @brief Dimension of VEV space
   *
   */
  const size_t dim;

  /**
   * @brief Number of times we let the profile relax without getting better.
   *
   */
  size_t NotBetterThreshold = 5;

  /**
   * @brief Bubble width
   *
   */
  double Lw;

  /**
   * @brief Try to solve the tunnel profile from -Lw * LwToSolve -> Lw *
   * LwToSolve
   *
   */
  double LwToSolve = 10;

  /**
   * @brief Max iterations for the relaxation method
   *
   */
  size_t itmax = 50;

  /**
   * @brief Calculates \f$ \frac{d\phi}{dz} \f$ at the edges and, if its too
   * large, increases the integration domain.
   *
   */
  const double dphiTreshold = 1.e-3;

  /**
   * @brief True and False Vacuum
   *
   */
  const std::vector<double> &TrueVacuum, &FalseVacuum;

  /**
   * @brief Potential
   *
   */
  const std::function<double(std::vector<double>)> &V;

  /**
   * @brief Gradient
   *
   */
  const std::function<std::vector<double>(std::vector<double>)> &dV;

  /**
   * @brief Hessian
   *
   */
  const std::function<std::vector<std::vector<double>>(std::vector<double>)>
      &Hessian;

  /**
   * @brief z positions corresponding to the path knots
   *
   */
  std::vector<double> z;

  /**
   * @brief Path that the solver takes
   *
   */
  MatDoub y;

  /**
   * @brief Scale of each field. Used for the solver
   *
   */
  VecDoub scalv;

  /**
   * @brief
   *
   */
  std::vector<tk::spline> splines;

  /**
   * @brief Generate the splines wih the fields
   *
   */
  void GenerateSplines();

  /**
   * @brief Calculates the vacuum profile at positon z
   *
   * @param z position
   */
  std::vector<double> GetVev(const double &zz, const int &diff = 0);

  /**
   * @brief Relocate the vacuum profile "center" to \f$ z = 0 \f$ but
   * putting
   * \f$ \max \left\{\frac{d\vec\phi}{dz}\right\}^2 \f$ at \f$ z = 0 \f$
   *
   */
  void CenterPath();

  /**
   * @brief Relocate the vacuum profile "center" to \f$ z = 0 \f$ but
   * putting
   * \f$ \max \left\{\frac{d\vec\phi}{dz}\right\}^2 \f$ at \f$ z = 0 \f$
   *
   * @param center calculated old center
   */
  void CenterPath(double &center);

  /**
   * @brief Calculate the order of the fields to pass to the numerical solver
   *
   * @param mode solver mode
   * @return VecInt reordering
   */
  VecInt Calcindexv();

  /**
   * @brief Calculate the bubble width using a rought approximation \f$ L_w =
   * \frac{v_c}{\sqrt{8 V_b}} \f$
   *
   * @return double
   */

  /**
   * @brief Calculate the bubble width using a rought approximation \f$ L_w =
   * \frac{v_c}{\sqrt{8 V_b}} \f$. Algorithm from BSMPTv2 on the straight line
   *
   * @param TrueVacuum True vacuum
   * @param FalseVacuum False vacuum
   * @param V Potential
   * @return double
   */
  static double
  CalculateWidth(const std::vector<double> &TrueVacuum,
                 const std::vector<double> &FalseVacuum,
                 const std::function<double(std::vector<double>)> &V);
  /**
   * @brief Create a Path object
   *
   * @param NumberOfSteps Number of steps
   */
  void CreatePath(const size_t &NumberOfSteps);

  /**
   * @brief Load path -> z, y
   *
   * @param z_In \f$ z \f$ positions
   * @param path_In field position at \f$ z \f$
   */
  void LoadPath(const std::vector<double> &z_In,
                const std::vector<std::vector<double>> &path_In);

  /**
   * @brief Calculate the vacuum profile
   *
   */
  void CalculateProfile();

  /**
   * @brief Construct a new Vacuum Profile object. Attemps to calculate \f$ L_w
   * \f$ for the first guess using the kink solution
   *
   * @param dim_In VEV dimension
   * @param TrueVacuum_In True Vacuum
   * @param FalseVacuum_In False Vacuum
   * @param V_In Potential
   * @param dV_In Potential Gradient
   * @param Hessian_In Potential Hessian
   */
  VacuumProfile(
      // Dimension of VEV space
      const size_t dim_In,
      // True
      const std::vector<double> &TrueVacuum_In,
      // False
      const std::vector<double> &FalseVacuum_In,
      // Potential
      const std::function<double(std::vector<double>)> &V_In,
      // Gradient
      const std::function<std::vector<double>(std::vector<double>)> &dV_In,
      // Hessian
      const std::function<std::vector<std::vector<double>>(std::vector<double>)>
          &Hessian_In);

  /**
   * @brief Construct a new Vacuum Profile object. Uses \f$ L_w \f$ for the
   * first path guess using the kink solution
   *
   * @param dim_In VEV dimension
   * @param TrueVacuum_In True Vacuum
   * @param FalseVacuum_In False Vacuum
   * @param V_In Potential
   * @param dV_In Potential Gradient
   * @param Hessian_In Potential Hessian
   */
  VacuumProfile(
      // Dimension of VEV space
      const size_t dim_In,
      // True
      const std::vector<double> &TrueVacuum_In,
      // False
      const std::vector<double> &FalseVacuum_In,
      // Potential
      const std::function<double(std::vector<double>)> &V_In,
      // Gradient
      const std::function<std::vector<double>(std::vector<double>)> &dV_In,
      // Hessian
      const std::function<std::vector<std::vector<double>>(std::vector<double>)>
          &Hessian_In,
      // Bubble width
      const double Lw_In,
      // Number of steps for the path
      const size_t NumberOfSteps = 1000);

  /**
   * @brief Construct a new Vacuum Profile object
   *
   * @param dim_In VEV dimension
   * @param TrueVacuum_In True Vacuum
   * @param FalseVacuum_In False Vacuum
   * @param V_In Potential
   * @param dV_In Potential Gradient
   * @param Hessian_In Potential Hessian
   * @param z_In List of the \f$ z \f$ position
   * @param path_In List of the field knots
   */
  VacuumProfile(
      // Dimension of VEV space
      const size_t dim_In,
      // True
      const std::vector<double> &TrueVacuum_In,
      // False
      const std::vector<double> &FalseVacuum_In,
      // Potential
      const std::function<double(std::vector<double>)> &V_In,
      // Gradient
      const std::function<std::vector<double>(std::vector<double>)> &dV_In,
      // Hessian
      const std::function<std::vector<std::vector<double>>(std::vector<double>)>
          &Hessian_In,
      // z knot list
      std::vector<double> &z_In,
      // field knot list
      std::vector<std::vector<double>> &path_In);
};

} // namespace VacuumProfileNS
} // namespace BSMPT
