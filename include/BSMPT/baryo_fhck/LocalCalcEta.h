#pragma once

#include <BSMPT/models/IncludeAllModels.h>
#include <memory>
#include <optional>
#include <vector>

namespace BSMPT
{
namespace Baryo
{
namespace FHCK
{

/**
 * @brief
 *
 */
struct LocalCalcEta
{
  /**
   * @brief model pointer
   *
   */
  const std::shared_ptr<Class_Potential_Origin> modelPointer;

  /**
   * @brief dimension
   *
   */
  const size_t dim;

  /**
   * @brief temperature
   *
   */
  const double T;

  /**
   * @brief True and false vacuum
   *
   */
  const std::vector<double> TrueVacuum, FalseVacuum;

  /**
   * @brief normal vector and reduced normal vector
   *
   */
  std::vector<double> n, reducedn;

  /**
   * @brief (dim-1) x dim matrix to project grad to plane
   *
   */
  std::vector<std::vector<double>> gradprojection;

  /**
   * @brief Save tunnel path
   *
   */
  std::vector<std::vector<double>> path;
  /**
   * @brief index to be reduced
   *
   */
  size_t index;

  /**
   * @brief Calculate \f$ L_w \f$ using simply a straight line.
   *
   */
  bool StraightLineApproximation = true;

  /**
   * @brief Bubble width
   *
   */
  std::optional<double> Lw;

  /**
   * @brief Phase top and bot quark mas
   *
   */
  std::optional<double> top_theta_sym, bot_theta_sym, top_theta_brk,
      bot_theta_brk;

  /**
   * @brief Project point into plane
   *
   * @param vev point
   * @return std::vector<double> projected point
   */
  std::vector<double> Reduce(std::vector<double> vev);

  /**
   * @brief Calculate point complete coordinates
   *
   * @param vev reduced point
   * @param d plane identification \f$\vec{n}\cdot\vec{\phi} = d\f$
   * @return std::vector<double> full vector
   */
  std::vector<double> FromPlane(std::vector<double> vev, const double &d);

  /**
   * @brief Generate `gradprojection`
   *
   */
  void GradProjection();

  /**
   * @brief Construct a new Local Calc Eta object
   *
   * @param modelPointer_input model pointer
   * @param T_In Temperature
   * @param TrueVacuum_In True Vacuum
   * @param FalseVacuum_In False Vacuum
   */
  LocalCalcEta(
      const std::shared_ptr<Class_Potential_Origin> &modelPointer_input,
      const double &T_In,
      const std::vector<double> &TrueVacuum_In,
      const std::vector<double> &FalseVacuum_In);

  /**
   * @brief Calculate bubble width \f$L_w\f$
   *
   */
  void CalculateLw();

  /**
   * @brief Calculate top and bottom phase \f$\theta\f$
   *
   */
  void CalculateTheta();
};

} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT
