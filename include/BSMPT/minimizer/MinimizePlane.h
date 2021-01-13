/*
 * MinimizePlane.h
 *
 *  Created on: Jan 16, 2019
 *      Author: basler
 */

#ifndef SRC_MINIMIZER_MINIMIZEPLANE_H_
#define SRC_MINIMIZER_MINIMIZEPLANE_H_

/**
 * @file
 */

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <gsl/gsl_vector_double.h> // for gsl_vector
#include <memory>                  // for shared_ptr
#include <vector>                  // for vector

namespace BSMPT
{

class Class_Potential_Origin;
namespace Minimizer
{

/**
 * @brief The MinimizePlaneReturn struct contains the value of the potential at
 * the minimum and the minimum
 */
struct MinimizePlaneReturn
{
  double PotVal;
  std::vector<double> Minimum;
};

/**
 * @brief The GSLPlaneReturn struct
 */
struct GSLPlaneReturn
{
  /**
   * @brief vc The value of the EW VEV in the minimum
   */
  double vc;
  /**
   * @brief PotVal The value of the potential in the minimum
   */
  double PotVal;
  /**
   * @brief Minimum The electroweak minimum
   */
  std::vector<double> Minimum;
  /**
   * @brief StatusFlag A flag which is false if something went wrong
   */
  bool StatusFlag;
};

/**
 * Struct which will store all necessary information to translate between the
 * nVEV coordinates of the potential and the nVEV -1 coordinates of the plane.
 */
struct PointerContainerMinPlane
{
  /**
   * Pointer to the model class
   */
  std::shared_ptr<Class_Potential_Origin> modelPointer;

  /**
   * Vacuum of the broken phase
   */
  std::vector<double> VEVBroken;
  /**
   * Vacuum of the symmetric phase
   */
  std::vector<double> VEVSymmetric;
  /**
   * Intersection point between the plane and the connection vector of the
   * symmetric and broken phase
   */
  std::vector<double> Point;
  /**
   * Connection vector between the broken and the symmmetric phase. This is also
   * the normal vector of the plane.
   */
  std::vector<double> normalvector;
  /**
   * Temerpature at which the point should be minimised
   */
  double Temp;
  /**
   * Finds the first component of the connection vector ! = 0. In this component
   * the value of the nVEV dimensional vacuum configuration will be calculated
   * from the condition (vev - point) * normalvector = 0
   */
  std::size_t Index;
  /**
   * Number of vevs
   */
  std::size_t nVEV;
};

/**
 * Calculates the minimum of a potential on a plane. For this the normal vector
 * of the plane is calculated as the connection vector between the symmetric and
 * the broken phase. At a given point the plane normal to the connection vector
 * is then calculated and the potential is minimised along this plane.
 * @param basepoint Parameter point at which the plane and the connection
 * between the symmetric and the broken minimum should be calculated
 * @param VEVSymmetric Symmetric minimum
 * @param VEVBroken Broken minimum
 * @param Model Decides which model should be used through FChoose
 * @param par Inputparameters for the parameterpoint
 * @param parCT Counterterm parameters for the parameterpoint
 * @param Temp Temperature at which the minimum should be calculated
 * @return MinimizePlaneReturn struct which has the minimum and the value of the
 * potential
 */
MinimizePlaneReturn
MinimizePlane(const std::vector<double> &basepoint,
              const std::vector<double> &VEVSymmetric,
              const std::vector<double> &VEVBroken,
              const ModelID::ModelIDs &Model,
              const std::vector<double> &par,
              const std::vector<double> &parCT,
              const double &Temp,
              const int &WhichMinimizer = WhichMinimizerDefault);
/**
 * Calculates the minimum of a potential on a plane. For this the normal vector
 * of the plane is calculated as the connection vector between the symmetric and
 * the broken phase. At a given point the plane normal to the connection vector
 * is then calculated and the potential is minimised along this plane.
 * @param basepoint Parameter point at which the plane and the connection
 * between the symmetric and the broken minimum should be calculated
 * @param VEVSymmetric Symmetric minimum
 * @param VEVBroken Broken minimum
 * @param modelPointer Pointer to the corresponding model and parameter point
 * through Class_Potential_Origin
 * @param Temp Temperature at which the minimum should be calculated
 * @return MinimizePlaneReturn struct which has the minimum and the value of the
 * potential
 */
MinimizePlaneReturn
MinimizePlane(const std::vector<double> &basepoint,
              const std::vector<double> &VEVSymmetric,
              const std::vector<double> &VEVBroken,
              const std::shared_ptr<Class_Potential_Origin> &modelPointer,
              const double &Temp,
              const int &WhichMinimizer = WhichMinimizerDefault);

/**
 * Transform from nVEV -1 coordinates used in the plane minimisation to the nVEV
 * coordinates used to evaluate the potential.
 * @param vMinTilde A point in the nVEV -1 dimensions which were used to
 * minimise the potential in the plane
 * @param params the PointerContainerMinPlane struct which containts the
 * potential and the index at which the vector has to be extended.
 * @return The point vMinTilde will be extended by the missing component to be
 * in the nVEV dimensional space which can be used to calculate the potential
 */
std::vector<double>
TransformCoordinates(const std::vector<double> &vMinTilde,
                     const struct PointerContainerMinPlane &params);

/**
 * Calculates the value of the effective potential at the vev v and temperature
 * p->Temp for the gsl interface
 * @param v VEV configuration at which potential should be evaluated
 * @param p Pointer to a PointerContainerMinPlane struct which contains the
 * model information
 */
double
GSL_VEFF_Minimize_Plane(const gsl_vector *v, void *p);

/**
 * Uses the GSL minimisation routines to find the next local minimum from a
 * given point
 * @param p GSL_params struct which containts the information about the model
 * and the potential
 * @param sol Vector to store the solution
 * @param start Starting point from where to look for the next local minimum
 * @returns The final status of the gsl minimization process.
 */
int
GSL_Minimize_Plane_From_S_gen_all(const struct PointerContainerMinPlane &p,
                                  std::vector<double> &sol,
                                  const std::vector<double> &start);

/**
 * Minimise the Potential from different random starting points and choose the
 * local minimum with the deepest potential value as the candidate for the
 * global minimum
 * @param params PointerContainerMinPlane which containts the model information
 * and potential
 * @param seed Seed for generating the random starting points
 * @param saveAllMinima Matrix in which all found local minima are saved
 * @param MinSol Number of local minima to look out for
 * @return A GSLPlaneReturn with a boolean being true if a candidate for the
 * global minimum is found and false otherwise and a vector in which the
 * candidate for the global minimum is stored
 */
GSLPlaneReturn
GSL_Minimize_Plane_gen_all(const struct PointerContainerMinPlane &params,
                           std::size_t seed,
                           std::vector<std::vector<double>> &saveAllMinima,
                           std::size_t MinSol);

/**
 * Minimise the Potential from different random starting points and choose the
 * local minimum with the deepest potential value as the candidate for the
 * global minimum
 * @param params PointerContainerMinPlane which containts the model information
 * and potential
 * @param seed Seed for generating the random starting points
 * @param MinSol Number of local minima to look out for
 * @return A GSLPlaneReturn with a boolean being true if a candidate for the
 * global minimum is found and false otherwise and a vector in which the
 * candidate for the global minimum is stored
 */
GSLPlaneReturn
GSL_Minimize_Plane_gen_all(const struct PointerContainerMinPlane &params,
                           std::size_t seed,
                           std::size_t MinSol);

/**
 * Minimise the Potential from different random starting points and look for 20
 * local minima. Choose the local minimum with the deepest potential value as
 * the candidate for the global minimum
 * @param params PointerContainerMinPlane which containts the model information
 * and potential
 * @param seed Seed for generating the random starting points
 * @return A GSLPlaneReturn with a boolean being true if a candidate for the
 * global minimum is found and false otherwise and a vector in which the
 * candidate for the global minimum is stored
 */
GSLPlaneReturn
GSL_Minimize_Plane_gen_all(const struct PointerContainerMinPlane &params,
                           std::size_t seed);

} // namespace Minimizer
} // namespace BSMPT
#endif /* SRC_MINIMIZER_MINIMIZEPLANE_H_ */
