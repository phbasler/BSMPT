#ifndef MINIMIZENLOPT_H
#define MINIMIZENLOPT_H

#include <memory>
#include <vector>
#include <nlopt.hpp>

/**
 *@file
 */

namespace BSMPT {
class Class_Potential_Origin;
namespace Minimizer {
struct PointerContainerMinPlane;
namespace LibNLOPT {

/**
 * @brief The ShareInformationNLOPT struct which is used to pass the model information to the NLopt minimization routine
 */
struct ShareInformationNLOPT{
    std::shared_ptr<Class_Potential_Origin> model;
    double Temp;
};

/**
 * @brief The NLOPTReturnType struct which is returned by the minimization
 */
struct NLOPTReturnType{
    std::vector<double> Minimum;
    double PotVal;
    nlopt::result NLOPResult;
};

/**
 * @brief NLOPTVEff
 * @param x VEV configuration to evaluate the effective potential
 * @param grad NULL as a derivative free method is used
 * @param data pointer to ShareInformationNLOPT struct
 * @return
 */
double NLOPTVEff(
        const std::vector<double>& x,
        std::vector<double>& grad,
        void* data);

/**
 * @brief MinimizeUsingNLOPT minimizes the effective potential using the NLopt LN_COBYLA algorithm
 * @param model pointer to the parameter point to be minimized
 * @param Temp Temperature at which the potential should be minimized
 * @return A ShareInformationNLOPT with the global minimum, the potential value and the nlopt::result of the minimization
 */
NLOPTReturnType MinimizeUsingNLOPT(
        const std::shared_ptr<Class_Potential_Origin>& model,
        const double& Temp);

/**
 * @brief MinimizePlaneUsingNLOPT minimizes the effective potential in a given plane using the NLopt LN_COBYLA algorithm
 * @param params The struct containing all necessary informations
 * @return A ShareInformationNLOPT with the global minimum, the potential value and the nlopt::result of the minimization
 */
NLOPTReturnType MinimizePlaneUsingNLOPT(
        const struct PointerContainerMinPlane& params);
/**
 * @brief NLOPTVEffPlane
 * @param x VEV configuration to evaluate the effective potential
 * @param grad NULL as a derivative free method is used
 * @param data pointer to PointerContainerMinPlane struct
 * @return
 */
double NLOPTVEffPlane(const std::vector<double>& x, std::vector<double>& grad, void* data);

}
}
}

#endif // MINIMIZENLOPT_H
