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
    ShareInformationNLOPT(const std::shared_ptr<Class_Potential_Origin>& modelIn, const double& TempIn):
        model{modelIn},Temp{TempIn}{}
};

/**
 * @brief The NLOPTReturnType struct which is returned by the minimization
 */
struct NLOPTReturnType{
    std::vector<double> Minimum;
    double PotVal;
    nlopt::result NLOPTResult;
    bool Success;
    NLOPTReturnType()
        : NLOPTReturnType(std::vector<double>(),0,nlopt::result(), false) {}
    NLOPTReturnType(const std::vector<double>& MinimumIn,
                    const double& PotValIn,
                    const nlopt::result& NLOPTResultIn,
                    bool SuccessIn):
        Minimum{MinimumIn}, PotVal{PotValIn}, NLOPTResult{NLOPTResultIn},Success{SuccessIn}{}
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

/**
 * @brief FindLocalMinimum finds the local minimum near the starting point
 * @param model parameter point to use
 * @param Start Starting point for the minimization
 * @param Temp temperature for the minimization
 * @return the local minimum
 */
NLOPTReturnType FindLocalMinimum(const std::shared_ptr<Class_Potential_Origin>& model,
                                     const std::vector<double>& Start,
                                     const double& Temp);

}
}
}

#endif // MINIMIZENLOPT_H
