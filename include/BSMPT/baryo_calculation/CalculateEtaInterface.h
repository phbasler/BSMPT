/*
 * CalculateEtaInterface.h
 *
 *  Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file
 */

#ifndef SRC_BARYO_CALCULATION_CALCULATEETAINTERFACE_H_
#define SRC_BARYO_CALCULATION_CALCULATEETAINTERFACE_H_

#include <BSMPT/baryo_calculation/Fluid_Type/gen_calc.h>
#include <BSMPT/baryo_calculation/Fluid_Type/gen_func_fluid.h>
#include <BSMPT/baryo_calculation/transport_equations.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <string>
#include <vector>

namespace BSMPT
{
namespace Baryo
{

struct AdditionalBaryoArgs
{
  const bool Used{true};

  AdditionalBaryoArgs(bool SetUsed = true) : Used{SetUsed} {}
};

/**
 * This is an interface class to call the different EWBG methods
 */
class CalculateEtaInterface
{
protected:
  /**
   * Bool vector to tell which EWBG method is used and which is not. At the
   * moment this is \n method_transport(0) --> top only included in transport
   * equations \n method_transport(1) --> top+bot included in transport
   * equations \n method_transport(2) --> top+bot+tau included in transport
   * equations \n method_transport(3) --> FH Ansatz with plasma velocities \n
   * method_transport(4) --> FH Ansatz with plasma velocities replaced by the
   * second derivatives
   */
  std::vector<bool> method_transport;

  /**
   * the bubble wall velocity
   */
  double vw;
  /**
   * The critical temperature
   */
  double TC;
  /**
   * The vev in the broken phase at the critical temperature
   */
  std::vector<double> vev_critical;
  /**
   * The vev in the symmetric phase at the critical temperature
   */
  std::vector<double> vev_symmetric;
  /**
   * modelPointer for the used parameter point
   */
  std::shared_ptr<Class_Potential_Origin> modelPointer;
  // /**
  //  * struct to pass all necessary informations for the transport methods to
  //  the respective classes. Also calculates the wall thickness LW
  //  */

  // GSL_integration_mubl GSL_integration_mubl_container;

  /**
   * Initialisation of Gam_M class for the numerical evaluation of Gamma_M
   */
  Calc_Gam_M Calc_Gam_inp;
  /**
   * Initialisation of S_CP class for the numerical evaluation of S_CP
   */
  Calc_Scp Calc_Scp_inp;
  /**
   * Initialisation of kappa class for the numerical evaluation of kappa
   */
  Calc_kappa_t Calc_kappa_inp;
  /**
   * The step size used in the fluid methods
   */
  const std::size_t n_step = 50;
  /**
   * The member instance of the Calc_eta class to calculate the fluid Ansatz
   */
  Calc_eta C_eta;
  /**
   * Set if the bot is treated as massive (1) or not (2)
   */
  int bot_mass_flag;
  /**
   * struct to pass all necessary informations for the transport methods to the
   * respective classes. Also calculates the wall thickness LW
   */
  GSL_integration_mubl GSL_integration_mubl_container;

public:
  /**
   * @brief CalculateEtaInterface Initialises the class with a config pair
   * @param config config.first sets the CalculateEtaInterface::method_transport
   * and second CalculateEtaInterface::bot_mass_flag
   */
  CalculateEtaInterface(const std::pair<std::vector<bool>, int> &config);

  /**
   * Initialises the class member and sets the
   * CalculateEtaInterface::method_transport and
   * CalculateEtaInterface::bot_mass_flag
   * @param method_input Sets the CalculateEtaInterface::method_transport member
   * @param bot_mass_flag_in Sets the CalculateEtaInterface::bot_mass_flag
   * member
   */
  CalculateEtaInterface(const std::vector<bool> &method_input,
                        const int &bot_mass_flag_in);

  /**
   * Initialises the class member and sets the
   * CalculateEtaInterface::method_transport and
   * CalculateEtaInterface::bot_mass_flag with the input given in the input file
   * @param file input file to get the settings
   */
  CalculateEtaInterface(const std::string &file);

  virtual ~CalculateEtaInterface();

  /**
   * Read in the input file to get the configuration for the class
   */
  std::pair<std::vector<bool>, int>
  ReadConfigFile(const std::string &file) const;

  /**
   * Returns a string with the labels of the used EWBG methods stored in
   * CalculateEtaInterface::method_transport
   */
  std::vector<std::string>
  legend() const;
  /**
   * Sets the numerical values needed for the calculation
   * @param vw_input Sets the wall velocity CalculateEtaInterface::vw
   * @param vev_critical_input Sets the vev in the broken phase
   * CalculateEtaInterface::vev_critical
   * @param vev_symmetric_input Sets the vev in the symmetric phase
   * CalculateEtaInterface::vev_symmetric
   * @param TC_input Sets the critical temperature CalculateEtaInterface::TC
   * @param modelPointer_input Sets the used parameter point
   * CalculateEtaInterface::modlePointer
   * * @param WhichMinimizer defines which minimizers should be used
   * @return The labels of the different EWBG methods
   */
  void
  setNumerics(const double &vw_input,
              std::vector<double> &vev_critical_input,
              std::vector<double> &vev_symmetric_input,
              const double &TC_input,
              std::shared_ptr<Class_Potential_Origin> &modelPointer_input,
              const int &WhichMinimizer = Minimizer::WhichMinimizerDefault);
  /**
   * Sets the numerical values needed for the calculation
   * @param vw_input Sets the wall velocity CalculateEtaInterface::vw
   * @param vev_critical_input Sets the vev in the broken phase
   * CalculateEtaInterface::vev_critical
   * @param vev_symmetric_input Sets the vev in the symmetric phase
   * CalculateEtaInterface::vev_symmetric
   * @param TC_input Sets the critical temperature CalculateEtaInterface::TC
   * @param modelPointer_input Sets the used parameter point
   * CalculateEtaInterface::modlePointer
   * @param AddBaryoArgs Instance of AdditionalBaryoArgs used for optional
   * parameters in the different methods
   * @param WhichMinimizer defines which minimizers should be used
   * @return The labels of the different EWBG methods
   */
  void
  setNumerics(const double &vw_input,
              std::vector<double> &vev_critical_input,
              std::vector<double> &vev_symmetric_input,
              const double &TC_input,
              std::shared_ptr<Class_Potential_Origin> &modelPointer_input,
              const AdditionalBaryoArgs &AddBaryoArgs,
              const int &WhichMinimizer = Minimizer::WhichMinimizerDefault);
  /**
   * Calculates all EWBG methods turned on in
   * CalculateEtaInterface::method_transport with the numerical values set in
   * CalculateEtaInterface::setNumerics()
   * @return The results of the different EWBG methods
   */
  std::vector<double>
  CalcEta();
  /**
   * Calls the CalculateEtaInterface::setNumerics() function and then calculates
   * the EWBG methods stored in CalculateEtaInterface::method_transport
   * @param vw_input Sets the wall velocity CalculateEtaInterface::vw
   * @param vev_critical_input Sets the vev in the broken phase
   * CalculateEtaInterface::vev_critical
   * @param vev_symmetric_input Sets the vev in the symmetric phase
   * CalculateEtaInterface::vev_symmetric
   * @param TC_input Sets the critical temperature CalculateEtaInterface::TC
   * @param modelPointer_input Sets the used parameter point
   * CalculateEtaInterface::modlePointer
   * @param WhichMinimizer defines which minimizers should be used
   * @return The results of the different EWBG methods
   */
  std::vector<double>
  CalcEta(const double &vw_input,
          std::vector<double> &vev_critical_input,
          std::vector<double> &vev_symmetric_input,
          const double &TC_input,
          std::shared_ptr<Class_Potential_Origin> &modelPointer_input,
          const int &WhichMinimizer = Minimizer::WhichMinimizerDefault);

  /**
   * Calls the CalculateEtaInterface::setNumerics() function and then calculates
   * the EWBG methods stored in CalculateEtaInterface::method_transport
   * @param vw_input Sets the wall velocity CalculateEtaInterface::vw
   * @param vev_critical_input Sets the vev in the broken phase
   * CalculateEtaInterface::vev_critical
   * @param vev_symmetric_input Sets the vev in the symmetric phase
   * CalculateEtaInterface::vev_symmetric
   * @param TC_input Sets the critical temperature CalculateEtaInterface::TC
   * @param modelPointer_input Sets the used parameter point
   * CalculateEtaInterface::modlePointer
   * @param AddBaryoArgs Instance of AdditionalBaryoArgs used for optional
   * parameters in the different methods
   * @param WhichMinimizer defines which minimizers should be used
   * @return The results of the different EWBG methods
   */
  std::vector<double>
  CalcEta(const double &vw_input,
          std::vector<double> &vev_critical_input,
          std::vector<double> &vev_symmetric_input,
          const double &TC_input,
          std::shared_ptr<Class_Potential_Origin> &modelPointer_input,
          const AdditionalBaryoArgs &AddBaryoArgs,
          const int &WhichMinimizer = Minimizer::WhichMinimizerDefault);

  /**
   * Sets the wall velocity CalculateEtaInterface::vw
   * @param vw_in Input value for the wall velocity CalculateEtaInterface::vw
   */
  void
  setvw(double vw_in);

  /**
   * @return Wall thickness LW calculated through GSL_integration_mubl::init()
   */
  double
  getLW() const;

  /**
   * @brief get_class_CalcGamM
   * @return Calc_Gam_inp
   */
  Calc_Gam_M
  get_class_CalcGamM() const;
  /**
   * @brief get_class_Scp
   * @return Calc_Scp_inp
   */
  Calc_Scp
  get_class_Scp() const;
  /**
   * @brief get_class_kappa
   * @return Calc_kappa_inp
   */
  Calc_kappa_t
  get_class_kappa() const;

  /**
   * @brief getSymmetricCPViolatingPhase_top
   */
  auto
  getSymmetricCPViolatingPhase_top() const
  {
    return GSL_integration_mubl_container.getSymmetricCPViolatingPhase_top();
  }

  /**
   * @brief getBrokenCPViolatingPhase_top
   */
  auto
  getBrokenCPViolatingPhase_top() const
  {
    return GSL_integration_mubl_container.getBrokenCPViolatingPhase_top();
  }

  /**
   * @brief getSymmetricCPViolatingPhase_bot
   */
  auto
  getSymmetricCPViolatingPhase_bot() const
  {
    return GSL_integration_mubl_container.getSymmetricCPViolatingPhase_bot();
  }

  /**
   * @brief getBrokenCPViolatingPhase_bot
   */
  auto
  getBrokenCPViolatingPhase_bot() const
  {
    return GSL_integration_mubl_container.getBrokenCPViolatingPhase_bot();
  }

  /**
   * @brief getSymmetricCPViolatingPhase_tau
   */
  auto
  getSymmetricCPViolatingPhase_tau() const
  {
    return GSL_integration_mubl_container.getSymmetricCPViolatingPhase_tau();
  }

  /**
   * @brief getBrokenCPViolatingPhase_tau
   */
  auto
  getBrokenCPViolatingPhase_tau() const
  {
    return GSL_integration_mubl_container.getBrokenCPViolatingPhase_tau();
  }

  /**
   * @brief set_transport_method calls the set_transport_method of the
   * underlying GSL_integration_mubl
   * @param method
   */
  void
  set_transport_method(TransportMethod method)
  {
    GSL_integration_mubl_container.set_transport_method(method);
  }

  /**
   * @brief getGSL_integration_mubl_container
   */
  auto
  getGSL_integration_mubl_container() const
  {
    return GSL_integration_mubl_container;
  }
};

} // namespace Baryo
} // namespace BSMPT

#endif /* SRC_BARYO_CALCULATION_CALCULATEETAINTERFACE_H_ */
