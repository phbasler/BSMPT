// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#include "stdexcept"
#include <BSMPT/baryo_calculation/CalculateEtaInterface.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <algorithm>
#include <boost/algorithm/string/predicate.hpp>
#include <fstream>

#include <BSMPT/baryo_calculation/Fluid_Type/bot_source.h>
#include <BSMPT/baryo_calculation/Fluid_Type/tau_source.h>
#include <BSMPT/baryo_calculation/Fluid_Type/top_source.h>

namespace BSMPT
{
namespace Baryo
{

std::pair<std::vector<bool>, int>
CalculateEtaInterface::ReadConfigFile(const std::string &file) const
{
  std::pair<std::vector<bool>, int> res;
  std::string line;
  std::ifstream configFile(file);
  if (not configFile.good())
  {
    std::string errmsg("Can not open file: " + file + " in function " +
                       __func__);
    throw std::runtime_error(errmsg);
  }
  std::vector<bool> method_transport_local;
  int bot_mass_flag_local = 0;
  while (std::getline(configFile, line))
  {
    if (boost::starts_with(line, "#") or line.empty()) continue;
    std::for_each(line.begin(), line.end(), [](char &c) { c = ::tolower(c); });
    if (boost::starts_with(line, "include"))
    {
      std::stringstream ss(line);
      std::vector<std::string> words;
      std::string tmp;
      while (ss >> tmp)
        words.push_back(tmp);
      if (words.size() < 2)
      {
        std::string errmsg = "One of the settings in the EWBG config files is "
                             "not configured correctly.";
        throw std::runtime_error(errmsg);
      }
      if (words.at(1) != "yes" and words.at(1) != "no")
      {
        std::string errmsg = "One of the settings for the EWBG config file is "
                             "not set to yes or no. Please change this.";
        throw std::runtime_error(errmsg);
      }
      method_transport_local.push_back(words.at(1) == "yes");
    }

    if (boost::starts_with(line, "massive"))
    {
      std::stringstream ss(line);
      std::vector<std::string> words;
      std::string tmp;
      while (ss >> tmp)
        words.push_back(tmp);
      if (words.size() < 2)
      {
        std::string errmsg = "The setting for the bottom mass flag is missing.";
        throw std::runtime_error(errmsg);
      }
      if (words.at(1) != "yes" and words.at(1) != "no")
      {
        std::string errmsg = "The setting for the bottom mass flag is not set "
                             "correctly. Please change it to yes or no.";
        throw std::runtime_error(errmsg);
      }
      if (words.at(1) == "yes") bot_mass_flag_local = 1;
    }
  }

  configFile.close();

  res.first  = method_transport_local;
  res.second = bot_mass_flag_local;

  return res;
}

CalculateEtaInterface::CalculateEtaInterface(
    const std::pair<std::vector<bool>, int> &config)
    : CalculateEtaInterface(config, GetSMConstants())
{
}

CalculateEtaInterface::CalculateEtaInterface(
    const std::pair<std::vector<bool>, int> &config,
    const ISMConstants &smConstants)
    : method_transport{config.first}
    , C_eta(smConstants)
    , bot_mass_flag{config.second}
    , SMConstants{smConstants}

{
  if (config.first.size() != 5)
  {
    std::string errmsg = "Warning: ";
    errmsg += __func__;
    errmsg += " expects a vector of length 5 but only received ";
    errmsg += std::to_string(config.first.size());
    errmsg += ".";
    throw std::runtime_error(errmsg);
  }
}

CalculateEtaInterface::CalculateEtaInterface(const std::string &file)
    : CalculateEtaInterface(file, GetSMConstants())
{
}

CalculateEtaInterface::CalculateEtaInterface(const std::string &file,
                                             const ISMConstants &smConstants)
    : CalculateEtaInterface(ReadConfigFile(file), smConstants)
{
}

CalculateEtaInterface::CalculateEtaInterface(
    const std::vector<bool> &method_input,
    const int &bot_mass_flag_in)
    : CalculateEtaInterface(method_input, bot_mass_flag_in, GetSMConstants())
{
}

CalculateEtaInterface::CalculateEtaInterface(
    const std::vector<bool> &method_input,
    const int &bot_mass_flag_in,
    const ISMConstants &smConstants)
    : CalculateEtaInterface(std::make_pair(method_input, bot_mass_flag_in),
                            smConstants)
{
}

CalculateEtaInterface::~CalculateEtaInterface()
{
  // TODO Auto-generated destructor stub
}

std::vector<std::string> CalculateEtaInterface::legend() const
{
  std::vector<std::string> etaLegend;
  if (method_transport.at(0)) etaLegend.push_back("eta_TopOnly");
  if (method_transport.at(1)) etaLegend.push_back("eta_TopBot");
  if (method_transport.at(2)) etaLegend.push_back("eta_TopBotTau");
  if (method_transport.at(3)) etaLegend.push_back("eta_PlasmaVelocities");
  if (method_transport.at(4))
    etaLegend.push_back("eta_PlasmaVelocitiesReplaced");

  return etaLegend;
}

void CalculateEtaInterface::setNumerics(
    const double &vw_input,
    const std::vector<double> &vev_critical_input,
    const std::vector<double> &vev_symmetric_input,
    const double &TC_input,
    std::shared_ptr<Class_Potential_Origin> &modelPointer_input,
    const int &WhichMinimizer)
{
  AdditionalBaryoArgs optionalArgs(false);
  setNumerics(vw_input,
              vev_critical_input,
              vev_symmetric_input,
              TC_input,
              modelPointer_input,
              optionalArgs,
              WhichMinimizer);
}

void CalculateEtaInterface::setNumerics(
    const double &vw_input,
    const std::vector<double> &vev_critical_input,
    const std::vector<double> &vev_symmetric_input,
    const double &TC_input,
    std::shared_ptr<Class_Potential_Origin> &modelPointer_input,
    const AdditionalBaryoArgs &AddBaryoArgs,
    const int &WhichMinimizer)
{
  vw            = vw_input;
  vev_critical  = vev_critical_input;
  vev_symmetric = vev_symmetric_input;
  TC            = TC_input;
  modelPointer  = modelPointer_input;
  if (modelPointer->get_Model() != ModelID::ModelIDs::C2HDM)
  {
    throw std::runtime_error(
        "Baryogenesis is only implemented for the C2HDM at the moment.");
  }
  GSL_integration_mubl_container.init(
      vw, vev_critical, vev_symmetric, TC, modelPointer, WhichMinimizer);
  if (AddBaryoArgs.Used)
  {
    // Work with the additional arguments
  }
}

void CalculateEtaInterface::setvw(double vw_in)
{
  vw = vw_in;
  GSL_integration_mubl_container.setvw(vw_in);
}

std::vector<double> CalculateEtaInterface::CalcEta(
    const double &vw_input,
    const std::vector<double> &vev_critical_input,
    const std::vector<double> &vev_symmetric_input,
    const double &TC_input,
    std::shared_ptr<Class_Potential_Origin> &modelPointer_input,
    const AdditionalBaryoArgs &AddBaryoArgs,
    const int &WhichMinimizer)
{
  setNumerics(vw_input,
              vev_critical_input,
              vev_symmetric_input,
              TC_input,
              modelPointer_input,
              AddBaryoArgs,
              WhichMinimizer);
  return CalcEta();
}

std::vector<double> CalculateEtaInterface::CalcEta(
    const double &vw_input,
    const std::vector<double> &vev_critical_input,
    const std::vector<double> &vev_symmetric_input,
    const double &TC_input,
    std::shared_ptr<Class_Potential_Origin> &modelPointer_input,
    const int &WhichMinimizer)
{
  AdditionalBaryoArgs OptionalArgs{false};
  setNumerics(vw_input,
              vev_critical_input,
              vev_symmetric_input,
              TC_input,
              modelPointer_input,
              OptionalArgs,
              WhichMinimizer);
  return CalcEta();
}

std::vector<double> CalculateEtaInterface::CalcEta()
{
  std::vector<double> eta;
  if (method_transport.at(0))
  {
    GSL_integration_mubl_container.set_transport_method(TransportMethod::top);
    top_source C_top(SMConstants);
    C_top.set_class(bot_mass_flag,
                    GSL_integration_mubl_container,
                    Calc_Gam_inp,
                    Calc_Scp_inp,
                    Calc_kappa_inp);
    auto arr_nL = set_up_nL_grid(n_step, GSL_integration_mubl_container, C_top);
    C_eta.set_class(arr_nL, TC, vw);
    eta.push_back(
        Nintegrate_eta(C_eta, 0, GSL_integration_mubl_container.getZMAX()));
  }
  if (method_transport.at(1))
  {
    GSL_integration_mubl_container.set_transport_method(
        TransportMethod::bottom);
    bot_source C_bot(SMConstants);
    C_bot.set_class(bot_mass_flag,
                    GSL_integration_mubl_container,
                    Calc_Gam_inp,
                    Calc_Scp_inp,
                    Calc_kappa_inp);
    auto arr_nL = set_up_nL_grid(n_step, GSL_integration_mubl_container, C_bot);
    C_eta.set_class(arr_nL, TC, vw);
    eta.push_back(
        Nintegrate_eta(C_eta, 0, GSL_integration_mubl_container.getZMAX()));
  }
  if (method_transport.at(2))
  {
    GSL_integration_mubl_container.set_transport_method(TransportMethod::tau);
    tau_source C_tau(SMConstants);
    C_tau.set_class(bot_mass_flag,
                    GSL_integration_mubl_container,
                    Calc_Gam_inp,
                    Calc_Scp_inp,
                    Calc_kappa_inp);
    auto arr_nL = set_up_nL_grid(n_step, GSL_integration_mubl_container, C_tau);
    C_eta.set_class(arr_nL, TC, vw);
    eta.push_back(
        Nintegrate_eta(C_eta, 0, GSL_integration_mubl_container.getZMAX()));
  }
  if (method_transport.at(3))
  {
    GSL_integration_mubl_container.setUseVelocityTransportEquations(true);
    eta.push_back(Integrate_mubl_interpolated(GSL_integration_mubl_container));
  }
  if (method_transport.at(4))
  {
    GSL_integration_mubl_container.setUseVelocityTransportEquations(false);
    eta.push_back(Integrate_mubl_interpolated(GSL_integration_mubl_container));
  }
  return eta;
}

double CalculateEtaInterface::getLW() const
{
  return GSL_integration_mubl_container.getLW();
}
Calc_Gam_M CalculateEtaInterface::get_class_CalcGamM() const
{
  return Calc_Gam_inp;
}
Calc_Scp CalculateEtaInterface::get_class_Scp() const
{
  return Calc_Scp_inp;
}
Calc_kappa_t CalculateEtaInterface::get_class_kappa() const
{
  return Calc_kappa_inp;
}

} // namespace Baryo
} // namespace BSMPT
