// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/baryo_calculation/Fluid_Type/bot_source.h>
#include <BSMPT/baryo_calculation/Fluid_Type/gen_calc.h>
#include <BSMPT/baryo_calculation/Fluid_Type/gen_func_fluid.h>
#include <BSMPT/baryo_calculation/Fluid_Type/tau_source.h>
#include <BSMPT/baryo_calculation/Fluid_Type/top_source.h>

/**
 * @file
 */
namespace BSMPT
{
namespace Baryo
{

std::pair<std::vector<double>, std::vector<double>>
set_up_nL_grid(std::size_t n_step,
               GSL_integration_mubl &container,
               boost::any const &classpointer)
{
  std::vector<double> arr_z(n_step + 1);
  std::vector<double> arr_nL(n_step + 1);

  double wall_factor = container.getZMAX();
  double zstart      = container.getZMAX();
  if (container.get_transport_method() == TransportMethod::top)
  {
    auto C_class = boost::any_cast<top_source>(&classpointer);
    if (not C_class)
    {
      std::string errmsg = "boost::any_cast failed @ setting to top_source\n";
      throw std::runtime_error(errmsg);
    }
    for (std::size_t i = 0; i <= n_step; i++)
    {
      double zend  = i * wall_factor / n_step;
      arr_z.at(i)  = zend;
      arr_nL.at(i) = C_class->Calc_nL(zstart, zend);
    }
  }
  if (container.get_transport_method() == TransportMethod::bottom)
  {
    auto C_class = boost::any_cast<bot_source>(&classpointer);
    if (not C_class)
    {
      std::string errmsg = "boost::any_cast failed @ setting to bot_source\n";
      throw std::runtime_error(errmsg);
    }
    for (std::size_t i = 0; i <= n_step; i++)
    {
      double zend  = i * wall_factor / n_step;
      arr_z.at(i)  = zend;
      arr_nL.at(i) = C_class->Calc_nL(zstart, zend);
    }
  }
  if (container.get_transport_method() == TransportMethod::tau)
  {
    auto C_class = boost::any_cast<tau_source>(&classpointer);
    if (not C_class)
    {
      std::string errmsg = "boost::any_cast failed @ setting to tau_source\n";
      throw std::runtime_error(errmsg);
    }
    for (std::size_t i = 0; i <= n_step; i++)
    {
      double zend  = i * wall_factor / n_step;
      arr_z.at(i)  = zend;
      arr_nL.at(i) = C_class->Calc_nL(zstart, zend);
    }
  }
  std::pair<std::vector<double>, std::vector<double>> res =
      std::make_pair(arr_z, arr_nL);
  return res;
}

} // namespace Baryo
} // namespace BSMPT
