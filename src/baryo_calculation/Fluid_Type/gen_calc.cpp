/*
 * gen_calc.cpp
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

#include <BSMPT/baryo_calculation/Fluid_Type/gen_calc.h>

/**
 * @file
 */
namespace BSMPT{
namespace Baryo{

std::pair<std::vector<double> , std::vector<double>> set_up_nL_grid(
    std::size_t n_step,
    GSL_integration_mubl& container, 
    boost::any const  & classpointer)
{
    std::vector<double> arr_z;
    std::vector<double> arr_nL;

    arr_z.resize(n_step);
    arr_nL.resize(n_step);
    double wall_factor = container.getZMAX();
    double zstart   =   container.getZMAX();
    if(container.get_transport_method()== TransportMethod::top)
    {
        auto C_class = boost::any_cast<top_source>(&classpointer);
        if(not C_class) 
        {
            std::string errmsg = "boost::any_cast failed @ setting to top_source\n";
            throw std::runtime_error(errmsg);
        }
        for(std::size_t i=0;i<=n_step;i++)
        {
            double zend     =   i*wall_factor/n_step;
            arr_z[i]        =   zend;
            arr_nL[i]       =   C_class->Calc_nL(zstart , zend);
        }
    }    
    if(container.get_transport_method()== TransportMethod::bottom)
    {
        auto C_class = boost::any_cast<bot_source>(&classpointer);
        if(not C_class) 
            {
                std::string errmsg = "boost::any_cast failed @ setting to bot_source\n";
                throw std::runtime_error(errmsg);
            }
        for(std::size_t i=0;i<=n_step;i++)
        {
            double zend     =   i*wall_factor/n_step;
            arr_z[i]        =   zend;
            arr_nL[i]       =   C_class->Calc_nL(zstart , zend);
        }
    }
    if(container.get_transport_method()== TransportMethod::tau)
    {
        auto C_class = boost::any_cast<tau_source>(&classpointer);
        if(not C_class) 
            {
                std::string errmsg = "boost::any_cast failed @ setting to tau_source\n";
                throw std::runtime_error(errmsg);
            }
        for(std::size_t i=0;i<=n_step;i++)
        {
            double zend     =   i*wall_factor/n_step;
            arr_z[i]        =   zend;
            arr_nL[i]       =   C_class->Calc_nL(zstart , zend);
        }
    }
    std::pair<std::vector<double> , std::vector<double>>  res  = std::make_pair(arr_z,arr_nL);
    return res;
}

}
}
