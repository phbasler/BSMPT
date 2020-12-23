/*
 * KfactorsinterpolatedGSL.cpp
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


#include <BSMPT/Kfactors/KfactorsinterpolatedGSL.h>
#include <BSMPT/Kfactors/Kfactors.h>
#include <BSMPT/Kfactors/Kfactors_grid/Kfunctions_grid.h>
#include <BSMPT/Kfactors/Kfactors_grid/KtildeInterpolation.h>

#include <iostream>


#include <boost/version.hpp>
#if BOOST_VERSION >= 107200
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#else
#include <boost/math/interpolators/cubic_b_spline.hpp>
#endif



namespace BSMPT {
namespace Kfactors{

#if BOOST_VERSION >= 107200
template<typename T>
using boost_cubic_b_spline = boost::math::interpolators::cardinal_cubic_b_spline<T>;
#else
template<typename T>
using boost_cubic_b_spline = boost::math::cubic_b_spline<T>;
#endif

/**
 * @brief GSLAcclType Type used by the interpolations of the Kfunctions
 */
using GSLAcclType = std::unique_ptr<gsl_interp_accel,decltype(&gsl_interp_accel_free)>;

double CalculateNorm1(const double& T)
{
    return -std::pow(M_PI,3)*std::pow(T,2) * 2.0/3.0;
}

double CalculateNorm2(const double& msquared, const double& T, const int& s)
{
//    static const double startingValue = 0;
//    static const double StepSize = 1;

    static const boost_cubic_b_spline<double> KtildeNormalisationFermionSpline(Data::KtildeNormFermion_grid.data(),
                                                           Data::KtildeNormFermion_grid.size()
                                                           ,0,1);

    static const boost_cubic_b_spline<double> KtildeNormalisationBosonSpline(Data::KtildeNormBoson_grid.data(),
                                                                                    Data::KtildeNormBoson_grid.size()
                                                                                    ,0,1);

    double norm = -1;
    if(s==1){
        norm = KtildeNormalisationFermionSpline(msquared/std::pow(T,2));
    }
    else{
        norm = KtildeNormalisationBosonSpline(msquared/std::pow(T,2));
    }
    norm *= std::pow(T,3);
    return norm;
}

std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK1fermionGrid()
{
    auto grid= std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)>(gsl_spline2d_alloc(gsl_interp2d_bicubic, Data::msg_size, Data::Tg_size),gsl_spline2d_free);
    std::vector<double> zav(Data::msg_size*Data::Tg_size);
    for(std::size_t i=0;i<Data::msg_size;i++){
        for(std::size_t j=0;j<Data::Tg_size;j++){
            zav.at(j*Data::msg_size+i) = Data::K1p[i][j];
        }
    }
    gsl_spline2d_init(grid.get(), Data::msg.data(), Data::Tg.data(), zav.data(), Data::msg_size, Data::Tg_size);
    return grid;
}
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK1bosonGrid()
{
    auto grid= std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)>(gsl_spline2d_alloc(gsl_interp2d_bicubic, Data::msg_size, Data::Tg_size),gsl_spline2d_free);
    std::vector<double> zav(Data::msg_size*Data::Tg_size);
    for(std::size_t i=0;i<Data::msg_size;i++){
        for(std::size_t j=0;j<Data::Tg_size;j++){
            zav.at(j*Data::msg_size+i) = Data::K1m[i][j];
        }
    }
    gsl_spline2d_init(grid.get(), Data::msg.data(), Data::Tg.data(), zav.data(), Data::msg_size, Data::Tg_size);
    return grid;
}
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK2fermionGrid()
{
    auto grid= std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)>(gsl_spline2d_alloc(gsl_interp2d_bicubic, Data::msg_size, Data::Tg_size),gsl_spline2d_free);
    std::vector<double> zav(Data::msg_size*Data::Tg_size);
    for(std::size_t i=0;i<Data::msg_size;i++){
        for(std::size_t j=0;j<Data::Tg_size;j++){
            zav.at(j*Data::msg_size+i) = Data::K2p[i][j];
        }
    }
    gsl_spline2d_init(grid.get(), Data::msg.data(), Data::Tg.data(), zav.data(), Data::msg_size, Data::Tg_size);
    return grid;
}
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK4fermionGrid()
{
    auto grid= std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)>(gsl_spline2d_alloc(gsl_interp2d_bicubic, Data::msg_size, Data::Tg_size),gsl_spline2d_free);
    std::vector<double> zav(Data::msg_size*Data::Tg_size);
    for(std::size_t i=0;i<Data::msg_size;i++){
        for(std::size_t j=0;j<Data::Tg_size;j++){
            zav.at(j*Data::msg_size+i) = Data::K4p[i][j];
        }
    }
    gsl_spline2d_init(grid.get(), Data::msg.data(), Data::Tg.data(), zav.data(), Data::msg_size, Data::Tg_size);
    return grid;
}
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK4bosonGrid()
{
    auto grid= std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)>(gsl_spline2d_alloc(gsl_interp2d_bicubic, Data::msg_size, Data::Tg_size),gsl_spline2d_free);
    std::vector<double> zav(Data::msg_size*Data::Tg_size);
    for(std::size_t i=0;i<Data::msg_size;i++){
        for(std::size_t j=0;j<Data::Tg_size;j++){
            zav.at(j*Data::msg_size+i) = Data::K4m[i][j];
        }
    }
    gsl_spline2d_init(grid.get(), Data::msg.data(), Data::Tg.data(), zav.data(), Data::msg_size, Data::Tg_size);
    return grid;
}
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK5fermionGrid()
{
    auto grid= std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)>(gsl_spline2d_alloc(gsl_interp2d_bicubic, Data::msg_size, Data::Tg_size),gsl_spline2d_free);
    std::vector<double> zav(Data::msg_size*Data::Tg_size);
    for(std::size_t i=0;i<Data::msg_size;i++){
        for(std::size_t j=0;j<Data::Tg_size;j++){
            zav.at(j*Data::msg_size+i) = Data::K5p[i][j];
        }
    }
    gsl_spline2d_init(grid.get(), Data::msg.data(), Data::Tg.data(), zav.data(), Data::msg_size, Data::Tg_size);
    return grid;
}
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK5bosonGrid()
{
    auto grid= std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)>(gsl_spline2d_alloc(gsl_interp2d_bicubic, Data::msg_size, Data::Tg_size),gsl_spline2d_free);
    std::vector<double> zav(Data::msg_size*Data::Tg_size);
    for(std::size_t i=0;i<Data::msg_size;i++){
        for(std::size_t j=0;j<Data::Tg_size;j++){
            zav.at(j*Data::msg_size+i) = Data::K5m[i][j];
        }
    }
    gsl_spline2d_init(grid.get(), Data::msg.data(), Data::Tg.data(), zav.data(), Data::msg_size, Data::Tg_size);
    return grid;
}
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK6fermionGrid()
{
    auto grid= std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)>(gsl_spline2d_alloc(gsl_interp2d_bicubic, Data::msg_size, Data::Tg_size),gsl_spline2d_free);
    std::vector<double> zav(Data::msg_size*Data::Tg_size);
    for(std::size_t i=0;i<Data::msg_size;i++){
        for(std::size_t j=0;j<Data::Tg_size;j++){
            zav.at(j*Data::msg_size+i) = Data::K6p[i][j];
        }
    }
    gsl_spline2d_init(grid.get(), Data::msg.data(), Data::Tg.data(), zav.data(), Data::msg_size, Data::Tg_size);
    return grid;
}
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK8fermionGrid()
{
    auto grid= std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)>(gsl_spline2d_alloc(gsl_interp2d_bicubic, Data::msg_size, Data::Tg_size),gsl_spline2d_free);
    std::vector<double> zav(Data::msg_size*Data::Tg_size);
    for(std::size_t i=0;i<Data::msg_size;i++){
        for(std::size_t j=0;j<Data::Tg_size;j++){
            zav.at(j*Data::msg_size+i) = Data::K8p[i][j];
        }
    }
    gsl_spline2d_init(grid.get(), Data::msg.data(), Data::Tg.data(), zav.data(), Data::msg_size, Data::Tg_size);
    return grid;
}
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK9fermionGrid()
{
    auto grid= std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)>(gsl_spline2d_alloc(gsl_interp2d_bicubic, Data::msg_size, Data::Tg_size),gsl_spline2d_free);
    std::vector<double> zav(Data::msg_size*Data::Tg_size);
    for(std::size_t i=0;i<Data::msg_size;i++){
        for(std::size_t j=0;j<Data::Tg_size;j++){
            zav.at(j*Data::msg_size+i) = Data::K9p[i][j];
        }
    }
    gsl_spline2d_init(grid.get(), Data::msg.data(), Data::Tg.data(), zav.data(), Data::msg_size, Data::Tg_size);
    return grid;
}


double K1fermion(double msquared, double T)
{
    static const std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> grid = initializeK1fermionGrid();
    GSLAcclType xacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    GSLAcclType yacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    double res = gsl_spline2d_eval(grid.get(),msquared,T,xacc.get(),yacc.get());
    return res;
}

double K1boson(double msquared, double T)
{
    static const std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> grid = initializeK1bosonGrid();
    GSLAcclType xacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    GSLAcclType yacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    double res = gsl_spline2d_eval(grid.get(),msquared,T,xacc.get(),yacc.get());
    return res;
}

double K2fermion(double msquared, double T)
{
    static const std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> grid = initializeK2fermionGrid();
    GSLAcclType xacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    GSLAcclType yacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    double res = gsl_spline2d_eval(grid.get(),msquared,T,xacc.get(),yacc.get());
    return res;
}

double K4fermion(double msquared, double T)
{
    static const std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> grid = initializeK4fermionGrid();
    GSLAcclType xacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    GSLAcclType yacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    double res = gsl_spline2d_eval(grid.get(),msquared,T,xacc.get(),yacc.get());
    return res;
}

double K4boson(double msquared, double T)
{
    static const std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> grid = initializeK4bosonGrid();
    GSLAcclType xacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    GSLAcclType yacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    double res = gsl_spline2d_eval(grid.get(),msquared,T,xacc.get(),yacc.get());
    return res;
}

double K5fermion(double msquared, double T)
{
    static const std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> grid = initializeK5fermionGrid();
    GSLAcclType xacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    GSLAcclType yacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    double res = gsl_spline2d_eval(grid.get(),msquared,T,xacc.get(),yacc.get());
    return res;
}

double K5boson(double msquared, double T)
{
    static const std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> grid = initializeK5bosonGrid();
    GSLAcclType xacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    GSLAcclType yacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    double res = gsl_spline2d_eval(grid.get(),msquared,T,xacc.get(),yacc.get());
    return res;
}

double K6fermion(double msquared, double T)
{
    static const std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> grid = initializeK6fermionGrid();
    GSLAcclType xacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    GSLAcclType yacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    double res = gsl_spline2d_eval(grid.get(),msquared,T,xacc.get(),yacc.get());
    return res;
}

double K8fermion(double msquared, double T)
{
    static const std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> grid = initializeK8fermionGrid();
    GSLAcclType xacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    GSLAcclType yacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    double res = gsl_spline2d_eval(grid.get(),msquared,T,xacc.get(),yacc.get());
    return res;
}

double K9fermion(double msquared, double T)
{
    static const std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> grid = initializeK9fermionGrid();
    GSLAcclType xacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    GSLAcclType yacc(gsl_interp_accel_alloc(),gsl_interp_accel_free);
    double res = gsl_spline2d_eval(grid.get(),msquared,T,xacc.get(),yacc.get());
    return res;
}



double K1fermion_normalized(double msquared, double T)
{
	double nom = K1fermion(msquared,T);
	double denom = CalculateNorm1(T);
	return nom/denom;
}


double K1boson_normalized(double msquared, double T)
{
	return K1boson(msquared,T)/CalculateNorm1(T);
}


double K2fermion_normalized(double msquared, double T)
{
	return K2fermion(msquared,T)/CalculateNorm1(T);
}



double K4fermion_normalized(double msquared, double T)
{
	return K4fermion(msquared,T)/CalculateNorm1(T);
}



double K4boson_normalized(double msquared, double T)
{
	return K4boson(msquared,T)/CalculateNorm1(T);
}


double K5fermion_normalized(double msquared, double T)
{
	return K5fermion(msquared,T)/CalculateNorm2(msquared,T,1);
}



double K5boson_normalized(double msquared, double T)
{
	return K5boson(msquared,T)/CalculateNorm2(msquared,T,-1);
}


double K6fermion_normalized(double msquared, double T)
{
	return K6fermion(msquared,T)/CalculateNorm2(msquared,T,1);
}



double K8fermion_normalized(double msquared, double T)
{
	return K8fermion(msquared,T)/CalculateNorm2(msquared,T,1);
}



double K9fermion_normalized(double msquared, double T)
{
	return K9fermion(msquared,T)/CalculateNorm1(T);
}



}
}
