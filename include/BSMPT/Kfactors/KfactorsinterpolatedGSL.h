/*
 * KfactorsinterpolatedGSL.h
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

#ifndef SRC_BARYO_CALCULATION_KFACTORSINTERPOLATEDGSL_H_
#define SRC_BARYO_CALCULATION_KFACTORSINTERPOLATEDGSL_H_


#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <cmath>
#include <boost/math/interpolators/cubic_b_spline.hpp>

/**
 * @file
 * These functions use the gsl_interp2d_bicubic interpolation routine from GSL to calculate the K-functions
 * defined in arXiv:hep-ph/0604159v2 Eq. (23). These are needed in the transport equations which are implemented
 * in transport_equations.
 * For this a the integrals (without their normalisation) were evaluated numerically in a grid and stored in the
 * header Kfunctions_grid.h
 */

namespace BSMPT{
namespace Kfactors {

/**
 * Calculates the norm for < > Integrals
 * @param T temperature at which the normalisation should be evaluated
 */
double CalculateNorm1(const double& T);

/**
 * Calculates the norm for [] Integrals
 * @param msquared m^2 value at which the normalisation should be evaluated
 * @param T temperature at which the normalisation should be evaluated
 * @param s switch for bosons (s=-1) and fermions (s=1)
 */
double CalculateNorm2(const double& msquared, const double& T, const int& s);

/**
 * @brief initializeK1fermionGrid
 * @return a unique pointer to the grid for the K_1 fermion interpolation
 */
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK1fermionGrid();
/**
 * @brief initializeK1bosonGrid
 * @return a unique pointer to the grid for the K_1 boson interpolation
 */
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK1bosonGrid();
/**
 * @brief initializeK2fermionGrid
 * @return a unique pointer to the grid for the K_2 fermion interpolation
 */
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK2fermionGrid();
/**
 * @brief initializeK4fermionGrid
 * @return a unique pointer to the grid for the K_4 fermion interpolation
 */
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK4fermionGrid();
/**
 * @brief initializeK4bosonGrid
 * @return a unique pointer to the grid for the K_4 boson interpolation
 */
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK4bosonGrid();
/**
 * @brief initializeK5fermionGrid
 * @return a unique pointer to the grid for the K_5 fermion interpolation
 */
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK5fermionGrid();
/**
 * @brief initializeK5bosonGrid
 * @return a unique pointer to the grid for the K_5 boson interpolation
 */
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK5bosonGrid();
/**
 * @brief initializeK6fermionGrid
 * @return a unique pointer to the grid for the K_6 fermion interpolation
 */
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK6fermionGrid();
/**
 * @brief initializeK8fermionGrid
 * @return a unique pointer to the grid for the K_8 fermion interpolation
 */
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK8fermionGrid();
/**
 * @brief initializeK9fermionGrid
 * @return a unique pointer to the grid for the K_9 fermion interpolation
 */
std::unique_ptr<gsl_spline2d,decltype(&gsl_spline2d_free)> initializeK9fermionGrid();
/**
 * Calculates the non normalised function K1 for fermions
 * @param msquared m^2 [GeV^2]
 * @param T temperature [GeV]
 */
double K1fermion(double msquared, double T);

/**
 * Calculates the non normalised function K1 for bosons
 * @param msquared m^2 [GeV^2]
 * @param T temperature [GeV]
 */
double K1boson(double msquared, double T);

/**
 * Calculates the non normalised function K2 for fermions
 * @param msquared m^2 [GeV^2]
 * @param T temperature [GeV]
 */
double K2fermion(double msquared, double T);
/**
 * Calculates the non normalised function K4 for fermions
 * @param msquared m^2 [GeV^2]
 * @param T temperature [GeV]
 */
double K4fermion(double msquared, double T);
/**
 * Calculates the non normalised function K4 for bosons
 * @param msquared m^2 [GeV^2]
 * @param T temperature [GeV]
 */
double K4boson(double msquared, double T);
/**
 * Calculates the non normalised function K5 for fermions
 * @param msquared m^2 [GeV^2]
 * @param T temperature [GeV]
 */
double K5fermion(double msquared, double T);
/**
 * Calculates the non normalised function K5 for bosons
 * @param msquared m^2 [GeV^2]
 * @param T temperature [GeV]
 */
double K5boson(double msquared, double T);
/**
 * Calculates the non normalised function K6 for fermions
 * @param msquared m^2 [GeV^2]
 * @param T temperature [GeV]
 */
double K6fermion(double msquared, double T);
/**
 * Calculates the non normalised function K8 for fermions
 * @param msquared m^2 [GeV^2]
 * @param T temperature [GeV]
 */
double K8fermion(double msquared, double T);
/**
 * Calculates the non normalised function K9 for fermions
 * @param msquared m^2 [GeV^2]
 * @param T temperature [GeV]
 */
double K9fermion(double msquared, double T);


/**
     * Calculates the normalised function K1 for fermions
     * @param msquared m^2 [GeV^2]
     * @param T temperature [GeV]
     */
double K1fermion_normalized(double msquared, double T);
/**
     * Calculates the normalised function K1 for bosons
     * @param msquared m^2 [GeV^2]
     * @param T temperature [GeV]
     */
double K1boson_normalized(double msquared, double T);
/**
     * Calculates the normalised function K2 for fermions
     * @param msquared m^2 [GeV^2]
     * @param T temperature [GeV]
     */
double K2fermion_normalized(double msquared, double T);
/**
     * Calculates the normalised function K4 for fermions
     * @param msquared m^2 [GeV^2]
     * @param T temperature [GeV]
     */
double K4fermion_normalized(double msquared, double T);
/**
     * Calculates the normalised function K4 for bosons
     * @param msquared m^2 [GeV^2]
     * @param T temperature [GeV]
     */
double K4boson_normalized(double msquared, double T);
/**
     * Calculates the normalised function K5 for fermions
     * @param msquared m^2 [GeV^2]
     * @param T temperature [GeV]
     */
double K5fermion_normalized(double msquared, double T);
/**
     * Calculates the normalised function K5 for bosons
     * @param msquared m^2 [GeV^2]
     * @param T temperature [GeV]
     */
double K5boson_normalized(double msquared, double T);
/**
     * Calculates the normalised function K6 for fermions
     * @param msquared m^2 [GeV^2]
     * @param T temperature [GeV]
     */
double K6fermion_normalized(double msquared, double T);
/**
     * Calculates the normalised function K8 for fermions
     * @param msquared m^2 [GeV^2]
     * @param T temperature [GeV]
     */
double K8fermion_normalized(double msquared, double T);
/**
     * Calculates the normalised function K9 for fermions
     * @param msquared m^2 [GeV^2]
     * @param T temperature [GeV]
     */
double K9fermion_normalized(double msquared, double T);



}
}
#endif /* SRC_BARYO_CALCULATION_KFACTORSINTERPOLATEDGSL_H_ */
