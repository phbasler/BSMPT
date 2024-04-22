/*
 * vev_spline.cpp
 *
 * cubic spline interpolation for path in VEV space using constant speed (in
 * particular = 1)
 *
 * We use the lib spline.h to make a spline in each direction, the parameter is
 * the linear length between the points (not important as long as it
 * monotonically increasing) Then we make a thoroughly integration (Simpson 3/8
 * rule) to convert from linear length to spline length, and vice-versa.
 * Interpolate using a spline
 *
 */

#include <BSMPT/utility/Logger.h> // for Logger Class
#include <BSMPT/utility/spline.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#ifndef INCLUDE_BSMPT_GravitationalWaves_GWUtils_ConstantVelocitySpline_H_
#define INCLUDE_BSMPT_GravitationalWaves_GWUtils_ConstantVelocitySpline_H_

namespace cvspline
{
double length_vector(std::vector<double> x0, int dim);

std::vector<std::vector<double>> transpose(std::vector<std::vector<double>> &A);

class cvspline
{
private:
  std::vector<double> list_x; // List of x for linear lengths division
  std::vector<double> list_l; // List of l for spline lengths division
  double lin_abs_deriv(double x);

  double Simpson_step(double t0, double t1);

public:
  // Attributes
  int dim;           // Dimension of the VEV space
  int num_points;    // Number of points in the path (knots)
  int num_inter;     // Number of point from x to l (and viceversa)
  tk::spline x_to_l; // Spline to convert from linear length to spline length
  tk::spline l_to_x; // Spline to convert from spline length to linear length
  float linL;        // Linear length of the (spline) path
  float L;           // True length of the (spline) path
  std::vector<double>
      lin_lengths; // We need a parameter, so I used the linear lengths (not the
                   // path spline length) as parameter, this is not optimal as
                   // the velocity will not be 1.
  std::vector<double>
      vev_position; // Position of the VEV list using spline length
  std::vector<tk::spline>
      splines; // Vector of each 2D splines (x = linear length (not important,
               // must be increasing), y = vev_i)
  std::vector<std::vector<double>>
      transposed_phi; // Transpose of the path given, easier for calculations
  std::vector<std::vector<double>> phipath; // List of VEVs paths

  cvspline(); // Default constructor

  cvspline(std::vector<std::vector<double>> phipath_in);

  cvspline(std::vector<std::vector<double>> phipath_in, int input_num_inter);

  void initialize();
  void add_point(double p, std::vector<double> new_vev, bool compile = true);

  std::vector<double> deriv(int order, double x);

  std::vector<double> dl(double l);

  std::vector<double> d2l(double l);

  std::vector<double> operator()(double l);

  void print_path();

  void save_path(std::string file_name, bool header = true);
};
} // namespace cvspline
#endif