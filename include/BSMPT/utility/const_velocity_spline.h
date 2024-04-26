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

#pragma once

#include <BSMPT/utility/Logger.h> // for Logger Class
#include <BSMPT/utility/spline.h>
#include <BSMPT/utility/utility.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

/**
 * @brief Constructs a spline \f$ s(l) \f$ with constant velocity, i.e.  \f$
 * \frac{ds(l)}{dl} \equiv 1\f$. \f$ l \f$  acts as the length alongisde the
 * spline. It works by using tk::spline to construct a cubic spline which
 * depends on another arbitrary parameter \f$ x \f$ and then construct an
 * additional spline to convert from \f$ x \to l \equiv l(x) \f$
 *
 */
class cvspline
{
private:
  /**
   * @brief List of x for linear lengths division
   *
   */
  std::vector<double> list_x;
  /**
   * @brief // List of l for spline lengths division
   *
   */
  std::vector<double> list_l;

  /**
   * @brief derivative in \f$ x \f$
   *
   * @param x
   * @return double
   */
  double lin_abs_deriv(double x);

  /**
   * @brief Integrates from \f$ t_0 \f$ to \f$ t_1 \f$ using a single Simposons'
   * 3/8 integration step.
   *
   * @param t0 lower limit
   * @param t1 upper limit
   * @return double value of the integral
   */
  double Simpson_step(double t0, double t1);

public:
  /**
   * @brief Dimension of the VEV space
   *
   */
  int dim;
  /**
   * @brief Number of points in the path (knots)
   *
   */
  int num_points;
  /**
   * @brief Number of point from x to l (and viceversa)
   *
   */
  int num_inter;
  /**
   * @brief Spline to convert from linear length to spline length
   *
   */
  tk::spline x_to_l;
  /**
   * @brief Spline to convert from spline length to linear length
   *
   */
  tk::spline l_to_x;
  /**
   * @brief Linear length of the (spline) path
   *
   */
  float linL;
  /**
   * @brief True length of the (spline) path
   *
   */
  float L;
  /**
   * @brief We need a parameter, so I used the linear lengths (not the path
   * spline length) as parameter, this is not optimal as the velocity will not
   * be 1.
   *
   */
  std::vector<double> lin_lengths;
  /**
   * @brief Position of the VEV list "lin_lengths" using spline length
   *
   */
  std::vector<double> vev_position;
  /**
   * @brief Vector of each 2D splines (x = linear length (not important,
   * must be increasing), y = vev_i)
   *
   */
  std::vector<tk::spline> splines;
  /**
   * @brief Transpose of the path given, easier for calculations
   *
   */
  std::vector<std::vector<double>> transposed_phi;
  /**
   * @brief // List of VEVs paths
   *
   */
  std::vector<std::vector<double>> phipath;
  /**
   * @brief Default constructor
   *
   */
  cvspline();
  /**
   * @brief Construct a new cvspline object
   *
   * @param phipath_in knots of the path
   */
  cvspline(std::vector<std::vector<double>> phipath_in);
  /**
   * @brief Construct a new cvspline object
   *
   * @param phipath_in knots of the path
   * @param input_num_inter number of interpolations between \f$ x \f$ and \f$ l
   * \f$
   */
  cvspline(std::vector<std::vector<double>> phipath_in, int input_num_inter);
  /**
   * @brief initialize the constant velocity spline using the given path
   *
   */
  void initialize();
  /**
   * @brief introduce a new point on the knot list
   *
   * @param p parameter of the point
   * @param new_vev value of the knot
   * @param compile initialize after point is added?
   */
  void add_point(double p, std::vector<double> new_vev, bool compile = true);
  /**
   * @brief derivative of the constant velocity spline in \f$ x \f$
   *
   * @param order order of the derivative
   * @param x where to calculate the derivative
   * @return std::vector<double> result
   */
  std::vector<double> deriv(int order, double x);
  /**
   * @brief first derivative of the constant velocity spline in \f$ l \f$
   *
   * @param l where to calculate the derivative
   * @return std::vector<double> result
   */
  std::vector<double> dl(double l);
  /**
   * @brief second derivative of the constant velocity spline in \f$ l \f$
   *
   * @param l where to calculate the derivative
   * @return std::vector<double> result
   */
  std::vector<double> d2l(double l);
  /**
   * @brief value of the constant velocity spline at \f$ l \f$
   *
   * @param l
   * @return std::vector<double>
   */
  std::vector<double> operator()(double l);
  /**
   * @brief print the current knots of the spline
   *
   */
  void print_path();
  /**
   * @brief save the knots of the splien into a file
   *
   * @param file_name name of the file
   * @param header print the names of the VEVs in the first column?
   */
  void save_path(std::string file_name, bool header = true);
};