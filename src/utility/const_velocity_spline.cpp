// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/*
 * spline.cpp
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

#include <BSMPT/utility/const_velocity_spline.h>

double cvspline::lin_abs_deriv(
    double x) // Length of vector speed of non constant velocity spline
{
  double r = 0;
  for (int i = 0; i < dim; i++)
  {
    r += std::pow(splines[i].deriv(1, x), 2); // Calculates the speed vector
  }
  return std::sqrt(r);
}
double cvspline::Simpson_step(
    double t0,
    double
        t1) // Function to take single integration step using Simpson 3/8 rule
            // https://en.wikipedia.org/wiki/Simpson%27s_rule#:~:text=Simpson%27s%20method.-,Simpson%27s%203/8%20rule%5Bedit%5D,-Simpson%27s%203/8
{
  double r = lin_abs_deriv(t0);
  r += 3 * lin_abs_deriv((2 * t0 + t1) / 3);
  r += 3 * lin_abs_deriv((t0 + 2 * t1) / 3);
  r += lin_abs_deriv(t1);

  return r * (t1 - t0) / 8;

  /* 1/3 Simpson Rule
  double r = lin_abs_deriv(t0);
  r += 4*lin_abs_deriv((t0+t1)/2);
  r += lin_abs_deriv(t1);
  r *= (t1-t0)/6;
  return r;
  */
}

cvspline::cvspline() = default; // Default constructor

cvspline::cvspline(const std::vector<std::vector<double>> &phipath_in)
    : num_inter(10000) // Constructor
{
  this->phipath = phipath_in;
  initialize();
}

cvspline::cvspline(
    const std::vector<std::vector<double>> &phipath_in,
    int input_num_inter) // Constructor with specified num of interpolations
{
  this->phipath = phipath_in;
  num_inter     = input_num_inter;
  initialize();
}

void cvspline::initialize() // Initalizes the class
{
  dim        = phipath[0].size();
  num_points = phipath.size();

  if (phipath.size() == 2) // Cannot fit a spline with two points (or less),
                           // adding a middle point if len = 2
  {
    std::vector<double> middle_vector1;
    std::vector<double> middle_vector2;

    for (int i = 0; i < dim; i++)
    {
      middle_vector1.push_back(phipath[0][i] * 2 / 3 +
                               phipath[1][i] * 1 /
                                   3); // Calculates middle point
    }
    // phipath.insert(phipath.begin() + 1, middle_vector); // Adds it to the
    // 2nd position

    // middle_vector.clear();

    for (int i = 0; i < dim; i++)
    {
      middle_vector2.push_back(phipath[0][i] * 1 / 3 +
                               phipath[1][i] * 2 /
                                   3); // Calculates middle point
    }
    phipath.insert(phipath.begin() + 1,
                   middle_vector1); // Adds it to the 2nd position
    phipath.insert(phipath.begin() + 2,
                   middle_vector2); // Adds it to the 3nd position

    num_points = phipath.size();
  }

  std::vector<double> direction;
  lin_lengths = {0}; // Reset the list
  for (std::size_t i = 1; i < phipath.size(); i++)
  {
    direction = phipath[i];
    for (int j = 0; j < dim; j++)
    {
      direction[j] -= phipath[i - 1][j];
    }
    lin_lengths.push_back(lin_lengths.back() + BSMPT::L2NormVector(direction));
  }

  linL =
      lin_lengths
          .back(); // Linear Length is the last element of the list of lengths

  transposed_phi = BSMPT::Transpose(
      phipath); // Transposes the path, helpfull when computing the spline
  splines = {}; // Reset splines list
  for (int i = 0; i < dim; i++)
  {
    tk::spline s(lin_lengths,
                 transposed_phi[i],
                 tk::spline::cspline,
                 false,
                 tk::spline::not_a_knot,
                 0.0,
                 tk::spline::not_a_knot,
                 0.0);
    splines.push_back(s);
  }

  list_l = {0}; // Reset the list
  list_x = {0}; // Reset the list

  for (int i = 1; i <= num_inter;
       i++) // Integrates the path using Simpson 3/8 rule to find the
            // correspondence between x and l
  {
    double step = linL / num_inter;
    list_l.push_back(list_l.back() +
                     Simpson_step(list_l.back(), list_l.back() + step));
    list_x.push_back(list_x.back() + step);
  }
  L = list_l.back(); // Set complet length of the path
                     // This method seems to be stable, l_to_x(x_to_l(value) =
                     // value
  x_to_l.set_points(list_x, list_l); // Compute x to l spline
  l_to_x.set_points(list_l, list_x); // Compute l to x spline

  vev_position = {}; // Reset the list
  for (int i = 0; i < num_points;
       i++) // Integrates the path using Simpson 3/8 rule to find the
            // correspondence between x and l
  {
    vev_position.push_back(x_to_l(lin_lengths[i]));
  }
}
void cvspline::add_point(const double &p,
                         const std::vector<double> &new_vev,
                         const bool &compile)
{
  // Function to add a new vev "new_vev" at the position "p"
  // We allow for p < 0 as an easy way to add a new true vacuum
  // The position "p" is going to be rewritten when the path is recompiled
  // It is necessary in order to disambiguate the order of the VEVs

  if (p < 0)
  {
    phipath.insert(phipath.begin(), new_vev); // Add VEV
    if (compile) initialize();                // Redo calculations
    return;
  }
  // We allow for p > L as an easy way to add a new false vacuum, no (good)
  // reason to do this tho
  if (p > L)
  {
    phipath.push_back(new_vev); // Add VEV
    if (compile) initialize();  // Redo calculations
    return;
  }

  int lower = 0;              // Index of lower bound
  int upper = phipath.size(); // Index of upper bound
  int lower_temp, upper_temp; // Floor/Ceil of the middle point

  while (upper - lower !=
         1) // When upper - lower = 1 then we add the point at that position
  {

    // Use the integer version of the bisection method
    lower_temp = int(floor((upper + lower) / 2)); // Floor of the middle point
    upper_temp = int(ceil((upper + lower) / 2));  // Ceiling of the middle point

    if ((vev_position[upper] == upper_temp) or
        (vev_position[lower] == lower_temp))
    {
      // Point already exists in the position "p", no point is added
      BSMPT::Logger::Write(BSMPT::LoggingLevel::TransitionDetailed,
                           "Point in that position already exists!");
      return;
    }
    if (vev_position[lower_temp] < p)
    {
      lower = lower_temp;
    }
    if (vev_position[upper_temp] > p)
    {
      upper = upper_temp;
    }
  }
  phipath.insert(phipath.begin() + upper, new_vev); // Add the point
  vev_position.insert(vev_position.begin() + upper,
                      p); // Update the "vev_positions" list, very importante
                          // to initialize() after all new VEVs are added.

  if (compile) initialize(); // Redo calculations
}
std::vector<double> cvspline::deriv(int order, double x)
{
  // Derivatives of spline with ( linear argument - speed with norm NOT equal
  // to 1)
  std::vector<double> r;
  for (int i = 0; i < dim; i++)
  {
    r.push_back(splines[i].deriv(order, x)); // Calculates the speed vector
  }
  return r;
}

std::vector<double> cvspline::dl(double l)
{
  // Vector speed of spline with size = 1 (this way is analytically)
  std::vector<double> r;
  double tL = 0;
  for (int i = 0; i < dim; i++)
  {
    r.push_back(splines[i].deriv(1, l_to_x(l))); // Calculates the speed vector
    tL += r.back() * r.back(); // Calculates the length of the speed vector
  }
  tL = std::sqrt(tL);
  for (int i = 0; i < dim; i++)
  {
    r[i] = r[i] / tL; // Normalizes speed vector
  }

  return r;
}
std::vector<double> cvspline::d2l(double l)
{
  // Normal acceleration (since v = 1 all acceleration is normal) of spline
  // The analytical form seems less stable than the numerical one,
  // i.e. numerical(| gamma'(s).gamma''(s) |) < | gamma'(s).gamma''(s) | which
  // should be zero We're still going to use the analytical version in the
  // meanwhile

  std::vector<double> rr, gammaPPx, gammaPl;
  double dldxSqrd, d2ld2x;

  double x = l_to_x(l);

  gammaPPx = deriv(2, x);                     // Gamma''(x)
  gammaPl  = dl(l);                           // Gamma'(l)
  dldxSqrd = std::pow(x_to_l.deriv(1, x), 2); // dldx(l)^2 = l'(x)^2
  d2ld2x   = x_to_l.deriv(2, x);              // d2ldx(l) = l''(x)

  for (int i = 0; i < dim; i++)
  {
    rr.push_back((gammaPPx[i] - gammaPl[i] * d2ld2x) / dldxSqrd);
  }

  return rr;
}
std::vector<double> cvspline::operator()(double l)
{
  // Calculates the VEV at the position l in spline length,
  l = l_to_x(l); // Convert from spline length to linear length.
  std::vector<double> r;
  for (int i = 0; i < dim; i++)
  {
    r.push_back(splines[i](l));
  }
  return r;
}

void cvspline::print_path()
{
  int wid = 15;

  std::stringstream ss;

  ss << std::string(15 * (dim + 2), '-') << std::endl;
  ss << std::setw(5) << std::right << "id" << std::setw(wid) << std::right
     << "l";

  for (int i = 0; i < dim; i++)
  {
    ss << std::setw(wid) << std::right << ("Field " + std::to_string(i))
       << std::setw(wid) << std::right;
  }

  ss << "\n" << std::string(15 * (dim + 2), '-') << std::endl;

  // ss << "-------------------------------\n";
  for (int i = 0; i < num_points; i++)
  {

    ss << std::setw(5) << std::right << i << std::setw(wid) << std::right
       << vev_position[i] << std::setw(wid) << std::right;
    for (int j = 0; j < dim; j++)
    {
      ss << phipath[i][j] << std::setw(wid) << std::right;
    }
    ss << " " << std::endl;
  }
  ss << std::string(15 * (dim + 2), '-') << std::endl;

  BSMPT::Logger::Write(BSMPT::LoggingLevel::BounceDetailed, ss.str());
}
void cvspline::save_path(const std::string &file_name, const bool &header)
{
  std::ofstream myfile;
  myfile.open(file_name);
  if (header)
  {
    for (int i = 0; i < dim; i++)
    {
      myfile << "field " << i << "\t";
    }
    myfile << "\n";
  }
  for (int i = 0; i < num_points; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      myfile << phipath[i][j] << "\t";
    }
    myfile << " "
           << "\n";
  }

  myfile.close();
}
