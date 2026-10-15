#pragma once

#include <BSMPT/utility/NumericalIntegration.h>
#include <iostream>
#include <memory>

namespace BSMPT
{
namespace Baryo
{
namespace FHCK
{

/**
 * @brief Types of Kernel functions
 *
 */
enum class KernelType
{
  K,
  D,
  Q,
  Qe,
  Q8o,
  Q9o,
  Rb,
  K4FH
};

/**
 * @brief Types of particles
 *
 */
enum class ParticleType
{
  LeftFermion,
  RightFermion,
  Boson
};

/**
 * @brief Function for integer powers. Better suitable than the std one
 *
 * @param x Value of the number
 * @param i Power to which we take the number
 */
double ipow(const double x, const int n);

/**
 * @brief Fermi/Dirac and Bose/Einstein distribution with their derivatives
 *
 * @param x Value at which we want to evaluate
 * @param s Spin statistic. 1 for fermions, -1 for bosons
 * @param diff Order of derivative {0,1,2}
 */
double f0(const double x, const double s, const int diff);

/**
 * @brief y-dependent integrand of the kernel
 *
 */
class KernelInty
{
private:
  /**
   * @brief Settings which adjust the kernel function F(f^(k), V, n, m)
   * see [2510.21915]. Structure = {0, 1, 2} -> V={1, V_s, V_h}
   *
   */
  const int k, structure, n, m;

  /**
   * @brief Wall velocity, wall gamma factor, spin-statistic, lower boundary of
   * the integral
   *
   */
  const double vw, gamw, s, x;

  /**
   * @brief Variables which occur during the integration over y
   *
   */
  double w, pwt;

public:
  /**
   * @brief Useful prefactor to save computation time
   *
   */
  double pre;

  /**
   * @brief Construct a new Kernel Integrand y object
   *
   * @param k_in derivative of f0
   * @param structure_in structure of V
   * @param n_in n variable of kernel see [2407.13639]
   * @param m_in m variable of kernel see [2407.13639]
   * @param vw_in Wall velocity
   * @param gamw_in Gamma factor of the wall
   * @param s_in Spin statistic of particle
   * @param x_in x = m / T
   */
  KernelInty(const int k_in,
             const int structure_in,
             const int n_in,
             const int m_in,
             const double vw_in,
             const double gamw_in,
             const double s_in,
             const double x_in)
      : k(k_in)
      , structure(structure_in)
      , n(n_in)
      , m(m_in)
      , vw(vw_in)
      , gamw(gamw_in)
      , s(s_in)
      , x(x_in) {};

  /**
   * @brief Sets all relevant variables of the class
   * @param  u Integration variable. Defined through w = x + (1 - u) / u. Takes
   * value between 0 and 1
   */
  void set_all(const double u);

  /**
   * @brief Operator which evaluates the integrand of the y-dependent part of
   * the kernel integrand
   * @param y Integration variable between [-1, 1]
   */
  double operator()(const double y);

  ~KernelInty() {};
};

/**
 * @brief w-dependent integrand of the kernel
 *
 */
class KernelIntw
{
private:
  /**
   * @brief Integrand depending on y.
   */
  KernelInty Integrand;

public:
  /**
   * @brief Construct a new Kernel Integrand w object
   *
   * @param k_in derivative of f0
   * @param structure_in structure of V
   * @param n_in n variable of kernel see [2407.13639]
   * @param m_in m variable of kernel see [2407.13639]
   * @param vw_in Wall velocity
   * @param gamw_in Gamma factor of the wall
   * @param s_in Spin statistic of particle
   * @param x_in x = m / T
   */
  KernelIntw(const int k_in,
             const int structure_in,
             const int n_in,
             const int m_in,
             const double vw_in,
             const double gamw_in,
             const double s_in,
             const double x_in)
      : Integrand(k_in, structure_in, n_in, m_in, vw_in, gamw_in, s_in, x_in) {
      };

  /**
   * @brief Evaluate the integrand depending on u between [0, 1].
   */
  double operator()(const double u);

  ~KernelIntw() {};
};

/**
 * @brief Separate class for Q9 kernels as they involve two normal kernels.
 * y-dependent part
 *
 */
class Q9KernelInty
{
private:
  /**
   * @brief Structure = {0, 1, 2} -> V={1, V_s, V_h}. l is the moment of the
   * kernel
   */
  const int structure, l;

  /**
   * @brief Wall velocity, wall gamma factor, spin-statistic, lower boundary of
   * the integral
   *
   */
  const double vw, gamw, s, x;

  /**
   * @brief Variables which occur during the integration over y
   *
   */
  double w, pwt;

public:
  /**
   * @brief Useful prefactors to save computation time
   *
   */
  double pre1, pre2;

  /**
   * @brief Construct a new Q9 Kernel Integrand y object
   *
   * @param structure_in structure of V
   * @param l_in Moment of the kernel
   * @param vw_in Wall velocity
   * @param gamw_in Gamma factor of the wall
   * @param vw_in Spin statistic of particle
   * @param x_in x = m / T
   */
  Q9KernelInty(const int structure_in,
               const int l_in,
               const double vw_in,
               const double gamw_in,
               const double s_in,
               const double x_in)
      : structure(structure_in)
      , l(l_in)
      , vw(vw_in)
      , gamw(gamw_in)
      , s(s_in)
      , x(x_in) {};

  /**
   * @brief Sets all relevant variables of the class
   * @param  u Integration variable. Defined through w = x + (1 - u) / u. Takes
   * value between 0 and 1
   */
  void set_all(const double u);

  /**
   * @brief operator which evaluates the integrand of the y-dependent part of
   * the kernel integrand
   * @param y Integration variable between [-1, 1]
   */
  double operator()(const double y);

  ~Q9KernelInty() {};
};

/**
 * @brief w-dependent integrand of the Q9 kernel
 *
 */
class Q9KernelIntw
{
private:
  /**
   * @brief Integrand depending on y.
   */
  Q9KernelInty Integrand;

public:
  /**
   * @brief Construct a new Q9 Kernel Integrand w object
   *
   * @param structure_in structure of V
   * @param l_in Moment of the kernel
   * @param vw_in Wall velocity
   * @param gamw_in Gamma factor of the wall
   * @param vw_in Spin statistic of particle
   * @param x_in x = m / T
   */
  Q9KernelIntw(const int structure_in,
               const int l_in,
               const double vw_in,
               const double gamw_in,
               const double s_in,
               const double x_in)
      : Integrand(structure_in, l_in, vw_in, gamw_in, s_in, x_in) {};

  /**
   * @brief Evaluate the integrand depending on u between [0, 1].
   */
  double operator()(const double u);
  ~Q9KernelIntw() {};
};

/**
 * @brief Integrand of the N0 kernel
 */
class N0Int
{
private:
  /**
   * @brief Wall velocity, wall gamma factor, spin-statistic, lower boundary of
   * the integral
   *
   */
  const double s, x;

public:
  /**
   * @brief Construct a new N0 Kernel Integrand object
   *
   * @param vw_in Spin statistic of particle
   * @param x_in x = m / T
   */
  N0Int(const double s_in, const double x_in) : s(s_in), x(x_in) {}

  /**
   * @brief Evaluate the integrand depending on u between [0, 1].
   */
  double operator()(const double u);
  ~N0Int() {};
};

/**
 * @brief Integrand of the Rbar kernel
 */
class RbarInt
{
private:
  /**
   * @brief Wall velocity, wall gamma factor, spin-statistic, lower boundary of
   * the integral
   *
   */
  const double vw, gamw, s, x;

public:
  /**
   * @brief Construct a new Rbar Kernel Integrand object
   *
   * @param vw_in Wall velocity
   * @param gamw_in Gamma factor of the wall
   * @param vw_in Spin statistic of particle
   * @param x_in x = m / T
   */
  RbarInt(const double vw_in,
          const double gamw_in,
          const double s_in,
          const double x_in)
      : vw(vw_in)
      , gamw(gamw_in)
      , s(s_in)
      , x(x_in) {};

  /**
   * @brief Evaluate the integrand depending on u between [0, 1].
   */
  double operator()(const double u);
  ~RbarInt() {};
};

/**
 * @brief Integrand of the K4 kernel
 */
class K4Int
{
private:
  /**
   * @brief spin-statistic, lower boundary of the integral
   *
   */
  const double s, x;

public:
  /**
   * @brief Construct a new K4 Kernel Integrand object
   *
   * @param vw_in Spin statistic of particle
   * @param x_in x = m / T
   */
  K4Int(const double s_in, const double x_in) : s(s_in), x(x_in) {}

  /**
   * @brief Evaluate the integrand depending on u between [0, 1].
   */
  double operator()(const double u);

  ~K4Int() {};
};

/**
 * @brief Generic Kernel function for any moment
 */
class Kernel
{
private:
  /**
   * @brief Moment of the Kernel
   */
  const int l;

  /**
   * @brief Structure = {0, 1, 2} -> V={1, V_s, V_h}
   */
  const int structure;

public:
  /**
   * @brief Construct a new Kernel object
   *
   * @param l_in l-th moment of the kernel
   * @param structure_in structure of V
   */
  Kernel(const int l_in, const int structure_in)
      : l(l_in)
      , structure(structure_in) {};

  /**
   * @brief Evaluates the needed kernel at a given x and wall velocity
   * @param KernelType Type of Kernel
   * @param ParticleType Type of Particle
   * @param x Value at which we want to evaluate the kernel x=m/T
   * @param vw Wall velocity at which we want to evaluate the kernel
   */
  double operator()(const KernelType Kern,
                    const ParticleType P,
                    const double x,
                    const double vw);
  ~Kernel() {};
};

} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT