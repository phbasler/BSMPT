#include <BSMPT/utility/NumericalDerivatives.h>

namespace BSMPT
{

std::vector<double>
NablaNumerical(const std::vector<double> &phi,
               const std::function<double(std::vector<double>)> &f,
               const double &eps,
               const int &dim)
{
  std::vector<double> result(dim);

  for (int i = 0; i < dim; i++)
  {
    std::vector<double> lp2 = phi;
    lp2[i] += 2 * eps;
    std::vector<double> lp1 = phi;
    lp1[i] += eps;
    std::vector<double> lm1 = phi;
    lm1[i] -= eps;
    std::vector<double> lm2 = phi;
    lm2[i] -= 2 * eps;
    result[i] = (-f(lp2) + 8 * f(lp1) - 8 * f(lm1) + f(lm2)) / (12 * eps);
  }
  return result;
}


std::vector<std::vector<double>> HessianNumerical_BI(
    const std::vector<double> &phi,
    const std::function<double(std::vector<double>)> &V,
    double eps,
    const int &dim)
{
  std::vector<std::vector<double>> result(dim, std::vector<double>(dim));
  for (int i = 0; i < dim; i++)
  {
    // https://en.wikipedia.org/wiki/Finite_difference
    for (int j = i; j < dim; j++)
    {
      double r = 0;

      if (i == j) eps /= 2;

      std::vector<double> xp = phi; // F(x+h, y+h)
      xp[i] += eps;
      xp[j] += eps;
      r += V(xp);

      xp = phi; //-F(x+h, y-h)
      xp[i] += eps;
      xp[j] -= eps;
      r -= V(xp);

      xp = phi; //-F(x-h, y+h)
      xp[i] -= eps;
      xp[j] += eps;
      r -= V(xp);

      xp = phi; // F(x-h, y-h)
      xp[i] -= eps;
      xp[j] -= eps;
      r += V(xp);

      result[i][j] = r / (4 * eps * eps);
      result[j][i] = r / (4 * eps * eps);

      if (i == j) eps *= 2;
    }
  }
  return result;
}

std::vector<std::vector<double>>
HessianNumerical_MT(
    const std::vector<double> &phi,
                 const std::function<double(std::vector<double>)> &f,
                 const double &eps,
                 const int &dim)
{
  std::vector<std::vector<double>> result(dim, std::vector<double>(dim));
  for (int i = 0; i < dim; i++)
  {
    // https://en.wikipedia.org/wiki/Finite_difference
    for (int j = 0; j < dim; j++)
    {
      double r = 0;

      std::vector<double> xp = phi; // F(x+h, y+h)
      xp[i] += eps;
      xp[j] += eps;
      r += f(xp);

      xp = phi; //-F(x+h, y-h)
      xp[i] += eps;
      xp[j] -= eps;
      r -= f(xp);

      xp = phi; //-F(x-h, y+h)
      xp[i] -= eps;
      xp[j] += eps;
      r -= f(xp);

      xp = phi; // F(x-h, y-h)
      xp[i] -= eps;
      xp[j] -= eps;
      r += f(xp);

      result[i][j] = r / (4 * eps * eps);
    }
  }
  return result;
}
}