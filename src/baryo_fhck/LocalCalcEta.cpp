#include <BSMPT/baryo_fhck/LocalCalcEta.h>
#include <BSMPT/minimum_tracer/minimum_tracer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/NumericalDerivatives.h>
#include <BSMPT/utility/utility.h>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

namespace BSMPT
{
namespace Baryo
{
namespace FHCK
{

std::vector<double> LocalCalcEta::FromPlane(std::vector<double> vev,
                                            const double &d)
{
  const double missing_element = (d - reducedn * vev) / n[index];
  vev.insert(vev.begin() + index, missing_element);
  return vev;
}

std::vector<double> LocalCalcEta::Reduce(std::vector<double> vev)
{
  vev.erase(vev.begin() + index);
  return vev;
}

void LocalCalcEta::GradProjection()
{
  // Generate identity matrix
  gradprojection =
      std::vector<std::vector<double>>(dim, std::vector<double>(dim, 0.0));
  for (size_t i = 0; i < dim; i++)
  {
    gradprojection.at(i).at(i)     = 1.;
    gradprojection.at(i).at(index) = -n.at(i) / n.at(index); // project
  }
  gradprojection.erase(gradprojection.begin() + index);
}

LocalCalcEta::LocalCalcEta(
    const std::shared_ptr<Class_Potential_Origin> &modelPointer_input,
    const double &T_In,
    const std::vector<double> &TrueVacuum_In,
    const std::vector<double> &FalseVacuum_In)
    : modelPointer(modelPointer_input)
    , dim(TrueVacuum_In.size())
    , T(T_In)
    , TrueVacuum(TrueVacuum_In)
    , FalseVacuum(FalseVacuum_In)
{

  std::stringstream ss;
  ss << "\n-------------- Eta Interface (Local Minimizer) --------------\n";

  n        = FalseVacuum - TrueVacuum;
  index    = std::distance(n.begin(),
                        std::max_element(n.begin(),
                                         n.end(),
                                         [](const int &a, const int &b)
                                         { return abs(a) < abs(b); }));
  n        = n / L2NormVector(n);
  reducedn = Reduce(n);
  ss << "Normal vector\t" << n << "\n";
  ss << "Index\t" << index << "\n";

  GradProjection();
  ss << "\nGrad projection matrix\n";
  for (const auto &i : gradprojection)
  {
    for (const auto &j : i)
      ss << j << "\t";
    ss << "\n";
  }
  Logger::Write(LoggingLevel::FHCK, ss.str());

  CalculateLw();

  CalculateTheta();
}

void LocalCalcEta::CalculateLw()
{
  double eps                                   = 0.1;
  double GradientThreshold                     = 1e-3;
  std::function<double(std::vector<double>)> V = [&](std::vector<double> vev)
  {
    // Potential wrapper
    std::vector<double> res = modelPointer->MinimizeOrderVEV(vev);
    return modelPointer->VEff(res, T);
  };
  std::function<std::vector<double>(std::vector<double>)> dV =
      [=](auto const &arg) { return NablaNumerical(arg, V, eps); };
  std::function<std::vector<std::vector<double>>(std::vector<double>)> Hessian =
      [=](auto const &arg) { return HessianNumerical(arg, V, eps); };

  // Start at the false vacuum
  std::vector<double> new_point, point = Reduce(FalseVacuum);
  MinimumTracer mintracer;

  size_t num_points = 200;

  double Vb          = -1;
  const double Vtrue = V(TrueVacuum);
  const double DistanceThreshold =
      L2NormVector(TrueVacuum - FalseVacuum) / num_points * 10.;

  for (size_t b = 0; b <= num_points; b++)
  {
    std::vector<double> basePoint =
        static_cast<double>(b) / num_points * TrueVacuum +
        (1. - static_cast<double>(b) / num_points) * FalseVacuum;
    const double d = basePoint * n;
    // Calculates two iteration when StraightLineApproximation
    if ((not StraightLineApproximation) or (b < 2))
      new_point = mintracer.LocateMinimum(
          point,
          [this, d, dV](auto const &vev)
          { return gradprojection * dV(FromPlane(vev, d)); },
          [this, d, Hessian](auto const &vev)
          {
            return gradprojection * Hessian(FromPlane(vev, d)) *
                   Transpose(gradprojection);
          },
          1e-2 * GradientThreshold * std::min(1., dim - 1.));

    if (L2NormVector(new_point - point) > DistanceThreshold)
    {
      // Diverged
      return;
    }

    path.push_back(FromPlane(new_point, d));

    if (StraightLineApproximation)
      Vb = std::max(Vb, V(basePoint) - Vtrue);
    else
      Vb = std::max(Vb, V(path.back()) - Vtrue);

    point = new_point;
  }

  if (StraightLineApproximation and (Vb > 0))
    Lw = L2NormVector(TrueVacuum - FalseVacuum) / sqrt(8 * Vb);

  if ((not StraightLineApproximation) and (Vb > 0) and
      (L2NormVector(TrueVacuum - path.back()) < dim * 1e-1) /* extra check */)
    Lw = L2NormVector(TrueVacuum - FalseVacuum) / sqrt(8 * Vb);

  std::stringstream ss;

  ss << "L2\t" << L2NormVector(TrueVacuum - FalseVacuum) << "\n";
  ss << "Vb\t" << Vb << "\n";
  ss << "Distance final point path to TrueVacuum \t"
     << L2NormVector(TrueVacuum - path.back()) << "\n";

  ss << "Lw =\t" << Lw.value_or(NAN);
  Logger::Write(LoggingLevel::FHCK, ss.str());
}

void LocalCalcEta::CalculateTheta()
{
  if (path.size() < 1 and not Lw.has_value())
  {
    Logger::Write(LoggingLevel::FHCK, "Masses phase could not be calculated");
    return;
  }

  Eigen::MatrixXcd MIJQuarks =
      modelPointer->QuarkMassMatrix(modelPointer->MinimizeOrderVEV(path.at(1)));
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> esQuark(MIJQuarks);

  double fac;

  // symmetruc

  std::complex<double> topmass =
      esQuark.eigenvalues()[esQuark.eigenvalues().size() - 1];
  fac = 2 * (topmass.real() >= 0) - 1;

  top_theta_sym = std::arg(topmass / fac);

  std::complex<double> botmass =
      esQuark.eigenvalues()[esQuark.eigenvalues().size() - 1 - 2];
  fac           = 2 * (botmass.real() >= 0) - 1;
  bot_theta_sym = std::arg(botmass / fac);

  std::stringstream ss;

  ss << "theta_sym(mtop) = " << top_theta_sym.value() << "\n";
  ss << "theta_sym(mbot) = " << bot_theta_sym.value() << "\n";

  // brk

  MIJQuarks =
      modelPointer->QuarkMassMatrix(modelPointer->MinimizeOrderVEV(TrueVacuum));
  esQuark = Eigen::ComplexEigenSolver<Eigen::MatrixXcd>(MIJQuarks);

  // symmetric

  topmass = esQuark.eigenvalues()[esQuark.eigenvalues().size() - 1];
  fac     = 2 * (topmass.real() >= 0) - 1;

  top_theta_brk = std::arg(topmass / fac);

  botmass       = esQuark.eigenvalues()[esQuark.eigenvalues().size() - 1 - 2];
  fac           = 2 * (botmass.real() >= 0) - 1;
  bot_theta_brk = std::arg(botmass / fac);

  ss << "theta_brk(mtop) = " << top_theta_brk.value() << "\n";
  ss << "theta_brk(mbot) = " << bot_theta_brk.value() << "\n";

  Logger::Write(LoggingLevel::FHCK, ss.str());
}

} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT