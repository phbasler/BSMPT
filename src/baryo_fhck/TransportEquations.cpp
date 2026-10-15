#include <BSMPT/baryo_fhck/TransportEquations.h>
#include <BSMPT/utility/asciiplotter/asciiplotter.h>
#include <fstream>
#include <sstream>
#define _USE_MATH_DEFINES
#include <cmath>
namespace BSMPT
{

namespace Baryo
{
namespace FHCK
{
TransportEquations::TransportEquations(
    const std::shared_ptr<TransportModel> &model_in,
    const double &Tstar_in,
    const std::vector<size_t> &moments_in)
    : transportmodel(model_in)
    , Tstar(Tstar_in)
    , moments(moments_in)
    , BAUeta(moments.size())
{
  std::stringstream ss;
  ss << "Moments to calculate = " << moments;
  Logger::Write(LoggingLevel::FHCK, ss.str());
  Initialize();
}

void TransportEquations::Initialize()
{
  transportmodel->Initialize();

  if (transportmodel->status != TransportModelStatus::Success)
  {
    Logger::Write(LoggingLevel::FHCK,
                  "Field profile calculation failed. Aborting...\n");
    return;
  }

  // Plot top mass
  NumberOfSteps = DefaultNumberOfSteps;
  LwMultiplier  = DefaultLwMultiplier;
  MakeDistribution(LwMultiplier * transportmodel->Lw, NumberOfSteps);
  transportmodel->GenerateFermionMass(zList, true);

  gamwall = 1. / std::sqrt(1. - transportmodel->vwall * transportmodel->vwall);

  Logger::Write(LoggingLevel::FHCK, "Building Kernels Interpolations...");

  BuildKernelInterpolation();

  Logger::Write(LoggingLevel::FHCK, "\033[92mSuccess.\033[0m");
}

void TransportEquations::GenerateIntegrationSpace()
{
  // These are only used to calculate the estimated NumberOfSteps and
  // LwMultiplier
  NumberOfSteps = DefaultNumberOfSteps;
  LwMultiplier  = DefaultLwMultiplier;

  MakeDistribution(LwMultiplier * transportmodel->Lw, NumberOfSteps);
  transportmodel->GenerateFermionMass(zList);
  Solution = MatDoub(nEqs, zList.size(), 0.);

  double HighestNegReEigenvalue, HighestImEigenvalue;
  CheckBoundary(HighestNegReEigenvalue, HighestImEigenvalue);

  stringstream ss;

  ss << "HighestNegRe  = " << HighestNegReEigenvalue << "\n";
  ss << "HighestNegEigenvalue  = " << HighestImEigenvalue;

  Logger::Write(LoggingLevel::FHCK, ss.str());

  LwMultiplier = std::ceil(std::log(LwMultiplierCutoff) /
                           (transportmodel->Lw * HighestNegReEigenvalue));

  NumberOfSteps =
      StepsPerCycle * std::ceil(2 * LwMultiplier * transportmodel->Lw *
                                HighestImEigenvalue / (2 * M_PI));

  // Upper bound on NumberOfSteps
  const double LowHighStepFactor =
      StepsPerCycle /
      StepsPerCycleHigh; // prevent same MaxNumberOfSteps for high and low cycle
  NumberOfSteps = std::min(
      NumberOfSteps, (size_t)std::ceil(MaxNumberOfSteps * LowHighStepFactor));
  // Lower bound on LwMultiplier
  LwMultiplier = std::max(LwMultiplier, MinLwMultiplier);

  Logger::Write(LoggingLevel::FHCK,
                "Calculated NumberOfSteps = " + std::to_string(NumberOfSteps));
  Logger::Write(LoggingLevel::FHCK,
                "Calculated LwMultiplier = " + std::to_string(LwMultiplier) +
                    "\n");

  MakeDistribution(LwMultiplier * transportmodel->Lw, NumberOfSteps);
  transportmodel->GenerateFermionMass(zList);
  Solution = MatDoub(nEqs, zList.size(), 0.);

  Logger::Write(LoggingLevel::FHCK,
                "Integration limits in z | " + std::to_string(zList.front()) +
                    " -> " + std::to_string(zList.back()));
}

void TransportEquations::InitializeMoment(const size_t &moment_in)
{
  moment = moment_in;

  nParticles = nFermions + nBosons;
  nEqs       = moment * nParticles;

  MtildeM.resize(nEqs, nEqs);
  MtildeP.resize(nEqs, nEqs);
  StildeM.resize(nEqs);
  StildeP.resize(nEqs);

  GenerateIntegrationSpace();
}

tk::spline TransportEquations::InterpolateKernel(const std::string &kernel_name,
                                                 const bool is_1D)
{
  double x, x_save, vw, y;
  std::vector<double> vw_vals, y_vals, x_vals, res_vals;

  std::string abscissaeFile = std::string(KERNEL_DIR) + "/x_vw.dat";
  std::string kernelFile = std::string(KERNEL_DIR) + "/" + kernel_name + ".dat";
  std::ifstream kernel(kernelFile);
  std::ifstream abscissae(abscissaeFile);
  if (is_1D)
  {
    while (kernel >> x >> y)
    {
      x_vals.push_back(x);
      res_vals.push_back(y);
    }
  }
  else
  {
    abscissae >> x >> vw;
    kernel >> y;
    x_save = x;
    vw_vals.push_back(vw);
    y_vals.push_back(y);

    while (abscissae >> x >> vw && kernel >> y)
    {
      if (x != x_save)
      {
        tk::spline spl(vw_vals, y_vals);
        x_vals.push_back(x_save);
        res_vals.push_back(spl(transportmodel->vwall));
        vw_vals.clear();
        y_vals.clear();
      }
      vw_vals.push_back(vw);
      y_vals.push_back(y);
      x_save = x;
    }
    tk::spline spl(vw_vals, y_vals);
    x_vals.push_back(x_save);
    res_vals.push_back(spl(transportmodel->vwall));
  }
  kernel.close();
  abscissae.close();

  return tk::spline(x_vals, res_vals);
}

void TransportEquations::BuildKernelInterpolation()
{
  std::string kernel_name;
  bool is_1d = true;

  K4FHf = InterpolateKernel("K4FH_f", true);
  K4FHb = InterpolateKernel("K4FH_b", true);
  Rbarf = InterpolateKernel("Rbar_f", false);
  Rbarb = InterpolateKernel("Rbar_b", false);

  // Pad first element with empty spline
  Qlf.push_back(tk::spline());
  Q8ol.push_back(tk::spline());
  Q9ol.push_back(tk::spline());
  Qlb.push_back(tk::spline());

  const size_t max_moment = *std::max_element(moments.begin(), moments.end());

  for (size_t l = 0; l <= max_moment; l++)
  {
    for (int type = 0; type <= 1; type++)
    {
      std::string suffix = (type == 0 ? "_f" : "_b");

      kernel_name = "D" + std::to_string(l) + suffix;
      type == 0 ? Dlf.push_back(InterpolateKernel(kernel_name, is_1d))
                : Dlb.push_back(InterpolateKernel(kernel_name, is_1d));

      if (l == 1)
        type == 0 ? Klf.push_back(tk::spline()) : Klb.push_back(tk::spline());
      else
      {
        kernel_name = "K" + std::to_string(l) + suffix;
        type == 0 ? Klf.push_back(InterpolateKernel(kernel_name, is_1d))
                  : Klb.push_back(InterpolateKernel(kernel_name, is_1d));
      }

      if (l != 0)
      {
        kernel_name = "Q" + std::to_string(l) + suffix;
        type == 0 ? Qlf.push_back(InterpolateKernel(kernel_name, is_1d))
                  : Qlb.push_back(InterpolateKernel(kernel_name, is_1d));

        if (type == 0)
        {
          kernel_name = "Q8o" + std::to_string(l) + suffix;
          Q8ol.push_back(InterpolateKernel(kernel_name, is_1d));

          kernel_name = "Q9o" + std::to_string(l) + suffix;
          Q9ol.push_back(InterpolateKernel(kernel_name, is_1d));
        }
      }
    }
    is_1d = false;
  }
}

MatDoub TransportEquations::CalculateCollisionMatrix(const double &xW,
                                                     VecDoub &FermionMasses,
                                                     VecDoub &BosonMasses)
{
  constexpr double KAPPABE = 0.685;
  constexpr double KAPPAFD = 0.912;
  if (FermionMasses.size() != 3)
    throw std::runtime_error(
        "Transport equation are only supported with 3 Weyl fermions.");

  if (BosonMasses.size() != 1)
    throw std::runtime_error(
        "Transport equation are only supported with 1 boson.");

  std::vector<double> mass = {sqrt(FermionMasses[0]),
                              sqrt(FermionMasses[1]),
                              sqrt(FermionMasses[2]),
                              sqrt(BosonMasses[0])};

  std::vector<ParticleType> ptype = {ParticleType::LeftFermion,
                                     ParticleType::RightFermion,
                                     ParticleType::LeftFermion,
                                     ParticleType::Boson};

  const double D0T = Dlf[0](mass[0]);
  const double D0B = Dlf[0](mass[2]);
  const double D0H = Dlb[0](mass[3]);
  const double D2T = Dlf[2](mass[0]);
  const double D2B = Dlf[2](mass[2]);
  const double D2H = Dlb[2](mass[3]);

  const std::vector<double> gammatot = {D2T / (D0T * 7.1) * Tstar,
                                        D2T / (D0T * 7.6) * Tstar,
                                        D2B / (D0B * 7.1) * Tstar,
                                        D2H / (D0H * 14.) * Tstar};

  const std::vector<double> rates = {0.26 * pow(mass[0], 2) * Tstar, // gammaM
                                     5.9e-3 * Tstar,                 // gammaY
                                     gammatot[3],                    // gammaW
                                     1.5 * xW * xW * Tstar,          // gammaH
                                     2.7e-3 * Tstar};                // gammaSS

  std::vector<double> trunc = {1., 1., 1., 1., 1.};

  const std::vector<std::vector<double>> tL_interactions = {
      {2., 1., 1., 0., (1. + 9. * D0T)},
      {-2., -1., 0., 0., (-1. + 9. * D0T)},
      {0., 0., -1., 0., (1. + 9. * D0B)},
      {0., 1., 0., 0., 0.}};
  const std::vector<std::vector<double>> tR_interactions = {
      {-2., -1., 0., 0., -(1. + 9. * D0T)},
      {2., 2., 0., 0., (1. - 9. * D0T)},
      {0., -1., 0., 0., -(1. + 9. * D0B)},
      {0., -2., 0., 0., 0.}};
  const std::vector<std::vector<double>> bL_interactions = {
      {0., 0., -1., 0., (1. + 9. * D0T)},
      {0., -1., 0., 0., (-1. + 9. * D0T)},
      {0., 1., 1., 0., (1. + 9. * D0B)},
      {0., 1., 0., 0., 0.}};
  const std::vector<std::vector<double>> h_interactions = {
      {0., 1.5, 0., 0., 0.},
      {0., -3., 0., 0., 0},
      {0., 1.5, 0., 0., 0.},
      {0., 3., 0., 1., 0.}};
  const std::vector<std::vector<std::vector<double>>> prtcl = {
      tL_interactions, tR_interactions, bL_interactions, h_interactions};

  MatDoub Gamma(nEqs, nEqs, 0.);
  double m, Kp, ul, kappa;
  for (size_t i = 0; i < prtcl.size(); i++)
  {
    m     = mass[i];
    kappa = (ptype[i] == ParticleType::Boson ? KAPPABE : KAPPAFD);
    for (size_t l = 0; l < moment; l++)
    {
      if (l == 0)
      {
        ul                      = 0.;
        trunc[rates.size() - 1] = 1.;
        Kp = (ptype[i] == ParticleType::Boson ? Klb[l](m) : Klf[l](m));
      }
      else if (l == 1)
      {
        ul                      = 1.;
        trunc[rates.size() - 1] = 1.;
        Kp                      = -transportmodel->vwall *
             (ptype[i] == ParticleType::Boson ? Klb[0](m) : Klf[0](m));
      }
      else
      {
        ul                      = 1.;
        trunc[rates.size() - 1] = 0.;
        Kp = (ptype[i] == ParticleType::Boson ? Klb[l](m) : Klf[l](m));
      }

      for (size_t j = 0; j < prtcl.size(); j++)
      {
        for (size_t k = 0; k < rates.size(); k++)
        {
          Gamma[moment * i + l][moment * j] +=
              prtcl[i][j][k] * rates[k] * trunc[k];
        }
        Gamma[moment * i + l][moment * j] *= kappa * Kp;
      }
      Gamma[moment * i + l][moment * i + l] -= kappa * ul * gammatot[i];
    }
  }

  return Gamma;
}

std::vector<double> TransportEquations::calc_Ri(const size_t &particle,
                                                const size_t &k)
{
  std::vector<double> Ri(moment, 0);

  if (transportmodel->truncationscheme == TruncationScheme::Variance)
  {
    // R_1, ..., R_n-1, -1
    double n = static_cast<double>(moment); // This conversion is NECESSARY!
    Ri.at(0) = -n * n * pow(-Solution[particle * moment + 1][k], moment - 1);
    for (size_t i = 2; i < moment; i++)
    {
      double nck = nChoosek(moment, i);
      Ri.at(0) += nck * (n - i) * Solution[particle * moment + i][k] *
                  pow(-1 * Solution[particle * moment + 1][k], moment - i - 1);
      Ri.at(i - 1) =
          -nck * pow(-Solution[particle * moment + 1][k], moment - i);
    }
  }
  else
  {
    // 0, 0, ..., const, -1
    Ri.at(moment - 2) = transportmodel->truncationR;
  }
  Ri.at(moment - 1) = -1;

  return Ri;
}

MatDoub TransportEquations::calc_Ainv(const double &x,
                                      const ParticleType &type,
                                      const size_t &particle,
                                      const size_t &k)
{
  MatDoub res(moment, moment, 0.);
  std::vector<double> Di(moment, 1), Ri = calc_Ri(particle, k);

  for (size_t i = 0; i < moment - 1; i++)
  {
    // 1, D1, D2, ..., Dn-1
    Di.at(i + 1) =
        (type == ParticleType::Boson ? Dlb[i + 1](x) : Dlf[i + 1](x));
    // off-diagonal "diagonal"
    res[i + 1][i] = 1;
  }

  // (-1)^momentInv(A)
  double Dn = (type == ParticleType::Boson ? Dlb[moment](x) : Dlf[moment](x));
  for (size_t i = 0; i < moment - 1; i++)
    Dn -= Ri.at(i) * Di.at(i + 1);

  // Tensor product
  for (std::size_t i = 0; i < moment; i++)
    for (std::size_t j = 0; j < moment; j++)
      res[i][j] += Di[i] * Ri[j] / Dn;

  return res;
}

MatDoub TransportEquations::calc_m2B(const double &x,
                                     const double &dx2,
                                     const ParticleType &type)
{
  MatDoub res(moment, moment, 0.);
  const double fRbar = (type == ParticleType::Boson ? Rbarb(x) : Rbarf(x));

  for (size_t l = 1; l <= moment; l++)
  {
    const double fQ   = (type == ParticleType::Boson ? Qlb[l](x) : Qlf[l](x));
    res[l - 1][l - 1] = fRbar * (double)(l - 1) * dx2;
    res[l - 1][0]     = gamwall * transportmodel->vwall * fQ * dx2;
  }
  return res;
}

VecDoub TransportEquations::calc_source(const double &x,
                                        const double &dx2,
                                        const double &dth,
                                        const double &d2th,
                                        const int &h)
{
  VecDoub res(moment);
  for (size_t i = 0; i < moment; i++)
    res[i] = -transportmodel->vwall * gamwall * h *
             ((dx2 * dth + x * x * d2th) * Q8ol[i + 1](x) -
              dx2 * x * x * dth * Q9ol[i + 1](x)); // Si
  return res;
}

void TransportEquations::InsertBlockDiagonal(MatDoub &full,
                                             MatDoub &sub,
                                             const size_t position)
{
  assert((full.cols() % sub.cols()) == 0);
  size_t start = position * sub.cols();
  for (size_t i = 0; i < sub.rows(); i++)
    for (size_t j = 0; j < sub.cols(); j++)
      full[start + i][start + j] = sub[i][j];
}

void TransportEquations::Equations(const double &z,
                                   MatDoub &Mtilde,
                                   VecDoub &Stilde,
                                   const size_t &k)
{
  MatDoub Ainverse(nEqs, nEqs, 0.); // Store A^-1

  VecDoub S(nEqs, 0.); // Source vector

  MatDoub m2B(nEqs, nEqs, 0.); // m2' B;

  Stilde.zero(); // A^-1 * Source vector
  Mtilde.zero(); // Store A^-1 * M

  // Mass vector
  const double xW = transportmodel->GetWRatio(z);
  VecDoub FermionMasses(nFermions);
  VecDoub BosonMasses(nBosons);

  for (size_t particle = 0; particle < transportmodel->involvedparticles.size();
       particle++)
  {
    ParticleType ptype = transportmodel->involvedparticles[particle];
    double x2, x2prime, thetaprime, theta2prime;
    if (ptype == ParticleType::Boson)
    {
      x2             = 0;
      x2prime        = 0;
      BosonMasses[0] = x2;
    }
    else
    {
      // Calculate fermionic masses
      transportmodel->GetFermionRatio(
          z, particle, x2, x2prime, thetaprime, theta2prime);
      FermionMasses[particle] = x2;

      // Source terms
      int h = (ptype == ParticleType::LeftFermion ? -1 : 1);
      VecDoub tempS(calc_source(sqrt(x2), x2prime, thetaprime, theta2prime, h));
      for (size_t i = 0; i < moment; i++)
        S[moment * particle + i] = tempS[i];
    }

    // Calculate A inverse for fermion
    MatDoub tempA(calc_Ainv(sqrt(x2), ptype, particle, k));
    InsertBlockDiagonal(Ainverse, tempA, particle);

    // Calculate m2'B for fermion
    MatDoub tempB(calc_m2B(sqrt(x2), x2prime, ptype));
    InsertBlockDiagonal(m2B, tempB, particle);
  }

  // Gamma = deltaC - m2' B
  const MatDoub CollisionMatrix =
      CalculateCollisionMatrix(xW, FermionMasses, BosonMasses);

  // Calculate M = A^-1 * Gamma ( = deltaC - m2'B)
  for (size_t i = 0; i < nEqs; i++)
    for (size_t j = 0; j < nEqs; j++)
      for (size_t l = 0; l < nEqs; l++)
        Mtilde[i][j] += Ainverse[i][l] * (CollisionMatrix[l][j] - m2B[l][j]);

  // Calculate Stilde = A^-1 * S
  for (size_t i = 0; i < nEqs; i++)
    for (size_t j = 0; j < moment * (nFermions); j++)
      Stilde[i] += Ainverse[i][j] * S[j];
}

void TransportEquations::MakeDistribution(const double amplitude,
                                          const size_t npoints)
{
  zList.clear();
  for (size_t i = 0; i < npoints; i++)
  {
    // linear -1 -> 1
    double u = ((i - (npoints - 1) / 2.)) / ((npoints - 1) / 2.);
    u        = pow(u, 3);
    u        = amplitude * atanh(u * tanh(2));
    zList.push_back(u);
  }
}

void TransportEquations::CheckBoundary()
{
  double HighestNegRe, HighestNegEigenvalue;
  CheckBoundary(HighestNegRe, HighestNegEigenvalue);
}

void TransportEquations::CheckBoundary(double &HighestNegReEigenvalue,
                                       double &HighestImEigenvalue)
{
  // Initialize vars
  HighestNegReEigenvalue = -1e100;
  HighestImEigenvalue    = 3;

  // Calculate asymptotic matrixs
  Equations(zList.front(), MtildeM, StildeM, 0);
  Equations(zList.back(), MtildeP, StildeP, NumberOfSteps - 1);

  double STildeLength(0);
  size_t NumberOfNonDecayingModes(0);
  std::stringstream ss;

  Eigen::MatrixXcd EigenMtildeM(nEqs, nEqs);
  Eigen::MatrixXcd EigenMtildeP(nEqs, nEqs);

  for (size_t i = 0; i < nEqs; i++)
  {
    STildeLength += pow(StildeM[i], 2);
    STildeLength += pow(StildeP[i], 2);
    for (size_t j = 0; j < nEqs; j++)
    {
      EigenMtildeM(i, j) = MtildeM[i][j];
      EigenMtildeP(i, j) = MtildeP[i][j];
    }
  }

  // Checking if S vector are small enough
  STildeLength = sqrt(STildeLength);
  STildeLength /= 2. * nEqs;

  if (STildeLength > STildeThreshold)
  {
    Status = FHCKStatus::SmallIntegrationRegion;
    ss << "\033[31mThe integration region is too small, increase "
          "LwMultiplier.\033[0m\n";
  }

  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> EigenSolverM(EigenMtildeM);
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> EigenSolverP(EigenMtildeP);

  const double EigenValueThreshold = -1e-10;
  // Number of modes that have to be set to zero
  for (auto ev : EigenSolverM.eigenvalues())
    if (-ev.real() >= EigenValueThreshold) /* minus sign z -> -Infinity */
      NumberOfNonDecayingModes++;
    else
    {
      HighestNegReEigenvalue =
          std::max(HighestNegReEigenvalue,
                   -ev.real()); /* minus sign because z -> -Infinity */
      HighestImEigenvalue = std::max(HighestImEigenvalue, std::abs(ev.imag()));
    }
  for (auto ev : EigenSolverP.eigenvalues())
    if (ev.real() >= EigenValueThreshold)
      NumberOfNonDecayingModes++;
    else
    {
      HighestNegReEigenvalue = std::max(HighestNegReEigenvalue, ev.real());
      HighestImEigenvalue = std::max(HighestImEigenvalue, std::abs(ev.imag()));
    }

  if (NumberOfNonDecayingModes > nEqs)
  {
    ss << " \033[31m\nToo many non-decaying modes (" << NumberOfNonDecayingModes
       << ") for the boundary conditions (" << nEqs
       << "). Impossible to satisfy "
          "the "
          "boundary "
          "conditions of μ = u = 0.\033[0m\n";
    Status = FHCKStatus::UnphysicalBoundary;
  }

  Logger::Write(LoggingLevel::FHCK, ss.str());
}

void TransportEquations::smatrix(const int k,
                                 const int k1,
                                 const int k2,
                                 const int jsf,
                                 VecInt &indexv,
                                 MatDoub &s,
                                 MatDoub &y)
{
  double temp, zlow, zhigh;
  if (transportmodel->VevProfile == VevProfileMode::FieldEquation)
  {
    zlow  = transportmodel->vacuumprofile->z.front();
    zhigh = transportmodel->vacuumprofile->z.back();
  }
  else
  {
    zlow  = -100. * transportmodel->Lw;
    zhigh = 100. * transportmodel->Lw;
  }
  MatDoub Mtilde(nEqs, nEqs);
  VecDoub Stilde(nEqs);
  s.zero(); // Set matrix s = 0
  if (k == k1)
  {
    // Boundary conditions mu = 0 and first u on first boundary
    for (size_t i = 0; i < nEqs / 2; i++)
    {
      // Sn at the first boundary
      s[nEqs / 2 + i][nEqs + i] = 1.0;
      // B0
      s[nEqs / 2 + i][jsf] = y[indexv[i]][0];
    }
  }
  else if (k > k2 - 1)
  {
    // Boundary conditions mu = 0 and first u on second boundary
    for (size_t i = 0; i < nEqs / 2; i++)
    {
      // Sn at the last boundary
      s[i][nEqs + i] = 1.0;
      // C0
      s[i][jsf] = y[indexv[i]][zList.size() - 1];
    }
  }
  else
  {
    double zk = (zList[k] + zList[k - 1]) / 2.;
    double dz = (zList[k] - zList[k - 1]);
    if (zk < zlow)
    {
      Mtilde = MtildeM;
      Stilde = StildeM;
    }
    else if (zk > zhigh)
    {
      Mtilde = MtildeP;
      Stilde = StildeP;
    }
    else
      Equations(zk, Mtilde, Stilde, k);
    for (size_t j = 0; j < nEqs; j++)
    {
      for (size_t n = 0; n < nEqs; n++)
      {
        // s matrix for the middle point
        // S_{j,n}
        // M[k] and S[k] are evaluated between k and k-1
        s[j][indexv[n]] = -Delta(j, n) - 0.5 * dz * (Mtilde[j][n]);
        // S_{j,N + n}
        s[j][nEqs + indexv[n]] = Delta(j, n) - 0.5 * dz * (Mtilde[j][n]);
      }
      //  Equations for E(k,k-1)
      temp = Stilde[j];
      for (size_t i = 0; i < nEqs; i++)
        temp += Mtilde[j][i] * (y[i][k] + y[i][k - 1]) / 2.;
      s[j][jsf] = y[j][k] - y[j][k - 1] - dz * temp;
    }
  }
}

void TransportEquations::SolveTransportEquation()
{
  std::fill(BAUeta.begin(), BAUeta.end(), NAN);
  if (transportmodel->status != TransportModelStatus::Success)
  {

    Logger::Write(LoggingLevel::FHCK,
                  "SolveTransportEquation() failed. No VEV profile available.");
    return;
  }

  Logger::Write(LoggingLevel::FHCK,
                "Lw = " + std::to_string(transportmodel->Lw));

  for (size_t ell = 0; ell < moments.size(); ell++)
  {
    stringstream ss;
    ss << "---------------- Calculating BAU for ℓ = " << moments.at(ell)
       << " ----------------";
    Logger::Write(LoggingLevel::FHCK, ss.str());
    ss.str("");

    double bauLowPrecision = NAN, bauHighPrecision = NAN;

    StepsPerCycle   = StepsPerCycleLow;
    bauLowPrecision = SolveTransportEquationEll(ell);

    if (not std::isnan(bauLowPrecision))
    {
      StepsPerCycle    = StepsPerCycleHigh;
      bauHighPrecision = SolveTransportEquationEll(ell);
    }

    const double Uncertainty = std::abs(bauHighPrecision / bauLowPrecision - 1);

    PrintTransportEquation(120, "tL", "mu");
    PrintTransportEquation(120, "tR", "mu");
    PrintTransportEquation(120, "bL", "mu");
    PrintTransportEquation(120, "h", "mu");

    ss << "BAU(low precision) = " << bauLowPrecision
       << "\nBAU(high precision) = " << bauHighPrecision
       << "\nUncertainty = " << Uncertainty;

    Logger::Write(LoggingLevel::FHCK, ss.str());

    if (not isnan(bauLowPrecision) and not isnan(bauHighPrecision) and
        Uncertainty < UncertaintyThreshold)
      BAUeta.at(ell) = bauHighPrecision;
  }
}

double TransportEquations::SolveTransportEquationEll(const size_t &ell)
{
  Status = FHCKStatus::NotSet;

  InitializeMoment(moments.at(ell));

  Equations(zList.front(), MtildeM, StildeM, 0);
  Equations(zList.back(), MtildeP, StildeP, NumberOfSteps - 1);

  CheckBoundary();
  if (Status != FHCKStatus::NotSet)
  {
    Logger::Write(
        LoggingLevel::FHCK,
        "\033[31mTransport equation unable to be computed. Status code : " +
            FHCKStatusToString.at(FHCKStatus::NotSet) + "\033[0m\n");
    return NAN;
  }

  size_t itmax = 10;
  double conv  = 1e-10;
  double slowc = 1.;
  VecDoub scalv(nEqs, 1);
  VecInt indexv(nEqs);

  // fix μ and first u
  for (size_t l = 0; l < moment /* moment= 2 + 4k */; l++)
    for (size_t particle = 0; particle < nParticles; particle++)
      indexv[l + particle * moment] = particle + l * (nParticles);

  int NB = nEqs / 2;

  size_t it;
  double MinError    = 1.e100;
  double MaxError    = 1e-100;
  size_t NotBetter   = 0;
  auto Best_Solution = Solution;

  for (it = 0; it < itmax; it++)
  {
    if (NotBetter >= NotBetterThreshold) break;
    try
    {
      RelaxOde solvde(1, conv, slowc, scalv, indexv, NB, Solution, *this);
      stringstream ss;
      ss << "[Transport equations] it = " << it << ". Error = " << solvde.err;
      Logger::Write(LoggingLevel::FHCK, ss.str());

      MaxError = std::max(MaxError, solvde.err);

      if (solvde.err < MinError)
      {
        MinError      = solvde.err;
        Best_Solution = Solution;
        NotBetter     = 0;
      }

      // Early exit in case error gets much smaller than the worst case
      if (solvde.err / MaxError < 1e-10) break;

      NotBetter++;
    }
    catch (const char *e)
    {
      Logger::Write(LoggingLevel::FHCK,
                    "\033[31mException caught in solvde: " + std::string(e) +
                        "\033[0m\n");
      return NAN;
    }
  }

  if (it == NotBetter)
  {
    Logger::Write(
        LoggingLevel::FHCK,
        "Relaxation transport equations failed! [did not converge once]");
    return NAN;
  }

  Solution = Best_Solution;

  CalculateBAU();

  return bau;
}

double TransportEquations::Gws(const double &z)
{
  const double h    = transportmodel->EWSBVEV(z);
  const double Gsph = 8.e-7 * Tstar;
  return Gsph * std::min(1., 1.7 * Tstar / Gsph * std::exp(-37. * h / Tstar));
}

double TransportEquations::WashoutFactor(const double &z)
{
  const double A   = 15. / 2.;
  const double fac = -A * nf / (2 * transportmodel->vwall * gamwall);

  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(10000);

  gsl_function F;
  F.function = [](double zz, void *params) -> double
  { return static_cast<TransportEquations *>(params)->Gws(zz); };
  F.params = this;

  double result, error;
  gsl_integration_qag(&F,
                      zList.front(),
                      z,
                      1e-10,
                      1e-10,
                      10000,
                      GSL_INTEG_GAUSS51,
                      workspace,
                      &result,
                      &error);

  gsl_integration_workspace_free(workspace);

  return std::exp(fac * result);
}

void TransportEquations::CalculateBAU()
{
  const double Nc = 3;

  double xt2, xb2, x2prime, thetaprime, theta2prime; // temporary vars
  // Weak spharelon rate
  double r;                // temporary variable to store the result
  std::vector<double> z;   // list of u positions
  std::vector<double> muB; // muB integrand at position u
  for (size_t i = 0; i < zList.size(); i++)
  {
    const double zi = zList.at(i);
    // Calculate the vev at z
    transportmodel->GetFermionRatio(
        zi, 0, xt2, x2prime, thetaprime, theta2prime);
    transportmodel->GetFermionRatio(
        zi, 2, xb2, x2prime, thetaprime, theta2prime);
    // Results
    r = 0;
    r += (1 + 4 * Dlf[0](sqrt(xt2))) / 2. * Solution[0][i];          // tL
    r += (1 + 4 * Dlf[0](sqrt(xb2))) / 2. * Solution[moment * 2][i]; // bL
    r += 2. * Dlf[0](sqrt(xt2)) * Solution[moment][i];               // tR
    r *= Gws(zi);
    r *= WashoutFactor(zi);
    // Save in list to pass to integrator
    z.push_back(zi);
    muB.push_back(r); // integrand
  }

  size_t n = z.size();
  double a = z.front(); // Lower bound of integration
  double b = z.back();  // Upper bound of integration

  // linear approximation
  gsl_interp_accel *acc_l = gsl_interp_accel_alloc();
  gsl_spline *spline_l    = gsl_spline_alloc(gsl_interp_linear, n);
  gsl_spline_init(spline_l, z.data(), muB.data(), n);
  double result_l = gsl_spline_eval_integ(spline_l, a, b, acc_l);
  gsl_spline_free(spline_l);
  gsl_interp_accel_free(acc_l);

  // cubic spline approximation
  gsl_interp_accel *acc_c = gsl_interp_accel_alloc();
  gsl_spline *spline_c    = gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline_init(spline_c, z.data(), muB.data(), n);
  double result_c = gsl_spline_eval_integ(spline_c, a, b, acc_c);
  gsl_spline_free(spline_c);
  gsl_interp_accel_free(acc_c);

  double error = abs(result_c - result_l);

  // Results
  const double prefactor = -45. * nf * Nc /
                           (4 * M_PI * M_PI * transportmodel->vwall * gamwall *
                            106.75 * Tstar); // TODO Fix gstas

  result_c *= prefactor;
  error *= prefactor;

  // Step 3: Output the result
  stringstream ss;
  ss << "eta = " << result_c << " with error " << error << std::endl;

  double Uncertainty = abs(error / result_c);

  if (Uncertainty > UncertaintyThreshold)
  {
    ss << "Calculation failed!\t" << Uncertainty << "\n";
    result_c = NAN;
  }

  // Save the result
  bau = result_c;

  ss << "eta/eta_obs = " << result_c / (8.7e-11);
  Logger::Write(LoggingLevel::FHCK, ss.str());
}

void TransportEquations::PrintTransportEquation(const int &size,
                                                const std::string &Particle,
                                                const std::string &MuOrU)
{
  AsciiPlotter Plot(Particle + " " + MuOrU, size, ceil(size / 3.));
  std::optional<int> ind;
  std::vector<double> z, y;

  if (Particle == "tL") ind = 0 * moment;
  if (Particle == "tR") ind = 1 * moment;
  if (Particle == "bL") ind = 2 * moment;
  if (Particle == "h") ind = 3 * moment;

  if (not ind.has_value()) throw("Invalid particle to plot the solution.");

  if (MuOrU == "u") ind = ind.value() + 1;

  size_t i_min_left = 0, i_min_right = zList.size() - 1;
  double max = -1.;
  for (size_t i = 0; i < zList.size(); i++)
    max = std::max(max, abs(Solution[ind.value()][i]));

  for (size_t i = 0; i < zList.size(); i++)
  {
    if (abs(Solution[ind.value()][i]) > max / 100.)
    {
      i_min_left = i;
      break;
    }
  }

  for (size_t i = zList.size(); i-- > 0;)
  {
    if (abs(Solution[ind.value()][i]) > max / 100.)
    {
      i_min_right = i;
      break;
    }
  }

  for (size_t i = i_min_left; i <= i_min_right; i++)
  {
    z.push_back(zList[i]);
    y.push_back(Solution[ind.value()][i]);
  }
  Plot.addPlot(z, y, "", '*');
  std::stringstream ss;
  Plot.show(ss);
  Logger::Write(LoggingLevel::FHCK, ss.str());
}
} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT
