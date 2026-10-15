/**
 * @file
 */

#include <BSMPT/vacuum_profile/vacuum_profile.h>

namespace BSMPT
{
namespace VacuumProfileNS
{

VecInt VacuumProfile::Calcindexv()
{
  VecInt indexv(2 * dim);
  for (size_t i = 0; i < 2 * dim; i++)
  {
    // 0 1 2 3 4 5...
    indexv[i] = i;
    // switch index of field and deriv
    // 0 1 2 3 4 5-> 3 4 5 1 2 3
    indexv[i] += dim * ((i < dim) * 2 - 1) *
                 /* reordeing only necessary for dirichlet */
                 (mode == ProfileSolverMode::Field);
  }
  return indexv;
}

void VacuumProfile::CreatePath(const size_t &NumberOfSteps)
{
  // temporary variables
  std::vector<std::vector<double>> path;
  std::vector<double> vev_k(dim), zpath;
  // Use the kink to calculate a first approximation
  Logger::Write(LoggingLevel::VacuumProfile,
                "Creating path | Lw = " + std::to_string(Lw) +
                    " and NumberOfSteps = " + std::to_string(NumberOfSteps));
  for (size_t k = 0; k < NumberOfSteps; k++)
  {
    // uniformly distributed between -10 Lw and 10 Lw
    double u        = (-1. + 2. * k / (NumberOfSteps - 1.));
    u               = pow(u, 3);
    u               = atanh(u * tanh(1));
    const double zk = u * Lw * LwToSolve;
    for (size_t i = 0; i < dim; i++)
      vev_k.at(i) =
          TrueVacuum.at(i) +
          (1. + tanh(zk / Lw)) / 2. * (FalseVacuum.at(i) - TrueVacuum.at(i));
    zpath.push_back(zk);
    path.push_back(vev_k);
  }
  LoadPath(zpath, path);
}

void VacuumProfile::LoadPath(const std::vector<double> &z_In,
                             const std::vector<std::vector<double>> &path_In)
{
  if (z_In.size() != path_In.size())
    throw("z and path have different lengths.");
  scalv = VecDoub(2 * dim, 1); // min of 1
  z     = z_In;
  y     = MatDoub(dim * 2, z.size(), 0.);
  for (size_t k = 0; k < path_In.size(); k++)
  {
    // load z position
    z[k] = z_In[k];
    for (size_t i = 0; i < dim; i++)
    {
      // load field
      y[dim + i][k] = path_In[k][i];
      // load field z-derivative
      if (k != 0 and k != path_In.size() - 1) // ignore first and last point
        y[i][k] = (path_In[k + 1][i] - path_In[k - 1][i]) /
                  (z_In[k + 1] - z_In[k - 1]);
      // calculate scale of quantities
      if (abs(y[i][k]) > scalv[i]) scalv[i] = abs(y[i][k]);
      if (abs(y[dim + i][k]) > scalv[dim + i])
        scalv[dim + i] = abs(y[dim + i][k]);
    }
  }
}

void VacuumProfile::CalculateProfile()
{
  size_t it;
  double conv  = 1e-10;
  double slowc = 1e-1;

  VecInt indexv = Calcindexv();
  std::stringstream ss;

  ss << "---------------------------------------\nStarting vacuum profile "
        "calculation\n---------------------------------------\n";

  ss << "True vacuum = " << TrueVacuum << "\n";
  ss << "False vacuum = " << FalseVacuum << "\n";
  ss << "VEV dimension = " << dim << "\n";
  ss << "Calculating profile in z ∈ [" << z.front() << ", " << z.back() << "]";
  Logger::Write(LoggingLevel::VacuumProfile, ss.str());
  ss.clear();

  Difeq_VacuumProfile difeq_vacuumprofile(
      mode, dim, z, TrueVacuum, FalseVacuum, V, dV, Hessian);

  double MinError  = 1.e100;
  size_t NotBetter = 0;
  auto Best_y      = y;

  double oldmu = 1e100;

  for (it = 0; it < itmax; it++)
  {
    CenterPath();
    if (NotBetter >= NotBetterThreshold) break;
    RelaxOde solvde(1, conv, slowc, scalv, indexv, dim, y, difeq_vacuumprofile);
    std::stringstream sss;
    sss << "[Relaxation Vacuum profile] it = " << it
        << ". Error = " << solvde.err << ". mu = " << difeq_vacuumprofile.mu
        << " | relative variation of mu = "
        << std::abs(difeq_vacuumprofile.mu / oldmu - 1.);

    oldmu = difeq_vacuumprofile.mu;

    Logger::Write(LoggingLevel::VacuumProfile, sss.str());
    if (solvde.err < MinError)
    {
      MinError  = solvde.err;
      Best_y    = y;
      NotBetter = 0;
    }
    // Check if edges are moving
    double dphi = 0;
    for (size_t i = 0; i < dim; i++)
      dphi +=
          pow(y[i][0] / scalv[i], 2) + pow(y[i][z.size() - 1] / scalv[i], 2);
    if (sqrt(dphi) / (2 * dim) > dphiTreshold)
    {

      if (LwToSolve > 100)
      {
        Logger::Write(LoggingLevel::VacuumProfile,
                      "\nVacuum profile calculation failed! [small domain]\n");
        status = VacuumProfileStatus::Failed;
        return;
      }
      LwToSolve += 10;
      // Generate new path. Keep roughly same number of points
      CreatePath(ceil(z.size() * LwToSolve / (LwToSolve - 10)));

      ss << "\n\033[31mVacuum profile small domain (dphi = "
         << sqrt(dphi) / (2 * dim) << "), increase LwToSolve to " << LwToSolve
         << "\033[0m";
      Logger::Write(LoggingLevel::VacuumProfile, ss.str());
      CalculateProfile();
      return;
    }

    NotBetter++;
  }

  if (it == NotBetter)
  {
    Logger::Write(LoggingLevel::VacuumProfile,
                  "Vacuum profile calculation failed! [did not converge once]");
    status = VacuumProfileStatus::Failed;
    return;
  }

  ss << "\nμ = " << difeq_vacuumprofile.mu << "\n";

  status = VacuumProfileStatus::Success;
  y      = Best_y;
  CenterPath();

  Logger::Write(LoggingLevel::VacuumProfile, ss.str());
}

double VacuumProfile::CalculateWidth(
    const std::vector<double> &TrueVacuum,
    const std::vector<double> &FalseVacuum,
    const std::function<double(std::vector<double>)> &V)
{
  size_t NumberPointsBarrier = 100;

  const double vc     = BSMPT::L2NormVector(TrueVacuum - FalseVacuum);
  const double Vtrue  = V(TrueVacuum);
  const double Vfalse = V(FalseVacuum);
  double Vb           = -1.;
  for (size_t k = 0; k < NumberPointsBarrier; k++)
  {
    const double Vk = V(TrueVacuum + (double)k / (NumberPointsBarrier - 1.) *
                                         (FalseVacuum - TrueVacuum));
    Vb              = std::max(Vb, std::max(abs(Vtrue - Vk), abs(Vfalse - Vk)));
  }

  const double Lw = vc / sqrt(8 * Vb);

  Logger::Write(LoggingLevel::VacuumProfile,
                "Estimating Lw using a rough approximation -> Lw = " +
                    std::to_string(Lw));

  return Lw;
}

void VacuumProfile::GenerateSplines()
{
  splines.clear();
  for (size_t i = 0; i < dim; i++)
  {
    std::vector<double> phi;
    for (size_t k = 0; k < z.size(); k++)
      phi.push_back(y[dim + i][k]);
    tk::spline s(z,
                 phi,
                 tk::spline::cspline,
                 false,
                 tk::spline::not_a_knot,
                 0.0,
                 tk::spline::not_a_knot,
                 0.0);
    splines.push_back(s);
  }
}

std::vector<double> VacuumProfile::GetVev(const double &zz, const int &diff)
{
  if (splines.size() != dim)
  {
    Logger::Write(LoggingLevel::VacuumProfile,
                  "Vacuum profile calculation failed. Do not call GetVev(z)");
    return std::vector<double>(dim, 0); // return 0, ..., 0
  }
  if (zz < z.front())
  {
    if (diff == 0)
      return TrueVacuum;
    else
      return std::vector<double>(dim, 0.);
  }
  if (zz > z.back())
  {
    if (diff == 0)
      return FalseVacuum;
    else
      return std::vector<double>(dim, 0.);
  }
  std::vector<double> r;
  for (size_t i = 0; i < dim; i++)
  {
    if (diff == 0)
    {
      r.push_back(splines.at(i)(zz));
    }
    else
    {
      r.push_back(splines.at(i).deriv(diff, zz));
    }
  }
  return r;
}

void VacuumProfile::CenterPath()
{
  double center;
  CenterPath(center);
}

void VacuumProfile::CenterPath(double &center)
{
  GenerateSplines();
  // locate maximum of dphi/dz
  double max_dphidz = -1;
  for (size_t k = 0; k < z.size(); k++)
  {
    double dphidz = 0;
    for (size_t i = 0; i < dim; i++)
      dphidz += pow(y[i][k], 2);
    if (dphidz > max_dphidz)
    {
      center     = z[k];
      max_dphidz = dphidz;
    }
  }
  std::vector<std::vector<double>> new_path;
  for (const double &zk : z)
  {
    new_path.push_back(GetVev(zk + center));
  }
  LoadPath(z, new_path);
  GenerateSplines();
}

VacuumProfile::VacuumProfile(
    // Dimension of VEV space
    const size_t dim_In,
    // True
    const std::vector<double> &TrueVacuum_In,
    // False
    const std::vector<double> &FalseVacuum_In,
    // Potential
    const std::function<double(std::vector<double>)> &V_In,
    // Gradient
    const std::function<std::vector<double>(std::vector<double>)> &dV_In,
    // Hessian
    const std::function<std::vector<std::vector<double>>(std::vector<double>)>
        &Hessian_In,
    // z knot list
    std::vector<double> &z_In,
    // field knot list
    std::vector<std::vector<double>> &path_In)
    : dim(dim_In)
    , TrueVacuum(TrueVacuum_In)
    , FalseVacuum(FalseVacuum_In)
    , V(V_In)
    , dV(dV_In)
    , Hessian(Hessian_In)
{
  LoadPath(z_In, path_In);
}

VacuumProfile::VacuumProfile(
    // Dimension of VEV space
    const size_t dim_In,
    // True
    const std::vector<double> &TrueVacuum_In,
    // False
    const std::vector<double> &FalseVacuum_In,
    // Potential
    const std::function<double(std::vector<double>)> &V_In,
    // Gradient
    const std::function<std::vector<double>(std::vector<double>)> &dV_In,
    // Hessian
    const std::function<std::vector<std::vector<double>>(std::vector<double>)>
        &Hessian_In,
    // Bubble width
    const double Lw_In,
    // Number of steps for the path
    const size_t NumberOfSteps)
    : dim(dim_In)
    , Lw(Lw_In)
    , TrueVacuum(TrueVacuum_In)
    , FalseVacuum(FalseVacuum_In)
    , V(V_In)
    , dV(dV_In)
    , Hessian(Hessian_In)
{
  CreatePath(NumberOfSteps);
}

VacuumProfile::VacuumProfile(
    // Dimension of VEV space
    const size_t dim_In,
    // True
    const std::vector<double> &TrueVacuum_In,
    // False
    const std::vector<double> &FalseVacuum_In,
    // Potential
    const std::function<double(std::vector<double>)> &V_In,
    // Gradient
    const std::function<std::vector<double>(std::vector<double>)> &dV_In,
    // Hessian
    const std::function<std::vector<std::vector<double>>(std::vector<double>)>
        &Hessian_In)
    : VacuumProfile(dim_In,
                    TrueVacuum_In,
                    FalseVacuum_In,
                    V_In,
                    dV_In,
                    Hessian_In,
                    CalculateWidth(TrueVacuum_In, FalseVacuum_In, V_In)) {};

} // namespace VacuumProfileNS
} // namespace BSMPT