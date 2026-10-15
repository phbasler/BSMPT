#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Approx = Catch::Approx;

#include <BSMPT/baryo_fhck/Kernels.h>
#include <BSMPT/baryo_fhck/TransportEquations.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/NumericalDerivatives.h>
#include <BSMPT/utility/utility.h>
#include <BSMPT/vacuum_profile/vacuum_profile.h>

TEST_CASE("Domain Wall lambda^4", "[vacuumprofile]")
{
  using namespace BSMPT;
  using namespace Baryo::FHCK;
  using namespace VacuumProfileNS;

  SetLogger({"--logginglevel::complete=true"});

  size_t dim                      = 1;
  std::vector<double> TrueVacuum  = {-1};
  std::vector<double> FalseVacuum = {1};

  double m   = 0;
  double lam = 1;

  std::function<double(std::vector<double>)> V = [=](auto const &arg)
  {
    return lam * pow(arg[0], 2) * (-2 + pow(arg[0], 2)) +
           m * (arg[0] - pow(arg[0], 3) / 3.);
  };
  std::function<std::vector<double>(std::vector<double>)> dV =
      [=](auto const &arg)
  {
    return std::vector<double>(
        {-((m - 4 * lam * arg[0]) * (-1 + pow(arg[0], 2)))});
  };
  std::function<std::vector<std::vector<double>>(std::vector<double>)> Hessian =
      [=](auto const &arg)
  {
    return std::vector<std::vector<double>>(
        {{-2 * m * arg[0] + 10 * lam * pow(arg[0], 2) +
          2 * lam * (-2 + pow(arg[0], 2))}});
  };

  VacuumProfile vacuumprofile(dim, TrueVacuum, FalseVacuum, V, dV, Hessian);

  vacuumprofile.CalculateProfile();

  REQUIRE(1 == 1);
}
TEST_CASE("Bubble profile lambda^4", "[vacuumprofile]")
{
  using namespace BSMPT;
  using namespace Baryo::FHCK;
  using namespace VacuumProfileNS;

  SetLogger({"--logginglevel::complete=true"});

  size_t dim                      = 1;
  std::vector<double> TrueVacuum  = {-1};
  std::vector<double> FalseVacuum = {1};

  double m   = 0.1;
  double lam = 1;

  std::function<double(std::vector<double>)> V = [=](auto const &arg)
  {
    return lam * pow(arg[0], 2) * (-2 + pow(arg[0], 2)) +
           m * (arg[0] - pow(arg[0], 3) / 3.);
  };
  std::function<std::vector<double>(std::vector<double>)> dV =
      [=](auto const &arg)
  {
    return std::vector<double>(
        {-((m - 4 * lam * arg[0]) * (-1 + pow(arg[0], 2)))});
  };
  std::function<std::vector<std::vector<double>>(std::vector<double>)> Hessian =
      [=](auto const &arg)
  {
    return std::vector<std::vector<double>>(
        {{-2 * m * arg[0] + 10 * lam * pow(arg[0], 2) +
          2 * lam * (-2 + pow(arg[0], 2))}});
  };

  VacuumProfile vacuumprofile(dim,
                              TrueVacuum,
                              FalseVacuum,
                              V,
                              dV,
                              Hessian,
                              1 / sqrt(3.) /* wrong on purpose*/,
                              10000);
  vacuumprofile.CalculateProfile();
  vacuumprofile.CenterPath();

  std::function<std::vector<double>(double)> solution = [=](auto const &arg)
  { return std::vector<double>({tanh(sqrt(2) * arg)}); };

  for (const auto z : std::vector<double>({{-0.5, -0.1, 0, 0.1, 0.5}}))
    REQUIRE(solution(z)[0] == Approx(vacuumprofile.GetVev(z)[0]).margin(0.001));
}

TEST_CASE("Test indexv - field", "[vacuumprofile]")
{
  using namespace BSMPT;
  using namespace Baryo::FHCK;
  using namespace VacuumProfileNS;

  size_t dim = 4;
  std::vector<double> TrueVacuum(dim, -1);
  std::vector<double> FalseVacuum(dim, 1);

  std::function<double(std::vector<double>)> V;
  std::function<std::vector<double>(std::vector<double>)> dV;
  std::function<std::vector<std::vector<double>>(std::vector<double>)> Hessian;

  VacuumProfile vacuumprofile(dim, TrueVacuum, FalseVacuum, V, dV, Hessian, 1);

  VecInt indexvField(8);

  indexvField[0] = 4;
  indexvField[1] = 5;
  indexvField[2] = 6;
  indexvField[3] = 7;
  indexvField[4] = 0;
  indexvField[5] = 1;
  indexvField[6] = 2;
  indexvField[7] = 3;

  vacuumprofile.mode = VacuumProfileNS::ProfileSolverMode::Field;
  VecInt indexv      = vacuumprofile.Calcindexv();
  for (size_t i = 0; i < 2 * dim; i++)
    CHECK(indexv[i] == indexvField[i]);
}

TEST_CASE("Test indexv - deriv", "[vacuumprofile]")
{
  using namespace BSMPT;
  using namespace Baryo::FHCK;
  using namespace VacuumProfileNS;

  size_t dim = 4;
  std::vector<double> TrueVacuum(dim, -1);
  std::vector<double> FalseVacuum(dim, 1);

  std::function<double(std::vector<double>)> V;
  std::function<std::vector<double>(std::vector<double>)> dV;
  std::function<std::vector<std::vector<double>>(std::vector<double>)> Hessian;

  VacuumProfile vacuumprofile(dim, TrueVacuum, FalseVacuum, V, dV, Hessian, 1);

  VecInt indexvDeriv(8);

  indexvDeriv[0] = 0;
  indexvDeriv[1] = 1;
  indexvDeriv[2] = 2;
  indexvDeriv[3] = 3;
  indexvDeriv[4] = 4;
  indexvDeriv[5] = 5;
  indexvDeriv[6] = 6;
  indexvDeriv[7] = 7;

  vacuumprofile.mode = VacuumProfileNS::ProfileSolverMode::Deriv;
  VecInt indexv      = vacuumprofile.Calcindexv();
  for (size_t i = 0; i < 2 * dim; i++)
    REQUIRE(indexv[i] == indexvDeriv[i]);
}

TEST_CASE("Test center path", "[vacuumprofile]")
{
  using namespace BSMPT;
  using namespace Baryo::FHCK;
  using namespace VacuumProfileNS;

  size_t dim                      = 1;
  std::vector<double> TrueVacuum  = {-1};
  std::vector<double> FalseVacuum = {1};

  double m   = 0;
  double lam = 1;

  std::function<double(std::vector<double>)> V = [=](auto const &arg)
  {
    return lam * pow(arg[0], 2) * (-2 + pow(arg[0], 2)) +
           m * (arg[0] - pow(arg[0], 3) / 3.);
  };
  std::function<std::vector<double>(std::vector<double>)> dV =
      [=](auto const &arg)
  {
    return std::vector<double>(
        {-((m - 4 * lam * arg[0]) * (-1 + pow(arg[0], 2)))});
  };
  std::function<std::vector<std::vector<double>>(std::vector<double>)> Hessian =
      [=](auto const &arg)
  {
    return std::vector<std::vector<double>>(
        {{-2 * m * arg[0] + 10 * lam * pow(arg[0], 2) +
          2 * lam * (-2 + pow(arg[0], 2))}});
  };

  VacuumProfile vacuumprofile(dim, TrueVacuum, FalseVacuum, V, dV, Hessian);

  vacuumprofile.CalculateProfile();

  double center;
  vacuumprofile.CenterPath(center); // Center path
  vacuumprofile.CenterPath(center); // Recenter does nothing

  REQUIRE(center == Approx(0.).margin(1e-2));
}

TEST_CASE("Test path small domain", "[vacuumprofile]")
{
  using namespace BSMPT;
  using namespace Baryo::FHCK;
  using namespace VacuumProfileNS;

  SetLogger({"--logginglevel::complete=true"});

  size_t dim                      = 1;
  std::vector<double> TrueVacuum  = {-1};
  std::vector<double> FalseVacuum = {1};

  double m   = 0;
  double lam = 1;

  std::function<double(std::vector<double>)> V = [=](auto const &arg)
  {
    return lam * pow(arg[0], 2) * (-2 + pow(arg[0], 2)) +
           m * (arg[0] - pow(arg[0], 3) / 3.);
  };
  std::function<std::vector<double>(std::vector<double>)> dV =
      [=](auto const &arg)
  {
    return std::vector<double>(
        {-((m - 4 * lam * arg[0]) * (-1 + pow(arg[0], 2)))});
  };
  std::function<std::vector<std::vector<double>>(std::vector<double>)> Hessian =
      [=](auto const &arg)
  {
    return std::vector<std::vector<double>>(
        {{-2 * m * arg[0] + 10 * lam * pow(arg[0], 2) +
          2 * lam * (-2 + pow(arg[0], 2))}});
  };

  VacuumProfile vacuumprofile(dim, TrueVacuum, FalseVacuum, V, dV, Hessian);

  vacuumprofile.LwToSolve = 2;
  vacuumprofile.CreatePath(1000);

  vacuumprofile.CalculateProfile();
  vacuumprofile.CenterPath();

  std::function<std::vector<double>(double)> solution = [=](auto const &arg)
  { return std::vector<double>({tanh(sqrt(2) * arg)}); };

  for (const auto z : std::vector<double>({{-0.5, -0.1, 0, 0.1, 0.5}}))
    REQUIRE(solution(z)[0] == Approx(vacuumprofile.GetVev(z)[0]).margin(0.001));
}

TEST_CASE("Test not enough iterations", "[vacuumprofile]")
{
  using namespace BSMPT;
  using namespace Baryo::FHCK;
  using namespace VacuumProfileNS;

  size_t dim                      = 1;
  std::vector<double> TrueVacuum  = {-1};
  std::vector<double> FalseVacuum = {1};

  double m   = 0;
  double lam = 1;

  std::function<double(std::vector<double>)> V = [=](auto const &arg)
  {
    return lam * pow(arg[0], 2) * (-2 + pow(arg[0], 2)) +
           m * (arg[0] - pow(arg[0], 3) / 3.);
  };
  std::function<std::vector<double>(std::vector<double>)> dV =
      [=](auto const &arg)
  {
    return std::vector<double>(
        {-((m - 4 * lam * arg[0]) * (-1 + pow(arg[0], 2)))});
  };
  std::function<std::vector<std::vector<double>>(std::vector<double>)> Hessian =
      [=](auto const &arg)
  {
    return std::vector<std::vector<double>>(
        {{-2 * m * arg[0] + 10 * lam * pow(arg[0], 2) +
          2 * lam * (-2 + pow(arg[0], 2))}});
  };

  VacuumProfile vacuumprofile(dim, TrueVacuum, FalseVacuum, V, dV, Hessian);

  vacuumprofile.itmax = 0;

  vacuumprofile.CalculateProfile();

  REQUIRE(vacuumprofile.status == VacuumProfileStatus::Failed);
}