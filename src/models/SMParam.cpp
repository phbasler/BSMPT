#include <BSMPT/models/SMparam.h>

namespace BSMPT
{
const ISMConstants GetSMConstants()
{
  ISMConstants SM;
  SM.C_Wolfenstein_lambda = 0.22537;
  SM.C_Wolfenstein_A      = 0.814;
  SM.C_Wolfenstein_rho    = 0.117;
  SM.C_Wolfenstein_eta    = 0.353;
  SM.theta12              = std::asin(SM.C_Wolfenstein_lambda);
  SM.theta23 =
      std::asin(SM.C_Wolfenstein_A * std::pow(SM.C_Wolfenstein_lambda, 2));
  SM.delta =
      std::arg(SM.C_Wolfenstein_A * std::pow(SM.C_Wolfenstein_lambda, 3) *
               (SM.C_Wolfenstein_rho + II * SM.C_Wolfenstein_eta));
  SM.theta13 = std::asin(
      std::abs(SM.C_Wolfenstein_A * std::pow(SM.C_Wolfenstein_lambda, 3) *
               (SM.C_Wolfenstein_rho + II * SM.C_Wolfenstein_eta)));
  SM.C_Vud = std::cos(SM.theta12) * std::cos(SM.theta13);
  SM.C_Vus = std::sin(SM.theta12) * std::cos(SM.theta13);
  SM.C_Vub = std::sin(SM.theta13) * std::exp(-SM.delta * II);
  SM.C_Vcd = -std::sin(SM.theta12) * std::cos(SM.theta23) -
             std::cos(SM.theta12) * std::sin(SM.theta23) *
                 std::sin(SM.theta13) * std::exp(II * SM.delta);
  SM.C_Vcs = std::cos(SM.theta12) * std::cos(SM.theta23) -
             std::sin(SM.theta12) * std::sin(SM.theta23) *
                 std::sin(SM.theta13) * std::exp(II * SM.delta);
  SM.C_Vcb = std::sin(SM.theta23) * std::cos(SM.theta13);
  SM.C_Vtd = std::sin(SM.theta12) * std::sin(SM.theta23) -
             std::cos(SM.theta12) * std::cos(SM.theta23) *
                 std::sin(SM.theta13) * std::exp(II * SM.delta);
  SM.C_Vts = -std::cos(SM.theta12) * std::sin(SM.theta23) -
             std::sin(SM.theta12) * std::cos(SM.theta23) *
                 std::sin(SM.theta13) * std::exp(II * SM.delta);
  SM.C_Vtb = std::cos(SM.theta23) * std::cos(SM.theta13);

  SM.C_MassW        = 80.385;
  SM.C_MassZ        = 91.1876;
  SM.C_MassSMHiggs  = 125.09;
  SM.C_MassUp       = 0.1;
  SM.C_MassDown     = 0.1;
  SM.C_MassStrange  = 0.1;
  SM.C_MassTop      = 172.5;
  SM.C_MassCharm    = 1.51;
  SM.C_MassBottom   = 4.92;
  SM.C_MassTau      = 1.77682;
  SM.C_MassMu       = 0.1056583715;
  SM.C_MassElectron = 0.510998928 * std::pow(10.0, -3.0);
  SM.C_GF           = 1.1663787 * 1e-5;

  SM.C_sinsquaredWeinberg =
      1 - (SM.C_MassW * SM.C_MassW) / (SM.C_MassZ * SM.C_MassZ);
  SM.C_vev0 = std::sqrt(1 / std::sqrt(2) * 1 / SM.C_GF);
  SM.C_g    = 2 * SM.C_MassW / SM.C_vev0;
  SM.C_gs   = 2 * std::sqrt(std::pow(SM.C_MassZ, 2) - std::pow(SM.C_MassW, 2)) /
            SM.C_vev0;
  SM.C_SMTriHiggs = 3 * SM.C_MassSMHiggs * SM.C_MassSMHiggs / (SM.C_vev0);

  return SM;
}

const std::complex<double> II(0, 1);

} // namespace BSMPT
