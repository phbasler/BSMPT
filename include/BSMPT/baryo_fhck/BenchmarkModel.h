#pragma once

#include <BSMPT/baryo_fhck/TransportModel.h>
#include <vector>

namespace BSMPT
{
namespace Baryo
{
namespace FHCK
{

/**
 * @brief Class based on the benchmark model provided in [2407.13639] Sec.7
 *
 */
class BenchmarkModel : public TransportModel
{
private:
  /**
   *@brief Input parameters of the benchmark model
   *
   */
  double vn, wn, Tn, LAM, Ls;

  /**
   *@brief Top yukawa coupling and weak coupling constant
   */
  const double yt = 1., g = 0.65;

public:
  /**
   * @brief Construct Benchmark model, where everything else is dependent on T
   *
   * @param Tn_in Nucleation temperature
   * @param vw_in Wall velocity
   * @param trunactionscheme_in Truncation scheme to be used
   * @param truncationR_in Truncation constant \f$ R \f$
   */
  BenchmarkModel(
      const double Tn_in,
      const double vw_in,
      const TruncationScheme &truncationscheme_in = TruncationScheme::MinusVw,
      const double &truncationR_in                = 0);

  /**
   * @brief Construct Benchmark model
   *
   * @param vn_in Vev of the Higgs like field
   * @param wn_in Vev of the second scalar
   * @param Tn_in Nucleation temperature
   * @param LAM_in New physics scale
   * @param Lw_in Thickness of the bubble wall
   * @param vw_in Wall velocity
   * @param trunactionscheme_in Truncation scheme to be used
   * @param truncationR_in Truncation constant \f$ R \f$
   */
  BenchmarkModel(
      const double vn_in,
      const double wn_in,
      const double Tn_in,
      const double LAM_in,
      const double Lw_in,
      const double Ls_in,
      const double vw_in,
      const TruncationScheme &truncationscheme_in = TruncationScheme::MinusVw,
      const double &truncationR_in                = 0);

  /**
   * @brief Override function to do nothing
   */
  void Initialize() override;

  /**
   * @brief Override function to do nothing
   *
   * @param zList z distribution where we interpolate the fermion mass
   * @param MakeTopMassPlot Make the plots or not
   */
  void GenerateFermionMass(const std::vector<double> &zList,
                           const bool &MakeTopMassPlot) override;

  /**
   * @brief Vev profile of the Higgs field
   *
   * @param z Distance from the bubble wall 
   * @param deriv Order of derivative with respect to z
   */
  double hvev(const double &z, const int &deriv);

  /**
   * @brief Vev profile of the second scalar field
   *
   * @param z Distance from the bubble wall 
   * @param deriv Order of derivative with respect to z
   */
  double svev(const double &z, const int &deriv);

  /**
   * @brief Get the Fermion Mass object Calculate the fermion mass and its
   * derivatives
   *
   * @param z distance from the bubble wall
   * @param fermion which fermion, 0 = most massive
   * @param x2 \f$ x^2 \f$
   * @param x2prime \f$ x'^2 \f$
   * @param thetaprime \f$ \theta' \f$
   * @param theta2prime \f$ \theta'' \f$
   */
  void GetFermionRatio(const double &z,
                      const size_t &fermion,
                      double &x2,
                      double &x2prime,
                      double &thetaprime,
                      double &theta2prime) override;

  /**
   * @brief Calculate the W boson mass/temperature.
   *
   * @param vev VEV
   * @return double W boson mass/temperature
   */
  double GetWRatio(const double &z) override;

  /**
   * @brief This function calculates the EW breaking VEV from hvev
   *
   * @param z distance from the bubble wall
   */
  double EWSBVEV(const double &z) override;
};
} // namespace FHCK
} // namespace Baryo
} // namespace BSMPT