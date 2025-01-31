// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
using Approx = Catch::Approx;

#include <BSMPT/baryo_calculation/CalculateEtaInterface.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/modeltests/ModelTestfunctions.h>

#include "C2HDM.h"

#include <fstream>
#include <optional>

const std::vector<double> example_point_C2HDM{/* lambda_1 = */ 3.29771,
                                              /* lambda_2 = */ 0.274365,
                                              /* lambda_3 = */ 4.71019,
                                              /* lambda_4 = */ -2.23056,
                                              /* Re(lambda_5) = */ -2.43487,
                                              /* Im(lambda_5) = */ 0.124948,
                                              /* Re(m_{12}^2) = */ 2706.86,
                                              /* tan(beta) = */ 4.64487,
                                              /* Yukawa Type = */ 1};

namespace
{
void writeBaryoConfigFile(
    const std::string &filename,
    const std::vector<std::optional<std::string>> &includeStrings,
    std::optional<int> massiveConfig)
{
  std::ofstream file(filename);
  if (includeStrings.at(0).has_value())
  {
    file << "VIA Ansatz only including the top quark in the transport "
            "equations "
         << "\n"
         << "Include: " << includeStrings.at(0).value() << "\n";
  }

  if (includeStrings.at(1).has_value())
  {
    file << "VIA Ansatz including the top and bottom quark in the transport "
            "equations"
         << "\n"
         << "Include: " << includeStrings.at(1).value() << "\n";
  }

  if (includeStrings.at(2).has_value())
  {
    file << "VIA Ansatz including the top and bottom quark and the tau lepton "
            "in "
            "the transport equations "
         << "\n"
         << "Include: " << includeStrings.at(2).value() << "\n";
  }

  if (massiveConfig.has_value())
  {
    file << "VIA Ansatz treating the bottom quark massive " << "\n"
         << "Massive: " << (massiveConfig.value() == 1 ? " yes " : " no ")
         << "\n";
  }

  if (includeStrings.at(3).has_value())
  {
    file << "FH Ansatz with the plasma velocities " << "\n"
         << "Include: " << includeStrings.at(3).value() << "\n";
  }

  if (includeStrings.at(4).has_value())
  {
    file << "FH Ansatz with the plasma velocities replaced through the second "
            "derivatives "
         << "\n"
         << "Include: " << includeStrings.at(4).value() << "\n";
  }
}

void CheckFileForConfig(const std::string &filename,
                        const std::pair<std::vector<bool>, int> &expectedConfig,
                        const std::string &expectedException)
{
  std::string exceptionMessage;
  std::pair<std::vector<bool>, int> config;
  try
  {
    BSMPT::Baryo::CalculateEtaInterface etaInterface(
        std::vector<bool>(5, true), 1, BSMPT::GetSMConstants());
    config = etaInterface.ReadConfigFile(filename);
  }
  catch (std::exception &e)
  {
    exceptionMessage = e.what();
  }

  if (expectedException != std::string())
  {
    REQUIRE(exceptionMessage == expectedException);
  }
  else
  {
    REQUIRE(expectedConfig.first == config.first);
    REQUIRE(expectedConfig.second == config.second);
  }
}

template <typename T>
void CompareValues(T expected, T result, double epsilon, double threshold)
{
  if (std::abs(expected) >= threshold)
  {
    REQUIRE(result == Approx(expected).epsilon(epsilon));
  }
  else
  {
    REQUIRE(std::abs(result) <= threshold);
  }
}

} // namespace

const Compare_C2HDM Expected;

TEST_CASE("Checking EWBG for C2HDM", "[c2hdm]")
{

  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM, SMConstants);
  modelPointer->initModel(example_point_C2HDM);

  const auto WhichMin = Minimizer::WhichMinimizerDefault;

  const auto EWPT = Expected.EWPTPerSetting.at(WhichMin);

  std::vector<double> vevsymmetricSolution, checksym, startpoint;
  for (const auto &el : EWPT.EWMinimum)
    startpoint.push_back(0.5 * el);
  vevsymmetricSolution = Minimizer::Minimize_gen_all(
      modelPointer, EWPT.Tc + 1, checksym, startpoint, WhichMin, true);
  auto min_expected = Expected.vevSymmetricPerSetting.at(WhichMin);
  REQUIRE(vevsymmetricSolution.size() == min_expected.size());
  for (std::size_t i{0}; i < vevsymmetricSolution.size(); ++i)
  {
    auto res      = std::abs(vevsymmetricSolution.at(i));
    auto expected = std::abs(min_expected.at(i));
    UNSCOPED_INFO("Current Option for Minimizer:\t" << WhichMin);
    UNSCOPED_INFO("This is the position:"
                  << i << "\tFound solution =" << vevsymmetricSolution.at(i)
                  << "\tExpected solution = " << min_expected.at(i));
    CompareValues(expected, res, 1e-4, 1e-4);
  }

  // Call: Calculation of eta in the different implemented approaches

  auto config =
      std::pair<std::vector<bool>, int>{std::vector<bool>(5, true), 1};
  Baryo::CalculateEtaInterface EtaInterface(config, SMConstants);
  const double testVW = Expected.testVW;

  REQUIRE(EWPT.vc > EWPT.Tc);
  auto eta      = EtaInterface.CalcEta(testVW,
                                  EWPT.EWMinimum,
                                  vevsymmetricSolution,
                                  EWPT.Tc,
                                  modelPointer,
                                  WhichMin);
  const auto LW = EtaInterface.getLW();
  REQUIRE(LW == Approx(Expected.LWPerSetting.at(WhichMin)).epsilon(1e-4));

  auto expectedEta = Expected.etaPerSetting.at(WhichMin);
  REQUIRE(eta.size() == expectedEta.size());
  UNSCOPED_INFO("Check eta components");
  const double etaThreshold = 1e-15;
  for (std::size_t i{0}; i < vevsymmetricSolution.size(); ++i)
  {
    auto res      = std::abs(eta.at(i));
    auto expected = std::abs(expectedEta.at(i));

    UNSCOPED_INFO("Current Option for Minimizer:\t" << WhichMin);
    UNSCOPED_INFO("This is the position:"
                  << i << "\tFound solution =" << eta.at(i)
                  << "\tExpected solution = " << expectedEta.at(i));
    CompareValues(expected, res, 1e-2, etaThreshold);
  }
}

TEST_CASE("Checking ReadConfig", "[baryo]")
{
  {
    std::string filename = "Readconfig_testA.txt";
    std::pair<std::vector<bool>, int> expectedConfig;
    expectedConfig.first  = std::vector<bool>{true, false, true, true, false};
    expectedConfig.second = 1;

    std::vector<std::optional<std::string>> includeValues{
        "yes", "no", "yes", "yes", "no"};
    std::optional<int> massiveConfig = 1;
    writeBaryoConfigFile(filename, includeValues, massiveConfig);
    CheckFileForConfig(filename, expectedConfig, std::string());
  }

  {
    std::string filename = "Readconfig_testB.txt";
    std::pair<std::vector<bool>, int> expectedConfig;
    expectedConfig.first  = std::vector<bool>{};
    expectedConfig.second = 0;

    std::vector<std::optional<std::string>> includeValues{
        "yes", "no", "yes", "yes", "someWrongvalue"};
    std::optional<int> massiveConfig = 1;
    std::string expectedException =
        "One of the settings for the EWBG config file is not set to yes or no. "
        "Please change this.";
    writeBaryoConfigFile(filename, includeValues, massiveConfig);
    CheckFileForConfig(filename, expectedConfig, expectedException);
  }
}
