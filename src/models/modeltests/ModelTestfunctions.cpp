// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/modeltests/ModelTestfunctions.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
#include <array>
#include <iomanip>

namespace
{
BSMPT::ModelTests::TestResults CheckFourPermutation(
    const std::vector<std::vector<std::vector<std::vector<double>>>> &Tensor,
    std::size_t i,
    std::size_t j,
    std::size_t k,
    std::size_t l)
{
  std::array<std::array<std::size_t, 4>, 24> PermutationOrders;
  PermutationOrders[0][0]  = i;
  PermutationOrders[0][1]  = j;
  PermutationOrders[0][2]  = k;
  PermutationOrders[0][3]  = l;
  PermutationOrders[1][0]  = i;
  PermutationOrders[1][1]  = j;
  PermutationOrders[1][2]  = l;
  PermutationOrders[1][3]  = k;
  PermutationOrders[2][0]  = i;
  PermutationOrders[2][1]  = k;
  PermutationOrders[2][2]  = j;
  PermutationOrders[2][3]  = l;
  PermutationOrders[3][0]  = i;
  PermutationOrders[3][1]  = k;
  PermutationOrders[3][2]  = l;
  PermutationOrders[3][3]  = j;
  PermutationOrders[4][0]  = i;
  PermutationOrders[4][1]  = l;
  PermutationOrders[4][2]  = j;
  PermutationOrders[4][3]  = k;
  PermutationOrders[5][0]  = i;
  PermutationOrders[5][1]  = l;
  PermutationOrders[5][2]  = k;
  PermutationOrders[5][3]  = j;
  PermutationOrders[6][0]  = j;
  PermutationOrders[6][1]  = i;
  PermutationOrders[6][2]  = k;
  PermutationOrders[6][3]  = l;
  PermutationOrders[7][0]  = j;
  PermutationOrders[7][1]  = i;
  PermutationOrders[7][2]  = l;
  PermutationOrders[7][3]  = k;
  PermutationOrders[8][0]  = j;
  PermutationOrders[8][1]  = k;
  PermutationOrders[8][2]  = i;
  PermutationOrders[8][3]  = l;
  PermutationOrders[9][0]  = j;
  PermutationOrders[9][1]  = k;
  PermutationOrders[9][2]  = l;
  PermutationOrders[9][3]  = i;
  PermutationOrders[10][0] = j;
  PermutationOrders[10][1] = l;
  PermutationOrders[10][2] = i;
  PermutationOrders[10][3] = k;
  PermutationOrders[11][0] = j;
  PermutationOrders[11][1] = l;
  PermutationOrders[11][2] = k;
  PermutationOrders[11][3] = i;
  PermutationOrders[12][0] = k;
  PermutationOrders[12][1] = i;
  PermutationOrders[12][2] = j;
  PermutationOrders[12][3] = l;
  PermutationOrders[13][0] = k;
  PermutationOrders[13][1] = i;
  PermutationOrders[13][2] = l;
  PermutationOrders[13][3] = j;
  PermutationOrders[14][0] = k;
  PermutationOrders[14][1] = j;
  PermutationOrders[14][2] = i;
  PermutationOrders[14][3] = l;
  PermutationOrders[15][0] = k;
  PermutationOrders[15][1] = j;
  PermutationOrders[15][2] = l;
  PermutationOrders[15][3] = i;
  PermutationOrders[16][0] = k;
  PermutationOrders[16][1] = l;
  PermutationOrders[16][2] = i;
  PermutationOrders[16][3] = j;
  PermutationOrders[17][0] = k;
  PermutationOrders[17][1] = l;
  PermutationOrders[17][2] = j;
  PermutationOrders[17][3] = i;
  PermutationOrders[18][0] = l;
  PermutationOrders[18][1] = i;
  PermutationOrders[18][2] = j;
  PermutationOrders[18][3] = k;
  PermutationOrders[19][0] = l;
  PermutationOrders[19][1] = i;
  PermutationOrders[19][2] = k;
  PermutationOrders[19][3] = j;
  PermutationOrders[20][0] = l;
  PermutationOrders[20][1] = j;
  PermutationOrders[20][2] = i;
  PermutationOrders[20][3] = k;
  PermutationOrders[21][0] = l;
  PermutationOrders[21][1] = j;
  PermutationOrders[21][2] = k;
  PermutationOrders[21][3] = i;
  PermutationOrders[22][0] = l;
  PermutationOrders[22][1] = k;
  PermutationOrders[22][2] = i;
  PermutationOrders[22][3] = j;
  PermutationOrders[23][0] = l;
  PermutationOrders[23][1] = k;
  PermutationOrders[23][2] = j;
  PermutationOrders[23][3] = i;

  for (std::size_t counter{1}; counter < 24; ++counter)
  {
    if (Tensor.at(PermutationOrders[0][0])
            .at(PermutationOrders[0][1])
            .at(PermutationOrders[0][2])
            .at(PermutationOrders[0][3]) !=
        Tensor.at(PermutationOrders[counter][0])
            .at(PermutationOrders[counter][1])
            .at(PermutationOrders[counter][2])
            .at(PermutationOrders[counter][3]))
    {
      return BSMPT::ModelTests::TestResults::Fail;
    }
  }
  return BSMPT::ModelTests::TestResults::Pass;
}
} // namespace

namespace BSMPT
{
namespace ModelTests
{
std::string TestResultsToString(const TestResults &input)
{
  std::string result;
  switch (input)
  {
  case TestResults::Pass:
  {
    result = "Pass";
    break;
  }
  case TestResults::Fail:
  {
    result = "Fail";
    break;
  }
  default:
  {
    result = "Invalid";
    break;
  }
  }
  return result;
}

TestResults CheckNumberOfCTParameters(const Class_Potential_Origin &point)
{
  auto result = TestResults::Fail;
  if (point.get_nParCT() == point.addLegendCT().size())
  {
    result = TestResults::Pass;
  }
  else
  {
    Logger::Write(LoggingLevel::Default,
                  "WARNING: The number of labels for the Counterterms does not "
                  "match the number of Counterterms."
                  " If you don't fix this, then your header will not match the "
                  "numerical output in the output file.");
  }
  return result;
}

TestResults CheckNumberOfVEVLabels(const Class_Potential_Origin &point)
{
  auto result = TestResults::Fail;

  if (point.get_nVEV() == point.addLegendVEV().size())
  {
    result = TestResults::Pass;
  }
  else
  {
    Logger::Write(LoggingLevel::Default,
                  "WARNING: The number of labels for the VEVs does not match "
                  "the number of VEVs defined in the model."
                  " If you don't fix this, then your header will not match the "
                  "numerical output in the output file.");
  }

  return result;
}

TestResults CheckLegendTemp(const Class_Potential_Origin &point)
{
  auto result = (point.addLegendTemp().size() == point.get_nVEV() + 3)
                    ? TestResults::Pass
                    : TestResults::Fail;
  if (result == TestResults::Fail)
  {
    Logger::Write(
        LoggingLevel::Default,
        "WARNING: The number of labels in addLegendTemp does not match the "
        "number of VEVs + 3."
        " If you don't fix this, then your header will not match the "
        "numerical output in the output file."
        " It is expected to be label for the critical temperature, the cirtical VEV, the ratio of VEV and temperature and the labels\
                       for the VEVs.");
  }
  return result;
}

TestResults CheckNumberOfTripleCouplings(const Class_Potential_Origin &point)
{
  std::size_t ExpectedTripleHiggs = 0;
  for (std::size_t i = 0; i < point.get_NHiggs(); i++)
  {
    for (std::size_t j = i; j < point.get_NHiggs(); j++)
    {
      for (std::size_t k = j; k < point.get_NHiggs(); k++)
      {
        ExpectedTripleHiggs++;
      }
    }
  }
  ExpectedTripleHiggs *= 3;
  auto result = (ExpectedTripleHiggs == point.addLegendTripleCouplings().size())
                    ? TestResults::Pass
                    : TestResults::Fail;
  if (result == TestResults::Fail)
  {
    Logger::Write(
        LoggingLevel::Default,
        "WARNING: The number of labels in addLegendTripleCouplings does "
        "not match the number of calculated Triple Higgs Couplings."
        " If you don't fix this, then your header will not match the "
        "numerical output in the output file.");
  }
  return result;
}

TestResults CheckGaugeBosonMasses(const Class_Potential_Origin &point)
{

  std::vector<double> gaugeMassesInput;
  gaugeMassesInput.push_back(0);
  gaugeMassesInput.push_back(pow(point.SMConstants.C_MassW, 2));
  gaugeMassesInput.push_back(pow(point.SMConstants.C_MassW, 2));
  gaugeMassesInput.push_back(pow(point.SMConstants.C_MassZ, 2));
  std::sort(gaugeMassesInput.begin(), gaugeMassesInput.end());
  auto GaugeMassCalculated = point.GaugeMassesSquared(
      point.MinimizeOrderVEV(point.get_vevTreeMin()), 0, 0);
  if (GaugeMassCalculated.size() != gaugeMassesInput.size())
  {
    return TestResults::Fail;
  }
  double sum{0};
  for (std::size_t i{0}; i < GaugeMassCalculated.size(); ++i)
  {
    if (gaugeMassesInput.at(i) == 0)
    {
      sum += std::abs(GaugeMassCalculated.at(i));
    }
    else
    {
      sum += (std::abs(GaugeMassCalculated.at(i) - gaugeMassesInput.at(i))) /
             gaugeMassesInput.at(i);
    }
  }
  auto result = sum > 1e-5 ? TestResults::Fail : TestResults::Pass;

  std::string prsize_tline1 = "The SM gauge boson masses squared are : ";
  std::string prsize_tline2 =
      "The calculated gauge boson masses squared are : ";
  auto maxlength = std::max(prsize_tline1.size(), prsize_tline2.size());
  auto addline1  = maxlength - prsize_tline1.size();
  auto addline2  = maxlength - prsize_tline2.size();

  std::stringstream ss;
  ss << prsize_tline1 << std::setw(addline1) << " ";
  for (auto x : gaugeMassesInput)
    ss << x << sep;
  ss << std::endl;
  ss << prsize_tline2 << std::setw(addline2);
  for (auto x : GaugeMassCalculated)
    ss << x << sep;
  ss << std::endl;

  ss << "The result is " << sum << std::endl;

  Logger::Write(LoggingLevel::ProgDetailed, ss.str());

  return result;
}

std::pair<TestResults, TestResults>
CheckFermionicMasses(const Class_Potential_Origin &point)
{
  auto result = std::make_pair(TestResults::Pass, TestResults::Pass);

  std::vector<double> leptonMassesInput, quarkMassesInput;
  leptonMassesInput.push_back(0);
  leptonMassesInput.push_back(0);
  leptonMassesInput.push_back(0);
  leptonMassesInput.push_back(pow(point.SMConstants.C_MassElectron, 2));
  leptonMassesInput.push_back(pow(-point.SMConstants.C_MassElectron, 2));
  leptonMassesInput.push_back(pow(point.SMConstants.C_MassMu, 2));
  leptonMassesInput.push_back(pow(-point.SMConstants.C_MassMu, 2));
  leptonMassesInput.push_back(pow(point.SMConstants.C_MassTau, 2));
  leptonMassesInput.push_back(pow(-point.SMConstants.C_MassTau, 2));

  quarkMassesInput.push_back(pow(point.SMConstants.C_MassUp, 2));
  quarkMassesInput.push_back(pow(-point.SMConstants.C_MassUp, 2));
  quarkMassesInput.push_back(pow(point.SMConstants.C_MassCharm, 2));
  quarkMassesInput.push_back(pow(-point.SMConstants.C_MassCharm, 2));
  quarkMassesInput.push_back(pow(point.SMConstants.C_MassTop, 2));
  quarkMassesInput.push_back(pow(-point.SMConstants.C_MassTop, 2));

  quarkMassesInput.push_back(pow(point.SMConstants.C_MassDown, 2));
  quarkMassesInput.push_back(pow(-point.SMConstants.C_MassDown, 2));
  quarkMassesInput.push_back(pow(point.SMConstants.C_MassStrange, 2));
  quarkMassesInput.push_back(pow(-point.SMConstants.C_MassStrange, 2));
  quarkMassesInput.push_back(pow(point.SMConstants.C_MassBottom, 2));
  quarkMassesInput.push_back(pow(-point.SMConstants.C_MassBottom, 2));

  if (point.get_NLepton() == 0)
  {
    for (const auto &el : leptonMassesInput)
      quarkMassesInput.push_back(el);
  }
  std::sort(leptonMassesInput.begin(), leptonMassesInput.end());
  std::sort(quarkMassesInput.begin(), quarkMassesInput.end());

  auto leptonMassCalculated = point.LeptonMassesSquared(
      point.MinimizeOrderVEV(point.get_vevTreeMin()), 0);
  auto quarkMassCalculated = point.QuarkMassesSquared(
      point.MinimizeOrderVEV(point.get_vevTreeMin()), 0);

  const double ZeroMass = 1e-5;
  if (point.get_NLepton() != 0)
  {
    std::string prsize_tline1 = "The SM lepton masses squared are : ";
    std::string prsize_tline2 = "The calculated lepton masses squared are : ";
    auto maxlength = std::max(prsize_tline1.size(), prsize_tline2.size());
    auto addline1  = maxlength - prsize_tline1.size();
    auto addline2  = maxlength - prsize_tline2.size();

    std::stringstream ss;
    ss << prsize_tline1 << std::setw(addline1) << " ";
    for (auto x : leptonMassesInput)
      ss << x << sep;
    ss << std::endl;
    ss << prsize_tline2 << std::setw(addline2);
    for (auto x : leptonMassCalculated)
      ss << x << sep;
    ss << std::endl;
    ss << std::endl;
    Logger::Write(LoggingLevel::ProgDetailed, ss.str());

    double sum{0};
    for (std::size_t i{0};
         i < std::min(leptonMassCalculated.size(), leptonMassesInput.size());
         ++i)
    {
      sum += std::abs(leptonMassCalculated.at(i) - leptonMassesInput.at(i));
    }
    if (leptonMassCalculated.size() > leptonMassesInput.size())
    {
      for (std::size_t i{leptonMassesInput.size()};
           i < leptonMassCalculated.size();
           ++i)
      {
        sum += std::abs(leptonMassCalculated.at(i));
      }
    }
    else if (leptonMassesInput.size() > leptonMassCalculated.size())
    {
      for (std::size_t i{leptonMassCalculated.size()};
           i < leptonMassesInput.size();
           ++i)
      {
        sum += std::abs(leptonMassesInput.at(i));
      }
    }
    if (sum > ZeroMass)
    {
      result.first = TestResults::Fail;
    }
  }

  std::string prsize_tline1 = "The SM quark masses squared are : ";
  std::string prsize_tline2 = "The calculated quark masses squared are : ";
  auto maxlength = std::max(prsize_tline1.size(), prsize_tline2.size());
  auto addline1  = maxlength - prsize_tline1.size();
  auto addline2  = maxlength - prsize_tline2.size();

  {
    std::stringstream ss;
    ss << prsize_tline1 << std::setw(addline1) << " ";
    for (auto x : quarkMassesInput)
      ss << x << sep;
    ss << std::endl;
    ss << prsize_tline2 << std::setw(addline2);
    for (auto x : quarkMassCalculated)
      ss << x << sep;
    ss << std::endl << std::endl;
    Logger::Write(LoggingLevel::ProgDetailed, ss.str());
  }

  double sum{0};
  for (std::size_t i{0}; i < quarkMassCalculated.size(); ++i)
  {
    sum += std::abs(quarkMassCalculated.at(i) - quarkMassesInput.at(i));
  }
  if (sum > ZeroMass)
  {
    result.second = TestResults::Fail;
  }

  return result;
}

TestResults CheckTreeLevelMin(const Class_Potential_Origin &point,
                              int WhichMinimizer)
{
  auto result = TestResults::Pass;

  std::vector<double> CalculatedHiggsVEV, CheckVector, start;
  for (const auto &x : point.get_vevTreeMin())
    start.push_back(0.5 * x);

  CalculatedHiggsVEV =
      Minimizer::Minimize_gen_all_tree_level(point.get_Model(),
                                             point.get_parStored(),
                                             point.get_parCTStored(),
                                             point.SMConstants,
                                             CheckVector,
                                             start,
                                             WhichMinimizer);
  std::string prsize_tline1 = "The given VEV configuration at tree-level is : ";
  std::string prsize_tline2 =
      "The calculated VEV configuration at tree-level is : ";
  auto maxlength = std::max(prsize_tline1.size(), prsize_tline2.size());
  auto addline1  = maxlength - prsize_tline1.size();
  auto addline2  = maxlength - prsize_tline2.size();

  std::stringstream ss;
  ss << prsize_tline1 << std::setw(addline1) << " ";
  for (auto x : point.get_vevTreeMin())
    ss << x << sep;
  ss << std::endl;
  ss << prsize_tline2 << std::setw(addline2);
  for (auto x : CalculatedHiggsVEV)
    ss << x << sep;
  ss << std::endl << std::endl;

  double sum{0};
  for (std::size_t i{0}; i < point.get_nVEV(); ++i)
  {
    sum += std::abs(std::abs(CalculatedHiggsVEV.at(i)) -
                    std::abs(point.get_vevTreeMin(i)));
  }
  if (sum > 0.5)
  {
    result = TestResults::Fail;
    ss << "Test failed with difference = " << sum;
  }

  Logger::Write(LoggingLevel::ProgDetailed, ss.str());
  return result;
}

TestResults CheckTadpoleRelations(const Class_Potential_Origin &point)
{
  auto result = TestResults::Pass;

  double SurviveTadpole = 0;
  auto transformedVEV   = point.MinimizeOrderVEV(point.get_vevTreeMin());
  for (std::size_t i = 0; i < point.get_NHiggs(); i++)
    SurviveTadpole += std::abs(point.VTree(transformedVEV, i + 1));

  if (SurviveTadpole > 1e-5)
  {
    Logger::Write(
        LoggingLevel::Default,
        "The given input parameter does not fulfill the tadpole relations and "
        "is not a minimum of the potential.\n This may happen if all your "
        "parameters are read in from an input file. Try applying the minimum "
        "conditions in the set_gen function.\n");
    result = TestResults::Fail;
  }
  return result;
}

TestResults CheckNLOMasses(const Class_Potential_Origin &point)
{
  using namespace Eigen;
  auto NHiggs  = point.get_NHiggs();
  auto vevTree = point.MinimizeOrderVEV(point.get_vevTreeMin());
  auto result  = TestResults::Pass;

  auto HesseWeinberg = point.WeinbergSecondDerivativeAsMatrixXd();

  std::vector<double> WeinbergNabla;
  WeinbergNabla = point.WeinbergFirstDerivative();

  VectorXd NablaWeinberg(NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    NablaWeinberg[i] = WeinbergNabla[i];
  }

  auto NablaVCT   = point.NablaVCT(vevTree);
  auto HesseVCT   = point.HessianCT(vevTree);
  auto MassMatrix = point.HiggsMassMatrix(vevTree);

  SelfAdjointEigenSolver<MatrixXd> esTree(MassMatrix, EigenvaluesOnly);
  SelfAdjointEigenSolver<MatrixXd> esNLO(MassMatrix + HesseVCT + HesseWeinberg,
                                         EigenvaluesOnly);

  std::vector<double> TreeMass, NLOMass;

  const double ZeroMass = 1e-5;
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    if (std::abs(esTree.eigenvalues()[i]) < ZeroMass)
      TreeMass.push_back(0);
    else
      TreeMass.push_back(esTree.eigenvalues()[i]);
    if (std::abs(esNLO.eigenvalues()[i]) < ZeroMass)
      NLOMass.push_back(0);
    else
      NLOMass.push_back(esNLO.eigenvalues()[i]);
  }

  std::stringstream ss;
  ss << "The higgs masses squared at LO | NLO are : " << std::endl;
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    ss << "m_i^2 = " << TreeMass[i] << " | " << NLOMass[i] << std::endl;
  }
  Logger::Write(LoggingLevel::ProgDetailed, ss.str());

  double sum{0.0};
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    double z = std::abs(std::abs(NLOMass[i]) - TreeMass.at(i));
    double n = std::max(std::abs(NLOMass[i]), std::abs(TreeMass[i]));
    if (n != 0) sum += z / n;
  }
  if (sum > 0.5)
  {
    result = TestResults::Fail;
  }
  return result;
}

TestResults CheckVTreeSimplified(const Class_Potential_Origin &point)
{
  auto result = TestResults::Pass;

  std::default_random_engine randGen(0);
  double ValSimplified = 0, ValTensor = 0, PotentialDifference = 0;
  std::vector<double> VevDummyConfigForMinimiser(point.get_nVEV());
  double RNDMax = 500;
  for (int i = 0; i < 10; i++)
  {
    for (std::size_t j = 0; j < point.get_nVEV(); j++)
      VevDummyConfigForMinimiser.at(j) =
          RNDMax *
          (-1 +
           2 * std::generate_canonical<double,
                                       std::numeric_limits<double>::digits>(
                   randGen));
    auto VevDummyConfigForPotential =
        point.MinimizeOrderVEV(VevDummyConfigForMinimiser);
    ValSimplified = point.VTree(VevDummyConfigForPotential);
    ValTensor     = point.VTree(
        VevDummyConfigForPotential, 0, true /* = ForceExplicitCalculation */);
    PotentialDifference += std::abs(ValTensor - ValSimplified) /
                           (std::abs(ValTensor) + std::abs(ValSimplified));
  }
  if (PotentialDifference > 1e-5)
  {
    Logger::Write(LoggingLevel::Default,
                  "You provided a simplified version of the tree-level "
                  "potential but it yields"
                  " different results for the same input compared to the "
                  "explicit calculation. "
                  "Recheck your implementation of the simplified tree-level "
                  "potential.");
    result = TestResults::Fail;
  }

  return result;
}

TestResults CheckVCounterSimplified(const Class_Potential_Origin &point)
{
  auto result = TestResults::Pass;
  std::default_random_engine randGen(0);
  double ValCTSimplified = 0, ValCTTensor = 0, PotentialCTDifference = 0;
  std::vector<double> VevDummyConfigForMinimiser(point.get_nVEV());
  double RNDMax = 500;
  for (int i = 0; i < 10; i++)
  {
    for (std::size_t j = 0; j < point.get_nVEV(); j++)
      VevDummyConfigForMinimiser.at(j) =
          RNDMax *
          (-1 +
           2 * std::generate_canonical<double,
                                       std::numeric_limits<double>::digits>(
                   randGen));
    auto VevDummyConfigForPotential =
        point.MinimizeOrderVEV(VevDummyConfigForMinimiser);
    ValCTSimplified = point.CounterTerm(VevDummyConfigForPotential);
    ValCTTensor     = point.CounterTerm(VevDummyConfigForPotential, 0, true);
    PotentialCTDifference +=
        std::abs(ValCTTensor - ValCTSimplified) /
        (std::abs(ValCTSimplified) + std::abs(ValCTTensor));
  }
  if (PotentialCTDifference > 1)
  {
    Logger::Write(LoggingLevel::Default,
                  "You provided a simplified version of the counterterm "
                  "potential but it yields"
                  " different results for the same input compared to the "
                  "explicit calculation. "
                  "Recheck your implementation of the simplified counterterm "
                  "potential.");
    result = TestResults::Fail;
  }
  return result;
}

TestResults CheckCKMUnitarity(const ISMConstants &SMConstants)
{
  using namespace Eigen;
  MatrixXcd VCKM(3, 3);
  VCKM(0, 0) = SMConstants.C_Vud;
  VCKM(0, 1) = SMConstants.C_Vus;
  VCKM(0, 2) = SMConstants.C_Vub;
  VCKM(1, 0) = SMConstants.C_Vcd;
  VCKM(1, 1) = SMConstants.C_Vcs;
  VCKM(1, 2) = SMConstants.C_Vcb;
  VCKM(2, 0) = SMConstants.C_Vtd;
  VCKM(2, 1) = SMConstants.C_Vts;
  VCKM(2, 2) = SMConstants.C_Vtb;

  double ZeroMass = std::pow(10, -5);
  auto norm       = (VCKM.adjoint() * VCKM - MatrixXcd::Identity(3, 3)).norm();
  auto result     = norm > ZeroMass ? TestResults::Fail : TestResults::Pass;
  if (result == TestResults::Fail)
  {

    Logger::Write(LoggingLevel::Default,
                  "Your CKM implementation is not unitary!");
    std::stringstream ss;
    ss << "Your CKM Matrix V is given by \n" << VCKM << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
    ss.clear();
    ss << "with adjoint(V)*V = \n" << VCKM.adjoint() * VCKM << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
    Logger::Write(LoggingLevel::Default,
                  "The norm deviating from 1 is " + std::to_string(norm));
  }
  return result;
}

TestResults
CheckCTConditionsFirstDerivative(const Class_Potential_Origin &point)
{
  using namespace Eigen;
  auto NHiggs  = point.get_NHiggs();
  auto vevTree = point.MinimizeOrderVEV(point.get_vevTreeMin());
  auto result  = TestResults::Pass;

  std::vector<double> WeinbergNabla;
  WeinbergNabla = point.WeinbergFirstDerivative();

  auto NablaVCT = point.NablaVCT(vevTree);

  for (std::size_t i{0}; i < NHiggs; ++i)
  {
    if (std::abs(WeinbergNabla.at(i) + NablaVCT(i)) > 1e-5)
    {
      result = TestResults::Fail;
      break;
    }
  }
  return result;
}

TestResults
CheckCTConditionsSecondDerivative(const Class_Potential_Origin &point)
{
  using namespace Eigen;
  auto NHiggs  = point.get_NHiggs();
  auto vevTree = point.MinimizeOrderVEV(point.get_vevTreeMin());

  auto HesseWeinberg = point.WeinbergSecondDerivativeAsMatrixXd();

  auto HesseVCT = point.HessianCT(vevTree);
  for (std::size_t i{0}; i < NHiggs; ++i)
  {
    for (std::size_t j{0}; j < NHiggs; ++j)
    {
      if (std::abs(HesseVCT(i, j) + HesseWeinberg(i, j)) > 1e-5)
      {
        std::stringstream ss;
        ss << "Failed at << (" << i << "," << j << ")"
           << " with CT = " << HesseVCT(i, j)
           << "; CW = " << HesseWeinberg(i, j)
           << "; CT+CW = " << HesseVCT(i, j) + HesseWeinberg(i, j) << std::endl;
        Logger::Write(LoggingLevel::ProgDetailed, ss.str());
        return TestResults::Fail;
      }
    }
  }
  return TestResults::Pass;
}

TestResults CheckCTIdentities(const Class_Potential_Origin &point)
{
  auto result     = TestResults::Pass;
  auto identities = point.GetCTIdentities();
  for (const auto &el : identities)
  {
    if (std::abs(el) > 1e-5)
    {
      result = TestResults::Fail;
      break;
    }
  }
  return result;
}

TestResults CheckCTNumber(const Class_Potential_Origin &point)
{
  return (point.calc_CT().size() == point.get_nParCT()) ? TestResults::Pass
                                                        : TestResults::Fail;
}

TestResults
CheckSymmetricTensorScalarSecond(const std::vector<std::vector<double>> &Tensor)
{
  if (Tensor.empty())
  {
    return TestResults::Pass;
  }
  const auto nRows = Tensor.size();
  const auto nCols = Tensor.at(0).size();

  if (nRows != nCols)
  {
    return TestResults::Fail;
  }

  for (const auto &rows : Tensor)
  {
    if (rows.size() != nCols)
    {
      return TestResults::Fail;
    }
  }

  for (std::size_t i{0}; i < nRows; ++i)
  {
    for (std::size_t j{i}; j < nCols; j++)
    {
      if (Tensor.at(i).at(j) != Tensor.at(j).at(i))
      {
        return TestResults::Fail;
      }
    }
  }

  return TestResults::Pass;
}
TestResults CheckSymmetricTensorScalarThird(
    const std::vector<std::vector<std::vector<double>>> &Tensor)
{
  if (Tensor.empty())
  {
    return TestResults::Pass;
  }

  const auto nDim = Tensor.size();
  for (const auto &first : Tensor)
  {
    if (first.size() != nDim)
    {
      return TestResults::Fail;
    }
    for (const auto &sec : first)
    {
      if (sec.size() != nDim)
      {
        return TestResults::Fail;
      }
    }
  }

  for (std::size_t i{0}; i < nDim; ++i)
  {
    for (std::size_t j{0}; j < nDim; ++j)
    {
      for (std::size_t k{0}; k < nDim; ++k)
      {
        if (Tensor.at(i).at(j).at(k) != Tensor.at(i).at(k).at(j))
        {
          return TestResults::Fail;
        }

        if (Tensor.at(i).at(j).at(k) != Tensor.at(j).at(i).at(k))
        {
          return TestResults::Fail;
        }

        if (Tensor.at(i).at(j).at(k) != Tensor.at(j).at(k).at(i))
        {
          return TestResults::Fail;
        }

        if (Tensor.at(i).at(j).at(k) != Tensor.at(k).at(i).at(j))
        {
          return TestResults::Fail;
        }

        if (Tensor.at(i).at(j).at(k) != Tensor.at(k).at(j).at(i))
        {
          return TestResults::Fail;
        }
      }
    }
  }

  return TestResults::Pass;
}
TestResults CheckSymmetricTensorScalarFourth(
    const std::vector<std::vector<std::vector<std::vector<double>>>> &Tensor)
{
  if (Tensor.empty())
  {
    return TestResults::Pass;
  }

  const auto nDim = Tensor.size();
  for (const auto &first : Tensor)
  {
    if (nDim != first.size())
    {
      return TestResults::Fail;
    }
    for (const auto &sec : first)
    {
      if (sec.size() != nDim)
      {
        return TestResults::Fail;
      }
      for (const auto &third : sec)
      {
        if (third.size() != nDim)
        {
          return TestResults::Fail;
        }
      }
    }
  }

  for (std::size_t i{0}; i < nDim; ++i)
  {
    for (std::size_t j{0}; j < nDim; ++j)
    {
      for (std::size_t k{0}; k < nDim; ++k)
      {
        for (std::size_t l{0}; l < nDim; ++l)
        {
          if (CheckFourPermutation(Tensor, i, j, k, l) == TestResults::Fail)
          {
            return TestResults::Fail;
          }
        }
      }
    }
  }
  return TestResults::Pass;
}

TestResults CheckSymmetricTensorLeptonsThird(
    const std::vector<std::vector<std::vector<std::complex<double>>>> &Tensor)
{
  if (Tensor.empty())
  {
    return TestResults::Pass;
  }

  const auto nDim = Tensor.size();
  for (const auto &rows : Tensor)
  {
    if (nDim != rows.size())
    {
      return TestResults::Fail;
    }
  }

  const auto nHiggs = Tensor.at(0).at(0).size();

  for (std::size_t i{0}; i < nDim; ++i)
  {
    for (std::size_t j{0}; j < nDim; ++j)
    {
      for (std::size_t k{0}; k < nHiggs; ++k)
      {
        if (Tensor.at(i).at(j).at(k) != Tensor.at(j).at(i).at(k))
        {
          return TestResults::Fail;
        }
      }
    }
  }

  return TestResults::Pass;
}
TestResults CheckSymmetricTensorQuarksThird(
    const std::vector<std::vector<std::vector<std::complex<double>>>> &Tensor)
{
  return CheckSymmetricTensorLeptonsThird(Tensor);
}

TestResults CheckSymmetricTensorLeptons(
    const std::vector<std::vector<std::complex<double>>> &Tensor)
{
  if (Tensor.empty())
  {
    return TestResults::Pass;
  }

  const auto nDim = Tensor.size();
  for (const auto &rows : Tensor)
  {
    if (nDim != rows.size())
    {
      return TestResults::Fail;
    }
  }

  for (std::size_t i{0}; i < nDim; ++i)
  {
    for (std::size_t j{0}; j < nDim; ++j)
    {
      if (Tensor.at(i).at(j) != Tensor.at(j).at(i))
      {
        return TestResults::Fail;
      }
    }
  }

  return TestResults::Pass;
}
TestResults CheckSymmetricTensorQuarks(
    const std::vector<std::vector<std::complex<double>>> &Tensor)
{
  return CheckSymmetricTensorLeptons(Tensor);
}
TestResults CheckSymmetricTensorGauge(
    const std::vector<std::vector<std::vector<std::vector<double>>>> &Tensor)
{
  if (Tensor.empty())
  {
    return TestResults::Pass;
  }

  const auto nGauge = Tensor.size();
  for (const auto &rows : Tensor)
  {
    if (rows.size() != nGauge)
    {
      return TestResults::Fail;
    }
  }

  const auto nHiggs = Tensor.at(0).at(0).size();

  if (nHiggs == 0)
  {
    return TestResults::Fail;
  }

  for (const auto &nRowsGauge : Tensor)
  {
    for (const auto &nRowsHiggs : nRowsGauge)
    {
      if (nRowsHiggs.size() != nHiggs)
      {
        return TestResults::Fail;
      }
      for (const auto &nColHiggs : nRowsHiggs)
      {
        if (nColHiggs.size() != nHiggs)
        {
          return TestResults::Fail;
        }
      }
    }
  }

  for (std::size_t a{0}; a < nGauge; ++a)
  {
    for (std::size_t b{0}; b < nGauge; ++b)
    {
      for (std::size_t i{0}; i < nHiggs; ++i)
      {
        for (std::size_t j{0}; j < nHiggs; ++j)
        {
          if (Tensor.at(a).at(b).at(i).at(j) != Tensor.at(a).at(b).at(j).at(i))
          {
            return TestResults::Fail;
          }

          if (Tensor.at(a).at(b).at(i).at(j) != Tensor.at(b).at(a).at(j).at(i))
          {
            return TestResults::Fail;
          }

          if (Tensor.at(a).at(b).at(i).at(j) != Tensor.at(b).at(a).at(i).at(j))
          {
            return TestResults::Fail;
          }
        }
      }
    }
  }

  return TestResults::Pass;
}

void CheckImplementation(const Class_Potential_Origin &point,
                         const int &WhichMinimizer)
{
  using std::pow;

  const std::string Pass = "Pass";
  const std::string Fail = "Fail";

  std::stringstream ss_model;
  ss_model << "The tested Model is the " << point.get_Model();
  Logger::Write(LoggingLevel::Default, ss_model.str());

  std::vector<std::string> TestNames, TestResults;

  TestNames.push_back("CT number/label match");
  TestResults.push_back(ModelTests::TestResultsToString(
      ModelTests::CheckNumberOfCTParameters(point)));

  TestNames.push_back("VEV number/label match");
  TestResults.push_back(ModelTests::TestResultsToString(
      ModelTests::CheckNumberOfVEVLabels(point)));

  TestNames.push_back("addLegendTemp number/label match");
  TestResults.push_back(
      ModelTests::TestResultsToString(ModelTests::CheckLegendTemp(point)));

  TestNames.push_back("addLegendTripleCouplings number/label match");
  TestResults.push_back(ModelTests::TestResultsToString(
      ModelTests::CheckNumberOfTripleCouplings(point)));

  TestNames.push_back("CKM matrix unitarity");
  TestResults.push_back(ModelTests::TestResultsToString(
      ModelTests::CheckCKMUnitarity(point.SMConstants)));

  Logger::Write(LoggingLevel::Default,
                "This function calculates the masses of the gauge bosons, "
                "fermions and Higgs boson and compares them "
                "with the parameters defined in SMparam.h.");

  TestNames.push_back("Matching gauge boson masses with SMparam.h");
  TestResults.push_back(ModelTests::TestResultsToString(
      ModelTests::CheckGaugeBosonMasses(point)));

  auto FermionCheck = ModelTests::CheckFermionicMasses(point);
  TestNames.push_back("Matching lepton masses with SMparam.h");
  TestResults.push_back(ModelTests::TestResultsToString(FermionCheck.first));
  TestNames.push_back("Matching quark masses with SMparam.h");
  TestResults.push_back(ModelTests::TestResultsToString(FermionCheck.second));

  TestNames.push_back("Correct EW Minimum");
  TestResults.push_back(ModelTests::TestResultsToString(
      ModelTests::CheckTreeLevelMin(point, WhichMinimizer)));

  TestNames.push_back("Tadpole relations fullfilled");
  TestResults.push_back(ModelTests::TestResultsToString(
      ModelTests::CheckTadpoleRelations(point)));

  TestNames.push_back("Matching Masses between NLO and tree-level");
  TestResults.push_back(
      ModelTests::TestResultsToString(ModelTests::CheckNLOMasses(point)));

  if (point.UseVTreeSimplified)
  {
    TestNames.push_back("Checking VTreeSimplified");
    TestResults.push_back(ModelTests::TestResultsToString(
        ModelTests::CheckVTreeSimplified(point)));
  }

  if (point.UseVCounterSimplified)
  {
    TestNames.push_back("Checking VCounterSimplified");
    TestResults.push_back(ModelTests::TestResultsToString(
        ModelTests::CheckVCounterSimplified(point)));
  }

  TestNames.push_back("Checking second derivative of CW+CT");
  TestResults.push_back(ModelTests::TestResultsToString(
      ModelTests::CheckCTConditionsSecondDerivative(point)));

  TestNames.push_back("Checking first derivative of CW+CT");
  TestResults.push_back(ModelTests::TestResultsToString(
      ModelTests::CheckCTConditionsFirstDerivative(point)));

  TestNames.push_back("Check symmetric properties of the scalar tensor Lij");
  TestResults.push_back(ModelTests::TestResultsToString(
      ModelTests::CheckSymmetricTensorScalarSecond(
          point.Get_Curvature_Higgs_L2())));

  TestNames.push_back("Check symmetric properties of the scalar tensor Lijk");
  TestResults.push_back(ModelTests::TestResultsToString(
      ModelTests::CheckSymmetricTensorScalarThird(
          point.Get_Curvature_Higgs_L3())));

  TestNames.push_back("Check symmetric properties of the scalar tensor Lijkl");
  TestResults.push_back(ModelTests::TestResultsToString(
      ModelTests::CheckSymmetricTensorScalarFourth(
          point.Get_Curvature_Higgs_L4())));

  TestNames.push_back(
      "Check symmetric properties of the gauge tensor in the C2HDM");
  TestResults.push_back(ModelTests::TestResultsToString(
      ModelTests::CheckSymmetricTensorGauge(point.Get_Curvature_Gauge_G2H2())));

  TestNames.push_back(
      "Check symmetric properties of the Lepton tensor in the C2HDM");
  TestResults.push_back(ModelTests::TestResultsToString(
      ModelTests::CheckSymmetricTensorLeptonsThird(
          point.Get_Curvature_Lepton_F2H1())));

  TestNames.push_back(
      "Check symmetric properties of the mass Lepton tensor in the C2HDM");
  TestResults.push_back(
      ModelTests::TestResultsToString(ModelTests::CheckSymmetricTensorLeptons(
          point.Get_Curvature_Lepton_F2())));

  TestNames.push_back(
      "Check symmetric properties of the Quark tensor in the C2HDM");
  TestResults.push_back(ModelTests::TestResultsToString(
      ModelTests::CheckSymmetricTensorQuarksThird(
          point.Get_Curvature_Quark_F2H1())));

  TestNames.push_back(
      "Check symmetric properties of the mass Quark tensor in the C2HDM");
  TestResults.push_back(ModelTests::TestResultsToString(
      ModelTests::CheckSymmetricTensorQuarks(point.Get_Curvature_Quark_F2())));

  if (TestNames.size() != TestResults.size())
  {
    std::stringstream ss;
    ss << "TestNames : " << std::endl << TestNames << std::endl;
    ss << "TestResults : " << std::endl << TestResults << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
    std::string errmsg{
        "Mismatch between the number of labels for the tests and the results."};
    errmsg += "TestNames.size() = " + std::to_string(TestNames.size());
    errmsg += " TestResults.size() = " + std::to_string(TestResults.size());
    throw std::runtime_error(errmsg);
  }

  std::size_t Passes{0};
  for (const auto &el : TestResults)
  {
    if (el == Pass) Passes++;
  }

  std::stringstream ss;
  ss << "\nTEST | Pass/Fail\n================================\n" << std::endl;
  auto maxsize = TestNames.at(0).size();
  for (const auto &el : TestNames)
  {
    if (el.size() > maxsize) maxsize = el.size();
  }
  maxsize += 5;
  for (std::size_t i{0}; i < TestNames.size(); ++i)
  {
    ss << std::setw(maxsize) << std::left << TestNames.at(i) << "| "
       << TestResults.at(i) << std::endl;
  }

  ss << Passes << " tests out of " << TestResults.size() << " passed.\n"
     << std::endl;
  if (Passes != TestResults.size())
  {
    ss << TestResults.size() - Passes
       << " tests failed. Please check and try again." << std::endl;
  }
  else
  {
    ss << "You're good to go!\n" << std::endl;
  }
  Logger::Write(LoggingLevel::Default, ss.str());
}

} // namespace ModelTests
} // namespace BSMPT
