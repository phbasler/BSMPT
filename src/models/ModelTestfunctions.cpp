// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/ModelTestfunctions.h>
#include <BSMPT/utility.h>
#include <iomanip>

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
    std::cerr << "WARNING: The number of labels for the Counterterms does not "
                 "match the number of Counterterms."
              << " If you don't fix this, then your header will not match the "
                 "numerical output in the output file."
              << std::endl;
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
    std::cerr << "WARNING: The number of labels for the VEVs does not match "
                 "the number of VEVs defined in the model."
              << " If you don't fix this, then your header will not match the "
                 "numerical output in the output file."
              << std::endl;
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
    std::cerr
        << "WARNING: The number of labels in addLegendTemp does not match the "
           "number of VEVs + 3."
        << " If you don't fix this, then your header will not match the "
           "numerical output in the output file."
        << " It is expected to be label for the critical temperature, the cirtical VEV, the ratio of VEV and temperature and the labels\
                       for the VEVs."
        << std::endl;
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
    std::cerr
        << "WARNING: The number of labels in addLegendTripleCouplings does "
           "not match the number of calculated Triple Higgs Couplings."
        << " If you don't fix this, then your header will not match the "
           "numerical output in the output file."
        << std::endl;
  }
  return result;
}

TestResults CheckGaugeBosonMasses(const Class_Potential_Origin &point)
{

  std::vector<double> gaugeMassesInput;
  gaugeMassesInput.push_back(0);
  gaugeMassesInput.push_back(pow(C_MassW, 2));
  gaugeMassesInput.push_back(pow(C_MassW, 2));
  gaugeMassesInput.push_back(pow(C_MassZ, 2));
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
    sum += std::abs(GaugeMassCalculated.at(i) - gaugeMassesInput.at(i));
  }
  auto result = sum > 1e-5 ? TestResults::Fail : TestResults::Pass;

  std::string prsize_tline1 = "The SM gauge boson masses squared are : ";
  std::string prsize_tline2 =
      "The calculated gauge boson masses squared are : ";
  auto maxlength = std::max(prsize_tline1.size(), prsize_tline2.size());
  auto addline1  = maxlength - prsize_tline1.size();
  auto addline2  = maxlength - prsize_tline2.size();

  std::cout << prsize_tline1 << std::setw(addline1) << " ";
  for (auto x : gaugeMassesInput)
    std::cout << x << sep;
  std::cout << std::endl;
  std::cout << prsize_tline2 << std::setw(addline2);
  for (auto x : GaugeMassCalculated)
    std::cout << x << sep;
  std::cout << std::endl;

  std::cout << std::endl;

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
  leptonMassesInput.push_back(pow(C_MassElectron, 2));
  leptonMassesInput.push_back(pow(-C_MassElectron, 2));
  leptonMassesInput.push_back(pow(C_MassMu, 2));
  leptonMassesInput.push_back(pow(-C_MassMu, 2));
  leptonMassesInput.push_back(pow(C_MassTau, 2));
  leptonMassesInput.push_back(pow(-C_MassTau, 2));

  quarkMassesInput.push_back(pow(C_MassUp, 2));
  quarkMassesInput.push_back(pow(-C_MassUp, 2));
  quarkMassesInput.push_back(pow(C_MassCharm, 2));
  quarkMassesInput.push_back(pow(-C_MassCharm, 2));
  quarkMassesInput.push_back(pow(C_MassTop, 2));
  quarkMassesInput.push_back(pow(-C_MassTop, 2));

  quarkMassesInput.push_back(pow(C_MassDown, 2));
  quarkMassesInput.push_back(pow(-C_MassDown, 2));
  quarkMassesInput.push_back(pow(C_MassStrange, 2));
  quarkMassesInput.push_back(pow(-C_MassStrange, 2));
  quarkMassesInput.push_back(pow(C_MassBottom, 2));
  quarkMassesInput.push_back(pow(-C_MassBottom, 2));

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

    std::cout << prsize_tline1 << std::setw(addline1) << " ";
    for (auto x : leptonMassesInput)
      std::cout << x << sep;
    std::cout << std::endl;
    std::cout << prsize_tline2 << std::setw(addline2);
    for (auto x : leptonMassCalculated)
      std::cout << x << sep;
    std::cout << std::endl;
    std::cout << std::endl;

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

  std::cout << prsize_tline1 << std::setw(addline1) << " ";
  for (auto x : quarkMassesInput)
    std::cout << x << sep;
  std::cout << std::endl;
  std::cout << prsize_tline2 << std::setw(addline2);
  for (auto x : quarkMassCalculated)
    std::cout << x << sep;
  std::cout << std::endl << std::endl;

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
                                             CheckVector,
                                             start,
                                             WhichMinimizer);
  std::string prsize_tline1 = "The given VEV configuration at tree-level is : ";
  std::string prsize_tline2 =
      "The calculated VEV configuration at tree-level is : ";
  auto maxlength = std::max(prsize_tline1.size(), prsize_tline2.size());
  auto addline1  = maxlength - prsize_tline1.size();
  auto addline2  = maxlength - prsize_tline2.size();

  std::cout << prsize_tline1 << std::setw(addline1) << " ";
  for (auto x : point.get_vevTreeMin())
    std::cout << x << sep;
  std::cout << std::endl;
  std::cout << prsize_tline2 << std::setw(addline2);
  for (auto x : CalculatedHiggsVEV)
    std::cout << x << sep;
  std::cout << std::endl << std::endl;

  double sum{0};
  for (std::size_t i{0}; i < point.get_nVEV(); ++i)
  {
    sum += std::abs(std::abs(CalculatedHiggsVEV.at(i)) -
                    std::abs(point.get_vevTreeMin(i)));
  }
  if (sum > 0.5)
  {
    result = TestResults::Fail;
  }

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
    std::cout << "The given input parameter does not fulfill the tadpole "
                 "relations and is not a minimum of the potential."
              << std::endl
              << "This may happen if all your parameters are read in from an "
                 "input file. Try applying the minimum conditions"
              << " in the set_gen function." << std::endl
              << std::endl;
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

  std::cout << "The higgs masses squared at LO | NLO are : " << std::endl;
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    std::cout << "m_i^2 = " << TreeMass[i] << " | " << NLOMass[i] << std::endl;
  }

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
    std::cout << "You provided a simplified version of the tree-level "
                 "potential but it yields"
              << " different results for the same input compared to the "
                 "explicit calculation. "
              << "Recheck your implementation of the simplified tree-level "
                 "potential."
              << std::endl;
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
    std::cout << "You provided a simplified version of the counterterm "
                 "potential but it yields"
              << " different results for the same input compared to the "
                 "explicit calculation. "
              << "Recheck your implementation of the simplified counterterm "
                 "potential."
              << std::endl;
    result = TestResults::Fail;
  }
  return result;
}

TestResults CheckCKMUnitarity()
{
  using namespace Eigen;
  MatrixXcd VCKM(3, 3);
  VCKM(0, 0) = C_Vud;
  VCKM(0, 1) = C_Vus;
  VCKM(0, 2) = C_Vub;
  VCKM(1, 0) = C_Vcd;
  VCKM(1, 1) = C_Vcs;
  VCKM(1, 2) = C_Vcb;
  VCKM(2, 0) = C_Vtd;
  VCKM(2, 1) = C_Vts;
  VCKM(2, 2) = C_Vtb;

  double ZeroMass = std::pow(10, -5);
  auto norm       = (VCKM.adjoint() * VCKM - MatrixXcd::Identity(3, 3)).norm();
  auto result     = norm > ZeroMass ? TestResults::Fail : TestResults::Pass;
  if (result == TestResults::Fail)
  {

    std::cerr << "Your CKM implementation is not unitary!" << std::endl;
    std::cerr << "Your CKM Matrix V is given by \n" << VCKM << std::endl;
    std::cerr << "with adjoint(V)*V = \n" << VCKM.adjoint() * VCKM << std::endl;
    std::cerr << "The norm deviating from 1 is " << norm << std::endl;
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
  auto result  = TestResults::Pass;

  auto HesseWeinberg = point.WeinbergSecondDerivativeAsMatrixXd();

  auto HesseVCT = point.HessianCT(vevTree);
  for (std::size_t i{0}; i < NHiggs; ++i)
  {
    for (std::size_t j{0}; j < NHiggs; ++j)
    {
      if (std::abs(HesseVCT(i, j) + HesseWeinberg(i, j)) > 1e-5)
      {
        result = TestResults::Fail;
        break;
      }
    }
  }
  return result;
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
} // namespace ModelTests
} // namespace BSMPT
