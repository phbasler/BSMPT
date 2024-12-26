// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <complex>
#include <string>
#include <vector>

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/SMparam.h>

namespace BSMPT
{
class Class_Potential_Origin;
namespace ModelTests
{

enum class TestResults
{
  Pass,
  Fail
};

std::string TestResultsToString(const TestResults &input);

TestResults CheckNumberOfCTParameters(const Class_Potential_Origin &point);
TestResults CheckNumberOfVEVLabels(const Class_Potential_Origin &point);
TestResults CheckLegendTemp(const Class_Potential_Origin &point);
TestResults CheckNumberOfTripleCouplings(const Class_Potential_Origin &point);
TestResults CheckGaugeBosonMasses(const Class_Potential_Origin &point);
std::pair<TestResults, TestResults>
CheckFermionicMasses(const Class_Potential_Origin &point);
TestResults CheckTreeLevelMin(const Class_Potential_Origin &point,
                              int WhichMinimizer);

TestResults CheckTadpoleRelations(const Class_Potential_Origin &point);
TestResults CheckNLOMasses(const Class_Potential_Origin &point);
TestResults CheckVTreeSimplified(const Class_Potential_Origin &point);
TestResults CheckVCounterSimplified(const Class_Potential_Origin &point);
TestResults
CheckCTConditionsFirstDerivative(const Class_Potential_Origin &point);
TestResults
CheckCTConditionsSecondDerivative(const Class_Potential_Origin &point);
TestResults CheckCTIdentities(const Class_Potential_Origin &point);
TestResults CheckCTNumber(const Class_Potential_Origin &point);

TestResults CheckCKMUnitarity(const ISMConstants &SMConstants);
TestResults CheckSymmetricTensorScalarSecond(
    const std::vector<std::vector<double>> &Tensor);
TestResults CheckSymmetricTensorScalarThird(
    const std::vector<std::vector<std::vector<double>>> &Tensor);
TestResults CheckSymmetricTensorScalarFourth(
    const std::vector<std::vector<std::vector<std::vector<double>>>> &Tensor);
TestResults CheckSymmetricTensorLeptonsThird(
    const std::vector<std::vector<std::vector<std::complex<double>>>> &Tensor);
TestResults CheckSymmetricTensorQuarksThird(
    const std::vector<std::vector<std::vector<std::complex<double>>>> &Tensor);
TestResults CheckSymmetricTensorLeptons(
    const std::vector<std::vector<std::complex<double>>> &Tensor);
TestResults CheckSymmetricTensorQuarks(
    const std::vector<std::vector<std::complex<double>>> &Tensor);
TestResults CheckSymmetricTensorGauge(
    const std::vector<std::vector<std::vector<std::vector<double>>>> &Tensor);

/**
 * Checks if the tensors are correctly implemented. For this the fermion,
 * quark and gauge boson masses are calculated and printed next to the values
 * defined in SMparah.h
 */
void CheckImplementation(
    const Class_Potential_Origin &point,
    const int &WhichMinimizer = Minimizer::WhichMinimizerDefault);

} // namespace ModelTests
} // namespace BSMPT
