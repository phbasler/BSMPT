// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * Model file for the SM + real singlet with the tree-level potential
 */

#pragma once

#include <string> // for string
#include <vector> // for vector

#include <BSMPT/models/ClassPotentialOrigin.h>
namespace BSMPT
{
namespace Models
{

/**
 * @brief The Class_RxSM class
 * Lagrangian from https://arxiv.org/pdf/1512.05355 Eq. (11)
 */
class Class_RxSM : public Class_Potential_Origin
{
public:
  Class_RxSM(const ISMConstants &smConstants);
  virtual ~Class_RxSM();

  // Choice of parameters of Lagrangian from https://arxiv.org/pdf/1512.05355 Eq. (11)
  double lambdaS, lambdaHS, vS;

  // Not an input parameter; lambda is fixed via the requirement of having
  // one of the Higgs bosons as the SM one with 125.09 GeV
  double lambda;

  // Not an input parameter; set to the SM value
  double vH;

  // Not input parameters; set through the tadpole equations
  double msq, mSsq;

  double alpha;

  bool UnbrokenSingletPhase;

  double dmsq, dlambda, dmSsq, dlambdaS, dlambdaHS, dT1, dT2, dT3, dT4, dT5;

  int pos_G1, pos_G2, pos_G0, pos_h, pos_H;
  int pos_h_SM, pos_h_H;

  void ReadAndSet(const std::string &linestr,
                  std::vector<double> &par) override;
  std::vector<std::string> addLegendCT() const override;
  std::vector<std::string> addLegendTemp() const override;
  std::vector<std::string> addLegendTripleCouplings() const override;
  std::vector<std::string> addLegendVEV() const override;

  void set_gen(const std::vector<double> &par) override;
  void set_CT_Pot_Par(const std::vector<double> &par) override;
  void write() const override;

  void AdjustRotationMatrix() override;
  void TripleHiggsCouplings() override;
  std::vector<double> calc_CT() const override;

  void SetCurvatureArrays() override;
  bool CalculateDebyeSimplified() override;
  bool CalculateDebyeGaugeSimplified() override;
  double VTreeSimplified(const std::vector<double> &v) const override;
  double VCounterSimplified(const std::vector<double> &v) const override;
  void Debugging(const std::vector<double> &input,
                 std::vector<double> &output) const override;
};

} // namespace Models
} // namespace BSMPT
