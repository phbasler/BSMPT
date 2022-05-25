// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
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
 * @brief The C2HDMSympy class
 * Template for implementing a new model
 */
class C2HDMSympy : public Class_Potential_Origin
{
public:
  C2HDMSympy();
  virtual ~C2HDMSympy();

  // Add here your parameters for the Lagrangian as well as for the counterterm
  // potential Add here your variables in which you will save the Debye
  // correction factors

  double Imlambda5, Imlambda6, Imlambda7, Imm12sq;
  double Relambda5, Relambda6, Relambda7, lambda1, lambda2, lambda3, lambda4;
  double m11sq, m22sq, Rem12sq;

  double dImlambda5, dImlambda6, dImlambda7, dImm12sq;
  double dRelambda5, dRelambda6, dRelambda7, dlambda1, dlambda2, dlambda3,
      dlambda4;
  double dm11sq, dm22sq, dRem12sq;

  double dT1, dT2, dT3, dT4, dT5, dT6, dT7, dT8;

  enum class THDMType
  {
    Invalid        = 0,
    TypeI          = 1,
    TypeII         = 2,
    LeptonSpecific = 3,
    Flipped        = 4
  };
  THDMType Type;

  double CTempC1 = 0, CTempC2 = 0, CTempCS = 0;
  double R_Hh_1 = 0, R_Hh_2 = 0, R_Hh_3 = 0, R_Hl_1 = 0, R_Hl_2 = 0, R_Hl_3 = 0,
         R_Hsm_1 = 0, R_Hsm_2 = 0, R_Hsm_3 = 0;

  double TanBeta = 0, C_CosBeta = 0, C_SinBeta = 0, C_CosBetaSquared = 0,
         C_SinBetaSquared = 0;

  void ReadAndSet(const std::string &linestr,
                  std::vector<double> &par) override;
  std::vector<std::string> addLegendCT() const override;
  std::vector<std::string> addLegendTemp() const override;
  std::vector<std::string> addLegendTripleCouplings() const override;
  std::vector<std::string> addLegendVEV() const override;

  void set_gen(const std::vector<double> &par) override;
  void set_CT_Pot_Par(const std::vector<double> &par) override;
  void write() const override;

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
