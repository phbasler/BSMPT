// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/IterativeLinearSolvers"
#include <BSMPT/models/SMparam.h> // for C_vev0, C_MassTop, C_g
#include <algorithm>              // for max, copy
#include <iostream>               // for operator<<, endl, basic_o...
#include <memory>                 // for allocator_traits<>::value...
#include <stddef.h>               // for std::size_t

#include <BSMPT/models/C2HDMUsingSympy.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
using namespace Eigen;

/**
 * @file
 * Template for adding a new model class
 */

namespace BSMPT
{
namespace Models
{
/**
 * Here you have to adjust NNeutralHiggs, NChargedHiggs, nPar (number of
 * Lagrangian parameters AFTER using the tadpole conditions), nParCT (number of
 * counterterms) as well as nVEV (number of VEVs for minimization)
 */
C2HDMSympy::C2HDMSympy()
{
  Model = ModelID::ModelIDs::C2HDMSympy; // global int constant which will be
                                         // used to tell the program which model
                                         // is called
  NNeutralHiggs = 4; // number of neutral Higgs bosons at T = 0
  NChargedHiggs = 4; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar   = 9;  // number of parameters in the tree-Level Lagrangian
  nParCT = 22; // number of parameters in the counterterm potential

  nVEV = 4; // number of VEVs to minimize the potential

  NHiggs = NNeutralHiggs + NChargedHiggs;

  VevOrder.resize(nVEV);
  // Here you have to tell which scalar field gets which VEV.
  VevOrder.at(0) = 2;
  VevOrder.at(1) = 4;
  VevOrder.at(2) = 6;
  VevOrder.at(3) = 7;

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = true;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = true;
}

C2HDMSympy::~C2HDMSympy()
{
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given input file
 */
std::vector<std::string> C2HDMSympy::addLegendCT() const
{
  std::vector<std::string> labels(nParCT);
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * input file
 */
std::vector<std::string> C2HDMSympy::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c");
  labels.push_back("omega_c");
  labels.push_back("omega_c/T_c");
  labels.push_back("omega_CB(T_c)");
  labels.push_back("omega_1(T_c)");
  labels.push_back("omega_2(T_c)");
  labels.push_back("omega_CP(T_c)");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * Higgs couplings. Use this to complement the legend of the given input file
 *
 */
std::vector<std::string> C2HDMSympy::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;
  particles.resize(NHiggs);
  // here you have to define the particle names in the vector particles

  particles[0] = ("G^+");
  particles[1] = ("G^-");
  particles[2] = ("H^+");
  particles[3] = ("H^-");
  particles[4] = ("G^0");
  particles[5] = ("h_1");
  particles[6] = ("h_2");
  particles[7] = ("h_3");

  std::string out = "Tree_";
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = i; j < NHiggs; j++)
    {
      for (std::size_t k = j; k < NHiggs; k++)
      {
        labels.push_back("Tree_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
        labels.push_back("CT_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
        labels.push_back("CW_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
      }
    }
  }

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs.
 * Use this to complement the legend of the given input file
 */
std::vector<std::string> C2HDMSympy::addLegendVEV() const
{
  std::vector<std::string> labels;
  // out = "Your VEV order";
  labels.push_back("omega_CB");
  labels.push_back("omega_1");
  labels.push_back("omega_2");
  labels.push_back("omega_CP");
  return labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void C2HDMSympy::ReadAndSet(const std::string &linestr,
                            std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k < 10; k++)
  {
    ss >> tmp;
    if (k == 1)
      Type = static_cast<THDMType>(tmp);
    else if (k == 2)
      lambda1 = tmp;
    else if (k == 3)
      lambda2 = tmp;
    else if (k == 4)
      lambda3 = tmp;
    else if (k == 5)
      lambda4 = tmp;
    else if (k == 6)
      Relambda5 = tmp;
    else if (k == 7)
      Imlambda5 = tmp;
    else if (k == 8)
      Rem12sq = tmp;
    else if (k == 9)
      TanBeta = tmp;
  }

  par[0] = lambda1;
  par[1] = lambda2;
  par[2] = lambda3;
  par[3] = lambda4;
  par[4] = Relambda5;
  par[5] = Imlambda5;
  par[6] = Rem12sq;
  par[7] = TanBeta;
  par[8] = static_cast<int>(Type);

  set_gen(par);
  return;
}

/**
 * Set Class Object as well as the VEV configuration
 */
void C2HDMSympy::set_gen(const std::vector<double> &p)
{
  scale            = C_vev0;
  lambda1          = p[0];
  lambda2          = p[1];
  lambda3          = p[2];
  lambda4          = p[3];
  Relambda5        = p[4];
  Imlambda5        = p[5];
  Rem12sq          = p[6];
  TanBeta          = p[7];
  Type             = static_cast<THDMType>(p[8]);
  C_CosBetaSquared = 1.0 / (1 + TanBeta * TanBeta);
  C_CosBeta        = sqrt(C_CosBetaSquared);
  C_SinBetaSquared = TanBeta * TanBeta * C_CosBetaSquared;
  C_SinBeta        = sqrt(C_SinBetaSquared);

  m11sq = Rem12sq * TanBeta -
          C_vev0 * C_vev0 * C_SinBetaSquared * (lambda4 + Relambda5 + lambda3) /
              0.2e1 -
          C_vev0 * C_vev0 * C_CosBetaSquared * lambda1 / 0.2e1;
  m22sq = Rem12sq * 1.0 / TanBeta -
          C_vev0 * C_vev0 * C_CosBetaSquared * (lambda4 + Relambda5 + lambda3) /
              0.2e1 -
          C_vev0 * C_vev0 * C_SinBetaSquared * lambda2 / 0.2e1;
  Imm12sq = C_vev0 * C_vev0 * TanBeta * C_CosBetaSquared * Imlambda5 * 0.5;

  Relambda6 = 0;
  Relambda7 = 0;
  Imlambda6 = 0;
  Imlambda7 = 0;

  double cb = 0;

  if (Type == THDMType::TypeI or
      Type == THDMType::LeptonSpecific) // Type I 2HDM or Lepton Specific
  {
    cb = std::sqrt(2) * C_MassBottom / (C_vev0 * C_SinBeta);
  }
  if (Type == THDMType::TypeII or
      Type == THDMType::Flipped) // Type II 2HDM or Flipped
  {
    cb = std::sqrt(2) * C_MassBottom / (C_vev0 * C_CosBeta);
  }
  CTempC1 = 1.0 / 48 *
            (12 * lambda1 + 8 * lambda3 + 4 * lambda4 +
             3 * (3 * C_g * C_g + C_gs * C_gs));
  double ct = std::sqrt(2) * C_MassTop / (C_vev0 * C_SinBeta);
  CTempC2   = 1.0 / 48 *
            (12 * lambda2 + 8 * lambda3 + 4 * lambda4 +
             3 * (3 * C_g * C_g + C_gs * C_gs) + 12 * ct * ct);

  if (Type == THDMType::TypeI or Type == THDMType::LeptonSpecific)
  {
    CTempC2 += 12.0 / 48.0 * cb * cb;
  }
  else
  {
    CTempC1 += 12.0 / 48.0 * cb * cb;
  }
  vevTreeMin.resize(nVEV);
  vevTreeMin[0] = 0;
  vevTreeMin[1] = C_vev0 * C_CosBeta;
  vevTreeMin[2] = C_vev0 * C_SinBeta;
  vevTreeMin[3] = 0;
  vevTree.resize(NHiggs);
  vevTree = MinimizeOrderVEV(vevTreeMin);

  double v1 = C_vev0 * C_CosBeta;
  double v2 = C_vev0 * C_SinBeta;

  double lambda1Calc =
      -Relambda5 * std::pow(v2, 2) / std::pow(v1, 2) - 3 * Relambda6 * v2 / v1 -
      Relambda7 * std::pow(v2, 3) / std::pow(v1, 3) +
      2 * Rem12sq * v2 / std::pow(v1, 3) -
      lambda3 * std::pow(v2, 2) / std::pow(v1, 2) -
      lambda4 * std::pow(v2, 2) / std::pow(v1, 2) - 2 * m11sq / std::pow(v1, 2);

  double lambda2Calc =
      -Relambda5 * std::pow(v1, 2) / std::pow(v2, 2) -
      Relambda6 * std::pow(v1, 3) / std::pow(v2, 3) - 3 * Relambda7 * v1 / v2 +
      2 * Rem12sq * v1 / std::pow(v2, 3) -
      lambda3 * std::pow(v1, 2) / std::pow(v2, 2) -
      lambda4 * std::pow(v1, 2) / std::pow(v2, 2) - 2 * m22sq / std::pow(v2, 2);

  double diffLambda1 = (lambda1Calc - lambda1) / lambda1Calc;
  if (std::abs(diffLambda1) > 1e-5)
  {
    std::string msg = "lambda1 = " + std::to_string(lambda1) +
                      "; lambda1Calc = " + std::to_string(lambda1Calc) +
                      "; diff = " + std::to_string(diffLambda1);
    throw std::runtime_error(msg.c_str());
  }

  double diffLambda2 = (lambda2Calc - lambda2) / lambda2Calc;
  if (std::abs(diffLambda2) > 1e-5)
  {
    throw std::runtime_error("lambda2 not equal");
  }

  if (!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the
 * entries of Curvature_Higgs_CT_L1 to Curvature_Higgs_CT_L4.
 */
void C2HDMSympy::set_CT_Pot_Par(const std::vector<double> &par)
{

  // Begin CT Order for set_CT_Pot_Par
  dlambda1   = par.at(0);
  dlambda2   = par.at(1);
  dlambda3   = par.at(2);
  dlambda4   = par.at(3);
  dRelambda5 = par.at(4);
  dImlambda5 = par.at(5);
  dRelambda6 = par.at(6);
  dImlambda6 = par.at(7);
  dRelambda7 = par.at(8);
  dImlambda7 = par.at(9);
  dm11sq     = par.at(10);
  dm22sq     = par.at(11);
  dRem12sq   = par.at(12);
  dImm12sq   = par.at(13);
  dT1        = par.at(14);
  dT2        = par.at(15);
  dT3        = par.at(16);
  dT4        = par.at(17);
  dT5        = par.at(18);
  dT6        = par.at(19);
  dT7        = par.at(20);
  dT8        = par.at(21);
  // End CT Order

  // Begin of Higgs CT curvature tensors
  Curvature_Higgs_CT_L1.at(0)                   = dT1;
  Curvature_Higgs_CT_L1.at(1)                   = dT2;
  Curvature_Higgs_CT_L1.at(2)                   = dT3;
  Curvature_Higgs_CT_L1.at(3)                   = dT4;
  Curvature_Higgs_CT_L1.at(4)                   = dT5;
  Curvature_Higgs_CT_L1.at(5)                   = dT6;
  Curvature_Higgs_CT_L1.at(6)                   = dT7;
  Curvature_Higgs_CT_L1.at(7)                   = dT8;
  Curvature_Higgs_CT_L2.at(0).at(0)             = dm11sq;
  Curvature_Higgs_CT_L2.at(0).at(2)             = -dRem12sq;
  Curvature_Higgs_CT_L2.at(0).at(3)             = dImm12sq;
  Curvature_Higgs_CT_L2.at(1).at(1)             = dm11sq;
  Curvature_Higgs_CT_L2.at(1).at(2)             = -dImm12sq;
  Curvature_Higgs_CT_L2.at(1).at(3)             = -dRem12sq;
  Curvature_Higgs_CT_L2.at(2).at(0)             = -dRem12sq;
  Curvature_Higgs_CT_L2.at(2).at(1)             = -dImm12sq;
  Curvature_Higgs_CT_L2.at(2).at(2)             = dm22sq;
  Curvature_Higgs_CT_L2.at(3).at(0)             = dImm12sq;
  Curvature_Higgs_CT_L2.at(3).at(1)             = -dRem12sq;
  Curvature_Higgs_CT_L2.at(3).at(3)             = dm22sq;
  Curvature_Higgs_CT_L2.at(4).at(4)             = dm11sq;
  Curvature_Higgs_CT_L2.at(4).at(6)             = -dRem12sq;
  Curvature_Higgs_CT_L2.at(4).at(7)             = dImm12sq;
  Curvature_Higgs_CT_L2.at(5).at(5)             = dm11sq;
  Curvature_Higgs_CT_L2.at(5).at(6)             = -dImm12sq;
  Curvature_Higgs_CT_L2.at(5).at(7)             = -dRem12sq;
  Curvature_Higgs_CT_L2.at(6).at(4)             = -dRem12sq;
  Curvature_Higgs_CT_L2.at(6).at(5)             = -dImm12sq;
  Curvature_Higgs_CT_L2.at(6).at(6)             = dm22sq;
  Curvature_Higgs_CT_L2.at(7).at(4)             = dImm12sq;
  Curvature_Higgs_CT_L2.at(7).at(5)             = -dRem12sq;
  Curvature_Higgs_CT_L2.at(7).at(7)             = dm22sq;
  Curvature_Higgs_CT_L4.at(0).at(0).at(0).at(0) = 3 * dlambda1;
  Curvature_Higgs_CT_L4.at(0).at(0).at(0).at(2) = 3 * dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(0).at(0).at(3) = -3 * dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(0).at(1).at(1) = dlambda1;
  Curvature_Higgs_CT_L4.at(0).at(0).at(1).at(2) = dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(0).at(1).at(3) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(0).at(2).at(0) = 3 * dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(0).at(2).at(1) = dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(0).at(2).at(2) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(0).at(2).at(3) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(0).at(3).at(0) = -3 * dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(0).at(3).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(0).at(3).at(2) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(0).at(3).at(3) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(0).at(4).at(4) = dlambda1;
  Curvature_Higgs_CT_L4.at(0).at(0).at(4).at(6) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(0).at(4).at(7) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(0).at(5).at(5) = dlambda1;
  Curvature_Higgs_CT_L4.at(0).at(0).at(5).at(6) = dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(0).at(5).at(7) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(0).at(6).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(0).at(6).at(5) = dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(0).at(6).at(6) = dlambda3;
  Curvature_Higgs_CT_L4.at(0).at(0).at(7).at(4) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(0).at(7).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(0).at(7).at(7) = dlambda3;
  Curvature_Higgs_CT_L4.at(0).at(1).at(0).at(1) = dlambda1;
  Curvature_Higgs_CT_L4.at(0).at(1).at(0).at(2) = dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(1).at(0).at(3) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(1).at(1).at(0) = dlambda1;
  Curvature_Higgs_CT_L4.at(0).at(1).at(1).at(2) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(1).at(1).at(3) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(1).at(2).at(0) = dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(1).at(2).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(1).at(2).at(2) = dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(1).at(2).at(3) = dRelambda5;
  Curvature_Higgs_CT_L4.at(0).at(1).at(3).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(1).at(3).at(1) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(1).at(3).at(2) = dRelambda5;
  Curvature_Higgs_CT_L4.at(0).at(1).at(3).at(3) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(2).at(0).at(0) = 3 * dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(2).at(0).at(1) = dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(2).at(0).at(2) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(2).at(0).at(3) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(2).at(1).at(0) = dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(2).at(1).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(2).at(1).at(2) = dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(2).at(1).at(3) = dRelambda5;
  Curvature_Higgs_CT_L4.at(0).at(2).at(2).at(0) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(2).at(2).at(1) = dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(2).at(2).at(2) = 3 * dRelambda7;
  Curvature_Higgs_CT_L4.at(0).at(2).at(2).at(3) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(0).at(2).at(3).at(0) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(2).at(3).at(1) = dRelambda5;
  Curvature_Higgs_CT_L4.at(0).at(2).at(3).at(2) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(0).at(2).at(3).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(0).at(2).at(4).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(2).at(4).at(6) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(2).at(4).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(2).at(5).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(2).at(5).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(2).at(5).at(7) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(2).at(6).at(4) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(2).at(6).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(2).at(6).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(0).at(2).at(7).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(2).at(7).at(5) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(2).at(7).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(0).at(3).at(0).at(0) = -3 * dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(3).at(0).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(3).at(0).at(2) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(3).at(0).at(3) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(3).at(1).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(3).at(1).at(1) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(3).at(1).at(2) = dRelambda5;
  Curvature_Higgs_CT_L4.at(0).at(3).at(1).at(3) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(3).at(2).at(0) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(3).at(2).at(1) = dRelambda5;
  Curvature_Higgs_CT_L4.at(0).at(3).at(2).at(2) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(0).at(3).at(2).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(0).at(3).at(3).at(0) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(3).at(3).at(1) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(3).at(3).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(0).at(3).at(3).at(3) = -3 * dImlambda7;
  Curvature_Higgs_CT_L4.at(0).at(3).at(4).at(4) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(3).at(4).at(6) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(3).at(4).at(7) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(3).at(5).at(5) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(3).at(5).at(6) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(3).at(5).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(3).at(6).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(3).at(6).at(5) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(3).at(6).at(6) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(0).at(3).at(7).at(4) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(3).at(7).at(5) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(3).at(7).at(7) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(0).at(4).at(0).at(4) = dlambda1;
  Curvature_Higgs_CT_L4.at(0).at(4).at(0).at(6) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(4).at(0).at(7) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(4).at(2).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(4).at(2).at(6) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(4).at(2).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(4).at(3).at(4) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(4).at(3).at(6) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(4).at(3).at(7) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(4).at(4).at(0) = dlambda1;
  Curvature_Higgs_CT_L4.at(0).at(4).at(4).at(2) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(4).at(4).at(3) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(4).at(6).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(4).at(6).at(2) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(4).at(6).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(4).at(7).at(0) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(4).at(7).at(2) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(4).at(7).at(3) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(5).at(0).at(5) = dlambda1;
  Curvature_Higgs_CT_L4.at(0).at(5).at(0).at(6) = dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(5).at(0).at(7) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(5).at(2).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(5).at(2).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(5).at(2).at(7) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(5).at(3).at(5) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(5).at(3).at(6) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(5).at(3).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(5).at(5).at(0) = dlambda1;
  Curvature_Higgs_CT_L4.at(0).at(5).at(5).at(2) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(5).at(5).at(3) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(5).at(6).at(0) = dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(5).at(6).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(5).at(6).at(3) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(5).at(7).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(5).at(7).at(2) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(5).at(7).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(6).at(0).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(6).at(0).at(5) = dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(6).at(0).at(6) = dlambda3;
  Curvature_Higgs_CT_L4.at(0).at(6).at(2).at(4) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(6).at(2).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(6).at(2).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(0).at(6).at(3).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(6).at(3).at(5) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(6).at(3).at(6) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(0).at(6).at(4).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(6).at(4).at(2) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(6).at(4).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(6).at(5).at(0) = dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(6).at(5).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(6).at(5).at(3) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(6).at(6).at(0) = dlambda3;
  Curvature_Higgs_CT_L4.at(0).at(6).at(6).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(0).at(6).at(6).at(3) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(0).at(7).at(0).at(4) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(7).at(0).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(7).at(0).at(7) = dlambda3;
  Curvature_Higgs_CT_L4.at(0).at(7).at(2).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(7).at(2).at(5) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(7).at(2).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(0).at(7).at(3).at(4) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(7).at(3).at(5) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(7).at(3).at(7) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(0).at(7).at(4).at(0) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(0).at(7).at(4).at(2) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(7).at(4).at(3) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(7).at(5).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(0).at(7).at(5).at(2) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(0).at(7).at(5).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(0).at(7).at(7).at(0) = dlambda3;
  Curvature_Higgs_CT_L4.at(0).at(7).at(7).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(0).at(7).at(7).at(3) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(1).at(0).at(0).at(1) = dlambda1;
  Curvature_Higgs_CT_L4.at(1).at(0).at(0).at(2) = dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(0).at(0).at(3) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(0).at(1).at(0) = dlambda1;
  Curvature_Higgs_CT_L4.at(1).at(0).at(1).at(2) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(0).at(1).at(3) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(0).at(2).at(0) = dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(0).at(2).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(0).at(2).at(2) = dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(0).at(2).at(3) = dRelambda5;
  Curvature_Higgs_CT_L4.at(1).at(0).at(3).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(0).at(3).at(1) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(0).at(3).at(2) = dRelambda5;
  Curvature_Higgs_CT_L4.at(1).at(0).at(3).at(3) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(1).at(0).at(0) = dlambda1;
  Curvature_Higgs_CT_L4.at(1).at(1).at(0).at(2) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(1).at(0).at(3) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(1).at(1).at(1) = 3 * dlambda1;
  Curvature_Higgs_CT_L4.at(1).at(1).at(1).at(2) = 3 * dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(1).at(1).at(3) = 3 * dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(1).at(2).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(1).at(2).at(1) = 3 * dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(1).at(2).at(2) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(1).at(2).at(3) = dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(1).at(3).at(0) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(1).at(3).at(1) = 3 * dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(1).at(3).at(2) = dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(1).at(3).at(3) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(1).at(4).at(4) = dlambda1;
  Curvature_Higgs_CT_L4.at(1).at(1).at(4).at(6) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(1).at(4).at(7) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(1).at(5).at(5) = dlambda1;
  Curvature_Higgs_CT_L4.at(1).at(1).at(5).at(6) = dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(1).at(5).at(7) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(1).at(6).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(1).at(6).at(5) = dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(1).at(6).at(6) = dlambda3;
  Curvature_Higgs_CT_L4.at(1).at(1).at(7).at(4) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(1).at(7).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(1).at(7).at(7) = dlambda3;
  Curvature_Higgs_CT_L4.at(1).at(2).at(0).at(0) = dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(2).at(0).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(2).at(0).at(2) = dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(2).at(0).at(3) = dRelambda5;
  Curvature_Higgs_CT_L4.at(1).at(2).at(1).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(2).at(1).at(1) = 3 * dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(2).at(1).at(2) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(2).at(1).at(3) = dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(2).at(2).at(0) = dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(2).at(2).at(1) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(2).at(2).at(2) = 3 * dImlambda7;
  Curvature_Higgs_CT_L4.at(1).at(2).at(2).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(1).at(2).at(3).at(0) = dRelambda5;
  Curvature_Higgs_CT_L4.at(1).at(2).at(3).at(1) = dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(2).at(3).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(1).at(2).at(3).at(3) = dImlambda7;
  Curvature_Higgs_CT_L4.at(1).at(2).at(4).at(4) = dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(2).at(4).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(2).at(4).at(7) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(2).at(5).at(5) = dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(2).at(5).at(6) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(2).at(5).at(7) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(2).at(6).at(4) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(2).at(6).at(5) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(2).at(6).at(6) = dImlambda7;
  Curvature_Higgs_CT_L4.at(1).at(2).at(7).at(4) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(2).at(7).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(2).at(7).at(7) = dImlambda7;
  Curvature_Higgs_CT_L4.at(1).at(3).at(0).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(3).at(0).at(1) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(3).at(0).at(2) = dRelambda5;
  Curvature_Higgs_CT_L4.at(1).at(3).at(0).at(3) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(3).at(1).at(0) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(3).at(1).at(1) = 3 * dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(3).at(1).at(2) = dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(3).at(1).at(3) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(3).at(2).at(0) = dRelambda5;
  Curvature_Higgs_CT_L4.at(1).at(3).at(2).at(1) = dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(3).at(2).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(1).at(3).at(2).at(3) = dImlambda7;
  Curvature_Higgs_CT_L4.at(1).at(3).at(3).at(0) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(3).at(3).at(1) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(3).at(3).at(2) = dImlambda7;
  Curvature_Higgs_CT_L4.at(1).at(3).at(3).at(3) = 3 * dRelambda7;
  Curvature_Higgs_CT_L4.at(1).at(3).at(4).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(3).at(4).at(6) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(3).at(4).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(3).at(5).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(3).at(5).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(3).at(5).at(7) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(3).at(6).at(4) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(3).at(6).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(3).at(6).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(1).at(3).at(7).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(3).at(7).at(5) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(3).at(7).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(1).at(4).at(1).at(4) = dlambda1;
  Curvature_Higgs_CT_L4.at(1).at(4).at(1).at(6) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(4).at(1).at(7) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(4).at(2).at(4) = dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(4).at(2).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(4).at(2).at(7) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(4).at(3).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(4).at(3).at(6) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(4).at(3).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(4).at(4).at(1) = dlambda1;
  Curvature_Higgs_CT_L4.at(1).at(4).at(4).at(2) = dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(4).at(4).at(3) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(4).at(6).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(4).at(6).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(4).at(6).at(3) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(4).at(7).at(1) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(4).at(7).at(2) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(4).at(7).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(5).at(1).at(5) = dlambda1;
  Curvature_Higgs_CT_L4.at(1).at(5).at(1).at(6) = dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(5).at(1).at(7) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(5).at(2).at(5) = dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(5).at(2).at(6) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(5).at(2).at(7) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(5).at(3).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(5).at(3).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(5).at(3).at(7) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(5).at(5).at(1) = dlambda1;
  Curvature_Higgs_CT_L4.at(1).at(5).at(5).at(2) = dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(5).at(5).at(3) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(5).at(6).at(1) = dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(5).at(6).at(2) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(5).at(6).at(3) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(5).at(7).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(5).at(7).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(5).at(7).at(3) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(6).at(1).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(6).at(1).at(5) = dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(6).at(1).at(6) = dlambda3;
  Curvature_Higgs_CT_L4.at(1).at(6).at(2).at(4) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(6).at(2).at(5) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(6).at(2).at(6) = dImlambda7;
  Curvature_Higgs_CT_L4.at(1).at(6).at(3).at(4) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(6).at(3).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(6).at(3).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(1).at(6).at(4).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(6).at(4).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(6).at(4).at(3) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(6).at(5).at(1) = dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(6).at(5).at(2) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(6).at(5).at(3) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(6).at(6).at(1) = dlambda3;
  Curvature_Higgs_CT_L4.at(1).at(6).at(6).at(2) = dImlambda7;
  Curvature_Higgs_CT_L4.at(1).at(6).at(6).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(1).at(7).at(1).at(4) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(7).at(1).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(7).at(1).at(7) = dlambda3;
  Curvature_Higgs_CT_L4.at(1).at(7).at(2).at(4) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(7).at(2).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(7).at(2).at(7) = dImlambda7;
  Curvature_Higgs_CT_L4.at(1).at(7).at(3).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(7).at(3).at(5) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(7).at(3).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(1).at(7).at(4).at(1) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(1).at(7).at(4).at(2) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(7).at(4).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(7).at(5).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(1).at(7).at(5).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(1).at(7).at(5).at(3) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(1).at(7).at(7).at(1) = dlambda3;
  Curvature_Higgs_CT_L4.at(1).at(7).at(7).at(2) = dImlambda7;
  Curvature_Higgs_CT_L4.at(1).at(7).at(7).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(0).at(0).at(0) = 3 * dRelambda6;
  Curvature_Higgs_CT_L4.at(2).at(0).at(0).at(1) = dImlambda6;
  Curvature_Higgs_CT_L4.at(2).at(0).at(0).at(2) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(0).at(0).at(3) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(0).at(1).at(0) = dImlambda6;
  Curvature_Higgs_CT_L4.at(2).at(0).at(1).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(2).at(0).at(1).at(2) = dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(0).at(1).at(3) = dRelambda5;
  Curvature_Higgs_CT_L4.at(2).at(0).at(2).at(0) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(0).at(2).at(1) = dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(0).at(2).at(2) = 3 * dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(0).at(2).at(3) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(0).at(3).at(0) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(0).at(3).at(1) = dRelambda5;
  Curvature_Higgs_CT_L4.at(2).at(0).at(3).at(2) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(0).at(3).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(0).at(4).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(2).at(0).at(4).at(6) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(0).at(4).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(0).at(5).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(2).at(0).at(5).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(0).at(5).at(7) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(0).at(6).at(4) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(0).at(6).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(0).at(6).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(0).at(7).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(0).at(7).at(5) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(0).at(7).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(1).at(0).at(0) = dImlambda6;
  Curvature_Higgs_CT_L4.at(2).at(1).at(0).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(2).at(1).at(0).at(2) = dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(1).at(0).at(3) = dRelambda5;
  Curvature_Higgs_CT_L4.at(2).at(1).at(1).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(2).at(1).at(1).at(1) = 3 * dImlambda6;
  Curvature_Higgs_CT_L4.at(2).at(1).at(1).at(2) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(1).at(1).at(3) = dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(1).at(2).at(0) = dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(1).at(2).at(1) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(1).at(2).at(2) = 3 * dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(1).at(2).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(1).at(3).at(0) = dRelambda5;
  Curvature_Higgs_CT_L4.at(2).at(1).at(3).at(1) = dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(1).at(3).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(1).at(3).at(3) = dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(1).at(4).at(4) = dImlambda6;
  Curvature_Higgs_CT_L4.at(2).at(1).at(4).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(1).at(4).at(7) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(1).at(5).at(5) = dImlambda6;
  Curvature_Higgs_CT_L4.at(2).at(1).at(5).at(6) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(1).at(5).at(7) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(1).at(6).at(4) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(1).at(6).at(5) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(1).at(6).at(6) = dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(1).at(7).at(4) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(1).at(7).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(1).at(7).at(7) = dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(2).at(0).at(0) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(2).at(0).at(1) = dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(2).at(0).at(2) = 3 * dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(2).at(0).at(3) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(2).at(1).at(0) = dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(2).at(1).at(1) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(2).at(1).at(2) = 3 * dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(2).at(1).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(2).at(2).at(0) = 3 * dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(2).at(2).at(1) = 3 * dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(2).at(2).at(2) = 3 * dlambda2;
  Curvature_Higgs_CT_L4.at(2).at(2).at(3).at(0) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(2).at(3).at(1) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(2).at(3).at(3) = dlambda2;
  Curvature_Higgs_CT_L4.at(2).at(2).at(4).at(4) = dlambda3;
  Curvature_Higgs_CT_L4.at(2).at(2).at(4).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(2).at(4).at(7) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(2).at(5).at(5) = dlambda3;
  Curvature_Higgs_CT_L4.at(2).at(2).at(5).at(6) = dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(2).at(5).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(2).at(6).at(4) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(2).at(6).at(5) = dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(2).at(6).at(6) = dlambda2;
  Curvature_Higgs_CT_L4.at(2).at(2).at(7).at(4) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(2).at(7).at(5) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(2).at(7).at(7) = dlambda2;
  Curvature_Higgs_CT_L4.at(2).at(3).at(0).at(0) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(3).at(0).at(1) = dRelambda5;
  Curvature_Higgs_CT_L4.at(2).at(3).at(0).at(2) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(3).at(0).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(3).at(1).at(0) = dRelambda5;
  Curvature_Higgs_CT_L4.at(2).at(3).at(1).at(1) = dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(3).at(1).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(3).at(1).at(3) = dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(3).at(2).at(0) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(3).at(2).at(1) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(3).at(2).at(3) = dlambda2;
  Curvature_Higgs_CT_L4.at(2).at(3).at(3).at(0) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(3).at(3).at(1) = dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(3).at(3).at(2) = dlambda2;
  Curvature_Higgs_CT_L4.at(2).at(4).at(0).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(2).at(4).at(0).at(6) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(4).at(0).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(4).at(1).at(4) = dImlambda6;
  Curvature_Higgs_CT_L4.at(2).at(4).at(1).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(4).at(1).at(7) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(4).at(2).at(4) = dlambda3;
  Curvature_Higgs_CT_L4.at(2).at(4).at(2).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(4).at(2).at(7) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(4).at(4).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(2).at(4).at(4).at(1) = dImlambda6;
  Curvature_Higgs_CT_L4.at(2).at(4).at(4).at(2) = dlambda3;
  Curvature_Higgs_CT_L4.at(2).at(4).at(6).at(0) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(4).at(6).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(4).at(6).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(4).at(7).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(4).at(7).at(1) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(4).at(7).at(2) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(5).at(0).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(2).at(5).at(0).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(5).at(0).at(7) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(5).at(1).at(5) = dImlambda6;
  Curvature_Higgs_CT_L4.at(2).at(5).at(1).at(6) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(5).at(1).at(7) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(5).at(2).at(5) = dlambda3;
  Curvature_Higgs_CT_L4.at(2).at(5).at(2).at(6) = dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(5).at(2).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(5).at(5).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(2).at(5).at(5).at(1) = dImlambda6;
  Curvature_Higgs_CT_L4.at(2).at(5).at(5).at(2) = dlambda3;
  Curvature_Higgs_CT_L4.at(2).at(5).at(6).at(0) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(5).at(6).at(1) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(5).at(6).at(2) = dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(5).at(7).at(0) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(5).at(7).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(5).at(7).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(6).at(0).at(4) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(6).at(0).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(6).at(0).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(6).at(1).at(4) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(6).at(1).at(5) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(6).at(1).at(6) = dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(6).at(2).at(4) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(6).at(2).at(5) = dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(6).at(2).at(6) = dlambda2;
  Curvature_Higgs_CT_L4.at(2).at(6).at(4).at(0) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(6).at(4).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(6).at(4).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(6).at(5).at(0) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(6).at(5).at(1) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(6).at(5).at(2) = dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(6).at(6).at(0) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(6).at(6).at(1) = dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(6).at(6).at(2) = dlambda2;
  Curvature_Higgs_CT_L4.at(2).at(7).at(0).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(7).at(0).at(5) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(7).at(0).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(7).at(1).at(4) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(7).at(1).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(7).at(1).at(7) = dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(7).at(2).at(4) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(7).at(2).at(5) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(7).at(2).at(7) = dlambda2;
  Curvature_Higgs_CT_L4.at(2).at(7).at(4).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(7).at(4).at(1) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(7).at(4).at(2) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(7).at(5).at(0) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(2).at(7).at(5).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(2).at(7).at(5).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(7).at(7).at(0) = dRelambda7;
  Curvature_Higgs_CT_L4.at(2).at(7).at(7).at(1) = dImlambda7;
  Curvature_Higgs_CT_L4.at(2).at(7).at(7).at(2) = dlambda2;
  Curvature_Higgs_CT_L4.at(3).at(0).at(0).at(0) = -3 * dImlambda6;
  Curvature_Higgs_CT_L4.at(3).at(0).at(0).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(3).at(0).at(0).at(2) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(0).at(0).at(3) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(0).at(1).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(3).at(0).at(1).at(1) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(3).at(0).at(1).at(2) = dRelambda5;
  Curvature_Higgs_CT_L4.at(3).at(0).at(1).at(3) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(0).at(2).at(0) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(0).at(2).at(1) = dRelambda5;
  Curvature_Higgs_CT_L4.at(3).at(0).at(2).at(2) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(0).at(2).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(0).at(3).at(0) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(0).at(3).at(1) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(0).at(3).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(0).at(3).at(3) = -3 * dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(0).at(4).at(4) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(3).at(0).at(4).at(6) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(0).at(4).at(7) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(0).at(5).at(5) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(3).at(0).at(5).at(6) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(0).at(5).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(0).at(6).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(0).at(6).at(5) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(0).at(6).at(6) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(0).at(7).at(4) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(0).at(7).at(5) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(0).at(7).at(7) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(1).at(0).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(3).at(1).at(0).at(1) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(3).at(1).at(0).at(2) = dRelambda5;
  Curvature_Higgs_CT_L4.at(3).at(1).at(0).at(3) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(1).at(1).at(0) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(3).at(1).at(1).at(1) = 3 * dRelambda6;
  Curvature_Higgs_CT_L4.at(3).at(1).at(1).at(2) = dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(1).at(1).at(3) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(1).at(2).at(0) = dRelambda5;
  Curvature_Higgs_CT_L4.at(3).at(1).at(2).at(1) = dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(1).at(2).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(1).at(2).at(3) = dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(1).at(3).at(0) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(1).at(3).at(1) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(1).at(3).at(2) = dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(1).at(3).at(3) = 3 * dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(1).at(4).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(3).at(1).at(4).at(6) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(1).at(4).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(1).at(5).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(3).at(1).at(5).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(1).at(5).at(7) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(1).at(6).at(4) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(1).at(6).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(1).at(6).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(1).at(7).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(1).at(7).at(5) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(1).at(7).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(2).at(0).at(0) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(2).at(0).at(1) = dRelambda5;
  Curvature_Higgs_CT_L4.at(3).at(2).at(0).at(2) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(2).at(0).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(2).at(1).at(0) = dRelambda5;
  Curvature_Higgs_CT_L4.at(3).at(2).at(1).at(1) = dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(2).at(1).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(2).at(1).at(3) = dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(2).at(2).at(0) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(2).at(2).at(1) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(2).at(2).at(3) = dlambda2;
  Curvature_Higgs_CT_L4.at(3).at(2).at(3).at(0) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(2).at(3).at(1) = dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(2).at(3).at(2) = dlambda2;
  Curvature_Higgs_CT_L4.at(3).at(3).at(0).at(0) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(3).at(0).at(1) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(3).at(0).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(3).at(0).at(3) = -3 * dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(3).at(1).at(0) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(3).at(1).at(1) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(3).at(1).at(2) = dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(3).at(1).at(3) = 3 * dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(3).at(2).at(0) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(3).at(2).at(1) = dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(3).at(2).at(2) = dlambda2;
  Curvature_Higgs_CT_L4.at(3).at(3).at(3).at(0) = -3 * dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(3).at(3).at(1) = 3 * dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(3).at(3).at(3) = 3 * dlambda2;
  Curvature_Higgs_CT_L4.at(3).at(3).at(4).at(4) = dlambda3;
  Curvature_Higgs_CT_L4.at(3).at(3).at(4).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(3).at(4).at(7) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(3).at(5).at(5) = dlambda3;
  Curvature_Higgs_CT_L4.at(3).at(3).at(5).at(6) = dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(3).at(5).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(3).at(6).at(4) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(3).at(6).at(5) = dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(3).at(6).at(6) = dlambda2;
  Curvature_Higgs_CT_L4.at(3).at(3).at(7).at(4) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(3).at(7).at(5) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(3).at(7).at(7) = dlambda2;
  Curvature_Higgs_CT_L4.at(3).at(4).at(0).at(4) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(3).at(4).at(0).at(6) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(4).at(0).at(7) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(4).at(1).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(3).at(4).at(1).at(6) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(4).at(1).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(4).at(3).at(4) = dlambda3;
  Curvature_Higgs_CT_L4.at(3).at(4).at(3).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(4).at(3).at(7) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(4).at(4).at(0) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(3).at(4).at(4).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(3).at(4).at(4).at(3) = dlambda3;
  Curvature_Higgs_CT_L4.at(3).at(4).at(6).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(4).at(6).at(1) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(4).at(6).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(4).at(7).at(0) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(4).at(7).at(1) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(4).at(7).at(3) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(5).at(0).at(5) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(3).at(5).at(0).at(6) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(5).at(0).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(5).at(1).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(3).at(5).at(1).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(5).at(1).at(7) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(5).at(3).at(5) = dlambda3;
  Curvature_Higgs_CT_L4.at(3).at(5).at(3).at(6) = dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(5).at(3).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(5).at(5).at(0) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(3).at(5).at(5).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(3).at(5).at(5).at(3) = dlambda3;
  Curvature_Higgs_CT_L4.at(3).at(5).at(6).at(0) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(5).at(6).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(5).at(6).at(3) = dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(5).at(7).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(5).at(7).at(1) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(5).at(7).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(6).at(0).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(6).at(0).at(5) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(6).at(0).at(6) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(6).at(1).at(4) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(6).at(1).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(6).at(1).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(6).at(3).at(4) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(6).at(3).at(5) = dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(6).at(3).at(6) = dlambda2;
  Curvature_Higgs_CT_L4.at(3).at(6).at(4).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(6).at(4).at(1) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(6).at(4).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(6).at(5).at(0) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(6).at(5).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(6).at(5).at(3) = dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(6).at(6).at(0) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(6).at(6).at(1) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(6).at(6).at(3) = dlambda2;
  Curvature_Higgs_CT_L4.at(3).at(7).at(0).at(4) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(7).at(0).at(5) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(7).at(0).at(7) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(7).at(1).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(7).at(1).at(5) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(7).at(1).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(7).at(3).at(4) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(7).at(3).at(5) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(7).at(3).at(7) = dlambda2;
  Curvature_Higgs_CT_L4.at(3).at(7).at(4).at(0) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(7).at(4).at(1) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(7).at(4).at(3) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(7).at(5).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(3).at(7).at(5).at(1) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(3).at(7).at(5).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(7).at(7).at(0) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(3).at(7).at(7).at(1) = dRelambda7;
  Curvature_Higgs_CT_L4.at(3).at(7).at(7).at(3) = dlambda2;
  Curvature_Higgs_CT_L4.at(4).at(0).at(0).at(4) = dlambda1;
  Curvature_Higgs_CT_L4.at(4).at(0).at(0).at(6) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(0).at(0).at(7) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(0).at(2).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(0).at(2).at(6) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(0).at(2).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(0).at(3).at(4) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(0).at(3).at(6) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(0).at(3).at(7) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(0).at(4).at(0) = dlambda1;
  Curvature_Higgs_CT_L4.at(4).at(0).at(4).at(2) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(0).at(4).at(3) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(0).at(6).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(0).at(6).at(2) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(0).at(6).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(0).at(7).at(0) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(0).at(7).at(2) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(0).at(7).at(3) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(1).at(1).at(4) = dlambda1;
  Curvature_Higgs_CT_L4.at(4).at(1).at(1).at(6) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(1).at(1).at(7) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(1).at(2).at(4) = dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(1).at(2).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(1).at(2).at(7) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(1).at(3).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(1).at(3).at(6) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(1).at(3).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(1).at(4).at(1) = dlambda1;
  Curvature_Higgs_CT_L4.at(4).at(1).at(4).at(2) = dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(1).at(4).at(3) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(1).at(6).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(1).at(6).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(1).at(6).at(3) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(1).at(7).at(1) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(1).at(7).at(2) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(1).at(7).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(2).at(0).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(2).at(0).at(6) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(2).at(0).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(2).at(1).at(4) = dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(2).at(1).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(2).at(1).at(7) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(2).at(2).at(4) = dlambda3;
  Curvature_Higgs_CT_L4.at(4).at(2).at(2).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(4).at(2).at(2).at(7) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(4).at(2).at(4).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(2).at(4).at(1) = dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(2).at(4).at(2) = dlambda3;
  Curvature_Higgs_CT_L4.at(4).at(2).at(6).at(0) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(2).at(6).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(2).at(6).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(4).at(2).at(7).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(2).at(7).at(1) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(2).at(7).at(2) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(4).at(3).at(0).at(4) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(3).at(0).at(6) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(3).at(0).at(7) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(3).at(1).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(3).at(1).at(6) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(3).at(1).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(3).at(3).at(4) = dlambda3;
  Curvature_Higgs_CT_L4.at(4).at(3).at(3).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(4).at(3).at(3).at(7) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(4).at(3).at(4).at(0) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(3).at(4).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(3).at(4).at(3) = dlambda3;
  Curvature_Higgs_CT_L4.at(4).at(3).at(6).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(3).at(6).at(1) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(3).at(6).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(4).at(3).at(7).at(0) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(3).at(7).at(1) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(3).at(7).at(3) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(4).at(4).at(0).at(0) = dlambda1;
  Curvature_Higgs_CT_L4.at(4).at(4).at(0).at(2) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(4).at(0).at(3) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(4).at(1).at(1) = dlambda1;
  Curvature_Higgs_CT_L4.at(4).at(4).at(1).at(2) = dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(4).at(1).at(3) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(4).at(2).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(4).at(2).at(1) = dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(4).at(2).at(2) = dlambda3;
  Curvature_Higgs_CT_L4.at(4).at(4).at(3).at(0) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(4).at(3).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(4).at(3).at(3) = dlambda3;
  Curvature_Higgs_CT_L4.at(4).at(4).at(4).at(4) = 3 * dlambda1;
  Curvature_Higgs_CT_L4.at(4).at(4).at(4).at(6) = 3 * dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(4).at(4).at(7) = -3 * dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(4).at(5).at(5) = dlambda1;
  Curvature_Higgs_CT_L4.at(4).at(4).at(5).at(6) = dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(4).at(5).at(7) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(4).at(6).at(4) = 3 * dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(4).at(6).at(5) = dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(4).at(6).at(6) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(4).at(6).at(7) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(4).at(7).at(4) = -3 * dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(4).at(7).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(4).at(7).at(6) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(4).at(7).at(7) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(5).at(4).at(5) = dlambda1;
  Curvature_Higgs_CT_L4.at(4).at(5).at(4).at(6) = dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(5).at(4).at(7) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(5).at(5).at(4) = dlambda1;
  Curvature_Higgs_CT_L4.at(4).at(5).at(5).at(6) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(5).at(5).at(7) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(5).at(6).at(4) = dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(5).at(6).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(5).at(6).at(6) = dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(5).at(6).at(7) = dRelambda5;
  Curvature_Higgs_CT_L4.at(4).at(5).at(7).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(5).at(7).at(5) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(5).at(7).at(6) = dRelambda5;
  Curvature_Higgs_CT_L4.at(4).at(5).at(7).at(7) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(6).at(0).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(6).at(0).at(2) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(6).at(0).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(6).at(1).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(6).at(1).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(6).at(1).at(3) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(6).at(2).at(0) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(6).at(2).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(6).at(2).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(4).at(6).at(3).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(6).at(3).at(1) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(6).at(3).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(4).at(6).at(4).at(4) = 3 * dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(6).at(4).at(5) = dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(6).at(4).at(6) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(6).at(4).at(7) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(6).at(5).at(4) = dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(6).at(5).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(6).at(5).at(6) = dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(6).at(5).at(7) = dRelambda5;
  Curvature_Higgs_CT_L4.at(4).at(6).at(6).at(4) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(6).at(6).at(5) = dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(6).at(6).at(6) = 3 * dRelambda7;
  Curvature_Higgs_CT_L4.at(4).at(6).at(6).at(7) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(4).at(6).at(7).at(4) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(6).at(7).at(5) = dRelambda5;
  Curvature_Higgs_CT_L4.at(4).at(6).at(7).at(6) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(4).at(6).at(7).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(4).at(7).at(0).at(0) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(7).at(0).at(2) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(7).at(0).at(3) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(7).at(1).at(1) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(7).at(1).at(2) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(7).at(1).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(7).at(2).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(7).at(2).at(1) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(7).at(2).at(2) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(4).at(7).at(3).at(0) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(7).at(3).at(1) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(7).at(3).at(3) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(4).at(7).at(4).at(4) = -3 * dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(7).at(4).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(7).at(4).at(6) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(7).at(4).at(7) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(7).at(5).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(4).at(7).at(5).at(5) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(4).at(7).at(5).at(6) = dRelambda5;
  Curvature_Higgs_CT_L4.at(4).at(7).at(5).at(7) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(7).at(6).at(4) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(7).at(6).at(5) = dRelambda5;
  Curvature_Higgs_CT_L4.at(4).at(7).at(6).at(6) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(4).at(7).at(6).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(4).at(7).at(7).at(4) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(4).at(7).at(7).at(5) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(4).at(7).at(7).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(4).at(7).at(7).at(7) = -3 * dImlambda7;
  Curvature_Higgs_CT_L4.at(5).at(0).at(0).at(5) = dlambda1;
  Curvature_Higgs_CT_L4.at(5).at(0).at(0).at(6) = dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(0).at(0).at(7) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(0).at(2).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(0).at(2).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(0).at(2).at(7) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(0).at(3).at(5) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(0).at(3).at(6) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(0).at(3).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(0).at(5).at(0) = dlambda1;
  Curvature_Higgs_CT_L4.at(5).at(0).at(5).at(2) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(0).at(5).at(3) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(0).at(6).at(0) = dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(0).at(6).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(0).at(6).at(3) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(0).at(7).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(0).at(7).at(2) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(0).at(7).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(1).at(1).at(5) = dlambda1;
  Curvature_Higgs_CT_L4.at(5).at(1).at(1).at(6) = dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(1).at(1).at(7) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(1).at(2).at(5) = dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(1).at(2).at(6) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(1).at(2).at(7) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(1).at(3).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(1).at(3).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(1).at(3).at(7) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(1).at(5).at(1) = dlambda1;
  Curvature_Higgs_CT_L4.at(5).at(1).at(5).at(2) = dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(1).at(5).at(3) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(1).at(6).at(1) = dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(1).at(6).at(2) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(1).at(6).at(3) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(1).at(7).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(1).at(7).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(1).at(7).at(3) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(2).at(0).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(2).at(0).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(2).at(0).at(7) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(2).at(1).at(5) = dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(2).at(1).at(6) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(2).at(1).at(7) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(2).at(2).at(5) = dlambda3;
  Curvature_Higgs_CT_L4.at(5).at(2).at(2).at(6) = dImlambda7;
  Curvature_Higgs_CT_L4.at(5).at(2).at(2).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(5).at(2).at(5).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(2).at(5).at(1) = dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(2).at(5).at(2) = dlambda3;
  Curvature_Higgs_CT_L4.at(5).at(2).at(6).at(0) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(2).at(6).at(1) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(2).at(6).at(2) = dImlambda7;
  Curvature_Higgs_CT_L4.at(5).at(2).at(7).at(0) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(2).at(7).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(2).at(7).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(5).at(3).at(0).at(5) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(3).at(0).at(6) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(3).at(0).at(7) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(3).at(1).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(3).at(1).at(6) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(3).at(1).at(7) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(3).at(3).at(5) = dlambda3;
  Curvature_Higgs_CT_L4.at(5).at(3).at(3).at(6) = dImlambda7;
  Curvature_Higgs_CT_L4.at(5).at(3).at(3).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(5).at(3).at(5).at(0) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(3).at(5).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(3).at(5).at(3) = dlambda3;
  Curvature_Higgs_CT_L4.at(5).at(3).at(6).at(0) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(3).at(6).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(3).at(6).at(3) = dImlambda7;
  Curvature_Higgs_CT_L4.at(5).at(3).at(7).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(3).at(7).at(1) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(3).at(7).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(5).at(4).at(4).at(5) = dlambda1;
  Curvature_Higgs_CT_L4.at(5).at(4).at(4).at(6) = dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(4).at(4).at(7) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(4).at(5).at(4) = dlambda1;
  Curvature_Higgs_CT_L4.at(5).at(4).at(5).at(6) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(4).at(5).at(7) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(4).at(6).at(4) = dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(4).at(6).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(4).at(6).at(6) = dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(4).at(6).at(7) = dRelambda5;
  Curvature_Higgs_CT_L4.at(5).at(4).at(7).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(4).at(7).at(5) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(4).at(7).at(6) = dRelambda5;
  Curvature_Higgs_CT_L4.at(5).at(4).at(7).at(7) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(5).at(0).at(0) = dlambda1;
  Curvature_Higgs_CT_L4.at(5).at(5).at(0).at(2) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(5).at(0).at(3) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(5).at(1).at(1) = dlambda1;
  Curvature_Higgs_CT_L4.at(5).at(5).at(1).at(2) = dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(5).at(1).at(3) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(5).at(2).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(5).at(2).at(1) = dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(5).at(2).at(2) = dlambda3;
  Curvature_Higgs_CT_L4.at(5).at(5).at(3).at(0) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(5).at(3).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(5).at(3).at(3) = dlambda3;
  Curvature_Higgs_CT_L4.at(5).at(5).at(4).at(4) = dlambda1;
  Curvature_Higgs_CT_L4.at(5).at(5).at(4).at(6) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(5).at(4).at(7) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(5).at(5).at(5) = 3 * dlambda1;
  Curvature_Higgs_CT_L4.at(5).at(5).at(5).at(6) = 3 * dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(5).at(5).at(7) = 3 * dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(5).at(6).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(5).at(6).at(5) = 3 * dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(5).at(6).at(6) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(5).at(6).at(7) = dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(5).at(7).at(4) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(5).at(7).at(5) = 3 * dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(5).at(7).at(6) = dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(5).at(7).at(7) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(6).at(0).at(0) = dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(6).at(0).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(6).at(0).at(3) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(6).at(1).at(1) = dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(6).at(1).at(2) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(6).at(1).at(3) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(6).at(2).at(0) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(6).at(2).at(1) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(6).at(2).at(2) = dImlambda7;
  Curvature_Higgs_CT_L4.at(5).at(6).at(3).at(0) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(6).at(3).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(6).at(3).at(3) = dImlambda7;
  Curvature_Higgs_CT_L4.at(5).at(6).at(4).at(4) = dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(6).at(4).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(6).at(4).at(6) = dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(6).at(4).at(7) = dRelambda5;
  Curvature_Higgs_CT_L4.at(5).at(6).at(5).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(6).at(5).at(5) = 3 * dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(6).at(5).at(6) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(6).at(5).at(7) = dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(6).at(6).at(4) = dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(6).at(6).at(5) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(6).at(6).at(6) = 3 * dImlambda7;
  Curvature_Higgs_CT_L4.at(5).at(6).at(6).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(5).at(6).at(7).at(4) = dRelambda5;
  Curvature_Higgs_CT_L4.at(5).at(6).at(7).at(5) = dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(6).at(7).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(5).at(6).at(7).at(7) = dImlambda7;
  Curvature_Higgs_CT_L4.at(5).at(7).at(0).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(7).at(0).at(2) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(7).at(0).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(7).at(1).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(7).at(1).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(7).at(1).at(3) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(7).at(2).at(0) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(7).at(2).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(7).at(2).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(5).at(7).at(3).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(7).at(3).at(1) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(7).at(3).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(5).at(7).at(4).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(7).at(4).at(5) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(7).at(4).at(6) = dRelambda5;
  Curvature_Higgs_CT_L4.at(5).at(7).at(4).at(7) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(7).at(5).at(4) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(5).at(7).at(5).at(5) = 3 * dRelambda6;
  Curvature_Higgs_CT_L4.at(5).at(7).at(5).at(6) = dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(7).at(5).at(7) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(7).at(6).at(4) = dRelambda5;
  Curvature_Higgs_CT_L4.at(5).at(7).at(6).at(5) = dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(7).at(6).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(5).at(7).at(6).at(7) = dImlambda7;
  Curvature_Higgs_CT_L4.at(5).at(7).at(7).at(4) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(5).at(7).at(7).at(5) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(5).at(7).at(7).at(6) = dImlambda7;
  Curvature_Higgs_CT_L4.at(5).at(7).at(7).at(7) = 3 * dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(0).at(0).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(6).at(0).at(0).at(5) = dImlambda6;
  Curvature_Higgs_CT_L4.at(6).at(0).at(0).at(6) = dlambda3;
  Curvature_Higgs_CT_L4.at(6).at(0).at(2).at(4) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(0).at(2).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(0).at(2).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(0).at(3).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(0).at(3).at(5) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(0).at(3).at(6) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(0).at(4).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(6).at(0).at(4).at(2) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(0).at(4).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(0).at(5).at(0) = dImlambda6;
  Curvature_Higgs_CT_L4.at(6).at(0).at(5).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(0).at(5).at(3) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(0).at(6).at(0) = dlambda3;
  Curvature_Higgs_CT_L4.at(6).at(0).at(6).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(0).at(6).at(3) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(1).at(1).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(6).at(1).at(1).at(5) = dImlambda6;
  Curvature_Higgs_CT_L4.at(6).at(1).at(1).at(6) = dlambda3;
  Curvature_Higgs_CT_L4.at(6).at(1).at(2).at(4) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(1).at(2).at(5) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(1).at(2).at(6) = dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(1).at(3).at(4) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(1).at(3).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(1).at(3).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(1).at(4).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(6).at(1).at(4).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(1).at(4).at(3) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(1).at(5).at(1) = dImlambda6;
  Curvature_Higgs_CT_L4.at(6).at(1).at(5).at(2) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(1).at(5).at(3) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(1).at(6).at(1) = dlambda3;
  Curvature_Higgs_CT_L4.at(6).at(1).at(6).at(2) = dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(1).at(6).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(2).at(0).at(4) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(2).at(0).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(2).at(0).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(2).at(1).at(4) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(2).at(1).at(5) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(2).at(1).at(6) = dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(2).at(2).at(4) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(2).at(2).at(5) = dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(2).at(2).at(6) = dlambda2;
  Curvature_Higgs_CT_L4.at(6).at(2).at(4).at(0) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(2).at(4).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(2).at(4).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(2).at(5).at(0) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(2).at(5).at(1) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(2).at(5).at(2) = dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(2).at(6).at(0) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(2).at(6).at(1) = dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(2).at(6).at(2) = dlambda2;
  Curvature_Higgs_CT_L4.at(6).at(3).at(0).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(3).at(0).at(5) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(3).at(0).at(6) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(3).at(1).at(4) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(3).at(1).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(3).at(1).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(3).at(3).at(4) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(3).at(3).at(5) = dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(3).at(3).at(6) = dlambda2;
  Curvature_Higgs_CT_L4.at(6).at(3).at(4).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(3).at(4).at(1) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(3).at(4).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(3).at(5).at(0) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(3).at(5).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(3).at(5).at(3) = dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(3).at(6).at(0) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(3).at(6).at(1) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(3).at(6).at(3) = dlambda2;
  Curvature_Higgs_CT_L4.at(6).at(4).at(0).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(6).at(4).at(0).at(2) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(4).at(0).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(4).at(1).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(6).at(4).at(1).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(4).at(1).at(3) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(4).at(2).at(0) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(4).at(2).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(4).at(2).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(4).at(3).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(4).at(3).at(1) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(4).at(3).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(4).at(4).at(4) = 3 * dRelambda6;
  Curvature_Higgs_CT_L4.at(6).at(4).at(4).at(5) = dImlambda6;
  Curvature_Higgs_CT_L4.at(6).at(4).at(4).at(6) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(4).at(4).at(7) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(4).at(5).at(4) = dImlambda6;
  Curvature_Higgs_CT_L4.at(6).at(4).at(5).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(6).at(4).at(5).at(6) = dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(4).at(5).at(7) = dRelambda5;
  Curvature_Higgs_CT_L4.at(6).at(4).at(6).at(4) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(4).at(6).at(5) = dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(4).at(6).at(6) = 3 * dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(4).at(6).at(7) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(4).at(7).at(4) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(4).at(7).at(5) = dRelambda5;
  Curvature_Higgs_CT_L4.at(6).at(4).at(7).at(6) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(4).at(7).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(5).at(0).at(0) = dImlambda6;
  Curvature_Higgs_CT_L4.at(6).at(5).at(0).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(5).at(0).at(3) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(5).at(1).at(1) = dImlambda6;
  Curvature_Higgs_CT_L4.at(6).at(5).at(1).at(2) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(5).at(1).at(3) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(5).at(2).at(0) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(5).at(2).at(1) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(5).at(2).at(2) = dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(5).at(3).at(0) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(5).at(3).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(5).at(3).at(3) = dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(5).at(4).at(4) = dImlambda6;
  Curvature_Higgs_CT_L4.at(6).at(5).at(4).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(6).at(5).at(4).at(6) = dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(5).at(4).at(7) = dRelambda5;
  Curvature_Higgs_CT_L4.at(6).at(5).at(5).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(6).at(5).at(5).at(5) = 3 * dImlambda6;
  Curvature_Higgs_CT_L4.at(6).at(5).at(5).at(6) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(5).at(5).at(7) = dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(5).at(6).at(4) = dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(5).at(6).at(5) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(5).at(6).at(6) = 3 * dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(5).at(6).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(5).at(7).at(4) = dRelambda5;
  Curvature_Higgs_CT_L4.at(6).at(5).at(7).at(5) = dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(5).at(7).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(5).at(7).at(7) = dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(6).at(0).at(0) = dlambda3;
  Curvature_Higgs_CT_L4.at(6).at(6).at(0).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(6).at(0).at(3) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(6).at(1).at(1) = dlambda3;
  Curvature_Higgs_CT_L4.at(6).at(6).at(1).at(2) = dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(6).at(1).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(6).at(2).at(0) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(6).at(2).at(1) = dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(6).at(2).at(2) = dlambda2;
  Curvature_Higgs_CT_L4.at(6).at(6).at(3).at(0) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(6).at(3).at(1) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(6).at(3).at(3) = dlambda2;
  Curvature_Higgs_CT_L4.at(6).at(6).at(4).at(4) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(6).at(4).at(5) = dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(6).at(4).at(6) = 3 * dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(6).at(4).at(7) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(6).at(5).at(4) = dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(6).at(5).at(5) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(6).at(6).at(5).at(6) = 3 * dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(6).at(5).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(6).at(6).at(4) = 3 * dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(6).at(6).at(5) = 3 * dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(6).at(6).at(6) = 3 * dlambda2;
  Curvature_Higgs_CT_L4.at(6).at(6).at(7).at(4) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(6).at(7).at(5) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(6).at(7).at(7) = dlambda2;
  Curvature_Higgs_CT_L4.at(6).at(7).at(4).at(4) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(7).at(4).at(5) = dRelambda5;
  Curvature_Higgs_CT_L4.at(6).at(7).at(4).at(6) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(7).at(4).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(7).at(5).at(4) = dRelambda5;
  Curvature_Higgs_CT_L4.at(6).at(7).at(5).at(5) = dImlambda5;
  Curvature_Higgs_CT_L4.at(6).at(7).at(5).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(7).at(5).at(7) = dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(7).at(6).at(4) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(7).at(6).at(5) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(7).at(6).at(7) = dlambda2;
  Curvature_Higgs_CT_L4.at(6).at(7).at(7).at(4) = dRelambda7;
  Curvature_Higgs_CT_L4.at(6).at(7).at(7).at(5) = dImlambda7;
  Curvature_Higgs_CT_L4.at(6).at(7).at(7).at(6) = dlambda2;
  Curvature_Higgs_CT_L4.at(7).at(0).at(0).at(4) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(7).at(0).at(0).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(7).at(0).at(0).at(7) = dlambda3;
  Curvature_Higgs_CT_L4.at(7).at(0).at(2).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(0).at(2).at(5) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(0).at(2).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(0).at(3).at(4) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(0).at(3).at(5) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(0).at(3).at(7) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(0).at(4).at(0) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(7).at(0).at(4).at(2) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(0).at(4).at(3) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(0).at(5).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(7).at(0).at(5).at(2) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(0).at(5).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(0).at(7).at(0) = dlambda3;
  Curvature_Higgs_CT_L4.at(7).at(0).at(7).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(0).at(7).at(3) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(1).at(1).at(4) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(7).at(1).at(1).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(7).at(1).at(1).at(7) = dlambda3;
  Curvature_Higgs_CT_L4.at(7).at(1).at(2).at(4) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(1).at(2).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(1).at(2).at(7) = dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(1).at(3).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(1).at(3).at(5) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(1).at(3).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(1).at(4).at(1) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(7).at(1).at(4).at(2) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(1).at(4).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(1).at(5).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(7).at(1).at(5).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(1).at(5).at(3) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(1).at(7).at(1) = dlambda3;
  Curvature_Higgs_CT_L4.at(7).at(1).at(7).at(2) = dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(1).at(7).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(2).at(0).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(2).at(0).at(5) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(2).at(0).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(2).at(1).at(4) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(2).at(1).at(5) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(2).at(1).at(7) = dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(2).at(2).at(4) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(2).at(2).at(5) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(2).at(2).at(7) = dlambda2;
  Curvature_Higgs_CT_L4.at(7).at(2).at(4).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(2).at(4).at(1) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(2).at(4).at(2) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(2).at(5).at(0) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(2).at(5).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(2).at(5).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(2).at(7).at(0) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(2).at(7).at(1) = dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(2).at(7).at(2) = dlambda2;
  Curvature_Higgs_CT_L4.at(7).at(3).at(0).at(4) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(3).at(0).at(5) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(3).at(0).at(7) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(3).at(1).at(4) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(3).at(1).at(5) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(3).at(1).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(3).at(3).at(4) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(3).at(3).at(5) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(3).at(3).at(7) = dlambda2;
  Curvature_Higgs_CT_L4.at(7).at(3).at(4).at(0) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(3).at(4).at(1) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(3).at(4).at(3) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(3).at(5).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(3).at(5).at(1) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(3).at(5).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(3).at(7).at(0) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(3).at(7).at(1) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(3).at(7).at(3) = dlambda2;
  Curvature_Higgs_CT_L4.at(7).at(4).at(0).at(0) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(7).at(4).at(0).at(2) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(4).at(0).at(3) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(4).at(1).at(1) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(7).at(4).at(1).at(2) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(4).at(1).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(4).at(2).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(4).at(2).at(1) =
      (1.0 / 2.0) * dRelambda5 - 1.0 / 2.0 * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(4).at(2).at(2) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(4).at(3).at(0) =
      -1.0 / 2.0 * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(4).at(3).at(1) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(4).at(3).at(3) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(4).at(4).at(4) = -3 * dImlambda6;
  Curvature_Higgs_CT_L4.at(7).at(4).at(4).at(5) = dRelambda6;
  Curvature_Higgs_CT_L4.at(7).at(4).at(4).at(6) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(4).at(4).at(7) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(4).at(5).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(7).at(4).at(5).at(5) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(7).at(4).at(5).at(6) = dRelambda5;
  Curvature_Higgs_CT_L4.at(7).at(4).at(5).at(7) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(4).at(6).at(4) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(4).at(6).at(5) = dRelambda5;
  Curvature_Higgs_CT_L4.at(7).at(4).at(6).at(6) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(4).at(6).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(4).at(7).at(4) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(4).at(7).at(5) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(4).at(7).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(4).at(7).at(7) = -3 * dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(5).at(0).at(0) = dRelambda6;
  Curvature_Higgs_CT_L4.at(7).at(5).at(0).at(2) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(5).at(0).at(3) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(5).at(1).at(1) = dRelambda6;
  Curvature_Higgs_CT_L4.at(7).at(5).at(1).at(2) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(5).at(1).at(3) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(5).at(2).at(0) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(5).at(2).at(1) = (1.0 / 2.0) * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(5).at(2).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(5).at(3).at(0) = -1.0 / 2.0 * dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(5).at(3).at(1) =
      (1.0 / 2.0) * dRelambda5 + (1.0 / 2.0) * dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(5).at(3).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(5).at(4).at(4) = dRelambda6;
  Curvature_Higgs_CT_L4.at(7).at(5).at(4).at(5) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(7).at(5).at(4).at(6) = dRelambda5;
  Curvature_Higgs_CT_L4.at(7).at(5).at(4).at(7) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(5).at(5).at(4) = -dImlambda6;
  Curvature_Higgs_CT_L4.at(7).at(5).at(5).at(5) = 3 * dRelambda6;
  Curvature_Higgs_CT_L4.at(7).at(5).at(5).at(6) = dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(5).at(5).at(7) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(5).at(6).at(4) = dRelambda5;
  Curvature_Higgs_CT_L4.at(7).at(5).at(6).at(5) = dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(5).at(6).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(5).at(6).at(7) = dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(5).at(7).at(4) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(5).at(7).at(5) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(5).at(7).at(6) = dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(5).at(7).at(7) = 3 * dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(6).at(4).at(4) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(6).at(4).at(5) = dRelambda5;
  Curvature_Higgs_CT_L4.at(7).at(6).at(4).at(6) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(6).at(4).at(7) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(6).at(5).at(4) = dRelambda5;
  Curvature_Higgs_CT_L4.at(7).at(6).at(5).at(5) = dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(6).at(5).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(6).at(5).at(7) = dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(6).at(6).at(4) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(6).at(6).at(5) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(6).at(6).at(7) = dlambda2;
  Curvature_Higgs_CT_L4.at(7).at(6).at(7).at(4) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(6).at(7).at(5) = dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(6).at(7).at(6) = dlambda2;
  Curvature_Higgs_CT_L4.at(7).at(7).at(0).at(0) = dlambda3;
  Curvature_Higgs_CT_L4.at(7).at(7).at(0).at(2) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(7).at(0).at(3) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(7).at(1).at(1) = dlambda3;
  Curvature_Higgs_CT_L4.at(7).at(7).at(1).at(2) = dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(7).at(1).at(3) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(7).at(2).at(0) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(7).at(2).at(1) = dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(7).at(2).at(2) = dlambda2;
  Curvature_Higgs_CT_L4.at(7).at(7).at(3).at(0) = -dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(7).at(3).at(1) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(7).at(3).at(3) = dlambda2;
  Curvature_Higgs_CT_L4.at(7).at(7).at(4).at(4) =
      -dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(7).at(4).at(5) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(7).at(4).at(6) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(7).at(4).at(7) = -3 * dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(7).at(5).at(4) = -dImlambda5;
  Curvature_Higgs_CT_L4.at(7).at(7).at(5).at(5) =
      dRelambda5 + dlambda3 + dlambda4;
  Curvature_Higgs_CT_L4.at(7).at(7).at(5).at(6) = dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(7).at(5).at(7) = 3 * dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(7).at(6).at(4) = dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(7).at(6).at(5) = dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(7).at(6).at(6) = dlambda2;
  Curvature_Higgs_CT_L4.at(7).at(7).at(7).at(4) = -3 * dImlambda7;
  Curvature_Higgs_CT_L4.at(7).at(7).at(7).at(5) = 3 * dRelambda7;
  Curvature_Higgs_CT_L4.at(7).at(7).at(7).at(7) = 3 * dlambda2;
  // End of Higgs CT curvature tensors
}

/**
 * console output of all Parameters
 */
void C2HDMSympy::write() const
{

  std::stringstream ss;
  ss.precision(std::numeric_limits<double>::max_digits10);

  double MSM = 0, MhUp = 0, MhDown = 0;

  ss << "scale = " << scale << std::endl;

  ss << "The parameters are :  \n";
  ss << "Model = " << Model << "\n";
  ss << "v1 = " << C_vev0 * C_CosBeta << "\n";
  ss << "v2 = " << C_vev0 * C_SinBeta << "\n";
  ss << "Type = " << static_cast<int>(Type) << "\n";

  ss << "tan(beta) = " << TanBeta << std::endl;
  ss << "Lambda1 = " << lambda1 << std::endl;
  ss << "Lambda2 = " << lambda2 << std::endl;
  ss << "Lambda3 = " << lambda3 << std::endl;
  ss << "Lambda4 = " << lambda4 << std::endl;
  ss << "Re(Lambda5) = " << Relambda5 << std::endl;
  ss << "Im(Lambda5) = " << Imlambda5 << std::endl;
  ss << "Re(m_12^2) = " << Rem12sq << std::endl;
  ss << "m_{11}^2 = " << m11sq << std::endl;
  ss << "m_{22}^2 = " << m22sq << std::endl;
  ss << "Im(m_{12}^2) = " << Imm12sq << std::endl;

  ss << "The counterterms are :\n";

  ss << "DL1 := " << dlambda1 << ";\n";
  ss << "DL2 := " << dlambda2 << ";\n";
  ss << "DL3 := " << dlambda3 << ";\n";
  ss << "DL4 := " << dlambda4 << ";\n";
  ss << "DRL5 := " << dRelambda5 << ";\n";
  ss << "DIL5 := " << dImlambda5 << ";\n";
  ss << "Du1 := " << dm11sq << ";\n";
  ss << "Du2 := " << dm22sq << ";\n";
  ss << "DRu3 := " << dRem12sq << ";\n";
  ss << "DIu3 := " << dImm12sq << ";\n";
  ss << "DT1 := " << dT1 << ";\n";
  ss << "DT2 := " << dT2 << ";\n";
  ss << "DT3 := " << dT3 << ";\n";
  ss << "DT4 := " << dT4 << ";\n";
  ss << "DT5 := " << dT5 << ";\n";
  ss << "DT6 := " << dT6 << ";\n";
  ss << "DT7 := " << dT7 << ";\n";
  ss << "DT8 := " << dT8 << ";\n";
  ss << "DIL6:= " << dImlambda6 << ";\n";

  if (CalcCouplingsdone)
  {
    MatrixXd HiggsRot(NHiggs, NHiggs);
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      for (std::size_t j = 0; j < NHiggs; j++)
      {
        HiggsRot(i, j) = HiggsRotationMatrix[i][j];
      }
    }

    int posMHCS1 = 0;
    int posN[3];
    int countposN              = 0;
    int posG0                  = 0;
    double testsum             = 0;
    const double ZeroThreshold = 1e-5;
    for (int i = 0; i < 3; i++)
    {
      //    			testsum = std::abs(HiggsRot(i,0)) + std::abs(HiggsRot(i,2));
      //    			if(testsum > ZeroThreshold) posG1 = i;
      //    			testsum = std::abs(HiggsRot(i,1)) + std::abs(HiggsRot(i,3));
      //    			if(testsum > ZeroThreshold) posG2 = i;
      testsum = std::abs(HiggsRot(i, 5)) + std::abs(HiggsRot(i, 7));
      if (testsum > ZeroThreshold) posG0 = i;
    }
    for (std::size_t i = 3; i < NHiggs; i++)
    {
      testsum = std::abs(HiggsRot(i, 0)) + std::abs(HiggsRot(i, 2));
      if (testsum > ZeroThreshold) posMHCS1 = i;
      //    			testsum = std::abs(HiggsRot(i,1)) + std::abs(HiggsRot(i,3));
      //    			if(testsum > ZeroThreshold) posMHCS2 = i;
      testsum = 0;
      for (int k = 4; k < 8; k++)
        testsum += std::abs(HiggsRot(i, k));
      if (testsum > ZeroThreshold)
      {
        posN[countposN] = i;
        countposN++;
      }
    }

    std::vector<double> HiggsMasses;
    HiggsMasses = HiggsMassesSquared(vevTree, 0);

    double NeutralHiggs[3];
    for (int i = 0; i < 3; i++)
    {
      NeutralHiggs[i] = HiggsMasses[posN[i]];
    }
    for (int i = 0; i < 3; i++)
    {
      if (std::sqrt(NeutralHiggs[i]) < 126 and std::sqrt(NeutralHiggs[i]) > 124)
        MSM = std::sqrt(NeutralHiggs[i]);
    }
    if (std::sqrt(NeutralHiggs[0]) == MSM)
    {
      MhUp   = std::sqrt(NeutralHiggs[2]);
      MhDown = std::sqrt(NeutralHiggs[1]);
    }
    else if (std::sqrt(NeutralHiggs[1]) == MSM)
    {
      MhUp   = std::sqrt(NeutralHiggs[2]);
      MhDown = std::sqrt(NeutralHiggs[0]);
    }
    else
    {
      MhUp   = std::sqrt(NeutralHiggs[1]);
      MhDown = std::sqrt(NeutralHiggs[0]);
    }

    if (MSM > MhUp)
    {
      double tmp = posN[1];
      posN[1]    = posN[2];
      posN[2]    = tmp;
    }
    if (MSM > MhDown)
    {
      double tmp = posN[0];
      posN[0]    = posN[1];
      posN[1]    = tmp;
    }

    MatrixXd NeutralMatrix(4, 4);
    for (int j = 0; j < 4; j++)
      NeutralMatrix(0, j) = HiggsRot(posG0, j + 4);
    for (int i = 1; i < 4; i++)
    {
      for (int j = 0; j < 4; j++)
        NeutralMatrix(i, j) = HiggsRot(posN[i - 1], j + 4);
    }

    double beta = std::atan(TanBeta);
    MatrixXd MassMixing(3, 3);
    for (int i = 0; i < 3; i++)
    {
      MassMixing(i, 0) = NeutralMatrix(i + 1, 0);
      MassMixing(i, 1) = NeutralMatrix(i + 1, 2);
      MassMixing(i, 2) = -std::sin(beta) * NeutralMatrix(i + 1, 1) +
                         std::cos(beta) * NeutralMatrix(i + 1, 3);
    }

    ss << "The mass spectrum is given by :\n";
    ss << "m_{H^+} = " << std::sqrt(HiggsMasses[posMHCS1]) << " GeV \n"
       << "m_{H_SM} = " << MSM << " GeV \n"
       << "m_{H_l} = " << MhDown << " GeV \n"
       << "m_{H_h} = " << MhUp << " GeV \n";
    ss << "The neutral mixing Matrix is given by :\n";
    bool IsNegative = MassMixing(0, 1) < 0;
    ss << "H_{SM} = " << MassMixing(0, 0) << " zeta_1 ";
    if (IsNegative)
      ss << "-";
    else
      ss << "+";
    ss << std::abs(MassMixing(0, 1)) << " zeta_2 ";
    IsNegative = MassMixing(0, 2) < 0;
    if (IsNegative)
      ss << "-";
    else
      ss << "+";
    ss << std::abs(MassMixing(0, 2)) << " zeta_3 \n"
       << "H_{l} = " << MassMixing(1, 0) << " zeta_1 ";
    IsNegative = MassMixing(1, 1) < 0;
    if (IsNegative)
      ss << "-";
    else
      ss << "+";
    ss << std::abs(MassMixing(1, 1)) << " zeta_2 ";
    IsNegative = MassMixing(1, 2) < 0;
    if (IsNegative)
      ss << "-";
    else
      ss << "+";
    ss << std::abs(MassMixing(1, 2)) << " zeta_3 \n"
       << "H_{h} = " << MassMixing(2, 0) << " zeta_1 ";
    IsNegative = MassMixing(2, 1) < 0;
    if (IsNegative)
      ss << "-";
    else
      ss << "+";
    ss << std::abs(MassMixing(2, 1)) << " zeta_2 ";
    IsNegative = MassMixing(2, 2) < 0;
    if (IsNegative)
      ss << "-";
    else
      ss << "+";
    ss << std::abs(MassMixing(2, 2)) << " zeta_3 \n";
  }

  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms. Here you need to work out the scheme and
 * implement the formulas.
 */
std::vector<double> C2HDMSympy::calc_CT() const
{

  std::vector<double> parCT;

  if (!SetCurvatureDone)
  {
    std::string retmes = __func__;
    retmes += " was called before SetCurvatureArrays()!\n";
    throw std::runtime_error(retmes);
  }
  if (!CalcCouplingsdone)
  {
    std::string retmes = __func__;
    retmes += " was called before CalculatePhysicalCouplings()!\n";
    throw std::runtime_error(retmes);
  }

  std::vector<double> WeinbergNabla, WeinbergHesse;
  WeinbergNabla = WeinbergFirstDerivative();
  WeinbergHesse = WeinbergSecondDerivative();

  VectorXd NablaWeinberg(NHiggs);
  MatrixXd HesseWeinberg(NHiggs, NHiggs), HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    NablaWeinberg[i] = WeinbergNabla[i];
    for (std::size_t j = 0; j < NHiggs; j++)
      HesseWeinberg(i, j) = WeinbergHesse.at(j * NHiggs + i);
  }

  double v1 = C_vev0 * C_CosBeta;
  double v2 = C_vev0 * C_SinBeta;

  // Here you have to use your formulae for the counterterm scheme
  // Begin CT Calculation
  // Begin CT Calculation
  parCT.push_back(2 * HesseWeinberg(3, 3) * std::pow(v2, 2) / std::pow(v1, 4) -
                  HesseWeinberg(4, 4) / std::pow(v1, 2) +
                  HesseWeinberg(5, 5) / std::pow(v1, 2) -
                  2 * HesseWeinberg(7, 7) * std::pow(v2, 2) /
                      std::pow(v1, 4)); // dlambda1
  parCT.push_back(2 * HesseWeinberg(3, 3) / std::pow(v2, 2) -
                  HesseWeinberg(6, 6) / std::pow(v2, 2) -
                  HesseWeinberg(7, 7) / std::pow(v2, 2)); // dlambda2
  parCT.push_back(-HesseWeinberg(4, 6) / (v1 * v2) +
                  HesseWeinberg(5, 7) / (v1 * v2)); // dlambda3
  parCT.push_back(0);                               // dlambda4
  parCT.push_back(-2 * HesseWeinberg(3, 3) / std::pow(v1, 2) +
                  2 * HesseWeinberg(7, 7) / std::pow(v1, 2)); // dRelambda5
  parCT.push_back(2 * HesseWeinberg(6, 7) / std::pow(v1, 2)); // dImlambda5
  parCT.push_back(0);                                         // dRelambda6
  parCT.push_back(HesseWeinberg(4, 7) / std::pow(v1, 2) +
                  HesseWeinberg(5, 6) / std::pow(v1, 2)); // dImlambda6
  parCT.push_back(0);                                     // dRelambda7
  parCT.push_back(0);                                     // dImlambda7
  parCT.push_back(-2 * HesseWeinberg(3, 3) * std::pow(v2, 2) / std::pow(v1, 2) +
                  (1.0 / 2.0) * HesseWeinberg(4, 4) +
                  (1.0 / 2.0) * HesseWeinberg(4, 6) * v2 / v1 -
                  3.0 / 2.0 * HesseWeinberg(5, 5) -
                  1.0 / 2.0 * HesseWeinberg(5, 7) * v2 / v1 +
                  2 * HesseWeinberg(7, 7) * std::pow(v2, 2) /
                      std::pow(v1, 2)); // dm11sq
  parCT.push_back(-2 * HesseWeinberg(3, 3) +
                  (1.0 / 2.0) * HesseWeinberg(4, 6) * v1 / v2 -
                  1.0 / 2.0 * HesseWeinberg(5, 7) * v1 / v2 +
                  (1.0 / 2.0) * HesseWeinberg(6, 6) +
                  (1.0 / 2.0) * HesseWeinberg(7, 7)); // dm22sq
  parCT.push_back(-2 * HesseWeinberg(3, 3) * v2 / v1 + HesseWeinberg(5, 7) +
                  2 * HesseWeinberg(7, 7) * v2 / v1); // dRem12sq
  parCT.push_back((1.0 / 2.0) * HesseWeinberg(4, 7) +
                  (3.0 / 2.0) * HesseWeinberg(5, 6) +
                  2 * HesseWeinberg(6, 7) * v2 / v1); // dImm12sq
  parCT.push_back(-NablaWeinberg(0));                 // dT1
  parCT.push_back(-NablaWeinberg(1));                 // dT2
  parCT.push_back(-NablaWeinberg(2));                 // dT3
  parCT.push_back(-NablaWeinberg(3));                 // dT4
  parCT.push_back(HesseWeinberg(5, 5) * v1 + HesseWeinberg(5, 7) * v2 -
                  NablaWeinberg(4)); // dT5
  parCT.push_back(HesseWeinberg(5, 6) * v2 +
                  HesseWeinberg(6, 7) * std::pow(v2, 2) / v1 -
                  NablaWeinberg(5)); // dT6
  parCT.push_back(HesseWeinberg(5, 7) * v1 + HesseWeinberg(7, 7) * v2 -
                  NablaWeinberg(6)); // dT7
  parCT.push_back(-HesseWeinberg(5, 6) * v1 - HesseWeinberg(6, 7) * v2 -
                  NablaWeinberg(7)); // dT8
  // End CT Calculation

  return parCT;
}

void C2HDMSympy::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsdone) CalculatePhysicalCouplings();

  std::vector<double> HiggsOrder(NHiggs);
  // Here you have to set the vector HiggsOrder. By telling e.g. HiggsOrder[0] =
  // 5 you always want your 6th lightest particle to be the first particle in
  // the vector (which has the index 5 because they are sorted by mass)

  MatrixXd HiggsRotSort(NHiggs, NHiggs);
  int posMHCS1 = 0, posMHCS2 = 0;
  int posN[3];
  int countposN = 0;
  int posG1 = 0, posG2 = 0, posG0 = 0;
  double testsum             = 0;
  const double ZeroThreshold = 1e-5;
  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
    }
  }
  for (int i = 0; i < 3; i++)
  {
    testsum = std::abs(HiggsRot(i, 0)) + std::abs(HiggsRot(i, 2));
    if (testsum > ZeroThreshold) posG1 = i;
    testsum = std::abs(HiggsRot(i, 1)) + std::abs(HiggsRot(i, 3));
    if (testsum > ZeroThreshold) posG2 = i;
    testsum = std::abs(HiggsRot(i, 5)) + std::abs(HiggsRot(i, 7));
    if (testsum > ZeroThreshold) posG0 = i;
  }
  for (std::size_t i = 3; i < NHiggs; i++)
  {
    testsum = std::abs(HiggsRot(i, 0)) + std::abs(HiggsRot(i, 2));
    if (testsum > ZeroThreshold) posMHCS1 = i;
    testsum = std::abs(HiggsRot(i, 1)) + std::abs(HiggsRot(i, 3));
    if (testsum > ZeroThreshold) posMHCS2 = i;
    testsum = 0;
    for (int k = 4; k < 8; k++)
      testsum += std::abs(HiggsRot(i, k));
    if (testsum > ZeroThreshold)
    {
      posN[countposN] = i;
      countposN++;
    }
  }

  if (false)
  {
    std::vector<double> HiggsMasses;
    double MhUp = 0, MhDown = 0, MSM = 0;
    HiggsMasses = HiggsMassesSquared(vevTree, 0);

    double NeutralHiggs[3];
    for (int i = 0; i < 3; i++)
    {
      NeutralHiggs[i] = HiggsMasses[posN[i]];
    }
    for (int i = 0; i < 3; i++)
    {
      if (std::sqrt(NeutralHiggs[i]) < 126 and std::sqrt(NeutralHiggs[i]) > 124)
        MSM = std::sqrt(NeutralHiggs[i]);
    }
    if (std::sqrt(NeutralHiggs[0]) == MSM)
    {
      MhUp   = std::sqrt(NeutralHiggs[2]);
      MhDown = std::sqrt(NeutralHiggs[1]);
    }
    else if (std::sqrt(NeutralHiggs[1]) == MSM)
    {
      MhUp   = std::sqrt(NeutralHiggs[2]);
      MhDown = std::sqrt(NeutralHiggs[0]);
    }
    else
    {
      MhUp   = std::sqrt(NeutralHiggs[1]);
      MhDown = std::sqrt(NeutralHiggs[0]);
    }

    if (MSM > MhUp)
    {
      double tmp = posN[1];
      posN[1]    = posN[2];
      posN[2]    = tmp;
    }
    if (MSM > MhDown)
    {
      double tmp = posN[0];
      posN[0]    = posN[1];
      posN[1]    = tmp;
    }
  }

  HiggsOrder[0] = posG1;
  HiggsOrder[1] = posG2;
  HiggsOrder[2] = posMHCS1;
  HiggsOrder[3] = posMHCS2;
  HiggsOrder[4] = posG0;
  for (int i = 5; i < 8; i++)
    HiggsOrder[i] = posN[i - 5];

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsRotSort.row(i) = HiggsRot.row(HiggsOrder[i]);
  }

  std::vector<double> TripleDeriv;
  TripleDeriv = WeinbergThirdDerivative();
  std::vector<std::vector<std::vector<double>>> GaugeBasis(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        GaugeBasis[i][j][k] =
            TripleDeriv.at(i + j * NHiggs + k * NHiggs * NHiggs);
      }
    }
  }

  TripleHiggsCorrectionsCWPhysical.resize(NHiggs);
  TripleHiggsCorrectionsTreePhysical.resize(NHiggs);
  TripleHiggsCorrectionsCTPhysical.resize(NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        TripleHiggsCorrectionsCWPhysical[i][j][k]   = 0;
        TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;
        TripleHiggsCorrectionsCTPhysical[i][j][k]   = 0;
        for (std::size_t l = 0; l < NHiggs; l++)
        {
          for (std::size_t m = 0; m < NHiggs; m++)
          {
            for (std::size_t n = 0; n < NHiggs; n++)
            {
              double RotFac =
                  HiggsRotSort(i, l) * HiggsRotSort(j, m) * HiggsRotSort(k, n);
              TripleHiggsCorrectionsCWPhysical[i][j][k] +=
                  RotFac * GaugeBasis[l][m][n];
              TripleHiggsCorrectionsTreePhysical[i][j][k] +=
                  RotFac * LambdaHiggs_3[l][m][n];
              TripleHiggsCorrectionsCTPhysical[i][j][k] +=
                  RotFac * LambdaHiggs_3_CT[l][m][n];
            }
          }
        }
      }
    }
  }
}

void C2HDMSympy::SetCurvatureArrays()
{
  /*
   *  Here you have to set the vectors
   *  Curvature_Higgs_L1,Curvature_Higgs_L2,Curvature_Higgs_L3,Curvature_Higgs_L4
   *  Curvature_Gauge_G2H2
   *  Curvature_Quark_F2H1, Curvature_Lepton_F2H1
   *  as described in the potential in the paper.
   */

  initVectors();
  SetCurvatureDone = true;
  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

  // Begin of Higgs curvature tensors
  {
    Curvature_Higgs_L2.at(0).at(0)             = m11sq;
    Curvature_Higgs_L2.at(0).at(2)             = -Rem12sq;
    Curvature_Higgs_L2.at(0).at(3)             = Imm12sq;
    Curvature_Higgs_L2.at(1).at(1)             = m11sq;
    Curvature_Higgs_L2.at(1).at(2)             = -Imm12sq;
    Curvature_Higgs_L2.at(1).at(3)             = -Rem12sq;
    Curvature_Higgs_L2.at(2).at(0)             = -Rem12sq;
    Curvature_Higgs_L2.at(2).at(1)             = -Imm12sq;
    Curvature_Higgs_L2.at(2).at(2)             = m22sq;
    Curvature_Higgs_L2.at(3).at(0)             = Imm12sq;
    Curvature_Higgs_L2.at(3).at(1)             = -Rem12sq;
    Curvature_Higgs_L2.at(3).at(3)             = m22sq;
    Curvature_Higgs_L2.at(4).at(4)             = m11sq;
    Curvature_Higgs_L2.at(4).at(6)             = -Rem12sq;
    Curvature_Higgs_L2.at(4).at(7)             = Imm12sq;
    Curvature_Higgs_L2.at(5).at(5)             = m11sq;
    Curvature_Higgs_L2.at(5).at(6)             = -Imm12sq;
    Curvature_Higgs_L2.at(5).at(7)             = -Rem12sq;
    Curvature_Higgs_L2.at(6).at(4)             = -Rem12sq;
    Curvature_Higgs_L2.at(6).at(5)             = -Imm12sq;
    Curvature_Higgs_L2.at(6).at(6)             = m22sq;
    Curvature_Higgs_L2.at(7).at(4)             = Imm12sq;
    Curvature_Higgs_L2.at(7).at(5)             = -Rem12sq;
    Curvature_Higgs_L2.at(7).at(7)             = m22sq;
    Curvature_Higgs_L4.at(0).at(0).at(0).at(0) = 3 * lambda1;
    Curvature_Higgs_L4.at(0).at(0).at(0).at(2) = 3 * Relambda6;
    Curvature_Higgs_L4.at(0).at(0).at(0).at(3) = -3 * Imlambda6;
    Curvature_Higgs_L4.at(0).at(0).at(1).at(1) = lambda1;
    Curvature_Higgs_L4.at(0).at(0).at(1).at(2) = Imlambda6;
    Curvature_Higgs_L4.at(0).at(0).at(1).at(3) = Relambda6;
    Curvature_Higgs_L4.at(0).at(0).at(2).at(0) = 3 * Relambda6;
    Curvature_Higgs_L4.at(0).at(0).at(2).at(1) = Imlambda6;
    Curvature_Higgs_L4.at(0).at(0).at(2).at(2) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(0).at(0).at(2).at(3) = -Imlambda5;
    Curvature_Higgs_L4.at(0).at(0).at(3).at(0) = -3 * Imlambda6;
    Curvature_Higgs_L4.at(0).at(0).at(3).at(1) = Relambda6;
    Curvature_Higgs_L4.at(0).at(0).at(3).at(2) = -Imlambda5;
    Curvature_Higgs_L4.at(0).at(0).at(3).at(3) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(0).at(0).at(4).at(4) = lambda1;
    Curvature_Higgs_L4.at(0).at(0).at(4).at(6) = Relambda6;
    Curvature_Higgs_L4.at(0).at(0).at(4).at(7) = -Imlambda6;
    Curvature_Higgs_L4.at(0).at(0).at(5).at(5) = lambda1;
    Curvature_Higgs_L4.at(0).at(0).at(5).at(6) = Imlambda6;
    Curvature_Higgs_L4.at(0).at(0).at(5).at(7) = Relambda6;
    Curvature_Higgs_L4.at(0).at(0).at(6).at(4) = Relambda6;
    Curvature_Higgs_L4.at(0).at(0).at(6).at(5) = Imlambda6;
    Curvature_Higgs_L4.at(0).at(0).at(6).at(6) = lambda3;
    Curvature_Higgs_L4.at(0).at(0).at(7).at(4) = -Imlambda6;
    Curvature_Higgs_L4.at(0).at(0).at(7).at(5) = Relambda6;
    Curvature_Higgs_L4.at(0).at(0).at(7).at(7) = lambda3;
    Curvature_Higgs_L4.at(0).at(1).at(0).at(1) = lambda1;
    Curvature_Higgs_L4.at(0).at(1).at(0).at(2) = Imlambda6;
    Curvature_Higgs_L4.at(0).at(1).at(0).at(3) = Relambda6;
    Curvature_Higgs_L4.at(0).at(1).at(1).at(0) = lambda1;
    Curvature_Higgs_L4.at(0).at(1).at(1).at(2) = Relambda6;
    Curvature_Higgs_L4.at(0).at(1).at(1).at(3) = -Imlambda6;
    Curvature_Higgs_L4.at(0).at(1).at(2).at(0) = Imlambda6;
    Curvature_Higgs_L4.at(0).at(1).at(2).at(1) = Relambda6;
    Curvature_Higgs_L4.at(0).at(1).at(2).at(2) = Imlambda5;
    Curvature_Higgs_L4.at(0).at(1).at(2).at(3) = Relambda5;
    Curvature_Higgs_L4.at(0).at(1).at(3).at(0) = Relambda6;
    Curvature_Higgs_L4.at(0).at(1).at(3).at(1) = -Imlambda6;
    Curvature_Higgs_L4.at(0).at(1).at(3).at(2) = Relambda5;
    Curvature_Higgs_L4.at(0).at(1).at(3).at(3) = -Imlambda5;
    Curvature_Higgs_L4.at(0).at(2).at(0).at(0) = 3 * Relambda6;
    Curvature_Higgs_L4.at(0).at(2).at(0).at(1) = Imlambda6;
    Curvature_Higgs_L4.at(0).at(2).at(0).at(2) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(0).at(2).at(0).at(3) = -Imlambda5;
    Curvature_Higgs_L4.at(0).at(2).at(1).at(0) = Imlambda6;
    Curvature_Higgs_L4.at(0).at(2).at(1).at(1) = Relambda6;
    Curvature_Higgs_L4.at(0).at(2).at(1).at(2) = Imlambda5;
    Curvature_Higgs_L4.at(0).at(2).at(1).at(3) = Relambda5;
    Curvature_Higgs_L4.at(0).at(2).at(2).at(0) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(0).at(2).at(2).at(1) = Imlambda5;
    Curvature_Higgs_L4.at(0).at(2).at(2).at(2) = 3 * Relambda7;
    Curvature_Higgs_L4.at(0).at(2).at(2).at(3) = -Imlambda7;
    Curvature_Higgs_L4.at(0).at(2).at(3).at(0) = -Imlambda5;
    Curvature_Higgs_L4.at(0).at(2).at(3).at(1) = Relambda5;
    Curvature_Higgs_L4.at(0).at(2).at(3).at(2) = -Imlambda7;
    Curvature_Higgs_L4.at(0).at(2).at(3).at(3) = Relambda7;
    Curvature_Higgs_L4.at(0).at(2).at(4).at(4) = Relambda6;
    Curvature_Higgs_L4.at(0).at(2).at(4).at(6) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(2).at(4).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(2).at(5).at(5) = Relambda6;
    Curvature_Higgs_L4.at(0).at(2).at(5).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(0).at(2).at(5).at(7) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(2).at(6).at(4) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(2).at(6).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(0).at(2).at(6).at(6) = Relambda7;
    Curvature_Higgs_L4.at(0).at(2).at(7).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(2).at(7).at(5) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(2).at(7).at(7) = Relambda7;
    Curvature_Higgs_L4.at(0).at(3).at(0).at(0) = -3 * Imlambda6;
    Curvature_Higgs_L4.at(0).at(3).at(0).at(1) = Relambda6;
    Curvature_Higgs_L4.at(0).at(3).at(0).at(2) = -Imlambda5;
    Curvature_Higgs_L4.at(0).at(3).at(0).at(3) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(0).at(3).at(1).at(0) = Relambda6;
    Curvature_Higgs_L4.at(0).at(3).at(1).at(1) = -Imlambda6;
    Curvature_Higgs_L4.at(0).at(3).at(1).at(2) = Relambda5;
    Curvature_Higgs_L4.at(0).at(3).at(1).at(3) = -Imlambda5;
    Curvature_Higgs_L4.at(0).at(3).at(2).at(0) = -Imlambda5;
    Curvature_Higgs_L4.at(0).at(3).at(2).at(1) = Relambda5;
    Curvature_Higgs_L4.at(0).at(3).at(2).at(2) = -Imlambda7;
    Curvature_Higgs_L4.at(0).at(3).at(2).at(3) = Relambda7;
    Curvature_Higgs_L4.at(0).at(3).at(3).at(0) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(0).at(3).at(3).at(1) = -Imlambda5;
    Curvature_Higgs_L4.at(0).at(3).at(3).at(2) = Relambda7;
    Curvature_Higgs_L4.at(0).at(3).at(3).at(3) = -3 * Imlambda7;
    Curvature_Higgs_L4.at(0).at(3).at(4).at(4) = -Imlambda6;
    Curvature_Higgs_L4.at(0).at(3).at(4).at(6) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(3).at(4).at(7) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(3).at(5).at(5) = -Imlambda6;
    Curvature_Higgs_L4.at(0).at(3).at(5).at(6) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(0).at(3).at(5).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(3).at(6).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(3).at(6).at(5) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(0).at(3).at(6).at(6) = -Imlambda7;
    Curvature_Higgs_L4.at(0).at(3).at(7).at(4) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(3).at(7).at(5) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(3).at(7).at(7) = -Imlambda7;
    Curvature_Higgs_L4.at(0).at(4).at(0).at(4) = lambda1;
    Curvature_Higgs_L4.at(0).at(4).at(0).at(6) = Relambda6;
    Curvature_Higgs_L4.at(0).at(4).at(0).at(7) = -Imlambda6;
    Curvature_Higgs_L4.at(0).at(4).at(2).at(4) = Relambda6;
    Curvature_Higgs_L4.at(0).at(4).at(2).at(6) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(4).at(2).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(4).at(3).at(4) = -Imlambda6;
    Curvature_Higgs_L4.at(0).at(4).at(3).at(6) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(4).at(3).at(7) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(4).at(4).at(0) = lambda1;
    Curvature_Higgs_L4.at(0).at(4).at(4).at(2) = Relambda6;
    Curvature_Higgs_L4.at(0).at(4).at(4).at(3) = -Imlambda6;
    Curvature_Higgs_L4.at(0).at(4).at(6).at(0) = Relambda6;
    Curvature_Higgs_L4.at(0).at(4).at(6).at(2) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(4).at(6).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(4).at(7).at(0) = -Imlambda6;
    Curvature_Higgs_L4.at(0).at(4).at(7).at(2) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(4).at(7).at(3) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(5).at(0).at(5) = lambda1;
    Curvature_Higgs_L4.at(0).at(5).at(0).at(6) = Imlambda6;
    Curvature_Higgs_L4.at(0).at(5).at(0).at(7) = Relambda6;
    Curvature_Higgs_L4.at(0).at(5).at(2).at(5) = Relambda6;
    Curvature_Higgs_L4.at(0).at(5).at(2).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(0).at(5).at(2).at(7) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(5).at(3).at(5) = -Imlambda6;
    Curvature_Higgs_L4.at(0).at(5).at(3).at(6) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(0).at(5).at(3).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(5).at(5).at(0) = lambda1;
    Curvature_Higgs_L4.at(0).at(5).at(5).at(2) = Relambda6;
    Curvature_Higgs_L4.at(0).at(5).at(5).at(3) = -Imlambda6;
    Curvature_Higgs_L4.at(0).at(5).at(6).at(0) = Imlambda6;
    Curvature_Higgs_L4.at(0).at(5).at(6).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(0).at(5).at(6).at(3) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(0).at(5).at(7).at(0) = Relambda6;
    Curvature_Higgs_L4.at(0).at(5).at(7).at(2) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(5).at(7).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(6).at(0).at(4) = Relambda6;
    Curvature_Higgs_L4.at(0).at(6).at(0).at(5) = Imlambda6;
    Curvature_Higgs_L4.at(0).at(6).at(0).at(6) = lambda3;
    Curvature_Higgs_L4.at(0).at(6).at(2).at(4) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(6).at(2).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(0).at(6).at(2).at(6) = Relambda7;
    Curvature_Higgs_L4.at(0).at(6).at(3).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(6).at(3).at(5) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(0).at(6).at(3).at(6) = -Imlambda7;
    Curvature_Higgs_L4.at(0).at(6).at(4).at(0) = Relambda6;
    Curvature_Higgs_L4.at(0).at(6).at(4).at(2) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(6).at(4).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(6).at(5).at(0) = Imlambda6;
    Curvature_Higgs_L4.at(0).at(6).at(5).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(0).at(6).at(5).at(3) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(0).at(6).at(6).at(0) = lambda3;
    Curvature_Higgs_L4.at(0).at(6).at(6).at(2) = Relambda7;
    Curvature_Higgs_L4.at(0).at(6).at(6).at(3) = -Imlambda7;
    Curvature_Higgs_L4.at(0).at(7).at(0).at(4) = -Imlambda6;
    Curvature_Higgs_L4.at(0).at(7).at(0).at(5) = Relambda6;
    Curvature_Higgs_L4.at(0).at(7).at(0).at(7) = lambda3;
    Curvature_Higgs_L4.at(0).at(7).at(2).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(7).at(2).at(5) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(7).at(2).at(7) = Relambda7;
    Curvature_Higgs_L4.at(0).at(7).at(3).at(4) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(7).at(3).at(5) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(7).at(3).at(7) = -Imlambda7;
    Curvature_Higgs_L4.at(0).at(7).at(4).at(0) = -Imlambda6;
    Curvature_Higgs_L4.at(0).at(7).at(4).at(2) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(7).at(4).at(3) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(7).at(5).at(0) = Relambda6;
    Curvature_Higgs_L4.at(0).at(7).at(5).at(2) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(0).at(7).at(5).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(0).at(7).at(7).at(0) = lambda3;
    Curvature_Higgs_L4.at(0).at(7).at(7).at(2) = Relambda7;
    Curvature_Higgs_L4.at(0).at(7).at(7).at(3) = -Imlambda7;
    Curvature_Higgs_L4.at(1).at(0).at(0).at(1) = lambda1;
    Curvature_Higgs_L4.at(1).at(0).at(0).at(2) = Imlambda6;
    Curvature_Higgs_L4.at(1).at(0).at(0).at(3) = Relambda6;
    Curvature_Higgs_L4.at(1).at(0).at(1).at(0) = lambda1;
    Curvature_Higgs_L4.at(1).at(0).at(1).at(2) = Relambda6;
    Curvature_Higgs_L4.at(1).at(0).at(1).at(3) = -Imlambda6;
    Curvature_Higgs_L4.at(1).at(0).at(2).at(0) = Imlambda6;
    Curvature_Higgs_L4.at(1).at(0).at(2).at(1) = Relambda6;
    Curvature_Higgs_L4.at(1).at(0).at(2).at(2) = Imlambda5;
    Curvature_Higgs_L4.at(1).at(0).at(2).at(3) = Relambda5;
    Curvature_Higgs_L4.at(1).at(0).at(3).at(0) = Relambda6;
    Curvature_Higgs_L4.at(1).at(0).at(3).at(1) = -Imlambda6;
    Curvature_Higgs_L4.at(1).at(0).at(3).at(2) = Relambda5;
    Curvature_Higgs_L4.at(1).at(0).at(3).at(3) = -Imlambda5;
    Curvature_Higgs_L4.at(1).at(1).at(0).at(0) = lambda1;
    Curvature_Higgs_L4.at(1).at(1).at(0).at(2) = Relambda6;
    Curvature_Higgs_L4.at(1).at(1).at(0).at(3) = -Imlambda6;
    Curvature_Higgs_L4.at(1).at(1).at(1).at(1) = 3 * lambda1;
    Curvature_Higgs_L4.at(1).at(1).at(1).at(2) = 3 * Imlambda6;
    Curvature_Higgs_L4.at(1).at(1).at(1).at(3) = 3 * Relambda6;
    Curvature_Higgs_L4.at(1).at(1).at(2).at(0) = Relambda6;
    Curvature_Higgs_L4.at(1).at(1).at(2).at(1) = 3 * Imlambda6;
    Curvature_Higgs_L4.at(1).at(1).at(2).at(2) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(1).at(1).at(2).at(3) = Imlambda5;
    Curvature_Higgs_L4.at(1).at(1).at(3).at(0) = -Imlambda6;
    Curvature_Higgs_L4.at(1).at(1).at(3).at(1) = 3 * Relambda6;
    Curvature_Higgs_L4.at(1).at(1).at(3).at(2) = Imlambda5;
    Curvature_Higgs_L4.at(1).at(1).at(3).at(3) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(1).at(1).at(4).at(4) = lambda1;
    Curvature_Higgs_L4.at(1).at(1).at(4).at(6) = Relambda6;
    Curvature_Higgs_L4.at(1).at(1).at(4).at(7) = -Imlambda6;
    Curvature_Higgs_L4.at(1).at(1).at(5).at(5) = lambda1;
    Curvature_Higgs_L4.at(1).at(1).at(5).at(6) = Imlambda6;
    Curvature_Higgs_L4.at(1).at(1).at(5).at(7) = Relambda6;
    Curvature_Higgs_L4.at(1).at(1).at(6).at(4) = Relambda6;
    Curvature_Higgs_L4.at(1).at(1).at(6).at(5) = Imlambda6;
    Curvature_Higgs_L4.at(1).at(1).at(6).at(6) = lambda3;
    Curvature_Higgs_L4.at(1).at(1).at(7).at(4) = -Imlambda6;
    Curvature_Higgs_L4.at(1).at(1).at(7).at(5) = Relambda6;
    Curvature_Higgs_L4.at(1).at(1).at(7).at(7) = lambda3;
    Curvature_Higgs_L4.at(1).at(2).at(0).at(0) = Imlambda6;
    Curvature_Higgs_L4.at(1).at(2).at(0).at(1) = Relambda6;
    Curvature_Higgs_L4.at(1).at(2).at(0).at(2) = Imlambda5;
    Curvature_Higgs_L4.at(1).at(2).at(0).at(3) = Relambda5;
    Curvature_Higgs_L4.at(1).at(2).at(1).at(0) = Relambda6;
    Curvature_Higgs_L4.at(1).at(2).at(1).at(1) = 3 * Imlambda6;
    Curvature_Higgs_L4.at(1).at(2).at(1).at(2) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(1).at(2).at(1).at(3) = Imlambda5;
    Curvature_Higgs_L4.at(1).at(2).at(2).at(0) = Imlambda5;
    Curvature_Higgs_L4.at(1).at(2).at(2).at(1) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(1).at(2).at(2).at(2) = 3 * Imlambda7;
    Curvature_Higgs_L4.at(1).at(2).at(2).at(3) = Relambda7;
    Curvature_Higgs_L4.at(1).at(2).at(3).at(0) = Relambda5;
    Curvature_Higgs_L4.at(1).at(2).at(3).at(1) = Imlambda5;
    Curvature_Higgs_L4.at(1).at(2).at(3).at(2) = Relambda7;
    Curvature_Higgs_L4.at(1).at(2).at(3).at(3) = Imlambda7;
    Curvature_Higgs_L4.at(1).at(2).at(4).at(4) = Imlambda6;
    Curvature_Higgs_L4.at(1).at(2).at(4).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(2).at(4).at(7) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(1).at(2).at(5).at(5) = Imlambda6;
    Curvature_Higgs_L4.at(1).at(2).at(5).at(6) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(2).at(5).at(7) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(2).at(6).at(4) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(2).at(6).at(5) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(2).at(6).at(6) = Imlambda7;
    Curvature_Higgs_L4.at(1).at(2).at(7).at(4) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(1).at(2).at(7).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(2).at(7).at(7) = Imlambda7;
    Curvature_Higgs_L4.at(1).at(3).at(0).at(0) = Relambda6;
    Curvature_Higgs_L4.at(1).at(3).at(0).at(1) = -Imlambda6;
    Curvature_Higgs_L4.at(1).at(3).at(0).at(2) = Relambda5;
    Curvature_Higgs_L4.at(1).at(3).at(0).at(3) = -Imlambda5;
    Curvature_Higgs_L4.at(1).at(3).at(1).at(0) = -Imlambda6;
    Curvature_Higgs_L4.at(1).at(3).at(1).at(1) = 3 * Relambda6;
    Curvature_Higgs_L4.at(1).at(3).at(1).at(2) = Imlambda5;
    Curvature_Higgs_L4.at(1).at(3).at(1).at(3) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(1).at(3).at(2).at(0) = Relambda5;
    Curvature_Higgs_L4.at(1).at(3).at(2).at(1) = Imlambda5;
    Curvature_Higgs_L4.at(1).at(3).at(2).at(2) = Relambda7;
    Curvature_Higgs_L4.at(1).at(3).at(2).at(3) = Imlambda7;
    Curvature_Higgs_L4.at(1).at(3).at(3).at(0) = -Imlambda5;
    Curvature_Higgs_L4.at(1).at(3).at(3).at(1) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(1).at(3).at(3).at(2) = Imlambda7;
    Curvature_Higgs_L4.at(1).at(3).at(3).at(3) = 3 * Relambda7;
    Curvature_Higgs_L4.at(1).at(3).at(4).at(4) = Relambda6;
    Curvature_Higgs_L4.at(1).at(3).at(4).at(6) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(3).at(4).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(1).at(3).at(5).at(5) = Relambda6;
    Curvature_Higgs_L4.at(1).at(3).at(5).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(3).at(5).at(7) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(3).at(6).at(4) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(3).at(6).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(3).at(6).at(6) = Relambda7;
    Curvature_Higgs_L4.at(1).at(3).at(7).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(1).at(3).at(7).at(5) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(3).at(7).at(7) = Relambda7;
    Curvature_Higgs_L4.at(1).at(4).at(1).at(4) = lambda1;
    Curvature_Higgs_L4.at(1).at(4).at(1).at(6) = Relambda6;
    Curvature_Higgs_L4.at(1).at(4).at(1).at(7) = -Imlambda6;
    Curvature_Higgs_L4.at(1).at(4).at(2).at(4) = Imlambda6;
    Curvature_Higgs_L4.at(1).at(4).at(2).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(4).at(2).at(7) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(1).at(4).at(3).at(4) = Relambda6;
    Curvature_Higgs_L4.at(1).at(4).at(3).at(6) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(4).at(3).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(1).at(4).at(4).at(1) = lambda1;
    Curvature_Higgs_L4.at(1).at(4).at(4).at(2) = Imlambda6;
    Curvature_Higgs_L4.at(1).at(4).at(4).at(3) = Relambda6;
    Curvature_Higgs_L4.at(1).at(4).at(6).at(1) = Relambda6;
    Curvature_Higgs_L4.at(1).at(4).at(6).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(4).at(6).at(3) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(4).at(7).at(1) = -Imlambda6;
    Curvature_Higgs_L4.at(1).at(4).at(7).at(2) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(1).at(4).at(7).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(1).at(5).at(1).at(5) = lambda1;
    Curvature_Higgs_L4.at(1).at(5).at(1).at(6) = Imlambda6;
    Curvature_Higgs_L4.at(1).at(5).at(1).at(7) = Relambda6;
    Curvature_Higgs_L4.at(1).at(5).at(2).at(5) = Imlambda6;
    Curvature_Higgs_L4.at(1).at(5).at(2).at(6) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(5).at(2).at(7) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(5).at(3).at(5) = Relambda6;
    Curvature_Higgs_L4.at(1).at(5).at(3).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(5).at(3).at(7) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(5).at(5).at(1) = lambda1;
    Curvature_Higgs_L4.at(1).at(5).at(5).at(2) = Imlambda6;
    Curvature_Higgs_L4.at(1).at(5).at(5).at(3) = Relambda6;
    Curvature_Higgs_L4.at(1).at(5).at(6).at(1) = Imlambda6;
    Curvature_Higgs_L4.at(1).at(5).at(6).at(2) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(5).at(6).at(3) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(5).at(7).at(1) = Relambda6;
    Curvature_Higgs_L4.at(1).at(5).at(7).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(5).at(7).at(3) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(6).at(1).at(4) = Relambda6;
    Curvature_Higgs_L4.at(1).at(6).at(1).at(5) = Imlambda6;
    Curvature_Higgs_L4.at(1).at(6).at(1).at(6) = lambda3;
    Curvature_Higgs_L4.at(1).at(6).at(2).at(4) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(6).at(2).at(5) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(6).at(2).at(6) = Imlambda7;
    Curvature_Higgs_L4.at(1).at(6).at(3).at(4) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(6).at(3).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(6).at(3).at(6) = Relambda7;
    Curvature_Higgs_L4.at(1).at(6).at(4).at(1) = Relambda6;
    Curvature_Higgs_L4.at(1).at(6).at(4).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(6).at(4).at(3) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(6).at(5).at(1) = Imlambda6;
    Curvature_Higgs_L4.at(1).at(6).at(5).at(2) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(6).at(5).at(3) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(6).at(6).at(1) = lambda3;
    Curvature_Higgs_L4.at(1).at(6).at(6).at(2) = Imlambda7;
    Curvature_Higgs_L4.at(1).at(6).at(6).at(3) = Relambda7;
    Curvature_Higgs_L4.at(1).at(7).at(1).at(4) = -Imlambda6;
    Curvature_Higgs_L4.at(1).at(7).at(1).at(5) = Relambda6;
    Curvature_Higgs_L4.at(1).at(7).at(1).at(7) = lambda3;
    Curvature_Higgs_L4.at(1).at(7).at(2).at(4) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(1).at(7).at(2).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(7).at(2).at(7) = Imlambda7;
    Curvature_Higgs_L4.at(1).at(7).at(3).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(1).at(7).at(3).at(5) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(7).at(3).at(7) = Relambda7;
    Curvature_Higgs_L4.at(1).at(7).at(4).at(1) = -Imlambda6;
    Curvature_Higgs_L4.at(1).at(7).at(4).at(2) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(1).at(7).at(4).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(1).at(7).at(5).at(1) = Relambda6;
    Curvature_Higgs_L4.at(1).at(7).at(5).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(1).at(7).at(5).at(3) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(1).at(7).at(7).at(1) = lambda3;
    Curvature_Higgs_L4.at(1).at(7).at(7).at(2) = Imlambda7;
    Curvature_Higgs_L4.at(1).at(7).at(7).at(3) = Relambda7;
    Curvature_Higgs_L4.at(2).at(0).at(0).at(0) = 3 * Relambda6;
    Curvature_Higgs_L4.at(2).at(0).at(0).at(1) = Imlambda6;
    Curvature_Higgs_L4.at(2).at(0).at(0).at(2) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(2).at(0).at(0).at(3) = -Imlambda5;
    Curvature_Higgs_L4.at(2).at(0).at(1).at(0) = Imlambda6;
    Curvature_Higgs_L4.at(2).at(0).at(1).at(1) = Relambda6;
    Curvature_Higgs_L4.at(2).at(0).at(1).at(2) = Imlambda5;
    Curvature_Higgs_L4.at(2).at(0).at(1).at(3) = Relambda5;
    Curvature_Higgs_L4.at(2).at(0).at(2).at(0) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(2).at(0).at(2).at(1) = Imlambda5;
    Curvature_Higgs_L4.at(2).at(0).at(2).at(2) = 3 * Relambda7;
    Curvature_Higgs_L4.at(2).at(0).at(2).at(3) = -Imlambda7;
    Curvature_Higgs_L4.at(2).at(0).at(3).at(0) = -Imlambda5;
    Curvature_Higgs_L4.at(2).at(0).at(3).at(1) = Relambda5;
    Curvature_Higgs_L4.at(2).at(0).at(3).at(2) = -Imlambda7;
    Curvature_Higgs_L4.at(2).at(0).at(3).at(3) = Relambda7;
    Curvature_Higgs_L4.at(2).at(0).at(4).at(4) = Relambda6;
    Curvature_Higgs_L4.at(2).at(0).at(4).at(6) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(0).at(4).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(2).at(0).at(5).at(5) = Relambda6;
    Curvature_Higgs_L4.at(2).at(0).at(5).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(0).at(5).at(7) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(0).at(6).at(4) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(0).at(6).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(0).at(6).at(6) = Relambda7;
    Curvature_Higgs_L4.at(2).at(0).at(7).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(2).at(0).at(7).at(5) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(0).at(7).at(7) = Relambda7;
    Curvature_Higgs_L4.at(2).at(1).at(0).at(0) = Imlambda6;
    Curvature_Higgs_L4.at(2).at(1).at(0).at(1) = Relambda6;
    Curvature_Higgs_L4.at(2).at(1).at(0).at(2) = Imlambda5;
    Curvature_Higgs_L4.at(2).at(1).at(0).at(3) = Relambda5;
    Curvature_Higgs_L4.at(2).at(1).at(1).at(0) = Relambda6;
    Curvature_Higgs_L4.at(2).at(1).at(1).at(1) = 3 * Imlambda6;
    Curvature_Higgs_L4.at(2).at(1).at(1).at(2) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(2).at(1).at(1).at(3) = Imlambda5;
    Curvature_Higgs_L4.at(2).at(1).at(2).at(0) = Imlambda5;
    Curvature_Higgs_L4.at(2).at(1).at(2).at(1) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(2).at(1).at(2).at(2) = 3 * Imlambda7;
    Curvature_Higgs_L4.at(2).at(1).at(2).at(3) = Relambda7;
    Curvature_Higgs_L4.at(2).at(1).at(3).at(0) = Relambda5;
    Curvature_Higgs_L4.at(2).at(1).at(3).at(1) = Imlambda5;
    Curvature_Higgs_L4.at(2).at(1).at(3).at(2) = Relambda7;
    Curvature_Higgs_L4.at(2).at(1).at(3).at(3) = Imlambda7;
    Curvature_Higgs_L4.at(2).at(1).at(4).at(4) = Imlambda6;
    Curvature_Higgs_L4.at(2).at(1).at(4).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(1).at(4).at(7) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(2).at(1).at(5).at(5) = Imlambda6;
    Curvature_Higgs_L4.at(2).at(1).at(5).at(6) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(1).at(5).at(7) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(1).at(6).at(4) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(1).at(6).at(5) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(1).at(6).at(6) = Imlambda7;
    Curvature_Higgs_L4.at(2).at(1).at(7).at(4) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(2).at(1).at(7).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(1).at(7).at(7) = Imlambda7;
    Curvature_Higgs_L4.at(2).at(2).at(0).at(0) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(2).at(2).at(0).at(1) = Imlambda5;
    Curvature_Higgs_L4.at(2).at(2).at(0).at(2) = 3 * Relambda7;
    Curvature_Higgs_L4.at(2).at(2).at(0).at(3) = -Imlambda7;
    Curvature_Higgs_L4.at(2).at(2).at(1).at(0) = Imlambda5;
    Curvature_Higgs_L4.at(2).at(2).at(1).at(1) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(2).at(2).at(1).at(2) = 3 * Imlambda7;
    Curvature_Higgs_L4.at(2).at(2).at(1).at(3) = Relambda7;
    Curvature_Higgs_L4.at(2).at(2).at(2).at(0) = 3 * Relambda7;
    Curvature_Higgs_L4.at(2).at(2).at(2).at(1) = 3 * Imlambda7;
    Curvature_Higgs_L4.at(2).at(2).at(2).at(2) = 3 * lambda2;
    Curvature_Higgs_L4.at(2).at(2).at(3).at(0) = -Imlambda7;
    Curvature_Higgs_L4.at(2).at(2).at(3).at(1) = Relambda7;
    Curvature_Higgs_L4.at(2).at(2).at(3).at(3) = lambda2;
    Curvature_Higgs_L4.at(2).at(2).at(4).at(4) = lambda3;
    Curvature_Higgs_L4.at(2).at(2).at(4).at(6) = Relambda7;
    Curvature_Higgs_L4.at(2).at(2).at(4).at(7) = -Imlambda7;
    Curvature_Higgs_L4.at(2).at(2).at(5).at(5) = lambda3;
    Curvature_Higgs_L4.at(2).at(2).at(5).at(6) = Imlambda7;
    Curvature_Higgs_L4.at(2).at(2).at(5).at(7) = Relambda7;
    Curvature_Higgs_L4.at(2).at(2).at(6).at(4) = Relambda7;
    Curvature_Higgs_L4.at(2).at(2).at(6).at(5) = Imlambda7;
    Curvature_Higgs_L4.at(2).at(2).at(6).at(6) = lambda2;
    Curvature_Higgs_L4.at(2).at(2).at(7).at(4) = -Imlambda7;
    Curvature_Higgs_L4.at(2).at(2).at(7).at(5) = Relambda7;
    Curvature_Higgs_L4.at(2).at(2).at(7).at(7) = lambda2;
    Curvature_Higgs_L4.at(2).at(3).at(0).at(0) = -Imlambda5;
    Curvature_Higgs_L4.at(2).at(3).at(0).at(1) = Relambda5;
    Curvature_Higgs_L4.at(2).at(3).at(0).at(2) = -Imlambda7;
    Curvature_Higgs_L4.at(2).at(3).at(0).at(3) = Relambda7;
    Curvature_Higgs_L4.at(2).at(3).at(1).at(0) = Relambda5;
    Curvature_Higgs_L4.at(2).at(3).at(1).at(1) = Imlambda5;
    Curvature_Higgs_L4.at(2).at(3).at(1).at(2) = Relambda7;
    Curvature_Higgs_L4.at(2).at(3).at(1).at(3) = Imlambda7;
    Curvature_Higgs_L4.at(2).at(3).at(2).at(0) = -Imlambda7;
    Curvature_Higgs_L4.at(2).at(3).at(2).at(1) = Relambda7;
    Curvature_Higgs_L4.at(2).at(3).at(2).at(3) = lambda2;
    Curvature_Higgs_L4.at(2).at(3).at(3).at(0) = Relambda7;
    Curvature_Higgs_L4.at(2).at(3).at(3).at(1) = Imlambda7;
    Curvature_Higgs_L4.at(2).at(3).at(3).at(2) = lambda2;
    Curvature_Higgs_L4.at(2).at(4).at(0).at(4) = Relambda6;
    Curvature_Higgs_L4.at(2).at(4).at(0).at(6) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(4).at(0).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(2).at(4).at(1).at(4) = Imlambda6;
    Curvature_Higgs_L4.at(2).at(4).at(1).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(4).at(1).at(7) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(2).at(4).at(2).at(4) = lambda3;
    Curvature_Higgs_L4.at(2).at(4).at(2).at(6) = Relambda7;
    Curvature_Higgs_L4.at(2).at(4).at(2).at(7) = -Imlambda7;
    Curvature_Higgs_L4.at(2).at(4).at(4).at(0) = Relambda6;
    Curvature_Higgs_L4.at(2).at(4).at(4).at(1) = Imlambda6;
    Curvature_Higgs_L4.at(2).at(4).at(4).at(2) = lambda3;
    Curvature_Higgs_L4.at(2).at(4).at(6).at(0) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(4).at(6).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(4).at(6).at(2) = Relambda7;
    Curvature_Higgs_L4.at(2).at(4).at(7).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(2).at(4).at(7).at(1) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(2).at(4).at(7).at(2) = -Imlambda7;
    Curvature_Higgs_L4.at(2).at(5).at(0).at(5) = Relambda6;
    Curvature_Higgs_L4.at(2).at(5).at(0).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(5).at(0).at(7) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(5).at(1).at(5) = Imlambda6;
    Curvature_Higgs_L4.at(2).at(5).at(1).at(6) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(5).at(1).at(7) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(5).at(2).at(5) = lambda3;
    Curvature_Higgs_L4.at(2).at(5).at(2).at(6) = Imlambda7;
    Curvature_Higgs_L4.at(2).at(5).at(2).at(7) = Relambda7;
    Curvature_Higgs_L4.at(2).at(5).at(5).at(0) = Relambda6;
    Curvature_Higgs_L4.at(2).at(5).at(5).at(1) = Imlambda6;
    Curvature_Higgs_L4.at(2).at(5).at(5).at(2) = lambda3;
    Curvature_Higgs_L4.at(2).at(5).at(6).at(0) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(5).at(6).at(1) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(5).at(6).at(2) = Imlambda7;
    Curvature_Higgs_L4.at(2).at(5).at(7).at(0) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(5).at(7).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(5).at(7).at(2) = Relambda7;
    Curvature_Higgs_L4.at(2).at(6).at(0).at(4) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(6).at(0).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(6).at(0).at(6) = Relambda7;
    Curvature_Higgs_L4.at(2).at(6).at(1).at(4) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(6).at(1).at(5) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(6).at(1).at(6) = Imlambda7;
    Curvature_Higgs_L4.at(2).at(6).at(2).at(4) = Relambda7;
    Curvature_Higgs_L4.at(2).at(6).at(2).at(5) = Imlambda7;
    Curvature_Higgs_L4.at(2).at(6).at(2).at(6) = lambda2;
    Curvature_Higgs_L4.at(2).at(6).at(4).at(0) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(6).at(4).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(6).at(4).at(2) = Relambda7;
    Curvature_Higgs_L4.at(2).at(6).at(5).at(0) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(6).at(5).at(1) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(6).at(5).at(2) = Imlambda7;
    Curvature_Higgs_L4.at(2).at(6).at(6).at(0) = Relambda7;
    Curvature_Higgs_L4.at(2).at(6).at(6).at(1) = Imlambda7;
    Curvature_Higgs_L4.at(2).at(6).at(6).at(2) = lambda2;
    Curvature_Higgs_L4.at(2).at(7).at(0).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(2).at(7).at(0).at(5) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(7).at(0).at(7) = Relambda7;
    Curvature_Higgs_L4.at(2).at(7).at(1).at(4) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(2).at(7).at(1).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(7).at(1).at(7) = Imlambda7;
    Curvature_Higgs_L4.at(2).at(7).at(2).at(4) = -Imlambda7;
    Curvature_Higgs_L4.at(2).at(7).at(2).at(5) = Relambda7;
    Curvature_Higgs_L4.at(2).at(7).at(2).at(7) = lambda2;
    Curvature_Higgs_L4.at(2).at(7).at(4).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(2).at(7).at(4).at(1) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(2).at(7).at(4).at(2) = -Imlambda7;
    Curvature_Higgs_L4.at(2).at(7).at(5).at(0) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(2).at(7).at(5).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(2).at(7).at(5).at(2) = Relambda7;
    Curvature_Higgs_L4.at(2).at(7).at(7).at(0) = Relambda7;
    Curvature_Higgs_L4.at(2).at(7).at(7).at(1) = Imlambda7;
    Curvature_Higgs_L4.at(2).at(7).at(7).at(2) = lambda2;
    Curvature_Higgs_L4.at(3).at(0).at(0).at(0) = -3 * Imlambda6;
    Curvature_Higgs_L4.at(3).at(0).at(0).at(1) = Relambda6;
    Curvature_Higgs_L4.at(3).at(0).at(0).at(2) = -Imlambda5;
    Curvature_Higgs_L4.at(3).at(0).at(0).at(3) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(3).at(0).at(1).at(0) = Relambda6;
    Curvature_Higgs_L4.at(3).at(0).at(1).at(1) = -Imlambda6;
    Curvature_Higgs_L4.at(3).at(0).at(1).at(2) = Relambda5;
    Curvature_Higgs_L4.at(3).at(0).at(1).at(3) = -Imlambda5;
    Curvature_Higgs_L4.at(3).at(0).at(2).at(0) = -Imlambda5;
    Curvature_Higgs_L4.at(3).at(0).at(2).at(1) = Relambda5;
    Curvature_Higgs_L4.at(3).at(0).at(2).at(2) = -Imlambda7;
    Curvature_Higgs_L4.at(3).at(0).at(2).at(3) = Relambda7;
    Curvature_Higgs_L4.at(3).at(0).at(3).at(0) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(3).at(0).at(3).at(1) = -Imlambda5;
    Curvature_Higgs_L4.at(3).at(0).at(3).at(2) = Relambda7;
    Curvature_Higgs_L4.at(3).at(0).at(3).at(3) = -3 * Imlambda7;
    Curvature_Higgs_L4.at(3).at(0).at(4).at(4) = -Imlambda6;
    Curvature_Higgs_L4.at(3).at(0).at(4).at(6) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(0).at(4).at(7) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(0).at(5).at(5) = -Imlambda6;
    Curvature_Higgs_L4.at(3).at(0).at(5).at(6) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(3).at(0).at(5).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(0).at(6).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(0).at(6).at(5) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(3).at(0).at(6).at(6) = -Imlambda7;
    Curvature_Higgs_L4.at(3).at(0).at(7).at(4) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(0).at(7).at(5) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(0).at(7).at(7) = -Imlambda7;
    Curvature_Higgs_L4.at(3).at(1).at(0).at(0) = Relambda6;
    Curvature_Higgs_L4.at(3).at(1).at(0).at(1) = -Imlambda6;
    Curvature_Higgs_L4.at(3).at(1).at(0).at(2) = Relambda5;
    Curvature_Higgs_L4.at(3).at(1).at(0).at(3) = -Imlambda5;
    Curvature_Higgs_L4.at(3).at(1).at(1).at(0) = -Imlambda6;
    Curvature_Higgs_L4.at(3).at(1).at(1).at(1) = 3 * Relambda6;
    Curvature_Higgs_L4.at(3).at(1).at(1).at(2) = Imlambda5;
    Curvature_Higgs_L4.at(3).at(1).at(1).at(3) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(3).at(1).at(2).at(0) = Relambda5;
    Curvature_Higgs_L4.at(3).at(1).at(2).at(1) = Imlambda5;
    Curvature_Higgs_L4.at(3).at(1).at(2).at(2) = Relambda7;
    Curvature_Higgs_L4.at(3).at(1).at(2).at(3) = Imlambda7;
    Curvature_Higgs_L4.at(3).at(1).at(3).at(0) = -Imlambda5;
    Curvature_Higgs_L4.at(3).at(1).at(3).at(1) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(3).at(1).at(3).at(2) = Imlambda7;
    Curvature_Higgs_L4.at(3).at(1).at(3).at(3) = 3 * Relambda7;
    Curvature_Higgs_L4.at(3).at(1).at(4).at(4) = Relambda6;
    Curvature_Higgs_L4.at(3).at(1).at(4).at(6) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(1).at(4).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(1).at(5).at(5) = Relambda6;
    Curvature_Higgs_L4.at(3).at(1).at(5).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(3).at(1).at(5).at(7) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(1).at(6).at(4) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(1).at(6).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(3).at(1).at(6).at(6) = Relambda7;
    Curvature_Higgs_L4.at(3).at(1).at(7).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(1).at(7).at(5) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(1).at(7).at(7) = Relambda7;
    Curvature_Higgs_L4.at(3).at(2).at(0).at(0) = -Imlambda5;
    Curvature_Higgs_L4.at(3).at(2).at(0).at(1) = Relambda5;
    Curvature_Higgs_L4.at(3).at(2).at(0).at(2) = -Imlambda7;
    Curvature_Higgs_L4.at(3).at(2).at(0).at(3) = Relambda7;
    Curvature_Higgs_L4.at(3).at(2).at(1).at(0) = Relambda5;
    Curvature_Higgs_L4.at(3).at(2).at(1).at(1) = Imlambda5;
    Curvature_Higgs_L4.at(3).at(2).at(1).at(2) = Relambda7;
    Curvature_Higgs_L4.at(3).at(2).at(1).at(3) = Imlambda7;
    Curvature_Higgs_L4.at(3).at(2).at(2).at(0) = -Imlambda7;
    Curvature_Higgs_L4.at(3).at(2).at(2).at(1) = Relambda7;
    Curvature_Higgs_L4.at(3).at(2).at(2).at(3) = lambda2;
    Curvature_Higgs_L4.at(3).at(2).at(3).at(0) = Relambda7;
    Curvature_Higgs_L4.at(3).at(2).at(3).at(1) = Imlambda7;
    Curvature_Higgs_L4.at(3).at(2).at(3).at(2) = lambda2;
    Curvature_Higgs_L4.at(3).at(3).at(0).at(0) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(3).at(3).at(0).at(1) = -Imlambda5;
    Curvature_Higgs_L4.at(3).at(3).at(0).at(2) = Relambda7;
    Curvature_Higgs_L4.at(3).at(3).at(0).at(3) = -3 * Imlambda7;
    Curvature_Higgs_L4.at(3).at(3).at(1).at(0) = -Imlambda5;
    Curvature_Higgs_L4.at(3).at(3).at(1).at(1) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(3).at(3).at(1).at(2) = Imlambda7;
    Curvature_Higgs_L4.at(3).at(3).at(1).at(3) = 3 * Relambda7;
    Curvature_Higgs_L4.at(3).at(3).at(2).at(0) = Relambda7;
    Curvature_Higgs_L4.at(3).at(3).at(2).at(1) = Imlambda7;
    Curvature_Higgs_L4.at(3).at(3).at(2).at(2) = lambda2;
    Curvature_Higgs_L4.at(3).at(3).at(3).at(0) = -3 * Imlambda7;
    Curvature_Higgs_L4.at(3).at(3).at(3).at(1) = 3 * Relambda7;
    Curvature_Higgs_L4.at(3).at(3).at(3).at(3) = 3 * lambda2;
    Curvature_Higgs_L4.at(3).at(3).at(4).at(4) = lambda3;
    Curvature_Higgs_L4.at(3).at(3).at(4).at(6) = Relambda7;
    Curvature_Higgs_L4.at(3).at(3).at(4).at(7) = -Imlambda7;
    Curvature_Higgs_L4.at(3).at(3).at(5).at(5) = lambda3;
    Curvature_Higgs_L4.at(3).at(3).at(5).at(6) = Imlambda7;
    Curvature_Higgs_L4.at(3).at(3).at(5).at(7) = Relambda7;
    Curvature_Higgs_L4.at(3).at(3).at(6).at(4) = Relambda7;
    Curvature_Higgs_L4.at(3).at(3).at(6).at(5) = Imlambda7;
    Curvature_Higgs_L4.at(3).at(3).at(6).at(6) = lambda2;
    Curvature_Higgs_L4.at(3).at(3).at(7).at(4) = -Imlambda7;
    Curvature_Higgs_L4.at(3).at(3).at(7).at(5) = Relambda7;
    Curvature_Higgs_L4.at(3).at(3).at(7).at(7) = lambda2;
    Curvature_Higgs_L4.at(3).at(4).at(0).at(4) = -Imlambda6;
    Curvature_Higgs_L4.at(3).at(4).at(0).at(6) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(4).at(0).at(7) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(4).at(1).at(4) = Relambda6;
    Curvature_Higgs_L4.at(3).at(4).at(1).at(6) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(4).at(1).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(4).at(3).at(4) = lambda3;
    Curvature_Higgs_L4.at(3).at(4).at(3).at(6) = Relambda7;
    Curvature_Higgs_L4.at(3).at(4).at(3).at(7) = -Imlambda7;
    Curvature_Higgs_L4.at(3).at(4).at(4).at(0) = -Imlambda6;
    Curvature_Higgs_L4.at(3).at(4).at(4).at(1) = Relambda6;
    Curvature_Higgs_L4.at(3).at(4).at(4).at(3) = lambda3;
    Curvature_Higgs_L4.at(3).at(4).at(6).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(4).at(6).at(1) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(4).at(6).at(3) = Relambda7;
    Curvature_Higgs_L4.at(3).at(4).at(7).at(0) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(4).at(7).at(1) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(4).at(7).at(3) = -Imlambda7;
    Curvature_Higgs_L4.at(3).at(5).at(0).at(5) = -Imlambda6;
    Curvature_Higgs_L4.at(3).at(5).at(0).at(6) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(3).at(5).at(0).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(5).at(1).at(5) = Relambda6;
    Curvature_Higgs_L4.at(3).at(5).at(1).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(3).at(5).at(1).at(7) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(5).at(3).at(5) = lambda3;
    Curvature_Higgs_L4.at(3).at(5).at(3).at(6) = Imlambda7;
    Curvature_Higgs_L4.at(3).at(5).at(3).at(7) = Relambda7;
    Curvature_Higgs_L4.at(3).at(5).at(5).at(0) = -Imlambda6;
    Curvature_Higgs_L4.at(3).at(5).at(5).at(1) = Relambda6;
    Curvature_Higgs_L4.at(3).at(5).at(5).at(3) = lambda3;
    Curvature_Higgs_L4.at(3).at(5).at(6).at(0) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(3).at(5).at(6).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(3).at(5).at(6).at(3) = Imlambda7;
    Curvature_Higgs_L4.at(3).at(5).at(7).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(5).at(7).at(1) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(5).at(7).at(3) = Relambda7;
    Curvature_Higgs_L4.at(3).at(6).at(0).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(6).at(0).at(5) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(3).at(6).at(0).at(6) = -Imlambda7;
    Curvature_Higgs_L4.at(3).at(6).at(1).at(4) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(6).at(1).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(3).at(6).at(1).at(6) = Relambda7;
    Curvature_Higgs_L4.at(3).at(6).at(3).at(4) = Relambda7;
    Curvature_Higgs_L4.at(3).at(6).at(3).at(5) = Imlambda7;
    Curvature_Higgs_L4.at(3).at(6).at(3).at(6) = lambda2;
    Curvature_Higgs_L4.at(3).at(6).at(4).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(6).at(4).at(1) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(6).at(4).at(3) = Relambda7;
    Curvature_Higgs_L4.at(3).at(6).at(5).at(0) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(3).at(6).at(5).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(3).at(6).at(5).at(3) = Imlambda7;
    Curvature_Higgs_L4.at(3).at(6).at(6).at(0) = -Imlambda7;
    Curvature_Higgs_L4.at(3).at(6).at(6).at(1) = Relambda7;
    Curvature_Higgs_L4.at(3).at(6).at(6).at(3) = lambda2;
    Curvature_Higgs_L4.at(3).at(7).at(0).at(4) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(7).at(0).at(5) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(7).at(0).at(7) = -Imlambda7;
    Curvature_Higgs_L4.at(3).at(7).at(1).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(7).at(1).at(5) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(7).at(1).at(7) = Relambda7;
    Curvature_Higgs_L4.at(3).at(7).at(3).at(4) = -Imlambda7;
    Curvature_Higgs_L4.at(3).at(7).at(3).at(5) = Relambda7;
    Curvature_Higgs_L4.at(3).at(7).at(3).at(7) = lambda2;
    Curvature_Higgs_L4.at(3).at(7).at(4).at(0) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(7).at(4).at(1) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(7).at(4).at(3) = -Imlambda7;
    Curvature_Higgs_L4.at(3).at(7).at(5).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(3).at(7).at(5).at(1) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(3).at(7).at(5).at(3) = Relambda7;
    Curvature_Higgs_L4.at(3).at(7).at(7).at(0) = -Imlambda7;
    Curvature_Higgs_L4.at(3).at(7).at(7).at(1) = Relambda7;
    Curvature_Higgs_L4.at(3).at(7).at(7).at(3) = lambda2;
    Curvature_Higgs_L4.at(4).at(0).at(0).at(4) = lambda1;
    Curvature_Higgs_L4.at(4).at(0).at(0).at(6) = Relambda6;
    Curvature_Higgs_L4.at(4).at(0).at(0).at(7) = -Imlambda6;
    Curvature_Higgs_L4.at(4).at(0).at(2).at(4) = Relambda6;
    Curvature_Higgs_L4.at(4).at(0).at(2).at(6) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(0).at(2).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(0).at(3).at(4) = -Imlambda6;
    Curvature_Higgs_L4.at(4).at(0).at(3).at(6) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(0).at(3).at(7) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(0).at(4).at(0) = lambda1;
    Curvature_Higgs_L4.at(4).at(0).at(4).at(2) = Relambda6;
    Curvature_Higgs_L4.at(4).at(0).at(4).at(3) = -Imlambda6;
    Curvature_Higgs_L4.at(4).at(0).at(6).at(0) = Relambda6;
    Curvature_Higgs_L4.at(4).at(0).at(6).at(2) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(0).at(6).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(0).at(7).at(0) = -Imlambda6;
    Curvature_Higgs_L4.at(4).at(0).at(7).at(2) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(0).at(7).at(3) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(1).at(1).at(4) = lambda1;
    Curvature_Higgs_L4.at(4).at(1).at(1).at(6) = Relambda6;
    Curvature_Higgs_L4.at(4).at(1).at(1).at(7) = -Imlambda6;
    Curvature_Higgs_L4.at(4).at(1).at(2).at(4) = Imlambda6;
    Curvature_Higgs_L4.at(4).at(1).at(2).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(4).at(1).at(2).at(7) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(4).at(1).at(3).at(4) = Relambda6;
    Curvature_Higgs_L4.at(4).at(1).at(3).at(6) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(1).at(3).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(1).at(4).at(1) = lambda1;
    Curvature_Higgs_L4.at(4).at(1).at(4).at(2) = Imlambda6;
    Curvature_Higgs_L4.at(4).at(1).at(4).at(3) = Relambda6;
    Curvature_Higgs_L4.at(4).at(1).at(6).at(1) = Relambda6;
    Curvature_Higgs_L4.at(4).at(1).at(6).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(4).at(1).at(6).at(3) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(1).at(7).at(1) = -Imlambda6;
    Curvature_Higgs_L4.at(4).at(1).at(7).at(2) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(4).at(1).at(7).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(2).at(0).at(4) = Relambda6;
    Curvature_Higgs_L4.at(4).at(2).at(0).at(6) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(2).at(0).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(2).at(1).at(4) = Imlambda6;
    Curvature_Higgs_L4.at(4).at(2).at(1).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(4).at(2).at(1).at(7) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(4).at(2).at(2).at(4) = lambda3;
    Curvature_Higgs_L4.at(4).at(2).at(2).at(6) = Relambda7;
    Curvature_Higgs_L4.at(4).at(2).at(2).at(7) = -Imlambda7;
    Curvature_Higgs_L4.at(4).at(2).at(4).at(0) = Relambda6;
    Curvature_Higgs_L4.at(4).at(2).at(4).at(1) = Imlambda6;
    Curvature_Higgs_L4.at(4).at(2).at(4).at(2) = lambda3;
    Curvature_Higgs_L4.at(4).at(2).at(6).at(0) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(2).at(6).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(4).at(2).at(6).at(2) = Relambda7;
    Curvature_Higgs_L4.at(4).at(2).at(7).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(2).at(7).at(1) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(4).at(2).at(7).at(2) = -Imlambda7;
    Curvature_Higgs_L4.at(4).at(3).at(0).at(4) = -Imlambda6;
    Curvature_Higgs_L4.at(4).at(3).at(0).at(6) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(3).at(0).at(7) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(3).at(1).at(4) = Relambda6;
    Curvature_Higgs_L4.at(4).at(3).at(1).at(6) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(3).at(1).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(3).at(3).at(4) = lambda3;
    Curvature_Higgs_L4.at(4).at(3).at(3).at(6) = Relambda7;
    Curvature_Higgs_L4.at(4).at(3).at(3).at(7) = -Imlambda7;
    Curvature_Higgs_L4.at(4).at(3).at(4).at(0) = -Imlambda6;
    Curvature_Higgs_L4.at(4).at(3).at(4).at(1) = Relambda6;
    Curvature_Higgs_L4.at(4).at(3).at(4).at(3) = lambda3;
    Curvature_Higgs_L4.at(4).at(3).at(6).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(3).at(6).at(1) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(3).at(6).at(3) = Relambda7;
    Curvature_Higgs_L4.at(4).at(3).at(7).at(0) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(3).at(7).at(1) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(3).at(7).at(3) = -Imlambda7;
    Curvature_Higgs_L4.at(4).at(4).at(0).at(0) = lambda1;
    Curvature_Higgs_L4.at(4).at(4).at(0).at(2) = Relambda6;
    Curvature_Higgs_L4.at(4).at(4).at(0).at(3) = -Imlambda6;
    Curvature_Higgs_L4.at(4).at(4).at(1).at(1) = lambda1;
    Curvature_Higgs_L4.at(4).at(4).at(1).at(2) = Imlambda6;
    Curvature_Higgs_L4.at(4).at(4).at(1).at(3) = Relambda6;
    Curvature_Higgs_L4.at(4).at(4).at(2).at(0) = Relambda6;
    Curvature_Higgs_L4.at(4).at(4).at(2).at(1) = Imlambda6;
    Curvature_Higgs_L4.at(4).at(4).at(2).at(2) = lambda3;
    Curvature_Higgs_L4.at(4).at(4).at(3).at(0) = -Imlambda6;
    Curvature_Higgs_L4.at(4).at(4).at(3).at(1) = Relambda6;
    Curvature_Higgs_L4.at(4).at(4).at(3).at(3) = lambda3;
    Curvature_Higgs_L4.at(4).at(4).at(4).at(4) = 3 * lambda1;
    Curvature_Higgs_L4.at(4).at(4).at(4).at(6) = 3 * Relambda6;
    Curvature_Higgs_L4.at(4).at(4).at(4).at(7) = -3 * Imlambda6;
    Curvature_Higgs_L4.at(4).at(4).at(5).at(5) = lambda1;
    Curvature_Higgs_L4.at(4).at(4).at(5).at(6) = Imlambda6;
    Curvature_Higgs_L4.at(4).at(4).at(5).at(7) = Relambda6;
    Curvature_Higgs_L4.at(4).at(4).at(6).at(4) = 3 * Relambda6;
    Curvature_Higgs_L4.at(4).at(4).at(6).at(5) = Imlambda6;
    Curvature_Higgs_L4.at(4).at(4).at(6).at(6) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(4).at(4).at(6).at(7) = -Imlambda5;
    Curvature_Higgs_L4.at(4).at(4).at(7).at(4) = -3 * Imlambda6;
    Curvature_Higgs_L4.at(4).at(4).at(7).at(5) = Relambda6;
    Curvature_Higgs_L4.at(4).at(4).at(7).at(6) = -Imlambda5;
    Curvature_Higgs_L4.at(4).at(4).at(7).at(7) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(4).at(5).at(4).at(5) = lambda1;
    Curvature_Higgs_L4.at(4).at(5).at(4).at(6) = Imlambda6;
    Curvature_Higgs_L4.at(4).at(5).at(4).at(7) = Relambda6;
    Curvature_Higgs_L4.at(4).at(5).at(5).at(4) = lambda1;
    Curvature_Higgs_L4.at(4).at(5).at(5).at(6) = Relambda6;
    Curvature_Higgs_L4.at(4).at(5).at(5).at(7) = -Imlambda6;
    Curvature_Higgs_L4.at(4).at(5).at(6).at(4) = Imlambda6;
    Curvature_Higgs_L4.at(4).at(5).at(6).at(5) = Relambda6;
    Curvature_Higgs_L4.at(4).at(5).at(6).at(6) = Imlambda5;
    Curvature_Higgs_L4.at(4).at(5).at(6).at(7) = Relambda5;
    Curvature_Higgs_L4.at(4).at(5).at(7).at(4) = Relambda6;
    Curvature_Higgs_L4.at(4).at(5).at(7).at(5) = -Imlambda6;
    Curvature_Higgs_L4.at(4).at(5).at(7).at(6) = Relambda5;
    Curvature_Higgs_L4.at(4).at(5).at(7).at(7) = -Imlambda5;
    Curvature_Higgs_L4.at(4).at(6).at(0).at(0) = Relambda6;
    Curvature_Higgs_L4.at(4).at(6).at(0).at(2) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(6).at(0).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(6).at(1).at(1) = Relambda6;
    Curvature_Higgs_L4.at(4).at(6).at(1).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(4).at(6).at(1).at(3) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(6).at(2).at(0) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(6).at(2).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(4).at(6).at(2).at(2) = Relambda7;
    Curvature_Higgs_L4.at(4).at(6).at(3).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(6).at(3).at(1) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(6).at(3).at(3) = Relambda7;
    Curvature_Higgs_L4.at(4).at(6).at(4).at(4) = 3 * Relambda6;
    Curvature_Higgs_L4.at(4).at(6).at(4).at(5) = Imlambda6;
    Curvature_Higgs_L4.at(4).at(6).at(4).at(6) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(4).at(6).at(4).at(7) = -Imlambda5;
    Curvature_Higgs_L4.at(4).at(6).at(5).at(4) = Imlambda6;
    Curvature_Higgs_L4.at(4).at(6).at(5).at(5) = Relambda6;
    Curvature_Higgs_L4.at(4).at(6).at(5).at(6) = Imlambda5;
    Curvature_Higgs_L4.at(4).at(6).at(5).at(7) = Relambda5;
    Curvature_Higgs_L4.at(4).at(6).at(6).at(4) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(4).at(6).at(6).at(5) = Imlambda5;
    Curvature_Higgs_L4.at(4).at(6).at(6).at(6) = 3 * Relambda7;
    Curvature_Higgs_L4.at(4).at(6).at(6).at(7) = -Imlambda7;
    Curvature_Higgs_L4.at(4).at(6).at(7).at(4) = -Imlambda5;
    Curvature_Higgs_L4.at(4).at(6).at(7).at(5) = Relambda5;
    Curvature_Higgs_L4.at(4).at(6).at(7).at(6) = -Imlambda7;
    Curvature_Higgs_L4.at(4).at(6).at(7).at(7) = Relambda7;
    Curvature_Higgs_L4.at(4).at(7).at(0).at(0) = -Imlambda6;
    Curvature_Higgs_L4.at(4).at(7).at(0).at(2) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(7).at(0).at(3) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(7).at(1).at(1) = -Imlambda6;
    Curvature_Higgs_L4.at(4).at(7).at(1).at(2) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(4).at(7).at(1).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(7).at(2).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(7).at(2).at(1) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(4).at(7).at(2).at(2) = -Imlambda7;
    Curvature_Higgs_L4.at(4).at(7).at(3).at(0) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(4).at(7).at(3).at(1) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(4).at(7).at(3).at(3) = -Imlambda7;
    Curvature_Higgs_L4.at(4).at(7).at(4).at(4) = -3 * Imlambda6;
    Curvature_Higgs_L4.at(4).at(7).at(4).at(5) = Relambda6;
    Curvature_Higgs_L4.at(4).at(7).at(4).at(6) = -Imlambda5;
    Curvature_Higgs_L4.at(4).at(7).at(4).at(7) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(4).at(7).at(5).at(4) = Relambda6;
    Curvature_Higgs_L4.at(4).at(7).at(5).at(5) = -Imlambda6;
    Curvature_Higgs_L4.at(4).at(7).at(5).at(6) = Relambda5;
    Curvature_Higgs_L4.at(4).at(7).at(5).at(7) = -Imlambda5;
    Curvature_Higgs_L4.at(4).at(7).at(6).at(4) = -Imlambda5;
    Curvature_Higgs_L4.at(4).at(7).at(6).at(5) = Relambda5;
    Curvature_Higgs_L4.at(4).at(7).at(6).at(6) = -Imlambda7;
    Curvature_Higgs_L4.at(4).at(7).at(6).at(7) = Relambda7;
    Curvature_Higgs_L4.at(4).at(7).at(7).at(4) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(4).at(7).at(7).at(5) = -Imlambda5;
    Curvature_Higgs_L4.at(4).at(7).at(7).at(6) = Relambda7;
    Curvature_Higgs_L4.at(4).at(7).at(7).at(7) = -3 * Imlambda7;
    Curvature_Higgs_L4.at(5).at(0).at(0).at(5) = lambda1;
    Curvature_Higgs_L4.at(5).at(0).at(0).at(6) = Imlambda6;
    Curvature_Higgs_L4.at(5).at(0).at(0).at(7) = Relambda6;
    Curvature_Higgs_L4.at(5).at(0).at(2).at(5) = Relambda6;
    Curvature_Higgs_L4.at(5).at(0).at(2).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(0).at(2).at(7) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(0).at(3).at(5) = -Imlambda6;
    Curvature_Higgs_L4.at(5).at(0).at(3).at(6) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(5).at(0).at(3).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(5).at(0).at(5).at(0) = lambda1;
    Curvature_Higgs_L4.at(5).at(0).at(5).at(2) = Relambda6;
    Curvature_Higgs_L4.at(5).at(0).at(5).at(3) = -Imlambda6;
    Curvature_Higgs_L4.at(5).at(0).at(6).at(0) = Imlambda6;
    Curvature_Higgs_L4.at(5).at(0).at(6).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(0).at(6).at(3) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(5).at(0).at(7).at(0) = Relambda6;
    Curvature_Higgs_L4.at(5).at(0).at(7).at(2) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(0).at(7).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(5).at(1).at(1).at(5) = lambda1;
    Curvature_Higgs_L4.at(5).at(1).at(1).at(6) = Imlambda6;
    Curvature_Higgs_L4.at(5).at(1).at(1).at(7) = Relambda6;
    Curvature_Higgs_L4.at(5).at(1).at(2).at(5) = Imlambda6;
    Curvature_Higgs_L4.at(5).at(1).at(2).at(6) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(1).at(2).at(7) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(1).at(3).at(5) = Relambda6;
    Curvature_Higgs_L4.at(5).at(1).at(3).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(1).at(3).at(7) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(1).at(5).at(1) = lambda1;
    Curvature_Higgs_L4.at(5).at(1).at(5).at(2) = Imlambda6;
    Curvature_Higgs_L4.at(5).at(1).at(5).at(3) = Relambda6;
    Curvature_Higgs_L4.at(5).at(1).at(6).at(1) = Imlambda6;
    Curvature_Higgs_L4.at(5).at(1).at(6).at(2) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(1).at(6).at(3) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(1).at(7).at(1) = Relambda6;
    Curvature_Higgs_L4.at(5).at(1).at(7).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(1).at(7).at(3) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(2).at(0).at(5) = Relambda6;
    Curvature_Higgs_L4.at(5).at(2).at(0).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(2).at(0).at(7) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(2).at(1).at(5) = Imlambda6;
    Curvature_Higgs_L4.at(5).at(2).at(1).at(6) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(2).at(1).at(7) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(2).at(2).at(5) = lambda3;
    Curvature_Higgs_L4.at(5).at(2).at(2).at(6) = Imlambda7;
    Curvature_Higgs_L4.at(5).at(2).at(2).at(7) = Relambda7;
    Curvature_Higgs_L4.at(5).at(2).at(5).at(0) = Relambda6;
    Curvature_Higgs_L4.at(5).at(2).at(5).at(1) = Imlambda6;
    Curvature_Higgs_L4.at(5).at(2).at(5).at(2) = lambda3;
    Curvature_Higgs_L4.at(5).at(2).at(6).at(0) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(2).at(6).at(1) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(2).at(6).at(2) = Imlambda7;
    Curvature_Higgs_L4.at(5).at(2).at(7).at(0) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(2).at(7).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(2).at(7).at(2) = Relambda7;
    Curvature_Higgs_L4.at(5).at(3).at(0).at(5) = -Imlambda6;
    Curvature_Higgs_L4.at(5).at(3).at(0).at(6) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(5).at(3).at(0).at(7) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(5).at(3).at(1).at(5) = Relambda6;
    Curvature_Higgs_L4.at(5).at(3).at(1).at(6) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(3).at(1).at(7) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(3).at(3).at(5) = lambda3;
    Curvature_Higgs_L4.at(5).at(3).at(3).at(6) = Imlambda7;
    Curvature_Higgs_L4.at(5).at(3).at(3).at(7) = Relambda7;
    Curvature_Higgs_L4.at(5).at(3).at(5).at(0) = -Imlambda6;
    Curvature_Higgs_L4.at(5).at(3).at(5).at(1) = Relambda6;
    Curvature_Higgs_L4.at(5).at(3).at(5).at(3) = lambda3;
    Curvature_Higgs_L4.at(5).at(3).at(6).at(0) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(5).at(3).at(6).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(3).at(6).at(3) = Imlambda7;
    Curvature_Higgs_L4.at(5).at(3).at(7).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(5).at(3).at(7).at(1) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(3).at(7).at(3) = Relambda7;
    Curvature_Higgs_L4.at(5).at(4).at(4).at(5) = lambda1;
    Curvature_Higgs_L4.at(5).at(4).at(4).at(6) = Imlambda6;
    Curvature_Higgs_L4.at(5).at(4).at(4).at(7) = Relambda6;
    Curvature_Higgs_L4.at(5).at(4).at(5).at(4) = lambda1;
    Curvature_Higgs_L4.at(5).at(4).at(5).at(6) = Relambda6;
    Curvature_Higgs_L4.at(5).at(4).at(5).at(7) = -Imlambda6;
    Curvature_Higgs_L4.at(5).at(4).at(6).at(4) = Imlambda6;
    Curvature_Higgs_L4.at(5).at(4).at(6).at(5) = Relambda6;
    Curvature_Higgs_L4.at(5).at(4).at(6).at(6) = Imlambda5;
    Curvature_Higgs_L4.at(5).at(4).at(6).at(7) = Relambda5;
    Curvature_Higgs_L4.at(5).at(4).at(7).at(4) = Relambda6;
    Curvature_Higgs_L4.at(5).at(4).at(7).at(5) = -Imlambda6;
    Curvature_Higgs_L4.at(5).at(4).at(7).at(6) = Relambda5;
    Curvature_Higgs_L4.at(5).at(4).at(7).at(7) = -Imlambda5;
    Curvature_Higgs_L4.at(5).at(5).at(0).at(0) = lambda1;
    Curvature_Higgs_L4.at(5).at(5).at(0).at(2) = Relambda6;
    Curvature_Higgs_L4.at(5).at(5).at(0).at(3) = -Imlambda6;
    Curvature_Higgs_L4.at(5).at(5).at(1).at(1) = lambda1;
    Curvature_Higgs_L4.at(5).at(5).at(1).at(2) = Imlambda6;
    Curvature_Higgs_L4.at(5).at(5).at(1).at(3) = Relambda6;
    Curvature_Higgs_L4.at(5).at(5).at(2).at(0) = Relambda6;
    Curvature_Higgs_L4.at(5).at(5).at(2).at(1) = Imlambda6;
    Curvature_Higgs_L4.at(5).at(5).at(2).at(2) = lambda3;
    Curvature_Higgs_L4.at(5).at(5).at(3).at(0) = -Imlambda6;
    Curvature_Higgs_L4.at(5).at(5).at(3).at(1) = Relambda6;
    Curvature_Higgs_L4.at(5).at(5).at(3).at(3) = lambda3;
    Curvature_Higgs_L4.at(5).at(5).at(4).at(4) = lambda1;
    Curvature_Higgs_L4.at(5).at(5).at(4).at(6) = Relambda6;
    Curvature_Higgs_L4.at(5).at(5).at(4).at(7) = -Imlambda6;
    Curvature_Higgs_L4.at(5).at(5).at(5).at(5) = 3 * lambda1;
    Curvature_Higgs_L4.at(5).at(5).at(5).at(6) = 3 * Imlambda6;
    Curvature_Higgs_L4.at(5).at(5).at(5).at(7) = 3 * Relambda6;
    Curvature_Higgs_L4.at(5).at(5).at(6).at(4) = Relambda6;
    Curvature_Higgs_L4.at(5).at(5).at(6).at(5) = 3 * Imlambda6;
    Curvature_Higgs_L4.at(5).at(5).at(6).at(6) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(5).at(5).at(6).at(7) = Imlambda5;
    Curvature_Higgs_L4.at(5).at(5).at(7).at(4) = -Imlambda6;
    Curvature_Higgs_L4.at(5).at(5).at(7).at(5) = 3 * Relambda6;
    Curvature_Higgs_L4.at(5).at(5).at(7).at(6) = Imlambda5;
    Curvature_Higgs_L4.at(5).at(5).at(7).at(7) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(5).at(6).at(0).at(0) = Imlambda6;
    Curvature_Higgs_L4.at(5).at(6).at(0).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(6).at(0).at(3) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(5).at(6).at(1).at(1) = Imlambda6;
    Curvature_Higgs_L4.at(5).at(6).at(1).at(2) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(6).at(1).at(3) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(6).at(2).at(0) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(6).at(2).at(1) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(6).at(2).at(2) = Imlambda7;
    Curvature_Higgs_L4.at(5).at(6).at(3).at(0) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(5).at(6).at(3).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(6).at(3).at(3) = Imlambda7;
    Curvature_Higgs_L4.at(5).at(6).at(4).at(4) = Imlambda6;
    Curvature_Higgs_L4.at(5).at(6).at(4).at(5) = Relambda6;
    Curvature_Higgs_L4.at(5).at(6).at(4).at(6) = Imlambda5;
    Curvature_Higgs_L4.at(5).at(6).at(4).at(7) = Relambda5;
    Curvature_Higgs_L4.at(5).at(6).at(5).at(4) = Relambda6;
    Curvature_Higgs_L4.at(5).at(6).at(5).at(5) = 3 * Imlambda6;
    Curvature_Higgs_L4.at(5).at(6).at(5).at(6) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(5).at(6).at(5).at(7) = Imlambda5;
    Curvature_Higgs_L4.at(5).at(6).at(6).at(4) = Imlambda5;
    Curvature_Higgs_L4.at(5).at(6).at(6).at(5) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(5).at(6).at(6).at(6) = 3 * Imlambda7;
    Curvature_Higgs_L4.at(5).at(6).at(6).at(7) = Relambda7;
    Curvature_Higgs_L4.at(5).at(6).at(7).at(4) = Relambda5;
    Curvature_Higgs_L4.at(5).at(6).at(7).at(5) = Imlambda5;
    Curvature_Higgs_L4.at(5).at(6).at(7).at(6) = Relambda7;
    Curvature_Higgs_L4.at(5).at(6).at(7).at(7) = Imlambda7;
    Curvature_Higgs_L4.at(5).at(7).at(0).at(0) = Relambda6;
    Curvature_Higgs_L4.at(5).at(7).at(0).at(2) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(7).at(0).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(5).at(7).at(1).at(1) = Relambda6;
    Curvature_Higgs_L4.at(5).at(7).at(1).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(7).at(1).at(3) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(7).at(2).at(0) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(7).at(2).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(5).at(7).at(2).at(2) = Relambda7;
    Curvature_Higgs_L4.at(5).at(7).at(3).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(5).at(7).at(3).at(1) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(5).at(7).at(3).at(3) = Relambda7;
    Curvature_Higgs_L4.at(5).at(7).at(4).at(4) = Relambda6;
    Curvature_Higgs_L4.at(5).at(7).at(4).at(5) = -Imlambda6;
    Curvature_Higgs_L4.at(5).at(7).at(4).at(6) = Relambda5;
    Curvature_Higgs_L4.at(5).at(7).at(4).at(7) = -Imlambda5;
    Curvature_Higgs_L4.at(5).at(7).at(5).at(4) = -Imlambda6;
    Curvature_Higgs_L4.at(5).at(7).at(5).at(5) = 3 * Relambda6;
    Curvature_Higgs_L4.at(5).at(7).at(5).at(6) = Imlambda5;
    Curvature_Higgs_L4.at(5).at(7).at(5).at(7) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(5).at(7).at(6).at(4) = Relambda5;
    Curvature_Higgs_L4.at(5).at(7).at(6).at(5) = Imlambda5;
    Curvature_Higgs_L4.at(5).at(7).at(6).at(6) = Relambda7;
    Curvature_Higgs_L4.at(5).at(7).at(6).at(7) = Imlambda7;
    Curvature_Higgs_L4.at(5).at(7).at(7).at(4) = -Imlambda5;
    Curvature_Higgs_L4.at(5).at(7).at(7).at(5) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(5).at(7).at(7).at(6) = Imlambda7;
    Curvature_Higgs_L4.at(5).at(7).at(7).at(7) = 3 * Relambda7;
    Curvature_Higgs_L4.at(6).at(0).at(0).at(4) = Relambda6;
    Curvature_Higgs_L4.at(6).at(0).at(0).at(5) = Imlambda6;
    Curvature_Higgs_L4.at(6).at(0).at(0).at(6) = lambda3;
    Curvature_Higgs_L4.at(6).at(0).at(2).at(4) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(0).at(2).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(0).at(2).at(6) = Relambda7;
    Curvature_Higgs_L4.at(6).at(0).at(3).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(6).at(0).at(3).at(5) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(6).at(0).at(3).at(6) = -Imlambda7;
    Curvature_Higgs_L4.at(6).at(0).at(4).at(0) = Relambda6;
    Curvature_Higgs_L4.at(6).at(0).at(4).at(2) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(0).at(4).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(6).at(0).at(5).at(0) = Imlambda6;
    Curvature_Higgs_L4.at(6).at(0).at(5).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(0).at(5).at(3) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(6).at(0).at(6).at(0) = lambda3;
    Curvature_Higgs_L4.at(6).at(0).at(6).at(2) = Relambda7;
    Curvature_Higgs_L4.at(6).at(0).at(6).at(3) = -Imlambda7;
    Curvature_Higgs_L4.at(6).at(1).at(1).at(4) = Relambda6;
    Curvature_Higgs_L4.at(6).at(1).at(1).at(5) = Imlambda6;
    Curvature_Higgs_L4.at(6).at(1).at(1).at(6) = lambda3;
    Curvature_Higgs_L4.at(6).at(1).at(2).at(4) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(1).at(2).at(5) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(1).at(2).at(6) = Imlambda7;
    Curvature_Higgs_L4.at(6).at(1).at(3).at(4) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(1).at(3).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(1).at(3).at(6) = Relambda7;
    Curvature_Higgs_L4.at(6).at(1).at(4).at(1) = Relambda6;
    Curvature_Higgs_L4.at(6).at(1).at(4).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(1).at(4).at(3) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(1).at(5).at(1) = Imlambda6;
    Curvature_Higgs_L4.at(6).at(1).at(5).at(2) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(1).at(5).at(3) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(1).at(6).at(1) = lambda3;
    Curvature_Higgs_L4.at(6).at(1).at(6).at(2) = Imlambda7;
    Curvature_Higgs_L4.at(6).at(1).at(6).at(3) = Relambda7;
    Curvature_Higgs_L4.at(6).at(2).at(0).at(4) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(2).at(0).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(2).at(0).at(6) = Relambda7;
    Curvature_Higgs_L4.at(6).at(2).at(1).at(4) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(2).at(1).at(5) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(2).at(1).at(6) = Imlambda7;
    Curvature_Higgs_L4.at(6).at(2).at(2).at(4) = Relambda7;
    Curvature_Higgs_L4.at(6).at(2).at(2).at(5) = Imlambda7;
    Curvature_Higgs_L4.at(6).at(2).at(2).at(6) = lambda2;
    Curvature_Higgs_L4.at(6).at(2).at(4).at(0) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(2).at(4).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(2).at(4).at(2) = Relambda7;
    Curvature_Higgs_L4.at(6).at(2).at(5).at(0) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(2).at(5).at(1) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(2).at(5).at(2) = Imlambda7;
    Curvature_Higgs_L4.at(6).at(2).at(6).at(0) = Relambda7;
    Curvature_Higgs_L4.at(6).at(2).at(6).at(1) = Imlambda7;
    Curvature_Higgs_L4.at(6).at(2).at(6).at(2) = lambda2;
    Curvature_Higgs_L4.at(6).at(3).at(0).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(6).at(3).at(0).at(5) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(6).at(3).at(0).at(6) = -Imlambda7;
    Curvature_Higgs_L4.at(6).at(3).at(1).at(4) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(3).at(1).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(3).at(1).at(6) = Relambda7;
    Curvature_Higgs_L4.at(6).at(3).at(3).at(4) = Relambda7;
    Curvature_Higgs_L4.at(6).at(3).at(3).at(5) = Imlambda7;
    Curvature_Higgs_L4.at(6).at(3).at(3).at(6) = lambda2;
    Curvature_Higgs_L4.at(6).at(3).at(4).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(6).at(3).at(4).at(1) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(3).at(4).at(3) = Relambda7;
    Curvature_Higgs_L4.at(6).at(3).at(5).at(0) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(6).at(3).at(5).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(3).at(5).at(3) = Imlambda7;
    Curvature_Higgs_L4.at(6).at(3).at(6).at(0) = -Imlambda7;
    Curvature_Higgs_L4.at(6).at(3).at(6).at(1) = Relambda7;
    Curvature_Higgs_L4.at(6).at(3).at(6).at(3) = lambda2;
    Curvature_Higgs_L4.at(6).at(4).at(0).at(0) = Relambda6;
    Curvature_Higgs_L4.at(6).at(4).at(0).at(2) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(4).at(0).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(6).at(4).at(1).at(1) = Relambda6;
    Curvature_Higgs_L4.at(6).at(4).at(1).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(4).at(1).at(3) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(4).at(2).at(0) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(4).at(2).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(4).at(2).at(2) = Relambda7;
    Curvature_Higgs_L4.at(6).at(4).at(3).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(6).at(4).at(3).at(1) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(4).at(3).at(3) = Relambda7;
    Curvature_Higgs_L4.at(6).at(4).at(4).at(4) = 3 * Relambda6;
    Curvature_Higgs_L4.at(6).at(4).at(4).at(5) = Imlambda6;
    Curvature_Higgs_L4.at(6).at(4).at(4).at(6) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(6).at(4).at(4).at(7) = -Imlambda5;
    Curvature_Higgs_L4.at(6).at(4).at(5).at(4) = Imlambda6;
    Curvature_Higgs_L4.at(6).at(4).at(5).at(5) = Relambda6;
    Curvature_Higgs_L4.at(6).at(4).at(5).at(6) = Imlambda5;
    Curvature_Higgs_L4.at(6).at(4).at(5).at(7) = Relambda5;
    Curvature_Higgs_L4.at(6).at(4).at(6).at(4) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(6).at(4).at(6).at(5) = Imlambda5;
    Curvature_Higgs_L4.at(6).at(4).at(6).at(6) = 3 * Relambda7;
    Curvature_Higgs_L4.at(6).at(4).at(6).at(7) = -Imlambda7;
    Curvature_Higgs_L4.at(6).at(4).at(7).at(4) = -Imlambda5;
    Curvature_Higgs_L4.at(6).at(4).at(7).at(5) = Relambda5;
    Curvature_Higgs_L4.at(6).at(4).at(7).at(6) = -Imlambda7;
    Curvature_Higgs_L4.at(6).at(4).at(7).at(7) = Relambda7;
    Curvature_Higgs_L4.at(6).at(5).at(0).at(0) = Imlambda6;
    Curvature_Higgs_L4.at(6).at(5).at(0).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(5).at(0).at(3) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(6).at(5).at(1).at(1) = Imlambda6;
    Curvature_Higgs_L4.at(6).at(5).at(1).at(2) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(5).at(1).at(3) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(5).at(2).at(0) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(5).at(2).at(1) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(6).at(5).at(2).at(2) = Imlambda7;
    Curvature_Higgs_L4.at(6).at(5).at(3).at(0) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(6).at(5).at(3).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(6).at(5).at(3).at(3) = Imlambda7;
    Curvature_Higgs_L4.at(6).at(5).at(4).at(4) = Imlambda6;
    Curvature_Higgs_L4.at(6).at(5).at(4).at(5) = Relambda6;
    Curvature_Higgs_L4.at(6).at(5).at(4).at(6) = Imlambda5;
    Curvature_Higgs_L4.at(6).at(5).at(4).at(7) = Relambda5;
    Curvature_Higgs_L4.at(6).at(5).at(5).at(4) = Relambda6;
    Curvature_Higgs_L4.at(6).at(5).at(5).at(5) = 3 * Imlambda6;
    Curvature_Higgs_L4.at(6).at(5).at(5).at(6) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(6).at(5).at(5).at(7) = Imlambda5;
    Curvature_Higgs_L4.at(6).at(5).at(6).at(4) = Imlambda5;
    Curvature_Higgs_L4.at(6).at(5).at(6).at(5) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(6).at(5).at(6).at(6) = 3 * Imlambda7;
    Curvature_Higgs_L4.at(6).at(5).at(6).at(7) = Relambda7;
    Curvature_Higgs_L4.at(6).at(5).at(7).at(4) = Relambda5;
    Curvature_Higgs_L4.at(6).at(5).at(7).at(5) = Imlambda5;
    Curvature_Higgs_L4.at(6).at(5).at(7).at(6) = Relambda7;
    Curvature_Higgs_L4.at(6).at(5).at(7).at(7) = Imlambda7;
    Curvature_Higgs_L4.at(6).at(6).at(0).at(0) = lambda3;
    Curvature_Higgs_L4.at(6).at(6).at(0).at(2) = Relambda7;
    Curvature_Higgs_L4.at(6).at(6).at(0).at(3) = -Imlambda7;
    Curvature_Higgs_L4.at(6).at(6).at(1).at(1) = lambda3;
    Curvature_Higgs_L4.at(6).at(6).at(1).at(2) = Imlambda7;
    Curvature_Higgs_L4.at(6).at(6).at(1).at(3) = Relambda7;
    Curvature_Higgs_L4.at(6).at(6).at(2).at(0) = Relambda7;
    Curvature_Higgs_L4.at(6).at(6).at(2).at(1) = Imlambda7;
    Curvature_Higgs_L4.at(6).at(6).at(2).at(2) = lambda2;
    Curvature_Higgs_L4.at(6).at(6).at(3).at(0) = -Imlambda7;
    Curvature_Higgs_L4.at(6).at(6).at(3).at(1) = Relambda7;
    Curvature_Higgs_L4.at(6).at(6).at(3).at(3) = lambda2;
    Curvature_Higgs_L4.at(6).at(6).at(4).at(4) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(6).at(6).at(4).at(5) = Imlambda5;
    Curvature_Higgs_L4.at(6).at(6).at(4).at(6) = 3 * Relambda7;
    Curvature_Higgs_L4.at(6).at(6).at(4).at(7) = -Imlambda7;
    Curvature_Higgs_L4.at(6).at(6).at(5).at(4) = Imlambda5;
    Curvature_Higgs_L4.at(6).at(6).at(5).at(5) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(6).at(6).at(5).at(6) = 3 * Imlambda7;
    Curvature_Higgs_L4.at(6).at(6).at(5).at(7) = Relambda7;
    Curvature_Higgs_L4.at(6).at(6).at(6).at(4) = 3 * Relambda7;
    Curvature_Higgs_L4.at(6).at(6).at(6).at(5) = 3 * Imlambda7;
    Curvature_Higgs_L4.at(6).at(6).at(6).at(6) = 3 * lambda2;
    Curvature_Higgs_L4.at(6).at(6).at(7).at(4) = -Imlambda7;
    Curvature_Higgs_L4.at(6).at(6).at(7).at(5) = Relambda7;
    Curvature_Higgs_L4.at(6).at(6).at(7).at(7) = lambda2;
    Curvature_Higgs_L4.at(6).at(7).at(4).at(4) = -Imlambda5;
    Curvature_Higgs_L4.at(6).at(7).at(4).at(5) = Relambda5;
    Curvature_Higgs_L4.at(6).at(7).at(4).at(6) = -Imlambda7;
    Curvature_Higgs_L4.at(6).at(7).at(4).at(7) = Relambda7;
    Curvature_Higgs_L4.at(6).at(7).at(5).at(4) = Relambda5;
    Curvature_Higgs_L4.at(6).at(7).at(5).at(5) = Imlambda5;
    Curvature_Higgs_L4.at(6).at(7).at(5).at(6) = Relambda7;
    Curvature_Higgs_L4.at(6).at(7).at(5).at(7) = Imlambda7;
    Curvature_Higgs_L4.at(6).at(7).at(6).at(4) = -Imlambda7;
    Curvature_Higgs_L4.at(6).at(7).at(6).at(5) = Relambda7;
    Curvature_Higgs_L4.at(6).at(7).at(6).at(7) = lambda2;
    Curvature_Higgs_L4.at(6).at(7).at(7).at(4) = Relambda7;
    Curvature_Higgs_L4.at(6).at(7).at(7).at(5) = Imlambda7;
    Curvature_Higgs_L4.at(6).at(7).at(7).at(6) = lambda2;
    Curvature_Higgs_L4.at(7).at(0).at(0).at(4) = -Imlambda6;
    Curvature_Higgs_L4.at(7).at(0).at(0).at(5) = Relambda6;
    Curvature_Higgs_L4.at(7).at(0).at(0).at(7) = lambda3;
    Curvature_Higgs_L4.at(7).at(0).at(2).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(0).at(2).at(5) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(0).at(2).at(7) = Relambda7;
    Curvature_Higgs_L4.at(7).at(0).at(3).at(4) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(0).at(3).at(5) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(0).at(3).at(7) = -Imlambda7;
    Curvature_Higgs_L4.at(7).at(0).at(4).at(0) = -Imlambda6;
    Curvature_Higgs_L4.at(7).at(0).at(4).at(2) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(0).at(4).at(3) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(0).at(5).at(0) = Relambda6;
    Curvature_Higgs_L4.at(7).at(0).at(5).at(2) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(0).at(5).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(0).at(7).at(0) = lambda3;
    Curvature_Higgs_L4.at(7).at(0).at(7).at(2) = Relambda7;
    Curvature_Higgs_L4.at(7).at(0).at(7).at(3) = -Imlambda7;
    Curvature_Higgs_L4.at(7).at(1).at(1).at(4) = -Imlambda6;
    Curvature_Higgs_L4.at(7).at(1).at(1).at(5) = Relambda6;
    Curvature_Higgs_L4.at(7).at(1).at(1).at(7) = lambda3;
    Curvature_Higgs_L4.at(7).at(1).at(2).at(4) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(7).at(1).at(2).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(7).at(1).at(2).at(7) = Imlambda7;
    Curvature_Higgs_L4.at(7).at(1).at(3).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(1).at(3).at(5) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(1).at(3).at(7) = Relambda7;
    Curvature_Higgs_L4.at(7).at(1).at(4).at(1) = -Imlambda6;
    Curvature_Higgs_L4.at(7).at(1).at(4).at(2) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(7).at(1).at(4).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(1).at(5).at(1) = Relambda6;
    Curvature_Higgs_L4.at(7).at(1).at(5).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(7).at(1).at(5).at(3) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(1).at(7).at(1) = lambda3;
    Curvature_Higgs_L4.at(7).at(1).at(7).at(2) = Imlambda7;
    Curvature_Higgs_L4.at(7).at(1).at(7).at(3) = Relambda7;
    Curvature_Higgs_L4.at(7).at(2).at(0).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(2).at(0).at(5) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(2).at(0).at(7) = Relambda7;
    Curvature_Higgs_L4.at(7).at(2).at(1).at(4) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(7).at(2).at(1).at(5) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(7).at(2).at(1).at(7) = Imlambda7;
    Curvature_Higgs_L4.at(7).at(2).at(2).at(4) = -Imlambda7;
    Curvature_Higgs_L4.at(7).at(2).at(2).at(5) = Relambda7;
    Curvature_Higgs_L4.at(7).at(2).at(2).at(7) = lambda2;
    Curvature_Higgs_L4.at(7).at(2).at(4).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(2).at(4).at(1) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(7).at(2).at(4).at(2) = -Imlambda7;
    Curvature_Higgs_L4.at(7).at(2).at(5).at(0) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(2).at(5).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(7).at(2).at(5).at(2) = Relambda7;
    Curvature_Higgs_L4.at(7).at(2).at(7).at(0) = Relambda7;
    Curvature_Higgs_L4.at(7).at(2).at(7).at(1) = Imlambda7;
    Curvature_Higgs_L4.at(7).at(2).at(7).at(2) = lambda2;
    Curvature_Higgs_L4.at(7).at(3).at(0).at(4) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(3).at(0).at(5) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(3).at(0).at(7) = -Imlambda7;
    Curvature_Higgs_L4.at(7).at(3).at(1).at(4) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(3).at(1).at(5) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(3).at(1).at(7) = Relambda7;
    Curvature_Higgs_L4.at(7).at(3).at(3).at(4) = -Imlambda7;
    Curvature_Higgs_L4.at(7).at(3).at(3).at(5) = Relambda7;
    Curvature_Higgs_L4.at(7).at(3).at(3).at(7) = lambda2;
    Curvature_Higgs_L4.at(7).at(3).at(4).at(0) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(3).at(4).at(1) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(3).at(4).at(3) = -Imlambda7;
    Curvature_Higgs_L4.at(7).at(3).at(5).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(3).at(5).at(1) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(3).at(5).at(3) = Relambda7;
    Curvature_Higgs_L4.at(7).at(3).at(7).at(0) = -Imlambda7;
    Curvature_Higgs_L4.at(7).at(3).at(7).at(1) = Relambda7;
    Curvature_Higgs_L4.at(7).at(3).at(7).at(3) = lambda2;
    Curvature_Higgs_L4.at(7).at(4).at(0).at(0) = -Imlambda6;
    Curvature_Higgs_L4.at(7).at(4).at(0).at(2) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(4).at(0).at(3) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(4).at(1).at(1) = -Imlambda6;
    Curvature_Higgs_L4.at(7).at(4).at(1).at(2) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(7).at(4).at(1).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(4).at(2).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(4).at(2).at(1) =
        (1.0 / 2.0) * Relambda5 - 1.0 / 2.0 * lambda4;
    Curvature_Higgs_L4.at(7).at(4).at(2).at(2) = -Imlambda7;
    Curvature_Higgs_L4.at(7).at(4).at(3).at(0) =
        -1.0 / 2.0 * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(4).at(3).at(1) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(4).at(3).at(3) = -Imlambda7;
    Curvature_Higgs_L4.at(7).at(4).at(4).at(4) = -3 * Imlambda6;
    Curvature_Higgs_L4.at(7).at(4).at(4).at(5) = Relambda6;
    Curvature_Higgs_L4.at(7).at(4).at(4).at(6) = -Imlambda5;
    Curvature_Higgs_L4.at(7).at(4).at(4).at(7) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(7).at(4).at(5).at(4) = Relambda6;
    Curvature_Higgs_L4.at(7).at(4).at(5).at(5) = -Imlambda6;
    Curvature_Higgs_L4.at(7).at(4).at(5).at(6) = Relambda5;
    Curvature_Higgs_L4.at(7).at(4).at(5).at(7) = -Imlambda5;
    Curvature_Higgs_L4.at(7).at(4).at(6).at(4) = -Imlambda5;
    Curvature_Higgs_L4.at(7).at(4).at(6).at(5) = Relambda5;
    Curvature_Higgs_L4.at(7).at(4).at(6).at(6) = -Imlambda7;
    Curvature_Higgs_L4.at(7).at(4).at(6).at(7) = Relambda7;
    Curvature_Higgs_L4.at(7).at(4).at(7).at(4) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(7).at(4).at(7).at(5) = -Imlambda5;
    Curvature_Higgs_L4.at(7).at(4).at(7).at(6) = Relambda7;
    Curvature_Higgs_L4.at(7).at(4).at(7).at(7) = -3 * Imlambda7;
    Curvature_Higgs_L4.at(7).at(5).at(0).at(0) = Relambda6;
    Curvature_Higgs_L4.at(7).at(5).at(0).at(2) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(5).at(0).at(3) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(5).at(1).at(1) = Relambda6;
    Curvature_Higgs_L4.at(7).at(5).at(1).at(2) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(7).at(5).at(1).at(3) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(5).at(2).at(0) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(5).at(2).at(1) = (1.0 / 2.0) * Imlambda5;
    Curvature_Higgs_L4.at(7).at(5).at(2).at(2) = Relambda7;
    Curvature_Higgs_L4.at(7).at(5).at(3).at(0) = -1.0 / 2.0 * Imlambda5;
    Curvature_Higgs_L4.at(7).at(5).at(3).at(1) =
        (1.0 / 2.0) * Relambda5 + (1.0 / 2.0) * lambda4;
    Curvature_Higgs_L4.at(7).at(5).at(3).at(3) = Relambda7;
    Curvature_Higgs_L4.at(7).at(5).at(4).at(4) = Relambda6;
    Curvature_Higgs_L4.at(7).at(5).at(4).at(5) = -Imlambda6;
    Curvature_Higgs_L4.at(7).at(5).at(4).at(6) = Relambda5;
    Curvature_Higgs_L4.at(7).at(5).at(4).at(7) = -Imlambda5;
    Curvature_Higgs_L4.at(7).at(5).at(5).at(4) = -Imlambda6;
    Curvature_Higgs_L4.at(7).at(5).at(5).at(5) = 3 * Relambda6;
    Curvature_Higgs_L4.at(7).at(5).at(5).at(6) = Imlambda5;
    Curvature_Higgs_L4.at(7).at(5).at(5).at(7) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(7).at(5).at(6).at(4) = Relambda5;
    Curvature_Higgs_L4.at(7).at(5).at(6).at(5) = Imlambda5;
    Curvature_Higgs_L4.at(7).at(5).at(6).at(6) = Relambda7;
    Curvature_Higgs_L4.at(7).at(5).at(6).at(7) = Imlambda7;
    Curvature_Higgs_L4.at(7).at(5).at(7).at(4) = -Imlambda5;
    Curvature_Higgs_L4.at(7).at(5).at(7).at(5) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(7).at(5).at(7).at(6) = Imlambda7;
    Curvature_Higgs_L4.at(7).at(5).at(7).at(7) = 3 * Relambda7;
    Curvature_Higgs_L4.at(7).at(6).at(4).at(4) = -Imlambda5;
    Curvature_Higgs_L4.at(7).at(6).at(4).at(5) = Relambda5;
    Curvature_Higgs_L4.at(7).at(6).at(4).at(6) = -Imlambda7;
    Curvature_Higgs_L4.at(7).at(6).at(4).at(7) = Relambda7;
    Curvature_Higgs_L4.at(7).at(6).at(5).at(4) = Relambda5;
    Curvature_Higgs_L4.at(7).at(6).at(5).at(5) = Imlambda5;
    Curvature_Higgs_L4.at(7).at(6).at(5).at(6) = Relambda7;
    Curvature_Higgs_L4.at(7).at(6).at(5).at(7) = Imlambda7;
    Curvature_Higgs_L4.at(7).at(6).at(6).at(4) = -Imlambda7;
    Curvature_Higgs_L4.at(7).at(6).at(6).at(5) = Relambda7;
    Curvature_Higgs_L4.at(7).at(6).at(6).at(7) = lambda2;
    Curvature_Higgs_L4.at(7).at(6).at(7).at(4) = Relambda7;
    Curvature_Higgs_L4.at(7).at(6).at(7).at(5) = Imlambda7;
    Curvature_Higgs_L4.at(7).at(6).at(7).at(6) = lambda2;
    Curvature_Higgs_L4.at(7).at(7).at(0).at(0) = lambda3;
    Curvature_Higgs_L4.at(7).at(7).at(0).at(2) = Relambda7;
    Curvature_Higgs_L4.at(7).at(7).at(0).at(3) = -Imlambda7;
    Curvature_Higgs_L4.at(7).at(7).at(1).at(1) = lambda3;
    Curvature_Higgs_L4.at(7).at(7).at(1).at(2) = Imlambda7;
    Curvature_Higgs_L4.at(7).at(7).at(1).at(3) = Relambda7;
    Curvature_Higgs_L4.at(7).at(7).at(2).at(0) = Relambda7;
    Curvature_Higgs_L4.at(7).at(7).at(2).at(1) = Imlambda7;
    Curvature_Higgs_L4.at(7).at(7).at(2).at(2) = lambda2;
    Curvature_Higgs_L4.at(7).at(7).at(3).at(0) = -Imlambda7;
    Curvature_Higgs_L4.at(7).at(7).at(3).at(1) = Relambda7;
    Curvature_Higgs_L4.at(7).at(7).at(3).at(3) = lambda2;
    Curvature_Higgs_L4.at(7).at(7).at(4).at(4) = -Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(7).at(7).at(4).at(5) = -Imlambda5;
    Curvature_Higgs_L4.at(7).at(7).at(4).at(6) = Relambda7;
    Curvature_Higgs_L4.at(7).at(7).at(4).at(7) = -3 * Imlambda7;
    Curvature_Higgs_L4.at(7).at(7).at(5).at(4) = -Imlambda5;
    Curvature_Higgs_L4.at(7).at(7).at(5).at(5) = Relambda5 + lambda3 + lambda4;
    Curvature_Higgs_L4.at(7).at(7).at(5).at(6) = Imlambda7;
    Curvature_Higgs_L4.at(7).at(7).at(5).at(7) = 3 * Relambda7;
    Curvature_Higgs_L4.at(7).at(7).at(6).at(4) = Relambda7;
    Curvature_Higgs_L4.at(7).at(7).at(6).at(5) = Imlambda7;
    Curvature_Higgs_L4.at(7).at(7).at(6).at(6) = lambda2;
    Curvature_Higgs_L4.at(7).at(7).at(7).at(4) = -3 * Imlambda7;
    Curvature_Higgs_L4.at(7).at(7).at(7).at(5) = 3 * Relambda7;
    Curvature_Higgs_L4.at(7).at(7).at(7).at(7) = 3 * lambda2;
  }
  // End of Higgs curvature tensors

  // Begin of Gauge interaction tensors
  {
    Curvature_Gauge_G2H2.at(0).at(0).at(0).at(0) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(0).at(0).at(1).at(1) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(0).at(0).at(2).at(2) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(0).at(0).at(3).at(3) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(0).at(0).at(4).at(4) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(0).at(0).at(5).at(5) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(0).at(0).at(6).at(6) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(0).at(0).at(7).at(7) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(0).at(3).at(0).at(4) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(0).at(3).at(1).at(5) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(0).at(3).at(2).at(6) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(0).at(3).at(3).at(7) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(0).at(3).at(4).at(0) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(0).at(3).at(5).at(1) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(0).at(3).at(6).at(2) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(0).at(3).at(7).at(3) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(1).at(1).at(0).at(0) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(1).at(1).at(1).at(1) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(1).at(1).at(2).at(2) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(1).at(1).at(3).at(3) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(1).at(1).at(4).at(4) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(1).at(1).at(5).at(5) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(1).at(1).at(6).at(6) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(1).at(1).at(7).at(7) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(1).at(3).at(0).at(5) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(1).at(3).at(1).at(4) = -1.0 / 2.0 * C_g * C_gs;
    Curvature_Gauge_G2H2.at(1).at(3).at(2).at(7) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(1).at(3).at(3).at(6) = -1.0 / 2.0 * C_g * C_gs;
    Curvature_Gauge_G2H2.at(1).at(3).at(4).at(1) = -1.0 / 2.0 * C_g * C_gs;
    Curvature_Gauge_G2H2.at(1).at(3).at(5).at(0) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(1).at(3).at(6).at(3) = -1.0 / 2.0 * C_g * C_gs;
    Curvature_Gauge_G2H2.at(1).at(3).at(7).at(2) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(2).at(2).at(0).at(0) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(2).at(2).at(1).at(1) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(2).at(2).at(2).at(2) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(2).at(2).at(3).at(3) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(2).at(2).at(4).at(4) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(2).at(2).at(5).at(5) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(2).at(2).at(6).at(6) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(2).at(2).at(7).at(7) =
        (1.0 / 2.0) * std::pow(C_g, 2);
    Curvature_Gauge_G2H2.at(2).at(3).at(0).at(0) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(2).at(3).at(1).at(1) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(2).at(3).at(2).at(2) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(2).at(3).at(3).at(3) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(2).at(3).at(4).at(4) = -1.0 / 2.0 * C_g * C_gs;
    Curvature_Gauge_G2H2.at(2).at(3).at(5).at(5) = -1.0 / 2.0 * C_g * C_gs;
    Curvature_Gauge_G2H2.at(2).at(3).at(6).at(6) = -1.0 / 2.0 * C_g * C_gs;
    Curvature_Gauge_G2H2.at(2).at(3).at(7).at(7) = -1.0 / 2.0 * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(0).at(0).at(4) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(0).at(1).at(5) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(0).at(2).at(6) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(0).at(3).at(7) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(0).at(4).at(0) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(0).at(5).at(1) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(0).at(6).at(2) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(0).at(7).at(3) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(1).at(0).at(5) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(1).at(1).at(4) = -1.0 / 2.0 * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(1).at(2).at(7) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(1).at(3).at(6) = -1.0 / 2.0 * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(1).at(4).at(1) = -1.0 / 2.0 * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(1).at(5).at(0) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(1).at(6).at(3) = -1.0 / 2.0 * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(1).at(7).at(2) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(2).at(0).at(0) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(2).at(1).at(1) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(2).at(2).at(2) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(2).at(3).at(3) = (1.0 / 2.0) * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(2).at(4).at(4) = -1.0 / 2.0 * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(2).at(5).at(5) = -1.0 / 2.0 * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(2).at(6).at(6) = -1.0 / 2.0 * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(2).at(7).at(7) = -1.0 / 2.0 * C_g * C_gs;
    Curvature_Gauge_G2H2.at(3).at(3).at(0).at(0) =
        (1.0 / 2.0) * std::pow(C_gs, 2);
    Curvature_Gauge_G2H2.at(3).at(3).at(1).at(1) =
        (1.0 / 2.0) * std::pow(C_gs, 2);
    Curvature_Gauge_G2H2.at(3).at(3).at(2).at(2) =
        (1.0 / 2.0) * std::pow(C_gs, 2);
    Curvature_Gauge_G2H2.at(3).at(3).at(3).at(3) =
        (1.0 / 2.0) * std::pow(C_gs, 2);
    Curvature_Gauge_G2H2.at(3).at(3).at(4).at(4) =
        (1.0 / 2.0) * std::pow(C_gs, 2);
    Curvature_Gauge_G2H2.at(3).at(3).at(5).at(5) =
        (1.0 / 2.0) * std::pow(C_gs, 2);
    Curvature_Gauge_G2H2.at(3).at(3).at(6).at(6) =
        (1.0 / 2.0) * std::pow(C_gs, 2);
    Curvature_Gauge_G2H2.at(3).at(3).at(7).at(7) =
        (1.0 / 2.0) * std::pow(C_gs, 2);
  }
  // End of Gauge interaction tensors

  double v1 = C_vev0 * C_CosBeta;
  double v2 = C_vev0 * C_SinBeta;

  auto Vub = C_Vub;
  auto Vus = C_Vus;
  auto Vud = C_Vud;
  auto Vcb = C_Vcb;
  auto Vcs = C_Vcs;
  auto Vcd = C_Vcd;
  auto Vtb = C_Vtb;
  auto Vts = C_Vts;
  auto Vtd = C_Vtd;

  auto LeptonTypeI = [&]() {
    // Begin of Lepton interaction tensors
    Curvature_Lepton_F2H1.at(0).at(3).at(2) = C_MassElectron / v2;
    Curvature_Lepton_F2H1.at(0).at(3).at(3) = C_MassElectron * II / v2;
    Curvature_Lepton_F2H1.at(1).at(4).at(2) = C_MassMu / v2;
    Curvature_Lepton_F2H1.at(1).at(4).at(3) = C_MassMu * II / v2;
    Curvature_Lepton_F2H1.at(2).at(5).at(2) = C_MassTau / v2;
    Curvature_Lepton_F2H1.at(2).at(5).at(3) = C_MassTau * II / v2;
    Curvature_Lepton_F2H1.at(3).at(0).at(2) = C_MassElectron / v2;
    Curvature_Lepton_F2H1.at(3).at(0).at(3) = C_MassElectron * II / v2;
    Curvature_Lepton_F2H1.at(3).at(6).at(6) = C_MassElectron / v2;
    Curvature_Lepton_F2H1.at(3).at(6).at(7) = C_MassElectron * II / v2;
    Curvature_Lepton_F2H1.at(4).at(1).at(2) = C_MassMu / v2;
    Curvature_Lepton_F2H1.at(4).at(1).at(3) = C_MassMu * II / v2;
    Curvature_Lepton_F2H1.at(4).at(7).at(6) = C_MassMu / v2;
    Curvature_Lepton_F2H1.at(4).at(7).at(7) = C_MassMu * II / v2;
    Curvature_Lepton_F2H1.at(5).at(2).at(2) = C_MassTau / v2;
    Curvature_Lepton_F2H1.at(5).at(2).at(3) = C_MassTau * II / v2;
    Curvature_Lepton_F2H1.at(5).at(8).at(6) = C_MassTau / v2;
    Curvature_Lepton_F2H1.at(5).at(8).at(7) = C_MassTau * II / v2;
    Curvature_Lepton_F2H1.at(6).at(3).at(6) = C_MassElectron / v2;
    Curvature_Lepton_F2H1.at(6).at(3).at(7) = C_MassElectron * II / v2;
    Curvature_Lepton_F2H1.at(7).at(4).at(6) = C_MassMu / v2;
    Curvature_Lepton_F2H1.at(7).at(4).at(7) = C_MassMu * II / v2;
    Curvature_Lepton_F2H1.at(8).at(5).at(6) = C_MassTau / v2;
    Curvature_Lepton_F2H1.at(8).at(5).at(7) = C_MassTau * II / v2;
    // End of Lepton interaction tensors
  };

  auto QuarkTypeI = [&]() {
    // Begin of Quark interaction tensors
    Curvature_Quark_F2H1.at(0).at(6).at(6)  = C_MassUp / v2;
    Curvature_Quark_F2H1.at(0).at(6).at(7)  = -C_MassUp * II / v2;
    Curvature_Quark_F2H1.at(0).at(9).at(2)  = -C_MassUp * conj(Vud) / v2;
    Curvature_Quark_F2H1.at(0).at(9).at(3)  = C_MassUp * II * conj(Vud) / v2;
    Curvature_Quark_F2H1.at(0).at(10).at(2) = -C_MassUp * conj(Vus) / v2;
    Curvature_Quark_F2H1.at(0).at(10).at(3) = C_MassUp * II * conj(Vus) / v2;
    Curvature_Quark_F2H1.at(0).at(11).at(2) = -C_MassUp * conj(Vub) / v2;
    Curvature_Quark_F2H1.at(0).at(11).at(3) = C_MassUp * II * conj(Vub) / v2;
    Curvature_Quark_F2H1.at(1).at(7).at(6)  = C_MassCharm / v2;
    Curvature_Quark_F2H1.at(1).at(7).at(7)  = -C_MassCharm * II / v2;
    Curvature_Quark_F2H1.at(1).at(9).at(2)  = -C_MassCharm * conj(Vcd) / v2;
    Curvature_Quark_F2H1.at(1).at(9).at(3)  = C_MassCharm * II * conj(Vcd) / v2;
    Curvature_Quark_F2H1.at(1).at(10).at(2) = -C_MassCharm * conj(Vcs) / v2;
    Curvature_Quark_F2H1.at(1).at(10).at(3) = C_MassCharm * II * conj(Vcs) / v2;
    Curvature_Quark_F2H1.at(1).at(11).at(2) = -C_MassCharm * conj(Vcb) / v2;
    Curvature_Quark_F2H1.at(1).at(11).at(3) = C_MassCharm * II * conj(Vcb) / v2;
    Curvature_Quark_F2H1.at(2).at(8).at(6)  = C_MassTop / v2;
    Curvature_Quark_F2H1.at(2).at(8).at(7)  = -C_MassTop * II / v2;
    Curvature_Quark_F2H1.at(2).at(9).at(2)  = -C_MassTop * conj(Vtd) / v2;
    Curvature_Quark_F2H1.at(2).at(9).at(3)  = C_MassTop * II * conj(Vtd) / v2;
    Curvature_Quark_F2H1.at(2).at(10).at(2) = -C_MassTop * conj(Vts) / v2;
    Curvature_Quark_F2H1.at(2).at(10).at(3) = C_MassTop * II * conj(Vts) / v2;
    Curvature_Quark_F2H1.at(2).at(11).at(2) = -C_MassTop * conj(Vtb) / v2;
    Curvature_Quark_F2H1.at(2).at(11).at(3) = C_MassTop * II * conj(Vtb) / v2;
    Curvature_Quark_F2H1.at(3).at(6).at(2)  = C_MassDown * Vud / v2;
    Curvature_Quark_F2H1.at(3).at(6).at(3)  = C_MassDown * II * Vud / v2;
    Curvature_Quark_F2H1.at(3).at(7).at(2)  = C_MassDown * Vcd / v2;
    Curvature_Quark_F2H1.at(3).at(7).at(3)  = C_MassDown * II * Vcd / v2;
    Curvature_Quark_F2H1.at(3).at(8).at(2)  = C_MassDown * Vtd / v2;
    Curvature_Quark_F2H1.at(3).at(8).at(3)  = C_MassDown * II * Vtd / v2;
    Curvature_Quark_F2H1.at(3).at(9).at(6)  = C_MassDown / v2;
    Curvature_Quark_F2H1.at(3).at(9).at(7)  = C_MassDown * II / v2;
    Curvature_Quark_F2H1.at(4).at(6).at(2)  = C_MassStrange * Vus / v2;
    Curvature_Quark_F2H1.at(4).at(6).at(3)  = C_MassStrange * II * Vus / v2;
    Curvature_Quark_F2H1.at(4).at(7).at(2)  = C_MassStrange * Vcs / v2;
    Curvature_Quark_F2H1.at(4).at(7).at(3)  = C_MassStrange * II * Vcs / v2;
    Curvature_Quark_F2H1.at(4).at(8).at(2)  = C_MassStrange * Vts / v2;
    Curvature_Quark_F2H1.at(4).at(8).at(3)  = C_MassStrange * II * Vts / v2;
    Curvature_Quark_F2H1.at(4).at(10).at(6) = C_MassStrange / v2;
    Curvature_Quark_F2H1.at(4).at(10).at(7) = C_MassStrange * II / v2;
    Curvature_Quark_F2H1.at(5).at(6).at(2)  = C_MassBottom * Vub / v2;
    Curvature_Quark_F2H1.at(5).at(6).at(3)  = C_MassBottom * II * Vub / v2;
    Curvature_Quark_F2H1.at(5).at(7).at(2)  = C_MassBottom * Vcb / v2;
    Curvature_Quark_F2H1.at(5).at(7).at(3)  = C_MassBottom * II * Vcb / v2;
    Curvature_Quark_F2H1.at(5).at(8).at(2)  = C_MassBottom * Vtb / v2;
    Curvature_Quark_F2H1.at(5).at(8).at(3)  = C_MassBottom * II * Vtb / v2;
    Curvature_Quark_F2H1.at(5).at(11).at(6) = C_MassBottom / v2;
    Curvature_Quark_F2H1.at(5).at(11).at(7) = C_MassBottom * II / v2;
    Curvature_Quark_F2H1.at(6).at(0).at(6)  = C_MassUp / v2;
    Curvature_Quark_F2H1.at(6).at(0).at(7)  = -C_MassUp * II / v2;
    Curvature_Quark_F2H1.at(6).at(3).at(2)  = C_MassDown * Vud / v2;
    Curvature_Quark_F2H1.at(6).at(3).at(3)  = C_MassDown * II * Vud / v2;
    Curvature_Quark_F2H1.at(6).at(4).at(2)  = C_MassStrange * Vus / v2;
    Curvature_Quark_F2H1.at(6).at(4).at(3)  = C_MassStrange * II * Vus / v2;
    Curvature_Quark_F2H1.at(6).at(5).at(2)  = C_MassBottom * Vub / v2;
    Curvature_Quark_F2H1.at(6).at(5).at(3)  = C_MassBottom * II * Vub / v2;
    Curvature_Quark_F2H1.at(7).at(1).at(6)  = C_MassCharm / v2;
    Curvature_Quark_F2H1.at(7).at(1).at(7)  = -C_MassCharm * II / v2;
    Curvature_Quark_F2H1.at(7).at(3).at(2)  = C_MassDown * Vcd / v2;
    Curvature_Quark_F2H1.at(7).at(3).at(3)  = C_MassDown * II * Vcd / v2;
    Curvature_Quark_F2H1.at(7).at(4).at(2)  = C_MassStrange * Vcs / v2;
    Curvature_Quark_F2H1.at(7).at(4).at(3)  = C_MassStrange * II * Vcs / v2;
    Curvature_Quark_F2H1.at(7).at(5).at(2)  = C_MassBottom * Vcb / v2;
    Curvature_Quark_F2H1.at(7).at(5).at(3)  = C_MassBottom * II * Vcb / v2;
    Curvature_Quark_F2H1.at(8).at(2).at(6)  = C_MassTop / v2;
    Curvature_Quark_F2H1.at(8).at(2).at(7)  = -C_MassTop * II / v2;
    Curvature_Quark_F2H1.at(8).at(3).at(2)  = C_MassDown * Vtd / v2;
    Curvature_Quark_F2H1.at(8).at(3).at(3)  = C_MassDown * II * Vtd / v2;
    Curvature_Quark_F2H1.at(8).at(4).at(2)  = C_MassStrange * Vts / v2;
    Curvature_Quark_F2H1.at(8).at(4).at(3)  = C_MassStrange * II * Vts / v2;
    Curvature_Quark_F2H1.at(8).at(5).at(2)  = C_MassBottom * Vtb / v2;
    Curvature_Quark_F2H1.at(8).at(5).at(3)  = C_MassBottom * II * Vtb / v2;
    Curvature_Quark_F2H1.at(9).at(0).at(2)  = -C_MassUp * conj(Vud) / v2;
    Curvature_Quark_F2H1.at(9).at(0).at(3)  = C_MassUp * II * conj(Vud) / v2;
    Curvature_Quark_F2H1.at(9).at(1).at(2)  = -C_MassCharm * conj(Vcd) / v2;
    Curvature_Quark_F2H1.at(9).at(1).at(3)  = C_MassCharm * II * conj(Vcd) / v2;
    Curvature_Quark_F2H1.at(9).at(2).at(2)  = -C_MassTop * conj(Vtd) / v2;
    Curvature_Quark_F2H1.at(9).at(2).at(3)  = C_MassTop * II * conj(Vtd) / v2;
    Curvature_Quark_F2H1.at(9).at(3).at(6)  = C_MassDown / v2;
    Curvature_Quark_F2H1.at(9).at(3).at(7)  = C_MassDown * II / v2;
    Curvature_Quark_F2H1.at(10).at(0).at(2) = -C_MassUp * conj(Vus) / v2;
    Curvature_Quark_F2H1.at(10).at(0).at(3) = C_MassUp * II * conj(Vus) / v2;
    Curvature_Quark_F2H1.at(10).at(1).at(2) = -C_MassCharm * conj(Vcs) / v2;
    Curvature_Quark_F2H1.at(10).at(1).at(3) = C_MassCharm * II * conj(Vcs) / v2;
    Curvature_Quark_F2H1.at(10).at(2).at(2) = -C_MassTop * conj(Vts) / v2;
    Curvature_Quark_F2H1.at(10).at(2).at(3) = C_MassTop * II * conj(Vts) / v2;
    Curvature_Quark_F2H1.at(10).at(4).at(6) = C_MassStrange / v2;
    Curvature_Quark_F2H1.at(10).at(4).at(7) = C_MassStrange * II / v2;
    Curvature_Quark_F2H1.at(11).at(0).at(2) = -C_MassUp * conj(Vub) / v2;
    Curvature_Quark_F2H1.at(11).at(0).at(3) = C_MassUp * II * conj(Vub) / v2;
    Curvature_Quark_F2H1.at(11).at(1).at(2) = -C_MassCharm * conj(Vcb) / v2;
    Curvature_Quark_F2H1.at(11).at(1).at(3) = C_MassCharm * II * conj(Vcb) / v2;
    Curvature_Quark_F2H1.at(11).at(2).at(2) = -C_MassTop * conj(Vtb) / v2;
    Curvature_Quark_F2H1.at(11).at(2).at(3) = C_MassTop * II * conj(Vtb) / v2;
    Curvature_Quark_F2H1.at(11).at(5).at(6) = C_MassBottom / v2;
    Curvature_Quark_F2H1.at(11).at(5).at(7) = C_MassBottom * II / v2;
    // End of Quark interaction tensors
  };

  auto LeptonTypeII = [&]() {
    // Begin of Lepton interaction tensors
    Curvature_Lepton_F2H1.at(0).at(3).at(4) = C_MassElectron / v1;
    Curvature_Lepton_F2H1.at(0).at(3).at(5) = C_MassElectron * II / v1;
    Curvature_Lepton_F2H1.at(1).at(4).at(4) = C_MassMu / v1;
    Curvature_Lepton_F2H1.at(1).at(4).at(5) = C_MassMu * II / v1;
    Curvature_Lepton_F2H1.at(2).at(5).at(4) = C_MassTau / v1;
    Curvature_Lepton_F2H1.at(2).at(5).at(5) = C_MassTau * II / v1;
    Curvature_Lepton_F2H1.at(3).at(0).at(4) = C_MassElectron / v1;
    Curvature_Lepton_F2H1.at(3).at(0).at(5) = C_MassElectron * II / v1;
    Curvature_Lepton_F2H1.at(3).at(6).at(0) = C_MassElectron / v1;
    Curvature_Lepton_F2H1.at(3).at(6).at(1) = C_MassElectron * II / v1;
    Curvature_Lepton_F2H1.at(4).at(1).at(4) = C_MassMu / v1;
    Curvature_Lepton_F2H1.at(4).at(1).at(5) = C_MassMu * II / v1;
    Curvature_Lepton_F2H1.at(4).at(7).at(0) = C_MassMu / v1;
    Curvature_Lepton_F2H1.at(4).at(7).at(1) = C_MassMu * II / v1;
    Curvature_Lepton_F2H1.at(5).at(2).at(4) = C_MassTau / v1;
    Curvature_Lepton_F2H1.at(5).at(2).at(5) = C_MassTau * II / v1;
    Curvature_Lepton_F2H1.at(5).at(8).at(0) = C_MassTau / v1;
    Curvature_Lepton_F2H1.at(5).at(8).at(1) = C_MassTau * II / v1;
    Curvature_Lepton_F2H1.at(6).at(3).at(0) = C_MassElectron / v1;
    Curvature_Lepton_F2H1.at(6).at(3).at(1) = C_MassElectron * II / v1;
    Curvature_Lepton_F2H1.at(7).at(4).at(0) = C_MassMu / v1;
    Curvature_Lepton_F2H1.at(7).at(4).at(1) = C_MassMu * II / v1;
    Curvature_Lepton_F2H1.at(8).at(5).at(0) = C_MassTau / v1;
    Curvature_Lepton_F2H1.at(8).at(5).at(1) = C_MassTau * II / v1;
    // End of Lepton interaction tensors
  };
  auto QuarkTypeII = [&]() {
    // Begin of Quark interaction tensors
    Curvature_Quark_F2H1.at(0).at(6).at(6)  = C_MassUp / v2;
    Curvature_Quark_F2H1.at(0).at(6).at(7)  = -C_MassUp * II / v2;
    Curvature_Quark_F2H1.at(0).at(9).at(0)  = C_MassDown * Vud / v1;
    Curvature_Quark_F2H1.at(0).at(9).at(1)  = C_MassDown * II * Vud / v1;
    Curvature_Quark_F2H1.at(0).at(10).at(0) = C_MassStrange * Vus / v1;
    Curvature_Quark_F2H1.at(0).at(10).at(1) = C_MassStrange * II * Vus / v1;
    Curvature_Quark_F2H1.at(0).at(11).at(0) = C_MassBottom * Vub / v1;
    Curvature_Quark_F2H1.at(0).at(11).at(1) = C_MassBottom * II * Vub / v1;
    Curvature_Quark_F2H1.at(1).at(7).at(6)  = C_MassCharm / v2;
    Curvature_Quark_F2H1.at(1).at(7).at(7)  = -C_MassCharm * II / v2;
    Curvature_Quark_F2H1.at(1).at(9).at(0)  = C_MassDown * Vcd / v1;
    Curvature_Quark_F2H1.at(1).at(9).at(1)  = C_MassDown * II * Vcd / v1;
    Curvature_Quark_F2H1.at(1).at(10).at(0) = C_MassStrange * Vcs / v1;
    Curvature_Quark_F2H1.at(1).at(10).at(1) = C_MassStrange * II * Vcs / v1;
    Curvature_Quark_F2H1.at(1).at(11).at(0) = C_MassBottom * Vcb / v1;
    Curvature_Quark_F2H1.at(1).at(11).at(1) = C_MassBottom * II * Vcb / v1;
    Curvature_Quark_F2H1.at(2).at(8).at(6)  = C_MassTop / v2;
    Curvature_Quark_F2H1.at(2).at(8).at(7)  = -C_MassTop * II / v2;
    Curvature_Quark_F2H1.at(2).at(9).at(0)  = C_MassDown * Vtd / v1;
    Curvature_Quark_F2H1.at(2).at(9).at(1)  = C_MassDown * II * Vtd / v1;
    Curvature_Quark_F2H1.at(2).at(10).at(0) = C_MassStrange * Vts / v1;
    Curvature_Quark_F2H1.at(2).at(10).at(1) = C_MassStrange * II * Vts / v1;
    Curvature_Quark_F2H1.at(2).at(11).at(0) = C_MassBottom * Vtb / v1;
    Curvature_Quark_F2H1.at(2).at(11).at(1) = C_MassBottom * II * Vtb / v1;
    Curvature_Quark_F2H1.at(3).at(6).at(6)  = -C_MassUp * conj(Vud) / v2;
    Curvature_Quark_F2H1.at(3).at(6).at(7)  = C_MassUp * II * conj(Vud) / v2;
    Curvature_Quark_F2H1.at(3).at(7).at(6)  = -C_MassCharm * conj(Vcd) / v2;
    Curvature_Quark_F2H1.at(3).at(7).at(7)  = C_MassCharm * II * conj(Vcd) / v2;
    Curvature_Quark_F2H1.at(3).at(8).at(6)  = -C_MassTop * conj(Vtd) / v2;
    Curvature_Quark_F2H1.at(3).at(8).at(7)  = C_MassTop * II * conj(Vtd) / v2;
    Curvature_Quark_F2H1.at(3).at(9).at(0)  = C_MassDown / v1;
    Curvature_Quark_F2H1.at(3).at(9).at(1)  = C_MassDown * II / v1;
    Curvature_Quark_F2H1.at(4).at(6).at(6)  = -C_MassUp * conj(Vus) / v2;
    Curvature_Quark_F2H1.at(4).at(6).at(7)  = C_MassUp * II * conj(Vus) / v2;
    Curvature_Quark_F2H1.at(4).at(7).at(6)  = -C_MassCharm * conj(Vcs) / v2;
    Curvature_Quark_F2H1.at(4).at(7).at(7)  = C_MassCharm * II * conj(Vcs) / v2;
    Curvature_Quark_F2H1.at(4).at(8).at(6)  = -C_MassTop * conj(Vts) / v2;
    Curvature_Quark_F2H1.at(4).at(8).at(7)  = C_MassTop * II * conj(Vts) / v2;
    Curvature_Quark_F2H1.at(4).at(10).at(0) = C_MassStrange / v1;
    Curvature_Quark_F2H1.at(4).at(10).at(1) = C_MassStrange * II / v1;
    Curvature_Quark_F2H1.at(5).at(6).at(6)  = -C_MassUp * conj(Vub) / v2;
    Curvature_Quark_F2H1.at(5).at(6).at(7)  = C_MassUp * II * conj(Vub) / v2;
    Curvature_Quark_F2H1.at(5).at(7).at(6)  = -C_MassCharm * conj(Vcb) / v2;
    Curvature_Quark_F2H1.at(5).at(7).at(7)  = C_MassCharm * II * conj(Vcb) / v2;
    Curvature_Quark_F2H1.at(5).at(8).at(6)  = -C_MassTop * conj(Vtb) / v2;
    Curvature_Quark_F2H1.at(5).at(8).at(7)  = C_MassTop * II * conj(Vtb) / v2;
    Curvature_Quark_F2H1.at(5).at(11).at(0) = C_MassBottom / v1;
    Curvature_Quark_F2H1.at(5).at(11).at(1) = C_MassBottom * II / v1;
    Curvature_Quark_F2H1.at(6).at(0).at(6)  = C_MassUp / v2;
    Curvature_Quark_F2H1.at(6).at(0).at(7)  = -C_MassUp * II / v2;
    Curvature_Quark_F2H1.at(6).at(3).at(6)  = -C_MassUp * conj(Vud) / v2;
    Curvature_Quark_F2H1.at(6).at(3).at(7)  = C_MassUp * II * conj(Vud) / v2;
    Curvature_Quark_F2H1.at(6).at(4).at(6)  = -C_MassUp * conj(Vus) / v2;
    Curvature_Quark_F2H1.at(6).at(4).at(7)  = C_MassUp * II * conj(Vus) / v2;
    Curvature_Quark_F2H1.at(6).at(5).at(6)  = -C_MassUp * conj(Vub) / v2;
    Curvature_Quark_F2H1.at(6).at(5).at(7)  = C_MassUp * II * conj(Vub) / v2;
    Curvature_Quark_F2H1.at(7).at(1).at(6)  = C_MassCharm / v2;
    Curvature_Quark_F2H1.at(7).at(1).at(7)  = -C_MassCharm * II / v2;
    Curvature_Quark_F2H1.at(7).at(3).at(6)  = -C_MassCharm * conj(Vcd) / v2;
    Curvature_Quark_F2H1.at(7).at(3).at(7)  = C_MassCharm * II * conj(Vcd) / v2;
    Curvature_Quark_F2H1.at(7).at(4).at(6)  = -C_MassCharm * conj(Vcs) / v2;
    Curvature_Quark_F2H1.at(7).at(4).at(7)  = C_MassCharm * II * conj(Vcs) / v2;
    Curvature_Quark_F2H1.at(7).at(5).at(6)  = -C_MassCharm * conj(Vcb) / v2;
    Curvature_Quark_F2H1.at(7).at(5).at(7)  = C_MassCharm * II * conj(Vcb) / v2;
    Curvature_Quark_F2H1.at(8).at(2).at(6)  = C_MassTop / v2;
    Curvature_Quark_F2H1.at(8).at(2).at(7)  = -C_MassTop * II / v2;
    Curvature_Quark_F2H1.at(8).at(3).at(6)  = -C_MassTop * conj(Vtd) / v2;
    Curvature_Quark_F2H1.at(8).at(3).at(7)  = C_MassTop * II * conj(Vtd) / v2;
    Curvature_Quark_F2H1.at(8).at(4).at(6)  = -C_MassTop * conj(Vts) / v2;
    Curvature_Quark_F2H1.at(8).at(4).at(7)  = C_MassTop * II * conj(Vts) / v2;
    Curvature_Quark_F2H1.at(8).at(5).at(6)  = -C_MassTop * conj(Vtb) / v2;
    Curvature_Quark_F2H1.at(8).at(5).at(7)  = C_MassTop * II * conj(Vtb) / v2;
    Curvature_Quark_F2H1.at(9).at(0).at(0)  = C_MassDown * Vud / v1;
    Curvature_Quark_F2H1.at(9).at(0).at(1)  = C_MassDown * II * Vud / v1;
    Curvature_Quark_F2H1.at(9).at(1).at(0)  = C_MassDown * Vcd / v1;
    Curvature_Quark_F2H1.at(9).at(1).at(1)  = C_MassDown * II * Vcd / v1;
    Curvature_Quark_F2H1.at(9).at(2).at(0)  = C_MassDown * Vtd / v1;
    Curvature_Quark_F2H1.at(9).at(2).at(1)  = C_MassDown * II * Vtd / v1;
    Curvature_Quark_F2H1.at(9).at(3).at(0)  = C_MassDown / v1;
    Curvature_Quark_F2H1.at(9).at(3).at(1)  = C_MassDown * II / v1;
    Curvature_Quark_F2H1.at(10).at(0).at(0) = C_MassStrange * Vus / v1;
    Curvature_Quark_F2H1.at(10).at(0).at(1) = C_MassStrange * II * Vus / v1;
    Curvature_Quark_F2H1.at(10).at(1).at(0) = C_MassStrange * Vcs / v1;
    Curvature_Quark_F2H1.at(10).at(1).at(1) = C_MassStrange * II * Vcs / v1;
    Curvature_Quark_F2H1.at(10).at(2).at(0) = C_MassStrange * Vts / v1;
    Curvature_Quark_F2H1.at(10).at(2).at(1) = C_MassStrange * II * Vts / v1;
    Curvature_Quark_F2H1.at(10).at(4).at(0) = C_MassStrange / v1;
    Curvature_Quark_F2H1.at(10).at(4).at(1) = C_MassStrange * II / v1;
    Curvature_Quark_F2H1.at(11).at(0).at(0) = C_MassBottom * Vub / v1;
    Curvature_Quark_F2H1.at(11).at(0).at(1) = C_MassBottom * II * Vub / v1;
    Curvature_Quark_F2H1.at(11).at(1).at(0) = C_MassBottom * Vcb / v1;
    Curvature_Quark_F2H1.at(11).at(1).at(1) = C_MassBottom * II * Vcb / v1;
    Curvature_Quark_F2H1.at(11).at(2).at(0) = C_MassBottom * Vtb / v1;
    Curvature_Quark_F2H1.at(11).at(2).at(1) = C_MassBottom * II * Vtb / v1;
    Curvature_Quark_F2H1.at(11).at(5).at(0) = C_MassBottom / v1;
    Curvature_Quark_F2H1.at(11).at(5).at(1) = C_MassBottom * II / v1;
    // End of Quark interaction tensors
  };

  auto LeptonFlipped = [&]() { LeptonTypeI(); };
  auto QuarkFlipped  = [&]() { QuarkTypeII(); };

  auto LeptonLS = [&]() { LeptonTypeII(); };
  auto QuarkLS  = [&]() { QuarkTypeI(); };

  switch (Type)
  {
  case THDMType::TypeI:
  {
    LeptonTypeI();
    QuarkTypeI();
    break;
  }
  case THDMType::TypeII:
  {
    LeptonTypeII();
    QuarkTypeII();
    break;
  }
  case THDMType::LeptonSpecific:
  {
    LeptonLS();
    QuarkLS();
    break;
  }
  case THDMType::Flipped:
  {
    LeptonFlipped();
    QuarkFlipped();
    break;
  }
  case THDMType::Invalid:
  default: throw std::runtime_error("Invalid 2HDM Type");
  }
}

bool C2HDMSympy::CalculateDebyeSimplified()
{
  return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs
   * mass matrix and implement your formula here and return true. The vector
   * is given by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool C2HDMSympy::CalculateDebyeGaugeSimplified()
{
  /*
   * Use this function if you calculated the Debye corrections to the gauge
   * mass matrix and implement your formula here and return true. The vector
   * is given by DebyeGauge[NGauge][NGauge]
   */

  return false;
}
double C2HDMSympy::VTreeSimplified(const std::vector<double> &v) const
{
  if (not UseVTreeSimplified) return 0;
  double res = 0;

  double wCB = v.at(2);
  double w1  = v.at(4);
  double w2  = v.at(6);
  double wCP = v.at(7);
  res        = -1.0 / 2.0 * Imlambda5 * std::pow(w1, 2) * w2 * wCP -
        1.0 / 2.0 * Imlambda6 * std::pow(w1, 3) * wCP -
        1.0 / 2.0 * Imlambda7 * w1 * std::pow(w2, 2) * wCP -
        1.0 / 2.0 * Imlambda7 * w1 * std::pow(wCB, 2) * wCP -
        1.0 / 2.0 * Imlambda7 * w1 * std::pow(wCP, 3) + Imm12sq * w1 * wCP +
        (1.0 / 4.0) * Relambda5 * std::pow(w1, 2) * std::pow(w2, 2) -
        1.0 / 4.0 * Relambda5 * std::pow(w1, 2) * std::pow(wCP, 2) +
        (1.0 / 2.0) * Relambda6 * std::pow(w1, 3) * w2 +
        (1.0 / 2.0) * Relambda7 * w1 * std::pow(w2, 3) +
        (1.0 / 2.0) * Relambda7 * w1 * w2 * std::pow(wCB, 2) +
        (1.0 / 2.0) * Relambda7 * w1 * w2 * std::pow(wCP, 2) -
        Rem12sq * w1 * w2 + (1.0 / 8.0) * lambda1 * std::pow(w1, 4) +
        (1.0 / 8.0) * lambda2 * std::pow(w2, 4) +
        (1.0 / 4.0) * lambda2 * std::pow(w2, 2) * std::pow(wCB, 2) +
        (1.0 / 4.0) * lambda2 * std::pow(w2, 2) * std::pow(wCP, 2) +
        (1.0 / 8.0) * lambda2 * std::pow(wCB, 4) +
        (1.0 / 4.0) * lambda2 * std::pow(wCB, 2) * std::pow(wCP, 2) +
        (1.0 / 8.0) * lambda2 * std::pow(wCP, 4) +
        (1.0 / 4.0) * lambda3 * std::pow(w1, 2) * std::pow(w2, 2) +
        (1.0 / 4.0) * lambda3 * std::pow(w1, 2) * std::pow(wCB, 2) +
        (1.0 / 4.0) * lambda3 * std::pow(w1, 2) * std::pow(wCP, 2) +
        (1.0 / 4.0) * lambda4 * std::pow(w1, 2) * std::pow(w2, 2) +
        (1.0 / 4.0) * lambda4 * std::pow(w1, 2) * std::pow(wCP, 2) +
        (1.0 / 2.0) * m11sq * std::pow(w1, 2) +
        (1.0 / 2.0) * m22sq * std::pow(w2, 2) +
        (1.0 / 2.0) * m22sq * std::pow(wCB, 2) +
        (1.0 / 2.0) * m22sq * std::pow(wCP, 2);

  return res;
}

double C2HDMSympy::VCounterSimplified(const std::vector<double> &v) const
{
  if (not UseVCounterSimplified) return 0;
  double res = 0;
  double wCB = v.at(2);
  double w1  = v.at(4);
  double w2  = v.at(6);
  double wCP = v.at(7);
  res        = -1.0 / 2.0 * dImlambda5 * std::pow(w1, 2) * w2 * wCP -
        1.0 / 2.0 * dImlambda6 * std::pow(w1, 3) * wCP -
        1.0 / 2.0 * dImlambda7 * w1 * std::pow(w2, 2) * wCP -
        1.0 / 2.0 * dImlambda7 * w1 * std::pow(wCB, 2) * wCP -
        1.0 / 2.0 * dImlambda7 * w1 * std::pow(wCP, 3) + dImm12sq * w1 * wCP +
        (1.0 / 4.0) * dRelambda5 * std::pow(w1, 2) * std::pow(w2, 2) -
        1.0 / 4.0 * dRelambda5 * std::pow(w1, 2) * std::pow(wCP, 2) +
        (1.0 / 2.0) * dRelambda6 * std::pow(w1, 3) * w2 +
        (1.0 / 2.0) * dRelambda7 * w1 * std::pow(w2, 3) +
        (1.0 / 2.0) * dRelambda7 * w1 * w2 * std::pow(wCB, 2) +
        (1.0 / 2.0) * dRelambda7 * w1 * w2 * std::pow(wCP, 2) -
        dRem12sq * w1 * w2 + dT3 * wCB + dT5 * w1 + dT7 * w2 + dT8 * wCP +
        (1.0 / 8.0) * dlambda1 * std::pow(w1, 4) +
        (1.0 / 8.0) * dlambda2 * std::pow(w2, 4) +
        (1.0 / 4.0) * dlambda2 * std::pow(w2, 2) * std::pow(wCB, 2) +
        (1.0 / 4.0) * dlambda2 * std::pow(w2, 2) * std::pow(wCP, 2) +
        (1.0 / 8.0) * dlambda2 * std::pow(wCB, 4) +
        (1.0 / 4.0) * dlambda2 * std::pow(wCB, 2) * std::pow(wCP, 2) +
        (1.0 / 8.0) * dlambda2 * std::pow(wCP, 4) +
        (1.0 / 4.0) * dlambda3 * std::pow(w1, 2) * std::pow(w2, 2) +
        (1.0 / 4.0) * dlambda3 * std::pow(w1, 2) * std::pow(wCB, 2) +
        (1.0 / 4.0) * dlambda3 * std::pow(w1, 2) * std::pow(wCP, 2) +
        (1.0 / 4.0) * dlambda4 * std::pow(w1, 2) * std::pow(w2, 2) +
        (1.0 / 4.0) * dlambda4 * std::pow(w1, 2) * std::pow(wCP, 2) +
        (1.0 / 2.0) * dm11sq * std::pow(w1, 2) +
        (1.0 / 2.0) * dm22sq * std::pow(w2, 2) +
        (1.0 / 2.0) * dm22sq * std::pow(wCB, 2) +
        (1.0 / 2.0) * dm22sq * std::pow(wCP, 2);

  return res;
}

void C2HDMSympy::Debugging(const std::vector<double> &input,
                           std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
