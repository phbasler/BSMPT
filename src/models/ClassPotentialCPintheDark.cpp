/*
 * ClassPotentialCPintheDark.cpp
 */
// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner, Jonas
// Müller and Lisa Biermann
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/IterativeLinearSolvers"
#include <BSMPT/models/SMparam.h> // for SMConstants.C_vev0, SMConstants.C_MassTop, SMConstants.C_g
#include <algorithm> // for max, copy
#include <iomanip>
#include <iostream> // for operator<<, endl, basic_o...
#include <memory>   // for allocator_traits<>::value...
#include <stddef.h> // for std::size_t

#include <BSMPT/models/ClassPotentialCPintheDark.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>

#include <optional>

using namespace Eigen;

namespace BSMPT
{
namespace Models
{

Class_Potential_CPintheDark::Class_Potential_CPintheDark(
    const ISMConstants &smConstants)
    : Class_Potential_Origin(smConstants)
{
  Model         = ModelID::ModelIDs::CPINTHEDARK;
  NNeutralHiggs = 5; // number of neutral Higgs bosons at T = 0
  NChargedHiggs = 4; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar = 13;   // number of parameters in the tree-Level Lagrangian AFTER using
               // tadpole equations
  nParCT = 19; // number of parameters in the counterterm potential

  nVEV = 5; // number of VEVs to minimize the potential

  NHiggs = NNeutralHiggs + NChargedHiggs;

  VevOrder.resize(nVEV);
  VevOrder[0] = 2; // omegaCB
  VevOrder[1] = 4; // omega1
  VevOrder[2] = 6; // omega2
  VevOrder[3] = 7; // omegaCP
  VevOrder[4] = 8; // omegaS

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = true;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = true;
}

Class_Potential_CPintheDark::~Class_Potential_CPintheDark()
{
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_Potential_CPintheDark::addLegendCT() const
{
  std::vector<std::string> labels;

  labels.push_back("dm11s");
  labels.push_back("dm22s");
  labels.push_back("dmSs");
  labels.push_back("dReA");
  labels.push_back("dImA");
  labels.push_back("dL1");
  labels.push_back("dL2");
  labels.push_back("dL3");
  labels.push_back("dL4");
  labels.push_back("dL5");
  labels.push_back("dL6");
  labels.push_back("dL7");
  labels.push_back("dL8");
  labels.push_back("dT1");
  labels.push_back("dT2");
  labels.push_back("dTCP");
  labels.push_back("dTCB");
  labels.push_back("dTS");
  labels.push_back("dImL5");

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * input file
 */
std::vector<std::string> Class_Potential_CPintheDark::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c");     // Label for the critical temperature
  labels.push_back("v_c");     // Label for the critical vev
  labels.push_back("v_c/T_c"); // Label for xi_c
  // out += "VEV order";
  labels.push_back("omega_{CB}(T_c)");
  labels.push_back("omega_1(T_c)");
  labels.push_back("omega_2(T_c)");
  labels.push_back("omega_{CP}(T_c)");
  labels.push_back("omega_S(T_c)");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * Higgs couplings. Use this to complement the legend of the given input file
 */
std::vector<std::string>
Class_Potential_CPintheDark::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;

  // mass basis
  particles.push_back("G^+");
  particles.push_back("G^-");
  particles.push_back("H^+");
  particles.push_back("H^-");
  particles.push_back("h");
  particles.push_back("G^0");
  particles.push_back("h_1");
  particles.push_back("h_2");
  particles.push_back("h_3");

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
std::vector<std::string> Class_Potential_CPintheDark::addLegendVEV() const
{
  std::vector<std::string> labels;
  // out = "Your VEV order";
  labels.push_back("omega_{CB}");
  labels.push_back("omega_1");
  labels.push_back("omega_2");
  labels.push_back("omega_{CP}");
  labels.push_back("omega_S");
  return labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_Potential_CPintheDark::ReadAndSet(const std::string &linestr,
                                             std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  double lreA = 0, limA = 0, lL1 = 0, lL2 = 0, lL3 = 0, lL4 = 0, lL5 = 0,
         lL6 = 0, lL7 = 0, lL8 = 0, lm11s = 0, lm22s = 0, lmSs = 0;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 27; k++)
  {
    ss >> tmp;
    if (k == 15)
      lL1 = tmp;
    else if (k == 16)
      lL2 = tmp;
    else if (k == 17)
      lL3 = tmp;
    else if (k == 18)
      lL4 = tmp;
    else if (k == 19)
      lL5 = tmp;
    else if (k == 20)
      lL6 = tmp;
    else if (k == 21)
      lL7 = tmp;
    else if (k == 22)
      lL8 = tmp;
    else if (k == 23)
      lreA = tmp; // Tr
    else if (k == 24)
      limA = tmp; // Ti
    else if (k == 25)
      lm11s = tmp;
    else if (k == 26)
      lm22s = tmp;
    else if (k == 27)
      lmSs = tmp;
  }

  par[0]  = lm11s;
  par[1]  = lm22s;
  par[2]  = lmSs;
  par[3]  = lreA;
  par[4]  = limA;
  par[5]  = lL1;
  par[6]  = lL2;
  par[7]  = lL3;
  par[8]  = lL4;
  par[9]  = lL5;
  par[10] = lL6;
  par[11] = lL7;
  par[12] = lL8;

  set_gen(par);
  return;
}

/**
 * Set Class Object as well as the VEV configuration
 */
void Class_Potential_CPintheDark::set_gen(const std::vector<double> &par)
{

  m11s = par[0];
  m22s = par[1];
  mSs  = par[2];
  ReA  = par[3];
  ImA  = par[4];
  L1   = par[5];
  L2   = par[6];
  L3   = par[7];
  L4   = par[8];
  L5   = par[9];
  L6   = par[10];
  L7   = par[11];
  L8   = par[12];

  // set vev
  v1 = SMConstants.C_vev0;

  scale = SMConstants.C_vev0; // renormalisation scale is set to the SM VEV

  vevTreeMin.resize(nVEV);
  vevTree.resize(NHiggs);
  // set the vector vevTreeMin. vevTree will then be set by the
  // function MinimizeOrderVEV
  vevTreeMin[0] = 0;
  vevTreeMin[1] = v1;
  vevTreeMin[2] = 0;
  vevTreeMin[3] = 0;
  vevTreeMin[4] = 0;

  vevTree = MinimizeOrderVEV(vevTreeMin);
  if (!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the
 * entries of Curvature_Higgs_CT_L1 to Curvature_Higgs_CT_L4.
 */
void Class_Potential_CPintheDark::set_CT_Pot_Par(const std::vector<double> &par)
{
  dm11s = par[0];
  dm22s = par[1];
  dmSs  = par[2];
  dReA  = par[3];
  dImA  = par[4];
  dL1   = par[5];
  dL2   = par[6];
  dL3   = par[7];
  dL4   = par[8];
  dL5   = par[9];
  dL6   = par[10];
  dL7   = par[11];
  dL8   = par[12];
  dTCB  = par[13];
  dT1   = par[14];
  dT2   = par[15];
  dTCP  = par[16];
  dTS   = par[17];
  dImL5 = par[18];

  // assign the non-zero entries
  Curvature_Higgs_CT_L1[2] = dTCB;
  Curvature_Higgs_CT_L1[4] = dT1;
  Curvature_Higgs_CT_L1[6] = dT2;
  Curvature_Higgs_CT_L1[7] = dTCP;
  Curvature_Higgs_CT_L1[8] = dTS;

  Curvature_Higgs_CT_L2[0][0] = dm11s;
  Curvature_Higgs_CT_L2[1][1] = dm11s;
  Curvature_Higgs_CT_L2[2][2] = dm22s;
  Curvature_Higgs_CT_L2[3][3] = dm22s;
  Curvature_Higgs_CT_L2[4][4] = dm11s;
  Curvature_Higgs_CT_L2[5][5] = dm11s;
  Curvature_Higgs_CT_L2[6][6] = dm22s;
  Curvature_Higgs_CT_L2[7][7] = dm22s;
  Curvature_Higgs_CT_L2[8][8] = dmSs;

  Curvature_Higgs_CT_L3[0][2][8] = dReA;
  Curvature_Higgs_CT_L3[0][3][8] = -dImA;
  Curvature_Higgs_CT_L3[0][8][2] = dReA;
  Curvature_Higgs_CT_L3[0][8][3] = -dImA;
  Curvature_Higgs_CT_L3[1][2][8] = dImA;
  Curvature_Higgs_CT_L3[1][3][8] = dReA;
  Curvature_Higgs_CT_L3[1][8][2] = dImA;
  Curvature_Higgs_CT_L3[1][8][3] = dReA;
  Curvature_Higgs_CT_L3[2][0][8] = dReA;
  Curvature_Higgs_CT_L3[2][1][8] = dImA;
  Curvature_Higgs_CT_L3[2][8][0] = dReA;
  Curvature_Higgs_CT_L3[2][8][1] = dImA;
  Curvature_Higgs_CT_L3[3][0][8] = -dImA;
  Curvature_Higgs_CT_L3[3][1][8] = dReA;
  Curvature_Higgs_CT_L3[3][8][0] = -dImA;
  Curvature_Higgs_CT_L3[3][8][1] = dReA;
  Curvature_Higgs_CT_L3[4][6][8] = dReA;
  Curvature_Higgs_CT_L3[4][7][8] = -dImA;
  Curvature_Higgs_CT_L3[4][8][6] = dReA;
  Curvature_Higgs_CT_L3[4][8][7] = -dImA;
  Curvature_Higgs_CT_L3[5][6][8] = dImA;
  Curvature_Higgs_CT_L3[5][7][8] = dReA;
  Curvature_Higgs_CT_L3[5][8][6] = dImA;
  Curvature_Higgs_CT_L3[5][8][7] = dReA;
  Curvature_Higgs_CT_L3[6][4][8] = dReA;
  Curvature_Higgs_CT_L3[6][5][8] = dImA;
  Curvature_Higgs_CT_L3[6][8][4] = dReA;
  Curvature_Higgs_CT_L3[6][8][5] = dImA;
  Curvature_Higgs_CT_L3[7][4][8] = -dImA;
  Curvature_Higgs_CT_L3[7][5][8] = dReA;
  Curvature_Higgs_CT_L3[7][8][4] = -dImA;
  Curvature_Higgs_CT_L3[7][8][5] = dReA;
  Curvature_Higgs_CT_L3[8][0][2] = dReA;
  Curvature_Higgs_CT_L3[8][0][3] = -dImA;
  Curvature_Higgs_CT_L3[8][1][2] = dImA;
  Curvature_Higgs_CT_L3[8][1][3] = dReA;
  Curvature_Higgs_CT_L3[8][2][0] = dReA;
  Curvature_Higgs_CT_L3[8][2][1] = dImA;
  Curvature_Higgs_CT_L3[8][3][0] = -dImA;
  Curvature_Higgs_CT_L3[8][3][1] = dReA;
  Curvature_Higgs_CT_L3[8][4][6] = dReA;
  Curvature_Higgs_CT_L3[8][4][7] = -dImA;
  Curvature_Higgs_CT_L3[8][5][6] = dImA;
  Curvature_Higgs_CT_L3[8][5][7] = dReA;
  Curvature_Higgs_CT_L3[8][6][4] = dReA;
  Curvature_Higgs_CT_L3[8][6][5] = dImA;
  Curvature_Higgs_CT_L3[8][7][4] = -dImA;
  Curvature_Higgs_CT_L3[8][7][5] = dReA;

  Curvature_Higgs_CT_L4[0][0][0][0] = 0.3e1 * dL1;
  Curvature_Higgs_CT_L4[0][0][1][1] = dL1;
  Curvature_Higgs_CT_L4[0][0][2][2] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[0][0][2][3] = -dImL5;
  Curvature_Higgs_CT_L4[0][0][3][2] = -dImL5;
  Curvature_Higgs_CT_L4[0][0][3][3] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[0][0][4][4] = dL1;
  Curvature_Higgs_CT_L4[0][0][5][5] = dL1;
  Curvature_Higgs_CT_L4[0][0][6][6] = dL3;
  Curvature_Higgs_CT_L4[0][0][7][7] = dL3;
  Curvature_Higgs_CT_L4[0][0][8][8] = dL7;
  Curvature_Higgs_CT_L4[0][1][0][1] = dL1;
  Curvature_Higgs_CT_L4[0][1][1][0] = dL1;
  Curvature_Higgs_CT_L4[0][1][2][2] = dImL5;
  Curvature_Higgs_CT_L4[0][1][2][3] = dL5;
  Curvature_Higgs_CT_L4[0][1][3][2] = dL5;
  Curvature_Higgs_CT_L4[0][1][3][3] = -dImL5;
  Curvature_Higgs_CT_L4[0][2][0][2] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[0][2][0][3] = -dImL5;
  Curvature_Higgs_CT_L4[0][2][1][2] = dImL5;
  Curvature_Higgs_CT_L4[0][2][1][3] = dL5;
  Curvature_Higgs_CT_L4[0][2][2][0] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[0][2][2][1] = dImL5;
  Curvature_Higgs_CT_L4[0][2][3][0] = -dImL5;
  Curvature_Higgs_CT_L4[0][2][3][1] = dL5;
  Curvature_Higgs_CT_L4[0][2][4][6] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][4][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][5][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][5][7] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][6][4] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][6][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][7][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][7][5] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][0][2] = -dImL5;
  Curvature_Higgs_CT_L4[0][3][0][3] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[0][3][1][2] = dL5;
  Curvature_Higgs_CT_L4[0][3][1][3] = -dImL5;
  Curvature_Higgs_CT_L4[0][3][2][0] = -dImL5;
  Curvature_Higgs_CT_L4[0][3][2][1] = dL5;
  Curvature_Higgs_CT_L4[0][3][3][0] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[0][3][3][1] = -dImL5;
  Curvature_Higgs_CT_L4[0][3][4][6] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][4][7] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][5][6] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][5][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][6][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][6][5] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][7][4] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][7][5] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][4][0][4] = dL1;
  Curvature_Higgs_CT_L4[0][4][2][6] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][4][2][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][4][3][6] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][4][3][7] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][4][4][0] = dL1;
  Curvature_Higgs_CT_L4[0][4][6][2] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][4][6][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][4][7][2] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][4][7][3] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][5][0][5] = dL1;
  Curvature_Higgs_CT_L4[0][5][2][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][5][2][7] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][5][3][6] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][5][3][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][5][5][0] = dL1;
  Curvature_Higgs_CT_L4[0][5][6][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][5][6][3] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][5][7][2] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][5][7][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][6][0][6] = dL3;
  Curvature_Higgs_CT_L4[0][6][2][4] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][6][2][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][6][3][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][6][3][5] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][6][4][2] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][6][4][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][6][5][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][6][5][3] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][6][6][0] = dL3;
  Curvature_Higgs_CT_L4[0][7][0][7] = dL3;
  Curvature_Higgs_CT_L4[0][7][2][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][7][2][5] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][7][3][4] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][7][3][5] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][7][4][2] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][7][4][3] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][7][5][2] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][7][5][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[0][7][7][0] = dL3;
  Curvature_Higgs_CT_L4[0][8][0][8] = dL7;
  Curvature_Higgs_CT_L4[0][8][8][0] = dL7;
  Curvature_Higgs_CT_L4[1][0][0][1] = dL1;
  Curvature_Higgs_CT_L4[1][0][1][0] = dL1;
  Curvature_Higgs_CT_L4[1][0][2][2] = dImL5;
  Curvature_Higgs_CT_L4[1][0][2][3] = dL5;
  Curvature_Higgs_CT_L4[1][0][3][2] = dL5;
  Curvature_Higgs_CT_L4[1][0][3][3] = -dImL5;
  Curvature_Higgs_CT_L4[1][1][0][0] = dL1;
  Curvature_Higgs_CT_L4[1][1][1][1] = 0.3e1 * dL1;
  Curvature_Higgs_CT_L4[1][1][2][2] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[1][1][2][3] = dImL5;
  Curvature_Higgs_CT_L4[1][1][3][2] = dImL5;
  Curvature_Higgs_CT_L4[1][1][3][3] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[1][1][4][4] = dL1;
  Curvature_Higgs_CT_L4[1][1][5][5] = dL1;
  Curvature_Higgs_CT_L4[1][1][6][6] = dL3;
  Curvature_Higgs_CT_L4[1][1][7][7] = dL3;
  Curvature_Higgs_CT_L4[1][1][8][8] = dL7;
  Curvature_Higgs_CT_L4[1][2][0][2] = dImL5;
  Curvature_Higgs_CT_L4[1][2][0][3] = dL5;
  Curvature_Higgs_CT_L4[1][2][1][2] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[1][2][1][3] = dImL5;
  Curvature_Higgs_CT_L4[1][2][2][0] = dImL5;
  Curvature_Higgs_CT_L4[1][2][2][1] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[1][2][3][0] = dL5;
  Curvature_Higgs_CT_L4[1][2][3][1] = dImL5;
  Curvature_Higgs_CT_L4[1][2][4][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][4][7] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][5][6] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][5][7] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][6][4] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][6][5] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][7][4] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][7][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][0][2] = dL5;
  Curvature_Higgs_CT_L4[1][3][0][3] = -dImL5;
  Curvature_Higgs_CT_L4[1][3][1][2] = dImL5;
  Curvature_Higgs_CT_L4[1][3][1][3] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[1][3][2][0] = dL5;
  Curvature_Higgs_CT_L4[1][3][2][1] = dImL5;
  Curvature_Higgs_CT_L4[1][3][3][0] = -dImL5;
  Curvature_Higgs_CT_L4[1][3][3][1] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[1][3][4][6] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][4][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][5][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][5][7] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][6][4] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][6][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][7][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][7][5] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][4][1][4] = dL1;
  Curvature_Higgs_CT_L4[1][4][2][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][4][2][7] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][4][3][6] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][4][3][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][4][4][1] = dL1;
  Curvature_Higgs_CT_L4[1][4][6][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][4][6][3] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][4][7][2] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][4][7][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][5][1][5] = dL1;
  Curvature_Higgs_CT_L4[1][5][2][6] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][5][2][7] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][5][3][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][5][3][7] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][5][5][1] = dL1;
  Curvature_Higgs_CT_L4[1][5][6][2] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][5][6][3] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][5][7][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][5][7][3] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][6][1][6] = dL3;
  Curvature_Higgs_CT_L4[1][6][2][4] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][6][2][5] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][6][3][4] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][6][3][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][6][4][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][6][4][3] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][6][5][2] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][6][5][3] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][6][6][1] = dL3;
  Curvature_Higgs_CT_L4[1][7][1][7] = dL3;
  Curvature_Higgs_CT_L4[1][7][2][4] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][7][2][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][7][3][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][7][3][5] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][7][4][2] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][7][4][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][7][5][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][7][5][3] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[1][7][7][1] = dL3;
  Curvature_Higgs_CT_L4[1][8][1][8] = dL7;
  Curvature_Higgs_CT_L4[1][8][8][1] = dL7;
  Curvature_Higgs_CT_L4[2][0][0][2] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[2][0][0][3] = -dImL5;
  Curvature_Higgs_CT_L4[2][0][1][2] = dImL5;
  Curvature_Higgs_CT_L4[2][0][1][3] = dL5;
  Curvature_Higgs_CT_L4[2][0][2][0] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[2][0][2][1] = dImL5;
  Curvature_Higgs_CT_L4[2][0][3][0] = -dImL5;
  Curvature_Higgs_CT_L4[2][0][3][1] = dL5;
  Curvature_Higgs_CT_L4[2][0][4][6] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][0][4][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][0][5][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][0][5][7] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][0][6][4] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][0][6][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][0][7][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][0][7][5] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][1][0][2] = dImL5;
  Curvature_Higgs_CT_L4[2][1][0][3] = dL5;
  Curvature_Higgs_CT_L4[2][1][1][2] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[2][1][1][3] = dImL5;
  Curvature_Higgs_CT_L4[2][1][2][0] = dImL5;
  Curvature_Higgs_CT_L4[2][1][2][1] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[2][1][3][0] = dL5;
  Curvature_Higgs_CT_L4[2][1][3][1] = dImL5;
  Curvature_Higgs_CT_L4[2][1][4][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][1][4][7] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][1][5][6] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][1][5][7] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][1][6][4] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][1][6][5] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][1][7][4] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][1][7][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][2][0][0] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[2][2][0][1] = dImL5;
  Curvature_Higgs_CT_L4[2][2][1][0] = dImL5;
  Curvature_Higgs_CT_L4[2][2][1][1] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[2][2][2][2] = 0.3e1 * dL2;
  Curvature_Higgs_CT_L4[2][2][3][3] = dL2;
  Curvature_Higgs_CT_L4[2][2][4][4] = dL3;
  Curvature_Higgs_CT_L4[2][2][5][5] = dL3;
  Curvature_Higgs_CT_L4[2][2][6][6] = dL2;
  Curvature_Higgs_CT_L4[2][2][7][7] = dL2;
  Curvature_Higgs_CT_L4[2][2][8][8] = dL8;
  Curvature_Higgs_CT_L4[2][3][0][0] = -dImL5;
  Curvature_Higgs_CT_L4[2][3][0][1] = dL5;
  Curvature_Higgs_CT_L4[2][3][1][0] = dL5;
  Curvature_Higgs_CT_L4[2][3][1][1] = dImL5;
  Curvature_Higgs_CT_L4[2][3][2][3] = dL2;
  Curvature_Higgs_CT_L4[2][3][3][2] = dL2;
  Curvature_Higgs_CT_L4[2][4][0][6] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][4][0][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][4][1][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][4][1][7] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][4][2][4] = dL3;
  Curvature_Higgs_CT_L4[2][4][4][2] = dL3;
  Curvature_Higgs_CT_L4[2][4][6][0] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][4][6][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][4][7][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][4][7][1] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][5][0][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][5][0][7] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][5][1][6] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][5][1][7] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][5][2][5] = dL3;
  Curvature_Higgs_CT_L4[2][5][5][2] = dL3;
  Curvature_Higgs_CT_L4[2][5][6][0] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][5][6][1] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][5][7][0] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][5][7][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][6][0][4] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][6][0][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][6][1][4] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][6][1][5] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][6][2][6] = dL2;
  Curvature_Higgs_CT_L4[2][6][4][0] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][6][4][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][6][5][0] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][6][5][1] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][6][6][2] = dL2;
  Curvature_Higgs_CT_L4[2][7][0][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][7][0][5] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][7][1][4] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][7][1][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][7][2][7] = dL2;
  Curvature_Higgs_CT_L4[2][7][4][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][7][4][1] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][7][5][0] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][7][5][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[2][7][7][2] = dL2;
  Curvature_Higgs_CT_L4[2][8][2][8] = dL8;
  Curvature_Higgs_CT_L4[2][8][8][2] = dL8;
  Curvature_Higgs_CT_L4[3][0][0][2] = -dImL5;
  Curvature_Higgs_CT_L4[3][0][0][3] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[3][0][1][2] = dL5;
  Curvature_Higgs_CT_L4[3][0][1][3] = -dImL5;
  Curvature_Higgs_CT_L4[3][0][2][0] = -dImL5;
  Curvature_Higgs_CT_L4[3][0][2][1] = dL5;
  Curvature_Higgs_CT_L4[3][0][3][0] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[3][0][3][1] = -dImL5;
  Curvature_Higgs_CT_L4[3][0][4][6] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][0][4][7] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][0][5][6] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][0][5][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][0][6][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][0][6][5] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][0][7][4] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][0][7][5] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][1][0][2] = dL5;
  Curvature_Higgs_CT_L4[3][1][0][3] = -dImL5;
  Curvature_Higgs_CT_L4[3][1][1][2] = dImL5;
  Curvature_Higgs_CT_L4[3][1][1][3] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[3][1][2][0] = dL5;
  Curvature_Higgs_CT_L4[3][1][2][1] = dImL5;
  Curvature_Higgs_CT_L4[3][1][3][0] = -dImL5;
  Curvature_Higgs_CT_L4[3][1][3][1] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[3][1][4][6] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][1][4][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][1][5][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][1][5][7] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][1][6][4] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][1][6][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][1][7][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][1][7][5] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][2][0][0] = -dImL5;
  Curvature_Higgs_CT_L4[3][2][0][1] = dL5;
  Curvature_Higgs_CT_L4[3][2][1][0] = dL5;
  Curvature_Higgs_CT_L4[3][2][1][1] = dImL5;
  Curvature_Higgs_CT_L4[3][2][2][3] = dL2;
  Curvature_Higgs_CT_L4[3][2][3][2] = dL2;
  Curvature_Higgs_CT_L4[3][3][0][0] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[3][3][0][1] = -dImL5;
  Curvature_Higgs_CT_L4[3][3][1][0] = -dImL5;
  Curvature_Higgs_CT_L4[3][3][1][1] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[3][3][2][2] = dL2;
  Curvature_Higgs_CT_L4[3][3][3][3] = 0.3e1 * dL2;
  Curvature_Higgs_CT_L4[3][3][4][4] = dL3;
  Curvature_Higgs_CT_L4[3][3][5][5] = dL3;
  Curvature_Higgs_CT_L4[3][3][6][6] = dL2;
  Curvature_Higgs_CT_L4[3][3][7][7] = dL2;
  Curvature_Higgs_CT_L4[3][3][8][8] = dL8;
  Curvature_Higgs_CT_L4[3][4][0][6] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][4][0][7] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][4][1][6] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][4][1][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][4][3][4] = dL3;
  Curvature_Higgs_CT_L4[3][4][4][3] = dL3;
  Curvature_Higgs_CT_L4[3][4][6][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][4][6][1] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][4][7][0] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][4][7][1] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][5][0][6] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][5][0][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][5][1][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][5][1][7] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][5][3][5] = dL3;
  Curvature_Higgs_CT_L4[3][5][5][3] = dL3;
  Curvature_Higgs_CT_L4[3][5][6][0] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][5][6][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][5][7][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][5][7][1] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][6][0][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][6][0][5] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][6][1][4] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][6][1][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][6][3][6] = dL2;
  Curvature_Higgs_CT_L4[3][6][4][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][6][4][1] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][6][5][0] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][6][5][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][6][6][3] = dL2;
  Curvature_Higgs_CT_L4[3][7][0][4] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][7][0][5] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][7][1][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][7][1][5] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][7][3][7] = dL2;
  Curvature_Higgs_CT_L4[3][7][4][0] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][7][4][1] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][7][5][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][7][5][1] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[3][7][7][3] = dL2;
  Curvature_Higgs_CT_L4[3][8][3][8] = dL8;
  Curvature_Higgs_CT_L4[3][8][8][3] = dL8;
  Curvature_Higgs_CT_L4[4][0][0][4] = dL1;
  Curvature_Higgs_CT_L4[4][0][2][6] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][0][2][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][0][3][6] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][0][3][7] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][0][4][0] = dL1;
  Curvature_Higgs_CT_L4[4][0][6][2] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][0][6][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][0][7][2] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][0][7][3] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][1][1][4] = dL1;
  Curvature_Higgs_CT_L4[4][1][2][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][1][2][7] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][1][3][6] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][1][3][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][1][4][1] = dL1;
  Curvature_Higgs_CT_L4[4][1][6][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][1][6][3] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][1][7][2] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][1][7][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][2][0][6] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][2][0][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][2][1][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][2][1][7] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][2][2][4] = dL3;
  Curvature_Higgs_CT_L4[4][2][4][2] = dL3;
  Curvature_Higgs_CT_L4[4][2][6][0] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][2][6][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][2][7][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][2][7][1] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][3][0][6] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][3][0][7] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][3][1][6] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][3][1][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][3][3][4] = dL3;
  Curvature_Higgs_CT_L4[4][3][4][3] = dL3;
  Curvature_Higgs_CT_L4[4][3][6][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][3][6][1] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][3][7][0] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][3][7][1] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][4][0][0] = dL1;
  Curvature_Higgs_CT_L4[4][4][1][1] = dL1;
  Curvature_Higgs_CT_L4[4][4][2][2] = dL3;
  Curvature_Higgs_CT_L4[4][4][3][3] = dL3;
  Curvature_Higgs_CT_L4[4][4][4][4] = 0.3e1 * dL1;
  Curvature_Higgs_CT_L4[4][4][5][5] = dL1;
  Curvature_Higgs_CT_L4[4][4][6][6] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[4][4][6][7] = -dImL5;
  Curvature_Higgs_CT_L4[4][4][7][6] = -dImL5;
  Curvature_Higgs_CT_L4[4][4][7][7] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[4][4][8][8] = dL7;
  Curvature_Higgs_CT_L4[4][5][4][5] = dL1;
  Curvature_Higgs_CT_L4[4][5][5][4] = dL1;
  Curvature_Higgs_CT_L4[4][5][6][6] = dImL5;
  Curvature_Higgs_CT_L4[4][5][6][7] = dL5;
  Curvature_Higgs_CT_L4[4][5][7][6] = dL5;
  Curvature_Higgs_CT_L4[4][5][7][7] = -dImL5;
  Curvature_Higgs_CT_L4[4][6][0][2] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][6][0][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][6][1][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][6][1][3] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][6][2][0] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][6][2][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][6][3][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][6][3][1] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][6][4][6] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[4][6][4][7] = -dImL5;
  Curvature_Higgs_CT_L4[4][6][5][6] = dImL5;
  Curvature_Higgs_CT_L4[4][6][5][7] = dL5;
  Curvature_Higgs_CT_L4[4][6][6][4] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[4][6][6][5] = dImL5;
  Curvature_Higgs_CT_L4[4][6][7][4] = -dImL5;
  Curvature_Higgs_CT_L4[4][6][7][5] = dL5;
  Curvature_Higgs_CT_L4[4][7][0][2] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][7][0][3] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][7][1][2] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][7][1][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][7][2][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][7][2][1] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][7][3][0] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][7][3][1] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[4][7][4][6] = -dImL5;
  Curvature_Higgs_CT_L4[4][7][4][7] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[4][7][5][6] = dL5;
  Curvature_Higgs_CT_L4[4][7][5][7] = -dImL5;
  Curvature_Higgs_CT_L4[4][7][6][4] = -dImL5;
  Curvature_Higgs_CT_L4[4][7][6][5] = dL5;
  Curvature_Higgs_CT_L4[4][7][7][4] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[4][7][7][5] = -dImL5;
  Curvature_Higgs_CT_L4[4][8][4][8] = dL7;
  Curvature_Higgs_CT_L4[4][8][8][4] = dL7;
  Curvature_Higgs_CT_L4[5][0][0][5] = dL1;
  Curvature_Higgs_CT_L4[5][0][2][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][0][2][7] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][0][3][6] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][0][3][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][0][5][0] = dL1;
  Curvature_Higgs_CT_L4[5][0][6][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][0][6][3] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][0][7][2] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][0][7][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][1][1][5] = dL1;
  Curvature_Higgs_CT_L4[5][1][2][6] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][1][2][7] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][1][3][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][1][3][7] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][1][5][1] = dL1;
  Curvature_Higgs_CT_L4[5][1][6][2] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][1][6][3] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][1][7][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][1][7][3] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][2][0][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][2][0][7] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][2][1][6] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][2][1][7] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][2][2][5] = dL3;
  Curvature_Higgs_CT_L4[5][2][5][2] = dL3;
  Curvature_Higgs_CT_L4[5][2][6][0] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][2][6][1] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][2][7][0] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][2][7][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][3][0][6] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][3][0][7] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][3][1][6] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][3][1][7] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][3][3][5] = dL3;
  Curvature_Higgs_CT_L4[5][3][5][3] = dL3;
  Curvature_Higgs_CT_L4[5][3][6][0] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][3][6][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][3][7][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][3][7][1] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][4][4][5] = dL1;
  Curvature_Higgs_CT_L4[5][4][5][4] = dL1;
  Curvature_Higgs_CT_L4[5][4][6][6] = dImL5;
  Curvature_Higgs_CT_L4[5][4][6][7] = dL5;
  Curvature_Higgs_CT_L4[5][4][7][6] = dL5;
  Curvature_Higgs_CT_L4[5][4][7][7] = -dImL5;
  Curvature_Higgs_CT_L4[5][5][0][0] = dL1;
  Curvature_Higgs_CT_L4[5][5][1][1] = dL1;
  Curvature_Higgs_CT_L4[5][5][2][2] = dL3;
  Curvature_Higgs_CT_L4[5][5][3][3] = dL3;
  Curvature_Higgs_CT_L4[5][5][4][4] = dL1;
  Curvature_Higgs_CT_L4[5][5][5][5] = 0.3e1 * dL1;
  Curvature_Higgs_CT_L4[5][5][6][6] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[5][5][6][7] = dImL5;
  Curvature_Higgs_CT_L4[5][5][7][6] = dImL5;
  Curvature_Higgs_CT_L4[5][5][7][7] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[5][5][8][8] = dL7;
  Curvature_Higgs_CT_L4[5][6][0][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][6][0][3] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][6][1][2] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][6][1][3] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][6][2][0] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][6][2][1] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][6][3][0] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][6][3][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][6][4][6] = dImL5;
  Curvature_Higgs_CT_L4[5][6][4][7] = dL5;
  Curvature_Higgs_CT_L4[5][6][5][6] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[5][6][5][7] = dImL5;
  Curvature_Higgs_CT_L4[5][6][6][4] = dImL5;
  Curvature_Higgs_CT_L4[5][6][6][5] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[5][6][7][4] = dL5;
  Curvature_Higgs_CT_L4[5][6][7][5] = dImL5;
  Curvature_Higgs_CT_L4[5][7][0][2] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][7][0][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][7][1][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][7][1][3] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][7][2][0] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][7][2][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][7][3][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][7][3][1] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[5][7][4][6] = dL5;
  Curvature_Higgs_CT_L4[5][7][4][7] = -dImL5;
  Curvature_Higgs_CT_L4[5][7][5][6] = dImL5;
  Curvature_Higgs_CT_L4[5][7][5][7] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[5][7][6][4] = dL5;
  Curvature_Higgs_CT_L4[5][7][6][5] = dImL5;
  Curvature_Higgs_CT_L4[5][7][7][4] = -dImL5;
  Curvature_Higgs_CT_L4[5][7][7][5] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[5][8][5][8] = dL7;
  Curvature_Higgs_CT_L4[5][8][8][5] = dL7;
  Curvature_Higgs_CT_L4[6][0][0][6] = dL3;
  Curvature_Higgs_CT_L4[6][0][2][4] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][0][2][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][0][3][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][0][3][5] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][0][4][2] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][0][4][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][0][5][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][0][5][3] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][0][6][0] = dL3;
  Curvature_Higgs_CT_L4[6][1][1][6] = dL3;
  Curvature_Higgs_CT_L4[6][1][2][4] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][1][2][5] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][1][3][4] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][1][3][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][1][4][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][1][4][3] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][1][5][2] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][1][5][3] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][1][6][1] = dL3;
  Curvature_Higgs_CT_L4[6][2][0][4] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][2][0][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][2][1][4] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][2][1][5] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][2][2][6] = dL2;
  Curvature_Higgs_CT_L4[6][2][4][0] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][2][4][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][2][5][0] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][2][5][1] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][2][6][2] = dL2;
  Curvature_Higgs_CT_L4[6][3][0][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][3][0][5] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][3][1][4] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][3][1][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][3][3][6] = dL2;
  Curvature_Higgs_CT_L4[6][3][4][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][3][4][1] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][3][5][0] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][3][5][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][3][6][3] = dL2;
  Curvature_Higgs_CT_L4[6][4][0][2] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][4][0][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][4][1][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][4][1][3] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][4][2][0] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][4][2][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][4][3][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][4][3][1] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][4][4][6] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[6][4][4][7] = -dImL5;
  Curvature_Higgs_CT_L4[6][4][5][6] = dImL5;
  Curvature_Higgs_CT_L4[6][4][5][7] = dL5;
  Curvature_Higgs_CT_L4[6][4][6][4] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[6][4][6][5] = dImL5;
  Curvature_Higgs_CT_L4[6][4][7][4] = -dImL5;
  Curvature_Higgs_CT_L4[6][4][7][5] = dL5;
  Curvature_Higgs_CT_L4[6][5][0][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][5][0][3] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][5][1][2] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][5][1][3] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][5][2][0] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][5][2][1] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][5][3][0] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][5][3][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[6][5][4][6] = dImL5;
  Curvature_Higgs_CT_L4[6][5][4][7] = dL5;
  Curvature_Higgs_CT_L4[6][5][5][6] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[6][5][5][7] = dImL5;
  Curvature_Higgs_CT_L4[6][5][6][4] = dImL5;
  Curvature_Higgs_CT_L4[6][5][6][5] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[6][5][7][4] = dL5;
  Curvature_Higgs_CT_L4[6][5][7][5] = dImL5;
  Curvature_Higgs_CT_L4[6][6][0][0] = dL3;
  Curvature_Higgs_CT_L4[6][6][1][1] = dL3;
  Curvature_Higgs_CT_L4[6][6][2][2] = dL2;
  Curvature_Higgs_CT_L4[6][6][3][3] = dL2;
  Curvature_Higgs_CT_L4[6][6][4][4] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[6][6][4][5] = dImL5;
  Curvature_Higgs_CT_L4[6][6][5][4] = dImL5;
  Curvature_Higgs_CT_L4[6][6][5][5] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[6][6][6][6] = 0.3e1 * dL2;
  Curvature_Higgs_CT_L4[6][6][7][7] = dL2;
  Curvature_Higgs_CT_L4[6][6][8][8] = dL8;
  Curvature_Higgs_CT_L4[6][7][4][4] = -dImL5;
  Curvature_Higgs_CT_L4[6][7][4][5] = dL5;
  Curvature_Higgs_CT_L4[6][7][5][4] = dL5;
  Curvature_Higgs_CT_L4[6][7][5][5] = dImL5;
  Curvature_Higgs_CT_L4[6][7][6][7] = dL2;
  Curvature_Higgs_CT_L4[6][7][7][6] = dL2;
  Curvature_Higgs_CT_L4[6][8][6][8] = dL8;
  Curvature_Higgs_CT_L4[6][8][8][6] = dL8;
  Curvature_Higgs_CT_L4[7][0][0][7] = dL3;
  Curvature_Higgs_CT_L4[7][0][2][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][0][2][5] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][0][3][4] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][0][3][5] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][0][4][2] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][0][4][3] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][0][5][2] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][0][5][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][0][7][0] = dL3;
  Curvature_Higgs_CT_L4[7][1][1][7] = dL3;
  Curvature_Higgs_CT_L4[7][1][2][4] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][1][2][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][1][3][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][1][3][5] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][1][4][2] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][1][4][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][1][5][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][1][5][3] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][1][7][1] = dL3;
  Curvature_Higgs_CT_L4[7][2][0][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][2][0][5] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][2][1][4] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][2][1][5] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][2][2][7] = dL2;
  Curvature_Higgs_CT_L4[7][2][4][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][2][4][1] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][2][5][0] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][2][5][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][2][7][2] = dL2;
  Curvature_Higgs_CT_L4[7][3][0][4] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][3][0][5] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][3][1][4] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][3][1][5] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][3][3][7] = dL2;
  Curvature_Higgs_CT_L4[7][3][4][0] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][3][4][1] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][3][5][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][3][5][1] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][3][7][3] = dL2;
  Curvature_Higgs_CT_L4[7][4][0][2] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][4][0][3] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][4][1][2] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][4][1][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][4][2][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][4][2][1] = -dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][4][3][0] = dL4 / 0.2e1 - dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][4][3][1] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][4][4][6] = -dImL5;
  Curvature_Higgs_CT_L4[7][4][4][7] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[7][4][5][6] = dL5;
  Curvature_Higgs_CT_L4[7][4][5][7] = -dImL5;
  Curvature_Higgs_CT_L4[7][4][6][4] = -dImL5;
  Curvature_Higgs_CT_L4[7][4][6][5] = dL5;
  Curvature_Higgs_CT_L4[7][4][7][4] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[7][4][7][5] = -dImL5;
  Curvature_Higgs_CT_L4[7][5][0][2] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][5][0][3] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][5][1][2] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][5][1][3] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][5][2][0] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][5][2][1] = dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][5][3][0] = -dImL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][5][3][1] = dL4 / 0.2e1 + dL5 / 0.2e1;
  Curvature_Higgs_CT_L4[7][5][4][6] = dL5;
  Curvature_Higgs_CT_L4[7][5][4][7] = -dImL5;
  Curvature_Higgs_CT_L4[7][5][5][6] = dImL5;
  Curvature_Higgs_CT_L4[7][5][5][7] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[7][5][6][4] = dL5;
  Curvature_Higgs_CT_L4[7][5][6][5] = dImL5;
  Curvature_Higgs_CT_L4[7][5][7][4] = -dImL5;
  Curvature_Higgs_CT_L4[7][5][7][5] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[7][6][4][4] = -dImL5;
  Curvature_Higgs_CT_L4[7][6][4][5] = dL5;
  Curvature_Higgs_CT_L4[7][6][5][4] = dL5;
  Curvature_Higgs_CT_L4[7][6][5][5] = dImL5;
  Curvature_Higgs_CT_L4[7][6][6][7] = dL2;
  Curvature_Higgs_CT_L4[7][6][7][6] = dL2;
  Curvature_Higgs_CT_L4[7][7][0][0] = dL3;
  Curvature_Higgs_CT_L4[7][7][1][1] = dL3;
  Curvature_Higgs_CT_L4[7][7][2][2] = dL2;
  Curvature_Higgs_CT_L4[7][7][3][3] = dL2;
  Curvature_Higgs_CT_L4[7][7][4][4] = dL3 + dL4 - dL5;
  Curvature_Higgs_CT_L4[7][7][4][5] = -dImL5;
  Curvature_Higgs_CT_L4[7][7][5][4] = -dImL5;
  Curvature_Higgs_CT_L4[7][7][5][5] = dL3 + dL4 + dL5;
  Curvature_Higgs_CT_L4[7][7][6][6] = dL2;
  Curvature_Higgs_CT_L4[7][7][7][7] = 0.3e1 * dL2;
  Curvature_Higgs_CT_L4[7][7][8][8] = dL8;
  Curvature_Higgs_CT_L4[7][8][7][8] = dL8;
  Curvature_Higgs_CT_L4[7][8][8][7] = dL8;
  Curvature_Higgs_CT_L4[8][0][0][8] = dL7;
  Curvature_Higgs_CT_L4[8][0][8][0] = dL7;
  Curvature_Higgs_CT_L4[8][1][1][8] = dL7;
  Curvature_Higgs_CT_L4[8][1][8][1] = dL7;
  Curvature_Higgs_CT_L4[8][2][2][8] = dL8;
  Curvature_Higgs_CT_L4[8][2][8][2] = dL8;
  Curvature_Higgs_CT_L4[8][3][3][8] = dL8;
  Curvature_Higgs_CT_L4[8][3][8][3] = dL8;
  Curvature_Higgs_CT_L4[8][4][4][8] = dL7;
  Curvature_Higgs_CT_L4[8][4][8][4] = dL7;
  Curvature_Higgs_CT_L4[8][5][5][8] = dL7;
  Curvature_Higgs_CT_L4[8][5][8][5] = dL7;
  Curvature_Higgs_CT_L4[8][6][6][8] = dL8;
  Curvature_Higgs_CT_L4[8][6][8][6] = dL8;
  Curvature_Higgs_CT_L4[8][7][7][8] = dL8;
  Curvature_Higgs_CT_L4[8][7][8][7] = dL8;
  Curvature_Higgs_CT_L4[8][8][0][0] = dL7;
  Curvature_Higgs_CT_L4[8][8][1][1] = dL7;
  Curvature_Higgs_CT_L4[8][8][2][2] = dL8;
  Curvature_Higgs_CT_L4[8][8][3][3] = dL8;
  Curvature_Higgs_CT_L4[8][8][4][4] = dL7;
  Curvature_Higgs_CT_L4[8][8][5][5] = dL7;
  Curvature_Higgs_CT_L4[8][8][6][6] = dL8;
  Curvature_Higgs_CT_L4[8][8][7][7] = dL8;
  Curvature_Higgs_CT_L4[8][8][8][8] = 0.6e1 * dL6;
}

/**
 * console output of all parameters
 */
void Class_Potential_CPintheDark::write() const
{
  std::stringstream ss;
  typedef std::numeric_limits<double> dbl;
  ss.precision(dbl::max_digits10);

  ss << "Model = " << Model << "\n";

  ss << "The parameters are : "
     << "\n";
  ss << "m11s = " << m11s << "\n";
  ss << "m22s = " << m22s << "\n";
  ss << "mSs = " << mSs << "\n";
  ss << "ReA = " << ReA << "\n";
  ss << "ImA = " << ImA << "\n";
  ss << "L1 = " << L1 << "\n";
  ss << "L2 = " << L2 << "\n";
  ss << "L3 = " << L3 << "\n";
  ss << "L4 = " << L4 << "\n";
  ss << "L5 = " << L5 << "\n";
  ss << "L6 = " << L6 << "\n";
  ss << "L7 = " << L7 << "\n";
  ss << "L8 = " << L8 << "\n";

  ss << "The counterterm parameters are : "
     << "\n";
  ss << "dm11s = " << dm11s << "\n";
  ss << "dm22s = " << dm22s << "\n";
  ss << "dmSs = " << dmSs << "\n";
  ss << "dReA = " << dReA << "\n";
  ss << "dImA = " << dImA << "\n";
  ss << "dL1 = " << dL1 << "\n";
  ss << "dL2 = " << dL2 << "\n";
  ss << "dL3 = " << dL3 << "\n";
  ss << "dL4 = " << dL4 << "\n";
  ss << "dL5 = " << dL5 << "\n";
  ss << "dL6 = " << dL6 << "\n";
  ss << "dL7 = " << dL7 << "\n";
  ss << "dL8 = " << dL8 << "\n";
  ss << "dTCB = " << dTCB << "\n";
  ss << "dT1 = " << dT1 << "\n";
  ss << "dT2 = " << dT2 << "\n";
  ss << "dTCP = " << dTCP << "\n";
  ss << "dTS = " << dTS << "\n";
  ss << "dImL5 = " << dImL5 << "\n";

  ss << "The scale is given by mu = " << scale << " GeV \n";

  std::vector<double> HiggsMasses;
  HiggsMasses = HiggsMassesSquared(vevTree, 0);

  ss << "The mass spectrum is given by :\n";
  ss << "m_{G^+} = " << std::sqrt(HiggsMasses[pos_Gp]) << " GeV \n"
     << "m_{G^-} = " << std::sqrt(HiggsMasses[pos_Gm]) << " GeV \n"
     << "m_{H^+} = " << std::sqrt(HiggsMasses[pos_Hp]) << " GeV \n"
     << "m_{H^-} = " << std::sqrt(HiggsMasses[pos_Hm]) << " GeV \n"
     << "m_{hSM} = " << std::sqrt(HiggsMasses[pos_HSM]) << " GeV \n"
     << "m_{G^0} = " << std::sqrt(HiggsMasses[pos_G0]) << " GeV \n"
     << "m_{h_1} = " << std::sqrt(HiggsMasses[pos_h1]) << " GeV \n"
     << "m_{h_2} = " << std::sqrt(HiggsMasses[pos_h2]) << " GeV \n"
     << "m_{h_3} = " << std::sqrt(HiggsMasses[pos_h3]) << " GeV \n";

  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms. Here you need to work out the scheme and
 * implement the formulas.
 */
std::vector<double> Class_Potential_CPintheDark::calc_CT() const
{
  std::vector<double> parCT;

  if (!SetCurvatureDone)
  {
    std::string retmes = __func__;
    retmes += " was called before SetCurvatureArrays()!\n";
    throw std::runtime_error(retmes);
  }
  if (!CalcCouplingsDone)
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

  // set the free parameters
  double L2ct = 0;
  double L3ct = 0;
  double L6ct = 0;
  double L7ct = 0;
  double L8ct = 0;

  // formulae for the counterterm scheme
  parCT.push_back(HesseWeinberg(4, 4) / 0.2e1 -
                  0.3e1 / 0.2e1 * HesseWeinberg(5, 5));              // dm11s
  parCT.push_back(-HesseWeinberg(3, 3) - pow(v1, 2) / 0.2e1 * L3ct); // dm22s
  parCT.push_back(-HesseWeinberg(8, 8) - pow(v1, 2) / 0.2e1 * L7ct); // dmSs
  parCT.push_back(-HesseWeinberg(6, 8) / v1);                        // dReA
  parCT.push_back(HesseWeinberg(7, 8) / v1);                         // dImA
  parCT.push_back((-HesseWeinberg(4, 4) + HesseWeinberg(5, 5)) *
                  pow(v1, -2)); // dL1
  parCT.push_back(L2ct);        // dL2
  parCT.push_back(L3ct);        // dL3
  parCT.push_back((0.2e1 * HesseWeinberg(3, 3) - HesseWeinberg(6, 6) -
                   HesseWeinberg(7, 7)) *
                  pow(v1, -2)); // dL4
  parCT.push_back((-HesseWeinberg(6, 6) + HesseWeinberg(7, 7)) *
                  pow(v1, -2));                                 // dL5
  parCT.push_back(L6ct);                                        // dL6
  parCT.push_back(L7ct);                                        // dL7
  parCT.push_back(L8ct);                                        // dL8
  parCT.push_back(-NablaWeinberg(2));                           // dTCB
  parCT.push_back(HesseWeinberg(5, 5) * v1 - NablaWeinberg(4)); // dT1
  parCT.push_back(-NablaWeinberg(6));                           // dT2
  parCT.push_back(-NablaWeinberg(7));                           // dTCP
  parCT.push_back(-NablaWeinberg(8));                           // dTS
  parCT.push_back(0.2e1 * HesseWeinberg(6, 7) * pow(v1, -2));   // dImL5

  return parCT;
}

/**
 * Ensures the correct rotation matrix convention
 */
void Class_Potential_CPintheDark::AdjustRotationMatrix()
{
  const double ZeroThreshold = 1e-5;

  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsDone) CalculatePhysicalCouplings();

  if (!CheckRotationMatrix()) // Check whether generically generated rotation
                              // matrix is proper rotation matrix
  {
    throw std::runtime_error("Error in rotation matrix.");
  }

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
    }
  }

  std::vector<double> HiggsMasses = HiggsMassesSquared(vevTree, 0);
  if (HiggsMasses.front() <= -ZeroThreshold)
  {
    std::stringstream ss;
    ss.precision(std::numeric_limits<double>::max_digits10);
    ss << "Warning, at least one negative mass squared in spectrum: "
       << HiggsMasses.front() << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
  }

  // CP in the Dark interaction basis
  // 0     1     2     3     4      5     6      7     8
  // rho1, eta1, rho2, eta2, zeta1, psi1, zeta2, psi2, rhoS
  const std::size_t pos_rho1 = 0, pos_eta1 = 1, pos_rho2 = 2, pos_eta2 = 3,
                    pos_zeta1 = 4, pos_psi1 = 5, pos_zeta2 = 6, pos_psi2 = 7,
                    pos_rhoS = 8;

  // Indices of mass eigenstates for rotation from interaction to mass basis
  std::optional<std::size_t> tpos_Gp, tpos_Gm, tpos_Hp, tpos_Hm, tpos_HSM,
                             tpos_G0, tpos_h1, tpos_h2, tpos_h3;

  // basis = {rho1, eta1, rho2, eta2, zeta1, psi1, zeta2, psi2, rhoS}
  // the rotation matrix is diagonal besides for the neutral dark scalars
  for (std::size_t i = 0; i < NHiggs;
       i++) // mass base index i corresponds to mass vector sorted in ascending
            // mass
  {
    bool hasZeroMass = std::abs(HiggsMasses[i]) < ZeroThreshold;
    if (std::abs(HiggsRot(i, pos_rho1)) > ZeroThreshold)
    {
      if (not tpos_Gp.has_value() and hasZeroMass)
      {
        tpos_Gp = i;
      }
      else
      {
        throw std::runtime_error("Error. Goldstone Gp not massless "
                                 "or not diagonal.");
      }
    }
    else if (std::abs(HiggsRot(i, pos_eta1)) > ZeroThreshold)
    {
      if (not tpos_Gm.has_value() and hasZeroMass)
      {
        tpos_Gm = i;
      }
      else
      {
        throw std::runtime_error("Error. Goldstone Gm not massless "
                                 "or not diagonal.");
      }
    }
    else if (std::abs(HiggsRot(i, pos_rho2)) > ZeroThreshold)
    {
      if (not tpos_Hp.has_value())
      {
        tpos_Hp = i;
      }
      else
      {
        throw std::runtime_error("Error. Charged Higgs Hp not diagonal.");
      }
    }
    else if (std::abs(HiggsRot(i, pos_eta2)) > ZeroThreshold)
    {
      if (not tpos_Hm.has_value())
      {
        tpos_Hm = i;
      }
      else
      {
        throw std::runtime_error("Error. Charged Higgs Hm not diagonal.");
      }
    }
    else if (std::abs(HiggsRot(i, pos_zeta1)) > ZeroThreshold)
    {
      if (not tpos_HSM.has_value())
      {
        tpos_HSM = i;
      }
      else
      {
        throw std::runtime_error("Error. Higgs HSM not diagonal.");
      }
    }
    else if (std::abs(HiggsRot(i, pos_psi1)) > ZeroThreshold)
    {
      if (not tpos_G0.has_value() and hasZeroMass)
      {
        tpos_G0 = i;
      }
      else
      {
        throw std::runtime_error("Error. Goldstone G0 not massless "
                                 "or not diagonal.");
      }
    }
    // the neutral dark scalars mix
    else if ((std::abs(HiggsRot(i, pos_zeta2)) +
              std::abs(HiggsRot(i, pos_psi2)) +
              std::abs(HiggsRot(i, pos_rhoS))) > ZeroThreshold)
    {
      // use that scalars are sorted by mass
      if (not tpos_h1.has_value())
      {
        tpos_h1 = i;
      }
      else if (not tpos_h2.has_value())
      {
        tpos_h2 = i;
      }
      else if (not tpos_h3.has_value())
      {
        tpos_h3 = i;
      }
      else
      {
        throw std::runtime_error("Error. Neutral submatrix mixing "
                                 "with other components.");
      }
    }
    else
    {
      throw std::runtime_error("Error. Invalid mixing matrix containing row "
                               "with all zeroes.");
    }
  }

  // Sanity check if all position indices are set
  if (not (tpos_Gp.has_value() and tpos_Gm.has_value() and tpos_Hp.has_value()
           and tpos_Hm.has_value() and tpos_HSM.has_value()
           and tpos_G0.has_value() and tpos_h1.has_value()
           and tpos_h2.has_value() and tpos_h3.has_value())
     )
  {
    throw std::runtime_error("Error. Not all position indices are set.");
  }

  pos_Gp = tpos_Gp.value();
  pos_Gm = tpos_Gm.value();
  pos_Hp = tpos_Hp.value();
  pos_Hm = tpos_Hm.value();
  pos_HSM = tpos_HSM.value();
  pos_G0 = tpos_G0.value();
  pos_h1 = tpos_h1.value();
  pos_h2 = tpos_h2.value();
  pos_h3 = tpos_h3.value();

  // check if all other elements of rotation matrix are zero
  bool zero_element = false;
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      if (not((j == pos_rho1 and i == pos_Gp) or
              (j == pos_eta1 and i == pos_Gm) or
              (j == pos_rho2 and i == pos_Hp) or
              (j == pos_eta2 and i == pos_Hm) or
              (j == pos_zeta1 and i == pos_HSM) or
              (j == pos_psi1 and i == pos_G0) or
              (j == pos_zeta2 and
               (i == pos_h1 or i == pos_h2 or i == pos_h3)) or
              (j == pos_psi2 and
               (i == pos_h1 or i == pos_h2 or i == pos_h3)) or
              (j == pos_rhoS and
               (i == pos_h1 or i == pos_h2 or i == pos_h3))))
      {
        zero_element = true;
      }

      if (zero_element and std::abs(HiggsRot(i, j)) > ZeroThreshold)
      {
        throw std::runtime_error("Error. Invalid rotation matrix detected.");
      }
      zero_element = false;
    }
  }

  MatrixXd HiggsRotFixed(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsRotFixed.row(i) = HiggsRot.row(i);
  }

  if (HiggsRotFixed(pos_Gp, pos_rho1) < 0) // Gp rho1 (= +1)
  {
    HiggsRotFixed.row(pos_Gp) *= -1;
  }
  if (HiggsRotFixed(pos_Gm, pos_eta1) < 0) // Gm eta1 (= +1)
  {
    HiggsRotFixed.row(pos_Gm) *= -1;
  }
  if (HiggsRotFixed(pos_Hp, pos_rho2) < 0) // Hp rho2 (= +1)
  {
    HiggsRotFixed.row(pos_Hp) *= -1;
  }
  if (HiggsRotFixed(pos_Hm, pos_eta2) < 0) // Hm eta2 (= +1)
  {
    HiggsRotFixed.row(pos_Hm) *= -1;
  }
  if (HiggsRotFixed(pos_HSM, pos_zeta1) < 0) // HSM zeta1 (= +1)
  {
    HiggsRotFixed.row(pos_HSM) *= -1;
  }
  if (HiggsRotFixed(pos_G0, pos_psi1) < 0) // G0 psi1 (= +1)
  {
    HiggsRotFixed.row(pos_G0) *= -1;
  }

  // Check dark neutral submatrix
  // Use the "ScannerS" criteria from arXiv:2007.02985 Eq. (6)
  // (since ScannerS uses the same parametrisation of the angles as BSMPT):
  // * (1) if R[1][1] < 0: h1 -> -h1 (i.e. multiply the h1 row with -1)
  // * (2) if R[3][3] < 0: h3 -> -h3 (i.e. multiply the h3 row with -1)
  // * (3) if det R < 0: h2 -> -h2 (i.e. multiply the h2 row with -1)

  if (HiggsRotFixed(pos_h1, pos_zeta2) < 0)
  // h1 zeta2 (condition (1) above, R11 < 0)
  {
    // if negative, flip sign of h1
    HiggsRotFixed.row(pos_h1) *= -1;
  }

  if (HiggsRotFixed(pos_h3, pos_rhoS) < 0)
  // h3 rhoS (condition (2) above, R33 < 0)
  {
    // if negative, flip sign of h3
    HiggsRotFixed.row(pos_h3) *= -1;
  }

  // Calculate the determinant AFTER flipping the signs for rows 1 and 3 above
  MatrixXd HiggsRotFixedNeutral(3, 3);
  HiggsRotFixedNeutral(0, 0) = HiggsRotFixed(pos_h1, pos_zeta2);
  HiggsRotFixedNeutral(0, 1) = HiggsRotFixed(pos_h1, pos_psi2);
  HiggsRotFixedNeutral(0, 2) = HiggsRotFixed(pos_h1, pos_rhoS);

  HiggsRotFixedNeutral(1, 0) = HiggsRotFixed(pos_h2, pos_zeta2);
  HiggsRotFixedNeutral(1, 1) = HiggsRotFixed(pos_h2, pos_psi2);
  HiggsRotFixedNeutral(1, 2) = HiggsRotFixed(pos_h2, pos_rhoS);

  HiggsRotFixedNeutral(2, 0) = HiggsRotFixed(pos_h3, pos_zeta2);
  HiggsRotFixedNeutral(2, 1) = HiggsRotFixed(pos_h3, pos_psi2);
  HiggsRotFixedNeutral(2, 2) = HiggsRotFixed(pos_h3, pos_rhoS);

  if (HiggsRotFixedNeutral.determinant() < 0)
  // condition (3) above, det(R) < 0
  {
    // if negative, flip sign of h2
    HiggsRotFixed.row(pos_h2) *= -1;
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRotationMatrixEnsuredConvention[i][j] = HiggsRotFixed(i, j);
    }
  }

  return;
}

// mass basis triple couplings
void Class_Potential_CPintheDark::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsDone) CalculatePhysicalCouplings();

  if (CalculatedTripleCopulings) return;
  CalculatedTripleCopulings = true;

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrixEnsuredConvention[i][j];
    }
  }

  std::vector<std::size_t> HiggsOrder(NHiggs);
  HiggsOrder[0] = pos_Gp;
  HiggsOrder[1] = pos_Gm;
  HiggsOrder[2] = pos_Hp;
  HiggsOrder[3] = pos_Hm;
  HiggsOrder[4] = pos_HSM;
  HiggsOrder[5] = pos_G0;
  HiggsOrder[6] = pos_h1;
  HiggsOrder[7] = pos_h2;
  HiggsOrder[8] = pos_h3;

  MatrixXd HiggsRotSort(NHiggs, NHiggs);
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

void Class_Potential_CPintheDark::SetCurvatureArrays()
{
  initVectors();
  SetCurvatureDone = true;
  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

  // assign the non-zero entries
  Curvature_Higgs_L2[0][0] = m11s;
  Curvature_Higgs_L2[1][1] = m11s;
  Curvature_Higgs_L2[2][2] = m22s;
  Curvature_Higgs_L2[3][3] = m22s;
  Curvature_Higgs_L2[4][4] = m11s;
  Curvature_Higgs_L2[5][5] = m11s;
  Curvature_Higgs_L2[6][6] = m22s;
  Curvature_Higgs_L2[7][7] = m22s;
  Curvature_Higgs_L2[8][8] = mSs;

  Curvature_Higgs_L3[0][2][8] = ReA;
  Curvature_Higgs_L3[0][3][8] = -ImA;
  Curvature_Higgs_L3[0][8][2] = ReA;
  Curvature_Higgs_L3[0][8][3] = -ImA;
  Curvature_Higgs_L3[1][2][8] = ImA;
  Curvature_Higgs_L3[1][3][8] = ReA;
  Curvature_Higgs_L3[1][8][2] = ImA;
  Curvature_Higgs_L3[1][8][3] = ReA;
  Curvature_Higgs_L3[2][0][8] = ReA;
  Curvature_Higgs_L3[2][1][8] = ImA;
  Curvature_Higgs_L3[2][8][0] = ReA;
  Curvature_Higgs_L3[2][8][1] = ImA;
  Curvature_Higgs_L3[3][0][8] = -ImA;
  Curvature_Higgs_L3[3][1][8] = ReA;
  Curvature_Higgs_L3[3][8][0] = -ImA;
  Curvature_Higgs_L3[3][8][1] = ReA;
  Curvature_Higgs_L3[4][6][8] = ReA;
  Curvature_Higgs_L3[4][7][8] = -ImA;
  Curvature_Higgs_L3[4][8][6] = ReA;
  Curvature_Higgs_L3[4][8][7] = -ImA;
  Curvature_Higgs_L3[5][6][8] = ImA;
  Curvature_Higgs_L3[5][7][8] = ReA;
  Curvature_Higgs_L3[5][8][6] = ImA;
  Curvature_Higgs_L3[5][8][7] = ReA;
  Curvature_Higgs_L3[6][4][8] = ReA;
  Curvature_Higgs_L3[6][5][8] = ImA;
  Curvature_Higgs_L3[6][8][4] = ReA;
  Curvature_Higgs_L3[6][8][5] = ImA;
  Curvature_Higgs_L3[7][4][8] = -ImA;
  Curvature_Higgs_L3[7][5][8] = ReA;
  Curvature_Higgs_L3[7][8][4] = -ImA;
  Curvature_Higgs_L3[7][8][5] = ReA;
  Curvature_Higgs_L3[8][0][2] = ReA;
  Curvature_Higgs_L3[8][0][3] = -ImA;
  Curvature_Higgs_L3[8][1][2] = ImA;
  Curvature_Higgs_L3[8][1][3] = ReA;
  Curvature_Higgs_L3[8][2][0] = ReA;
  Curvature_Higgs_L3[8][2][1] = ImA;
  Curvature_Higgs_L3[8][3][0] = -ImA;
  Curvature_Higgs_L3[8][3][1] = ReA;
  Curvature_Higgs_L3[8][4][6] = ReA;
  Curvature_Higgs_L3[8][4][7] = -ImA;
  Curvature_Higgs_L3[8][5][6] = ImA;
  Curvature_Higgs_L3[8][5][7] = ReA;
  Curvature_Higgs_L3[8][6][4] = ReA;
  Curvature_Higgs_L3[8][6][5] = ImA;
  Curvature_Higgs_L3[8][7][4] = -ImA;
  Curvature_Higgs_L3[8][7][5] = ReA;

  Curvature_Higgs_L4[0][0][0][0] = 0.3e1 * L1;
  Curvature_Higgs_L4[0][0][1][1] = L1;
  Curvature_Higgs_L4[0][0][2][2] = L3 + L4 + L5;
  Curvature_Higgs_L4[0][0][3][3] = L3 + L4 - L5;
  Curvature_Higgs_L4[0][0][4][4] = L1;
  Curvature_Higgs_L4[0][0][5][5] = L1;
  Curvature_Higgs_L4[0][0][6][6] = L3;
  Curvature_Higgs_L4[0][0][7][7] = L3;
  Curvature_Higgs_L4[0][0][8][8] = L7;
  Curvature_Higgs_L4[0][1][0][1] = L1;
  Curvature_Higgs_L4[0][1][1][0] = L1;
  Curvature_Higgs_L4[0][1][2][3] = L5;
  Curvature_Higgs_L4[0][1][3][2] = L5;
  Curvature_Higgs_L4[0][2][0][2] = L3 + L4 + L5;
  Curvature_Higgs_L4[0][2][1][3] = L5;
  Curvature_Higgs_L4[0][2][2][0] = L3 + L4 + L5;
  Curvature_Higgs_L4[0][2][3][1] = L5;
  Curvature_Higgs_L4[0][2][4][6] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][2][5][7] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][2][6][4] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][2][7][5] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][3][0][3] = L3 + L4 - L5;
  Curvature_Higgs_L4[0][3][1][2] = L5;
  Curvature_Higgs_L4[0][3][2][1] = L5;
  Curvature_Higgs_L4[0][3][3][0] = L3 + L4 - L5;
  Curvature_Higgs_L4[0][3][4][7] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[0][3][5][6] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][3][6][5] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][3][7][4] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[0][4][0][4] = L1;
  Curvature_Higgs_L4[0][4][2][6] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][4][3][7] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[0][4][4][0] = L1;
  Curvature_Higgs_L4[0][4][6][2] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][4][7][3] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[0][5][0][5] = L1;
  Curvature_Higgs_L4[0][5][2][7] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][5][3][6] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][5][5][0] = L1;
  Curvature_Higgs_L4[0][5][6][3] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][5][7][2] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][6][0][6] = L3;
  Curvature_Higgs_L4[0][6][2][4] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][6][3][5] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][6][4][2] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][6][5][3] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][6][6][0] = L3;
  Curvature_Higgs_L4[0][7][0][7] = L3;
  Curvature_Higgs_L4[0][7][2][5] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][7][3][4] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[0][7][4][3] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[0][7][5][2] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[0][7][7][0] = L3;
  Curvature_Higgs_L4[0][8][0][8] = L7;
  Curvature_Higgs_L4[0][8][8][0] = L7;
  Curvature_Higgs_L4[1][0][0][1] = L1;
  Curvature_Higgs_L4[1][0][1][0] = L1;
  Curvature_Higgs_L4[1][0][2][3] = L5;
  Curvature_Higgs_L4[1][0][3][2] = L5;
  Curvature_Higgs_L4[1][1][0][0] = L1;
  Curvature_Higgs_L4[1][1][1][1] = 0.3e1 * L1;
  Curvature_Higgs_L4[1][1][2][2] = L3 + L4 - L5;
  Curvature_Higgs_L4[1][1][3][3] = L3 + L4 + L5;
  Curvature_Higgs_L4[1][1][4][4] = L1;
  Curvature_Higgs_L4[1][1][5][5] = L1;
  Curvature_Higgs_L4[1][1][6][6] = L3;
  Curvature_Higgs_L4[1][1][7][7] = L3;
  Curvature_Higgs_L4[1][1][8][8] = L7;
  Curvature_Higgs_L4[1][2][0][3] = L5;
  Curvature_Higgs_L4[1][2][1][2] = L3 + L4 - L5;
  Curvature_Higgs_L4[1][2][2][1] = L3 + L4 - L5;
  Curvature_Higgs_L4[1][2][3][0] = L5;
  Curvature_Higgs_L4[1][2][4][7] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][2][5][6] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[1][2][6][5] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[1][2][7][4] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][3][0][2] = L5;
  Curvature_Higgs_L4[1][3][1][3] = L3 + L4 + L5;
  Curvature_Higgs_L4[1][3][2][0] = L5;
  Curvature_Higgs_L4[1][3][3][1] = L3 + L4 + L5;
  Curvature_Higgs_L4[1][3][4][6] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][3][5][7] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][3][6][4] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][3][7][5] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][4][1][4] = L1;
  Curvature_Higgs_L4[1][4][2][7] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][4][3][6] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][4][4][1] = L1;
  Curvature_Higgs_L4[1][4][6][3] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][4][7][2] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][5][1][5] = L1;
  Curvature_Higgs_L4[1][5][2][6] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[1][5][3][7] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][5][5][1] = L1;
  Curvature_Higgs_L4[1][5][6][2] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[1][5][7][3] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][6][1][6] = L3;
  Curvature_Higgs_L4[1][6][2][5] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[1][6][3][4] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][6][4][3] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][6][5][2] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[1][6][6][1] = L3;
  Curvature_Higgs_L4[1][7][1][7] = L3;
  Curvature_Higgs_L4[1][7][2][4] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][7][3][5] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][7][4][2] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][7][5][3] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[1][7][7][1] = L3;
  Curvature_Higgs_L4[1][8][1][8] = L7;
  Curvature_Higgs_L4[1][8][8][1] = L7;
  Curvature_Higgs_L4[2][0][0][2] = L3 + L4 + L5;
  Curvature_Higgs_L4[2][0][1][3] = L5;
  Curvature_Higgs_L4[2][0][2][0] = L3 + L4 + L5;
  Curvature_Higgs_L4[2][0][3][1] = L5;
  Curvature_Higgs_L4[2][0][4][6] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][0][5][7] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][0][6][4] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][0][7][5] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][1][0][3] = L5;
  Curvature_Higgs_L4[2][1][1][2] = L3 + L4 - L5;
  Curvature_Higgs_L4[2][1][2][1] = L3 + L4 - L5;
  Curvature_Higgs_L4[2][1][3][0] = L5;
  Curvature_Higgs_L4[2][1][4][7] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][1][5][6] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[2][1][6][5] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[2][1][7][4] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][2][0][0] = L3 + L4 + L5;
  Curvature_Higgs_L4[2][2][1][1] = L3 + L4 - L5;
  Curvature_Higgs_L4[2][2][2][2] = 0.3e1 * L2;
  Curvature_Higgs_L4[2][2][3][3] = L2;
  Curvature_Higgs_L4[2][2][4][4] = L3;
  Curvature_Higgs_L4[2][2][5][5] = L3;
  Curvature_Higgs_L4[2][2][6][6] = L2;
  Curvature_Higgs_L4[2][2][7][7] = L2;
  Curvature_Higgs_L4[2][2][8][8] = L8;
  Curvature_Higgs_L4[2][3][0][1] = L5;
  Curvature_Higgs_L4[2][3][1][0] = L5;
  Curvature_Higgs_L4[2][3][2][3] = L2;
  Curvature_Higgs_L4[2][3][3][2] = L2;
  Curvature_Higgs_L4[2][4][0][6] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][4][1][7] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][4][2][4] = L3;
  Curvature_Higgs_L4[2][4][4][2] = L3;
  Curvature_Higgs_L4[2][4][6][0] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][4][7][1] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][5][0][7] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][5][1][6] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[2][5][2][5] = L3;
  Curvature_Higgs_L4[2][5][5][2] = L3;
  Curvature_Higgs_L4[2][5][6][1] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[2][5][7][0] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][6][0][4] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][6][1][5] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[2][6][2][6] = L2;
  Curvature_Higgs_L4[2][6][4][0] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][6][5][1] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[2][6][6][2] = L2;
  Curvature_Higgs_L4[2][7][0][5] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][7][1][4] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][7][2][7] = L2;
  Curvature_Higgs_L4[2][7][4][1] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][7][5][0] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[2][7][7][2] = L2;
  Curvature_Higgs_L4[2][8][2][8] = L8;
  Curvature_Higgs_L4[2][8][8][2] = L8;
  Curvature_Higgs_L4[3][0][0][3] = L3 + L4 - L5;
  Curvature_Higgs_L4[3][0][1][2] = L5;
  Curvature_Higgs_L4[3][0][2][1] = L5;
  Curvature_Higgs_L4[3][0][3][0] = L3 + L4 - L5;
  Curvature_Higgs_L4[3][0][4][7] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[3][0][5][6] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][0][6][5] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][0][7][4] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[3][1][0][2] = L5;
  Curvature_Higgs_L4[3][1][1][3] = L3 + L4 + L5;
  Curvature_Higgs_L4[3][1][2][0] = L5;
  Curvature_Higgs_L4[3][1][3][1] = L3 + L4 + L5;
  Curvature_Higgs_L4[3][1][4][6] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][1][5][7] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][1][6][4] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][1][7][5] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][2][0][1] = L5;
  Curvature_Higgs_L4[3][2][1][0] = L5;
  Curvature_Higgs_L4[3][2][2][3] = L2;
  Curvature_Higgs_L4[3][2][3][2] = L2;
  Curvature_Higgs_L4[3][3][0][0] = L3 + L4 - L5;
  Curvature_Higgs_L4[3][3][1][1] = L3 + L4 + L5;
  Curvature_Higgs_L4[3][3][2][2] = L2;
  Curvature_Higgs_L4[3][3][3][3] = 0.3e1 * L2;
  Curvature_Higgs_L4[3][3][4][4] = L3;
  Curvature_Higgs_L4[3][3][5][5] = L3;
  Curvature_Higgs_L4[3][3][6][6] = L2;
  Curvature_Higgs_L4[3][3][7][7] = L2;
  Curvature_Higgs_L4[3][3][8][8] = L8;
  Curvature_Higgs_L4[3][4][0][7] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[3][4][1][6] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][4][3][4] = L3;
  Curvature_Higgs_L4[3][4][4][3] = L3;
  Curvature_Higgs_L4[3][4][6][1] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][4][7][0] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[3][5][0][6] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][5][1][7] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][5][3][5] = L3;
  Curvature_Higgs_L4[3][5][5][3] = L3;
  Curvature_Higgs_L4[3][5][6][0] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][5][7][1] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][6][0][5] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][6][1][4] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][6][3][6] = L2;
  Curvature_Higgs_L4[3][6][4][1] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][6][5][0] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][6][6][3] = L2;
  Curvature_Higgs_L4[3][7][0][4] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[3][7][1][5] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][7][3][7] = L2;
  Curvature_Higgs_L4[3][7][4][0] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[3][7][5][1] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[3][7][7][3] = L2;
  Curvature_Higgs_L4[3][8][3][8] = L8;
  Curvature_Higgs_L4[3][8][8][3] = L8;
  Curvature_Higgs_L4[4][0][0][4] = L1;
  Curvature_Higgs_L4[4][0][2][6] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][0][3][7] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[4][0][4][0] = L1;
  Curvature_Higgs_L4[4][0][6][2] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][0][7][3] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[4][1][1][4] = L1;
  Curvature_Higgs_L4[4][1][2][7] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][1][3][6] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][1][4][1] = L1;
  Curvature_Higgs_L4[4][1][6][3] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][1][7][2] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][2][0][6] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][2][1][7] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][2][2][4] = L3;
  Curvature_Higgs_L4[4][2][4][2] = L3;
  Curvature_Higgs_L4[4][2][6][0] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][2][7][1] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][3][0][7] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[4][3][1][6] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][3][3][4] = L3;
  Curvature_Higgs_L4[4][3][4][3] = L3;
  Curvature_Higgs_L4[4][3][6][1] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][3][7][0] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[4][4][0][0] = L1;
  Curvature_Higgs_L4[4][4][1][1] = L1;
  Curvature_Higgs_L4[4][4][2][2] = L3;
  Curvature_Higgs_L4[4][4][3][3] = L3;
  Curvature_Higgs_L4[4][4][4][4] = 0.3e1 * L1;
  Curvature_Higgs_L4[4][4][5][5] = L1;
  Curvature_Higgs_L4[4][4][6][6] = L3 + L4 + L5;
  Curvature_Higgs_L4[4][4][7][7] = L3 + L4 - L5;
  Curvature_Higgs_L4[4][4][8][8] = L7;
  Curvature_Higgs_L4[4][5][4][5] = L1;
  Curvature_Higgs_L4[4][5][5][4] = L1;
  Curvature_Higgs_L4[4][5][6][7] = L5;
  Curvature_Higgs_L4[4][5][7][6] = L5;
  Curvature_Higgs_L4[4][6][0][2] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][6][1][3] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][6][2][0] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][6][3][1] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][6][4][6] = L3 + L4 + L5;
  Curvature_Higgs_L4[4][6][5][7] = L5;
  Curvature_Higgs_L4[4][6][6][4] = L3 + L4 + L5;
  Curvature_Higgs_L4[4][6][7][5] = L5;
  Curvature_Higgs_L4[4][7][0][3] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[4][7][1][2] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][7][2][1] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[4][7][3][0] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[4][7][4][7] = L3 + L4 - L5;
  Curvature_Higgs_L4[4][7][5][6] = L5;
  Curvature_Higgs_L4[4][7][6][5] = L5;
  Curvature_Higgs_L4[4][7][7][4] = L3 + L4 - L5;
  Curvature_Higgs_L4[4][8][4][8] = L7;
  Curvature_Higgs_L4[4][8][8][4] = L7;
  Curvature_Higgs_L4[5][0][0][5] = L1;
  Curvature_Higgs_L4[5][0][2][7] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][0][3][6] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][0][5][0] = L1;
  Curvature_Higgs_L4[5][0][6][3] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][0][7][2] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][1][1][5] = L1;
  Curvature_Higgs_L4[5][1][2][6] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[5][1][3][7] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][1][5][1] = L1;
  Curvature_Higgs_L4[5][1][6][2] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[5][1][7][3] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][2][0][7] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][2][1][6] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[5][2][2][5] = L3;
  Curvature_Higgs_L4[5][2][5][2] = L3;
  Curvature_Higgs_L4[5][2][6][1] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[5][2][7][0] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][3][0][6] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][3][1][7] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][3][3][5] = L3;
  Curvature_Higgs_L4[5][3][5][3] = L3;
  Curvature_Higgs_L4[5][3][6][0] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][3][7][1] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][4][4][5] = L1;
  Curvature_Higgs_L4[5][4][5][4] = L1;
  Curvature_Higgs_L4[5][4][6][7] = L5;
  Curvature_Higgs_L4[5][4][7][6] = L5;
  Curvature_Higgs_L4[5][5][0][0] = L1;
  Curvature_Higgs_L4[5][5][1][1] = L1;
  Curvature_Higgs_L4[5][5][2][2] = L3;
  Curvature_Higgs_L4[5][5][3][3] = L3;
  Curvature_Higgs_L4[5][5][4][4] = L1;
  Curvature_Higgs_L4[5][5][5][5] = 0.3e1 * L1;
  Curvature_Higgs_L4[5][5][6][6] = L3 + L4 - L5;
  Curvature_Higgs_L4[5][5][7][7] = L3 + L4 + L5;
  Curvature_Higgs_L4[5][5][8][8] = L7;
  Curvature_Higgs_L4[5][6][0][3] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][6][1][2] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[5][6][2][1] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[5][6][3][0] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][6][4][7] = L5;
  Curvature_Higgs_L4[5][6][5][6] = L3 + L4 - L5;
  Curvature_Higgs_L4[5][6][6][5] = L3 + L4 - L5;
  Curvature_Higgs_L4[5][6][7][4] = L5;
  Curvature_Higgs_L4[5][7][0][2] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][7][1][3] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][7][2][0] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][7][3][1] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[5][7][4][6] = L5;
  Curvature_Higgs_L4[5][7][5][7] = L3 + L4 + L5;
  Curvature_Higgs_L4[5][7][6][4] = L5;
  Curvature_Higgs_L4[5][7][7][5] = L3 + L4 + L5;
  Curvature_Higgs_L4[5][8][5][8] = L7;
  Curvature_Higgs_L4[5][8][8][5] = L7;
  Curvature_Higgs_L4[6][0][0][6] = L3;
  Curvature_Higgs_L4[6][0][2][4] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][0][3][5] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][0][4][2] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][0][5][3] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][0][6][0] = L3;
  Curvature_Higgs_L4[6][1][1][6] = L3;
  Curvature_Higgs_L4[6][1][2][5] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[6][1][3][4] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][1][4][3] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][1][5][2] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[6][1][6][1] = L3;
  Curvature_Higgs_L4[6][2][0][4] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][2][1][5] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[6][2][2][6] = L2;
  Curvature_Higgs_L4[6][2][4][0] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][2][5][1] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[6][2][6][2] = L2;
  Curvature_Higgs_L4[6][3][0][5] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][3][1][4] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][3][3][6] = L2;
  Curvature_Higgs_L4[6][3][4][1] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][3][5][0] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][3][6][3] = L2;
  Curvature_Higgs_L4[6][4][0][2] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][4][1][3] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][4][2][0] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][4][3][1] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][4][4][6] = L3 + L4 + L5;
  Curvature_Higgs_L4[6][4][5][7] = L5;
  Curvature_Higgs_L4[6][4][6][4] = L3 + L4 + L5;
  Curvature_Higgs_L4[6][4][7][5] = L5;
  Curvature_Higgs_L4[6][5][0][3] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][5][1][2] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[6][5][2][1] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[6][5][3][0] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[6][5][4][7] = L5;
  Curvature_Higgs_L4[6][5][5][6] = L3 + L4 - L5;
  Curvature_Higgs_L4[6][5][6][5] = L3 + L4 - L5;
  Curvature_Higgs_L4[6][5][7][4] = L5;
  Curvature_Higgs_L4[6][6][0][0] = L3;
  Curvature_Higgs_L4[6][6][1][1] = L3;
  Curvature_Higgs_L4[6][6][2][2] = L2;
  Curvature_Higgs_L4[6][6][3][3] = L2;
  Curvature_Higgs_L4[6][6][4][4] = L3 + L4 + L5;
  Curvature_Higgs_L4[6][6][5][5] = L3 + L4 - L5;
  Curvature_Higgs_L4[6][6][6][6] = 0.3e1 * L2;
  Curvature_Higgs_L4[6][6][7][7] = L2;
  Curvature_Higgs_L4[6][6][8][8] = L8;
  Curvature_Higgs_L4[6][7][4][5] = L5;
  Curvature_Higgs_L4[6][7][5][4] = L5;
  Curvature_Higgs_L4[6][7][6][7] = L2;
  Curvature_Higgs_L4[6][7][7][6] = L2;
  Curvature_Higgs_L4[6][8][6][8] = L8;
  Curvature_Higgs_L4[6][8][8][6] = L8;
  Curvature_Higgs_L4[7][0][0][7] = L3;
  Curvature_Higgs_L4[7][0][2][5] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][0][3][4] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[7][0][4][3] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[7][0][5][2] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][0][7][0] = L3;
  Curvature_Higgs_L4[7][1][1][7] = L3;
  Curvature_Higgs_L4[7][1][2][4] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][1][3][5] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][1][4][2] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][1][5][3] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][1][7][1] = L3;
  Curvature_Higgs_L4[7][2][0][5] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][2][1][4] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][2][2][7] = L2;
  Curvature_Higgs_L4[7][2][4][1] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][2][5][0] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][2][7][2] = L2;
  Curvature_Higgs_L4[7][3][0][4] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[7][3][1][5] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][3][3][7] = L2;
  Curvature_Higgs_L4[7][3][4][0] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[7][3][5][1] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][3][7][3] = L2;
  Curvature_Higgs_L4[7][4][0][3] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[7][4][1][2] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][4][2][1] = -L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][4][3][0] = L4 / 0.2e1 - L5 / 0.2e1;
  Curvature_Higgs_L4[7][4][4][7] = L3 + L4 - L5;
  Curvature_Higgs_L4[7][4][5][6] = L5;
  Curvature_Higgs_L4[7][4][6][5] = L5;
  Curvature_Higgs_L4[7][4][7][4] = L3 + L4 - L5;
  Curvature_Higgs_L4[7][5][0][2] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][5][1][3] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][5][2][0] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][5][3][1] = L4 / 0.2e1 + L5 / 0.2e1;
  Curvature_Higgs_L4[7][5][4][6] = L5;
  Curvature_Higgs_L4[7][5][5][7] = L3 + L4 + L5;
  Curvature_Higgs_L4[7][5][6][4] = L5;
  Curvature_Higgs_L4[7][5][7][5] = L3 + L4 + L5;
  Curvature_Higgs_L4[7][6][4][5] = L5;
  Curvature_Higgs_L4[7][6][5][4] = L5;
  Curvature_Higgs_L4[7][6][6][7] = L2;
  Curvature_Higgs_L4[7][6][7][6] = L2;
  Curvature_Higgs_L4[7][7][0][0] = L3;
  Curvature_Higgs_L4[7][7][1][1] = L3;
  Curvature_Higgs_L4[7][7][2][2] = L2;
  Curvature_Higgs_L4[7][7][3][3] = L2;
  Curvature_Higgs_L4[7][7][4][4] = L3 + L4 - L5;
  Curvature_Higgs_L4[7][7][5][5] = L3 + L4 + L5;
  Curvature_Higgs_L4[7][7][6][6] = L2;
  Curvature_Higgs_L4[7][7][7][7] = 0.3e1 * L2;
  Curvature_Higgs_L4[7][7][8][8] = L8;
  Curvature_Higgs_L4[7][8][7][8] = L8;
  Curvature_Higgs_L4[7][8][8][7] = L8;
  Curvature_Higgs_L4[8][0][0][8] = L7;
  Curvature_Higgs_L4[8][0][8][0] = L7;
  Curvature_Higgs_L4[8][1][1][8] = L7;
  Curvature_Higgs_L4[8][1][8][1] = L7;
  Curvature_Higgs_L4[8][2][2][8] = L8;
  Curvature_Higgs_L4[8][2][8][2] = L8;
  Curvature_Higgs_L4[8][3][3][8] = L8;
  Curvature_Higgs_L4[8][3][8][3] = L8;
  Curvature_Higgs_L4[8][4][4][8] = L7;
  Curvature_Higgs_L4[8][4][8][4] = L7;
  Curvature_Higgs_L4[8][5][5][8] = L7;
  Curvature_Higgs_L4[8][5][8][5] = L7;
  Curvature_Higgs_L4[8][6][6][8] = L8;
  Curvature_Higgs_L4[8][6][8][6] = L8;
  Curvature_Higgs_L4[8][7][7][8] = L8;
  Curvature_Higgs_L4[8][7][8][7] = L8;
  Curvature_Higgs_L4[8][8][0][0] = L7;
  Curvature_Higgs_L4[8][8][1][1] = L7;
  Curvature_Higgs_L4[8][8][2][2] = L8;
  Curvature_Higgs_L4[8][8][3][3] = L8;
  Curvature_Higgs_L4[8][8][4][4] = L7;
  Curvature_Higgs_L4[8][8][5][5] = L7;
  Curvature_Higgs_L4[8][8][6][6] = L8;
  Curvature_Higgs_L4[8][8][7][7] = L8;
  Curvature_Higgs_L4[8][8][8][8] = 0.6e1 * L6;

  Curvature_Gauge_G2H2[0][0][0][0] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][4][4] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][5][5] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][6][6] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][7][7] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][0][4] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][1][5] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][2][6] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][3][7] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][4][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][5][1] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][6][2] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][7][3] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][0][0] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][4][4] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][5][5] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][6][6] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][7][7] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][0][5] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][1][4] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][2][7] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][3][6] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][4][1] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][5][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][6][3] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][7][2] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][0][0] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][4][4] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][5][5] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][6][6] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][7][7] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][0][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][1][1] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][2][2] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][3][3] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][4][4] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][5][5] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][6][6] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][7][7] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][0][4] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][1][5] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][2][6] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][3][7] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][4][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][5][1] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][6][2] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][7][3] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][0][5] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][1][4] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][2][7] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][3][6] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][4][1] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][5][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][6][3] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][7][2] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][0][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][1][1] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][2][2] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][3][3] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][4][4] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][5][5] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][6][6] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][7][7] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][3][0][0] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][1][1] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][2][2] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][3][3] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][4][4] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][5][5] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][6][6] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][7][7] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;

  std::complex<double> V11, V12, V13, V21, V22, V23, V31, V32, V33;
  V11 = SMConstants.C_Vud;
  V12 = SMConstants.C_Vus;
  V13 = SMConstants.C_Vub;
  V21 = SMConstants.C_Vcd;
  V22 = SMConstants.C_Vcs;
  V23 = SMConstants.C_Vcb;
  V31 = SMConstants.C_Vtd;
  V32 = SMConstants.C_Vts;
  V33 = SMConstants.C_Vtb;

  Curvature_Quark_F2H1[0][1][4]  = 0.1e1 / v1 * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[0][1][5]  = -II / v1 * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[0][3][0]  = 0.1e1 / v1 * SMConstants.C_MassDown * V11;
  Curvature_Quark_F2H1[0][3][1]  = II / v1 * SMConstants.C_MassDown * V11;
  Curvature_Quark_F2H1[0][7][0]  = V12 / v1 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[0][7][1]  = II * V12 / v1 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[0][11][0] = V13 / v1 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[0][11][1] = II * V13 / v1 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[1][0][4]  = 0.1e1 / v1 * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[1][0][5]  = -II / v1 * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[1][2][0] =
      -0.1e1 / v1 * SMConstants.C_MassUp * std::conj(V11);
  Curvature_Quark_F2H1[1][2][1] =
      II / v1 * SMConstants.C_MassUp * std::conj(V11);
  Curvature_Quark_F2H1[1][6][0] =
      -0.1e1 / v1 * SMConstants.C_MassUp * std::conj(V12);
  Curvature_Quark_F2H1[1][6][1] =
      II / v1 * SMConstants.C_MassUp * std::conj(V12);
  Curvature_Quark_F2H1[1][10][0] =
      -0.1e1 / v1 * SMConstants.C_MassUp * std::conj(V13);
  Curvature_Quark_F2H1[1][10][1] =
      II / v1 * SMConstants.C_MassUp * std::conj(V13);
  Curvature_Quark_F2H1[2][1][0] =
      -0.1e1 / v1 * SMConstants.C_MassUp * std::conj(V11);
  Curvature_Quark_F2H1[2][1][1] =
      II / v1 * SMConstants.C_MassUp * std::conj(V11);
  Curvature_Quark_F2H1[2][3][4] = 0.1e1 / v1 * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[2][3][5] = II / v1 * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[2][5][0] =
      -0.1e1 / v1 * SMConstants.C_MassCharm * std::conj(V21);
  Curvature_Quark_F2H1[2][5][1] =
      II / v1 * SMConstants.C_MassCharm * std::conj(V21);
  Curvature_Quark_F2H1[2][9][0] =
      -0.1e1 / v1 * SMConstants.C_MassTop * std::conj(V31);
  Curvature_Quark_F2H1[2][9][1] =
      II / v1 * SMConstants.C_MassTop * std::conj(V31);
  Curvature_Quark_F2H1[3][0][0]  = 0.1e1 / v1 * SMConstants.C_MassDown * V11;
  Curvature_Quark_F2H1[3][0][1]  = II / v1 * SMConstants.C_MassDown * V11;
  Curvature_Quark_F2H1[3][2][4]  = 0.1e1 / v1 * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][2][5]  = II / v1 * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][4][0]  = V21 / v1 * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][4][1]  = II * V21 / v1 * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][8][0]  = 0.1e1 / v1 * SMConstants.C_MassDown * V31;
  Curvature_Quark_F2H1[3][8][1]  = II / v1 * SMConstants.C_MassDown * V31;
  Curvature_Quark_F2H1[4][3][0]  = V21 / v1 * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[4][3][1]  = II * V21 / v1 * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[4][5][4]  = 0.1e1 / v1 * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[4][5][5]  = -II / v1 * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[4][7][0]  = V22 / v1 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][7][1]  = II * V22 / v1 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][11][0] = V23 / v1 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[4][11][1] = II * V23 / v1 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][2][0] =
      -0.1e1 / v1 * SMConstants.C_MassCharm * std::conj(V21);
  Curvature_Quark_F2H1[5][2][1] =
      II / v1 * SMConstants.C_MassCharm * std::conj(V21);
  Curvature_Quark_F2H1[5][4][4] = 0.1e1 / v1 * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[5][4][5] = -II / v1 * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[5][6][0] =
      -0.1e1 / v1 * SMConstants.C_MassCharm * std::conj(V22);
  Curvature_Quark_F2H1[5][6][1] =
      II / v1 * SMConstants.C_MassCharm * std::conj(V22);
  Curvature_Quark_F2H1[5][10][0] =
      -0.1e1 / v1 * SMConstants.C_MassCharm * std::conj(V23);
  Curvature_Quark_F2H1[5][10][1] =
      II / v1 * SMConstants.C_MassCharm * std::conj(V23);
  Curvature_Quark_F2H1[6][1][0] =
      -0.1e1 / v1 * SMConstants.C_MassUp * std::conj(V12);
  Curvature_Quark_F2H1[6][1][1] =
      II / v1 * SMConstants.C_MassUp * std::conj(V12);
  Curvature_Quark_F2H1[6][5][0] =
      -0.1e1 / v1 * SMConstants.C_MassCharm * std::conj(V22);
  Curvature_Quark_F2H1[6][5][1] =
      II / v1 * SMConstants.C_MassCharm * std::conj(V22);
  Curvature_Quark_F2H1[6][7][4] = 0.1e1 / v1 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[6][7][5] = II / v1 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[6][9][0] =
      -0.1e1 / v1 * SMConstants.C_MassTop * std::conj(V32);
  Curvature_Quark_F2H1[6][9][1] =
      II / v1 * SMConstants.C_MassTop * std::conj(V32);
  Curvature_Quark_F2H1[7][0][0]  = V12 / v1 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[7][0][1]  = II * V12 / v1 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[7][4][0]  = V22 / v1 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[7][4][1]  = II * V22 / v1 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[7][6][4]  = 0.1e1 / v1 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[7][6][5]  = II / v1 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[7][8][0]  = V32 / v1 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[7][8][1]  = II * V32 / v1 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[8][3][0]  = 0.1e1 / v1 * SMConstants.C_MassDown * V31;
  Curvature_Quark_F2H1[8][3][1]  = II / v1 * SMConstants.C_MassDown * V31;
  Curvature_Quark_F2H1[8][7][0]  = V32 / v1 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[8][7][1]  = II * V32 / v1 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[8][9][4]  = 0.1e1 / v1 * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[8][9][5]  = -II / v1 * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[8][11][0] = V33 / v1 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[8][11][1] = II * V33 / v1 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[9][2][0] =
      -0.1e1 / v1 * SMConstants.C_MassTop * std::conj(V31);
  Curvature_Quark_F2H1[9][2][1] =
      II / v1 * SMConstants.C_MassTop * std::conj(V31);
  Curvature_Quark_F2H1[9][6][0] =
      -0.1e1 / v1 * SMConstants.C_MassTop * std::conj(V32);
  Curvature_Quark_F2H1[9][6][1] =
      II / v1 * SMConstants.C_MassTop * std::conj(V32);
  Curvature_Quark_F2H1[9][8][4] = 0.1e1 / v1 * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[9][8][5] = -II / v1 * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[9][10][0] =
      -0.1e1 / v1 * SMConstants.C_MassTop * std::conj(V33);
  Curvature_Quark_F2H1[9][10][1] =
      II / v1 * SMConstants.C_MassTop * std::conj(V33);
  Curvature_Quark_F2H1[10][1][0] =
      -0.1e1 / v1 * SMConstants.C_MassUp * std::conj(V13);
  Curvature_Quark_F2H1[10][1][1] =
      II / v1 * SMConstants.C_MassUp * std::conj(V13);
  Curvature_Quark_F2H1[10][5][0] =
      -0.1e1 / v1 * SMConstants.C_MassCharm * std::conj(V23);
  Curvature_Quark_F2H1[10][5][1] =
      II / v1 * SMConstants.C_MassCharm * std::conj(V23);
  Curvature_Quark_F2H1[10][9][0] =
      -0.1e1 / v1 * SMConstants.C_MassTop * std::conj(V33);
  Curvature_Quark_F2H1[10][9][1] =
      II / v1 * SMConstants.C_MassTop * std::conj(V33);
  Curvature_Quark_F2H1[10][11][4] = 0.1e1 / v1 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[10][11][5] = II / v1 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[11][0][0]  = V13 / v1 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[11][0][1]  = II * V13 / v1 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[11][4][0]  = V23 / v1 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[11][4][1]  = II * V23 / v1 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[11][8][0]  = V33 / v1 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[11][8][1]  = II * V33 / v1 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[11][10][4] = 0.1e1 / v1 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[11][10][5] = II / v1 * SMConstants.C_MassBottom;

  Curvature_Lepton_F2H1[0][1][4] = 0.1e1 / v1 * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[0][1][5] = II / v1 * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][0][4] = 0.1e1 / v1 * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][0][5] = II / v1 * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][6][0] = 0.1e1 / v1 * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][6][1] = II / v1 * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[2][3][4] = 0.1e1 / v1 * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[2][3][5] = II / v1 * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][2][4] = 0.1e1 / v1 * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][2][5] = II / v1 * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][7][0] = 0.1e1 / v1 * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][7][1] = II / v1 * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[4][5][4] = 0.1e1 / v1 * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[4][5][5] = II / v1 * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][4][4] = 0.1e1 / v1 * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][4][5] = II / v1 * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][8][0] = 0.1e1 / v1 * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][8][1] = II / v1 * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[6][1][0] = 0.1e1 / v1 * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[6][1][1] = II / v1 * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[7][3][0] = 0.1e1 / v1 * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[7][3][1] = II / v1 * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[8][5][0] = 0.1e1 / v1 * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[8][5][1] = II / v1 * SMConstants.C_MassTau;
}

bool Class_Potential_CPintheDark::CalculateDebyeSimplified()
{
  return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool Class_Potential_CPintheDark::CalculateDebyeGaugeSimplified()
{
  /*
   * Use this function if you calculated the Debye corrections to the gauge mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeGauge[NGauge][NGauge]
   */
  return false;
}
double
Class_Potential_CPintheDark::VTreeSimplified(const std::vector<double> &v) const
{
  double res = 0;
  if (not UseVTreeSimplified) return 0;

  double vcb, v_1, v_2, vcp, vs;

  vcb = v[2];
  v_1 = v[4];
  v_2 = v[6];
  vcp = v[7];
  vs  = v[8];

  double C22 = std::pow(vcb, 2) + std::pow(v_2, 2) + std::pow(vcp, 2);

  res += 0.5 * m11s * std::pow(v_1, 2);
  res += 0.5 * m22s * C22;
  res += 0.5 * mSs * std::pow(vs, 2);
  res += ReA * v_1 * v_2 * vs;
  res += -ImA * v_1 * vcp * vs;
  res += 1.0 / 8.0 * L1 * std::pow(v_1, 4);
  res += 1.0 / 8.0 * L2 * std::pow(C22, 2);
  res += 1.0 / 4.0 * L3 * std::pow(v_1, 2) * C22;
  res +=
      1.0 / 4.0 * L4 * std::pow(v_1, 2) * (std::pow(v_2, 2) + std::pow(vcp, 2));
  res +=
      1.0 / 4.0 * L5 * std::pow(v_1, 2) * (std::pow(v_2, 2) - std::pow(vcp, 2));
  res += 1.0 / 4.0 * L6 * std::pow(vs, 4);
  res += 1.0 / 4.0 * L7 * std::pow(v_1, 2) * std::pow(vs, 2);
  res += 1.0 / 4.0 * L8 * C22 * std::pow(vs, 2);

  return res;
}

double Class_Potential_CPintheDark::VCounterSimplified(
    const std::vector<double> &v) const
{
  double res = 0;
  if (not UseVCounterSimplified) return 0;

  double vcb, v_1, v_2, vcp, vs;

  vcb = v[2];
  v_1 = v[4];
  v_2 = v[6];
  vcp = v[7];
  vs  = v[8];

  double C22 = std::pow(vcb, 2) + std::pow(v_2, 2) + std::pow(vcp, 2);

  res += 0.5 * dm11s * std::pow(v_1, 2);
  res += 0.5 * dm22s * C22;
  res += 0.5 * dmSs * std::pow(vs, 2);
  res += dReA * v_1 * v_2 * vs;
  res += -dImA * v_1 * vcp * vs;
  res += 1.0 / 8.0 * dL1 * std::pow(v_1, 4);
  res += 1.0 / 8.0 * dL2 * std::pow(C22, 2);
  res += 1.0 / 4.0 * dL3 * std::pow(v_1, 2) * C22;
  res += 1.0 / 4.0 * dL4 * std::pow(v_1, 2) *
         (std::pow(v_2, 2) + std::pow(vcp, 2));
  res += 1.0 / 4.0 * dL5 * std::pow(v_1, 2) *
         (std::pow(v_2, 2) - std::pow(vcp, 2));
  res += 1.0 / 4.0 * dL6 * std::pow(vs, 4);
  res += 1.0 / 4.0 * dL7 * std::pow(v_1, 2) * std::pow(vs, 2);
  res += 1.0 / 4.0 * dL8 * C22 * std::pow(vs, 2);
  res += dTCB * vcb;
  res += dT1 * v_1;
  res += dT2 * v_2;
  res += dTCP * vcp;
  res += dTS * vs;
  res += -0.5 * dImL5 * std::pow(v_1, 2) * v_2 * vcp;

  return res;
}

void Class_Potential_CPintheDark::Debugging(const std::vector<double> &input,
                                            std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
