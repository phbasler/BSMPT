// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/models/ClassPotentialCxSM.h>
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
Class_CxSM::Class_CxSM()
{
  Model = ModelID::ModelIDs::CXSM; // global int constant which will be used to
                                   // tell the program which model is called
  NNeutralHiggs = 4;               // number of neutral Higgs bosons at T = 0
  NChargedHiggs = 2; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar   = 9 + 3; // number of parameters in the tree-Level Lagrangian
  nParCT = 9 + 6; // number of parameters in the counterterm potential

  nVEV = 3; // number of VEVs to minimize the potential

  NHiggs = NNeutralHiggs + NChargedHiggs;

  VevOrder.resize(nVEV);
  // Here you have to tell which scalar field gets which VEV.
  VevOrder[0] = 3;
  VevOrder[1] = 4;
  VevOrder[2] = 5;

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_CxSM::~Class_CxSM()
{
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_CxSM::addLegendCT() const
{
  std::vector<std::string> labels;
  labels.push_back("dms");
  labels.push_back("dlambda");
  labels.push_back("ddelta2");
  labels.push_back("db2");
  labels.push_back("dd2");
  labels.push_back("dReb1");
  labels.push_back("dImb1");
  labels.push_back("dRea1");
  labels.push_back("dIma1");
  labels.push_back("dT1");
  labels.push_back("dT2");
  labels.push_back("dT3");
  labels.push_back("dT4");
  labels.push_back("dT5");
  labels.push_back("dT6");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * input file
 */
std::vector<std::string> Class_CxSM::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c");
  labels.push_back("v_c");
  labels.push_back("omega_c/T_c");
  labels.push_back("omega_c");
  labels.push_back("omega_sc");
  labels.push_back("omega_ac");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * Higgs couplings. Use this to complement the legend of the given input file
 *
 */
std::vector<std::string> Class_CxSM::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;
  particles.resize(NHiggs);
  // here you have to define the particle names in the vector particles

  particles[0] = "G+";
  particles[1] = "G-";
  particles[2] = "G0";
  particles[3] = "H1";
  particles[4] = "H2";
  particles[5] = "H3";

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
std::vector<std::string> Class_CxSM::addLegendVEV() const
{
  std::vector<std::string> labels;
  labels.push_back("omega");
  labels.push_back("omega_s");
  labels.push_back("omega_a");
  return labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_CxSM::ReadAndSet(const std::string &linestr,
                            std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 12; k++)
  {
    ss >> tmp;

    if (k == 1)
      par[0] = tmp; // v_h
    else if (k == 2)
      par[1] = tmp; // v_s
    else if (k == 3)
      par[2] = tmp; // v_a
    else if (k == 4)
      par[3] = tmp; // m^2
    else if (k == 5)
      par[4] = tmp; // lambda
    else if (k == 6)
      par[5] = tmp; // delta_2
    else if (k == 7)
      par[6] = tmp; // b2
    else if (k == 8)
      par[7] = tmp; // d2
    else if (k == 9)
      par[8] = tmp; // Reb1;
    else if (k == 10)
      par[9] = tmp; // Imb1
    else if (k == 11)
      par[10] = tmp; // Rea1;
    else if (k == 12)
      par[11] = tmp; // Ima1
  }

  set_gen(par);
  return;
}

void Class_CxSM::set_gen(const std::vector<double> &par)
{

  vh     = par[0];
  vs     = par[1];
  va     = par[2];
  msq    = par[3];
  lambda = par[4];
  delta2 = par[5];
  b2     = par[6];
  d2     = par[7];
  Reb1   = par[8];
  Imb1   = par[9];
  Rea1   = par[10];
  Ima1   = par[11];

  if (va != 0 and vs != 0)
  {
    Rea1 = (-sqrt(0.2e1) * d2 * pow(va, 0.3e1) * vs -
            sqrt(0.2e1) * d2 * va * pow(vs, 0.3e1) -
            sqrt(0.2e1) * delta2 * va * vh * vh * vs + 0.4e1 * Ima1 * vs +
            Imb1 * sqrt(0.2e1) * va * va + Imb1 * sqrt(0.2e1) * vs * vs -
            0.2e1 * sqrt(0.2e1) * b2 * va * vs) /
           va / 0.4e1;
    Reb1 = 0.1e1 / va *
           (d2 * pow(va, 0.3e1) + d2 * va * vs * vs + delta2 * va * vh * vh -
            0.4e1 * Ima1 * sqrt(0.2e1) - 0.2e1 * Imb1 * vs + 0.2e1 * b2 * va) /
           0.2e1;
    msq = -delta2 * va * va / 0.2e1 - delta2 * vs * vs / 0.2e1 -
          lambda * vh * vh / 0.2e1;
  }
  else if (va != 0 and vs == 0)
  {
    Rea1 = sqrt(0.2e1) * Imb1 * va / 0.4e1;
    Reb1 = 0.1e1 / va *
           (d2 * pow(va, 0.3e1) + delta2 * va * vh * vh -
            0.4e1 * Ima1 * sqrt(0.2e1) + 0.2e1 * b2 * va) /
           0.2e1;
    msq = -delta2 * va * va / 0.2e1 - lambda * vh * vh / 0.2e1;
  }
  else if (va == 0 and vs != 0)
  {

    Ima1 = -sqrt(0.2e1) * Imb1 * vs / 0.4e1;
    Reb1 = -0.1e1 / vs *
           (d2 * pow(vs, 0.3e1) + delta2 * vh * vh * vs +
            0.4e1 * Rea1 * sqrt(0.2e1) + 0.2e1 * b2 * vs) /
           0.2e1;
    msq = -delta2 * vs * vs / 0.2e1 - lambda * vh * vh / 0.2e1;
  }
  else // va == 0 and vs == 0
  {
    Ima1 = 0;
    Rea1 = 0;
    msq  = -lambda * vh * vh / 0.2e1;
  }

  scale = vh;

  vevTreeMin.resize(nVEV);
  vevTree.resize(NHiggs);

  // Here you have to set the vector vevTreeMin. The vector vevTree will then be
  // set by the function MinimizeOrderVEV
  vevTreeMin[0] = vh;
  vevTreeMin[1] = vs;
  vevTreeMin[2] = va;

  vevTree = MinimizeOrderVEV(vevTreeMin);
  if (!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the
 * entries of Curvature_Higgs_CT_L1 to Curvature_Higgs_CT_L4.
 */
void Class_CxSM::set_CT_Pot_Par(const std::vector<double> &par)
{

  dmsq    = par[0];
  dlambda = par[1];
  ddelta2 = par[2];
  db2     = par[3];
  dd2     = par[4];
  dReb1   = par[5];
  dImb1   = par[6];
  dRea1   = par[7];
  dIma1   = par[8];
  dT1     = par[9];
  dT2     = par[10];
  dT3     = par[11];
  dT4     = par[12];
  dT5     = par[13];
  dT6     = par[14];

  Curvature_Higgs_CT_L1[0] = dT1;
  Curvature_Higgs_CT_L1[1] = dT2;
  Curvature_Higgs_CT_L1[2] = dT3;
  Curvature_Higgs_CT_L1[3] = dT4;
  Curvature_Higgs_CT_L1[4] = sqrt(0.2e1) * dRea1 + dT5;
  Curvature_Higgs_CT_L1[5] = -sqrt(0.2e1) * dIma1 + dT6;

  Curvature_Higgs_CT_L2[0][0] = dmsq / 0.2e1;
  Curvature_Higgs_CT_L2[1][1] = dmsq / 0.2e1;
  Curvature_Higgs_CT_L2[2][2] = dmsq / 0.2e1;
  Curvature_Higgs_CT_L2[3][3] = dmsq / 0.2e1;
  Curvature_Higgs_CT_L2[4][4] = db2 / 0.2e1 + dReb1 / 0.2e1;
  Curvature_Higgs_CT_L2[4][5] = -dImb1 / 0.2e1;
  Curvature_Higgs_CT_L2[5][4] = -dImb1 / 0.2e1;
  Curvature_Higgs_CT_L2[5][5] = db2 / 0.2e1 - dReb1 / 0.2e1;

  Curvature_Higgs_CT_L4[0][0][0][0] = 0.3e1 / 0.2e1 * dlambda;
  Curvature_Higgs_CT_L4[0][0][1][1] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][0][2][2] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][0][3][3] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][0][4][4] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[0][0][5][5] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[0][1][0][1] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][1][1][0] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][0][2] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][2][0] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][0][3] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][3][0] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][4][0][4] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[0][4][4][0] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[0][5][0][5] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[0][5][5][0] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[1][0][0][1] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][0][1][0] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][1][0][0] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][1][1][1] = 0.3e1 / 0.2e1 * dlambda;
  Curvature_Higgs_CT_L4[1][1][2][2] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][1][3][3] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][1][4][4] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[1][1][5][5] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][1][2] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][2][1] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][1][3] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][3][1] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][4][1][4] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[1][4][4][1] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[1][5][1][5] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[1][5][5][1] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[2][0][0][2] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][0][2][0] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][1][1][2] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][1][2][1] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][2][0][0] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][2][1][1] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][2][2][2] = 0.3e1 / 0.2e1 * dlambda;
  Curvature_Higgs_CT_L4[2][2][3][3] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][2][4][4] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[2][2][5][5] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[2][3][2][3] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][3][3][2] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][4][2][4] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[2][4][4][2] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[2][5][2][5] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[2][5][5][2] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[3][0][0][3] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[3][0][3][0] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[3][1][1][3] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[3][1][3][1] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[3][2][2][3] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[3][2][3][2] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[3][3][0][0] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[3][3][1][1] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[3][3][2][2] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[3][3][3][3] = 0.3e1 / 0.2e1 * dlambda;
  Curvature_Higgs_CT_L4[3][3][4][4] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[3][3][5][5] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[3][4][3][4] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[3][4][4][3] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[3][5][3][5] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[3][5][5][3] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[4][0][0][4] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[4][0][4][0] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[4][1][1][4] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[4][1][4][1] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[4][2][2][4] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[4][2][4][2] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[4][3][3][4] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[4][3][4][3] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[4][4][0][0] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[4][4][1][1] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[4][4][2][2] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[4][4][3][3] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[4][4][4][4] = 0.3e1 / 0.2e1 * dd2;
  Curvature_Higgs_CT_L4[4][4][5][5] = dd2 / 0.2e1;
  Curvature_Higgs_CT_L4[4][5][4][5] = dd2 / 0.2e1;
  Curvature_Higgs_CT_L4[4][5][5][4] = dd2 / 0.2e1;
  Curvature_Higgs_CT_L4[5][0][0][5] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[5][0][5][0] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[5][1][1][5] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[5][1][5][1] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[5][2][2][5] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[5][2][5][2] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[5][3][3][5] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[5][3][5][3] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[5][4][4][5] = dd2 / 0.2e1;
  Curvature_Higgs_CT_L4[5][4][5][4] = dd2 / 0.2e1;
  Curvature_Higgs_CT_L4[5][5][0][0] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[5][5][1][1] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[5][5][2][2] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[5][5][3][3] = ddelta2 / 0.2e1;
  Curvature_Higgs_CT_L4[5][5][4][4] = dd2 / 0.2e1;
  Curvature_Higgs_CT_L4[5][5][5][5] = 0.3e1 / 0.2e1 * dd2;
}

/**
 * console output of all Parameters
 */
void Class_CxSM::write() const
{
  std::stringstream ss;
  ss << "The parameters are : " << std::endl;
  ss << "\tlambda = " << lambda << std::endl
     << "\tm^2 = " << msq << std::endl
     << "\tdelta_2 = " << delta2 << std::endl
     << "\tb2 = " << b2 << std::endl
     << "\td2 = " << d2 << std::endl
     << "\tReb1 = " << Reb1 << "\n"
     << "\tImb1 = " << Imb1 << "\n"
     << "\tRea1 = " << Rea1 << "\n"
     << "\tIma1 = " << Ima1 << "\n"
     << "\tvh = " << vh << "\n"
     << "\tva = " << va << "\n"
     << "\tvs = " << vs << std::endl;

  ss << "The counterterm parameters are : " << std::endl;
  ss << "\tdlambda = " << dlambda << std::endl
     << "\tdm^2 = " << dmsq << std::endl
     << "\tddelta_2 = " << ddelta2 << std::endl
     << "\tdb2 = " << db2 << "\n"
     << "\tdd2 = " << dd2 << "\n"
     << "\tdReb1 = " << dReb1 << "\n"
     << "\tdImb1 = " << dImb1 << "\n"
     << "\tdRea1 = " << dRea1 << "\n"
     << "\tdIma1 = " << dIma1 << "\n"
     << "\tdT1 = " << dT1 << "\n"
     << "\tdT2 = " << dT2 << "\n"
     << "\tdT3 = " << dT3 << "\n"
     << "\tdT4 = " << dT4 << "\n"
     << "\tdT5 = " << dT5 << "\n"
     << "\tdT6 = " << dT6 << std::endl;

  ss << "The scale is given by mu = " << scale << " GeV " << std::endl;

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
    }
  }

  int posN[3];
  posN[0]         = 3;
  posN[1]         = 4;
  posN[2]         = 5;
  int posGCharged = 0, posG0 = 0;
  double testsum             = 0;
  const double ZeroThreshold = 1e-5;
  for (int i = 0; i < 3; i++)
  {
    testsum = std::abs(HiggsRot(i, 0)) + std::abs(HiggsRot(i, 1));
    if (testsum > ZeroThreshold and posGCharged == 0)
    {
      posGCharged = i;
    }
    testsum = std::abs(HiggsRot(i, 2));
    if (testsum > ZeroThreshold) posG0 = i;
  }

  std::vector<double> HiggsMasses;
  HiggsMasses = HiggsMassesSquared(vevTree, 0);

  double MhUp = 0, MhDown = 0, MSM = 0;
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

  MatrixXd NeutralMatrix(3, 3);
  for (int j = 0; j < 3; j++)
  {
    for (int i = 0; i < 3; i++)
      NeutralMatrix(i, j) = HiggsRot(posN[i], j + 3);
  }

  ss << "The mass spectrum is given by :\n";
  ss << "m_{G^+}^2 = " << HiggsMasses[posGCharged] << " GeV^2 \n"
     << "m_{G^0}^2 = " << HiggsMasses[posG0] << " GeV^2 \n"
     << "m_{H_SM} = " << MSM << " GeV \n"
     << "m_{H_l} = " << MhDown << " GeV \n"
     << "m_{H_h} = " << MhUp << " GeV \n";
  ss << "The neutral mixing Matrix is given by :\n";
  bool IsNegative = NeutralMatrix(0, 1) < 0;
  ss << "H_{SM} = " << NeutralMatrix(0, 0) << " h ";
  if (IsNegative)
    ss << "-";
  else
    ss << "+";
  ss << std::abs(NeutralMatrix(0, 1)) << " a ";
  IsNegative = NeutralMatrix(0, 2) < 0;
  if (IsNegative)
    ss << "-";
  else
    ss << "+";
  ss << std::abs(NeutralMatrix(0, 2)) << " s \n"
     << "H_{l} = " << NeutralMatrix(1, 0) << " h ";
  IsNegative = NeutralMatrix(1, 1) < 0;
  if (IsNegative)
    ss << "-";
  else
    ss << "+";
  ss << std::abs(NeutralMatrix(1, 1)) << " a ";
  IsNegative = NeutralMatrix(1, 2) < 0;
  if (IsNegative)
    ss << "-";
  else
    ss << "+";
  ss << std::abs(NeutralMatrix(1, 2)) << " s \n"
     << "H_{h} = " << NeutralMatrix(2, 0) << " h ";
  IsNegative = NeutralMatrix(2, 1) < 0;
  if (IsNegative)
    ss << "-";
  else
    ss << "+";
  ss << std::abs(NeutralMatrix(2, 1)) << " a ";
  IsNegative = NeutralMatrix(2, 2) < 0;
  if (IsNegative)
    ss << "-";
  else
    ss << "+";
  ss << std::abs(NeutralMatrix(2, 2)) << " s \n";
  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms. Here you need to work out the scheme and
 * implement the formulas.
 */
std::vector<double> Class_CxSM::calc_CT() const
{
  std::vector<double> parCT;
  parCT.resize(nParCT);

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

  double ZeroCut = 1e-5;
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    if (std::abs(NablaWeinberg[i]) <= ZeroCut) NablaWeinberg[i] = 0;
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      if (std::abs(HesseWeinberg(i, j)) <= ZeroCut) HesseWeinberg(i, j) = 0;
    }
  }

  // Here you have to use your formulas for the counterterm scheme

  if (va != 0 and vs != 0 and Ima1 != 0)
  {
    // You can either choose if Ima1 or Imb1 has a CT, in case Ima1 != 0 it will
    // get the CT, otherwise Imb1
    parCT[0] = ((va * va + vs * vs) * HesseWeinberg(3, 5) -
                3 * vh * va * (HesseWeinberg(2, 2) - HesseWeinberg(3, 3) / 3)) /
               vh / va;
    parCT[1] =
        (2 * HesseWeinberg(2, 2) - 2 * HesseWeinberg(3, 3)) * pow(vh, -2);
    parCT[2] = -2 / va / vh * HesseWeinberg(3, 5);
    parCT[3] = ((2 * va * va + 2 * vs * vs) * HesseWeinberg(4, 5) -
                (-HesseWeinberg(3, 5) * vh +
                 va * (HesseWeinberg(4, 4) + HesseWeinberg(5, 5))) *
                    vs) /
               vs / va;
    parCT[4] = -2 / vs / va * HesseWeinberg(4, 5);
    parCT[5] = HesseWeinberg(5, 5) + vs / va * HesseWeinberg(4, 5) -
               va / vs * HesseWeinberg(4, 5) - HesseWeinberg(4, 4);
    parCT[6] = 0;
    parCT[7] = (-HesseWeinberg(4, 5) * vs * vs +
                va * (vs * HesseWeinberg(4, 4) - NablaWeinberg(4))) *
               sqrt(2) / va / 2;
    parCT[8] = sqrt(2) *
               (va * va * HesseWeinberg(4, 5) -
                vs * (HesseWeinberg(5, 5) * va - NablaWeinberg(5))) /
               vs / 2;
    parCT[9]  = -NablaWeinberg(0);
    parCT[10] = -NablaWeinberg(1);
    parCT[11] = -NablaWeinberg(2);
    parCT[12] = HesseWeinberg(2, 2) * vh - NablaWeinberg(3);
    parCT[13] = 0;
    parCT[14] = 0;
  }
  else if (va != 0 and vs != 0)
  {
    parCT[0] = ((va * va + vs * vs) * HesseWeinberg(3, 5) -
                3 * vh * (HesseWeinberg(2, 2) - HesseWeinberg(3, 3) / 3) * va) /
               vh / va;
    parCT[1] =
        (2 * HesseWeinberg(2, 2) - 2 * HesseWeinberg(3, 3)) * pow(vh, -2);
    parCT[2] = -2 / va / vh * HesseWeinberg(3, 5);
    parCT[3] = (2 * HesseWeinberg(4, 5) * vs + HesseWeinberg(5, 5) * va -
                HesseWeinberg(4, 4) * va + HesseWeinberg(3, 5) * vh -
                2 * NablaWeinberg(5)) /
               va;
    parCT[4] = (-2 * HesseWeinberg(4, 5) * vs - 2 * HesseWeinberg(5, 5) * va +
                2 * NablaWeinberg(5)) /
               (va * va + vs * vs) / va;
    parCT[5] = ((-pow(va, 3) - vs * vs * va) * HesseWeinberg(4, 4) +
                (-vs * va * va + pow(vs, 3)) * HesseWeinberg(4, 5) +
                2 * HesseWeinberg(5, 5) * va * vs * vs +
                NablaWeinberg(5) * (va * va - vs * vs)) /
               (va * va + vs * vs) / va;
    parCT[6] = (2 * va * va * HesseWeinberg(4, 5) -
                2 * vs * (HesseWeinberg(5, 5) * va - NablaWeinberg(5))) /
               (va * va + vs * vs);
    parCT[7] = (HesseWeinberg(4, 5) * (va * va - vs * vs) +
                HesseWeinberg(4, 4) * va * vs - HesseWeinberg(5, 5) * va * vs -
                NablaWeinberg(4) * va + NablaWeinberg(5) * vs) *
               sqrt(2) / va / 2;
    parCT[8]  = 0;
    parCT[9]  = -NablaWeinberg(0);
    parCT[10] = -NablaWeinberg(1);
    parCT[11] = -NablaWeinberg(2);
    parCT[12] = HesseWeinberg(2, 2) * vh - NablaWeinberg(3);
    parCT[13] = 0;
    parCT[14] = 0;
  }
  else if (va != 0 and vs == 0)
  {
    parCT[0] = (-3 * HesseWeinberg(2, 2) * vh + HesseWeinberg(3, 3) * vh +
                va * HesseWeinberg(3, 5)) /
               vh;
    parCT[1] =
        (2 * HesseWeinberg(2, 2) - 2 * HesseWeinberg(3, 3)) * pow(vh, -2);
    parCT[2] = -2 / va / vh * HesseWeinberg(3, 5);
    parCT[3] = (HesseWeinberg(3, 5) * vh - HesseWeinberg(4, 4) * va +
                HesseWeinberg(5, 5) * va - 2 * NablaWeinberg(5)) /
               va;
    parCT[4] =
        (-2 * HesseWeinberg(5, 5) * va + 2 * NablaWeinberg(5)) * pow(va, -3);
    parCT[5]  = (-HesseWeinberg(4, 4) * va + NablaWeinberg(5)) / va;
    parCT[6]  = 2 * HesseWeinberg(4, 5);
    parCT[7]  = -sqrt(2) * (-va * HesseWeinberg(4, 5) + NablaWeinberg(4)) / 2;
    parCT[8]  = 0;
    parCT[9]  = -NablaWeinberg(0);
    parCT[10] = -NablaWeinberg(1);
    parCT[11] = -NablaWeinberg(2);
    parCT[12] = HesseWeinberg(2, 2) * vh - NablaWeinberg(3);
    parCT[13] = 0;
    parCT[14] = 0;
  }
  else if (va == 0 and vs != 0)
  {
    parCT[0] = (-3 * HesseWeinberg(2, 2) * vh + HesseWeinberg(3, 3) * vh +
                vs * HesseWeinberg(3, 4)) /
               vh;
    parCT[1] =
        (2 * HesseWeinberg(2, 2) - 2 * HesseWeinberg(3, 3)) * pow(vh, -2);
    parCT[2] = -2 / vh / vs * HesseWeinberg(3, 4);
    parCT[3] = (HesseWeinberg(3, 4) * vh + vs * HesseWeinberg(4, 4) -
                HesseWeinberg(5, 5) * vs - 2 * NablaWeinberg(4)) /
               vs;
    parCT[4] =
        (-2 * vs * HesseWeinberg(4, 4) + 2 * NablaWeinberg(4)) * pow(vs, -3);
    parCT[5]  = (HesseWeinberg(5, 5) * vs - NablaWeinberg(4)) / vs;
    parCT[6]  = 2 * HesseWeinberg(4, 5);
    parCT[7]  = 0;
    parCT[8]  = -sqrt(2) * (HesseWeinberg(4, 5) * vs - NablaWeinberg(5)) / 2;
    parCT[9]  = -NablaWeinberg(0);
    parCT[10] = -NablaWeinberg(1);
    parCT[11] = -NablaWeinberg(2);
    parCT[12] = HesseWeinberg(2, 2) * vh - NablaWeinberg(3);
    parCT[13] = 0;
    parCT[14] = 0;
  }
  else
  {
    // case va == 0 and vs == 0
    parCT[0] = -3 * HesseWeinberg(2, 2) + HesseWeinberg(3, 3);
    parCT[1] =
        (2 * HesseWeinberg(2, 2) - 2 * HesseWeinberg(3, 3)) * pow(vh, -2);
    parCT[2]  = 0;
    parCT[3]  = -HesseWeinberg(4, 4) - HesseWeinberg(5, 5);
    parCT[4]  = 0;
    parCT[5]  = HesseWeinberg(5, 5) - HesseWeinberg(4, 4);
    parCT[6]  = 2 * HesseWeinberg(4, 5);
    parCT[7]  = -sqrt(2) * NablaWeinberg(4) / 2;
    parCT[8]  = sqrt(2) * NablaWeinberg(5) / 2;
    parCT[9]  = -NablaWeinberg(0);
    parCT[10] = -NablaWeinberg(1);
    parCT[11] = -NablaWeinberg(2);
    parCT[12] = HesseWeinberg(2, 2) * vh - NablaWeinberg(3);
    parCT[13] = 0;
    parCT[14] = 0;
  }
  for (std::size_t i = 0; i < nParCT; i++)
  {
    if (std::abs(parCT[i]) <= ZeroCut) parCT[i] = 0;
  }

  return parCT;
}

void Class_CxSM::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsdone) CalculatePhysicalCouplings();

  std::vector<double> HiggsOrder(NHiggs);
  // Here you have to set the vector HiggsOrder. By telling e.g. HiggsOrder[0] =
  // 5 you always want your 6th lightest particle to be the first particle in
  // the vector (which has the index 5 because they are sorted by mass)

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
    }
  }

  std::size_t posGp = 0, posGm = 0, posG0 = 0;
  std::size_t posH1 = 0, posH2 = 0, posH3 = 0;
  const double ZeroThreshold = 1e-5;

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    // the rotation matrix is diagonal besides for the neutral scalars
    if (std::abs(HiggsRot(i, 0)) > ZeroThreshold)
      posGp = i;
    else if (std::abs(HiggsRot(i, 1)) > ZeroThreshold)
      posGm = i;
    else if (std::abs(HiggsRot(i, 2)) > ZeroThreshold)
      posG0 = i;

    // the neutral scalars mix
    if ((std::abs(HiggsRot(i, 3)) + std::abs(HiggsRot(i, 4)) +
         std::abs(HiggsRot(i, 5))) > ZeroThreshold)
    {
      // use that scalars are sorted by mass
      if (posH1 == 0)
      {
        posH1 = i;
      }
      else
      {
        if (posH2 == 0)
        {
          posH2 = i;
        }
        else
        {
          posH3 = i;
        }
      }
    }
  }

  // mass order: Gp, Gm, G0, H1, H2, H3
  HiggsOrder[0] = posGp;
  HiggsOrder[1] = posGm;
  HiggsOrder[2] = posG0;
  HiggsOrder[3] = posH1;
  HiggsOrder[4] = posH2;
  HiggsOrder[5] = posH3;

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

  MatrixXd HiggsRotSort(NHiggs, NHiggs);

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsRotSort.row(i) = HiggsRot.row(HiggsOrder[i]);
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

void Class_CxSM::SetCurvatureArrays()
{
  /*
   *  Here you have to set the vectors
   *  Curvature_Higgs_L1,Curvature_Higgs_L2,Curvature_Higgs_L3,Curvature_Higgs_L4
   *  Curvature_Gauge_G2H2
   *  Curvature_Quark_F2H1, Curvature_Lepton_F2H1
   *  as described in the potential in the paper.
   */

  initVectors();

  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

  Curvature_Higgs_L1[4] = Rea1 * sqrt(0.2e1);
  Curvature_Higgs_L1[5] = -Ima1 * sqrt(0.2e1);

  Curvature_Higgs_L2[0][0] = msq / 0.2e1;
  Curvature_Higgs_L2[1][1] = msq / 0.2e1;
  Curvature_Higgs_L2[2][2] = msq / 0.2e1;
  Curvature_Higgs_L2[3][3] = msq / 0.2e1;
  Curvature_Higgs_L2[4][4] = Reb1 / 0.2e1 + b2 / 0.2e1;
  Curvature_Higgs_L2[4][5] = -Imb1 / 0.2e1;
  Curvature_Higgs_L2[5][4] = -Imb1 / 0.2e1;
  Curvature_Higgs_L2[5][5] = -Reb1 / 0.2e1 + b2 / 0.2e1;

  Curvature_Higgs_L4[0][0][0][0] = 0.3e1 / 0.2e1 * lambda;
  Curvature_Higgs_L4[0][0][1][1] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][0][2][2] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][0][3][3] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][0][4][4] = delta2 / 0.2e1;
  Curvature_Higgs_L4[0][0][5][5] = delta2 / 0.2e1;
  Curvature_Higgs_L4[0][1][0][1] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][1][1][0] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][2][0][2] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][2][2][0] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][3][0][3] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][3][3][0] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][4][0][4] = delta2 / 0.2e1;
  Curvature_Higgs_L4[0][4][4][0] = delta2 / 0.2e1;
  Curvature_Higgs_L4[0][5][0][5] = delta2 / 0.2e1;
  Curvature_Higgs_L4[0][5][5][0] = delta2 / 0.2e1;
  Curvature_Higgs_L4[1][0][0][1] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][0][1][0] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][1][0][0] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][1][1][1] = 0.3e1 / 0.2e1 * lambda;
  Curvature_Higgs_L4[1][1][2][2] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][1][3][3] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][1][4][4] = delta2 / 0.2e1;
  Curvature_Higgs_L4[1][1][5][5] = delta2 / 0.2e1;
  Curvature_Higgs_L4[1][2][1][2] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][2][2][1] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][3][1][3] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][3][3][1] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][4][1][4] = delta2 / 0.2e1;
  Curvature_Higgs_L4[1][4][4][1] = delta2 / 0.2e1;
  Curvature_Higgs_L4[1][5][1][5] = delta2 / 0.2e1;
  Curvature_Higgs_L4[1][5][5][1] = delta2 / 0.2e1;
  Curvature_Higgs_L4[2][0][0][2] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][0][2][0] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][1][1][2] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][1][2][1] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][2][0][0] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][2][1][1] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][2][2][2] = 0.3e1 / 0.2e1 * lambda;
  Curvature_Higgs_L4[2][2][3][3] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][2][4][4] = delta2 / 0.2e1;
  Curvature_Higgs_L4[2][2][5][5] = delta2 / 0.2e1;
  Curvature_Higgs_L4[2][3][2][3] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][3][3][2] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][4][2][4] = delta2 / 0.2e1;
  Curvature_Higgs_L4[2][4][4][2] = delta2 / 0.2e1;
  Curvature_Higgs_L4[2][5][2][5] = delta2 / 0.2e1;
  Curvature_Higgs_L4[2][5][5][2] = delta2 / 0.2e1;
  Curvature_Higgs_L4[3][0][0][3] = lambda / 0.2e1;
  Curvature_Higgs_L4[3][0][3][0] = lambda / 0.2e1;
  Curvature_Higgs_L4[3][1][1][3] = lambda / 0.2e1;
  Curvature_Higgs_L4[3][1][3][1] = lambda / 0.2e1;
  Curvature_Higgs_L4[3][2][2][3] = lambda / 0.2e1;
  Curvature_Higgs_L4[3][2][3][2] = lambda / 0.2e1;
  Curvature_Higgs_L4[3][3][0][0] = lambda / 0.2e1;
  Curvature_Higgs_L4[3][3][1][1] = lambda / 0.2e1;
  Curvature_Higgs_L4[3][3][2][2] = lambda / 0.2e1;
  Curvature_Higgs_L4[3][3][3][3] = 0.3e1 / 0.2e1 * lambda;
  Curvature_Higgs_L4[3][3][4][4] = delta2 / 0.2e1;
  Curvature_Higgs_L4[3][3][5][5] = delta2 / 0.2e1;
  Curvature_Higgs_L4[3][4][3][4] = delta2 / 0.2e1;
  Curvature_Higgs_L4[3][4][4][3] = delta2 / 0.2e1;
  Curvature_Higgs_L4[3][5][3][5] = delta2 / 0.2e1;
  Curvature_Higgs_L4[3][5][5][3] = delta2 / 0.2e1;
  Curvature_Higgs_L4[4][0][0][4] = delta2 / 0.2e1;
  Curvature_Higgs_L4[4][0][4][0] = delta2 / 0.2e1;
  Curvature_Higgs_L4[4][1][1][4] = delta2 / 0.2e1;
  Curvature_Higgs_L4[4][1][4][1] = delta2 / 0.2e1;
  Curvature_Higgs_L4[4][2][2][4] = delta2 / 0.2e1;
  Curvature_Higgs_L4[4][2][4][2] = delta2 / 0.2e1;
  Curvature_Higgs_L4[4][3][3][4] = delta2 / 0.2e1;
  Curvature_Higgs_L4[4][3][4][3] = delta2 / 0.2e1;
  Curvature_Higgs_L4[4][4][0][0] = delta2 / 0.2e1;
  Curvature_Higgs_L4[4][4][1][1] = delta2 / 0.2e1;
  Curvature_Higgs_L4[4][4][2][2] = delta2 / 0.2e1;
  Curvature_Higgs_L4[4][4][3][3] = delta2 / 0.2e1;
  Curvature_Higgs_L4[4][4][4][4] = 0.3e1 / 0.2e1 * d2;
  Curvature_Higgs_L4[4][4][5][5] = d2 / 0.2e1;
  Curvature_Higgs_L4[4][5][4][5] = d2 / 0.2e1;
  Curvature_Higgs_L4[4][5][5][4] = d2 / 0.2e1;
  Curvature_Higgs_L4[5][0][0][5] = delta2 / 0.2e1;
  Curvature_Higgs_L4[5][0][5][0] = delta2 / 0.2e1;
  Curvature_Higgs_L4[5][1][1][5] = delta2 / 0.2e1;
  Curvature_Higgs_L4[5][1][5][1] = delta2 / 0.2e1;
  Curvature_Higgs_L4[5][2][2][5] = delta2 / 0.2e1;
  Curvature_Higgs_L4[5][2][5][2] = delta2 / 0.2e1;
  Curvature_Higgs_L4[5][3][3][5] = delta2 / 0.2e1;
  Curvature_Higgs_L4[5][3][5][3] = delta2 / 0.2e1;
  Curvature_Higgs_L4[5][4][4][5] = d2 / 0.2e1;
  Curvature_Higgs_L4[5][4][5][4] = d2 / 0.2e1;
  Curvature_Higgs_L4[5][5][0][0] = delta2 / 0.2e1;
  Curvature_Higgs_L4[5][5][1][1] = delta2 / 0.2e1;
  Curvature_Higgs_L4[5][5][2][2] = delta2 / 0.2e1;
  Curvature_Higgs_L4[5][5][3][3] = delta2 / 0.2e1;
  Curvature_Higgs_L4[5][5][4][4] = d2 / 0.2e1;
  Curvature_Higgs_L4[5][5][5][5] = 0.3e1 / 0.2e1 * d2;

  Curvature_Gauge_G2H2[0][0][0][0] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][1][1] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][2][2] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][3][3] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][0][3] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][1][2] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][2][1] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][3][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][0][0] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][1][1] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][2][2] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][3][3] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][0][2] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][1][3] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][2][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][3][1] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][0][0] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][1][1] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][2][2] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][3][3] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][0][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][1][1] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][2][2] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][3][3] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][0][3] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][1][2] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][2][1] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][3][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][0][2] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][1][3] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][2][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][3][1] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][0][0] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][1][1] = C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][2][2] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][3][3] = -C_gs * C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][3][0][0] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][1][1] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][2][2] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][3][3] = C_gs * C_gs / 0.2e1;

  MatrixXcd YIJQc1(NQuarks, NQuarks), YIJQc2(NQuarks, NQuarks),
      YIJQc2OI(NQuarks, NQuarks), YIJQg0(NQuarks, NQuarks),
      YIJQg0OI(NQuarks, NQuarks), YIJQh1(NQuarks, NQuarks),
      YIJQh2(NQuarks, NQuarks), YIJQh3(NQuarks, NQuarks);

  std::complex<double> V11, V12, V13, V21, V22, V23, V31, V32, V33;
  V11 = C_Vud;
  V12 = C_Vus;
  V13 = C_Vub;
  V21 = C_Vcd;
  V22 = C_Vcs;
  V23 = C_Vcb;
  V31 = C_Vtd;
  V32 = C_Vts;
  V33 = C_Vtb;

  std::complex<double> VC11, VC12, VC13, VC21, VC22, VC23, VC31, VC32, VC33;
  VC11 = std::conj(C_Vud);
  VC12 = std::conj(C_Vus);
  VC13 = std::conj(C_Vub);
  VC21 = std::conj(C_Vcd);
  VC22 = std::conj(C_Vcs);
  VC23 = std::conj(C_Vcb);
  VC31 = std::conj(C_Vtd);
  VC32 = std::conj(C_Vts);
  VC33 = std::conj(C_Vtb);

  std::complex<double> II(0, 1);

  Curvature_Lepton_F2H1[0][1][2] = II / vh * C_MassElectron;
  Curvature_Lepton_F2H1[0][1][3] = 0.1e1 / vh * C_MassElectron;
  Curvature_Lepton_F2H1[1][0][2] = II / vh * C_MassElectron;
  Curvature_Lepton_F2H1[1][0][3] = 0.1e1 / vh * C_MassElectron;
  Curvature_Lepton_F2H1[1][6][0] = 0.1e1 / vh * C_MassElectron;
  Curvature_Lepton_F2H1[1][6][1] = II / vh * C_MassElectron;
  Curvature_Lepton_F2H1[2][3][2] = II / vh * C_MassMu;
  Curvature_Lepton_F2H1[2][3][3] = 0.1e1 / vh * C_MassMu;
  Curvature_Lepton_F2H1[3][2][2] = II / vh * C_MassMu;
  Curvature_Lepton_F2H1[3][2][3] = 0.1e1 / vh * C_MassMu;
  Curvature_Lepton_F2H1[3][7][0] = 0.1e1 / vh * C_MassMu;
  Curvature_Lepton_F2H1[3][7][1] = II / vh * C_MassMu;
  Curvature_Lepton_F2H1[4][5][2] = II / vh * C_MassTau;
  Curvature_Lepton_F2H1[4][5][3] = 0.1e1 / vh * C_MassTau;
  Curvature_Lepton_F2H1[5][4][2] = II / vh * C_MassTau;
  Curvature_Lepton_F2H1[5][4][3] = 0.1e1 / vh * C_MassTau;
  Curvature_Lepton_F2H1[5][8][0] = 0.1e1 / vh * C_MassTau;
  Curvature_Lepton_F2H1[5][8][1] = II / vh * C_MassTau;
  Curvature_Lepton_F2H1[6][1][0] = 0.1e1 / vh * C_MassElectron;
  Curvature_Lepton_F2H1[6][1][1] = II / vh * C_MassElectron;
  Curvature_Lepton_F2H1[7][3][0] = 0.1e1 / vh * C_MassMu;
  Curvature_Lepton_F2H1[7][3][1] = II / vh * C_MassMu;
  Curvature_Lepton_F2H1[8][5][0] = 0.1e1 / vh * C_MassTau;
  Curvature_Lepton_F2H1[8][5][1] = II / vh * C_MassTau;

  Curvature_Quark_F2H1[0][6][2]  = -II / vh * C_MassUp;
  Curvature_Quark_F2H1[0][6][3]  = 0.1e1 / vh * C_MassUp;
  Curvature_Quark_F2H1[0][9][0]  = -0.1e1 / vh * C_MassUp * conj(V11);
  Curvature_Quark_F2H1[0][9][1]  = II / vh * C_MassUp * conj(V11);
  Curvature_Quark_F2H1[0][10][0] = -0.1e1 / vh * C_MassUp * conj(V12);
  Curvature_Quark_F2H1[0][10][1] = II / vh * C_MassUp * conj(V12);
  Curvature_Quark_F2H1[0][11][0] = -0.1e1 / vh * C_MassUp * conj(V13);
  Curvature_Quark_F2H1[0][11][1] = II / vh * C_MassUp * conj(V13);
  Curvature_Quark_F2H1[1][7][2]  = -II / vh * C_MassCharm;
  Curvature_Quark_F2H1[1][7][3]  = 0.1e1 / vh * C_MassCharm;
  Curvature_Quark_F2H1[1][9][0]  = -0.1e1 / vh * C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[1][9][1]  = II / vh * C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[1][10][0] = -0.1e1 / vh * C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[1][10][1] = II / vh * C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[1][11][0] = -0.1e1 / vh * C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[1][11][1] = II / vh * C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[2][8][2]  = -II / vh * C_MassTop;
  Curvature_Quark_F2H1[2][8][3]  = 0.1e1 / vh * C_MassTop;
  Curvature_Quark_F2H1[2][9][0]  = -0.1e1 / vh * C_MassTop * conj(V31);
  Curvature_Quark_F2H1[2][9][1]  = II / vh * C_MassTop * conj(V31);
  Curvature_Quark_F2H1[2][10][0] = -0.1e1 / vh * C_MassTop * conj(V32);
  Curvature_Quark_F2H1[2][10][1] = II / vh * C_MassTop * conj(V32);
  Curvature_Quark_F2H1[2][11][0] = -0.1e1 / vh * C_MassTop * conj(V33);
  Curvature_Quark_F2H1[2][11][1] = II / vh * C_MassTop * conj(V33);
  Curvature_Quark_F2H1[3][6][0]  = 0.1e1 / vh * C_MassDown * V11;
  Curvature_Quark_F2H1[3][6][1]  = II / vh * C_MassDown * V11;
  Curvature_Quark_F2H1[3][7][0]  = V21 / vh * C_MassDown;
  Curvature_Quark_F2H1[3][7][1]  = II * V21 / vh * C_MassDown;
  Curvature_Quark_F2H1[3][8][0]  = 0.1e1 / vh * C_MassDown * V31;
  Curvature_Quark_F2H1[3][8][1]  = II / vh * C_MassDown * V31;
  Curvature_Quark_F2H1[3][9][2]  = II / vh * C_MassDown;
  Curvature_Quark_F2H1[3][9][3]  = 0.1e1 / vh * C_MassDown;
  Curvature_Quark_F2H1[4][6][0]  = 0.1e1 / vh * C_MassStrange * V12;
  Curvature_Quark_F2H1[4][6][1]  = II / vh * C_MassStrange * V12;
  Curvature_Quark_F2H1[4][7][0]  = V22 / vh * C_MassStrange;
  Curvature_Quark_F2H1[4][7][1]  = II * V22 / vh * C_MassStrange;
  Curvature_Quark_F2H1[4][8][0]  = 0.1e1 / vh * C_MassStrange * V32;
  Curvature_Quark_F2H1[4][8][1]  = II / vh * C_MassStrange * V32;
  Curvature_Quark_F2H1[4][10][2] = II / vh * C_MassStrange;
  Curvature_Quark_F2H1[4][10][3] = 0.1e1 / vh * C_MassStrange;
  Curvature_Quark_F2H1[5][6][0]  = V13 / vh * C_MassBottom;
  Curvature_Quark_F2H1[5][6][1]  = II / vh * C_MassBottom * V13;
  Curvature_Quark_F2H1[5][7][0]  = V23 / vh * C_MassBottom;
  Curvature_Quark_F2H1[5][7][1]  = II / vh * C_MassBottom * V23;
  Curvature_Quark_F2H1[5][8][0]  = V33 / vh * C_MassBottom;
  Curvature_Quark_F2H1[5][8][1]  = II / vh * C_MassBottom * V33;
  Curvature_Quark_F2H1[5][11][2] = II / vh * C_MassBottom;
  Curvature_Quark_F2H1[5][11][3] = 0.1e1 / vh * C_MassBottom;
  Curvature_Quark_F2H1[6][0][2]  = -II / vh * C_MassUp;
  Curvature_Quark_F2H1[6][0][3]  = 0.1e1 / vh * C_MassUp;
  Curvature_Quark_F2H1[6][3][0]  = 0.1e1 / vh * C_MassDown * V11;
  Curvature_Quark_F2H1[6][3][1]  = II / vh * C_MassDown * V11;
  Curvature_Quark_F2H1[6][4][0]  = 0.1e1 / vh * C_MassStrange * V12;
  Curvature_Quark_F2H1[6][4][1]  = II / vh * C_MassStrange * V12;
  Curvature_Quark_F2H1[6][5][0]  = V13 / vh * C_MassBottom;
  Curvature_Quark_F2H1[6][5][1]  = II / vh * C_MassBottom * V13;
  Curvature_Quark_F2H1[7][1][2]  = -II / vh * C_MassCharm;
  Curvature_Quark_F2H1[7][1][3]  = 0.1e1 / vh * C_MassCharm;
  Curvature_Quark_F2H1[7][3][0]  = V21 / vh * C_MassDown;
  Curvature_Quark_F2H1[7][3][1]  = II * V21 / vh * C_MassDown;
  Curvature_Quark_F2H1[7][4][0]  = V22 / vh * C_MassStrange;
  Curvature_Quark_F2H1[7][4][1]  = II * V22 / vh * C_MassStrange;
  Curvature_Quark_F2H1[7][5][0]  = V23 / vh * C_MassBottom;
  Curvature_Quark_F2H1[7][5][1]  = II / vh * C_MassBottom * V23;
  Curvature_Quark_F2H1[8][2][2]  = -II / vh * C_MassTop;
  Curvature_Quark_F2H1[8][2][3]  = 0.1e1 / vh * C_MassTop;
  Curvature_Quark_F2H1[8][3][0]  = 0.1e1 / vh * C_MassDown * V31;
  Curvature_Quark_F2H1[8][3][1]  = II / vh * C_MassDown * V31;
  Curvature_Quark_F2H1[8][4][0]  = 0.1e1 / vh * C_MassStrange * V32;
  Curvature_Quark_F2H1[8][4][1]  = II / vh * C_MassStrange * V32;
  Curvature_Quark_F2H1[8][5][0]  = V33 / vh * C_MassBottom;
  Curvature_Quark_F2H1[8][5][1]  = II / vh * C_MassBottom * V33;
  Curvature_Quark_F2H1[9][0][0]  = -0.1e1 / vh * C_MassUp * conj(V11);
  Curvature_Quark_F2H1[9][0][1]  = II / vh * C_MassUp * conj(V11);
  Curvature_Quark_F2H1[9][1][0]  = -0.1e1 / vh * C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[9][1][1]  = II / vh * C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[9][2][0]  = -0.1e1 / vh * C_MassTop * conj(V31);
  Curvature_Quark_F2H1[9][2][1]  = II / vh * C_MassTop * conj(V31);
  Curvature_Quark_F2H1[9][3][2]  = II / vh * C_MassDown;
  Curvature_Quark_F2H1[9][3][3]  = 0.1e1 / vh * C_MassDown;
  Curvature_Quark_F2H1[10][0][0] = -0.1e1 / vh * C_MassUp * conj(V12);
  Curvature_Quark_F2H1[10][0][1] = II / vh * C_MassUp * conj(V12);
  Curvature_Quark_F2H1[10][1][0] = -0.1e1 / vh * C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[10][1][1] = II / vh * C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[10][2][0] = -0.1e1 / vh * C_MassTop * conj(V32);
  Curvature_Quark_F2H1[10][2][1] = II / vh * C_MassTop * conj(V32);
  Curvature_Quark_F2H1[10][4][2] = II / vh * C_MassStrange;
  Curvature_Quark_F2H1[10][4][3] = 0.1e1 / vh * C_MassStrange;
  Curvature_Quark_F2H1[11][0][0] = -0.1e1 / vh * C_MassUp * conj(V13);
  Curvature_Quark_F2H1[11][0][1] = II / vh * C_MassUp * conj(V13);
  Curvature_Quark_F2H1[11][1][0] = -0.1e1 / vh * C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[11][1][1] = II / vh * C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[11][2][0] = -0.1e1 / vh * C_MassTop * conj(V33);
  Curvature_Quark_F2H1[11][2][1] = II / vh * C_MassTop * conj(V33);
  Curvature_Quark_F2H1[11][5][2] = II / vh * C_MassBottom;
  Curvature_Quark_F2H1[11][5][3] = 0.1e1 / vh * C_MassBottom;

  SetCurvatureDone = true;
}

bool Class_CxSM::CalculateDebyeSimplified()
{
  return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool Class_CxSM::CalculateDebyeGaugeSimplified()
{
  /*
   * Use this function if you calculated the Debye corrections to the gauge mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeGauge[NGauge][NGauge]
   */

  return false;
}
double Class_CxSM::VTreeSimplified(const std::vector<double> &v) const
{
  (void)v;

  return 0;
}

double Class_CxSM::VCounterSimplified(const std::vector<double> &v) const
{
  (void)v;
  if (not UseVCounterSimplified) return 0;
  return 0;
}

void Class_CxSM::Debugging(const std::vector<double> &input,
                           std::vector<double> &output) const
{

  (void)input;
  (void)output;
  Logger::Write(LoggingLevel::Debug, std::string("Start ") + __func__);

  bool Debug = true;

  if (not Debug) return;

  write();

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

  Logger::Write(LoggingLevel::Debug, "Start ");

  // Calculate Tree Level Masses
  std::vector<double> LOMasses(NHiggs), NLOMasses(NHiggs);
  LOMasses = HiggsMassesSquared(vevTree, 0);

  MatrixXd TreeMassMatrix(NHiggs, NHiggs), CTMassMatrix(NHiggs, NHiggs),
      CWMassMatrix(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      TreeMassMatrix(i, j) = 0;
      CTMassMatrix(i, j)   = 0;
      CWMassMatrix(i, j)   = 0;
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      TreeMassMatrix(i, j) = Curvature_Higgs_L2[i][j];
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        TreeMassMatrix(i, j) += Curvature_Higgs_L3[i][j][k] * vevTree[k];
        for (std::size_t l = 0; l < NHiggs; l++)
          TreeMassMatrix(i, j) +=
              0.5 * Curvature_Higgs_L4[i][j][k][l] * vevTree[k] * vevTree[l];
      }

      CTMassMatrix(i, j) = Curvature_Higgs_CT_L2[i][j];
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        CTMassMatrix(i, j) += Curvature_Higgs_CT_L3[i][j][k] * vevTree[k];
        for (std::size_t l = 0; l < NHiggs; l++)
          CTMassMatrix(i, j) +=
              0.5 * Curvature_Higgs_CT_L4[i][j][k][l] * vevTree[k] * vevTree[l];
      }
    }
  }

  CWMassMatrix = HesseWeinberg;

  SelfAdjointEigenSolver<MatrixXd> es(
      TreeMassMatrix + CTMassMatrix + CWMassMatrix, EigenvaluesOnly);

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    NLOMasses[i] = es.eigenvalues()[i];
  }

  Logger::Write(LoggingLevel::Debug, "Print LO | NLO Mass^2 ");
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    Logger::Write(LoggingLevel::Debug,
                  std::to_string(LOMasses[i]) + " | " +
                      std::to_string(NLOMasses[i]));
  }

  MatrixXd HesseCT(NHiggs, NHiggs);
  VectorXd NablaCT(NHiggs);

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HesseCT(i, j) = Curvature_Higgs_CT_L2[i][j];
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        HesseCT(i, j) += Curvature_Higgs_CT_L3[i][j][k] * vevTree[k];
        for (std::size_t l = 0; l < NHiggs; l++)
          HesseCT(i, j) +=
              0.5 * Curvature_Higgs_CT_L4[i][j][k][l] * vevTree[k] * vevTree[l];
      }
    }
  }

  MatrixXd Add(NHiggs, NHiggs);
  Add = HesseCT + HesseWeinberg;
  for (std::size_t i = 0; i < NHiggs; i++)
    for (std::size_t j = 0; j < NHiggs; j++)
      if (std::abs(Add(i, j)) <= 1e-5) Add(i, j) = 0;

  std::stringstream ss;
  ss << "Hesse CT  = \n" << HesseCT << "\n\n\n CT + CW = " << Add << std::endl;

  ss << "Hesse CW = \n" << HesseWeinberg << std::endl;
  Logger::Write(LoggingLevel::Debug, ss.str());

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    NablaCT[i] = Curvature_Higgs_CT_L1[i];
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      NablaCT[i] += Curvature_Higgs_CT_L2[i][j] * vevTree[j];
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        NablaCT[i] +=
            0.5 * Curvature_Higgs_CT_L3[i][j][k] * vevTree[j] * vevTree[k];
        for (std::size_t l = 0; l < NHiggs; l++)
        {
          NablaCT[i] += 1.0 / 6.0 * Curvature_Higgs_CT_L4[i][j][k][l] *
                        vevTree[j] * vevTree[k] * vevTree[l];
        }
      }
    }
  }

  VectorXd AddNabla(NHiggs);
  AddNabla = NablaWeinberg + NablaCT;
  for (std::size_t i = 0; i < NHiggs; i++)
    if (std::abs(AddNabla[i]) <= 1e-5) AddNabla[i] = 0;
  ss.clear();
  ss << "NablaCT  = " << NablaCT.transpose() << "\n"
     << "NablaCW = " << NablaWeinberg.transpose() << "\n"
     << "AddNabla = " << AddNabla.transpose() << "\n";

  ss << "dT4 = " << HesseWeinberg(0, 0) * vh - NablaWeinberg(3) << std::endl;
  ss << "ahhsI = "
     << std::sqrt(2) *
            (HesseWeinberg(4, 5) * va * va -
             vs * (HesseWeinberg(5, 5) * va - NablaWeinberg(5))) /
            (vh * vh * vs)
     << std::endl;

  ss << "HT = " << HesseWeinberg(3, 5) * vs - HesseWeinberg(3, 4) * va
     << std::endl
     << "HT rel =  "
     << (HesseWeinberg(3, 5) * vs - HesseWeinberg(3, 4) * va) /
            (HesseWeinberg(3, 5) * vs + HesseWeinberg(3, 4) * va)
     << std::endl;

  ss << "End " << __func__ << std::endl;
  Logger::Write(LoggingLevel::Debug, ss.str());
}
} // namespace Models
} // namespace BSMPT
