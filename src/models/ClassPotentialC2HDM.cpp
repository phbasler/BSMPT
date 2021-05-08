// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/models/ClassPotentialC2HDM.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/utility.h>
using namespace Eigen;

namespace BSMPT
{
namespace Models
{
Class_Potential_C2HDM::Class_Potential_C2HDM()
{
  // TODO Auto-generated constructor stub
  Model         = ModelID::ModelIDs::C2HDM;
  NNeutralHiggs = 4;
  NChargedHiggs = 4;

  nPar   = 9;
  nParCT = 15;

  nVEV = 4;
  if (!IncludeChargeBreakingVEV) nVEV = 3;

  NHiggs  = NNeutralHiggs + NChargedHiggs;
  NLepton = 9;
  NQuarks = 12;
  NGauge  = 4;

  VevOrder.resize(nVEV);
  if (IncludeChargeBreakingVEV)
  {
    VevOrder[0] = 2; // Charge breaking
    VevOrder[1] = 4; // v1
    VevOrder[2] = 6; // v2
    VevOrder[3] = 7; // CP breaking
  }
  else
  {
    VevOrder[0] = 4; // v1
    VevOrder[1] = 6; // v2
    VevOrder[2] = 7; // CP breaking
  }

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = true;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = true;
}

Class_Potential_C2HDM::~Class_Potential_C2HDM()
{
  // TODO Auto-generated destructor stub
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given inputfile
 */
std::vector<std::string> Class_Potential_C2HDM::addLegendCT() const
{
  std::vector<std::string> labels;
  labels.push_back("Dm11sq");
  labels.push_back("Dm22sq");
  labels.push_back("Dim_m12sq");
  labels.push_back("Dre_m12sq");
  labels.push_back("DL1");
  labels.push_back("DL2");
  labels.push_back("DL3");
  labels.push_back("DL4");
  labels.push_back("Dre_L5");
  labels.push_back("Dim_L5");
  labels.push_back("DTCharged");
  labels.push_back("DT1");
  labels.push_back("DT2");
  labels.push_back("DT3");
  labels.push_back("Dim_L6");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * inputfile
 */
std::vector<std::string> Class_Potential_C2HDM::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c");
  labels.push_back("omega_c");
  labels.push_back("omega_c/T_c");
  if (IncludeChargeBreakingVEV) labels.push_back("omega_CB(T_c)");
  labels.push_back("omega_1(T_c)");
  labels.push_back("omega_2(T_c)");
  labels.push_back("omega_CP(T_c)");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * higgs couplings. Use this to complement the legend of the given inputfile
 *
 *
 */
std::vector<std::string> Class_Potential_C2HDM::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;
  particles.resize(NHiggs);
  particles[0] = ("G^+");
  particles[1] = ("G^-");
  particles[2] = ("H^+");
  particles[3] = ("H^-");
  particles[4] = ("G^0");
  if (UseHsmNotationInTripleHiggs)
  {
    particles[5] = ("h_SM");
    particles[6] = ("h_l");
    particles[7] = ("h_H");
  }
  else
  {
    particles[5] = ("h_1");
    particles[6] = ("h_2");
    particles[7] = ("h_3");
  }

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
 * returns a list of labels which tells the user the chronological order of the
 * VEVs. Use this to complement the legend of the given inputfile
 */
std::vector<std::string> Class_Potential_C2HDM::addLegendVEV() const
{
  std::vector<std::string> labels;
  if (IncludeChargeBreakingVEV) labels.push_back("omega_CB");
  labels.push_back("omega_1");
  labels.push_back("omega_2");
  labels.push_back("omega_CP");
  return labels;
}

void Class_Potential_C2HDM::ReadAndSet(const std::string &linestr,
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
      Type = static_cast<int>(tmp);
    else if (k == 2)
      L1 = tmp;
    else if (k == 3)
      L2 = tmp;
    else if (k == 4)
      L3 = tmp;
    else if (k == 5)
      L4 = tmp;
    else if (k == 6)
      RL5 = tmp;
    else if (k == 7)
      IL5 = tmp;
    else if (k == 8)
      RealMMix = tmp;
    else if (k == 9)
      TanBeta = tmp;
  }

  par[0] = L1;
  par[1] = L2;
  par[2] = L3;
  par[3] = L4;
  par[4] = RL5;
  par[5] = IL5;
  par[6] = RealMMix;
  par[7] = TanBeta;
  par[8] = Type;

  set_gen(par);
  return;
}

void Class_Potential_C2HDM::set_gen(const std::vector<double> &p)
{

  //	double *p = (double *)par;
  scale            = C_vev0;
  L1               = p[0];
  L2               = p[1];
  L3               = p[2];
  L4               = p[3];
  RL5              = p[4];
  IL5              = p[5];
  RealMMix         = p[6];
  TanBeta          = p[7];
  beta             = std::atan(TanBeta);
  Type             = static_cast<int>(p[8]);
  C_CosBetaSquared = 1.0 / (1 + TanBeta * TanBeta);
  C_CosBeta        = sqrt(C_CosBetaSquared);
  C_SinBetaSquared = TanBeta * TanBeta * C_CosBetaSquared;
  C_SinBeta        = sqrt(C_SinBetaSquared);

  u1 = RealMMix * TanBeta -
       C_vev0 * C_vev0 * C_SinBetaSquared * (L4 + RL5 + L3) / 0.2e1 -
       C_vev0 * C_vev0 * C_CosBetaSquared * L1 / 0.2e1;
  u2 = RealMMix * 1.0 / TanBeta -
       C_vev0 * C_vev0 * C_CosBetaSquared * (L4 + RL5 + L3) / 0.2e1 -
       C_vev0 * C_vev0 * C_SinBetaSquared * L2 / 0.2e1;
  Iu3 = C_vev0 * C_vev0 * TanBeta * C_CosBetaSquared * IL5 * 0.5;

  double cb = 0;

  if (Type == 1 or Type == 3) // Type I 2HDM or Lepton Specific
  {
    cb = std::sqrt(2) * C_MassBottom / (C_vev0 * C_SinBeta);
  }
  if (Type == 2 or Type == 4) // Type II 2HDM or Flipped
  {
    cb = std::sqrt(2) * C_MassBottom / (C_vev0 * C_CosBeta);
  }
  CTempC1 = 1.0 / 48 *
            (12 * L1 + 8 * L3 + 4 * L4 + 3 * (3 * C_g * C_g + C_gs * C_gs));
  double ct = std::sqrt(2) * C_MassTop / (C_vev0 * C_SinBeta);
  CTempC2   = 1.0 / 48 *
            (12 * L2 + 8 * L3 + 4 * L4 + 3 * (3 * C_g * C_g + C_gs * C_gs) +
             12 * ct * ct);

  if (Type == 1 or Type == 3)
  {
    CTempC2 += 12.0 / 48.0 * cb * cb;
  }
  else
  {
    CTempC1 += 12.0 / 48.0 * cb * cb;
  }
  vevTreeMin.resize(nVEV);
  if (IncludeChargeBreakingVEV)
  {
    vevTreeMin[0] = 0;
    vevTreeMin[1] = C_vev0 * C_CosBeta;
    vevTreeMin[2] = C_vev0 * C_SinBeta;
    vevTreeMin[3] = 0;
  }
  else
  {
    vevTreeMin[0] = C_vev0 * C_CosBeta;
    vevTreeMin[1] = C_vev0 * C_SinBeta;
    vevTreeMin[2] = 0;
  }
  vevTree.resize(NHiggs);
  vevTree = MinimizeOrderVEV(vevTreeMin);
  if (!SetCurvatureDone) SetCurvatureArrays();
}

void Class_Potential_C2HDM::set_CT_Pot_Par(const std::vector<double> &p)
{
  //	double *p = (double *)par;

  Du1CT     = p[0];
  Du2CT     = p[1];
  DIu3CT    = p[2];
  DRu3CT    = p[3];
  DL1CT     = p[4];
  DL2CT     = p[5];
  DL3CT     = p[6];
  DL4CT     = p[7];
  DRL5CT    = p[8];
  DIL5CT    = p[9];
  DTCharged = p[10];
  DT1       = p[11];
  DT2       = p[12];
  DT3       = p[13];
  DIL6CT    = p[14];

  Curvature_Higgs_CT_L1[2] = DTCharged;
  Curvature_Higgs_CT_L1[4] = DT1;
  Curvature_Higgs_CT_L1[6] = DT2;
  Curvature_Higgs_CT_L1[7] = DT3;

  Curvature_Higgs_CT_L2[0][0] = Du1CT;
  Curvature_Higgs_CT_L2[0][2] = -DRu3CT;
  Curvature_Higgs_CT_L2[0][3] = DIu3CT;
  Curvature_Higgs_CT_L2[1][1] = Du1CT;
  Curvature_Higgs_CT_L2[1][2] = -DIu3CT;
  Curvature_Higgs_CT_L2[1][3] = -DRu3CT;
  Curvature_Higgs_CT_L2[2][0] = -DRu3CT;
  Curvature_Higgs_CT_L2[2][1] = -DIu3CT;
  Curvature_Higgs_CT_L2[2][2] = Du2CT;
  Curvature_Higgs_CT_L2[3][0] = DIu3CT;
  Curvature_Higgs_CT_L2[3][1] = -DRu3CT;
  Curvature_Higgs_CT_L2[3][3] = Du2CT;
  Curvature_Higgs_CT_L2[4][4] = Du1CT;
  Curvature_Higgs_CT_L2[4][6] = -DRu3CT;
  Curvature_Higgs_CT_L2[4][7] = DIu3CT;
  Curvature_Higgs_CT_L2[5][5] = Du1CT;
  Curvature_Higgs_CT_L2[5][6] = -DIu3CT;
  Curvature_Higgs_CT_L2[5][7] = -DRu3CT;
  Curvature_Higgs_CT_L2[6][4] = -DRu3CT;
  Curvature_Higgs_CT_L2[6][5] = -DIu3CT;
  Curvature_Higgs_CT_L2[6][6] = Du2CT;
  Curvature_Higgs_CT_L2[7][4] = DIu3CT;
  Curvature_Higgs_CT_L2[7][5] = -DRu3CT;
  Curvature_Higgs_CT_L2[7][7] = Du2CT;

  {
    Curvature_Higgs_CT_L4[0][0][0][0] = 3 * DL1CT;

    Curvature_Higgs_CT_L4[0][0][0][3] = -3 * DIL6CT;

    Curvature_Higgs_CT_L4[0][0][1][1] = DL1CT;

    Curvature_Higgs_CT_L4[0][0][1][2] = DIL6CT;

    Curvature_Higgs_CT_L4[0][0][2][1] = DIL6CT;

    Curvature_Higgs_CT_L4[0][0][2][2] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[0][0][2][3] = -DIL5CT;

    Curvature_Higgs_CT_L4[0][0][3][0] = -3 * DIL6CT;

    Curvature_Higgs_CT_L4[0][0][3][2] = -DIL5CT;

    Curvature_Higgs_CT_L4[0][0][3][3] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[0][0][4][4] = DL1CT;

    Curvature_Higgs_CT_L4[0][0][4][7] = -DIL6CT;

    Curvature_Higgs_CT_L4[0][0][5][5] = DL1CT;

    Curvature_Higgs_CT_L4[0][0][5][6] = DIL6CT;

    Curvature_Higgs_CT_L4[0][0][6][5] = DIL6CT;

    Curvature_Higgs_CT_L4[0][0][6][6] = DL3CT;

    Curvature_Higgs_CT_L4[0][0][7][4] = -DIL6CT;

    Curvature_Higgs_CT_L4[0][0][7][7] = DL3CT;

    Curvature_Higgs_CT_L4[0][1][0][1] = DL1CT;

    Curvature_Higgs_CT_L4[0][1][0][2] = DIL6CT;

    Curvature_Higgs_CT_L4[0][1][1][0] = DL1CT;

    Curvature_Higgs_CT_L4[0][1][1][3] = -DIL6CT;

    Curvature_Higgs_CT_L4[0][1][2][0] = DIL6CT;

    Curvature_Higgs_CT_L4[0][1][2][2] = DIL5CT;

    Curvature_Higgs_CT_L4[0][1][2][3] = DRL5CT;

    Curvature_Higgs_CT_L4[0][1][3][1] = -DIL6CT;

    Curvature_Higgs_CT_L4[0][1][3][2] = DRL5CT;

    Curvature_Higgs_CT_L4[0][1][3][3] = -DIL5CT;

    Curvature_Higgs_CT_L4[0][2][0][1] = DIL6CT;

    Curvature_Higgs_CT_L4[0][2][0][2] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[0][2][0][3] = -DIL5CT;

    Curvature_Higgs_CT_L4[0][2][1][0] = DIL6CT;

    Curvature_Higgs_CT_L4[0][2][1][2] = DIL5CT;

    Curvature_Higgs_CT_L4[0][2][1][3] = DRL5CT;

    Curvature_Higgs_CT_L4[0][2][2][0] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[0][2][2][1] = DIL5CT;

    Curvature_Higgs_CT_L4[0][2][3][0] = -DIL5CT;

    Curvature_Higgs_CT_L4[0][2][3][1] = DRL5CT;

    Curvature_Higgs_CT_L4[0][2][4][6] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][2][4][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][2][5][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][2][5][7] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][2][6][4] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][2][6][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][2][7][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][2][7][5] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][3][0][0] = -3 * DIL6CT;

    Curvature_Higgs_CT_L4[0][3][0][2] = -DIL5CT;

    Curvature_Higgs_CT_L4[0][3][0][3] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[0][3][1][1] = -DIL6CT;

    Curvature_Higgs_CT_L4[0][3][1][2] = DRL5CT;

    Curvature_Higgs_CT_L4[0][3][1][3] = -DIL5CT;

    Curvature_Higgs_CT_L4[0][3][2][0] = -DIL5CT;

    Curvature_Higgs_CT_L4[0][3][2][1] = DRL5CT;

    Curvature_Higgs_CT_L4[0][3][3][0] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[0][3][3][1] = -DIL5CT;

    Curvature_Higgs_CT_L4[0][3][4][4] = -DIL6CT;

    Curvature_Higgs_CT_L4[0][3][4][6] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][3][4][7] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][3][5][5] = -DIL6CT;

    Curvature_Higgs_CT_L4[0][3][5][6] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][3][5][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][3][6][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][3][6][5] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][3][7][4] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][3][7][5] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][4][0][4] = DL1CT;

    Curvature_Higgs_CT_L4[0][4][0][7] = -DIL6CT;

    Curvature_Higgs_CT_L4[0][4][2][6] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][4][2][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][4][3][4] = -DIL6CT;

    Curvature_Higgs_CT_L4[0][4][3][6] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][4][3][7] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][4][4][0] = DL1CT;

    Curvature_Higgs_CT_L4[0][4][4][3] = -DIL6CT;

    Curvature_Higgs_CT_L4[0][4][6][2] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][4][6][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][4][7][0] = -DIL6CT;

    Curvature_Higgs_CT_L4[0][4][7][2] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][4][7][3] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][5][0][5] = DL1CT;

    Curvature_Higgs_CT_L4[0][5][0][6] = DIL6CT;

    Curvature_Higgs_CT_L4[0][5][2][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][5][2][7] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][5][3][5] = -DIL6CT;

    Curvature_Higgs_CT_L4[0][5][3][6] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][5][3][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][5][5][0] = DL1CT;

    Curvature_Higgs_CT_L4[0][5][5][3] = -DIL6CT;

    Curvature_Higgs_CT_L4[0][5][6][0] = DIL6CT;

    Curvature_Higgs_CT_L4[0][5][6][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][5][6][3] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][5][7][2] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][5][7][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][6][0][5] = DIL6CT;

    Curvature_Higgs_CT_L4[0][6][0][6] = DL3CT;

    Curvature_Higgs_CT_L4[0][6][2][4] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][6][2][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][6][3][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][6][3][5] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][6][4][2] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][6][4][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][6][5][0] = DIL6CT;

    Curvature_Higgs_CT_L4[0][6][5][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][6][5][3] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][6][6][0] = DL3CT;

    Curvature_Higgs_CT_L4[0][7][0][4] = -DIL6CT;

    Curvature_Higgs_CT_L4[0][7][0][7] = DL3CT;

    Curvature_Higgs_CT_L4[0][7][2][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][7][2][5] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][7][3][4] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][7][3][5] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][7][4][0] = -DIL6CT;

    Curvature_Higgs_CT_L4[0][7][4][2] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][7][4][3] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][7][5][2] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[0][7][5][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[0][7][7][0] = DL3CT;

    Curvature_Higgs_CT_L4[1][0][0][1] = DL1CT;

    Curvature_Higgs_CT_L4[1][0][0][2] = DIL6CT;

    Curvature_Higgs_CT_L4[1][0][1][0] = DL1CT;

    Curvature_Higgs_CT_L4[1][0][1][3] = -DIL6CT;

    Curvature_Higgs_CT_L4[1][0][2][0] = DIL6CT;

    Curvature_Higgs_CT_L4[1][0][2][2] = DIL5CT;

    Curvature_Higgs_CT_L4[1][0][2][3] = DRL5CT;

    Curvature_Higgs_CT_L4[1][0][3][1] = -DIL6CT;

    Curvature_Higgs_CT_L4[1][0][3][2] = DRL5CT;

    Curvature_Higgs_CT_L4[1][0][3][3] = -DIL5CT;

    Curvature_Higgs_CT_L4[1][1][0][0] = DL1CT;

    Curvature_Higgs_CT_L4[1][1][0][3] = -DIL6CT;

    Curvature_Higgs_CT_L4[1][1][1][1] = 3 * DL1CT;

    Curvature_Higgs_CT_L4[1][1][1][2] = 3 * DIL6CT;

    Curvature_Higgs_CT_L4[1][1][2][1] = 3 * DIL6CT;

    Curvature_Higgs_CT_L4[1][1][2][2] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[1][1][2][3] = DIL5CT;

    Curvature_Higgs_CT_L4[1][1][3][0] = -DIL6CT;

    Curvature_Higgs_CT_L4[1][1][3][2] = DIL5CT;

    Curvature_Higgs_CT_L4[1][1][3][3] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[1][1][4][4] = DL1CT;

    Curvature_Higgs_CT_L4[1][1][4][7] = -DIL6CT;

    Curvature_Higgs_CT_L4[1][1][5][5] = DL1CT;

    Curvature_Higgs_CT_L4[1][1][5][6] = DIL6CT;

    Curvature_Higgs_CT_L4[1][1][6][5] = DIL6CT;

    Curvature_Higgs_CT_L4[1][1][6][6] = DL3CT;

    Curvature_Higgs_CT_L4[1][1][7][4] = -DIL6CT;

    Curvature_Higgs_CT_L4[1][1][7][7] = DL3CT;

    Curvature_Higgs_CT_L4[1][2][0][0] = DIL6CT;

    Curvature_Higgs_CT_L4[1][2][0][2] = DIL5CT;

    Curvature_Higgs_CT_L4[1][2][0][3] = DRL5CT;

    Curvature_Higgs_CT_L4[1][2][1][1] = 3 * DIL6CT;

    Curvature_Higgs_CT_L4[1][2][1][2] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[1][2][1][3] = DIL5CT;

    Curvature_Higgs_CT_L4[1][2][2][0] = DIL5CT;

    Curvature_Higgs_CT_L4[1][2][2][1] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[1][2][3][0] = DRL5CT;

    Curvature_Higgs_CT_L4[1][2][3][1] = DIL5CT;

    Curvature_Higgs_CT_L4[1][2][4][4] = DIL6CT;

    Curvature_Higgs_CT_L4[1][2][4][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][2][4][7] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][2][5][5] = DIL6CT;

    Curvature_Higgs_CT_L4[1][2][5][6] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][2][5][7] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][2][6][4] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][2][6][5] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][2][7][4] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][2][7][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][3][0][1] = -DIL6CT;

    Curvature_Higgs_CT_L4[1][3][0][2] = DRL5CT;

    Curvature_Higgs_CT_L4[1][3][0][3] = -DIL5CT;

    Curvature_Higgs_CT_L4[1][3][1][0] = -DIL6CT;

    Curvature_Higgs_CT_L4[1][3][1][2] = DIL5CT;

    Curvature_Higgs_CT_L4[1][3][1][3] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[1][3][2][0] = DRL5CT;

    Curvature_Higgs_CT_L4[1][3][2][1] = DIL5CT;

    Curvature_Higgs_CT_L4[1][3][3][0] = -DIL5CT;

    Curvature_Higgs_CT_L4[1][3][3][1] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[1][3][4][6] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][3][4][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][3][5][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][3][5][7] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][3][6][4] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][3][6][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][3][7][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][3][7][5] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][4][1][4] = DL1CT;

    Curvature_Higgs_CT_L4[1][4][1][7] = -DIL6CT;

    Curvature_Higgs_CT_L4[1][4][2][4] = DIL6CT;

    Curvature_Higgs_CT_L4[1][4][2][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][4][2][7] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][4][3][6] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][4][3][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][4][4][1] = DL1CT;

    Curvature_Higgs_CT_L4[1][4][4][2] = DIL6CT;

    Curvature_Higgs_CT_L4[1][4][6][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][4][6][3] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][4][7][1] = -DIL6CT;

    Curvature_Higgs_CT_L4[1][4][7][2] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][4][7][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][5][1][5] = DL1CT;

    Curvature_Higgs_CT_L4[1][5][1][6] = DIL6CT;

    Curvature_Higgs_CT_L4[1][5][2][5] = DIL6CT;

    Curvature_Higgs_CT_L4[1][5][2][6] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][5][2][7] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][5][3][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][5][3][7] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][5][5][1] = DL1CT;

    Curvature_Higgs_CT_L4[1][5][5][2] = DIL6CT;

    Curvature_Higgs_CT_L4[1][5][6][1] = DIL6CT;

    Curvature_Higgs_CT_L4[1][5][6][2] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][5][6][3] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][5][7][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][5][7][3] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][6][1][5] = DIL6CT;

    Curvature_Higgs_CT_L4[1][6][1][6] = DL3CT;

    Curvature_Higgs_CT_L4[1][6][2][4] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][6][2][5] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][6][3][4] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][6][3][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][6][4][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][6][4][3] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][6][5][1] = DIL6CT;

    Curvature_Higgs_CT_L4[1][6][5][2] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][6][5][3] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][6][6][1] = DL3CT;

    Curvature_Higgs_CT_L4[1][7][1][4] = -DIL6CT;

    Curvature_Higgs_CT_L4[1][7][1][7] = DL3CT;

    Curvature_Higgs_CT_L4[1][7][2][4] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][7][2][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][7][3][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][7][3][5] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][7][4][1] = -DIL6CT;

    Curvature_Higgs_CT_L4[1][7][4][2] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][7][4][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][7][5][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[1][7][5][3] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[1][7][7][1] = DL3CT;

    Curvature_Higgs_CT_L4[2][0][0][1] = DIL6CT;

    Curvature_Higgs_CT_L4[2][0][0][2] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[2][0][0][3] = -DIL5CT;

    Curvature_Higgs_CT_L4[2][0][1][0] = DIL6CT;

    Curvature_Higgs_CT_L4[2][0][1][2] = DIL5CT;

    Curvature_Higgs_CT_L4[2][0][1][3] = DRL5CT;

    Curvature_Higgs_CT_L4[2][0][2][0] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[2][0][2][1] = DIL5CT;

    Curvature_Higgs_CT_L4[2][0][3][0] = -DIL5CT;

    Curvature_Higgs_CT_L4[2][0][3][1] = DRL5CT;

    Curvature_Higgs_CT_L4[2][0][4][6] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][0][4][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][0][5][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][0][5][7] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][0][6][4] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][0][6][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][0][7][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][0][7][5] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][1][0][0] = DIL6CT;

    Curvature_Higgs_CT_L4[2][1][0][2] = DIL5CT;

    Curvature_Higgs_CT_L4[2][1][0][3] = DRL5CT;

    Curvature_Higgs_CT_L4[2][1][1][1] = 3 * DIL6CT;

    Curvature_Higgs_CT_L4[2][1][1][2] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[2][1][1][3] = DIL5CT;

    Curvature_Higgs_CT_L4[2][1][2][0] = DIL5CT;

    Curvature_Higgs_CT_L4[2][1][2][1] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[2][1][3][0] = DRL5CT;

    Curvature_Higgs_CT_L4[2][1][3][1] = DIL5CT;

    Curvature_Higgs_CT_L4[2][1][4][4] = DIL6CT;

    Curvature_Higgs_CT_L4[2][1][4][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][1][4][7] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][1][5][5] = DIL6CT;

    Curvature_Higgs_CT_L4[2][1][5][6] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][1][5][7] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][1][6][4] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][1][6][5] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][1][7][4] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][1][7][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][2][0][0] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[2][2][0][1] = DIL5CT;

    Curvature_Higgs_CT_L4[2][2][1][0] = DIL5CT;

    Curvature_Higgs_CT_L4[2][2][1][1] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[2][2][2][2] = 3 * DL2CT;

    Curvature_Higgs_CT_L4[2][2][3][3] = DL2CT;

    Curvature_Higgs_CT_L4[2][2][4][4] = DL3CT;

    Curvature_Higgs_CT_L4[2][2][5][5] = DL3CT;

    Curvature_Higgs_CT_L4[2][2][6][6] = DL2CT;

    Curvature_Higgs_CT_L4[2][2][7][7] = DL2CT;

    Curvature_Higgs_CT_L4[2][3][0][0] = -DIL5CT;

    Curvature_Higgs_CT_L4[2][3][0][1] = DRL5CT;

    Curvature_Higgs_CT_L4[2][3][1][0] = DRL5CT;

    Curvature_Higgs_CT_L4[2][3][1][1] = DIL5CT;

    Curvature_Higgs_CT_L4[2][3][2][3] = DL2CT;

    Curvature_Higgs_CT_L4[2][3][3][2] = DL2CT;

    Curvature_Higgs_CT_L4[2][4][0][6] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][4][0][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][4][1][4] = DIL6CT;

    Curvature_Higgs_CT_L4[2][4][1][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][4][1][7] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][4][2][4] = DL3CT;

    Curvature_Higgs_CT_L4[2][4][4][1] = DIL6CT;

    Curvature_Higgs_CT_L4[2][4][4][2] = DL3CT;

    Curvature_Higgs_CT_L4[2][4][6][0] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][4][6][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][4][7][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][4][7][1] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][5][0][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][5][0][7] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][5][1][5] = DIL6CT;

    Curvature_Higgs_CT_L4[2][5][1][6] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][5][1][7] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][5][2][5] = DL3CT;

    Curvature_Higgs_CT_L4[2][5][5][1] = DIL6CT;

    Curvature_Higgs_CT_L4[2][5][5][2] = DL3CT;

    Curvature_Higgs_CT_L4[2][5][6][0] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][5][6][1] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][5][7][0] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][5][7][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][6][0][4] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][6][0][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][6][1][4] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][6][1][5] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][6][2][6] = DL2CT;

    Curvature_Higgs_CT_L4[2][6][4][0] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][6][4][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][6][5][0] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][6][5][1] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][6][6][2] = DL2CT;

    Curvature_Higgs_CT_L4[2][7][0][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][7][0][5] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][7][1][4] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][7][1][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][7][2][7] = DL2CT;

    Curvature_Higgs_CT_L4[2][7][4][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][7][4][1] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][7][5][0] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[2][7][5][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[2][7][7][2] = DL2CT;

    Curvature_Higgs_CT_L4[3][0][0][0] = -3 * DIL6CT;

    Curvature_Higgs_CT_L4[3][0][0][2] = -DIL5CT;

    Curvature_Higgs_CT_L4[3][0][0][3] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[3][0][1][1] = -DIL6CT;

    Curvature_Higgs_CT_L4[3][0][1][2] = DRL5CT;

    Curvature_Higgs_CT_L4[3][0][1][3] = -DIL5CT;

    Curvature_Higgs_CT_L4[3][0][2][0] = -DIL5CT;

    Curvature_Higgs_CT_L4[3][0][2][1] = DRL5CT;

    Curvature_Higgs_CT_L4[3][0][3][0] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[3][0][3][1] = -DIL5CT;

    Curvature_Higgs_CT_L4[3][0][4][4] = -DIL6CT;

    Curvature_Higgs_CT_L4[3][0][4][6] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][0][4][7] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][0][5][5] = -DIL6CT;

    Curvature_Higgs_CT_L4[3][0][5][6] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][0][5][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][0][6][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][0][6][5] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][0][7][4] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][0][7][5] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][1][0][1] = -DIL6CT;

    Curvature_Higgs_CT_L4[3][1][0][2] = DRL5CT;

    Curvature_Higgs_CT_L4[3][1][0][3] = -DIL5CT;

    Curvature_Higgs_CT_L4[3][1][1][0] = -DIL6CT;

    Curvature_Higgs_CT_L4[3][1][1][2] = DIL5CT;

    Curvature_Higgs_CT_L4[3][1][1][3] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[3][1][2][0] = DRL5CT;

    Curvature_Higgs_CT_L4[3][1][2][1] = DIL5CT;

    Curvature_Higgs_CT_L4[3][1][3][0] = -DIL5CT;

    Curvature_Higgs_CT_L4[3][1][3][1] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[3][1][4][6] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][1][4][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][1][5][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][1][5][7] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][1][6][4] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][1][6][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][1][7][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][1][7][5] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][2][0][0] = -DIL5CT;

    Curvature_Higgs_CT_L4[3][2][0][1] = DRL5CT;

    Curvature_Higgs_CT_L4[3][2][1][0] = DRL5CT;

    Curvature_Higgs_CT_L4[3][2][1][1] = DIL5CT;

    Curvature_Higgs_CT_L4[3][2][2][3] = DL2CT;

    Curvature_Higgs_CT_L4[3][2][3][2] = DL2CT;

    Curvature_Higgs_CT_L4[3][3][0][0] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[3][3][0][1] = -DIL5CT;

    Curvature_Higgs_CT_L4[3][3][1][0] = -DIL5CT;

    Curvature_Higgs_CT_L4[3][3][1][1] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[3][3][2][2] = DL2CT;

    Curvature_Higgs_CT_L4[3][3][3][3] = 3 * DL2CT;

    Curvature_Higgs_CT_L4[3][3][4][4] = DL3CT;

    Curvature_Higgs_CT_L4[3][3][5][5] = DL3CT;

    Curvature_Higgs_CT_L4[3][3][6][6] = DL2CT;

    Curvature_Higgs_CT_L4[3][3][7][7] = DL2CT;

    Curvature_Higgs_CT_L4[3][4][0][4] = -DIL6CT;

    Curvature_Higgs_CT_L4[3][4][0][6] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][4][0][7] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][4][1][6] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][4][1][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][4][3][4] = DL3CT;

    Curvature_Higgs_CT_L4[3][4][4][0] = -DIL6CT;

    Curvature_Higgs_CT_L4[3][4][4][3] = DL3CT;

    Curvature_Higgs_CT_L4[3][4][6][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][4][6][1] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][4][7][0] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][4][7][1] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][5][0][5] = -DIL6CT;

    Curvature_Higgs_CT_L4[3][5][0][6] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][5][0][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][5][1][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][5][1][7] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][5][3][5] = DL3CT;

    Curvature_Higgs_CT_L4[3][5][5][0] = -DIL6CT;

    Curvature_Higgs_CT_L4[3][5][5][3] = DL3CT;

    Curvature_Higgs_CT_L4[3][5][6][0] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][5][6][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][5][7][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][5][7][1] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][6][0][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][6][0][5] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][6][1][4] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][6][1][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][6][3][6] = DL2CT;

    Curvature_Higgs_CT_L4[3][6][4][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][6][4][1] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][6][5][0] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][6][5][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][6][6][3] = DL2CT;

    Curvature_Higgs_CT_L4[3][7][0][4] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][7][0][5] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][7][1][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][7][1][5] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][7][3][7] = DL2CT;

    Curvature_Higgs_CT_L4[3][7][4][0] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][7][4][1] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][7][5][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[3][7][5][1] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[3][7][7][3] = DL2CT;

    Curvature_Higgs_CT_L4[4][0][0][4] = DL1CT;

    Curvature_Higgs_CT_L4[4][0][0][7] = -DIL6CT;

    Curvature_Higgs_CT_L4[4][0][2][6] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][0][2][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][0][3][4] = -DIL6CT;

    Curvature_Higgs_CT_L4[4][0][3][6] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][0][3][7] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][0][4][0] = DL1CT;

    Curvature_Higgs_CT_L4[4][0][4][3] = -DIL6CT;

    Curvature_Higgs_CT_L4[4][0][6][2] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][0][6][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][0][7][0] = -DIL6CT;

    Curvature_Higgs_CT_L4[4][0][7][2] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][0][7][3] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][1][1][4] = DL1CT;

    Curvature_Higgs_CT_L4[4][1][1][7] = -DIL6CT;

    Curvature_Higgs_CT_L4[4][1][2][4] = DIL6CT;

    Curvature_Higgs_CT_L4[4][1][2][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][1][2][7] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][1][3][6] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][1][3][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][1][4][1] = DL1CT;

    Curvature_Higgs_CT_L4[4][1][4][2] = DIL6CT;

    Curvature_Higgs_CT_L4[4][1][6][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][1][6][3] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][1][7][1] = -DIL6CT;

    Curvature_Higgs_CT_L4[4][1][7][2] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][1][7][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][2][0][6] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][2][0][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][2][1][4] = DIL6CT;

    Curvature_Higgs_CT_L4[4][2][1][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][2][1][7] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][2][2][4] = DL3CT;

    Curvature_Higgs_CT_L4[4][2][4][1] = DIL6CT;

    Curvature_Higgs_CT_L4[4][2][4][2] = DL3CT;

    Curvature_Higgs_CT_L4[4][2][6][0] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][2][6][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][2][7][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][2][7][1] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][3][0][4] = -DIL6CT;

    Curvature_Higgs_CT_L4[4][3][0][6] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][3][0][7] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][3][1][6] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][3][1][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][3][3][4] = DL3CT;

    Curvature_Higgs_CT_L4[4][3][4][0] = -DIL6CT;

    Curvature_Higgs_CT_L4[4][3][4][3] = DL3CT;

    Curvature_Higgs_CT_L4[4][3][6][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][3][6][1] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][3][7][0] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][3][7][1] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][4][0][0] = DL1CT;

    Curvature_Higgs_CT_L4[4][4][0][3] = -DIL6CT;

    Curvature_Higgs_CT_L4[4][4][1][1] = DL1CT;

    Curvature_Higgs_CT_L4[4][4][1][2] = DIL6CT;

    Curvature_Higgs_CT_L4[4][4][2][1] = DIL6CT;

    Curvature_Higgs_CT_L4[4][4][2][2] = DL3CT;

    Curvature_Higgs_CT_L4[4][4][3][0] = -DIL6CT;

    Curvature_Higgs_CT_L4[4][4][3][3] = DL3CT;

    Curvature_Higgs_CT_L4[4][4][4][4] = 3 * DL1CT;

    Curvature_Higgs_CT_L4[4][4][4][7] = -3 * DIL6CT;

    Curvature_Higgs_CT_L4[4][4][5][5] = DL1CT;

    Curvature_Higgs_CT_L4[4][4][5][6] = DIL6CT;

    Curvature_Higgs_CT_L4[4][4][6][5] = DIL6CT;

    Curvature_Higgs_CT_L4[4][4][6][6] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[4][4][6][7] = -DIL5CT;

    Curvature_Higgs_CT_L4[4][4][7][4] = -3 * DIL6CT;

    Curvature_Higgs_CT_L4[4][4][7][6] = -DIL5CT;

    Curvature_Higgs_CT_L4[4][4][7][7] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[4][5][4][5] = DL1CT;

    Curvature_Higgs_CT_L4[4][5][4][6] = DIL6CT;

    Curvature_Higgs_CT_L4[4][5][5][4] = DL1CT;

    Curvature_Higgs_CT_L4[4][5][5][7] = -DIL6CT;

    Curvature_Higgs_CT_L4[4][5][6][4] = DIL6CT;

    Curvature_Higgs_CT_L4[4][5][6][6] = DIL5CT;

    Curvature_Higgs_CT_L4[4][5][6][7] = DRL5CT;

    Curvature_Higgs_CT_L4[4][5][7][5] = -DIL6CT;

    Curvature_Higgs_CT_L4[4][5][7][6] = DRL5CT;

    Curvature_Higgs_CT_L4[4][5][7][7] = -DIL5CT;

    Curvature_Higgs_CT_L4[4][6][0][2] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][6][0][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][6][1][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][6][1][3] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][6][2][0] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][6][2][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][6][3][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][6][3][1] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][6][4][5] = DIL6CT;

    Curvature_Higgs_CT_L4[4][6][4][6] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[4][6][4][7] = -DIL5CT;

    Curvature_Higgs_CT_L4[4][6][5][4] = DIL6CT;

    Curvature_Higgs_CT_L4[4][6][5][6] = DIL5CT;

    Curvature_Higgs_CT_L4[4][6][5][7] = DRL5CT;

    Curvature_Higgs_CT_L4[4][6][6][4] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[4][6][6][5] = DIL5CT;

    Curvature_Higgs_CT_L4[4][6][7][4] = -DIL5CT;

    Curvature_Higgs_CT_L4[4][6][7][5] = DRL5CT;

    Curvature_Higgs_CT_L4[4][7][0][0] = -DIL6CT;

    Curvature_Higgs_CT_L4[4][7][0][2] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][7][0][3] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][7][1][1] = -DIL6CT;

    Curvature_Higgs_CT_L4[4][7][1][2] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][7][1][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][7][2][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][7][2][1] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][7][3][0] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[4][7][3][1] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[4][7][4][4] = -3 * DIL6CT;

    Curvature_Higgs_CT_L4[4][7][4][6] = -DIL5CT;

    Curvature_Higgs_CT_L4[4][7][4][7] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[4][7][5][5] = -DIL6CT;

    Curvature_Higgs_CT_L4[4][7][5][6] = DRL5CT;

    Curvature_Higgs_CT_L4[4][7][5][7] = -DIL5CT;

    Curvature_Higgs_CT_L4[4][7][6][4] = -DIL5CT;

    Curvature_Higgs_CT_L4[4][7][6][5] = DRL5CT;

    Curvature_Higgs_CT_L4[4][7][7][4] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[4][7][7][5] = -DIL5CT;

    Curvature_Higgs_CT_L4[5][0][0][5] = DL1CT;

    Curvature_Higgs_CT_L4[5][0][0][6] = DIL6CT;

    Curvature_Higgs_CT_L4[5][0][2][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][0][2][7] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][0][3][5] = -DIL6CT;

    Curvature_Higgs_CT_L4[5][0][3][6] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][0][3][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][0][5][0] = DL1CT;

    Curvature_Higgs_CT_L4[5][0][5][3] = -DIL6CT;

    Curvature_Higgs_CT_L4[5][0][6][0] = DIL6CT;

    Curvature_Higgs_CT_L4[5][0][6][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][0][6][3] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][0][7][2] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][0][7][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][1][1][5] = DL1CT;

    Curvature_Higgs_CT_L4[5][1][1][6] = DIL6CT;

    Curvature_Higgs_CT_L4[5][1][2][5] = DIL6CT;

    Curvature_Higgs_CT_L4[5][1][2][6] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][1][2][7] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][1][3][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][1][3][7] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][1][5][1] = DL1CT;

    Curvature_Higgs_CT_L4[5][1][5][2] = DIL6CT;

    Curvature_Higgs_CT_L4[5][1][6][1] = DIL6CT;

    Curvature_Higgs_CT_L4[5][1][6][2] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][1][6][3] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][1][7][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][1][7][3] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][2][0][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][2][0][7] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][2][1][5] = DIL6CT;

    Curvature_Higgs_CT_L4[5][2][1][6] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][2][1][7] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][2][2][5] = DL3CT;

    Curvature_Higgs_CT_L4[5][2][5][1] = DIL6CT;

    Curvature_Higgs_CT_L4[5][2][5][2] = DL3CT;

    Curvature_Higgs_CT_L4[5][2][6][0] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][2][6][1] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][2][7][0] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][2][7][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][3][0][5] = -DIL6CT;

    Curvature_Higgs_CT_L4[5][3][0][6] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][3][0][7] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][3][1][6] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][3][1][7] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][3][3][5] = DL3CT;

    Curvature_Higgs_CT_L4[5][3][5][0] = -DIL6CT;

    Curvature_Higgs_CT_L4[5][3][5][3] = DL3CT;

    Curvature_Higgs_CT_L4[5][3][6][0] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][3][6][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][3][7][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][3][7][1] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][4][4][5] = DL1CT;

    Curvature_Higgs_CT_L4[5][4][4][6] = DIL6CT;

    Curvature_Higgs_CT_L4[5][4][5][4] = DL1CT;

    Curvature_Higgs_CT_L4[5][4][5][7] = -DIL6CT;

    Curvature_Higgs_CT_L4[5][4][6][4] = DIL6CT;

    Curvature_Higgs_CT_L4[5][4][6][6] = DIL5CT;

    Curvature_Higgs_CT_L4[5][4][6][7] = DRL5CT;

    Curvature_Higgs_CT_L4[5][4][7][5] = -DIL6CT;

    Curvature_Higgs_CT_L4[5][4][7][6] = DRL5CT;

    Curvature_Higgs_CT_L4[5][4][7][7] = -DIL5CT;

    Curvature_Higgs_CT_L4[5][5][0][0] = DL1CT;

    Curvature_Higgs_CT_L4[5][5][0][3] = -DIL6CT;

    Curvature_Higgs_CT_L4[5][5][1][1] = DL1CT;

    Curvature_Higgs_CT_L4[5][5][1][2] = DIL6CT;

    Curvature_Higgs_CT_L4[5][5][2][1] = DIL6CT;

    Curvature_Higgs_CT_L4[5][5][2][2] = DL3CT;

    Curvature_Higgs_CT_L4[5][5][3][0] = -DIL6CT;

    Curvature_Higgs_CT_L4[5][5][3][3] = DL3CT;

    Curvature_Higgs_CT_L4[5][5][4][4] = DL1CT;

    Curvature_Higgs_CT_L4[5][5][4][7] = -DIL6CT;

    Curvature_Higgs_CT_L4[5][5][5][5] = 3 * DL1CT;

    Curvature_Higgs_CT_L4[5][5][5][6] = 3 * DIL6CT;

    Curvature_Higgs_CT_L4[5][5][6][5] = 3 * DIL6CT;

    Curvature_Higgs_CT_L4[5][5][6][6] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[5][5][6][7] = DIL5CT;

    Curvature_Higgs_CT_L4[5][5][7][4] = -DIL6CT;

    Curvature_Higgs_CT_L4[5][5][7][6] = DIL5CT;

    Curvature_Higgs_CT_L4[5][5][7][7] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[5][6][0][0] = DIL6CT;

    Curvature_Higgs_CT_L4[5][6][0][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][6][0][3] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][6][1][1] = DIL6CT;

    Curvature_Higgs_CT_L4[5][6][1][2] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][6][1][3] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][6][2][0] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][6][2][1] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][6][3][0] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][6][3][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][6][4][4] = DIL6CT;

    Curvature_Higgs_CT_L4[5][6][4][6] = DIL5CT;

    Curvature_Higgs_CT_L4[5][6][4][7] = DRL5CT;

    Curvature_Higgs_CT_L4[5][6][5][5] = 3 * DIL6CT;

    Curvature_Higgs_CT_L4[5][6][5][6] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[5][6][5][7] = DIL5CT;

    Curvature_Higgs_CT_L4[5][6][6][4] = DIL5CT;

    Curvature_Higgs_CT_L4[5][6][6][5] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[5][6][7][4] = DRL5CT;

    Curvature_Higgs_CT_L4[5][6][7][5] = DIL5CT;

    Curvature_Higgs_CT_L4[5][7][0][2] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][7][0][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][7][1][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][7][1][3] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][7][2][0] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][7][2][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][7][3][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[5][7][3][1] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[5][7][4][5] = -DIL6CT;

    Curvature_Higgs_CT_L4[5][7][4][6] = DRL5CT;

    Curvature_Higgs_CT_L4[5][7][4][7] = -DIL5CT;

    Curvature_Higgs_CT_L4[5][7][5][4] = -DIL6CT;

    Curvature_Higgs_CT_L4[5][7][5][6] = DIL5CT;

    Curvature_Higgs_CT_L4[5][7][5][7] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[5][7][6][4] = DRL5CT;

    Curvature_Higgs_CT_L4[5][7][6][5] = DIL5CT;

    Curvature_Higgs_CT_L4[5][7][7][4] = -DIL5CT;

    Curvature_Higgs_CT_L4[5][7][7][5] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[6][0][0][5] = DIL6CT;

    Curvature_Higgs_CT_L4[6][0][0][6] = DL3CT;

    Curvature_Higgs_CT_L4[6][0][2][4] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][0][2][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][0][3][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][0][3][5] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][0][4][2] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][0][4][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][0][5][0] = DIL6CT;

    Curvature_Higgs_CT_L4[6][0][5][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][0][5][3] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][0][6][0] = DL3CT;

    Curvature_Higgs_CT_L4[6][1][1][5] = DIL6CT;

    Curvature_Higgs_CT_L4[6][1][1][6] = DL3CT;

    Curvature_Higgs_CT_L4[6][1][2][4] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][1][2][5] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][1][3][4] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][1][3][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][1][4][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][1][4][3] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][1][5][1] = DIL6CT;

    Curvature_Higgs_CT_L4[6][1][5][2] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][1][5][3] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][1][6][1] = DL3CT;

    Curvature_Higgs_CT_L4[6][2][0][4] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][2][0][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][2][1][4] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][2][1][5] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][2][2][6] = DL2CT;

    Curvature_Higgs_CT_L4[6][2][4][0] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][2][4][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][2][5][0] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][2][5][1] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][2][6][2] = DL2CT;

    Curvature_Higgs_CT_L4[6][3][0][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][3][0][5] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][3][1][4] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][3][1][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][3][3][6] = DL2CT;

    Curvature_Higgs_CT_L4[6][3][4][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][3][4][1] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][3][5][0] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][3][5][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][3][6][3] = DL2CT;

    Curvature_Higgs_CT_L4[6][4][0][2] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][4][0][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][4][1][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][4][1][3] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][4][2][0] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][4][2][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][4][3][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][4][3][1] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][4][4][5] = DIL6CT;

    Curvature_Higgs_CT_L4[6][4][4][6] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[6][4][4][7] = -DIL5CT;

    Curvature_Higgs_CT_L4[6][4][5][4] = DIL6CT;

    Curvature_Higgs_CT_L4[6][4][5][6] = DIL5CT;

    Curvature_Higgs_CT_L4[6][4][5][7] = DRL5CT;

    Curvature_Higgs_CT_L4[6][4][6][4] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[6][4][6][5] = DIL5CT;

    Curvature_Higgs_CT_L4[6][4][7][4] = -DIL5CT;

    Curvature_Higgs_CT_L4[6][4][7][5] = DRL5CT;

    Curvature_Higgs_CT_L4[6][5][0][0] = DIL6CT;

    Curvature_Higgs_CT_L4[6][5][0][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][5][0][3] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][5][1][1] = DIL6CT;

    Curvature_Higgs_CT_L4[6][5][1][2] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][5][1][3] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][5][2][0] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][5][2][1] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][5][3][0] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[6][5][3][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[6][5][4][4] = DIL6CT;

    Curvature_Higgs_CT_L4[6][5][4][6] = DIL5CT;

    Curvature_Higgs_CT_L4[6][5][4][7] = DRL5CT;

    Curvature_Higgs_CT_L4[6][5][5][5] = 3 * DIL6CT;

    Curvature_Higgs_CT_L4[6][5][5][6] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[6][5][5][7] = DIL5CT;

    Curvature_Higgs_CT_L4[6][5][6][4] = DIL5CT;

    Curvature_Higgs_CT_L4[6][5][6][5] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[6][5][7][4] = DRL5CT;

    Curvature_Higgs_CT_L4[6][5][7][5] = DIL5CT;

    Curvature_Higgs_CT_L4[6][6][0][0] = DL3CT;

    Curvature_Higgs_CT_L4[6][6][1][1] = DL3CT;

    Curvature_Higgs_CT_L4[6][6][2][2] = DL2CT;

    Curvature_Higgs_CT_L4[6][6][3][3] = DL2CT;

    Curvature_Higgs_CT_L4[6][6][4][4] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[6][6][4][5] = DIL5CT;

    Curvature_Higgs_CT_L4[6][6][5][4] = DIL5CT;

    Curvature_Higgs_CT_L4[6][6][5][5] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[6][6][6][6] = 3 * DL2CT;

    Curvature_Higgs_CT_L4[6][6][7][7] = DL2CT;

    Curvature_Higgs_CT_L4[6][7][4][4] = -DIL5CT;

    Curvature_Higgs_CT_L4[6][7][4][5] = DRL5CT;

    Curvature_Higgs_CT_L4[6][7][5][4] = DRL5CT;

    Curvature_Higgs_CT_L4[6][7][5][5] = DIL5CT;

    Curvature_Higgs_CT_L4[6][7][6][7] = DL2CT;

    Curvature_Higgs_CT_L4[6][7][7][6] = DL2CT;

    Curvature_Higgs_CT_L4[7][0][0][4] = -DIL6CT;

    Curvature_Higgs_CT_L4[7][0][0][7] = DL3CT;

    Curvature_Higgs_CT_L4[7][0][2][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][0][2][5] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][0][3][4] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][0][3][5] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][0][4][0] = -DIL6CT;

    Curvature_Higgs_CT_L4[7][0][4][2] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][0][4][3] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][0][5][2] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][0][5][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][0][7][0] = DL3CT;

    Curvature_Higgs_CT_L4[7][1][1][4] = -DIL6CT;

    Curvature_Higgs_CT_L4[7][1][1][7] = DL3CT;

    Curvature_Higgs_CT_L4[7][1][2][4] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][1][2][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][1][3][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][1][3][5] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][1][4][1] = -DIL6CT;

    Curvature_Higgs_CT_L4[7][1][4][2] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][1][4][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][1][5][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][1][5][3] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][1][7][1] = DL3CT;

    Curvature_Higgs_CT_L4[7][2][0][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][2][0][5] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][2][1][4] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][2][1][5] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][2][2][7] = DL2CT;

    Curvature_Higgs_CT_L4[7][2][4][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][2][4][1] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][2][5][0] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][2][5][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][2][7][2] = DL2CT;

    Curvature_Higgs_CT_L4[7][3][0][4] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][3][0][5] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][3][1][4] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][3][1][5] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][3][3][7] = DL2CT;

    Curvature_Higgs_CT_L4[7][3][4][0] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][3][4][1] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][3][5][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][3][5][1] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][3][7][3] = DL2CT;

    Curvature_Higgs_CT_L4[7][4][0][0] = -DIL6CT;

    Curvature_Higgs_CT_L4[7][4][0][2] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][4][0][3] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][4][1][1] = -DIL6CT;

    Curvature_Higgs_CT_L4[7][4][1][2] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][4][1][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][4][2][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][4][2][1] = -DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][4][3][0] = DL4CT / 2. - DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][4][3][1] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][4][4][4] = -3 * DIL6CT;

    Curvature_Higgs_CT_L4[7][4][4][6] = -DIL5CT;

    Curvature_Higgs_CT_L4[7][4][4][7] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[7][4][5][5] = -DIL6CT;

    Curvature_Higgs_CT_L4[7][4][5][6] = DRL5CT;

    Curvature_Higgs_CT_L4[7][4][5][7] = -DIL5CT;

    Curvature_Higgs_CT_L4[7][4][6][4] = -DIL5CT;

    Curvature_Higgs_CT_L4[7][4][6][5] = DRL5CT;

    Curvature_Higgs_CT_L4[7][4][7][4] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[7][4][7][5] = -DIL5CT;

    Curvature_Higgs_CT_L4[7][5][0][2] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][5][0][3] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][5][1][2] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][5][1][3] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][5][2][0] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][5][2][1] = DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][5][3][0] = -DIL5CT / 2.;

    Curvature_Higgs_CT_L4[7][5][3][1] = DL4CT / 2. + DRL5CT / 2.;

    Curvature_Higgs_CT_L4[7][5][4][5] = -DIL6CT;

    Curvature_Higgs_CT_L4[7][5][4][6] = DRL5CT;

    Curvature_Higgs_CT_L4[7][5][4][7] = -DIL5CT;

    Curvature_Higgs_CT_L4[7][5][5][4] = -DIL6CT;

    Curvature_Higgs_CT_L4[7][5][5][6] = DIL5CT;

    Curvature_Higgs_CT_L4[7][5][5][7] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[7][5][6][4] = DRL5CT;

    Curvature_Higgs_CT_L4[7][5][6][5] = DIL5CT;

    Curvature_Higgs_CT_L4[7][5][7][4] = -DIL5CT;

    Curvature_Higgs_CT_L4[7][5][7][5] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[7][6][4][4] = -DIL5CT;

    Curvature_Higgs_CT_L4[7][6][4][5] = DRL5CT;

    Curvature_Higgs_CT_L4[7][6][5][4] = DRL5CT;

    Curvature_Higgs_CT_L4[7][6][5][5] = DIL5CT;

    Curvature_Higgs_CT_L4[7][6][6][7] = DL2CT;

    Curvature_Higgs_CT_L4[7][6][7][6] = DL2CT;

    Curvature_Higgs_CT_L4[7][7][0][0] = DL3CT;

    Curvature_Higgs_CT_L4[7][7][1][1] = DL3CT;

    Curvature_Higgs_CT_L4[7][7][2][2] = DL2CT;

    Curvature_Higgs_CT_L4[7][7][3][3] = DL2CT;

    Curvature_Higgs_CT_L4[7][7][4][4] = DL3CT + DL4CT - DRL5CT;

    Curvature_Higgs_CT_L4[7][7][4][5] = -DIL5CT;

    Curvature_Higgs_CT_L4[7][7][5][4] = -DIL5CT;

    Curvature_Higgs_CT_L4[7][7][5][5] = DL3CT + DL4CT + DRL5CT;

    Curvature_Higgs_CT_L4[7][7][6][6] = DL2CT;

    Curvature_Higgs_CT_L4[7][7][7][7] = 3 * DL2CT;
  }

  return;
}

/**
 * Console-Output of all Parameters
 */
void Class_Potential_C2HDM::write() const
{
  typedef std::numeric_limits<double> dbl;
  std::cout.precision(dbl::max_digits10);

  double MSM = 0, MhUp = 0, MhDown = 0;

  std::cout << "scale = " << scale << std::endl;

  std::cout << "The parameters are :  \n";
  std::cout << "Model = " << Model << "\n";
  std::cout << "v1 = " << C_vev0 * C_CosBeta << "\n";
  std::cout << "v2 = " << C_vev0 * C_SinBeta << "\n";
  std::cout << "Type = " << Type << "\n";

  std::cout << "beta = " << beta << std::endl;
  std::cout << "tan(beta) = " << TanBeta << std::endl;
  std::cout << "Lambda1 = " << L1 << std::endl;
  std::cout << "Lambda2 = " << L2 << std::endl;
  std::cout << "Lambda3 = " << L3 << std::endl;
  std::cout << "Lambda4 = " << L4 << std::endl;
  std::cout << "Re(Lambda5) = " << RL5 << std::endl;
  std::cout << "Im(Lambda5) = " << IL5 << std::endl;
  std::cout << "Re(m_12^2) = " << RealMMix << std::endl;
  std::cout << "m_{11}^2 = " << u1 << std::endl;
  std::cout << "m_{22}^2 = " << u2 << std::endl;
  std::cout << "Im(m_{12}^2) = " << Iu3 << std::endl;

  std::cout << "The counterterms are :\n";

  std::cout << "DL1 := " << DL1CT << ";\n";
  std::cout << "DL2 := " << DL2CT << ";\n";
  std::cout << "DL3 := " << DL3CT << ";\n";
  std::cout << "DL4 := " << DL4CT << ";\n";
  std::cout << "DRL5 := " << DRL5CT << ";\n";
  std::cout << "DIL5 := " << DIL5CT << ";\n";
  std::cout << "Du1 := " << Du1CT << ";\n";
  std::cout << "Du2 := " << Du2CT << ";\n";
  std::cout << "DRu3 := " << DRu3CT << ";\n";
  std::cout << "DIu3 := " << DIu3CT << ";\n";
  std::cout << "DT1 := " << DT1 << ";\n";
  std::cout << "DT2 := " << DT2 << ";\n";
  std::cout << "DT3:= " << DT3 << ";\n";
  std::cout << "DIL6:= " << DIL6CT << ";\n";

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
    int countposN  = 0;
    int posG0      = 0;
    double testsum = 0;
    for (int i = 0; i < 3; i++)
    {
      //    			testsum = std::abs(HiggsRot(i,0)) + std::abs(HiggsRot(i,2));
      //    			if(testsum != 0) posG1 = i;
      //    			testsum = std::abs(HiggsRot(i,1)) + std::abs(HiggsRot(i,3));
      //    			if(testsum != 0) posG2 = i;
      testsum = std::abs(HiggsRot(i, 5)) + std::abs(HiggsRot(i, 7));
      if (testsum != 0) posG0 = i;
    }
    for (std::size_t i = 3; i < NHiggs; i++)
    {
      testsum = std::abs(HiggsRot(i, 0)) + std::abs(HiggsRot(i, 2));
      if (testsum != 0) posMHCS1 = i;
      //    			testsum = std::abs(HiggsRot(i,1)) + std::abs(HiggsRot(i,3));
      //    			if(testsum != 0) posMHCS2 = i;
      testsum = 0;
      for (int k = 4; k < 8; k++)
        testsum += std::abs(HiggsRot(i, k));
      if (testsum != 0)
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
      // std::cout << NeutralHiggs[i] << "\t" << std::sqrt(NeutralHiggs[i]) <<
      // std::endl;
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

    MatrixXd MassMixing(3, 3);
    for (int i = 0; i < 3; i++)
    {
      MassMixing(i, 0) = NeutralMatrix(i + 1, 0);
      MassMixing(i, 1) = NeutralMatrix(i + 1, 2);
      MassMixing(i, 2) = -std::sin(beta) * NeutralMatrix(i + 1, 1) +
                         std::cos(beta) * NeutralMatrix(i + 1, 3);
    }

    std::cout << "The mass spectrum is given by :\n";
    std::cout << "m_{H^+} = " << std::sqrt(HiggsMasses[posMHCS1]) << " GeV \n"
              << "m_{H_SM} = " << MSM << " GeV \n"
              << "m_{H_l} = " << MhDown << " GeV \n"
              << "m_{H_h} = " << MhUp << " GeV \n";
    std::cout << "The neutral mixing Matrix is given by :\n";
    bool IsNegative = MassMixing(0, 1) < 0;
    std::cout << "H_{SM} = " << MassMixing(0, 0) << " zeta_1 ";
    if (IsNegative)
      std::cout << "-";
    else
      std::cout << "+";
    std::cout << std::abs(MassMixing(0, 1)) << " zeta_2 ";
    IsNegative = MassMixing(0, 2) < 0;
    if (IsNegative)
      std::cout << "-";
    else
      std::cout << "+";
    std::cout << std::abs(MassMixing(0, 2)) << " zeta_3 \n"
              << "H_{l} = " << MassMixing(1, 0) << " zeta_1 ";
    IsNegative = MassMixing(1, 1) < 0;
    if (IsNegative)
      std::cout << "-";
    else
      std::cout << "+";
    std::cout << std::abs(MassMixing(1, 1)) << " zeta_2 ";
    IsNegative = MassMixing(1, 2) < 0;
    if (IsNegative)
      std::cout << "-";
    else
      std::cout << "+";
    std::cout << std::abs(MassMixing(1, 2)) << " zeta_3 \n"
              << "H_{h} = " << MassMixing(2, 0) << " zeta_1 ";
    IsNegative = MassMixing(2, 1) < 0;
    if (IsNegative)
      std::cout << "-";
    else
      std::cout << "+";
    std::cout << std::abs(MassMixing(2, 1)) << " zeta_2 ";
    IsNegative = MassMixing(2, 2) < 0;
    if (IsNegative)
      std::cout << "-";
    else
      std::cout << "+";
    std::cout << std::abs(MassMixing(2, 2)) << " zeta_3 \n";
  }
}

/**
 * Calculates the counterterms in the 2HDM (CP violating as well as CP
 * conserving)
 */
std::vector<double> Class_Potential_C2HDM::calc_CT() const
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

  double v1 = C_vev0 * C_CosBeta;
  double v2 = C_vev0 * C_SinBeta;

  VectorXd NablaWeinberg(8);
  MatrixXd HesseWeinberg(8, 8), HiggsRot(8, 8);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    NablaWeinberg[i] = WeinbergNabla[i];
    for (std::size_t j = 0; j < NHiggs; j++)
      HesseWeinberg(i, j) = WeinbergHesse.at(j * NHiggs + i);
  }

  double freepar =
      -(2 * HesseWeinberg(3, 3)) / std::pow(v1, 2) +
      (2 * HesseWeinberg(7, 7)) /
          std::pow(v1,
                   2); // Free CT parameter is chosen such that DL4CT is zero
  // Du1CT
  parCT.push_back(freepar * std::pow(v2, 2) + HesseWeinberg(4, 4) / 2. +
                  (v2 * HesseWeinberg(4, 6)) / (2. * v1) -
                  (3 * HesseWeinberg(5, 5)) / 2. -
                  (v2 * HesseWeinberg(5, 7)) / (2. * v1));
  // Du2CT
  parCT.push_back(freepar * std::pow(v1, 2) +
                  (v1 * HesseWeinberg(4, 6)) / (2. * v2) -
                  (v1 * HesseWeinberg(5, 7)) / (2. * v2) +
                  HesseWeinberg(6, 6) / 2. - (3 * HesseWeinberg(7, 7)) / 2.);
  // DIu3CT
  parCT.push_back(HesseWeinberg(4, 7) / 2. +
                  (v2 * HesseWeinberg(6, 7)) / (2. * v1) +
                  (3 * NablaWeinberg(5)) / (2. * v2));
  // DRu3CT
  parCT.push_back(freepar * v1 * v2 + HesseWeinberg(5, 7));
  // DL1CT
  parCT.push_back(-((freepar * std::pow(v2, 2)) / std::pow(v1, 2)) -
                  HesseWeinberg(4, 4) / std::pow(v1, 2) +
                  HesseWeinberg(5, 5) / std::pow(v1, 2));
  // DL2CT
  parCT.push_back(-((freepar * std::pow(v1, 2)) / std::pow(v2, 2)) -
                  HesseWeinberg(6, 6) / std::pow(v2, 2) +
                  HesseWeinberg(7, 7) / std::pow(v2, 2));
  // DL3CT
  parCT.push_back(-freepar + HesseWeinberg(1, 3) / (v1 * v2) -
                  HesseWeinberg(3, 3) / std::pow(v1, 2) -
                  HesseWeinberg(4, 6) / (v1 * v2) +
                  HesseWeinberg(7, 7) / std::pow(v1, 2));
  // DL4CT
  parCT.push_back(freepar + (2 * HesseWeinberg(3, 3)) / std::pow(v1, 2) -
                  (2 * HesseWeinberg(7, 7)) / std::pow(v1, 2));
  // DRL5CT
  parCT.push_back(freepar);
  // DIL5CT
  parCT.push_back((2 * HesseWeinberg(6, 7)) / std::pow(v1, 2));
  // DTCharged
  parCT.push_back(-NablaWeinberg(2));
  // DT1
  parCT.push_back(v1 * HesseWeinberg(5, 5) + v2 * HesseWeinberg(5, 7) -
                  NablaWeinberg(4));
  // DT2
  parCT.push_back(v1 * HesseWeinberg(1, 3) + v2 * HesseWeinberg(3, 3) -
                  NablaWeinberg(6));
  // DT3
  parCT.push_back(-((v1 * NablaWeinberg(5)) / v2) - NablaWeinberg(7));
  // CTImLam6
  parCT.push_back(-(HesseWeinberg(4, 5) / (v1 * v2)) -
                  (v2 * HesseWeinberg(6, 7)) / std::pow(v1, 3));

  // std::vector<double> CRelations{

  //     v1 * HCW(0, 3) + v1 * HCW(5, 6) + v2 * HCW(6, 7),

  //     HCW(0, 0) - HCW(1, 1),

  //     v1 * HCW(0, 0) + v2 * HCW(0, 2) - v1 * HCW(5, 5) - v2 * HCW(5, 7),

  //     HCW(0, 3) + HCW(1, 2),

  //     HCW(0, 2) - HCW(1, 3),

  //     v1 * HCW(0, 2) + v2 * HCW(2, 2) - v1 * HCW(5, 7) - v2 * HCW(7, 7),

  //     v2 * HCW(0, 3) - v1 * HCW(4, 5) - v2 * HCW(4, 7),

  //     HCW(2, 2) - HCW(3, 3),

  //     v2 * HCW(0, 3) + NablaWeinberg(5),

  // };

  return parCT;
}

void Class_Potential_C2HDM::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsdone) CalculatePhysicalCouplings();

  if (CalculatedTripleCopulings) return;
  CalculatedTripleCopulings = true;

  std::vector<double> TripleDeriv;
  TripleDeriv = WeinbergThirdDerivative();
  std::vector<std::vector<std::vector<double>>> GaugeBasis(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));
  //    double GaugeBasis[NHiggs][NHiggs][NHiggs];
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

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
    }
  }
  //    std::cout << "HiggsRot = \n" << HiggsRot2 << "\n" << std::endl;
  MatrixXd HiggsRotSort(NHiggs, NHiggs);
  int posMHCS1 = 0, posMHCS2 = 0;
  int posN[3];
  int countposN = 0;
  int posG1 = 0, posG2 = 0, posG0 = 0;
  double testsum = 0;
  for (int i = 0; i < 3; i++)
  {
    testsum = std::abs(HiggsRot(i, 0)) + std::abs(HiggsRot(i, 2));
    if (testsum != 0) posG1 = i;
    testsum = std::abs(HiggsRot(i, 1)) + std::abs(HiggsRot(i, 3));
    if (testsum != 0) posG2 = i;
    testsum = std::abs(HiggsRot(i, 5)) + std::abs(HiggsRot(i, 7));
    if (testsum != 0) posG0 = i;
  }
  for (std::size_t i = 3; i < NHiggs; i++)
  {
    testsum = std::abs(HiggsRot(i, 0)) + std::abs(HiggsRot(i, 2));
    if (testsum != 0) posMHCS1 = i;
    testsum = std::abs(HiggsRot(i, 1)) + std::abs(HiggsRot(i, 3));
    if (testsum != 0) posMHCS2 = i;
    testsum = 0;
    for (int k = 4; k < 8; k++)
      testsum += std::abs(HiggsRot(i, k));
    if (testsum != 0)
    {
      posN[countposN] = i;
      countposN++;
    }
  }

  if (UseHsmNotationInTripleHiggs)
  {
    std::vector<double> HiggsMasses;
    double MhUp = 0, MhDown = 0, MSM = 0;
    HiggsMasses = HiggsMassesSquared(vevTree, 0);

    double NeutralHiggs[3];
    for (int i = 0; i < 3; i++)
    {
      NeutralHiggs[i] = HiggsMasses[posN[i]];
      // std::cout << NeutralHiggs[i] << "\t" << std::sqrt(NeutralHiggs[i]) <<
      // std::endl;
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

  std::vector<double> HiggsOrder(NHiggs);
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

  PosSM = 5;

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
              //  			  double RotFac =
              //  (HiggsRot(i,l)*HiggsRot(j,m)*HiggsRot(k,n)).real();
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

void Class_Potential_C2HDM::SetCurvatureArrays()
{
  initVectors();

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    Curvature_Higgs_L1[i] = 0;
    HiggsVev[i]           = 0;
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
        Curvature_Higgs_L3[i][j][k] = 0;
    }
  }

  HiggsVev[4] = C_vev0 * C_CosBeta;
  HiggsVev[6] = C_vev0 * C_SinBeta;

  Curvature_Higgs_L2[0][0] = u1;
  Curvature_Higgs_L2[0][1] = 0;
  Curvature_Higgs_L2[0][2] = -RealMMix;
  Curvature_Higgs_L2[0][3] = Iu3;
  Curvature_Higgs_L2[0][4] = 0;
  Curvature_Higgs_L2[0][5] = 0;
  Curvature_Higgs_L2[0][6] = 0;
  Curvature_Higgs_L2[0][7] = 0;
  Curvature_Higgs_L2[1][0] = 0;
  Curvature_Higgs_L2[1][1] = u1;
  Curvature_Higgs_L2[1][2] = -Iu3;
  Curvature_Higgs_L2[1][3] = -RealMMix;
  Curvature_Higgs_L2[1][4] = 0;
  Curvature_Higgs_L2[1][5] = 0;
  Curvature_Higgs_L2[1][6] = 0;
  Curvature_Higgs_L2[1][7] = 0;
  Curvature_Higgs_L2[2][0] = -RealMMix;
  Curvature_Higgs_L2[2][1] = -Iu3;
  Curvature_Higgs_L2[2][2] = u2;
  Curvature_Higgs_L2[2][3] = 0;
  Curvature_Higgs_L2[2][4] = 0;
  Curvature_Higgs_L2[2][5] = 0;
  Curvature_Higgs_L2[2][6] = 0;
  Curvature_Higgs_L2[2][7] = 0;
  Curvature_Higgs_L2[3][0] = Iu3;
  Curvature_Higgs_L2[3][1] = -RealMMix;
  Curvature_Higgs_L2[3][2] = 0;
  Curvature_Higgs_L2[3][3] = u2;
  Curvature_Higgs_L2[3][4] = 0;
  Curvature_Higgs_L2[3][5] = 0;
  Curvature_Higgs_L2[3][6] = 0;
  Curvature_Higgs_L2[3][7] = 0;
  Curvature_Higgs_L2[4][0] = 0;
  Curvature_Higgs_L2[4][1] = 0;
  Curvature_Higgs_L2[4][2] = 0;
  Curvature_Higgs_L2[4][3] = 0;
  Curvature_Higgs_L2[4][4] = u1;
  Curvature_Higgs_L2[4][5] = 0;
  Curvature_Higgs_L2[4][6] = -RealMMix;
  Curvature_Higgs_L2[4][7] = Iu3;
  Curvature_Higgs_L2[5][0] = 0;
  Curvature_Higgs_L2[5][1] = 0;
  Curvature_Higgs_L2[5][2] = 0;
  Curvature_Higgs_L2[5][3] = 0;
  Curvature_Higgs_L2[5][4] = 0;
  Curvature_Higgs_L2[5][5] = u1;
  Curvature_Higgs_L2[5][6] = -Iu3;
  Curvature_Higgs_L2[5][7] = -RealMMix;
  Curvature_Higgs_L2[6][0] = 0;
  Curvature_Higgs_L2[6][1] = 0;
  Curvature_Higgs_L2[6][2] = 0;
  Curvature_Higgs_L2[6][3] = 0;
  Curvature_Higgs_L2[6][4] = -RealMMix;
  Curvature_Higgs_L2[6][5] = -Iu3;
  Curvature_Higgs_L2[6][6] = u2;
  Curvature_Higgs_L2[6][7] = 0;
  Curvature_Higgs_L2[7][0] = 0;
  Curvature_Higgs_L2[7][1] = 0;
  Curvature_Higgs_L2[7][2] = 0;
  Curvature_Higgs_L2[7][3] = 0;
  Curvature_Higgs_L2[7][4] = Iu3;
  Curvature_Higgs_L2[7][5] = -RealMMix;
  Curvature_Higgs_L2[7][6] = 0;
  Curvature_Higgs_L2[7][7] = u2;

  Curvature_Higgs_L4[0][0][0][0] = 3 * L1;
  Curvature_Higgs_L4[0][0][1][1] = L1;
  Curvature_Higgs_L4[0][0][2][2] = L3 + L4 + RL5;
  Curvature_Higgs_L4[0][0][2][3] = -IL5;
  Curvature_Higgs_L4[0][0][3][3] = L3 + L4 - RL5;
  Curvature_Higgs_L4[0][0][4][4] = L1;
  Curvature_Higgs_L4[0][0][5][5] = L1;
  Curvature_Higgs_L4[0][0][6][6] = L3;
  Curvature_Higgs_L4[0][0][7][7] = L3;
  Curvature_Higgs_L4[0][1][2][2] = IL5;
  Curvature_Higgs_L4[0][1][2][3] = RL5;
  Curvature_Higgs_L4[0][1][3][3] = -IL5;
  Curvature_Higgs_L4[0][2][4][6] = L4 / 0.2e1 + RL5 / 0.2e1;
  Curvature_Higgs_L4[0][2][4][7] = -IL5 / 0.2e1;
  Curvature_Higgs_L4[0][2][5][6] = IL5 / 0.2e1;
  Curvature_Higgs_L4[0][2][5][7] = L4 / 0.2e1 + RL5 / 0.2e1;
  Curvature_Higgs_L4[0][3][4][6] = -IL5 / 0.2e1;
  Curvature_Higgs_L4[0][3][4][7] = L4 / 0.2e1 - RL5 / 0.2e1;
  Curvature_Higgs_L4[0][3][5][6] = -L4 / 0.2e1 + RL5 / 0.2e1;
  Curvature_Higgs_L4[0][3][5][7] = -IL5 / 0.2e1;
  Curvature_Higgs_L4[1][1][1][1] = 3 * L1;
  Curvature_Higgs_L4[1][1][2][2] = L3 + L4 - RL5;
  Curvature_Higgs_L4[1][1][2][3] = IL5;
  Curvature_Higgs_L4[1][1][3][3] = L3 + L4 + RL5;
  Curvature_Higgs_L4[1][1][4][4] = L1;
  Curvature_Higgs_L4[1][1][5][5] = L1;
  Curvature_Higgs_L4[1][1][6][6] = L3;
  Curvature_Higgs_L4[1][1][7][7] = L3;
  Curvature_Higgs_L4[1][2][4][6] = IL5 / 0.2e1;
  Curvature_Higgs_L4[1][2][4][7] = -L4 / 0.2e1 + RL5 / 0.2e1;
  Curvature_Higgs_L4[1][2][5][6] = L4 / 0.2e1 - RL5 / 0.2e1;
  Curvature_Higgs_L4[1][2][5][7] = IL5 / 0.2e1;
  Curvature_Higgs_L4[1][3][4][6] = L4 / 0.2e1 + RL5 / 0.2e1;
  Curvature_Higgs_L4[1][3][4][7] = -IL5 / 0.2e1;
  Curvature_Higgs_L4[1][3][5][6] = IL5 / 0.2e1;
  Curvature_Higgs_L4[1][3][5][7] = L4 / 0.2e1 + RL5 / 0.2e1;
  Curvature_Higgs_L4[2][2][2][2] = 3 * L2;
  Curvature_Higgs_L4[2][2][3][3] = L2;
  Curvature_Higgs_L4[2][2][4][4] = L3;
  Curvature_Higgs_L4[2][2][5][5] = L3;
  Curvature_Higgs_L4[2][2][6][6] = L2;
  Curvature_Higgs_L4[2][2][7][7] = L2;
  Curvature_Higgs_L4[3][3][3][3] = 3 * L2;
  Curvature_Higgs_L4[3][3][4][4] = L3;
  Curvature_Higgs_L4[3][3][5][5] = L3;
  Curvature_Higgs_L4[3][3][6][6] = L2;
  Curvature_Higgs_L4[3][3][7][7] = L2;
  Curvature_Higgs_L4[4][4][4][4] = 3 * L1;
  Curvature_Higgs_L4[4][4][5][5] = L1;
  Curvature_Higgs_L4[4][4][6][6] = L3 + L4 + RL5;
  Curvature_Higgs_L4[4][4][6][7] = -IL5;
  Curvature_Higgs_L4[4][4][7][7] = L3 + L4 - RL5;
  Curvature_Higgs_L4[4][5][6][6] = IL5;
  Curvature_Higgs_L4[4][5][6][7] = RL5;
  Curvature_Higgs_L4[4][5][7][7] = -IL5;
  Curvature_Higgs_L4[5][5][5][5] = 3 * L1;
  Curvature_Higgs_L4[5][5][6][6] = L3 + L4 - RL5;
  Curvature_Higgs_L4[5][5][6][7] = IL5;
  Curvature_Higgs_L4[5][5][7][7] = L3 + L4 + RL5;
  Curvature_Higgs_L4[6][6][6][6] = 3 * L2;
  Curvature_Higgs_L4[6][6][7][7] = L2;
  Curvature_Higgs_L4[7][7][7][7] = 3 * L2;

  for (std::size_t k1 = 0; k1 < NHiggs; k1++)
  {
    for (std::size_t k2 = k1; k2 < NHiggs; k2++)
    {
      for (std::size_t k3 = k2; k3 < NHiggs; k3++)
      {
        for (std::size_t k4 = k3; k4 < NHiggs; k4++)
        {
          Curvature_Higgs_L4[k2][k3][k4][k1] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k3][k4][k1][k2] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k4][k1][k2][k3] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k2][k1][k3][k4] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k4][k2][k1][k3] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k3][k4][k2][k1] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k1][k3][k4][k2] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k3][k2][k1][k4] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k4][k3][k2][k1] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k1][k4][k3][k2] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k2][k1][k4][k3] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k4][k2][k3][k1] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k1][k4][k2][k3] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k3][k1][k4][k2] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k2][k3][k1][k4] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k1][k3][k2][k4] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k4][k1][k3][k2] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k2][k4][k1][k3] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k3][k2][k4][k1] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k1][k2][k4][k3] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k3][k1][k2][k4] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k4][k3][k1][k2] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
          Curvature_Higgs_L4[k2][k4][k3][k1] =
              Curvature_Higgs_L4[k1][k2][k3][k4];
        }
      }
    }
  }

  for (std::size_t a = 0; a < NGauge; a++)
  {
    for (std::size_t b = 0; b < NGauge; b++)
    {
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        for (std::size_t j = 0; j < NHiggs; j++)
          Curvature_Gauge_G2H2[a][b][i][j] = 0;
      }
    }
  }

  Curvature_Gauge_G2H2[0][0][0][0] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][1][1] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][2][2] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][3][3] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][4][4] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][5][5] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][6][6] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][7][7] = C_g * C_g / 0.2e1;

  Curvature_Gauge_G2H2[0][3][0][4] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][1][5] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][2][6] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][3][7] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][4][0] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][5][1] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][6][2] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][7][3] = C_g * C_gs / 0.2e1;

  Curvature_Gauge_G2H2[1][1][0][0] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][1][1] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][2][2] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][3][3] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][4][4] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][5][5] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][6][6] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][7][7] = C_g * C_g / 0.2e1;

  Curvature_Gauge_G2H2[1][3][0][5] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][1][4] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][2][7] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][3][6] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][4][1] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][5][0] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][6][3] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][7][2] = C_g * C_gs / 0.2e1;

  Curvature_Gauge_G2H2[2][2][0][0] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][1][1] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][2][2] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][3][3] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][4][4] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][5][5] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][6][6] = C_g * C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][7][7] = C_g * C_g / 0.2e1;

  Curvature_Gauge_G2H2[2][3][0][0] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][1][1] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][2][2] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][3][3] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][4][4] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][5][5] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][6][6] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][7][7] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][0][4] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][1][5] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][2][6] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][3][7] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][4][0] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][5][1] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][6][2] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][7][3] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][0][5] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][1][4] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][2][7] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][3][6] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][4][1] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][5][0] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][6][3] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][7][2] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][0][0] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][1][1] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][2][2] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][3][3] = C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][4][4] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][5][5] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][6][6] = -C_g * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][7][7] = -C_g * C_gs / 0.2e1;

  Curvature_Gauge_G2H2[3][3][0][0] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][1][1] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][2][2] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][3][3] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][4][4] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][5][5] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][6][6] = C_gs * C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][7][7] = C_gs * C_gs / 0.2e1;

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

  MatrixXcd YIJR2(NQuarks, NQuarks), YIJE2(NQuarks, NQuarks),
      YIJS2(NQuarks, NQuarks), YIJP2(NQuarks, NQuarks), YIJRD(NQuarks, NQuarks),
      YIJED(NQuarks, NQuarks), YIJSD(NQuarks, NQuarks), YIJPD(NQuarks, NQuarks);
  MatrixXcd YIJRL(NLepton, NLepton), YIJEL(NLepton, NLepton),
      YIJSL(NLepton, NLepton), YIJPL(NLepton, NLepton);
  YIJR2 = MatrixXcd::Zero(NQuarks, NQuarks);
  YIJE2 = MatrixXcd::Zero(NQuarks, NQuarks);
  YIJS2 = MatrixXcd::Zero(NQuarks, NQuarks);
  YIJP2 = MatrixXcd::Zero(NQuarks, NQuarks);
  YIJRD = MatrixXcd::Zero(NQuarks, NQuarks);
  YIJED = MatrixXcd::Zero(NQuarks, NQuarks);
  YIJSD = MatrixXcd::Zero(NQuarks, NQuarks);
  YIJPD = MatrixXcd::Zero(NQuarks, NQuarks);
  YIJRL = MatrixXcd::Zero(NLepton, NLepton);
  YIJEL = MatrixXcd::Zero(NLepton, NLepton);
  YIJSL = MatrixXcd::Zero(NLepton, NLepton);
  YIJPL = MatrixXcd::Zero(NLepton, NLepton);

  std::complex<double> II(0, 1);

  double v1 = C_vev0 * C_CosBeta;
  double v2 = C_vev0 * C_SinBeta;
  double vL = v2;
  double vD = v2;
  if (Type == 2)
  {
    vL = v1;
    vD = v1;
  }
  else if (Type == 3)
    vL = v1;
  else if (Type == 4)
    vD = v1;

  YIJR2(0, 9)  = -std::conj(V11) * C_MassUp / v2;
  YIJR2(0, 10) = -std::conj(V12) * C_MassUp / v2;
  YIJR2(0, 11) = -std::conj(V13) * C_MassUp / v2;

  YIJR2(1, 9)  = -std::conj(V21) * C_MassCharm / v2;
  YIJR2(1, 10) = -std::conj(V22) * C_MassCharm / v2;
  YIJR2(1, 11) = -std::conj(V23) * C_MassCharm / v2;

  YIJR2(2, 9)  = -std::conj(V31) * C_MassTop / v2;
  YIJR2(2, 10) = -std::conj(V32) * C_MassTop / v2;
  YIJR2(2, 11) = -std::conj(V33) * C_MassTop / v2;

  YIJS2(0, 6) = C_MassUp / v2;
  YIJS2(1, 7) = C_MassCharm / v2;
  YIJS2(2, 8) = C_MassTop / v2;

  YIJSD(3, 9)  = C_MassDown / vD;
  YIJSD(4, 10) = C_MassStrange / vD;
  YIJSD(5, 11) = C_MassBottom / vD;

  YIJRD(3, 6) = V11 * C_MassDown / vD;
  YIJRD(3, 7) = V21 * C_MassDown / vD;
  YIJRD(3, 8) = V31 * C_MassDown / vD;
  YIJRD(4, 6) = V12 * C_MassStrange / vD;
  YIJRD(4, 7) = V22 * C_MassStrange / vD;
  YIJRD(4, 8) = V32 * C_MassStrange / vD;
  YIJRD(5, 6) = V13 * C_MassBottom / vD;
  YIJRD(5, 7) = V23 * C_MassBottom / vD;
  YIJRD(5, 8) = V33 * C_MassBottom / vD;

  YIJRL(1, 6) = C_MassElectron / vL;
  YIJRL(3, 7) = C_MassMu / vL;
  YIJRL(5, 8) = C_MassTau / vL;

  YIJSL(0, 1) = C_MassElectron / vL;
  YIJSL(2, 3) = C_MassMu / vL;
  YIJSL(4, 5) = C_MassTau / vL;

  for (std::size_t i = 0; i < NQuarks; i++)
  {
    for (std::size_t j = 0; j < i; j++)
    {
      YIJR2(i, j) = YIJR2(j, i);
      YIJS2(i, j) = YIJS2(j, i);
      YIJRD(i, j) = YIJRD(j, i);
      YIJSD(i, j) = YIJSD(j, i);
    }
  }
  for (std::size_t i = 0; i < NLepton; i++)
  {
    for (std::size_t j = 0; j < i; j++)
    {
      YIJRL(i, j) = YIJRL(j, i);
      YIJSL(i, j) = YIJSL(j, i);
    }
  }

  YIJP2 = std::complex<double>(-1, 0) * II * YIJS2;
  YIJE2 = std::complex<double>(-1, 0) * II * YIJR2;

  YIJPD = II * YIJSD;
  YIJED = II * YIJRD;

  YIJPL = II * YIJSL;
  YIJEL = II * YIJRL;

  for (std::size_t i = 0; i < NQuarks; i++)
  {
    for (std::size_t j = 0; j < NQuarks; j++)
    {
      Curvature_Quark_F2H1[i][j][0] = 0;
      Curvature_Quark_F2H1[i][j][1] = 0;
      Curvature_Quark_F2H1[i][j][2] = YIJR2(i, j);
      Curvature_Quark_F2H1[i][j][3] = YIJE2(i, j);
      Curvature_Quark_F2H1[i][j][4] = 0;
      Curvature_Quark_F2H1[i][j][5] = 0;
      Curvature_Quark_F2H1[i][j][6] = YIJS2(i, j);
      Curvature_Quark_F2H1[i][j][7] = YIJP2(i, j);

      if (Type == 1 or Type == 3)
      {
        Curvature_Quark_F2H1[i][j][2] += YIJRD(i, j);
        Curvature_Quark_F2H1[i][j][3] += YIJED(i, j);
        Curvature_Quark_F2H1[i][j][6] += YIJSD(i, j);
        Curvature_Quark_F2H1[i][j][7] += YIJPD(i, j);
      }
      else
      {
        Curvature_Quark_F2H1[i][j][0] += YIJRD(i, j);
        Curvature_Quark_F2H1[i][j][1] += YIJED(i, j);
        Curvature_Quark_F2H1[i][j][4] += YIJSD(i, j);
        Curvature_Quark_F2H1[i][j][5] += YIJPD(i, j);
      }
    }
  }

  for (std::size_t i = 0; i < NLepton; i++)
  {
    for (std::size_t j = 0; j < NLepton; j++)
    {
      if (Type == 1 or Type == 4)
      {
        Curvature_Lepton_F2H1[i][j][0] = 0;
        Curvature_Lepton_F2H1[i][j][1] = 0;
        Curvature_Lepton_F2H1[i][j][2] = YIJRL(i, j);
        Curvature_Lepton_F2H1[i][j][3] = YIJEL(i, j);
        Curvature_Lepton_F2H1[i][j][4] = 0;
        Curvature_Lepton_F2H1[i][j][5] = 0;
        Curvature_Lepton_F2H1[i][j][6] = YIJSL(i, j);
        Curvature_Lepton_F2H1[i][j][7] = YIJPL(i, j);
      }
      else
      {
        Curvature_Lepton_F2H1[i][j][2] = 0;
        Curvature_Lepton_F2H1[i][j][3] = 0;
        Curvature_Lepton_F2H1[i][j][0] = YIJRL(i, j);
        Curvature_Lepton_F2H1[i][j][1] = YIJEL(i, j);
        Curvature_Lepton_F2H1[i][j][6] = 0;
        Curvature_Lepton_F2H1[i][j][7] = 0;
        Curvature_Lepton_F2H1[i][j][4] = YIJSL(i, j);
        Curvature_Lepton_F2H1[i][j][5] = YIJPL(i, j);
      }
    }
  }

  SetCurvatureDone = true;
}

bool Class_Potential_C2HDM::CalculateDebyeSimplified()
{
  double cb = 0;

  if (Type == 1 or Type == 3) // Type I 2HDM oder Lepton Specific
  {
    cb = std::sqrt(2) * C_MassBottom / (C_vev0 * C_SinBeta);
  }
  if (Type == 2 or Type == 4) // Type II 2HDM oder Flipped
  {
    cb = std::sqrt(2) * C_MassBottom / (C_vev0 * C_CosBeta);
  }
  CTempC1 = 1.0 / 48 *
            (12 * L1 + 8 * L3 + 4 * L4 + 3 * (3 * C_g * C_g + C_gs * C_gs));
  double ct = std::sqrt(2) * C_MassTop / (C_vev0 * C_SinBeta);
  CTempC2   = 1.0 / 48 *
            (12 * L2 + 8 * L3 + 4 * L4 + 3 * (3 * C_g * C_g + C_gs * C_gs) +
             12 * ct * ct);

  if (Type == 1 or Type == 3)
  {
    CTempC2 += 12.0 / 48.0 * cb * cb;
  }
  else
  {
    CTempC1 += 12.0 / 48.0 * cb * cb;
  }

  DebyeHiggs[0][0] = CTempC1;
  DebyeHiggs[1][1] = CTempC1;
  DebyeHiggs[2][2] = CTempC2;
  DebyeHiggs[3][3] = CTempC2;
  DebyeHiggs[4][4] = CTempC1;
  DebyeHiggs[5][5] = CTempC1;
  DebyeHiggs[6][6] = CTempC2;
  DebyeHiggs[7][7] = CTempC2;

  return true;
}

bool Class_Potential_C2HDM::CalculateDebyeGaugeSimplified()
{
  DebyeGauge[0][0] = 2 * C_g * C_g;
  DebyeGauge[1][1] = 2 * C_g * C_g;
  DebyeGauge[2][2] = 2 * C_g * C_g;
  DebyeGauge[3][3] = 2 * C_gs * C_gs;

  return true;
}

double
Class_Potential_C2HDM::VTreeSimplified(const std::vector<double> &v) const
{
  double res = 0;
  if (not UseVTreeSimplified) return 0;

  double v1, v2, vcp, vcb;

  vcb = v[2];
  v1  = v[4];
  v2  = v[6];
  vcp = v[7];

  double C22 = v2 * v2 + vcp * vcp + vcb * vcb;

  res += 0.5 * u1 * v1 * v1;
  res += 0.5 * u2 * C22;
  res += -RealMMix * v1 * v2;
  res += std::pow(v1, 4) / 8.0 * L1;
  res += 1.0 / 8.0 * std::pow(C22, 2) * L2;
  res += L3 / 4.0 * v1 * v1 * C22;
  res += L4 / 4.0 * v1 * v1 * (vcp * vcp + v2 * v2);
  res += RL5 / 4.0 * v1 * v1 * (-vcp * vcp + v2 * v2);

  res += -IL5 / 2.0 * v1 * v1 * v2 * vcp;
  res += Iu3 * v1 * vcp;

  return res;
}

double
Class_Potential_C2HDM::VCounterSimplified(const std::vector<double> &v) const
{
  if (not UseVCounterSimplified) return 0;
  double res = 0;

  double v1, v2, vcp, vcb;

  //	std::cout << "v.size() = " << v.size() << std::endl;

  vcb = v[2];
  v1  = v[4];
  v2  = v[6];
  vcp = v[7];

  //	std::cout << vcb << "\t" << v1 << "\t" << v2 << "\t" << vcp << std::endl;

  double C22 = v2 * v2 + vcp * vcp + vcb * vcb;

  res += 0.5 * Du1CT * v1 * v1;
  res += 0.5 * Du2CT * C22;
  res += -DRu3CT * v1 * v2;
  res += std::pow(v1, 4) / 8.0 * DL1CT;
  res += 1.0 / 8.0 * std::pow(C22, 2) * DL2CT;
  res += DL3CT / 4.0 * v1 * v1 * C22;
  res += DL4CT / 4.0 * v1 * v1 * (vcp * vcp + v2 * v2);
  res += DRL5CT / 4.0 * v1 * v1 * (-vcp * vcp + v2 * v2);
  res += -DIL5CT / 2.0 * v1 * v1 * v2 * vcp;
  res += DIu3CT * v1 * vcp;
  res += DT1 * v1;
  res += DT2 * v2;
  res += DT3 * vcp;
  res += DTCharged * vcb;
  res += -0.5 * DIL6CT * std::pow(v1, 3) * vcp; // Additional CT

  return res;
}

void Class_Potential_C2HDM::Debugging(const std::vector<double> &input,
                                      std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
