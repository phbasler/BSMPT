// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/models/ClassPotentialC2HDM.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>

#include <optional>

using namespace Eigen;

namespace BSMPT
{
namespace Models
{
Class_Potential_C2HDM::Class_Potential_C2HDM(const ISMConstants &smConstants)
    : Class_Potential_Origin(smConstants)
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
  scale            = SMConstants.C_vev0;
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
       SMConstants.C_vev0 * SMConstants.C_vev0 * C_SinBetaSquared *
           (L4 + RL5 + L3) / 0.2e1 -
       SMConstants.C_vev0 * SMConstants.C_vev0 * C_CosBetaSquared * L1 / 0.2e1;
  u2 = RealMMix * 1.0 / TanBeta -
       SMConstants.C_vev0 * SMConstants.C_vev0 * C_CosBetaSquared *
           (L4 + RL5 + L3) / 0.2e1 -
       SMConstants.C_vev0 * SMConstants.C_vev0 * C_SinBetaSquared * L2 / 0.2e1;
  Iu3 = SMConstants.C_vev0 * SMConstants.C_vev0 * TanBeta * C_CosBetaSquared *
        IL5 * 0.5;

  double cb = 0;

  if (Type == 1 or Type == 3) // Type I 2HDM or Lepton Specific
  {
    cb = std::sqrt(2) * SMConstants.C_MassBottom /
         (SMConstants.C_vev0 * C_SinBeta);
  }
  if (Type == 2 or Type == 4) // Type II 2HDM or Flipped
  {
    cb = std::sqrt(2) * SMConstants.C_MassBottom /
         (SMConstants.C_vev0 * C_CosBeta);
  }
  CTempC1 = 1.0 / 48 *
            (12 * L1 + 8 * L3 + 4 * L4 +
             3 * (3 * SMConstants.C_g * SMConstants.C_g +
                  SMConstants.C_gs * SMConstants.C_gs));
  double ct =
      std::sqrt(2) * SMConstants.C_MassTop / (SMConstants.C_vev0 * C_SinBeta);
  CTempC2 = 1.0 / 48 *
            (12 * L2 + 8 * L3 + 4 * L4 +
             3 * (3 * SMConstants.C_g * SMConstants.C_g +
                  SMConstants.C_gs * SMConstants.C_gs) +
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
    vevTreeMin[1] = SMConstants.C_vev0 * C_CosBeta;
    vevTreeMin[2] = SMConstants.C_vev0 * C_SinBeta;
    vevTreeMin[3] = 0;
  }
  else
  {
    vevTreeMin[0] = SMConstants.C_vev0 * C_CosBeta;
    vevTreeMin[1] = SMConstants.C_vev0 * C_SinBeta;
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
  std::stringstream ss;
  ss.precision(std::numeric_limits<double>::max_digits10);

  ss << "scale = " << scale << std::endl;

  ss << "The parameters are :  \n";
  ss << "Model = " << Model << "\n";
  ss << "v1 = " << SMConstants.C_vev0 * C_CosBeta << "\n";
  ss << "v2 = " << SMConstants.C_vev0 * C_SinBeta << "\n";
  ss << "Type = " << Type << "\n";

  ss << "beta = " << beta << std::endl;
  ss << "tan(beta) = " << TanBeta << std::endl;
  ss << "Lambda1 = " << L1 << std::endl;
  ss << "Lambda2 = " << L2 << std::endl;
  ss << "Lambda3 = " << L3 << std::endl;
  ss << "Lambda4 = " << L4 << std::endl;
  ss << "Re(Lambda5) = " << RL5 << std::endl;
  ss << "Im(Lambda5) = " << IL5 << std::endl;
  ss << "Re(m_12^2) = " << RealMMix << std::endl;
  ss << "m_{11}^2 = " << u1 << std::endl;
  ss << "m_{22}^2 = " << u2 << std::endl;
  ss << "Im(m_{12}^2) = " << Iu3 << std::endl;

  ss << "The counterterms are :\n";

  ss << "DL1 := " << DL1CT << ";\n";
  ss << "DL2 := " << DL2CT << ";\n";
  ss << "DL3 := " << DL3CT << ";\n";
  ss << "DL4 := " << DL4CT << ";\n";
  ss << "DRL5 := " << DRL5CT << ";\n";
  ss << "DIL5 := " << DIL5CT << ";\n";
  ss << "Du1 := " << Du1CT << ";\n";
  ss << "Du2 := " << Du2CT << ";\n";
  ss << "DRu3 := " << DRu3CT << ";\n";
  ss << "DIu3 := " << DIu3CT << ";\n";
  ss << "DT1 := " << DT1 << ";\n";
  ss << "DT2 := " << DT2 << ";\n";
  ss << "DT3:= " << DT3 << ";\n";
  ss << "DIL6:= " << DIL6CT << ";\n";

  if (CalcCouplingsDone)
  {
    MatrixXd HiggsRot(NHiggs, NHiggs);
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      for (std::size_t j = 0; j < NHiggs; j++)
      {
        HiggsRot(i, j) = HiggsRotationMatrixEnsuredConvention[i][j];
      }
    }

    std::vector<double> HiggsMasses;
    HiggsMasses = HiggsMassesSquared(vevTree, 0);

    MatrixXd NeutralMatrix(4, 4);
    for (int j = 0; j < 4; j++)
    {
      NeutralMatrix(0, j) = HiggsRot(pos_G0, j + 4);
      NeutralMatrix(1, j) = HiggsRot(pos_h_SM, j + 4);
      NeutralMatrix(2, j) = HiggsRot(pos_h_l, j + 4);
      NeutralMatrix(3, j) = HiggsRot(pos_h_H, j + 4);
    }

    MatrixXd MassMixing(3, 3);
    for (int i = 0; i < 3; i++)
    {
      MassMixing(i, 0) = NeutralMatrix(i + 1, 0);
      MassMixing(i, 1) = NeutralMatrix(i + 1, 2);
      MassMixing(i, 2) = -std::sin(beta) * NeutralMatrix(i + 1, 1) +
                         std::cos(beta) * NeutralMatrix(i + 1, 3);
    }

    ss << "The mass spectrum is given by :\n"
       << "m_{H^+} = " << std::sqrt(HiggsMasses[pos_Hp]) << " GeV \n"
       << "m_{H_SM} = " << std::sqrt(HiggsMasses[pos_h_SM]) << " GeV \n"
       << "m_{H_l} = " << std::sqrt(HiggsMasses[pos_h_l]) << " GeV \n"
       << "m_{H_h} = " << std::sqrt(HiggsMasses[pos_h_H]) << " GeV \n";

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
  if (!CalcCouplingsDone)
  {
    std::string retmes = __func__;
    retmes += " was called before CalculatePhysicalCouplings()!\n";
    throw std::runtime_error(retmes);
  }
  std::vector<double> WeinbergNabla, WeinbergHesse;
  WeinbergNabla = WeinbergFirstDerivative();
  WeinbergHesse = WeinbergSecondDerivative();

  double v1 = SMConstants.C_vev0 * C_CosBeta;
  double v2 = SMConstants.C_vev0 * C_SinBeta;

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

/**
 * Ensures the correct rotation matrix convention
 */
void Class_Potential_C2HDM::AdjustRotationMatrix()
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

  // C2HDM interaction basis
  // rho1, eta1, rho2, eta2, zeta1, psi1, zeta2, psi2
  const std::size_t pos_rho1 = 0, pos_eta1 = 1, pos_rho2 = 2, pos_eta2 = 3,
                    pos_zeta1 = 4, pos_psi1 = 5, pos_zeta2 = 6, pos_psi2 = 7;

  // Indices of mass eigenstates for rotation from interaction to mass basis.
  // Using temporary optional variables, later being set to the instance variables
  // pos_G0, pos_Gp, etc.
  std::optional<std::size_t> tpos_G0, tpos_Gp, tpos_Gm, tpos_Hp, tpos_Hm,
                             tpos_h1, tpos_h2, tpos_h3;

  for (std::size_t i = 0; i < NHiggs; i++)
  // mass base index i corresponds to mass vector sorted in ascending mass
  {
    // Goldstones have zero mass in the Landau gauge
    bool hasZeroMass = std::abs(HiggsMasses[i]) < ZeroThreshold;
    bool hasPosPsi12  = std::abs(HiggsRot(i, pos_psi1))
                        + std::abs(HiggsRot(i, pos_psi2)) > ZeroThreshold;
    bool hasPosZeta12 = std::abs(HiggsRot(i, pos_zeta1))
                        + std::abs(HiggsRot(i, pos_zeta2)) > ZeroThreshold;
    // Charged submatrix.
    // Make use of the fact that there is no mixing between the rho1,2 and
    // the eta1,2 states, otherwise this part with if/else if would not work
    if (std::abs(HiggsRot(i, pos_rho1)) + std::abs(HiggsRot(i, pos_rho2)) >
        ZeroThreshold)
    {
      if (not tpos_Gp.has_value() and hasZeroMass)
      {
        tpos_Gp = i;
      }
      else if (not tpos_Hp.has_value())
      {
        tpos_Hp = i;
      }
      else
      {
        throw std::runtime_error("Error. Charged submatrix Gp/Hp mixing with "
                                 "other components.");
      }
    }
    else if (std::abs(HiggsRot(i, pos_eta1)) + std::abs(HiggsRot(i, pos_eta2)) >
             ZeroThreshold)
    {
      if (not tpos_Gm.has_value() and hasZeroMass)
      {
        tpos_Gm = i;
      }
      else if (not tpos_Hm.has_value())
      {
        tpos_Hm = i;
      }
      else
      {
        throw std::runtime_error("Error. Charged submatrix Gm/Hm mixing with "
                                 "other components.");
      }
    }
    // Neutral CP-mixed submatrix
    else if (hasPosPsi12 or hasPosZeta12)
    {
      // Goldstone is massless and has no zeta1,2 componenta
      if (not tpos_G0.has_value() and hasZeroMass and not hasPosZeta12)
      {
        tpos_G0 = i;
      }
      // use that mh1 < mh2 < mh3
      else if (not tpos_h1.has_value())
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
        throw std::runtime_error("Error. Neutral submatrix mixing with "
                                 "other components.");
      }
    }
    else
    {
      throw std::runtime_error("Error. Invalid mixing matrix containing row "
                               "with all zeroes.");
    }
  }

  // Sanity check if all position indices are set
  if (not (tpos_G0.has_value() and tpos_Gp.has_value() and tpos_Gm.has_value()
           and tpos_Hp.has_value() and tpos_Hp.has_value()
           and tpos_h1.has_value() and tpos_h2.has_value()
           and tpos_h3.has_value())
     )
  {
    throw std::runtime_error("Error. Not all position indices are set.");
  }

  pos_G0 = tpos_G0.value();
  pos_Gp = tpos_Gp.value();
  pos_Gm = tpos_Gm.value();
  pos_Hp = tpos_Hp.value();
  pos_Hm = tpos_Hm.value();
  pos_h1 = tpos_h1.value();
  pos_h2 = tpos_h2.value();
  pos_h3 = tpos_h3.value();

  // check if all other elements of rotation matrix are zero
  bool zero_element = false;
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      if (not((j == pos_rho1 and (i == pos_Gp or i == pos_Hp)) or
              (j == pos_rho2 and (i == pos_Gp or i == pos_Hp)) or
              (j == pos_eta1 and (i == pos_Gm or i == pos_Hm)) or
              (j == pos_eta2 and (i == pos_Gm or i == pos_Hm)) or
              (j == pos_psi1 and (i == pos_G0 or i == pos_h1 or
                                  i == pos_h2 or i == pos_h3)) or
              (j == pos_psi2 and (i == pos_G0 or i == pos_h1 or
                                  i == pos_h2 or i == pos_h3)) or
              (j == pos_zeta1 and
              (i == pos_h1 or i == pos_h2 or i == pos_h3)) or
              (j == pos_zeta2 and
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

  // Determine the additional indices for the SM-like
  // and lighter/heavier Higgses.
  // Due to the masses being ordered, we will always have
  //  HiggsMasses[pos_h1] <= HiggsMasses[pos_h2] <= HiggsMasses[pos_h3]
  double diff1 =
      std::abs(std::sqrt(HiggsMasses[pos_h1]) - SMConstants.C_MassSMHiggs);
  double diff2 =
      std::abs(std::sqrt(HiggsMasses[pos_h2]) - SMConstants.C_MassSMHiggs);
  double diff3 =
      std::abs(std::sqrt(HiggsMasses[pos_h3]) - SMConstants.C_MassSMHiggs);
  if (diff1 < diff2 and diff1 < diff3)
  {
    pos_h_SM = pos_h1;
    pos_h_l  = pos_h2;
    pos_h_H  = pos_h3;
  }
  else if (diff2 < diff1 and diff2 < diff3)
  {
    pos_h_l  = pos_h1;
    pos_h_SM = pos_h2;
    pos_h_H  = pos_h3;
  }
  else
  {
    pos_h_l  = pos_h1;
    pos_h_H  = pos_h2;
    pos_h_SM = pos_h3;
  }

  // Steps:
  // (1) Rotate mass matrix from interaction to semi-interaction basis
  //     (i.e. interaction basis with neutral Goldstone rotated out):
  //
  //     From interaction basis
  //       0     1     2     3     4      5      6      7
  //       rho1, eta1, rho2, eta2, zeta1, psi1,  zeta2, psi2
  //     to semi-interaction basis (Goldstone rotated out)
  //       0     1     2     3     4      5      6      7
  //       G^0,  rho1, eta1, rho2, eta2,  zeta1, zeta2, zeta3
  //
  // (2) Diagonalise mass matrix in semi-interaction basis
  //     -> obtain rotation matrix from semi-interaction to mass basis
  //
  // (3) Neutral part of rotation matrix from step 2 must have the 3x3 form
  //     as in arXiv:1803.02846 Eq. (3.91), so this is the one to check for
  //     R11 > 0, R33 > 0, det > 0 (see arXiv:2007.02985 Eq. (6))

  // Find position of neutral Goldstone:
  // * Mass eigenvalues are ordered from smallest to largest
  //   => First three rows correspond to Goldstones
  //      (= massless in Landau gauge)
  // * Charged and neutral Goldstones do not mix
  //   => Look for row which has psi1 and psi2 mixing components =/= 0,
  //      and the rest = 0

  // Matrix to "rotate out" the neutral Goldstone boson, see arXiv:1803.02846
  // Eq. (3.89)
  MatrixXd RotGoldstone(NHiggs, NHiggs);
  RotGoldstone.row(0) << 0., 0., 0., 0., 0., C_CosBeta, 0., C_SinBeta;
  RotGoldstone.row(1) << 1., 0., 0., 0., 0., 0., 0., 0.;
  RotGoldstone.row(2) << 0., 1., 0., 0., 0., 0., 0., 0.;
  RotGoldstone.row(3) << 0., 0., 1., 0., 0., 0., 0., 0.;
  RotGoldstone.row(4) << 0., 0., 0., 1., 0., 0., 0., 0.;
  RotGoldstone.row(5) << 0., 0., 0., 0., 1., 0., 0., 0.;
  RotGoldstone.row(6) << 0., 0., 0., 0., 0., 0., 1., 0.;
  RotGoldstone.row(7) << 0., 0., 0., 0., 0., -C_SinBeta, 0., C_CosBeta;

  // Swap rows to ensure that G0 is always in the first row
  // for the following manipulations
  MatrixXd MoveGoldstoneFirst(NHiggs, NHiggs);
  MoveGoldstoneFirst.setIdentity(NHiggs, NHiggs);
  if (pos_G0 != 0)
  {
    MoveGoldstoneFirst(0, 0)           = 0.;
    MoveGoldstoneFirst(pos_G0, pos_G0) = 0.;
    MoveGoldstoneFirst(0, pos_G0)      = 1.;
    MoveGoldstoneFirst(pos_G0, 0)      = 1.;
  }

  // Compute rotation matrix from the "semi-interaction" (with G0 rotated out)
  // to the mass basis, to get the same rotation as in arXiv:1803.02846 Eqs.
  // (3.90)-(3.91)
  MatrixXd RotGoldstoneMassBasis(NHiggs, NHiggs);
  RotGoldstoneMassBasis =
      MoveGoldstoneFirst * HiggsRot * RotGoldstone.transpose();

  // Semi-interaction basis (neutral Goldstone rotated out)
  // G^0 == 0 (not used), rho1, eta1, rho2, eta2, zeta1, zeta2, zeta3
  const std::size_t pos_si_G0 = 0, pos_si_rho1 = 1, pos_si_eta1 = 2,
                    pos_si_rho2 = 3, pos_si_eta2 = 4, pos_si_zeta1 = 5,
                    pos_si_zeta2 = 6, pos_si_zeta3 = 7;

  double row1 = 0.0, col1 = 0.0;
  // Sum only over index starting from 1 (i.e. don't include the (0,0) element,
  // i.e. the upper left element, of the matrix); the sum should be very small
  for (std::size_t i = 1; i < NHiggs; i++)
  {
    row1 += std::abs(RotGoldstoneMassBasis(0, i));
    col1 += std::abs(RotGoldstoneMassBasis(i, 0));
  }

  // Consistency check that the Goldstone was rotated out properly:
  // first row/column should contain only zeroes except for the upper left
  // element
  if (std::abs(std::abs(RotGoldstoneMassBasis(0, 0)) - 1.0) > ZeroThreshold or
      std::abs(row1) > ZeroThreshold or std::abs(col1) > ZeroThreshold)
  {
    throw std::runtime_error("Error. Something went wrong after rotating "
                             "out the neutral Goldstone.");
  }

  // Indices of mass eigenstates for rotation from semi-interaction to mass
  // basis; position of neutral Goldstone is fixed to 0, see above
  // int pos_si_G1 = -1, pos_si_G2 = -1, pos_si_H1 = -1, pos_si_H2 = -1;
  // int pos_si_h1 = -1, pos_si_h2 = -1, pos_si_h3 = -1;

  std::optional<std::size_t> tpos_si_Gp, tpos_si_Gm, tpos_si_Hp, tpos_si_Hm,
                             tpos_si_h1, tpos_si_h2, tpos_si_h3;

  // Start with i = 1, i.e. skip over the neutral Goldstone
  for (std::size_t i = 1; i < NHiggs; i++)
  // mass base index i corresponds to mass vector sorted in ascending mass
  {
    bool hasZeroMass = std::abs(HiggsMasses[i]) < ZeroThreshold;
    // Charged submatrices
    // Check if the field with index i has a rho1 or rho2 component;
    // since Goldstone mass is zero (Landau gauge), it appears before the
    // charged Higgs
    if (std::abs(RotGoldstoneMassBasis(i, pos_si_rho1)) +
            std::abs(RotGoldstoneMassBasis(i, pos_si_rho2)) >
        ZeroThreshold)
    // use that 0 = mGpm < mHpm
    {
      if (not tpos_si_Gp.has_value() and hasZeroMass)
      {
        tpos_si_Gp = i;
      }
      else if (not tpos_si_Hp.has_value())
      {
        tpos_si_Hp = i;
      }
      else
      {
        throw std::runtime_error("Error. Charged submatrix Gp/Hp_si mixing "
                                 "with other components.");
      }
    }
    else if (std::abs(RotGoldstoneMassBasis(i, pos_si_eta1)) +
                 std::abs(RotGoldstoneMassBasis(i, pos_si_eta2)) >
             ZeroThreshold)
    // use that 0 = mGpm < mHpm
    {
      if (not tpos_si_Gm.has_value() and hasZeroMass)
      {
        tpos_si_Gm = i;
      }
      else if (not tpos_si_Hm.has_value())
      {
        tpos_si_Hm = i;
      }
      else
      {
        throw std::runtime_error("Error. Charged submatrix Gm/Hm_si mixing "
                                 "with other components.");
      }
    }
    // Neutral submatrix (mixed CP-even and CP-odd states);
    // neutral Goldstone already rotated out
    else if (std::abs(RotGoldstoneMassBasis(i, pos_si_zeta1)) +
                 std::abs(RotGoldstoneMassBasis(i, pos_si_zeta2)) +
                 std::abs(RotGoldstoneMassBasis(i, pos_si_zeta3)) >
             ZeroThreshold)
    // use that mh1 < mh2 < mh3
    {
      if (not tpos_si_h1.has_value())
      {
        tpos_si_h1 = i;
      }
      else if (not tpos_si_h2.has_value())
      {
        tpos_si_h2 = i;
      }
      else if (not tpos_si_h3.has_value())
      {
        tpos_si_h3 = i;
      }
      else
      {
        throw std::runtime_error("Error. Neutral submatrix _si mixing "
                                 "with other components.");
      }
    }
    else
    {
      throw std::runtime_error("Error. Invalid mixing matrix _si containing "
                               "row with all zeroes.");
    }
  }

  // Check if all position indices are set
  if (not (tpos_si_Gp.has_value() and tpos_si_Gm.has_value()
           and tpos_si_Hp.has_value() and tpos_si_Hm.has_value()
           and tpos_si_h1.has_value() and tpos_si_h2.has_value()
           and tpos_si_h3.has_value())
     )
  {
    throw std::runtime_error("Error. Not all position indices are set.");
  }

  std::size_t pos_si_Gp = tpos_si_Gp.value();
  std::size_t pos_si_Gm = tpos_si_Gm.value();
  std::size_t pos_si_Hp = tpos_si_Hp.value();
  std::size_t pos_si_Hm = tpos_si_Hm.value();
  std::size_t pos_si_h1 = tpos_si_h1.value();
  std::size_t pos_si_h2 = tpos_si_h2.value();
  std::size_t pos_si_h3 = tpos_si_h3.value();

  // Check if all other elements of rotation matrix are zero
  zero_element = false;
  // Start with i, j = 1, skip neutral Goldstone
  for (std::size_t i = 1; i < NHiggs; i++)
  {
    for (std::size_t j = 1; j < NHiggs; j++)
    {
      if (not((j == pos_si_rho1 and (i == pos_si_Gp or i == pos_si_Hp)) or
              (j == pos_si_eta1 and (i == pos_si_Gm or i == pos_si_Hm)) or
              (j == pos_si_rho2 and (i == pos_si_Gp or i == pos_si_Hp)) or
              (j == pos_si_eta2 and (i == pos_si_Gm or i == pos_si_Hm)) or
              (j == pos_si_zeta1 and
               (i == pos_si_h1 or i == pos_si_h2 or i == pos_si_h3)) or
              (j == pos_si_zeta2 and
               (i == pos_si_h1 or i == pos_si_h2 or i == pos_si_h3)) or
              (j == pos_si_zeta3 and
               (i == pos_si_h1 or i == pos_si_h2 or i == pos_si_h3))))
      {
        zero_element = true;
      }

      if (zero_element and
          std::abs(RotGoldstoneMassBasis(i, j)) > ZeroThreshold)
      {
        throw std::runtime_error("Error. Invalid rotation matrix detected.");
      }
      zero_element = false;
    }
  }

  MatrixXd HiggsRotFixed(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsRotFixed.row(i) = RotGoldstoneMassBasis.row(i);
  }

  // Neutral Goldstone; flip sign if its element, which should be "1", is
  // negative
  if (HiggsRotFixed(pos_si_G0, pos_si_G0) < 0) // G0 G0 (+1)
  {
    HiggsRotFixed.row(pos_si_G0) *= -1;
  }

  // charged submatrix
  if (HiggsRotFixed(pos_si_Gp, pos_si_rho1) < 0) // Gp rho1 (+ cos(beta))
  {
    HiggsRotFixed.row(pos_si_Gp) *= -1;
  }
  if (HiggsRotFixed(pos_si_Gm, pos_si_eta1) < 0) // Gm eta1 (+ cos(beta))
  {
    HiggsRotFixed.row(pos_si_Gm) *= -1;
  }
  if (HiggsRotFixed(pos_si_Hp, pos_si_rho2) < 0) // Hp rho2 (+ cos(beta))
  {
    HiggsRotFixed.row(pos_si_Hp) *= -1;
  }
  if (HiggsRotFixed(pos_si_Hm, pos_si_eta2) < 0) // Hm eta2 (+ cos(beta))
  {
    HiggsRotFixed.row(pos_si_Hm) *= -1;
  }

  // Check neutral submatrix
  // Use the "ScannerS" criteria from arXiv:2007.02985 Eq. (6)
  // (since they use the same parametrisation of the angles as BSMPT):
  // * (1) if R[1][1] < 0: h1 -> -h1 (i.e. multiply the h1 row with -1)
  // * (2) if R[3][3] < 0: h3 -> -h3 (i.e. multiply the h3 row with -1)
  // * (3) if det R < 0: h2 -> -h2 (i.e. multiply the h2 row with -1)

  // check neutral, CP-even submatrix
  if (HiggsRotFixed(pos_si_h1, pos_si_zeta1) < 0)
  // h1 zeta1 (condition (1) above, R11 < 0)
  {
    // if negative, flip sign of h1
    HiggsRotFixed.row(pos_si_h1) *= -1;
  }

  if (HiggsRotFixed(pos_si_h3, pos_si_zeta3) < 0)
  // h3 zeta3 (condition (2) above, R33 < 0)
  {
    // if negative, flip sign of h3
    HiggsRotFixed.row(pos_si_h3) *= -1;
  }

  // Calculate the determinant AFTER flipping the signs for rows 1 and 3 above
  MatrixXd HiggsRotFixedNeutral(3, 3);
  HiggsRotFixedNeutral(0, 0) = HiggsRotFixed(pos_si_h1, pos_si_zeta1);
  HiggsRotFixedNeutral(0, 1) = HiggsRotFixed(pos_si_h1, pos_si_zeta2);
  HiggsRotFixedNeutral(0, 2) = HiggsRotFixed(pos_si_h1, pos_si_zeta3);

  HiggsRotFixedNeutral(1, 0) = HiggsRotFixed(pos_si_h2, pos_si_zeta1);
  HiggsRotFixedNeutral(1, 1) = HiggsRotFixed(pos_si_h2, pos_si_zeta2);
  HiggsRotFixedNeutral(1, 2) = HiggsRotFixed(pos_si_h2, pos_si_zeta3);

  HiggsRotFixedNeutral(2, 0) = HiggsRotFixed(pos_si_h3, pos_si_zeta1);
  HiggsRotFixedNeutral(2, 1) = HiggsRotFixed(pos_si_h3, pos_si_zeta2);
  HiggsRotFixedNeutral(2, 2) = HiggsRotFixed(pos_si_h3, pos_si_zeta3);

  if (HiggsRotFixedNeutral.determinant() < 0)
  // condition (3) above, det(R) < 0
  {
    // if negative, flip sign of h2
    HiggsRotFixed.row(pos_si_h2) *= -1;
  }

  // Undo Goldstone rotation and flip that made G0 the first element
  MatrixXd HiggsRotFixedGoldstone(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsRotFixedGoldstone.row(i) = HiggsRotFixed.row(i);
  }

  HiggsRotFixed = MoveGoldstoneFirst * HiggsRotFixedGoldstone * RotGoldstone;

  // Extract the fixed mixing angles
  double sina2 = HiggsRotFixedGoldstone(pos_si_h1, pos_si_zeta3); // +sin(a2)
  double cosa2 = std::sqrt(1.0 - sina2 * sina2);
  alpha1       = std::asin(HiggsRotFixedGoldstone(pos_si_h1, pos_si_zeta2) /
                     cosa2); // +sin(a1) cos(a2)
  alpha2       = std::asin(sina2);
  alpha3       = std::asin(HiggsRotFixedGoldstone(pos_si_h2, pos_si_zeta3) /
                     cosa2); // +cos(a2) sin(a3)

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRotationMatrixEnsuredConvention[i][j] = HiggsRotFixed(i, j);
    }
  }

  return;
}

void Class_Potential_C2HDM::TripleHiggsCouplings()
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
  HiggsOrder[4] = pos_G0;
  if (UseHsmNotationInTripleHiggs)
  {
    HiggsOrder[5] = pos_h_SM;
    HiggsOrder[6] = pos_h_l;
    HiggsOrder[7] = pos_h_H;
  }
  else
  {
    HiggsOrder[5] = pos_h1;
    HiggsOrder[6] = pos_h2;
    HiggsOrder[7] = pos_h3;
  }

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

  HiggsVev[4] = SMConstants.C_vev0 * C_CosBeta;
  HiggsVev[6] = SMConstants.C_vev0 * C_SinBeta;

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

  Curvature_Gauge_G2H2[0][0][0][0] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][4][4] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][5][5] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][6][6] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][7][7] = SMConstants.C_g * SMConstants.C_g / 0.2e1;

  Curvature_Gauge_G2H2[0][3][0][4] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][1][5] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][2][6] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][3][7] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][4][0] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][5][1] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][6][2] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[0][3][7][3] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;

  Curvature_Gauge_G2H2[1][1][0][0] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][4][4] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][5][5] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][6][6] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][7][7] = SMConstants.C_g * SMConstants.C_g / 0.2e1;

  Curvature_Gauge_G2H2[1][3][0][5] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][1][4] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][2][7] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][3][6] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][4][1] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][5][0] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][6][3] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[1][3][7][2] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;

  Curvature_Gauge_G2H2[2][2][0][0] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][4][4] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][5][5] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][6][6] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][7][7] = SMConstants.C_g * SMConstants.C_g / 0.2e1;

  Curvature_Gauge_G2H2[2][3][0][0] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][1][1] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][2][2] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][3][3] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][4][4] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][5][5] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][6][6] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][7][7] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][0][4] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][1][5] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][2][6] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][3][7] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][4][0] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][5][1] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][6][2] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][0][7][3] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][0][5] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][1][4] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][2][7] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][3][6] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][4][1] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][5][0] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][6][3] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][1][7][2] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][0][0] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][1][1] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][2][2] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][3][3] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][4][4] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][5][5] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][6][6] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][7][7] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;

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

  double v1 = SMConstants.C_vev0 * C_CosBeta;
  double v2 = SMConstants.C_vev0 * C_SinBeta;
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

  YIJR2(0, 9)  = -std::conj(V11) * SMConstants.C_MassUp / v2;
  YIJR2(0, 10) = -std::conj(V12) * SMConstants.C_MassUp / v2;
  YIJR2(0, 11) = -std::conj(V13) * SMConstants.C_MassUp / v2;

  YIJR2(1, 9)  = -std::conj(V21) * SMConstants.C_MassCharm / v2;
  YIJR2(1, 10) = -std::conj(V22) * SMConstants.C_MassCharm / v2;
  YIJR2(1, 11) = -std::conj(V23) * SMConstants.C_MassCharm / v2;

  YIJR2(2, 9)  = -std::conj(V31) * SMConstants.C_MassTop / v2;
  YIJR2(2, 10) = -std::conj(V32) * SMConstants.C_MassTop / v2;
  YIJR2(2, 11) = -std::conj(V33) * SMConstants.C_MassTop / v2;

  YIJS2(0, 6) = SMConstants.C_MassUp / v2;
  YIJS2(1, 7) = SMConstants.C_MassCharm / v2;
  YIJS2(2, 8) = SMConstants.C_MassTop / v2;

  YIJSD(3, 9)  = SMConstants.C_MassDown / vD;
  YIJSD(4, 10) = SMConstants.C_MassStrange / vD;
  YIJSD(5, 11) = SMConstants.C_MassBottom / vD;

  YIJRD(3, 6) = V11 * SMConstants.C_MassDown / vD;
  YIJRD(3, 7) = V21 * SMConstants.C_MassDown / vD;
  YIJRD(3, 8) = V31 * SMConstants.C_MassDown / vD;
  YIJRD(4, 6) = V12 * SMConstants.C_MassStrange / vD;
  YIJRD(4, 7) = V22 * SMConstants.C_MassStrange / vD;
  YIJRD(4, 8) = V32 * SMConstants.C_MassStrange / vD;
  YIJRD(5, 6) = V13 * SMConstants.C_MassBottom / vD;
  YIJRD(5, 7) = V23 * SMConstants.C_MassBottom / vD;
  YIJRD(5, 8) = V33 * SMConstants.C_MassBottom / vD;

  YIJRL(1, 6) = SMConstants.C_MassElectron / vL;
  YIJRL(3, 7) = SMConstants.C_MassMu / vL;
  YIJRL(5, 8) = SMConstants.C_MassTau / vL;

  YIJSL(0, 1) = SMConstants.C_MassElectron / vL;
  YIJSL(2, 3) = SMConstants.C_MassMu / vL;
  YIJSL(4, 5) = SMConstants.C_MassTau / vL;

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
    cb = std::sqrt(2) * SMConstants.C_MassBottom /
         (SMConstants.C_vev0 * C_SinBeta);
  }
  if (Type == 2 or Type == 4) // Type II 2HDM oder Flipped
  {
    cb = std::sqrt(2) * SMConstants.C_MassBottom /
         (SMConstants.C_vev0 * C_CosBeta);
  }
  CTempC1 = 1.0 / 48 *
            (12 * L1 + 8 * L3 + 4 * L4 +
             3 * (3 * SMConstants.C_g * SMConstants.C_g +
                  SMConstants.C_gs * SMConstants.C_gs));
  double ct =
      std::sqrt(2) * SMConstants.C_MassTop / (SMConstants.C_vev0 * C_SinBeta);
  CTempC2 = 1.0 / 48 *
            (12 * L2 + 8 * L3 + 4 * L4 +
             3 * (3 * SMConstants.C_g * SMConstants.C_g +
                  SMConstants.C_gs * SMConstants.C_gs) +
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
  DebyeGauge[0][0] = 2 * SMConstants.C_g * SMConstants.C_g;
  DebyeGauge[1][1] = 2 * SMConstants.C_g * SMConstants.C_g;
  DebyeGauge[2][2] = 2 * SMConstants.C_g * SMConstants.C_g;
  DebyeGauge[3][3] = 2 * SMConstants.C_gs * SMConstants.C_gs;

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

  vcb = v[2];
  v1  = v[4];
  v2  = v[6];
  vcp = v[7];

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
