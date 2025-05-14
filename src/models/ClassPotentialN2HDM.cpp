// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
  @file
 */

#include <BSMPT/models/ClassPotentialN2HDM.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
#include <iomanip>

#include <optional>

using namespace Eigen;

namespace BSMPT
{
namespace Models
{

Class_Potential_N2HDM::Class_Potential_N2HDM(const ISMConstants &smConstants)
    : Class_Potential_Origin(smConstants)
{
  // TODO Auto-generated constructor stub
  Model         = ModelID::ModelIDs::N2HDM;
  NNeutralHiggs = 5;
  NChargedHiggs = 4;

  nPar   = 12;
  nParCT = 15;

  nVEV = 5;

  NHiggs  = NNeutralHiggs + NChargedHiggs;
  NGauge  = 4;
  NLepton = 9;
  NQuarks = 12;

  VevOrder.resize(nVEV);
  VevOrder[0] = 1; // Charged
  VevOrder[1] = 5; // CP breaking
  VevOrder[2] = 6; // v1
  VevOrder[3] = 7; // v2
  VevOrder[4] = 8; // vs

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_Potential_N2HDM::~Class_Potential_N2HDM()
{
  // TODO Auto-generated destructor stub
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given inputfile
 */
std::vector<std::string> Class_Potential_N2HDM::addLegendCT() const
{
  std::vector<std::string> labels;
  labels.push_back("Du11sq");
  labels.push_back("Du22sq");
  labels.push_back("Dm12sq");
  labels.push_back("DL1");
  labels.push_back("DL2");
  labels.push_back("DL3");
  labels.push_back("DL4");
  labels.push_back("DL5");
  labels.push_back("Dussq");
  labels.push_back("DL6");
  labels.push_back("DL7");
  labels.push_back("DL8");
  labels.push_back("DT1");
  labels.push_back("DT2");
  labels.push_back("DTs");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * inputfile
 */
std::vector<std::string> Class_Potential_N2HDM::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c");
  labels.push_back("omega_c");
  labels.push_back("omega_c/T_c");

  labels.push_back("omega_{CB}(T_c)");
  labels.push_back("omega_{CP}(T_c)");
  labels.push_back("omega_1(T_c)");
  labels.push_back("omega_2(T_c)");
  labels.push_back("omega_s(T_c)");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * higgs couplings. Use this to complement the legend of the given inputfile
 *
 * Beware, this is not implemented yet!
 */
std::vector<std::string> Class_Potential_N2HDM::addLegendTripleCouplings() const
{
  std::vector<std::string> particles;
  std::vector<std::string> labels;

  particles.push_back("G^+");
  particles.push_back("G^-");
  particles.push_back("H^+");
  particles.push_back("H^-");
  particles.push_back("G^0");
  particles.push_back("A");
  particles.push_back("h_SM");
  particles.push_back("h_l");
  particles.push_back("h_H");
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
 * Use this to complement the legend of the given inputfile
 */
std::vector<std::string> Class_Potential_N2HDM::addLegendVEV() const
{
  std::vector<std::string> labels;
  labels.push_back("omega_{CB}");
  labels.push_back("omega_{CP}");
  labels.push_back("omega_1");
  labels.push_back("omega_2");
  labels.push_back("omega_s");
  return labels;
}

void Class_Potential_N2HDM::ReadAndSet(const std::string &linestr,
                                       std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;
  double lu3 = 0, ll1 = 0, ll2 = 0, ll3 = 0, ll4 = 0, ll5 = 0, ll6 = 0, ll7 = 0,
         ll8 = 0, lvs = 0, ltanbeta = 0;
  double ltype = 0;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 12; k++)
  {
    ss >> tmp;
    if (k == 1)
      ltype = tmp;
    else if (k == 2)
      ll1 = tmp;
    else if (k == 3)
      ll2 = tmp;
    else if (k == 4)
      ll3 = tmp;
    else if (k == 5)
      ll4 = tmp;
    else if (k == 6)
      ll5 = tmp;
    else if (k == 7)
      ll6 = tmp;
    else if (k == 8)
      ll7 = tmp;
    else if (k == 9)
      ll8 = tmp;
    else if (k == 10)
      lvs = tmp;
    else if (k == 11)
      ltanbeta = tmp;
    else if (k == 12)
      lu3 = tmp;
  }
  par[8]  = ltanbeta;
  par[10] = lu3;
  par[9]  = lvs;
  par[11] = ltype;
  par[0]  = ll1;
  par[1]  = ll2;
  par[2]  = ll3;
  par[3]  = ll4;
  par[4]  = ll5;
  par[5]  = ll6;
  par[6]  = ll7;
  par[7]  = ll8;

  set_gen(par);
  return;
}

void Class_Potential_N2HDM::set_gen(const std::vector<double> &p)
{

  L1       = p[0];
  L2       = p[1];
  L3       = p[2];
  L4       = p[3];
  RL5      = p[4];
  NL6      = p[5];
  NL7      = p[6];
  NL8      = p[7];
  TanBeta  = p[8];
  Nvs      = p[9];
  RealMMix = p[10];
  scale    = SMConstants.C_vev0;
  Type     = p[11];

  C_CosBetaSquared = 1.0 / (1 + TanBeta * TanBeta);
  C_CosBeta        = sqrt(C_CosBetaSquared);
  C_SinBetaSquared = TanBeta * TanBeta * C_CosBetaSquared;
  C_SinBeta        = sqrt(C_SinBetaSquared);
  beta             = std::atan(TanBeta);

  double v1 = SMConstants.C_vev0 * C_CosBeta;
  double v2 = SMConstants.C_vev0 * C_SinBeta;

  u1 = RealMMix * TanBeta - 0.5 * v1 * v1 * L1 -
       0.5 * v2 * v2 * (L3 + L4 + RL5) - 0.5 * Nvs * Nvs * NL7;
  u2 = RealMMix / TanBeta - 0.5 * v2 * v2 * L2 -
       0.5 * v1 * v1 * (L3 + L4 + RL5) - 0.5 * Nvs * Nvs * NL8;
  Nus = -0.5 * NL6 * Nvs * Nvs - 0.5 * NL7 * v1 * v1 - 0.5 * NL8 * v2 * v2;

  vevTreeMin.resize(nVEV);
  vevTreeMin[0] = 0;
  vevTreeMin[1] = 0;
  vevTreeMin[2] = SMConstants.C_vev0 * C_CosBeta;
  vevTreeMin[3] = SMConstants.C_vev0 * C_SinBeta;
  vevTreeMin[4] = Nvs;
  vevTree.resize(NHiggs);
  vevTree = MinimizeOrderVEV(vevTreeMin);
}

void Class_Potential_N2HDM::set_CT_Pot_Par(const std::vector<double> &p)
{
  //	double *p = (double *)par;

  Du1CT  = p[0];
  Du2CT  = p[1];
  DRu3CT = p[2];
  DL1CT  = p[3];
  DL2CT  = p[4];
  DL3CT  = p[5];
  DL4CT  = p[6];
  DRL5CT = p[7];
  NDus   = p[8];
  NDL6   = p[9];
  NDL7   = p[10];
  NDL8   = p[11];
  DT1    = p[12];
  DT2    = p[13];
  NDTS   = p[14];

  Curvature_Higgs_CT_L1[0] = 0;
  Curvature_Higgs_CT_L1[1] = DTCharged;
  Curvature_Higgs_CT_L1[2] = 0;
  Curvature_Higgs_CT_L1[3] = 0;
  Curvature_Higgs_CT_L1[4] = 0;
  Curvature_Higgs_CT_L1[5] = DT3;
  Curvature_Higgs_CT_L1[6] = DT1;
  Curvature_Higgs_CT_L1[7] = DT2;
  Curvature_Higgs_CT_L1[8] = NDTS;

  Curvature_Higgs_CT_L2[0][0] = Du1CT;
  Curvature_Higgs_CT_L2[0][1] = -DRu3CT;
  Curvature_Higgs_CT_L2[1][1] = Du2CT;
  Curvature_Higgs_CT_L2[2][2] = Du1CT;
  Curvature_Higgs_CT_L2[2][3] = -DRu3CT;
  Curvature_Higgs_CT_L2[3][3] = Du2CT;
  Curvature_Higgs_CT_L2[4][4] = Du1CT;
  Curvature_Higgs_CT_L2[4][5] = -DRu3CT;
  Curvature_Higgs_CT_L2[5][5] = Du2CT;
  Curvature_Higgs_CT_L2[6][6] = Du1CT;
  Curvature_Higgs_CT_L2[6][7] = -DRu3CT;
  Curvature_Higgs_CT_L2[7][7] = Du2CT;
  Curvature_Higgs_CT_L2[8][8] = NDus;

  Curvature_Higgs_CT_L4[0][0][0][0] = 3 * DL1CT;
  Curvature_Higgs_CT_L4[0][0][1][1] = DL3CT + DL4CT + DRL5CT;
  Curvature_Higgs_CT_L4[0][0][2][2] = DL1CT;
  Curvature_Higgs_CT_L4[0][0][3][3] = DL3CT + DL4CT - DRL5CT;
  Curvature_Higgs_CT_L4[0][0][4][4] = DL1CT;
  Curvature_Higgs_CT_L4[0][0][5][5] = DL3CT;
  Curvature_Higgs_CT_L4[0][0][6][6] = DL1CT;
  Curvature_Higgs_CT_L4[0][0][7][7] = DL3CT;
  Curvature_Higgs_CT_L4[0][0][8][8] = NDL7;
  Curvature_Higgs_CT_L4[0][1][2][3] = DRL5CT;
  Curvature_Higgs_CT_L4[0][1][4][5] = DL4CT / 0.2e1 + DRL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[0][1][6][7] = DL4CT / 0.2e1 + DRL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][4][7] = -DL4CT / 0.2e1 + DRL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][5][6] = DL4CT / 0.2e1 - DRL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[1][1][1][1] = 3 * DL2CT;
  Curvature_Higgs_CT_L4[1][1][2][2] = DL3CT + DL4CT - DRL5CT;
  Curvature_Higgs_CT_L4[1][1][3][3] = DL2CT;
  Curvature_Higgs_CT_L4[1][1][4][4] = DL3CT;
  Curvature_Higgs_CT_L4[1][1][5][5] = DL2CT;
  Curvature_Higgs_CT_L4[1][1][6][6] = DL3CT;
  Curvature_Higgs_CT_L4[1][1][7][7] = DL2CT;
  Curvature_Higgs_CT_L4[1][1][8][8] = NDL8;
  Curvature_Higgs_CT_L4[1][2][4][7] = DL4CT / 0.2e1 - DRL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][5][6] = -DL4CT / 0.2e1 + DRL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[2][2][2][2] = 3 * DL1CT;
  Curvature_Higgs_CT_L4[2][2][3][3] = DL3CT + DL4CT + DRL5CT;
  Curvature_Higgs_CT_L4[2][2][4][4] = DL1CT;
  Curvature_Higgs_CT_L4[2][2][5][5] = DL3CT;
  Curvature_Higgs_CT_L4[2][2][6][6] = DL1CT;
  Curvature_Higgs_CT_L4[2][2][7][7] = DL3CT;
  Curvature_Higgs_CT_L4[2][2][8][8] = NDL7;
  Curvature_Higgs_CT_L4[2][3][4][5] = DL4CT / 0.2e1 + DRL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[2][3][6][7] = DL4CT / 0.2e1 + DRL5CT / 0.2e1;
  Curvature_Higgs_CT_L4[3][3][3][3] = 3 * DL2CT;
  Curvature_Higgs_CT_L4[3][3][4][4] = DL3CT;
  Curvature_Higgs_CT_L4[3][3][5][5] = DL2CT;
  Curvature_Higgs_CT_L4[3][3][6][6] = DL3CT;
  Curvature_Higgs_CT_L4[3][3][7][7] = DL2CT;
  Curvature_Higgs_CT_L4[3][3][8][8] = NDL8;
  Curvature_Higgs_CT_L4[4][4][4][4] = 3 * DL1CT;
  Curvature_Higgs_CT_L4[4][4][5][5] = DL3CT + DL4CT + DRL5CT;
  Curvature_Higgs_CT_L4[4][4][6][6] = DL1CT;
  Curvature_Higgs_CT_L4[4][4][7][7] = DL3CT + DL4CT - DRL5CT;
  Curvature_Higgs_CT_L4[4][4][8][8] = NDL7;
  Curvature_Higgs_CT_L4[4][5][6][7] = DRL5CT;
  Curvature_Higgs_CT_L4[5][5][5][5] = 3 * DL2CT;
  Curvature_Higgs_CT_L4[5][5][6][6] = DL3CT + DL4CT - DRL5CT;
  Curvature_Higgs_CT_L4[5][5][7][7] = DL2CT;
  Curvature_Higgs_CT_L4[5][5][8][8] = NDL8;
  Curvature_Higgs_CT_L4[6][6][6][6] = 3 * DL1CT;
  Curvature_Higgs_CT_L4[6][6][7][7] = DL3CT + DL4CT + DRL5CT;
  Curvature_Higgs_CT_L4[6][6][8][8] = NDL7;
  Curvature_Higgs_CT_L4[7][7][7][7] = 3 * DL2CT;
  Curvature_Higgs_CT_L4[7][7][8][8] = NDL8;
  Curvature_Higgs_CT_L4[8][8][8][8] = 3 * NDL6;

  for (std::size_t k1 = 0; k1 < NHiggs; k1++)
  {
    for (std::size_t k2 = k1; k2 < NHiggs; k2++)
    {
      Curvature_Higgs_CT_L2[k2][k1] = Curvature_Higgs_CT_L2[k1][k2];
      for (std::size_t k3 = k2; k3 < NHiggs; k3++)
      {
        for (std::size_t k4 = k3; k4 < NHiggs; k4++)
        {
          Curvature_Higgs_CT_L4[k2][k3][k4][k1] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k3][k4][k1][k2] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k4][k1][k2][k3] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k2][k1][k3][k4] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k4][k2][k1][k3] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k3][k4][k2][k1] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k1][k3][k4][k2] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k3][k2][k1][k4] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k4][k3][k2][k1] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k1][k4][k3][k2] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k2][k1][k4][k3] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k4][k2][k3][k1] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k1][k4][k2][k3] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k3][k1][k4][k2] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k2][k3][k1][k4] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k1][k3][k2][k4] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k4][k1][k3][k2] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k2][k4][k1][k3] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k3][k2][k4][k1] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k1][k2][k4][k3] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k3][k1][k2][k4] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k4][k3][k1][k2] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
          Curvature_Higgs_CT_L4[k2][k4][k3][k1] =
              Curvature_Higgs_CT_L4[k1][k2][k3][k4];
        }
      }
    }
  }

  return;
}

/**
 * Console-Output of all Parameters
 */
void Class_Potential_N2HDM::write() const
{
  std::stringstream ss;
  typedef std::numeric_limits<double> dbl;
  ss.precision(dbl::max_digits10);
  ss << std::setprecision(16);

  ss << "scale = " << scale << "\n";

  ss << "The parameters are :   \n";
  ss << "Model = " << Model << "\n";
  ss << "Renorm Scale = " << scale << "\n";
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
  ss << "Re(m_12^2) = " << RealMMix << std::endl;
  ss << "m_{11}^2 = " << u1 << std::endl;
  ss << "m_{22}^2 = " << u2 << std::endl;
  ss << "L6 = " << NL6 << "\n";
  ss << "L7 = " << NL7 << "\n";
  ss << "L8 = " << NL8 << "\n";
  ss << "m_s^2 = " << Nus << "\n";
  ss << "v_s = " << Nvs << "\n";

  ss << "The counterterms are :\n";

  ss << "DL1 := " << DL1CT << ";\n";
  ss << "DL2 := " << DL2CT << ";\n";
  ss << "DL3 := " << DL3CT << ";\n";
  ss << "DL4 := " << DL4CT << ";\n";
  ss << "DRL5 := " << DRL5CT << ";\n";

  ss << "Du1 := " << Du1CT << ";\n";
  ss << "Du2 := " << Du2CT << ";\n";
  ss << "DRu3 := " << DRu3CT << ";\n";

  ss << "DT1 := " << DT1 << ";\n";
  ss << "DT2 := " << DT2 << ";\n";
  ss << "DT3:= " << DT3 << ";\n";

  ss << "Dus := " << NDus << ";\n";
  ss << "DL6 := " << NDL6 << ";\n";
  ss << "DL7 := " << NDL7 << ";\n";
  ss << "DL8 := " << NDL8 << ";\n";
  ss << "DTS := " << NDTS << ";\n";

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

  MatrixXd NeutralMixing(3, 3);
  for (int i = 0; i < 3; i++)
  {
    NeutralMixing(0, i) = HiggsRot(pos_h_SM, i + 6);
    NeutralMixing(1, i) = HiggsRot(pos_h_l, i + 6);
    NeutralMixing(2, i) = HiggsRot(pos_h_H, i + 6);
  }

  ss << "The mass spectrum is given by :\n";
  ss << "m_{H^+} = " << std::sqrt(HiggsMasses[pos_Hp]) << " GeV \n"
     << "m_A = " << std::sqrt(HiggsMasses[pos_A]) << " GeV \n"
     << "m_{H_SM} = " << std::sqrt(HiggsMasses[pos_h_SM]) << " GeV \n"
     << "m_{H_l} = " << std::sqrt(HiggsMasses[pos_h_l]) << " GeV \n"
     << "m_{H_h} = " << std::sqrt(HiggsMasses[pos_h_H]) << " GeV \n";

  ss << "The neutral mixing Matrix is given by :\n";
  bool IsNegative = NeutralMixing(0, 1) < 0;
  ss << "H_{SM} = " << NeutralMixing(0, 0) << " zeta_1 ";
  if (IsNegative)
    ss << "-";
  else
    ss << "+";
  ss << std::abs(NeutralMixing(0, 1)) << " zeta_2 ";
  IsNegative = NeutralMixing(0, 2) < 0;
  if (IsNegative)
    ss << "-";
  else
    ss << "+";
  ss << std::abs(NeutralMixing(0, 2)) << " zeta_3 \n"
     << "H_{l} = " << NeutralMixing(1, 0) << " zeta_1 ";
  IsNegative = NeutralMixing(1, 1) < 0;
  if (IsNegative)
    ss << "-";
  else
    ss << "+";
  ss << std::abs(NeutralMixing(1, 1)) << " zeta_2 ";
  IsNegative = NeutralMixing(1, 2) < 0;
  if (IsNegative)
    ss << "-";
  else
    ss << "+";
  ss << std::abs(NeutralMixing(1, 2)) << " zeta_3 \n"
     << "H_{h} = " << NeutralMixing(2, 0) << " zeta_1 ";
  IsNegative = NeutralMixing(2, 1) < 0;
  if (IsNegative)
    ss << "-";
  else
    ss << "+";
  ss << std::abs(NeutralMixing(2, 1)) << " zeta_2 ";
  IsNegative = NeutralMixing(2, 2) < 0;
  if (IsNegative)
    ss << "-";
  else
    ss << "+";
  ss << std::abs(NeutralMixing(2, 2)) << " zeta_3 \n";
  Logger::Write(LoggingLevel::Default, ss.str());
}

std::vector<double> Class_Potential_N2HDM::calc_CT() const
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

  VectorXd NablaWeinberg(NHiggs);
  MatrixXd HesseWeinberg(NHiggs, NHiggs), HiggsRot(NHiggs, NHiggs);

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    NablaWeinberg(i) = WeinbergNabla[i];
    for (std::size_t j = 0; j < NHiggs; j++)
      HesseWeinberg(i, j) = WeinbergHesse.at(j * NHiggs + i);
  }

  double tL4  = 0; // Value of dL4
  double tDTS = 0; // Value of DTS

  double VB0[15], VB1[15], VB2[15];
  double VSol[15];
  if (Nvs != 0)
  {
    VB0[0] =
        (double)((Nvs * HesseWeinberg(6, 8) - HesseWeinberg(2, 3) * v2 +
                  2 * HesseWeinberg(4, 4) * v1 + HesseWeinberg(6, 7) * v2 -
                  5 * HesseWeinberg(2, 2) * v1 + HesseWeinberg(6, 6) * v1) /
                 v1) /
        0.2e1;
    VB0[1] = (double)((Nvs * HesseWeinberg(7, 8) * v2 +
                       HesseWeinberg(7, 7) * v2 * v2 -
                       3 * HesseWeinberg(5, 5) * v2 * v2 -
                       HesseWeinberg(2, 3) * v1 * v2 +
                       5 * v1 * v1 * HesseWeinberg(4, 4) +
                       v1 * v2 * HesseWeinberg(6, 7) -
                       5 * HesseWeinberg(2, 2) * v1 * v1) *
                      (double)pow((double)v2, (double)(-2))) /
             0.2e1;
    VB0[2] = (HesseWeinberg(2, 3) * v2 + HesseWeinberg(4, 4) * v1 -
              HesseWeinberg(2, 2) * v1) /
             v2;
    VB0[3] =
        (-HesseWeinberg(4, 4) + 2 * HesseWeinberg(2, 2) - HesseWeinberg(6, 6)) *
        (double)pow((double)v1, (double)(-2));
    VB0[4] = (-HesseWeinberg(7, 7) * v2 * v2 + HesseWeinberg(5, 5) * v2 * v2 -
              2 * v1 * v1 * HesseWeinberg(4, 4) +
              2 * HesseWeinberg(2, 2) * v1 * v1) *
             (double)pow((double)v2, (double)(-4));
    VB0[5] = (HesseWeinberg(2, 2) * v1 + HesseWeinberg(2, 3) * v2 -
              HesseWeinberg(4, 4) * v1 - HesseWeinberg(6, 7) * v2) /
             v1 * (double)pow((double)v2, (double)(-2));
    VB0[6] = 0;
    VB0[7] = (2 * HesseWeinberg(4, 4) - 2 * HesseWeinberg(2, 2)) *
             (double)pow((double)v2, (double)(-2));
    VB0[8] = (double)((HesseWeinberg(8, 8) * Nvs + v2 * HesseWeinberg(7, 8) +
                       v1 * HesseWeinberg(6, 8) - 3 * NablaWeinberg(8)) /
                      Nvs) /
             0.2e1;
    VB0[9] = -(double)pow((double)Nvs, (double)(-3)) *
             (HesseWeinberg(8, 8) * Nvs - NablaWeinberg(8));
    VB0[10] = -1 / v1 / Nvs * HesseWeinberg(6, 8);
    VB0[11] = -1 / v2 / Nvs * HesseWeinberg(7, 8);
    VB0[12] =
        HesseWeinberg(2, 2) * v1 + HesseWeinberg(2, 3) * v2 - NablaWeinberg(6);
    VB0[13] = (HesseWeinberg(0, 0) * v1 * v1 + HesseWeinberg(2, 3) * v1 * v2 -
               v1 * v1 * HesseWeinberg(4, 4) + HesseWeinberg(5, 5) * v2 * v2 -
               NablaWeinberg(7) * v2) /
              v2;
    VB0[14] = 0;

    VB1[0]  = v2 * v2;
    VB1[1]  = v1 * v1;
    VB1[2]  = v1 * v2;
    VB1[3]  = -pow(v1, -0.2e1) * v2 * v2;
    VB1[4]  = -pow(v2, -0.2e1) * v1 * v1;
    VB1[5]  = -1;
    VB1[6]  = 1;
    VB1[7]  = 1;
    VB1[8]  = 0;
    VB1[9]  = 0;
    VB1[10] = 0;
    VB1[11] = 0;
    VB1[12] = 0;
    VB1[13] = 0;
    VB1[14] = 0;

    VB2[0]  = 0;
    VB2[1]  = 0;
    VB2[2]  = 0;
    VB2[3]  = 0;
    VB2[4]  = 0;
    VB2[5]  = 0;
    VB2[6]  = 0;
    VB2[7]  = 0;
    VB2[8]  = -0.3e1 / 0.2e1 / Nvs;
    VB2[9]  = pow(Nvs, -0.3e1);
    VB2[10] = 0;
    VB2[11] = 0;
    VB2[12] = 0;
    VB2[13] = 0;
    VB2[14] = 1;

    for (int i = 0; i < 15; i++)
    {
      VSol[i] = VB0[i] + tL4 * VB1[i] + tDTS * VB2[i];
    }
  }

  //    Du1CT = VSol[0];
  //    Du2CT = VSol[1];
  //    DRu3CT = VSol[2];
  //    DL1CT = VSol[3];
  //    DL2CT = VSol[4];
  //    DL3CT = VSol[5];
  //    DL4CT = VSol[6];
  //    DRL5CT = VSol[7];
  //    NDus = VSol[8];
  //    NDL6 = VSol[9];
  //    NDL7 = VSol[10];
  //    NDL8 = VSol[11];
  //    DT1 = VSol[12];
  //    DT2 = VSol[13];
  //    NDTS = VSol[14];

  //  	for(int i=0;i<15;i++) par[i] = VSol[i];
  for (const auto &x : VSol)
    parCT.push_back(x);

  //  	double Identities[5];
  //  	Identities[0] = (HesseWeinberg(2, 3) - HesseWeinberg(0, 1));
  //  	Identities[1] = (HesseWeinberg(0, 0) - HesseWeinberg(2, 2));
  //  	Identities[2] = (HesseWeinberg(0, 0) * v1 * v1 + v2 * v2 *
  //  HesseWeinberg(5, 5) - v1 * v1 * HesseWeinberg(4, 4) - v2 * v2 *
  //  HesseWeinberg(3, 3)); 	Identities[3] = (HesseWeinberg(0, 0) * v1 -
  //  HesseWeinberg(4, 4) * v1 + HesseWeinberg(2, 3) * v2 - HesseWeinberg(4, 5)
  //  * v2); 	Identities[4] = (HesseWeinberg(0, 0) * v1 * v1 + v2 * v2 *
  //  HesseWeinberg(5, 5) - v1 * v1 * HesseWeinberg(4, 4) - v2 * v2 *
  //  HesseWeinberg(1, 1));

  return parCT;
}

/**
 * Ensures the correct rotation matrix convention
 */
void Class_Potential_N2HDM::AdjustRotationMatrix()
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

  // N2HDM interaction basis
  // 0     1     2     3     4     5     6      7      8
  // rho1, rho2, eta1, eta2, psi1, psi2, zeta1, zeta2, rhoS
  const std::size_t pos_rho1 = 0, pos_rho2 = 1, pos_eta1 = 2, pos_eta2 = 3,
                    pos_psi1 = 4, pos_psi2 = 5, pos_zeta1 = 6, pos_zeta2 = 7,
                    pos_rhoS = 8;

  // Unlike in the C2HDM, there is no CP-mixing, so we do not have to rotate out
  // the Goldstone

  // Indices of mass eigenstates for rotation from interaction to mass basis
  std::optional<std::size_t> tpos_G0, tpos_Gp, tpos_Gm, tpos_Hp, tpos_Hm,
                             tpos_h1, tpos_h2, tpos_h3, tpos_A;


  for (std::size_t i = 0; i < NHiggs; i++)
  // mass base index i corresponds to mass vector sorted in ascending mass
  {
    // Goldstones have zero mass in the Landau gauge
    bool hasZeroMass = std::abs(HiggsMasses[i]) < ZeroThreshold;
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
    // Neutral CP-odd submatrix
    else if (std::abs(HiggsRot(i, pos_psi1)) + std::abs(HiggsRot(i, pos_psi2)) >
             ZeroThreshold)
    {
      if (not tpos_G0.has_value() and hasZeroMass)
      {
        tpos_G0 = i;
      }
      else if (not tpos_A.has_value())
      {
        tpos_A = i;
      }
      else
      {
        throw std::runtime_error("Error. Neutral CP-odd submatrix mixing with "
                                 "other components.");
      }
    }
    // Neutral CP-even submatrix
    else if (std::abs(HiggsRot(i, pos_zeta1)) +
                 std::abs(HiggsRot(i, pos_zeta2)) +
                 std::abs(HiggsRot(i, pos_rhoS)) > ZeroThreshold)
    {
      // use that mh1 < mh2 < mh3
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
        throw std::runtime_error("Error. Neutral CP-even submatrix mixing with"
                                 " other components.");
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
           and tpos_h3.has_value() and tpos_A.has_value())
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
  pos_A = tpos_A.value();

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
              (j == pos_psi1 and (i == pos_G0 or i == pos_A)) or
              (j == pos_psi2 and (i == pos_G0 or i == pos_A)) or
              (j == pos_zeta1 and
               (i == pos_h1 or i == pos_h2 or i == pos_h3)) or
              (j == pos_zeta2 and
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

  MatrixXd HiggsRotFixed(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsRotFixed.row(i) = HiggsRot.row(i);
  }

  // Neutral CP-odd submatrix
  if (HiggsRotFixed(pos_G0, pos_psi1) < 0) // G0 psi1 (+ cos(beta))
  {
    HiggsRotFixed.row(pos_G0) *= -1;
  }
  if (HiggsRotFixed(pos_A, pos_psi2) < 0) // A psi2 (+ cos(beta))
  {
    HiggsRotFixed.row(pos_A) *= -1;
  }

  // charged submatrix
  if (HiggsRotFixed(pos_Gp, pos_rho1) < 0) // Gp rho1 (+ cos(beta))
  {
    HiggsRotFixed.row(pos_Gp) *= -1;
  }
  if (HiggsRotFixed(pos_Gm, pos_eta1) < 0) // Gm eta1 (+ cos(beta))
  {
    HiggsRotFixed.row(pos_Gm) *= -1;
  }
  if (HiggsRotFixed(pos_Hp, pos_rho2) < 0) // Hp rho2 (+ cos(beta))
  {
    HiggsRotFixed.row(pos_Hp) *= -1;
  }
  if (HiggsRotFixed(pos_Hm, pos_eta2) < 0) // Hm eta2 (+ cos(beta))
  {
    HiggsRotFixed.row(pos_Hm) *= -1;
  }

  // Check neutral CP-even submatrix
  // Use the "ScannerS" criteria from arXiv:2007.02985 Eq. (6)
  // (since ScannerS uses the same parametrisation of the angles as BSMPT):
  // * (1) if R[1][1] < 0: h1 -> -h1 (i.e. multiply the h1 row with -1)
  // * (2) if R[3][3] < 0: h3 -> -h3 (i.e. multiply the h3 row with -1)
  // * (3) if det R < 0: h2 -> -h2 (i.e. multiply the h2 row with -1)

  // check neutral, CP-even submatrix
  if (HiggsRotFixed(pos_h1, pos_zeta1) < 0)
  // h1 zeta1 (condition (1) above, R11 < 0)
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
  HiggsRotFixedNeutral(0, 0) = HiggsRotFixed(pos_h1, pos_zeta1);
  HiggsRotFixedNeutral(0, 1) = HiggsRotFixed(pos_h1, pos_zeta2);
  HiggsRotFixedNeutral(0, 2) = HiggsRotFixed(pos_h1, pos_rhoS);

  HiggsRotFixedNeutral(1, 0) = HiggsRotFixed(pos_h2, pos_zeta1);
  HiggsRotFixedNeutral(1, 1) = HiggsRotFixed(pos_h2, pos_zeta2);
  HiggsRotFixedNeutral(1, 2) = HiggsRotFixed(pos_h2, pos_rhoS);

  HiggsRotFixedNeutral(2, 0) = HiggsRotFixed(pos_h3, pos_zeta1);
  HiggsRotFixedNeutral(2, 1) = HiggsRotFixed(pos_h3, pos_zeta2);
  HiggsRotFixedNeutral(2, 2) = HiggsRotFixed(pos_h3, pos_rhoS);

  if (HiggsRotFixedNeutral.determinant() < 0)
  // condition (3) above, det(R) < 0
  {
    // if negative, flip sign of h2
    HiggsRotFixed.row(pos_h2) *= -1;
  }

  // Extract the fixed mixing angles
  double sina2 = HiggsRotFixed(pos_h1, pos_rhoS); // +sin(a2)
  double cosa2 = std::sqrt(1.0 - sina2 * sina2);
  alpha1 =
      std::asin(HiggsRotFixed(pos_h1, pos_zeta2) / cosa2); // +sin(a1) cos(a2)
  alpha2 = std::asin(sina2);
  alpha3 =
      std::asin(HiggsRotFixed(pos_h2, pos_rhoS) / cosa2); // +cos(a2) sin(a3)

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRotationMatrixEnsuredConvention[i][j] = HiggsRotFixed(i, j);
    }
  }

  return;
}

void Class_Potential_N2HDM::TripleHiggsCouplings()
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
  HiggsOrder[5] = pos_A;
  HiggsOrder[6] = pos_h_SM;
  HiggsOrder[7] = pos_h_l;
  HiggsOrder[8] = pos_h_H;

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

void Class_Potential_N2HDM::SetCurvatureArrays()
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

  HiggsVev[6]              = SMConstants.C_vev0 * C_CosBeta;
  HiggsVev[7]              = SMConstants.C_vev0 * C_SinBeta;
  HiggsVev[8]              = Nvs;
  Curvature_Higgs_L2[0][0] = u1;
  Curvature_Higgs_L2[0][1] = -RealMMix;
  Curvature_Higgs_L2[0][2] = 0;
  Curvature_Higgs_L2[0][3] = 0;
  Curvature_Higgs_L2[0][4] = 0;
  Curvature_Higgs_L2[0][5] = 0;
  Curvature_Higgs_L2[0][6] = 0;
  Curvature_Higgs_L2[0][7] = 0;
  Curvature_Higgs_L2[0][8] = 0;
  Curvature_Higgs_L2[1][0] = -RealMMix;
  Curvature_Higgs_L2[1][1] = u2;
  Curvature_Higgs_L2[1][2] = 0;
  Curvature_Higgs_L2[1][3] = 0;
  Curvature_Higgs_L2[1][4] = 0;
  Curvature_Higgs_L2[1][5] = 0;
  Curvature_Higgs_L2[1][6] = 0;
  Curvature_Higgs_L2[1][7] = 0;
  Curvature_Higgs_L2[1][8] = 0;
  Curvature_Higgs_L2[2][0] = 0;
  Curvature_Higgs_L2[2][1] = 0;
  Curvature_Higgs_L2[2][2] = u1;
  Curvature_Higgs_L2[2][3] = -RealMMix;
  Curvature_Higgs_L2[2][4] = 0;
  Curvature_Higgs_L2[2][5] = 0;
  Curvature_Higgs_L2[2][6] = 0;
  Curvature_Higgs_L2[2][7] = 0;
  Curvature_Higgs_L2[2][8] = 0;
  Curvature_Higgs_L2[3][0] = 0;
  Curvature_Higgs_L2[3][1] = 0;
  Curvature_Higgs_L2[3][2] = -RealMMix;
  Curvature_Higgs_L2[3][3] = u2;
  Curvature_Higgs_L2[3][4] = 0;
  Curvature_Higgs_L2[3][5] = 0;
  Curvature_Higgs_L2[3][6] = 0;
  Curvature_Higgs_L2[3][7] = 0;
  Curvature_Higgs_L2[3][8] = 0;
  Curvature_Higgs_L2[4][0] = 0;
  Curvature_Higgs_L2[4][1] = 0;
  Curvature_Higgs_L2[4][2] = 0;
  Curvature_Higgs_L2[4][3] = 0;
  Curvature_Higgs_L2[4][4] = u1;
  Curvature_Higgs_L2[4][5] = -RealMMix;
  Curvature_Higgs_L2[4][6] = 0;
  Curvature_Higgs_L2[4][7] = 0;
  Curvature_Higgs_L2[4][8] = 0;
  Curvature_Higgs_L2[5][0] = 0;
  Curvature_Higgs_L2[5][1] = 0;
  Curvature_Higgs_L2[5][2] = 0;
  Curvature_Higgs_L2[5][3] = 0;
  Curvature_Higgs_L2[5][4] = -RealMMix;
  Curvature_Higgs_L2[5][5] = u2;
  Curvature_Higgs_L2[5][6] = 0;
  Curvature_Higgs_L2[5][7] = 0;
  Curvature_Higgs_L2[5][8] = 0;
  Curvature_Higgs_L2[6][0] = 0;
  Curvature_Higgs_L2[6][1] = 0;
  Curvature_Higgs_L2[6][2] = 0;
  Curvature_Higgs_L2[6][3] = 0;
  Curvature_Higgs_L2[6][4] = 0;
  Curvature_Higgs_L2[6][5] = 0;
  Curvature_Higgs_L2[6][6] = u1;
  Curvature_Higgs_L2[6][7] = -RealMMix;
  Curvature_Higgs_L2[6][8] = 0;
  Curvature_Higgs_L2[7][0] = 0;
  Curvature_Higgs_L2[7][1] = 0;
  Curvature_Higgs_L2[7][2] = 0;
  Curvature_Higgs_L2[7][3] = 0;
  Curvature_Higgs_L2[7][4] = 0;
  Curvature_Higgs_L2[7][5] = 0;
  Curvature_Higgs_L2[7][6] = -RealMMix;
  Curvature_Higgs_L2[7][7] = u2;
  Curvature_Higgs_L2[7][8] = 0;
  Curvature_Higgs_L2[8][0] = 0;
  Curvature_Higgs_L2[8][1] = 0;
  Curvature_Higgs_L2[8][2] = 0;
  Curvature_Higgs_L2[8][3] = 0;
  Curvature_Higgs_L2[8][4] = 0;
  Curvature_Higgs_L2[8][5] = 0;
  Curvature_Higgs_L2[8][6] = 0;
  Curvature_Higgs_L2[8][7] = 0;
  Curvature_Higgs_L2[8][8] = Nus;

  Curvature_Higgs_L4[0][0][0][0] = 3 * L1;
  Curvature_Higgs_L4[0][0][0][1] = 0;
  Curvature_Higgs_L4[0][0][0][2] = 0;
  Curvature_Higgs_L4[0][0][0][3] = 0;
  Curvature_Higgs_L4[0][0][0][4] = 0;
  Curvature_Higgs_L4[0][0][0][5] = 0;
  Curvature_Higgs_L4[0][0][0][6] = 0;
  Curvature_Higgs_L4[0][0][0][7] = 0;
  Curvature_Higgs_L4[0][0][0][8] = 0;
  Curvature_Higgs_L4[0][0][1][1] = L3 + L4 + RL5;
  Curvature_Higgs_L4[0][0][1][2] = 0;
  Curvature_Higgs_L4[0][0][1][3] = 0;
  Curvature_Higgs_L4[0][0][1][4] = 0;
  Curvature_Higgs_L4[0][0][1][5] = 0;
  Curvature_Higgs_L4[0][0][1][6] = 0;
  Curvature_Higgs_L4[0][0][1][7] = 0;
  Curvature_Higgs_L4[0][0][1][8] = 0;
  Curvature_Higgs_L4[0][0][2][2] = L1;
  Curvature_Higgs_L4[0][0][2][3] = 0;
  Curvature_Higgs_L4[0][0][2][4] = 0;
  Curvature_Higgs_L4[0][0][2][5] = 0;
  Curvature_Higgs_L4[0][0][2][6] = 0;
  Curvature_Higgs_L4[0][0][2][7] = 0;
  Curvature_Higgs_L4[0][0][2][8] = 0;
  Curvature_Higgs_L4[0][0][3][3] = L3 + L4 - RL5;
  Curvature_Higgs_L4[0][0][3][4] = 0;
  Curvature_Higgs_L4[0][0][3][5] = 0;
  Curvature_Higgs_L4[0][0][3][6] = 0;
  Curvature_Higgs_L4[0][0][3][7] = 0;
  Curvature_Higgs_L4[0][0][3][8] = 0;
  Curvature_Higgs_L4[0][0][4][4] = L1;
  Curvature_Higgs_L4[0][0][4][5] = 0;
  Curvature_Higgs_L4[0][0][4][6] = 0;
  Curvature_Higgs_L4[0][0][4][7] = 0;
  Curvature_Higgs_L4[0][0][4][8] = 0;
  Curvature_Higgs_L4[0][0][5][5] = L3;
  Curvature_Higgs_L4[0][0][5][6] = 0;
  Curvature_Higgs_L4[0][0][5][7] = 0;
  Curvature_Higgs_L4[0][0][5][8] = 0;
  Curvature_Higgs_L4[0][0][6][6] = L1;
  Curvature_Higgs_L4[0][0][6][7] = 0;
  Curvature_Higgs_L4[0][0][6][8] = 0;
  Curvature_Higgs_L4[0][0][7][7] = L3;
  Curvature_Higgs_L4[0][0][7][8] = 0;
  Curvature_Higgs_L4[0][0][8][8] = NL7;
  Curvature_Higgs_L4[0][1][1][1] = 0;
  Curvature_Higgs_L4[0][1][1][2] = 0;
  Curvature_Higgs_L4[0][1][1][3] = 0;
  Curvature_Higgs_L4[0][1][1][4] = 0;
  Curvature_Higgs_L4[0][1][1][5] = 0;
  Curvature_Higgs_L4[0][1][1][6] = 0;
  Curvature_Higgs_L4[0][1][1][7] = 0;
  Curvature_Higgs_L4[0][1][1][8] = 0;
  Curvature_Higgs_L4[0][1][2][2] = 0;
  Curvature_Higgs_L4[0][1][2][3] = RL5;
  Curvature_Higgs_L4[0][1][2][4] = 0;
  Curvature_Higgs_L4[0][1][2][5] = 0;
  Curvature_Higgs_L4[0][1][2][6] = 0;
  Curvature_Higgs_L4[0][1][2][7] = 0;
  Curvature_Higgs_L4[0][1][2][8] = 0;
  Curvature_Higgs_L4[0][1][3][3] = 0;
  Curvature_Higgs_L4[0][1][3][4] = 0;
  Curvature_Higgs_L4[0][1][3][5] = 0;
  Curvature_Higgs_L4[0][1][3][6] = 0;
  Curvature_Higgs_L4[0][1][3][7] = 0;
  Curvature_Higgs_L4[0][1][3][8] = 0;
  Curvature_Higgs_L4[0][1][4][4] = 0;
  Curvature_Higgs_L4[0][1][4][5] = L4 / 0.2e1 + RL5 / 0.2e1;
  Curvature_Higgs_L4[0][1][4][6] = 0;
  Curvature_Higgs_L4[0][1][4][7] = 0;
  Curvature_Higgs_L4[0][1][4][8] = 0;
  Curvature_Higgs_L4[0][1][5][5] = 0;
  Curvature_Higgs_L4[0][1][5][6] = 0;
  Curvature_Higgs_L4[0][1][5][7] = 0;
  Curvature_Higgs_L4[0][1][5][8] = 0;
  Curvature_Higgs_L4[0][1][6][6] = 0;
  Curvature_Higgs_L4[0][1][6][7] = L4 / 0.2e1 + RL5 / 0.2e1;
  Curvature_Higgs_L4[0][1][6][8] = 0;
  Curvature_Higgs_L4[0][1][7][7] = 0;
  Curvature_Higgs_L4[0][1][7][8] = 0;
  Curvature_Higgs_L4[0][1][8][8] = 0;
  Curvature_Higgs_L4[0][2][2][2] = 0;
  Curvature_Higgs_L4[0][2][2][3] = 0;
  Curvature_Higgs_L4[0][2][2][4] = 0;
  Curvature_Higgs_L4[0][2][2][5] = 0;
  Curvature_Higgs_L4[0][2][2][6] = 0;
  Curvature_Higgs_L4[0][2][2][7] = 0;
  Curvature_Higgs_L4[0][2][2][8] = 0;
  Curvature_Higgs_L4[0][2][3][3] = 0;
  Curvature_Higgs_L4[0][2][3][4] = 0;
  Curvature_Higgs_L4[0][2][3][5] = 0;
  Curvature_Higgs_L4[0][2][3][6] = 0;
  Curvature_Higgs_L4[0][2][3][7] = 0;
  Curvature_Higgs_L4[0][2][3][8] = 0;
  Curvature_Higgs_L4[0][2][4][4] = 0;
  Curvature_Higgs_L4[0][2][4][5] = 0;
  Curvature_Higgs_L4[0][2][4][6] = 0;
  Curvature_Higgs_L4[0][2][4][7] = 0;
  Curvature_Higgs_L4[0][2][4][8] = 0;
  Curvature_Higgs_L4[0][2][5][5] = 0;
  Curvature_Higgs_L4[0][2][5][6] = 0;
  Curvature_Higgs_L4[0][2][5][7] = 0;
  Curvature_Higgs_L4[0][2][5][8] = 0;
  Curvature_Higgs_L4[0][2][6][6] = 0;
  Curvature_Higgs_L4[0][2][6][7] = 0;
  Curvature_Higgs_L4[0][2][6][8] = 0;
  Curvature_Higgs_L4[0][2][7][7] = 0;
  Curvature_Higgs_L4[0][2][7][8] = 0;
  Curvature_Higgs_L4[0][2][8][8] = 0;
  Curvature_Higgs_L4[0][3][3][3] = 0;
  Curvature_Higgs_L4[0][3][3][4] = 0;
  Curvature_Higgs_L4[0][3][3][5] = 0;
  Curvature_Higgs_L4[0][3][3][6] = 0;
  Curvature_Higgs_L4[0][3][3][7] = 0;
  Curvature_Higgs_L4[0][3][3][8] = 0;
  Curvature_Higgs_L4[0][3][4][4] = 0;
  Curvature_Higgs_L4[0][3][4][5] = 0;
  Curvature_Higgs_L4[0][3][4][6] = 0;
  Curvature_Higgs_L4[0][3][4][7] = -L4 / 0.2e1 + RL5 / 0.2e1;
  Curvature_Higgs_L4[0][3][4][8] = 0;
  Curvature_Higgs_L4[0][3][5][5] = 0;
  Curvature_Higgs_L4[0][3][5][6] = L4 / 0.2e1 - RL5 / 0.2e1;
  Curvature_Higgs_L4[0][3][5][7] = 0;
  Curvature_Higgs_L4[0][3][5][8] = 0;
  Curvature_Higgs_L4[0][3][6][6] = 0;
  Curvature_Higgs_L4[0][3][6][7] = 0;
  Curvature_Higgs_L4[0][3][6][8] = 0;
  Curvature_Higgs_L4[0][3][7][7] = 0;
  Curvature_Higgs_L4[0][3][7][8] = 0;
  Curvature_Higgs_L4[0][3][8][8] = 0;
  Curvature_Higgs_L4[0][4][4][4] = 0;
  Curvature_Higgs_L4[0][4][4][5] = 0;
  Curvature_Higgs_L4[0][4][4][6] = 0;
  Curvature_Higgs_L4[0][4][4][7] = 0;
  Curvature_Higgs_L4[0][4][4][8] = 0;
  Curvature_Higgs_L4[0][4][5][5] = 0;
  Curvature_Higgs_L4[0][4][5][6] = 0;
  Curvature_Higgs_L4[0][4][5][7] = 0;
  Curvature_Higgs_L4[0][4][5][8] = 0;
  Curvature_Higgs_L4[0][4][6][6] = 0;
  Curvature_Higgs_L4[0][4][6][7] = 0;
  Curvature_Higgs_L4[0][4][6][8] = 0;
  Curvature_Higgs_L4[0][4][7][7] = 0;
  Curvature_Higgs_L4[0][4][7][8] = 0;
  Curvature_Higgs_L4[0][4][8][8] = 0;
  Curvature_Higgs_L4[0][5][5][5] = 0;
  Curvature_Higgs_L4[0][5][5][6] = 0;
  Curvature_Higgs_L4[0][5][5][7] = 0;
  Curvature_Higgs_L4[0][5][5][8] = 0;
  Curvature_Higgs_L4[0][5][6][6] = 0;
  Curvature_Higgs_L4[0][5][6][7] = 0;
  Curvature_Higgs_L4[0][5][6][8] = 0;
  Curvature_Higgs_L4[0][5][7][7] = 0;
  Curvature_Higgs_L4[0][5][7][8] = 0;
  Curvature_Higgs_L4[0][5][8][8] = 0;
  Curvature_Higgs_L4[0][6][6][6] = 0;
  Curvature_Higgs_L4[0][6][6][7] = 0;
  Curvature_Higgs_L4[0][6][6][8] = 0;
  Curvature_Higgs_L4[0][6][7][7] = 0;
  Curvature_Higgs_L4[0][6][7][8] = 0;
  Curvature_Higgs_L4[0][6][8][8] = 0;
  Curvature_Higgs_L4[0][7][7][7] = 0;
  Curvature_Higgs_L4[0][7][7][8] = 0;
  Curvature_Higgs_L4[0][7][8][8] = 0;
  Curvature_Higgs_L4[0][8][8][8] = 0;
  Curvature_Higgs_L4[1][1][1][1] = 3 * L2;
  Curvature_Higgs_L4[1][1][1][2] = 0;
  Curvature_Higgs_L4[1][1][1][3] = 0;
  Curvature_Higgs_L4[1][1][1][4] = 0;
  Curvature_Higgs_L4[1][1][1][5] = 0;
  Curvature_Higgs_L4[1][1][1][6] = 0;
  Curvature_Higgs_L4[1][1][1][7] = 0;
  Curvature_Higgs_L4[1][1][1][8] = 0;
  Curvature_Higgs_L4[1][1][2][2] = L3 + L4 - RL5;
  Curvature_Higgs_L4[1][1][2][3] = 0;
  Curvature_Higgs_L4[1][1][2][4] = 0;
  Curvature_Higgs_L4[1][1][2][5] = 0;
  Curvature_Higgs_L4[1][1][2][6] = 0;
  Curvature_Higgs_L4[1][1][2][7] = 0;
  Curvature_Higgs_L4[1][1][2][8] = 0;
  Curvature_Higgs_L4[1][1][3][3] = L2;
  Curvature_Higgs_L4[1][1][3][4] = 0;
  Curvature_Higgs_L4[1][1][3][5] = 0;
  Curvature_Higgs_L4[1][1][3][6] = 0;
  Curvature_Higgs_L4[1][1][3][7] = 0;
  Curvature_Higgs_L4[1][1][3][8] = 0;
  Curvature_Higgs_L4[1][1][4][4] = L3;
  Curvature_Higgs_L4[1][1][4][5] = 0;
  Curvature_Higgs_L4[1][1][4][6] = 0;
  Curvature_Higgs_L4[1][1][4][7] = 0;
  Curvature_Higgs_L4[1][1][4][8] = 0;
  Curvature_Higgs_L4[1][1][5][5] = L2;
  Curvature_Higgs_L4[1][1][5][6] = 0;
  Curvature_Higgs_L4[1][1][5][7] = 0;
  Curvature_Higgs_L4[1][1][5][8] = 0;
  Curvature_Higgs_L4[1][1][6][6] = L3;
  Curvature_Higgs_L4[1][1][6][7] = 0;
  Curvature_Higgs_L4[1][1][6][8] = 0;
  Curvature_Higgs_L4[1][1][7][7] = L2;
  Curvature_Higgs_L4[1][1][7][8] = 0;
  Curvature_Higgs_L4[1][1][8][8] = NL8;
  Curvature_Higgs_L4[1][2][2][2] = 0;
  Curvature_Higgs_L4[1][2][2][3] = 0;
  Curvature_Higgs_L4[1][2][2][4] = 0;
  Curvature_Higgs_L4[1][2][2][5] = 0;
  Curvature_Higgs_L4[1][2][2][6] = 0;
  Curvature_Higgs_L4[1][2][2][7] = 0;
  Curvature_Higgs_L4[1][2][2][8] = 0;
  Curvature_Higgs_L4[1][2][3][3] = 0;
  Curvature_Higgs_L4[1][2][3][4] = 0;
  Curvature_Higgs_L4[1][2][3][5] = 0;
  Curvature_Higgs_L4[1][2][3][6] = 0;
  Curvature_Higgs_L4[1][2][3][7] = 0;
  Curvature_Higgs_L4[1][2][3][8] = 0;
  Curvature_Higgs_L4[1][2][4][4] = 0;
  Curvature_Higgs_L4[1][2][4][5] = 0;
  Curvature_Higgs_L4[1][2][4][6] = 0;
  Curvature_Higgs_L4[1][2][4][7] = L4 / 0.2e1 - RL5 / 0.2e1;
  Curvature_Higgs_L4[1][2][4][8] = 0;
  Curvature_Higgs_L4[1][2][5][5] = 0;
  Curvature_Higgs_L4[1][2][5][6] = -L4 / 0.2e1 + RL5 / 0.2e1;
  Curvature_Higgs_L4[1][2][5][7] = 0;
  Curvature_Higgs_L4[1][2][5][8] = 0;
  Curvature_Higgs_L4[1][2][6][6] = 0;
  Curvature_Higgs_L4[1][2][6][7] = 0;
  Curvature_Higgs_L4[1][2][6][8] = 0;
  Curvature_Higgs_L4[1][2][7][7] = 0;
  Curvature_Higgs_L4[1][2][7][8] = 0;
  Curvature_Higgs_L4[1][2][8][8] = 0;
  Curvature_Higgs_L4[1][3][3][3] = 0;
  Curvature_Higgs_L4[1][3][3][4] = 0;
  Curvature_Higgs_L4[1][3][3][5] = 0;
  Curvature_Higgs_L4[1][3][3][6] = 0;
  Curvature_Higgs_L4[1][3][3][7] = 0;
  Curvature_Higgs_L4[1][3][3][8] = 0;
  Curvature_Higgs_L4[1][3][4][4] = 0;
  Curvature_Higgs_L4[1][3][4][5] = 0;
  Curvature_Higgs_L4[1][3][4][6] = 0;
  Curvature_Higgs_L4[1][3][4][7] = 0;
  Curvature_Higgs_L4[1][3][4][8] = 0;
  Curvature_Higgs_L4[1][3][5][5] = 0;
  Curvature_Higgs_L4[1][3][5][6] = 0;
  Curvature_Higgs_L4[1][3][5][7] = 0;
  Curvature_Higgs_L4[1][3][5][8] = 0;
  Curvature_Higgs_L4[1][3][6][6] = 0;
  Curvature_Higgs_L4[1][3][6][7] = 0;
  Curvature_Higgs_L4[1][3][6][8] = 0;
  Curvature_Higgs_L4[1][3][7][7] = 0;
  Curvature_Higgs_L4[1][3][7][8] = 0;
  Curvature_Higgs_L4[1][3][8][8] = 0;
  Curvature_Higgs_L4[1][4][4][4] = 0;
  Curvature_Higgs_L4[1][4][4][5] = 0;
  Curvature_Higgs_L4[1][4][4][6] = 0;
  Curvature_Higgs_L4[1][4][4][7] = 0;
  Curvature_Higgs_L4[1][4][4][8] = 0;
  Curvature_Higgs_L4[1][4][5][5] = 0;
  Curvature_Higgs_L4[1][4][5][6] = 0;
  Curvature_Higgs_L4[1][4][5][7] = 0;
  Curvature_Higgs_L4[1][4][5][8] = 0;
  Curvature_Higgs_L4[1][4][6][6] = 0;
  Curvature_Higgs_L4[1][4][6][7] = 0;
  Curvature_Higgs_L4[1][4][6][8] = 0;
  Curvature_Higgs_L4[1][4][7][7] = 0;
  Curvature_Higgs_L4[1][4][7][8] = 0;
  Curvature_Higgs_L4[1][4][8][8] = 0;
  Curvature_Higgs_L4[1][5][5][5] = 0;
  Curvature_Higgs_L4[1][5][5][6] = 0;
  Curvature_Higgs_L4[1][5][5][7] = 0;
  Curvature_Higgs_L4[1][5][5][8] = 0;
  Curvature_Higgs_L4[1][5][6][6] = 0;
  Curvature_Higgs_L4[1][5][6][7] = 0;
  Curvature_Higgs_L4[1][5][6][8] = 0;
  Curvature_Higgs_L4[1][5][7][7] = 0;
  Curvature_Higgs_L4[1][5][7][8] = 0;
  Curvature_Higgs_L4[1][5][8][8] = 0;
  Curvature_Higgs_L4[1][6][6][6] = 0;
  Curvature_Higgs_L4[1][6][6][7] = 0;
  Curvature_Higgs_L4[1][6][6][8] = 0;
  Curvature_Higgs_L4[1][6][7][7] = 0;
  Curvature_Higgs_L4[1][6][7][8] = 0;
  Curvature_Higgs_L4[1][6][8][8] = 0;
  Curvature_Higgs_L4[1][7][7][7] = 0;
  Curvature_Higgs_L4[1][7][7][8] = 0;
  Curvature_Higgs_L4[1][7][8][8] = 0;
  Curvature_Higgs_L4[1][8][8][8] = 0;
  Curvature_Higgs_L4[2][2][2][2] = 3 * L1;
  Curvature_Higgs_L4[2][2][2][3] = 0;
  Curvature_Higgs_L4[2][2][2][4] = 0;
  Curvature_Higgs_L4[2][2][2][5] = 0;
  Curvature_Higgs_L4[2][2][2][6] = 0;
  Curvature_Higgs_L4[2][2][2][7] = 0;
  Curvature_Higgs_L4[2][2][2][8] = 0;
  Curvature_Higgs_L4[2][2][3][3] = L3 + L4 + RL5;
  Curvature_Higgs_L4[2][2][3][4] = 0;
  Curvature_Higgs_L4[2][2][3][5] = 0;
  Curvature_Higgs_L4[2][2][3][6] = 0;
  Curvature_Higgs_L4[2][2][3][7] = 0;
  Curvature_Higgs_L4[2][2][3][8] = 0;
  Curvature_Higgs_L4[2][2][4][4] = L1;
  Curvature_Higgs_L4[2][2][4][5] = 0;
  Curvature_Higgs_L4[2][2][4][6] = 0;
  Curvature_Higgs_L4[2][2][4][7] = 0;
  Curvature_Higgs_L4[2][2][4][8] = 0;
  Curvature_Higgs_L4[2][2][5][5] = L3;
  Curvature_Higgs_L4[2][2][5][6] = 0;
  Curvature_Higgs_L4[2][2][5][7] = 0;
  Curvature_Higgs_L4[2][2][5][8] = 0;
  Curvature_Higgs_L4[2][2][6][6] = L1;
  Curvature_Higgs_L4[2][2][6][7] = 0;
  Curvature_Higgs_L4[2][2][6][8] = 0;
  Curvature_Higgs_L4[2][2][7][7] = L3;
  Curvature_Higgs_L4[2][2][7][8] = 0;
  Curvature_Higgs_L4[2][2][8][8] = NL7;
  Curvature_Higgs_L4[2][3][3][3] = 0;
  Curvature_Higgs_L4[2][3][3][4] = 0;
  Curvature_Higgs_L4[2][3][3][5] = 0;
  Curvature_Higgs_L4[2][3][3][6] = 0;
  Curvature_Higgs_L4[2][3][3][7] = 0;
  Curvature_Higgs_L4[2][3][3][8] = 0;
  Curvature_Higgs_L4[2][3][4][4] = 0;
  Curvature_Higgs_L4[2][3][4][5] = L4 / 0.2e1 + RL5 / 0.2e1;
  Curvature_Higgs_L4[2][3][4][6] = 0;
  Curvature_Higgs_L4[2][3][4][7] = 0;
  Curvature_Higgs_L4[2][3][4][8] = 0;
  Curvature_Higgs_L4[2][3][5][5] = 0;
  Curvature_Higgs_L4[2][3][5][6] = 0;
  Curvature_Higgs_L4[2][3][5][7] = 0;
  Curvature_Higgs_L4[2][3][5][8] = 0;
  Curvature_Higgs_L4[2][3][6][6] = 0;
  Curvature_Higgs_L4[2][3][6][7] = L4 / 0.2e1 + RL5 / 0.2e1;
  Curvature_Higgs_L4[2][3][6][8] = 0;
  Curvature_Higgs_L4[2][3][7][7] = 0;
  Curvature_Higgs_L4[2][3][7][8] = 0;
  Curvature_Higgs_L4[2][3][8][8] = 0;
  Curvature_Higgs_L4[2][4][4][4] = 0;
  Curvature_Higgs_L4[2][4][4][5] = 0;
  Curvature_Higgs_L4[2][4][4][6] = 0;
  Curvature_Higgs_L4[2][4][4][7] = 0;
  Curvature_Higgs_L4[2][4][4][8] = 0;
  Curvature_Higgs_L4[2][4][5][5] = 0;
  Curvature_Higgs_L4[2][4][5][6] = 0;
  Curvature_Higgs_L4[2][4][5][7] = 0;
  Curvature_Higgs_L4[2][4][5][8] = 0;
  Curvature_Higgs_L4[2][4][6][6] = 0;
  Curvature_Higgs_L4[2][4][6][7] = 0;
  Curvature_Higgs_L4[2][4][6][8] = 0;
  Curvature_Higgs_L4[2][4][7][7] = 0;
  Curvature_Higgs_L4[2][4][7][8] = 0;
  Curvature_Higgs_L4[2][4][8][8] = 0;
  Curvature_Higgs_L4[2][5][5][5] = 0;
  Curvature_Higgs_L4[2][5][5][6] = 0;
  Curvature_Higgs_L4[2][5][5][7] = 0;
  Curvature_Higgs_L4[2][5][5][8] = 0;
  Curvature_Higgs_L4[2][5][6][6] = 0;
  Curvature_Higgs_L4[2][5][6][7] = 0;
  Curvature_Higgs_L4[2][5][6][8] = 0;
  Curvature_Higgs_L4[2][5][7][7] = 0;
  Curvature_Higgs_L4[2][5][7][8] = 0;
  Curvature_Higgs_L4[2][5][8][8] = 0;
  Curvature_Higgs_L4[2][6][6][6] = 0;
  Curvature_Higgs_L4[2][6][6][7] = 0;
  Curvature_Higgs_L4[2][6][6][8] = 0;
  Curvature_Higgs_L4[2][6][7][7] = 0;
  Curvature_Higgs_L4[2][6][7][8] = 0;
  Curvature_Higgs_L4[2][6][8][8] = 0;
  Curvature_Higgs_L4[2][7][7][7] = 0;
  Curvature_Higgs_L4[2][7][7][8] = 0;
  Curvature_Higgs_L4[2][7][8][8] = 0;
  Curvature_Higgs_L4[2][8][8][8] = 0;
  Curvature_Higgs_L4[3][3][3][3] = 3 * L2;
  Curvature_Higgs_L4[3][3][3][4] = 0;
  Curvature_Higgs_L4[3][3][3][5] = 0;
  Curvature_Higgs_L4[3][3][3][6] = 0;
  Curvature_Higgs_L4[3][3][3][7] = 0;
  Curvature_Higgs_L4[3][3][3][8] = 0;
  Curvature_Higgs_L4[3][3][4][4] = L3;
  Curvature_Higgs_L4[3][3][4][5] = 0;
  Curvature_Higgs_L4[3][3][4][6] = 0;
  Curvature_Higgs_L4[3][3][4][7] = 0;
  Curvature_Higgs_L4[3][3][4][8] = 0;
  Curvature_Higgs_L4[3][3][5][5] = L2;
  Curvature_Higgs_L4[3][3][5][6] = 0;
  Curvature_Higgs_L4[3][3][5][7] = 0;
  Curvature_Higgs_L4[3][3][5][8] = 0;
  Curvature_Higgs_L4[3][3][6][6] = L3;
  Curvature_Higgs_L4[3][3][6][7] = 0;
  Curvature_Higgs_L4[3][3][6][8] = 0;
  Curvature_Higgs_L4[3][3][7][7] = L2;
  Curvature_Higgs_L4[3][3][7][8] = 0;
  Curvature_Higgs_L4[3][3][8][8] = NL8;
  Curvature_Higgs_L4[3][4][4][4] = 0;
  Curvature_Higgs_L4[3][4][4][5] = 0;
  Curvature_Higgs_L4[3][4][4][6] = 0;
  Curvature_Higgs_L4[3][4][4][7] = 0;
  Curvature_Higgs_L4[3][4][4][8] = 0;
  Curvature_Higgs_L4[3][4][5][5] = 0;
  Curvature_Higgs_L4[3][4][5][6] = 0;
  Curvature_Higgs_L4[3][4][5][7] = 0;
  Curvature_Higgs_L4[3][4][5][8] = 0;
  Curvature_Higgs_L4[3][4][6][6] = 0;
  Curvature_Higgs_L4[3][4][6][7] = 0;
  Curvature_Higgs_L4[3][4][6][8] = 0;
  Curvature_Higgs_L4[3][4][7][7] = 0;
  Curvature_Higgs_L4[3][4][7][8] = 0;
  Curvature_Higgs_L4[3][4][8][8] = 0;
  Curvature_Higgs_L4[3][5][5][5] = 0;
  Curvature_Higgs_L4[3][5][5][6] = 0;
  Curvature_Higgs_L4[3][5][5][7] = 0;
  Curvature_Higgs_L4[3][5][5][8] = 0;
  Curvature_Higgs_L4[3][5][6][6] = 0;
  Curvature_Higgs_L4[3][5][6][7] = 0;
  Curvature_Higgs_L4[3][5][6][8] = 0;
  Curvature_Higgs_L4[3][5][7][7] = 0;
  Curvature_Higgs_L4[3][5][7][8] = 0;
  Curvature_Higgs_L4[3][5][8][8] = 0;
  Curvature_Higgs_L4[3][6][6][6] = 0;
  Curvature_Higgs_L4[3][6][6][7] = 0;
  Curvature_Higgs_L4[3][6][6][8] = 0;
  Curvature_Higgs_L4[3][6][7][7] = 0;
  Curvature_Higgs_L4[3][6][7][8] = 0;
  Curvature_Higgs_L4[3][6][8][8] = 0;
  Curvature_Higgs_L4[3][7][7][7] = 0;
  Curvature_Higgs_L4[3][7][7][8] = 0;
  Curvature_Higgs_L4[3][7][8][8] = 0;
  Curvature_Higgs_L4[3][8][8][8] = 0;
  Curvature_Higgs_L4[4][4][4][4] = 3 * L1;
  Curvature_Higgs_L4[4][4][4][5] = 0;
  Curvature_Higgs_L4[4][4][4][6] = 0;
  Curvature_Higgs_L4[4][4][4][7] = 0;
  Curvature_Higgs_L4[4][4][4][8] = 0;
  Curvature_Higgs_L4[4][4][5][5] = L3 + L4 + RL5;
  Curvature_Higgs_L4[4][4][5][6] = 0;
  Curvature_Higgs_L4[4][4][5][7] = 0;
  Curvature_Higgs_L4[4][4][5][8] = 0;
  Curvature_Higgs_L4[4][4][6][6] = L1;
  Curvature_Higgs_L4[4][4][6][7] = 0;
  Curvature_Higgs_L4[4][4][6][8] = 0;
  Curvature_Higgs_L4[4][4][7][7] = L3 + L4 - RL5;
  Curvature_Higgs_L4[4][4][7][8] = 0;
  Curvature_Higgs_L4[4][4][8][8] = NL7;
  Curvature_Higgs_L4[4][5][5][5] = 0;
  Curvature_Higgs_L4[4][5][5][6] = 0;
  Curvature_Higgs_L4[4][5][5][7] = 0;
  Curvature_Higgs_L4[4][5][5][8] = 0;
  Curvature_Higgs_L4[4][5][6][6] = 0;
  Curvature_Higgs_L4[4][5][6][7] = RL5;
  Curvature_Higgs_L4[4][5][6][8] = 0;
  Curvature_Higgs_L4[4][5][7][7] = 0;
  Curvature_Higgs_L4[4][5][7][8] = 0;
  Curvature_Higgs_L4[4][5][8][8] = 0;
  Curvature_Higgs_L4[4][6][6][6] = 0;
  Curvature_Higgs_L4[4][6][6][7] = 0;
  Curvature_Higgs_L4[4][6][6][8] = 0;
  Curvature_Higgs_L4[4][6][7][7] = 0;
  Curvature_Higgs_L4[4][6][7][8] = 0;
  Curvature_Higgs_L4[4][6][8][8] = 0;
  Curvature_Higgs_L4[4][7][7][7] = 0;
  Curvature_Higgs_L4[4][7][7][8] = 0;
  Curvature_Higgs_L4[4][7][8][8] = 0;
  Curvature_Higgs_L4[4][8][8][8] = 0;
  Curvature_Higgs_L4[5][5][5][5] = 3 * L2;
  Curvature_Higgs_L4[5][5][5][6] = 0;
  Curvature_Higgs_L4[5][5][5][7] = 0;
  Curvature_Higgs_L4[5][5][5][8] = 0;
  Curvature_Higgs_L4[5][5][6][6] = L3 + L4 - RL5;
  Curvature_Higgs_L4[5][5][6][7] = 0;
  Curvature_Higgs_L4[5][5][6][8] = 0;
  Curvature_Higgs_L4[5][5][7][7] = L2;
  Curvature_Higgs_L4[5][5][7][8] = 0;
  Curvature_Higgs_L4[5][5][8][8] = NL8;
  Curvature_Higgs_L4[5][6][6][6] = 0;
  Curvature_Higgs_L4[5][6][6][7] = 0;
  Curvature_Higgs_L4[5][6][6][8] = 0;
  Curvature_Higgs_L4[5][6][7][7] = 0;
  Curvature_Higgs_L4[5][6][7][8] = 0;
  Curvature_Higgs_L4[5][6][8][8] = 0;
  Curvature_Higgs_L4[5][7][7][7] = 0;
  Curvature_Higgs_L4[5][7][7][8] = 0;
  Curvature_Higgs_L4[5][7][8][8] = 0;
  Curvature_Higgs_L4[5][8][8][8] = 0;
  Curvature_Higgs_L4[6][6][6][6] = 3 * L1;
  Curvature_Higgs_L4[6][6][6][7] = 0;
  Curvature_Higgs_L4[6][6][6][8] = 0;
  Curvature_Higgs_L4[6][6][7][7] = L3 + L4 + RL5;
  Curvature_Higgs_L4[6][6][7][8] = 0;
  Curvature_Higgs_L4[6][6][8][8] = NL7;
  Curvature_Higgs_L4[6][7][7][7] = 0;
  Curvature_Higgs_L4[6][7][7][8] = 0;
  Curvature_Higgs_L4[6][7][8][8] = 0;
  Curvature_Higgs_L4[6][8][8][8] = 0;
  Curvature_Higgs_L4[7][7][7][7] = 3 * L2;
  Curvature_Higgs_L4[7][7][7][8] = 0;
  Curvature_Higgs_L4[7][7][8][8] = NL8;
  Curvature_Higgs_L4[7][8][8][8] = 0;
  Curvature_Higgs_L4[8][8][8][8] = 3 * NL6;

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
  Curvature_Gauge_G2H2[0][0][0][1] = 0;
  Curvature_Gauge_G2H2[0][0][0][2] = 0;
  Curvature_Gauge_G2H2[0][0][0][3] = 0;
  Curvature_Gauge_G2H2[0][0][0][4] = 0;
  Curvature_Gauge_G2H2[0][0][0][5] = 0;
  Curvature_Gauge_G2H2[0][0][0][6] = 0;
  Curvature_Gauge_G2H2[0][0][0][7] = 0;
  Curvature_Gauge_G2H2[0][0][0][8] = 0;
  Curvature_Gauge_G2H2[0][0][1][0] = 0;
  Curvature_Gauge_G2H2[0][0][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][1][2] = 0;
  Curvature_Gauge_G2H2[0][0][1][3] = 0;
  Curvature_Gauge_G2H2[0][0][1][4] = 0;
  Curvature_Gauge_G2H2[0][0][1][5] = 0;
  Curvature_Gauge_G2H2[0][0][1][6] = 0;
  Curvature_Gauge_G2H2[0][0][1][7] = 0;
  Curvature_Gauge_G2H2[0][0][1][8] = 0;
  Curvature_Gauge_G2H2[0][0][2][0] = 0;
  Curvature_Gauge_G2H2[0][0][2][1] = 0;
  Curvature_Gauge_G2H2[0][0][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][2][3] = 0;
  Curvature_Gauge_G2H2[0][0][2][4] = 0;
  Curvature_Gauge_G2H2[0][0][2][5] = 0;
  Curvature_Gauge_G2H2[0][0][2][6] = 0;
  Curvature_Gauge_G2H2[0][0][2][7] = 0;
  Curvature_Gauge_G2H2[0][0][2][8] = 0;
  Curvature_Gauge_G2H2[0][0][3][0] = 0;
  Curvature_Gauge_G2H2[0][0][3][1] = 0;
  Curvature_Gauge_G2H2[0][0][3][2] = 0;
  Curvature_Gauge_G2H2[0][0][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][3][4] = 0;
  Curvature_Gauge_G2H2[0][0][3][5] = 0;
  Curvature_Gauge_G2H2[0][0][3][6] = 0;
  Curvature_Gauge_G2H2[0][0][3][7] = 0;
  Curvature_Gauge_G2H2[0][0][3][8] = 0;
  Curvature_Gauge_G2H2[0][0][4][0] = 0;
  Curvature_Gauge_G2H2[0][0][4][1] = 0;
  Curvature_Gauge_G2H2[0][0][4][2] = 0;
  Curvature_Gauge_G2H2[0][0][4][3] = 0;
  Curvature_Gauge_G2H2[0][0][4][4] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][4][5] = 0;
  Curvature_Gauge_G2H2[0][0][4][6] = 0;
  Curvature_Gauge_G2H2[0][0][4][7] = 0;
  Curvature_Gauge_G2H2[0][0][4][8] = 0;
  Curvature_Gauge_G2H2[0][0][5][0] = 0;
  Curvature_Gauge_G2H2[0][0][5][1] = 0;
  Curvature_Gauge_G2H2[0][0][5][2] = 0;
  Curvature_Gauge_G2H2[0][0][5][3] = 0;
  Curvature_Gauge_G2H2[0][0][5][4] = 0;
  Curvature_Gauge_G2H2[0][0][5][5] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][5][6] = 0;
  Curvature_Gauge_G2H2[0][0][5][7] = 0;
  Curvature_Gauge_G2H2[0][0][5][8] = 0;
  Curvature_Gauge_G2H2[0][0][6][0] = 0;
  Curvature_Gauge_G2H2[0][0][6][1] = 0;
  Curvature_Gauge_G2H2[0][0][6][2] = 0;
  Curvature_Gauge_G2H2[0][0][6][3] = 0;
  Curvature_Gauge_G2H2[0][0][6][4] = 0;
  Curvature_Gauge_G2H2[0][0][6][5] = 0;
  Curvature_Gauge_G2H2[0][0][6][6] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][6][7] = 0;
  Curvature_Gauge_G2H2[0][0][6][8] = 0;
  Curvature_Gauge_G2H2[0][0][7][0] = 0;
  Curvature_Gauge_G2H2[0][0][7][1] = 0;
  Curvature_Gauge_G2H2[0][0][7][2] = 0;
  Curvature_Gauge_G2H2[0][0][7][3] = 0;
  Curvature_Gauge_G2H2[0][0][7][4] = 0;
  Curvature_Gauge_G2H2[0][0][7][5] = 0;
  Curvature_Gauge_G2H2[0][0][7][6] = 0;
  Curvature_Gauge_G2H2[0][0][7][7] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][7][8] = 0;
  Curvature_Gauge_G2H2[0][0][8][0] = 0;
  Curvature_Gauge_G2H2[0][0][8][1] = 0;
  Curvature_Gauge_G2H2[0][0][8][2] = 0;
  Curvature_Gauge_G2H2[0][0][8][3] = 0;
  Curvature_Gauge_G2H2[0][0][8][4] = 0;
  Curvature_Gauge_G2H2[0][0][8][5] = 0;
  Curvature_Gauge_G2H2[0][0][8][6] = 0;
  Curvature_Gauge_G2H2[0][0][8][7] = 0;
  Curvature_Gauge_G2H2[0][0][8][8] = 0;
  Curvature_Gauge_G2H2[0][1][0][0] = 0;
  Curvature_Gauge_G2H2[0][1][0][1] = 0;
  Curvature_Gauge_G2H2[0][1][0][2] = 0;
  Curvature_Gauge_G2H2[0][1][0][3] = 0;
  Curvature_Gauge_G2H2[0][1][0][4] = 0;
  Curvature_Gauge_G2H2[0][1][0][5] = 0;
  Curvature_Gauge_G2H2[0][1][0][6] = 0;
  Curvature_Gauge_G2H2[0][1][0][7] = 0;
  Curvature_Gauge_G2H2[0][1][0][8] = 0;
  Curvature_Gauge_G2H2[0][1][1][0] = 0;
  Curvature_Gauge_G2H2[0][1][1][1] = 0;
  Curvature_Gauge_G2H2[0][1][1][2] = 0;
  Curvature_Gauge_G2H2[0][1][1][3] = 0;
  Curvature_Gauge_G2H2[0][1][1][4] = 0;
  Curvature_Gauge_G2H2[0][1][1][5] = 0;
  Curvature_Gauge_G2H2[0][1][1][6] = 0;
  Curvature_Gauge_G2H2[0][1][1][7] = 0;
  Curvature_Gauge_G2H2[0][1][1][8] = 0;
  Curvature_Gauge_G2H2[0][1][2][0] = 0;
  Curvature_Gauge_G2H2[0][1][2][1] = 0;
  Curvature_Gauge_G2H2[0][1][2][2] = 0;
  Curvature_Gauge_G2H2[0][1][2][3] = 0;
  Curvature_Gauge_G2H2[0][1][2][4] = 0;
  Curvature_Gauge_G2H2[0][1][2][5] = 0;
  Curvature_Gauge_G2H2[0][1][2][6] = 0;
  Curvature_Gauge_G2H2[0][1][2][7] = 0;
  Curvature_Gauge_G2H2[0][1][2][8] = 0;
  Curvature_Gauge_G2H2[0][1][3][0] = 0;
  Curvature_Gauge_G2H2[0][1][3][1] = 0;
  Curvature_Gauge_G2H2[0][1][3][2] = 0;
  Curvature_Gauge_G2H2[0][1][3][3] = 0;
  Curvature_Gauge_G2H2[0][1][3][4] = 0;
  Curvature_Gauge_G2H2[0][1][3][5] = 0;
  Curvature_Gauge_G2H2[0][1][3][6] = 0;
  Curvature_Gauge_G2H2[0][1][3][7] = 0;
  Curvature_Gauge_G2H2[0][1][3][8] = 0;
  Curvature_Gauge_G2H2[0][1][4][0] = 0;
  Curvature_Gauge_G2H2[0][1][4][1] = 0;
  Curvature_Gauge_G2H2[0][1][4][2] = 0;
  Curvature_Gauge_G2H2[0][1][4][3] = 0;
  Curvature_Gauge_G2H2[0][1][4][4] = 0;
  Curvature_Gauge_G2H2[0][1][4][5] = 0;
  Curvature_Gauge_G2H2[0][1][4][6] = 0;
  Curvature_Gauge_G2H2[0][1][4][7] = 0;
  Curvature_Gauge_G2H2[0][1][4][8] = 0;
  Curvature_Gauge_G2H2[0][1][5][0] = 0;
  Curvature_Gauge_G2H2[0][1][5][1] = 0;
  Curvature_Gauge_G2H2[0][1][5][2] = 0;
  Curvature_Gauge_G2H2[0][1][5][3] = 0;
  Curvature_Gauge_G2H2[0][1][5][4] = 0;
  Curvature_Gauge_G2H2[0][1][5][5] = 0;
  Curvature_Gauge_G2H2[0][1][5][6] = 0;
  Curvature_Gauge_G2H2[0][1][5][7] = 0;
  Curvature_Gauge_G2H2[0][1][5][8] = 0;
  Curvature_Gauge_G2H2[0][1][6][0] = 0;
  Curvature_Gauge_G2H2[0][1][6][1] = 0;
  Curvature_Gauge_G2H2[0][1][6][2] = 0;
  Curvature_Gauge_G2H2[0][1][6][3] = 0;
  Curvature_Gauge_G2H2[0][1][6][4] = 0;
  Curvature_Gauge_G2H2[0][1][6][5] = 0;
  Curvature_Gauge_G2H2[0][1][6][6] = 0;
  Curvature_Gauge_G2H2[0][1][6][7] = 0;
  Curvature_Gauge_G2H2[0][1][6][8] = 0;
  Curvature_Gauge_G2H2[0][1][7][0] = 0;
  Curvature_Gauge_G2H2[0][1][7][1] = 0;
  Curvature_Gauge_G2H2[0][1][7][2] = 0;
  Curvature_Gauge_G2H2[0][1][7][3] = 0;
  Curvature_Gauge_G2H2[0][1][7][4] = 0;
  Curvature_Gauge_G2H2[0][1][7][5] = 0;
  Curvature_Gauge_G2H2[0][1][7][6] = 0;
  Curvature_Gauge_G2H2[0][1][7][7] = 0;
  Curvature_Gauge_G2H2[0][1][7][8] = 0;
  Curvature_Gauge_G2H2[0][1][8][0] = 0;
  Curvature_Gauge_G2H2[0][1][8][1] = 0;
  Curvature_Gauge_G2H2[0][1][8][2] = 0;
  Curvature_Gauge_G2H2[0][1][8][3] = 0;
  Curvature_Gauge_G2H2[0][1][8][4] = 0;
  Curvature_Gauge_G2H2[0][1][8][5] = 0;
  Curvature_Gauge_G2H2[0][1][8][6] = 0;
  Curvature_Gauge_G2H2[0][1][8][7] = 0;
  Curvature_Gauge_G2H2[0][1][8][8] = 0;
  Curvature_Gauge_G2H2[0][2][0][0] = 0;
  Curvature_Gauge_G2H2[0][2][0][1] = 0;
  Curvature_Gauge_G2H2[0][2][0][2] = 0;
  Curvature_Gauge_G2H2[0][2][0][3] = 0;
  Curvature_Gauge_G2H2[0][2][0][4] = 0;
  Curvature_Gauge_G2H2[0][2][0][5] = 0;
  Curvature_Gauge_G2H2[0][2][0][6] = 0;
  Curvature_Gauge_G2H2[0][2][0][7] = 0;
  Curvature_Gauge_G2H2[0][2][0][8] = 0;
  Curvature_Gauge_G2H2[0][2][1][0] = 0;
  Curvature_Gauge_G2H2[0][2][1][1] = 0;
  Curvature_Gauge_G2H2[0][2][1][2] = 0;
  Curvature_Gauge_G2H2[0][2][1][3] = 0;
  Curvature_Gauge_G2H2[0][2][1][4] = 0;
  Curvature_Gauge_G2H2[0][2][1][5] = 0;
  Curvature_Gauge_G2H2[0][2][1][6] = 0;
  Curvature_Gauge_G2H2[0][2][1][7] = 0;
  Curvature_Gauge_G2H2[0][2][1][8] = 0;
  Curvature_Gauge_G2H2[0][2][2][0] = 0;
  Curvature_Gauge_G2H2[0][2][2][1] = 0;
  Curvature_Gauge_G2H2[0][2][2][2] = 0;
  Curvature_Gauge_G2H2[0][2][2][3] = 0;
  Curvature_Gauge_G2H2[0][2][2][4] = 0;
  Curvature_Gauge_G2H2[0][2][2][5] = 0;
  Curvature_Gauge_G2H2[0][2][2][6] = 0;
  Curvature_Gauge_G2H2[0][2][2][7] = 0;
  Curvature_Gauge_G2H2[0][2][2][8] = 0;
  Curvature_Gauge_G2H2[0][2][3][0] = 0;
  Curvature_Gauge_G2H2[0][2][3][1] = 0;
  Curvature_Gauge_G2H2[0][2][3][2] = 0;
  Curvature_Gauge_G2H2[0][2][3][3] = 0;
  Curvature_Gauge_G2H2[0][2][3][4] = 0;
  Curvature_Gauge_G2H2[0][2][3][5] = 0;
  Curvature_Gauge_G2H2[0][2][3][6] = 0;
  Curvature_Gauge_G2H2[0][2][3][7] = 0;
  Curvature_Gauge_G2H2[0][2][3][8] = 0;
  Curvature_Gauge_G2H2[0][2][4][0] = 0;
  Curvature_Gauge_G2H2[0][2][4][1] = 0;
  Curvature_Gauge_G2H2[0][2][4][2] = 0;
  Curvature_Gauge_G2H2[0][2][4][3] = 0;
  Curvature_Gauge_G2H2[0][2][4][4] = 0;
  Curvature_Gauge_G2H2[0][2][4][5] = 0;
  Curvature_Gauge_G2H2[0][2][4][6] = 0;
  Curvature_Gauge_G2H2[0][2][4][7] = 0;
  Curvature_Gauge_G2H2[0][2][4][8] = 0;
  Curvature_Gauge_G2H2[0][2][5][0] = 0;
  Curvature_Gauge_G2H2[0][2][5][1] = 0;
  Curvature_Gauge_G2H2[0][2][5][2] = 0;
  Curvature_Gauge_G2H2[0][2][5][3] = 0;
  Curvature_Gauge_G2H2[0][2][5][4] = 0;
  Curvature_Gauge_G2H2[0][2][5][5] = 0;
  Curvature_Gauge_G2H2[0][2][5][6] = 0;
  Curvature_Gauge_G2H2[0][2][5][7] = 0;
  Curvature_Gauge_G2H2[0][2][5][8] = 0;
  Curvature_Gauge_G2H2[0][2][6][0] = 0;
  Curvature_Gauge_G2H2[0][2][6][1] = 0;
  Curvature_Gauge_G2H2[0][2][6][2] = 0;
  Curvature_Gauge_G2H2[0][2][6][3] = 0;
  Curvature_Gauge_G2H2[0][2][6][4] = 0;
  Curvature_Gauge_G2H2[0][2][6][5] = 0;
  Curvature_Gauge_G2H2[0][2][6][6] = 0;
  Curvature_Gauge_G2H2[0][2][6][7] = 0;
  Curvature_Gauge_G2H2[0][2][6][8] = 0;
  Curvature_Gauge_G2H2[0][2][7][0] = 0;
  Curvature_Gauge_G2H2[0][2][7][1] = 0;
  Curvature_Gauge_G2H2[0][2][7][2] = 0;
  Curvature_Gauge_G2H2[0][2][7][3] = 0;
  Curvature_Gauge_G2H2[0][2][7][4] = 0;
  Curvature_Gauge_G2H2[0][2][7][5] = 0;
  Curvature_Gauge_G2H2[0][2][7][6] = 0;
  Curvature_Gauge_G2H2[0][2][7][7] = 0;
  Curvature_Gauge_G2H2[0][2][7][8] = 0;
  Curvature_Gauge_G2H2[0][2][8][0] = 0;
  Curvature_Gauge_G2H2[0][2][8][1] = 0;
  Curvature_Gauge_G2H2[0][2][8][2] = 0;
  Curvature_Gauge_G2H2[0][2][8][3] = 0;
  Curvature_Gauge_G2H2[0][2][8][4] = 0;
  Curvature_Gauge_G2H2[0][2][8][5] = 0;
  Curvature_Gauge_G2H2[0][2][8][6] = 0;
  Curvature_Gauge_G2H2[0][2][8][7] = 0;
  Curvature_Gauge_G2H2[0][2][8][8] = 0;
  Curvature_Gauge_G2H2[0][3][0][0] = 0;
  Curvature_Gauge_G2H2[0][3][0][1] = 0;
  Curvature_Gauge_G2H2[0][3][0][2] = 0;
  Curvature_Gauge_G2H2[0][3][0][3] = 0;
  Curvature_Gauge_G2H2[0][3][0][4] = 0;
  Curvature_Gauge_G2H2[0][3][0][5] = 0;
  Curvature_Gauge_G2H2[0][3][0][6] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][0][7] = 0;
  Curvature_Gauge_G2H2[0][3][0][8] = 0;
  Curvature_Gauge_G2H2[0][3][1][0] = 0;
  Curvature_Gauge_G2H2[0][3][1][1] = 0;
  Curvature_Gauge_G2H2[0][3][1][2] = 0;
  Curvature_Gauge_G2H2[0][3][1][3] = 0;
  Curvature_Gauge_G2H2[0][3][1][4] = 0;
  Curvature_Gauge_G2H2[0][3][1][5] = 0;
  Curvature_Gauge_G2H2[0][3][1][6] = 0;
  Curvature_Gauge_G2H2[0][3][1][7] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][1][8] = 0;
  Curvature_Gauge_G2H2[0][3][2][0] = 0;
  Curvature_Gauge_G2H2[0][3][2][1] = 0;
  Curvature_Gauge_G2H2[0][3][2][2] = 0;
  Curvature_Gauge_G2H2[0][3][2][3] = 0;
  Curvature_Gauge_G2H2[0][3][2][4] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][2][5] = 0;
  Curvature_Gauge_G2H2[0][3][2][6] = 0;
  Curvature_Gauge_G2H2[0][3][2][7] = 0;
  Curvature_Gauge_G2H2[0][3][2][8] = 0;
  Curvature_Gauge_G2H2[0][3][3][0] = 0;
  Curvature_Gauge_G2H2[0][3][3][1] = 0;
  Curvature_Gauge_G2H2[0][3][3][2] = 0;
  Curvature_Gauge_G2H2[0][3][3][3] = 0;
  Curvature_Gauge_G2H2[0][3][3][4] = 0;
  Curvature_Gauge_G2H2[0][3][3][5] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][3][6] = 0;
  Curvature_Gauge_G2H2[0][3][3][7] = 0;
  Curvature_Gauge_G2H2[0][3][3][8] = 0;
  Curvature_Gauge_G2H2[0][3][4][0] = 0;
  Curvature_Gauge_G2H2[0][3][4][1] = 0;
  Curvature_Gauge_G2H2[0][3][4][2] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][4][3] = 0;
  Curvature_Gauge_G2H2[0][3][4][4] = 0;
  Curvature_Gauge_G2H2[0][3][4][5] = 0;
  Curvature_Gauge_G2H2[0][3][4][6] = 0;
  Curvature_Gauge_G2H2[0][3][4][7] = 0;
  Curvature_Gauge_G2H2[0][3][4][8] = 0;
  Curvature_Gauge_G2H2[0][3][5][0] = 0;
  Curvature_Gauge_G2H2[0][3][5][1] = 0;
  Curvature_Gauge_G2H2[0][3][5][2] = 0;
  Curvature_Gauge_G2H2[0][3][5][3] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][5][4] = 0;
  Curvature_Gauge_G2H2[0][3][5][5] = 0;
  Curvature_Gauge_G2H2[0][3][5][6] = 0;
  Curvature_Gauge_G2H2[0][3][5][7] = 0;
  Curvature_Gauge_G2H2[0][3][5][8] = 0;
  Curvature_Gauge_G2H2[0][3][6][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][6][1] = 0;
  Curvature_Gauge_G2H2[0][3][6][2] = 0;
  Curvature_Gauge_G2H2[0][3][6][3] = 0;
  Curvature_Gauge_G2H2[0][3][6][4] = 0;
  Curvature_Gauge_G2H2[0][3][6][5] = 0;
  Curvature_Gauge_G2H2[0][3][6][6] = 0;
  Curvature_Gauge_G2H2[0][3][6][7] = 0;
  Curvature_Gauge_G2H2[0][3][6][8] = 0;
  Curvature_Gauge_G2H2[0][3][7][0] = 0;
  Curvature_Gauge_G2H2[0][3][7][1] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][7][2] = 0;
  Curvature_Gauge_G2H2[0][3][7][3] = 0;
  Curvature_Gauge_G2H2[0][3][7][4] = 0;
  Curvature_Gauge_G2H2[0][3][7][5] = 0;
  Curvature_Gauge_G2H2[0][3][7][6] = 0;
  Curvature_Gauge_G2H2[0][3][7][7] = 0;
  Curvature_Gauge_G2H2[0][3][7][8] = 0;
  Curvature_Gauge_G2H2[0][3][8][0] = 0;
  Curvature_Gauge_G2H2[0][3][8][1] = 0;
  Curvature_Gauge_G2H2[0][3][8][2] = 0;
  Curvature_Gauge_G2H2[0][3][8][3] = 0;
  Curvature_Gauge_G2H2[0][3][8][4] = 0;
  Curvature_Gauge_G2H2[0][3][8][5] = 0;
  Curvature_Gauge_G2H2[0][3][8][6] = 0;
  Curvature_Gauge_G2H2[0][3][8][7] = 0;
  Curvature_Gauge_G2H2[0][3][8][8] = 0;
  Curvature_Gauge_G2H2[1][0][0][0] = 0;
  Curvature_Gauge_G2H2[1][0][0][1] = 0;
  Curvature_Gauge_G2H2[1][0][0][2] = 0;
  Curvature_Gauge_G2H2[1][0][0][3] = 0;
  Curvature_Gauge_G2H2[1][0][0][4] = 0;
  Curvature_Gauge_G2H2[1][0][0][5] = 0;
  Curvature_Gauge_G2H2[1][0][0][6] = 0;
  Curvature_Gauge_G2H2[1][0][0][7] = 0;
  Curvature_Gauge_G2H2[1][0][0][8] = 0;
  Curvature_Gauge_G2H2[1][0][1][0] = 0;
  Curvature_Gauge_G2H2[1][0][1][1] = 0;
  Curvature_Gauge_G2H2[1][0][1][2] = 0;
  Curvature_Gauge_G2H2[1][0][1][3] = 0;
  Curvature_Gauge_G2H2[1][0][1][4] = 0;
  Curvature_Gauge_G2H2[1][0][1][5] = 0;
  Curvature_Gauge_G2H2[1][0][1][6] = 0;
  Curvature_Gauge_G2H2[1][0][1][7] = 0;
  Curvature_Gauge_G2H2[1][0][1][8] = 0;
  Curvature_Gauge_G2H2[1][0][2][0] = 0;
  Curvature_Gauge_G2H2[1][0][2][1] = 0;
  Curvature_Gauge_G2H2[1][0][2][2] = 0;
  Curvature_Gauge_G2H2[1][0][2][3] = 0;
  Curvature_Gauge_G2H2[1][0][2][4] = 0;
  Curvature_Gauge_G2H2[1][0][2][5] = 0;
  Curvature_Gauge_G2H2[1][0][2][6] = 0;
  Curvature_Gauge_G2H2[1][0][2][7] = 0;
  Curvature_Gauge_G2H2[1][0][2][8] = 0;
  Curvature_Gauge_G2H2[1][0][3][0] = 0;
  Curvature_Gauge_G2H2[1][0][3][1] = 0;
  Curvature_Gauge_G2H2[1][0][3][2] = 0;
  Curvature_Gauge_G2H2[1][0][3][3] = 0;
  Curvature_Gauge_G2H2[1][0][3][4] = 0;
  Curvature_Gauge_G2H2[1][0][3][5] = 0;
  Curvature_Gauge_G2H2[1][0][3][6] = 0;
  Curvature_Gauge_G2H2[1][0][3][7] = 0;
  Curvature_Gauge_G2H2[1][0][3][8] = 0;
  Curvature_Gauge_G2H2[1][0][4][0] = 0;
  Curvature_Gauge_G2H2[1][0][4][1] = 0;
  Curvature_Gauge_G2H2[1][0][4][2] = 0;
  Curvature_Gauge_G2H2[1][0][4][3] = 0;
  Curvature_Gauge_G2H2[1][0][4][4] = 0;
  Curvature_Gauge_G2H2[1][0][4][5] = 0;
  Curvature_Gauge_G2H2[1][0][4][6] = 0;
  Curvature_Gauge_G2H2[1][0][4][7] = 0;
  Curvature_Gauge_G2H2[1][0][4][8] = 0;
  Curvature_Gauge_G2H2[1][0][5][0] = 0;
  Curvature_Gauge_G2H2[1][0][5][1] = 0;
  Curvature_Gauge_G2H2[1][0][5][2] = 0;
  Curvature_Gauge_G2H2[1][0][5][3] = 0;
  Curvature_Gauge_G2H2[1][0][5][4] = 0;
  Curvature_Gauge_G2H2[1][0][5][5] = 0;
  Curvature_Gauge_G2H2[1][0][5][6] = 0;
  Curvature_Gauge_G2H2[1][0][5][7] = 0;
  Curvature_Gauge_G2H2[1][0][5][8] = 0;
  Curvature_Gauge_G2H2[1][0][6][0] = 0;
  Curvature_Gauge_G2H2[1][0][6][1] = 0;
  Curvature_Gauge_G2H2[1][0][6][2] = 0;
  Curvature_Gauge_G2H2[1][0][6][3] = 0;
  Curvature_Gauge_G2H2[1][0][6][4] = 0;
  Curvature_Gauge_G2H2[1][0][6][5] = 0;
  Curvature_Gauge_G2H2[1][0][6][6] = 0;
  Curvature_Gauge_G2H2[1][0][6][7] = 0;
  Curvature_Gauge_G2H2[1][0][6][8] = 0;
  Curvature_Gauge_G2H2[1][0][7][0] = 0;
  Curvature_Gauge_G2H2[1][0][7][1] = 0;
  Curvature_Gauge_G2H2[1][0][7][2] = 0;
  Curvature_Gauge_G2H2[1][0][7][3] = 0;
  Curvature_Gauge_G2H2[1][0][7][4] = 0;
  Curvature_Gauge_G2H2[1][0][7][5] = 0;
  Curvature_Gauge_G2H2[1][0][7][6] = 0;
  Curvature_Gauge_G2H2[1][0][7][7] = 0;
  Curvature_Gauge_G2H2[1][0][7][8] = 0;
  Curvature_Gauge_G2H2[1][0][8][0] = 0;
  Curvature_Gauge_G2H2[1][0][8][1] = 0;
  Curvature_Gauge_G2H2[1][0][8][2] = 0;
  Curvature_Gauge_G2H2[1][0][8][3] = 0;
  Curvature_Gauge_G2H2[1][0][8][4] = 0;
  Curvature_Gauge_G2H2[1][0][8][5] = 0;
  Curvature_Gauge_G2H2[1][0][8][6] = 0;
  Curvature_Gauge_G2H2[1][0][8][7] = 0;
  Curvature_Gauge_G2H2[1][0][8][8] = 0;
  Curvature_Gauge_G2H2[1][1][0][0] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][0][1] = 0;
  Curvature_Gauge_G2H2[1][1][0][2] = 0;
  Curvature_Gauge_G2H2[1][1][0][3] = 0;
  Curvature_Gauge_G2H2[1][1][0][4] = 0;
  Curvature_Gauge_G2H2[1][1][0][5] = 0;
  Curvature_Gauge_G2H2[1][1][0][6] = 0;
  Curvature_Gauge_G2H2[1][1][0][7] = 0;
  Curvature_Gauge_G2H2[1][1][0][8] = 0;
  Curvature_Gauge_G2H2[1][1][1][0] = 0;
  Curvature_Gauge_G2H2[1][1][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][1][2] = 0;
  Curvature_Gauge_G2H2[1][1][1][3] = 0;
  Curvature_Gauge_G2H2[1][1][1][4] = 0;
  Curvature_Gauge_G2H2[1][1][1][5] = 0;
  Curvature_Gauge_G2H2[1][1][1][6] = 0;
  Curvature_Gauge_G2H2[1][1][1][7] = 0;
  Curvature_Gauge_G2H2[1][1][1][8] = 0;
  Curvature_Gauge_G2H2[1][1][2][0] = 0;
  Curvature_Gauge_G2H2[1][1][2][1] = 0;
  Curvature_Gauge_G2H2[1][1][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][2][3] = 0;
  Curvature_Gauge_G2H2[1][1][2][4] = 0;
  Curvature_Gauge_G2H2[1][1][2][5] = 0;
  Curvature_Gauge_G2H2[1][1][2][6] = 0;
  Curvature_Gauge_G2H2[1][1][2][7] = 0;
  Curvature_Gauge_G2H2[1][1][2][8] = 0;
  Curvature_Gauge_G2H2[1][1][3][0] = 0;
  Curvature_Gauge_G2H2[1][1][3][1] = 0;
  Curvature_Gauge_G2H2[1][1][3][2] = 0;
  Curvature_Gauge_G2H2[1][1][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][3][4] = 0;
  Curvature_Gauge_G2H2[1][1][3][5] = 0;
  Curvature_Gauge_G2H2[1][1][3][6] = 0;
  Curvature_Gauge_G2H2[1][1][3][7] = 0;
  Curvature_Gauge_G2H2[1][1][3][8] = 0;
  Curvature_Gauge_G2H2[1][1][4][0] = 0;
  Curvature_Gauge_G2H2[1][1][4][1] = 0;
  Curvature_Gauge_G2H2[1][1][4][2] = 0;
  Curvature_Gauge_G2H2[1][1][4][3] = 0;
  Curvature_Gauge_G2H2[1][1][4][4] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][4][5] = 0;
  Curvature_Gauge_G2H2[1][1][4][6] = 0;
  Curvature_Gauge_G2H2[1][1][4][7] = 0;
  Curvature_Gauge_G2H2[1][1][4][8] = 0;
  Curvature_Gauge_G2H2[1][1][5][0] = 0;
  Curvature_Gauge_G2H2[1][1][5][1] = 0;
  Curvature_Gauge_G2H2[1][1][5][2] = 0;
  Curvature_Gauge_G2H2[1][1][5][3] = 0;
  Curvature_Gauge_G2H2[1][1][5][4] = 0;
  Curvature_Gauge_G2H2[1][1][5][5] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][5][6] = 0;
  Curvature_Gauge_G2H2[1][1][5][7] = 0;
  Curvature_Gauge_G2H2[1][1][5][8] = 0;
  Curvature_Gauge_G2H2[1][1][6][0] = 0;
  Curvature_Gauge_G2H2[1][1][6][1] = 0;
  Curvature_Gauge_G2H2[1][1][6][2] = 0;
  Curvature_Gauge_G2H2[1][1][6][3] = 0;
  Curvature_Gauge_G2H2[1][1][6][4] = 0;
  Curvature_Gauge_G2H2[1][1][6][5] = 0;
  Curvature_Gauge_G2H2[1][1][6][6] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][6][7] = 0;
  Curvature_Gauge_G2H2[1][1][6][8] = 0;
  Curvature_Gauge_G2H2[1][1][7][0] = 0;
  Curvature_Gauge_G2H2[1][1][7][1] = 0;
  Curvature_Gauge_G2H2[1][1][7][2] = 0;
  Curvature_Gauge_G2H2[1][1][7][3] = 0;
  Curvature_Gauge_G2H2[1][1][7][4] = 0;
  Curvature_Gauge_G2H2[1][1][7][5] = 0;
  Curvature_Gauge_G2H2[1][1][7][6] = 0;
  Curvature_Gauge_G2H2[1][1][7][7] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][7][8] = 0;
  Curvature_Gauge_G2H2[1][1][8][0] = 0;
  Curvature_Gauge_G2H2[1][1][8][1] = 0;
  Curvature_Gauge_G2H2[1][1][8][2] = 0;
  Curvature_Gauge_G2H2[1][1][8][3] = 0;
  Curvature_Gauge_G2H2[1][1][8][4] = 0;
  Curvature_Gauge_G2H2[1][1][8][5] = 0;
  Curvature_Gauge_G2H2[1][1][8][6] = 0;
  Curvature_Gauge_G2H2[1][1][8][7] = 0;
  Curvature_Gauge_G2H2[1][1][8][8] = 0;
  Curvature_Gauge_G2H2[1][2][0][0] = 0;
  Curvature_Gauge_G2H2[1][2][0][1] = 0;
  Curvature_Gauge_G2H2[1][2][0][2] = 0;
  Curvature_Gauge_G2H2[1][2][0][3] = 0;
  Curvature_Gauge_G2H2[1][2][0][4] = 0;
  Curvature_Gauge_G2H2[1][2][0][5] = 0;
  Curvature_Gauge_G2H2[1][2][0][6] = 0;
  Curvature_Gauge_G2H2[1][2][0][7] = 0;
  Curvature_Gauge_G2H2[1][2][0][8] = 0;
  Curvature_Gauge_G2H2[1][2][1][0] = 0;
  Curvature_Gauge_G2H2[1][2][1][1] = 0;
  Curvature_Gauge_G2H2[1][2][1][2] = 0;
  Curvature_Gauge_G2H2[1][2][1][3] = 0;
  Curvature_Gauge_G2H2[1][2][1][4] = 0;
  Curvature_Gauge_G2H2[1][2][1][5] = 0;
  Curvature_Gauge_G2H2[1][2][1][6] = 0;
  Curvature_Gauge_G2H2[1][2][1][7] = 0;
  Curvature_Gauge_G2H2[1][2][1][8] = 0;
  Curvature_Gauge_G2H2[1][2][2][0] = 0;
  Curvature_Gauge_G2H2[1][2][2][1] = 0;
  Curvature_Gauge_G2H2[1][2][2][2] = 0;
  Curvature_Gauge_G2H2[1][2][2][3] = 0;
  Curvature_Gauge_G2H2[1][2][2][4] = 0;
  Curvature_Gauge_G2H2[1][2][2][5] = 0;
  Curvature_Gauge_G2H2[1][2][2][6] = 0;
  Curvature_Gauge_G2H2[1][2][2][7] = 0;
  Curvature_Gauge_G2H2[1][2][2][8] = 0;
  Curvature_Gauge_G2H2[1][2][3][0] = 0;
  Curvature_Gauge_G2H2[1][2][3][1] = 0;
  Curvature_Gauge_G2H2[1][2][3][2] = 0;
  Curvature_Gauge_G2H2[1][2][3][3] = 0;
  Curvature_Gauge_G2H2[1][2][3][4] = 0;
  Curvature_Gauge_G2H2[1][2][3][5] = 0;
  Curvature_Gauge_G2H2[1][2][3][6] = 0;
  Curvature_Gauge_G2H2[1][2][3][7] = 0;
  Curvature_Gauge_G2H2[1][2][3][8] = 0;
  Curvature_Gauge_G2H2[1][2][4][0] = 0;
  Curvature_Gauge_G2H2[1][2][4][1] = 0;
  Curvature_Gauge_G2H2[1][2][4][2] = 0;
  Curvature_Gauge_G2H2[1][2][4][3] = 0;
  Curvature_Gauge_G2H2[1][2][4][4] = 0;
  Curvature_Gauge_G2H2[1][2][4][5] = 0;
  Curvature_Gauge_G2H2[1][2][4][6] = 0;
  Curvature_Gauge_G2H2[1][2][4][7] = 0;
  Curvature_Gauge_G2H2[1][2][4][8] = 0;
  Curvature_Gauge_G2H2[1][2][5][0] = 0;
  Curvature_Gauge_G2H2[1][2][5][1] = 0;
  Curvature_Gauge_G2H2[1][2][5][2] = 0;
  Curvature_Gauge_G2H2[1][2][5][3] = 0;
  Curvature_Gauge_G2H2[1][2][5][4] = 0;
  Curvature_Gauge_G2H2[1][2][5][5] = 0;
  Curvature_Gauge_G2H2[1][2][5][6] = 0;
  Curvature_Gauge_G2H2[1][2][5][7] = 0;
  Curvature_Gauge_G2H2[1][2][5][8] = 0;
  Curvature_Gauge_G2H2[1][2][6][0] = 0;
  Curvature_Gauge_G2H2[1][2][6][1] = 0;
  Curvature_Gauge_G2H2[1][2][6][2] = 0;
  Curvature_Gauge_G2H2[1][2][6][3] = 0;
  Curvature_Gauge_G2H2[1][2][6][4] = 0;
  Curvature_Gauge_G2H2[1][2][6][5] = 0;
  Curvature_Gauge_G2H2[1][2][6][6] = 0;
  Curvature_Gauge_G2H2[1][2][6][7] = 0;
  Curvature_Gauge_G2H2[1][2][6][8] = 0;
  Curvature_Gauge_G2H2[1][2][7][0] = 0;
  Curvature_Gauge_G2H2[1][2][7][1] = 0;
  Curvature_Gauge_G2H2[1][2][7][2] = 0;
  Curvature_Gauge_G2H2[1][2][7][3] = 0;
  Curvature_Gauge_G2H2[1][2][7][4] = 0;
  Curvature_Gauge_G2H2[1][2][7][5] = 0;
  Curvature_Gauge_G2H2[1][2][7][6] = 0;
  Curvature_Gauge_G2H2[1][2][7][7] = 0;
  Curvature_Gauge_G2H2[1][2][7][8] = 0;
  Curvature_Gauge_G2H2[1][2][8][0] = 0;
  Curvature_Gauge_G2H2[1][2][8][1] = 0;
  Curvature_Gauge_G2H2[1][2][8][2] = 0;
  Curvature_Gauge_G2H2[1][2][8][3] = 0;
  Curvature_Gauge_G2H2[1][2][8][4] = 0;
  Curvature_Gauge_G2H2[1][2][8][5] = 0;
  Curvature_Gauge_G2H2[1][2][8][6] = 0;
  Curvature_Gauge_G2H2[1][2][8][7] = 0;
  Curvature_Gauge_G2H2[1][2][8][8] = 0;
  Curvature_Gauge_G2H2[1][3][0][0] = 0;
  Curvature_Gauge_G2H2[1][3][0][1] = 0;
  Curvature_Gauge_G2H2[1][3][0][2] = 0;
  Curvature_Gauge_G2H2[1][3][0][3] = 0;
  Curvature_Gauge_G2H2[1][3][0][4] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][0][5] = 0;
  Curvature_Gauge_G2H2[1][3][0][6] = 0;
  Curvature_Gauge_G2H2[1][3][0][7] = 0;
  Curvature_Gauge_G2H2[1][3][0][8] = 0;
  Curvature_Gauge_G2H2[1][3][1][0] = 0;
  Curvature_Gauge_G2H2[1][3][1][1] = 0;
  Curvature_Gauge_G2H2[1][3][1][2] = 0;
  Curvature_Gauge_G2H2[1][3][1][3] = 0;
  Curvature_Gauge_G2H2[1][3][1][4] = 0;
  Curvature_Gauge_G2H2[1][3][1][5] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][1][6] = 0;
  Curvature_Gauge_G2H2[1][3][1][7] = 0;
  Curvature_Gauge_G2H2[1][3][1][8] = 0;
  Curvature_Gauge_G2H2[1][3][2][0] = 0;
  Curvature_Gauge_G2H2[1][3][2][1] = 0;
  Curvature_Gauge_G2H2[1][3][2][2] = 0;
  Curvature_Gauge_G2H2[1][3][2][3] = 0;
  Curvature_Gauge_G2H2[1][3][2][4] = 0;
  Curvature_Gauge_G2H2[1][3][2][5] = 0;
  Curvature_Gauge_G2H2[1][3][2][6] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][2][7] = 0;
  Curvature_Gauge_G2H2[1][3][2][8] = 0;
  Curvature_Gauge_G2H2[1][3][3][0] = 0;
  Curvature_Gauge_G2H2[1][3][3][1] = 0;
  Curvature_Gauge_G2H2[1][3][3][2] = 0;
  Curvature_Gauge_G2H2[1][3][3][3] = 0;
  Curvature_Gauge_G2H2[1][3][3][4] = 0;
  Curvature_Gauge_G2H2[1][3][3][5] = 0;
  Curvature_Gauge_G2H2[1][3][3][6] = 0;
  Curvature_Gauge_G2H2[1][3][3][7] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][3][8] = 0;
  Curvature_Gauge_G2H2[1][3][4][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][4][1] = 0;
  Curvature_Gauge_G2H2[1][3][4][2] = 0;
  Curvature_Gauge_G2H2[1][3][4][3] = 0;
  Curvature_Gauge_G2H2[1][3][4][4] = 0;
  Curvature_Gauge_G2H2[1][3][4][5] = 0;
  Curvature_Gauge_G2H2[1][3][4][6] = 0;
  Curvature_Gauge_G2H2[1][3][4][7] = 0;
  Curvature_Gauge_G2H2[1][3][4][8] = 0;
  Curvature_Gauge_G2H2[1][3][5][0] = 0;
  Curvature_Gauge_G2H2[1][3][5][1] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][5][2] = 0;
  Curvature_Gauge_G2H2[1][3][5][3] = 0;
  Curvature_Gauge_G2H2[1][3][5][4] = 0;
  Curvature_Gauge_G2H2[1][3][5][5] = 0;
  Curvature_Gauge_G2H2[1][3][5][6] = 0;
  Curvature_Gauge_G2H2[1][3][5][7] = 0;
  Curvature_Gauge_G2H2[1][3][5][8] = 0;
  Curvature_Gauge_G2H2[1][3][6][0] = 0;
  Curvature_Gauge_G2H2[1][3][6][1] = 0;
  Curvature_Gauge_G2H2[1][3][6][2] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][6][3] = 0;
  Curvature_Gauge_G2H2[1][3][6][4] = 0;
  Curvature_Gauge_G2H2[1][3][6][5] = 0;
  Curvature_Gauge_G2H2[1][3][6][6] = 0;
  Curvature_Gauge_G2H2[1][3][6][7] = 0;
  Curvature_Gauge_G2H2[1][3][6][8] = 0;
  Curvature_Gauge_G2H2[1][3][7][0] = 0;
  Curvature_Gauge_G2H2[1][3][7][1] = 0;
  Curvature_Gauge_G2H2[1][3][7][2] = 0;
  Curvature_Gauge_G2H2[1][3][7][3] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][7][4] = 0;
  Curvature_Gauge_G2H2[1][3][7][5] = 0;
  Curvature_Gauge_G2H2[1][3][7][6] = 0;
  Curvature_Gauge_G2H2[1][3][7][7] = 0;
  Curvature_Gauge_G2H2[1][3][7][8] = 0;
  Curvature_Gauge_G2H2[1][3][8][0] = 0;
  Curvature_Gauge_G2H2[1][3][8][1] = 0;
  Curvature_Gauge_G2H2[1][3][8][2] = 0;
  Curvature_Gauge_G2H2[1][3][8][3] = 0;
  Curvature_Gauge_G2H2[1][3][8][4] = 0;
  Curvature_Gauge_G2H2[1][3][8][5] = 0;
  Curvature_Gauge_G2H2[1][3][8][6] = 0;
  Curvature_Gauge_G2H2[1][3][8][7] = 0;
  Curvature_Gauge_G2H2[1][3][8][8] = 0;
  Curvature_Gauge_G2H2[2][0][0][0] = 0;
  Curvature_Gauge_G2H2[2][0][0][1] = 0;
  Curvature_Gauge_G2H2[2][0][0][2] = 0;
  Curvature_Gauge_G2H2[2][0][0][3] = 0;
  Curvature_Gauge_G2H2[2][0][0][4] = 0;
  Curvature_Gauge_G2H2[2][0][0][5] = 0;
  Curvature_Gauge_G2H2[2][0][0][6] = 0;
  Curvature_Gauge_G2H2[2][0][0][7] = 0;
  Curvature_Gauge_G2H2[2][0][0][8] = 0;
  Curvature_Gauge_G2H2[2][0][1][0] = 0;
  Curvature_Gauge_G2H2[2][0][1][1] = 0;
  Curvature_Gauge_G2H2[2][0][1][2] = 0;
  Curvature_Gauge_G2H2[2][0][1][3] = 0;
  Curvature_Gauge_G2H2[2][0][1][4] = 0;
  Curvature_Gauge_G2H2[2][0][1][5] = 0;
  Curvature_Gauge_G2H2[2][0][1][6] = 0;
  Curvature_Gauge_G2H2[2][0][1][7] = 0;
  Curvature_Gauge_G2H2[2][0][1][8] = 0;
  Curvature_Gauge_G2H2[2][0][2][0] = 0;
  Curvature_Gauge_G2H2[2][0][2][1] = 0;
  Curvature_Gauge_G2H2[2][0][2][2] = 0;
  Curvature_Gauge_G2H2[2][0][2][3] = 0;
  Curvature_Gauge_G2H2[2][0][2][4] = 0;
  Curvature_Gauge_G2H2[2][0][2][5] = 0;
  Curvature_Gauge_G2H2[2][0][2][6] = 0;
  Curvature_Gauge_G2H2[2][0][2][7] = 0;
  Curvature_Gauge_G2H2[2][0][2][8] = 0;
  Curvature_Gauge_G2H2[2][0][3][0] = 0;
  Curvature_Gauge_G2H2[2][0][3][1] = 0;
  Curvature_Gauge_G2H2[2][0][3][2] = 0;
  Curvature_Gauge_G2H2[2][0][3][3] = 0;
  Curvature_Gauge_G2H2[2][0][3][4] = 0;
  Curvature_Gauge_G2H2[2][0][3][5] = 0;
  Curvature_Gauge_G2H2[2][0][3][6] = 0;
  Curvature_Gauge_G2H2[2][0][3][7] = 0;
  Curvature_Gauge_G2H2[2][0][3][8] = 0;
  Curvature_Gauge_G2H2[2][0][4][0] = 0;
  Curvature_Gauge_G2H2[2][0][4][1] = 0;
  Curvature_Gauge_G2H2[2][0][4][2] = 0;
  Curvature_Gauge_G2H2[2][0][4][3] = 0;
  Curvature_Gauge_G2H2[2][0][4][4] = 0;
  Curvature_Gauge_G2H2[2][0][4][5] = 0;
  Curvature_Gauge_G2H2[2][0][4][6] = 0;
  Curvature_Gauge_G2H2[2][0][4][7] = 0;
  Curvature_Gauge_G2H2[2][0][4][8] = 0;
  Curvature_Gauge_G2H2[2][0][5][0] = 0;
  Curvature_Gauge_G2H2[2][0][5][1] = 0;
  Curvature_Gauge_G2H2[2][0][5][2] = 0;
  Curvature_Gauge_G2H2[2][0][5][3] = 0;
  Curvature_Gauge_G2H2[2][0][5][4] = 0;
  Curvature_Gauge_G2H2[2][0][5][5] = 0;
  Curvature_Gauge_G2H2[2][0][5][6] = 0;
  Curvature_Gauge_G2H2[2][0][5][7] = 0;
  Curvature_Gauge_G2H2[2][0][5][8] = 0;
  Curvature_Gauge_G2H2[2][0][6][0] = 0;
  Curvature_Gauge_G2H2[2][0][6][1] = 0;
  Curvature_Gauge_G2H2[2][0][6][2] = 0;
  Curvature_Gauge_G2H2[2][0][6][3] = 0;
  Curvature_Gauge_G2H2[2][0][6][4] = 0;
  Curvature_Gauge_G2H2[2][0][6][5] = 0;
  Curvature_Gauge_G2H2[2][0][6][6] = 0;
  Curvature_Gauge_G2H2[2][0][6][7] = 0;
  Curvature_Gauge_G2H2[2][0][6][8] = 0;
  Curvature_Gauge_G2H2[2][0][7][0] = 0;
  Curvature_Gauge_G2H2[2][0][7][1] = 0;
  Curvature_Gauge_G2H2[2][0][7][2] = 0;
  Curvature_Gauge_G2H2[2][0][7][3] = 0;
  Curvature_Gauge_G2H2[2][0][7][4] = 0;
  Curvature_Gauge_G2H2[2][0][7][5] = 0;
  Curvature_Gauge_G2H2[2][0][7][6] = 0;
  Curvature_Gauge_G2H2[2][0][7][7] = 0;
  Curvature_Gauge_G2H2[2][0][7][8] = 0;
  Curvature_Gauge_G2H2[2][0][8][0] = 0;
  Curvature_Gauge_G2H2[2][0][8][1] = 0;
  Curvature_Gauge_G2H2[2][0][8][2] = 0;
  Curvature_Gauge_G2H2[2][0][8][3] = 0;
  Curvature_Gauge_G2H2[2][0][8][4] = 0;
  Curvature_Gauge_G2H2[2][0][8][5] = 0;
  Curvature_Gauge_G2H2[2][0][8][6] = 0;
  Curvature_Gauge_G2H2[2][0][8][7] = 0;
  Curvature_Gauge_G2H2[2][0][8][8] = 0;
  Curvature_Gauge_G2H2[2][1][0][0] = 0;
  Curvature_Gauge_G2H2[2][1][0][1] = 0;
  Curvature_Gauge_G2H2[2][1][0][2] = 0;
  Curvature_Gauge_G2H2[2][1][0][3] = 0;
  Curvature_Gauge_G2H2[2][1][0][4] = 0;
  Curvature_Gauge_G2H2[2][1][0][5] = 0;
  Curvature_Gauge_G2H2[2][1][0][6] = 0;
  Curvature_Gauge_G2H2[2][1][0][7] = 0;
  Curvature_Gauge_G2H2[2][1][0][8] = 0;
  Curvature_Gauge_G2H2[2][1][1][0] = 0;
  Curvature_Gauge_G2H2[2][1][1][1] = 0;
  Curvature_Gauge_G2H2[2][1][1][2] = 0;
  Curvature_Gauge_G2H2[2][1][1][3] = 0;
  Curvature_Gauge_G2H2[2][1][1][4] = 0;
  Curvature_Gauge_G2H2[2][1][1][5] = 0;
  Curvature_Gauge_G2H2[2][1][1][6] = 0;
  Curvature_Gauge_G2H2[2][1][1][7] = 0;
  Curvature_Gauge_G2H2[2][1][1][8] = 0;
  Curvature_Gauge_G2H2[2][1][2][0] = 0;
  Curvature_Gauge_G2H2[2][1][2][1] = 0;
  Curvature_Gauge_G2H2[2][1][2][2] = 0;
  Curvature_Gauge_G2H2[2][1][2][3] = 0;
  Curvature_Gauge_G2H2[2][1][2][4] = 0;
  Curvature_Gauge_G2H2[2][1][2][5] = 0;
  Curvature_Gauge_G2H2[2][1][2][6] = 0;
  Curvature_Gauge_G2H2[2][1][2][7] = 0;
  Curvature_Gauge_G2H2[2][1][2][8] = 0;
  Curvature_Gauge_G2H2[2][1][3][0] = 0;
  Curvature_Gauge_G2H2[2][1][3][1] = 0;
  Curvature_Gauge_G2H2[2][1][3][2] = 0;
  Curvature_Gauge_G2H2[2][1][3][3] = 0;
  Curvature_Gauge_G2H2[2][1][3][4] = 0;
  Curvature_Gauge_G2H2[2][1][3][5] = 0;
  Curvature_Gauge_G2H2[2][1][3][6] = 0;
  Curvature_Gauge_G2H2[2][1][3][7] = 0;
  Curvature_Gauge_G2H2[2][1][3][8] = 0;
  Curvature_Gauge_G2H2[2][1][4][0] = 0;
  Curvature_Gauge_G2H2[2][1][4][1] = 0;
  Curvature_Gauge_G2H2[2][1][4][2] = 0;
  Curvature_Gauge_G2H2[2][1][4][3] = 0;
  Curvature_Gauge_G2H2[2][1][4][4] = 0;
  Curvature_Gauge_G2H2[2][1][4][5] = 0;
  Curvature_Gauge_G2H2[2][1][4][6] = 0;
  Curvature_Gauge_G2H2[2][1][4][7] = 0;
  Curvature_Gauge_G2H2[2][1][4][8] = 0;
  Curvature_Gauge_G2H2[2][1][5][0] = 0;
  Curvature_Gauge_G2H2[2][1][5][1] = 0;
  Curvature_Gauge_G2H2[2][1][5][2] = 0;
  Curvature_Gauge_G2H2[2][1][5][3] = 0;
  Curvature_Gauge_G2H2[2][1][5][4] = 0;
  Curvature_Gauge_G2H2[2][1][5][5] = 0;
  Curvature_Gauge_G2H2[2][1][5][6] = 0;
  Curvature_Gauge_G2H2[2][1][5][7] = 0;
  Curvature_Gauge_G2H2[2][1][5][8] = 0;
  Curvature_Gauge_G2H2[2][1][6][0] = 0;
  Curvature_Gauge_G2H2[2][1][6][1] = 0;
  Curvature_Gauge_G2H2[2][1][6][2] = 0;
  Curvature_Gauge_G2H2[2][1][6][3] = 0;
  Curvature_Gauge_G2H2[2][1][6][4] = 0;
  Curvature_Gauge_G2H2[2][1][6][5] = 0;
  Curvature_Gauge_G2H2[2][1][6][6] = 0;
  Curvature_Gauge_G2H2[2][1][6][7] = 0;
  Curvature_Gauge_G2H2[2][1][6][8] = 0;
  Curvature_Gauge_G2H2[2][1][7][0] = 0;
  Curvature_Gauge_G2H2[2][1][7][1] = 0;
  Curvature_Gauge_G2H2[2][1][7][2] = 0;
  Curvature_Gauge_G2H2[2][1][7][3] = 0;
  Curvature_Gauge_G2H2[2][1][7][4] = 0;
  Curvature_Gauge_G2H2[2][1][7][5] = 0;
  Curvature_Gauge_G2H2[2][1][7][6] = 0;
  Curvature_Gauge_G2H2[2][1][7][7] = 0;
  Curvature_Gauge_G2H2[2][1][7][8] = 0;
  Curvature_Gauge_G2H2[2][1][8][0] = 0;
  Curvature_Gauge_G2H2[2][1][8][1] = 0;
  Curvature_Gauge_G2H2[2][1][8][2] = 0;
  Curvature_Gauge_G2H2[2][1][8][3] = 0;
  Curvature_Gauge_G2H2[2][1][8][4] = 0;
  Curvature_Gauge_G2H2[2][1][8][5] = 0;
  Curvature_Gauge_G2H2[2][1][8][6] = 0;
  Curvature_Gauge_G2H2[2][1][8][7] = 0;
  Curvature_Gauge_G2H2[2][1][8][8] = 0;
  Curvature_Gauge_G2H2[2][2][0][0] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][0][1] = 0;
  Curvature_Gauge_G2H2[2][2][0][2] = 0;
  Curvature_Gauge_G2H2[2][2][0][3] = 0;
  Curvature_Gauge_G2H2[2][2][0][4] = 0;
  Curvature_Gauge_G2H2[2][2][0][5] = 0;
  Curvature_Gauge_G2H2[2][2][0][6] = 0;
  Curvature_Gauge_G2H2[2][2][0][7] = 0;
  Curvature_Gauge_G2H2[2][2][0][8] = 0;
  Curvature_Gauge_G2H2[2][2][1][0] = 0;
  Curvature_Gauge_G2H2[2][2][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][1][2] = 0;
  Curvature_Gauge_G2H2[2][2][1][3] = 0;
  Curvature_Gauge_G2H2[2][2][1][4] = 0;
  Curvature_Gauge_G2H2[2][2][1][5] = 0;
  Curvature_Gauge_G2H2[2][2][1][6] = 0;
  Curvature_Gauge_G2H2[2][2][1][7] = 0;
  Curvature_Gauge_G2H2[2][2][1][8] = 0;
  Curvature_Gauge_G2H2[2][2][2][0] = 0;
  Curvature_Gauge_G2H2[2][2][2][1] = 0;
  Curvature_Gauge_G2H2[2][2][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][2][3] = 0;
  Curvature_Gauge_G2H2[2][2][2][4] = 0;
  Curvature_Gauge_G2H2[2][2][2][5] = 0;
  Curvature_Gauge_G2H2[2][2][2][6] = 0;
  Curvature_Gauge_G2H2[2][2][2][7] = 0;
  Curvature_Gauge_G2H2[2][2][2][8] = 0;
  Curvature_Gauge_G2H2[2][2][3][0] = 0;
  Curvature_Gauge_G2H2[2][2][3][1] = 0;
  Curvature_Gauge_G2H2[2][2][3][2] = 0;
  Curvature_Gauge_G2H2[2][2][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][3][4] = 0;
  Curvature_Gauge_G2H2[2][2][3][5] = 0;
  Curvature_Gauge_G2H2[2][2][3][6] = 0;
  Curvature_Gauge_G2H2[2][2][3][7] = 0;
  Curvature_Gauge_G2H2[2][2][3][8] = 0;
  Curvature_Gauge_G2H2[2][2][4][0] = 0;
  Curvature_Gauge_G2H2[2][2][4][1] = 0;
  Curvature_Gauge_G2H2[2][2][4][2] = 0;
  Curvature_Gauge_G2H2[2][2][4][3] = 0;
  Curvature_Gauge_G2H2[2][2][4][4] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][4][5] = 0;
  Curvature_Gauge_G2H2[2][2][4][6] = 0;
  Curvature_Gauge_G2H2[2][2][4][7] = 0;
  Curvature_Gauge_G2H2[2][2][4][8] = 0;
  Curvature_Gauge_G2H2[2][2][5][0] = 0;
  Curvature_Gauge_G2H2[2][2][5][1] = 0;
  Curvature_Gauge_G2H2[2][2][5][2] = 0;
  Curvature_Gauge_G2H2[2][2][5][3] = 0;
  Curvature_Gauge_G2H2[2][2][5][4] = 0;
  Curvature_Gauge_G2H2[2][2][5][5] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][5][6] = 0;
  Curvature_Gauge_G2H2[2][2][5][7] = 0;
  Curvature_Gauge_G2H2[2][2][5][8] = 0;
  Curvature_Gauge_G2H2[2][2][6][0] = 0;
  Curvature_Gauge_G2H2[2][2][6][1] = 0;
  Curvature_Gauge_G2H2[2][2][6][2] = 0;
  Curvature_Gauge_G2H2[2][2][6][3] = 0;
  Curvature_Gauge_G2H2[2][2][6][4] = 0;
  Curvature_Gauge_G2H2[2][2][6][5] = 0;
  Curvature_Gauge_G2H2[2][2][6][6] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][6][7] = 0;
  Curvature_Gauge_G2H2[2][2][6][8] = 0;
  Curvature_Gauge_G2H2[2][2][7][0] = 0;
  Curvature_Gauge_G2H2[2][2][7][1] = 0;
  Curvature_Gauge_G2H2[2][2][7][2] = 0;
  Curvature_Gauge_G2H2[2][2][7][3] = 0;
  Curvature_Gauge_G2H2[2][2][7][4] = 0;
  Curvature_Gauge_G2H2[2][2][7][5] = 0;
  Curvature_Gauge_G2H2[2][2][7][6] = 0;
  Curvature_Gauge_G2H2[2][2][7][7] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][7][8] = 0;
  Curvature_Gauge_G2H2[2][2][8][0] = 0;
  Curvature_Gauge_G2H2[2][2][8][1] = 0;
  Curvature_Gauge_G2H2[2][2][8][2] = 0;
  Curvature_Gauge_G2H2[2][2][8][3] = 0;
  Curvature_Gauge_G2H2[2][2][8][4] = 0;
  Curvature_Gauge_G2H2[2][2][8][5] = 0;
  Curvature_Gauge_G2H2[2][2][8][6] = 0;
  Curvature_Gauge_G2H2[2][2][8][7] = 0;
  Curvature_Gauge_G2H2[2][2][8][8] = 0;
  Curvature_Gauge_G2H2[2][3][0][0] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][0][1] = 0;
  Curvature_Gauge_G2H2[2][3][0][2] = 0;
  Curvature_Gauge_G2H2[2][3][0][3] = 0;
  Curvature_Gauge_G2H2[2][3][0][4] = 0;
  Curvature_Gauge_G2H2[2][3][0][5] = 0;
  Curvature_Gauge_G2H2[2][3][0][6] = 0;
  Curvature_Gauge_G2H2[2][3][0][7] = 0;
  Curvature_Gauge_G2H2[2][3][0][8] = 0;
  Curvature_Gauge_G2H2[2][3][1][0] = 0;
  Curvature_Gauge_G2H2[2][3][1][1] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][1][2] = 0;
  Curvature_Gauge_G2H2[2][3][1][3] = 0;
  Curvature_Gauge_G2H2[2][3][1][4] = 0;
  Curvature_Gauge_G2H2[2][3][1][5] = 0;
  Curvature_Gauge_G2H2[2][3][1][6] = 0;
  Curvature_Gauge_G2H2[2][3][1][7] = 0;
  Curvature_Gauge_G2H2[2][3][1][8] = 0;
  Curvature_Gauge_G2H2[2][3][2][0] = 0;
  Curvature_Gauge_G2H2[2][3][2][1] = 0;
  Curvature_Gauge_G2H2[2][3][2][2] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][2][3] = 0;
  Curvature_Gauge_G2H2[2][3][2][4] = 0;
  Curvature_Gauge_G2H2[2][3][2][5] = 0;
  Curvature_Gauge_G2H2[2][3][2][6] = 0;
  Curvature_Gauge_G2H2[2][3][2][7] = 0;
  Curvature_Gauge_G2H2[2][3][2][8] = 0;
  Curvature_Gauge_G2H2[2][3][3][0] = 0;
  Curvature_Gauge_G2H2[2][3][3][1] = 0;
  Curvature_Gauge_G2H2[2][3][3][2] = 0;
  Curvature_Gauge_G2H2[2][3][3][3] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][3][4] = 0;
  Curvature_Gauge_G2H2[2][3][3][5] = 0;
  Curvature_Gauge_G2H2[2][3][3][6] = 0;
  Curvature_Gauge_G2H2[2][3][3][7] = 0;
  Curvature_Gauge_G2H2[2][3][3][8] = 0;
  Curvature_Gauge_G2H2[2][3][4][0] = 0;
  Curvature_Gauge_G2H2[2][3][4][1] = 0;
  Curvature_Gauge_G2H2[2][3][4][2] = 0;
  Curvature_Gauge_G2H2[2][3][4][3] = 0;
  Curvature_Gauge_G2H2[2][3][4][4] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][4][5] = 0;
  Curvature_Gauge_G2H2[2][3][4][6] = 0;
  Curvature_Gauge_G2H2[2][3][4][7] = 0;
  Curvature_Gauge_G2H2[2][3][4][8] = 0;
  Curvature_Gauge_G2H2[2][3][5][0] = 0;
  Curvature_Gauge_G2H2[2][3][5][1] = 0;
  Curvature_Gauge_G2H2[2][3][5][2] = 0;
  Curvature_Gauge_G2H2[2][3][5][3] = 0;
  Curvature_Gauge_G2H2[2][3][5][4] = 0;
  Curvature_Gauge_G2H2[2][3][5][5] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][5][6] = 0;
  Curvature_Gauge_G2H2[2][3][5][7] = 0;
  Curvature_Gauge_G2H2[2][3][5][8] = 0;
  Curvature_Gauge_G2H2[2][3][6][0] = 0;
  Curvature_Gauge_G2H2[2][3][6][1] = 0;
  Curvature_Gauge_G2H2[2][3][6][2] = 0;
  Curvature_Gauge_G2H2[2][3][6][3] = 0;
  Curvature_Gauge_G2H2[2][3][6][4] = 0;
  Curvature_Gauge_G2H2[2][3][6][5] = 0;
  Curvature_Gauge_G2H2[2][3][6][6] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][6][7] = 0;
  Curvature_Gauge_G2H2[2][3][6][8] = 0;
  Curvature_Gauge_G2H2[2][3][7][0] = 0;
  Curvature_Gauge_G2H2[2][3][7][1] = 0;
  Curvature_Gauge_G2H2[2][3][7][2] = 0;
  Curvature_Gauge_G2H2[2][3][7][3] = 0;
  Curvature_Gauge_G2H2[2][3][7][4] = 0;
  Curvature_Gauge_G2H2[2][3][7][5] = 0;
  Curvature_Gauge_G2H2[2][3][7][6] = 0;
  Curvature_Gauge_G2H2[2][3][7][7] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[2][3][7][8] = 0;
  Curvature_Gauge_G2H2[2][3][8][0] = 0;
  Curvature_Gauge_G2H2[2][3][8][1] = 0;
  Curvature_Gauge_G2H2[2][3][8][2] = 0;
  Curvature_Gauge_G2H2[2][3][8][3] = 0;
  Curvature_Gauge_G2H2[2][3][8][4] = 0;
  Curvature_Gauge_G2H2[2][3][8][5] = 0;
  Curvature_Gauge_G2H2[2][3][8][6] = 0;
  Curvature_Gauge_G2H2[2][3][8][7] = 0;
  Curvature_Gauge_G2H2[2][3][8][8] = 0;
  Curvature_Gauge_G2H2[3][0][0][0] = 0;
  Curvature_Gauge_G2H2[3][0][0][1] = 0;
  Curvature_Gauge_G2H2[3][0][0][2] = 0;
  Curvature_Gauge_G2H2[3][0][0][3] = 0;
  Curvature_Gauge_G2H2[3][0][0][4] = 0;
  Curvature_Gauge_G2H2[3][0][0][5] = 0;
  Curvature_Gauge_G2H2[3][0][0][6] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][0][7] = 0;
  Curvature_Gauge_G2H2[3][0][0][8] = 0;
  Curvature_Gauge_G2H2[3][0][1][0] = 0;
  Curvature_Gauge_G2H2[3][0][1][1] = 0;
  Curvature_Gauge_G2H2[3][0][1][2] = 0;
  Curvature_Gauge_G2H2[3][0][1][3] = 0;
  Curvature_Gauge_G2H2[3][0][1][4] = 0;
  Curvature_Gauge_G2H2[3][0][1][5] = 0;
  Curvature_Gauge_G2H2[3][0][1][6] = 0;
  Curvature_Gauge_G2H2[3][0][1][7] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][1][8] = 0;
  Curvature_Gauge_G2H2[3][0][2][0] = 0;
  Curvature_Gauge_G2H2[3][0][2][1] = 0;
  Curvature_Gauge_G2H2[3][0][2][2] = 0;
  Curvature_Gauge_G2H2[3][0][2][3] = 0;
  Curvature_Gauge_G2H2[3][0][2][4] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][2][5] = 0;
  Curvature_Gauge_G2H2[3][0][2][6] = 0;
  Curvature_Gauge_G2H2[3][0][2][7] = 0;
  Curvature_Gauge_G2H2[3][0][2][8] = 0;
  Curvature_Gauge_G2H2[3][0][3][0] = 0;
  Curvature_Gauge_G2H2[3][0][3][1] = 0;
  Curvature_Gauge_G2H2[3][0][3][2] = 0;
  Curvature_Gauge_G2H2[3][0][3][3] = 0;
  Curvature_Gauge_G2H2[3][0][3][4] = 0;
  Curvature_Gauge_G2H2[3][0][3][5] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][3][6] = 0;
  Curvature_Gauge_G2H2[3][0][3][7] = 0;
  Curvature_Gauge_G2H2[3][0][3][8] = 0;
  Curvature_Gauge_G2H2[3][0][4][0] = 0;
  Curvature_Gauge_G2H2[3][0][4][1] = 0;
  Curvature_Gauge_G2H2[3][0][4][2] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][4][3] = 0;
  Curvature_Gauge_G2H2[3][0][4][4] = 0;
  Curvature_Gauge_G2H2[3][0][4][5] = 0;
  Curvature_Gauge_G2H2[3][0][4][6] = 0;
  Curvature_Gauge_G2H2[3][0][4][7] = 0;
  Curvature_Gauge_G2H2[3][0][4][8] = 0;
  Curvature_Gauge_G2H2[3][0][5][0] = 0;
  Curvature_Gauge_G2H2[3][0][5][1] = 0;
  Curvature_Gauge_G2H2[3][0][5][2] = 0;
  Curvature_Gauge_G2H2[3][0][5][3] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][5][4] = 0;
  Curvature_Gauge_G2H2[3][0][5][5] = 0;
  Curvature_Gauge_G2H2[3][0][5][6] = 0;
  Curvature_Gauge_G2H2[3][0][5][7] = 0;
  Curvature_Gauge_G2H2[3][0][5][8] = 0;
  Curvature_Gauge_G2H2[3][0][6][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][6][1] = 0;
  Curvature_Gauge_G2H2[3][0][6][2] = 0;
  Curvature_Gauge_G2H2[3][0][6][3] = 0;
  Curvature_Gauge_G2H2[3][0][6][4] = 0;
  Curvature_Gauge_G2H2[3][0][6][5] = 0;
  Curvature_Gauge_G2H2[3][0][6][6] = 0;
  Curvature_Gauge_G2H2[3][0][6][7] = 0;
  Curvature_Gauge_G2H2[3][0][6][8] = 0;
  Curvature_Gauge_G2H2[3][0][7][0] = 0;
  Curvature_Gauge_G2H2[3][0][7][1] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][7][2] = 0;
  Curvature_Gauge_G2H2[3][0][7][3] = 0;
  Curvature_Gauge_G2H2[3][0][7][4] = 0;
  Curvature_Gauge_G2H2[3][0][7][5] = 0;
  Curvature_Gauge_G2H2[3][0][7][6] = 0;
  Curvature_Gauge_G2H2[3][0][7][7] = 0;
  Curvature_Gauge_G2H2[3][0][7][8] = 0;
  Curvature_Gauge_G2H2[3][0][8][0] = 0;
  Curvature_Gauge_G2H2[3][0][8][1] = 0;
  Curvature_Gauge_G2H2[3][0][8][2] = 0;
  Curvature_Gauge_G2H2[3][0][8][3] = 0;
  Curvature_Gauge_G2H2[3][0][8][4] = 0;
  Curvature_Gauge_G2H2[3][0][8][5] = 0;
  Curvature_Gauge_G2H2[3][0][8][6] = 0;
  Curvature_Gauge_G2H2[3][0][8][7] = 0;
  Curvature_Gauge_G2H2[3][0][8][8] = 0;
  Curvature_Gauge_G2H2[3][1][0][0] = 0;
  Curvature_Gauge_G2H2[3][1][0][1] = 0;
  Curvature_Gauge_G2H2[3][1][0][2] = 0;
  Curvature_Gauge_G2H2[3][1][0][3] = 0;
  Curvature_Gauge_G2H2[3][1][0][4] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][0][5] = 0;
  Curvature_Gauge_G2H2[3][1][0][6] = 0;
  Curvature_Gauge_G2H2[3][1][0][7] = 0;
  Curvature_Gauge_G2H2[3][1][0][8] = 0;
  Curvature_Gauge_G2H2[3][1][1][0] = 0;
  Curvature_Gauge_G2H2[3][1][1][1] = 0;
  Curvature_Gauge_G2H2[3][1][1][2] = 0;
  Curvature_Gauge_G2H2[3][1][1][3] = 0;
  Curvature_Gauge_G2H2[3][1][1][4] = 0;
  Curvature_Gauge_G2H2[3][1][1][5] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][1][6] = 0;
  Curvature_Gauge_G2H2[3][1][1][7] = 0;
  Curvature_Gauge_G2H2[3][1][1][8] = 0;
  Curvature_Gauge_G2H2[3][1][2][0] = 0;
  Curvature_Gauge_G2H2[3][1][2][1] = 0;
  Curvature_Gauge_G2H2[3][1][2][2] = 0;
  Curvature_Gauge_G2H2[3][1][2][3] = 0;
  Curvature_Gauge_G2H2[3][1][2][4] = 0;
  Curvature_Gauge_G2H2[3][1][2][5] = 0;
  Curvature_Gauge_G2H2[3][1][2][6] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][2][7] = 0;
  Curvature_Gauge_G2H2[3][1][2][8] = 0;
  Curvature_Gauge_G2H2[3][1][3][0] = 0;
  Curvature_Gauge_G2H2[3][1][3][1] = 0;
  Curvature_Gauge_G2H2[3][1][3][2] = 0;
  Curvature_Gauge_G2H2[3][1][3][3] = 0;
  Curvature_Gauge_G2H2[3][1][3][4] = 0;
  Curvature_Gauge_G2H2[3][1][3][5] = 0;
  Curvature_Gauge_G2H2[3][1][3][6] = 0;
  Curvature_Gauge_G2H2[3][1][3][7] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][3][8] = 0;
  Curvature_Gauge_G2H2[3][1][4][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][4][1] = 0;
  Curvature_Gauge_G2H2[3][1][4][2] = 0;
  Curvature_Gauge_G2H2[3][1][4][3] = 0;
  Curvature_Gauge_G2H2[3][1][4][4] = 0;
  Curvature_Gauge_G2H2[3][1][4][5] = 0;
  Curvature_Gauge_G2H2[3][1][4][6] = 0;
  Curvature_Gauge_G2H2[3][1][4][7] = 0;
  Curvature_Gauge_G2H2[3][1][4][8] = 0;
  Curvature_Gauge_G2H2[3][1][5][0] = 0;
  Curvature_Gauge_G2H2[3][1][5][1] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][5][2] = 0;
  Curvature_Gauge_G2H2[3][1][5][3] = 0;
  Curvature_Gauge_G2H2[3][1][5][4] = 0;
  Curvature_Gauge_G2H2[3][1][5][5] = 0;
  Curvature_Gauge_G2H2[3][1][5][6] = 0;
  Curvature_Gauge_G2H2[3][1][5][7] = 0;
  Curvature_Gauge_G2H2[3][1][5][8] = 0;
  Curvature_Gauge_G2H2[3][1][6][0] = 0;
  Curvature_Gauge_G2H2[3][1][6][1] = 0;
  Curvature_Gauge_G2H2[3][1][6][2] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][6][3] = 0;
  Curvature_Gauge_G2H2[3][1][6][4] = 0;
  Curvature_Gauge_G2H2[3][1][6][5] = 0;
  Curvature_Gauge_G2H2[3][1][6][6] = 0;
  Curvature_Gauge_G2H2[3][1][6][7] = 0;
  Curvature_Gauge_G2H2[3][1][6][8] = 0;
  Curvature_Gauge_G2H2[3][1][7][0] = 0;
  Curvature_Gauge_G2H2[3][1][7][1] = 0;
  Curvature_Gauge_G2H2[3][1][7][2] = 0;
  Curvature_Gauge_G2H2[3][1][7][3] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][7][4] = 0;
  Curvature_Gauge_G2H2[3][1][7][5] = 0;
  Curvature_Gauge_G2H2[3][1][7][6] = 0;
  Curvature_Gauge_G2H2[3][1][7][7] = 0;
  Curvature_Gauge_G2H2[3][1][7][8] = 0;
  Curvature_Gauge_G2H2[3][1][8][0] = 0;
  Curvature_Gauge_G2H2[3][1][8][1] = 0;
  Curvature_Gauge_G2H2[3][1][8][2] = 0;
  Curvature_Gauge_G2H2[3][1][8][3] = 0;
  Curvature_Gauge_G2H2[3][1][8][4] = 0;
  Curvature_Gauge_G2H2[3][1][8][5] = 0;
  Curvature_Gauge_G2H2[3][1][8][6] = 0;
  Curvature_Gauge_G2H2[3][1][8][7] = 0;
  Curvature_Gauge_G2H2[3][1][8][8] = 0;
  Curvature_Gauge_G2H2[3][2][0][0] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][0][1] = 0;
  Curvature_Gauge_G2H2[3][2][0][2] = 0;
  Curvature_Gauge_G2H2[3][2][0][3] = 0;
  Curvature_Gauge_G2H2[3][2][0][4] = 0;
  Curvature_Gauge_G2H2[3][2][0][5] = 0;
  Curvature_Gauge_G2H2[3][2][0][6] = 0;
  Curvature_Gauge_G2H2[3][2][0][7] = 0;
  Curvature_Gauge_G2H2[3][2][0][8] = 0;
  Curvature_Gauge_G2H2[3][2][1][0] = 0;
  Curvature_Gauge_G2H2[3][2][1][1] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][1][2] = 0;
  Curvature_Gauge_G2H2[3][2][1][3] = 0;
  Curvature_Gauge_G2H2[3][2][1][4] = 0;
  Curvature_Gauge_G2H2[3][2][1][5] = 0;
  Curvature_Gauge_G2H2[3][2][1][6] = 0;
  Curvature_Gauge_G2H2[3][2][1][7] = 0;
  Curvature_Gauge_G2H2[3][2][1][8] = 0;
  Curvature_Gauge_G2H2[3][2][2][0] = 0;
  Curvature_Gauge_G2H2[3][2][2][1] = 0;
  Curvature_Gauge_G2H2[3][2][2][2] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][2][3] = 0;
  Curvature_Gauge_G2H2[3][2][2][4] = 0;
  Curvature_Gauge_G2H2[3][2][2][5] = 0;
  Curvature_Gauge_G2H2[3][2][2][6] = 0;
  Curvature_Gauge_G2H2[3][2][2][7] = 0;
  Curvature_Gauge_G2H2[3][2][2][8] = 0;
  Curvature_Gauge_G2H2[3][2][3][0] = 0;
  Curvature_Gauge_G2H2[3][2][3][1] = 0;
  Curvature_Gauge_G2H2[3][2][3][2] = 0;
  Curvature_Gauge_G2H2[3][2][3][3] = SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][3][4] = 0;
  Curvature_Gauge_G2H2[3][2][3][5] = 0;
  Curvature_Gauge_G2H2[3][2][3][6] = 0;
  Curvature_Gauge_G2H2[3][2][3][7] = 0;
  Curvature_Gauge_G2H2[3][2][3][8] = 0;
  Curvature_Gauge_G2H2[3][2][4][0] = 0;
  Curvature_Gauge_G2H2[3][2][4][1] = 0;
  Curvature_Gauge_G2H2[3][2][4][2] = 0;
  Curvature_Gauge_G2H2[3][2][4][3] = 0;
  Curvature_Gauge_G2H2[3][2][4][4] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][4][5] = 0;
  Curvature_Gauge_G2H2[3][2][4][6] = 0;
  Curvature_Gauge_G2H2[3][2][4][7] = 0;
  Curvature_Gauge_G2H2[3][2][4][8] = 0;
  Curvature_Gauge_G2H2[3][2][5][0] = 0;
  Curvature_Gauge_G2H2[3][2][5][1] = 0;
  Curvature_Gauge_G2H2[3][2][5][2] = 0;
  Curvature_Gauge_G2H2[3][2][5][3] = 0;
  Curvature_Gauge_G2H2[3][2][5][4] = 0;
  Curvature_Gauge_G2H2[3][2][5][5] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][5][6] = 0;
  Curvature_Gauge_G2H2[3][2][5][7] = 0;
  Curvature_Gauge_G2H2[3][2][5][8] = 0;
  Curvature_Gauge_G2H2[3][2][6][0] = 0;
  Curvature_Gauge_G2H2[3][2][6][1] = 0;
  Curvature_Gauge_G2H2[3][2][6][2] = 0;
  Curvature_Gauge_G2H2[3][2][6][3] = 0;
  Curvature_Gauge_G2H2[3][2][6][4] = 0;
  Curvature_Gauge_G2H2[3][2][6][5] = 0;
  Curvature_Gauge_G2H2[3][2][6][6] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][6][7] = 0;
  Curvature_Gauge_G2H2[3][2][6][8] = 0;
  Curvature_Gauge_G2H2[3][2][7][0] = 0;
  Curvature_Gauge_G2H2[3][2][7][1] = 0;
  Curvature_Gauge_G2H2[3][2][7][2] = 0;
  Curvature_Gauge_G2H2[3][2][7][3] = 0;
  Curvature_Gauge_G2H2[3][2][7][4] = 0;
  Curvature_Gauge_G2H2[3][2][7][5] = 0;
  Curvature_Gauge_G2H2[3][2][7][6] = 0;
  Curvature_Gauge_G2H2[3][2][7][7] =
      -SMConstants.C_g * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][2][7][8] = 0;
  Curvature_Gauge_G2H2[3][2][8][0] = 0;
  Curvature_Gauge_G2H2[3][2][8][1] = 0;
  Curvature_Gauge_G2H2[3][2][8][2] = 0;
  Curvature_Gauge_G2H2[3][2][8][3] = 0;
  Curvature_Gauge_G2H2[3][2][8][4] = 0;
  Curvature_Gauge_G2H2[3][2][8][5] = 0;
  Curvature_Gauge_G2H2[3][2][8][6] = 0;
  Curvature_Gauge_G2H2[3][2][8][7] = 0;
  Curvature_Gauge_G2H2[3][2][8][8] = 0;
  Curvature_Gauge_G2H2[3][3][0][0] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][0][1] = 0;
  Curvature_Gauge_G2H2[3][3][0][2] = 0;
  Curvature_Gauge_G2H2[3][3][0][3] = 0;
  Curvature_Gauge_G2H2[3][3][0][4] = 0;
  Curvature_Gauge_G2H2[3][3][0][5] = 0;
  Curvature_Gauge_G2H2[3][3][0][6] = 0;
  Curvature_Gauge_G2H2[3][3][0][7] = 0;
  Curvature_Gauge_G2H2[3][3][0][8] = 0;
  Curvature_Gauge_G2H2[3][3][1][0] = 0;
  Curvature_Gauge_G2H2[3][3][1][1] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][1][2] = 0;
  Curvature_Gauge_G2H2[3][3][1][3] = 0;
  Curvature_Gauge_G2H2[3][3][1][4] = 0;
  Curvature_Gauge_G2H2[3][3][1][5] = 0;
  Curvature_Gauge_G2H2[3][3][1][6] = 0;
  Curvature_Gauge_G2H2[3][3][1][7] = 0;
  Curvature_Gauge_G2H2[3][3][1][8] = 0;
  Curvature_Gauge_G2H2[3][3][2][0] = 0;
  Curvature_Gauge_G2H2[3][3][2][1] = 0;
  Curvature_Gauge_G2H2[3][3][2][2] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][2][3] = 0;
  Curvature_Gauge_G2H2[3][3][2][4] = 0;
  Curvature_Gauge_G2H2[3][3][2][5] = 0;
  Curvature_Gauge_G2H2[3][3][2][6] = 0;
  Curvature_Gauge_G2H2[3][3][2][7] = 0;
  Curvature_Gauge_G2H2[3][3][2][8] = 0;
  Curvature_Gauge_G2H2[3][3][3][0] = 0;
  Curvature_Gauge_G2H2[3][3][3][1] = 0;
  Curvature_Gauge_G2H2[3][3][3][2] = 0;
  Curvature_Gauge_G2H2[3][3][3][3] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][3][4] = 0;
  Curvature_Gauge_G2H2[3][3][3][5] = 0;
  Curvature_Gauge_G2H2[3][3][3][6] = 0;
  Curvature_Gauge_G2H2[3][3][3][7] = 0;
  Curvature_Gauge_G2H2[3][3][3][8] = 0;
  Curvature_Gauge_G2H2[3][3][4][0] = 0;
  Curvature_Gauge_G2H2[3][3][4][1] = 0;
  Curvature_Gauge_G2H2[3][3][4][2] = 0;
  Curvature_Gauge_G2H2[3][3][4][3] = 0;
  Curvature_Gauge_G2H2[3][3][4][4] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][4][5] = 0;
  Curvature_Gauge_G2H2[3][3][4][6] = 0;
  Curvature_Gauge_G2H2[3][3][4][7] = 0;
  Curvature_Gauge_G2H2[3][3][4][8] = 0;
  Curvature_Gauge_G2H2[3][3][5][0] = 0;
  Curvature_Gauge_G2H2[3][3][5][1] = 0;
  Curvature_Gauge_G2H2[3][3][5][2] = 0;
  Curvature_Gauge_G2H2[3][3][5][3] = 0;
  Curvature_Gauge_G2H2[3][3][5][4] = 0;
  Curvature_Gauge_G2H2[3][3][5][5] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][5][6] = 0;
  Curvature_Gauge_G2H2[3][3][5][7] = 0;
  Curvature_Gauge_G2H2[3][3][5][8] = 0;
  Curvature_Gauge_G2H2[3][3][6][0] = 0;
  Curvature_Gauge_G2H2[3][3][6][1] = 0;
  Curvature_Gauge_G2H2[3][3][6][2] = 0;
  Curvature_Gauge_G2H2[3][3][6][3] = 0;
  Curvature_Gauge_G2H2[3][3][6][4] = 0;
  Curvature_Gauge_G2H2[3][3][6][5] = 0;
  Curvature_Gauge_G2H2[3][3][6][6] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][6][7] = 0;
  Curvature_Gauge_G2H2[3][3][6][8] = 0;
  Curvature_Gauge_G2H2[3][3][7][0] = 0;
  Curvature_Gauge_G2H2[3][3][7][1] = 0;
  Curvature_Gauge_G2H2[3][3][7][2] = 0;
  Curvature_Gauge_G2H2[3][3][7][3] = 0;
  Curvature_Gauge_G2H2[3][3][7][4] = 0;
  Curvature_Gauge_G2H2[3][3][7][5] = 0;
  Curvature_Gauge_G2H2[3][3][7][6] = 0;
  Curvature_Gauge_G2H2[3][3][7][7] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][7][8] = 0;
  Curvature_Gauge_G2H2[3][3][8][0] = 0;
  Curvature_Gauge_G2H2[3][3][8][1] = 0;
  Curvature_Gauge_G2H2[3][3][8][2] = 0;
  Curvature_Gauge_G2H2[3][3][8][3] = 0;
  Curvature_Gauge_G2H2[3][3][8][4] = 0;
  Curvature_Gauge_G2H2[3][3][8][5] = 0;
  Curvature_Gauge_G2H2[3][3][8][6] = 0;
  Curvature_Gauge_G2H2[3][3][8][7] = 0;
  Curvature_Gauge_G2H2[3][3][8][8] = 0;

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
      Curvature_Quark_F2H1[i][j][1] = YIJR2(i, j);
      Curvature_Quark_F2H1[i][j][2] = 0;
      Curvature_Quark_F2H1[i][j][3] = YIJE2(i, j);
      Curvature_Quark_F2H1[i][j][4] = 0;
      Curvature_Quark_F2H1[i][j][5] = YIJP2(i, j);
      Curvature_Quark_F2H1[i][j][6] = 0;
      Curvature_Quark_F2H1[i][j][7] = YIJS2(i, j);
      Curvature_Quark_F2H1[i][j][8] = 0;

      if (Type == 1 or Type == 3)
      {
        Curvature_Quark_F2H1[i][j][1] += YIJRD(i, j);
        Curvature_Quark_F2H1[i][j][3] += YIJED(i, j);
        Curvature_Quark_F2H1[i][j][7] += YIJSD(i, j);
        Curvature_Quark_F2H1[i][j][5] += YIJPD(i, j);
      }
      else
      {
        Curvature_Quark_F2H1[i][j][0] += YIJRD(i, j);
        Curvature_Quark_F2H1[i][j][2] += YIJED(i, j);
        Curvature_Quark_F2H1[i][j][6] += YIJSD(i, j);
        Curvature_Quark_F2H1[i][j][4] += YIJPD(i, j);
      }
    }
  }

  for (std::size_t i = 0; i < NLepton; i++)
  {
    for (std::size_t j = 0; j < NLepton; j++)
    {
      Curvature_Lepton_F2H1[i][j][8] = 0;
      if (Type == 1 or Type == 4)
      {
        Curvature_Lepton_F2H1[i][j][0] = 0;
        Curvature_Lepton_F2H1[i][j][1] = YIJRL(i, j);
        Curvature_Lepton_F2H1[i][j][2] = 0;
        Curvature_Lepton_F2H1[i][j][3] = YIJEL(i, j);
        Curvature_Lepton_F2H1[i][j][4] = 0;
        Curvature_Lepton_F2H1[i][j][5] = YIJPL(i, j);
        Curvature_Lepton_F2H1[i][j][6] = 0;
        Curvature_Lepton_F2H1[i][j][7] = YIJSL(i, j);
      }
      else
      {
        Curvature_Lepton_F2H1[i][j][0] = YIJRL(i, j);
        Curvature_Lepton_F2H1[i][j][1] = 0;
        Curvature_Lepton_F2H1[i][j][2] = YIJEL(i, j);
        Curvature_Lepton_F2H1[i][j][3] = 0;
        Curvature_Lepton_F2H1[i][j][4] = YIJPL(i, j);
        Curvature_Lepton_F2H1[i][j][5] = 0;
        Curvature_Lepton_F2H1[i][j][6] = YIJSL(i, j);
        Curvature_Lepton_F2H1[i][j][7] = 0;
      }
    }
  }

  SetCurvatureDone = true;
}

bool Class_Potential_N2HDM::CalculateDebyeSimplified()
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

  CTempC1 += 1.0 / 48.0 * (2 * NL7);
  CTempC2 += 1.0 / 48.0 * (2 * NL8);
  CTempCS = 1.0 / 6.0 * (NL7 + NL8) + 1.0 / 8.0 * NL6;

  DebyeHiggs[0][0] = CTempC1;
  DebyeHiggs[2][2] = CTempC1;
  DebyeHiggs[1][1] = CTempC2;
  DebyeHiggs[3][3] = CTempC2;
  DebyeHiggs[4][4] = CTempC1;
  DebyeHiggs[6][6] = CTempC1;
  DebyeHiggs[5][5] = CTempC2;
  DebyeHiggs[7][7] = CTempC2;
  DebyeHiggs[8][8] = CTempCS;

  return true;
}

bool Class_Potential_N2HDM::CalculateDebyeGaugeSimplified()
{
  DebyeGauge[0][0] = 2 * SMConstants.C_g * SMConstants.C_g;
  DebyeGauge[1][1] = 2 * SMConstants.C_g * SMConstants.C_g;
  DebyeGauge[2][2] = 2 * SMConstants.C_g * SMConstants.C_g;
  DebyeGauge[3][3] = 2 * SMConstants.C_gs * SMConstants.C_gs;

  return true;
}

double
Class_Potential_N2HDM::VTreeSimplified(const std::vector<double> &v) const
{
  (void)v;
  double res = 0;

  return res;
}

double
Class_Potential_N2HDM::VCounterSimplified(const std::vector<double> &v) const
{
  (void)v;
  if (not UseVCounterSimplified) return 0;
  double res = 0;
  return res;
}

void Class_Potential_N2HDM::Debugging(const std::vector<double> &input,
                                      std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
