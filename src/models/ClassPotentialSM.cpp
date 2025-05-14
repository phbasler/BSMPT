// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/models/ClassPotentialSM.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
using namespace Eigen;

/**
 * @file
 * Implementation of the SM
 */

namespace BSMPT
{
namespace Models
{

Class_SM::Class_SM(const ISMConstants &smConstants)
    : Class_Potential_Origin(smConstants)
{
  Model         = ModelID::ModelIDs::SM;
  NNeutralHiggs = 2; // number of neutral Higgs bosons at T = 0
  NChargedHiggs = 2; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar   = 2; // number of parameters in the tree-Level Lagrangian
  nParCT = 6; // number of parameters in the counterterm potential

  nVEV = 1; // number of VEVs to minimize the potential

  NHiggs = NNeutralHiggs + NChargedHiggs;

  VevOrder.resize(nVEV);
  VevOrder[0] = 2;

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_SM::~Class_SM()
{
}

std::vector<std::string> Class_SM::addLegendCT() const
{
  std::vector<std::string> labels;
  labels.push_back("dmuSq");
  labels.push_back("dlambda");
  labels.push_back("dT1");
  labels.push_back("dT2");
  labels.push_back("dT3");
  labels.push_back("dT4");
  return labels;
}

std::vector<std::string> Class_SM::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c");
  labels.push_back("v_c");
  labels.push_back("omega_c/T_c");
  labels.push_back("omega_c");
  return labels;
}

std::vector<std::string> Class_SM::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;
  particles.resize(NHiggs);

  particles[0] = "G+";
  particles[1] = "G-";
  particles[2] = "G0";
  particles[3] = "H";

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

std::vector<std::string> Class_SM::addLegendVEV() const
{
  std::vector<std::string> labels;
  labels.push_back("omega");
  return labels;
}

void Class_SM::ReadAndSet(const std::string &linestr, std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 2; k++)
  {
    ss >> tmp;

    if (k == 1)
      par[0] = tmp; // muSq
    else if (k == 2)
      par[1] = tmp; // lambda
  }

  set_gen(par);
  return;
}

void Class_SM::set_gen(const std::vector<double> &par)
{
  v0 = SMConstants.C_vev0;

  // uncomment if you want to load muSq from data
  // muSq = par[0];

  // muSq calculated from SM Higgs Mass directly
  (void)par;
  muSq = -std::pow(SMConstants.C_MassSMHiggs, 0.2e1) / 2;

  lambda = -muSq * std::pow(v0, -0.2e1);

  scale = v0;

  vevTreeMin.resize(nVEV);
  vevTree.resize(NHiggs);

  vevTreeMin[0] = v0;

  vevTree = MinimizeOrderVEV(vevTreeMin);
  if (!SetCurvatureDone) SetCurvatureArrays();
}

void Class_SM::set_CT_Pot_Par(const std::vector<double> &par)
{
  dmuSq   = par[0];
  dlambda = par[1];
  dT1     = par[2];
  dT2     = par[3];
  dT3     = par[4];
  dT4     = par[5];

  Curvature_Higgs_CT_L1[0] = dT1;
  Curvature_Higgs_CT_L1[1] = dT2;
  Curvature_Higgs_CT_L1[2] = dT3;
  Curvature_Higgs_CT_L1[3] = dT4;

  Curvature_Higgs_CT_L2[0][0] = dmuSq;
  Curvature_Higgs_CT_L2[1][1] = dmuSq;
  Curvature_Higgs_CT_L2[2][2] = dmuSq;
  Curvature_Higgs_CT_L2[3][3] = dmuSq;

  Curvature_Higgs_CT_L4[0][0][0][0] = 6 * dlambda;
  Curvature_Higgs_CT_L4[0][0][1][1] = 2 * dlambda;
  Curvature_Higgs_CT_L4[0][0][2][2] = 2 * dlambda;
  Curvature_Higgs_CT_L4[0][0][3][3] = 2 * dlambda;
  Curvature_Higgs_CT_L4[1][1][1][1] = 6 * dlambda;
  Curvature_Higgs_CT_L4[1][1][2][2] = 2 * dlambda;
  Curvature_Higgs_CT_L4[1][1][3][3] = 2 * dlambda;
  Curvature_Higgs_CT_L4[2][2][2][2] = 6 * dlambda;
  Curvature_Higgs_CT_L4[2][2][3][3] = 2 * dlambda;
  Curvature_Higgs_CT_L4[3][3][3][3] = 6 * dlambda;

  sym4Dim(Curvature_Higgs_CT_L4, NHiggs, NHiggs, NHiggs, NHiggs);
}

void Class_SM::write() const
{
  std::stringstream ss;
  typedef std::numeric_limits<double> dbl;
  ss.precision(dbl::max_digits10);

  ss << "The parameters are : "
     << "\n"
     << "\tmuSq = " << muSq << "\n"
     << "\tlambda = " << lambda << "\n"
     << "\tv0 = " << v0 << "\n";

  ss << "The counterterm parameters are : "
     << "\n";
  ss << "\tdmuSq = " << dmuSq << "\n"
     << "\tdlambda = " << dlambda << "\n"
     << "\tdT1 = " << dT1 << "\n"
     << "\tdT2 = " << dT2 << "\n"
     << "\tdT3 = " << dT3 << "\n"
     << "\tdT4 = " << dT4 << "\n";

  ss << "The scale is given by mu = " << scale << " GeV "
     << "\n";

  std::vector<double> HiggsMasses;
  HiggsMasses = HiggsMassesSquared(vevTree, 0);

  ss << "The mass spectrum is given by :\n";
  ss << "m_{G^+}^2 = " << HiggsMasses[0] << " GeV^2 \n"
     << "m_{G^-}^2 = " << HiggsMasses[1] << " GeV^2 \n"
     << "m_{G^0}^2 = " << HiggsMasses[2] << " GeV^2 \n"
     << "m_{H_SM} = " << std::sqrt(HiggsMasses[3]) << " GeV \n";

  Logger::Write(LoggingLevel::Default, ss.str());
}

std::vector<double> Class_SM::calc_CT() const
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

  parCT.push_back(HesseWeinberg(2, 2) / 2 -
                  0.3e1 / 0.2e1 * HesseWeinberg(3, 3)); // dmuSq
  parCT.push_back((-HesseWeinberg(2, 2) + HesseWeinberg(3, 3)) * pow(v0, -2) /
                  2);                                           // dlambda
  parCT.push_back(-NablaWeinberg(0));                           // dT1
  parCT.push_back(-NablaWeinberg(1));                           // dT2
  parCT.push_back(HesseWeinberg(3, 3) * v0 - NablaWeinberg(2)); // dT3
  parCT.push_back(-NablaWeinberg(3));                           // dT4

  return parCT;
}

void Class_SM::AdjustRotationMatrix()
{
}

void Class_SM::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsDone) CalculatePhysicalCouplings();

  if (CalculatedTripleCopulings) return;
  CalculatedTripleCopulings = true;

  std::vector<std::size_t> HiggsOrder(NHiggs);

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsOrder[i] = i;
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

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
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

void Class_SM::SetCurvatureArrays()
{
  initVectors();

  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

  Curvature_Higgs_L2[0][0] = muSq;
  Curvature_Higgs_L2[1][1] = muSq;
  Curvature_Higgs_L2[2][2] = muSq;
  Curvature_Higgs_L2[3][3] = muSq;

  Curvature_Higgs_L4[0][0][0][0] = 6 * lambda;
  Curvature_Higgs_L4[0][0][1][1] = 2 * lambda;
  Curvature_Higgs_L4[0][0][2][2] = 2 * lambda;
  Curvature_Higgs_L4[0][0][3][3] = 2 * lambda;
  Curvature_Higgs_L4[1][1][1][1] = 6 * lambda;
  Curvature_Higgs_L4[1][1][2][2] = 2 * lambda;
  Curvature_Higgs_L4[1][1][3][3] = 2 * lambda;
  Curvature_Higgs_L4[2][2][2][2] = 6 * lambda;
  Curvature_Higgs_L4[2][2][3][3] = 2 * lambda;
  Curvature_Higgs_L4[3][3][3][3] = 6 * lambda;

  sym4Dim(Curvature_Higgs_L4, NHiggs, NHiggs, NHiggs, NHiggs);

  Curvature_Gauge_G2H2[0][0][0][0] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][0][2] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][1][3] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][2][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][3][1] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][0][0] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][0][3] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][1][2] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][2][1] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][3][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][0][0] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][0][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][1][1] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][2][2] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][3][3] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][0][2] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][1][3] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][2][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][3][1] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][0][3] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][1][2] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][2][1] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][3][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][0][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][1][1] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][2][2] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][3][3] =
      -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][3][0][0] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][1][1] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][2][2] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][3][3] =
      SMConstants.C_gs * SMConstants.C_gs / 0.2e1;

  MatrixXcd YIJQc1(NQuarks, NQuarks), YIJQc2(NQuarks, NQuarks),
      YIJQc2OI(NQuarks, NQuarks), YIJQg0(NQuarks, NQuarks),
      YIJQg0OI(NQuarks, NQuarks), YIJQh1(NQuarks, NQuarks),
      YIJQh2(NQuarks, NQuarks), YIJQh3(NQuarks, NQuarks);

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

  std::complex<double> VC11, VC12, VC13, VC21, VC22, VC23, VC31, VC32, VC33;
  VC11 = std::conj(SMConstants.C_Vud);
  VC12 = std::conj(SMConstants.C_Vus);
  VC13 = std::conj(SMConstants.C_Vub);
  VC21 = std::conj(SMConstants.C_Vcd);
  VC22 = std::conj(SMConstants.C_Vcs);
  VC23 = std::conj(SMConstants.C_Vcb);
  VC31 = std::conj(SMConstants.C_Vtd);
  VC32 = std::conj(SMConstants.C_Vts);
  VC33 = std::conj(SMConstants.C_Vtb);

  Curvature_Lepton_F2H1[0][1][2] = 0.1e1 / v0 * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[0][1][3] = II / v0 * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][0][2] = 0.1e1 / v0 * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][0][3] = II / v0 * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][6][0] = 0.1e1 / v0 * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][6][1] = II / v0 * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[2][3][2] = 0.1e1 / v0 * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[2][3][3] = II / v0 * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][2][2] = 0.1e1 / v0 * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][2][3] = II / v0 * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][7][0] = 0.1e1 / v0 * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][7][1] = II / v0 * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[4][5][2] = 0.1e1 / v0 * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[4][5][3] = II / v0 * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][4][2] = 0.1e1 / v0 * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][4][3] = II / v0 * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][8][0] = 0.1e1 / v0 * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][8][1] = II / v0 * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[6][1][0] = 0.1e1 / v0 * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[6][1][1] = II / v0 * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[7][3][0] = 0.1e1 / v0 * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[7][3][1] = II / v0 * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[8][5][0] = 0.1e1 / v0 * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[8][5][1] = II / v0 * SMConstants.C_MassTau;

  Curvature_Quark_F2H1[0][6][2] = 0.1e1 / v0 * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[0][6][3] = -II / v0 * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[0][9][0] =
      -0.1e1 / v0 * SMConstants.C_MassUp * conj(V11);
  Curvature_Quark_F2H1[0][9][1] = II / v0 * SMConstants.C_MassUp * conj(V11);
  Curvature_Quark_F2H1[0][10][0] =
      -0.1e1 / v0 * SMConstants.C_MassUp * conj(V12);
  Curvature_Quark_F2H1[0][10][1] = II / v0 * SMConstants.C_MassUp * conj(V12);
  Curvature_Quark_F2H1[0][11][0] =
      -0.1e1 / v0 * SMConstants.C_MassUp * conj(V13);
  Curvature_Quark_F2H1[0][11][1] = II / v0 * SMConstants.C_MassUp * conj(V13);
  Curvature_Quark_F2H1[1][7][2]  = 0.1e1 / v0 * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[1][7][3]  = -II / v0 * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[1][9][0] =
      -0.1e1 / v0 * SMConstants.C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[1][9][1] = II / v0 * SMConstants.C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[1][10][0] =
      -0.1e1 / v0 * SMConstants.C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[1][10][1] =
      II / v0 * SMConstants.C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[1][11][0] =
      -0.1e1 / v0 * SMConstants.C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[1][11][1] =
      II / v0 * SMConstants.C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[2][8][2] = 0.1e1 / v0 * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[2][8][3] = -II / v0 * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[2][9][0] =
      -0.1e1 / v0 * SMConstants.C_MassTop * conj(V31);
  Curvature_Quark_F2H1[2][9][1] = II / v0 * SMConstants.C_MassTop * conj(V31);
  Curvature_Quark_F2H1[2][10][0] =
      -0.1e1 / v0 * SMConstants.C_MassTop * conj(V32);
  Curvature_Quark_F2H1[2][10][1] = II / v0 * SMConstants.C_MassTop * conj(V32);
  Curvature_Quark_F2H1[2][11][0] =
      -0.1e1 / v0 * SMConstants.C_MassTop * conj(V33);
  Curvature_Quark_F2H1[2][11][1] = II / v0 * SMConstants.C_MassTop * conj(V33);
  Curvature_Quark_F2H1[3][6][0]  = 0.1e1 / v0 * SMConstants.C_MassDown * V11;
  Curvature_Quark_F2H1[3][6][1]  = II / v0 * SMConstants.C_MassDown * V11;
  Curvature_Quark_F2H1[3][7][0]  = V21 / v0 * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][7][1]  = II * V21 / v0 * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][8][0]  = 0.1e1 / v0 * SMConstants.C_MassDown * V31;
  Curvature_Quark_F2H1[3][8][1]  = II / v0 * SMConstants.C_MassDown * V31;
  Curvature_Quark_F2H1[3][9][2]  = 0.1e1 / v0 * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][9][3]  = II / v0 * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[4][6][0]  = 0.1e1 / v0 * SMConstants.C_MassStrange * V12;
  Curvature_Quark_F2H1[4][6][1]  = II / v0 * SMConstants.C_MassStrange * V12;
  Curvature_Quark_F2H1[4][7][0]  = V22 / v0 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][7][1]  = II * V22 / v0 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][8][0]  = 0.1e1 / v0 * SMConstants.C_MassStrange * V32;
  Curvature_Quark_F2H1[4][8][1]  = II / v0 * SMConstants.C_MassStrange * V32;
  Curvature_Quark_F2H1[4][10][2] = 0.1e1 / v0 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][10][3] = II / v0 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[5][6][0]  = V13 / v0 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][6][1]  = II / v0 * SMConstants.C_MassBottom * V13;
  Curvature_Quark_F2H1[5][7][0]  = V23 / v0 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][7][1]  = II / v0 * SMConstants.C_MassBottom * V23;
  Curvature_Quark_F2H1[5][8][0]  = V33 / v0 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][8][1]  = II / v0 * SMConstants.C_MassBottom * V33;
  Curvature_Quark_F2H1[5][11][2] = 0.1e1 / v0 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][11][3] = II / v0 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[6][0][2]  = 0.1e1 / v0 * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[6][0][3]  = -II / v0 * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[6][3][0]  = 0.1e1 / v0 * SMConstants.C_MassDown * V11;
  Curvature_Quark_F2H1[6][3][1]  = II / v0 * SMConstants.C_MassDown * V11;
  Curvature_Quark_F2H1[6][4][0]  = 0.1e1 / v0 * SMConstants.C_MassStrange * V12;
  Curvature_Quark_F2H1[6][4][1]  = II / v0 * SMConstants.C_MassStrange * V12;
  Curvature_Quark_F2H1[6][5][0]  = V13 / v0 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[6][5][1]  = II / v0 * SMConstants.C_MassBottom * V13;
  Curvature_Quark_F2H1[7][1][2]  = 0.1e1 / v0 * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[7][1][3]  = -II / v0 * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[7][3][0]  = V21 / v0 * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[7][3][1]  = II * V21 / v0 * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[7][4][0]  = V22 / v0 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[7][4][1]  = II * V22 / v0 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[7][5][0]  = V23 / v0 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[7][5][1]  = II / v0 * SMConstants.C_MassBottom * V23;
  Curvature_Quark_F2H1[8][2][2]  = 0.1e1 / v0 * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[8][2][3]  = -II / v0 * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[8][3][0]  = 0.1e1 / v0 * SMConstants.C_MassDown * V31;
  Curvature_Quark_F2H1[8][3][1]  = II / v0 * SMConstants.C_MassDown * V31;
  Curvature_Quark_F2H1[8][4][0]  = 0.1e1 / v0 * SMConstants.C_MassStrange * V32;
  Curvature_Quark_F2H1[8][4][1]  = II / v0 * SMConstants.C_MassStrange * V32;
  Curvature_Quark_F2H1[8][5][0]  = V33 / v0 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[8][5][1]  = II / v0 * SMConstants.C_MassBottom * V33;
  Curvature_Quark_F2H1[9][0][0] =
      -0.1e1 / v0 * SMConstants.C_MassUp * conj(V11);
  Curvature_Quark_F2H1[9][0][1] = II / v0 * SMConstants.C_MassUp * conj(V11);
  Curvature_Quark_F2H1[9][1][0] =
      -0.1e1 / v0 * SMConstants.C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[9][1][1] = II / v0 * SMConstants.C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[9][2][0] =
      -0.1e1 / v0 * SMConstants.C_MassTop * conj(V31);
  Curvature_Quark_F2H1[9][2][1] = II / v0 * SMConstants.C_MassTop * conj(V31);
  Curvature_Quark_F2H1[9][3][2] = 0.1e1 / v0 * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[9][3][3] = II / v0 * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[10][0][0] =
      -0.1e1 / v0 * SMConstants.C_MassUp * conj(V12);
  Curvature_Quark_F2H1[10][0][1] = II / v0 * SMConstants.C_MassUp * conj(V12);
  Curvature_Quark_F2H1[10][1][0] =
      -0.1e1 / v0 * SMConstants.C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[10][1][1] =
      II / v0 * SMConstants.C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[10][2][0] =
      -0.1e1 / v0 * SMConstants.C_MassTop * conj(V32);
  Curvature_Quark_F2H1[10][2][1] = II / v0 * SMConstants.C_MassTop * conj(V32);
  Curvature_Quark_F2H1[10][4][2] = 0.1e1 / v0 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[10][4][3] = II / v0 * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[11][0][0] =
      -0.1e1 / v0 * SMConstants.C_MassUp * conj(V13);
  Curvature_Quark_F2H1[11][0][1] = II / v0 * SMConstants.C_MassUp * conj(V13);
  Curvature_Quark_F2H1[11][1][0] =
      -0.1e1 / v0 * SMConstants.C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[11][1][1] =
      II / v0 * SMConstants.C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[11][2][0] =
      -0.1e1 / v0 * SMConstants.C_MassTop * conj(V33);
  Curvature_Quark_F2H1[11][2][1] = II / v0 * SMConstants.C_MassTop * conj(V33);
  Curvature_Quark_F2H1[11][5][2] = 0.1e1 / v0 * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[11][5][3] = II / v0 * SMConstants.C_MassBottom;

  SetCurvatureDone = true;
}

bool Class_SM::CalculateDebyeSimplified()
{
  return false;
}

bool Class_SM::CalculateDebyeGaugeSimplified()
{
  return false;
}
double Class_SM::VTreeSimplified(const std::vector<double> &v) const
{
  double vtmp = v[2];

  double res = 0;

  res += 0.5 * muSq * std::pow(vtmp, 2);
  res += 0.25 * lambda * std::pow(vtmp, 4);

  return res;
}

double Class_SM::VCounterSimplified(const std::vector<double> &v) const
{
  if (not UseVCounterSimplified) return 0;

  double vtmp = v[2];
  double res  = 0;

  res += 0.5 * dmuSq * std::pow(vtmp, 2);
  res += 0.25 * dlambda * std::pow(vtmp, 4);
  res += dT3 * vtmp;

  return res;
}

void Class_SM::Debugging(const std::vector<double> &input,
                         std::vector<double> &output) const
{
  (void)input;
  (void)output;
}
} // namespace Models
} // namespace BSMPT
