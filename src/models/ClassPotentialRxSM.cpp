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
#include <cstddef>               // for std::size_t

#include <BSMPT/models/ClassPotentialRxSM.h>
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
Class_RxSM::Class_RxSM(const ISMConstants &smConstants)
    : Class_Potential_Origin(smConstants)
{
  Model = ModelID::ModelIDs::RXSM; // global int constant which will be used to
                                   // tell the program which model is called
  NNeutralHiggs = 3; // number of neutral Higgs bosons at T = 0
  NChargedHiggs = 2; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar   = 4; // number of independent input parameters (in the tree-Level Lagrangian)
  nParCT = 10; // number of parameters in the counterterm potential

  nVEV = 2; // number of VEVs to minimize the potential

  NHiggs = NNeutralHiggs + NChargedHiggs;

  VevOrder.resize(nVEV);
  // Here you have to tell which scalar field gets which VEV.
  VevOrder[0] = 3;
  VevOrder[1] = 4;

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_RxSM::~Class_RxSM()
{
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_RxSM::addLegendCT() const
{
  std::vector<std::string> labels;
  labels.push_back("dmsq");
  labels.push_back("dlambda");
  labels.push_back("dmSsq");
  labels.push_back("dlambdaS");
  labels.push_back("dlambdaHS");
  labels.push_back("dT1");
  labels.push_back("dT2");
  labels.push_back("dT3");
  labels.push_back("dT4");
  labels.push_back("dT5");

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * input file
 */
std::vector<std::string> Class_RxSM::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c"); // Label for the critical temperature
  labels.push_back("omega_c"); // Label for the critical vev
  labels.push_back("omega_c/T_c");
  labels.push_back("omega_H(T_c)");
  labels.push_back("omega_S(T_c)");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * Higgs couplings. Use this to complement the legend of the given input file
 *
 */
std::vector<std::string> Class_RxSM::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;
  particles.resize(NHiggs);
  // here you have to define the particle names in the vector particles

  particles[0] = "G+";
  particles[1] = "G-";
  particles[2] = "G0";
  particles[3] = "h_SM";
  particles[4] = "h_H";

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
std::vector<std::string> Class_RxSM::addLegendVEV() const
{
  std::vector<std::string> labels;
  // out = "Your VEV order";
  labels.push_back("omega_H");
  labels.push_back("omega_S");
  return labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_RxSM::ReadAndSet(const std::string &linestr,
                                std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  double lambdaSIn{0}, lambdaHSIn{0}, vSIn{0}, MassSIn{0};

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 4; ++k)
  {
    ss >> tmp;
    if (k == 1)
      lambdaSIn = tmp;
    else if (k == 2)
      lambdaHSIn = tmp;
    else if (k == 3)
      vSIn = tmp;
    else if (k == 4)
      MassSIn = tmp;
  }
  par[0] = lambdaSIn;
  par[1] = lambdaHSIn;
  par[2] = vSIn;
  par[3] = MassSIn;

  set_gen(par); // This you have to call so that everything will be set
  return;

}

/**
 * Set Class Object as well as the VEV configuration
 */
void Class_RxSM::set_gen(const std::vector<double> &par)
{
  lambdaS      = par[0];
  lambdaHS     = par[1];
  vS           = par[2];
  double MassS = par[3]; // this input is only used if vS == 0

  const double ZeroThreshold = 1e-5;

  UnbrokenSingletPhase = false;
  if (std::abs(vS) < ZeroThreshold)
  {
    UnbrokenSingletPhase = true;
    vS = 0.;
  }

  vH = SMConstants.C_vev0;

  double mh2 = SMConstants.C_MassSMHiggs*SMConstants.C_MassSMHiggs;

  // Fix lambda via the condition that one of the Higgses must be SM-like
  // Computed with Mathematica from diagonalising the 2x2 mass matrix and
  // then solving for lambda when fixing one mass to mh = 125.09 GeV

  // Be careful when the denominator = 0, i.e. when vS*vS*lambdaS == 3*mh2;
  // this is the point where the SM-like state switches from being the
  // heavier or lighter one to the corresponding other one
  if (std::abs(vS*vS*lambdaS - 3*mh2) < ZeroThreshold)
  {
    throw std::runtime_error("Invalid parameter point, lambda is divergent.");
  }

  lambda = 2*mh2/vH/vH + 6*vS*vS*lambdaHS*lambdaHS/(vS*vS*lambdaS - 3*mh2);

  // Tree-level tadpole equations
  msq  = -lambda*vH*vH/2.  - lambdaHS*vS*vS;

  if (UnbrokenSingletPhase)
  {
    mSsq = MassS*MassS - lambdaHS*vH*vH/2.;
  }
  else
  {
    mSsq = -lambdaS*vS*vS/6. - lambdaHS*vH*vH/2.;
  }

  scale = vH; // Renormalisation scale is set to the SM VEV

  vevTreeMin.resize(nVEV);
  vevTree.resize(NHiggs);
  vevTreeMin[0] = vH;
  vevTreeMin[1] = vS;
  vevTree = MinimizeOrderVEV(vevTreeMin);

  if (!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the
 * entries of Curvature_Higgs_CT_L1 to Curvature_Higgs_CT_L4.
 */
void Class_RxSM::set_CT_Pot_Par(const std::vector<double> &par)
{

  dmsq      = par[0];
  dlambda   = par[1];
  dmSsq     = par[2];
  dlambdaS  = par[3];
  dlambdaHS = par[4];
  dT1       = par[5];
  dT2       = par[6];
  dT3       = par[7];
  dT4       = par[8];
  dT5       = par[9];

  Curvature_Higgs_CT_L1[0] = dT1;
  Curvature_Higgs_CT_L1[1] = dT2;
  Curvature_Higgs_CT_L1[2] = dT3;
  Curvature_Higgs_CT_L1[3] = dT4;
  Curvature_Higgs_CT_L1[4] = dT5;

  Curvature_Higgs_CT_L2[0][0] = dmsq / 0.2e1;
  Curvature_Higgs_CT_L2[1][1] = dmsq / 0.2e1;
  Curvature_Higgs_CT_L2[2][2] = dmsq / 0.2e1;
  Curvature_Higgs_CT_L2[3][3] = dmsq / 0.2e1;
  Curvature_Higgs_CT_L2[4][4] = dmSsq;

  Curvature_Higgs_CT_L4[0][0][0][0] = 0.3e1 / 0.2e1 * dlambda;
  Curvature_Higgs_CT_L4[0][0][1][1] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][0][2][2] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][0][3][3] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][0][4][4] = dlambdaHS;
  Curvature_Higgs_CT_L4[0][1][0][1] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][1][1][0] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][0][2] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][2][2][0] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][0][3] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][3][3][0] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[0][4][0][4] = dlambdaHS;
  Curvature_Higgs_CT_L4[0][4][4][0] = dlambdaHS;
  Curvature_Higgs_CT_L4[1][0][0][1] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][0][1][0] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][1][0][0] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][1][1][1] = 0.3e1 / 0.2e1 * dlambda;
  Curvature_Higgs_CT_L4[1][1][2][2] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][1][3][3] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][1][4][4] = dlambdaHS;
  Curvature_Higgs_CT_L4[1][2][1][2] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][2][2][1] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][1][3] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][3][3][1] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[1][4][1][4] = dlambdaHS;
  Curvature_Higgs_CT_L4[1][4][4][1] = dlambdaHS;
  Curvature_Higgs_CT_L4[2][0][0][2] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][0][2][0] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][1][1][2] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][1][2][1] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][2][0][0] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][2][1][1] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][2][2][2] = 0.3e1 / 0.2e1 * dlambda;
  Curvature_Higgs_CT_L4[2][2][3][3] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][2][4][4] = dlambdaHS;
  Curvature_Higgs_CT_L4[2][3][2][3] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][3][3][2] = dlambda / 0.2e1;
  Curvature_Higgs_CT_L4[2][4][2][4] = dlambdaHS;
  Curvature_Higgs_CT_L4[2][4][4][2] = dlambdaHS;
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
  Curvature_Higgs_CT_L4[3][3][4][4] = dlambdaHS;
  Curvature_Higgs_CT_L4[3][4][3][4] = dlambdaHS;
  Curvature_Higgs_CT_L4[3][4][4][3] = dlambdaHS;
  Curvature_Higgs_CT_L4[4][0][0][4] = dlambdaHS;
  Curvature_Higgs_CT_L4[4][0][4][0] = dlambdaHS;
  Curvature_Higgs_CT_L4[4][1][1][4] = dlambdaHS;
  Curvature_Higgs_CT_L4[4][1][4][1] = dlambdaHS;
  Curvature_Higgs_CT_L4[4][2][2][4] = dlambdaHS;
  Curvature_Higgs_CT_L4[4][2][4][2] = dlambdaHS;
  Curvature_Higgs_CT_L4[4][3][3][4] = dlambdaHS;
  Curvature_Higgs_CT_L4[4][3][4][3] = dlambdaHS;
  Curvature_Higgs_CT_L4[4][4][0][0] = dlambdaHS;
  Curvature_Higgs_CT_L4[4][4][1][1] = dlambdaHS;
  Curvature_Higgs_CT_L4[4][4][2][2] = dlambdaHS;
  Curvature_Higgs_CT_L4[4][4][3][3] = dlambdaHS;
  Curvature_Higgs_CT_L4[4][4][4][4] = dlambdaS;
}

/**
 * console output of all Parameters
 */
void Class_RxSM::write() const
{

  std::stringstream ss;
  ss << "Model = " << Model << std::endl;

  ss << "The parameters are : " << std::endl;
  ss << "lambda    = " << lambda << " (fixed via requirement of "
     << "SM Higgs mass == 125.09 GeV)" << std::endl
     << "lambda_S  = " << lambdaS << std::endl
     << "lambda_HS = " << lambdaHS << std::endl
     << "v_S       = " << vS << std::endl
     << "v_H       = " << vH << " (fixed to SM value)" << std::endl
     << "m^2       = " << msq << " (via tadpole eqs.)" << std::endl
     << "ms^2      = " << mSsq << " (via tadpole eqs.)" << std::endl;

  ss << "The counterterm parameters are : " << std::endl;
  ss << "dm^2       = " << dmsq << std::endl
     << "dlambda    = " << dlambda << std::endl
     << "dms^2      = " << dmSsq << std::endl
     << "dlambda_S  = " << dlambdaS << std::endl
     << "dlambda_HS = " << dlambdaHS << std::endl
     << "dT1        = " << dT1 << std::endl
     << "dT2        = " << dT2 << std::endl
     << "dT3        = " << dT3 << std::endl
     << "dT4        = " << dT4 << std::endl
     << "dT5        = " << dT5 << std::endl;

  ss << "The scale is given by mu = " << scale << " GeV " << std::endl;

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

  ss << "The mass spectrum is given by :\n"
     << "m_{G^+}^2 = " << HiggsMasses[pos_G1] << " GeV \n"
     << "m_{G^-}^2 = " << HiggsMasses[pos_G2] << " GeV \n"
     << "m_{G^0}^2 = " << HiggsMasses[pos_G0] << " GeV \n"
     << "m_h_SM    = " << std::sqrt(HiggsMasses[pos_h_SM]) << " GeV \n"
     << "m_h_H     = " << std::sqrt(HiggsMasses[pos_h_H]) << " GeV \n";

  ss << "The neutral mixing Matrix is given by :\n";
  ss << "h_SM = " << HiggsRot(pos_h_SM, 3) << " zeta_1 ";
  bool IsNegative = HiggsRot(pos_h_SM, 4) < 0;
  if (IsNegative)
    ss << "- ";
  else
    ss << "+ ";
  ss << std::abs(HiggsRot(pos_h_SM, 4)) << " zeta_S\n"
     << "h_H  = " << HiggsRot(pos_h_H, 3) << " zeta_1 ";
  IsNegative = HiggsRot(pos_h_H, 4) < 0;
  if (IsNegative)
    ss << "- ";
  else
    ss << "+ ";
  ss << std::abs(HiggsRot(pos_h_H, 4)) << " zeta_S" << std::endl;
  ss << "The mixing angle is: alpha = " << alpha << std::endl;

  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms. Here you need to work out the scheme and
 * implement the formulas.
 */
std::vector<double> Class_RxSM::calc_CT() const
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

  if (UnbrokenSingletPhase)
  {
    // Free parameters t4 and t5 chosen such that dlambdaS = dlambdaHS = 0

    // dmsq
    parCT.push_back(-3*HesseWeinberg(2, 2) + HesseWeinberg(3, 3));

    // dlambda
    parCT.push_back((2*HesseWeinberg(2, 2) - 2*HesseWeinberg(3, 3))/vH/vH);

    // dmSsq
    parCT.push_back(-HesseWeinberg(4, 4));

    // dlambdaS
    parCT.push_back(0);

    // dlambdaHS
    parCT.push_back(0);

    // dT1
    parCT.push_back(-NablaWeinberg(0));

    // dT2
    parCT.push_back(-NablaWeinberg(1));

    // dT3
    parCT.push_back(-NablaWeinberg(2));

    // dT4
    parCT.push_back(HesseWeinberg(2, 2)*vH - NablaWeinberg(3));

    // dT5
    parCT.push_back(-NablaWeinberg(4));
  }
  else
  {
    // Free parameter t chosen such that all tadpole CTs vanish

    // dmsq
    parCT.push_back((-3*HesseWeinberg(2, 2)*vH
                     + HesseWeinberg(3, 3)*vH
                     + HesseWeinberg(3, 4)*vS)/vH);

    // dlambda
    parCT.push_back((2*HesseWeinberg(2, 2)
                     - 2*HesseWeinberg(3, 3))/vH/vH);

    // dmSsq
    // parCT.push_back((HesseWeinberg(3, 4)*vH/2
    //                  - HesseWeinberg(4, 4)*vS)/vS);
    parCT.push_back((HesseWeinberg(3, 4)*vH
                     + HesseWeinberg(4, 4)*vS
                     - 3*NablaWeinberg(4))/vS/2);

    // dlambdaS
    // parCT.push_back(0);
    parCT.push_back((-3*HesseWeinberg(4, 4)*vS
                     + 3*NablaWeinberg(4))*pow(vS, -3));

    // dlambdaHS
    parCT.push_back(-HesseWeinberg(3, 4)/vH/vS);

    // dT1
    parCT.push_back(-NablaWeinberg(0));

    // dT2
    parCT.push_back(-NablaWeinberg(1));

    // dT3
    parCT.push_back(-NablaWeinberg(2));

    // dT4
    parCT.push_back(HesseWeinberg(2, 2)*vH - NablaWeinberg(3));

    // dT5
    // parCT.push_back(HesseWeinberg(4, 4)*vS - NablaWeinberg(4));
    parCT.push_back(0);
  }

  return parCT;
}

/**
 * Ensures the correct rotation matrix convention
 */
void Class_RxSM::AdjustRotationMatrix()
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

  // initialize position indices (new initialization for each point in multiline
  // files)
  pos_G1 = -1, pos_G2 = -1, pos_G0 = -1, pos_h = -1, pos_H = -1;

  // Doublet: Phi1 = 1/Sqrt(2)*{rho1 + I*eta1, zeta1 + vH + I*psi1}
  // Singlet: S = zetaS + vS
  // higgsbasis = {rho1, eta1, psi1, zeta1, zetaS}

  // interaction basis
  // rho1, eta1, psi1, zeta1, zetaS
  int pos_rho1 = 0, pos_eta1 = 1, pos_psi1 = 2,
      pos_zeta1 = 3, pos_zetaS = 4;

  for (std::size_t i = 0; i < NHiggs;
       i++) // mass base index i corresponds to mass vector sorted in ascending
            // mass
  {
    // The Goldstones are all on the diagonal, there is no mixing
    if (std::abs(HiggsRot(i, pos_rho1)) > ZeroThreshold)
    {
      pos_G1 = i;
    }
    else if (std::abs(HiggsRot(i, pos_eta1)) > ZeroThreshold)
    {
      pos_G2 = i;
    }
    else if (std::abs(HiggsRot(i, pos_psi1)) > ZeroThreshold)
    {
      pos_G0 = i;
    }
    else if (std::abs(HiggsRot(i, pos_zeta1)) + std::abs(HiggsRot(i, pos_zetaS)) >
        ZeroThreshold) // use that mh < mH
    {
      if (pos_h == -1)
      {
        pos_h = i;
      }
      else
      {
        pos_H = i;
      }
    }
  }

  // check if all position indices are set
  if (pos_G1 == -1 or pos_G2 == -1 or pos_G0 == -1 or
      pos_h == -1 or pos_H == -1)
  {
    throw std::runtime_error("Error. Not all position indices are set.");
  }

  // check if all other elements of rotation matrix are zero
  bool zero_element = false;
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      int ii = int(i);
      int jj = int(j);
      if (not((jj == pos_rho1 and ii == pos_G1) or
              (jj == pos_eta1 and ii == pos_G2) or
              (jj == pos_psi1 and ii == pos_G0) or
              (jj == pos_zeta1 and (ii == pos_h or ii == pos_H)) or
              (jj == pos_zetaS and (ii == pos_h or ii == pos_H))))
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
  // and exotic Higgses
  pos_h_SM = -1, pos_h_H = -1;

  std::vector<double> HiggsMasses;
  HiggsMasses = HiggsMassesSquared(vevTree, 0);

  // Due to the masses being ordered, we will always have
  //  HiggsMasses[pos_h] <= HiggsMasses[pos_H]
  double diff1 = std::abs(std::sqrt(HiggsMasses[pos_h])
                          - SMConstants.C_MassSMHiggs);
  double diff2 = std::abs(std::sqrt(HiggsMasses[pos_H])
                          - SMConstants.C_MassSMHiggs);
  if (diff1 < diff2)
  {
    pos_h_SM = pos_h;
    pos_h_H = pos_H;
  }
  else
  {
    pos_h_H = pos_h;
    pos_h_SM = pos_H;
  }

  MatrixXd HiggsRotFixed(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsRotFixed.row(i) = HiggsRot.row(i);
  }

  // charged submatrix
  if (HiggsRotFixed(pos_G1, pos_rho1) < 0) // G1 rho1 (= +1)
  {
    HiggsRotFixed.row(pos_G1) *= -1;
  }
  if (HiggsRotFixed(pos_G2, pos_eta1) < 0) // G2 eta1 (= +1)
  {
    HiggsRotFixed.row(pos_G2) *= -1;
  }

  // check neutral, CP-odd submatrix
  if (HiggsRotFixed(pos_G0, pos_psi1) < 0) // G0 psi1 (= +1)
  {
    HiggsRotFixed.row(pos_G0) *= -1;
  }

  // // check neutral, CP-even submatrix
  if (HiggsRotFixed(pos_h, pos_zeta1) < 0) // h zeta1 (+ cos(alpha))
  {
    // if negative, rotate h
    HiggsRotFixed.row(pos_h) *= -1; // h
  }
  if (HiggsRotFixed(pos_H, pos_zetaS) < 0) // H zetaS (+ cos(alpha))
  {
    // if negative, rotate H
    HiggsRotFixed.row(pos_H) *= -1; // H
  }

  // Extract the fixed mixing angle
  alpha = std::asin(HiggsRotFixed(pos_h, pos_zetaS)); // h zetaS (+ sin(alpha))

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRotationMatrixEnsuredConvention[i][j] = HiggsRotFixed(i, j);
    }
  }

  return;
}

void Class_RxSM::TripleHiggsCouplings()
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

  std::vector<double> HiggsOrder(NHiggs);
  HiggsOrder[0] = pos_G1;
  HiggsOrder[1] = pos_G2;
  HiggsOrder[2] = pos_G0;
  HiggsOrder[3] = pos_h_SM;
  HiggsOrder[4] = pos_h_H;

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

void Class_RxSM::SetCurvatureArrays()
{
  initVectors();
  SetCurvatureDone = true;
  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

  Curvature_Higgs_L2[0][0] = msq / 0.2e1;
  Curvature_Higgs_L2[1][1] = msq / 0.2e1;
  Curvature_Higgs_L2[2][2] = msq / 0.2e1;
  Curvature_Higgs_L2[3][3] = msq / 0.2e1;
  Curvature_Higgs_L2[4][4] = mSsq;

  Curvature_Higgs_L4[0][0][0][0] = 0.3e1 / 0.2e1 * lambda;
  Curvature_Higgs_L4[0][0][1][1] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][0][2][2] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][0][3][3] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][0][4][4] = lambdaHS;
  Curvature_Higgs_L4[0][1][0][1] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][1][1][0] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][2][0][2] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][2][2][0] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][3][0][3] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][3][3][0] = lambda / 0.2e1;
  Curvature_Higgs_L4[0][4][0][4] = lambdaHS;
  Curvature_Higgs_L4[0][4][4][0] = lambdaHS;
  Curvature_Higgs_L4[1][0][0][1] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][0][1][0] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][1][0][0] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][1][1][1] = 0.3e1 / 0.2e1 * lambda;
  Curvature_Higgs_L4[1][1][2][2] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][1][3][3] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][1][4][4] = lambdaHS;
  Curvature_Higgs_L4[1][2][1][2] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][2][2][1] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][3][1][3] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][3][3][1] = lambda / 0.2e1;
  Curvature_Higgs_L4[1][4][1][4] = lambdaHS;
  Curvature_Higgs_L4[1][4][4][1] = lambdaHS;
  Curvature_Higgs_L4[2][0][0][2] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][0][2][0] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][1][1][2] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][1][2][1] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][2][0][0] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][2][1][1] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][2][2][2] = 0.3e1 / 0.2e1 * lambda;
  Curvature_Higgs_L4[2][2][3][3] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][2][4][4] = lambdaHS;
  Curvature_Higgs_L4[2][3][2][3] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][3][3][2] = lambda / 0.2e1;
  Curvature_Higgs_L4[2][4][2][4] = lambdaHS;
  Curvature_Higgs_L4[2][4][4][2] = lambdaHS;
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
  Curvature_Higgs_L4[3][3][4][4] = lambdaHS;
  Curvature_Higgs_L4[3][4][3][4] = lambdaHS;
  Curvature_Higgs_L4[3][4][4][3] = lambdaHS;
  Curvature_Higgs_L4[4][0][0][4] = lambdaHS;
  Curvature_Higgs_L4[4][0][4][0] = lambdaHS;
  Curvature_Higgs_L4[4][1][1][4] = lambdaHS;
  Curvature_Higgs_L4[4][1][4][1] = lambdaHS;
  Curvature_Higgs_L4[4][2][2][4] = lambdaHS;
  Curvature_Higgs_L4[4][2][4][2] = lambdaHS;
  Curvature_Higgs_L4[4][3][3][4] = lambdaHS;
  Curvature_Higgs_L4[4][3][4][3] = lambdaHS;
  Curvature_Higgs_L4[4][4][0][0] = lambdaHS;
  Curvature_Higgs_L4[4][4][1][1] = lambdaHS;
  Curvature_Higgs_L4[4][4][2][2] = lambdaHS;
  Curvature_Higgs_L4[4][4][3][3] = lambdaHS;
  Curvature_Higgs_L4[4][4][4][4] = lambdaS;

  Curvature_Gauge_G2H2[0][0][0][0] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][0][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][0][3] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][1][2] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][2][1] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[0][3][3][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][0][0] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][1][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][0][2] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][1][3] = -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][2][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[1][3][3][1] = -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][0][0] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][1][1] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][2][2] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][2][3][3] = SMConstants.C_g * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][0][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][1][1] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][2][2] = -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[2][3][3][3] = -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][0][3] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][1][2] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][2][1] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][0][3][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][0][2] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][1][3] = -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][2][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][1][3][1] = -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][0][0] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][1][1] = SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][2][2] = -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][2][3][3] = -SMConstants.C_gs * SMConstants.C_g / 0.2e1;
  Curvature_Gauge_G2H2[3][3][0][0] = SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][1][1] = SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][2][2] = SMConstants.C_gs * SMConstants.C_gs / 0.2e1;
  Curvature_Gauge_G2H2[3][3][3][3] = SMConstants.C_gs * SMConstants.C_gs / 0.2e1;

  Curvature_Lepton_F2H1[0][1][2] = II / vH * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[0][1][3] = 0.1e1 / vH * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][0][2] = II / vH * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][0][3] = 0.1e1 / vH * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][6][0] = 0.1e1 / vH * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[1][6][1] = II / vH * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[2][3][2] = II / vH * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[2][3][3] = 0.1e1 / vH * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][2][2] = II / vH * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][2][3] = 0.1e1 / vH * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][7][0] = 0.1e1 / vH * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[3][7][1] = II / vH * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[4][5][2] = II / vH * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[4][5][3] = 0.1e1 / vH * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][4][2] = II / vH * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][4][3] = 0.1e1 / vH * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][8][0] = 0.1e1 / vH * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[5][8][1] = II / vH * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[6][1][0] = 0.1e1 / vH * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[6][1][1] = II / vH * SMConstants.C_MassElectron;
  Curvature_Lepton_F2H1[7][3][0] = 0.1e1 / vH * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[7][3][1] = II / vH * SMConstants.C_MassMu;
  Curvature_Lepton_F2H1[8][5][0] = 0.1e1 / vH * SMConstants.C_MassTau;
  Curvature_Lepton_F2H1[8][5][1] = II / vH * SMConstants.C_MassTau;


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

  Curvature_Quark_F2H1[0][6][2] = -II / vH * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[0][6][3] = 0.1e1 / vH * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[0][9][0] = -0.1e1 / vH * SMConstants.C_MassUp * conj(V11);
  Curvature_Quark_F2H1[0][9][1] = II / vH * SMConstants.C_MassUp * conj(V11);
  Curvature_Quark_F2H1[0][10][0] = -0.1e1 / vH * SMConstants.C_MassUp * conj(V12);
  Curvature_Quark_F2H1[0][10][1] = II / vH * SMConstants.C_MassUp * conj(V12);
  Curvature_Quark_F2H1[0][11][0] = -0.1e1 / vH * SMConstants.C_MassUp * conj(V13);
  Curvature_Quark_F2H1[0][11][1] = II / vH * SMConstants.C_MassUp * conj(V13);
  Curvature_Quark_F2H1[1][7][2] = -II / vH * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[1][7][3] = 0.1e1 / vH * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[1][9][0] = -0.1e1 / vH * SMConstants.C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[1][9][1] = II / vH * SMConstants.C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[1][10][0] = -0.1e1 / vH * SMConstants.C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[1][10][1] = II / vH * SMConstants.C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[1][11][0] = -0.1e1 / vH * SMConstants.C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[1][11][1] = II / vH * SMConstants.C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[2][8][2] = -II / vH * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[2][8][3] = 0.1e1 / vH * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[2][9][0] = -0.1e1 / vH * SMConstants.C_MassTop * conj(V31);
  Curvature_Quark_F2H1[2][9][1] = II / vH * SMConstants.C_MassTop * conj(V31);
  Curvature_Quark_F2H1[2][10][0] = -0.1e1 / vH * SMConstants.C_MassTop * conj(V32);
  Curvature_Quark_F2H1[2][10][1] = II / vH * SMConstants.C_MassTop * conj(V32);
  Curvature_Quark_F2H1[2][11][0] = -0.1e1 / vH * SMConstants.C_MassTop * conj(V33);
  Curvature_Quark_F2H1[2][11][1] = II / vH * SMConstants.C_MassTop * conj(V33);
  Curvature_Quark_F2H1[3][6][0] = V11 / vH * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][6][1] = II * V11 / vH * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][7][0] = V21 / vH * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][7][1] = II * V21 / vH * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][8][0] = V31 / vH * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][8][1] = II * V31 / vH * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][9][2] = II / vH * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[3][9][3] = 0.1e1 / vH * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[4][6][0] = 0.1e1 / vH * SMConstants.C_MassStrange * V12;
  Curvature_Quark_F2H1[4][6][1] = II / vH * SMConstants.C_MassStrange * V12;
  Curvature_Quark_F2H1[4][7][0] = V22 / vH * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][7][1] = II * V22 / vH * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][8][0] = 0.1e1 / vH * SMConstants.C_MassStrange * V32;
  Curvature_Quark_F2H1[4][8][1] = II / vH * SMConstants.C_MassStrange * V32;
  Curvature_Quark_F2H1[4][10][2] = II / vH * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[4][10][3] = 0.1e1 / vH * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[5][6][0] = V13 / vH * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][6][1] = II / vH * SMConstants.C_MassBottom * V13;
  Curvature_Quark_F2H1[5][7][0] = V23 / vH * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][7][1] = II / vH * SMConstants.C_MassBottom * V23;
  Curvature_Quark_F2H1[5][8][0] = V33 / vH * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][8][1] = II / vH * SMConstants.C_MassBottom * V33;
  Curvature_Quark_F2H1[5][11][2] = II / vH * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[5][11][3] = 0.1e1 / vH * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[6][0][2] = -II / vH * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[6][0][3] = 0.1e1 / vH * SMConstants.C_MassUp;
  Curvature_Quark_F2H1[6][3][0] = V11 / vH * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[6][3][1] = II * V11 / vH * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[6][4][0] = 0.1e1 / vH * SMConstants.C_MassStrange * V12;
  Curvature_Quark_F2H1[6][4][1] = II / vH * SMConstants.C_MassStrange * V12;
  Curvature_Quark_F2H1[6][5][0] = V13 / vH * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[6][5][1] = II / vH * SMConstants.C_MassBottom * V13;
  Curvature_Quark_F2H1[7][1][2] = -II / vH * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[7][1][3] = 0.1e1 / vH * SMConstants.C_MassCharm;
  Curvature_Quark_F2H1[7][3][0] = V21 / vH * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[7][3][1] = II * V21 / vH * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[7][4][0] = V22 / vH * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[7][4][1] = II * V22 / vH * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[7][5][0] = V23 / vH * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[7][5][1] = II / vH * SMConstants.C_MassBottom * V23;
  Curvature_Quark_F2H1[8][2][2] = -II / vH * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[8][2][3] = 0.1e1 / vH * SMConstants.C_MassTop;
  Curvature_Quark_F2H1[8][3][0] = V31 / vH * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[8][3][1] = II * V31 / vH * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[8][4][0] = 0.1e1 / vH * SMConstants.C_MassStrange * V32;
  Curvature_Quark_F2H1[8][4][1] = II / vH * SMConstants.C_MassStrange * V32;
  Curvature_Quark_F2H1[8][5][0] = V33 / vH * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[8][5][1] = II / vH * SMConstants.C_MassBottom * V33;
  Curvature_Quark_F2H1[9][0][0] = -0.1e1 / vH * SMConstants.C_MassUp * conj(V11);
  Curvature_Quark_F2H1[9][0][1] = II / vH * SMConstants.C_MassUp * conj(V11);
  Curvature_Quark_F2H1[9][1][0] = -0.1e1 / vH * SMConstants.C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[9][1][1] = II / vH * SMConstants.C_MassCharm * conj(V21);
  Curvature_Quark_F2H1[9][2][0] = -0.1e1 / vH * SMConstants.C_MassTop * conj(V31);
  Curvature_Quark_F2H1[9][2][1] = II / vH * SMConstants.C_MassTop * conj(V31);
  Curvature_Quark_F2H1[9][3][2] = II / vH * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[9][3][3] = 0.1e1 / vH * SMConstants.C_MassDown;
  Curvature_Quark_F2H1[10][0][0] = -0.1e1 / vH * SMConstants.C_MassUp * conj(V12);
  Curvature_Quark_F2H1[10][0][1] = II / vH * SMConstants.C_MassUp * conj(V12);
  Curvature_Quark_F2H1[10][1][0] = -0.1e1 / vH * SMConstants.C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[10][1][1] = II / vH * SMConstants.C_MassCharm * conj(V22);
  Curvature_Quark_F2H1[10][2][0] = -0.1e1 / vH * SMConstants.C_MassTop * conj(V32);
  Curvature_Quark_F2H1[10][2][1] = II / vH * SMConstants.C_MassTop * conj(V32);
  Curvature_Quark_F2H1[10][4][2] = II / vH * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[10][4][3] = 0.1e1 / vH * SMConstants.C_MassStrange;
  Curvature_Quark_F2H1[11][0][0] = -0.1e1 / vH * SMConstants.C_MassUp * conj(V13);
  Curvature_Quark_F2H1[11][0][1] = II / vH * SMConstants.C_MassUp * conj(V13);
  Curvature_Quark_F2H1[11][1][0] = -0.1e1 / vH * SMConstants.C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[11][1][1] = II / vH * SMConstants.C_MassCharm * conj(V23);
  Curvature_Quark_F2H1[11][2][0] = -0.1e1 / vH * SMConstants.C_MassTop * conj(V33);
  Curvature_Quark_F2H1[11][2][1] = II / vH * SMConstants.C_MassTop * conj(V33);
  Curvature_Quark_F2H1[11][5][2] = II / vH * SMConstants.C_MassBottom;
  Curvature_Quark_F2H1[11][5][3] = 0.1e1 / vH * SMConstants.C_MassBottom;

}

bool Class_RxSM::CalculateDebyeSimplified()
{
  return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool Class_RxSM::CalculateDebyeGaugeSimplified()
{
  /*
   * Use this function if you calculated the Debye corrections to the gauge mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeGauge[NGauge][NGauge]
   */

  return false;
}
double Class_RxSM::VTreeSimplified(const std::vector<double> &v) const
{
  (void)v;
  if (not UseVTreeSimplified) return 0;
  double res = 0;
  return res;
}

double Class_RxSM::VCounterSimplified(const std::vector<double> &v) const
{
  (void)v;
  if (not UseVCounterSimplified) return 0;
  double res = 0;
  return res;
}

void Class_RxSM::Debugging(const std::vector<double> &input,
                               std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
