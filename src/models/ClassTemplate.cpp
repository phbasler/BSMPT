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

#include <BSMPT/models/ClassTemplate.h>
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
Class_Template::Class_Template(const ISMConstants &smConstants)
    : Class_Potential_Origin(smConstants)
{
  Model =
      ModelID::ModelIDs::TEMPLATE; // global int constant which will be used to
                                   // tell the program which model is called
  NNeutralHiggs = 1;               // number of neutral Higgs bosons at T = 0
  NChargedHiggs = 0; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar   = 2; // number of parameters in the tree-Level Lagrangian
  nParCT = 3; // number of parameters in the counterterm potential

  nVEV = 1; // number of VEVs to minimize the potential

  NHiggs = NNeutralHiggs + NChargedHiggs;

  VevOrder.resize(nVEV);
  // Here you have to tell which scalar field gets which VEV.
  VevOrder[0] = 0;

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_Template::~Class_Template()
{
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_Template::addLegendCT() const
{
  std::vector<std::string> labels;
  labels.push_back("dT");
  labels.push_back("dlambda");
  labels.push_back("dmsquared");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * input file
 */
std::vector<std::string> Class_Template::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c"); // Label for the critical temperature
  labels.push_back("v_c"); // Label for the critical vev
  labels.push_back(
      "v_c/T_c"); // Label for v_c/T_c, you could use xi_c also for example
  // out += "Your VEV order"; // Now you have to put the label for your vevs
  labels.push_back("omega");
  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * Higgs couplings. Use this to complement the legend of the given input file
 *
 */
std::vector<std::string> Class_Template::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;
  particles.resize(NHiggs);
  // here you have to define the particle names in the vector particles

  particles[0] = "H";

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
std::vector<std::string> Class_Template::addLegendVEV() const
{
  std::vector<std::string> labels;
  // out = "Your VEV order";
  labels.push_back("omega");
  return labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_Template::ReadAndSet(const std::string &linestr,
                                std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  double lms{0}, llambda{0};

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 2; k++)
  {
    ss >> tmp;
    if (k == 1)
      lms = tmp;
    else if (k == 2)
      llambda = tmp;
  }
  par[0] = lms;
  par[1] = llambda;

  set_gen(par); // This you have to call so that everything will be set
  return;
}

/**
 * Set Class Object as well as the VEV configuration
 */
void Class_Template::set_gen(const std::vector<double> &par)
{
  ms     = par[0]; // Class member is set accordingly to the input parameters
  lambda = par[1]; // Class member is set accordingly to the input parameters
  g      = SMConstants.C_g; // SM SU (2) gauge coupling --> SMparam .h
  yt     = std::sqrt(2) / SMConstants.C_vev0 *
       SMConstants.C_MassTop; // Top Yukawa coupling --> SMparam .h
  scale = SMConstants.C_vev0; // Renormalisation scale is set to the SM VEV
  vevTreeMin.resize(nVEV);
  vevTree.resize(NHiggs);
  // Here you have to set the vector vevTreeMin. The vector vevTree will then be
  // set by the function MinimizeOrderVEV
  vevTreeMin[0] = SMConstants.C_vev0;
  vevTree       = MinimizeOrderVEV(vevTreeMin);
  if (!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the
 * entries of Curvature_Higgs_CT_L1 to Curvature_Higgs_CT_L4.
 */
void Class_Template::set_CT_Pot_Par(const std::vector<double> &par)
{

  dT      = par[0];
  dlambda = par[1];
  dms     = par[2];

  Curvature_Higgs_CT_L1[0]          = dT;
  Curvature_Higgs_CT_L2[0][0]       = dms;
  Curvature_Higgs_CT_L4[0][0][0][0] = dlambda;
}

/**
 * console output of all Parameters
 */
void Class_Template::write() const
{

  std::stringstream ss;
  ss << "Model = " << Model << std::endl;

  ss << "The parameters are : " << std::endl;
  ss << "lambda = " << lambda << std::endl << "m^2 = " << ms << std::endl;

  ss << "The counterterm parameters are : " << std::endl;
  ss << "dT = " << dT << std::endl
     << "dlambda = " << dlambda << std::endl
     << "dm^2 = " << dms << std::endl;

  ss << "The scale is given by mu = " << scale << " GeV " << std::endl;
  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms. Here you need to work out the scheme and
 * implement the formulas.
 */
std::vector<double> Class_Template::calc_CT() const
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

  // Here you have to use your formulae for the counterterm scheme
  double t = 0;
  parCT.push_back(t); // dT
  parCT.push_back(3.0 * t / std::pow(SMConstants.C_vev0, 3) +
                  3.0 / std::pow(SMConstants.C_vev0, 3) * NablaWeinberg(0) -
                  3.0 / std::pow(SMConstants.C_vev0, 2) *
                      HesseWeinberg(0, 0)); // dlambda
  parCT.push_back(-3.0 / (2 * SMConstants.C_vev0) * NablaWeinberg(0) +
                  1.0 / 2.0 * HesseWeinberg(0, 0) -
                  3.0 * t / (2 * SMConstants.C_vev0)); // dms

  return parCT;
}

/**
 * Ensures the correct rotation matrix convention
 */
void Class_Template::AdjustRotationMatrix()
{
  // Here you implement the rotation matrix convention of your model
  // and define HiggsRotationMatrixEnsuredConvention, use then
  // HiggsRotationMatrixEnsuredConvention in TripleHiggsCouplings
}

void Class_Template::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsDone) CalculatePhysicalCouplings();

  if (CalculatedTripleCopulings) return;
  CalculatedTripleCopulings = true;

  std::vector<double> HiggsOrder(NHiggs);
  // Here you have to set the vector HiggsOrder. By telling e.g. HiggsOrder[0] =
  // 5 you always want your 6th lightest particle to be the first particle in
  // the vector (which has the index 5 because they are sorted by mass)

  // example for keeping the mass order
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

void Class_Template::SetCurvatureArrays()
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

  Curvature_Higgs_L2[0][0]       = ms;
  Curvature_Higgs_L4[0][0][0][0] = lambda;

  Curvature_Gauge_G2H2[0][0][0][0] = 4 * std::pow(g, 2);

  Curvature_Quark_F2H1[1][0][0] = yt;
  Curvature_Quark_F2H1[0][1][0] = yt;
}

bool Class_Template::CalculateDebyeSimplified()
{
  return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool Class_Template::CalculateDebyeGaugeSimplified()
{
  /*
   * Use this function if you calculated the Debye corrections to the gauge mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeGauge[NGauge][NGauge]
   */

  return false;
}
double Class_Template::VTreeSimplified(const std::vector<double> &v) const
{
  if (not UseVTreeSimplified) return 0;
  double res = 0;

  double vIn = v[0];
  res = 0.5 * ms * std::pow(vIn, 2) + 1.0 / 24.0 * lambda * std::pow(vIn, 4);

  return res;
}

double Class_Template::VCounterSimplified(const std::vector<double> &v) const
{
  if (not UseVCounterSimplified) return 0;
  double res = 0;
  double vIn = v[0];
  res = 0.5 * dms * std::pow(vIn, 2) + 1.0 / 24.0 * dlambda * std::pow(vIn, 4) +
        dT * vIn;
  return res;
}

void Class_Template::Debugging(const std::vector<double> &input,
                               std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
