// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gsl/gsl_sf_gamma.h>
#include <iomanip>
#include <random>

#include "Eigen/Dense"

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_zeta.h>

#include <BSMPT/ThermalFunctions/NegativeBosonSpline.h>
#include <BSMPT/ThermalFunctions/ThermalFunctions.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/modeltests/ModelTestfunctions.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
using namespace Eigen;

namespace BSMPT
{
Class_Potential_Origin::Class_Potential_Origin()
    : Class_Potential_Origin(GetSMConstants())
{
}

Class_Potential_Origin::Class_Potential_Origin(const ISMConstants &smConstants)
    : SMConstants{smConstants}
    , scale{SMConstants.C_vev0}

{
  // TODO Auto-generated constructor stub
}

Class_Potential_Origin::~Class_Potential_Origin()
{
  // TODO Auto-generated destructor stub
}

/**
 * This will call set_gen(par), SetCurvatureArrays, set_CT_Pot_Par(parCT),
 * CalculateDebye() as well as CalculateDebyeGauge()
 */
void Class_Potential_Origin::set_All(const std::vector<double> &par,
                                     const std::vector<double> &parCT)
{

  set_gen(par);
  if (!SetCurvatureDone) SetCurvatureArrays();
  set_CT_Pot_Par(parCT);
  CalculateDebye();
  CalculateDebyeGauge();
}

void Class_Potential_Origin::Prepare_Triple()
{
  for (std::size_t a = 0; a < NHiggs; a++)
  {
    for (std::size_t b = 0; b < NHiggs; b++)
    {
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        LambdaHiggs_3_CT[a][b][i] = Curvature_Higgs_CT_L3[a][b][i];
        for (std::size_t j = 0; j < NHiggs; j++)
        {
          LambdaHiggs_3_CT[a][b][i] +=
              Curvature_Higgs_CT_L4[a][b][i][j] * HiggsVev[j];
        }
      }
    }
  }
}

double Class_Potential_Origin::FCW(double MassSquared) const
{
  double res = 0;
  double x;
  double Boarder = std::pow(10, -200);
  if (std::isnan(MassSquared))
    x = Boarder;
  else if (std::abs(MassSquared) < Boarder)
    x = Boarder;
  else
    x = std::abs(MassSquared);
  res = std::log(x) - 2 * std::log(scale);
  return res;
}

double
Class_Potential_Origin::CWTerm(double MassSquared, double cb, int diff) const
{

  if (std::abs(MassSquared) < C_threshold) return 0;
  double LogTerm = 0, PotVal = 0;
  LogTerm = FCW(MassSquared);
  if (diff == 0)
    PotVal =
        1.0 / (64 * M_PI * M_PI) * MassSquared * MassSquared * (LogTerm - cb);
  else if (diff > 0)
  {
    PotVal = 1.0 / (32 * M_PI * M_PI) * MassSquared * (LogTerm - cb + 0.5);
  }
  return PotVal;
}

double Class_Potential_Origin::boson(double MassSquared,
                                     double Temp,
                                     double cb,
                                     int diff) const
{
  double res = 0;
  if (diff >= 0) res = CWTerm(std::abs(MassSquared), cb, diff);
  if (Temp == 0) return res;
  double Ratio = MassSquared / std::pow(Temp, 2);
  if (diff == 0)
  {
    res += std::pow(Temp, 4) / (2 * std::pow(M_PI, 2)) *
           ThermalFunctions::JbosonInterpolated(Ratio);
  }
  else if (diff == 1)
  {
    res += std::pow(Temp, 2) / (2 * std::pow(M_PI, 2)) *
           ThermalFunctions::JbosonNumericalIntegration(Ratio, 1);
  }
  else if (diff == -1)
  {
    res += 1.0 / (2 * std::pow(M_PI, 2)) *
           (4 * std::pow(Temp, 3) *
                ThermalFunctions::JbosonNumericalIntegration(Ratio, 0) -
            2 * Temp * MassSquared *
                ThermalFunctions::JbosonNumericalIntegration(Ratio, 1));
  }
  return res;
}

double
Class_Potential_Origin::fermion(double MassSquared, double Temp, int diff) const
{
  double res = 0;
  if (diff >= 0) res = CWTerm(std::abs(MassSquared), C_CWcbFermion, diff);
  double Ratio = MassSquared / std::pow(Temp, 2);
  if (Temp == 0) return res;
  if (diff == 0)
  {
    res += std::pow(Temp, 4) / (2 * std::pow(M_PI, 2)) *
           ThermalFunctions::JfermionInterpolated(Ratio);
  }
  else if (diff == 1)
  {
    res += std::pow(Temp, 2) / (2 * std::pow(M_PI, 2)) *
           ThermalFunctions::JfermionNumericalIntegration(Ratio, 1);
  }
  else if (diff == -1)
  {
    res += 1.0 / (2 * std::pow(M_PI, 2)) *
           (4 * std::pow(Temp, 3) *
                ThermalFunctions::JfermionNumericalIntegration(Ratio) -
            2 * Temp * MassSquared *
                ThermalFunctions::JfermionNumericalIntegration(Ratio, 1));
  }
  return res;
}

std::vector<double> Class_Potential_Origin::FirstDerivativeOfEigenvalues(
    const Ref<MatrixXcd> M,
    const Ref<MatrixXcd> MDiff) const
{
  std::vector<double> res;
  const std::size_t nRows = M.rows();
  const std::size_t nCols = M.cols();

  const double EVThres = std::pow(10, -6);

  if (nCols != nRows)
  {
    throw std::runtime_error("ERROR ! M needs to be an quadratic Matrix for "
                             "calculating the derivatives !\n");
  }

  const std::size_t nSize = nRows;

  SelfAdjointEigenSolver<MatrixXcd> es;
  es.compute(M);

  std::vector<std::complex<double>> Eigenvalues(nSize);
  std::vector<std::complex<double>> Derivatives(nSize);
  std::vector<double> AlreadyCalculated(
      nSize); // Array to check which EVs already been calculated.
  for (std::size_t i = 0; i < nSize; i++)
  {
    Eigenvalues[i] = es.eigenvalues()[i];
    if (std::abs(Eigenvalues[i]) < EVThres)
    {
      Eigenvalues[i] = 0;
    }
  }

  std::vector<std::vector<double>> Mapping(nSize, std::vector<double>(nSize));

  for (std::size_t i = 0; i < nSize; i++)
  {
    AlreadyCalculated[i] = -1;
    for (std::size_t j = i; j < nSize; j++)
    {
      if (std::abs(Eigenvalues[i] - Eigenvalues[j]) > EVThres)
      {
        Mapping[i][j] = 0;
      }
      else
      {
        Mapping[i][j] = 1;
      }
    }
  }
  for (std::size_t i = 1; i < nSize; i++)
  {
    for (std::size_t j = 0; j < i; j++)
      Mapping[i][j] = Mapping[j][i];
  }

  for (std::size_t p = 0; p < nSize; p++)
  {
    if (AlreadyCalculated[p] == -1)
    {
      std::size_t NumOfReps = 0;
      for (std::size_t i = p + 1; i < nSize; i++)
      {
        NumOfReps += Mapping[p][i];
      }
      if (NumOfReps == 0)
      {
        VectorXcd v(nSize);
        v = es.eigenvectors().col(p);
        // Derivatives[p] = (v.transpose()*MDiff*v).value();
        Derivatives[p]       = (v.adjoint() * MDiff * v).value();
        AlreadyCalculated[p] = 1;
      }
      else
      {
        MatrixXcd Phi(nSize, NumOfReps + 1);
        std::size_t helpCol = 0;
        MatrixXcd MXWork(NumOfReps + 1, NumOfReps + 1);

        for (std::size_t i = p; i < nSize; i++)
        {
          if (Mapping[p][i] == 1)
          {
            Phi.col(helpCol) = es.eigenvectors().col(i);
            helpCol++;
          }
        }
        MXWork = Phi.adjoint() * MDiff * Phi;
        SelfAdjointEigenSolver<MatrixXcd> esWork;
        esWork.compute(MXWork);
        helpCol = 0;
        for (std::size_t i = p; i < nSize; i++)
        {
          if (Mapping[p][i] == 1)
          {
            AlreadyCalculated[i] = 1;
            Derivatives[i]       = esWork.eigenvalues()[helpCol];
            helpCol++;
          }
        }
      }
    }
  }
  for (std::size_t i = 0; i < nSize; i++)
  {
    if (std::abs(Derivatives[i]) < EVThres) Derivatives[i] = 0;
  }

  for (std::size_t i = 0; i < nSize; i++)
  {
    res.push_back(Eigenvalues[i].real());
  }
  for (std::size_t i = 0; i < nSize; i++)
    res.push_back(Derivatives[i].real());
  return res;
}

double Class_Potential_Origin::fbaseTri(double MassSquaredA,
                                        double MassSquaredB,
                                        double MassSquaredC) const
{
  double res  = 0;
  double mas  = MassSquaredA;
  double mbs  = MassSquaredB;
  double mcs  = MassSquaredC;
  double LogA = 0, LogB = 0, LogC = 0;
  double Thres = 1e-8;
  if (std::abs(mas) < Thres) mas = 0;
  if (std::abs(mbs) < Thres) mbs = 0;
  if (std::abs(mcs) < Thres) mcs = 0;
  if (std::abs(mas - mbs) < Thres) mas = mbs;
  if (std::abs(mas - mcs) < Thres) mas = mcs;
  if (std::abs(mbs - mcs) < Thres) mbs = mcs;

  if (mas != 0) LogA = std::log(mas) - 2 * std::log(scale);
  if (mbs != 0) LogB = std::log(mbs) - 2 * std::log(scale);
  if (mcs != 0) LogC = std::log(mcs) - 2 * std::log(scale);

  std::size_t C = 1;
  if (mas == 0 and mbs == 0 and mcs == 0)
    res = 0;
  else if (mas != 0 and mbs == 0 and mcs == 0)
  {
    C   = 2;
    res = 1.0 / mas * (LogA - 1);
  }
  else if (mas == 0 and mbs != 0 and mcs == 0)
  {
    C   = 3;
    res = (LogB - 1) / mbs;
  }
  else if (mas == 0 and mbs == 0 and mcs != 0)
  {
    C   = 4;
    res = (LogC - 1) / mcs;
  }
  else if (mas == mbs and mas != 0 and mas != mcs and mcs != 0)
  {
    C   = 6;
    res = (mbs - mcs + mcs * std::log(mcs / mbs)) / std::pow(mbs - mcs, 2);
  }
  else if (mas == mcs and mas != 0 and mas != mbs and mbs != 0)
  {
    C   = 7;
    res = (mbs * log(mbs / mcs) - mbs + mcs) / std::pow(mbs - mcs, 2);
  }
  else if (mbs == mcs and mas != 0 and mbs != mas and mbs != 0)
  {
    C   = 8;
    res = (mas * std::log(mas / mcs) - mas + mcs) / std::pow(mas - mcs, 2);
  }
  else if (mas == mbs and mas == mcs and mas != 0)
  {
    C   = 9;
    res = 1.0 / (2 * mcs);
  }
  else if (mas == mbs and mas != mcs and mas != 0 and mcs == 0)
  {
    C   = 10;
    res = 1.0 / mbs;
  }
  else if (mas == mcs and mas != mbs and mas != 0 and mbs == 0)
  {
    C   = 11;
    res = 1.0 / mas;
  }
  else if (mbs == mcs and mbs != 0 and mbs != mas and mas == 0)
  {
    C   = 12;
    res = 1.0 / mbs;
  }
  else
  {
    C   = 5;
    res = mas * LogA / ((mas - mbs) * (mas - mcs)) +
          mbs * LogB / ((mbs - mas) * (mbs - mcs));
    res += mcs * LogC / ((mcs - mas) * (mcs - mbs));
  }

  if (std::isnan(res) or std::isinf(res))
  {
    std::string throwstring = "Found nan at line = ";
    throwstring += std::to_string(InputLineNumber);
    throwstring += " in function ";
    throwstring += __func__;
    throwstring += "\n";
    std::stringstream ss;
    ss << "Found nan at line = " << InputLineNumber << " in function "
       << __func__ << std::endl;
    ss << mas << sep << mbs << sep << mcs << sep << res << sep << C
       << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str(), __FILE__, __LINE__);
    throw std::runtime_error(throwstring.c_str());
  }

  return res;
}

double Class_Potential_Origin::fbaseFour(double MassSquaredA,
                                         double MassSquaredB,
                                         double MassSquaredC,
                                         double MassSquaredD) const
{

  double res  = 0;
  double mas  = MassSquaredA;
  double mbs  = MassSquaredB;
  double mcs  = MassSquaredC;
  double mds  = MassSquaredD;
  double LogA = 0, LogB = 0, LogC = 0, LogD = 0;
  double thresZero = 1e-6;
  double thresDeg  = 1e-6;
  if (std::abs(mas) < thresZero) mas = 0;
  if (std::abs(mbs) < thresZero) mbs = 0;
  if (std::abs(mcs) < thresZero) mcs = 0;
  if (std::abs(mds) < thresZero) mds = 0;
  if (std::abs(mas - mbs) < thresDeg) mbs = mas;
  if (std::abs(mas - mcs) < thresDeg) mcs = mas;
  if (std::abs(mas - mds) < thresDeg) mds = mas;
  if (std::abs(mbs - mcs) < thresDeg) mcs = mbs;
  if (std::abs(mbs - mds) < thresDeg) mds = mbs;
  if (std::abs(mcs - mds) < thresDeg) mds = mcs;

  if (mas != 0) LogA = std::log(mas) - 2 * std::log(scale);
  if (mbs != 0) LogB = std::log(mbs) - 2 * std::log(scale);
  if (mcs != 0) LogC = std::log(mcs) - 2 * std::log(scale);
  if (mds != 0) LogD = std::log(mds) - 2 * std::log(scale);

  std::size_t C = 0;

  // all masses are zero
  if (mas == 0 and mbs == 0 and mcs == 0 and mds == 0)
  {
    C = 1;
    // f0000
    res = 0;
  }

  // one mass is non-zero, the other ones are zero
  else if (mas != 0 and mbs == 0 and mcs == 0 and mds == 0)
  {
    C = 2;
    // fa000
    res = LogA / (mas * mas);
  }
  else if (mas == 0 and mbs != 0 and mcs == 0 and mds == 0)
  {
    C = 3;
    // f0b00
    res = LogB / (mbs * mbs);
  }
  else if (mas == 0 and mbs == 0 and mcs != 0 and mds == 0)
  {
    C = 4;
    // f00c0
    res = LogC / (mcs * mcs);
  }
  else if (mas == 0 and mbs == 0 and mcs == 0 and mds != 0)
  {
    C = 5;
    // f000d
    res = LogD / (mds * mds);
  }

  // two masses are non-zero, the other masses are zero
  // 1) they are equal
  else if (mas == mbs and mas != 0 and mcs == 0 and mds == 0)
  {
    C = 6;
    // faa00
    res = (2. - LogA) / (mas * mas);
  }
  else if (mas == mcs and mas != 0 and mbs == 0 and mds == 0)
  {
    C = 7;
    // fa0a0
    res = (2. - LogA) / (mas * mas);
  }
  else if (mas == mds and mas != 0 and mbs == 0 and mcs == 0)
  {
    C = 8;
    // fa00a
    res = (2. - LogA) / (mas * mas);
  }
  else if (mbs == mcs and mbs != 0 and mas == 0 and mds == 0)
  {
    C = 9;
    // f0bb0
    res = (2. - LogB) / (mbs * mbs);
  }
  else if (mbs == mds and mbs != 0 and mas == 0 and mcs == 0)
  {
    C = 10;
    // f0b0b
    res = (2. - LogB) / (mbs * mbs);
  }
  else if (mcs == mds and mcs != 0 and mas == 0 and mbs == 0)
  {
    C = 11;
    // f00cc
    res = (2. - LogC) / (mcs * mcs);
  }

  // 2) they are non-equal
  else if (mas != 0 and mbs != 0 and mas != mbs and mcs == 0 and mds == 0)
  {
    C = 12;
    // fab00
    res = 1. / (mas * mbs) +
          (mbs * LogA - mas * LogB) / (mas * mbs * (mas - mbs));
  }
  else if (mas != 0 and mcs != 0 and mas != mcs and mbs == 0 and mds == 0)
  {
    C = 13;
    // fa0c0
    res = 1. / (mas * mcs) +
          (mcs * LogA - mas * LogC) / (mas * mcs * (mas - mcs));
  }
  else if (mas != 0 and mds != 0 and mas != mds and mbs == 0 and mcs == 0)
  {
    C = 14;
    // fa00d
    res = 1. / (mas * mds) +
          (mds * LogA - mas * LogD) / (mas * mds * (mas - mds));
  }
  else if (mbs != 0 and mcs != 0 and mbs != mcs and mas == 0 and mds == 0)
  {
    C = 15;
    // f0bc0
    res = 1. / (mbs * mcs) +
          (mcs * LogB - mbs * LogC) / (mbs * mcs * (mbs - mcs));
  }
  else if (mbs != 0 and mds != 0 and mbs != mds and mas == 0 and mcs == 0)
  {
    C = 16;
    // f0b0d
    res = 1. / (mbs * mds) +
          (mds * LogB - mbs * LogD) / (mbs * mds * (mbs - mds));
  }
  else if (mcs != 0 and mds != 0 and mcs != mds and mas == 0 and mbs == 0)
  {
    C = 17;
    // f00cd
    res = 1. / (mcs * mds) +
          (mds * LogC - mcs * LogD) / (mcs * mds * (mcs - mds));
  }

  // three masses are non-zero, the remaining mass is zero
  // 1) the three masses are equal
  else if (mas == 0 and mbs == mcs and mbs == mds and mbs != 0)
  {
    C = 18;
    // f0bbb
    res = -1. / (2. * mbs * mbs);
  }
  else if (mbs == 0 and mas == mcs and mas == mds and mas != 0)
  {
    C = 19;
    // fa0aa
    res = -1. / (2. * mas * mas);
  }
  else if (mcs == 0 and mas == mbs and mas == mds and mas != 0)
  {
    C = 20;
    // faa0a
    res = -1. / (2. * mas * mas);
  }
  else if (mds == 0 and mas == mbs and mas == mcs and mas != 0)
  {
    C = 21;
    // faaa0
    res = -1. / (2. * mas * mas);
  }

  // 2) two of the three non-zero masses are equal
  else if (mas == 0 and mbs != 0 and mcs != 0 and mds != 0)
  {
    if (mbs == mcs and mds != mbs)
    {
      C = 22;
      // f0bbd
      res = (mbs - mds + mbs * LogD - mbs * LogB) /
            (mbs * (mbs - mds) * (mbs - mds));
    }
    else if (mbs == mds and mcs != mbs)
    {
      C = 23;
      // f0bcb
      res = (mbs - mcs + mbs * LogC - mbs * LogB) /
            (mbs * (mbs - mcs) * (mbs - mcs));
    }
    else if (mcs == mds and mbs != mcs)
    {
      C = 24;
      // f0bcc
      res = (mcs - mbs + mcs * LogB - mcs * LogC) /
            (mcs * (mbs - mcs) * (mbs - mcs));
    }
  }
  else if (mas != 0 and mbs == 0 and mcs != 0 and mds != 0)
  {
    if (mas == mcs and mds != mas)
    {
      C = 25;
      // fa0ad
      res = (mas - mds + mas * LogD - mas * LogA) /
            (mas * (mas - mds) * (mas - mds));
    }
    else if (mas == mds and mcs != mas)
    {
      C = 26;
      // fa0ca
      res = (mas - mcs + mas * LogC - mas * LogA) /
            (mas * (mas - mcs) * (mas - mcs));
    }
    else if (mcs == mds and mas != mcs)
    {
      C = 27;
      // fa0cc
      res = (mcs - mas + mcs * LogA - mcs * LogC) /
            (mcs * (mcs - mas) * (mcs - mas));
    }
  }
  else if (mas != 0 and mbs != 0 and mcs == 0 and mds != 0)
  {
    if (mas == mbs and mds != mas)
    {
      C = 28;
      // faa0d
      res = (mas - mds + mas * LogD - mas * LogA) /
            (mas * (mas - mds) * (mas - mds));
    }
    else if (mas == mds and mbs != mas)
    {
      C = 29;
      // fab0a
      res = (mas - mbs + mas * LogB - mas * LogA) /
            (mas * (mas - mbs) * (mas - mbs));
    }
    else if (mbs == mds and mas != mbs)
    {
      C = 30;
      // fab0b
      res = (mbs - mas + mbs * LogA - mbs * LogB) /
            (mbs * (mas - mbs) * (mas - mbs));
    }
  }
  else if (mas != 0 and mbs != 0 and mcs != 0 and mds == 0)
  {
    if (mas == mbs and mcs != mas)
    {
      C = 31;
      // faac0
      res = (mas - mcs + mas * LogC - mas * LogA) /
            (mas * (mas - mcs) * (mas - mcs));
    }
    else if (mas == mcs and mbs != mas)
    {
      C = 32;
      // faba0
      res = (mas - mbs + mas * LogB - mas * LogA) /
            (mas * (mas - mbs) * (mas - mbs));
    }
    else if (mbs == mcs and mas != mbs)
    {
      C = 33;
      // fabb0
      res = (mbs - mas + mbs * LogA - mbs * LogB) /
            (mbs * (mas - mbs) * (mas - mbs));
    }
  }

  // 3) all three non-zero masses are different
  else if (mas == 0 and mbs != mcs and mbs != mds and mcs != mds and
           mbs != 0 and mcs != 0 and mds != 0)
  {
    C = 34;
    // f0bcd
    res = LogB / (mbs - mcs) / (mbs - mds) + LogC / (mcs - mbs) / (mcs - mds) +
          LogD / (mds - mbs) / (mds - mcs);
  }
  else if (mbs == 0 and mas != mcs and mas != mds and mcs != mds and
           mas != 0 and mcs != 0 and mds != 0)
  {
    C = 35;
    // fa0cd
    res = LogA / (mas - mcs) / (mas - mds) + LogC / (mcs - mas) / (mcs - mds) +
          LogD / (mds - mas) / (mds - mcs);
  }
  else if (mcs == 0 and mas != mbs and mas != mds and mbs != mds and
           mas != 0 and mbs != 0 and mds != 0)
  {
    C = 36;
    // fab0d
    res = LogA / (mas - mbs) / (mas - mds) + LogB / (mbs - mas) / (mbs - mds) +
          LogD / (mds - mas) / (mds - mbs);
  }
  else if (mds == 0 and mas != mbs and mas != mcs and mbs != mcs and
           mas != 0 and mbs != 0 and mcs != 0)
  {
    C = 37;
    // fabc0
    res = LogA / (mas - mbs) / (mas - mcs) + LogB / (mbs - mas) / (mbs - mcs) +
          LogC / (mcs - mas) / (mcs - mbs);
  }

  // all four masses are non-zero
  // 1) all four masses are equal
  else if (mas == mbs and mbs == mcs and mcs == mds and mas != 0)
  {
    C = 38;
    // faaaa
    res = -1 / (6 * mas * mas);
  }
  // 2) only two of the masses are equal
  // 2.1) remaining two are not equal
  else if (mas == mbs and mcs != mds and mas != mcs and mas != mds and
           mas != 0 and mcs != 0 and mds != 0)
  {
    C = 39;
    // faacd
    res = (mcs * (mas - mds) * (mas - mds) * LogC -
           mds * (mas - mcs) * (mas - mcs) * LogD +
           (mcs - mds) *
               ((mas - mcs) * (mas - mds) + (-mas * mas + mcs * mds) * LogA)) /
          ((mas - mcs) * (mas - mcs) * (mas - mds) * (mas - mds) * (mcs - mds));
  }
  else if (mas == mcs and mbs != mds and mas != mbs and mas != mds and
           mas != 0 and mbs != 0 and mds != 0)
  {
    C = 40;
    // fabad
    res = (mbs * (mas - mds) * (mas - mds) * LogB -
           mds * (mas - mbs) * (mas - mbs) * LogD +
           (mbs - mds) *
               ((mas - mbs) * (mas - mds) + (-mas * mas + mbs * mds) * LogA)) /
          ((mas - mbs) * (mas - mbs) * (mas - mds) * (mas - mds) * (mbs - mds));
  }
  else if (mas == mds and mbs != mcs and mas != mbs and mas != mcs and
           mas != 0 and mbs != 0 and mcs != 0)
  {
    C = 41;
    // fabca
    res = (mbs * (mas - mcs) * (mas - mcs) * LogB -
           mcs * (mas - mbs) * (mas - mbs) * LogC +
           (mbs - mcs) *
               ((mas - mbs) * (mas - mcs) + (-mas * mas + mbs * mcs) * LogA)) /
          ((mas - mbs) * (mas - mbs) * (mas - mcs) * (mas - mcs) * (mbs - mcs));
  }
  else if (mbs == mcs and mas != mds and mbs != mas and mbs != mds and
           mbs != 0 and mas != 0 and mds != 0)
  {
    C = 42;
    // fabbd
    res = (mas * (mbs - mds) * (mbs - mds) * LogA -
           mds * (mas - mbs) * (mas - mbs) * LogD +
           (mas - mds) *
               (-(mas - mbs) * (mbs - mds) + (-mbs * mbs + mas * mds) * LogB)) /
          ((mas - mbs) * (mas - mbs) * (mbs - mds) * (mbs - mds) * (mas - mds));
  }
  else if (mbs == mds and mas != mcs and mbs != mas and mbs != mcs and
           mbs != 0 and mas != 0 and mcs != 0)
  {
    C = 43;
    // fabcb
    res = (mas * (mbs - mcs) * (mbs - mcs) * LogA -
           mcs * (mas - mbs) * (mas - mbs) * LogC +
           (mas - mcs) *
               (-(mas - mbs) * (mbs - mcs) + (-mbs * mbs + mas * mcs) * LogB)) /
          ((mas - mbs) * (mas - mbs) * (mbs - mcs) * (mbs - mcs) * (mas - mcs));
  }
  else if (mcs == mds and mas != mbs and mcs != mas and mcs != mbs and
           mcs != 0 and mcs != 0 and mds != 0)
  {
    C = 44;
    // fabcc
    res = (mas * (mbs - mcs) * (mbs - mcs) * LogA -
           mbs * (mas - mcs) * (mas - mcs) * LogB +
           (mas - mbs) *
               ((mas - mcs) * (mbs - mcs) + (-mcs * mcs + mas * mbs) * LogC)) /
          ((mas - mcs) * (mas - mcs) * (mbs - mcs) * (mbs - mcs) * (mas - mbs));
  }
  // 2.2) remaining two are also equal
  else if (mas == mbs and mcs == mds and mas != mcs and mas != 0 and mcs != 0)
  {
    C = 45;
    // faacc
    res = (2. * (mas - mcs) + (mas + mcs) * (LogC - LogA)) /
          ((mas - mcs) * (mas - mcs) * (mas - mcs));
  }
  else if (mas == mcs and mbs == mds and mas != mbs and mas != 0 and mbs != 0)
  {
    C = 46;
    // fabab
    res = (2. * (mas - mbs) + (mas + mbs) * (LogB - LogA)) /
          ((mas - mbs) * (mas - mbs) * (mas - mbs));
  }
  else if (mas == mds and mbs == mcs and mas != mbs and mas != 0 and mbs != 0)
  {
    C = 47;
    // fabba
    res = (2. * (mas - mbs) + (mas + mbs) * (LogB - LogA)) /
          ((mas - mbs) * (mas - mbs) * (mas - mbs));
  }

  // 3) three of the masses are equal
  else if (mas == mbs and mas == mcs and mas != 0 and mas != mds and mds != 0)
  {
    C = 48;
    // faaad
    res = (-mas * mas + mds * mds + 2. * mas * mds * (LogA - LogD)) /
          (2. * mas * (mas - mds) * (mas - mds) * (mas - mds));
  }
  else if (mas == mbs and mas == mds and mas != 0 and mas != mcs and mcs != 0)
  {
    C = 49;
    // faaca
    res = (-mas * mas + mcs * mcs + 2. * mas * mcs * (LogA - LogC)) /
          (2. * mas * (mas - mcs) * (mas - mcs) * (mas - mcs));
  }
  else if (mas == mcs and mas == mds and mas != 0 and mas != mbs and mbs != 0)
  {
    C = 50;
    // fabaa
    res = (-mas * mas + mbs * mbs + 2. * mas * mbs * (LogA - LogB)) /
          (2. * mas * (mas - mbs) * (mas - mbs) * (mas - mbs));
  }
  else if (mbs == mcs and mbs == mds and mbs != 0 and mas != mbs and mas != 0)
  {
    C = 51;
    // fabbb
    res = (-mbs * mbs + mas * mas + 2. * mbs * mas * (LogB - LogA)) /
          (2. * mbs * (mbs - mas) * (mbs - mas) * (mbs - mas));
  }

  // 4) all four masses are non-equal
  else
  {
    C = 52;
    // fabcd
    res = mas * LogA / ((mas - mbs) * (mas - mcs) * (mas - mds));
    res += mbs * LogB / ((mbs - mas) * (mbs - mcs) * (mbs - mds));
    res += mcs * LogC / ((mcs - mas) * (mcs - mbs) * (mcs - mds));
    res += mds * LogD / ((mds - mas) * (mds - mbs) * (mds - mcs));
  }

  if (std::isnan(res) or std::isinf(res))
  {
    std::string throwstring = "Found nan at line = ";
    throwstring += std::to_string(InputLineNumber);
    throwstring += " in function ";
    throwstring += __func__;
    throwstring += "\n";
    std::stringstream ss;
    // maximize the numerical precision
    typedef std::numeric_limits<double> dbl;
    ss << std::setprecision(dbl::max_digits10);
    ss << "Found nan at line = " << InputLineNumber << " in function "
       << __func__ << std::endl;
    ss << mas << sep << mbs << sep << mcs << sep << mds << sep << C << sep
       << res << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str(), __FILE__, __LINE__);
    throw std::runtime_error(throwstring.c_str());
  }

  return res;
}

double Class_Potential_Origin::fbase(double MassSquaredA,
                                     double MassSquaredB) const
{
  double res  = 0;
  double LogA = 0;
  if (MassSquaredA == 0 and MassSquaredB == 0) return 1;
  double ZB = std::pow(10, -5);
  if (MassSquaredA != 0) LogA = std::log(MassSquaredA) - 2 * std::log(scale);
  if (std::abs(MassSquaredA - MassSquaredB) > ZB)
  {
    double LogB = 0;
    if (MassSquaredB != 0) LogB = std::log(MassSquaredB) - 2 * std::log(scale);
    if (MassSquaredA == 0)
      res = LogB;
    else if (MassSquaredB == 0)
      res = LogA;
    else
      res = (LogA * MassSquaredA - LogB * MassSquaredB) /
            (MassSquaredA - MassSquaredB);
  }
  else
  {
    res = 1 + LogA;
  }
  return res;
}

std::vector<double>
Class_Potential_Origin::SecondDerivativeOfEigenvaluesNonRepeated(
    const Eigen::Ref<Eigen::MatrixXd> M,
    const Eigen::Ref<Eigen::MatrixXd> MDiffX,
    const Eigen::Ref<Eigen::MatrixXd> MDiffY,
    const Eigen::Ref<Eigen::MatrixXd> MDiffXY) const
{
  std::vector<double> res;
  const std::size_t nRows = M.rows();
  const std::size_t nCols = M.cols();

  const std::size_t EVThres = std::pow(10, -6);

  if (nCols != nRows)
  {
    throw std::runtime_error("ERROR ! M needs to be an quadratic Matrix for "
                             "calculating the derivatives !\n");
  }

  const std::size_t nSize = nRows;

  SelfAdjointEigenSolver<MatrixXd> es;
  es.compute(M);

  std::vector<double> Eigenvalues(nSize);
  for (std::size_t i = 0; i < nSize; i++)
    Eigenvalues[i] = es.eigenvalues()[i];
  for (std::size_t i = 0; i < nSize - 1; i++)
  {
    if (std::abs(Eigenvalues[i] - Eigenvalues[i + 1]) < EVThres)
    {
      Logger::Write(LoggingLevel::Default, "ERROR ! repeated eigenvalues.");
    }
  }

  std::vector<std::vector<double>> Deriv(nSize, std::vector<double>(4));
  VectorXd v(nSize);
  MatrixXd C(nSize, nSize), E(nSize, nSize), Identity(nSize, nSize);
  Identity = MatrixXd::Identity(nSize, nSize);
  VectorXd vDiffX(nSize), vDiffY(nSize);

  for (std::size_t i = 0; i < nSize; i++)
  {
    Deriv[i][0] = Eigenvalues[i];
    v           = es.eigenvectors().col(i);
    Deriv[i][1] = v.transpose() * MDiffX * v;
    Deriv[i][2] = v.transpose() * MDiffY * v;

    C = (M - Deriv[i][0] * Identity).transpose() *
            (M - Deriv[i][0] * Identity) +
        v * v.transpose();
    E = (M - Deriv[i][0] * Identity).transpose() *
        (MDiffX - Deriv[i][1] * Identity);

    vDiffX = C.colPivHouseholderQr().solve(-E * v);

    E = (M - Deriv[i][0] * Identity).transpose() *
        (MDiffY - Deriv[i][2] * Identity);
    vDiffY = C.colPivHouseholderQr().solve(-E * v);

    Deriv[i][3] = v.transpose() * MDiffXY * v;
    Deriv[i][3] += v.transpose() * (MDiffX - Deriv[i][1] * Identity) * vDiffY;
    Deriv[i][3] += v.transpose() * (MDiffY - Deriv[i][2] * Identity) * vDiffX;
  }

  for (const auto &x : Eigenvalues)
    res.push_back(x);

  for (std::size_t i = 0; i < nSize; i++)
  {
    for (std::size_t j = 0; j < 4; j++)
      res.push_back(Deriv[i][j]);
  }
  return res;
}

// Sanity check to make sure HiggsRotationMatrix is a proper rotation
// matrix, i.e. its inverse should correspond to its transpose, and its
// determinant should be +1 or -1
bool Class_Potential_Origin::CheckRotationMatrix()
{
  MatrixXd mat(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      mat(i, j) = HiggsRotationMatrix[i][j];
    }
  }

  double precision = 1e-10;

  bool AbsDetIsOne   = almost_the_same(std::abs(mat.determinant()), 1.,
                                       precision);
  bool InvEqTrans = true;

  auto inv    = mat.inverse();
  auto transp = mat.transpose();

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      if (!almost_the_same(inv(i, j), transp(i, j), precision))
      {
        InvEqTrans = false;
        break;
      }
    }
  }

  if (AbsDetIsOne and InvEqTrans)
  {
    return true;
  }
  return false;
}

void Class_Potential_Origin::CalculatePhysicalCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  const double ZeroMass = std::pow(10, -5);
  MatrixXd MassHiggs(NHiggs, NHiggs), MassGauge(NGauge, NGauge);
  MatrixXcd MassQuark(NQuarks, NQuarks), MassLepton(NLepton, NLepton);
  MassHiggs  = MatrixXd::Zero(NHiggs, NHiggs);
  MassGauge  = MatrixXd::Zero(NGauge, NGauge);
  MassQuark  = MatrixXcd::Zero(NQuarks, NQuarks);
  MassLepton = MatrixXcd::Zero(NLepton, NLepton);

  MassSquaredGauge.resize(NGauge);
  MassSquaredHiggs.resize(NHiggs);
  MassSquaredQuark.resize(NQuarks);
  MassSquaredLepton.resize(NLepton);
  HiggsRotationMatrix.resize(NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsRotationMatrix[i].resize(NHiggs);

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      MassHiggs(i, j) += Curvature_Higgs_L2[i][j];
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        MassHiggs(i, j) += Curvature_Higgs_L3[i][j][k] * HiggsVev[k];
        for (std::size_t l = 0; l < NHiggs; l++)
          MassHiggs(i, j) +=
              0.5 * Curvature_Higgs_L4[i][j][k][l] * HiggsVev[k] * HiggsVev[l];
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
        {
          MassGauge(a, b) += 0.5 * Curvature_Gauge_G2H2[a][b][i][j] *
                             HiggsVev[i] * HiggsVev[j];
        }
      }
    }
  }

  MatrixXcd MIJQuarks = QuarkMassMatrix(HiggsVev);

  MassQuark = MIJQuarks.conjugate() * MIJQuarks;

  MatrixXcd MIJLeptons = LeptonMassMatrix(HiggsVev);

  MassLepton = MIJLeptons.conjugate() * MIJLeptons;

  MatrixXd HiggsRot(NHiggs, NHiggs), GaugeRot(NGauge, NGauge),
      QuarkRot(NQuarks, NQuarks), LepRot(NLepton, NLepton);
  HiggsRot = MatrixXd::Identity(NHiggs, NHiggs);
  GaugeRot = MatrixXd::Identity(NGauge, NGauge);
  QuarkRot = MatrixXd::Identity(NQuarks, NQuarks);
  LepRot   = MatrixXd::Identity(NLepton, NLepton);

  SelfAdjointEigenSolver<MatrixXd> es;

  es.compute(MassHiggs);
  HiggsRot = es.eigenvectors().transpose();
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      if (std::abs(HiggsRot(i, j)) < std::pow(10, -10)) HiggsRot(i, j) = 0;
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    MassSquaredHiggs[i] = es.eigenvalues()[i];
    if (std::abs(MassSquaredHiggs[i]) < ZeroMass) MassSquaredHiggs[i] = 0;
  }

  es.compute(MassGauge);
  GaugeRot = es.eigenvectors().transpose();

  for (std::size_t i = 0; i < NGauge; i++)
  {
    MassSquaredGauge[i] = es.eigenvalues()[i];
    if (std::abs(MassSquaredGauge[i]) < ZeroMass) MassSquaredGauge[i] = 0;
  }

  SelfAdjointEigenSolver<MatrixXcd> esQuark(MassQuark);
  QuarkRot = esQuark.eigenvectors().transpose().real();
  for (std::size_t i = 0; i < NQuarks; i++)
    MassSquaredQuark[i] = esQuark.eigenvalues().real()[i];

  SelfAdjointEigenSolver<MatrixXcd> esLepton(MassLepton);
  LepRot = esLepton.eigenvectors().transpose().real();
  for (std::size_t i = 0; i < NLepton; i++)
    MassSquaredLepton[i] = esLepton.eigenvalues().real()[i];

  for (std::size_t a = 0; a < NGauge; a++)
  {
    for (std::size_t b = 0; b < NGauge; b++)
    {
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        LambdaGauge_3[a][b][i] = 0;
        for (std::size_t j = 0; j < NHiggs; j++)
          LambdaGauge_3[a][b][i] +=
              Curvature_Gauge_G2H2[a][b][i][j] * HiggsVev[j];
      }
    }
  }
  for (std::size_t a = 0; a < NHiggs; a++)
  {
    for (std::size_t b = 0; b < NHiggs; b++)
    {
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        LambdaHiggs_3[a][b][i] = Curvature_Higgs_L3[a][b][i];

        for (std::size_t j = 0; j < NHiggs; j++)
        {
          LambdaHiggs_3[a][b][i] +=
              Curvature_Higgs_L4[a][b][i][j] * HiggsVev[j];
        }
      }
    }
  }

  for (std::size_t i = 0; i < NQuarks; i++)
  {
    for (std::size_t j = 0; j < NQuarks; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        LambdaQuark_3[i][j][k] = 0;
        for (std::size_t l = 0; l < NQuarks; l++)
        {
          LambdaQuark_3[i][j][k] +=
              conj(Curvature_Quark_F2H1[i][l][k]) * MIJQuarks(l, j);
          LambdaQuark_3[i][j][k] +=
              conj(MIJQuarks(i, l)) * Curvature_Quark_F2H1[l][j][k];
        }
        for (std::size_t m = 0; m < NHiggs; m++)
        {
          LambdaQuark_4[i][j][k][m] = 0;
          for (std::size_t l = 0; l < NQuarks; l++)
          {
            LambdaQuark_4[i][j][k][m] += conj(Curvature_Quark_F2H1[i][l][k]) *
                                         Curvature_Quark_F2H1[l][j][m];
            LambdaQuark_4[i][j][k][m] += conj(Curvature_Quark_F2H1[i][l][m]) *
                                         Curvature_Quark_F2H1[l][j][k];
          }
        }
      }
    }
  }

  for (std::size_t i = 0; i < NLepton; i++)
  {
    for (std::size_t j = 0; j < NLepton; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        LambdaLepton_3[i][j][k] = 0;
        for (std::size_t l = 0; l < NLepton; l++)
        {
          LambdaLepton_3[i][j][k] +=
              conj(Curvature_Lepton_F2H1[i][l][k]) * MIJLeptons(l, j);
          LambdaLepton_3[i][j][k] +=
              conj(MIJLeptons(i, l)) * Curvature_Lepton_F2H1[l][j][k];
        }
        for (std::size_t m = 0; m < NHiggs; m++)
        {
          LambdaLepton_4[i][j][k][m] = 0;
          for (std::size_t l = 0; l < NLepton; l++)
          {
            LambdaLepton_4[i][j][k][m] += conj(Curvature_Lepton_F2H1[i][l][k]) *
                                          Curvature_Lepton_F2H1[l][j][m];
            LambdaLepton_4[i][j][k][m] += conj(Curvature_Lepton_F2H1[i][l][m]) *
                                          Curvature_Lepton_F2H1[l][j][k];
          }
        }
      }
    }
  }

  // Rotate and save std::size_to corresponding vectors

  Couplings_Higgs_Quartic.resize(NHiggs);
  Couplings_Higgs_Triple.resize(NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    Couplings_Higgs_Quartic[i].resize(NHiggs);
    Couplings_Higgs_Triple[i].resize(NHiggs);
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      Couplings_Higgs_Quartic[i][j].resize(NHiggs);
      Couplings_Higgs_Triple[i][j].resize(NHiggs);
      for (std::size_t k = 0; k < NHiggs; k++)
        Couplings_Higgs_Quartic[i][j][k].resize(NHiggs);
    }
  }

  Couplings_Gauge_Higgs_22.resize(NGauge);
  Couplings_Gauge_Higgs_21.resize(NGauge);
  for (std::size_t a = 0; a < NGauge; a++)
  {
    Couplings_Gauge_Higgs_22[a].resize(NGauge);
    Couplings_Gauge_Higgs_21[a].resize(NGauge);
    for (std::size_t b = 0; b < NGauge; b++)
    {
      Couplings_Gauge_Higgs_22[a][b].resize(NHiggs);
      Couplings_Gauge_Higgs_21[a][b].resize(NHiggs);
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        Couplings_Gauge_Higgs_22[a][b][i].resize(NHiggs);
      }
    }
  }

  Couplings_Quark_Higgs_22.resize(NQuarks);
  Couplings_Quark_Higgs_21.resize(NQuarks);
  for (std::size_t a = 0; a < NQuarks; a++)
  {
    Couplings_Quark_Higgs_22[a].resize(NQuarks);
    Couplings_Quark_Higgs_21[a].resize(NQuarks);
    for (std::size_t b = 0; b < NQuarks; b++)
    {
      Couplings_Quark_Higgs_22[a][b].resize(NHiggs);
      Couplings_Quark_Higgs_21[a][b].resize(NHiggs);
      for (std::size_t i = 0; i < NHiggs; i++)
        Couplings_Quark_Higgs_22[a][b][i].resize(NHiggs);
    }
  }

  Couplings_Lepton_Higgs_22.resize(NLepton);
  Couplings_Lepton_Higgs_21.resize(NLepton);
  for (std::size_t a = 0; a < NLepton; a++)
  {
    Couplings_Lepton_Higgs_22[a].resize(NLepton);
    Couplings_Lepton_Higgs_21[a].resize(NLepton);
    for (std::size_t b = 0; b < NLepton; b++)
    {
      Couplings_Lepton_Higgs_22[a][b].resize(NHiggs);
      Couplings_Lepton_Higgs_21[a][b].resize(NHiggs);
      for (std::size_t i = 0; i < NHiggs; i++)
        Couplings_Lepton_Higgs_22[a][b][i].resize(NHiggs);
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        Couplings_Higgs_Triple[i][j][k] = 0;
        for (std::size_t is = 0; is < NHiggs; is++)
        {
          for (std::size_t js = 0; js < NHiggs; js++)
          {
            for (std::size_t ks = 0; ks < NHiggs; ks++)
            {
              Couplings_Higgs_Triple[i][j][k] +=
                  HiggsRot(i, is) * HiggsRot(j, js) * HiggsRot(k, ks) *
                  LambdaHiggs_3[is][js][ks];
            }
          }
        }
        for (std::size_t l = 0; l < NHiggs; l++)
        {
          Couplings_Higgs_Quartic[i][j][k][l] = 0;
          for (std::size_t is = 0; is < NHiggs; is++)
          {
            for (std::size_t js = 0; js < NHiggs; js++)
            {
              for (std::size_t ks = 0; ks < NHiggs; ks++)
              {
                for (std::size_t ls = 0; ls < NHiggs; ls++)
                {
                  Couplings_Higgs_Quartic[i][j][k][l] +=
                      HiggsRot(i, is) * HiggsRot(j, js) * HiggsRot(k, ks) *
                      HiggsRot(l, ls) * Curvature_Higgs_L4[is][js][ks][ls];
                }
              }
            }
          }
        }
      }
    }
  }

  // Gauge Rot
  for (std::size_t a = 0; a < NGauge; a++)
  {
    for (std::size_t b = 0; b < NGauge; b++)
    {
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        Couplings_Gauge_Higgs_21[a][b][i] = 0;
        for (std::size_t as = 0; as < NGauge; as++)
        {
          for (std::size_t bs = 0; bs < NGauge; bs++)
          {
            for (std::size_t is = 0; is < NHiggs; is++)
              Couplings_Gauge_Higgs_21[a][b][i] +=
                  GaugeRot(a, as) * GaugeRot(b, bs) * HiggsRot(i, is) *
                  LambdaGauge_3[as][bs][is];
          }
        }
        for (std::size_t j = 0; j < NHiggs; j++)
        {
          Couplings_Gauge_Higgs_22[a][b][i][j] = 0;
          for (std::size_t as = 0; as < NGauge; as++)
          {
            for (std::size_t bs = 0; bs < NGauge; bs++)
            {
              for (std::size_t is = 0; is < NHiggs; is++)
              {
                for (std::size_t js = 0; js < NHiggs; js++)
                {
                  double RotFac = GaugeRot(a, as) * GaugeRot(b, bs) *
                                  HiggsRot(i, is) * HiggsRot(j, js);
                  Couplings_Gauge_Higgs_22[a][b][i][j] +=
                      RotFac * Curvature_Gauge_G2H2[as][bs][is][js];
                }
              }
            }
          }
        }
      }
    }
  }

  // Quark

  for (std::size_t a = 0; a < NQuarks; a++)
  {
    for (std::size_t b = 0; b < NQuarks; b++)
    {
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        Couplings_Quark_Higgs_21[a][b][i] = 0;
        for (std::size_t as = 0; as < NQuarks; as++)
        {
          for (std::size_t bs = 0; bs < NQuarks; bs++)
          {
            for (std::size_t is = 0; is < NHiggs; is++)
            {
              double RotFac =
                  QuarkRot(a, as) * QuarkRot(b, bs) * HiggsRot(i, is);
              Couplings_Quark_Higgs_21[a][b][i] +=
                  RotFac * LambdaQuark_3[as][bs][is];
            }
          }
        }
        for (std::size_t j = 0; j < NHiggs; j++)
        {
          Couplings_Quark_Higgs_22[a][b][i][j] = 0;
          for (std::size_t as = 0; as < NQuarks; as++)
          {
            for (std::size_t bs = 0; bs < NQuarks; bs++)
            {
              for (std::size_t is = 0; is < NHiggs; is++)
              {
                for (std::size_t js = 0; js < NHiggs; js++)
                {
                  double RotFac = QuarkRot(a, as) * QuarkRot(b, bs) *
                                  HiggsRot(i, is) * HiggsRot(j, js);
                  Couplings_Quark_Higgs_22[a][b][i][j] +=
                      RotFac * LambdaQuark_4[as][bs][is][js];
                }
              }
            }
          }
        }
      }
    }
  }

  // Lepton

  for (std::size_t a = 0; a < NLepton; a++)
  {
    for (std::size_t b = 0; b < NLepton; b++)
    {
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        Couplings_Lepton_Higgs_21[a][b][i] = 0;
        for (std::size_t as = 0; as < NLepton; as++)
        {
          for (std::size_t bs = 0; bs < NLepton; bs++)
          {
            for (std::size_t is = 0; is < NHiggs; is++)
            {
              double RotFac = LepRot(a, as) * LepRot(b, bs) * HiggsRot(i, is);
              Couplings_Lepton_Higgs_21[a][b][i] +=
                  RotFac * LambdaLepton_3[as][bs][is];
            }
          }
        }
        for (std::size_t j = 0; j < NHiggs; j++)
        {
          Couplings_Lepton_Higgs_22[a][b][i][j] = 0;
          for (std::size_t as = 0; as < NLepton; as++)
          {
            for (std::size_t bs = 0; bs < NLepton; bs++)
            {
              for (std::size_t is = 0; is < NHiggs; is++)
              {
                for (std::size_t js = 0; js < NHiggs; js++)
                {
                  double RotFac = LepRot(a, as) * LepRot(b, bs) *
                                  HiggsRot(i, is) * HiggsRot(j, js);
                  Couplings_Lepton_Higgs_22[a][b][i][j] +=
                      RotFac * LambdaLepton_4[as][bs][is][js];
                }
              }
            }
          }
        }
      }
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRotationMatrix[i][j] = HiggsRot(i, j);
    }
  }

  CalcCouplingsDone = true;

  return;
}

std::vector<double> Class_Potential_Origin::WeinbergFirstDerivative() const
{
  std::vector<double> res;
  if (!CalcCouplingsDone)
  {
    //        CalculatePhysicalCouplings();
    std::string retmes = __func__;
    retmes += " tries to use Physical couplings but they are not initialised.";
    throw std::runtime_error(retmes);
  }
  const double NumZero = std::pow(10, -10);
  VectorXd FirstDeriv(NHiggs), FirstDerivGauge(NHiggs), FirstDerivHiggs(NHiggs),
      FirstDerivQuark(NHiggs), FirstDerivLepton(NHiggs);
  FirstDeriv       = VectorXd::Zero(NHiggs);
  FirstDerivGauge  = VectorXd::Zero(NHiggs);
  FirstDerivHiggs  = VectorXd::Zero(NHiggs);
  FirstDerivQuark  = VectorXd::Zero(NHiggs);
  FirstDerivLepton = VectorXd::Zero(NHiggs);
  double epsilon   = 1.0 / (16 * M_PI * M_PI);

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t a = 0; a < NGauge; a++)
    {
      if (MassSquaredGauge[a] != 0)
      {
        FirstDerivGauge[i] +=
            MassSquaredGauge[a] * Couplings_Gauge_Higgs_21[a][a][i] *
            (std::log(MassSquaredGauge[a] / std::pow(scale, 2)) - C_CWcbGB +
             0.5);
      }
    }

    for (std::size_t a = 0; a < NHiggs; a++)
    {
      if (MassSquaredHiggs[a] != 0)
      {

        FirstDerivHiggs[i] +=
            MassSquaredHiggs[a] * Couplings_Higgs_Triple[a][a][i] *
            (std::log(MassSquaredHiggs[a] / std::pow(scale, 2)) - C_CWcbHiggs +
             0.5);
      }
    }
    for (std::size_t a = 0; a < NQuarks; a++)
    {
      if (MassSquaredQuark[a] != 0)
      {
        double Coup = Couplings_Quark_Higgs_21[a][a][i].real();
        FirstDerivQuark[i] +=
            MassSquaredQuark[a] * Coup *
            (std::log(MassSquaredQuark[a] / std::pow(scale, 2)) -
             C_CWcbFermion + 0.5);
      }
    }
    for (std::size_t a = 0; a < NLepton; a++)
    {
      if (MassSquaredLepton[a] != 0)
      {
        double Coup = Couplings_Lepton_Higgs_21[a][a][i].real();
        FirstDerivLepton[i] +=
            MassSquaredLepton[a] * Coup *
            (std::log(MassSquaredLepton[a] / std::pow(scale, 2)) -
             C_CWcbFermion + 0.5);
      }
    }
  }
  FirstDerivGauge *= 1.5;
  FirstDerivHiggs *= 0.5;
  FirstDerivQuark *= -3;
  FirstDerivLepton *= -1;

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
  }

  FirstDeriv = HiggsRot.transpose() * (FirstDerivGauge + FirstDerivHiggs +
                                       FirstDerivQuark + FirstDerivLepton);
  FirstDeriv *= epsilon;

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    if (std::abs(FirstDeriv[i]) < NumZero) FirstDeriv[i] = 0;
  }

  for (std::size_t i = 0; i < NHiggs; i++)
    res.push_back(FirstDeriv[i]);

  return res;
}

Eigen::MatrixXd
Class_Potential_Origin::WeinbergSecondDerivativeAsMatrixXd() const
{
  if (!CalcCouplingsDone)
  {
    //        CalculatePhysicalCouplings();
    std::string retmes = __func__;
    retmes += " tries to use Physical couplings but they are not initialised.";
    throw std::runtime_error(retmes);
  }
  const double NumZero = std::pow(10, -10);
  MatrixXd GaugePart(NHiggs, NHiggs), HiggsPart(NHiggs, NHiggs),
      QuarkPart(NHiggs, NHiggs), LeptonPart(NHiggs, NHiggs);
  GaugePart  = MatrixXd::Zero(NHiggs, NHiggs);
  HiggsPart  = MatrixXd::Zero(NHiggs, NHiggs);
  QuarkPart  = MatrixXd::Zero(NHiggs, NHiggs);
  LeptonPart = MatrixXd::Zero(NHiggs, NHiggs);

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t a = 0; a < NGauge; a++)
      {
        for (std::size_t b = 0; b < NGauge; b++)
        {
          double Coup1 = Couplings_Gauge_Higgs_21[a][b][i];
          double Coup2 = Couplings_Gauge_Higgs_21[b][a][j];
          double Br =
              fbase(MassSquaredGauge[a], MassSquaredGauge[b]) - C_CWcbGB + 0.5;
          GaugePart(i, j) += Coup1 * Coup2 * Br;
        }
        if (MassSquaredGauge[a] != 0)
        {
          GaugePart(i, j) +=
              MassSquaredGauge[a] * Couplings_Gauge_Higgs_22[a][a][i][j] *
              (std::log(MassSquaredGauge[a] / std::pow(scale, 2)) - C_CWcbGB +
               0.5);
        }
      }

      for (std::size_t a = 0; a < NHiggs; a++)
      {
        for (std::size_t b = 0; b < NHiggs; b++)
        {
          double Coup1 = Couplings_Higgs_Triple[a][b][i];
          double Coup2 = Couplings_Higgs_Triple[b][a][j];
          double Br    = fbase(MassSquaredHiggs[a], MassSquaredHiggs[b]) -
                      C_CWcbHiggs + 0.5;
          HiggsPart(i, j) += Coup1 * Coup2 * Br;
        }
        if (MassSquaredHiggs[a] != 0)
        {
          HiggsPart(i, j) +=
              MassSquaredHiggs[a] * Couplings_Higgs_Quartic[a][a][i][j] *
              (std::log(MassSquaredHiggs[a] / std::pow(scale, 2)) -
               C_CWcbHiggs + 0.5);
        }
      }

      for (std::size_t a = 0; a < NQuarks; a++)
      {
        for (std::size_t b = 0; b < NQuarks; b++)
        {
          double Coup = (Couplings_Quark_Higgs_21[a][b][i] *
                         Couplings_Quark_Higgs_21[b][a][j])
                            .real();
          double Br = fbase(MassSquaredQuark[a], MassSquaredQuark[b]) -
                      C_CWcbFermion + 0.5;
          QuarkPart(i, j) += Coup * Br;
        }
        if (MassSquaredQuark[a] != 0)
        {
          double Coup = Couplings_Quark_Higgs_22[a][a][i][j].real();
          QuarkPart(i, j) +=
              MassSquaredQuark[a] * Coup *
              (std::log(MassSquaredQuark[a] / std::pow(scale, 2)) -
               C_CWcbFermion + 0.5);
        }
      }

      for (std::size_t a = 0; a < NLepton; a++)
      {
        for (std::size_t b = 0; b < NLepton; b++)
        {
          double Coup = (Couplings_Lepton_Higgs_21[a][b][i] *
                         Couplings_Lepton_Higgs_21[b][a][j])
                            .real();
          double Br = fbase(MassSquaredLepton[a], MassSquaredLepton[b]) -
                      C_CWcbFermion + 0.5;
          LeptonPart(i, j) += Coup * Br;
        }
        if (MassSquaredLepton[a] != 0)
        {
          double Coup = Couplings_Lepton_Higgs_22[a][a][i][j].real();
          LeptonPart(i, j) +=
              Coup * MassSquaredLepton[a] *
              (std::log(MassSquaredLepton[a] / std::pow(scale, 2)) -
               C_CWcbFermion + 0.5);
        }
      }
    }
  }

  HiggsPart *= 0.5;
  GaugePart *= 1.5;
  QuarkPart *= -3;
  LeptonPart *= -1;

  MatrixXd Storage(NHiggs, NHiggs);
  Storage = HiggsPart + GaugePart + QuarkPart + LeptonPart;

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      if (std::abs(Storage(i, j)) < NumZero) Storage(i, j) = 0;
    }
  }
  MatrixXd ResMatrix;
  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
    }
  }

  ResMatrix =
      0.5 * HiggsRot.transpose() * (Storage + Storage.transpose()) * HiggsRot;
  double epsilon = 1.0 / (16.0 * M_PI * M_PI);
  ResMatrix *= epsilon;

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      if (std::abs(ResMatrix(i, j)) < NumZero) ResMatrix(i, j) = 0;
    }
  }
  return ResMatrix;
}
std::vector<double> Class_Potential_Origin::WeinbergSecondDerivative() const
{

  auto ResMatrix = WeinbergSecondDerivativeAsMatrixXd();
  std::vector<double> res;
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      res.push_back(ResMatrix(j, i));
    }
  }

  return res;
}

std::vector<double> Class_Potential_Origin::WeinbergThirdDerivative() const
{

  if (not CalcCouplingsDone)
  {
    std::string retmes = __func__;
    retmes += " tries to use Physical couplings but they are not initialised.";
    throw std::runtime_error(retmes);
  }
  const double NumZero = std::pow(10, -10);
  double epsilon       = 1.0 / (16.0 * M_PI * M_PI);

  std::vector<double> res;

  std::vector<std::vector<std::vector<std::complex<double>>>> restmp;
  std::vector<std::vector<std::vector<std::complex<double>>>> QuarkPart;
  std::vector<std::vector<std::vector<std::complex<double>>>> LeptonPart;
  std::vector<std::vector<std::vector<std::complex<double>>>> QuarkPartSym;
  std::vector<std::vector<std::vector<std::complex<double>>>> LeptonPartSym;
  restmp.resize(NHiggs);
  QuarkPart.resize(NHiggs);
  LeptonPart.resize(NHiggs);
  QuarkPartSym.resize(NHiggs);
  LeptonPartSym.resize(NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    restmp[i].resize(NHiggs);
    QuarkPart[i].resize(NHiggs);
    LeptonPart[i].resize(NHiggs);
    QuarkPartSym[i].resize(NHiggs);
    LeptonPartSym[i].resize(NHiggs);
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      restmp[i][j].resize(NHiggs);
      QuarkPart[i][j].resize(NHiggs);
      LeptonPart[i][j].resize(NHiggs);
      QuarkPartSym[i][j].resize(NHiggs);
      LeptonPartSym[i][j].resize(NHiggs);
    }
  }

  std::vector<std::vector<std::vector<double>>> resGaugeBase(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));
  std::vector<std::vector<std::vector<double>>> Higgspart(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));

  std::vector<std::vector<std::vector<double>>> GaugePart(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));

  std::vector<std::vector<std::vector<double>>> HiggspartSym(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));

  std::vector<std::vector<std::vector<double>>> GaugePartSym(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        Higgspart[i][j][k] = 0;
        for (std::size_t a = 0; a < NHiggs; a++)
        {
          for (std::size_t b = 0; b < NHiggs; b++)
          {
            for (std::size_t c = 0; c < NHiggs; c++)
            {
              double f1 = fbaseTri(MassSquaredHiggs[a],
                                   MassSquaredHiggs[b],
                                   MassSquaredHiggs[c]);
              double f2 = Couplings_Higgs_Triple[a][b][i];
              double f3 = Couplings_Higgs_Triple[b][c][j];
              double f4 = Couplings_Higgs_Triple[c][a][k];
              Higgspart[i][j][k] += 2 * f1 * f2 * f3 * f4;
            }
            double f1 = Couplings_Higgs_Quartic[a][b][i][j];
            double f2 = Couplings_Higgs_Triple[b][a][k];
            double f3 = fbase(MassSquaredHiggs[a], MassSquaredHiggs[b]) -
                        C_CWcbHiggs + 0.5;
            Higgspart[i][j][k] += 3 * f1 * f2 * f3;
          }
        }

        GaugePart[i][j][k] = 0;
        for (std::size_t a = 0; a < NGauge; a++)
        {
          for (std::size_t b = 0; b < NGauge; b++)
          {
            for (std::size_t c = 0; c < NGauge; c++)
            {
              double f1 = fbaseTri(MassSquaredGauge[a],
                                   MassSquaredGauge[b],
                                   MassSquaredGauge[c]);
              double f2 = Couplings_Gauge_Higgs_21[a][b][i];
              double f3 = Couplings_Gauge_Higgs_21[b][c][j];
              double f4 = Couplings_Gauge_Higgs_21[c][a][k];
              GaugePart[i][j][k] += 2 * f1 * f2 * f3 * f4;
            }
            double f1 = Couplings_Gauge_Higgs_22[a][b][i][j];
            double f2 = Couplings_Gauge_Higgs_21[b][a][k];
            double f3 = fbase(MassSquaredGauge[a], MassSquaredGauge[b]) -
                        C_CWcbGB + 0.5;
            GaugePart[i][j][k] += 3 * f1 * f2 * f3;
          }
        }

        QuarkPart[i][j][k] = 0;
        for (std::size_t a = 0; a < NQuarks; a++)
        {
          for (std::size_t b = 0; b < NQuarks; b++)
          {
            for (std::size_t c = 0; c < NQuarks; c++)
            {
              std::complex<double> f1 = fbaseTri(MassSquaredQuark[a],
                                                 MassSquaredQuark[b],
                                                 MassSquaredQuark[c]);
              std::complex<double> f2 = Couplings_Quark_Higgs_21[a][b][i];
              std::complex<double> f3 = Couplings_Quark_Higgs_21[b][c][j];
              std::complex<double> f4 = Couplings_Quark_Higgs_21[c][a][k];
              QuarkPart[i][j][k] += 2.0 * f1 * f2 * f3 * f4;
            }
            std::complex<double> f1 = Couplings_Quark_Higgs_22[a][b][i][j];
            std::complex<double> f2 = Couplings_Quark_Higgs_21[b][a][k];
            std::complex<double> f3 =
                fbase(MassSquaredQuark[a], MassSquaredQuark[b]) -
                C_CWcbFermion + 0.5;
            QuarkPart[i][j][k] += 3.0 * f1 * f2 * f3;
          }
        }
        LeptonPart[i][j][k] = 0;
        for (std::size_t a = 0; a < NLepton; a++)
        {
          for (std::size_t b = 0; b < NLepton; b++)
          {
            for (std::size_t c = 0; c < NLepton; c++)
            {
              std::complex<double> f1 = fbaseTri(MassSquaredLepton[a],
                                                 MassSquaredLepton[b],
                                                 MassSquaredLepton[c]);
              std::complex<double> f2 = Couplings_Lepton_Higgs_21[a][b][i];
              std::complex<double> f3 = Couplings_Lepton_Higgs_21[b][c][j];
              std::complex<double> f4 = Couplings_Lepton_Higgs_21[c][a][k];
              LeptonPart[i][j][k] += 2.0 * f1 * f2 * f3 * f4;
            }
            std::complex<double> f1 = Couplings_Lepton_Higgs_22[a][b][i][j];
            std::complex<double> f2 = Couplings_Lepton_Higgs_21[b][a][k];
            std::complex<double> f3 =
                fbase(MassSquaredLepton[a], MassSquaredLepton[b]) -
                C_CWcbFermion + 0.5;
            LeptonPart[i][j][k] += 3.0 * f1 * f2 * f3;
          }
        }
      }
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        HiggspartSym[i][j][k] = Higgspart[i][j][k] + Higgspart[i][k][j];
        HiggspartSym[i][j][k] += Higgspart[j][i][k] + Higgspart[j][k][i];
        HiggspartSym[i][j][k] += Higgspart[k][i][j] + Higgspart[k][j][i];
        HiggspartSym[i][j][k] *= 1.0 / 6.0;

        GaugePartSym[i][j][k] = GaugePart[i][j][k] + GaugePart[i][k][j];
        GaugePartSym[i][j][k] += GaugePart[j][i][k] + GaugePart[j][k][i];
        GaugePartSym[i][j][k] += GaugePart[k][i][j] + GaugePart[k][j][i];
        GaugePartSym[i][j][k] *= 1.0 / 6.0;

        QuarkPartSym[i][j][k] = QuarkPart[i][j][k] + QuarkPart[i][k][j];
        QuarkPartSym[i][j][k] += QuarkPart[j][i][k] + QuarkPart[j][k][i];
        QuarkPartSym[i][j][k] += QuarkPart[k][i][j] + QuarkPart[k][j][i];
        QuarkPartSym[i][j][k] *= 1.0 / 6.0;

        LeptonPartSym[i][j][k] = LeptonPart[i][j][k] + LeptonPart[i][k][j];
        LeptonPartSym[i][j][k] += LeptonPart[j][i][k] + LeptonPart[j][k][i];
        LeptonPartSym[i][j][k] += LeptonPart[k][i][j] + LeptonPart[k][j][i];
        LeptonPartSym[i][j][k] *= 1.0 / 6.0;
      }
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        restmp[i][j][k] = 0.5 * HiggspartSym[i][j][k];
        restmp[i][j][k] += 1.5 * GaugePartSym[i][j][k];
        restmp[i][j][k] += -1.0 * LeptonPartSym[i][j][k];
        restmp[i][j][k] += -3.0 * QuarkPartSym[i][j][k];
      }
    }
  }

  for (std::size_t l = 0; l < NHiggs; l++)
  {
    for (std::size_t m = 0; m < NHiggs; m++)
    {
      for (std::size_t n = 0; n < NHiggs; n++)
      {
        resGaugeBase[l][m][n] = 0;
        for (std::size_t i = 0; i < NHiggs; i++)
        {
          for (std::size_t j = 0; j < NHiggs; j++)
          {
            for (std::size_t k = 0; k < NHiggs; k++)
            {
              double RotFac = HiggsRotationMatrix[i][l] *
                              HiggsRotationMatrix[j][m] *
                              HiggsRotationMatrix[k][n];
              resGaugeBase[l][m][n] += RotFac * restmp[i][j][k].real();
            }
          }
        }
        resGaugeBase[l][m][n] *= epsilon;
        if (std::abs(resGaugeBase[l][m][n]) < NumZero)
          resGaugeBase[l][m][n] = 0;
      }
    }
  }

  for (std::size_t l = 0; l < NHiggs; l++)
  {
    for (std::size_t m = 0; m < NHiggs; m++)
    {
      for (std::size_t n = 0; n < NHiggs; n++)
      {
        res.push_back(resGaugeBase[l][m][n]);
      }
    }
  }

  return res;
}

std::vector<double> Class_Potential_Origin::WeinbergForthDerivative() const
{

  if (not CalcCouplingsDone)
  {
    std::string retmes = __func__;
    retmes += " tries to use Physical couplings but they are not initialised.";
    throw std::runtime_error(retmes);
  }

  const double NumZero = std::pow(10, -10);
  double epsilon       = 1.0 / (16.0 * M_PI * M_PI);

  std::vector<double> res;

  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>
      restmp;
  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>
      QuarkPart;
  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>
      LeptonPart;
  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>
      QuarkPartSym;
  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>
      LeptonPartSym;
  std::vector<std::vector<std::vector<std::vector<double>>>> resGaugeBase;
  std::vector<std::vector<std::vector<std::vector<double>>>> HiggsPart;
  std::vector<std::vector<std::vector<std::vector<double>>>> GaugePart;
  std::vector<std::vector<std::vector<std::vector<double>>>> HiggsPartSym;
  std::vector<std::vector<std::vector<std::vector<double>>>> GaugePartSym;

  restmp.resize(NHiggs);
  QuarkPart.resize(NHiggs);
  LeptonPart.resize(NHiggs);
  QuarkPartSym.resize(NHiggs);
  LeptonPartSym.resize(NHiggs);
  resGaugeBase.resize(NHiggs);
  HiggsPart.resize(NHiggs);
  GaugePart.resize(NHiggs);
  HiggsPartSym.resize(NHiggs);
  GaugePartSym.resize(NHiggs);

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    restmp[i].resize(NHiggs);
    QuarkPart[i].resize(NHiggs);
    LeptonPart[i].resize(NHiggs);
    QuarkPartSym[i].resize(NHiggs);
    LeptonPartSym[i].resize(NHiggs);
    resGaugeBase[i].resize(NHiggs);
    HiggsPart[i].resize(NHiggs);
    GaugePart[i].resize(NHiggs);
    HiggsPartSym[i].resize(NHiggs);
    GaugePartSym[i].resize(NHiggs);
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      restmp[i][j].resize(NHiggs);
      QuarkPart[i][j].resize(NHiggs);
      LeptonPart[i][j].resize(NHiggs);
      QuarkPartSym[i][j].resize(NHiggs);
      LeptonPartSym[i][j].resize(NHiggs);
      resGaugeBase[i][j].resize(NHiggs);
      HiggsPart[i][j].resize(NHiggs);
      GaugePart[i][j].resize(NHiggs);
      HiggsPartSym[i][j].resize(NHiggs);
      GaugePartSym[i][j].resize(NHiggs);
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        restmp[i][j][k].resize(NHiggs);
        QuarkPart[i][j][k].resize(NHiggs);
        LeptonPart[i][j][k].resize(NHiggs);
        QuarkPartSym[i][j][k].resize(NHiggs);
        LeptonPartSym[i][j][k].resize(NHiggs);
        resGaugeBase[i][j][k].resize(NHiggs);
        HiggsPart[i][j][k].resize(NHiggs);
        GaugePart[i][j][k].resize(NHiggs);
        HiggsPartSym[i][j][k].resize(NHiggs);
        GaugePartSym[i][j][k].resize(NHiggs);
      }
    }
  }

  for (std::size_t i1 = 0; i1 < NHiggs; i1++)
  {
    for (std::size_t i2 = 0; i2 < NHiggs; i2++)
    {
      for (std::size_t i3 = 0; i3 < NHiggs; i3++)
      {
        for (std::size_t i4 = 0; i4 < NHiggs; i4++)
        {
          HiggsPart[i1][i2][i3][i4] = 0;
          for (std::size_t a = 0; a < NHiggs; a++)
          {
            for (std::size_t b = 0; b < NHiggs; b++)
            {
              double f1 = Couplings_Higgs_Quartic[a][b][i1][i4];
              double f2 = Couplings_Higgs_Quartic[b][a][i2][i3];
              double f3 = fbase(MassSquaredHiggs[a], MassSquaredHiggs[b]) -
                          C_CWcbHiggs + 0.5;

              HiggsPart[i1][i2][i3][i4] += f1 * f2 * f3;

              for (std::size_t c = 0; c < NHiggs; c++)
              {
                for (std::size_t d = 0; d < NHiggs; d++)
                {
                  double f11 = fbaseFour(MassSquaredHiggs[a],
                                         MassSquaredHiggs[b],
                                         MassSquaredHiggs[c],
                                         MassSquaredHiggs[d]);
                  double f12 = Couplings_Higgs_Triple[a][b][i4];
                  double f13 = Couplings_Higgs_Triple[b][c][i1];
                  double f14 = Couplings_Higgs_Triple[c][d][i2];
                  double f15 = Couplings_Higgs_Triple[d][a][i3];

                  HiggsPart[i1][i2][i3][i4] +=
                      2.0 * f11 * f12 * f13 * f14 * f15;
                }

                double f21 = fbaseTri(MassSquaredHiggs[a],
                                      MassSquaredHiggs[b],
                                      MassSquaredHiggs[c]);
                double f22 = Couplings_Higgs_Triple[a][b][i4];
                double f23 = Couplings_Higgs_Quartic[b][c][i1][i2];
                double f24 = Couplings_Higgs_Triple[c][a][i3];

                HiggsPart[i1][i2][i3][i4] += 4.0 * f21 * f22 * f23 * f24;
              }
            }
          }

          GaugePart[i1][i2][i3][i4] = 0;
          for (std::size_t a = 0; a < NGauge; a++)
          {
            for (std::size_t b = 0; b < NGauge; b++)
            {
              double f1 = Couplings_Gauge_Higgs_22[a][b][i1][i4];
              double f2 = Couplings_Gauge_Higgs_22[b][a][i2][i3];
              double f3 = fbase(MassSquaredGauge[a], MassSquaredGauge[b]) -
                          C_CWcbGB + 0.5;

              GaugePart[i1][i2][i3][i4] += f1 * f2 * f3;

              for (std::size_t c = 0; c < NGauge; c++)
              {
                for (std::size_t d = 0; d < NGauge; d++)
                {
                  double f11 = fbaseFour(MassSquaredGauge[a],
                                         MassSquaredGauge[b],
                                         MassSquaredGauge[c],
                                         MassSquaredGauge[d]);
                  double f12 = Couplings_Gauge_Higgs_21[a][b][i4];
                  double f13 = Couplings_Gauge_Higgs_21[b][c][i1];
                  double f14 = Couplings_Gauge_Higgs_21[c][d][i2];
                  double f15 = Couplings_Gauge_Higgs_21[d][a][i3];

                  GaugePart[i1][i2][i3][i4] +=
                      2.0 * f11 * f12 * f13 * f14 * f15;
                }

                double f21 = fbaseTri(MassSquaredGauge[a],
                                      MassSquaredGauge[b],
                                      MassSquaredGauge[c]);
                double f22 = Couplings_Gauge_Higgs_21[a][b][i4];
                double f23 = Couplings_Gauge_Higgs_22[b][c][i1][i2];
                double f24 = Couplings_Gauge_Higgs_21[c][a][i3];

                GaugePart[i1][i2][i3][i4] += 4.0 * f21 * f22 * f23 * f24;
              }
            }
          }

          QuarkPart[i1][i2][i3][i4] = 0;
          for (std::size_t a = 0; a < NQuarks; a++)
          {
            for (std::size_t b = 0; b < NQuarks; b++)
            {
              std::complex<double> f1 = Couplings_Quark_Higgs_22[a][b][i1][i4];
              std::complex<double> f2 = Couplings_Quark_Higgs_22[b][a][i2][i3];
              std::complex<double> f3 =
                  fbase(MassSquaredQuark[a], MassSquaredQuark[b]) -
                  C_CWcbFermion + 0.5;

              QuarkPart[i1][i2][i3][i4] += f1 * f2 * f3;

              for (std::size_t c = 0; c < NQuarks; c++)
              {
                for (std::size_t d = 0; d < NQuarks; d++)
                {
                  std::complex<double> f11 = fbaseFour(MassSquaredQuark[a],
                                                       MassSquaredQuark[b],
                                                       MassSquaredQuark[c],
                                                       MassSquaredQuark[d]);
                  std::complex<double> f12 = Couplings_Quark_Higgs_21[a][b][i4];
                  std::complex<double> f13 = Couplings_Quark_Higgs_21[b][c][i1];
                  std::complex<double> f14 = Couplings_Quark_Higgs_21[c][d][i2];
                  std::complex<double> f15 = Couplings_Quark_Higgs_21[d][a][i3];

                  QuarkPart[i1][i2][i3][i4] +=
                      2.0 * f11 * f12 * f13 * f14 * f15;
                }

                std::complex<double> f21 = fbaseTri(MassSquaredQuark[a],
                                                    MassSquaredQuark[b],
                                                    MassSquaredQuark[c]);
                std::complex<double> f22 = Couplings_Quark_Higgs_21[a][b][i4];
                std::complex<double> f23 =
                    Couplings_Quark_Higgs_22[b][c][i1][i2];
                std::complex<double> f24 = Couplings_Quark_Higgs_21[c][a][i3];

                QuarkPart[i1][i2][i3][i4] += 4.0 * f21 * f22 * f23 * f24;
              }
            }
          }

          LeptonPart[i1][i2][i3][i4] = 0;
          for (std::size_t a = 0; a < NLepton; a++)
          {
            for (std::size_t b = 0; b < NLepton; b++)
            {
              std::complex<double> f1 = Couplings_Lepton_Higgs_22[a][b][i1][i4];
              std::complex<double> f2 = Couplings_Lepton_Higgs_22[b][a][i2][i3];
              std::complex<double> f3 =
                  fbase(MassSquaredLepton[a], MassSquaredLepton[b]) -
                  C_CWcbFermion + 0.5;

              LeptonPart[i1][i2][i3][i4] += f1 * f2 * f3;

              for (std::size_t c = 0; c < NLepton; c++)
              {
                for (std::size_t d = 0; d < NLepton; d++)
                {
                  std::complex<double> f11 = fbaseFour(MassSquaredLepton[a],
                                                       MassSquaredLepton[b],
                                                       MassSquaredLepton[c],
                                                       MassSquaredLepton[d]);
                  std::complex<double> f12 =
                      Couplings_Lepton_Higgs_21[a][b][i4];
                  std::complex<double> f13 =
                      Couplings_Lepton_Higgs_21[b][c][i1];
                  std::complex<double> f14 =
                      Couplings_Lepton_Higgs_21[c][d][i2];
                  std::complex<double> f15 =
                      Couplings_Lepton_Higgs_21[d][a][i3];

                  LeptonPart[i1][i2][i3][i4] +=
                      2.0 * f11 * f12 * f13 * f14 * f15;
                }

                std::complex<double> f21 = fbaseTri(MassSquaredLepton[a],
                                                    MassSquaredLepton[b],
                                                    MassSquaredLepton[c]);
                std::complex<double> f22 = Couplings_Lepton_Higgs_21[a][b][i4];
                std::complex<double> f23 =
                    Couplings_Lepton_Higgs_22[b][c][i1][i2];
                std::complex<double> f24 = Couplings_Lepton_Higgs_21[c][a][i3];

                LeptonPart[i1][i2][i3][i4] += 4.0 * f21 * f22 * f23 * f24;
              }
            }
          }
        }
      }
    }
  }

  for (std::size_t i1 = 0; i1 < NHiggs; i1++)
  {
    for (std::size_t i2 = 0; i2 < NHiggs; i2++)
    {
      for (std::size_t i3 = 0; i3 < NHiggs; i3++)
      {
        for (std::size_t i4 = 0; i4 < NHiggs; i4++)
        {
          HiggsPartSym[i1][i2][i3][i4] =
              HiggsPart[i1][i2][i3][i4] + HiggsPart[i1][i2][i4][i3] +
              HiggsPart[i1][i3][i2][i4] + HiggsPart[i1][i3][i4][i2] +
              HiggsPart[i1][i4][i2][i3] + HiggsPart[i1][i4][i3][i2];
          HiggsPartSym[i1][i2][i3][i4] +=
              HiggsPart[i2][i1][i3][i4] + HiggsPart[i2][i1][i4][i3] +
              HiggsPart[i2][i3][i1][i4] + HiggsPart[i2][i3][i4][i1] +
              HiggsPart[i2][i4][i1][i3] + HiggsPart[i2][i4][i3][i1];
          HiggsPartSym[i1][i2][i3][i4] +=
              HiggsPart[i3][i1][i2][i4] + HiggsPart[i3][i1][i4][i2] +
              HiggsPart[i3][i2][i1][i4] + HiggsPart[i3][i2][i4][i1] +
              HiggsPart[i3][i4][i1][i2] + HiggsPart[i3][i4][i2][i1];
          HiggsPartSym[i1][i2][i3][i4] +=
              HiggsPart[i4][i1][i2][i3] + HiggsPart[i4][i1][i3][i2] +
              HiggsPart[i4][i2][i1][i3] + HiggsPart[i4][i2][i3][i1] +
              HiggsPart[i4][i3][i1][i2] + HiggsPart[i4][i3][i2][i1];
          HiggsPartSym[i1][i2][i3][i4] *= 1.0 / 24.0;

          GaugePartSym[i1][i2][i3][i4] =
              GaugePart[i1][i2][i3][i4] + GaugePart[i1][i2][i4][i3] +
              GaugePart[i1][i3][i2][i4] + GaugePart[i1][i3][i4][i2] +
              GaugePart[i1][i4][i2][i3] + GaugePart[i1][i4][i3][i2];
          GaugePartSym[i1][i2][i3][i4] +=
              GaugePart[i2][i1][i3][i4] + GaugePart[i2][i1][i4][i3] +
              GaugePart[i2][i3][i1][i4] + GaugePart[i2][i3][i4][i1] +
              GaugePart[i2][i4][i1][i3] + GaugePart[i2][i4][i3][i1];
          GaugePartSym[i1][i2][i3][i4] +=
              GaugePart[i3][i1][i2][i4] + GaugePart[i3][i1][i4][i2] +
              GaugePart[i3][i2][i1][i4] + GaugePart[i3][i2][i4][i1] +
              GaugePart[i3][i4][i1][i2] + GaugePart[i3][i4][i2][i1];
          GaugePartSym[i1][i2][i3][i4] +=
              GaugePart[i4][i1][i2][i3] + GaugePart[i4][i1][i3][i2] +
              GaugePart[i4][i2][i1][i3] + GaugePart[i4][i2][i3][i1] +
              GaugePart[i4][i3][i1][i2] + GaugePart[i4][i3][i2][i1];
          GaugePartSym[i1][i2][i3][i4] *= 1.0 / 24.0;

          QuarkPartSym[i1][i2][i3][i4] =
              QuarkPart[i1][i2][i3][i4] + QuarkPart[i1][i2][i4][i3] +
              QuarkPart[i1][i3][i2][i4] + QuarkPart[i1][i3][i4][i2] +
              QuarkPart[i1][i4][i2][i3] + QuarkPart[i1][i4][i3][i2];
          QuarkPartSym[i1][i2][i3][i4] +=
              QuarkPart[i2][i1][i3][i4] + QuarkPart[i2][i1][i4][i3] +
              QuarkPart[i2][i3][i1][i4] + QuarkPart[i2][i3][i4][i1] +
              QuarkPart[i2][i4][i1][i3] + QuarkPart[i2][i4][i3][i1];
          QuarkPartSym[i1][i2][i3][i4] +=
              QuarkPart[i3][i1][i2][i4] + QuarkPart[i3][i1][i4][i2] +
              QuarkPart[i3][i2][i1][i4] + QuarkPart[i3][i2][i4][i1] +
              QuarkPart[i3][i4][i1][i2] + QuarkPart[i3][i4][i2][i1];
          QuarkPartSym[i1][i2][i3][i4] +=
              QuarkPart[i4][i1][i2][i3] + QuarkPart[i4][i1][i3][i2] +
              QuarkPart[i4][i2][i1][i3] + QuarkPart[i4][i2][i3][i1] +
              QuarkPart[i4][i3][i1][i2] + QuarkPart[i4][i3][i2][i1];
          QuarkPartSym[i1][i2][i3][i4] *= 1.0 / 24.0;

          LeptonPartSym[i1][i2][i3][i4] =
              LeptonPart[i1][i2][i3][i4] + LeptonPart[i1][i2][i4][i3] +
              LeptonPart[i1][i3][i2][i4] + LeptonPart[i1][i3][i4][i2] +
              LeptonPart[i1][i4][i2][i3] + LeptonPart[i1][i4][i3][i2];
          LeptonPartSym[i1][i2][i3][i4] +=
              LeptonPart[i2][i1][i3][i4] + LeptonPart[i2][i1][i4][i3] +
              LeptonPart[i2][i3][i1][i4] + LeptonPart[i2][i3][i4][i1] +
              LeptonPart[i2][i4][i1][i3] + LeptonPart[i2][i4][i3][i1];
          LeptonPartSym[i1][i2][i3][i4] +=
              LeptonPart[i3][i1][i2][i4] + LeptonPart[i3][i1][i4][i2] +
              LeptonPart[i3][i2][i1][i4] + LeptonPart[i3][i2][i4][i1] +
              LeptonPart[i3][i4][i1][i2] + LeptonPart[i3][i4][i2][i1];
          LeptonPartSym[i1][i2][i3][i4] +=
              LeptonPart[i4][i1][i2][i3] + LeptonPart[i4][i1][i3][i2] +
              LeptonPart[i4][i2][i1][i3] + LeptonPart[i4][i2][i3][i1] +
              LeptonPart[i4][i3][i1][i2] + LeptonPart[i4][i3][i2][i1];
          LeptonPartSym[i1][i2][i3][i4] *= 1.0 / 24.0;
        }
      }
    }
  }

  for (std::size_t i1 = 0; i1 < NHiggs; i1++)
  {
    for (std::size_t i2 = 0; i2 < NHiggs; i2++)
    {
      for (std::size_t i3 = 0; i3 < NHiggs; i3++)
      {
        for (std::size_t i4 = 0; i4 < NHiggs; i4++)
        {
          restmp[i1][i2][i3][i4] = 3.0 * 0.5 * HiggsPartSym[i1][i2][i3][i4];
          restmp[i1][i2][i3][i4] += 3.0 * 1.5 * GaugePartSym[i1][i2][i3][i4];
          restmp[i1][i2][i3][i4] +=
              3.0 * (-1.0) * LeptonPartSym[i1][i2][i3][i4];
          restmp[i1][i2][i3][i4] += 3.0 * (-3.0) * QuarkPartSym[i1][i2][i3][i4];
        }
      }
    }
  }

  for (std::size_t j1 = 0; j1 < NHiggs; j1++)
  {
    for (std::size_t j2 = 0; j2 < NHiggs; j2++)
    {
      for (std::size_t j3 = 0; j3 < NHiggs; j3++)
      {
        for (std::size_t j4 = 0; j4 < NHiggs; j4++)
        {
          resGaugeBase[j1][j2][j3][j4] = 0;

          for (std::size_t i1 = 0; i1 < NHiggs; i1++)
          {
            for (std::size_t i2 = 0; i2 < NHiggs; i2++)
            {
              for (std::size_t i3 = 0; i3 < NHiggs; i3++)
              {
                for (std::size_t i4 = 0; i4 < NHiggs; i4++)
                {
                  double RotFac = HiggsRotationMatrix[i1][j1] *
                                  HiggsRotationMatrix[i2][j2] *
                                  HiggsRotationMatrix[i3][j3] *
                                  HiggsRotationMatrix[i4][j4];
                  resGaugeBase[j1][j2][j3][j4] +=
                      RotFac * restmp[i1][i2][i3][i4].real();
                }
              }
            }
          }

          resGaugeBase[j1][j2][j3][j4] *= epsilon;
          if (std::abs(resGaugeBase[j1][j2][j3][j4]) < NumZero)
            resGaugeBase[j1][j2][j3][j4] = 0;
        }
      }
    }
  }

  for (std::size_t j1 = 0; j1 < NHiggs; j1++)
  {
    for (std::size_t j2 = 0; j2 < NHiggs; j2++)
    {
      for (std::size_t j3 = 0; j3 < NHiggs; j3++)
      {
        for (std::size_t j4 = 0; j4 < NHiggs; j4++)
        {
          res.push_back(resGaugeBase[j4][j3][j2][j1]);
        }
      }
    }
  }

  return res;
}

MatrixXd Class_Potential_Origin::HiggsMassMatrix(const std::vector<double> &v,
                                                 double Temp,
                                                 int diff) const
{
  MatrixXd res(NHiggs, NHiggs);
  if (v.size() != nVEV and v.size() != NHiggs)
  {
    std::string ErrorString =
        std::string("You have called ") + std::string(__func__) +
        std::string(
            " with an invalid vev configuration. Your vev is of dimension ") +
        std::to_string(v.size()) + std::string(" and it should be ") +
        std::to_string(NHiggs) + std::string(".");
    throw std::runtime_error(ErrorString);
  }
  if (v.size() == nVEV and nVEV != NHiggs)
  {
    std::stringstream ss;
    ss << __func__
       << " is being called with a wrong sized vev configuration. It "
          "has the dimension of "
       << nVEV << " while it should have " << NHiggs
       << ". For now this is transformed but please fix this to reduce "
          "the runtime."
       << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
    std::vector<double> Transformedv;
    Transformedv = MinimizeOrderVEV(v);
    res          = HiggsMassMatrix(Transformedv, Temp, diff);
    return res;
  }
  if (!SetCurvatureDone)
  {
    //        SetCurvatureArrays();
    throw std::runtime_error(
        "SetCurvatureDone is not set. The Model is not initiliased correctly");
  }

  if (diff == 0)
  {
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      for (std::size_t j = i; j < NHiggs; j++)
      {
        res(i, j) = Curvature_Higgs_L2[i][j];
        for (std::size_t k = 0; k < NHiggs; k++)
        {
          res(i, j) += Curvature_Higgs_L3[i][j][k] * v[k];
          for (std::size_t l = 0; l < NHiggs; l++)
          {
            res(i, j) += 0.5 * Curvature_Higgs_L4[i][j][k][l] * v[k] * v[l];
          }
        }

        if (Temp != 0)
        {
          res(i, j) += DebyeHiggs[i][j] * std::pow(Temp, 2);
        }
      }
    }
    for (std::size_t i{1}; i < NHiggs; ++i)
    {
      for (std::size_t j{0}; j < i; ++j)
      {
        res(i, j) = res(j, i);
      }
    }
  }
  else if (static_cast<size_t>(diff) <= NHiggs and diff > 0)
  {
    std::size_t x0 = diff - 1;
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      for (std::size_t j = 0; j < NHiggs; j++)
      {
        res(i, j) = Curvature_Higgs_L3[i][j][x0];
        for (std::size_t k = 0; k < NHiggs; k++)
        {
          res(i, j) += Curvature_Higgs_L4[i][j][x0][k] * v[k];
        }
      }
    }
  }
  else if (diff == -1)
  {
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      for (std::size_t j = 0; j < NHiggs; j++)
      {
        res(i, j) = 2 * DebyeHiggs[i][j] * Temp;
      }
    }
  }
  return res;
}

std::vector<double>
Class_Potential_Origin::HiggsMassesSquared(const std::vector<double> &v,
                                           const double &Temp,
                                           const int &diff) const
{
  std::vector<double> res;

  auto MassMatrix = HiggsMassMatrix(v, Temp, diff);

  double ZeroMass = std::pow(10, -5);

  if (diff == 0 and res.size() == 0)
  {
    SelfAdjointEigenSolver<MatrixXd> es(MassMatrix, EigenvaluesOnly);
    const auto EV = es.eigenvalues();
    for (std::size_t i{0}; i < NHiggs; ++i)
    {
      if (std::abs(EV[i]) < ZeroMass)
      {
        res.push_back(0);
      }
      else
      {
        res.push_back(EV[i]);
      }
    }
  }
  else if (diff == 0 and res.size() == NHiggs)
  {
    SelfAdjointEigenSolver<MatrixXd> es(MassMatrix, EigenvaluesOnly);
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      double tmp = es.eigenvalues()[i];
      if (std::abs(tmp) < ZeroMass) tmp = 0;
      res[i] = tmp;
    }
  }
  else if (diff == 0 and res.size() != 0 and res.size() != NHiggs)
  {
    Logger::Write(LoggingLevel::Debug,
                  std::string("Something went wrong in ") + __func__ + ".\n" +
                      __func__ + "Is calculating the mass for " +
                      std::to_string(NHiggs) +
                      "fields but the resolution vector has a size of " +
                      std::to_string(res.size()) + ". This should be zero or " +
                      std::to_string(NHiggs));
  }
  else if (static_cast<std::size_t>(diff) <= NHiggs and diff > 0)
  {
    MatrixXd Diff(NHiggs, NHiggs);
    std::size_t x0 = diff - 1;
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      for (std::size_t j = 0; j < NHiggs; j++)
      {
        Diff(i, j) = Curvature_Higgs_L3[i][j][x0];
        for (std::size_t k = 0; k < NHiggs; k++)
        {
          Diff(i, j) += Curvature_Higgs_L4[i][j][x0][k] * v[k];
        }
      }
    }
    MatrixXcd MassCast(NHiggs, NHiggs);
    MassCast = MassMatrix;
    MatrixXcd DiffCast(NHiggs, NHiggs);
    DiffCast = Diff;
    res      = FirstDerivativeOfEigenvalues(MassCast, DiffCast);
  }
  else if (diff == -1)
  {
    MatrixXd Diff(NHiggs, NHiggs);
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      for (std::size_t j = 0; j < NHiggs; j++)
      {
        Diff(i, j) = 2 * DebyeHiggs[i][j] * Temp;
      }
    }

    MatrixXcd MassCast(NHiggs, NHiggs);
    MassCast = MassMatrix;
    MatrixXcd DiffCast(NHiggs, NHiggs);
    DiffCast = Diff;
    res      = FirstDerivativeOfEigenvalues(MassCast, DiffCast);
  }

  return res;
}

std::vector<double>
Class_Potential_Origin::GaugeMassesSquared(const std::vector<double> &v,
                                           const double &Temp,
                                           const int &diff) const
{
  std::vector<double> res;
  if (v.size() != nVEV and v.size() != NHiggs)
  {
    std::string ErrorString =
        std::string("You have called ") + std::string(__func__) +
        std::string(
            " with an invalid vev configuration. Your vev is of dimension ") +
        std::to_string(v.size()) + std::string(" and it should be ") +
        std::to_string(NHiggs) + std::string(".");
    throw std::runtime_error(ErrorString);
  }
  if (v.size() == nVEV and nVEV != NHiggs)
  {
    std::stringstream ss;
    ss << __func__
       << " is being called with a wrong sized vev configuration. It "
          "has the dimension of "
       << nVEV << " while it should have " << NHiggs
       << ". For now this is transformed but please fix this to reduce "
          "the runtime."
       << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
    std::vector<double> Transformedv;
    Transformedv = MinimizeOrderVEV(v);
    res          = GaugeMassesSquared(Transformedv, Temp, diff);
    return res;
  }
  if (!SetCurvatureDone)
  {
    //        SetCurvatureArrays();
    std::string retmes = __func__;
    retmes += "was called while the model was not initialised correctly.\n";
    throw std::runtime_error(retmes);
  }
  MatrixXd MassMatrix(NGauge, NGauge);
  double ZeroMass = std::pow(10, -5);
  for (std::size_t a = 0; a < NGauge; a++)
  {
    for (std::size_t b = 0; b < NGauge; b++)
    {
      MassMatrix(a, b) = 0;
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        for (std::size_t j = 0; j < NHiggs; j++)
          MassMatrix(a, b) +=
              0.5 * Curvature_Gauge_G2H2[a][b][i][j] * v.at(i) * v.at(j);
      }

      if (Temp != 0)
      {
        MassMatrix(a, b) += DebyeGauge[a][b] * std::pow(Temp, 2);
      }
    }
  }

  for (std::size_t a{1}; a < NGauge; ++a)
  {
    for (std::size_t b{0}; b < a; ++b)
    {
      MassMatrix(a, b) = MassMatrix(b, a);
    }
  }

  if (diff == 0)
  {

    SelfAdjointEigenSolver<MatrixXd> es(MassMatrix, EigenvaluesOnly);
    for (std::size_t i = 0; i < NGauge; i++)
    {
      double tmp = es.eigenvalues()[i];
      if (std::abs(tmp) < ZeroMass)
        res.push_back(0);
      else
        res.push_back(tmp);
    }
  }
  else if (diff > 0 and static_cast<size_t>(diff) <= NHiggs)
  {
    std::size_t i = diff - 1;
    MatrixXd Diff(NGauge, NGauge);
    Diff = MatrixXd::Zero(NGauge, NGauge);
    for (std::size_t a = 0; a < NGauge; a++)
    {
      for (std::size_t b = 0; b < NGauge; b++)
      {
        for (std::size_t j = 0; j < NHiggs; j++)
          Diff(a, b) += Curvature_Gauge_G2H2[a][b][i][j] * v[j];
      }
    }
    MatrixXcd MassCast(NGauge, NGauge);
    MassCast = MassMatrix;
    MatrixXcd DiffCast(NGauge, NGauge);
    DiffCast = Diff;
    res      = FirstDerivativeOfEigenvalues(MassCast, DiffCast);
  }
  else if (diff == -1)
  {
    MatrixXd Diff(NGauge, NGauge);
    for (std::size_t i = 0; i < NGauge; i++)
    {
      for (std::size_t j = 0; j < NGauge; j++)
      {
        Diff(i, j) = 2 * DebyeGauge[i][j] * Temp;
      }
    }

    MatrixXcd MassCast(NGauge, NGauge);
    MassCast = MassMatrix;
    MatrixXcd DiffCast(NGauge, NGauge);
    DiffCast = Diff;
    res      = FirstDerivativeOfEigenvalues(MassCast, DiffCast);
  }

  return res;
}

std::vector<double>
Class_Potential_Origin::QuarkMassesSquared(const std::vector<double> &v,
                                           const int &diff) const
{
  std::vector<double> res;
  MatrixXcd MassMatrix(NQuarks, NQuarks), MIJ(NQuarks, NQuarks);
  MIJ             = QuarkMassMatrix(v);
  double ZeroMass = std::pow(10, -10);

  MassMatrix = MIJ.conjugate() * MIJ;

  if (diff <= 0) // No temperature dependent part here
  {
    SelfAdjointEigenSolver<MatrixXcd> es(MassMatrix, EigenvaluesOnly);
    for (std::size_t i = 0; i < NQuarks; i++)
    {
      double tmp = es.eigenvalues().real()[i];
      if (std::abs(tmp) < ZeroMass)
        res.push_back(0);
      else
        res.push_back(tmp);
    }
  }
  else if (static_cast<size_t>(diff) <= NHiggs)
  {
    std::size_t m = diff - 1;
    MatrixXcd Diff(NQuarks, NQuarks);
    Diff = MatrixXcd::Zero(NQuarks, NQuarks);
    for (std::size_t a = 0; a < NQuarks; a++)
    {
      for (std::size_t b = 0; b < NQuarks; b++)
      {
        for (std::size_t i = 0; i < NQuarks; i++)
        {
          Diff(a, b) += std::conj(Curvature_Quark_F2H1[a][i][m]) * MIJ(i, b);
          Diff(a, b) += std::conj(MIJ(a, i)) * Curvature_Quark_F2H1[i][b][m];
        }
      }
    }

    res = FirstDerivativeOfEigenvalues(MassMatrix, Diff);

    for (std::size_t j = 0; j < res.size(); j++)
    {
      if (std::isnan(res.at(j)))
      {
        std::stringstream ss;
        ss << "MassMatrix = \n"
           << MassMatrix << "\nDiff = \n"
           << Diff << std::endl;
        ss << "Fermion Masses : ";
        for (std::size_t i = 0; i < NQuarks; i++)
          ss << std::sqrt(std::abs(res.at(i))) << sep;
        ss << std::endl;
        ss << "VEV fields : ";
        for (std::size_t i = 0; i < v.size(); i++)
          ss << v.at(i) << sep;
        ss << std::endl;

        for (std::size_t l = 0; l < NHiggs; l++)
        {

          ss << "Curvature_Quark * v an Higgs  =  :" << l << "\n";
          for (std::size_t a = 0; a < NQuarks; a++)
          {
            for (std::size_t i = 0; i < NQuarks; i++)
            {
              ss << Curvature_Quark_F2H1[a][i][l] * v[l] << sep;
            }
            ss << std::endl;
          }
          ss << "conj Curvature_Quark an Higgs = :" << l << "\n";
          for (std::size_t a = 0; a < NQuarks; a++)
          {
            for (std::size_t i = 0; i < NQuarks; i++)
            {
              ss << std::conj(Curvature_Quark_F2H1[a][i][l]) * v[l] << sep;
            }
            ss << std::endl;
          }
        }

        Logger::Write(LoggingLevel::Debug, ss.str());

        std::string retmessage = "Nan found in ";
        retmessage += __func__;
        retmessage += " at deriv number ";
        retmessage += std::to_string(j);
        retmessage += " and m = ";
        retmessage += std::to_string(m);
        throw std::runtime_error(retmessage);
      }
    }
  }

  return res;
}

std::vector<double>
Class_Potential_Origin::LeptonMassesSquared(const std::vector<double> &v,
                                            const int &diff) const
{
  std::vector<double> res;
  MatrixXcd MassMatrix(NLepton, NLepton), MIJ(NLepton, NLepton);
  double ZeroMass = std::pow(10, -10);
  MIJ             = LeptonMassMatrix(v);

  MassMatrix = MIJ.conjugate() * MIJ;

  if (diff <= 0) // no temperature part here
  {
    SelfAdjointEigenSolver<MatrixXcd> es(MassMatrix, EigenvaluesOnly);
    for (std::size_t i = 0; i < NLepton; i++)
    {
      double tmp = es.eigenvalues().real()[i];
      if (std::abs(tmp) < ZeroMass)
        res.push_back(0);
      else
        res.push_back(tmp);
    }
  }
  else if (static_cast<size_t>(diff) <= NHiggs)
  {

    auto k         = diff - 1;
    MatrixXcd Diff = MatrixXcd::Zero(NLepton, NLepton);
    for (std::size_t I{0}; I < NLepton; ++I)
    {
      for (std::size_t J{0}; J < NLepton; ++J)
      {
        for (std::size_t L{0}; L < NLepton; ++L)
        {
          Diff(I, J) += std::conj(Curvature_Lepton_F2H1[I][L][k]) * MIJ(L, J);
          Diff(I, J) += std::conj(MIJ(I, L)) * Curvature_Lepton_F2H1[L][J][k];
        }
      }
    }

    res = FirstDerivativeOfEigenvalues(MassMatrix, Diff);

    for (std::size_t j = 0; j < res.size(); j++)
    {
      if (std::isnan(res.at(j)))
      {
        std::stringstream ss;
        ss << "MassMatrix = \n"
           << MassMatrix << "\nDiff = \n"
           << Diff << std::endl;
        ss << "Fermion Masses : ";
        for (std::size_t i = 0; i < NLepton; i++)
          ss << std::sqrt(std::abs(res.at(i))) << sep;
        ss << std::endl;
        ss << "VEV fields : ";
        for (std::size_t i = 0; i < v.size(); i++)
          ss << v.at(i) << sep;
        ss << std::endl;

        for (std::size_t l = 0; l < NHiggs; l++)
        {

          ss << "Curvature_Lepton * v an Higgs  =  :" << l << "\n";
          for (std::size_t a = 0; a < NLepton; a++)
          {
            for (std::size_t i = 0; i < NLepton; i++)
            {
              ss << Curvature_Lepton_F2H1[a][i][l] * v[l] << sep;
            }
            ss << std::endl;
          }
          ss << "conj Curvature_Lepton an Higgs = :" << l << "\n";
          for (std::size_t a = 0; a < NLepton; a++)
          {
            for (std::size_t i = 0; i < NLepton; i++)
            {
              ss << std::conj(Curvature_Lepton_F2H1[a][i][l]) * v[l] << sep;
            }
            ss << std::endl;
          }
        }

        Logger::Write(LoggingLevel::Debug, ss.str());

        std::string retmessage = "Nan found in ";
        retmessage += __func__;
        retmessage += " at deriv number ";
        retmessage += std::to_string(j);
        retmessage += " and m = ";
        retmessage += std::to_string(diff - 1);
        throw std::runtime_error(retmessage);
      }
    }
  }

  return res;
}

double Class_Potential_Origin::VTree(const std::vector<double> &v,
                                     int diff,
                                     bool ForceExplicitCalculation) const
{
  double res = 0;

  if (not ForceExplicitCalculation)
  {
    res = VTreeSimplified(v);
    if (UseVTreeSimplified and diff == 0)
    {

      return res;
    }
  }
  res = 0;

  if (diff == 0)
  {
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      if (v[i] != 0)
      {
        res += Curvature_Higgs_L1[i] * v[i];
        for (std::size_t j = 0; j < NHiggs; j++)
        {
          if (v[j] != 0)
          {
            res += 0.5 * Curvature_Higgs_L2[i][j] * v[i] * v[j];
            for (std::size_t k = 0; k < NHiggs; k++)
            {
              res +=
                  1.0 / 6.0 * Curvature_Higgs_L3[i][j][k] * v[i] * v[j] * v[k];
              for (std::size_t l = 0; l < NHiggs; l++)
              {
                res += 1.0 / 24.0 * Curvature_Higgs_L4[i][j][k][l] * v[i] *
                       v[j] * v[k] * v[l];
              }
            }
          }
        }
      }
    }
  }
  else if (diff > 0 and static_cast<size_t>(diff) <= NHiggs)
  {
    std::size_t i = diff - 1;
    res           = Curvature_Higgs_L1[i];
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      res += Curvature_Higgs_L2[i][j] * v[j];
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        res += 0.5 * Curvature_Higgs_L3[i][j][k] * v[j] * v[k];
        for (std::size_t l = 0; l < NHiggs; l++)
        {
          res +=
              1.0 / 6.0 * Curvature_Higgs_L4[i][j][k][l] * v[j] * v[k] * v[l];
        }
      }
    }
  }

  return res;
}

double Class_Potential_Origin::CounterTerm(const std::vector<double> &v,
                                           int diff,
                                           bool ForceExplicitCalculation) const
{

  double res = 0;
  if (not ForceExplicitCalculation and UseVCounterSimplified)
  {
    res = VCounterSimplified(v);
    if (UseVCounterSimplified and diff == 0) return res;
  }

  res = 0;
  if (diff == 0)
  {
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      res += Curvature_Higgs_CT_L1[i] * v[i];
      for (std::size_t j = 0; j < NHiggs; j++)
      {
        res += 0.5 * Curvature_Higgs_CT_L2[i][j] * v[i] * v[j];
        for (std::size_t k = 0; k < NHiggs; k++)
        {
          res +=
              1.0 / 6.0 * Curvature_Higgs_CT_L3[i][j][k] * v[i] * v[j] * v[k];
          for (std::size_t l = 0; l < NHiggs; l++)
          {
            res += 1.0 / 24.0 * Curvature_Higgs_CT_L4[i][j][k][l] * v[i] *
                   v[j] * v[k] * v[l];
          }
        }
      }
    }
  }
  else if (diff > 0 and static_cast<size_t>(diff) <= NHiggs)
  {
    std::size_t i = diff - 1;
    res           = Curvature_Higgs_CT_L1[i];
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      res += Curvature_Higgs_CT_L2[i][j] * v[j];
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        res += 0.5 * Curvature_Higgs_CT_L3[i][j][k] * v[j] * v[k];
        for (std::size_t l = 0; l < NHiggs; l++)
        {
          res += 1.0 / 6.0 * Curvature_Higgs_CT_L4[i][j][k][l] * v[j] * v[k] *
                 v[l];
        }
      }
    }
  }

  return res;
}

double Class_Potential_Origin::VEff(const std::vector<double> &v,
                                    double Temp,
                                    int diff,
                                    int Order) const
{
  if (v.size() != nVEV and v.size() != NHiggs)
  {
    std::string ErrorString =
        std::string("You have called ") + std::string(__func__) +
        std::string(
            " with an invalid vev configuration. Your vev is of dimension ") +
        std::to_string(v.size()) + std::string(" and it should be ") +
        std::to_string(NHiggs) + std::string(".");
    throw std::runtime_error(ErrorString);
  }
  if (v.size() == nVEV and nVEV != NHiggs)
  {
    std::stringstream ss;
    ss << __func__
       << " is being called with a wrong sized vev configuration. It "
          "has the dimension of "
       << nVEV << " while it should have " << NHiggs
       << ". For now this is transformed but please fix this to reduce "
          "the runtime."
       << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
    std::vector<double> Transformedv;
    Transformedv = MinimizeOrderVEV(v);
    return VEff(Transformedv, Temp, diff);
  }

  double resOut = 0;
  resOut        = VTree(v, diff);
  if (Order != 0 and not UseTreeLevel)
  {
    resOut += CounterTerm(v, diff);
    resOut += V1Loop(v, Temp, diff);
  }
  // for(std::size_t i=0;i<NHiggs;i++) resOut +=
  // DebyeHiggs[i][i]*0.5*std::pow(v.at(i),2)*std::pow(Temp,2);
  return resOut;
}

double Class_Potential_Origin::V1Loop(const std::vector<double> &v,
                                      double Temp,
                                      int diff) const
{
  double res = 0;

  /**
   * In case of diff != 0 the mass vectors will directly have the derivatives of
   * the masses
   */
  std::vector<double> HiggsMassesVec, QuarkMassesVec, GaugeMassesVec,
      LeptonMassesVec, HiggsMassesZeroTempVec, GaugeMassesZeroTempVec;
  HiggsMassesVec         = HiggsMassesSquared(v, Temp, diff);
  GaugeMassesVec         = GaugeMassesSquared(v, Temp, diff);
  GaugeMassesZeroTempVec = GaugeMassesSquared(v, 0, diff);
  QuarkMassesVec         = QuarkMassesSquared(v, diff);
  LeptonMassesVec        = LeptonMassesSquared(v, diff);

  if (diff == 0)
  {
    if (C_UseParwani)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
        res += boson(HiggsMassesVec[k], Temp, C_CWcbHiggs, 0);
      for (std::size_t k = 0; k < NGauge; k++)
        res += boson(GaugeMassesVec[k], Temp, C_CWcbGB, 0);
      for (std::size_t k = 0; k < NGauge; k++)
        res += 2 * boson(GaugeMassesZeroTempVec[k], Temp, C_CWcbGB, 0);
      for (std::size_t k = 0; k < NQuarks; k++)
        res += -6 * fermion(QuarkMassesVec[k], Temp, 0);
      for (std::size_t k = 0; k < NLepton; k++)
        res += -2 * fermion(LeptonMassesVec[k], Temp, 0);
    }
    else
    {
      HiggsMassesZeroTempVec = HiggsMassesSquared(v, 0);
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        res += boson(HiggsMassesZeroTempVec[k], Temp, C_CWcbHiggs, 0);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        res += 3 * boson(GaugeMassesZeroTempVec[k], Temp, C_CWcbGB, 0);
      }
      double AddContQuark = 0;
      for (std::size_t k = 0; k < NQuarks; k++)
      {
        AddContQuark += -2 * fermion(QuarkMassesVec[k], Temp, 0);
      }
      res += NColour * AddContQuark;
      for (std::size_t k = 0; k < NLepton; k++)
      {
        res += -2 * fermion(LeptonMassesVec[k], Temp, 0);
      }

      double VDebye = 0;
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        if (HiggsMassesVec[k] > 0) VDebye += std::pow(HiggsMassesVec[k], 1.5);
        if (HiggsMassesZeroTempVec[k] > 0)
          VDebye += -std::pow(HiggsMassesZeroTempVec[k], 1.5);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        if (GaugeMassesVec[k] > 0) VDebye += std::pow(GaugeMassesVec[k], 1.5);
        if (GaugeMassesZeroTempVec[k] > 0)
          VDebye += -std::pow(GaugeMassesZeroTempVec[k], 1.5);
      }

      VDebye *= -Temp / (12 * M_PI);
      res += VDebye;
    }
  }
  else if (diff > 0 and static_cast<size_t>(diff) <= NHiggs)
  {
    if (C_UseParwani)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        res += HiggsMassesVec.at(k + NHiggs) *
               boson(HiggsMassesVec.at(k), Temp, C_CWcbHiggs, diff);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        res += GaugeMassesVec.at(k + NGauge) *
               boson(GaugeMassesVec.at(k), Temp, C_CWcbHiggs, diff);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        res += 2 * GaugeMassesZeroTempVec.at(k + NGauge) *
               boson(GaugeMassesZeroTempVec.at(k), Temp, C_CWcbHiggs, diff);
      }
      for (std::size_t k = 0; k < NQuarks; k++)
      {
        res += -6 * QuarkMassesVec.at(k + NQuarks) *
               fermion(QuarkMassesVec[k], Temp, diff);
      }
      for (std::size_t k = 0; k < NLepton; k++)
      {
        res += -2 * LeptonMassesVec.at(k + NLepton) *
               fermion(LeptonMassesVec[k], Temp, diff);
      }
    }
    else
    {
      HiggsMassesZeroTempVec = HiggsMassesSquared(v, 0, diff);
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        res += HiggsMassesZeroTempVec.at(k + NHiggs) *
               boson(HiggsMassesZeroTempVec[k], Temp, C_CWcbHiggs, diff);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        res += 3 * GaugeMassesZeroTempVec.at(k + NGauge) *
               boson(GaugeMassesZeroTempVec[k], Temp, C_CWcbGB, diff);
      }
      double AddContQuark = 0;
      for (std::size_t k = 0; k < NQuarks; k++)
      {
        AddContQuark += -2 * QuarkMassesVec.at(k + NQuarks) *
                        fermion(QuarkMassesVec[k], Temp, diff);
      }
      for (std::size_t k = 0; k < NColour; k++)
        res += AddContQuark;
      for (std::size_t k = 0; k < NLepton; k++)
      {
        res += -2 * LeptonMassesVec.at(k + NLepton) *
               fermion(LeptonMassesVec[k], Temp, diff);
      }

      double VDebye = 0;
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        if (HiggsMassesVec[k] > 0)
        {
          VDebye += 1.5 * HiggsMassesVec.at(k + NHiggs) *
                    std::pow(HiggsMassesVec.at(k), 0.5);
        }
        if (HiggsMassesZeroTempVec[k] > 0)
        {
          VDebye += -1.5 * HiggsMassesZeroTempVec.at(k + NHiggs) *
                    std::pow(HiggsMassesZeroTempVec[k], 0.5);
        }
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        if (GaugeMassesVec[k] > 0)
        {
          VDebye += 1.5 * GaugeMassesVec.at(k + NGauge) *
                    std::pow(GaugeMassesVec[k], 0.5);
        }
        if (GaugeMassesZeroTempVec[k] > 0)
        {
          VDebye += -1.5 * GaugeMassesZeroTempVec.at(k + NGauge) *
                    std::pow(GaugeMassesZeroTempVec[k], 0.5);
        }
      }

      VDebye *= -Temp / (12 * M_PI);
      res += VDebye;
    }
  }
  else if (diff == -1)
  {
    if (C_UseParwani)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        res += HiggsMassesVec.at(k + NHiggs) *
               boson(HiggsMassesVec[k], Temp, C_CWcbHiggs, 0);
        res += boson(HiggsMassesVec[k], Temp, C_CWcbHiggs, -1);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        res += GaugeMassesVec.at(k + NGauge) *
               boson(GaugeMassesVec[k], Temp, C_CWcbGB, 0);
        res += boson(GaugeMassesVec.at(k), Temp, C_CWcbGB, -1);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        res += 2 * boson(GaugeMassesZeroTempVec[k], Temp, C_CWcbGB, -1);
      }
      for (std::size_t k = 0; k < NQuarks; k++)
      {
        res += -6 * fermion(QuarkMassesVec[k], Temp, -1);
      }
      for (std::size_t k = 0; k < NLepton; k++)
      {
        res += -2 * fermion(LeptonMassesVec[k], Temp, -1);
      }
    }
    else
    {
      HiggsMassesZeroTempVec = HiggsMassesSquared(v, 0, diff);
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        res += boson(HiggsMassesZeroTempVec[k], Temp, C_CWcbHiggs, -1);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        res += 3 * boson(GaugeMassesZeroTempVec[k], Temp, C_CWcbGB, -1);
      }
      double AddContQuark = 0;
      for (std::size_t k = 0; k < NQuarks; k++)
        AddContQuark += -2 * fermion(QuarkMassesVec[k], Temp, -1);
      for (std::size_t k = 0; k < NColour; k++)
        res += AddContQuark;
      for (std::size_t k = 0; k < NLepton; k++)
        res += -2 * fermion(LeptonMassesVec[k], Temp, -1);

      double VDebye = 0;
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        if (HiggsMassesVec[k] > 0) VDebye += std::pow(HiggsMassesVec[k], 1.5);
        if (HiggsMassesZeroTempVec[k] > 0)
          VDebye += -std::pow(HiggsMassesZeroTempVec[k], 1.5);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        if (GaugeMassesVec[k] > 0) VDebye += std::pow(GaugeMassesVec[k], 1.5);
        if (GaugeMassesZeroTempVec[k] > 0)
          VDebye += -std::pow(GaugeMassesZeroTempVec[k], 1.5);
      }

      VDebye *= -1.0 / (12 * M_PI);
      res += VDebye;
    }
  }

  return res;
}

void Class_Potential_Origin::CalculateDebye(bool forceCalculation)
{
  if (!SetCurvatureDone) SetCurvatureArrays();

  bool Calculate = forceCalculation or not CalculateDebyeSimplified();
  if (Calculate)
  {
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      for (std::size_t j = i; j < NHiggs; j++)
      {
        DebyeHiggs[i][j] = 0;
        for (std::size_t k = 0; k < NHiggs; k++)
        {
          DebyeHiggs[i][j] += 0.5 * Curvature_Higgs_L4[i][j][k][k] / 12.0;
        }
        for (std::size_t k = 0; k < NGauge; k++)
        {
          DebyeHiggs[i][j] += 3 * 0.5 * Curvature_Gauge_G2H2[k][k][i][j] / 12.0;
        }

        for (std::size_t a = 0; a < NQuarks; a++)
        {
          for (std::size_t b = 0; b < NQuarks; b++)
          {
            double tmp = 0.5 * (std::conj(Curvature_Quark_F2H1[a][b][j]) *
                                    Curvature_Quark_F2H1[a][b][i] +
                                std::conj(Curvature_Quark_F2H1[a][b][i]) *
                                    Curvature_Quark_F2H1[a][b][j])
                                   .real();
            DebyeHiggs[i][j] += 6.0 / 24.0 * tmp;
          }
        }

        for (std::size_t a = 0; a < NLepton; a++)
        {
          for (std::size_t b = 0; b < NLepton; b++)
          {
            double tmp = 0.5 * (std::conj(Curvature_Lepton_F2H1[a][b][j]) *
                                    Curvature_Lepton_F2H1[a][b][i] +
                                std::conj(Curvature_Lepton_F2H1[a][b][i]) *
                                    Curvature_Lepton_F2H1[a][b][j])
                                   .real();
            DebyeHiggs[i][j] += 2.0 / 24.0 * tmp;
          }
        }

        //	            if(i==j) DebyeHiggs[i][j] *= 0.5;
      }
    }

    for (std::size_t i = 0; i < NHiggs; i++)
    {
      for (std::size_t j = i; j < NHiggs; j++)
      {
        if (std::abs(DebyeHiggs[i][j]) <= 1e-5) DebyeHiggs[i][j] = 0;
      }
    }

    for (std::size_t i = 0; i < NHiggs; i++)
    {
      for (std::size_t j = 0; j < i; j++)
      {
        DebyeHiggs[i][j] = DebyeHiggs[j][i];
      }
    }
  }
}

void Class_Potential_Origin::CalculateDebyeGauge()
{
  for (std::size_t i = 0; i < NGauge; i++)
  {
    for (std::size_t j = 0; j < NGauge; j++)
      DebyeGauge[i][j] = 0;
  }

  bool Done = CalculateDebyeGaugeSimplified();
  if (Done) return;

  std::size_t nGaugeHiggs = 0;

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    if (Curvature_Gauge_G2H2[0][0][i][i] != 0)
    {
      nGaugeHiggs++;
    }
  }
  for (std::size_t i = 0; i < NGauge; i++)
  {
    double GaugeFac = 0;
    for (std::size_t k = 0; k < NHiggs; k++)
    {
      GaugeFac += Curvature_Gauge_G2H2[i][i][k][k];
    }
    GaugeFac *= 1.0 / nGaugeHiggs;
    DebyeGauge[i][i] = 2.0 / 3.0 * (nGaugeHiggs / 8.0 + 5) * GaugeFac;
  }

  for (std::size_t i = 0; i < NGauge; i++)
  {
    for (std::size_t j = 0; j < NGauge; j++)
    {
      if (std::abs(DebyeGauge[i][j]) <= 1e-5) DebyeGauge[i][j] = 0;
    }
  }
}

void Class_Potential_Origin::initVectors()
{
  using vec2 = std::vector<std::vector<double>>;
  using vec3 = std::vector<std::vector<std::vector<double>>>;
  using vec4 = std::vector<std::vector<std::vector<std::vector<double>>>>;

  using vec1Complex = std::vector<std::complex<double>>;
  using vec2Complex = std::vector<std::vector<std::complex<double>>>;
  using vec3Complex =
      std::vector<std::vector<std::vector<std::complex<double>>>>;
  using vec4Complex =
      std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>;

  Curvature_Higgs_L1 = std::vector<double>(NHiggs, 0);
  Curvature_Higgs_L2 = vec2{NHiggs, std::vector<double>(NHiggs, 0)};
  Curvature_Higgs_L3 =
      vec3{NHiggs, vec2{NHiggs, std::vector<double>(NHiggs, 0)}};
  Curvature_Higgs_L4 =
      vec4{NHiggs, vec3{NHiggs, vec2{NHiggs, std::vector<double>(NHiggs, 0)}}};

  Curvature_Higgs_CT_L1 = std::vector<double>(NHiggs, 0);
  Curvature_Higgs_CT_L2 = vec2{NHiggs, std::vector<double>(NHiggs, 0)};
  Curvature_Higgs_CT_L3 =
      vec3{NHiggs, vec2{NHiggs, std::vector<double>(NHiggs, 0)}};
  Curvature_Higgs_CT_L4 =
      vec4{NHiggs, vec3{NHiggs, vec2{NHiggs, std::vector<double>(NHiggs, 0)}}};

  VEVSymmetric = std::vector<double>(NHiggs, 0);

  DebyeHiggs = vec2{NHiggs, std::vector<double>(NHiggs, 0)};

  LambdaHiggs_3    = vec3{NHiggs, vec2{NHiggs, std::vector<double>(NHiggs, 0)}};
  LambdaHiggs_3_CT = vec3{NHiggs, vec2{NHiggs, std::vector<double>(NHiggs, 0)}};

  Curvature_Gauge_G2H2 =
      vec4{NGauge, vec3{NGauge, vec2{NHiggs, std::vector<double>(NHiggs, 0)}}};
  DebyeGauge    = vec2{NGauge, std::vector<double>(NGauge, 0)};
  LambdaGauge_3 = vec3{NGauge, vec2{NGauge, std::vector<double>(NHiggs, 0)}};

  Curvature_Lepton_F2 = vec2Complex{NLepton, vec1Complex(NLepton, 0)};
  Curvature_Lepton_F2H1 =
      vec3Complex{NLepton, vec2Complex{NLepton, vec1Complex(NHiggs, 0)}};
  LambdaLepton_3 =
      vec3Complex{NLepton, vec2Complex{NLepton, vec1Complex(NHiggs, 0)}};
  LambdaLepton_4 = vec4Complex{
      NLepton,
      vec3Complex{NLepton, vec2Complex{NHiggs, vec1Complex(NHiggs, 0)}}};

  Curvature_Quark_F2 = vec2Complex{NQuarks, vec1Complex(NQuarks, 0)};
  Curvature_Quark_F2H1 =
      vec3Complex{NQuarks, vec2Complex{NQuarks, vec1Complex(NHiggs, 0)}};
  LambdaQuark_3 =
      vec3Complex{NQuarks, vec2Complex{NQuarks, vec1Complex(NHiggs, 0)}};
  LambdaQuark_4 = vec4Complex{
      NQuarks,
      vec3Complex{NQuarks, vec2Complex{NHiggs, vec1Complex(NHiggs, 0)}}};

  HiggsVev = std::vector<double>(NHiggs, 0);

  HiggsRotationMatrixEnsuredConvention =
      std::vector<std::vector<double>>{NHiggs, std::vector<double>(NHiggs, 0)};
}

void Class_Potential_Origin::sym2Dim(
    std::vector<std::vector<double>> &Tensor2Dim,
    std::size_t Nk1,
    std::size_t Nk2)
{
  for (std::size_t k1 = 0; k1 < Nk1; k1++)
  {
    for (std::size_t k2 = k1; k2 < Nk2; k2++)
    {
      Tensor2Dim[k2][k1] = Tensor2Dim[k1][k2];
    }
  }
}

void Class_Potential_Origin::sym3Dim(
    std::vector<std::vector<std::vector<double>>> &Tensor3Dim,
    std::size_t Nk1,
    std::size_t Nk2,
    std::size_t Nk3)
{
  for (std::size_t k1 = 0; k1 < Nk1; k1++)
  {
    for (std::size_t k2 = k1; k2 < Nk2; k2++)
    {
      for (std::size_t k3 = k2; k3 < Nk3; k3++)
      {
        Tensor3Dim[k1][k3][k2] = Tensor3Dim[k1][k2][k3];
        Tensor3Dim[k2][k1][k3] = Tensor3Dim[k1][k2][k3];
        Tensor3Dim[k2][k3][k1] = Tensor3Dim[k1][k2][k3];
        Tensor3Dim[k3][k1][k2] = Tensor3Dim[k1][k2][k3];
        Tensor3Dim[k3][k2][k1] = Tensor3Dim[k1][k2][k3];
      }
    }
  }
}

void Class_Potential_Origin::sym4Dim(
    std::vector<std::vector<std::vector<std::vector<double>>>> &Tensor4Dim,
    std::size_t Nk1,
    std::size_t Nk2,
    std::size_t Nk3,
    std::size_t Nk4)
{
  for (std::size_t k1 = 0; k1 < Nk1; k1++)
  {
    for (std::size_t k2 = k1; k2 < Nk2; k2++)
    {
      for (std::size_t k3 = k2; k3 < Nk3; k3++)
      {
        for (std::size_t k4 = k3; k4 < Nk4; k4++)
        {
          Tensor4Dim[k1][k2][k4][k3] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k1][k3][k2][k4] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k1][k3][k4][k2] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k1][k4][k2][k3] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k1][k4][k3][k2] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k2][k1][k3][k4] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k2][k1][k4][k3] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k2][k3][k1][k4] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k2][k3][k4][k1] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k2][k4][k1][k3] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k2][k4][k3][k1] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k3][k1][k2][k4] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k3][k1][k4][k2] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k3][k2][k1][k4] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k3][k2][k4][k1] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k3][k4][k1][k2] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k3][k4][k2][k1] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k4][k1][k2][k3] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k4][k1][k3][k2] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k4][k2][k1][k3] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k4][k2][k3][k1] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k4][k3][k1][k2] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k4][k3][k2][k1] = Tensor4Dim[k1][k2][k3][k4];
        }
      }
    }
  }
}

void Class_Potential_Origin::resetbools()
{
  SetCurvatureDone          = false;
  CalcCouplingsDone         = false;
  CalculatedTripleCopulings = false;
  parStored.clear();
  parCTStored.clear();
}

bool Class_Potential_Origin::CheckNLOVEV(const std::vector<double> &v) const
{
  // std::vector<double> vPotential;
  double MaxDiff           = 0;
  double AllowedDifference = 1;
  for (std::size_t i = 0; i < nVEV; i++)
  {
    double tmp = std::abs(std::abs(v[i]) - std::abs(vevTreeMin[i]));
    if (tmp > MaxDiff) MaxDiff = tmp;
  }

  return (MaxDiff < AllowedDifference);
}

double Class_Potential_Origin::EWSBVEV(const std::vector<double> &v) const
{
  double res = 0;
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    double checkgauge = 0;
    for (std::size_t j = 0; j < NGauge; j++)
    {
      checkgauge += std::abs(Curvature_Gauge_G2H2[j][j][i][i]);
    }
    if (checkgauge != 0) res += std::pow(v.at(i), 2);
  }
  res = std::sqrt(res);

  if (res <= 0.5)
  {
    return 0;
  }
  return res;
}

void Class_Potential_Origin::SetEWVEVZero(std::vector<double> &sol) const
{
  int count = 0;
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    double checkgauge = 0;
    for (std::size_t j = 0; j < NGauge; j++)
    {
      checkgauge += std::abs(Curvature_Gauge_G2H2[j][j][i][i]);
    }
    if (checkgauge > 1e-10 and i == VevOrder[count]) // true for doublet vevs
    {
      sol.at(count) = 0;
      count += 1;
    }
  }
}

void Class_Potential_Origin::setUseIndexCol(std::string legend)
{
  UseIndexCol = legend.rfind(sep, 0) == 0;
}

bool Class_Potential_Origin::getUseIndexCol()
{
  return UseIndexCol;
}

void Class_Potential_Origin::FindSignSymmetries()
{
  SignSymmetries.clear();
  std::vector<double> testvev, testvevPotential;
  for (std::size_t i = 0; i < nVEV; i++)
    testvev.push_back(i + 1);
  testvevPotential      = MinimizeOrderVEV(testvev);
  double referenceValue = VEff(testvevPotential);

  std::vector<double> vevdummy, vevdummyPotential;
  double VEffDummy;

  // Fill a dummy vector with a certain amout of -1 and look at all possible
  // permutations of it
  for (std::size_t countNegative = 1; countNegative <= nVEV; countNegative++)
  {
    std::vector<double> tmpSymmetry;
    for (std::size_t i = 0; i < countNegative; i++)
      tmpSymmetry.push_back(-1);
    for (std::size_t i = countNegative; i < nVEV; i++)
      tmpSymmetry.push_back(1);

    do
    {
      vevdummy.clear();
      for (std::size_t i = 0; i < nVEV; i++)
        vevdummy.push_back(tmpSymmetry.at(i) * testvev.at(i));
      vevdummyPotential = MinimizeOrderVEV(vevdummy);
      VEffDummy         = VEff(vevdummyPotential);
      if (std::abs(VEffDummy - referenceValue) <=
          1e-3 * std::abs(referenceValue))
        SignSymmetries.push_back(tmpSymmetry);
    } while (std::next_permutation(tmpSymmetry.begin(), tmpSymmetry.end()));
  }
}

void Class_Potential_Origin::SetUseTreeLevel(bool val)
{
  UseTreeLevel = val;
}

std::pair<std::vector<double>, std::vector<double>>
Class_Potential_Origin::initModel(std::string linestr)
{
  std::vector<double> par(nPar), parCT(nParCT);
  resetbools();
  ReadAndSet(linestr, par);
  parCT = initModel(par);

  parStored   = par;
  parCTStored = parCT;

  std::pair<std::vector<double>, std::vector<double>> res;
  res.first  = par;
  res.second = parCT;

  return res;
}

std::vector<double>
Class_Potential_Origin::initModel(const std::vector<double> &par)
{
  std::vector<double> parCT(nParCT);
  resetbools();
  set_gen(par);
  CalculatePhysicalCouplings();
  parCT = calc_CT();
  set_CT_Pot_Par(parCT);
  CalculateDebye();
  CalculateDebyeGauge();

  AdjustRotationMatrix();

  parStored   = par;
  parCTStored = parCT;

  return parCT;
}

std::vector<double> Class_Potential_Origin::resetScale(const double &newScale)
{
  scale      = newScale;
  auto parCT = calc_CT();
  set_CT_Pot_Par(parCT);

  parCTStored = parCT;

  return parCT;
}

Eigen::MatrixXcd
Class_Potential_Origin::QuarkMassMatrix(const std::vector<double> &v) const
{
  MatrixXcd MIJ(NQuarks, NQuarks);
  if (v.size() != nVEV and v.size() != NHiggs)
  {
    std::string ErrorString =
        std::string("You have called ") + std::string(__func__) +
        std::string(
            " with an invalid vev configuration. Your vev is of dimension ") +
        std::to_string(v.size()) + std::string(" and it should be ") +
        std::to_string(NHiggs) + std::string(".");
    throw std::runtime_error(ErrorString);
  }
  if (v.size() == nVEV and nVEV != NHiggs)
  {
    std::stringstream ss;
    ss << __func__
       << " is being called with a wrong sized vev configuration. It "
          "has the dimension of "
       << nVEV << " while it should have " << NHiggs
       << ". For now this is transformed but please fix this to reduce "
          "the runtime."
       << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
    std::vector<double> Transformedv;
    Transformedv = MinimizeOrderVEV(v);
    MIJ          = QuarkMassMatrix(Transformedv);
    return MIJ;
  }
  if (!SetCurvatureDone)
  {
    std::string retmes = __func__;
    retmes += " is called before SetCurvatureArrays() is called. \n";
    throw std::runtime_error(retmes);
  }

  MIJ = MatrixXcd::Zero(NQuarks, NQuarks);

  for (std::size_t i = 0; i < NQuarks; i++)
  {
    for (std::size_t j = 0; j < NQuarks; j++)
    {
      MIJ(i, j) = Curvature_Quark_F2[i][j];
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        MIJ(i, j) += Curvature_Quark_F2H1[i][j][k] * v[k];
      }
    }
  }

  return MIJ;
}

std::vector<std::complex<double>>
Class_Potential_Origin::QuarkMasses(const std::vector<double> &v) const
{
  std::vector<std::complex<double>> res;
  double ZeroMass = std::pow(10, -10);

  auto MIJ = QuarkMassMatrix(v);

  ComplexEigenSolver<MatrixXcd> es(MIJ, false);

  for (std::size_t i = 0; i < NQuarks; i++)
  {
    auto tmp = es.eigenvalues()[i];
    if (std::abs(tmp) < ZeroMass)
      res.push_back(0);
    else
      res.push_back(tmp);
  }

  return res;
}

MatrixXcd
Class_Potential_Origin::LeptonMassMatrix(const std::vector<double> &v) const
{
  MatrixXcd res = MatrixXcd::Zero(NLepton, NLepton);
  if (v.size() != nVEV and v.size() != NHiggs)
  {
    std::string ErrorString =
        std::string("You have called ") + std::string(__func__) +
        std::string(
            " with an invalid vev configuration. Your vev is of dimension ") +
        std::to_string(v.size()) + std::string(" and it should be ") +
        std::to_string(NHiggs) + std::string(".");
    throw std::runtime_error(ErrorString);
  }
  if (v.size() == nVEV and nVEV != NHiggs)
  {
    std::stringstream ss;
    ss << __func__
       << " is being called with a wrong sized vev configuration. It "
          "has the dimension of "
       << nVEV << " while it should have " << NHiggs
       << ". For now this is transformed but please fix this to reduce "
          "the runtime."
       << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
    std::vector<double> Transformedv;
    Transformedv = MinimizeOrderVEV(v);
    res          = LeptonMassMatrix(Transformedv);
    return res;
  }
  if (!SetCurvatureDone)
  {
    std::string retmes = __func__;
    retmes += " is called before SetCurvatureArrays();\n";
    throw std::runtime_error(retmes);
  }

  for (std::size_t i = 0; i < NLepton; i++)
  {
    for (std::size_t j = 0; j < NLepton; j++)
    {
      res(i, j) = Curvature_Lepton_F2[i][j];
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        res(i, j) += Curvature_Lepton_F2H1[i][j][k] * v[k];
      }
    }
  }

  return res;
}

std::vector<std::complex<double>>
Class_Potential_Origin::LeptonMasses(const std::vector<double> &v) const
{
  std::vector<std::complex<double>> res;

  MatrixXcd MIJ = LeptonMassMatrix(v);

  ComplexEigenSolver<MatrixXcd> es(MIJ, false);
  double ZeroMass = 1e-10;
  for (std::size_t i = 0; i < NLepton; i++)
  {
    auto tmp = es.eigenvalues()[i];
    if (std::abs(tmp) < ZeroMass)
      res.push_back(0);
    else
      res.push_back(tmp);
  }
  return res;
}

double Class_Potential_Origin::CalculateRatioAlpha(
    const std::vector<double> &vev_symmetric,
    const std::vector<double> &vev_broken,
    const double &Temp) const
{
  (void)vev_symmetric;
  (void)vev_broken;
  (void)Temp;
  //  double res                          = 0;
  //  double PotentialSymmetricPhaseValue = VEff(vev_symmetric, Temp, 0);
  //  double PotentialSymmetricPhaseDeriv = VEff(vev_symmetric, Temp, -1);
  //  double PotentialBrokenPhaseValue    = VEff(vev_broken, Temp, 0);
  //  double PotentialBrokenPhaseDeriv    = VEff(vev_broken, Temp, -1);
  //  res = -(PotentialBrokenPhaseValue - PotentialSymmetricPhaseValue) +
  //        Temp * (PotentialBrokenPhaseDeriv - PotentialSymmetricPhaseDeriv);
  // TODO:: Unfinished!
  throw std::runtime_error("The CalculateRatioAlpha function is not finished "
                           "and you should not be using it!");
}

std::vector<double> Class_Potential_Origin::MinimizeOrderVEV(
    const std::vector<double> &vevMinimizer) const
{
  std::vector<double> vevFunction;

  std::size_t count = 0;
  for (std::size_t i = 0; i < NHiggs; ++i)
  {
    if (i == VevOrder[count])
    {
      vevFunction.push_back(vevMinimizer.at(count));
      count++;
    }
    else
      vevFunction.push_back(0);
  }
  return vevFunction;
}

Eigen::VectorXd
Class_Potential_Origin::NablaVCT(const std::vector<double> &v) const
{
  VectorXd result(NHiggs);
  for (std::size_t i{0}; i < NHiggs; ++i)
  {
    result(i) = Curvature_Higgs_CT_L1[i];
    for (std::size_t j{0}; j < NHiggs; ++j)
    {
      result(i) += Curvature_Higgs_CT_L2[i][j] * v.at(j);
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        result(i) += 0.5 * Curvature_Higgs_CT_L3[i][j][k] * v.at(j) * v.at(k);
        for (std::size_t l = 0; l < NHiggs; l++)
        {
          result(i) += 1.0 / 6.0 * Curvature_Higgs_CT_L4[i][j][k][l] * v.at(j) *
                       v.at(k) * v.at(l);
        }
      }
    }
  }
  return result;
}

Eigen::MatrixXd
Class_Potential_Origin::HessianCT(const std::vector<double> &v) const
{
  Eigen::MatrixXd result(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      result(i, j) = Curvature_Higgs_CT_L2[i][j];
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        result(i, j) += Curvature_Higgs_CT_L3[i][j][k] * v.at(k);
        for (std::size_t l = 0; l < NHiggs; l++)
        {
          result(i, j) +=
              0.5 * Curvature_Higgs_CT_L4[i][j][k][l] * v.at(k) * v.at(l);
        }
      }
    }
  }
  return result;
}

std::vector<double> Class_Potential_Origin::GetCTIdentities() const
{
  return std::vector<double>();
}

} // namespace BSMPT
