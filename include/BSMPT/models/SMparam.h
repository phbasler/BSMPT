// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * Standard Model input parameters
 */
#ifndef SMPARAM_H_
#define SMPARAM_H_

#include <cmath>
#include <complex>

namespace BSMPT
{

/**
 * @brief The ISMConstants struct containing all necessary SM constants.
 */
struct ISMConstants
{
  /**
   * @brief The lambda parameter in the Wolfenstein parametrisation of the
   * CKM-Matrix LHCHXSWG-INT-2015-006
   */
  double C_Wolfenstein_lambda;

  /**
   * @brief The A parameter in the Wolfenstein parametrisation of the CKM-Matrix
   * LHCHXSWG-INT-2015-006
   */
  double C_Wolfenstein_A;

  /**
   * @brief The rho parameter in the Wolfenstein parametrisation of the
   * CKM-Matrix LHCHXSWG-INT-2015-006
   */
  double C_Wolfenstein_rho;
  /**
   * @brief The eta parameter in the Wolfenstein parametrisation of the
   * CKM-Matrix LHCHXSWG-INT-2015-006
   */
  double C_Wolfenstein_eta;

  /**
   * @brief The theta_12 mixing angle in the CKM matrix calculated by the
   * Wolfenstein parameters
   */
  double theta12;
  /**
   * @brief The theta_23 mixing angle in the CKM matrix calculated by the
   * Wolfenstein parameters
   */
  double theta23;
  /**
   * @brief The CP-violating angle in the CKM matrix calculated by the
   * Wolfenstein parameters
   */
  double delta;
  /**
   * @brief The theta_13 mixing angle in the CKM matrix calculated by the
   * Wolfenstein parameters
   */
  double theta13;

  /**
   * @brief ud element of the CKM matrix
   */
  std::complex<double> C_Vud;
  /**
   * @brief us element of the CKM matrix
   */
  std::complex<double> C_Vus;
  /**
   * @brief ub element of the CKM matrix
   */
  std::complex<double> C_Vub;

  /**
   * @brief cd element of the CKM matrix
   */
  std::complex<double> C_Vcd;
  /**
   * @brief cs element of the CKM matrix
   */
  std::complex<double> C_Vcs;
  /**
   * @brief cb element of the CKM matrix
   */
  std::complex<double> C_Vcb;

  /**
   * @brief td element of the CKM matrix
   */
  std::complex<double> C_Vtd;
  /**
   * @brief ts element of the CKM matrix
   */
  std::complex<double> C_Vts;
  /**
   * @brief tb element of the CKM matrix
   */
  std::complex<double> C_Vtb;

  /**
   * @brief Mass of the W-Boson
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassW;
  /**
   * @brief Mass of the Z-Boson
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassZ;
  /**
   * @brief Mass of the Higgs Boson
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassSMHiggs;

  /**
   * @brief Mass of the up quark
   * Unit: GeV
   * https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG
   */
  double C_MassUp;
  /**
   * @brief Mass of the down quark
   * Unit: GeV
   * https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG
   */
  double C_MassDown;
  /**
   * @brief Mass of the strange quark
   * Unit: GeV
   * https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG
   */
  double C_MassStrange;
  /**
   * @brief Mass of the top quark in the OS Scheme
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassTop;
  /**
   * @brief Mass of the charm quark in the OS Scheme
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassCharm;
  /**
   * @brief Mass of the bottom quark in the OS Scheme
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassBottom;

  /**
   * @brief Mass of the tau lepton
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassTau;
  /**
   * @brief Mass of the muon
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassMu;
  /**
   * @brief Mass of the electron
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassElectron;

  /**
   * @brief Fermi constant
   * Unit: GeV^{-2}
   * LHCHXSWG-INT-2015-006
   */
  double C_GF;

  /**
   * @brief sin^2(theta_Weinberg) derived from the W- and Z-Boson masses
   */
  double C_sinsquaredWeinberg;

  /**
   * @brief Vacuum expectation value of the SM derived through the Fermi
   * constant
   */
  double C_vev0;
  /**
   * @brief gauge coupling of the U(2)_L with the SM Higgs doublett, derived
   * through the W-Boson mass and the SM VEV Unit: GeV
   */
  double C_g;
  /**
   * @brief gauge coupling of the U(1) with the SM Higgs doublett, derived
   * through the W- and Z-Boson masses and the SM VEV
   */
  double C_gs;

  /**
   * @brief Trilinear coupling between three SM Higgs Boson, calculated as the
   * third derivative of the SM Higgs Potential Unit: GeV
   */
  double C_SMTriHiggs;
};

/**
 * @brief imaginary number i
 */
const std::complex<double> II(0, 1);

/**
 * @brief GetSMConstants returns a set of SM contants as indicated by the
 * sources described for each parameter.
 * @return The SM Constants used by default in BSMPT
 */
const ISMConstants GetSMConstants();

} // namespace BSMPT

#endif /* SMPARAM_H_ */
