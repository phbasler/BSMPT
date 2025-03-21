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
  double C_Wolfenstein_lambda{0};

  /**
   * @brief The A parameter in the Wolfenstein parametrisation of the CKM-Matrix
   * LHCHXSWG-INT-2015-006
   */
  double C_Wolfenstein_A{0};

  /**
   * @brief The rho parameter in the Wolfenstein parametrisation of the
   * CKM-Matrix LHCHXSWG-INT-2015-006
   */
  double C_Wolfenstein_rho{0};
  /**
   * @brief The eta parameter in the Wolfenstein parametrisation of the
   * CKM-Matrix LHCHXSWG-INT-2015-006
   */
  double C_Wolfenstein_eta{0};

  /**
   * @brief The theta_12 mixing angle in the CKM matrix calculated by the
   * Wolfenstein parameters
   */
  double theta12{0};
  /**
   * @brief The theta_23 mixing angle in the CKM matrix calculated by the
   * Wolfenstein parameters
   */
  double theta23{0};
  /**
   * @brief The CP-violating angle in the CKM matrix calculated by the
   * Wolfenstein parameters
   */
  double delta{0};
  /**
   * @brief The theta_13 mixing angle in the CKM matrix calculated by the
   * Wolfenstein parameters
   */
  double theta13{0};

  /**
   * @brief ud element of the CKM matrix
   */
  std::complex<double> C_Vud{0, 0};
  /**
   * @brief us element of the CKM matrix
   */
  std::complex<double> C_Vus{0, 0};
  /**
   * @brief ub element of the CKM matrix
   */
  std::complex<double> C_Vub{0, 0};

  /**
   * @brief cd element of the CKM matrix
   */
  std::complex<double> C_Vcd{0, 0};
  /**
   * @brief cs element of the CKM matrix
   */
  std::complex<double> C_Vcs{0, 0};
  /**
   * @brief cb element of the CKM matrix
   */
  std::complex<double> C_Vcb{0, 0};

  /**
   * @brief td element of the CKM matrix
   */
  std::complex<double> C_Vtd{0, 0};
  /**
   * @brief ts element of the CKM matrix
   */
  std::complex<double> C_Vts{0, 0};
  /**
   * @brief tb element of the CKM matrix
   */
  std::complex<double> C_Vtb{0, 0};

  /**
   * @brief Mass of the W-Boson
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassW{0};
  /**
   * @brief Mass of the Z-Boson
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassZ{0};
  /**
   * @brief Mass of the Higgs Boson
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassSMHiggs{0};

  /**
   * @brief Mass of the up quark
   * Unit: GeV
   * https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG
   */
  double C_MassUp{0};
  /**
   * @brief Mass of the down quark
   * Unit: GeV
   * https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG
   */
  double C_MassDown{0};
  /**
   * @brief Mass of the strange quark
   * Unit: GeV
   * https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG
   */
  double C_MassStrange{0};
  /**
   * @brief Mass of the top quark in the OS Scheme
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassTop{0};
  /**
   * @brief Mass of the charm quark in the OS Scheme
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassCharm{0};
  /**
   * @brief Mass of the bottom quark in the OS Scheme
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassBottom{0};

  /**
   * @brief Mass of the tau lepton
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassTau{0};
  /**
   * @brief Mass of the muon
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassMu{0};
  /**
   * @brief Mass of the electron
   * Unit: GeV
   * LHCHXSWG-INT-2015-006
   */
  double C_MassElectron{0};

  /**
   * @brief Fermi constant
   * Unit: GeV^{-2}
   * LHCHXSWG-INT-2015-006
   */
  double C_GF{0};

  /**
   * @brief sin^2(theta_Weinberg) derived from the W- and Z-Boson masses
   */
  double C_sinsquaredWeinberg{0};

  /**
   * @brief Vacuum expectation value of the SM derived through the Fermi
   * constant
   */
  double C_vev0{0};
  /**
   * @brief gauge coupling of the U(2)_L with the SM Higgs doublett, derived
   * through the W-Boson mass and the SM VEV Unit: GeV
   */
  double C_g{0};
  /**
   * @brief gauge coupling of the U(1) with the SM Higgs doublett, derived
   * through the W- and Z-Boson masses and the SM VEV
   */
  double C_gs{0};

  /**
   * @brief Trilinear coupling between three SM Higgs Boson, calculated as the
   * third derivative of the SM Higgs Potential Unit: GeV
   */
  double C_SMTriHiggs{0};

  /**
   * @brief speed of sound in the bag model
   */
  const double Csound = 0.5773502691896258; // 1/sqrt(3)

  /**
   * @brief reduced Planck mass = MPl / (8 Pi)
   */
  const double MPl = 2.4e18;
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
