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
 * @brief imaginary number i
 */
const std::complex<double> II(0, 1);

// CKM Matrix

// const double C_Vts = 0.0404;
// const double C_Vtd = 0.00867;
// const double C_Vtb = std::sqrt(1-C_Vts*C_Vts - C_Vtd*C_Vtd); //0.9991;
// const double C_Vcb = 0.0412;
// const double C_Vcd = 0.22520;
// const double C_Vcs = std::sqrt(1-C_Vcb*C_Vcb-C_Vcd*C_Vcd);//0.97344;
// const double C_Vub = 0.00351;
// const double C_Vus = 0.22534;
// const double C_Vud = std::sqrt(1-C_Vub*C_Vub-C_Vus*C_Vus);//0.97427;

// CKM Matrix as unitary

// const std::complex<double> C_Vts=0;
// const std::complex<double> C_Vtd = 0;
// const std::complex<double> C_Vtb = 1; //0.9991;
// const std::complex<double> C_Vcb = 0;
// const std::complex<double> C_Vcd = 0;
// const std::complex<double> C_Vcs = 1;//0.97344;
// const std::complex<double> C_Vub = 0;
// const std::complex<double> C_Vus = 0;
// const std::complex<double> C_Vud = 1; //0.97427;

/* Here is an example of the CKM Matrix given by the standard parameters. The
 * elements V11, V23 and V33 are real in this parametrisation and are calculated
 * by the other elements and the unitarity conditions. If the unitarity at
 * numerical precision is not given you will end up with massive charged
 * Goldstone bosons
 */

/**
 * @brief The lambda parameter in the Wolfenstein parametrisation of the
 * CKM-Matrix LHCHXSWG-INT-2015-006
 */
const double C_Wolfenstein_lambda = 0.22537;
/**
 * @brief The A parameter in the Wolfenstein parametrisation of the CKM-Matrix
 * LHCHXSWG-INT-2015-006
 */
const double C_Wolfenstein_A = 0.814;
/**
 * @brief The rho parameter in the Wolfenstein parametrisation of the CKM-Matrix
 * LHCHXSWG-INT-2015-006
 */
const double C_Wolfenstein_rho = 0.117;
/**
 * @brief The eta parameter in the Wolfenstein parametrisation of the CKM-Matrix
 * LHCHXSWG-INT-2015-006
 */
const double C_Wolfenstein_eta = 0.353;

/**
 * @brief The theta_12 mixing angle in the CKM matrix calculated by the
 * Wolfenstein parameters
 */
const double theta12 = std::asin(C_Wolfenstein_lambda);
/**
 * @brief The theta_23 mixing angle in the CKM matrix calculated by the
 * Wolfenstein parameters
 */
const double theta23 =
    std::asin(C_Wolfenstein_A * std::pow(C_Wolfenstein_lambda, 2));
/**
 * @brief The CP-violating angle in the CKM matrix calculated by the Wolfenstein
 * parameters
 */
const double delta =
    std::arg(C_Wolfenstein_A * std::pow(C_Wolfenstein_lambda, 3) *
             (C_Wolfenstein_rho + II * C_Wolfenstein_eta));
/**
 * @brief The theta_13 mixing angle in the CKM matrix calculated by the
 * Wolfenstein parameters
 */
const double theta13 =
    std::asin(std::abs(C_Wolfenstein_A * std::pow(C_Wolfenstein_lambda, 3) *
                       (C_Wolfenstein_rho + II * C_Wolfenstein_eta)));

/**
 * @brief ud element of the CKM matrix
 */
const std::complex<double> C_Vud = std::cos(theta12) * std::cos(theta13);
/**
 * @brief us element of the CKM matrix
 */
const std::complex<double> C_Vus = std::sin(theta12) * std::cos(theta13);
/**
 * @brief ub element of the CKM matrix
 */
const std::complex<double> C_Vub = std::sin(theta13) * std::exp(-delta * II);

/**
 * @brief cd element of the CKM matrix
 */
const std::complex<double> C_Vcd = -std::sin(theta12) * std::cos(theta23) -
                                   std::cos(theta12) * std::sin(theta23) *
                                       std::sin(theta13) * std::exp(II * delta);
/**
 * @brief cs element of the CKM matrix
 */
const std::complex<double> C_Vcs = std::cos(theta12) * std::cos(theta23) -
                                   std::sin(theta12) * std::sin(theta23) *
                                       std::sin(theta13) * std::exp(II * delta);
/**
 * @brief cb element of the CKM matrix
 */
const std::complex<double> C_Vcb = std::sin(theta23) * std::cos(theta13);

/**
 * @brief td element of the CKM matrix
 */
const std::complex<double> C_Vtd = std::sin(theta12) * std::sin(theta23) -
                                   std::cos(theta12) * std::cos(theta23) *
                                       std::sin(theta13) * std::exp(II * delta);
/**
 * @brief ts element of the CKM matrix
 */
const std::complex<double> C_Vts = -std::cos(theta12) * std::sin(theta23) -
                                   std::sin(theta12) * std::cos(theta23) *
                                       std::sin(theta13) * std::exp(II * delta);
/**
 * @brief tb element of the CKM matrix
 */
const std::complex<double> C_Vtb = std::cos(theta23) * std::cos(theta13);

/**
 * @brief Mass of the W-Boson
 * Unit: GeV
 * LHCHXSWG-INT-2015-006
 */
const double C_MassW = 80.385;
/**
 * @brief Mass of the Z-Boson
 * Unit: GeV
 * LHCHXSWG-INT-2015-006
 */
const double C_MassZ = 91.1876;
/**
 * @brief Mass of the Higgs Boson
 * Unit: GeV
 * LHCHXSWG-INT-2015-006
 */
const double C_MassSMHiggs = 125.09;

/**
 * @brief Mass of the up quark
 * Unit: GeV
 * https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG
 */
const double C_MassUp = 0.13;
/**
 * @brief Mass of the down quark
 * Unit: GeV
 * https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG
 */
const double C_MassDown = 0.14;
/**
 * @brief Mass of the strange quark
 * Unit: GeV
 * https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG
 */
const double C_MassStrange = 0.15;
/**
 * @brief Mass of the top quark in the OS Scheme
 * Unit: GeV
 * LHCHXSWG-INT-2015-006
 */
const double C_MassTop = 172.5;
/**
 * @brief Mass of the charm quark in the OS Scheme
 * Unit: GeV
 * LHCHXSWG-INT-2015-006
 */
const double C_MassCharm = 1.51;
/**
 * @brief Mass of the bottom quark in the OS Scheme
 * Unit: GeV
 * LHCHXSWG-INT-2015-006
 */
const double C_MassBottom = 4.92;

/**
 * @brief Mass of the tau lepton
 * Unit: GeV
 * LHCHXSWG-INT-2015-006
 */
const double C_MassTau = 1.77682;
/**
 * @brief Mass of the muon
 * Unit: GeV
 * LHCHXSWG-INT-2015-006
 */
const double C_MassMu = 0.1056583715;
/**
 * @brief Mass of the electron
 * Unit: GeV
 * LHCHXSWG-INT-2015-006
 */
const double C_MassElectron = 0.510998928 * std::pow(10.0, -3.0);

/**
 * @brief Fermi constant
 * Unit: GeV^{-2}
 * LHCHXSWG-INT-2015-006
 */
const double C_GF = 1.1663787 * 1e-5;

/**
 * @brief sin^2(theta_Weinberg) derived from the W- and Z-Boson masses
 */
const double C_sinsquaredWeinberg =
    1 - (C_MassW * C_MassW) / (C_MassZ * C_MassZ);

/**
 * @brief Vacuum expectation value of the SM derived through the Fermi constant
 */
const double C_vev0 = std::sqrt(1 / std::sqrt(2) * 1 / C_GF);
/**
 * @brief gauge coupling of the U(2)_L with the SM Higgs doublett, derived
 * through the W-Boson mass and the SM VEV Unit: GeV
 */
const double C_g = 2 * C_MassW / C_vev0;
/**
 * @brief gauge coupling of the U(1) with the SM Higgs doublett, derived through
 * the W- and Z-Boson masses and the SM VEV
 */
const double C_gs =
    2 * std::sqrt(std::pow(C_MassZ, 2) - std::pow(C_MassW, 2)) / C_vev0;

/**
 * @brief Trilinear coupling between three SM Higgs Boson, calculated as the
 * third derivative of the SM Higgs Potential Unit: GeV
 */
const double C_SMTriHiggs = 3 * C_MassSMHiggs * C_MassSMHiggs / (C_vev0);

} // namespace BSMPT

#endif /* SMPARAM_H_ */
