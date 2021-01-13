/*
 * SMParam.h
 *
 *  Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file Standard Model input parameters
 */
#ifndef SMPARAM_H_
#define SMPARAM_H_

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

const double theta12 = 13.04 / 180.0 * M_PI;
const double theta13 = 0.201 / 180.0 * M_PI;
const double theta23 = 2.38 / 180.0 * M_PI;
const double delta   = 1.20;

const double s12 = std::sin(theta12);
const double s13 = std::sin(theta13);
const double s23 = std::sin(theta23);
const double c12 = std::cos(theta12);
const double c13 = std::cos(theta13);
const double c23 = std::cos(theta23);

const std::complex<double> imagNumber(0, 1);
const std::complex<double> C_Vus = s12 * c13;
const std::complex<double> C_Vub = s13 * exp(-imagNumber * delta);
const std::complex<double> C_Vcd =
    s12 * c23 - c12 * s23 * s13 * exp(imagNumber * delta);
const std::complex<double> C_Vcs =
    -c12 * c23 - s12 * s23 * s13 * exp(imagNumber * delta);
const std::complex<double> C_Vtd =
    s12 * s23 - c12 * c23 * s13 * exp(imagNumber * delta);
const std::complex<double> C_Vts =
    -c12 * s23 - s12 * c23 * s13 * exp(imagNumber * delta);

const std::complex<double> C_Vud = c12 * c13;
const std::complex<double> C_Vtb = c23 * c13;
const std::complex<double> C_Vcb = s23 * c13;

const double C_MassW       = 80.385;
const double C_MassZ       = 91.1876;
const double C_MassSMHiggs = 125.09;

const double C_MassUp      = 0.1;
const double C_MassDown    = 0.1;
const double C_MassStrange = 0.1;
const double C_MassTop     = 172.5;
const double C_MassCharm   = 1.51;
const double C_MassBottom  = 4.92;

const double C_MassTau      = 1.77682;
const double C_MassMu       = 0.1056583715;
const double C_MassElectron = 0.510998928 * std::pow(10.0, -3.0);

const double C_GF    = 1.1663787 * 1e-5;
const double alphaEW = 1.0 / 128.862;

const double C_alpha_S = 0.119;
// const double C_sinsquaredWeinberg = 0.23126;
const double C_sinsquaredWeinberg =
    1 - (C_MassW * C_MassW) / (C_MassZ * C_MassZ);

const double C_vev0 = std::sqrt(1 / std::sqrt(2) * 1 / C_GF);
const double C_g    = 2 * C_MassW / C_vev0;
const double C_gs =
    2 * std::sqrt(std::pow(C_MassZ, 2) - std::pow(C_MassW, 2)) / C_vev0;
const double C_ElectricCharge = C_g * std::sqrt(C_sinsquaredWeinberg);

const double C_SMTriHiggs = 3 * C_MassSMHiggs * C_MassSMHiggs / (C_vev0);

#endif /* SMPARAM_H_ */
