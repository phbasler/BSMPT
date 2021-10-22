// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * \mainpage BSMPT - A tool for calculating the temperature dependent  effective
potential in various extended Higgs models
 *
 * \section intro_sec Introduction
*BSMPT - Beyond the Standard Model Phase Transitions:
*BSMPT is a C++ tool BSMPT for the calculation of the strength of the
*   electroweak phase transition in extended Higgs sectors. For this the
*   loop-corrected effective potential at finite temperature is calculated
*   including the daisy resummation of the bosonic masses. The program
*   calculates the vacuum expectation value (VEV) v of the potential as a
*   function of the temperature, and in particular the critical VEV v_c at the
*   temperature T_c where the phase transition takes place. Furthermore, the
*   loop-corrected trilinear Higgs self-couplings are provided at T=0. The
*   renormalization scheme is chosen such that the loop-corrected
*   masses and mixing angles remain at their tree-level input values.
*   For details, see: https://arxiv.org/abs/1803.02846.
*   The code itself can be downloaded at https://github.com/phbasler/BSMPT.

 *  \section Citation
 *  If you use this program for your work please cite
https://arxiv.org/abs/1803.02846 and https://arxiv.org/abs/2007.01725

 * \section Installation
 * \subsection Dependencies
1. [GSL](http://www.gnu.org/software/gsl/) library: The code assume GSL is
installed in PATH
2. [Eigen3](http://eigen.tuxfamily.org) : If eigen3 is not found automatically,
you need to give the path to your eigen3 installation.  To install eigen3 you go
to your downloaded eigen3 folder and install it through

        mkdir build && cd build
        cmake -DCMAKE_INSTALL_PREFIX=/path/to/installedEigen  ..
        make install
After that you can use `-DEigen3_Dir=/path/to/installedEigen/share/eigen3/cmake`
to link Eigen3

3. [libcmaes](https://github.com/beniz/libcmaes): Additionally to GSL you should
either use libcmaes or NLopt. If libcmaes is not found automatically you can set
the path with

    `-Dlibcmaes_ROOT=/path/to/cmaes`

   If cmaes is not installed then it will be installed in your build directory.
For more details on the libcmaes installation, e.g. possible dependencies, visit
their [wiki](https://github.com/CMA-ES/libcmaes/wiki). If you don't want to
install or use it, you can set `-DUseLibCMAES=OFF`

4. [NLopt](https://nlopt.readthedocs.io/en/latest/): With
`-DNLopt_DIR=/path/to/installedNLopt/lib/cmake/nlopt` you need to tell where
NLopt is installed. If you do not then NLopt will not be used.
5. [Boost](https://www.boost.org/) : It should be found automatically, if not
you can use `-DBOOST_ROOT=/path/to/boost`

\subsection build
With the dependencies and options you can build the programm with

        mkdir build && cd build
        cmake (OPTIONS from Dependencies) ..
        cmake --build . -j


Note to Mac Users: You have to use the g++ compiler as clang does not support
OpenMP. If you get the error "missing libcmaes_config.h" during the compilation,
please check if the file is in build/libcmaes-0.9.5. If so, copy it to
/path/to/libcmaes/include/libcmaes


   \section executables_sec Executables
 *  We briefly give an overview of the available executables and their function.
*We illustrate the call of the executables for the CP-violationg 2HDM, with the
*input parameter values given in 'example/C2HDM_Input.dat.' For details, see
*also the section 'Executables' of the manual. Each executable can be called
with
* the option '--help' to get information on possible input parameters.
 *
 *  \subsection BSMPT
 *  BSMPT calculates the EWPT for every parameter point in the input file and
 *  gives out the results of those parameter points for which \f$v_c/T_c\f$ > 1.
 *  To find these points a bisection method is used with the temperature
starting
 *  between 0 and 300 GeV. The executable is called through:
 *
 *  	./bin/BSMPT Model Inputfile Outputfile LineStart LineEnd
 *
 *  This will call the specific model to be used, identified through 'Model',
and
 *  calculate the EWPT for each parameter point (corresponding to one line)
between
 *  'LineStart' and 'LineEnd'.
 *
 *	For our example the command
 *
 *  	./bin/BSMPT c2hdm example/C2HDM_Input.dat example/test_BSMPT.dat 2 2
 *
 *  will calculate for the C2HDM the EWPT for one
 *  parameter point given in line 2 in C2HDM_Input.dat. This will generate the
output file example/test_BSMPT.dat 2 2
 *  which can be compared with the already available file
example/C2HDM_Input.dat_BSMPT.
 *
 *
 * \subsection CalcCT
 * This will calculate the counterterms. In the output file the information on
 * the input parameter point is given and the counterterms are added at the end
 * of the line.
 *
 * It is called through the command
 *
 *  	./bin/CalcCT Model Inputfile Outputfile LineStart LineEnd
 *
 * For the C2HDM example this reads
 *
 * 		./bin/CalcCT c2hdm example/C2HDM_Input.dat example/test_CalcCT.dat 2 2
 *
 * which will generate the output file example/test_CalcCT.dat. This can be
compared with the already available file
 * example/C2HDM_Input.dat_CalcCT.
 *
 *  \subsection NLOVEV
 *
 * This calculates the VEV at 1-loop order at vanishing temperature in the
 * effective potential approach. This can be used to investigate the vacuum
 * stability of the model. It is called through
 *
 *  	./bin/NLOVEV Model Inputfile Outputfile LineStart LineEnd
 *
 * and for the C2HDM example it is given by
 *
 *  	./bin/NLOVEV c2hdm example/C2HDM_Input.dat example/test_NLOVEV.dat 2 2
 *
 * where the result is written into the file example/test_NLOVEV.dat which
 * can be compared with the already available file
example/C2HDM_Input.dat_NLOVEV.
 *
 * \subsection VEVEVO
 * This program calculates the evolution of the vacuum expecation value of a
given point with the temperature.
 * It is called through
 *
 *	  ./bin/VEVEVO Model Inputfile Outputfile Line Tempstart Tempstep Tempend
 *
 * where 'Tempstart' is the starting value of the temperature which increases
with 'Tempstep' until 'Tempend'.
 *
 * For our C2HDM example this would be
 *
 *  	./bin/VEVEVO c2hdm example/C2HDM_Input.dat example/test_VEVEVO.dat 2 0 5
        150
 *
 * where the result for the NLO VEV is given in example/test_VEVEVO.dat as
 * function of the temperature in the interval between 0 and 150 GeV in steps of
 * 5 GeV. This can be compared with the already available file
example/C2HDM_Input.dat_vevevo.
 *
 *
 *  \subsection TripleHiggsCouplingNLO
 *
 * This program calculates the trilinear Higgs self-couplings at NLO at zero
 * temperature. It is called through
 *
 * 	 	./bin/TripleHiggsNLO Model Inputfile Outputfile LineStart LineEnd
 *
 * The C2HDM example is called through
 *
 * 	 	./bin/TripleHiggsNLO c2hdm example/C2HDM_Input.dat
example/test_TripleHiggsCouplingNLO.dat 2 2
 *
 * with the result given in example/test_TripleHiggsNLO.dat which
 * can be compared with the already available file
example/C2HDM_Input.dat_TripleHiggsCouplingNLO .
 *
 * \subsection CalculateEWBG
 *
 * This program calculates the difference between baryons and anti-baryons
normalized to the photon density generated through the EWPT.
 *  Please beware that this is only tested for the C2HDM so far and the general
implementation is future work. It is called through
 *
 *  	./bin/CalculateEWBG c2hdm Inputfile Outputfile LineStart LineEnd
        config_file

 *
 * An example is given for the example/C2HDM_Input.dat parameter point through
 *
 *      ./bin/CalculateEWBG c2hdm example/C2HDM_Input.dat example/test_EWBG.dat
        2 2 example/EWBG_config.txt
 *
 * with the result given in example/test_EWBG.dat which can be compared with
the already available file example/C2HDM_Input.dat_EWBG.
 *
 *
 * \subsection PlotEWBG_vw
 * This executable varies the wall velocity of a given parameter point and
calculates the EWBG for each velocity.
 *
 * \subsection PlotEWBG_nL
 * This executable calculates the left-handed fermion density in front of the
wall as a function of the distance to the bubble wall.
 *
 *
 *
 */

#ifndef CLASSPOTENTIALORIGIN_H_
#define CLASSPOTENTIALORIGIN_H_

#include "Eigen/Eigenvalues"
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/SMparam.h>
#include <iostream>
#include <vector>

namespace BSMPT
{

/**
 * @brief C_UseParwani Use the Parwani Method instead of Arnold-Espinosa
 */
const bool C_UseParwani = false;

/**
 * @brief C_PT Lower threshold to stop the EWPT calculation
 */
const double C_PT = 0;

/**
 * @brief C_threshold threshold to check if a mass is numerically zero
 */
const double C_threshold = 1e-4;
/**
 * @brief C_CWcbFermion constant used in the CW potential for fermions
 */
static constexpr double C_CWcbFermion = 1.5;
/**
 * @brief C_CWcbGB constant used in the CW potential for gauge bosons
 */
const double C_CWcbGB = 5.0 / 6.0;
/**
 * @brief C_CWcbHiggs constant used in the CW potential for Higgs bosons
 */
const double C_CWcbHiggs = 1.5;

/**
 * @brief The Class_Potential_Origin class
 * Base class for all models. This class contains all numerical calculations on
 * the tensors and the inherited classes only have to set them.
 */
class Class_Potential_Origin
{
protected:
  /**
   * @brief UseTreeLevel Enforces VEff to only use the tree-level potential
   */
  bool UseTreeLevel = false;

  /**
   * MSBar renormalization scale
   */
  double scale = C_vev0;

  /**
   * Number of Lagrange parameters in the Higgs Tree-Level potential
   */
  std::size_t nPar = 0;
  /**
   * Number of Lagrange parameters in the Higgs Counterterm potential
   */
  std::size_t nParCT = 0;

  /**
   * Vector to store the parameter of the Potential
   */
  std::vector<double> parStored;

  /**
   * Vector to store the parameter of the Counterterm Potential
   */
  std::vector<double> parCTStored;

  /**
   * Parameter to differentiate the models
   */
  ModelID::ModelIDs Model = ModelID::ModelIDs::NotSet;
  /**
   * Number of neutral Higgs bosons
   */
  std::size_t NNeutralHiggs = 0;
  /**
   * Number of charged Higgs bosons
   */
  std::size_t NChargedHiggs = 0;
  /**
   * Number of all Higgs particles. This will define the size of your Higgs mass
   * matrix.
   */
  std::size_t NHiggs = NNeutralHiggs + NChargedHiggs;
  /**
   * Number of gauge bosons. Do not change this in the current version as we
   * only investigate extended Higgs sectors. If you want to extend the other
   * sectors as well the Debye corrections have to be calculated by hand.
   */
  std::size_t NGauge = 4;
  /**
   * Number of quarks. Do not change this in the current version as we only
   * investigate extended Higgs sectors. If you want to extend the other sectors
   * as well the Debye corrections have to be calculated by hand.
   */
  std::size_t NQuarks = 12;

  /**
   * @brief NColour Number of colours of the quarks.
   * Do not change this in the current version as we only investigate extended
   * Higgs sectors. If you want to extend the other sectors as well the Debye
   * corrections have to be calculated by hand.
   */
  std::size_t NColour = 3;

  /**
   * Number of leptons. Do not change this in the current version as we only
   * investigate extended Higgs sectors. If you want to extend the other sectors
   * as well the Debye corrections have to be calculated by hand.
   */
  std::size_t NLepton = 9;

  /**
   * Number of VEVs you want to minimize.
   */
  std::size_t nVEV = 0;

  /**
   * Storage of the symmetric VEV in all Higgs fields
   */
  std::vector<double> VEVSymmetric;
  /**
   * Storage of the tree-level VEV in all Higgs fields
   */
  std::vector<double> vevTree;
  /**
   * Storage of the tree-Level VEV in the configuration for the minimizer
   */
  std::vector<double> vevTreeMin;
  /**
   * Storage of the contributions of the Coleman-Weinberg potential to the
   * triple Higgs couplings in the gauge basis
   */
  std::vector<double> TripleHiggsCorrectionsCW;
  /**
   * Storage of the contributions of the Coleman-Weinberg potential to the
   * triple Higgs couplings in the mass basis
   */
  std::vector<std::vector<std::vector<double>>>
      TripleHiggsCorrectionsCWPhysical;
  /**
   * Storage of the contributions of the tree-level potential to the triple
   * Higgs couplings in the mass basis
   */
  std::vector<std::vector<std::vector<double>>>
      TripleHiggsCorrectionsTreePhysical;
  /**
   * Storage of the contributions of the counterterm potential to the triple
   * Higgs couplings in the mass basis
   */
  std::vector<std::vector<std::vector<double>>>
      TripleHiggsCorrectionsCTPhysical;

  /**
   * @brief SetCurvatureDone Used to check if the tensors are set
   */
  bool SetCurvatureDone = false;
  /**
   * @brief CalcCouplingsdone Used to check if CalculatePhysicalCouplings has
   * already been called
   */
  bool CalcCouplingsdone = false;
  /**
   * @brief CalculatedTripleCopulings Used to check if TripleHiggsCouplings has
   * already been called
   */
  bool CalculatedTripleCopulings = false;

  /**
   * @brief InputLineNumber Used for the error message in fbaseTri
   */
  int InputLineNumber = -1;

  /**
   * Variable to check if the input file has an index column or not
   */
  bool UseIndexCol = false;

  /**
   * Storage of the Tree-Level Higgs VEV
   */
  std::vector<double> HiggsVev;
  /**
   * L_{(S)}^{i}
   */
  std::vector<double> Curvature_Higgs_L1;
  /**
   * L_{(S)}^{ij}
   */
  std::vector<std::vector<double>> Curvature_Higgs_L2;
  /**
   * L_{(S)}^{ijk}
   */
  std::vector<std::vector<std::vector<double>>> Curvature_Higgs_L3;
  /**
   * L_{(S)}^{ijkl}
   */
  std::vector<std::vector<std::vector<std::vector<double>>>> Curvature_Higgs_L4;
  /**
   * L_{(S),CT}^{i}
   */
  std::vector<double> Curvature_Higgs_CT_L1;
  /**
   * L_{(S),CT}^{ij}
   */
  std::vector<std::vector<double>> Curvature_Higgs_CT_L2;
  /**
   * L_{(S),CT}^{ijk}
   */
  std::vector<std::vector<std::vector<double>>> Curvature_Higgs_CT_L3;
  /**
   * L_{(S),CT}^{ijkl}
   */
  std::vector<std::vector<std::vector<std::vector<double>>>>
      Curvature_Higgs_CT_L4;
  /**
   * G^{abij}
   */
  std::vector<std::vector<std::vector<std::vector<double>>>>
      Curvature_Gauge_G2H2;
  /**
   * Y^{IJk} for Quarks
   */
  std::vector<std::vector<std::vector<std::complex<double>>>>
      Curvature_Quark_F2H1;
  /**
   * Y^{IJ} for Quarks
   */
  std::vector<std::vector<std::complex<double>>> Curvature_Quark_F2;
  /**
   * Y^{IJk} for Leptons
   */
  std::vector<std::vector<std::vector<std::complex<double>>>>
      Curvature_Lepton_F2H1;
  /**
   * Y^{IJ} for Leptons
   */
  std::vector<std::vector<std::complex<double>>> Curvature_Lepton_F2;
  /**
   * @brief MassSquaredHiggs Stores the masses of the Higgs Bosons calculated in
   * CalculatePhysicalCouplings
   */
  std::vector<double> MassSquaredHiggs;
  /**
   * @brief MassSquaredGauge Stores the masses of the gauge Bosons calculated in
   * CalculatePhysicalCouplings
   */
  std::vector<double> MassSquaredGauge;
  /**
   * @brief MassSquaredQuark Stores the masses of the quarks calculated in
   * CalculatePhysicalCouplings
   */
  std::vector<double> MassSquaredQuark;
  /**
   * @brief MassSquaredLepton Stores the masses of the leptons calculated in
   * CalculatePhysicalCouplings
   */
  std::vector<double> MassSquaredLepton;
  /**
   * Storage of the Higgs rotation matrix for the Higgs mass matrix at the
   * tree-level Vacuum
   */
  std::vector<std::vector<double>> HiggsRotationMatrix;
  /**
   * @brief Couplings_Higgs_Quartic Stores the quartic Higgs couplings in the
   * mass base
   */
  std::vector<std::vector<std::vector<std::vector<double>>>>
      Couplings_Higgs_Quartic;
  /**
   * @brief Couplings_Higgs_Triple Stores the triple Higgs couplings in the mass
   * base
   */
  std::vector<std::vector<std::vector<double>>> Couplings_Higgs_Triple;
  /**
   * @brief Couplings_Gauge_Higgs_22 Stores the couplings between two Higgs and
   * two gauge bosons in the mass base
   */
  std::vector<std::vector<std::vector<std::vector<double>>>>
      Couplings_Gauge_Higgs_22;
  /**
   * @brief Couplings_Gauge_Higgs_21 Stores the coupling between two gauge and
   * one Higgs boson in the mass base
   */
  std::vector<std::vector<std::vector<double>>> Couplings_Gauge_Higgs_21;
  /**
   * @brief Couplings_Quark_Higgs_22 Stores the couplings between two quarks and
   * two Higgs bosons in the mass base
   */
  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>
      Couplings_Quark_Higgs_22;
  /**
   * @brief Couplings_Quark_Higgs_21 Stores the couplings between two quarks and
   * one Higgs boson in the mass base
   */
  std::vector<std::vector<std::vector<std::complex<double>>>>
      Couplings_Quark_Higgs_21;
  /**
   * @brief Couplings_Lepton_Higgs_22 Stores the couplings between two leptons
   * and two Higgs bosons in the mass base
   */
  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>
      Couplings_Lepton_Higgs_22;
  /**
   * @brief Couplings_Quark_Higgs_21 Stores the couplings between two leptons
   * and one Higgs boson in the mass base
   */
  std::vector<std::vector<std::vector<std::complex<double>>>>
      Couplings_Lepton_Higgs_21;

  /**
   * @brief LambdaGauge_3 Stores the Lambda_{(G)}^{abi} tensor
   */
  std::vector<std::vector<std::vector<double>>> LambdaGauge_3;
  /**
   * @brief LambdaHiggs_3 Stores the Lambda_{(S)}^{ijk} tensor
   */
  std::vector<std::vector<std::vector<double>>> LambdaHiggs_3;
  /**
   * @brief LambdaHiggs_3 Stores the Lambda_{(S)}^{ijk} tensor for the
   * counterterm parameters
   */
  std::vector<std::vector<std::vector<double>>> LambdaHiggs_3_CT;
  /**
   * @brief LambdaQuark_3 Stores the Lambda_{(F)}^{IJk} tensor for quarks ,
   * describing the derivative of Lambda_{(F)}^{IJ} w.r.t. the Higgs field k
   */
  std::vector<std::vector<std::vector<std::complex<double>>>> LambdaQuark_3;
  /**
   * @brief LambdaLepton_3 Stores the Lambda_{(F)}^{IJk} tensor for leptons ,
   * describing the derivative of Lambda_{(F)}^{IJ} w.r.t. the Higgs field k
   */
  std::vector<std::vector<std::vector<std::complex<double>>>> LambdaLepton_3;
  /**
   * @brief LambdaQuark_4 Stores the Lambda_{(F)}^{IJkm} tensor for quarks ,
   * describing the derivative of Lambda_{(F)}^{IJ} w.r.t. the Higgs fields k
   * and m
   */
  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>
      LambdaQuark_4;
  /**
   * @brief LambdaLepton_4 Stores the Lambda_{(F)}^{IJkm} tensor for leptons ,
   * describing the derivative of Lambda_{(F)}^{IJ} w.r.t. the Higgs fields k
   * and m
   */
  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>
      LambdaLepton_4;

  /**
   * @brief DebyeHiggs Stores the debye corrections to the mass matrix of the
   * Higgs bosons
   */
  std::vector<std::vector<double>> DebyeHiggs;
  /**
   * @brief DebyeGauge Stores the debye corrections to the mass matrix of the
   * gauge bosons
   */
  std::vector<std::vector<double>> DebyeGauge;
  /**
   * @brief VevOrder Stores the matching order used in MinimizeOrderVEV, set in
   * the constructor of the model
   */
  std::vector<std::size_t> VevOrder;

public:
  Class_Potential_Origin();
  virtual ~Class_Potential_Origin();

  /**
   * @brief Get_Curvature_Higgs_L2 allows const ref access to the L_{(S)}^{ij}
   * Tensor
   * @return
   */
  const auto &Get_Curvature_Higgs_L2() const { return Curvature_Higgs_L2; }

  /**
   * @brief Get_Curvature_Higgs_L3 allows const ref access to the L_{(S)}^{ijk}
   * Tensor
   * @return
   */
  const auto &Get_Curvature_Higgs_L3() const { return Curvature_Higgs_L3; }

  /**
   * @brief Get_Curvature_Higgs_L4 allows const ref access to the L_{(S)}^{ijkl}
   * Tensor
   * @return
   */
  const auto &Get_Curvature_Higgs_L4() const { return Curvature_Higgs_L4; }

  /**
   * @brief Get_Curvature_Gauge_G2H2 allows const ref access to the G^{abij}
   * Tensor
   * @return
   */
  const auto &Get_Curvature_Gauge_G2H2() const { return Curvature_Gauge_G2H2; }

  /**
   * @brief Get_Curvature_Lepton_F2H1 allows const ref access to the Y^{IJk}
   * Tensor for Leptons
   * @return
   */
  const auto &Get_Curvature_Lepton_F2H1() const
  {
    return Curvature_Lepton_F2H1;
  }

  /**
   * @brief Get_Curvature_Lepton_F2 allows const ref access to the Y^{IJ}
   * Tensor for Leptons
   * @return
   */
  const auto &Get_Curvature_Lepton_F2() const { return Curvature_Lepton_F2; }

  /**
   * @brief Get_Curvature_Quark_F2H1 allows const ref access to the Y^{IJk}
   * Tensor for Quarks
   * @return
   */
  const auto &Get_Curvature_Quark_F2H1() const { return Curvature_Quark_F2H1; }

  /**
   * @brief Get_Curvature_Quark_F2 allows const ref access to the Y^{IJ}
   * Tensor for Quarks
   * @return
   */
  const auto &Get_Curvature_Quark_F2() const { return Curvature_Quark_F2; }

  /**
   * @brief get_scale
   * @return the MSBar renormalisation scale
   */
  double get_scale() const { return scale; }

  /**
   * @brief set_scale sets the MSBar renormalisation scale to scale_new
   * @param scale_new
   */
  void set_scale(double scale_new) { scale = scale_new; }
  /**
   * @brief get_nPar
   * @return nPar
   */
  std::size_t get_nPar() const { return nPar; }
  /**
   * @brief get_nParCT
   * @return nParCT
   */
  std::size_t get_nParCT() const { return nParCT; }
  /**
   * @brief get_nVEV
   * @return nVEV
   */
  std::size_t get_nVEV() const { return nVEV; }
  /**
   * @brief get_vevTreeMin
   * @return vevTreeMin
   */
  std::vector<double> get_vevTreeMin() const { return vevTreeMin; }
  /**
   * @brief get_vevTreeMin
   * @param k
   * @return vevTreeMin.at(k)
   */
  double get_vevTreeMin(const std::size_t &k) const { return vevTreeMin.at(k); }
  /**
   * @brief get_parStored
   * @return parStored
   */
  std::vector<double> get_parStored() const { return parStored; }
  /**
   * @brief get_parCTStored
   * @return parCTStored
   */
  std::vector<double> get_parCTStored() const { return parCTStored; }
  /**
   * @brief get_NGauge
   * @return NGauge
   */
  std::size_t get_NGauge() const { return NGauge; }
  /**
   * @brief get_NQuarks
   * @return NQuarks
   */
  std::size_t get_NQuarks() const { return NQuarks; }
  /**
   * @brief get_NHiggs
   * @return NHiggs
   */
  std::size_t get_NHiggs() const { return NHiggs; }
  /**
   * @brief get_NLepton
   * @return NLepton
   */
  std::size_t get_NLepton() const { return NLepton; }
  /**
   * @brief get_Model
   * @return ModelID of the Model
   */
  ModelID::ModelIDs get_Model() const { return Model; }
  /**
   * @brief set_InputLineNumber
   * @param InputLineNumber_in value to set InputLineNumber
   */
  void set_InputLineNumber(int InputLineNumber_in)
  {
    InputLineNumber = InputLineNumber_in;
  }
  /**
   * @brief get_TripleHiggsCorrectionsTreePhysical
   * @param i
   * @param j
   * @param k
   * @return TripleHiggsCorrectionsTreePhysical[i][j][k]
   */
  double get_TripleHiggsCorrectionsTreePhysical(std::size_t i,
                                                std::size_t j,
                                                std::size_t k) const
  {
    return TripleHiggsCorrectionsTreePhysical.at(i).at(j).at(k);
  }
  /**
   * @brief get_TripleHiggsCorrectionsCTPhysical
   * @param i
   * @param j
   * @param k
   * @return TripleHiggsCorrectionsCTPhysical[i][j][k]
   */
  double get_TripleHiggsCorrectionsCTPhysical(std::size_t i,
                                              std::size_t j,
                                              std::size_t k) const
  {
    return TripleHiggsCorrectionsCTPhysical.at(i).at(j).at(k);
  }
  /**
   * @brief get_TripleHiggsCorrectionsCWPhysical
   * @param i
   * @param j
   * @param k
   * @return TripleHiggsCorrectionsCWPhysical[i][j][k]
   */
  double get_TripleHiggsCorrectionsCWPhysical(std::size_t i,
                                              std::size_t j,
                                              std::size_t k) const
  {
    return TripleHiggsCorrectionsCWPhysical.at(i).at(j).at(k);
  }

  /**
   * Sets the UseIndexCol var
   */
  void setUseIndexCol(std::string legend);

  /**
   * Initializes all vectors needed for the calculations.
   */
  void initVectors();

  /**
   * Calculates the effective potential and its derivatives.
   * @param v vev configuration at which the potential should be evaluated
   * @param Temp temperature at which the potential should be evaluated
   * @param diff Switch for the derivative of the potential. Default is 0 for
   * the value of the potential
   * @param Order 0 returns the tree level potential and 1 the NLO potential.
   * Default value is the NLO potential
   */
  double VEff(const std::vector<double> &v,
              double Temp = 0,
              int diff    = 0,
              int Order   = 1) const;
  /**
   * Calculates the tree-level potential and its derivatives.
   * @param v the configuration of all VEVs at which the potential should be
   * calculated
   * @param diff 0 returns the potential and i!= 0 returns the derivative of the
   * potential w.r.t v_i
   * @param ForceExplicitCalculation Calculate the tensors directly from the
   * tensors even if VTreeSimplified() is given
   */
  double VTree(const std::vector<double> &v,
               int diff                      = 0,
               bool ForceExplicitCalculation = false) const;
  /**
   * Calculates the counterterm potential and its derivatives
   * @param v the configuration of all VEVs at which the potential should be
   * calculated
   * @param diff 0 returns the potential and i!= 0 returns the derivative of the
   * potential w.r.t v_i
   * @param ForceExplicitCalculation Calculate the tensors directly from the
   * tensors even if VCounterSimplified() is given
   */
  double CounterTerm(const std::vector<double> &v,
                     int diff                      = 0,
                     bool ForceExplicitCalculation = false) const;
  /**
   * Calculates the Coleman-Weinberg and temperature-dependent 1-loop part of
   * the effective potential and its derivatives.
   * @param diff 0: No derivative, > 0 derivative w.r.t to the Higgs field; -1 :
   * derivative w.r.t to the temperature
   * @param v the configuration of all VEVs at which the potential should be
   * calculated
   * @param Temp the temperature at which the potential should be evaluated
   * @return the value of the one-loop part of the effective potential
   */
  double V1Loop(const std::vector<double> &v, double Temp, int diff) const;
  /**
   * This function calculates the EW breaking VEV from all contributing field
   * configurations.
   */
  double EWSBVEV(const std::vector<double> &v) const;

  /**
   * Reads the string linestr and sets the parameter point
   */
  virtual void ReadAndSet(const std::string &linestr,
                          std::vector<double> &par) = 0;
  /**
   * Adds the name of the counterterms to the legend of the output file. This
   * has to be specified in the model file.
   */
  virtual std::vector<std::string> addLegendCT() const = 0;
  /**
   * Adds the name of the VEVs, together with the critical VEV, the critical
   * temperature and the ratio, to the legend of the output file. This has to be
   * specified in the model file.
   */
  virtual std::vector<std::string> addLegendTemp() const = 0;
  /**
   * Adds the name of the triple Higgs couplings at T=0 to the legend of the
   * output file. This has to be specified in the model file.
   */
  virtual std::vector<std::string> addLegendTripleCouplings() const = 0;
  /**
   * Adds the name of the VEVs to the legend of the output file. This has to be
   * specified in the model file.
   */
  virtual std::vector<std::string> addLegendVEV() const = 0;

  /**
   * Reads the Lagrangian parameters from the vector 'par' and sets them to the
   * model parameters. This also sets the VEV configuration as well as the
   * scale. This has to be specified in the model file.
   */
  virtual void set_gen(const std::vector<double> &par) = 0;
  /**
   * Reads the counterterm parameters from the vector 'par' and sets them to the
   * model parameters as well as the Tensors \p Curvature_Higgs_CT_L1, \p
   * Curvature_Higgs_CT_L2, \p Curvature_Higgs_CT_L3 , \p Curvature_Higgs_CT_L4
   * . This has to be specified in the model file.
   */
  virtual void set_CT_Pot_Par(const std::vector<double> &par) = 0;
  /**
   * This will produce a terminal output of all the model parameters.
   * This has to be specified in the model file.
   */
  virtual void write() const = 0;

  void set_All(const std::vector<double> &par,
               const std::vector<double> &parCT);

  /**
   * This will set all the tensors needed to describe the tree-level Lagrangian
   *  except for the counterterms in the potential.
   * This has to be specified in the model file.
   */
  virtual void SetCurvatureArrays() = 0;
  /**
    Calculates all triple and quartic couplings in the physical basis
 */
  void CalculatePhysicalCouplings();
  /**
   * Calculates the first derivative of the Coleman-Weinberg potential evaluated
   * at the tree-level minimum.
   */
  std::vector<double> WeinbergFirstDerivative() const;
  /**
   * Calculates the second derivative of the Coleman-Weinberg potential at the
   * tree-level minimum.
   */
  std::vector<double> WeinbergSecondDerivative() const;
  /**
   * Calculates the second derivative of the Coleman-Weinberg potential at the
   * tree-level minimum.
   */
  Eigen::MatrixXd WeinbergSecondDerivativeAsMatrixXd() const;
  /**
   * Calculates the third derivative of the Coleman-Weinberg potential at the
   * tree-level minimum.
   */
  std::vector<double> WeinbergThirdDerivative() const;
  /**
   * Calculates the forth derivative of the Coleman-Weinberg potential at the
   * tree-level minimum.
   */
  std::vector<double> WeinbergForthDerivative() const;
  /**
   * Calculates the Debye corrections to the Higgs mass matrix.
   * If you can provide CalculateDebyeSimplified() with the Matrix as this will
   * reduce the runtime.
   */
  void CalculateDebye();
  /**
   * Calculates the Debye corrections to the gauge sector. By using
   * CalculateDebyeGaugeSimplified() the runtime can be reduced.
   */
  void CalculateDebyeGauge();
  /**
   * Sets a tensor needed to calculate the contribution of the counterterm
   * potential to the triple Higgs couplings.
   */
  void Prepare_Triple();

  /**
   * You can give the explicit Debye corrections to the Higgs mass matrix with
   * this function if you know it. Otherwise it is calculated via
   * CalculateDebye(). If you know the corrections use this and let the function
   * return true, this will save you a lot of computing time.
   */
  virtual bool CalculateDebyeSimplified() = 0;

  /**
   * You can give the explicit Debye corrections to the gauge boson mass matrix
   * with this function if you know it. Otherwise it is calculated via
   * CalculateDebyeGauge(). If you know the corrections use this and let the
   * function return true, this will save you a lot of computing time.
   */
  virtual bool CalculateDebyeGaugeSimplified() = 0;

  /**
   * You can give the explicit form of your tree-level potential here. This
   * speeds up the computation time.
   */
  virtual double VTreeSimplified(const std::vector<double> &v) const = 0;
  /**
   * @brief UseVTreeSimplified Decides wether VTreeSimplified will be used or
   * not. VTreeSimplified returns 0 if UseVTreeSimplified is false Set in
   * constructor of the implemented models
   */
  bool UseVTreeSimplified = false;
  /**
   * You can give the explicit form of your counterterm potential here. This
   * speeds up the computation time.
   */
  virtual double VCounterSimplified(const std::vector<double> &v) const = 0;
  /**
   * @brief UseVCounterSimplified Decides wether VCounterSimplified will be used
   * or not. VCounterSimplified returns 0 if UseVCounterSimplified is false Set
   * in constructor of the implemented models
   */
  bool UseVCounterSimplified = false;

  /**
   * Calculates the Higgs mass matrix and saves all eigenvalues
   * @param v the configuration of all VEVs at which the eigenvalues should be
   * evaluated
   * @param Temp The temperature at which the Debye corrected masses should be
   * calculated
   * @param diff 0 returns the masses and i!=0 returns the derivative of m^2
   * w.r.t v_i, i = -1 returns the derivative w.r.t. to the temperature
   * @return Vector in which the eigenvalues m^2 of the mass matrix will be
   * stored
   */
  std::vector<double> HiggsMassesSquared(const std::vector<double> &v,
                                         const double &Temp = 0,
                                         const int &diff    = 0) const;

  /**
   * @brief HiggsMassMatrix calculates the Higgs mass matrix
   * @param v the configuration of all VEVs at which the Mass Matrix should be
   * evaluated
   * @param Temp The temperature at which the Debye corrected masses should be
   * calculated
   * @param diff 0 returns the masses and i!=0 returns the derivative the Mass
   * Matrix w.r.t v_i, i = -1 returns the derivative w.r.t. to the temperature
   * @return
   */
  Eigen::MatrixXd HiggsMassMatrix(const std::vector<double> &v,
                                  double Temp = 0,
                                  int diff    = 0) const;

  /**
   * Calculates the gauge mass matrix and saves all eigenvalues
   * @param v the configuration of all VEVs at which the eigenvalues should be
   * evaluated
   * @param Temp The temperature at which the Debye corrected masses should be
   * calculated
   * @param diff 0 returns the masses and i!=0 returns the derivative of m^2
   * w.r.t v_i, -1 returns the derivative w.r.t. the temperature
   * @return Vector in which the eigenvalues m^2 of the mass matrix will be
   * stored
   */
  std::vector<double> GaugeMassesSquared(const std::vector<double> &v,
                                         const double &Temp = 0,
                                         const int &diff    = 0) const;
  /**
   * Calculates the quark mass matrix and saves all eigenvalues, this assumes
   * the same masses for different colours.
   * @param v the configuration of all VEVs at which the eigenvalues should be
   * evaluated
   * @param diff 0 returns the masses and i!=0 returns the derivative of m^2
   * w.r.t v_i
   * @return Vector in which the eigenvalues m^2 of the mass matrix will be
   * stored
   */
  std::vector<double> QuarkMassesSquared(const std::vector<double> &v,
                                         const int &diff = 0) const;
  /**
   * Calculates the lepton mass matrix and saves all eigenvalues
   * @param v the configuration of all VEVs at which the eigenvalues should be
   * evaluated
   * @param diff 0 returns the masses and i!=0 returns the derivative of m^2
   * w.r.t v_i
   * @return Vector in which the eigenvalues m^2 of the mass matrix will be
   * stored
   */
  std::vector<double> LeptonMassesSquared(const std::vector<double> &v,
                                          const int &diff = 0) const;

  /**
   * Calculates the quark mass matrix and saves all eigenvalues, this assumes
   * the same masses for different colours.
   * @param v the configuration of all VEVs at which the eigenvalues should be
   * evaluated
   * @return Vector in which the complex eigenvalues m of the mass matrix will
   * be stored
   */
  std::vector<std::complex<double>>
  QuarkMasses(const std::vector<double> &v) const;
  /**
   * @brief QuarkMassMatrix calculates the Mass Matrix for the Quarks of the
   * form $ M^{IJ} = Y^{IJ} + Y^{IJk} v_k $
   * @param v the configuration of all VEVs at which the matrix should be
   * calculated
   * @return the Mass Matrix for the Quarks of the form $ M^{IJ} = Y^{IJ} +
   * Y^{IJk} v_k $
   */
  Eigen::MatrixXcd QuarkMassMatrix(const std::vector<double> &v) const;
  /**
   * Calculates the quark mass matrix and saves all eigenvalues, this assumes
   * the same masses for different colours.
   * @param v the configuration of all VEVs at which the eigenvalues should be
   * evaluated
   * @return Vector in which the complex eigenvalues m of the mass matrix will
   * be stored
   */
  std::vector<std::complex<double>>
  LeptonMasses(const std::vector<double> &v) const;
  /**
   * @brief LeptonMassMatrix calculates the Mass Matrix for the Leptons of the
   * form $ M^{IJ} = Y^{IJ} + Y^{IJk} v_k $
   * @param v the configuration of all VEVs at which the matrix should be
   * calculated
   * @return the Mass Matrix for the Leptons of the form $ M^{IJ} = Y^{IJ} +
   * Y^{IJk} v_k $
   */
  Eigen::MatrixXcd LeptonMassMatrix(const std::vector<double> &v) const;

  /**
   * Calculates the triple Higgs couplings at NLO in the mass basis.
   *
   * Use the vector TripleHiggsCorrectionsCWPhysical to save your couplings and
   * set the nTripleCouplings to the number of couplings you want as output.
   */
  virtual void TripleHiggsCouplings() = 0;
  /**
   * Calculates the function f^1 needed for the derivatives of the Coleman
   * Weinberg potential.
   */
  double fbase(double MassSquaredA, double MassSquaredB) const;
  /**
   * Calculates the function f^2 needed for the 3rd derivatives of the Coleman
   * Weinberg potential.
   */
  double
  fbaseTri(double MassSquaredA, double MassSquaredB, double MassSquaredC) const;
  /**
   * Calculates the function f needed for the 4th derivatives of the Coleman
   * Weinberg potential.
   */
  double fbaseFour(double MassSquaredA,
                   double MassSquaredB,
                   double MassSquaredC,
                   double MassSquaredD) const;

  /**
   * Calculates the counterterm parameters. Here you need to work out the scheme
   * and implement the formulas. This has to be specified in the model file.
   */
  virtual std::vector<double> calc_CT() const = 0;

  /**
   * Calculates the Coleman-Weinberg contribution of a particle with m^2 =
   * MassSquared and the constant scheme-dependent parameter cb as well as the
   * first derivative with respect to m^2.
   */
  double CWTerm(double MassSquared, double cb, int diff = 0) const;
  /**
   * Calculates Re(log(MassSquared)) and returns 0 if the argument is too small
   * as this function is only called with an (m^2)^n in front of it.
   */
  double FCW(double MassSquared) const;
  /**
   * @brief Calculation of the bosonic std::size_tegral + Coleman-Weinberg
   * potential __without__ taking d.o.f. std::size_to account
   *
   * @param MassSquared = m^2 of the particle
   * @param Temp = Temperature at which the Debye masses and std::size_tegrals
   * should be calculated
   * @param cb = Parameter of the renormalisation-Scheme in the Coleman-Weinberg
   * potential
   * @param diff: 0 returns the value of the std::size_tegral and diff >0 the
   * derivative w.r.t. m^2 and diff = -1 w.r.t. Temp
   *
   *
   */
  double boson(double MassSquared, double Temp, double cb, int diff = 0) const;
  /**
   * Deprecated version of boson() as present in the v1.X release. Still here
   * for legacy reasons
   */
  [[deprecated(
      "Use this only if you want to calculate the effective potential with the "
      "exact same routine used in BSMPT v1.x. Use boson otherwise.")]] double
  boson_legacy(double MassSquared, double Temp, double cb, int diff = 0) const;
  /**
   * @brief Calculation of the fermionic std::size_tegral + Coleman-Weinberg
   * potential __without__ taking d.o.f. std::size_to account
   *
   * @param MassSquared = m^2 of the particle
   * @param Temp = temperature at which the std::size_tegrals should be
   * calculated
   * @param diff :  0 = Value of the potential, diff > 0 returns the derivative
   * w.r.t. m^2 and diff=-1 w.r.t Temp
   */
  double fermion(double MassSquared, double Temp, int diff = 0) const;
  /**
   * Deprecated version of fermion() as present in the v1.X release. Still here
   * for legacy reasons
   */
  [[deprecated(
      "Use this only if you want to calculate the effective potential with the "
      "exact same routine used in BSMPT v1.x. Use fermion otherwise.")]] double
  fermion_legacy(double Mass, double Temp, int diff = 0) const;
  /**
   * Calculates the large m^2/T^2 approximation to order n for the
   * temperature-dependent std::size_tegrals.
   */
  [[deprecated("Use this only if you want to calculate the effective potential "
               "with the exact same routine used in BSMPT v1.x. "
               "Otherwise use ThermalFunctions::JInterpolatedHigh.")]] double
  Vl(double MassSquared, double Temp, int n, int diff = 0) const;
  /**
   * Calculates the small m^2/T^2 approximation of the fermionic temperature
   * dependent std::size_tegral to the n-th order.
   * @param MassSquared m^2
   * @param Temp temperature
   * @param n the order to which expand the expansion
   * @param diff 0 = value; >0 derivative w.r.t.
   */
  [[deprecated(
      "Use this only if you want to calculate the effective potential with the "
      "exact same routine used in BSMPT v1.x. "
      "Otherwise use ThermalFunctions::JfermionInterpolatedLow.")]] double
  Vsf(double MassSquared, double Temp, int n, int diff = 0) const;

  /**
   * Calculates the small m^2/T^2 approximation of the bosonic temperature
   * dependent std::size_tegral to the n-th order.
   */
  [[deprecated(
      "Use this only if you want to calculate the effective potential with the "
      "exact same routine used in BSMPT v1.x. "
      "Otherwise use ThermalFunctions::JbosonInterpolatedLow.")]] double
  Vsb(double MassSquaredIn, double Temp, int n, int diff = 0) const;

  /**
   * Resets all bools. Needed if you want to deal with multiple posize_ts one
   * after another with the same pointer.
   */
  void resetbools();

  /**
   * This will convert the vector VEVminimizer with the VEVs used by the
   * minimizer to a vector with the complete VEV configuration VEVFunction used
   * by the Class functions. The mapping is done with the VevOrder vector which
   * has to be set in the constructor of your model
   * @param VEVminimizer vector with the VEVs, used as inputs for the minimizers
   * @return vector with all VEVs, filled with zeroes for fields which do not
   * have a VEV
   */
  std::vector<double>
  MinimizeOrderVEV(const std::vector<double> &VEVminimizer) const;

  /**
   * Calculates the first derivatives of the eigenvalues of a given matrix
   * @param M : the original matrix
   * @param MDiff : the element-wise first derivative of the matrix M with
   * respect to the parameter you want to consider
   * @return vector containing the mass eigenvalues and then the derivatives in
   * the same order
   */
  std::vector<double>
  FirstDerivativeOfEigenvalues(const Eigen::Ref<Eigen::MatrixXcd> M,
                               const Eigen::Ref<Eigen::MatrixXcd> MDiff) const;
  /**
   * This function calculates the second derivatives of all eigenvalues.
   * The matrix must not have a repeated eigenvalue for this!
   * @param M : the original matrix
   * @param MDiffX : the element-wise first derivative of the matrix M with
   * respect to the first parameter you want to consider
   * @param MDiffY : the element-wise first derivative of the matrix M with
   * respect to the second parameter you want to consider
   * @param MDiffXY : the element-wise second derivative of the matrix M with
   * respect to both parameters you want to consider
   * @return the mass eigenvalues in the vector and then the derivatives in the
   * same order
   */
  std::vector<double> SecondDerivativeOfEigenvaluesNonRepeated(
      const Eigen::Ref<Eigen::MatrixXd> M,
      const Eigen::Ref<Eigen::MatrixXd> MDiffX,
      const Eigen::Ref<Eigen::MatrixXd> MDiffY,
      const Eigen::Ref<Eigen::MatrixXd> MDiffXY) const;

  /**
   * This function will check if the VEV at NLO is still close enough to the LO
   * VEV. v has to be given in the configuration of the minimizer.
   */
  bool CheckNLOVEV(const std::vector<double> &v) const;

  /**
   * This is a possible debugging function.
   */
  virtual void Debugging(const std::vector<double> &input,
                         std::vector<double> &output) const = 0;

  /**
   * Checks if the tensors are correctly implemented. For this the fermion,
   * quark and gauge boson masses are calculated and printed next to the values
   * defined in SMparah.h
   */
  void CheckImplementation(
      const int &WhichMinimizer = Minimizer::WhichMinimizerDefault) const;

  /**
   * Find all possible sign combinations of the vevs under which the potential
   * is invariant
   */
  std::vector<std::vector<double>> SignSymmetries;
  /**
   * @brief FindSignSymmetries checks for all possible sign changes in the VEV
   * components and checks for all possible Z2 symmetries
   */
  void FindSignSymmetries();

  /**
   * Set the parameter UseTreeLevel to the input
   */
  void SetUseTreeLevel(bool val);

  /**
   * Gets the parameter line as an Input and calls
   * resetbools, ReadAndSet, calc_CT, set_CT_Pot_Par,CalculateDebye and
   * CalculateDebyeGauge
   */
  std::pair<std::vector<double>, std::vector<double>>
  initModel(std::string linestr);

  /**
   * Gets the parameter in a vector in the same way as in set_gen and calls
   * resetbools, set_gen, calc_CT, set_CT_Pot_Par,CalculateDebye and
   * CalculateDebyeGauge
   * @param par Parameters to define the parameter point. Same input as in
   * set_gen()
   * @return Vector with the CT
   */
  std::vector<double> initModel(const std::vector<double> &par);

  /**
   * @brief resetScale changes the MSBar scale to newScale
   * @param newScale the new Scale in GeV
   * @return the CT at the new scale
   */
  std::vector<double> resetScale(const double &newScale);

  /**
   * Calculate the ratio of the latent heat w.r.t. the energy density of the
   * plasma background Formula taken from
   */

  double CalculateRatioAlpha(const std::vector<double> &vev_symmetric,
                             const std::vector<double> &vev_broken,
                             const double &Temp) const;

  /**
   * @brief NablaVCT
   * @return
   */
  Eigen::VectorXd NablaVCT(const std::vector<double> &v) const;
  /**
   * @brief HessianWeinberg
   * @return
   */
  Eigen::MatrixXd HessianCT(const std::vector<double> &v) const;

  /**
   * @brief GetCTIdentities
   * @return vector with the identities required to vanish for the CT
   * calculations. Returns an empty vector if not set otherwise.
   */
  virtual std::vector<double> GetCTIdentities() const;
};

} // namespace BSMPT
#endif /* CLASSPOTENTIALORIGIN_H_ */
