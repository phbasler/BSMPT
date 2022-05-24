<!--
SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller

SPDX-License-Identifier: GPL-3.0-or-later
-->

Program: BSMPT version 2.4.0

Released by: Philipp Basler and Lisa Biermann and Margarete Mühlleitner and Jonas Müller

[!["GitHub Discussions"](https://img.shields.io/badge/%20GitHub-%20Discussions-gray.svg?longCache=true&logo=github&colorB=purple)](https://github.com/phbasler/BSMPT/discussions)
[![Ubuntu unit tests](https://github.com/phbasler/BSMPT/actions/workflows/test.yml/badge.svg?branch=master)](https://github.com/phbasler/BSMPT/actions/workflows/test.yml)
[![Mac unit tests](https://github.com/phbasler/BSMPT/actions/workflows/test-mac.yml/badge.svg?branch=master)](https://github.com/phbasler/BSMPT/actions/workflows/test-mac.yml)
[![Windows unit tests](https://github.com/phbasler/BSMPT/actions/workflows/windows_unit_tests.yml/badge.svg?branch=master)](https://github.com/phbasler/BSMPT/actions/workflows/windows_unit_tests.yml)
[![codecov master](https://codecov.io/gh/phbasler/BSMPT/branch/master/graph/badge.svg?token=LDGNQTADB5)](https://codecov.io/gh/phbasler/BSMPT)
[![Documentation](https://img.shields.io/badge/Documentation-master-success)][DoxygenLink]
[![Benchmarks](https://img.shields.io/badge/Benchmark-master-success)](https://phbasler.github.io/BSMPT/benchmarks/)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/phbasler/bsmpt/graphs/commit-activity)
[![GitHub license](https://img.shields.io/github/license/phbasler/bsmpt.svg)](https://github.com/phbasler/BSMPT/blob/master/LICENSE.md)
[![Latest release](https://badgen.net/github/release/phbasler/bsmpt)](https://github.com/phbasler/bsmpt/releases)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/phbasler/BSMPT.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/phbasler/BSMPT/context:cpp)



Manual: version 2.0

BSMPT - Beyond the Standard Model Phase Transitions:
The C++ program package BSMPT calculates the strength of the electroweak phase transition in extended Higgs sectors. For this the loop-corrected effective potential at finite temperature is calculated including the daisy resummation of the bosonic masses.
The program computes the vacuum expectation value (VEV) \f$ v \f$ of the potential
as a function of the temperature, and in particular the critical VEV
\f$v_c\f$ at the temperature \f$T_c\f$ where the phase transition takes
place. 
In addition, the loop-corrected trilinear Higgs self-couplings are
provided. We apply an 'on-shell' renormalization scheme in the sense
that the loop-corrected masses and mixing angles are required to be
equal to their tree-level input values. This allows for efficient
scans in the parameter space of the models. 

The models implemented so far are

  - CP-conserving 2-Higgs-Doublet Models (R2HDM)
  - CP-violating  2-Higgs-Doublet Models (C2HDM)
  - Next-to-Minimal 2HDM (N2HDM)
  - CP in the Dark ([arXiv 1807.10322](https://arxiv.org/abs/1807.10322), [arXiv 2204.13425](https://arxiv.org/abs/2204.13425))
  - Complex Singlet Extension (CxSM)

The code is structured such that users can add their own models.



The program package can be downloaded at:
https://github.com/phbasler/BSMPT

The documentation of the code is provided at [https://phbasler.github.io/BSMPT/documentation][DoxygenLink].

Sample input and output files are provided in the directory 'example'.

Modifications and corrected bugs are reported in the file 'Changelog.md'.



For additional information, comments, complaints or suggestions please e-mail
to:  bsmpt@lists.kit.edu, open a corresponding [issue](https://github.com/phbasler/BSMPT/issues) or start a [discussion](https://github.com/phbasler/BSMPT/discussions).

---

### Citation:
If you use this program for your work, please cite 

  - [1803.02846](https://arxiv.org/abs/1803.02846)
  - [2007.01725](https://arxiv.org/abs/2007.01725)

## Installation:

### Dependencies
BSMPT uses cmake which will try to find the following libraries. If they are not found automatically, you can use the -DpackageName_ROOT option of cmake to tell it where to look for the package.

1. [GSL library](https://www.gnu.org/software/gsl/). 
2. [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page)  
3. [libcmaes](https://github.com/CMA-ES/libcmaes): Additionally to GSL you should either use libcmaes or NLopt. If libcmaes is installed through cmake BSMPT should find it automatically, otherwise you can point it to the install direction with
    `-Dlibcmaes_ROOT=/path/to/cmaes`  
    
    If cmaes is not installed then it will be installed in your build directory. For more details on the libcmaes installation, e.g. possible dependencies, visit their [wiki](https://github.com/CMA-ES/libcmaes/wiki). If you don't want to install or use it, you can set `-DUseLibCMAES=OFF` 
    
4. [NLopt](https://nlopt.readthedocs.io/en/latest/): If you do not want to use NLopt, you can set `-DUseNLopt=OFF`
5. [Boost](https://www.boost.org/) : This is optional and only required for the Baryogenesis calculations. If it is not found then these parts are not available.

### Alternative Install method 
If you have [conan](https://conan.io/) installed, then you can set the `-DUseConan=ON` flag for cmake and it will download boost, gsl, eigen3 and NLopt (if UseNLopt was not turned off) through [conancenter](https://conan.io/center/).

### build
With the dependencies and options you can build the programm with
  
        mkdir build && cd build  
        cmake (OPTIONS from Dependencies) ..  
        cmake --build . -j  
        cmake --build . -j -t doc
    

The `-t doc` will use doxygen to create the online help in build/html which can be opened locally.


### Unit tests
After compiling the code call `ctest -j` in the build folder to run some checks. Here the NLO VEV and EWPT for the R2HDM, C2HDM and N2HDM example points will be calculated and compared to the expected results. 

---

## How to add a new model (for further details, also see the manual):

To add a new model you have to modify/create five files  

1. Go to include/BSMPT/models and copy ClassTemplate.h to YourModel.h. Adjust the Class_Template name to your new model. For step 5 I will assume that your class is named Class_YourModel. 

2. Go to src/models and copy ClassTemplate.cpp to YourModel.cpp, and change the Class_Template class name in the file to your model name. Also follow the instructions in here and in the manual to set up your new model. 

3. For your model to compile you have to open src/models/CMakeLists.txt and add ${header_path}/YourModel.h in the set(header enviroment as well as YourModel.cpp in the set(src enviroment)

4. In include/BSMPT/models/IncludeAllModels.h you need to add a new entry in the ModelIDs enum above the `stop` entry which is different from the ones already in the enum, e.g. YourModel. Additionally, you have to create a new entry in the `const std::unordered_map<std::string,ModelIDs> ModelNames` map in the same file and add a new line with {"YourModelName",ModelIDs::YourModel} , the matching will be done automatically.
Then you can call your model with `./binary YourModelName ...` .

5. In src/models/IncludeAllModels.cpp you have to add `#include <BSMPT/models/YourModel.h>` to the include list. Also to actually call your model you have to extend the FChoose function. For this you add a new case to the switch statement, which reads

        case ModelIDs::YourModel: return std::make_unique<Class_YourModel>(); break;



You can use the Test executable to detect possible errors in your implementation. If the Test executable does not show you an error, but something is still wrong, contact us at bsmpt@lists.kit.edu

Also contact us if you have a custom model for BSMPT v1.x and you have trouble converting it to the new notation.

## Executables
BSMPT provides multiple executables. Here we give a quick overview of them. For every executable you can call them with the `--help` option to get an overview of possible input parameters.
For the following examples the C2HDM with the example/C2HDM_Input.dat file is used.

### BSMPT
BSMPT calculates the EWPT for every parameter point in the input file and gives out the results of those parameter points for which \f$v_c/T_c\f$ > 1. 
To find these points a bisection method is used with the temperature starting between 0 and 300 GeV. The executable is called through:

        ./bin/BSMPT Model Inputfile Outputfile LineStart LineEnd

This will call the specific model to be used, identified through 'Model', and calculate the EWPT for each parameter point (corresponding to one line) between 'LineStart' and 'LineEnd'.

For our example the command

    	./bin/BSMPT c2hdm example/C2HDM_Input.dat example/test_BSMPT.dat 2 2

will calculate for the C2HDM the EWPT for one parameter point given in line 2 in C2HDM_Input.dat. This will generate the output file example/test_BSMPT.dat 2 2 which can be compared with the already available file example/C2HDM_Input.dat_BSMPT.

### CalcCT
This will calculate the counterterms. In the output file the information on the input parameter point is given and the counterterms are added at the end of the line.

It is called through the command

    	./bin/CalcCT Model Inputfile Outputfile LineStart LineEnd

For the C2HDM example this reads

    	./bin/CalcCT c2hdm example/C2HDM_Input.dat example/test_CalcCT.dat 2 2

which will generate the output file example/test_CalcCT.dat. This can be compared with the already available file example/C2HDM_Input.dat_CalcCT.

### NLOVEV
This calculates the VEV at 1-loop order at vanishing temperature in the effective potential approach. This can be used to investigate the vacuum stability of the model. It is called through

    	./bin/NLOVEV Model Inputfile Outputfile LineStart LineEnd

and for the C2HDM example it is given by

    	./bin/NLOVEV c2hdm example/C2HDM_Input.dat example/test_NLOVEV.dat 2 2

where the result is written into the file example/test_NLOVEV.dat which can be compared with the already available file example/C2HDM_Input.dat_NLOVEV.

### VEVEVO
This program calculates the evolution of the vacuum expecation value of a given point with the temperature. It is called through

    	./bin/VEVEVO Model Inputfile Outputfile Line Tempstart Tempstep Tempend

where 'Tempstart' is the starting value of the temperature which increases with 'Tempstep' until 'Tempend'.

For our C2HDM example this would be

     	./bin/VEVEVO c2hdm example/C2HDM_Input.dat example/test_VEVEVO.dat 2 0 5 150

where the result for the NLO VEV is given in example/test_VEVEVO.dat as function of the temperature in the interval between 0 and 150 GeV in steps of 5 GeV. This can be compared with the already available file example/C2HDM_Input.dat_vevevo.

### TripleHiggsCouplingNLO
This program calculates the trilinear Higgs self-couplings at NLO at zero temperature. It is called through

    	./bin/TripleHiggsNLO Model Inputfile Outputfile LineStart LineEnd

The C2HDM example is called through

    	./bin/TripleHiggsNLO c2hdm example/C2HDM_Input.dat example/test_TripleHiggsCouplingNLO.dat 2 2

with the result given in example/test_TripleHiggsNLO.dat which can be compared with the already available file example/C2HDM_Input.dat_TripleHiggsCouplingNLO .

### CalculateEWBG
This program calculates the difference between baryons and anti-baryons
normalized to the photon density generated through the EWPT.
Please beware that this is only tested for the C2HDM so far and the general
implementation is future work. It is called through

    	./bin/CalculateEWBG c2hdm Inputfile Outputfile LineStart LineEnd config_file

 An example is given for the example/C2HDM_Input.dat parameter point through

    	./bin/CalculateEWBG c2hdm example/C2HDM_Input.dat example/test_EWBG.dat 2 2 example/EWBG_config.txt

with the result given in example/test_EWBG.dat which can be compared with
the already available file example/C2HDM_Input.dat_EWBG.

### PlotEWBG_vw
This executable varies the wall velocity of a given parameter point and
calculates the EWBG for each velocity.

### PlotEWBG_nL
This executable calculates the left-handed fermion density in front of the
wall as a function of the distance to the bubble wall.

### Generate the C++ code for a model
We provide currently two methods to generate the tensors and calculate the counter terms for a new model.

1. At tools/ModelGeneration/Maple we provide the maple Worksheet CreateModel.mw which you can use to implement your model and get the tensors. 
2. At tools/ModelGeneration/sympy we provide a setup using only python3 with sympy (at least version 1.10!, if your packet manager only has an older installed, e.g. ubuntu 20.04 only has v1.6, then you have to install v1.10 or up with pip). Here we provide two examples, SM.py and G2HDM.py which both implement two different models and use the ModelGenerator.py module to calculate the tensors and CT. You can get the CT using `python3 SM.py --show ct` and the tensors by calling `python3 SM.py --show tensors`. If your counterterms don't have a unique solution, then the solution space will be shown to you and you have to add additional equations until you have a unique solution (e.g. the G2HDM example).

[DoxygenLink]: https://phbasler.github.io/BSMPT/documentation
