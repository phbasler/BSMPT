<!--
SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller

SPDX-License-Identifier: GPL-3.0-or-later
-->

Program: BSMPT version 2.3.3

Released by: Philipp Basler and Margarete Mühlleitner and Jonas Müller

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
  - CP in the dark ([arXiv 1807.10322](https://arxiv.org/abs/1807.10322))
  - Complex Singlet Extension (CxSM)

The code is structured such that users can add their own models.



The program package can be downloaded at:
https://github.com/phbasler/BSMPT

A short documentation is given at: [DoxygenLink]

Sample input and output files are provided in the directory 'example'.

Modifications and corrected bugs are reported in the file 'Changelog.md'.



For additional information, comments, complaints or suggestions please e-mail
to:  bsmpt@lists.kit.edu or create a corresponding Issue.

---

### Citation:
If you use this program for your work, please cite [1803.02846](https://arxiv.org/abs/1803.02846) and [2007.01725](https://arxiv.org/abs/2007.01725)

##Installation:

### Dependencies
1. GSL library: The code assume GSL is installed in PATH
2. Eigen3: If eigen3 is not found automatically, you need to give the path to your eigen3 installation.  To install eigen3 you go to your downloaded eigen3 folder and install it through   
  
        mkdir build && cd build  
        cmake -DCMAKE_INSTALL_PREFIX=/path/to/installedEigen  ..  
        make install  
After that you can use `-DEigen3_Dir=/path/to/installedEigen/share/eigen3/cmake` to link Eigen3
  
3. libcmaes: Additionally to GSL you should either use libcmaes or NLopt. If libcmaes is installed through cmake BSMPT should find it automatically, otherwise you can point it to the install direction with
    `-Dlibcmaes_ROOT=/path/to/cmaes`  
    
    If cmaes is not installed then it will be installed in your build directory. For more details on the libcmaes installation, e.g. possible dependencies, visit their [wiki](https://github.com/CMA-ES/libcmaes/wiki). If you don't want to install or use it, you can set `-DUseLibCMAES=OFF` 
    
4. [NLopt](https://nlopt.readthedocs.io/en/latest/): If NLopt is installed through your packet manager cmake will find it automatically. Otherwise, with `-DNLopt_DIR=/path/to/installedNLopt/lib/cmake/nlopt` you can tell where NLopt is installed. If you do not want to use NLopt, you can set `-DUseNLopt=OFF`
5. [Boost](https://www.boost.org/) : It should be found automatically, if not you can use `-DBOOST_ROOT=/path/to/boost`

### Alternative Install method 
If you have [conan](https://conan.io/) installed, then you can set the `-DUseConan=ON` flag for cmake and it will download boost, gsl, eigen3 and NLopt (if UseNLopt was not turned off) through [conancenter](https://conan.io/center/).

### build
With the dependencies and options you can build the programm with
  
        mkdir build && cd build  
        cmake (OPTIONS from Dependencies) ..  
        cmake --build . -j  
        cmake --build . -j -t doc
    

The make doc will use doxygen to create the online help in build/html which can be opened locally.


Note to Mac Users: You have to use the g++ compiler as clang does not support OpenMP. If you get the error "missing libcmaes_config.h" during the compilation, please check if the file is in build/libcmaes-0.9.5. If so, copy it to /path/to/libcmaes/include/libcmaes


### Unit tests
After calling `make` in the build directory you can call `ctest`or `./bin/GenericTests` to run some checks. Here the NLO VEV and EWPT for the R2HDM, C2HDM and N2HDM example points will be calculated and compared to the expected results. 

---

##How to add a new model (for further details, also see the manual):

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

[DoxygenLink]: https://phbasler.github.io/BSMPT/documentation