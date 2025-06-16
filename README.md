<!--
SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller

SPDX-License-Identifier: GPL-3.0-or-later
-->

Program: BSMPT version 3.1.2

Released by: Philipp Basler, Lisa Biermann, Margarete Mühlleitner, Jonas Müller, Rui Santos and João Viana

[![GitHub Discussions](https://img.shields.io/badge/%20GitHub-%20Discussions-gray.svg?longCache=true&logo=github&colorB=purple)](https://github.com/phbasler/BSMPT/discussions)
[![Unit tests](https://github.com/phbasler/BSMPT/actions/workflows/test.yml/badge.svg?branch=master)](https://github.com/phbasler/BSMPT/actions/workflows/test.yml)
[![codecov master](https://codecov.io/gh/phbasler/BSMPT/branch/master/graph/badge.svg?token=LDGNQTADB5)](https://codecov.io/gh/phbasler/BSMPT)
[![Documentation](https://img.shields.io/badge/Documentation-master-success)][DoxygenLink]
[![Benchmarks](https://img.shields.io/badge/Benchmark-master-success)](https://phbasler.github.io/BSMPT/benchmarks/)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/phbasler/bsmpt/graphs/commit-activity)
[![GitHub license](https://img.shields.io/github/license/phbasler/bsmpt.svg)](https://github.com/phbasler/BSMPT/blob/master/LICENSE.md)
[![Latest release](https://badgen.net/github/release/phbasler/bsmpt)](https://github.com/phbasler/bsmpt/releases)


Manual: version 3.0

BSMPT - Beyond the Standard Model Phase Transitions:

The C++ program package BSMPT allows for the detailed study of (multi-step) phase transitions between temperature-dependent minima in the one-loop daisy-resummed finite-temperature effective potential.

The program tracks temperature-dependent minima, calculates the bounce solution, the characteristic temperatures and gravitational wave signals of first-order phase transitions.
The code also allows to derive the loop-corrected trilinear Higgs self-couplings and provides the computation of the baryon asymmetry for the CP-violating 2-Higgs Doublet Model (C2HDM).

We apply an 'on-shell' renormalization scheme in the sense
that the loop-corrected masses and mixing angles are required to be
equal to their tree-level input values. This allows for efficient
scans in the parameter space of the models. 

The models implemented so far are

  - Standard Model (SM)
  - CP-conserving 2-Higgs-Doublet Model (R2HDM)
  - CP-violating 2-Higgs-Doublet Model (C2HDM)
  - Next-to-Minimal 2HDM (N2HDM)
  - CP in the Dark ([arXiv 1807.10322](https://arxiv.org/abs/1807.10322), [arXiv 2204.13425](https://arxiv.org/abs/2204.13425))
  - Complex Singlet Extension (CxSM)

The code is structured such that users can add their own models.



The program package can be downloaded at:
https://github.com/phbasler/BSMPT

The documentation of the code is provided at [https://phbasler.github.io/BSMPT/documentation][DoxygenLink].

Sample input and output files are provided in the directory 'example'.

Modifications and corrected bugs are reported in the file 'Changelog.md'.



For additional information, comments, complaints or suggestions open a corresponding [issue](https://github.com/phbasler/BSMPT/issues) or start a [discussion](https://github.com/phbasler/BSMPT/discussions). For non-public matters please send an e-mail to bsmpt@lists.kit.edu.

---

### Citation:
If you use this program for your work, please cite 

  - [1803.02846](https://arxiv.org/abs/1803.02846)
  - [2007.01725](https://arxiv.org/abs/2007.01725)
  - [2404.19037](https://arxiv.org/abs/2404.19037)

## Installation:
BSMPT uses cmake and [conan 2](https://conan.io/) for its installation which can be installed through pip with `pip3 install cmake conan`.
In addition, you need a `C` and `C++` compiler installed.

### build - simple
If you want the default installation of BSMPT, you can then use the `Build.py` script.
The script `Build.py` installs the necessary conan profiles for your operating system, handles the dependencies and compiles `BSMPT` with its default settings in release mode. You can execute it with

```bash
python3 Build.py
```

### Dependencies
BSMPT uses cmake and [conan 2](https://conan.io/) for dependencies. The used dependencies are:

1. [GSL library](https://www.gnu.org/software/gsl/). 
2. [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page)  
3. [libcmaes](https://github.com/CMA-ES/libcmaes): Additionally to GSL you should either use libcmaes or NLopt. For more details on the libcmaes installation, e.g. possible dependencies, visit their [wiki](https://github.com/CMA-ES/libcmaes/wiki). If you don't want to install or use it, you can set `--options UseLibCMAES=False` when using the detailed build, as described below.
4. [NLopt](https://nlopt.readthedocs.io/en/latest/): If you do not want to use NLopt, you can set `--options UseNLopt=False` when using the detailed build, as described below.
5. [Boost](https://www.boost.org/) : This is optional and only required for the Baryogenesis calculations. In order to compile the Baryogenesis calculation, set `--options CompileBaryo=True` when using the detailed build, as described below.


### build - detailed
We provide the script `Setup.py` which installs conan profiles for your operating system and runs `conan install` to download the dependencies. If you want to use other profiles feel free to execute `conan install` with your profile manually or add it to the script.

You can build the code with

```bash
python3 Setup.py  
cmake --preset ${profile}  
cmake --build --preset ${profile} -j  
cmake --build --preset ${profile} -j -t doc
```
    

The `-t doc` will use doxygen to create the online help in build/html which can be opened locally.
The `${profile}` parameter depends on your operating system. After running the `Setup` script you can call `cmake --list-presets` to show the found presets.

The script `Setup.py` can take several optional arguments, run `python3 Setup.py -h` or `python3 Setup.py --help` to display them.


### Unit tests
After compiling the code call `ctest --preset ${profile} -j` in the root folder to run some checks. Here the NLO VEV and EWPT for the R2HDM, C2HDM and N2HDM example points will be calculated and compared to the expected results. 


### Development
Most modern IDEs support cmake profiles. After running the `Setup.py` script you can open the root folder in an IDE of your choice (e.g. VSCode with cmake extension) and it will recognise the cmake profile.


Code from the following repositories is used in BSMPT:
- [`ttk592::spline`](https://github.com/ttk592/spline) by Tino Kluge, a C++ cubic spline interpolation library 
- [`AsciiPlotter`](https://github.com/joehood/asciiplotter) by Joe Hood, for ASCII plots in the terminal 

---

## How to include BSMPT in another program

You can call the `Setup.py` script with the `--create` option. This will generate a local conan package of BSMPT which then can be used in your program.
If your code is already a conan project you can add BSMPT just as a requirement. If you have a pure CMake project you can use the [conan-provider](https://github.com/conan-io/cmake-conan/) to include it.

---

## How to add a new model:

To add a new model, you have to modify/create five files (for further details, also consult the manual):

1. Go to `include/BSMPT/models` and copy `ClassTemplate.h` to `YourModel.h`. Adjust the name of the class `Class_Template` to `Class_YourModel`.

2. Go to `src/models` and copy `ClassTemplate.cpp` to `YourModel.cpp`, and again change `Class_Template` to `Class_YourModel`. Also, follow the instructions in this file and in the manual to set up your new model. 

3. For your model to compile, you have to open `src/models/CMakeLists.txt` and add `${header_path}/YourModel.h` as well as `YourModel.cpp` to the listed headers and source files.

4. In `include/BSMPT/utility/ModelIDs.h` you need to add a new entry in the `enum class ModelIDs` above the `stop` entry which is different from the already defined `ModelIDs`, e.g. `YourModel`. Additionally, you have to create a new entry in the `const std::unordered_map<std::string, ModelIDs> ModelNames` map in the same file and add a new line with `{"YourModelName",ModelIDs::YourModel}`.

5. In `src/models/IncludeAllModels.cpp` you have to add `#include <BSMPT/models/YourModel.h>` to the include list. Also, to be able to call your model, you have to extend the `FChoose` function. For this you add a new case to the switch statement, which reads

        case ModelIDs::YourModel: return std::make_unique<Class_YourModel>(); break;

### Generate the C++ code for a model
We provide currently three methods to generate the tensors and calculate the counter terms for a new model.

1. At `tools/ModelGeneration/sympy` we provide a setup using only `python3` with `sympy` (at least version 1.10!, if your packet manager only has an older installed, e.g. ubuntu 20.04 only has v1.6, then you have to install v1.10 or up with pip). Here we provide two examples, `SM.py` and `G2HDM.py` (generic 2HDM) which both use the `ModelGenerator.py` module to calculate the tensors and CT. You can get the CT using `python3 SM.py --show ct` and the tensors by calling `python3 SM.py --show tensors`. If your counterterms don't have a unique solution, then the solution space will be shown to you and you have to add additional equations until you have a unique solution (e.g. in the G2HDM example). To show the simplified tree-level and counterterm potentials, you can use `python3 SM.py --show treeSimpl` and `python3 SM.py --show CTSimpl`.
2. At `tools/ModelGeneration/Mathematica` we provide a Mathematica framework to implement your model.
3. At `tools/ModelGeneration/Maple` we provide the maple Worksheet `CreateModel.mw` which you can use to implement your model and get the tensors.


You can use the Test executable to detect possible errors in your implementation. If the Test executable does not show you an error, but something is still wrong, contact us at bsmpt@lists.kit.edu

Also contact us if you have a custom model for BSMPT v1.x and you have trouble converting it to the new notation.

## Executables
BSMPT provides multiple executables. Here we give a quick overview of them. 
Every executable can be called with the `--help` option to see how it can be run and to get an overview of all its required and optional arguments.
Also, consult the manuals ([BSMPTv1](https://arxiv.org/abs/1803.02846), [BSMPTv2](https://arxiv.org/abs/2007.01725), [BSMPTv3](https://arxiv.org/abs/24XX.XXXXX])) for more details on the executables and their input parameters.

Note, that every executable has the option to set the `--json=/path/to/your/file.json` which contains a json string with the parameters you can set through the CLI. This can be useful if you want to store the parameters you used for a given call. Please beware that all paths in the json file are considered relative to the current working directory and not to the location of the json file. Examples can be found in `example/JSON`. If you want to be sure to have the correct output file we recommend using absolute paths.


In BSMPTv3, the following four executables are added:


### MinimaTracer
MinimaTracer tracks temperature-dependent local minima in a user-defined temperature interval.

### CalcTemps
CalcTemps identifies regions of coexisting minima, calculates the bounce solutions and characteristic temperature scales (critical, nucleation, percolation and completion temperature) of first-order phase transitions. Based on that, we report a transition history for the point.

### CalcGW
CalcGW expands CalcTemps by the additional calculation of the gravitational waves spectra sourced by first-order phase transitions.

### PotPlotter
PotPlotter calculates user-defined data grids that can be used for the visualization of multi-dimensional potential contours.


The following executables were released with BSMPTv1 and BSMPTv2:


### Test
Test checks the model implementation for a provided parameter point. Some of the performed tests are e.g.: matching fermion masses and tree-level electroweak minimum with SM, tadpole relations, matching scalar masses between tree-level and NLO and symmetries of the coupling-tensors. The number of passed/failed tests is reported.

### BSMPT
BSMPT calculates the strength of a single-step electroweak phase transition (EWPT), defined as the ratio of the vacuum expectation value (VEV) at the critical temperature \f$v_c\f$ over the critical temperature \f$T_c\f$, based on finding a discontinuity in the electroweak VEV of the temperature-dependent global minimum.

### CalcCT
CalcCT calculates the (finite) counterterms for the 'on-shell' renormalization scheme.

### NLOVEV
NLOVEV calculates the zero-temperature VEV at one-loop order. This can be used to investigate the vacuum stability of the model.

### VEVEVO
VEVEVO calculates the evolution of the global minimum of a given point in a user-specified temperature range.

### TripleHiggsCouplingNLO
TripleHiggsCouplingNLO calculates the trilinear Higgs self-couplings at NLO at zero temperature.

### CalculateEWBG
CalculateEWBG calculates the difference between baryons and anti-baryons
normalized to the photon density generated through the EWPT.
Please beware that this is only tested for the C2HDM so far and the general
implementation is future work.

### PlotEWBG_vw
PlotEWBG_vw varies the wall velocity of a given parameter point and
calculates the baryon asymmetry for each velocity.

### PlotEWBG_nL
PlotEWBG_nL calculates the left-handed fermion density in front of the
wall as a function of the distance to the bubble wall.



[DoxygenLink]: https://phbasler.github.io/BSMPT/documentation
