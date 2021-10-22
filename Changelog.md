<!--
SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller

SPDX-License-Identifier: GPL-3.0-or-later
-->

# Changelog

# 2021/10/20: Release of v2.3.2
- Fixed wrong linker part if BSMPT::Minimizer was used in an external cmake package 
- Included tests to check if the Tensors are symmetric
- Use symmetric tensors to update performance 
- Bumped catch2 to 2.3.17

# 2021/09/23: Release of v2.3.1
- Fixed a bug in the CxSM which would cause Re(a1)=Im(a1)=0 if read a point from a file
- Added CxSM example point
- Added CxSM unit tests
- Fixed a race condition in the thermal interpolations causing segfaults on mac
- Bumped libcmaes to 0.10

# 2021/08/30: Release of v2.3
- Changed cout output to a Logger which provides different levels. How to set them is shown under the --help option.
- CP in the Dark is now added as a model. The implementation follows the conventions of [1807.10322], for further details on the chosen renormalization scheme we refer to our respective publication.

# 2021/05/06: Release of v2.2 
- Unit tests are provided for the models 'ch2dm', 'n2hdm' and 'r2hdm'
- Fixed the renormalisation scheme for the C2HDM 
- Added the option to turn off the multithreading in the minimization
- Removed the following binaries: 'RenormScale' and 'EWBGRenormScale'

# 2020/01/03: Release of v2.1
- Changed to multithreading in the minimizer
- Included a fix for deprecating cubic splines starting in boost 1.72
- Fixed an out-of-bounds error which caused a crash on mac
- Provided some first unit tests with the ctest environment
- Updated the cmaes inclusion to their new cmake setup
- Updated the NLopt inclusion to use their delivered FindCmake instead of our own


# 2020/07/06: Release of v2.0

## New Physics

With this version, we provide the calculation for the electroweak baryogenesis in the C2HDM with multiple, different approaches.
For this several smaller libraries have been implemented:
### Baryo
This library contains the numerical evaluation of the different methods to calculate the electroweak baryogenesis. The different approaches are described in the manual.
### WallThickness
The thickness of the bubble wall is calculated as described in the manual.
### ThermalFunctions
The thermal integrals J_+ and J_- have been moved to their own library. The old implementations are still in the main code and can be accessed for legacy code if wished. 
### Kfactors
During the calculation of the transport equations, the K-functions are needed. This library provides a numerical integration for those but also the bicubic spline interpolation which is used in the numerical evaluation of the transport equations. The data points for the interpolation are given in 'include/BSMPT/Kfactors/Kfactors_grid/Kfunctions_grid.h' which has a size of 148mb.
### Minimizer
NLopt is now a possible option for the minimization. You can set the default settings in Minimizer.h 
## Models
You CAN NOT call the models anymore by their number. You have to call them through the string set in IncludeAllModels.h
	
## BSMPT as a package
It is now possible to include BSMPT as a library into your program through cmake. In your cmake file you can use find_package(BSMPT) after compiling BSMPT.

## Changes in Test
* Extends the Test binary to check if the minimum conditions are fulfilled
* If a simplified tree-level or counterterm potential is given it is compared at random points with the potential calculated through the tensor structures to check for possible errors in the simplified potential.

## Changes in the Installation Routine
Due to changes in the cmake interfaces, Eigen is now included through -DEigen3_DIR. The instructions in the Readme are updated.


## Differences in how to include a new model
Due to the restructuring of the Code, the source file for your model has to be put in include/BSMPT/models/YourModel.h and the source file in src/models/YourModel.cpp. In src/models you have to put ${header_path}/YourModel.h in the set(header list and YourModel.cpp in the set(src list. 
	
    
# 2019/12/11: Release of v1.1.2
Many thanks for Andrew Fowlie, Csaba Balazs, Peter Athron and Yang Zhang for pointing out a small bug in the calculation of the bosonic thermal corrections. The finite shift delta_- introduced in Eq. (2.43) of the manual was not carried over in the case m^2 < 0, introducing a non-continuity in the function. This is now fixed.

#2019/10/08: Release of v1.1.1
v1.1 did not calculate the tree-level minimum for comparison with the input but the NLO minimum, this is now fixed.
The code can now work with data samples where the first column is an index column without a label, assuming the columns are tab separated.
Cleaned up the interface of the minimisation functions, this changed nothing at the results.
Introduced the initModel function which will handle ReadAndSet from the line, setting the parameters, calculating the Counterterm parameters and setting them as well as preparing everything necessary in the background. This was cleanup of the interface and did not change the numerics. 

# 2019/10/07: Release of v1.1
It is now possible to call the models with the names 'ch2dm', 'n2hdm', 'r2hdm' instead of the numbers. Additionally added the Test binary which provides a possibility to check the implemented models.

# 2019/09/27
Fixed a bug introduced with the changes of 2019/08/15

# 2019/08/05
Fixed a small bug which occurred if the Gauge fields were already given in a diagonal basis

# 2018/11/06
Updated the manual with the cmake installation and fixed some typos

# 2018/09/25
Fixed a small bug which set all singlet VEVs to 0 for T > T_C. This did not effect BSMPT but only VEVEVO for plotting for T > T_C in the N2HDM.
# 2018/07/11
Fixed a bug in the calculation of the triple Higgs couplings. Furthermore, the terminal output was changed into a more readable format. Thanks to Jonas Wittbrodt the installation is now done through a cmake file where the procedure is described in the README.
# 2018/07/02
The notation of omega_c and T_c was swapped in the N2HDM model file and in the example/N2HDM_Input.dat_BSMPT output file. This is fixed now.
# 2018/06/25
Fixed a mistake in eq (2.30) in the manual and the corresponding formula in the code. Thanks to Peter Athron for noticing this. Additionally, the numbers of minimizations were increased if only the GSL minimizer is used. Before it searched for 20 local minima, now it searches for 50. 
# 2018/04/23
Fixed a wrong sign in the Yukawa couplings to the CP-odd Higgs fields
# 2018/03/30
Corrected a small bug in the minimizer which caused a segfault if vevsolTmp > 0.5 and ModifiedVEVVectorDim has dimension 0
# 2018/03/07: v1.0
Release 
