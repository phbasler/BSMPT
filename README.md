Program: BSMPT version 1.1.2

Released by: Philipp Basler and Margarete Mühlleitner

Manual: version 1.0 

BSMPT - Beyond the Standard Model Phase Transitions:
The C++ program package BSMPT calculates for models with extended
Higgs sectors the loop-corrected effective potential at finite temperature
including the daisy resummation of the bosonic masses. The program
computes the vacuum expectation value (VEV) \f$ v \f$ of the potential
as a function of the temperature, and in particular the critical VEV
\f$v_c\f$ at the temperature \f$T_c\f$ where the phase transition takes
place. In addition, the loop-corrected trilinear Higgs self-couplings are
provided. We apply an 'on-shell' renormalization scheme in the sense
that the loop-corrected masses and mixing angles are required to be
equal to their tree-level input values. This allows for efficient
scans in the parameter space of the models. The models implemented so far
are the CP-conserving and the CP-violating 2-Higgs-Doublet Models (2HDM) and the
Next-to-Minimal 2HDM (N2HDM). The program structure is such that the
user can easily implement further models.

The program package can be downloaded at:
https://github.com/phbasler/BSMPT

A short documentation is given at: https://phbasler.github.io/BSMPT/

Sample input and output files are provided in the directory 'example'.

Modifications and corrected bugs are reported in the file 'Changelog.md'.



For additional information, comments, complaints or suggestions please e-mail
to:  bsmpt@lists.kit.edu

---
   

##Installation:


1) The program requires the GSL library which the code assumes to be installed in PATH. 

2) To compile the program type `mkdir build && cd build` and there `cmake ..` where this part can be modified with

    CC=CCompiler
    
    CXX=C++Compiler
    
    EIGEN3_ROOT=/path/to/eigen3  
        This is only necessary if eigen3 is not installed through your packet mananger. If it is not go to http://eigen.tuxfamily.org and download and unzip it to /path/to/eigen3.
    
    CMAES_ROOT=/path/to/libcmaes 
        If you do not give this option the default path will be the 'tools' folder in your BSMPT directory. For downloading libcmaes you need either wget or curl installed on your system.
    
    
For example a complete call would be `CC=gcc-7 CXX=g++-7 EIGEN3_ROOT=/path/to/eigen3 CMAES_ROOT=/path/to/libcmaes cmake ..` . After this execute `make` to compile the executables.

 
---
  
##How to add a new model (for further details, also see the manual):

1) Go to the file IncludeAllModels.h and rename the variable
   'C_ModelTemplate' to the variable name with which the new model shall
   be selected by the program.

2) In the file `src/models/CMakeLists.txt` rename `ClassTemplate.cpp` with the name of your new model file.

3) Go to IncludeAllModels.cpp and add the header for your model.  Also add

``` c++    
	  else if(choice == C_ModelTemplate)
     {
       return std::unique_ptr<Class_Potential_Origin> { new Class_Template };
     }
```


   to the function Fchoose. 

4) Adjust the functions in ClassTemplate.cpp as needed for the new model.

5) Have fun!


