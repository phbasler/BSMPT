2019/10/08: Release of v1.1.1:
v1.1 did not calculate the tree-level minimum for comparison with the input but the NLO minimum, this is now fixed.
The code can now work with data samples where the first column is an index column without a label, assuming the columns are tab separated.
Cleaned up the interface of the minimisation functions, this changed nothing at the results.
Introduced the initModel function which will handle ReadAndSet from the line, setting the parameters, calculating the Counterterm parameters and setting them as well as preparing everything necessary in the background. This was cleanup of the interface and did not change the numerics. 
2019/10/07: Release of v1.1. It is now possible to call the models with the names 'ch2dm', 'n2hdm', 'r2hdm' instead of the numbers. Additionally added the Test binary which provides a possibility to check the implemented models.
2019/09/27: Fixed a bug introduced with the changes of 2019/08/15
2019/08/05: Fixed a small bug which occurred if the Gauge fields were already given in a diagonal basis
2018/11/06: Updated the manual with the cmake installation and fixed some typos
2018/09/25: Fixed a small bug which set all singlet VEVs to 0 for T > T_C. This did not effect BSMPT but only VEVEVO for plotting for T > T_C in the N2HDM.
2018/07/11: Fixed a bug in the calculation of the triple Higgs couplings. Furthermore, the terminal output was changed into a more readable format. Thanks to Jonas Wittbrodt the installation is now done through a cmake file where the procedure is described in the README.
2018/07/02: The notation of omega_c and T_c was swapped in the N2HDM model file and in the example/N2HDM_Input.dat_BSMPT output file. This is fixed now.
2018/06/25 : Fixed a mistake in eq (2.30) in the manual and the corresponding formula in the code. Thanks to Peter Athron for noticing this. Additionally, the numbers of minimizations were increased if only the GSL minimizer is used. Before it searched for 20 local minima, now it searches for 50. 
2018/04/23: Fixed a wrong sign in the Yukawa couplings to the CP-odd Higgs fields
2018/03/30: Corrected a small bug in the minimizer which caused a segfault if vevsolTmp > 0.5 and ModifiedVEVVectorDim has dimension 0
2018/03/07: v1.0: Release 
