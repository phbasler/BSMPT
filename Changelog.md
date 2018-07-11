2018/07/11 : Fixed a bug in the calculation of the triple Higgs couplings. Furthermore the terminal output was changed into a more readable format. Thanks to Jonas Wittbrodt the installation is now done through a cmake file where the procedure is described in the README.
2018/07/02 : The notation of omega_c and T_c was swapped in the N2HDM model file and in the example/N2HDM_Input.dat_BSMPT output file. This is fixed now.
2018/06/25 : Fixed a mistake in eq (2.30) in the manual and the corresponding formula in the code. Thanks to Peter Athron for noticing this. Additionally the numbers of minimizations was increased if only the GSL minimizer is used. Before it searched for 20 local minima, now it searches for 50. 
2018/04/23 : Fixed a wrong sign in the yukawa couplings to the CP-odd higgs fields
2018/03/30 : Corrected a small bug in the minimizer which caused a seg fault if vevsolTmp > 0.5 and ModifiedVEVVectorDim has dimension 0
2018/03/07 : v1.0 : Relase 
