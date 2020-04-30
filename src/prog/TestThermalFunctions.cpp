/*
 * TestThermalFunctions.cpp
 *
 *  Created on: 02.12.2019
 *      Author: phil
 */

#include <BSMPT/models/IncludeAllModels.h>
#include <iostream>
#include <BSMPT/ThermalFunctions/ThermalFunctions.h>
#include <BSMPT/ThermalFunctions/NegativeBosonSpline.h>
#include <bits/exception.h>                              // for exception
#include <math.h>                                        // for abs, M_PI
#include <stdlib.h>                                      // for size_t, EXIT...
#include <iomanip>                                       // for operator<<
#include <memory>                                        // for unique_ptr
#include <BSMPT/models/ClassPotentialOrigin.h>           // for Class_Potent...
#include <BSMPT/utility.h>

using namespace std;
using namespace BSMPT;
using namespace BSMPT::ThermalFunctions;

double JB(double msquared,double T){
    std::shared_ptr<Class_Potential_Origin> modelPointer = ModelID::FChoose(ModelID::ModelIDs::C2HDM);
	double x = msquared/(T*T);
	double xprev = 0, xnext = 0, fprev = 0, fnext = 0, PotVal = 0;
	if(x<0){
		if(-x >= C_NegLine-2) {
			xprev = NegLinearInt[C_NegLine-2][0];
			xnext = NegLinearInt[C_NegLine-1][0];
			fprev = NegLinearInt[C_NegLine-2][1];
			fnext = NegLinearInt[C_NegLine-1][1];
		}
		else{
			size_t pos = (size_t (-x));
			xprev = NegLinearInt[pos][0];
			fprev = NegLinearInt[pos][1];
			xnext = NegLinearInt[pos+1][0];
			fnext = NegLinearInt[pos+1][1];
		}

		PotVal = (fnext-fprev)/(xnext-xprev) * (x-xprev) + fprev;
		PotVal+= C_BosonShift;
		PotVal *= std::pow(T,4)/(2*M_PI*M_PI);
	}
	else{
		double xbp = C_BosonTheta;
		if ( x < xbp) {
			PotVal += modelPointer->Vsb(msquared, T, 3) + C_BosonShift*std::pow(T,4)/(2*M_PI*M_PI);
		} else if (x >= xbp) {
			PotVal += modelPointer->Vl(msquared, T, 3);
		}
	}

	return PotVal;
}


int main() try{

//	std::ofstream outfile("TestJB.tsv");
//	double Temp = 1e2;
//	outfile << "x\tJB" << std::endl;
//	double stepsize = 0.1;
//	for(int i=-100;i<=100;i++){
//		double x = i*stepsize;
//		double res = JB(x*std::pow(Temp,2),Temp);
//		outfile << x << sep << res << std::endl;
//	}
//	outfile.close();

    auto modelPointer = ModelID::FChoose(ModelID::ModelIDs::C2HDM);


	double StepSize =2;
	std::cout << "Boson check : " << std::endl;
	for(int i = 0;i<=5;i++)
	{
		double x = i*StepSize;
		double Integral = 1.0/(2*std::pow(M_PI,2)) * JbosonNumericalIntegration(x,0);
		double cb = 0.5;
		double OldVersion = modelPointer->boson_legacy(x,1,cb) - modelPointer->CWTerm(std::abs(x),cb);
		double NewVersion = modelPointer->boson(x,1,cb) - modelPointer->CWTerm(std::abs(x),cb);
		double DeltaIntOld = Integral-OldVersion, DeltaIntNew = Integral-NewVersion, DeltaOldNew = OldVersion-NewVersion;
		std::cout << std::scientific << std::setprecision(10) << "Jb(" << x << ") = "
				<< Integral
                << sep << " Old = " << OldVersion
                << sep << " New = " << NewVersion
                << sep << " Delta(Int-Old) = " << DeltaIntOld
                << sep << " Delta(Int-New) = " << DeltaIntNew
                << sep << " Delta(Old-New) = " << DeltaOldNew
				<< std::endl;
	}
	std::cout << "Fermion check : " << std::endl;
	for(int i = 0;i<=5;i++)
	{
		double x = i*StepSize;
		double Integral = 1.0/(2*std::pow(M_PI,2)) * JfermionNumericalIntegration(x,0);
        double OldVersion = modelPointer->fermion_legacy(x,1) - modelPointer->CWTerm(std::abs(x),C_CWcbFermion);
		double NewVersion = modelPointer->fermion(x,1) - modelPointer->CWTerm(std::abs(x),C_CWcbFermion);
		double DeltaIntOld = Integral-OldVersion, DeltaIntNew = Integral-NewVersion, DeltaOldNew = OldVersion-NewVersion;
		std::cout << std::scientific << std::setprecision(10) << "Jf(" << x << ") = "
				<< Integral
                << sep << " Old = " << OldVersion
                << sep << " New = " << NewVersion
                << sep << " Delta(Int-Old) = " << DeltaIntOld
                << sep << " Delta(Int-New) = " << DeltaIntNew
                << sep << " Delta(Old-New) = " << DeltaOldNew
				<< std::endl;
	}

	return EXIT_SUCCESS;
}

catch(exception& e){
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
}
