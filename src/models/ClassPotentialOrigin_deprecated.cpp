/*
 * ClassPotentialOrigin.cpp
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




#include <iomanip>
#include <gsl/gsl_sf_gamma.h>
#include <random>

#include "Eigen/Dense"

#include <gsl/gsl_sf_zeta.h>

#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/ThermalFunctions/ThermalFunctions.h>
#include <BSMPT/ThermalFunctions/NegativeBosonSpline.h>
#include <BSMPT/utility.h>
using namespace Eigen;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"


namespace BSMPT{


double Class_Potential_Origin::Vl(double MassSquared, double Temp, int n,int diff) const
{

    double Mass = std::sqrt((MassSquared));
    double PotVal;
    PotVal = 0;
    if(diff == 0)
      {
        for(int l=0; l<=n; l++)
            {
                PotVal += 1/(std::pow(2,l)*gsl_sf_fact(l)) * gsl_sf_gamma(2.5+l)/gsl_sf_gamma(2.5-l) * std::pow(Temp/Mass,l);
            }
            PotVal *= -std::exp(-Mass/Temp) * std::pow(Mass*Temp/(2*M_PI),1.5)*Temp;
      }
    else if(diff >0){
        for(int l=0;l<=n;l++)
          {
        PotVal += 1/(std::pow(2,l)*gsl_sf_fact(l)) * gsl_sf_gamma(2.5+l)/gsl_sf_gamma(2.5-l)*std::pow(Temp/Mass,l)*(2*Temp*l+2*Mass-3*Temp);
          }
        PotVal *= 1.0/(16*M_PI) * (Temp/Mass) * std::exp(-Mass/Temp) * std::sqrt(2*Mass*Temp/M_PI);
    }
    else if(diff == -1){
        double Sum = 0, SumDerivative = 0;
        for(int l=0; l<=n; l++)
        {
            Sum += 1/(std::pow(2,l)*gsl_sf_fact(l)) * gsl_sf_gamma(2.5+l)/gsl_sf_gamma(2.5-l) * std::pow(Temp/Mass,l);
            SumDerivative += 1/(std::pow(2,l)*gsl_sf_fact(l)) * gsl_sf_gamma(2.5+l)/gsl_sf_gamma(2.5-l) * std::pow(Temp/Mass,l) * l/Temp;
        }
        PotVal = std::pow(Temp,2) *Sum*std::sqrt(MassSquared)/Temp+5*Sum*std::pow(Temp,2)/2-2*SumDerivative*MassSquared;
        PotVal *= -std::sqrt(2)*std::pow(MassSquared,2)*std::exp(-std::sqrt(MassSquared)/Temp);
        PotVal *= 1.0/(4*std::pow(M_PI,1.5)*std::sqrt(Temp) * std::pow(MassSquared,2.5));

    }

    return PotVal;
}

double Class_Potential_Origin::Vsf(double MassSquared, double Temp, int n, int diff) const
{
    double PotVal=0;
//	std::cout << scale << std::endl;
    double logM = FCW(MassSquared) + std::log(scale*scale);
//	std::cout << logM << std::endl;
    double cf = 1.5+2*std::log(4*M_PI)-2*ThermalFunctions::C_euler_gamma-2*std::log(4);
    if(diff == 0)
      {
        PotVal = -7*M_PI*M_PI/(720.0)*std::pow(Temp,4);
        PotVal += MassSquared*Temp*Temp/(48.0);
        PotVal += std::pow(MassSquared,2)/(64*M_PI*M_PI)*(logM-std::log(Temp*Temp)- cf);
        for(int l=2; l<=n; l++)
        {
            PotVal -= MassSquared*Temp*Temp/2.0 * std::pow(-MassSquared/(4*M_PI*M_PI*Temp*Temp),l)*gsl_sf_doublefact(2*l-3) * gsl_sf_zeta(2*l-1) *1/(gsl_sf_doublefact(2*l)) * 1.0/(l+1) * ( std::pow(2,2*l-1) - 1 );
        }
      }
    else if(diff > 0){
        PotVal = Temp*Temp/(48.0);
        PotVal += MassSquared*(logM-std::log(Temp*Temp) - cf)/(32.0*M_PI*M_PI);
        PotVal += MassSquared/(64.0*M_PI*M_PI);
        for(int l=2;l<=n;l++)
          {
            double Kl = gsl_sf_doublefact(2*l-3) * gsl_sf_zeta(2*l-1) *1/(gsl_sf_doublefact(2*l)) * 1.0/(l+1) * ( std::pow(2,2*l-1) - 1 );
            PotVal += -0.5*Temp*Temp*std::pow(-MassSquared/(4.0*M_PI*M_PI*Temp*Temp),l)*Kl*(l+1);
          }
    }
    else if(diff == -1 ){
        double Sum = 0, SumDerivative = 0;
        for(int l=2;l<=n;l++){
            double Kl = gsl_sf_doublefact(2*l-3) * gsl_sf_zeta(2*l-1) *1/(gsl_sf_doublefact(2*l)) * 1.0/(l+1) * ( std::pow(2,2*l-1) - 1 );
            Sum += std::pow(-MassSquared/std::pow(Temp,2) *1.0/(4*std::pow(M_PI,2)),l) *Kl;
            SumDerivative+= -2/Temp * std::pow(-MassSquared/(4*std::pow(M_PI*Temp,2)),l )*Kl;
        }
        PotVal = -7*std::pow(Temp,3)*std::pow(M_PI,2)/180.0;
        PotVal += Temp*MassSquared/24.0;
        PotVal += -Temp*MassSquared *Sum;
        PotVal += -std::pow(MassSquared,2)/(32*Temp*std::pow(M_PI,2));
        PotVal += std::pow(MassSquared,2) *SumDerivative /Temp;

    }

    return PotVal;
}

double Class_Potential_Origin::Vsb(double MassSquaredIn, double Temp, int n, int diff) const
{
    double MassSquared = (MassSquaredIn);
    double PotVal=0;
    double LogTerm = FCW(MassSquared) + std::log(scale*scale);
    double cb = 1.5+2*std::log(4*M_PI)-2*ThermalFunctions::C_euler_gamma;
    if(diff == 0)
      {
        PotVal = -M_PI*M_PI/90.0*std::pow(Temp,4);
        PotVal += MassSquared*Temp*Temp/(24.0);
        PotVal -= std::pow(MassSquared,1.5)*Temp/(12.0*M_PI);
        PotVal -= std::pow(MassSquared,2.0)/(64*M_PI*M_PI)*(LogTerm - std::log(Temp*Temp)-cb);
        for(int l=2; l<=n; l++)
        {
            PotVal += MassSquared*Temp*Temp*0.5*std::pow(-MassSquared/(4*M_PI*M_PI*Temp*Temp),l)*gsl_sf_doublefact(2*l-3) * gsl_sf_zeta(2*l-1) *1/(gsl_sf_doublefact(2*l)) * 1.0/(l+1);
        }
      }
    else if(diff > 0){
        PotVal = Temp*Temp/(24.0);
        PotVal += -1.0/8.0*std::sqrt(MassSquared)*Temp/M_PI;
        PotVal += - MassSquared*(LogTerm - std::log(Temp*Temp)-cb)/(32.0*M_PI*M_PI);
        PotVal += -MassSquared/(64.0*M_PI*M_PI);
        for(int l=2;l<=n;l++)
          {
        double Kl= gsl_sf_doublefact(2*l-3) * gsl_sf_zeta(2*l-1) *1/(gsl_sf_doublefact(2*l)) * 1.0/(l+1);
        PotVal += 0.5*Temp*Temp*std::pow(-MassSquared/(4*M_PI*M_PI*Temp*Temp),l)*Kl*(l+1);

          }
    }
    else if(diff == -1){
        double Sum = 0, SumDerivative = 0;
        for(int l=2;l<=n;l++){
            double Factorial = gsl_sf_doublefact(2*l-3) * gsl_sf_zeta(2*l-1) *1/(gsl_sf_doublefact(2*l)) * 1.0/(l+1);
            Sum += std::pow(-MassSquared/(4*std::pow(M_PI*Temp,2)),l)*Factorial;
            SumDerivative += -2/Temp * std::pow(-MassSquared/(4*std::pow(M_PI*Temp,2)),l)*Factorial;
        }
        PotVal = -2*std::pow(Temp,3)*std::pow(M_PI,2)/45.0;
        PotVal +=  Temp *MassSquared *Sum;
        PotVal += MassSquared*Temp/12.0;
        if(MassSquared > 0) PotVal += -MassSquared *std::sqrt(MassSquared);
        if(MassSquared > 0) PotVal += -std::pow(Temp,3)*std::log(MassSquared/std::pow(Temp,2))/(16*M_PI*M_PI);
        PotVal +=  std::pow(Temp,3)*cb/(16*M_PI*M_PI);
        if(MassSquared > 0) PotVal += std::pow(MassSquared,2)/(12*M_PI*std::sqrt(MassSquared));
        PotVal += std::pow(MassSquared,2)*SumDerivative/Temp;
        PotVal += std::pow(Temp,3)/(32*M_PI*M_PI);
    }


    return PotVal;
}

double Class_Potential_Origin::boson_legacy(double MassSquared, double Temp, double cb, int diff) const
{
    double resPotVal=0;
    long double PotVal = 0;
    long double PotCW = CWTerm(std::abs(MassSquared),cb, diff);
    if(Temp == 0)
    {

        return (double) PotCW;
    }
    if(diff==0)
      {
        if(MassSquared < 0 )
        {


            double xprev,fprev,xnext,fnext;
            double x = MassSquared/(Temp*Temp);
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
            PotVal += ThermalFunctions::C_BosonShift;
            PotVal *= std::pow(Temp,4)/(2*M_PI*M_PI);

        }
        else{
            double xbp = ThermalFunctions::C_BosonTheta;
            if (MassSquared / (Temp * Temp) < xbp) {
                PotVal += Vsb(MassSquared, Temp, 3) + ThermalFunctions::C_BosonShift*std::pow(Temp,4)/(2*M_PI*M_PI);
            } else if (MassSquared / (Temp * Temp) >= xbp) {
                PotVal += Vl(MassSquared, Temp, 3);
            }
        }
    }
    else if(diff >0){
        if(std::abs(MassSquared/std::pow(Temp,2))<= C_threshold) return Temp*Temp/(24.0);
        else if(MassSquared >= 0)
        {
            double xbp = ThermalFunctions::C_BosonTheta;
            if (MassSquared / (Temp * Temp) < xbp) {
                PotVal += Vsb(MassSquared, Temp, 3,diff);
            }
            else if (MassSquared / (Temp * Temp) >= xbp) {
                PotVal += Vl(MassSquared, Temp, 3,diff);
            }
        }
        else{
            double x = -MassSquared/std::pow(Temp,2);
            bool NegSpline = true;
            if(NegSpline)
                {
                    PotVal =(std::log(2)/4.0 - 0.691643) * std::sqrt(std::abs(x));
                    PotVal *= Temp*Temp/(2*M_PI*M_PI);
                }
            else{
            PotVal = 0;
            std::cerr << "This is not implemented yet! You called the derivative of the bosonic"
                << " size_tegral for negative m^2" << std::endl;
            }
        }
    }
    resPotVal = (double) (PotVal+PotCW);
    return resPotVal;

}

double Class_Potential_Origin::fermion_legacy(double MassSquared, double Temp, int diff) const
{
    long double PotVal = 0;
    long double PotCW = CWTerm(MassSquared,C_CWcbFermion,diff);
    double resPotVal= (double) PotCW;
    if(Temp == 0) return resPotVal;
    double x = MassSquared/std::pow(Temp,2);
    if(diff == 0)
    {
        if(x < ThermalFunctions::C_FermionTheta)
        {
            PotVal = -(Vsf(MassSquared,Temp,4)+ThermalFunctions::C_FermionShift*std::pow(Temp,4)/(2*M_PI*M_PI));
        }
        else if( x >= ThermalFunctions::C_FermionTheta)
        {
            PotVal = -(Vl(MassSquared,Temp,3));
        }
    }
    else{
        if(x < ThermalFunctions::C_FermionTheta)
        {
            PotVal = -(Vsf(MassSquared,Temp,4,diff));
            if(diff == -1){
                PotVal += -ThermalFunctions::C_FermionShift * 4*std::pow(Temp,3)/(2*M_PI*M_PI);
            }
        }
        else if( x >= ThermalFunctions::C_FermionTheta)
        {
            PotVal = -(Vl(MassSquared,Temp,3,diff));
        }
    }


    resPotVal = (double) (PotVal+PotCW);
    return resPotVal;
}



}

#pragma GCC diagnostic pop
