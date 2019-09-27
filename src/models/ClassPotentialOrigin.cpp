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

#include "ClassPotentialOrigin.h"
#include "IncludeAllModels.h"
using namespace Eigen;


Class_Potential_Origin::Class_Potential_Origin ()
{
  // TODO Auto-generated constructor stub

}

Class_Potential_Origin::~Class_Potential_Origin ()
{
  // TODO Auto-generated destructor stub
}

/**
   * This will call set_gen(par), SetCurvatureArrays, set_CT_Pot_Par(parCT), CalculateDebye() as well as
   * CalculateDebyeGauge()
   */
void Class_Potential_Origin::set_All(const std::vector<double>& par,const std::vector<double>& parCT)
{


  set_gen(par);
  if(!SetCurvatureDone) SetCurvatureArrays();
  set_CT_Pot_Par(parCT);
  CalculateDebye();
  CalculateDebyeGauge();

}

/**
 * Sets a tensor needed to calculate the contribution of the counterterm potential to the triple Higgs couplings.
 */
void Class_Potential_Origin::Prepare_Triple()
{
	for(int a=0;a<NHiggs;a++)
	{
		for(int b=0;b<NHiggs;b++)
		{
			for(int i=0;i<NHiggs;i++)
			{
				LambdaHiggs_3_CT[a][b][i] = Curvature_Higgs_CT_L3[a][b][i];
				for(int j=0;j<NHiggs;j++)
				{
					LambdaHiggs_3_CT[a][b][i] += Curvature_Higgs_CT_L4[a][b][i][j]*HiggsVev[j];
				}
			}
		}
	}

}

/**
 * Calculates the small m^2/T^2 approximation of the bosonic temperature dependent integral to the n-th order.
 */
double Class_Potential_Origin::Vsb(double MassSquaredIn, double Temp, int n, int diff)
{
	double MassSquared = (MassSquaredIn);
	double PotVal;
	double LogTerm = FCW(MassSquared) + std::log(scale*scale);
	double cb = 1.5+2*std::log(4*M_PI)-2*C_euler_gamma;
	double cf= cb - 2*std::log(4);
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
	else{
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


	return PotVal;
}

/**
 * Calculates Re(log(MassSquared)) and returns 0 if the argument is too small as this function is only called with an (m^2)^n
 * in front of it.
 */
double Class_Potential_Origin::FCW(double MassSquared)
{
	bool Reg = false;
	double res=0;
	double x;
	double Boarder=std::pow(10,-200);
	if(std::isnan(MassSquared)) x = Boarder;
	else if(std::abs(MassSquared) < Boarder) x = Boarder;
	else x = std::abs(MassSquared);
	res = std::log(x) - 2*std::log(scale);
	return res;

}


/**
 * Calculates the Coleman-Weinberg contribution of a particle with m^2 =  MassSquared and the constant scheme-dependent parameter cb
 * as well as the first derivative with respect to m^2.
 */
double Class_Potential_Origin::CWTerm(double MassSquared, double cb, int diff)
{

    if(std::abs(MassSquared) < C_threshold) return 0;
   double LogTerm=0, PotVal=0;
   LogTerm = FCW(MassSquared);
   if(diff == 0) PotVal = 1.0/(64*M_PI*M_PI) * MassSquared*MassSquared*(LogTerm-cb);
   else{
    PotVal = 1.0/(32*M_PI*M_PI)*MassSquared*(LogTerm-cb+0.5);

   }
   return PotVal;
}


/**
 * Calculates the large m^2/T^2 approximation to order n for the temperature-dependent integrals.
 */
double Class_Potential_Origin::Vl(double MassSquared, double Temp, int n,int diff)
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
	else if(diff != 0){
	    for(int l=0;l<=n;l++)
	      {
		PotVal += 1/(std::pow(2,l)*gsl_sf_fact(l)) * gsl_sf_gamma(2.5+l)/gsl_sf_gamma(2.5-l)*std::pow(Temp/Mass,l)*(2*Temp*l+2*Mass-3*Temp);
	      }
	    PotVal *= 1.0/(16*M_PI) * (Temp/Mass) * std::exp(-Mass/Temp) * std::sqrt(2*Mass*Temp/M_PI);
	}

	return PotVal;
}


/**
 * Calculates the small m^2/T^2 approximation of the fermionic temperature dependent integral to the n-th order.
 */
double Class_Potential_Origin::Vsf(double MassSquared, double Temp, int n, int diff)
{
	double PotVal=0;
//	std::cout << scale << std::endl;
	double logM = FCW(MassSquared) + std::log(scale*scale);
//	std::cout << logM << std::endl;
	double cf = 1.5+2*std::log(4*M_PI)-2*C_euler_gamma-2*std::log(4);
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
	else{
	    PotVal = Temp*Temp/(48.0);
	    PotVal += MassSquared*(logM-std::log(Temp*Temp) - cf)/(32.0*M_PI*M_PI);
	    PotVal += MassSquared/(64.0*M_PI*M_PI);
	    for(int l=2;l<=n;l++)
	      {
		double Kl = gsl_sf_doublefact(2*l-3) * gsl_sf_zeta(2*l-1) *1/(gsl_sf_doublefact(2*l)) * 1.0/(l+1) * ( std::pow(2,2*l-1) - 1 );
		PotVal += -0.5*Temp*Temp*std::pow(-MassSquared/(4.0*M_PI*M_PI*Temp*Temp),l)*Kl*(l+1);
	      }
	}

	return PotVal;
}


/**
 * @brief Calculation of the bosonic integral + Coleman-Weinberg potential __without__ taking d.o.f. into account
 *
 * @param MassSquared = m^2 of the particle
 * @param Temp = Temperature at which the Debye masses and integrals should be calculated
 * @param cb = Parameter of the renormalisation-Scheme in the Coleman-Weinberg potential
 * @param diff: 0 returns the value of the integral and diff != 0 the derivative w.r.t. m^2
 *
 * Not implemented yet:: Calculate the derivative for m^2 < 0
 *
 */
double Class_Potential_Origin::boson(double MassSquared, double Temp, double cb, int diff=0)
{
	bool NegativIntegral=true;
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
				int pos = (int (-x));
				xprev = NegLinearInt[pos][0];
				fprev = NegLinearInt[pos][1];
				xnext = NegLinearInt[pos+1][0];
				fnext = NegLinearInt[pos+1][1];
			}

			PotVal = (fnext-fprev)/(xnext-xprev) * (x-xprev) + fprev;
			PotVal *= std::pow(Temp,4)/(2*M_PI*M_PI);

        }
        else{
            double xbp = C_BosonTheta;
            if (MassSquared / (Temp * Temp) < xbp) {
                PotVal += Vsb(MassSquared, Temp, 3) + C_BosonShift*std::pow(Temp,4)/(2*M_PI*M_PI);
            } else if (MassSquared / (Temp * Temp) >= xbp) {
                PotVal += Vl(MassSquared, Temp, 3);
            }
        }
    }
    else{
        if(std::abs(MassSquared/std::pow(Temp,2))<= C_threshold) return Temp*Temp/(24.0);
        else if(MassSquared >= 0)
        {
            double xbp = C_BosonTheta;
	    if (MassSquared / (Temp * Temp) < xbp) {
		PotVal += Vsb(MassSquared, Temp, 3,1);
	    } else if (MassSquared / (Temp * Temp) >= xbp) {
		PotVal += Vl(MassSquared, Temp, 3,1);
	    }
        }
        else{
            double x = -MassSquared/std::pow(Temp,2);
            double a = -MassSquared;
            bool NegSpline = true;
            if(NegSpline)
            	{
            		PotVal =(std::log(2)/4.0 - 0.691643) * std::sqrt(std::abs(x));
            		PotVal *= Temp*Temp/(2*M_PI*M_PI);
            	}
            else{
        	PotVal = 0;
        	std::cerr << "This is not implemented yet! You called the derivative of the bosonic"
        	    << " Integral for negative m^2" << std::endl;
            }
        }
    }
    resPotVal = (double) (PotVal+PotCW);
    return resPotVal;

}

/**
 * @brief Calculation of the fermionic integral + Coleman-Weinberg potential __without__ taking d.o.f. into account
 *
 * @param MassSquared = m^2 of the particle
 * @param Temp = temperature at which the integrals should be calculated
 * @param diff :  0 = Value of the potential, i != 0 returns the derivative w.r.t. m^2
 */

double Class_Potential_Origin::fermion(double MassSquared, double Temp, int diff)
{
    long double PotVal = 0;
    long double PotCW = CWTerm(MassSquared,C_CWcbFermion,diff);
    double resPotVal= (double) PotCW;
    if(Temp == 0) return resPotVal;
    double x = MassSquared/std::pow(Temp,2);
    if(diff == 0)
    {
	if(x < C_FermionTheta)
	{
	    PotVal = -(Vsf(MassSquared,Temp,4)+C_FermionShift*std::pow(Temp,4)/(2*M_PI*M_PI));
	}
	else if( x >= C_FermionTheta)
	{
	    PotVal = -(Vl(MassSquared,Temp,3));
	}
    }
    else{
	if(x < C_FermionTheta)
	{
	    PotVal = -(Vsf(MassSquared,Temp,4,1));
	}
	else if( x >= C_FermionTheta)
	{
	    PotVal = -(Vl(MassSquared,Temp,3,1));
	}
    }


    resPotVal = (double) (PotVal+PotCW);
    return resPotVal;
}


/**
 * Calculates the first derivatives of the eigenvalues of a given matrix
 * @param M : the original matrix
 * @param MDiff : the element-wise first derivative of the matrix M with respect to the parameter you want to consider
 * @param res : writes the mass eigenvalues in the vector and then the derivatives in the same order
 */
void Class_Potential_Origin::FirstDerivativeOfEigenvalues(const Ref<MatrixXd> M, const Ref<MatrixXd> MDiff, std::vector<double> &res)
{
	const int nRows = M.rows();
	const int nCols = M.cols();

	const int EVThres = std::pow(10,-6);

	if(nCols != nRows) {
		std::cout << "ERROR ! M needs to be an quadratic Matrix for calculating the derivatives !\n";
		return ;
	}

	const int nSize = nRows;

	SelfAdjointEigenSolver<MatrixXd> es;
	es.compute(M);

	double Eigenvalues[nSize];
	double Derivatives[nSize];
	double AlreadyCalculated[nSize]; //Array to check which EVs already been calculated.
	for(int i=0;i<nSize;i++) {
        Eigenvalues[i] = es.eigenvalues()[i];
        if(std::abs(Eigenvalues[i]) < EVThres) Eigenvalues[i] = 0;
    }

	double Mapping[nSize][nSize];

	for(int i=0;i<nSize;i++)
	{
		AlreadyCalculated[i] = -1;
		for(int j=i;j<nSize;j++)
		{
			if(std::abs(Eigenvalues[i]-Eigenvalues[j]) > EVThres)
			{
				Mapping[i][j] = 0;
			}
			else{
				Mapping[i][j] = 1;
			}
		}
	}
	for(int i=1;i<nSize;i++)
	{
		for(int j=0;j<i;j++) Mapping[i][j] = Mapping[j][i];
	}

	for(int p=0;p<nSize;p++)
	{
		if(AlreadyCalculated[p] == -1)
		{
			int NumOfReps = 0;
			for(int i=p+1;i<nSize;i++)
			{
				NumOfReps += Mapping[p][i];
			}
			if(NumOfReps == 0)
			{
				VectorXd v(nSize);
				v = es.eigenvectors().col(p);
				Derivatives[p] = (v.transpose()*MDiff*v).value();
				AlreadyCalculated[p] = 1;
			}
			else{
				MatrixXd Phi(nSize,NumOfReps+1);
				int helpCol = 0;
				MatrixXd MXWork(NumOfReps+1,NumOfReps+1);

				for(int i=p;i<nSize;i++)
				{
					if(Mapping[p][i] == 1)
					{
						Phi.col(helpCol) = es.eigenvectors().col(i);
						helpCol++;
					}
				}
				MXWork = Phi.transpose()*MDiff*Phi;
				SelfAdjointEigenSolver<MatrixXd> esWork(MXWork,EigenvaluesOnly);
				helpCol = 0;
				for(int i=p;i<nSize;i++)
				{
					if(Mapping[p][i] == 1)
					{
						AlreadyCalculated[i] = 1;
						Derivatives[i] = esWork.eigenvalues()[helpCol];
						helpCol++;
					}
				}

			}
		}

	}
    for(int i=0;i<nSize;i++)
    {
        if(std::abs(Derivatives[i]) < EVThres) Derivatives[i] = 0;
    }

	for(int i=0;i<nSize;i++)
	  {
	    res.push_back(Eigenvalues[i]);
	  }
	  for(int i=0;i<nSize;i++) res.push_back(Derivatives[i]);





}


/**
 * Calculates the function f^2 needed for the 3rd derivatives of the Coleman Weinberg potential.
 */
double Class_Potential_Origin::fbaseTri(double MassSquaredA, double MassSquaredB, double MassSquaredC)
{
    double res = 0;
    double mas = MassSquaredA;
    double mbs = MassSquaredB;
    double mcs = MassSquaredC;
    double LogA = 0, LogB = 0, LogC = 0;
    double Thres = 1e-8;
    if(std::abs(mas) < Thres) mas = 0;
    if(std::abs(mbs) < Thres) mbs = 0;
    if(std::abs(mcs) < Thres) mcs = 0;
    if(std::abs(mas-mbs) < Thres) mas = mbs;
    if(std::abs(mas-mcs) < Thres) mas = mcs;
    if(std::abs(mbs-mcs) < Thres) mbs = mcs;

    if(mas != 0) LogA = std::log(mas) - 2*std::log(scale);
    if(mbs != 0) LogB = std::log(mbs) - 2*std::log(scale);
    if(mcs != 0) LogC = std::log(mcs) - 2*std::log(scale);

    int C = 1;
    if(mas == 0 and mbs == 0 and mcs == 0) res = 0;
    else if(mas != 0 and mbs == 0 and mcs == 0)
      {
	C=2;
	res = 1.0/mas*(std::log(mas) - 1);
	res = 1.0/mas*(LogA-1);
      }
    else if(mas == 0 and mbs != 0 and mcs == 0)
      {
	C=3;
	res = 1.0/mbs*(std::log(mbs)-1);
	res = (LogB-1)/mbs;
      }
    else if(mas == 0 and mbs == 0 and mcs != 0)
      {
	C=4;
	res = 1.0/mcs*(std::log(mcs)-1);
	res = (LogC-1)/mcs;
      }
    else if(mas == mbs and mas != 0 and mas != mcs and mcs != 0)
      {
	C=6;
	res = (mbs-mcs+mcs*std::log(mcs/mbs))/std::pow(mbs-mcs,2);
      }
    else if(mas == mcs and mas != 0 and mas != mbs and mbs != 0 )
      {
	C=7;
	res = (mbs*log(mbs/mcs) - mbs + mcs)/std::pow(mbs-mcs,2);
      }
    else if(mbs == mcs and mas != 0 and mbs != mas and mbs != 0 )
      {
	C=8;
	res = (mas*std::log(mas/mcs) - mas + mcs)/std::pow(mas-mcs,2);
      }
    else if(mas == mbs and mas == mcs and mas != 0)
      {
	C=9;
	res = 1.0/(2*mcs);
      }
    else if(mas == mbs and mas != mcs and mas != 0 and mcs == 0 )
      {
	C=10;
	res = 1.0/mbs;
      }
    else if(mas == mcs and mas != mbs and mas != 0 and mbs == 0)
      {
	C=11;
	res = 1.0/mas;
      }
    else if(mbs == mcs and mbs != 0 and mbs != mas and mas == 0 )
      {
	C=12;
	res = 1.0/mbs;
      }
    else
      {
    	C=5;
    	res = mas*LogA/((mas-mbs)*(mas-mcs)) + mbs*LogB/((mbs-mas)*(mbs-mcs));
    	res += mcs*LogC/((mcs-mas)*(mcs-mbs));
      }

    if(std::isnan(res) or std::isinf(res))   std::cout << mas << "\t" << mbs << "\t" << mcs << "\t" << res << "\t" << C << std::endl;



    return res;
}

/**
 * Calculates the function f^1 needed for the derivatives of the Coleman Weinberg potential.
 */
double Class_Potential_Origin::fbase(double MassSquaredA,double MassSquaredB)
{
    double res = 0;
    double LogA = 0;
    if(MassSquaredA == 0 and MassSquaredB == 0 ) return 1;
    double ZB = std::pow(10,-5);
    if(MassSquaredA != 0) LogA = std::log(MassSquaredA) - 2*std::log(scale);
    if(std::abs(MassSquaredA-MassSquaredB) > ZB)
    {
        double LogB = 0;
        if(MassSquaredB != 0)LogB = std::log(MassSquaredB) - 2*std::log(scale);
        if(MassSquaredA == 0) res = LogB;
        else if(MassSquaredB == 0 ) res = LogA;
        else res = (LogA*MassSquaredA - LogB*MassSquaredB)/(MassSquaredA-MassSquaredB);
    }
    else{
        res = 1 + LogA;
    }
    return res;
}

/**
 * This function calculates the second derivatives of all eigenvalues.
 * The matrix must not have a repeated eigenvalue for this!
 * @param M : the original matrix
 * @param MDiffX : the element-wise first derivative of the matrix M with respect to the first parameter you want to consider
 * @param MDiffY : the element-wise first derivative of the matrix M with respect to the second parameter you want to consider
 * @param MDiffXY : the element-wise second derivative of the matrix M with respect to both parameters you want to consider
 * @param res : writes the mass eigenvalues in the vector and then the derivatives in the same order
 */
void Class_Potential_Origin::SecondDerivativeOfEigenvaluesNonRepeated(const Eigen::Ref<Eigen::MatrixXd> M,const Eigen::Ref<Eigen::MatrixXd> MDiffX,const Eigen::Ref<Eigen::MatrixXd> MDiffY,const Eigen::Ref<Eigen::MatrixXd> MDiffXY, std::vector<double> &res)
{
  const int nRows = M.rows();
  const int nCols = M.cols();

  const int EVThres = std::pow(10,-6);

  if(nCols != nRows) {
	  std::cout << "ERROR ! M needs to be an quadratic Matrix for calculating the derivatives !\n";
	  return ;
  }

  const int nSize = nRows;

  SelfAdjointEigenSolver<MatrixXd> es;
  es.compute(M);

  double Eigenvalues[nSize];
  double Derivatives[nSize];
  for(int i=0;i<nSize;i++) Eigenvalues[i] = es.eigenvalues()[i];
  for(int i=0;i<nSize-1;i++)
    {
      if(std::abs(Eigenvalues[i]-Eigenvalues[i+1]) < EVThres) {
	  std::cerr << "ERROR ! repeated eigenvalues. \n";
      }
    }

  double Deriv[nSize][4];
  VectorXd v(nSize);
  MatrixXd C(nSize,nSize),E(nSize,nSize),Identity(nSize,nSize);
  Identity = MatrixXd::Identity(nSize,nSize);
  VectorXd vDiffX(nSize),vDiffY(nSize);

  for(int i=0;i<nSize;i++)
    {
      Deriv[i][0] = Eigenvalues[i];
      v = es.eigenvectors().col(i);
      Deriv[i][1] = v.transpose()*MDiffX*v;
      Deriv[i][2] = v.transpose()*MDiffY*v;

      C = (M-Deriv[i][0]*Identity).transpose()*(M-Deriv[i][0]*Identity) + v*v.transpose();
      E = (M-Deriv[i][0]*Identity).transpose()*(MDiffX-Deriv[i][1]*Identity);

      vDiffX = C.colPivHouseholderQr().solve(-E*v);

      E = (M-Deriv[i][0]*Identity).transpose()*(MDiffY-Deriv[i][2]*Identity);
      vDiffY = C.colPivHouseholderQr().solve(-E*v);

      Deriv[i][3] = v.transpose()*MDiffXY*v;
      Deriv[i][3] += v.transpose()*(MDiffX-Deriv[i][1]*Identity)*vDiffY;
      Deriv[i][3] += v.transpose()*(MDiffY-Deriv[i][2]*Identity)*vDiffX;

    }

  for(int i=0;i<nSize;i++)
    {
      for(int j=0;j<4;j++) res.push_back(Deriv[i][j]);
    }
}

/**
    Calculates all triple and quartic couplings in the physical basis
 */
void Class_Potential_Origin::CalculatePhysicalCouplings()
{
    bool Debug = false;
    if(!SetCurvatureDone) SetCurvatureArrays();
    if(Debug) std::cout << "Debug turned on in " << __func__ << std::endl;
    if(Debug) std::cout << "NH = " << NHiggs << "\tnG = " << NGauge << "\tnQ = " << NQuarks << "\tnL = " << NLepton << std::endl;
    const double ZeroMass = std::pow(10,-5);
    MatrixXd MassHiggs(NHiggs,NHiggs),MassGauge(NGauge,NGauge);
    MatrixXcd MassQuark(NQuarks,NQuarks),MassLepton(NLepton,NLepton);
    MassHiggs = MatrixXd::Zero(NHiggs,NHiggs);
    MassGauge = MatrixXd::Zero(NGauge,NGauge);
    MassQuark = MatrixXcd::Zero(NQuarks,NQuarks);
    MassLepton = MatrixXcd::Zero(NLepton,NLepton);

    MassSquaredGauge.resize(NGauge);
    MassSquaredHiggs.resize(NHiggs);
    MassSquaredQuark.resize(NQuarks);
    MassSquaredLepton.resize(NLepton);
    HiggsRotationMatrix.resize(NHiggs);
    for(int i=0;i<NHiggs;i++) HiggsRotationMatrix[i].resize(NHiggs);

    if(Debug) std::cout << "Setup done" << std::endl;


    for(int i=0;i<NHiggs;i++)
    {
        for(int j=0;j<NHiggs;j++)
        {
            MassHiggs(i,j) += Curvature_Higgs_L2[i][j];
            for(int k=0;k<NHiggs;k++)
            {
                MassHiggs(i,j) += Curvature_Higgs_L3[i][j][k]*HiggsVev[k];
                for(int l=0;l<NHiggs;l++) MassHiggs(i,j) += 0.5*Curvature_Higgs_L4[i][j][k][l]*HiggsVev[k]*HiggsVev[l];
            }
        }
    }

    for(int a=0;a<NGauge;a++)
    {
        for(int b=0;b<NGauge;b++)
        {
            for(int i=0;i<NHiggs;i++)
            {
                for(int j=0;j<NHiggs;j++)
                {
                    MassGauge(a,b) += 0.5*Curvature_Gauge_G2H2[a][b][i][j]*HiggsVev[i]*HiggsVev[j];
                }
            }
        }
    }

    if(Debug) std::cout << "GaugeMassMatrix done " << std::endl;

    MatrixXcd MIJQuarks(NQuarks,NQuarks);
    MIJQuarks = MatrixXcd::Zero(NQuarks,NQuarks);
    for(int a=0;a<NQuarks;a++)
    {
        for(int b=0;b<NQuarks;b++)
        {
            // MIJQuarks(a,b) = 0;
            for(int k=0;k<NHiggs;k++) MIJQuarks(a,b) += Curvature_Quark_F2H1[a][b][k]*HiggsVev[k];
        }
    }

    MassQuark = MIJQuarks.conjugate()*MIJQuarks;

    if(Debug) std::cout << "QuarkMassMatrix done " << std::endl;

    MatrixXcd MIJLeptons(NLepton,NLepton);
    MIJLeptons=MatrixXcd::Zero(NLepton,NLepton);
    for(int a=0;a<NLepton;a++)
    {
        for(int b=0;b<NLepton;b++)
        {
            // MIJLeptons(a,b) = 0;
            for(int k=0;k<NHiggs;k++) MIJLeptons(a,b) += Curvature_Lepton_F2H1[a][b][k]*HiggsVev[k];
        }
    }

    MassLepton = MIJLeptons.conjugate()*MIJLeptons;

    if(Debug) {
    	std::cout << "LeptonMassMatrix done " << std::endl;
    	std::cout << "MassGauge =  \n" << MassGauge << std::endl;
    }



    if(Debug) std::cout << "Set" << std::endl;


    MatrixXd HiggsRot(NHiggs,NHiggs),GaugeRot(NGauge,NGauge),QuarkRot(NQuarks,NQuarks),LepRot(NLepton,NLepton);
    HiggsRot = MatrixXd::Identity(NHiggs,NHiggs);
    GaugeRot = MatrixXd::Identity(NGauge,NGauge);
    QuarkRot = MatrixXd::Identity(NQuarks,NQuarks);
    LepRot = MatrixXd::Identity(NLepton,NLepton);

    SelfAdjointEigenSolver<MatrixXd> es;

	es.compute(MassHiggs);
	HiggsRot = es.eigenvectors().transpose();
	for(int i=0;i<NHiggs;i++)
	{
		for(int j=0;j<NHiggs;j++)
		{
			if(std::abs(HiggsRot(i,j)) < std::pow(10,-10)) HiggsRot(i,j) = 0;
		}
	}

    for(int i=0;i<NHiggs;i++)
    {
        MassSquaredHiggs[i] = es.eigenvalues()[i];
        if(std::abs(MassSquaredHiggs[i]) < ZeroMass) MassSquaredHiggs[i] = 0;
    }


	es.compute(MassGauge);
	GaugeRot = es.eigenvectors().transpose();

    for(int i=0;i<NGauge;i++)
    {
        MassSquaredGauge[i] = es.eigenvalues()[i];
        if(std::abs(MassSquaredGauge[i]) < ZeroMass) MassSquaredGauge[i] = 0;
    }


	SelfAdjointEigenSolver<MatrixXcd> esQuark(MassQuark);
	QuarkRot = esQuark.eigenvectors().transpose().real();
	for(int i=0;i<NQuarks;i++) MassSquaredQuark[i] = esQuark.eigenvalues().real()[i];

	SelfAdjointEigenSolver<MatrixXcd> esLepton(MassLepton);
	LepRot = esLepton.eigenvectors().transpose().real();
	for(int i=0;i<NLepton;i++) MassSquaredLepton[i] = esLepton.eigenvalues().real()[i];

    if(Debug) std::cout << "Calculated Masses and Rotation matrices " << std::endl;
    if(Debug) {
     std::cout << "Gauge Masses";
     for(int i=0;i<NGauge;i++) std::cout << "\t" << MassSquaredGauge[i];
     std::cout << std::endl;
     }





    // Higgs kopplungen aus tensoren zusammensetzen

    for(int a=0;a<NGauge;a++)
    {
        for(int b=0;b<NGauge;b++)
        {
            for(int i=0;i<NHiggs;i++)
            {
                LambdaGauge_3[a][b][i] = 0;
                for(int j=0;j<NHiggs;j++) LambdaGauge_3[a][b][i] += Curvature_Gauge_G2H2[a][b][i][j]*HiggsVev[j];
            }
        }
    }
    if(Debug){
    	std::cout << "LambdaGauge_3 done" << std::endl;
    }
    for(int a=0;a<NHiggs;a++)
    {
        for(int b=0;b<NHiggs;b++)
        {
            for(int i=0;i<NHiggs;i++)
            {
                LambdaHiggs_3[a][b][i] = Curvature_Higgs_L3[a][b][i];

                for(int j=0;j<NHiggs;j++) {
                	LambdaHiggs_3[a][b][i] += Curvature_Higgs_L4[a][b][i][j]*HiggsVev[j];

                }
            }
        }
    }
    if(Debug){
        	std::cout << "LambdaHiggs_3 done" << std::endl;
        }

    for(int i=0;i<NQuarks;i++)
    {
        for(int j=0;j<NQuarks;j++)
        {
            for(int k=0;k<NHiggs;k++)
            {
                LambdaQuark_3[i][j][k] = 0;
                for(int l=0;l<NQuarks;l++){
                    LambdaQuark_3[i][j][k] += conj(Curvature_Quark_F2H1[i][l][k])*MIJQuarks(l,j);
                    LambdaQuark_3[i][j][k] += conj(MIJQuarks(i,l)) * Curvature_Quark_F2H1[l][j][k];
                }
                for(int m=0;m<NHiggs;m++)
                {
                    LambdaQuark_4[i][j][k][m] = 0;
                    for(int l=0;l<NQuarks;l++)
                    {
                        LambdaQuark_4[i][j][k][m] += conj(Curvature_Quark_F2H1[i][l][k])*Curvature_Quark_F2H1[l][j][m];
                        LambdaQuark_4[i][j][k][m] += conj(Curvature_Quark_F2H1[i][l][m])*Curvature_Quark_F2H1[l][j][k];
                    }
                }
            }
        }
    }
    if(Debug){
            	std::cout << "LambdaQuark_4 done" << std::endl;
            }


    for(int i=0;i<NLepton;i++)
    {
        for(int j=0;j<NLepton;j++)
        {
            for(int k=0;k<NHiggs;k++)
            {
                LambdaLepton_3[i][j][k] = 0;
                for(int l=0;l<NLepton;l++){
                    LambdaLepton_3[i][j][k] += conj(Curvature_Lepton_F2H1[i][l][k])*MIJLeptons(l,j);
                    LambdaLepton_3[i][j][k] += conj(MIJLeptons(i,l)) * Curvature_Lepton_F2H1[l][j][k];
                }
                for(int m=0;m<NHiggs;m++)
                {
                    LambdaLepton_4[i][j][k][m] = 0;
                    for(int l=0;l<NLepton;l++)
                    {
                        LambdaLepton_4[i][j][k][m] += conj(Curvature_Lepton_F2H1[i][l][k])*Curvature_Lepton_F2H1[l][j][m];
                        LambdaLepton_4[i][j][k][m] += conj(Curvature_Lepton_F2H1[i][l][m])*Curvature_Lepton_F2H1[l][j][k];
                    }
                }
            }
        }
    }
    if(Debug){
		std::cout << "LambdaLepton_4 done" << std::endl;
	}

    // Rotate and save into corresponding vectors

    Couplings_Higgs_Quartic.resize(NHiggs);
    Couplings_Higgs_Triple.resize(NHiggs);
    for(int i=0;i<NHiggs;i++) {
        Couplings_Higgs_Quartic[i].resize(NHiggs);
        Couplings_Higgs_Triple[i].resize(NHiggs);
        for(int j=0;j<NHiggs;j++){
            Couplings_Higgs_Quartic[i][j].resize(NHiggs);
            Couplings_Higgs_Triple[i][j].resize(NHiggs);
            for(int k=0;k<NHiggs;k++) Couplings_Higgs_Quartic[i][j][k].resize(NHiggs);
        }
    }



    Couplings_Gauge_Higgs_22.resize(NGauge);
    Couplings_Gauge_Higgs_21.resize(NGauge);
    for(int a=0;a<NGauge;a++)
    {
        Couplings_Gauge_Higgs_22[a].resize(NGauge);
        Couplings_Gauge_Higgs_21[a].resize(NGauge);
        for(int b=0;b<NGauge;b++)
        {
            Couplings_Gauge_Higgs_22[a][b].resize(NHiggs);
            Couplings_Gauge_Higgs_21[a][b].resize(NHiggs);
            for(int i=0;i<NHiggs;i++)
            {
                Couplings_Gauge_Higgs_22[a][b][i].resize(NHiggs);
            }
        }
    }



    Couplings_Quark_Higgs_22.resize(NQuarks);
    Couplings_Quark_Higgs_21.resize(NQuarks);
    for(int a=0;a<NQuarks;a++)
    {
        Couplings_Quark_Higgs_22[a].resize(NQuarks);
        Couplings_Quark_Higgs_21[a].resize(NQuarks);
        for(int b=0;b<NQuarks;b++) {
            Couplings_Quark_Higgs_22[a][b].resize(NHiggs);
            Couplings_Quark_Higgs_21[a][b].resize(NHiggs);
            for(int i=0;i<NHiggs;i++) Couplings_Quark_Higgs_22[a][b][i].resize(NHiggs);
        }
    }




    Couplings_Lepton_Higgs_22.resize(NLepton);
    Couplings_Lepton_Higgs_21.resize(NLepton);
    for(int a=0;a<NLepton;a++)
    {
        Couplings_Lepton_Higgs_22[a].resize(NLepton);
        Couplings_Lepton_Higgs_21[a].resize(NLepton);
        for(int b=0;b<NLepton;b++)
        {
            Couplings_Lepton_Higgs_22[a][b].resize(NHiggs);
            Couplings_Lepton_Higgs_21[a][b].resize(NHiggs);
            for(int i=0;i<NHiggs;i++) Couplings_Lepton_Higgs_22[a][b][i].resize(NHiggs);
        }
    }

    if(Debug){
    	std::cout << "Couplings_Lepton_Higgs_ done" << std::endl;
    }



  for(int i=0;i<NHiggs;i++)
	{
		for(int j=0;j<NHiggs;j++)
		{
			for(int k=0;k<NHiggs;k++)
			{
				Couplings_Higgs_Triple[i][j][k] = 0;
				for(int is=0;is<NHiggs;is++)
				{
					for(int js=0;js<NHiggs;js++)
					{
						for(int ks=0;ks<NHiggs;ks++)
						{
							Couplings_Higgs_Triple[i][j][k] += HiggsRot(i,is)*HiggsRot(j,js)*HiggsRot(k,ks)*LambdaHiggs_3[is][js][ks];
						}
					}
				}
				for(int l=0;l<NHiggs;l++)
				{
					Couplings_Higgs_Quartic[i][j][k][l] = 0;
					for(int is=0;is<NHiggs;is++)
					{
						for(int js=0;js<NHiggs;js++)
						{
							for(int ks=0;ks<NHiggs;ks++)
							{
								for(int ls=0;ls<NHiggs;ls++)
								{
									Couplings_Higgs_Quartic[i][j][k][l] += HiggsRot(i,is)*HiggsRot(j,js)*HiggsRot(k,ks)*HiggsRot(l,ls)*Curvature_Higgs_L4[is][js][ks][ls];
								}
							}
						}
					}
				}
			}
		}
	}




    // Gauge Rot
	for(int a=0;a<NGauge;a++)
	{
		for(int b=0;b<NGauge;b++)
		{
			for(int i=0;i<NHiggs;i++)
			{
				Couplings_Gauge_Higgs_21[a][b][i] = 0;
				for(int as=0;as<NGauge;as++)
				{
					for(int bs=0;bs<NGauge;bs++)
					{
						for(int is=0;is<NHiggs;is++)
						Couplings_Gauge_Higgs_21[a][b][i] += GaugeRot(a,as)*GaugeRot(b,bs)*HiggsRot(i,is)*LambdaGauge_3[as][bs][is];
					}
				}
				for(int j=0;j<NHiggs;j++)
				{
					Couplings_Gauge_Higgs_22[a][b][i][j] = 0;
					for(int as=0;as<NGauge;as++)
					{
						for(int bs=0;bs<NGauge;bs++)
						{
							for(int is=0;is<NHiggs;is++)
							{
								for(int js=0;js<NHiggs;js++){
									double RotFac = GaugeRot(a,as)*GaugeRot(b,bs)*HiggsRot(i,is)*HiggsRot(j,js);
									Couplings_Gauge_Higgs_22[a][b][i][j] += RotFac*Curvature_Gauge_G2H2[as][bs][is][js];
								}
							}

						}
					}
				}
			}
		}
	}


    // Quark

	for(int a=0;a<NQuarks;a++)
	{
		for(int b=0;b<NQuarks;b++)
		{
			for(int i=0;i<NHiggs;i++)
			{
				Couplings_Quark_Higgs_21[a][b][i] = 0;
				for(int as = 0;as<NQuarks;as++)
				{
					for(int bs=0;bs<NQuarks;bs++)
					{
						for(int is=0;is<NHiggs;is++)
						{
							double RotFac = QuarkRot(a,as)*QuarkRot(b,bs)*HiggsRot(i,is);
							Couplings_Quark_Higgs_21[a][b][i] += RotFac*LambdaQuark_3[as][bs][is];
						}
					}
				}
				for(int j=0;j<NHiggs;j++)
				{
					Couplings_Quark_Higgs_22[a][b][i][j] = 0;
					for(int as=0;as<NQuarks;as++)
					{
						for(int bs=0;bs<NQuarks;bs++)
						{
							for(int is=0;is<NHiggs;is++)
							{
								for(int js=0;js<NHiggs;js++)
								{
									double RotFac = QuarkRot(a,as)*QuarkRot(b,bs)*HiggsRot(i,is)*HiggsRot(j,js);
									Couplings_Quark_Higgs_22[a][b][i][j] += RotFac*LambdaQuark_4[as][bs][is][js];
								}
							}
						}
					}
				}
			}
		}
	}


    // Lepton


	for(int a=0;a<NLepton;a++)
	{
		for(int b=0;b<NLepton;b++)
		{
			for(int i=0;i<NHiggs;i++)
			{
				Couplings_Lepton_Higgs_21[a][b][i] = 0;
				for(int as = 0;as<NLepton;as++)
				{
					for(int bs=0;bs<NLepton;bs++)
					{
						for(int is=0;is<NHiggs;is++)
						{
							double RotFac = LepRot(a,as)*LepRot(b,bs)*HiggsRot(i,is);
							Couplings_Lepton_Higgs_21[a][b][i] += RotFac*LambdaLepton_3[as][bs][is];
						}
					}
				}
				for(int j=0;j<NHiggs;j++)
				{
					Couplings_Lepton_Higgs_22[a][b][i][j] = 0;
					for(int as=0;as<NLepton;as++)
					{
						for(int bs=0;bs<NLepton;bs++)
						{
							for(int is=0;is<NHiggs;is++)
							{
								for(int js=0;js<NHiggs;js++)
								{
									double RotFac = LepRot(a,as)*LepRot(b,bs)*HiggsRot(i,is)*HiggsRot(j,js);
									Couplings_Lepton_Higgs_22[a][b][i][j] += RotFac*LambdaLepton_4[as][bs][is][js];
								}
							}
						}
					}
				}
			}
		}
	}

    for(int i=0;i<NHiggs;i++)
    {
        for(int j=0;j<NHiggs;j++)
        {
            HiggsRotationMatrix[i][j] = HiggsRot(i,j);
        }
    }


    CalcCouplingsdone = true;

    if(Debug)
    {
        std::cout << "HiggsRot = \n" << HiggsRot << "\n" << std::endl;
        std::cout << "Higgsmassen : \n";
        for(int i=0;i<NHiggs;i++) std::cout << MassSquaredHiggs[i] << "\t";
        std::cout << std::endl;
        std::cout << "GaugeRot = \n" << GaugeRot << "\n" << std::endl;
        std::cout << "LepRot = \n" << LepRot << "\n" << std::endl;
        std::cout << "QuarkRot = \n" << QuarkRot << "\n" << std::endl;


        std::cout << "HiggsMassMatrix \n" << MassHiggs << std::endl;


        std::cout << "HiggsVev = ";
        for(int i=0;i<NHiggs;i++) std::cout << "\t" << HiggsVev[i];
        std::cout << std::endl;
    }

    if(Debug) std::cout << "\nEnd of Debug output of " << __func__ << std::endl;

    return;


}

/**
 * Calculates the first derivative of the Coleman-Weinberg potential.
 */
void Class_Potential_Origin::WeinbergFirstDerivative(std::vector<double>& res){
    bool Debug = false;
    if(!CalcCouplingsdone) CalculatePhysicalCouplings();
    const double NumZero = std::pow(10,-10);
    if(Debug) std::cout << "Debug turned on in" << __func__ << std::endl;
    VectorXd FirstDeriv(NHiggs),FirstDerivGauge(NHiggs),FirstDerivHiggs(NHiggs),FirstDerivQuark(NHiggs),FirstDerivLepton(NHiggs);
    FirstDeriv = VectorXd::Zero(NHiggs);
    FirstDerivGauge = VectorXd::Zero(NHiggs);
    FirstDerivHiggs = VectorXd::Zero(NHiggs);
    FirstDerivQuark = VectorXd::Zero(NHiggs);
    FirstDerivLepton = VectorXd::Zero(NHiggs);
    double epsilon = 1.0/(16*M_PI*M_PI);

    for(int i=0;i<NHiggs;i++)
    {
        for(int a=0;a<NGauge;a++)
        {
            if(MassSquaredGauge[a] != 0)
            {
                FirstDerivGauge[i] += MassSquaredGauge[a]*Couplings_Gauge_Higgs_21[a][a][i]*(std::log(MassSquaredGauge[a]/std::pow(scale,2))-C_CWcbGB + 0.5 );
            }
        }

        for(int a = 0;a<NHiggs;a++)
        {
            if(MassSquaredHiggs[a] != 0)
            {

                FirstDerivHiggs[i] += MassSquaredHiggs[a]*Couplings_Higgs_Triple[a][a][i]*(std::log(MassSquaredHiggs[a]/std::pow(scale,2)) - C_CWcbHiggs + 0.5);
            }
        }
        for(int a=0;a<NQuarks;a++)
        {
            if(MassSquaredQuark[a] != 0)
            {
                double Coup = Couplings_Quark_Higgs_21[a][a][i].real();
                FirstDerivQuark[i] += MassSquaredQuark[a]*Coup*(std::log(MassSquaredQuark[a]/std::pow(scale,2)) - C_CWcbFermion + 0.5);
            }
        }
        for(int a=0;a<NLepton;a++)
        {
            if(MassSquaredLepton[a] != 0)
            {
                double Coup = Couplings_Lepton_Higgs_21[a][a][i].real();
                FirstDerivLepton[i] += MassSquaredLepton[a]*Coup*(std::log(MassSquaredLepton[a]/std::pow(scale,2)) - C_CWcbFermion + 0.5);
            }
        }
    }
    FirstDerivGauge *= 1.5;
    FirstDerivHiggs *= 0.5;
    FirstDerivQuark *= -3;
    FirstDerivLepton *= -1;

    MatrixXd HiggsRot(NHiggs,NHiggs);
    for(int i=0;i<NHiggs;i++)
    {
        for(int j=0;j<NHiggs;j++) HiggsRot(i,j) = HiggsRotationMatrix[i][j];
    }

    FirstDeriv = HiggsRot.transpose()*(FirstDerivGauge+FirstDerivHiggs+FirstDerivQuark+FirstDerivLepton);
    FirstDeriv *= epsilon;

    for(int i=0;i<NHiggs;i++)
    {
        if(std::abs(FirstDeriv[i]) < NumZero) FirstDeriv[i] = 0;
    }

    for(int i=0;i<NHiggs;i++)    res.push_back(FirstDeriv[i]);

    if(Debug) {
        std::cout << "WeinbergFirstDerivative = \t" << FirstDeriv.transpose() << "\n" << std::endl;
    }



}

/**
 * Calculates the second derivative of the Coleman-Weinberg potential.
 */
void Class_Potential_Origin::WeinbergSecondDerivative(std::vector<double>& res)
{
    bool Debug = false;
    if(!CalcCouplingsdone) CalculatePhysicalCouplings();
    const double NumZero = std::pow(10,-10);
    if(Debug) std::cout << "Debug turned on in " << __func__ << std::endl;
    MatrixXd GaugePart(NHiggs,NHiggs),HiggsPart(NHiggs,NHiggs),QuarkPart(NHiggs,NHiggs),LeptonPart(NHiggs,NHiggs);
    GaugePart = MatrixXd::Zero(NHiggs,NHiggs);
    HiggsPart = MatrixXd::Zero(NHiggs,NHiggs);
    QuarkPart = MatrixXd::Zero(NHiggs,NHiggs);
    LeptonPart = MatrixXd::Zero(NHiggs,NHiggs);

    for(int i=0;i<NHiggs;i++)
    {
        for(int j=0;j<NHiggs;j++)
        {
            for(int a=0;a<NGauge;a++)
            {
                for(int b=0;b<NGauge;b++)
                {
                    double Coup1 = Couplings_Gauge_Higgs_21[a][b][i];
                    double Coup2 = Couplings_Gauge_Higgs_21[b][a][j];
                    double Br = fbase(MassSquaredGauge[a],MassSquaredGauge[b]) - C_CWcbGB + 0.5;
                    GaugePart(i,j) += Coup1*Coup2*Br;
                }
                if(MassSquaredGauge[a] != 0)
                {
                    GaugePart(i,j) += MassSquaredGauge[a]*Couplings_Gauge_Higgs_22[a][a][i][j]*(std::log(MassSquaredGauge[a]/std::pow(scale,2)) - C_CWcbGB + 0.5);
                }
            }


            for(int a=0;a<NHiggs;a++)
            {
                for(int b=0;b<NHiggs;b++)
                {
                    double Coup1 = Couplings_Higgs_Triple[a][b][i];
                    double Coup2 = Couplings_Higgs_Triple[b][a][j];
                    double Br = fbase(MassSquaredHiggs[a],MassSquaredHiggs[b]) - C_CWcbHiggs + 0.5;
                    HiggsPart(i,j) += Coup1*Coup2*Br;
                }
                if(MassSquaredHiggs[a] != 0)
                {
                    HiggsPart(i,j) += MassSquaredHiggs[a]*Couplings_Higgs_Quartic[a][a][i][j]*(std::log(MassSquaredHiggs[a]/std::pow(scale,2)) - C_CWcbHiggs + 0.5);
                }
            }

            for(int a=0;a<NQuarks;a++)
            {
                for(int b=0;b<NQuarks;b++)
                {
                    double Coup = (Couplings_Quark_Higgs_21[a][b][i] * Couplings_Quark_Higgs_21[b][a][j]).real();
                    double Br = fbase(MassSquaredQuark[a],MassSquaredQuark[b]) - C_CWcbFermion + 0.5;
                    QuarkPart(i,j) += Coup*Br;
                }
                if(MassSquaredQuark[a] != 0)
                {
                    double Coup = Couplings_Quark_Higgs_22[a][a][i][j].real();
                    QuarkPart(i,j) += MassSquaredQuark[a]*Coup*(std::log(MassSquaredQuark[a]/std::pow(scale,2)) - C_CWcbFermion + 0.5);

                }
            }

            for(int a=0;a<NLepton;a++)
            {
                for(int b=0;b<NLepton;b++)
                {
                    double Coup = (Couplings_Lepton_Higgs_21[a][b][i]*Couplings_Lepton_Higgs_21[b][a][j]).real();
                    double Br = fbase(MassSquaredLepton[a],MassSquaredLepton[b]) - C_CWcbFermion + 0.5;
                    LeptonPart(i,j) += Coup*Br;
                }
                if(MassSquaredLepton[a] != 0)
                {
                    double Coup = Couplings_Lepton_Higgs_22[a][a][i][j].real();
                    LeptonPart(i,j) += Coup*MassSquaredLepton[a]*(std::log(MassSquaredLepton[a]/std::pow(scale,2)) - C_CWcbFermion + 0.5);
                }
            }
        }
    }

    HiggsPart *= 0.5;
    GaugePart *= 1.5;
    QuarkPart*= -3;
    LeptonPart *= -1;

    MatrixXd Storage(NHiggs,NHiggs);
    Storage = HiggsPart+GaugePart+QuarkPart+LeptonPart;
    MatrixXd ResMatrix;
    MatrixXd HiggsRot(NHiggs,NHiggs);
    for(int i=0;i<NHiggs;i++)
    {
        for(int j=0;j<NHiggs;j++)
        {
            HiggsRot(i,j) = HiggsRotationMatrix[i][j];
        }
    }

    ResMatrix = 0.5*HiggsRot.transpose()*(Storage+Storage.transpose())*HiggsRot;
    double epsilon = 1.0/(16.0*M_PI*M_PI);
    ResMatrix *= epsilon;

    for(int i=0;i<NHiggs;i++)
    {
        for(int j=0;j<NHiggs;j++)
        {
            if(std::abs(ResMatrix(i,j)) < NumZero) ResMatrix(i,j) = 0;
        }
    }

    for(int i=0;i<NHiggs;i++)
    {
        for(int j=0;j<NHiggs;j++)
        {
            res.push_back(ResMatrix(j,i));
        }
    }




}

/**
 * Calculates the third derivative of the Coleman-Weinberg potential.
 */
void Class_Potential_Origin::WeinbergThirdDerivative(std::vector<double>& res){

	bool Debug = false;
	if(Debug) std::cout << "Debug turned on in " << __func__ << std::endl;
  const double NumZero = std::pow(10,-10);
  double epsilon = 1.0/(16.0*M_PI*M_PI);
//  std::complex<double> restmp[NHiggs][NHiggs][NHiggs];
//  std::complex<double> QuarkPart[NHiggs][NHiggs][NHiggs];
//  std::complex<double> LeptonPart[NHiggs][NHiggs][NHiggs];
//  std::complex<double> QuarkPartSym[NHiggs][NHiggs][NHiggs];
//  std::complex<double> LeptonPartSym[NHiggs][NHiggs][NHiggs];

  std::vector<std::vector<std::vector< std::complex<double>>>> restmp;
  std::vector<std::vector<std::vector< std::complex<double>>>> QuarkPart;
  std::vector<std::vector<std::vector< std::complex<double>>>> LeptonPart;
  std::vector<std::vector<std::vector< std::complex<double>>>> QuarkPartSym;
  std::vector<std::vector<std::vector< std::complex<double>>>> LeptonPartSym;
  restmp.resize(NHiggs);
  QuarkPart.resize(NHiggs);
  LeptonPart.resize(NHiggs);
  QuarkPartSym.resize(NHiggs);
  LeptonPartSym.resize(NHiggs);
  if(Debug) std::cout << "Setup 1D " << std::endl;
  for(int i=0;i<NHiggs;i++)
  {
	  restmp[i].resize(NHiggs);
	  QuarkPart[i].resize(NHiggs);
	  LeptonPart[i].resize(NHiggs);
	  QuarkPartSym[i].resize(NHiggs);
	  LeptonPartSym[i].resize(NHiggs);
	  if(Debug) std::cout << "Setup 2D " << std::endl;
	  for(int j=0;j<NHiggs;j++) {
		  restmp[i][j].resize(NHiggs);
		  QuarkPart[i][j].resize(NHiggs );
		  LeptonPart[i][j].resize(NHiggs);
		  QuarkPartSym[i][j].resize(NHiggs);
		  LeptonPartSym[i][j].resize(NHiggs);
		  if(Debug) std::cout << "Setup 3D " << std::endl;
	  }
  }

  double resGaugeBase[NHiggs][NHiggs][NHiggs];
  double Higgspart[NHiggs][NHiggs][NHiggs];


  double GaugePart[NHiggs][NHiggs][NHiggs];

  double HiggspartSym[NHiggs][NHiggs][NHiggs];

  double GaugePartSym[NHiggs][NHiggs][NHiggs];

  if(Debug) std::cout << "Setup done " << std::endl;


  for(int i=0;i<NHiggs;i++)
    {
      for(int j=0;j<NHiggs;j++)
	{
	  for(int k=0;k<NHiggs;k++)
	    {
	      Higgspart[i][j][k] = 0;
	      for(int a=0;a<NHiggs;a++)
		{
		  for(int b=0;b<NHiggs;b++)
		    {
		      for(int c=0;c<NHiggs;c++)
			{
			  double f1 = fbaseTri(MassSquaredHiggs[a],MassSquaredHiggs[b],MassSquaredHiggs[c]);
			  double f2 = Couplings_Higgs_Triple[a][b][i];
			  double f3 = Couplings_Higgs_Triple[b][c][j];
			  double f4 = Couplings_Higgs_Triple[c][a][k];
			  Higgspart[i][j][k] += 2*f1*f2*f3*f4;
			}
		      double f1 = Couplings_Higgs_Quartic[a][b][i][j];
		      double f2 = Couplings_Higgs_Triple[b][a][k];
		      double f3 = fbase(MassSquaredHiggs[a],MassSquaredHiggs[b]) - C_CWcbHiggs + 0.5;
		      Higgspart[i][j][k] += 3*f1*f2*f3;
		    }
		}

	      GaugePart[i][j][k] = 0;
	      for(int a=0;a<NGauge;a++)
		{
		  for(int b=0;b<NGauge;b++)
		    {
		      for(int c=0;c<NGauge;c++)
			{
			  double f1 = fbaseTri(MassSquaredGauge[a],MassSquaredGauge[b],MassSquaredGauge[c]);
			  double f2 = Couplings_Gauge_Higgs_21[a][b][i];
			  double f3 = Couplings_Gauge_Higgs_21[b][c][j];
			  double f4 = Couplings_Gauge_Higgs_21[c][a][k];
			  GaugePart[i][j][k] += 2*f1*f2*f3*f4;
			}
		      double f1 = Couplings_Gauge_Higgs_22[a][b][i][j];
		      double f2 = Couplings_Gauge_Higgs_21[b][a][k];
		      double f3 = fbase(MassSquaredGauge[a],MassSquaredGauge[b]) - C_CWcbGB + 0.5;
		      GaugePart[i][j][k] += 3*f1*f2*f3;
		    }
		}

	      QuarkPart[i][j][k] = 0;
	      for(int a=0;a<NQuarks;a++)
		{
		  for(int b=0;b<NQuarks;b++)
		    {
		      for(int c=0;c<NQuarks;c++)
			{
			  std::complex<double> f1 = fbaseTri(MassSquaredQuark[a],MassSquaredQuark[b],MassSquaredQuark[c]);
			  std::complex<double> f2 = Couplings_Quark_Higgs_21[a][b][i];
			  std::complex<double> f3 = Couplings_Quark_Higgs_21[b][c][j];
			  std::complex<double> f4 = Couplings_Quark_Higgs_21[c][a][k];
			  QuarkPart[i][j][k] += 2.0*f1*f2*f3*f4;
			}
		      std::complex<double> f1 = Couplings_Quark_Higgs_22[a][b][i][j];
		      std::complex<double> f2 = Couplings_Quark_Higgs_21[b][a][k];
		      std::complex<double> f3 = fbase(MassSquaredQuark[a],MassSquaredQuark[b]) - C_CWcbFermion + 0.5;
		      QuarkPart[i][j][k] += 3.0*f1*f2*f3;
		    }
		}
	      LeptonPart[i][j][k] = 0;
	      for(int a=0;a<NLepton;a++)
		{
		  for(int b=0;b<NLepton;b++)
		    {
		      for(int c=0;c<NLepton;c++)
			{
			  std::complex<double> f1 = fbaseTri(MassSquaredLepton[a],MassSquaredLepton[b],MassSquaredLepton[c]);
			  std::complex<double> f2 = Couplings_Lepton_Higgs_21[a][b][i];
			  std::complex<double> f3 = Couplings_Lepton_Higgs_21[b][c][j];
			  std::complex<double> f4 = Couplings_Lepton_Higgs_21[c][a][k];
			  LeptonPart[i][j][k] += 2.0*f1*f2*f3*f4;
			}
		      std::complex<double> f1 = Couplings_Lepton_Higgs_22[a][b][i][j];
		      std::complex<double> f2 = Couplings_Lepton_Higgs_21[b][a][k];
		      std::complex<double> f3 = fbase(MassSquaredLepton[a],MassSquaredLepton[b]) - C_CWcbFermion + 0.5;
		      LeptonPart[i][j][k] += 3.0*f1*f2*f3;
		    }
		}
	    }
	}
    }

  for(int i=0;i<NHiggs;i++)
    {
      for(int j=0;j<NHiggs;j++)
	{
	  for(int k=0;k<NHiggs;k++)
	    {
	      HiggspartSym[i][j][k] = Higgspart[i][j][k] + Higgspart[i][k][j];
	      HiggspartSym[i][j][k] += Higgspart[j][i][k] + Higgspart[j][k][i];
	      HiggspartSym[i][j][k] += Higgspart[k][i][j] + Higgspart[k][j][i];
	      HiggspartSym[i][j][k] *= 1.0/6.0;

	      GaugePartSym[i][j][k] = GaugePart[i][j][k] + GaugePart[i][k][j];
	      GaugePartSym[i][j][k] += GaugePart[j][i][k] + GaugePart[j][k][i];
	      GaugePartSym[i][j][k] += GaugePart[k][i][j] + GaugePart[k][j][i];
	      GaugePartSym[i][j][k] *= 1.0/6.0;

	      QuarkPartSym[i][j][k] = QuarkPart[i][j][k] + QuarkPart[i][k][j];
	      QuarkPartSym[i][j][k] += QuarkPart[j][i][k] + QuarkPart[j][k][i];
	      QuarkPartSym[i][j][k] += QuarkPart[k][i][j] + QuarkPart[k][j][i];
	      QuarkPartSym[i][j][k] *= 1.0/6.0;

	      LeptonPartSym[i][j][k] = LeptonPart[i][j][k] + LeptonPart[i][k][j];
	      LeptonPartSym[i][j][k] += LeptonPart[j][i][k] + LeptonPart[j][k][i];
	      LeptonPartSym[i][j][k] += LeptonPart[k][i][j] + LeptonPart[k][j][i];
	      LeptonPartSym[i][j][k] *= 1.0/6.0;
	    }
	}
    }

  for(int i=0;i<NHiggs;i++)
  {
    for(int j=0;j<NHiggs;j++)
    {
      for(int k=0;k<NHiggs;k++)
	{
	  restmp[i][j][k] = 0.5*HiggspartSym[i][j][k];
	  restmp[i][j][k] += 1.5*GaugePartSym[i][j][k];
	  restmp[i][j][k] += -1.0*LeptonPartSym[i][j][k];
	  restmp[i][j][k] += -3.0*QuarkPartSym[i][j][k];
	}
    }
  }

  for(int l=0;l<NHiggs;l++)
    {
      for(int m=0;m<NHiggs;m++)
      {
        for(int n=0;n<NHiggs;n++)
  	{
  	  resGaugeBase[l][m][n] = 0;
  	  for(int i=0;i<NHiggs;i++)
  	    {
  	      for(int j=0;j<NHiggs;j++)
  		{
  		  for(int k=0;k<NHiggs;k++)
  		    {
  		      double RotFac = HiggsRotationMatrix[i][l]*HiggsRotationMatrix[j][m]*HiggsRotationMatrix[k][n];
  		      resGaugeBase[l][m][n] += RotFac*restmp[i][j][k].real();
  		    }
  		}
  	    }
  	  resGaugeBase[l][m][n] *= epsilon;
  	  if(std::abs(resGaugeBase[l][m][n]) < NumZero) resGaugeBase[l][m][n] = 0;
  	}
      }
    }

  for(int l=0;l<NHiggs;l++)
    {
      for(int m=0;m<NHiggs;m++)
	{
	  for(int n=0;n<NHiggs;n++)
	    {
	      res.push_back(resGaugeBase[l][m][n]);
	    }
	}
    }


}

/**
 * Calculates the Higgs mass matrix and saves all eigenvalues
 * @param res Vector in which the eigenvalues m^2 of the mass matrix will be stored
 * @param v the configuration of all VEVs at which the eigenvalues should be evaluated
 * @param Temp The temperature at which the Debye corrected masses should be calculated
 * @param diff 0 returns the masses and i!=0 returns the derivative of m^2 w.r.t v_i
 */
void Class_Potential_Origin::HiggsMassesSquared(std::vector<double>& res, const std::vector<double>& v, double Temp, int diff)
{
    if(!SetCurvatureDone) SetCurvatureArrays();
    MatrixXd MassMatrix(NHiggs,NHiggs);
    double ZeroMass = std::pow(10,-5);
    for(int i=0;i<NHiggs;i++)
    {
        for(int j=0;j<NHiggs;j++)
        {
            MassMatrix(i,j) = Curvature_Higgs_L2[i][j];
            for(int k=0;k<NHiggs;k++)
            {
                MassMatrix(i,j) += Curvature_Higgs_L3[i][j][k]*v[k];
                for(int l=0;l<NHiggs;l++) MassMatrix(i,j) += 0.5*Curvature_Higgs_L4[i][j][k][l] * v[k]*v[l];
            }

            if(Temp != 0)
            {
                MassMatrix(i,j) += DebyeHiggs[i][j]*std::pow(Temp,2);
            }
        }
    }

//    std::cout << MassMatrix << std::endl;

    if(diff == 0 and res.size() == 0)
    {
        SelfAdjointEigenSolver<MatrixXd> es(MassMatrix,EigenvaluesOnly);
        for(int i =0;i<NHiggs;i++)
        {
            double tmp = es.eigenvalues()[i];
            if(std::abs(tmp) < ZeroMass ) res.push_back(0);
            else res.push_back(tmp);
        }
    }
    else if(diff == 0 and res.size() == NHiggs){
    	SelfAdjointEigenSolver<MatrixXd> es(MassMatrix,EigenvaluesOnly);
		for(int i =0;i<NHiggs;i++)
		{
			double tmp = es.eigenvalues()[i];
			if(std::abs(tmp) < ZeroMass ) tmp = 0;
			res[i] = tmp;
		}
    }
    else if(diff == 0 and res.size()!= 0 and res.size() != NHiggs){
    	std::cout << "Something went wrong in " << __func__ << std::endl;
    	std::cout << __func__ << "Is calculating the mass for " << NHiggs << "fields but the resolution vector has a size of "
    			<< res.size() << ". This should be zero or " << NHiggs << std::endl;
    }
    else if(diff <= nVEV)
    {
        MatrixXd Diff(NHiggs,NHiggs);
        int x0 = diff -1;
        for(int i=0;i<NHiggs;i++)
        {
            for(int j=0;j<NHiggs;j++)
            {
                Diff(i,j) = Curvature_Higgs_L3[i][j][x0];
                for(int k=0;k<NHiggs;k++)
                {
                    Diff(i,j) += Curvature_Higgs_L4[i][j][x0][k]*v[k];
                }
            }
        }
        FirstDerivativeOfEigenvalues(MassMatrix,Diff,res);
    }


}

/**
 * Calculates the gauge mass matrix and saves all eigenvalues
 * @param res Vector in which the eigenvalues m^2 of the mass matrix will be stored
 * @param v the configuration of all VEVs at which the eigenvalues should be evaluated
 * @param Temp The temperature at which the Debye corrected masses should be calculated
 * @param diff 0 returns the masses and i!=0 returns the derivative of m^2 w.r.t v_i
 */
void Class_Potential_Origin::GaugeMassesSquared(std::vector<double>& res,const std::vector<double>& v, double Temp,int diff)
{
    if(!SetCurvatureDone) SetCurvatureArrays();
    MatrixXd MassMatrix(NGauge,NGauge);
    double ZeroMass = std::pow(10,-5);
    for(int i=0;i<NGauge;i++)
    {
        for(int j=0;j<NGauge;j++)
        {
            MassMatrix(i,j) = 0;
            for(int k=0;k<NHiggs;k++)
            {
                for(int l=0;l<NHiggs;l++) MassMatrix(i,j) += 0.5*Curvature_Gauge_G2H2[i][j][k][l] * v[k]*v[l];
            }
            if(Temp != 0)
            {
                MassMatrix(i,j) += DebyeGauge[i][j]*std::pow(Temp,2);
            }
        }
    }

    if(diff == 0)
    {

        SelfAdjointEigenSolver<MatrixXd> es(MassMatrix,EigenvaluesOnly);
        for(int i =0;i<NGauge;i++)
        {
            double tmp = es.eigenvalues()[i];
            if(std::abs(tmp) < ZeroMass ) res.push_back(0);
            else res.push_back(tmp);
        }
    }

    else if(diff <= nVEV)
    {
        int i = diff -1;
        MatrixXd Diff(NGauge,NGauge);
        Diff = MatrixXd::Zero(NGauge,NGauge);
        for(int a=0;a<NGauge;a++)
        {
            for(int b=0;b<NGauge;b++)
            {
                for(int j=0;j<NHiggs;j++) Diff(a,b) += Curvature_Gauge_G2H2[a][b][i][j]*v[j];
            }
        }
        FirstDerivativeOfEigenvalues(MassMatrix,Diff,res);
    }


}

/**
 * Calculates the quark mass matrix and saves all eigenvalues, this assumes the same masses
 *  for different colours.
 * @param res Vector in which the eigenvalues m^2 of the mass matrix will be stored
 * @param v the configuration of all VEVs at which the eigenvalues should be evaluated
 * @param Temp The temperature at which the Debye corrected masses should be calculated
 * @param diff 0 returns the masses and i!=0 returns the derivative of m^2 w.r.t v_i
 */
void Class_Potential_Origin::QuarkMassesSquared(std::vector<double>& res, const std::vector<double>& v, int diff)
{
    if(!SetCurvatureDone) {SetCurvatureArrays(); std::cout << "Reset of Set Curvature " << std::endl;};
    MatrixXcd MassMatrix(NQuarks,NQuarks),MIJ(NQuarks,NQuarks);
    MIJ = MatrixXcd::Zero(NQuarks,NQuarks);
    double ZeroMass = std::pow(10,-10);
    for(int i=0;i<NQuarks;i++)
    {
        for(int j=0;j<NQuarks;j++)
        {

            for(int k=0;k<NHiggs;k++)
            {
                MIJ(i,j) += Curvature_Quark_F2H1[i][j][k]*v[k];
            }
        }
    }

    MassMatrix = MIJ.conjugate()*MIJ;

    if(diff == 0)
    {
        SelfAdjointEigenSolver<MatrixXcd> es(MassMatrix,EigenvaluesOnly);
        for(int i =0;i<NQuarks;i++)
        {
            double tmp = es.eigenvalues().real()[i];
            if(std::abs(tmp) < ZeroMass ) res.push_back(0);
            else res.push_back(tmp);
        }
    }




}

/**
 * Calculates the lepton mass matrix and saves all eigenvalues
 * @param res Vector in which the eigenvalues m^2 of the mass matrix will be stored
 * @param v the configuration of all VEVs at which the eigenvalues should be evaluated
 * @param Temp The temperature at which the Debye corrected masses should be calculated
 * @param diff 0 returns the masses and i!=0 returns the derivative of m^2 w.r.t v_i
 */
void Class_Potential_Origin::LeptonMassesSquared(std::vector<double>& res, const std::vector<double>& v, int diff)
{
    if(!SetCurvatureDone) SetCurvatureArrays();
    MatrixXcd MassMatrix(NLepton,NLepton),MIJ(NLepton,NLepton);
    MIJ = MatrixXcd::Zero(NLepton,NLepton);
    double ZeroMass = std::pow(10,-10);
    for(int i=0;i<NLepton;i++)
    {
        for(int j=0;j<NLepton;j++)
        {

            for(int k=0;k<NHiggs;k++)
            {
                MIJ(i,j) += Curvature_Lepton_F2H1[i][j][k]*v[k];
            }
        }
    }

    MassMatrix = MIJ.conjugate()*MIJ;

    SelfAdjointEigenSolver<MatrixXcd> es(MassMatrix,EigenvaluesOnly);
    for(int i =0;i<NLepton;i++)
    {
        double tmp = es.eigenvalues().real()[i];
        if(std::abs(tmp) < ZeroMass ) res.push_back(0);
        else res.push_back(tmp);
    }


}

/**
 * Calculates the tree-level potential and its derivatives.
 * @param v the configuration of all VEVs at which the potential should be calculated
 * @param diff 0 returns the potential and i!= 0 returns the derivative of the potential w.r.t v_i
 */
double Class_Potential_Origin::VTree(const std::vector<double>& v, int diff)
{
    double res = VTreeSimplified(v) ;
    if(UseVTreeSimplified)
    {

    	return res;
    }
    res = 0;


    if(diff == 0)
    {
        for(int i=0;i<NHiggs;i++)
        {
        	if(v[i] != 0)
        	{
				res += Curvature_Higgs_L1[i]*v[i];
				for(int j=0;j<NHiggs;j++)
				{
					if(v[j] !=0)
					{
						res += 0.5*Curvature_Higgs_L2[i][j]*v[i]*v[j];
						for(int k=0;k<NHiggs;k++)
						{
							res += 1.0/6.0 * Curvature_Higgs_L3[i][j][k]*v[i]*v[j]*v[k];
							for(int l=0;l<NHiggs;l++)
							{
								res += 1.0/24.0*Curvature_Higgs_L4[i][j][k][l]*v[i]*v[j]*v[k]*v[l];
							}
						}
					}
				}
        	}
        }
    }
    else if(diff <= NHiggs)
    {
        int i = diff -1;
        res = Curvature_Higgs_L1[i];
        for(int j=0;j<NHiggs;j++)
        {
            res += Curvature_Higgs_L2[i][j]*v[j];
            for(int k=0;k<NHiggs;k++)
            {
                res += 0.5*Curvature_Higgs_L3[i][j][k]*v[j]*v[k];
                for(int l=0;l<NHiggs;l++)
                {
                    res += 1.0/6.0*Curvature_Higgs_L4[i][j][k][l]*v[j]*v[k]*v[l];
                }
            }
        }
    }

    return res;

}

/**
 * Calculates the counterterm potential and its derivatives
 * @param v the configuration of all VEVs at which the potential should be calculated
 * @param diff 0 returns the potential and i!= 0 returns the derivative of the potential w.r.t v_i
 */
double Class_Potential_Origin::CounterTerm(const std::vector<double>& v, int diff)
{
    double res = VCounterSimplified(v) ;
    if(UseVCounterSimplified) return res;
//    std::cout << "res Sim = " << res << std::endl;





    res =0;
    if(diff == 0)
    {
        for(int i=0;i<NHiggs;i++)
        {
            res += Curvature_Higgs_CT_L1[i]*v[i];
            for(int j=0;j<NHiggs;j++)
            {
                res += 0.5*Curvature_Higgs_CT_L2[i][j]*v[i]*v[j];
                for(int k=0;k<NHiggs;k++)
                {
                    res += 1.0/6.0 * Curvature_Higgs_CT_L3[i][j][k]*v[i]*v[j]*v[k];
                    for(int l=0;l<NHiggs;l++)
                    {
                        res += 1.0/24.0*Curvature_Higgs_CT_L4[i][j][k][l]*v[i]*v[j]*v[k]*v[l];
                    }
                }
            }
        }
    }
    else if(diff <= NHiggs)
    {
        int i = diff -1;
        res = Curvature_Higgs_CT_L1[i];
        for(int j=0;j<NHiggs;j++)
        {
            res += Curvature_Higgs_CT_L2[i][j]*v[j];
            for(int k=0;k<NHiggs;k++)
            {
                res += 0.5*Curvature_Higgs_CT_L3[i][j][k]*v[j]*v[k];
                for(int l=0;l<NHiggs;l++)
                {
                    res += 1.0/6.0*Curvature_Higgs_CT_L4[i][j][k][l]*v[j]*v[k]*v[l];
                }
            }
        }
    }

    return res;
}


/**
 * Calculates the effective potential and its derivatives.
 */
double Class_Potential_Origin::VEff(const std::vector<double>& v,
		double Temp=0, int diff=0) {
		double resOut = 0;
		resOut = VTree(v,diff);
		resOut+= CounterTerm(v,diff);
		resOut+= V1Loop(v,Temp,diff);
		return resOut;
}


/**
 * Calculates the Coleman-Weinberg and temperature-dependent 1-loop part of the effective potential and its derivatives.
 */
double Class_Potential_Origin::V1Loop(const std::vector<double>& v, double Temp, int diff)
{
	double res = 0;
	double HiggsMasses[NHiggs];

    std::vector<double> HiggsMassesVec,QuarkMassesVec,GaugeMassesVec,LeptonMassesVec,HiggsMassesZeroTempVec,GaugeMassesZeroTempVec;
    HiggsMassesSquared(HiggsMassesVec,v,Temp);
    GaugeMassesSquared(GaugeMassesVec,v,Temp);
    GaugeMassesSquared(GaugeMassesZeroTempVec,v,0);
    QuarkMassesSquared(QuarkMassesVec,v);
    LeptonMassesSquared(LeptonMassesVec,v);


	if(diff == 0)
	{
        if(C_UseParwani)
        {
            for(int k=0;k<NHiggs;k++) res += boson(HiggsMassesVec[k],Temp,C_CWcbHiggs,0);
            for(int k=0;k<NGauge;k++) res += boson(GaugeMassesVec[k],Temp,C_CWcbGB,0);
            for(int k=0;k<NGauge;k++) res += 2*boson(GaugeMassesZeroTempVec[k],Temp,C_CWcbGB,0);
            for(int k=0;k<NQuarks;k++) res += -6*fermion(QuarkMassesVec[k],Temp,0);
            for(int k=0;k<NLepton;k++) res += -2*fermion(LeptonMassesVec[k],Temp,0);
        }
        else{
            HiggsMassesSquared(HiggsMassesZeroTempVec,v,0);
            for(int k=0;k<NHiggs;k++) res += boson(HiggsMassesZeroTempVec[k],Temp,C_CWcbHiggs,0);
            for(int k=0;k<NGauge;k++) res += 3*boson(GaugeMassesZeroTempVec[k],Temp,C_CWcbGB,0);
            double AddContQuark=0;
            for(int k=0;k<NQuarks;k++) AddContQuark += -2*fermion(QuarkMassesVec[k],Temp,0);
            for(int k=0;k<NColour;k++) res += AddContQuark;
            for(int k=0;k<NLepton;k++) res += -2*fermion(LeptonMassesVec[k],Temp,0);


            double VDebay = 0;
            for(int k=0;k<NHiggs;k++)
            {
                if( HiggsMassesVec[k] > 0 ) VDebay += std::pow(HiggsMassesVec[k],1.5);
                if( HiggsMassesZeroTempVec[k] > 0) VDebay += -std::pow(HiggsMassesZeroTempVec[k],1.5);
            }
            for(int k=0;k<NGauge;k++)
            {
                if(GaugeMassesVec[k] > 0 ) VDebay += std::pow(GaugeMassesVec[k],1.5);
                if(GaugeMassesZeroTempVec[k] > 0 ) VDebay += -std::pow(GaugeMassesZeroTempVec[k],1.5);
            }

            VDebay *= -Temp/(12*M_PI);
            res += VDebay;
        }
	}


        return res;




}

/**
 * Calculates the Debye corrections to the Higgs mass matrix.
 * If you can provide CalculateDebyeSimplified() with the Matrix as this will reduce the runtime.
 */
void Class_Potential_Origin::CalculateDebye()
{
    bool Debug = false;



    if(!SetCurvatureDone) SetCurvatureArrays();
    if(Debug) std::cout << "Debug turned on in " << __func__ << std::endl;

    bool Done = CalculateDebyeSimplified();
    if(Debug and Done) std::cout << "Done = true " << std::endl;
    if(Debug and !Done) std::cout << "Done = false " << std::endl;
    if(!Done)
      {
	for(int i=0;i<NHiggs;i++)
	    {
	        for(int j=i;j<NHiggs;j++)
	        {
	            DebyeHiggs[i][j] = 0;
	            for(int k=0;k<NHiggs;k++)
	            {
	                DebyeHiggs[i][j] += 0.5*Curvature_Higgs_L4[i][j][k][k]/12.0;
	            }
	            for(int k=0;k<NGauge;k++)
	            {
	                DebyeHiggs[i][j] += 3*0.5*Curvature_Gauge_G2H2[k][k][i][j]/12.0;
	            }

	            for(int a=0;a<NQuarks;a++)
	            {
	                for(int b=0;b<NQuarks;b++)
	                {
	                    double tmp = 0.5*(std::conj(Curvature_Quark_F2H1[a][b][j])*Curvature_Quark_F2H1[a][b][i] + std::conj(Curvature_Quark_F2H1[a][b][i])*Curvature_Quark_F2H1[a][b][j]).real();
	                    DebyeHiggs[i][j] += 6.0/24.0*tmp;
	                }
	            }

	            for(int a=0;a<NLepton;a++)
	            {
	                for(int b=0;b<NLepton;b++)
	                {
	                    double tmp = 0.5*(std::conj(Curvature_Lepton_F2H1[a][b][j])*Curvature_Lepton_F2H1[a][b][i] + std::conj(Curvature_Lepton_F2H1[a][b][i])*Curvature_Lepton_F2H1[a][b][j]).real();
	                    DebyeHiggs[i][j] += 2.0/24.0*tmp;
	                }
	            }

//	            if(i==j) DebyeHiggs[i][j] *= 0.5;
	        }
	    }

		for(int i=0;i<NHiggs;i++){
			for(int j=i;j<NHiggs;j++) {
				if(std::abs(DebyeHiggs[i][j])<=1e-5) DebyeHiggs[i][j] = 0;
			}
		}


	    for(int i=0;i<NHiggs;i++)
	    {
	        for(int j=0;j<i;j++) {
	        	DebyeHiggs[i][j] = DebyeHiggs[j][i];
	        }
	    }
      }


    int nHiggsGauge = 0;
    for(int i=0;i<NHiggs;i++)
    {
        if(Curvature_Gauge_G2H2[0][0][i][i] != 0) nHiggsGauge++;
    }

    if(Debug)
    {
        std::cout << "Debye : " << std::endl;
        for(int i=0;i<NHiggs;i++)
        {
            for(int j=0;j<NHiggs;j++) std::cout << DebyeHiggs[i][j] << "\t";
            std::cout << std::endl;
        }
        std::cout << "NHG = " << nHiggsGauge << std::endl;
    }


}

/**
 * Calculates the Debye corrections to the gauge sector. By using CalculateDebyeGaugeSimplified() the runtime can be reduced.
 */
void Class_Potential_Origin::CalculateDebyeGauge(){
  bool Debug = false;
  if(Debug) std::cout << "Debug turned on in " << __func__ << std::endl;
    for(int i=0;i<NGauge;i++){
        for(int j=0;j<NGauge;j++) DebyeGauge[i][j] = 0;
    }

    bool Done = CalculateDebyeGaugeSimplified();
    if(Done) return ;

    int nGaugeHiggs = 0;



    for(int i=0;i<NHiggs;i++) {
	if(Curvature_Gauge_G2H2[0][0][i][i] != 0){
	    nGaugeHiggs++;
	}
    }
    for(int i=0;i<NGauge;i++)
      {
    	double GaugeFac = 0;
    	for(int k=0;k<NHiggs;k++)
    	{
    		GaugeFac += Curvature_Gauge_G2H2[i][i][k][k];
    	}
    	GaugeFac*= 1.0/nGaugeHiggs;
    	if(Debug)std::cout << "gaugefac = " << GaugeFac << std::endl;
    	DebyeGauge[i][i] = 2.0/3.0*(nGaugeHiggs/8.0 + 5)*GaugeFac;
      }

    for(int i=0;i<NGauge;i++){
    	for(int j=0;j<NGauge;j++)
    	{
    		if(std::abs(DebyeGauge[i][j]) <=1e-5) DebyeGauge[i][j] = 0;
    	}
    }

    if(Debug)
      {
	for(int i=0;i<NGauge;i++)
	  {
	    for(int j=0;j<NGauge;j++) std::cout << DebyeGauge[i][j] << "\t";
	    std::cout << std::endl;
	  }
	std::cout << std::endl;
      }
}

/**
 * Initializes all vectors needed for the calculations.
 */
void Class_Potential_Origin::initVectors(){
bool Debug = false;

	VEVSymmetric.resize(NHiggs);

    Curvature_Higgs_L1.resize(NHiggs);
  Curvature_Higgs_L2.resize(NHiggs);
  Curvature_Higgs_L3.resize(NHiggs);
  Curvature_Higgs_L4.resize(NHiggs);

  Curvature_Higgs_CT_L1.resize(NHiggs);
  Curvature_Higgs_CT_L2.resize(NHiggs);
  Curvature_Higgs_CT_L3.resize(NHiggs);
  Curvature_Higgs_CT_L4.resize(NHiggs);

  LambdaHiggs_3.resize(NHiggs);
  LambdaHiggs_3_CT.resize(NHiggs);


  DebyeHiggs.resize(NHiggs);
  for(int i=0;i<NHiggs;i++)
    {
	  VEVSymmetric[i] = 0;

      DebyeHiggs[i].resize(NHiggs);
      Curvature_Higgs_L1[i] = 0;
      Curvature_Higgs_L2[i].resize(NHiggs);
      Curvature_Higgs_L3[i].resize(NHiggs);
      Curvature_Higgs_L4[i].resize(NHiggs);

      Curvature_Higgs_CT_L1[i] = 0;
      Curvature_Higgs_CT_L2[i].resize(NHiggs);
      Curvature_Higgs_CT_L3[i].resize(NHiggs);
      Curvature_Higgs_CT_L4[i].resize(NHiggs);

      LambdaHiggs_3[i].resize(NHiggs);
      LambdaHiggs_3_CT[i].resize(NHiggs);

      for(int j=0;j<NHiggs;j++)
	{
        DebyeHiggs[i][j] = 0;

        Curvature_Higgs_CT_L2[i][j] = 0;
        Curvature_Higgs_CT_L3[i][j].resize(NHiggs);
        Curvature_Higgs_CT_L4[i][j].resize(NHiggs);

        Curvature_Higgs_L2[i][j] = 0;
        Curvature_Higgs_L3[i][j].resize(NHiggs);
	Curvature_Higgs_L4[i][j].resize(NHiggs);

	LambdaHiggs_3[i][j].resize(NHiggs);
	LambdaHiggs_3_CT[i][j].resize(NHiggs);


	  for(int k=0;k<NHiggs;k++)
	    {
	      Curvature_Higgs_L3[i][j][k] = 0;
	      Curvature_Higgs_L4[i][j][k].resize(NHiggs);

	      Curvature_Higgs_CT_L3[i][j][k] = 0;
	      Curvature_Higgs_CT_L4[i][j][k].resize(NHiggs);

	      for(int l=0;l<NHiggs;l++) {
		  Curvature_Higgs_L4[i][j][k][l] = 0;
		  Curvature_Higgs_CT_L4[i][j][k][l] = 0;
	      }

	    }
	}
    }
    Curvature_Gauge_G2H2.resize(NGauge);
    DebyeGauge.resize(NGauge);
    LambdaGauge_3.resize(NGauge);
    for(int a=0;a<NGauge;a++)
    {
        DebyeGauge[a].resize(NGauge);
        Curvature_Gauge_G2H2[a].resize(NGauge);
        LambdaGauge_3[a].resize(NGauge);
        for(int b=0;b<NGauge;b++)
        {
            DebyeGauge[a][b] = 0;
            Curvature_Gauge_G2H2[a][b].resize(NHiggs);
            LambdaGauge_3[a][b].resize(NHiggs);
            for(int i=0;i<NHiggs;i++)
            {
                Curvature_Gauge_G2H2[a][b][i].resize(NHiggs);
                for(int j=0;j<NHiggs;j++) Curvature_Gauge_G2H2[a][b][i][j] = 0;
            }
        }
    }
    Curvature_Lepton_F2H1.resize(NLepton);
    LambdaLepton_3.resize(NLepton);
    LambdaLepton_4.resize(NLepton);
    for(int i=0;i<NLepton;i++)
    {
        Curvature_Lepton_F2H1[i].resize(NLepton);
        LambdaLepton_3[i].resize(NLepton);
	LambdaLepton_4[i].resize(NLepton);
        for(int j=0;j<NLepton;j++)
        {
            Curvature_Lepton_F2H1[i][j].resize(NHiggs);
            LambdaLepton_3[i][j].resize(NHiggs);
	    LambdaLepton_4[i][j].resize(NHiggs);
            for(int l=0;l<NHiggs;l++) {
        	Curvature_Lepton_F2H1[i][j][l] = 0;
        	LambdaLepton_4[i][j][l].resize(NHiggs);
            }
        }
    }
    Curvature_Quark_F2H1.resize(NQuarks);
    LambdaQuark_3.resize(NQuarks);
    LambdaQuark_4.resize(NQuarks);
    for(int i=0;i<NQuarks;i++)
    {
        Curvature_Quark_F2H1[i].resize(NQuarks);
        LambdaQuark_3[i].resize(NQuarks);
	LambdaQuark_4[i].resize(NQuarks);

        for(int j=0;j<NQuarks;j++) {
            Curvature_Quark_F2H1[i][j].resize(NHiggs);
            LambdaQuark_3[i][j].resize(NHiggs);
	    LambdaQuark_4[i][j].resize(NHiggs);
            for(int l=0;l<NHiggs;l++) {
        	Curvature_Quark_F2H1[i][j][l] = 0;
        	LambdaQuark_4[i][j][l].resize(NHiggs);
            }
        }

    }

    HiggsVev.resize(NHiggs);
    for(int i=0;i<NHiggs;i++) HiggsVev[i] = 0;
  if(Debug) std::cout << "Resize done " << std::endl;
}


/**
 * Resets all bools. Needed if you want to deal with multiple points one after another with the same pointer.
 */
void Class_Potential_Origin::resetbools()
{
  SetCurvatureDone = false;
  CalcCouplingsdone=false;
  CalculatedTripleCopulings = false;
}


bool Class_Potential_Origin::CheckNLOVEV(const std::vector<double>& v)
{
	//std::vector<double> vPotential;
	double MaxDiff = 0;
	double AllowedDifference = 1;
	for(int i=0;i<nVEV;i++)
	{
		double tmp=std::abs(std::abs(v[i])-std::abs(vevTreeMin[i]));
		if(tmp > MaxDiff) MaxDiff = tmp;
	}

	return (MaxDiff < AllowedDifference);
}



double Class_Potential_Origin::EWSBVEV(std::vector<double> v)
{
  double res=0;
//  for(int i=0;i<NHiggs-1;i++) res += std::pow(v.at(i),2);
//  std::cout << "res = " << res << std::endl;
//  res=0;
  for(int i=0;i<NHiggs;i++)
  {
	  double checkgauge=0;
	  for(int j=0;j<NGauge;j++)
	  {
		  checkgauge += std::abs(Curvature_Gauge_G2H2[j][j][i][i]);
	  }
	  if(checkgauge != 0) res += std::pow(v.at(i),2);
  }
//  std::cout << "Res v2 =  " << res << std::endl;
  res = std::sqrt(res);

  if(res <= 0.5)
    {
//	  std::cout << nVEV << std::endl;
      ModifiedVEVVectorDim.resize(nVEV);
      for(int i=0;i<nVEV;i++)
      {
    	  double checkgauge=0;
    	  int FieldNumber = VevOrder[i];
//    	  std::cout << "i = " << i << "\t" << FieldNumber << std::endl;
		  for(int j=0;j<NGauge;j++)
		  {
			  checkgauge += std::abs(Curvature_Gauge_G2H2[j][j][FieldNumber][FieldNumber]);
		  }
		  if(checkgauge != 0) {
			  ModifiedVEVVectorDim[i] = 0;
		  }
		  else{
			  ModifiedVEVVectorDim[i] = v.at(FieldNumber);
		  }
      }
      res = 0;
    }
  return res;
}
