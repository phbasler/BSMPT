/*
 * ClassPotentialR2HDM.h
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

#ifndef SRC_CLASSPOTENTIALR2HDM_H_
#define SRC_CLASSPOTENTIALR2HDM_H_

#include "ClassPotentialOrigin.h"
using namespace Eigen;

class Class_Potential_R2HDM : public Class_Potential_Origin
{
public:
  Class_Potential_R2HDM ();
  virtual  ~Class_Potential_R2HDM ();

  long double L1=0,L2=0,L3=0,L4=0,RL5=0,RealMMix=0,u1=0,u2=0;
  long double DL1CT=0,DL2CT=0,DL3CT=0,DL4CT=0,DRL5CT=0,Du2CT=0,Du1CT=0,DRu3CT=0;
  double DT1=0,DT2=0,DT3=0,DTCharged=0;
  double TanBeta=0,C_CosBeta=0,C_SinBeta=0,C_CosBetaSquared=0,C_SinBetaSquared=0;
  double beta=0;
  long double Mh=0,MH=0,MA=0,MHP=0 , alpha = 0;
  int Type=0;
  double CTempC1=0,CTempC2=0,CTempCS=0;




  void ReadAndSet(const std::string& linestr, std::vector<double>& par);
  std::string addLegendCT();
  std::string addLegendTemp();
  std::string addLegendTripleCouplings();
  std::string addLegendVEV();

  void set_gen(const std::vector<double>& par);
  void set_CT_Pot_Par(const std::vector<double>& par);
  void write();

//  double WeinbergCurvature(std::vector<double>& res,bool CalcTriCouplings=false);
  void TripleHiggsCouplings();
  void calc_CT(std::vector<double>& par);

    void MinimizeOrderVEV(const std::vector<double>& vevminimizer, std::vector<double>& vevFunction);
  void SetCurvatureArrays();

  bool CalculateDebyeSimplified();
  bool CalculateDebyeGaugeSimplified();
  double VTreeSimplified(const std::vector<double>& v);
  double VCounterSimplified(const std::vector<double>& v);
  void Debugging(const std::vector<double>& input, std::vector<double>& output);
};

#endif /* SRC_CLASSPOTENTIALR2HDM_H_ */
