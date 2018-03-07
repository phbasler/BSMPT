/*
 * ClassPotentialC2HDM.h
 *
 *  Created on: Mar 10, 2017
 *      Author: basler
 */

#ifndef SRC_CLASSPOTENTIALC2HDM_H_
#define SRC_CLASSPOTENTIALC2HDM_H_

#include "ClassPotentialOrigin.h"

class Class_Potential_C2HDM : public Class_Potential_Origin
{
public:
  Class_Potential_C2HDM ();
  virtual  ~Class_Potential_C2HDM ();

  long double L1=0,L2=0,L3=0,L4=0,RL5=0,RealMMix=0,u1=0,u2=0;
  long double IL5=0,Iu3=0;
  long double DL1CT=0,DL2CT=0,DL3CT=0,DL4CT=0,DRL5CT=0,Du2CT=0,Du1CT=0,DRu3CT=0;
  long double DIL5CT=0,DIu3CT=0;
  double DT1=0,DT2=0,DT3=0,DTCharged=0;
  double TanBeta=0,C_CosBeta=0,C_SinBeta=0,C_CosBetaSquared=0,C_SinBetaSquared=0;
  double beta=0;
  long double M1=0,M2=0,M3=0,alpha1=0,alpha2=0,alpha3=0,MHC=0;
  long double MSM=0,MhUp=0,MhDown=0;
  int Type=0;
  double CTempC1=0,CTempC2=0,CTempCS=0;
  double R_Hh_1=0,R_Hh_2=0,R_Hh_3=0,R_Hl_1=0,R_Hl_2=0,R_Hl_3=0,R_Hsm_1=0,R_Hsm_2=0,R_Hsm_3=0;






  void ReadAndSet(const std::string& linestr, std::vector<double>& par);
  std::string addLegendCT();
  std::string addLegendTemp();
  std::string addLegendTripleCouplings();
  std::string addLegendVEV();

  void set_gen(const std::vector<double>& par);
  void set_CT_Pot_Par(const std::vector<double>& par);
  void write();

  void TripleHiggsCouplings();
  void calc_CT(std::vector<double>& par);


  void MinimizeOrderVEV(const std::vector<double>& vevminimizer, std::vector<double>& vevFunction);

  void SetCurvatureArrays();
  bool CalculateDebyeSimplified();
  bool CalculateDebyeGaugeSimplified();
  double VTreeSimplified(const std::vector<double>& v);
  double VCounterSimplified(const std::vector<double>& v);
  void Debugging(const std::vector<double>& input, std::vector<double>& output);

  bool IncludeChargeBreakingVEV = true;
  bool CTAlternative=false;
  double DiffDelta=0;

  int PosSM=-1;


};

#endif /* SRC_CLASSPOTENTIALC2HDM_H_ */
