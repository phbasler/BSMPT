/*
 * ClassPotentialCN2HDM.h
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

/**
  * @file
  * Model file for the CP-violating 2HDM + complex singlet
 */

#ifndef SRC_CLASSCN2HDM_H_
#define SRC_CLASSCN2HDM_H_



#include <BSMPT/models/ClassPotentialOrigin.h>

namespace BSMPT {
namespace Models {


/**
 * @brief The Class_Potential_CN2HDM class
 * Another take at the CP-violating 2HDM + complex singlet
 *
 * Scalar lagrangian:
 * \f$ -L_S = m112 \Phi_1^\dagger \Phi_1  + m222 \Phi_2^\dagger \Phi_2 + \frac{L1}{2} (\Phi_1^\dagger \Phi_1)^2 + \frac{L2}{2} (\Phi_2^\dagger \Phi_2)^\dagger + L3 \Phi_1^\dagger \Phi_1 \Phi_2^\dagger \Phi_2 + L4 \Phi_1^\dagger \Phi_2 \Phi_2^\dagger \Phi_1 + \frac{ReL5 + i ImL5}{2} (\Phi_1^\dagger \Phi_2)^2 + \frac{ReL5 - i ImL5}{2} (\Phi_2^\dagger \Phi_1)^2 + \frac{L6}{8} + |S|^4 + \frac{L7}{2} \Phi_1^\dagger \Phi_1 |S|^2 + \frac{L8}{2} \Phi_2^\dagger \Phi_2 |S|^2 + \frac{1}{2} ms2 |S|^2 - (Rem122 + i Imm122) \Phi_1^\dagger \Phi_2 - (Rem122 - i Imm122) \Phi_2^\dagger \Phi_1 - \frac{mDM2}{4} (S^2 + (S^*)^2) \f$
 *
 * Yukawa lagrangian:
 *
 * \f$ -L_Y = \overline{U}_L V \text{diag}(y_d,y_s,y_b) D_R \Phi_d^+  + \overline{D}_L \text{diag}(y_d,y_s,y_b) D_R \Phi_d^0 + \overline{U}_L \text{diag}(y_u, y_c, y_t) U_R (\Phi_2^0)^\ast - \overline{D}_L V^\dagger \text{diag}(y_u, y_c, y_t) U_R \Phi_2^- + \overline{E}_L \text{diag}(y_e, y_\mu, y_\tau) E_R \Phi_l^0 +\overline{\nu}_L \text{diag}(y_e, y_\mu, y_\tau) E_R \Phi_l^+ + c.c.  \f$
 *
 * Gauge Lagrangian:
 *
 * \f$ L_G = (D_\mu \Phi_1)^\dagger (D^\mu \Phi_1) + (D_\mu \Phi_2)^\dagger (D^\mu \Phi_2)\f$
 *
 * with
 * * \f$ \Phi_1 = \begin{pmatrix} \Phi_1^+ \\  \Phi_1^0 \end{pmatrix} = 1/\sqrt{2} \begin{pmatrix} higgsbase[0] + I*higgsbase[1] \\  higgsbase[2] + I *higgsbase[3] \end{pmatrix}\,, \f$
 * * \f$ \Phi_2 = \begin{pmatrix} \Phi_2^+ \\  \Phi_2^0 \end{pmatrix} = 1/\sqrt{2} \begin{pmatrix} higgsbase[4] + I*higgsbase[5] \\  higgsbase[6] + I *higgsbase[7] \end{pmatrix} \f$
 * * \f$ S = 1/\sqrt{2}( higgsbase[8] + I* higgsbase[9]) \f$
 * * \f$ D_\mu = -I C\_{}g/2 * W_\mu^a \sigma^a -I C\_{}gs/2 B_\mu = -\frac{I}{2} \begin{pmatrix} C\_{}gs B + C\_{}g W3  & C\_{}g (W1 -I W2) \\ C\_{}g (W1 +I W2) & C\_{}gs B - C\_{}g W3 \end{pmatrix} \f$
 *
 * with the gauge base = [W1,W2,W3,B]
 * * \f$ \overline{U}_L = \begin{pmatrix} \overline{u}_L &  \overline{c}_L & \overline{t}_L \end{pmatrix}\,, \f$
 * * \f$ U_R = \begin{pmatrix} u_R \\ c_R \\ t_R \end{pmatrix}\,, \f$
 * * \f$ \overline{D}_L = \begin{pmatrix} \overline{d}_L & \overline{s}_L & \overline{b}_L \end{pmatrix}\,, \f$
 * * \f$ D_R = \begin{pmatrix} d_R \\ s_R \\ b_R \end{pmatrix}\,, \f$
 * * \f$ \overline{E}_L = \begin{pmatrix} \overline{e}_L & \overline{\mu}_L & \overline{\tau}_L \end{pmatrix}\,, \f$
 * * \f$ E_R = \begin{pmatrix} e_R \\ \mu_R \\ \tau_R \end{pmatrix}\,, \f$
 * * \f$ \overline{\nu}_L = \begin{pmatrix} \overline{\nu}_{e,L} & \overline{\nu}_{\mu,L} & \overline{\nu}_{\tau,L} \end{pmatrix} \f$
 *
 * The basis of the quarks given by \f$ [u_R, c_R, t_R, d_R, s_R, b_R, \overline{u}_L, \overline{c}_L, \overline{t}_L, \overline{d}_L, \overline{s}_L, \overline{b}_L ] \f$
 * and for the leptons \f$ [\overline{e}_L, e_R, \overline{\mu}_L , \mu_R, \overline{\tau}_L, \tau_R, \overline{\nu}_{e,L}, \overline{\nu}_{\mu,L}, \overline{\nu}_{\tau,L} ] \f$
 * The Yukawa types defined by
 *
 * * Type I : \f$ \Phi_d = \Phi_l = \Phi_2 \f$
 * * Type II: \f$ \Phi_d = \Phi_l = \Phi_1 \f$
 * * Type  LS(=3): \f$ \Phi_d = \Phi_2\,, \Phi_l = \Phi_1 \f$
 * * Type FL(=4): \f$ \Phi_d = \Phi_1 \,,\Phi_l = \Phi_2 \f$
 */
class Class_Potential_CN2HDM : public Class_Potential_Origin
{
public:
	Class_Potential_CN2HDM ();
  virtual
  ~Class_Potential_CN2HDM ();



  double Rem122=0,Imm122=0,L1=0,L2=0,L3=0,L4=0,ReL5=0,ImL5=0,L6=0,L7=0,L8=0,ms2=0,m112=0,m222=0,mds=0,vs=0,TanBeta=0, mDM2=0;
  int Type = 0;

  std::vector<std::string>
  paramsName{"Rem122","Imm122","m112","m222","ms2","mds","L1","L2","L3","L4","ReL5","ImL5","L6","L7","L8"};

  double v1=0,v2=0;

  double C_CosBetaSquared = 0, C_CosBeta = 0, C_SinBetaSquared = 0, C_SinBeta = 0;

  double
  dRem122=0,dImm122=0,dL1=0,dL2=0,dL3=0,dL4=0,dReL5=0,dImL5=0,dL6=0,dL7=0,dL8=0,dms2=0,dm112=0,dm222=0,dmDM2=0;
  double dT1=0,dT2=0,dT3=0,dT4=0,dT5=0,dT6=0,dT7=0,dT8=0,dT9=0,dT10=0;

  const bool UseDMVEV = true; //Allow the DM candidate to develop a VEV





  void ReadAndSet(const std::string& linestr, std::vector<double>& par) override;
  std::vector<std::string> addLegendCT() const override;
  std::vector<std::string> addLegendTemp() const override;
  std::vector<std::string> addLegendTripleCouplings() const override;
  std::vector<std::string> addLegendVEV() const override;


  void set_gen(const std::vector<double>& par) override;
  void set_CT_Pot_Par(const std::vector<double>& par) override;
  void write() const override;

  void TripleHiggsCouplings() override;
  std::vector<double> calc_CT() const override;


  void SetCurvatureArrays() override;
  bool CalculateDebyeSimplified() override;
  bool CalculateDebyeGaugeSimplified() override;
  double VTreeSimplified(const std::vector<double>& v) const override;
  double VCounterSimplified(const std::vector<double>& v) const override;
  void Debugging(const std::vector<double>& input, std::vector<double>& output) const override;
};

}
}
#endif /* SRC_CLASSCN2HDM_H_ */
