/*
 * ClassPotentialCxSM.h
 *
 *  Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller

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
  * Model file for the SM + complex singlet with the tree-level potential
  */

#pragma once

#include <BSMPT/models/ClassPotentialOrigin.h>

namespace BSMPT {
namespace Models {

/**
 * @brief The Class_CxSM class
 * Implementation of the CxSM as shown in the manual of BSMPT v2
 *
 * * \f$ -L_S = msq/2 H^\dagger H + (\lambda/4) (H^\dagger H)^2 + (\delta2 / 2) (H^\dagger H) | S | ^2 + ( b2 / 2) | S |^2
 *         + ( d2 / 4) |S|^4 + ( ( b1 / 4) S^2 + a1 S + c.c.) \f$
 * * \f$ -L_Y = \overline{U}_L V \text{diag}(y_d,y_s,y_b) D_R H^+  + \overline{D}_L \text{diag}(y_d,y_s,y_b) D_R H^0 + \overline{U}_L \text{diag}(y_u, y_c, y_t) U_R (H^0)^\ast - \overline{D}_L V^\dagger \text{diag}(y_u, y_c, y_t) U_R H^- + \overline{E}_L \text{diag}(y_e, y_\mu, y_\tau) E_R H^0 +\overline{\nu}_L \text{diag}(y_e, y_\mu, y_\tau) E_R H^+ + c.c.  \f$
 * * \f$ L_G = (D_\mu H)^\dagger (D^\mu H) \f$
 *
 * with
 *
 * * \f$ H = \begin{pmatrix} H^+ \\  H^0 \end{pmatrix} = 1/\sqrt{2} \begin{pmatrix} higgsbase[0] + I*higgsbase[1] \\  higgsbase[3] + I *higgsbase[2] \end{pmatrix}\,, \f$
 * * \f$ S = 1/\sqrt{2}(higgsfield[4] + I higgsfield[5]) \f$
 * * \f$ D_\mu = -I C\_{}g/2 * W_\mu^a \sigma^a -I C\_{}gs/2 B_\mu = -\frac{I}{2} \begin{pmatrix} C\_{}gs B + C\_{}g W3  & C\_{}g (W1 -I W2) \\ C\_{}g (W1 +I W2) & C\_{}gs B - C\_{}g W3 \end{pmatrix} \f$
 * with the gauge base = [W1,W2,W3,B]
 *
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
 *
 */
class Class_CxSM : public Class_Potential_Origin
{
public:
  Class_CxSM ();
  virtual
  ~Class_CxSM ();


  // Add here your parameters for the Lagrangian as well as for the counterterm potential
  // Add here your variables in which you will save the Debye correction factors

  double msq, lambda, delta2, b2, d2, Reb1 , Imb1,  Rea1, Ima1;

  double dmsq, dlambda, ddelta2, db2, dd2, dReb1, dImb1, dRea1,  dIma1, dT1, dT2, dT3, dT4, dT5, dT6;

  double g1 = C_gs;
  double g2 = C_g;

  double vh, vs, va;




  void ReadAndSet(const std::string& linestr, std::vector<double>& par) override;
  std::vector<std::string> addLegendCT() const override;
  std::vector<std::string> addLegendTemp() const override;
  std::vector<std::string> addLegendTripleCouplings() const override;
  std::vector<std::string> addLegendVEV() const override;


  /**
   * @brief set_gen
   * @param par[0] = vh
   * @param par[1] = vs
   * @param par[2] = va
   * @param par[3] = msq
   * @param par[4] = lambda
   * @param par[5] = delta2
   * @param par[6] = b2
   * @param par[7] = d2
   * @param par[8] = Reb1
   * @param par[9] = Imb1
   * @param par[10] = Rea1
   * @param par[11] = Ima1
   */
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
