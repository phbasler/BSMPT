// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#ifndef SRC_CLASSPOTENTIALC2HDM_H_
#define SRC_CLASSPOTENTIALC2HDM_H_

#include <BSMPT/models/ClassPotentialOrigin.h>
namespace BSMPT
{
namespace Models
{

/**
 * @brief The Class_Potential_C2HDM class
 * Implementation of the CP-violating 2HDM as given in the manual
 *
 *
 *
 * * \f$ -L_S = u1 \Phi_1^\dagger \Phi_1 + u2 \Phi_2^\dagger \Phi_2 +
 * \frac{L1}{2} (\Phi_1^\dagger \Phi_1)^2 + \frac{L2}{2} (\Phi_2^\dagger
 * \Phi_2)^2 + L3 \Phi_1^\dagger \Phi_1 \Phi_2^\dagger \Phi_2 + L4
 * \Phi_1^\dagger \Phi_2 \Phi_2^\dagger \Phi_1 + [ \frac{RL5+I*IL5}{2}
 * (\Phi_1^\dagger \Phi_2)^2 - (RealMMix + I* Iu3) \Phi_1^\dagger \Phi_2 ] \f$
 * * \f$ -L_Y = \overline{U}_L V \text{diag}(y_d,y_s,y_b) D_R \Phi_d^+  +
 * \overline{D}_L \text{diag}(y_d,y_s,y_b) D_R \Phi_d^0 + \overline{U}_L
 * \text{diag}(y_u, y_c, y_t) U_R (\Phi_2^0)^\ast - \overline{D}_L V^\dagger
 * \text{diag}(y_u, y_c, y_t) U_R \Phi_2^- + \overline{E}_L \text{diag}(y_e,
 * y_\mu, y_\tau) E_R \Phi_l^0 +\overline{\nu}_L \text{diag}(y_e, y_\mu, y_\tau)
 * E_R \Phi_l^+ + c.c.  \f$
 * * \f$ L_G = (D_\mu \Phi_1)^\dagger (D^\mu \Phi_1) + (D_\mu \Phi_2)^\dagger
 * (D^\mu \Phi_2)\f$
 *
 * with
 *
 * * \f$ \Phi_1 = \begin{pmatrix} \Phi_1^+ \\  \Phi_1^0 \end{pmatrix} =
 * 1/\sqrt{2} \begin{pmatrix} higgsbase[0] + I*higgsbase[1] \\  higgsbase[4] + I
 * *higgsbase[5] \end{pmatrix}\,, \f$
 * * \f$ \Phi_2 = \begin{pmatrix} \Phi_2^+ \\  \Phi_2^0 \end{pmatrix} =
 * 1/\sqrt{2} \begin{pmatrix} higgsbase[2] + I*higgsbase[3] \\  higgsbase[6] + I
 * *higgsbase[7] \end{pmatrix}  \f$
 * * \f$ D_\mu = -I C\_{}g/2 * W_\mu^a \sigma^a -I C\_{}gs/2 B_\mu =
 * -\frac{I}{2} \begin{pmatrix} C\_{}gs B + C\_{}g W3  & C\_{}g (W1 -I W2) \\
 * C\_{}g (W1 +I W2) & C\_{}gs B - C\_{}g W3 \end{pmatrix} \f$ with the gauge
 * base = [W1,W2,W3,B]
 *
 * * \f$ \overline{U}_L = \begin{pmatrix} \overline{u}_L &  \overline{c}_L &
 * \overline{t}_L \end{pmatrix}\,, \f$
 * * \f$ U_R = \begin{pmatrix} u_R \\ c_R \\ t_R \end{pmatrix}\,, \f$
 * * \f$ \overline{D}_L = \begin{pmatrix} \overline{d}_L & \overline{s}_L &
 * \overline{b}_L \end{pmatrix}\,, \f$
 * * \f$ D_R = \begin{pmatrix} d_R \\ s_R \\ b_R \end{pmatrix}\,, \f$
 * * \f$ \overline{E}_L = \begin{pmatrix} \overline{e}_L & \overline{\mu}_L &
 * \overline{\tau}_L \end{pmatrix}\,, \f$
 * * \f$ E_R = \begin{pmatrix} e_R \\ \mu_R \\ \tau_R \end{pmatrix}\,, \f$
 * * \f$ \overline{\nu}_L = \begin{pmatrix} \overline{\nu}_{e,L} &
 * \overline{\nu}_{\mu,L} & \overline{\nu}_{\tau,L} \end{pmatrix} \f$
 *
 * The basis of the quarks given by \f$ [u_R, c_R, t_R, d_R, s_R, b_R,
 * \overline{u}_L, \overline{c}_L, \overline{t}_L, \overline{d}_L,
 * \overline{s}_L, \overline{b}_L ] \f$ and for the leptons \f$ [\overline{e}_L,
 * e_R, \overline{\mu}_L , \mu_R, \overline{\tau}_L, \tau_R,
 * \overline{\nu}_{e,L}, \overline{\nu}_{\mu,L}, \overline{\nu}_{\tau,L} ] \f$
 *
 * The Yukawa types defined by
 *
 * * Type I : \f$ \Phi_d = \Phi_l = \Phi_2 \f$
 * * Type II: \f$ \Phi_d = \Phi_l = \Phi_1 \f$
 * * Type  LS(=3): \f$ \Phi_d = \Phi_2\,, \Phi_l = \Phi_1 \f$
 * * Type FL(=4): \f$ \Phi_d = \Phi_1 \,,\Phi_l = \Phi_2 \f$
 */
class Class_Potential_C2HDM : public Class_Potential_Origin
{
public:
  Class_Potential_C2HDM(const ISMConstants &smConstants);
  virtual ~Class_Potential_C2HDM() override;

  bool UseHsmNotationInTripleHiggs = false;

  double L1 = 0, L2 = 0, L3 = 0, L4 = 0, RL5 = 0, RealMMix = 0, u1 = 0, u2 = 0;
  double IL5 = 0, Iu3 = 0;
  double DL1CT = 0, DL2CT = 0, DL3CT = 0, DL4CT = 0, DRL5CT = 0, Du2CT = 0,
         Du1CT = 0, DRu3CT = 0;
  double DIL5CT = 0, DIu3CT = 0;
  double DT1 = 0, DT2 = 0, DT3 = 0, DTCharged = 0;
  double DIL6CT  = 0;
  double TanBeta = 0, C_CosBeta = 0, C_SinBeta = 0, C_CosBetaSquared = 0,
         C_SinBetaSquared = 0;
  double beta             = 0;
  double M1 = 0, M2 = 0, M3 = 0, alpha1 = 0, alpha2 = 0, alpha3 = 0, MHC = 0;
  int Type       = 0;
  double CTempC1 = 0, CTempC2 = 0, CTempCS = 0;
  double R_Hh_1 = 0, R_Hh_2 = 0, R_Hh_3 = 0, R_Hl_1 = 0, R_Hl_2 = 0, R_Hl_3 = 0,
         R_Hsm_1 = 0, R_Hsm_2 = 0, R_Hsm_3 = 0;

  std::size_t pos_G0, pos_Gp, pos_Gm, pos_Hp, pos_Hm, pos_h1, pos_h2, pos_h3;
  std::size_t pos_h_SM, pos_h_l, pos_h_H;

  void ReadAndSet(const std::string &linestr,
                  std::vector<double> &par) override;
  std::vector<std::string> addLegendCT() const override;
  std::vector<std::string> addLegendTemp() const override;
  std::vector<std::string> addLegendTripleCouplings() const override;
  std::vector<std::string> addLegendVEV() const override;

  /**
   * Set the numerical values for the Lagrange parameters
   * @param par[0] = lambda_1
   * @param par[1] = lambda_2
   * @param par[2] = lambda_3
   * @param par[3] = lambda_4
   * @param par[4] = Re(lambda_5)
   * @param par[5] = Im(lambda_5)
   * @param par[6] = Re(m_{12}^2)
   * @param par[7] = tan(beta)
   * @param par[8] = Yukawa Type
   */
  void set_gen(const std::vector<double> &par) override;
  void set_CT_Pot_Par(const std::vector<double> &par) override;
  void write() const override;

  void AdjustRotationMatrix() override;
  void TripleHiggsCouplings() override;
  std::vector<double> calc_CT() const override;

  void SetCurvatureArrays() override;
  bool CalculateDebyeSimplified() override;
  bool CalculateDebyeGaugeSimplified() override;
  double VTreeSimplified(const std::vector<double> &v) const override;
  double VCounterSimplified(const std::vector<double> &v) const override;
  void Debugging(const std::vector<double> &input,
                 std::vector<double> &output) const override;

  bool IncludeChargeBreakingVEV = true;
  bool CTAlternative            = false;
  double DiffDelta              = 0;

  int PosSM = -1;
};

} // namespace Models
} // namespace BSMPT
#endif /* SRC_CLASSPOTENTIALC2HDM_H_ */
