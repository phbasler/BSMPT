// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SRC_CLASSPOTENTIALRN2HDM_H_
#define SRC_CLASSPOTENTIALRN2HDM_H_

#include <BSMPT/models/ClassPotentialOrigin.h>

/**
 * @file
 */

namespace BSMPT
{
namespace Models
{

/**
 * @brief The Class_Potential_RN2HDM class
 * Implementation of the real N2HDM, as shown in the manual
 *
 * * \f$ -L_S = u1 \Phi_1^\dagger \Phi_1 + u2 \Phi_2^\dagger \Phi_2 +
 * \frac{L1}{2} (\Phi_1^\dagger \Phi_1)^2 + \frac{L2}{2} (\Phi_2^\dagger
 * \Phi_2)^2 + L3 \Phi_1^\dagger \Phi_1 \Phi_2^\dagger \Phi_2 + L4
 * \Phi_1^\dagger \Phi_2 \Phi_2^\dagger \Phi_1 + [ \frac{RL5}{2} (\Phi_1^\dagger
 * \Phi_2)^2 - RealMMix \Phi_1^\dagger \Phi_2 ] + Nus/2 S^2 + NL6/8 S^4 + NL7/2
 * \Phi_1^\dagger \Phi_1 S^2 + NL8/2 \Phi_2^\dagger \Phi_2 S^2 \f$
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
 * 1/\sqrt{2} \begin{pmatrix} higgsbase[0] + I*higgsbase[2] \\  higgsbase[6] + I
 * *higgsbase[4] \end{pmatrix}\,, \f$
 * * \f$ \Phi_2 = \begin{pmatrix} \Phi_2^+ \\  \Phi_2^0 \end{pmatrix} =
 * 1/\sqrt{2} \begin{pmatrix} higgsbase[1] + I*higgsbase[3] \\  higgsbase[7] + I
 * *higgsbase[5] \end{pmatrix}  \f$
 * * \f$ S = higgsbase[8] \f$
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
 *
 */
class Class_Potential_RN2HDM : public Class_Potential_Origin
{
public:
  Class_Potential_RN2HDM(const ISMConstants &smConstants);
  virtual ~Class_Potential_RN2HDM();

  double L1 = 0, L2 = 0, L3 = 0, L4 = 0, RL5 = 0, RealMMix = 0, u1 = 0, u2 = 0;
  double DL1CT = 0, DL2CT = 0, DL3CT = 0, DL4CT = 0, DRL5CT = 0, Du2CT = 0,
         Du1CT = 0, DRu3CT = 0;
  double DT1 = 0, DT2 = 0, DT3 = 0;
  double TanBeta = 0, C_CosBeta = 0, C_SinBeta = 0, C_CosBetaSquared = 0,
         C_SinBetaSquared = 0;
  double beta             = 0;
  int Type                = 0;
  double CTempC1 = 0, CTempC2 = 0, CTempCS = 0;
  double alpha1 = 0, alpha2 = 0, alpha3 = 0;
  double MSM = 0, MhUp = 0, MhDown = 0;

  double Nus = 0, NL6 = 0, NL7 = 0, NL8 = 0, Nvs = 0;
  double NDus = 0, NDL6 = 0, NDL7 = 0, NDL8 = 0, NDvs = 0, NDTS = 0;
  double DTCharged = 0;

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
   * @param par[4] = lambda_5
   * @param par[5] = lambda_6
   * @param par[6] = lambda_7
   * @param par[7] = lambda_8
   * @param par[8] = tan(beta)
   * @param par[9] = v_s
   * @param par[10] = m_{12}^2
   * @param par[11] = Yukawa Type
   */
  void set_gen(const std::vector<double> &par) override;
  void set_CT_Pot_Par(const std::vector<double> &par) override;
  void write() const override;

  void TripleHiggsCouplings() override;
  std::vector<double> calc_CT() const override;

  void SetCurvatureArrays() override;
  bool CalculateDebyeSimplified() override;
  bool CalculateDebyeGaugeSimplified() override;
  double VTreeSimplified(const std::vector<double> &v) const override;
  double VCounterSimplified(const std::vector<double> &v) const override;
  void Debugging(const std::vector<double> &input,
                 std::vector<double> &output) const override;
};

} // namespace Models
} // namespace BSMPT
#endif /* SRC_CLASSPOTENTIALRN2HDM_H_ */
