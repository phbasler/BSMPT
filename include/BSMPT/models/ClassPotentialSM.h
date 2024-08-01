// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * Model file for the SM
 */

#pragma once

#include <BSMPT/models/ClassPotentialOrigin.h>

namespace BSMPT
{
namespace Models
{

/**
 * @brief The Class_SM class
 * Implementation of the Standard Model
 *
 *
 * * \f$ -L_S = muSq \Phi^\dagger \Phi + lambda (\Phi^\dagger \Phi)^2 \f$
 * * \f$ -L_Y = \overline{U}_L V \text{diag}(y_d,y_s,y_b) D_R \Phi^+  +
 * \overline{D}_L \text{diag}(y_d,y_s,y_b) D_R \Phi^0 + \overline{U}_L
 * \text{diag}(y_u, y_c, y_t) U_R (\Phi^0)^\ast - \overline{D}_L V^\dagger
 * \text{diag}(y_u, y_c, y_t) U_R \Phi^- + \overline{E}_L \text{diag}(y_e,
 * y_\mu, y_\tau) E_R \Phi^0 +\overline{\nu}_L \text{diag}(y_e, y_\mu, y_\tau)
 * E_R \Phi^+ + c.c.  \f$
 * * \f$ L_G = (D_\mu \Phi)^\dagger (D^\mu \Phi) \f$
 *
 * with
 *
 * * \f$ \Phi = \begin{pmatrix} \Phi^+ \\  \Phi^0 \end{pmatrix} =
 * 1/\sqrt{2} \begin{pmatrix} higgsbase[0] + I*higgsbase[1] \\  higgsbase[2] + I
 * *higgsbase[3] \end{pmatrix}\,, \f$
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
 */
class Class_SM : public Class_Potential_Origin
{
public:
  Class_SM(const ISMConstants &smConstants);
  virtual ~Class_SM();

  double muSq, lambda;

  double dmuSq, dlambda, dT1, dT2, dT3, dT4;

  double v0;

  void ReadAndSet(const std::string &linestr,
                  std::vector<double> &par) override;
  std::vector<std::string> addLegendCT() const override;
  std::vector<std::string> addLegendTemp() const override;
  std::vector<std::string> addLegendTripleCouplings() const override;
  std::vector<std::string> addLegendVEV() const override;

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
