.. _n2hdm:

N2HDM
==============

This implementation was based on [`1612.01309 <https://arxiv.org/abs/1612.01309>`_] and [`1912.10477 <https://arxiv.org/abs/1912.10477>`_]. 

On the scalar sector there are two doublets 

.. math::
   \Phi_1=\frac{1}{\sqrt{2}}\binom{\rho_1+\mathrm{i} \eta_1}{\zeta_1+\omega_1+\mathrm{i} \psi_1}, \quad \Phi_2=\frac{1}{\sqrt{2}}\binom{\rho_2+\omega_{\mathrm{CB}}+\mathrm{i} \eta_2}{\zeta_2+\omega_2+\mathrm{i}\left(\psi_2+\omega_{\mathrm{CP}}\right)}, \quad \Phi_S=\zeta_3+\omega_S

where we impose two softly-broken :math:`\mathbb{Z}_2` symmetries. One :math:`\mathbb{Z}_2` which takes

.. math::
   \begin{align}
   \Phi_{1} &\to -\Phi_{1}\,,\\
   \Phi_{2} &\to \Phi_{2}\,,\\
   \Phi_S &\to \Phi_S\,,
   \end{align}

and another :math:`\mathbb{Z}_2` which takes

.. math::
   \begin{align}
   \Phi_{1} &\to \Phi_{1}\,,\\
   \Phi_{2} &\to \Phi_{2}\,,\\
   \Phi_S &\to -\Phi_S\,.
   \end{align}

The potential is given by

.. math::
   \begin{aligned}V_{\mathrm{N} 2 \mathrm{HDM}}= & m_{11}^2 \Phi_1^{\dagger} \Phi_1+m_{22}^2 \Phi_2^{\dagger} \Phi_2-m_{12}^2\left(\Phi_1^{\dagger} \Phi_2+\text { h.c. }\right)+\frac{\lambda_1}{2}\left(\Phi_1^{\dagger} \Phi_1\right)^2+\frac{\lambda_2}{2}\left(\Phi_2^{\dagger} \Phi_2\right)^2 \\& +\lambda_3 \Phi_1^{\dagger} \Phi_1 \Phi_2^{\dagger} \Phi_2+\lambda_4 \Phi_1^{\dagger} \Phi_2 \Phi_2^{\dagger} \Phi_1+\frac{\lambda_5}{2}\left(\left(\Phi_1^{\dagger} \Phi_2\right)^2+\text { h.c. }\right) \\& +\frac{1}{2} m_S^2 \Phi_S^2+\frac{\lambda_6}{8} \Phi_S^4+\lambda_7\left(\Phi_1^{\dagger} \Phi_1\right) \Phi_S^2+\lambda_8\left(\Phi_2^{\dagger} \Phi_2\right) \Phi_S^2,\end{aligned}

where all coupling constants are real to preserve CP in the model at :math:`T = 0`.

The :math:`T=0` vacuum is chosen to be

.. math::
   \left.\left\langle\Phi_1\right\rangle\right|_{T=0}=\frac{1}{\sqrt{2}}\binom{0}{v_1},\left.\quad\left\langle\Phi_2\right\rangle\right|_{T=0}=\frac{1}{\sqrt{2}}\binom{0}{v_2},\left.\quad\langle S\rangle\right|_{T=0}=v_S

The Yukawa types defined by
   * Type I : :math:`\Phi_d = \Phi_l = \Phi_2`
   * Type II: :math:`\Phi_d = \Phi_l = \Phi_1`
   * Type  LS(=3): :math:`\Phi_d = \Phi_2\,, \Phi_l = \Phi_1`
   * Type FL(=4): :math:`\Phi_d = \Phi_1 \,,\Phi_l = \Phi_2`
