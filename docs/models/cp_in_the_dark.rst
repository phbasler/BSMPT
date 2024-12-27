.. _cp_in_the_dark:

CP in the dark
==============

This implementation was based on [`1807.10322 <https://arxiv.org/abs/1807.10322>`_] and [`2204.13425 <https://arxiv.org/abs/2204.13425>`_]. 

On the scalar sector there are two doublets 

.. math::
   \Phi_1=\frac{1}{\sqrt{2}}\binom{\rho_1+\mathrm{i} \eta_1}{\zeta_1+\omega_1+\mathrm{i} \psi_1}, \quad \Phi_2=\frac{1}{\sqrt{2}}\binom{\rho_2+\omega_{\mathrm{CB}}+\mathrm{i} \eta_2}{\zeta_2+\omega_2+\mathrm{i}\left(\psi_2+\omega_{\mathrm{CP}}\right)}, \quad \Phi_S=\zeta_3+\omega_S

where we impose a softly-broken :math:`\mathbb{Z}_2` symmetry which takes

.. math::
   \begin{align}
   \Phi_{1} &\to -\Phi_{1}\,,\\
   \Phi_{2} &\to \Phi_{2}\,,\\
   \Phi_S &\to -\Phi_S\,,
   \end{align}

The potential is given by

.. math::
    \begin{aligned}V_\text{CP in the Dark} &= m_{11}^2\left|\Phi_1\right|^2+m_{22}^2\left|\Phi_2\right|^2+\frac{1}{2} m_S^2 \Phi_S^2+\left(A \Phi_1^{\dagger} \Phi_2 \Phi_S+\text { h.c. }\right) \\& +\frac{1}{2} \lambda_1\left|\Phi_1\right|^4+\frac{1}{2} \lambda_2\left|\Phi_2\right|^4+\lambda_3\left|\Phi_1\right|^2\left|\Phi_2\right|^2+\lambda_4\left|\Phi_1^{\dagger} \Phi_2\right|^2+\frac{1}{2} \lambda_5\left[\left(\Phi_1^{\dagger} \Phi_2\right)^2+\text { h.c. }\right] \\& +\frac{1}{4} \lambda_6 \Phi_S^4+\frac{1}{2} \lambda_7\left|\Phi_1\right|^2 \Phi_S^2+\frac{1}{2} \lambda_8\left|\Phi_2\right|^2 \Phi_S^2,\end{aligned}

where all coupling constants are real except :math:`A`, whose complex phase is directly linked with CP violation.

The :math:`T=0` vacuum is chosen to be

.. math::
   \left.\left\langle\Phi_1\right\rangle\right|_{T=0}=\frac{1}{\sqrt{2}}\binom{0}{v_1},\left.\quad\left\langle\Phi_2\right\rangle\right|_{T=0}=\frac{1}{\sqrt{2}}\binom{0}{v_2},\left.\quad\langle S\rangle\right|_{T=0}=v_S

The Yukawa types defined by
   * Type I : :math:`\Phi_d = \Phi_l = \Phi_2`
   * Type II: :math:`\Phi_d = \Phi_l = \Phi_1`
   * Type  LS(=3): :math:`\Phi_d = \Phi_2\,, \Phi_l = \Phi_1`
   * Type FL(=4): :math:`\Phi_d = \Phi_1 \,,\Phi_l = \Phi_2`
