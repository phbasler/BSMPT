.. _c2hdm:

C2HDM
==============

This implementation was based on [`2108.03580 <https://arxiv.org/abs/2108.03580>`_]. 

On the scalar sector there are two doublets 

.. math::
   \Phi_1=\frac{1}{\sqrt{2}}\binom{\rho_1+\mathrm{i} \eta_1}{\zeta_1+\omega_1+\mathrm{i} \psi_1}, \quad \Phi_2=\frac{1}{\sqrt{2}}\binom{\rho_2+\omega_{\mathrm{CB}}+\mathrm{i} \eta_2}{\zeta_2+\omega_2+\mathrm{i}\left(\psi_2+\omega_{\mathrm{CP}}\right)}

where we impose a softly-broken :math:`\mathbb{Z}_2` which takes

.. math::
   \begin{align}
   \Phi_{1} &\to -\Phi_{1}\,,\\
   \Phi_{2} &\to \Phi_{2}.
   \end{align}

The potential is given by

.. math::
   \begin{aligned}V_{\text {C2HDM}}= & m_{11}^2 \Phi_1^{\dagger} \Phi_1+m_{22}^2 \Phi_2^{\dagger} \Phi_2-\left[m_{12}^2 \Phi_1^{\dagger} \Phi_2+\text { h.c. }\right]+\frac{1}{2} \lambda_1\left(\Phi_1^{\dagger} \Phi_1\right)^2+\frac{1}{2} \lambda_2\left(\Phi_2^{\dagger} \Phi_2\right)^2 \\& +\lambda_3\left(\Phi_1^{\dagger} \Phi_1\right)\left(\Phi_2^{\dagger} \Phi_2\right)+\lambda_4\left(\Phi_1^{\dagger} \Phi_2\right)\left(\Phi_2^{\dagger} \Phi_1\right)+\left[\frac{1}{2} \lambda_5\left(\Phi_1^{\dagger} \Phi_2\right)^2+\text { h.c. }\right] .\end{aligned}

where all parameters are real except :math:`m_{12}` and :math:`\lambda_5` which, if non-real, explicitly break CP in the model at :math:`T = 0`.

The :math:`T=0` vacuum is chosen to be

.. math::
   \left.\left\langle\Phi_1\right\rangle\right|_{T=0}=\frac{1}{\sqrt{2}}\binom{0}{v_1},\left.\quad\left\langle\Phi_2\right\rangle\right|_{T=0}=\frac{1}{\sqrt{2}}\binom{0}{v_2}\,.

The Yukawa types defined by
   * Type I : :math:`\Phi_d = \Phi_l = \Phi_2`
   * Type II: :math:`\Phi_d = \Phi_l = \Phi_1`
   * Type  LS(=3): :math:`\Phi_d = \Phi_2\,, \Phi_l = \Phi_1`
   * Type FL(=4): :math:`\Phi_d = \Phi_1 \,,\Phi_l = \Phi_2`
