.. _cxsm:

CxSM
==============

This implementation was based on [`1512.05355 <https://arxiv.org/abs/1512.05355>`_]. 

On the scalar sector there are a doublet and a complex scalar field

.. math::
   H=\frac{1}{\sqrt{2}}\binom{G^{+}}{v+h+i G^0} \quad \text { and } \quad \mathbb{S}=\frac{1}{\sqrt{2}}\left[v_S+s+i\left(v_A+a\right)\right]

The potential is given by

.. math::
   V_{\mathrm{CxSM}}=\frac{m^2}{2} H^{\dagger} H+\frac{\lambda}{4}\left(H^{\dagger} H\right)^2+\frac{\delta_2}{2} H^{\dagger} H|\mathbb{S}|^2+\frac{b_2}{2}|\mathbb{S}|^2+\frac{d_2}{4}|\mathbb{S}|^4+\left(\frac{b_1}{4} \mathbb{S}^2+a_1 \mathbb{S}+c . c .\right)

where all parameters are real except :math:`a_1` and :math:`b_1` (one of the phases can be absorbed by a phase rotation).

The :math:`T=0` vacuum is chosen to be

.. math::
   \left.\left\langle\Phi_1\right\rangle\right|_{T=0}=\frac{1}{\sqrt{2}}\binom{0}{v_1},\left.\quad\left\langle\mathbb{S}\right\rangle\right|_{T=0}=\frac{1}{\sqrt{2}}(v_s+i v_a)\,.

