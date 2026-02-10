.. _description:

######################
Theoretical Background
######################

The seismic hazard assessment conducted by the Ha3Py algorithm employs a hybrid approach.
The maximum magnitude :math:`m_{max}` is estimated separately using various methods.
In contrast other parameters of the magnitude exceed probability,
such as :math:`\beta` (the equivalent of the b-value in the Gutenberg-Richter (G-R) magnitude distribution)
and the :math:`\lambda` coefficient of the Poisson events occurrence probability,
are estimated independently using maximum-likelihood methods.
Nevertheless, the estimation errors for all parameters are assessed using the maximum likelihood method.

The earthquake recurrence parameters and their estimation
#########################################################

The sought area-characteristic earthquake recurrence parameters (ERP)
:math:`\mathbf{\theta}=\left( \lambda, \beta, m_{max} \right)` are estimated by the maximum likelihood method.
The chosen assessment procedure requires knowledge of the likelihood function of the parameters
:math:`\mathbf{\theta}=\left( \lambda, \beta, m_{max} \right)`.
Following the multiplicative property of the likelihood function (:cite:t:`Rao1973`),
the joint likelihood :math:`\mathcal{L}\left(\mathbf{\theta}\right)` based on prehistoric,
historic, and complete parts of the catalogue, is of the form:

.. math::
    :label: eq_a1

    \mathcal{L}\left(\mathbf{\theta}\right)=
    \mathcal{L}_P{\left(\mathbf{\theta}\right)}^{w_P}\times
    \mathcal{L}_H{\left(\mathbf{\theta}\right)}^{w_H}\times
    \prod_{i=1}^{N_C}{\mathcal{L}_C^{\left(i\right)}{\left(\mathbf{\theta}\right)}^{w_{Ci}}},

where :math:`\mathcal{L}_P{\left(\mathbf{\theta}\right)}`, :math:`\mathcal{L}_H{\left(\mathbf{\theta}\right)}`,
and :math:`\mathcal{L}_C^{\left(i\right)}{\left(\mathbf{\theta}\right)}`
denote the likelihood functions based on prehistoric, historic,
and complete sections of the catalogue :numref:`(Fig. %s) <fig_g1>`,
:math:`N_C` is the number of complete catalogues
and :math:`w_P,w_H, w_{Ci}` are weights of prehistoric (P), historic (H) and complete (Ci) sub-catalogues.
The form of the likelihood functions :math:`\mathcal{L}_P{\left(\mathbf{\theta}\right)}`,
:math:`\mathcal{L}_H{\left(\mathbf{\theta}\right)}`,
and :math:`\mathcal{L}_C^{\left(i\right)}{\left(\mathbf{\theta}\right)}` is provided,
e.g. by :cite:t:`Kijko_atal_2016` and :cite:t:`Smit_etal_2019`.


.. _description_m_max:

Maximum magnitude assessment methods
####################################

There is no universally accepted method for estimating the value of :math:`m_{max}`.
The presented software evaluates :math:`m_{max}` by nine different methods.
(For details, :cite:t:`Kijko_2004`; :cite:t:`KijkoSingh2011`; :cite:t:`VermeulenKijko2017` and :cite:t:`Kijko2025`).
However, the assessment of this parameter is only possible,
if information about it is available and provided by seismic event catalogues or other independent sources.

Primitive Procedure
===================

.. math::
    :label: eq_b1

    m_{max} = m_{max}^{obs} + 0.5


Robson-Whitlock Procedure
=========================

:cite:t:`Robson_and_Whitlock_1964` showed that, under very general conditions,
when magnitudes of :math:`n` events are arranged in ascending order
:math:`m_1 \le m_2 \le .... \le m_{n-1} \le m_{max}^{obs}`,
the estimation of :math:`m_{max}` takes the form:

.. math::
    :label: eq_b2

    \widehat{m}_{max} = m_{max}^{obs}+\left( m_{max}^{obs} - m_{n-1} \right)

with variance

.. math::
    \text{VAR}\left( \widehat{m}_{max} \right)=
    5\sigma_M^2 + \left( m_{max}^{obs} - m_{n-1} \right)^2,

where :math:`\sigma_M^2` denotes the standard error
of the largest observed magnitude :math:`m_{max}^{obs}` determination.

Robson-Whitlock-Cooke Procedure
===============================

:cite:t:`Cooke_1979` showed that a reduction in the mean squared error of the Robson-Whitlock estimator is possible
when some information about the shape of the upper tail
of the probability distribution function :math:`f_M \left( m \right)` is known.
Assuming that the observed magnitudes are sampled from a distribution that is truncated,
as in the double truncated G-R relation (:cite:t:`GutenbergRichter1956`),
the improved version of the Robson-Whitlock estimator is

.. math::
    :label: eq_b3

    \widehat{m}_{max} = m_{max}^{obs}+0.5\left( m_{max}^{obs} - m_{n-1} \right)

with the variance

.. math::
    \text{VAR}\left( \widehat{m}_{max} \right)=
    1.5\sigma_M^2 + 0.25\left( m_{max}^{obs} - m_{n-1} \right)^2


Note that the *Primitive*, *Robson-Whitlock*, and *Robson-Whitlock-Cooke* procedures
do not require specification of the magnitude distribution.

Gibowicz-Kijko Procedure
========================

The procedure relies on the properties of end-point estimators of the uniform distribution.
Since the values of the magnitude distribution :math:`F_M \left( m \right)` follow a uniform distribution,
the following relation can be anticipated

.. math::
    :label: eq_b3a

    F_M \left( m_{max}^{obs} \right)= \frac{n}{n+1} .

Since :math:`F_M \left( m \right)` depends on :math:`m_{max}`,
we can derive the estimator of the :math:`\widehat{m}_{max}` by solving equation :eq:`eq_b3a`.

Tate-Pisarenko Procedure
========================

The Tate-Pisarenko (and following Kijko-Sellevoll) assess :math:`m_{max}` by solving the equation

.. math::
    :label: eq_b4

    \widehat{m}_{max} = m_{max}^{obs}+\Delta,

where :math:`\Delta` depends on :math:`m_{max}`.
In the case of applying the Tate-Pisarenko procedure,
the correction factor :math:`\Delta` takes the form

.. math::
    :label: eq_b5

    \Delta=\frac{1}{nf_M\left( m_{max}^{obs} \right)},

where the probability density function of the magnitude distribution :math:`f_M\left( m \right)`
depends on :math:`m_{max}`. The approximate variance of the Tate-Pisarenko estimator is of the form
(:cite:t:`Kijko_and_Graham_1998`; :cite:t:`Kijko_2004`)

.. math::
    \text{VAR}\left( \widehat{m}_{max} \right)=
    \sigma_M^2+\frac{n+1}{n^3f_M^2\left( m_{max}^{obs} \right)}

where :math:`\sigma_M^2` denotes the standard error of the largest observed magnitude determination.

Kijko-Sellevoll procedure
=========================

Similar to the previous Tate-Pisarenko procedure :eq:`eq_b4`,
the Kijko-Sellevoll procedure assesses :math:`m_{max}` by solving the equation

.. math::
    \widehat{m}_{max} = m_{max}^{obs}+\Delta,

where

.. math::
    :label: eq_b6

    \Delta=\int_{m_{min}}^{m_{max}}F_{M}\left(m\right)^{n}dm,

The approximate variance of the Kijko-Sellevoll estimator of :math:`m_{max}`
for the frequency-magnitude G-R distribution is of the form

.. math::
    \text{VAR}\left( \widehat{m}_{max} \right) = \sigma_M^2+\Delta^2

Procedure based on the largest few earthquakes
==============================================

The procedure based on the largest few earthquakes is a special case of the Kijko-Sellevoll procedure.
The non-parametric magnitude distribution model is applied,
using only the largest earthquakes for model building.

Bayesian maximum magnitude assessment methods
#############################################

The Bayesian :math:`m_{max}` assessment methods incorporate any relevant information about :math:`m_{max}`.
This information is external and can come from geology, tectonics, or the seismicity of similar regions.
The idea of :math:`m_{max}` Bayesian estimation was first described by :cite:`Cornell1994`.
It combines two sources of information:

*   The information of the :math:`m_{max}` prior distribution is described as :math:`\pi\left(m_{max}\right)`.
    The information is based on observations and comes from external sources.
    Following :cite:`Coppersmith1994`, it is assumed that :math:`\pi\left(m_{max}\right)`
    has the form of a double-truncated Gaussian distribution
    :math:`\pi\left(m_{max}\right)=\mathcal{N}\left(m_{max}^{prior},\sigma_{m_{max}^{prior}}, m_{max}^L, m_{max}^U \right)`,
    where :math:`m_{\max}^U` is the largest magnitude that might ever happen,
    and :math:`m_{max}^L` is the magnitude that we are confident has occurred.
    Usually, it is the maximum observed magnitude.
*   The likelihood function :math:`\mathcal{L}\left(\mathbf{m}|m_{max}\right)` of observed earthquake magnitudes.

Both the prior knowledge about :math:`m_{max}` and the knowledge derived from the observed magnitudes :math:`\mathbf{m}`,
which is known as the posterior distribution of :math:`m_{max}`, are summarised in the form

.. math::
    :label: eq_b7

    p_{m_{max}}\left(m_{max}|\mathbf{m}\right)=
    \begin{cases}
    0 & : m_{max} < m_{\max}^L \\
    C\cdot \pi\left(m_{max}\right)\mathcal{L}\left(\mathbf{m}|m_{max}\right)
    & : m_{\max}^L \leqslant m_{max} \leqslant m_{\max}^U \\
    0 & : m_{max} > m_{\max}^U
    \end{cases},

where :math:`C` is a normalising constant

.. math::
    C=1/\int_{m_{max}^L}^{m_{max}^U}
    \pi\left(m_{max}\right)L\left(\mathbf{m}|m_{max}\right)dm_{max}.

Based on :eq:`eq_b7`, three Bayesian analogues of the maximum likelihood (ML) point estimators are used:

*   The maximum posterior estimate (MAP) value

    .. math::
        p_{m_{max}}\left({\widehat{m}}_{max}^{posterior}|\mathbf{m}\right)=
        \text{maximum}

*   The posterior mean (PM) value (expected value)

    .. math::
        {\widehat{m}}_{max}^{posterior}=
        \int{\zeta\ p_{m_{max}}\left(\zeta\right)}d\zeta

*   The posterior median is defined by the solution of the equation

    .. math::
        \int_{m_{\max}^L}^{{\widehat{m}}_{max}^{posterior}}
        {p_{m_{max}}\left( \zeta \right)}d\zeta=\frac{1}{2}

:cite:`Cornell1994` proposed

.. math::
    :label: eq_b8

    \mathcal{L}\left(\mathbf{m}|m_{max}\right)=
    \prod_{i=1}^{n}{f_M\left(m_i|m_{max}\right)}.

However, if :math:`f_M\left(m|m_{max}\right)` is the G-R distribution,
the MAP that fulfills :eq:`eq_b8` gives the maximum observed magnitude :math:`m_{\max}^{obs}` as
the solution of :eq:`eq_b7`, which is not expected as the maximum possible magnitude :math:`m_{max}`.
Therefore, several techniques were applied to correct the :cite:`Cornell1994` procedure.

Bayesian :math:`m_{max}` assessment based on the shift of the likelihood function
=================================================================================

In this method :math:`m_{max}^L={\widehat{m}}_{max}`,
where :math:`{\widehat{m}}_{max}` can be assessed by any of those described in the
:ref:`Maximum magnitude assessment methods` section.
The likelihood function is

.. math::
    :label: eq_b9

    \mathcal{L}\left(\mathbf{m}|m_{max}\right)=
    \prod_{i=1}^{n}{f_M\left(m_i | m_{max}\right)}.

Bayesian :math:`m_{max}` assessment based on the Gaussian distribution
======================================================================

In this method :math:`m_{max}^L=m_{max}^{obs}`,
and the likelihood function is

.. math::
    :label: eq_b10

    \mathcal{L}\left(\mathbf{m}|m_{max}\right)=
    \mathcal{N}\left({\widehat{m}}_{max},\sqrt{\text{VAR}\left( \widehat{m}_{max} \right)} \right)

where :math:`{\widehat{m}}_{max}` and :math:`\text{VAR}\left( \widehat{m}_{max} \right)`
can be assessed by any of methods described in the :ref:`Maximum magnitude assessment methods` section.

Bayesian :math:`m_{max}` assessment based on the Fiducial :math:`m_{max}` distribution
======================================================================================

The method assumes that the database information on :math:`m_{max}`
(in our case, seismic event catalogue)
is expressed in the form of the fiducial distribution :cite:`Pisarenko1991`

.. math::
    F_{M_{max}}^{FID}(m_{max})=
    1-\ \left[F_M\left(m_{max}^{obs}|m_{max} \right)\right]^n

where :math:`n` is the number of earthquakes of magnitude :math:`m\geq m_C`.
The probability density function is

.. math::
    :label: eq_b11

    f_{M_{max}}^{FID} \left( m_{max} \right)=
    n\left[F_M\left(m_{max}^{obs}|m_{max} \right)\right]^{n-1}
    \frac{\partial S_M\left(m_{max}^{obs}|m_{max} \right)}{\partial m_{max}}

and the likelihood :math:`\mathcal{L}\left(\mathbf{m}|m_{max}\right)=f_{M_{max}}^{FID} \left( m_{max} \right)`.
The equation :eq:`eq_b11` is written in the form that uses methods defined in
the magnitude distribution classes (see :ref:`classes_rev`).

Estimation uncertainty assessment
#################################

.. math::
    \mathbf{COV_\Theta} = \mathbf{H_\mathcal{L}}^{-1},

where :math:`\mathbf{\Theta}` is the vector estimated by the maximum likelihood coefficients,
and :math:`\mathbf{H_\mathcal{L}}` is the Hessian of the likelihood logarithm

.. math::
    \mathbf{H_\mathcal{L}}=
    \begin{bmatrix}
    \frac{\partial^2 \ln\left(\mathcal{L}\right)}{\partial \Theta_1 \partial \Theta_1} &
    \cdots &
    \frac{\partial^2 \ln\left(\mathcal{L}\right)}{\partial \Theta_1 \partial \Theta_K} \\
    \vdots & \ddots & \vdots \\
    \frac{\partial^2 \ln\left(\mathcal{L}\right)}{\partial \Theta_K \partial \Theta_1} &
    \cdots &
    \frac{\partial^2 \ln\left(\mathcal{L}\right)}{\partial \Theta_K \partial \Theta_K} \end{bmatrix}

The standard deviation of the earthquake occurrence probability coefficients is

.. math::
    \sigma_\mathbf{\Theta} = \sqrt{\text tr\left( \mathbf{COV_\Theta} \right)}

* Example:
    When the  coefficients are :math:`\beta,\lambda`,

    .. math::
        \mathbf{COV}_{\beta,\lambda}=
        \left[\begin{matrix}-\frac{\partial^{2}\ln\left(\mathcal{L}\right)}{\partial\lambda^{2}} &
        -\frac{\partial^{2}\ln\left(\mathcal{L}\right)}{\partial\lambda\partial\beta} \\
        -\frac{\partial^{2}\ln\left(\mathcal{L}\right)}{\partial\lambda\partial\beta} &
        -\frac{\partial^{2}\ln\left(\mathcal{L}\right)}{\partial\beta^{2}}
        \end{matrix}\right]^{-1},

    where

    .. math::
        \sigma_\lambda = \mathbf{COV}_{\beta,\lambda}[1,1]

        \sigma_\beta = \mathbf{COV}_{\beta,\lambda}[2,2]

Return period computing methods
###############################

The return period (:math:`T_R \left( m \right)`) is a function of magnitude.
It is defined as the inverse of the annual probability of not occurrence,
which is the survival function of the magnitude distribution.

.. math::
    :label: eq_c1

    T_{R}\left(m,\lambda_{0},F_{M}\right)=\frac{1}{S_M^{max}\left(m\right)}

Sometimes, the annual probability of not-occurrence is described as :math:`\lambda`
(Not to confuse with the lambda notation as the occurrence probability coefficient).
Then, e.g., the return period (:math:`T_R \left( m \right)`) is

.. math::

    T_R\left(m \right)=\frac{1}{\lambda\left(m\right)}

where :math:`\lambda\left( m \right)=\lambda_0S_M\left( m|m_0 \right)`,
:math:`m_0` is the assumed minimum magnitude,
:math:`\lambda_0` is the annual probability of not-occurrence corresponding
to the minimum magnitude,
and :math:`S_{M}\left(m|m_{0}\right)=1-F_{M}\left(m|m_{0}\right)` is the survival function
of the :ref:`Magnitude distribution <api_md>`.

The return period uncertainty is

.. math::
    :label: eq_c2

    \left(\nabla_\mathbf\Theta T_{R}\right)^{T}\mathbf{COV_\Theta}\nabla_\mathbf\Theta T_{R}

where

.. math::
    \frac{\partial T_R\left( m \right)}{\partial\Theta_i}=
    \frac{-1}{\left( S_M^{max}\left(m\right) \right)^2}
    \frac{\partial S_M^{max}\left(m\right)}{\partial\Theta_i}

The gradient of the survival function of magnitude distribution
depends on not occurrence distribution models
(see specific applications in :ref:`Seismic event occurrence probability <api_ed>` section).
In the typical simplified case, when the magnitude and not-occurrence distribution coefficients are
:math:`\lambda_0`, :math:`\beta` and :math:`m_{max}`, gradient of the return period is

.. math::
    :label: eq_c3

    \frac{\partial T_R\left( m \right)}{\partial\lambda_0}=
    \frac{\partial\left(\frac{1}{\lambda_{0}\cdot S_{M}\left(m|m_{0}\right)}\right)}
    {\partial\lambda_{0}}=\frac{-1}{\lambda_{0}^{2}\cdot S_{M}\left(m|m_{0}\right)},

    \frac{\partial T_{R}\left(m\right)}{\partial\beta}=
    \frac{\partial\left(\frac{1}{\lambda_{0}\cdot S_{M}\left(m|m_{0}\right)}\right)}
    {\partial\beta}=\frac{-1}{\lambda\cdot S_{M}^{2}\left(m|m_{0}\right)}
    \frac{\partial S_{M}\left(m|m_{0}\right)}{\partial\beta},

    \frac{\partial T_{R}\left(m\right)}{\partial m_{max}}=
    \frac{\partial\left(\frac{1}{\lambda_{0}\cdot S_{M}\left(m|m_{0}\right)}\right)}
    {\partial m_{max}}=\frac{-1}{\lambda\cdot S_{M}^{2}\left(m|m_{0}\right)}
    \frac{\partial S_{M}\left(m|m_{0}\right)}{\partial m_{max}},

where gradients of survival function :math:`\partial S_{M}\left(m|m_{0}\right) / \partial m_{max}`
are defined with magnitude distribution classes.

Synthetic catalogues
####################

Simulation of two types of synthetic catalogues is possible:

* Each seismic event is specified by its origin time and magnitude
* Each seismic event is specified by its magnitude

The first method should be applied for simulating paleo and historical catalogues.
In contrast, the second method creates complete catalogues.

The catalogue, which contains both the origin times and magnitudes created by the simulation,
with event times, can be used for all likelihood functions.
The catalogue only includes seismic event magnitudes can be used by the likelihood function.

Two methods can calculate the catalog of earthquakes containing origin times and magnitudes:

:Extreme events simulation: This method divides the catalog period
    into roughly equal time intervals with random margins,
    during which it counts a random number of earthquakes
    with random magnitudes following the specified occurrence probability.
    From these events, only one earthquake with extreme magnitude is selected and recorded in the catalog.
:Full simulation incremental: This method is iterative with incremental event time.
    First, the time is set to the catalog's starting time.
    In each iteration, the algorithm generates random periods of non-occurrence of the seismic event
    drawn from a defined probability of occurrence.
    The time is incremented by the period and set as the event time.
    The earthquake with a random magnitude based on the defined magnitude distribution is recorded in the catalog
    and the time for the next iteration is set to the event time.
    This process repeats until the end of the catalog is reached.

The catalog of earthquakes containing only origin magnitudes without origin time can be calculated by one method:

:Full simulation without date: The method generates a random number of earthquakes for the catalog period,
    assigns random magnitudes based on the specified occurrence probabilities,
    and then records all events with their magnitudes into the catalog.