"""
Gamma compound Poisson events occurrence probability
----------------------------------------------------

:copyright:
    Jan Wiszniowski <jwisz@igf.edu.pl>,
    Andrzej Kijko <andrzej.kijko@up.ac.za>
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
:version 0.0.1:
    2025-01-01
"""

from abc import ABC
import numpy as np
from scipy.special import factorial, gamma, loggamma
from ha3py.LambdaOccurrence import LambdaOccurrence
from ha3py.BaseOccurrence import get_event_occurrence_parameters


class PoissonGammaCompoundOccurrence(LambdaOccurrence, ABC):
    r"""
    It is the combination of Poisson distribution with the gamma distribution
    to create the Poisson-gamma compound distribution to obtain the probability
    to observe :math:`n` seismic events, within a time interval :math:`t`,
    for temporal varying seismic activity :math:`\lambda` (Benjamin, 1968), as follows

    .. math::
        p_n(n|\lambda,t)=
        \frac{\Gamma\left( n+q_\lambda \right)}{n!\Gamma\left( q_\lambda \right)}
        \left( \frac{p_\lambda}{t+p_\lambda} \right)^{q_\lambda}
        \left( \frac{t}{t+p_\lambda} \right)^n

    where in which :math:`\Gamma()` is the gamma function, :math:`q_\lambda={\lambda}^2/{\sigma_\lambda}^2`
    and :math:`p_\lambda=q_\lambda/\lambda` are the constant parameters of gamma distribution,
    and :math:`\lambda` means the mean value of the activity rate.
    This calculation is performed as

        .. math::
            p_n(n|\lambda,t)= \exp\left( sum  \right),

        where

        .. math::
            sum = \ln \Gamma\left( n+q_\lambda \right)
            - \ln \Gamma\left( n+1 \right)
            - \ln \Gamma\left( q_\lambda \right)

            + q_\lambda\ln\left( \frac{p_\lambda}{t+p_\lambda} \right)
            + n\ln\left( \frac{t}{t+p_\lambda} \right).

    The PDF of Poisson-gamma compound distribution of not exceeding magnitude :math:`m`
    in time :math:`t` is

    .. math::
        f_M^{max}\left( m|\lambda,t \right)=
        \frac{\lambda t q_\lambda f_M\left( m \right)F_M^{max}\left( m|\lambda,t \right)}
        {q_\lambda+\lambda t S_M\left( m \right)}

    The CDF of Poisson-gamma compound distribution of not exceeding magnitude :math:`m`
    in time :math:`t` is

    .. math::
        F_M^{max}\left( m|\lambda,t \right)
        \left[ \frac{q_\lambda}{q_\lambda+\lambda t S_M\left( m \right)} \right]^{q_\lambda}

    Gradients of the of non occurrence of event in the time probability survive function are:

    .. math::
        \frac{\partial S_M^{max}\left(m | t\right)}{\partial\lambda}=
        \left(1+\frac{\lambda t \left\{S_M\left(m \right)\right\} }{q_\lambda} \right)^{-q_\lambda-1}
        t\left\{S_M\left(m \right)\right\}

        \frac{\partial S_M^{max}\left(m | t\right)}{\partial\theta_i}=
        t\lambda\left(1+\frac{ \lambda t\left\{S_M\left(m \right)\right\}}{q_\lambda}\right)^{-q_\lambda-1}
        \left\{\frac{\partial S_M\left(m \right)}{\partial\theta_i}\right\}

    Gradients of the :math:`n` events in the time in the time probability survive function are:

    Poisson-gamma compound coefficient are: :math:`\lambda` ('lambda') and magnitude distribution coefficient.
    One constant coefficient is :math:`q_\lambda` ('q_lambda'). The class overwrite the OccurrenceBase methods:

    `_grad_sf_magnitude_distribution`:

    `_d_grad_sf_magnitude_distribution`:

    `_grad_sf`:

    `_d_grad_sf`:

    `d_mean`:

    """

    def __init__(self, configuration, **kwargs):
        self.q_lambda = configuration.get('q_lambda', 100.0)
        super().__init__(configuration, 'Poisson-gamma compound', **kwargs)

    def factor(self, n):
        result = 1.0
        for idx in range(n):
            result *= (self.q_lambda+n-1) / n
        return result

    # def _d_pmf1(self, k, t):
    #     p_lambda = self.q_lambda / self.lamb
    #     den = t + p_lambda
    #
    #         component1 = loggamma(k + self.q_lambda) - loggamma(self.q_lambda) - loggamma(k + 1)
    #         component2 = np.log(p_lambda / den) * self.q_lambda + np.log(t / den) * k
    #         return np.exp(component1 + component2)
    #     else:
    #         return (p_lambda / den) ** self.q_lambda

    def _d_pmf(self, n, t):
        r"""
        Function

        .. math::
            p_n(n|\lambda,t)=
            \frac{\Gamma\left( n+q_\lambda \right)}{n!\Gamma\left( q_\lambda \right)}
            \left( \frac{p_\lambda}{t+p_\lambda} \right)^{q_\lambda}
            \left( \frac{t}{t+p_\lambda} \right)^n

        is realized as

        .. math::
            p_n(n|\lambda,t)= \exp\left[
            \ln \Gamma\left( n+q_\lambda \right)
            - \ln \Gamma\left( n+1 \right)
            - \ln \Gamma\left( q_\lambda \right)
            + q_\lambda\ln\left( \frac{p_\lambda}{t+p_\lambda} \right)
            + n\ln\left( \frac{t}{t+p_\lambda} \right) \right]

        :param n:
        :type n:
        :param t:
        :type t:
        :return:
        :rtype:
        """
        p_lambda = self.q_lambda / self.lamb
        den = t + p_lambda
        component1 = loggamma(n + self.q_lambda) - loggamma(self.q_lambda) - loggamma(n + 1)
        component2 = np.log(p_lambda / den) * self.q_lambda + np.log(t / den) * n
        return np.exp(component1 + component2)

    def ln_d_pmf(self, n, *args):
        t = get_event_occurrence_parameters(*args, default=1.0)
        p_lambda = self.q_lambda / self.lamb
        den = t + p_lambda
        component1 = loggamma(n + self.q_lambda) - loggamma(self.q_lambda) - loggamma(n + 1)
        component2 = np.log(p_lambda / den) * self.q_lambda + np.log(t / den) * n
        return component1 + component2

    def _d_cdf(self, n, t):
        cum_sum = 0.0
        for idx in range(n + 1):
            cum_sum += self._d_pmf(idx, t)
            if cum_sum >= 1.0:
                return 1.0
        return cum_sum

    def _d_rvs(self, t):
        u = np.random.uniform()
        cum_sum = 0.0
        n = 0
        while True:
            cum_sum += self._d_pmf(n, t)
            if cum_sum > u:
                return n
            n += 1
            if n > 10000:
                print("WARNING! Random many events, n > 10000")
                return n

    def _l_cdf(self, m, t):
        r"""
        .. math::
            F_M^{max}\left( m|\lambda,t \right)
            \left[ \frac{q_\lambda}{q_\lambda+\lambda t S_M\left( m \right)} \right]^{q_\lambda}

        :param self:
        :type self:
        :param m:
        :type m: float
        :param t:
        :type t: float
        :return:
        :rtype:
        """

        return ((self.q_lambda + self.lamb * t * self.magnitude_distribution.sf(m)) / self.q_lambda) ** (-self.q_lambda)

    def _l_pdf(self, m, t):
        r"""
        .. math::
            f_M^{max}\left( m|\lambda,t \right)=
            \frac{\lambda t q_\lambda f_M\left( m \right)F_M^{max}\left( m|\lambda,t \right)}
            {q_\lambda+\lambda t S_M\left( m \right)}

        :param m:
        :type m: float
        :param t:
        :type t: float
        :return:
        :rtype:
        """
        nominator = np.multiply(self.lamb * t * self.q_lambda * self.magnitude_distribution.pdf(m), self._l_cdf(m, t))
        return np.divide(nominator, self.q_lambda + self.lamb * t * self.magnitude_distribution.sf(m))

    # def _t_cdf(self, t, m):
    #     pass

    def _const_coefficients(self):
        return ['q_lambda']

    def _grad_sf_magnitude_distribution(self, m, t):
        n = t * self.lamb
        return n * (1.0 + n * self.magnitude_distribution.sf(m) / self.q_lambda) ** (-1.0 - self.q_lambda)

    def _grad_sf(self, m, t):
        d_lambda = t * (1.0 + t * self.lamb * self.magnitude_distribution.sf(m) / self.q_lambda) ** (
                    -1.0 - self.q_lambda)
        d_lambda = np.multiply(self.magnitude_distribution.sf(m), d_lambda)
        return {'lambda': d_lambda}

    def _d_grad_sf_magnitude_distribution(self, n, m, t):
        r"""

        :param m:
        :param n:
        :param t:
        :return:
        """
        raise Exception(f"Undefined _d_grad_sf_magnitude_distribution in the {self.local_name} class")

    def _d_grad_sf(self, n, t):
        r"""

        :param m:
        :param n:
        :param t:
        :return:
        """
        raise Exception(f"Undefined _d_grad_sf in the {self.local_name} class")
