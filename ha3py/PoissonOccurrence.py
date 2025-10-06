"""
Poisson events occurrence probability
-------------------------------------

..
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
from ha3py.LambdaOccurrence import LambdaOccurrence
from scipy.stats import poisson
from ha3py.BaseOccurrence import get_event_occurrence_parameters


class PoissonOccurrence(LambdaOccurrence, ABC):
    r"""
    The PDF of Poisson distribution of not exceeding magnitude :math:`m`
    in time :math:`t` is

    .. math::
        f_M^{max}\left( m|\lambda,\mathbf{\Theta_\beta},t \right)=
        \lambda t f_m\left( m|\mathbf{\Theta_\beta} \right)
        \exp\left( \lambda t S_M\left( m|\mathbf{\Theta_\beta} \right) \right).

    The CDF of Poisson distribution of not exceeding magnitude :math:`m`
    in time :math:`t` is

    .. math::
        F_M^{max}\left( m|\lambda,\mathbf{\Theta_\beta},t \right)=
        1 - e^{-\lambda t S_M\left( m|\mathbf{\Theta_\beta} \right)}.

    Gradients of the of non occurrence of event in the time probability survive function are:

    .. math::
        \frac{\partial S_M^{max}\left(m | t\right)}{\partial\lambda}=
        \exp\left(-\lambda t \left\{S_M\left( m|\mathbf{\Theta_\beta} \right)\right\} \right)
        t\left\{S_M\left( m |\mathbf{\Theta_\beta} \right) \right\},

        \frac{\partial S_M^{max}\left(m | t\right)}{\partial\theta_i}=
        \lambda t \exp\left(-\lambda t \left\{S_M\left( m|\mathbf{\Theta_\beta} \right)\right\} \right)
        \left\{\frac{\partial S_M\left(m \right)}{\partial\theta_i}\right\}.

    """
    def __init__(self, configuration, **kwargs):
        super().__init__(configuration, 'Poisson', **kwargs)

    def _d_cdf(self, n, t):
        return poisson.cdf(n, self.lamb * t)

    def _d_pmf(self, n, t):
        wyn = poisson.pmf(n, self.lamb * t)
        l = self.lamb * t
        return wyn

    def ln_d_pmf(self, n, *args):
        t = get_event_occurrence_parameters(*args, default=1.0)
        return poisson.logpmf(n, self.lamb * t)

    def _d_rvs(self, t):
        return poisson.rvs(self.lamb * t)

    def _l_cdf(self, m, t):
        r"""

        .. math::
            F_M^{max}\left( m|\lambda,t \right)=
            1-\exp\left(-\lambda t S_M\left( m \right)\right)

        :param m:
        :type m:
        :param t:
        :type t:
        :return:
        :rtype:
        """
        lambda_t = self.lamb * t
        return np.exp(-lambda_t * self.magnitude_distribution.sf(m))

    def _l_pdf(self, m, t):
        r"""

        .. math::
            f_M^{max}\left( m|\lambda,t \right)=
            \lambda t f_M\left( m \right) \exp\left(-\lambda t S_M\left( m \right)\right)

        :param m:
        :type m:
        :param t:
        :type t:
        :return:
        :rtype:
        """
        lambda_t = self.lamb * t
        return np.multiply(lambda_t * self.magnitude_distribution.pdf(m),
                           np.exp(-lambda_t * self.magnitude_distribution.sf(m)))

    def time_cdf(self, t, m):
        if t < 0:
            return 0.0
        l = self.lamb
        sf = self.magnitude_distribution.sf(m)
        fm = self.magnitude_distribution.pdf(m)
        ex = np.exp(-sf * t * l)
        return -(-fm + fm*ex + sf*t*fm*l*ex) / fm

    def _const_coefficients(self):
        return []

    def _grad_sf_magnitude_distribution(self, m, t):
        r"""

        .. math::

            \frac{\partial S_M^{max}\left(m | t\right)}{\partial S_M\left( m|\mathbf{\Theta_\beta} \right)}=
            \lambda t \exp\left(-\lambda t S_M\left( m|\mathbf{\Theta_\beta} \right) \right)

        :param m:
        :param t:
        :return:
        """
        n = t * self.lamb
        return n * np.exp(-n * self.magnitude_distribution.sf(m))

    def _grad_sf(self, m, t):
        d_lambda = t * np.multiply(np.exp(-t * self.lamb * self.magnitude_distribution.sf(m)),
                               self.magnitude_distribution.sf(m))
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

        :param n:
        :param t:
        :return:
        """
        raise Exception(f"Undefined _d_grad_sf in the {self.local_name} class")
