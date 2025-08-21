r"""
Delta (:math:`\Delta`) calculation classes
------------------------------------------

:copyright:
    Jan Wiszniowski <jwisz@igf.edu.pl>,
    Andrzej Kijko <andrzej.kijko@up.ac.za>
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
:version 0.0.1:
    2025-01-01

The delta classes define the :math:`\Delta` calculation methods,
They are applied in a few :math:`m_{max}` estimation algorithms,
and realises an object-oriented approach to this issue.
It provides a very flexible approach to :math:`m_{max}` estimation,
allowing the assessment of multiple :math:`m_{max}` estimation algorithms
in various combinations.
The base :math:`\Delta` is

There exist two formulas for :math:`\Delta` calculation:

* based on Tate-Pisarenko theory,
* based on Kijko-Sellevoll theory.
"""

import scipy.integrate as integrate
from abc import ABC, abstractmethod
from ha3py.get_magnitude_distribution import get_magnitude_distribution
from ha3py.utils import HaPyException


class BaseDelta(ABC):
    r"""
    BaseDelta :math:`\Delta` calculation class.
    """
    def __init__(self, parameters, magnitude_distribution=None, m_max=None, m_max_obs=None):
        if magnitude_distribution:
            self.magnitude_distribution = magnitude_distribution
        else:
            self.magnitude_distribution = get_magnitude_distribution(parameters, m_max=m_max)
        if m_max_obs is None:
            self.m_max_obs = parameters.get('m_max_obs', 10.0)
        else:
            self.m_max_obs = m_max_obs

    def delta(self, n=None, time=1.0, annual_lambda=1.0):
        r"""
        Calculates the delta value. Instead of using `delta` the object can be called oneself. E.g::

            delta_object = KijkoSellevollDelta(params)
            delta = delta_object(t=123.0, annual_lambda=0.37)

        :param n: number if events. If n is unset then it is determined as t and :math:`\lambda`
             :math:`n=t\lambda`
        :param time: the time duration in years
        :param annual_lambda: annual occurrence - :math:`\lambda` value
        :return: the :math:`\Delta` value - result of the virtual function _method_delta(m, n)

        """
        if n is None:
            n = time * annual_lambda
        return self._delta(n)

    def __call__(self, n=None, time=1.0, annual_lambda=1.0):
        return self.delta(n=n, time=time, annual_lambda=annual_lambda)

    @property
    def m_max(self):
        """It is the minimum magnitude"""
        return self.magnitude_distribution.m_max

    @m_max.setter
    def m_max(self, val):
        self.magnitude_distribution.m_max = val

    @abstractmethod
    def _delta(self, n):
        r"""
        The abstract :math:`\Delta` calculation method.
        :param n: number of events
        """
        raise Exception('Undefined')


class KijkoSellevoll(BaseDelta):
    r"""
    Kijko-Sellevoll :math:`\Delta` calculation class

    .. math::
        \Delta =\int_{m_{min}}^{m_{max}}F_M\left( m \right)^ndm

    The integration is performed numerically.
    """

    def __init__(self, configuration, magnitude_distribution=None, m_max=None, m_max_obs=None):
        super().__init__(configuration, magnitude_distribution=magnitude_distribution, m_max=m_max, m_max_obs=m_max_obs)

    def _delta(self, n):
        r"""
        The Kijko-Sellevoll :math:`\Delta` calculation method.

        :param n: Number of events, It can be float value
        :return: the :math:`\Delta` value

        .. math::
            \Delta =\int_{m_{min}}^{m_{max}}F_M\left( m \right)^ndm

        The integration is performed numerically.
        """
        # return integrate.quad(lambda x: exp(log(self.cdf(x))*n), self.m_min, m)
        delta = integrate.quad(lambda x: self.magnitude_distribution.cdf(x) ** n,
                               self.magnitude_distribution.m_min, self.magnitude_distribution.m_max)
        return delta[0]


class TatePisarenko(BaseDelta):
    r"""
    Tate-Pisarenko :math:`\Delta` calculation class

    .. math::
        \Delta =\frac{1}{nf_M\left( m_{max}^{obs} \right)}

    """

    def __init__(self, configuration, magnitude_distribution=None, m_max=None, m_max_obs=None):
        super().__init__(configuration, magnitude_distribution=magnitude_distribution, m_max=m_max, m_max_obs=m_max_obs)

    def _delta(self, n):
        r"""
        The Tate-Pisarenko :math:`\Delta` calculation method.

        :param n: Number of events, It can be float value
        :return: the :math:`\Delta` value

        .. math::
        \Delta =\frac{1}{nf_M\left( m_{max}^{obs} \right)}

        """
        # Modification of pdf_m for stability ???
        if self.m_max_obs < self.magnitude_distribution.m_max - 0.025:
            pdf_m = self.magnitude_distribution.pdf(self.m_max_obs)
        else:
            pdf_m = self.magnitude_distribution.pdf(self.magnitude_distribution.m_max - 0.025)
        return 1.0 / n / pdf_m


def get_delta(configuration, magnitude_distribution=None, m_max=None):
    delta = configuration.get('delta', 'Kijko-Sellevoll')
    if delta == 'Kijko-Sellevoll':
        return KijkoSellevoll(configuration, magnitude_distribution=magnitude_distribution, m_max=m_max)
    elif delta == 'Tate-Pisarenko':
        return TatePisarenko(configuration, magnitude_distribution=magnitude_distribution, m_max=m_max)
    else:
        raise HaPyException('Unknown delta computation')
