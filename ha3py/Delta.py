r"""
Delta (:math:`\Delta`) calculation classes
------------------------------------------

..
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
from ha3py.constant_values import EPS


class BaseDelta(ABC):
    r"""
    BaseDelta :math:`\Delta` calculation class.
    """

    def __init__(self, name, parameters, magnitude_distribution=None, m_max=None, m_max_obs=None):
        self.name = name
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

        :param n: number if events. If n is unset
            then it is determined based on t and :math:`\lambda`: :math:`n=t\lambda`
        :type n: float
        :param time: the time duration in years
        :type time: float
        :param annual_lambda: annual occurrence - :math:`\lambda` value
        :type annual_lambda: float
        :return: the :math:`\Delta` value - result of the virtual function _method_delta(n) for :math:`m=m_{max}`

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

    def exist_solution(self, n=None, time=1.0, annual_lambda=1.0):
        if n is None:
            n = time * annual_lambda
        return self._exist_solution(n)

    def _exist_solution(self, n):
        r"""
        The abstract :math:`\Delta` calculation method.
        :param n: number of events
        """
        return n > 0.0


class KijkoSellevoll(BaseDelta):
    r"""
    Kijko-Sellevoll :math:`\Delta` calculation class is:

    .. math::
        \Delta =\int_{m_{min}}^{m_{max}}F_M\left( m | m_{max} \right)^ndm

    The integration is performed numerically.
    """

    def __init__(self, configuration, magnitude_distribution=None, m_max=None, m_max_obs=None):
        super().__init__('Kijko-Sellevoll', configuration, magnitude_distribution=magnitude_distribution,
                         m_max=m_max, m_max_obs=m_max_obs)

    def _delta(self, n):
        r"""
        The Kijko-Sellevoll :math:`\Delta` calculation method.

        :param n: Number of events, It can be float value
        :return: the :math:`\Delta` value

        .. math::
            \Delta =\int_{m_{min}}^{m_{max}}F_M\left( m | m_{max} \right)^ndm

        The integration is performed numerically.
        """
        # return integrate.quad(lambda x: exp(log(self.cdf(x))*n), self.m_min, m)
        delta = integrate.quad(lambda x: self.magnitude_distribution.cdf(x) ** n,
                               self.magnitude_distribution.m_min, self.magnitude_distribution.m_max)
        return delta[0]

    def _exist_solution(self, n):
        r""""
        Test whether exist Tate-Pisarenko solution condition:

        .. math::
            ln\left(n\right)\geq\beta\left(m_{max}^{obs}{-m}_{min}\right)-0.58.

        """
        return True
        # left = log(n)
        # right = self.magnitude_distribution.beta * (self.m_max_obs - self.magnitude_distribution.m_min) - 0.5772156649
        # return left >= right


class KijkoSellevollMmaxObs(BaseDelta):
    r"""
    Kijko-Sellevoll :math:`\Delta` calculation class is:

    .. math::
        \Delta =\int_{m_{min}}^{m_{max}^{obs}}F_M\left( m | m_{max} \right)^ndm

    The integration is performed numerically.
    """

    def __init__(self, configuration, magnitude_distribution=None, m_max=None, m_max_obs=None):
        super().__init__('Kijko-Sellevoll for m_max_obs', configuration,
                         magnitude_distribution=magnitude_distribution, m_max=m_max, m_max_obs=m_max_obs)

    def _delta(self, n):
        r"""
        The Kijko-Sellevoll simplified :math:`\Delta` calculation method.

        :param n: Number of events, It can be float value
        :return: the :math:`\Delta` value

        .. math::
            \Delta =\int_{m_{min}}^{m_{max}^{obs}}F_M\left( m | m_{max} \right)^ndm

        The integration is performed numerically.
        """
        # return integrate.quad(lambda x: exp(log(self.cdf(x))*n), self.m_min, m)
        delta = integrate.quad(lambda x: self.magnitude_distribution.cdf(x) ** n,
                               self.magnitude_distribution.m_min, self.m_max_obs)
        return delta[0]


class TatePisarenko(BaseDelta):
    r"""
    Tate-Pisarenko :math:`\Delta` calculation class is:

    .. math::
        \Delta =\frac{1}{nf_M\left( m_{max} | m_{max} \right)}

    """

    def __init__(self, configuration, magnitude_distribution=None, m_max=None, m_max_obs=None):
        super().__init__('Tate-Pisarenko', configuration, magnitude_distribution=magnitude_distribution,
                         m_max=m_max, m_max_obs=m_max_obs)

    def _delta(self, n):
        r"""
        The Tate-Pisarenko :math:`\Delta` calculation method.

        :param n: Number of events, It can be float value
        :return: the :math:`\Delta` value

        .. math::
        \Delta =\frac{1}{nf_M\left( m_{max} | m_{max} \right)}

        """
        pdf_m = self.magnitude_distribution.pdf(self.magnitude_distribution.m_max - EPS)
        return 1.0 / n / pdf_m

    def _exist_solution(self, n):
        r""""
        Test whether exist Tate-Pisarenko solution condition:

        .. math::
            m_{max}^{obs}{-m}_{min}\le\frac{\ln{\left(n\right)}}{\beta}-\frac{n-1}{n\beta}

        """
        return True
        # left = self.m_max_obs - self.magnitude_distribution.m_min
        # right = log(n) / self.magnitude_distribution.beta - (n - 1) / n / self.magnitude_distribution.beta
        # return left <= right


class TatePisarenkoMmaxObs(BaseDelta):
    r"""
    Tate-Pisarenko simplified :math:`\Delta` calculation class is:

    .. math::
        \Delta =\frac{1}{nf_M\left( m_{max}^{obs} | m_{max} \right)}

    Please note that the original
    Tate-Pisarenko :math:`\Delta=\frac{1}{nf_M\left( m_{max}^{obs} | m_{max}^{obs} \right)}` calculation
    is called directly without solving the equation :math:`m_{max}=m_{max}^{obs}+\Delta(m_{max})`,
    since :math:`\Delta(m_{max})=\text{constant}`.

    """

    def __init__(self, configuration, magnitude_distribution=None, m_max=None, m_max_obs=None):
        super().__init__('Tate-Pisarenko m_max_obs', configuration, magnitude_distribution=magnitude_distribution,
                         m_max=m_max, m_max_obs=m_max_obs)

    def _delta(self, n):
        r"""
        The Tate-Pisarenko :math:`\Delta` calculation method.

        :param n: Number of events, It can be float value
        :return: the :math:`\Delta` value

        .. math::
        \Delta =\frac{1}{nf_M\left( m_{max}^{obs} \right)}

        """
        # Modification of pdf_m for stability ???
        # pdf_m = self.magnitude_distribution.pdf(self.m_max - 0.025)
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
    if delta == 'Kijko-Sellevoll m_max_obs':
        return KijkoSellevollMmaxObs(configuration, magnitude_distribution=magnitude_distribution, m_max=m_max)
    elif delta == 'Tate-Pisarenko m_max_obs':
        return TatePisarenkoMmaxObs(configuration, magnitude_distribution=magnitude_distribution, m_max=m_max)
    else:
        raise HaPyException('Unknown delta computation')
