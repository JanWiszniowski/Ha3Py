"""
BaseDelta class of magnitude distributions
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

"""

from scipy.stats import rv_continuous
from abc import ABC, abstractmethod
from ha3py.constant_values import EPS2


class BaseMagnitudeDistribution(rv_continuous, ABC):
    r"""
    BaseDelta magnitude distribution class manages :math:`m_{min}` and :math:`m_{max}`.
    The base magnitude distribution is a derived class from the
    `SciPy generic continuous random variable class.
    <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rv_continuous.html>`_

    :param configuration: dictionary of Ha3Py params params
        Required items in the params dictionary are unless they are define in the constructor:

        * m_min
        * m_max_current
        * m_max (required if m_max_current is missing in the params dictionary)

    :param name: Magnitude distribution prompt
    :type name: str
    :param long_name: Magnitude distribution long prompt
    :type long_name: str
    :param m_min: Minimum value of the magnitude distribution
        If missing, the maximum magnitude is taken from configuration
    :type m_min: float
    :param m_max: Maximum value of the magnitude distribution.
        If missing, the maximum magnitude is taken from configuration
    :type m_max: float

    Classes derived from the MagnitudeDistribution classes define exact magnitude distribution,
    e.g. Gutenberg-Richter magnitude distribution. They must define methods:

        * _prepare - preparation of probability distribution computation,
        * _parameters - return the list of all probability distribution parameters,
        * _grad_sf - return the survive function gradients of all probability distribution parameters (see grad_sf),
        * _const_coefficients - return list of names of defined unestimated coefficients,
        * _coefficient_names - return list of names of coefficients including m_max as a last coefficient,
        * _coefficient_values -  return list of values of coefficients including m_max,
        * _pdf - return probability density function of magnitude(s),
        * _cdf - return cumulative distribution function of magnitude(s),

    Required params if they are not define in the constructor:

        * m_min,
        * m_max_current,
        * m_max (required if 'm_max_current' is missing in the params dictionary)

   """

    def __init__(self, configuration, name, long_name=None, m_min=None, m_max=None):
        """

        :param configuration:
        :param name:
        :param long_name:
        :param m_min:
        :param m_max:
        """
        r"""
        Required params:
            m_min
            m_max_current
        Optional gr_parameters:
            m_max (if m_max_current is missing in the params dictionary)
        """
        if m_min is not None:
            self._m_min = m_min
        else:
            self._m_min = configuration.get('m_min', 0.0)
        if m_max is not None:
            self._m_max = m_max
        elif 'm_max_current' in configuration:
            self._m_max = configuration['m_max_current']
        else:
            self._m_max = configuration.get('m_max', 10.0)
        if long_name is None:
            long_name = name
        rv_continuous.__init__(self, a=self._m_min - EPS2, b=self._m_max + EPS2, name=name, longname=long_name)
        self._prepare()

    @property
    def m_min(self):
        """It is the minimum magnitude"""
        return self._m_min

    @m_min.setter
    def m_min(self, val):
        self._m_min = val
        self.rv_continuous.a = val - EPS2
        self._prepare()

    @m_min.getter
    def m_min(self):
        return self._m_min

    @property
    def m_max(self):
        """It is the minimum magnitude"""
        return self._m_max

    @m_max.setter
    def m_max(self, val):
        self._m_max = val
        self.rv_continuous.b = val + EPS2
        self._prepare()

    @m_max.getter
    def m_max(self):
        return self._m_max

    @abstractmethod
    def _prepare(self):
        raise Exception(f"Undefined _prepare in the {self.name} class")

    @abstractmethod
    def _grad_sf(self, m):
        raise Exception(f"Undefined _grad_sf in the {self.name} class")

    @abstractmethod
    def _coefficient_names(self):
        raise Exception(f"Undefined _coefficient_names in the {self.name} class")

    @abstractmethod
    def _coefficient_values(self):
        raise Exception(f"Undefined _coefficient_values in the {self.name} class")

    @property
    def coefficients(self):
        """They are magnitude distribution parameters"""
        return self._coefficient_values()

    @property
    def coefficient_names(self):
        """They are magnitude distribution parameters"""
        return self._coefficient_names()

    @abstractmethod
    def _const_coefficients(self):
        raise Exception(f"Undefined _const_coefficients in the {self.name} class")

    @property
    def const_coefficients(self):
        """They are magnitude distribution parameters"""
        return self._const_coefficients()

    def grad_sf(self, m, coefficient_name=None):
        r"""
        Compute gradients of magnitude distribution survive function parameters

        .. math::
            \frac{\partial S_M\left( m \right)}{\partial x_i}, i=1,...

        where a survive function :math:`S_M \left( m \right) = 1 - F_M \left( m \right)`
        and :math:`x_i, i = 1,...` are the magnitude distribution parameters .

        :param coefficient_name: The coefficient name for which gradient is calculated.
            If coefficient_name is empty method return dictionary of gradients of all coefficients.
        :type coefficient_name: str
        :param m: magnitude distribution of the survive function
        :type m: float
        :return: Dictionary of magnitude distribution parameters names and their gradients
            or value of gradient of the coefficient.
        :rtype: dict

        The magnitude distribution coefficients depend on the magnitude distribution
        """
        if coefficient_name is None:
            return self._grad_sf(m)
        # coefficient_names = self.coefficient_names
        grad = self._grad_sf(m)
        return grad[coefficient_name]
        # for idx, name in enumerate(coefficient_names):
        #     if name == coefficient_name:
        #         return grad[idx]
        # return None
