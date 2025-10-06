"""
BaseDelta class of magnitude distributions
------------------------------------------

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

    :param parameters: dictionary of Ha3Py params params
        Required items in the params dictionary are unless they are define in the constructor:

        * m_min
        * m_max_current
        * m_max (required if m_max_current is missing in the params dictionary)

    :param name: magnitude distribution prompt
    :param long_name:  magnitude distribution long prompt
    :param m_min: minimum value of the magnitude distribution
    :param m_max: maximum value of the magnitude distribution

    Classes derived from the MagnitudeDistribution classes define exact magnitude distribution,
    e.g. Gutenberg-Richter magnitude distribution. They must define methods:

    Required params if they are not define in the contractor:

    * m_min
    * m_max_current
    * m_max (required if 'm_max_current' is missing in the params dictionary)

    * _prepare - preparation of probability distribution computation
    * _parameters - return the list of all probability distribution Parameters
    * _grad_sf - return the survive function gradients of all probability distribution Parameters
        (see grad_sf)

   """
    def __init__(self, parameters, name, long_name=None, m_min=None, m_max=None):
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
            self._m_min = parameters.get('m_min', 0.0)
        if m_max is not None:
            self._m_max = m_max
        elif 'm_max_current' in parameters:
            self._m_max = parameters['m_max_current']
        else:
            self._m_max = parameters.get('m_max', 10.0)
        if long_name is None:
            long_name = name
        rv_continuous.__init__(self, a=self._m_min-EPS2, b=self._m_max+EPS2, name=name, longname=long_name)
        self._prepare()

    @property
    def m_min(self):
        """It is the minimum magnitude"""
        return self._m_min

    @m_min.setter
    def m_min(self, val):
        self._m_min = val
        self.a = val - EPS2
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
        self.b = val + EPS2
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
        """They are magnitude distribution Parameters"""
        return self._coefficient_values()

    @property
    def coefficient_names(self):
        """They are magnitude distribution Parameters"""
        return self._coefficient_names()

    @abstractmethod
    def _const_coefficients(self):
        raise Exception(f"Undefined _const_coefficients in the {self.name} class")

    @property
    def const_coefficients(self):
        """They are magnitude distribution Parameters"""
        return self._const_coefficients()

    def grad_sf(self, m, coefficient_name=None):
        r"""
        Compute gradients of magnitude distribution survive function Parameters

        .. math::
            \frac{\partial S_M\left( m \right)}{\partial x_i}, i=1,...

        where a survive function :math:`S_M \left( m \right) = 1 - F_M \left( m \right)`
        and :math:`x_i, i = 1,...` are the magnitude distribution Parameters .

        :param coefficient_name:
        :type coefficient_name:
        :param m: magnitude distribution of the survive function
        :return: dictionary of magnitude distribution Parameters names and their gradients.
                 The magnitude distribution Parameters depend on the magnitude distribution
        """
        if coefficient_name is None:
            return self._grad_sf(m)
        coefficient_names = self.coefficient_names
        grad = self._grad_sf(m)
        return grad[coefficient_name]
        # for idx, name in enumerate(coefficient_names):
        #     if name == coefficient_name:
        #         return grad[idx]
        # return None
