"""
Base classes of events occurrence probabilities
-----------------------------------------------

..
    :copyright:
        Jan Wiszniowski <jwisz@igf.edu.pl>,
        Andrzej Kijko <andrzej.kijko@up.ac.za>
    :license:
        GNU Lesser General Public License, Version 3
        (https://www.gnu.org/copyleft/lesser.html)
    :version 0.0.1:
        2025-02-01

"""

from scipy.stats import rv_continuous
import numpy as np
from abc import ABC, abstractmethod
from ha3py.get_magnitude_distribution import get_magnitude_distribution
from ha3py.constant_values import EPS2


def get_event_occurrence_parameters(*args, default=0.0):
    if len(args) > 0:
        return args[0]
    else:
        return default


class OccurrenceBase(rv_continuous, ABC):
    r"""
    The base class of event occurrence.
    It defines two probabilities:

    *   The discrete probability :math:`p_n(n|t)` of occurrence of :math:`n` event
        equal or greater than the minimum magnitude in the period :math:`t`

        **methods:** d_pmf, d_cdf, ln_d_pmf, d_grad_sf,

    *   The continues probability that in the period :math:`time` magnitude of none event exceed the value :math:`m`

        **methods:** pdf, cdf, grad_sf, and other methods of the SciPy continues probability (rv_continuous).

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict
    :param name: The name of the derived class object
    :type name: str
    :param kwargs: Optional parameters are: magnitude_distribution (magnitude distribution object),
        m_min (float), m_max (float), theta (list or np.array).
        They are used for the object creation. If they are missing,
        appropriate objects are created based on the configuration
    :type kwargs: dict


    **Methods:**


    `_cdf`:
        Local definition of the cumulative distribution function :math:`F_M^{max}\left( m|t, \mathbf{\Theta} \right)`
        of continues probability that in the period :math:`t` magnitude of none event exceed the value :math:`m`,
        where coefficients :math:`\mathbf{\Theta}`, which depend on the specific realisation of the distribution,
        must be defined before the method call. E.g, If you wanted to compute occurrence
        for magnitudes :math:`m` or greater, you must set :math:`m_{min}=m` earlier.

        It modifies the `rv_continuous.pdf` and others rv_continuous's methods that requre the pdf.
        Derived classes shouldn't overwrite this method but rather define the abstract `_l_cdf`.
    `_pdf`:
        Local definition of the probability distribution function :math:`f_M^{max}\left( m|t, \mathbf{\Theta} \right)`
        of continues probability that in the period :math:`t` magnitude of none event exceed the value :math:`m`,
        where coefficients :math:`\mathbf{\Theta}`, which depend on the specific realisation of the distribution,
        must be defined before the method call. E.g, If you wanted to compute occurrence
        for magnitudes :math:`m` or greater, you must set :math:`m_{min}=m` earlier.
        It modifies the `rv_continuous.pdf` and others rv_continuous's methods that requre the pdf.
        Derived classes shouldn't overwrite this method but rather define the abstract `_l_pdf`.
    `d_cdf`:
        Definition of the cumulative distribution function :math:`F_n\left( n|t, \mathbf{\Theta} \right)`
        of discrete probability of occurrence of :math:`n` events equal or greater than the minimum magnitude
        in the period :math:`t`, where coefficients :math:`\mathbf{\Theta}`,
        which depend on the specific realisation of the distribution,
        must be defined before the method call. E.g, If you wanted to compute occurrence
        for magnitudes :math:`m` or greater, you must set :math:`m_{min}=m` earlier.

        Derived classes shouldn't overwrite this method but rather define the abstract `_d_cdf`.
    `d_pmf`:
        Definition of the discrete probability mass function :math:`p_n(n|t, \mathbf{\Theta} )`
        of occurrence of :math:`n` events equal or greater than the minimum magnitude in the period :math:`t`
        where coefficients :math:`\mathbf{\Theta}`, which depend on the specific realisation of the distribution,
        must be defined before the method call. E.g, If you wanted to compute occurrence
        for magnitudes :math:`m` or greater, you must set :math:`m_{min}=m` earlier.
        Derived classes shouldn't overwrite this method but rather define the abstract `_d_cdf`.
    `d_expected`:
        Definition of the expected value of discrete probability :math:`p_n(n|t)`
        of occurrence of :math:`n` event in the period :math:`t`:
    `ln_d_pmf`:
        The function compute the natural probability logarithm :math:`\ln\left[ p_n(n|t) \right]`
        of occurrence of :math:`n` events in the period :math:`t` probability mass function.
    `grad_sf`:
        Compute gradients of non occurrence of event in the time survive function (see `_cdf`, `_pdf`).

        .. math::

            \frac{\partial S_M^{max}\left(m | t\right)}{\partial\theta_{\lambda i}}

        and

        .. math::

            \frac{\partial S_M^{max}\left(m | t\right)}{\partial\theta_{\beta i}} =
            \frac{\partial S_M^{max}\left(m | t\right)}{\partial S_M\left( m|\mathbf{\Theta_\beta} \right)}
            \frac{\partial S_M\left(m | \mathbf{\Theta_\beta} \right)}{\partial\theta_{\beta i}},

        where :math:`S_M^{max}\left(m | t\right) = 1-F_M^{max}\left(m | t\right)`,
        :math:`\theta_{\beta i} \in \Theta_\beta` are magnitude distribution coefficients,
        and :math:`\theta_{\lambda i} \in \Theta_\lambda` are the occurrence probability own coefficients.

    `d_grad_sf`:
        Compute gradients of the discrete probability that, in the time :math:`t`, :math:`n`
        appears with magnitude :math:`m` or greater survive function (see `d_cdf`, `d_pmf`).

        .. math::

            \frac{\partial S_n\left(n | t\right)}{\partial\theta_{\lambda i}}

        and

        .. math::

            \frac{\partial S_n\left(n | t\right)}{\partial\theta_{\beta i}} =
            \frac{\partial S_n\left(n | t\right)}{\partial S_M\left( m|\mathbf{\Theta_\beta} \right)}
            \frac{\partial S_M\left(m | \mathbf{\Theta_\beta} \right)}{\partial\theta_{\beta i}},

        where :math:`S_n\left(n | t\right) = 1-F_n^{max}\left(n | t\right)`,
        :math:`\theta_{\beta i} \in \Theta_\beta` are magnitude distribution coefficients,
        and :math:`\theta_{\lambda i} \in \Theta_\lambda` are the occurrence probability own coefficients.


    **Properties:**

    `coefficients`:  (read only)
        Returns the list of all the occurrence probability estimated coefficients values :math:`\mathbf{\Theta}`.
        Both the continues probability :math:`f_M^{max}\left( m|t, \mathbf{\Theta} \right)`
        that in the period :math:`t` magnitude of none event exceed the value :math:`m` and
        the discrete probability :math:`p_n(n|t, \mathbf{\Theta} )`
        of occurrence of :math:`n` event in the period :math:`t` have the same coefficients.
        First in the list are own occurrence probability coefficients :math:`\Theta_\lambda`,
        next are magnitude distribution coefficients :math:`\Theta_\beta`.
        The list ends with the maximum magnitude
    `coefficient_names`: (read only)
        Returns the list of all the occurrence probability estimated coefficients names.
        Both the continues probability :math:`f_M^{max}\left( m|t, \mathbf{\Theta} \right)`
        that in the period :math:`t` magnitude of none event exceed the value :math:`m` and
        the discrete probability :math:`p_n(n|t, \mathbf{\Theta} )`
        of occurrence of :math:`n` event in the period :math:`t` have the same coefficients.
        First in the list are own occurrence probability coefficients,
        next are magnitude distribution coefficients.
        The list ends with the maximum magnitude
    `const_coefficients`:  (read only)
        Returns the names of all the occurrence probability constant coefficients.
        These coefficients must be defined (some values can be default) but are non-estimated
    `m_min`:
        The substitution sets the minimum magnitude :math:`m_{min}` of occurrence probability
        and modifies the probability.
    `m_max`:
        The substitution sets the maximum magnitude :math:`m_{max}` of occurrence probability
        and modifies the probability.

    **Required definition of abstract methods:**

    All the following methods must be defined in derived classes

    `_l_cdf`:
        Returns in the derived class
        the cumulative distribution function :math:`F_M^{max}\left( m|t \right)`
        of continues probability that in the period :math:`t` magnitude of none event exceed the value :math:`m`

    `_l_pdf`:
        Returns in the derived class
        the probability distribution function :math:`f_M^{max}\left( m|t \right)`
        of continues probability that in the period :math:`t` magnitude of none event exceed the value :math:`m`,

    `_d_cdf`:
        Returns in the derived class
        the cumulative distribution function :math:`F_n\left( n|t \right)`
        of discrete probability of occurrence of :math:`n` event equal or greater than the minimum magnitude
        in the period :math:`t`.

    `_d_pmf`:
        Returns in the derived class
        the discrete probability :math:`p_n(n|t )` of occurrence of :math:`n` event equal or greater
        than the minimum magnitude in the period :math:`t` mass function
    `_coefficient_names`:
        Returns the list of own occurrence probability estimated coefficients :math:`\mathbf{\Theta_{\lambda}}` names
        excluding the magnitude distribution coefficients :math:`\mathbf{\Theta_{\beta}}`.
    `_coefficient_values`:
        Returns the list of own occurrence probability estimated coefficients :math:`\mathbf{\Theta_{\lambda}}`values
        excluding the values of magnitude distribution coefficients :math:`\mathbf{\Theta_{\beta}}`.
    `_const_coefficients`:
        Returns the names of own occurrence probability constant (non-estimated) coefficients
        excluding the names of constant magnitude distribution coefficients.
    `_grad_sf_magnitude_distribution(m, t)`:
        Returns the gradient of the of non occurrence of event in the time probability survive function
        with respect to magnitude distribution  survive function:

        .. math::

            \frac{\partial S_M^{max}\left(m | t\right)}{\partial S_M\left( m|\mathbf{\Theta_\beta} \right)}

    `_d_grad_sf_magnitude_distribution(n, m, t)`:
        Return the gradient of the discrete probability that, in the time interval :math:`t`, :math:`n`
        appears with magnitude :math:`m` or greater survive function
        with respect to magnitude distribution  survive function:

        .. math::

            \frac{\partial S_n\left(n | m, t\right)}{\partial S_n\left( m|\mathbf{\Theta_\beta} \right)}

    `_grad_sf(m, t)`:
        Return gradients of the non occurrence of event in the time survive function
        with respect to the occurrence probability own coefficients

        .. math::

            \frac{\partial S_M^{max}\left(n | t\right)}{\partial\theta}

    `_d_grad_sf(n, t)`:
        Return gradients of the discrete probability that, in the time :math:`t`, :math:`n`
        appears with magnitude :math:`m` or greater survive function
        with respect to the occurrence probability own coefficients

        .. math::

            \frac{\partial S_n\left(n | t\right)}{\partial\theta}

    `d_mean`:

    """

    def __init__(self, configuration, name, **kwargs):
        """

        :param configuration:
        :type configuration:
        :param name:
        :type name:
        :param kwargs:
        :type kwargs:
        """
        self.local_name = name
        magnitude_distribution = kwargs.get('magnitude_distribution')
        if magnitude_distribution is None:
            magnitude_distribution = get_magnitude_distribution(configuration,
                                                                m_min=kwargs.get('m_min'),
                                                                m_max=kwargs.get('m_max'),
                                                                theta=kwargs.get('theta'))
        self.magnitude_distribution = magnitude_distribution
        rv_continuous.__init__(self, name=name, longname=f'{name}_md_{self.magnitude_distribution.name}',
                               a=magnitude_distribution.a, shapes='t')

    # -------------------------------------------------------------
    # Non occurrence probability of magnitude in time

    @abstractmethod
    def _l_cdf(self, m, t):
        pass

    def _cdf(self, m, *args):
        r"""


        :param m:
        :param args:
        :return:
        """
        t = get_event_occurrence_parameters(*args, default=1.0)
        return self._l_cdf(m, t)

    @abstractmethod
    def _l_pdf(self, m, t):
        pass

    def _pdf(self, m, *args):
        r"""


        :param m:
        :param args:
        :return:
        """
        t = get_event_occurrence_parameters(*args, default=1.0)
        return self._l_pdf(m, t)

    # -------------------------------------------------------------
    # Probability of occurrence od n event

    @abstractmethod
    def _d_cdf(self, n, t):
        pass

    def d_cdf(self, n, *args):
        r"""
        The function compute the cumulate distribution function :math:`F_n(n|t)` of occurrence of :math:`n`
        event in the period :math:`t`. If you wanted to compute occurrence for magnitudes :math:`m` or greater,
        you must set :math:`m_{min}=m` earlier.

        :param n: Number of events
        :type n: int
        :param args: Optional period value. If no argument missing, one year is assumed.
        :type args: list
        :return: The CDF for :math:`n` events
        :rtype: float

        """
        t = get_event_occurrence_parameters(*args, default=1.0)
        return self._d_cdf(n, t)

    @abstractmethod
    def _d_pmf(self, n, t):
        pass

    def d_pmf(self, n, *args):
        """
        The function compute the  probability mass function :math:`p_n(n|t)` of occurrence of :math:`n`
        events in the period :math:`t`. If you wanted to compute occurrence for magnitudes :math:`m` or greater,
        you must set :math:`m_{min}=m` earlier.

        :param n: Number of events
        :type n: int
        :param args: Optional period value. If no argument missing, one year is assumed.
        :type args: list
        :return: The probability of :math:`n` event occurrence in the period :math:`time`
        :rtype: float

        """
        t = get_event_occurrence_parameters(*args, default=1.0)
        return self._d_pmf(n, t)

    @abstractmethod
    def _d_mean(self, t):
        pass

    def d_mean(self, *args):
        """

        :return: The
        :rtype: float
        """
        t = get_event_occurrence_parameters(*args, default=1.0)
        return self._d_mean(t)

    def d_expected(self, t):
        """
        The function compute the expected value
        of occurrence of :math:`n` events in the period :math:`time` having magnitudes :math:`m` or greater.
        It is realised by :math:``
        :param t:
        :return:
        """
        x = 0.0
        p = 100.0
        n = 0
        s = 0.0
        while p > 0.00001 * x:
            p = self._d_pmf(n, t)
            s += n * p
            n += 1
            if p > x:
                x = p
        return s

    def ln_d_pmf(self, n, *args):
        r"""
        The function compute the natural probability logarithm :math:`\ln\left [ p_n(n|t) \right ]`
        of occurrence of :math:`n` events in the period :math:`t`  probability mass function.
        If you wanted to compute occurrence for magnitudes :math:`m` or greater,
        you must set :math:`m_{min}=m` earlier.

        :param n:
        :param args:
        :return:
        """
        return np.log(self.d_pmf(n,*args))

    @abstractmethod
    def _d_rvs(self, *args):
        pass

    def d_rvs(self, *args):
        t = get_event_occurrence_parameters(*args, default=1.0)
        return self._d_rvs(t)
    # -------------------------------------------------------------
    # coefficient

    @abstractmethod
    def _coefficient_names(self):
        raise Exception(f'Undefined _coefficient_names in the {self.local_name} class')

    @abstractmethod
    def _coefficient_values(self):
        raise Exception(f'Undefined _coefficient_values in the {self.local_name} class')

    @property
    def coefficients(self):
        r"""
        The coefficients function gets the values of current event occurrence parameters.
        They are given in order: first are own event occurrence parameters e.g. :math:`\lambda`,
        next are magnitude distribution parameters e.g. :math:`\beta`,
        and at the end is maximum value parameters :math:`m_{max}`.

        :return: The list of current event occurrence parameters
        :rtype: list(float)
        """
        return self._coefficient_values() + self.magnitude_distribution.coefficients

    @property
    def coefficient_names(self):
        """
        The coefficient_names function gets the names of variable event occurrence parameters.
        The names correspond to the coefficient values of the coefficients function.
        They are given in order: first are own event occurrence parameters names e.g. '*lambda*',
        next are magnitude distribution parameters names e.g. '*beta*',
        and at the end is maximum value parameters '*m_max*'.

        :return: The list of current event occurrence parameters names
        :rtype: list(str)
        """
        return self._coefficient_names() + self.magnitude_distribution.coefficient_names

    @abstractmethod
    def _const_coefficients(self):
        raise Exception(f'Undefined _const_coefficients in the {self.local_name} class')

    @property
    def const_coefficients(self):
        """
        The const_coefficients gets the values of constant (not estimated) event occurrence parameters

        :return: The list of constant event occurrence parameters
        :rtype: list(float)
        """
        """They are magnitude distribution Parameters"""
        return self.magnitude_distribution.const_coefficients + self._const_coefficients()

    # -------------------------------------------------------------
    # gradients

    @abstractmethod
    def _grad_sf_magnitude_distribution(self, m, t):
        raise Exception(f'Undefined _grad_sf_magnitude_distribution in the {self.local_name} class')

    @abstractmethod
    def _d_grad_sf_magnitude_distribution(self, n, m, t):
        raise Exception(f'Undefined _d_grad_sf_magnitude_distribution in the {self.local_name} class')

    @abstractmethod
    def _grad_sf(self, m, t):
        raise Exception(f'Undefined _grad_sf in the {self.local_name} class')

    @abstractmethod
    def _d_grad_sf(self, n, t):
        raise Exception(f'Undefined _d_grad_sf in the {self.local_name} class')

    def grad_sf(self, m, t):
        r"""
        Compute gradients of event non occurrence time survive function
        for all event non occurrence time probability coefficients
        including magnitude distribution coefficients and maxim magnitude

        .. math::
            \frac{\partial S_M\left( time | m \right)}{\partial \theta_i}, i=1,...

        where a survive function :math:`S_M^{max} \left( time | m \right) = 1 - F_M^{max} \left( time | m \right)`
        and :math:`\theta_i, i = 1,...` are the magnitude non occurrence time probability coefficients.

        :param t: The time event does not exceed the magnitude :math:`m`
        :type t: float
        :param m: The magnitude that will be not exceeded in the time :math:`time`
        :type m: float
        :return: The dictionary of variable events occurrence parameters names and their gradients.
                 The events occurrence parameters depend on the events occurrence
                 and magnitude distribution objects
        :rtype: dict
        """
        general_md_grad = self._grad_sf_magnitude_distribution(m, t)
        particular_md_grad = self.magnitude_distribution.grad_sf(m)
        for key in particular_md_grad.keys():
            particular_md_grad[key] *= general_md_grad
        return self._grad_sf(m, t) | particular_md_grad

    def d_grad_sf(self, n, m, t):
        r"""
        Compute gradients of the probability that, in the time interval :math:`time`, :math:`n`
        appears with magnitude :math:`m` or greater.

        .. math::
            \frac{\partial S_n\left( n | m,time \right)}{\partial \theta_i}, i=1,...

        where a survive function :math:`S_M \left( m \right) = 1 - F_M \left( m \right)`
        and :math:`x_i, i = 1,...` are the magnitude distribution Parameters .

        :param n: The number of events
        :type n: int
        :param t: The time event does not exceed the magnitude :math:`m`
        :type t: float
        :param m: The magnitude that will be not exceeded in the time :math:`time`
        :type m: float
        :return: The dictionary of variable events occurrence parameters names and their gradients.
                 The events occurrence parameters depend on the events occurrence
                 and magnitude distribution objects
        :rtype: dict
        """
        general_md_grad = self._d_grad_sf_magnitude_distribution(n, m, t)
        particular_md_grad = self.magnitude_distribution.grad_sf(m)
        for key in particular_md_grad.keys():
            particular_md_grad[key] *= general_md_grad
        return self._d_grad_sf(n, t) | particular_md_grad

    # -------------------------------------------------------------
    # m_min

    @property
    def m_min(self):
        """It is the minimum magnitude"""
        return self.magnitude_distribution.m_min

    @abstractmethod
    def _set_m_min(self, val):
        raise Exception(f'Undefined _set_m_min in the {self.local_name} class')

    @m_min.setter
    def m_min(self, val):
        self.a = val - EPS2
        self._set_m_min(val)
        self.magnitude_distribution.m_min = val

    @m_min.getter
    def m_min(self):
        return self.magnitude_distribution.m_min

    # -------------------------------------------------------------
    # m_max

    @property
    def m_max(self):
        """It is the minimum magnitude"""
        return self.magnitude_distribution.m_max

    @m_max.setter
    def m_max(self, val):
        self.magnitude_distribution.m_max = val

    @m_max.getter
    def m_max(self):
        return self.magnitude_distribution.m_max
