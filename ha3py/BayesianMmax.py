"""
Module base classes of bayesian maximum magnitude likelihood
------------------------------------------------------------

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
import numpy as np
from scipy.stats import rv_continuous
from scipy.stats import truncnorm
import scipy.integrate as integrate
from ha3py.m_max_utils import non_bayesian_m_max_estimation
from abc import ABC, abstractmethod
from ha3py.constant_values import EPS2
from ha3py.utils import HaPyException


def init_bayesian_m_max(configuration, magnitude_distribution=None, m_max_pair=None):
    if m_max_pair is None:
        m_max, sd_m_max = non_bayesian_m_max_estimation(configuration,
                                                        magnitude_distribution=magnitude_distribution)
    else:
        m_max, sd_m_max = m_max_pair
    prior_m_max = configuration['prior_m_max']
    sd_prior_m_max = configuration['sd_prior_m_max']
    if m_max > prior_m_max:
        print(f"Prior m_max ({prior_m_max}) is smaller than m_max estimated from the catalog ({m_max})")
        # raise HaPyException('Prior m_max error')
        exit(-1)
    return m_max, sd_m_max, prior_m_max, sd_prior_m_max

class BayesianBase(ABC):
    r"""
    The base class likelihood probability for most Bayesian methods.
    It estimates the non-bayesian :math:`\hat{m}_{max}` and
    defines the truncated normal distribution (:math:`\pi\left( m_{max} \right)`) of apriori :math:`m_{max}`.
    The derived classes must define the `_likelihood` method

    :param configuration: General configuration container,
        which is the dictionary of all parameters required for Ha3Py modules
        and results of all computations.
    :type configuration: dict

    """
    def __init__(self, configuration, magnitude_distribution=None, m_max_pair=None):
        """
        Initialisation

        :param configuration: General configuration container,
            which is the dictionary of all parameters required for Ha3Py modules
            and results of all computations.
        :type configuration: dict

        """
        self.m_max, self.sd_m_max, self.prior_m_max, self.sd_prior_m_max = init_bayesian_m_max(
            configuration, magnitude_distribution=magnitude_distribution, m_max_pair=m_max_pair)
        # self.m_max, self.sd_m_max = non_bayesian_m_max_estimation(configuration,
        #                                                           magnitude_distribution=magnitude_distribution)
        # self.prior_m_max = configuration['prior_m_max']
        # self.sd_prior_m_max = configuration['sd_prior_m_max']
        self.m_max_u = configuration.get('upper_m_max', 9.5)
        # self.m_max_l = configuration.get('lower_m_max', configuration['m_max_obs'])
        self.m_max_l = configuration.get('lower_m_max', self.m_max)
        a = (self.m_max_l - self.prior_m_max) / self.sd_prior_m_max
        b = (self.m_max_u - self.prior_m_max) / self.sd_prior_m_max
        self.pi = truncnorm(a, b, loc=self.prior_m_max, scale=self.sd_prior_m_max)

        # rv_continuous.__init__(self, a=self.m_max_l-EPS2, b=self.m_max_u+EPS2, name='BbS',
        #                        longname='Bayesian_by_shift')
        self.den = 1
        # self.den = self.cdf(self.m_max_u) - self.cdf(self.m_max_l)
        den =  integrate.quad(lambda x: self.pdf(x), self.m_max_l, self.m_max_u)
        self.den = den[0]

    def pdf(self, m, *args):
        """

        :param m:
        :type m:
        :return:
        :rtype:
        """
        # print(f"{m}: {self._likelihood(m)}, {self.pi.pdf(m)}")
        return np.multiply(self._likelihood(m), self.pi.pdf(m)) / self.den

    @abstractmethod
    def _likelihood(self, m):
        """

        :param m:
        :type m:
        :return:
        :rtype:
        """
        pass
